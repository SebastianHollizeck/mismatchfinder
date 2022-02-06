import datetime
from logging import debug, error, info, basicConfig, getLogger
from multiprocessing import Process
from sys import getsizeof
from collections import defaultdict

import pysam
import re
from numpy import array, quantile, sort, mean

from mismatchfinder.results.Results import MismatchCandidates
from mismatchfinder.core.EndMotives import EndMotives
from mismatchfinder.utils.Misc import countLowerCase
from mismatchfinder.utils.Output import outputExists

# TODO: remove after profiling
# from pyinstrument import Profiler


class BamScanner(Process):
    """Class which scans a Bam for TMB and other interesting features"""

    def __init__(
        self,
        semaphore,
        lock,
        results,
        bamFilePath,
        referenceFile,
        minMQ,
        minBQ,
        minAvgBQ,
        maxMisMatchesPerRead,
        minMisMatchesPerRead,
        maxMisMatchesPerFragment,
        minMisMatchesPerFragment,
        fragmentLengthIntervals,
        filterSecondaries=True,
        onlyOverlap=False,
        strictOverlap=False,
        blackList=None,
        whiteList=None,
        germObj=None,
        afCutOff=0,
        germlineRequirePass=False,
        outFileRoot=None,
        kmer=4,
        log=None,
        writeEvidenceBam=False,
        writeEvidenceReadPairs=False,
        overwrite=False,
    ):
        super(BamScanner, self).__init__()

        self.semaphore = semaphore
        self.results = results
        self.lock = lock

        self.minMQ = minMQ
        self.minBQ = minBQ
        self.minAvgBQ = minAvgBQ

        self.maxMisMatchesPerRead = maxMisMatchesPerRead
        self.minMisMatchesPerRead = minMisMatchesPerRead

        self.maxMisMatchesPerFragment = maxMisMatchesPerFragment
        self.minMisMatchesPerFragment = minMisMatchesPerFragment

        self.fragmentLengthIntervals = fragmentLengthIntervals

        self.filterSecondaries = filterSecondaries

        self.onlyOverlap = onlyOverlap
        self.strictOverlap = strictOverlap

        self.bamFilePath = bamFilePath
        self.referenceFile = referenceFile

        self.blackList = blackList
        self.whiteList = whiteList
        self.germObj = germObj
        self.afCutOff = afCutOff
        self.germlineRequirePass = germlineRequirePass

        self.outFileRoot = outFileRoot
        self.writeEvidenceBam = writeEvidenceBam
        self.writeEvidenceReadPairs = writeEvidenceReadPairs
        self.overwrite = overwrite

        self.endMotives = EndMotives(kmer)

        self.logger = log
        self.logLevel = log.level

        # this pattern will match the MDString of reads with at least 1, up to a maximum of
        # maxMisMatchesPerRead mismatches in a read (for both SBS and DBS)
        self.MDstrPattern = re.compile(
            r"\d+([ACGT](0[ACGT])?\d+){"
            + str(minMisMatchesPerRead)
            + ","
            + str(maxMisMatchesPerRead)
            + "}",
        )

    # this function gets all sites of mismatches from any mapped read in a bam
    # use qualThreshold, bedObj and minMQ to exclude reads and mismatches
    # qualThreshold is for the base quality
    # bedObj is to discard reads mapping to those blacklisted areas in the bed
    # minMQ directly ignores reads which have a mappingquality lower
    def getMutationSites(self):

        # state your purpose ;)
        self.logger.info("Checking reads for mismatches")
        # store the found sites of mutations and how often they occured
        mutSites = {}

        # store counts for stats later
        nLowQualReads = 0
        nNoMisMatchReads = 0
        nReads = 0
        nBlackListedReads = 0
        nMisMatches = 0
        nAlignedBases = 0
        nAlignedReads = 0
        nSecondaryHits = 0
        nFragSize = 0
        nAnalysedBases = 0

        # store the fragment lengths for later
        fragLengths = []

        # get the time we started with this
        startTime = datetime.datetime.now()
        currTime = datetime.datetime.now()

        # we work directly on the bam without iterator creation
        bamFile = pysam.AlignmentFile(
            self.bamFilePath, "r", reference_filename=self.referenceFile
        )

        # possibly write the evidence reads to a file as well
        if not self.outFileRoot is None and self.writeEvidenceBam:
            evBamPath = f"{self.outFileRoot.parent}/{self.outFileRoot.name}_{self.bamFilePath.name}_evidence.bam"
            evidenceBam = pysam.AlignmentFile(evBamPath, "wb", template=bamFile)
            self.logger.info(f"Writing evidence reads to {evBamPath}")
        else:
            evidenceBam = None

        # store the reads for which we do not have a partner yet
        readCache = {}
        # cache the wrapper for black and whitelist
        isReadInRegionOfInterest = self.isReadInRegionOfInterest

        # Profiling revealed, the most time is spent in isReadInRegionOfInterest closely followed
        # by the mean function
        # prof = Profiler()
        # prof.start()
        for read in bamFile.fetch(until_eof=True):
            nReads += 1
            # we check every 10K reads how uch time has passed and if it has been more than 30
            # we give an update
            if nReads % 100000 == 0:

                now = datetime.datetime.now()
                # time since the last output
                currDelta = now - currTime

                # if it has been 10 seconds we last said something
                if currDelta.total_seconds() > 10:

                    # time since start
                    deltaTime = now - startTime

                    readsPerSec = nReads / deltaTime.total_seconds()

                    # we really just care about the general aread so we round to the next hundred
                    readsPerSec = int(round(readsPerSec, -2))

                    # give some info how far we are already through the bam
                    self.logger.info(
                        f"Read through {nReads} reads - processing {readsPerSec:9d} reads per second"
                    )
                    currTime = now

            # we only want proper reads and no secondaries. We can be pretty lenient here, because
            # we check for the mismatch itself if the sequencing quality is high enough later.
            qname = read.query_name
            reads = []
            if read.is_qcfail or read.is_secondary or read.is_supplementary:
                # in this case, we dont want to do anything with the read
                nLowQualReads += 1
            elif (
                read.is_unmapped
                or read.is_duplicate
                or read.mapping_quality < self.minMQ
                or mean(read.query_qualities) < self.minAvgBQ
            ):
                # in this case, we care about the info, that this end of the fragment is not used
                # we add in this None value, so we have nothing stuck in the case where one does
                # not pass the filter but the mate does not
                if read.is_paired:
                    if not qname in readCache:
                        # we only bother to add this, if the other read has a chance of an overlap
                        if read.reference_id == read.next_reference_id:
                            readCache[qname] = None
                            continue
                    else:
                        # in this case, this IS the mate and if this is a proper read, we need to
                        # analyse it
                        mate = readCache[qname]
                        if not mate is None:
                            reads = [mate]
                        del readCache[qname]

                nLowQualReads += 1
            else:
                # get the blacklist status
                readWhiteListed = isReadInRegionOfInterest(read)

                # if the read is paired, there is the chance of a
                # overlap, so we instead store the read until we find the mate
                # so that we can build a consensus
                if read.is_paired and read.reference_id == read.next_reference_id:
                    if qname not in readCache:
                        readCache[qname] = read
                        continue
                    else:
                        mate = readCache[qname]
                        # see the blacklist staus of the mate
                        mateWhiteListed = isReadInRegionOfInterest(mate)

                        # if the mate is high quality its not none and both of the reads need
                        # to be in the region of interest for an overlap to be within the region
                        # of interest and if the consensus part isnt in the region of interest it
                        # does not change anyhing and we can save computational time
                        if not mate is None and (readWhiteListed and mateWhiteListed):
                            read1, read2, interLen = makeConsensusRead(
                                read,
                                mate,
                                onlyOverlap=self.onlyOverlap,
                                strict=self.strictOverlap,
                            )
                            # if there was no overlap, but onlyOverlap was true, the reads will
                            # both be None, so we can skip them
                            if read1 is None:
                                reads = []
                            else:
                                reads = [read1, read2]
                                # we store how big the region of analysis was (this is only relevant if we have onlyOverlap enabled)
                                nAnalysedBases += interLen
                        else:

                            # at this point only one of the two reads can be in the region of
                            # interest
                            if readWhiteListed:
                                # in this case we do only analyse this read by itself
                                reads = [read]
                            if mateWhiteListed:
                                # in this case we do only analyse the mate
                                reads = [mate]

                            # this means one of the reads  was high quality, but was not in a
                            # region of interest
                            nBlackListedReads += 1

                        del readCache[qname]

                elif readWhiteListed:
                    # if its not a pair, we just use it as is if it is in the region of interest
                    reads = [read]
                else:
                    # this means the read was high quality, but either not whitelisted or was
                    # blacklisted
                    nBlackListedReads += 1

            # prefilter the reads in case we only want to analyse the reads if BOTH reads pass
            # all filters
            scanList = []
            for r in reads:

                # even if the fragment doesnt have any mismatches it is important to
                # store the fragment length of this read for fragment size statistics, but
                # only for the first read, because otherwise we just have everything twice
                if not r.is_paired or r.is_read1:
                    fragLengths.append(abs(r.template_length))
                # single end reads will not have a template length
                # that is meaningfull for us, but its good to have
                # the info that it was single end in the final
                # result as
                # (number discordant reads == number of reads)

                # add the amount of bases of this read that were aligned
                nAlignedBases += r.query_alignment_length
                nAlignedReads += 1

                # we also care about the endmotives of the reads, so we store those
                self.endMotives.count(r)

                if not hasNMisMatches(r, self.MDstrPattern):
                    nNoMisMatchReads += 1
                elif self.filterSecondaries and hasSecondaryMatches(r):
                    nSecondaryHits += 1
                elif not readInFragmentLength(r, self.fragmentLengthIntervals):
                    nFragSize += 1
                else:
                    scanList.append(r)

            # if we only want to look at the overlap and want to be strict about it, we stop if
            # one of the reads got filtered out
            # we could probably save LOTS of time if we were to only look at one of the two here,
            # but until i can recalculate the MD tag of a read, i am stuck with checking both
            if self.onlyOverlap and self.strictOverlap and len(scanList) != 2:
                continue

            # now we execute the mismatch finding for the reads that we selected (read + mate)
            # to enable the check for mismatches per fragment, we save the mismatches instead of
            # evaluating right away
            tmpMisMatches = []
            nFragMisMatches = 0
            tmpPerReadUsage = []
            for r in scanList:

                # get all mismatches in this read
                readMisMatches = self.scanAlignedSegment(r)

                # because we only checked with the MD string so far, but didnt really actually check
                # if we then end up with the right amount of mismatches, we do that here.
                nReadMisMatches = len(readMisMatches)

                # we really only need to check, if we still have enough, as the MDstr already gave
                # us the information that it will not be MORE than that. (this is slightly more
                # efficient and we need every help we can get)
                if nReadMisMatches >= self.minMisMatchesPerRead:
                    tmpMisMatches += readMisMatches
                    # just saving this should be faster than calculating the length again for the
                    # joined list (and still correct)
                    nFragMisMatches += nReadMisMatches
                    # store in the per read, so we can write the evidence only for the reads we used
                    tmpPerReadUsage.append(True)
                else:
                    tmpPerReadUsage.append(False)

            # then we check if between the two reads (or even just the single if the other one was
            # discard) we have enough mismatches to keep this in the analysis
            # in contrast to the perRead check, we can definitly have too many here and need to
            # check both boundaries
            if (
                nFragMisMatches >= self.minMisMatchesPerFragment
                and nFragMisMatches <= self.maxMisMatchesPerFragment
            ):
                # store the mismatches and keep a record how often each was found
                for mm in tmpMisMatches:
                    if mm in mutSites:
                        mutSites[mm] += 1
                    else:
                        mutSites[mm] = 1
                    # we also store the amount of mismatches found, so we can calculate
                    # the mismatches per read which should be stable between samples
                    nMisMatches += 1

                # write the evidence bam to file, for the reads we used
                if not evidenceBam is None:
                    for r, use in zip(scanList, tmpPerReadUsage):
                        # only write the read if it was used, or if we want the pair info
                        if use or self.writeEvidenceReadPairs:
                            evidenceBam.write(r)

        # at this point, there should be no more reads in our read storage, or we did something
        # wrong (or the bam is truncated in some test case)
        if len(readCache) != 0:
            self.logger.error(
                f"Found reads left over in the cache after the main loop, this is a bug or shows a problem with orphaned reads in the input\nListing read names:"
            )
            counter = 0
            for read in readCache:
                self.logger.error(read.query_name)
                counter += 1
                if counter == 10:
                    self.logger.error(f"and {len(readCache)-10} more")

        # we are done so we update the status as well
        self.logger.info(
            f"Read through 100.00% of reads in {(datetime.datetime.now()-startTime).total_seconds()/60:.1f} minutes"
        )

        bamFile.close()
        if not evidenceBam is None:
            evidenceBam.close()

        # TODO remove after profiling
        # prof.stop()
        # print(prof.output_text(unicode=True, color=True))
        # prof.output_html()

        # did we have an issue with reads?
        if nAlignedBases == 0:
            self.logger.error(
                f"Could not detect any reads from Bam ({self.bamFilePath}). Further analysis is not possible\nLowQualReads: {nLowQualReads}\nNoMisMatchReads: {nNoMisMatchReads}\nBlacklistedReads: {nBlackListedReads}"
            )
            raise Exception(f"No reads left for {self.bamFilePath}")

        # sort the fragment lengths so we can get the median and the sides of the distribution for
        # estimation of fusions or similar (because those would have higher fragment lengths as
        # expected )
        fragLengths = array(sort(fragLengths))

        lenFragLengths = len(fragLengths)
        # we get the middle of the array, which should be the median of the distribution
        median = fragLengths[int(lenFragLengths / 2)]

        # find the index, at which the read distance gets to normal amounts
        lowSizeIdx = 0
        for i in range(0, lenFragLengths):
            if fragLengths[i] < 35:
                lowSizeIdx = i
            else:
                break

        # find the index, at which the read distance gets to normal amounts
        highSizeIdx = lenFragLengths
        for i in range((lenFragLengths - 1), 0, -1):
            if fragLengths[i] > 1500:
                highSizeIdx = i
            else:
                break

        # amount of discordant reads is the amount of lower than 35 plus the amount of reads longer
        # than 1500 bases
        nDiscordantFragments = lowSizeIdx + (lenFragLengths - highSizeIdx)

        # get the fragment length distribution
        quantileRange = [0.01, 0.05, 0.15, 0.25, 0.5, 0.75, 0.85, 0.95, 0.99]
        fragLenQuantilesAr = quantile(fragLengths, quantileRange)
        fragLenQuantiles = {}
        for i in range(0, len(quantileRange)):
            fragLenQuantiles[quantileRange[i]] = fragLenQuantilesAr[i]

        # if we have only overLap, we need to only look at the intersection rather than all of
        # the aligned bases, but we stored that already in the variant, so we overwrite it if
        # its not overlap only
        if not self.onlyOverlap:
            nAnalysedBases = nAlignedBases

        # return a dict of the counts we made
        return MismatchCandidates(
            mutSites=mutSites,
            nReads=nReads,
            nLowQualReads=nLowQualReads,
            nNoMisMatchReads=nNoMisMatchReads,
            nBlackListedReads=nBlackListedReads,
            nAlignedReads=nAlignedReads,
            fragmentSizeQuantiles=fragLenQuantiles,
            nDiscordantFragments=nDiscordantFragments,
            nMisMatches=nMisMatches,
            nAnalysedBases=nAnalysedBases,
            bam=self.bamFilePath.name,
        )

    # this is a wrapper for the check for black and whitelist
    def isReadInRegionOfInterest(self, read):
        if read is None:
            # if we get a None, we have to assume its blacklisted
            return False
        # if we have a blacklist and the read fits, then we have a blacklisted read
        if not self.blackList is None and self.blackList.isReadWithinRegion(read):
            return False
        # if we have a whitelist and the read does NOT overlap, then its blacklisted
        elif not self.whiteList is None and not self.whiteList.isReadWithinRegion(read):
            return False
        else:
            # well otherwise its not blacklisted and can be used
            return True

    def run(self):

        # for some reason, we need to reset the logLevel here otherwise it resets to WARNING
        # maybe it has to do with using "spawn" as a start method
        self.logger.setLevel(self.logLevel)

        # first we check if we actually need to create output
        if not self.overwrite and outputExists(self.outFileRoot, self.bamFilePath):
            self.logger.info(
                f"Skipping scan of {self.bamFilePath.name} found previous results"
            )
        else:
            self.logger.info(f"Starting scan of {self.bamFilePath.name}")
            # initiate bam object

            # execute
            # first get all possible sites
            mutCands = self.getMutationSites()

            # filter out germline mismatches
            mutCands.checkGermlineStatus(
                self.germObj,
                afCutOff=self.afCutOff,
                discard=True,
                requirePass=self.germlineRequirePass,
            )

            # count each context (strand agnosticly, so if the rreverse complement is already found we
            # instead count the reverse complement)
            mutCands.countContexts()

            if not self.outFileRoot is None:
                # we request only one thread to write to the files at one time
                with self.lock:
                    self.logger.info(f"Writing result counts to file")
                    mutCands.writeSBSToFile(self.outFileRoot, self.bamFilePath)
                    mutCands.writeDBSToFile(self.outFileRoot, self.bamFilePath)
                    mutCands.writeStatsToFile(self.outFileRoot, self.bamFilePath)
                    self.endMotives.writeFreqsToFile(self.outFileRoot, self.bamFilePath)

                # dont need a lock here, as we create a new file for each bam
                mutCands.writeSitesToFile(self.outFileRoot, self.bamFilePath)
            else:
                self.logger.info(f"No output file root found so no output written")

            # sadly it takes VERY long to actually put things into the result queue if the sites
            # copied as well... the issus is the pickling. So we really need to work with this within
            # this process
            mutCands.mutSites = None
            mutCands.fragmentSizeQuantiles = None
            ## TODO: write out mutSites for "panel of normals"
            ## TODO: fit a density function for the fragmentSizes?
            ## TODO: maybe also "delete" the field instead of set to None
            self.logger.info(f"Finished scan of {self.bamFilePath.name}")

        self.logger.debug("Releasing requested resource lock")
        # release the block for resources again
        self.semaphore.release()

        # explicitly destroy the semaphore and lock
        del self.semaphore
        del self.lock

    # read through one read and find all mismatches and their contexts
    # can visualize the contexts on the reads
    def scanAlignedSegment(self, AlignedSegment, vis=False):

        # this is where we will store the sites we find in this read
        mutations = []

        # we do this here, so we only do it once instead of in the loop
        alignedRefSequence = AlignedSegment.get_reference_sequence()
        referencePositions = array(
            AlignedSegment.get_reference_positions(full_length=True)
        )
        # this is a store to see where in the reference sequence we are, because its not possible to
        # calculate back from the aligned pairs
        # we start with -1 so that we can increase it in the beginning and not in the end of the
        # loop
        rpos = -1
        query_qualities = AlignedSegment.query_qualities
        # loop through the read
        for (readPos, contigPos, seq) in AlignedSegment.get_aligned_pairs(
            with_seq=True, matches_only=False
        ):
            # if the contigPos is None, that means we have an insertion (or softclipping) which we
            # skip without increasing the reference sequence pos
            if contigPos is None:
                continue

            # if readPos is none, that means its a deletion, which means the refpos need to be
            # advanced and we can skip this, because its not a mismatch
            rpos += 1
            if readPos is None:
                continue

            # if it is lower case it symbolises a mismatch
            if seq.islower():

                # because we might have deleted this by creating a consensus, we check if the query_sequence is equal to what we actually have
                # @TODO: proper fix for MD tag when building consensus, because
                # that would remove the need of this tag alternatively, build
                # our own class, which holds all the info we need and we can
                # manipulate
                if seq.upper() == AlignedSegment.query_sequence[readPos]:
                    # # need to fix this, because the ref is calculated from
                    # the MD string which we kinda fucked up by changing
                    # the sequence
                    tmpRef = list(alignedRefSequence)

                    # mdStr = AlignedSegment.get_tag("MD")
                    # we need actually cant really do this properly, because there isnt a good way
                    # to get back to the position in the template if there is a deletion
                    tmpRef[rpos] = seq.upper()

                    alignedRefSequence = "".join(tmpRef)
                    # but then we just skip this corrected snp
                    continue

                # offset the soft cliping at the beginning
                readPos = readPos - AlignedSegment.query_alignment_start
                # offset the reference location (this could probably be substituted by rpos but i
                # dont have enough time to debug it, so this stays)
                templatePos = contigPos - AlignedSegment.reference_start

                # if everything is right, we get the trinucl context of the mismatch in the reference
                # and the query
                altContext = AlignedSegment.query_alignment_sequence[
                    readPos - 1 : readPos + 2
                ]
                refContext = alignedRefSequence[templatePos - 1 : templatePos + 2]

                # if either of those contexts is not the full 3 nucleotides, we will skip to the next
                if len(altContext) != 3 or len(refContext) != 3:
                    continue

                # now we check if what we have is actually from a continous stretch, or if there is an
                # indel hidden within the snp
                # as we use the full length refpositions list we need to add the softclip positions back
                refIndex = readPos + AlignedSegment.query_alignment_start
                # we want an increment of 1 from the position left to middle and one more again for
                # the position to the right
                if (
                    referencePositions[refIndex - 1] != referencePositions[refIndex] - 1
                    or referencePositions[refIndex + 1]
                    != referencePositions[refIndex] + 1
                ):
                    continue

                # now that we have a proper context defined we can go about and check what type of variant this

                # if we have more than 1 lower case char in our context, that means we ae at the
                # beginning of a multi mismatch region and we should combine the two into a multi
                # variant thing
                # so here we just move on to the next position if we find there is a mismatch in the
                # next position
                # this will also take care of trinucleotide variants
                if refContext[2].islower():
                    continue
                elif refContext[0].islower():
                    # this means its a DBS
                    misMatchClass = 2
                    # in this case, we check if the average of the two quals is higher than the
                    # min bq
                    # if the quality of the base is too low, we drop this
                    if (
                        query_qualities[refIndex] + query_qualities[refIndex - 1]
                    ) / 2 < self.minBQ:
                        continue
                else:
                    # and this is a SBS
                    misMatchClass = 1
                    # if the quality of the base is too low, we drop this
                    if query_qualities[refIndex] < self.minBQ:
                        continue

                # however we would get the last position of a 3bp mismatch as a dinucleotide change,
                # which we do not want so we also check if there are more than 2 variants within 4
                # positions
                if (
                    countLowerCase(
                        alignedRefSequence[templatePos - 2 : templatePos + 2]
                    )
                    > 2
                ):
                    continue

                # that would leave us with only single or dinucleotide variants which are shifted to
                # the left e.g.
                # alt: TTG     alt: CTG
                # ref: ccG     ref: CcG
                # which we can easily convert to the classes that we want which are the C>T changes in
                # the dinucleotide context CC,CT,TC because those are specific for melanoma
                mut = (
                    AlignedSegment.reference_name,
                    contigPos,
                    refContext,
                    altContext,
                    misMatchClass,
                )

                # just so we can check what kind of variants we actually put into this we can print
                # the contexts for easy readability (this does not put the context at the same
                # position if there are indels in the read but at the position it is taken from)
                if vis:
                    # print the query (the actual read) on top
                    self.logger.debug(AlignedSegment.query_alignment_sequence)
                    line = ""
                    for i in range(0, readPos - 1):
                        line += " "
                    line += altContext
                    self.logger.debug(line)
                    # print the template (the reference below)
                    self.logger.debug(alignedRefSequence, addTime=False)
                    line = ""
                    for i in range(0, readPos - 1):
                        line += " "
                    line += refContext
                    self.logger.debug(line)

                # we do one last check if the actual mismatch is in the region
                if (
                    not self.blackList is None
                    and self.blackList.isMisMatchWithinRegion(mut)
                ):
                    # we do nothing here, because the mismatch is outside the region of interest
                    append = False
                    debug(
                        f"read: {AlignedSegment.query_name} is not blacklisted, but mut {mut} is"
                    )

                elif self.whiteList is None or self.whiteList.isMisMatchWithinRegion(
                    mut
                ):
                    # add to the found mutations now that everything is sorted out
                    mutations.append(mut)

        return mutations


# this is a simple check if a MDstring contains a certain number of mismatches
# we make this a method, so we can use the same precompiled pattern
def hasNMisMatches(read, pattern):
    # if there is no MD tag we cant do this, so we just say it does not have mismatches
    try:
        mdStr = read.get_tag("MD")
    except KeyError:
        # this occures only in unmapped reads so we cant actually say if there are mismatches
        return False

    # use the mdstr pattern, which allows a set number of mismatches, if the pattern doesnt
    # match it will return none
    if re.fullmatch(pattern, mdStr) is None:
        return False
    else:
        return True


# we might need to discard read, which also map to different locations
def hasSecondaryMatches(read):
    # if there is no XA tag (KeyError), then there is no secondary
    # but if it doesnt fail, we have secondary hits
    try:
        xa = read.get_tag("XA")
        return True
    except KeyError:
        return False


# calculates the consensus of two reads and can also set the non overlapping part to base quality 0
# so that it will not be analysed further
def makeConsensusRead(read1, read2, onlyOverlap=False, strict=False):

    # we can just check, if there are overlaps from just the start and end pos, which is significant
    # ly faster than using the md str
    # we use the equals, because the end is non inclusive
    if (
        read1.reference_end <= read2.reference_start
        or read1.reference_start >= read2.reference_end
    ):
        if onlyOverlap:
            return (None, None, None)
        return (read1, read2, abs(read1.template_length))

    # get the reference positions to see if read1 and 2 are aligned to the same area
    read1RefPos = read1.get_reference_positions(full_length=True)
    read2RefPos = read2.get_reference_positions(full_length=True)
    # get the intersection of the reference positions of the two reads
    inter = set(read1RefPos).intersection(read2RefPos)

    # we do this only if there is an intersection, so that we speed things up if there isnt
    if len(inter) > 0:
        read1Seq = list(read1.query_sequence)
        read1Quals = read1.query_qualities
        read1IndDict = dict((k, i) for i, k in enumerate(read1RefPos))

        read2Seq = list(read2.query_sequence)
        read2Quals = read2.query_qualities
        read2IndDict = dict((k, i) for i, k in enumerate(read2RefPos))
    else:
        # if there is no intersection, we just return the reads unchanged (this shouldnt happen,
        # because we checked the template length before, but sure)
        if onlyOverlap:
            return (None, None, none)
        return (read1, read2, abs(read1.template_length))

    # go through all of the overlaps and decide which of the reads is better
    for pos in inter:
        # we might get a None, if there was softclipping, so we discard those
        if pos is None:
            continue
        read1IntPos = read1IndDict[pos]
        read2IntPos = read2IndDict[pos]
        # we only really care if there is a difference in the sequence
        if read1Seq[read1IntPos] != read2Seq[read2IntPos]:

            # if the strict method is used, overlaps that dont show the same base
            # are discarded (quality set to 0)
            if strict:
                read1Quals[read1IntPos] = 0
                read2Quals[read2IntPos] = 0
                continue

            # we see which has the higher quality and then take that as the ground truth, but also
            # we try to avoid counting the overlap twice by reducing the quality of the other to 0
            # but also, because they did not agree, we slightly reduce the quality of the kept read
            if read1Quals[read1IntPos] > read2Quals[read2IntPos]:
                read2Seq[read2IntPos] = read1Seq[read1IntPos]
                read1Quals[read1IntPos] = read1Quals[read1IntPos] - (
                    int(read2Quals[read2IntPos] / 2)
                )
                read2Quals[read2IntPos] = 0
            elif read1Quals[read1IntPos] < read2Quals[read2IntPos]:
                read1Seq[read1IntPos] = read2Seq[read2IntPos]
                read2Quals[read2IntPos] = read2Quals[read2IntPos] - (
                    int(read1Quals[read1IntPos] / 2)
                )
                read1Quals[read1IntPos] = 0
            else:
                # this is the case where both are likely, in this case we just
                # use the reference that one of them hopefully has?

                # this is the ref base of read1, we can use either, because we
                # know they overlap here, but we need to make it upper, because
                # if read1 is the one with the mismatch we get a lowercase

                # we kinda use this to check if an md string is present
                readForBase = None
                try:
                    read1.get_tag("MD")
                    readForBase = read1
                except KeyError:
                    # now we try the seond read
                    try:
                        read2.get_tag("MD")
                        readForBase = read2
                    except KeyError:
                        # well, we tried without an md str we cant build a consensus
                        return (read1, read2)

                for (readPos, contigPos, refBase) in readForBase.get_aligned_pairs(
                    with_seq=True, matches_only=True
                ):
                    if contigPos == pos:
                        refBase = refBase.upper()
                        break

                if read1Seq[read1IntPos] == refBase:
                    read2Seq[read2IntPos] = refBase
                    read2Quals[read2IntPos] = read1Quals[read1IntPos]
                elif read2Seq[read2IntPos] == refBase:
                    read1Seq[read1IntPos] = refBase
                    read1Quals[read1IntPos] = read2Quals[read2IntPos]
                else:
                    # at this point we just take whatever they say
                    pass
        else:
            # this is a quick and ugly fix to not count mismatches from the two reads corresponding
            # to one dna fragment twice. If both reads are reference, we dont care anyways, but if
            # both show a mismatch, we only want to count one, we also sum the readQuals up,
            # because thats what mpileup does as well
            if read1Quals[read1IntPos] >= read2Quals[read2IntPos]:
                read1Quals[read1IntPos] = (
                    read1Quals[read1IntPos] + read2Quals[read2IntPos]
                )
                read2Quals[read2IntPos] = 0
            elif read1Quals[read1IntPos] < read2Quals[read2IntPos]:
                read2Quals[read2IntPos] = (
                    read1Quals[read1IntPos] + read2Quals[read2IntPos]
                )
                read1Quals[read1IntPos] = 0

    # if we only made a consensus, we report the normal length
    fLen = abs(read1.template_length)
    if onlyOverlap:
        # here we set every position that was not in the intersection to BQ 0
        for k, intPos in enumerate(read1RefPos):
            if not intPos in inter:
                read1Quals[k] = 0
        for k, intPos in enumerate(read2RefPos):
            if not intPos in inter:
                read2Quals[k] = 0
        # we only report the length of the intersection
        fLen = len(inter)

    # finally we have to create new reads
    read1New = read1
    read1New.query_sequence = "".join(read1Seq)
    read1New.query_qualities = read1Quals
    # same for read2
    read2New = read2
    read2New.query_sequence = "".join(read2Seq)
    read2New.query_qualities = read2Quals

    return (read1New, read2New, fLen)


def readInFragmentLength(read, fragmentLengthIntervals):
    # check all the fragmentlength intervals we have if this fits
    fragSize = abs(read.template_length)
    for min, max in fragmentLengthIntervals:
        if fragSize >= min or fragSize <= max:
            return True

    # we return fals if none of them worked
    return False
