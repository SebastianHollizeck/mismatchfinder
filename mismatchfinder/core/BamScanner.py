import datetime
from logging import debug, error, info, basicConfig, getLogger
from multiprocessing import Process
from sys import getsizeof
from collections import defaultdict

import pysam
from numpy import array, quantile, sort, mean

from mismatchfinder.results.Results import MismatchCandidates
from mismatchfinder.core.EndMotives import EndMotives
from mismatchfinder.utils.Misc import countLowerCase

# TODO: remove after profiling
from pyinstrument import Profiler


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
        blackList=None,
        whiteList=None,
        germObj=None,
        outFileRoot=None,
        kmer=4,
        log=None,
    ):
        super(BamScanner, self).__init__()

        self.semaphore = semaphore
        self.results = results
        self.lock = lock

        self.minMQ = minMQ
        self.minBQ = minBQ
        self.minAvgBQ = minAvgBQ
        self.maxMisMatchesPerRead = maxMisMatchesPerRead

        self.bamFilePath = bamFilePath
        self.referenceFile = referenceFile

        self.blackList = blackList
        self.whiteList = whiteList
        self.germObj = germObj

        self.outFileRoot = outFileRoot

        self.endMotives = EndMotives(kmer)

        self.logger = log
        self.logLevel = log.level

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

        # store the fragment lengths for later
        fragLengths = []

        # get the time we started with this
        startTime = datetime.datetime.now()
        currTime = datetime.datetime.now()

        # we work directly on the bam without iterator creation
        bamFile = pysam.AlignmentFile(
            self.bamFilePath, "r", reference_filename=self.referenceFile
        )

        # store the reads for which we do not have a partner yet
        readCache = {}
        # cache the wrapper for black and whitelist
        isReadInRegionOfInterest = self.isReadInRegionOfInterest

        # @TODO: disable this when done with profiling
        prof = Profiler()
        prof.start()
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
                        f"Read through {nReads} reads - processing {readsPerSec:6d} reads per second"
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
                            reads = [readCache[qname]]
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
                            read1, read2 = makeConsensusRead(read, mate)
                            reads = [read1, read2]
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

            # now we execute the mismatch finding for the reads that we selected (read + mate)
            for r in reads:

                if not hasMisMatches(r):
                    nNoMisMatchReads += 1
                else:
                    # get all mismatches in this read
                    tmpMisMatches = self.scanAlignedSegment(r)
                    # store the mismatches and keep a record how often each was found
                    for mm in tmpMisMatches:
                        if mm in mutSites:
                            mutSites[mm] += 1
                        else:
                            mutSites[mm] = 1
                    # we also store the amount of mismatches found, so we can calculate
                    # the mismatches per read which should be stable between samples
                    nMisMatches += len(tmpMisMatches)

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

        # at this point, there should be no more reads in our read storage, or we did something
        # wrong (or the bam is truncated in some test case)
        if len(readCache) != 0:
            self.logger.error(
                f"Found reads left over in the cache after the main loop, this is a bug {readCache}"
            )

        # we are done so we update the status as well
        self.logger.info(
            f"Read through 100.00% of reads in {(datetime.datetime.now()-startTime).total_seconds()/60:.1f} minutes"
        )
        # TODO remove after profiling
        prof.stop()
        print(prof.output_text(unicode=True, color=True))
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

        # we get the middle of the array, which should be the median of the distribution
        median = fragLengths[int(len(fragLengths) / 2)]

        # find the index, at which the read distance gets to normal amounts
        lowSizeIdx = -1
        for i in range(0, len(fragLengths)):
            if fragLengths[i] < 35:
                lowSizeIdx = i
            else:
                break

        # find the index, at which the read distance gets to normal amounts
        highSizeIdx = -1
        for i in range((len(fragLengths) - 1), 0, -1):
            if fragLengths[i] > 1500:
                highSizeIdx = i
            else:
                break

        # amount of discordant reads is the amount of lower than 35 plus the amount of reads longer than
        # 1500 bases
        nDiscordantReads = lowSizeIdx + (nAlignedReads - highSizeIdx)

        # get the fragment length distribution
        quantileRange = [0.01, 0.05, 0.15, 0.25, 0.5, 0.75, 0.85, 0.95, 0.99]
        fragLenQuantilesAr = quantile(fragLengths, quantileRange)
        fragLenQuantiles = {}
        for i in range(0, len(quantileRange)):
            fragLenQuantiles[quantileRange[i]] = fragLenQuantilesAr[i]

        # return a dict of the counts we made
        return MismatchCandidates(
            mutSites=mutSites,
            nReads=nReads,
            nLowQualReads=nLowQualReads,
            nNoMisMatchReads=nNoMisMatchReads,
            nBlackListedReads=nBlackListedReads,
            nAlignedReads=nAlignedReads,
            fragmentSizeQuantiles=fragLenQuantiles,
            nDiscordantReads=nDiscordantReads,
            nMisMatches=nMisMatches,
            nAlignedBases=nAlignedBases,
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

        self.logger.info(f"Starting scan of {self.bamFilePath.name}")
        # initiate bam object

        # execute
        # first get all possible sites
        mutCands = self.getMutationSites()

        # filter out germline mismatches
        mutCands.checkGermlineStatus(self.germObj)

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

        self.logger.debug("Releasing requested resource lock")
        # release the block for resources again
        self.semaphore.release()
        self.logger.info(f"Finished scan of {self.bamFilePath.name}")

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
        refIndDict = dict(
            (k, i) for i, k in enumerate(AlignedSegment.get_reference_positions())
        )

        # loop through the read
        for (readPos, contigPos, seq) in AlignedSegment.get_aligned_pairs(
            with_seq=True, matches_only=True
        ):

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
                    mappedPos = readPos - AlignedSegment.query_alignment_start
                    mdStr = AlignedSegment.get_tag("MD")
                    debug(
                        f"correcting pos {mappedPos} of read {AlignedSegment.qname} with md {mdStr} from {seq} to {seq.upper()} after read error correction with mate"
                    )
                    tmpRef[refIndDict[contigPos]] = seq.upper()

                    alignedRefSequence = "".join(tmpRef)
                    # but then we just skipr this corrected snp
                    continue

                # we really only want high quality mismatches
                # we need to do this before we do the softclip adjustment
                # because the base qualities are not soft clipped
                qual = AlignedSegment.query_qualities[readPos]

                # if the quality of the base is too low, we drop this
                if qual < self.minBQ:
                    continue

                # offset the soft cliping at the beginning
                readPos = readPos - AlignedSegment.query_alignment_start
                # offset the reference location
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
                else:
                    # and this is a SBS
                    misMatchClass = 1

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
                    debug(f"read: {read} is not blacklisted, but mut {mut} is")

                elif self.whiteList is None or self.whiteList.isMisMatchWithinRegion(
                    mut
                ):
                    # add to the found mutations now that everything is sorted out
                    mutations.append(mut)
                    # if we go over the allowed amount of mismatches in this read, we return an
                    # empty set of mismatches instead
                    if len(mutations) > self.maxMisMatchesPerRead:
                        return []

        return mutations


# this is a simple check if a MDstring contains any mismatches
def hasMisMatches(read):
    # if there is no MD tag we cant do this, so we just say it does not have mismatches
    try:
        mdStr = read.get_tag("MD")
    except KeyError:
        # this occures only in unmapped reads so we cant actually say if there are mismatches
        return False

    try:
        # if the whole string is just an int it tells us all positions were matches and no
        # mismatches
        mdInt = int(mdStr)
        return False
    except ValueError:
        # but if there is a ValueError we know there is more in there than that (could be an indel
        # as well, but we have to look closer to actually check that)
        return True


def makeConsensusRead(read1, read2):

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
        # if there is no intersection, we just return the reads unchanged
        return (read1, read2)

    # go through all of the overlaps and decide which of the reads is better
    for pos in inter:
        # we might get a None, if there was softclipping, so we discard those
        if pos is None:
            continue
        read1IntPos = read1IndDict[pos]
        read2IntPos = read2IndDict[pos]
        # we only really care if there is a difference in the sequence
        if read1Seq[read1IntPos] != read2Seq[read2IntPos]:

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
            # both show a mismatch, we only want to count one
            if read1Quals[read1IntPos] >= read2Quals[read2IntPos]:
                read2Quals[read2IntPos] = 0
            elif read1Quals[read1IntPos] < read2Quals[read2IntPos]:
                read1Seq[read1IntPos] = read2Seq[read2IntPos]
                read1Quals[read1IntPos] = 0

    # finally we have to create new reads
    read1New = read1
    read1New.query_sequence = "".join(read1Seq)
    read1New.query_qualities = read1Quals
    # same for read2
    read2New = read2
    read2New.query_sequence = "".join(read2Seq)
    read2New.query_qualities = read2Quals

    return (read1New, read2New)
