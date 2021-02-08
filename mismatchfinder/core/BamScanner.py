import datetime
from logging import debug, error, info, basicConfig, getLogger
from multiprocessing import Process
from sys import getsizeof

import pysam
from numpy import array, quantile, sort

from mismatchfinder.results.Results import MismatchCandidates
from mismatchfinder.core.EndMotives import EndMotives
from mismatchfinder.utils.Misc import countLowerCase


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

        # store the fragment lengths for later
        fragLengths = []

        # get the time we started with this
        startTime = datetime.datetime.now()
        currTime = datetime.datetime.now()

        # we work directly on the bam without iterator creation
        bamFile = pysam.AlignmentFile(
            self.bamFilePath, "r", reference_filename=self.referenceFile
        )

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

            # we only want proper reads and no secondaries. We can be pretty lenient here, because we
            # check for the mismatch itself if the sequencing quality is high enough later.
            if (
                read.is_duplicate
                or read.is_qcfail
                or read.is_secondary
                or read.is_supplementary
                or read.is_unmapped
                or read.mapping_quality < self.minMQ
            ):
                nLowQualReads += 1
                continue
            else:

                # if we did get a blacklist and the read actually is within that region, we
                # discard it and keep a count on how many we discarded
                if not self.blackList is None and self.blackList.isWithinRegion(read):
                    nBlackListedReads += 1

                else:
                    # only analyse read in the whitelist region if there is one
                    if self.whiteList is None or self.whiteList.isReadWithinRegion(
                        read
                    ):
                        # we only do our mismacth analysis if the read actually does have mismatches
                        if not hasMisMatches(read):
                            nNoMisMatchReads += 1
                        else:
                            # get all mismatches in this read
                            tmpMisMatches = self.scanAlignedSegment(read, self.minBQ)
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
                        if read.is_read1:
                            fragLengths.append(abs(read.template_length))

                        # add the amount of bases of this read that were aligned
                        nAlignedBases += read.query_alignment_length

                        # we also care about the endmotives of the reads, so we store those
                        self.endMotives.count(read)
                    else:
                        nBlackListedReads += 1

        # we are done so we update the status as well

        self.logger.info(
            f"Read through 100.00% of reads in {(datetime.datetime.now()-startTime).total_seconds()/60:.1f} minutes"
        )

        # did we have an issue with reads?
        if len(fragLengths) == 0:
            self.logger.error(
                f"Could not detect any reads from Bam ({self.bamFilePath}). Further analysis is not possible\nLowQualReads: {nLowQualReads}\nNoMisMatchReads: {nNoMisMatchReads}\nBlacklistedReads: {nBlackListedReads}"
            )
            raise Exception(f"No reads left for {self.bamFilePath}")

        # we use fragLengths here, because it already contains all aligned reads
        nAlignedReads = len(fragLengths)

        # sort the fragment lengths so we can get the median and the sides of the distribution for
        # estimation of fusions or similar (because those would have higher fragment lengths as
        # expected )
        fragLengths = array(sort(fragLengths))

        # we get the middle of the array, which should be the median of the distribution
        median = fragLengths[int(nAlignedReads / 2)]

        # find the index, at which the read distance gets to normal amounts
        lowSizeIdx = -1
        for i in range(0, nAlignedReads):
            if fragLengths[i] < 35:
                lowSizeIdx = i
            else:
                break

        # find the index, at which the read distance gets to normal amounts
        highSizeIdx = -1
        for i in range(nAlignedReads - 1, 0, -1):
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
    def scanAlignedSegment(self, AlignedSegment, qualThreshold=21, vis=False):

        # this is where we will store the sites we find in this read
        mutations = []

        # we do this here, so we only do it once instead of in the loop
        alignedRefSequence = AlignedSegment.get_reference_sequence()
        referencePositions = AlignedSegment.get_reference_positions(full_length=True)

        # loop through the read
        for (readPos, contigPos, seq) in AlignedSegment.get_aligned_pairs(
            with_seq=True, matches_only=True
        ):

            # if it is lower case it symbolises a mismatch
            if seq.islower():
                # offset the soft cliping at the beginning
                readPos = readPos - AlignedSegment.query_alignment_start
                # offset the reference location
                templatePos = contigPos - AlignedSegment.reference_start

                # we really only want high quality mismatches
                qual = AlignedSegment.query_qualities[readPos]
                if qual < qualThreshold:
                    continue

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
