from multiprocessing import Process
from logging import basicConfig, debug, DEBUG
from numpy import array, sort, quantile
import pysam
import datetime

#
# basicConfig(level=DEBUG, format="(%(threadName)-9s) %(message)s")
#


class BamScanner(Process):
    """Class which scans a Bam for TMB and other interesting features"""

    def __init__(
        self, bamFile, referenceFile, minMQ, minBQ, blackList=None, whiteList=None
    ):
        super(BamScanner, self).__init__()
        self.minMQ = minMQ
        self.minBQ = minBQ
        self.bamFile = pysam.AlignmentFile(
            bamFile, "r", require_index=True, reference_filename=referenceFile
        )

        self.blackList = blackList
        self.whiteList = whiteList

        basicConfig(level=DEBUG, format="%(processName)-10s  %(message)s")

    # this function gets all sites of mismatches from any mapped read in a bam
    # use qualThreshold, bedObj and minMQ to exclude reads and mismatches
    # qualThreshold is for the base quality
    # bedObj is to discard reads mapping to those blacklisted areas in the bed
    # minMQ directly ignores reads which have a mappingquality lower
    def getMutationSites(self):

        # state your purpose ;)
        debug("Checking reads for mismatches")

        # store the found sites of mutations and how often they occured
        mutSites = {}

        # store counts for stats later
        nLowQualReads = 0
        nNoMisMReads = 0
        nReads = 0
        nBlackListed = 0
        nNonWhiteListed = 0
        nMisMatches = 0
        nAlignedBases = 0

        # store the fragment lengths for later
        fragLengths = []

        totalIndexReads = self.bamFile.mapped + self.bamFile.unmapped

        # only a bam index can tell you how many reads are mapped and unmapped, a cram index doesnt
        # tell you
        if self.bamFile.is_bam:
            debug(
                f"Found {totalIndexReads} reads in the index (mapped: {bamFile.mapped}; unmapped: {bamFile.unmapped})"
            )
        # get the time we started with this
        startTime = datetime.datetime.now()

        # we work directly on the bam without iterator creation
        for read in self.bamFile:
            nReads += 1
            # report every 100 reads
            if nReads % 1000000 == 0:
                # time so far
                currTime = datetime.datetime.now()
                deltaTime = currTime - startTime
                readsPerSec = nReads / deltaTime.total_seconds()
                # we really just care about the general aread so we round to the next hundred
                readsPerSec = int(round(readsPerSec, -2))

                # because only for a bam we know how many reads there are in total, we only print
                # percentages there
                if self.bamFile.is_bam:
                    debug(
                        f"Read through {(nReads/totalIndexReads):3.2%} of reads {readsPerSec:6d} reads per second"
                    )
                else:
                    debug(
                        f"Read through {nReads} reads - processing {readsPerSec:6d} reads per second"
                    )

            # we only want proper reads and no secondaries. We can be pretty lenient here, because we
            # check for the mismatch itself if the sequencing quality is high enough later.
            if (
                read.is_duplicate
                or read.is_qcfail
                or read.is_secondary
                or read.is_unmapped
                or read.mapping_quality < self.minMQ
            ):
                nLowQualReads += 1
                continue
            else:
                # it is faster to check if there are mismatches before checking if they align to
                # either white or blacklisted regions
                if not hasMisMatches(read):
                    nNoMisMReads += 1
                else:

                    # if the read is aligned to any low complexity area, we just cant deal with it
                    if not self.blackList is None and self.blackList.isWithinRegion(
                        read
                    ):
                        nBlackListed += 1

                    else:
                        # only analyse read in the whitelist region if there is one
                        if self.whiteList is None or self.whiteList.isWithinRegion(
                            read
                        ):
                            # get all mismatches in this read
                            tmpMisMatches = scanAlignedSegment(read, self.minBQ)
                            # store the mismatches and keep a record how often each was found
                            for mm in tmpMisMatches:
                                if mm in mutSites:
                                    mutSites[mm] += 1
                                else:
                                    mutSites[mm] = 1
                            # we also store the amount of mismatches found, so we can calculate
                            # the mismatches per read which should be stable between samples
                            nMisMatches += len(tmpMisMatches)
                            # store the fragment length of this read for fragment size statistics
                            fragLengths.append(abs(read.template_length))

                            # add the amount of bases of this read that were aligned
                            nAlignedBases += read.query_alignment_length
                        else:
                            nNonWhiteListed += 1

        # we are done so we update the status as well

        debug("Read through 100.00% of reads ")

        # did we have an issue with reads?
        if len(fragLengths) == 0:
            debug(
                f"Could not detect any reads from Bam ({self.bamFile}). Further analysis is not possible"
            )
            exit(1)

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
        return {
            "sites": mutSites,
            "nLowQualReads": nLowQualReads,
            "nNoMisMatchReads": nNoMisMReads,
            "nBlackListedReads": nBlackListed,
            "nAlignedReads": nAlignedReads,
            "nReads": nReads,
            "nAlignedBases": nAlignedBases,
            "fragSizeQuantiles": fragLenQuantiles,
            "nMisMatches": nMisMatches,
            "nSites": len(mutSites),
            "misMatchesPerRead": nMisMatches / nAlignedReads,
            "nDiscordantReads": nDiscordantReads,
        }

    def run(self):
        self.getMutationSites()


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


# read through one read and find all mismatches and their contexts
# can visualize the contexts on the reads
def scanAlignedSegment(AlignedSegment, qualThreshold=21, vis=False):

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
                or referencePositions[refIndex + 1] != referencePositions[refIndex] + 1
            ):
                continue

            # now that we have a proper context defined we can go about and check what type of variant this

            # if we have more than 1 lower case char in our context, that means we ae at the
            # beginning of a multi mismatch region and we should combine the two into a multi
            # variant thing
            # so here we just move on to the next position if we find there is a mismatch in the
            # next position
            # this will also take care of trinucleotide variants
            if refContext[2:].islower():
                continue

            # however we would get the last position of a 3bp mismatch as a dinucleotide change,
            # which we do not want so we also check if there are more than 2 variants within 4
            # positions
            if (
                countLowerCase(alignedRefSequence[templatePos - 2 : templatePos + 2])
                > 2
            ):
                continue

            # that would leave us with only single or dinucleotide variants which are shifted to
            # the left e.g.
            # alt: TTG     alt: CTG
            # ref: ccG     ref: CcG
            # which we can easily convert to the classes that we want which are the C>T changes in
            # the dinucleotide context CC,CT,TC because those are specific for melanoma
            mut = (AlignedSegment.reference_name, contigPos, refContext, altContext)

            # just so we can check what kind of variants we actually put into this we can print
            # the contexts for easy readability (this does not put the context at the same
            # position if there are indels in the read but at the position it is taken from)
            if vis:
                # print the query (the actual read) on top
                printLog(AlignedSegment.query_alignment_sequence, addTime=False)
                for i in range(0, readPos - 1):
                    printLog(" ", end="", addTime=False)
                printLog(altContext, addTime=False)
                # print the template (the the refernece below)
                printLog(alignedRefSequence, addTime=False)
                for i in range(0, readPos - 1):
                    printLog(" ", end="", addTime=False)
                printLog(refContext, addTime=False)
                printLog("\n", addTime=False)

            # add to the found mutations now that everything is sorted out
            mutations.append(mut)

    return mutations


# this function counts how many lower case chars are in a string
def countLowerCase(string):
    return sum(1 for c in string if c.islower())
