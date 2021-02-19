import datetime
from logging import debug, error, info, basicConfig, getLogger
from multiprocessing import Process
from sys import getsizeof
from collections import defaultdict

import pysam
from numpy import array, quantile, sort

from mismatchfinder.results.Results import MismatchCandidates
from mismatchfinder.core.EndMotives import EndMotives
from mismatchfinder.utils.Misc import countLowerCase
from mismatchfinder.core.Fragment import Fragment


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
        nAlignedReads = 0
        nDiscordantReads = 0

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
        read_dict = defaultdict(lambda: [None, None])

        for read in bamFile.fetch(until_eof=True):
            nReads += 1
            # we check every 10K reads how uch time has passed and if it has been more than 30
            # we give an update
            if nReads % 1000 == 0:

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

            # we build proper fragments of DNA from the reads, if they are paired in sequencing
            # or we discard them due to quality issues, but we still discard any supplemenary or
            # secondary reads
            if read.is_duplicate or read.is_secondary or read.is_supplementary:
                nLowQualReads += 1
                continue

            # store the read until we find the second one, but only if they are aligned to the
            # chromosome, otherwise we just continue
            qname = read.query_name
            if read.next_reference_id == read.reference_id:
                if qname not in read_dict:
                    read_dict[qname] = read
                    continue
                else:
                    # now that we have access to both reads, we can decide, how we build this
                    # fragment
                    mate = read_dict[qname]
                    del read_dict[qname]
                    r1U = isReadUsable(read, minMQ=self.minMQ)
                    r2U = isReadUsable(mate, minMQ=self.minMQ)
                    if r1U and r2U:
                        frag = Fragment(read, mate)
                    elif r1U:
                        frag = Fragment(read)
                        nLowQualReads += 1
                    elif r2U:
                        frag = Fragment(mate)
                        nLowQualReads += 1
                    else:
                        # both of the reads are bad, so we just continue
                        nLowQualReads += 2
                        continue

            else:
                # if this is single end sequencing, we just build the fragment right away
                frag = Fragment(read)
                # but if it is paired, then this is a discordant read
                if read.is_paired:
                    nDiscordantReads += 1

            # now we only need to check if there is a blacklist or whitelist, which restricts us further
            if not self.blackList is None and self.blackList.isFragmentInRegion(frag):
                if fragment.is_paired:
                    nBlackListedReads += 2
                else:
                    nBlackListedReads += 1

            else:
                # analyse if the fragment overlaps with the whitelist
                if self.whiteList is None or self.whiteList.isFragmentInRegion(frag):
                    # get the mutations the fragment contains
                    tmpMM = frag.getMismatches(self.minBQ)

                    # check if the actual mismatch is in the white/blacklist
                    for mm in tmpMM:
                        if (
                            not self.blackList is None
                            and self.blackList.isMisMatchWithinRegion(mm)
                        ):
                            continue
                        if (
                            self.whiteList is None
                            or self.whiteList.isMisMatchWithinRegion(mm)
                        ):
                            # increase the count for this mismatch if we found it already
                            if mm in mutSites:
                                mutSites[mm] += 1
                            else:
                                mutSites[mm] = 1
                            nMisMatches += 1

                    nAlignedBases += len(frag.querySeq)
                    if frag.is_paired:
                        nAlignedReads += 2
                    else:
                        nAlignedReads += 1

                    # TODO: add endmotif analysis to the fragments

                    fragLengths.append(frag.fragment_end - frag.fragment_start)
                    del frag
                else:
                    nBlackListedReads += 1

        self.logger.info(
            f"Read through 100.00% of reads in {(datetime.datetime.now()-startTime).total_seconds()/60:.1f} minutes"
        )

        # did we have an issue with reads?
        if nAlignedBases == 0:
            self.logger.error(
                f"Could not detect any reads from Bam ({self.bamFilePath}). Further analysis is not possible\nLowQualReads: {nLowQualReads}\nNoMisMatchReads: {nNoMisMatchReads}\nBlackListedReads: {nBlackListedReads}"
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


def isReadUsable(read, minMQ):

    return (
        read.mapping_quality >= minMQ
        and not read.is_qcfail
        and not read.is_unmapped
        and not read.is_duplicate
    )
