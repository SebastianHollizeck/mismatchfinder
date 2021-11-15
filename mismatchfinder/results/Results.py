from copy import deepcopy
from logging import debug, info

from pandas import DataFrame, concat

from mismatchfinder.core.GermlineObject import ChromosomeCache
from mismatchfinder.utils.Misc import (
    DBSrefs,
    convertToDBSTable,
    convertToSBSTable,
    countLowerCase,
    reverseComplement,
)


def convertToPandasDataFrame(resultQueue):
    if resultQueue.empty():
        return None
    rowDfs = []
    # go through all results and convert them to a pandas dataframe
    while not resultQueue.empty():
        res = resultQueue.get()
        rowDfs.append(res)

    df = concat(rowDfs)
    return df


class MismatchCandidates(object):
    """This class contains all the information of mismatched found in a sample"""

    def __init__(
        self,
        mutSites,
        nReads,
        nNoMisMatchReads,
        nBlackListedReads,
        nLowQualReads,
        nAnalysedBases,
        fragmentSizeQuantiles,
        nMisMatches,
        nDiscordantReads,
        nAlignedBases,
        bam,
    ):

        super(MismatchCandidates, self).__init__(),
        self.bam = bam
        self.mutSites = mutSites
        self.nReads = nReads
        self.nNoMisMatchReads = nNoMisMatchReads
        self.nBlackListedReads = nBlackListedReads
        self.nLowQualReads = nLowQualReads
        self.nAlignedReads = nAlignedReads
        self.fragmentSizeQuantiles = fragmentSizeQuantiles
        self.nMisMatches = nMisMatches
        self.nSites = len(self.mutSites)
        self.nDiscordantReads = nDiscordantReads
        self.nAnalysedBases = nAnalysedBases

    def convertToPandasDataFrameRow(self):
        fields = deepcopy(self.__dict__)
        # remove fields which have objects or similar as value
        try:
            fields.pop("mutSites")
        except KeyError:
            # this is acceptable because they are actually removed for pickling the results
            debug("mitSites field was already deleted")
        try:
            fields.pop("fragmentSizeQuantiles")
        except KeyError:
            # this is acceptable because they are actually removed for pickling the results
            debug("fragmentSizeQuantiles field was already deleted")
        try:
            fields.pop("SBScontexts")
        except KeyError:
            # this is acceptable because they are actually removed for pickling the results
            debug("SBScontexts field was already deleted")
        try:
            fields.pop("DBScontexts")
        except KeyError:
            # this is acceptable because they are actually removed for pickling the results
            debug("DBScontexts field was already deleted")
        try:
            fields.pop("chrHash")
        except KeyError:
            # this is acceptable because they are actually removed for pickling the results
            debug("chrHash field was already deleted")
        try:
            fields.pop("contextHash")
        except KeyError:
            # this is acceptable because they are actually removed for pickling the results
            debug("contextHash field was already deleted")

        # use the rest to make a dataframe
        df = DataFrame.from_records([fields], index=[self.bam])
        return df

    # just so we can actually print the object
    def __str__(self):
        return self.convertToPandasDataFrameRow().__str__()

    # this function aggregates the found mismatches and their contexts into counts
    # please use a SORTED list of mismatches, to avoid unnecessary caching when using germline check
    # without the check it doesnt really matter
    def countContexts(self):
        # just tell me what you are doing
        info("Counting context occurrences in all mismatches")

        resSBS = {}
        resDBS = {}

        for (chr, pos, refContext, altContext, misMatchClass) in self.mutSites:

            # get the key before we play around with the contexts
            occurrences = self.mutSites[
                (chr, pos, refContext, altContext, misMatchClass)
            ]

            # this means we have only one mismatch
            if misMatchClass == 1:
                # we need to convert the mismatch into the right strand (reverse complement, if the
                # mutated base is not either a C or a T)
                if not refContext[1] == "c" and not refContext[1] == "t":
                    refContext = reverseComplement(refContext)
                    altContext = reverseComplement(altContext)

                # now that the variants are normalised, we convert them into counts
                key = refContext.upper() + ">" + altContext[1]
                if not key in resSBS:
                    resSBS[key] = occurrences
                else:
                    resSBS[key] += occurrences
            elif misMatchClass == 2:
                # DBS can only happen in the first two bases, because we shifted the mismatches
                # so we can discard the last char of the context and convert it to upper, as
                # we know that the are the changes
                refContext = refContext[:2].upper()
                altContext = altContext[:2].upper()
                # there is only a limited amount of ref doublets, that cannot be explained through
                # their reverse complement so we look for those and count them

                if not refContext in DBSrefs:
                    refContext = reverseComplement(refContext)
                    altContext = reverseComplement(altContext)

                if not refContext in DBSrefs:
                    debug(
                        f"Contexts ref:{reverseComplement(refContext)} alt:{reverseComplement(altContext)} not found in DBS context database... skipping"
                    )
                    continue

                key = refContext + ">" + altContext

                if not key in resDBS:
                    resDBS[key] = occurrences
                else:
                    resDBS[key] += occurrences

        # need to have a newline for that as well
        info("Checked 100.00% of mismatches               ")
        self.SBScontexts = convertToSBSTable(resSBS)
        self.DBScontexts = convertToDBSTable(resDBS)
        debug("Done counting contexts")

    def checkGermlineStatus(self, germlineObj, discard=False, confidenceThreshold=2):

        if not germlineObj is None:
            info(f"Checking germline status for {len(self.mutSites)} sites")
            self.nSomaticMisMatches = 0
            self.nSomaticMisMatchSites = 0
            self.nConfidentSomaticMisMatchSites = 0
        else:
            info("No germline resource found, so no check will be performed")
            self.nSomaticMisMatches = None
            self.nSomaticMisMatchSites = None
            self.nConfidentSomaticMisMatchSites = None
            return

        nShiftFail = 0

        res = {}

        cache = None
        # we go through all sites, for caching purposes its important here to have a sorted list
        # otherwise we have to refresh the cache with every change of the chromosome
        for (chr, pos, refContext, altContext, misMatchClass) in sorted(self.mutSites):

            # first we need to cache the data from the chromosome that the site is on
            if cache is None or cache.chr != chr:
                cache = ChromosomeCache(germlineObj, chr)
                # if there is no cache for this chromosome we just get none back
                if cache is None:
                    # in this case we cannot check the germline status and just assume everything is
                    # somatic, but we can discard sites on chromosomes when we dont know the status
                    if discard:
                        continue

            occurrences = self.mutSites[
                (chr, pos, refContext, altContext, misMatchClass)
            ]
            # get the index of all germline variants with a trinucleotide at this position
            varIdxs = cache.findOverlaps(pos)

            # we do not really need the start and stop here, we could just use the idx, but
            # NCLS gives the full info anyways
            for (start, stop, idx) in varIdxs:
                # this is the position in the contexts we are talking about (again pysam is 0
                # based which makes this part a tad easier)
                offset = start - pos

                # we can ignore this position if it is upper case in the ref, because that would
                # mean we didnt find a variant there
                if refContext[offset].isupper():
                    continue

                # if the ref doesnt have length 1 we are dealing with an indel, and we dont care
                ref = cache.refs[idx]
                if len(ref) != 1:
                    continue

                # we can have several alts, so we check all of them
                for alt in cache.alts[idx]:
                    # again if this is not of length 1 we have an indel and we do not really
                    # want to deal with it
                    if len(alt) != 1:
                        continue
                    # now here we need to check if the alt matches the alt that we found in the
                    # reads
                    if altContext[offset].upper() == alt:
                        # this means this is a germline variant and we just make it an upper case
                        # to keep the variant but show its not somatic
                        refContext = (
                            f"{refContext[:offset]}{ref}{refContext[offset+1:]}"
                        )
                        # we can also stop iterating through the alts, because only one ref and alt
                        # combination will match with the observed variant
                        break
            # with us changing around things for germline variants, we need to discard the
            # cases, that now dont have variants anymore (at least no somatic ones)
            misMatchClass = countLowerCase(refContext)
            if misMatchClass == 0:
                continue
            else:
                if refContext[1].isupper():
                    nShiftFail += 1

            res[(chr, pos, refContext, altContext, misMatchClass)] = occurrences
            # we want to count the double mismatches as double the single mismatches
            # at least at the moment
            self.nSomaticMisMatches += misMatchClass * occurrences
            self.nSomaticMisMatchSites += 1
            if occurrences >= confidenceThreshold:
                self.nConfidentSomaticMisMatchSites += 1

        info(f"Found {nShiftFail} trinucleotides without a center variant")
        info(f"Found {self.nSomaticMisMatchSites} somatic sites")
        self.mutSites = res

        # we really need to clean up after ourselves because memory foot prints are horrendous
        del cache

    def writeSBSToFile(self, outFileRoot, bamFilePath):

        file = outFileRoot.parent / (outFileRoot.name + "_SBScontexts.tsv")

        # we append this, as we do not want to overwrite if there are several threads writing to it
        debug(f"Writing SBS to file: {file}")
        with open(file, "a") as outFH:
            # because we only want one row, we need newline instead of delimiter
            joined = "\t".join(str(x) for x in self.SBScontexts)
            outFH.write(f"{bamFilePath.name}\t{joined}\n")

    def writeDBSToFile(self, outFileRoot, bamFilePath):

        file = outFileRoot.parent / (outFileRoot.name + "_DBScontexts.tsv")

        # we append this, as we do not want to overwrite if there are several threads writing to it
        debug(f"Writing DBS to file: {file}")
        with open(file, "a") as outFH:
            # because we only want one row, we need newline instead of delimiter
            joined = "\t".join(str(x) for x in self.DBScontexts)
            outFH.write(f"{bamFilePath.name}\t{joined}\n")

    def writeStatsToFile(self, outFileRoot, bamFilePath):

        file = outFileRoot.parent / (outFileRoot.name + "_stats.tsv")
        stats = self.convertToPandasDataFrameRow()
        # we append this, as we do not want to overwrite if there are several threads writing to it
        debug(f"Writing stats to file: {file}")
        with open(file, "a") as outFH:

            # if the position in the file is still 0 we need to write the header, otherwise we just
            # add the line
            if outFH.tell() == 0:
                debug(f"Writing header to file {file}")
                stats.to_csv(outFH, sep="\t", index=False)
            else:
                stats.to_csv(outFH, sep="\t", header=False, index=False)

    def cleanUpForPickle(self):
        fields = self.__dict__

        if "mutSites" in fields:
            del self.mutSites
        if "fragmentSizeQuantiles" in fields:
            del self.fragmentSizeQuantiles
        if "SBScontexts" in fields:
            del self.SBScontexts
        if "DBScontexts" in fields:
            del self.DBScontexts

    def writeSitesToFile(self, outFileRoot, bamFilePath):
        # we put the bamname into the file name
        file = outFileRoot.parent / (f"{outFileRoot.name}_{bamFilePath.name}_sites.tsv")
        debug(f"Spilling all sites to file: {file}")
        with open(file, "w") as outFH:
            for (chr, pos, refContext, altContext, misMatchClass) in self.mutSites:
                count = self.mutSites[(chr, pos, refContext, altContext, misMatchClass)]
                outFH.write(
                    f"{chr}\t{pos}\t{refContext}\t{altContext}\t{count}\t{misMatchClass}\n"
                )
