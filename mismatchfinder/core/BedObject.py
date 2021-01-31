#!/usr/bin/env python3

from ncls import NCLS
from numpy import array, int64


class BedObject(object):
    """Create a queryable object from a bed file."""

    def __init__(self, bedFile):
        super(BedObject, self).__init__()

        # This function builds an index tree from the bed file to have a fast check if a read falls within
        # a specified area or not.
        # the hard work is done by NCLS (https://github.com/biocore-ntnu/ncls) which is also used by the
        # pyranges module

        starts = []
        ends = []
        currChr = None
        self.__ncls = {}

        with open(bedFile) as f:
            for line in f:
                # break the line into fields
                lineArray = line.strip().split()
                # if the chromosome is still the same, or we do this the first time, we append
                if currChr == lineArray[0] or currChr is None:
                    # this is not changing anything but for the first time (when currChr is None), but
                    # thats fine as this is neither time consuming, nor the bottle neck, its just not
                    # pretty
                    currChr = lineArray[0]

                    # add the starts and stops to the list
                    starts.append(int(lineArray[1]))
                    ends.append(int(lineArray[2]))

                else:
                    # convert to array with dtype (ncls needs that)
                    starts = array(starts, dtype=int64)
                    ends = array(ends, dtype=int64)

                    # create the data structure (third column is ids... which could be anything, but
                    # needs to be a number )
                    tmpNcls = NCLS(starts, ends, starts)
                    # store the data structure under its chromosome name
                    self.__ncls[currChr] = tmpNcls

                    # reset all the things for the next chromosome (and initialise it while we are
                    # already at it)
                    currChr = lineArray[0]
                    starts = [lineArray[1]]
                    ends = [lineArray[2]]

    @classmethod
    def parseFile(cls, file):
        """Takes care of returning None if the file is empty or a proper bed object if thefile is a bed file"""
        if file == "" or file is None:
            return None
        else:
            return BedObject(file)

    # check if the read is aligned to any blacklisted area (again NCLS does the heavy lifting)
    def isReadWithinRegion(self, read):
        try:
            # the idea is to retrieve the NCLS object which correcponds to the contig the read mapped
            # to (read.reference_name), then ask NCLS to find an overlap between the area where the
            # read aligned to (reference_start-reference_end) and the already stored bed file
            #
            bedIt = list(
                self.__ncls[read.reference_name].find_overlap(
                    read.reference_start, read.reference_end
                )
            )
        except KeyError:
            # in this case, there was no entry on the bedfile for this contig, so we assume it is not
            # blacklisted
            return False

        # now that that worked, we just need to see if we have any overlaps at all. and if we do, we
        # can confidently say that this read is blacklisted
        if len(bedIt) > 0:
            return True
        else:
            return False

    # this is very similar to the read check, but we check again if the mismatch itself is also in
    # the region, because even though a read overlaps with a region, the mismatch could be outside
    def isMisMatchWithinRegion(self, mismatch):

        chr, pos, refContext, altContext, misMatchClass = mismatch
        # here we calculate the starting pos of the mismatch, if its a snp, the start pos is the
        # center variant, which is pos +2 -1, where if its a double, we have the actual pos +2 -2
        internalMisMatchStartPos = pos + 2 - misMatchClass
        # the end pos is always the pos after the center pos
        internalMisMatchEndPos = pos + 2

        print(
            f"startPos: {internalMisMatchStartPos} and endPos: {internalMisMatchEndPos}"
        )
        try:
            bedIt = list(
                self.__ncls[chr].find_overlap(
                    internalMisMatchStartPos, internalMisMatchEndPos
                )
            )
        except KeyError:
            return false

        # and now we decide if there was an overlap
        if len(bedIt) > 0:
            return True
        else:
            return False
