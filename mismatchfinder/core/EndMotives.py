from itertools import product
from pandas import DataFrame

from mismatchfinder.utils.Misc import reverseComplement, countLowerCase

from logging import debug


class EndMotives(object):
    def __init__(self, kmer=4):

        super(EndMotives, self).__init__()

        self.kmer = kmer
        self.counts = EndMotives.createEndMotivesDict(kmer)
        self.n = 0
        self.freqs = None

    def count(self, read):

        endSeq = EndMotives.get5primeReadEnd(read, self.kmer)

        try:
            self.counts[endSeq] = self.counts[endSeq] + 1
        except KeyError:
            debug(f"Unknown read end: {endSeq}")
            self.counts["unknown"] = self.counts["unknown"] + 1

        self.n = self.n + 1

    @classmethod
    def createEndMotivesDict(cls, kmer=4):
        debug(f"Creating possible end motives of length {kmer}")
        bases = ["A", "T", "C", "G"]

        motiveCombs = ["".join(p) for p in product(bases, bases, repeat=2)]

        motiveCountDict = dict.fromkeys(motiveCombs, 0)

        debug(f"Created {len(motiveCountDict)} different end motives")

        # we also add a count for how often we found an end that we did not expect
        motiveCountDict["unknown"] = 0
        return motiveCountDict

    @classmethod
    def get5primeReadEnd(cls, read, length):

        # we only need the reference sequence to get the sequence
        alignedQuerySequence = read.query_alignment_sequence

        # we need the downstream end for if its a forward read, and the upstream one for the
        # reverse read
        if read.is_reverse:
            endSequence = alignedQuerySequence[len(alignedQuerySequence) - length :]
            endSequence = reverseComplement(endSequence)
        else:
            endSequence = alignedQuerySequence[0:length]

        return endSequence

    def getFrequencies(self, refresh=False):

        debug("Converting counts to frequencies")
        if self.freqs is None or refresh:
            # we create a frequency dictionary by dividing every value by the total amount of ends
            self.freqs = {k: v / self.n for k, v in self.counts.items()}

        return self.freqs

    def writeFreqsToFile(self, outFileRoot, bamFilePath):

        file = outFileRoot.parent / (outFileRoot.name + "_endMotives.tsv")

        with open(file, "a") as outFH:

            df = DataFrame([self.getFrequencies()])

            # add the bam column at the beginning
            df.insert(0, "bam", bamFilePath)

            # if the position in the file is still 0 we need to write the header, otherwise we just
            # add the line
            if outFH.tell() == 0:
                debug(f"Writing header to file {file}")
                df.to_csv(outFH, sep="\t", index=False)
            else:
                df.to_csv(outFH, sep="\t", header=False, index=False)
