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

        (downStreamEnd, upStreamEnd) = EndMotives.get5primeFragmentEnds(read, self.kmer)

        self.counts[downStreamEnd] = self.counts[downStreamEnd] +1
        self.counts[upStreamEnd] = self.counts[upStreamEnd] +1

        self.n = self.n+2



    @classmethod
    def createEndMotivesDict(cls, kmer=4):
        debug(f"Creating possible end motives of length {kmer}")
        bases=["A", "T", "C", "G"]

        motiveCombs = [ "".join(p) for p in product(bases, bases, repeat=2)]

        motiveCountDict = dict.fromkeys(motiveCombs, 0)

        debug(f"Created {len(motiveCountDict)} different end motives")

        #we also add a count for how often we had a mismatch in the ends
        motiveCountDict["mismatch"] = 0
        return(motiveCountDict)

    @classmethod
    def get5primeFragmentEnds(cls, read, length):

        #we only need the reference sequence to get the sequence
        alignedRefSequence = read.get_reference_sequence()

        downStream = alignedRefSequence[0:length]
        upStream = alignedRefSequence[len(alignedRefSequence)-length:]

        #we do want the reverse complement of the downstream sequence to get the 5' of both ends
        downStream = reverseComplement(downStream)

        #now we check if there is a mismatch in the end sequences and set the and accordingly
        if(countLowerCase(downStream) != 0):
            downStream = "mismatch"
        if(countLowerCase(upStream) != 0):
            upStream ="mismatch"

        return (downStream, upStream)


    def getFrequencies(self, refresh=False):

        debug("Converting counts to frequencies")
        if self.freqs is None or refresh:
            #we create a frequency dictionary by dividing every value by the total amount of ends
            self.freqs = {k: v / self.n for k, v in self.counts.items()}

        return(self.freqs)

    def writeFreqsToFile(self, outFileRoot, bamFilePath):

        file = outFileRoot.parent / (outFileRoot.name + "_endMotives.tsv")

        with open(file, "a") as outFH:

            df = DataFrame([self.getFrequencies()])

            #add the bam column at the beginning
            df.insert(0, "bam", bamFilePath)

            # if the position in the file is still 0 we need to write the header, otherwise we just
            # add the line
            if outFH.tell() == 0:
                debug(f"Writing header to file {file}")
                df.to_csv(outFH, sep="\t", index=False)
            else:
                df.to_csv(outFH, sep="\t", header=False, index=False)
