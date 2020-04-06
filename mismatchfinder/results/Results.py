from pandas import DataFrame, concat
from numpy import array


class Results(object):
    """This class contains all results accumulated after analysing a bam."""

    sample = None

    nMisMatches = None
    nSites = None
    nReads = None
    nAlignedReads = None
    nAlignedBases = None

    nSomMisMatches = None
    nSomSites = None

    nSomConfidentSites = None

    def __init__(
        self,
        sample,
        nMisMatches,
        nSites,
        nReads,
        nAlignedReads,
        nAlignedBases,
        nSomMisMatches,
        nSomSites,
        nSomConfidentSites,
        nDiscordantReads,
    ):
        super(Results, self).__init__()
        self.sample = sample
        self.nMisMatches = nMisMatches
        self.nSites = nSites
        self.nReads = nReads
        self.nAlignedReads = nAlignedReads
        self.nAlignedBases = nAlignedBases
        self.nSomMisMatches = nSomMisMatches
        self.nSomSites = nSomSites
        self.nSomConfidentSites = nSomConfidentSites
        self.nDiscordantReads = nDiscordantReads

    # just so we can actually print the object
    def __str__(self):
        return self.convertToPandasDataFrameRow().__str__()

    def convertToPandasDataFrameRow(self):
        df = DataFrame(
            data={
                "Sample": [self.sample],
                "nReads": [self.nReads],
                "nAlignedReads": [self.nAlignedReads],
                "nAlignedBases": [self.nAlignedBases],
                "nMisMatches": [self.nMisMatches],
                "nSites": [self.nSites],
                "nSomMisMatches": [self.nSomMisMatches],
                "nSomSites": [self.nSomSites],
                "nSomConfidentSites": [self.nSomConfidentSites],
                "nDiscordantReads": [self.nDiscordantReads],
            }
        )
        return df


def convertToPandasDataFrame(resultQueue):
    if resultQueue.empty():
        return None
    rowDfs = []
    # go through all results and convert them to a pandas dataframe
    while not resultQueue.empty():
        res = resultQueue.get()
        resDf = res.convertToPandasDataFrameRow()
        rowDfs.append(resDf)

    df = concat(rowDfs, ignore_index=True)
    return df
