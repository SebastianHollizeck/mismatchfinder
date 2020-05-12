from logging import debug, error

import matplotlib
import matplotlib.pyplot as plt
# import seaborn as sns
from numpy import array, eye, float64, hstack, matmul, ones, sum
from pandas import DataFrame, read_csv
from quadprog import solve_qp

from mismatchfinder import ext

matplotlib.use("Agg")

try:
    import importlib.resources as pkg_resources
except ImportError:
    debug("falling back to backported importlib module")
    import importlib_resources as pkg_resources


class Signature(object):
    def __init__(self, pandas, type):
        super(Signature, self).__init__()

        self.type = type
        self.data = pandas.values
        self.df = pandas
        self.nSigs = len(self.df.columns)
        # this is not actually required, as we only use 96 different contexts, but in the future we
        # might not
        self.nTypes = len(self.df.index)

        # normalising the signature input just to be save
        self.data = self.data / self.data.sum(axis=0)[None, :]

        # create the input matrix from the reference signatures
        self.G = matmul(self.data.T, self.data)

        # set up restrictions
        # all weights need to be >0
        # and all weights need to sum up to 1
        self.C = hstack(
            (ones((self.nSigs, 1), dtype=float64), eye(self.nSigs, dtype=float64))
        )
        self.b = array([1.0] + [0.0] * self.nSigs, dtype=float64)

    def whichSignatures(self, counts):

        if counts.shape[1] != self.nTypes:
            error(
                "The counts provided do not match the loaded signature, you might need to load a different signature"
            )
            return None

        normalisedCounts = counts / counts.sum(axis=1)[:, None]

        # create the input vector of observed counts
        D = matmul(normalisedCounts, self.data)
        debug("Estimating contributions of signatures to counts")
        weights = array([solve_qp(self.G, d, self.C, self.b, meq=1)[0] for d in D])

        # floats are a bit wonky sometimes, so we check that everything is still > 0 and renormalise
        weights[weights < 0] = 0
        weights = weights / weights.sum(axis=1)[:, None]
        # # TODO: Add "unexplained" percentage
        return weights

    @classmethod
    def loadSBSSignaturesFromFile(cls, file=None):
        if file is None:
            debug("Using default SBS signatures")
            with pkg_resources.path(ext, "sigProfiler_SBS_signatures.csv") as path:
                file = path


        debug(f"Loading reference signature file {file}")
        with open(file, "r") as sigFH:

            prelimSigs = read_csv(sigFH, header=0, sep=",", index_col=0)

            # and here we just discard them
            return Signature(prelimSigs, "SBS")

    @classmethod
    def loadDBSSignaturesFromFile(cls, file=None):

        if file is None:
            debug("Using default DBS signatures")
            with pkg_resources.path(ext, "sigProfiler_DBS_signatures.csv") as path:
                file = path

        debug(f"Loading reference signature file {file}")
        with open(file, "r") as sigFH:

            prelimSigs = read_csv(sigFH, header=0, sep=",", index_col=0)

            # and here we just discard them
            return Signature(prelimSigs, "DBS")

    def analyseCountsFile(self, file):

        debug(f"Reading in context counts file {file}")
        with open(file, mode="r") as countFH:
            countsTable = read_csv(countFH, header=0, sep="\t", index_col=0)

        sigWeights = self.whichSignatures(countsTable.values)

        df = DataFrame(
            data=sigWeights, index=countsTable.index, columns=self.df.columns
        )
        return df

#
# def plotSignatures(weights):
#
#     items = []
#     sigNames = weights.columns
#     for bam, r in weights.iterrows():
#         for i, sig in enumerate(weights.columns):
#             items.append({"bam": bam, "signature": sigNames[i], "weight": r[sig]})
#
#     df = DataFrame(items)
#
#     print(df)
#
#     # Create plots and output to file
#     sns.set_style("whitegrid")
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
#     sns.boxplot(x="signature", y="weight", data=df, ax=ax1)
#     sns.boxplot(x="signature", y="weight", data=df, ax=ax2)
#     plt.suptitle("testplot", y=1)
#
#     plt.tight_layout()
#     plt.savefig("test.png")
