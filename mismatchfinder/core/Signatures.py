from logging import debug, error

import matplotlib.pyplot as plt
from numpy import array, eye, float64, hstack, matmul, ones, sum
from pandas import DataFrame, read_csv
from quadprog import solve_qp

from mismatchfinder import ext

try:
    import importlib.resources as pkg_resources
except ImportError:
    debug("falling back to backported importlib module")
    import importlib_resources as pkg_resources


class Signature(object):

    def __init__(self, pandas, type):
        """
        Note: dont use this constructor, but instead use the class methods 'loadSignaturesFromFile'
        to ensure this is initialised the right way.
        This class stores the signatures and allows the decontruction of counts

        :param pandas the signature dataframe in pandas notations
        :param type flag to store if this is an SBS or DBS signature
        """

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

    def whichSignaturesQP(self, counts):
        """
        Deconstucts the counts for the trinucleotide contexts into the weight of each of the signatures saved in this object
        This method is based on the math from the paper: 'Decomposition of mutational context signatures using quadratic programming methods'

        :param counts The numpy array of counts with each instance being a row and each context being a column
        """

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
    def loadSignaturesFromFile(cls, file=None, type="SBS"):
        """
        Method to create a signature object from a file.
        There are some default signatures available for usage if no file is specified
        :param file the file to create the signatures from (csv with signatures in the columns and the trinucleotides in the rows with header as well as rownames) will load default signatures if no file given
        :param type sets the type of the signatures and is responsible to decide which decault signatures to load if no file is given (SBS/DBS)
        """

        if file is None:
            if type != "SBS" and type != "DBS":
                error(f"No signature associated with default value: {type}")
                exit(1)

            debug(f"Using default {type} signatures")
            with pkg_resources.path(ext, f"sigProfiler_{type}_signatures.csv") as path:
                file = path

        debug(f"Loading reference signature file {file}")
        with open(file, "r") as sigFH:

            prelimSigs = read_csv(sigFH, header=0, sep=",", index_col=0)

            # and here we just discard them
            return Signature(prelimSigs, type)


    def analyseCountsFile(self, file):

        debug(f"Reading in context counts file {file}")
        with open(file, mode="r") as countFH:
            countsTable = read_csv(countFH, header=0, sep="\t", index_col=0)

        sigWeights = self.whichSignaturesQP(countsTable.values)

        df = DataFrame(
            data=sigWeights, index=countsTable.index, columns=self.df.columns
        )
        return df
