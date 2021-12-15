from logging import debug, error

import matplotlib.pyplot as plt
from numpy import (
    apply_along_axis,
    argmin,
    array,
    eye,
    float64,
    full,
    hstack,
    inf,
    isinf,
    matmul,
    nonzero,
    ones,
    sum,
    zeros,
    sqrt,
    copy,
    isnan,
)
from pandas import DataFrame, read_csv
from quadprog import solve_qp
from scipy.optimize import minimize_scalar

from mismatchfinder import ext
from mismatchfinder.utils.Misc import normaliseCounts

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
                f"The counts provided ({counts.shape}) do not match the loaded signature ({self.df.shape}), you might need to load a different signature"
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

    def analyseCountsFile(self, file, method="ILM", oligoCounts=None, flatNorm=True):

        debug(f"Reading in context counts file {file}")
        with open(file, mode="r") as countFH:
            countsTable = read_csv(countFH, header=0, sep="\t", index_col=0)
            # check if we have issues
            if countsTable.isna().values.any():
                error(f"File {file} could not be parsed")
                raise Exception(f"File {file} could not be parsed")
            elif len(countsTable.index) == 0:
                error(f"File {file} did not contain any data")
                raise Exception(f"File {file} did not contain any data")
            # @TODO: make this work, because its kinda important
            # elif (countsTable.values == 0).all:
            #     error(f"File {file} did not contain any non zero counts")
            #     return DataFrame()

            if not oligoCounts is None:
                # normalise the counts, by dividing each count by the occurance in the counts
                countsTable = normaliseCounts(
                    countsDf=countsTable, contextCountDf=oligoCounts, flatNorm=flatNorm
                )

        if method == "QP":
            sigWeights = self.whichSignaturesQP(countsTable.values)
        elif method == "ILM":
            sigWeights = apply_along_axis(
                self.whichSignaturesILM, axis=1, arr=countsTable.values
            )
        else:
            error(f"Signature deconstruction method {method} not implemented")
            raise Exception(f"Signature deconstruction method {method} not implemented")

        df = DataFrame(
            data=sigWeights, index=countsTable.index, columns=self.df.columns
        )
        return df

    def whichSignaturesILM(
        self, counts, limitSigs=None, errorThreshold=1e-3, limitIterations=1000
    ):

        if len(counts) != self.nTypes:
            error(
                "The counts provided do not match the loaded signature, you might need to load a different signature"
            )
            return None

        normalisedCounts = counts / sum(counts)

        # if we didnt get a limit, we just use all of them
        if limitSigs is None:
            limitSigs = self.nSigs

        absentSigs = self.calculateAbsentSignatures(normalisedCounts)
        debug(f"Will ignore these signatures: {absentSigs}")
        # get only the indexes
        ignoreIdxs = [a["sigIdx"] for a in absentSigs]

        weights, currError = self.seedWeights(normalisedCounts, ignoreIdxs)

        iteration = 0
        errorDiff = inf
        while errorDiff > errorThreshold and iteration < limitIterations:
            iteration += 1
            debug(f"Optimisation iteration {iteration}")

            weights, newError = self.updateWeights(
                counts=normalisedCounts,
                weights=weights,
                limitSigs=limitSigs,
                ignoreSigs=ignoreIdxs,
            )
            # see how much we improved upon the last time
            if not isinf(currError):
                errorDiff = (currError - newError) / currError
                debug(f"reduced error by {errorDiff}")

            # update our state to the new weights
            currError = newError

            # if the error is 0 we might as well stop optimising
            if currError == 0:
                debug(f"Found a perfect solution")
                break

        # normalise the weights
        weights = weights / sum(weights)
        ## TODO: put in a filtering?
        # where(weights < 0.1, weights, 0)

        return weights

    def calculateAbsentSignatures(self, counts):

        # first we want to know which contexts basically never show up in the normalised counts
        # if all contexts would show equally, we would expect each to have 0.0104 so this means
        # this signature is occurs 5 times less frequent than the average
        contextsNotPresent = counts < 0.002

        # this is the place to store signatures we find are absent
        signaturesAbsent = []

        # we go through all signatures (in the columns of the dataframe)
        for i in range(self.df.shape[1]):
            contextFractions = self.df.iloc[:, i]
            # then look at their fractions affecting this contexts
            for j in range(len(contextFractions)):
                cf = contextFractions[j]
                # if this signature has a peak here
                if cf > 0.5:
                    # but the context is not present in the counts
                    if contextsNotPresent[j]:
                        # we can safely ignore the signature
                        signaturesAbsent.append(
                            {
                                "sigIdx": i,
                                "sig": self.df.columns[i],
                                "contextIdx:": j,
                                "context": self.df.index[j],
                                "cf": cf,
                            }
                        )
                    # if this peak happend, we dont need to look at any others anymore for this
                    # signature
                    break

        return signaturesAbsent

    def seedWeights(self, counts, ignoreSigs=None):
        debug(f"Seeding initial weights")
        # if we dont ignore anything, then we just transform this to an empty list
        if ignoreSigs is None:
            ignoreSigs = []

        error = inf
        weights = []
        for i in range(self.nSigs):
            if i not in ignoreSigs:
                currWeights = zeros(self.nSigs)
                currWeights[i] = 1
                currError = self.getError(counts, currWeights)
                if currError < error:
                    weights = currWeights
                    error = currError

        return weights, error

    def getError(self, counts, weights):
        reconstProf = self.reconstructProfile(weights)

        error = sum(sqrt((counts - reconstProf) ** 2))

        return error

    def reconstructProfile(self, weights):
        weights = weights / sum(weights)

        return weights.dot(self.data.T)

    def updateWeights(self, counts, weights, limitSigs, ignoreSigs):

        if ignoreSigs is None:
            ignoreSigs = []

        # see how many signatures are already used in this solution
        nUsedSigs = len([weight for weight in weights if weight != 0])

        # the error with the current weights
        currError = self.getError(counts, weights)

        # calculate which weights we are allowed to change
        # if we still have less signatures used than the limit, then we can play with any of them
        # as we only add one more signature each iteration
        if nUsedSigs < limitSigs:
            availSigs = range(self.nSigs)
            debug("adding more signatures possible")
        else:
            # if we already are at the limit, then we only adjust those which are already used
            availSigs = nonzero(weights)[0]
            debug(
                "we couldnt add any more signatures, so we need to adjust the ones we have"
            )

        # now we discard those signatures which we already know we can ignore
        availSigs = [i for i in availSigs if i not in ignoreSigs]

        # create a matrix to store optimisation results
        v = zeros((self.nSigs, self.nSigs))

        # store new errors per signature
        newErrors = full(self.nSigs, inf)

        # go through all indices that are available
        for i in availSigs:

            # for signature i find the weight that minimises the error
            def minimise(x):
                mWeights = copy(weights)
                mWeights[i] += x
                return self.getError(counts, mWeights)

            # minimise in respect to the weight > 0 and weight < 1
            errorMinimiser = minimize_scalar(
                minimise, bounds=(-weights[i], 1), method="bounded"
            ).x
            # save the result
            v[i, i] = errorMinimiser
            newWeights = weights + v[i]

            # find out how well we did
            newErrors[i] = self.getError(counts, newWeights)

        # check which of the signatures can be adjusted to minimise the error
        minNewError = min(newErrors)

        # update the weights if we find a more optimal solution
        if minNewError < currError:
            minIdx = argmin(newErrors)
            weights[minIdx] += v[minIdx, minIdx]
            currError = minNewError
        else:
            debug(f"No new improvement found")

        return weights, currError


if __name__ == "__main__":

    from logging import basicConfig

    basicConfig(
        level="DEBUG",
        format="%(asctime)s %(levelname)-8s %(processName)-10s  %(message)s",
    )

    sig = Signature.loadSignaturesFromFile(type="SBS")

    sig.analyseCountsFile(
        "/Users/hollizecksebastian/Documents/workspace/MisMatchFinder/tests/SBScounts.tsv"
    )
    # building artificial counts
    # mix = zeros(67)
    # mix[0] = mix[2] = mix[5] = mix[9] = 0.25
    #
    # # get the counts how they look like
    # counts = sig.data.dot(mix)
    #
    # print(len(counts))
    # print(counts)
    #
    # # deconstruct
    # weights = sig.whichSignaturesILM(counts, limitSigs=4)
    #
    # error = sum(sqrt((weights - mix) ** 2))
    # print(error)
    # print(weights)
