import unittest
import pathlib

from mismatchfinder.core.Signatures import Signature
import numpy as np


def calculateError(one, two):
    return sum(np.sqrt((one - two) ** 2))


class sigDecompTest(unittest.TestCase):
    def test_loadDefaultSBSSignature(self):
        self.assertIsInstance(Signature.loadSignaturesFromFile(type="SBS"), Signature)

    def test_loadDefaultDBSSignature(self):
        self.assertIsInstance(Signature.loadSignaturesFromFile(type="DBS"), Signature)

    def test_analyseCountFile(self):
        countFile = pathlib.Path(__file__).parent / "SBScounts.tsv"
        sig = Signature.loadSignaturesFromFile(type="SBS")
        weights = sig.analyseCountsFile(countFile).values.flatten()

        # building the truth
        mix = np.zeros(67)
        mix[0] = mix[2] = mix[5] = mix[9] = 0.25

        # compare true weights with deconstructed weights
        error = calculateError(weights, mix)

        # error should be very small
        self.assertTrue(error < 0.0001)

    def test_analyseCountFileFail(self):
        countFile = pathlib.Path(__file__).parent / "falseStructure.tsv"
        sig = Signature.loadSignaturesFromFile(type="SBS")
        try:
            weights = sig.analyseCountsFile(countFile)
        except Exception as e:
            self.assertEqual(str(e), f"File {countFile} could not be parsed")

    def test_decompILM(self):
        sig = Signature.loadSignaturesFromFile(type="SBS")

        # building artificial counts
        mix = np.zeros(67)
        mix[0] = mix[2] = mix[5] = mix[9] = 0.25

        # get the counts how they look like
        counts = sig.data.dot(mix) * 100000

        weights = sig.whichSignaturesILM(counts)
        # compare true weights with deconstructed weights
        error = calculateError(weights, mix)

        # error should be very small
        self.assertTrue(error < 0.0001)

    def test_decompQP(self):
        sig = Signature.loadSignaturesFromFile(type="SBS")

        # building artificial counts
        mix = np.zeros(67)
        mix[0] = mix[2] = mix[5] = mix[9] = 0.25

        # get the counts how they look like
        counts = sig.data.dot(mix) * 100000

        # QP actually wants a matrix
        counts = counts.reshape((1, 96))

        # deconstruct
        weights = sig.whichSignaturesQP(counts)

        weights = weights.flatten()

        # compare true weights with deconstructed weights
        error = calculateError(weights, mix)

        # error should be very small
        self.assertTrue(error < 0.0001)


if __name__ == "__main__":
    unittest.main()
