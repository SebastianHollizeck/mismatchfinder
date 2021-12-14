from itertools import product
from logging import debug
from collections import defaultdict

from ncls import NCLS
from numpy import arange, array, int64, sum, empty, divide
from pysam import FastaFile
from pyranges import read_bed, from_dict


# this function counts how many lower case chars are in a string
def countLowerCase(string):
    return sum(1 for c in string if c.islower())


# this function builds an ncls index from chromosomal sites for fast query for certain locations
def buildNCLSindex(sites):

    starts = array(sites, dtype=int64)
    ends = array(starts + 1, dtype=int64)
    idxs = arange(len(starts))

    index = NCLS(starts, ends, idxs)
    return index


# the complement for each normal base to be used in reverse complement
complement = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
}

# returns the reverse complement of the bases put in (will fail for non standard bases)
def reverseComplement(dna):
    bases = list(dna)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = "".join(bases)
    return bases


# generate the order of the counts to work with the order from cosmic we need all possible
# of acgt first three times always with a c in the middle and then three times always with a T in
# the middle to have all possibilities first C>A then C>G then C>T
SBSorder = [a + "C" + b + ">A" for a, b in product(["A", "C", "G", "T"], repeat=2)]
SBSorder += [a + "C" + b + ">G" for a, b in product(["A", "C", "G", "T"], repeat=2)]
SBSorder += [a + "C" + b + ">T" for a, b in product(["A", "C", "G", "T"], repeat=2)]
SBSorder += [a + "T" + b + ">A" for a, b in product(["A", "C", "G", "T"], repeat=2)]
SBSorder += [a + "T" + b + ">C" for a, b in product(["A", "C", "G", "T"], repeat=2)]
SBSorder += [a + "T" + b + ">G" for a, b in product(["A", "C", "G", "T"], repeat=2)]


def convertToSBSTable(countDict):
    debug(f"Converting dictionary with counts to fixed array in SANGER style for SBS")
    res = []
    for entry in SBSorder:
        if entry in countDict:
            res.append(countDict[entry])
        else:
            res.append(0)
    debug("Finished with SBS")
    res = array(res)
    return res


DBSrefs = set(["AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT"])
DBSorder = [
    "AC>CA",
    "AC>CG",
    "AC>CT",
    "AC>GA",
    "AC>GG",
    "AC>GT",
    "AC>TA",
    "AC>TG",
    "AC>TT",
    "AT>CA",
    "AT>CC",
    "AT>CG",
    "AT>GA",
    "AT>GC",
    "AT>TA",
    "CC>AA",
    "CC>AG",
    "CC>AT",
    "CC>GA",
    "CC>GG",
    "CC>GT",
    "CC>TA",
    "CC>TG",
    "CC>TT",
    "CG>AT",
    "CG>GC",
    "CG>GT",
    "CG>TA",
    "CG>TC",
    "CG>TT",
    "CT>AA",
    "CT>AC",
    "CT>AG",
    "CT>GA",
    "CT>GC",
    "CT>GG",
    "CT>TA",
    "CT>TC",
    "CT>TG",
    "GC>AA",
    "GC>AG",
    "GC>AT",
    "GC>CA",
    "GC>CG",
    "GC>TA",
    "TA>AT",
    "TA>CG",
    "TA>CT",
    "TA>GC",
    "TA>GG",
    "TA>GT",
    "TC>AA",
    "TC>AG",
    "TC>AT",
    "TC>CA",
    "TC>CG",
    "TC>CT",
    "TC>GA",
    "TC>GG",
    "TC>GT",
    "TG>AA",
    "TG>AC",
    "TG>AT",
    "TG>CA",
    "TG>CC",
    "TG>CT",
    "TG>GA",
    "TG>GC",
    "TG>GT",
    "TT>AA",
    "TT>AC",
    "TT>AG",
    "TT>CA",
    "TT>CC",
    "TT>CG",
    "TT>GA",
    "TT>GC",
    "TT>GG",
]


def convertToDBSTable(countDict):
    debug(f"Converting dictionary with counts to fixed array in SANGER style for DBS")
    res = []
    for entry in DBSorder:
        if entry in countDict:
            res.append(countDict[entry])
        else:
            res.append(0)
    debug(f"Finished with DBS")
    res = array(res)
    return res


def countContexts(fastaFilePath, whiteListBed=None, blackListBed=None):
    debug(f"Starting to count contexts of nucleotides in {fastaFilePath}")

    triNucCounts = defaultdict(int)
    diNucCounts = defaultdict(int)
    # open the fastaFile
    with FastaFile(fastaFilePath) as fastaFile:

        # if we do not have a whitelist to start out, we make one from the fasta, which includes
        # everything
        if whiteListBed is None:
            wlObj = from_dict(
                {
                    "Chromosome": fastaFile.references,
                    "Start": [1] * fastaFile.nreferences,
                    "End": fastaFile.lengths,
                }
            )
        else:
            wlObj = read_bed(whiteListBed)
            wlObj = wlObj.merge()

        # if we have a blacklist, we subtract that from the whitelist, otherwise we leave it how
        # it is
        if not blackListBed is None:
            blObj = read_bed(blackListBed)
            blObj = blObj.merge()
            wlObj = wlObj.subtract(blObj)
            # shouldnt need to merge again here, as we only have less ranges than before

        # while we could use the get_fasta function from pyranges, it needs another
        # dependency (pyfaidx) and is slower (from my preliminary testing)
        # i terate over all chromosomes and each of the ranges
        for chr, df in wlObj:
            for region in df.iterrows():
                seq = fastaFile.fetch(
                    reference=chr, start=region["Start"], end=region["End"]
                )

                for i in range(len(seq) - 2):
                    diNucCounts[seq[i : i + 2]] += 1
                    triNucCounts[seq[i : i + 3]] += 1

    return (diNucCounts, triNucCounts)


def normaliseCounts(countsDf, contextCountDf):
    debug("Normalising counts with reference context frequency")
    # go over all columns and get the corresponding counts from the contexts, then make a vector from those counts in the same order
    countsLength = len(countsDf.columns)
    # create an empty array
    weights = empty(countsLength)
    for i in range(countsLength):
        # get the mutation type
        col = countsDf.columns[i]
        # snip the mutation off
        refContext = col[:3]
        # we also need the reverse complement of the context as we collapse the counts
        refContextRevComp = reverseComplement(refContext)

        weights[i] = contextCountDf[refContext] + contextCountDf[refContextRevComp]
        debug(f"{i}: {refContext}/{refContextRevComp} -> {weights[i]}")

    # we create the weights, but multiply them, so that we dont come into double precision range
    weights = divide(weights, (sum(weights) / 100))
    debug(f"Final weights: {weights}")

    # finally the normalisation (which is pretty easy)
    return divide(countsDf, weights)
