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
            # we cast this to string, because pyranges wants string and we use the Path type
            wlObj = read_bed(str(whiteListBed))
            wlObj = wlObj.merge()

        # if we have a blacklist, we subtract that from the whitelist, otherwise we leave it how
        # it is
        if not blackListBed is None:
            # we cast this to string, because pyranges wants string and we use the Path type
            blObj = read_bed(str(blackListBed))
            blObj = blObj.merge()
            wlObj = wlObj.subtract(blObj)
            # shouldnt need to merge again here, as we only have less ranges than before

        # while we could use the get_fasta function from pyranges, it needs another
        # dependency (pyfaidx) and is slower (from my preliminary testing)
        # i terate over all chromosomes and each of the ranges
        for chr, df in wlObj:
            # iterrows has to return the index, even though we dont use it
            for idx, region in df.iterrows():
                seq = fastaFile.fetch(
                    reference=chr, start=region["Start"], end=region["End"]
                )

                for i in range(len(seq) - 2):
                    diNucCounts[seq[i : i + 2]] += 1
                    triNucCounts[seq[i : i + 3]] += 1
            debug(f"contect frequency analysis complete for chromsome {chr}")

    return (diNucCounts, triNucCounts)


# these counts are used to generate a weighting for normalisation and they were generated with
# Biostrings on the BSgenome.GRCh38
nucleotideCountsGRCh38 = {
    # diNucCounts
    "AA": 287025139,
    "AC": 148150331,
    "AG": 205752406,
    "AT": 226225785,
    "CA": 212880749,
    "CC": 151236932,
    "CG": 29401795,
    "CT": 205524144,
    "GA": 175847498,
    "GC": 124732844,
    "GG": 152432158,
    "GT": 148502457,
    "TA": 191400248,
    "TC": 174923630,
    "TG": 213928532,
    "TT": 289690054,
    # triNuc
    "AAA": 112465943,
    "AAC": 43532050,
    "AAG": 58439928,
    "AAT": 72587151,
    "ACA": 59305516,
    "ACC": 33784390,
    "ACG": 7584302,
    "ACT": 47476086,
    "AGA": 65552680,
    "AGC": 41073623,
    "AGG": 51723263,
    "AGT": 47402783,
    "ATA": 60308591,
    "ATC": 39076747,
    "ATG": 53548035,
    "ATT": 73292370,
    "CAA": 55220609,
    "CAC": 44001434,
    "CAG": 59791771,
    "CAT": 53866888,
    "CCA": 53293160,
    "CCC": 38036593,
    "CCG": 8026845,
    "CCT": 51880303,
    "CGA": 6511692,
    "CGC": 7021552,
    "CGG": 8229568,
    "CGT": 7638969,
    "CTA": 37666053,
    "CTC": 49481013,
    "CTG": 59039769,
    "CTT": 59337262,
    "GAA": 58990420,
    "GAC": 27737004,
    "GAG": 49560877,
    "GAT": 39559024,
    "GCA": 42481943,
    "GCC": 34497599,
    "GCG": 7078395,
    "GCT": 40674873,
    "GGA": 46022042,
    "GGC": 34474720,
    "GGG": 38148838,
    "GGT": 33786518,
    "GTA": 33265786,
    "GTC": 27466578,
    "GTG": 44578403,
    "GTT": 43191653,
    "TAA": 60348082,
    "TAC": 32879810,
    "TAG": 37959659,
    "TAT": 60212654,
    "TCA": 57800075,
    "TCC": 44918305,
    "TCG": 6712244,
    "TCT": 65492835,
    "TGA": 57760931,
    "TGC": 42162935,
    "TGG": 54330453,
    "TGT": 59674158,
    "TTA": 60159779,
    "TTC": 58899235,
    "TTG": 56762262,
    "TTT": 113868707,
}


def normaliseCounts(countsDf, contextCountDf, flatNorm=True):
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

        # we generate the weight, by comparing the count we had in our analysis region to the counts
        # in the whole genome and adjust them that way (or just a flat normalisation)
        if flatNorm:
            baseline = nucleotideCountsGRCh38[refContext]
            +nucleotideCountsGRCh38[refContextRevComp]
        else:
            baseline = 1
        weights[i] = (
            contextCountDf[refContext] + contextCountDf[refContextRevComp] / baseline
        )
        debug(f"{i}: {refContext}/{refContextRevComp} -> {weights[i]}")

    # we create the weights, but multiply them, so that we dont come into double precision range
    weights = divide(weights, (sum(weights) / 100))
    debug(f"Final weights: {weights}")

    # finally the normalisation (which is pretty easy)
    return divide(countsDf, weights)
