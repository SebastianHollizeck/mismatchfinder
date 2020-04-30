from numpy import array, int64, arange
from ncls import NCLS

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
