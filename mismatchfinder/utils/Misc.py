import linecache
import resource
import tracemalloc
from itertools import product
from logging import debug

from ncls import NCLS
from numpy import arange, array, int64


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


# with this function we limit the available memory for the process to 2Gb, which should make it so
# the garbage collection is stricter. THis is a requirement on high ram machines to not overload
# the system with multiple processes
def limitMemory():
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    debug(f"Found resource limits: (soft={soft},hard={hard}) ")
    debug(f"Setting resource limits to soft: 2 and hard: {hard}")
    resource.setrlimit(resource.RLIMIT_AS, (2048, 4096))


def display_top(snapshot, key_type="lineno", limit=10):
    snapshot = snapshot.filter_traces(
        (
            tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
            tracemalloc.Filter(False, "<unknown>"),
        )
    )
    top_stats = snapshot.statistics(key_type)

    debug("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        debug(
            "#%s: %s:%s: %.1f KiB"
            % (index, frame.filename, frame.lineno, stat.size / 1024)
        )
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            debug("    %s" % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        debug("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    debug("Total allocated size: %.1f KiB" % (total / 1024))
