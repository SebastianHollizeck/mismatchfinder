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
