from zarr import open_group
from ncls import NCLS
from numpy import array, int64, arange
from mismatchfinder.utils.Misc import countLowerCase, buildNCLSindex
from logging import debug, info, error
from sys import exit


class GermlineObject(object):
    """docstring for GermlineObject."""

    def __init__(self, zarrRootFolder):
        super(GermlineObject, self).__init__()
        self.zarrObj = open_group(zarrRootFolder, mode="r")

    # This is a very small function, but it has a heavy backbone attached. this zarr storage was
    # precomputed via scikit-allel, to contain all of gnomad, which then can be loaded with almost
    # no computational cost involved.
    # if anything fails in this, we return an empty object, which signals downstream methods to not
    # do any tests for germline status
    # the hg38 version of gnomad was converted with
    # # import sys, allel, pysam
    # #
    # # zarr_path= '/data/reference/dawson_labs/gnomad/3/zarr'
    # # vcf_path = '/data/databases/gnomAD/3/gnomad.genomes.r3.0.sites.vcf.bgz'
    # #
    # # varsFH = pysam.VariantFile(vcf_path)
    # # chrs = list(varsFH.header.contigs)
    # #
    # # for chr in chrs:
    # #     try:
    # #         allel.vcf_to_zarr(vcf_path, zarr_path, group=f"{chr}", fields='*', log=sys.stdout,
    # #               tabix="tabix", region=f"{chr}")
    # #     except ValueError:
    # #         # this means that this already exists, which is fine, we just continue with the
    # #         # next chr
    # #         print(f"chr {chr} was already found in zarr storage")
    #
    @classmethod
    def parseZarrRoot(cls, dir):
        # if we do not get a dir we just return None
        if dir is None or dir == "":
            return None
        try:
            return GermlineObject(dir)
        except Exception as e:
            # we do not print anything here, because this is just to get None for an empty input
            # if the zarr thing does not exist
            return None

    def loadChromosomeIndex(self, chr):
        debug(f"loading index for chr {chr}")
        # build the NCLS index instead of what allel does, because the caching is
        # slower but the query is SO much faster
        try:
            return buildNCLSindex(self.zarrObj[f"{chr}/variants/POS"][...])
        except Exception as e:
            debug("Exception {0}".format(e))
            debug(f"No variants found in germline resource for chr {chr}")
            raise KeyError(f"No variants found for chr {chr}")

    def loadChromsomeAlts(self, chr):
        # this loads the zarr object into memory, which is possible because we only load
        # a small part of the whole thing
        return self.zarrObj[f"{chr}/variants/ALT"][...]

    def loadChromsomeRefs(self, chr):
        return self.zarrObj[f"{chr}/variants/REF"][...]

    def checkGermlineStatus(self, candMisMatches, discard=False):

        if len(candMisMatches) == 0:
            debug(f"we didnt find any sites, so no germline status to be checked")
            debug(f"Checking germline status for {len(candMisMatches)} mismatches")
        nShiftFail = 0

        res = {}

        cache = None
        # we go through all sites, for caching purposes its important here to have a sorted list
        # otherwise we have to refresh the cache with every change of the chromosome
        for (chr, pos, refContext, altContext) in sorted(candMisMatches):

            # first we need to cache the data from the chromosome that the site is on
            if cache is None or cache.chr != chr:
                cache = ChromosomeCache(self, chr)
                # if there is no cache for this chromosome we just get none back
                if cache is None:
                    # in this case we cannot check the germline status and just assume everything is
                    # somatic, but we can discard sites on chromosomes when we dont know the status
                    if discard:
                        continue

            occurrences = candMisMatches[(chr, pos, refContext, altContext)]
            # get the index of all germline variants with a trinucleotide at this position
            varIdxs = cache.findOverlaps(pos)

            # we do not really need the start and stop here, we could just use the idx, but
            # NCLS gives the full info anyways
            for (start, stop, idx) in varIdxs:
                # this is the position in the contexts we are talking about (again pysam is 0
                # based which makes this part a tad easier)
                offset = start - pos

                # we can ignore this position if it is upper case in the ref, because that would
                # mean we didnt find a variant there
                if refContext[offset].isupper():
                    continue

                # if the ref doesnt have length 1 we are dealing with an indel, and we dont care
                ref = cache.refs[idx]
                if len(ref) != 1:
                    continue

                # we can have several alts, so we check all of them
                for alt in cache.alts[idx]:
                    # again if this is not of length 1 we have an indel and we do not really
                    # want to deal with it
                    if len(alt) != 1:
                        continue
                    # now here we need to check if the alt matches the alt that we found in the
                    # reads
                    if altContext[offset].upper() == alt:
                        # this means this is a germline variant and we just make it an upper case
                        # to keep the variant but show its not somatic
                        refContext = (
                            f"{refContext[:offset]}{ref}{refContext[offset+1:]}"
                        )
                        # we can also stop iterating through the alts, because only one ref and alt
                        # combination will match with the observed variant
                        break
            # with us changing around things for germline variants, we need to discard the
            # cases, that now dont have variants anymore (at least no somatic ones)
            if countLowerCase(refContext) == 0:
                continue
            else:
                if refContext[1].isupper():
                    nShiftFail += 1
            res[(chr, pos, refContext, altContext)] = occurrences

        info(f"Found {nShiftFail} trinucleotides without a center variant")
        return res


class ChromosomeCache(object):
    """Very simple object that just stores the cached information from the zarr object"""

    def __init__(self, germlineObj, chr):
        super(ChromosomeCache, self).__init__()
        debug(f"creating cache for chr {chr}")
        try:
            self.index = germlineObj.loadChromosomeIndex(chr)
            self.refs = germlineObj.loadChromsomeRefs(chr)
            self.alts = germlineObj.loadChromsomeAlts(chr)
            self.chr = chr
        except KeyError:
            # if this fails it means we couldnt load the index, because it doesnt exist in the cache
            debug(f"Couldnt find variants for chr {chr} in the cache")
            self.chr = None
        except Exception as e:
            eStr = getattr(e, "message", repr(e))
            error(
                f"Unknown exception {eStr} '{str(e)}' when trying to cache from zarr storage"
            )
            exit(1)

    def findOverlaps(self, pos, loff=0, roff=2):
        # because pysam gives 0 based indexes, we want pos:pos+1 to get all variants in the
        # first two positions of the trinucleotide context, which is enough, as we left
        # shifted variants and didnt allow a variant in the third position
        # getting the index which corresponds to the chromosomal position
        try:
            res = self.index.find_overlap(pos + loff, pos + roff)

        except AttributeError:
            # this is when there are no germline variants reported for this site
            res = []
        except Exception as e:
            eStr = getattr(e, "message", repr(e))
            error(f"Unknown exception {eStr} '{str(e)}' when querying cache")
            exit(1)
        return res
