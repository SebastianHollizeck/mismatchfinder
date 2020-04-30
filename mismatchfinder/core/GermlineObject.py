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
        self.__zarrObj = open_group(zarrRootFolder, mode="r")

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
            return GermlineObject(dir.as_posix())
        except Exception as e:
            # we do not print anything here, because this is just to get None for an empty input
            # if the zarr thing does not exist
            return None

    def loadChromosomeIndex(self, chr):
        # build the NCLS index instead of what allel does, because the caching is
        # slower but the query is SO much faster
        return buildNCLSindex(self.__zarrObj[f"{chr}/variants/POS"][...])

    def loadChromsomeAlts(self, chr):
        # this loads the zarr object into memory, which is possible because we only load
        # a small part of the whole thing
        return self.__zarrObj[f"{chr}/variants/ALT"][...]

    def loadChromsomeRefs(self, chr):
        return self.__zarrObj[f"{chr}/variants/REF"][...]


class ChromosomeCache(object):
    """Very simple object that just stores the cached information from the zarr object"""

    def __init__(self, germlineObj, chr):
        super(ChromosomeCache, self).__init__()
        debug(f"creating cache for chr {chr}")
        try:
            self.index = germlineObj.loadChromosomeIndex(chr)
            self.refs = germlineObj.loadChromsomeRefs(chr)
            self.alts = germlineObj.loadChromsomeAlts(chr)

        except KeyError:
            # if this fails it means we couldnt load the index, because it doesnt exist in the cache
            debug(f"Zarr storage did not contain chr {chr}")

        except Exception as e:
            eStr = getattr(e, "message", repr(e))
            error(
                f"Unknown exception {eStr} '{str(e)}' when trying to cache from zarr storage"
            )
            exit(1)
        # we set this so we dont try to cache again for the same chromosome
        self.chr = chr

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
