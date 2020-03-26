#!/usr/bin/env python3

from argparse import ArgumentParser
from os.path import isfile, isdir
from pysam import AlignmentFile


class InputParser(object):
    """Class to parse any inputs"""

    def __init__(self):
        super(InputParser, self).__init__()

        # get the parsing done here
        parser = ArgumentParser()

        parser.add_argument(
            "bams", help="the aligned bam file (needs index)", metavar="BAM", nargs="+"
        )
        parser.add_argument(
            "-T",
            "--reference",
            help="the reference to use (only required for cram input)",
        )
        parser.add_argument(
            "-b", "--blackList", help="the black list bed file to filter"
        )
        parser.add_argument(
            "-w", "--whiteList", help="the white list bed file to filter"
        )
        parser.add_argument(
            "-g", "--germline", help="source of germline variants to filter"
        )
        parser.add_argument(
            "-q",
            "--baseQuality",
            help="min base quality required to be counted as a variant",
            type=int,
            default=25,
        )
        parser.add_argument(
            "-Q",
            "--mappingQuality",
            help="min mapping quality required for a read to be included",
            type=int,
            default=25,
        )

        parser.add_argument(
            "-t",
            "--threads",
            help="Amount of additional threads to use for decompression of the bam",
            type=int,
            default=1,
        )

        parser.add_argument(
            "-v",
            "--verbose",
            help="flag to enable verbose status output (this will use carriage return to write to the same line again, so dont use it in redirects)",
            action="count",
            default=0,
        )

        parser.add_argument(
            "-O",
            "--outputType",
            help="the format to output the result",
            choices=["json", "R"],
            default="json",
        )
        params = parser.parse_args()

        ############################################################################################
        ###                                    sanity check                                      ###
        ############################################################################################

        # we better check if the bam files actually exist
        self.bamFiles = []
        for bam in params.bams:
            # we might actually waive this exception (maybe with an option) if just one file is not
            # found? then again we dont want to just do half an analysis.
            if not isfile(bam):
                raise Exception("Could not find bam file: " + bam)
            else:
                # this will be done multiple times, so you can combine bams and crams in the
                # analysis without any issues
                with AlignmentFile(bam, "r", require_index=True) as tFile:
                    if tFile.is_cram:
                        # we need to check for reference here, because crams need one
                        if params.reference is None or not isfile(params.reference):
                            raise Exception("CRAMs need a reference")
                        else:
                            self.referenceFile = params.reference
                    else:
                        self.referenceFile = None
                self.bamFiles.append(bam)

        self.blackListFile = None
        # we really only need to check if the file exists, if a file was actually given to us
        if not params.blackList is None:
            if not isfile(params.blackList):
                raise Exception(
                    "Could not find black list bed file: " + params.blackList
                )
            else:
                self.blackListFile = params.blackList

        self.whiteListFile = None
        if not params.whiteList is None:
            if not isfile(params.whiteList):
                raise Exception(
                    "Could not find whitelist bed file: " + params.whiteList
                )
            else:
                self.whiteListFile = params.whiteList

        self.germlineFile = ""
        # we really only need to check if the file exists, if a file was actually given to us
        if not params.germline is None:
            if not isdir(params.germline):
                raise Exception(
                    "Could not find germline zarr folder: " + params.germline
                )
            else:
                self.germlineFile = params.germline

        # some info if weird values get specified
        if params.baseQuality < 0:
            print("base quality needs to be a positive integer")
            exit(1)
        elif params.baseQuality < 20:
            print(
                "lower base quality migh significantly affect the calculations",
                file=sys.stderr,
            )
        self.minBQ = params.baseQuality

        # same here, if weird values are set, we just tell and continue unless its completly wrong
        if params.mappingQuality < 0:
            raise Exception("base quality needs to be a positive integer or 0")

        elif params.mappingQuality > 60:
            print(
                "BWA caps mapping quality at 60, there will be no reads to consider if MQ>60",
                file=sys.stderr,
            )
        self.minMQ = params.mappingQuality

        # dont need to check anything here, because the parser already does everything for us
        self.verbose = params.verbose
