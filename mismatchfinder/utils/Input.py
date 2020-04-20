#!/usr/bin/env python3

from argparse import ArgumentParser, REMAINDER
from os.path import isfile, isdir, abspath
from pysam import AlignmentFile
from logging import debug, info, error, basicConfig


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
            "-O",
            "--outputType",
            help="the format to output the result",
            choices=["json", "R"],
            default="json",
        )

        parser.add_argument(
            "-v",
            "--verbosity",
            help="Level of messages to be displayed",
            choices=["debug", "info", "error"],
            default="info",
        )

        parser.add_argument(
            "-n",
            "--normals",
            help="the normal bams that should be used for plots",
            nargs="+",
            metavar="BAM",
            default=[],
        )
        params = parser.parse_args()

        # set up the logging before we do anything else
        # now we tell people how much they want to know

        basicConfig(
            level=params.verbosity.upper(), format="%(processName)-10s  %(message)s"
        )

        ############################################################################################
        ###                                    sanity check                                      ###
        ############################################################################################

        # find out if we have a reference
        if params.reference is None or not isfile(params.reference):
            self.referenceFile = None
        else:
            self.referenceFile = params.reference

        # we better check if the bam files actually exist
        self.bamFiles = []
        for bam in params.bams:
            debug(f"Checking alignment input file: {bam}")
            # we might actually waive this exception (maybe with an option) if just one file is not
            # found? then again we dont want to just do half an analysis.
            if not isfile(bam):
                raise Exception("Could not find bam file: " + bam)
            else:
                # this will be done multiple times, so you can combine bams and crams in the
                # analysis without any issues
                with AlignmentFile(bam, "r", require_index=True) as tFile:
                    if tFile.is_cram and self.referenceFile is None:
                        raise Exception("CRAMs need a reference")
                self.bamFiles.append(abspath(bam))

        # we do the same test for the normals
        self.normals = []
        for bam in params.normals:
            debug(f"Checking alignment input file: {bam}")
            # we might actually waive this exception (maybe with an option) if just one file is not
            # found? then again we dont want to just do half an analysis.
            if not isfile(bam):
                raise Exception("Could not find bam file: " + bam)
            else:
                # this will be done multiple times, so you can combine bams and crams in the
                # analysis without any issues
                with AlignmentFile(bam, "r", require_index=True) as tFile:
                    if tFile.is_cram and self.referenceFile is None:
                        raise Exception("CRAMs need a reference")
                self.normals.append(abspath(bam))

        self.blackListFile = None
        # we really only need to check if the file exists, if a file was actually given to us
        if not params.blackList is None:
            debug(f"Checking blacklist input file: {params.blackList}")
            if not isfile(params.blackList):
                raise Exception(
                    "Could not find black list bed file: " + params.blackList
                )
            else:
                self.blackListFile = abspath(params.blackList)

        self.whiteListFile = None
        if not params.whiteList is None:
            debug(f"Checking whitelist input file: {params.whiteList}")
            if not isfile(params.whiteList):
                raise Exception(
                    "Could not find whitelist bed file: " + params.whiteList
                )
            else:
                self.whiteListFile = abspath(params.whiteList)

        self.germlineFile = ""
        # we really only need to check if the file exists, if a file was actually given to us
        if not params.germline is None:
            debug(f"Checking germline input file: {params.germline}")
            if not isdir(params.germline):
                raise Exception(
                    "Could not find germline zarr folder: " + params.germline
                )
            else:
                self.germlineFile = abspath(params.germline)

        # some info if weird values get specified
        if params.baseQuality < 0:
            error("base quality needs to be a positive integer")
            exit(1)
        elif params.baseQuality < 20:
            info(
                f"Selected low base quality ({params.baseQuality}) might significantly affect results"
            )
        self.minBQ = params.baseQuality

        # same here, if weird values are set, we just tell and continue unless its completly wrong
        if params.mappingQuality < 0:
            raise Exception("base quality needs to be a positive integer or 0")

        elif params.mappingQuality > 60:
            info(
                "BWA caps mapping quality at 60, there will be no reads to consider if MQ>60"
            )
        self.minMQ = params.mappingQuality

        # dont need to check anything here, because the parser already does everything for us
        self.verbosity = params.verbosity.upper()
        self.threads = params.threads
