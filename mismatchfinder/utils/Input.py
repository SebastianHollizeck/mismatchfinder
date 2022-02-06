from argparse import ArgumentParser, ArgumentTypeError
from logging import basicConfig, debug, error, info, getLogger, root
from pathlib import Path
from sys import exit

from pysam import AlignmentFile

# stole this from stackoverflow
# (https://stackoverflow.com/questions/12116685/how-can-i-require-my-python-scripts-argument-to-be-a-float-between-0-0-1-0-usin)
def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


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
            "-Q",
            "--baseQuality",
            help="min base quality required to be counted as a variant [default: %(default)s]",
            type=int,
            default=25,
        )
        parser.add_argument(
            "--minAverageBaseQuality",
            help="min average base quality through the read to consider it in analysis [default: %(default)s]",
            type=int,
            default=20,
        )
        parser.add_argument(
            "-q",
            "--mappingQuality",
            help="min mapping quality required for a read to be included [default: %(default)s]",
            type=int,
            default=37,
        )
        parser.add_argument(
            "--maxMisMatchesPerRead",
            help="maximum mismatches a read is allowed before it is discarded [default: %(default)s]",
            type=int,
            default=7,
        )
        parser.add_argument(
            "--minMisMatchesPerRead",
            help="minimum mismatches a read needs before it is analysed [default: %(default)s]",
            type=int,
            default=1,
        )
        parser.add_argument(
            "--maxMisMatchesPerFragment",
            help="maximum mismatches a pair of reads is allowed before it is discarded [default: %(default)s]",
            type=int,
            default=7,
        )
        parser.add_argument(
            "--minMisMatchesPerFragment",
            help="minimum mismatches a pair of read needs before it is analysed [default: %(default)s]",
            type=int,
            default=1,
        )
        parser.add_argument(
            "--maxFragmentLengths",
            help="comma seperated list of maximum fragment lengths for a read to be considered for analysis [default: %(default)s] -1 corresponds to unlimited; will be checked with --minFragmentLengths for proper non overlapping intervals",
            type=str,
            default="144,325",
        )
        parser.add_argument(
            "--minFragmentLengths",
            help="comma seperated list of minimum fragment lengths for a read to be considered for analysis [default: %(default)s]; will be checked with --maxFragmentLengths for proper non overlapping intervals",
            type=str,
            default="74,240",
        )
        parser.add_argument(
            "--onlyOverlap",
            help="only consider the overlapping part of the paired end fragment [default: %(default)s]",
            action="store_true",
        )
        parser.set_defaults(onlyOverlap=False)
        parser.add_argument(
            "--strictOverlap",
            help="when enabled and the reads do not agree in the overlap, the mismatch will be discarded [default: %(default)s]",
            action="store_true",
        )
        parser.set_defaults(strictOverlap=False)

        parser.add_argument(
            "-t",
            "--threads",
            help="Amount of additional threads to use for decompression of the bam",
            type=int,
            default=1,
        )

        # parser.add_argument(
        #     "-O",
        #     "--outputType",
        #     help="the format to output the result",
        #     choices=["json", "R"],
        #     default="json",
        # )

        parser.add_argument(
            "-o",
            "--outFileRoot",
            help="the file to write output to (will be created if not present)",
            required=True,
            default=None,
        )
        parser.add_argument(
            "--overwrite",
            help="wether to overwrite previous results [default: %(default)s]",
            action="store_true",
        )
        parser.set_defaults(overwrite=False)
        parser.add_argument(
            "--writeEvidenceBam",
            help="wether to output the reads, which were used to calculate the mismatches as a bam [default: %(default)s]",
            action="store_true",
        )
        parser.set_defaults(writeEvidenceBam=False)
        parser.add_argument(
            "--writeEvidenceAsReadPairs",
            help="wether to output both reads, which were used to calculate the mismatches when writing the evidence bam [default: %(default)s]",
            action="store_true",
        )
        parser.set_defaults(writeEvidenceAsReadPairs=True)

        parser.add_argument(
            "-v",
            "--verbosity",
            help="Level of messages to be displayed",
            choices=["debug", "info", "error"],
            default="info",
        )

        parser.add_argument(
            "-m",
            "--method",
            help="The method to use to deconstruct signatures QP: Quadratic Programming ILM: Iterative Linear [default: %(default)s]",
            choices=["QP", "ILM"],
            default="QP",
        )

        parser.add_argument(
            "--germlineAFCutOff",
            help="minimum allele frequency for a variant in the germline resource to be used to filter [default: %(default)s]",
            type=restricted_float,
            default=0,
        )

        parser.add_argument(
            "--germlineRequirePass",
            help="Flag wheter the germline variant needs to be a PASS filter value to be used for filtering",
            action="store_true",
        )
        # we set this to false, because an issue in the germline calling can hint at a problem we
        # might have with the detection as well
        parser.set_defaults(germlineRequirePass=False)

        parser.add_argument(
            "--normaliseCounts",
            help="Flag to enable normalisation based on the occurrences of the (di/tri)-nucleotides",
            action="store_true",
        )
        # we set this to false, because an issue in the germline calling can hint at a problem we
        # might have with the detection as well
        parser.set_defaults(normaliseCounts=False)

        parser.add_argument(
            "--flatNormalisation",
            help="Flag to enable an even normalisation on the number of the contexts in the genome (default is to use the fraction of found contexts vs the total number in the genome)",
            action="store_true",
        )
        # we set this to false, because an issue in the germline calling can hint at a problem we
        # might have with the detection as well
        parser.set_defaults(flatNormalisation=False)

        # parser.add_argument(
        #     "-n",
        #     "--normals",
        #     help="the normal bams that should be used for plots",
        #     nargs="+",
        #     metavar="BAM",
        #     default=[],
        # )

        params = parser.parse_args()

        # set up the logging before we do anything else
        # now we tell people how much they want to know
        root.setLevel(params.verbosity.upper())

        ############################################################################################
        ###                                    sanity check                                      ###
        ############################################################################################

        # find out if we have a reference
        try:
            rFile = Path(params.reference)
            if rFile.is_file():
                self.referenceFile = params.reference
            else:
                self.referenceFile = None
        except TypeError:
            self.referenceFile = None

        # we better check if the bam files actually exist
        self.bamFiles = []
        for bam in params.bams:
            debug(f"Checking alignment input file: {bam}")
            # we might actually waive this exception (maybe with an option) if just one file is not
            # found? then again we dont want to just do half an analysis.
            bam = Path(bam)
            if not bam.is_file():
                raise Exception("Could not find bam file: " + bam)
            else:
                # this will be done multiple times, so you can combine bams and crams in the
                # analysis without any issues
                # we could theoretically do this with just file endings and it would probably be
                # faster, but then we also would need to check for the index
                with AlignmentFile(bam, "r") as tFile:
                    if tFile.is_cram and self.referenceFile is None:
                        raise Exception("CRAMs need a reference")
                self.bamFiles.append(bam)

        # we do the same test for the normals
        self.normals = []
        # for bam in params.normals:
        #     debug(f"Checking alignment input file: {bam}")
        #     # we might actually waive this exception (maybe with an option) if just one file is not
        #     # found? then again we dont want to just do half an analysis.
        #     bam = Path(bam)
        #     if not bam.is_file():
        #         raise Exception(f"Could not find bam file: {bam}")
        #     else:
        #         # this will be done multiple times, so you can combine bams and crams in the
        #         # analysis without any issues
        #         with AlignmentFile(bam, "r") as tFile:
        #             if tFile.is_cram and self.referenceFile is None:
        #                 raise Exception("CRAMs need a reference")
        #         self.normals.append(bam)

        self.blackListFile = None
        # we really only need to check if the file exists, if a file was actually given to us
        if not params.blackList is None:
            debug(f"Checking blacklist input file: {params.blackList}")

            blFile = Path(params.blackList)
            if not blFile.is_file():
                raise Exception(
                    f"Could not find black list bed file: {params.blackList}"
                )
            else:
                self.blackListFile = blFile

        self.whiteListFile = None
        if not params.whiteList is None:
            debug(f"Checking whitelist input file: {params.whiteList}")

            wlFile = Path(params.whiteList)
            if not wlFile.is_file():
                raise Exception(
                    f"Could not find whitelist bed file: {params.whiteList}"
                )
            else:
                self.whiteListFile = wlFile

        self.germlineFile = ""
        # we really only need to check if the file exists, if a file was actually given to us
        if not params.germline is None:
            debug(f"Checking germline input file: {params.germline}")

            glFile = Path(params.germline)
            if not glFile.is_dir():
                raise Exception(
                    f"Could not find germline zarr folder: {params.germline}"
                )
            else:
                self.germlineFile = glFile

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

        oFile = None
        if not params.outFileRoot is None:
            # I REALLY WANT TO TEST HERE IF THE directory allows write acces, as it would suck to run through the whole process only to not be able to write the results
            oFile = Path(params.outFileRoot)
        # we also want the None, so we do this outside of the if
        self.outFileRoot = oFile

        # dont need to check anything here, because the parser already does everything for us
        self.verbosity = params.verbosity.upper()
        self.threads = params.threads
        self.method = params.method

        self.minAverageBaseQuality = params.minAverageBaseQuality

        self.maxMisMatchesPerRead = params.maxMisMatchesPerRead
        self.minMisMatchesPerRead = params.minMisMatchesPerRead

        self.maxMisMatchesPerFragment = params.maxMisMatchesPerFragment
        self.minMisMatchesPerFragment = params.minMisMatchesPerFragment

        self.onlyOverlap = params.onlyOverlap
        self.strictOverlap = params.strictOverlap

        # if both strict overlap and only overlap are active, we actually double the
        # base quality, which we require, as the base quality of two agreeing reads
        # gets added up
        # TODO: think if we should add this in here, or if we should just leave it
        # to the user to specify a higher base quality

        self.writeEvidenceBam = params.writeEvidenceBam
        self.writeEvidenceReadPairs = params.writeEvidenceAsReadPairs
        self.overwrite = params.overwrite

        # here we check and parse the fragment length intervals
        minFragList = params.minFragLengths.split(",")
        maxFragList = params.maxFragLengths.split(",")

        if len(minFragList) != len(maxFragList):
            error(
                "Length of minimum and maximum fragment sizes does not match\n--minFragmentLengths and --maxFragmentLengths need to have the same length"
            )
            exit(1)

        self.fragmentLengthIntervals = []
        for min, max in zip(minFragList, maxFragList):
            try:
                minNum = int(min)
                maxNum = int(max)

                if maxNum == -1:
                    maxNum = float("inf")

                if minNum < 0 or minNum >= maxNum:
                    raise Exception("why you do the stupid?")

                self.fragmentIntervals.append((min, max))

            except Exception as e:
                error(
                    f"Specified fragment size interval is no proper size min:{min} max:{max}"
                )

        self.afCutOff = params.germlineAFCutOff
        self.germlineRequirePass = params.germlineRequirePass

        self.normaliseCounts = params.normaliseCounts

        self.flatNormalisation = params.flatNormalisation
