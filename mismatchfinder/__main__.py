#!/usr/bin/env python3


from logging import debug, info, basicConfig, getLogger

# set the logging up
basicConfig(format="%(asctime)s %(processName)-13s %(levelname)-6s %(message)s",)

from multiprocessing import Lock, Semaphore, SimpleQueue, set_start_method

# until this bug is fixed, we have this in here
import warnings

warnings.filterwarnings("ignore")

from .core.BamScanner import BamScanner
from .core.BedObject import BedObject
from .core.GermlineObject import GermlineObject
from .core.Signatures import Signature
from .utils.Input import InputParser
from .utils.Output import createOutputFiles


# This will run everything.


def main():

    # start with parameter Parsing
    inputs = InputParser()

    # generate the blacklist
    debug("Parsing blacklist file")
    blackList = BedObject.parseFile(inputs.blackListFile)

    # generate white list
    whiteList = BedObject.parseFile(inputs.whiteListFile)

    # generate germline object
    germline = GermlineObject.parseZarrRoot(inputs.germlineFile)

    createOutputFiles(inputs.outFileRoot)

    # create a block to analyse only the specified amount of processes in parallel
    semaphore = Semaphore(inputs.threads)
    # create a lock for file write access to ensure everything is save
    fhLock = Lock()

    # store all created processes
    processes = []

    # store the results of the parallel processing
    tumourResults = SimpleQueue()
    normalResults = SimpleQueue()

    logger = getLogger()

    # spawn processes to analyse each bam individually
    for bam in inputs.bamFiles:
        p = BamScanner(
            bamFilePath=bam,
            referenceFile=inputs.referenceFile,
            minBQ=inputs.minBQ,
            minMQ=inputs.minMQ,
            blackList=blackList,
            whiteList=whiteList,
            semaphore=semaphore,
            lock=fhLock,
            results=tumourResults,
            germObj=germline,
            outFileRoot=inputs.outFileRoot,
            log=logger,
        )
        # we request the resources here, but return them inside the thread, once it is done
        semaphore.acquire()
        # then we start the process
        p.start()
        processes.append(p)

    # spawn processes to analyse each bam individually
    for bam in inputs.normals:
        p = BamScanner(
            bamFilePath=bam,
            referenceFile=inputs.referenceFile,
            minBQ=inputs.minBQ,
            minMQ=inputs.minMQ,
            blackList=blackList,
            whiteList=whiteList,
            semaphore=semaphore,
            lock=fhLock,
            results=normalResults,
            germObj=germline,
            outFileRoot=inputs.outFileRoot,
            log=logger,
        )
        # we request the resources here, but return them inside the thread, once it is done
        semaphore.acquire()
        # then we start the process
        p.start()
        processes.append(p)

    # TODO: check if we need to consume the queue here already to reduce the memory footprint, or
    # if the memory requirement is actually from the zarr storage being cached so often

    debug("Waiting for all parallel processes to finish")
    # wait for all processes to finish before we continue
    for p in processes:
        p.join()
        p.close()

    del processes

    # destroy the semaphore after we are done
    semaphore = None
    fhLock = None

    # read the files of the parallel processed back in and analyse the results
    SBSsig = Signature.loadSignaturesFromFile(type="SBS")
    SBSFile = inputs.outFileRoot.parent / (inputs.outFileRoot.name + "_SBScontexts.tsv")
    SBSweights = SBSsig.analyseCountsFile(SBSFile, method=inputs.method)
    SBSweightsFile = inputs.outFileRoot.parent / (
        inputs.outFileRoot.name + "_SBSweights.tsv"
    )
    debug(f"Writing SBS signature weights to {SBSweightsFile}")
    SBSweights.to_csv(SBSweightsFile, sep="\t")

    DBSsig = Signature.loadSignaturesFromFile(type="DBS")
    DBSFile = inputs.outFileRoot.parent / (inputs.outFileRoot.name + "_DBScontexts.tsv")
    DBSweights = DBSsig.analyseCountsFile(DBSFile, method=inputs.method)
    DBSweightsFile = inputs.outFileRoot.parent / (
        inputs.outFileRoot.name + "_DBSweights.tsv"
    )
    debug(f"Writing DBS signature weights to {DBSweightsFile}")
    DBSweights.to_csv(DBSweightsFile, sep="\t")

    info("FINISHED")


if __name__ == "__main__":
    # we specify how to create the subprocesses
    # we use spawn to keep the memory footprint of the process smaller, even though creationg takes
    # a bit longer. But computation of the process will be the limiting factor anyways
    set_start_method("spawn")

    # and then run the whole thing
    main()
