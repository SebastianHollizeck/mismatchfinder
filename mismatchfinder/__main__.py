#!/usr/bin/env python3


from logging import debug, info
from multiprocessing import Semaphore, SimpleQueue

from .core.BamScanner import BamScanner
from .core.BedObject import BedObject
from .core.GermlineObject import GermlineObject
from .results.Results import convertToPandasDataFrame
from .utils.Input import InputParser
from .utils.Output import createOutputFiles, writeStatsFile

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

    # store all created processes
    processes = []

    # store the results of the parallel processing
    tumourResults = SimpleQueue()
    normalResults = SimpleQueue()

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
            results=tumourResults,
            germObj=germline,
            outFileRoot=inputs.outFileRoot,
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
            results=normalResults,
            germObj=germline,
            outFileRoot=inputs.outFileRoot,
        )
        # we request the resources here, but return them inside the thread, once it is done
        semaphore.acquire()
        # then we start the process
        p.start()
        processes.append(p)

    debug("waiting for join")
    # wait for all processes to finish before we continue
    for p in processes:
        p.join()

    statsTable = convertToPandasDataFrame(tumourResults)
    writeStatsFile(statsTable, outFileRoot=inputs.outFileRoot)

    # output results as just the data or even with plot
    # plotStats(
    #     convertToPandasDataFrame(tumourResults),
    #     normalDf=convertToPandasDataFrame(normalResults),
    # )
    info("FINISHED")


if __name__ == "__main__":
    main()
