#!/usr/bin/env python3


from .utils.Input import InputParser
from .results.Results import convertToPandasDataFrame
from .utils.Output import plotStats
from .core.BedObject import BedObject
from .core.GermlineObject import GermlineObject
from .core.BamScanner import BamScanner
from multiprocessing import Semaphore, SimpleQueue
from logging import basicConfig, debug, DEBUG

# This will run everything.


def main():
    basicConfig(level=DEBUG, format="%(processName)-10s  %(message)s")
    # start with parameter Parsing
    inputs = InputParser()

    # generate the blacklist
    blackList = BedObject.parseFile(inputs.blackListFile)

    # generate white list
    whiteList = BedObject.parseFile(inputs.whiteListFile)

    # generate germline object
    germline = GermlineObject.parseZarrRoot(inputs.germlineFile)

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
            bamFile=bam,
            referenceFile=inputs.referenceFile,
            minBQ=inputs.minBQ,
            minMQ=inputs.minMQ,
            blackList=blackList,
            whiteList=whiteList,
            semaphore=semaphore,
            results=tumourResults,
            germObj=germline,
        )
        # this shedules the running, but will be blocked by the semaphore inside the actual process
        # to not overload the system
        p.start()
        processes.append(p)

    # spawn processes to analyse each bam individually
    for bam in inputs.normals:
        p = BamScanner(
            bamFile=bam,
            referenceFile=inputs.referenceFile,
            minBQ=inputs.minBQ,
            minMQ=inputs.minMQ,
            blackList=blackList,
            whiteList=whiteList,
            semaphore=semaphore,
            results=normalResults,
            germObj=germline,
        )
        # this shedules the running, but will be blocked by the semaphore inside the actual process
        # to not overload the system
        p.start()
        processes.append(p)

    # wait for all processes to finish before we continue
    for p in processes:
        p.join()

    # output results as just the data or even with plot
    # plotStats(
    #     convertToPandasDataFrame(tumourResults),
    #     normalDf=convertToPandasDataFrame(normalResults),
    # )
    # print(tumourResults.get())


if __name__ == "__main__":
    main()
