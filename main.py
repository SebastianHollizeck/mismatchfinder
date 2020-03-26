#!/usr/bin/env python3


from utils.Input import InputParser
from core.BedObject import BedObject
from core.BamScanner import BamScanner

# This will run everything.

# start with parameter Parsing
inputs = InputParser()

# generate the blacklist
blackList = BedObject(inputs.blackListFile)

# generate white list
whiteList = BedObject(inputs.whiteListFile)

for bam in inputs.bamFiles:
    print(bam)
    t = BamScanner(
        bamFile=bam,
        referenceFile=inputs.referenceFile,
        minBQ=inputs.minBQ,
        minMQ=inputs.minMQ,
        blackList=blackList,
        whiteList=whiteList,
    )
    t.start()
