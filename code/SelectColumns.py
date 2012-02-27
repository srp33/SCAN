import os, sys, glob
import utilities

inFilePath = sys.argv[1]
columnIndices = [int(x) for x in sys.argv[2].split(",")]
outFilePath = sys.argv[3]

lineCount = 0
outFile = open(outFilePath, 'w')
out = ""

for line in file(inFilePath):
    lineCount += 1
    lineItems = line.rstrip().split("\t")
    outItems = [lineItems[i] for i in columnIndices]

    out += "\t".join(outItems) + "\n"

    if lineCount % 100000 == 0:
        print lineCount
        outFile.write(out)
        out = ""

outFile.write(out)
outFile.close()
