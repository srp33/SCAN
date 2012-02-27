import os, sys, glob
from operator import itemgetter, attrgetter
import utilities

inFilePath1 = sys.argv[1]
inFilePath2 = sys.argv[2]
joinColumnIndex1 = int(sys.argv[3])
joinColumnIndex2 = int(sys.argv[4])
outFilePath = sys.argv[5]
outputJoinColumn = sys.argv[6] == "True"

outFile = open(outFilePath, 'w')
outLines = []
outLineCount = 0

inFile2 = open(inFilePath2)
lineItems2 = inFile2.readline().rstrip().split("\t")

for line in file(inFilePath1):
    lineItems1 = line.rstrip().split("\t")

    key1 = lineItems1[joinColumnIndex1]
    key2 = lineItems2[joinColumnIndex2]

    if key1.isdigit():
        key1 = int(key1)
        key2 = int(key2)

    if key1 == key2:
        out1 = list(lineItems1)
        out2 = list(lineItems2)
        del out1[joinColumnIndex1]
        del out2[joinColumnIndex2]

        out = out1 + out2
        if outputJoinColumn:
            out.insert(0, str(key1))

        outLines.append("\t".join(out))
        outLineCount += 1

        if len(outLines) == 100000:
            outFile.write("\n".join(outLines) + "\n")
            outLines = []
            print outLineCount

        lineItems2 = inFile2.readline().rstrip().split("\t")
    elif key1 > key2:
        lineItems2 = inFile2.readline().rstrip().split("\t")

    if len(lineItems2) < len(lineItems1):
        break

if len(outLines) > 0:
    outFile.write("\n".join(outLines) + "\n")
    print outLineCount

inFile2.close()
outFile.close()
