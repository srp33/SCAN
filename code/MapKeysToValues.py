import os, sys, glob
import utilities

inFilePath = sys.argv[1]
numHeaderRows = int(sys.argv[2])
keyColumnIndex = int(sys.argv[3])
valueColumnIndex = int(sys.argv[4])
outFilePath = sys.argv[5]

inFile = open(inFilePath)

for i in range(numHeaderRows):
    inFile.readline()

lineCount = 0
keyValueDict = {}
for line in inFile:
    lineCount += 1
    if lineCount % 10000 == 0:
        print lineCount
    lineItems = line.rstrip().split("\t")
    key = lineItems[keyColumnIndex]
    value = lineItems[valueColumnIndex]
    keyValueDict[key] = keyValueDict.setdefault(key, []) + [value]

inFile.close()

outData = []
for key in sorted(keyValueDict.keys()):
    outData.append([key, ",".join(keyValueDict[key])])

utilities.writeMatrixToFile(outData, outFilePath)
