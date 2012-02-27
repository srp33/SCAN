import os, sys, glob
import utilities

inFilePath = sys.argv[1]
splitColumnIndex = int(sys.argv[2])
splitText = sys.argv[3].decode('string-escape')
outFilePath = sys.argv[4]

data = utilities.readMatrixFromFile(inFilePath)

outData = []

for row in data:
    for splitItem in row[splitColumnIndex].split(splitText):
        newRow = list(row)
        newRow[splitColumnIndex] = splitItem
        outData.append(newRow)

utilities.writeMatrixToFile(outData, outFilePath)
