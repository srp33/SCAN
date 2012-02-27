import os, sys, glob
import utilities
from operator import itemgetter, attrgetter

inFilePath = sys.argv[1]
columnIndex = int(sys.argv[2])
reverse = sys.argv[3] == "reverse"
numHeaderRows = int(sys.argv[4])

data = utilities.readMatrixFromFile(inFilePath)

headerRows = []
for i in range(numHeaderRows):
    headerRows.append(data.pop(0))

data.sort(key=itemgetter(columnIndex), reverse=reverse)

outFile = open(inFilePath, 'w')
for row in headerRows:
    outFile.write("\t".join(row) + "\n")

for row in data:
    outFile.write("\t".join([str(x) for x in row]) + "\n")

outFile.close()
