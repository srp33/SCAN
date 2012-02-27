import os, sys, glob
import utilities

inFilePath = sys.argv[1]
filterColumnIndex = int(sys.argv[2])
filterValuesFilePathOrValue = sys.argv[3]
numHeaderRows = int(sys.argv[4])
outFilePath = sys.argv[5]

print "Reading data from %s" % inFilePath
data = utilities.readMatrixFromFile(inFilePath)

headerRows = []
for i in range(numHeaderRows):
    headerRows.append(data.pop(0))

if os.path.exists(filterValuesFilePathOrValue):
    filterValues = [x.rstrip() for x in file(filterValuesFilePathOrValue)]
else:
    filterValues = filterValuesFilePathOrValue.split(",")
    if len(filterValues) == 1:
        print "Warning: No file %s exists." % filterValues[0]

print "Filtering data"
values = [row[filterColumnIndex] for row in data]

valueIndexDict = {}
for x in enumerate(values):
    index = x[0]
    value = x[1]
    valueIndexDict[value] = index

indicesToKeep = [valueIndexDict[value] for value in filterValues if valueIndexDict.has_key(value)]
missingFilterValues = [value for value in filterValues if not valueIndexDict.has_key(value)]

data = [data[i] for i in indicesToKeep]
print "%i rows remain after filtering" % len(data)

print "Outputting data"
outFile = open(outFilePath, 'w')
for row in headerRows:
    outFile.write("\t".join(row) + "\n")
for row in data:
    outFile.write("\t".join(row) + "\n")
outFile.close()
