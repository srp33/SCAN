import os, sys, glob
import utilities

inFilePath = sys.argv[1]
columnIndex = int(sys.argv[2])

data = utilities.readMatrixFromFile(inFilePath)
values = [x[columnIndex] for x in data]
unique = list(set(values))

#countDict = {}
#for x in data:
#    value = x[columnIndex]
#    countDict[value] = countDict.setdefault(value, 0) + 1
#
#for key in countDict:
#    count = countDict[key]
#    if count > 1:
#        print key, count
#        for x in data:
#            if x[columnIndex] == key:
#                print x
#        exit()

print "%i unique values" % len(unique)
