import os, sys, glob
import utilities

inFilePath = sys.argv[1]
filterExpression = sys.argv[2]
numHeaderRows = int(sys.argv[3])
outFilePath = sys.argv[4]

inFile = open(inFilePath)
outFile = open(outFilePath, 'w')

def filterAndWrite(data):
    data = filter(lambda x: eval(filterExpression), data)

    for row in data:
        outFile.write("\t".join(row) + "\n")

    return len(data)

for i in range(numHeaderRows):
    outFile.write(inFile.readline())

linesRead = 0
linesWritten = 0
data = []

for line in inFile:
    data.append(line.rstrip().split("\t"))
    linesRead += 1

    if len(data) == 100000:
        print "Data rows read: %i" % linesRead
        linesWritten += filterAndWrite(data)
        data = []

linesWritten += filterAndWrite(data)

outFile.close()
inFile.close()

print "Data rows read: %i" % linesRead
print "Data rows after filtering: %i" % linesWritten
