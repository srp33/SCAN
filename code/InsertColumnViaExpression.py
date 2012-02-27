import os, sys, glob
import utilities

inFilePath = sys.argv[1]
expression = sys.argv[2]
outFilePath = sys.argv[3]

outFile = open(outFilePath, 'w')
data = []
outCount = 0

for line in file(inFilePath):
    lineItems = line.rstrip().split("\t")
    data.append(lineItems)
    outCount += 1

    if len(data) == 100000:
        data = map(lambda x: eval("x + [" + expression + "]"), data)
        outFile.write("\n".join(["\t".join(x) for x in data]) + "\n")
        data = []
        print outCount

if len(data) > 0:
    data = map(lambda x: eval("x + [" + expression + "]"), data)
    outFile.write("\n".join(["\t".join(x) for x in data]) + "\n")
    print outCount

outFile.close()
