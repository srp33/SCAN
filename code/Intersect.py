import os, sys, glob

inFilePath1 = sys.argv[1]
inFilePath2 = sys.argv[2]
outFilePath = sys.argv[3]

values1 = set([x.rstrip() for x in file(inFilePath1)])
values2 = set([x.rstrip() for x in file(inFilePath2)])

intersection = list(values1 & values2)

outFile = open(outFilePath, 'w')
outFile.write("\n".join(intersection))
outFile.close()
