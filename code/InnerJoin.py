import os, sys, glob
import utilities

inFilePath1 = sys.argv[1]
inFilePath2 = sys.argv[2]
joinColumnIndex1 = int(sys.argv[3])
joinColumnIndex2 = int(sys.argv[4])
outFilePath = sys.argv[5]
outputJoinColumn = sys.argv[6] == "True"

data1 = utilities.readMatrixFromFile(inFilePath1)
data2 = utilities.readMatrixFromFile(inFilePath2)

dataDict1 = {}
dataDict2 = {}

for row in data1:
    if dataDict1.has_key(row[joinColumnIndex1]):
        dataDict1[row[joinColumnIndex1]].append(row)
    else:
        dataDict1[row[joinColumnIndex1]] = [row]

for row in data2:
    if dataDict2.has_key(row[joinColumnIndex2]):
        dataDict2[row[joinColumnIndex2]].append(row)
    else:
        dataDict2[row[joinColumnIndex2]] = [row]

commonKeys = sorted(list(set(dataDict1.keys()) & set(dataDict2.keys())))

outData = []
for joinValue in commonKeys:
    for row1 in dataDict1[joinValue]:
        for row2 in dataDict2[joinValue]:
            row1Out = [row1[i] for i in range(len(row1)) if i != joinColumnIndex1]
            row2Out = [row2[i] for i in range(len(row2)) if i != joinColumnIndex2]

            out = []
            if outputJoinColumn:
                out.append(joinValue)

            out += row1Out
            out += row2Out

            outData.append(out)

utilities.writeMatrixToFile(outData, outFilePath)
