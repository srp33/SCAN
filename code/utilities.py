import glob, os, posix, sys, math, collections
import scipy
from scipy.stats import *
from operator import itemgetter, attrgetter
import itertools
from random import uniform, sample
import numpy
from collections import defaultdict

def printMatrix(data):
    for x in data:
        print x
    print ""

def getProbes(probeTabFilePath):
    probes = []

    probeTabFile = open(probeTabFilePath)
    headerItems = [x.lower() for x in probeTabFile.readline().rstrip().split("\t")]

    for line in probeTabFile:
        lineItems = line.rstrip().split("\t")
        if headerItems.count("probe set name") > 0:
            probeset = lineItems[headerItems.index("probe set name")]
        else:
            if headerItems.count("probe set id") > 0:
                probeset = lineItems[headerItems.index("probe set id")]
            else:
                print "No probe set name or probe set id column in %s" % probeTabFilePath

        probeX = lineItems[headerItems.index("probe x")]
        probeY = lineItems[headerItems.index("probe y")]
        probe = probeset + "#" + probeX + "_" + probeY
        probes.append(probe)

    return probes

def getProbesetProbesDict(probes):
    probesetProbesDict = {}

    for probe in probes:
        probeset = probe[:probe.find("#")]
        probesetProbesDict[probeset] = probesetProbesDict.setdefault(probeset, []) + [probe]

    return probesetProbesDict

def getPatientIDs(normDirPath, normFileSuffix):
    ids = []

    #print normDirPath + "*" + normFileSuffix
    #sys.exit(0)
    for filePath in glob.glob(normDirPath + "*" + normFileSuffix):
        ids.append(filePath.replace(normDirPath, "").replace(normFileSuffix, ""))

    ids.sort()
    return ids

def readScalarFromFile(filePath):
    return readMatrixFromFile(filePath)[0][0]

def writeScalarToFile(x, filePath):
    outFile = open(filePath, 'w')
    outFile.write(x)
    outFile.close()

def readVectorFromFile(filePath):
    return [line.rstrip() for line in file(filePath)]

def writeVectorToFile(data, filePath):
    outFile = open(filePath, 'w')
    for x in data:
        outFile.write(str(x) + "\n")
    outFile.close()

def readMatrixFromFile(filePath, numLines=None):
    matrix = []
    for line in file(filePath):
        if numLines != None and len(matrix) >= numLines:
            break

        matrix.append(line.rstrip().split("\t"))

        if len(matrix) % 100000 == 0:
            print len(matrix)

    return matrix

def writeMatrixToFile(x, filePath):
    outFile = open(filePath, 'w')
    for y in x:
        outFile.write("\t".join(y) + "\n")
    outFile.close()

def readTextFromFile(filePath):
    text = ""

    for line in file(filePath):
        text += line

    return text

def calculateMean(values):
    return sum(values) / len(values)

def calculateVarianceMean(values):
    mu = calculateMean(values)
    diffValues = [(x - mu)**2 for x in values]
    return calculateMean(diffValues) / (len(diffValues) - 1)

def calculateWeightedMean(values, weights):
    if len(values) != len(weights):
        print "When calculating a weighted mean, the values must be the same length as the weights."
        raise

def calculateStandardDeviation(values):
    xbar = calculateMean(values)
    residuals = [x - xbar for x in values]
    residualsSquared = [x**2 for x in residuals]
    return math.sqrt(sum(residualsSquared) / (len(values) - 1))

def calculateZscore(x):
    mean = calculateMean(x)
    standardDeviation = calculateStandardDeviation(x)
    return [(y - mean) / standardDeviation for y in x]

def calculateTrimmedMean(values, trimProportion=0.10):
    if values == None or len(values) == 0:
        return None

    values = sorted([float(x) for x in values])

    if len(values) < 3:
        return calculateMean(values)
    elif len(values) == 3:
        return values[1]
    elif len(values) == 4:
        return calculateMean(values[1:3])
    elif len(values) == 5:
        return calculateMean(values[1:4])

    values = scipy.stats.trimboth(values, trimProportion)

    return float(calculateMean(values))

def calculateEuclideanDistance(xList, yList):
    zSum = 0.0

    for i in range(len(xList)):
        x = xList[i]
        y = yList[i]
        z = math.pow(x - y, 2)
        zSum += z

    return math.sqrt(zSum)

def calculateCorrelationCoefficient(xList, yList):
    return numpy.corrcoef(xList, yList)[0,1]

def calculateTTest(xList, yList):
    return ttest_ind(xList, yList, 0)[1]

def calculateMedian(values):
  sortedValues = sorted(values)

  if len(sortedValues) % 2 == 1:
      return sortedValues[(len(sortedValues)+1)/2-1]
  else:
      lower = sortedValues[len(sortedValues)/2-1]
      upper = sortedValues[len(sortedValues)/2]
      return (float(lower + upper)) / 2

def calculateFoldChange(values1, values2):
    overallMin = min(min(values1), min(values2))

    values1 = [x - overallMin + 1 for x in values1]
    values2 = [x - overallMin + 1 for x in values2]

    mean1 = calculateMean(values1)
    mean2 = calculateMean(values2)

    return mean1 / mean2

def calculateAbsoluteFoldChange(values1, values2):
    overallMin = min(min(values1), min(values2))

    values1 = [x - overallMin + 1 for x in values1]
    values2 = [x - overallMin + 1 for x in values2]

    mean1 = calculateMean(values1)
    mean2 = calculateMean(values2)

    ratioA = mean1 / mean2
    ratioB = mean2 / mean1

    return min(ratioA, ratioB)

def getNormalizedProbes(normFilePath):
    print "Getting normalized probes"
    return [line.split(" ")[0] for line in file(normFilePath)]

def getKeyProbeDict(filePath):
    keyProbeDict = {}

    for line in file(filePath):
        lineItems = line.rstrip().split("\t")
        key = lineItems[0]

        if len(lineItems) > 1:
            fileProbes = [x for x in lineItems[1].split(",") if x != ""]

            if len(fileProbes) > 0:
                keyProbeDict[key] = keyProbeDict.setdefault(key, []) + fileProbes

    return keyProbeDict

def getTranscriptProbeDict(filePath, normFilePath):
    normalizedProbes = set(getNormalizedProbes(normFilePath))

    print "Getting transcript-probe dictionary"
    transcriptProbeDict = {}
    for line in file(filePath):
        lineItems = line.rstrip().split("\t")
        transcript = lineItems[0]
        probes = lineItems[1].split(",")
        probes = list(set(probes) & normalizedProbes)

        transcriptProbeDict[transcript] = probes

    return transcriptProbeDict

def getPatientsKeyValuesDict(sourceDir, patientIDs, fileSuffix, dataValueIndex, keys=None):
    patientsKeyValuesDict = collections.defaultdict(dict)

    if len(patientIDs) == 0:
        return patientsKeyValuesDict

    keyLineIndicesDict = {}
    lineCount = 0

    for line in file(sourceDir + patientIDs[0] + fileSuffix):
        key = line.rstrip().split("\t")[0]
        keyLineIndicesDict[key] = lineCount

        lineCount += 1
        #if lineCount % 100000 == 0:
        #    print "Parsing file line indices: %i" % lineCount

    #print "Creating key line indices list from dict"
    if keys == None:
        keyLineIndices = [(key, keyLineIndicesDict[key]) for key in keyLineIndicesDict.keys()]
    else:
        keyLineIndices = [(key, keyLineIndicesDict[key]) for key in keys if key in keyLineIndicesDict.keys()]

    #print "Sorting key line indices"
    keyLineIndices.sort(key=itemgetter(1))

    patientFileHandles = {}
    for patientID in patientIDs:
        patientFileHandles[patientID] = open(checkDirPath(sourceDir) + patientID + fileSuffix)

    for patientID in patientIDs:
        #print patientID
        patientFile = open(checkDirPath(sourceDir) + patientID + fileSuffix)

        previousLineIndex = 0
        for keyLineIndex in keyLineIndices:
            for i in range(previousLineIndex, keyLineIndex[1]):
                patientFile.readline()
            previousLineIndex = keyLineIndex[1] + 1

            lineItems = patientFile.readline().rstrip().split("\t")
            patientsKeyValuesDict[patientID][lineItems[0]] = lineItems[dataValueIndex]

        patientFile.close()

    return patientsKeyValuesDict

def getPatientKeyValuesDict(filePath, dataColumnIndex, probes=None):
    probeValues = {}

    for line in file(filePath):
        lineItems = line.rstrip().split("\t")
        probe = lineItems[0]
        value = lineItems[dataColumnIndex]

        probeValues[probe] = value

    if not probes:
        return probeValues
    else:
        modProbeValues = {}
        for probe in probes:
            modProbeValues[probe] = probeValues[probe]
        return modProbeValues

def savePatientKeyValuesDict(patientDict, outFilePath):
    outFile = open(outFilePath, 'w')

    for key in sorted(patientDict.keys()):
        outFile.write("%s\t%s\n" % (key, patientDict[key]))

    outFile.close()

#This is slower
#def getPatientProbeValuesDict(sourceDir, dataColumnIndex, patientID, probes=None):
#    probeValues = {}
#
#    for line in file(sourceDir + patientID + ".norm.txt"):
#        lineItems = line.rstrip().split(" ")
#
#        if probes == None or lineItems[0] in probes:
#            probeValues[lineItems[0]] = lineItems[dataColumnIndex]
#
#    return probeValues

def checkDirPath(dirPath):
    if not os.path.exists(dirPath):
        posix.mkdir(dirPath)

    if not dirPath.endswith("/"):
        dirPath = dirPath + "/"

    return dirPath

def lastIndexOf(theList, value):
    return len(theList) - 1 - theList[::-1].index(value)

def getTranscriptGeneDict(filePath):
    transcriptGeneDict = {}

    for line in file(filePath):
        lineItems = line.rstrip().split("\t")
        transcript = lineItems[0]

        gene = lineItems[1]
        if len(lineItems) == 3:
            gene = lineItems[2]

        transcriptGeneDict[transcript] = gene

    return transcriptGeneDict

def getGeneTranscriptDict(filePath):
    geneTranscriptDict = {}

    for line in file(filePath):
        lineItems = line.rstrip().split("\t")
        transcript = lineItems[0]

        gene = lineItems[1]
        if len(lineItems) == 3:
            gene = lineItems[2]

        geneTranscriptDict[gene] = geneTranscriptDict.setdefault(gene, []) + [transcript]

    return geneTranscriptDict

def transposeMatrix(x):
    transposed = zip(*x)

    for i in range(len(transposed)):
        transposed[i] = list(transposed[i])

    return transposed

# Copied from: http://code.activestate.com/recipes/491268-ordering-and-ranking-for-lists/
def order(x, NoneIsLast = True, decreasing = False):
    """
    Returns the ordering of the elements of x. The list
    [ x[j] for j in order(x) ] is a sorted version of x.

    Missing values in x are indicated by None. If NoneIsLast is true,
    then missing values are ordered to be at the end.
    Otherwise, they are ordered at the beginning.
    """
    omitNone = False
    if NoneIsLast == None:
        NoneIsLast = True
        omitNone = True

    n  = len(x)
    ix = range(n)
    if None not in x:
        ix.sort(reverse = decreasing, key = lambda j : x[j])
    else:
        # Handle None values properly.
        def key(i, x = x):
            elem = x[i]
            # Valid values are True or False only.
            if decreasing == NoneIsLast:
                return not(elem is None), elem
            else:
                return elem is None, elem
        ix = range(n)
        ix.sort(key=key, reverse=decreasing)

    if omitNone:
        n = len(x)
        for i in range(n-1, -1, -1):
            if x[ix[i]] == None:
                n -= 1
        return ix[:n]
    return ix

# Copied from: http://code.activestate.com/recipes/491268-ordering-and-ranking-for-lists/
def rankSmart(x, NoneIsLast=True, decreasing = False, ties = "first"):
    """
    Returns the ranking of the elements of x. The position of the first
    element in the original vector is rank[0] in the sorted vector.

    Missing values are indicated by None.  Calls the order() function.
    Ties are NOT averaged by default. Choices are:
                 "first" "average" "min" "max" "random" "average"
    """
    omitNone = False
    if NoneIsLast == None:
        NoneIsLast = True
        omitNone = True
    O = order(x, NoneIsLast = NoneIsLast, decreasing = decreasing)
    R = O[:]
    n = len(O)
    for i in range(n):
        R[O[i]] = i
    if ties == "first" or ties not in ["first", "average", "min", "max", "random"]:
        return R

    blocks     = []
    isnewblock = True
    newblock   = []
    for i in range(1,n) :
        if x[O[i]] == x[O[i-1]]:
            if i-1 not in newblock:
                newblock.append(i-1)
            newblock.append(i)
        else:
            if len(newblock) > 0:
                blocks.append(newblock)
                newblock = []
    if len(newblock) > 0:
        blocks.append(newblock)

    for i, block  in enumerate(blocks):
        # Don't process blocks of None values.
        if x[O[block[0]]] == None:
            continue
        if ties == "average":
            s = 0.0
            for j in block:
                s += j
            s /= float(len(block))
            for j in block:
                R[O[j]] = s
        elif ties == "min":
            s = min(block)
            for j in block:
                R[O[j]] = s
        elif ties == "max":
            s =max(block)
            for j in block:
                R[O[j]] = s
        elif ties == "random":
            s = sample([O[i] for i in block], len(block))
            for i,j in enumerate(block):
                R[O[j]] = s[i]
        else:
            for i,j in enumerate(block):
                R[O[j]] = j
    if omitNone:
        R = [ R[j] for j in range(n) if x[j] != None]
    return R

# The following function came from http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
def rank2(a):
    n = len(a)
    ivec=rank_simple(a)
    svec=[a[rank] for rank in ivec]
    sumranks = 0
    dupcount = 0
    newarray = [0]*n
    for i in xrange(n):
        sumranks += i
        dupcount += 1
        if i==n-1 or svec[i] != svec[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in xrange(i-dupcount+1,i+1):
                newarray[ivec[j]] = averank
            sumranks = 0
            dupcount = 0

    return newarray

def globFilesSortedByModTime(pattern):
    def getModifiedTime(filename):
        return os.stat(filename).st_mtime

    return sorted(glob.glob(pattern), key=getModifiedTime)

## From http://stackoverflow.com/questions/34518/natural-sorting-algorithm
def naturalSort(x, reverse=False):
    def natural_key(s):
        return tuple(
            int(''.join(chars)) if isdigit else ''.join(chars)
            for isdigit, chars in itertools.groupby(s, str.isdigit)
        )

    return sorted(x, key=natural_key, reverse=reverse)

def getItemFrequencyMap(x):
    d = defaultdict(int)
    for item in x:
        d[item] += 1

    return d
