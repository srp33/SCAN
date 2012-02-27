import glob, os, posix, shutil, sys
import utilities
import scipy.stats
from collections import defaultdict

def calculateVarianceMean(values):
    return utilities.calculateVarianceMean(values)

def calculateMean(values):
    return utilities.calculateMean(values)

def calculateTrimmedMean(values):
    return utilities.calculateTrimmedMean(values)

inFilePath = sys.argv[1]
dataColumnIndex = int(sys.argv[2])
keyProbeFilePath = sys.argv[3]
minNumProbesPer = int(sys.argv[4])
summarizeFunction = getattr(sys.modules[__name__], sys.argv[5])
outFilePath = sys.argv[6]

print "Getting key/probe dict from %s" % keyProbeFilePath
keyProbeDict = utilities.getKeyProbeDict(keyProbeFilePath)

print "Removing keys with few probes"
keysToRemove = []
for key in keyProbeDict.keys():
    if len(keyProbeDict[key]) < minNumProbesPer:
        keysToRemove.append(key)

if len(keysToRemove) > 0:
    print "Removing %i keys because they have few probes" % len(keysToRemove)
    for key in keysToRemove:
        del keyProbeDict[key]

print "Building reverse key/value dict"
reverseKeyProbeDict = {}
for key, probes in keyProbeDict.iteritems():
    probes = tuple(sorted(probes))
    reverseKeyProbeDict[probes] = reverseKeyProbeDict.setdefault(probes, []) + [key]

print "Combining duplicate keys"
keyProbeDict2 = {}
for key in reverseKeyProbeDict.keys():
    valueList = reverseKeyProbeDict[key]
    ### Note: If multiple keys have the same probes, we just pick the first one in the list
    modValue = valueList[0]
    keyProbeDict2[modValue] = key

print "Reading data from %s" % inFilePath
probeValuesDict = utilities.getPatientKeyValuesDict(inFilePath, dataColumnIndex)

outFile = open(outFilePath, 'w')

for key in sorted(keyProbeDict2.keys()):
    probeValues = [float(probeValuesDict[probe]) for probe in keyProbeDict2[key] if probeValuesDict.has_key(probe)]

    if len(probeValues) == 0:
        continue

    summaryValue = summarizeFunction(probeValues)
#    print key
#    print summarizeFunction
#    print probeValues
#    print summaryValue
    outFile.write("%s\t%.9f\n" % (key, summaryValue))

outFile.close()
