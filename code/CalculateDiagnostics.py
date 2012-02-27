import os, sys, glob, math
import utilities

inFilePattern = sys.argv[1]
outFilePath = sys.argv[2]

outFile = open(outFilePath, 'w')
outFile.write("\t".join(["Sample", "MeanProbCode", "R2BackgroundActual", "R2PercentTotalVariation", "PercentVariationBackground", "PreSignalNoise", "PostSignalNoise", "PostPreSignalNoise"]) + "\n")

for inFilePath in sorted(glob.glob(inFilePattern)):
    fileName = os.path.basename(inFilePath)
    print "Processing %s" % fileName

    normalized = []
    probCodes = []
    raw = []
    m1 = []
    m2 = []
    percentTotalVariation = []

    for line in file(inFilePath):
        lineItems = line.rstrip().split("\t")

        pProbe = float(lineItems[2])
        m1Probe = float(lineItems[4])
        m2Probe = float(lineItems[5])

        normalized.append(float(lineItems[1]))
        probCodes.append(pProbe)
        raw.append(float(lineItems[3]))
        m1.append(m1Probe)
        m2.append(m2Probe)
        percentTotalVariation.append(((1 - pProbe) * m1Probe) + (pProbe * m2Probe))

    meanProbCode = utilities.calculateMean(probCodes)
    r2BackgroundActual = math.pow(utilities.calculateCorrelationCoefficient(raw, m1), 2)
    r2PercentTotalVariation = math.pow(utilities.calculateCorrelationCoefficient(raw, percentTotalVariation), 2)
    percentVariationBackground = r2BackgroundActual / r2PercentTotalVariation
    preSignalNoise = (r2PercentTotalVariation - r2BackgroundActual) / r2BackgroundActual
    postSignalNoise = preSignalNoise / (1 - r2PercentTotalVariation)
    postPreSignalNoise = postSignalNoise / preSignalNoise

    outFile.write("\t".join([str(x) for x in [fileName, meanProbCode, r2BackgroundActual, r2PercentTotalVariation, percentVariationBackground, preSignalNoise, postSignalNoise, postPreSignalNoise]]) + "\n")
    outFile.flush()

outFile.close()
