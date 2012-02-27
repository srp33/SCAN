import glob,os,sys,time,posix
import utilities
from mycel import MyCEL
from normalize import *

def getMetadata(probeSequenceTabFilePath):
    coord2probe = {}
    probe2seq = {}
    probe2type = {}

    probeSequenceTabFile = open(probeSequenceTabFilePath)

    for line in probeSequenceTabFile:
        lineItems = line.rstrip().split("\t")
        probe_id = lineItems[0]
        probe_coor_x = int(lineItems[1])
        probe_coor_y = int(lineItems[2])
        probe_seq = lineItems[3]
        probe_type = lineItems[4]

        coord2probe[(probe_coor_x, probe_coor_y)] = probe_id
        probe2seq[probe_id] = probe_seq
        probe2type[probe_id] = probe_type

    probeSequenceTabFile.close()

    return coord2probe, probe2seq, probe2type

celFilePath = sys.argv[1]
celFileName = os.path.basename(celFilePath)
probeSequenceTabFilePath = sys.argv[2]
outFilePath = sys.argv[3]

outputRawIntensities = False
if len(sys.argv) >= 5:
    outputRawIntensities = sys.argv[4] == "True"

cel = MyCEL()
norm = Normalize()

print "Reading annotations from %s " % probeSequenceTabFilePath
coord2probe, probe2seq, probe2type = getMetadata(probeSequenceTabFilePath)

print "Reading " + celFilePath
probe2intensity = cel.read_cel(celFilePath, coord2probe)

print "Normalizing to " + outFilePath
normValues = norm.normalize(probe2intensity, probe2seq, probe2type, outFilePath, outputRawIntensities)
