import glob,os,sys,time,posix
import utilities
from numpy import *
from numpy.linalg import *
from em import EM
from mycel import MyCEL

class Normalize(object):
    def __init__(self):
        self.code= {'A':0,'a':0,'G':1,'g':1,'C':2,'C':2,'T':3,'t':3}
        self.N = 50000

    def _encode(self,seq):
        return [self.code[s] for s in seq]

    def _sig(self,y,m):
        resid=y-m
        s2=sqrt(dot(resid,resid)/size(y))
        return s2

    def _sample(self,total,start=1):
        interval = total/self.N
        if interval <= 1:
            interval = 1
        return [i for i in range(start,total,interval)]

    def _quantile_normalize(self,x):
        xr = sort(x,axis=0)
        xm = mean(x,axis=1)
        return xm[argsort(x,axis=0)]

    def _design_matrix(self,PMProbe):
        x = zeros((PMProbe.shape[0],80), 'f')
        x[:,0] = sum(PMProbe == self.code['T'],1).astype('float32')
        j = 1
        for ibase in 'ACG':
            x[:, j:j+25] = (PMProbe == self.code[ibase])
            j += 25
        for ibase in 'ACGT':
            count = sum(PMProbe == self.code[ibase], 1).astype('float32')
            x[:, j] = count**2
            j += 1
        return x

    def normalize(self, probe_int, probe_seq, probe_type, outFilePath, outputExtraModelData=False):
        probe_ids = sorted(list(set(probe_seq.keys()) & set(probe_int.keys())))

        pseq = []
        pint = []
        ptype = []
        for p in probe_ids:
            pseq.append(probe_seq[p])
            pint.append(probe_int[p])
            ptype.append(probe_type[p])

        mx = self._design_matrix(array([self._encode(v) for v in pseq]))
        random.seed(1)
        randomNoise = array([random.random() / 10000000 for x in range(len(pint))])
        my = log2(array((pint)))
        my = my + randomNoise # This avoids a problem where too many identical values results in some bins having no values

        model = EM()
        samplingProbeIndices = self._sample(len(probe_ids))
        mainProbeIndices = [i for i in range(len(ptype)) if ptype[i] == "M"]
        samplingProbeIndices = list(set(samplingProbeIndices) & set(mainProbeIndices))

        nbins = 25
        b1, b2 = model.EM_vMix(my[samplingProbeIndices,],mx[samplingProbeIndices,], bins=nbins)
        m1 = dot(mx, b1)
        m2 = dot(mx, b2)

        binsize = 5000
        nGroups = int(ceil(size(my) / binsize))
        index = argsort(m1)
        y_norm = zeros(size(my), 'f')

        for i in arange(nGroups):
            tmp = index[(binsize * i):min([binsize * i + binsize, size(my)])]
            tmpSd = self._sig(my[tmp], m1[tmp])
            y_norm[tmp] = ((my[tmp] - m1[tmp]) / tmpSd).tolist()

        model.assign_bin(m1, bins=nbins)
        gam = model.vresp(my, mx)

        print "Outputting to %s" % outFilePath
        outFile = file(outFilePath, 'w')

        outLines = []
        for i, pid in enumerate(probe_ids):
            out = pid + "\t%.9f\t%.9f" % (y_norm[i], gam[i, 1])

            if outputExtraModelData:
                out += "\t%.9f\t%.9f\t%.9f" % (my[i], m1[i], m2[i])

            outLines.append(out)

            if len(outLines) % 100000 == 0:
                outFile.write("\n".join(outLines) + "\n")
                outLines = []

        if len(outLines) > 0:
            outFile.write("\n".join(outLines) + "\n")

        outFile.close()
