# E-M model to estimate model parameters.

from numpy import *
import numpy
from numpy.linalg import *

class EM:
    def __init__(self):
        pass

    def dnorm(self,y,m,s2):
        fy=1/(sqrt(2*pi*s2))*exp(-1/(2*s2)*(y-m)**2)
        return fy

    def vresp(self,y,X):
        vars0=numpy.array([self.s1[i] for i in self.bin])
        L0=self.dnorm(y,dot(X,self.b1),vars0)
        vars1=numpy.array([self.s2[i] for i in self.bin])
        L1=self.dnorm(y,dot(X,self.b2),vars1)

        gam=zeros((shape(y)[0],2),'f')
        gam[:,0]=self.p[0]*L0/(self.p[0]*L0+self.p[1]*L1)
        gam[:,1]=1-gam[:,0]

        return gam

    def beta(self,y,X,gam):
        sqgam=sqrt(gam)
        Xw=transpose(sqgam*transpose(X))
        yw=sqgam*y
        beta=dot(dot(inv(dot(transpose(Xw), Xw)),transpose(Xw)), yw)
        return beta

    def vbeta(self,y,X,gam,s2):
        vars=numpy.sqrt(numpy.array([s2[i] for i in self.bin]))
        sqgam=sqrt(gam)
        Xw=transpose(1/vars*sqgam*transpose(X))
        yw=1/vars*sqgam*y
        beta=dot(dot(inv(dot(transpose(Xw), Xw)),transpose(Xw)), yw)
        return beta

    def vsig(self,y,X,b,gam,bins):
        s2=numpy.zeros(bins)+1

        for i in arange(bins):
            ystar=y[self.bin==i]
            Xstar=X[self.bin==i,:]
            gamstar=gam[self.bin==i]+.01
            resid=ystar-dot(Xstar,b)
            s2[i]=dot(resid*gamstar,resid)/sum(gamstar)

        return s2

    def assign_bin(self, y, bins):
        sy=numpy.sort(y)
        n=len(y)
        quans=[sy[int(n*i/bins)] for i in range(1,bins)]
        self.bin = numpy.zeros(n,'i')
        for i in range(bins-1):
            self.bin[y>quans[i]]+=1

    def EM_vMix(self,y,X,p=.5,bins=10,conv=.01):
        print "Starting EM"
        self.assign_bin(y,bins)

        self.gam = zeros((len(y),2),'f')
        quan=sort(y)[int(p*len(y))-1]
        self.gam[:,0]=y<=quan
        self.gam[:,1]=y>quan

        self.p =mean(self.gam,axis=0) #numpy
        self.b1=self.beta(y,X,self.gam[:,0])
        self.b2=self.beta(y,X,self.gam[:,1])
        self.s1=self.vsig(y,X,self.b1,self.gam[:,0],bins)
        self.s2=self.vsig(y,X,self.b2,self.gam[:,1],bins)

        theta_old=concatenate((self.p,self.b1,self.s1,self.b2,self.s2))
        it=0
        c=1
        while c>conv and it<100:
            # Expectation Step:
            self.gam=self.vresp(y,X)

            #M-Step
            self.p = mean(self.gam,axis=0) #numpy
            self.b1=self.vbeta(y,X,self.gam[:,0],self.s1)
            self.assign_bin(dot(X,self.b1),bins) ## assign bins based on background expected value
            self.b2=self.vbeta(y,X,self.gam[:,1],self.s2)
            self.s1=self.vsig(y,X,self.b1,self.gam[:,0],bins)
            self.s2=self.vsig(y,X,self.b2,self.gam[:,1],bins)

            c = max(abs(concatenate((self.p,self.b1,self.s1,self.b2,self.s2))-theta_old)/theta_old)
            theta_old=concatenate((self.p,self.b1,self.s1,self.b2,self.s2))
            it+=1

            print "Iteration %i, c = %.6f" % (it, c)

        print 'Converged in', it, 'iterations. Proportion of background probes:', self.p[0]
        return self.b1, self.b2

'''
from em import *
import numpy

y=numpy.array([2.643254, 4.033662, 6.312421, 4.539054, 6.416255, 2.535445, 4.730555, 6.530505, 9.602606, 11.02814, 3.210855, 6.368144, 7.803356, 6.077369, 7.816515, 8.289626, 9.64034, 8.603153, 6.890807, 10.25251])

X=numpy.array([[1, 1, 1], [1, 2, 2], [1, 3, 3], [1, 4, 4], [1, 5, 5], [1, 6, 1], [1, 7, 2], [1, 8, 3], [1, 9, 4], [1, 10, 5], [1, 1, 1], [1, 2, 2], [1, 3, 3], [1, 4, 4], [1, 5, 5], [1, 6, 1], [1, 7, 2], [1, 8, 3], [1, 9, 4], [1, 10, 5]])

example=EM()

gam1=example.EM_Mix(y,X)
gam2=example.EM_uMix(y,X)
gam3=example.EM_vMix(y,X)

gam1
gam2
gam3


#debug:
self=example
bins=3;p=.5
self.assign_bin(y,bins)

self.gam = zeros((len(y),2),'f')
quan=sort(y)[int(p*len(y))-1]
self.gam[:,0]=y<=quan
self.gam[:,1]=y>quan
        
self.p =mean(self.gam,axis=0) #numpy
self.b1=self.beta(y,X,self.gam[:,0])
self.b2=self.beta(y,X,self.gam[:,1])
self.s1=self.vsig(y,X,self.b1,self.gam[:,0],bins)
self.s2=self.vsig(y,X,self.b2,self.gam[:,1],bins)
        
theta_old=concatenate((self.p,self.b1,self.s1,self.b2,self.s2))

self.gam=self.vresp(y,X)
      
            #M-Step
self.p =mean(self.gam,axis=0) #numpy
self.b1=self.vbeta(y,X,self.gam[:,0],self.s1)
self.assign_bin(dot(X,self.b1),bins) ## assign bins based on background expected value
self.b2=self.vbeta(y,X,self.gam[:,1],self.s2)
self.s1=self.vsig(y,X,self.b1,self.gam[:,0],bins)
self.s2=self.vsig(y,X,self.b2,self.gam[:,1],bins)
'''
