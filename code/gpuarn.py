import ConfigParser, os
import logging
import random
import numpy as np
from time import clock
from numpy import array as nparray
from functools import partial
from bitstring import *
from math import exp
import sys
from sys import stdout
from subprocess import call
from utils import *
from utils.bitstrutils import *
from arn import *
import gpukernel
from pycuda import gpuarray, driver as drv


try:
    import matplotlib.pyplot as plt
    from matplotlib import collections, legend
except ImportError:
    logging.warning('matplotlib not found.')
    logging.warning('use of silentmode = 0 will result in an error...')

log = logging.getLogger(__name__)

TFACTORS = 0
STRUCTS = 1


def iterate(arnet,samplerate, simtime, silentmode, simstep,delta,**kwargs):
    s = clock()
    time = 0
    nump = arnet.numtf
    numsamples = 1
    numrec = arnet.numrec
    numeff = arnet.numeff


    cuweights_e = arnet.gpu_eweights
    cuweights_i = arnet.gpu_iweights

    #print kwargs['inputs']
    #cuccs = gpuarray.to_gpu(arnet.ccs)
    while time < simtime:
        #for i in range(numrec):
        #    arnet.ccs[nump+i] = kwargs['inputs'][i]

        _update(arnet.proteins,arnet.ccs,cuweights_e,cuweights_i,delta, **kwargs)
        #normalize ccs, ignoring receptors
        #print arnet.ccs[:nump]
        #cuccs.get(arnet.ccs)
        totparcels = (arnet.ccs[:nump].tolist() +
                      arnet.ccs[nump+numrec:].tolist())
        arnet.ccs /= sum(totparcels)
        #cuccs.get(arnet.ccs)
        #enforce again input/receptor cc
        #for i in range(numrec):
        #    arnet.ccs[nump+i] = kwargs['inputs'][i]

        if time % int(simtime*samplerate) == 0:
            log.debug('TIME: '+ str(time))
            arnet.cchistory = np.column_stack((arnet.cchistory,
                                               arnet.ccs[:nump]))
            arnet.receptorhist = np.column_stack((arnet.receptorhist,
                                                  arnet.ccs[nump:nump+numrec]))
            arnet.effectorhist = np.column_stack((arnet.effectorhist,
                                                  arnet.ccs[nump+numrec:]))
        time+=simstep


    print 'Elapsed time: %f sec.' % (clock()-s,)
    return arnet


def _update(proteins, cuccs, exciteweights, inhibitweights,delta,**kwargs):
    #print ccs.shape
    #print cuccs.shape
    esignals = _getSignalArrayCUDA(cuccs[:len(exciteweights)],
                                   exciteweights,**kwargs)
    isignals = _getSignalArrayCUDA(cuccs[:len(inhibitweights)],
                                   inhibitweights,**kwargs)
    deltas = esignals - isignals
    #print deltas.shape
    #print cuccs.shape
    deltas *= np.float32(delta)
    #kwargs['esignals'] *= np.float32(delta)
    deltas *= cuccs
    #kwargs['esignals'] = kwargs['esignals'] * cuccs
    #total = sum(ccs)+sum(deltas)
    cuccs += deltas
    #ccs/=total

def _getSignalArrayCUDA(ccs, weightstable, **kwargs):
    #print ccs.shape
    #print weightstable.shape[1]
    out = gpuarray.to_gpu(np.zeros(weightstable.shape[1], dtype=np.float32))
    cuccs = gpuarray.to_gpu(ccs)
    #print signals.shape
    #fun = gpukernel.getkernel(16,weightstable.shape[1])
    kwargs['gpudot'](out , cuccs, weightstable,
                     np.int32(weightstable.shape[0]),
                     np.int32(weightstable.shape[1]))
    return np.float32(1.0/len(ccs)) * out.get()
    #return np.float32(1.0/len(ccs)) *

def _getprotein(idx, code, bind_size, gene_size, protein_size):
    signature = BitStream(bin=
                        applymajority(code[bind_size*3:bind_size*3+gene_size],
                                      protein_size))
    #EXTENDED version - Weak linkage (needs double size gene/proteins)
    #p = [code[:self.bind_size],
     #       code[self.bind_size:self.bind_size*2],
      #      signature[0:self.bind_size],
       #     signature[self.bind_size:self.protein_size]]
    #ORIGINAL version
    p = [idx,
         code[:bind_size],
         code[bind_size:bind_size*2],
         signature,
         signature]
    log.debug(p)
    return p

def _getweights(bindings, bindingsize, beta, **kwargs):
    weights = bindings - bindingsize
    weights *= beta
    return np.exp(weights)

#organized in columns for the target equation
def getbindings(bindtype, proteins, match_threshold,**kwargs):
    return nparray([[XORmatching(p[3],otherps[1+bindtype],match_threshold)
                    for otherps in proteins]
                    for p in proteins],dtype=np.float32)
class ARNetwork:
    def __init__(self, gcode, config, **kwargs):
        self.code = gcode
        self.simtime = config.getint('default','simtime')

        promfun = bindparams(config, buildpromlist)
        productsfun = bindparams(config, buildproducts)
        self.promlist = promfun(gcode)
        self.proteins = productsfun( gcode, self.promlist)

        self.effectors=[]
        self.effectorproms = promfun(gcode, promoter='00000000')
        if self.effectorproms:
            #print 'EFFECTORS:', self.effectorproms
            self.effectors = productsfun(gcode,self.effectorproms)

        self.receptors=[]
        self.receptorproms = promfun(gcode, promoter='11111111')
        if self.receptorproms:
            #print 'RECEPTORS:', self.receptorproms
            self.receptors = productsfun(gcode,self.receptorproms)

        pbindfun = bindparams(config, getbindings)
        weightsfun = bindparams(config, _getweights)

        prob = kwargs['problem']
        self.numtf = len(self.proteins)
        self.numeff = min(len(self.effectors),prob.nout)
        self.effectors = self.effectors[:self.numeff]
        self.numrec = min(len(self.receptors),prob.ninp)
        self.receptors = self.receptors[:self.numrec]
        self.ccs = []
        if self.promlist:
            self.ccs = nparray([1.0/(self.numtf+self.numeff+self.numrec)]*
                               (self.numtf+self.numeff+self.numrec),
                               dtype=np.float32)
            self._initializehistory()
            self._initializebindings(pbindfun)
            self._initializeweights(weightsfun)
            self.dot_ = gpukernel.getkernel(22,self.eweights.shape[1])
            self.esignals = gpuarray.to_gpu(np.zeros(self.eweights.shape[1],
                                                dtype=np.float32))
            self.isignals = gpuarray.to_gpu(np.zeros(self.iweights.shape[1],
                                                dtype=np.float32))
            for i in range(len(self.proteins)):
                self.proteins[i].append(self.ccs[i])

        self.simfun = bindparams(config,iterate)
        self.delta = config.getfloat('default','delta')

    def _initializebindings(self, pbindfun):
        self.ebindings = pbindfun(0, self.proteins  + self.receptors +
                                  self.effectors)
        self.ibindings = pbindfun(1, self.proteins + self.receptors +
                                  self.effectors)
        #effectors and dummy receptors do not regulate
        if self.effectors or self.receptors:
                self.ebindings = self.ebindings[:self.numtf+self.numrec,:]
                #print 'ebinds: ', self.ebindings.shape
                self.ibindings = self.ibindings[:self.numtf+self.numrec,:]
                #print 'ibinds: ', self.ibindings.shape

    def _initializeweights(self, weightsfun):
        self.eweights = weightsfun(self.ebindings)
        #print 'ebinds: ', self.eweights.shape
        self.iweights = weightsfun(self.ibindings)
        #print 'ebinds: ', self.iweights.shape
        self.gpu_eweights = gpuarray.to_gpu(self.eweights)
        self.gpu_iweights = gpuarray.to_gpu(self.iweights)

    def _initializehistory(self):
        self.cchistory=nparray(self.ccs[:self.numtf])
        self.receptorhist=nparray(self.ccs[self.numtf:self.numtf+self.numrec])
        self.effectorhist=nparray(self.ccs[self.numtf+self.numrec:])

    def reset(self):
        self._initializehistory()
        self.ccs = nparray([1.0/(self.numtf+self.numeff+self.numrec)]*
                           (self.numtf+self.numeff+self.numrec))

    def __str__(self):
        return str(self.proteins)

    def simulate(self, *inputs):
        if not inputs:
            inputs = [.0]*self.numrec
        #print self.numtf
        #print self.numrec
        #print str(len(self.ccs) - (self.numtf+self.numrec))
        if self.simtime > 0:
            self.simfun(self,inputs = inputs,
                        gpudot = self.dot_,
                        esignals = self.esignals,
                        isignals = self.isignals)
            for i in range(self.numtf):
                self.proteins[i][-1] = self.ccs[i]
            #for i in range(len(self.effectors)):
               # self.effectors[i][-1] = self.ccs[i+self.numtf+self.numrec]

    def stepsimulate(self, proteins, ccs):
         _updatenonorm(proteins, ccs, self.eweights, self.iweights, self.delta)
         return ccs

    def nstepsim(self, n = 1000, *inputs):
        self.simfun(self, simtime = n, numeffectors = self.numeff,
                    inputs = inputs)
        for i in range(len(self.proteins)):
            self.proteins[i][-1] = self.ccs[i]

################################################################################
### Test                                                                ########
################################################################################

if __name__ == '__main__':
        _usecuda_ = True
        class Problem:
            pass
        p = Problem()
        p.ninp=3
        p.nout = 1
        arnconfigfile = '../configfiles/arnsimlong.cfg'
        saveto = 'genome.save'
        log.setLevel(logging.DEBUG)
        cfg = ConfigParser.ConfigParser()
        cfg.readfp(open(arnconfigfile))
        proteins=[]
        nump = 0
        prob_inp=[.0,.0,.0]
        try:
            f = open(sys.argv[1], 'r')
            genome = BitStream(bin=f.readline())
            arnet = ARNetwork(genome, cfg, problem = p)
            saveto = sys.argv[1]
        except:
            while nump < 200 or nump > 512 or numeff == 0 or not arnet.receptors:
                genome = BitStream(float=random.random(), length=32)
                for i in range(cfg.getint('default','initdm')):
                    genome = dm_event(genome,
                                      .02)

                    arnet = ARNetwork(genome, cfg, problem = p)
                    nump = len(arnet.promlist)
                    print nump
                    numeff = len(arnet.effectors)

        #f = open(saveto,'w')
        #f.write(genome.bin)
        #f.close
        #exit(0)
        #s = clock()
        arnet.simulate()
        #print clock() - s

        if not cfg.getint('default', 'silentmode'):
            displayARNresults(arnet.proteins, arnet.cchistory,
                              cfg.getfloat( 'default','samplerate'), temp=0)
            extralabels = ['R']*arnet.numrec + ['E']*arnet.numeff
            struct_prots = arnet.receptors + arnet.effectors
            hist = np.vstack((arnet.receptorhist,arnet.effectorhist))
            displayARNresults(struct_prots, hist,
                              cfg.getfloat( 'default','samplerate'),
                              temp = 1, extralabels = extralabels)

        for p in arnet.proteins: print p
        print 'effectors: ', arnet.effectors
        f = open(saveto,'w')
        f.write(genome.bin)
        f.close
        #print genome.bin
