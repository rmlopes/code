import ConfigParser, os
import logging
import random
import numpy as np
from numpy import array as nparray
from functools import partial
from bitstring import *
from math import exp
import sys
from sys import stdout
from subprocess import call
from utils import *
from utils.bitstrutils import *
from time import clock
from arn import bindparams, generatechromo, buildpromlist, \
    buildproducts, getbindings, _getweights, _getSignalArray
from extendedarn import displayARNresults
import copy

log = logging.getLogger(__name__)


INPUT_SIGNATURES = [BitStream(bin=('0'*32)),
                    BitStream(bin=('0'*16 + '1'*16)),
                    BitStream(bin=('1'*16 + '0'*16)),
                    BitStream(bin=('1'*32))]

def build_customproducts(signatures = INPUT_SIGNATURES):
    products = []
    i = 0
    for sig in signatures:
        products.append([i, sig, sig, sig, sig])
        i += 1
    return products


def iterate(arnet,samplerate, simtime, silentmode, simstep,delta,**kwargs):
    #s = clock()
    time = 1
    nump = arnet.numtf
    numsamples = 1
    numrec = arnet.numrec
    numeff = arnet.numeff

    #print kwargs['inputs']
    while time <= simtime:
        for i in range(numrec):
            arnet.ccs[nump+i] = kwargs['inputs'][i]

        _update(arnet.proteins,arnet.ccs,arnet.eweights,arnet.iweights,
                delta, numtf = nump)
        #normalize ccs, ignoring effectors
        #print arnet.ccs[:nump]
        totparcels = arnet.ccs[:nump+numrec].tolist()
        #receptors ccs is not modified here
        arnet.ccs[:nump] /= sum(totparcels)
        #normalize outputs
        #FIXME: should this be done?
        #for i in range(numeff):
        #    if arnet.ccs[-1-i] > 1.0:
        #        arnet.ccs[-1-i] = 1.0
        #    elif arnet.ccs[-1-i] < .0:
        #        arnet.ccs[-1-i] = .0
        if numeff > 1:
            arnet.ccs[nump+numrec:] /= sum(arnet.ccs[nump+numrec:])

        if time % int(simtime*samplerate) == 0:
            log.debug('TIME: '+ str(time))
            arnet.updatehistory()
        time+=simstep

    #print 'Elapsed time: %f sec.' % (clock()-s,)
    return arnet


def _update(proteins, ccs, exciteweights, inhibitweights,delta,**kwargs):
    #print ccs.shape
    deltas = (_getSignalArray(ccs[:len(exciteweights)],exciteweights) -
              _getSignalArray(ccs[:len(inhibitweights)],inhibitweights))
    #print deltas.shape
    deltas *= delta

    #NOTE: previous output concentration shall not be accounted for
    numreg = exciteweights.shape[0]
    #NOTE: not touching receptors ccs
    deltas[:kwargs['numtf']] *= ccs[:kwargs['numtf']]
    ccs[:kwargs['numtf']] += deltas[:kwargs['numtf']]
    ccs[numreg:] += deltas[numreg:]
    #ccs/=total

class ARNetwork:
    def __init__(self, gcode, config, **kwargs):
        self.code = gcode
        self.simtime = config.getint('default','simtime')

        promfun = bindparams(config, buildpromlist)
        productsfun = bindparams(config, buildproducts)
        self.promlist = promfun(gcode)
        self.proteins = productsfun( gcode, self.promlist)

        self.effectors=[]
        self.effectorproms = promfun(gcode, promoter='11111111')
        if self.effectorproms:
            #print 'EFFECTORS:', self.effectorproms
            self.effectors = productsfun(gcode,self.effectorproms)

        self.receptors= build_customproducts()
        #self.receptorproms = promfun(gcode, promoter='11111111')
        #if self.receptorproms:
            #print 'RECEPTORS:', self.receptorproms
            #self.receptors = productsfun(gcode,self.receptorproms)

        pbindfun = bindparams(config, getbindings)
        weightsfun = bindparams(config, _getweights)

        prob = kwargs['problem']
        self.numtf = len(self.proteins)
        #self.numeff = min(len(self.effectors),prob.nout)
        self.numeff = len(self.effectors)
        #self.effectors = self.effectors[:self.numeff]
        self.numrec = min(4,prob.ninp)
        self.receptors = self.receptors[:self.numrec]
        self.ccs = []
        if self.promlist:
            self.ccs = nparray([1.0/(self.numtf+self.numeff+self.numrec)]*
                               (self.numtf+self.numeff+self.numrec))
            self._initializehistory()
            self._initializebindings(pbindfun)
            self._initializeweights(weightsfun)
            for i in range(len(self.proteins)):
                self.proteins[i].append(self.ccs[i])

        self.simfun = bindparams(config,iterate)
        self.delta = config.getfloat('default','delta')
        self.output_idx = 0

    def _initializebindings(self, pbindfun):
        self.ebindings = pbindfun(0, self.proteins  + self.receptors +
                                  self.effectors)
        self.ibindings = pbindfun(1, self.proteins + self.receptors +
                                  self.effectors)

        #effectors do not regulate
        if self.effectors:
                self.ebindings = self.ebindings[:self.numtf+self.numrec,:]
                #print 'ebinds: ', self.ebindings.shape
                self.ibindings = self.ibindings[:self.numtf+self.numrec,:]
                #print 'ibinds: ', self.ibindings.shape

    def _initializeweights(self, weightsfun):
        self.eweights = weightsfun(self.ebindings)
        #print 'ebinds: ', self.eweights.shape
        self.iweights = weightsfun(self.ibindings)
        #print 'ebinds: ', self.iweights.shape

    def _initializehistory(self):
        self.cchistory=nparray(self.ccs[:self.numtf])
        self.receptorhist=nparray(self.ccs[self.numtf:self.numtf+self.numrec])
        self.effectorhist=nparray(self.ccs[self.numtf+self.numrec:])

    def updatehistory(self):
        nump = self.numtf
        self.cchistory = np.column_stack((self.cchistory,
                                          self.ccs[:nump]))
        self.receptorhist = np.column_stack((self.receptorhist,
                                             self.ccs[nump:nump+self.numrec]))
        self.effectorhist = np.column_stack((self.effectorhist,
                                             self.ccs[nump+self.numrec:]))

    def reset(self, cc_state = []):
        if len(cc_state) == 0:
            self.ccs = nparray([1.0/(self.numtf+self.numeff+self.numrec)]*
                           (self.numtf+self.numeff+self.numrec))
        else:
            #NOTE: although some defend the use of [:] instead of
            #deepcopy, even if it is a list of floats would only
            #shallow copy, resulting in inconsistency
            self.ccs = copy.deepcopy(cc_state)
            #print 'reseting state to: ', cc_state
        self._initializehistory()
        #FIXME: pickled ccs are not correct
        #print cc_state
        #self.updatehistory()
        #self.nstepsim(self.simtime,*[.0,.0,.0,.0])
        #print self.ccs
        #if len(cc_state) == 0:

        #else:
           # self.ccs = nparray(cc_state[:])


    def __str__(self):
        return str(self.proteins)

    def simulate(self, *inputs):
        if not inputs:
            inputs = [.0]*self.numrec
        #print self.numtf
        #print self.numrec
        #print str(len(self.ccs) - (self.numtf+self.numrec))
        if self.simtime > 0:
            self.simfun(self,inputs = inputs)
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
        arnconfigfile = '../configfiles/arnsim-miguel.cfg'
        class Problem:
            pass
        p = Problem()
        p.ninp=4
        p.nout = 1
        log.setLevel(logging.DEBUG)
        cfg = ConfigParser.ConfigParser()
        cfg.readfp(open(arnconfigfile))
        proteins=[]
        nump = 0
        prob_inp=[.0,.0,.0,.0]
        try:
            f = open(sys.argv[1], 'r')
            genome = BitStream(bin=f.readline())
            arnet = ARNetwork(genome, cfg, problem = p)
        except:
            while nump < 4 or nump > 32 or numeff == 0 or not arnet.receptors:
                genome = BitStream(float=random.random(), length=32)
                for i in range(cfg.getint('default','initdm')):
                    genome = dm_event(genome,
                                      .02)

                    arnet = ARNetwork(genome, cfg, problem = p)
                    nump = len(arnet.promlist)
                    numeff = len(arnet.effectors)

        arnet.simulate()
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
        print 'receptors: ', arnet.receptors
        #f = open('genome.save','w')
        #f.write(genome.bin)
        #f.close
        #print genome.bin
