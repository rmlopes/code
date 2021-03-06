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
import arn
from arn import bindparams, generatechromo, \
    buildproducts, getbindings, _getweights, _getSignalArray
from extendedarn import displayARNresults
import copy
from bitarray import *

log = logging.getLogger(__name__)

#transcription factors aux index
TF = 0
#products aux index
PR = 1
'''
INPUT_SIGNATURES = [BitStream(bin=('0'*32)),
                    BitStream(bin=('1'*16 + '0'*16)),
                    BitStream(bin=('0'*16 + '1'*16)),
                    BitStream(bin=('1'*32)),
                    BitStream(bin=('0111'*8)),
                    BitStream(bin=('0011'*8)),
                    BitStream(bin=('0001'*8)),
                    BitStream(bin=('0101'*8)),
                    BitStream(bin=('00001111'*4)),
                    BitStream(bin=('11110000'*4)),
                    BitStream(bin=('0000000011111111'*2)),
                    BitStream(bin=('1111111100000000'*2))]
'''
INPUT_SIGNATURES = [32*bitarray('0'),
                    16*bitarray('1') + 16*bitarray('0'),
                    16*bitarray('0') + 16*bitarray('1'),
                    32*bitarray('1'),
                    8*bitarray('0111'),
                    8*bitarray('0011'),
                    8*bitarray('0001'),
                    8*bitarray('0101'),
                    4*bitarray('00001111'),
                    4*bitarray('11110000'),
                    2*bitarray('0000000011111111'),
                    2*bitarray('1111111100000000')]

def _get_all(promoter, genome, excite_offset,genesize):
    promsize = len(promoter)
    plist = genome.search(bitarray(promoter))
    plist = filter( lambda index:
                     int(excite_offset) <= index <  (genome.length()-(int(genesize)+promsize )),
                     plist)
    return plist

def _filteroverlapping(plist, genesize):
    #plist of tuples (prom, type)
    return reduce(lambda indxlst, indx:
                  indxlst + [indx] if indx[0]-indxlst[-1][0] >=  genesize + 32 + 64 else indxlst,
                  plist[1:],
                  plist[:1])

def _filteroverlaptf(plist, efflist):
    return filter(lambda x: not any([(e-256 < x < e + 256) for e in efflist]),plist)

def buildpromlist(genome, excite_offset, genesize,
                  promoter, productprom, **kwargs):
    plist1 = _get_all(promoter,genome, excite_offset, genesize)
    plist2 = _get_all(productprom,genome, excite_offset, genesize)#[:kwargs['nout']]
    alltogether = zip(plist1, [TF]*len(plist1)) + zip(plist2, [PR]*len(plist2))
    alltogether.sort(key=lambda x: x[0])
    alltogether = _filteroverlapping(alltogether, genesize)
    return ([i[0] for i in alltogether if i[1] == TF],
            [i[0] for i in alltogether if i[1] == PR])

def build_customproducts(signatures = INPUT_SIGNATURES):
    products = []
    i = 0
    for sig in signatures:
        products.append([i, sig, sig, sig, sig])
        i += 1
    return products


def iterate(arnet,samplerate, simtime, silentmode, simstep,delta,**kwargs):
    time = 1
    nump = arnet.numtf
    numsamples = 1
    numrec = arnet.numrec
    numeff = arnet.numeff

    while time <= simtime:
        for i in range(numrec):
            arnet.ccs[nump+i] = kwargs['inputs'][i]

        _update(arnet.proteins,arnet.ccs,arnet.eweights,arnet.iweights,
                delta, numtf = nump)
        #normalize ccs, ignoring effectors
        totparcels = arnet.ccs[:nump+numrec].tolist()
        #receptors ccs is not modified here
        arnet.ccs[:nump] /= sum(totparcels)
        #normalize outputs
        #FIXME: should this be done?
        for i in range(numeff):
            if arnet.ccs[-1-i] > 1.0:
                arnet.ccs[-1-i] = 1.0
            if arnet.ccs[-1-i] < .0:
                arnet.ccs[-1-i] = .0
        if numeff > 1:
            arnet.ccs[nump+numrec:] /= sum(arnet.ccs[nump+numrec:])

        if time % int(simtime*samplerate) == 0:
            log.debug('TIME: '+ str(time))
            arnet.updatehistory()
        time+=simstep
    return arnet


def _update(proteins, ccs, exciteweights, inhibitweights,delta,**kwargs):
    deltas = (_getSignalArray(ccs[:len(exciteweights)],exciteweights) -
              _getSignalArray(ccs[:len(inhibitweights)],inhibitweights))

    deltas *= delta
    #NOTE: previous output concentration shall not be accounted for
    numreg = exciteweights.shape[0]
    #NOTE: not touching receptors ccs
    deltas[:kwargs['numtf']] *= ccs[:kwargs['numtf']]
    ccs[:kwargs['numtf']] += deltas[:kwargs['numtf']]
    ccs[numreg:] += deltas[numreg:]

class ARNetwork(arn.ARNetwork):
    def __init__(self, gcode, config, **kwargs):
        self.code = gcode
        self.simtime = config.getint('default','simtime')
        prob = kwargs['problem']

        promfun = bindparams(config, buildpromlist)
        productsfun = bindparams(config, buildproducts)
        #TODO: add nout to the parameters file or bring the problem instance here?
        self.promlist,self.effectorproms = \
                        promfun(gcode,
                                productprom='11111111',
                                nout = prob.nout)
        self.proteins = productsfun( gcode, self.promlist)

        self.effectors=[]
        if self.effectorproms:
            self.effectors = productsfun(gcode,self.effectorproms)

        self.receptors= build_customproducts()
        self.receptorproms = [r[0] for r in self.receptors]

        pbindfun = bindparams(config, getbindings)
        weightsfun = bindparams(config, _getweights)

        self.numtf = len(self.proteins)
        if prob.nout == 1:
            self.numeff = len(self.effectors)
        else:
            self.numeff = min(len(self.effectors),prob.nout)

        self.effectorproms = self.effectorproms[:self.numeff]
        self.effectors = self.effectors[:self.numeff]
        self.numrec = min(len(self.receptorproms),prob.ninp)
        self.receptorproms = self.receptorproms[:self.numrec]
        self.receptors = self.receptors[:self.numrec]
        self.ccs = []
        if self.promlist or self.effectorproms:
            self.ccs = nparray([1.0/(self.numtf+self.numeff+self.numrec)]*
                               (self.numtf+self.numeff+self.numrec))
            self._initializehistory()
            self._initializebindings(pbindfun)
            self._initializeweights(weightsfun)
            for i in range(len(self.proteins)):
                self.proteins[i].append(self.ccs[i])

        self.simfun = bindparams(config,iterate)

        self.delta = config.getfloat('default','delta')
        log.debug(self.snapshot())
        self.output_idx = 0

    def _initializebindings(self, pbindfun):
        self.ebindings = pbindfun(0, self.proteins  + self.receptors +
                                  self.effectors)
        self.ibindings = pbindfun(1, self.proteins + self.receptors +
                                  self.effectors)

        #effectors do not regulate
        if self.effectors:
                self.ebindings = self.ebindings[:self.numtf+self.numrec,:]
                self.ibindings = self.ibindings[:self.numtf+self.numrec,:]

    def _initializeweights(self, weightsfun):
        self.eweights = weightsfun(self.ebindings)
        self.iweights = weightsfun(self.ibindings)

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
        self._initializehistory()

    def __len__(self):
        return len(self.proteins)+len(self.effectors)

    def __str__(self):
        return str(self.proteins)

    def simulate(self, *inputs):
        if not inputs:
            inputs = [.0]*self.numrec
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

    def printgraph(self):
        '''In order to be used as a Phenotype object'''
        displayARNresults(self.proteins, self.cchistory, temp=0)
        extralabels = ['R']*self.numrec + ['E']*self.numeff
        struct_prots = self.receptors + self.effectors
        hist = np.vstack((self.receptorhist,self.effectorhist))
        displayARNresults(struct_prots, hist,
                          temp = 1, extralabels = extralabels)
        return "ARN simulation chart was printed!"


    def snapshot(self):
        s = 'digraph best {\nordering = out;\n'
        shape = 'hexagon'
        labelidx = 0
        for outprom in self.effectorproms:
            s += '%i [label="%s",shape="hexagon"];\n' % (outprom, labelidx)
            for e,h,i in zip(self.ebindings[:,self.numtf+self.numrec+labelidx],
                             self.ibindings[:,self.numtf+self.numrec+labelidx],
                             range(self.ebindings.shape[0])):
                if e > 0:
                    s += '%i -> %i [dir=back,weight=%i];\n' % \
                             (outprom, (self.promlist + self.receptorproms)[i],e)
                if h > 0:
                    s += '%i -> %i [dir=back,style=dotted,weight=%i];\n' % \
                             (outprom, (self.promlist + self.receptorproms)[i],h)
            labelidx += 1

        for tf in self.promlist:
            s += '%i [label="%s"];\n' % (tf, labelidx )
            for e,h,i in zip(self.ebindings[:,labelidx-self.numeff],
                             self.ibindings[:,labelidx-self.numeff],
                             range(self.ebindings.shape[0])):
                if e > 0:
                    s += '%i -> %i [dir=back,weight=%i];\n' % \
                             (tf, (self.promlist + self.receptorproms)[i], e)
                if h > 0:
                    s += '%i -> %i [dir=back,style=dotted,weight=%i];\n' % \
                             (tf, (self.promlist + self.receptorproms)[i],h)
            labelidx += 1

        for rec in self.receptorproms:
            s += '%i [label="%s",shape="rectangle"];\n' % (rec, labelidx)
            labelidx += 1
        s += '}'
        return s

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
