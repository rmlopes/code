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
try:
    import matplotlib.pyplot as plt
    from matplotlib import collections, legend
except ImportError:
    logging.warning('matplotlib not found.')
    logging.warning('use of silentmode = 0 will result in an error...')

log = logging.getLogger(__name__)


def bindparams(config,fun):
    '''Binds the ARN configuration file parameters to a function.'''
    return partial(fun,
                   bindingsize = config.getint('default','bindingsize'),
                   proteinsize = config.getint('default','proteinsize'),
                   genesize = config.getint('default','genesize'),
                   promoter = config.get('default','promoter'),
                   excite_offset = config.getint('default','excite_offset'),
                   match_threshold = config.getint('default','match_threshold'),
                   beta = config.getfloat('default','beta'),
                   delta = config.getfloat('default','delta'),
                   samplerate = config.getfloat('default','samplerate'),
                   simtime = config.getint('default','simtime'),
                   simstep = config.getint('default','simstep'),
                   silentmode = config.getboolean('default','silentmode'),
                   initdm = config.getint('default','initdm'),
                   mutratedm = config.getfloat('default','mutratedm'))

def generatechromo(initdm, mutratedm, genesize,
                   promoter, excite_offset, **bindargs):
    '''
    Default function to generate an ARN chromosome.
    To be used with bindparams.
    '''
    valid = False
    while not 48 > valid >= 4:
        genome = BitStream(float=random.random(),length=32);
        for i in range(0,initdm):
            genome = dm_event(genome, mutratedm)
        promlist = buildpromlist(genome, excite_offset, genesize, promoter)
        valid = len(promlist)

    return genome

TFACTORS = 0
STRUCTS = 1

def displayARNresults(proteins, ccs,
                      samplerate=.001, temp = 0,extralabels=None,**kwargs):
    log.warning('Plotting simulation results for ' +
                str(len(proteins)) + ' genes/proteins')
    #plt.figure(kwargs['figure'])
    plt.clf()
    xx = nparray(range(ccs.shape[1]))
    if extralabels:
        for i in range(len(proteins)):
            plt.plot(xx, ccs[i],label="%s%i"%(extralabels[i],proteins[i][0],))
        plt.legend()

    else:
        for i in range(len(proteins)):
            plt.plot(xx, ccs[i])
    plt.savefig('ccoutput_' + str(temp) + '.png')
    #plt.show()
    #call(["open",'ccoutput_' + str(temp) + '.png'])

def buildpromlist(genome, excite_offset, genesize, promoter,**kwargs):
    gene_index = genome.findall(BitStream(bin=promoter))
    promsize = len(promoter)
    promlist = filter( lambda index:
                       int(excite_offset) <= index <  (genome.length-
                                                       (int(genesize)+
                                                        promsize )),
                       gene_index)
    proms = reduce(lambda indxlst, indx:
                   indxlst + [indx] if(indx-indxlst[-1] >= promsize
                             + (32-promsize)) else indxlst,
                   promlist,
                   [0])
    return proms[1:]

def buildproducts(genome, promlist, excite_offset, promoter,
                  genesize, bindingsize, proteinsize, **kwargs):
     log.debug("Building ARN with " + str(len(promlist)) + " genes")
    #each protein is
    #[protein_index(=prom_index), e-bind, h-bind,
    # bind-signature, function-signature ]
     proteins = list()
     for pidx in promlist:
         proteins.append(_getprotein(pidx,
                                     genome[pidx-excite_offset:pidx+genesize+len(promoter)],
                                     bindingsize,
                                     genesize,
                                     proteinsize))
     return proteins


#organized in columns for the target equation
def getbindings(bindtype, proteins, match_threshold,**kwargs):
    return nparray([[XORmatching(p[3],otherps[1+bindtype],match_threshold)
                         for otherps in proteins]
                        for p in proteins],dtype=float);

def iterate(arnet,samplerate, simtime, silentmode, simstep,delta,**kwargs):
    #s = clock()
    time = 0
    nump = arnet.numtf
    numsamples = 1
    numrec = arnet.numrec
    numeff = arnet.numeff

    #print kwargs['inputs']
    while time < simtime:
        for i in range(numrec):
            arnet.ccs[nump+i] = kwargs['inputs'][i]

        _update(arnet.proteins,arnet.ccs,arnet.eweights,arnet.iweights,delta)
        #normalize ccs, ignoring receptors
        #print arnet.ccs[:nump]
        totparcels = (arnet.ccs[:nump].tolist() +
                      arnet.ccs[nump+numrec:].tolist())
        arnet.ccs /= sum(totparcels)
        #enforce again input/receptor cc
        for i in range(numrec):
            arnet.ccs[nump+i] = kwargs['inputs'][i]

        if time % int(simtime*samplerate) == 0:
            log.debug('TIME: '+ str(time))
            arnet.cchistory = np.column_stack((arnet.cchistory,
                                               arnet.ccs[:nump]))
            arnet.receptorhist = np.column_stack((arnet.receptorhist,
                                                  arnet.ccs[nump:nump+numrec]))
            arnet.effectorhist = np.column_stack((arnet.effectorhist,
                                                  arnet.ccs[nump+numrec:]))
        time+=simstep

    #print 'Elapsed time: %f sec.' % (clock()-s,)
    return arnet


def _update(proteins, ccs, exciteweights, inhibitweights,delta):
    #print ccs.shape
    deltas = (_getSignalArray(ccs[:len(exciteweights)],exciteweights) -
              _getSignalArray(ccs[:len(inhibitweights)],inhibitweights))
    #print deltas.shape
    deltas *= delta
    deltas *= ccs
    #total = sum(ccs)+sum(deltas)
    ccs+=deltas
    #ccs/=total

def _getSignalArray(ccs, weightstable):
    return 1.0/len(ccs) * np.dot(ccs,weightstable)

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
                               (self.numtf+self.numeff+self.numrec))
            self._initializehistory()
            self._initializebindings(pbindfun)
            self._initializeweights(weightsfun)
            for i in range(len(self.proteins)):
                self.proteins[i].append(self.ccs[i])

        self.simfun = bindparams(config,iterate)
        self.delta = config.getfloat('default','delta')

    def _initializebindings(self, pbindfun):
        self.ebindings = pbindfun(0, self.proteins  + self.receptors +
                                  self.effectors)
        self.ibindings = pbindfun(1, self.proteins + self.receptors +
                                  self.effectors)
        # testing:non-dummy effectors also regulate
        #effectors and dummy receptors do not regulate
        #if self.effectors or self.receptors:
         #       self.ebindings = self.ebindings[:self.numtf+self.numrec,:]
                #print 'ebinds: ', self.ebindings.shape
          #      self.ibindings = self.ibindings[:self.numtf+self.numrec,:]
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

    def reset(self, cc_state = None):
        self._initializehistory()
        if len(cc_state)==0:
            self.ccs = nparray([1.0/(self.numtf+self.numeff+self.numrec)]*
                               (self.numtf+self.numeff+self.numrec))
        else:
            self.ccs = nparray(cc_state[:])


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
        arnconfigfile = '../configfiles/arnsimlong.cfg'
        class Problem:
            pass
        p = Problem()
        p.ninp=3
        p.nout = 1
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
        f = open('genome.save','w')
        f.write(genome.bin)
        f.close
        #print genome.bin
