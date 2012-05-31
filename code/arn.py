import ConfigParser, os
import logging
import random
import numpy as np      
from numpy import array as nparray
from functools import partial
from bitstring import *
from math import exp
from sys import stdout
from subprocess import call
from utils import *
from utils.bitstrutils import *

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
    while not 30 > valid >= 4:
        genome = BitStream(float=random.random(),length=32);
        for i in range(0,initdm):
            genome = dm_event(genome, mutratedm)
        promlist = buildpromlist(genome, excite_offset, genesize, promoter)
        valid = len(promlist)

    return genome

def displayARNresults(proteins, ccs, step):
    log.warning('Plotting simulation results for ' + 
                str(len(proteins)) + ' genes/proteins')
    plt.clf()
    xx=[i*step for i in range(len(ccs[0]))]
    for i in range(len(proteins)):
        plt.plot(xx, ccs[i],label=str(proteins[i][0]))
    plt.savefig('ccoutput.png')
    call(["open", "ccoutput.png"])

def buildpromlist(genome, excite_offset, genesize, promoter,**kwargs):
    gene_index = genome.findall(BitStream(bin=promoter))
    promsize = len(promoter)
    promlist = filter( lambda index: 
                       int(excite_offset) <= index <  (genome.length-(int(genesize)+promsize )), 
                       gene_index)
    proms = reduce(lambda indxlst, indx: 
                   indxlst + [indx] if(indx-indxlst[-1] >= 32) else indxlst, 
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
   
def iterate(proteins, ccs, excite_weights, inhibit_weights,
            samplerate, simtime, silentmode, simstep,delta,**kwargs):
    time = 0
    cchistory=nparray(ccs)
    while time < simtime:
        _update(proteins,ccs,excite_weights,inhibit_weights,delta)
        if(not(silentmode) and 
           (time % (simtime*samplerate) == 0)):
            log.debug('TIME: '+ str(time))
            for p in proteins:
                cchistory = np.column_stack((cchistory,ccs))                        
        time+=simstep
            
    if not silentmode:
        displayARNresults(proteins, cchistory,simstep)        
    for i in range(len(proteins)):
        proteins[i].append(ccs[i])

    
def _update(proteins, ccs, exciteweights, inhibitweights,delta):
    deltas = (_getSignalArray(ccs,exciteweights) - 
              _getSignalArray(ccs,inhibitweights))
    deltas *= delta
    deltas *= ccs
    total = sum(ccs)+sum(deltas)
    ccs+=deltas
    ccs/=total

def _getSignalArray(ccs, weightstable):
    return 1.0/len(ccs) * np.dot(ccs,weightstable)            
   
def _getprotein(idx, code, bind_size, gene_size, protein_size):
    signature = BitStream(bin=applymajority(code[bind_size*3:bind_size*3+gene_size], protein_size))
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
    def __init__(self, gcode, config):
        self.code = gcode
        self.simtime = config.getint('default','simtime')

        promfun = bindparams(config, buildpromlist)
        productsfun = bindparams(config, buildproducts)
        self.promlist = promfun(gcode)
        self.proteins = productsfun( gcode, self.promlist) 
        
        pbindfun = bindparams(config, getbindings)
        weightsfun = bindparams(config, _getweights)
        nump = len(self.proteins)
        if self.promlist:
            self.ebindings = pbindfun(0,self.proteins)
            self.ibindings = pbindfun(1,self.proteins)
            self.eweights = weightsfun(self.ebindings)
            self.iweights = weightsfun(self.ibindings)
            self.ccs=nparray([1.0/nump]*nump)
        self.simfun = bindparams(config,iterate)

    def __str__(self):
        return str(self.proteins)
           
    def simulate(self):
        if self.simtime > 0:
            self.simfun(self.proteins, self.ccs, self.eweights, self.iweights)

###########################################################################
### Test                                                                ###
###########################################################################

if __name__ == '__main__': 
    arnconfigfile = 'arnsim.cfg'
    log.setLevel(logging.DEBUG)
    cfg = ConfigParser.ConfigParser()
    cfg.readfp(open('test.cfg'))
    proteins=[]
    nump = 0
    while nump < 3:
        genome = BitStream(float=random.random(), length=32)
        for i in range(cfg.getint('default','initdm')):
            genome = dm_event(genome, 
                              .02)
            
        arnet = ARNetwork(genome, 'arnsim.cfg')
        nump = len(arnet.promlist)
        
    for p in arnet.proteins: print p
    arnet.simulate()



###########################################################################
### Other Helpers                                                       ###
###########################################################################
#deprecated
def generatechromoepi(init_dm, dm_mutrate,**bindargs):
    valid = False
    #initdm = random.gauss(float(init_dm),1.0)
    while not 30 > valid >= 4:
        genome = BitStream(float=random.random(),length=32);
        for i in range(0,int(init_dm)):
            genome = dm_event(genome, dm_mutrate)
        promlist = buildpromlist(genome, bindargs['excite_offset'], 
                                 bindargs['genesize'], bindargs['promoter'])
        valid = len(promlist)

    proteins = buildproducts(genome, promlist, 
                             bindargs['excite_offset'], 
                             len(bindargs['promoter']),
                             bindargs['genesize'], 
                             bindargs['bindingsize'], 
                             bindargs['proteinsize'])
    cromatines = [1.0/float(valid)]*valid
    return (genome,dict(zip(promlist,proteins)),
            dict(zip(promlist,cromatines)),10000,0)


def buildpromlistEPI(genome, excite_offset, genesize, promoter):
    gene_index = genome.findall(promoter)
    promsize = len(promoter)
    promlist = filter( lambda index: 
                       excite_offset <= index <  (genome.length-(genesize+promsize )), 
                       gene_index)
    proms = reduce(lambda indxlst, indx: 
                   indxlst + [indx] if(indx-indxlst[-1] >= genesize+3*32) else indxlst, 
                   promlist,
                   [0])
    return proms[1:]

#individual is (genome,proteindict,[epigenome])
#deprecated
def simulate(individual, bindingsize, proteinsize, genesize, promoter,
             excite_offset, match_threshold, beta, delta, 
             samplerate, simtime, simstep, silentmode):
    #genome,proteins,epig,lf,inactive = individual
    #MODIFIED
    #promlist = epig.keys()#buildpromlist(genome, excite_offset,genesize,promoter)
    promlist = individual[1].keys()
    proteins = individual[1].values()
    threshold = match_threshold
   
    if promlist and simtime > 0:
            excite_weights = getweights(0, proteins, threshold, 
                                        bindingsize, beta)
            inhibit_weights = getweights(1, proteins, threshold, 
                                         bindingsize, beta)
            initccs = [1.0/len(proteins)]*len(proteins)
            ccs=nparray(initccs)
            iterate(proteins, ccs, excite_weights, inhibit_weights,
                    samplerate,simtime, silentmode,simstep, delta)

    return proteins

#organized in columns for the target equation
#deprecated
def getweights(bindtype, proteins, threshold, bindsize,beta):
    bindings = getbindings(bindtype, proteins, threshold, bindsize)
    bindings -= bindsize
    bindings *= beta
    return np.exp(bindings)

