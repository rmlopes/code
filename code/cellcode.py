import ConfigParser, os
import logging
import random
import numpy as np
import extendedarn as arn
from evodevo import Agent
from numpy import array as nparray
from functools import partial
from bitstring import *
from math import exp
import sys
from sys import stdout
from subprocess import call
from utils import *
from utils.bitstrutils import *
from extendedarn import *


log = logging.getLogger(__name__)

### Problem base to use with ReNCoDe
class CellProb:
        def __init__(self, evaluate, numins, numouts):
                self.eval_ = evaluate
                self.ninp = numins
                self.nout = numouts
                self.print_ = printcell

def printcell(arn):
        return arn.code.bin

class Cell(Agent):
        phenotype = None
        genotype = None
        fitness = None
        def __init__(self, config, gcode = None, parentfit = 1e4, **kwargs):
            Agent.__init__(self, parentfit)
            generator = arn.bindparams(config, arn.generatechromo)
            if gcode == None:
                    gcode = generator()

            self.genotype = arn.ARNetwork(gcode, config, **kwargs)
            #because now the phenotype is expressed at
            #evaluatiuon time
            self.phenotype = self.genotype
            #self.phenotype = arn.ARNetwork(gcode,config)
            while self.phenotype.numeff == 0 or self.phenotype.numrec == 0:
                gcode = generator()
                self.genotype = arn.ARNetwork(gcode, config, **kwargs)
                self.phenotype = self.genotype

            #initialize phenotype
            self.phenotype.nstepsim(config.getint('default','simtime'),
                                    *nparray([0.0,.0,.0]))
            self.fitness = 1e4

        def __str__(self):
            return "### Agent ###\n%s\n%s: %f" % (self.arn,self.circuit,
                                                  self.fitness)

        def setARN(self, arn):
                self.arn = arn

def plotindividual(arnet, **kwargs):
        displayARNresults(arnet.proteins, arnet.cchistory,
                          kwargs['samplerate'], temp=0,
                          figure=TFACTORS)
        extralabels = ['R']*arnet.numrec + ['E']*arnet.numeff
        struct_prots = arnet.receptors + arnet.effectors
        hist = np.vstack((arnet.receptorhist,arnet.effectorhist))
        displayARNresults(struct_prots, hist, kwargs['samplerate'],
                          temp = 1, extralabels = extralabels,
                          figure=STRUCTS)

def getbinaryoutput(arn, **kwargs):
        index = 0
        gradientsum = []
        leftlim =  -(1.0 / kwargs['samplerate'])
        output = np.sum(np.gradient(arn.effectorhist[index][leftlim:]))
        return 0 if output > 0 else 1

def getoutputp0p1(arn, **kwargs):
        index = 0
        leftlim =  -(1.0 / kwargs['samplerate'])
        g1 = np.sum(np.gradient(arn.effectorhist[index][leftlim:]))

        try:
                g2 = np.sum(np.gradient(arn.effectorhist[index + 1][leftlim:]))
        except IndexError: return 0

        return 1 if g2 >= g1 else 0

#TODO: move this inside the problem
def evaluatecircuit(phenotype, test = False, **kwargs):
        mapfun = getbinaryoutput
        try:
                mapfun = kwargs['mapfun']
        except KeyError: pass
        n = 3
        ok=0
        intinps = range(pow(2,n))
        if not test:
                intinps = intinps + intinps
        #random.shuffle(intinps)
        try:
                if kwargs['shuffle']:
                        random.shuffle(intinps)
        except KeyError: pass

        for i in intinps:
                inputs = BitStream(uint = i, length = n)
            #print inputs.bin
                normalized = nparray([float(inputs.bin[i])
                                      for i in range(n)])
                normalized *= .1
                phenotype.nstepsim(kwargs['simtime'],*normalized)
                out = mapfun(phenotype, **kwargs)
                        #print 'OUT: ', out
                if out == inputs[1+inputs[0]]:
                        ok += 1

        #print 'SILENT: ', kwargs['silentmode']
        if not kwargs['silentmode']:
            plotindividual(phenotype,**kwargs)
        return len(intinps) - ok

def evaluatewithreset(phenotype, test = False, **kwargs):
        mapfun = getbinaryoutput
        try:
                mapfun = kwargs['mapfun']
        except KeyError: pass
        n = 3
        ok=0
        intinps = range(pow(2,n))

        initstate = phenotype.ccs

        for i in intinps:
                inputs = BitStream(uint = i, length = n)
            #print inputs.bin
                normalized = nparray([float(inputs.bin[i])
                                      for i in range(n)])
                normalized *= .1
                phenotype.nstepsim(kwargs['simtime'],*normalized)
                out = mapfun(phenotype, **kwargs)
                        #print 'OUT: ', out
                if out == inputs[1+inputs[0]]:
                        ok += 1
                phenotype.reset(initstate)

        #print 'SILENT: ', kwargs['silentmode']
        if not kwargs['silentmode']:
            plotindividual(phenotype,**kwargs)
        return len(intinps) - ok

if __name__ == '__main__':
        arnconfigfile = '../configfiles/arnsim.cfg'
        log.setLevel(logging.DEBUG)
        cfg = ConfigParser.ConfigParser()
        cfg.readfp(open(arnconfigfile))

        proteins=[]
        nump = 0
        try:
            f = open(sys.argv[1], 'r')
            genome = BitStream(bin=f.readline())
            arnet = ARNetwork(genome, cfg)
        except:
            while nump < 4 or nump > 32 or numeff == 0 or not arnet.receptors:
                genome = BitStream(float=random.random(), length=32)
                for i in range(cfg.getint('default','initdm')):
                    genome = dm_event(genome,
                                      .02)

                    arnet = ARNetwork(genome, cfg, numinputs = 3,numoutputs=1)
                    nump = len(arnet.promlist)
                    numeff = len(arnet.effectors)

        offspring = None
        themother = Cell(cfg,genome )
        eval_ = bindparams(cfg, evaluatecircuit)
        plot_ = bindparams(cfg, plotindividual)

        pop = [(themother, eval_(themother)),
               (None,0)]

        while pop[0][1] > 0:
                offspring = Cell(cfg, bitflipmutation(
                        pop[0][0].genotype.code,.01))
                while not offspring.phenotype.effectors:
                    offspring = Cell(cfg, bitflipmutation(
                            pop[0][0].genotype,.01))
                pop[1] = (offspring,eval_(offspring))
                pop.sort(key = lambda x: x[1])
                print pop

                #for p in arnet.proteins: print p


        f = open('genome.save','w')
        f.write(pop[0][0].genotype.code.bin)
        f.close
        #print genome.bin
        plot_(pop[0][0].genotype)

def buildcircuit(agent, problem, **kwargs):
        #DELETE THIS: NO NEED TO BUILD A CIRCUIT, BEHAVIOR WILL BE THE
        #DEVELOPMENT OF THE REGULATORY NETWORK INTERACTING WITH THE
        #ENVIRONMENT THROUGH RECEPTORS AND EFFECTORS, DURING FITNESS EVALUATION
        """Returns the circuit to be fed into the evaluation function"""
        arn = agent.genotype
        if not arn.promlist:
                return []

        arn.nstepsim()

        #map signatures to bitwise operators
        #cumsum of ccs to probabilities

        return circuit