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

### Problem base to use with CellCoDe
class CellProb:
        def __init__(self, evaluate, numins, numouts):
                self.eval_ = evaluate
                self.ninp = numins
                self.nout = numouts
                self.print_ = printcell

def printcell(cell):
        return (cell.input_weights, cell.hidden_weights, cell.output_weights)

def nn(inputs, input_weights, hidden_weights, output_weights):
        if len(inputs) > input_weights.shape[0]:
            inputs = inputs[:input_weights.shape[0]]
        hidenum = hidden_weights.shape[0]
        outnum = output_weights.shape[1]
        #print inputs.shape, ' ', input_weights.shape, ' ', hidden_weights.shape, ' ', output_weights.shape

        K,Z,Y = [],[],[]
        for i in range(hidenum):
            K.append(sigmoid(inputs, input_weights[:,i]))
        for i in range(hidenum):
            Z.append(sigmoid(K, hidden_weights[:,i]))
        for i in range(outnum):
            Y.append(sigmoid(Z,output_weights[:,i]))
        return Y

def sigmoid(inputs,weights):
        #print len(inputs), ' ', len(weights)
        return 1.0 / (1.0 + exp(-np.dot(inputs, weights)))

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
            while (self.genotype.numeff == 0 or
                   self.genotype.numrec == 0 or
                   self.genotype.numtf == 0):
                gcode = generator()
                self.genotype = arn.ARNetwork(gcode, config, **kwargs)

            #initialize phenotype
            weights = self.genotype.eweights - self.genotype.iweights
            nump = self.genotype.numtf
            numinp = self.genotype.numrec
            numout = self.genotype.numeff
            class Pheno(): pass
            self.phenotype = Pheno()
            self.phenotype.input_weights = weights[nump:nump+numinp,:nump]
            #print self.phenotype.input_weights
            self.phenotype.hidden_weights = weights[:nump,:nump]
            #print self.phenotype.hidden_weights
            self.phenotype.output_weights = weights[:nump,nump+numinp:]
            #print self.phenotype.output_weights
            self.fitness = 1e9

        def __str__(self):
            return "### Agent ###\n%s\n%s: %f" % (self.arn,self.circuit,
                                                  self.fitness)

        def setARN(self, arn):
                self.arn = arn


#TODO: move this inside the problem
def evaluatecircuit(phenotype, test = False, **kwargs):
        n = 3
        ok=0
        intinps = range(pow(2,n))
        #if not test:
         #       intinps = intinps + intinps
        #random.shuffle(intinps)
        try:
                if kwargs['shuffle']:
                        random.shuffle(intinps)
        except KeyError: pass

        for i in intinps:
                inputs = BitStream(uint = i, length = n)
            #print inputs.bin
                #not normalized only floated
                normalized = nparray([float(inputs.bin[i])
                                      for i in range(n)])
                out = nn(normalized,
                         phenotype.input_weights,
                         phenotype.hidden_weights,
                         phenotype.output_weights)
                print 'OUT: ', out
                out = (0 if out[0] < .5 else 1)
                print 'OUT: ', out
                if out == inputs[1+inputs[0]]:
                        ok += 1

        #print 'SILENT: ', kwargs['silentmode']
        #if not kwargs['silentmode']:
         #   plotindividual(phenotype,**kwargs)
        return len(intinps) - ok

if __name__ == '__main__':
        arnconfigfile = '../configfiles/arnsim.cfg'
        log.setLevel(logging.DEBUG)
        cfg = ConfigParser.ConfigParser()
        cfg.readfp(open(arnconfigfile))

        proteins=[]
        p = CellProb(evaluatecircuit, 3, 1)
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

                    arnet = ARNetwork(genome, cfg, problem = p)
                    nump = len(arnet.promlist)
                    numeff = len(arnet.effectors)

        offspring = None
        themother = Cell(cfg,genome,problem = p )
        eval_ = bindparams(cfg, evaluatecircuit)
        #plot_ = bindparams(cfg, plotindividual)

        pop = [(themother, eval_(themother.phenotype)),
               (None,0)]

        while pop[0][1] > 0:
                offspring = Cell(cfg, bitflipmutation(
                        pop[0][0].genotype.code,.01), problem = p)
                while not offspring.genotype.effectors:
                    offspring = Cell(cfg, bitflipmutation(
                            pop[0][0].genotype,.01))
                pop[1] = (offspring,eval_(offspring.phenotype))
                pop.sort(key = lambda x: x[1])
                print pop

                #for p in arnet.proteins: print p


        f = open('genome.save','w')
        f.write(pop[0][0].genotype.code.bin)
        f.close
        #print genome.bin
        #plot_(pop[0][0].genotype)
