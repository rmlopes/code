import logging
from functools import partial
from bitstring import BitStream
import math
from arn import *
import rencode
from evodevo import Problem, Agent
from utils import *
from utils.mathlogic import *
from utils.bitstrutils import *
import copy

log = logging.getLogger(__name__)

### Problem base to use with ReNCoDe
class Prob():
        def __init__(self, evaluate, grammar, startcodon, **kwargs):
            try:
                pp = kwargs['printf']
            except KeyError:
                self.print_ = (lambda x: x.circuit)
            self.eval_ = evaluate
            self.grammar = grammar
            self.start = startcodon


### Agent model to use with this CoDe module
class GEARNetAgent(Agent):
        genotype = None
        phenotype = None
        fitness = None
        def __init__(self, config, problem, gcode = None, parent = None):
                Agent.__init__(self, parent)
                generator = bindparams(config, self.generate)
                if gcode == None:
                    numtfs = 0
                    while not numtfs:
                        arnet = ARNetwork(generator(),config, problem=problem)
                        numtfs = arnet.numtf
                else:
                    arnet = ARNetwork(gcode,config, problem=problem)
                self.genotype = arnet
                self.genotype.simulate()
                #print arnet.promlist
                self.phenotype = Phenotype(arnet,problem)
                self.problem = problem
                self.fitness = 1e9

        def __str__(self):
                return "### Agent ###\n%s\n%s: %f" % (self.arn,self.phenotype,
                                                      self.fitness)

        def pickled(self):
            return self.genotype.code.bin

        def print_(self):
            return self.phenotype.circuit

class DMAgent(GEARNetAgent):
        def __init__(self, config, problem, gcode = None, parent = None):
                self.generate = generatechromo
                GEARNetAgent.__init__(self, config, problem, gcode, parent)

class RndAgent(GEARNetAgent):
        def __init__(self, config, problem, gcode = None, parent = None):
            self.generate = partial(generatechromo_rnd,
                    genomesize = 32 * pow(2,config.getint('default','initdm')))
            GEARNetAgent.__init__(self, config, problem, gcode, parent)

def replace( intnum, key, g):
    i = intnum % len(g[key])
    log.debug("Returning rule " + str(i) + " for key " +
              key + " -> "+str( g[key][i]))
    return g[key][i]

class Phenotype:
    def __init__(self, arnet, problem):
        self.problem = problem
        self.arnet = arnet
        self.int_sequence = self.extract_int_seq()
        #print self.int_sequence
        self.circuit = reduce(lambda l,m: l+m,
                              self.translate(copy.deepcopy(self.int_sequence),
                                             problem.grammar,[problem.start],self.int_sequence))
        print self.circuit

    def __len__(self):
        return len(self.int_sequence)

    def __call__(self, *inputs):
        #log.info(self.circuit)
        if not self.isvalid():
            return 1e6
        return eval(self.circuit)

    def isvalid(self):
        for key in self.problem.grammar.iterkeys():
            if self.circuit.find(key) >= 0:
                return False
        return True

    def __eq__(self,other):
        return self.circuit == other.circuit

    def printgraph(self):
        return self.circuit

    def extract_int_seq(self):
        sz = 8
        functions = [[p[0],
                      [p[4][i*sz:sz+i*sz].uint for i in range(len(p[4])/sz)],
                      #[p[4][i*sz:sz+i*sz].count(1) for i in range(len(p[4])/sz)],
                      p[5]] for p in self.arnet.proteins]
        functions.sort(key=lambda function: function[2],reverse=True)
        log.debug(functions)
        return [el for f in functions for el in f[1]]

    def translate(self, int_sequence, grm, tree ,originalseq=[], wrap=1):
        if len(tree) == 0:
            return []
        #print "top el: ", tree[0]
        #print grm.iterkeys()
        if tree[0] in grm.iterkeys():
            if len(int_sequence) == 0:
                if wrap:
                    int_sequence.extend(originalseq)
                    wrap -= 1
                else:
                    return tree
            x = int_sequence.pop(0)
            #print tree
            tree = replace(x,tree[0],grm) + tree[1:]
            return self.translate(int_sequence, grm,
                                  tree,originalseq,wrap)
        else:
            return [tree[0]] + self.translate(int_sequence, grm,
                                            tree[1:],originalseq,wrap)



###########################################################################
### Test                                                                ###
###########################################################################

if __name__ == '__main__':
        import ConfigParser
        import random
        from bitstring import BitStream
        from bitstringutils import dm_event
        from arn import ARNetwork
        from evodevo import Problem

        log.setLevel(logging.DEBUG)
        cfg = ConfigParser.ConfigParser()
        cfg.readfp(open('test-rencode.cfg'))
        arncfg = ConfigParser.ConfigParser()
        arncfg.readfp(open('arnsim-rencode2.cfg'))
        proteins=[]
        nump = 0
        while nump < 3:
                genome = BitStream(float=random.random(), length=32)
                for i in range(5):
                        genome = dm_event(genome,
                                          .02)

                arnet = ARNetwork(genome, arncfg)
                nump = len(arnet.promlist)
        for p in arnet.proteins: print p
        prob = Problem(defaultnodemap,regressionfun, None)
        circuit = buildcircuit(arnet, prob,False)
        _printcircuit(circuit)
        print printdotcircuit(circuit, prob.labels)
