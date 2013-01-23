import logging
from functools import partial
from bitstring import BitStream
import math
import arn
import rencode
from arnmiguel import *
from evodevo import Problem, Agent
from utils import *
from utils.bitstrutils import *

log = logging.getLogger(__name__)

def printrencode2(phenotype, **kwargs):
    return rencode.printdotcircuit(
            phenotype.circuits[phenotype.arnet.output_idx],**kwargs)

### Agent model to use with this CoDe module
class ARNGPAgent(Agent):
        genotype = None
        phenotype = None
        fitness = None
        def __init__(self, config, problem, gcode = None, parentfit = 1e4):
                Agent.__init__(self, parentfit)
                generator = bindparams(config, self.generate)
                if gcode == None:
                    numprods = 0
                    while not numprods:
                        arnet = ARNetwork(generator(),config, problem=problem)
                        numprods = arnet.numeff
                else:
                    arnet = ARNetwork(gcode,config, problem=problem)
                self.genotype = arnet
                self.phenotype = Phenotype(self.genotype,problem)
                self.problem = problem
                self.fitness = 1e9

        def __str__(self):
                return "### Agent ###\n%s\n%s: %f" % (self.arn,self.phenotype,
                                                      self.fitness)

        def pickled(self):
            return (self.genotype.code.bin,self.phenotype.output_idx)

class DMAgent(ARNGPAgent):
        def __init__(self, config, problem, gcode = None, parentfit = 1e4):
                self.generate = arn.generatechromo
                ARNGPAgent.__init__(self, config, problem, gcode, parentfit)

class RndAgent(ARNGPAgent):
        def __init__(self, config, problem, gcode = None, parentfit = 1e4):
            self.generate = partial(arn.generatechromo_rnd,
                    genomesize = 32 * pow(2,config.getint('default','initdm')))
            ARNGPAgent.__init__(self, config, problem, gcode, parentfit)

class Phenotype:
    def __init__(self, arnet, problem):
        self.problem = problem
        self.products = arnet.effectorproms
        self.circuits = []
        self.graph = arnet.ebindings - arnet.ibindings
        self.arnet = arnet
        if problem.nout == 1:
            for i in range(len(self.products)):
                arnet.output_idx = i
                queue = self.getinputs(arnet.numtf+arnet.numrec+i)
                #print "effector inputs: ", queue
                fset = problem.funs
                if not queue:
                        fset = problem.terms
                circuit = [(arnet.effectors[i][0],
                            self.problem.nodemap_(arnet.effectors[i][4],fset),
                            [q[0] for q in queue])]
                self.circuits.append(self.buildcircuit(queue, circuit))
                #print self.circuits[-1]
                #print problem.print_(self)
        else:
            queue=[]
            circuit=[]
            for i in range(arnet.numeff):
                inps = self.getinputs(arnet.numtf+arnet.numrec+i)
                fset = problem.funs
                if not inps:
                        fset = problem.terms
                circuit.append((arnet.effectors[i][0],
                                self.problem.nodemap_(arnet.effectors[i][4],
                                                      fset),
                                [inp[0] for inp in inps]))
                queue.extend(inps)
            queue.sort(key=lambda x: x[1], reverse=True)
            self.circuits.append(self.buildcircuit(queue,circuit))

        #map receptors to inputs
        #map tfs to functions
        #map output(s) to functions
        #several options to get the connections
        #start by using match strength using positive/negative to
        #give the direction

    def getcircuit(self, index=0):
        return self.circuits[index]

    def buildcircuit(self, queue, circuit=[], force_inputs=False):
        if not queue:
            if len(circuit) > 1 or force_inputs:
                circuit.extend(self.getinputnodes())
            return circuit

        pnext = queue.pop(0)
        blacklist = [c[0] for c in circuit]
        if pnext[0] in blacklist:
            return self.buildcircuit(queue, circuit)
            #break;
        elif pnext[0] in self.arnet.receptorproms:
            return self.buildcircuit(queue, circuit, force_inputs=True)
            #break;
        else:
            blacklist.append(pnext[0])
            index = self.arnet.promlist.index(pnext[0])
            signature = self.arnet.proteins[index][4]
            fun = self.problem.nodemap_(signature,self.problem.funs)
            inputs = self.getinputs(index)
            inputs = filter(lambda x: x[0] not in blacklist, inputs)
            ar = self.problem.arity[fun]
            if ar > 0:
                    inputs = inputs[:ar]
            if not inputs:
                #print signature, math.tanh(signature.int)
                fun = str(math.tanh(signature.int))
            circuit.append((pnext[0],fun,[i[0] for i in inputs]))

            queue.extend(inputs)
            queue.sort(key=lambda x: x[1], reverse=True)
            return self.buildcircuit(queue,circuit)

    def getinputs(self, index):
        inputs=[]
        rlist = self.arnet.receptorproms# [r[0] for r in self.arnet.receptors]
        pmap = zip(self.arnet.promlist+rlist,self.graph[:,index])
        for p,w in pmap:
                if w > 0:
                        inputs.append((p,w))
        inputs.sort(key=lambda x: x[1], reverse = True)
        return inputs

    def getinputnodes(self):
        return map(lambda x: (x[0],
                              self.problem.terms[x[0]],
                              []),
                   self.arnet.receptors)


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
