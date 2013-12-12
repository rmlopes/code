import logging
from functools import partial
from bitstring import BitStream
import math
import arn
import rencode
from rencode import evaluatecircuit, nnlikefun, mergefun,P
from arnmiguel import *
from evodevo import Problem, Agent
from utils import *
from utils.bitstrutils import *
from utils.mathlogic import *

log = logging.getLogger(__name__)

#evaluatecircuit = partial(rencode.evaluatecircuit,
 #                         nout = 1)

def printrencode2(phenotype, **kwargs):
    return rencode.printdotcircuit(
            phenotype.circuits[phenotype.output_idx],**kwargs)


def printmultiplecircuit(phenotype, labels=None, arnet = None):
    circuit = phenotype.getcircuit(phenotype.output_idx)
    if not arnet:
        arnet = phenotype.arnet
    s = 'digraph best {\nordering = out;\n'
    for c in circuit:
        shape='oval'
        if c[0] in arnet.effectorproms:
            shape='hexagon'
        elif c[0] in arnet.receptorproms:
            shape='rectangle'
        s += '%i [label="%s",shape="%s"];\n' % (c[0], c[1], shape)
                #else labels[c[1]])
        for inp in c[2]:
            aux = "dir=back"
            if inp < 0:
                aux += ",style=dotted"
            s += '%i -> %i [%s];\n' % (c[0],abs(inp),aux)
    s += '}'
    return s

### Agent model to use with this CoDe module
class ARNGPAgent(Agent):
        genotype = None
        phenotype = None
        fitness = None
        def __init__(self, config, problem, gcode = None, parent = None):
                Agent.__init__(self, parent)
                generator = bindparams(config, self.generate)
                if gcode == None:
                    numprods = 0
                    numtfs = 0
                    while not numprods or not numtfs:
                        arnet = ARNetwork(generator(),config, problem=problem)
                        numprods = arnet.numeff
                        numtfs = arnet.numtf
                else:
                    arnet = ARNetwork(gcode,config, problem=problem)
                self.genotype = arnet
                #print arnet.promlist
                self.phenotype = Phenotype(arnet,problem)
                self.problem = problem
                self.fitness = 1e9

        def __str__(self):
                return "### Agent ###\n%s\n%s: %f" % (self.arn,self.phenotype,
                                                      self.fitness)

        def pickled(self):
            return (self.genotype.code.bin,self.phenotype.output_idx)

        def print_(self):
            return printmultiplecircuit(self.phenotype)

class DMAgent(ARNGPAgent):
        def __init__(self, config, problem, gcode = None, parent = None):
                self.generate = arn.generatechromo
                ARNGPAgent.__init__(self, config, problem, gcode, parent)

class RndAgent(ARNGPAgent):
        def __init__(self, config, problem, gcode = None, parent = None):
            self.generate = partial(arn.generatechromo_rnd,
                    genomesize = 32 * pow(2,config.getint('default','initdm')))
            ARNGPAgent.__init__(self, config, problem, gcode, parent)

class Phenotype:
    def __init__(self, arnet, problem):
        self.funskel = nnlikefun
        self.problem = problem
        self.products = arnet.effectorproms
        self.circuits = []
        self.graph = arnet.ebindings - arnet.ibindings
        self.arnet = arnet
        self.output_idx = 0
        if problem.nout == 1:
            for i in range(len(self.products)):
                self.output_idx = i
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
            outputs = self.selectoutputs(problem.nout, problem)
            for o in outputs:
                #slice-out the index of the output
                circuit.append(o[:3])
                queue.extend(self.getinputs(arnet.numtf+arnet.numrec+o[3]))
            #for i in range(arnet.numeff):
             #   inps = self.getinputs(arnet.numtf+arnet.numrec+i)
              #  fset = problem.funs
              #  if not inps:
              #          fset = problem.terms
              #  circuit.append((arnet.effectors[i][0],
              #                  self.problem.nodemap_(arnet.effectors[i][4],
              #                                        fset),
              #                  [inp[0] for inp in inps]))
              #  queue.extend(inps)
            queue.sort(key=lambda x: x[1], reverse=True)
            #print "Building circuit with this: "
            #print self.printgraph()
            self.circuits.append(self.buildcircuit(queue,circuit))
            #self.circuitstr = map(lambda x: compile(x,'<string>','eval'),
             #                     evaluatecircuit(self.circuits[0], mergefun, dict(),
              #                                    *problem.terms, nout=problem.nout))
            #self.replaceinp = re.compile('inputs')


        #map receptors to inputs
        #map tfs to functions
        #map output(s) to functions
        #several options to get the connections
        #start by using match strength using positive/negative to
        #give the direction
    def __len__(self, index = 0):
        return len(self.circuits[index])

    def __call__(self, *inputs, **kwargs):
     #   #print self.circuitstr
        #vardict = dict(globals())
      #  globals()['inputs']=inputs
        #print vardict
        #return map(lambda x: eval(x,globals(),{'inputs':inputs}), self.circuitstr)
        #print r
        #return r
        return evaluatecircuit(self.getcircuit(self.output_idx),
                               self.funskel,dict(), *inputs, nout=self.problem.nout)

    def __eq__(self, other):
        if len(self.circuits) != len(other.circuits):
            return False
        for i in range(len(self.circuits)):
            if self.circuits[i] != other.circuits[i]:
                return False
        return True

    def selectoutputs(self,n, problem):
        '''Selects the output nodes by number of input edges.'''
        temp = list()
        for i in range(self.arnet.numeff):
            inps = self.getinputs(self.arnet.numtf+self.arnet.numrec+i)
            fset = problem.funs
            if not inps:
                fset = problem.terms
            temp.append((self.arnet.effectors[i][0],
                         self.problem.nodemap_(self.arnet.effectors[i][4],
                                               fset),
                         [inp[0] for inp in inps],
                         i))
        #temp.sort(key=lambda x: len(x[2]), reverse=True)
        return temp[:n]

    def getcircuit(self, index=0):
        return self.circuits[index]

    def buildcircuit(self, queue, circuit=[], force_inputs=False):
        if not queue:
            #if len(circuit) > 1 or force_inputs:
            circuit.extend(self.getinputnodes())
            return circuit

        pnext = queue.pop(0)
        blacklist = [c[0] for c in circuit]
        if pnext[0] in blacklist or pnext[0] in self.arnet.receptorproms:
            return self.buildcircuit(queue, circuit)
            #break;
        #elif pnext[0] in self.arnet.receptorproms:
         #   return self.buildcircuit(queue, circuit, force_inputs=True)
            #break;
        else:
            blacklist.append(pnext[0])
            index = self.arnet.promlist.index(pnext[0])
            signature = self.arnet.proteins[index][4]
            fun = self.problem.nodemap_(signature,self.problem.funs)
            #print "NUUUUTS:", self.printgraph()
            inputs = self.getinputs(index)
            #print self.arnet.promlist
            #print 'index is',index
            #print inputs
            inputs = filter(lambda x: x[0] not in blacklist, inputs)
            ar = self.problem.arity[fun]
            if ar > 0:
                    inputs = inputs[:ar]
            if not inputs:
                #print signature, math.tanh(signature.int)
                #fun = str(math.tanh(signature.int))
                fun = str(random.choice([0.0,1.0]))
                #fun = self.problem.nodemap_(signature,self.problem.terms)
            circuit.append((pnext[0],fun,[i[0] for i in inputs]))

            queue.extend(inputs)
            queue.sort(key=lambda x: x[1], reverse=True)
            return self.buildcircuit(queue,circuit)

    def getinputs(self, index):
        inputs=[]
        #rlist = self.arnet.receptorproms# [r[0] for r in self.arnet.receptors]
        #pmap = zip(self.arnet.promlist+rlist,self.graph[:,index])
        for e,i in zip(self.graph[:,index],
                       range(self.graph.shape[0])):
                       #range(len(self.arnet.promlist))):
            if e > 0:
                inputs.append(((self.arnet.promlist + self.arnet.receptorproms)[i],e))
        inputs.sort(key=lambda x: x[1], reverse = True)
        return inputs

    def getinputnodes(self):
        return map(lambda x: (x[0],
                              self.problem.terms[x[0]],
                              []),
                   self.arnet.receptors)

    def printgraph(self):
        s = 'digraph best {\nordering = out;\n'
        shape = 'hexagon'
        arnet = self.arnet
        labelidx = 0
        #print "CODE"
        #print self.graph
        for outprom in arnet.effectorproms:
            s += '%i [label="%s",shape="hexagon"];\n' % (outprom, labelidx)
            for e,i in zip(self.graph[:,arnet.numtf+arnet.numrec+labelidx],
                           range(self.graph.shape[0])):
                if e > 0:
                    s += '%i -> %i [dir=back,weight=%i];\n' % \
                             (outprom, (arnet.promlist + arnet.receptorproms)[i],e)
            labelidx += 1

        for tf,idx in zip(arnet.promlist,range(arnet.numtf)):
            s += '%i [label="%s"];\n' % (tf, labelidx)
            for e,i in zip(self.graph[:,idx],
                           range(self.graph.shape[0])):
                if e > 0:
                    s += '%i -> %i [dir=back,weight=%i];\n' % \
                             (tf, (arnet.promlist + arnet.receptorproms)[i],e)

            labelidx += 1

        for rec in arnet.receptorproms:
            s += '%i [label="%s",shape="rectangle"];\n' % (rec, labelidx)
            labelidx += 1
        s += '}'
        return s



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
