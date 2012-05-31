import logging
from functools import partial
from bitstring import BitStream
import math
import arn
from evodevo import Problem, Agent
from utils import *
from utils.bitstrutils import *

log = logging.getLogger(__name__)

### Problem base to use with ReNCoDe
class ReNCoDeProb(Problem):
	#read fun set from config file??
	funs = [ 'add_', 'sub_', 'mul_', 'div_']
	terms = [ 'inputs[0]' ]
	labels = {'add_':'+', 'sub_':'-', 'mul_':'*', 'div_':'/',
		  'inputs[0]':'x'}
	arity = {}#{'add_':0, 'sub_':0, 'mul_':0, 'div_':0}
	feedback = False
	def __init__(self, evaluate):
		Problem.__init__(self, evaluate, defaultnodemap, 
				 printdotcircuit)

### Agent model to use with this CoDe module
class ReNCoDeAgent(Agent):
	genotype = None
	phenotype = None
	fitness = None
	def __init__(self, config, gcode = None):
		Agent.__init__(self)
		generator = arn.bindparams(config, arn.generatechromo)
		if gcode == None:
			gcode = generator()
			
		self.genotype = arn.ARNetwork(gcode,config)
		self.phenotype = []
		self.fitness = 1e4
	
	def __str__(self):
		return "### Agent ###\n%s\n%s: %f" % (self.arn,self.circuit,
						      self.fitness)

def regressionfun(mapped, node_inputs, inputs ):
	if not node_inputs:
		return eval(mapped)
	mainmod = __import__('__main__')
	if len(node_inputs) == 1:
		return getattr(mainmod, mapped)(node_inputs[0])
	return reduce(lambda m,n: getattr(mainmod, mapped)(m,n),
		      node_inputs)

def mergefun(mapped, node_inputs, inputs):
	if not node_inputs:
		return mapped
	result = mapped + '(' 
	if node_inputs:
		result += reduce(lambda m, n: m + ',' + n, node_inputs)
	result += ')'
	return result
	
def defaultnodemap(signature, mappingset):
	if len(mappingset) < 2:
		return mappingset[0]
	index = BitStream(bin=applymajority(
			signature,
			int(math.ceil(math.log(len(mappingset),2)))))
	intindex = index.uint
	if intindex >= len(mappingset):
		intindex -= len(mappingset)
	return mappingset[intindex]
	
def evaluatecircuit(circuit, circuitmap, resultdict, *inputs):
	for i in range(len(circuit)):
		inputvalues =  [resultdict[abs(v)] for v in circuit[-1-i][2]]
		result = circuitmap(circuit[-1-i][1], 
				    inputvalues,
				    inputs)
		resultdict[circuit[-1-i][0]] = result
	return result

def buildcircuit(arn, problem, **kwargs):
	"""Returns circuit to be fed into evaluation function"""
	if not arn.promlist:
		return []
	graph = arn.ebindings - arn.ibindings
	_cleanpairs(graph)
	promlist = [(p[0], 
		     _getinputlist(
				arn.promlist, 
				graph[:,arn.promlist.index(p[0])].tolist()))
		    for p in arn.proteins]
	
	promlist.sort(key=lambda x: len(x[1]),reverse=True)
	circuit = []
	pdict = dict(zip(arn.promlist,arn.proteins))
	_buildcircuit(circuit, problem, pdict, dict(promlist), 
		      graph, [promlist[0][0]], [], [])
	return circuit

def _buildcircuit(circuit, problem, proteindict, inputdict, graph,
		  pqueue, secondqueue, blacklist):
	if not pqueue:
		if not secondqueue: return circuit
		else: 
			secondqueue.sort(key=lambda x: len(inputdict[x]),
					 reverse=True)
			pqueue.append(secondqueue.pop(0))
			
	next = pqueue.pop(0)
	pnext = proteindict[next]
	blacklist.append(next)
	inputs = []
	for i in inputdict[next]:
		if i in blacklist:
			if problem.feedback:
				inputs.append(-i)
		else:
			inputs.append(i)
	if inputs:
		fun = problem.nodemap_(pnext[4],problem.funs)
		ar = problem.arity[fun]
		if ar > 0:
			inputs = inputs[:ar]
	else:
		fun = problem.nodemap_(pnext[4],problem.terms)

	circuit.append((next, fun , inputs))

	secondqueue.extend([inp for inp in circuit[-1][2]
			    if  inp not in (pqueue+secondqueue+blacklist) 
			    and (inp >= 0)]) 
	
	pqueue, secondqueue = _mergequeues(pqueue, secondqueue, inputdict)
	pqueue.sort(key=lambda x: len(inputdict[x]), reverse = True)

	return _buildcircuit(circuit, problem, proteindict, inputdict, 
			     graph, pqueue,secondqueue, blacklist)	

def _cleanpairs(matrix):
	s = len(matrix)
	for i in range(s):
		for j in range(i):
			if matrix[i][j] >= matrix[j][i]:
				matrix[j][i] = 0
			else:
				matrix[i][j] = 0

def _mergequeues(q1, q2,inputdict):
	disjunction = q2[:]
	dependent = False
	for e in q2:
		if e in [p 
			 for pq_el in q1
			 for p in inputdict[pq_el]]:
			dependent = True
		if e in [p 
			 for pq_el in q2
			 for p in inputdict[pq_el]
			 if pq_el != e]:
			dependent = True
		
		if dependent:
			disjunction.remove(e)
	q1.extend(disjunction)
	for e in disjunction:
		q2.remove(e)
	return q1,q2

def _getinputlist(promlist, weights):
	"""Returns the relevant inputs given the weights."""
	inputs=[]
	pmap = zip(promlist,weights)
	for p,w in pmap:
		if w > 0:
			inputs.append(p)
	return inputs

def printcircuit(circuit):
	s = ''
	for c in circuit:
		s += "%i [%s]: %s\n" % (c[0], c[1],
				      reduce(lambda m,n: "%s %s " % (m,n),c[2],
					     ""))
	return s[:-2]

def printdotcircuit(circuit, labels=None):
	s = 'digraph best {\nordering = out;\n'
	for c in circuit:
		s += '%i [label="%s"];\n' % (c[0], c[1] if not labels 
					     else labels[c[1]])
		for inp in c[2]:
			aux = "dir=back"
			if inp < 0:
				aux += ",style=dotted"
			s += '%i -> %i [%s];\n' % (c[0],abs(inp),aux)
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
	cfg.readfp(open('test.cfg'))
	arncfg = ConfigParser.ConfigParser()
	arncfg.readfp(open('arnsim.cfg'))
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
