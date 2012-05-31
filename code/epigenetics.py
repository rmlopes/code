import logging
import functools
from arn import bindparams, getbindings
from grammar import *
from protectedmath import protdiv
from operator import add, sub, mul
from bitstring import BitStream
from bitstringutils import applymajority
import math
from numpy import subtract
from gp_ant import *

log = logging.getLogger(__name__)



def cleanpairs(matrix):
	s = len(matrix)
	for i in range(s):
		for j in range(i):
			if matrix[i][j] >= matrix[j][i]:
				matrix[j][i] = 0
			else:
				matrix[i][j] = 0

def _buildcircuitfunctions(individual, **kwargs):
	"""
	Returns circuit to be fed into evaluation function
	"""
   	sz=32
	genome, protdict, epigenome,lf,inactivenum = individual
	#print individual
	epigenome = dict(epigenome)
	#print 'EPI: ', epigenome
        #remove proteins from deactivated genes
	epithr = (1.0/float(len(protdict)))*kwargs['epithreshold']
	oldnum = len(protdict)
	#print epigenome
	#print proteins
	proteins = filter(lambda x: x[0] in epigenome.keys() and 
			  epigenome[x[0]] >= epithr, 
			  protdict.values())
	print 'Removed %i deactivated/new proteins' % (oldnum - len(proteins))
        #excite weights
        eweights = getbindings(0, proteins, kwargs['match_threshold'], 
                             kwargs['bindingsize'])
	iweights = getbindings(1, proteins, kwargs['match_threshold'], 
			      kwargs['bindingsize'])	
	eweights -= iweights
	
	cleanpairs(eweights)
	proteinorder = [p[0] for p in proteins]
	promlist = [(p[0],1, filter(lambda x: x[0] != p[0],
				    _getinputlist(
					proteinorder, 
					eweights[:][proteinorder.index(p[0])].tolist()))) 
		    for p in proteins]
	#print 'PROTEINS:'
	#for p in proteins: print str(p[0]) + ' ' + str(p[5])
	#print 'PROMLIST: ',promlist
	#promlist.sort(key=lambda x: len(x[2]), reverse=True)
	log.debug(promlist)
	#print "ORDERED PLIST:"
	#for p in promlist: print p

	circuitset = []
	#for p in filter(lambda x: x[5] >= 1.0/len(proteins), proteins):
	#	circuit=[]
	#	proteindict = dict(zip(promlist,proteins))
	#	circuitset.append(_buildcircuit( circuit, [(p[0],1)], proteindict, promlist, eweights))
	circuit=[]

	proteindict = dict(zip(proteinorder,proteins))
	circuitset.append(_buildcircuit( circuit, [random.choice(promlist)], proteindict, proteinorder, eweights,[],[],kwargs['feedback']))
	#i=0
	#for c in circuitset:
	#	print "CIRCUIT ",i,':'
	#	_printcircuit(c)
		#for n in c: print n
	#	i+=1
	
	log.debug(circuit)
	return circuit

def _buildcircuit(circuit, pqueue, proteindict, promlist, 
		  weights, blacklist,sec_queue, feedback):
	if not pqueue:
		if not sec_queue: return circuit
		else: 
			sec_queue.sort(key=lambda x: len(x[2]),reverse=True)
			pqueue.append(sec_queue.pop(0))
			
	next = pqueue.pop(0)#(id,w,inputs)
	pnext = proteindict[next[0]]
	blacklist.append(next[0])
	#print 'BLACKLIST:', blacklist
	#queueblacklist = [p[0] for p in pqueue]
	#print 'p[0] after = %s', first
	#adiciona no ao circuito (32-bit sign.) com lista de inputs
	#print 'Wcol: ',weights[:][promlist.index(next[0])].tolist()
	inputs = []
	for i in next[2]:
		if i[0] in blacklist:
			if feedback:
				inputs.append((-i[0],i[1]))
		else:
			inputs.append(i)
	circuit.append(([next[:1],pnext[4]], inputs))

	queueblacklist = [q[0] for q in pqueue]
	queueblacklist.extend([q[0] for q in sec_queue])
	sec_queue.extend([nodetuple+
		       (nofilterfeedbackinputs(nodetuple,promlist,
					       weights, blacklist),)
		       for nodetuple in circuit[-1][1]
		       if  nodetuple[0] not in queueblacklist+blacklist and 
			  nodetuple[0] > 0 ]) 

	pqueue, sec_queue = mergequeues(pqueue, sec_queue, blacklist)
	#no sort

	return _buildcircuit(circuit, pqueue, proteindict, promlist, 
			     weights,blacklist,sec_queue, feedback)

def mergequeues(q1, q2, blacklist):
	disjunction = q2[:]
	dependent = False
	for e in q2:
		if e[0] in [p[0] 
		       for pq_el in q1
		       for p in pq_el[2]]:
			dependent = True
		if e[0] in [p[0] 
		       for pq_el in q2
		       for p in pq_el[2]]:
			dependent = True
		
		if dependent:
			disjunction.remove(e)
	q1.extend(disjunction)
	for e in disjunction:
		q2.remove(e)
	return q1,q2

def nofilterfeedbackinputs(nodetuple, promlist, weights, blacklist):
	ilist = filter( lambda x: x[0] != nodetuple[0],
			_getinputlist(promlist, 
				      weights[:][promlist.index(nodetuple[0])].tolist()))
	result = []
	return ilist

#TODO: should be trimmed to return only the promotor (id)
def _getinputlist(promlist, weights):
	"""
	Returns the relevant inputs given the weights.
	"""
	#print 'W: ', weights
	inputs=[]
	pmap = zip(promlist,weights)
	for p,w in pmap:
		if w > 0:
			inputs.append((p,w))
	return inputs

	
def getmapping(arncfg, fback = False):
	"""
	Returns the mapping (from genome and expression to circuit)
	"""
	return functools.partial( bindparams(arncfg,_buildcircuitfunctions),
				  feedback = fback) 

def evaluateregression(circuit, target, inputs):
	#print 'Evaluating...'
	#_printcircuit(circuit)
	errors = [abs(target(inp) - 
		      _evaluatecircuit(circuit, _defaultmap,_regressionfun,inp)) 
		  for inp in inputs]
	sum_ = sum(errors)
	return 1e4 if math.isinf(sum_) else sum_

def evaluate_ant(circuit):
	dmap = functools.partial(_defaultmap, funs=ANT_FSET, terms=ANT_TSET)
	#_printcircuit(circuit, funs=ANT_FSET, terms=ANT_TSET)
	routine = _evaluatecircuit(circuit, dmap, _antfun)
	#print routine
	ant.runstring(routine,True)
	return 89 - ant.eaten
	
			
#FIXME
def _evaluatecircuit(circuit, mapfun, domainfun, *inputs):
	resultdict = dict()
	#print circuit
	for i in range(len(circuit)):
		inputvalues = _getinputvalues(circuit[-1-i][1],resultdict) 
		result = domainfun(mapfun,
				   circuit[-1-i][0][1], 
				   inputvalues,
				   inputs)
		#result = _genericfunction(circuit[-1-i][0][1], 
		#			  inputvalues,
		#			  mapfun,
		#			  domainfun,
		#			  inputs)
		#print inputvalues
		#if not circuit[-1-i][1] and not noeval:
			#result = eval(result)
		#circuit[-1-i] += (result,)
		resultdict[circuit[-1-i][0][0][0]] = result
		#print resultdict
	return result

def _getinputvalues( inputlist, results):
	return [results[v[0]] for v in inputlist]
	
	

FUNCTION_SET = [ 'add', 'sub', 'mul', 'protdiv','add','sub','mul','protdiv']
TERMINAL_SET = [ 'inputs[0]','inputs[0]','inputs[0]','inputs[0]',
		 'inputs[0]','inputs[0]','inputs[0]','inputs[0]']
	
def _regressionfun(mapfun, signature, node_inputs, inputs ):
	if not node_inputs:
		str_ = mapfun(signature,node_inputs)
		#print str_
		return eval(str_)
	return reduce(lambda m,n: globals()[mapfun(signature,node_inputs)](m,n),
		      node_inputs)

def _antfun(mapfun, signature, node_inputs, inputs):
	if not node_inputs:
		return mapfun(signature,node_inputs)
	result = mapfun(signature,node_inputs) + '(' 
	if node_inputs:
		result += reduce(lambda m,n: m+','+n, node_inputs)
	result += ')'
	#print result
	return result
	
def _defaultmap(signature,inputs, funs = FUNCTION_SET, terms = TERMINAL_SET):
	index = BitStream(bin=applymajority(
			signature,
			int(math.log(len(funs),2))))
	#print index.bin
	if not inputs:
		#print 'term = ', TERMINAL_SET[signature.uint%len(TERMINAL_SET)]
		#return TERMINAL_SET[signature.uint%len(TERMINAL_SET)]
		return terms[index.uint]
	#fun = FUNCTION_SET[signature.uint%len(FUNCTION_SET)]
	return funs[index.uint]
	

def _printcircuit(circuit, funs=FUNCTION_SET, terms=TERMINAL_SET):
	print 'Circuit: '
	for c in circuit:
		print "%i [%s]: %s" % (c[0][0][0],
				       _defaultmap(c[0][1],
						   c[1],
						   funs,
						   terms),
				       reduce(lambda m,n: "%s %s " % (m,n[0]) ,
					      c[1],
					      ""))

def _printdotcircuit(circuit, funs=FUNCTION_SET, terms=TERMINAL_SET, labels=[]):
	s = 'digraph best {\nordering = out;\n'
	for c in circuit:
		nodemap = _defaultmap(c[0][1], c[1], funs, terms)
		s += '%i [label="%s"];\n' % (c[0][0][0],
					     nodemap 
					     if not labels else labels[nodemap])
		for inp in c[1]:
			s += '%i -> %i [dir=back];\n' % (c[0][0][0],inp[0])
	s += '}'
	return s

ANTLABELS = {'ant.move_forward':'M',
	     'ant.if_food_ahead':'IFA', 
	     'progN':'progn',
	     'ant.turn_left':'L',
	     'ant.turn_right':'R'}

ANT_FSET = ['ant.if_food_ahead', 'progN', 'ant.if_food_ahead','progN',
	    'ant.if_food_ahead', 'progN', 'ant.if_food_ahead','progN']
ANT_TSET = ['ant.move_forward','ant.turn_left','ant.turn_right','ant.move_forward',
	    'ant.move_forward','ant.turn_left','ant.turn_right','ant.move_forward']
#def _antfunction(signature, inputs):
#	index = BitStream(bin=applymajority(
#			signature,
#			int(math.log(len(FUNCTION_SET),2))))
#	if not inputs:
#		return ANT_TSET[index.uint]
#	fun = ANT_FSET[index.uint]
#	return fun + '(' + reduce(lambda m,n: m+','+n, inputs) + ')'
				    
	
