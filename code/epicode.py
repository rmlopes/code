import logging
from functools import partial
from bitstring import BitStream
from utils import *
from utils.bitstrutils import *
from rencode import ReNCoDeAgent, recursivebuild, cleanpairs

log = logging.getLogger(__name__)

### Agent model to use with this CoDe module
class EpiCoDeAgent(ReNCoDeAgent):
        epigenome = None
	def __init__(self, config, gcode = None):
		ReNCoDeAgent.__init__(self,config,gcode)
		self.epigenome = BitStream(bin='1'*len(self.genotype.promlist))

	def __str__(self):
		return "### Agent ###\n%s\n%s\n%s\nFitness: %f" % (
			self.arn,self.circuit,self.epigenome,
			self.fitness)


class LocalSearch:
	def __init__(self, numiters, codefun):
		self.numiters = numiters
		self.codefun = codefun
		
	def search(self, agent,problem):
		for i in range(self.numiters):
			idx = random.randint(0,len(agent.epigenome)-1)
			agent.epigenome.invert(idx)
			newphenotype = self.codefun(agent,problem)
			if newphenotype != agent.phenotype:
				newfitness = problem.eval_(newphenotype)
				if newfitness < agent.fitness:
					#log.debug('Successful local search')
					agent.fitness = newfitness
					agent.phenotype = newphenotype
	

def buildcircuit(agent, problem, **kwargs):
	"""Returns circuit to be fed into evaluation function"""
	circuit = []
	arn = agent.genotype
	epig = agent.epigenome.bin
	if not arn.promlist:
		return []
	epipromlist = filter(lambda x: int(x[1]), zip(arn.promlist,epig))
	if epipromlist:
		newproms,garbage = zip(*epipromlist)

		graph = arn.ebindings - arn.ibindings
		cleanpairs(graph)
		promlist = [(p, 
			     _getinputlist(
					newproms, 
					graph[:,arn.promlist.index(p)].tolist(),
					epig))
			    for p in newproms]

		promlist.sort(key=lambda x: len(x[1]),reverse=True)
		pdict = dict(zip(newproms,
				 filter(lambda x: x[0] in newproms, arn.proteins)))
		recursivebuild(circuit, problem, pdict, dict(promlist), 
			      graph, [promlist[0][0]], [], [])
	return circuit


def _getinputlist(promlist, weights, epigenome):
	"""Returns the relevant inputs given the weights."""
	inputs=[]
	epigweights = zip(weights,epigenome)
	pmap = zip(promlist,[w[0] for w in epigweights if w[1]])
	for p,w in pmap:
		if w > 0:
			inputs.append(p)
	return inputs
