import ConfigParser
import os
import logging
import logging.config
from operator import itemgetter, attrgetter
import random
from functools import partial
from time import clock, sleep
import copy
from bitstring import BitStream
from utils import *
from utils.bitstrutils import *
from numpy import cumsum
from abc import *

log = logging.getLogger(__name__)

evologhead = str('Best\tNumProteins\tNumFunctions\tNumInactive\t'+
		 'GeneSize\tEvaluations\tAvgBest\tAvgNumProteins\t'+
		 'AvgNumFunctions\tAvgGeneSize')

class Agent:
	__metaclass__ = ABCMeta
	@abstractproperty
	def genotype(self): pass
	@abstractproperty
	def phenotype(self): pass
	@abstractproperty
	def fitness(self): pass
	

#Problem base for htis module
class Problem:
	__metaclass__ = ABCMeta
	@abstractproperty
	def funs(self): pass
	@abstractproperty
	def terms(self): pass
	@abstractproperty
	def labels(self): pass
	@abstractproperty
	def arity(self): pass

	def __init__(self, evaluate, nodemap, print_):
		self.nodemap_ = nodemap
		self.eval_ = evaluate	
		self.print_ = partial(print_, labels = self.labels)
		self.arity = dict(zip(self.funs,[0]*len(self.funs)))
		
### Main class of the library
class EvoDevoWorkbench:
	evolog = logging.getLogger('evolution')
	circuitlog = logging.getLogger('circuit')
	arnlog = logging.getLogger('arntofile')
	def __init__(self, configfname, problem, codefun, agentclass):
		config = ConfigParser.ConfigParser()
		config.readfp(open(configfname))
		
		logging.config.fileConfig(config.get('default','logconf'))
    		log.info('Starting evolution...')
		self.evolog.critical(evologhead)

		self.problem = problem
		self.codefun = codefun
		self.popsize = config.getint('default','popsize')
		self.parentpsize = config.getint('default','parentpopsize')
		self.maxiters = config.getint('default','numiters')
		self.popratio = self.popsize / self.parentpsize
		self.numevals = 0
		
		opnames = config.get('default','operators')
		oprates = config.get('default','oprates')
		self.opargs = config.get('default','opargs').split(',')
		self.ops_, self.oprates = _initialize_ops(opnames,oprates)
		log.debug(self.ops_)
		log.debug(self.oprates)
		arncfg = config.get('default','arnconf')
		self.arnconfig = ConfigParser.ConfigParser()
		self.arnconfig.readfp(open(arncfg))
		self.agentclass = partial(agentclass, config = self.arnconfig)
		log.info('Initializing population...')
		self.population = [self.agentclass() 
				   for p in range(self.popsize)]
		self.best = self.agentclass()
		self.parents = []
		mutrate = config.getfloat('default','mutrate')
		self.mutate_ = partial(bitflipmutation,
				       mutrate = mutrate)
		self.itercount = 0
		
	def step(self):
		log.info('Mapping population to circuit...')
		for i in self.population:
			i.phenotype = self.codefun(i.genotype,self.problem)
		
		log.info('Evaluating population...')
		for i in self.population:
			i.fitness = self.problem.eval_(i.phenotype)
		self.numevals += len(self.population)
		
		log.info('Selecting parents...')
		#there is a bug in functools.partial in Py2.6
		#cannot use deepcopy with partial'ed objects
		self.parents.extend(self.population[:])
		self.parents.sort(key=lambda a: a.fitness)
		del self.parents[self.parentpsize:]
	
		log.info('Generating next population...')
		#This will have to be refactored to use crossover!
		op_ = _selectop(self.ops_,self.oprates)
		gcodeclass = self.parents[0].genotype.code.__class__
		self.population = [self.agentclass( 
				gcode = self.mutate_(
					op_(gcodeclass(p.genotype.code),
					    *self.opargs))) 
				   for p in self.parents
				   for i in range(self.popratio)]
		
	def run(self, terminate = (lambda x,y: x <= 1e-3 or y <= 0)):
		start = clock()
		while(not(terminate(self.best.fitness,
				    self.maxiters-self.itercount))):
			log.info('--- Iteration #%i ---', self.itercount)
			self.step()
			log.debug(self.parents[0].fitness)
			log.debug(self.best.fitness)
			if self.parents[0].fitness < self.best.fitness:
				#self.best = copy.deepcopy(self.parents[0])
				self.best = self.parents[0]
				self.circuitlog.critical(
					self.problem.print_(self.best.phenotype))
				log.info('Best:\n%s',
					 self.problem.print_(self.best.phenotype))
					 
			self._logiteration()
			self.itercount += 1
		print str(clock()-start)+ "sec."
		return self.best

	def _logiteration(self):
		log.info('Logging...')
		tolog = [self.best.fitness]
		tolog.append(len(self.best.genotype.promlist))
		tolog.append(len(self.best.phenotype))
		tolog.append(tolog[-2] - tolog[-1])
		tolog.append(len(self.best.genotype.code))
		tolog.append(self.numevals)
		tolog.append(sum([p.fitness 
				  for p in self.parents])/self.parentpsize)
		tolog.append(sum([len(p.genotype.promlist) 
				  for p in self.parents])/self.parentpsize)
		tolog.append(sum([len(p.phenotype) 
				  for p in self.parents])/self.parentpsize)
		tolog.append(sum([len(p.genotype.code) 
				  for p in self.parents])/self.parentpsize)
		self.evolog.critical(
			reduce(lambda m,n: str(m)+'\t'+str(n), tolog))

		
def _idle(gcode, *args):
	return gcode

def _initialize_ops(opnames, rates):
	if not opnames or not rates:
		ops, ratescum = [_idle],[1.0]
	else:
		opnames = opnames.split(',')
		mainmod = __import__('__main__')
		try:
			ops = [getattr(mainmod, x) for x in opnames]
		except:
			log.warning('Operator not found, loading aborted...')
			log.debug('__main__ module contents: %s', dir(mainmod))
			return [_idle],[1.0]
		ops.append(_idle)
		rates = rates.split(',')
		rates = [float(r) for r in rates]
		rates.append(1.0 - sum(rates))
		ratescum = cumsum(rates)
	return ops,ratescum
		

def _selectop(ops, rates):
	rnd = random.random()
	for i in range(len(rates)):
		if rnd < rates[i]:
			return ops[i]
###########################################################################
### Test                                                                ###
###########################################################################

if __name__ == '__main__':
	from rencode import *
	p = Problem(defaultnodemap,regressionfun, lambda x,y: 1)
	edw = EvoDevoWorkbench('configfiles/test.cfg',p,buildcircuit,Agent)
	edw.run()

