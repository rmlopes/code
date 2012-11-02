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
                 'AvgNumFunctions\tAvgGeneSize\tMutRate')

class Agent:
        __metaclass__ = ABCMeta
        parentfit = 1e4
        @abstractproperty
        def genotype(self): pass
        @abstractproperty
        def phenotype(self): pass
        @abstractproperty
        def fitness(self): pass

        def __init__(self, pfit):
                self.parentfit = pfit


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



def register_adf(adfcount,adf,problem):
        name = 'adf%i'%(adfcount,)
        mainmod = __import__('__main__')
        setattr(mainmod, name , adf)
        problem.terms = problem.terms +[name]
        log.debug('Registered %s', name)

### Main class of the library
class EvoDevoWorkbench:
        evolog = logging.getLogger('evolution')
        circuitlog = logging.getLogger('circuit')
        arnlog = logging.getLogger('arntofile')

        #TODO: adapt rencode to work with  the
        #buildfunction (codefun) out of here
        def __init__(self, configfname, problem, agentclass):
                config = ConfigParser.ConfigParser()
                config.readfp(open(configfname))

                logging.config.fileConfig(config.get('default','logconf'))
                log.info('Setting up evolutionary workbench...')
                self.evolog.critical(evologhead)
                mainmod = __import__('__main__')

                self.problem = problem
                #self.codefun = codefun
                self.popsize = config.getint('default','popsize')
                self.parentpsize = config.getint('default','parentpopsize')
                self.maxiters = config.getint('default','numiters')
                self.popratio = self.popsize / self.parentpsize

                opnames = config.get('default','operators')
                oprates = config.get('default','oprates')
                self.opargs = config.get('default','opargs').split(',')
                self.ops_, self.oprates = _initialize_ops(opnames,oprates)
                log.debug(self.ops_)
                log.debug(self.oprates)
                arncfg = config.get('default','arnconf')
                self.arnconfig = ConfigParser.ConfigParser()
                self.arnconfig.readfp(open(arncfg))
                self.agentclass = partial(agentclass,
                                          config = self.arnconfig,
                                          problem = problem )
                self.mutrate = config.getfloat('default','mutrate')
                self.orig_mutrate = self.mutrate
                self.mutate_ = partial(bitflipmutation,
                                       mutrate = self.mutrate)
                self.improves = 0
                self.tempevals = 0
                self.adfcount = 0

                self.localsearch = config.get('default','localsearch')
                if self.localsearch:
                        log.info('Initializing local search holder')
                        self.localsearch = getattr(mainmod,
                                                   self.localsearch)(5,codefun)

                self.basicadf = config.get('default','adf')
                if self.basicadf:
                        log.info('Initializing multiplex adf skeleton')
                        self.basicadf = getattr(mainmod,
                                                   self.basicadf)

                self.interactive = config.get('default','gui')
                if self.interactive:
                        log.info('Interactive mode enabled. Initialiazing GUI.')
                        _guiclass = getattr(mainmod,
                                           self.interactive)
                        self.gui = _guiclass(self.popsize,problem=self.problem)

                self.numevals = None
                self.population = None
                self.parents = None
                self.best = None
                self.itercount = None

        def step(self):
                log.info('Evaluating population...')
                if not self.interactive:
                    for i in self.population:
                        i.fitness = self.problem.eval_(i)
                else:
                    self.gui.pop = self.population
                    self.gui.unpause()
                    for i in range(self.popsize):
                        self.population[i].fitness = 1 if i in self.gui.selected else 10
                    if not self.gui._running:
                        self.population[self.gui.selected[0]].fitness = 0



                self.numevals += len(self.population)

                #self.adaptmutrate()

                if self.localsearch:
                        log.info('Performing local search...')
                        for i in self.population:
                                self.numevals += self.localsearch.search(
                                        i,self.problem)

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
                if not self.interactive:
                    self.population = [self.agentclass(
                                gcode = self.mutate_(
                                        op_(gcodeclass(p.genotype.code),
                                            *self.opargs)),
                                parentfit = p.fitness)
                                   for p in self.parents
                                   for i in range(self.popratio)]
                else:
                    self.population = []
                    for p in self.parents:
                        self.population.append(p)
                        for i in range(self.popratio - 1):
                            self.population.append(
                                self.agentclass(
                                    gcode = self.mutate_(
                                        op_(gcodeclass(p.genotype.code),
                                            *self.opargs)),
                                    parentfit = p.fitness))

        def run(self, terminate = (lambda x,y: x <= 1e-3 or y <= 0)):
                mainmod = __import__('__main__')
                start = clock()
                self.numevals = 0
                self.itercount = 0
                log.info('Initializing population...')
                self.population = [self.agentclass()
                                   for p in range(self.popsize)]
                self.best = self.agentclass()
                self.parents = []
                if self.adfcount:
                        self.problem.terms = filter(lambda t: t[:3]!='adf',
                                                    self.problem.terms)
                        self.adfcount = 0

                while(not(terminate(self.best.fitness,
                                    self.maxiters-self.itercount))):
                        log.info('--- Iteration #%i ---', self.itercount)
                        self.step()
                        log.debug(self.parents[0].fitness)
                        log.debug(self.best.fitness)
                        if self.parents[0].fitness < self.best.fitness:
                                #self.best = copy.deepcopy(self.parents[0])
                                self.best = self.parents[0]
                                if self.basicadf:
                                        register_adf(self.adfcount,
                                                     partial(self.basicadf,
                                                             circuit = self.best.phenotype),
                                                     self.problem)
                                        self.adfcount += 1
                                self.circuitlog.critical(
                                        self.problem.print_(self.best.phenotype))
                                log.info('Best:\n%s',
                                         self.problem.print_(self.best.phenotype))

                        self._logiteration()
                        self.itercount += 1
                print str(clock()-start)+ "sec."
                return self.best

        def adaptmutrate(self):
                improves = [p.fitness <= p.parentfit for p in self.population]
                self.improves += improves.count(True)
                self.tempevals += self.popsize


                if self.itercount % 10 == 0:
                        if self.improves >= .2*self.tempevals:
                                self.mutrate /= 0.85
                        else:
                                self.mutrate *= 0.85
                        self.improves = 0
                        self.tempevals = 0

                if self.mutrate < .5*self.orig_mutrate:
                        self.mutrate *= 2

                #print 'mutrate = %f' % (self.mutrate,)

                self.mutate_ = partial(bitflipmutation,
                                       mutrate = self.mutrate)


        def _logiteration(self):
                log.info('Logging...')
                tolog = [self.best.fitness]
                tolog.append(len(self.best.genotype.promlist))
                tolog.append(0)#len(self.best.phenotype))
                tolog.append(0)#tolog[-2] - tolog[-1])
                tolog.append(len(self.best.genotype.code))
                tolog.append(self.numevals)
                tolog.append(sum([p.fitness
                                  for p in self.parents])/self.parentpsize)
                tolog.append(sum([len(p.genotype.promlist)
                                  for p in self.parents])/self.parentpsize)
                tolog.append(0)#sum([len(p.phenotype)
                                 # for p in self.parents])/self.parentpsize)
                tolog.append(sum([len(p.genotype.code)
                                  for p in self.parents])/self.parentpsize)
                tolog.append(self.mutrate)
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
