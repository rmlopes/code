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
from utils.filestuff import DefaultRunLog
from numpy import cumsum
from abc import *
import cPickle as pickle
import sys

log = logging.getLogger(__name__)

class Agent:
        __metaclass__ = ABCMeta
        parentfit = 1e4
        @abstractproperty
        def genotype(self): pass
        @abstractproperty
        def phenotype(self): pass
        @abstractproperty
        def fitness(self): pass

        def __init__(self, parent = None):
                self.parent = parent
                #for backward compatibility
                if parent:
                        self.parentfit = parent.fitness
                #name, (total, neutral)
                self.oplog = dict()
                #total, neutral
                self.mutlog = [0,0]

        def traceops(self, opname, neutral):
            log.debug("tracing use of operator %s (neutral=%s)"%(opname,
                                                                 neutral))
            try:
                self.oplog[opname][0] += 1
                if neutral:
                    self.oplog[opname][1] += 1
            except KeyError:
                self.oplog[opname] = [1, int(neutral)]

        def tracemuts(self, total, neutral):
            log.debug("tracing %i total mutations (%i neutral)"%(total,
                                                                  neutral))
            self.mutlog[0] += total
            self.mutlog[1] += neutral

        def pickled(self): return self

        def print_(self): return str(self)

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

        def __init__(self, evaluate):
                #self.nodemap_ = nodemap
                self.eval_ = evaluate
                #self.print_ = partial(print_, labels = self.labels)
                self.arity = dict(zip(self.funs,[0]*len(self.funs)))


def register_adf(adfcount,adf,problem):
        name = 'adf%i'%(adfcount,)
        mainmod = __import__('__main__')
        setattr(mainmod, name , adf)
        problem.terms = problem.terms +[name]
        log.debug('Registered %s', name)

### Main class of the library
class EvoDevoWorkbench:
        def __init__(self, loadedconfig, problem, logmanager = DefaultRunLog,
                     **kwargs):
                config, self.arnconfig = loadedconfig
                logging.config.fileConfig(config.get('default','logconf'),
                                          disable_existing_loggers=False)
                log.info('Setting up evolutionary workbench...')
                self.runlog = logmanager(log)
                self.trace =  config.getboolean('default',
                                                'ancestortrace')
                mainmod = __import__('__main__')

                self.problem = problem
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

                self.xover_ = config.get('default','xover')
                if self.xover_:
                    log.info('Initializing crossover operator...')
                    self.xover_ = getattr(mainmod,self.xover_)
                    self.xrate = config.getfloat('default','xrate')

                #arncfg = config.get('default','arnconf')
                #self.arnconfig = ConfigParser.ConfigParser()
                #self.arnconfig.readfp(open(arncfg))
                aclass = config.get('default','agent').split('.')
                mod = aclass[0] + "." + aclass[1]
                self.device = __import__(mod,fromlist=aclass[2])
                self.problem.eval_ = partial(self.problem.eval_,
                                             device = self.device)
                log.info("CoDe module: %s"%(mod,))
                log.info("Agent class: %s"%(aclass[2],))
                self.agentclass = partial(getattr(self.device, aclass[2]),
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
                        i.fitness = self.problem.eval_(i.phenotype)
                else:
                    self.gui.pop = self.population
                    self.gui.unpause()
                    self.parents = []
                    for i in self.gui.selected: #[-self.parentpsize]:
                        self.population[i].fitness -= 1
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

                #TODO: clean this code
                offsprings = []
                if self.xover_:
                    while len(offsprings)/float(self.popsize)\
                    < self.xrate:
                        breeders = random.sample(self.parents,2)
                        offsprings.extend(self._reproduct_agents(*breeders))
                    log.debug('Created %i offsprings by xover.'
                              %(len(offsprings),))

                if not self.interactive:
                    mutants = [self._create_mutant(p)
                               for p in self.parents
                               for i in range(self.popratio)]
                    self.population = (offsprings +
                                       mutants[:(self.popsize-len(offsprings))])
                else:
                    if len(offsprings) < self.popsize:
                        mutants = [self._create_mutant(p)
                                   for p in self.parents
                                   for i in range(self.popratio-1)]
                    else:
                        mutants = []
                    self.population = (self.parents + offsprings +
                                       mutants[:self.popsize -
                                            (self.parentpsize+len(offsprings))])

        def _reproduct_agents(self, p1, p2):
            return [self.agentclass(gcode = self.mutate_(o),
                                    parent = p1,
                                    problem = self.problem)
                    for o in self.xover_(p1,
                                         p2)]

        def _create_mutant(self, parent):
            op_ = _selectop(self.ops_,self.oprates)
            gcodeclass = parent.genotype.code.__class__
            gcode = op_(gcodeclass(parent.genotype.code),
                        *self.opargs,
                        arnet = parent.genotype)
            if self.trace and op_.__name__ != '_idle':
                agent = self.agentclass(gcode = gcode,
                                        parent = parent,
                                        problem = self.problem)
                if agent.phenotype == parent.phenotype:
                    neutraloperator = True
                else:
                    neutraloperator = False

            mut_gcode = self.mutate_(gcodeclass(gcode))
            mutagent = self.agentclass(gcode = mut_gcode,
                                       parent = parent,
                                       problem = self.problem)
            if self.trace:
                if op_.__name__ == '_idle':
                    agent = parent
                mutmask = mut_gcode ^ agent.genotype.code
                idxmut = filter(lambda l: l[0] == '1',
                                zip(mutmask.bin, range(len(mutmask.bin))))
                totalmut = len(idxmut)
                if agent.phenotype == mutagent.phenotype:
                    neutralmut = totalmut
                else:
                    #print idxmut
                    #print agent.genotype.promlist
                    #FIXME: constants
                    neutralmut = 0
                    try:
                        promlist = agent.genotype.promlist + \
                                   agent.genotype.effectorproms
                    except:
                        promlist = agent.genotype.promlist

                    for i in idxmut:
                        n = True
                        for p in promlist:
                            if  p - 88 <= i[1] < p + 168:
                                n = False
                                break
                        if n:
                            neutralmut += 1
                if op_.__name__ != '_idle':
                    mutagent.traceops(op_.__name__, neutraloperator)
                mutagent.tracemuts(totalmut, neutralmut)

            return mutagent

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
                                #NOTE: cannot pickle on python 2.6
                                        #(issue 1398)
                                #self.circuitlog.critical(
                                 #       pickle.dumps(self.best,2))
                                self.runlog.dumpcircuit( self.best )
                                log.info('Best:\n%s',
                                         self.best.print_())
                                log.debug(self.best.genotype.snapshot())
                                log.debug(self.best.phenotype.printgraph())

                        self.runlog.step(self.parents, numevals=self.numevals)
                        self.itercount += 1
                print str(clock()-start)+ "sec."
                if self.trace:
                    self.runlog.dumpancestortree(self.best)
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





def _idle(gcode, *args, **kwargs):
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
