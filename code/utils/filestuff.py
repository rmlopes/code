from logging.handlers import RotatingFileHandler
import logging
import cPickle as pickle

log = logging.getLogger(__name__)

class WriteRotateFileHandler(RotatingFileHandler):
        def __init__(self, filename):
                RotatingFileHandler.__init__(self, filename, 'w',0,0)

        def emit(self, record):
            self.doRollover()
            RotatingFileHandler.emit(self, record)

class DefaultRunLog:
    evologhead = str('Best\tNumProteins\tNumFunctions\tNumInactive\t'+
                     'GeneSize\tEvaluations\tAvgBest\tAvgNumProteins\t'+
                     'AvgNumFunctions\tAvgGeneSize\tMutRate')
    def __init__(self, mainlogger):
        self.main = mainlogger
        self.evolog = logging.getLogger('evolution')
        self.circuitlog = logging.getLogger('circuit')
        self.arnlog = logging.getLogger('arntofile')
        self.validatelog = logging.getLogger('validationfile')
        self.evolog.critical(self.evologhead)

    def step(self, parents, **kwargs):
        log.info('Logging...')
        best = parents[0]
        parentpsize = float(len(parents))
        tolog = [best.fitness]
        tolog.append(len(best.genotype.promlist))
        tolog.append(len(best.phenotype))
        tolog.append(tolog[-2] - tolog[-1])
        tolog.append(len(best.genotype.code))
        tolog.append(kwargs['numevals'])
        tolog.append(sum([p.fitness
                          for p in parents])/parentpsize)
        tolog.append(sum([len(p.genotype.promlist)
                          for p in parents])/parentpsize)
        tolog.append(sum([len(p.phenotype)
                          for p in parents])/parentpsize)
        tolog.append(sum([len(p.genotype.code)
                          for p in parents])/parentpsize)
        #tolog.append(self.mutrate)
        self.evolog.critical(
                reduce(lambda m,n: str(m)+'\t'+str(n), tolog))

    def dumpcircuit(self, best):
        self.circuitlog.critical(pickle.dumps(best.pickled()))

    def dumpancestortree(self, best):
        ancestorlog = logging.getLogger('ancestortrace')
        ancestorlog.critical(pickle.dumps(best.pickled()))
        self._print_ancestors(best, ancestorlog)
        print best.mutlog, best.oplog
        mutshares = self._sum_mutations(best)
        mutshares = filter(lambda s: s>-1, mutshares)
        print mutshares
        if mutshares:
            print "AVG_NEUTRAL_MUT:", sum(mutshares) / float(len(mutshares))
        shares = self._sum_neutralshare(best)
        print shares
        print "AVG_NEUTRAL_GENOME:", sum(shares) / float(len(shares))
        opshares = self._sum_ops(best)
        if opshares:
            for op, sh in opshares.items():
                print "AVG_NEUTRAL_"+op+":"+str(sh[1]/float(sh[0]))

    def _print_ancestors(self, ind, alog):
        parent = ind.parent
        if not parent:
            return
        else:
            alog.critical(pickle.dumps(parent.pickled()))
            self._print_ancestors(parent, alog)

    def _sum_mutations(self, ind):
        parent = ind.parent
        if not parent:
            if ind.mutlog[0] != 0:
                return [ind.mutlog[1]/float(ind.mutlog[0])]
            else:
                return [-1.0]
        else:
            if ind.mutlog[0] != 0:
                l = [ind.mutlog[1]/float(ind.mutlog[0])]
            else:
                l = [-1.0]
            l.extend(self._sum_mutations(parent))
            return l

    def _sum_ops(self, ind):
        parent = ind.parent
        if not parent:
            return ind.oplog
        else:
            p_oplog = self._sum_ops(parent)
            for k,v in p_oplog.items():
                try:
                    ind.oplog[k][0] += v[0]
                    ind.oplog[k][1] += v[1]
                except KeyError:
                    ind.oplog[k] = v
            return ind.oplog

    def _sum_neutralshare(self,ind):
        parent = ind.parent
        if not parent:
            return [ind.genotype.getneutralshare()]
        else:
            #neutralsum = self._sum_neutralshare(parent)
            l =  [ind.genotype.getneutralshare()]
            l.extend(self._sum_neutralshare(parent))
            return l


#if __name__ == '__main__':
       # wrfh = WriteRotateFileHandler('../../results/code_circuit_0.save')
