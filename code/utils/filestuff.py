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
        muts = self._sum_mutations(best)
        if muts[0] != 0:
                print "Avg. Neutral mutations ratio: ", muts[1]/float(muts[0])
        shares = self._sum_neutralshare(best)
        print shares
        print "Avg. Neutral portion of the genome", sum(shares) / len(shares)
        print self._sum_ops(best)

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
            return ind.mutlog
        else:
            pmuts = self._sum_mutations(parent)
            return [ind.mutlog[0]+pmuts[0], ind.mutlog[1] + pmuts[1]]

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
