import sys
from code.utils.config import parsecmd, loadconfig
from code.evodevo import *
from code.operators import *
from code.rencode2 import *
from code.rencode import ReNCoDeProb, evaluatecircuit, regressionfun
from code.utils.mathlogic import *
from code.utils.kfolds import *
from code.extendedarn import bindparams, displayARNresults
import logging

log = logging.getLogger(__name__)


class EnergyProb(ReNCoDeProb):
        #extrafuns = ['exp_','min_','max_','log_','tanh_','sin_','cos_','sinh_','cosh_','tan_']
        funs = ['add_','mul_','div_', 'sub_']#'sqrt_','cos_','reciprocalsum','log_']
        terms = ['inputs[0]']#,'tanh_','sin_','cos_','sinh_','cosh_','tan_']
        arity={}
        labels={}
        def __init__(self, evaluate, **kwargs):
                self.labels = None
                ReNCoDeProb.__init__(self,evaluate,**kwargs)
                self.terms.extend(["inputs[%i]"%i 
                                   for i in range(1,8)])
                #self.funs.extend(self.extrafuns)
                self.arity.update(zip(self.funs,[0]*len(self.funs)))
                #print self.arity
                self.nout = 2
                self.ninp = 8

def evaluate(phenotype, test=False, **kwargs):
    y1error = 0.0
    y2error = 0.0
    mainmod = __import__('__main__')
    src = getattr(mainmod,TESTSET) if test else getattr(mainmod,TRAINSET)


    for r in src:
        inputs = r[:8]
        outputs = r[8:]
        out = phenotype(*inputs)
        #error += abs(out[0]-outputs[0]) + abs(out[1]-outputs[1])
        
        y1error += abs(out[0]-outputs[0])
        y2error += abs(out[1]-outputs[1])

        
    phenotype.parcialfit = (y1error/len(src),y2error/len(src))
    return sum(phenotype.parcialfit)

if __name__ == '__main__':
    #log.setLevel(logging.DEBUG)
    #random.seed(1225123451234)
    #import pandas
    #bedata = pandas.read_excel('datafiles/ENB2012_data.xlsx')
    #enb= bedata.values
    import pickle
    enb = pickle.load(open('datafiles/ENB2012_data.pkl'))
    #random.shuffle(enb)
    #thethird = int(.7*len(enb))
    #traindata = enb[:thethird]
    #testdata = enb[thethird:]
    evalf = partial(evaluate)
    cfg = loadconfig(parsecmd())
    p = EnergyProb(evalf)
    edw = EvoDevoWorkbench(cfg,p)

    genvalues = []
    folds = createfolds(enb)
    for i in range(len(folds)):
        setcurrent(i, folds)
        edw.run(terminate = (lambda x,y: x < 1e-6 or y <= 0))
        gr = evaluate(edw.best.phenotype,test=True)
        genvalues.append((gr,edw.best.phenotype.parcialfit))

    edw.runlog.validatelog.critical(genvalues)
    for r in  genvalues:
        print r

    #pr.disable()
    #s = StringIO.StringIO()
    #sortby = 'cumulative'
    #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    #ps.print_stats()
    #print s.getvalue()

    #f = open('genome.save','w')
    #f.write(edw.best.phenotype.code.bin)
    #zzf.close

    #plot_ = bindparams(edw.arnconfig, plotindividual)
    #plot_(edw.best.phenotype)
    #genresult = test(ewd.best.circuit, range(2,20))
    #print 'Generalization: '
    #print genresult
