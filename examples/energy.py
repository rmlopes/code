import sys
from code.utils.config import parsecmd, loadconfig
from code.evodevo import *
from code.operators import *
from code.rencode2 import *
from code.rencode import ReNCoDeProb, evaluatecircuit, regressionfun
from code.utils.mathlogic import *
from code.extendedarn import bindparams, displayARNresults
import logging

log = logging.getLogger(__name__)


class EnergyProb(ReNCoDeProb):
        #extrafuns = ['exp_','min_','max_','log_','tanh_','sin_','cos_','sinh_','cosh_','tan_']
        funs = ['add_','mul_','sqrt_','cos_','reciprocalsum','log_']
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

def evaluate(phenotype, src, **kwargs):
    error = 0.0


    for r in src:
        inputs = r[:8]
        outputs = r[8:]
        out = phenotype(*inputs)
        error += abs(out[0]-outputs[0]) + abs(out[1]-outputs[1])

    return error

if __name__ == '__main__':
    #log.setLevel(logging.DEBUG)
    #random.seed(1234)
    import pandas
    bedata = pandas.read_excel('datafiles/ENB2012_data.xlsx')
    enb= bedata.values
    random.shuffle(enb)
    thethird = int(.7*len(enb))
    traindata = enb[:thethird]
    testdata = enb[thethird:]

    evalf = partial(evaluate, src =  traindata)
    #mapfun = getoutputp0p1)
    cfg = loadconfig(parsecmd())
    p = EnergyProb(evalf)
    edw = EvoDevoWorkbench(cfg,p)
    #p.eval_ = bindparams(edw.arnconfig, p.eval_)
    edw.run(terminate = (lambda x,y: x < 1e-6 or y <= 0))
    genv = evaluate(edw.best.phenotype,testdata)
    edw.runlog.validatelog.critical(genv)    
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
