import sys
from code.evodevo import *
from code.operators import *
from code.rencode import *
from code.utils.mathlogic import *
#from code.epicode import LocalSearch,EpiCoDeAgent
from code.utils.config import *

class Hexagons(ReNCoDeProb):
    #funs = ['add_','sub_','mul_','div_']
    terms = ['inputs[0]']
    feedback = True
    def __init__(self, evaluate):
        ReNCoDeProb.__init__(self,evaluate)
        self.labels['inputs[0]'] = '1'


def evaluate(ind, nelems = 100, **kwargs):
    circuit = ind.getcircuit()
    if not any([cnodeinput < 0
                for c in circuit
                for cnodeinput in c[2]]):
        return nelems + 1
    ok = 0
    resultdict = dict(zip([c[0] for c in circuit],
                          [1]*len(circuit)))
    #print nelems
    results = []
    expected = [ (2*n*(2*n - 1))/2 for n in range(1,nelems+1)]
    for cur in range(nelems):
        #n = cur + 1
        results.append(evaluatecircuit(circuit, regressionfun,
                                 resultdict, *[1]))
        #expected = (2*n*(2*n - 1))/2.0
        #print n, expected
        #print result, expected
        #if result == expected:
           # ok += 1
        #else:
          #  break

    #print results
    suc = map(lambda t: 1 if t[0]==t[1] else 0, zip(expected,results))
    return  nelems - sum(suc)

if __name__ == '__main__':
    p  = Hexagons(evaluate)
    #print sys.argv
    #random.seed(1234*int(os.getenv('SGE_TASK_ID')))
    cfg = loadconfig(parsecmd())
    edw = EvoDevoWorkbench(cfg,p)
    edw.run()
    genresult = evaluate(edw.best.phenotype, 1000)
    print 'Generalization: ', genresult
    edw.runlog.validatelog.critical(genresult)
