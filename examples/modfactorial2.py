import sys
from code.evodevo import *
from code.operators import *
from code.rencode import *
from code.utils.mathlogic import *
#from code.epicode import LocalSearch,EpiCoDeAgent
from code.utils.config import *

class Squares(ReNCoDeProb):
    funs = ['add_','sub_']
    terms = ['inputs[0]']
    feedback = True
    def __init__(self, evaluate):
        ReNCoDeProb.__init__(self,evaluate)
        self.labels['inputs[0]'] = '1'


def evaluate(ind, nelems = 10, k = 2, step = 1, **kwargs):
    circuit = ind.getcircuit()
    if not any([cnodeinput < 0
                for c in circuit
                for cnodeinput in c[2]]):
        return nelems + 1
    ok = 0
    resultdict = dict(zip([c[0] for c in circuit],
                          [1]*len(circuit)))
    outputs = []
    for cur in range(nelems):
        result = evaluatecircuit(circuit, regressionfun,
                                 resultdict, *[1.0])
        if cur < step:
            expected = 1
        else:
            expected = k * outputs[-step]
        #print result, expected
        if result == expected:
            ok += 1
            outputs.append(result)
        else:
            break

    return  nelems - ok

if __name__ == '__main__':
    eval_ = partial(evaluate,
                    k = 2,
                    step = 3)
    p  = Squares(eval_)
    #print sys.argv
    #random.seed(1234*int(os.getenv('SGE_TASK_ID')))
    cfg = loadconfig(parsecmd())
    edw = EvoDevoWorkbench(cfg,p)
    edw.run()
    genresult = evaluate(edw.best.phenotype, nelems = 100)
    print 'Generalization: ', genresult
    edw.runlog.validatelog.critical(genresult)
