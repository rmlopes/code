import sys
from code.evodevo import *
from code.operators import *
from code.cellcode import *
from code.utils.mathlogic import *
from code.extendedarn import bindparams, displayARNresults

if __name__ == '__main__':
    #p  = BooleanProb(evaluate)
    evalf = evaluatecircuit
    #mapfun = getoutputp0p1)
    p = CellProb(evalf, 3, 1)
    edw = EvoDevoWorkbench(sys.argv[1],p,Cell)
    p.eval_ = bindparams(edw.arnconfig, p.eval_)
    edw.run()

    f = open('genome.save','w')
    f.write(edw.best.phenotype.code.bin)
    f.close

    #plot_ = bindparams(edw.arnconfig, plotindividual)
    #plot_(edw.best.phenotype)
    #genresult = test(ewd.best.circuit, range(2,20))
    #print 'Generalization: '
    #print genresult
