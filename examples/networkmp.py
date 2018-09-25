import sys
from code.evodevo import *
from code.operators import *
from code.cellnetwork import *
from code.utils.mathlogic import *
from code.extendedarn import bindparams, displayARNresults
from code.utils.config import parsecmd, loadconfig

def testmp(phenotype, intinps):
    #FIXME: n = log_2(len(intinps))
    n = 3
    ok = 0
    for i in intinps:
        inputs = BitStream(uint = i, length = n)
        #print inputs.bin
        normalized = nparray([float(inputs.bin[i])
                              for i in range(n)])
        normalized *= .1
        #normalized += 0.05
        phenotype.simulate(phenotype.simtime,*normalized)
        out = phenotype.getbinaryoutput()
        if out == inputs[1+inputs[0]]:
            ok += 1
    return ok

def evaluatecircuit(phenotype, test = False, **kwargs):
        n = 3
        intinps = range(pow(2,n))
        repeated = list()
        if not test:
            #while len(repeated) < 16:
             #   shuffledinps = copy.deepcopy(intinps)# + intinps[:] + intinps[:] + intinps[:]
              #  random.shuffle(shuffledinps)
               # repeated.extend(shuffledinps[:])
            repeated = intinps + intinps
        else:
            repeated = intinps
        try:
                if kwargs['shuffle']:
                    print 'Shuffling input cases...'
                    random.shuffle(intinps)
        except KeyError: pass

        fit = testmp(phenotype, repeated)
        #if not kwargs['silentmode']:
         #   plotindividual(phenotype,**kwargs)
        return len(repeated) - fit

if __name__ == '__main__':
    #p  = BooleanProb(evaluate)
    evalf = evaluatecircuit
    #mapfun = getoutputp0p1)
    p = CellProb(evalf, 3, 1)
    cfg = loadconfig(parsecmd())
    #edw = EvoDevoWorkbench(cfg,p)

    edw = EvoDevoWorkbench(cfg,p)
    p.eval_ = bindparams(edw.arnconfig, p.eval_)
    edw.run()

    #f = open('genome.save','w')
    #f.write(edw.best.phenotype.code.bin)
    #f.close
    #best =
    #plot_ = bindparams(edw.arnconfig, plotindividual)
    #plot_(edw.best.phenotype)
    #genresult = test(ewd.best.circuit, range(2,20))
    #print 'Generalization: '
    #print genresult
