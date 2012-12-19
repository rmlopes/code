import sys
from code.evodevo import *
from code.operators import *
from code.cellcode import *
from code.utils.mathlogic import *
from code.extendedarn import bindparams, displayARNresults

def testmp(phenotype, intinps):
    #FIXME: n = log_2(len(intinps))
    n = 3
    ok = 0
    outidx = phenotype.output_idx
    for i in intinps:
        inputs = BitStream(uint = i, length = n)
        #print inputs.bin
        normalized = nparray([float(inputs.bin[i])
                              for i in range(n)])
        normalized *= .1
        phenotype.nstepsim(phenotype.simtime,*normalized)
        out = (phenotype.effectorhist[outidx][-1] -
               phenotype.effectorhist[outidx][-2])
        #print 'OUT: ', out
        out = 0 if out <= 0 else 1
        if out == inputs[1+inputs[0]]:
            ok += 1
    return ok

#TODO: move this inside the problem
def evaluatecircuit(phenotype, test = False, **kwargs):
        n = 3
        intinps = range(pow(2,n))
        if not test:
                intinps = intinps[:]# + intinps[:]
        #random.shuffle(intinps)
        try:
                if kwargs['shuffle']:
                    print 'Shuffling input cases...'
                    random.shuffle(intinps)
        except KeyError: pass

        bestfit = 0
        bestout = 0
        orig_state = phenotype.ccs[:]
        if not test:
            for eff in range(phenotype.numeff):
                phenotype.reset(orig_state)
                phenotype.output_idx = eff
                ok = testmp(phenotype, intinps)
                if ok > bestfit:
                    bestfit = ok
                    bestout = eff
                    print 'best output index is now ',eff
            phenotype.output_idx = bestout
        else:
            print 'output index is ', phenotype.output_idx
            bestfit = testmp(phenotype, intinps)
        if not kwargs['silentmode']:
            plotindividual(phenotype,**kwargs)
        return len(intinps) - bestfit

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
