import sys,os
from code.evodevo import *
from code.operators import *
from code.cellcode import *
from code.utils.mathlogic import *
from code.extendedarn import bindparams, displayARNresults
from code.utils.config import parsecmd, loadconfig

def testmp2(phenotype, intinps):
    #FIXME: n = log_2(len(intinps))
    n = 3
    ok = 0
    #outidx = phenotype.output_idx
    for i in intinps:
        inputs = BitStream(uint = i, length = n)
        #print inputs.bin
        normalized = nparray([float(inputs.bin[i])
                              for i in range(n)])
        normalized *= .15
        phenotype.nstepsim(phenotype.simtime,*normalized)
        diffs = (phenotype.effectorhist[:,-1] -
                 phenotype.effectorhist[:,-2])
        indexdiff = zip(phenotype.effectorproms,diffs)
        indexdiff.sort(key=lambda x: x[1],reverse=True)
        winner = indexdiff[0][0]
        #print 'winner is: ', winner
        if winner < len(phenotype.code) / 4.0:
            out = reduce(and_, inputs)
            #out = 0
        elif len(phenotype.code) / 4.0 <= winner < len(phenotype.code) / 2.0:
            out = reduce(nand, inputs)
        elif len(phenotype.code) / 2.0 <= winner < len(phenotype.code) * .75:
            out = reduce(or_, inputs)
        else:
            out = reduce(nor, inputs)
            #out = 1
        #print 'OUT: ', out
        #out = 0 if out <= 0 else 1
        if out == inputs[1+inputs[0]]:
            ok += 1
    return ok

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
            intinps = intinps[:] #+ intinps[:]
        else:
            random.shuffle(intinps)
        try:
                if kwargs['shuffle']:
                    print 'Shuffling input cases...'
                    random.shuffle(intinps)
        except KeyError: pass

        bestfit = 0
        bestout = 0
        orig_state = copy.deepcopy(phenotype.ccs)
        #print 'ORIGINAL: ',orig_state
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

def evaluatecircuit2(phenotype, test = False, **kwargs):
        n = 3
        intinps = range(pow(2,n))
        repeated = list()
        if not test:
            while len(repeated) < 16:
                shuffledinps = copy.deepcopy(intinps)# + intinps[:] + intinps[:] + intinps[:]
                random.shuffle(shuffledinps)
                repeated.extend(shuffledinps[:])
                #repeated.extend(intinps)
        else:
            repeated = intinps
        try:
                if kwargs['shuffle']:
                    print 'Shuffling input cases...'
                    random.shuffle(intinps)
        except KeyError: pass
        try:
                if kwargs['forceinputs']:
                    repeated = kwargs['forceinputs']
        except KeyError: pass

        fit = testmp2(phenotype, repeated)
        #if not kwargs['silentmode']:
         #   plotindividual(phenotype,**kwargs)
        return len(repeated) - fit

if __name__ == '__main__':
    #p  = BooleanProb(evaluate)
    evalf = evaluatecircuit2
    #mapfun = getoutputp0p1)
    p = CellProb(evalf, 3, 1)
    cfg = loadconfig(parsecmd())
    #edw = EvoDevoWorkbench(cfg,p)

    edw = EvoDevoWorkbench(cfg,p)
    p.eval_ = bindparams(edw.arnconfig, p.eval_)
    edw.run()



    f = open(os.environ['JOB_NAME'] + '_' + os.environ['SGE_TASK_ID']+'.save','w')
    f.write(edw.best.phenotype.code.bin)
    f.close
    #best =
    #plot_ = bindparams(edw.arnconfig, plotindividual)
    #plot_(edw.best.phenotype)
    #genresult = test(ewd.best.circuit, range(2,20))
    #print 'Generalization: '
    #print genresult
