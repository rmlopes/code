import sys
from code.evodevo import *
from code.operators import *
from code.rencode2 import *
from code.rencode import ReNCoDeProb, evaluatecircuit, regressionfun
from code.utils.mathlogic import *
from code.extendedarn import bindparams, displayARNresults

class BooleanProb(ReNCoDeProb):
        labels = {'and_':'AND',
                  'or_':'OR',
                  'nand':'NAND',
                  'nor':'NOR',
                  'inputs[0]':'INP'}
        funs = [ 'and_', 'or_', 'nand', 'nor','and_','or_','nand','nor']
        feedback = True
        def __init__(self, nbits, evalf):
                self.ninp = nbits
                self.nout = 1
                self.terms = ['inputs[%i]'%(idx,) for idx in range(nbits)]
                ReNCoDeProb.__init__(self,evalf, printf=printrencode2)

def testmp(phenotype, intinps):
    #FIXME: n = log_2(len(intinps))
    n = int(math.log(len(intinps),2))
    address = int(math.floor(math.log(n,2)))
    #print address
    ok = 0
    outidx = phenotype.output_idx
    for i in intinps:
        inputs = BitStream(uint = i, length = n)
        #print inputs.bin
        #normalized = nparray([float(inputs.bin[i])
         #                     for i in range(n)])
        #normalized *= .1
        #phenotype.nstepsim(phenotype.simtime,*normalized)
        #out = (phenotype.effectorhist[outidx][-1] -
        #       phenotype.effectorhist[outidx][-2])
        try:
                out = evaluatecircuit(phenotype.getcircuit(outidx),
                              regressionfun, dict(),
                              *inputs.bin)
        except:
                print phenotype.getcircuit(outidx)
                print inputs.bin
        #print 'OUT: ', out
        #out = 0 if out <= 0 else 1
        index = BitStream(bin=inputs.bin[:address]).uint
        #print index, inputs.bin[:address]
        #assert(index == inputs[0])
        if out == inputs[address+index]:
            #print out, inputs
            ok += 1
    return ok

def evaluate(phenotype, test = False, nbits = 3, **kwargs):
        n = nbits
        intinps = range(pow(2,n))
        if not test:
                intinps = intinps[:]
        #random.shuffle(intinps)
        try:
                if kwargs['shuffle']:
                    print 'Shuffling input cases...'
                    random.shuffle(intinps)
        except KeyError: pass

        bestfit = 0
        bestout = 0
        #print 'ORIGINAL: ',orig_state
        if not test:
            for eff in range(phenotype.arnet.numeff):
                phenotype.output_idx = eff
                ok = testmp(phenotype, intinps)
                if ok > bestfit:
                    bestfit = ok
                    bestout = eff
                    #print 'best output index is now ',eff
            phenotype.output_idx = bestout
        else:
            print 'output index is ', phenotype.output_idx
            bestfit = testmp(phenotype, intinps)
        return len(intinps) - bestfit

if __name__ == '__main__':
    #p  = BooleanProb(evaluate)
    evalf = partial(evaluate, nbits=11)
    #mapfun = getoutputp0p1)
    p = BooleanProb(11,evalf)
    edw = EvoDevoWorkbench(sys.argv[1],p,RndAgent)
    #p.eval_ = bindparams(edw.arnconfig, p.eval_)
    edw.run()

    #f = open('genome.save','w')
    #f.write(edw.best.phenotype.code.bin)
    #zzf.close

    #plot_ = bindparams(edw.arnconfig, plotindividual)
    #plot_(edw.best.phenotype)
    #genresult = test(ewd.best.circuit, range(2,20))
    #print 'Generalization: '
    #print genresult
