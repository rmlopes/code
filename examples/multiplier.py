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

class BooleanProb(ReNCoDeProb):
        labels = {'and_':'AND',
                  'or_':'OR',
                  'nand':'NAND',
                  'nor':'NOR',
                  'inputs[0]':'INP'}
        funs = [ 'and_', 'or_', 'nand', 'nor','and_','or_','nand','nor']
        feedback = False
        def __init__(self, nbits, evalf, **kwargs):
                self.ninp = nbits*2
                self.nout = nbits*2
                self.terms = ['inputs[%i]'%(idx,) for idx in range(nbits*2)]
                ReNCoDeProb.__init__(self,evalf, **kwargs)

def testadder(phenotype, intinps):
        #FIXME: n = log_2(len(intinps))
    n = int(math.log(len(intinps),2))/2
    #print address
    ok = 0
    for i in intinps:
        inputs = BitStream(uint = i, length = n*2)
        a = inputs[:n]
        b = inputs[n:]
        #print inputs.bin
        #normalized = nparray([float(inputs.bin[i])
         #                     for i in range(n)])
        #normalized *= .1
        #phenotype.nstepsim(phenotype.simtime,*normalized)
        #out = (phenotype.effectorhist[outidx][-1] -
        #       phenotype.effectorhist[outidx][-2])
        try:
                out = list(evaluatecircuit(phenotype.getcircuit(),
                                      regressionfun, dict(),
                                      *inputs.bin, nout = 4))
        except:
                print phenotype.getcircuit()
                print inputs.bin
                return -1
        #print 'OUT: ', out
        #out = 0 if out <= 0 else 1
        out = reduce(lambda x,k: x + ('1' if k else '0'), out, '')
        if a.uint * b.uint == BitStream(bin = out).uint:
        #result = int(out0) * 2 + int(out1) * 1
        #if result == inputs.count(1):
            #print out, inputs
            ok += 1
    return ok

def evaluate(phenotype, test = False, nbits = 2, **kwargs):
        intinps = range(pow(2,nbits*2))
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
        bestfit = testadder(phenotype, intinps)
        if bestfit == -1:
                return 1e6
        return len(intinps) - bestfit

if __name__ == '__main__':
    #log.setLevel(logging.DEBUG)
    #random.seed(1234*int(os.getenv('SGE_TASK_ID')))
    #p  = BooleanProb(evaluate)
    evalf = partial(evaluate, nbits=2)
    #mapfun = getoutputp0p1)
    cfg = loadconfig(parsecmd())
    p = BooleanProb(2,evalf)
    edw = EvoDevoWorkbench(cfg,p)
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
