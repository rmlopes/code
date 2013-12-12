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
        feedback = True
        def __init__(self, nbits, evalf, **kwargs):
                self.ninp = nbits
                self.nout = 2
                self.terms = ['inputs[%i]'%(idx,) for idx in range(nbits)]
                ReNCoDeProb.__init__(self,evalf, **kwargs)

def testadder(phenotype, intinps):
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
        #try:
        out = phenotype(*inputs.bin)
                #evaluatecircuit(phenotype.getcircuit(outidx),
                #                                     regressionfun,
                # mergefun,
                #                                    dict(),
                #                                   *inputs.bin, nout = 2)
        #except:
         #       print phenotype.getcircuit(outidx)
          #      print inputs.bin
        r = BitStream(bin='%i%i'%(int(out[0]),int(out[1])))
	expected = BitStream(uint = inputs.count(1), length = n-1)
	ok += (expected ^ r).count(1)
        #print 'OUT: ', out
        #out = 0 if out <= 0 else 1
        #out = reduce(lambda x,k: x + ('1' if k else '0'), out, '')
        #result = BitStream(bin=out).uint
        #if result == inputs.count(1):
            #print out, inputs
        #    ok += 1
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
        r = testadder(phenotype, intinps)
        #return bestfit
        #r = len(intinps) - bestfit
        if r == 0:
            return 0
        else:
            #numtf = float(phenotype.arnet.numtf) + 1.0
            numeff = phenotype.arnet.numeff
            if numeff > 1:
                p = 1.0 - 2/float(numeff)#.05 * abs(numeff-2)
            else:
                p = 1.0#

            return r + .45*p #(.45 / float(numtf)) + (.45 * p)

if __name__ == '__main__':
    #log.setLevel(logging.DEBUG)
    random.seed(1234)
    #p  = BooleanProb(evaluate)
    evalf = partial(evaluate, nbits=3)
    #mapfun = getoutputp0p1)
    cfg = loadconfig(parsecmd())
    p = BooleanProb(3,evalf)
    edw = EvoDevoWorkbench(cfg,p)
    #p.eval_ = bindparams(edw.arnconfig, p.eval_)
    edw.run(terminate = (lambda x,y: x ==0 or y <= 0))
    #edw.run()

    #f = open('genome.save','w')
    #f.write(edw.best.phenotype.code.bin)
    #zzf.close

    #plot_ = bindparams(edw.arnconfig, plotindividual)
    #plot_(edw.best.phenotype)
    #genresult = test(ewd.best.circuit, range(2,20))
    #print 'Generalization: '
    #print genresult
