import sys
from code.evodevo import *
from code.operators import *
from code.cellcode import *
from code.utils.mathlogic import *
from code.extendedarn import bindparams, displayARNresults

def test(phen, rangelist ):
    ok = 0
    for i in rangelist:
        r = evaluate(phen,i)
        if r == 0:
            ok += 1
    return len(rangelist) - ok

trainstream = '0101100111001100010111010111000010110100100010101'

def evaluate(phenotype, nbits = 3, **kwargs):
        mapfun = getbinaryoutput
        try:
                mapfun = kwargs['mapfun']
        except KeyError: pass

        ok = 0
        #for cur in range(pow(2,nbits)):
        #t = BitStream(uint=cur,length=nbits)
        t = BitStream(bin=trainstream)
        for i in range(len(t.bin)):
            #print i
            norm = [int(t.bin[i])*0.1]
            phenotype.nstepsim(kwargs['simtime'],*norm)
            #print 'OUT: ', out
            out = mapfun(phenotype, **kwargs)
            if bool(out) == bool(not (t[:i+1].count(1)%2)):
                ok += 1
            #print '%s: %i (%i)' % (t[:i+1].bin, out, ok)

        #return  pow(2,nbits) - ok
        return len(trainstream) - ok

if __name__ == '__main__':
    p = CellProb(evaluate, 1, 1)
    edw = EvoDevoWorkbench(sys.argv[1],p,Cell)
    p.eval_ = bindparams(edw.arnconfig, p.eval_)
    edw.run()

    f = open('genome.save','w')
    f.write(edw.best.phenotype.code.bin)
    f.close

    #plot_ = bindparams(edw.arnconfig, plotindividual)
    #plot_(edw.best.phenotype)
    genresult = test(edw.best.phenotype, range(2,20))
    print 'Generalization: '
    print genresult
