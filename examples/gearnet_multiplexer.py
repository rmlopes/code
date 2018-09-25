from symbregression import *
import code.gearnet as gnet
import code.grammars.grammardb as gdb
from code.operators import *
from code.utils.mathlogic import *

def testmp(phenotype, intinps):
    #FIXME: n = log_2(len(intinps))
    n = int(math.log(len(intinps),2))
    address = int(math.floor(math.log(n,2)))
    #print address
    ok = 0
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
        out = phenotype(*inputs)
        #except:
         ##       print phenotype
           #     print inputs.bin
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
        ok = testmp(phenotype, intinps)
        return len(intinps) - ok

if __name__ == '__main__':
    #p  = BooleanProb(evaluate)
    evalf = partial(evaluate, nbits=11)
    p = gnet.Prob(evalf,gdb._boolean, 'start')
    cfg = loadconfig(parsecmd())
    edw = EvoDevoWorkbench(cfg,p)
    edw.run()

    #f = open('genome.save','w')
    #f.write(edw.best.phenotype.code.bin)
    #zzf.close

    #plot_ = bindparams(edw.arnconfig, plotindividual)
    #plot_(edw.best.phenotype)
    #genresult = test(ewd.best.circuit, range(2,20))
    #print 'Generalization: '
    #print genresult
