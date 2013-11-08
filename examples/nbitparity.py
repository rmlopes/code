import sys
from code.evodevo import *
from code.operators import *
from code.rencode import *
from code.utils.mathlogic import *
from code.utils.config import *

class BooleanProb(ReNCoDeProb):
	labels = {'and_':'AND',
                  'or_':'OR',
                  'nand':'NAND',
                  'nor':'NOR',
                  'inputs[0]':'INP'}
        funs = [ 'and_', 'or_', 'nand', 'nor','and_','or_','nand','nor']
        feedback = True


def test(circuit, rangelist ):
    ok = 0
    for i in rangelist:
        r = evaluate(circuit,i)
        if r == 0:
            ok += 1

    return len(rangelist) - ok

def evaluate(pheno, nbits = 3, **kwargs):
        circuit = pheno.getcircuit()
	if not any([cnodeinput < 0
		    for c in circuit
		    for cnodeinput in c[2]]):
		return pow(2,nbits)
	ok = 0
	for cur in range(pow(2,nbits)):
		t = BitStream(uint=cur,length=nbits)
		resultdict = dict(zip([c[0] for c in circuit],[1]*len(circuit)))
		for i in list(t.bin[:]):
                    result = evaluatecircuit(circuit, regressionfun,
                                             resultdict, i)
                if int(result) == int(t.count(1)%2):
			ok += 1

        return  pow(2,nbits) - ok

if __name__ == '__main__':
    p  = BooleanProb(evaluate)
    cfg = loadconfig(parsecmd())
    edw = EvoDevoWorkbench(cfg,p)
    edw.run()
    genresult = test(edw.best.phenotype, range(2,20))
    print 'Generalization: ', genresult
    edw.runlog.validatelog.critical(genresult)
