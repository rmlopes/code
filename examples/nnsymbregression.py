from code.operators import *
from code.utils import *
from code.utils.mathlogic import *
from code.evodevo import EvoDevoWorkbench
from code.arncode import *
import logging
import sys
import math

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def kozapolynomial(inp):
    return pow(inp,2) - 2*pow(inp,4) + pow(inp,6)

def quarticpolynomial(inp):
    return pow(inp,4) + pow(inp,3) + pow(inp,2) + inp

def evaluate(circuit, target, inputs, **kwargs):
    errors = [abs(target(inp) - nn([inp],
                                   circuit.input_weights,
                                   circuit.hidden_weights,
                                   circuit.output_weights)[0])
              for inp in inputs]
    try:
        sum_ = sum(errors)
    except:
        log.warning('Invalid individual: overflow error...')
        return 100

    return 100 if math.isinf(sum_) else sum_


if __name__ == '__main__':
    evalfun = partial(evaluate,
                      target=kozapolynomial,
                      inputs=list(drange(-1,1.1,.1)))
    p = CellProb(evalfun, 1, 1)
    edw = EvoDevoWorkbench(sys.argv[1],p,Cell)
    p.eval_ = bindparams(edw.arnconfig, p.eval_)
    edw.run()

    f = open('genome.save','w')
    f.write(edw.best.phenotype.code.bin)
    f.close