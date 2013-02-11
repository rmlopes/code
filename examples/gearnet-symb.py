from code.operators import *
from code.utils import *
from code.utils.mathlogic import *
from code.evodevo import EvoDevoWorkbench
from code.gearnet import *
from code.grammars import grammardb
import logging
import sys


log = logging.getLogger(__name__)

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def kozapolynomial(inp):
    return pow(inp,2) - 2*pow(inp,4) + pow(inp,6)

def quarticpolynomial(inp):
    return pow(inp,4) + pow(inp,3) + pow(inp,2) + inp

def evaluate(phenotype, target, inputs, grm = grammardb._base_grammar):
    circuit = phenotype.circuit
    #print 'duh'
    #print circuit
    for k in grm.iterkeys():
        #print 'key:', k
        if circuit.find(k) != -1:
            return 1e9

    errors = [abs(target(inp) -
                  eval(circuit))
              for inp in inputs]
    try:
        sum_ = sum(errors)
    except:
        log.warning('Invalid individual: overflow error...')
        return 100

    return 100 if math.isinf(sum_) else sum_


if __name__ == '__main__':
    log.setLevel(logging.DEBUG)
    evalfun = partial(evaluate,
                      target=kozapolynomial,
                      inputs=list(drange(-1,1.1,.1)),
                      grm = grammardb._fixed_grammar_b)
    p = Prob(evalfun, grammardb._fixed_grammar_b, ['expr'])

    try:
        agentclass = getattr(__import__('__main__'),sys.argv[2])
    except:
        print 'Loading default agent class'
        agentclass = DMAgent

    edw = EvoDevoWorkbench(sys.argv[1],p,agentclass)
    edw.run()
    print edw.best.phenotype.circuit
