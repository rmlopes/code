from code.operators import *
from code.utils import *
from code.utils.mathlogic import *
from code.evodevo import EvoDevoWorkbench
from code.rencode import *
import logging
import sys

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def kozapolynomial(inp):
    return pow(inp,2) - 2*pow(inp,4) + pow(inp,6)

def quarticpolynomial(inp):
    return pow(inp,4) + pow(inp,3) + pow(inp,2) + inp

def Keijzer6(inp):
    return sum([i for i in range(1,inp+1,1)])

def evaluate(circuit, target, inputs):
    if len(circuit) < 4:
        return 100
    errors = [abs(target(inp) -
                  evaluatecircuit(circuit, regressionfun,dict(),inp))
              for inp in inputs]
    try:
        sum_ = sum(errors)
    except:
        log.warning('Invalid individual: overflow error...')
        return 100

    return 100 if math.isinf(sum_) else sum_

def evaluatennlike(circuit, target, inputs):
    if len(circuit) < 4:
        return 1e6
    errors = [abs(target(inp) -
                  evaluatecircuit(circuit, nnlikefun,dict(),inp))
              for inp in inputs]
    try:
        sum_ = sum(errors)
    except:
        log.warning('Invalid individual: overflow error...')
        return 1e6

    return 1e6 if math.isinf(sum_) else sum_


class KeijzerProb(ReNCoDeProb):
        #extrafuns = ['exp_','min_','max_','log_','tanh_','sin_','cos_','sinh_','cosh_','tan_']
        funs = ['add_','mul_','sqrt_','cos_','reciprocalsum','log_']
        terms = ['inputs[0]']#,'tanh_','sin_','cos_','sinh_','cosh_','tan_']
        def __init__(self, evaluate, **kwargs):
                self.labels = None
                ReNCoDeProb.__init__(self,evaluate, **kwargs)
                #self.terms.extend(["inputs[%i]"%i
                 #                  for i in range(1,numfeat)])
                #self.funs.extend(self.extrafuns)
                self.arity.update(zip(self.funs,[0]*len(self.funs)))
                print self.arity

if __name__ == '__main__':
    #evalfun = partial(evaluate,
     #                 target=kozapolynomial,
      #                inputs=list(drange(-1,1.1,.1)))
    evalfun = partial(evaluatennlike,
                      target=Keijzer6,
                      inputs=list(range(1,51,1)))
    #p = ReNCoDeProb(evalfun)
    p = KeijzerProb(evalfun)
    try:
        agentclass = getattr(__import__('__main__'),sys.argv[2])
    except:
        print 'Loading default agent class'
        agentclass = DMAgent

    edw = EvoDevoWorkbench(sys.argv[1],p,agentclass)
    edw.run()
