import code
from code.operators import *
from code.utils import *
from code.utils.mathlogic import *
from code.evodevo import EvoDevoWorkbench,  Problem
from code.utils.config import parsecmd, loadconfig
from code.rencode import ReNCoDeProb, defaultnodemap
import logging

logging.basicConfig()

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
    return sum([1.0/i for i in range(1,inp+1,1)])

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

#device is the CoDe module, deduced from the agent class full name
def evaluatekeijzer(circuit, target, inputs, device, printY = False):
        if len(circuit) == 0:
            return 1e6
        #mod = __import__(code, fromlist=['evaluatecircuit'])
        #agentclass = getattr(mod, 'evaluatecircuit')

        ty_tuples = map(lambda x: (target(x),
                                   device.evaluatecircuit(circuit,
                                                          device.nnlikefun,
                                                          dict(),x)),
                               inputs)
        targets, outputs = zip(*ty_tuples)

        if printY:
                print "T: ", targets
                print "Y: ", outputs
        avgtarget = sum(targets)/float(len(targets))
        avgoutput = sum(outputs)/float(len(outputs))

        try:
            bdenom = sum([pow(out - avgoutput,2)
                      for out in outputs])
        except OverflowError:
            return 1e6

        if bdenom == 0:
                return 1e6
        b = float(sum([(t - avgtarget) * (out - avgoutput)
                 for t, out in ty_tuples])) / bdenom
        if b == 1:
            return 1e6
        a = avgtarget - b * avgoutput
        if a == 0:
            return 1e6

        shape_mse = 1.0/len(targets) * sum([pow(a + b*y - t,2)
                                      for t,y in ty_tuples])
        return shape_mse + \
        sum([pow(y - t,2) for t,y in ty_tuples])/float(len(ty_tuples))

def wrapevaluate(circuit, target, inputs, device, printY = False, test = False):
    if device.__name__ == 'code.rencode':
        return evaluatekeijzer(circuit, target, inputs, device, printY)
    else:
        if not test:
            bestfit = 1e6
            outidx = 0
            for i in range(len(circuit.products)):
                fit = evaluatekeijzer(circuit.getcircuit(i),
                                      target, inputs, device, printY)
                if fit < bestfit:
                    bestfit = fit
                    outidx = i
            circuit.output_idx = outidx
        else:
            print 'output index is ', circuit.output_idx
            bestfit = evaluatekeijzer(circuit.getcircuit(circuit.output_idx),
                                      target, inputs, device, printY)
        return bestfit

class KeijzerProb(ReNCoDeProb):
        #extrafuns = ['exp_','min_','max_','log_','tanh_','sin_','cos_','sinh_','cosh_','tan_']
        funs = ['add_','mul_','sqrt_','cos_','reciprocalsum','log_']
        terms = ['inputs[0]']#,'tanh_','sin_','cos_','sinh_','cosh_','tan_']
        arity={}
        labels={}
        def __init__(self, evaluate, **kwargs):
                self.labels = None
                ReNCoDeProb.__init__(self,evaluate,**kwargs)
                #self.terms.extend(["inputs[%i]"%i
                 #                  for i in range(1,numfeat)])
                #self.funs.extend(self.extrafuns)
                self.arity.update(zip(self.funs,[0]*len(self.funs)))
                #print self.arity
                self.nout = 1
                self.ninp = 1
                #self.nodemap_ = defaultnodemap

if __name__ == '__main__':
    #log.setLevel(logging.INFO)
    #evalfun = partial(evaluate,
     #                 target=kozapolynomial,
      #                inputs=list(drange(-1,1.1,.1)))

    ##fromlist is needed for classes. Functions may be called directly
    #mod = __import__(sys.argv[3], fromlist=[sys.argv[2]])
    #agentclass = getattr(mod, sys.argv[2])

    evalfun = partial(wrapevaluate,
                      target=Keijzer6,
                      inputs=list(range(1,51,1)))
    p = KeijzerProb(evalfun)

    cfg = loadconfig(parsecmd())
    edw = EvoDevoWorkbench(cfg,p)
    edw.run()
    print wrapevaluate(edw.best.phenotype, target=Keijzer6,
                       inputs=list(range(1,121,1)),
                       device=edw.device,
                       printY = True,
                       test = True)