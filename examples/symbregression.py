import code
from code.operators import *
from code.utils import *
from code.utils.mathlogic import *
from code.evodevo import EvoDevoWorkbench,  Problem
from code.utils.config import parsecmd, loadconfig
from code.rencode import ReNCoDeProb, defaultnodemap
from code.grammars import grammardb
import logging
import numpy as np
import math

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

def Keijzer6(inp):
    return sum([1.0/i for i in range(1,inp+1,1)])

def evaluate(circuit, target, inputs, test = False,**kwargs):
    if len(circuit) < 4:
        return 1e6

    try:
        tytuples = map(lambda x: (target(x),
                              circuit(x)),
                   inputs)

        targets, outputs = zip(*tytuples)

        if not test:
            mse =  sum(abs(np.array(targets)-np.array(outputs)))
        else:
            mse = nrms_(0,1,zip(targets,outputs))
    except OverflowError:
        log.warning('Invalid individual: overflow error...')
        return 1e6

    return 1e6 if math.isinf(mse) else mse

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

def mse_orig(a,b,tytuples):
    return (1.0/len(tytuples)) * sum([pow((a + b*y) - t,2)
                                   for t,y in tytuples])

def mse_(a,b,tytuples):
    return sum([abs((a + b*y) - t)
                for t,y in tytuples])

def nrms_(a,b,tytuples):
    mse = mse_orig(a,b,tytuples)
    n = len(tytuples)
    targets, outputs = zip(*tytuples)
    return 100 * sqrt(mse) /(max(targets)-min(targets))


def evaluatekeijzer(circuit, target, inputs, printY = False):
        if len(circuit) == 0:
            return 1e6
        #mod = __import__(code, fromlist=['evaluatecircuit'])
        #agentclass = getattr(mod, 'evaluatecircuit')
        try:
            ty_tuples = map(lambda x: (target(x),
                                   circuit(x)),
                        inputs)
            targets, outputs = zip(*ty_tuples)

            if printY:
                print "T: ", targets
                print "Y: ", outputs
                a = circuit.a
                b = circuit.b

            avgtarget = sum(targets)/float(len(targets))
            avgoutput = sum(outputs)/float(len(outputs))

            if not printY:
                bdenom = sum([pow(out - avgoutput,2)
                          for out in outputs])


                if bdenom == 0:
                    return 1e6
                b = float(sum([(t - avgtarget) * (out - avgoutput)
                               for t, out in ty_tuples])) / bdenom
                #if b == 1:
                 #   return 1e6
                a = avgtarget - b * avgoutput
                #if a == 0:
                #    return 1e6
                mse =  mse_(a,b,ty_tuples)
            else:
                mse = nrms_(a,b,ty_tuples)

            #mse = 1.0/len(targets) * sum([pow(a + b*y - t,2)
             #                  for t,y in ty_tuples])
            #n = len(ty_tuples)
            #nrms = 100 * sqrt((n/(n-1))*shape_mse) / numpy.std(targets)
            if printY:
                print 'a= ',circuit.a
                print 'b= ',circuit.b
            else:
                circuit.a = a
                circuit.b = b
        except OverflowError:
            return 1e6

        return mse


def wrapevaluate(circuit, target, inputs, device, test = False):
    if device.__name__ == 'code.rencode' or device.__name__ == 'code.gearnet':
        return evaluatekeijzer(circuit, target, inputs, test)
    else:
        if not test:
            bestfit = 1e6
            outidx = 0
            for i in range(len(circuit.products)):
                circuit.output_idx = i
                fit = evaluatekeijzer(circuit,target, inputs, test)
                if fit < bestfit:
                    bestfit = fit
                    outidx = i
            circuit.output_idx = outidx
        else:
            print 'output index is ', circuit.output_idx
            bestfit = evaluatekeijzer(circuit,target, inputs, test)
        return bestfit


if __name__ == '__main__':
    #log.setLevel(logging.INFO)
    evalfun = partial(evaluate,
                      target=kozapolynomial,
                      inputs=list(drange(-1,1.1,.1)))

    ##fromlist is needed for classes. Functions may be called directly
    #mod = __import__(sys.argv[3], fromlist=[sys.argv[2]])
    #agentclass = getattr(mod, sys.argv[2])

    p = ReNCoDeProb(evalfun)

    cfg = loadconfig(parsecmd())
    edw = EvoDevoWorkbench(cfg,p)

    edw.run(terminate = (lambda x,y: x <= 1e-3 or y <= 0))
    print wrapevaluate(edw.best.phenotype,
                       target=kozapolynomial,
                       inputs=list(drange(-1,1.1,.1)),
                       device=edw.device,
                       test = True)
