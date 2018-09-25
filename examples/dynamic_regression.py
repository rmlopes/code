from symbregression import quarticpolynomial, kozapolynomial,drange
from code.evodevo import *
from code.gearnet import *
from code.utils.config import *
import code.grammars.grammardb as gdb
from code.operators import *

def evaluate(circuit, targetlist, inputs, **kwargs):
    if len(circuit) < 4:
        return 1e6

    #if kwargs['itercount'] % 100 == 0:
     #   targetpointer += 1

    targetpointer = kwargs['itercount'] / 25
    target = targetlist[targetpointer%len(targetlist)]
    if kwargs['itercount'] % 25 == 0:
        print target
    try:
        tytuples = map(lambda x: (target(x),
                              circuit(x)),
                   inputs)

        targets, outputs = zip(*tytuples)

        #if not test:
        mse =  sum(abs(np.array(targets)-np.array(outputs)))
        #else:
         #   mse = nrms_(0,1,zip(targets,outputs))
    except OverflowError:
        log.warning('Invalid individual: overflow error...')
        return 1e6

    return 1e6 if math.isinf(mse) else mse



if __name__ == '__main__':
    #log.setLevel(logging.INFO)
    evalfun = partial(evaluate,
                      targetlist=[quarticpolynomial,kozapolynomial],
                      inputs=list(drange(-1,1.1,.1)))

    ##Fromlist is needed for classes. Functions may be called directly
    #mod = __import__(sys.argv[3], fromlist=[sys.argv[2]])
    #agentclass = getattr(mod, sys.argv[2])

    p = Prob(evalfun, gdb._galapagos, 'start')

    cfg = loadconfig(parsecmd())
    edw = EvoDevoWorkbench(cfg,p)

    edw.run(terminate = (lambda x,y: x <= 1e-3 or y <= 0))
    #print wrapevaluate(edw.best.phenotype,
     #                  target=kozapolynomial,
      #                 inputs=list(drange(-1,1.1,.1)),
       #                device=edw.device,
        #               test = True)
