from symbregression import *
import code.gearnet as gnet
import code.grammars.grammardb as gdb

if __name__ == '__main__':
    #log.setLevel(logging.INFO)
    evalfun = partial(evaluate,
                      target=quarticpolynomial,
                      inputs=list(drange(-1,1.1,.1)))

    ##Fromlist is needed for classes. Functions may be called directly
    #mod = __import__(sys.argv[3], fromlist=[sys.argv[2]])
    #agentclass = getattr(mod, sys.argv[2])

    p = gnet.Prob(evalfun,gdb._galapagos, 'start')

    cfg = loadconfig(parsecmd())
    edw = EvoDevoWorkbench(cfg,p)

    edw.run(terminate = (lambda x,y: x <= 1e-3 or y <= 0))
    #print wrapevaluate(edw.best.phenotype,
     #                  target=kozapolynomial,
      #                 inputs=list(drange(-1,1.1,.1)),
       #                device=edw.device,
        #               test = True)
