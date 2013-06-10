from symbregression import *

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
                #for gearnet
                self.grammar = grammardb._fixed_grammar_b
                self.start = ['expr']

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

    edw.run(terminate = (lambda x,y: x <= 1e-6 or y <= 0))
    edw.best.fitness = wrapevaluate(edw.best.phenotype, target=Keijzer6,
                                    inputs=list(range(1,121,1)),
                                    device=edw.device,
                                    test = True)
    edw.runlog.step(edw.parents,numevals=edw.numevals)
