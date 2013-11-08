import sys
from code.evodevo import *
from code.operators import *
from code.rencode import *
from code.utils.mathlogic import *
#from code.epicode import LocalSearch,EpiCoDeAgent
from code.utils.config import *

class Fibonacci(ReNCoDeProb):
	terms = ['inputs[0]','inputs[1]']
        feedback = True
	def __init__(self, evaluate):
		ReNCoDeProb.__init__(self,evaluate)
		self.labels['inputs[0]'] = '0'
		self.labels['inputs[1]'] = '1'


def evaluate(ind, nelems = 10, **kwargs):
        circuit = ind.getcircuit()
	if not any([cnodeinput < 0
		    for c in circuit
		    for cnodeinput in c[2]]):
		return nelems + 1
	ok = 0
	resultdict = dict(zip([c[0] for c in circuit],
			      [1]*len(circuit)))
	last1 = 0
	last2 = 0
	for cur in range(nelems):
		result = evaluatecircuit(circuit, regressionfun,
					 resultdict, *[0,1])
		expected = last1 + last2
		#print result, expected
		if cur == 0 and result == 0:
			ok += 1
			last1 = 1
		elif cur==1 and result == expected:
			ok += 1
			last2 = 1
		elif cur > 1:
			if result == expected:
				ok += 1
			last2 = last1
			last1 = expected
	return  nelems - ok

if __name__ == '__main__':
    p  = Fibonacci(evaluate)
    #print sys.argv
    random.seed(1234*int(os.getenv('SGE_TASK_ID')))
    cfg = loadconfig(parsecmd())
    edw = EvoDevoWorkbench(cfg,p)
    edw.run()
    genresult = evaluate(edw.best.phenotype, 100)
    print 'Generalization: ', genresult
    edw.runlog.validatelog.critical(genresult)
