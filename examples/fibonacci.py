import sys
from code.evodevo import *
from code.operators import *
from code.rencode import *
from code.utils.mathlogic import *


class Fibonacci(ReNCoDeProb):
	terms = ['inputs[0]','inputs[1]']
        feedback = True
	def __init__(self, evaluate):
		ReNCoDeProb.__init__(self,evaluate)
		self.labels['inputs[0]'] = '0'
		self.labels['inputs[1]'] = '1'
    

def evaluate(circuit, nelems = 10):
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
    print sys.argv
    edw = EvoDevoWorkbench(sys.argv[1],p,buildcircuit,ReNCoDeAgent)
    edw.run()
    genresult = evaluate(edw.best.circuit, 100)
    print 'Generalization: '
    print genresult

