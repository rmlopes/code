#from code import operators
from code.evodevo import *
from code.operators import *
from math import *
from code.utils.mathlogic import *
from code.rencode import *
from random import sample
#import matplotlib.mlab as mlab
from numpy import array as nparray
import numpy
import logging

log = logging.getLogger('mirex')

numclasses = 5
allclasses = nparray([int(c)
                      for c in open('datafiles/MIREXclasses.txt').readlines()])
allfeatures =nparray([map(lambda t: float(t), l.split(','))
                      for l in open('datafiles/FMnorm.txt').readlines()])

projected = numpy.load('datafiles/projectedfeat-01.npy')
defaultdiscard = (0,0,1,1)

def sigmoid(x):
    try:
        return 1 / (1 + exp(-x))
    except OverflowError:
        return 0

def evaluate(circuit, test = False, relax = False):
    if len(circuit) < 4 and not relax:
        return defaultdiscard

    mainmod = __import__('__main__')
    #workingset = globals()['testset'] if test else globals()['trainset']
    workingset = getattr(mainmod,
                         'testset') if test else getattr(mainmod,'trainset')
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    for c,feats,origc in workingset:
        #try:
        result = evaluatecircuit(circuit, nnlikefun,
                                 dict(), *feats)
        if isnan(result) or isinf(result):
            return defaultdiscard
        #except:
        #   return defaultdiscard
        try:
            result = 1 if tanh(result) > 0 else 0
        except OverflowError:
            return defaultdiscard
        #if test:
        #print c, result
        if c and result :
            tp += 1
        elif c and not result:
            fn += 1
        elif not c and not result:
            tn += 1
        elif not c and result:
            fp += 1

            #print sample(results,10)
    return tp, tn, fp, fn

def fmeas_eval(circuit, test = False):
    results = evaluate(circuit, test)
    tp, tn, fp, fn = results
    if test:
        log.critical("%i\t%i\t%i\t%i", *results)
    try:
        f = fmeasure(tp, fp, fn)
    except ZeroDivisionError:
        f = -1

    return (1 - f)

def evalaccuracy(circuit):
    results = evaluate(circuit)
    tp, tn, fp, fn = results

    return 1 - (tp + tn) / float(tp+tn+fp+fn)

globalconfusion = nparray([[0]*numclasses
                           for i in range(numclasses)])

def getfalsepos(confusion,target):
    return sum([confusion[i][target] for i in range(len(confusion))
                if i != target])

def getfalseneg(confusion,target):
    return sum([confusion[target][i] for i in range(len(confusion))
                if i != target])

def combinedeval(classifiers, testset):
    '''Highest value identifies the class'''
    rl = list(range(1,len(classifiers)+1))
    results = dict(zip(rl,
                       [0,0,0,0,0]))
    confusion = nparray([[0] * len(classifiers)
                         for i in range(len(classifiers))])
    for binc,feats,annot in testset:
        c, indx = annot
        outputs = [sigmoid(evaluatecircuit(circuit, nnlikefun, dict(), *feats))
                   for circuit in classifiers]

        zipouts = zip(range(1,6),outputs)
        pred = max(zipouts, key=lambda x: x[1])[0]

        confusion[c-1][pred-1] += 1
        globalconfusion[c-1][pred-1] += 1

    print confusion

    return map(lambda x: 0 if isnan(x) else x,
               [fmeasure(confusion[c][c],
                    getfalsepos(confusion,c),
                    getfalseneg(confusion,c))
                for c in range(len(classifiers))])

def runfold(i, datasets, edw):
    print "Running fold %i"%i
    mainmod = __import__('__main__')
    testset = []
    circuits=[]
    for samples in datasets:
        trainset = []
        for j in range(len(samples)):
            if j != i:
                trainset.extend(samples[j])
            if j == i:
                testset.extend(samples[j])
        setattr(mainmod,'trainset',trainset)
        edw.run()
        circuits.append(edw.best.phenotype)

    used = set()
    auxts = []
    for i in range(len(testset)):
        if testset[i][2][1] not in used:
            used.add(testset[i][2][1])
            auxts.append(testset[i])

    return combinedeval(circuits, auxts)



def createdataset(data, target):
    zipped1 = map(lambda x: (1,x[1],(x[0],x[2])) if x[0] == target else (0,x[1],(x[0],x[2])),
                  data)
    cset = filter(lambda x: x[0] == 1, zipped1)
    otherset = filter(lambda x: x[0] == 0, zipped1)
    workset = cset + random.sample(otherset, len(cset))
    random.shuffle(workset)
    foldsize = int(len(workset)/numfolds)
    folds = []
    for i in range(numfolds-1):
        start = i*foldsize
        folds.append(workset[start:start+foldsize])
    folds.append(workset[(numfolds-1)*foldsize:])
    return folds

if __name__ == '__main__':
    import sys
    try:
        if sys.argv[2]  == 'pca': #nopca
            zipped = zip(allclasses, projected, range(len(allclasses)))
        elif sys.argv[2] == 'original':
            zipped = zip(allclasses, allfeatures, range(len(allclasses)))

        numfolds = int(sys.argv[3])
        print "NumFeatures: %i "% (len(zipped[0][1]))

        fit = sys.argv[4]
    except:
        print "usage: python -m examples.mirex <config> <feat:pca|original> <numfolds> <fit:acc|f1>"
        print "example: python -m examples.mirex test.cfg pca 3"
        exit(0)

    datasets = []
    for i in range(1,6):
        datasets.append(createdataset(zipped,i))

    if fit == 'f1':
        problem  = ClassifProb(fmeas_eval,len(zipped[0][1]))
    else:
        problem = ClassifProb(evalaccuracy,len(zipped[0][1]))

    #problem.funs = ['add_','sub_','mul_','div_','min_','max_']
    edw = EvoDevoWorkbench(sys.argv[1],problem,buildcircuit,ReNCoDeAgent)

    vresults = 0
    for i in range(numfolds):
        v = runfold(i, datasets, edw)
        cv = sum(v)/len(v)
        log.critical("%f\t%f\t%f\t%f\t%f\t%f"% (tuple(v)+(cv,)))
        vresults += cv

    print 'CONFUSION:'
    for i in range(5):
        print "%i\t%i\t%i\t%i\t%i" % tuple(globalconfusion[i,:])

    globalresults =  map(lambda x: 0 if isnan(x) else x,
                         [fmeasure(globalconfusion[c][c],
                                   getfalsepos(globalconfusion,c),
                                   getfalseneg(globalconfusion,c))
                          for c in range(len(classifiers))])

    print globalresults
    log.critical("%f", sum(globalresults)/len(globalresults))
    log.critical("%f", vresults/numfolds)
