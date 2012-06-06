#from code import operators
import code
from code.operators import *
from math import *
from code.utils.mathlogic import *
from code.rencode import *
from code.evodevo import *
from random import sample
import matplotlib.mlab as mlab
from numpy import array as nparray
import numpy

numclasses = 5
allclasses = nparray([int(c) 
                      for c in open('datafiles/MIREXclasses.txt').readlines()])
allfeatures =nparray([map(lambda t: float(t), l.split(',')) 
                      for l in open('datafiles/FMnorm.txt').readlines()])

pcafeats = mlab.PCA(allfeatures)

#print pcafeats.fracs
#print pcafeats.Y[0,:]

projected = pcafeats.project(allfeatures,0.03)
print projected[:1]
numpy.save('datafiles/projectedfeat-03', projected)
