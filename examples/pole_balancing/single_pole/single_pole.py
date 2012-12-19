# ******************************************************** #
# Adapted from neat-py.examples.pole_balancing.single_pole #
# @url http://code.google.com/p/neat-python/               #
# Rui L. Lopes (c)2012                                     #
# ******************************************************** #
# Single pole balancing experiment                         #
# ******************************************************** #
#from neat import config, population, chromosome, genome2, visualize
#from neat import nn
import cPickle as pickle
import math, random
import numpy as np
from code.evodevo import *
from code.operators import *
from code.cellcode import *
from code.utils.mathlogic import *
from code.extendedarn import bindparams, displayARNresults
import sys

rstate = random.getstate()
save = open('rstate','w')
pickle.dump(rstate, save)
save.close()

#dumped = open('rstate','r')
#rstate = pickle.load(dumped)
#random.setstate(rstate)
#dumped.close()
GRAVITY = 9.8
MASSCART = 1.0
MASSPOLE = 0.1
TOTAL_MASS = (MASSPOLE + MASSCART)
LENGTH = 0.5    # actually half the pole's length
POLEMASS_LENGTH = (MASSPOLE * LENGTH)
FORCE_MAG = 10.0
TAU = 0.02  # seconds between state updates
FOURTHIRDS = 1.3333333333333
TWELVE_DEGREES = 0.2094384 #radians
ONEANDHALF_DEGREES = 0.02618 #radians

def cart_pole(net_output, x, x_dot, theta, theta_dot):
    ''' Directly copied from Stanley's C++ source code '''
    if net_output > .0:
        force = FORCE_MAG
    else:
        force = -FORCE_MAG

    costheta = math.cos(theta)
    sintheta = math.sin(theta)

    temp = (force + POLEMASS_LENGTH * theta_dot * theta_dot * sintheta)/ TOTAL_MASS

    thetaacc = (GRAVITY*sintheta - costheta*temp)\
               /(LENGTH * (FOURTHIRDS - MASSPOLE * costheta * costheta/TOTAL_MASS))

    xacc  = (temp - POLEMASS_LENGTH * thetaacc * costheta) / TOTAL_MASS

    #Update the four state variables, using Euler's method
    x         += TAU * x_dot
    x_dot     += TAU * xacc
    theta     += TAU * theta_dot
    theta_dot += TAU * thetaacc

    return x, x_dot, theta, theta_dot

def evaluate_individual(phenotype, test = None, **kwargs):

    try:
        outidx = phenotype.output_idx
    except AttributeError:
        outidx = 0
   #radians
    #Nicolau et al. 2010
    if not test:
        num_steps = 120000 #10**5
    else:
        num_steps = 1000

    #for chromo in population:

        #net = nn.create_phenotype(chromo)
        # initial conditions (as used by Stanley)
    if not test:
        x         = random.randint(0, 4799)/1000.0 - 2.4
        x_dot     = random.randint(0, 1999)/1000.0 - 1.0
        theta     = random.randint(0,  399)/1000.0 - 0.2
        theta_dot = random.randint(0, 2999)/1000.0 - 1.5
        print "Initial conditions:"
        print "%2.4f   %2.4f   %2.4f   %2.4f" %(x, x_dot, theta, theta_dot)
    else:
        x, x_dot, theta, theta_dot = test
        #x = 0.0
        #x_dot = 0.0
        #theta = 0.0
        #theta_dot = 0.0

    fitness = 0
    leftlim =  -(1.0 / kwargs['samplerate'])
    last = .0
    for trials in xrange(num_steps):

        # maps into [0,1]
        # RL-NOTE: x_dot and theta_dot may go outside the defined
        # boundaries: how to normalize then?
        # modified to match Nicolau et al 2010 description.
        #inputs = np.array([(x + 2.4)/4.8,
        #                   (x_dot + 0.75)/1.5,
        #                   (theta + twelve_degrees)/0.41,
        #                   (theta_dot + 1.0)/2.0],
        #                  dtype = float)
        inputs = np.array([(x + 2.4)/4.8,
                           (x_dot + 1.0)/2.0,
                           #(theta + twelve_degrees)/0.41,
                           (theta + TWELVE_DEGREES)/( 2*TWELVE_DEGREES),
                           (theta_dot + 1.5)/3.0],
                          dtype = float)
        # maps into cc [0,.1]
        inputs *= .1
        for i in range(len(inputs)):
                if inputs[i] > .1:
                        inputs[i] = .1
                elif inputs[i] < 0:
                        inputs[i] = .0

        # a normalizacao so acontece para estas condicoes iniciais
        # nada garante que a evolucao do sistema leve a outros
        # valores de x, x_dot e etc...
        #print inputs
        phenotype.nstepsim(kwargs['simtime'],*inputs)
        #output =
        # np.sum(np.gradient(phenotype.effectorhist[outidx][leftlim-1:]))
        #NOTE: only valid for samplerate =1 .0
        output = (phenotype.effectorhist[outidx][-1] -
                  phenotype.effectorhist[outidx][-2])
        #print output
        # Apply action to the simulated cart-pole
        #if no cc change then repeat last
        if abs(output) <= 1e-20:
            output = last
        x, x_dot, theta, theta_dot = \
            cart_pole(output, x, x_dot, theta, theta_dot)

        # Check for failure.  If so, return steps
        # the number of steps indicates the fitness: higher = better
        fitness += 1
        last = output
        #RL-NOTE: original does not check for the speed boundaries
        # Problem description in Nicoulau et al. 2010 shows closed
        # intervals (>= should be >)
        if (abs(x) >= 2.4 or abs(theta) >= TWELVE_DEGREES):
            #or abs(x_dot) > 1 or abs(theta_dot) > 1.5):
            # the cart/pole has run/inclined out of the limits
            break

    #Fitness as defined in (Nicolau et al., 2010)
    #adapted to minimize untill zero
    #plotindividual(phenotype,**kwargs)
    return num_steps/float(fitness) - 1

if __name__ == "__main__":
    evalf = partial(evaluate_individual)
    p = CellProb(evalf, 4, 1)
    edw = EvoDevoWorkbench(sys.argv[1],p,Cell)
    p.eval_ = bindparams(edw.arnconfig, p.eval_)
    edw.run()

    f = open('genome.save','w')
    f.write(edw.best.phenotype.code.bin)
    f.close