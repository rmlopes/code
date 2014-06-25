#!/usr/bin/python
import sys
import os
import cPickle as pickle
import pydot
from pydot import graph_from_dot_data

treefname = sys.argv[1]
tfile = open(treefname,'r')
lines= tfile.readlines()
filename = treefname.split('.')[0]
for i in range(0,len(lines),2):
    dotstr = pickle.loads(lines[i]+lines[i+1])
    g = graph_from_dot_data(dotstr)
    g.set_splines('ortho')
    print dotstr
    g.write_png(filename + "_" + str((len(lines)-i)/2) + '.png')
