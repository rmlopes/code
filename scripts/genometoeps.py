#!/usr/bin/python
import sys
import os
import cPickle as pickle
#from pydot import graph_from_dot_data

fname = sys.argv[1]
tfile = open(fname,'r')
#lines= tfile.readlines()
#outfile = sys.argv[2]

dotstr = pickle.load(tfile)
print dotstr
#g = graph_from_dot_data(dotstr)
#g.write_eps(outfile + '.eps')
