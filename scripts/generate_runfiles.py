#!/usr/bin/env python
"""Generate the run scripts (SGE) for a set of experiments
"""

import csv
import sys
import argparse

body = '''#!/bin/bash
#$ -N %s
#$ -cwd
#$ -t 1-%s

python -m %s %s configfiles/baseline.cfg
'''

options = {'dm5':['--initdm 5',
                  '-o True'],
           'rnd5':['--initdm 5',
                   '--agent code.rencode.RndAgent']}


parser = argparse.ArgumentParser()
parser.add_argument("name", help="a base name for the experiment set.")
parser.add_argument("module", help="the module of the problem to be run.")
parser.add_argument("numruns", help="the number of runs for each experiment")
args = parser.parse_args()

f = open('run-'+args.name,'w')
f.write(body%(args.name, args.numruns,args.module, ''))
f.close()

for k,options in options.items():
    f = open('run-'+args.name + '-' + k, 'w')
    optionstr=""
    for pair in options:
        optionstr += pair + " "
    f.write(body%(args.name+"-"+k, args.numruns, args.module, optionstr))
f.close()