#!/usr/bin/env python
"""Generate the run scripts (SGE) for a set of experiments
"""

import csv
import sys
import argparse
import experiments
body = '''#!/bin/bash
#$ -N %s
#$ -cwd
#$ -t 1-%s

python -m examples.%s %s configfiles/baseline.cfg
'''

options = experiments.xover#experiments.ops#representation

parser = argparse.ArgumentParser()
parser.add_argument("name", help="a list of base name for the experiment sets, separated by ':'.")
parser.add_argument("module", help="the list of modules of the problems to be run, separated by ':'.")
parser.add_argument("numruns", help="the number of runs for each experiment")
args = parser.parse_args()

nametokens = args.name.split(':')
moduletokens = args.module.split(':')
assert(len(nametokens)==len(moduletokens))

for name, module in zip(nametokens, moduletokens):
    #f = open('run-'+args.name,'w')
    #f.write(body%(args.name, args.numruns,args.module, ''))
    #f.close()

    for k,optlist in options.items():
        f = open('run-'+name + '-' + k, 'w')
        optionstr=""
        for pair in optlist:
            optionstr += pair + " "
        f.write(body%(name+"-"+k, args.numruns, module, optionstr))
        f.close()
