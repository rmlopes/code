#!/usr/bin/python
# Filename: group-data.py

import sys
import os
import argparse
import string
#FIXME import from code
evologhead = str('Best\tNumProteins\tNumFunctions\tNumInactive\t'+
                 'GeneSize\tEvaluations\tAvgBest\tAvgNumProteins\t'+
                 'AvgNumFunctions\tAvgGeneSize\tGeneralization\t'+
                 'NeutralMut\tNeutralBits')

parser = argparse.ArgumentParser()
parser.add_argument("dir",
                    help="the directory containing the evolutionary logs.")
parser.add_argument("-p","--pattern", default="_evolution_",
                    help="the pattern used to distinguish the \
                    evolutionary run log files.")
parser.add_argument("-g","--include_generalization", default=True,
                    help="the module of the problem to be run.")
args = parser.parse_args()

files = filter(lambda x: x.find(args.pattern) != -1,
               os.listdir(args.dir))
newlines=[]
headers = "Experiment\tRunIdx\t" + evologhead + "\n"
sys.stdout.write(headers)
for f in files:
    tokens = f.split('_')
    try:
        label = tokens[0].split('-',1)[1]
        label = string.replace(label,"-nooverlap","F")
    except:
        label = 'dm5'
    idx = tokens[2]
    fh = open(args.dir+f,'r')
    lines = fh.readlines()
    towrite = label+"\t"+idx+"\t"+lines[-1][:-1]

    if args.include_generalization:
        lenidx = len(tokens[-1])
        genfile = filter(lambda x: x.find(tokens[0]) != -1 and \
                         x[-lenidx-1:].find("."+tokens[-1]) != -1,
                         os.listdir(args.dir))
        genfh = open(args.dir+genfile[0],'r')
        genlines = genfh.readlines()
        if genlines:
            genresult = genlines[-1][:-1] + "\t"
            neutrals = filter(lambda x: x.find("Avg. Neutral") != -1,
                              genlines)
            if len(neutrals) > 1:
                neutralmuts = neutrals[0].split(" ")[-1][:-1]
                neutralbits = neutrals[1].split(" ")[-1]
            else:
                neutralmuts = '0.0'
                neutralbits = neutrals[0].split(" ")[-1]

            genresult += "\t" + neutralmuts +"\t"
            genresult += neutralbits

        else:
            genresult = '1e6\n'

        towrite += "\t" + genresult

    #to filter failed experiments
    if lines[-1].find('Best') == -1:
            sys.stdout.write(towrite)
