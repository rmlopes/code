#!/usr/bin/python
# Filename: group-data.py

import sys
import os
import argparse
import string
#FIXME import from code
evologhead = str('Best\tNumProteins\tNumFunctions\tNumInactive\t'+
                 'GeneSize\tEvaluations\tAvgBest\tAvgNumProteins\t'+
                 'AvgNumFunctions\tAvgGeneSize\tValidation\t')#+
                 #'NeutralMut\tNeutralBits')

class FilenameParsers:
    class representation:
        headers = "Problem\tExperiment\tRunIdx\t" + evologhead + "\n"
        def parsefilename(self,filename):
            class fn: pass
            f = fn()
            tokens = filename.split('_')
            try:
                f.problem =  tokens[0].split('-',1)[0]
                f.label = tokens[0].split('-',1)[1]
                #label = string.replace(label,"-nooverlap","F")
            except:
                f.label = 'dm5'
            f.idx = tokens[-1]
            return f

        def tostr(self,s):
            return s.problem+"\t"+s.label+"\t"+s.idx

    class operators:
        headers = "Problem\tOp\tOpSize\tRates\tRunIdx\t" + evologhead + "\n"
        def parsefilename(self,filename):
            class fn: pass
            f = fn()
            tokens = filename.split('_')
            f.problem =  tokens[0].split('-',1)[0]
            f.ops = tokens[0].split('-',1)[1]
            if f.ops == "GCD":
                f.size = "len(gene)"
                f.rates = tokens[-3]
            elif f.ops == 'rnd7F':
                f.size = '0'
                f.rates = '0'
            else:
                f.size = tokens[1]
                f.rates = tokens[-3]
            f.idx = tokens[-1]
            return f

        def tostr(self,s):
            return s.problem+"\t"+s.ops+"\t"+s.size +"\t"+s.rates+"\t"+s.idx

    class xover:
        headers = "Problem\tOp\tRate\tRunIdx\t" + evologhead + "\n"
        def parsefilename(self,filename):
            class fn: pass
            f = fn()
            tokens = filename.split('_')
            f.problem =  tokens[0].split('-',1)[0]
            f.op = tokens[0].split('-')[1]

            if f.op == "rnd7F":
                f.rate = '0.0'
            else:
                f.rate = tokens[0].split('-')[2]

            f.idx = tokens[-1].split('.')[0]
            return f

        def tostr(self,s):
            return s.problem+"\t"+s.op+"\t"+s.rate+"\t"+s.idx

    class feedback:
        headers = "Problem\tRunIdx\t" + evologhead + "\n"
        def parsefilename(self,filename):
            class fn: pass
            f = fn()
            tokens = filename.split('_')
            f.problem =  tokens[0]#.split('-',1)[0]
            #f.op = tokens[0].split('-')[1]

            #if f.op == "rnd7F":
            #    f.rate = '0.0'
            #else:
             #   f.rate = tokens[0].split('-')[2]

            f.idx = tokens[-1]
            return f

        def tostr(self,s):
            return s.problem+"\t"+s.idx

parser = argparse.ArgumentParser()
parser.add_argument("dir",
                    help="the directory containing the evolutionary logs.")
parser.add_argument("-p","--pattern", default="_evolution_",
                    help="the pattern used to distinguish the \
                    evolutionary run log files.")
parser.add_argument("-fp", "--fnameparser", default="representation",
                    help="the parser to read the filename and separate the experiments.")
parser.add_argument("-g","--validation", default="_validation_",
                    help="the pattern used to distinguish the validation result files.")
args = parser.parse_args()

files = filter(lambda x: x.find(args.pattern) != -1,
               os.listdir(args.dir))
newlines=[]

absfnparse = getattr(FilenameParsers,args.fnameparser)
fnparser = absfnparse()
sys.stdout.write(absfnparse.headers)
for f in files:
    parsed = fnparser.parsefilename(f)
    fh = open(args.dir+f,'r')
    lines = fh.readlines()
    valfh = open(args.dir+ f.replace("evolution","validation")+'.save','r')
    vlines = valfh.readlines()
    if vlines:
        v = vlines[-1]
    else:
        v = '-1\n'
    towrite = fnparser.tostr(parsed) +"\t"+lines[-1][:-1]+'\t'+v
    sys.stdout.write(towrite)

'''
Search in error and output stream files.

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
'''
    #to filter failed experiments
    #if lines[-1].find('Best') == -1:
