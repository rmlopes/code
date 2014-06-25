#!/usr/bin/python
import sys
import os
import cPickle as pickle
import pydot
from pydot import graph_from_dot_data
import re
import argparse

_NLABEL = ' -  %s  - \n %s '

parser = argparse.ArgumentParser()
parser.add_argument("file",
                    help="the ancestors file.")
parser.add_argument("-png",type=bool,
                    help="use to create the frame images  for the anim.")
parser.add_argument("-gif",type=bool,
                    help="use to create the gif from the frame images.")

def split_nodes_edges(gstr):
    nodes = []
    edges = []
    tokens = gstr.split('\n')
    for line in tokens[2:-1]:
        if '->' in line:
            edges.append(line)
        else:
            nodes.append(line)
    return (nodes,edges)

def getclusters(nodes):
    clusters = []
    nodes = map(lambda x:int(x),nodes)
    nodes.sort()
    #print nodes
    for n in nodes:
        open_ = filter(lambda x: x not in [n
                                           for c in clusters
                                           for n in c],
                       nodes)
        cluster = filter(lambda x: abs(int(x) - int(n)) < 256,
                         open_)
        if len(cluster) > 1:
            clusters.append(cluster)
        elif len(cluster) == 1:
            for c in clusters:
                if any(map(lambda x: abs(x - cluster[0]) < 256, c)):
                    c.append(cluster[0])
    return clusters

def setcolor(c, g):
    g.obj_dict['attributes']['color'] = c
    g.obj_dict['attributes']['fontcolor'] = c
    g.obj_dict['attributes']['labelfontcolor'] = c
    g.get_node('node')[0].obj_dict['attributes']['color'] = c
    g.get_node('node')[0].obj_dict['attributes']['fontcolor'] = c
    g.get_node('edge')[0].obj_dict['attributes']['color'] = c

def convert(fname):
    cmd = 'convert -density 300x300 %s.ps %s.png' % (fname,fname)
    os.system(cmd)
    os.system('rm %s*.ps' % fname)


############ START ##############
args = parser.parse_args()

treefname = args.file
tfile = open(treefname,'r')
lines= tfile.readlines()
filename = treefname.split('.')[-2].split('/')[-1]


#load ancestors (old to new)
graphs = []
for i in range(0,len(lines),2):
    dotstr = pickle.loads(lines[i]+lines[i+1])
    graphs.append(dotstr)

graphs = graphs[::-1]
graphs = map(split_nodes_edges, graphs)

#build layered graph
layers = map(lambda x: 'l%02d' % x, range(len(graphs)))

nodesglobal = dict()
edgesglobal = dict()

#keys = layers
outputmap = dict()
inputmap = dict()
labelmap = dict()

#to read the nodes
labelp = re.compile(r'\"(.+?)\"')
nump = re.compile(r'[0-9]*')

for g,l in zip(graphs, layers):
    lstr = '[layer="' + l
    #lg += 'node [layer="%s"];\n' % l
    #print g[0]
    outputmap[l] = int(g[0][0].split(' ')[0])
    inputmap[l] = []
    labelmap[l] = dict()
    for n in g[0]:
        idx = nump.search(n).group()
        label = labelp.findall(n)
        labelmap[l][idx] = label[0]
        if 'M' in n or 'L' in n or 'R' in n or 'input' in n:
            inputmap[l].append(idx)
        if idx not in nodesglobal.keys():
            nodesglobal[idx] = lstr
        else:
            nodesglobal[idx] = "%s%s%s" % (nodesglobal[idx],',',l)
    #print g[0]
    #print labelmap[l]

lg = 'digraph layered {\ngraph[rotate=90,size="10,10"];\nordering = out;\nsplines="ortho"\nlayers="%s"\n' % reduce(lambda x,y: x+':'+y, layers)
#strnodes = map(lambda x: '%s %s"];\n' % (x[0],x[1]), nodesglobal.items())
strnodes = map(lambda x: '%s;\n' % x, nodesglobal.keys())
lg += 'node[layer="all",shape=box];\n'
lg += reduce(lambda x,y: x + y , strnodes)

clusters = getclusters(nodesglobal.keys())
for cluster,idx in zip(clusters,range(len(clusters))):
    cstr = 'subgraph cluster%i{\ngraph[style=dotted];\nedge[style=invis];\n' % idx
    rankstr = ''
    reducerule = '\n'
    if cluster[-1]-cluster[0] < 256 and len(cluster)>2:
        rankstr = 'rank=same\n'
        reducerule = ' -> '
    cstr += rankstr
    cstr += reduce(lambda x,y: str(x) + reducerule + str(y), cluster)
    lg += cstr + '\n}\n'


lg += 'edge[style=solid];\n'
for g,l in zip(graphs, layers):
    lstr = 'layer="' + l
    #lg += 'node [layer="%s"];\n' % l
    for e in g[1]:
        if e not in edgesglobal.keys():
            edgesglobal[e] = lstr
        else:
            edgesglobal[e] = "%s%s%s" % (edgesglobal[e],',',l)

stredges = map(lambda x: '%s,%s"];\n' % (x[0][:-2],x[1]), edgesglobal.items())
lg += reduce(lambda x,y: x + y ,stredges)
lg += '}'

#print lg
g = graph_from_dot_data(lg)
g.write('%s_anim.dot' % filename)


#build animation
animcmd = 'convert '
#set initial labels to obtain correct layout
for n in filter(lambda n: n.get_name() not in ['graph','edge'], g.get_nodes()):
    n.obj_dict['attributes']['label'] = _NLABEL % (n.get_name(),'_')
#layers = split(':',g.get_layers())
#stylep = re.compile(r'(,*?invis?)?')

for l in layers:
    print '### %s ###' % l
    loutname = "%s_%s" % (filename, l)

    if args.png:
        g.obj_dict['attributes']['layerselect'] = l
        g.get_node(str(outputmap[l]))[0].obj_dict['attributes']['style'] = 'diagonals'
        for idx in inputmap[l]:
            g.get_node(str(idx))[0].obj_dict['attributes']['style'] = 'dashed'
        for n in g.get_nodes():
            if n.get_name() in labelmap[l].keys():#[gr[0] for gr in]n.obj_dict['attributes']['layer']:
                n.obj_dict['attributes']['label'] = _NLABEL  % (n.get_name(), labelmap[l][n.get_name()])
            else:
                if n.get_name()not in ['graph','node','edge']:
                    n.obj_dict['attributes']['style'] = 'invis'

        setcolor('black',g)
        g.write_ps("%s_out.ps" % loutname)
        setcolor('grey',g)
        g.write_ps('%s_mask.ps' % loutname)

        for n in g.get_nodes():
            if n.get_name() not in ['graph','edge']:
                n.obj_dict['attributes']['style'] = "solid"

        convert('%s_out' % loutname)
        convert('%s_mask' %  loutname)

    animcmd += '-dispose none -delay 0 %s_mask.png ' % loutname
    animcmd += '-dispose previous -delay 300 %s_out.png ' % loutname

animcmd += '-loop 1 -rotate 90  %s.gif' % loutname
print animcmd

if args.gif:
    os.system(animcmd)
exit(0)
