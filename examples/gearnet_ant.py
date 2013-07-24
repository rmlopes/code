from code.operators import *
from code.evodevo import EvoDevoWorkbench
from code.gearnet import *
from gp_ant import *
from code.utils.config import *
import code.grammars.grammardb as gdb

def evaluateant(agent, **kwargs):
        #agent.funskel = mergefun
        #if len(agent) < 3:
         #       return 90
        #routine = agent()
        #print routine
        #if len(routine) < 1000:
        routine = agent.circuit
        #print routine
        if agent.isvalid() and len(routine)<1000:
           ant.runstring(routine,True)  
        else:
           return 1e6
        #a = agent()
        #print a
        #
        #else:
        #        return 100
        return 89 - ant.eaten


if __name__ == '__main__':
        p = Prob(evaluateant, gdb._artificial_ant, 'code')
        cfg = loadconfig(parsecmd())
        edw = EvoDevoWorkbench(cfg,p)
        edw.run()
