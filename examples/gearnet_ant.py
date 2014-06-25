from code.operators import *
from code.evodevo import EvoDevoWorkbench
from code.gearnet import *
from gp_ant import *
from code.utils.config import *
import code.grammars.grammardb as gdb

#from galapagos
def python_filter(txt):
    """ Create correct python syntax.
    We use {: and :} as special open and close brackets, because
    it's not possible to specify indentation correctly in a BNF
    grammar without this type of scheme."""
    indent_level = 0
    tmp = txt[:]
    i = 0
    while i < len(tmp):
      tok = tmp[i:i+2]
      if tok == "{:":
        indent_level += 1
      elif tok == ":}":
        indent_level -= 1
      tabstr = "\n" + "  " * indent_level
      if tok == "{:" or tok == ":}" or tok == "\\n":
        tmp = tmp.replace(tok, tabstr, 1)
      i += 1
    # Strip superfluous blank lines.
    txt = "\n".join([line for line in tmp.split("\n") if line.strip() != ""])
    return txt

def evaluateant(agent, **kwargs):
        #agent.funskel = mergefun
        #if len(agent) < 3:
         #       return 90
        #routine = agent()
        #print routine
        #if len(routine) < 1000:
        agent.circuit = python_filter(agent.circuit)
        routine = agent.circuit
        #print routine
        if agent.isvalid():
           ant.runstring(routine)
        else:
           return 1e6
        #a = agent()
        #print a
        #
        #else:
        #        return 100
        return 89 - ant.eaten


if __name__ == '__main__':
        p = Prob(evaluateant, gdb._ant_grammar, 'code')
        cfg = loadconfig(parsecmd())
        edw = EvoDevoWorkbench(cfg,p)
        edw.run()
        edw.runlog.validatelog.critical(edw.best.moves)
