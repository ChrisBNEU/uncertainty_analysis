import os
import sys
sys.path.append("/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/")
from rmg_gua.gua_cantera.Spinning_basket_reactor.sbrTest import TestSBR


tsbr = TestSBR()
tsbr.test_change_reactions()
tsbr.test_change_be()
tsbr.teatdown()


