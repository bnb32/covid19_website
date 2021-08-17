import sys
sys.path.append('../../')
from codebase.sir import run_comp_and_plot
from codebase.sir import fetch

run_comp_and_plot("New York",n_days=1,detection_rate=0.1,piq=0.9,pai=0.6,rai=0.65,refit=True)
