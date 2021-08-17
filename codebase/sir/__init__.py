from codebase import get_logger
from datetime import date, timedelta
from os import environ, path
from tempfile import gettempdir
import json
import numpy as np
import matplotlib.pyplot as plt
import urllib.request
from scipy import stats
import codebase.sir.models
import codebase.sir.fetch
import codebase.sir.display
import codebase.sir.estimation

logger = get_logger()
vals = models.params()

def run_comp_and_plot(state,county=None,
                     n_days=30,
                     n_ramp=3,
                     Tea=0.5,
                     Tiq=0.5,
                     Tai=5.0,
                     Tir=14.0,
                     Tid=20.0,
                     Tih=1.0,
                     piq=0.9,
                     pai=0.6,
                     rai=0.4,
                     Sd=0.6,
                     Sd_period=0,
                     Sd_delay=12,
                     detection_rate=0.2,
                     data_days=30,
                     n_substeps=10,
                     refit=False):

    logger.info('Running run_comp_and_plot')
    
    model = models.SEAIQHRD(state,county,
                            n_days,n_substeps,n_ramp,
                            Tea,Tai,Tir,Tiq,Tid,
                            Tih,pai,piq,rai,Sd,
                            Sd_delay,detection_rate,
                            data_days,refit)
    '''

    model = models.SIR(state,county,
                       n_days,n_substeps,n_ramp,
                       Tir,Sd,Sd_delay,detection_rate,
                       data_days,refit)
    
    '''
    
    #get data and initial values
    model.initialize()    
    model.define_parameters()
    
    #get model parameters
    model.update_parameters()
    print("beta,Sd,Sd_delay,n_ramp: %s,%s,%s,%s" %(model.params.beta,model.params.Sd,model.params.Sd_delay,model.params.n_ramp))

    #run model
    model.simulate()

    #error on data
    err=model.calculate_error()
    print("Data difference: %s" %err)
    
    #plot results
    #comp_plot_path = display.plot_comparison_SIR(model)
    comp_plot_path = display.plot_comparison_SEAIQHRD(model)
    
    return comp_plot_path
