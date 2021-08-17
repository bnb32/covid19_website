from codebase import get_logger
from datetime import date, timedelta
from os import environ, path
from tempfile import gettempdir
import numpy as np
import matplotlib.pyplot as plt
from codebase.sir import estimation

logger=get_logger()

def plot_doubling_trend(dates,data,state):
    days=[i for i in range(len(dates))]
    plt.plot(days,data)
    tick_step=max((1,int(len(days)/3)-1))
    plt.xticks(days[0::tick_step],dates[0::tick_step],rotation=30)

    plt.title(state)
    plt.grid(True)
    plt.ylabel("Doubling Time (days)")
    plt.xlabel("Date")
    plt.tight_layout()
    plt.savefig('Td_trend.png')
    plt.clf()


def plot_comparison_SIR(model):
    params=model.params
    init_vals=model.init_vals
    data=model.data
    sim=model.results
    plot_name = path.join(gettempdir(), 'comp_plot.png')
    dates=data['dates']
    data_length=len(dates)

    dates=[str((dates[0]+timedelta(days=i)).strftime('%m-%d-%y')) for i in range(len(sim['I']))]
    dat_days=[i for i in range(len(data['I_ref']))]
    sim_days=[i for i in range(len(sim['I']))]
    
    
    plt.scatter(dat_days,data['I_ref'],color=(0,0,1),marker='x',label='Infections-Data')
    
    plt.plot(sim_days,sim['I'],color=(0,0,1),label='Infections-Sim')
 
    tick_step=max((1,int(len(sim_days)/10)-1))
    plt.xticks(sim_days[0::tick_step],dates[0::tick_step],rotation=90)
    
    if params.county==None:
        plt.title(params.state+' (R0=%s, Rf=%s)'%(float('%.3g' %(params.R0)),float('%.3g' %(params.Rf))))
    else:    
        plt.title(params.state+'-'+params.county+' (R0=%s, Rf=%s)'%(float('%.3g' %(params.R0)),float('%.3g' %(params.Rf))))
    
    ax=plt.gca()
    ax.axvline(x=data_length-1,color='k')    
    plt.grid(True)
    
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(plot_name)
    logger.info("Saved plot as {0}".format(plot_name))
    
    plt.clf()

    return plot_name

def plot_comparison_SEAIQHRD(model):
    params=model.params
    init_vals=model.init_vals
    data=model.data
    sim=model.results
    plot_name = path.join(gettempdir(), 'comp_plot.png')
    dates=data['dates']
    data_length=len(dates)

    dates=[str((dates[0]+timedelta(days=i)).strftime('%m-%d-%y')) for i in range(len(sim['I']))]
    dat_days=[i for i in range(len(data['I_ref']))]
    sim_days=[i for i in range(len(sim['I']))]
    
    #plt.scatter(dat_days,data['I']+data['Q']+data['A'],color=(0,0,1),marker='x', label='Infected-Data')
    #plt.plot(sim_days,sim['I']+sim['Q']+sim['A'],color=(0,0,1),label='Infected-Sim')
    
    plt.scatter(dat_days,data['I']+data['Q'],color=(0,1,0),marker='x', label='I-Data')
    plt.plot(sim_days,sim['I']+sim['Q'],color=(0,1,0),label='I-Sim')
    
    plt.scatter(dat_days,data['A'],color=(0,0,1),marker='x', label='A-Data')
    
    plt.plot(sim_days,sim['A'],color=(0,0,1),label='A-Sim')

    plt.scatter(dat_days,data['H'],color=(1,0,1),marker='x',label='H-Data')
    plt.plot(sim_days,sim['H'],color=(1,0,1),label='H-Sim')
    
    plt.scatter(dat_days,data['D'],color=(1,0,0),marker='x',label='D-Data')
    
    plt.plot(sim_days,sim['D'],color=(1,0,0),label='D-Sim')
    tick_step=max((1,int(len(sim_days)/10)-1))
    plt.xticks(sim_days[0::tick_step],dates[0::tick_step],rotation=90)
    
    if params.county==None:
        plt.title(params.state+' (R0=%s, Rf=%s)'%(float('%.3g' %(params.R0)),float('%.3g' %(params.Rf))))
    else:    
        plt.title(params.state+'-'+params.county+' (R0=%s, Rf=%s)'%(float('%.3g' %(params.R0)),float('%.3g' %(params.Rf))))
    
    ax=plt.gca()
    ax.axvline(x=data_length-1,color='k')    
    plt.grid(True)
    
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(plot_name)
    logger.info("Saved plot as {0}".format(plot_name))
    
    plt.clf()

    return plot_name
