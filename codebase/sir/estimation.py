import numpy as np
import math
from datetime import date,timedelta,datetime
from scipy import stats
from scipy.optimize import curve_fit
from codebase.sir import fetch


def exp_func(t,a,b):
    return a*np.exp(b*t)

def beta_interp(t,beta,Sd,t0,n_ramp):
    if t<t0:
        return beta
    elif t0<=t<=t0+n_ramp:
        return (1-(t-t0)*Sd/(n_ramp))*beta
    else:
        return (1-Sd)*beta
    
def exp_decay_func(t,a,beta,Sd,t0,n_ramp):
    
    return [a*np.exp(i*beta_interp(i,beta,Sd,t0,n_ramp)) for i in t]

def exp_lin_func(t,a,b,c):
    return a*np.exp(b*t)+c*t

def logistic_func(t,a,b,c,d,e):
    return a/(d*np.exp(-b*(t-e))+c)

class initial_values:
    pass

def cases_to_beta(data,params,init_vals):
    
    tmp=data['I_ref']
    e_idx=params.Sd_delay 
    popt,pconv=curve_fit(exp_func,[i for i in range(len(tmp[:e_idx]))],tmp[:e_idx])
    
    params.R0=(popt[1]/params.gamma+1)
    params.beta=params.R0*params.gamma

    return params

def cases_to_intercepts(data,params):
    vals=initial_values()
    e_idx=params.Sd_delay
    
    vals={}
    days=[i for i in range(len(data['dates'][:e_idx]))]
    for k,v in data.items():
        if k in ['S','E','A','I','Q','H','R','D']:
            try:
                #popt,pconv=curve_fit(exp_func,days,
                #                     v[:e_idx])
                #intercept=popt[0]

            #except:
                r=stats.linregress(days,
                [np.log(x+1) for x in v[:e_idx]])
                intercept=np.exp(r.intercept)
            except:
                pass
            
            vals[k]=intercept
    return vals

def doubling_trend_to_Sd(Td_array):
    days=[i for i in range(len(Td_array))]
    r=stats.linregress(days,Td_array)
    Sd=r.slope/(1+r.slope)
    return Sd 

