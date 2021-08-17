import numpy as np
import math
from datetime import date,timedelta,datetime
from scipy import stats
from codebase import get_logger
from codebase.sir import estimation
from codebase.sir import fetch
import itertools
from os import path, system
from tempfile import gettempdir
import json
from multiprocessing import Pool
from scipy.optimize import curve_fit

logger = get_logger()
parallel_execute=False
#parallel_execute=True


class params:
    pass

class initial_values:
    pass

class SIR:
    def __init__(self,state,county,
                 n_days,n_substeps,
                 n_ramp,Tir,
                 Sd,Sd_delay,
                 detection_rate,
                 data_days,refit):
        
        self.params=params()
        self.params.refit=refit
        self.params.state=state
        self.params.county=county
        self.params.n_days=n_days
        self.params.n_substeps=n_substeps
        self.params.dt=1.0/n_substeps
        self.params.Tir=Tir
        self.params.Sd=Sd
        self.params.Sd_delay=Sd_delay
        self.params.n_ramp=n_ramp
        self.params.detection_rate=detection_rate
        self.params.data_days=data_days
        self.params.search_steps=30
        self.params.descent_steps=500
        self.params.eps=1e-8

        self.init_vals=initial_values()
        self.init_vals.N=fetch.location_to_population(state,county)
        order_date=fetch.order_dates[state]
        if order_date: 
            self.params.Sd_date=datetime.strptime(order_date,"%Y-%m-%d").date()
        else:
            self.params.Sd_date=datetime.today().date()

        
        self.results={}
        self.data={}

    def initialize(self):
        self.data=fetch.location_to_data(self.params)
        self.data['I_ref']=self.data['I_act']
        self.params.data_days=len(self.data['dates'])
        self.params.Sd_delay=min(((self.params.Sd_date-self.data['dates'][0]).days,self.params.data_days))

        self.params.n_days=self.params.data_days+self.params.n_days
        
        vals=estimation.cases_to_intercepts(self.data,self.params)   
        v=self.init_vals
        v.I=vals['I']
        v.R=vals['D']+vals['R']
        v.S=v.N-v.I-v.R
        self.init_vals=v
        return self.init_vals

    def update_parameters(self):
       
        p=self.params
        
        self.params.i_decay=1.0/self.params.Tir
        self.params.r_growth=self.params.i_decay
        self.params.gamma=p.i_decay

        self.params=estimation.cases_to_beta(self.data,self.params,self.init_vals)
        self.param_fit()
        #self.grad_descent_param_fit()
        self.params.R0=self.params.beta/self.params.gamma
        
        self.params.Rf=self.params.R0*(1-self.params.Sd)
        self.params.beta_array=get_beta_ramp(self.params)
        
        return self.params

    def rk_step(self,v,N,beta):
        p=self.params
        dt=p.dt
        rk_step=[0.5,0.5,1,0]
        Sk=[];Ik=[];Rk=[];
    
        s,i,r=v.S,v.I,v.R
        for t in range(4):
    
            sk=-beta*s*(i)/N*dt
            ik=(beta*s*(i)/N-i*p.i_decay)*dt
            rk=(i*p.r_growth)*dt
            
            s=sk*rk_step[t]+v.S
            i=ik*rk_step[t]+v.I
            r=rk*rk_step[t]+v.R
            
            Sk.append(sk)
            Ik.append(ik)
            Rk.append(rk)
    
        vn=initial_values()
        stencil=[1,2,2,1]
        vn.S=v.S+1.0/6.0*np.dot(Sk,stencil)
        vn.I=v.I+1.0/6.0*np.dot(Ik,stencil)
        vn.R=v.R+1.0/6.0*np.dot(Rk,stencil)
    
        if vn.S<0.0: vn.S=0.0
        if vn.I<0.0: vn.I=0.0
        if vn.R<0.0: vn.R=0.0
    
        return vn

    def param_iteration(self,params,param_combos,i):
        p=param_combos[i]
        Sd=p[0]
        Sd_delay=p[1]
        n_ramp=p[2]

        beta=params.beta
        params.Sd=Sd
        params.Sd_delay=Sd_delay
        params.n_ramp=n_ramp
        beta_array=get_beta_ramp(params)
        self.run_model(params.data_days,
                       beta_array)
        err=self.calculate_error()
        #print("Sd,Sd_delay,n_ramp,err: %s,%s,%s,%s"%(Sd,Sd_delay,n_ramp,err))
        return err
        #error_queue.put({i:err})

    def param_fit(self):
        return self.grid_search_param_fit()
    
    def grid_search_param_fit(self):
        params=self.params
        v=self.init_vals
         
        Sd_delay_min=max((params.Sd_delay-2,0))
        Sd_delay_max=min((params.Sd_delay+2,params.data_days))

        Sd_list=np.linspace(0.1,0.9,params.search_steps)
        Sd_delay_list=np.arange(Sd_delay_min,Sd_delay_max)
        n_ramp_list=np.arange(1,25)
        param_combos=list(itertools.product(*[Sd_list,Sd_delay_list,n_ramp_list]))
        
        errs=[0]*len(param_combos)
        for i in range(len(errs)):
            #job=Process(target=self.param_iteration, args=(params,param_combos,i))
            #job.start()
            #errs[i]=error_queue.get()[i]
            #job.join()
            errs[i]=self.param_iteration(params,param_combos,i)
        
        idx=errs.index(min(errs))
        self.params.Sd=param_combos[idx][0]
        self.params.Sd_delay=param_combos[idx][1]
        self.params.n_ramp=param_combos[idx][2]
        return self.params

    def grad_descent_param_fit(self):
        params=self.params
        v=self.init_vals

        p0={
            'Sd':0.5,
            'Sd_delay':params.Sd_delay,
            'n_ramp':10}
        p1=p0.copy()
        alpha={
               'Sd':1e-6,
               'Sd_delay':1e-8,
               'n_ramp':1e-8}
        
        dp={}; J0={}; J1={};
        for k in p0:
            dp[k]=1e-3
            J0[k]=1
            J1[k]=0

        cost0=1
        cost1=0
        p_increase=0.05
        p_decrease=0.1
        step=1

        while np.abs(cost1-cost0)>params.eps:   
           cost_diff=cost1-cost0
           #print("cost diff: %s"%(cost_diff))
           #print("cost: %s"%(cost0))
           
           if cost_diff<0:
               for k in alpha:
                   alpha[k]*=(1+p_increase)
               p0=p1
           else:
               for k in alpha:    
                   alpha[k]*=p_decrease
           
           params.Sd=p0['Sd']
           params.Sd_delay=p0['Sd_delay']
           params.n_ramp=p0['n_ramp']
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           
           cost0=self.calculate_error()

           for k in p0:
               J0[k]=cost0
               tmp=p0.copy()
               tmp[k]+=dp[k]
               params.Sd=tmp['Sd']
               params.Sd_delay=tmp['Sd_delay']
               params.n_ramp=tmp['n_ramp']
               beta_array=get_beta_ramp(params)
               self.run_model(params.data_days,
                              beta_array)
               J1[k]=self.calculate_error()
           
           p1=grad_descent(p0,J0,J1,dp,alpha)
           params.Sd=p1['Sd']
           params.Sd_delay=p1['Sd_delay']
           params.n_ramp=p1['n_ramp']
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           cost1=self.calculate_error()
           step+=1

        self.params.Sd=p0['Sd']
        self.params.Sd_delay=p0['Sd_delay']
        self.params.n_ramp=p0['n_ramp']
        return self.params

    def simulate(self):
        return self.run_model(self.params.n_days,
                              self.params.beta_array)

    def run_model(self,n_days,beta_array):
        params=self.params
        v=self.init_vals
        s,i,r=[v.S],[v.I],[v.R]
        N=v.S+v.I+v.R
        for day in range(n_days):
            beta=beta_array[day]
            for t in range(params.n_substeps):
                v=self.rk_step(v,N,beta)
            s.append(v.S)
            i.append(v.I)
            r.append(v.R)
        self.results={'S':np.array(s),
                      'I':np.array(i),
                      'R':np.array(r),
                     }
        return self.results

    def calculate_error(self):
        errs=[]
        #errs.append(cost_func(self.data['I_ref'],self.results['I']))
        errs.append(cost_func(vec_diff(self.data['I_ref']),vec_diff(self.results['I'])))
        '''
        for k in self.data:
            if k in self.results:
                errs.append(cost_func(self.data[k],
                                      self.results[k]))
        '''
        return np.sum(errs) 

class SEAIQHRD:
    def __init__(self,state,county,
                 n_days,n_substeps,
                 n_ramp,
                 Tea,Tai,Tir,
                 Tiq,Tid,
                 Tih,pai,piq,rai,
                 Sd,Sd_delay,
                 detection_rate,
                 data_days,refit):
        
        self.params=params()
        self.params.refit=refit
        self.params.beta_decay=0.0
        self.params.state=state
        self.params.county=county
        self.params.n_days=n_days
        self.params.n_substeps=n_substeps
        self.params.dt=1.0/n_substeps
        self.params.Tea=Tea
        self.params.Tai=Tai
        self.params.Tir=Tir
        self.params.Tiq=Tiq
        self.params.Tid=Tid
        self.params.Tih=Tih
        self.params.rai=rai
        self.params.Sd=Sd
        self.params.Sd_delay=Sd_delay
        self.params.n_ramp=n_ramp
        self.params.detection_rate=detection_rate
        self.params.data_days=data_days
        self.params.piq=piq
        self.params.pai=pai
        self.params.rai=1.0/self.params.pai-1
        self.params.search_steps=15
        self.params.descent_steps=500
        self.params.eps=1e-8

        self.init_vals=initial_values()
        self.init_vals.N=fetch.location_to_population(state,county)
        order_date=fetch.order_dates[state]
        if order_date: 
            self.params.Sd_date=datetime.strptime(order_date,"%Y-%m-%d").date()
        else:
            self.params.Sd_date=datetime.today().date()
        
        self.results={}
        self.data={}

    def initialize(self):
        self.data=fetch.location_to_data(self.params)
        self.rates=fetch.data_to_rates(self.data,self.params)

                
    def update_data(self):
        self.data=fetch.refine_data(self.data,self.params)
        self.data['I_ref']=self.data['I']+self.data['A']+self.data['Q']
        self.params.data_days=len(self.data['dates'])
        self.params.Sd_delay=min(((self.params.Sd_date-self.data['dates'][0]).days,self.params.data_days))
        self.params.n_days=self.params.data_days+self.params.n_days
        
        vals=estimation.cases_to_intercepts(self.data,self.params)
        v=self.init_vals
        v.E=0
        v.A=vals['A']
        v.I=vals['I']
        v.D=vals['D']
        v.R=vals['R']
        v.H=vals['H']
        v.Q=vals['Q']
        v.S=v.N-v.A-v.I-v.Q-v.E-v.R-v.D-v.H
        self.init_vals=v

    
    def define_parameters(self):

        p=self.params
        dt=p.dt
        
        p.pid=self.rates.pid
        p.pqd=self.rates.pqd
        p.phd=self.rates.phd
        p.pih=self.rates.pih
        p.pqh=self.rates.pqh
        p.Tar=p.Tai+p.Tir
        p.Tqh=max((p.Tih-p.Tiq,dt))
        p.Th=(p.Tih+p.Tqh)/2.0
        p.Tqr=max((p.Tir-p.Tiq,dt))
        p.Tqd=max((p.Tid-p.Tiq,dt))
        p.Thr=max((p.Tir-p.Th,dt))
        p.Thd=max((p.Tid-p.Th,dt))
        p.pir=(1-p.piq-p.pid-p.pih)
        p.pqr=(1-p.pqd-p.pqh)
        p.par=(1-p.pai)
        p.phr=(1-p.phd)
                
        p.e_decay=1.0/p.Tea
        
        p.a_growth=p.e_decay
        p.a_recover=p.par/p.Tar
        p.a_decay=p.pai/p.Tai+p.a_recover

        p.i_growth=p.pai/p.Tai
        p.i_recover=p.pir/p.Tir
        p.i_death=p.pid/p.Tid
        p.i_decay=p.pih/p.Tih+p.piq/p.Tiq+p.i_recover+p.i_death
        
        p.q_growth=p.piq/p.Tiq
        p.q_recover=p.pqr/p.Tqr
        p.q_death=p.pqd/p.Tqd
        p.q_decay=p.q_recover+p.q_death+p.pqh/p.Tqh
        
        p.h_growth_q=p.pqh/p.Tqh
        p.h_growth_i=p.pih/p.Tih
        p.h_recover=p.phr/p.Thr
        p.h_death=p.phd/p.Thd
        p.h_decay=p.h_death+p.h_recover
        
        p.r_growth_q=p.q_recover
        p.r_growth_i=p.i_recover
        p.r_growth_a=p.a_recover
        p.r_growth_h=p.h_recover
        
        p.d_growth_q=p.q_death
        p.d_growth_i=p.i_death
        p.d_growth_h=p.h_death
        p.gamma=(p.a_decay*p.i_decay)/(p.i_growth+p.i_decay)
        
        self.params=p
        return self.params

    def update_parameters(self):
        self.update_data()
        self.params=estimation.cases_to_beta(self.data,
                                             self.params,
                                             self.init_vals)
        if not path.exists(fetch.params_path):
            system('rm -rf %s/param_data*'%(gettempdir()))
            outfile=open(fetch.params_path,'w')
        
            self.param_fit()
            params_dict={self.params.state:
                         {'beta':self.params.beta,
                          'Sd':self.params.Sd,
                          'Sd_delay':float(self.params.Sd_delay),
                          'n_ramp':float(self.params.n_ramp)}}
            json.dump(params_dict,outfile) 
            logger.info("Saved parameter file: %s" %fetch.params_path)
            #logger.info(params_dict)

        else:
            json_file=open(fetch.params_path,'r')
            params_dict=json.load(json_file)
            if self.params.state not in params_dict or self.params.refit:
                self.param_fit()
                params_dict[self.params.state]={
                          'beta':self.params.beta,
                          'Sd':self.params.Sd,
                          'Sd_delay':float(self.params.Sd_delay),
                          'n_ramp':float(self.params.n_ramp)}
                
                outfile=open(fetch.params_path,'w')
                json.dump(params_dict,outfile) 
                logger.info("Refit: %s"%self.params.refit)
                logger.info("Saved parameter file: %s" %fetch.params_path)
                #logger.info(params_dict)



            else:
                self.params.beta=params_dict[self.params.state]['beta']
                self.params.Sd=params_dict[self.params.state]['Sd']
                self.params.Sd_delay=params_dict[self.params.state]['Sd_delay']
                self.params.n_ramp=params_dict[self.params.state]['n_ramp']

        self.params.R0=self.params.beta/self.params.gamma
        
        self.params.Rf=self.params.R0*(1-self.params.Sd)
        self.params.beta_array=get_beta_ramp(self.params)
        
        return self.params

    def rk_step(self,v,N,beta):
        p=self.params
        dt=p.dt
        rk_step=[0.5,0.5,1,0]
        Sk=[];Ek=[];Ak=[];Ik=[];Qk=[];Hk=[];Rk=[];Dk=[]
    
        s,e,a,i,q,h,r,d=v.S,v.E,v.A,v.I,v.Q,v.H,v.R,v.D
        for t in range(4):
    
            div=N-q-h-d
            sk=-beta*s*(i+a)/div*dt
            ek=(beta*s*(i+a)/div-e*p.e_decay)*dt
            ak=(e*p.a_growth-p.a_decay*a)*dt
            ik=(a*p.i_growth-p.i_decay*i)*dt
            qk=(i*p.q_growth-p.q_decay*q)*dt
            
            h_growth=q*p.h_growth_q+i*p.h_growth_i
            hk=(h_growth-p.h_decay*h)*dt
            
            r_growth=q*p.r_growth_q+i*p.r_growth_i+a*p.r_growth_a+h*p.r_growth_h
            rk=r_growth*dt
            
            d_growth=q*p.d_growth_q+i*p.d_growth_i+h*p.d_growth_h
            dk=d_growth*dt
            
            s=sk*rk_step[t]+v.S
            e=ek*rk_step[t]+v.E
            a=ak*rk_step[t]+v.A
            i=ik*rk_step[t]+v.I
            q=qk*rk_step[t]+v.Q
            h=hk*rk_step[t]+v.H
            r=rk*rk_step[t]+v.R
            d=dk*rk_step[t]+v.D
            
            Sk.append(sk)
            Ek.append(ek)
            Ak.append(ak)
            Ik.append(ik)
            Qk.append(qk)
            Hk.append(hk)
            Rk.append(rk)
            Dk.append(dk)
    
        vn=initial_values()
        stencil=[1,2,2,1]
        vn.S=v.S+1.0/6.0*np.dot(Sk,stencil)
        vn.E=v.E+1.0/6.0*np.dot(Ek,stencil)
        vn.A=v.A+1.0/6.0*np.dot(Ak,stencil)
        vn.I=v.I+1.0/6.0*np.dot(Ik,stencil)
        vn.Q=v.Q+1.0/6.0*np.dot(Qk,stencil)
        vn.H=v.H+1.0/6.0*np.dot(Hk,stencil)
        vn.R=v.R+1.0/6.0*np.dot(Rk,stencil)
        vn.D=v.D+1.0/6.0*np.dot(Dk,stencil)
    
        if vn.S<0.0: vn.S=0.0
        if vn.E<0.0: vn.E=0.0
        if vn.A<0.0: vn.A=0.0
        if vn.I<0.0: vn.I=0.0
        if vn.Q<0.0: vn.Q=0.0
        if vn.H<0.0: vn.H=0.0
        if vn.R<0.0: vn.R=0.0
        if vn.D<0.0: vn.D=0.0
    
        return vn

    def param_fit(self):
        return self.curve_fit_param_fit()

    def grad_descent_param_fit(self):
        params=self.params
        v=self.init_vals

        p0={
            'Sd':0.5,
            'Sd_delay':params.Sd_delay,
            'n_ramp':5}
        p1=p0.copy()
        alpha={
               'Sd':1e-6,
               'Sd_delay':1e-8,
               'n_ramp':1e-8}
        
        dp={}; J0={}; J1={};
        for k in p0:
            dp[k]=1e-3
            J0[k]=1
            J1[k]=0

        cost0=1
        cost1=0
        p_increase=0.05
        p_decrease=0.5
        step=1

        while np.abs(cost1-cost0)>params.eps and step<params.descent_steps:   
           cost_diff=cost1-cost0
           #print("cost diff: %s"%(cost_diff))
           #print("cost: %s"%(cost0))
          
           if cost_diff<0 and step>1:
               for k in alpha:
                   alpha[k]*=(1+p_increase)
               p0=p1
               cost0=cost1
           elif step==1:
               params.Sd=p0['Sd']
               params.Sd_delay=p0['Sd_delay']
               params.n_ramp=p0['n_ramp']
               beta_array=get_beta_ramp(params)
               self.run_model(params.data_days,
                              beta_array)
               cost0=self.calculate_error()
           else:    
               for k in alpha:    
                   alpha[k]*=p_decrease
           
           for k in p0:
               J0[k]=cost0
               tmp=p0.copy()
               tmp[k]+=dp[k]
               params.Sd=tmp['Sd']
               params.Sd_delay=tmp['Sd_delay']
               
               #why does this work?
               '''
               if k=='Sd_delay':
                   params.n_ramp=tmp['Sd_delay']
               else:    
                   params.n_ramp=tmp['n_ramp']
               #####????
               '''
               params.n_ramp=tmp['n_ramp'] 
               beta_array=get_beta_ramp(params)
               self.run_model(params.data_days,
                              beta_array)
               J1[k]=self.calculate_error()
           
           p1=grad_descent(p0,J0,J1,dp,alpha)
           
           params.Sd=p1['Sd']
           params.Sd_delay=p1['Sd_delay']
           params.n_ramp=p1['n_ramp']
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           cost1=self.calculate_error()
           step+=1

        self.params.Sd=p0['Sd']
        self.params.Sd_delay=p0['Sd_delay']
        self.params.n_ramp=p0['n_ramp']
        return self.params

    def param_iteration(self,p):
        
        Sd=p[0]
        Sd_delay=p[1]
        n_ramp=p[2]
        beta=p[3]

        self.params.beta=beta
        self.params.Sd=Sd
        self.params.Sd_delay=Sd_delay
        self.params.n_ramp=n_ramp
        beta_array=get_beta_ramp(self.params)
        self.run_model(self.params.data_days,
                       beta_array)
        
        return self.calculate_error()
        #print("beta,Sd,Sd_delay,n_ramp,err: %s,%s,%s,%s,%s"%(beta,Sd,Sd_delay,n_ramp,err))

    def grid_search_param_fit(self):
        params=self.params
        v=self.init_vals
         
        Sd_delay_min=min((20,params.Sd_delay-2))
        Sd_delay_max=min((params.Sd_delay+2,params.data_days))
        beta_list=np.linspace(0.1,1.5,params.search_steps)
        Sd_list=np.linspace(0.0,1.0,params.search_steps)
        Sd_delay_list=np.arange(Sd_delay_min,Sd_delay_max)
        n_ramp_list=np.arange(5,10)
        combos=[Sd_list,Sd_delay_list,n_ramp_list,beta_list]
        param_combos=list(itertools.product(*combos))
        
        if parallel_execute:
            pool=Pool()
            errs=pool.map(self.param_iteration,param_combos)
            pool.close()
        else:
            errs=[self.param_iteration(p) for p in param_combos]

        idx=errs.index(min(errs))
        self.params.Sd=param_combos[idx][0]
        self.params.Sd_delay=param_combos[idx][1]
        self.params.n_ramp=param_combos[idx][2]
        self.params.beta=param_combos[idx][3]
        return self.params

    def curve_fit_param_fit(self):
        days=[i for i in range(self.params.data_days)]
        dat=self.data['A']#+self.data['I']
        ddat=vec_diff(dat)
        popt,pconv=curve_fit(self.model_fit_func,
                             days,dat,method='trf',
                             bounds=([0.0,0.1,10,5],
                             [1.5,1.0,self.params.data_days,15]))
        
        self.params.beta=popt[0]
        self.params.Sd=popt[1]
        self.params.Sd_delay=popt[2]
        self.params.n_ramp=popt[3]
        return self.params

    def simulate(self):
        return self.run_model(self.params.n_days,
                       self.params.beta_array)

    def model_fit_func(self,t,beta,Sd,Sd_delay,n_ramp):
        self.params.beta=beta
        self.params.Sd=Sd
        self.params.Sd_delay=Sd_delay
        self.params.n_ramp=n_ramp
        beta_array=get_beta_ramp(self.params)
        results=self.run_model(len(t)-1,beta_array)
        return results['A']#+results['I']#+results['Q']

    def run_model(self,n_days,beta_array):
        params=self.params
        v=self.init_vals
        s,e,a,i,q,h,r,d=[v.S],[v.E],[v.A],[v.I],[v.Q],[v.H],[v.R],[v.D]
        N=v.S+v.E+v.A+v.I+v.Q+v.H+v.R+v.D
        dt=params.dt
        for day in range(n_days):
            beta_decay=(1-params.beta_decay)**day
            beta=beta_decay*beta_array[day]
            for t in range(params.n_substeps):
                v=self.rk_step(v,N,beta)
            s.append(v.S)
            e.append(v.E)
            a.append(v.A)
            i.append(v.I)
            q.append(v.Q)
            h.append(v.H)
            r.append(v.R)
            d.append(v.D)
        self.results={'S':np.array(s),
                     'E':np.array(e),
                     'A':np.array(a),
                     'I':np.array(i),
                     'Q':np.array(q),
                     'H':np.array(h),
                     'R':np.array(r),
                     'D':np.array(d)
                     }
        return self.results

    def calculate_error(self):
        errs=[]
        dat=self.data['I']+self.data['A']
        sim=self.results['I']+self.results['A']
        errs.append(cost_func(dat,sim))
        errs.append(cost_func(vec_diff(dat),vec_diff(sim)))
        '''
        for k in self.data:
            if k in self.results:
                errs.append(cost_func(self.data[k],
                            self.results[k]))
                ddat=vec_diff(self.data[k])
                dsim=vec_diff(self.results[k])
                errs.append(cost_func(ddat,dsim))
        '''
        return np.sum(errs) 

def vec_diff(v):
    diff=[v[i+1]-v[i] for i in range(len(v)-1)]
    return diff

def get_beta_ramp(p):
    beta_array=[]
    for i in range(p.n_days):
        if i<p.Sd_delay:
            tmp=p.beta
        elif p.Sd_delay<=i<p.Sd_delay+p.n_ramp:
            tmp=p.beta*(1-p.Sd*(float((i-p.Sd_delay)/p.n_ramp)))
        else:
            tmp=p.beta*(1-p.Sd)
        beta_array.append(tmp)

    return beta_array    

def get_seasonal_beta_variation(date,peak_date,amp):
    var=amp*math.cos(2*math.pi(date-peak_date).days/365.0)
    return var

def cost_func(v0,v1):
    avg_step=3
    tmp0=[np.mean(v0[i:i+avg_step]) for i in range(len(v0)-avg_step+1)]
    tmp1=v1[:len(tmp0)]
    min_len=min((len(tmp0),len(tmp1)))
    return np.sqrt(np.sum([((tmp1[i]-tmp0[i]))**2 for i in range(min_len)]))/(min_len)

def grad_descent(p0,J0,J1,dp,alpha):
    new_p={}
    for k in p0:
        new_p[k]=p0[k]-alpha[k]*(J1[k]-J0[k])/dp[k]
    return new_p    

