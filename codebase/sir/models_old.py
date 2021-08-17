import numpy as np
import math
from datetime import date,timedelta,datetime
from scipy import stats
from codebase.sir import estimation
from codebase.sir import fetch

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
                 data_days,
                 s_idx,e_idx):
        
        self.params=params()
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
        self.params.s_idx=s_idx
        self.params.e_idx=e_idx
        self.params.search_steps=20
        self.params.descent_steps=500

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
        self.data['I_ref']=self.data['I_current']
        self.params.data_days=len(self.data['dates'])
        self.params.Sd_delay=min(((self.params.Sd_date-self.data['dates'][0]).days,self.params.data_days))

        self.params.n_days=self.params.data_days+self.params.n_days
        
        vals=estimation.cases_to_intercepts(self.data,self.params)   
        v=self.init_vals
        v.I=vals['I_current']
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
        #self.grid_search_param_fit()
        self.grad_descent_param_fit()
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
        if vn.I<0.5: vn.I=0.0
        if vn.R<0.0: vn.R=0.0
    
        return vn

    def grad_descent_param_fit(self):
        params=self.params
        v=self.init_vals

        p0=[0.5,params.Sd_delay,10]
        p1=p0
        dp=[1e-3,1e-3,1e-3]
        alpha=[1e-6,1e-8,1e-8]

        cost0=1
        cost1=0
        step=1
        
        while step<params.descent_steps and (cost1-cost0)!=0:   
           cost_diff=cost1-cost0
           print("cost diff: %s"%(cost_diff))
           if cost_diff<0:
               alpha=[x*1.05 for x in alpha]
               p0=p1
           else:
               alpha=[x*0.5 for x in alpha]

           params.Sd=p0[0]
           params.Sd_delay=p0[1]
           params.n_ramp=p0[2]
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           
           cost0=self.calculate_error()
           print("cost: %s"%(cost0))
           
           J0=[cost0]*len(p0)
           J1=[0]*len(p0)

           params.Sd=p0[0]+dp[0]
           params.Sd_delay=p0[1]
           params.n_ramp=p0[2]
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           J1[0]=self.calculate_error()
           
           params.Sd=p0[0]
           params.Sd_delay=p0[1]+dp[1]
           params.n_ramp=p0[2]
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           J1[1]=self.calculate_error()
           
           params.Sd=p0[0]
           params.Sd_delay=p0[1]
           params.n_ramp=p0[2]+dp[2]
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           J1[2]=self.calculate_error()

           p1=grad_descent(p0,J0,J1,dp,alpha)
           params.Sd=p1[0]
           params.Sd_delay=p1[1]
           params.n_ramp=p1[2]
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           cost1=self.calculate_error()

           step+=1
           #print(p0)

        self.params.Sd=p0[0]
        self.params.Sd_delay=p0[1]
        self.params.n_ramp=p0[2]
        return self.params

    def grid_search_param_fit(self):
        params=self.params
        v=self.init_vals
        N=v.S+v.I+v.R
         
        err=np.inf

        fits={'Sd':0.0,'n_ramp':0,'Sd_delay':0.0}

        Sd_delay_min=max((params.Sd_delay-2,0))
        Sd_delay_max=params.Sd_delay+2

        for Sd in np.linspace(0.1,0.9,params.search_steps):
            for Sd_delay in np.arange(Sd_delay_min,Sd_delay_max):
                for n_ramp in np.arange(1,20):
                    params.Sd=Sd
                    params.Sd_delay=Sd_delay
                    params.n_ramp=n_ramp
                    beta_array=get_beta_ramp(params)
                    self.run_model(params.data_days,
                                   beta_array)
                    tmp=self.calculate_error()
                    if tmp<err: 
                       #print("error,beta,Sd,Sd_delay,n_ramp: %s,%s,%s,%s,%s"%(tmp,beta,Sd,Sd_delay,n_ramp))
                       err=tmp
                       fits['Sd']=Sd
                       fits['n_ramp']=n_ramp
                       fits['Sd_delay']=Sd_delay
        
        self.params.n_ramp=fits['n_ramp']
        self.params.Sd=fits['Sd']
        self.params.Sd_delay=fits['Sd_delay']
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
        dv0=[self.data['I_ref'][i+1]-self.data['I_ref'][i] for i in range(self.params.data_days-1)]
        dv1=[self.results['I'][i+1]-self.results['I'][i] for i in range(self.params.data_days-1)]
        v0=self.data['I_ref']
        v1=self.results['I']
        
        val_cost=cost_func(v0,v1)
        deriv_cost=cost_func(dv0,dv1)
        
        return val_cost+deriv_cost

class SEAIQHRD:
    def __init__(self,state,county,
                 n_days,n_substeps,
                 n_ramp,
                 Tea,Tai,Tir,
                 Tiq,Tid,Tqh,
                 Tih,pai,piq,
                 Sd,Sd_delay,
                 detection_rate,
                 data_days,
                 s_idx,e_idx):
        
        self.params=params()
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
        self.params.Tqh=Tqh
        self.params.Sd=Sd
        self.params.Sd_delay=Sd_delay
        self.params.n_ramp=n_ramp
        self.params.detection_rate=detection_rate
        self.params.data_days=data_days
        self.params.piq=piq
        self.params.pai=pai
        self.params.s_idx=s_idx
        self.params.e_idx=e_idx
        self.params.search_steps=20
        self.params.descent_steps=500

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
        self.data=fetch.refine_data(self.data,self.params)
        self.data['I_ref']=self.data['I_current']-self.data['H']
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
        self.rates=fetch.data_to_rates(self.data,self.params)

    def update_parameters(self):
       
        p=self.params
        p.pid=self.rates.pid
        p.pqd=self.rates.pqd
        p.phd=self.rates.phd
        p.pih=self.rates.pih
        p.pqh=self.rates.pqh
        dt=p.dt
        p.Tar=p.Tai+p.Tir
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
        
        self.params=estimation.cases_to_beta(self.data,p,self.init_vals)
        self.grad_descent_param_fit()
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
        if vn.A<0.5: vn.A=0.0
        if vn.I<0.5: vn.I=0.0
        if vn.Q<0.5: vn.Q=0.0
        if vn.H<0.0: vn.H=0.0
        if vn.R<0.0: vn.R=0.0
        if vn.D<0.0: vn.D=0.0
    
        return vn


    def grad_descent_param_fit(self):
        params=self.params
        v=self.init_vals

        p0=[0.5,params.Sd_delay,10]
        p1=p0
        dp=[1e-3,1e-3,1e-3]
        alpha=[1e-10,1e-10,1e-10]

        cost0=1
        cost1=0
        step=1
        
        while np.abs(cost1-cost0)>1e-6:   
           cost_diff=cost1-cost0
           print("cost diff: %s"%(cost_diff))
           if cost_diff<0:
               alpha=[x*1.05 for x in alpha]
               p0=p1
           else:
               alpha=[x*0.5 for x in alpha]

           params.Sd=p0[0]
           params.Sd_delay=p0[1]
           params.n_ramp=p0[2]
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           
           cost0=self.calculate_error()
           print("cost: %s"%(cost0))
           
           J0=[cost0]*len(p0)
           J1=[0]*len(p0)

           params.Sd=p0[0]+dp[0]
           params.Sd_delay=p0[1]
           params.n_ramp=p0[2]
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           J1[0]=self.calculate_error()
           
           params.Sd=p0[0]
           params.Sd_delay=p0[1]+dp[1]
           params.n_ramp=p0[2]
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           J1[1]=self.calculate_error()
           
           params.Sd=p0[0]
           params.Sd_delay=p0[1]
           params.n_ramp=p0[2]+dp[2]
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           J1[2]=self.calculate_error()

           p1=grad_descent(p0,J0,J1,dp,alpha)
           params.Sd=p1[0]
           params.Sd_delay=p1[1]
           params.n_ramp=p1[2]
           beta_array=get_beta_ramp(params)
           self.run_model(params.data_days,
                          beta_array)
           cost1=self.calculate_error()

           step+=1
           #print(p0)

        self.params.Sd=p0[0]
        self.params.Sd_delay=p0[1]
        self.params.n_ramp=p0[2]
        return self.params

    def simulate(self):
        return self.run_model(self.params.n_days,
                       self.params.beta_array)

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
        I_sim=self.results['I']+self.results['A']+self.results['Q']+self.results['E']
        dv0=[self.data['I_ref'][i+1]-self.data['I_ref'][i] for i in range(self.params.data_days-1)]
        dv1=[I_sim[i+1]-I_sim[i] for i in range(self.params.data_days-1)]
        v0=self.data['I_ref']
        v1=I_sim
        
        val_cost=cost_func(v0,v1)
        deriv_cost=cost_func(dv0,dv1)
        
        return val_cost#+deriv_cost
 

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
    min_len=min((len(v0),len(v1)))
    return np.sqrt(np.sum([((v1[i]-v0[i]))**2 for i in range(min_len)]))/(min_len)

def grad_descent(p0,J0,J1,dp,alpha):
    new_p=[0]*len(p0)
    for i in range(len(p0)): 
        new_p[i]=p0[i]-alpha[i]*(J1[i]-J0[i])/dp[i]
    return new_p    
