from codebase import get_logger
from datetime import date, timedelta, datetime
from os import environ, path, system
import json
import numpy as np
import urllib.request
from scipy import stats
from tempfile import gettempdir

logger = get_logger()

pop_path = path.join(environ['LAMBDA_TASK_ROOT'], 'pop.json')
states_path = path.join(environ['LAMBDA_TASK_ROOT'], 'us_states.json')

cases_path = path.join(gettempdir(), 'case_data_%s-%s-%s.json'%(date.today().day,date.today().month,date.today().year))

params_path = path.join(gettempdir(), 'params_%s-%s-%s.json'%(date.today().day,date.today().month,date.today().year))

orders_path = path.join(environ['LAMBDA_TASK_ROOT'], 'sah_order_dates.json')

covid_uri="https://covidtracking.com/api/states/daily?state={0}"

class initial_values:
    pass

with open(pop_path, 'rt') as pop_json:
    pop_data = json.load(pop_json)

with open(states_path, 'rt') as states_json:
    us_states = json.load(states_json)

with open(orders_path, 'rt') as orders_json:
    order_dates = json.load(orders_json)


def download_data_file(case_file):    
    url="https://covidtracking.com/api/v1/states/daily.json"
    if not path.exists(case_file):
        system('rm -rf %s/case_data*'%(gettempdir()))
        logger.info("Downloading case data: %s"%(case_file))
        v=urllib.request.urlopen(url).read()
        j=json.loads(v)
        with open(case_file, 'w') as outfile:
            json.dump(j,outfile)
    else:
        logger.info("Case data already downloaded: %s"%(case_file))
 
download_data_file(cases_path)
with open(cases_path, 'rt') as data_json:
    case_data = json.load(data_json)       

def latlon_to_place(api_key, lat, lon):
    uri = 'https://maps.googleapis.com/maps/api/geocode/json?latlng={Lat},{Lon}&key={ApiKey}'.format(
        ApiKey=api_key,
        Lat=lat,
        Lon=lon,
    )
    
    logger.info('Retrieving GeoCoding for ({0}, {1})...'.format(lat, lon))
    
    v=urllib.request.urlopen(uri).read()
    logger.info(v)
    j=json.loads(v)
    
    components=j['results'][0]['address_components']
    country=state=county=None
    for c in components:
        if "country" in c['types']:
            country=c['long_name']
        if "administrative_area_level_1" in c['types']:
            state=c['long_name']
        if "administrative_area_level_2" in c['types']:
            county=c['long_name']

    return county,state,country        

def format_date(d):
    return datetime.strptime(str(d),'%Y%m%d').date()

def location_to_population(state,county=None):
    if county==None: county=state
    N=pop_data[state][county] 
    return N
    
def location_to_data(params):
    
    abbrev=us_states[params.state]
    data={}
    
    I_total=np.array([case_data[i].get("positive",0) for i in range(len(case_data)) if case_data[i]['state']==abbrev])
        
    I_new=np.array([case_data[i].get("positiveIncrease",0) for i in range(len(case_data)) if case_data[i]['state']==abbrev])

    D=np.array([case_data[i].get("death",0) for i in range(len(case_data)) if case_data[i]['state']==abbrev])
    
    R=np.array([case_data[i].get("recovered",0) for i in range(len(case_data)) if case_data[i]['state']==abbrev])
    
    H_total=np.array([case_data[i].get("hospitalizedCumulative",0) for i in range(len(case_data)) if case_data[i]['state']==abbrev])
    
    H=np.array([case_data[i].get("hospitalizedCurrently",0) for i in range(len(case_data)) if case_data[i]['state']==abbrev])
     
    dates=np.array([case_data[i].get("date",0) for i in range(len(case_data)) if case_data[i]['state']==abbrev])
    dates=[format_date(x) for x in dates]

    I_total[np.where(I_total==None)[0]]=0
    data['I_cum']=I_total[::-1]
    
    I_new[np.where(I_new==None)[0]]=0
    data['I_new']=I_new[::-1]
    
    D[np.where(D==None)[0]]=0
    data['D']=D[::-1]
    
    R[np.where(R==None)[0]]=0
    data['R']=R[::-1]
    
    H[np.where(H==None)[0]]=0
    data['H']=H[::-1]
    
    H_total[np.where(H_total==None)[0]]=0
    data['H_cum']=H_total[::-1]
    
    data['I_act']=data['I_cum']-data['R']-data['D']
    
    data['I']=data['I_act']

    data['dates']=dates[::-1]  
    return data

def refine_data(data,params):    

    #approx detected active cases
    tmp=(data['I_act']-data['H'])

    data['I_act']=data['I_act']/params.detection_rate
    data['I_cum']=data['I_cum']/params.detection_rate
    
    data['I_new']=data['I_new']/params.detection_rate
    
    data['R']=data['R']/params.detection_rate
    
    #approx undetected cases
    data['I_total']=tmp/params.detection_rate
    
    data['A']=(data['I_total'])*params.rai/(1+params.rai)

    data['I']=np.array([0]*len(data['A']))
    data['Q']=np.array([0]*len(data['A']))

    for n in range(1,len(data['A'])):
        data['I'][n]=data['I'][n-1]*(1-params.i_decay)+params.i_growth*data['A'][n-1]
        if data['I'][n]<0: data['I'][n]=0
        data['Q'][n]=data['Q'][n-1]*(1-params.q_decay)+params.q_growth*data['I'][n-1]
        if data['Q'][n]<0: data['Q'][n]=0
        
    return data

def data_to_initial_values(data,params):
    v=initial_values()   

    #initialize compartments
    v.E=0
    v.Q=data['Q'][-1]
    v.R=data['R'][-1]
    v.I=data['I'][-1]
    v.A=data['A'][-1]
    v.H=data['H'][-1]
    v.D=data['D'][-1]
    v.N=params.N
    v.S=v.N-v.A-v.I-v.Q-v.E-v.R-v.D-v.H

    return v

def data_to_rates(data,params):
    
    #death rate of symptomatic
    r=initial_values()
    
    min_death_rate=0.001
    r.death_rate=max((params.detection_rate*data['D'][-1]/(data['I_cum'][-1]/(1+params.rai)),min_death_rate))
    r.pid=r.pqd=r.phd=r.death_rate

    #hosp rate of symptomatic
    min_hosp_rate=0.001
    r.hosp_rate=max((params.detection_rate*data['H_cum'][-1]/(data['I_cum'][-1]/(1+params.rai)),min_hosp_rate))
    
    r.pih=r.pqh=r.hosp_rate

    return r

