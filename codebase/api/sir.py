from boto3.session import Session
from base64 import b64encode
from codebase.sir.fetch import latlon_to_place
from codebase.sir import run_comp_and_plot
from os import environ, remove

aws = Session()
ssm = aws.client('ssm')

google_key_path = '/{0}/google-api-key'.format(environ['APPLICATION'])
google_key = ssm.get_parameter(Name=google_key_path)['Parameter']['Value']

def get_plot(lat,lon,**args):
  county, state, country = latlon_to_place(google_key,lat,lon)
  
  args['state']=state
  plot_path = run_comp_and_plot(**args)
  

  with open(plot_path, 'rb') as plot_png:
    plot_b64 = b64encode(plot_png.read()).decode('utf-8')

  remove(plot_path)
  return plot_b64
  
