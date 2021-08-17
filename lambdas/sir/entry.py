from codebase import get_logger, json_handler
from codebase.api.responses import binary_response, client_error, critical_error, server_error
from codebase.api.exceptions import ServerlessClientException, ServerlessException
from codebase.api.sir import get_plot
from os import environ, path
import json
import logging

logger = get_logger()

def lambda_handler(event, context):
  logger.debug('Request ' + json.dumps(event, default=json_handler))
  
  try:
    lat = event['pathParameters']['lat']
    lon = event['pathParameters']['lon']
    n_days = event['pathParameters']['n_days']
    detection_rate = event['pathParameters']['detection_rate']
    piq = event['pathParameters']['piq']
    pai = event['pathParameters']['pai']
    refit = event['pathParameters']['refit']
    
    if refit=="True": refit=True
    if refit=="False": refit=False

    args={
          'n_days':int(n_days),
          'piq':float(piq),
          'pai':float(pai),
          'detection_rate':float(detection_rate),
          'refit':refit,
          }
    
    return binary_response(
      Body=get_plot(lat,lon,**args),
      ContentType='image/png',
    )

  except ServerlessClientException as e:
    return client_error(e)
  
  except ServerlessException as e:
    return server_error(e)
  
  except Exception as e:
    return critical_error(e)
