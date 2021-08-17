import csv
import json
import environment as env

infile=env.csv_pop_file
outfile=env.population_file

def csv_to_json(infile,outfile):
    csvFile=open(infile,'r')
    j=dict()
    csvReader=csv.DictReader(csvFile)
    for l in csvReader:
        state=l['State']
        county=l['County']
        pop=int(l['Population'])
        if state not in j:
            j[state]={}
            j[state][county]=pop
        else:
            j[state][county]=pop
    jsonFile=open(outfile,'w')
    jsonFile.write(json.dumps(j))

csv_to_json(infile,outfile)
