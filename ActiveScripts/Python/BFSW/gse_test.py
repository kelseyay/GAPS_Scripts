#ssh -p 55227 -L 44555:localhost:44555 gse@gamma1.ssl.berkeley.edu
#Password is th@tonE. You know. 
#Need to source the python environment!
#source /home/kelsey/bfsw_env/bin/activate
#Need to connect to gse7!
#python3 gse_test.py

#import uproot
import argparse
import numpy as np
from datetime import datetime
from pybfsw.gse.gsequery import GSEQuery #It's a filepath
from tqdm import tqdm
from pybfsw.gse import gsequery

t1 =  1766869931
#1767981712
t2 = 1766870031
#1767981812

print("Hello world!")
query = GSEQuery(project = "gaps") #This will fail if gse7 isn't being spoken to properly
#data = query.time_query1("@asictemp_l0r0m0", t1, t2)

#Third time's the charm
print(query.time_query3("@asictemp_l0r1m5", t1,t2))

#data = query.get_latest_tracker_temp(1, 1, 1)
#print(data)

#Can't remember how to get all of the names!
#data = query.time_query1("@asictemp_l0r0m0", t1, t2)
#data = query.get_latest_tracker_temp(1, 0, 0)
#data = query.
#print(data)

#transfer functions should be in
