import uproot
import argparse
import numpy as np
from datetime import datetime
from pybfsw.gse.gsequery import GSEQuery
from tqdm import tqdm

query = GSEQuery(project='gaps')

table_dtypes = [('layer', 'u1'),
               ('row', 'u1'),
               ('module', 'u1'),
               ('channel', 'u1'),
               ('adc', 'u2'),
               ('code', 'u1'),
               ('evtid', 'u4'),
               ('evttime', 'u8'),
               ('gcutime','u4')]

parser = argparse.ArgumentParser(description = 'exports tracker db data to root file')
parser.add_argument('-t', '--times', nargs=2, required=True, type=float, help = 'unix time stamp for start and end times of extraction')
parser.add_argument('-n', '--name', help = 'name of output file')

args = parser.parse_args()

def n_packets(t1,t2):
    sql = f'select count(gcutime) from gapstrackerpacket where gcutime >= {t1} and gcutime < {t2}'
    return query.dbi.query(sql)[0][0]

def n_hits(t1,t2):
    table = query.tracker_query1a(t1, t2)
    return len(table)
def empty_tree(n):
    return np.zeros(shape = n, dtype = table_dtypes)
# we've defined a basket to contain ~ 100M MB of data or ~5,000,000 events
# we can compute an average event rate and query for those 5 Mevents at a time

t_start, t_stop = args.times
t_mid = (t_stop+t_start)/2

total_n_packets = n_packets(t_start, t_stop)

print(f'There are {total_n_packets} packets between {t_start} and {t_stop}')

dt = 500/(total_n_packets/(t_stop-t_start))

ratio = n_hits(t_mid-dt, t_mid+dt)/n_packets(t_mid-dt, t_mid+dt) # for speed we take a small section of data and figure out the ratio of packets to hits 

print(f'There are {ratio} hits per packets in the middle {2*dt} seconds of this data')

event_rate = ratio*total_n_packets/(t_stop-t_start)

print(f'The approximate trigger rate is {event_rate} hits/sec')

dt = 5e6/event_rate

intervals = np.arange(t_start, t_stop+1, dt)

file = uproot.recreate(args.name)
file.mktree('tree', table_dtypes)

for t_from, t_to in tqdm(zip(intervals[:-1], intervals[1:]), total=len(intervals)-1):
    table = query.tracker_query1a(t_from, t_to)
    tree = empty_tree(len(table))
    for i, (layer,row,module,channel,adc,evt_code,rowid,time,evtid,evttime) in enumerate(table):
        tree['layer'][i] = layer-128
        tree['row'][i] = row
        tree['module'][i] = module
        tree['channel'][i] = channel
        tree['adc'][i] = adc
        tree['code'][i] = evt_code
        tree['evtid'][i] = evtid
        tree['evttime'][i] = evttime
        tree['gcutime'][i] = time
    file['tree'].extend(tree)

file.close()