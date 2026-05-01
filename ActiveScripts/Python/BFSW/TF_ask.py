#Testing Gabriel's TF instructions

import argparse
import numpy as np
from datetime import datetime
from pybfsw.gse.gsequery import GSEQuery #It's a filepath
from tqdm import tqdm
from pybfsw.gse import gsequery

def query_calibrations(query, timestamps):
    table = []
    for ts in timestamps:
        table += query.get_tracker_cal_data(ts)
    dac_pulse_heights = np.unique([row[12] for row in table])
    calibration_data = np.zeros((7,6,6,32,len(cal_pulse_heights),2))
    for row in table:
        time = row[2]
        l, r, m = row[7:10]
        dac = row[12]
        i = np.where(cal_pulse_heights==dac)[0][0]
        for ch in range(32):
            calibration_data[l,r,m,ch,i,0] = row[5*ch + 17]/8
            calibration_data[l,r,m,ch,i,1] = row[5*ch + 16]/8

    return dac_pulse_heights, calibration_data
