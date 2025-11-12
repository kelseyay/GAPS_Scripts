import numpy as np
import os

output_dir = "output" #Fix later!
pfile = open(os.path.join(output_dir, "best_pulsed_channel.txt"),"r")

with open(os.path.join(output_dir, "best_pulsed_channel.json"),'w') as jfile:
    jfile.write("{\n")

with open(os.path.join(output_dir, "best_pulsed_channel.json"),'a') as jfile:
    for i, x in enumerate(pfile):
        tarray = x.strip().split(' ')
        masked = [int(tarray[3]),int(tarray[4])]
        print(masked)
        byt_str = ''
        for ch in range(32):
            if ch in masked:
                byt_str += '1'
            else:
                byt_str += '0'
        byt_str = byt_str[::-1]
        byt_str = "\t\"" + str(int(tarray[0])) + str(int(tarray[1])) + str(int(tarray[2])) + "\": \"" + str(hex(int(byt_str,2))) + "\","
        #print(byt_str)
        jfile.write(byt_str)
jfile.close()

with open(os.path.join(output_dir, "best_pulsed_channel.json"),'rb+') as jfile:
    jfile.seek(-1,os.SEEK_END)
    jfile.truncate()

with open(os.path.join(output_dir, "best_pulsed_channel.json"),'a') as jfile:
    jfile.write("\n}")
