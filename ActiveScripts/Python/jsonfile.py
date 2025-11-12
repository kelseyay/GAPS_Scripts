import numpy as np
import os

output_dir = "output" #Fix later!
pfile = open(os.path.join(output_dir, "RankedChannels.txt"),"r")
mfile = open("JsonMask.txt","r")
pmax = 20
pcount = 0

nlayers = 7
nrows = 6
nmods = 6
nch = 32

byt_str = []
masked_ch = []

for i, x in enumerate(mfile):
    if i == 0:
        continue
    tarray = x.strip().split('\t')
    #print(tarray[0])
    masked_ch.append(tarray[0])

print(masked_ch)

for l in range(nlayers):
    for r in range(nrows):
        for m in range(nmods):
            byt_str.append([str(l) + str(r) + str(m), '0'*nch]) #Start by making an array with each lrm and '0'*32 string which will be turned into hexadec

for i, x in enumerate(pfile):
    if i == 0 or i == 1:
        continue
    tarray = x.strip().split('\t')
    if(str(tarray[0])+str(tarray[1])+str(tarray[2])) in masked_ch:
        continue
    pulsed = -1
    pulsed = int(tarray[3])
    #print(str(tarray[0]) + str(tarray[1]) + str(tarray[2]))

    for n in range(nlayers*nrows*nmods):
        if str(tarray[0]) + str(tarray[1]) + str(tarray[2]) == byt_str[n][0]:
            if pulsed > -1 and pcount < pmax:
                print(byt_str[n][0] + " pls pulse " + str([pulsed]) )
                byt_str[n][1] = byt_str[n][1][:(nch-pulsed-1)] + str(1) + byt_str[n][1][(nch-pulsed):]
                pcount = pcount +1
                print(byt_str[n][1])
            #print(byt_str[n][0] + " pls pulse " + str([pulsed]) )


with open(os.path.join(output_dir, "RankedChannels.json"),'w') as jfile:
    jfile.write("{\n")
jfile.close()

with open(os.path.join(output_dir, "RankedChannels.json"),'a') as jfile:
    for n in range(nlayers*nrows*nmods-1):
        #print("lrm: " + byt_str[n][0] + " bin: " + byt_str[n][1])
        #print("Associated hex: " + str(hex(int(byt_str[n][1],2)))) #The 2 means base 2!
        jfile.write("\t\"" + byt_str[n][0] + "\": \"" + str(hex(int(byt_str[n][1],2))) + "\",\n" )
    n = nlayers*nrows*nmods-1
    jfile.write("\t\"" + byt_str[n][0] + "\": \"" + str(hex(int(byt_str[n][1],2))) + "\"\n" ) #No comma on the last element

jfile.close()

with open(os.path.join(output_dir, "RankedChannels.json"),'a') as jfile:
    jfile.write("\n}")

jfile.close()
