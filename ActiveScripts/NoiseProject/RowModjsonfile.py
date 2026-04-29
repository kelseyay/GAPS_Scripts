import numpy as np
import os

output_dir = "output"  # Fix later!
pfile = open(os.path.join(output_dir, "RankedMods.txt"), "r")

nlayers = 7
nrows = 6
nmods = 6
nch = 32

pmax = 20
pcount = 0

byt_str = []
for l in range(nlayers):
    for r in range(nrows):
        byt_str.append([str(l) + str(r), '0'*nmods]) #Start by making an array with each lr and '0'*6 string

with open(os.path.join(output_dir, "RankedMods.json"), 'w') as jfile:
    jfile.write("{\n")

for i, x in enumerate(pfile):
    if i == 0 or i == 1:
        continue #Skip the first two lines
    tarray = x.strip().split('\t')
    pulsed1 = -1
    pulsed1 = int(tarray[3])
    pulsed2 = -1
    pulsed2 = int(tarray[4])
    #print(str(tarray[0]) + str(tarray[1]) + str(tarray[2]))

    for n in range(nlayers*nrows):
        if str(tarray[0]) + str(tarray[1]) == byt_str[n][0]:
            if (pulsed1 > -1 or pulsed2 > -1) and pcount < pmax:
                print(byt_str[n][0] + " pls pulse mod " + str(tarray[2]) + " strips " + str([pulsed1]) + " and/or " + str([pulsed2]) )
                byt_str[n][1] = byt_str[n][1][:(nmods-int(tarray[2])-1)] + str(1) + byt_str[n][1][(nmods-int(tarray[2])):] #Replace a 0 with a 1 at the module
                pcount = pcount +1
                print(byt_str[n][1])
            #print(byt_str[n][0] + " pls pulse " + str([pulsed]) )


with open(os.path.join(output_dir, "RankedMods.json"),'w') as jfile:
    jfile.write("{\n")
jfile.close()

with open(os.path.join(output_dir, "RankedMods.json"),'a') as jfile:
    for n in range(nlayers*nrows-1):
        print("lr: " + byt_str[n][0] + " bin: " + byt_str[n][1])
        #print("Associated hex: " + str(hex(int(byt_str[n][1],2)))) #The 2 means base 2!
        #jfile.write("\t\"" + byt_str[n][0] + "\": \"" + str(hex(int(byt_str[n][1],2))) + "\",\n" )
        jfile.write("\t\"" + byt_str[n][0] + "\": \"" + byt_str[n][1] + "\",\n" ) #No comma on the last element
    n = nlayers*nrows-1
    jfile.write("\t\"" + byt_str[n][0] + "\": \"" + byt_str[n][1] + "\"\n" ) #No comma on the last element
    #jfile.write("\t\"" + byt_str[n][0] + "\": \"" + str(hex(int(byt_str[n][1],2))) + "\"\n" ) #No comma on the last element

jfile.close()

with open(os.path.join(output_dir, "RankedMods.json"),'a') as jfile:
    jfile.write("\n}")

jfile.close()
