#!/usr/bin/env python
import numpy as np
import sys

def pulldata(logname):
    data={}
    with open("log.production") as f:
        lines=f.readlines()
        flag=0
        keys=[]
        for line in lines:
            if "Loop" in line:
                flag=0
                print("stop")
            if flag == 1 and "colvars" not in line:
                for key in keys:
                    data[key].append(float(line.strip().split()[loc[key]]))
            if "Step Time" in line:
                flag=1
                data={}
                loc={}
                keys=line.strip().split()
                count = 0
                for key in keys:
                    data[key]=[]
                    loc[key]=count
                    count+=1
                print("start")

    for key in data:
        data[key].pop()
        np.savetxt("%s_init.out"%key,np.c_[data[key]])
        if key == "Volume":
            L = np.array(data[key])**(1./3.)
            np.savetxt("L.dat",np.c_[L])

logfile=str(sys.argv[1])
pulldata(logfile)

