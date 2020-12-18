#!/usr/local/bin/python

import sys
import os
import glob
import math
import numpy as np
sys.path.append("/mnt/home/mrenzo/codes/python_stuff/plotFunc/")
from utilsLib import *


def setup_one_model(M1, M2, PERIOD, ROOT, RUNFILE):
    headerline = "export OMP_NUM_THREADS=4 && export MESA_DIR=/mnt/home/mrenzo/codes/mesa/mesa_15140/mesa15140 && export MESASDK_ROOT=/mnt/home/mrenzo/codes/mesa/mesa_15140/mesasdk && source $MESASDK_ROOT/bin/mesasdk_init.sh"
    backline = " && ./clean && ./mk && ./rn 2>&1 | tee output"+"\n" # missing: move to destination
    folder_name = ROOT+"/"+f"M1_{M1:.0f}_M2_{M2:.0f}_P_{PERIOD:.0f}/"
    print("----------")
    print(folder_name)
    print("----------")
    # check if it exists
    content = checkFolder(folder_name)
    if not content:
        go_on = 'Y'
    else:
        print(str(folder_name), "is not empty")
        print(content)
        go_on = input("Go on anyways? [Y/y]")
    if (go_on =='Y') or (go_on =='y'):
        os.system('mkdir -p '+folder_name)
        os.chdir(folder_name)
        # # copy stuff
        copy = 'cp -r '+TEMPLATE+'/* ./'
        os.system(copy)
        with open(RUNFILE, "a") as F:
            F.writelines(headerline+" && cd "+folder_name+backline)
        # modify the inlists
        os.system('perl -pi.back -e \'s/MASS1/'+str(M1)+'/g;\' inlist_binary')
        os.system('perl -pi.back -e \'s/MASS2/'+str(M2)+'/g;\' inlist_binary')
        os.system('perl -pi.back -e \'s/PERIOD/'+str(PERIOD)+'/g;\' inlist_binary')
        os.system('perl -pi.back -e \'s/MASS1/'+str(M1)+'/g;\' inlist1')
        os.system('perl -pi.back -e \'s/MASS2/'+str(M2)+'/g;\' inlist2')
    os.chdir(ROOT)


PERIODS = [50, 75, 100] # days
M1 = [18, 19,20,21,22,23,24,25]
M2 = [15, 16,17,18,19,20]

# these will all become sys.argv if needed
TEMPLATE = "/mnt/home/mrenzo/Templates/zeta_oph/MESA_setup/binary/"
ROOT = "/mnt/home/mrenzo/ceph/RUNS/ZETA_OPH/Z_0.01/"
RUNFILE=ROOT+'/runFile.txt' # to use with disBatch.py

count_models = 0
for m1 in M1:
    for m2 in M2:
        if m2>=m1:
            # skip this combination
            continue
        else:
            for period in PERIODS:
                setup_one_model(m1, m2, period, ROOT, RUNFILE)
                count_models+=1
print(str(count_models)+" to be started")
