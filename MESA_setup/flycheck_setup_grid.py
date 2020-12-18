#!/usr/local/bin/python

import sys
import os
import glob
import math
import numpy as np
sys.path.append("/mnt/home/mrenzo/codes/python_stuff/plotFunc/")
from utilsLib import *


def setup_one_model(M1, M2, PERIOD, ROOT, RUNFILE):
    headerline = "export OMP_NUM_THREADS=7 && export MESA_DIR=/mnt/home/mrenzo/codes/mesa/mesa_15140/mea15140 && export MESASDK_ROOT=/mnt/home/mrenzo/codes/mesa/mesa_15140/mesasdk && source $MESASDK_ROOT/bin/mesasdk_init.sh"
    backline = " && ./clean && ./mk && ./rn 2>&1 | tee output"+"\n" # missing: move to destination
    folder_name = f"M1_{M1:.0f}_M2_{M2:.0f}_P_{PERIOD:.0f}/"
    print("----------")
    print(folder_name)
    print("----------")
    os.system('mkdir -p '+ROOT+'/'+folder_name)
    os.chdir(ROOT+'/'+folder_name)
    # # copy stuff
    copy = 'cp -r '+TEMPLATE+'/* ./'
    os.system(copy)
    with open(RUNFILE, "a") as F:
        F.writelines(headerline+" && cd "+ROOT+'/'+folder_name+backline)
    # modify the inlists
    os.system('perl -pi.back -e \'s/MASS1/'+str(M1)+'/g;\' inlist_binary')
    os.system('perl -pi.back -e \'s/MASS2/'+str(M2)+'/g;\' inlist_binary')
    os.system('perl -pi.back -e \'s/PERIOD/'+str(PERIOD)+'/g;\' inlist_binary')
    os.system('perl -pi.back -e \'s/MASS1/'+str(M1)+'/g;\' inlist1')
    os.system('perl -pi.back -e \'s/MASS2/'+str(M2)+'/g;\' inlist2')
    os.chdir(ROOT)


PERIODS = [50, 75, 100] # days
M1 = [18,19,20,21,22,23,24,25]
M2 = [15,16,17,18,19,20]

# these will all become sys.argv if needed
TEMPLATE = "/mnt/home/mrenzo/Templates/zeta_oph/MESA_setup/binary/"
ROOT = "/mnt/home/mrenzo/ceph/RUNS/ZETA_OPH/Z_0.01/"
RUNFILE=ROOT+'/runFile.txt' # to use with disBatch.py

for m1 in M1:
    
        content = checkFolder(ROOT)
        if not content:
            go_on_root = "Y"
        else:
            print(str(ROOT), "is not empty")
            print(content)
            go_on_root = input("Go on anyways? [Y/y]")
            if (go_on_root != 'Y' and go_on_root != 'y'): sys.exit()
            content = checkFolder(DESTINATION)
        if not content:
            go_on_dest = "Y"
        else:
            print(str(DESTINATION), "is not empty")
            print(content)
            go_on_dest = input("Go on anyways? [Y/y]")
            if (go_on_dest != 'Y' and go_on_dest != 'y'): sys.exit()
        if (((go_on_root == 'Y') or (go_on_root == 'y')) and ((go_on_dest == 'Y') or (go_on_dest == 'y'))):
            # now set up stuff
            description="binary semiconv and thermohaline description"+str(sem)+", q=0.8, small grid, mesa testing"
            # save template in destination
            os.system('tar -czf template.tar.xz '+TEMPLATE+' && mv template.tar.xz '+DESTINATION)
            # now write the RUNFILE
            with open(RUNFILE,"a") as F:
