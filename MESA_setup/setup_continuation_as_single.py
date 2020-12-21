#!/usr/local/bin/python

import sys
import os
import glob
import math
import numpy as np
sys.path.append("/mnt/home/mrenzo/codes/python_stuff/plotFunc/")
from utilsLib import *

target_params = [(25,17,75),
                 (25,19,50),
                 (24,17,75),
                 (25,15,100),
                 (25,19,75),
                 (25,18,100),
                 (25,19,100),
                 (25,16,100),
                 (25,18,75)]

global ROOT
ROOT = "/mnt/home/mrenzo/ceph/RUNS/ZETA_OPH/Z_0.01/"
TEMPLATE = "/mnt/home/mrenzo/Templates/zeta_oph/MESA_setup/continue_as_single/"
folders = glob.glob(ROOT+'/*/')


def get_params(folder):
    fname = folder.split('/')[-2]
    array = fname.split('_')
    M1 = int(array[1])
    M2 = int(array[3])
    P = int(array[-1].rstrip('/'))
    return (M1, M2, P)



def setup_one_model(folder, RUNFILE=ROOT+"run_as_single_to_TAMS.txt"):
    headerline = "export OMP_NUM_THREADS=4 && export MESA_DIR=/mnt/home/mrenzo/codes/mesa/mesa_15140/mesa15140 && export MESASDK_ROOT=/mnt/home/mrenzo/codes/mesa/mesa_15140/mesasdk && source $MESASDK_ROOT/bin/mesasdk_init.sh"
    backline = " && ./clean && ./mk && ./rn 2>&1 | tee output"+"\n" # missing: move to destination
    folder_name = folder+"/accretor_to_TAMS/"
    print("----------")
    print(folder_name)
    print("----------")
    # get values
    tup = get_params(folder)
    M1 = tup[0]
    M2 = tup[1]
    P = tup[2]
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
        os.system('perl -pi.back -e \'s/MASS2/'+str(M2)+'/g;\' inlist_extra')
        os.system('perl -pi.back -e \'s/PERIOD/'+str(P)+'/g;\' inlist_extra')
    os.chdir(ROOT)


for f in folders:
    init_cond = get_params(f)
    if init_cond in target_params:
        print("will run continuation of "+f)
        setup_one_model(f, RUNFILE=ROOT+"run_as_single_to_TAMS.txt")
