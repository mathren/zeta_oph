#!/usr/local/bin/python

import sys
import os
import glob
import math
import numpy as np
sys.path.append("/mnt/home/mrenzo/codes/python_stuff/plotFunc/")
from utilsLib import *

def get_params(folder):
    array = folder.split('_')
    M1 = int(array[1])
    M2 = int(array[3])
    P = int(array[-1].rstrip('/'))
    return (M1, M2, P)

target_params = [(25,17,75),
                 (25,19,50),
                 (24,17,75),
                 (25,15,100),
                 (25,19,75),
                 (25,18,100),
                 (25,19,100),
                 (25,16,100),
                 (25,18,75)]


root = "/mnt/home/mrenzo/ceph/RUNS/ZETA_OPH/Z_0.01"
folders = glob.glob(root+'/*/')
print(folders)
for f in folders:
    init_cond = get_params(f)
    print(init_cond)
    if init_cond in target_params:
        print("will run continuation of "+f)
        
# def setup_one_model(M1, M2, PERIOD, ROOT, RUNFILE):
#     headerline = "export OMP_NUM_THREADS=4 && export MESA_DIR=/mnt/home/mrenzo/codes/mesa/mesa_15140/mesa15140 && export MESASDK_ROOT=/mnt/home/mrenzo/codes/mesa/mesa_15140/mesasdk && source $MESASDK_ROOT/bin/mesasdk_init.sh"
#     backline = " && ./clean && ./mk && ./rn 2>&1 | tee output"+"\n" # missing: move to destination
#     folder_name = ROOT+"/"+f"M1_{M1:.0f}_M2_{M2:.0f}_P_{PERIOD:.0f}/"
#     print("----------")
#     print(folder_name)
#     print("----------")
#     # check if it exists
#     content = checkFolder(folder_name)
#     if not content:
#         go_on = 'Y'
#     else:
#         print(str(folder_name), "is not empty")
#         print(content)
#         go_on = input("Go on anyways? [Y/y]")
#     if (go_on =='Y') or (go_on =='y'):
#         os.system('mkdir -p '+folder_name)
#         os.chdir(folder_name)
#         # # copy stuff
#         copy = 'cp -r '+TEMPLATE+'/* ./'
#         os.system(copy)
#         with open(RUNFILE, "a") as F:
#             F.writelines(headerline+" && cd "+folder_name+backline)
#         # modify the inlists
#         os.system('perl -pi.back -e \'s/MASS1/'+str(M1)+'/g;\' inlist_binary')
#         os.system('perl -pi.back -e \'s/MASS2/'+str(M2)+'/g;\' inlist_binary')
#         os.system('perl -pi.back -e \'s/PERIOD/'+str(PERIOD)+'/g;\' inlist_binary')
#         os.system('perl -pi.back -e \'s/MASS1/'+str(M1)+'/g;\' inlist1')
#         os.system('perl -pi.back -e \'s/MASS2/'+str(M2)+'/g;\' inlist2')
#     os.chdir(ROOT)
