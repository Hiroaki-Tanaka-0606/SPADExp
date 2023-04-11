# SPADExp make_VPS_database.py

import sys
import h5py
import os
import re
from datetime import datetime
import numpy as np

if len(sys.argv)<2:
    print("Usage: %s DataDir" % sys.argv[0])
    sys.exit(0)
    
dataDir=sys.argv[1]
print("Data directory is %s" % dataDir)

VPSDir=os.path.join(dataDir, "VPS")
print("Target directory is %s" % VPSDir)

ffList=os.listdir(VPSDir)

VPSPathList=[]
VPSNameList=[]
for ff in ffList:
    path=os.path.join(VPSDir, ff)
    re_result=re.findall(r"^(.*)\.vps$", ff)
    if os.path.isfile(path) and len(re_result)==1:
        VPSPathList.append(path)
        VPSNameList.append(re_result[0])

with h5py.File("../Pseudopotentials.hdf5", "w") as h5f:
    h5f.attrs.create("datetime", datetime.now().isoformat(" "))
    
    for i, VPSPath in enumerate(VPSPathList):
        print("Read %s" % VPSPath)
        VPSName=VPSNameList[i]
        VPSG=h5f.create_group(VPSName)

        num_grid=0
        num_vps_proj=0
        num_vps=0
        cutoff_list=None
        cutoff_max=-1
        j_depend=False
        j_digit=1
        l_list=None # [l]
        E_list=None # [l][j]
        x=None # 1st column log(r) 
        r=None # 2nd column 
        V_loc=None # [r]
        V_nonloc=None # [l][j][r]
        np_set1=False # for x, r, V_loc (<=> num_grid)
        np_set2=False # for l_list, E_list, V_nonloc (<=> num_grid and num_vps_proj)
        np_set3=False # for cutoff_list (<=> num_vps)
        with open(VPSPath) as f:
            while True:
                line_raw=f.readline()
                if len(line_raw)==0:
                    print("End of file")
                    break
                line=line_raw.strip()
                line_sp=line.split()
                if re.findall(r"^grid\.num\.output", line):
                    num_grid=int(line_sp[1])
                    print("Number of grid is %d" % num_grid)
                    x=np.zeros((num_grid,), dtype=float)
                    r=np.zeros((num_grid,), dtype=float)
                    V_loc=np.zeros((num_grid,), dtype=float)
                    np_set1=True
                elif re.findall(r"^j\.dependent\.pseudo\.potentials", line):
                    j_depend=line_sp[1]=="on"
                    j_digit=2 if j_depend else 1
                    print("j-dependence is %s (%s)" % (line_sp[1], "true" if j_depend else "false"))
                elif re.findall(r"^number\.vps", line):
                    num_vps=int(line_sp[1])
                    print("Number of vps is %d" % num_vps)
                    cutoff_list=np.zeros((num_vps,), dtype=float)
                    np_set3=True
                elif re.findall(r"^<pseudo.NandL", line):
                    if np_set3==False:
                        print("Error: num_vps not set")
                        break
                    for i in range(num_vps):
                        line_sp=f.readline().strip().split()
                        cutoff_list[i]=float(line_sp[3])
                        if cutoff_list[i]>cutoff_max:
                            cutoff_max=cutoff_list[i]
                    line=f.readline().strip()
                    if re.findall(r"pseudo.NandL>", line):
                        print("<pseudo.NandL> loaded")
                        print("Maximum cutoff is %.1f Bohr" % cutoff_max)
                    else:
                        print("Error in <pseudo.NandL>")
                        break
                    
                elif re.findall(r"^<project.energies", line):
                    num_vps_proj=int(f.readline().strip().split()[0])
                    print("Number of vps projector is %d" % num_vps_proj)
                    if np_set1==False:
                        print("Error: num_grid not set")
                        break
                    l_list=np.zeros((num_vps_proj,), dtype=int)
                    E_list=np.zeros((num_vps_proj, j_digit), dtype=float)
                    V_nonloc=np.zeros((num_vps_proj, j_digit, num_grid), dtype=float)
                    np_set2=True
                    for i in range(num_vps_proj):
                        line_sp=f.readline().strip().split()
                        l_list[i]=int(line_sp[0])
                        for j in range(j_digit):
                            E_list[i][j]=float(line_sp[j+1])
                    line=f.readline().strip()
                    if re.findall(r"^project.energies>", line):
                        print("<project.energies> loaded")
                    else:
                        print("Error in <project.energies>")
                        break
                    
                elif re.findall(r"^<Pseudo.Potentials", line):
                    if np_set1==False or np_set2==False:
                        print("Error: something is not set")
                        break
                    for i in range(num_grid):
                        line_sp=f.readline().strip().split()
                        x[i]=float(line_sp[0])
                        r[i]=float(line_sp[1])
                        V_loc[i]=float(line_sp[2])
                        for l in range(num_vps_proj):
                            for j in range(j_digit):
                                V_nonloc[l][j][i]=line_sp[3+l*j_digit+j]
                    line=f.readline().strip()
                    if re.findall(r"Pseudo.Potentials>", line):
                        print("<Pseudo.Potentials> loaded")
                        # print("%.3f" % V_loc[0])
                        # for l in range(num_vps_proj):
                        # for j in range(j_digit):
                        # print("%.3e" % V_nonloc[l][j][0])
                    else:
                        print("Error in <Pseudo.Potentials>")
                        break
        VPSG.attrs.create("num_vps", num_vps)
        VPSG.attrs.create("num_vps_proj", num_vps_proj)
        VPSG.attrs.create("num_grid", num_grid)
        VPSG.attrs.create("j_dependent", j_depend)
        VPSG.attrs.create("r", r)
        VPSG.attrs.create("max_cutoff", cutoff_max)
        VPSG.create_dataset("orbital_angular_momenta", data=l_list)
        VPSG.create_dataset("project_energies", data=E_list)
        VPSG.create_dataset("local_potential", data=V_loc)
        VPSG.create_dataset("nonlocal_potentials", data=V_nonloc)
        
    
        
