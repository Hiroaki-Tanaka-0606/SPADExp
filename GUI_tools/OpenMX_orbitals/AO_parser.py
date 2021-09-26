# AO_parser
# read .ao files (from VPS)

import re
import numpy as np
import math

class AO_parser:
    def __init__(self, filePath):
        # variables
        self.grid_output=0
        self.maxL_pao=0
        self.max_N=0
        self.AOs=None
        self.x=None
        self.r=None
        self.OrthNorm=None
        
        # grid.num.output, maxL.pao, max_N
        # read atomic.orbitals
        with open(filePath, "r") as f:
            # grid.num.output
            while True:
                line=f.readline()
                re_result=re.findall(r"grid\.num\.output\s*(\d+)", line)
                if len(re_result)>0:
                    self.grid_output=int(re_result[0])
                    print(("grid.num.output is {:d}").format(self.grid_output))
                    break
                if len(line)==0:
                    print("Error: cannot find grid.num.output")
                    return

            # maxL.pao
            while True:
                line=f.readline()
                re_result=re.findall(r"maxL\.pao\s*(\d+)", line)
                if len(re_result)>0:
                    self.maxL_pao=int(re_result[0])
                    print(("maxL.pao is {:d}").format(self.maxL_pao))
                    break
                if len(line)==0:
                    print("Error: cannot find maxL.pao")
                    return
                
            # max.N
            while True:
                line=f.readline()
                re_result=re.findall(r"max\.N\s*(\d+)", line)
                if len(re_result)>0:
                    self.max_N=int(re_result[0])
                    print(("max.N is {:d}").format(self.max_N))
                    break
                if len(line)==0:
                    print("Error: cannot find max.N")
                    return

            self.OrthNorm=np.zeros((self.maxL_pao+1, self.max_N, self.max_N))

            # atomic orbitals
            self.x=np.zeros((self.grid_output,))
            self.r=np.zeros((self.grid_output,))
            self.AOs=np.zeros((self.maxL_pao+1, self.max_N, self.grid_output))
            for l in range(0, self.maxL_pao+1):
                while True:
                    line=f.readline()
                    re_result=re.findall(r"^<atomic\.orbitals\.L="+str(l), line)
                    if len(re_result)>0:
                        break
                    if len(line)==0:
                        print("Error: cannot find atomic orbitals")
                        return
                for i in range(0, self.grid_output):
                    line_arr=f.readline().split()
                    try:
                        self.x[i]=float(line_arr[0])
                        self.r[i]=float(line_arr[1])
                        for j in range(0, min(self.max_N, len(line_arr)-2)):
                            self.AOs[l][j][i]=float(line_arr[j+2])
                    except Exception as e:
                        print(e)
                        return
                line=f.readline()
                re_result=re.findall(r"atomic.orbitals.L="+str(l)+">", line)
                if len(re_result)>0:
                    print(("Reading atomic orbitals L= {0:d} finished").format(l))

            # normalization
            for l in range(0, self.maxL_pao+1):
                for i in range(0, self.max_N-l):
                    norm=0.0
                    for k in range(0, self.grid_output-1):
                        norm+=math.pow(self.AOs[l][i][k],2)*math.pow(self.r[k], 2)*(self.r[k+1]-self.r[k])
                        
                    self.AOs[l][i]/=math.sqrt(norm)
                    
            # Orthonormalization check or orbitals
            for l in range(0, self.maxL_pao+1):
                for i in range(0, self.max_N):
                    for j in range(0, i+1):
                        norm=0.0
                        for k in range(0, self.grid_output-1):
                            norm+=self.AOs[l][i][k]*self.AOs[l][j][k]*math.pow(self.r[k], 2)*(self.r[k+1]-self.r[k])
                        self.OrthNorm[l][i][j]=norm
                        self.OrthNorm[l][j][i]=norm
