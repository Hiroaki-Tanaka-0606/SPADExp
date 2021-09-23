# PAO_parser
# read .pao files (before or after optimized)

import re
import numpy as np

class PAO_parser:
    def __init__(self, filePath, optimized):
        # variables
        self.number_optpao=0
        self.grid_output=0
        self.maxL_pao=0
        self.num_pao=0
        self.PAO_Lmax=0
        self.PAO_Mul=0
        self.CCoes=None
        self.PAOs=None
        self.x=None
        self.r=None
        
        # read contraction coefficients, grid.num.output, maxL.pao, num.pao
        # read PAO.Lmax(=maxL.pao), PAO.Mul(=num.pao), pseudo.atomic.orbitals
        with open(filePath, "r") as f:
            # number.optpao
            while optimized:
                line=f.readline()
                re_result=re.findall(r"^number\.optpao\s*(\d+)", line)
                if len(re_result)>0:
                    self.number_optpao=int(re_result[0])
                    print(("number.optpao is {0:d}").format(self.number_optpao))
                    break
                re_result=re.findall(r"^\s*Input file\s*$", line)
                if len(re_result)>0:
                    print("Warning: number.optpao is not found")
                    break
                
            ccoes=[]
            # Contraction.coefficients
            for cindex in range(1, self.number_optpao+1):
                while len(re.findall(r"^<Contraction.coefficients"+str(cindex), line))==0:
                    line=f.readline()
                num_rows=int(f.readline())
                for lindex in range(0, num_rows):
                    line=f.readline()
                    re_result=re.findall(r"^\s+Atom=\s*\d+\s*"+\
                                         r"L=\s*(\d+)\s*"+\
                                         r"Mul=\s*\d+\s*"+\
                                         r"p=\s*(\d+)\s*"+\
                                         r"([0-9\.\-]+)", line)
                    if len(re_result)==0:
                        print("Error in parsing Contraction.coefficients")
                        return 
                    ccoes.append([int(re_result[0][0]),\
                                  int(re_result[0][1]),\
                                  float(re_result[0][2])])
                line=f.readline()
                if len(re.findall(r"^Contraction.coefficients"+str(cindex)+">", line))==1:
                    print(("Reading contraction coefficients {0:d} finished").format(cindex))
                else:
                    print("Error in reading Contraction.coefficients")
                    return

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

            # num.pao
            while True:
                line=f.readline()
                re_result=re.findall(r"num\.pao\s*(\d+)", line)
                if len(re_result)>0:
                    self.num_pao=int(re_result[0])
                    print(("num.pao is {:d}").format(self.num_pao))
                    break
                if len(line)==0:
                    print("Error: cannot find num.pao")
                    return
        
            # PAO.Lmax
            while True:
                line=f.readline()
                re_result=re.findall(r"PAO.Lmax\s*(\d+)", line)
                if len(re_result)>0:
                    self.PAO_Lmax=int(re_result[0])
                    print(("PAO.Lmax is {:d}").format(self.PAO_Lmax))
                    if self.PAO_Lmax!=self.maxL_pao:
                        print("Error: maxL.pao!=PAO.Lmax")
                        return
                    break
                if len(line)==0:
                    print("Error: cannot find PAO.Lmax")
                    return

            # PAO.Mul
            while True:
                line=f.readline()
                re_result=re.findall(r"PAO\.Mul\s*(\d+)", line)
                if len(re_result)>0:
                    self.PAO_Mul=int(re_result[0])
                    print(("PAO.Mul is {:d}").format(self.PAO_Mul))
                    if self.PAO_Mul!=self.num_pao:
                        print("Error: num.pao!=PAO_Mul")
                        return
                    break
                if len(line)==0:
                    print("Error: cannot find PAO_Mul")
                    return

            # put contraction coefficients in numpy
            self.CCoes=np.zeros((self.PAO_Lmax, self.PAO_Mul, self.PAO_Mul))
            Mul_indices=[]
            for i in range(0,self.PAO_Lmax):
                Mul_indices.append(-1)
                for j in range(0, self.PAO_Mul):
                    self.CCoes[i][j][j]=1
                    
            for ccoe in ccoes:
                if ccoe[1]==0:
                    Mul_indices[ccoe[0]]+=1
                i=Mul_indices[ccoe[0]]
                self.CCoes[ccoe[0]][i][ccoe[1]]=ccoe[2]

            # for mat in self.CCoes:
            # print(mat)

            # pseudo atomic orbitals
            self.x=np.zeros((self.grid_output,))
            self.r=np.zeros((self.grid_output,))
            self.PAOs=np.zeros((self.PAO_Lmax+1, self.PAO_Mul, self.grid_output))
            for l in range(0, self.PAO_Lmax+1):
                while True:
                    line=f.readline()
                    re_result=re.findall(r"^<pseudo\.atomic\.orbitals\.L="+str(l), line)
                    if len(re_result)>0:
                        break
                    if len(line)==0:
                        print("Error: cannot find pseudo atomic orbitals")
                        return
                for i in range(0, self.grid_output):
                    line_arr=f.readline().split()
                    try:
                        self.x[i]=float(line_arr[0])
                        self.r[i]=float(line_arr[1])
                        for j in range(0, self.PAO_Mul):
                            self.PAOs[l][j][i]=float(line_arr[j+2])
                    except Exception as e:
                        print(e)
                        return
                line=f.readline()
                re_result=re.findall(r"pseudo.atomic.orbitals.L="+str(l)+">", line)
                if len(re_result)>0:
                    print(("Reading pseudo atomic orbitals L= {0:d} finished").format(l))
                    
            
