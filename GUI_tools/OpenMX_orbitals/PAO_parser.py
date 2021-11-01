# PAO_parser
# read .pao files (before or after optimized)

import re
import numpy as np
import math
import Config

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
        self.OrthNorm=None
        self.valence_min=[]
        
        
        # read contraction coefficients, grid.num.output, number.vps, <pseudo.NandL>
        # read maxL.pao, num.pao
        # read PAO.Lmax(=maxL.pao), PAO.Mul(=num.pao), pseudo.atomic.orbitals
        with open(filePath, "r") as f:
            # number.optpao
            while optimized:
                line=f.readline()
                re_result=re.findall(r"number\.optpao\s*(\d+)", line)
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
                while len(re.findall(r"<Contraction.coefficients"+str(cindex), line))==0:
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
                if len(re.findall(r"Contraction.coefficients"+str(cindex)+">", line))==1:
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

            # number.vps
            number_vps=-1
            while True:
                line=f.readline()
                re_result=re.findall(r"number\.vps\s*(\d+)", line)
                if len(re_result)>0:
                    number_vps=int(re_result[0])
                    print(("number.vps is {:d}").format(number_vps))
                    break
                if len(line)==0:
                    print("Error: cannot find number.vps")
                    return

            # <pseudo.NandL>
            pseudo_NL=[]
            maxL=-1
            
            while True:
                line=f.readline()
                re_result=re.findall(r"<pseudo.NandL", line)
                if len(re_result)>0:
                    for i in range(0, number_vps):
                        line=f.readline()
                        re_result=re.findall(r"\s*\d+\s*(\d+)\s*(\d+)", line)
                        if len(re_result)>0:
                            pseudo_NL.append([int(re_result[0][0]), int(re_result[0][1])])
                            if maxL<int(re_result[0][1]):
                                maxL=int(re_result[0][1])
                    line=f.readline()
                    re_result=re.findall("pseudo.NandL>", line)
                    if len(re_result)==0:
                        print("Error in pseudo.NandL")
                        return
                    for l in range(0, maxL+1):
                        self.valence_min.append(-1)
                    for data in pseudo_NL:
                        if self.valence_min[data[1]]<0 or self.valence_min[data[1]]>data[0]:
                            self.valence_min[data[1]]=data[0]
                    break
                        

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
            self.CCoes=np.zeros((self.PAO_Lmax+1, self.PAO_Mul, self.PAO_Mul))
            Mul_indices=[]
            for i in range(0,self.PAO_Lmax+1):
                Mul_indices.append(-1)
                for j in range(0, self.PAO_Mul):
                    self.CCoes[i][j][j]=1
                    
            for ccoe in ccoes:
                if ccoe[1]==0:
                    Mul_indices[ccoe[0]]+=1
                i=Mul_indices[ccoe[0]]
                self.CCoes[ccoe[0]][i][ccoe[1]]=ccoe[2]

            self.OrthNorm=np.zeros((self.PAO_Lmax+1, self.PAO_Mul, self.PAO_Mul))

            # for mat in self.CCoes:
            # print(mat)

            # pseudo atomic orbitals
            self.x=np.zeros((self.grid_output,))
            self.r=np.zeros((self.grid_output,))
            self.PAOs=np.zeros((self.PAO_Lmax+1, self.PAO_Mul, self.grid_output))
            for l in range(0, self.PAO_Lmax+1):
                while True:
                    line=f.readline()
                    re_result=re.findall(r"<pseudo\.atomic\.orbitals\.L="+str(l), line)
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

            # Gram-Schmidt orthonormalization
            if optimized:
                for l, CCoe_mat in enumerate(self.CCoes):
                    for i in range(1, self.PAO_Mul):
                        for j in range(0, i):
                            # orthogonalization
                            norm=np.dot(CCoe_mat[i], CCoe_mat[j])
                            CCoe_mat[i]-=norm*CCoe_mat[j]
                        
                        # normalization
                        norm2=np.dot(CCoe_mat[i], CCoe_mat[i])
                        CCoe_mat[i]/=math.sqrt(norm2)

                        
                    # print(("Orthonormalized contraction coefficients for L= {0:d}").format(l))
                    # for i in range(0, self.PAO_Mul):
                    #     for j in range(0, self.PAO_Mul):
                    #         print(("{0:7.4f} ").format(CCoe_mat[i][j]), end="")
                    #     
                    #     print("")
                    # print("")


            # Orthonormalization check of orbitals
            for l in range(0, self.PAO_Lmax+1):
                for i in range(0, self.PAO_Mul):
                    for j in range(0, i+1):
                        norm=0.0
                        for k in range(0, self.grid_output-1):
                            norm+=self.PAOs[l][i][k]*self.PAOs[l][j][k]*math.pow(self.r[k], 2)*(self.r[k+1]-self.r[k])
                        self.OrthNorm[l][i][j]=norm
                        self.OrthNorm[l][j][i]=norm


def calcContraction(after, before, matrix):
    Lmax=after.PAO_Lmax
    Mul1=after.PAO_Mul
    Mul2=before.PAO_Mul
    grid=after.grid_output

    for l in range(0, Lmax+1):
        for i in range(0, Mul1):
            for j in range(0, Mul2):
                norm1=0.0
                norm2=0.0
                for k in range(0, grid-1):
                    norm1+=after.PAOs[l][i][k]*before.PAOs[l][j][k]*math.pow(before.r[k],2)*(before.r[k+1]-before.r[k])
                    norm2+=before.PAOs[l][j][k]*before.PAOs[l][j][k]*math.pow(before.r[k],2)*(before.r[k+1]-before.r[k])
                matrix[l][i][j]=norm1/norm2

            
def reproducePAO(before, ccoes, reproduced):
    Lsize=len(ccoes)
    Mul1=len(ccoes[0])
    Mul2=before.PAO_Mul

    for l in range(0, Lsize):
        for i in range(0, Mul1):
            for j in range(0, Mul2):
                reproduced[l][i]+=before.PAOs[l][j]*ccoes[l][i][j]


def normCheck(PAO, AO):
    Lmax=PAO.PAO_Lmax
    valence_min=PAO.valence_min
    Mul=PAO.PAO_Mul
    grid=PAO.grid_output
    
    for Mul_PAO in range(0, Mul):
        for l in range(0, Lmax+1):
            Mul_AO=Mul_PAO+(valence_min[l]-(l+1) if len(valence_min)>l else 0)
            norm1=0.0
            norm2=0.0
            for k in range(round(grid*0.9), grid-1):
                norm1+=PAO.PAOs[l][Mul_PAO][k]*AO.AOs[l][Mul_AO][k]*math.pow(PAO.r[k],2)*(PAO.r[k+1]-PAO.r[k])
                norm2+=PAO.PAOs[l][Mul_PAO][k]*PAO.PAOs[l][Mul_PAO][k]*math.pow(PAO.r[k],2)*(PAO.r[k+1]-PAO.r[k])
            if norm1/norm2<Config.invert_criterion:
                print(("{0:d}{1:s} in AO is inverted, norm={2:.3f}").format(Mul_AO+l+1, Config.azimuthal[l], norm1/norm2))
                AO.AOs[l][Mul_AO]*=-1
        

            
def reproduceAO(AO, PAO, ccoes, reproduced):
    Lsize=len(ccoes)
    Mul1=len(ccoes[0])
    Mul2=PAO.PAO_Mul

    for l in range(0, Lsize):
        for i in range(0, Mul1):
            for j in range(0, Mul2):
                Mul_AO=j+(PAO.valence_min[l]-(l+1) if len(PAO.valence_min)>l else 0)
                reproduced[l][i]+=AO.AOs[l][Mul_AO]*ccoes[l][i][j]
