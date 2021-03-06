# SPADExp GUI
# class LCAO

import numpy as np
import h5py
import math

import Config
from lib import physical_tools as pt
import time
    
class LCAO:
    def __init__(self):
        self.filePath="" # file path to the OpenMX output (.hdf5)
        # /Input/UnitCells
        self.unit="" # unit of the unit cells (Ang or AU)
        self.BandCell=None # real space unit cell for the band dispersion
        self.AtomCell=None # real space unit cell for the crystal structure
        # /Input/Atomic.Species
        self.Atom_specs={} # key: atom label, value: [PAO, Pseudopotential, Orbits]
        # /Input/Atoms.SpeciesAndCoordinates
        self.Atoms=[] # list of atom labels
        self.Atom_unit="" # unit for atom positions (Ang, AU, or FRAC)
        self.Atom_coordinates=None # list of [x, y, z] or [p, q, r] (depends on unit)
        self.numAtoms=0 # number of atoms
        # calculated from /Input/Atomc.SpeciesAndCoordinates
        self.Atom_au=None # list of atom coordinates (atomic unit)
        # /Input/Kpath
        ## 2D: kx=0, 1, ..., numPnts_kx-1
        ## 3D: (kx, ky)=(0, 0), (1, 0), ..., (numPnts_kx-1, 0), (0, 1), ..., (numPnts_kx-1, numPnts_ky-1)
        ## see OpenMX_tools/preproc_main.cpp for the order of kpoints
        self.Dimension=0 # dimension in the reciprocal space (1 or 2)
        self.Curved=False # whether the band path is curved or not
        self.Origin_frac=None # the origin for the band path selection
        self.Xvector_frac=None # X vector in the fractional coordinate
        self.Xrange=None # X range
        self.numPnts_kx=1 # number of points along the kx
        self.Yvector_frac=None # Y vector in the fractional coordinate
        self.Yrange=None # Y range
        self.numPnts_ky=1 # number of points along the ky
        self.Kpath=None # list of the coordinates (fractional)
        # calculated from /Input/Kpath
        self.RecCell=None # reciprocal unit cell
        self.Xvector=None # X vector in the reciprocal coordinate
        self.Xlength=0 # length of the X vector
        self.Yvector=None # Y vector in the reciprocal coordinate
        self.Ylength=0 # length of the Y vector
        self.numPnts_k=0 # total number of points
        self.dXlength=0 # distance between two points aligned along kx
        self.dYlength=0 # distance between two points aligned along ky
        self.Kpath_au=None # list of the coordinates (atomic unit)
        # /Output
        self.Spin="" # spin polarization: off, on, nc
        self.Spin_i=0 # index corresponding to the spin polarization: 0, 1, 2 respectively
        self.EF_Eh=0 # The Fermi energy in units of the Hartree energy (27.2 eV)
        self.Band=None # band dispersion, used in off or nc
        self.BandUp=None # band dispersion, used in on
        self.BandDn=None # band dispersion, used in on
        self.numBands=0 # number of bands in the data
        self.LCAO_atoms=[] # labels "{:d}_{:s}"
        self.LCAO=[] # LCAO coefficients [atom index][orbital index](numpy[k][band][m][Re, Im])
        self.LCAO_labels=[] # LCAO labels, related to orbital index, probably alphabetic order
        self.LCAO_index=[] # [i] is inversion of self.LCAO_labels[i]
        

    def open(self, filePath):
        self.filePath=filePath

        with h5py.File(self.filePath, "r") as f:
            # /Input/UnitCells
            self.unit=f["Input"]["UnitCells"].attrs["Unit"]
            self.BandCell=np.array(f["Input"]["UnitCells"]["Bands"])
            self.AtomCell=np.array(f["Input"]["UnitCells"]["Atoms"])
            # /Input/Atomic.Species
            self.Atom_specs={}
            Atom_labels=f["Input"]["Atomic.Species"]["Labels"]
            PAOs=f["Input"]["Atomic.Species"]["PAOs"]
            PPs=f["Input"]["Atomic.Species"]["Pseudopotentials"]
            Orbits=f["Input"]["Atomic.Species"]["Orbits"]
            for i, atom_label in enumerate(Atom_labels):
                self.Atom_specs[Atom_labels[i].decode("utf_8")]=[\
                    PAOs[i].decode("utf_8"),\
                    PPs[i].decode("utf_8"),\
                    Orbits[i].decode("utf_8")]
            # /Input/Atoms.SpeciesAndCoordinates
            Atoms_h5=f["Input"]["Atoms.SpeciesAndCoordinates"]["Labels"]
            self.Atoms=[]
            for at in Atoms_h5:
                self.Atoms.append(at.decode("utf_8"))
            self.Atom_unit=f["Input"]["Atoms.SpeciesAndCoordinates"].attrs["Unit"]
            self.Atom_coordinates=np.array(f["Input"]["Atoms.SpeciesAndCoordinates"]["Coordinates"])
            self.numAtoms=f["Input"]["Atoms.SpeciesAndCoordinates"].attrs["Length"]
            # calculation from /Input/Atoms.SpeciesAndCoordinates
            self.AtomCell_au=self.AtomCell.copy()
            if self.unit.lower()=="ang":
                self.AtomCell_au/=Config.au_ang
            self.Atom_au=self.Atom_coordinates.copy()
            if self.Atom_unit.lower()=="au":
                pass
            elif self.Atom_unit.lower()=="ang":
                self.Atom_au/=Config.au_ang
            elif self.Atom_unit.lower()=="frac":
                for k in range(0, self.Atom_au.shape[0]):
                    for i in range(0, 3):
                        self.Atom_au[k][i]=0.0
                        for j in range(0, 3):
                            self.Atom_au[k][i]+=self.AtomCell_au[j][i]*self.Atom_coordinates[k][j]
            # /Input/Kpath
            self.Dimension=int(f["Input"]["Kpath"].attrs["Dimension"])
            self.Curved=f["Input"]["Kpath"].attrs["Curved"]
            self.Origin_frac=f["Input"]["Kpath"].attrs["Origin"]
            self.Xvector_frac=np.array(f["Input"]["Kpath"].attrs["Xvector"])
            self.Xrange=np.array(f["Input"]["Kpath"].attrs["Xrange"])
            self.numPnts_kx=int(f["Input"]["Kpath"].attrs["Xcount"])
            if self.Dimension==2:
                self.Yvector_frac=np.array(f["Input"]["Kpath"].attrs["Yvector"])
                self.Yrange=np.array(f["Input"]["Kpath"].attrs["Yrange"])
                self.numPnts_ky=int(f["Input"]["Kpath"].attrs["Ycount"])
            self.Kpath=np.array(f["Input"]["Kpath"]["Coordinates"])
            # calculation from /Input/Kpath
            ## The reciprocal unit cell (a.u.^-1)
            self.RecCell=np.zeros((3, 3))
            op=np.cross(self.BandCell[1], self.BandCell[2])
            det=np.inner(self.BandCell[0], op)
            ### 0
            self.RecCell[0]=2*math.pi*op/det
            ### 1
            op=np.cross(self.BandCell[2], self.BandCell[0])
            self.RecCell[1]=2*math.pi*op/det
            ### 2
            op=np.cross(self.BandCell[0], self.BandCell[1])
            self.RecCell[2]=2*math.pi*op/det
            if self.unit.lower()=="ang":
                self.RecCell*=Config.au_ang
            ## X and Y vectors
            self.Xvector=np.zeros((3,))
            for i in range(0, 3):
                self.Xvector+=self.Xvector_frac[i]*self.RecCell[i]
            self.Xlength=math.sqrt(np.inner(self.Xvector, self.Xvector))
            if self.Dimension==2:
                self.Yvector=np.zeros((3,))
                for i in range(0, 3):
                    self.Yvector+=self.Yvector_frac[i]*self.RecCell[i]
                self.Ylength=math.sqrt(np.inner(self.Yvector, self.Yvector))
            ## total number of points
            self.numPnts_k=self.numPnts_kx*self.numPnts_ky
            ## dx and dy
            self.dx_length=self.Xlength*(self.Xrange[1]-self.Xrange[0])/(self.numPnts_kx-1)
            if self.Dimension==2:
                self.dy_length=self.Ylength*(self.Yrange[1]-self.Yrange[0])/(self.numPnts_ky-1)
            ## K points
            self.Kpath_au=np.zeros(self.Kpath.shape)
            for k in range(0, self.Kpath_au.shape[0]):
                for i in range(0, 3):
                    for j in range(0, 3):
                        self.Kpath_au[k][i]+=self.RecCell[j][i]*self.Kpath[k][j]                
            # /Output
            self.Spin=str(f["Output"].attrs["Spin"])
            self.EF_Eh=float(f["Output"].attrs["EF_Eh"])
            ## band dispersion
            if self.Spin.lower()=="off":
                self.Spin_i=0
                self.Band=np.array(f["Output"]["Band"])
                self.numBands=self.Band.shape[1]
            elif self.Spin.lower()=="on":
                self.Spin_i=1
                self.BandUp=np.array(f["Output"]["BandUp"])
                self.BandDn=np.array(f["Output"]["BandDn"])
                self.numBands=self.BandUp.shape[1]
            elif self.Spin.lower()=="nc":
                self.Spin_i=2
                self.Band=np.array(f["Output"]["Band"])
                self.numBands=self.Band.shape[1]
            ## LCAO coefficients
            self.LCAO_atoms=[]
            self.LCAO=[]
            self.LCAO_labels=[]
            self.LCAO_index=[]
            for atom in f["Output"]["LCAO"].keys():
                self.LCAO_atoms.append(str(atom))
            for atom in self.LCAO_atoms:
                lcao_tmp=[]
                lcao_label_tmp=[]
                for lcao_key in f["Output"]["LCAO"][atom].keys():
                    lcao_label_tmp.append(lcao_key)
                    lcao_tmp.append(np.array(f["Output"]["LCAO"][atom][lcao_key]))
                self.LCAO.append(lcao_tmp)
                self.LCAO_labels.append(lcao_label_tmp)
                
            for LCAO_label1 in self.LCAO_labels:
                inv_LCAO_label={}
                for j, LCAO_label2 in enumerate(LCAO_label1):
                    inv_LCAO_label[LCAO_label2]=j
                self.LCAO_index.append(inv_LCAO_label)

class Wfn:
    def __init__(self, specs):
        # specs = [PAO, Pseudopotential, Orbits], same as Atom_specs[label] in LCAO
        self.PAO=specs[0]
        self.PP=specs[1]
        self.Orbits_str=specs[2] # like "s2p2d1"
        self.length=0 # data length
        self.r=None # numpy[length], values of r
        self.dr=None # numpy[length], dr[i]=r[i+1]-r[i], dr[-1]=0
        self.Wfn=[] # list of [r][PAO or AO], the order follws Orbits

        # decomposition of Orbits_str
        self.Orbits=[] # like "s0", "s1", "p0", "p1", "d0", "d1"
        self.l=[] # like 0, 0, 1, 1, 2, 2
        lList={"s": 0, "p": 1, "d": 2, "f": 3}
        i=0
        while i<len(self.Orbits_str):
            orbit_label=self.Orbits_str[i]
            mul=int(self.Orbits_str[i+1])
            for j in range(0, mul):
                self.Orbits.append(("{0:s}{1:1d}").format(orbit_label, j))
                self.l.append(lList[orbit_label])
            i+=2

        groupName=self.PAO+"&"+self.PP
        # load PAO and AO
        with h5py.File(Config.PAO_and_AO, "r") as f:
            g=f[groupName]
            self.length=g.attrs["length"]
            self.r=np.array(g.attrs["r"])
            self.dr=np.zeros(self.r.shape)
            for orbit in self.Orbits:
                self.Wfn.append(np.array(g[orbit]))
            for i in range(0, self.length-1):
                self.dr[i]=self.r[i+1]-self.r[i]
            
class PAD:
    def __init__(self, LCAO, Wfns):
        self.LCAO=LCAO
        self.Wfns=Wfns
        self.initialStates_i=0 # 0->PAO, 1->AO
        self.finalStates_i=0 # 0->Plane wave, 1->Calculated
        self.finalStates_step=0.0
        self.Y_coeff=[0, 0, 0]
        self.configured=False
        self.Gaunt=np.zeros((4,9,5,11)) # [l][m+l][lp][mp+lp]
        for lp in range(0, 5):
            for mp in range(-lp, lp+1):
                for l in range(0, 4):
                    for m in range(-l, l+1):
                        self.Gaunt[l][m+l][lp][mp+lp]=pt.Gaunt(lp, mp, l, m)

    def setSystem(self, initialStates_i, finalStates_i, finalStates_step, Y_coeff):
        self.initialStates_i=initialStates_i
        self.finalStates_i=finalStates_i
        self.finalStates_step=finalStates_step
        for i in range (0,3):
            self.Y_coeff[i]=Y_coeff[i]
        self.configured=True

    def calcFinalState(self, wfn_final, l, k, r):
        if self.configured==False:
            return
        
        if self.finalStates_i==0:
            # plane wave: spherical bessel function * r
            for i, ri in enumerate(r):
                wfn_final[i]=ri*pt.spBessel(l, k*ri)
    
    def calc(self, ik, ib, useUp=True):
        # ik: kpoint
        # ib: band
        # ia: atom
        # io: orbital of initial states
        # iL: LCAO index (found)
        # il: LCAO index (searching)
        # l: azimuthal quantum number of the initial state
        # dl: +1 or -1, l'=l+dl
        # lp: l', azimutlal quantum number of the final state
        # mpl: m+1 (0, 1, 2, ..., 2l)
        # m: magnetic quantum number of the initial state, m=mp1-l
        # jp1: j+1 (0, 1, 2)
        # j: magnetic quantum number of the perturbation, j=jp1-1

        dls=[-1, 1]
        # dls=[1]
        m1jlp=[1, -1j, -1, 1j, 1]

        k_au=self.LCAO.Kpath_au[ik]
        k_length=math.sqrt(np.inner(k_au, k_au))
        ret=0.0
        ret2=0.0
        # prepare spherical harmonics
        Ylm_k=pt.sphericalHarmonics(k_au)
        
        # for debug
        # if ik==0:
        #     print(k_au)
        #     print(Ylm_k)

        # time1=time.time()
        for ia in range(0, self.LCAO.numAtoms):
            # time2=time.time()
            # print(("calc {0:d}").format(ia), time2-time1) # 0.0035 sec 
            # time1=time2
            atom_label=self.LCAO.Atoms[ia]
            wfnObj=self.Wfns[atom_label]
            kt=np.inner(k_au, self.LCAO.Atom_au[ia])
            # print(kt)
            atom_phase=math.cos(kt)-math.sin(kt)*1j
            # prepare final states, atom dependent in case of calculated final states
            finalStates=np.zeros((5, wfnObj.length))
            for lp in range(0, 5):
                self.calcFinalState(finalStates[lp], lp, k_length, wfnObj.r)
                    
            for io, orbit_label in enumerate(wfnObj.Orbits):
                if self.LCAO.Spin_i==1:
                    # collinear spin
                    # either of Up or Dn is used (depending on argument useUp)
                    if useUp==True and orbit_label[2:]=="Dn":
                        continue
                    if useUp==False and orbit_label[2:]=="Up":
                        continue

                    if useUp==True:
                        orbit_label+="Up"
                    else:
                        orbit_label+="Dn"
                        
                l=wfnObj.l[io]
                iL=self.LCAO.LCAO_index[ia][orbit_label]
                LCAO_iL=self.LCAO.LCAO[ia][iL][ik][ib]
                LCAO_iL_conv=np.zeros((7,), dtype=complex)
                LCAO_iL_conv2=np.zeros((7,), dtype=complex) # for noncollinear spin
                # conversion of LCAO coefficients
                if l==0:
                    LCAO_iL_conv[0]=LCAO_iL[0][0]+LCAO_iL[0][1]*1j
                    if self.LCAO.Spin_i==2:
                        LCAO_iL_conv2[0]=LCAO_iL[0][2]+LCAO_iL[0][3]*1j
                elif l==1:
                    p1=LCAO_iL[0][0]+LCAO_iL[0][1]*1j
                    p2=LCAO_iL[1][0]+LCAO_iL[1][1]*1j
                    p3=LCAO_iL[2][0]+LCAO_iL[2][1]*1j
                    pt.convertLCAO_p(p1, p2, p3, LCAO_iL_conv)
                    if self.LCAO.Spin_i==2:
                        p1=LCAO_iL[0][2]+LCAO_iL[0][3]*1j
                        p2=LCAO_iL[1][2]+LCAO_iL[1][3]*1j
                        p3=LCAO_iL[2][2]+LCAO_iL[2][3]*1j
                        pt.convertLCAO_p(p1, p2, p3, LCAO_iL_conv2)                        
                elif l==2:
                    d1=LCAO_iL[0][0]+LCAO_iL[0][1]*1j
                    d2=LCAO_iL[1][0]+LCAO_iL[1][1]*1j
                    d3=LCAO_iL[2][0]+LCAO_iL[2][1]*1j
                    d4=LCAO_iL[3][0]+LCAO_iL[3][1]*1j
                    d5=LCAO_iL[4][0]+LCAO_iL[4][1]*1j
                    pt.convertLCAO_d(d1, d2, d3, d4, d5, LCAO_iL_conv)
                    if self.LCAO.Spin_i==2:                        
                        d1=LCAO_iL[0][2]+LCAO_iL[0][3]*1j
                        d2=LCAO_iL[1][2]+LCAO_iL[1][3]*1j
                        d3=LCAO_iL[2][2]+LCAO_iL[2][3]*1j
                        d4=LCAO_iL[3][2]+LCAO_iL[3][3]*1j
                        d5=LCAO_iL[4][2]+LCAO_iL[4][3]*1j
                        pt.convertLCAO_d(d1, d2, d3, d4, d5, LCAO_iL_conv2)
                elif l==3:
                    f1=LCAO_iL[0][0]+LCAO_iL[0][1]*1j
                    f2=LCAO_iL[1][0]+LCAO_iL[1][1]*1j
                    f3=LCAO_iL[2][0]+LCAO_iL[2][1]*1j
                    f4=LCAO_iL[3][0]+LCAO_iL[3][1]*1j
                    f5=LCAO_iL[4][0]+LCAO_iL[4][1]*1j
                    f6=LCAO_iL[5][0]+LCAO_iL[5][1]*1j
                    f7=LCAO_iL[6][0]+LCAO_iL[6][1]*1j
                    pt.convertLCAO_f(f1, f2, f3, f4, f5, f6, f7, LCAO_iL_conv)
                    if self.LCAO.Spin_i==3:                        
                        f1=LCAO_iL[0][2]+LCAO_iL[0][3]*1j
                        f2=LCAO_iL[1][2]+LCAO_iL[1][3]*1j
                        f3=LCAO_iL[2][2]+LCAO_iL[2][3]*1j
                        f4=LCAO_iL[3][2]+LCAO_iL[3][3]*1j
                        f5=LCAO_iL[4][2]+LCAO_iL[4][3]*1j
                        f6=LCAO_iL[5][2]+LCAO_iL[5][3]*1j
                        f7=LCAO_iL[6][2]+LCAO_iL[6][3]*1j
                        pt.convertLCAO_f(f1, f2, f3, f4, f5, f6, f7, LCAO_iL_conv2)
                wfn_initial=wfnObj.Wfn[io][:][self.initialStates_i]
                for dl in dls:
                    lp=l+dl
                    if lp<0:
                        continue
                    radialPart=(wfn_initial*finalStates[lp]*wfnObj.r*wfnObj.dr).sum()

                    coeffs=np.zeros((7,), dtype=complex)
                    for mpl in range(0, 2*l+1):
                        m=mpl-l

                        jp1St=max(-1, -(m+lp))+1 # Min of j+1
                        jp1En=min(1, lp-m)+2     # Max of j+1+1
                        mpjplpSt=jp1St-1+m+lp    # Min of m+j+lp
                        mpjplpEn=jp1En-1+m+lp    # Max of m+j+lp+1
                        # for jp1 in range(0, 3):
                        #    j=jp1-1
                        #    if m+j<-lp or m+j>lp:
                        #        continue
                        #    gaunt_coeff=self.Gaunt[lp][l][m+j+lp][m+l]
                        #    Ylpmp=Ylm_k[lp][m+j+lp]
                        #    ret2+=Ylpmp*self.Y_coeff[jp1]*gaunt_coeff
                        # ret1+=ret2*LCAO_iL_conv[mpl]

                        coeffs[mpl]=(Ylm_k[lp][mpjplpSt:mpjplpEn]*self.Gaunt[l][mpl][lp][mpjplpSt:mpjplpEn]*self.Y_coeff[jp1St:jp1En]).sum()

                    # print(coeffs)
                    ret+=(LCAO_iL_conv*coeffs).sum()*m1jlp[lp]*radialPart*atom_phase
                    if self.LCAO.Spin_i==2:
                        ret2+=(LCAO_iL_conv2*coeffs).sum()*m1jlp[lp]*radialPart*atom_phase
        # print(ret)
        if self.LCAO.Spin_i!=2:
            return ret.real**2+ret.imag**2
        else:
            return ret.real**2+ret.imag**2+ret2.real**2+ret2.imag**2


class Elements():
    def __init__(self):
        self.filePath=Config.elements_file # list of elements

        with open(self.filePath) as f:
            self.lines=f.readlines()

        count=len(self.lines)
        self.labels=[]
        self.radii=np.zeros((count, 3), dtype=float)
        self.colors=np.zeros((count, 3), dtype=float)

        index=0
        for line in self.lines:
            line_sp=line.split()
            self.labels.append(line_sp[1])
            self.radii[index][0]=float(line_sp[2])
            self.radii[index][1]=float(line_sp[3])
            self.radii[index][2]=float(line_sp[4])
            # self.colors[index][0]=round(255*float(line_sp[5]))
            # self.colors[index][1]=round(255*float(line_sp[6]))
            # self.colors[index][2]=round(255*float(line_sp[7]))
            self.colors[index][0]=float(line_sp[5])
            self.colors[index][1]=float(line_sp[6])
            self.colors[index][2]=float(line_sp[7])
            
            index+=1

        # print(self.labels)
        # print(self.radii)
        # print(self.colors)

class Dispersion():
    def __init__(self):
        self.filePath="" # file path to the SPADExp output (.hdf5)
        self.Dimension=0 # reciprocal space dimension (1 or 2)
        self.Dispersion=None # dispersion (2D or 3D)
        self.Offset=None # offset [kx, E] or [kx, ky, E]
        self.Delta=None # [delta_kx, delta_E] or [delta_kx, delta_ky, delta_E]
        self.Size=None #  [numPnts_kx, numPnts_E], [numPnts_kx, numPnts_ky, numPnts_E]
        self.Theta=0.0 # polarization angle theta
        self.Phi=0.0 # polarization angle phi
        self.Xvector=None # kx vector (a.u. orthonormal coordinate)
        self.Yvector=None # ky vector (a.u. orthonormal coordinate)
        self.Weighting=False # weighting
        self.AtomCell=None # unit cell for atoms
        self.Atom_label=None # atom labels
        self.Atom_coord=None # atom coordinates
        self.Atom_weighting=None # weighting
    
    def open(self, filePath):
        self.filePath=filePath

        with h5py.File(self.filePath, "r") as f:
            self.Dimension=f.attrs["Dimension"]
            self.Dispersion=np.array(f["Dispersion"])
            self.Offset=np.array(f.attrs["Offset"])
            self.Delta=np.array(f.attrs["Delta"])
            self.Size=np.array(f.attrs["Size"])
            self.Theta=f.attrs["Theta"]
            self.Phi=f.attrs["Phi"]
            self.Xvector=np.array(f.attrs["Xvector"])
            if self.Dimension==2:
                self.Yvector=np.array(f.attrs["Yvector"])
                
            self.AtomCell=np.array(f["Atoms"]["UnitCell"])
            self.Atom_label=[]
            Atoms_h5=f["Atoms"]["Labels"]
            for at in Atoms_h5:
                self.Atom_label.append(at.decode("utf_8"))
                
            self.Atom_coord=np.array(f["Atoms"]["Coordinates"])

            if "Weighting" in f.attrs.keys():
                self.Weighting=f.attrs["Weighting"]
            else:
                self.Weighting=False
            
            if self.Weighting==True:
                self.Atom_weighting=np.array(f["Atoms"]["Weighting"])
                
            
            
            
