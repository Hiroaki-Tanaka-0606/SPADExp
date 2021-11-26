# calcPSF GUI
# class LCAO

import numpy as np
import h5py
import math

import Config
from lib import physical_tools as pt
    
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
            self.Atom_au=self.Atom_coordinates.copy()
            if self.Atom_unit.lower()=="au":
                pass
            elif self.Atom_unit.lower()=="ang":
                self.Atom_au/=Config.au_ang
            elif self.Atom_unit.lower()=="frac":
                atomCell_au=self.AtomCell.copy()
                if self.unit.lower()=="ang":
                    atomCell_au/=Config.au_ang
                for k in range(0, self.Atom_au.shape[0]):
                    for i in range(0, 3):
                        self.Atom_au[k][i]=0.0
                        for j in range(0, 3):
                            self.Atom_au[k][i]+=atomCell_au[j][i]*self.Atom_coordinates[k][j]
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
            ## The reciprocal unit cell
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
            recCell_au=self.RecCell.copy()
            if self.unit.lower()=="ang":
                recCell_au*=Config.au_ang
            self.Kpath_au=np.zeros(self.Kpath.shape)
            for k in range(0, self.Kpath_au.shape[0]):
                for i in range(0, 3):
                    for j in range(0, 3):
                        self.Kpath_au[k][i]+=recCell_au[j][i]*self.Kpath[k][j]                
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


class Wfn:
    def __init__(self, specs):
        # specs = [PAO, Pseudopotential, Orbits], same as Atom_specs[label] in LCAO
        self.PAO=specs[0]
        self.PP=specs[1]
        self.Orbits_str=specs[2] # like "s2p2d1"
        self.length=0 # data length
        self.r=None # numpy[length], values of r
        self.Wfn=[] # list of [r][AO or PAO], the order follws Orbits

        # decomposition of Orbits_str
        self.Orbits=[] # like "s0", "s1", "p0", "p1", "d0", "d1"
        i=0
        while i<len(self.Orbits_str):
            orbit_label=self.Orbits_str[i]
            mul=int(self.Orbits_str[i+1])
            for j in range(0, mul):
                self.Orbits.append(("{0:s}{1:1d}").format(orbit_label, j))
            i+=2

        groupName=self.PAO+"&"+self.PP
        # load PAO and AO
        with h5py.File(Config.PAO_and_AO, "r") as f:
            g=f[groupName]
            self.length=g.attrs["length"]
            self.r=np.array(g.attrs["r"])
            for orbit in self.Orbits:
                self.Wfn.append(np.array(g[orbit]))
            
class PSF:
    def __init__(self, LCAO, Wfns):
        self.LCAO=LCAO
        self.Wfns=Wfns
        self.initialStates_i=0 # 0->PAO, 1->AO
        self.finalStates_i=0 # 0->Plane wave, 1->Calculated
        self.finalStates_step=0.0
        self.Y_coeff=[0, 0, 0]
        self.configured=False
        self.Gaunt=np.zeros((5,11,4,9)) # [lp][mp+lp][l][m+l]
        for lp in range(0, 5):
            for mp in range(-lp, lp+1):
                for l in range(0, 4):
                    for m in range(-l, l+1):
                        self.Gaunt[lp][mp+lp][l][m+l]=pt.Gaunt(lp, mp, l, m)

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
    
    def calc(self, ik, ib):
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

        k_au=self.LCAO.Kpath_au[ik]
        k_length=math.sqrt(np.inner(k_au, k_au))
        ret=0.0
        # prepare spherical harmonics
        Ylm_k=np.zeros((5,11), dtype=complex)
        for lp in range(0, 5):
            for mp in range(-lp, lp+1):
                Ylm_k[lp][mp+lp]=pt.sphericalHarmonics(lp, mp, k_au)
        # for debug
        # if ik==14825:
        #     print(k_au)
        #     print(Ylm_k)
        
        for ia in range(0, self.LCAO.numAtoms):
            atom_label=self.LCAO.Atoms[ia]
            wfnObj=self.Wfns[atom_label]
            kt=np.inner(k_au, self.LCAO.Atom_au[ia])
            atom_phase=math.cos(kt)-math.sin(kt)*1j
            # prepare final states, atom dependent in case of calculated final states
            finalStates=np.zeros((5, wfnObj.length))
            for lp in range(0, 5):
                self.calcFinalState(finalStates[lp], lp, k_length, wfnObj.r)
                    
            for io, orbit_label in enumerate(wfnObj.Orbits):
                l=0
                if orbit_label[0]=="s":
                    pass
                elif orbit_label[0]=="p":
                    l=1
                elif orbit_label[0]=="d":
                    l=2
                elif orbit_label[0]=="f":
                    l=3
                else:
                    print("Error: invalid azimuthal quantum label")
                    return 0.0
                iL=0
                iL_found=False
                for il, LCAO_label in enumerate(self.LCAO.LCAO_labels[ia]):
                    if LCAO_label==orbit_label:
                        iL=il
                        iL_found=True
                        break
                if iL_found==False:
                    return 0.0
                LCAO_iL=self.LCAO.LCAO[ia][iL][ik][ib]
                LCAO_iL_conv=[]
                # conversion of LCAO coefficients
                if l==0:
                    LCAO_iL_conv.append(LCAO_iL[0][0]+LCAO_iL[0][1]*1j)
                elif l==1:
                    p1=LCAO_iL[0][0]+LCAO_iL[0][1]*1j
                    p2=LCAO_iL[1][0]+LCAO_iL[1][1]*1j
                    p3=LCAO_iL[2][0]+LCAO_iL[2][1]*1j
                    LCAO_iL_conv=pt.convertLCAO_p(p1, p2, p3)
                elif l==2:
                    d1=LCAO_iL[0][0]+LCAO_iL[0][1]*1j
                    d2=LCAO_iL[1][0]+LCAO_iL[1][1]*1j
                    d3=LCAO_iL[2][0]+LCAO_iL[2][1]*1j
                    d4=LCAO_iL[3][0]+LCAO_iL[3][1]*1j
                    d5=LCAO_iL[4][0]+LCAO_iL[4][1]*1j
                    LCAO_iL_conv=pt.convertLCAO_d(d1, d2, d3, d4, d5)
                elif l==3:
                    f1=LCAO_iL[0][0]+LCAO_iL[0][1]*1j
                    f2=LCAO_iL[1][0]+LCAO_iL[1][1]*1j
                    f3=LCAO_iL[2][0]+LCAO_iL[2][1]*1j
                    f4=LCAO_iL[3][0]+LCAO_iL[3][1]*1j
                    f5=LCAO_iL[4][0]+LCAO_iL[4][1]*1j
                    f6=LCAO_iL[5][0]+LCAO_iL[5][1]*1j
                    f7=LCAO_iL[6][0]+LCAO_iL[6][1]*1j
                    LCAO_iL_conv=pt.convertLCAO_f(f1, f2, f3, f4, f5, f6, f7)
                    
                wfn_initial=wfnObj.Wfn[io][:][self.initialStates_i]
                for dl in dls:
                    lp=l+dl
                    if lp<0:
                        continue
                    radialPart=pt.radialIntegral(wfn_initial, finalStates[lp], wfnObj.r)
                    for mpl in range(0, 2*l+1):
                        m=mpl-l
                        for jp1 in range(0, 3):
                            j=jp1-1
                            if m+j<-lp or m+j>lp:
                                continue
                            gaunt_coeff=self.Gaunt[lp][m+j+lp][l][m+l]
                            Ylpmp=Ylm_k[lp][m+j+lp]
                            ret+=((-1j)**lp)*Ylpmp*atom_phase*self.Y_coeff[jp1]*LCAO_iL_conv[mpl]*gaunt_coeff*radialPart
                        
                
                    
        return ret.real**2+ret.imag**2
