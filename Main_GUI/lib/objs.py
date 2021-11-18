# calcPSF GUI
# class LCAO

import numpy as np
import h5py
import math
    
class LCAO:
    def __init__(self):
        self.filePath="" # file path to the OpenMX output (.hdf5)
        # /Input/UnitCells
        self.unit="" # unit of the unit cells (Ang or AU)
        self.BandCell=None # real space unit cell for the band dispersion
        self.AtomCell=None # real space unit cell for the crystal structure
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
        self.Kpath=None # list of the coordinates
        # calculated from /Input/Kpath
        self.RecCell=None # reciprocal unit cell
        self.Xvector=None # X vector in the reciprocal coordinate
        self.Xlength=0 # length of the X vector
        self.Yvector=None # Y vector in the reciprocal coordinate
        self.Ylength=0 # length of the Y vector
        self.numPnts_k=0 # total number of points
        # /Output
        self.Spin="" # spin polarization: off, on, nc
        self.Spin_i=0 # index corresponding to the spin polarization: 0, 1, 2 respectively
        self.EF_Eh=0 # The Fermi energy in units of the Hartree energy (27.2 eV)
        self.Band=None # band dispersion, used in off or nc
        self.BandUp=None # band dispersion, used in on
        self.BandDn=None # band dispersion, used in on
        self.numBands=0 # number of bands in the data
        self.Atoms=[] # labels for atoms
        self.LCAO=[] # LCAO coefficients [atom index][orbital index](numpy[k][band][l][Re, Im])
        self.LCAO_labels=[] # LCAO labels, related to orbital index
        

    def open(self, filePath):
        self.filePath=filePath

        with h5py.File(self.filePath, "r") as f:
            # /Input/UnitCells
            self.unit=str(f["Input"]["UnitCells"].attrs["Unit"])
            self.BandCell=np.array(f["Input"]["UnitCells"]["Bands"])
            self.AtomCell=np.array(f["Input"]["UnitCells"]["Atoms"])
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
            for atom in f["Output"]["LCAO"].keys():
                self.Atoms.append(str(atom))
            for atom in self.Atoms:
                lcao_tmp=[]
                lcao_label_tmp=[]
                for lcao_key in f["Output"]["LCAO"][atom].keys():
                    lcao_label_tmp.append(lcao_key)
                    lcao_tmp.append(np.array(f["Output"]["LCAO"][atom][lcao_key]))
                self.LCAO.append(lcao_tmp)
                self.LCAO_labels.append(lcao_label_tmp)
