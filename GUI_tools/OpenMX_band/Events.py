# Events

from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import numpy as np
import re
import h5py
import math
from scipy.stats import norm

import Config

def openFile(win):
    global Band
    global BandUp
    global BandDn
    global BandCell
    global RecCell
    global Dimension
    global Spin
    global Spin_i
    global EF_Eh
    global Xrange
    global Yrange
    global Curved

    currentFile=win.filePath.text()
    selectedFile, _filter=QtGui.QFileDialog.getOpenFileName(caption="Open file", directory=currentFile)
    if selectedFile!="":
        # valid
        win.filePath.setText(selectedFile)

        with h5py.File(selectedFile, "r") as f:
            unit=str(f["Input"]["UnitCells"].attrs["Unit"])
            win.plot.setLabel(axis="bottom", text=("Wavevector ({0:s}^-1)").format(unit))
            BandCell=np.array(f["Input"]["UnitCells"]["Bands"])
            Dimension=int(f["Input"]["Kpath"].attrs["Dimension"])
            Spin=str(f["Output"].attrs["Spin"])
            EF_Eh=float(f["Output"].attrs["EF_Eh"])
            Xrange=np.array(f["Input"]["Kpath"].attrs["Xrange"])
            if Dimension==2:
                Yrange=np.array(f["Input"]["Kpath"].attrs["Yrange"])

            if Spin.lower()=="off":
                Spin_i=0
                Band=np.array(f["Output"]["Band"])
            elif Spin.lower()=="on":
                Spin_i=1
                BandUp=np.array(f["Output"]["BandUp"])
                BandDn=np.array(f["Output"]["BandDn"])
            elif Spin.lower()=="nc":
                Spin_i=2
                Band=np.array(f["Output"]["Band"])
            Curved=f["Input"]["Kpath"].attrs["Curved"]


        print(("Dimension: {0:d}").format(Dimension))
        print(("Spin: {0:s}").format(Spin))
        print(("EF: {0:.3f} eV").format(Config.Eh*EF_Eh))
        print("X range:")
        print(Xrange)
        if Dimension==2:
            print("Y range:")
            print(Yrange)
        print(("Curved: {0:s}").format(str(Curved)))

        # calculate the reciprocal unit cell
        RecCell=np.zeros((3, 3))
        op=np.cross(BandCell[1], BandCell[2])
        det=np.inner(BandCell[0], op)
        
        ## 0
        RecCell[0]=2*math.pi*op/det
        
        ## 1
        op=np.cross(BandCell[2], BandCell[0])
        RecCell[1]=2*math.pi*op/det

        ## 2
        op=np.cross(BandCell[0], BandCell[1])
        RecCell[2]=2*math.pi*op/det

        print("Reciprocal unit cell:")
        print(RecCell)

def appendDispersion(i, n, EMin, EPixel, tailProfile):
    global Dispersion
    global Band
    global BandUp
    global BandDn
    global EF_Eh
    global Spin_i

    ret=0
    Esize=Dispersion.shape[1]
    tailSize=tailProfile.shape[0]
    if Spin_i==1:
        continueFlag=False
        
        eigen=(BandUp[i][n]-EF_Eh)*Config.Eh
        eigen_index=round((eigen-EMin)/EPixel)
        if eigen_index-tailSize>=Esize:
            pass
        else:
            continueFlag=True
            for j in range(-tailSize+1, tailSize):
                if eigen_index+j>=0 and eigen_index+j<Esize:
                    Dispersion[i][eigen_index+j]+=tailProfile[abs(j)]
        
        eigen=(BandDn[i][n]-EF_Eh)*Config.Eh
        eigen_index=round((eigen-EMin)/EPixel)
        if eigen_index-tailSize>=Esize:
            pass
        else:
            continueFlag=True
            for j in range(-tailSize+1, tailSize):
                if eigen_index+j>=0 and eigen_index+j<Esize:
                    Dispersion[i][eigen_index+j]+=tailProfile[abs(j)]

        return continueFlag

    else:
        eigen=(Band[i][n]-EF_Eh)*Config.Eh
        eigen_index=round((eigen-EMin)/EPixel)
        if eigen_index-tailSize>=Esize:
            return False
        for j in range(-tailSize+1, tailSize):
            if eigen_index+j>=0 and eigen_index+j<Esize:
                Dispersion[i][eigen_index+j]+=tailProfile[abs(j)]
        return True
        
        

def plot(win):
    global Dispersion
    global Band
    global BandUp
    global BandDn
    global RecCell
    global Spin_i
    global EF_Eh
    global Dimension
    global Xrange
    global Yrange

    EMin=float(win.EMin.text())
    EMax=float(win.EMax.text())
    dE=float(win.dE.text())
    EPixel=float(win.EPixel.text())

    tailIndex=math.floor(dE*Config.sigma_max/EPixel)
    tailProfile=np.zeros((tailIndex+1,))
    for i in range(0, tailIndex+1):
        tailProfile[i]=norm.pdf(i*EPixel, loc=0, scale=dE)
    print(tailProfile)

    numPnts_E=math.ceil((EMax-EMin)/EPixel+1)
    if numPnts_E<0:
        print("Energy range error")
        return

    print(("{0:d} points along the energy").format(numPnts_E))

    numPnts_k=0
    numBands=0
    if Spin_i==1:
        numPnts_k=BandUp.shape[0]
        numBands=BandUp.shape[1]
    else:
        numPnts_k=Band.shape[0]
        numBands=Band.shape[1]

    print(("{0:d} points along the wavevector(s)").format(numPnts_k))

    if Dimension==1:
        Dispersion=np.zeros((numPnts_k, numPnts_E))
        for i in range(0, numPnts_k):
            for n in range(0, numBands):
                if appendDispersion(i, n, EMin, EPixel, tailProfile):
                    continue
                else:
                    break

        xLength=math.sqrt(np.inner(RecCell[0], RecCell[0]))
        dx_length=xLength*(Xrange[1]-Xrange[0])/(numPnts_k-1)
        tr=QtGui.QTransform()
        tr.translate(xLength*Xrange[0],EMin)
        tr.scale(dx_length, EPixel)
        win.img.setTransform(tr)


    win.img.setImage(Dispersion)
            




