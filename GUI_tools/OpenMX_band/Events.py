# Events

from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import numpy as np
import re
import h5py
import math
from scipy.stats import norm

import pyqtgraph.opengl as gl

import Config


BandCell=None
RecCell=None
Band=None
BandUp=None
BandDn=None
Dimension=0
Spin=""
Spin_i=0
EF_Eh=0
Dispersion=None
Atoms=[]
LCAO=[]
LCAO_labels=[]
numPnts_kx=0
numPnts_ky=0
Xrange=None
Yrange=None
Xvector=None
Yvector=None
Xlength=0
Ylength=0

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
    global Xvector
    global Yvector
    global Xlength
    global Ylength
    global Curved
    global Atoms
    global LCAO
    global LCAOUp
    global LCAODn
    global LCAO_labels
    global numPnts_kx
    global numPnts_ky

    currentFile=win.filePath.text()
    selectedFile, _filter=QtGui.QFileDialog.getOpenFileName(caption="Open file", directory=currentFile)
    if selectedFile!="":
        # valid
        win.filePath.setText(selectedFile)

        with h5py.File(selectedFile, "r") as f:
            unit=str(f["Input"]["UnitCells"].attrs["Unit"])
            win.plot.setLabel(axis="bottom", text=("Wavevector ({0:s}^-1)").format(unit))
            win.plotEx.setLabel(axis="bottom", text=("Wavevector ({0:s}^-1)").format(unit))
            win.plotEy.setLabel(axis="left", text=("Wavevector ({0:s}^-1)").format(unit))
            win.plotxy.setLabel(axis="left", text=("Wavevector ({0:s}^-1)").format(unit))
            win.plotxy.setLabel(axis="bottom", text=("Wavevector ({0:s}^-1)").format(unit))
            BandCell=np.array(f["Input"]["UnitCells"]["Bands"])
            Dimension=int(f["Input"]["Kpath"].attrs["Dimension"])
            Spin=str(f["Output"].attrs["Spin"])
            EF_Eh=float(f["Output"].attrs["EF_Eh"])
            Xrange=np.array(f["Input"]["Kpath"].attrs["Xrange"])
            Xvector_frac=np.array(f["Input"]["Kpath"].attrs["Xvector"])
            numPnts_kx=int(f["Input"]["Kpath"].attrs["Xcount"])
            numPnts_ky=1
            if Dimension==2:
                Yrange=np.array(f["Input"]["Kpath"].attrs["Yrange"])
                Yvector_frac=np.array(f["Input"]["Kpath"].attrs["Yvector"])
                numPnts_ky=int(f["Input"]["Kpath"].attrs["Ycount"])
            numPnts_k=numPnts_kx*numPnts_ky
            
            Curved=f["Input"]["Kpath"].attrs["Curved"]

            if Spin.lower()=="off":
                Spin_i=0
                Band=np.array(f["Output"]["Band"])
                win.UpButton.setCheckable(False)
                win.DnButton.setCheckable(False)
                numBands=Band.shape[1]
            elif Spin.lower()=="on":
                Spin_i=1
                BandUp=np.array(f["Output"]["BandUp"])
                BandDn=np.array(f["Output"]["BandDn"])
                win.UpButton.setCheckable(True)
                win.DnButton.setCheckable(True)
                win.UpButton.setChecked(True)
                numPnts_k=BandUp.shape[0]
                numBands=BandUp.shape[1]
            elif Spin.lower()=="nc":
                Spin_i=2
                Band=np.array(f["Output"]["Band"])
                win.UpButton.setCheckable(False)
                win.DnButton.setCheckable(False)
                numPnts_k=Band.shape[0]
                numBands=Band.shape[1]

            
            win.kxIndex.setMaximum(numPnts_kx-1)
            win.kyIndex.setMaximum(numPnts_ky-1)
            win.bIndex.setMaximum(numBands-1)

            Atoms=[]
            LCAO=[]
            LCAO_labels=[]
            for atom in f["Output"]["LCAO"].keys():
                Atoms.append(str(atom))
            win.Atom.clear()
            for atom in Atoms:
                win.Atom.addItem(atom)
                lcao_tmp=[]
                lcao_label_tmp=[]
                for lcao_key in f["Output"]["LCAO"][atom].keys():
                    lcao_label_tmp.append(lcao_key)
                    lcao_tmp.append(np.array(f["Output"]["LCAO"][atom][lcao_key]))
                LCAO.append(lcao_tmp)
                LCAO_labels.append(lcao_label_tmp)

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

        Xvector=np.zeros((3,))
        for i in range(0, 3):
            Xvector+=Xvector_frac[i]*RecCell[i]
        Xlength=math.sqrt(np.inner(Xvector, Xvector))
        print("X vector:")
        print(Xvector)
        if Dimension==2:
            Yvector=np.zeros((3,))
            for i in range(0, 3):
                Yvector+=Yvector_frac[i]*RecCell[i]
            Ylength=math.sqrt(np.inner(Yvector, Yvector))
            print("Y vector:")
            print(Yvector)





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
        

def appendDispersion3(ix, iy, n, EMin, EPixel, tailProfile):
    global Dispersion
    global Band
    global BandUp
    global BandDn
    global EF_Eh
    global Spin_i
    global numPnts_kx

    ret=0
    Esize=Dispersion.shape[2]
    tailSize=tailProfile.shape[0]

    i=ix+iy*numPnts_kx
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
                    Dispersion[ix][iy][eigen_index+j]+=tailProfile[abs(j)]
        
        eigen=(BandDn[i][n]-EF_Eh)*Config.Eh
        eigen_index=round((eigen-EMin)/EPixel)
        if eigen_index-tailSize>=Esize:
            pass
        else:
            continueFlag=True
            for j in range(-tailSize+1, tailSize):
                if eigen_index+j>=0 and eigen_index+j<Esize:
                    Dispersion[ix][iy][eigen_index+j]+=tailProfile[abs(j)]

        return continueFlag

    else:
        eigen=(Band[i][n]-EF_Eh)*Config.Eh
        eigen_index=round((eigen-EMin)/EPixel)
        if eigen_index-tailSize>=Esize:
            return False
        for j in range(-tailSize+1, tailSize):
            if eigen_index+j>=0 and eigen_index+j<Esize:
                Dispersion[ix][iy][eigen_index+j]+=tailProfile[abs(j)]
        return True
                

def plot(win):
    # for 2D plot
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
    global Xlength
    global Ylength

    EMin=float(win.EMin.text())
    EMax=float(win.EMax.text())
    dE=float(win.dE.text())
    EPixel=float(win.EPixel.text())

    tailIndex=math.floor(dE*Config.sigma_max/EPixel)
    tailProfile=np.zeros((tailIndex+1,))
    for i in range(0, tailIndex+1):
        tailProfile[i]=norm.pdf(i*EPixel, loc=0, scale=dE)

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

    print(("{0:d} points along the kx").format(numPnts_k))

    if Dimension==1:
        Dispersion=np.zeros((numPnts_k, numPnts_E))
        for i in range(0, numPnts_k):
            for n in range(0, numBands):
                if appendDispersion(i, n, EMin, EPixel, tailProfile):
                    continue
                else:
                    break

        dx_length=Xlength*(Xrange[1]-Xrange[0])/(numPnts_k-1)
        tr=QtGui.QTransform()
        tr.translate(Xlength*Xrange[0]-dx_length/2,EMin-EPixel/2)
        tr.scale(dx_length, EPixel)
        win.img.setTransform(tr)

    win.img.setImage(Dispersion)
    drawCursor(win)


def makeDispersion3(win):
    # for 3D plot
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
    global Xlength
    global Ylength
    global numPnts_kx
    global numPnts_ky

    EMin=float(win.EMin.text())
    EMax=float(win.EMax.text())
    dE=float(win.dE.text())
    EPixel=float(win.EPixel.text())

    tailIndex=math.floor(dE*Config.sigma_max/EPixel)
    tailProfile=np.zeros((tailIndex+1,))
    for i in range(0, tailIndex+1):
        tailProfile[i]=norm.pdf(i*EPixel, loc=0, scale=dE)

    numPnts_E=math.ceil((EMax-EMin)/EPixel+1)
    if numPnts_E<0:
        print("Energy range error")
        return

    print(("{0:d} points along the energy").format(numPnts_E))

    win.eIndex.setMaximum(numPnts_E-1)

    numBands=0
    if Spin_i==1:
        numBands=BandUp.shape[1]
    else:
        numBands=Band.shape[1]

    print(("{0:d} points along the kx").format(numPnts_kx))
    print(("{0:d} points along the ky").format(numPnts_ky))

    if Dimension==2:
        Dispersion=np.zeros((numPnts_kx, numPnts_ky, numPnts_E))
        for i in range(0, numPnts_kx):
            for j in range(0, numPnts_ky):
                for n in range(0, numBands):
                    if appendDispersion3(i, j, n, EMin, EPixel, tailProfile):
                        continue
                    else:
                        break

        dx_length=Xlength*(Xrange[1]-Xrange[0])/(numPnts_kx-1)
        dy_length=Ylength*(Yrange[1]-Yrange[0])/(numPnts_ky-1)

        tr_x=QtGui.QTransform()
        tr_x.translate(Xlength*Xrange[0]-dx_length/2,EMin-EPixel/2)
        tr_x.scale(dx_length, EPixel)
        win.imgEx.setTransform(tr_x)

        tr_y=QtGui.QTransform()
        tr_y.translate(EMin-EPixel/2,Ylength*Yrange[0]-dy_length/2)
        tr_y.rotate(-90)
        tr_y.scale(-dy_length, EPixel)
        win.imgEy.setTransform(tr_y)

        tr_E=QtGui.QTransform()
        tr_E.translate(Xlength*Xrange[0]-dx_length/2,Ylength*Yrange[0]-dy_length/2)
        tr_E.scale(dx_length, dy_length)
        win.imgxy.setTransform(tr_E)
    else:
        return

    Maxpoint=Dispersion.max()

    win.Cube=np.zeros((numPnts_kx, numPnts_ky, numPnts_E, 4))
    win.Cube[:,:,:,0]=255
    win.Cube[:,:,:,1]=255
    win.Cube[:,:,:,2]=255
    win.Cube[:,:,:,3]=Dispersion/Maxpoint*100
    
    win.bandCube=gl.GLVolumeItem(win.Cube)
    win.bandCube.scale(dx_length, dx_length, EPixel)    
    win.bandCube.translate(Xlength*Xrange[0]-dx_length/2,Ylength*Yrange[0]-dy_length/2,0)
    win.plot3D.clear()
    win.plot3D.addItem(win.bandCube)
    plot3(win)

def plot3(win):

    global Dispersion
    Maxpoint=Dispersion.max()

    kx=win.kxIndex.value()
    ky=win.kyIndex.value()
    ei=win.eIndex.value()
    win.imgEx.setImage(Dispersion[:,ky,:])
    win.imgEy.setImage(Dispersion[kx,:,:])
    win.imgxy.setImage(Dispersion[:,:,ei])

    win.Cube[:,:,:,0]=255
    win.Cube[:,:,:,1]=255
    win.Cube[:,:,:,2]=255
    win.Cube[:,:,:,3]=Dispersion/Maxpoint*100

    win.Cube[kx,:,:,0]=Config.pen1[0]
    win.Cube[kx,:,:,1]=Config.pen1[1]
    win.Cube[kx,:,:,2]=Config.pen1[2]
    win.Cube[kx,:,:,3]=Config.gridAlpha

    win.Cube[:,ky,:,0]=Config.pen2[0]
    win.Cube[:,ky,:,1]=Config.pen2[1]
    win.Cube[:,ky,:,2]=Config.pen2[2]
    win.Cube[:,ky,:,3]=Config.gridAlpha

    win.Cube[:,:,ei,0]=Config.pen3[0]
    win.Cube[:,:,ei,1]=Config.pen3[1]
    win.Cube[:,:,ei,2]=Config.pen3[2]
    win.Cube[:,:,ei,3]=Config.gridAlpha

    win.bandCube.setData(win.Cube)

    ExMax=Dispersion[:,ky,:].max()
    EyMax=Dispersion[kx,:,:].max()
    xyMax=Dispersion[:,:,ei].max()
    Maxpoint=max(ExMax, EyMax, xyMax)
    win.imgEx.setLevels([0,Maxpoint])
    win.imgEy.setLevels([0,Maxpoint])
    win.imgxy.setLevels([0,Maxpoint])

    drawCursor3(win)
            
def drawCursor(win):
    global BandUp
    global BandDn
    global Band
    global Spin_i
    global RecCell
    global Xrange
    global Xlength
    global EF_Eh
    
    k=win.kxIndex.value()
    b=win.bIndex.value()
    
    numPnts_k=0
    numBands=0
    UseUp=False
    if Spin_i==1:
        numPnts_k=BandUp.shape[0]
        numBands=BandUp.shape[1]
        if win.UpButton.isChecked():
            UseUp=True
        elif win.DnButton.isChecked():
            UseUp=False
        else:
            print("None of Up and Dn is checked")
            return
    else:
        numPnts_k=Band.shape[0]
        numBands=Band.shape[1]

    dx_length=Xlength*(Xrange[1]-Xrange[0])/(numPnts_k-1)
    
    if 0<=k and k<numPnts_k and 0<=b and b<numBands:
        k_value=Xlength*Xrange[0]+dx_length*k
        b_value=0
        if Spin_i==1:
            if UseUp:
                b_value=BandUp[k][b]
            else:
                b_value=BandDn[k][b]
        else:
            b_value=Band[k][b]

        b_value=(b_value-EF_Eh)*Config.Eh
        win.vLine.setPos(k_value)
        win.hLine.setPos(b_value)
        win.kxValue.setText(("({0:.3f})").format(k_value))
        win.bValue.setText(("({0:.3f})").format(b_value))
        
    else:
        print("Index error")
        return

    makeLCAOTable(win)

            
def drawCursor3(win):
    global BandUp
    global BandDn
    global Band
    global Spin_i
    global RecCell
    global Xrange
    global Xlength
    global Yrange
    global Ylength
    global EF_Eh
    global numPnts_kx
    global numPnts_ky

    EMin=float(win.EMin.text())
    EPixel=float(win.EPixel.text())
    
    kx=win.kxIndex.value()
    ky=win.kyIndex.value()
    b=win.bIndex.value()
    ei=win.eIndex.value()
    numPnts_k=numPnts_kx*numPnts_ky
    numBands=0
    UseUp=False
    if Spin_i==1:
        numBands=BandUp.shape[1]
        if win.UpButton.isChecked():
            UseUp=True
        elif win.DnButton.isChecked():
            UseUp=False
        else:
            print("None of Up and Dn is checked")
            return
    else:
        numBands=Band.shape[1]

    dx_length=Xlength*(Xrange[1]-Xrange[0])/(numPnts_kx-1)
    dy_length=Ylength*(Yrange[1]-Yrange[0])/(numPnts_ky-1)
    k=kx+ky*numPnts_kx

    if 0<=k and k<numPnts_k and 0<=b and b<numBands:
        kx_value=Xlength*Xrange[0]+dx_length*kx
        ky_value=Ylength*Yrange[0]+dy_length*ky
        e_value=EMin+EPixel*ei
        b_value=0
        if Spin_i==1:
            if UseUp:
                b_value=BandUp[k][b]
            else:
                b_value=BandDn[k][b]
        else:
            b_value=Band[k][b]

        b_value=(b_value-EF_Eh)*Config.Eh
        win.vLineEx.setPos(kx_value)
        win.hLineEx.setPos(e_value)
        win.vLineEy.setPos(e_value)
        win.hLineEy.setPos(ky_value)
        win.vLinexy.setPos(kx_value)
        win.hLinexy.setPos(ky_value)
        win.bLineEx.setPos(b_value)
        win.bLineEy.setPos(b_value)
        win.kxValue.setText(("({0:.3f})").format(kx_value))
        win.kyValue.setText(("({0:.3f})").format(ky_value))
        win.bValue.setText(("({0:.3f})").format(b_value))
        win.eValue.setText(("({0:.3f})").format(e_value))
        
    else:
        print("Index error")
        return

    makeLCAOTable(win)

            
def makeLCAOTable(win):
    global BandUp
    global BandDn
    global Band
    global Spin_i
    global RecCell
    global Xrange
    global EF_Eh
    global Atoms
    global LCAO_labels
    global LCAO
    global numPnts_kx
    global numPnts_ky
    global Dimension
    
    kx=win.kxIndex.value()
    ky=win.kyIndex.value()
    k=0
    if Dimension==1:
        k=kx
    elif Dimension==2:
        k=kx+ky*numPnts_kx

    b=win.bIndex.value()
    at=win.Atom.currentIndex()
    
    numPnts_k=numPnts_kx*numPnts_ky
    numBands=0
    UseUp=False
    if Spin_i==1:
        numBands=BandUp.shape[1]
        if win.UpButton.isChecked():
            UseUp=True
        elif win.DnButton.isChecked():
            UseUp=False
        else:
            print("None of Up and Dn is checked")
            return
    else:
        numBands=Band.shape[1]

    if len(LCAO_labels)!=len(Atoms):
        return

    win.LCAOTable.setRowCount(0)
    if Spin_i==2:
        win.LCAOTable.setColumnCount(4)
        win.LCAOTable.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem("Up (raw)"))
        win.LCAOTable.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem("Up (calc)"))
        win.LCAOTable.setHorizontalHeaderItem(2, QtGui.QTableWidgetItem("Dn (raw)"))
        win.LCAOTable.setHorizontalHeaderItem(3, QtGui.QTableWidgetItem("Dn (calc)"))
    else:
        win.LCAOTable.setColumnCount(2)
        if Spin_i==1:
            win.LCAOTable.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem(("{0:s} (raw)").format("Up" if UseUp else "Dn")))
            win.LCAOTable.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem(("{0:s} (calc)").format("Up" if UseUp else "Dn")))
        else:
            win.LCAOTable.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem("raw"))
            win.LCAOTable.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem("calc"))

                                              

    currentRow=0
    if 0<=k and k<numPnts_k and 0<=b and b<numBands and 0<=at and at<len(Atoms):
        orbitLabels=["s", "p", "d", "f"]
        for orbitLabel in orbitLabels:
            mul_index=0
            while(True):
                LCAO_label=("{0:s}{1:1d}").format(orbitLabel, mul_index)
                if Spin_i==1:
                    if UseUp:
                        LCAO_label+="Up"
                    else:
                        LCAO_label+="Dn"
                    
                LCAO_found=False
                for i, LCAO_label_ref in enumerate(LCAO_labels[at]):
                    if LCAO_label==LCAO_label_ref:
                        LCAO_found=True

                        LCAO_disp=LCAO[at][i][k][b]
                        numSpin=1
                        if Spin_i==2:
                            numSpin=2
                        if orbitLabel=="s":
                            currentRow+=1
                        elif orbitLabel=="p":
                            currentRow+=3
                        elif orbitLabel=="d":
                            currentRow+=5
                        elif orbitLabel=="f":
                            currentRow+=7
                        win.LCAOTable.setRowCount(currentRow)
                        for s in range(0, numSpin):
                            s2=s*2
                            if orbitLabel=="s":
                                # s orbital: nothing to calculate
                                item1=QtGui.QTableWidgetItem(("s ({0:.2f}, {1:.2f})").format(LCAO_disp[0][0+s2], LCAO_disp[0][1+s2]))
                                win.LCAOTable.setItem(currentRow-1, 0+s2, item1)
                                item2=QtGui.QTableWidgetItem(("s ({0:.2f}, {1:.2f})").format(LCAO_disp[0][0+s2], LCAO_disp[0][1+s2]))
                                win.LCAOTable.setItem(currentRow-1, 1+s2, item2)
                                head=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                win.LCAOTable.setVerticalHeaderItem(currentRow-1, head)
                                
                            if orbitLabel=="p":
                                # p orbital: px (cos P), py (sin P), pz (1)
                                px=LCAO_disp[0][0+s2]+LCAO_disp[0][1+s2]*1j
                                py=LCAO_disp[1][0+s2]+LCAO_disp[1][1+s2]*1j
                                pz=LCAO_disp[2][0+s2]+LCAO_disp[2][1+s2]*1j
                                pm1=(px-1j*py)/math.sqrt(2)
                                pp0=pz
                                pp1=(px+1j*py)/math.sqrt(2)
                                # table
                                item1=QtGui.QTableWidgetItem(("px ({0:.2f}, {1:.2f})").format(LCAO_disp[0][0+s2], LCAO_disp[0][1+s2]))
                                item2=QtGui.QTableWidgetItem(("py ({0:.2f}, {1:.2f})").format(LCAO_disp[1][0+s2], LCAO_disp[1][1+s2]))
                                item3=QtGui.QTableWidgetItem(("pz ({0:.2f}, {1:.2f})").format(LCAO_disp[2][0+s2], LCAO_disp[2][1+s2]))
                                win.LCAOTable.setItem(currentRow-3, 0+s2, item1)
                                win.LCAOTable.setItem(currentRow-2, 0+s2, item2)
                                win.LCAOTable.setItem(currentRow-1, 0+s2, item3)
                                item1=QtGui.QTableWidgetItem(("p(-1) ({0:.2f}, {1:.2f})").format(pm1.real, pm1.imag))
                                item2=QtGui.QTableWidgetItem(("p(+0) ({0:.2f}, {1:.2f})").format(pp0.real, pp0.imag))
                                item3=QtGui.QTableWidgetItem(("p(+1) ({0:.2f}, {1:.2f})").format(pp1.real, pp1.imag))
                                win.LCAOTable.setItem(currentRow-3, 1+s2, item1)
                                win.LCAOTable.setItem(currentRow-2, 1+s2, item2)
                                win.LCAOTable.setItem(currentRow-1, 1+s2, item3)
                                head1=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head2=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head3=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                win.LCAOTable.setVerticalHeaderItem(currentRow-3, head1)
                                win.LCAOTable.setVerticalHeaderItem(currentRow-2, head2)
                                win.LCAOTable.setVerticalHeaderItem(currentRow-1, head3)
                                
                            if orbitLabel=="d":
                                # d orbital: d3z^2-r^2 (1), dx^2-y^2 (cos 2P), xy (sin 2P), xz (cos P), yz (sin P)
                                d3z2r2=LCAO_disp[0][0+s2]+LCAO_disp[0][1+s2]*1j
                                dx2y2 =LCAO_disp[1][0+s2]+LCAO_disp[1][1+s2]*1j
                                dxy   =LCAO_disp[2][0+s2]+LCAO_disp[2][1+s2]*1j
                                dxz   =LCAO_disp[3][0+s2]+LCAO_disp[3][1+s2]*1j
                                dyz   =LCAO_disp[4][0+s2]+LCAO_disp[4][1+s2]*1j
                                dm2=(dx2y2-1j*dxy)/math.sqrt(2)
                                dm1=(dxz-1j*dyz)/math.sqrt(2)
                                dp0=d3z2r2
                                dp1=(dxz+1j*dyz)/math.sqrt(2)
                                dp2=(dx2y2+1j*dxy)/math.sqrt(2)
                                # table
                                item1=QtGui.QTableWidgetItem(("d3z2r2 ({0:.2f}, {1:.2f})").format(LCAO_disp[0][0+s2], LCAO_disp[0][1+s2]))
                                item2=QtGui.QTableWidgetItem(("dx2y2  ({0:.2f}, {1:.2f})").format(LCAO_disp[1][0+s2], LCAO_disp[1][1+s2]))
                                item3=QtGui.QTableWidgetItem(("dxy    ({0:.2f}, {1:.2f})").format(LCAO_disp[2][0+s2], LCAO_disp[2][1+s2]))
                                item4=QtGui.QTableWidgetItem(("dxz    ({0:.2f}, {1:.2f})").format(LCAO_disp[3][0+s2], LCAO_disp[3][1+s2]))
                                item5=QtGui.QTableWidgetItem(("dyz    ({0:.2f}, {1:.2f})").format(LCAO_disp[4][0+s2], LCAO_disp[4][1+s2]))
                                win.LCAOTable.setItem(currentRow-5, 0+s2, item1)
                                win.LCAOTable.setItem(currentRow-4, 0+s2, item2)
                                win.LCAOTable.setItem(currentRow-3, 0+s2, item3)
                                win.LCAOTable.setItem(currentRow-2, 0+s2, item4)
                                win.LCAOTable.setItem(currentRow-1, 0+s2, item5)
                                item1=QtGui.QTableWidgetItem(("d(-2) ({0:.2f}, {1:.2f})").format(dm2.real, dm2.imag))
                                item2=QtGui.QTableWidgetItem(("d(-1) ({0:.2f}, {1:.2f})").format(dm1.real, dm1.imag))
                                item3=QtGui.QTableWidgetItem(("d(+0) ({0:.2f}, {1:.2f})").format(dp0.real, dp0.imag))
                                item4=QtGui.QTableWidgetItem(("d(+1) ({0:.2f}, {1:.2f})").format(dp1.real, dp1.imag))
                                item5=QtGui.QTableWidgetItem(("d(+2) ({0:.2f}, {1:.2f})").format(dp2.real, dp2.imag))
                                win.LCAOTable.setItem(currentRow-5, 1+s2, item1)
                                win.LCAOTable.setItem(currentRow-4, 1+s2, item2)
                                win.LCAOTable.setItem(currentRow-3, 1+s2, item3)
                                win.LCAOTable.setItem(currentRow-2, 1+s2, item4)
                                win.LCAOTable.setItem(currentRow-1, 1+s2, item5)
                                head1=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head2=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head3=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head4=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head5=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                win.LCAOTable.setVerticalHeaderItem(currentRow-5, head1)
                                win.LCAOTable.setVerticalHeaderItem(currentRow-4, head2)
                                win.LCAOTable.setVerticalHeaderItem(currentRow-3, head3)
                                win.LCAOTable.setVerticalHeaderItem(currentRow-2, head4)
                                win.LCAOTable.setVerticalHeaderItem(currentRow-1, head5)

                            if orbitLabel=="f":
                                # f orbital: f5z23r2 (1), f5xy2xr2 (cos P), f5yz2yr2 (sin P),
                                #            fzx2zy2 (cos 2P), fxyz (sin 2P), fx33xy2 (cos 3P), f3yx2y3 (sin 3P)
                                f5z23r2 =LCAO_disp[0][0+s2]+LCAO_disp[0][1+s2]*1j
                                f5xy2xr2=LCAO_disp[1][0+s2]+LCAO_disp[1][1+s2]*1j
                                f5yz2yr2=LCAO_disp[2][0+s2]+LCAO_disp[2][1+s2]*1j
                                fzx2zy2 =LCAO_disp[3][0+s2]+LCAO_disp[3][1+s2]*1j
                                fxyz    =LCAO_disp[4][0+s2]+LCAO_disp[4][1+s2]*1j
                                fx33xy2 =LCAO_disp[5][0+s2]+LCAO_disp[5][1+s2]*1j
                                f3yx2y3 =LCAO_disp[6][0+s2]+LCAO_disp[6][1+s2]*1j
                                fm3=(fx33xy2-1j*f3yx2y3)/math.sqrt(2)
                                fm2=(fzx2zy2-1j*fxyz)/math.sqrt(2)
                                fm1=(f5xy2xr2-1j*f5yz2yr2)/math.sqrt(2)
                                fp0=f5z23r2
                                fp1=(f5xy2xr2+1j*f5yz2yr2)/math.sqrt(2)
                                fp2=(fzx2zy2+1j*fxyz)/math.sqrt(2)
                                fp3=(fx33xy2+1j*f3yx2y3)/math.sqrt(2)
                                # table
                                item1=QtGui.QTableWidgetItem(("f5z23r2  ({0:.2f}, {1:.2f})").format(LCAO_disp[0][0+s2], LCAO_disp[0][1+s2]))
                                item2=QtGui.QTableWidgetItem(("f5xy2xr2 ({0:.2f}, {1:.2f})").format(LCAO_disp[1][0+s2], LCAO_disp[1][1+s2]))
                                item3=QtGui.QTableWidgetItem(("f5yz2yr2 ({0:.2f}, {1:.2f})").format(LCAO_disp[2][0+s2], LCAO_disp[2][1+s2]))
                                item4=QtGui.QTableWidgetItem(("fzx2zy2  ({0:.2f}, {1:.2f})").format(LCAO_disp[3][0+s2], LCAO_disp[3][1+s2]))
                                item5=QtGui.QTableWidgetItem(("fxyz     ({0:.2f}, {1:.2f})").format(LCAO_disp[4][0+s2], LCAO_disp[4][1+s2]))
                                item6=QtGui.QTableWidgetItem(("fx33xy2  ({0:.2f}, {1:.2f})").format(LCAO_disp[5][0+s2], LCAO_disp[5][1+s2]))
                                item7=QtGui.QTableWidgetItem(("f3yx2y3  ({0:.2f}, {1:.2f})").format(LCAO_disp[6][0+s2], LCAO_disp[6][1+s2]))
                                win.LCAOTable.setItem(currentRow-7, 0+s2, item1)
                                win.LCAOTable.setItem(currentRow-6, 0+s2, item2)
                                win.LCAOTable.setItem(currentRow-5, 0+s2, item3)
                                win.LCAOTable.setItem(currentRow-4, 0+s2, item4)
                                win.LCAOTable.setItem(currentRow-3, 0+s2, item5)
                                win.LCAOTable.setItem(currentRow-2, 0+s2, item6)
                                win.LCAOTable.setItem(currentRow-1, 0+s2, item7)
                                item1=QtGui.QTableWidgetItem(("f(-3) ({0:.2f}, {1:.2f})").format(fm3.real, fm3.imag))
                                item2=QtGui.QTableWidgetItem(("f(-2) ({0:.2f}, {1:.2f})").format(fm2.real, fm2.imag))
                                item3=QtGui.QTableWidgetItem(("f(-1) ({0:.2f}, {1:.2f})").format(fm1.real, fm1.imag))
                                item4=QtGui.QTableWidgetItem(("f(+0) ({0:.2f}, {1:.2f})").format(fp0.real, fp0.imag))
                                item5=QtGui.QTableWidgetItem(("f(+1) ({0:.2f}, {1:.2f})").format(fp1.real, fp1.imag))
                                item6=QtGui.QTableWidgetItem(("f(+2) ({0:.2f}, {1:.2f})").format(fp2.real, fp2.imag))
                                item7=QtGui.QTableWidgetItem(("f(+3) ({0:.2f}, {1:.2f})").format(fp3.real, fp3.imag))
                                win.LCAOTable.setItem(currentRow-7, 1+s2, item1)
                                win.LCAOTable.setItem(currentRow-6, 1+s2, item2)
                                win.LCAOTable.setItem(currentRow-5, 1+s2, item3)
                                win.LCAOTable.setItem(currentRow-4, 1+s2, item4)
                                win.LCAOTable.setItem(currentRow-3, 1+s2, item5)
                                win.LCAOTable.setItem(currentRow-2, 1+s2, item6)
                                win.LCAOTable.setItem(currentRow-1, 1+s2, item7)
                                head1=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head2=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head3=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head4=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head5=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head6=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head7=QtGui.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                win.LCAOTable.setVerticalHeaderItem(currentRow-7, head1)
                                win.LCAOTable.setVerticalHeaderItem(currentRow-6, head2)
                                win.LCAOTable.setVerticalHeaderItem(currentRow-5, head3)
                                win.LCAOTable.setVerticalHeaderItem(currentRow-4, head4)
                                win.LCAOTable.setVerticalHeaderItem(currentRow-3, head5)
                                win.LCAOTable.setVerticalHeaderItem(currentRow-2, head6)
                                win.LCAOTable.setVerticalHeaderItem(currentRow-1, head7)
                                                                                                
                        mul_index+=1
                        break
                if LCAO_found==False:
                    break
        
