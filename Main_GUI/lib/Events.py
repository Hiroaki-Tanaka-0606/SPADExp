# Events

from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import numpy as np
import re
import h5py
import math
from scipy.stats import norm

import pyqtgraph.opengl as gl
from lib import physical_tools as pt
from lib import objs

from datetime import datetime
import time

import Config

Dispersion=None # dispersion calculated in the specified energy range (numpy array)

def openFile(win, LCAO, Wfns):
    currentFile=win.filePath.text()
    selectedFile, _filter=QtWidgets.QFileDialog.getOpenFileName(caption="Open file", directory=currentFile)
    if selectedFile!="":
        # valid
        win.filePath.setText(selectedFile)
        LCAO.open(selectedFile)
        print("Finished loading LCAO file")

        Wfns.clear()
        for key, value in LCAO.Atom_specs.items():
            Wfns[key]=objs.Wfn(value)
        print("Finished loading AO and PAO file")
            
        win.kxIndex.setMaximum(LCAO.numPnts_kx-1)
        win.kyIndex.setMaximum(LCAO.numPnts_ky-1)
        win.bIndex.setMaximum(LCAO.numBands-1)
        
        win.Atom.clear()
        for atom in LCAO.LCAO_atoms:
            win.Atom.addItem(atom)
        makeOrbitalList(win, LCAO, Wfns)

        print(("Dimension: {0:d}").format(LCAO.Dimension))
        print(("Spin: {0:s}").format(LCAO.Spin))
        print(("EF: {0:.3f} eV").format(Config.Eh*LCAO.EF_Eh))
        print("X range:")
        print(LCAO.Xrange)
        if LCAO.Dimension==2:
            print("Y range:")
            print(LCAO.Yrange)
        print(("Curved: {0:s}").format(str(LCAO.Curved)))

        print("Reciprocal unit cell:")
        print(LCAO.RecCell)
        print("X vector:")
        print(LCAO.Xvector)
        if LCAO.Dimension==2:
            print("Y vector:")
            print(LCAO.Yvector)
        print(" in unit of a.u.^-1")

        if LCAO.Spin.lower()=="on":
            win.UpButton.setCheckable(True)
            win.DnButton.setCheckable(True)
            win.UpButton.setChecked(True)
        else:
            win.UpButton.setCheckable(False)
            win.DnButton.setCheckable(False)

        print("Atoms:")
        for i, atom in enumerate(LCAO.Atoms):
            print(("{0:d} {1:s} {2:s} {3:s}").format(i, atom, LCAO.Atom_specs[atom][0], LCAO.Atom_specs[atom][1]))
            

        print("Coordinates of atoms:")
        print(LCAO.Atom_coordinates)
        print((" in unit of {0:s}").format(LCAO.Atom_unit))

def appendDispersion(i, n, EMin, EPixel, tailProfile, LCAO, plotPAD, PADobj):
    global Dispersion

    ret=0
    Esize=Dispersion.shape[1]
    tailSize=tailProfile.shape[0]
    if LCAO.Spin_i==1:
        continueFlag=False
        # spin Up
        eigen=(LCAO.BandUp[i][n]-LCAO.EF_Eh)*Config.Eh
        eigen_index=round((eigen-EMin)/EPixel)
        if eigen_index-tailSize>=Esize:
            pass
        else:
            continueFlag=True
            PAD=1.0
            if plotPAD==True and eigen_index+tailSize-1>=0:
                PAD=PADobj.calc(i, n, useUp=True)
                
            for j in range(-tailSize+1, tailSize):
                if eigen_index+j>=0 and eigen_index+j<Esize:
                    Dispersion[i][eigen_index+j]+=tailProfile[abs(j)]*PAD
        # spin Dn
        eigen=(LCAO.BandDn[i][n]-LCAO.EF_Eh)*Config.Eh
        eigen_index=round((eigen-EMin)/EPixel)
        if eigen_index-tailSize>=Esize:
            pass
        else:
            continueFlag=True
            PAD=1.0
            if plotPAD==True and eigen_index+tailSize-1>=0:
                PAD=PADobj.calc(i, n, useUp=False)
                
            for j in range(-tailSize+1, tailSize):
                if eigen_index+j>=0 and eigen_index+j<Esize:
                    Dispersion[i][eigen_index+j]+=tailProfile[abs(j)]*PAD

        return continueFlag
    else:
        eigen=(LCAO.Band[i][n]-LCAO.EF_Eh)*Config.Eh
        eigen_index=round((eigen-EMin)/EPixel)
        if eigen_index-tailSize>=Esize:
            return False
        PAD=1.0
        if plotPAD==True and eigen_index+tailSize-1>=0:
            PAD=PADobj.calc(i, n)
            
        for j in range(-tailSize+1, tailSize):
            if eigen_index+j>=0 and eigen_index+j<Esize:
                Dispersion[i][eigen_index+j]+=tailProfile[abs(j)]*PAD
        return True
        

def appendDispersion3(ix, iy, n, EMin, EPixel, tailProfile, LCAO, plotPAD, PADobj):
    global Dispersion

    ret=0
    Esize=Dispersion.shape[2]
    tailSize=tailProfile.shape[0]

    i=ix+iy*LCAO.numPnts_kx
    if LCAO.Spin_i==1:
        continueFlag=False
        
        eigen=(LCAO.BandUp[i][n]-LCAO.EF_Eh)*Config.Eh
        eigen_index=round((eigen-EMin)/EPixel)
        if eigen_index-tailSize>=Esize:
            pass
        else:
            continueFlag=True
            PAD=1.0
            if plotPAD==True and eigen_index+tailSize-1>=0:
                PAD=PADobj.calc(i, n, useUp=True)
                
            for j in range(-tailSize+1, tailSize):
                if eigen_index+j>=0 and eigen_index+j<Esize:
                    Dispersion[ix][iy][eigen_index+j]+=tailProfile[abs(j)]*PAD
        
        eigen=(LCAO.BandDn[i][n]-LCAO.EF_Eh)*Config.Eh
        eigen_index=round((eigen-EMin)/EPixel)
        if eigen_index-tailSize>=Esize:
            pass
        else:
            continueFlag=True
            PAD=1.0
            if plotPAD==True and eigen_index+tailSize-1>=0:
                PAD=PADobj.calc(i, n, useUp=False)
                
            for j in range(-tailSize+1, tailSize):
                if eigen_index+j>=0 and eigen_index+j<Esize:
                    Dispersion[ix][iy][eigen_index+j]+=tailProfile[abs(j)]*PAD

        return continueFlag

    else:
        eigen=(LCAO.Band[i][n]-LCAO.EF_Eh)*Config.Eh
        eigen_index=round((eigen-EMin)/EPixel)
        if eigen_index-tailSize>=Esize:
            return False
        PAD=1.0
        if plotPAD==True and eigen_index+tailSize-1>=0:
            # time1=time.time()
            PAD=PADobj.calc(i, n) 
            # time2=time.time()
            # print("PAD", time2-time1) about 0.013 sec
            
        for j in range(-tailSize+1, tailSize):
            if eigen_index+j>=0 and eigen_index+j<Esize:
                Dispersion[ix][iy][eigen_index+j]+=tailProfile[abs(j)]*PAD
                
        return True
                

def genTailProfile(EPixel, dE):
    tailIndex=math.floor(dE*Config.sigma_max/EPixel)
    ret=np.zeros((tailIndex+1,))
    for i in range(0, tailIndex+1):
        ret[i]=norm.pdf(i*EPixel, loc=0, scale=dE)

    return ret


def plot(win, LCAO, PADobj):
    # for 2D plot
    global Dispersion

    EMin=float(win.EMin.text())
    EMax=float(win.EMax.text())
    dE=float(win.dE.text())
    EPixel=float(win.EPixel.text())

    tailProfile=genTailProfile(EPixel, dE)

    numPnts_E=math.ceil((EMax-EMin)/EPixel+1)
    if numPnts_E<0:
        print("Energy range error")
        return

    plotPAD=False
    initialStates_i=1 # 0->PAO, 1->AO
    finalStates_i=0 # 0->Plane wave, 1->Calculated
    polarization_i=0 # 0->Linear, 1->Right circular, 2->Left circular
    finalStates_step=0.0
    theta=0.0
    phi=0.0
    Y_coeff=[0, 0, 0] # coeffcients of operators r Y_{1,m} [m=-1, m=0, m=1]
    if win.plotDispersion.isChecked():
        print("Plot the band dispersion")
    elif win.plotPAD.isChecked():
        plotPAD=True
        print("Plot the PAD")
        # initial state
        if win.AOButton.isChecked():
            print("Initial state: atomic orbital")
        elif win.PAOButton.isChecked():
            print("Initial state: pseudo-atomic orbital")
            initialStates_i=0
        else:
            print("Error: AO or PAO should be selected")
            return        
        # final state
        if win.PWButton.isChecked():
            print("Final state: plane wave")
        elif win.CalcButton.isChecked():
            print("Final state: calculated wavefunction")
            finalStates_i=1
            finalStates_step=float(win.finalStates_step.text())
        else:
            print("Error: final state is not selected")
            return
        # polarization
        if win.linear.isChecked():
            print("Polarization: linear")
        elif win.rCircular.isChecked():
            print("Polarization: right circular")
            polarization_i=1
        elif win.lCircular.isChecked():
            print("Polarization: left circular")
            polarization_i=2
        else:
            print("Error: polarization is not selected")
            return
        # angle
        theta=float(win.theta.text())
        phi=float(win.phi.text())
        print(("Angle: theta={0:.2f} deg, phi={1:.2f} deg").format(theta, phi))
            
        pt.calcOperatorCoeff(Y_coeff, polarization_i, theta, phi)
        PADobj.setSystem(initialStates_i, finalStates_i, finalStates_step, Y_coeff)
            
    else:
        print("Error: Band dispersion or PAD should be checked")
        return
    

    print(("{0:d} points along the energy").format(numPnts_E))
    print(("{0:d} points along the kx").format(LCAO.numPnts_kx))

    if LCAO.Dimension!=1:
        print("Dimension error")
        return
    
    Dispersion=np.zeros((LCAO.numPnts_kx, numPnts_E))
    
    for i in range(0, LCAO.numPnts_kx):
        print(("Calculating k = {0:6d}").format(i))
        for n in range(0, LCAO.numBands):
            if appendDispersion(i, n, EMin, EPixel, tailProfile, LCAO, plotPAD, PADobj):
                continue
            else:
                break

    tr=QtGui.QTransform()
    tr.translate(LCAO.Xlength*LCAO.Xrange[0]-LCAO.dx_length/2,EMin-EPixel/2)
    tr.scale(LCAO.dx_length, EPixel)
    win.img.setTransform(tr)
    win.img.setImage(Dispersion)

    maxpoint=Dispersion.max()
    print(("Max point: {0:8.3e}").format(maxpoint))
    win.bar.setLevels((0, maxpoint))
            
def makeDispersion3(win, LCAO, PADobj):
    # for 3D plot
    global Dispersion

    EMin=float(win.EMin.text())
    EMax=float(win.EMax.text())
    dE=float(win.dE.text())
    EPixel=float(win.EPixel.text())

    tailProfile=genTailProfile(EPixel, dE)

    numPnts_E=math.ceil((EMax-EMin)/EPixel+1)
    if numPnts_E<0:
        print("Energy range error")
        return

    plotPAD=False
    initialStates_i=1 # 0->PAO, 1->AO
    finalStates_i=0 # 0->Plane wave, 1->Calculated
    polarization_i=0 # 0->Linear, 1->Right circular, 2->Left circular
    finalStates_step=0.0
    theta=0.0
    phi=0.0
    Y_coeff=[0, 0, 0] # coeffcients of operators r Y_{1,m} [m=-1, m=0, m=1]
    if win.plotDispersion.isChecked():
        print("Plot the band dispersion")
    elif win.plotPAD.isChecked():
        plotPAD=True
        print("Plot the PAD")
        # initial state
        if win.AOButton.isChecked():
            print("Initial state: atomic orbital")
        elif win.PAOButton.isChecked():
            print("Initial state: pseudo-atomic orbital")
            initialStates_i=0
        else:
            print("Error: PAO or AO should be selected")
            return        
        # final state
        if win.PWButton.isChecked():
            print("Final state: plane wave")
        elif win.CalcButton.isChecked():
            print("Final state: calculated wavefunction")
            finalStates_i=1
            finalStates_step=float(win.finalStates_step.text())
        else:
            print("Error: final state is not selected")
            return
        # polarization
        if win.linear.isChecked():
            print("Polarization: linear")
        elif win.rCircular.isChecked():
            print("Polarization: right circular")
            polarization_i=1
        elif win.lCircular.isChecked():
            print("Polarization: left circular")
            polarization_i=2
        else:
            print("Error: polarization is not selected")
            return
        # angle
        theta=float(win.theta.text())
        phi=float(win.phi.text())
        print(("Angle: theta={0:.2f} deg, phi={1:.2f} deg").format(theta, phi))
            
        pt.calcOperatorCoeff(Y_coeff, polarization_i, theta, phi)
        PADobj.setSystem(initialStates_i, finalStates_i, finalStates_step, Y_coeff)
            
    else:
        print("Error: Band dispersion or PAD should be checked")
        return
    
    print(("{0:d} points along the energy").format(numPnts_E))

    win.eIndex.setMaximum(numPnts_E-1)

    print(("{0:d} points along the kx").format(LCAO.numPnts_kx))
    print(("{0:d} points along the ky").format(LCAO.numPnts_ky))

    if LCAO.Dimension!=2:
        print("Dimension error")
        return
    
    Dispersion=np.zeros((LCAO.numPnts_kx, LCAO.numPnts_ky, numPnts_E))
    # time1=time.time()
    for i in range(0, LCAO.numPnts_kx):
        print(("Calculating kx = {0:6d}").format(i))
        for j in range(0, LCAO.numPnts_ky):
            # print(("  Calculating ky = {0:6d}").format(j))
            for n in range(0, LCAO.numBands):
                # time2=time.time()
                # print(time2-time1) about 0.013 sec
                # time1=time2
                if appendDispersion3(i, j, n, EMin, EPixel, tailProfile, LCAO, plotPAD, PADobj):
                    continue
                else:
                    break


    tr_x=QtGui.QTransform()
    tr_x.translate(LCAO.Xlength*LCAO.Xrange[0]-LCAO.dx_length/2,EMin-EPixel/2)
    tr_x.scale(LCAO.dx_length, EPixel)
    win.imgEx.setTransform(tr_x)

    tr_y=QtGui.QTransform()
    tr_y.translate(EMin-EPixel/2,LCAO.Ylength*LCAO.Yrange[0]-LCAO.dy_length/2)
    tr_y.rotate(-90)
    tr_y.scale(-LCAO.dy_length, EPixel)
    win.imgEy.setTransform(tr_y)

    tr_E=QtGui.QTransform()
    tr_E.translate(LCAO.Xlength*LCAO.Xrange[0]-LCAO.dx_length/2,LCAO.Ylength*LCAO.Yrange[0]-LCAO.dy_length/2)
    tr_E.scale(LCAO.dx_length, LCAO.dy_length)
    win.imgxy.setTransform(tr_E)

    Maxpoint=Dispersion.max()

    win.Cube=np.zeros((LCAO.numPnts_kx, LCAO.numPnts_ky, numPnts_E, 4))
    win.Cube[:,:,:,0]=255
    win.Cube[:,:,:,1]=255
    win.Cube[:,:,:,2]=255
    win.Cube[:,:,:,3]=Dispersion/Maxpoint*100
    
    win.bandCube=gl.GLVolumeItem(win.Cube)
    win.bandCube.scale(LCAO.dx_length, LCAO.dx_length, EPixel)    
    win.bandCube.translate(LCAO.Xlength*LCAO.Xrange[0]-LCAO.dx_length/2,LCAO.Ylength*LCAO.Yrange[0]-LCAO.dy_length/2,0)
    win.plot3D.clear()
    win.plot3D.addItem(win.bandCube)

def plot3(win, LCAO):

    global Dispersion

    kx=win.kxIndex.value()
    ky=win.kyIndex.value()
    ei=win.eIndex.value()
    ExMax=Dispersion[:,ky,:].max()
    EyMax=Dispersion[kx,:,:].max()
    xyMax=Dispersion[:,:,ei].max()
    MaxPoint=max(ExMax, EyMax, xyMax)
    win.imgEx.setImage(Dispersion[:,ky,:])
    win.imgEy.setImage(Dispersion[kx,:,:])
    win.imgxy.setImage(Dispersion[:,:,ei])

    Maxpoint=Dispersion.max()
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
            
def drawCursor(win, LCAO):
    
    k=win.kxIndex.value()
    b=win.bIndex.value()
    
    UseUp=False
    if LCAO.Spin_i==1:
        if win.UpButton.isChecked():
            UseUp=True
        elif win.DnButton.isChecked():
            UseUp=False
        else:
            print("None of Up and Dn is checked")
            return
    
    if 0<=k and k<LCAO.numPnts_kx and 0<=b and b<LCAO.numBands:
        k_value=LCAO.Xlength*LCAO.Xrange[0]+LCAO.dx_length*k
        b_value=0
        if LCAO.Spin_i==1:
            if UseUp:
                b_value=LCAO.BandUp[k][b]
            else:
                b_value=LCAO.BandDn[k][b]
        else:
            b_value=LCAO.Band[k][b]

        b_value=(b_value-LCAO.EF_Eh)*Config.Eh
        win.vLine.setPos(k_value)
        win.hLine.setPos(b_value)
        win.kxValue.setText(("({0:.3f})").format(k_value))
        win.bValue.setText(("({0:.3f})").format(b_value))
        
    else:
        print("Index error")
        return
            
def drawCursor3(win, LCAO):
    EMin=float(win.EMin.text())
    EPixel=float(win.EPixel.text())
    
    kx=win.kxIndex.value()
    ky=win.kyIndex.value()
    b=win.bIndex.value()
    ei=win.eIndex.value()
    UseUp=False
    if LCAO.Spin_i==1:
        if win.UpButton.isChecked():
            UseUp=True
        elif win.DnButton.isChecked():
            UseUp=False
        else:
            print("None of Up and Dn is checked")
            return

    k=kx+ky*LCAO.numPnts_kx

    if 0<=k and k<LCAO.numPnts_k and 0<=b and b<LCAO.numBands:
        kx_value=LCAO.Xlength*LCAO.Xrange[0]+LCAO.dx_length*kx
        ky_value=LCAO.Ylength*LCAO.Yrange[0]+LCAO.dy_length*ky
        e_value=EMin+EPixel*ei
        b_value=0
        if LCAO.Spin_i==1:
            if UseUp:
                b_value=LCAO.BandUp[k][b]
            else:
                b_value=LCAO.BandDn[k][b]
        else:
            b_value=LCAO.Band[k][b]

        b_value=(b_value-LCAO.EF_Eh)*Config.Eh
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

def makeOrbitalList(win, LCAO, Wfns):
    win.orbitalToPlot.clear()

    at=win.Atom.currentIndex()
    at_label=LCAO.Atoms[at]
    orbits=Wfns[at_label].Orbits
    for orbit in orbits:
        win.orbitalToPlot.addItem(orbit)


def plotOrbital(win, LCAO, Wfns, PADobj):
    at=win.Atom.currentIndex()
    at_label=LCAO.Atoms[at]
    orb=win.orbitalToPlot.currentIndex()
    
    kx=win.kxIndex.value()
    ky=win.kyIndex.value()
    k=0
    if LCAO.Dimension==1:
        k=kx
    elif LCAO.Dimension==2:
        k=kx+ky*LCAO.numPnts_kx
    
    k_au=LCAO.Kpath_au[k]
    k_length=math.sqrt(np.inner(k_au, k_au))
        
    wfn=Wfns[at_label].Wfn[orb]
    wfn_finalp1=np.zeros((Wfns[at_label].length,))
    wfn_finalm1=np.zeros((Wfns[at_label].length,))
    orbit_label=Wfns[at_label].Orbits[orb]
    l=0
    if orbit_label[0]=="s":
        pass
    elif orbit_label[0]=="p":
        l=1
    elif orbit_label[0]=="d":
        l=2
    elif orbit_label[0]=="f":
        l=3
 
    r=Wfns[at_label].r
    PADobj.calcFinalState(wfn_finalp1, l+1, k_length, r)
    if l-1>=0:
        PADobj.calcFinalState(wfn_finalm1, l-1, k_length, r)

    win.wfnPlot.clear()

    win.wfnPlot.plot(y=wfn[:][0], x=r, name="PAO", pen=Config.pen_PAO)
    win.wfnPlot.plot(y=wfn[:][1], x=r, name="AO", pen=Config.pen_AO)
    win.wfnPlot.plot(y=wfn_finalp1, x=r, name="Final (l+1)", pen=Config.pen_finalp1)
    if l-1>=0:
        win.wfnPlot.plot(y=wfn_finalm1, x=r, name="Final (l-1)", pen=Config.pen_finalm1)
    

    
def makeLCAOTable(win, LCAO):
    
    kx=win.kxIndex.value()
    ky=win.kyIndex.value()
    k=0
    if LCAO.Dimension==1:
        k=kx
    elif LCAO.Dimension==2:
        k=kx+ky*LCAO.numPnts_kx

    b=win.bIndex.value()
    at=win.Atom.currentIndex()
    
    UseUp=False
    if LCAO.Spin_i==1:
        if win.UpButton.isChecked():
            UseUp=True
        elif win.DnButton.isChecked():
            UseUp=False
        else:
            print("None of Up and Dn is checked")
            return
        
    if len(LCAO.LCAO_labels)!=len(LCAO.Atoms):
        return

    win.LCAOTable.setRowCount(0)
    if LCAO.Spin_i==2:
        win.LCAOTable.setColumnCount(4)
        win.LCAOTable.setHorizontalHeaderItem(0, QtWidgets.QTableWidgetItem("Up (raw)"))
        win.LCAOTable.setHorizontalHeaderItem(1, QtWidgets.QTableWidgetItem("Up (calc)"))
        win.LCAOTable.setHorizontalHeaderItem(2, QtWidgets.QTableWidgetItem("Dn (raw)"))
        win.LCAOTable.setHorizontalHeaderItem(3, QtWidgets.QTableWidgetItem("Dn (calc)"))
    else:
        win.LCAOTable.setColumnCount(2)
        if LCAO.Spin_i==1:
            win.LCAOTable.setHorizontalHeaderItem(0, QtWidgets.QTableWidgetItem(("{0:s} (raw)").format("Up" if UseUp else "Dn")))
            win.LCAOTable.setHorizontalHeaderItem(1, QtWidgets.QTableWidgetItem(("{0:s} (calc)").format("Up" if UseUp else "Dn")))
        else:
            win.LCAOTable.setHorizontalHeaderItem(0, QtWidgets.QTableWidgetItem("raw"))
            win.LCAOTable.setHorizontalHeaderItem(1, QtWidgets.QTableWidgetItem("calc"))

    currentRow=0
    if 0<=k and k<LCAO.numPnts_k and 0<=b and b<LCAO.numBands and 0<=at and at<len(LCAO.Atoms):
        orbitLabels=["s", "p", "d", "f"]
        for orbitLabel in orbitLabels:
            mul_index=0
            while(True):
                LCAO_label=("{0:s}{1:1d}").format(orbitLabel, mul_index)
                if LCAO.Spin_i==1:
                    if UseUp:
                        LCAO_label+="Up"
                    else:
                        LCAO_label+="Dn"
                    
                LCAO_found=False
                for i, LCAO_label_ref in enumerate(LCAO.LCAO_labels[at]):
                    if LCAO_label==LCAO_label_ref:
                        LCAO_found=True

                        LCAO_disp=LCAO.LCAO[at][i][k][b]
                        numSpin=1
                        if LCAO.Spin_i==2:
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
                        LCAO_conv=np.zeros((7,), dtype=complex)
                        for s in range(0, numSpin):
                            s2=s*2
                            # see OpenMX/source/AngularF.c for the order of the spherical harmonics
                            if orbitLabel=="s":
                                # s orbital: nothing to calculate
                                item1=QtWidgets.QTableWidgetItem(("s ({0:.3f}, {1:.3f})").format(LCAO_disp[0][0+s2], LCAO_disp[0][1+s2]))
                                win.LCAOTable.setItem(currentRow-1, 0+s2, item1)
                                item2=QtWidgets.QTableWidgetItem(("s ({0:.3f}, {1:.3f})").format(LCAO_disp[0][0+s2], LCAO_disp[0][1+s2]))
                                win.LCAOTable.setItem(currentRow-1, 1+s2, item2)
                                head=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                win.LCAOTable.setVerticalHeaderItem(currentRow-1, head)
                                
                            if orbitLabel=="p":
                                # p orbital: px (cos P), py (sin P), pz (1)
                                px=LCAO_disp[0][0+s2]+LCAO_disp[0][1+s2]*1j
                                py=LCAO_disp[1][0+s2]+LCAO_disp[1][1+s2]*1j
                                pz=LCAO_disp[2][0+s2]+LCAO_disp[2][1+s2]*1j
                                pt.convertLCAO_p(px, py, pz, LCAO_conv)
                                # table
                                item1=QtWidgets.QTableWidgetItem(("px ({0:.3f}, {1:.3f})").format(LCAO_disp[0][0+s2], LCAO_disp[0][1+s2]))
                                item2=QtWidgets.QTableWidgetItem(("py ({0:.3f}, {1:.3f})").format(LCAO_disp[1][0+s2], LCAO_disp[1][1+s2]))
                                item3=QtWidgets.QTableWidgetItem(("pz ({0:.3f}, {1:.3f})").format(LCAO_disp[2][0+s2], LCAO_disp[2][1+s2]))
                                win.LCAOTable.setItem(currentRow-3, 0+s2, item1)
                                win.LCAOTable.setItem(currentRow-2, 0+s2, item2)
                                win.LCAOTable.setItem(currentRow-1, 0+s2, item3)
                                item1=QtWidgets.QTableWidgetItem(("p(-1) ({0:.3f}, {1:.3f})").format(LCAO_conv[0].real, LCAO_conv[0].imag))
                                item2=QtWidgets.QTableWidgetItem(("p(+0) ({0:.3f}, {1:.3f})").format(LCAO_conv[1].real, LCAO_conv[1].imag))
                                item3=QtWidgets.QTableWidgetItem(("p(+1) ({0:.3f}, {1:.3f})").format(LCAO_conv[2].real, LCAO_conv[2].imag))
                                win.LCAOTable.setItem(currentRow-3, 1+s2, item1)
                                win.LCAOTable.setItem(currentRow-2, 1+s2, item2)
                                win.LCAOTable.setItem(currentRow-1, 1+s2, item3)
                                head1=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head2=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head3=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
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
                                pt.convertLCAO_d(d3z2r2, dx2y2, dxy, dxz, dyz, LCAO_conv)
                                # table
                                item1=QtWidgets.QTableWidgetItem(("d3z2r2 ({0:.3f}, {1:.3f})").format(LCAO_disp[0][0+s2], LCAO_disp[0][1+s2]))
                                item2=QtWidgets.QTableWidgetItem(("dx2y2  ({0:.3f}, {1:.3f})").format(LCAO_disp[1][0+s2], LCAO_disp[1][1+s2]))
                                item3=QtWidgets.QTableWidgetItem(("dxy    ({0:.3f}, {1:.3f})").format(LCAO_disp[2][0+s2], LCAO_disp[2][1+s2]))
                                item4=QtWidgets.QTableWidgetItem(("dxz    ({0:.3f}, {1:.3f})").format(LCAO_disp[3][0+s2], LCAO_disp[3][1+s2]))
                                item5=QtWidgets.QTableWidgetItem(("dyz    ({0:.3f}, {1:.3f})").format(LCAO_disp[4][0+s2], LCAO_disp[4][1+s2]))
                                win.LCAOTable.setItem(currentRow-5, 0+s2, item1)
                                win.LCAOTable.setItem(currentRow-4, 0+s2, item2)
                                win.LCAOTable.setItem(currentRow-3, 0+s2, item3)
                                win.LCAOTable.setItem(currentRow-2, 0+s2, item4)
                                win.LCAOTable.setItem(currentRow-1, 0+s2, item5)
                                item1=QtWidgets.QTableWidgetItem(("d(-2) ({0:.3f}, {1:.3f})").format(LCAO_conv[0].real, LCAO_conv[0].imag))
                                item2=QtWidgets.QTableWidgetItem(("d(-1) ({0:.3f}, {1:.3f})").format(LCAO_conv[1].real, LCAO_conv[1].imag))
                                item3=QtWidgets.QTableWidgetItem(("d(+0) ({0:.3f}, {1:.3f})").format(LCAO_conv[2].real, LCAO_conv[2].imag))
                                item4=QtWidgets.QTableWidgetItem(("d(+1) ({0:.3f}, {1:.3f})").format(LCAO_conv[3].real, LCAO_conv[3].imag))
                                item5=QtWidgets.QTableWidgetItem(("d(+2) ({0:.3f}, {1:.3f})").format(LCAO_conv[4].real, LCAO_conv[4].imag))
                                win.LCAOTable.setItem(currentRow-5, 1+s2, item1)
                                win.LCAOTable.setItem(currentRow-4, 1+s2, item2)
                                win.LCAOTable.setItem(currentRow-3, 1+s2, item3)
                                win.LCAOTable.setItem(currentRow-2, 1+s2, item4)
                                win.LCAOTable.setItem(currentRow-1, 1+s2, item5)
                                head1=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head2=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head3=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head4=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head5=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
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
                                pt.convertLCAO_f(f5z23r2, f5xy2xr2, f5yz2yr2, fzx2zy2, fxyz, fx33xy2, f3yx2y3, LCAO_conv)
                                # table
                                item1=QtWidgets.QTableWidgetItem(("f5z23r2  ({0:.3f}, {1:.3f})").format(LCAO_disp[0][0+s2], LCAO_disp[0][1+s2]))
                                item2=QtWidgets.QTableWidgetItem(("f5xy2xr2 ({0:.3f}, {1:.3f})").format(LCAO_disp[1][0+s2], LCAO_disp[1][1+s2]))
                                item3=QtWidgets.QTableWidgetItem(("f5yz2yr2 ({0:.3f}, {1:.3f})").format(LCAO_disp[2][0+s2], LCAO_disp[2][1+s2]))
                                item4=QtWidgets.QTableWidgetItem(("fzx2zy2  ({0:.3f}, {1:.3f})").format(LCAO_disp[3][0+s2], LCAO_disp[3][1+s2]))
                                item5=QtWidgets.QTableWidgetItem(("fxyz     ({0:.3f}, {1:.3f})").format(LCAO_disp[4][0+s2], LCAO_disp[4][1+s2]))
                                item6=QtWidgets.QTableWidgetItem(("fx33xy2  ({0:.3f}, {1:.3f})").format(LCAO_disp[5][0+s2], LCAO_disp[5][1+s2]))
                                item7=QtWidgets.QTableWidgetItem(("f3yx2y3  ({0:.3f}, {1:.3f})").format(LCAO_disp[6][0+s2], LCAO_disp[6][1+s2]))
                                win.LCAOTable.setItem(currentRow-7, 0+s2, item1)
                                win.LCAOTable.setItem(currentRow-6, 0+s2, item2)
                                win.LCAOTable.setItem(currentRow-5, 0+s2, item3)
                                win.LCAOTable.setItem(currentRow-4, 0+s2, item4)
                                win.LCAOTable.setItem(currentRow-3, 0+s2, item5)
                                win.LCAOTable.setItem(currentRow-2, 0+s2, item6)
                                win.LCAOTable.setItem(currentRow-1, 0+s2, item7)
                                item1=QtWidgets.QTableWidgetItem(("f(-3) ({0:.3f}, {1:.3f})").format(LCAO_conv[0].real, LCAO_conv[0].imag))
                                item2=QtWidgets.QTableWidgetItem(("f(-2) ({0:.3f}, {1:.3f})").format(LCAO_conv[1].real, LCAO_conv[1].imag))
                                item3=QtWidgets.QTableWidgetItem(("f(-1) ({0:.3f}, {1:.3f})").format(LCAO_conv[2].real, LCAO_conv[2].imag))
                                item4=QtWidgets.QTableWidgetItem(("f(+0) ({0:.3f}, {1:.3f})").format(LCAO_conv[3].real, LCAO_conv[3].imag))
                                item5=QtWidgets.QTableWidgetItem(("f(+1) ({0:.3f}, {1:.3f})").format(LCAO_conv[4].real, LCAO_conv[4].imag))
                                item6=QtWidgets.QTableWidgetItem(("f(+2) ({0:.3f}, {1:.3f})").format(LCAO_conv[5].real, LCAO_conv[5].imag))
                                item7=QtWidgets.QTableWidgetItem(("f(+3) ({0:.3f}, {1:.3f})").format(LCAO_conv[6].real, LCAO_conv[6].imag))
                                win.LCAOTable.setItem(currentRow-7, 1+s2, item1)
                                win.LCAOTable.setItem(currentRow-6, 1+s2, item2)
                                win.LCAOTable.setItem(currentRow-5, 1+s2, item3)
                                win.LCAOTable.setItem(currentRow-4, 1+s2, item4)
                                win.LCAOTable.setItem(currentRow-3, 1+s2, item5)
                                win.LCAOTable.setItem(currentRow-2, 1+s2, item6)
                                win.LCAOTable.setItem(currentRow-1, 1+s2, item7)
                                head1=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head2=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head3=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head4=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head5=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head6=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
                                head7=QtWidgets.QTableWidgetItem(("{0:s}{1:1d}").format(orbitLabel, mul_index))
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
# real-space image
# in unit of Ang, so conversion is necessary
def makeRealSpace(win, LCAO, Elements):
    # ia: atom
    # ie: element
    boundaryA=win.boundaryA.value()
    boundaryB=win.boundaryB.value()
    boundaryC=win.boundaryC.value()
    win.realSpace.clear()
    # print(LCAO.atomCell_au)
    
    au_ang=Config.au_ang

    # unit cell
    pts=np.zeros((5,3))
    # base plane
    pts[1]=LCAO.AtomCell_au[0]*au_ang
    pts[2]=(LCAO.AtomCell_au[0]+LCAO.AtomCell_au[1])*au_ang
    pts[3]=LCAO.AtomCell_au[1]*au_ang
    
    plt=gl.GLLinePlotItem(pos=pts)
    win.realSpace.addItem(plt)
    
    #top plane
    pts2=pts.copy()
    for i in range(0, 5):
        pts2[i]+=LCAO.AtomCell_au[2]*au_ang
    plt=gl.GLLinePlotItem(pos=pts2)
    win.realSpace.addItem(plt)

    # vertical pillars
    for i in range(0, 4):
        pts3=np.zeros((2, 3))
        for j in range(0, 3):
            pts3[0][j]=pts[i][j]
            pts3[1][j]=pts[i][j]+LCAO.AtomCell_au[2][j]*au_ang
        plt=gl.GLLinePlotItem(pos=pts3)
        win.realSpace.addItem(plt)

    # kx and ky vector (direction)
    kx=np.zeros((2,3))
    kx[1]=LCAO.Xvector*Config.reciprocal_coeff/au_ang
    kx_color=np.zeros((2,4))
    kx_color[:,0]=Config.pen_kx[0]
    kx_color[:,1]=Config.pen_kx[1]
    kx_color[:,2]=Config.pen_kx[2]
    kx_color[:,3]=1.0
    plt=gl.GLLinePlotItem(pos=kx, color=kx_color, width=Config.reciprocal_axis_width)
    win.realSpace.addItem(plt)
    if LCAO.Dimension==2:
        ky=np.zeros((2,3))
        ky[1]=LCAO.Yvector*Config.reciprocal_coeff/au_ang
        ky_color=np.zeros((2,4))
        ky_color[:,0]=Config.pen_ky[0]
        ky_color[:,1]=Config.pen_ky[1]
        ky_color[:,2]=Config.pen_ky[2]
        ky_color[:,3]=1.0
        plt=gl.GLLinePlotItem(pos=ky, color=ky_color, width=Config.reciprocal_axis_width)
        win.realSpace.addItem(plt)

    # Polarization
    pol=np.zeros((2,3))
    theta=math.radians(float(win.theta.text()))
    phi=math.radians(float(win.phi.text()))
    pol[1][0]=math.sin(theta)*math.cos(phi)
    pol[1][1]=math.sin(theta)*math.sin(phi)
    pol[1][2]=math.cos(theta)
    pol[1]*=Config.polarization_length
    pol_color=np.zeros((2,4))
    pol_color[:,0]=Config.pen_pol[0]
    pol_color[:,1]=Config.pen_pol[1]
    pol_color[:,2]=Config.pen_pol[2]
    pol_color[:,3]=1
    plt=gl.GLLinePlotItem(pos=pol, color=pol_color, width=Config.polarization_width)
    win.realSpace.addItem(plt)
    
    # atoms
    for ia, atom_label in enumerate(LCAO.Atoms):
        # print(atom_label)
        el_name=""
        el_index=-1
        for ie, el_label in enumerate(Elements.labels):
            el_match=re.findall(r"^"+el_label, atom_label)
            if len(el_match)>0 and len(el_match[0])>len(el_name):
                el_name=el_match[0]
                el_index=ie

        if el_index==-1:
            print(("Error: atom {0:s} not found").format(atom_label))
            el_index=Config.not_found_element

        r=Elements.radii[el_index][Config.radius_index]*Config.radius_coeff
        color=Elements.colors[el_index]
            
        md=gl.MeshData.sphere(rows=10, cols=20, radius=r)
        meshcolor=np.zeros((md.faceCount(), 4), dtype=float)
        meshcolor[:,0]=color[0]
        meshcolor[:,1]=color[1]
        meshcolor[:,2]=color[2]
        meshcolor[:,3]=1.0
        # print(meshcolor)

        md.setFaceColors(meshcolor)

        for iA in range(0, boundaryA):
            for iB in range(0, boundaryB):
                for iC in range(0, boundaryC):
                    mi=gl.GLMeshItem(meshdata=md, smooth=False)
                    coordinate=LCAO.Atom_au[ia].copy()
                    coordinate+=iA*LCAO.AtomCell_au[0]
                    coordinate+=iB*LCAO.AtomCell_au[1]
                    coordinate+=iC*LCAO.AtomCell_au[2]
                    
                    mi.translate(coordinate[0]*au_ang, coordinate[1]*au_ang, coordinate[2]*au_ang)
                    win.realSpace.addItem(mi)



def export(win, LCAO, PADobj):
    currentFile=win.filePath.text()
    selectedFile, _filter=QtWidgets.QFileDialog.getSaveFileName(caption="Open file", directory=currentFile)
    if selectedFile!="":
        with h5py.File(selectedFile, "w") as f:
            f.attrs.create("Datetime", datetime.now().isoformat(" "))
            f.create_dataset("Dispersion", data=Dispersion)
            f.attrs.create("Dimension", LCAO.Dimension)
            f.attrs.create("Weighting", False)

            if LCAO.Dimension==1:
                offset=[LCAO.Xlength*LCAO.Xrange[0], float(win.EMin.text())]
                f.attrs.create("Offset", offset)
                delta=[LCAO.dx_length, float(win.EPixel.text())]
                f.attrs.create("Delta", delta)
                f.attrs.create("Xvector", LCAO.Xvector)
            else:
                offset=[LCAO.Xlength*LCAO.Xrange[0], LCAO.Ylength*LCAO.Yrange[0], float(win.EMin.text())]
                f.attrs.create("Offset", offset)
                delta=[LCAO.dx_length, LCAO.dy_length, float(win.EPixel.text())]
                f.attrs.create("Delta", delta)
                f.attrs.create("Xvector", LCAO.Xvector)
                f.attrs.create("Yvector", LCAO.Yvector)
            
            size=Dispersion.shape
            f.attrs.create("Size", size)
            dE=float(win.dE.text())
            f.attrs.create("dE", dE)

            initialStates=(["PAO", "AO"])[PADobj.initialStates_i]
            finalStates=(["PW", "Calc"])[PADobj.finalStates_i]
            f.attrs.create("Initial_state", initialStates)
            f.attrs.create("Final_state", finalStates)

            polarization=""
            if win.linear.isChecked():
                polarization="Linear"
            elif win.rCircular.isChecked():
                polarization="RCircular"
            elif win.lCircular.isChecked():
                polarzation="LCircular"
            f.attrs.create("Polarization", polarization)

            theta=float(win.theta.text())
            phi=float(win.phi.text())
            f.attrs.create("Theta", theta)
            f.attrs.create("Phi", phi)

            atomG=f.create_group("Atoms")
            atomG.create_dataset("Labels", data=LCAO.Atoms)
            atomG.create_dataset("Coordinates", data=LCAO.Atom_au)
            atomG.create_dataset("UnitCell", data=LCAO.AtomCell_au)

    print("Export finished")
