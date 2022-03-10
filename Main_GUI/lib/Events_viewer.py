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

def openFile(win, Disp, Elements):
    currentFile=win.filePath.text()
    selectedFile, _filter=QtGui.QFileDialog.getOpenFileName(caption="Open file", directory=currentFile)
    if selectedFile!="":
        # valid
        win.filePath.setText(selectedFile)
        Disp.open(selectedFile)
        print("Finished loading dispersion file")
        
        win.kxIndex.setMaximum(Disp.Size[0]-1)
        if Disp.Dimension==1:
            win.eIndex.setMaximum(Disp.Size[1]-1)
        elif Disp.Dimension==2:
            win.kyIndex.setMaximum(Disp.Size[1]-1)
            win.eIndex.setMaximum(Disp.Size[2]-1)
        else:
            print("Dimension error")
            return
        
        if Disp.Dimension==1:
            plot(win, Disp)
            win.graphTab.setCurrentIndex(0)
        elif Disp.Dimension==2:
            plot3(win, Disp)
            win.graphTab.setCurrentIndex(1)

        makeRealSpace(win, Disp, Elements)

def plot(win, Disp):
    tr=QtGui.QTransform()
    tr.translate(Disp.Offset[0]-Disp.Delta[0]/2,Disp.Offset[1]-Disp.Delta[1]/2)
    tr.scale(Disp.Delta[0], Disp.Delta[1])
    win.img.setTransform(tr)
    win.img.setImage(Disp.Dispersion)

    maxpoint=Disp.Dispersion.max()
    print(("Max point: {0:8.3e}").format(maxpoint))
    win.bar.setLevels((0, maxpoint))

def plot3(win, Disp):

    tr_x=QtGui.QTransform()
    tr_x.translate(Disp.Offset[0]-Disp.Delta[0]/2,Disp.Offset[2]-Disp.Delta[2]/2)
    tr_x.scale(Disp.Delta[0], Disp.Delta[2])
    win.imgEx.setTransform(tr_x)

    tr_y=QtGui.QTransform()
    tr_y.translate(Disp.Offset[2]-Disp.Delta[2]/2, Disp.Offset[1]-Disp.Delta[1]/2)
    tr_y.rotate(-90)
    tr_y.scale(-Disp.Delta[1], Disp.Delta[2])
    win.imgEy.setTransform(tr_y)

    tr_E=QtGui.QTransform()
    tr_E.translate(Disp.Offset[0]-Disp.Delta[0]/2, Disp.Offset[1]-Disp.Delta[1]/2)
    tr_E.scale(Disp.Delta[0], Disp.Delta[1])
    win.imgxy.setTransform(tr_E)

    kx=win.kxIndex.value()
    ky=win.kyIndex.value()
    ei=win.eIndex.value()
    win.imgEx.setImage(Disp.Dispersion[:,ky,:])
    win.imgEy.setImage(Disp.Dispersion[kx,:,:])
    win.imgxy.setImage(Disp.Dispersion[:,:,ei])

    Maxpoint=Disp.Dispersion.max()
    win.Cube=np.zeros((Disp.Dispersion.shape[0], Disp.Dispersion.shape[1], Disp.Dispersion.shape[2], 4))
    win.Cube[:,:,:,0]=255
    win.Cube[:,:,:,1]=255
    win.Cube[:,:,:,2]=255
    win.Cube[:,:,:,3]=Disp.Dispersion/Maxpoint*100

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

    win.bandCube=gl.GLVolumeItem(win.Cube)
    win.bandCube.setData(win.Cube)    
    win.bandCube.scale(Disp.Delta[0], Disp.Delta[1], Disp.Delta[2])    
    win.bandCube.translate(Disp.Offset[0]-Disp.Delta[0]/2, Disp.Offset[1]-Disp.Delta[1]/2, Disp.Offset[2]-Disp.Delta[2]/2)
    win.plot3D.clear()
    win.plot3D.addItem(win.bandCube)

    drawCursor3(win, Disp)
            
            
def drawCursor3(win, Disp):
    
    kx=win.kxIndex.value()
    ky=win.kyIndex.value()
    ei=win.eIndex.value()

    kx_value=Disp.Offset[0]+Disp.Delta[0]*kx
    ky_value=Disp.Offset[1]+Disp.Delta[1]*ky
    e_value=Disp.Offset[2]+Disp.Delta[2]*ei
    
    win.vLineEx.setPos(kx_value)
    win.hLineEx.setPos(e_value)
    win.vLineEy.setPos(e_value)
    win.hLineEy.setPos(ky_value)
    win.vLinexy.setPos(kx_value)
    win.hLinexy.setPos(ky_value)
    win.kxValue.setText(("({0:.3f})").format(kx_value))
    win.kyValue.setText(("({0:.3f})").format(ky_value))
    win.eValue.setText(("({0:.3f})").format(e_value))
    
    
# real-space image
# in unit of Ang, so conversion is necessary
def makeRealSpace(win, Disp, Elements):
    # ia: atom
    # ie: element
    boundaryA=win.boundaryA.value()
    boundaryB=win.boundaryB.value()
    boundaryC=win.boundaryC.value()

    enableWeighting=win.enableWeight.isChecked()
    
    win.realSpace.clear()
    # print(LCAO.atomCell_au)
    
    au_ang=Config.au_ang

    # unit cell
    pts=np.zeros((5,3))
    # base plane
    pts[1]=Disp.AtomCell[0]*au_ang
    pts[2]=(Disp.AtomCell[0]+Disp.AtomCell[1])*au_ang
    pts[3]=Disp.AtomCell[1]*au_ang
    
    plt=gl.GLLinePlotItem(pos=pts)
    win.realSpace.addItem(plt)
    
    #top plane
    pts2=pts.copy()
    for i in range(0, 5):
        pts2[i]+=Disp.AtomCell[2]*au_ang
    plt=gl.GLLinePlotItem(pos=pts2)
    win.realSpace.addItem(plt)

    # vertical pillars
    for i in range(0, 4):
        pts3=np.zeros((2, 3))
        for j in range(0, 3):
            pts3[0][j]=pts[i][j]
            pts3[1][j]=pts[i][j]+Disp.AtomCell[2][j]*au_ang
        plt=gl.GLLinePlotItem(pos=pts3)
        win.realSpace.addItem(plt)

    # kx and ky vector (direction)
    kx=np.zeros((2,3))
    kx[1]=Disp.Xvector*Config.reciprocal_coeff/au_ang
    kx_color=np.zeros((2,4))
    kx_color[:,0]=Config.pen_kx[0]
    kx_color[:,1]=Config.pen_kx[1]
    kx_color[:,2]=Config.pen_kx[2]
    kx_color[:,3]=1.0
    plt=gl.GLLinePlotItem(pos=kx, color=kx_color, width=Config.reciprocal_axis_width)
    win.realSpace.addItem(plt)
    if Disp.Dimension==2:
        ky=np.zeros((2,3))
        ky[1]=Disp.Yvector*Config.reciprocal_coeff/au_ang
        ky_color=np.zeros((2,4))
        ky_color[:,0]=Config.pen_ky[0]
        ky_color[:,1]=Config.pen_ky[1]
        ky_color[:,2]=Config.pen_ky[2]
        ky_color[:,3]=1.0
        plt=gl.GLLinePlotItem(pos=ky, color=ky_color, width=Config.reciprocal_axis_width)
        win.realSpace.addItem(plt)
    # print(Disp.Xvector)
    # print(Disp.Yvector)

    # Polarization
    pol=np.zeros((2,3))
    theta=math.radians(Disp.Theta)
    phi=math.radians(Disp.Phi)
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
    for ia, atom_label in enumerate(Disp.Atom_label):
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
        if enableWeighting and Disp.Weighting:
            meshcolor[:,3]=Disp.Atom_weighting[ia]*(1.0-Config.Weighting_offset)+Config.Weighting_offset
            # print(Disp.Atom_weighting[ia]*(1.0-Config.Weighting_offset)+Config.Weighting_offset)
        # print(meshcolor)

        md.setFaceColors(meshcolor)

        for iA in range(0, boundaryA):
            for iB in range(0, boundaryB):
                for iC in range(0, boundaryC):
                    mi=gl.GLMeshItem(meshdata=md, smooth=False)
                    coordinate=Disp.Atom_coord[ia].copy()
                    coordinate+=iA*Disp.AtomCell[0]
                    coordinate+=iB*Disp.AtomCell[1]
                    coordinate+=iC*Disp.AtomCell[2]
                    
                    mi.translate(coordinate[0]*au_ang, coordinate[1]*au_ang, coordinate[2]*au_ang)
                    mi.setGLOptions("additive")
                    win.realSpace.addItem(mi)

