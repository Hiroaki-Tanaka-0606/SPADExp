# Events
# clicked, index changed etc.

import Config
import PAO_parser as PAOp
import AO_parser as AOp
import os
import re
import subprocess
import shutil
from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import numpy as np
import h5py
from datetime import datetime

def changeAtom(win):
    # QComboBox atomNumber is changed

    Z=win.atomNumber.currentIndex()

    win.paoFile.clear()
    for paoName in Config.paoArr[Z]:
        win.paoFile.addItem(paoName)

    win.vpsFile.clear()
    for vpsName in Config.vpsArr[Z]:
        win.vpsFile.addItem(vpsName)

    changeAnalysis(win)

def clearAll(win):
    win.orbitalTable.setRowCount(0)
    win.orbitalTable.setColumnCount(0)
    win.matrix.setRowCount(0)
    win.matrix.setColumnCount(0)
    win.orbitalGraph.clear()
    
def changeAnalysis(win):
    # QComboBox analysisType is changed

    global PAO_after
    global PAO_before
    global PAO_fromVPS
    global PAO_reproduced
    global AO_fromVPS
    global Calc_CCoes
    global Selected_orbitals
    analysisIndex=win.analysisType.currentIndex()
    Z=win.atomNumber.currentIndex()
    paoIndex=win.paoFile.currentIndex()
    vpsIndex=win.vpsFile.currentIndex()
    if 0<=paoIndex and paoIndex<len(Config.paoArr[Z]):
        paoName=Config.paoArr[Z][paoIndex]
    else:
        return
        
    if 0<=vpsIndex and vpsIndex<len(Config.vpsArr[Z]):
        vpsName=Config.vpsArr[Z][vpsIndex]
    else:
        return
    
    if analysisIndex==1: # PAO: Before and after optimization
        win.matrixType.clear()
        # check PAO files (before and after optimization) exist
        paoPath_after=os.path.join(Config.dirPath_PAO, paoName)
        paoPath_before=os.path.join(Config.dirPath_PAOfromPAO, paoName)

        if os.path.isfile(paoPath_after):
            print(("OK: {0:s} after optimization exists").format(paoName))
        else:
            print(("Not OK: {0:s} after optimization does not exist").format(paoName))
            clearAll(win)
            return
            
        if os.path.isfile(paoPath_before):
            print(("OK: {0:s} before optimization exists").format(paoName))
        else:
            print(("Not OK: {0:s} before optimization does not exist").format(paoName))
            clearAll(win)
            return

        print(("Read {0:s}").format(paoPath_after))
        PAO_after=PAOp.PAO_parser(paoPath_after, True)
        print(("Read {0:s}").format(paoPath_before))
        PAO_before=PAOp.PAO_parser(paoPath_before, False)
        Calc_CCoes=np.zeros((PAO_after.PAO_Lmax+1, PAO_after.PAO_Mul, PAO_after.PAO_Mul))
        PAOp.calcContraction(PAO_after, PAO_before, Calc_CCoes)
        
        Selected_orbitals=np.zeros((PAO_after.PAO_Lmax+1, PAO_after.PAO_Mul))

        # construct orbital table
        win.orbitalTable.setRowCount(PAO_after.PAO_Lmax+1)
        win.orbitalTable.setColumnCount(PAO_after.PAO_Mul)
        for l in range(0, PAO_after.PAO_Lmax+1):
            win.orbitalTable.setVerticalHeaderItem(l, QtGui.QTableWidgetItem(Config.azimuthal[l]))
            
        for p in range(0, PAO_after.PAO_Mul):
            for l in range(0, PAO_after.PAO_Lmax+1):
                item=QtGui.QTableWidgetItem(Config.azimuthal[l]+str(p))
                item.setBackground(Config.unselected_item)
                win.orbitalTable.setItem(l, p, item)

        # prepare matrixType combobox
        for l in range(0, PAO_after.PAO_Lmax+1):
            win.matrixType.addItem(("Contraction coefficients for L= {0:d}").format(l))
        for l in range(0, PAO_after.PAO_Lmax+1):
            win.matrixType.addItem(("Orthonormalization matrix (before) for L= {0:d}").format(l))
        for l in range(0, PAO_after.PAO_Lmax+1):
            win.matrixType.addItem(("Orthonormalization matrix (after) for L= {0:d}").format(l))

    elif analysisIndex==2: # PAO: Before optimization and from VPS
        win.matrixType.clear()
        # check PAO files (before optimization and from VPS) exist
        paoPath_before=os.path.join(Config.dirPath_PAOfromPAO, paoName)
        outPaoName=re.sub(r"\.vps$", ".pao", vpsName)
        paoPath_fromVPS=os.path.join(Config.dirPath_PAOfromVPS, outPaoName)

        if os.path.isfile(paoPath_before):
            print(("OK: {0:s} before optimization exists").format(paoName))
        else:
            print(("Not OK: {0:s} before optimization does not exist").format(paoName))
            clearAll(win)
            return
            
        if os.path.isfile(paoPath_fromVPS):
            print(("OK: {0:s} from VPS exists").format(outPaoName))
        else:
            print(("Not OK: {0:s} from VPS does not exist").format(outPaoName))
            clearAll(win)
            return

        print(("Read {0:s}").format(paoPath_before))
        PAO_before=PAOp.PAO_parser(paoPath_before, True)
        print(("Read {0:s}").format(paoPath_fromVPS))
        PAO_fromVPS=PAOp.PAO_parser(paoPath_fromVPS, False)

        Calc_CCoes=np.zeros((PAO_before.PAO_Lmax+1, PAO_before.PAO_Mul, PAO_fromVPS.PAO_Mul))
        PAO_reproduced=np.zeros((PAO_before.PAO_Lmax+1, PAO_before.PAO_Mul, PAO_before.grid_output))

        # compatibility check
        if PAO_before.grid_output != PAO_fromVPS.grid_output or\
           abs(PAO_before.x[0]-PAO_fromVPS.x[0])>0.01 or\
           abs(PAO_before.x[-1]-PAO_fromVPS.x[-1])>0.01:
            print("Error: grid mismatch")
            return
        
        
        PAOp.calcContraction(PAO_before, PAO_fromVPS, Calc_CCoes)
        Lmax=max(PAO_before.PAO_Lmax, PAO_fromVPS.PAO_Lmax)
        Mul=max(PAO_before.PAO_Mul, PAO_fromVPS.PAO_Mul)
        Selected_orbitals=np.zeros((Lmax+1, Mul))

        PAOp.reproducePAO(PAO_fromVPS, Calc_CCoes, PAO_reproduced)

        # construct orbital table
        win.orbitalTable.setRowCount(Lmax+1)
        win.orbitalTable.setColumnCount(Mul)
        for l in range(0, Lmax+1):
            win.orbitalTable.setVerticalHeaderItem(l, QtGui.QTableWidgetItem(Config.azimuthal[l]))
            
        for p in range(0, Mul):
            for l in range(0, Lmax+1):
                item=QtGui.QTableWidgetItem(Config.azimuthal[l]+str(p))
                item.setBackground(Config.unselected_item)
                win.orbitalTable.setItem(l, p, item)

        # prepare matrixType combobox
        for l in range(0, PAO_before.PAO_Lmax+1):
            win.matrixType.addItem(("Contraction coefficients for L= {0:d}").format(l))
        for l in range(0, PAO_before.PAO_Lmax+1):
            win.matrixType.addItem(("Orthonormalization matrix (before) for L= {0:d}").format(l))
        for l in range(0, PAO_fromVPS.PAO_Lmax+1):
            win.matrixType.addItem(("Orthonormalization matrix (from VPS) for L= {0:d}").format(l))

    elif analysisIndex==3: # PAO: After optimization and from VPS
        win.matrixType.clear()
        # check PAO files (after optimization and from VPS) exist
        paoPath_after=os.path.join(Config.dirPath_PAO, paoName)
        outPaoName=re.sub(r"\.vps$", ".pao", vpsName)
        paoPath_fromVPS=os.path.join(Config.dirPath_PAOfromVPS, outPaoName)

        if os.path.isfile(paoPath_after):
            print(("OK: {0:s} after optimization exists").format(paoName))
        else:
            print(("Not OK: {0:s} after optimization does not exist").format(paoName))
            clearAll(win)
            return
            
        if os.path.isfile(paoPath_fromVPS):
            print(("OK: {0:s} from VPS exists").format(outPaoName))
        else:
            print(("Not OK: {0:s} from VPS does not exist").format(outPaoName))
            clearAll(win)
            return

        print(("Read {0:s}").format(paoPath_after))
        PAO_after=PAOp.PAO_parser(paoPath_after, True)
        print(("Read {0:s}").format(paoPath_fromVPS))
        PAO_fromVPS=PAOp.PAO_parser(paoPath_fromVPS, False)

        Calc_CCoes=np.zeros((PAO_after.PAO_Lmax+1, PAO_after.PAO_Mul, PAO_fromVPS.PAO_Mul))
        PAO_reproduced=np.zeros((PAO_after.PAO_Lmax+1, PAO_after.PAO_Mul, PAO_after.grid_output))

        # compatibility check
        if PAO_after.grid_output != PAO_fromVPS.grid_output or\
           abs(PAO_after.x[0]-PAO_fromVPS.x[0])>0.01 or\
           abs(PAO_after.x[-1]-PAO_fromVPS.x[-1])>0.01:
            print("Error: grid mismatch")
            return
        
        
        PAOp.calcContraction(PAO_after, PAO_fromVPS, Calc_CCoes)
        Lmax=max(PAO_after.PAO_Lmax, PAO_fromVPS.PAO_Lmax)
        Mul=max(PAO_after.PAO_Mul, PAO_fromVPS.PAO_Mul)
        Selected_orbitals=np.zeros((Lmax+1, Mul))

        PAOp.reproducePAO(PAO_fromVPS, Calc_CCoes, PAO_reproduced)

        # construct orbital table
        win.orbitalTable.setRowCount(Lmax+1)
        win.orbitalTable.setColumnCount(Mul)
        for l in range(0, Lmax+1):
            win.orbitalTable.setVerticalHeaderItem(l, QtGui.QTableWidgetItem(Config.azimuthal[l]))
            
        for p in range(0, Mul):
            for l in range(0, Lmax+1):
                item=QtGui.QTableWidgetItem(Config.azimuthal[l]+str(p))
                item.setBackground(Config.unselected_item)
                win.orbitalTable.setItem(l, p, item)

        # prepare matrixType combobox
        for l in range(0, PAO_after.PAO_Lmax+1):
            win.matrixType.addItem(("Contraction coefficients for L= {0:d}").format(l))
        for l in range(0, PAO_after.PAO_Lmax+1):
            win.matrixType.addItem(("Orthonormalization matrix (after) for L= {0:d}").format(l))
        for l in range(0, PAO_fromVPS.PAO_Lmax+1):
            win.matrixType.addItem(("Orthonormalization matrix (from VPS) for L= {0:d}").format(l))

    elif analysisIndex==4: # PAO and AO from VPS
        win.matrixType.clear()
        # check PAO file from VPS exists
        outPaoName=re.sub(r"\.vps$", ".pao", vpsName)
        paoPath_fromVPS=os.path.join(Config.dirPath_PAOfromVPS, outPaoName)
        # check AO file from VPS exists
        outAoName=re.sub(r"\.vps$", ".ao", vpsName)
        aoPath_fromVPS=os.path.join(Config.dirPath_AOfromVPS, outAoName)

        if os.path.isfile(paoPath_fromVPS):
            print(("OK: {0:s} from VPS exists").format(outPaoName))
        else:
            print(("Not OK: {0:s} from VPS does not exist").format(outPaoName))
            clearAll(win)
            return
            
        if os.path.isfile(aoPath_fromVPS):
            print(("OK: {0:s} from VPS exists").format(outAoName))
        else:
            print(("Not OK: {0:s} from VPS does not exist").format(outAoName))
            clearAll(win)
            return

        print(("Read {0:s}").format(paoPath_fromVPS))
        PAO_fromVPS=PAOp.PAO_parser(paoPath_fromVPS, False)
        print(("Read {0:s}").format(aoPath_fromVPS))
        AO_fromVPS=AOp.AO_parser(aoPath_fromVPS)

        # compatibility check
        if PAO_fromVPS.grid_output != AO_fromVPS.grid_output or\
           abs(PAO_fromVPS.x[0]-AO_fromVPS.x[0])>0.01 or\
           abs(PAO_fromVPS.x[-1]-AO_fromVPS.x[-1])>0.01:
            print("Error: grid mismatch")
            return

        PAOp.normCheck(PAO_fromVPS, AO_fromVPS)
        
        # construct orbital table
        win.orbitalTable.setRowCount(AO_fromVPS.maxL_pao+1)
        win.orbitalTable.setColumnCount(AO_fromVPS.max_N)
        Selected_orbitals=np.zeros((AO_fromVPS.maxL_pao+1, AO_fromVPS.max_N))
        for l in range(0, AO_fromVPS.maxL_pao+1):
            win.orbitalTable.setVerticalHeaderItem(l, QtGui.QTableWidgetItem(Config.azimuthal[l]))

            
        for n in range(1, AO_fromVPS.max_N+1):
            for l in range(0, AO_fromVPS.maxL_pao+1):
                item_str=str(n)+Config.azimuthal[l]
                if l>=n:
                    item_str=""
                elif n<(PAO_fromVPS.valence_min[l] if len(PAO_fromVPS.valence_min)>l else l+1):
                    item_str="("+item_str+")"
                elif n>=(PAO_fromVPS.valence_min[l] if len(PAO_fromVPS.valence_min)>l else l+1)+PAO_fromVPS.PAO_Mul:
                    item_str="["+ item_str+"]"
                    
                item=QtGui.QTableWidgetItem(item_str)
                item.setBackground(Config.unselected_item)
                win.orbitalTable.setItem(l, n-1, item)

        # prepare matrixType combobox
        for l in range(0, PAO_fromVPS.PAO_Lmax+1):
            win.matrixType.addItem(("Orthonormalization matrix (PAO) for L= {0:d}").format(l))
        for l in range(0, AO_fromVPS.maxL_pao+1):
            win.matrixType.addItem(("Orthonormalization matrix (AO) for L= {0:d}").format(l))

    elif analysisIndex==5: # PAO and AO after optimization
        win.matrixType.clear()
        # check PAO file after optimization exists
        paoPath_after=os.path.join(Config.dirPath_PAO, paoName)
        # check PAO file from VPS exists
        outPaoName=re.sub(r"\.vps$", ".pao", vpsName)
        paoPath_fromVPS=os.path.join(Config.dirPath_PAOfromVPS, outPaoName)
        # check AO file from VPS exists
        outAoName=re.sub(r"\.vps$", ".ao", vpsName)
        aoPath_fromVPS=os.path.join(Config.dirPath_AOfromVPS, outAoName)

        if os.path.isfile(paoPath_after):
            print(("OK: {0:s} after optimization exists").format(paoName))
        else:
            print(("Not OK: {0:s} after optimization does not exist").format(paoName))
            clearAll(win)
            return
            
        if os.path.isfile(paoPath_fromVPS):
            print(("OK: {0:s} from VPS exists").format(outPaoName))
        else:
            print(("Not OK: {0:s} from VPS does not exist").format(outPaoName))
            clearAll(win)
            return
            
        if os.path.isfile(aoPath_fromVPS):
            print(("OK: {0:s} from VPS exists").format(outAoName))
        else:
            print(("Not OK: {0:s} from VPS does not exist").format(outAoName))
            clearAll(win)
            return

        print(("Read {0:s}").format(paoPath_after))
        PAO_after=PAOp.PAO_parser(paoPath_after, True)
        print(("Read {0:s}").format(paoPath_fromVPS))
        PAO_fromVPS=PAOp.PAO_parser(paoPath_fromVPS, False)
        print(("Read {0:s}").format(aoPath_fromVPS))
        AO_fromVPS=AOp.AO_parser(aoPath_fromVPS)

        # compatibility check
        if PAO_after.grid_output != PAO_fromVPS.grid_output or\
           abs(PAO_after.x[0]-PAO_fromVPS.x[0])>0.01 or\
           abs(PAO_after.x[-1]-PAO_fromVPS.x[-1])>0.01 or\
           PAO_fromVPS.grid_output != AO_fromVPS.grid_output or\
           abs(PAO_fromVPS.x[0]-AO_fromVPS.x[0])>0.01 or\
           abs(PAO_fromVPS.x[-1]-AO_fromVPS.x[-1])>0.01:
            print("Error: grid mismatch")
            return

        Calc_CCoes=np.zeros((PAO_after.PAO_Lmax+1, PAO_after.PAO_Mul, PAO_fromVPS.PAO_Mul))
        PAO_reproduced=np.zeros((PAO_after.PAO_Lmax+1, PAO_after.PAO_Mul, PAO_after.grid_output))
        PAOp.normCheck(PAO_fromVPS, AO_fromVPS)

        PAOp.calcContraction(PAO_after, PAO_fromVPS, Calc_CCoes)
        Selected_orbitals=np.zeros((PAO_after.PAO_Lmax, PAO_after.PAO_Mul))

        PAOp.reproduceAO(AO_fromVPS, PAO_fromVPS, Calc_CCoes, PAO_reproduced)
        
        # construct orbital table
        win.orbitalTable.setRowCount(PAO_after.PAO_Lmax)
        win.orbitalTable.setColumnCount(PAO_after.PAO_Mul)
        
        for l in range(0, PAO_fromVPS.PAO_Lmax+1):
            win.orbitalTable.setVerticalHeaderItem(l, QtGui.QTableWidgetItem(Config.azimuthal[l]))

            
        for p in range(0, PAO_after.PAO_Mul):
            for l in range(0, PAO_after.PAO_Lmax+1):
                item_str=Config.azimuthal[l]+str(p)
                item=QtGui.QTableWidgetItem(item_str)
                item.setBackground(Config.unselected_item)
                win.orbitalTable.setItem(l, p, item)

def selectOrbital(win, row, column):
    global Selected_orbitals
    if Selected_orbitals[row][column]==0:
        Selected_orbitals[row][column]=1
        win.orbitalTable.item(row, column).setBackground(Config.selected_item)
    else:
        Selected_orbitals[row][column]=0
        win.orbitalTable.item(row, column).setBackground(Config.unselected_item)
    drawOrbitalGraph(win)

def drawOrbitalGraph(win):
    global Selected_orbitals
    global PAO_after
    global PAO_before
    global PAO_fromVPS
    global AO_fromVPS
    
    radialType=win.radialType.checkedButton().text()[0]
    analysisIndex=win.analysisType.currentIndex()
    if analysisIndex==1: # PAO: Before and after optimization
        win.orbitalGraph.clear()
        for l in range(0, len(Selected_orbitals)):
            for mul in range(0, len(Selected_orbitals[l])):
                if Selected_orbitals[l][mul]==1:
                    if radialType=="P":
                        win.orbitalGraph.plot(y=PAO_after.PAOs[l][mul]*PAO_after.r,\
                                              x=PAO_after.r,\
                                              pen=Config.pen_after,\
                                              name=Config.azimuthal[l]+str(mul)+" after")
                        win.orbitalGraph.plot(y=PAO_before.PAOs[l][mul]*PAO_before.r,\
                                              x=PAO_before.r, \
                                              pen=Config.pen_before,\
                                              name=Config.azimuthal[l]+str(mul)+" before")
                    else:
                        win.orbitalGraph.plot(y=PAO_after.PAOs[l][mul], \
                                              x=PAO_after.r, \
                                              pen=Config.pen_after, \
                                              name=Config.azimuthal[l]+str(mul)+" after")
                        win.orbitalGraph.plot(y=PAO_before.PAOs[l][mul], \
                                              x=PAO_before.r, \
                                              pen=Config.pen_before, \
                                              name=Config.azimuthal[l]+str(mul)+" before")
    elif analysisIndex==2: # PAO: Before optimization and from VPS
        win.orbitalGraph.clear()
        for l in range(0, len(Selected_orbitals)):
            for mul in range(0, len(Selected_orbitals[l])):
                if Selected_orbitals[l][mul]==1:
                    if radialType=="P":
                        if l<len(PAO_before.PAOs) and mul<len(PAO_before.PAOs[l]):
                            win.orbitalGraph.plot(y=PAO_before.PAOs[l][mul]*PAO_before.r,\
                                                  x=PAO_before.r,\
                                                  pen=Config.pen_before,\
                                                  name=Config.azimuthal[l]+str(mul)+" before")
                            win.orbitalGraph.plot(y=PAO_reproduced[l][mul]*PAO_before.r,\
                                                  x=PAO_before.r,\
                                                  pen=pg.mkPen(color=Config.pen_reproduced,style=QtCore.Qt.DashLine),\
                                                  name=Config.azimuthal[l]+str(mul)+" before (reproduced)")
                        if l<len(PAO_fromVPS.PAOs) and mul<len(PAO_fromVPS.PAOs[l]):
                            win.orbitalGraph.plot(y=PAO_fromVPS.PAOs[l][mul]*PAO_fromVPS.r,\
                                                  x=PAO_fromVPS.r, \
                                                  pen=Config.pen_fromVPS,\
                                                  name=Config.azimuthal[l]+str(mul)+" from VPS")
                    else:
                        if l<len(PAO_before.PAOs) and mul<len(PAO_before.PAOs[l]):
                            win.orbitalGraph.plot(y=PAO_before.PAOs[l][mul],\
                                                  x=PAO_before.r,\
                                                  pen=Config.pen_before,\
                                                  name=Config.azimuthal[l]+str(mul)+" before")
                            win.orbitalGraph.plot(y=PAO_reproduced[l][mul],\
                                                  x=PAO_before.r,\
                                                  pen=pg.mkPen(color=Config.pen_reproduced,style=QtCore.Qt.DashLine),\
                                                  name=Config.azimuthal[l]+str(mul)+" before (reproduced)")
                        if l<len(PAO_fromVPS.PAOs) and mul<len(PAO_fromVPS.PAOs[l]):
                            win.orbitalGraph.plot(y=PAO_fromVPS.PAOs[l][mul],\
                                                  x=PAO_fromVPS.r, \
                                                  pen=Config.pen_fromVPS,\
                                                  name=Config.azimuthal[l]+str(mul)+" from VPS")
    elif analysisIndex==3: # PAO: After optimization and from VPS
        win.orbitalGraph.clear()
        for l in range(0, len(Selected_orbitals)):
            for mul in range(0, len(Selected_orbitals[l])):
                if Selected_orbitals[l][mul]==1:
                    # print(PAO_after.PAOs[l][mul]-PAO_reproduced[l][mul])
                    if radialType=="P":
                        if l<len(PAO_after.PAOs) and mul<len(PAO_after.PAOs[l]):
                            win.orbitalGraph.plot(y=PAO_after.PAOs[l][mul]*PAO_after.r,\
                                                  x=PAO_after.r,\
                                                  pen=Config.pen_after,\
                                                  name=Config.azimuthal[l]+str(mul)+" after")
                            win.orbitalGraph.plot(y=PAO_reproduced[l][mul]*PAO_after.r,\
                                                  x=PAO_after.r,\
                                                  pen=pg.mkPen(color=Config.pen_reproduced,style=QtCore.Qt.DashLine),\
                                                  name=Config.azimuthal[l]+str(mul)+" after (reproduced)")
                        if l<len(PAO_fromVPS.PAOs) and mul<len(PAO_fromVPS.PAOs[l]):
                            win.orbitalGraph.plot(y=PAO_fromVPS.PAOs[l][mul]*PAO_fromVPS.r,\
                                                  x=PAO_fromVPS.r, \
                                                  pen=Config.pen_fromVPS,\
                                                  name=Config.azimuthal[l]+str(mul)+" from VPS")
                    else:
                        if l<len(PAO_after.PAOs) and mul<len(PAO_after.PAOs[l]):
                            win.orbitalGraph.plot(y=PAO_after.PAOs[l][mul],\
                                                  x=PAO_after.r,\
                                                  pen=Config.pen_after,\
                                                  name=Config.azimuthal[l]+str(mul)+" after")
                            win.orbitalGraph.plot(y=PAO_reproduced[l][mul],\
                                                  x=PAO_after.r,\
                                                  pen=pg.mkPen(color=Config.pen_reproduced,style=QtCore.Qt.DashLine),\
                                                  name=Config.azimuthal[l]+str(mul)+" after (reproduced)")
                        if l<len(PAO_fromVPS.PAOs) and mul<len(PAO_fromVPS.PAOs[l]):
                            win.orbitalGraph.plot(y=PAO_fromVPS.PAOs[l][mul],\
                                                  x=PAO_fromVPS.r, \
                                                  pen=Config.pen_fromVPS,\
                                                  name=Config.azimuthal[l]+str(mul)+" from VPS")

    elif analysisIndex==4: # PAO and AO from VPS
        win.orbitalGraph.clear()
        for l in range(0, len(Selected_orbitals)):
            for mul in range(0, len(Selected_orbitals[l])):
                if Selected_orbitals[l][mul]==1 and l<(mul+1):
                    mul_PAO=mul+1-(PAO_fromVPS.valence_min[l] if len(PAO_fromVPS.valence_min)>l else l+1)
                    mul_AO=mul+1-(l+1)
                    if radialType=="P":
                        if mul_PAO<PAO_fromVPS.PAO_Mul and mul_PAO>=0:
                            win.orbitalGraph.plot(y=PAO_fromVPS.PAOs[l][mul_PAO]*PAO_fromVPS.r,\
                                                  x=PAO_fromVPS.r,\
                                                  pen=Config.pen_fromVPS,\
                                                  name=str(mul+1)+Config.azimuthal[l]+" ("+str(mul_PAO+1)+Config.azimuthal[l]+" in PAO)")
                        if mul_AO>=0:
                            win.orbitalGraph.plot(y=AO_fromVPS.AOs[l][mul_AO]*AO_fromVPS.r,\
                                                  x=AO_fromVPS.r, \
                                                  pen=pg.mkPen(color=Config.pen_fromVPS, style=QtCore.Qt.DashLine),\
                                                  name=str(mul+1)+Config.azimuthal[l]+" AO")
                    else:
                        if mul_PAO<PAO_fromVPS.PAO_Mul and mul_PAO>=0:
                            win.orbitalGraph.plot(y=PAO_fromVPS.PAOs[l][mul_PAO],\
                                                  x=PAO_fromVPS.r,\
                                                  pen=Config.pen_fromVPS,\
                                                  name=str(mul+1)+Config.azimuthal[l]+" ("+str(mul_PAO+1)+Config.azimuthal[l]+" in PAO)")
                        if mul_AO>=0:
                            win.orbitalGraph.plot(y=AO_fromVPS.AOs[l][mul_AO],\
                                                  x=AO_fromVPS.r, \
                                                  pen=pg.mkPen(color=Config.pen_fromVPS, style=QtCore.Qt.DashLine),\
                                                  name=str(mul+1)+Config.azimuthal[l]+" AO")

    elif analysisIndex==5: # PAO and AO after optimization
        win.orbitalGraph.clear()
        for l in range(0, len(Selected_orbitals)):
            for mul in range(0, len(Selected_orbitals[l])):
                if Selected_orbitals[l][mul]==1:
                    if radialType=="P":
                        win.orbitalGraph.plot(y=PAO_after.PAOs[l][mul]*PAO_after.r,\
                                              x=PAO_after.r,\
                                              pen=Config.pen_after,\
                                              name=Config.azimuthal[l]+str(mul)+" PAO")
        
                        win.orbitalGraph.plot(y=PAO_reproduced[l][mul]*PAO_after.r,\
                                              x=PAO_after.r, \
                                              pen=Config.pen_fromVPS,\
                                              name=Config.azimuthal[l]+str(mul)+" AO")
                    else:
                        win.orbitalGraph.plot(y=PAO_after.PAOs[l][mul],\
                                              x=PAO_after.r,\
                                              pen=Config.pen_after,\
                                              name=Config.azimuthal[l]+str(mul)+" PAO")
        
                        win.orbitalGraph.plot(y=PAO_reproduced[l][mul],\
                                              x=PAO_after.r, \
                                              pen=Config.pen_fromVPS,\
                                              name=Config.azimuthal[l]+str(mul)+" AO")
        
def changeMatrix(win):
    # QComboBox matrixType is changed

    analysisIndex=win.analysisType.currentIndex()
    matrixIndex=win.matrixType.currentIndex()

    global PAO_after
    global PAO_before
    global Calc_CCoes
    
    if analysisIndex==1: # PAO: Before and after optimization
        if not ("PAO_after" in globals()) or not ("PAO_before" in globals()):
            return
        win.matrix.setRowCount(PAO_after.PAO_Mul)
        win.matrix.setColumnCount(PAO_after.PAO_Mul)
        if matrixIndex>=0 and matrixIndex<=PAO_after.PAO_Lmax:
            # Contraction coefficients
            l=matrixIndex
            for i in range(0, PAO_after.PAO_Mul):
                for j in range(0, PAO_after.PAO_Mul):
                    item=QtGui.QTableWidgetItem(("{0:7.4f}\n{1:7.4f}").format(PAO_after.CCoes[l][i][j], Calc_CCoes[l][i][j]))
                    win.matrix.setItem(j, i, item)
                    error=PAO_after.CCoes[l][i][j]-Calc_CCoes[l][i][j]
                    if abs(error)>Config.Contraction_threshold:
                        item.setBackground(Config.error_item)
        elif matrixIndex>=PAO_after.PAO_Lmax+1 and matrixIndex<=PAO_after.PAO_Lmax*2+1:
            # Orthonormalization matrix (before)
            l=matrixIndex-PAO_after.PAO_Lmax-1
            for i in range(0, PAO_after.PAO_Mul):
                for j in range(0, PAO_after.PAO_Mul):
                    item=QtGui.QTableWidgetItem(("{0:7.4f}").format(PAO_before.OrthNorm[l][i][j]))
                    win.matrix.setItem(j, i, item)
                    error=PAO_before.OrthNorm[l][i][j]
                    if i==j:
                        error-=1
                    if abs(error)>Config.OrthNorm_threshold:
                        item.setBackground(Config.error_item)
        elif matrixIndex>=PAO_after.PAO_Lmax*2+2 and matrixIndex<=PAO_after.PAO_Lmax*3+2:
            # Orthonormalization matrix (after)
            l=matrixIndex-PAO_after.PAO_Lmax*2-2
            for i in range(0, PAO_after.PAO_Mul):
                for j in range(0, PAO_after.PAO_Mul):
                    item=QtGui.QTableWidgetItem(("{0:7.4f}").format(PAO_after.OrthNorm[l][i][j]))
                    win.matrix.setItem(j, i, item)
                    error=PAO_after.OrthNorm[l][i][j]
                    if i==j:
                        error-=1
                    if abs(error)>Config.OrthNorm_threshold:
                        item.setBackground(Config.error_item)

    elif analysisIndex==2: # PAO: Before optimization and from VPS
        if not ("PAO_fromVPS" in globals()) or not ("PAO_before" in globals()):
            return
        if matrixIndex>=0 and matrixIndex<=PAO_before.PAO_Lmax:
            # Contraction coefficients
            l=matrixIndex
            win.matrix.setRowCount(PAO_fromVPS.PAO_Mul)
            win.matrix.setColumnCount(PAO_before.PAO_Mul)
            for i in range(0, PAO_before.PAO_Mul):
                for j in range(0, PAO_fromVPS.PAO_Mul):
                    item=QtGui.QTableWidgetItem(("{0:7.4f}").format(Calc_CCoes[l][i][j]))
                    win.matrix.setItem(j, i, item)
        elif matrixIndex>=PAO_before.PAO_Lmax+1 and matrixIndex<=PAO_before.PAO_Lmax*2+1:
            # Orthonormalization matrix (before)
            l=matrixIndex-PAO_before.PAO_Lmax-1
            win.matrix.setColumnCount(PAO_before.PAO_Mul)
            win.matrix.setRowCount(PAO_before.PAO_Mul)
            for i in range(0, PAO_before.PAO_Mul):
                for j in range(0, PAO_before.PAO_Mul):
                    item=QtGui.QTableWidgetItem(("{0:7.4f}").format(PAO_before.OrthNorm[l][i][j]))
                    win.matrix.setItem(j, i, item)
                    error=PAO_before.OrthNorm[l][i][j]
                    if i==j:
                        error-=1
                    if abs(error)>Config.OrthNorm_threshold:
                        item.setBackground(Config.error_item)
        elif matrixIndex>=PAO_before.PAO_Lmax*2+2 and matrixIndex<=PAO_before.PAO_Lmax*2+PAO_fromVPS.PAO_Lmax+2:
            # Orthonormalization matrix (from VPS)
            l=matrixIndex-PAO_before.PAO_Lmax*2-2
            win.matrix.setColumnCount(PAO_fromVPS.PAO_Mul)
            win.matrix.setRowCount(PAO_fromVPS.PAO_Mul)
            for i in range(0, PAO_fromVPS.PAO_Mul):
                for j in range(0, PAO_fromVPS.PAO_Mul):
                    item=QtGui.QTableWidgetItem(("{0:7.4f}").format(PAO_fromVPS.OrthNorm[l][i][j]))
                    win.matrix.setItem(j, i, item)
                    error=PAO_fromVPS.OrthNorm[l][i][j]
                    if i==j:
                        error-=1
                    if abs(error)>Config.OrthNorm_threshold:
                        item.setBackGround(Config.error_item)

    elif analysisIndex==3: # PAO: After optimization and from VPS
        if not ("PAO_fromVPS" in globals()) or not ("PAO_after" in globals()):
            return
        if matrixIndex>=0 and matrixIndex<=PAO_after.PAO_Lmax:
            # Contraction coefficients
            l=matrixIndex
            win.matrix.setRowCount(PAO_fromVPS.PAO_Mul)
            win.matrix.setColumnCount(PAO_after.PAO_Mul)
            for i in range(0, PAO_after.PAO_Mul):
                for j in range(0, PAO_fromVPS.PAO_Mul):
                    item=QtGui.QTableWidgetItem(("{0:7.4f}\n({1:7.4f})").format(Calc_CCoes[l][i][j],(PAO_after.CCoes[l][i][j] if j<PAO_after.PAO_Mul else 0)))
                    win.matrix.setItem(j, i, item)
                    error=(PAO_after.CCoes[l][i][j] if j<PAO_after.PAO_Mul else 0)-Calc_CCoes[l][i][j]
                    if abs(error)>Config.Contraction_threshold:
                        item.setBackground(Config.error_item)
        elif matrixIndex>=PAO_after.PAO_Lmax+1 and matrixIndex<=PAO_after.PAO_Lmax*2+1:
            # Orthonormalization matrix (after)
            l=matrixIndex-PAO_after.PAO_Lmax-1
            win.matrix.setColumnCount(PAO_after.PAO_Mul)
            win.matrix.setRowCount(PAO_after.PAO_Mul)
            for i in range(0, PAO_after.PAO_Mul):
                for j in range(0, PAO_after.PAO_Mul):
                    item=QtGui.QTableWidgetItem(("{0:7.4f}").format(PAO_after.OrthNorm[l][i][j]))
                    win.matrix.setItem(j, i, item)
                    error=PAO_after.OrthNorm[l][i][j]
                    if i==j:
                        error-=1
                    if abs(error)>Config.OrthNorm_threshold:
                        item.setBackground(Config.error_item)
        elif matrixIndex>=PAO_after.PAO_Lmax*2+2 and matrixIndex<=PAO_after.PAO_Lmax*2+PAO_fromVPS.PAO_Lmax+2:
            # Orthonormalization matrix (from VPS)
            l=matrixIndex-PAO_after.PAO_Lmax*2-2
            win.matrix.setColumnCount(PAO_fromVPS.PAO_Mul)
            win.matrix.setRowCount(PAO_fromVPS.PAO_Mul)
            for i in range(0, PAO_fromVPS.PAO_Mul):
                for j in range(0, PAO_fromVPS.PAO_Mul):
                    item=QtGui.QTableWidgetItem(("{0:7.4f}").format(PAO_fromVPS.OrthNorm[l][i][j]))
                    win.matrix.setItem(j, i, item)
                    error=PAO_fromVPS.OrthNorm[l][i][j]
                    if i==j:
                        error-=1
                    if abs(error)>Config.OrthNorm_threshold:
                        item.setBackGround(Config.error_item)

    elif analysisIndex==4: # PAO and AO from VPS
        if not ("PAO_fromVPS" in globals()) or not ("AO_fromVPS" in globals()):
            return
        if matrixIndex>=0 and matrixIndex<=PAO_fromVPS.PAO_Lmax:
            # Orthonormalization matrix (PAO from VPS)
            l=matrixIndex-PAO_fromVPS.PAO_Lmax-1
            win.matrix.setColumnCount(PAO_fromVPS.PAO_Mul)
            win.matrix.setRowCount(PAO_fromVPS.PAO_Mul)
            for i in range(0, PAO_fromVPS.PAO_Mul):
                for j in range(0, PAO_fromVPS.PAO_Mul):
                    item=QtGui.QTableWidgetItem(("{0:7.4f}").format(PAO_fromVPS.OrthNorm[l][i][j]))
                    win.matrix.setItem(j, i, item)
                    error=PAO_fromVPS.OrthNorm[l][i][j]
                    if i==j:
                        error-=1
                    if abs(error)>Config.OrthNorm_threshold:
                        item.setBackground(Config.error_item)
        elif matrixIndex>=PAO_fromVPS.PAO_Lmax+1 and matrixIndex<=PAO_fromVPS.PAO_Lmax+AO_fromVPS.maxL_pao+1:
            # Orthonormalization matrix (AO from VPS)
            l=matrixIndex-PAO_fromVPS.PAO_Lmax-1
            win.matrix.setColumnCount(AO_fromVPS.max_N)
            win.matrix.setRowCount(AO_fromVPS.max_N)
            for i in range(0, AO_fromVPS.max_N):
                for j in range(0, AO_fromVPS.max_N):
                    item=QtGui.QTableWidgetItem(("{0:7.4f}").format(AO_fromVPS.OrthNorm[l][i][j]))
                    win.matrix.setItem(j, i, item)
                    error=AO_fromVPS.OrthNorm[l][i][j]
                    if i==j:
                        error-=1
                    if abs(error)>Config.OrthNorm_threshold:
                        item.setBackground(Config.error_item)

    elif analysisIndex==5: # PAO and AO after optimization
        win.matrix.setColumnCount(0)
        win.matrix.setRowCount(0)
        
                                                                    
                

    
def performCalculation_PAO(win):
    # QPuthButton calcButton_PAO is clicked

    Z=win.atomNumber.currentIndex()
    paoIndex=win.paoFile.currentIndex()
    vpsIndex=win.vpsFile.currentIndex()
    paoName=Config.paoArr[Z][paoIndex]
    vpsName=Config.vpsArr[Z][vpsIndex]

    pwd=os.getcwd()
    
    # 1: PAO before optimization
    print("Calculation 1: PAO before optimization")
    paoPath_after=os.path.join(Config.dirPath_PAO, paoName)
    inPaoName=re.sub(r"\.pao$", ".inp", paoName)
    inPaoPath=os.path.join(Config.dirPath_inPAO4PAO, inPaoName)
    outPaoPath=os.path.join(Config.dirPath_PAOfromPAO, paoName)
    sysName=""
    with open(paoPath_after, "r") as f1:
        with open(inPaoPath, "w") as f2:
            line=""
            while len(re.findall(r"^\s*Input file\s*$", line))==0:
                line=f1.readline()
            while line[0] != "#":
                line=f1.readline()
            while len(re.findall(r"^\*{32}", line))==0 and len(line)>0:
                re_result=re.findall(r"System\.Name\s*(\S*).*$", line)
                if len(re_result)>0:
                    sysName=re_result[0]
                    print(("System name is {0:s}").format(sysName))
                line=re.sub(r"ocupied","occupied",line)
                f2.write(line)
                line=f1.readline()

    if len(sysName)==0:
        print("Error: System.Name not found")
        return

    os.chdir(Config.dirPath_inPAO4PAO)
    adpack_proc=subprocess.Popen([Config.adpack, inPaoName, "-nt", "1"], stdout=subprocess.PIPE, text=True)
    while True:
        line=adpack_proc.stdout.readline().rstrip()
        if line:
            print(line)
        if (not line) and adpack_proc.poll() is not None:
            break
    shutil.move(sysName+".pao", outPaoPath)

    os.chdir(pwd)
    changeAnalysis(win)

def performCalculation_VPS(win):
    # QPuthButton calcButton_VPS is clicked

    Z=win.atomNumber.currentIndex()
    paoIndex=win.paoFile.currentIndex()
    vpsIndex=win.vpsFile.currentIndex()
    paoName=Config.paoArr[Z][paoIndex]
    vpsName=Config.vpsArr[Z][vpsIndex]

    pwd=os.getcwd()
    
    # 1: PAO from VPS
    print("Calculation 1: PAO from VPS")
    vpsPath=os.path.join(Config.dirPath_VPS, vpsName)
    inVpsName=re.sub(r"\.vps$", ".inp", vpsName)
    inVpsPath=os.path.join(Config.dirPath_inVPS4PAO, inVpsName)
    outPaoName=re.sub(r"\.vps$", ".pao", vpsName)
    outPaoPath=os.path.join(Config.dirPath_PAOfromVPS, outPaoName)
    sysName=""
    with open(vpsPath, "r") as f1:
        with open(inVpsPath, "w") as f2:
            line=""
            while len(re.findall(r"^\s*Input file\s*$", line))==0:
                line=f1.readline()
                if len(line)==0:
                    print("Error: cannot find input parameters")
                    return
            line=f1.readline()
            line=f1.readline()
            # while line[0] != "#":
            #     line=f1.readline()
            while len(re.findall(r"^\*{32}", line))==0 and len(line)>0:
                re_result=re.findall(r"System\.Name\s*(\S*).*$", line)
                if len(re_result)>0:
                    sysName=re_result[0]
                    print(("System name is {0:s}").format(sysName))
                    
                re_result=re.findall(r"(calc\.type\s*)\S*(.*)$", line)
                if len(re_result)>0:
                    line=re_result[0][0]+"pao"+re_result[0][1]+"\n"
                    
                re_result=re.findall(r"(search\.UpperE\s*)(\S*)(.*)$", line)
                if len(re_result)>0:
                    line=re_result[0][0]+("{0:.4f}").format(float(re_result[0][1])*Config.search_UpperE_coeff)+re_result[0][2]+"\n"
                    
                re_result=re.findall(r"(num\.of\.partition\s*)(\S*)(.*)$", line)
                if len(re_result)>0:
                    line=re_result[0][0]+("{0:d}").format(round(int(re_result[0][1])*Config.num_of_partition_coeff))+re_result[0][2]+"\n"
                    
                line=re.sub(r"ocupied","occupied",line)
                f2.write(line)
                line=f1.readline()

    if len(sysName)==0:
        print("Error: System.Name not found")
        return

    os.chdir(Config.dirPath_inVPS4PAO)
    adpack_proc=subprocess.Popen([Config.adpack, inVpsName, "-nt", "1"], stdout=subprocess.PIPE, text=True)
    while True:
    # while False:
        line=adpack_proc.stdout.readline().rstrip()
        if line:
            print(line)
        if (not line) and adpack_proc.poll() is not None:
            break
    shutil.move(sysName+".pao", outPaoPath)

    # 2: AO from VPS
    print("Calculation 2: AO from VPS")
    vpsPath=os.path.join(Config.dirPath_VPS, vpsName)
    inVpsName=re.sub(r"\.vps$", ".inp", vpsName)
    inVpsPath=os.path.join(Config.dirPath_inVPS4AO, inVpsName)
    outAoName=re.sub(r"\.vps$", ".ao", vpsName)
    outAoPath=os.path.join(Config.dirPath_AOfromVPS, outAoName)
    sysName=""
    num_vps=-1
    valence_min=[]
    maxL_pao=-1
    num_pao=-1
    with open(vpsPath, "r") as f1:
        with open(inVpsPath, "w") as f2:
            line=""
            while len(re.findall(r"^\s*Input file\s*$", line))==0:
                line=f1.readline()
                if len(line)==0:
                    print("Error: cannot find input parameters")
                    return
            line=f1.readline()
            line=f1.readline()
            # while line[0] != "#":
            #    line=f1.readline()
            while len(re.findall(r"^\*{32}", line))==0 and len(line)>0:
                re_result=re.findall(r"System\.Name\s*(\S*).*$", line)
                if len(re_result)>0:
                    sysName=re_result[0]
                    print(("System name is {0:s}").format(sysName))
                    
                re_result=re.findall(r"(calc\.type\s*)\S*(.*)$", line)
                if len(re_result)>0:
                    line=re_result[0][0]+"all"+re_result[0][1]+"\n"
                    
                re_result=re.findall(r"(search\.UpperE\s*)(\S*)(.*)$", line)
                if len(re_result)>0:
                    line=re_result[0][0]+("{0:.4f}").format(float(re_result[0][1])*Config.search_UpperE_coeff)+re_result[0][2]+"\n"
                    
                re_result=re.findall(r"(num\.of\.partition\s*)(\S*)(.*)$", line)
                if len(re_result)>0:
                    line=re_result[0][0]+("{0:d}").format(round(int(re_result[0][1])*Config.num_of_partition_coeff))+re_result[0][2]+"\n"

                re_result=re.findall(r"number\.vps\s*(\S*).*$", line)
                if len(re_result)>0:
                    num_vps=int(re_result[0])
                    
                re_result=re.findall(r"^<pseudo\.NandL", line)
                if len(re_result)>0:
                    if num_vps<=0:
                        print("Error: number.vps not found")
                        return
                    f2.write(line)
                    maxL=-1
                    pseudo_NL=[]
                    for i in range(0, num_vps):
                        line=f1.readline()
                        re_result=re.findall(r"^\s*\S*\s*(\S*)\s*(\S*)", line)
                        if len(re_result)>0:
                            pseudo_NL.append([int(re_result[0][0]), int(re_result[0][1])])
                            if maxL<int(re_result[0][1]):
                                maxL=int(re_result[0][1])
                        f2.write(line)
                    # print(pseudo_NL)
                    line=f1.readline()
                    re_result=re.findall(r"^pseudo\.NandL>", line)
                    if len(re_result)==0:
                        print("Error in pseudo.NandL")
                        return
                    for l in range(0, maxL+1):
                        valence_min.append(-1)
                    for data in pseudo_NL:
                        if valence_min[data[1]]<0 or valence_min[data[1]]>data[0]:
                            valence_min[data[1]]=data[0]
                    # print(valence_min)

                re_result=re.findall(r"^maxL\.pao\s*(\S*)", line)
                if len(re_result)>0:
                    maxL_pao=int(re_result[0])
                re_result=re.findall(r"^num\.pao\s*(\S*)", line)
                if len(re_result)>0:
                    num_pao=int(re_result[0])
                    # valence_min, maxL_pao, and num_pao is loaded
                    maxN=-1
                    for l in range(0, maxL_pao+1):
                        valence_min_n=(valence_min[l] if len(valence_min)>l else l+1)+num_pao-1
                        if maxN<valence_min_n:
                            maxN=valence_min_n
                    f2.write(("max.N {0:d} # added \n").format(maxN))
                                
                    
                    
                line=re.sub(r"ocupied","occupied",line)
                f2.write(line)
                line=f1.readline()

    if len(sysName)==0:
        print("Error: System.Name not found")
        return

    os.chdir(Config.dirPath_inVPS4AO)
    adpack_proc=subprocess.Popen([Config.adpack, inVpsName, "-nt", "1"], stdout=subprocess.PIPE, text=True)
    while True:
        line=adpack_proc.stdout.readline().rstrip()
        if line:
            print(line)
        if (not line) and adpack_proc.poll() is not None:
            break
    shutil.move(sysName+".ao", outAoPath)

    os.chdir(pwd)
    changeAnalysis(win)

def outputToHdf5(win):
    # output PAO and AO after optimization
    # check whether "PAO and AO after optimization" is selected
    global PAO_after
    global PAO_reproduced
    
    analysisIndex=win.analysisType.currentIndex()
    if analysisIndex!=5:
        print("Error: 'PAO and AO after optimization' must be selected")
        return

    Z=win.atomNumber.currentIndex()
    paoIndex=win.paoFile.currentIndex()
    vpsIndex=win.vpsFile.currentIndex()
    paoName=Config.paoArr[Z][paoIndex]
    vpsName=Config.vpsArr[Z][vpsIndex]

    paoPrefix=re.sub(r"\.pao$", "", paoName)
    vpsPrefix=re.sub(r"\.vps$", "", vpsName)

    groupName=paoPrefix+"&"+vpsPrefix
    
    # PAO_after.PAOs contains PAO after opt
    # PAO_reproduced contains AO after opt
    
    # data in the HDF5 file
    # PAO and AO are P(r) (= rR(r)) format
    # (groupName) /s0 = [PAO, AO]
    #             /s1 = [PAO, AO]
    #             ...
    #             /r (attribute)
    #             /log_r (attribute)
    #             /datetime (attribute)
    #             /length (attribute) <- len(r) (= len(log_r) = len(PAO) = len(AO))

    with h5py.File(Config.hdf5File, "a") as f:
        # print(list(f.keys()))
        if groupName in list(f.keys()):
            print(("Group {0:s} exists, overwrite it").format(groupName))
            del f[groupName]

        f.create_group(groupName)
        lSize=PAO_after.PAOs.shape[0]
        mulSize=PAO_after.PAOs.shape[1]
        rSize=PAO_after.r.shape[0]
        PAO_and_AO=np.zeros((2,rSize))
        for l in range(0, lSize):
            for mul in range(0, mulSize):
                datasetName=("{0:s}{1:d}").format(Config.azimuthal[l], mul)
                for i in range(0, rSize):
                    PAO_and_AO[0][i]=PAO_after.PAOs[l][mul][i]*PAO_after.r[i]
                    PAO_and_AO[1][i]=PAO_reproduced[l][mul][i]*PAO_after.r[i]
                f[groupName].create_dataset(datasetName, data=PAO_and_AO)
        f[groupName].attrs.create("length", rSize)
        f[groupName].attrs.create("r", PAO_after.r)
        f[groupName].attrs.create("log_r", PAO_after.x)
        f[groupName].attrs.create("datetime", datetime.now().isoformat(" "))
        f[groupName].attrs.create("Z", Z)

        print(("Output to the group {0:s} finished.").format(groupName))


def createDatabase(win):
    # calculate AO after opt and save data to the HDF5 file

    for db in Config.database_list:
        print(db)
        win.analysisType.setCurrentIndex(0)
        Z=db[0]
        paoName=db[1]+".pao"
        vpsName=db[2]+".vps"
        win.atomNumber.setCurrentIndex(Z)
        paoCount=win.paoFile.count()
        vpsCount=win.vpsFile.count()
        pao_found=False
        vps_found=False
        for p in range(0, paoCount):
            if win.paoFile.itemText(p)==paoName:
                win.paoFile.setCurrentIndex(p)
                pao_found=True
                break
        if pao_found==False:
            print(("Error: {0:s} not found").format(paoName))
            return

        for v in range(0, vpsCount):
            if win.vpsFile.itemText(v)==vpsName:
                win.vpsFile.setCurrentIndex(v)
                vps_found=True
                break
        if vps_found==False:
            print(("Error: {0:s} not found").format(vpsName))
            return

        performCalculation_PAO(win)
        performCalculation_VPS(win)

        win.analysisType.setCurrentIndex(5)
        outputToHdf5(win)


        





        
            

