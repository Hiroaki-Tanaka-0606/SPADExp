# Events
# clicked, index changed etc.

import Config
import PAO_parser as PAOp
import os
import re
import subprocess
import shutil
from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import numpy as np

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
    win.orbitalGraph.clear()
    
def changeAnalysis(win):
    # QComboBox analysisType is changed

    analysisIndex=win.analysisType.currentIndex()
    Z=win.atomNumber.currentIndex()
    paoIndex=win.paoFile.currentIndex()
    vpsIndex=win.vpsFile.currentIndex()
    paoName=Config.paoArr[Z][paoIndex]
    vpsName=Config.vpsArr[Z][vpsIndex]
    
    if analysisIndex==1: # PAO: Before and after optimization
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
        

        global PAO_after
        global PAO_before
        print(("Read {0:s}").format(paoPath_after))
        PAO_after=PAOp.PAO_parser(paoPath_after, True)
        print(("Read {0:s}").format(paoPath_before))
        PAO_before=PAOp.PAO_parser(paoPath_before, False)
        
        global Selected_orbitals
        Selected_orbitals=np.zeros((PAO_after.PAO_Lmax+1, PAO_after.PAO_Mul))

        # construct orbital table
        win.orbitalTable.setRowCount(PAO_after.PAO_Lmax+1)
        win.orbitalTable.setColumnCount(PAO_after.PAO_Mul)
        for l in range(0, PAO_after.PAO_Lmax+1):
            win.orbitalTable.setVerticalHeaderItem(l, QtGui.QTableWidgetItem(Config.azimuthal[l]))
            
        for n in range(1, PAO_after.PAO_Mul+1):
            for l in range(0, PAO_after.PAO_Lmax+1):
                item=QtGui.QTableWidgetItem(str(n)+Config.azimuthal[l])
                item.setBackground(Config.unselected_item)
                win.orbitalTable.setItem(l, n-1, item)


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
    
    analysisIndex=win.analysisType.currentIndex()
    if analysisIndex==1: # PAO: Before and after optimization
        win.orbitalGraph.clear()
        radialType=win.radialType.checkedButton().text()[0]
        for l in range(0, len(Selected_orbitals)):
            for mul in range(0, len(Selected_orbitals[l])):
                if Selected_orbitals[l][mul]==1:
                    if radialType=="P":
                        win.orbitalGraph.plot(y=PAO_after.PAOs[l][mul]*PAO_after.r, x=PAO_after.r, pen=Config.pen_after, name=str(mul+1)+Config.azimuthal[l]+" after")
                        win.orbitalGraph.plot(y=PAO_before.PAOs[l][mul]*PAO_before.r, x=PAO_before.r, pen=Config.pen_before, name=str(mul+1)+Config.azimuthal[l]+" before")
                    else:
                        win.orbitalGraph.plot(y=PAO_after.PAOs[l][mul], x=PAO_after.r, pen=Config.pen_after, name=str(mul+1)+Config.azimuthal[l]+" after")
                        win.orbitalGraph.plot(y=PAO_before.PAOs[l][mul], x=PAO_before.r, pen=Config.pen_before, name=str(mul+1)+Config.azimuthal[l]+" before")
        
def performCalculation(win):
    # QPuthButton calcButton is clicked

    Z=win.atomNumber.currentIndex()
    paoIndex=win.paoFile.currentIndex()
    vpsIndex=win.vpsFile.currentIndex()
    paoName=Config.paoArr[Z][paoIndex]
    vpsName=Config.vpsArr[Z][vpsIndex]

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
                re_result=re.findall(r"^System.Name\s*(\S*).*$", line)
                if len(re_result)>0:
                    sysName=re_result[0]
                    print(("System name is {0:s}").format(sysName))
                line=re.sub(r"ocupied","occupied",line)
                f2.write(line)
                line=f1.readline()

    if len(sysName)==0:
        print("Error: System.Name not found")
        return

    pwd=os.getcwd()
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
