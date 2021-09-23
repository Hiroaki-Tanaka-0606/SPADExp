# OpenMX orbitals
# Analysis tool for pao (pseudo atomic orbitals) and ao (all-electron atomic orbitals)


# PyQt5 and PyqQtGraph are necessary
from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import os
import re

import Config
import Events

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setWindowTitle("OpenMX orbitals")
        
        vbox=QtGui.QVBoxLayout()
        vbox.setContentsMargins(*Config.ContentsMargins)
        vbox.setAlignment(QtCore.Qt.AlignTop)


        mainWidget=QtGui.QWidget()
        mainWidget.setLayout(vbox)
        self.setCentralWidget(mainWidget)

        # Row 1: working directory
        row1=QtGui.QHBoxLayout()
        row1.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row1)

        dirPath=QtGui.QLabel()
        row1.addWidget(dirPath)
        dirPath.setText(("Directory: {0:s}").format(Config.workingDirectory))

        # Row 2: atomic number, PAO and VPS files
        row2=QtGui.QHBoxLayout()
        row2.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row2)

        ## atomic number
        zLabel=QtGui.QLabel("Atomic number")
        row2.addWidget(zLabel)
        self.atomNumber=QtGui.QComboBox()
        row2.addWidget(self.atomNumber)
        for z, el in enumerate(Config.el_symbol):
            self.atomNumber.addItem(("{0:d}: {1:s}").format(z, el))

        ## PAO file
        pLabel=QtGui.QLabel("PAO file")
        row2.addWidget(pLabel)
        self.paoFile=QtGui.QComboBox()
        row2.addWidget(self.paoFile)

        ## VPS file
        vLabel=QtGui.QLabel("VPS file")
        row2.addWidget(vLabel)
        self.vpsFile=QtGui.QComboBox()
        row2.addWidget(self.vpsFile)

        # Row 3: analysis type & calculation button
        row3=QtGui.QHBoxLayout()
        row3.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row3)

        mLabel=QtGui.QLabel("Analysis type")
        row3.addWidget(mLabel)

        self.analysisType=QtGui.QComboBox()
        for analysisType in Config.analysisTypes:
            self.analysisType.addItem(analysisType)
        row3.addWidget(self.analysisType)

        self.calcButton=QtGui.QPushButton("Calculate")
        row3.addWidget(self.calcButton)

        # Row 4: select orbitals
        row4=QtGui.QHBoxLayout()
        row4.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row4)
        
        self.orbitalTable=QtGui.QTableWidget()
        row4.addWidget(self.orbitalTable)
        self.orbitalTable.horizontalHeader().setSectionResizeMode(QtGui.QHeaderView.ResizeToContents)
        self.orbitalTable.verticalHeader().setSectionResizeMode(QtGui.QHeaderView.Stretch)
        self.orbitalTable.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)

        # Row 5: graph prefs
        row5=QtGui.QHBoxLayout()
        row5.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row5)

        fLabel=QtGui.QLabel("Radial function type")
        row5.addWidget(fLabel)

        self.radialType=QtGui.QButtonGroup()
        radial_r=QtGui.QRadioButton("R(r)")
        self.radialType.addButton(radial_r)
        row5.addWidget(radial_r)
        radial_p=QtGui.QRadioButton("P(r)=r R(r)")
        radial_p.setChecked(True)
        self.radialType.addButton(radial_p)
        row5.addWidget(radial_p)

        # Row 6: graph
        row6=QtGui.QHBoxLayout()
        self.orbitalGraph=pg.PlotWidget()
        self.orbitalGraph.addLegend()
        row6.addWidget(self.orbitalGraph)
        vbox.addLayout(row6)


app=QtGui.QApplication([])
# win: MainWindow object
win=MainWindow()
font=QtGui.QFont()
font.setPixelSize(Config.fontSize_normal)
font.setFamilies(Config.fontFamilies)
win.setFont(font)


# search for "VPS" and "PAO" directories in Config.workingDirectory
print(("Search for PAO and VPS dirs in {0:s}").format(Config.workingDirectory))
Config.dirPath_PAO=os.path.join(Config.workingDirectory, "PAO")
Config.dirPath_VPS=os.path.join(Config.workingDirectory, "VPS")
Config.dirPath_inPAO4PAO=os.path.join(Config.workingDirectory, "PAO_input_for_PAO")
Config.dirPath_PAOfromPAO=os.path.join(Config.workingDirectory, "PAO_from_PAO")

for dirPath in [Config.dirPath_PAO, Config.dirPath_VPS]:
    if os.path.isdir(dirPath):
        print(("OK: {0:s} is a directory").format(dirPath))
    else:
        print(("Not OK: {0:s} is not a directory or does not exist").format(dirPath))

for dirPath in [Config.dirPath_inPAO4PAO, Config.dirPath_PAOfromPAO]:
    if os.path.isdir(dirPath):
        print(("Directory {0:s} exists").format(dirPath))
    else:
        print(("Directory {0:s} does not exist, make it").format(dirPath))
        os.mkdir(dirPath)
        
Config.paoArr=[]
Config.vpsArr=[]
for el in Config.el_symbol:
    Config.paoArr.append([])
    Config.vpsArr.append([])
for paoName in os.listdir(Config.dirPath_PAO):
    el_index=-1
    el_name=""
    for i, el in enumerate(Config.el_symbol):
        el_match=re.findall("^("+el+r").*\.pao$", paoName)
        if len(el_match)>0 and len(el_match[0])>len(el_name):
            el_name=el_match[0]
            el_index=i
    if len(el_name)>0:
        Config.paoArr[el_index].append(paoName)
    else:
        print(("Error: file {0:s} does not match to any element").format(paoName))

for vpsName in os.listdir(Config.dirPath_VPS):
    el_index=-1
    el_name=""
    for i, el in enumerate(Config.el_symbol):
        el_match=re.findall("^("+el+r").*\.vps$", vpsName)
        if len(el_match)>0 and len(el_match[0])>len(el_name):
            el_name=el_match[0]
            el_index=i
    if len(el_name)>0:
        Config.vpsArr[el_index].append(vpsName)
    else:
        print(("Error: file {0:s} does not match to any element").format(vpsName))

# E (Z=0) uses Kr (Z=36)
for Kr_pao in Config.paoArr[36]:
    Config.paoArr[0].append(Kr_pao)
        
for paoList in Config.paoArr:
    paoList.sort()

# objects
PAO_after=None
PAO_before=None
Selected_orbitals=None

# Event
win.atomNumber.currentIndexChanged.connect(lambda: Events.changeAtom(win))
Events.changeAtom(win)

win.paoFile.currentIndexChanged.connect(lambda: Events.changeAnalysis(win))
win.vpsFile.currentIndexChanged.connect(lambda: Events.changeAnalysis(win))
win.analysisType.currentIndexChanged.connect(lambda: Events.changeAnalysis(win))
win.calcButton.clicked.connect(lambda: Events.performCalculation(win))
win.orbitalTable.cellClicked.connect(lambda row, column: Events.selectOrbital(win, row, column))
for b in win.radialType.buttons():
    b.clicked.connect(lambda: Events.drawOrbitalGraph(win))

pg.setConfigOptions(antialias=True)
win.show()
app.exec_()

        
