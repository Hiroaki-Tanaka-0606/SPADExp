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
        
        vbox=QtWidgets.QVBoxLayout()
        vbox.setContentsMargins(*Config.ContentsMargins)
        vbox.setAlignment(QtCore.Qt.AlignTop)


        mainWidget=QtWidgets.QWidget()
        mainWidget.setLayout(vbox)
        self.setCentralWidget(mainWidget)

        # Row 1: working directory
        row1=QtWidgets.QHBoxLayout()
        row1.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row1)

        dirPath=QtWidgets.QLabel()
        row1.addWidget(dirPath)
        dirPath.setText(("Directory: {0:s}").format(Config.workingDirectory))

        # Row 2: atomic number, PAO and VPS files
        row2=QtWidgets.QHBoxLayout()
        row2.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row2)

        ## atomic number
        zLabel=QtWidgets.QLabel("Atomic number")
        row2.addWidget(zLabel)
        self.atomNumber=QtWidgets.QComboBox()
        self.atomNumber.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        row2.addWidget(self.atomNumber)
        for z, el in enumerate(Config.el_symbol):
            self.atomNumber.addItem(("{0:d}: {1:s}").format(z, el))

        ## PAO file
        pLabel=QtWidgets.QLabel("PAO file")
        row2.addWidget(pLabel)
        self.paoFile=QtWidgets.QComboBox()
        self.paoFile.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        row2.addWidget(self.paoFile)

        ## VPS file
        vLabel=QtWidgets.QLabel("VPS file")
        row2.addWidget(vLabel)
        self.vpsFile=QtWidgets.QComboBox()
        self.vpsFile.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        row2.addWidget(self.vpsFile)

        # Row 3: analysis type & calculation button
        row3=QtWidgets.QHBoxLayout()
        row3.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row3)

        mLabel=QtWidgets.QLabel("Analysis type")
        row3.addWidget(mLabel)

        self.analysisType=QtWidgets.QComboBox()
        for analysisType in Config.analysisTypes:
            self.analysisType.addItem(analysisType)
        row3.addWidget(self.analysisType)

        self.calcButton_PAO=QtWidgets.QPushButton("Calculate for PAO")
        row3.addWidget(self.calcButton_PAO)
        
        self.calcButton_VPS=QtWidgets.QPushButton("Calculate for VPS")
        row3.addWidget(self.calcButton_VPS)

        self.outToHdf5Button=QtWidgets.QPushButton("Output to the hdf5 file")
        row3.addWidget(self.outToHdf5Button)

        self.databaseButton=QtWidgets.QPushButton("Create database")
        row3.addWidget(self.databaseButton)


        # Row 4: select orbitals
        row4=QtWidgets.QHBoxLayout()
        row4.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row4)
        
        self.orbitalTable=QtWidgets.QTableWidget()
        row4.addWidget(self.orbitalTable)
        self.orbitalTable.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        self.orbitalTable.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.orbitalTable.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)

        # Row 5: graph and matrix prefs
        row5=QtWidgets.QHBoxLayout()
        row5.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row5)

        fLabel=QtWidgets.QLabel("Radial function type")
        row5.addWidget(fLabel)

        self.radialType=QtWidgets.QButtonGroup()
        radial_r=QtWidgets.QRadioButton("R(r)")
        self.radialType.addButton(radial_r)
        row5.addWidget(radial_r)
        radial_p=QtWidgets.QRadioButton("P(r)=r R(r)")
        radial_p.setChecked(True)
        self.radialType.addButton(radial_p)
        row5.addWidget(radial_p)

        mLabel=QtWidgets.QLabel("Matrix")
        row5.addWidget(mLabel)
        self.matrixType=QtWidgets.QComboBox()
        self.matrixType.addItem("The choices are presented after an analysis is selected")
        row5.addWidget(self.matrixType)

        # Row 6: graph and matrix
        row6=QtWidgets.QHBoxLayout()
        self.orbitalGraph=pg.PlotWidget()
        self.orbitalGraph.addLegend()
        row6.addWidget(self.orbitalGraph)
        vbox.addLayout(row6)

        self.matrix=QtWidgets.QTableWidget()
        row6.addWidget(self.matrix)
        self.matrix.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        self.matrix.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        self.matrix.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)


app=QtWidgets.QApplication([])
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
Config.dirPath_inVPS4PAO=os.path.join(Config.workingDirectory, "VPS_input_for_PAO")
Config.dirPath_PAOfromVPS=os.path.join(Config.workingDirectory, "PAO_from_VPS")
Config.dirPath_inVPS4AO=os.path.join(Config.workingDirectory, "VPS_input_for_AO")
Config.dirPath_AOfromVPS=os.path.join(Config.workingDirectory, "AO_from_VPS")

for dirPath in [Config.dirPath_PAO, Config.dirPath_VPS]:
    if os.path.isdir(dirPath):
        print(("OK: {0:s} is a directory").format(dirPath))
    else:
        print(("Not OK: {0:s} is not a directory or does not exist").format(dirPath))

for dirPath in [Config.dirPath_inPAO4PAO, Config.dirPath_PAOfromPAO, \
                Config.dirPath_inVPS4PAO, Config.dirPath_PAOfromVPS,\
                Config.dirPath_inVPS4AO,  Config.dirPath_AOfromVPS,]:
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


# Event
win.atomNumber.currentIndexChanged.connect(lambda: Events.changeAtom(win))
Events.changeAtom(win)

win.paoFile.currentIndexChanged.connect(lambda: Events.changeAnalysis(win))
win.vpsFile.currentIndexChanged.connect(lambda: Events.changeAnalysis(win))
win.analysisType.currentIndexChanged.connect(lambda: Events.changeAnalysis(win))
win.calcButton_PAO.clicked.connect(lambda: Events.performCalculation_PAO(win))
win.calcButton_VPS.clicked.connect(lambda: Events.performCalculation_VPS(win))
win.outToHdf5Button.clicked.connect(lambda: Events.outputToHdf5(win))
win.databaseButton.clicked.connect(lambda: Events.createDatabase(win))
win.orbitalTable.cellClicked.connect(lambda row, column: Events.selectOrbital(win, row, column))
for b in win.radialType.buttons():
    b.clicked.connect(lambda: Events.drawOrbitalGraph(win))
win.matrixType.currentIndexChanged.connect(lambda: Events.changeMatrix(win))

pg.setConfigOptions(antialias=True)
win.show()
app.exec_()

        
