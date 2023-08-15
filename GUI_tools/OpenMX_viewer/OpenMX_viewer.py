# OpenMX viewer
# Viewer for vps (pseudopotential) and pao ((optimized) pseudo atomic orbitals)


# PyQt5 and PyqQtGraph are necessary
from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg

import Config
import Events

blockNames=[]
blockData=[]

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setWindowTitle("OpenMX viewer")
        
        vbox=QtWidgets.QVBoxLayout()
        vbox.setContentsMargins(*Config.ContentsMargins)
        vbox.setAlignment(QtCore.Qt.AlignTop)


        mainWidget=QtWidgets.QWidget()
        mainWidget.setLayout(vbox)
        self.setCentralWidget(mainWidget)

        # Row 1: Open file button
        row1=QtWidgets.QHBoxLayout()
        row1.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row1)
        self.openFileButton=QtWidgets.QPushButton("Open file")
        row1.addWidget(self.openFileButton)

        self.filePath=QtWidgets.QLabel()
        row1.addWidget(self.filePath)

        # Row 2: select data block
        row2=QtWidgets.QHBoxLayout()
        row2.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row2)
        label2A=QtWidgets.QLabel("Data block")
        row2.addWidget(label2A)
        self.selectDataBlock=QtWidgets.QComboBox()
        self.selectDataBlock.addItem("---- The choices are displayed after the file is opened ----")
        row2.addWidget(self.selectDataBlock)

        # Row 3: select x axis and y axes
        row3=QtWidgets.QHBoxLayout()
        row3.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row3)
        # Column 3A: x axis
        self.xCol=QtWidgets.QVBoxLayout()
        self.xCol.setAlignment(QtCore.Qt.AlignTop)
        row3.addLayout(self.xCol)
        label3AA=QtWidgets.QLabel("X axis")
        self.xCol.addWidget(label3AA)
        self.xAxis=QtWidgets.QButtonGroup()
        # Column 3B: y axes
        self.yCol=QtWidgets.QVBoxLayout()
        self.yCol.setAlignment(QtCore.Qt.AlignTop)
        row3.addLayout(self.yCol)
        label3BA=QtWidgets.QLabel("Y axes")
        self.yCol.addWidget(label3BA)
        self.yAxes=QtWidgets.QButtonGroup()
        self.yAxes.setExclusive(False)

        # Row 4: graph
        row4=QtWidgets.QVBoxLayout()
        row4.setAlignment(QtCore.Qt.AlignTop)
        row3.addLayout(row4)

        self.graphPlot=pg.PlotWidget()
        row4.addWidget(self.graphPlot)

        # Row 5: export to hdf5
        row5=QtWidgets.QVBoxLayout()
        row5.setAlignment(QtCore.Qt.AlignTop)
        vbox.addLayout(row5)

        self.exportHDF5=QtWidgets.QPushButton("Export to HDF5")
        row5.addWidget(self.exportHDF5)


app=QtWidgets.QApplication([])
# win: MainWindow object
win=MainWindow()
font=QtGui.QFont()
font.setPixelSize(Config.fontSize_normal)
font.setFamilies(Config.fontFamilies)
win.setFont(font)
win.openFileButton.clicked.connect(lambda: Events.openFile(win))
win.selectDataBlock.currentIndexChanged.connect(lambda: Events.setDataBlock(win))
win.xAxis.buttonToggled.connect(lambda: Events.setGraph(win))
win.yAxes.buttonToggled.connect(lambda: Events.setGraph(win))
win.exportHDF5.clicked.connect(lambda: Events.exportData(win))
pg.setConfigOptions(antialias=True)
win.show()
app.exec_()

        
