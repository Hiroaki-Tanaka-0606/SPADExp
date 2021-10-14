# OpenMX band
# Viewer for band dispersion


# PyQt5 and PyqQtGraph are necessary
from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg

import Config
import Events

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

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setWindowTitle("OpenMX band dispersion")

        
        font=QtGui.QFont()
        font.setFamilies(Config.fontFamilies)
        font.setPixelSize(Config.fontSize_normal)
        
        vbox=QtGui.QVBoxLayout()
        vbox.setContentsMargins(*Config.ContentsMargins)
        vbox.setAlignment(QtCore.Qt.AlignTop)


        mainWidget=QtGui.QWidget()
        mainWidget.setLayout(vbox)
        self.setCentralWidget(mainWidget)

        # Row 1: Open file button
        row1=QtGui.QHBoxLayout()
        row1.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row1)
        self.openFileButton=QtGui.QPushButton("Open file")
        row1.addWidget(self.openFileButton)

        self.filePath=QtGui.QLabel()
        row1.addWidget(self.filePath)

        # Row 2: configuration
        row2=QtGui.QHBoxLayout()
        row2.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row2)

        ## Energy min and max
        label2A=QtGui.QLabel("E-EF (min and max, eV)")
        row2.addWidget(label2A)
        self.EMin=QtGui.QLineEdit()
        self.EMin.setText("-6.0")
        self.EMax=QtGui.QLineEdit()
        self.EMax.setText("+2.0")
        row2.addWidget(self.EMin)
        row2.addWidget(self.EMax)

        ## Gauss width 
        label2B=QtGui.QLabel("dE (eV)")
        row2.addWidget(label2B)
        self.dE=QtGui.QLineEdit()
        self.dE.setText("0.1")
        row2.addWidget(self.dE)

        ## pixel size (along the energy)
        label2C=QtGui.QLabel("Pixel (eV)")
        row2.addWidget(label2C)
        self.EPixel=QtGui.QLineEdit()
        self.EPixel.setText("0.05")
        row2.addWidget(self.EPixel)

        ## plot button
        self.plotButton=QtGui.QPushButton("Plot")
        row2.addWidget(self.plotButton)


        # Row 3: graph
        row3=QtGui.QVBoxLayout()
        row3.setAlignment(QtCore.Qt.AlignTop)
        vbox.addLayout(row3)
        
        self.plot=pg.PlotWidget()
        self.img=pg.ImageItem()
        row3.addWidget(self.plot)
        self.plot.addItem(self.img)
        
        labelStyle={"font-size":str(Config.fontSize_normal)+"px", "color": "white"}
        self.plot.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.plot.getAxis("bottom").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plot.getAxis("left").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plot.getAxis("bottom").setPen((255,255,255))
        self.plot.getAxis("left").setPen((255,255,255))
        self.plot.getAxis("bottom").setTextPen((255,255,255))
        self.plot.getAxis("left").setTextPen((255,255,255))
        self.plot.getAxis("bottom").setLabel(**labelStyle)
        self.plot.getAxis("left").setLabel(**labelStyle)
        self.plot.setLabel(axis="left", text="E-EF (eV)")

        self.plot.showGrid(x=True, y=True, alpha=1.0)
        self.plot.getAxis("bottom").setZValue(1)
        self.plot.getAxis("left").setZValue(1)



app=QtGui.QApplication([])
# win: MainWindow object
win=MainWindow()
font=QtGui.QFont()
font.setPixelSize(Config.fontSize_normal)
font.setFamilies(Config.fontFamilies)
win.setFont(font)

win.openFileButton.clicked.connect(lambda: Events.openFile(win))
win.plotButton.clicked.connect(lambda: Events.plot(win))

pg.setConfigOptions(antialias=True)
win.show()
app.exec_()

        
