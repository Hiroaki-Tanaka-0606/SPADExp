# OpenMX band
# Viewer for band dispersion


# PyQt5 and PyqQtGraph are necessary
from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import pyqtgraph.opengl as gl


import Config
import Events

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setWindowTitle("OpenMX band dispersion")

        
        font=QtGui.QFont()
        font.setFamilies(Config.fontFamilies)
        font.setPixelSize(Config.fontSize_normal)
        
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

        # Row 2: configuration
        row2=QtWidgets.QHBoxLayout()
        row2.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row2)

        ## Energy min and max
        label2A=QtWidgets.QLabel("E-EF (min and max, eV)")
        row2.addWidget(label2A)
        self.EMin=QtWidgets.QLineEdit()
        self.EMin.setText("-6.0")
        self.EMax=QtWidgets.QLineEdit()
        self.EMax.setText("+2.0")
        row2.addWidget(self.EMin)
        row2.addWidget(self.EMax)

        ## Gauss width 
        label2B=QtWidgets.QLabel("dE (eV)")
        row2.addWidget(label2B)
        self.dE=QtWidgets.QLineEdit()
        self.dE.setText("0.1")
        row2.addWidget(self.dE)

        ## pixel size (along the energy)
        label2C=QtWidgets.QLabel("Pixel (eV)")
        row2.addWidget(label2C)
        self.EPixel=QtWidgets.QLineEdit()
        self.EPixel.setText("0.05")
        row2.addWidget(self.EPixel)

        ## plot button
        self.plotButton=QtWidgets.QPushButton("Plot")
        row2.addWidget(self.plotButton)


        # Row 3 left: graph (2D / 3D)
        row3=QtWidgets.QHBoxLayout()
        row3.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row3)

        self.graphTab=QtWidgets.QTabWidget()
        row3.addWidget(self.graphTab, 3)
        
        ## 2D
        self.plot=pg.PlotWidget()
        self.img=pg.ImageItem()
        # row3.addWidget(self.plot)
        self.graphTab.addTab(self.plot, "2D")
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

        self.vLine=pg.InfiniteLine(angle=90, movable=False, pen=Config.pen1)
        self.plot.addItem(self.vLine, ignoreBounds=True)
        self.hLine=pg.InfiniteLine(angle=0, movable=False, pen=Config.pen1)
        self.plot.addItem(self.hLine, ignoreBounds=True)

        ## 3D
        ## the four panels are placed like the following
        ## [3D] | [Ey] [xy]
        ##      |      [Ex]

        row3_3d=QtWidgets.QHBoxLayout()
        row3_3d.setAlignment(QtCore.Qt.AlignLeft)

        tabWidget=QtWidgets.QWidget()
        tabWidget.setLayout(row3_3d)
        self.graphTab.addTab(tabWidget, "3D")
        
        self.plot3D=gl.GLViewWidget()
        self.plot3D.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.plot3D.setMinimumSize(Config.plot3D_minWidth, Config.plot3D_minHeight)
        row3_3d.addWidget(self.plot3D)

        self.plot3=pg.GraphicsLayoutWidget()
        self.plot3.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        row3_3d.addWidget(self.plot3)

        self.plotEy=self.plot3.addPlot()
        self.plotxy=self.plot3.addPlot()
        self.plot3.nextRow()
        self.plot3.nextColumn()
        self.plotEx=self.plot3.addPlot()

        self.imgEx=pg.ImageItem()
        self.imgEy=pg.ImageItem()
        self.imgxy=pg.ImageItem()
        self.plotEx.addItem(self.imgEx)
        self.plotEy.addItem(self.imgEy)
        self.plotxy.addItem(self.imgxy)

        self.plotEx.setLabel(axis="left", text="E-EF (eV)")
        self.plotEy.setLabel(axis="bottom", text="E-EF (eV)")
        
        self.plotEx.getAxis("bottom").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotEx.getAxis("left").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotEx.getAxis("bottom").setPen((255,255,255))
        self.plotEx.getAxis("left").setPen((255,255,255))
        self.plotEx.getAxis("bottom").setTextPen((255,255,255))
        self.plotEx.getAxis("left").setTextPen((255,255,255))
        self.plotEx.getAxis("bottom").setLabel(**labelStyle)
        self.plotEx.getAxis("left").setLabel(**labelStyle)
        self.plotEx.setLabel(axis="left", text="E-EF (eV)")
        self.plotEx.showGrid(x=True, y=True, alpha=1.0)
        self.plotEx.getAxis("bottom").setZValue(1)
        self.plotEx.getAxis("left").setZValue(1)

        
        self.plotEy.getAxis("bottom").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotEy.getAxis("left").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotEy.getAxis("bottom").setPen((255,255,255))
        self.plotEy.getAxis("left").setPen((255,255,255))
        self.plotEy.getAxis("bottom").setTextPen((255,255,255))
        self.plotEy.getAxis("left").setTextPen((255,255,255))
        self.plotEy.getAxis("bottom").setLabel(**labelStyle)
        self.plotEy.getAxis("left").setLabel(**labelStyle)
        self.plotEy.setLabel(axis="bottom", text="E-EF (eV)")
        self.plotEy.showGrid(x=True, y=True, alpha=1.0)
        self.plotEy.getAxis("bottom").setZValue(1)
        self.plotEy.getAxis("left").setZValue(1)

        
        self.plotxy.getAxis("bottom").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotxy.getAxis("left").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotxy.getAxis("bottom").setPen((255,255,255))
        self.plotxy.getAxis("left").setPen((255,255,255))
        self.plotxy.getAxis("bottom").setTextPen((255,255,255))
        self.plotxy.getAxis("left").setTextPen((255,255,255))
        self.plotxy.getAxis("bottom").setLabel(**labelStyle)
        self.plotxy.getAxis("left").setLabel(**labelStyle)
        self.plotxy.showGrid(x=True, y=True, alpha=1.0)
        self.plotxy.getAxis("bottom").setZValue(1)
        self.plotxy.getAxis("left").setZValue(1)

        
        self.vLineEx=pg.InfiniteLine(angle=90, movable=False, pen=Config.pen1)
        self.plotEx.addItem(self.vLineEx, ignoreBounds=True)
        self.hLineEx=pg.InfiniteLine(angle=0, movable=False, pen=Config.pen1)
        self.plotEx.addItem(self.hLineEx, ignoreBounds=True)

        self.vLineEy=pg.InfiniteLine(angle=90, movable=False, pen=Config.pen1)
        self.plotEy.addItem(self.vLineEy, ignoreBounds=True)
        self.hLineEy=pg.InfiniteLine(angle=0, movable=False, pen=Config.pen1)
        self.plotEy.addItem(self.hLineEy, ignoreBounds=True)

        self.vLinexy=pg.InfiniteLine(angle=90, movable=False, pen=Config.pen1)
        self.plotxy.addItem(self.vLinexy, ignoreBounds=True)
        self.hLinexy=pg.InfiniteLine(angle=0, movable=False, pen=Config.pen1)
        self.plotxy.addItem(self.hLinexy, ignoreBounds=True)

        self.bLineEx=pg.InfiniteLine(angle=0, movable=False, pen=Config.pen2)
        self.plotEx.addItem(self.bLineEx, ignoreBounds=True)
        self.bLineEy=pg.InfiniteLine(angle=90, movable=False, pen=Config.pen2)
        self.plotEy.addItem(self.bLineEy, ignoreBounds=True)

        def changeKBIndices(e):
            if e.key()==QtCore.Qt.Key_Down:
                self.bIndex.setValue(self.bIndex.value()-1)
            elif e.key()==QtCore.Qt.Key_Up:
                self.bIndex.setValue(self.bIndex.value()+1)
            elif e.key()==QtCore.Qt.Key_Right:
                self.kxIndex.setValue(self.kxIndex.value()+1)
            elif e.key()==QtCore.Qt.Key_Left:
                self.kxIndex.setValue(self.kxIndex.value()-1)

        def changeKXYIndices(e):
            if e.key()==QtCore.Qt.Key_Down:
                self.kyIndex.setValue(self.kyIndex.value()-1)
            elif e.key()==QtCore.Qt.Key_Up:
                self.kyIndex.setValue(self.kyIndex.value()+1)
            elif e.key()==QtCore.Qt.Key_Right:
                self.kxIndex.setValue(self.kxIndex.value()+1)
            elif e.key()==QtCore.Qt.Key_Left:
                self.kxIndex.setValue(self.kxIndex.value()-1)
            elif e.key()==QtCore.Qt.Key_PageUp:
                self.eIndex.setValue(self.eIndex.value()+1)
            elif e.key()==QtCore.Qt.Key_PageDown:
                self.eIndex.setValue(self.eIndex.value()-1)
            elif e.key()==QtCore.Qt.Key_Home:
                self.bIndex.setValue(self.bIndex.value()+1)
            elif e.key()==QtCore.Qt.Key_End:
                self.bIndex.setValue(self.bIndex.value()-1)
            else:
                return
            Events.plot3(win)
                
        self.plot.keyPressEvent=changeKBIndices
        self.plot3.keyPressEvent=changeKXYIndices

        # Row 3 right: properties
        vbox2=QtWidgets.QVBoxLayout()
        vbox2.setAlignment(QtCore.Qt.AlignTop)
        row3.addLayout(vbox2, 1)

        ## row 1: k index
        row3r1=QtWidgets.QHBoxLayout()
        row3r1.setAlignment(QtCore.Qt.AlignLeft)
        vbox2.addLayout(row3r1)
        
        label3r1A=QtWidgets.QLabel("kx index")
        row3r1.addWidget(label3r1A)
        self.kxIndex=QtWidgets.QSpinBox()
        self.kxIndex.setSingleStep(1)
        self.kxIndex.setMinimum(0)
        row3r1.addWidget(self.kxIndex)

        self.kxValue=QtWidgets.QLabel()
        row3r1.addWidget(self.kxValue)

        label3r1B=QtWidgets.QLabel("ky index")
        row3r1.addWidget(label3r1B)
        self.kyIndex=QtWidgets.QSpinBox()
        self.kyIndex.setSingleStep(1)
        self.kyIndex.setMinimum(0)
        row3r1.addWidget(self.kyIndex)

        self.kyValue=QtWidgets.QLabel()
        row3r1.addWidget(self.kyValue)

        ## row 2: band index and Constant energy index
        row3r2=QtWidgets.QHBoxLayout()
        row3r2.setAlignment(QtCore.Qt.AlignLeft)
        vbox2.addLayout(row3r2)
        
        label3r2A=QtWidgets.QLabel("Band index")
        row3r2.addWidget(label3r2A)
        self.bIndex=QtWidgets.QSpinBox()
        self.bIndex.setSingleStep(1)
        self.bIndex.setMinimum(0)
        row3r2.addWidget(self.bIndex)
        self.bValue=QtWidgets.QLabel()
        row3r2.addWidget(self.bValue)

        label3r2B=QtWidgets.QLabel("Energy index")
        row3r2.addWidget(label3r2B)
        self.eIndex=QtWidgets.QSpinBox()
        self.eIndex.setSingleStep(1)
        self.eIndex.setMinimum(0)
        row3r2.addWidget(self.eIndex)
        self.eValue=QtWidgets.QLabel()
        row3r2.addWidget(self.eValue)

        ## row 3: Up or Dn
        row3r3=QtWidgets.QHBoxLayout()
        row3r3.setAlignment(QtCore.Qt.AlignLeft)
        vbox2.addLayout(row3r3)

        self.UpDn=QtWidgets.QButtonGroup()
        self.UpButton=QtWidgets.QRadioButton("Up")
        self.DnButton=QtWidgets.QRadioButton("Dn")
        row3r3.addWidget(self.UpButton)
        row3r3.addWidget(self.DnButton)
        self.UpDn.addButton(self.UpButton)
        self.UpDn.addButton(self.DnButton)

        ## row 4: Atom
        row3r4=QtWidgets.QHBoxLayout()
        row3r4.setAlignment(QtCore.Qt.AlignLeft)
        vbox2.addLayout(row3r4)

        label3r4A=QtWidgets.QLabel("Atom")
        row3r4.addWidget(label3r4A)

        self.Atom=QtWidgets.QComboBox()
        self.Atom.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        row3r4.addWidget(self.Atom)

        ## row 5: LCAO table
        self.LCAOTable=QtWidgets.QTableWidget()
        vbox2.addWidget(self.LCAOTable)
        



app=QtWidgets.QApplication([])
# win: MainWindow object
win=MainWindow()
font=QtGui.QFont()
font.setPixelSize(Config.fontSize_normal)
font.setFamilies(Config.fontFamilies)
win.setFont(font)

win.openFileButton.clicked.connect(lambda: Events.openFile(win))

def plotEvent(win):
    if Events.Dimension==1:
        Events.plot(win)
        win.graphTab.setCurrentIndex(0)
    elif Events.Dimension==2:
        Events.makeDispersion3(win)
        win.graphTab.setCurrentIndex(1)
    else:
        print("Dimension error")
        return

def cursorEvent(win):
    if Events.Dimension==1:
        Events.drawCursor(win)
    elif Events.Dimension==2:
        Events.drawCursor3(win)
    else:
        print("Dimension error")
        return

win.plotButton.clicked.connect(lambda: plotEvent(win))

win.kxIndex.valueChanged.connect(lambda: cursorEvent(win))
win.kyIndex.valueChanged.connect(lambda: cursorEvent(win))
win.bIndex.valueChanged.connect(lambda: cursorEvent(win))
win.eIndex.valueChanged.connect(lambda: cursorEvent(win))
win.UpButton.clicked.connect(lambda: cursorEvent(win))
win.DnButton.clicked.connect(lambda: cursorEvent(win))
win.Atom.currentIndexChanged.connect(lambda: Events.makeLCAOTable(win))

pg.setConfigOptions(antialias=True)
win.show()
app.exec_()

        
