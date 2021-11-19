# calcPSF GUI

# PyQt5 and PyqQtGraph are necessary
from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import pyqtgraph.opengl as gl


import Config
from lib import objs, Events

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setWindowTitle("calcPSF GUI")

        
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

        # Row 2: configuration (energy)
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

        # Row 3: configuration (final state)
        row3=QtGui.QHBoxLayout()
        row3.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row3)

        self.finalState=QtGui.QButtonGroup()
        self.PWButton=QtGui.QRadioButton("Plane wave")
        self.CalcButton=QtGui.QRadioButton("Calculated")
        row3.addWidget(self.PWButton)
        row3.addWidget(self.CalcButton)
        self.finalState.addButton(self.PWButton)
        self.finalState.addButton(self.CalcButton)
        self.PWButton.setChecked(True)

        label3A=QtGui.QLabel("Calc. step of the final states (a.u.^-1)")
        row3.addWidget(label3A)
        self.finalStates_step=QtGui.QLineEdit()
        self.finalStates_step.setText("0.01")
        row3.addWidget(self.finalStates_step)

        # Row 4: configuration (polarization)
        row4=QtGui.QHBoxLayout()
        row4.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row4)

        self.polarization=QtGui.QButtonGroup()
        self.linear=QtGui.QRadioButton("Linear")
        self.rCircular=QtGui.QRadioButton("Right circular")
        self.lCircular=QtGui.QRadioButton("Left circular")
        row4.addWidget(self.linear)
        row4.addWidget(self.rCircular)
        row4.addWidget(self.lCircular)
        self.polarization.addButton(self.linear)
        self.polarization.addButton(self.rCircular)
        self.polarization.addButton(self.lCircular)
        self.linear.setChecked(True)

        label4A=QtGui.QLabel("Theta (deg)")
        row4.addWidget(label4A)
        self.theta=QtGui.QLineEdit()
        self.theta.setText("0.0")
        row4.addWidget(self.theta)
        label4B=QtGui.QLabel("Phi (deg)")
        row4.addWidget(label4B)
        self.phi=QtGui.QLineEdit()
        self.phi.setText("0.0")
        row4.addWidget(self.phi)

        # Row 5: configuration (plot)
        row5=QtGui.QHBoxLayout()
        row5.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row5)

        self.dataToPlot=QtGui.QButtonGroup()
        self.plotDispersion=QtGui.QRadioButton("Band dispersion")
        self.plotPSF=QtGui.QRadioButton("Photoemission structure factor (PSF)")
        row5.addWidget(self.plotDispersion)
        row5.addWidget(self.plotPSF)
        self.dataToPlot.addButton(self.plotDispersion)
        self.dataToPlot.addButton(self.plotPSF)
        self.plotDispersion.setChecked(True)
        
        self.plotButton=QtGui.QPushButton("Plot")
        row5.addWidget(self.plotButton)


        # Row 6 left: graph (2D / 3D)
        row6=QtGui.QHBoxLayout()
        row6.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row6)

        self.graphTab=QtGui.QTabWidget()
        row6.addWidget(self.graphTab, 3)
        
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

        row6_3d=QtGui.QHBoxLayout()
        row6_3d.setAlignment(QtCore.Qt.AlignLeft)

        tabWidget=QtGui.QWidget()
        tabWidget.setLayout(row6_3d)
        self.graphTab.addTab(tabWidget, "3D")
        
        self.plot3D=gl.GLViewWidget()
        self.plot3D.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.plot3D.setMinimumSize(Config.plot3D_minWidth, Config.plot3D_minHeight)
        row6_3d.addWidget(self.plot3D)

        self.plot3=pg.GraphicsLayoutWidget()
        self.plot3.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        row6_3d.addWidget(self.plot3)

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
            Events.plot3(win, LCAO)
                
        self.plot.keyPressEvent=changeKBIndices
        self.plot3.keyPressEvent=changeKXYIndices

        # Row 6 right: properties
        vbox2=QtGui.QVBoxLayout()
        vbox2.setAlignment(QtCore.Qt.AlignTop)
        row6.addLayout(vbox2, 1)

        ## row 1: k index
        row6r1=QtGui.QHBoxLayout()
        row6r1.setAlignment(QtCore.Qt.AlignLeft)
        vbox2.addLayout(row6r1)
        
        label6r1A=QtGui.QLabel("kx index")
        row6r1.addWidget(label6r1A)
        self.kxIndex=QtGui.QSpinBox()
        self.kxIndex.setSingleStep(1)
        self.kxIndex.setMinimum(0)
        row6r1.addWidget(self.kxIndex)

        self.kxValue=QtGui.QLabel()
        row6r1.addWidget(self.kxValue)

        label6r1B=QtGui.QLabel("ky index")
        row6r1.addWidget(label6r1B)
        self.kyIndex=QtGui.QSpinBox()
        self.kyIndex.setSingleStep(1)
        self.kyIndex.setMinimum(0)
        row6r1.addWidget(self.kyIndex)

        self.kyValue=QtGui.QLabel()
        row6r1.addWidget(self.kyValue)

        ## row 2: band index and Constant energy index
        row6r2=QtGui.QHBoxLayout()
        row6r2.setAlignment(QtCore.Qt.AlignLeft)
        vbox2.addLayout(row6r2)
        
        label6r2A=QtGui.QLabel("Band index")
        row6r2.addWidget(label6r2A)
        self.bIndex=QtGui.QSpinBox()
        self.bIndex.setSingleStep(1)
        self.bIndex.setMinimum(0)
        row6r2.addWidget(self.bIndex)
        self.bValue=QtGui.QLabel()
        row6r2.addWidget(self.bValue)

        label6r2B=QtGui.QLabel("Energy index")
        row6r2.addWidget(label6r2B)
        self.eIndex=QtGui.QSpinBox()
        self.eIndex.setSingleStep(1)
        self.eIndex.setMinimum(0)
        row6r2.addWidget(self.eIndex)
        self.eValue=QtGui.QLabel()
        row6r2.addWidget(self.eValue)

        ## row 3: Up or Dn
        row6r3=QtGui.QHBoxLayout()
        row6r3.setAlignment(QtCore.Qt.AlignLeft)
        vbox2.addLayout(row6r3)

        self.UpDn=QtGui.QButtonGroup()
        self.UpButton=QtGui.QRadioButton("Up")
        self.DnButton=QtGui.QRadioButton("Dn")
        row6r3.addWidget(self.UpButton)
        row6r3.addWidget(self.DnButton)
        self.UpDn.addButton(self.UpButton)
        self.UpDn.addButton(self.DnButton)

        ## row 4: Atom
        row6r4=QtGui.QHBoxLayout()
        row6r4.setAlignment(QtCore.Qt.AlignLeft)
        vbox2.addLayout(row6r4)

        label6r4A=QtGui.QLabel("Atom")
        row6r4.addWidget(label6r4A)

        self.Atom=QtGui.QComboBox()
        self.Atom.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        row6r4.addWidget(self.Atom)

        ## row 5: LCAO table
        self.LCAOTable=QtGui.QTableWidget()
        vbox2.addWidget(self.LCAOTable)
        



app=QtGui.QApplication([])
# win: MainWindow object
win=MainWindow()
font=QtGui.QFont()
font.setPixelSize(Config.fontSize_normal)
font.setFamilies(Config.fontFamilies)
win.setFont(font)

LCAO=objs.LCAO()
win.openFileButton.clicked.connect(lambda: Events.openFile(win, LCAO))

def plotEvent():
    if LCAO.Dimension==1:
        Events.plot(win, LCAO)
        win.graphTab.setCurrentIndex(0)
    elif LCAO.Dimension==2:
        Events.makeDispersion3(win, LCAO)
        win.graphTab.setCurrentIndex(1)
    else:
        print("Dimension error")
        return

def cursorEvent():
    if LCAO.Dimension==1:
        Events.drawCursor(win, LCAO)
    elif LCAO.Dimension==2:
        Events.drawCursor3(win, LCAO)
    else:
        print("Dimension error")
        return

win.plotButton.clicked.connect(plotEvent)

win.kxIndex.valueChanged.connect(cursorEvent)
win.kyIndex.valueChanged.connect(cursorEvent)
win.bIndex.valueChanged.connect(cursorEvent)
win.eIndex.valueChanged.connect(cursorEvent)
win.UpButton.clicked.connect(cursorEvent)
win.DnButton.clicked.connect(cursorEvent)
win.Atom.currentIndexChanged.connect(lambda: Events.makeLCAOTable(win, LCAO))

pg.setConfigOptions(antialias=True)
win.show()
app.exec_()

        
