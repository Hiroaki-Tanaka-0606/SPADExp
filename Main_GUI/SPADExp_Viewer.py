# SPADExp Viewer

# PyQt5 and PyqQtGraph are necessary
from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import pyqtgraph.opengl as gl


import Config
from lib import Events_viewer as Events
from lib import objs

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setWindowTitle("SPADExp Viewer")

        
        font=QtGui.QFont()
        font.setFamilies(Config.fontFamilies)
        font.setPixelSize(Config.fontSize_normal)

        bFont=QtGui.QFont(font)
        bFont.setBold(True)
        
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

        # Row 6 left: real-space image
        row6=QtGui.QHBoxLayout()
        row6.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row6)

        row6l=QtGui.QVBoxLayout()
        row6l.setAlignment(QtCore.Qt.AlignTop)
        row6.addLayout(row6l)

        ## Row 1: Boundaries
        row6l1=QtGui.QHBoxLayout()
        row6l1.setAlignment(QtCore.Qt.AlignLeft)
        row6l.addLayout(row6l1)

        label6l1A=QtGui.QLabel("Boundaries (a, b, c)")
        label6l1A.setFont(bFont)
        row6l1.addWidget(label6l1A)
        self.boundaryA=QtGui.QSpinBox()
        self.boundaryA.setMinimum(1)
        self.boundaryA.setValue(1)
        row6l1.addWidget(self.boundaryA)
        self.boundaryB=QtGui.QSpinBox()
        self.boundaryB.setMinimum(1)
        self.boundaryB.setValue(1)
        row6l1.addWidget(self.boundaryB)
        self.boundaryC=QtGui.QSpinBox()
        self.boundaryC.setMinimum(1)
        self.boundaryC.setValue(1)
        row6l1.addWidget(self.boundaryC)

        ## Row 2: Weighting
        self.enableWeight=QtGui.QCheckBox("Enable weighting")
        row6l.addWidget(self.enableWeight)

        
        ## Row 3: plot button
        row6l3=QtGui.QHBoxLayout()
        row6l3.setAlignment(QtCore.Qt.AlignLeft)
        row6l.addLayout(row6l3)
        self.realSpacePlot=QtGui.QPushButton("Draw")
        row6l3.addWidget(self.realSpacePlot)
        
        ## Row 4: GLView
        row6l4=QtGui.QHBoxLayout()
        row6l4.setAlignment(QtCore.Qt.AlignLeft)
        row6l.addLayout(row6l4)
        self.realSpace=gl.GLViewWidget()
        self.realSpace.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.realSpace.setMinimumSize(Config.plot3D_minWidth, Config.plot3D_minHeight)
        row6l4.addWidget(self.realSpace)

        # Row 6 right: graph (2D / 3D)
        row6r=QtGui.QVBoxLayout()
        row6r.setAlignment(QtCore.Qt.AlignTop)
        row6.addLayout(row6r)

        ## Row6-1: indices button
        row6r1=QtGui.QHBoxLayout()
        row6r1.setAlignment(QtCore.Qt.AlignLeft)
        row6r.addLayout(row6r1)

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

        label6r2B=QtGui.QLabel("Energy index")
        row6r1.addWidget(label6r2B)
        self.eIndex=QtGui.QSpinBox()
        self.eIndex.setSingleStep(1)
        self.eIndex.setMinimum(0)
        row6r1.addWidget(self.eIndex)
        self.eValue=QtGui.QLabel()
        row6r1.addWidget(self.eValue)
        
        ## Row6-2: graph
        self.graphTab=QtGui.QTabWidget()
        row6r.addWidget(self.graphTab, 3)

        ## 2D
        self.plot=pg.PlotWidget()
        self.img=pg.ImageItem()
        # row3.addWidget(self.plot)
        self.plot.addItem(self.img)
        self.graphTab.addTab(self.plot, "2D")

        cmap=pg.colormap.get("CET-L9")
        self.bar=pg.ColorBarItem(colorMap=cmap)
        self.bar.setImageItem(self.img)
                            

        
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
        self.plot.setLabel(axis="left", text="Wavevector (a.u.^-1)")

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
        self.barEx=pg.ColorBarItem(colorMap=cmap, values=(0,1))
        self.barEy=pg.ColorBarItem(colorMap=cmap, values=(0,1))
        self.barxy=pg.ColorBarItem(colorMap=cmap, values=(0,1))
        self.barEx.setImageItem(self.imgEx)
        self.barEy.setImageItem(self.imgEy)
        self.barxy.setImageItem(self.imgxy)
        
        self.plotEx.getAxis("bottom").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotEx.getAxis("left").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotEx.getAxis("bottom").setPen((255,255,255))
        self.plotEx.getAxis("left").setPen((255,255,255))
        self.plotEx.getAxis("bottom").setTextPen((255,255,255))
        self.plotEx.getAxis("left").setTextPen((255,255,255))
        self.plotEx.getAxis("bottom").setLabel(**labelStyle)
        self.plotEx.getAxis("left").setLabel(**labelStyle)
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
        
        self.plotEx.setLabel(axis="left", text="E-EF (eV)")
        self.plotEx.setLabel(axis="bottom", text="kx (a.u.^-1)")
        
        self.plotEy.setLabel(axis="bottom", text="E-EF (eV)")
        self.plotEy.setLabel(axis="left", text="ky (a.u.^-1)")

        self.plotxy.setLabel(axis="left", text="ky (a.u.^-1)")
        self.plotxy.setLabel(axis="bottom", text="kx (a.u.^-1)")

        
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
            else:
                return
            
        self.plot3.keyPressEvent=changeKXYIndices
    

app=QtGui.QApplication([])
# win: MainWindow object
win=MainWindow()
font=QtGui.QFont()
font.setPixelSize(Config.fontSize_normal)
font.setFamilies(Config.fontFamilies)
win.setFont(font)
Elements=objs.Elements()
Disp=objs.Dispersion()

win.openFileButton.clicked.connect(lambda: Events.openFile(win, Disp, Elements))
win.realSpacePlot.clicked.connect(lambda: Events.makeRealSpace(win, Disp, Elements))

def cursorEvent():
    if Disp.Dimension==2:
        Events.plot3(win, Disp)
        Events.drawCursor3(win, Disp)
    elif Disp.Dimension==1:
        pass
    else:
        print("Dimension error")
        return

win.kxIndex.valueChanged.connect(cursorEvent)
win.kyIndex.valueChanged.connect(cursorEvent)
win.eIndex.valueChanged.connect(cursorEvent)

pg.setConfigOptions(antialias=True)
win.show()
app.exec_()

        
