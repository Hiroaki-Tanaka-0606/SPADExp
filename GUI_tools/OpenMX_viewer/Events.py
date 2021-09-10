# Events

from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import numpy as np
import re

def openFile(win):
    currentFile=win.filePath.text()
    selectedFile, _filter=QtGui.QFileDialog.getOpenFileName(caption="Open file", directory=currentFile)
    if selectedFile!="":
        # valid
        win.filePath.setText(selectedFile)

        in_the_block=False
        currentBlock=""
        currentBlockData=[]
        currentMaxItems=0
        line_number=0
        global blockNames
        global blockData
        blockNames=[]
        blockData=[]
        with open(selectedFile, "r") as f:
            while True:
                line=f.readline()
                if not line:
                    break
                line_number+=1
                # find a row starting from '<' (beginning of a block)
                re_result=re.findall(r"^<(\S+)", line)
                if len(re_result)>0:
                    if in_the_block:
                        print(("Error in row {0:d}: a block started in the block").format(line_number))
                        break
                    if re_result[0][-1]==">":
                        continue
                    in_the_block=True
                    currentBlock=re_result[0]
                    blockNames.append(currentBlock)
                    currentBlockData=[]
                    currentMaxItems=0
                    continue

                re_result=re.findall("^"+currentBlock+">", line)
                if len(re_result)>0:
                    in_the_block=False
                    validItems=0
                    for lineData in currentBlockData:
                        if len(lineData)==currentMaxItems:
                            validItems+=1
                    np_data=np.zeros((currentMaxItems, validItems))
                    i=0
                    j=0
                    for lineData in currentBlockData:
                        j=0
                        if len(lineData)==currentMaxItems:
                            for value in lineData:
                                np_data[j][i]=value
                                j+=1
                            i+=1
                    blockData.append(np_data)
                    continue

                if in_the_block:
                    items=line.split()
                    currentLineData=[]
                    for item in items:
                        item_float=0.0
                        try:
                            item_float=float(item)
                            currentLineData.append(item_float)
                        except:
                            pass
                    if currentMaxItems<len(currentLineData):
                        currentMaxItems=len(currentLineData)
                    currentBlockData.append(currentLineData)

            # update datablock list                
            win.selectDataBlock.clear()
            for blockName in blockNames:
                win.selectDataBlock.addItem(blockName)
                
def setDataBlock(win):
    index=win.selectDataBlock.currentIndex()

    # X axis
    while len(win.xAxis.buttons())>0:
        win.xCol.removeWidget(win.xAxis.buttons()[0])
        
    for i in range(len(blockData[index])):
        xButton=QtGui.QRadioButton(str(i+1))
        win.xAxis.addButton(xButton)
        win.xCol.addWidget(xButton)

    # Y axes
    while len(win.yAxes.buttons())>0:
        win.yCol.removeWidget(win.yAxes.buttons()[0])

    for i in range(len(blockData[index])):
        yButton=QtGui.QCheckBox(str(i+1))
        win.yAxes.addButton(yButton)
        win.yCol.addWidget(yButton)


def setGraph(win):
    blockIndex=win.selectDataBlock.currentIndex()
    xIndex=-1
    if win.xAxis.checkedButton()!=None:
        xIndex=int(win.xAxis.checkedButton().text())-1
        
    yIndices=[]
    for yButton in win.yAxes.buttons():
        if yButton.isChecked():
            yIndices.append(int(yButton.text())-1)

    if xIndex==-1 or len(yIndices)<1:
        return
        
    win.graphPlot.clear()
    for yIndex in yIndices:
        win.graphPlot.plot(y=blockData[blockIndex][yIndex], x=blockData[blockIndex][xIndex])