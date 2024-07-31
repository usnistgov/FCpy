# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 20:38:45 2022

@author: Evan Groopman
@affiliation: National Institute of Standards and Technology (NIST)
@email: evan.groopman@nist.gov
"""

import os
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import Qt
from PyQt5 import uic
import time
import numpy as np

#local import
try:
    from FCpy.FCpy import FC
except ModuleNotFoundError:
    try:
        import FC
    except ModuleNotFoundError:
        from FCpy import FC

#convenience loader function
def loadUI_N(parentPath, relativePath, parent= None):
    """parentPath will be the full path name of __file__"""
    directory = os.path.dirname(parentPath)
    path = os.path.join(directory, relativePath)
    uic.loadUi(path, parent)

#results table for GUI
class TableModel(QtCore.QAbstractTableModel):
    """simple table model"""
    def __init__(self, data, header, index):
        super(TableModel, self).__init__()
        self._data = data #numpy array
        self._header = header
        self._index = index

    def data(self, index, role):
        if role == Qt.DisplayRole:
            if isinstance(self._data, list):
                value = self._data[index.row()][index.column()]
            elif isinstance(self._data, np.ndarray):
                value = self._data[index.row(), index.column()]
            return str(value)

    def rowCount(self, index):
        return np.shape(self._data)[0]

    def columnCount(self, index):
        return np.shape(self._data)[1]

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self._header[section])

            if orientation == Qt.Vertical:
                return str(self._index[section])
    
    def setData(self, index, value, role= Qt.DisplayRole):
        if isinstance(self._data, list):
            self._data[index.row()][index.column()] = value
        elif isinstance(self._data, np.ndarray):
            self._data[index.row(), index.column()] = value
        return True

    def updateData(self, data):
        self._data = data
        self.dataChanged.emit(QtCore.QModelIndex(), QtCore.QModelIndex())
    
    def updateHeader(self, header):
        self._header = header
        self.headerDataChanged.emit(Qt.Horizontal, 0, 0)
        
    def updateIndex(self, index):
        self._index = index
        self.headerDataChanged.emit(Qt.Vertical, 0, 0)

    
class FC_GUI(QtWidgets.QMainWindow):
    
    def __init__(self, parent=None, default_dir= None):
        super(FC_GUI, self).__init__(parent)
        loadUI_N(__file__, "UI" + os.path.sep +"FC_UI.ui", self)
        
        self._signals_and_slots()
        
        #add dummy results to the results table to start
        results = np.array([[0,0],[0,0],[0,0]])
        self.resultsModel = TableModel(results, ['CI Low', 'CI High'], ['68.3%', '95%', '99.7%'])
        self.resultsTable.setModel(self.resultsModel)
        

    def _signals_and_slots(self):
        self.b_calculate.clicked.connect(self._calculate)
    
    def _calculate(self):
        N = int(self.sb_n0.value())
        b = float(self.sb_background.value())
        t = float(self.sb_time.value())
        corr = self.cb_correction.isChecked()
        manualConf = self.cb_manualConf.isChecked()
        conf = float(self.sb_conf.value())
        
        go = True
        if corr:
            # if N is high and corr checked, pop up verification box
            if N > 35:
                mbox = QtWidgets.QMessageBox
                ret = mbox.question(self,'Confirm', 'N=%d calculation with\nthe correction\could take a long time.\nAre you sure?' % N, mbox.Yes | mbox.No)
                if ret == mbox.No:
                    go = False
                    if manualConf:
                        self.resultsModel.updateData(np.array([[0,0]]))
                        self.resultsModel.updateIndex(['%.1f' % round(conf*100,1) + '%'])
                    else:
                        self.resultsModel.updateData(np.array([[0,0],[0,0],[0,0]]))
                        self.resultsModel.updateIndex(['68.3%', '95%', '99.7%'])
        
        if go:
            if manualConf:
                startTime = time.time()
                CIlow, CIhigh = FC.FC_poisson(N, b, t, conf=conf, useCorrection=corr) #use defualt tolerance 5E-4
                calcTime = (time.time() - startTime)*1000
                self.resultsModel.updateData(np.array([[CIlow, CIhigh]]).round(3))
                #find number of decimals to display
                
                value = round(conf*100,3)
                if value == int(value):
                    form = '%d'
                    value = int(conf*100)
                else:
                    ndigits = len(str(value).split('.')[-1])
                    form = '%.' + '%d' % ndigits + 'f'
                    value = round(conf*100,ndigits)

                self.resultsModel.updateIndex([form % value + '%'])
                
            else:
                startTime = time.time()
                CIs = FC.FC_poisson123(N, b, t, useCorrection=corr) #use defualt tolerance
                calcTime = (time.time() - startTime)*1000
                
                self.resultsModel.updateData(np.array(CIs).round(3))
                self.resultsModel.updateIndex(['68.3%', '95%', '99.7%'])
            
            msg = 'Calculation time: %.1f ms' % calcTime
            self.statusbar.showMessage(msg)

if __name__ == '__main__':
    import platform
    import sys
    
    # for launching from ipython
    try:
        if __IPYTHON__:
            ipython = get_ipython()
            ipython.magic('reset_selective -f win') #clear any previous invocations of the program
    except NameError:
        ipython = None
    
        
    #this conditional is needed for Spyder to launch in IPython
    #prevents kernel from dying... strange bug, but it works.
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance()
    app.setOrganizationName("National Institute of Standards and Technology")
    app.setOrganizationDomain("evan.groopman@nist.gov")
    app.setApplicationName("FC Poisson Calculations")
    win = FC_GUI()
    win.show()
    
    #for starting on macs
    if ipython is None or platform.platform().startswith('Darwin'):
        sys.exit(app.exec_())