# -*- coding: utf-8 -*-
"""
NIST-developed software is provided by NIST as a public service. You may use, 
copy, and distribute copies of the software in any medium, provided that you 
keep intact this entire notice. You may improve, modify, and create derivative 
works of the software or any portion of the software, and you may copy and 
distribute such modifications or works. Modified works should carry a notice 
stating that you changed the software and should note the date and nature of 
any such change. Please explicitly acknowledge the National Institute of 
Standards and Technology as the source of the software. 

NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY 
OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, 
INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. 
NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL 
BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST 
DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE 
OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, 
RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

You are solely responsible for determining the appropriateness of using and 
distributing the software and you assume all risks associated with its use, 
including but not limited to the risks and costs of program errors, compliance 
with applicable laws, damage to or loss of data, programs or equipment, and the
unavailability or interruption of operation. This software is not intended to 
be used in any situation where a failure could cause risk of injury or damage 
to property. The software developed by NIST employees is not subject to 
copyright protection within the United States.

@author: Evan Groopman
@affiliation: National Institute of Standards and Technology
@email: evan.groopman@nist.gov
"""
# import multiprocessing
import os, sys
from PyQt5 import QtWidgets, QtCore, QtGui

import FCpy

try:
    #for running from command line or pyinstaller executable
    from FCpy.FC_GUI import FC_GUI
except ModuleNotFoundError:
    try:
        #for running in IPython (e.g., Spyder)
        from FCpy.FCpy.FC_GUI import FC_GUI
    except ModuleNotFoundError:
        from FC_GUI import FC_GUI


try:
    from ctypes import windll  # Only exists on Windows.
    myappid = 'NIST.FCpy.04112024'
    windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
except ImportError:
    pass

if __name__ == '__main__':
    # multiprocessing.freeze_support()
    import platform
    
    try:
        if __IPYTHON__:
            ipython = get_ipython()
            ipython.run_line_magic('reset_selective', '-f win') #clear any previous invocations of the program
    except NameError:
        ipython = None
    
        
    #this conditional is needed for Spyder to launch in IPython
    #prevents kernel from dying... strange bug, but it works.
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance()
    
    app.setOrganizationName("National Institute of Standards and Technology")
    app.setOrganizationDomain("evan.groopman@nist.gov")
    app.setApplicationName("FC CI Calculator")
    #app.setWindowIcon(QtGui.QIcon(os.path.join(os.path.dirname(__file__), "APMview/UI/icons/icon.ico")))
    win = FC_GUI()
    win.show()
    
    #for starting on macs
    if ipython:
        pass
    elif ipython is None or platform.platform().startswith('Darwin'):
        sys.exit(app.exec_())
    else:
        app.exec_()