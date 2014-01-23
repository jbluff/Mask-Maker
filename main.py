# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:24:56 2014

@author: jzb4
"""

import gdspy


from gdspy_subroutines import * 
from qubits import *

params = {'finger_width' : 450,
          'antenna_funnel_length': 20*um}

device_number = 1
qubit_cell = gdspy.Cell('device[%d]'%device_number)

plgs = transmon(params)

add_plgs(qubit_cell, plgs)



gdspy.LayoutViewer()


