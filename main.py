# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:24:56 2014

@author: jzb4

Note:  gdspy generally doesn't like spyder.  Make sure you run this with the
setting that it runs in it's own, new, python interpreter.
(I think python works better than ipython, as well.)

Things you need:
gdspy installed.
python files:
- gdspy_subroutines
- junctions
- qubits
- resonators
- main (this file)


"""

import gdspy


from gdspy_subroutines import *
import qubits
import junctions


device_number = 1

params = {'finger_width' : 450,
          'antenna_funnel_length': 20*um,
          'junction_type': 'dolan',
          'device_number':device_number}


qubit_cell = gdspy.Cell('device[%d]'%device_number)


plgs = qubits.transmon(params)
#plgs.add(two_inch_wafer())
wafer_plgs = two_inch_wafer()

add_plgs(qubit_cell, plgs)

gdspy.LayoutViewer(pattern={'default':7})
