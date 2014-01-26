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
wafer_name = 'AA15'

dice_kerf = 250*um
params = {'finger_width' : 450,
          'antenna_funnel_length': 20*um,
          'dice_width':3.2*mm,
          'dice_length':6.7*mm,
          'junction_type': 'dolan',
          'device_number':device_number,
          'wafer_name':wafer_name}


qubit_cell = gdspy.Cell('%s[%d]' % (wafer_name,device_number))


#plgs = qubits.transmon(params)
plgs = two_inch_wafer()

idcs = populate_wafer(params['dice_length']+dice_kerf,params['dice_width']+dice_kerf)

dx = params['dice_length']
dy = params['dice_width']
dxp = dx+dice_kerf
dyp = dy+dice_kerf

for idx in idcs:
    plg = gdspy.Rectangle((-dx/2,-dy/2),(dx/2,dy/2))
    translate(plg,(idx[0]*dxp,idx[1]*dyp))
    plgs.add(plg)

add_plgs(qubit_cell, plgs)
gdspy.LayoutViewer(pattern={'default':7})
