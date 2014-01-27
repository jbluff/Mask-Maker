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

wafer_name = 'AA15'

dice_kerf = 250*um
dx = 6.7*mm
dy = 3.2*mm

layers = {'coarse_layer':1,
          'fine_layer':2,
          'undercut_layer':3,
          'label_layer':14,
          'fine_floating_fields_layer':6,
          'coarse_floating_fields_layer':5}

params = {'finger_width' : 450,
          'antenna_funnel_length': 20*um,
          'dice_width':dy,
          'dice_length':dx,
          'junction_type': 'dolan',
          'wafer_name':wafer_name}

params.update(layers)
param_list = [params,]*51

'''Draw the wafer'''
main_cell = gdspy.Cell('%s' % (wafer_name,))
plgs = two_inch_wafer()
add_plgs_to_cell(main_cell, plgs)

'''Figure out our dice locations'''
dxp = dx + dice_kerf
dyp = dy + dice_kerf
location_idcs = populate_wafer(dxp,dyp)


'''Draw the qubits (dices)'''
for device_number,params in enumerate(param_list):
    #print '%d/%d' % (device_number, len(param_list))
    params['device_number'] = device_number

    '''Draw the qubit in its own cell'''
    qubit_cell = gdspy.Cell('%s[%d]' % (wafer_name,device_number))
    plgs = qubits.transmon(params)
    add_plgs_to_cell(qubit_cell, plgs)

    '''Add the qubit cell to the main cell'''
    location_idx = location_idcs[device_number]
    cell_ref = gdspy.CellReference(qubit_cell,
                                   origin=(location_idx[0]*dxp,
                                           location_idx[1]*dyp))
    add_plgs_to_cell(main_cell, cell_ref)


'''Draw row labels, dice marks'''
label_plgs = dice_marks_and_labels(location_idcs,dxp,dyp,params['coarse_layer'])
add_plgs_to_cell(main_cell, label_plgs)


'''Add the plassys alignment aid'''
pam = plassys_alignment_mark(params['coarse_layer'])
add_plgs_to_cell(main_cell,pam)

save_gds(wafer_name)

#gdspy.LayoutViewer(pattern={'default':7})
