# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 17:58:54 2014

@author: rsl
"""

import gdspy
reload(gdspy)
from gdspy_subroutines import *
import qubits
reload(qubits)
import junctions

wafer_name = 'test1'

layers = {'coarse_layer':1,
          'fine_layer':2,
          'undercut_layer':3,
          'label_layer':14,
          'fine_floating_fields_layer':6,
          'coarse_floating_fields_layer':5}


'''first the junctions'''
dx = 0.5*mm
dy = 0.5*mm

def params_width(width):
    params= {'finger_width' : width,
          'antenna_funnel_length': 20*um,
          'antenna_length':50*um,
          'dice_width':dy,
          'dice_length':dx,
          'junction_type': 'dolan',
          'wafer_name':wafer_name,
          'test_devices':False,
          'litho_label':False,
          'gds_label':False,
          'short':False}
    params.update(layers)
    return params

max_width = 400
min_width = 100
num_widths = 31

num_reps = 5

param_list = [params_width(width) for width in np.linspace(min_width,
                                                           max_width,
                                                           num_widths)]*num_reps

'''Draw the wafer'''
gdspy.Cell.cell_dict.clear()
main_cell = gdspy.Cell('%s' % (wafer_name,))

plgs = two_inch_wafer()
add_plgs_to_cell(main_cell, plgs)

'''Figure out our dice locations'''
yoffset = 30
yspacing = 1.5

xs = np.arange(-num_widths,num_widths,2)+1
ys = yoffset-np.arange(num_reps)*yspacing
#location_idcs = [(x,35-y*2) for y in range(num_reps) for x in xs ]
location_idcs = [(x,y) for y in ys for x in xs ]
'''Draw the qubits '''
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
                                   origin=(location_idx[0]*dx,
                                           (location_idx[1])*dy))
    add_plgs_to_cell(main_cell, cell_ref)

    '''Add some labels'''
    if location_idx[1] == max(ys):
        lbl = text_plg('%d' % (params['finger_width'],),
                       size=100*um,
                       layer=layers['coarse_layer'])
        translate(lbl,((location_idx[0]-0.3)*dx,(location_idx[1]+1)*dy))
        add_plgs_to_cell(main_cell,lbl)

    if location_idx[1] == min(ys):
        lbl = text_plg('%d' % (params['finger_width'],),
                       size=100*um,
                       layer=layers['coarse_layer'])
        translate(lbl,((location_idx[0]-0.3)*dx,(location_idx[1]-1)*dy))
        add_plgs_to_cell(main_cell,lbl)

    if location_idx[0] == min(xs):
        lbl = text_plg('%d' % (-(location_idx[1]-yoffset)/yspacing,),
                       size=100*um,
                       layer=layers['coarse_layer'])
        translate(lbl,((location_idx[0]-1)*dx,(location_idx[1])*dy))
        add_plgs_to_cell(main_cell,lbl)

    if location_idx[0] == max(xs):
        lbl = text_plg('%d' % (-(location_idx[1]-yoffset)/yspacing,),
                       size=100*um,
                       layer=layers['coarse_layer'])
        translate(lbl,((location_idx[0]+1)*dx,(location_idx[1])*dy))
        add_plgs_to_cell(main_cell,lbl)

#==============================================================================
# Resonator example
#==============================================================================

save_gds(wafer_name)
#gdspy.gds_print('thing.gds',unit=1.0e-9, precision=1.0e-9)
#gdspy.LayoutViewer(pattern={'default':7})