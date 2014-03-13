# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 17:10:51 2014

@author: jzb4
"""

import gdspy
import junctions
import copy
from gdspy_subroutines import *

'''Default parameters for transmons'''

transmon_ps = {'junction_type' : 'dolan',
                'separation': 100*um,
                'antenna_length': 500*um,
                'antenna_width': 250*um,
                'dice_width':3.2*mm,
                'dice_length':6.7*mm,
                'final_trace_width':10000,
                'antenna_funnel_length':10*um,
                'edge_rounding':10*um,
                'short':True,

                'coarse_layer':1,
                'undercut_layer':3,
                'label_layer':14,
                'fine_floating_field_layer':6,

                'test_devices':True,
                'litho_label':True,
                'gds_label':True,

                'just_qubit': False #useful for HFSS
                }

                 #Later we'll deal with rotation angle with beziers.
def transmon(params):

    plgs = poly_list()

    '''Update our parameter array'''
    transmon_ps.update(params)

    junction_plgs, junction_length,ret_ps = junctions.junction(transmon_ps)
    plgs.add(junction_plgs)

    sep = ret_ps['separation']
    al = ret_ps['antenna_length']
    aw = ret_ps['antenna_width']

    tw = ret_ps['final_trace_width']
    afl = ret_ps['antenna_funnel_length']
    er = ret_ps['edge_rounding']

    dw = ret_ps['dice_width']
    dl = ret_ps['dice_length']

    cl = ret_ps['coarse_layer']
    ul = ret_ps['undercut_layer']
    fl = ret_ps['fine_layer']
    ll = ret_ps['label_layer']
    fffl = ret_ps['fine_floating_field_layer']

    just_qubit = ret_ps['just_qubit']
    td = ret_ps['test_devices']
    lithl = ret_ps['litho_label']
    gl = ret_ps['gds_label']

    jl = junction_length

    '''Draw the main antenna'''
    left_trace = gdspy.Rectangle((+jl/2,tw/2),(sep/2-afl/2,-tw/2),layer=cl)
    right_trace = gdspy.Rectangle((+jl/2,tw/2),(sep/2-afl/2,-tw/2),layer=cl)
    right_trace.rotate(np.pi)

    plgs.add([left_trace,right_trace])

    left_funnel = funnel(tw,aw,layer = cl,length = afl)
    translate(left_funnel, (sep/2-afl/2,0))
    right_funnel = funnel(tw,aw,layer = cl,length = afl)
    right_funnel.rotate(np.pi)
    translate(right_funnel, (-sep/2+afl/2,0))

    plgs.add([left_funnel,right_funnel])

    right_plgs = []
    right_pad = gdspy.Rectangle((sep/2+afl/2,aw/2),(sep/2+al-er,-aw/2),layer=cl)
    right_upper_corner = round_plg((sep/2+al-er,aw/2-er),er,
                                     initial_angle = 0,
                                     final_angle = np.pi/2,layer=cl)

    right_lower_corner = copy.deepcopy(right_upper_corner)
    xflip(right_lower_corner)
    right_edge_plg = gdspy.Rectangle((sep/2+al-er,aw/2-er),(sep/2+al,-aw/2+er),layer=cl)
    right_plgs = [right_pad, right_upper_corner,right_lower_corner,right_edge_plg]

    left_plgs = copy.deepcopy(right_plgs)
    yflip(left_plgs)

    plgs.add(right_plgs)
    plgs.add(left_plgs)

    if ret_ps['short']:
        plgs.add(gdspy.Rectangle((-afl-sep/2,aw/2),(afl+sep/2,aw/2-2*um),layer=cl))

    if not just_qubit:
        '''Dice size indicator'''
        plgs.add(gdspy.Rectangle((-dl/2,-dw/2),(dl/2,dw/2),layer=ll))

        '''Floating field'''
        ffs = 500*um #Size of the floating field box.
        plgs.add(gdspy.Rectangle((-ffs/2,-ffs/2),(ffs/2,ffs/2),layer=fffl))

        '''Test junctions.
        Three test junctions lets you see visual if the device is upside
        down, so I'm a fan of that configuration.'''

        if td:
            tj_height = 250*um
            tj_width = 250*um
            tj_xdim = (-tj_width/2,-tj_height/2)
            tj_ydim = (tj_width/2,tj_height/2)

            tj_bounding_box = poly_list([
                gdspy.Rectangle(tj_xdim,tj_ydim,layer=cl),
                gdspy.Rectangle(tj_xdim,tj_ydim,layer=fl),
                gdspy.Rectangle(tj_xdim,tj_ydim,layer=ul),
                gdspy.Rectangle(tj_xdim,tj_ydim,layer=fffl)
                                ])

            tj_plgs = plg_bool(plgs,tj_bounding_box,'int')
            tj_plgs2 = copy.deepcopy(tj_plgs)
            tj_plgs3 = copy.deepcopy(tj_plgs)

            tj_x1 = -dl/2 + tj_width/2 + 100*um
            tj_x3 = dl/2 - tj_width/2 - 100*um
            tj_y1 = -dw/2 + tj_height/2 + 100*um
            tj_y2 = dw/2 - tj_height/2 - 100*um

            translate(tj_plgs,(tj_x1,tj_y1))
            translate(tj_plgs2,(tj_x1,tj_y2))
            translate(tj_plgs3,(tj_x3,tj_y2))

            plgs.add(tj_plgs)
            plgs.add(tj_plgs2)
            plgs.add(tj_plgs3)


        '''Litho label  (actually printed)'''
        if lithl:
            label_text = '%s[%s]'%(ret_ps['wafer_name'],ret_ps['device_number'])
            size = 150*um
            text_plgs = text_plg(label_text,size=size,angle=np.pi/2,layer=cl)

            text_bb = bounding_box(text_plgs)
            y_center = text_bb[1][1] - text_bb[0][1]

            translate(text_plgs,(-dl/2+size+150*um,-y_center/2))
            plgs.add(text_plgs)


        '''GDS label (not actually exposed.)
        I generally keep the finger width printed on the GDS for reference
        purposes.  Any other swept variable also makes sense.'''
        if gl:
            label_text = '%dnm'%(ret_ps['finger_width'])
            size = 250*um
            text_plgs = text_plg(label_text,size=size,layer=ll)

            text_bb = bounding_box(text_plgs)
            y_center = text_bb[1][1] - text_bb[0][1]

            translate(text_plgs,(-dl/3.0+size,-y_center/2))
            plgs.add(text_plgs)



        #print y_center

        #plgs.add(gdspy.Rectangle(text_bb[0],text_bb[1]))

    return plgs

def parallel_capacitor_transmon(params):
    #These have a better name.
    pass

def rotated_transmon(params):
    pass


def fluxnonium(params):
    pass

