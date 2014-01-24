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
                'coarse_layer':1,
                'undercut_layer':3,
                'label_layer':14,
                'final_trace_width':10000,
                'antenna_funnel_length':10*um,
                'edge_rounding':10*um
                }

                 #Later we'll deal with rotation angle with beziers.
def transmon(params):

    plgs = poly_list()

    transmon_ps.update(params)

    junction_plgs, junction_length = junctions.junction(transmon_ps)
    #plgs.extend(junction_plgs)
    plgs.add(junction_plgs)

    sep = transmon_ps['separation']
    al = transmon_ps['antenna_length']
    aw = transmon_ps['antenna_width']
    dw = transmon_ps['dice_width']
    dl = transmon_ps['dice_length']
    cl = transmon_ps['coarse_layer']
    ul = transmon_ps['undercut_layer']
    ll = transmon_ps['label_layer']
    tw = transmon_ps['final_trace_width']
    afl = transmon_ps['antenna_funnel_length']
    er = transmon_ps['edge_rounding']
    jl = junction_length

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

    if dl and dw:
        plgs.add(gdspy.Rectangle((-dl/2,-dw/2),(dl/2,dw/2),layer=10))

        #Delete this.
        #plgs = poly_list(gdspy.Rectangle((-dl/2,-dw/2),(dl/2,dw/2),layer=10))
    return plgs

def rotated_transmon(params):
    pass


def fluxnonium(params):
    pass

