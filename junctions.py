# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:24:06 2014

@author: jzb4
"""

import gdspy

from gdspy_subroutines import *

'''These are the fault parameters.  They'll be overwritten if the passed params
dict contains entries with the same key.'''

dolan_ps = {  'finger_width' : 400,
              'bridge_width' : 500,
              'bridge_height' : 600,
              'needle_width' : 200,
              'finger_length' : 10000,
              'undercut_layer': 3,
              'fine_layer' : 2,
              'inner_trace_width': 1000,
              'inner_trace_length':10000,
              'final_trace_width':10000,
              'inner_funnel_length':500,
              'outer_funnel_length':2000,
              'fine_coarse_overlap':200}

def junction(params):
    #plgs = poly_list([])  #Polygons to be returned
    junction_length = 0  #Length of the returned polygon

    if params['junction_type'] == 'dolan':
        dolan_ps.update(params)


        '''Some convenient shorthands:
        Bridge is the undercut region, touched by finger and needle with given widths.
        The finger width defines the resistance.  They extend out finger and needle =
        length, which are here the same, and funnel out to final_trace_width,
        which is often one micron for some length, then out again to a final trace width '''

        fw = dolan_ps['finger_width']
        fl = dolan_ps['finger_length']
        bw = dolan_ps['bridge_width']
        bh = dolan_ps['bridge_height']
        nw = dolan_ps['needle_width']

        itw = dolan_ps['inner_trace_width']
        itl = dolan_ps['inner_trace_length']
        ftw = dolan_ps['final_trace_width']

        ifl = dolan_ps['inner_funnel_length']
        ofl = dolan_ps['outer_funnel_length']

        ul = dolan_ps['undercut_layer']
        ml = dolan_ps['fine_layer'] #Main layer

        fco = dolan_ps['fine_coarse_overlap']

        bridge = gdspy.Rectangle((-bw/2, -bh/2), (bw/2, bh/2), layer=ul)

        finger = gdspy.Rectangle((0,0),(fl,fw), layer =ml)
        translate(finger,(bw/2,-fw/2))

        needle = gdspy.Rectangle((0,0),(-fl,nw), layer =ml)
        translate(needle,(-bw/2,-nw/2))

        #plgs.extend([bridge, finger, needle])
        plgs = poly_list([bridge, finger, needle])

        finger_funnel = funnel(fw,itw,layer = ml,length = ifl)
        translate(finger_funnel,(bw/2+fl,0))

        needle_funnel = funnel(nw,itw,layer = ml,length = ifl)
        needle_funnel.rotate(np.pi)
        translate(needle_funnel,(-bw/2-fl,0))

        #plgs.extend([finger_funnel, needle_funnel])
        plgs.add([finger_funnel, needle_funnel])

        inner_trace_fs = gdspy.Rectangle((bw/2+fl+ifl,-itw/2),(bw/2+fl+itl,itw/2),layer=ml)
        inner_trace_ns = gdspy.Rectangle((-(bw/2+fl+ifl),-itw/2),(-(bw/2+fl+itl),itw/2),layer=ml)

        #plgs.extend([inner_trace_fs,inner_trace_ns])
        plgs.add([inner_trace_fs,inner_trace_ns])

        outer_funnel_fs = funnel(itw,ftw,layer = ml,length = ofl)
        translate(outer_funnel_fs,(bw/2+fl+itl,0))

        outer_funnel_ns = funnel(itw,ftw,layer = ml,length = ofl)
        outer_funnel_ns.rotate(np.pi)
        translate(outer_funnel_ns,(-(bw/2+fl+itl),0))

        plgs.add([outer_funnel_fs, outer_funnel_ns])

        #plgs.extend([outer_funnel_fs, outer_funnel_ns])

        #add fine_coarse overlap!!
        #fco =

        junction_length = bw + 2*fl + 2*itl +2*ofl


    if params['junction_type'] == 'dolan_squid':
        pass

    if params['junction_type'] == 'bridge_free':
        pass

    if params['junction_type'] == 'bridge_free_squid':
        pass



    return plgs, junction_length



