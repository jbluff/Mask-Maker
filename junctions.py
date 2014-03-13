# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:24:06 2014

@author: jzb4
"""

import gdspy

from gdspy_subroutines import *

'''These are the default parameters.  They'll be overwritten if the passed params
dict contains entries with the same key.'''

dolan_ps = {  'finger_width' : 400,
              'bridge_extra_width' : 200,
              'bridge_length' : 460,
              'finger_extra_width' : 200,
              'finger_length' : 4000,
              'undercut_layer': 3,
              'fine_layer' : 2,
              'inner_trace_width': 1000,
              'inner_trace_length':9500,
              'final_trace_width':10000,
              'inner_funnel_length':500,
              'outer_funnel_length':2000,
              'fine_coarse_overlap':200}

def junction(params):
    plgs = poly_list()  #Polygons to be returned
    junction_length = 0  #Length of the returned plgs (to align qubit)

    if params['junction_type'] == 'dolan':
        dolan_ps.update(params)


        '''Some convenient shorthands:
        Bridge is the undercut region, touched by finger and needle with given widths.
        The finger width defines the resistance.  They extend out finger and needle =
        length, which are here the same, and funnel out to final_trace_width,
        which is often one micron for some length, then out again to a final trace width '''

        fw = dolan_ps['finger_width']
        fl = dolan_ps['finger_length']
        bew = dolan_ps['bridge_extra_width']
        bl = dolan_ps['bridge_length']
        few = dolan_ps['finger_extra_width']

        itw = dolan_ps['inner_trace_width']
        itl = dolan_ps['inner_trace_length']
        ftw = dolan_ps['final_trace_width']

        ifl = dolan_ps['inner_funnel_length']
        ofl = dolan_ps['outer_funnel_length']

        ul = dolan_ps['undercut_layer']
        ml = dolan_ps['fine_layer'] #Main layer

        fco = dolan_ps['fine_coarse_overlap']

        bridge = gdspy.Rectangle((-bl/2,-(fw+few+bew)/2),
                                 (bl/2,(fw+few+bew)/2),
                                 layer=ul)

        finger = gdspy.Rectangle((0,0),(-fl,fw), layer =ml)
        translate(finger,(-bl/2,-fw/2))

        wider_finger = gdspy.Rectangle((0,0),(fl,fw+few), layer =ml)
        translate(wider_finger,(bl/2,-(fw+few)/2))

        plgs.add([bridge, finger, wider_finger])

        finger_funnel = funnel(fw,itw,layer = ml,length = ifl)
        finger_funnel.rotate(np.pi)
        translate(finger_funnel,(-bl/2-fl,0))

        wider_funnel = funnel(fw+few,itw,layer = ml,length = ifl)
        translate(wider_funnel,(+bl/2+fl,0))

        plgs.add([finger_funnel, wider_funnel])

        inner_trace_fs = gdspy.Rectangle((bl/2+fl+ifl,-itw/2),(bl/2+fl+itl,itw/2),layer=ml)
        inner_trace_ns = gdspy.Rectangle((-(bl/2+fl+ifl),-itw/2),(-(bl/2+fl+itl),itw/2),layer=ml)

        plgs.add([inner_trace_fs,inner_trace_ns])

        outer_funnel_fs = funnel(itw,ftw,layer = ml,length = ofl)
        translate(outer_funnel_fs,(bl/2+fl+itl,0))

        outer_funnel_ns = funnel(itw,ftw,layer = ml,length = ofl)
        outer_funnel_ns.rotate(np.pi)
        translate(outer_funnel_ns,(-(bl/2+fl+itl),0))

        plgs.add([outer_funnel_fs, outer_funnel_ns])


        #add fine_coarse overlap!!
        #fco =

        junction_length = bl + 2*fl + 2*itl +2*ofl


    if params['junction_type'] == 'dolan_squid':
        pass

    if params['junction_type'] == 'bridge_free':
        pass

    if params['junction_type'] == 'bridge_free_squid':
        pass



    return plgs, junction_length, dolan_ps
    '''
    Returning dolan_ps never overwrites things we set in qubit_ps or the original
    parameters.  It's a way of making sure that later in the program we know
    What the default parameters in earlier steps were.  This can be abused,
    but it lets us use the apprpriate default parameters, rather than have to
    move around parameters which aren't applicable, i.e. having the default
    bridge_free parameters passed around when you're making dolan qubits.
    '''



