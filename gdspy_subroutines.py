# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:44:42 2014

@author: jzb4
"""

import gdspy
#import matplotlib.pyplot as plt
import numpy as np
import math

um = 1000
mm = 1000000
mil = 25400


def flatten(l):
    '''Flatten lists of lists (arbitrary depth)'''
    out = []
    for item in l:
        if isinstance(item, (list, tuple, np.ndarray)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out
  
def add_plgs(cell,plg_list):
    plg_list = flatten(plg_list)
    for plg in plg_list:
        cell.add(plg)
        
def translate(plgs,delta):
    '''translate a plg array (dx,dy)'''
    plgs = flatten([plgs])
    for plg in plgs:
        d = list(delta)
        pts = plg.points
        for point in pts:
            point += d
        plg.points = pts
    
def xflip(plgs):
    '''Flip a plg array across the x axis'''
    plgs = flatten([plgs])
    for plg in plgs:
        pts = plg.points
        for point in pts:
            point[1] *= -1
        plg.points = pts   
   
def yflip(plgs):
    '''Flip a plg array across the y axis'''
    plgs = flatten([plgs])
    for plg in plgs:
        pts = plg.points
        for point in pts:
            point[0] *= -1
        plg.points = pts   
        
'''Funnel subroutines, taken from https://gist.github.com/Alquimista/1274149'''
def binomial(i, n):
    """Binomial coefficient"""
    return math.factorial(n) / float(
        math.factorial(i) * math.factorial(n - i))
 
def bernstein(t, i, n):
    """Bernstein polynom"""
    return binomial(i, n) * (t ** i) * ((1 - t) ** (n - i))
 
def bezier(t, points):
    """Calculate coordinate of a point in the bezier curve"""
    n = len(points) - 1
    x = y = 0
    for i, pos in enumerate(points):
        bern = bernstein(t, i, n)
        x += pos[0] * bern
        y += pos[1] * bern
    return x, y
 
def bezier_curve_range(n, points):
    """Range of points in a curve bezier"""
    for i in xrange(n):
        t = i / float(n - 1)
        yield bezier(t, points)
        
        
        
def funnel(start_width,end_width,layer=0, length = 500,steepness = 0.5,steps =50):
    sw = start_width
    ew = end_width
    
    steps = 50
    x = np.zeros(steps*2+1)
    y = np.zeros(steps*2+1)
    
    control_points = [(0,sw/2),
                      (steepness*length,sw/2),
                      ((1-steepness)*length,ew/2),
                      (length,ew/2)]
    
    #Draw the first curve
    for idx,pt in enumerate(bezier_curve_range(steps, control_points)):
        x[idx] = pt[0]
        y[idx] = pt[1]

    #Draw the second curve (backwards, to keep the points in the right order)
    for idx,pt in enumerate(bezier_curve_range(steps, control_points)):
        x[2*steps-idx-1] = pt[0]
        y[2*steps-idx-1] = -pt[1]
    
    #Finish the plg.
    x[-1] = x[0]
    y[-1] = y[0]

    pts = list(np.zeros(2*steps+1))
    for idx in range(steps*2+1):
        pts[idx] = (x[idx],y[idx])
        
    plg = gdspy.Polygon(pts,layer=layer)
    return plg
        
        
        
def round_to_poly(rnd):
    '''Convert GDSPY Round objects (useless) into GDSPY polygon objects'''
    layer = rnd.layers[0]
    plgs = rnd.polygons
    pts = plgs[0]
    outplg = gdspy.Polygon(pts,layer=layer)
    return outplg

        

    