# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:44:42 2014

@author: jzb4
"""

import gdspy
import matplotlib.pyplot as plt
import numpy as np
import math

um = 1000
mm = 1000000
mil = 25400

#==============================================================================
# Operations on lists, cells
#==============================================================================
def flatten(l):
    '''Flatten lists of lists (arbitrary depth)'''
    out = []
    for item in l:
        if isinstance(item, (list, tuple, np.ndarray)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out


class poly_list:
    '''a poly_list is a dictionary, indexed by layer, that contains a gdspy
    polygonSet object for each layer.  However, adding polygons
    or other polygon_lists doesn't require such considerations.
    
    This is the polygon object of choice, and subroutines will sometimes force 
    you into using them.'''

    def __init__(self,item_list=None):
        self.dict = {}
        if item_list:
            self.add(item_list)

    def __getitem__(self,idx):
        '''You can index the polylist by layer'''
        return self.dict[idx]
#
    def __setitem__(self,idx,item):
        self.dict[idx] = item

    def __iter__(self):
        for x in self.dict:
            yield x

    def has_key(self,x):
        return self.dict.has_key(x)

    def extend(self,item):
        self.add(item)

    def add(self,items):
        '''Adding polygons, other poly_list objects, or lists of those items
        to a poly_list'''
        items = flatten([items,])
        for item in items:
            if isinstance(item,poly_list):
                for layer in item.dict:
    
                    if self.dict.has_key(layer):
                        self.dict[layer].polygons.extend(item[layer].polygons)
                    else:
                        self.dict[layer] = item[layer]
    
            else:
                layer = item.layer
                points = item.points
                if self.dict.has_key(layer):
                    self.dict[layer].polygons.append(points)
                else:
                    self.dict[layer] = gdspy.PolygonSet([points,],layer=layer)
    pass



def add_plgs(cell,plg_list):
    #Layers and Datatypes need to be as long as the number of plgs

    if isinstance(plg_list, poly_list):
        for layer in plg_list.dict:
            plg_list[layer].layers = [layer,]*len(plg_list[layer].polygons)
            plg_list[layer].datatypes = [0,]*len(plg_list[layer].polygons)
            cell.add(plg_list[layer])

#==============================================================================
# Polygon transformations and operations
#==============================================================================
'''There are times when raw polygons are better than poly_lists, so the basic
transformations should work on both, I think.'''

def translate(plgs,delta):
    '''translate a plg array or polylist by (dx,dy)'''
    plgs = flatten([plgs])
    for plg in plgs:
        if isinstance(plg,poly_list):
            for layer in plg.dict:
                for plgn in plg[layer].polygons:
                    for point in plgn:
                        point += delta
        else:           
            pts = plg.points
            for point in pts:
                point += delta
            plg.points = pts


def xflip(plgs):
    '''Flip a poly_list, plg, or combination across the x axis'''
    plgs = flatten([plgs])
    for plg in plgs:
        if isinstance(plg,poly_list):
            for layer in plg.dict:
                for plgn in plg[layer].polygons:
                    for point in plgn:
                        point[1] *= -1
        else:           
            pts = plg.points
            for point in pts:
                point[0] *= -1
            plg.points = pts

def yflip(plgs):
    '''Flip a poly_list, plg, or combination across the y axis'''
    plgs = flatten([plgs])
    for plg in plgs:
        if isinstance(plg,poly_list):
            for layer in plg.dict:
                for plgn in plg[layer].polygons:
                    for point in plgn:
                        point[0] *= -1
        else:           
            pts = plg.points
            for point in pts:
                point[0] *= -1
            plg.points = pts


def bounding_box(plgs):
    '''Find the extend of a poly_list'''
    if not isinstance(plgs,poly_list):
        plgs = poly_list(plgs)
    
    min_x = 1e10
    max_x = -1e10
    
    min_y = 1e10
    max_y = -1e10
    
    for layer in plgs:
        for plg_pts in plgs[layer].polygons:
            for pt in plg_pts:
                if pt[0] > max_x: max_x = pt[0]
                if pt[0] < min_x: min_x = pt[0]
                if pt[1] > max_y: max_y = pt[1]
                if pt[1] < min_y: min_y = pt[1]
        
    return (min_x,min_y),(max_x,max_y)

#==============================================================================
# Bezier code
#==============================================================================
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


#==============================================================================
# Shapes
#==============================================================================

def round_plg(*args,**kwargs):
    '''This is a wrapper for the gdspy.Round object, which is basically useless.'''
    rnd =  gdspy.Round(*args,**kwargs)
    layer = rnd.layers[0]
    plgs = rnd.polygons
    pts = plgs[0]
    outplg = poly_list(gdspy.Polygon(pts,layer=layer))
    return outplg
    
def text_plg(*args,**kwargs):
    '''This is a wrapper for the gdspy.Text object, which is basically useless.'''
    txt =  gdspy.Text(*args,**kwargs)
    layer = txt.layers[0]
    plg_pts = txt.polygons
    outplgs = poly_list()
    for pts in plg_pts:
        new_plg = gdspy.Polygon(pts,layer=layer)
        outplgs.add(new_plg)
    return outplgs
    
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


def wafer(diameter,flat_length,layer,num_points=100):
    angle = np.tan(flat_length/diameter)
    radius = diameter/2
    pts = np.empty([num_points,2])

    angles = np.linspace(angle,2*np.pi-angle,num_points-1)-np.pi/2
    for i,angle in enumerate(angles):
        pts[i,0] = radius*np.cos(angle)
        pts[i,1] = radius*np.sin(angle)

    pts[num_points-1,0] = pts[0,0]
    pts[num_points-1,1] = pts[0,1]

    plg = gdspy.Polygon(pts,layer=layer)
    return plg

def two_inch_wafer(layer=15):
    num_points = 100
    diameter = 50.8 * mm
    flat_length = 16.0 * mm
    outer_plg = wafer(diameter,flat_length,layer,num_points)

    inner_plg = wafer(diameter-4*mm,(diameter-4*mm)/diameter*flat_length,
                      layer,num_points)

    return poly_list([inner_plg,outer_plg])



#==============================================================================
# Clipping
#==============================================================================
'''Try to use boolean operations in an efficient way, this isn't very fast code'''

def plg_bool(plgsa,plgsb,operation,**kwargs):
    '''take two poly_lists and perform layer-wise boolean operations.

    This is not grade-A code.  It doesn't support subtracting a contained object.'''
    if not isinstance(plgsa,poly_list) or not isinstance(plgsb,poly_list):
        print 'both objects need to be poly lists!'
        pass

    if operation == 'int':
        f = lambda a,b: a and b
    elif operation == 'union':
        f = lambda a,b: a or b
    elif operation == 'sub':
        f = lambda a, b: a and not b

    out_plgs = poly_list()

    for layer in plgsa:
        if plgsb.has_key(layer):
#            print 'found two in the same layer'
#            print 'a:', plgsa[layer]
#            print 'b:', plgsb[layer]
            ret = gdspy.boolean([plgsa[layer], plgsb[layer]], f, eps=1e-10)
#            print 'ret:', ret
            ret.layers = [layer,]*len(ret.polygons)
            ret.datatypes = [0,]*len(ret.polygons)
#            print 'ret:', ret
            out_plgs[layer] = ret


    return out_plgs

def populate_wafer(dice_x, dice_y, radius = 25.4*mm - 2.0*mm):
    '''
    1) See how many columns can fit.
    
    2) For each # of columns:
        2a) fill the column top to bottom with devices.
        2b) go through and see which devices have (at least) one corner outside
            the bounding circle.
            remove these.
            (later I want to delete a row if it only has one device)
        2c) count the remaining devices.
        
    3) use the number of columns that fit the most.
    
    The return format is the translation indices from the center of the gds
    to the center of the dice.  
    '''
    max_columns = int((2*radius)//dice_x)
    max_rows = int((2*radius)//dice_y)

    idcs = []
    num_devs = []
    
    for num_cols in range(1,max_columns+1):
        valid_idcs =[]
        num_devices = num_cols*max_rows
        for col in range(num_cols):
            col -= num_cols/2.0 -0.5
            left_x = (col -0.5)*dice_x
            right_x = (col + 0.5)*dice_x
            
            for row in range(max_rows):
                row -= max_rows/2.0 - 0.5
                valid = True
                bottom_y = (row - 0.5)*dice_y
                top_y = (row + 0.5)*dice_y
                
                pts = [(left_x,bottom_y),
                       (left_x,top_y),
                        (right_x,bottom_y),
                        (right_x,top_y)]
                
                for pt in pts:
                    ext = pt[0]**2+pt[1]**2
                    if ext > radius**2:
                        valid = False
                        num_devices -= 1
                        break
                    
                if valid:
                    valid_idcs.append((col,row))
                    
                
        idcs.append(valid_idcs)
        num_devs.append(num_devices)
        
    opt_idx = np.argmax(num_devs)
    return idcs[opt_idx]