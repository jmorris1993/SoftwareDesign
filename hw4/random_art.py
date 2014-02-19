# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 11:34:57 2014

@author: pruvolo
"""

# you do not have to use these particular modules, but they may help
from random import randint
import Image

def build_random_function(min_depth, max_depth):
    if max_depth<=0:
        xylist=[["x"],["y"]]
        return xylist[randint(0,1)]
    
    a=build_random_function(min_depth-1, max_depth-1)
    b=build_random_function(min_depth-1, max_depth-1)
    funclist=[["prod",a,b],["sin_pi",a],["cos_pi",a],["x"],["y"]]
        
    if min_depth>=0:
        return funclist[randint(0,2)]
    else:
        return funclist[randint(0,4)]
    
    
    
    
#print build_random_function(1,5)

from math import *
def evaluate_random_function(f, x, y):
    
    
    if f[0]=='x':
            return x
            
    elif f[0]=='y':
            return y
    
    elif f[0]=='sin_pi':
            return sin(pi*evaluate_random_function(f[1],x,y))
            
    elif f[0]=='cos_pi':
            return cos(pi*evaluate_random_function(f[1],x,y))
            
    elif f[0]=='prod':
            return evaluate_random_function(f[1],x,y )* evaluate_random_function(f[2],x,y)



#print evaluate_random_function(build_random_function(1,5),3,6)



    
def remap_interval(val, input_interval_start, input_interval_end, output_interval_start, output_interval_end):
    """ Maps the input value that is in the interval [input_interval_start, input_interval_end]
        to the output interval [output_interval_start, output_interval_end].  The mapping
        is an affine one (i.e. output = input*c + b).
    
        TODO: please fill out the rest of this docstring
    """
    
    scaling=(output_interval_end-output_interval_start)/float(input_interval_end-input_interval_start)
    finalval=(val-input_interval_start)*scaling+output_interval_start
    return finalval
    
#remap_interval(0,-1,1,0,255)


import time
    
def finalimage(IMAGESIZE,mindepth,maxdepth):   
    redfunc= build_random_function(mindepth,maxdepth)
    bluefunc= build_random_function(mindepth,maxdepth)
    greenfunc= build_random_function(mindepth,maxdepth)
    im = Image.new("RGB",(IMAGESIZE,IMAGESIZE))
 
    g=evaluate_random_function(greenfunc,200,108)
    green=remap_interval(g,-1,1,0,255)
    
    
    for x in range(0,IMAGESIZE):
        for y in range(0,IMAGESIZE):
            a=remap_interval(x,0,IMAGESIZE-1,-1,1)
            b=remap_interval(y,0,IMAGESIZE-1,-1,1)
            r=evaluate_random_function(redfunc,a,b)
            g=evaluate_random_function(greenfunc,a,b)
            b=evaluate_random_function(bluefunc,a,b)
            red=remap_interval(r,-1,1,0,255)
            green=remap_interval(g,-1,1,0,255)
            blue=remap_interval(b,-1,1,0,255)
            im.putpixel((x,y),(int(red),int(green),int(blue)))
    im.save("./image"+str(time.time())+".bmp")
        
finalimage(350,1,6)







