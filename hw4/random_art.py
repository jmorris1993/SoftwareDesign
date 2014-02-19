# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 11:34:57 2014

@author: pruvolo
"""

# you do not have to use these particular modules, but they may help
from random import randint
import Image

def build_random_function(min_depth, max_depth):
    """ builds the random function using recursion. The if and else statements are
    for making sure that the depth is greater than min_depth and less than max_depth)
    """
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
    """evaluates the random function created above, using recursion. The if/elif statments 
    are for evaluating different operations.
    """
    
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
    
        The mapping includes multiplying the interval by a scaler, and then shifting the interval.
        The scaling variable was calculated by taking the ratio of the magnitude of the two intervals.
        The final value is calculated using the scaler and then shifting it by the value of the output_interval_start.
    """
    
    scaling=(output_interval_end-output_interval_start)/float(input_interval_end-input_interval_start)
    finalval=(val-input_interval_start)*scaling+output_interval_start
    return finalval
    
#remap_interval(0,-1,1,0,255)


import time
    
def finalimage(IMAGESIZE,mindepth,maxdepth): 
    """
    This function takes in the image size, and the min and max depth of the function you want to generate. 
    The 4 variables at the top builds the 3 functions for each color, and the blank image.
    Then, the for loop is used for x,y coordinates of each pixel. a and b remaps the interval of 0 to 349 
    to -1 to 1. Then those values are used to evaluate the random functions generated at the top for each color.
    And finally, the generated values for the r,g,b are mapped to fit between 0 and 255.
    These red,green and blue values are put onto the pixel, and this will repeat for all pixels.
    Then, the image is saved in the current folder for us to enjoy!
    """
    redfunc= build_random_function(mindepth,maxdepth)
    bluefunc= build_random_function(mindepth,maxdepth)
    greenfunc= build_random_function(mindepth,maxdepth)
    im = Image.new("RGB",(IMAGESIZE,IMAGESIZE))
    
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







