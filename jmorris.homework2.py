SoftwareDesign
==============

The base repository for Olin College's Software Design Spring 2014\

SoftDes
Julian Morris
Homework 2 


"""HW 2 Exercise 2.4"""

"""part 1"""
"""The volume of a sphere with radius r is 4/3 π r3.
What is the volume of a sphere with radius 5?"""

import math

r=5.0
v=4.0/3.0*math.pi*r**3.0

print v

"""part 2"""
"""Suppose the cover price of a book is $24.95, 
but bookstores get a 40% discount. 
Shipping costs $3 for the first copy and 75 cents for each additional copy. 
What is the total wholesale cost for 60 copies?"""

c=24.95*0.6
e=c+3+59*(c+0.75)

print e

"""part 3"""
"""If I leave my house at 6:52 am and 
run 1 mile at an easy pace (8:15 per mile), 
then 3 miles at tempo (7:12 per mile) and 1 mile at easy pace again,
what time do I get home for breakfast?"""

"""
using unum
"""

from unum.units import *
runtime=(8*min+15*s)+3*(7*min+12*s)+(8*min+15*s)
minutes=(52*min)+runtime
total=((6*h)+minutes)

print total

"""
without using unum
"""

minutes=((8*2+7*3)+(52)+(15.0/60*2)+(12.0/60*3))
finish_minutes=minutes%60
finish_hour=6+int(minutes/60)

print '%d:%dam' %(finish_hour, finish_minutes)


"""
Exercise 3.5 part 1
"""

def twice(f):
    f()
    f()
    
def four_times(f):
    f()
    f()
    f()
    f()
    
def pm():
    print '+ - - - -' ,
    
def plusminus():
    twice(pm)
    print '+'

def v():
    print '|        ', 

def vertical():
    twice(v)
    print '|'
    
def sets():
    plusminus()
    four_times(vertical)

def result():
    twice(sets)
    plusminus()

result()


"""
Exercise 3.5 part 2
"""

def twice2(f):
    f()
    f()
    
def four_times2(f):
    f()
    f()
    f()
    f()
    
def pm2():
    print '+ - - - -' ,
    
def plusminus2():
    four_times2(pm2)
    print '+'

def v2():
    print '|        ', 

def vertical2():
    four_times2(v2)
    print '|'

def section2():
    plusminus2()
    four_times2(vertical2)

def result2():
    four_times2(section2)
    plusminus2()
    
result2()

"""
Exercise 5.3 part 1
"""


def check_fermat(a,b,c,n):
    if n>2:
        if a**n+b**n==c**n:
            print "Holy smokes, Fermat was wrong!"
        else:
            print "No, that doesn't work."
        
check_femat(2,3,4,5)

"""
Exercise 5.3 part 2
"""
    
def check_input():
    a = int(raw_input("a: "))
    b = int(raw_input("b: "))
    c = int(raw_input("c: "))
    n = int(raw_input("n: "))
    print check_fermat(a, b, c, n)
    

check_input()

"""
Exercise 6.1. Write a compare function that returns 1 if x>y,
 0 if x==y, and -1 if x<y.
"""

def compare(x,y):
    if x>y:
        return 1
    elif x==y:
        return 0
    elif x<y:
        return -1
    
print compare(2,2)
