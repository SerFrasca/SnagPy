    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

import numpy as np
from scipy.fft import fft, ifft
import cmath as cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
import time
import SERV

pi=cm.pi
deg2rad=pi/180

# class gd2  -----------------------------------------------
"""  
gd2 is the class-container for 2-D data in SnagPy.
The attributes are:
# > y      the 2-dim ordinate (basic data)
# > n      
# > m
# > ini    initial abscissa (used in type 1 gds)
# > dx     sampling step (used in type 1 gds)
# > ini2
# > dx2
# > x      abscissas (used in type 2 gds)
# > typ    determine type 1 (virtual abscissa) or type 2 (real abscissa)
# > capt   caption (a string)
# > cont   a control variable (in the case of a bsd, a special structure)
# > unc    ordinate uncertainty (typically not used)
# > uncx   abscissa uncertainty (used for other)
#

Any gd2 can be modified by edit_gd2. 
A gd2 can be created giving just the m and n parameters, 
or transform a 2-dimensional array to a gd2, or give the other information. 
Addition and multiplication are overloaded for gd2 objects.
"""

class gd2:    # gd2 creation
    def __init__(self,y,**gdpar): # y ordinate or n and m
        if isinstance(y,int):
            y=np.zeros(y)
        self.y=y
        if 'ini' in gdpar:
            self.ini=gdpar['ini']
        else:
            self.ini=0
        if 'dx' in gdpar:
            self.dx=gdpar['dx']
        else:
            self.dx=1
        if 'x' in gdpar:
            self.x=gdpar['x']
            self.typ=2
        else:
            self.x=[]
            self.typ=1
        if 'y' in gdpar:
            self.y=y
        self.n=len(y) 
        if 'capt' in gdpar:
            self.capt=gdpar['capt']
        else:
           self.capt='' 
        if 'cont' in gdpar:
            self.cont=gdpar['cont']
        else:
            self.cont=0 
        if 'unc' in gdpar:
            self.unc=gdpar['unc']
        else:
            self.unc=0
        if 'uncx' in gdpar:
            self.uncx=gdpar['uncx']
        else:
            self.uncx=0
        self.label='SnagPy gd2 created on '+time.asctime()

    def __add__(self,other):
        outgd2=copy.deepcopy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd2.y=outgd2.y+other
        else:
            outgd2.y=outgd2.y+other.y
        return outgd2   

    def __radd__(self,other):
        outgd2=copy.deepcopy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd2.y=outgd2.y+other
        return outgd2   

    def __mul__(self,other):
        outgd2=copy.deepcopy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd2.y=outgd2.y*other
        else:
            outgd2.y=outgd2.y
            other.y
        return outgd2  

    def __rmul__(self,other):
        outgd2=copy.deepcopy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd2.y=outgd2.y*other
        return outgd2   
  

def edit_gd2(ingd2,**gdpar): # 'new'-1  -> new object
    if 'new' in gdpar:
        outgd2=copy.deepcopy(ingd2)
    else:
        outgd2=ingd2
    if 'x' in gdpar:
        outgd2.x=gdpar['x']
        outgd2.typ=2
    if 'y' in gdpar:
        outgd2.y=gdpar['y']
        outgd2.n=len(outgd2.y)
    if 'ini' in gdpar:
        outgd2.ini=gdpar['ini']
    if 'dx' in gdpar:
        print(gdpar)
        outgd2.dx=gdpar["dx"]
    if 'capt' in gdpar:
        outgd2.capt=gdpar['capt']
    if 'cont' in gdpar:
        outgd2.x=gdpar['cont']
    if 'unc' in gdpar:
        outgd2.unc=gdpar['unc']
    if 'uncx' in gdpar:
        outgd2.uncx=gdpar['uncx']
    return outgd2

