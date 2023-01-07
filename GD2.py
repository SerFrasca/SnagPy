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
# > y      the ordinate (basic data)
# > n      the length
# > ini    initial abscissa (used in type 1 gds)
# > dx     sampling step (used in type 1 gds)
# > x      abscissas (used in type 2 gds)
# > typ    determine type 1 (virtual abscissa) or type 2 (real abscissa)
# > capt   caption (a string)
# > cont   a control variable (in the case of a bsd, a special structure)
# > unc    ordinate uncertainty (typically not used)
# > uncx   abscissa uncertainty (used for other)
#

Any gd can be modified by edit_gd. 
A gd can be created giving just the number of samples, 
or transform a 1-dimensional array to a gd, or give the other information. 
Addition and multiplication are overloaded for gd objects.
"""

class gd2:    # gd2 creation
    def __init__(self,y,**gdpar): # y ordinate or n
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
        self.label='SnagPy gd created on '+time.asctime()

    def __add__(self,other):
        outgd=copy.copy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd.y=outgd.y+other
        else:
            outgd.y=outgd.y+other.y
        return outgd   

    def __radd__(self,other):
        outgd=copy.copy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd.y=outgd.y+other
        return outgd   

    def __mul__(self,other):
        outgd=copy.copy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd.y=outgd.y*other
        else:
            outgd.y=outgd.y
            other.y
        return outgd   

    def __rmul__(self,other):
        outgd=copy.copy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd.y=outgd.y*other
        return outgd   
  

def edit_gd(ingd,**gdpar): # 'new'-1  -> new object
    if 'new' in gdpar:
        outgd=copy.copy(ingd)
    else:
        outgd=ingd
    if 'x' in gdpar:
        outgd.x=gdpar['x']
        outgd.typ=2
    if 'y' in gdpar:
        outgd.y=gdpar['y']
        outgd.n=len(outgd.y)
    if 'ini' in gdpar:
        outgd.ini=gdpar['ini']
    if 'dx' in gdpar:
        print(gdpar)
        outgd.dx=gdpar["dx"]
    if 'capt' in gdpar:
        outgd.capt=gdpar['capt']
    if 'cont' in gdpar:
        outgd.x=gdpar['cont']
    if 'unc' in gdpar:
        outgd.unc=gdpar['unc']
    if 'uncx' in gdpar:
        outgd.uncx=gdpar['uncx']
    return outgd


def x_gd(ingd):
    if ingd.typ == 1 :
        x=ingd.ini+np.arange(ingd.n)*ingd.dx
    else:
        x=ingd.x
    x=SERV.deshape(x)

    return x        
