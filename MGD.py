# Copyright (C) 2023  Sergio Frasca
#  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
        Module MGD

Multi-ordinate GD
'''

from scipy.fft import fft, ifft
import cmath as cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
import time
import numpy as np
import BASIC,SERV

pi=cm.pi
deg2rad=pi/180

# class gd  -----------------------------------------------

class mgd: 
    '''
            mgd creation

    mgd is a single abscissa, multiple ordinate class-container in SnagPy.

    The attributes are:

    > y      the ordinates (nxm array)
    > n      the length (number of rows)
    > m      number of channels (columns)
    > ini    initial abscissa (used in type 1 gds)
    > dx     sampling step (used in type 1 gds)
    > x      abscissas (used in type 2 gds)
    > typ    determine type 1 (virtual abscissa) or type 2 (real abscissa)
    > capt   caption (a string)
    > cont   a control variable (in the case of a bsd, a special structure)

    '''
    def __init__(self,y,**gdpar): # y matrix or (n,m) tuple
        if isinstance(y,tuple):
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
        if 'titl' in gdpar:
            self.titl=gdpar['titl']
        else:
           self.titl=tuple(np.repeat('')) 
        if 'scal' in gdpar:
            self.scal=gdpar['scal']
        else:
           self.scal=tuple(np.repeat(1))
        if 'capt' in gdpar:
            self.capt=gdpar['capt']
        else:
           self.capt='' 
        if 'cont' in gdpar:
            self.cont=gdpar['cont']
        else:
            self.cont=0 
        self.label='SnagPy mgd created on '+time.asctime()
