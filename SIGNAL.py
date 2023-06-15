    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
        Module SIGNAL

Signal analysis

Sections:

> Signals                   -> dummy_sig
> Operations                -> dummy_oper
> Filters                   -> dummy_filt
> Particular transform      -> dummy_transf
> Processes                 -> dummy_proc

'''

def sections():
    sec=[
    ]
    return sec

import scipy.signal as sip
import scipy.interpolate as sinp
import numpy as np
import copy

# Signals ---------------------

def dummy_sig():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def chirp():
    pass


def gausspulse():
    pass


def sawtooth():
    pass


def square():
    pass


# Operations --------------------

def dummy_oper():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def detrend(ingd):
    pass

def resample(ingd,newdx):
    pass

def relminmax(ingd):
    pass

def cubicspline():
    pass

def convol():
    pass

def correl():
    pass


# Filters ----------------------

def dummy_filt():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def Filter(ingd,a,b,zi):  # not to be confused with "filter" function
    pass

def FiltFilt(ingd,a,b,norm=1):
    '''
    Apply a digital filter forward and backward to a signal

    ingd    gd or 1-D array
    a       denominator coefficient vector of the filter (AR or IIR coefficients)
    b       denominator coefficient vector of the filter (MA or FIR coefficients)
    norm    normalization type (0 -> no, 1 -> max DC)
    '''
    
    if isinstance(ingd,np.ndarray):
        icgd=0
        y=ingd
    else: 
        icgd=1
        y=ingd.y

    out=sip.filtfilt(b,a,y)
    if norm == 1:
        outdc=sip.filtfilt(b,a,np.ones(100000))
        oo=max(outdc)
        out=out/oo

    if icgd == 1:
        icgd=copy.deepcopy(icgd)
        icgd.y=out
        out=icgd

    return out


def multi_sour_adapt():
    '''
    '''


def bode():
    pass


# Particular transform

def dummy_transf():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def hilbert(ingd):
    pass

def radon(ingd):
    pass

def hough(ingd):
    pass


    

# Processes --------------------

def dummy_proc():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def normal_proc():
    pass

def poisson_proc():
    pass
