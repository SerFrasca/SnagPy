# Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
    
'''
      Module BSD
    BSD functions

A bsd is a standard gd with peculiar cont structure (primary bsd).
It is very useful for narrow band signal detection.

It is a type 1 gd that contains complex sampled data with information on the gravitational wave signal in a band. 
The bsd have a particular cont structure.
There are basic or “primary” bsd that are stored and collected in a particular folder structure, the BSD collection.
Starting from a BSD collection, one can construct other “secondary” bsd.

In secondary bsd should be an oper structure (or list) that contains information and parameters used to generate it. 

If the starting bsd was secondary, the oper structure contains inside the original oper structure.  And so on, recursively.

A full-band bsd is a particular type of bsd, that has the entire (Nyquist) band. In this case some parameters can be missing. 
It has the half sampling frequency of the real data.

A bsd is created starting from the GWA sampled data with Doppler information.
It is then "decorated" with some basic analysis.

This module allows general bsd management. (in Matlab are more than 200 functions)

The module BSD_NC is focalized to the non-coherent detection procedures.


Sections:
> File access                   -> dummy_file
> BSD basic                     -> dummy_basic
> BSD services                  -> dummy_service
> Create BSD                    -> dummy_create
> BSD decoration                -> dummy_deco
> BSD collection access         -> dummy_access
> BSD Lego                      -> dummy_lego
> BSD 5-vect                    -> dummy_5vec
> BSD PSS                       -> dummy_pss
> BSD other searches            -> dummy_other
> BSD time events               -> dummy_time
> BSD GUI                       -> dummy_gui
> BSD analysis procedures       -> dummy_analys
> Non-BSD useful functions      -> dummy_usef


'''

def sections():
    sec=[
    ]
    return sec

import GD,ML_PY

# File access ----------------------------

def dummy_file():
    '''
    The bsds can be stored in files as:
    - HDF5
    - mat v.7
    - mat v7.3
    - pickle
    '''

def bsd_lm7(fil):
    g=ML_PY.gd_lm7(fil)
    cont=g.cont
    names=cont.dtype.names
    cval=cont[0][0]
    nval=len(names)
    dt=[]
    for ii in range(nval):
        aa=cval[ii]
        dt.append(len(aa.dtype))

    return names,cval,dt




# BSD basic ----------------------------

def dummy_basic():
    '''
    Basic bsd operation
    '''

def bsd_type(bin):
    '''
    type and format of a bsd (or gd)
    '''


def bsd_real(bin):
    '''
    Converts to real data (sampling rate is doubled)
    '''

def bsd_freqshift(bin,dfr,tref=57279):
    '''
    Frequency shift for bsd. Shift is circular.
    
    in    input bsd
    dfr   frequency shift
    tref  reference time (def 14 sep 2015)

    '''


def enlarge(bin,enl=2):
    '''
    Enlarge frequency domain

    
    in    input bsd
    enl   frequency domain enlargement

    '''





# BSD services ----------------------------

def dummy_service():
    '''
    Here are the basic services of a bsd. 
    They are peculiar classes:
    > BSD_basic
    > BSD_doppler
    > BSD_subsp
    > BSD_clean
    > BSD_interv
    > BSD_peaks
    > BSD_sky
    '''

class BSD_basic:
    '''
    BSD basic parameters.

    The attributes are:

    > t0
    > inifr
    > bandw
    > ant
    > run
    > cal
    > Tfft
    > tcreation
    > durcreation

    '''

    def __init__(self,bas_dict):
        self.t0=bas_dict['t0']
        self.inifr=bas_dict['inifr']
        self.bandw=bas_dict['bandw']
        self.ant=bas_dict['ant']
        self.run=bas_dict['run']
        self.cal=bas_dict['cal']
        self.Tfft=bas_dict['Tfft']
        self.tcreation=bas_dict['tcreation']
        self.durcreation=bas_dict['durcreation']


class BSD_doppler:
    '''
    BSD Doppler data.

    The attributes are:

    > v_eq is the detector velocity components in unit of c, with the time in s
      since the beginning of data
    > p_eq is the detector position components in a SSB system in unit light seconds,
      with the time in s since the beginning of data

    '''

    def __init__(self,v_eq,p_eq):
        self.v_eq=v_eq
        self.p_eq=p_eq


class BSD_subsp:
    '''
    BSD short spectrum.

    The attributes are:

    > 
    '''


# Create BSD ----------------------------

def dummy_create():
    '''
    '''


# BSD decoration ----------------------------

def dummy_deco():
    '''
    '''


class deco_tfclean:
    '''
    Cleaning decoration class
    
     - tf filter
     - persistence
     - noise adaptive filter

       tfstr         tf structure as created by bsd_peakmap
       anabasic      procedure control structure
        .pers_thr  persistence threshold (0 -> no threshold; typical 0.3)
        .tfh_df    time-frequency histogram df
        .tfh_dt    time-frequency histogram dt (days) def 0.5
        .tfh_pht   time-frequency histogram time phase (hours) def 8 local time
        .tfh_thr   time_frequency threshold (0 -> no threshold)
        .noplot    = 0  plot (default)
    '''


    class deco_shortsp:
        '''
        '''


    class deco_sky:
        '''
        '''


# BSD collection access ----------------------------

def dummy_access():
    '''
    '''



# BSD Lego ----------------------------

def dummy_lego():
    '''
    '''



# BSD 5-vect ----------------------------

def dummy_5vec():
    '''
    '''



# BSD PSS ----------------------------

def dummy_pss():
    '''
    '''



# BSD other searches ----------------------------

def dummy_other():
    '''
    '''



# BSD time events ----------------------------

def dummy_time():
    '''
    '''







# BSD GUI ----------------------------

def dummy_gui():
    '''
    '''







# BSD analysis procedures ----------------------------

def dummy_analys():
    '''
    '''







# Non-BSD useful functions ----------------------------

def dummy_usef():
    '''
    '''



