    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
        Module DS
As a gd is a “group of data” of determined (known and not too big) length,
a ds (Data Stream) is used to handle sampled data (in time domain) with unknown 
(or very big) length, of which one has at a given moment just a chunk.
        DS operation
The basic idea is the client-server metaphor: the client asks for chunks of data, defining at the beginning the modalities and then calling iteratively the server. 
  There are three fundamental operations:
•	Initial Setting
•	ds Servicing
•	Client Processing

Data stream server:
 - frame
 - HDF5
 - Simulation

'''
import numpy as np

Ntot=0
Ntot1=0

class ds:
    '''
    baselen   typical user max length (multiple of 4)
    dt        sampling time
    dsconf    ds configuration:
                start index   typical nny
                shift         typical ceil(baselen/2) or baselen
    '''
    def __init__(self,baselen,dt,dsconf):
    
        baselen4=np.ceil(baselen/4)
        if baselen4*4 > baselen:
            baselen=baselen4*4
            print('*** New baselen = ',baselen)
        ny=baselen*2
        baselen1=3*baselen4
        self.y=np.zeros(ny)
        self.ny=ny
        self.dt=dt
    


def ds_shift_charge():
    pass



def ds_filter():
    pass


def ds_simul():
    pass


def ds_read_frame():
    pass


def ds_hdf5():
    pass


def ds_bsd():
    pass


def ds_pows():
    pass


def ds_spectrogram():
    pass


def ds_stat():
    pass


def ds_quality():
    pass