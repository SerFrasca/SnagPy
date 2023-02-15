    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

# Data stream server:
#  - frame
#  - HDF5
#  - Simulation

import numpy as np

Ntot=0
Ntot1=0

class ds:
    def __init__(self,baselen,dt,dsconf):
#  baselen   typical user max length (multiple of 4)
#  dt        sampling time
#  dsconf    ds configuration:
#             start index   typical nny
#             shift         typical ceil(baselen/2) or baselen
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