    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

import ASTROTIME

# CW sources  ---------------------------------------

class cw():
    def __init__(self):
        self.name='crab';
        self.a=83.633218
        self.d=22.01446361111;
        self.v_a=0.0   # marcs/y
        self.v_d=0.0   # marcs/y
        self.pepoch=5.463200000022558e+04
        self.f0=59.473551262000001
        self.df0=-7.433700409800000e-10
        self.ddf0=2.359316462680000e-20
        self.fepoch=54936

        self.t00=ASTROTIME.v2mjd([2000,1,1,0,0,0])
        self.eta=-0.7667
        self.iota=62.165
        self.siota=0.8531
        self.psi=35.155
        self.spsi=0.0908
        self.h=1.400000000000000e-24
        self.snr=1
        self.coord=0
        self.chphase=0

        self.dfrsim=-0.2
        self.ephfile='/storage/targeted/band_recos/crab_ephemeris_file_O1.txt'
        self.ephstarttime=57249
        self.ephendtime=57403



# Antennas ------------------------- 

class Antenna():
    def __init__(self):
        self.name='virgo'
        self.lat=43.6314133
        self.long=10.5044968
        self.azim=199.4326
        self.incl=0
        self.height=4
        self.whour=+1
        self.shour=+2
        self.type=2


# 5-vect --------------------------------------




# Other ---------------------------------------


