# copy and paste helper snippets

# initial importing ------------------

import numpy as np
import SNAGPY,GD,GD2,BASIC,ML_PY,SERV,STAT,SIGNAL,ASTROTIME,GD2,GWDATA
import BSD,GUISNAG,GUISNAG_QT,GUISNAG_TK,GWDATA
import importlib,copy
# importlib.reload(Module)
from importlib import reload

# useful services -------------

from BASIC import tic,toc
tic()

from BASIC import var,Var

# global sh_edited_data

# symbols -----------------

snagpy_p,gwdata_p,examples_p,exper_p,pers_p,projects_p,\
antennas_p,cwinj_o2_p,cwinj_o3_p,cwsour_p,\
table_ligoh_p,table_ligol_p,table_virgo_p,table_kagra_p,de440s_p\
    =GWDATA.set_symbols()

table_par=[1261872018.0, 1577490618.0, 600.0]   # start, stop in gps, step in s
# Rectangular equatorial: posx,posy,posz in light seconds , velx/C vely/C velz/C deinstein
# Time is: TDT given as mjd and gpstime [s]. Leap seconds: 37 last added January 2017

virgo,ligol,ligoh,kagra=GWDATA.antennas(antennas_p)
