# copy and paste helper snippets

# initial importing ------------------

import numpy as np
import GD,BASIC,ML_PY,SERV,STAT,SIGNAL,ASTROTIME,GD2,GWDATA,BSD,GUISNAG
import importlib
# importlib.reload(Module)


# useful services -------------

from BASIC import tic,toc
tic()


# symbols -----------------

snagpy_p,gwdata_p,examples_p,exper_p,pers_p,projects_p,\
antennas_p,cwinj_o2_p,cwinj_o3_p,cwsour_p,\
table_ligoh_p,table_ligol_p,table_virgo_p,table_kagra_p,de440s_p\
    =GWDATA.set_symbols()

table_par=[1261872018.0, 1577490618.0, 600.0]

virgo,ligol,ligoh,kagra=GWDATA.antennas(antennas_p)