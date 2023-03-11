# Cut and paste only the piece of your interest

# Inizialization

import numpy as np
import GD,GD2,BASIC,ML_PY,SERV,STAT,SIGNAL,ASTROTIME,GD2,GWDATA
import BSD,GUISNAG,GWDATA
import importlib
# importlib.reload(Module)

# useful services -------------

from BASIC import tic,toc
tic()

from BASIC import var

#-----------------------------------------------

# charge gd from a mat file v7

g=ML_PY.gd_lm7('L_C02_20170104_0060_0070_O2_tfstr.mat')
sp,tt,ff=STAT.gd_spectrogram(g,1000)
Ini,Fin,dim,Inan=BASIC.findnodata(g)

# Not used

def explore_hdf5(fil):
# explores hdf5 frame data
    gwData=h5py.File(fil, 'r')
    Groups=list(gwData.keys())
    print('Groups ',Groups)
    for it in Groups:
        dset=gwData[it]
        ldset=list(dset)
        print('Group ',it,': ',ldset)
        for it2 in ldset:
            print('Dataset ',it2,': ',it2 )
