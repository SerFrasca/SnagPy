    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

import GD,ML_PY

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



def mult_gd_par():
    pass
