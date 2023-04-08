    # Copyright (C) 2023  Sergio Frasca, Edoardo Giancarli
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
        Module ML_PY

This module deals with operations with Matlab. 
The most important part is the management of the mat files. 
Some mat files could not be read because contain Matlab tables: 
in such case the solution is saving the Matlab table 
as a csv file (by the writetable internal function
or table2cvs Snag function).

Sections:
> read - write v.7 format mat file      -> dummy_v7
> read - write v.7.3 format mat file    -> dummy_v73
> read - write csv                      -> dummy_csv
'''
import numpy as np
import scipy.io as sio
import csv
import mat73
import GD,SERV,BASIC

# read - write v.7 format mat file -------------

def dummy_v7():
    '''
    '''

def loadmat7(fil):
# Generic load
    dat=sio.loadmat(fil)

    return dat



def gd_lm7(fil):
    mat=sio.loadmat(fil)
    ll=list(mat.keys())
    nam=ll[3]
    arstr=mat[nam]
    dic=BASIC.arrstruct2dict(arstr)
    x_=dic['x']
    y=dic['y']
    y=np.squeeze(y)
    n_=dic['n']
    typ=dic['type']
    ini_=dic['ini']
    dx_=dic['dx']
    capt_=dic['capt']
    cont_=dic['cont']
    unc_=dic['unc']
    uncx_=dic['uncx']
    
    if typ == 1:
        gdout=GD.gd(y,ini=ini_,dx=dx_,typ=1,capt=capt_,cont=cont_,unc=unc_,uncx=uncx_)
    else:
        gdout=GD.gd(y,x=x_,ini=ini_,dx=dx_,typ=2,capt=capt_,cont=cont_,unc=unc_,uncx=uncx_)

    return gdout


def gdsavedicmat7(ingd):
    '''
    save a gd transformed to a dictionary
    '''
    nam=SERV.retrieve_name(ingd)
    outdic=GD.gd2dict(ingd)
    filnam=nam+'_dic.mat'
    sio.savemat(filnam,outdic)

    return 'File '+filnam+' created'


def gdloaddicmat7(fil):
# load a gd transformed to a dictionary
    indic=sio.loadmat(fil)
    outgd=GD.dict2gd(indic)

    return outgd


# read - write v.7.3 format mat file ---------------------

def dummy_v73():
    '''
    '''

def loadmat73(fil):
    '''
    Generic load
    '''
    dat=mat73.loadmat(fil)
    obj=list(dat.keys())

    return dat


def gd_lm73(fil):
    '''
    load a file containing a gd
    '''
    dat=mat73.loadmat(fil)
    obgd=list(dat.keys())
    obgd=obgd[0]
    gddic=dat[obgd]
    y=gddic['y']
    ini_=gddic['ini']
    dx_=gddic['dx']
    typ=gddic['type']
    if typ == 2:
        x_=gddic['x']
    capt_=gddic['capt']
    cont_=gddic['cont']
    unc_=gddic['unc']
    uncx_=gddic['uncx']

    if typ == 1:
        gdout=GD.gd(y,ini=ini_,dx=dx_,typ=1,capt=capt_,cont=cont_,unc=unc_,uncx=uncx_)
    else:
        gdout=GD.gd(y,x=x_,ini=ini_,dx=dx_,typ=2,capt=capt_,cont=cont_,unc=unc_,uncx=uncx_)
    
    return gdout



# read - write csv --------------------------

def dummy_csv():
    '''
    '''

def csv_read(filn):
    with open(filn, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            print(', '.join(row))


# BSD --------------------------------
