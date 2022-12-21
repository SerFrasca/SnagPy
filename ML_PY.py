import numpy as np
import scipy.io as sio
import h5py
import csv
import mat73
import GD,SERV,BASIC

'''
This module deals with operations with Matlab. 
The most important part is the management of the mat files. 
Some mat files could not be read because contain Matlab tables: 
in such case the solution is saving the Matlab table 
as a csv file (by the writetable internal function
or table2cvs Snag function).
'''
# read - write v.7 format mat file -------------

def loadmat7(fil):
# Generic load
    dat=sio.loadmat(fil)

    return dat


def gdloadmat7(fil):
# load a file containing a gd
    matc=sio.loadmat(fil)
    ke=list(matc.keys())
    gdnam=ke[3]
    mc=matc[gdnam]
    base=mc.base
    tup=mc.item()
    x_=tup[0]
 #   if len(x_) == 1:
 #       x_=[]
    x_=np.reshape(x_,len(x_))
    y=tup[1]
    y=np.reshape(y,len(y))
    if isinstance(y[1],complex):
        print('y is complex')
    else:
        print('y is real') 
    n_=tup[2]
    n_=np.reshape(n_,1)
    typ=tup[3]
    ini_=tup[4]
    ini_=np.reshape(ini_,1)
    dx_=tup[5]
    dx_=np.reshape(dx_,1)
    capt_=tup[6][0]
    cont_=tup[7]
    unc_=tup[8]
    uncx_=tup[9]
    
    if typ == 1:
        gdout=GD.gd(y,ini=ini_,dx=dx_,typ=1,capt=capt_,cont=cont_,unc=unc_,uncx=uncx_)
    else:
        gdout=GD.gd(y,x=x_,ini=ini_,dx=dx_,typ=2,capt=capt_,cont=cont_,unc=unc_,uncx=uncx_)

    return gdout



def gdsavemat7(fil,ingd):   ### NON FUNZIONA
    mdic={'ingd':ingd,'label':'gd by SnagPy'}
    sio.savemat(fil,mdic)



def gdsavedicmat(ingd):
# save a gd transformed to a dictionary
    nam=SERV.retrieve_name(ingd)
    outdic=GD.gd2dict(ingd)
    filnam=nam+'_dic.mat'
    sio.savemat(filnam,outdic)

    return 'File '+filnam+' created'



def gdloaddicmat(fil):
# load a gd transformed to a dictionary
    indic=sio.loadmat(fil)
    outgd=GD.dict2gd(indic)

    return outgd


# by Edoardo Giancarli ---------------------------


def mat_to_dict(pathfil):
# converts the data from MATLAB (v7 format) to dictionary
    
    """
    Conversion from MATLAB data file to dict.
    Parameters: 

    path : (str) path of your MATLAB data file
        
    data_dict: (dict) dict from MATLAB data file     
    """
    
    # SciPy reads in structures as structured NumPy arrays of dtype object
    # The size of the array is the size of the structure array, not the number-
    #   -elements in any particular field. The shape defaults to 2-dimensional.
    # For convenience make a dictionary of the data using the names from dtypes
    # Since the structure has only one element, but is 2-D, index it at [0, 0]
    
    mat = sio.loadmat(pathfil)                                 # load mat-file
    path,filnam,ext=BASIC.path_fil_ext(pathfil)
    mdata = mat[filnam]                                              # variable in mat file
    mdtype = mdata.dtype                                          # dtypes of structures are "unsized objects"
    data_dict = {n: mdata[n][0, 0] for n in mdtype.names}         # express mdata as a dict
        
    # n = data_dict['n'][0, 0]
    # data_dict['y'] = list(data_dict['y'][i][0] for i in range(n))
    y = data_dict['y']
    y = y.reshape(len(y))
    data_dict['y'] = y                               # perc of total zero data in y (data from the Obs run)
        
    cont = data_dict['cont']
    cont_dtype = cont.dtype
    cont_dict = {u: cont[str(u)] for u in cont_dtype.names}       # cont in data_dict is a structured array, I converted it in a dict
    data_dict['cont'] = cont_dict                                 # now we have a fully accessible dict
    
    return data_dict



def mat2gd(pathfil):                
# converts the data from MATLAB (v7 format) to a gd

    data_dict=mat_to_dict(pathfil)
    outgd=GD.dict2gd(data_dict)
    return outgd
 


# read - write v.7.3 format mat file ---------------------

def loadmat73(fil):
# Generic load
    dat=mat73.loadmat(fil)
    obj=list(dat.keys())

    return obj


def gdloadmat73(fil):
# load a file containing a gd
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

def csv_read(filn):
    with open(filn, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            print(', '.join(row))