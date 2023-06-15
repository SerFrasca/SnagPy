# Copyright (C) 2023  Sergio Frasca
#  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
        Module BASIC
Base functions for SnagPy

Sections:

> Small apps                                            -> dummy_small
> load & save dictionary: text, csv, json and pickle    -> dummy_gendic
> HDF5                                                  -> dummy_hdf5
> List                                                  -> dummy_list
> Dictionaries                                          -> dummy_dict
> simple dict                                           -> dummy_simpd
> Numpy arrays                                          -> dummy_numarr
> Numpy structures                                      -> dummy_numstr
> SnagTable                                             -> dummy_sntab
> ArrayTable                                            -> dummy_artab
> system                                                -> dummy_sys
> Graphic                                               -> dummy_graph
> Documentation                                         -> dummy_doc
'''

def sections():
    sec=[
        'small','gendic','hdf5','list','dict','simpd','numarr',
        'numstr','sntab','artab','sys','graph','doc'
    ]
    return sec

import sys
import csv
import json
import pickle
import h5py
import silx.io.dictdump as silx
import sys
import time
import copy
import inspect
# import os.path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors
#import matplotlib.colors as mcolorsNotFoundError
import GD
import GD2
import os,importlib

"""
   To reimport a module:

 import importlib
 importlib.reload(Module)

 Matlab "exist" function:
 if var in locals():

"""

# Small apps ---------------------------------------

def dummy_small():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def envir_var(var):
    '''
    value of environment variables
      var    name of the environment variable (a string. ex.: BASIC.envir_var('SNAGPY_PATH'))
    '''
    AA = dict(os.environ)
    try:
        BB = AA[var]
    except:
        print(var+' variable not present')
        BB = var+' variable not present'

    return BB


def var(v,typ=0):
    '''
    variable analysis
      v     input variable
      typ   =1 methods, =2 attributes, =3 methods and attributes
    '''
    if isinstance(v,np.ndarray):
        print('Numpy array')
        arr(v,typ=typ)
        return
    
    print('\n', type(v))
    print('id : ', id(v), '\n')
    print(v, '\n')

    siz=sys.getsizeof(v)
    print('Memory occupancy ',siz,' bytes \n')
    print('Recursive occupancy (ESTIMATION):')

    rsiz=len(pickle.dumps(v))
    print('   ',rsiz,'bytes')
    print('   ',rsiz/1024,'kbytes')
    print('   ',rsiz/1024**2,'Mbytes')
    print('   ',rsiz/1024**3,'Gbytes')

    if typ > 0:
        di = dir(v)

        if typ == 1 or typ == 3:
            i=0
            print('\n   Methods\n')
            for it in di:
                if hasattr(v, it) and callable(getattr(v, it)):
                    i+=1
                    print(i,it)

        if typ == 2 or typ == 3:
            i=0
            print('\n   Attributes\n')
            for it in di:
                if hasattr(v, it) and not callable(getattr(v, it)):
                    i+=1
                    print(i,it)



def Var(v):    
    '''
    Complete variable analysis
      v     input variabl
    '''
    var(v,typ=3)



def arr(v,typ=0):
    '''
    array analysis
      v     input variable
      typ   =1 methods, =2 attributes, =3 methods and attributes
    '''

    print('\n', type(v))
    print('id : ', id(v), '\n')
    print(v, '\n')
    print('dtype:',v.dtype)

    print('shape:',v.shape)
    print('itemsize:',v.itemsize,'bytes')

    siz=sys.getsizeof(v)
    print('Memory occupancy ',siz,' bytes \n')
    print('Recursive occupancy (ESTIMATION):')

    rsiz=len(pickle.dumps(v))
    print('   ',rsiz,'bytes')
    print('   ',rsiz/1024,'kbytes')
    print('   ',rsiz/1024**2,'Mbytes')
    print('   ',rsiz/1024**3,'Gbytes')

    if typ > 0:
        di = dir(v)

        if typ == 1 or typ == 3:
            i=0
            print('\n   Methods\n')
            for it in di:
                if hasattr(v, it) and callable(getattr(v, it)):
                    i+=1
                    print(i,it)

        if typ == 2 or typ == 3:
            i=0
            print('\n   Attributes\n')
            for it in di:
                if hasattr(v, it) and not callable(getattr(v, it)):
                    i+=1
                    print(i,it)




def Exec(file):   # exec(open('.\test_file').read())
    # just to see how to use exec.
    # Alternatively, launch at the command prompt python -i D:\OneDrive\SF\_Prog\Python\SnagPy\starting.py
    a = "exec(open('"
    # b = "').read(), globals(), locals())"
    b = "').read())"

    exec(a+file+b)

# aaastart="exec(open('"
# bbbend="').read())"


def isa(arg, typ=0):
    '''
    Type of argument.
    if typ is absent (or = 0) out is
    '''
    out = 'Unrecognized object'
    if isinstance(arg, int):
        out = 1
    if isinstance(arg, float):
        out = 2
    if isinstance(arg, complex):
        out = 3
    if isinstance(arg, str):
        out = 4
    if isinstance(arg, list):
        out = 11
    if isinstance(arg, tuple):
        out = 12
    if isinstance(arg, dict):
        out = 13
    if isinstance(arg, set):
        out = 14
    if isinstance(arg, np.ndarray):
        out = 15
    if isinstance(arg, GD.gd):
        out = 21
    if isinstance(arg, GD2.gd2):
        out = 22
    elif typ != 0:
        out = isinstance(arg, typ)

    return out


def dims(arr):
    '''
    Dimensions for arrays or gds

    arr     array or gd or gd2
    '''
    if isinstance(arr,GD.gd):
        arr=arr.y
    if isinstance(arr,GD2.gd2):
        arr=arr.y

    dim=arr.shape
    if len(dim) == 1:
        nr=1
        nc=dim[0]
    else:
        [nr,nc]=dim

    return nr,nc



def byte2str(dat):
    '''
    converts byte object or list of b.o. to str or list of str
    '''

    if isinstance(dat, bytes):
        dat = dat.decode("utf-8")
    elif isinstance(dat, list) or isinstance(dat, np.ndarray):
        ii = 0
        for it in dat:
            if isinstance(it, bytes):
                it = it.decode("utf-8")
                dat[ii] = it
            ii += 1

    return dat


def tic():
    '''
    similar to Matlab tic
    '''
    global ttic
    ttic = time.time()
    return ttic


def toc(tticin=0):
    '''
    similar to Matlab toc
    '''
    if tticin == 0:
        return time.time()-ttic
    else:
        return time.time()-tticin


def num_var_odd(lin):
    '''
    analyze the string lin to identifies numbers or variable names
    1 num, 2 var, 3 odd
    '''
    res = 0
    lin1 = lin.strip(' \n')
    try:
        lin0 = float(lin1)
    except:
      #       ValueError
        lin2 = lin1.replace('_', '')
        res = 3
        if lin2.isalnum():
            if lin2[0].isalpha():
                res = 2
        lin0 = lin1
    else:
        res = 1

    return res, lin0


def path_fil_ext(ffil):
    '''
    identifies in a file name with path
    the path, the file name without extension, the extension
    '''
    path, a2 = os.path.split(ffil)
    filnam, ext = os.path.splitext(a2)

    return path, filnam, ext


def listfil(pat,indir=1,opat=0):
    '''
    lists files in folders
    pat    path of the folder
    indir  =1 no subfolders, =0 also in subfolders
    opat   output path (together to files)
    '''

    print('opat = ',opat)
    res = []
    dirs=[]

    for path in os.listdir(pat):
        if os.path.isfile(os.path.join(pat, path)):
            if opat == 0:
                res.append(path)
            else:
                ap=opat+os.sep+path
                res.append(ap)
        else:
            dirs.append(path)

    if indir == 0:
        print('*** dirs',dirs)
        for it in dirs:
            os.chdir(it)
            print('+++>>',it)
            res1=listfil('.',opat=it)
            for it1 in res1:
                res.append(it1)
            os.chdir('..')

    return res


def ind_from_inidx(xx, ini, dx, typ=2):
    '''
    index from ini and dx
     xx      value
     ini,dx  parameters
     typ     1 floor, 2 round, 3 ceil
    '''
    aa = (xx-ini)/dx
    if typ == 2:
        ind = np.round(aa)
    elif typ < 2:
        ind = np.floor(aa)
    else:
        ind = np.ceil(aa)
    ind = int(ind)

    return ind


def ind_from_x(xx, x):
    '''
    index from x  (for type 2 gds)
     xx      value or sorted array or list
     x       abscissas
    '''

    lx=len(x)
    n=len(xx)
    ind=np.zeros(n)
    ii=0
    for i in range(n):
        if xx[i] <= x[0]:
            ind[i]=0
            ii+=1
        if xx[i] >= x[lx-1]:
            ind[i]=lx-1
    i=ii+1
    for k in range(lx):
        if k < lx-1:
            dx=(x[k+1]-x[k])/2
        if xx[i] >= x[k]-dx and xx[i] < x[k]+dx:
            ind[i]=k
            if i < n-1:
                i+=1

    return ind

# load & save dictionary: text, csv, json and pickle ---------------

def dummy_gendic():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def dict2text(dic, fil):
    '''
    write a dictionary on a text file
    fil  file without extension
    '''
    fil1 = fil+'.txt'
    f = open(fil1, 'w')
    f.write(str(dic))
    f.close()


def dict2csv(dic, fil):
    '''
    write a dictionary on a csv file
    fil  file without extension
    '''
    fil1 = fil+'.csv'
    w = csv.writer(open(fil1, 'w', delimiter=';'))
    for key, val in dic.items():
        w.writerow([key, val])


def dict2json(dic, fil):
    '''
    write a dictionary on a json file
    fil  file without extension
     NOT WORKING WITH NUMPY DATA
    '''
    fil1 = fil+'.json'
    json1 = json.dumps(dic)
    f = open(fil1, 'w')
    f.write(json1)
    f.close()


def dict2pkl(dic, fil):
    '''
    write a dictionary on a pikle file
    fil  file without extension
    '''
    fil1 = fil+'.pkl'
    f = open(fil1, 'wb')
    pickle.dump(dic, f)
    f.close()


def pkl2dict(fil):
    '''
    read a dictionary from a pickle file
    '''

    f = open(fil, 'rb')
    dic = pickle.load(f)
    f.close()

    return dic


# HDF5 -------------------------

def dummy_hdf5():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun


def list_baskeys_hdf5(fil):
    f = h5py.File(fil, 'r')
    li = list(f.keys())
    print(li)


def write_hdf5_s(dic, fil):
    '''
    writes in a hdf5 file a simple dictionary
     dats   simple dictionary containing the data
    '''
    N = len(dic)
    kk = list(dic.keys())
    vv = list(dic.values())
    with h5py.File(fil, 'w') as f:
        for i in range(N):
            nam = kk[i]
            val = vv[i]
            try:
                valsh = val.shape
                valdt = val.dtype
                dset = f.create_dataset(
                    nam, valsh, valdt, data=val, compression='gzip')
            except:
                dset = f.create_dataset(nam, data=val)


def read_hdf5_s(fil):
    '''
    writes in a hdf5 file a simple dictionary
     dats   simple dictionary containing the data
    '''
    dicout = {}
    with h5py.File(fil, 'r') as f:
        lk = list(f.keys())
        for i in range(len(lk)):
            dset = f[lk[i]]
    #        print(lk[i],type(dset))
            try:
                ll = len(dset)
                dset = dset[0:ll]
            except:
                dset = dset[()]
            dset = byte2str(dset)
            dicout[lk[i]] = dset

    return dicout


def write_dic2hdf5(dic, fil):
    '''
    writes a dictionary in a hdf5 file
    '''
    silx.dicttoh5(dic, fil, h5path='/', mode='a',
                  overwrite_data=False, create_dataset_args=None)


def read_hdf52dic(fil):
    '''
    reads a dictionary in a hdf5 file, transforming data in numpy data
    '''
    dic = silx.h5todict(fil, path='/')

    return dic


def read_hdf5_part(fil, dataset, slice):
    '''
    partial read of HDF5 file
      fil      file (with path)
      dataset  dataset (string)
      slice    string (ex.: '[1:20,[6,7]]')
    '''

    hf = h5py.File(fil, 'r')
    ds = hf.get(dataset)
    out = eval("hf[dataset]"+slice)

    return out



def explore_hdf5(fil):
    '''
    explores hdf5 frame data
    '''
    dic=read_hdf52dic(fil)
    Keys,Dics=expl_dict(dic)
    show_dict_struct(Keys, Dics)

    return Keys,Dics


# List ----------------------------

def dummy_list():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun


def list2string(lis):
    strin = ''
    for x in lis:
        strin += x

    return strin


def show_list(lis, sort=1):
    '''
    shows a list
    sort can be put = 1 only if the elements are all numbers or all strings
    '''
    print('    ')
    if sort > 0:
        lis.sort()
    k = 0
    for i in lis:
        k += 1
        print('> ', k, ' - ', i)
    print('    ')


def list_len(lis):
    '''
    Checks if lis is a list of lists and gives the length of the lists
    '''
    le=len(lis)
    lens=np.zeros(le)
    ii=0
    for l in lis:
        if isinstance(l,list):
            lens[ii]=len(l)
        else:
            print('*** the element',ii,'is not a list')
            return 0
        ii+=1

    return lens


def array_list_len(lis):
    '''
    Checks if lis is a list of arrays and gives the length of the arrays
    '''
    le=len(lis)
    lens=np.zeros(le,dtype=int)
    ii=0
    for l in lis:
        if isinstance(l,np.ndarray):
            lens[ii]=len(l)
        else:
            print('*** the element',ii,'is not an array')
            return 0
        ii+=1

    return lens


# Dictionary --------------------------

def dummy_dict():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun


def show_dict(dict):
    print('    ')
    for keys, values in dict.items():
        if isinstance(values, float):
            print(f'{keys:15} ==> {values:15f}')
        elif isinstance(values, int):
            print(f'{keys:15} ==> {values:15d}')
        elif isinstance(values, str):
            print(f'{keys:15} ==> {values:15}')
        else:
            print(f'{keys:15} ==>   xxx')
    print('    ')


def expl_dict(dic, Keys=[], Dics=[]):
    '''
    def expl_dict(dic,dicname='-0',Keys=[],Dics=[]):
        explores dictionary structure
        dic     dictionary to explore
        Keys    inizialization of tne names of the variables
        Val         "     "    of the values of the variables
        Dics        "     "    of the names of the structures
    '''

    try:
        kkeys = tuple(dic.keys())
    except:
        print(' *** Not a dictionary')
        return Keys, Dics
    Keys.append(kkeys)
    if len(Dics) == 0:
        Dics.append('start')
    nval = len(kkeys)
    print(nval, kkeys)

    nstr = 0
    for ii in range(nval):
        iikeys = kkeys[ii]
        # print('*** ',ii,nval)
        # print(' *** ',iikeys)
        val = dic[iikeys]
        if isinstance(val, dict):
            # print(' >>>1 ',type(val))
            # print(' >>>2 ',ii,iikeys)
            Dics.append(iikeys)
            k1, d1 = expl_dict(val, Keys, Dics)
            nstr += 1

    return Keys, Dics


def show_dict_struct(Keys, Dics):
    '''
    displays the dictionary structure explored by expl_array
    '''
    Rows = []

    for i in range(len(Keys)):
        lis = [i, Dics[i], Keys[i]]
        Rows.append(lis)
        print(i, Dics[i], Keys[i])

    print(' - - - - - - - - - - - - - - - - - - - - ')

    Gen = []
    for i in range(len(Rows)-1, 0, -1):
        it = Rows[i][1]
        gen = [it]
        for ii in range(i-1, 0, -1):
            if it in Rows[ii][2]:
                it = Rows[ii][1]
                gen.append(it)
        Gen.append(gen)

    for ii in range(len(Gen)):
        Gen[ii].reverse()

    Gen.reverse()

    for i in range(len(Keys)):
        if i == 0:
            print(' start ', Keys[i])
        else:
            print(Gen[i-1], Keys[i])

    return Rows, Gen


def show_dict_struct_2(Keys, Dics):
    '''
    displays the dictionary structure explored by expl_dict
    '''

    for i in range(len(Keys)):
        print(i, Dics[i], Keys[i])


def dict_extract_2(dic, listkeys):
    '''
    extracts a value from a dictionary using the list of subsequent keys
      dic       the dictionary
      listkeys  list of the subsequent keys in the tree
    '''
    aa = dic
    try:
        for key in listkeys:
            aa = aa[key]

        return aa
    except:
        print(listkeys, 'not available')


# simple dict ----------------------

def dummy_simpd():
    '''
    A simple dictionary is something that describes something like
    a simple table as:

    Antenna     =  "Virgo"
    Long        =  12.5
    Lat         =  42
    Azimut      =  30
    Altro       =  [1, 23, 4.5]
    CC          =  (11+4.2j)

    '''

    clas=[

    ]
    fun=[
        
    ]
    return clas,fun


class simp_dict:
    '''
    A simple_dict is a "decorated" dictionary, with some limitations:
    - the keys are always strings
    - the values are never dictionaries (or simp_dict)

    Attributes:

    > dic       dictionary
    > keys
    > values
    > len       length
    > name
    > capt
    '''
    def __init__(self,dic,**gdpar): 
        self.dict=dic
        keys=dic.keys()
        values=list(dic.values())
        N=len(dic)
        self.len=N
        self.keys=keys
        self.values=values
        flavor=np.zeros(N,dtype=int)
        for i in range(N):
            flavor[i]=isa(values[i])
        self.flavor=flavor
        self.level=0
        ok=1
        for i in keys:
            if not isinstance(i,str):
                ok=0
                print('key ',i,' not a string')
        ii=0
        values=list(values)
        for i in values:
            if isinstance(i,simp_dict):
                self.level=1+i.level
            if isinstance(i,dict):
            # if isinstance(i,dict):
                ok=0
                print('value ',ii,' not permitted')
            ii+=1

        if ok == 0:
            print(' *** simp_dict not correct')

        self.ok=ok

        if 'name' in gdpar:
            self.name=gdpar['name']
        else:
            self.name='simp_dict'

        if 'capt' in gdpar:
            self.capt=gdpar['capt']
        else:
            self.capt='no capt'



def eq_interpr(line):
    '''
    line equation interpretation for simple tables
    '''
    k = line.find('=')
    lin1 = line[0:k].strip()
    lin2 = line[k+1:].strip(' \n')

    if k == -1:
        lin1 = ''
        lin2 = ''
        typ = -1
    else:
        typ = 1

    if len(lin1) > 0:
        if lin1[0] == '#':
            typ = 0

    kk = lin2.find(',')

    if kk != -1:
        typ = 2

    return lin1, lin2, typ


def simp2sdict(file):
    '''
    puts the content of a simple table to a simple dictionary
    '''
    tot = '{'
    f = open(file, 'r')
    lin = f.readline()
    while lin:
        if len(lin) > 0:
            lin1, lin2, typ = eq_interpr(lin)
            if typ > 0:
                tot = tot+'"'+lin1+'":'
            if typ == 1:
                tot = tot+lin2+','
            if typ == 2:
                tot = tot+lin2+','
        lin = f.readline()

    tot = tot[0:-1]+'}'
    tot

    dic=eval(tot)
    path,filnam,ext=path_fil_ext(file)
    sdict=simp_dict(dic,name=filnam)

    return sdict


def show_simp(sdic, spac=10, file=0, tit=''):
    '''
    shows a simple dictionary
     spac is the key field length
     file is the name of a desired output fil (if any, def no)
    '''
    print(tit)
    print(sdic.name)
    print(sdic.capt)
    print(sdic.level)
    print(sdic.ok)

    for key, value in sdic.items():
        print(key.ljust(spac), ' = ', value)
    if file != 0:
        stdout0 = sys.stdout
        sys.stdout = open(file, 'w')
        for key, value in sdic.items():
            print(key.ljust(spac), ' = ', value)
        sys.stdout.close()
        sys.stdout = stdout0


# numpy arrays -----------------------------------

def dummy_numarr():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun


def array_rowcol(arr):
    '''
    number of row, columns for np.arrays or gd or gd2
    '''
    if isinstance(arr,GD.gd):
        arr=arr.y
    if isinstance(arr,GD2.gd2):
        arr=arr.y
    dsh=arr.shape
    if len(dsh) == 1:
        nr = 1
        nc = len(arr)
    else:
        nr = dsh[0]
        nc = dsh[1]
    return nr,nc


# numpy structures -------------------------------------

def dummy_numstr():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def expl_array(arr, Names=[], Cval=[], StName=[]):
    '''
    explores array structure
     Names   inizialization of tne names of the variables
     Cval        "     "    of the values of the variables
     StName      "     "    of the names of the structures
    '''

    names = arr.dtype.names
    try:
        cval = arr[0][0]
    except:
        print(' *** Not a structured array')
        return Names, Cval, StName
    Names.append(names)
    Cval.append(cval)
    if len(StName) == 0:
        StName.append('start')
    nval = len(names)
    print(names)
    dt = []
    for ii in range(nval):
        aa = cval[ii]
        dt.append(len(aa.dtype))

    nstr = 0
    for ii in range(nval):
        if dt[ii] > 0:
            stname = names[ii]
            StName.append(stname)
            n1, c1, st1 = expl_array(cval[ii], Names, Cval, StName)
            nstr += 1

    return Names, Cval, StName


def expl_array_2(arr, arrname='-0', Names=[], Cval=[], StName=[]):
    '''
    explores array structure
     arr     numpy ndarray to explore
     arrname array name as string (optional)
     Names   inizialization of tne names of the variables
     Cval        "     "    of the values of the variables
     StName      "     "    of the names of the structures
    '''

    if arrname == '-0':
        arrname = 'start'
    names = arr.dtype.names
    # cval=arr[0][0]
    try:
        cval = arr[0][0]
    except:
        print(' *** Not a structured array')
        return Names, Cval, StName
    Names.append(names)
    Cval.append(cval)
    if len(StName) == 0:
        StName.append(arrname)
    nval = len(names)
    print(names)
    dt = []
    for ii in range(nval):
        aa = cval[ii]
        dt.append(len(aa.dtype))

    nstr = 0
    for ii in range(nval):
        if dt[ii] > 0:
            stname = names[ii]
            StName.append(stname)
            n1, c1, st1 = expl_array(cval[ii], Names, Cval, StName)
            nstr += 1

    return Names, Cval, StName


def show_array_struct_2(Names, StName):
    '''
    displays the array structure explored by expl_array
    '''

    for i in range(len(Names)):
        print(i, StName[i], Names[i])


def show_array_struct(Names, StName):
    '''
    displays the array structure explored by expl_array
    '''
    Rows = []

    for i in range(len(Names)):
        lis = [i, StName[i], Names[i]]
        Rows.append(lis)
        print(i, StName[i], Names[i])

    print(' - - - - - - - - - - - - - - - - - - - - ')

    Gen = []
    for i in range(len(Rows)-1, 0, -1):
        it = Rows[i][1]
        gen = [it]
        for ii in range(i-1, 0, -1):
            if it in Rows[ii][2]:
                it = Rows[ii][1]
                gen.append(it)
        Gen.append(gen)

    for ii in range(len(Gen)):
        Gen[ii].reverse()

    Gen.reverse()

    for i in range(len(Names)):
        if i == 0:
            print(' start ', Names[i])
        else:
            print(Gen[i-1], Names[i])

    return Rows, Gen


def array_extract(kstr, var, Names, Cval, outdic=0):
    '''
    extracts a value from an array structure
      kstr          number of structure (as in show_array_struct)
      var           as it is shown by show_array_struct
      Names & Cval  as produced by expl_array
      outdic        -1 output dictionary
    '''

    n1 = Names[kstr]
    pos = n1.index(var)

    out = Cval[kstr][pos]
    try:
        if outdic == 1:
            out = arrstruct2dict(out)
    except:
        out = {var: out.squeeze()}

    return out


def arrstruct2dict(arstr):
    '''
    dictionary from a numpy array structure
    '''
    asdtyp = arstr.dtype
    dict = {n: arstr[n][0, 0] for n in asdtyp.names}

    return dict


def val_from_key(st, kk):
    '''
    value from key for dictionary or structured array
       st   structured array or dictionary
       kk   key (string)
    '''

    if isinstance(st, dict):
        val = st[kk]
    else:
        dic = arrstruct2dict(st)
        val = dic[kk]
        val = val.squeeze()

    return val


# SnagTable --------------------------

def dummy_sntab():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

class snag_table:
    '''
    simple Matlab-like table management
    
     data    full data (with titles, list of lists or tuples, all of the same length)
     nr      number of rows
     nc      number of colums
     tup     = 1 list of tuples, = 0 list of lists
     titles  first row of data
     name    name
     capt    caption
     cont    meta data or control variable (typically a structure or dictionary)
    '''

    def __init__(self, data, name='snag_table',capt='table', cont=0, tup=1):
        self.data = data
        self.nr = len(data)-1
        self.nc = len(data[1])
        flavor=np.zeros([self.nr,self.nc],dtype=int)
        for i in range(self.nr): 
            for j in range(self.nc):
                flavor[i][j]=isa(data[i+1][j])
        self.flavor=flavor
        self.tup = tup
        self.titles = data[0]
        self.name = name
        self.capt = capt
        self.cont = cont


def snag_table_show(st):
    print(st.titles)
    for i in range(1, st.nr+1):
        print(i-1, st.data[i])


def extr_st_col(st, which):
    '''
     st     snag table
     which  the number of the column or the title
    '''
    data = st.data
    titles = st.titles
    nr = st.nr
    print('nr ', nr)
    try:
        if isinstance(which, str):
            which = titles.index(which)
            print('column ', which)
    except:
        print(which, ' erroneous key')
        print('Only ', titles)
        return 0

    outcol = []
    for i in range(1, nr+1):
        outcol.append(data[i][which])

    if not isinstance(outcol[0], str):
        outcol = np.array(outcol)

    return outcol


def extr_st_rows(st, which):
    '''
    Creates a new table with a subset of rows
     st     snag table
     which  a list with the selected indices (ex.: [0,19,20,101])
    '''

    dataout = []
    data = st.data
    dataout.append(data[0])
    ii = 0

    for i in which:
        dataout.append(data[which[ii]+1])

    outst = snag_table(dataout)

    return outst


def extr_st_row(st, k):
    '''
    Creates a new table with a subset of rows
     st     snag table
     k      the index of the row
    '''

    dataout = []
    data = st.data
    dataout.append(data[k+1])

    return dataout[0]


def extr_st_dict(st, which):
    '''
     st     snag table
     which  a 2-values list, with the numbers of the columns uses for keys and values
    '''
    ke = extr_st_col(st, which[0])
    va = extr_st_col(st, which[1])

    dic = dict(zip(ke, va))

    return dic


def csv2list(fil, tup=1):
    '''
    reads data from csv file
     tup     = 1 list of tuples, = 0 list of lists
    '''
    lis = []
    with open(fil, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            print(row)
            row1 = decode_simp_list(row)
            if tup == 1:
                row1 = tuple(row1)
            lis.append(row1)

    return lis


def csv2st(fil, tup=1):
    '''
    creates a snag_table from csv file
     tup     = 1 list of tuples, = 0 list of lists
    '''
    data = csv2list(fil, tup=1)
    path,filnam,ext=path_fil_ext(fil)
    st = snag_table(data, tup=tup, name=filnam)

    return st


def st2csv(st, fil):
    '''
    store a snag_table in a csv file
    fil  file without extension
    '''
    data = st.data
    fil1 = fil+'.csv'
    w = csv.writer(open(fil1, 'w', newline=''), delimiter=';')
    for it in data:
        w.writerow(it)


def decode_simp_list(lis):
    '''
    decodes strings in numbers in simple lists
    '''

    out = []

    for elem in lis:
        try:
            elem1 = int(elem)
        except:
            try:
                elem1 = float(elem)
            except:
                elem1 = elem
        out.append(elem1)

    return out

# ArrayTable --------------------------

def dummy_artab():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

class array_table:
    '''
    container for huge numerical table
    
     titles  names of the columns
     data    full data (nr x nc array)
     nr      number of rows
     nc      number of colums
     capt    caption
     cont    meta data or control variable (typically a structure or dictionary)
    '''

    def __init__(self, titles, data):
        self.titles = titles
        self.data = data
        self.nr = len(data)
        self.nc = len(data[1])
        self.capt = 'table'
        self.cont = 0


def text2array(fil, noline, nr, items):
    '''
    to read a table of floats
     fil     file
     noline  number of lenes to jump
     nr      number of output rows
     items   column to output (e.g. [0,1,3,7])
    '''

    nc = len(items)
    f = open(fil, 'r')
    arr = np.zeros((nr, nc))

    for i in range(noline):
        f.readline()

    for i in range(nr):
        lin = f.readline()
        lis = lin.split()
        num = []
        for k in range(len(lis)):
            num.append(float(lis[k]))

        num1 = []
        for k in items:
            num1.append(num[k])

        arr[i] = num1

    return arr


def array_table_to_dict(at, fil=''):
    '''
     at    array_table
     fil   output file (HDF5)
    '''

    atdic = {'capt': at.capt, 'titles': at.titles, 'nr': at.nr, 'nc': at.nc, 'array': at.data,
             'cont': at.cont}

    if fil != '':
        write_hdf5_s(atdic, fil)

    return atdic


# system -------------------------

def dummy_sys():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun


def func_in_module(modul):
    '''
    creates a list of the functions in a module
    the module should be imported
    '''
    res1 = []
    for i in dir(modul):
        if type(getattr(modul, i)).__name__ == "function":
            res1.append(getattr(modul, i))

    res = []
    for i in res1:
        res.append(i.__name__)

    return res


def deshape(inarr):
    '''
    Eliminates "shape" (and pletoric parentheses) in array
    reducing to simple 1-D array
    It is similar to squeeze
    '''
    outarr = inarr
    aa = inarr.shape
    if len(aa) == 1:
        return outarr
    elif aa[0] == 1:
        outarr = inarr[0]
    elif aa[1] == 1:
        outarr = outarr.transpose()
        outarr = outarr[0]
    else:
        n = aa[0]*aa[1]
        outarr = outarr.reshape(n)

    return outarr


def ana_module(modu, filout):
    '''
    Modules analysis
    
    modu    path with the name of the module
            (if installed, see module_name.__file__)
             e.g. BASIC.ana_module(GD.__file__,'prova.out'))
    filout  output file
    '''

    fo = open(filout, 'w')
    path, filnam, ext = path_fil_ext(modu)
    fo.write('         Module '+filnam+'\n\n')

    f = open(modu, 'r')
    nchap = 0
    nfun = 0
    ncla = 0
    ncom = 0
    lins = 0
    comon = 0
    npass = 0
    while True:
        line = f.readline()
        if line == '':
            strin1 = '\n {0}+{1} functions  {2} classes  {3} comment lines  {4} total lines\n'.format(
                nfun-npass, npass, ncla, ncom, lins)
            fo.write(strin1)
            fo.close()
            return
        lins += 1
        if ('"""' in line or "'''" in line) and not 'line' in line:
            if comon == 1:
                comon = 0
            else:
                comon = 1
        if comon == 1:
          #  fo.write(line[0:50]+'\n')
            fo.write(line)
            ncom += 1

        if 'pass' in line and len(line.strip()) == 4:
            npass += 1

        if '---' in line and 'line' not in line:
            nchap += 1
            line = line.rstrip()
            fo.write('\n'+'Section '+line[1:]+'\n\n')
            line = ' '

        if line[0:3] == 'def':
            nfun += 1
            line = line.rstrip()
            strin1 = '\n{}   >fun {} row {}\n::\n'.format(line, nfun, lins)
            fo.write(strin1)
        if line[0:5] == 'class':
            ncla += 1
            line = line.rstrip()
            strin1 = '{}   >{} row {}\n::\n'.format(line, ncla, lins)
            fo.write(strin1)
        if line[0] == '#':
            line = line.rstrip()
            fo.write(line[1:]+'\n')
            ncom += 1
        if line[0:6] == 'import':
            line = line.rstrip()
            fo.write(line+'\n')


def all_modules(pack_path, filout):
    '''
    All modules synthetic analysis
    
    pack_path   main path of the package
    filout      output file
    
    ex.: BASIC.all_modules('D:\\OneDrive\\SF\\_Prog\\Python\\SnagPy\\','prova.out')
    '''

    mod_list = ['GD',
                'GD2',
                'MGD',
                'DS',
                'BASIC',
                'SERV',
                'ML_PY',
                'STAT',
                'SIGNAL',
                'IMAGE',
                'ASTROTIME',
                'GWDATA',
                'PSS',
                'BSD',
                'GWOTH',
                'GUISNAG',
                'EXT_PACK',
                'MAN_SUPER',
                'PARGPU',
                'WEB_SNAG',
                'FANCY_FIG',
                'DEEPSNAG']
    # 'PERS//SF']

    print(pack_path, filout)

    fo = open(filout, 'w')

    for modu in mod_list:
        funcs = []
        fo.write('\n\n__________________________________________\n')
        fo.write('\n'+'         Module '+modu+'\n')

        print(pack_path+modu+'.py')

        f = open(pack_path+modu+'.py', 'r')
        nchap = 0
        nfun = 0
        ncla = 0
        lins = 0
        npass = 0
        chapon = 0
        while True:
            line = f.readline()
            if line == '':
                strin = list2string(funcs)
                fo.write('Functions: '+strin)
                funcs = []
                # strin1='\n {0}+{1} functions  {2} classes  {3} comment lines  {4} total lines\n'.format(nfun-npass,npass,ncla,ncom,lins)
                # fo.write(strin1)
                # fo.close()
                break
            lins += 1

            if 'pass' in line and len(line.strip()) == 4:
                npass += 1

            if '---' in line and 'line' not in line:
                if chapon == 1:
                    strin = list2string(funcs)
                    fo.write('Functions: '+strin)
                funcs = []
                chapon = 1
                nchap += 1
                line = line.rstrip()
                fo.write('\n\n'+'Section '+line[1:]+'\n')
                line = ' '

            if line[0:3] == 'def':
                nfun += 1
                k = line.find('(')
                funcs.append(line[4:k]+', ')
                # line=line.rstrip()
                # strin1='\n{}   >fun {} row {}\n::\n'.format(line,nfun,lins)
                # fo.write(strin1)
            if line[0:5] == 'class':
                ncla += 1
                k = line.find(':')
                fo.write(line[0:k]+'\n')
                # line=line.rstrip()
                # strin1='{}   >{} row {}\n::\n'.format(line,ncla,lins)
                # fo.write(strin1)

    fo.close()

    
def list_modules(pack_path, filout):
    '''
    All modules synthetic analysis
    
    pack_path   main path of the package
    filout      output file
    
    ex.: BASIC.all_modules('D:\\OneDrive\\SF\\_Prog\\Python\\SnagPy\\','prova.out')
    '''

    mod_list = ['GD',
                'GD2',
                'MGD',
                'DS',
                'BASIC',
                'SERV',
                'ML_PY',
                'STAT',
                'SIGNAL',
                'IMAGE',
                'ASTROTIME',
                'GWDATA',
                'PSS',
                'BSD',
                'GWOTH',
                'GUISNAG',
                'EXT_PACK',
                'MAN_SUPER',
                'PARGPU',
                'WEB_SNAG',
                'FANCY_FIG',
                'DEEPSNAG']
    # 'PERS//SF']

    print(pack_path, filout)

    fo = open(filout, 'w')

    for modu in mod_list:
        funcs = []
        fo.write('\n\n__________________________________________\n')
        fo.write('\n'+'         Module '+modu+'\n')

        print(pack_path+modu+'.py')

        f = open(pack_path+modu+'.py', 'r')
        nchap = 0
        nfun = 0
        ncla = 0
        lins = 0
        npass = 0
        chapon = 0
        while True:
            line = f.readline()
            if line == '':
                strin = list2string(funcs)
                # fo.write('Functions: '+strin)
                funcs = []
                # strin1='\n {0}+{1} functions  {2} classes  {3} comment lines  {4} total lines\n'.format(nfun-npass,npass,ncla,ncom,lins)
                # fo.write(strin1)
                # fo.close()
                break
            lins += 1

            if 'pass' in line and len(line.strip()) == 4:
                npass += 1

            if '---' in line and 'line' not in line:
                if chapon == 1:
                    strin = list2string(funcs)
                    # fo.write('Functions: '+strin)
                funcs = []
                chapon = 1
                nchap += 1
                line = line.rstrip()
                fo.write('\n\n'+'Section '+line[1:]+'\n')
                line = ' '

            if line[0:3] == 'def':
                nfun += 1
                k = line.find('(')
                fo.write(line[4:k]+'\n')
                # line=line.rstrip()
                # strin1='\n{}   >fun {} row {}\n::\n'.format(line,nfun,lins)
                # fo.write(strin1)
            if line[0:5] == 'class':
                ncla += 1
                k = line.find(':')
                fo.write(line[0:k]+'\n')
                # line=line.rstrip()
                # strin1='{}   >{} row {}\n::\n'.format(line,ncla,lins)
                # fo.write(strin1)

    fo.close()



# Graphic ------------------------

def dummy_graph():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun


def fig_dim():
    '''
    figure dimensions
    outputs xlim and ylim
    '''
    axes = plt.gca()
    xx = axes.get_xlim()
    yy = axes.get_ylim()

    return xx, yy


def inter_grid(xx, lev=1):
    '''
    to define grids
    '''
    d = xx[1]-xx[0]
    n10 = int(np.log10(d))
    d10 = 10**(n10-lev)
    x = np.ceil(xx[0]/d10)*d10
    print(d, n10, d10)
    gr = []
    while x <= xx[1]:
        gr.append(x)
        x += d10

    return gr


def plot_colortable(colors, sort_colors=True, emptycols=0):

    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12

    # Sort colors by hue, saturation, value and name.
    if sort_colors is True:
        by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(color))),
                         name)
                        for name, color in colors.items())
        names = [name for hsv, name in by_hsv]
    else:
        names = list(colors)

    n = len(names)
    ncols = 4 - emptycols
    nrows = n // ncols + int(n % ncols > 0)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + 2 * margin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-margin)/height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        text_pos_x = cell_width * col + swatch_width + 7

        ax.text(text_pos_x, y, name, fontsize=14,
                horizontalalignment='left',
                verticalalignment='center')

        ax.add_patch(
            Rectangle(xy=(swatch_start_x, y-9), width=swatch_width,
                      height=18, facecolor=colors[name], edgecolor='0.7')
        )

    return fig


def base_colors():
    plot_colortable(mcolors.BASE_COLORS, sort_colors=False, emptycols=1)


def tab_palette():
    plot_colortable(mcolors.TABLEAU_COLORS, sort_colors=False, emptycols=2)


def css_colors():
    plot_colortable(mcolors.CSS4_COLORS)


# Documentation -------------------------

def dummy_doc():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun


def SnagH(modul,sect=0,func=0):
    '''
    SnagPy Help on module, section and function

        modul   Module
        sect    Module section
        func    function

    You can also use help() function
    '''



def list_of_func(modul):
    '''
    creates a list of the functions in a module with addresses
    the module should be imported
    '''
    list_of_functions = inspect.getmembers(modul, inspect.isfunction)

    return list_of_functions