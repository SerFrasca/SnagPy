    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

import sys
import csv
import json
import pickle
import h5py
import hdfdict 
import sys
import inspect
import os.path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors
import GD,GD2

"""
   To reimport a module:

 import importlib
 importlib.reload(Module)

 Matlab "exist" function:
 if var in locals():

"""

# Small apps ---------------------------------------

def Exec(file):   # exec(open('test_file').read())
    a="exec(open('"
    b="').read())"

    exec(a+file+b) 



def isa(arg,typ=0):
# Type of argument. 
# if typ is absent (or = 0) out is
    out='Unrecognized object'
    if isinstance(arg,int):
        out=1
    if isinstance(arg,float):
        out=2
    if isinstance(arg,complex):
        out=3
    if isinstance(arg,str):
        out=4
    if isinstance(arg,list):
        out=11
    if isinstance(arg,tuple):
        out=12
    if isinstance(arg,dict):
        out=13
    if isinstance(arg,set):
        out=14
    if isinstance(arg,np.ndarray):
        out=15
    if isinstance(arg,GD.gd):
        out=21
    if isinstance(arg,GD2.gd2):
        out=22
    elif typ != 0:
        out=isinstance(arg,typ)

    return out



# load & save dictionary: text, csv, json and pickle ---------------

def dict2text(dic,fil):
# write a dictionary on a text file
# fil  file without extension
    fil1=fil+'.txt'
    f=open(fil1,'w')
    f.write(str(dic))
    f.close()


def dict2csv(dic,fil):
# write a dictionary on a csv file
# fil  file without extension
    fil1=fil+'.csv'
    w=csv.writer(open(fil1,'w'))
    for key, val in dic.items():
        w.writerow([key,val])


def dict2json(dic,fil):  
# write a dictionary on a json file
# fil  file without extension
#  NOT WORKING WITH NUMPY DATA
    fil1=fil+'.json'
    json1=json.dumps(dic)
    f=open(fil1,'w')
    f.write(json1)
    f.close()


def dict2pkl(dic,fil):  
# write a dictionary on a pikle file
# fil  file without extension
    fil1=fil+'.pkl'
    f=open(fil1,'wb')
    pickle.dump(dic,f)
    f.close()


# HDF5 -------------------------

def list_baskeys_hdf5(fil):
    f=h5py.File(fil,'r')
    li=list(f.keys())
    print(li)

def read_hdf5(fil):
    out=hdfdict.load(fil)

    return out


def write_hdf5(dic,fil):
    hdfdict.dump(dic,fil)


# List ----------------------------

def list2string(lis):
    strin=''
    for x in lis:
        strin+=x

    return strin

def show_list(lis,sort=1):
# shows a list
# sort can be put = 1 only if the elements are all numbers or all strings
    print('    ')
    if sort > 0:
        lis.sort()
    k=0
    for i in lis:
        k+=1
        print('> ',k,' - ',i)
    print('    ')


# Dictionary --------------------------

def show_dict(dict):
    print('    ')
    for keys,values in dict.items():
        if isinstance(values,float):
            print(f'{keys:15} ==> {values:15f}')
        elif isinstance(values,int):
            print(f'{keys:15} ==> {values:15d}')
        elif isinstance(values,str):
            print(f'{keys:15} ==> {values:15}')
        else:
            print(f'{keys:15} ==>   xxx')
    print('    ')



def expl_dict(dic,Keys=[],Dics=[]):
   # def expl_dict(dic,dicname='-0',Keys=[],Dics=[]):
# explores dictionary structure
#  dic     dictionary to explore
#  Keys    inizialization of tne names of the variables
#  Val         "     "    of the values of the variables
#  Dics        "     "    of the names of the structures

    try:
        kkeys=tuple(dic.keys())
    except:
        print(' *** Not a dictionary')
        return Keys,Dics
    Keys.append(kkeys)
    if len(Dics) == 0:
        Dics.append('start')
    nval=len(kkeys)
    print(nval,kkeys)

    nstr=0
    for ii in range(nval):
        iikeys=kkeys[ii]
        # print('*** ',ii,nval)
        # print(' *** ',iikeys)
        val=dic[iikeys]
        if isinstance(val,dict):
            # print(' >>>1 ',type(val))
            # print(' >>>2 ',ii,iikeys)
            Dics.append(iikeys)
            k1,d1=expl_dict(val,Keys,Dics)
            nstr+=1

    return Keys,Dics



def show_dict_struct(Keys,Dics):
# displays the dictionary structure explored by expl_array 
    Rows=[]

    for i in range(len(Keys)):
        lis=[i,Dics[i],Keys[i]]
        Rows.append(lis)
        print(i,Dics[i],Keys[i])

    print(' - - - - - - - - - - - - - - - - - - - - ')

    Gen=[]
    for i in range(len(Rows)-1,0,-1):
        it=Rows[i][1]
        gen=[it]
        for ii in range(i-1,0,-1):
            if it in Rows[ii][2]:
                it=Rows[ii][1]
                gen.append(it)
        Gen.append(gen)

    for ii in range(len(Gen)):
        Gen[ii].reverse()

    Gen.reverse()

    for i in range(len(Keys)):
        if i == 0:
            print(' start ',Keys[i])    
        else:
            print(Gen[i-1],Keys[i])

    return Rows,Gen



def show_dict_struct_2(Keys,Dics):
# displays the dictionary structure explored by expl_dict

    for i in range(len(Keys)):
        print(i,Dics[i],Keys[i])



def dict_extract_2(dic,listkeys):
# extracts a value from a dictionary using the list of subsequent keys
#   dic       the dictionary
#   listkeys  list of the subsequent keys in the tree
    aa=dic
    try:
        for key in listkeys:
            aa=aa[key]

        return aa
    except:
            print(listkeys, 'not available')



# simple dict ----------------------
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

def eq_interpr(line):
# line equation interpretation for simple tables
    k=line.find('=')
    lin1=line[0:k].strip()
    lin2=line[k+1:-1].strip(' \n')

    if k == -1:
        lin1=''
        lin2=''
        typ=-1
    else:
        typ=1

    if len(lin1) > 0:
        if lin1[0] == '#':
            typ=0
    
    kk=lin2.find(',')

    if kk != -1:
        typ=2 

    return lin1,lin2,typ


def simp2dict(file):
# puts the content of a simple table to a simple dictionary
    tot='{'
    f=open(file,'r')
    lin=f.readline()
    while lin:
        if len(lin) > 0:
            lin1,lin2,typ=eq_interpr(lin)
            if typ > 0:
                tot=tot+'"'+lin1+'":'
            if typ == 1:
              tot=tot+lin2+','
            if typ == 2:
              tot=tot+lin2+','
        lin=f.readline()

    tot=tot[0:-1]+'}'
    tot

    return eval(tot)


def show_simp(sdic,spac=10,file=0):  
# shows a simple dictionary
#  spac is the key field length
#  file is the name of a desired output fil (if any, def no)
    for key, value in sdic.items():
        print(key.ljust(spac),' = ',value)
    if file != 0:
        stdout0=sys.stdout
        sys.stdout=open(file,'w')
        for key, value in sdic.items():
            print(key.ljust(spac),' = ',value)
        sys.stdout.close()
        sys.stdout=stdout0



# numpy structures -------------------------------------

def expl_array(arr,Names=[],Cval=[],StName=[]):
# explores array structure
#  Names   inizialization of tne names of the variables
#  Cval        "     "    of the values of the variables
#  StName      "     "    of the names of the structures

    names=arr.dtype.names
    try:
        cval=arr[0][0]
    except:
        print(' *** Not a structured array')
        return Names,Cval,StName
    Names.append(names)
    Cval.append(cval)
    if len(StName) == 0:
        StName.append('start')
    nval=len(names)
    print(names)
    dt=[]
    for ii in range(nval):
        aa=cval[ii]
        dt.append(len(aa.dtype))

    nstr=0
    for ii in range(nval):
        if dt[ii] > 0:
            stname=names[ii]
            StName.append(stname)
            n1,c1,st1=expl_array(cval[ii],Names,Cval,StName)
            nstr+=1

    return Names,Cval,StName



def expl_array_2(arr,arrname='-0',Names=[],Cval=[],StName=[]):
# explores array structure
#  arr     numpy ndarray to explore
#  arrname array name as string (optional)
#  Names   inizialization of tne names of the variables
#  Cval        "     "    of the values of the variables
#  StName      "     "    of the names of the structures

    if arrname == '-0':
        arrname='start'
    names=arr.dtype.names
    # cval=arr[0][0]
    try:
        cval=arr[0][0]
    except:
        print(' *** Not a structured array')
        return Names,Cval,StName
    Names.append(names)
    Cval.append(cval)
    if len(StName) == 0:
        StName.append(arrname)
    nval=len(names)
    print(names)
    dt=[]
    for ii in range(nval):
        aa=cval[ii]
        dt.append(len(aa.dtype))

    nstr=0
    for ii in range(nval):
        if dt[ii] > 0:
            stname=names[ii]
            StName.append(stname)
            n1,c1,st1=expl_array(cval[ii],Names,Cval,StName)
            nstr+=1

    return Names,Cval,StName



def show_array_struct_2(Names,StName):
# displays the array structure explored by expl_array 

    for i in range(len(Names)):
        print(i,StName[i],Names[i])



def show_array_struct(Names,StName):
# displays the array structure explored by expl_array 
    Rows=[]

    for i in range(len(Names)):
        lis=[i,StName[i],Names[i]]
        Rows.append(lis)
        print(i,StName[i],Names[i])

    print(' - - - - - - - - - - - - - - - - - - - - ')

    Gen=[]
    for i in range(len(Rows)-1,0,-1):
        it=Rows[i][1]
        gen=[it]
        for ii in range(i-1,0,-1):
            if it in Rows[ii][2]:
                it=Rows[ii][1]
                gen.append(it)
        Gen.append(gen)

    for ii in range(len(Gen)):
        Gen[ii].reverse()

    Gen.reverse()

    for i in range(len(Names)):
        if i == 0:
            print(' start ',Names[i])    
        else:
            print(Gen[i-1],Names[i])

    return Rows,Gen


def array_extract(kstr,var,Names,Cval,outdic=0):
# extracts a value from an array structure 
#   kstr          number of structure (as in show_array_struct)
#   var           as it is shown by show_array_struct
#   Names & Cval  as produced by expl_array
#   outdic        -1 output dictionary

    n1=Names[kstr]
    pos=n1.index(var)

    out=Cval[kstr][pos]
    try:
        if outdic == 1:
            out=arrstruct2dict(out)
    except:
        out={var: out.squeeze()}

    return out



def arrstruct2dict(arstr):
# dictionary from a numpy array structure
    asdtyp=arstr.dtype
    dict = {n: arstr[n][0, 0] for n in asdtyp.names}

    return dict


# SnagTable --------------------------

class snag_table:
# simple Matlab-like table management
#
#  data    full data (with titles, list of lists, all of the same length)
#  nr      number of rows
#  nc      number of colums
#  titles  first row of data
#  capt    caption
#  cont    control variable (typicalli a structure or dictionary)

    def __init__(self,data):
 #   def __init__(self,data,**stpar):
        self.data=data
        self.nr=len(data)-1
        self.nc=len(data[1])
        self.titles=data[0]
        self.capt='table'
        self.cont=0

def extr_st_col(st,which):
#  st     snag table
#  which  the number of the column or the title
    data=st.data
    titles=st.titles
    nr=st.nr
    print('nr ',nr)
    try:
        if isinstance(which,str):
            which=titles.index(which)
            print('column ',which)
    except:
        print(which,' erroneous key')
        print('Only ',titles)
        return 0

    outcol=[]
    for i in range(1,nr+1):
        outcol.append(data[i][which])

    if not isinstance(outcol[0],str):
        outcol=np.array(outcol)

    return outcol


def extr_st_rows(st,which):
# Creates a new table with a subset of rows
#  st     snag table
#  which  a list with the selected indices (ex.: [0,19,20,101])

    dataout=[]
    data=st.data
    dataout.append(data[0])
    ii=0

    for i in which:
        dataout.append(data[which[ii]+1])

    outst=snag_table(dataout)

    return outst



def csv2list(fil):
# reads dictionary from csv file
    lis=[]
    with open(fil, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            print(row)
            row1=decode_simp_list(row)
            lis.append(row1)

    return lis


def decode_simp_list(lis):
# decodes strings in numbers in simple lists

    out=[]

    for elem in lis:
        try:
            elem1=int(elem)
        except:
            try:
                elem1=float(elem)
            except:
                elem1=elem
        out.append(elem1)

    return out


# system -------------------------

def func_in_module(modul): 
# creates a list of the functions in a module
# the module should be imported
    res1=[]
    for i in dir(modul):
        if type(getattr(modul, i)).__name__ == "function":
            res1.append(getattr(modul, i))
    
    res=[]
    for i in res1:
        res.append(i.__name__)

    return res


def list_of_func(modul):
# creates a list of the functions in a module with addresses
# the module should be imported
        list_of_functions = inspect.getmembers(modul, inspect.isfunction)

        return list_of_functions



def deshape(inarr):
# Eliminates "shape" (and pletoric parentheses) in array
# reducing to simple 1-D array
    outarr=inarr
    aa=inarr.shape
    if len(aa) == 1:
        return outarr 
    elif aa[0] == 1:
        outarr=inarr[0]
    elif aa[1] == 1:
        outarr=outarr.transpose()
        outarr=outarr[0]
    else:
        n=aa[0]*aa[1]
        outarr=outarr.reshape(n)

    return outarr


def num_var_odd(lin):  
# analyze the string lin to identifies numbers or variable names
# 1 num, 2 var, 3 odd
    res=0
    lin1=lin.strip(' \n')
    try:
        lin0=float(lin1)
    except:
  #       ValueError
        lin2=lin1.replace('_','')
        res=3
        if lin2.isalnum():
            if lin2[0].isalpha():
                res=2
        lin0=lin1
    else:
        res=1

    return res,lin0



def path_fil_ext(ffil):
# identifies in a file name with path 
# the path, the file name without extension, the extension
    path,a2=os.path.split(ffil)
    filnam,ext=os.path.splitext(a2)

    return path,filnam,ext



def ana_module(modu,filout):  
# Modules analysis
#
# modu    path with the name of the module 
#         (if installed, see module_name.__file__)
#          e.g. BASIC.ana_module(GD.__file__,'prova.out'))
# filout  output file 

    fo=open(filout,'w')
    path,filnam,ext=path_fil_ext(modu)
    fo.write('         Module '+filnam+'\n\n')

    f=open(modu,'r')
    nchap=0
    nfun=0
    ncla=0
    ncom=0
    lins=0
    comon=0
    npass=0
    while True:
        line=f.readline()
        if line == '':
            strin1='\n {0}+{1} functions  {2} classes  {3} comment lines  {4} total lines\n'.format(nfun-npass,npass,ncla,ncom,lins)
            fo.write(strin1)
            fo.close()
            return
        lins+=1
        if ('"""' in line or "'''" in line) and not 'line' in line:
            if comon == 1:
                comon=0
            else:
                comon=1
        if comon == 1:
          #  fo.write(line[0:50]+'\n')
            fo.write(line)
            ncom+=1

        if 'pass' in line and len(line.strip()) == 4:
            npass+=1

        if '---' in line and 'line' not in line:
            nchap+=1
            line=line.rstrip()
            fo.write('\n'+'Section '+line[1:]+'\n\n')
            line=' '

        if line[0:3] == 'def':
            nfun+=1
            line=line.rstrip()
            strin1='\n{}   >fun {} row {}\n::\n'.format(line,nfun,lins)
            fo.write(strin1)
        if line[0:5] == 'class':
            ncla+=1
            line=line.rstrip()
            strin1='{}   >{} row {}\n::\n'.format(line,ncla,lins)
            fo.write(strin1)
        if line[0] == '#':
            line=line.rstrip()
            fo.write(line[1:]+'\n')
            ncom+=1
        if line[0:6] == 'import':
            line=line.rstrip()
            fo.write(line+'\n')


def all_modules(pack_path,filout): 
# All modules synthetic analysis
#
# pack_path   main path of the package 
# filout      output file 

    mod_list=['GD',
    'GD2',
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
    'PARGPU',
    'WEB_SNAG']
    # 'DEEPSNAG',
    # 'PERS//SF']

    fo=open(filout,'w')

    for modu in mod_list:
        fo.write('\n\n__________________________________________\n')
        fo.write('\n'+'         Module '+modu)

        f=open(pack_path+modu+'.py','r')
        nchap=0
        nfun=0
        ncla=0
        lins=0
        npass=0
        chapon=0
        while True:
            line=f.readline()
            if line == '':
                strin=list2string(funcs)
                fo.write('Functions: '+strin)
                # strin1='\n {0}+{1} functions  {2} classes  {3} comment lines  {4} total lines\n'.format(nfun-npass,npass,ncla,ncom,lins)
                # fo.write(strin1)
                # fo.close()
                break
            lins+=1
 
            if 'pass' in line and len(line.strip()) == 4:
                npass+=1

            if '---' in line and 'line' not in line:
                if chapon == 1:
                    strin=list2string(funcs)
                    fo.write('Functions: '+strin)
                funcs=[]
                chapon=1
                nchap+=1
                line=line.rstrip()
                fo.write('\n\n'+'Section '+line[1:]+'\n')
                line=' '

            if line[0:3] == 'def':
                nfun+=1
                k=line.find('(')
                funcs.append(line[4:k]+', ')
                # line=line.rstrip()
                # strin1='\n{}   >fun {} row {}\n::\n'.format(line,nfun,lins)
                # fo.write(strin1)
            if line[0:5] == 'class':
                ncla+=1
                k=line.find(':')
                fo.write(line[0:k]+'\n')
                # line=line.rstrip()
                # strin1='{}   >{} row {}\n::\n'.format(line,ncla,lins)
                # fo.write(strin1)

    fo.close()




# Graphic ------------------------

def fig_dim():
# figure dimensions
# outputs xlim and ylim
    axes = plt.gca()
    xx=axes.get_xlim()
    yy=axes.get_ylim()

    return xx,yy


def inter_grid(xx,lev=1):
    d=xx[1]-xx[0]
    n10=int(np.log10(d))
    d10=10**(n10-lev)
    x=np.ceil(xx[0]/d10)*d10
    print(d,n10,d10)
    gr=[]
    while x <= xx[1]:
        gr.append(x)
        x+=d10

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


