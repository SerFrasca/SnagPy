import sys
import csv
import json
import pickle
import sys
import inspect
import os.path

"""
   To reimport a module:

 import importlib
 importlib.reload(Module)

"""

def Exec(file):   # exec(open('test_file').read())
    a="exec(open('"
    b="').read())"

    exec(a+file+b) 


# load & save dictionary: text, csv, json and pickle

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



# simple dict ------------
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
        if '"""' in line or "'''" in line:
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

        if '---' in line:
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

 