    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

import numpy as np
import ASTROTIME,BASIC
import os

Pi=np.pi

# General -----------------------------------------

def set_symbols(snagpy_p=0):
# Symbols to folders, files and data
#  snagpypath   the path to snagpy, as a string (es.: 'D:\\OneDrive\\SF\\_Prog\\Python\\SnagPy')
#               if absent, taken by the environment variable SNAGPY_PATH

    sep=os.sep

    if snagpy_p == 0:
        snagpy_p=BASIC.envir_var('SNAGPY_PATH')
    
    gwdata_p=snagpy_p+sep+'GWdata'
    examples_p=snagpy_p+sep+'Examples'
    exper_p=snagpy_p+sep+'Exper'
    pers_p=snagpy_p+sep+'Pers'
    projects_p=snagpy_p+sep+'Projects'

    antennas_p=gwdata_p+sep+'Antennas.csv'
    cwinj_o2_p=gwdata_p+sep+'cwinj_O2.csv'
    cwinj_o3_p=gwdata_p+sep+'cwinj_O3.csv'
    cwsour_p=gwdata_p+sep+'cwsour.csv'
    table_ligoh_p=gwdata_p+sep+'table_ligoh.hdf5'
    table_ligol_p=gwdata_p+sep+'table_ligol.hdf5'
    table_virgo_p=gwdata_p+sep+'table_virgo.hdf5'
    table_kagra_p=gwdata_p+sep+'table_kagra.hdf5'
    de440s_p=gwdata_p+sep+'de440s.bsp'

    return snagpy_p,gwdata_p,examples_p,exper_p,pers_p,projects_p,\
        antennas_p,cwinj_o2_p,cwinj_o3_p,cwsour_p,\
        table_ligoh_p,table_ligol_p,table_virgo_p,table_kagra_p,de440s_p


# Doppler data --------------------------------

def extr_doppler(tab,tin,tfi,table_par):
# extract data from doppler tables (in hdf5 format)
#   tab        Doppler table (ex.: table_virgo_p)
#   tin,tfi    times (vect as [a,m,d,h,m,s])
#   table_par  as defined in the symbols or starting.py

    gtin=ASTROTIME.now(tin,form='gps')
    gtfi=ASTROTIME.now(tfi,form='gps')

    ini=int((gtin-table_par[0])/table_par[2])
    ifi=int((gtfi-table_par[0])/table_par[2])+1

    out=BASIC.read_hdf5_part(tab,'array','['+str(ini)+':'+str(ifi)+']')

    return out
    

# CW sources ---------------------------------------

def ligo2virgo_cw_table(ligotab,capttab='',run=''):
# creates a Virgo format cw list from a Ligo table
    ligtab=BASIC.csv2list(ligotab)
    virgotab=[]
    title=[]
    lin=[]
    t00=ASTROTIME.v2mjd([2000,1,1,0,0,0])

    if run != '':
        title=('run','name','a','d','v_a','v_d','fepoch','f0','df0','ddf0','pepoch',
        't00','eps','eta','psi','h')
    else:
        title=('name','a','d','v_a','v_d','fepoch','f0','df0','ddf0','pepoch',
        't00','eps','eta','psi','h')

    virgotab.append(title)
    N=len(ligtab)

    for i in range(1,N):
        lline=ligtab[i]
        print(lline)
        ii=i-1
        a=lline[16]*180/Pi
        d=lline[15]*180/Pi
        v_a=0   # marcs/y
        v_d=0   # marcs/y
        fepoch=0 # to be corrected after gps2mjd
        f0=lline[2]
        df0=lline[3]
        ddf0=0
        pepoch=fepoch
        eps=1
        cosi=np.cos(lline[11])
        eta=-2*cosi/(1+cosi**2)
        psi=lline[12]*180/Pi
        h=lline[8]*np.sqrt(1+6*cosi**2+cosi**4)/2

        if run != '':
            lin=(run,str(ii),a,d,v_a,v_d,fepoch,f0,df0,ddf0,pepoch,t00,eps,eta,psi,h)
        else:
            lin=(str(ii),a,d,v_a,v_d,fepoch,f0,df0,ddf0,pepoch,t00,eps,eta,psi,h)

        virgotab.append(lin)

        vt=BASIC.snag_table(virgotab,capt=capttab+run)

        return vt







# Antennas ------------------------- 

def anten(ant,anten_tab=0):
# extract data for an antenna (a dictionary)
#  ant         antenna name ('virgo','ligol',...)
#  anten_tab   antenna table (def possible if the environment variable is set)
    if anten_tab == 0:
        anten_tab=BASIC.envir_var('ANT_TAB')

    lis=BASIC.csv2list(anten_tab)
    st=BASIC.snag_table(lis)
    tit=list(st.titles)
    tit[0]='Antenna'
    col=BASIC.extr_st_col(st,0)
    try:
        n=col.index(ant)
    except:
        print(ant+' not found in table'+anten_tab)
        return 0

    antlis=list(BASIC.extr_st_row(st,n))

    antdic=dict(zip(tit,antlis))

    return antdic



def antennas(anten_tab=0):
# creates dictionarys for each antenna
#  typical call: virgo,ligol,ligoh,kagra=GWDATA.antennas()
    virgo=anten('virgo')
    ligol=anten('ligol')
    ligoh=anten('ligoh')
    kagra=anten('kagra')

    return virgo,ligol,ligoh,kagra



# 5-vect --------------------------------------




# Other ---------------------------------------

def conv_jpl_data(datin,noline,nr,items):
# converts text data (typically table_xxxx.dat) to hdf5 format
#   datin    complete path
#   noline   number of lenes to jump (typ. 4)
#   nr       number of output rows   (typ. 526032)
#   items    column to output (e.g. [0,1,3,7], typ. [1,2,3,4,5,6,7,8])

    sep=os.sep

    path,fil,ext=BASIC.path_fil_ext(datin)

    arr=BASIC.text2array(path+sep+fil+ext,noline,nr,items)
    at=BASIC.array_table(['gpst','x','y','z','vx','vy','vz','einst'],arr)
    dic=BASIC.array_table_to_dict(at,fil+'.hdf5')

    return dic