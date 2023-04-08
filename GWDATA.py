    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
        Module GWDATA
Basic GW data management

Sections:

> General       -> dummy_gen
> GW data       -> dummy_data
> Doppler data  -> dummy_doppl
> CW sources    -> dummy_cw
> Antennas      -> dummy_ant
> 5-vect        -> dummy_5vec
> Other         -> dummy_oth
'''

import numpy as np
import ASTROTIME,BASIC
import os
import h5py

Pi=np.pi

# General -----------------------------------------

def dummy_gen():
    '''
    '''

def set_symbols(snagpy_p=0):
    '''
    Symbols to folders, files and data
    snagpypath   the path to snagpy, as a string (es.: 'D:\\OneDrive\\SF\\_Prog\\Python\\SnagPy')
                if absent, taken by the environment variable SNAGPY_PATH
    '''

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



# GW data --------------------------------------

def dummy_data():
    '''
    '''

def explore_gw_hdf5(fil):
    '''
    explores hdf5 frame data
    '''
    dic=BASIC.read_hdf52dic(fil)
    # Keys,Dics=BASIC.expl_dict(dic)

    meta=dic['meta']
    print('                meta')
    print('Description    : ',str(meta['Description']))
    print('DescriptionURL : ',str(meta['DescriptionURL']))
    print('Detector       : ',str(meta['Detector']))
    print('Duration       : ',meta['Duration'])
    print('FrameType      : ',str(meta['FrameType']))
    print('GPSstart       : ',meta['Description'])
    print('Observatory    : ',str(meta['Observatory']))
    print('StrainChannel  : ',str(meta['StrainChannel']))
    print('Type           : ',str(meta['Type']))
    print('UTCstart       : ',str(meta['UTCstart']))

    strai=dic['strain']
    l=list(strai.keys())
    print('           ')
    print('                strain')
    for it in strai:
        print(it,' <-> ',strai[it])

    qual=dic['quality']

    return meta,qual



def read_gw_hdf5(fil):
    '''
    reads hdf5 frame data
    '''
    dic=BASIC.read_hdf52dic(fil)
    meta=dic['meta']
    qual=dic['quality']
    strai=dic['strain']
    strain=strai['Strain']

    return meta,strain,qual



# Doppler data --------------------------------

def dummy_doppl():
    '''
    '''

def extr_doppler(tab,tin,tfi,table_par):
    '''
    extract data from doppler tables (in hdf5 format)

     tab        Doppler table (ex.: table_virgo_p)
     tin,tfi    times (vect as [a,m,d,h,m,s])
     table_par  as defined in the symbols or starting.py
    '''

    gtin=ASTROTIME.now(tin,form='gps')
    gtfi=ASTROTIME.now(tfi,form='gps')

    ini=int((gtin-table_par[0])/table_par[2])
    ifi=int((gtfi-table_par[0])/table_par[2])+1

    out=BASIC.read_hdf5_part(tab,'array','['+str(ini)+':'+str(ifi)+']')
    t=out[:,0]
    tmjd=ASTROTIME.t_conv(t,'gps','mjd')

    return out,tmjd


# CW sources ---------------------------------------

def dummy_cw():
    '''
    '''

def ligo2virgo_cw_table(ligotab,capttab='',run=''):
    '''
    creates a Virgo format cw list from a Ligo table
    '''
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

def dummy_ant():
    '''
    '''

def anten(ant,anten_tab=0):
    '''
    extract data for an antenna (a dictionary)
     ant         antenna name ('virgo','ligol',...)
     anten_tab   antenna table (def possible if the environment variable is set)
    '''
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
    '''
    creates dictionarys for each antenna
    typical call: virgo,ligol,ligoh,kagra=GWDATA.antennas()
    '''
    virgo=anten('virgo')
    ligol=anten('ligol')
    ligoh=anten('ligoh')
    kagra=anten('kagra')

    return virgo,ligol,ligoh,kagra



# 5-vect --------------------------------------

def dummy_5vec():
    '''
    '''




# Other ---------------------------------------

def dummy_oth():
    '''
    '''

def conv_jpl_data(datin,noline,nr,items):
    '''
    converts text data (typically table_xxxx.dat) to hdf5 format
        datin    complete path
        noline   number of lenes to jump (typ. 4)
        nr       number of output rows   (typ. 526032)
        items    column to output (e.g. [0,1,3,7], typ. [1,2,3,4,5,6,7,8])
    '''

    sep=os.sep

    path,fil,ext=BASIC.path_fil_ext(datin)

    arr=BASIC.text2array(path+sep+fil+ext,noline,nr,items)
    at=BASIC.array_table(['gpst','x','y','z','vx','vy','vz','einst'],arr)
    dic=BASIC.array_table_to_dict(at,fil+'.hdf5')

    return dic