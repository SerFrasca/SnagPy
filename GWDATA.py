    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

import numpy as np
import ASTROTIME,BASIC

Pi=np.pi

# CW sources  ---------------------------------------

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


