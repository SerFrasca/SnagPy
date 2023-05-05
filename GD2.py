    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
'''
        Module GD2

Class and functions for GD2 management

 Sections:

 > class gd2                    -> dummy_gd2
 > gd2 - dictionary management  -> dummy_dict
 > gd2 display                  -> dummy_disp
 > modification function        -> dummy_mod
 > map functions                -> dummy_map
'''

def sections():
    sec=[
    ]
    return sec

import numpy as np
from scipy.fft import fft2, ifft2
import cmath
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import copy
import time
import SERV,BASIC,GD

pi=cmath.pi
deg2rad=pi/180

# class gd2  -----------------------------------------------

def dummy_gd2():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

class gd2:    # gd2 creation
    """  
    gd2 is the class-container for 2-D data in SnagPy.
    The horizontal abscissa is the principal one and can be real or virtual; it
    has a dimension n/m, the number of columns. Typically is the time.
    The secondary abscissa is the vertical one, so it is the first of the array;
    it is always virtual and has dimension m, the namber of rows. Typically it
    is the frequency.
    The attributes are:
    > y      the 2-dim ordinate (basic data)
    > n      total number of values
    > m      secondary abscissa dimension
    > ini    initial abscissa (used in type 1 gd2s)
    > dx     sampling step (used in type 1 gd2s)
    > ini2
    > dx2
    > x      abscissas (used in type 2 gds)
    > typ    determine type 1 (virtual abscissa) or type 2 (real abscissa)
    > capt   caption (a string)
    > cont   a control variable (in the case of a bsd, a special structure)
    > unc    ordinate uncertainty (typically not used)
    > uncx   abscissa uncertainty (used for other)
    

    Any gd2 can be modified by edit_gd2. 
    A gd2 can be created giving just the m and n parameters, 
    or transform a 2-dimensional array to a gd2, or give the other information. 
    Addition and multiplication are overloaded for gd2 objects.
    """
    def __init__(self,y,**gdpar): # y ordinate or n and m
        if isinstance(y,int):
            y=np.zeros(y)
        self.y=y
        if 'ini' in gdpar:
            self.ini=gdpar['ini']
        else:
            self.ini=0
        if 'ini2' in gdpar:
            self.ini2=gdpar['ini2']
        else:
            self.ini2=0
        if 'dx' in gdpar:
            self.dx=gdpar['dx']
        else:
            self.dx=1
        if 'dx2' in gdpar:
            self.dx2=gdpar['dx2']
        else:
            self.dx2=1
        if 'x' in gdpar:
            self.x=gdpar['x']
            if len(self.x) > 0:
                self.typ=2
        else:
            self.x=[]
            self.typ=1
        if 'y' in gdpar:
            self.y=y
        A=y.shape
        self.m=A[0]
        self.n=A[0]*A[1] 
        if 'capt' in gdpar:
            self.capt=gdpar['capt']
        else:
           self.capt='' 
        if 'cont' in gdpar:
            self.cont=gdpar['cont']
        else:
            self.cont=0 
        if 'm' in gdpar:
            self.m=gdpar['m']
        if 'ini2' in gdpar:
            self.ini2=gdpar['ini2']
        if 'dx2' in gdpar:
            self.dx2=gdpar['dx2']
        self.label='SnagPy gd2 created on '+time.asctime()

    def __add__(self,other):
        outgd2=copy.deepcopy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd2.y=outgd2.y+other
        else:
            outgd2.y=outgd2.y+other.y
        return outgd2   

    def __radd__(self,other):
        outgd2=copy.deepcopy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd2.y=outgd2.y+other
        return outgd2   

    def __mul__(self,other):
        outgd2=copy.deepcopy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd2.y=outgd2.y*other
        else:
            outgd2.y=outgd2.y
            other.y
        return outgd2  

    def __rmul__(self,other):
        outgd2=copy.deepcopy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd2.y=outgd2.y*other
        return outgd2   
  

def edit_gd2(ingd2,**gdpar): # 'new'-1  -> new object
    if 'new' in gdpar:
        outgd2=copy.deepcopy(ingd2)
    else:
        outgd2=ingd2
    if 'x' in gdpar:
        outgd2.x=gdpar['x']
        outgd2.typ=2
    if 'y' in gdpar:
        outgd2.y=gdpar['y']
        outgd2.n=len(outgd2.y)
    if 'ini' in gdpar:
        outgd2.ini=gdpar['ini']
    if 'dx' in gdpar:
        print(gdpar)
        outgd2.dx=gdpar["dx"]
    if 'capt' in gdpar:
        outgd2.capt=gdpar['capt']
    if 'cont' in gdpar:
        outgd2.x=gdpar['cont']
    if 'unc' in gdpar:
        outgd2.unc=gdpar['unc']
    if 'uncx' in gdpar:
        outgd2.uncx=gdpar['uncx']
    return outgd2


def x_gd2(ingd2):
    '''
    gd2 main abscissa (real or realized) (operates also on arrays)
    '''

    if not isinstance(ingd2,np.ndarray):
        if ingd2.typ == 1 :
            x=ingd2.ini+np.arange(ingd2.n/ingd2.m)*ingd2.dx
        else:
            x=ingd2.x
        x=BASIC.deshape(x)
    else:
        A=ingd2.shape
        x=np.arange(A[0])

    return x


def x2_gd2(ingd2):
    '''
    gd2 secondary abscissa (real or realized) (operates also on arrays)
    '''
    if not isinstance(ingd2,np.ndarray):
        x2=ingd2.ini2+np.arange(ingd2.m)*ingd2.dx2
        x2=BASIC.deshape(x2)
    else:
        A=ingd2.shape
        x2=np.arange(A(1))

    return x2       


def zero_nan_gd2(ingd2,v=0.):
    pass


# gd2 - dictionary management ------------------------

def dummy_dict():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def dict2gd2(dicin):
    pass


def gd22dict(ingd):
    pass


# gd2 display -----------------------------

def dummy_disp():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def show_gd2(ingd2):
    print('type   ',ingd2.typ)
    print('n      ',ingd2.n)
    print('m      ',ingd2.m)
    print('n/m    ',ingd2.n/ingd2.m)
    print('ini    ',ingd2.ini)
    print('dx     ',ingd2.dx)
    print('ini2   ',ingd2.ini2)
    print('dx2    ',ingd2.dx2)
    print('y      ',ingd2.y)
    print('capt   ',ingd2.capt)
    print('cont   ',ingd2.cont)
    if hasattr(ingd2,'label'):
        print('  ',ingd2.label)


# modification function -----------------------------------------------

def dummy_mod():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def modif_gd2(ingd2,fun,par1=1,par2=0.1,par3=0):
    '''
    gd modification by a function fun

    fun  ex.:'abs','real','imag','angle','log10','xlog10','loglog10'
    '''
    outgd2=copy.deepcopy(ingd2)
    y=ingd2.y
    if fun == 'abs':
        y=np.abs(y)
    if fun == 'real':
        y=np.real(y)
    if fun == 'imag':
        y=np.imag(y)
    if fun == 'angle':
        y=np.angle(y,deg=True)
    if fun == 'log10':
        y=np.log10(y)
        outgd2.capt=outgd2.capt+' log10'
    # if fun == 'xlog10':
    #     x=x_gd2(ingd2)
    #     x=np.imag(x)
    #     outgd2.x=x
    #     outgd2.typ=2
    #     outgd2.capt=outgd2.capt+' xlog10'
    # if fun == 'loglog10':
    #     x=x_gd2(ingd2)
    #     x=np.imag(x)
    #     outgd2.x=x
    #     outgd2.typ=2
    #     y=np.log10(y)
    #     outgd2.capt=outgd2.capt+' loglog10'

    outgd2.y=y
    return outgd2     


def rota_gd2(ingd2,n):
    pass



def fft_gd2(ingd2,fif=1):  
    '''
    fft2 of the ordinate of a gd2
    fif = -1 ifft
    '''
    outgd2=ingd2
    y=ingd2.y
    dx=1/(ingd2.n*ingd2.dx)
    if fif == 1:
        y=fft2(y)
        outgd2.capt='fft of '+outgd2.capt
    else:
        y=ifft2(y)
        outgd2.capt='ifft of '+outgd2.capt

    outgd2.ini=0
    outgd2.dx=dx
    outgd2.y=y

    return outgd2


def im_cut(ingd2,iout,jout,typ=1):
    '''
    image cut
    ingd2   input gd2 or matrix
    iout    range x1 out
    jout    range x2 out
    typ     1 range in samples, 2 range in x
    '''
    if isinstance(ingd2,np.ndarray):
        ini1=0
        dx1=1
        ini2=0
        dx2=1
        yy=ingd2
        icgd2=0
    else:
        ini1=ingd2.ini
        dx1=ingd2.dx
        ini2=ingd2.ini2
        dx2=ingd2.dx2
        yy=ingd2.y
        icgd2=1

    # print(iout[0],iout[1],ini1,dx1)
    if typ == 1:
        iiout0=iout[0]
        iiout1=iout[1]
        jjout0=jout[0]
        jjout1=jout[1]
    else:
        iiout0=BASIC.ind_from_inidx(iout[0],ini1,dx1)
        iiout1=BASIC.ind_from_inidx(iout[1],ini1,dx1)
        jjout0=BASIC.ind_from_inidx(jout[0],ini2,dx2)
        jjout1=BASIC.ind_from_inidx(jout[1],ini2,dx2)

    ini=ingd2.ini+(iiout0-1)*ingd2.dx
    ini2=ingd2.ini2+(jjout0-1)*ingd2.dx
    maxy,maxx=yy.shape

    if iiout1 > maxx:
        iiout1=maxx
    if jjout1 > maxy:
        jjout1=maxy

    print(iout,jout)
    print(iiout0,iiout1,jjout0,jjout1)
    # cut=yy[np.ix_(range(iiout0,iiout1),range(jjout0,jjout1))]
    cut=yy[jjout0:jjout1,iiout0:iiout1]
    print(cut.shape)

    if icgd2 == 1:
        cut=gd2(cut,dx=ingd2.dx,ini=ini,dx2=ingd2.dx2,ini2=ini2)

    return cut




def im_reduce(ingd2,nr,nc,typ):
    '''
    image reduction
    '''
    pass


def gd2_stat_nz(ingd2):
    if isinstance(ingd2,np.ndarray):
        y=ingd2
    else:
        y=ingd2.y
    aa=y.shape
    ny=aa[0]*aa[1]
    mask=np.nzero(y)
    yy=y(mask).flatten
    nyy=len(yy)
    stat=GD.stat_gd(yy,nbins=100)
    st1={'tot_el':ny,'nzero_el':nyy}
    stat_nz=dict(st1,**stat)

    return stat_nz


# map functions -----------------------------------------------

def dummy_map():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def newfig2(siz=1):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    return fig,ax


def grey_map(ingd2,fun='none',cmap='cool',modif=0,modpar=1,vmin=0,vmax=0,
             MH=0,gridlin=0.5,gridstyl='--',gridcol='y',norm='linear',alpha=1,
             figsize=1,colorbar=1):
    ''' 
    grey map of 2-D array or gd2

    ingd2       a gd2 or 2D array
    fun         function applied to ingd2 ('abs', 'log', 'sqrt')
    cmap        colormap ('viridis','plasma','cool','hot')
                 see https://matplotlib.org/stable/tutorials/colors/colormaps.html
    modif       cmap modification (0 no (def), 1,...)
    vmin,vmax   range of grey to be considered 
                (if vmin=vmax not considered)
    MH          helper dict
    gridlin,gridstyl,gridcol    grid control
    norm        scale "linear", "log", "symlog", "logit" 
    alpha       transparency (0->1)
    '''

    if isinstance(MH,dict):
        if 'gridlin' in MH:
            gridlin=MH['gridlin']
        if 'gridstyl' in MH:
            gridstyl=MH['gridstyl']
        if 'gridcol' in MH:
            gridcol=MH['gridcol']
        if 'norm' in MH:
            norm=MH['norm']
        if 'alpha' in MH:
            alpha=MH['alpha']

    if isinstance(ingd2,np.ndarray):
        y=ingd2
        aa=y.shape
        inix=0
        finx=aa[1]-1
        iniy=0
        finy=aa[0]-1
    else:
        y=ingd2.y
        inix=ingd2.ini
        finx=inix+(ingd2.n/ingd2.m-1)*ingd2.dx
        iniy=ingd2.ini2
        finy=iniy+(ingd2.m-1)*ingd2.dx2
        
    if fun == 'abs':
        y=abs(y)   
    if fun == 'log':
        y=np.log10(abs(y))  
    if fun == 'sqrt':
        y=np.sqrt(abs(y))
    ext=(inix,finx,iniy,finy)
    
    siz=figsize
    if not isinstance(siz,list):
        siz1=siz
        siz=[siz1,siz1]
    hor=6.4*siz[0]
    ver=4.8*siz[1]
    if colorbar == 1:
        hor=hor*1.2
    fig = plt.figure(figsize=[hor,ver])
    # fig = plt.figure()
    plt.ion()
    ax = fig.add_subplot(111)
    if modif != 0:
        cmap=modif_colormap(cmap,typ=modif,par=modpar)
    if vmin == vmax:
        im = ax.imshow(y, interpolation='none', aspect='auto', origin='lower',
            cmap=cmap,alpha=alpha,extent=ext)
    else:
        im = ax.imshow(y, interpolation='none',aspect='auto',origin='lower',
            cmap=cmap,alpha=alpha,extent=ext,vmin=vmin,vmax=vmax)
    if colorbar == 1:
        plt.colorbar(im)

    ax.grid(which='major', color=gridcol, linestyle=gridstyl, linewidth=gridlin)
    plt.show()

    return im,ax


def map_helper(cmap='cool',norm='linear',alpha=1,
        gridlin=0.5,gridstyl='--',gridcol='y'):
    """ 
    map helper 

        norm   scale "linear", "log", "symlog", "logit" 
        vmin
        vmax
        alpha   transparency (0->1) 

    OUTPUT:
        MH      a dictionary
    """

    MH={'cmap':cmap,'norm':norm,'alpha':alpha,
        'gridlin':gridlin,'gridstyl':gridstyl,'gridcol':gridcol}

    return MH


def post_map(tit,xlab,ylab):
    '''
    inserts title and h and v labels
    '''
    plt.title(tit)
    plt.xlabel(xlab)
    plt.ylabel(ylab)


def scat_gd2(ingd2,fun):
    pass


def post_fig2(tit,xlab,ylab):
# inserts title and h and v labels
    plt.title(tit)
    plt.xlabel(xlab)
    plt.ylabel(ylab)


def modif_colormap(cmap,typ=1,par=1):
    '''
    Creates a modified color map.

    cmap    color map (ex: 'viridis')
    typ     modification type
            = 1 black for lower values
    par     type parameters
    '''

    nres=256
    ccmap=mpl.colormaps[cmap].resampled(nres)
    newcolors = ccmap(np.linspace(0, 1, nres))
    black=np.array([0.1, 0.1, 0.1, 0.5])
    white=np.array([1, 1, 1, 1])
    newcolors[:par,:]=black
    newcolors[nres-par:nres,:]=white
    modcmap=ListedColormap(newcolors,name='cm_modif')

    return modcmap


def list_cmaps():
    cmaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
    BASIC.show_list(cmaps)