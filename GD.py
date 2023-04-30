    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
 .           Module GD

Class and functions for GD management

 Sections:

 > class gd 	                -> dummy_gd
 > gd - dictionary management   -> dummy_gd_dic
 > gd display                   -> dummy_gd_disp
 > set functions                -> dummy_set
 > modification functions       -> dummy_mod
 > simple functions             -> dummy_simp
 > plot functions               -> dummy_plot
 > calc functions               -> dummy_calc
'''

def sections():
    sec=[
        'gd','gd_dic','gd_disp','set','mod','plot','calc'
    ]
    return sec

from scipy.fft import fft, ifft
import cmath as cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
import time
import numpy as np
import BASIC,SERV

pi=cm.pi
deg2rad=pi/180

# class gd  -----------------------------------------------

def dummy_gd():
    '''
    '''
    clas=['gd'
    ]
    fun=[
        'x_gd','div_gd','zero_nan_gd'
    ]

class gd:  
    '''
            gd creation

    gd is the basic class-container in SnagPy.

    The attributes are:

    > y      the ordinate (basic data)
    > n      the length
    > ini    initial abscissa (used in type 1 gds)
    > dx     sampling step (used in type 1 gds)
    > x      abscissas (used in type 2 gds)
    > typ    determine type 1 (virtual abscissa) or type 2 (real abscissa)
    > capt   caption (a string)
    > cont   a control variable (in the case of a bsd, a special structure)
    > unc    ordinate uncertainty (typically not used)
    > uncx   abscissa uncertainty (used for other)

    Any gd can be modified by edit_gd. 
    A gd can be created giving just the number of samples, 
    or transform a 1-dimensional array to a gd, or give the other information. 
    Addition and multiplication are overloaded for gd objects, 
    so you can write, e.g.,  30.5*gd1+gd2*gd3

    Note that gd1 = gd2  is just a different name for gd1, but
    gd1 = gd2 + 0  is a new gd
    Any copy of a gd should be a deepcopy
    '''
    def __init__(self,y,**gdpar): # y ordinate or n
        if isinstance(y,int):
            y=np.zeros(y)
        self.y=y
        if 'ini' in gdpar:
            self.ini=gdpar['ini']
        else:
            self.ini=0.
        if 'dx' in gdpar:
            self.dx=gdpar['dx']
        else:
            dx=1.
            self.dx=dx
            # self.dx=dx.squeeze()
        if 'x' in gdpar:
            self.x=gdpar['x']
            if len(self.x) > 0:
                self.typ=2
        else:
            self.x=[]
            self.typ=1
        if 'y' in gdpar:
            self.y=y
        self.n=len(y) 
        if 'capt' in gdpar:
            self.capt=gdpar['capt']
        else:
           self.capt='' 
        if 'cont' in gdpar:
            self.cont=gdpar['cont']
        else:
            self.cont=0 
        if 'unc' in gdpar:
            self.unc=gdpar['unc']
        else:
            self.unc=0
        if 'uncx' in gdpar:
            self.uncx=gdpar['uncx']
        else:
            self.uncx=0
        self.label='SnagPy gd created on '+time.asctime()

    def __add__(self,other):
        outgd=copy.deepcopy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd.y=outgd.y+other
        else:
            outgd.y=outgd.y+other.y
        return outgd   

    def __radd__(self,other):
        outgd=copy.deepcopy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd.y=outgd.y+other
        return outgd   

    def __mul__(self,other):
        outgd=copy.deepcopy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd.y=outgd.y*other
        else:
            outgd.y=outgd.y*other.y
        return outgd   

    def __rmul__(self,other):
        outgd=copy.deepcopy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd.y=outgd.y*other
        return outgd   
  

    def edit_gd(ingd,**gdpar): 
        """ 
        changes parameters of a gd or copy a gd
        'new'-1  -> new object
        in gdpar any couple varname=varvalue 
        
        Output:
            outgd   edited gd
        """

        if 'new' in gdpar:
            outgd=copy.deepcopy(ingd)
        else:
            outgd=ingd
        if 'x' in gdpar:
            outgd.x=gdpar['x']
            outgd.typ=2
        if 'y' in gdpar:
            outgd.y=gdpar['y']
            outgd.n=len(outgd.y)
        if 'ini' in gdpar:
            outgd.ini=gdpar['ini']
        if 'dx' in gdpar:
            print(gdpar)
            outgd.dx=gdpar["dx"]
        if 'capt' in gdpar:
            outgd.capt=gdpar['capt']
        if 'cont' in gdpar:
            outgd.x=gdpar['cont']
        if 'unc' in gdpar:
            outgd.unc=gdpar['unc']
        if 'uncx' in gdpar:
            outgd.uncx=gdpar['uncx']
        return outgd


def x_gd(ingd):
    '''
    gd abscissa (real or realized) (operates also on arrays)
    
    Output:
        x   gd abscissa
    '''
    if not isinstance(ingd,np.ndarray):
        if ingd.typ == 1 :
            x=ingd.ini+np.arange(ingd.n)*ingd.dx
        else:
            x=ingd.x
        x=BASIC.deshape(x)
    else:
        x=np.arange(len(ingd))

    return x        


def div_gd(gd1,gd2):
    '''
    gd division (operates also on arrays)
    '''
    try:
        y2=gd2.y
        outgd=copy.deepcopy(gd2)
        print('gd2 object')
    except:
        y2=gd2
        print('gd2 numeric')
    try:
        y1=gd1.y
        outgd=copy.deepcopy(gd1)
        print('gd1 object')
    except:
        y1=gd1
        print('gd1 numeric')

    y1=np.array(y1)    # strange, but fundamental
    y2=np.array(y2)

    outgd.y=y1/y2
    outgd.capt='divide'
    return outgd


def zero_nan_gd(ingd,v=0.):
    '''
    substitutes NaN values with 0 or v, for gds or arrays
    '''
    if isinstance(ingd,np.ndarray):
        iy=np.argwhere(np.isnan(ingd))
        iy=BASIC.deshape(iy)
        ingd[iy]=v
        return ingd
    else:
        iy=np.argwhere(np.isnan(ingd.y))
        iy=BASIC.deshape(iy)
        ingd.y[iy]=v
        return ingd


# gd - dictionary management -------------------------

def dummy_gd_dic():
    '''
    '''
    clas=[

    ]
    fun=['dummy_gd_dic','dict2gd','gd2dict',
        
    ]
    return clas,fun

def dict2gd(dicin):
    y=dicin['y']
    if isinstance(y[1],complex):
        print('y is complex')
    else:
        print('y is real') 
    outgd=gd(y,ini=dicin['ini'],dx=dicin['dx'],x=dicin['x'],
    capt=dicin['capt'],cont=dicin['cont'],unc=dicin['unc'],uncx=dicin['uncx'])
    outgd.typ=dicin['type']
    return outgd


def gd2dict(ingd):
    outdic={'ini':ingd.ini,'dx':ingd.dx,'x':ingd.x,'typ':ingd.typ,'y':ingd.y,
    'capt':ingd.capt,'cont':ingd.cont,'unc':ingd.unc,'uncx':ingd.uncx}
    return outdic


# gd display -----------------------------

def dummy_gd_disp():
    '''
    '''
    clas=[

    ]
    fun=['show_gd', 
        
    ]
    return clas,fun

def show_gd(ingd):
    print('type   ',ingd.typ)
    print('n      ',ingd.n)
    print('ini    ',ingd.ini)
    print('dx     ',ingd.dx)
    print('y      ',ingd.y)
    print('capt   ',ingd.capt)
    print('cont   ',ingd.cont)
    if hasattr(ingd,'label'):
        print('  ',ingd.label)



# set functions -----------------------------------------------

def dummy_set():
    '''
    To create gds with certain simple signals or random series
    '''
    clas=[

    ]
    fun=['set_gd','rand_gd',
        
    ]
    return clas,fun

def set_gd(ingd,fun,par1=1,par2=0.1,par3=0):
    '''
    ingd a gd or an integer that is the length of the new gd
    fun is the function and can be:
        'delt','step','ramp','cos','sin','cexp','exp','power','rect'
    par1,pae2,par3  are parameters for fun
    '''

    if isinstance(ingd,int):
        ingd=gd(np.zeros(ingd))

    x=x_gd(ingd)
    outgd=ingd

    if fun == 'delt': # par1 Amp, par2 position 
        Amp=par1
        pos=par2
        n=outgd.n
        y=np.zeros(n)
        y[pos]=1

        y=y*Amp

    if fun == 'step': # par1 Amp, par2 position 
        Amp=par1
        pos=par2
        n=outgd.n
        y=np.zeros(n)
        y[pos:n]=1

        y=y*Amp

    if fun == 'ramp': # par1 first value, par2 step (positive or negative)
        y=np.arange(ingd.n)*par2+par1

    if fun == 'cos': # par1 Amp, par2 fr, par3 ph  
        Amp=par1
        per=1/par2
        ph=par3
        x=(x%per)*2*pi/per+ph*deg2rad
        y=np.cos(x)*Amp
        
    if fun == 'sin': # par1 Amp, par2 fr, par3 ph  
        Amp=par1
        per=1/par2
        ph=par3
        x=(x%per)*2*pi/per+ph*deg2rad
        y=np.sin(x)*Amp
        
    if fun == 'cexp': # par1 Amp, par2 fr, par3 ph  
        Amp=par1
        per=1/par2
        ph=par3
        x=(x%per)*2*pi/per+ph*deg2rad
        y=np.exp(1j*x)*Amp
        
    if fun == 'exp': # par1 Amp, par2 tau 
        Amp=par1
        tau=par2
        x=(x/tau)
        y=np.exp(x)*Amp
    
    if fun == 'power': # par1 Amp, par2 exponent, par3 scale 
        Amp=par1
        expon=par2
        scal=par3
        x=(x*scal)
        y=np.x**expon*Amp
    
    if fun == 'rect': # par1 Amp, par2 length on, par3 length off
                      # if Amp < 0, begins witn off
        Amp=par1
        lon=par2
        loff=par3
        n=outgd.n
        y=np.zeros(n)
        if Amp > 0:
            ion=np.arange(0,n,lon+loff)
            nn=len(ion)
            for i in ion:
                y[i:i+lon]=1
        if Amp < 0:
            Amp=-Amp
            ion=np.arange(0,n,lon+loff)
            nn=len(ion)
            for i in ion:
                y[i+loff:i+lon+loff]=1

        y=y*Amp

    outgd.y=y
    return outgd



def rand_gd(ingd,dist,par1=0,par2=1,par3=0):
    ''' 
    Random numbers

    ingd       a gd or an integer that is the length of the new gd
    dist       is the distribution ('norm','unif',...)
    par1,par2  are the parameters of the distribution
    par3       seed
    '''

    if isinstance(ingd,int):
        ingd=gd(np.zeros(ingd))

    outgd=ingd
    if par3:
        rng=np.random.default_rng(par3)
    else:
        rng=np.random.default_rng()

    if dist == 'unif': # par1 min, par2 max, par3 seed 
        y=np.random.rand(ingd.n)*(par2-par1)+par1

    if dist == 'norm': # par1 mu, par2 sigma, par3 seed 
        y=np.random.randn(ingd.n)*par2+par1

    if dist == 'cnorm': # par1 mu, par2 sigma, par3 seed 
        y=(np.random.randn(ingd.n)+1j*np.random.randn(ingd.n))*par2+par1

    if dist == 'cauchy': # par1 mu, par2 scale, par3 seed 
        y=np.random.standard_cauchy(ingd.n)*par2+par1

    outgd.y=y
    return outgd     


# modification functions -----------------------------------------------

def dummy_mod():
    '''
    '''
    clas=[

    ]
    fun=['modif_gd','rota_gd','resamp_gd','parallel_gd'
        
    ]
    return clas,fun


def modif_gd(ingd,fun,par1=1,par2=0.1,par3=0):
    '''
    gd modification by a function fun

    fun  ex.:'abs','real','imag','angle','log10','xlog10','loglog10'
    '''
    outgd=copy.deepcopy(ingd)
    y=ingd.y
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
        outgd.capt=outgd.capt+' log10'
    if fun == 'xlog10':
        x=x_gd(ingd)
        x=np.imag(x)
        outgd.x=x
        outgd.typ=2
        outgd.capt=outgd.capt+' xlog10'
    if fun == 'loglog10':
        x=x_gd(ingd)
        x=np.imag(x)
        outgd.x=x
        outgd.typ=2
        y=np.log10(y)
        outgd.capt=outgd.capt+' loglog10'

    outgd.y=y
    return outgd     


def rota_gd(ingd,n):
    '''
    circular shift of a gd ordinate
    n is a positive or negative integer
    '''

    outgd=ingd
    outgd.capt=outgd.capt+'- rota'
    y=outgd.y

    y=list(y) 

    y = y[n:] + y[:n]

    outgd.y=np.array(y)

    return outgd  



def resamp_gd(ingd,dx,nmax=0):
    '''
    gd resampling with frequency domain filter
    dx new sampling step
    nmax if > 0, impose fft length
    '''
    N0=ingd.n
    DX0=ingd.dx
    T0=N0*DX0
    DF0=1/T0
    FRMAX0=1/DX0

    if dx == DX0:
        return ingd

    if nmax > 0:
        if nmax < N0:
            print('*** ERROR: nmax ',nmax,' < ',N0,' (gd length)')
            raise RuntimeError
        else:
            y=np.zeros(nmax)
    else:
        y=np.zeros(N0)

    me=np.mean(ingd.y)

    y[0:N0-1]=ingd.y.deepcopy()-me

    if isinstance(y,complex):
        iccompl=1
    else:
        iccompl=0
        if N0 % 2 == 1:
            N0=N0+1
            y=np.append(y,-me)

    enh=DX0/dx
    N1=round(enh*N0/2)*2
    DN=N1-N0

    #  fr0=1/(dx*(2-iccompl))
    y=fft(y)

    Y=np.zeros(N1)

    if iccompl == 0:
        if DN > 0:
            Y[0:N0/2]=y[0:N0/2]
        else:
            Y[0:N1/2]=y[0:N1/2]

        Y[N1:N1/2:-1]=Y[1:N1/2]
    else:
        if DN > 0:
            Y[0:N0]=y
        else:
            Y=y[0:N1]

    Y=ifft(Y)

    out_gd=gd(Y,ini=ingd.ini,dx=ingd.dx,capt=ingd.capt+' resampled')

    return out_gd
    



def parallel_gd(ingd1,ingd2,**gdpar):
    pass


# simple functions

def dummy_simp():
    '''
    '''
    clas=[

    ]
    fun=['stat_gd','fft_gd'
        
    ]
    return clas,fun


def stat_gd(ingd,nbins=20):
    '''
    Simple statistics (parameters and histogram)
    nbins number of bins of the output histogram (a gd). If = 0, no hist
    producs a simple dictionary with median, mean, std, skew and kurt

    ingd    gd or array

    OUTPUT:
        stat,Hist   
        or
        stat
    '''

    if isinstance(ingd,np.ndarray):
        y=ingd
    else:
        y=ingd.y
    N=len(y)
    median=np.nanmedian(y)
    mean=np.nanmean(y)
 #   std=np.sqrt(sum((y-mean)**2)/(N-1)) 
    std=np.nanstd(y)
    skew=sum((y-mean)**3)/(N*std**3)
    kurt=sum((y-mean)**4)/(N*std**4)-3
    mi=min(y)
    ma=max(y)

    stat={'N': N, 'mean': mean, 'median': median, 'stdev': std, 'Skewness': skew, 
          'kurtosis': kurt, 'min': mi, 'max': ma}
    if nbins > 0:
        hist,bin_edge=np.histogram(y,bins=nbins)
        dbin=bin_edge[1]-bin_edge[0]
        Hist=gd(hist,ini=bin_edge[0]+dbin/2,dx=dbin)

        return stat,Hist
    else:
        return stat


def fft_gd(ingd,fif=1):
    '''  
    fft of the ordianate of a gd
    fif = -1 ifft
    '''
    outgd=ingd
    y=ingd.y
    dx=1/(ingd.n*ingd.dx)
    if fif == 1:
        y=fft(y)
        outgd.capt='fft of '+outgd.capt
    else:
        y=ifft(y)
        outgd.capt='ifft of '+outgd.capt

    outgd.ini=0
    outgd.dx=dx
    outgd.y=y

    return outgd


# plot functions -----------------------------------------------

def dummy_plot():
    '''
                    Main plotting procedure
    To plot data (typically a gd, but also for arrays) the main procedure 
    consists in four steps:
    1) newfig      that can define the dimension of the window. If the 
                parameter siz is a two element list, it defines the relative 
                enhancement of the horizontal and vertical dimension of the
                window, if it is just a number it applies to both h and v.
                Example: newfig([1.33,1])
    2) plot_helper defines various aspects of the graph, e.g.:
                mode   : normal plot, steps, scatter plot
                scale  : linear or logatithmic scale for each axis.
                grid   : grid lines
                fmt    : format (color and texture)
                linewid: line width
                It creates a dictionary called P_H.
    3) plot_gd     general plotting function; the imput are:
                ingd  the input gd or array
                P_H   the output of plot_helper; if absent, the default
    4) post_plot   defines title and labels.

    The plot can be modified by the functions xlog, ylog, xlin, ylin,
    xlim, ylim and others.

    ----------------------------------------------------------------

    A different procedure, more synthetic is using simply the GD function "plot".
    Among the arguments you can define the title, the x and y labels, the log or lin, 
    the mode, the color, if the plot is in a new or in the old figure and other.
    Also a label can be attached and can be used in the legend with GD.legend(...).

    '''
    clas=[

    ]
    fun=['newfig','plot_gd','plot_helper','post_plot','c_plot_gd','plot','legend','close_fig'
        
    ]
    return clas,fun


def newfig(siz=1):
    '''
    creates a new figure. siz is a two-element list 
    containing the enhancement of the horizontal and vertical dimensions
    '''    
    if not isinstance(siz,list):
        siz1=siz
        siz=[siz1,siz1]
    hor=6.4*siz[0]
    ver=4.8*siz[1]
    plt.figure(figsize=[hor,ver])
    plt.ion()
    plt.grid(True)
    plt.show()



def plot_gd(ingd,x=[],P_H=0):
    '''
    ingd   input gd or array
    x      changed x
    P_H    output of plot_helper (can be defaulted)
    '''

    if len(x) > 0:
        xx=x
        chx=1
    else:
        chx=0

    if isinstance(P_H,int):
        print("default P_H")
        P_H=plot_helper()
    if isinstance(ingd,np.ndarray):
        x=np.arange(len(ingd))
        y=ingd
    else:
        x=x_gd(ingd)
        y=ingd.y

    if chx == 1:
        x=xx

    sca=P_H['scale']
    if sca[0:2] == 'lo':
        plt.xscale('log') 
    if sca[2:4] == 'lo':
        plt.yscale('log') 
    if sca[0:2] == 'li':
        plt.xscale('linear') 
    if sca[2:4] == 'li':
        plt.yscale('linear')

    if P_H['mode'] == 'plot':      
        plt.plot(x,y)
    elif P_H['mode'] == 'step': 
        plt.step(x,y,where='mid')
    elif P_H['mode'] == 'scat': 
        plt.scatter(x,y)

    gridon()



def plot_helper(mode='plot',scale='lili',grid='norm',absc='norm',fmt='b',linewid=1.):
    '''
    auxiliary input for plot_gd
    mode    'plot','step','scat'
    scale   'lili','lilo','loli','lolo'
    grid    'norm','dens','no'
    absc    'norm','min','hour','days','week','mjd','date'
    fmt     
    linewid
    '''

    P_H={'mode':mode,'scale':scale,'grid':grid,'absc':absc,fmt:'fmt',
    'linewid':linewid,'marker':'o','marksiz':10}

    return P_H



def post_plot(tit,xlab,ylab):
    '''
    inserts title and h and v labels
    '''
    plt.title(tit)
    plt.xlabel(xlab)
    plt.ylabel(ylab)



def c_plot_gd(typ,ingd):
    ''' 
    plot for complex array
    typ = 0 i vs r, 1 real, 2 imag, 3 abs, 4 angle
    '''
    if isinstance(ingd,np.ndarray):
        x=np.arange(len(ingd))
        y=ingd
    else:
        x=x_gd(ingd)
        y=ingd.y
    plt.ion()
    if typ == 0:
        x=np.real(y)
        y=np.imag(y)
        plt.plot(x,y)
    if typ == 1:
        y=np.real(y)
        plt.plot(x,y)
    if typ == 2:
        y=np.imag(y)
        plt.plot(x,y)
    if typ == 3:
        y=np.abs(y)
        plt.plot(x,y)
    if typ == 4:
        y=SERV.atan3(y)
        plt.plot(x,y)
 #   plt.show() 
    plt.grid()


def plot(ingd,**plotpar):
    '''
    ingd   input gd or array
    x      changed x

    keywords:
        more=1          no newfig
        mode='step'
        mode='scatter'
        xlog=1
        ylog=1
        loglog=1
        plmodif         plot modifier (es.: 'r+', 
                        see https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html)
        x=array         abscissa modification
        fmt='r','b',... (for step)

        title='Something'  put a title
        xlab='Something'   put a xlabel
        ylab='Something'   put a ylabel

        label='Something'  sets a label to the line
    '''

    if 'more' not in plotpar:
        print('New Fig')
        newfig()
    else:
        print('No New Fig')

    mode='plot'
    if 'mode' in plotpar:
        mode=plotpar['mode']
        
    if 'plmodif' in plotpar:
        plmodif=plotpar['plmodif']
    else:
        plmodif=''
        
    if 'fmt' in plotpar:
        fmt=plotpar['fmt']
    else:
        fmt=''

    plt.xscale('linear')
    plt.yscale('linear')
    if 'xlog' in plotpar:
        plt.xscale('log')
    if 'ylog' in plotpar:
        plt.yscale('log')
    if 'loglog' in plotpar:
        plt.xscale('log')
        plt.yscale('log')

    if 'x' in plotpar:
        xx=plotpar['x']
        chx=1
    else:
        chx=0

    if isinstance(ingd,np.ndarray):
        x=np.arange(len(ingd))
        y=ingd
    else:
        x=x_gd(ingd)
        y=ingd.y

    if chx == 1:
        x=xx

    if 'label' in plotpar:
        lab=plotpar['label']
    else:
        lab=''

    if mode == 'plot':      
        plt.plot(x,y,plmodif,label=lab)
    elif mode == 'step': 
        plt.step(x,y,fmt,where='mid',label=lab)
    elif mode == 'scat': 
        plt.scatter(x,y,label=lab)

    gridon()

    if 'title' in plotpar:
        tit=plotpar['title']
        plt.title(tit)
      
    if 'xlab' in plotpar:
        xlab=plotpar['xlab']
        plt.xlabel(xlab)
      
    if 'ylab' in plotpar:
        ylab=plotpar['ylab']
        plt.ylabel(ylab)
      

def legend(loc='lower right'):
    '''
    Legend creation using labels set by plot
    For the location, use 'best' or any position; default 'lower right'
    see https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html#matplotlib.pyplot.legend
    '''
    plt.legend(loc=loc)


def close_fig():
    plt.close()


def gridon():
    plt.grid(True)


def xlog():
    plt.xscale('log')


def ylog():
    plt.yscale('log')


def xlin():
    plt.xscale('linear')


def ylin():
    plt.yscale('linear')


def xlim(mi,ma):
    plt.xlim([mi,ma])


def ylim(mi,ma):
    plt.ylim([mi,ma])


def fig_limits():
    l,r=plt.xlim()
    xx=[l,r]
    d,u=plt.ylim()
    yy=[d,u]

    return xx,yy



def ioff():
    '''
    interactive off (use plt.show() to show the figure)
    '''
    plt.ioff()    

def ion():
# interactive on 
    plt.ion()  
    plt.show()


def holdoff():
    plt.clf()



# calc functions -------------------------------

def dummy_calc():
    '''
    '''
    clas=[

    ]
    fun=['minmax_gd'
        
    ]
    return clas,fun

def minmax_gd(ingd):
    '''
    Min and max for abscissa and ordinate
    '''
    x=x_gd(ingd)
    mima_x=[min(x),max(x)] 
    mima_y=[min(ingd.y),max(ingd.y)]
    
    return mima_x,mima_y 
