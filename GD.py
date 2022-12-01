import numpy as np
from scipy.fft import fft, ifft
import cmath as cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
import time
import SERV

pi=cm.pi
deg2rad=pi/180

# class gd  -----------------------------------------------

class gd:
    def __init__(self,y,**gdpar): # y ordinate or n
        if isinstance(y,int):
            y=np.zeros(y)
        self.y=y
        if 'ini' in gdpar:
            self.ini=gdpar['ini']
        else:
            self.ini=0
        if 'dx' in gdpar:
            self.dx=gdpar['dx']
        else:
            self.dx=1
        if 'x' in gdpar:
            self.x=gdpar['x']
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
        outgd=copy.copy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd.y=outgd.y+other
        else:
            outgd.y=outgd.y+other.y
        return outgd   

    def __radd__(self,other):
        outgd=copy.copy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd.y=outgd.y+other
        return outgd   

    def __mul__(self,other):
        outgd=copy.copy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd.y=outgd.y*other
        else:
            outgd.y=outgd.y
            other.y
        return outgd   

    def __rmul__(self,other):
        outgd=copy.copy(self)
        if isinstance(other,int):
            other=float(other)
        if isinstance(other,float) or isinstance(other,complex):
            outgd.y=outgd.y*other
        return outgd   
  

def edit_gd(ingd,**gdpar): # 'new'-1  -> new object
    if 'new' in gdpar:
        outgd=copy.copy(ingd)
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
    if ingd.typ == 1 :
        x=ingd.ini+np.arange(ingd.n)*ingd.dx
    else:
        x=ingd.x
    x=SERV.deshape(x)

    return x        


def dict2gd(dicin):
    outgd=gd(dicin['y'],ini=dicin['ini'],dx=dicin['dx'],x=dicin['x'],
    capt=dicin['capt'],cont=dicin['cont'],unc=dicin['unc'],uncx=dicin['uncx'])
    outgd.typ=dicin['type']
    return outgd


def gd2dict(ingd):
    outdic={'ini':ingd.ini,'dx':ingd.dx,'x':ingd.x,'typ':ingd.typ,'y':ingd.y,
    'capt':ingd.capt,'cont':ingd.cont,'unc':ingd.unc,'uncx':ingd.uncx}
    return outdic


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

def set_gd(ingd,fun,par1=1,par2=0.1,par3=0):
    # ingd a gd or an integer that is the length of the new gd
    # fun is the function

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
    # ingd a gd or an integer that is the length of the new gd
    # dist is the distribution
    # par3 seed

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

    if dist == 'cauchy': # par1 mu, par2 scale, par3 seed 
        y=np.random.standard_cauchy(ingd.n)*par2+par1

    outgd.y=y
    return outgd     


# modification function -----------------------------------------------

def modif_gd(ingd,fun,par1=1,par2=0.1,par3=0):
    outgd=copy.copy(ingd)
    y=ingd.y
    if fun == 'abs':
        y=np.abs(y)
    if fun == 'real':
        y=np.real(y)
    if fun == 'imag':
        y=np.imag(y)
    if fun == 'imag':
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

    outgd=ingd
    outgd.capt=outgd.capt+'- rota'
    y=outgd.y

    y=list(y) 

    y = y[n:] + y[:n]

    outgd.y=np.array(y)

    return outgd  


# simple functions

def stat_gd(ingd,nbins=20):
    y=ingd.y
    N=len(y)
    median=np.nanmedian(y)
    mean=np.nanmean(y)
 #   std=np.sqrt(sum((y-mean)**2)/(N-1)) 
    std=np.nanstd(y)
    skew=sum((y-mean)**3)/(N*std**3)
    kurt=sum((y-mean)**4)/(N*std**4)-3

    stat={'mean': mean, 'median': median, 'stdev': std, 'Skewness': skew, 'kurtosis': kurt}
    hist,bin_edge=np.histogram(y,bins=nbins)
    dbin=bin_edge[1]-bin_edge[0]
    Hist=gd(hist,ini=bin_edge[0]+dbin/2,dx=dbin)

    return stat,Hist


def fft_gd(ingd,fif=1):  # fif = -1 ifft
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

def plot_gd(ingd):
    x=x_gd(ingd)
    plt.ion()
    plt.plot(x,ingd.y)
    plt.show() 
    gridon()


def semilogx_gd(ingd):
    x=x_gd(ingd)
    plt.ion()
    plt.semilogx(x,ingd.y)
    plt.show() 
    gridon()


def semilogy_gd(ingd):
    x=x_gd(ingd)
    plt.ion()
    plt.semilogy(x,ingd.y)
    plt.show() 
    gridon()


def loglog_gd(ingd):
    x=x_gd(ingd)
    plt.ion()
    plt.loglog(x,ingd.y)
    plt.show()
    gridon()


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


# calc functions -------------------------------

def minmax_gd(ingd):
    x=x_gd(ingd)
    mima_x=[min(x),max(x)] 
    mima_y=[min(ingd.y),max(ingd.y)]
    
    return mima_x,mima_y 
