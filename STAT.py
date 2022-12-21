import numpy as np
from scipy.fft import fft, ifft
from scipy import stats
from scipy import signal
import scipy as sc
import matplotlib.pyplot as plt
import GD 


# Zero management ----------------------------

def ana_zero(indat,mode=0,eps=1.e-6):   
   # a zero period is at least 2 samples
   # mode 0 only percentage of zeros
   #      1 zero function (1 for zeros, 0 elsewhere)
   #      2 zero intervals
   #      3 one  intervals
   # eps to the lowest meaningful value, as percentage of maximum

   if isinstance(indat,gd):
      y=indat.y
      ty='gd'
   elif isinstance(indat,np.ndarray):
      y=indat
      ty='array'
   else:
      y=np.array(indat)
      ty='other'

   n=len(y)
   zer=np.zeros(n)
   yy=abs(y)
   ma=max(yy)
   eps1=eps*ma

   for i in range(n-1):
 #     if y[i] or y[i+1] == 0:
      if yy[i] <= eps1 and yy[i+1] <= eps1:
         zer[i]=1

   if abs(y[n-1]) <= 0:
      zer[n-1]=1

   perc=sum(zer)/n

   if mode == 0:
      return perc
   if mode == 1:
      if ty == 'gd':
         zer=edit_gd(indat,y=zer)
      else:
         zer=gd(zer)
      return perc,zer
   if mode == 2:
      di=np.diff(zer)
      N=int(np.ceil(np.sum(abs(di))/2))
      inds=np.zeros((2,N))
      ini=np.argwhere(di == 1)
      fin=np.argwhere(di == -1)
      N=len(ini)
      inds[0]=ini.reshape(len(ini))
      inds[1]=fin.reshape(len(fin))
      if len(ini) > len(fin):
         inds[1][N-1]=n-1
      if ty == 'gd':
         xzer=indat.ini+indat.dx*inds
         return xzer
      else:
         return inds
   if mode == 3:
      di=np.diff(1-zer)
      N=int(np.ceil(np.sum(abs(di))/2))
      inds=np.zeros((2,N))
      ini=np.argwhere(di == 1)
      fin=np.argwhere(di == -1)
      N=len(ini)
      inds[0]=ini.reshape(len(ini))
      inds[1]=fin.reshape(len(fin))
      if len(ini) > len(fin):
         inds[1][N-1]=n-1
      if ty == 'gd':
         xone=indat.ini+indat.dx*inds
         return xone
      else:
         return inds
      
   return perc,zer


# def zeroing(ingd,**zpar):
#    if 'xinter' in zpar:
#       pass
#    if 'inter' in zpar:
#       pass
#    if 'zgd' in zpar:
#       pass




# Histograms ----------------------------





# Power Spectra ----------------------------

def gd_pows(ingd,npiece=1,res=1,shift=1,nobias=1,notrend=1,window=2,singleb=1,center=1):
# Standard power spectrum estimation
# npiece    number of pieces (without shift)
# res       resoltion (minimal 1)
# shift     for interlaced pieces (1 no interlace, 0.5 one half shift)
# nobias    =1 subtract bias
# notrend   =1 subtract trend
# window    -0 no, 1 bartlett, 2 hann (default), 3 flat top cosine
# singleb   =1 single side band (for real data, normalized for single band)
# center    =1 0 frequency at center, no single side band

   N=ingd.n
   dx=ingd.dx
   y=ingd.y
   if isinstance(y[1],complex):
      singleb=0

   if npiece < 1:
      print(' *** Minimal number of pieces is 1')
      npiece=1

   if res < 1:
      print(' *** Minimal resolution is 1')
      res=1

   if shift > 1:
      print(' *** Maximal shift is 1')
      shift=1

   mea=0
   if nobias == 1:
      mea=np.mean(y)
      y=y-mea

   if notrend == 1:
      y=signal.detrend(y)

   n=int(N/npiece)
   sh=int(shift*n)
   
   win=np.ones(n)
   if window == 1:
      win=signal.windows.bartlett(n)
      win=win*np.sqrt(n/sum(win**2))
   if window == 2:
      win=signal.windows.hann(n)
      win=win*np.sqrt(n/sum(win**2))
   if window == 3:
      win=signal.windows.tukey(n,alpha=0.5)
      win=win*np.sqrt(n/sum(win**2))

   rang=range(0,N-n+1,sh)
   lS=int(np.ceil(res*n/2.)*2)
   lS2=int(lS/(1+singleb))
   # print(n,res)
   # print(lS,lS2)
   S=np.zeros(lS)
   df=1/(lS*dx)

   npie=0

   for i in rang:
      i1=i
      i2=i1+n
      npie+=1
      yy=y[i1:i2]*win
      YY=np.zeros(lS)
      if isinstance(yy[1],complex):
         YY=YY+1j*YY
      YY[0:n]=yy
      YY=np.abs(fft(YY))**2
      print(np.shape(S))
      print(np.shape(YY))
      S+=YY[0:lS2]

   S=S*dx/npie
   print('N,n,sh,npie,lS',N,n,sh,npie,lS)
   S=GD.gd(S,dx=df)

   return S


def gd_welch(ingd,lenfrac,res=1,shift=1,notrend=1,win='hann',singleb=1):
   y=ingd.y
   dx=ingd.dx
   fs=1/dx
   n=ingd.n
   nperseg=int(n*lenfrac/2)*2
   if shift < 1:
      noverlap=int(nperseg*shift)
   else:
      noverlap=None
   nfft=int(nperseg*res/2)*2
   f,s=signal.welch(y,fs=fs,nperseg=nperseg,noverlap=noverlap,
   nfft=nfft,window=win)
   s=GD.gd(s,dx=f[1])

   return s

   
def cross_pw(ingd1,ingd2):
   pass

def coher(ingd1,ingd2):
   pass

def lombscargle(ingd):
   pass

def pulse_spectrum(ingd):
   pass

def stft(ingd):
   pass



# Spectrograms ----------------------------

def spectrogram(ingd):
   pass




# Period analysis ------------------------

def gd_period():
   pass


def gd_worm():
   pass
