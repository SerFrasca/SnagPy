   # Copyright (C) 2023  Sergio Frasca, Riccardo Felicetti
   #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

import numpy as np
from scipy.fft import fft, ifft
from scipy import stats
from scipy import signal
import scipy as sc
import matplotlib.pyplot as plt
import GD,BASIC,ASTROTIME,GD2


# Zero management ----------------------------

def ana_zero(indat,mode=0,eps=1.e-6):   
   # a zero period is at least 2 samples
   # mode 0 only percentage of zeros
   #      1 zero function (1 for zeros, 0 elsewhere)
   #      2 zero intervals
   #      3 one  intervals
   # eps to the lowest meaningful value, as percentage of maximum
   #    = 0 no enlarged 0 output

   if isinstance(indat,GD.gd):
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

   if eps > 0:
      y1=y*zer

   perc=sum(zer)/n

   if mode == 0:
      return perc
   if mode == 1:
      if ty == 'gd':
         zer=GD.edit_gd(indat,y=zer)
      else:
         zer=GD.gd(zer)
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
      
   if eps == 0:
      return perc,zer
   else:
      if ty == 'gd':
         outdat=indat+0
         outdat.y=y1
      else:
         outdat=y1
      return perc,zer,outdat


# def zeroing(ingd,**zpar):
#    if 'xinter' in zpar:
#       pass
#    if 'inter' in zpar:
#       pass
#    if 'zgd' in zpar:
#       pass


def stat_nozero(indat):
# statistics on non-zero data

   if isinstance(indat,GD.gd):
      indat=indat.y
   if isinstance(indat,GD2.gd2):
      indat=indat.y
      
   indat=indat.flatten()

   out=np.where(indat != 0)
   out=out[0]

   stat=GD.stat_gd(out,nbins=0)

   return stat
   


# Histograms and Parameters ----------------------------

def CW_histogram(x,w=[],ini=[],step=1,n=[],enl=10,typ='tri',edg=0,verb=1):
# Convolutional weighted histogram
#  x     values
#  w     wheights
#  ini   initial value (def the minimum)
#  step  histogram step
#  n     number of bin (def such to cover all the x range)
#  enl   enlargement factor
#  typ   'tri' triangular, 'gau' gaussian, 'rec' rectangular, if array -> personal fun
#  edg   = 1, out values in the edge bins
#  verb  verbosity
   if verb > 1:
      BASIC.tic()

   x=np.array(x)
   N=len(x)
   dic={'N':N}
   mea=np.mean(x)
   med=np.median(x)
   std=np.std(x)
   ske=stats.skew(x)
   kur=stats.kurtosis(x)
   dic['xmean']=mea
   dic['xmedian']=med
   dic['xstd']=std
   dic['xskew']=ske
   dic['xkurt']=kur

   w=np.array(w)

   if ini == []:
      ini=min(x)
   if n == []:
      ma=max(x)
      n=np.ceil((ma-ini)/step)

   if len(w) == 0:
      w=np.ones(len(x))
   else:
      M=sum(w)
      wmea=np.mean(w)
      wmed=np.median(w)
      wstd=np.std(w)
      dic['M']=M
      dic['wmean']=wmea
      dic['wmedian']=wmed
      dic['wstd']=wstd

   sstep=step/enl
   N=round(n*enl)
   base=np.zeros(N)
   y=np.round((x-ini)/sstep)
   # print(n,N,ini,sstep)
   # print(len(w))

   for i in range(len(x)):
      # print(x[i])
      yi=int(y[i])
      if yi >=0 and yi < N:
         # print(i,yi)
         base[yi]=base[yi]+w[i]

   nenl=np.round(enl)
   nenl2=2*nenl-1
   # print(nenl,nenl2)

   if typ == 'tri':
      fun=np.arange(nenl2)
      for i in range(1,nenl):
         # print(i,i+nenl,fun[i+nenl],nenl2-fun[i+nenl])
         fun[i+nenl-1]=fun[nenl-i-1]
      # print(fun)
      rit=nenl-1
      dic['type']='Triangular'
   elif typ == 'gau':
      fun=np.zeros(nenl*6-1)
      for i in range(0,nenl*6-1):
         fun[i]=np.exp(-(i-nenl*3+1)**2/(2*nenl**2))
      rit=3*nenl-1
      dic['type']='Gaussian'
   elif typ == 'rec':
      fun=np.ones(nenl)
      rit=(nenl-1)/2
      dic['type']='Rectangular'
   else:
      dic['type']='Personal'
      fun=typ

   fun=fun*nenl/np.sum(np.abs(fun))

   out=signal.convolve(base,fun)
   out=GD.gd(out,ini=ini-rit/enl,dx=sstep)
   fun=GD.gd(fun,ini=-rit)

   mu,sig=param_from_hist(out)
   dic['mu_hist']=mu
   dic['sig_hist']=sig

   if verb > 1:
      mu,sig=param_from_hist(fun)
      dic['mu_fun']=mu
      dic['sig_fun']=sig
      dic['tictoc']=BASIC.toc()
      GD.newfig()
      GD.plot_gd(out)

   BASIC.show_simp(dic)

   return out,dic,fun


def param_from_hist(hist):
# estimate population parameters from histogram
#  hist   gd containing the histogram
   x=GD.x_gd(hist)
   y=hist.y
   dx=hist.dx
   y=y/(sum(y)*dx)
   mu=np.sum(x*y)*dx
   sig=np.sqrt(np.sum(y*(x-mu)**2)*dx)

   return mu,sig


# Power Spectra ----------------------------

def gd_pows(ingd,npiece=1,res=1,shift=1,nobias=1,notrend=1,window=2,singleb=1,
sqr=0,center=1):
# Standard power spectrum estimation
# npiece    number of pieces (without shift)
# res       resoltion (minimal 1)
# shift     for interlaced pieces (1 no interlace, 0.5 one half shift)
# nobias    =1 subtract bias
# notrend   =1 subtract trend
# window    -0 no, 1 bartlett, 2 hann (default), 3 flat top cosine
# singleb   =1 single side band (for real data, normalized for single band)
# sqr       =1 square root of spectrum
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
   S=np.zeros(lS2)
   df=1/(lS*dx)

   npie=0

   for i in rang:
      i1=i
      i2=i1+n
      npie+=1
      # print(y.shape)
      # print(win.shape)
      yy=y[i1:i2]*win
      YY=np.zeros(lS)
      if isinstance(yy[1],complex):
         YY=YY+1j*YY
      YY[0:n]=yy
      YY=np.abs(fft(YY))**2
  #    print(np.shape(S))
  #    print(np.shape(YY))
      S+=YY[0:lS2]

   coef=dx/(npie*lS2)
   S=S*coef
   if sqr == 1:
      S=np.sqrt(S)
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

def spectrogram_gd(ingd,l,zenh=2,shif=0.5,win='tuckey'):
# Spectrogram
#
#  ingd    input gd or array
#  l       base length (pieces length)
#  zenh    zero padding enhancement fraction (>= 1)
#  shif    shift fraction (<= 1)
#  wind    window ('tuckey','hann','no')

   if isinstance(ingd,np.ndarray):
      dx=1
      N=len(ingd)
   else:
      dx=ingd.dx
      N=ingd.n
      ingd=ingd.y

   nfft=int(l*zenh)
   noverl=int(shif*l)
   ff,tt,spec=signal.spectrogram(ingd,fs=1/dx,nfft=nfft,nperseg=l,noverlap=noverl)
   ff=ff.squeeze()

   print('N,nfft,noverl,nfft,nperseg,noverl,lent,lenf',N,nfft,noverl,nfft,l,noverl,len(tt),len(ff))

   return spec,tt,ff



# Period analysis ------------------------

def gd_period(ingd,per,nbin=48,nharm=5,ph=0,preproc=1):
# period analysis
#  ingd     gd or array
#  per      period (numeric or physical (only for gds))
#           Physical period: 'day', 'sid', 'week'
#  nbin     number of bins in the period
#  nharm    number of interesting harmonics
#  ph       phase (in deg)
#  preproc  pre-processing 0 nothing, 1 abs, 2 square
#
# The gd sampling time must be in s, if physical periods are considered

   if isinstance(ingd,np.ndarray):
      ingd=GD.gd(ingd)
      cc='arr'
   else:
      cc='gd'

   y=ingd.y
   if preproc == 1:
      y=np.abs(y)
   if preproc == 2:
      y=np.abs(y)
      y=y*y
   x=GD.x_gd(ingd)
   N=ingd.n

   ini=0
   dx=360/nbin

   if isinstance(per,str):
      if cc == 'arr':
         per=float(input('No physical period for arrays, give me a number '))
      else:
         cont=ingd.cont
         try:
            t0=BASIC.val_from_key(cont,'t0')
            week,day,sid,locsid=ASTROTIME.mjd_phase(t0)
            if per == 'week':
               per=7*86400
               dx=7/nbin
               if ph == 0:
                  ph=week
            if per == 'day':
               per=86400
               dx=24/nbin
               if ph == 0:
                  ph=day
            if per == 'sid':
               per=0.9972695663290856*86400 # 86164.090530833/86400 as 2000/1/1
               dx=24/nbin
               if ph == 0:
                  ph=sid
            if per == 'locsid':
               per=0.9972695663290856*86400 # 86164.090530833/86400 as 2000/1/1
               dx=24/nbin
               if ph == 0:
                  ph=locsid
         except:
            print('*** This is not a dated gd: no t0 defined')
            t0=0
            if per == 'week':
               per=7*86400
               dx=7/nbin
            if per == 'day':
               per=86400
            if per == 'sid':
               per=0.9972695663290856*86400 # 86164.090530833/86400 as 2000/1/1
               dx=24/nbin
            if per == 'locsid':
               per=0.9972695663290856*86400 # 86164.090530833/86400 as 2000/1/1
               dx=24/nbin
         
   ii=np.int16(((x/per)%1)*nbin)
   period=np.zeros(nbin)
   iperiod=np.zeros(nbin)

   for i in range(N):
      kk=ii[i]
      period[kk]+=y[i]
      iperiod[kk]+=np.abs(np.sign(y[i]))

   period=period/(iperiod+0.1)
   meanp=np.mean(period)

   harm=0
   perclean=0
   f=np.fft.fft(period-meanp)
   harm=f[0:nharm+1]
   f[nharm+1:nbin-nharm]=0
   perclean=np.fft.ifft(f)+meanp
   win=iperiod
   period=GD.gd(period,dx=dx)
   perclean=GD.gd(perclean,dx=dx)
   win=GD.gd(win,dx=dx)

   return period,meanp,harm,perclean,win
   


def gd_tperiod(ingd,per,nbin=(20,48),nharm=5,ph=0,preproc=1):
# period analysis with time variation
#  ingd     gd or array
#  per      period (numeric or physical (only for gds))
#           Physical period: 'day', 'sid', 'week'
#  nbin     number of bins in the oservation time [0] and in the period [1]
#  nharm    number of interesting harmonics
#  ph       phase (in deg)
#  preproc  pre-processing 0 nothing, 1 abs, 2 square
#
# The gd sampling time must be in s, if physical periods are considered

   tnbin=nbin[0]
   pnbin=nbin[1]
   if isinstance(ingd,np.ndarray):
      ingd=GD.gd(ingd)
      cc='arr'
   else:
      cc='gd'

   y=ingd.y
   if preproc == 1:
      y=np.abs(y)
   if preproc == 2:
      y=np.abs(y)
      y=y*y
   x=GD.x_gd(ingd)
   N=ingd.n

   ini=0
   dx=360/pnbin

   if isinstance(per,str):
      if cc == 'arr':
         per=float(input('No physical period for arrays, give me a number '))
      else:
         cont=ingd.cont
         try:
            t0=BASIC.val_from_key(cont,'t0')
            week,day,sid,locsid=ASTROTIME.mjd_phase(t0)
            if per == 'week':
               per=7*86400
               dx=7/pnbin
               if ph == 0:
                  ph=week
            if per == 'day':
               per=86400
               dx=24/pnbin
               if ph == 0:
                  ph=day
            if per == 'sid':
               per=0.9972695663290856*86400 # 86164.090530833/86400 as 2000/1/1
               dx=24/pnbin
               if ph == 0:
                  ph=sid
            if per == 'locsid':
               per=0.9972695663290856*86400 # 86164.090530833/86400 as 2000/1/1
               dx=24/pnbin
               if ph == 0:
                  ph=locsid
         except:
            print('*** This is not a dated gd: no t0 defined')
            t0=0
            if per == 'week':
               per=7*86400
               dx=7/pnbin
            if per == 'day':
               per=86400
            if per == 'sid':
               per=0.9972695663290856*86400 # 86164.090530833/86400 as 2000/1/1
               dx=24/pnbin
            if per == 'locsid':
               per=0.9972695663290856*86400 # 86164.090530833/86400 as 2000/1/1
               dx=24/pnbin
         
   ii=np.int16(((x/per)%1)*pnbin)
   tperiod=np.zeros([tnbin,pnbin])
   iperiod=np.zeros([tnbin,pnbin])

   for i in range(N):
      jj=int(tnbin*i/N)
      kk=ii[i]
      tperiod[jj,kk]+=y[i]
      iperiod[jj,kk]+=np.abs(np.sign(y[i]))

   tperiod=tperiod/(iperiod+0.1)
   meanp=np.mean(tperiod)

   harm=0
   tperclean=0
   # f=np.fft.fft2(tperiod-meanp)
   # harm=f[0:nharm+1]
   # f[nharm+1:nbin-nharm]=0
   # tperclean=np.fft.ifft(f)+meanp
   twin=iperiod
   tperiod=GD.gd(tperiod,dx=dx)
   # tperclean=GD.gd(tperclean,dx=dx)
   twin=GD.gd(twin,dx=dx)

   return tperiod,meanp,harm,tperclean,twin
   

def gd_worm():
   pass


# From Scipy -------------------------

dist_continu = [d for d in dir(stats) if
                isinstance(getattr(stats, d), stats.rv_continuous)]
dist_discrete = [d for d in dir(stats) if
                 isinstance(getattr(stats, d), stats.rv_discrete)]

# The main public methods for continuous RVs are:
#
# rvs: Random Variates
# pdf: Probability Density Function
# cdf: Cumulative Distribution Function
# sf: Survival Function (1-CDF)
# ppf: Percent Point Function (Inverse of CDF)
# isf: Inverse Survival Function (Inverse of SF)
# stats: Return mean, variance, (Fisher’s) skew, or (Fisher’s) kurtosis
# moment: non-central moments of the distribution
#
# https://docs.scipy.org/doc/scipy/tutorial/stats/discrete.html
# https://docs.scipy.org/doc/scipy/tutorial/stats/continuous.html

def rand_data(N,typ,par1=0,par2=1):
# Random numbers
#  N     how many
#  typ   type
#  pars  parameters
#
#    typ      pars
#   norm     mu,sigma

   if typ == 'norm': # par1 = mean, par2 = std
      dat=stats.norm.rvs(size=N)*par2+par1
   if typ == 'unif': # par1 = mean, par2 = length
      dat=(stats.uniform.rvs(size=N)+par1-0.5)*par2
   if typ == 'chi2': # par1 = dof
      df=par1
      if df == 0:
         df=4
      dat=stats.chi2.rvs(df,size=N)
   if typ == 'binom': # par1 = n, par2= p
      n=par1
      p=par2
      if n == 0:
         n=20
         p=0.5
      dat=stats.binom.rvs(n,p,size=N)

   return dat
