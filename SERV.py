    # Copyright (C) 2023  Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
        Module SERV

Service functions

Sections:

> Service computational routines    -> dummy_comp
> Intervals                         -> dummy_inter
> no-data management                -> dummy_nodata

'''

def sections():
    sec=[
    ]
    return sec

import numpy as np
import GD,GD2,BASIC,STAT,SIGNAL
import copy,time

# Service computational routines ---------------------------------

def dummy_comp():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

def rota(inp,n):
# circular shift of an array or a list
# n is a positive or negative integer

    isarr=0
    if isinstance(inp,np.ndarray):
        isarr=1
        inp=list(inp)        

    inp = inp[n:] + inp[:n]

    if isarr:
        inp=np.array(inp)

    return inp  


def thresh(inp,t1,t2):
#  threshold function (0 if out, 1 if in)
#     a=thresh(b,t1,t2)
#  if t1 < t2  is in the interval, else it is out
#  (used, e.g., in atan3)
    if t1 < t2:
        a1=np.floor((np.sign(inp-t1)+2)/2)
        a=a1*np.floor((np.sign(t2-inp)+2)/2)
    else:
        a1=np.floor((np.sign(t2-inp)+1)/2)
        a=a1+np.floor((np.sign(inp-t1)+1)/2)

    return a


def atan3(inp):
# operates on complex arrays: computes the phases in a domain larger than 360 deg
# computes atan3 (in number of turns; 1 turn = 360 degrees)
    pi2=np.pi*2
    n=len(inp)
    b=np.angle(inp)/pi2
    dif=np.zeros(n)
    dif[1:n]=np.diff(b)
    s=rota(b,-1)
    s[0]=0
    s=np.sign(s)

    cdif=thresh(np.abs(dif),0.5,2)
    cdif=-cdif*s+dif
    
    cdif1=rota(cdif,1)
    cdif2=rota(cdif,-1)
    i= cdif1*cdif2 > 0 & cdif1*cdif < 0
    cdif[i]=(cdif1(i)+cdif2(i))/2

    a=np.cumsum(cdif)

    return a


def mask(n,n1,**pars): 
# Mask function (e.g. for frequency filters)
#  n=length, n1 cut or list cuts, 
#    ** lenw (def 3), mode 'symfr'
    mas=np.ones(n)
    if isinstance(n1,int):
        nw=1
        n1=[n1]
    nw=len(n1)
    if 'lenw' in pars:
        lenw=pars('lenw')
    else:
        lenw=3
    ii=np.arange(1,lenw+1)   
    w=(np.cos(ii*np.pi/(lenw+1))+1)/2
    type(w)
    type(mas)
    w1=np.flip(w)
    fin=0

    for i in range(0,nw,2):
        len1=n1[i]-fin
        if i == nw-1:
            fin=n
        else:
            fin=n1[i+1]
            if n-fin > lenw:
                mas[n1[i+1]:n1[i+1]+lenw]=w1

        if len1 > lenw:
            mas[n1[i]-lenw:n1[i]]=w[0:lenw]

        mas[n1[i]+1:fin]=0  
        if 'mode' in pars:
            mode=pars['mode']
            if mode == 'symfr':
                nn=int(np.floor(n/2)-1)
                mas[n:n-nn-1:-1]=mas[1:nn+1]

    return mas


def ifr(ingd,fr):
# Frequency sample index for a given frequency
#  ingd   the gd to be transformed
#  fr     the requested frequency
    n=ingd.n
    dx=ingd.dx
    frmax=1/dx
    i=np.round(fr*n/frmax)
    if i >= n:
        print(' *** OUT OF THE BAND: i = ',i,' max ',n-1)

    return i


def range_shift(ini,fin,shift,lseg):
#  range with shift and lseg
# 
    n=0
    inis=[]
    inis.append(ini)
    fins=[]
    fins.append(inis[n]+lseg)
    while fins[n] <= fin:
        n+=1
        inis.append(inis[n-1]+shift)
        fins.append(inis[n]+lseg)

    inis=inis[0:n]
    fins=fins[0:n]
    return n,inis,fins


def vec_ccdot(in1,in2):
# vectorial complex conjugate dot product
#  in1,in2   equal shape arrays (nxm)
# output n array

    aa=in1.shape
    n=aa[0]

    out=np.zeros(n)

    for i in range(n):
        iin1=in1[i]
        iin2=in2[i]
        out[i]=np.vdot(iin1,iin2)

    return out


# Intervals --------------------------------

def dummy_inter():
    '''
    Time (or other) interval management
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun

class intervals():
    """
    Intervals class -------------------------

            intervals creation

    To define intervals in a 1-D array

    The attributes are:

    > lar    array length
    > nar    second dimension of the array (def=1)
    > cover  total coverage (number of 1s, per row)
    > typ    1 or 2 as GDs
    > x0     initial abscissa value (def=0)
    > dx     abscissa step (def=1)
    > x      abscissa (in the case of typ=2)
    > ini    interval init index
    > fin    interval end index (excluded value) 
    > xini    interval abscissa init
    > xfin    interval abscissa end (excluded value) 
    > label  e.g. 'hole', 'data', 'good',...; def ''
    """

    def __init__(self,lar,**gdpar): # y ordinate or n
        self.lar=lar
        self.typ=1
        if 'nar' in gdpar:
            self.nar=gdpar['nar']
        else:
            self.nar=1
        if 'x0' in gdpar:
            self.x0=gdpar['x0']
        else:
            self.x0=0
        if 'dx' in gdpar:
            self.dx=gdpar['dx']
        else:
            self.dx=1
        if 'x' in gdpar:
            self.x=np.array(gdpar['x'])
            if len(self.x) > 0:
                self.typ=2
        else:
            self.x=[]
            self.typ=1
        if 'ini' in gdpar:
            ini=gdpar['ini']
            self.ini=ini
            # if ini.ndim == 2:
            #     aa=ini.shape
            #     self.nar=aa[1]           
        else:
            self.ini=[]
        if 'fin' in gdpar:
            self.fin=gdpar['fin']
        else:
            self.fin=[]
        if 'xini' in gdpar:
            xini=gdpar['xini']
            self.xini=xini           
        else:
            self.xini=[]
        if 'xfin' in gdpar:
            self.xfin=gdpar['xfin']
        else:
            self.xfin=[]
        if 'label' in gdpar:
            self.label=gdpar['label']
        else:
            self.label=[]

        if self.nar == 1:
            self.ini=np.array(self.ini)
            self.fin=np.array(self.fin)

        self.cover=cover_calc(self)
  
    def set_data(self):
        self.label='data'

    def set_hole(self):
        self.label='hole'

    def invert(self):
        inver=invert_interv(self)
        return inver
    
    def check(self):
        ch=check_interv(self)
        print(ch)

    def x_(self):
        x=x_interv(self)
        return x
    
    def show(self):
        show_interv(self)

    def mask(self):
        mask=interv2mask(self)
        return mask

  
def x2ind(interv,x):
    if isinstance(ind,np.ndarray):
        ind=BASIC.ind_from_inidx(x, interv.x0, interv.dx)
    elif isinstance(x,list):
        nn=len(x)
        ind=[]
        for i in range(nn):
            ind0=BASIC.ind_from_inidx(x[i], interv.x0, interv.dx)
            ind.append(ind0)
    else:
        print('*** Error: x is ',type(x))
        return False
    
    return ind


def ind2x(interv,ind):
    if isinstance(ind,np.ndarray):
        if interv.typ == 1:
            x=interv.x0+ind*interv.dx
        else:
            print('not yet implemented')
            return False
    elif isinstance(ind,list):
        nn=len(ind)
        x=[[]]*nn
        if interv.typ == 1:
            for i in range(nn):
                xx=interv.x0+ind[i]*interv.dx
                x[i]=xx
        else:
            print('not yet implemented')
            return False
    else:
        print('*** Error: ind is ',type(ind))
        return False

    return x


def full_ind2x(interv):
    xini=ind2x(interv,interv.ini)
    xfin=ind2x(interv,interv.fin)
    interv.xini=xini
    interv.xfin=xfin

    return interv


def cover_calc(interv):
# computes cover
    nar=interv.nar
    ini=interv.ini
    fin=interv.fin
    
    if nar == 1:
        cover=sum(fin-ini)
    else:
        cover=np.zeros(nar,dtype=int)
        for i in range(nar):
            cover[i]=sum(fin[i]-ini[i])
        
    return cover


def mask2interv(mask, cc=1):
    '''
    mask (1-0 array) intervals
    m      mask
    cc     =0 0 intervals
            =1 1 intervals
    '''

    nar,lar=BASIC.array_rowcol(mask)
    n = lar
    interv=[]

    if cc == 0:
        mask = (1-mask)

    dm = np.diff(mask)
    if nar == 1:
        inz = np.nonzero(np.array(dm))[0][0]
        minnz = dm[inz]
        if minnz == 1:
            dm = np.insert(dm, 0, 0)
        else:
            dm = np.insert(dm, 0, 1)

        ini = np.argwhere(dm == 1)
        fin = np.argwhere(dm == -1)
        if len(ini) > len(fin):
            fin = np.append(fin, n)

        ini=ini.squeeze()
        fin=fin.squeeze()
    else:
        ini=[]
        fin=[]
        for i in range(nar):
            dm0=dm[i]
            inz = np.nonzero(dm0)[0][0]
            minnz = dm0[inz]
            if minnz == 1:
                dm0 = np.insert(dm0, 0, 0)
            else:
                dm0 = np.insert(dm0, 0, 1)

            ini0 = np.argwhere(dm0 == 1)
            fin0 = np.argwhere(dm0 == -1)
            if len(ini0) > len(fin0):
                fin0 = np.append(fin0, n)

            ini0=ini0.squeeze()
            fin0=fin0.squeeze()

            ini.append(ini0)
            fin.append(fin0)
    
    interv=intervals(n,ini=ini,fin=fin)
    interv.nar=nar

    return interv


def interv2mask(interv):
    lar=interv.lar
    nar=interv.nar
    if nar == 1:
        mask=np.zeros(lar,dtype='int8')
    else:
        mask=np.zeros([nar,lar],dtype='int8')
    if nar == 1:
        nin=len(interv.ini)
        for i in range(nin):
            mask[interv.ini[i]:interv.fin[i]]=1
    else:
        for j in range(nar):
            nin=len(interv.ini[j])
            for i in range(nin):
                mask[j,interv.ini[j][i]:interv.fin[j][i]]=1
    
    return mask


def interv2dict(interv):
    intdic={'lar':interv.lar,'nar':interv.nar,'ini':interv.ini,'fin':interv.fin,
           'x0':interv.x0,'dx':interv.dx,'x':interv.x,'typ':interv.typ}
    return intdic


def dict2interv(intdic):
    interv=intervals(intdic['lar'],nar=intdic['nar'],x0=intdic['x0'],dx=intdic['dx'],x=intdic['x'],
                     ini=intdic['ini'],fin=intdic['fin'])
    return interv


def x_interv(interv):
# intervals abscissa (real or realized)

    if isinstance(interv,dict):
        interv=dict2interv(interv)

    if interv.typ == 1 :
        x=interv.x0+np.arange(interv.lar)*interv.dx
    else:
        x=interv.x

    x=BASIC.deshape(x)

    return x 


def show_interv(interv,verb=1):
    '''
    show intervals attributes
    '''
    out=check_interv(interv)
    nar=interv.nar
    cover=interv.cover
    if nar > 1:
        verb=0
    if not out:
        print('*** Wrong interval')
        return
    if nar == 1:
        nint=len(interv.ini)
    else:
        nint=np.zeros(interv.nar,dtype=int)
        for i in range(interv.nar):
            nint[i]=len(interv.ini[i])

    print('\n')
    print('lar    =',interv.lar)
    print('nar    =',nar)
    print('cover  =',cover)
    print('x0     =',interv.x0)
    print('dx     =',interv.dx)
    print('x      =',interv.x)
    print('typ    =',interv.typ)
    print('N.int  =',nint)
    time.sleep(5)
    print('ini    =',interv.ini)
    print('fin    =',interv.fin)

    if verb == 1:
        ini=interv.ini
        fin=interv.fin
        print('   Intervals')
        for i in range(nint):
            print(i,' - ',ini[i],'<->',fin[i],' -->',fin[i]-ini[i])


def show_interv_2(interv):
    '''
    show intervals - short version
    '''
    ini=interv.ini
    fin=interv.fin
    nar=interv.nar

    print('   Intervals')
    for i in range(nar):
        print(i,' - ',ini[i],'<->',fin[i],' -->',fin[i]-ini[i])


def check_interv(inter):
# checks intervals object
    out=True
    try:
        if isinstance(inter,dict):
            inter=dict2interv(inter)
        lar=inter.lar
        nar=inter.nar
        ini=inter.ini
        fin=inter.fin
        typ=inter.typ
        print('lar =',lar)
        print('nar =',nar)
        if nar == 1:
            if len(ini) != len(fin):
                out=False
                print('ini and fin different length',len(ini),len(fin))
            kel=0
            for el in ini:
                if el > lar:
                    out=False
                    print('ini[',kel,'] = ',el,'> ',lar)
                kel+=1
            kel=0
            for el in fin:
                if el > lar:
                    out=False
                    print('fin[',kel,'] = ',el,'> ',lar)
                kel+=1
        else:
            for i in range(nar):
                if len(ini[i]) != len(fin[i]):
                    out=False
                    print('row',i,'-> ini and fin different length',len(ini),len(fin))
                kel=0
                for el in ini[i]:
                    if el > lar:
                        out=False
                        print('row',i,'-> ini[',kel,'] = ',el,'> ',lar)
                    kel+=1
                kel=0
                for el in fin[i]:
                    if el > lar:
                        out=False
                        print('row',i,'-> fin[',kel,'] = ',el,'> ',lar)
                    kel+=1

        if typ == 2:
            if len(inter.x) != lar:
                out=False
                print('length(x) =',len(inter.x))
    except:
        out=False
        print('intervals generic error')

    return out


def invert_interv(interv):
    out=check_interv(interv)
    if not out:
        print('*** Wrong input interval')
        return out
    if isinstance(interv,dict):
        inter=dict2interv(interv)
        conv=1
    else:
        conv=0
    outint=copy.deepcopy(interv)
    lar=interv.lar
    nar=interv.nar
    ini=interv.ini
    fin=interv.fin

    if nar == 1:
        nin=len(ini)

        inio=copy.copy(fin)
        fino=copy.copy(ini)
        for i in range(nin):
            if i < nin-1:
                fino[i]=ini[i+1]
            else:
                fino[i]=lar

        if ini[0] > 0:
            inio=np.insert(inio,0,0)
            fino=np.insert(fino,0,ini[0])
        non=len(inio)
        if fin[nin-1] == lar:
            inio=np.delete(inio,non-1)
            fino=np.delete(fino,non-1)

        outint.ini=inio
        outint.fin=fino
    else:
        for i in range(nar):
            nin=len(ini[i])

            inio=copy.copy(fin[i])
            fino=copy.copy(ini[i])
            for j in range(nin):
                if j < nin-1:
                    fino[j]=ini[i][j+1]
                else:
                    fino[j]=lar

            if ini[i][0] > 0:
                inio=np.insert(inio,0,0)
                fino=np.insert(fino,0,ini[i][0])
            non=len(inio)
            if fin[i][nin-1] == lar:
                inio=np.delete(inio,non-1)
                fino=np.delete(fino,non-1)

            outint.ini[i]=inio
            outint.fin[i]=fino

    if interv.label == 'data':
        interv.label='hole'
    if interv.label == 'hole':
        interv.label='data'

    return outint


def interv_and(*inter):
    '''
    AND of undefined number of arguments
    --- only for 1-D intervals ---
    '''

    Nint=len(inter)
    print(Nint,'arguments')

    for i in range(Nint):
        if inter[i].nar > 1:
            print('*** Error: multidimensional array')
            return False
        if inter[i].lar != inter[0].lar:
            print('*** Error: different lar: not compatible')
            return False

    if Nint == 1:
        return inter[0]
    else:
        mas=inter[0].mask()
        for i in range(1,Nint):
            ma=inter[i].mask()
            mas=mas*ma
    
    outint=mask2interv(mas)

    return outint


def interv_or(*inter):
    '''
    OR of undefined number of arguments
        --- only for 1-D intervals ---
    '''

    Nint=len(inter)
    print(Nint,'arguments')

    for i in range(Nint):
        if inter[i].nar > 1:
            print('*** Error: multidimensional array')
            return False
        if inter[i].lar != inter[0].lar:
            print('*** Error: different lar: not compatible')
            return False

    if Nint == 1:
        return inter[0]
    else:
        mas=inter[0].mask()
        for i in range(1,Nint):
            ma=inter[i].mask()
            mas=mas+ma
    
        mas=np.sign(mas)

    outint=mask2interv(mas)

    return outint


def interv_coverage(*inter):
    '''
    Coverage function of undefined number of arguments
    --- only for 1-D intervals ---
    '''

    Nint=len(inter)
    print(Nint,'arguments')

    for i in range(Nint):
        if inter[i].nar > 1:
            print('*** Error: multidimensional array')
            return False
        if inter[i].lar != inter[0].lar:
            print('*** Error: different lar: not compatible')
            return False

    if Nint == 1:
        return inter[0]
    else:
        cover=inter[0].mask()
        for i in range(1,Nint):
            ma=inter[i].mask()
            cover=cover+ma

    return cover


def collect_interv(*inter):
    '''
    Collects a number of 1-D intervals to create a 2-D intervals
    '''
    nar=len(inter)
    lar=inter[0].lar
    ini=[]
    fin=[]
    for i in range(nar):
        ini.append(inter[i].ini)
        fin.append(inter[i].fin)

    outinter=intervals(lar,nar=nar,ini=ini,fin=fin)

    return outinter


def sel_interv(ingd,interv):
    '''
    select data in the intervals
    ingd    gd, gd2 or np array
    '''
    lar=interv.lar
    nar=interv.nar
    if isinstance(ingd,GD.gd):
        ingd=ingd.y
    if isinstance(ingd,GD2.gd2):
        ingd=ingd.y

    y=ingd

    if nar == 1:
        nint=len(interv.ini)
        yy=np.array([])
        for i in range(nint):
            yy.append(y[interv.ini[i]:interv.fin[i]])
    else:
        yy=[]
        for k in range(nar):
            nint=len(interv.ini[k])
            yy0=np.array([])
            for i in range(nint):
                yy0.append(y[k][interv.ini[k][i]:interv.fin[k][i]])
                yy.append(yy0)

    return yy


def interv_reduc(inter,mi,ma=0):
    '''
    Interval reduction, based on the length

    mi   minimal interval length
    ma   maximal interval length; ma=0 -> no limit
    '''

    lar=inter.lar
    if ma == 0:
        ma=lar+1

    nin=len(inter.ini)
    ini=np.zeros(nin)
    fin=np.zeros(nin)
    redinter=copy.deepcopy(inter)
    ii=0

    for i in range(nin):
        w=inter.fin[i]-inter.ini[i]
        if w >= mi and w <= ma:
            ini[ii]=inter.ini[i]
            fin[ii]=inter.fin[i]
            ii+=1

    redinter.ini=ini[0:ii]
    redinter.fin=fin[0:ii]

    return redinter

    
def win_interv(inter,lwin,verb=1):
    '''
    Spectral window creation for intervalled data

    inter   intervals object
    lwin    window length parameter
    '''

    mask=interv2mask(inter)
    
    w=np.exp(-1/lwin)
    a=[1,-w]
    b=1
    win=SIGNAL.FiltFilt(mask,a,b)
    win=win*mask
    smask=np.sum(mask)
    swin=np.sum(win)
    print('smask/swin =',smask/swin)
    win=win*np.sqrt(smask/swin)

    if verb == 1:
        GD.newfig()
        GD.plot_gd(mask)
        GD.plot_gd(win)
        GD.post_plot('Windows','','')

        sp=STAT.gd_pows(mask)
        fsp=STAT.gd_pows(win)
        
        GD.newfig()
        GD.plot_gd(sp)
        GD.plot_gd(fsp)
        GD.post_plot('Windows spectra','','')

    return win



# no-data management -------------------------

def dummy_nodata():
    '''
    '''
    clas=[

    ]
    fun=[
        
    ]
    return clas,fun


def data_interv(dat,eps=1.e-5,lenz=2,dataon=1):
    '''
    finds data intervals

        dat     gd, gd2 or np array
        eps     relative level of no data
        lenz    min length of zero for no data
        dataon  = 0 -> holes, = 1 -> data
    '''

    if isinstance(dat, np.ndarray) == False:
        dat = dat.y
    dat=abs(dat)
    dsh = dat.shape
    if len(dsh) == 1:
        nr = 1
        nc = len(dat)
    else:
        nr = dsh[0]
        nc = dsh[1]

    lar=nc
    nar=nr

    dim = (nr, nc)
    yy = dat.flatten()
    yy = yy[~np.isnan(yy)]
    yy = np.abs(yy)
    mma = np.max(yy)
    # eps1a=eps*mma/np.sqrt(nr)
    eps1a = eps*mma
    if mma == 0:
        return print('*** All null values')

    Inan = [None]*nr
    Ini = [None]*nr
    Fin = [None]*nr

    print('nr=', nr)

    for i in range(nr):
        if nr == 1:
            y = dat
        else:
            y = dat[i]

        zer = np.zeros(nc)
        iii = np.argwhere(~np.isnan(y))
        yy = abs(y[iii])
        ma = max(yy)
        if ma < mma:
            ma = mma
        eps1 = eps*ma
        inan = np.argwhere(np.isnan(y))
        zer[inan] = 1
        Inan[i] = inan

        y[inan] = 0

        for k in range(nc-lenz+1):
            hol=1
            if y[k] <= eps1:
                for iz in range(lenz-1):
                    if y[k+1+iz] > eps1:
                        hol=hol*0
                if hol == 1:
                    for iz in range(lenz):
                        zer[k+iz] = 1

        print(len(zer),sum(zer))
        interv = mask2interv(zer)
        Ini[i] = interv.ini
        Fin[i] = interv.fin

    if nr == 1:
        Ini=Ini[0]
        Fin=Fin[0]
        Inan=Inan[0]

    interv=intervals(nc,nar=nr,ini=Ini,fin=Fin)
    interv.label='hole'

    if dataon == 1:
        interv=invert_interv(interv)
        interv.label='data'

    interv.Inan=Inan

    return interv
    

def show_nodata(Ini, Fin, dim, dat=0):
    '''
    Shows (and sets) nodata
     Ini,Fin,dim   output of findnodata
     dat           array or gd or gd2 to set
    '''
    nr = dim[0]
    nc = dim[1]
    Nint = 0
    Noel = 0
    Nint = np.zeros(nr)
    Noel = np.zeros(nr)

    for i in range(nr):
        Nint[i] = len(Ini[i])
        Noel[i] = np.sum(Fin[i]-Ini[i])
        print('row ', i, '>', Nint[i], ' intervals ',
              Noel[i], ' absent elements')

    nint = np.sum(Nint)
    noel = np.sum(Noel)

    print('Total number of intervals and absent elements', nint, noel)

    return Nint, Noel


def stat_interv(inter,low=150):
    '''
    statistical analysis of an intervals object
    '''
    lar=inter.lar
    nar=inter.nar
    ini=inter.ini
    fin=inter.fin
    w=fin-ini
    hist=STAT.loghist(w,nbins=100)

    inter1=invert_interv(inter)
    ini1=inter1.ini
    fin1=inter1.fin
    w1=fin1-ini1
    hist1=STAT.loghist(w1,nbins=100)
    
    GD.newfig()
    P_H=GD.plot_helper(scale='lolo',mode='step')
    GD.plot_gd(hist,P_H=P_H)
    GD.plot_gd(hist1,P_H=P_H)
    GD.post_plot('Data and Hole interval length','len','')
    
    rich=copy.deepcopy(hist)
    rich.y=rich.y*rich.x
    rich1=copy.deepcopy(hist1)
    rich1.y=rich1.y*rich1.x
    GD.newfig()
    P_H=GD.plot_helper(scale='lolo',mode='step')
    GD.plot_gd(rich,P_H=P_H)
    GD.plot_gd(rich1,P_H=P_H)
    GD.post_plot('Data and Hole interval richness','len','')

    stat,Hist=GD.stat_gd(w,100)
    BASIC.show_simp(stat)
    
    stat,Hist=GD.stat_gd(w1,100)
    BASIC.show_simp(stat)
