import numpy as np

'''
Service computational routines
'''

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
