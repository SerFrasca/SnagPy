# examp_interv_3.py
# script for win_interp

from scipy import stats
import numpy as np
import SERV,GD
import STAT,SIGNAL

q1 = 0.0001
q0 = q1*2

p1 = 1-q1
p0 = 1-q0

N = 1000000

mask = np.zeros(N)
p = stats.uniform.rvs(size=N)

stat = 1

for i in range(N):
    if stat == 1:
        if p[i] > p1:
            stat = 0
    else:
        if p[i] > p0:
            stat = 1
    mask[i] = stat

inter = SERV.mask2interv(mask)
mask1=1-mask

sp = STAT.gd_pows(mask)

GD.newfig()
GD.plot_gd(sp)

w=np.exp(-1/10)
a=[1,-w]
b=1
fask=SIGNAL.FiltFilt(mask,a,b)
fsp = STAT.gd_pows(fask)

