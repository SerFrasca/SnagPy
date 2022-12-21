import GD,STAT

A=GD.gd(100000,dx=0.001)
A1=GD.set_gd(A,'sin',par2=100)
id(A)
id(A1)
B=GD.gd(100000,dx=0.001)
B1=GD.rand_gd(B,'norm')
C=A+B*0.2
#GD.plot_gd(C)
s=STAT.gd_welch(C,0.1,shift=0.5)
GD.semilogy_gd(s)