'''
Simple intervals operations
'''
import SERV,GD

A=SERV.intervals(20,ini=[3,12],fin=[7,16])
SERV.show_interv(A)

B=SERV.intervals(20,ini=[0,11],fin=[6,16])
SERV.show_interv_2(B)

C=SERV.intervals(20,ini=[3,12],fin=[7,20])
SERV.show_interv_2(C)

AA=SERV.invert_interv(A)
SERV.show_interv_2(AA)
BB=SERV.invert_interv(B)
SERV.show_interv_2(BB)
CC=SERV.invert_interv(C)
SERV.show_interv_2(CC)

D=SERV.interv_and(A,B,C)
SERV.show_interv_2(D)
E=SERV.interv_or(A,B,C)
SERV.show_interv_2(E)

cover=SERV.interv_coverage(A,B,C)
GD.newfig()
GD.plot_gd(cover)

I2D=SERV.collect_interv(A,B,C)
SERV.check_interv(I2D)
SERV.show_interv(I2D)

mask=A=SERV.intervals(20,ini=[3,12],fin=[7,16])
SERV.show_interv(A)

B=SERV.intervals(20,ini=[0,11],fin=[6,16])
SERV.show_interv_2(B)

C=SERV.intervals(20,ini=[3,12],fin=[7,20])
SERV.show_interv_2(C)

AA=SERV.invert_interv(A)
SERV.show_interv_2(AA)
BB=SERV.invert_interv(B)
SERV.show_interv_2(BB)
CC=SERV.invert_interv(C)
SERV.show_interv_2(CC)

D=SERV.interv_and(A,B,C)
SERV.show_interv_2(D)
E=SERV.interv_or(A,B,C)
SERV.show_interv_2(E)

cover=SERV.interv_coverage(A,B,C)
GD.newfig()
GD.plot_gd(cover)

I2D=SERV.collect_interv(A,B,C)
SERV.check_interv(I2D)
SERV.show_interv(I2D)

mask=SERV.interv2mask(I2D)
