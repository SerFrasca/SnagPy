'''
More complex interval operations
'''
import ML_PY,STAT,SERV,GD,GD2

# charge gd from a mat file v7

g=ML_PY.gd_lm7('L_C02_20170104_0060_0070_O2_tfstr.mat')
sp,tt,ff=STAT.gd_spectrogram(g,1000)

# intervals

intervg=SERV.data_interv(g)

maskg=SERV.interv2mask(intervg)
GD.newfig() 
GD.plot_gd(maskg) 

intervsp=SERV.data_interv(sp)

masksp=SERV.interv2mask(intervsp)
GD2.grey_map(masksp)

