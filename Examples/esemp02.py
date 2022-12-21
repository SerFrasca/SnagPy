import ML_PY,STAT,GD

g7=ML_PY.gdloadmat7('L_C02_20170104_0060_0070_O2_tfstr.mat')

s7=STAT.gd_pows(g7,res=4,npiece=4)

GD.semilogy_gd(s7)
