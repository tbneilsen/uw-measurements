#%%
from TimeGate_UnderwaterTank import uwsoundspeed
from readLogFile import readLogFile

filename = 'ID001_001log.txt'
location = '/home/byu.local/sh747/underwater/uw-measurements-tank/2021/2021-11-12B/2021-11-12_scan3'+'/'

_,_,tempWater,_,_,_,_,_,depth,_,_,_,_,_,_ = readLogFile(filename,location)

D=depth/2
S=0.03
model='Garrett'
c=uwsoundspeed(D,tempWater,S,model)

print(c)
# %%
