#!/usr/bin/python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import random
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.client import FDSNNoDataException
from obspy import UTCDateTime
from obspy import read
from sys import exit, argv
import tqdm
from glob import glob
import os, datetime

try:
    NET=argv[0]
    STA=argv[1]
    LOC=argv[2]
    CHAN=argv[3]
except:
    NET,STA,LOC,CHAN="FR", "CURIE", "00", "HHZ"
    NET,STA,LOC,CHAN="AM", "R8F35", "00", "SHZ" # BG
    NET,STA,LOC,CHAN="AM", "RFD43", "00", "EHZ" # RB
    NET,STA,LOC,CHAN="AM", "RBFD5", "00", "SHZ" # NB
    NET,STA,LOC,CHAN="AM", "R9F1B", "00", "SHZ"

nslc = "{}.{}.{}.{}".format(NET, STA, LOC, CHAN)
# make sure that wildcard characters are not in nslc
nslc = nslc.replace("*", "").replace("?", "")

# ---------------------
headers = ['rms', 'Datetime']
dtypes = {'rms': 'float', 'Datetime': 'str'}
parse_dates = ['Datetime']

#fn_pattern="./csv/*_{}.csv".format(nslc)
fn_pattern="./csv/*_%s_%s_%s_%s.csv" % (NET,STA,LOC,CHAN)
for csvfile in glob(fn_pattern):
    print ("csvfile=",csvfile)
    if os.path.isfile(csvfile):
        try:
            df=pd.read_csv(csvfile, sep=',', header=0, names=headers, dtype=dtypes, parse_dates=parse_dates)
            Q=np.quantile(df.rms, 0.90)
            df.rms=np.clip(df.rms, 0, Q)
        except:
            raise()
    try:
        dfALL=pd.concat([dfALL, df], ignore_index=True)
    except:
        dfALL=df

dfALL['Datetime']=pd.to_datetime(dfALL['Datetime'], unit='s')
dfALL=dfALL.sort_values(by=['Datetime'])
fig=plt.figure()
ax=plt.axes()
##ax.plot(stD[0].times('matplotlib'), stD[0].data, 'k-', label="Trace", linewidth=0.2) 
ax.plot(dfALL.Datetime, dfALL.rms, 'r-', label="%s"%nslc, linewidth=0.5) 
fig.autofmt_xdate()
plt.legend(loc=0)
plt.grid()
plt.show()

