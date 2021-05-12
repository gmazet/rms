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
import os, datetime

try:
    NET=argv[0]
    STA=argv[1]
    LOC=argv[2]
    CHAN=argv[3]
    PROVIDER =argv[4]
except:
    NET,STA,LOC,CHAN="FR", "CURIE", "00", "HHZ"
    NET,STA,LOC,CHAN="AM", "R8F35", "00", "SHZ" # BG
    NET,STA,LOC,CHAN="AM", "RFD43", "00", "EHZ" # RB
    NET,STA,LOC,CHAN="AM", "RBFD5", "00", "SHZ" # NB
    NET,STA,LOC,CHAN="AM", "R9F1B", "00", "SHZ"

try:
    client = Client(PROVIDER)
except:
    client = Client("RASPISHAKE")

BEGTIME = "2021-04-05T00:00:00.000"
BEGTIME = "2021-01-01T00:00:00.000"
DURATION_HOURS=24*100
DURATION=3600.0*DURATION_HOURS-1
WINLEN_SEC=3600
STEP_SEC=int(WINLEN_SEC)

FMIN,FMAX=3.0,15.0

t1 = UTCDateTime(BEGTIME)
t2=t1+DURATION

datelist = pd.date_range(t1.datetime, min(t2, UTCDateTime()).datetime, freq="D")
nslc = "{}.{}.{}.{}".format(NET, STA, LOC, CHAN)
# make sure that wildcard characters are not in nslc
nslc = nslc.replace("*", "").replace("?", "")

print ("")
pbar = tqdm.tqdm(datelist)
for day in pbar:
    datestr = day.strftime("%Y-%m-%d")
    fn = "./data/{}_{}.mseed".format(datestr, nslc)
    if day != UTCDateTime().datetime and os.path.isfile(fn):
        continue
    else:
        pbar.set_description("Fetching %s" % fn)
        try: 
            st = client.get_waveforms(NET, STA, LOC, CHAN, UTCDateTime(day), UTCDateTime(day)+86400, attach_response=True)
            st.merge(fill_value='interpolate')
            sr=st[0].stats.sampling_rate
            if (sr>=50.0):
                print ("Decimate...")
                st.decimate(2, strict_length=False, no_filter=True)
        except FDSNNoDataException:
            pbar.set_description("No data on FDSN server for %s" % fn)
            continue
        st.write(fn)

inv = client.get_stations(network=NET, station=STA, location=LOC, channel=CHAN, level='RESP')
print (inv[0][0][0].response)
#inv.plot_response(0.001, output="DISP")
#print ("Sensitivity:")
#for f in [1.0,3.0,5.0,15.0]:
    #print (inv[0][0][0].response._get_overall_sensitivity_and_gain(frequency=f, output='DISP'))

SENSIB_DISP=1.0/inv[0][0][0].response._get_overall_sensitivity_and_gain(frequency=5.0, output='DISP')[1]*1e+9
print ("SENSIB_DISP=",SENSIB_DISP)


print ("")
force_reprocess=0
pbar = tqdm.tqdm(datelist)
for day in pbar:
    datestr = day.strftime("%Y-%m-%d")
    fn_in = "./data/{}_{}.mseed".format(datestr, nslc)
    pbar.set_description("Removing response %s" % fn_in)
    if not os.path.isfile(fn_in):
        continue
    stall = read(fn_in, headonly=True)
    for mseedid in list(set([tr.id for tr in stall])):
        fn_out = "./data/{}_{}_disp.mseed".format(datestr, mseedid)
        if os.path.isfile(fn_out) and not force_reprocess:
            continue
        st = read(fn_in, sourcename=mseedid)
        stD = st.copy()
        #data=stD[0].data*SENSIB_DISP
        #stD[0].data=data
        stD.remove_response(output='DISP', inventory=inv)
        stD.write(fn_out, format='MSEED')

    del stall

print ("")
pbar = tqdm.tqdm(datelist)
for day in pbar:
    datestr = day.strftime("%Y-%m-%d")
    fn_in = "./data/{}_{}_disp.mseed".format(datestr, nslc)
    #fn_in = "./data/{}_{}.mseed".format(datestr, nslc)
    csvfile="./csv/%s_%s_%s_%s_%s.csv" % (datestr,NET,STA,LOC,CHAN)
    pbar.set_description("Reading displacement files %s" % fn_in)
    if not os.path.isfile(fn_in):
        continue
    stall = read(fn_in, headonly=True)
    for mseedid in list(set([tr.id for tr in stall])):
        stD = read(fn_in, sourcename=mseedid)

    #stD.merge()
    day2=day+1
    ##stD.trim(UTCDateTime(day),UTCDateTime(day2+))
    sr=stD[0].stats.sampling_rate

    print (stD)

    print ("filter...")
    #stD.detrend()
    stD.filter(type='bandpass',freqmin=FMIN, freqmax=FMAX)

    print ("Max displacement=",max(abs(stD[0].data))*1e+6," µm")
    
    #st.plot() 
    #stD.plot() 

    print ("Build dataframe ...")
    #d = {'rtime': stD[0].times('utcdatetime') , 'rms': stD[0].data**2}
    d = { 'rms': stD[0].data**2}
    df = pd.DataFrame(data=d)
    #df['Datetime']=stD[0].times('utcdatetime')
    #print (df)

    print ("Rolling mean ...")
    WINLEN=int(sr*WINLEN_SEC)
    df2=df.rolling(WINLEN, center=True, min_periods=int(WINLEN/2)).mean().apply(np.sqrt)
    print ("...")
    df2['Datetime']=stD[0].times('timestamp')
    #df2['Datetime']=stD[0].times('utcdatetime')
    #print (df2)

    STEP=int(sr*STEP_SEC)
    df2=df2[int(STEP/2)::STEP]
    print (df2.iloc[0]['Datetime'])
    print (pd.to_datetime(df2.iloc[0]['Datetime'], unit='s'))
    print (pd.to_datetime(df2.iloc[-1]['Datetime'], unit='s'))
    print ("len(df2)=",len(df2))


    """
    t3=t1+86400-3600
    datelist = pd.date_range(t1.datetime,t3.datetime , freq="60min")
    print (df2)
    print (datelist)
    print (len(df2), len(datelist))
    df2['Datetime']=datelist
    """

    #print (df2)

    df2.to_csv(csvfile, index=False)

    del(stall)
    del(stD)

