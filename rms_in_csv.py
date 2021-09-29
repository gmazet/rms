#!/usr/bin/python3
# -*- coding: utf-8 -*-
from ppsd_utils import *
from obspy.signal import PPSD

dirname= (os.path.dirname(__file__))
basename = os.path.basename(__file__)
print ("============================================")
print ("start ",basename)

# ----------------------------
#OUTPUT="VEL" # DISP, VEL or ACC
#OUTPUT="DISP" # DISP, VEL or ACC
#PROCESSED_DATA_DIR="%s/%s" % (DATADIR,OUTPUT)
# ----------------------------

print ("")

try:
    begtime=argv[1]
    durationH=int(argv[2])
    network = argv[3]
    station=argv[4]
    location=argv[5]
    channel=argv[6]
    FMIN=float(argv[7])
    FMAX=float(argv[8])
    winlen_sec=int(argv[9])
    OUTPUT=argv[10]
except:
    print ('Usage: %s "YYYY-MM-DDTHH:mm:ss.s" <nb hours> <network> <station> <location> <channel> <freq min> <freq max> <window length (seconds)>')
    exit()


cfgfile="%s/rms.cfg"%dirname
SDSROOT, DATADIR, XMLDIR, TELSITE_XMLDIR = readcfg(cfgfile)

start = UTCDateTime(begtime)
end=start+durationH*3600
print ("start=",start,"end=",end)
step_sec=int(winlen_sec/2.0)

# ---------------------------------------------

datelist = pd.date_range(start.datetime, min(end, UTCDateTime()).datetime, freq="D")
nslc = "{}.{}.{}.{}".format(network, station, location, channel)
# make sure that wildcard characters are not in nslc
nslc = nslc.replace("*", "").replace("?", "")

# ---------------------------
print ("")
pbar = tqdm.tqdm(datelist)
for day in pbar:
    print ("day",day)
    YYYY = day.strftime("%Y")
    MM = day.strftime("%m")
    DD = day.strftime("%d")
    datestr = day.strftime("%Y-%m-%d")
    processeddatadir="{}/{}/{}/{}/{}".format(DATADIR,OUTPUT,YYYY,MM,DD)
    fn_in = "{}/{}_{}.mseed".format(processeddatadir,datestr, nslc)
    #fn_in = "{}/{}_{}.mseed".format(PROCESSED_DATA_DIR,datestr, nslc)
    #csvfile="%s/csv/%s_%s_%.1f_%.1fHz.csv" % (DATADIR,datestr,nslc,FMIN,FMAX)
    csvdatadir="{}/csv/{}/{}/{}".format(DATADIR,YYYY,MM,DD)
    pathlib.Path(csvdatadir).mkdir(parents=True, exist_ok=True)
    csvfile="%s/%s_%s_%.1f_%.1fHz.csv" % (csvdatadir,datestr,nslc,FMIN,FMAX)
    pbar.set_description("Reading displacement (or velocity) files %s" % fn_in)
    if not os.path.isfile(fn_in):
        continue
    stall = read(fn_in, headonly=True)
    for mseedid in list(set([tr.id for tr in stall])):
        stD = read(fn_in, sourcename=mseedid)

    #stD.merge()
    day2=(day+timedelta(days=1))

    sr=stD[0].stats.sampling_rate

    print ("filter...")
    #stD.detrend()
    stD.filter(type='bandpass',freqmin=FMIN, freqmax=FMAX, corners=4)
    stD.trim(UTCDateTime(day)-600,UTCDateTime(day2)+600)
    print (stD)

    if (OUTPUT=="DISP"):
        print ("Max displacement=",max(abs(stD[0].data))*1e+9," nm")
    if (OUTPUT=="VEL"):
        print ("Max velocity=",max(abs(stD[0].data))*1e+9," nm/s")
    if (OUTPUT=="ACC"):
        print ("Max acceleration=",max(abs(stD[0].data))*1e+9," nm/s/s")

    print ("Build dataframe ...")
    #d = {'rtime': stD[0].times('utcdatetime') , 'rms': stD[0].data**2}
    d = { 'rms': stD[0].data**2}
    df = pd.DataFrame(data=d)
    #df['Datetime']=stD[0].times('utcdatetime')
    #print (df)

    print ("Rolling rms ...")
    winlen_samp=int(sr*winlen_sec)
    df2=df.rolling(winlen_samp, center=True, min_periods=int(winlen_samp/2)).mean().apply(np.sqrt)
    print ("...")
    df2['Datetime']=stD[0].times('timestamp')
    #df2['Datetime']=stD[0].times('utcdatetime')
    #print (df2)

    step=int(sr*step_sec)
    #df2=df2[int(step/2)::step]
    df2=df2[step::step]
    #print (df2.iloc[0]['Datetime'])
    #print (pd.to_datetime(df2.iloc[0]['Datetime'], unit='s'))
    #print (pd.to_datetime(df2.iloc[-1]['Datetime'], unit='s'))
    #print ("len(df2)=",len(df2))

    df2.to_csv(csvfile, index=False)

    del(stall)
    del(stD)


