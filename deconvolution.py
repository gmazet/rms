#!/usr/bin/python3
# -*- coding: utf-8 -*-
from ppsd_utils import *

basename = os.path.basename(__file__)
print ("============================================")
print ("start ",basename)

# ----------------------------------
force_reprocess=0
#OUTPUT="VEL" # DISP, VEL or ACC
OUTPUT="DISP" # DISP, VEL or ACC
# ----------------------------------

print ("")

try:
    begtime=argv[1]
    durationH=int(argv[2])
    network = argv[3]
    station=argv[4]
    location=argv[5]
    channel=argv[6]
    OUTPUT=argv[7]
except:
    print ('Usage: %s "YYYY-MM-DDTHH:mm:ss.s" <nb hours> <network> <station> <location> <channel> <DISP or VEL or ACC>' % basename)
    exit()

PROCESSED_DATA_DIR="%s/%s" % (DATADIR,OUTPUT)

start = UTCDateTime(begtime)
end=start+durationH*3600
print ("start=",start,"end=",end)
#pre_filt = [0.05, 0.1, 30, 35]
#pre_filt = [0.5, 0.7, 30, 32]
pre_filt = [0.01, 0.05, 28, 30]

# ---------------------------------------------
datelist = pd.date_range(start.datetime, min(end, UTCDateTime()).datetime, freq="D")
nslc = "{}.{}.{}.{}".format(network, station, location, channel)
# make sure that wildcard characters are not in nslc
nslc = nslc.replace("*", "").replace("?", "")

print ("")
client=SDSClient(SDSROOT)
try:
    fclient=FDSNClient("RESIF")
except:
    pass
pbar = tqdm.tqdm(datelist)
for day in pbar:
    YYYY = day.strftime("%Y")
    MM = day.strftime("%m")
    DD = day.strftime("%d")
    datestr = day.strftime("%Y-%m-%d")
    outdir="{}/rawdata/{}/{}/{}".format(DATADIR,YYYY,MM,DD)
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
    fn = "{}/{}_{}.mseed".format(outdir,datestr, nslc)
    if day != UTCDateTime().datetime and os.path.isfile(fn):
        continue
    else:
        pbar.set_description("Fetching %s" % fn)
        try:
            st = client.get_waveforms(network, station, location, channel, UTCDateTime(day)-1800, UTCDateTime(day)+86400+1800, attach_response=True)
            st.merge(fill_value='interpolate')
            #except FDSNNoDataException:
        except:
            print ("Can't find data in SDS. Try FDSN request")
            pass

        if (not st):
            try:
                st = fclient.get_waveforms(network, station, location, channel, UTCDateTime(day)-1800, UTCDateTime(day)+86400+1800, attach_response=True)
                st.merge(fill_value='interpolate')
            except FDSNNoDataException:
                pbar.set_description("No data on FDSN server for %s" % fn)
                continue

        print (st)
        sr=st[0].stats.sampling_rate
        #if (sr>=100.0):
            #print ("Decimate...")
            #st.decimate(2, strict_length=False, no_filter=True)
        #elif (sr>=50.0):
            #print ("Decimate...")
            #st.decimate(2, strict_length=False, no_filter=True)
        st.write(fn)


try:
    inv=read_inventory('%s/%s.%s.%s.%s.xml' % (XMLDIR,network,station,location, channel))
except:
    try:
        inv=read_inventory('%s/%s.%s.%s.%s.xml' % (TELSITE_XMLDIR,network,station,location, channel))
    except:
        try:
            inv=fclient.get_stations(network=network, sta=station, loc=location, channel=channel, level="response")
        except:
            raise()

try:
    print (inv[0][0][0].response)
    SENSIB_DISP=1.0/inv[0][0][0].response._get_overall_sensitivity_and_gain(frequency=5.0, output='DISP')[1]*1e+9
    SENSIB_VEL=1.0/inv[0][0][0].response._get_overall_sensitivity_and_gain(frequency=5.0, output='VEL')[1]*1e+9
    print ("SENSIB_DISP= %.3f"%SENSIB_DISP," nm/count at 5.0Hz")
    print ("SENSIB_VEL = %.3f"%SENSIB_VEL," nm/s/count at 5.0Hz")

    SENSIB=1.0/inv[0][0][0].response._get_overall_sensitivity_and_gain(frequency=5.0, output=OUTPUT)[1]
except:
    pass

print ("")
pbar = tqdm.tqdm(datelist)
for day in pbar:
    YYYY = day.strftime("%Y")
    MM = day.strftime("%m")
    DD = day.strftime("%d")
    datestr = day.strftime("%Y-%m-%d")
    rawdatadir="{}/rawdata/{}/{}/{}".format(DATADIR,YYYY,MM,DD)

    fn_in = "{}/{}_{}.mseed".format(rawdatadir,datestr, nslc)
    pbar.set_description("Removing response %s" % fn_in)
    if not os.path.isfile(fn_in):
        continue
    stall = read(fn_in, headonly=True)
    for mseedid in list(set([tr.id for tr in stall])):
        #fn_out = "{}/{}_{}.mseed".format(PROCESSED_DATA_DIR,datestr, mseedid)
        processeddatadir="{}/{}/{}/{}/{}".format(DATADIR,OUTPUT,YYYY,MM,DD)
        pathlib.Path(processeddatadir).mkdir(parents=True, exist_ok=True)
        #fn_out = "{}/{}_{}.mseed".format(PROCESSED_DATA_DIR,datestr, mseedid)
        fn_out = "{}/{}_{}.mseed".format(processeddatadir,datestr, mseedid)
        if os.path.isfile(fn_out) and not force_reprocess:
            continue
        st = read(fn_in, sourcename=mseedid)
        stD = st.copy()
        ##data=stD[0].data*SENSIB
        ##stD[0].data=data
        ###stD.remove_response(output='DISP', inventory=inv)
        stD.remove_response(inventory=inv,pre_filt=pre_filt,zero_mean=True,taper=True, output=OUTPUT,water_level=60,plot=False)
        stD.write(fn_out, format='MSEED')

    del stall

