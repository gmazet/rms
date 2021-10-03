#!/usr/bin/python3
# -*- coding: utf-8 -*-
from ppsd_utils import *

basename = os.path.basename(__file__)
print ("============================================")
print ("start ",basename)

# ----------------------------------
force_reprocess=0
# ----------------------------------

print ("")

try:
    begtime=argv[1]
    durationH=int(argv[2])
    network = argv[3]
    station=argv[4]
    location=argv[5]
    channel=argv[6]
except:
    print ('Usage: %s "YYYY-MM-DDTHH:mm:ss.s" <nb hours> <network> <station> <location> <channel> <DISP or VEL or ACC>' % basename)
    exit()

cfgfile="./rms.cfg"
SDSROOT, DATADIR, XMLDIR, TELSITE_XMLDIR = readcfg(cfgfile)

print (SDSROOT)

start = UTCDateTime(begtime)
end=start+durationH*3600
print ("start=",start,"end=",end)

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
    print ("ERROR: Can't connect FDSN client to RESIF WS")
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

