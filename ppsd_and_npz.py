#!/usr/bin/python3
# -*- coding: utf-8 -*-

from ppsd_utils import *
from obspy.signal import PPSD

basename = os.path.basename(__file__)
print ("============================================")
print ("start ",basename)

# -----------------------------------
# For remove response
#pre_filt = [0.05, 0.1, 30, 35]
#pre_filt = [0.5, 0.7, 30, 32]
pre_filt = [0.01, 0.05, 28, 30]

# PPSD
freqmin=0.02 # seconds
freqmax=10.0 # seconds

force_reprocess=0
# -----------------------------------

try:
    begtime=argv[1]
    durationH=int(argv[2])
    network = argv[3]
    station=argv[4]
    location=argv[5]
    channel=argv[6]
except:
    print ('Usage: %s "YYYY-MM-DDTHH:mm:ss.s" <nb hours> <network> <station> <location> <channel>')
    exit()

start = UTCDateTime(begtime)
end=start+durationH*3600
print ("start=",start,"end=",end)

cfgfile="./rms.cfg"
SDSROOT, DATADIR, XMLDIR, TELSITE_XMLDIR = readcfg(cfgfile)

try:
    fclient=FDSNClient("RESIF")
except:
    pass

# ---------------------------------------------

datelist = pd.date_range(start.datetime, min(end, UTCDateTime()).datetime, freq="D")
nslc = "{}.{}.{}.{}".format(network, station, location, channel)
nslc = nslc.replace("*", "").replace("?", "")

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

#if (nslc=="RD.SFTF..HHZ"):
#    print ("Use paz files instead of xml")
#    inv=paz_T120
#    print(inv)

"""
try:
    print (inv[0][0][0].response)
    SENSIB_DISP=1.0/inv[0][0][0].response._get_overall_sensitivity_and_gain(frequency=5.0, output='DISP')[1]*1e+9
    SENSIB_VEL=1.0/inv[0][0][0].response._get_overall_sensitivity_and_gain(frequency=5.0, output='VEL')[1]*1e+9
    print ("SENSIB_DISP= %.3f"%SENSIB_DISP," nm/count at 5.0Hz")
    print ("SENSIB_VEL = %.3f"%SENSIB_VEL," nm/s/count at 5.0Hz")

    SENSIB=1.0/inv[0][0][0].response._get_overall_sensitivity_and_gain(frequency=5.0, output=OUTPUT)[1]
except:
    pass
"""


# PPSD
print ("")
pbar = tqdm.tqdm(datelist)
for day in pbar:
    print ("day",day)

    YYYY = day.strftime("%Y")
    MM = day.strftime("%m")
    DD = day.strftime("%d")
    datestr = day.strftime("%Y-%m-%d")
    rawdatadir="{}/rawdata/{}/{}/{}".format(DATADIR,YYYY,MM,DD)

    fn_in = "{}/{}_{}.mseed".format(rawdatadir,datestr, nslc)
    pbar.set_description("PPSD %s" % fn_in)
    if not os.path.isfile(fn_in):
        continue
    stall = read(fn_in, headonly=True)
    for mseedid in list(set([tr.id for tr in stall])):
        npzdatadir="{}/npz/{}/{}/{}".format(DATADIR,YYYY,MM,DD)
        pathlib.Path(npzdatadir).mkdir(parents=True, exist_ok=True)
        #fn_out = "{}/npz/{}_{}.npz".format(DATADIR,datestr, nslc)
        fn_out = "{}/{}_{}.npz".format(npzdatadir,datestr, nslc)

        if os.path.isfile(fn_out) and not force_reprocess:
            continue
        st = read(fn_in, sourcename=mseedid)
        trace=st[0]
        print ("ppsd...")
        ppsd = PPSD(trace.stats, inv,
                ppsd_length=120, overlap=0.5,
                period_smoothing_width_octaves=0.025,
                period_step_octaves=0.0125,
                period_limits=(freqmin, freqmax),
                db_bins=(-200, 20, 0.25))
                #period_limits=(0.008, 50),
        print ("add traces ...")
        ppsd.add(st)
        ppsd.save_npz(fn_out[:-4])
        del st, ppsd
    del stall

