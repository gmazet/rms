#!/usr/bin/python3
# -*- coding: utf-8 -*-

from obspy.core.inventory.inventory import read_inventory
from obspy import UTCDateTime, read
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx
from matplotlib import cm

from sys import exit,argv
import tqdm
import warnings

import pandas as pd
import numpy as np
from glob import glob

cmap1 = cm.magma
cmap1 = cm.plasma
cmap2 = cm.inferno

#SDSROOT="/cea/dsku/SDS_buffer_AF"
#XMLDIR="/dase/AF-OP/gm171646/Documents/LDG/Etudes/Eoliennes/signaux/xmlresponses"
#DATADIR="/cea/dsku/SDS_muru/SDS_muru/rms"

cfgfile="./rms.cfg"
SDSROOT, DATADIR, XMLDIR, TELSITE_XMLDIR = readcfg(cfgfile)

try:
    begtime=argv[1]
    durationH=int(argv[2])
    network = argv[3]
    station=argv[4]
    location=argv[5]
    channel=argv[6]
except:
    begtime = "2021-03-08T00:00:00.000"
    durationH=96
    network = "RD"
    station = "CABF"
    location = ""
    channel = "SHZ"

start = UTCDateTime(begtime)
end=start+durationH*3600
print ("start=",start,"end=",end)

# ---------------------------------------------

datelist = pd.date_range(start.datetime, min(end, UTCDateTime()).datetime, freq="D")
nslc = "{}.{}.{}.{}".format(network, station, location, channel)
# make sure that wildcard characters are not in nslc
nslc = nslc.replace("*", "").replace("?", "")

ppsds = {}
pbar = tqdm.tqdm(datelist)
for day in pbar:
    YYYY = day.strftime("%Y")
    MM = day.strftime("%m")
    DD = day.strftime("%d")
    datestr = day.strftime("%Y-%m-%d")
    npzdatadir="{}/npz/{}/{}/{}".format(DATADIR,YYYY,MM,DD)
    fn_pattern = "{}/{}_{}.npz".format(npzdatadir,datestr, nslc)
    pbar.set_description("Reading %s" % fn_pattern)
    for fn in glob(fn_pattern):
        mseedid = fn.replace(".npz", "").split("_")[-1]
        if mseedid not in ppsds:
            ppsds[mseedid] = PPSD.load_npz(fn)#, allow_pickle=True)
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ppsds[mseedid].add_npz(fn)#, allow_pickle=True)

#print (ppsds.items())
for mseedid, ppsd in ppsds.items():
    print (mseedid, nslc)
#ppsd.plot(max_percentage=5, show=True, xaxis_frequency=True, cmap=cmap1, show_coverage=True, period_lim=(0.1,25), filename="./figures/{}_spectrum.png".format(nslc))
ppsd.plot(max_percentage=10, show=True, xaxis_frequency=True, cmap=pqlx, show_coverage=True, period_lim=(0.1,25))


#[ppsd.plot_temporal(0.10) for mseedid, ppsd in ppsds.items()]
#ppsd.plot_spectrogram(clim=(-160,-100), cmap=cmap2, filename="./figures/{}_spectrogram.png".format(nslc))

