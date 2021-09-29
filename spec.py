#!/usr/bin/python3
# -*- coding: utf-8 -*-

from ppsd_utils import *

dirname= (os.path.dirname(__file__))
basename = os.path.basename(__file__)
print ("============================================")
print ("start ",basename)

from obspy.core.inventory.inventory import read_inventory
from obspy import UTCDateTime, read
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx
from matplotlib import cm
import matplotlib.pyplot as plt

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

cfgfile="%s/rms.cfg"%dirname
SDSROOT, DATADIR, XMLDIR, TELSITE_XMLDIR = readcfg(cfgfile)

try:
    begtime=argv[1]
    durationH=int(argv[2])
    network = argv[3]
    station=argv[4]
    location=argv[5]
    channel=argv[6]
    freqmax=float(argv[7])
except:
    raise()

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
#ppsd.plot(max_percentage=10, show=True, xaxis_frequency=True, cmap=pqlx, show_coverage=True, period_lim=(0.1,25), filename="./figures/{}_spectrum.png".format(nslc))
#fig=ppsd.plot(max_percentage=10, show=False, xaxis_frequency=True, cmap=pqlx, show_coverage=True, period_lim=(0.1,25))

from obspy.signal import spectral_estimation as spe
(nhnm_T, nhnm_psd) = spe.get_nhnm()
(nlnm_T, nlnm_psd) = spe.get_nlnm()
(nhnm_fq, nlnm_fq) = (np.reciprocal(nhnm_T), np.reciprocal(nlnm_T))

fig= plt.figure()
ax = fig.gca()
ax.semilogx(nhnm_fq, nhnm_psd, 'k--')
ax.semilogx(nlnm_fq, nlnm_psd, 'k--')

val50 = ppsd.get_percentile(50)
fq50 = np.reciprocal(val50[0])
amp50 = val50[1]

val5=ppsd.get_percentile(5)
val95=ppsd.get_percentile(95)
fq5, amp5 = np.reciprocal(val5[0]), val5[1]
fq95, amp95 = np.reciprocal(val95[0]), val95[1]

#Q=np.quantile(amp, 0.90)
#print (Q)
#amp=np.clip(amp, 0, Q)

ax.semilogx(fq50, amp50, 'ro', markersize=1, label='Median')
ax.semilogx(fq5, amp5, 'o', color='lightgrey', markersize=0.5, label='5th percentile')
ax.semilogx(fq95, amp95, 'o', color='grey', markersize=0.5, label='95th percentile')
#ax.semilogx(fq, ppsd, 'o', cmap='pqlx', markersize=0.1)

ax.set_ylim(-175, -85)
ax.set_xlim(0.1, freqmax)

ax.set(xlabel='Fréquence (Hz)')
ax.set(ylabel=r'PSD Acceleration 20log($m/s^2/\sqrt{Hz}$) dB')
ax.set_title("{} {} -> {}".format(nslc,begtime[:10],str(end)[:10]), fontsize=10)

plt.legend(loc='lower left')
plt.grid()

print (ax,ax.get_ylim(), ax.get_xlim())
#ax.set_xlim(0.1,25)
#ax.set_ylim(-175, -85)
#plt.show()
plt.savefig("{}/figures/{}_spectrum.png".format(dirname,nslc))

#[ppsd.plot_temporal(0.10) for mseedid, ppsd in ppsds.items()]
#ppsd.plot_spectrogram(clim=(-160,-100), cmap=cmap2, filename="./figures/{}_spectrogram.png".format(nslc))

