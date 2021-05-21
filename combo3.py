#!/usr/bin/python3

# -*- coding: utf-8 -*-
from ppsd_utils import *
from obspy.signal import PPSD

basename = os.path.basename(__file__)
print ("============================================")
print ("start ",basename)

import matplotlib
matplotlib.rcParams['agg.path.chunksize'] = 10000
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import cm

from obspy.imaging.cm import pqlx, obspy_sequential
from obspy.imaging.util import _set_xaxis_obspy_dates

cfgfile="./rms.cfg"
SDSROOT, DATADIR, XMLDIR, TELSITE_XMLDIR = readcfg(cfgfile)

# ---------------------------------
#OUTPUT="VEL"
#YMAX=25.0 # nm/s in rms 

OUTPUT="DISP"
YMAX1=7.0 # nm in signal 
YMAX2=0.8 # nm in rms 
cmap2 = cm.inferno

DECIM=1 # 6 = every 1 hour; 3 = every 1/2 hour; 1 = every 10 min
STEP=DECIM
QUANTILE=0.995

per_lap =0.9
#clip=[0.1,1.0]
clip=[0.6,0.99] # dbscale=True
clip=[0.01,0.5] # dbscale=False

specfreqmax=20.0
specfreqmin=0.5
# ---------------------------------

def plot_spectrogram_and_rms(ppsd, trace, df, cmap=cmap2, clim=None, grid=True,
                         filename=None, show=True):
        """
        Plot the temporal evolution of the PSD in a spectrogram-like plot.

        .. note::
            For example plots see the :ref:`Obspy Gallery <gallery>`.

        :type cmap: :class:`matplotlib.colors.Colormap`
        :param cmap: Specify a custom colormap instance. If not specified, then
            the default ObsPy sequential colormap is used.
        :type clim: list
        :param clim: Minimum/maximum dB values for lower/upper end of colormap.
            Specified as type ``float`` or ``None`` for no clipping on one end
            of the scale (e.g. ``clim=[-150, None]`` for a lower limit of
            ``-150`` dB and no clipping on upper end).
        :type grid: bool
        :param grid: Enable/disable grid in histogram plot.
        :type filename: str
        :param filename: Name of output file
        :type show: bool
        :param show: Enable/disable immediately showing the plot.
        """
        import matplotlib.pyplot as plt

        #fig, ax = plt.subplots()
        print ("Plot spectro...")
        #fig, axes = plt.subplots(2)
        #fig.set_size_inches(8, 5, forward=True)
        fig=plt.figure(figsize=[8,6])
        ax11 = fig.add_axes([0.1, 0.05, 0.85, 0.16]) # rms
        ax12 = fig.add_axes([0.1, 0.21, 0.85, 0.16]) # rms
        ax2 = fig.add_axes([0.1, 0.42, 0.85, 0.53]) # spectro

        #ax=ax2

        quadmeshes = []
        yedges = ppsd.period_xedges
        # Test GMR pour mettre le spectro en frequences...
        yedges = [1.0/per for per in ppsd.period_xedges]
        #yedges.reverse()
        #print (yedges[0], yedges[-1])
        ###exit()

        for times, psds in ppsd._get_gapless_psd():
            #print ("=======",times, ppsd.step)
            xedges = [t.matplotlib_date for t in times] + \
                [(times[-1] + ppsd.step).matplotlib_date]
            meshgrid_x, meshgrid_y = np.meshgrid(xedges, yedges)
            data = np.array(psds).T

            quadmesh = ax2.pcolormesh(meshgrid_x, meshgrid_y, data, cmap=cmap,
                                     zorder=-1)
            quadmeshes.append(quadmesh)

        if clim is None:
            cmin = min(qm.get_clim()[0] for qm in quadmeshes)
            cmax = max(qm.get_clim()[1] for qm in quadmeshes)
            clim = (cmin, cmax)

        for quadmesh in quadmeshes:
            quadmesh.set_clim(*clim)

        if grid:
            ax2.grid()

        #cb = plt.colorbar(quadmesh, ax=ax2)
        #cb.set_label('Amplitude [$m^2/s^4/Hz$] [dB]')

        #if ppsd.special_handling is None:
            #cb.ax2.set_ylabel('Amplitude [$m^2/s^4/Hz$] [dB]')
        #else:
            #cb.ax2.set_ylabel('Amplitude [dB]')
        ax2.set_ylabel('Fréquence [Hz]')

        fig.autofmt_xdate()
        _set_xaxis_obspy_dates(ax2)

        ax2.set_yscale("log")
        ax2.set_xlim(ppsd.times_processed[0].matplotlib_date,
                    (ppsd.times_processed[-1] + ppsd.step).matplotlib_date)
        #ax2.set_ylim(yedges[0], yedges[-1])
        try:
            ax2.set_facecolor('0.8')
        # mpl <2 has different API for setting Axes background color
        except AttributeError:
            ax2.set_axis_bgcolor('0.8')

        # Plot trace
        print ("Plot trace...", FMIN,FMAX)
        trace.filter(type='bandpass',freqmin=FMIN, freqmax=FMAX, corners=4)
        npts = trace.stats.npts
        dt = trace.stats.delta
        t = np.linspace(0, dt * npts, npts)
        t1,t2=trace.times("matplotlib")[0], trace.times("matplotlib")[-1]
        #trace.plot(automerge=False,fig=fig, marker=',', color='k', equal_scale=False)
        ax11.plot(trace.times("matplotlib"), trace.data, 'k', linewidth=0.4)
        #print (trace)
        print ("Max trace:",max(abs(trace.data)))
        #print (t)
        ax11.set_xlim(t1, t2)


        print ("Plot RMS...")
        ##ax11.plot(df['Datetime'], df['rms'], 'b-', linewidth=1.0)
        ax12.plot(df['Datetime'], df['rms'], 'b-', linewidth=0.8)

        _set_xaxis_obspy_dates(ax11)

        ax11.set_xlim(ppsd.times_processed[0].matplotlib_date,
                    (ppsd.times_processed[-1] + ppsd.step).matplotlib_date)
        ax12.set_xlim(ppsd.times_processed[0].matplotlib_date,
                    (ppsd.times_processed[-1] + ppsd.step).matplotlib_date)

        Q1=np.quantile(abs(trace.data), 0.98)
        #YMAX11=min (max(max(abs(trace.data))*1.1, YMAX1), 50)
        YMAX11=min (max(Q1*1.2, YMAX1), 50)
        ax11.set_ylim(-YMAX11,YMAX11) # signal

        Q2=np.quantile(df.rms*1.1, 0.9)
        YMAX22=max(Q2, YMAX2)
        ax12.set_ylim(0,YMAX22) # rms

        #YMAX11=max(abs(df['rms']))*1.05
        #if (station != "HYF"):
        #    ax11.set_ylim(-YMAX,YMAX) # nm/s
        #    ax12.set_ylim(-YMAX,YMAX) # nm/s
        #else:
        #    ax11.set_ylim(-YMAX11,YMAX11) # nm/s
        #    ax12.set_ylim(-YMAX11,YMAX11) # nm/s

        #ax11.set_yscale("log")
        #ax12.set_yscale("log")

        #bbox = {'fc': '0.8', 'pad': 0}
        props={'ha': 'left', 'va': 'center'}
        ax11.text(ax11.get_xlim()[1], ((ax11.get_ylim()[0]+ax11.get_ylim()[1])/2), FREQ1.replace("_","-").replace("Hz"," Hz"), props, rotation=270, color='r')
        ax12.text(ax12.get_xlim()[1], ((ax12.get_ylim()[0]+ax12.get_ylim()[1])/2), FREQ2.replace("_","-").replace("Hz"," Hz"), props, rotation=270, color='r')

        if (OUTPUT=="DISP"):
            ax11.set_ylabel('Displ (nm)')
            ax12.set_ylabel('RMS (nm)')
        if (OUTPUT=="VEL"):
            ax11.set_ylabel('Vel (nm/s)')
            ax12.set_ylabel('RMS (nm/s)')

        ax11.grid()
        ax12.xaxis.set_visible(False)
        ax12.grid()

        ax2.set_ylim(specfreqmax, specfreqmin) ## 1Hz max, pour ne pas voir le pic de houle
        ax2.set_title("%s" % (nslc))
        #ax11.set_title("Filtre %.1f-%.1f Hz" % (FMIN,FMAX),fontsize=10)


        #fig.tight_layout()
        #fig.subplots_adjust(top=0.95)
        ax2.invert_yaxis()

        if filename is not None:
            plt.savefig(filename)
            plt.close()
        elif show:
            plt.draw()
            plt.show()
        else:
            plt.draw()
            return fig


try:
    begtime=argv[1]
    durationH=int(argv[2])
    network = argv[3]
    station=argv[4]
    location=argv[5]
    channel=argv[6]
    FMIN=float(argv[7])
    FMAX=float(argv[8])
    FREQ1=argv[9]
    FREQ2=argv[10]
except:
    begtime = "2021-03-08T00:00:00.000"
    durationH=96
    network = "RD"
    station = "CABF"
    location = ""
    channel = "SHZ"
    FMIN,FMAX=3.0,15.0
    FREQ1="2.0_10.0Hz"
    FREQ2="1.0_6.0Hz"
    FREQ2=FREQ1

start = UTCDateTime(begtime)
DURATION=durationH*3600
end=start+DURATION
wlen=int(end-start)/60 # pour spectrogramme

t3=start+DURATION/100
t4=end-DURATION/100

XMLDIR="/dase/AF-OP/gm171646/Documents/LDG/Etudes/Eoliennes/signaux/xmlresponses"

SDSROOT="/cea/dsku/SDS_buffer_AF"
#client=SDSClient(SDSROOT)
#fclient=FDSNClient("RESIF")

##DATADIR="/cea/dsku/SDS_muru/SDS_muru/rms/%s"%OUTPUT
DATADIR="/cea/dsku/SDS_muru/SDS_muru/rms"

# ----------------------------
#OUTPUT="VEL" # DISP, VEL or ACC
OUTPUT="DISP" # DISP, VEL or ACC
PROCESSED_DATA_DIR="%s/%s" % (DATADIR,OUTPUT)
# ----------------------------


nslc = "{}.{}.{}.{}".format(network, station, location, channel)
nslc = nslc.replace("*", "").replace("?", "")

datelist = pd.date_range(start.datetime, min(end, UTCDateTime()).datetime, freq="D")

pbar = tqdm.tqdm(datelist)
for day in pbar:
    day2=(day+timedelta(days=1))
    YYYY = day.strftime("%Y")
    MM = day.strftime("%m")
    DD = day.strftime("%d")
    datestr = day.strftime("%Y-%m-%d")
    datadir="{}/{}/{}/{}/{}".format(DATADIR,OUTPUT,YYYY,MM,DD)

    fn_in = "{}/{}_{}.mseed".format(datadir,datestr, nslc)
    pbar.set_description("Reading deconvolved data %s" % fn_in)
    if not os.path.isfile(fn_in):
        continue
    stall = read(fn_in, headonly=True)
    for mseedid in list(set([tr.id for tr in stall])):
        st = read(fn_in, sourcename=mseedid)

    sr=st[0].stats.sampling_rate
    st.trim(UTCDateTime(day.to_pydatetime()), UTCDateTime(day2.to_pydatetime()))
    if (sr>=100.0):
        print ("Decimate...")
        st.decimate(4, strict_length=False, no_filter=True)
    elif (sr>=50.0):
        print ("Decimate...")
        st.decimate(2, strict_length=False, no_filter=True)

    st[0].data=st[0].data*1e+9 # m -> nm

    try:
        stALL+=st
    except:
        stALL=st

    del(st)

stALL.merge(fill_value='interpolate')
print (stALL)
#print (max(stALL[0].data))

trace=stALL[0]


# ---------------------------------
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

# ---------------------------------------------
# Reead rms files
# ---------------------------------------------
headers = ['rms', 'Datetime']
dtypes = {'rms': 'float', 'Datetime': 'str'}
parse_dates = ['Datetime']

print ("%s/csv/*_%s_%s.csv" % (DATADIR,nslc,FREQ1))
#for fn in glob("%s/csv/*_%s_%s.csv" % (DATADIR,nslc,FREQ1)):

pbar = tqdm.tqdm(datelist)
for day in pbar:
    YYYY = day.strftime("%Y")
    MM = day.strftime("%m")
    DD = day.strftime("%d")
    datestr = day.strftime("%Y-%m-%d")

    csvdatadir="{}/csv/{}/{}/{}".format(DATADIR,YYYY,MM,DD)
    fn_pattern = "{}/{}_{}_{}.csv".format(csvdatadir,datestr, nslc,FREQ1)
    pbar.set_description("Reading %s" % fn_pattern)
    for fn in glob(fn_pattern):
        df=pd.read_csv(fn, sep=',', header=0, names=headers, dtype=dtypes, parse_dates=parse_dates)
        Q=np.quantile(df.rms, QUANTILE)
        df.rms=np.clip(df.rms, 0, Q)
        df.rms=df.rms*1e+9 #nm
        try:
            dfALL=pd.concat([dfALL, df], ignore_index=True)
        except:
            dfALL=df

try:
    dfALL['Datetime']=pd.to_datetime(dfALL['Datetime'], unit='s')
    dfALL=dfALL.sort_values(by=['Datetime'])
except:
    raise()

dfRMS1=dfALL.rolling(DECIM, center=True, min_periods=STEP).mean() # to get rid of calibration/earthquakes effects
dfRMS1['Datetime']=dfALL['Datetime']
dfRMS1['Datetime']=pd.to_datetime(dfRMS1['Datetime'], unit='s')

del dfALL
del df

print ("----")

"""
# ---------------------------------------------
# Reead rms files
# ---------------------------------------------
for fn in glob("%s/csv/*_%s_%s.csv" % (DATADIR,nslc,FREQ2)):
    df=pd.read_csv(fn, sep=',', header=0, names=headers, dtype=dtypes, parse_dates=parse_dates)
    Q=np.quantile(df.rms, QUANTILE)
    df.rms=np.clip(df.rms, 0, Q)
    df.rms=df.rms*1e+9 #nm
    try:
        dfALL=pd.concat([dfALL, df], ignore_index=True)
    except:
        dfALL=df

try:
    dfALL['Datetime']=pd.to_datetime(dfALL['Datetime'], unit='s')
    dfALL=dfALL.sort_values(by=['Datetime'])
except:
    raise()

dfRMS2=dfALL.rolling(DECIM, center=True, min_periods=STEP).mean() # to get rid of calibration/earthquakes effects
dfRMS2['Datetime']=dfALL['Datetime']
dfRMS2['Datetime']=pd.to_datetime(dfRMS2['Datetime'], unit='s')
dfRMS2=dfRMS2.rename(columns = {'rms':'rms2'})

dfRMS3=pd.merge(dfRMS2, dfRMS1, on=["Datetime"])

# ---

print ("====")
dfRMS3=dfRMS3[::STEP]
print (dfRMS3)

"""

# ---------------------------------
for mseedid, ppsd in ppsds.items():
    print (mseedid, nslc)

#print ("Max=",max(tr.data))

plot_spectrogram_and_rms(ppsd, trace, dfRMS1, filename="./figures/%s_combo_signal.png"%nslc, show=False,clim=[-160, -120])
#plt.show()
exit()



