#!/usr/bin/python3

# -*- coding: utf-8 -*-
from obspy.geodetics import (gps2dist_azimuth, degrees2kilometers)
from obspy.clients.filesystem.sds import Client as SDSClient
from obspy.clients.fdsn import Client as FDSNClient
from obspy import UTCDateTime, read
from obspy.clients.fdsn.client import FDSNNoDataException
from obspy.core.inventory.inventory import read_inventory

from datetime import timedelta
import tqdm
import warnings
import numpy as np
import pandas as pd
from glob import glob

from sys import exit,argv
import os

SDSROOT="/cea/dsku/SDS_buffer_AF"

XMLDIR="/cea/data/outils/dase/seiscomp/xml"
TELSITE_XMLDIR="/cea/data/outils/dase/noise_analysis/telsite/invfiles"

DATADIR="/cea/dsku/SDS_muru/SDS_muru/rms"

font = {'family': 'Tahoma',
        'color':  'black',
        'weight': 'normal',
        'size': 10,
        }

POLES_AND_ZEROS={
'L22':[[-8.886 + 8.886j, -8.886 - 8.886j],[0j, 0j]],
'T120':[[-0.03859+ 0.03649j, -0.03859- 0.03649j, -32.55+0j, -142+0j, -364+404j, -364-404j, -1260+0j, -4900+5200j, -4900-5200j, -7100+1700j, -7100-1700j], [0j, 0j, -31.63+0j, -160+0j ,-350+0j, -3177+0j]],
'ZM500':[[-4.44 - 4.44j, -4.44 + 4.44j], [0j, 0j]]
}

paz_L22 = {
    'poles': POLES_AND_ZEROS['L22'][0],
    'zeros': POLES_AND_ZEROS['L22'][1],
    'gain': 1.0,
    'sensitivity': 0.3341e+9} # nm
    #'sensitivity': 2.099e+9} # nm/s
paz_T120 = {
    'poles': POLES_AND_ZEROS['T120'][0],
    'zeros': POLES_AND_ZEROS['T120'][1],
    'gain': 8.31871e+17,
    'sensitivity': 0.3118e+9} # nm
    #'sensitivity': 5.10399e+08} # nm
    #'sensitivity': 1.9593e+9} # nm/s


