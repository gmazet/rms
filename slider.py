#!/usr/bin/python3
# -*- coding: utf-8 -*-
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

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
            #Q=np.quantile(df.rms, 0.95)
            #df.rms=np.clip(df.rms, 0, Q)
        except:
            raise()
    try:
        dfALL=pd.concat([dfALL, df], ignore_index=True)
    except:
        dfALL=df

dfALL['Datetime']=pd.to_datetime(dfALL['Datetime'], unit='s')
dfALL=dfALL.sort_values(by=['Datetime'])


# Load data
#df = pd.read_csv( "https://raw.githubusercontent.com/plotly/datasets/master/finance-charts-apple.csv")
#df = pd.read_csv( "./csv/2021-04-10_AM_R9F1B_00_SHZ.csv", sep=',')
#df.columns = [col.replace("AAPL.", "") for col in df.columns]
#df['dt']=df.Date + 'T' +  df['Time UTC']
#print df['dt']
#df['Datetime']=pd.to_datetime(df['Datetime'], unit='s')
print (dfALL.columns)

# Create figure
fig = go.Figure()
#fig.add_trace( go.Scatter(x=list(df.Datetime), y=list(df.rms), mode='markers') )
#fig=px.histogram(df,x="Magnitude", nbins=5)
fig=px.line(dfALL,x="Datetime", y="rms")

# Set title
fig.update_layout( title_text="RMS")

"""
# Add range slider
fig.update_layout(
    xaxis=dict(
        rangeselector=dict(
            buttons=list([
                dict(count=1,
                     label="1m",
                     step="month",
                     stepmode="backward"),
                dict(count=6,
                     label="6m",
                     step="month",
                     stepmode="backward"),
                dict(count=1,
                     label="YTD",
                     step="year",
                     stepmode="todate"),
                dict(count=1,
                     label="1y",
                     step="year",
                     stepmode="backward"),
                dict(step="all")
            ])
        ),
        rangeslider=dict(
            visible=True
        ),
        type="date"
    )
)
"""

fig.show()
