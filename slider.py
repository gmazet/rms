# -*- coding: utf-8 -*-
import plotly.graph_objects as go
import plotly.express as px

import pandas as pd

# Load data
#df = pd.read_csv( "https://raw.githubusercontent.com/plotly/datasets/master/finance-charts-apple.csv")
df = pd.read_csv( "./csv/2021-04-10_AM_R9F1B_00_SHZ.csv", sep=',')
#df.columns = [col.replace("AAPL.", "") for col in df.columns]
#df['dt']=df.Date + 'T' +  df['Time UTC']
#print df['dt']
df['Datetime']=pd.to_datetime(df['Datetime'], unit='s')
print (df.columns)

# Create figure
fig = go.Figure()
fig.add_trace( go.Scatter(x=list(df.Datetime), y=list(df.rms), mode='markers') )
#fig=px.histogram(df,x="Magnitude", nbins=5)

# Set title
fig.update_layout( title_text="RMS")

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

fig.show()
