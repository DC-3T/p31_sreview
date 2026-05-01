"""
Script for generating a summary 31P-MRS figure illustrating static and
dynamic data acquired in skeletal muscle. Author: Donnie Cameron
"""

import math

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from skimage import io


# Begin by initialising a subplot figure comprising four figures
fig1 = make_subplots(rows=2, cols=2, column_widths=[0.4, 0.6],
                     specs=[[{}, {"type": "scatter3d"}],
                            [{}, {"secondary_y": True}]],
                     horizontal_spacing=0.075,  # Reduce horizontal space
                     vertical_spacing=0.025  # Reduce vertical space
                    )

# Add figure 1a - static axial image of thigh
img = io.imread('imgs/fig1a.jpg')
fig1.add_trace(go.Image(z=img), row=1, col=1)
fig1.update_xaxes(visible=False, row=1, col=1)
fig1.update_yaxes(visible=False, row=1, col=1)
fig1.update_layout(showlegend=False, margin=dict(r=5, l=5, b=5, t=5))


# Add figure 1b - 3D line plot of spectra
df2 = pd.read_csv('data/p31_ex_fig1.csv')  # Load spectra

time = 0
x_v = []
for i in range(np.shape(df2)[1]):
    x_v_temp = [1 * time] * np.shape(df2)[0]
    x_v_temp.extend([None])
    x_v.extend(x_v_temp)
    time += 9

y_v = []
for i in range(np.shape(df2)[1]):
    y_v_temp = list(range(np.shape(df2)[0]))
    y_v_temp.extend([None])
    y_v.extend(y_v_temp)

df2.loc[-1] = [None] * np.shape(df2)[1]
sigs = [c for c in df2 if c.startswith('Signal')]

df_fig1b = pd.melt(df2, value_vars=sigs, value_name='signal')

df_fig1b.insert(1, 'spec_no', x_v)
df_fig1b.insert(2, 'point', y_v)
df_fig1b_short = df_fig1b[0 : 2048 * 44]

# Convert x-axis from pts to ppm
df_fig1b_short["point"] = df_fig1b_short["point"].sub(1024.)
df_fig1b_short["point"] = df_fig1b_short["point"].mul(2994./2048.)
df_fig1b_short["point"] = df_fig1b_short["point"].div(51.7)

# Plot 3D scatter plot of spectra
fig1.add_trace(
    go.Scatter3d(x=df_fig1b_short["point"], y=df_fig1b_short["spec_no"],
                 z=df_fig1b_short["signal"], mode='lines',
                 line=dict(color=df_fig1b_short["signal"].fillna(0), width=7,
                           colorscale='RdBu_r', cmin=3E7, cmax=6E8)),
    row=1, col=2)

fig1.update_layout(
    scene=dict(
        bgcolor='rgba(0,0,0,0)',
        xaxis=dict(title="Frequency, ppm", range=[-18.5, 8.5],
                   backgroundcolor="rgba(0, 0, 0,0)", showgrid=False,
                   zerolinecolor="white",),
        yaxis=dict(title="Time, s", backgroundcolor="rgba(0, 0, 0,0)",
                   showgrid=False, zerolinecolor="white",),
        zaxis=dict(visible=False, backgroundcolor="rgba(0, 0, 0,0)",
                   showgrid=False, zerolinecolor="white",),
    aspectratio=dict(
        x=0.7, y=1, z=0.25),
camera=dict(
            up=dict(
                x=math.cos(math.pi/2),
                y=0,
                z=math.sin(math.pi/2)
            ),
            center=dict(
                x=-0.035,
                y=-0.035,
                z=-0.3
            ),
            eye=dict(
                x=0.65,
                y=-0.65,
                z=0.52,
            ))),
    width=1000,
    height=700,
    margin=dict(r=0, l=0, b=5, t=0),
    dragmode='orbit')


# Add figure 1c - plot of single spectrum
df3 = pd.read_csv('data/p31_rest_fig1.csv')  # Load spectrum data

# Convert x-axis from pts to ppm
df3["point"] = list(range(2048))
df3["point"] = df3["point"].sub(1024.)
df3["point"] = df3["point"].mul(2994./2048.)
df3["point"] = df3["point"].div(51.7)

# Plot spectrum
fig1.add_trace(
    go.Scatter(x=df3["point"], y=df3["fft(real)"], name="Spectrum",
               line=dict(color="#33334c")), secondary_y=False, row=2, col=1)
fig1.update_layout(showlegend=False, margin=dict(r=2, l=5, b=5, t=5))

# Add axis titles
fig1.update_xaxes(title="Frequency, ppm", range=[15, -22], row=2, col=1)
fig1.update_yaxes(title="Signal intensity, a.u.", row=2, col=1)

# Label PCr
fig1.add_trace(
    go.Scatter(x=[-2], y=[1.2E9], mode="text", name="Text", text="PCr",
               textfont=dict(size=14, color="#33334c"),
               textposition="top center"), row=2, col=1)

# Label PDE
fig1.add_trace(
    go.Scatter(x=[3.0], y=[70E6], mode="text", name="Text", text="PDE",
               textfont=dict(size=14, color="#33334c"),
               textposition="top center"), row=2, col=1)

# Label Pi
fig1.add_trace(
    go.Scatter(x=[5.8], y=[114E6], mode="text", name="Text", text="Pi",
               textfont=dict(size=14, color="#33334c"),
               textposition="top center"), row=2, col=1)

# Label gamma ATP
fig1.add_trace(
    go.Scatter(x=[-3], y=[170E6], mode="text", name="Text",
               text="\u03b3" + "ATP", textfont=dict(size=14, color="#33334c"),
               textposition="top center"), row=2, col=1)

# Label alpha ATP
fig1.add_trace(
    go.Scatter(x=[-8.2], y=[160E6], mode="text", name="Text",
               text="\u03b1" + "ATP", textfont=dict(size=14, color="#33334c"),
               textposition="top center"), row=2, col=1)

# Label beta ATP
fig1.add_trace(
    go.Scatter(x=[-15.8], y=[140E6], mode="text", name="Text",
               text="\u03B2" + "ATP", textfont=dict(size=14, color="#33334c"),
               textposition="top center"), row=2, col=1)


# Add figure 1d - time course of 31P metabolites / pH
df = pd.read_csv('data/p31_time_course.csv')  # Load example 31P time course

# Add traces
# Box over exercise period
fig1.add_trace(
    go.Scatter(x=[57, 57, 117, 117], y=[0, 3.25E7, 3.25E7, 0], fill="toself",
               fillcolor="rgba(210, 210, 210, 0.4)",
               mode="text"), row=2, col=2)

# Add phosphocreatine timecourse
fig1.add_trace(
     go.Scatter(x=df['time_s'] * 3, y=df["pcr_amp"], name="yaxis data",
                line=dict(color="#33334c")), row=2, col=2, secondary_y=False)

# Add inorganic phosphate timecourse
fig1.add_trace(
    go.Scatter(x=df['time_s'] * 3, y=df["pi_amp"], name="yaxis data",
               line=dict(color="#dd3d2d")), row=2, col=2,secondary_y=False)

# Add pH timecourse
fig1.add_trace(
    go.Scatter(x=df['time_s'] * 3, y=df["ph"], name="yaxis2 data",
               line=dict(color="#6ea6cd")), row=2, col=2, secondary_y=True)

# Add PCr label
fig1.add_trace(
    go.Scatter(x=[370], y=[29E6], mode="text", name="Text", text="PCr",
               textfont=dict(size=18, color="#33334c"),
               textposition="top center"), row=2, col=2)

# Add Pi label
fig1.add_trace(
    go.Scatter(x=[370], y=[3.5E6], mode="text", name="Text", text="Pi",
               textfont=dict(size=18, color="#dd3d2d"),
               textposition="top center"), row=2, col=2)

# Add pH label
fig1.add_trace(
    go.Scatter(x=[370], y=[18.5E6], mode="text", name="Text", text="pH",
               textfont=dict(size=18, color="#6ea6cd"),
               textposition="top center"), row=2, col=2)

# Use Scattergl to plot arrow on top of other plots
fig1.add_trace(
    go.Scattergl(x=[57, 57], y=[29.35E6, 16.2E6], mode="markers+lines",
                 marker_symbol=["arrow-bar-up", "arrow-bar-down"],
                 marker_size=[12, 12], marker_line_color="#33334c",
                 marker_line_width=2, marker_color="#33334c",
                 line_color="#33334c"), row=2, col=2)

# Add PCr drop label
fig1.add_trace(
    go.Scatter(x=[28.5], y=[22E6], mode="text", name="Text", text="PCr drop",
               textfont=dict(size=14, color="#33334c"),
               textposition="top center"), row=2, col=2)

# Add tauPCr label
fig1.add_trace(
    go.Scatter(x=[144], y=[25E6], mode="text", name="Text",
               text="\u03c4" + "<sub>PCr</sub>",
               textfont=dict(size=16, color="#33334c"),
               textposition="top center"), row=2, col=2)

# Add end-exercise line
fig1.add_trace(
    go.Scatter(x=[117, 117], y=[0, 3.25E7], mode="lines", line_color="black",
               line_dash="dashdot"),row=2, col=2)

# Add end-exercise label
fig1.add_annotation(
    dict(x=0.576, y=0.035, xref="paper", yref="paper", text="end exercise",
    font=dict(size=14, color="black"), textangle=90))

# Set x-axis title
fig1.update_xaxes(title_text="Time, seconds", row=2, col=2)

# Set y-axes titles
fig1.update_yaxes(title_text="Signal intensity, a.u.", range=[0, 3.25E7],
                  secondary_y=False, row=2, col=2)
fig1.update_yaxes(title_text="pH", title_font_color="#6ea6cd",
                  tickfont_color="#6ea6cd",
                 range=[6.75, 7.25], secondary_y=True, row=2, col=2)
fig1.update_layout(showlegend=False, margin=dict(r=5, l=5, b=5, t=5))


# Adding annotations to set labels for all subplots
annotations = [{'font': {'size': 30, 'family': 'Arial', 'color': 'white'},
                'showarrow': False,
                'text': '<b>a</b>',
                'x': 0.025,
                'xanchor': 'center',
                'xref': 'paper',
                'y': 0.94,
                'yanchor': 'bottom',
                'yref': 'paper'},
               {'font': {'size': 30, 'family': 'Arial', 'color': 'black'},
                'showarrow': False,
                'text': '<b>b</b>',
                'x': 0.445,
                'xanchor': 'center',
                'xref': 'paper',
                'y': 0.94,
                'yanchor': 'bottom',
                'yref': 'paper'},
               {'font': {'size': 30, 'family': 'Arial', 'color': 'black'},
                'showarrow': False,
                'text': '<b>c</b>',
                'x': 0.025,
                'xanchor': 'center',
                'xref': 'paper',
                'y': 0.43,
                'yanchor': 'bottom',
                'yref': 'paper'},
               {'font': {'size': 30, 'family': 'Arial',  'color': 'black'},
                'showarrow': False,
                'text': '<b>d</b>',
                'x': 0.445,
                'xanchor': 'center',
                'xref': 'paper',
                'y': 0.43,
                'yanchor': 'bottom',
                'yref': 'paper'}
               ]
fig1.update_layout(annotations=annotations)

fig1.show()

# Save image to png and interacitve html
pio.write_image(fig1, "fig1_complete.png", format="png",
                width=1000, height=700, scale=3)
fig1.write_html("fig1.html")

