# 3D TERRAIN MODELLING
import numpy as np
import math
import plotly.offline as go_offline
import plotly.graph_objects as go
import csv
from plotly.subplots import make_subplots
from copy import deepcopy
import pandas as pd
import plotly.express as px
import statsmodels.api as sm
import os
import chart_studio.plotly as py
import ipywidgets as widgets
import time
from plotly.offline import init_notebook_mode, plot
import matplotlib.pyplot as plt
from matplotlib import animation, rc

import mpl_toolkits.mplot3d.axes3d as p3
from IPython.display import HTML, Image
import subprocess
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import mpl_toolkits.mplot3d.axes3d as p3
from IPython.display import HTML


#--------------------------------------------------  Setup  ----------------------------------------------------------

#--General--
surveydata = 'InjectionAnalysis/Input/survey_data.csv'             #for terrain only interpolation
terrainCSV = 'InjectionAnalysis/General/TerrainInterpolation.csv'
WellsCSV = 'InjectionAnalysis/General/WellsInfo.csv'

wells = [[13000, 10000, 1000000], [13600, 9900, 1000000], [13200, 10400, 1500000]]
                                                                        #lat, lon, volume (in cubic meters)
ellipsoidratio = 200                                                      #ratio of x (or y) to z in ellipsoid eqn
n = 50                                             # number of interpolation point for x and y axis (i.e. resolution)
percentagewood = 25                                 # fraction of total volume which is woodchips within slurry

#--Evacuation--
valve = False                                       #skips the first log function in the piecewise funciton
maxrate = 50                                       #maximum rate at which water can flow through evacuating tube
firstslope = 7
secondslope = 7                                     #slope: slope of logarithmic function

#--Injection--
totaltime = 1000
timescaleinjection = 25
timescaleevacuation = 200


# READING AND PARSING THE DATA
file = open(surveydata, 'r')
lines = file.readlines()
n_line = len(lines)
x = []
y = []
z = []
for i in range(1, n_line):
    split_line = lines[i].split(",")
    xyz_t = []
    x.append(float(split_line[1].rstrip()))
    y.append(float(split_line[0].rstrip()))
    z.append(float(split_line[2].rstrip()))

#Saving well data to csv
fields = ['Well #', 'Latitude', 'Longitude', 'Volume of Slurry']
rows = [[0 for x in range(4)] for x in range(len(wells))]
curr = 0
for i in range(len(wells)):
    rows[i][0] = i
    rows[i][1] = wells[i][0]
    rows[i][2] = wells[i][1]
    rows[i][3] = wells[i][2]
    curr += 1

with open(WellsCSV, 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(fields)
    csvwriter.writerows(rows)

#--------------------------------------------------  Helper Methods  ------------------------------------------------

# DISTANCE FUNCTION
def distance(x1, y1, x2, y2):
    d = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    return d

# TO FIND COEFFICIENTS TO ELLIPSOID EQUATION
def findCoefficients(volume):
    c = ((3 * volume) / (4 * math.pi * (ellipsoidratio ** (2)))) ** (1/3)
    return (ellipsoidratio * c), c

# TO GET ALL A VALUES OF THE WELLS
def getAs():
    fin = []
    for well in wells:
        aaa, ccc = findCoefficients(well[2])
        fin.append(aaa)
    return fin

#TO FIND HEIGHT AT GIVEN COORDINATES (x,y) IN ELLIPSOID DEFINED BY a,b,c
def findHeight(aa, bb, cc, xx, yy):
    height = 2 * (cc * ((aa ** 2) * ((bb ** 2) - (yy ** 2)) - ((bb ** 2) * (xx ** 2))) ** (1/2)) / (aa * bb)
    return height

# TO FIND AND ADD HEIGHT OF ELLIPSOID AT INTERPOLATED POINTS
def updateInterpolationPoints(x_list, y_list, z_list, well, aa, bb, cc):
    yiterator = 0
    xiterator = 0
    z_new = z_list.copy()
    for xcoor in x_list:
        yiterator = 0
        for ycoor in y_list:
            if abs(distance(well[0], well[1], xcoor, ycoor)) < aa:
                z_new[yiterator][xiterator] += findHeight(aa, bb, cc, xcoor - well[0], ycoor - well[1])
            yiterator += 1
        xiterator += 1
    return z_new

# TO UPDATE INTERPOLATION MAP FOR EVERY WELL
def updateWithAllWells(x_list, y_list, z_list, limit):
    z_new = z_list.copy()
    for well in wells[0:limit]:
        aaa, ccc = findCoefficients(well[2])
        z_new = updateInterpolationPoints(x_list, y_list, z_list, well, aaa, aaa, ccc)
    return z_new

# TO UPDATE INTERPOLATION MAP FOR SPECIFIC WELL
def updateGivenWell(x_list, y_list, z_list, inttt):
    z_new = z_list.copy()
    well = wells[inttt - 1]
    aaa, ccc = findCoefficients(well[2])
    z_new = updateInterpolationPoints(x_list, y_list, z_list, well, aaa, aaa, ccc)
    return z_new

# TO UPDATE INTERPOLATION MAP FOR SPECIFIC WELL WITH C
def updateWellWithC(x_list, y_list, z_list, cc, inttt):
    z_new = z_list.copy()
    well = wells[inttt - 1]
    aaa, ccc = findCoefficients(well[2])
    z_new = updateInterpolationPoints(x_list, y_list, z_list, well, aaa, aaa, cc)
    return z_new

# TO UPDATE INTERPOLATION MAP FOR SPECIFIC WELL WITH PARAMETERS
def updateWellWithAC(x_list, y_list, z_list, aa, cc, inttt):
    z_new = z_list.copy()
    well = wells[inttt - 1]
    z_new = updateInterpolationPoints(x_list, y_list, z_list, well, aa, aa, cc)
    return z_new

#TO FIND ARRAY OF COORDINATES THAT NEED TO BE TESTED
def findCoordinates(aa, bb, cc):
    fin = []
    for count in range(0, math.ceil(aa) + 1, math.floor((math.ceil(aa) + 1) / ((math.ceil(aa) + 1) / 5))):
        for countt in range(0, math.ceil(bb) + 1, math.floor((math.ceil(bb) + 1) / ((math.ceil(bb) + 1) / 5))):
            if distance(count, countt, 0, 0) < aa:
                fin.append([count, countt])
    finn = fin.copy()
    for ele in fin:
        finn.append([(-1 * ele[0]), ele[1]])
    finnn = finn.copy()
    for ele in fin:
        finnn.append([(-1 * ele[0]), (-1 * ele[1])])
    finnnn = []
    for ele in finnn:
        if ele not in finnnn:
            finnnn.append([ele[0], ele[1], findHeight(aa, bb, cc, ele[0], ele[1])])
    return finnnn

# GENERATING COORDINATES FROM WELLS + ELLIPSOIDS
def realCoordinates(samples, well):
    for ele in samples:
        ele[0] += well[0]
        ele[1] += well[1]
    return samples

# GENERATE AN ELLIPSOID
def getOneEllipsoid(wellnum):
    well = wells[wellnum]
    aaa, ccc = findCoefficients(well[2])
    samples = findCoordinates(aaa, aaa, ccc)  # thus assuming that x and y dimensions of ellipsoid are equal -> geocores all equal
    return samples

# GENERATE ALL ELLIPSOIDS
def generateAllEllipsoids():
    fin = []
    for ele in wells:
        aaa, ccc = findCoefficients(ele[2])
        samples = findCoordinates(aaa, aaa, ccc) #thus assuming that x and y dimensions of ellipsoid are equal -> geocores all equal
        samples = realCoordinates(samples, ele)
        fin.append(samples)
    return fin

# GENERATE A, C OVER TIME OF ACTIVE INJECTION
def generatePumpingParameters(well):
    aa, cc = findCoefficients(well[2])
    fina = []
    finc = []
    for iter in range(totaltime):
        fina.append(iter * (aa / totaltime))
        finc.append(iter * (cc / totaltime))
    return fina, finc

# TO GET THE RATE AT A TIME T IN SECOND LOG FUNCTION
def getRateOne(a):
    fin = []
    end = math.exp(maxrate / a) - 1
    for inn in range(math.floor(end)):
        val = (a * math.log(inn + 1))
        fin.append(val)
    return fin, end

# TO GET THE RATE AT A TIME T IN SECOND LOG FUNCTION
def getRateTwo(a, b):
    fin = []
    start = -1 * math.exp(maxrate / a) + b + 1
    end = b
    for inn in range(math.floor(end - start)):
        val = (a * math.log((-1 * (start + inn) + 1 + b)))
        fin.append(val)
    return fin, end, start

# TO FIND VOLUME OF FIRST LOG PART
def getVolumeOne(a):
    fin, end = getRateOne(a)
    sum = 0
    for ele in fin:
        sum += abs((ele * end/len(fin)))
    return sum

# TO FIND VOLUME OF SECOND LOG PART
def getVolumeTwo(a, b):
    fin, end, start = getRateTwo(a, b)
    sum = 0
    for ele in fin:
        sum += abs((ele * (end - start)/len(fin)))
    return sum

# TO FIND NECESSARY VOLUME OF MIDDLE PART
def findMiddleLength(v1, v2, total):
    neededvol = total - (v1 + v2)
    length = neededvol / maxrate
    return length

def findMiddleLength1(v2, total):
    neededvol = total - v2
    length = neededvol / maxrate
    return length

# TO GENERATE SERIES OF RATES THROUGHOUT EVACUATION
def generateRatesOverTime(currentvolume):
    volumewater = currentvolume - (currentvolume * (percentagewood / 100))
    fin = []
    if not valve:
        tentative = getVolumeTwo(secondslope, 1000)
        mid = findMiddleLength1(tentative, volumewater)
        for _ in range(math.floor(mid)):
            fin.append(maxrate)
        ar, end, start = getRateTwo(secondslope, mid)
        for ele in ar:
            fin.append(ele)
        return fin
    else:
        tent = getVolumeOne(firstslope)
        tentative = getVolumeTwo(secondslope, 1000)
        mid = findMiddleLength1(tentative + tent, volumewater)
        ar, end = getRateOne(firstslope)
        for elem in ar:
            fin.append(elem)
        for _ in range(math.floor(mid)):
            fin.append(maxrate)
        ar, end, start = getRateTwo(secondslope, mid)
        for elem in ar:
            fin.append(elem)
        return fin

# TO CONVERT AND GENERATE A VOLUME FROM RATES
def generateVolumeOverTime(currentvolume):
    fin = generateRatesOverTime(currentvolume)
    runningtotal = 0
    finn = []
    for ele in fin:
        temp = ele
        finn.append(currentvolume - runningtotal)
        runningtotal += temp
    return finn

# TO CONVERT AND GENERATE C FROM RATES
def generateCOverTime(currentvolume):
    aa, cc = findCoefficients(currentvolume)
    fin = generateRatesOverTime(currentvolume)
    runningtotal = 0
    finn = []
    for ele in fin:
        temp = ele
        finn.append(cc * ((currentvolume - runningtotal) / currentvolume))
        runningtotal += temp
    return finn

# TO GENERATE RATES, VOLUME, AND C ARRAYS FOR EVERY WELL
def generateAllEvacuationData():
    finC = []
    finVolume = []
    finRates = []
    for ele in wells:
        finC.append(generateCOverTime(ele[2]))
        finVolume.append(generateVolumeOverTime(ele[2]))
        finRates.append(generateRatesOverTime(ele[2]))
    return finRates, finVolume, finC

# CREATING IDW FUNCTION
def idw_npoint(xz, yz, n_point, p):
    r = 10  # block radius iteration distance
    nf = 0
    while nf <= n_point:  # will stop when np reaching at least n_point
        x_block = []
        y_block = []
        z_block = []
        r += 10  # add 10 unit each iteration
        xr_min = xz - r
        xr_max = xz + r
        yr_min = yz - r
        yr_max = yz + r
        for i in range(len(x)):
            # condition to test if a point is within the block
            if ((x[i] >= xr_min and x[i] <= xr_max) and (y[i] >= yr_min and y[i] <= yr_max)):
                x_block.append(x[i])
                y_block.append(y[i])
                z_block.append(z[i])
        nf = len(x_block)  # calculate number of point in the block

    # calculate weight based on distance and p value
    w_list = []
    for j in range(len(x_block)):
        d = distance(xz, yz, x_block[j], y_block[j])
        if d > 0:
            w = 1 / (d ** p)
            w_list.append(w)
            z0 = 0
        else:
            w_list.append(0)  # if meet this condition, it means d<=0, weight is set to 0

    # check if there is 0 in weight list
    w_check = 0 in w_list
    if w_check == True:
        idx = w_list.index(0)  # find index for weight=0
        z_idw = z_block[idx]  # set the value to the current sample value
    else:
        wt = np.transpose(w_list)
        z_idw = np.dot(z_block, wt) / sum(w_list)  # idw calculation using dot product
    return z_idw

#--------------------------------------------------  Terrain Estimation Portion  -------------------------------------

# POPULATE INTERPOLATION POINTS
x_min = min(x)
x_max = max(x)
y_min = min(y)
y_max = max(y)
w = x_max - x_min  # width
h = y_max - y_min  # length
wn = w / n  # x interval
hn = h / n  # y interval

# list to store interpolation point and elevation
y_init = y_min
x_init = x_min
x_idw_list = []
y_idw_list = []
z_head = []
for i in range(n):
    xz = x_init + wn * i
    yz = y_init + hn * i
    y_idw_list.append(yz)
    x_idw_list.append(xz)
    z_idw_list = []
    for j in range(n):
        xz = x_init + wn * j
        z_idw = idw_npoint(xz, yz, 5, 1.5)  # min. point=5, p=1.5
        z_idw_list.append(z_idw)
    z_head.append(z_idw_list)

# if not os.path.exists('InjectionAnalysis/General'):
#     os.makedirs('InjectionAnalysis/General')
# # CREATING 3D TERRAIN
# fig = go.Figure()
# fig.add_trace(go.Surface(z=z_head, x=x_idw_list, y=y_idw_list))
# fig.update_layout(
#     scene=dict(aspectratio=dict(x=2, y=2, z=0.5), xaxis=dict(range=[x_min, x_max], ), yaxis=dict(range=[y_min, y_max])))
# fig.update_layout(title_text="Terrain Interpolation")
# go_offline.plot(fig, filename='InjectionAnalysis/General/TerrainInterpolation.html',
#                 validate=True, auto_open=False)
#
# #Saving to a csv file for future processing
# fields = ['Latitude', 'Longitude', 'Elevation']
# rows = [[0 for x in range(3)] for x in range(n*n)]
# curr = 0
# for i in range(n):
#     for j in range(n):
#         rows[curr][0] = y_idw_list[j]
#         rows[curr][1] = x_idw_list[n-i-1]
#         rows[curr][2] = z_head[j][n-i-1]
#         curr += 1
#
# with open(terrainCSV, 'w') as csvfile:
#     csvwriter = csv.writer(csvfile)
#     csvwriter.writerow(fields)
#     csvwriter.writerows(rows)
#
# #--------------------------------------------------  Graphics Generation Portion  -------------------------------------
#
# #---Wells---
#
#
# for ele in range(len(wells)):
#     z_head1 = deepcopy(z_head)
#     z_head1 = updateGivenWell(x_idw_list, y_idw_list, z_head1, ele + 1)
#     if not os.path.exists('InjectionAnalysis/IndividualWells/Well #' + str(ele)):
#         os.makedirs('InjectionAnalysis/IndividualWells/Well #' + str(ele))
#
#     fig = go.Figure()
#     fig.add_trace(go.Surface(z=z_head1, x=x_idw_list, y=y_idw_list))
#     fig.update_layout(
#         scene=dict(aspectratio=dict(x=2, y=2, z=0.5), xaxis=dict(range=[x_min, x_max], ),
#                    yaxis=dict(range=[y_min, y_max])))
#     fig.update_layout(title_text="Terrain Estimation After Initial Injection in Well #" + str(ele))
#     go_offline.plot(fig, filename='InjectionAnalysis/IndividualWells/Well #' + str(ele) +
#                                   '/TerrainEstimation.html', validate=True, auto_open=False)
#
# # #----Multiple Injections----
#
# if not os.path.exists('InjectionAnalysis/General'):
#     os.makedirs('InjectionAnalysis/General')
#
# fig = make_subplots(
#     rows=2, cols=2,
#     specs=[[{'type': 'surface'}, {'type': 'surface'}], [{'type': 'surface'}, {'type': 'surface'}]]
#     , subplot_titles=("Original Terrain", "Post-Injection (Single Well)", "Post-Injection (Double Well)",
#                       "Post-Injection (Triple Well)"))
#
# # adding surfaces to subplots.
# fig.add_trace(
#     go.Surface(x=x_idw_list, y=y_idw_list, z=z_head, colorscale='earth', showscale=False),
#     row=1, col=1)
#
# z_first = deepcopy(z_head)
# z_first = updateGivenWell(x_idw_list, y_idw_list, z_first, 1)
#
# fig.add_trace(
#     go.Surface(x=x_idw_list, y=y_idw_list, z=z_first, colorscale='RdBu', showscale=False),
#     row=1, col=2)
#
# z_first = updateGivenWell(x_idw_list, y_idw_list, z_first, 2)
#
# fig.add_trace(
#     go.Surface(x=x_idw_list, y=y_idw_list, z=z_first, colorscale='Viridis', showscale=False),
#     row=2, col=1)
#
# z_first = updateGivenWell(x_idw_list, y_idw_list, z_first, 3)
#
# fig.add_trace(
#     go.Surface(x=x_idw_list, y=y_idw_list, z=z_first, showscale=False),
#     row=2, col=2)
#
# fig.update_layout(
#     scene=dict(xaxis=dict(range=[x_min, x_max], ), yaxis=dict(range=[y_min, y_max])))
# fig.update_scenes(aspectratio=dict(x=2, y=2, z=0.5))
#
# fig.update_layout(
#     title_text='Slurry Injection Simulation',
#
# )
#
# go_offline.plot(fig, filename='InjectionAnalysis/General/First3Wells.html', validate=True,
#                 auto_open=False)
#
# # #----Injection Over Time----
#
# allrates, allvolume, allC = generateAllEvacuationData()
#
# for ele in range(len(allrates)):
#     if not os.path.exists('InjectionAnalysis/IndividualWells/Well #' + str(ele)):
#         os.makedirs('InjectionAnalysis/IndividualWells/Well #' + str(ele))
#     ttt = allrates[ele]
#     fig = px.line(x=range(len(ttt)), y=ttt)
#     fig.update_layout(
#         xaxis_title="Time",
#         yaxis_title="Rate",
#     )
#     fig.update_layout(title_text="Rate of Water Evacuation in Well #" + str(ele))
#     go_offline.plot(fig, filename='InjectionAnalysis/IndividualWells/Well #' + str(ele) +
#                                   '/RateData.html',
#                     validate=True, auto_open=False)
#
#     ttt = allvolume[ele]
#     fig = px.line(x=range(len(ttt)), y=ttt)
#     fig.update_layout(
#         xaxis_title="Time",
#         yaxis_title="Rate",
#     )
#     fig.update_layout(title_text="Water Volume in Well #" + str(ele))
#     go_offline.plot(fig, filename='InjectionAnalysis/IndividualWells/Well #' + str(ele) +
#                                   '/VolumeData.html',
#                     validate=True, auto_open=False)
#
#     ttt = allC[ele]
#     fig = px.line(x=range(len(ttt)), y=ttt)
#     fig.update_layout(
#         xaxis_title="Time",
#         yaxis_title="Rate",
#     )
#     fig.update_layout(title_text="Ellipsoidal Height of Well #" + str(ele))
#     go_offline.plot(fig, filename='InjectionAnalysis/IndividualWells/Well #' + str(ele) +
#                                   '/CData.html',
#                     validate=True, auto_open=False)
#
# # -- Map of Measured -> Interpolated Points in Terrain --
#
# ymapvisual = []
# xmapvisual = []
# ii = 0
# for elemm in y_idw_list:
#     for iii in range(len(x_idw_list)):
#         ymapvisual.append(y_idw_list[iii])
#     for iii in range(len(y_idw_list)):
#         xmapvisual.append(x_idw_list[ii])
#     ii += 1
#
# fig = go.Figure()
#
# fig.add_trace(go.Scatter(x=xmapvisual, y=ymapvisual,
#                     mode='markers',
#                     name='Interpolated Points'))
#
# fig.add_trace(go.Scatter(x=x, y=y,
#                     mode='markers',
#                     name='Measured Points'))
#
# xmapvisuall = []
# ymapvisuall = []
# rrr = []
# for eleee in getAs():
#     rrr.append(eleee / 2)
# for well in wells:
#     xmapvisuall.append(well[0])
#     ymapvisuall.append(well[1])
# fig.add_trace(go.Scatter(
#     x=xmapvisuall,
#     y=ymapvisuall,
#     mode='markers',
#     name='Affected Area',
#     marker=dict(size=rrr)
# ))
#
# fig.update_layout(
#         xaxis_title="X",
#         yaxis_title="Y",
#     )
# fig.update_layout(title_text="Important Markers")
#
# go_offline.plot(fig, filename='InjectionAnalysis/General/ScatterOfTerrain.html',
#                     validate=True, auto_open=False)
#
# # # -- Animation of Wells Over Time--
#
# init_notebook_mode()
# fig = go.Figure()
#
# fig_dict = {
#     "data": [],
#     "layout": {},
#     "frames": []
# }
#
# sliders_dict = {
#     "active": 0,
#     "yanchor": "top",
#     "xanchor": "left",
#     "currentvalue": {
#         "font": {"size": 20},
#         "prefix": "Time:",
#         "visible": True,
#         "xanchor": "right"
#     },
#     "transition": {"duration": 300, "easing": "cubic-in-out"},
#     "pad": {"b": 10, "t": 50},
#     "len": 0.9,
#     "x": 0.1,
#     "y": 0,
#     "steps": []
# }
#
# Covertime = generateCOverTime(wells[0][2])
# nco = []
# Aoverinjection, Coverinjection = generatePumpingParameters(wells[0])
# aaa = []
# print(Aoverinjection)
# for eleme in range(len(Aoverinjection)):
#     if (eleme % timescaleinjection == 0):
#         aaa.append(Aoverinjection[eleme])
# for eleme in range(len(Aoverinjection)):
#     if (eleme % timescaleinjection == 0):
#         nco.append(Coverinjection[eleme])
# for eleme in range(len(Covertime)):
#     if (eleme % timescaleevacuation == 0):
#         nco.append(Covertime[eleme])
#
# z_tnm = np.zeros((math.floor(len(nco)),n,n), dtype=float)
#
# for t in range(math.floor(len(nco))):
#     z_nn = deepcopy(z_head)
#     if (t < len(aaa)):
#         z_nn = updateWellWithAC(x_idw_list, y_idw_list, z_nn, aaa[t], nco[t], 1)
#     else:
#         z__nn = updateWellWithC(x_idw_list, y_idw_list, z_nn, nco[t], 1)
#     tsunami = np.array(z_nn)
#     z_tnm[t] = tsunami
#     slider_step = {"args": [
#         [tsunami],
#         {"frame": {"duration": 300, "redraw": False},
#          "mode": "immediate",
#          "transition": {"duration": 300}}
#     ],
#         "label": t,
#         "method": "animate"}
#     sliders_dict["steps"].append(slider_step)
#
# # Create FigureWidget and add surface trace
# #surface = fig.add_surface(z=z_bath,name='bath')
#
# # ----- change part ------
# bath = go.Surface(z=z_head)
# fig = go.Figure(data=[bath,bath])
# fig.add_trace(
#     go.Surface(x=x_idw_list, y=y_idw_list, z=z_head, colorscale='RdBu', showscale=False))
#
#
# # Set axis ranges to fixed values to keep them from resetting during animation
# fig.update_layout(
#     scene=dict(aspectratio=dict(x=2, y=2, z=0.5), xaxis=dict(range=[x_min, x_max], ), yaxis=dict(range=[y_min, y_max])))
#
#
# frames = []
# for i in range(math.floor(len(nco))):
#     frames.append(go.Frame(data=[{'type': 'surface', 'z': z_tnm[i], 'name':'tsunami'}]))
#
# fig.frames = frames
#
# fig.layout.updatemenus = [
#     {
#         "buttons": [
#             {
#                 "args": [None, {"frame": {"duration": 300, "redraw": False},
#                                 "fromcurrent": True, "transition": {"duration": 300,
#                                                                     "easing": "quadratic-in-out"}}],
#                 "label": "Play",
#                 "method": "animate"
#             },
#             {
#                 "args": [[None], {"frame": {"duration": 0, "redraw": False},
#                                   "mode": "immediate",
#                                   "transition": {"duration": 0}}],
#                 "label": "Pause",
#                 "method": "animate"
#             }
#         ],
#         "direction": "left",
#         "pad": {"r": 10, "t": 87},
#         "showactive": False,
#         "type": "buttons",
#         "x": 0.1,
#         "xanchor": "right",
#         "y": 0,
#         "yanchor": "top"
#     }
# ]
# fig.layout.sliders = [sliders_dict]
#
# plot(fig)
#
# go_offline.plot(fig, filename='InjectionAnalysis/General/animation.html',
#                     validate=True, auto_open=False)

# ------------------------------------------- DASHBOARD MOCK UP (unfinished) -----------------------------------------

#
# Covertime = generateCOverTime(wells[0][2])
# nco = []
# Aoverinjection, Coverinjection = generatePumpingParameters(wells[0])
# aaa = []
#
# for eleme in range(len(Aoverinjection)):
#     if (eleme % timescaleinjection == 0):
#         aaa.append(Aoverinjection[eleme])
# for eleme in range(len(Aoverinjection)):
#     if (eleme % timescaleinjection == 0):
#         nco.append(Coverinjection[eleme])
# for eleme in range(len(Covertime)):
#     if (eleme % timescaleevacuation == 0):
#         nco.append(Covertime[eleme])
# z_tnm = np.zeros((math.floor(len(nco)),n,n), dtype=float)
# tsunami = []
# for t in range(math.floor(len(nco))):
#     z_nn = deepcopy(z_head)
#     if (t < len(aaa)):
#         z_nn = updateWellWithAC(x_idw_list, y_idw_list, z_nn, aaa[t], nco[t], 1)
#     else:
#         z_nn = updateWellWithC(x_idw_list, y_idw_list, z_nn, nco[t], 1)
#     tsunami = np.array(z_nn)
#     z_tnm[t] = tsunami
#
# new_x_idw_list = []
# for length in range(len(x_idw_list)):
#     new_x_idw_list.append(x_idw_list)
# new_y_idw_list = []
# for length in range(len(y_idw_list)):
#     new_y_idw_list.append(y_idw_list)
# #
# # N = 500 # Meshsize
# # fps = 10 # frame per sec
# # frn = math.floor(len(nco)) # frame number of the animation
# #
# # fig = plt.figure()
# # ax = Axes3D(fig)
# #
# # def init():
# #     # Plot the surface.
# #     ax.plot_surface(new_x_idw_list, new_y_idw_list, z_tnm[0], cmap="magma",
# #                     linewidth=0, antialiased=False)
# #     return fig
# #
# # def update_plot(frame_number, zarray, plot):
# #     plot[0].remove()
# #     plot[0] = ax.plot_surface(new_x_idw_list, new_y_idw_list, z_tnm[frame_number], cmap="magma")
# #     # ax.view_init(elev=10, azim=frame_number * 4)
# #     return fig
# #
# #
# # plot = [ax.plot_surface(new_x_idw_list, new_y_idw_list, z_tnm[0], color='0.75', rstride=1, cstride=1)]
# # ax.set_zlim(598,614)
# # ax.set_aspect('auto')
# # ani = animation.FuncAnimation(fig, update_plot, init_func=init, fargs=(z_tnm, plot), interval=1000/fps)
# #
# # plt.show()
#
#
# X = np.arange(-5, 5, 0.25)
# Y = np.arange(-5, 5, 0.25)
# X, Y = np.meshgrid(X, Y)
# R = np.sqrt(X**2 + Y**2)
# Z = np.sin(R)
#
# # Create a figure and a 3D Axes
# fig = plt.figure()
# ax = Axes3D(fig)
#
# def init():
#     # Plot the surface.
#     ax.plot_surface(np.array(new_x_idw_list), np.array(new_y_idw_list), np.array(z_tnm[0]), cmap="magma",
#                     linewidth=0, antialiased=False)
#     return fig,
#
# def animate(i):
#     # elevation angle : -180 deg to 180 deg
#     ax.view_init(elev=(i-45)*4, azim=10)
#     return fig,
#
# # Animate
# ani = animation.FuncAnimation(fig, animate, init_func=init,
#                               frames=90, interval=50, blit=True)








