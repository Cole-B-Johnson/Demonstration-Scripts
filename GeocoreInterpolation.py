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
import pickle

#--------------------------------------------------  Setup  ----------------------------------------------------------

#--General--
terraindata = 'InjectionAnalysis/Input/survey_data.csv'  #surface terrain data that the
                                                                                    # geocores
                                                                                    # were extracted from

surveydata = ['InjectionAnalysis/Input/Geocore_Survey1.csv',
              'InjectionAnalysis/Input/Geocore_Survey2.csv',
              'InjectionAnalysis/Input/Geocore_Survey3.csv',
              'InjectionAnalysis/Input/Geocore_Survey4.csv',
              'InjectionAnalysis/Input/Geocore_Survey5.csv',
              'InjectionAnalysis/Input/Geocore_Survey6.csv']      #one csv file per geocore

n = 5                                             # number of interpolation point for x and y axis (i.e. resolution)
geocoredata = []
distinguishthresholdvoid = 3
distinguishthresholdpart = 3

#--------------------------------------------------  Classes  ------------------------------------------------

# GEOCORE CLASS FOR STORAGE
class Geocore:
    voidfraction = {}
    chemistry = {}
    particlesize = {}
    latitude = 0
    longitude = 0
    breakdown = []
    continuousparticlesize = {}
    continuousvoid = []
    averageprofile = {}     # for every different makeup of material throughout geocore, average of part and void

    def __init__(self, c, v, p, lat, lon):
        self.voidfraction = v
        self.chemistry = c
        self.particlesize = p
        self.latitude = lat
        self.longitude = lon
        self.breakdown = generateBreakdown(self)
        self.continuousparticlesize = continuousparticlesize(self.particlesize)
        self.continuousvoid = continuousvoid(self.voidfraction)
        self.averageprofile = generateaverageprofile(self)

    def __str__(self):
        return ("Void Fraction: " + str(self.voidfraction) + "\nChemistry: " + str(self.chemistry) + "\nParticle Size: "
                + str(self.particlesize) + "\nLatitude: " + str(self.latitude) + "\nLongitude: " + str(self.longitude)
                + "\nBreakdown: " + str(self.breakdown))

#--------------------------------------------------  Helper Methods  ------------------------------------------------

# TO GENERATE GEOCORE OBJECT FOR USE
def generateGeocore(filenum):
    file = open(surveydata[filenum])
    lines = file.readlines()
    n_line = len(lines)
    chem = {}
    void = {}
    particlesize = {}
    firstline = lines[1].split(",")
    lat = float(firstline[4].rstrip())
    lon = float(firstline[5].rstrip())
    chem[int(firstline[0])] = firstline[1].rstrip().split(";")
    void[int(firstline[0])] = float(firstline[2].rstrip())
    temppart = firstline[3].rstrip()
    ttemppart = {}
    for ele in temppart.split(";"):
        ttemppart[ele.split(":")[0].lstrip().rstrip()] = float(ele.split(":")[1].lstrip().rstrip())
    particlesize[int(firstline[0])] = ttemppart
    for i in range(2, n_line):
        split_line = lines[i].split(",")
        depth = split_line[0]
        tempchem = split_line[1].rstrip()
        tempvoid = split_line[2].rstrip()
        temppart = split_line[3].rstrip()

        chem[int(depth)] = tempchem.split(";")
        void[int(depth)] = float(tempvoid)
        ttemppart = {}
        for ele in temppart.split(";"):
            ttemppart[ele.split(":")[0].lstrip().rstrip()] = float(ele.split(":")[1].lstrip().rstrip())
        particlesize[int(depth)] = ttemppart
    gc = Geocore(chem, void, particlesize, lat, lon)
    return gc

# HELPER TO GENERATE CONTINUOUS DICTIONARY FOR VOID
def continuousvoid(void):
    contvoid = np.zeros(math.ceil(max(void.keys())) + 1)
    for ele in void.keys():
        contvoid[ele] = void[ele]
    prevele = 0
    for ele in range(len(contvoid)):
        if (contvoid[ele] == 0):
            continue
        else:
            for element in range(ele):
                contvoid[element] = contvoid[ele]
            break
    for ele in range(len(contvoid)):
        if (contvoid[ele] == 0):
            continue
        else:
            for element in range(prevele, ele):
                contvoid[element] = contvoid[prevele] + ((contvoid[ele] - contvoid[prevele]) *
                                                         abs((element - (prevele)) / (ele - (prevele))))
            prevele = ele
    return contvoid

# HELPER TO GENERATE CONTINUOUS DICTIONARY FOR PARTICLE SIZE
def continuousparticlesize(partsize):
    contpart = []
    for ele in range(math.ceil(max(partsize.keys())) + 1):
        if (ele in partsize.keys()):
            contpart.append(partsize[ele])
        else:
            contpart.append({})
    prevele = 0
    for ele in range(len(contpart)):
        if (contpart[ele] == {}):
            prevele += 1
            continue
        else:
            for element in range(ele):
                contpart[element] = contpart[ele]
            prevele += 1
            break
    for ele in range(prevele, len(contpart)):
        if (contpart[ele] == {}):
            continue
        else:
            for element in range(prevele, ele):
                for key in contpart[ele].keys():
                    if (key in contpart[prevele]):
                        contpart[element][key] = contpart[prevele][key] + ((contpart[ele][key] -
                                                                            contpart[prevele][key]) *
                                                            abs((element - (prevele)) / (ele - (prevele))))
                    else:
                        if (element > prevele):
                            contpart[element][key] = ((contpart[ele][key]) *
                                                                           abs((element - (prevele)) / (
                                                                                       ele - (prevele))))
                for key in contpart[prevele].keys():
                    if key not in contpart[ele]:
                        contpart[element][key] = contpart[prevele][key] - ((contpart[prevele][key]) *
                                                  abs((element - (prevele)) / (
                                                          ele - (prevele))))
            prevele = ele
    for ele in range(len(contpart)):
        if contpart[ele] == {}:
            contpart[ele] = contpart[ele - 1]
    return contpart                 #if compound doesnt exist in previous or future sample, then it is linearly interpd

# TO MESS WITH HOW TO DISTINGUISH LAYERS
def distinguishlayers(void, part):
    return ((void > distinguishthresholdvoid) or (part > distinguishthresholdpart))

# TO GENERATE BREAKDOWNS OF A GEOCORE
def generateBreakdown(geocore):
    layercounter = 0
    finalgeocore = []
    prevpart = {}
    prevvoid = 0
    contvoid = continuousvoid(geocore.voidfraction)
    contfrac = continuousparticlesize(geocore.particlesize)
    totaldepth = len(contvoid)
    for ele in range(totaldepth):
        partsizedistinguisher = 0
        for composite in contfrac[ele]:
            if (composite in prevpart.keys()):
                partsizedistinguisher += abs(contfrac[ele][composite] - prevpart[composite])
            else:
                partsizedistinguisher += abs(contfrac[ele][composite])
        prevpart = contfrac[ele]
        voiddistinguisher = abs(contvoid[ele] - prevvoid)
        prevvoid = contvoid[ele]
        finalgeocore.append(layercounter)
        if (distinguishlayers(voiddistinguisher, partsizedistinguisher)):
            layercounter += 1
    return finalgeocore

# TO GENERATE AVERAGE PROFILE FOR EVERY DIFFERENT BREAKDOWN IN A GEOCORE (dictionary of tuples)
def generateaverageprofile(geocore):
    breakdown = geocore.breakdown
    part = geocore.continuousparticlesize
    void = geocore.continuousvoid
    listofdepths = [0]
    curr = breakdown[0]
    averageprofile = {}
    for ele in range(len(breakdown)):
        if breakdown[ele] != curr:
            listofdepths.append(ele)
            curr = breakdown[ele]

    if (listofdepths[0] != 0):
        partcounter = {}
        voidcounter = 0
        for ele in range(0, listofdepths[0]):
            for key in part[ele].keys():
                if key in partcounter.keys():
                    partcounter[key] += part[ele][key]
                else:
                    partcounter[key] = part[ele][key]
            voidcounter += void[ele]
        voidcounter /= listofdepths[0]
        for ele in partcounter.keys():
            partcounter[ele] /= listofdepths[0]
        averageprofile[0] = (voidcounter, partcounter)

    for elem in range(1, len(listofdepths)):
        partcounter = {}
        voidcounter = 0
        for ele in range(listofdepths[elem - 1], listofdepths[elem]):
            for key in part[ele].keys():
                if key in partcounter.keys():
                    partcounter[key] += part[ele][key]
                else:
                    partcounter[key] = part[ele][key]
            voidcounter += void[ele]
        voidcounter /= (listofdepths[elem] - listofdepths[elem - 1])
        for ele in partcounter.keys():
            partcounter[ele] /= (listofdepths[elem] - listofdepths[elem - 1])
        averageprofile[listofdepths[elem]] = (voidcounter, partcounter)

    partcounter = {}
    voidcounter = 0
    for ele in range(listofdepths[len(listofdepths) - 1], len(part)):
        for key in part[ele].keys():
            if key in partcounter.keys():
                partcounter[key] += part[ele][key]
            else:
                partcounter[key] = part[ele][key]
        voidcounter += void[ele]
    voidcounter /= (len(part) - listofdepths[len(listofdepths) - 1])
    for ele in partcounter.keys():
        partcounter[ele] /= (len(part) - listofdepths[len(listofdepths) - 1])
    averageprofile[listofdepths[len(listofdepths) - 1]] = (voidcounter, partcounter)
    return averageprofile

# GIVEN A PROFILE, GENERATE ALL DEPTHS IN A GEOCORE THAT MATCH MATERIAL PROFILE
def generateAllDepthsMatching(profile, geocore):
    fin = []
    void, part = profile
    avgmap = geocore.averageprofile
    for depth in avgmap.keys():
        summ = 0
        for ele in avgmap[depth][1].keys():
            if ele in part:
                summ += abs(part[ele] - avgmap[depth][1][ele])
            else:
                summ += abs(avgmap[depth][1][ele])
        for ele in part.keys():
            if ele not in avgmap[depth][1].keys():
                summ += abs(part[ele])
        if not distinguishlayers(abs(void - avgmap[depth][0]), abs(summ)):
            fin.append(depth)
    return fin

# GENERATE ALL GEOCORES
def generateAllGeocores():
    for ele in range(len(surveydata)):
        geocoredata.append(generateGeocore(ele))
    return geocoredata

# GENERATE DICTIONARY IN FORM OF {PROFILE: [[LAT, LON, DEPTH]...], ...} FOR ALL PROFS AND ALL MATCHING LOCATIONS
def generateProfileDictionary(allgeocores = geocoredata):
    profiles = {}
    for geocore in allgeocores:
        for profs in geocore.averageprofile.values():
            xxx = [(profs[0])]
            yyy = []
            for key in profs[1]:
                yyy.append((key, profs[1][key]))
            xxx.append(tuple(yyy))
            profiles[tuple(xxx)] = []
    for geocore in allgeocores:                         #taking each geocore,
        for profs in geocore.averageprofile.values():   #and for each material profile contained,
            for key in geocore.averageprofile.keys():
                if geocore.averageprofile[key] == profs: #I'm finding the depth that the material profile is present
                                                            #I then need to check if it is approx.
                                                            # already present in the dictionary
                    skip = False
                    for materialalreadypresent in profiles.keys():
                        partdiff = 0
                        voiddiff = abs(geocore.averageprofile[key][0] - materialalreadypresent[0])
                        for tuplee in materialalreadypresent[1]:
                            if tuplee[0] in geocore.averageprofile[key][1]:
                                partdiff += abs(tuplee[1] - geocore.averageprofile[key][1][tuplee[0]])
                            else:
                                partdiff += tuplee[1]

                        if not distinguishlayers(voiddiff, partdiff):
                            profiles[materialalreadypresent].append([geocore.latitude, geocore.longitude, key])
                            skip = True

                    if not skip:
                        xxx = [(profs[0])]                      #If it is, then I'm appending it to the existing material
                                                                # profile
                        yyy = []
                        for k in profs[1]:                      #If not, then I am adding the material to the dictionary
                            yyy.append((k, profs[1][k]))
                        xxx.append(tuple(yyy))
                        profiles[tuple(xxx)].append([geocore.latitude, geocore.longitude, key])
    for materialpresent in profiles.keys():
        if profiles[materialpresent] == []:
            del profiles[materialpresent]       #Then i need to remove all keys in the dictionary that
                                                # point to an empty list just for efficiency
    for materialpresent in profiles.keys():
        lat = []
        lon = []
        removelist = []
        for location in profiles[materialpresent]:
            if (location[0] in lat) and (location[1] in lon):
                removelist.append(location)             # removing any residual duplicates at diff depths
            else:
                lat.append(location[0])
                lon.append(location[1])
        for ele in removelist:
            profiles[materialpresent].remove(ele)
    templist = []
    dellist = []
    depthlist = []
    for materialpresent in profiles.keys():
        if (len(profiles[materialpresent]) == 1):
            depthlist.append(profiles[materialpresent][0][2])
            if (((profiles[materialpresent][0][2] - 1) in depthlist) or (profiles[materialpresent][0][2] + 1) in depthlist):
                dellist.append(materialpresent)
        if (profiles[materialpresent] in templist):
            dellist.append(materialpresent)
        else:
            templist.append(profiles[materialpresent])
    for ele in dellist:
        del profiles[ele]
    return profiles

# DISTANCE FUNCTION
def distance(x1, y1, x2, y2):
    d = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    return d

# CREATING IDW FUNCTION
def idw_npoint(xz, yz, n_point, p, x, y, z):
    r = 10  # block radius iteration distance
    nf = 0
    counter = 0
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
        counter += 1
        if (counter > 1000):
            break

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

#--------------------------------------------------  Subterranean Estimation Portion  ----------------------------------

generateAllGeocores()
geocores = generateProfileDictionary()
listofx = []
listofy = []
listofz = []

# READING AND PARSING THE TERRAIN DATA
file = open(terraindata, 'r')
lines = file.readlines()
n_line = len(lines)
xterrain = []
yterrain = []
zterrain = []
for i in range(1, n_line):
    split_line = lines[i].split(",")
    xyz_t = []
    xterrain.append(float(split_line[1].rstrip()))
    yterrain.append(float(split_line[0].rstrip()))
    zterrain.append(float(split_line[2].rstrip()))

xt_min = min(xterrain)
xt_max = max(xterrain)
yt_min = min(yterrain)
yt_max = max(yterrain)
w = xt_max - xt_min  # width
h = yt_max - yt_min  # length
wn = w / n  # x interval
hn = h / n  # y interval

# list to store interpolation point and elevation
yt_init = yt_min
xt_init = xt_min
x_terrain = []
y_terrain = []
z_terrain = []
for i in range(n):
    xz = xt_init + wn * i
    yz = yt_init + hn * i
    y_terrain.append(yz)
    x_terrain.append(xz)
    z_idw_list = []
    for j in range(n):
        xz = xt_init + wn * j
        z_idw = idw_npoint(xz, yz, 5, 1.5, xterrain, yterrain, zterrain)  # min. point=5, p=1.5
        z_idw_list.append(z_idw)
    z_terrain.append(z_idw_list)

finnn = []
for setofpoints in geocores.values():
    # POPULATE INTERPOLATION POINTS
    tf = True
    for ele in setofpoints:
        if (ele[2] != 0):
            tf = False
    if (tf == True):
        for ele in geocores.keys():
            if geocores[ele] == setofpoints:
                finnn.append(ele)
        continue
    x = []
    y = []
    z = []
    for ele in setofpoints:
        x.append(ele[0])
        y.append(ele[1])
        z.append(ele[2])

    x_min = min(xterrain)
    x_max = max(xterrain)
    y_min = min(yterrain)
    y_max = max(yterrain)
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
            z_idw = idw_npoint(xz, yz, 5, 1.5, x, y, z)  # min. point=5, p=1.5
            z_idw_list.append(z_idw)
        z_new_list = []
        for ele in range(len(z_idw_list)):
            z_new_list.append(z_terrain[i][ele] - z_idw_list[ele])
        z_head.append(z_new_list)
    listofx.append(x_idw_list)
    listofy.append(y_idw_list)
    listofz.append(z_head)
for ele in finnn:
    del geocores[ele]


if not os.path.exists('InjectionAnalysis/Geocore'):
    os.makedirs('InjectionAnalysis/Geocore')
# CREATING 3D TERRAIN
fig = go.Figure()
fig.add_trace(go.Surface(z=z_terrain, x=x_terrain, y=y_terrain, colorscale='earth'))
for ele in range(len(listofz)):
    fig.add_trace(go.Surface(z=listofz[ele], x=listofx[ele], y=listofy[ele]))
for geoccore in geocoredata:
    fig.add_scatter3d(x=[geoccore.latitude, geoccore.latitude], y=[geoccore.longitude, geoccore.longitude],
                      z=[540, 610], mode='lines',
                      marker=dict(size=2,
                                  colorscale='Reds'))

fig.update_layout(
    scene=dict(aspectratio=dict(x=2, y=2, z=1), xaxis=dict(range=[min(xterrain), max(xterrain)], ), yaxis=dict(range=
                                                                                        [min(yterrain),
                                                                                         max(yterrain)])))
fig.update_layout(title_text="Terrain Interpolation")
go_offline.plot(fig, filename='InjectionAnalysis/Geocore/SubterraneanInterpolation.html',
                validate=True, auto_open=False)

#--------------------------------------------------  Saving Data Portion  ---------------------------------------------

if not os.path.exists('InjectionAnalysis/General/SubterraneanStructures'):
    os.makedirs('InjectionAnalysis/General/SubterraneanStructures')
listofz.insert(0, z_terrain)
#Saving to a csv file for future processing
fields = ['Structures', 'Latitude', 'Longitude', 'Elevation']
listtosave = [list(geocores), listofx, listofy, listofz]
file_name = "sample.pkl"

open_file = open('InjectionAnalysis/General/SubterraneanStructures/Structures.pkl', "wb")
pickle.dump(listtosave, open_file)
open_file.close()

