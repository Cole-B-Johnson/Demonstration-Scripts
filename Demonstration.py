import numpy as np
import plotly.graph_objects as go
from scipy.interpolate import griddata as gd

def readSurveyDataFromCSV(filename):
    surveydata = filename
    file = open(surveydata, 'r')
    lines = file.readlines()
    n_line = len(lines)
    x = []
    y = []
    z = []
    for i in range(1, n_line):
        split_line = lines[i].split(",")
        x.append(float(split_line[1].rstrip()))
        y.append(float(split_line[0].rstrip()))
        z.append(float(split_line[2].rstrip()))
    return x, y, z

def distance(x1, y1, x2, y2):
    d = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    return d

def idw_xpoint(x, y, z, xz, yz, n_point, p):
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


def interpolateSurface(x, y, z, extrapolation_interval=30, mode='cubic', resolutionfororiginal=50):
    # modes: 'nearest', 'linear', 'cubic', 'original'
    if mode != 'original':
        converted = []
        for ele in range(len(z)):
            converted.append([x[ele], y[ele], z[ele]])
        converted = np.array(converted)
        extrapolation_spots = get_plane(min(y), max(y), min(x), max(x), extrapolation_interval)
        return nearest_analysis(extrapolation_spots, converted, mode, x, y)
    else:
        n = resolutionfororiginal
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
                z_idw = idw_xpoint(x, y, z, xz, yz, 5, 1.5)  # min. point=5, p=1.5
                z_idw_list.append(z_idw)
            z_head.append(z_idw_list)
        return x_idw_list, y_idw_list, z_head

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

def nearest_analysis(extrapolation_spots, converted, mode, x, y):
    bot_extra = extrapolation(converted, extrapolation_spots)
    gridx_top, gridy_top, gridz_top = interpolation(bot_extra, mode, x, y)
    newx = []
    newy = []
    for ele in gridx_top:
        newx.append(ele[0])
    for ele in gridy_top[0]:
        newy.append(ele)
    return newx, newy, gridz_top

def nearest_neighbor_interpolation(data, x, y, p=0.5):
    n = len(data)
    vals = np.zeros((n, 2), dtype=np.float64)
    distance = lambda x1, x2, y1, y2: (x2 - x1)**2 + (y2 - y1)**2
    for i in range(n):
        vals[i, 0] = data[i, 2] / (distance(data[i, 0], x, data[i, 1], y))**p
        vals[i, 1] = 1          / (distance(data[i, 0], x, data[i, 1], y))**p
    z = np.sum(vals[:, 0]) / np.sum(vals[:, 1])
    return z

def get_plane(xl, xu, yl, yu, i):
    xx = np.arange(xl, xu, i)
    yy = np.arange(yl, yu, i)
    extrapolation_spots = np.zeros((len(xx) * len(yy), 2))
    count = 0
    for i in xx:
        for j in yy:
            extrapolation_spots[count, 0] = i
            extrapolation_spots[count, 1] = j
            count += 1
    return extrapolation_spots

def extrapolation(data, extrapolation_spots):
    new_points = np.zeros((len(extrapolation_spots), 3))
    new_points[:, 0] = extrapolation_spots[:, 0]
    new_points[:, 1] = extrapolation_spots[:, 1]
    for i in range(len(extrapolation_spots)):
        new_points[i, 2] = nearest_neighbor_interpolation(data,
                                                          extrapolation_spots[i, 0],
                                                          extrapolation_spots[i, 1], p=2)
    combined = np.concatenate((data, new_points))
    return combined

def interpolation(data, mode, x, y):
    gridx, gridy = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]
    gridz = gd(data[:, :2],data[:, 2], (gridx, gridy),
               method=mode)
    return gridx, gridy, gridz
import pickle
import math
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import plotly.offline as go_offline

# ------------------------------------------------- LOADING IN DATA FROM SUBTERRANEAN ANALYSIS ---------------------

open_file = open('InjectionAnalysis/General/SubterraneanStructures/Structures.pkl', "rb")
[composition, listofx, listofy, listofz] = pickle.load(open_file)
fakex = listofx
fakey = listofy
fakez = listofz
listofx = listofx[0]
listofy = listofy[0]
open_file.close()
vvoidfractions = []
voidfractions = []
particlesizes = []
for elem in composition:
    vvoidfractions.append(elem[0])
for elem in composition:
    voidfractions.append(elem[0] / max(vvoidfractions))  # REMOVE MAX STATEMENT WHEN ACTUALLY USING DATA
    particlesizes.append(elem[1])

# ------------------------------------------------- TUNING PARAMETERS FOR PENETRATION RESISTANCE ---------------------

k = 12  # for steepness of logistic curve for void fraction resistance
thresholdsresistance = {1: .5, 2: 1, 5: 2}  # clay, silt, and sand, respectively: {threshold: resistance, ...}
voidweight = 2
particleweight = 1
initialresistance = .1  # the resistance of the surface of the terrain as given by the penetrometer
maxdepth = 100  # the max depth that will be drilled to throughout injection

# ------------------------------------------------- TO GET RESISTANCE VALUES ------------------------------------------


# LOGISTIC CURVE CALCULATOR
def logisticCurve(voidfraction):
    weight = 1 / (1 + math.e ** (-1 * k * (voidfraction - .5)))
    return weight


# TO CALCULATE WEIGHT OF RESISTANCE FROM VOID FRACTION
def voidResistance(voidfraction):
    logg = logisticCurve(voidfraction)
    return 1 - logg


# TO GET FRACTION OF PARTICLE SIZES IN DISTRIBUTION THAT ARE BELOW THRESHOLD
def fractionOfParticles(threshold, particles):
    count = 0
    for particlesize in particles:
        if particlesize[1] >= threshold:
            count += 1
    return count / len(particles)


# TO GET OVERALL CHARACTERIZATION OF PARTICLE SIZE DISTRIBUTION
def particleCharacterization(particles):
    top = 0
    final = 0
    for threshold in list(thresholdsresistance):
        fin = fractionOfParticles(threshold, particles)
        if top < fin:
            top = fin
            final = threshold
    return final


# TO GET RESISTANCE FROM PARTICLE SIZE
def particleResistance(particles):
    return thresholdsresistance[particleCharacterization(particles)]


def getResistance(void, particles):
    return (voidweight * voidResistance(void)) + (particleweight * particleResistance(particles))


def getAllResistance():
    fin = [initialresistance]
    for ele in range(len(voidfractions)):
        fin.append(getResistance(voidfractions[ele], particlesizes[ele]))
    return fin


# ------------------------------------------- TO GET 3D MODEL OF UNDERGROUND ------------------------------------------

resistancelayers = getAllResistance()
numLayers = len(resistancelayers)

maxterrain = 0
for elem in listofz:
    for eleme in elem:
        for elemee in eleme:
            if elemee > maxterrain:
                maxterrain = elemee
basearray = [[[0 for elee in range(len(listofx))] for eleee in range(len(listofy))] for eleeee in
             range(math.ceil(maxterrain - maxdepth))]

# GET RESISTANCE OF GIVEN ELEMENT IN 3D MATRIX
def getResistanceAbove(coordinates):
    x, y, z = coordinates
    posofx = 0
    posofy = 0
    for ele in range(len(listofx)):
        if listofx[ele] == x:
            posofx = ele
        if listofy[ele] == y:
            posofy = ele

    correspondingz = {}
    for ele in range(len(listofz)):
        correspondingz[listofz[ele][posofx][posofy]] = resistancelayers[ele]
    maxx = max(list(correspondingz))
    for ele in correspondingz.keys():
        if z > ele:
            continue
        else:
            if maxx > ele:
                maxx = ele
    return correspondingz[maxx]


# SET ALL ELEMENTS IN BASEARRAY TO RESISTANCE MOST DIRECTLY ABOVE
def setEntireMatrix():
    for colentry in range(len(basearray)):
        for row in range(len(basearray[colentry])):
            for length in range(len(basearray[colentry][row])):
                basearray[colentry][row][length] = getResistanceAbove((listofx[row], listofy[length],
                                                                       maxdepth + colentry))
    return basearray

# ------------------------------------------- TO FIND OPTIMAL DRILLING PATH ------------------------------------------






















# ------------------------------------------------------ GRAPHICS GENERATION ------------------------------------------
def mockup():
    fig = go.Figure()
    for ele in range(1, len(fakex)):
        fig.add_trace(go.Surface(z=fakez[ele], x=fakex[ele], y=fakey[ele], colorscale='viridis'))
    fig.update_traces(opacity=.8, selector=dict(type='surface'))
    fig.add_trace(go.Surface(z=fakez[0], x=fakex[0], y=fakey[0]))

    fig.add_scatter3d(x=[12800, 12800, 12805, 12810, 12820, 12845, 13380, 13390, 13400, 13700],
                      y=[9600, 9600, 9605, 9610, 9620, 9645, 10180, 10190, 10200, 10500],
                      z=[610, 600, 598, 597, 596, 594.5, 555.7, 555.3, 555, 555], line=dict(color='red', width=10),
                      marker=dict(size=2, colorscale='Reds'))


    fig.update_layout(
        scene=dict(aspectratio=dict(x=2, y=2, z=1),
                   xaxis=dict(range=[min(fakex[0]), max(fakex[0])]), yaxis=dict(range=[min(fakey[0]),max(fakey[0])])))
    fig.update_layout(title_text="Drilling Optimization")
    go_offline.plot(fig, filename='InjectionAnalysis/Geocore/DrillingOptimization.html',
                    validate=True, auto_open=False)
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
def idw_xpoint(xz, yz, n_point, p, x, y, z):
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
        z_idw = idw_xpoint(xz, yz, 5, 1.5, xterrain, yterrain, zterrain)  # min. point=5, p=1.5
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
            z_idw = idw_xpoint(xz, yz, 5, 1.5, x, y, z)  # min. point=5, p=1.5
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

#--------------------------------------------------  Terrain Estimation Portion  -------------------------------------

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
