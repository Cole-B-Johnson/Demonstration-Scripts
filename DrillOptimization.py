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
