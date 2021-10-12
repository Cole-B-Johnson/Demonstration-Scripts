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
                z_idw = idw_npoint(x, y, z, xz, yz, 5, 1.5)  # min. point=5, p=1.5
                z_idw_list.append(z_idw)
            z_head.append(z_idw_list)
        return x_idw_list, y_idw_list, z_head

def idw_npoint(x, y, z, xz, yz, n_point, p):
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