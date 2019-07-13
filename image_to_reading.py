#!/usr/bin/env python

from osgeo import gdal
import numpy as np
import math
import os

pos = {}

with open("tram-pos.csv", 'r') as positions:
    flines = positions.readlines()
    first = True
    for l in flines:
        if first:
            first = False
            continue
        l = l.split(',')
        pos[l[0]] = (float(l[1]), float(l[2]))

outer_folder = {"Temp": ["MOD06", "MOD07", "MYD06", "MYD07", "VIIRS"],
                "Pressure": ["MOD06", "MOD07", "MYD06", "MYD07"],
                "WaterVapor": ["MOD05", "MOD07", "MYD05", "MYD07"]}        

pix_pos = {}
obs = pos.keys()

for od in os.listdir('.'):
    if od not in outer_folder:
        continue
    pref = outer_folder[od]
    data = open(od+"_data.csv", 'w+')
    infos = [["" for x in range(len(obs))] for y in range(31)]
    colname = 'name,day'

    for d in os.listdir(od):
        if d not in pref:
            continue
        current_folder = od + "/" + d
        print current_folder
        colname += ','
        colname += d
        for day in range(1, 32):
            suffix = str("%03d" % day)
            name = current_folder + "/test_im_" + suffix + ".tif"
            im = gdal.Open(name, gdal.GA_ReadOnly)
            #geo[0] = ulx, geo[1] = s1, geo[3] = uly, geo[5] = s2
            band = im.GetRasterBand(1)
            a = band.ReadAsArray().astype(np.float)

            for k in range(len(obs)):
                key = obs[k]
                if key not in pix_pos.keys():
                    geo = im.GetGeoTransform()
                    coord = pos[key]
                    pix_pos[key] = (int(math.floor((coord[0] - geo[3]) / geo[5])),
                                    int(math.floor((coord[1] - geo[0]) / geo[1])))
                infos[day - 1][k] += ","
                infos[day - 1][k] += str(a[pix_pos[key][0]][pix_pos[key][1]])

    data.write(colname + '\n')
    for i in range(len(obs)):
        tr = obs[i]
        for j in range(1, 32):
            trd = tr + ',' + str(j)
            data.write(trd + infos[j-1][i] + '\n')

print "Done"
data.close()
