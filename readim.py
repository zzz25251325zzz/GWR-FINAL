#!/usr/bin/env python

import os
import fnmatch
from osgeo import gdal
import numpy as np
import time

outer_folder = {"Temp": [("MOD06", 1, "MOD06_L2.A2014"), ("MOD07", 1, "MOD07_L2.A2014"), ("MYD06", 1, "MYD06_L2.A2014"), ("MYD07", 1, "MYD07_L2.A2014"), ("VIIRS", 0, "GMTCO-VLSTO_npp_d201401")],
                "Pressure": [("MOD06", 1, "MOD06_L2.A2014"), ("MOD07", 1, "MOD07_L2.A2014"), ("MYD06", 1, "MYD06_L2.A2014"), ("MYD07", 1, "MYD07_L2.A2014")],
                "WaterVapor": [("MOD05", 1, "MOD05_L2.A2014"), ("MOD07", 1, "MOD07_L2.A2014"), ("MYD05", 1, "MYD05_L2.A2014"), ("MYD07", 1, "MYD07_L2.A2014")]}

for od in os.listdir('.'): 
    if od not in outer_folder:
        continue
    pref = outer_folder[od]
    for d in os.listdir(od):
        cor = []
        for p in pref:
            if d == p[0]:
                cor = p
                break
        if len(cor) == 0:
            continue

        current_folder = od + "/" + d    
        print current_folder
        for day in range(1, 32):
            pics = []
            unknown = -9999
            driver = None
            dtype = None
            refim = None
            
            name_suff = str("%03d" % day)
            if (cor[1] == 1):
                suffix = str("%03d" % day)
            else:
                suffix = str("%02d" % day)

            for f in os.listdir(current_folder):
                name = cor[2] + suffix + "*.tif"
                if fnmatch.fnmatch(f, name):
                    im = gdal.Open(current_folder + "/" + f, gdal.GA_ReadOnly)
                    driver = im.GetDriver()
                    band = im.GetRasterBand(1)
                    dtype = band.DataType
                    refim = im
                    unknown = band.GetNoDataValue()
                    a = band.ReadAsArray().astype(np.float)
                    pics.append(a)

            assemble = np.asarray(pics)

            count = np.sum(assemble != unknown, axis = 0)
            (h, w) = count.shape
            total = np.sum(assemble, axis = 0)
            aver = np.ones(count.shape) * unknown
            maxnum = len(pics)
            u = (count == 0)
            nu = (count != 0)
            total = np.multiply(nu, total)
            count = count + u * maxnum
            aver = u * unknown + nu * (total - (maxnum - count) * unknown) / count

            outname = current_folder + "/test_im_" + name_suff + ".tif"
            output = driver.Create(outname, w, h, 1, dtype)
            output.SetGeoTransform(refim.GetGeoTransform())
            output.SetProjection(refim.GetProjection())
            outband = output.GetRasterBand(1)
            outband.SetNoDataValue(unknown)
            outband.WriteArray(aver)
            outband.FlushCache()

print "Done"
