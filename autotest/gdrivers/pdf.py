#!/usr/bin/env python
###############################################################################
# $Id$
#
# Project:  GDAL/OGR Test Suite
# Purpose:  PDF Testing.
# Author:   Even Rouault <even dot rouault at mines dash paris dot org>
#
###############################################################################
# Copyright (c) 2010, Even Rouault <even dot rouault at mines dash paris dot org>
# 
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
###############################################################################

import os
import sys
import gdal
import ogr

sys.path.append( '../pymod' )

import gdaltest
import ogrtest

###############################################################################
# Test OGC best practice geospatial PDF

def pdf_online_1():

    try:
        if gdal.GetDriverByName('PDF') is None:
            return 'skip'
    except:
        return 'skip'

    if not gdaltest.download_file('http://www.agc.army.mil/GeoPDFgallery/Imagery/Cherrydale_eDOQQ_1m_0_033_R1C1.pdf', 'Cherrydale_eDOQQ_1m_0_033_R1C1.pdf'):
        return 'skip'

    try:
        os.stat('tmp/cache/Cherrydale_eDOQQ_1m_0_033_R1C1.pdf')
    except:
        return 'skip'

    ds = gdal.Open('tmp/cache/Cherrydale_eDOQQ_1m_0_033_R1C1.pdf')
    if ds is None:
        return 'fail'

    if ds.RasterXSize != 620:
        gdaltest.post_reason('bad dimensions')
        return 'fail'

    gt = ds.GetGeoTransform()
    wkt = ds.GetProjectionRef()

    expected_gt = (-77.112328333299999, 1.8333311999999999e-05, 0.0, 38.897842488372, -0.0, -1.8333311999999999e-05)
    for i in range(6):
        if abs(gt[i] - expected_gt[i]) > 1e-15:
            gdaltest.post_reason('bad geotransform')
            print(gt)
            print(expected_gt)
            return 'fail'

    expected_wkt = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]]'
    if wkt != expected_wkt:
        gdaltest.post_reason('bad WKT')
        print(wkt)
        return 'fail'

    cs = ds.GetRasterBand(1).Checksum()
    if cs == 0:
        gdaltest.post_reason('bad checksum')
        return 'fail'

    return 'success'

###############################################################################

def pdf_online_2():

    try:
        if gdal.GetDriverByName('PDF') is None:
            return 'skip'
    except:
        return 'skip'

    try:
        os.stat('tmp/cache/Cherrydale_eDOQQ_1m_0_033_R1C1.pdf')
    except:
        return 'skip'

    ds = gdal.Open('PDF:1:tmp/cache/Cherrydale_eDOQQ_1m_0_033_R1C1.pdf')
    if ds is None:
        return 'fail'

    gt = ds.GetGeoTransform()
    wkt = ds.GetProjectionRef()

    expected_gt = (-77.112328333299999, 1.8333311999999999e-05, 0.0, 38.897842488372, -0.0, -1.8333311999999999e-05)
    for i in range(6):
        if abs(gt[i] - expected_gt[i]) > 1e-15:
            print(gt)
            print(expected_gt)
            return 'fail'

    expected_wkt = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]]'
    if wkt != expected_wkt:
        print(wkt)
        return 'fail'

    return 'success'

###############################################################################
# Test Adobe style geospatial pdf

def pdf_1():

    try:
        if gdal.GetDriverByName('PDF') is None:
            return 'skip'
    except:
        return 'skip'

    gdal.SetConfigOption('GDAL_PDF_DPI', '200')
    ds = gdal.Open('data/adobe_style_geospatial.pdf')
    gdal.SetConfigOption('GDAL_PDF_DPI', None)
    if ds is None:
        return 'fail'

    gt = ds.GetGeoTransform()
    wkt = ds.GetProjectionRef()

    expected_gt = (333274.61654367246, 31.764802242655662, 0.0, 4940391.7593506984, 0.0, -31.794745501708238)
    for i in range(6):
        if abs(gt[i] - expected_gt[i]) > 1e-8:
            gdaltest.post_reason('bad geotransform')
            print(gt)
            print(expected_gt)
            return 'fail'

    expected_wkt = 'PROJCS["WGS_1984_UTM_Zone_20N",GEOGCS["GCS_WGS_1984",DATUM["WGS_1984",SPHEROID["WGS_84",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-63.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]'
    if wkt != expected_wkt:
        gdaltest.post_reason('bad WKT')
        print(wkt)
        return 'fail'

    cs = ds.GetRasterBand(1).Checksum()
    if cs != 17740:
        gdaltest.post_reason('bad checksum')
        print(cs)
        return 'fail'

    neatline = ds.GetMetadataItem('NEATLINE')
    got_geom = ogr.CreateGeometryFromWkt(neatline)
    expected_geom = ogr.CreateGeometryFromWkt('POLYGON ((338304.150125828920864 4896673.639421294443309,338304.177293475600891 4933414.799376524984837,382774.271384406310972 4933414.546264361590147,382774.767329963855445 4896674.273581005632877,338304.150125828920864 4896673.639421294443309))')

    if ogrtest.check_feature_geometry(got_geom, expected_geom) != 0:
        gdaltest.post_reason('bad neatline')
        print(neatline)
        return 'fail'

    return 'success'

gdaltest_list = [
    pdf_online_1,
    pdf_online_2,
    pdf_1 ]


if __name__ == '__main__':

    gdaltest.setup_run( 'PDF' )

    gdaltest.run_tests( gdaltest_list )

    gdaltest.summarize()

