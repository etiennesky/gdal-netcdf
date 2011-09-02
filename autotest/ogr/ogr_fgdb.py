#!/usr/bin/env python
###############################################################################
# $Id$
#
# Project:  GDAL/OGR Test Suite
# Purpose:  FGDB driver testing.
# Author:   Even Rouault <even dot rouault at mines dash paris dot org>
#
###############################################################################
# Copyright (c) 2011, Even Rouault <even dot rouault at mines dash paris dot org>
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
import string
import shutil

sys.path.append( '../pymod' )

import gdaltest
import ogrtest
import ogr
import osr

###############################################################################
# Test if driver is available

def ogr_fgdb_init():

    ogrtest.fgdb_drv = None

    try:
        ogrtest.fgdb_drv = ogr.GetDriverByName('FileGDB')
    except:
        pass

    if ogrtest.fgdb_drv is None:
        return 'skip'

    try:
        shutil.rmtree("tmp/test.gdb")
    except:
        pass

    return 'success'

###############################################################################
# Write and read back various geometry types

def ogr_fgdb_1():
    if ogrtest.fgdb_drv is None:
        return 'skip'

    srs = osr.SpatialReference()
    srs.SetFromUserInput("WGS84")

    ds = ogrtest.fgdb_drv.CreateDataSource("tmp/test.gdb")

    datalist = [ [ "point", ogr.wkbPoint, "POINT (1 2)" ],
                 [ "multipoint", ogr.wkbMultiPoint, "MULTIPOINT (1 2,3 4)" ],
                 [ "linestring", ogr.wkbLineString, "LINESTRING (1 2,3 4)", "MULTILINESTRING ((1 2,3 4))" ],
                 [ "multilinestring", ogr.wkbMultiLineString, "MULTILINESTRING ((1 2,3 4))" ],
                 [ "polygon", ogr.wkbPolygon, "POLYGON ((0 0,0 1,1 1,1 0,0 0))", "MULTIPOLYGON (((0 0,0 1,1 1,1 0,0 0)))" ],
                 [ "multipolygon", ogr.wkbMultiPolygon, "MULTIPOLYGON (((0 0,0 1,1 1,1 0,0 0)))" ],
                 [ "point25D", ogr.wkbPoint25D, "POINT (1 2 3)" ],
                 [ "multipoint25D", ogr.wkbMultiPoint25D, "MULTIPOINT (1 2 -10,3 4 -20)" ],
                 [ "linestring25D", ogr.wkbLineString25D, "LINESTRING (1 2 -10,3 4 -20)", "MULTILINESTRING ((1 2 -10,3 4 -20))" ],
                 [ "multilinestring25D", ogr.wkbMultiLineString25D, "MULTILINESTRING ((1 2 -10,3 4 -20))" ],
                 [ "polygon25D", ogr.wkbPolygon25D, "POLYGON ((0 0 -10,0 1 -10,1 1 -10,1 0 -10,0 0 -10))", "MULTIPOLYGON (((0 0 -10,0 1 -10,1 1 -10,1 0 -10,0 0 -10)))" ],
                 [ "multipolygon25D", ogr.wkbMultiPolygon25D, "MULTIPOLYGON (((0 0 -10,0 1 -10,1 1 -10,1 0 -10,0 0 -10)))" ],
               ]

    for data in datalist:
        lyr = ds.CreateLayer(data[0], geom_type = data[1], srs = srs)
        lyr.CreateField(ogr.FieldDefn("str", ogr.OFTString))
        lyr.CreateField(ogr.FieldDefn("int", ogr.OFTInteger))
        lyr.CreateField(ogr.FieldDefn("real", ogr.OFTReal))
        feat = ogr.Feature(lyr.GetLayerDefn())
        feat.SetGeometry(ogr.CreateGeometryFromWkt(data[2]))
        feat.SetField("str", "foo_\xc3\xa9")
        feat.SetField("int", 123)
        feat.SetField("real", 4.56)
        lyr.CreateFeature(feat)

    for data in datalist:
        lyr = ds.GetLayerByName(data[0])
        if lyr.GetSpatialRef().IsSame(srs) != 1:
            print(lyr.GetSpatialRef())
            return 'fail'
        feat = lyr.GetNextFeature()
        try:
            expected_wkt = data[3]
        except:
            expected_wkt = data[2]
        if feat.GetGeometryRef().ExportToWkt() != expected_wkt:
            feat.DumpReadable()
            return 'fail'

    ds = None

    return 'success'

###############################################################################
# Run test_ogrsf

def ogr_fgdb_2():
    if ogrtest.fgdb_drv is None:
        return 'skip'

    import test_cli_utilities
    if test_cli_utilities.get_test_ogrsf_path() is None:
        return 'skip'

    ret = gdaltest.runexternal(test_cli_utilities.get_test_ogrsf_path() + ' -ro tmp/test.gdb')

    if ret.find('INFO') == -1 or ret.find('ERROR') != -1:
        print(ret)
        return 'fail'

    return 'success'

###############################################################################
# Run ogr2ogr

def ogr_fgdb_3():
    if ogrtest.fgdb_drv is None:
        return 'skip'

    import test_cli_utilities
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'
    if test_cli_utilities.get_test_ogrsf_path() is None:
        return 'skip'

    try:
        shutil.rmtree("tmp/poly.gdb")
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -f filegdb tmp/poly.gdb data/poly.shp -nlt MULTIPOLYGON -a_srs None')

    ret = gdaltest.runexternal(test_cli_utilities.get_test_ogrsf_path() + ' tmp/poly.gdb')
    #print ret

    if ret.find('INFO') == -1 or ret.find('ERROR') != -1:
        print(ret)
        return 'fail'

    return 'success'

###############################################################################
# Test delete layer

def ogr_fgdb_4():
    if ogrtest.fgdb_drv is None:
        return 'skip'

    for j in range(2):

        # Create a layer
        ds = ogr.Open("tmp/test.gdb", update = 1)
        srs = osr.SpatialReference()
        srs.SetFromUserInput("WGS84")
        lyr = ds.CreateLayer("layer_to_remove", geom_type = ogr.wkbPoint, srs = srs)
        lyr.CreateField(ogr.FieldDefn("str", ogr.OFTString))
        feat = ogr.Feature(lyr.GetLayerDefn())
        feat.SetGeometry(ogr.CreateGeometryFromWkt('POINT(2 49)'))
        feat.SetField("str", "foo")
        feat = None
        lyr = None

        if j == 1:
            ds = None
            ds = ogr.Open("tmp/test.gdb", update = 1)

        # Delete it
        for i in range(ds.GetLayerCount()):
            if ds.GetLayer(i).GetName() == 'layer_to_remove':
                ds.DeleteLayer(i)
                break

        # Check it no longer exists
        lyr = ds.GetLayerByName('layer_to_remove')
        ds = None

        if lyr is not None:
            gdaltest.post_reason('failed at iteration %d' % j)
            return 'fail'

    return 'success'

###############################################################################
# Cleanup

def ogr_fgdb_cleanup():
    if ogrtest.fgdb_drv is None:
        return 'skip'

    shutil.rmtree("tmp/test.gdb")

    try:
        shutil.rmtree("tmp/poly.gdb")
    except:
        pass

    return 'success'

gdaltest_list = [
    ogr_fgdb_init,
    ogr_fgdb_1,
    ogr_fgdb_2,
    ogr_fgdb_3,
    ogr_fgdb_4,
    ogr_fgdb_cleanup,
    ]

if __name__ == '__main__':

    gdaltest.setup_run( 'ogr_fgdb' )

    gdaltest.run_tests( gdaltest_list )

    gdaltest.summarize()



