#!/usr/bin/env python
###############################################################################
# $Id$
#
# Project:  GDAL/OGR Test Suite
# Purpose:  ogr2ogr testing
# Author:   Even Rouault <even dot rouault @ mines-paris dot org>
# 
###############################################################################
# Copyright (c) 2008, Even Rouault <even dot rouault @ mines-paris dot org>
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

import sys
import os
import shutil

sys.path.append( '../pymod' )
sys.path.append( '../ogr' )

import gdal
import ogr
import osr
import gdaltest
import ogrtest
import test_cli_utilities

###############################################################################
# Simple test

def test_ogr2ogr_1():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/poly.shp ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        return 'fail'

    feat0 = ds.GetLayer(0).GetFeature(0)
    if feat0.GetFieldAsDouble('AREA') != 215229.266:
        print(feat0.GetFieldAsDouble('AREA'))
        gdaltest.post_reason('Did not get expected value for field AREA')
        return 'fail'
    if feat0.GetFieldAsString('PRFEDEA') != '35043411':
        print(feat0.GetFieldAsString('PRFEDEA'))
        gdaltest.post_reason('Did not get expected value for field PRFEDEA')
        return 'fail'
        
    ds.Destroy()
    
    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -sql

def test_ogr2ogr_2():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/poly.shp ../ogr/data/poly.shp -sql "select * from poly"')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -spat

def test_ogr2ogr_3():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/poly.shp ../ogr/data/poly.shp -spat 479609 4764629 479764 4764817')

    ds = ogr.Open('tmp/poly.shp')
    if ogrtest.have_geos():
        if ds is None or ds.GetLayer(0).GetFeatureCount() != 4:
            return 'fail'
    else:
        if ds is None or ds.GetLayer(0).GetFeatureCount() != 5:
            return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -where

def test_ogr2ogr_4():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/poly.shp ../ogr/data/poly.shp -where "EAS_ID=171"')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 1:
        return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -append

def test_ogr2ogr_5():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/poly.shp ../ogr/data/poly.shp')
    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -update -append tmp/poly.shp ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 20:
        return 'fail'
        
    feat10 = ds.GetLayer(0).GetFeature(10)
    if feat10.GetFieldAsDouble('AREA') != 215229.266:
        print(feat10.GetFieldAsDouble('AREA'))
        gdaltest.post_reason('Did not get expected value for field AREA')
        return 'fail'
    if feat10.GetFieldAsString('PRFEDEA') != '35043411':
        print(feat10.GetFieldAsString('PRFEDEA'))
        gdaltest.post_reason('Did not get expected value for field PRFEDEA')
        return 'fail'
        
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -overwrite

def test_ogr2ogr_6():

    import ogr_pg

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'
    if test_cli_utilities.get_ogrinfo_path() is None:
        return 'skip'

    ogr_pg.ogr_pg_1()
    if gdaltest.pg_ds is None:
        return 'skip'
    gdaltest.pg_ds.Destroy()

    gdaltest.runexternal(test_cli_utilities.get_ogrinfo_path() + ' PG:"' + gdaltest.pg_connection_string + '" -sql "DELLAYER:tpoly"')

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -f PostgreSQL PG:"' + gdaltest.pg_connection_string + '" ../ogr/data/poly.shp -nln tpoly')
    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -update -overwrite -f PostgreSQL PG:"' + gdaltest.pg_connection_string + '" ../ogr/data/poly.shp -nln tpoly')

    ds = ogr.Open('PG:' + gdaltest.pg_connection_string)
    if ds is None or ds.GetLayerByName('tpoly').GetFeatureCount() != 10:
        return 'fail'
    ds.Destroy()

    gdaltest.runexternal(test_cli_utilities.get_ogrinfo_path() + ' PG:"' + gdaltest.pg_connection_string + '" -sql "DELLAYER:tpoly"')

    return 'success'

###############################################################################
# Test -gt

def test_ogr2ogr_7():

    import ogr_pg

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'
    if test_cli_utilities.get_ogrinfo_path() is None:
        return 'skip'

    ogr_pg.ogr_pg_1()
    if gdaltest.pg_ds is None:
        return 'skip'
    gdaltest.pg_ds.Destroy()

    gdaltest.runexternal(test_cli_utilities.get_ogrinfo_path() + ' PG:"' + gdaltest.pg_connection_string + '" -sql "DELLAYER:tpoly"')

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -f PostgreSQL PG:"' + gdaltest.pg_connection_string + '" ../ogr/data/poly.shp -nln tpoly -gt 1')

    ds = ogr.Open('PG:' + gdaltest.pg_connection_string)
    if ds is None or ds.GetLayerByName('tpoly').GetFeatureCount() != 10:
        return 'fail'
    ds.Destroy()

    gdaltest.runexternal(test_cli_utilities.get_ogrinfo_path() + ' PG:"' + gdaltest.pg_connection_string + '" -sql "DELLAYER:tpoly"')

    return 'success'

###############################################################################
# Test -t_srs

def test_ogr2ogr_8():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -t_srs EPSG:4326 tmp/poly.shp ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/poly.shp')
    if str(ds.GetLayer(0).GetSpatialRef()).find('1984') == -1:
        return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -a_srs

def test_ogr2ogr_9():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -a_srs EPSG:4326 tmp/poly.shp ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/poly.shp')
    if str(ds.GetLayer(0).GetSpatialRef()).find('1984') == -1:
        return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -select

def test_ogr2ogr_10():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -select AREA tmp/poly.shp ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/poly.shp')
    if ds.GetLayer(0).GetLayerDefn().GetFieldCount() != 1:
        return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -lco

def test_ogr2ogr_11():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -lco SHPT=POLYGONZ tmp/poly.shp ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/poly.shp')
    if ds.GetLayer(0).GetLayerDefn().GetGeomType() != ogr.wkbPolygon25D:
        return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -nlt

def test_ogr2ogr_12():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -nlt POLYGON25D tmp/poly.shp ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/poly.shp')
    if ds.GetLayer(0).GetLayerDefn().GetGeomType() != ogr.wkbPolygon25D:
        return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Add explicit source layer name

def test_ogr2ogr_13():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/poly.shp ../ogr/data/poly.shp poly')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -segmentize

def test_ogr2ogr_14():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -segmentize 100 tmp/poly.shp ../ogr/data/poly.shp poly')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        return 'fail'
    feat = ds.GetLayer(0).GetNextFeature()
    if feat.GetGeometryRef().GetGeometryRef(0).GetPointCount() != 36:
        return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -overwrite with a shapefile

def test_ogr2ogr_15():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/poly.shp ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        return 'fail'
    ds.Destroy()

    # Overwrite
    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -overwrite tmp ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    return 'success'

###############################################################################
# Test -fid

def test_ogr2ogr_16():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -fid 8 tmp/poly.shp ../ogr/data/poly.shp')

    src_ds = ogr.Open('../ogr/data/poly.shp')
    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 1:
        return 'fail'
    src_feat = src_ds.GetLayer(0).GetFeature(8)
    feat = ds.GetLayer(0).GetNextFeature()
    if feat.GetField("EAS_ID") != src_feat.GetField("EAS_ID"):
        return 'fail'
    ds.Destroy()
    src_ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    return 'success'

###############################################################################
# Test -progress

def test_ogr2ogr_17():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    ret = gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -progress tmp/poly.shp ../ogr/data/poly.shp')
    if ret.find('0...10...20...30...40...50...60...70...80...90...100 - done.') == -1:
        return 'fail'

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    return 'success'

###############################################################################
# Test -wrapdateline

def test_ogr2ogr_18():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'
        
    if ogrtest.have_geos() is 0:
        return 'skip'

    try:
        os.stat('tmp/wrapdateline_src.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/wrapdateline_src.shp')
    except:
        pass
        
    try:
        os.stat('tmp/wrapdateline_dst.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/wrapdateline_dst.shp')
    except:
        pass
        
    ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('tmp/wrapdateline_src.shp')
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32660);
    lyr = ds.CreateLayer('wrapdateline_src', srs = srs)
    feat = ogr.Feature(lyr.GetLayerDefn())
    geom = ogr.CreateGeometryFromWkt('POLYGON((700000 4000000,800000 4000000,800000 3000000,700000 3000000,700000 4000000))')
    feat.SetGeometryDirectly(geom)
    lyr.CreateFeature(feat)
    feat.Destroy()
    ds.Destroy()
    
    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -wrapdateline -t_srs EPSG:4326 tmp/wrapdateline_dst.shp tmp/wrapdateline_src.shp')
    
    expected_wkt = 'MULTIPOLYGON (((179.222391385437419 36.124095832129363,180.0 36.10605558800065,180.0 27.090340569400169,179.017505655195095 27.107979523625211,179.222391385437419 36.124095832129363)),((-180.0 36.10605558800065,-179.667822828781084 36.098349195413753,-179.974688335419557 27.089886143076747,-180.0 27.090340569400169,-180.0 36.10605558800065)))'
    expected_geom = ogr.CreateGeometryFromWkt(expected_wkt)
    ds = ogr.Open('tmp/wrapdateline_dst.shp')
    lyr = ds.GetLayer(0)
    feat = lyr.GetNextFeature()
    ret = ogrtest.check_feature_geometry(feat, expected_geom)
    feat.Destroy()
    expected_geom.Destroy()
    ds.Destroy()
    
    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/wrapdateline_src.shp')
    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/wrapdateline_dst.shp')

    if ret == 0:
        return 'success'
    else:
        return 'fail'

###############################################################################
# Test -clipsrc

def test_ogr2ogr_19():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'
        
    if not ogrtest.have_geos():
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/poly.shp ../ogr/data/poly.shp -clipsrc spat_extent -spat 479609 4764629 479764 4764817')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 4:
        return 'fail'
        
    if ds.GetLayer(0).GetExtent() != (479609, 479764, 4764629, 4764817):
        print(ds.GetLayer(0).GetExtent())
        gdaltest.post_reason('unexpected extent')
        return 'fail'
        
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'
    
###############################################################################
# Test correct remap of fields when laundering to Shapefile format
# Test that the data is going into the right field
# FIXME: Any field is skipped if a subsequent field with same name is found.

def test_ogr2ogr_20():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    expected_fields = [ 'a',
                        'A_1',
                        'a_1_2',
                        'aaaaaAAAAA',
                        'aAaaaAAA_1',
                        'aaaaaAAAAB',
                        'aaaaaAAA_2',
                        'aaaaaAAA_3',
                        'aaaaaAAA_4',
                        'aaaaaAAA_5',
                        'aaaaaAAA_6',
                        'aaaaaAAA_7',
                        'aaaaaAAA_8',
                        'aaaaaAAA_9',
                        'aaaaaAAA10' ]
    expected_data = [ '1',
                      '2',
                      '3',
                      '4',
                      '5',
                      '6',
                      '7',
                      '8',
                      '9',
                      '10',
                      '11',
                      '12',
                      '13',
                      '14',
                      '15' ]

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp data/Fields.csv')
    
    ds = ogr.Open('tmp/Fields.dbf')

    if ds is None:
        return 'fail'
    layer_defn = ds.GetLayer(0).GetLayerDefn()
    if layer_defn.GetFieldCount() != 15:
        gdaltest.post_reason('Unexpected field count: ' + str(ds.GetLayer(0).GetLayerDefn().GetFieldCount()) )
        ds.Destroy()
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/Fields.dbf')
        return 'fail'

    error_occured = False
    feat = ds.GetLayer(0).GetNextFeature()    
    for i in range( layer_defn.GetFieldCount() ):
        if layer_defn.GetFieldDefn( i ).GetNameRef() != expected_fields[i]:
            print('Expected ', expected_fields[i],',but got',layer_defn.GetFieldDefn( i ).GetNameRef())
            error_occured = True
        if feat.GetFieldAsString(i) != expected_data[i]:
            print('Expected the value ', expected_data[i],',but got',feat.GetFieldAsString(i))
            error_occured = True

    ds.Destroy()
    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/Fields.dbf')

    if error_occured:
        return 'fail'

    return 'success'
    
###############################################################################
# Test ogr2ogr when the output driver has already created the fields
# at dataset creation (#3247)

def test_ogr2ogr_21():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.remove('tmp/testogr2ogr21.gtm')
    except:
        pass
        
    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() +
        ' -f GPSTrackMaker tmp/testogr2ogr21.gtm data/dataforogr2ogr21.csv ' +
        '-sql "SELECT comment, name FROM dataforogr2ogr21" -nlt POINT')
    ds = ogr.Open('tmp/testogr2ogr21.gtm')

    if ds is None:
        return 'fail'
    layer_defn = ds.GetLayer(0).GetLayerDefn()
    lyr = ds.GetLayer(0)
    feat = lyr.GetNextFeature()
    if feat.GetFieldAsString('name') != 'NAME' or \
       feat.GetFieldAsString('comment') != 'COMMENT':
        print(feat.GetFieldAsString('name'))
        print(feat.GetFieldAsString('comment'))
        ds.Destroy()
        os.remove('tmp/testogr2ogr21.gtm')
        return 'fail'

    ds.Destroy()
    os.remove('tmp/testogr2ogr21.gtm')

    return 'success'

    
###############################################################################
# Test ogr2ogr when the output driver delays the destination layer defn creation (#3384)

def test_ogr2ogr_22():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'
        
    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() +
        ' -f "MapInfo File" tmp/testogr2ogr22.mif data/dataforogr2ogr21.csv ' +
        '-sql "SELECT comment, name FROM dataforogr2ogr21" -nlt POINT')
    ds = ogr.Open('tmp/testogr2ogr22.mif')

    if ds is None:
        return 'fail'
    layer_defn = ds.GetLayer(0).GetLayerDefn()
    lyr = ds.GetLayer(0)
    feat = lyr.GetNextFeature()
    if feat.GetFieldAsString('name') != 'NAME' or \
       feat.GetFieldAsString('comment') != 'COMMENT':
        print(feat.GetFieldAsString('name'))
        print(feat.GetFieldAsString('comment'))
        ds.Destroy()
        ogr.GetDriverByName('MapInfo File').DeleteDataSource('tmp/testogr2ogr22.mif')
        return 'fail'

    ds.Destroy()
    ogr.GetDriverByName('MapInfo File').DeleteDataSource('tmp/testogr2ogr22.mif')

    return 'success'

###############################################################################
# Same as previous but with -select

def test_ogr2ogr_23():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'
        
    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() +
        ' -f "MapInfo File" tmp/testogr2ogr23.mif data/dataforogr2ogr21.csv ' +
        '-sql "SELECT comment, name FROM dataforogr2ogr21" -select comment,name -nlt POINT')
    ds = ogr.Open('tmp/testogr2ogr23.mif')

    if ds is None:
        return 'fail'
    layer_defn = ds.GetLayer(0).GetLayerDefn()
    lyr = ds.GetLayer(0)
    feat = lyr.GetNextFeature()
    if feat.GetFieldAsString('name') != 'NAME' or \
       feat.GetFieldAsString('comment') != 'COMMENT':
        print(feat.GetFieldAsString('name'))
        print(feat.GetFieldAsString('comment'))
        ds.Destroy()
        ogr.GetDriverByName('MapInfo File').DeleteDataSource('tmp/testogr2ogr23.mif')
        return 'fail'

    ds.Destroy()
    ogr.GetDriverByName('MapInfo File').DeleteDataSource('tmp/testogr2ogr23.mif')

    return 'success'

###############################################################################
# Test -clipsrc with WKT geometry (#3530)

def test_ogr2ogr_24():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'
        
    if not ogrtest.have_geos():
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/poly.shp ../ogr/data/poly.shp -clipsrc "POLYGON((479609 4764629,479609 4764817,479764 4764817,479764 4764629,479609 4764629))"')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 4:
        return 'fail'
        
    if ds.GetLayer(0).GetExtent() != (479609, 479764, 4764629, 4764817):
        print(ds.GetLayer(0).GetExtent())
        gdaltest.post_reason('unexpected extent')
        return 'fail'
        
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -clipsrc with clip from external datasource

def test_ogr2ogr_25():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'
        
    if not ogrtest.have_geos():
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass
    
    f = open('tmp/clip.csv', 'wt')
    f.write('foo,WKT\n')
    f.write('foo,"POLYGON((479609 4764629,479609 4764817,479764 4764817,479764 4764629,479609 4764629))"\n')
    f.close()

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/poly.shp ../ogr/data/poly.shp -clipsrc tmp/clip.csv -clipsrcwhere foo=\'foo\'')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 4:
        return 'fail'
        
    if ds.GetLayer(0).GetExtent() != (479609, 479764, 4764629, 4764817):
        print(ds.GetLayer(0).GetExtent())
        gdaltest.post_reason('unexpected extent')
        return 'fail'
        
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    os.remove('tmp/clip.csv')

    return 'success'

###############################################################################
# Test -clipdst with WKT geometry (#3530)

def test_ogr2ogr_26():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'
        
    if not ogrtest.have_geos():
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/poly.shp ../ogr/data/poly.shp -clipdst "POLYGON((479609 4764629,479609 4764817,479764 4764817,479764 4764629,479609 4764629))"')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 4:
        return 'fail'
        
    if ds.GetLayer(0).GetExtent() != (479609, 479764, 4764629, 4764817):
        print(ds.GetLayer(0).GetExtent())
        gdaltest.post_reason('unexpected extent')
        return 'fail'
        
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')

    return 'success'

###############################################################################
# Test -clipdst with clip from external datasource

def test_ogr2ogr_27():
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'
        
    if not ogrtest.have_geos():
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass
    
    f = open('tmp/clip.csv', 'wt')
    f.write('foo,WKT\n')
    f.write('foo,"POLYGON((479609 4764629,479609 4764817,479764 4764817,479764 4764629,479609 4764629))"\n')
    f.close()

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -nlt MULTIPOLYGON tmp/poly.shp ../ogr/data/poly.shp -clipdst tmp/clip.csv -clipdstsql "SELECT * from clip"')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 4:
        return 'fail'
        
    if ds.GetLayer(0).GetExtent() != (479609, 479764, 4764629, 4764817):
        print(ds.GetLayer(0).GetExtent())
        gdaltest.post_reason('unexpected extent')
        return 'fail'
        
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    os.remove('tmp/clip.csv')

    return 'success'


###############################################################################
# Test -wrapdateline on linestrings

def test_ogr2ogr_28():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/wrapdateline_src.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/wrapdateline_src.shp')
    except:
        pass
        
    try:
        os.stat('tmp/wrapdateline_dst.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/wrapdateline_dst.shp')
    except:
        pass
        
    ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('tmp/wrapdateline_src.shp')
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326);
    lyr = ds.CreateLayer('wrapdateline_src', srs = srs)
    feat = ogr.Feature(lyr.GetLayerDefn())
    geom = ogr.CreateGeometryFromWkt('LINESTRING(160 0,165 1,170 2,175 3,177 4,-177 5,-175 6,-170 7,-177 8,177 9,170 10)')
    feat.SetGeometryDirectly(geom)
    lyr.CreateFeature(feat)
    feat.Destroy()
    ds.Destroy()
    
    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -wrapdateline tmp/wrapdateline_dst.shp tmp/wrapdateline_src.shp')
    
    expected_wkt = 'MULTILINESTRING ((160 0,165 1,170 2,175 3,177 4,180 4.5),(-180 4.5,-177 5,-175 6,-170 7,-177 8,-180 8.5),(180 8.5,177 9,170 10))'
    expected_geom = ogr.CreateGeometryFromWkt(expected_wkt)
    ds = ogr.Open('tmp/wrapdateline_dst.shp')
    lyr = ds.GetLayer(0)
    feat = lyr.GetNextFeature()
    ret = ogrtest.check_feature_geometry(feat, expected_geom)
    feat.Destroy()
    expected_geom.Destroy()
    ds.Destroy()
    
    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/wrapdateline_src.shp')
    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/wrapdateline_dst.shp')

    if ret == 0:
        return 'success'
    else:
        return 'fail'

###############################################################################
# Test -wrapdateline on polygons

def test_ogr2ogr_29():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    if ogrtest.have_geos() is 0:
        return 'skip'

    for i in range(2):
        try:
            os.stat('tmp/wrapdateline_src.shp')
            ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/wrapdateline_src.shp')
        except:
            pass
            
        try:
            os.stat('tmp/wrapdateline_dst.shp')
            ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/wrapdateline_dst.shp')
        except:
            pass
            
        ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('tmp/wrapdateline_src.shp')
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326);
        lyr = ds.CreateLayer('wrapdateline_src', srs = srs)
        feat = ogr.Feature(lyr.GetLayerDefn())
        
        if i == 0:
            geom = ogr.CreateGeometryFromWkt('POLYGON((179 40,179.5 40,-179.5 40,-179 40,-170 40,-165 40,-165 30,-170 30,-179 30,-179.5 30,179.5 30,179 30,179 40))')
        else:
            geom = ogr.CreateGeometryFromWkt('POLYGON((-165 30,-170 30,-179 30,-179.5 30,179.5 30,179 30,179 40,179.5 40,-179.5 40,-179 40,-170 40,-165 40,-165 30))')
        feat.SetGeometry(geom)
        lyr.CreateFeature(feat)
        feat.Destroy()
        ds.Destroy()
        
        gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -wrapdateline tmp/wrapdateline_dst.shp tmp/wrapdateline_src.shp')
        
        expected_wkt = 'MULTIPOLYGON (((180 30,179.5 30.0,179 30,179 40,179.5 40.0,180 40,180 30)),((-180 40,-179.5 40.0,-179 40,-170 40,-165 40,-165 30,-170 30,-179 30,-179.5 30.0,-180 30,-180 40)))'
        expected_geom = ogr.CreateGeometryFromWkt(expected_wkt)
        ds = ogr.Open('tmp/wrapdateline_dst.shp')
        lyr = ds.GetLayer(0)
        feat = lyr.GetNextFeature()
        ret = ogrtest.check_feature_geometry(feat, expected_geom)
        if ret != 0:
            print('src is : %s' % geom.ExportToWkt())
            print('got    : %s' % feat.GetGeometryRef().ExportToWkt())
    
        feat.Destroy()
        expected_geom.Destroy()
        ds.Destroy()
        
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/wrapdateline_src.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/wrapdateline_dst.shp')
    
        if ret != 0:
            return 'fail'

    return 'success'

###############################################################################
# Test -splitlistfields option

def test_ogr2ogr_30():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        ds = ogr.Open('../ogr/data/testlistfields.gml')
    except:
        ds = None

    if ds is None:
        return 'skip'
    ds = None

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -splitlistfields tmp/test_ogr2ogr_30.shp ../ogr/data/testlistfields.gml')

    ds = ogr.Open('tmp/test_ogr2ogr_30.shp')
    if ds is None:
        return 'fail'
    lyr = ds.GetLayer(0)
    feat = lyr.GetNextFeature()

    if feat.GetField('attrib11') != 'value1' or \
       feat.GetField('attrib12') != 'value2' or \
       feat.GetField('attrib2') != 'value3' or \
       feat.GetField('attrib31') != 4 or \
       feat.GetField('attrib32') != 5 or \
       feat.GetField('attrib41') != 6.1 or \
       feat.GetField('attrib42') != 7.1:
        gdaltest.post_reason('did not get expected attribs')
        feat.DumpReadable()
        return 'fail'

    ds = None
    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_30.shp')

    return 'success'

###############################################################################
# Test that -overwrite work if the output file doesn't yet exist (#3825)

def test_ogr2ogr_31():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/poly.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -overwrite tmp/poly.shp ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        return 'fail'
    ds.Destroy()

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/poly.shp')
    return 'success'

###############################################################################
# Test that -append/-overwrite to a single-file shapefile work without specifying -nln

def test_ogr2ogr_32():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/test_ogr2ogr_32.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_32.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/test_ogr2ogr_32.shp ../ogr/data/poly.shp')
    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -append tmp/test_ogr2ogr_32.shp ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/test_ogr2ogr_32.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 20:
        gdaltest.post_reason('-append failed')
        return 'fail'
    ds = None

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -overwrite tmp/test_ogr2ogr_32.shp ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/test_ogr2ogr_32.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        gdaltest.post_reason('-overwrite failed')
        return 'fail'
    ds = None

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_32.shp')
    return 'success'

###############################################################################
# Test -explodecollections

def test_ogr2ogr_33():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/test_ogr2ogr_33_src.csv')
        ogr.GetDriverByName('CSV').DeleteDataSource('tmp/test_ogr2ogr_33_src.csv')
    except:
        pass

    try:
        os.stat('tmp/test_ogr2ogr_33_dst.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_33_dst.shp')
    except:
        pass

    f = open('tmp/test_ogr2ogr_33_src.csv', 'wt')
    f.write('foo,WKT\n')
    f.write('bar,"MULTIPOLYGON (((10 10,10 11,11 11,11 10,10 10)),((100 100,100 200,200 200,200 100,100 100),(125 125,175 125,175 175,125 175,125 125)))"\n')
    f.write('baz,"POLYGON ((0 0,0 1,1 1,1 0,0 0))"\n')
    f.close()

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -explodecollections tmp/test_ogr2ogr_33_dst.shp tmp/test_ogr2ogr_33_src.csv -select foo')

    ds = ogr.Open('tmp/test_ogr2ogr_33_dst.shp')
    lyr = ds.GetLayer(0)
    if lyr.GetFeatureCount() != 3:
        gdaltest.post_reason('-explodecollections failed')
        print(lyr.GetFeatureCount())
        return 'fail'

    feat = lyr.GetFeature(0)
    if feat.GetField("foo") != 'bar':
        feat.DumpReadable()
        return 'fail'
    if feat.GetGeometryRef().ExportToWkt() != 'POLYGON ((10 10,10 11,11 11,11 10,10 10))':
        feat.DumpReadable()
        return 'fail'

    feat = lyr.GetFeature(1)
    if feat.GetField("foo") != 'bar':
        feat.DumpReadable()
        return 'fail'
    if feat.GetGeometryRef().ExportToWkt() != 'POLYGON ((100 100,100 200,200 200,200 100,100 100),(125 125,175 125,175 175,125 175,125 125))':
        feat.DumpReadable()
        return 'fail'

    feat = lyr.GetFeature(2)
    if feat.GetField("foo") != 'baz':
        feat.DumpReadable()
        return 'fail'
    if feat.GetGeometryRef().ExportToWkt() != 'POLYGON ((0 0,0 1,1 1,1 0,0 0))':
        feat.DumpReadable()
        return 'fail'

    ds = None

    ogr.GetDriverByName('CSV').DeleteDataSource('tmp/test_ogr2ogr_33_src.csv')
    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_33_dst.shp')
    return 'success'

###############################################################################
# Test 'ogr2ogr someDirThatDoesNotExist src.shp -nln someDirThatDoesNotExist'
# This should result in creating a someDirThatDoesNotExist directory with
# someDirThatDoesNotExist.shp/dbf/shx inside this directory

def test_ogr2ogr_34():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/test_ogr2ogr_34_dir')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_34_dir')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/test_ogr2ogr_34_dir ../ogr/data/poly.shp -nln test_ogr2ogr_34_dir')

    ds = ogr.Open('tmp/test_ogr2ogr_34_dir/test_ogr2ogr_34_dir.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        gdaltest.post_reason('initial shapefile creation failed')
        return 'fail'
    ds = None

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -append tmp/test_ogr2ogr_34_dir ../ogr/data/poly.shp -nln test_ogr2ogr_34_dir')

    ds = ogr.Open('tmp/test_ogr2ogr_34_dir/test_ogr2ogr_34_dir.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 20:
        gdaltest.post_reason('-append failed')
        return 'fail'
    ds = None

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -overwrite tmp/test_ogr2ogr_34_dir ../ogr/data/poly.shp -nln test_ogr2ogr_34_dir')

    ds = ogr.Open('tmp/test_ogr2ogr_34_dir/test_ogr2ogr_34_dir.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        gdaltest.post_reason('-overwrite failed')
        return 'fail'
    ds = None

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_34_dir')
    return 'success'

###############################################################################
# Test 'ogr2ogr someDirThatDoesNotExist src.shp'

def test_ogr2ogr_35():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/test_ogr2ogr_35_dir')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_35_dir')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/test_ogr2ogr_35_dir ../ogr/data/poly.shp ')

    ds = ogr.Open('tmp/test_ogr2ogr_35_dir/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        gdaltest.post_reason('initial shapefile creation failed')
        return 'fail'
    ds = None

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -append tmp/test_ogr2ogr_35_dir ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/test_ogr2ogr_35_dir/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 20:
        gdaltest.post_reason('-append failed')
        return 'fail'
    ds = None

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -overwrite tmp/test_ogr2ogr_35_dir ../ogr/data/poly.shp')

    ds = ogr.Open('tmp/test_ogr2ogr_35_dir/poly.shp')
    if ds is None or ds.GetLayer(0).GetFeatureCount() != 10:
        gdaltest.post_reason('-overwrite failed')
        return 'fail'
    ds = None

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_35_dir')
    return 'success'

###############################################################################
# Test ogr2ogr -zfield

def test_ogr2ogr_36():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/test_ogr2ogr_36.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_36.shp')
    except:
        pass

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/test_ogr2ogr_36.shp ../ogr/data/poly.shp -zfield EAS_ID')

    ds = ogr.Open('tmp/test_ogr2ogr_36.shp')
    feat = ds.GetLayer(0).GetNextFeature()
    wkt = feat.GetGeometryRef().ExportToWkt()
    ds = None

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_36.shp')

    if wkt.find(' 168,') == -1:
        return 'fail'

    return 'success'

###############################################################################
# Test 'ogr2ogr someDirThatDoesNotExist.shp dataSourceWithMultipleLayer'

def test_ogr2ogr_37():

    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    try:
        os.stat('tmp/test_ogr2ogr_37_dir.shp')
        ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_37_dir.shp')
    except:
        pass

    try:
        os.mkdir('tmp/test_ogr2ogr_37_src')
    except:
        pass
    shutil.copy('../ogr/data/poly.shp', 'tmp/test_ogr2ogr_37_src')
    shutil.copy('../ogr/data/poly.shx', 'tmp/test_ogr2ogr_37_src')
    shutil.copy('../ogr/data/poly.dbf', 'tmp/test_ogr2ogr_37_src')
    shutil.copy('../ogr/data/testpoly.shp', 'tmp/test_ogr2ogr_37_src')
    shutil.copy('../ogr/data/testpoly.shx', 'tmp/test_ogr2ogr_37_src')
    shutil.copy('../ogr/data/testpoly.dbf', 'tmp/test_ogr2ogr_37_src')

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' tmp/test_ogr2ogr_37_dir.shp tmp/test_ogr2ogr_37_src')

    ds = ogr.Open('tmp/test_ogr2ogr_37_dir.shp')
    if ds is None or ds.GetLayerCount() != 2:
        return 'fail'
    ds = None

    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_37_src')
    ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('tmp/test_ogr2ogr_37_dir.shp')

    return 'success'

gdaltest_list = [
    test_ogr2ogr_1,
    test_ogr2ogr_2,
    test_ogr2ogr_3,
    test_ogr2ogr_4,
    test_ogr2ogr_5,
    test_ogr2ogr_6,
    test_ogr2ogr_7,
    test_ogr2ogr_8,
    test_ogr2ogr_9,
    test_ogr2ogr_10,
    test_ogr2ogr_11,
    test_ogr2ogr_12,
    test_ogr2ogr_13,
    test_ogr2ogr_14,
    test_ogr2ogr_15,
    test_ogr2ogr_16,
    test_ogr2ogr_17,
    test_ogr2ogr_18,
    test_ogr2ogr_19,
    test_ogr2ogr_20,
    test_ogr2ogr_21,
    test_ogr2ogr_22,
    test_ogr2ogr_23,
    test_ogr2ogr_24,
    test_ogr2ogr_25,
    test_ogr2ogr_26,
    test_ogr2ogr_27,
    test_ogr2ogr_28,
    test_ogr2ogr_29,
    test_ogr2ogr_30,
    test_ogr2ogr_31,
    test_ogr2ogr_32,
    test_ogr2ogr_33,
    test_ogr2ogr_34,
    test_ogr2ogr_35,
    test_ogr2ogr_36,
    test_ogr2ogr_37 ]
    
if __name__ == '__main__':

    gdaltest.setup_run( 'test_ogr2ogr' )

    gdaltest.run_tests( gdaltest_list )

    gdaltest.summarize()

