#!/usr/bin/env python
###############################################################################
# $Id$
#
# Project:  GDAL/OGR Test Suite
# Purpose:  WFS driver testing.
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
import string

sys.path.append( '../pymod' )

import gdaltest
import ogrtest
import ogr
import osr
import gdal
import webserver

###############################################################################
# Test underlying OGR drivers
#

def ogr_wfs_init():

    gdaltest.wfs_drv = None

    try:
        gdaltest.wfs_drv = ogr.GetDriverByName('WFS')
    except:
        pass
        
    if gdaltest.wfs_drv is None:
        return 'skip'
    
    gdaltest.have_gml_reader = 0
    gdaltest.geoserver_wfs = None
    gdaltest.deegree_wfs = None
    gdaltest.ionic_wfs = None

    try:
        gml_ds = ogr.Open( 'data/ionic_wfs.gml' )
    except:
        gml_ds = None

    if gml_ds is None:
        if gdal.GetLastErrorMsg().find('Xerces') != -1:
            return 'skip'
        else:
            gdaltest.post_reason( 'failed to open test file.' )
            return 'skip'

    gdaltest.have_gml_reader = 1

    return 'success'

###############################################################################
# Test reading a MapServer WFS server

def ogr_wfs_mapserver():
    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        return 'skip'

    if gdaltest.gdalurlopen('http://www2.dmsolutions.ca/cgi-bin/mswfs_gmap') is None:
        print('cannot open URL')
        return 'skip'

    ds = ogr.Open('WFS:http://www2.dmsolutions.ca/cgi-bin/mswfs_gmap')
    if ds is None:
        gdaltest.post_reason('did not managed to open WFS datastore')
        return 'skip'

    if ds.GetLayerCount() != 2:
        gdaltest.post_reason('did not get expected layer count')
        print(ds.GetLayerCount())
        return 'fail'

    lyr = ds.GetLayer(0)
    if lyr.GetName() != 'park':
        gdaltest.post_reason('did not get expected layer name')
        print(lyr.GetName())
        return 'fail'

    sr = lyr.GetSpatialRef()
    sr2 = osr.SpatialReference()
    sr2.ImportFromEPSG(42304)
    if not sr.IsSame(sr2):
        gdaltest.post_reason('did not get expected SRS')
        print(sr)
        return 'fail'

    feat_count = lyr.GetFeatureCount()
    if feat_count != 46:
        gdaltest.post_reason('did not get expected feature count')
        print(feat_count)
        return 'fail'

    feat = lyr.GetNextFeature()
    geom = feat.GetGeometryRef()
    geom_wkt = geom.ExportToWkt()
    if geom_wkt.find("POLYGON ((389366.84375 3791519.75") == -1:
        gdaltest.post_reason('did not get expected feature')
        feat.DumpReadable()
        return 'fail'

    return 'success'


###############################################################################
# Test reading a GeoServer WFS server

def ogr_wfs_geoserver():
    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        return 'skip'

    if gdaltest.gdalurlopen('http://demo.opengeo.org/geoserver/wfs') is None:
        print('cannot open URL')
        gdaltest.geoserver_wfs = False
        return 'skip'
    gdaltest.geoserver_wfs = True

    ds = ogr.Open('WFS:http://demo.opengeo.org/geoserver/wfs?TYPENAME=za:za_points')
    if ds is None:
        gdaltest.post_reason('did not managed to open WFS datastore')
        return 'fail'

    if ds.GetLayerCount() != 1:
        gdaltest.post_reason('did not get expected layer count')
        print(ds.GetLayerCount())
        return 'fail'

    lyr = ds.GetLayer(0)
    if lyr.GetName() != 'za:za_points':
        gdaltest.post_reason('did not get expected layer name')
        print(lyr.GetName())
        return 'fail'

    sr = lyr.GetSpatialRef()
    sr2 = osr.SpatialReference()
    sr2.ImportFromEPSG(4326)
    if not sr.IsSame(sr2):
        gdaltest.post_reason('did not get expected SRS')
        print(sr)
        return 'fail'

    feat_count = lyr.GetFeatureCount()
    if feat_count < 14000:
        gdaltest.post_reason('did not get expected feature count')
        print(feat_count)
        return 'fail'

    if not lyr.TestCapability(ogr.OLCFastFeatureCount):
        gdaltest.post_reason('did not get OLCFastFeatureCount')
        return 'fail'

    ds = ogr.Open('WFS:http://demo.opengeo.org/geoserver/wfs?TYPENAME=za:za_points&MAXFEATURES=10')
    if ds is None:
        print('server perhaps overloaded')
        return 'skip'
    lyr = ds.GetLayer(0)
    gdal.ErrorReset()
    feat = lyr.GetNextFeature()

    # This error message is generally the sign of a server in a broken state
    if feat is None and gdal.GetLastErrorMsg().find('<ows:ExceptionText>org.geoserver.platform.ServiceException') != -1:
        print('server probably in a broken state')
        # Disable it for wfs-t test
        gdaltest.geoserver_wfs = False
        return 'skip'
    
    geom = feat.GetGeometryRef()
    geom_wkt = geom.ExportToWkt()
    if feat.GetField('name') != 'Alexander Bay' or \
       geom_wkt.find("POINT (16.4827778 -28.5947222)") == -1:
        gdaltest.post_reason('did not get expected feature (1)')
        feat.DumpReadable()
        return 'fail'

    # Same with VERSION=1.0.0
    ds = ogr.Open('WFS:http://demo.opengeo.org/geoserver/wfs?TYPENAME=za:za_points&MAXFEATURES=10&VERSION=1.0.0')
    if ds is None:
        print('server perhaps overloaded')
        return 'skip'
    lyr = ds.GetLayer(0)
    feat = lyr.GetNextFeature()
    geom = feat.GetGeometryRef()
    geom_wkt = geom.ExportToWkt()
    if feat.GetField('name') != 'Alexander Bay' or \
       ogrtest.check_feature_geometry(feat,'POINT (16.4827778 -28.5947222)',
                                      max_error = 0.000000001 ) != 0:
        gdaltest.post_reason('did not get expected feature (2)')
        feat.DumpReadable()
        return 'fail'

    # Test attribute filter
    ds = ogr.Open("WFS:http://demo.opengeo.org/geoserver/wfs?TYPENAME=za:za_points")
    if ds is None:
        print('server perhaps overloaded')
        return 'skip'
    lyr = ds.GetLayer(0)
    lyr.SetAttributeFilter("type is not null and name >= 'W' and type LIKE 'att%%ion'")
    feat_count = lyr.GetFeatureCount()
    if feat_count != 1:
        gdaltest.post_reason('did not get expected feature count after SetAttributeFilter (1)')
        print(feat_count)
        return 'fail'
    feat = lyr.GetNextFeature()
    if feat.GetField('gml_id') != 'za_points.2839':
        gdaltest.post_reason('did not get expected feature (3)')
        feat.DumpReadable()
        return 'fail'

    if False:
        # This GeoServer version doesn't understand <GmlObjectId>
        lyr.SetAttributeFilter("gml_id = 'za_points.2839'")
        feat_count = lyr.GetFeatureCount()
        if feat_count != 1:
            gdaltest.post_reason('did not get expected feature count after SetAttributeFilter (2)')
            print(feat_count)
            return 'fail'
        feat = lyr.GetNextFeature()
        if feat.GetField('gml_id') != 'za_points.2839':
            gdaltest.post_reason('did not get expected feature (4)')
            feat.DumpReadable()
            return 'fail'

    return 'success'

###############################################################################
# Test reading a GeoServer WFS server with OUTPUTFORMAT=json

def ogr_wfs_geoserver_json():
    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        if gdaltest.geoserver_wfs is None:
            if gdaltest.gdalurlopen('http://demo.opengeo.org/geoserver/wfs') is None:
                gdaltest.geoserver_wfs = False
                print('cannot open URL')
                return 'skip'
            gdaltest.geoserver_wfs = True

    if gdaltest.geoserver_wfs != True:
        return 'skip'

    ds = ogr.Open('WFS:http://demo.opengeo.org/geoserver/wfs?TYPENAME=za:za_points&MAXFEATURES=10&OUTPUTFORMAT=json')
    if ds is None:
        gdaltest.post_reason('did not managed to open WFS datastore')
        return 'fail'

    if ds.GetLayerCount() != 1:
        gdaltest.post_reason('did not get expected layer count')
        print(ds.GetLayerCount())
        return 'fail'

    lyr = ds.GetLayer(0)
    if lyr.GetName() != 'za:za_points':
        gdaltest.post_reason('did not get expected layer name')
        print(lyr.GetName())
        return 'fail'

    feat_count = lyr.GetFeatureCount()
    if feat_count != 10:
        gdaltest.post_reason('did not get expected feature count')
        print(feat_count)
        return 'fail'

    if not lyr.TestCapability(ogr.OLCFastFeatureCount):
        gdaltest.post_reason('did not get OLCFastFeatureCount')
        return 'fail'

    feat = lyr.GetNextFeature()
    geom = feat.GetGeometryRef()
    geom_wkt = geom.ExportToWkt()
    if feat.GetField('name') != 'Alexander Bay' or \
       ogrtest.check_feature_geometry(feat,'POINT (16.4827778 -28.5947222)',
                                      max_error = 0.000000001 ) != 0:
        gdaltest.post_reason('did not get expected feature')
        feat.DumpReadable()
        return 'fail'

    return 'success'


###############################################################################
# Test reading a GeoServer WFS server with OUTPUTFORMAT=SHAPE-ZIP

def ogr_wfs_geoserver_shapezip():
    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        if gdaltest.geoserver_wfs is None:
            if gdaltest.gdalurlopen('http://demo.opengeo.org/geoserver/wfs') is None:
                gdaltest.geoserver_wfs = False
                print('cannot open URL')
                return 'skip'
            gdaltest.geoserver_wfs = True

    if gdaltest.geoserver_wfs != True:
        return 'skip'

    ds = ogr.Open('WFS:http://demo.opengeo.org/geoserver/wfs?TYPENAME=za:za_points&MAXFEATURES=10&OUTPUTFORMAT=SHAPE-ZIP')
    if ds is None:
        gdaltest.post_reason('did not managed to open WFS datastore')
        return 'fail'

    if ds.GetLayerCount() != 1:
        gdaltest.post_reason('did not get expected layer count')
        print(ds.GetLayerCount())
        return 'fail'

    lyr = ds.GetLayer(0)
    if lyr.GetName() != 'za:za_points':
        gdaltest.post_reason('did not get expected layer name')
        print(lyr.GetName())
        return 'fail'

    feat_count = lyr.GetFeatureCount()
    if feat_count != 10:
        gdaltest.post_reason('did not get expected feature count')
        print(feat_count)
        return 'fail'

    if not lyr.TestCapability(ogr.OLCFastFeatureCount):
        gdaltest.post_reason('did not get OLCFastFeatureCount')
        return 'fail'

    feat = lyr.GetNextFeature()
    geom = feat.GetGeometryRef()
    geom_wkt = geom.ExportToWkt()
    if feat.GetField('name') != 'Alexander Bay' or \
       ogrtest.check_feature_geometry(feat,'POINT (16.4827778 -28.5947222)',
                                      max_error = 0.000000001 ) != 0:
        gdaltest.post_reason('did not get expected feature')
        feat.DumpReadable()
        return 'fail'

    return 'success'

###############################################################################
# Test reading a Deegree WFS server

def ogr_wfs_deegree():
    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        return 'skip'

    if gdaltest.gdalurlopen('http://demo.deegree.org/deegree-wfs/services') is None:
        gdaltest.deegree_wfs = False
        print('cannot open URL')
        return 'skip'
    gdaltest.deegree_wfs = True

    ds = ogr.Open("WFS:http://demo.deegree.org/deegree-wfs/services?MAXFEATURES=10")
    if ds is None:
        gdaltest.post_reason('did not managed to open WFS datastore')
        return 'fail'

    lyr = ds.GetLayerByName('app:Springs')
    if lyr.GetName() != 'app:Springs':
        gdaltest.post_reason('did not get expected layer name')
        print(lyr.GetName())
        return 'fail'

    sr = lyr.GetSpatialRef()
    sr2 = osr.SpatialReference()
    sr2.ImportFromEPSG(26912)
    if not sr.IsSame(sr2):
        gdaltest.post_reason('did not get expected SRS')
        print(sr)
        return 'fail'

    feat = lyr.GetNextFeature()
    geom = feat.GetGeometryRef()
    geom_wkt = geom.ExportToWkt()
    if feat.GetField('objectid') != 1 or \
       ogrtest.check_feature_geometry(feat,'POINT (558750.703125 4402882.05)',
                                      max_error = 0.000000001 ) != 0:
        gdaltest.post_reason('did not get expected feature')
        feat.DumpReadable()
        return 'fail'

    # Test attribute filter
    ds = ogr.Open("WFS:http://demo.deegree.org/deegree-wfs/services")
    lyr = ds.GetLayerByName('app:Springs')
    lyr.SetAttributeFilter('objectid = 9 or objectid = 100 or (objectid >= 20 and objectid <= 30 and objectid != 27)')
    feat_count = lyr.GetFeatureCount()
    if feat_count != 12:
        gdaltest.post_reason('did not get expected feature count after SetAttributeFilter')
        print(feat_count)
        return 'fail'

    # Test attribute filter with gml_id
    lyr.SetAttributeFilter("gml_id = 'SGID024_Springs30' or gml_id = 'SGID024_Springs100'")
    feat_count = lyr.GetFeatureCount()
    if feat_count != 2:
        gdaltest.post_reason('did not get expected feature count after SetAttributeFilter (2)')
        print(feat_count)
        return 'fail'
    return 'success'

###############################################################################
# Run test_ogrsf

def ogr_wfs_test_ogrsf():
    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        return 'skip'

    if gdaltest.deegree_wfs != True:
        return 'skip'

    import test_cli_utilities
    if test_cli_utilities.get_test_ogrsf_path() is None:
        return 'skip'

    ret = gdaltest.runexternal(test_cli_utilities.get_test_ogrsf_path() + ' -ro "WFS:http://demo.deegree.org/deegree-wfs/services?MAXFEATURES=10" app:Springs')

    if ret.find('INFO') == -1 or ret.find('ERROR') != -1:
        print(ret)
        return 'fail'

    return 'success'

###############################################################################
# Test reading a local fake WFS server

def ogr_wfs_fake_wfs_server():
    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        return 'skip'

    (process, port) = webserver.launch()
    if port == 0:
        return 'skip'

    ds = ogr.Open("WFS:http://127.0.0.1:%d/fakewfs" % port)
    if ds is None:
        gdaltest.post_reason('did not managed to open WFS datastore')
        webserver.server_stop(process, port)
        return 'fail'

    lyr = ds.GetLayerByName('rijkswegen')
    if lyr.GetName() != 'rijkswegen':
        gdaltest.post_reason('did not get expected layer name')
        print(lyr.GetName())
        webserver.server_stop(process, port)
        return 'fail'

    sr = lyr.GetSpatialRef()
    sr2 = osr.SpatialReference()
    sr2.ImportFromEPSG(28992)
    if not sr.IsSame(sr2):
        gdaltest.post_reason('did not get expected SRS')
        print(sr)
        webserver.server_stop(process, port)
        return 'fail'

    feat = lyr.GetNextFeature()
    if feat.GetField('MPLength') != '33513.' or \
       ogrtest.check_feature_geometry(feat,'MULTILINESTRING ((154898.65286 568054.62753,160108.36082 566076.78094,164239.254332 563024.70188,170523.31535 561231.219583,172676.42256 559253.37299,175912.80562 557459.89069,180043.699132 553508.779495,183294.491306 552250.182732))',
                                      max_error = 0.00001 ) != 0:
        gdaltest.post_reason('did not get expected feature')
        feat.DumpReadable()
        webserver.server_stop(process, port)
        return 'fail'

    webserver.server_stop(process, port)

    return 'success'

###############################################################################
# Test CreateFeature() / UpdateFeature() / DeleteFeature() (WFS-T)

def ogr_wfs_geoserver_wfst():
    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        return 'skip'

    if gdaltest.geoserver_wfs != True:
        return 'skip'

    ds = ogr.Open('WFS:http://demo.opengeo.org/geoserver/wfs', update = 1)
    if ds is None:
        return 'fail'

    lyr = ds.GetLayerByName('za:za_points')
    geom = ogr.CreateGeometryFromWkt('POINT(0 89.5)')
    feat = ogr.Feature(lyr.GetLayerDefn())
    feat.SetGeometry(geom)
    feat.SetField('name', 'name_set_by_ogr_wfs_8_test')
    feat.SetField('type', 'type_set_by_ogr_wfs_8_test')
    if lyr.CreateFeature(feat) != 0:
        gdaltest.post_reason('cannot create feature')
        return 'fail'

    print('Feature %d created !' % feat.GetFID())

    feat.SetField('name', 'name_modified_by_ogr_wfs_8_test')
    if lyr.SetFeature(feat) != 0:
        gdaltest.post_reason('cannot update feature')
        return 'fail'
    print('Feature %d updated !' % feat.GetFID())
    
    if lyr.DeleteFeature(feat.GetFID()) != 0:
        gdaltest.post_reason('could not delete feature')
        return 'fail'

    print('Feature %d deleted !' % feat.GetFID())

    # Test transactions
    if lyr.StartTransaction() != 0:
        gdaltest.post_reason('CommitTransaction() failed')
        return 'fail'
        
    geom = ogr.CreateGeometryFromWkt('POINT(0 89.5)')
    feat = ogr.Feature(lyr.GetLayerDefn())
    feat.SetGeometry(geom)
    feat.SetField('name', 'name_set_by_ogr_wfs_8_test')
    feat.SetField('type', 'type_set_by_ogr_wfs_8_test')
    if lyr.CreateFeature(feat) != 0:
        gdaltest.post_reason('cannot create feature')
        return 'fail'
    geom = ogr.CreateGeometryFromWkt('POINT(0 89.5)')
    feat = ogr.Feature(lyr.GetLayerDefn())
    feat.SetGeometry(geom)
    feat.SetField('name', 'name_set_by_ogr_wfs_8_test_2')
    feat.SetField('type', 'type_set_by_ogr_wfs_8_test_2')
    if lyr.CreateFeature(feat) != 0:
        gdaltest.post_reason('cannot create feature')
        return 'fail'

    if lyr.CommitTransaction() != 0:
        gdaltest.post_reason('CommitTransaction() failed')
        return 'fail'

    # Retrieve inserted features
    print('Retrieving created features gml:id')
    sql_lyr = ds.ExecuteSQL("SELECT _LAST_INSERTED_FIDS_ FROM za:za_points");
    feat = sql_lyr.GetNextFeature()
    while feat is not None:
        gml_id = feat.GetFieldAsString(0)
        print('Feature %s has been created in transaction !' % gml_id)
        feat = sql_lyr.GetNextFeature()
    feat = None
    count = sql_lyr.GetFeatureCount()
    ds.ReleaseResultSet(sql_lyr)

    if count != 2:
        gdaltest.post_reason('did not get expected feature count')
        return 'fail'

    # Delete a bunch of features
    print('Deleting created features')
    sql_lyr = ds.ExecuteSQL("DELETE FROM za:za_points WHERE name = 'name_set_by_ogr_wfs_8_test' OR name = 'name_set_by_ogr_wfs_8_test_2'")
    ds.ReleaseResultSet(sql_lyr)

    return 'success'


###############################################################################
# Test CreateFeature() / UpdateFeature() / DeleteFeature() with expected
# failure due to server not allowing insert & delete

def ogr_wfs_deegree_wfst():

    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        return 'skip'

    if gdaltest.gdalurlopen('http://testing.deegree.org/deegree-wfs/services') is None:
        print('cannot open URL')
        return 'skip'

    ds = ogr.Open('WFS:http://testing.deegree.org/deegree-wfs/services', update = 1)
    if ds is None:
        return 'fail'

    lyr = ds.GetLayerByName('app:CountyBoundaries_edited')
    geom = ogr.CreateGeometryFromWkt('POINT(2 49)')
    feat = ogr.Feature(lyr.GetLayerDefn())
    feat.SetGeometry(geom)
    feat.SetField('name', 'nameSetByOGR')
    feat.SetField('fips', '10')
    feat.SetField('feature_id', '123456')
    feat.SetField('objectid', '7890123')
    feat.SetField('shape_area', 12.34)
    feat.SetField('shape_len', 56.78)

    ret = lyr.CreateFeature(feat)
    if ret != 0:
        print('expected fail on CreateFeature')

    ret = lyr.DeleteFeature(1)
    if ret != 0:
        print('expected fail on DeleteFeature')

    feat = lyr.GetFeature(10)
    ret = lyr.SetFeature(feat)
    if ret != 0:
        print('expected fail on SetFeature')

    return 'success'

###############################################################################
# Test CreateFeature() / UpdateFeature() / DeleteFeature() on a WFS 1.0.0 server

def ogr_wfs_ionic_wfst():

    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        return 'skip'

    if gdaltest.gdalurlopen('http://webservices.ionicsoft.com/ionicweb/wfs/BOSTON_ORA') is None:
        print('cannot open URL')
        gdaltest.ionic_wfs = False
        return 'skip'
    gdaltest.ionic_wfs = True

    ds = ogr.Open('WFS:http://webservices.ionicsoft.com/ionicweb/wfs/BOSTON_ORA', update = 1)
    if ds is None:
        if gdal.GetLastErrorMsg().find('HTTP error code : 403') != -1:
            gdaltest.ionic_wfs = False
            return 'skip'
        return 'fail'

    lyr = ds.GetLayerByName('wfs:BUSINESS')
    geom = ogr.CreateGeometryFromWkt('POINT(234000 890000)')
    feat = ogr.Feature(lyr.GetLayerDefn())
    feat.SetGeometry(geom)
    feat.SetField('NAME', 'nameSetByOGR')
    feat.SetField('TOTAL_EMPLOYEES', '10')

    ret = lyr.CreateFeature(feat)
    if ret != 0:
        print('fail on CreateFeature')
        return 'fail'

    gmlid = feat.GetField('gml_id')

    ret = lyr.SetFeature(feat)
    if ret != 0:
        print('fail on SetFeature')
        return 'fail'

    ds.ExecuteSQL("DELETE FROM wfs:BUSINESS WHERE gml_id = '%s'" % gmlid)

    return 'success'

###############################################################################
# Test ExecuteSQL() where SQL should be turned into PROPERTYNAME and FILTER parameters

def ogr_wfs_ionic_sql():

    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        return 'skip'

    if gdaltest.ionic_wfs != True:
        return 'skip'

    ds = ogr.Open('WFS:http://webservices.ionicsoft.com/ionicweb/wfs/BOSTON_ORA')
    if ds is None:
        return 'fail'

    lyr = ds.ExecuteSQL("SELECT name FROM 'wfs:BUSINESS' WHERE total_employees = 105")
    count = lyr.GetFeatureCount()

    ds.ReleaseResultSet(lyr)

    if count != 1:
        return 'fail'

    return 'success'

###############################################################################
# Test opening a datasource from a XML description file
# The following test should issue 0 WFS http request

def ogr_wfs_xmldescriptionfile():

    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        return 'skip'

    ds = ogr.Open('data/testwfs.xml')
    lyr = ds.GetLayer(0)
    feature_defn = lyr.GetLayerDefn()
    index = feature_defn.GetFieldIndex('name')
    sr = lyr.GetSpatialRef()
    ds = None

    if index != 1:
        print(index)
        return 'fail'

    wkt = sr.ExportToWkt()
    if wkt.find('WGS 84') == -1:
        print(wkt)
        return 'fail'

    return 'success'

###############################################################################
# Test opening a datastore which only support GML 3.2.1 output

def ogr_wfs_deegree_gml321():

    if gdaltest.wfs_drv is None:
        return 'skip'
    if not gdaltest.have_gml_reader:
        return 'skip'

    if gdaltest.gdalurlopen('http://deegree3-testing.deegree.org:80/deegree-inspire-node/services') is None:
        print('cannot open URL')
        return 'skip'

    ds = ogr.Open('WFS:http://deegree3-testing.deegree.org:80/deegree-inspire-node/services?MAXFEATURES=10')
    if ds is None:
        return 'fail'

    lyr = ds.GetLayerByName("ad:Address")
    count = lyr.GetFeatureCount()
    if count != 10:
        print(count)
        return 'fail'

    lyr = ds.GetLayerByName("au:AdministrativeBoundary")
    count = lyr.GetFeatureCount()
    if count != 0:
        print(count)
        return 'fail'

    return 'success'

gdaltest_list = [ 
    ogr_wfs_init,
    ogr_wfs_mapserver,
    ogr_wfs_geoserver,
    ogr_wfs_geoserver_json,
    ogr_wfs_geoserver_shapezip,
    ogr_wfs_deegree,
    ogr_wfs_test_ogrsf,
    ogr_wfs_fake_wfs_server,
    ogr_wfs_geoserver_wfst,
    #ogr_wfs_deegree_wfst,
    ogr_wfs_ionic_wfst,
    ogr_wfs_ionic_sql,
    ogr_wfs_xmldescriptionfile,
    ogr_wfs_deegree_gml321,
    ]

if __name__ == '__main__':

    gdaltest.setup_run( 'ogr_wfs' )

    gdaltest.run_tests( gdaltest_list )

    gdaltest.summarize()

