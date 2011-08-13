#!/usr/bin/env python
###############################################################################
# $Id$
#
# Project:  GDAL/OGR Test Suite
# Purpose:  GML Reading Driver testing.
# Author:   Frank Warmerdam <warmerdam@pobox.com>
# 
###############################################################################
# Copyright (c) 2006, Frank Warmerdam <warmerdam@pobox.com>
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
import shutil

###############################################################################
# Test reading geometry and attribute from ionic wfs gml file.
#

def ogr_gml_1():

    gdaltest.have_gml_reader = 0

    try:
        gml_ds = ogr.Open( 'data/ionic_wfs.gml' )
    except:
        gml_ds = None

    if gml_ds is None:
        if gdal.GetLastErrorMsg().find('Xerces') != -1:
            return 'skip'
        else:
            gdaltest.post_reason( 'failed to open test file.' )
            return 'fail'

    gdaltest.have_gml_reader = 1

    if gml_ds.GetLayerCount() != 1:
        gdaltest.post_reason( 'wrong number of layers' )
        return 'fail'
    
    lyr = gml_ds.GetLayerByName('GEM')
    feat = lyr.GetNextFeature()

    if feat.GetField('Name') != 'Aartselaar':
        gdaltest.post_reason( 'Wrong name field value' )
        return 'fail'

    wkt = 'POLYGON ((44038 511549,44015 511548,43994 511522,43941 511539,43844 511514,43754 511479,43685 511521,43594 511505,43619 511452,43645 511417,4363 511387,437 511346,43749 511298,43808 511229,43819 511205,4379 511185,43728 511167,43617 511175,43604 511151,43655 511125,43746 511143,43886 511154,43885 511178,43928 511186,43977 511217,4404 511223,44008 511229,44099 51131,44095 511335,44106 51135,44127 511379,44124 511435,44137 511455,44105 511467,44098 511484,44086 511499,4407 511506,44067 511535,44038 511549))'
    
    if ogrtest.check_feature_geometry( feat, wkt):
        return 'fail'

    feat = lyr.GetNextFeature()
    if feat is not None:
        gdaltest.post_reason( 'got unexpected feature.' )
        return 'fail'

    return 'success'

###############################################################################
# Do the same test somewhere without a .gfs file.

def ogr_gml_2():
    if not gdaltest.have_gml_reader:
        return 'skip'

    # copy gml file (but not .gfs file)
    open('tmp/ionic_wfs.gml','w').write(open('data/ionic_wfs.gml').read())
    
    gml_ds = ogr.Open( 'tmp/ionic_wfs.gml' )    

    if gml_ds.GetLayerCount() != 1:
        gdaltest.post_reason( 'wrong number of layers' )
        return 'fail'
    
    lyr = gml_ds.GetLayerByName('GEM')
    feat = lyr.GetNextFeature()

    if feat.GetField('Name') != 'Aartselaar':
        gdaltest.post_reason( 'Wrong name field value' )
        return 'fail'

    wkt = 'POLYGON ((44038 511549,44015 511548,43994 511522,43941 511539,43844 511514,43754 511479,43685 511521,43594 511505,43619 511452,43645 511417,4363 511387,437 511346,43749 511298,43808 511229,43819 511205,4379 511185,43728 511167,43617 511175,43604 511151,43655 511125,43746 511143,43886 511154,43885 511178,43928 511186,43977 511217,4404 511223,44008 511229,44099 51131,44095 511335,44106 51135,44127 511379,44124 511435,44137 511455,44105 511467,44098 511484,44086 511499,4407 511506,44067 511535,44038 511549))'
    
    if ogrtest.check_feature_geometry( feat, wkt):
        return 'fail'

    feat = lyr.GetNextFeature()
    if feat is not None:
        gdaltest.post_reason( 'got unexpected feature.' )
        return 'fail'

    return 'success'

###############################################################################
# Similar test for RNF style line data.

def ogr_gml_3():
    if not gdaltest.have_gml_reader:
        return 'skip'

    gml_ds = ogr.Open( 'data/rnf_eg.gml' )    

    if gml_ds.GetLayerCount() != 1:
        gdaltest.post_reason( 'wrong number of layers' )
        return 'fail'
    
    lyr = gml_ds.GetLayerByName('RoadSegment')
    feat = lyr.GetNextFeature()

    if feat.GetField('ngd_id') != 817792:
        gdaltest.post_reason( 'Wrong ngd_id field value' )
        return 'fail'

    if feat.GetField('type') != 'HWY':
        gdaltest.post_reason( 'Wrong type field value' )
        return 'fail'

    wkt = 'LINESTRING (-63.500411040289066 46.240122507771368,-63.501009714909742 46.240344881690326,-63.502170462373471 46.241041855639622,-63.505862621395394 46.24195250605576,-63.506719184531178 46.242002742901576,-63.507197272602212 46.241931577811606,-63.508403092799554 46.241752283460158,-63.509946573455622 46.241745397977233)'
    
    if ogrtest.check_feature_geometry( feat, wkt):
        return 'fail'

    feat = lyr.GetNextFeature()
    if feat is not None:
        gdaltest.post_reason( 'got unexpected feature.' )
        return 'fail'

    return 'success'

###############################################################################
# Test of read GML file with UTF-8 BOM indicator.
# Test also support for nested GML elements (#3680)

def ogr_gml_4():
    if not gdaltest.have_gml_reader:
        return 'skip'

    gml_ds = ogr.Open( 'data/bom.gml' )    

    if gml_ds.GetLayerCount() != 1:
        gdaltest.post_reason( 'wrong number of layers' )
        return 'fail'

    lyr = gml_ds.GetLayerByName('CartographicText')

    if lyr.GetFeatureCount() != 3:
        gdaltest.post_reason( 'wrong number of features' )
        return 'fail'

    # Test 1st feature
    feat = lyr.GetNextFeature()

    if feat.GetField('featureCode') != 10198:
        gdaltest.post_reason( 'Wrong featureCode field value' )
        return 'fail'

    if feat.GetField('anchorPosition') != 8:
        gdaltest.post_reason( 'Wrong anchorPosition field value' )
        return 'fail'

    wkt = 'POINT (347243.85 461299.5)'

    if ogrtest.check_feature_geometry( feat, wkt):
        return 'fail'

    # Test 2nd feature
    feat = lyr.GetNextFeature()

    if feat.GetField('featureCode') != 10069:
        gdaltest.post_reason( 'Wrong featureCode field value' )
        return 'fail'

    wkt = 'POINT (347251.45 461250.85)'

    if ogrtest.check_feature_geometry( feat, wkt):
        return 'fail'

    return 'success'


###############################################################################
# Test of read GML file that triggeered bug #2349

def ogr_gml_5():

    if not gdaltest.have_gml_reader:
        return 'skip'

    gml_ds = ogr.Open( 'data/ticket_2349_test_1.gml' )

    lyr = gml_ds.GetLayerByName('MyPolyline')

    lyr.SetAttributeFilter( 'height > 300' )

    lyr.GetNextFeature()

    return 'success'

###############################################################################
# Test of various FIDs (various prefixes and lengths) (Ticket#1017) 
def ogr_gml_6():

    if not gdaltest.have_gml_reader: 
        return 'skip' 

    files = ['test_point1', 'test_point2', 'test_point3', 'test_point4'] 
    fids = [] 

    for filename in files: 
        fids[:] = [] 
        gml_ds = ogr.Open( 'data' + os.sep + filename + '.gml' ) 
        lyr = gml_ds.GetLayer() 
        feat = lyr.GetNextFeature() 
        while feat is not None: 
            if ( feat.GetFID() < 0 ) or ( feat.GetFID() in fids ): 
                os.remove( 'data' + os.sep + filename + '.gfs' ) 
                gdaltest.post_reason( 'Wrong FID value' ) 
                return 'fail' 
            fids.append(feat.GetFID()) 
            feat = lyr.GetNextFeature() 
        os.remove( 'data' + os.sep + filename + '.gfs' ) 

    return 'success'

###############################################################################
# Test of colon terminated prefixes for attribute values (Ticket#2493)

def ogr_gml_7():

    if not gdaltest.have_gml_reader:
        return 'skip'

    gdal.SetConfigOption('GML_EXPOSE_FID', 'FALSE')
    gml_ds = ogr.Open( 'data/test_point.gml' )
    gdal.SetConfigOption('GML_EXPOSE_FID', None)
    lyr = gml_ds.GetLayer()
    ldefn = lyr.GetLayerDefn()

    # Test fix for #2969
    if lyr.GetFeatureCount() != 5:
        gdaltest.post_reason( 'Bad feature count' )
        return 'fail'

    try:
        ldefn.GetFieldDefn(0).GetFieldTypeName
    except:
        return 'skip'

    if ldefn.GetFieldDefn(0).GetFieldTypeName(ldefn.GetFieldDefn(0).GetType())\
       != 'Real':
        return 'fail'
    if ldefn.GetFieldDefn(1).GetFieldTypeName(ldefn.GetFieldDefn(1).GetType())\
       != 'Integer':
        return 'fail'
    if ldefn.GetFieldDefn(2).GetFieldTypeName(ldefn.GetFieldDefn(2).GetType())\
       != 'String':
        return 'fail'

    return 'success'

###############################################################################
# Test a GML file with some non-ASCII UTF-8 content that triggered a bug (Ticket#2948)

def ogr_gml_8():

    if not gdaltest.have_gml_reader:
        return 'skip'

    gml_ds = ogr.Open( 'data/utf8.gml' )
    lyr = gml_ds.GetLayer()
    feat = lyr.GetNextFeature()
    if feat.GetFieldAsString('name') != '\xc4\x80liamanu':
        print(feat.GetFieldAsString('name'))
        return 'fail'

    gml_ds.Destroy()

    return 'success'

###############################################################################
# Test writing invalid UTF-8 content in a GML file (ticket #2971)

def ogr_gml_9():

    if not gdaltest.have_gml_reader:
        return 'skip'

    drv = ogr.GetDriverByName('GML')
    ds = drv.CreateDataSource('tmp/broken_utf8.gml')
    lyr = ds.CreateLayer('test')
    lyr.CreateField(ogr.FieldDefn('test', ogr.OFTString))

    dst_feat = ogr.Feature( lyr.GetLayerDefn() )
    dst_feat.SetField('test', '\x80bad')

    # Avoid the warning
    gdal.PushErrorHandler('CPLQuietErrorHandler')
    ret = lyr.CreateFeature( dst_feat )
    gdal.PopErrorHandler()

    if ret != 0:
        gdaltest.post_reason('CreateFeature failed.')
        return 'fail'

    dst_feat.Destroy()
    ds.Destroy()

    ds = ogr.Open('tmp/broken_utf8.gml')
    lyr = ds.GetLayerByName('test')
    feat = lyr.GetNextFeature()

    if feat.GetField('test') != '?bad':
        gdaltest.post_reason('Unexpected content.')
        return 'fail'

    feat.Destroy();
    ds.Destroy()

    os.remove('tmp/broken_utf8.gml')
    os.remove('tmp/broken_utf8.xsd')

    return 'success'

###############################################################################
# Test writing different data types in a GML file (ticket #2857)
# TODO: Add test for other data types as they are added to the driver.

def ogr_gml_10():

    if not gdaltest.have_gml_reader:
        return 'skip'

    drv = ogr.GetDriverByName('GML')
    ds = drv.CreateDataSource('tmp/fields.gml')
    lyr = ds.CreateLayer('test')
    field_defn = ogr.FieldDefn('string', ogr.OFTString)
    field_defn.SetWidth(100)
    lyr.CreateField(field_defn)
    lyr.CreateField(ogr.FieldDefn('date', ogr.OFTDate))
    field_defn = ogr.FieldDefn('real', ogr.OFTReal)
    field_defn.SetWidth(4)
    field_defn.SetPrecision(2)
    lyr.CreateField(field_defn)
    lyr.CreateField(ogr.FieldDefn('float', ogr.OFTReal))
    field_defn = ogr.FieldDefn('integer', ogr.OFTInteger)
    field_defn.SetWidth(5)
    lyr.CreateField(field_defn)

    dst_feat = ogr.Feature( lyr.GetLayerDefn() )
    dst_feat.SetField('string', 'test string of length 24')
    dst_feat.SetField('date', '2003/04/22')
    dst_feat.SetField('real', 12.34)
    dst_feat.SetField('float', 1234.5678)
    dst_feat.SetField('integer', '1234')

    ret = lyr.CreateFeature( dst_feat )

    if ret != 0:
        gdaltest.post_reason('CreateFeature failed.')
        return 'fail'

    dst_feat.Destroy()
    ds.Destroy()

    ds = ogr.Open('tmp/fields.gml')
    lyr = ds.GetLayerByName('test')
    feat = lyr.GetNextFeature()

    if feat.GetFieldDefnRef(feat.GetFieldIndex('string')).GetType() != ogr.OFTString:
        gdaltest.post_reason('String type is reported wrong. Got ' + str(feat.GetFieldDefnRef(feat.GetFieldIndex('string')).GetType()))
        return 'fail'
    if feat.GetFieldDefnRef(feat.GetFieldIndex('date')).GetType() != ogr.OFTString:
        gdaltest.post_reason('Date type is not reported as OFTString. Got ' + str(feat.GetFieldDefnRef(feat.GetFieldIndex('date')).GetType()))
        return 'fail'
    if feat.GetFieldDefnRef(feat.GetFieldIndex('real')).GetType() != ogr.OFTReal:
        gdaltest.post_reason('Real type is reported wrong. Got ' + str(feat.GetFieldDefnRef(feat.GetFieldIndex('real')).GetType()))
        return 'fail'
    if feat.GetFieldDefnRef(feat.GetFieldIndex('float')).GetType() != ogr.OFTReal:
        gdaltest.post_reason('Float type is not reported as OFTReal. Got ' + str(feat.GetFieldDefnRef(feat.GetFieldIndex('float')).GetType()))
        return 'fail'
    if feat.GetFieldDefnRef(feat.GetFieldIndex('integer')).GetType() != ogr.OFTInteger:
        gdaltest.post_reason('Integer type is reported wrong. Got ' + str(feat.GetFieldDefnRef(feat.GetFieldIndex('integer')).GetType()))
        return 'fail'

    if feat.GetField('string') != 'test string of length 24':
        gdaltest.post_reason('Unexpected string content.' + feat.GetField('string') )
        return 'fail'
    if feat.GetField('date') != '2003/04/22':
        gdaltest.post_reason('Unexpected string content.' + feat.GetField('date') )
        return 'fail'
    if feat.GetFieldAsDouble('real') != 12.34:
        gdaltest.post_reason('Unexpected real content.')
        return 'fail'
    if feat.GetField('float') != 1234.5678:
        gdaltest.post_reason('Unexpected float content.')
        return 'fail'
    if feat.GetField('integer') != 1234:
        gdaltest.post_reason('Unexpected integer content.')
        return 'fail'

    if lyr.GetLayerDefn().GetFieldDefn(lyr.GetLayerDefn().GetFieldIndex('string')).GetWidth() != 100:
        gdaltest.post_reason('Unexpected width of string field.')
        return 'fail'
    if lyr.GetLayerDefn().GetFieldDefn(lyr.GetLayerDefn().GetFieldIndex('real')).GetWidth() != 4:
        gdaltest.post_reason('Unexpected width of real field.')
        return 'fail'
    if lyr.GetLayerDefn().GetFieldDefn(lyr.GetLayerDefn().GetFieldIndex('real')).GetPrecision() != 2:
        gdaltest.post_reason('Unexpected precision of real field.')
        return 'fail'
    if lyr.GetLayerDefn().GetFieldDefn(lyr.GetLayerDefn().GetFieldIndex('integer')).GetWidth() != 5:
        gdaltest.post_reason('Unexpected width of integer field.')
        return 'fail'

    feat.Destroy();
    ds.Destroy()

    os.remove('tmp/fields.gml')
    os.remove('tmp/fields.xsd')

    return 'success'

###############################################################################
# Test reading a geometry element specified with <GeometryElementPath>

def ogr_gml_11():

    if not gdaltest.have_gml_reader:
        return 'skip'

    # Make sure the .gfs file is more recent that the .gml one
    try:
        gml_mtime = os.stat('data/testgeometryelementpath.gml').st_mtime
        gfs_mtime = os.stat('data/testgeometryelementpath.gfs').st_mtime
        touch_gfs = gfs_mtime <= gml_mtime
    except:
        touch_gfs = True
    if touch_gfs:
        print('Touching .gfs file')
        f = open('data/testgeometryelementpath.gfs', 'rb+')
        data = f.read(1)
        f.seek(0, 0)
        f.write(data)
        f.close()

    ds = ogr.Open('data/testgeometryelementpath.gml')
    lyr = ds.GetLayer(0)
    if lyr.GetGeometryColumn() != 'location1container|location1':
        gdaltest.post_reason('did not get expected geometry column name')
        return 'fail'

    feat = lyr.GetNextFeature()
    if feat.GetField('attrib1') != 'attrib1_value':
        gdaltest.post_reason('did not get expected value for attrib1')
        return 'fail'
    if feat.GetField('attrib2') != 'attrib2_value':
        gdaltest.post_reason('did not get expected value for attrib2')
        return 'fail'
    geom = feat.GetGeometryRef()
    if geom.ExportToWkt() != 'POINT (3 50)':
        gdaltest.post_reason('did not get expected geometry')
        return 'fail'
    ds = None
    return 'success'

###############################################################################
# Test reading a virtual GML file

def ogr_gml_12():

    if not gdaltest.have_gml_reader:
        return 'skip'

    ds = ogr.Open('/vsizip/data/testgeometryelementpath.zip/testgeometryelementpath.gml')
    lyr = ds.GetLayer(0)
    if lyr.GetGeometryColumn() != 'location1container|location1':
        gdaltest.post_reason('did not get expected geometry column name')
        return 'fail'

    feat = lyr.GetNextFeature()
    if feat.GetField('attrib1') != 'attrib1_value':
        gdaltest.post_reason('did not get expected value for attrib1')
        return 'fail'
    if feat.GetField('attrib2') != 'attrib2_value':
        gdaltest.post_reason('did not get expected value for attrib2')
        return 'fail'
    geom = feat.GetGeometryRef()
    if geom.ExportToWkt() != 'POINT (3 50)':
        gdaltest.post_reason('did not get expected geometry')
        return 'fail'
    ds = None
    return 'success'

###############################################################################
# Test reading GML with StringList, IntegerList and RealList fields

def ogr_gml_13():

    if not gdaltest.have_gml_reader:
        return 'skip'
    
    ds = ogr.Open('data/testlistfields.gml')
    lyr = ds.GetLayer(0)
    feat = lyr.GetNextFeature()
    if feat.GetFieldAsStringList(feat.GetFieldIndex('attrib1')) != ['value1','value2']:
        gdaltest.post_reason('did not get expected value for attrib1')
        return 'fail'
    if feat.GetField(feat.GetFieldIndex('attrib2')) != 'value3':
        gdaltest.post_reason('did not get expected value for attrib2')
        return 'fail'
    if feat.GetFieldAsIntegerList(feat.GetFieldIndex('attrib3')) != [4,5]:
        gdaltest.post_reason('did not get expected value for attrib3')
        return 'fail'
    if feat.GetFieldAsDoubleList(feat.GetFieldIndex('attrib4')) != [6.1,7.1]:
        gdaltest.post_reason('did not get expected value for attrib4')
        return 'fail'
    ds = None
    return 'success'

###############################################################################
# Test xlink resolution

def ogr_gml_14():

    if not gdaltest.have_gml_reader:
        return 'skip'

    # We need CURL for xlink resolution, and a sign that Curl is available
    # is the availability of the WMS driver
    try:
        gdaltest.wms_drv = gdal.GetDriverByName( 'WMS' )
    except:
        gdaltest.wms_drv = None
    if gdaltest.wms_drv is None:
        return 'skip'


    files = [ 'xlink1.gml', 'xlink2.gml', 'expected1.gml', 'expected2.gml' ]
    for file in files:
        if not gdaltest.download_file('http://download.osgeo.org/gdal/data/gml/' + file, file ):
            return 'skip'

    gdal.SetConfigOption( 'GML_SKIP_RESOLVE_ELEMS', 'NONE' )
    gdal.SetConfigOption( 'GML_SAVE_RESOLVED_TO', 'tmp/cache/xlink1resolved.gml' )
    gml_ds = ogr.Open( 'tmp/cache/xlink1.gml' )
    gml_ds = None
    gdal.SetConfigOption( 'GML_SKIP_RESOLVE_ELEMS', 'gml:directedNode' )
    gdal.SetConfigOption( 'GML_SAVE_RESOLVED_TO', 'tmp/cache/xlink2resolved.gml' )
    gml_ds = ogr.Open( 'tmp/cache/xlink1.gml' )
    gml_ds = None
    gdal.SetConfigOption( 'GML_SKIP_RESOLVE_ELEMS', 'ALL' )

    try:
        fp = open( 'tmp/cache/xlink1resolved.gml', 'r' )
        text = fp.read()
        fp.close()
        os.remove( 'tmp/cache/xlink1resolved.gml' )
        fp = open( 'tmp/cache/expected1.gml', 'r' )
        expectedtext = fp.read()
        fp.close()
    except:
        return 'fail'

    if text != expectedtext:
        print('Problem with file 1')
        return 'fail'

    try:
        fp = open( 'tmp/cache/xlink2resolved.gml', 'r' )
        text = fp.read()
        fp.close()
        os.remove( 'tmp/cache/xlink2resolved.gml' )
        fp = open( 'tmp/cache/expected2.gml', 'r' )
        expectedtext = fp.read()
        fp.close()
    except:
        return 'fail'

    if text != expectedtext:
        print('Problem with file 2')
        return 'fail'

    return 'success'

###############################################################################
# Run test_ogrsf

def ogr_gml_15():

    if not gdaltest.have_gml_reader:
        return 'skip'

    import test_cli_utilities
    if test_cli_utilities.get_test_ogrsf_path() is None:
        return 'skip'

    ret = gdaltest.runexternal(test_cli_utilities.get_test_ogrsf_path() + ' -ro data/test_point.gml')

    if ret.find('INFO') == -1 or ret.find('ERROR') != -1:
        print(ret)
        return 'fail'

    return 'success'

###############################################################################
# Read CityGML generic attributes

def ogr_gml_16():

    if not gdaltest.have_gml_reader:
        return 'skip'

    ds = ogr.Open('data/citygml.gml')
    lyr = ds.GetLayer(0)
    feat = lyr.GetNextFeature()

    if feat.GetField('Name_') != 'aname' or \
       feat.GetField('a_int_attr') != 2 or \
       feat.GetField('a_double_attr') != 3.45:
        feat.DumpReadable()
        gdaltest.post_reason('did not get expected values')
        return 'fail'

    return 'success'

###############################################################################
# Read layer SRS for WFS 1.0.0 return

def ogr_gml_17():

    if not gdaltest.have_gml_reader:
        return 'skip'

    ds = ogr.Open('data/gnis_pop_100.gml')
    lyr = ds.GetLayer(0)
    sr = lyr.GetSpatialRef()
    got_wkt = sr.ExportToWkt()
    if got_wkt.find('GEOGCS["WGS 84"') == -1:
        gdaltest.post_reason('did not get expected SRS')
        print(got_wkt)
        return 'fail'

    feat = lyr.GetNextFeature()
    geom = feat.GetGeometryRef()
    got_wkt = geom.ExportToWkt()
    if got_wkt != 'POINT (2.09 34.12)':
        gdaltest.post_reason('did not get expected geometry')
        print(got_wkt)
        return 'fail'

    return 'success'

###############################################################################
# Read layer SRS for WFS 1.1.0 return

def ogr_gml_18():

    if not gdaltest.have_gml_reader:
        return 'skip'

    ds = ogr.Open('data/gnis_pop_110.gml')
    lyr = ds.GetLayer(0)
    sr = lyr.GetSpatialRef()
    got_wkt = sr.ExportToWkt()
    if got_wkt.find('GEOGCS["WGS 84"') == -1:
        gdaltest.post_reason('did not get expected SRS')
        print(got_wkt)
        return 'fail'

    feat = lyr.GetNextFeature()
    geom = feat.GetGeometryRef()
    got_wkt = geom.ExportToWkt()
    if got_wkt != 'POINT (2.09 34.12)':
        gdaltest.post_reason('did not get expected geometry')
        print(got_wkt)
        return 'fail'

    return 'success'

###############################################################################
# Read layer SRS for WFS 1.1.0 return, but without trying to restore
# (long, lat) order. So we should get EPSGA:4326 and (lat, long) order

def ogr_gml_19():

    if not gdaltest.have_gml_reader:
        return 'skip'

    try:
        os.remove( 'data/gnis_pop_110.gfs' )
    except:
        pass

    gdal.SetConfigOption('GML_INVERT_AXIS_ORDER_IF_LAT_LONG', 'NO')
    ds = ogr.Open('data/gnis_pop_110.gml')
    gdal.SetConfigOption('GML_INVERT_AXIS_ORDER_IF_LAT_LONG', None)

    lyr = ds.GetLayer(0)
    sr = lyr.GetSpatialRef()
    got_wkt = sr.ExportToWkt()
    if got_wkt.find('GEOGCS["WGS 84"') == -1 or \
       got_wkt.find('AXIS["Latitude",NORTH],AXIS["Longitude",EAST]') == -1:
        gdaltest.post_reason('did not get expected SRS')
        print(got_wkt)
        return 'fail'

    feat = lyr.GetNextFeature()
    geom = feat.GetGeometryRef()
    got_wkt = geom.ExportToWkt()
    if got_wkt != 'POINT (34.12 2.09)':
        gdaltest.post_reason('did not get expected geometry')
        print(got_wkt)
        return 'fail'

    return 'success'

###############################################################################
# Test parsing a .xsd where the type definition is before its reference

def ogr_gml_20():

    if not gdaltest.have_gml_reader:
        return 'skip'

    try:
        os.remove( 'data/archsites.gfs' )
    except:
        pass

    ds = ogr.Open('data/archsites.gml')
    lyr = ds.GetLayer(0)
    ldefn = lyr.GetLayerDefn()

    try:
        ldefn.GetFieldDefn(0).GetFieldTypeName
    except:
        return 'skip'

    idx = ldefn.GetFieldIndex("gml_id")
    if idx == -1:
        gdaltest.post_reason('did not get expected column "gml_id"')
        return 'fail'

    idx = ldefn.GetFieldIndex("cat")
    fddefn = ldefn.GetFieldDefn(idx)
    if fddefn.GetFieldTypeName(fddefn.GetType()) != 'Integer':
        gdaltest.post_reason('did not get expected column type for col "cat"')
        return 'fail'
    idx = ldefn.GetFieldIndex("str1")
    fddefn = ldefn.GetFieldDefn(idx)
    if fddefn.GetFieldTypeName(fddefn.GetType()) != 'String':
        gdaltest.post_reason('did not get expected column type for col "str1"')
        return 'fail'

    if lyr.GetGeometryColumn() != 'the_geom':
        gdaltest.post_reason('did not get expected geometry column name')
        return 'fail'

    if ldefn.GetGeomType() != ogr.wkbPoint:
        gdaltest.post_reason('did not get expected geometry type')
        return 'fail'

    ds = None

    try:
        os.stat('data/archsites.gfs')
        gdaltest.post_reason('did not expected .gfs -> XSD parsing failed')
        return 'fail'
    except:
        return 'success'

###############################################################################
# Test writing GML3

def ogr_gml_21():

    if not gdaltest.have_gml_reader:
        return 'skip'

    # Create GML3 file
    sr = osr.SpatialReference()
    sr.ImportFromEPSG(4326)

    ds = ogr.GetDriverByName('GML').CreateDataSource('tmp/gml_21.gml', options = ['FORMAT=GML3'] )
    lyr = ds.CreateLayer('firstlayer', srs = sr)
    lyr.CreateField(ogr.FieldDefn('string_field', ogr.OFTString))

    feat = ogr.Feature(lyr.GetLayerDefn())
    geom = ogr.CreateGeometryFromWkt('POINT (2 49)')
    feat.SetGeometry(geom)
    lyr.CreateFeature(feat)

    feat = ogr.Feature(lyr.GetLayerDefn())
    feat.SetField(0, 'foo')
    geom = ogr.CreateGeometryFromWkt('POINT (3 48)')
    feat.SetGeometry(geom)
    lyr.CreateFeature(feat)

    ds = None

    # Reopen the file
    ds = ogr.Open('tmp/gml_21.gml')
    lyr = ds.GetLayer(0)
    feat = lyr.GetNextFeature()
    if feat.GetGeometryRef().ExportToWkt() != 'POINT (2 49)':
        gdaltest.post_reason('did not get expected geometry')
        return 'fail'
    ds = None

    # Test that .gml and .xsd are identical to what is expected
    f1 = open('tmp/gml_21.gml', 'rt')
    f2 = open('data/expected_gml_21.gml', 'rt')
    line1 = f1.readline()
    line2 = f2.readline()
    while line1 != '':
        line1 = line1.strip()
        line2 = line2.strip()
        if line1 != line2:
            gdaltest.post_reason('.gml file not identical to expected')
            print(open('tmp/gml_21.gml', 'rt').read())
            return 'fail'
        line1 = f1.readline()
        line2 = f2.readline()
    f1.close()
    f2.close()

    f1 = open('tmp/gml_21.xsd', 'rt')
    f2 = open('data/expected_gml_21.xsd', 'rt')
    line1 = f1.readline()
    line2 = f2.readline()
    while line1 != '':
        line1 = line1.strip()
        line2 = line2.strip()
        if line1 != line2:
            gdaltest.post_reason('.xsd file not identical to expected')
            print(open('tmp/gml_21.xsd', 'rt').read())
            return 'fail'
        line1 = f1.readline()
        line2 = f2.readline()
    f1.close()
    f2.close()

    return 'success'

###############################################################################
# Read a OpenLS DetermineRouteResponse document

def ogr_gml_22():

    if not gdaltest.have_gml_reader:
        return 'skip'

    ds = ogr.Open('data/paris_typical_strike_demonstration.xml')
    lyr = ds.GetLayerByName('RouteGeometry')
    if lyr is None:
        gdaltest.post_reason('cannot find RouteGeometry')
        return 'fail'
    lyr = ds.GetLayerByName('RouteInstruction')
    if lyr is None:
        gdaltest.post_reason('cannot find RouteInstruction')
        return 'fail'
    count = lyr.GetFeatureCount()
    if count != 9:
        gdaltest.post_reason('did not get expected feature count')
        print(count)
        return 'fail'

    ds = None
    return 'success'

###############################################################################
# Test that use SRS defined in global gml:Envelope if no SRS is set for any
# feature geometry

def ogr_gml_23():

    if not gdaltest.have_gml_reader:
        return 'skip'

    try:
        os.remove( 'tmp/global_geometry.gfs' )
    except:
        pass

    shutil.copy('data/global_geometry.xml', 'tmp/global_geometry.xml')

    # Here we use only the .xml file
    ds = ogr.Open('tmp/global_geometry.xml')

    lyr = ds.GetLayer(0)
    sr = lyr.GetSpatialRef()
    got_wkt = sr.ExportToWkt()
    if got_wkt.find('GEOGCS["WGS 84"') == -1 or \
       got_wkt.find('AXIS["Latitude",NORTH],AXIS["Longitude",EAST]') != -1:
        gdaltest.post_reason('did not get expected SRS')
        print(got_wkt)
        return 'fail'

    feat = lyr.GetNextFeature()
    geom = feat.GetGeometryRef()
    got_wkt = geom.ExportToWkt()
    if got_wkt != 'POINT (2 49)':
        gdaltest.post_reason('did not get expected geometry')
        print(got_wkt)
        return 'fail'

    extent = lyr.GetExtent()
    if extent != (2.0, 3.0, 49.0, 50.0):
        gdaltest.post_reason('did not get expected layer extent')
        print(extent)
        return 'fail'

    return 'success'

###############################################################################
# Test that use SRS defined in global gml:Envelope if no SRS is set for any
# feature geometry

def ogr_gml_24():

    if not gdaltest.have_gml_reader:
        return 'skip'

    try:
        os.remove( 'data/global_geometry.gfs' )
    except:
        pass

    # Here we use only the .xml file and the .xsd file
    ds = ogr.Open('data/global_geometry.xml')

    lyr = ds.GetLayer(0)

    # Because we read the .xsd, we (currently) don't find the SRS
    
    #sr = lyr.GetSpatialRef()
    #got_wkt = sr.ExportToWkt()
    #if got_wkt.find('GEOGCS["WGS 84"') == -1 or \
    #   got_wkt.find('AXIS["Latitude",NORTH],AXIS["Longitude",EAST]') != -1:
    #    gdaltest.post_reason('did not get expected SRS')
    #    print(got_wkt)
    #    return 'fail'

    feat = lyr.GetNextFeature()
    geom = feat.GetGeometryRef()
    got_wkt = geom.ExportToWkt()
    if got_wkt != 'POINT (2 49)':
        gdaltest.post_reason('did not get expected geometry')
        print(got_wkt)
        return 'fail'

    extent = lyr.GetExtent()
    if extent != (2.0, 3.0, 49.0, 50.0):
        gdaltest.post_reason('did not get expected layer extent')
        print(extent)
        return 'fail'

    return 'success'

###############################################################################
# Test fixes for #3934 and #3935

def ogr_gml_25():

    if not gdaltest.have_gml_reader:
        return 'skip'

    try:
        os.remove( 'data/curveProperty.gfs' )
    except:
        pass

    ds = ogr.Open('data/curveProperty.xml')

    lyr = ds.GetLayer(0)

    feat = lyr.GetNextFeature()
    geom = feat.GetGeometryRef()
    got_wkt = geom.ExportToWkt()
    if got_wkt != 'POLYGON ((14 21,6 21,6 21,6 9,6 9,14 9,14 9,22 9,22 9,22 21,22 21,14 21))':
        gdaltest.post_reason('did not get expected geometry')
        print(got_wkt)
        return 'fail'

    return 'success'

###############################################################################
# Test writing and reading 3D geoms (GML2)

def ogr_gml_26():

    if not gdaltest.have_gml_reader:
        return 'skip'

    import test_cli_utilities
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -f GML tmp/ogr_gml_26.gml data/poly.shp -zfield eas_id')

    f = open('tmp/ogr_gml_26.gml', 'rt')
    content = f.read()
    f.close()
    if content.find("<gml:coord><gml:X>478315.53125</gml:X><gml:Y>4762880.5</gml:Y><gml:Z>158</gml:Z></gml:coord>") == -1:
        return 'fail'

    ds = ogr.Open('tmp/ogr_gml_26.gml')

    lyr = ds.GetLayer(0)

    if lyr.GetGeomType() != ogr.wkbPolygon25D:
        return 'fail'

    ds = None

    return 'success'

###############################################################################
# Test writing and reading 3D geoms (GML3)

def ogr_gml_27():

    if not gdaltest.have_gml_reader:
        return 'skip'

    import test_cli_utilities
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -f GML tmp/ogr_gml_27.gml data/poly.shp -zfield eas_id -dsco FORMAT=GML3')

    f = open('tmp/ogr_gml_27.gml', 'rt')
    content = f.read()
    f.close()
    if content.find("<gml:lowerCorner>478315.53125 4762880.5 158</gml:lowerCorner>") == -1:
        return 'fail'

    ds = ogr.Open('tmp/ogr_gml_27.gml')

    lyr = ds.GetLayer(0)

    if lyr.GetGeomType() != ogr.wkbPolygon25D:
        return 'fail'

    ds = None

    return 'success'

###############################################################################
# Test writing and reading layers of type wkbNone (#4154)

def ogr_gml_28():

    if not gdaltest.have_gml_reader:
        return 'skip'

    import test_cli_utilities
    if test_cli_utilities.get_ogr2ogr_path() is None:
        return 'skip'

    gdaltest.runexternal(test_cli_utilities.get_ogr2ogr_path() + ' -f GML tmp/ogr_gml_28.gml data/idlink.dbf')

    # Try with .xsd
    ds = ogr.Open('tmp/ogr_gml_28.gml')
    lyr = ds.GetLayer(0)
    if lyr.GetGeomType() != ogr.wkbNone:
        return 'fail'
    ds = None

    os.unlink('tmp/ogr_gml_28.xsd')

    ds = ogr.Open('tmp/ogr_gml_28.gml')
    lyr = ds.GetLayer(0)
    if lyr.GetGeomType() != ogr.wkbNone:
        return 'fail'
    ds = None

    # Try with .gfs
    ds = ogr.Open('tmp/ogr_gml_28.gml')
    lyr = ds.GetLayer(0)
    if lyr.GetGeomType() != ogr.wkbNone:
        return 'fail'
    ds = None

    return 'success'

###############################################################################
# Test reading FME GMLs

def ogr_gml_29():

    if not gdaltest.have_gml_reader:
        return 'skip'

    ds = ogr.Open('data/testfmegml.gml')

    expected_results = [ [ ogr.wkbMultiPoint, 'MULTIPOINT (2 49)' ],
                         [ ogr.wkbMultiPolygon, 'MULTIPOLYGON (((2 49,3 49,3 50,2 50,2 49)))'],
                         [ ogr.wkbMultiLineString, 'MULTILINESTRING ((2 49,3 50))'],
                       ]

    for j in range(len(expected_results)):
        lyr = ds.GetLayer(j)
        if lyr.GetGeomType() != expected_results[j][0]:
            gdaltest.post_reason('layer %d, did not get expected layer geometry type' % j)
            return 'fail'
        for i in range(2):
            feat = lyr.GetNextFeature()
            geom = feat.GetGeometryRef()
            got_wkt = geom.ExportToWkt()
            if got_wkt != expected_results[j][1]:
                gdaltest.post_reason('layer %d, did not get expected geometry' % j)
                print(got_wkt)
                return 'fail'

    ds = None

    return 'success'
    
###############################################################################
#  Cleanup

def ogr_gml_cleanup():
    if not gdaltest.have_gml_reader:
        return 'skip'
    
    gdaltest.clean_tmp()
    
    return ogr_gml_clean_files()


def ogr_gml_clean_files():
    try:
        os.remove( 'data/bom.gfs' )
    except:
        pass
    try:
        os.remove( 'data/utf8.gfs' )
    except:
        pass
    try:
        os.remove( 'data/ticket_2349_test_1.gfs' )
    except:
        pass
    try:
        os.remove( 'data/citygml.gfs' )
    except:
        pass
    try:
        os.remove( 'data/gnis_pop_100.gfs' )
    except:
        pass
    try:
        os.remove( 'data/gnis_pop_110.gfs' )
    except:
        pass
    try:
        os.remove( 'data/paris_typical_strike_demonstration.gfs' )
    except:
        pass
    try:
        os.remove( 'data/global_geometry.gfs' )
    except:
        pass
    try:
        os.remove( 'tmp/global_geometry.gfs' )
    except:
        pass
    try:
        os.remove( 'tmp/global_geometry.xml' )
    except:
        pass
    try:
        os.remove( 'data/curveProperty.gfs' )
    except:
        pass
    try:
        os.remove( 'tmp/ogr_gml_26.gml' )
        os.remove( 'tmp/ogr_gml_26.xsd' )
    except:
        pass
    try:
        os.remove( 'tmp/ogr_gml_27.gml' )
        os.remove( 'tmp/ogr_gml_27.xsd' )
    except:
        pass
    try:
        os.remove( 'tmp/ogr_gml_28.gml' )
        os.remove( 'tmp/ogr_gml_28.gfs' )
    except:
        pass

    files = os.listdir('data')
    for filename in files:
        if len(filename) > 13 and filename[-13:] == '.resolved.gml':
            os.unlink('data/' + filename)

    return 'success'

gdaltest_list = [ 
    ogr_gml_clean_files,
    ogr_gml_1,
    ogr_gml_2,
    ogr_gml_3,
    ogr_gml_4,
    ogr_gml_5,
    ogr_gml_6,
    ogr_gml_7,
    ogr_gml_8,
    ogr_gml_9,
    ogr_gml_10,
    ogr_gml_11,
    ogr_gml_12,
    ogr_gml_13,
    ogr_gml_14,
    ogr_gml_15,
    ogr_gml_16,
    ogr_gml_17,
    ogr_gml_18,
    ogr_gml_19,
    ogr_gml_20,
    ogr_gml_21,
    ogr_gml_22,
    ogr_gml_23,
    ogr_gml_24,
    ogr_gml_25,
    ogr_gml_26,
    ogr_gml_27,
    ogr_gml_28,
    ogr_gml_29,
    ogr_gml_cleanup ]

if __name__ == '__main__':

    gdaltest.setup_run( 'ogr_gml_read' )

    gdaltest.run_tests( gdaltest_list )

    gdaltest.summarize()

