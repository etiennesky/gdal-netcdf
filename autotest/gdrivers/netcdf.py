#!/usr/bin/env python
###############################################################################
# $Id$
#
# Project:  GDAL/OGR Test Suite
# Purpose:  Test NetCDF driver support.
# Author:   Frank Warmerdam <warmerdam@pobox.com>
# 
###############################################################################
# Copyright (c) 2007, Frank Warmerdam <warmerdam@pobox.com>
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
import osr

sys.path.append( '../pymod' )

import gdaltest

###############################################################################
# Netcdf Functions
###############################################################################

###############################################################################
# Get netcdf version and test for supported files

def netcdf_setup():

    gdaltest.netcdf_drv_version = None
    gdaltest.netcdf_drv_has_nc2 = False
    gdaltest.netcdf_drv_has_nc4 = False
    gdaltest.netcdf_drv = gdal.GetDriverByName( 'NETCDF' )

    if gdaltest.netcdf_drv is None:
        print('NOTICE: netcdf not supported, skipping checks')
        return 'skip'

    #ugly hack to get netcdf version with 'ncdump', available in netcdf v3 avnd v4
    try:
        (ret, err) = gdaltest.runexternal_out_and_err('ncdump -h')
#        (ret, err) = gdaltest.runexternal_out_and_err('LD_LIBRARY_PATH=/home/src/netcdf-test/usr/lib:$LD_LIBRARY_PATH /home/src/netcdf-test/usr/bin/ncdump -h')
    except:
        #nothing is supported as ncdump not found
        print('NOTICE: netcdf version not found')
        return 'success'
    
    i = err.find('netcdf library version ')
    #version not found
    if i == -1:
        print('NOTICE: netcdf version not found')
        return 'success'

    #netcdf library version "3.6.3" of Dec 22 2009 06:10:17 $
    #netcdf library version 4.1.1 of Mar  4 2011 12:52:19 $
    v = err[ i+23 : ]
    v = v[ 0 : v.find(' ') ]
    v = v.strip('"');

    gdaltest.netcdf_drv_version = v

    #for version 3, assume nc2 supported and nc4 unsupported
    if v[0] == '3':
        gdaltest.netcdf_drv_has_nc2 = True
        gdaltest.netcdf_drv_has_nc4 = False

    # for version 4, use nc-config to test
    elif v[0] == '4':
    
        #check if netcdf library has nc2 (64-bit) support
        #this should be done at configure time by gdal like in cdo
        try:
            ret = gdaltest.runexternal('nc-config --has-nc2')
        except:
            gdaltest.netcdf_drv_has_nc2 = False
        else:
            #should test this on windows
            if ret.rstrip() == 'yes':
                gdaltest.netcdf_drv_has_nc2 = True
                    
        #check if netcdf library has nc4 support
        #this should be done at configure time by gdal like in cdo
        try:
            ret = gdaltest.runexternal('nc-config --has-nc4')
        except:
            gdaltest.netcdf_drv_has_nc4 = False
        else:
            #should test this on windows
            if ret.rstrip() == 'yes':
                gdaltest.netcdf_drv_has_nc4 = True

    print('NOTICE: using netcdf version ' + gdaltest.netcdf_drv_version+'  has_nc2: '+str(gdaltest.netcdf_drv_has_nc2)+'  has_nc4: '+str(gdaltest.netcdf_drv_has_nc4))
 
    return 'success'

###############################################################################
#test integrity of file copy

def netcdf_test_file_copy( ifile, tmpfile, driver_name ):

    #print 'netcdf_test_file_copy( '+ifile+', '+tmpfile+', '+driver_name+' )'

    #load driver(s)
    if gdaltest.netcdf_drv is None:
        return 'skip'
    driver = gdal.GetDriverByName( driver_name )
    if driver is None:
        return 'skip'

    #create copy
    src_ds = gdal.Open( ifile )
    if src_ds is None:
        gdaltest.post_reason( 'Failed to open file '+ifile )
        return 'fail'
    dst_ds = driver.CreateCopy( tmpfile, src_ds )
    if dst_ds is None:
        gdaltest.post_reason( 'Failed to create file '+tmpfile )
        return 'fail'
    dst_ds = None
    dst_ds = gdal.Open( tmpfile )

    #do some tests
    #print str(src_ds.GetGeoTransform())+'-'+str(dst_ds.GetGeoTransform())
    if src_ds.GetGeoTransform() != dst_ds.GetGeoTransform():
        gdaltest.post_reason( 'Incorrect geotransform, got '+str(dst_ds.GetGeoTransform())+' for '+tmpfile )
        return 'fail'
    
    #get projection in PROJ.4 and WKT format
    srs = osr.SpatialReference()
    src_wkt = src_ds.GetProjection()
    srs.ImportFromWkt(src_wkt)
    src_proj4 = srs.ExportToProj4()
    dst_wkt = dst_ds.GetProjection()
    srs.ImportFromWkt(src_wkt)
    dst_proj4 = srs.ExportToProj4()

    #check projection in WKT format - don't cause test to fail, as it can be too stringent
    #print src_wkt+'-'+dst_wkt
    if src_wkt != dst_wkt:
        print('WARNING: Possibly incorrect projection in file '+tmpfile+', got '+dst_wkt+' - will also check PROJ.4 string')
        #return 'fail'

    #check projection in PROJ.4 format - allow this test to fail
    #print src_proj4+'-'+dst_proj4
    if src_proj4 != dst_proj4:
        gdaltest.post_reason( 'Incorrect projection, got '+dst_proj4+' for '+tmpfile )
        return 'fail'

    if src_ds.GetRasterBand(1).Checksum() != dst_ds.GetRasterBand(1).Checksum():
        gdaltest.post_reason( 'Incorrect checksum, got'+dst_ds.GetRasterBand(1).Checksum()+' for '+tmpfile )
        return 'fail'

    src_ds = None
    dst_ds = None
    #don't cleanup as we could use the tmpfile later, cleanup after all tests finished
    #gdaltest.clean_tmp()

    return 'success'

###############################################################################
# Netcdf Tests
###############################################################################

###############################################################################
# Perform simple read test.

def netcdf_1():

    #setup netcdf environment
    netcdf_setup()

    if gdaltest.netcdf_drv is None:
        return 'skip'

    tst = gdaltest.GDALTest( 'NetCDF', 'NETCDF:"data/bug636.nc":tas', 1, 31621,
                             filename_absolute = 1 )

    # We don't want to gum up the test stream output with the
    # 'Warning 1: No UNIDATA NC_GLOBAL:Conventions attribute' message.
    gdal.PushErrorHandler( 'CPLQuietErrorHandler' )
    result = tst.testOpen()
    gdal.PopErrorHandler()

    return result

###############################################################################
# Verify a simple createcopy operation.  We can't do the trivial gdaltest
# operation because the new file will only be accessable via subdatasets!

def netcdf_2():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    src_ds = gdal.Open( 'data/byte.tif' )
    
    base_ds = gdaltest.netcdf_drv.CreateCopy( 'tmp/netcdf2.nc', src_ds)

    tst = gdaltest.GDALTest( 'NetCDF', 'tmp/netcdf2.nc',
                             1, 4672,
                             filename_absolute = 1 )

    wkt = """PROJCS["NAD27 / UTM zone 11N",
    GEOGCS["NAD27",
        DATUM["North_American_Datum_1927",
            SPHEROID["Clarke 1866",6378206.4,294.9786982139006,
                AUTHORITY["EPSG","7008"]],
            AUTHORITY["EPSG","6267"]],
        PRIMEM["Greenwich",0],
        UNIT["degree",0.0174532925199433],
        AUTHORITY["EPSG","4267"]],
    PROJECTION["Transverse_Mercator"],
    PARAMETER["latitude_of_origin",0],
    PARAMETER["central_meridian",-117],
    PARAMETER["scale_factor",0.9996],
    PARAMETER["false_easting",500000],
    PARAMETER["false_northing",0],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    AUTHORITY["EPSG","26711"]]"""

    result = tst.testOpen( check_prj = wkt )

    if result != 'success':
        return result

    base_ds = None
    gdaltest.clean_tmp()

    return 'success'

###############################################################################

def netcdf_3():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    ds = gdal.Open( 'data/sombrero.grd' )
    bnd = ds.GetRasterBand(1)
    minmax = bnd.ComputeRasterMinMax()

    if abs(minmax[0] - (-0.675758)) > 0.000001 or abs(minmax[1] - 1.0) > 0.000001:
        gdaltest.post_reason( 'Wrong min or max.' )
        return 'fail'

    bnd = None
    ds = None

    return 'success'
    
###############################################################################
# In #2582 5dimensional files were causing problems.  Verify use ok.

def netcdf_4():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    tst = gdaltest.GDALTest( 'NetCDF',
                             'NETCDF:data/foo_5dimensional.nc:temperature',
                             3, 1218, filename_absolute = 1 )

    # We don't want to gum up the test stream output with the
    # 'Warning 1: No UNIDATA NC_GLOBAL:Conventions attribute' message.
    gdal.PushErrorHandler( 'CPLQuietErrorHandler' )
    result = tst.testOpen()
    gdal.PopErrorHandler()

    return result
    
###############################################################################
# In #2583 5dimensional files were having problems unrolling the highest
# dimension - check handling now on band 7.

def netcdf_5():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    tst = gdaltest.GDALTest( 'NetCDF',
                             'NETCDF:data/foo_5dimensional.nc:temperature',
                             7, 1227, filename_absolute = 1 )

    # We don't want to gum up the test stream output with the
    # 'Warning 1: No UNIDATA NC_GLOBAL:Conventions attribute' message.
    gdal.PushErrorHandler( 'CPLQuietErrorHandler' )
    result = tst.testOpen()
    gdal.PopErrorHandler()

    return result

###############################################################################
#ticket #3324 check spatial reference reading for cf-1.4 lambert conformal
#1 standard parallel.
def netcdf_6():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    ds = gdal.Open( 'data/cf_lcc1sp.nc' )
    prj = ds.GetProjection( )

    sr = osr.SpatialReference( )
    sr.ImportFromWkt( prj )
    lat_origin = sr.GetProjParm( 'latitude_of_origin' )

    if lat_origin != 25:
        gdaltest.post_reason( 'Latitude of origin does not match expected:\n%f' 
                              % lat_origin )
        return 'fail'

    ds = None

    return 'success'
    
###############################################################################
#ticket #3324 check spatial reference reading for cf-1.4 lambert conformal
#2 standard parallels.
def netcdf_7():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    ds = gdal.Open( 'data/cf_lcc2sp.nc' )
    prj = ds.GetProjection( )

    sr = osr.SpatialReference( )
    sr.ImportFromWkt( prj )
    std_p1 = sr.GetProjParm( 'standard_parallel_1' )
    std_p2 = sr.GetProjParm( 'standard_parallel_2' )

    if std_p1 != 33.0 or std_p2 != 45.0:
        gdaltest.post_reason( 'Standard Parallels do not match expected:\n%f,%f' 
                              % ( std_p1, std_p2 ) )
        return 'fail'

    ds = None
    sr = None

    return 'success'
    
###############################################################################
#check for cf convention read of albers equal area
def netcdf_8():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    ds = gdal.Open( 'data/cf_aea2sp_invf.nc' )
    prj = ds.GetProjection( )

    if prj != 'PROJCS["unnamed",GEOGCS["unknown",DATUM["unknown",SPHEROID["Spheroid",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",29.83333333333334],PARAMETER["standard_parallel_2",45.83333333333334],PARAMETER["latitude_of_center",37.5],PARAMETER["longitude_of_center",-96],PARAMETER["false_easting",0],PARAMETER["false_northing",0]]':

        gdaltest.post_reason( 'Projection does not match expected:\n%s' % ( prj ) )
        return 'fail'

    ds = None

    return 'success'
    
###############################################################################
#check to see if projected systems default to wgs84 if no spheroid def
def netcdf_9():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    ds = gdal.Open( 'data/cf_no_sphere.nc' )

    prj = ds.GetProjection( )

    sr = osr.SpatialReference( )
    sr.ImportFromWkt( prj )
    spheroid = sr.GetAttrValue( 'SPHEROID' )

    if spheroid != 'WGS 84':
        gdaltest.post_reason( 'Incorrect spheroid read from file\n%s' 
                              % ( spheroid ) )
        return 'fail'

    ds = None
    sr = None

    return 'success'
    
###############################################################################
#check if km pixel size makes it through to gt
def netcdf_10():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    ds = gdal.Open( 'data/cf_no_sphere.nc' )

    prj = ds.GetProjection( )

    gt = ds.GetGeoTransform( )

    if gt != (-1897186.0290038721, 
               5079.3608398440065, 
               0.0, 
               2674684.0244560046, 
               0.0, 
               -5079.4721679684635):

        gdaltest.post_reason( 'Incorrect geotransform' )
        return 'fail'

    ds = None

    return 'success'
    
###############################################################################
#check if ll gets caught in km pixel size check
def netcdf_11():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    ds = gdal.Open( 'data/cf_geog.nc' )

    gt = ds.GetGeoTransform( )

    if gt != (-0.5, 1.0, 0.0, 10.5, 0.0, -1.0):

        gdaltest.post_reason( 'Incorrect geotransform' )
        return 'fail'

    ds = None

    return 'success'

###############################################################################
#check for scale/offset set/get.
def netcdf_12():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    ds = gdal.Open( 'data/scale_offset.nc' )

    scale = ds.GetRasterBand( 1 ).GetScale();
    offset = ds.GetRasterBand( 1 ).GetOffset()

    if scale != 0.01 or offset != 1.5:
        gdaltest.post_reason( 'Incorrect scale(%f) or offset(%f)' % ( scale, offset ) )
        return 'fail'

    ds = None

    return 'success'

###############################################################################
#check for scale/offset = 1.0/0.0 if no scale or offset is available
def netcdf_13():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    ds = gdal.Open( 'data/no_scale_offset.nc' )

    scale = ds.GetRasterBand( 1 ).GetScale();
    offset = ds.GetRasterBand( 1 ).GetOffset()

    if scale != 1.0 or offset != 0.0:
        gdaltest.post_reason( 'Incorrect scale or offset' )
        return 'fail'

    ds = None

    return 'success'

###############################################################################
#check for scale/offset for two variables
def netcdf_14():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    ds = gdal.Open( 'NETCDF:data/two_vars_scale_offset.nc:z' )

    scale = ds.GetRasterBand( 1 ).GetScale();
    offset = ds.GetRasterBand( 1 ).GetOffset()

    if scale != 0.01 or offset != 1.5:
        gdaltest.post_reason( 'Incorrect scale(%f) or offset(%f)' % ( scale, offset ) )
        return 'fail'
    
    ds = None

    ds = gdal.Open( 'NETCDF:data/two_vars_scale_offset.nc:q' )

    scale = ds.GetRasterBand( 1 ).GetScale();
    offset = ds.GetRasterBand( 1 ).GetOffset()

    scale = ds.GetRasterBand( 1 ).GetScale();
    offset = ds.GetRasterBand( 1 ).GetOffset()

    if scale != 0.1 or offset != 2.5:
        gdaltest.post_reason( 'Incorrect scale(%f) or offset(%f)' % ( scale, offset ) )
        return 'fail'

    return 'success'

###############################################################################
#check support for netcdf-2 (64 bit)
def netcdf_15():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    if gdaltest.netcdf_drv_has_nc2:
        ds = gdal.Open( 'data/trmm-nc2.nc' )
        if ds is None:
            return 'fail'
        else:
            ds = None
            return 'success'
    else:
        return 'skip'

    return 'success'
        
###############################################################################
#check support for netcdf-4
def netcdf_16():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    ifile = 'data/trmm-nc4.nc'
                
    if gdaltest.netcdf_drv_has_nc4:

        # test with Open()
        ds = gdal.Open( ifile )
        if ds is None:
            return 'fail'
        else:
            name = ds.GetDriver().GetDescription()
            ds = None
            #return fail if did not open with the netCDF driver (i.e. HDF5Image)
            if name != 'netCDF':
                return 'fail'

        # test with Identify()
        name = gdal.IdentifyDriver( ifile ).GetDescription()
        if name != 'netCDF':
            return 'fail'

    else:
        return 'skip'

    return 'success'

###############################################################################
#check support for netcdf-4 - make sure hdf5 is not read by netcdf driver
def netcdf_17():

    if gdaltest.netcdf_drv is None:
        return 'skip'

    ifile = 'data/u8be.h5'

    #skip test if Hdf5 is not enabled
    if gdal.GetDriverByName( 'HDF5' ) is None and \
            gdal.GetDriverByName( 'HDF5Image' ) is None:
        return 'skip'
    
    if gdaltest.netcdf_drv_has_nc4:

        #test with Open()
        ds = gdal.Open( ifile )
        if ds is None:
            return 'fail'
        else:
            name = ds.GetDriver().GetDescription()
            ds = None
                #return fail if opened with the netCDF driver
            if name == 'netCDF':
                return 'fail'

        # test with Identify()
        name = gdal.IdentifyDriver( ifile ).GetDescription()
        if name == 'netCDF':
            return 'fail'

    else:
        return 'skip'
    
    return 'success'

     
###############################################################################

gdaltest_list = [
    netcdf_1,
    netcdf_2,
    netcdf_3,
    netcdf_4,
    netcdf_5,
    netcdf_6,
    netcdf_7,
    netcdf_8,
    netcdf_9,
    netcdf_10, 
    netcdf_11,
    netcdf_12,
    netcdf_13,
    netcdf_14,
    netcdf_15,
    netcdf_16,
    netcdf_17
 ]

if __name__ == '__main__':

    gdaltest.setup_run( 'netcdf' )

    gdaltest.run_tests( gdaltest_list )

    #make sure we cleanup
    gdaltest.clean_tmp()

    gdaltest.summarize()

