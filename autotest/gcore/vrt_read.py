#!/usr/bin/env python
###############################################################################
# $Id$
#
# Project:  GDAL/OGR Test Suite
# Purpose:  Test basic read support for a all datatypes from a VRT file.
# Author:   Frank Warmerdam <warmerdam@pobox.com>
# 
###############################################################################
# Copyright (c) 2003, Frank Warmerdam <warmerdam@pobox.com>
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

sys.path.append( '../pymod' )

import gdaltest
import gdal

###############################################################################
# When imported build a list of units based on the files available.

gdaltest_list = []

init_list = [ \
    ('byte.vrt', 1, 4672, None),
    ('int16.vrt', 1, 4672, None),
    ('uint16.vrt', 1, 4672, None),
    ('int32.vrt', 1, 4672, None),
    ('uint32.vrt', 1, 4672, None),
    ('float32.vrt', 1, 4672, None),
    ('float64.vrt', 1, 4672, None),
    ('cint16.vrt', 1, 5028, None),
    ('cint32.vrt', 1, 5028, None),
    ('cfloat32.vrt', 1, 5028, None),
    ('cfloat64.vrt', 1, 5028, None),
    ('msubwinbyte.vrt', 2, 2699, None),
    ('utmsmall.vrt', 1, 50054, None),
    ('byte_nearest_50pct.vrt', 1, 1192, None),
    ('byte_averaged_50pct.vrt', 1, 1152, None),
    ('byte_nearest_200pct.vrt', 1, 18784, None),
    ('byte_averaged_200pct.vrt', 1, 18784, None)]

###############################################################################
# The VRT references a non existing TIF file

def vrt_read_1():

    gdal.PushErrorHandler('CPLQuietErrorHandler')
    ds = gdal.Open('data/idontexist.vrt')
    gdal.PopErrorHandler()

    if ds is None:
        return 'success'

    return 'fail'

###############################################################################
# The VRT references a non existing TIF file, but using the proxy pool dataset API (#2837)

def vrt_read_2():

    ds = gdal.Open('data/idontexist2.vrt')
    if ds is None:
        return 'fail'

    gdal.PushErrorHandler('CPLQuietErrorHandler')
    cs = ds.GetRasterBand(1).Checksum()
    gdal.PopErrorHandler()

    if cs != 0:
        return 'fail'

    ds.GetMetadata()
    ds.GetRasterBand(1).GetMetadata()
    ds.GetGCPs()

    ds = None
    return 'success'

###############################################################################
# Test init of band data in case of cascaded VRT (ticket #2867)

def vrt_read_3():

    driver_tif = gdal.GetDriverByName("GTIFF")

    output_dst = driver_tif.Create( 'tmp/test_mosaic1.tif', 100, 100, 3, gdal.GDT_Byte)
    output_dst.GetRasterBand(1).Fill(255)
    output_dst = None

    output_dst = driver_tif.Create( 'tmp/test_mosaic2.tif', 100, 100, 3, gdal.GDT_Byte)
    output_dst.GetRasterBand(1).Fill(127)
    output_dst = None
    
    ds = gdal.Open('data/test_mosaic.vrt')
    # A simple Checksum() cannot detect if the fix works or not as
    # Checksum() reads line per line, and we must use IRasterIO() on multi-line request
    data = ds.GetRasterBand(1).ReadRaster(90,0,20,100)
    import struct
    got = struct.unpack('B' * 20*100, data)
    for i in range(100):
        if got[i*20 + 9 ] != 255:
            gdaltest.post_reason('at line %d, did not find 255' % i)
            return 'fail'
    ds = None
    
    driver_tif.Delete('tmp/test_mosaic1.tif')
    driver_tif.Delete('tmp/test_mosaic2.tif')

    return 'success'


###############################################################################
# Test complex source with complex data (#3977)

def vrt_read_4():

    try:
        import numpy as np
    except:
        return 'skip'

    data = np.zeros((1, 1), np.complex64)
    data[0, 0] = 1. + 3.j

    drv = gdal.GetDriverByName('GTiff')
    ds = drv.Create("/vsimem/test.tif", 1, 1, 1, gdal.GDT_CFloat32)
    ds.GetRasterBand(1).WriteArray(data)
    ds = None

    complex_xml = '''<VRTDataset rasterXSize="1" rasterYSize="1">
  <VRTRasterBand dataType="CFloat32" band="1">
    <ComplexSource>
      <SourceFilename relativeToVRT="1">/vsimem/test.tif</SourceFilename>
      <SourceBand>1</SourceBand>
      <ScaleOffset>3</ScaleOffset>
      <ScaleRatio>2</ScaleRatio>
    </ComplexSource>
  </VRTRasterBand>
</VRTDataset>
'''

    ds = gdal.Open(complex_xml)
    scaleddata = ds.GetRasterBand(1).ReadAsArray()
    ds = None

    gdal.Unlink("/vsimem/test.tif")

    if scaleddata[0, 0].real != 5.0 or scaleddata[0, 0].imag != 9.0:
        gdaltest.post_reason('did not get expected value')
        print('scaleddata[0, 0]: %f %f' % (scaleddata[0, 0].real, scaleddata[0, 0].imag))
        return 'fail'

    return 'success'

for item in init_list:
    ut = gdaltest.GDALTest( 'VRT', item[0], item[1], item[2] )
    if ut is None:
        print( 'VRT tests skipped' )
        sys.exit()
    gdaltest_list.append( (ut.testOpen, item[0]) )
    
gdaltest_list.append( vrt_read_1 )
gdaltest_list.append( vrt_read_2 )
gdaltest_list.append( vrt_read_3 )
gdaltest_list.append( vrt_read_4 )

if __name__ == '__main__':

    gdaltest.setup_run( 'vrt_read' )

    gdaltest.run_tests( gdaltest_list )

    gdaltest.summarize()

