#!/usr/bin/env python
###############################################################################
# $Id$
#
# Project:  GDAL/OGR Test Suite
# Purpose:  Test COMPD_CS support.
# Author:   Frank Warmerdam <warmerdam@pobox.com>
# 
###############################################################################
# Copyright (c) 2010, Frank Warmerdam <warmerdam@pobox.com>
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
try:
    from osgeo import osr
except:
    import osr

example_compd_wkt = 'COMPD_CS["OSGB36 / British National Grid + ODN",PROJCS["OSGB 1936 / British National Grid",GEOGCS["OSGB 1936",DATUM["OSGB_1936",SPHEROID["Airy 1830",6377563.396,299.3249646,AUTHORITY["EPSG",7001]],TOWGS84[375,-111,431,0,0,0,0],AUTHORITY["EPSG",6277]],PRIMEM["Greenwich",0,AUTHORITY["EPSG",8901]],UNIT["DMSH",0.0174532925199433,AUTHORITY["EPSG",9108]],AXIS["Lat",NORTH],AXIS["Long",EAST],AUTHORITY["EPSG",4277]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",49],PARAMETER["central_meridian",-2],PARAMETER["scale_factor",0.999601272],PARAMETER["false_easting",400000],PARAMETER["false_northing",-100000],UNIT["metre_1",1,AUTHORITY["EPSG",9001]],AXIS["E",EAST],AXIS["N",NORTH],AUTHORITY["EPSG",27700]],VERT_CS["Newlyn",VERT_DATUM["Ordnance Datum Newlyn",2005,AUTHORITY["EPSG",5101]],UNIT["metre_2",1,AUTHORITY["EPSG",9001]],AXIS["Up",UP],AUTHORITY["EPSG",5701]],AUTHORITY["EPSG",7405]]'


###############################################################################
# Test parsing and a few operations on a compound coordinate system.

def osr_compd_1():

    srs = osr.SpatialReference()
    srs.ImportFromWkt( example_compd_wkt )

    if not srs.IsProjected():
        gdaltest.post_reason( 'Projected COMPD_CS not recognised as projected.')
        return 'fail'

    if srs.IsGeographic():
        gdaltest.post_reason( 'projected COMPD_CS misrecognised as geographic.')
        return 'fail'

    if srs.IsLocal():
        gdaltest.post_reason( 'projected COMPD_CS misrecognised as local.')
        return 'fail'

    expected_proj4 = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.999601272 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=375,-111,431,0,0,0,0 +units=m +vunits=m +no_defs '
    got_proj4 = srs.ExportToProj4()

    if expected_proj4 != got_proj4:
        print( 'Got:      %s' % got_proj4 )
        print( 'Expected: %s' % expected_proj4 )
        gdaltest.post_reason( 'did not get expected proj.4 translation of compd_cs' )
        return 'fail'
    
    if srs.GetLinearUnitsName() != 'metre_1':
        gdaltest.post_reason( 'Did not get expected linear units.' )
        return 'false'

    return 'success'

###############################################################################
# Test SetFromUserInput()

def osr_compd_2():

    srs = osr.SpatialReference()
    srs.SetFromUserInput( example_compd_wkt )

    if not srs.IsProjected():
        gdaltest.post_reason( 'Projected COMPD_CS not recognised as projected.')
        return 'fail'

    return 'success'

###############################################################################
# Test expansion of compound coordinate systems from EPSG definition.

def osr_compd_3():

    srs = osr.SpatialReference()
    srs.ImportFromEPSG( 7401 )

    exp_wkt = """COMPD_CS["NTF (Paris) / France II + NGF Lallemand",
    PROJCS["NTF (Paris) / France II (deprecated)",
        GEOGCS["NTF (Paris)",
            DATUM["Nouvelle_Triangulation_Francaise_Paris",
                SPHEROID["Clarke 1880 (IGN)",6378249.2,293.4660212936265,
                    AUTHORITY["EPSG","7011"]],
                TOWGS84[-168,-60,320,0,0,0,0],
                AUTHORITY["EPSG","6807"]],
            PRIMEM["Paris",2.33722917,
                AUTHORITY["EPSG","8903"]],
            UNIT["grad",0.01570796326794897,
                AUTHORITY["EPSG","9105"]],
            AUTHORITY["EPSG","4807"]],
        UNIT["metre",1,
            AUTHORITY["EPSG","9001"]],
        PROJECTION["Lambert_Conformal_Conic_1SP"],
        PARAMETER["latitude_of_origin",52],
        PARAMETER["central_meridian",0],
        PARAMETER["scale_factor",0.99987742],
        PARAMETER["false_easting",600000],
        PARAMETER["false_northing",2200000],
        AUTHORITY["EPSG","27582"],
        AXIS["X",EAST],
        AXIS["Y",NORTH]],
    VERT_CS["NGF Lallemand height",
        VERT_DATUM["Nivellement General de la France - Lallemand",2005,
            AUTHORITY["EPSG","5118"]],
        UNIT["metre",1,
            AUTHORITY["EPSG","9001"]],
        AXIS["Up",UP],
        AUTHORITY["EPSG","5719"]],
    AUTHORITY["EPSG","7401"]]"""
    wkt = srs.ExportToPrettyWkt() 
    if gdaltest.equal_srs_from_wkt( exp_wkt, wkt ) == 0:
        gdaltest.post_reason( 'did not get expected compound cs for EPSG:7401')
        return 'fail'
    
    return 'success'

###############################################################################
# Test expansion of GCS+VERTCS compound coordinate system.

def osr_compd_4():

    srs = osr.SpatialReference()
    srs.ImportFromEPSG( 7400 )

    exp_wkt = """COMPD_CS["NTF (Paris) + NGF IGN69 height",
    GEOGCS["NTF (Paris)",
        DATUM["Nouvelle_Triangulation_Francaise_Paris",
            SPHEROID["Clarke 1880 (IGN)",6378249.2,293.4660212936265,
                AUTHORITY["EPSG","7011"]],
            TOWGS84[-168,-60,320,0,0,0,0],
            AUTHORITY["EPSG","6807"]],
        PRIMEM["Paris",2.33722917,
            AUTHORITY["EPSG","8903"]],
        UNIT["grad",0.01570796326794897,
            AUTHORITY["EPSG","9105"]],
        AUTHORITY["EPSG","4807"]],
    VERT_CS["NGF IGN69 height",
        VERT_DATUM["Nivellement General de la France - IGN69",2005,
            AUTHORITY["EPSG","5119"]],
        UNIT["metre",1,
            AUTHORITY["EPSG","9001"]],
        AXIS["Up",UP],
        AUTHORITY["EPSG","5720"]],
    AUTHORITY["EPSG","7400"]]"""
    wkt = srs.ExportToPrettyWkt() 

    if gdaltest.equal_srs_from_wkt( exp_wkt, wkt ) == 0:
        gdaltest.post_reason( 'did not get expected compound cs for EPSG:7400')
        return 'fail'
    
    return 'success'

###############################################################################
# Test that compound coordinate systems with grid shift files are
# expanded properly and converted to PROJ.4 format with the grids.

def osr_compd_5():

    srs = osr.SpatialReference()
    srs.SetFromUserInput( 'EPSG:26911+5703' )

    exp_wkt = """COMPD_CS["NAD83 / UTM zone 11N + VERT_CS",
    PROJCS["NAD83 / UTM zone 11N",
        GEOGCS["NAD83",
            DATUM["North_American_Datum_1983",
                SPHEROID["GRS 1980",6378137,298.257222101,
                    AUTHORITY["EPSG","7019"]],
                TOWGS84[0,0,0,0,0,0,0],
                AUTHORITY["EPSG","6269"]],
            PRIMEM["Greenwich",0,
                AUTHORITY["EPSG","8901"]],
            UNIT["degree",0.0174532925199433,
                AUTHORITY["EPSG","9122"]],
            AUTHORITY["EPSG","4269"]],
        UNIT["metre",1,
            AUTHORITY["EPSG","9001"]],
        PROJECTION["Transverse_Mercator"],
        PARAMETER["latitude_of_origin",0],
        PARAMETER["central_meridian",-117],
        PARAMETER["scale_factor",0.9996],
        PARAMETER["false_easting",500000],
        PARAMETER["false_northing",0],
        AUTHORITY["EPSG","26911"],
        AXIS["Easting",EAST],
        AXIS["Northing",NORTH]],
    VERT_CS["NAVD88 height",
        VERT_DATUM["North American Vertical Datum 1988",2005,
            AUTHORITY["EPSG","5103"],
            EXTENSION["PROJ4_GRIDS","g2003conus.gtx,g2003alaska.gtx,g2003h01.gtx,g2003p01.gtx"]],
        UNIT["metre",1,
            AUTHORITY["EPSG","9001"]],
        AXIS["Up",UP],
        AUTHORITY["EPSG","5703"]]]"""
    wkt = srs.ExportToPrettyWkt() 

    if gdaltest.equal_srs_from_wkt( exp_wkt, wkt ) == 0:
        return 'fail'

    if wkt.find('g2003conus.gtx') == -1:
        gdaltest.post_reason( 'Did not get PROJ4_GRIDS EXTENSION node' )
        return 'fail'

    exp_proj4 = '+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +geoidgrids=g2003conus.gtx,g2003alaska.gtx,g2003h01.gtx,g2003p01.gtx +vunits=m +no_defs '
    proj4 = srs.ExportToProj4()
    if proj4 != exp_proj4:
        gdaltest.post_reason( 'Did not get expected proj.4 string, got:' + proj4 )
        return 'fail'
    
    return 'success'

###############################################################################
# Test conversion from PROJ.4 to WKT including vertical units.

def osr_compd_6():

    srs = osr.SpatialReference()
    srs.SetFromUserInput( '+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +geoidgrids=g2003conus.gtx,g2003alaska.gtx,g2003h01.gtx,g2003p01.gtx +vunits=us-ft +no_defs ' )

    exp_wkt = """COMPD_CS["UTM Zone 11, Northern Hemisphere + Unnamed Vertical Datum",
    PROJCS["UTM Zone 11, Northern Hemisphere",
        GEOGCS["GRS 1980(IUGG, 1980)",
            DATUM["unknown",
                SPHEROID["GRS80",6378137,298.257222101],
                TOWGS84[0,0,0,0,0,0,0]],
            PRIMEM["Greenwich",0],
            UNIT["degree",0.0174532925199433]],
        PROJECTION["Transverse_Mercator"],
        PARAMETER["latitude_of_origin",0],
        PARAMETER["central_meridian",-117],
        PARAMETER["scale_factor",0.9996],
        PARAMETER["false_easting",500000],
        PARAMETER["false_northing",0],
        UNIT["Meter",1]],
    VERT_CS["Unnamed",
        VERT_DATUM["Unnamed",2005,
            EXTENSION["PROJ4_GRIDS","g2003conus.gtx,g2003alaska.gtx,g2003h01.gtx,g2003p01.gtx"]],
        AXIS["Up",UP],
        UNIT["Foot_US",0.3048006096012192]]]"""
            
    wkt = srs.ExportToPrettyWkt() 

    if gdaltest.equal_srs_from_wkt( exp_wkt, wkt ) == 0:
        return 'fail'

    if wkt.find('g2003conus.gtx') == -1:
        gdaltest.post_reason( 'Did not get PROJ4_GRIDS EXTENSION node' )
        return 'fail'

    exp_proj4 = '+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +geoidgrids=g2003conus.gtx,g2003alaska.gtx,g2003h01.gtx,g2003p01.gtx +vunits=us-ft +no_defs '
    proj4 = srs.ExportToProj4()
    if proj4 != exp_proj4:
        gdaltest.post_reason( 'Did not get expected proj.4 string, got:' + proj4 )
        return 'fail'
    
    return 'success'

gdaltest_list = [ 
    osr_compd_1,
    osr_compd_2,
    osr_compd_3,
    osr_compd_4,
    osr_compd_5,
    osr_compd_6,
    None ]

if __name__ == '__main__':

    gdaltest.setup_run( 'osr_compd' )

    gdaltest.run_tests( gdaltest_list )

    gdaltest.summarize()

