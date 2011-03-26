#!/usr/bin/env python
###############################################################################
# $Id$
#
# Project:  GDAL/OGR Test Suite
# Purpose:  Test WEBP driver
# Author:   Even Rouault, <even dot rouault at mines dash paris dot org>
# 
###############################################################################
# Copyright (c) 2011, Even Rouault
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

sys.path.append( '../pymod' )

import gdaltest

###############################################################################
# Test if WEBP driver is present

def webp_1():

    try:
        gdaltest.webp_drv = gdal.GetDriverByName( 'WEBP' )
    except:
        gdaltest.webp_drv = None
        return 'skip'

    return 'success'

###############################################################################
# Open() test

def webp_2():

    if gdaltest.webp_drv is None:
        return 'skip'

    tst = gdaltest.GDALTest( 'WEBP', 'rgbsmall.webp', 1, 21498 )
    return tst.testOpen()

###############################################################################
# CreateCopy() test

def webp_3():

    if gdaltest.webp_drv is None:
        return 'skip'

    tst = gdaltest.GDALTest( 'WEBP', 'rgbsmall.tif', 1, 21498, options = ['QUALITY=80'] )
    return tst.testCreateCopy( vsimem = 1, check_minmax = False )

gdaltest_list = [
    webp_1,
    webp_2,
    webp_3 ]

if __name__ == '__main__':

    gdaltest.setup_run( 'webp' )

    gdaltest.run_tests( gdaltest_list )

    gdaltest.summarize()

