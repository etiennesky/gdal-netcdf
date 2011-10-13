/******************************************************************************
 * $Id: netcdfdataset.h 23197 2011-10-07 00:34:31Z pds $
 *
 * Project:  netCDF read/write Driver
 * Purpose:  GDAL bindings over netCDF library.
 * Author:   Frank Warmerdam, warmerdam@pobox.com
 *
 ******************************************************************************
 * Copyright (c) 2004, Frank Warmerdam
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ****************************************************************************/

#ifndef _NETCDFDATASET_H_INCLUDED_
#define _NETCDFATASET_H_INCLUDED_

#include <float.h>
#include "gdal_pam.h"
#include "gdal_priv.h"
#include "gdal_frmts.h"
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "netcdf.h"


/************************************************************************/
/* ==================================================================== */
/*			     netCDFDataset				*/
/* ==================================================================== */
/************************************************************************/

#define MAX_STR_LEN            8192

/* Useful defines of CF-1 convention standard variables related to mapping
 * & projection - see http://cf-pcmdi.llnl.gov/ */
#define AEA                    "albers_conical_equal_area"
#define AE                     "azimuthal_equidistant"
#define CEA                    "cylindrical_equal_area"
#define LAEA                   "lambert_azimuthal_equal_area"
#define LCEA                   "lambert_cylindrical_equal_area"
#define L_C_CONIC              "lambert_conformal_conic"
#define TM                     "transverse_mercator"
#define GRD_MAPPING_NAME       "grid_mapping_name"
#define GRD_MAPPING            "grid_mapping"
#define COORDINATES            "coordinates"
#define LONLAT                 "lon lat"
#define LATITUDE_LONGITUDE     "latitude_longitude"
#define MERCATOR               "mercator"
#define ORTHOGRAPHIC           "orthographic"
#define POLAR_STEREO           "polar_stereographic"
#define STEREO                 "stereographic"

#define STD_PARALLEL           "standard_parallel"
#define STD_PARALLEL_1         "standard_parallel_1"
#define STD_PARALLEL_2         "standard_parallel_2"
#define CENTRAL_MERIDIAN       "central_meridian"
#define LONG_CENTRAL_MERIDIAN  "longitude_of_central_meridian"
#define LON_PROJ_ORIGIN        "longitude_of_projection_origin"
#define LAT_PROJ_ORIGIN        "latitude_of_projection_origin"
#define SCALE_FACTOR_ORIGIN    "scale_factor_at_projection_origin"
#define PROJ_X_ORIGIN          "projection_x_coordinate_origin"
#define PROJ_Y_ORIGIN          "projection_y_coordinate_origin"
#define EARTH_SHAPE            "GRIB_earth_shape"
#define EARTH_SHAPE_CODE       "GRIB_earth_shape_code"
//this is not CF but there are two possible translations SCALE_FACTOR_MERIDIAN and  SCALE_FACTOR_ORIGIN
#define SCALE_FACTOR           "scale_factor" 
#define SCALE_FACTOR_MERIDIAN  "scale_factor_at_central_meridian"
#define VERT_LONG_FROM_POLE    "straight_vertical_longitude_from_pole"
#define FALSE_EASTING          "false_easting"
#define FALSE_NORTHING         "false_northing"
#define EARTH_RADIUS           "earth_radius"
#define INVERSE_FLATTENING     "inverse_flattening"
#define LONG_PRIME_MERIDIAN    "longitude_of_prime_meridian"
#define SEMI_MAJOR_AXIS        "semi_major_axis"
#define SEMI_MINOR_AXIS        "semi_minor_axis"

#define STD_NAME               "standard_name"
#define LNG_NAME               "long_name"
#define UNITS                  "units"
#define AXIS                   "axis"
#define BOUNDS                 "bounds"
#define ORIG_AXIS              "original_units"

#define GDALNBDIM  2

/* netcdf file types, as in libcdi/cdo */
#define NCDF_FILETYPE_NONE            0   /* Not a netCDF file */
#define NCDF_FILETYPE_NC              1   /* File type netCDF                     */
#define NCDF_FILETYPE_NC2             2   /* File type netCDF version 2 (64-bit)  */
#define NCDF_FILETYPE_NC4             3   /* File type netCDF version 4           */
#define NCDF_FILETYPE_NC4C            4   /* File type netCDF version 4 (classic) - not used yet */
/* File type HDF5, not supported here (lack of netCDF-4 support or extension is not .nc or .nc4 */
#define NCDF_FILETYPE_HDF5            5   
#define NCDF_FILETYPE_UNKNOWN         10  /* Filetype not determined (yet) */

/* new defs */
#define NCDF_DIMNAME_X "x"
#define NCDF_DIMNAME_Y "y"
#define NCDF_DIMNAME_LON "lon"
#define NCDF_DIMNAME_LAT "lat"
#define NCDF_CONVENTIONS_CF "CF-1.5"
#define NCDF_GDAL GDALVersionInfo("--version")
#define NCDF_SPATIAL_REF "spatial_ref"

#define NCDF_


/* Following are a series of mappings from CF-1 convention parameters
 * for each projection, to the equivalent in OGC WKT used internally by
 * GDAL.
 * See: http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.5/apf.html
 */

/* A struct allowing us to map from GDAL projection attributes (OGC WKT),
   and NetCDF ones (CF-1) */
typedef struct {
    const char *NCDF_ATT;
    const char *GDAL_ATT; 
    // TODO: mappings may need default values, like scale factor?
    //double defval;
} oNetcdfSRS_PP;

//Grid mapping attributes
//
//earth_radius
//false_easting 	
//false_northing 	
//grid_mapping_name 	
//grid_north_pole_latitude
//grid_north_pole_longitude
//inverse_flattening
//latitude_of_projection_origin 
//longitude_of_central_meridian 
//longitude_of_prime_meridian
//longitude_of_projection_origin
//north_pole_grid_longitude 
//perspective_point_height	
//scale_factor_at_central_meridian 
//scale_factor_at_projection_origin 
//semi_major_axis
//semi_minor_axis
//standard_parallel 	
//straight_vertical_longitude_from_pole 	


// default mappings, for the generic case
/* These 'generic' mappings are based on what was previously in the  
   poNetCDFSRS struct. They will be used as a fallback in case none 
   of the others match (ie you are exporting a projection that has 
   no CF-1 equivalent). 
   They are not used for known CF-1 projections since there is not a 
   unique 2-way projection-independent 
   mapping between OGC WKT params and CF-1 ones: it varies per-projection. 
*/ 

static const oNetcdfSRS_PP poGenericMappings[] = {
    {SCALE_FACTOR, SRS_PP_SCALE_FACTOR },    
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1 }, 
    {STD_PARALLEL_2, SRS_PP_STANDARD_PARALLEL_2 }, 
    {LONG_CENTRAL_MERIDIAN, SRS_PP_CENTRAL_MERIDIAN }, 
    {LONG_CENTRAL_MERIDIAN, SRS_PP_LONGITUDE_OF_CENTER }, 
    {LON_PROJ_ORIGIN, SRS_PP_LONGITUDE_OF_ORIGIN },  
    //Multiple mappings to LAT_PROJ_ORIGIN 
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN },  
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_CENTER },  
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },   
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },        
    {NULL, NULL },
};

// Albers equal area 
//
// grid_mapping_name = albers_conical_equal_area
// WKT: Albers_Conic_Equal_Area
// ESPG:9822 
//
// Map parameters:
//
//    * standard_parallel - There may be 1 or 2 values.
//    * longitude_of_central_meridian
//    * latitude_of_projection_origin
//    * false_easting
//    * false_northing
//
static const oNetcdfSRS_PP poAEAMappings[] = {
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1},
    {STD_PARALLEL_2, SRS_PP_STANDARD_PARALLEL_2},
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_CENTER},
    {LONG_CENTRAL_MERIDIAN, SRS_PP_LONGITUDE_OF_CENTER},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Azimuthal equidistant
//
// grid_mapping_name = azimuthal_equidistant
// WKT: Azimuthal_Equidistant
//
// Map parameters:
//
//    * longitude_of_projection_origin
//    * latitude_of_projection_origin
//    * false_easting
//    * false_northing
//
static const oNetcdfSRS_PP poAEMappings[] = {
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_CENTER},
    {LON_PROJ_ORIGIN, SRS_PP_LONGITUDE_OF_CENTER},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Lambert azimuthal equal area
//
// grid_mapping_name = lambert_azimuthal_equal_area
// WKT: Lambert_Azimuthal_Equal_Area
//
// Map parameters:
//
//    * longitude_of_projection_origin
//    * latitude_of_projection_origin
//    * false_easting
//    * false_northing
//
static const oNetcdfSRS_PP poLAEAMappings[] = {
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_CENTER},
    {LON_PROJ_ORIGIN, SRS_PP_LONGITUDE_OF_CENTER},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Lambert conformal
//
// grid_mapping_name = lambert_conformal_conic
// WKT: Lambert_Conformal_Conic_1SP / Lambert_Conformal_Conic_2SP
//
// Map parameters:
//
//    * standard_parallel - There may be 1 or 2 values.
//    * longitude_of_central_meridian
//    * latitude_of_projection_origin
//    * false_easting
//    * false_northing
//
// See http://www.remotesensing.org/geotiff/proj_list/lambert_conic_conformal_1sp.html 

// Lambert conformal conic - 1SP
static const oNetcdfSRS_PP poLCC1SPMappings[] = {
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1},
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN},
    {LONG_CENTRAL_MERIDIAN, SRS_PP_CENTRAL_MERIDIAN},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Lambert conformal conic - 2SP
static const oNetcdfSRS_PP poLCC2SPMappings[] = {
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1},
    {STD_PARALLEL_2, SRS_PP_STANDARD_PARALLEL_2},
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN},
    {LONG_CENTRAL_MERIDIAN, SRS_PP_CENTRAL_MERIDIAN},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Lambert cylindrical equal area
//
// grid_mapping_name = lambert_cylindrical_equal_area
// WKT: Cylindrical_Equal_Area
// EPSG:9834 (Spherical) and EPSG:9835 
//
// Map parameters:
//
//    * longitude_of_central_meridian
//    * either standard_parallel or scale_factor_at_projection_origin
//    * false_easting
//    * false_northing
//
// NB: CF-1 specifies a 'scale_factor_at_projection' alternative  
//  to std_parallel ... but no reference to this in EPSG/remotesensing.org 
//  ignore for now. 
//
static const oNetcdfSRS_PP poLCEAMappings[] = {
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1},
    {LONG_CENTRAL_MERIDIAN, SRS_PP_CENTRAL_MERIDIAN},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Latitude-Longitude
//
// grid_mapping_name = latitude_longitude
//
// Map parameters:
//
//    * None
//
// NB: handled as a special case - !isProjected()


// Mercator
//
// grid_mapping_name = mercator
// WKT: Mercator_1SP / Mercator_2SP
//
// Map parameters:
//
//    * longitude_of_projection_origin
//    * either standard_parallel or scale_factor_at_projection_origin
//    * false_easting
//    * false_northing

// Mercator 1 Standard Parallel (EPSG:9804) 
static const oNetcdfSRS_PP poM1SPMappings[] = {
    {LON_PROJ_ORIGIN, SRS_PP_CENTRAL_MERIDIAN},
    //LAT_PROJ_ORIGIN is always equator (0) in CF-1 
    {SCALE_FACTOR_ORIGIN, SRS_PP_SCALE_FACTOR},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Mercator 2 Standard Parallel
static const oNetcdfSRS_PP poM2SPMappings[] = {
    {LON_PROJ_ORIGIN, SRS_PP_CENTRAL_MERIDIAN},
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1},
    //From best understanding of this projection, only  
 	// actually specify one SP - it is the same N/S of equator. 
    //{STD_PARALLEL_2, SRS_PP_LATITUDE_OF_ORIGIN}, 
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Orthographic
// grid_mapping_name = orthographic
// WKT: Orthographic
//
// Map parameters:
//
//    * longitude_of_projection_origin
//    * latitude_of_projection_origin
//    * false_easting
//    * false_northing
//
static const oNetcdfSRS_PP poOrthoMappings[] = {
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN},
    {LON_PROJ_ORIGIN, SRS_PP_CENTRAL_MERIDIAN},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 }; 

// Polar stereographic
//
// grid_mapping_name = polar_stereographic
// WKT: Polar_Stereographic
//
// Map parameters:
//
//    * straight_vertical_longitude_from_pole
//    * latitude_of_projection_origin - Either +90. or -90.
//    * Either standard_parallel or scale_factor_at_projection_origin
//    * false_easting
//    * false_northing

/* TODO: CF-1 says 'standard_parallel' may replace scale factor
   TODO: In CF-1, LAT_PROJ_ORIGIN is either +90 or -90. Not sure
   how this maps to OGC stuff?
   eg do we have to check LAT_PROJ_ORIGIN is +90 or -90?
   Or in this case, is "std_parallel" in CF-1 really the same as
   LAT_OF_ORIGIN in WKT, and we fill in LAT_OF_ORIGIN with +90 or 
   -90? (LAT_PROJ_ORIGIN set to 90 if STD_PARALLEL > 0, else -90)
*/

static const oNetcdfSRS_PP poPSmappings[] = {
    //CF-1 says 'either' of the following attribs can be used to
    //define the projection, for now read and set both
    {STD_PARALLEL_1, SRS_PP_LATITUDE_OF_ORIGIN},
    {SCALE_FACTOR_ORIGIN, SRS_PP_SCALE_FACTOR},  
    {VERT_LONG_FROM_POLE, SRS_PP_CENTRAL_MERIDIAN},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
};

// Rotated Pole
//
// grid_mapping_name = rotated_latitude_longitude
// WKT: N/A
//
// Map parameters:
//
//    * grid_north_pole_latitude
//    * grid_north_pole_longitude
//    * north_pole_grid_longitude - This parameter is optional (default is 0.).

/* TODO: No GDAL equivalent of rotated pole? Doesn't seem to have an EPSG
   code or WKT ... so unless some advanced proj4 features can be used 
   seems to rule out.
   see GDAL bug #4285 for a possible fix or workaround
*/

// Stereographic
//
// grid_mapping_name = stereographic
// WKT: Stereographic (and/or Oblique_Stereographic??)
//
// Map parameters:
//
//    * longitude_of_projection_origin
//    * latitude_of_projection_origin
//    * scale_factor_at_projection_origin
//    * false_easting
//    * false_northing
//
// NB: see bug#4267 Stereographic vs. Oblique_Stereographic
//
static const oNetcdfSRS_PP poStMappings[] = {
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN},
    {LON_PROJ_ORIGIN, SRS_PP_CENTRAL_MERIDIAN},
    {SCALE_FACTOR_ORIGIN, SRS_PP_SCALE_FACTOR},  
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
  };

// Transverse Mercator
//
// grid_mapping_name = transverse_mercator
// WKT: Transverse_Mercator
//
// Map parameters:
//
//    * scale_factor_at_central_meridian
//    * longitude_of_central_meridian
//    * latitude_of_projection_origin
//    * false_easting
//    * false_northing
//
static const oNetcdfSRS_PP poTMMappings[] = {
    {SCALE_FACTOR_MERIDIAN, SRS_PP_SCALE_FACTOR},  
    {LONG_CENTRAL_MERIDIAN, SRS_PP_CENTRAL_MERIDIAN},
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
  };

// Vertical perspective
//
// grid_mapping_name = vertical_perspective
// WKT: ???
//
// Map parameters:
//
//    * latitude_of_projection_origin
//    * longitude_of_projection_origin
//    * perspective_point_height
//    * false_easting
//    * false_northing
//
// TODO: see how to map this to OGR


/* Mappings for various projections, including netcdf and GDAL projection names 
   and corresponding oNetcdfSRS_PP mapping struct. 
   A NULL mappings value means that the projection is not included in the CF
   standard and the generic mapping (poGenericMappings) will be used. */
typedef struct {
    const char *NCDF_SRS;
    const char *GDAL_SRS; 
    const oNetcdfSRS_PP* mappings;
} oNetcdfSRS_PT;

static const oNetcdfSRS_PT poNetcdfSRS_PT[] = {
    {AEA, SRS_PT_ALBERS_CONIC_EQUAL_AREA, poAEAMappings },
    {AE, SRS_PT_AZIMUTHAL_EQUIDISTANT, poAEMappings },
    {"cassini_soldner", SRS_PT_CASSINI_SOLDNER, NULL },
    {LCEA, SRS_PT_CYLINDRICAL_EQUAL_AREA, poLCEAMappings },
    {"eckert_iv", SRS_PT_ECKERT_IV, NULL },      
    {"eckert_vi", SRS_PT_ECKERT_VI, NULL },  
    {"equidistant_conic", SRS_PT_EQUIDISTANT_CONIC, NULL },
    {"equirectangular", SRS_PT_EQUIRECTANGULAR, NULL },
    {"gall_stereographic", SRS_PT_GALL_STEREOGRAPHIC, NULL },
    {"geostationary_satellite", SRS_PT_GEOSTATIONARY_SATELLITE, NULL },
    {"goode_homolosine", SRS_PT_GOODE_HOMOLOSINE, NULL },
    {"gnomonic", SRS_PT_GNOMONIC, NULL },
    {"hotine_oblique_mercator", SRS_PT_HOTINE_OBLIQUE_MERCATOR, NULL },
    {"hotine_oblique_mercator_2P", 
     SRS_PT_HOTINE_OBLIQUE_MERCATOR_TWO_POINT_NATURAL_ORIGIN, NULL },
    {"laborde_oblique_mercator", SRS_PT_LABORDE_OBLIQUE_MERCATOR, NULL },
    {L_C_CONIC, SRS_PT_LAMBERT_CONFORMAL_CONIC_1SP, poLCC1SPMappings },
    {L_C_CONIC, SRS_PT_LAMBERT_CONFORMAL_CONIC_2SP, poLCC2SPMappings },
    {LAEA, SRS_PT_LAMBERT_AZIMUTHAL_EQUAL_AREA, poLAEAMappings },
    {MERCATOR, SRS_PT_MERCATOR_1SP, poM1SPMappings },
    {MERCATOR, SRS_PT_MERCATOR_2SP, poM2SPMappings },
    {"miller_cylindrical", SRS_PT_MILLER_CYLINDRICAL, NULL },
    {"mollweide", SRS_PT_MOLLWEIDE, NULL },
    {"new_zealand_map_grid", SRS_PT_NEW_ZEALAND_MAP_GRID, NULL },
    /* for now map to STEREO, see bug #4267 */
    /* {"oblique_stereographic", SRS_PT_OBLIQUE_STEREOGRAPHIC, NULL },  */
    {STEREO, SRS_PT_OBLIQUE_STEREOGRAPHIC, poStMappings }, 
    {"orthographic", SRS_PT_ORTHOGRAPHIC, poOrthoMappings },
    {POLAR_STEREO, SRS_PT_POLAR_STEREOGRAPHIC, poPSmappings },
    {"polyconic", SRS_PT_POLYCONIC, NULL },
    {"robinson", SRS_PT_ROBINSON, NULL }, 
    {"sinusoidal", SRS_PT_SINUSOIDAL, NULL },  
    {STEREO, SRS_PT_STEREOGRAPHIC, poStMappings },
    {"swiss_oblique_cylindrical", SRS_PT_SWISS_OBLIQUE_CYLINDRICAL, NULL },
    {TM, SRS_PT_TRANSVERSE_MERCATOR, poTMMappings },
    {"TM_south_oriented", SRS_PT_TRANSVERSE_MERCATOR_SOUTH_ORIENTED, NULL },
    {NULL, NULL, NULL },
};


class netCDFRasterBand;

class netCDFDataset : public GDALPamDataset
{
    CPLString    osSubdatasetName;
    int          bTreatAsSubdataset;

    double      adfGeoTransform[6];
    char        **papszSubDatasets;
    char        **papszGeolocation;
    CPLString    osFilename;
    int          *panBandDimPos;         // X, Y, Z postion in array
    int          *panBandZLev;
    char         *pszProjection;
    int          bGotGeoTransform;
    double       rint( double );

    double       FetchCopyParm( const char *pszGridMappingValue, 
                                const char *pszParm, double dfDefault );

    char **      FetchStandardParallels( const char *pszGridMappingValue );

    static int IdentifyFileType( GDALOpenInfo *, bool );

  public:
    int           cdfid;
    char         **papszMetadata;
    char          papszDimName[NC_MAX_NAME][1024];
    int          *paDimIds;
    size_t        xdim, ydim;
    int           nDimXid, nDimYid;
    bool          bBottomUp;
    int           nFileType;

    netCDFDataset( );
    ~netCDFDataset( );
    
    static int Identify( GDALOpenInfo * );
    static GDALDataset *Open( GDALOpenInfo * );

    CPLErr      SafeStrcat(char**, char*, size_t*);
    CPLErr      ReadAttributes( int, int );

    CPLErr 	GetGeoTransform( double * );    

    const char * GetProjectionRef();

    char ** GetMetadata( const char * );

    void  CreateSubDatasetList( );

    void  SetProjection( int );

};

#endif
