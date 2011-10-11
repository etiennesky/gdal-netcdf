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
#define LCEA                   "lambert_cylindrical_equal_area"
#define L_C_CONIC              "lambert_conformal_conic"
#define TM                     "transverse_mercator"
#define LAEA                   "lambert_azimuthal_equal_area"
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
#define SCALE_FACTOR           "scale_factor_at_central_meridian" //this has to go
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

/*original mappings (both SRS_PT and SRS_PP) */

typedef struct {
    const char *netCDFSRS;
    const char *SRS; 
} oNetcdfSRS;

static const oNetcdfSRS poNetcdfSRS[] = {
    {"albers_conical_equal_area", SRS_PT_ALBERS_CONIC_EQUAL_AREA },
    {"azimuthal_equidistant", SRS_PT_AZIMUTHAL_EQUIDISTANT },
    {"cassini_soldner", SRS_PT_CASSINI_SOLDNER },
    {"lambert_cylindrical_equal_area", SRS_PT_CYLINDRICAL_EQUAL_AREA },
    {"eckert_iv", SRS_PT_ECKERT_IV },      
    {"eckert_vi", SRS_PT_ECKERT_VI },  
    {"equidistant_conic", SRS_PT_EQUIDISTANT_CONIC },
    {"equirectangular", SRS_PT_EQUIRECTANGULAR },
    {"gall_stereographic", SRS_PT_GALL_STEREOGRAPHIC },
    {"geostationary_satellite", SRS_PT_GEOSTATIONARY_SATELLITE },
    {"goode_homolosine", SRS_PT_GOODE_HOMOLOSINE },
    {"gnomonic", SRS_PT_GNOMONIC },
    {"hotine_oblique_mercator", SRS_PT_HOTINE_OBLIQUE_MERCATOR},
    {"hotine_oblique_mercator_2P", 
     SRS_PT_HOTINE_OBLIQUE_MERCATOR_TWO_POINT_NATURAL_ORIGIN},
    {"laborde_oblique_mercator", SRS_PT_LABORDE_OBLIQUE_MERCATOR },
    {"lambert_conformal_conic1", SRS_PT_LAMBERT_CONFORMAL_CONIC_1SP },
    {"lambert_conformal_conic", SRS_PT_LAMBERT_CONFORMAL_CONIC_2SP },
    {"lambert_azimuthal_equal_area", SRS_PT_LAMBERT_AZIMUTHAL_EQUAL_AREA },
    {"mercator_1sp", SRS_PT_MERCATOR_1SP },
    {"mercator_2sp", SRS_PT_MERCATOR_2SP },
    {"miller_cylindrical", SRS_PT_MILLER_CYLINDRICAL },
    {"mollweide", SRS_PT_MOLLWEIDE },
    {"new_zealand_map_grid", SRS_PT_NEW_ZEALAND_MAP_GRID },
    {"oblique_stereographic", SRS_PT_OBLIQUE_STEREOGRAPHIC }, 
    {"orthographic", SRS_PT_ORTHOGRAPHIC },
    {"polar_stereographic", SRS_PT_POLAR_STEREOGRAPHIC },
    {"polyconic", SRS_PT_POLYCONIC },
    {"robinson", SRS_PT_ROBINSON }, 
    {"sinusoidal", SRS_PT_SINUSOIDAL },  
    {"stereographic", SRS_PT_STEREOGRAPHIC },
    {"swiss_oblique_cylindrical", SRS_PT_SWISS_OBLIQUE_CYLINDRICAL},
    {"transverse_mercator", SRS_PT_TRANSVERSE_MERCATOR },
    {"TM_south_oriented", SRS_PT_TRANSVERSE_MERCATOR_SOUTH_ORIENTED },

    {LONG_CENTRAL_MERIDIAN, SRS_PP_CENTRAL_MERIDIAN },
    {SCALE_FACTOR, SRS_PP_SCALE_FACTOR },   
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1 },
    {STD_PARALLEL_2, SRS_PP_STANDARD_PARALLEL_2 },
    {"longitude_of_central_meridian", SRS_PP_LONGITUDE_OF_CENTER },
    {"longitude_of_projection_origin", SRS_PP_LONGITUDE_OF_ORIGIN }, 
    {"latitude_of_projection_origin", SRS_PP_LATITUDE_OF_ORIGIN }, 
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },       
    {NULL, NULL },
 };


/*new SRS_PT mappings */
typedef struct {
    const char *NCDF_SRS;
    const char *GDAL_SRS; 
} oNetcdf_SRS_PT;

static const oNetcdf_SRS_PT poNetcdf_SRS_PT[] = {
    {"albers_conical_equal_area", SRS_PT_ALBERS_CONIC_EQUAL_AREA },
    {"azimuthal_equidistant", SRS_PT_AZIMUTHAL_EQUIDISTANT },
    {"cassini_soldner", SRS_PT_CASSINI_SOLDNER },
    {"lambert_cylindrical_equal_area", SRS_PT_CYLINDRICAL_EQUAL_AREA },
    {"eckert_iv", SRS_PT_ECKERT_IV },      
    {"eckert_vi", SRS_PT_ECKERT_VI },  
    {"equidistant_conic", SRS_PT_EQUIDISTANT_CONIC },
    {"equirectangular", SRS_PT_EQUIRECTANGULAR },
    {"gall_stereographic", SRS_PT_GALL_STEREOGRAPHIC },
    {"geostationary_satellite", SRS_PT_GEOSTATIONARY_SATELLITE },
    {"goode_homolosine", SRS_PT_GOODE_HOMOLOSINE },
    {"gnomonic", SRS_PT_GNOMONIC },
    {"hotine_oblique_mercator", SRS_PT_HOTINE_OBLIQUE_MERCATOR},
    {"hotine_oblique_mercator_2P", 
     SRS_PT_HOTINE_OBLIQUE_MERCATOR_TWO_POINT_NATURAL_ORIGIN},
    {"laborde_oblique_mercator", SRS_PT_LABORDE_OBLIQUE_MERCATOR },
    /* {"lambert_conformal_conic1", SRS_PT_LAMBERT_CONFORMAL_CONIC_1SP }, */
    {"lambert_conformal_conic", SRS_PT_LAMBERT_CONFORMAL_CONIC_1SP },
    {"lambert_conformal_conic", SRS_PT_LAMBERT_CONFORMAL_CONIC_2SP },
    {"lambert_azimuthal_equal_area", SRS_PT_LAMBERT_AZIMUTHAL_EQUAL_AREA },
    /* {"mercator_1sp", SRS_PT_MERCATOR_1SP }, */
    /* {"mercator_2sp", SRS_PT_MERCATOR_2SP }, */
    {"mercator", SRS_PT_MERCATOR_1SP },
    {"mercator", SRS_PT_MERCATOR_2SP },
    {"miller_cylindrical", SRS_PT_MILLER_CYLINDRICAL },
    {"mollweide", SRS_PT_MOLLWEIDE },
    {"new_zealand_map_grid", SRS_PT_NEW_ZEALAND_MAP_GRID },
    //PDS: GDAL seems to currently treat oblique_stereographic same as stereographic
    //{"oblique_stereographic", SRS_PT_OBLIQUE_STEREOGRAPHIC }, 
    {STEREO, SRS_PT_OBLIQUE_STEREOGRAPHIC }, 
    {"orthographic", SRS_PT_ORTHOGRAPHIC },
    {POLAR_STEREO, SRS_PT_POLAR_STEREOGRAPHIC },
    {"polyconic", SRS_PT_POLYCONIC },
    {"robinson", SRS_PT_ROBINSON }, 
    {"sinusoidal", SRS_PT_SINUSOIDAL },  
    {STEREO, SRS_PT_STEREOGRAPHIC },
    {"swiss_oblique_cylindrical", SRS_PT_SWISS_OBLIQUE_CYLINDRICAL},
    {"transverse_mercator", SRS_PT_TRANSVERSE_MERCATOR },
    {"TM_south_oriented", SRS_PT_TRANSVERSE_MERCATOR_SOUTH_ORIENTED },
    {NULL, NULL },
};

/*new SRS_PP mappings */
typedef struct {
    const char *NCDF_SRS;
    const char *GDAL_SRS; 
    double defval;
} oNetcdf_SRS_PP;

static const oNetcdf_SRS_PP poNetcdf_SRS_PP[] = {
    {SCALE_FACTOR, SRS_PP_SCALE_FACTOR },   
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1 },
    {STD_PARALLEL_2, SRS_PP_STANDARD_PARALLEL_2 },
    {LONG_CENTRAL_MERIDIAN, SRS_PP_CENTRAL_MERIDIAN },
    //Another test to get "Longitude origin" to override "Longitude_center"
    {LON_PROJ_ORIGIN, SRS_PP_LONGITUDE_OF_CENTER }, 
    //{LONG_CENTRAL_MERIDIAN, SRS_PP_LONGITUDE_OF_CENTER },
    {LON_PROJ_ORIGIN, SRS_PP_LONGITUDE_OF_ORIGIN }, 
    //PDS Test change: adding extra lat_of_proj_origin mapping to NetCDF
    //2nd of these is what's needed for AEA and several others ...
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN }, 
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_CENTER }, 
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },       
     {NULL, NULL },
};
/*original SRS_PP mappings */
static const oNetcdf_SRS_PP poNetcdf_SRS_PP1[] = {
    {LONG_CENTRAL_MERIDIAN, SRS_PP_CENTRAL_MERIDIAN },
    {SCALE_FACTOR, SRS_PP_SCALE_FACTOR },   
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1 },
    {STD_PARALLEL_2, SRS_PP_STANDARD_PARALLEL_2 },
    {"longitude_of_central_meridian", SRS_PP_LONGITUDE_OF_CENTER },
    {"longitude_of_projection_origin", SRS_PP_LONGITUDE_OF_ORIGIN }, 
    {"latitude_of_projection_origin", SRS_PP_LATITUDE_OF_ORIGIN }, 
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },       
    {NULL, NULL },
};

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

//Albers equal area 
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
static const oNetcdfSRS_PP poAZEQMappings[] = {
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_CENTER},
    {LON_PROJ_ORIGIN, SRS_PP_LONGITUDE_OF_CENTER},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Lambert azimuthal equal area
static const oNetcdfSRS_PP poLAZEQMappings[] = {
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_CENTER},
    {LON_PROJ_ORIGIN, SRS_PP_LONGITUDE_OF_CENTER},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Lambert conformal conic - 1SP
static const oNetcdfSRS_PP poLC_1SPMappings[] = {
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1},
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN},
    {LONG_CENTRAL_MERIDIAN, SRS_PP_CENTRAL_MERIDIAN},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Lambert conformal conic - 2SP
static const oNetcdfSRS_PP poLC_2SPMappings[] = {
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1},
    {STD_PARALLEL_2, SRS_PP_STANDARD_PARALLEL_2},
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN},
    {LONG_CENTRAL_MERIDIAN, SRS_PP_CENTRAL_MERIDIAN},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Lambert cylindrical equal area
// Also a "scale_factor_at_projection' CF-1 alternative to std_parallel
//  ... but not in OGC WKT. Perhaps need translation formula on import?
static const oNetcdfSRS_PP poLCEAMappings[] = {
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1},
    {LONG_CENTRAL_MERIDIAN, SRS_PP_CENTRAL_MERIDIAN},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Mercator 1 Standard Parallel
static const oNetcdfSRS_PP poM_1SPMappings[] = {
    {LON_PROJ_ORIGIN, SRS_PP_CENTRAL_MERIDIAN},
    {SCALE_FACTOR_ORIGIN, SRS_PP_SCALE_FACTOR},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Mercator 2 Standard Parallel
static const oNetcdfSRS_PP poM_2SPMappings[] = {
    {LON_PROJ_ORIGIN, SRS_PP_CENTRAL_MERIDIAN},
    {STD_PARALLEL_1, SRS_PP_STANDARD_PARALLEL_1},
    // ? Not entirely sure about the std_parallel_2, needs testing ...
    {STD_PARALLEL_2, SRS_PP_LATITUDE_OF_ORIGIN},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 };

// Orthographic
static const oNetcdfSRS_PP poOrthoMappings[] = {
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN},
    {LON_PROJ_ORIGIN, SRS_PP_CENTRAL_MERIDIAN},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
 }; 

// Polar stereographic
// TODO: not sure about the mappings for this one, see comment in 
//   NCDFWriteProjAttribs()

// Rotated Pole
// TODO: No GDAL equivalent of rotated pole? Doesn't seem to have an EPSG
//  code or WKT ... so unless some advanced proj4 features can be used 
//  seems to rule out.

// Stereographic
// TODO: Issues of how GDAL handles STEREOGRAPHIC vs OBLIQUE_STEREOGRAPHIC
// (GDAL seems to create Oblique_stereographic when you request stereographic
//   using Proj4, not entirely sure they're different projections)
//  Haven't been able to create a GDAL regular Stereographic to test yet.
static const oNetcdfSRS_PP poStMappings[] = {
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN},
    {LON_PROJ_ORIGIN, SRS_PP_CENTRAL_MERIDIAN},
    {SCALE_FACTOR_ORIGIN, SRS_PP_SCALE_FACTOR},  
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
  };

// Transverse Mercator
static const oNetcdfSRS_PP poTMMappings[] = {
    {SCALE_FACTOR_MERIDIAN, SRS_PP_SCALE_FACTOR},  
    {LONG_CENTRAL_MERIDIAN, SRS_PP_CENTRAL_MERIDIAN},
    {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN},
    {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
    {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
    {NULL, NULL}
  };

/* ET here the meta-mapping with all mappings, mapped by SRS_PT */

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
