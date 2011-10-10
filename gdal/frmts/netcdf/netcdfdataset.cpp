/******************************************************************************
 * $Id: netcdfdataset.cpp 23197 2011-10-07 00:34:31Z pds $
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

#include "netcdfdataset.h"
#include "cpl_error.h"
CPL_CVSID("$Id: netcdfdataset.cpp 23197 2011-10-07 00:34:31Z pds $");

void NCDFAddHistory(int fpImage, const char *pszAddHist, const char *pszOldHist);

/* PDS: Internal Utility function declarations */
void NCDFWriteProjAttribs(const OGR_SRSNode *poPROJCS,
                            const char* pszProjection,
                            const int fpImage, const int NCDFVarID);

void NCDFWriteProjAttribsFromMappings(const OGR_SRSNode *poPROJCS,
                                        const oNetcdfSRS_PP mappings[],
                                        const int fpImage, const int NCDFVarID);

const char* GetProjParamVal(const OGR_SRSNode *poPROJCS, 
                            const char* findParamStr);

/************************************************************************/
/* ==================================================================== */
/*                         netCDFRasterBand                             */
/* ==================================================================== */
/************************************************************************/

class netCDFRasterBand : public GDALPamRasterBand
{
    nc_type nc_datatype;
    int         nZId;
    int         nZDim;
    int		nLevel;
    int         nBandXPos;
    int         nBandYPos;
    int         *panBandZPos;
    int         *panBandZLev;
    int         bNoDataSet;
    double      dfNoDataValue;
    double      dfScale;
    double      dfOffset;
    CPLErr	    CreateBandMetadata( ); 
    
  public:

    netCDFRasterBand( netCDFDataset *poDS, 
		      int nZId, 
		      int nZDim,
		      int nLevel, 
		      int *panBandZLen,
		      int *panBandPos, 
		      int nBand );
    ~netCDFRasterBand( );
    virtual double          GetNoDataValue( int * );
    virtual CPLErr          SetNoDataValue( double );
    virtual double          GetOffset( int * );
    virtual CPLErr          SetOffset( double );
    virtual double          GetScale( int * );
    virtual CPLErr          SetScale( double );
    virtual CPLErr IReadBlock( int, int, void * );

};

/************************************************************************/ 
/*                             GetOffset()                              */ 
/************************************************************************/ 
double netCDFRasterBand::GetOffset( int *pbSuccess ) 
{ 
    if( pbSuccess != NULL ) 
        *pbSuccess = TRUE; 
	 
    return dfOffset; 
}

/************************************************************************/ 
/*                             SetOffset()                              */ 
/************************************************************************/ 
CPLErr netCDFRasterBand::SetOffset( double dfNewOffset ) 
{ 
    dfOffset = dfNewOffset; 
    return CE_None; 
}

/************************************************************************/ 
/*                              GetScale()                              */ 
/************************************************************************/ 
double netCDFRasterBand::GetScale( int *pbSuccess ) 
{ 
    if( pbSuccess != NULL ) 
        *pbSuccess = TRUE; 
    return dfScale; 
}

/************************************************************************/ 
/*                              SetScale()                              */ 
/************************************************************************/ 
CPLErr netCDFRasterBand::SetScale( double dfNewScale )  
{ 
    dfScale = dfNewScale; 
    return CE_None; 
} 

/************************************************************************/
/*                            GetMetadata()                             */
/************************************************************************/
char **netCDFDataset::GetMetadata( const char *pszDomain )
{
    if( pszDomain != NULL && EQUALN( pszDomain, "SUBDATASETS", 11 ) )
        return papszSubDatasets;
    else
        return GDALDataset::GetMetadata( pszDomain );
}

/************************************************************************/
/*                          GetProjectionRef()                          */
/************************************************************************/

const char * netCDFDataset::GetProjectionRef()
{
    if( bGotGeoTransform )
        return pszProjection;
    else
        return GDALPamDataset::GetProjectionRef();
}

/************************************************************************/
/*                           GetNoDataValue()                           */
/************************************************************************/

double netCDFRasterBand::GetNoDataValue( int * pbSuccess )

{
    if( pbSuccess )
        *pbSuccess = bNoDataSet;

    if( bNoDataSet )
        return dfNoDataValue;
    else
        return GDALPamRasterBand::GetNoDataValue( pbSuccess );
}

/************************************************************************/
/*                           SetNoDataValue()                           */
/************************************************************************/

CPLErr netCDFRasterBand::SetNoDataValue( double dfNoData )

{
    bNoDataSet = TRUE;
    dfNoDataValue = dfNoData;

    return CE_None;
}

/************************************************************************/
/*                         ~netCDFRasterBand()                          */
/************************************************************************/

netCDFRasterBand::~netCDFRasterBand()
{
    if( panBandZPos ) 
        CPLFree( panBandZPos );
    if( panBandZLev )
        CPLFree( panBandZLev );
}

/************************************************************************/
/*                         CreateBandMetadata()                         */
/************************************************************************/

CPLErr netCDFRasterBand::CreateBandMetadata( ) 
{
    char     szVarName[NC_MAX_NAME];
    char     szMetaName[NC_MAX_NAME];
    char     szMetaTemp[MAX_STR_LEN];
    int      nd;
    int      i,j;
    int      Sum  = 1;
    int      Taken = 0;
    int      result = 0;
    int      status;
    int      nVarID = -1;
    int      nDims;
    size_t   start[1];
    size_t   count[1];
    char     szTemp[NC_MAX_NAME];
    const char *pszValue;

    nc_type nVarType;
    netCDFDataset *poDS;

    poDS = (netCDFDataset *) this->poDS;
/* -------------------------------------------------------------------- */
/*      Compute all dimensions from Band number and save in Metadata    */
/* -------------------------------------------------------------------- */
    nc_inq_varname( poDS->cdfid, nZId, szVarName );
    nc_inq_varndims( poDS->cdfid, nZId, &nd );
/* -------------------------------------------------------------------- */
/*      Compute multidimention band position                            */
/*                                                                      */
/* BandPosition = (Total - sum(PastBandLevels) - 1)/sum(remainingLevels)*/
/* if Data[2,3,4,x,y]                                                   */
/*                                                                      */
/*  BandPos0 = (nBand ) / (3*4)                                         */
/*  BandPos1 = (nBand - BandPos0*(3*4) ) / (4)                          */
/*  BandPos2 = (nBand - BandPos0*(3*4) ) % (4)                          */
/* -------------------------------------------------------------------- */

    sprintf( szMetaName,"NETCDF_VARNAME");
    sprintf( szMetaTemp,"%s",szVarName);
    SetMetadataItem( szMetaName, szMetaTemp );
    if( nd == 3 ) {
        Sum *= panBandZLev[0];
    }

    for( i=0; i < nd-2 ; i++ ) {
        if( i != nd - 2 -1 ) {
            Sum = 1;
            for( j=i+1; j < nd-2; j++ ) {
                Sum *= panBandZLev[j];
            }
            result = (int) ( ( nLevel-Taken) / Sum );
        }
        else {
            result = (int) ( ( nLevel-Taken) % Sum );
        }
        
        strcpy(szVarName, poDS->papszDimName[poDS->paDimIds[
                       panBandZPos[i]]] );

        sprintf( szMetaName,"NETCDF_DIMENSION_%s",  szVarName );

        status=nc_inq_varid(poDS->cdfid,  
                            szVarName,
                            &nVarID );

/* -------------------------------------------------------------------- */
/*      Try to uppercase the first letter of the variable               */
/* -------------------------------------------------------------------- */

        if( status != NC_NOERR ) {
            szVarName[0]=(char) toupper(szVarName[0]);
            status=nc_inq_varid(poDS->cdfid,  
                                szVarName,
                                &nVarID );
        }

        status = nc_inq_vartype( poDS->cdfid, nVarID, &nVarType );

        nDims = 0;
        status = nc_inq_varndims( poDS->cdfid, nVarID, &nDims );

        if( nDims == 1 ) {
            count[0]=1;
            start[0]=result;
            switch( nVarType ) {
                case NC_SHORT:
                    short sData;
                    status =  nc_get_vara_short( poDS->cdfid, nVarID, 
                                                 start,
                                                 count, &sData );
                    sprintf( szMetaTemp,"%d", sData );
                    break;
                case NC_INT:
                    int nData;
                    status =  nc_get_vara_int( poDS->cdfid, nVarID, 
                                               start,
                                               count, &nData );
                    sprintf( szMetaTemp,"%d", nData );
                    break;
                case NC_FLOAT:
                    float fData;
                    status =  nc_get_vara_float( poDS->cdfid, nVarID, 
                                                 start,
                                                 count, &fData );
                    sprintf( szMetaTemp,"%f", fData );
                    break;
                case NC_DOUBLE:
                    double dfData;
                    status =  nc_get_vara_double( poDS->cdfid, nVarID, 
                                                  start,
                                                  count, &dfData);
                    sprintf( szMetaTemp,"%.16g", dfData );
                    break;
                default:
                    break;
            }
        }
        else
            sprintf( szMetaTemp,"%d", result+1);
	
        SetMetadataItem( szMetaName, szMetaTemp );

/* -------------------------------------------------------------------- */
/*      Fetch dimension units                                           */
/* -------------------------------------------------------------------- */

        strcpy( szTemp, szVarName );
        strcat( szTemp, "#units" );
        pszValue = CSLFetchNameValue(poDS->papszMetadata, szTemp);
        if( pszValue != NULL ) {
            if( EQUAL( pszValue, "T") ) { 
                strcpy( szTemp, szVarName );
                strcat( szTemp, "#original_units" );
                pszValue = CSLFetchNameValue(poDS->papszMetadata, szTemp);
                strcpy( szTemp, "NETCDF_");
                strcat( szTemp, szVarName );
                strcat( szTemp, "_original_units" );
                SetMetadataItem( szTemp, pszValue );
            }
            else {
                strcpy( szTemp, "NETCDF_");
                strcat( szTemp, szVarName  );
                strcat( szTemp, "_units" );
                SetMetadataItem( szTemp, pszValue );
            }
        }
        Taken += result * Sum;
    }

/* -------------------------------------------------------------------- */
/*      Get all other metadata                                          */
/* -------------------------------------------------------------------- */

    int nAtt=0;
    nc_type  atttype=NC_NAT;
    size_t   attlen;
    float fval;
    double dval;
    int ival;

    nc_inq_varnatts( poDS->cdfid, nZId, &nAtt );
    for( i=0; i < nAtt ; i++ ) {
    	status = nc_inq_attname( poDS->cdfid, nZId, 
    				 i, szTemp);
    	status = nc_inq_att( poDS->cdfid, nZId, 
    			     szTemp, &atttype, &attlen);
    	if(strcmp(szTemp,_FillValue) ==0) continue;
    	sprintf( szMetaTemp,"%s",szTemp);
    	switch( atttype ) {
    	case NC_CHAR:
    	    char *fillc;
    	    fillc = (char *) CPLCalloc( attlen+1, sizeof(char) );
    	    status=nc_get_att_text( poDS->cdfid, nZId,
    				    szTemp, fillc );
    	    SetMetadataItem( szMetaTemp, fillc );
    	    CPLFree(fillc);
    	    break;
    	case NC_INT:
    	    status = nc_get_att_int( poDS->cdfid, nZId,
    				     szTemp, &ival );
    	    sprintf( szTemp,"%d",ival);
    	    SetMetadataItem( szMetaTemp, szTemp );
    	    break;
    	case NC_FLOAT:
    	    status = nc_get_att_float( poDS->cdfid, nZId,
    				       szTemp, &fval );
    	    sprintf( szTemp,"%f",fval);
    	    SetMetadataItem( szMetaTemp, szTemp );
    	    break;
    	case NC_DOUBLE:
    	    status = nc_get_att_double( poDS->cdfid, nZId,
    					szTemp, &dval );
    	    sprintf( szTemp,"%.16g",dval);
    	    SetMetadataItem( szMetaTemp, szTemp );
    	    break;
    	default:
    	    break;
    	}
    }

    return CE_None;
}

/************************************************************************/
/*                          netCDFRasterBand()                          */
/************************************************************************/

netCDFRasterBand::netCDFRasterBand( netCDFDataset *poDS, 
                                    int nZId, 
                                    int nZDim,
                                    int nLevel, 
                                    int *panBandZLev, 
                                    int *panBandDimPos, 
                                    int nBand)

{
    double   dfNoData;
    int      bNoDataSet = FALSE;
    nc_type  vartype=NC_NAT;
    nc_type  atttype=NC_NAT;
    size_t   attlen;
    int      status;
    char     szNoValueName[MAX_STR_LEN];


    this->panBandZPos = NULL;
    this->panBandZLev = NULL;
    this->poDS = poDS;
    this->nBand = nBand;
    this->nZId = nZId;
    this->nZDim = nZDim;
    this->nLevel = nLevel;
    this->nBandXPos = panBandDimPos[0];
    this->nBandYPos = panBandDimPos[1];

/* -------------------------------------------------------------------- */
/*      Take care of all other dimmensions                              */
/* ------------------------------------------------------------------ */
    if( nZDim > 2 ) {
        this->panBandZPos = 
            (int *) CPLCalloc( nZDim-1, sizeof( int ) );
        this->panBandZLev = 
            (int *) CPLCalloc( nZDim-1, sizeof( int ) );

        for ( int i=0; i < nZDim - 2; i++ ){
            this->panBandZPos[i] = panBandDimPos[i+2];
            this->panBandZLev[i] = panBandZLev[i];
        }
    }
    CreateBandMetadata();
    bNoDataSet    = FALSE;
    dfNoDataValue = -9999.0;

    nBlockXSize   = poDS->GetRasterXSize( );
    nBlockYSize   = 1;

/* -------------------------------------------------------------------- */
/*      Get the type of the "z" variable, our target raster array.      */
/* -------------------------------------------------------------------- */
    if( nc_inq_var( poDS->cdfid, nZId, NULL, &nc_datatype, NULL, NULL,
                    NULL ) != NC_NOERR ){
        CPLError( CE_Failure, CPLE_AppDefined, 
                  "Error in nc_var_inq() on 'z'." );
        return;
    }

    if( (nc_datatype == NC_BYTE) ) 
        eDataType = GDT_Byte;
    else if( (nc_datatype == NC_UBYTE) ) 
        eDataType = GDT_Byte;
    else if( (nc_datatype == NC_CHAR) ) 
        eDataType = GDT_Byte;
    else if( nc_datatype == NC_SHORT )
        eDataType = GDT_Int16;
    else if( nc_datatype == NC_INT )
        eDataType = GDT_Int32;
    else if( nc_datatype == NC_FLOAT )
        eDataType = GDT_Float32;
    else if( nc_datatype == NC_DOUBLE )
        eDataType = GDT_Float64;
    else
    {
        if( nBand == 1 )
            CPLError( CE_Warning, CPLE_AppDefined, 
                      "Unsupported netCDF datatype (%d), treat as Float32.", 
                      (int) nc_datatype );
        eDataType = GDT_Float32;
    }
/* -------------------------------------------------------------------- */
/*      Find out what is No Data for this variable                      */
/* -------------------------------------------------------------------- */

    status = nc_inq_att( poDS->cdfid, nZId, 
                         _FillValue, &atttype, &attlen);

/* -------------------------------------------------------------------- */
/*      Look for either Missing_Value or _FillValue attributes          */
/* -------------------------------------------------------------------- */

    if( status == NC_NOERR ) {
        strcpy(szNoValueName, _FillValue );
    }
    else {
        status = nc_inq_att( poDS->cdfid, nZId, 
                             "missing_value", &atttype, &attlen );
        if( status == NC_NOERR ) {

            strcpy( szNoValueName, "missing_value" );
        }
    }

    nc_inq_vartype( poDS->cdfid, nZId, &vartype );

    if( status == NC_NOERR ) {
        switch( atttype ) {
            case NC_CHAR:
                char *fillc;
                fillc = (char *) CPLCalloc( attlen+1, sizeof(char) );
                status=nc_get_att_text( poDS->cdfid, nZId,
                                        szNoValueName, fillc );
                dfNoData = atof( fillc );
                CPLFree(fillc);
                break;
            case NC_SHORT:
                short sNoData;
                status = nc_get_att_short( poDS->cdfid, nZId,
                                           szNoValueName, &sNoData );
                dfNoData = (double) sNoData;
                break;
            case NC_INT:
                int nNoData;
                status = nc_get_att_int( poDS->cdfid, nZId,
                                         szNoValueName, &nNoData );
                dfNoData = (double) nNoData;
                break;
            case NC_FLOAT:
                float fNoData;
                status = nc_get_att_float( poDS->cdfid, nZId,
                                           szNoValueName, &fNoData );
                dfNoData = (double) fNoData;
                break;
            case NC_DOUBLE:
                status = nc_get_att_double( poDS->cdfid, nZId,
                                            szNoValueName, &dfNoData );
                break;
            default:
                break;
        }
        status = nc_get_att_double( poDS->cdfid, nZId, 
                                    szNoValueName, &dfNoData );
	
    } else {
        switch( vartype ) {
            case NC_BYTE:
                /* don't do default fill-values for bytes, too risky */
                dfNoData = 0.0;
                /* should print a warning as users might not be expecting this */
                /* CPLError(CE_Warning, 1,"GDAL netCDF driver is setting default NoData value to 0.0 for NC_BYTE data\n"); */
               break;
            case NC_CHAR:
                dfNoData = NC_FILL_CHAR;
                break;
            case NC_SHORT:
                dfNoData = NC_FILL_SHORT;
                break;
            case NC_INT:
                dfNoData = NC_FILL_INT;
                break;
            case NC_FLOAT:
                dfNoData = NC_FILL_FLOAT;
                break;
            case NC_DOUBLE:
                dfNoData = NC_FILL_DOUBLE;
                break;
            default:
                dfNoData = 0.0;
                break;
        }
	    bNoDataSet = TRUE;
    }
    SetNoDataValue( dfNoData );

    /* -------------------------------------------------------------------- */
    /* Attempt to fetch the scale_factor and add_offset attributes for the  */
    /* variable and set them.  If these values are not available, set       */
    /* offset to 0 and scale to 1                                           */
    /* -------------------------------------------------------------------- */
    double dfOff = 0.0; 
    double dfScale = 1.0; 
    if ( nc_inq_attid ( poDS->cdfid, nZId, "add_offset", NULL) == NC_NOERR ) { 
        nc_get_att_double( poDS->cdfid, nZId, "add_offset", &dfOff );
    }
    if ( nc_inq_attid ( poDS->cdfid, nZId, 
			"scale_factor", NULL) == NC_NOERR ) { 
	nc_get_att_double( poDS->cdfid, nZId, "scale_factor", &dfScale ); 
    }
    SetOffset( dfOff ); 
    SetScale( dfScale ); 
}

/************************************************************************/
/*                             IReadBlock()                             */
/************************************************************************/

CPLErr netCDFRasterBand::IReadBlock( int nBlockXOff, int nBlockYOff,
                                     void * pImage )

{
    int    nErr=-1;
    int    cdfid = ( ( netCDFDataset * ) poDS )->cdfid;
    size_t start[ MAX_NC_DIMS ];
    size_t edge[ MAX_NC_DIMS ];
    char   pszName[ MAX_STR_LEN ];
    int    i,j;
    int    Sum=-1;
    int    Taken=-1;
    int    nd;

    *pszName='\0';
    memset( start, 0, sizeof( start ) );
    memset( edge,  0, sizeof( edge )  );
    nc_inq_varndims ( cdfid, nZId, &nd );

/* -------------------------------------------------------------------- */
/*      Locate X, Y and Z position in the array                         */
/* -------------------------------------------------------------------- */
	
    start[nBandXPos] = 0;          // x dim can move arround in array
    // check y order
    if( ( ( netCDFDataset *) poDS )->bBottomUp ) {
        start[nBandYPos] = ( ( netCDFDataset * ) poDS )->ydim - 1 - nBlockYOff;
    } else {
        start[nBandYPos] = nBlockYOff; // y
    }
        
    edge[nBandXPos] = nBlockXSize; 
    edge[nBandYPos] = 1;

    if( nd == 3 ) {
        start[panBandZPos[0]]  = nLevel;     // z
        edge [panBandZPos[0]]  = 1;
    }

/* -------------------------------------------------------------------- */
/*      Compute multidimention band position                            */
/*                                                                      */
/* BandPosition = (Total - sum(PastBandLevels) - 1)/sum(remainingLevels)*/
/* if Data[2,3,4,x,y]                                                   */
/*                                                                      */
/*  BandPos0 = (nBand ) / (3*4)                                         */
/*  BandPos1 = (nBand - (3*4) ) / (4)                                   */
/*  BandPos2 = (nBand - (3*4) ) % (4)                                   */
/* -------------------------------------------------------------------- */
    if (nd > 3) 
    {
        Taken = 0;
        for( i=0; i < nd-2 ; i++ ) 
        {
            if( i != nd - 2 -1 ) {
                Sum = 1;
                for( j=i+1; j < nd-2; j++ ) {
                    Sum *= panBandZLev[j];
                }
                start[panBandZPos[i]] = (int) ( ( nLevel-Taken) / Sum );
                edge[panBandZPos[i]] = 1;
            } else {
                start[panBandZPos[i]] = (int) ( ( nLevel-Taken) % Sum );
                edge[panBandZPos[i]] = 1;
            }
            Taken += start[panBandZPos[i]] * Sum;
        }
    }

    // printf("TMP ET read block\n\n\n\n");
    // int shuffle_in, deflate_in, deflate_level;
    // nErr = nc_inq_var_deflate(cdfid, nZId, &shuffle_in,
    //                           &deflate_in, &deflate_level);
    // // if ( nErr == NC_NOERR ) {
    //     CPLDebug( "GDAL_netCDF", 
    //               "nc_inq_var_deflate: shuffle=%d deflate=%d deflate_level=%d\n",
    //               shuffle_in, deflate_in, deflate_level);
    // // }

    if( eDataType == GDT_Byte )
        nErr = nc_get_vara_uchar( cdfid, nZId, start, edge, 
                                  (unsigned char *) pImage );
    else if( eDataType == GDT_Int16 )
        nErr = nc_get_vara_short( cdfid, nZId, start, edge, 
                                  (short int *) pImage );
    else if( eDataType == GDT_Int32 )
    {
        if( sizeof(long) == 4 )
            nErr = nc_get_vara_long( cdfid, nZId, start, edge, 
                                     (long *) pImage );
        else
            nErr = nc_get_vara_int( cdfid, nZId, start, edge, 
                                    (int *) pImage );
    }
    else if( eDataType == GDT_Float32 ){
        nErr = nc_get_vara_float( cdfid, nZId, start, edge, 
                                  (float *) pImage );
        for( i=0; i<nBlockXSize; i++ ){
            if( CPLIsNan( ( (float *) pImage )[i] ) )
                ( (float *)pImage )[i] = (float) dfNoDataValue;
        }
    }
    else if( eDataType == GDT_Float64 ){
        nErr = nc_get_vara_double( cdfid, nZId, start, edge, 
                                   (double *) pImage );
        for( i=0; i<nBlockXSize; i++ ){
            if( CPLIsNan( ( (double *) pImage)[i] ) ) 
                ( (double *)pImage )[i] = dfNoDataValue;
        }

    }

    if( nErr != NC_NOERR )
    {
        CPLError( CE_Failure, CPLE_AppDefined, 
                  "netCDF scanline fetch failed: %s", 
                  nc_strerror( nErr ) );
        return CE_Failure;
    }
    else
        return CE_None;
}

/************************************************************************/
/* ==================================================================== */
/*				netCDFDataset				*/
/* ==================================================================== */
/************************************************************************/

/************************************************************************/
/*                           netCDFDataset()                            */
/************************************************************************/

netCDFDataset::netCDFDataset()

{
    papszMetadata    = NULL;	
    papszSubDatasets = NULL;
    bGotGeoTransform = FALSE;
    pszProjection    = NULL;
    cdfid            = 0;
    bBottomUp        = FALSE;
    nFileType        = NCDF_FILETYPE_NONE;
}


/************************************************************************/
/*                           ~netCDFDataset()                           */
/************************************************************************/

netCDFDataset::~netCDFDataset()

{
    FlushCache();

    CSLDestroy( papszMetadata );
    CSLDestroy( papszSubDatasets );

    CPLFree( pszProjection );

    if( cdfid ) 
        nc_close( cdfid );
}

/************************************************************************/
/*                           FetchCopyParm()                            */
/************************************************************************/

double netCDFDataset::FetchCopyParm( const char *pszGridMappingValue, 
                                     const char *pszParm, double dfDefault )

{
    char         szTemp[ MAX_NC_NAME ];
    const char  *pszValue;

    strcpy(szTemp,pszGridMappingValue);
    strcat( szTemp, "#" );
    strcat( szTemp, pszParm );
    pszValue = CSLFetchNameValue(papszMetadata, szTemp);

    if( pszValue )
    {
        return CPLAtofM(pszValue);
    }
    else
        return dfDefault;
}

/************************************************************************/
/*                           FetchStandardParallels()                   */
/************************************************************************/

char** netCDFDataset::FetchStandardParallels( const char *pszGridMappingValue )
{
    char         szTemp[ MAX_NC_NAME ];
    const char   *pszValue;
    char         **papszValues = NULL;
    //cf-1.0 tags
    strcpy( szTemp,pszGridMappingValue );
    strcat( szTemp, "#" );
    strcat( szTemp, STD_PARALLEL );
    pszValue = CSLFetchNameValue( papszMetadata, szTemp );
    if( pszValue != NULL )
        papszValues = CSLTokenizeString2( pszValue, ",", CSLT_STRIPLEADSPACES |
                                          CSLT_STRIPENDSPACES );
    //try gdal tags
    else
    {
        strcpy( szTemp, pszGridMappingValue );
        strcat( szTemp, "#" );
        strcat( szTemp, STD_PARALLEL_1 );

        pszValue = CSLFetchNameValue( papszMetadata, szTemp );
	
        if ( pszValue != NULL )
            papszValues = CSLAddString( papszValues, pszValue );
				    
        strcpy( szTemp,pszGridMappingValue );
        strcat( szTemp, "#" );
        strcat( szTemp, STD_PARALLEL_2 );

        pszValue = CSLFetchNameValue( papszMetadata, szTemp );
	
        if( pszValue != NULL )	
            papszValues = CSLAddString( papszValues, pszValue );
    }
    
    return papszValues;
}

/************************************************************************/
/*                           SetProjection()                            */
/************************************************************************/
void netCDFDataset::SetProjection( int var )
{
/* -------------------------------------------------------------------- */
/*      Set Projection                                                  */
/* -------------------------------------------------------------------- */

    size_t       start[2], edge[2];
    int          status;
    unsigned int i;
    const char   *pszValue;
    int          nVarProjectionID;
    char         szVarName[ MAX_NC_NAME ];
    char         szTemp[ MAX_NC_NAME ];
    char         szGridMappingName[ MAX_NC_NAME ];
    char         szGridMappingValue[ MAX_NC_NAME ];

    double       dfStdP1=0.0;
    double       dfStdP2=0.0;
    double       dfCenterLat;
    double       dfCenterLon;
    double       dfScale;
    double       dfFalseEasting;
    double       dfFalseNorthing;
    double       dfCentralMeridian;
    double       dfEarthRadius;
    double       dfInverseFlattening;
    double       dfLonPrimeMeridian;
    double       dfSemiMajorAxis;
    double       dfSemiMinorAxis;
    
    int          bGotGeogCS = FALSE;
    int          bGotCfSRS = FALSE;
    int          bGotGdalSRS = FALSE;

    OGRSpatialReference oSRS;
    int          nVarDimXID = -1;
    int          nVarDimYID = -1;
    double       *pdfXCoord;
    double       *pdfYCoord;
    char         szDimNameX[ MAX_NC_NAME ];
    char         szDimNameY[ MAX_NC_NAME ];
    int          nSpacingBegin;
    int          nSpacingMiddle;
    int          nSpacingLast;

    const char *pszWKT;
    const char *pszGeoTransform;
    char **papszGeoTransform=NULL;
    //char         *pszProjectionGDAL = NULL;

    netCDFDataset * poDS;
    poDS = this;

/* -------------------------------------------------------------------- */
/*      Get x/y range information.                                      */
/* -------------------------------------------------------------------- */

    poDS->adfGeoTransform[0] = 0.0;
    poDS->adfGeoTransform[1] = 1.0;
    poDS->adfGeoTransform[2] = 0.0;
    poDS->adfGeoTransform[3] = 0.0;
    poDS->adfGeoTransform[4] = 0.0;
    poDS->adfGeoTransform[5] = 1.0;
    poDS->pszProjection = NULL;
    

/* -------------------------------------------------------------------- */
/*      Look for grid_mapping metadata                                  */
/* -------------------------------------------------------------------- */

    strcpy( szGridMappingValue, "" );
    strcpy( szGridMappingName, "" );

    nc_inq_varname(  cdfid, var, szVarName );
    strcpy(szTemp,szVarName);
    strcat(szTemp,"#");
    strcat(szTemp,GRD_MAPPING);
    pszValue = CSLFetchNameValue(poDS->papszMetadata, szTemp);
    if( pszValue ) {
        strcpy(szGridMappingName,szTemp);
        strcpy(szGridMappingValue,pszValue);
    }

/* -------------------------------------------------------------------- */
/*      Look for dimension: lon                                         */
/* -------------------------------------------------------------------- */

    memset( szDimNameX, '\0', sizeof( char ) * MAX_NC_NAME );
    memset( szDimNameY, '\0', sizeof( char ) * MAX_NC_NAME );

    for( i = 0; (i < strlen( poDS->papszDimName[ poDS->nDimXid ] )  && 
                 i < 3 ); i++ ) {
        szDimNameX[i]=(char)tolower( ( poDS->papszDimName[poDS->nDimXid] )[i] );
    }
    szDimNameX[3] = '\0';
    for( i = 0; (i < strlen( poDS->papszDimName[ poDS->nDimYid ] )  && 
                 i < 3 ); i++ ) {
        szDimNameY[i]=(char)tolower( ( poDS->papszDimName[poDS->nDimYid] )[i] );
    }
    szDimNameY[3] = '\0';

/* -------------------------------------------------------------------- */
/*      Read grid_mapping information and set projections               */
/* -------------------------------------------------------------------- */

    if( !( EQUAL(szGridMappingName,"" ) ) ) {
        nc_inq_varid( cdfid, szGridMappingValue, &nVarProjectionID );
        poDS->ReadAttributes( cdfid, nVarProjectionID );
    
        strcpy( szTemp, szGridMappingValue );
        strcat( szTemp, "#" );
        strcat( szTemp, GRD_MAPPING_NAME );
        pszValue = CSLFetchNameValue(poDS->papszMetadata, szTemp);

        if( pszValue != NULL ) {

/* -------------------------------------------------------------------- */
/*      Check for datum/spheroid information                            */
/* -------------------------------------------------------------------- */
            dfEarthRadius = 
                poDS->FetchCopyParm( szGridMappingValue, 
                                     EARTH_RADIUS, 
                                     -1.0 );

            dfLonPrimeMeridian = 
                poDS->FetchCopyParm( szGridMappingValue,
                                     LONG_PRIME_MERIDIAN, 
                                     0.0 );

            dfInverseFlattening = 
                poDS->FetchCopyParm( szGridMappingValue, 
                                     INVERSE_FLATTENING, 
                                     -1.0 );
	    
            dfSemiMajorAxis = 
                poDS->FetchCopyParm( szGridMappingValue, 
                                     SEMI_MAJOR_AXIS, 
                                     -1.0 );
	    
            dfSemiMinorAxis = 
                poDS->FetchCopyParm( szGridMappingValue, 
                                     SEMI_MINOR_AXIS, 
                                     -1.0 );

            //see if semi-major exists if radius doesn't
            if( dfEarthRadius < 0.0 )
                dfEarthRadius = dfSemiMajorAxis;
	    
            //if still no radius, check old tag
            if( dfEarthRadius < 0.0 )
                dfEarthRadius = poDS->FetchCopyParm( szGridMappingValue, 
                                                     "spherical_earth_radius_meters",
                                                     -1.0 );

            //has radius value
            if( dfEarthRadius > 0.0 ) {
                //check for inv_flat tag
                if( dfInverseFlattening < 0.0 ) {
                    //no inv_flat tag, check for semi_minor
                    if( dfSemiMinorAxis < 0.0 ) {
                        //no way to get inv_flat, use sphere
                        oSRS.SetGeogCS( "unknown", 
                                        NULL, 
                                        "Sphere", 
                                        dfEarthRadius, 0.0 );
                        bGotGeogCS = TRUE;
                    }
                    else {
                        if( dfSemiMajorAxis < 0.0 )
                            dfSemiMajorAxis = dfEarthRadius;
                        //set inv_flat using semi_minor/major
                        dfInverseFlattening = 
                            1.0 / ( dfSemiMajorAxis - dfSemiMinorAxis ) / dfSemiMajorAxis;
                        oSRS.SetGeogCS( "unknown", 
                                        NULL, 
                                        "Spheroid", 
                                        dfEarthRadius, dfInverseFlattening );
                        bGotGeogCS = TRUE;
                    }
                }
                else {
                    oSRS.SetGeogCS( "unknown", 
                                    NULL, 
                                    "Spheroid", 
                                    dfEarthRadius, dfInverseFlattening );
                    bGotGeogCS = TRUE;
                }
            
             }
            //no radius, set as wgs84 as default?
            else {
                // This would be too indiscrimant.  But we should set
                // it if we know the data is geographic.
                //oSRS.SetWellKnownGeogCS( "WGS84" );
            }
	    		
/* -------------------------------------------------------------------- */
/*      Transverse Mercator                                             */
/* -------------------------------------------------------------------- */

            if( EQUAL( pszValue, TM ) ) {

                // dfScale = 
                //     poDS->FetchCopyParm( szGridMappingValue, 
                //                          SCALE_FACTOR_ORIGIN, 1.0 );
                dfScale = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         SCALE_FACTOR_MERIDIAN, 1.0 );

                dfCenterLon = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LONG_CENTRAL_MERIDIAN, 0.0 );

                dfCenterLat = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LAT_PROJ_ORIGIN, 0.0 );

                dfFalseEasting = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_EASTING, 0.0 );

                dfFalseNorthing = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_NORTHING, 0.0 );

                CPLDebug( "GDAL_netCDF", 
                          "SetTM( %g, %g, %g, %g, %g)\n",
                          dfCenterLat, 
                          dfCenterLon,
                          dfScale,
                          dfFalseEasting,
                          dfFalseNorthing  );

                bGotCfSRS = TRUE;
                oSRS.SetTM( dfCenterLat, 
                            dfCenterLon,
                            dfScale,
                            dfFalseEasting,
                            dfFalseNorthing );

                if( !bGotGeogCS )
                    oSRS.SetWellKnownGeogCS( "WGS84" );
            }

/* -------------------------------------------------------------------- */
/*      Albers Equal Area                                               */
/* -------------------------------------------------------------------- */

            if( EQUAL( pszValue, AEA ) ) {
                char **papszStdParallels = NULL;
		
                dfCenterLon = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LONG_CENTRAL_MERIDIAN, 0.0 );

                dfCenterLat = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LAT_PROJ_ORIGIN, 0.0 );

                // dfScale = 
                //     poDS->FetchCopyParm( szGridMappingValue, 
                //                          SCALE_FACTOR, 1.0 );

                dfFalseEasting = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_EASTING, 0.0 );

                dfFalseNorthing = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_NORTHING, 0.0 );
		
                papszStdParallels = 
                    FetchStandardParallels( szGridMappingValue );

                if( papszStdParallels != NULL ) {
		  
                    if ( CSLCount( papszStdParallels ) == 1 ) {
                        dfStdP1 = CPLAtofM( papszStdParallels[0] );
                        dfStdP2 = dfStdP1;
                    }
		
                    else if( CSLCount( papszStdParallels ) == 2 ) {
                        dfStdP1 = CPLAtofM( papszStdParallels[0] );
                        dfStdP2 = CPLAtofM( papszStdParallels[1] );
                    }
                }
                //old default
                else {
                    dfStdP1 = 
                        poDS->FetchCopyParm( szGridMappingValue, 
                                             STD_PARALLEL_1, 0.0 );

                    dfStdP2 = 
                        poDS->FetchCopyParm( szGridMappingValue, 
                                             STD_PARALLEL_2, 0.0 );
                }

                dfCenterLat = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LAT_PROJ_ORIGIN, 0.0 );

                // PDS
                // Ideas for doing this more cleverly if we had defined
                //  in a big array
                //    int* res; = (Allocate to correct length)
                //for each mapped value
                /*
                    //TODO: upgrade this func to handle STD_PARALLEL
                    //Special cases
                    val = poDS->FetchCopyParm( szGridMappingValue, 
                                         mappings[iMap].NCDF_ATT,
                                         mappings[iMap].defVal );
                    res[iMap] = val
                    oSRS.SetACEA( res[0], res[1], res[2], res[3], res[4],
                        res[5], res[6] );
                //   clean up res array
                */                         

                bGotCfSRS = TRUE;
                oSRS.SetACEA( dfStdP1, dfStdP2, dfCenterLat, dfCenterLon,
                              dfFalseEasting, dfFalseNorthing );


                if( !bGotGeogCS )
                    oSRS.SetWellKnownGeogCS( "WGS84" );

                CSLDestroy( papszStdParallels );
            }

/* -------------------------------------------------------------------- */
/*      Cylindrical Equal Area                                          */
/* -------------------------------------------------------------------- */

            else if( EQUAL( pszValue, CEA ) || EQUAL( pszValue, LCEA ) ) {
                dfStdP1 = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         STD_PARALLEL_1, 0.0 );
                dfCentralMeridian = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LONG_CENTRAL_MERIDIAN, 0.0 );

                dfFalseEasting = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_EASTING, 0.0 );

                dfFalseNorthing = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_NORTHING, 0.0 );
		
                bGotCfSRS = TRUE;
                oSRS.SetCEA( dfStdP1, dfCentralMeridian,
                             dfFalseEasting, dfFalseNorthing );

                if( !bGotGeogCS )
                    oSRS.SetWellKnownGeogCS( "WGS84" );
		
            }

/* -------------------------------------------------------------------- */
/*      lambert_azimuthal_equal_area                                    */
/* -------------------------------------------------------------------- */
            else if( EQUAL( pszValue, LAEA ) ) {
                dfCenterLon = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LON_PROJ_ORIGIN, 0.0 );

                dfCenterLat = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LAT_PROJ_ORIGIN, 0.0 );

                dfFalseEasting = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_EASTING, 0.0 );

                dfFalseNorthing = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_NORTHING, 0.0 );

                /*
                  dfLonOrig =
                  poDS->FetchCopyParm( szGridMappingValue, 
                  LON_PROJ_ORIGIN, 0.0 );

                  dfLatOrig =
                  poDS->FetchCopyParm( szGridMappingValue, 
                  LAT_PROJ_ORIGIN, 0.0 );

                  dfScaleFactorOrig = 
                  poDS->FetchCopyParm( szGridMappingValue,
                  SCALE_FACTOR_ORIGN, 0.0 );

                  dfProjXOrig =
                  poDS->FetchCopyParm( szGridMappingValue,
                  PROJ_X_ORIGIN, 0.0 );

                  dfProjYOrig =
                  poDS->FetchCopyParm( szGridMappingValue,
                  PROJ_Y_ORIGIN, 0.0 );

                  dfFalseEasting = 
                  poDS->FetchCopyParm( szGridMappingValue, 
                  FALSE_EASTING, 0.0 );

                  dfFalseNorthing = 
                  poDS->FetchCopyParm( szGridMappingValue, 
                  FALSE_NORTHING, 0.0 );
                */
                oSRS.SetProjCS( "LAEA (WGS84) " );
		
                bGotCfSRS = TRUE;
                oSRS.SetLAEA( dfCenterLat, dfCenterLon,
                              dfFalseEasting, dfFalseNorthing );

                if( !bGotGeogCS )
                    oSRS.SetWellKnownGeogCS( "WGS84" );
		
            }

/* -------------------------------------------------------------------- */
/*      Azimuthal Equidistant                                           */
/* -------------------------------------------------------------------- */
            else if( EQUAL( pszValue, AE ) ) {
                dfCenterLon = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LON_PROJ_ORIGIN, 0.0 );

                dfCenterLat = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LAT_PROJ_ORIGIN, 0.0 );

                dfFalseEasting = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_EASTING, 0.0 );

                dfFalseNorthing = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_NORTHING, 0.0 );

                bGotCfSRS = TRUE;
                oSRS.SetAE( dfCenterLat, dfCenterLon,
                            dfFalseEasting, dfFalseNorthing );

                if( !bGotGeogCS )
                    oSRS.SetWellKnownGeogCS( "WGS84" );
		
            }

/* -------------------------------------------------------------------- */
/*      Lambert conformal conic                                         */
/* -------------------------------------------------------------------- */
            else if( EQUAL( pszValue, L_C_CONIC ) ) {
		
                char **papszStdParallels = NULL;
		
                dfCenterLon = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LONG_CENTRAL_MERIDIAN, 0.0 );

                dfCenterLat = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LAT_PROJ_ORIGIN, 0.0 );

                /* this is not CF!!! */
                dfScale = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         SCALE_FACTOR, 1.0 );

                dfFalseEasting = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_EASTING, 0.0 );

                dfFalseNorthing = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_NORTHING, 0.0 );
		
                papszStdParallels = 
                    FetchStandardParallels( szGridMappingValue );

                if( papszStdParallels != NULL ) {
		  
                   if ( CSLCount( papszStdParallels ) == 1 ) {
                        dfStdP1 = CPLAtofM( papszStdParallels[0] );
                        dfStdP2 = dfStdP1;
                        /* should use dfStdP1 and dfStdP2 instead of dfScale */
                        oSRS.SetLCC1SP( dfCenterLat, dfCenterLon, dfScale, 
                                        dfFalseEasting, dfFalseNorthing );
                    }
		
                    else if( CSLCount( papszStdParallels ) == 2 ) {
                        dfStdP1 = CPLAtofM( papszStdParallels[0] );
                        dfStdP2 = CPLAtofM( papszStdParallels[1] );
                        oSRS.SetLCC( dfStdP1, dfStdP2, dfCenterLat, dfCenterLon,
                                     dfFalseEasting, dfFalseNorthing );
                    }
                }
                //old default
                else {
                    dfStdP1 = 
                        poDS->FetchCopyParm( szGridMappingValue, 
                                             STD_PARALLEL_1, 0.0 );

                    dfStdP2 = 
                        poDS->FetchCopyParm( szGridMappingValue, 
                                             STD_PARALLEL_2, 0.0 );

                    oSRS.SetLCC( dfStdP1, dfStdP2, dfCenterLat, dfCenterLon,
                                 dfFalseEasting, dfFalseNorthing );
                }				

                bGotCfSRS = TRUE;
                if( !bGotGeogCS )
                    oSRS.SetWellKnownGeogCS( "WGS84" );

                CSLDestroy( papszStdParallels );
            }
		
/* -------------------------------------------------------------------- */
/*      Is this Latitude/Longitude Grid explicitly                      */
/* -------------------------------------------------------------------- */
	    
            else if ( EQUAL ( pszValue, LATITUDE_LONGITUDE ) ) {
                bGotCfSRS = TRUE;
                if( !bGotGeogCS )
                    oSRS.SetWellKnownGeogCS( "WGS84" );
            }
/* -------------------------------------------------------------------- */
/*      Mercator                                                        */
/* -------------------------------------------------------------------- */
		  
            else if ( EQUAL ( pszValue, MERCATOR ) ) {
                dfCenterLon = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LON_PROJ_ORIGIN, 0.0 );
	      
                dfCenterLat = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LAT_PROJ_ORIGIN, 0.0 );

                dfScale = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         SCALE_FACTOR_ORIGIN,
                                         1.0 );

                dfFalseEasting = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_EASTING, 0.0 );

                dfFalseNorthing = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_NORTHING, 0.0 );

                bGotCfSRS = TRUE;

                oSRS.SetMercator( dfCenterLat, dfCenterLon, dfScale, 
                                  dfFalseEasting, dfFalseNorthing );

                if( !bGotGeogCS )
                    oSRS.SetWellKnownGeogCS( "WGS84" );
            }

/* -------------------------------------------------------------------- */
/*      Orthographic                                                    */
/* -------------------------------------------------------------------- */
		  
            else if ( EQUAL ( pszValue, ORTHOGRAPHIC ) ) {
                dfCenterLon = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LON_PROJ_ORIGIN, 0.0 );
	      
                dfCenterLat = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LAT_PROJ_ORIGIN, 0.0 );

                dfFalseEasting = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_EASTING, 0.0 );

                dfFalseNorthing = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_NORTHING, 0.0 );

                bGotCfSRS = TRUE;

                oSRS.SetOrthographic( dfCenterLat, dfCenterLon, 
                                      dfFalseEasting, dfFalseNorthing );

                if( !bGotGeogCS )
                    oSRS.SetWellKnownGeogCS( "WGS84" );
            }

/* -------------------------------------------------------------------- */
/*      Polar Stereographic                                             */
/* -------------------------------------------------------------------- */
		  
            else if ( EQUAL ( pszValue, POLAR_STEREO ) ) {

                dfCenterLon = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LON_PROJ_ORIGIN, 0.0 );

                dfCenterLat = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LAT_PROJ_ORIGIN, 0.0 );
		
                dfScale = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         SCALE_FACTOR_ORIGIN, 
                                         1.0 );

                dfFalseEasting = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_EASTING, 0.0 );

                dfFalseNorthing = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_NORTHING, 0.0 );

                bGotCfSRS = TRUE;
                oSRS.SetPS( dfCenterLat, dfCenterLon, dfScale, 
                            dfFalseEasting, dfFalseNorthing );

                if( !bGotGeogCS )
                    oSRS.SetWellKnownGeogCS( "WGS84" );
            }

/* -------------------------------------------------------------------- */
/*      Stereographic                                                   */
/* -------------------------------------------------------------------- */
		  
            else if ( EQUAL ( pszValue, STEREO ) ) {
	        
                dfCenterLon = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LON_PROJ_ORIGIN, 0.0 );
	      
                dfCenterLat = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         LAT_PROJ_ORIGIN, 0.0 );

                dfScale = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         SCALE_FACTOR_ORIGIN,
                                         1.0 );

                dfFalseEasting = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_EASTING, 0.0 );

                dfFalseNorthing = 
                    poDS->FetchCopyParm( szGridMappingValue, 
                                         FALSE_NORTHING, 0.0 );

                bGotCfSRS = TRUE;
                oSRS.SetStereographic( dfCenterLat, dfCenterLon, dfScale, 
                                       dfFalseEasting, dfFalseNorthing );

                if( !bGotGeogCS )
                    oSRS.SetWellKnownGeogCS( "WGS84" );
            }
  
/* -------------------------------------------------------------------- */
/*      Is this Latitude/Longitude Grid, default                        */
/* -------------------------------------------------------------------- */
	    
        } else if( EQUAL( szDimNameX,"lon" ) ) {
            oSRS.SetWellKnownGeogCS( "WGS84" );

        } else {
            // This would be too indiscrimant.  But we should set
            // it if we know the data is geographic.
            //oSRS.SetWellKnownGeogCS( "WGS84" );
        }
    }
/* -------------------------------------------------------------------- */
/*      Read projection coordinates                                     */
/* -------------------------------------------------------------------- */

    nc_inq_varid( cdfid, poDS->papszDimName[nDimXid], &nVarDimXID );
    nc_inq_varid( cdfid, poDS->papszDimName[nDimYid], &nVarDimYID );
    
    if( ( nVarDimXID != -1 ) && ( nVarDimYID != -1 ) ) {
        pdfXCoord = (double *) CPLCalloc( xdim, sizeof(double) );
        pdfYCoord = (double *) CPLCalloc( ydim, sizeof(double) );
    
/* -------------------------------------------------------------------- */
/*      Is pixel spacing is uniform accross the map?                    */
/* -------------------------------------------------------------------- */
        start[0] = 0;
        edge[0]  = xdim;
	
        status = nc_get_vara_double( cdfid, nVarDimXID, 
                                     start, edge, pdfXCoord);
        edge[0]  = ydim;
        status = nc_get_vara_double( cdfid, nVarDimYID, 
                                     start, edge, pdfYCoord);

/* -------------------------------------------------------------------- */
/*      Check Longitude                                                 */
/* -------------------------------------------------------------------- */

        nSpacingBegin   = (int) poDS->rint((pdfXCoord[1]-pdfXCoord[0]) * 1000);
	
        nSpacingMiddle  = (int) poDS->rint((pdfXCoord[xdim / 2] - 
                                            pdfXCoord[(xdim / 2) + 1]) * 1000);
	
        nSpacingLast    = (int) poDS->rint((pdfXCoord[xdim - 2] - 
                                            pdfXCoord[xdim-1]) * 1000);
	
        if( ( abs( nSpacingBegin )  ==  abs( nSpacingLast )     )  &&
            ( abs( nSpacingBegin )  ==  abs( nSpacingMiddle )   ) &&
            ( abs( nSpacingMiddle ) ==  abs( nSpacingLast )     ) ) {

/* -------------------------------------------------------------------- */
/*      Longitude is equaly spaced, check lattitde                      */
/* -------------------------------------------------------------------- */
            nSpacingBegin   = (int) poDS->rint((pdfYCoord[1]-pdfYCoord[0]) * 
                                               1000); 
	    
            nSpacingMiddle  = (int) poDS->rint((pdfYCoord[ydim / 2] - 
                                                pdfYCoord[(ydim / 2) + 1]) * 
                                               1000);
	    
            nSpacingLast    = (int) poDS->rint((pdfYCoord[ydim - 2] - 
                                                pdfYCoord[ydim-1]) * 
                                               1000);

		    
/* -------------------------------------------------------------------- */
/*   For Latitude  we allow an error of 0.1 degrees for gaussion        */
/*   gridding                                                           */
/* -------------------------------------------------------------------- */

            if((( abs( abs(nSpacingBegin) - abs(nSpacingLast) ) )   < 100 ) &&
               (( abs( abs(nSpacingBegin) -  abs(nSpacingMiddle) ) ) < 100 ) &&
               (( abs( abs(nSpacingMiddle) - abs(nSpacingLast) ) )   < 100) ) {

                if( ( abs( nSpacingBegin )  !=  abs( nSpacingLast )     )  ||
                    ( abs( nSpacingBegin )  !=  abs( nSpacingMiddle )   ) ||
                    ( abs( nSpacingMiddle ) !=  abs( nSpacingLast )     ) ) {
		    
                    CPLError(CE_Warning, 1,"Latitude grid not spaced evenly.\nSeting projection for grid spacing is within 0.1 degrees threshold.\n");

                }
/* -------------------------------------------------------------------- */
/*      We have gridded data s we can set the Gereferencing info.       */
/* -------------------------------------------------------------------- */

/* -------------------------------------------------------------------- */
/*      Enable GeoTransform                                             */
/* -------------------------------------------------------------------- */
                /* ----------------------------------------------------------*/
                /*    In the following "actual_range" and "node_offset"      */
                /*    are attributes used by netCDF files created by GMT.    */
                /*    If we find them we know how to proceed. Else, use      */
                /*    the original algorithm.                                */
                /* --------------------------------------------------------- */
                double	dummy[2], xMinMax[2], yMinMax[2];
                int	node_offset = 0;

                poDS->bGotGeoTransform = TRUE;

                nc_get_att_int (cdfid, NC_GLOBAL, "node_offset", &node_offset);

                if (!nc_get_att_double (cdfid, nVarDimXID, "actual_range", dummy)) {
                    xMinMax[0] = dummy[0];		
                    xMinMax[1] = dummy[1];
                }
                else {
                    xMinMax[0] = pdfXCoord[0];
                    xMinMax[1] = pdfXCoord[xdim-1];
                    node_offset = 0;
                }

                if (!nc_get_att_double (cdfid, nVarDimYID, "actual_range", dummy)) {
                    yMinMax[0] = dummy[0];		
                    yMinMax[1] = dummy[1];
                }
                else {
                    yMinMax[0] = pdfYCoord[0];	
                    yMinMax[1] = pdfYCoord[ydim-1];
                    node_offset = 0;
                }

                 /* for CF-1 conventions, assume bottom first */
                // if( ( EQUAL( szDimNameY, "lat" ) || EQUAL( szDimNameY, "y" ) )
                //     && pdfYCoord[0] < pdfYCoord[1] )
                //     poDS->bBottomUp = TRUE;

                /* Check for Bottom-up for all files (not just CF) */
                /* Also, don't check for dimname "lon" or "y", because CF doesn't specify valid dimnames */
                if ( pdfYCoord[0] < pdfYCoord[1] )
                    poDS->bBottomUp = TRUE;

                /* Check for reverse order of y-coordinate */
                if ( yMinMax[0] > yMinMax[1] ) {
                    dummy[0] = yMinMax[1];
                    dummy[1] = yMinMax[0];
                    yMinMax[0] = dummy[0];
                    yMinMax[1] = dummy[1];
                }
                CPLDebug( "GDAL_netCDF", "bBottomUp = %d\n", poDS->bBottomUp );

                /* ----------------------------------------------------------*/
                /*    Many netcdf files are weather files distributed        */
                /*    in km for the x/y resolution.  This isn't perfect,     */
                /*    but geotransforms can be terribly off if this isn't    */
                /*    checked and accounted for.  Maybe one more level of    */
                /*    checking (grid_mapping_value#GRIB_param_Dx, or         */
                /*    x#grid_spacing), but those are not cf tags.            */
                /*    Have to change metadata value if change Create() to    */
                /*    write cf tags                                          */
                /* ----------------------------------------------------------*/
                
                //check units for x and y, expand to other values 
                //and conversions.
                if( oSRS.IsProjected( ) ) {
                    strcpy( szTemp, "x" );
                    strcat( szTemp, "#units" );
                    pszValue = CSLFetchNameValue( poDS->papszMetadata, 
                                                  szTemp );
                    if( pszValue != NULL ) {
                        if( EQUAL( pszValue, "km" ) ) {
                            xMinMax[0] = xMinMax[0] * 1000;
                            xMinMax[1] = xMinMax[1] * 1000;
                        }
                    }
                    strcpy( szTemp, "y" );
                    strcat( szTemp, "#units" );
                    pszValue = CSLFetchNameValue( poDS->papszMetadata, 
                                                  szTemp );
                    if( pszValue != NULL ) {
                        if( EQUAL( pszValue, "km" ) ) {
                            yMinMax[0] = yMinMax[0] * 1000;
                            yMinMax[1] = yMinMax[1] * 1000;
                        }
                    }
                }

                poDS->adfGeoTransform[0] = xMinMax[0];
                poDS->adfGeoTransform[2] = 0;
                poDS->adfGeoTransform[3] = yMinMax[1];
                poDS->adfGeoTransform[4] = 0;
                poDS->adfGeoTransform[1] = ( xMinMax[1] - xMinMax[0] ) / 
                    ( poDS->nRasterXSize + (node_offset - 1) );
                poDS->adfGeoTransform[5] = ( yMinMax[0] - yMinMax[1] ) / 
                    ( poDS->nRasterYSize + (node_offset - 1) );

/* -------------------------------------------------------------------- */
/*     Compute the center of the pixel                                  */
/* -------------------------------------------------------------------- */
                if ( !node_offset ) {	// Otherwise its already the pixel center
                    poDS->adfGeoTransform[0] -= (poDS->adfGeoTransform[1] / 2);
                    poDS->adfGeoTransform[3] -= (poDS->adfGeoTransform[5] / 2);
                }

                oSRS.exportToWkt( &(poDS->pszProjection) );
                // oSRS.exportToPrettyWkt( &(poDS->pszProjection) );
		    
            } 
        }

        CPLFree( pdfXCoord );
        CPLFree( pdfYCoord );
    }


/* -------------------------------------------------------------------- */
/*      Is this a netCDF file created by GDAL?                          */
/*      NOTE: this should perhaps be skipped if info already found      */
/* -------------------------------------------------------------------- */
    if( !EQUAL( szGridMappingValue, "" )  ) {
        strcpy( szTemp,szGridMappingValue);
        strcat( szTemp, "#" );
        strcat( szTemp, NCDF_SPATIAL_REF);
        pszWKT = CSLFetchNameValue(poDS->papszMetadata, szTemp);
	
        if( pszWKT != NULL ) {
	    
/* -------------------------------------------------------------------- */
/*      Compare CRS obtained from CF attributes and GDAL WKT            */
/*      If possible use the more complete GDAL WKT                      */
/* -------------------------------------------------------------------- */
           // pszProjectionGDAL = CPLStrdup( pszWKT );
            if ( ! bGotCfSRS ) {   
                bGotGdalSRS = TRUE;
                poDS->pszProjection = CPLStrdup( pszWKT );
            }
            else { /* use the SRS from GDAL if it doesn't conflict with the one from CF */
                char *pszProjectionGDAL = CPLStrdup( pszWKT );
                OGRSpatialReference oSRSGDAL;
                oSRSGDAL.importFromWkt( &pszProjectionGDAL );
                /* set datum to unknown, bug #4281 */
                if ( oSRSGDAL.GetAttrNode( "DATUM" ) )
                    oSRSGDAL.GetAttrNode( "DATUM" )->GetChild(0)->SetValue( "unknown" );
                if ( oSRS.IsSame(&oSRSGDAL) ) {
                    printf("ARE SAME, using GDAL WKT\n");
                    bGotGdalSRS = TRUE;
                    poDS->pszProjection = CPLStrdup( pszWKT );
                }
                // else {  /* we got a GDAL SRS but will not use it */
                //     printf("are not same...");
                //     printf("geog same: %d\n",oSRS.IsSameGeogCS(&oSRSGDAL));
                // }
            }

/* -------------------------------------------------------------------- */
/*      Look for GeoTransform Array, if not found previously            */
/* -------------------------------------------------------------------- */
            if ( !bGotGeoTransform ) {

            strcpy(szTemp,szGridMappingValue);
            strcat( szTemp, "#" );
            strcat( szTemp, "GeoTransform");
	    
            pszGeoTransform = CSLFetchNameValue(poDS->papszMetadata, szTemp);	    

            if( pszGeoTransform != NULL ) {
                CPLDebug( "GDAL_netCDF", "looking for geotransform, got %s\n", pszGeoTransform );
                papszGeoTransform = CSLTokenizeString2( pszGeoTransform,
                                                        " ", 
                                                        CSLT_HONOURSTRINGS );
                poDS->bGotGeoTransform   = TRUE;
		
                poDS->adfGeoTransform[0] = atof( papszGeoTransform[0] );
                poDS->adfGeoTransform[1] = atof( papszGeoTransform[1] );
                poDS->adfGeoTransform[2] = atof( papszGeoTransform[2] );
                poDS->adfGeoTransform[3] = atof( papszGeoTransform[3] );
                poDS->adfGeoTransform[4] = atof( papszGeoTransform[4] );
                poDS->adfGeoTransform[5] = atof( papszGeoTransform[5] );
/* -------------------------------------------------------------------- */
/*      Look for corner array values                                    */
/* -------------------------------------------------------------------- */
            } else {
                double dfNN=0.0, dfSN=0.0, dfEE=0.0, dfWE=0.0;
                CPLDebug( "GDAL_netCDF", "looking for geotransform values\n" );
                strcpy(szTemp,szGridMappingValue);
                strcat( szTemp, "#" );
                strcat( szTemp, "Northernmost_Northing");
                pszValue = CSLFetchNameValue(poDS->papszMetadata, szTemp);
		
                if( pszValue != NULL ) {
                    dfNN = atof( pszValue );
                }
                strcpy(szTemp,szGridMappingValue);
                strcat( szTemp, "#" );
                strcat( szTemp, "Southernmost_Northing");
                pszValue = CSLFetchNameValue(poDS->papszMetadata, szTemp);
		
                if( pszValue != NULL ) {
                    dfSN = atof( pszValue );
                }
		
                strcpy(szTemp,szGridMappingValue);
                strcat( szTemp, "#" );
                strcat( szTemp, "Easternmost_Easting");
                pszValue = CSLFetchNameValue(poDS->papszMetadata, szTemp);
		
                if( pszValue != NULL ) {
                    dfEE = atof( pszValue );
                }
		
                strcpy(szTemp,szGridMappingValue);
                strcat( szTemp, "#" );
                strcat( szTemp, "Westernmost_Easting");
                pszValue = CSLFetchNameValue(poDS->papszMetadata, szTemp);
		
                if( pszValue != NULL ) {
                    dfWE = atof( pszValue );
                }
		
                adfGeoTransform[0] = dfWE;
                adfGeoTransform[1] = (dfEE - dfWE) / 
                    ( poDS->GetRasterXSize() - 1 );
                adfGeoTransform[2] = 0.0;
                adfGeoTransform[3] = dfNN;
                adfGeoTransform[4] = 0.0;
                adfGeoTransform[5] = (dfSN - dfNN) / 
                    ( poDS->GetRasterYSize() - 1 );
/* -------------------------------------------------------------------- */
/*     Compute the center of the pixel                                  */
/* -------------------------------------------------------------------- */
                adfGeoTransform[0] = dfWE
                    - (adfGeoTransform[1] / 2);

                adfGeoTransform[3] = dfNN
                    - (adfGeoTransform[5] / 2);


                bGotGeoTransform = TRUE;
            }
            CSLDestroy( papszGeoTransform );
            }
        }
    }
    
/* -------------------------------------------------------------------- */
/*     Search for Well-known GeogCS if only got CF WKT                  */
/* -------------------------------------------------------------------- */
    if ( bGotCfSRS && ! bGotGdalSRS ) {
        /* ET - could use a more exhaustive method by scanning all EPSG codes in data/gcs.csv */
        /* as proposed by Even in the gdal-dev mailing list "help for comparing two WKT" */
        /* this code could be contributed to a new function */
        /* OGRSpatialReference * OGRSpatialReference::FindMatchingGeogCS( const OGRSpatialReference *poOther ) */ 
        const char *pszWKGCSList[] = { "WGS84", "WGS72", "NAD27", "NAD83" };
        char *pszWKGCS = NULL;
        oSRS.exportToPrettyWkt( &pszWKGCS );
        for( size_t i=0; i<sizeof(pszWKGCSList)/8; i++ ) {
            pszWKGCS = CPLStrdup( pszWKGCSList[i] );
            OGRSpatialReference oSRSTmp;
            oSRSTmp.SetWellKnownGeogCS( pszWKGCSList[i] );
            /* set datum to unknown, bug #4281 */
            if ( oSRSTmp.GetAttrNode( "DATUM" ) )
                oSRSTmp.GetAttrNode( "DATUM" )->GetChild(0)->SetValue( "unknown" );
            oSRSTmp.GetRoot()->StripNodes( "AXIS" );
            oSRSTmp.GetRoot()->StripNodes( "AUTHORITY" );
            oSRSTmp.GetRoot()->StripNodes( "EXTENSION" );
            oSRSTmp.exportToPrettyWkt( &pszWKGCS );
            if ( oSRS.IsSameGeogCS(&oSRSTmp) ) {
                oSRS.SetWellKnownGeogCS( pszWKGCSList[i] );
                oSRS.exportToWkt( &(poDS->pszProjection) );
                // oSRS.exportToPrettyWkt( &pszWKGCS );
                // printf("TMP ET set WKT to\n[%s]\n",pszWKGCS);
            }
        }
    }
}
/************************************************************************/
/*                          GetGeoTransform()                           */
/************************************************************************/

CPLErr netCDFDataset::GetGeoTransform( double * padfTransform )

{
    memcpy( padfTransform, adfGeoTransform, sizeof(double) * 6 );
    if( bGotGeoTransform )
        return CE_None;
    else
        return GDALPamDataset::GetGeoTransform( padfTransform );;
}

/************************************************************************/
/*                                rint()                                */
/************************************************************************/

double netCDFDataset::rint( double dfX)
{
    if( dfX > 0 ) {
        int nX = (int) (dfX+0.5);
        if( nX % 2 ) {
            double dfDiff = dfX - (double)nX;
            if( dfDiff == -0.5 )
                return double( nX-1 );
        }
        return double( nX );
    } else {
        int nX= (int) (dfX-0.5);
        if( nX % 2 ) {
            double dfDiff = dfX - (double)nX;
            if( dfDiff == 0.5 )
                return double(nX+1);
        }
        return double(nX);
    }
}

/************************************************************************/
/*                        ReadAttributes()                              */
/************************************************************************/
CPLErr netCDFDataset::SafeStrcat(char** ppszDest, char* pszSrc, size_t* nDestSize)
{
    /* Reallocate the data string until the content fits */
    while(*nDestSize < (strlen(*ppszDest) + strlen(pszSrc) + 1)) {
        (*nDestSize) *= 2;
        *ppszDest = (char*) CPLRealloc((void*) *ppszDest, *nDestSize);
    }
    strcat(*ppszDest, pszSrc);
    
    return CE_None;
}

CPLErr netCDFDataset::ReadAttributes( int cdfid, int var)

{
    char    szAttrName[ NC_MAX_NAME ];
    char    szVarName [ NC_MAX_NAME ];
    char    szMetaName[ NC_MAX_NAME * 2 ];
    char    *pszMetaTemp = NULL;
    size_t  nMetaTempSize;
    nc_type nAttrType;
    size_t  nAttrLen, m;
    int     nbAttr;
    char    szTemp[ MAX_STR_LEN ];

    nc_inq_varnatts( cdfid, var, &nbAttr );
    if( var == NC_GLOBAL ) {
        strcpy( szVarName,"NC_GLOBAL" );
    }
    else {
        nc_inq_varname(  cdfid, var, szVarName );
    }

    for( int l=0; l < nbAttr; l++) {
	
        nc_inq_attname( cdfid, var, l, szAttrName);
        sprintf( szMetaName, "%s#%s", szVarName, szAttrName  );
        nc_inq_att( cdfid, var, szAttrName, &nAttrType, &nAttrLen );
	
        /* Allocate guaranteed minimum size */
        nMetaTempSize = nAttrLen + 1;
        pszMetaTemp = (char *) CPLCalloc( nMetaTempSize, sizeof( char ));
        *pszMetaTemp = '\0';
	
        switch (nAttrType) {
            case NC_CHAR:
                nc_get_att_text( cdfid, var, szAttrName, pszMetaTemp );
                pszMetaTemp[nAttrLen]='\0';
                break;
            case NC_SHORT:
                short *psTemp;
                psTemp = (short *) CPLCalloc( nAttrLen, sizeof( short ) );
                nc_get_att_short( cdfid, var, szAttrName, psTemp );
                for(m=0; m < nAttrLen-1; m++) {
                    sprintf( szTemp, "%hd, ", psTemp[m] );
                    SafeStrcat(&pszMetaTemp, szTemp, &nMetaTempSize);
                }
                sprintf( szTemp, "%hd", psTemp[m] );
                SafeStrcat(&pszMetaTemp, szTemp, &nMetaTempSize);
                CPLFree(psTemp);
                break;
            case NC_INT:
                int *pnTemp;
                pnTemp = (int *) CPLCalloc( nAttrLen, sizeof( int ) );
                nc_get_att_int( cdfid, var, szAttrName, pnTemp );
                for(m=0; m < nAttrLen-1; m++) {
                    sprintf( szTemp, "%d, ", pnTemp[m] );
                    SafeStrcat(&pszMetaTemp, szTemp, &nMetaTempSize);
                }
        	    sprintf( szTemp, "%d", pnTemp[m] );
        	    SafeStrcat(&pszMetaTemp, szTemp, &nMetaTempSize);
                CPLFree(pnTemp);
                break;
            case NC_FLOAT:
                float *pfTemp;
                pfTemp = (float *) CPLCalloc( nAttrLen, sizeof( float ) );
                nc_get_att_float( cdfid, var, szAttrName, pfTemp );
                for(m=0; m < nAttrLen-1; m++) {
                    sprintf( szTemp, "%f, ", pfTemp[m] );
                    SafeStrcat(&pszMetaTemp, szTemp, &nMetaTempSize);
                }
        	    sprintf( szTemp, "%f", pfTemp[m] );
        	    SafeStrcat(&pszMetaTemp,szTemp, &nMetaTempSize);
                CPLFree(pfTemp);
                break;
            case NC_DOUBLE:
                double *pdfTemp;
                pdfTemp = (double *) CPLCalloc(nAttrLen, sizeof(double));
                nc_get_att_double( cdfid, var, szAttrName, pdfTemp );
                for(m=0; m < nAttrLen-1; m++) {
                    sprintf( szTemp, "%.16g, ", pdfTemp[m] );
                    SafeStrcat(&pszMetaTemp, szTemp, &nMetaTempSize);
                }
        	    sprintf( szTemp, "%.16g", pdfTemp[m] );
        	    SafeStrcat(&pszMetaTemp, szTemp, &nMetaTempSize);
                CPLFree(pdfTemp);
                break;
            default:
                break;
        }

        papszMetadata = CSLSetNameValue(papszMetadata, 
                                        szMetaName, 
                                        pszMetaTemp);
        CPLFree(pszMetaTemp);
    }
	

    return CE_None;

}


/************************************************************************/
/*                netCDFDataset::CreateSubDatasetList()                 */
/************************************************************************/
void netCDFDataset::CreateSubDatasetList( )
{

    char         szDim[ MAX_NC_NAME ];
    char         szTemp[ MAX_NC_NAME ];
    char         szType[ MAX_NC_NAME ];
    char         szName[ MAX_NC_NAME ];
    char         szVarStdName[ MAX_NC_NAME ];
    int          nDims;
    int          nVar;
    int          nVarCount;
    int          i;
    nc_type      nVarType;
    int          *ponDimIds;
    size_t       nDimLen;
    int          nSub;
    nc_type      nAttype;
    size_t       nAttlen;

    netCDFDataset 	*poDS;
    poDS = this;

    nSub=1;
    nc_inq_nvars ( cdfid, &nVarCount );
    for ( nVar = 0; nVar < nVarCount; nVar++ ) {

        nc_inq_varndims ( cdfid, nVar, &nDims );
        if( nDims >= 2 ) {
            ponDimIds = (int *) CPLCalloc( nDims, sizeof( int ) );
            nc_inq_vardimid ( cdfid, nVar, ponDimIds );
	    
/* -------------------------------------------------------------------- */
/*      Create Sub dataset list                                         */
/* -------------------------------------------------------------------- */
            szDim[0]='\0';
            for( i = 0; i < nDims; i++ ) {
                nc_inq_dimlen ( cdfid, ponDimIds[i], &nDimLen );
                sprintf(szTemp, "%d", (int) nDimLen);
                strcat(szTemp,  "x" );
                strcat(szDim,   szTemp);
            }

            nc_inq_vartype( cdfid, nVar, &nVarType );
/* -------------------------------------------------------------------- */
/*      Get rid of the last "x" character                               */
/* -------------------------------------------------------------------- */
            szDim[strlen(szDim) - 1] = '\0';
            switch( nVarType ) {
		
                case NC_BYTE:
                case NC_CHAR:
                    strcpy(szType, "8-bit character");
                    break;

                case NC_SHORT: 
                    strcpy(szType, "8-bit integer");
                    break;
                case NC_INT:
                    strcpy(szType, "16-bit integer");
                    break;
                case NC_FLOAT:
                    strcpy(szType, "32-bit floating-point");
                    break;
                case NC_DOUBLE:
                    strcpy(szType, "64-bit floating-point");
                    break;

                default:
                    break;
            }
            nc_inq_varname  ( cdfid, nVar, szName);
            nc_inq_att( cdfid, nVar, "standard_name", &nAttype, &nAttlen);
            if( nc_get_att_text ( cdfid, nVar, "standard_name", 
                                  szVarStdName ) == NC_NOERR ) {
                szVarStdName[nAttlen] = '\0';
            }
            else {
                strcpy( szVarStdName, szName );
            }
    
            sprintf( szTemp, "SUBDATASET_%d_NAME", nSub) ;
	    
            poDS->papszSubDatasets =
                CSLSetNameValue( poDS->papszSubDatasets, szTemp,
                                 CPLSPrintf( "NETCDF:\"%s\":%s",
                                             poDS->osFilename.c_str(),
                                             szName) ) ;
	    
            sprintf(  szTemp, "SUBDATASET_%d_DESC", nSub++ );

            poDS->papszSubDatasets =
                CSLSetNameValue( poDS->papszSubDatasets, szTemp,
                                 CPLSPrintf( "[%s] %s (%s)", 
                                             szDim,
                                             szVarStdName,
                                             szType ) );
            CPLFree(ponDimIds);
        }
    }

}
    
/************************************************************************/
/*                              IdentifyFileType()                      */
/************************************************************************/

int netCDFDataset::IdentifyFileType( GDALOpenInfo * poOpenInfo, bool bCheckHDF5 = TRUE )

{
/* -------------------------------------------------------------------- */
/*      Does this appear to be a netcdf file?                           */
/*      Note: proper care should be done at configure to detect which   */
/*        netcdf versions are supported (nc, nc2, nc4), as does CDO     */
/*      http://www.unidata.ucar.edu/software/netcdf/docs/faq.html#fv1_5 */
/* -------------------------------------------------------------------- */
    /* This should be extended to detect filetype with NETCDF: syntax */
    if( EQUALN(poOpenInfo->pszFilename,"NETCDF:",7) )
        return NCDF_FILETYPE_UNKNOWN;
    if ( poOpenInfo->nHeaderBytes < 4 )
        return NCDF_FILETYPE_NONE;
    if ( EQUALN((char*)poOpenInfo->pabyHeader,"CDF\001",4) )
        return NCDF_FILETYPE_NC;
    else if ( EQUALN((char*)poOpenInfo->pabyHeader,"CDF\002",4) )
        return NCDF_FILETYPE_NC2;
    else if ( EQUALN((char*)poOpenInfo->pabyHeader,"\211HDF\r\n\032\n",8) ) {
        /* Make sure this driver only opens netCDF-4 files and not other HDF5 files */
        /* This check should be relaxed, but there is no clear way to make a difference */
        /* If user really wants to open with this driver, use NETCDF:file.nc:var format */
        if ( TRUE == bCheckHDF5 ) { /* Only check if asked for */
            const char* pszExtension = CPLGetExtension( poOpenInfo->pszFilename );
            if ( ! ( EQUAL( pszExtension, "nc")  || EQUAL( pszExtension, "nc4") ) )
                return NCDF_FILETYPE_HDF5;
        }
        /* Requires netcdf v4, should also test for netCDF-4 support compiled in */
        /* This test could be done at configure like in CDO */
        /* Anyway, if support is not built-in the file will not open, and should fall-back to HDF5 */
        if ( nc_inq_libvers()[0]=='4' ) 
            return NCDF_FILETYPE_NC4;
        else
            return NCDF_FILETYPE_HDF5; 
    }

    return NCDF_FILETYPE_NONE;
} 

/************************************************************************/
/*                              Identify()                              */
/************************************************************************/

int netCDFDataset::Identify( GDALOpenInfo * poOpenInfo )

{
    if( EQUALN(poOpenInfo->pszFilename,"NETCDF:",7) ) {
        return TRUE;
    }
    int nTmpFileType = IdentifyFileType( poOpenInfo );
    if( NCDF_FILETYPE_NONE == nTmpFileType ||
        NCDF_FILETYPE_HDF5 == nTmpFileType ||
        NCDF_FILETYPE_UNKNOWN == nTmpFileType )
        return FALSE;
    else
        return TRUE;
} 

/************************************************************************/
/*                                Open()                                */
/************************************************************************/

GDALDataset *netCDFDataset::Open( GDALOpenInfo * poOpenInfo )
    
{
    int          j;
    unsigned int k;
    int          nd;
    int          cdfid, dim_count, var, var_count;
    int          i = 0;
    size_t       lev_count;
    size_t       nTotLevCount = 1;
    int          nDim = 2;
    int          status;
    int          nDimID;
    char         attname[NC_MAX_NAME];
    int          ndims, nvars, ngatts, unlimdimid;
    int          nCount=0;
    int          nVarID=-1;
    int          nTmpFileType=NCDF_FILETYPE_NONE;

/* -------------------------------------------------------------------- */
/*      Does this appear to be a netcdf file?                           */
/* -------------------------------------------------------------------- */
    if( ! EQUALN(poOpenInfo->pszFilename,"NETCDF:",7) ) {
        nTmpFileType = IdentifyFileType( poOpenInfo );
        /* Note: not calling Identify() directly, because I want to have the file type */
        /* Duplicating the hdf5 test (also in Identify()) */
        if( NCDF_FILETYPE_NONE == nTmpFileType ||
            NCDF_FILETYPE_HDF5 == nTmpFileType ||
            NCDF_FILETYPE_UNKNOWN == nTmpFileType )
            return NULL;
    }

    netCDFDataset 	*poDS;
    poDS = new netCDFDataset();

/* -------------------------------------------------------------------- */
/*      Disable PAM, at least temporarily. See bug #4244                */
/* -------------------------------------------------------------------- */
    poDS->nPamFlags |= GPF_DISABLED;


    poDS->SetDescription( poOpenInfo->pszFilename );
    
/* -------------------------------------------------------------------- */
/*       Check if filename start with NETCDF: tag                       */
/* -------------------------------------------------------------------- */
    if( EQUALN( poOpenInfo->pszFilename,"NETCDF:",7) )
    {
        char **papszName =
            CSLTokenizeString2( poOpenInfo->pszFilename,
                                ":", CSLT_HONOURSTRINGS|CSLT_PRESERVEESCAPES );
        
        /* -------------------------------------------------------------------- */
        /*    Check for drive name in windows NETCDF:"D:\...                    */
        /* -------------------------------------------------------------------- */
        if ( CSLCount(papszName) == 4 &&
             strlen(papszName[1]) == 1 &&
             (papszName[2][0] == '/' || papszName[2][0] == '\\') )
        {
            poDS->osFilename = papszName[1];
            poDS->osFilename += ':';
            poDS->osFilename += papszName[2];
            poDS->osSubdatasetName = papszName[3];
            poDS->bTreatAsSubdataset = TRUE;
            CSLDestroy( papszName );
        }
        else if( CSLCount(papszName) == 3 )
        {
            poDS->osFilename = papszName[1];
            poDS->osSubdatasetName = papszName[2];
            poDS->bTreatAsSubdataset = TRUE;
            CSLDestroy( papszName );
    	}
        else
        {
            CSLDestroy( papszName );
            delete poDS;
            CPLError( CE_Failure, CPLE_AppDefined,
                      "Failed to parse NETCDF: prefix string into expected three fields." );
            return NULL;
        }
        /* Identify filetype from real file, with bCheckHDF5=FALSE */ 
        GDALOpenInfo* poOpenInfo2 = new GDALOpenInfo(poDS->osFilename.c_str(), GA_ReadOnly );
        poDS->nFileType = IdentifyFileType( poOpenInfo2, FALSE );
        delete poOpenInfo2;
        if( NCDF_FILETYPE_NONE == poDS->nFileType ||
            NCDF_FILETYPE_UNKNOWN == poDS->nFileType ) {
            delete poDS;
            return NULL;
        }        
    }
    else 
    {
        poDS->osFilename = poOpenInfo->pszFilename;
        poDS->bTreatAsSubdataset = FALSE;
        poDS->nFileType = nTmpFileType;
    }

/* -------------------------------------------------------------------- */
/*      Try opening the dataset.                                        */
/* -------------------------------------------------------------------- */
    if( nc_open( poDS->osFilename, NC_NOWRITE, &cdfid ) != NC_NOERR ) {
        delete poDS;
        return NULL;
    }

/* -------------------------------------------------------------------- */
/*      Is this a real netCDF file?                                     */
/* -------------------------------------------------------------------- */
    status = nc_inq(cdfid, &ndims, &nvars, &ngatts, &unlimdimid);
    if( status != NC_NOERR ) {
        delete poDS;
        return NULL;
    }
    
/* -------------------------------------------------------------------- */
/*      Confirm the requested access is supported.                      */
/* -------------------------------------------------------------------- */
    if( poOpenInfo->eAccess == GA_Update )
    {
        CPLError( CE_Failure, CPLE_NotSupported, 
                  "The NETCDF driver does not support update access to existing"
                  " datasets.\n" );
        nc_close( cdfid );
        delete poDS;
        return NULL;
    }
    
/* -------------------------------------------------------------------- */
/*      Does the request variable exist?                                */
/* -------------------------------------------------------------------- */
    if( poDS->bTreatAsSubdataset )
    {
        status = nc_inq_varid( cdfid, poDS->osSubdatasetName, &var);
        if( status != NC_NOERR ) {
            CPLError( CE_Warning, CPLE_AppDefined, 
                      "%s is a netCDF file, but %s is not a variable.",
                      poOpenInfo->pszFilename, 
                      poDS->osSubdatasetName.c_str() );
            
            nc_close( cdfid );
            delete poDS;
            return NULL;
        }
    }

    if( nc_inq_ndims( cdfid, &dim_count ) != NC_NOERR || dim_count < 2 )
    {
        CPLError( CE_Warning, CPLE_AppDefined, 
                  "%s is a netCDF file, but not in GMT configuration.",
                  poOpenInfo->pszFilename );

        nc_close( cdfid );
        delete poDS;
        return NULL;
    }

    CPLDebug( "GDAL_netCDF", "dim_count = %d\n", dim_count );

    if( (status = nc_get_att_text( cdfid, NC_GLOBAL, "Conventions",
                                   attname )) != NC_NOERR ) {
        CPLError( CE_Warning, CPLE_AppDefined, 
                  "No UNIDATA NC_GLOBAL:Conventions attribute");
        /* note that 'Conventions' is always capital 'C' in CF spec*/
    }


/* -------------------------------------------------------------------- */
/*      Create band information objects.                                */
/* -------------------------------------------------------------------- */
    if ( nc_inq_nvars ( cdfid, &var_count) != NC_NOERR )
    {
        delete poDS;
        return NULL;
    }    
    
    CPLDebug( "GDAL_netCDF", "var_count = %d\n", var_count );

/* -------------------------------------------------------------------- */
/*      Create a corresponding GDALDataset.                             */
/*      Create Netcdf Subdataset if filename as NETCDF tag              */
/* -------------------------------------------------------------------- */
    poDS->cdfid = cdfid;

    poDS->ReadAttributes( cdfid, NC_GLOBAL );	

/* -------------------------------------------------------------------- */
/*  Verify if only one variable has 2 dimensions                        */
/* -------------------------------------------------------------------- */
    for ( j = 0; j < nvars; j++ ) {

        nc_inq_varndims ( cdfid, j, &ndims );
        if( ndims >= 2 ) {
            nVarID=j;
            nCount++;
        }
    }

/* -------------------------------------------------------------------- */
/*      We have more than one variable with 2 dimensions in the         */
/*      file, then treat this as a subdataset container dataset.        */
/* -------------------------------------------------------------------- */
    if( (nCount > 1) && !poDS->bTreatAsSubdataset )
    {
        poDS->CreateSubDatasetList();
        poDS->SetMetadata( poDS->papszMetadata );
        poDS->TryLoadXML();
        return( poDS );
    }

/* -------------------------------------------------------------------- */
/*      If we are not treating things as a subdataset, then capture     */
/*      the name of the single available variable as the subdataset.    */
/* -------------------------------------------------------------------- */
    if( !poDS->bTreatAsSubdataset ) // nCount must be 1!
    {
        char szVarName[NC_MAX_NAME];

        nc_inq_varname( cdfid, nVarID, szVarName);
        poDS->osSubdatasetName = szVarName;
    }

/* -------------------------------------------------------------------- */
/*      Open the NETCDF subdataset NETCDF:"filename":subdataset         */
/* -------------------------------------------------------------------- */
    var=-1;
    nc_inq_varid( cdfid, poDS->osSubdatasetName, &var);
    nd = 0;
    nc_inq_varndims ( cdfid, var, &nd );

    poDS->paDimIds = (int *)CPLCalloc(nd, sizeof( int ) );
    poDS->panBandDimPos = ( int * ) CPLCalloc( nd, sizeof( int ) );

    nc_inq_vardimid( cdfid, var, poDS->paDimIds );
	
/* -------------------------------------------------------------------- */
/*      Check fi somebody tried to pass a variable with less than 2D    */
/* -------------------------------------------------------------------- */
    if ( nd < 2 ) {
        CPLFree( poDS->paDimIds );
        CPLFree( poDS->panBandDimPos );
        delete poDS;
        return NULL;
    }

/* -------------------------------------------------------------------- */
/*      CF-1 Convention                                                 */
/*      dimensions to appear in the relative order T, then Z, then Y,   */
/*      then X  to the file. All other dimensions should, whenever      */
/*      possible, be placed to the left of the spatiotemporal           */
/*      dimensions.                                                     */
/* -------------------------------------------------------------------- */

/* -------------------------------------------------------------------- */
/*      Get X dimensions information                                    */
/* -------------------------------------------------------------------- */
    poDS->nDimXid = poDS->paDimIds[nd-1];
    nc_inq_dimlen ( cdfid, poDS->nDimXid, &poDS->xdim );
    poDS->nRasterXSize = poDS->xdim;

/* -------------------------------------------------------------------- */
/*      Get Y dimension information                                     */
/* -------------------------------------------------------------------- */
    poDS->nDimYid = poDS->paDimIds[nd-2];
    nc_inq_dimlen ( cdfid, poDS->nDimYid, &poDS->ydim );
    poDS->nRasterYSize = poDS->ydim;


    for( j=0,k=0; j < nd; j++ ){
        if( poDS->paDimIds[j] == poDS->nDimXid ){ 
            poDS->panBandDimPos[0] = j;         // Save Position of XDim
            k++;
        }
        if( poDS->paDimIds[j] == poDS->nDimYid ){
            poDS->panBandDimPos[1] = j;         // Save Position of YDim
            k++;
        }
    }
/* -------------------------------------------------------------------- */
/*      X and Y Dimension Ids were not found!                           */
/* -------------------------------------------------------------------- */
    if( k != 2 ) {
        CPLFree( poDS->paDimIds );
        CPLFree( poDS->panBandDimPos );
        return NULL;
    }
	    
/* -------------------------------------------------------------------- */
/*      Read Metadata for this variable                                 */
/* -------------------------------------------------------------------- */
    poDS->ReadAttributes( cdfid, var );
	
/* -------------------------------------------------------------------- */
/*      Read Metadata for each dimension                                */
/* -------------------------------------------------------------------- */
    
    for( j=0; j < dim_count; j++ ){
        nc_inq_dimname( cdfid, j, poDS->papszDimName[j] );
        status = nc_inq_varid( cdfid, poDS->papszDimName[j], &nDimID );
        if( status == NC_NOERR ) {
            poDS->ReadAttributes( cdfid, nDimID );
        }
    }

    poDS->SetProjection( var );
    poDS->SetMetadata( poDS->papszMetadata );

/* -------------------------------------------------------------------- */
/*      Create bands                                                    */
/* -------------------------------------------------------------------- */
    poDS->panBandZLev = (int *)CPLCalloc( nd-2, sizeof( int ) );
    
    nTotLevCount = 1;
    if ( dim_count > 2 ) {
        nDim=2;
        for( j=0; j < nd; j++ ){
            if( ( poDS->paDimIds[j] != poDS->nDimXid ) && 
                ( poDS->paDimIds[j] != poDS->nDimYid ) ){
                nc_inq_dimlen ( cdfid, poDS->paDimIds[j], &lev_count );
                nTotLevCount *= lev_count;
                poDS->panBandZLev[ nDim-2 ] = lev_count;
                poDS->panBandDimPos[ nDim++ ] = j;  //Save Position of ZDim
            }
        }
    }
    i=0;

    for ( unsigned int lev = 0; lev < nTotLevCount ; lev++ ) {
        char ** papszToken;
        papszToken=NULL;

        netCDFRasterBand *poBand =
            new netCDFRasterBand(poDS, var, nDim, lev,
                                 poDS->panBandZLev, poDS->panBandDimPos, i+1 );

        poDS->SetBand( i+1, poBand );
        i++;
    } 

    CPLFree( poDS->paDimIds );
    CPLFree( poDS->panBandDimPos );
    CPLFree( poDS->panBandZLev );
    
    poDS->nBands = i;

    // Handle angular geographic coordinates here

/* -------------------------------------------------------------------- */
/*      Initialize any PAM information.                                 */
/* -------------------------------------------------------------------- */
    if( poDS->bTreatAsSubdataset )
    {
        poDS->SetPhysicalFilename( poDS->osFilename );
        poDS->SetSubdatasetName( poDS->osSubdatasetName );
    }
    
    poDS->TryLoadXML();

    if( poDS->bTreatAsSubdataset )
        poDS->oOvManager.Initialize( poDS, ":::VIRTUAL:::" );
    else
        poDS->oOvManager.Initialize( poDS, poDS->osFilename );

    return( poDS );
}


/************************************************************************/
/*                            CopyMetadata()                            */
/*                                                                      */
/*      Create a copy of metadata for NC_GLOBAL or a variable           */
/************************************************************************/

void CopyMetadata( void  *poDS, int fpImage, int CDFVarID ) {

    char       **papszMetadata;
    char       **papszFieldData;
    const char *pszField;
    char       szMetaName[ MAX_STR_LEN ];
    char       szMetaValue[ MAX_STR_LEN ];
    char       szTemp[ MAX_STR_LEN ];
    int        nDataLength;
    int        nItems;
    int        bCopyItem;

    if( CDFVarID == NC_GLOBAL ) {
        papszMetadata = GDALGetMetadata( (GDALDataset *) poDS,"");
    } else {
        papszMetadata = GDALGetMetadata( (GDALRasterBandH) poDS, NULL );
    }

    nItems = CSLCount( papszMetadata );             
    
    for(int k=0; k < nItems; k++ ) {
        bCopyItem = TRUE;
        pszField = CSLGetField( papszMetadata, k );
        papszFieldData = CSLTokenizeString2 (pszField, "=", 
                                             CSLT_HONOURSTRINGS );
        if( papszFieldData[1] != NULL ) {
            strcpy( szMetaName,  papszFieldData[ 0 ] );
            strcpy( szMetaValue, papszFieldData[ 1 ] );

            /* Fix various fixes with metadata translation */ 
            if( CDFVarID == NC_GLOBAL ) {
                /* Remove NC_GLOBAL prefix for netcdf global Metadata */ 
                if( strncmp( szMetaName, "NC_GLOBAL#", 10 ) == 0 ) {
                    strcpy( szTemp, szMetaName+10 );
                    strcpy( szMetaName, szTemp );
                } 
                /* GDAL Metadata renamed as GDAL-[meta] */
                else if ( strstr( szMetaName, "#" ) == NULL ) {
                    strcpy( szTemp, "GDAL_" );
                    strcat( szTemp, szMetaName );
                    strcpy( szMetaName, szTemp );
                }
                /* Keep time, lev and depth information for safe-keeping */
                /* Time and vertical coordinate handling need improvements */
                else if( strncmp( szMetaName, "time#", 5 ) == 0 ) {
                    szMetaName[4] = '-';
                }
                else if( strncmp( szMetaName, "lev#", 4 ) == 0 ) {
                    szMetaName[3] = '-';
                }
                else if( strncmp( szMetaName, "depth#", 6 ) == 0 ) {
                    szMetaName[5] = '-';
                }
                /* Only copy data without # (previously all data was copied)  */
                if ( strstr( szMetaName, "#" ) != NULL ) {   
                    bCopyItem = FALSE;
                }
                /* netCDF attributes do not like the '#' character. */
                // for( unsigned int h=0; h < strlen( szMetaName ) -1 ; h++ ) {
                //     if( szMetaName[h] == '#' ) szMetaName[h] = '-'; 
                // }
            }
            else {
                if ( strncmp( szMetaName, "NETCDF_VARNAME", 14) == 0 ) 
                    bCopyItem = FALSE;
            }

            if ( bCopyItem ) {
                nDataLength = strlen( szMetaValue );
                nc_put_att_text( fpImage, 
                                 CDFVarID, 
                                 szMetaName,
                                 nDataLength,
                                 szMetaValue );                
            }
	    
        }
        CSLDestroy( papszFieldData );
    }

}

/************************************************************************/
/*                             CreateCopy()                             */
/************************************************************************/

/*
Driver options:

WRITELONLAT=yes/no (default: yes for geographic, no for projected)
WRITEGDALTAGS=yes/no (default: no)
TYPELONLAT=float/double (default: double for geographic, float for projected)
COMPRESS=NONE/DEFLATE/PACKED (default: NONE)
ZLEVEL=[1-9] (default: 6)
FILETYPE=NC/NC2/NC4 (COMPRESS=DEFLATE sets FILETYPE=NC4)

Processing steps:

1 error checking
2 create dataset
3 def dims and write metadata
4 write projection info
5 write variables: data, projection var, metadata
6 close dataset
7 write pam 
*/

//#include "netcdf-tmp.h"

static GDALDataset*
NCDFCreateCopy2( const char * pszFilename, GDALDataset *poSrcDS, 
                int bStrict, char ** papszOptions, 
                GDALProgressFunc pfnProgress, void * pProgressData )

{
    int  nBands = poSrcDS->GetRasterCount();
    int  nXSize = poSrcDS->GetRasterXSize();
    int  nYSize = poSrcDS->GetRasterYSize();
    int  nLonSize=0, nLatSize=0; 
    int  iBand;

    int  anBandDims[ NC_MAX_DIMS ];
    int  anBandMap[  NC_MAX_DIMS ];

    int  bWriteGeoTransform = FALSE;
    int  bWriteCFLonLat = FALSE;
    char pszNetcdfProjection[ NC_MAX_NAME ];
    int bWriteGDALTags = FALSE;

    const char *pszValue;
    int nFileType, nCompress, nZLevel;

    char   szTemp[ MAX_STR_LEN ];

    if (nBands == 0)
    {
        CPLError( CE_Failure, CPLE_NotSupported, 
                  "NetCDF driver does not support source dataset with zero band.\n");
        return NULL;
    }

    for( iBand=1; iBand <= nBands; iBand++ )
    {
        GDALRasterBand *poSrcBand = poSrcDS->GetRasterBand( iBand );
        GDALDataType eDT = poSrcBand->GetRasterDataType();
        if (eDT == GDT_Unknown || GDALDataTypeIsComplex(eDT))
        {
            CPLError( CE_Failure, CPLE_NotSupported, 
                      "NetCDF driver does not support source dataset with band of complex type.");
            return NULL;
        }
    }

    if( !pfnProgress( 0.0, NULL, pProgressData ) )
        return NULL;

/* -------------------------------------------------------------------- */
/*      Get Projection ref for netCDF data CF-1 Convention              */
/* -------------------------------------------------------------------- */

    OGRSpatialReference oSRS;
    char *pszWKT = (char *) poSrcDS->GetProjectionRef();
    if( pszWKT != NULL )
        oSRS.importFromWkt( &pszWKT );
    char *pszProj4Defn = NULL;
    oSRS.exportToProj4( &pszProj4Defn );

/* -------------------------------------------------------------------- */
/*      Process options.                                                */
/* -------------------------------------------------------------------- */

    /* Filetype and compression */
    nFileType = NCDF_FILETYPE_NC;
    pszValue = CSLFetchNameValue( papszOptions, "FILETYPE" );
    if ( pszValue != NULL ) {
        if ( EQUAL( pszValue, "NC2" ) )
            nFileType = NCDF_FILETYPE_NC2;
        else if ( EQUAL( pszValue, "NC4" ) )
            nFileType = NCDF_FILETYPE_NC4;
        else if ( EQUAL( pszValue, "NC4C" ) )
            CPLError( CE_Warning, CPLE_NotSupported,
                      "WARNING: Filetype NC4C not supported");
            // nFileType = NCDF_FILETYPE_NC4C;
    }

    nCompress = NCDF_COMPRESS_NONE;
    pszValue = CSLFetchNameValue( papszOptions, "COMPRESS" );
    if ( pszValue != NULL ) {
        if ( EQUAL( pszValue, "PACKED" ) ) {
            CPLError(  CE_Warning,  CPLE_NotSupported,
                       "WARNING: PACKED compression not supported yet." );
            // nCompress = NCDF_COMPRESS_PACKED;
        }
        else if ( EQUAL( pszValue, "DEFLATE" ) ) {
            nCompress = NCDF_COMPRESS_DEFLATE;
            if ( nFileType != NCDF_FILETYPE_NC4 ) {
                CPLError( CE_Warning,  CPLE_None,
                          "NOTICE: Filetype set to NC4 because compression is set to DEFLATE." );
                nFileType = NCDF_FILETYPE_NC4;
            }
        }
        else if ( EQUAL( pszValue, "SZIP" ) ) {
            CPLError(  CE_Warning,  CPLE_NotSupported,
                       "WARNING: SZIP compression not supported by netcdf." );
        }
    }

    nZLevel = NCDF_DEFLATE_LEVEL;
    pszValue = CSLFetchNameValue( papszOptions, "ZLEVEL" );
    if( pszValue  != NULL )
    {
        nZLevel =  atoi( pszValue );
        if (!(nZLevel >= 1 && nZLevel <= 9))
        {
            CPLError( CE_Warning, CPLE_IllegalArg, 
                    "ZLEVEL=%s value not recognised, ignoring.",
                    pszValue );
            nZLevel = NCDF_DEFLATE_LEVEL;
        }
    }
    /* Check for proper NC4 support, should also check for DEFLATE support */
    if (  ( nFileType == NCDF_FILETYPE_NC4 ||
            nFileType == NCDF_FILETYPE_NC4C ) &&
          nc_inq_libvers()[0]!='4' ) {
        CPLError( CE_Warning, CPLE_IllegalArg, 
                  "NC4 Filetype not supported, falling back to NC Filetype, without DEFLATE compression.\nPlease install netcdf-4 (using netcdf version %s).",
                  nc_inq_libvers() );
        nFileType = NCDF_FILETYPE_NC;
        if ( nCompress == NCDF_COMPRESS_DEFLATE )
            nCompress = NCDF_COMPRESS_NONE;
    }

    printf("TMP ET options filetype=%d compress=%d zlevel=%d\n",nFileType, nCompress, nZLevel );

/* -------------------------------------------------------------------- */
/*      Create the dataset.                                             */
/* -------------------------------------------------------------------- */

    int fpImage;
    int status;
    int nXDimID = 0;
    int nYDimID = 0;
    int nLonDimID = 0;
    int nLatDimID = 0;
    GDALDataType eDT;
    nc_type eLonLatType = NC_DOUBLE;
    int bBottomUp = FALSE;
    
    // status = nc_create( pszFilename, NC_CLOBBER,  &fpImage );
    int nCMode = NC_CLOBBER;
    if ( nFileType == NCDF_FILETYPE_NC2 )
        nCMode = NC_CLOBBER|NC_64BIT_OFFSET;
    else if ( nFileType == NCDF_FILETYPE_NC4 )
        nCMode = NC_CLOBBER|NC_NETCDF4;
    else if ( nFileType == NCDF_FILETYPE_NC4C )
        nCMode = NC_CLOBBER|NC_NETCDF4|NC_CLASSIC_MODEL;

    status = nc_create( pszFilename, nCMode,  &fpImage );

    if( status != NC_NOERR )
    {
        CPLError( CE_Failure, CPLE_OpenFailed, 
                  "Unable to create netCDF file %s (Error code %d): %s .\n", 
                  pszFilename, status, nc_strerror(status) );
        return NULL;
    }


    if( oSRS.IsProjected() ) 
    {
        nLonSize = nXSize * nYSize;
        nLatSize = nXSize * nYSize;
        status = nc_def_dim( fpImage, NCDF_DIMNAME_X, nXSize, &nXDimID );
        CPLDebug( "GDAL_netCDF", "status nc_def_dim %s = %d\n", NCDF_DIMNAME_X, status );   
        status = nc_def_dim( fpImage, NCDF_DIMNAME_Y, nYSize, &nYDimID );
        CPLDebug( "GDAL_netCDF", "status nc_def_dim %s = %d\n", NCDF_DIMNAME_Y, status );
        anBandDims[0] = nYDimID;
        anBandDims[1] = nXDimID;
        CPLDebug( "GDAL_netCDF", "nYDimID = %d\n", nXDimID );
        CPLDebug( "GDAL_netCDF", "nXDimID = %d\n", nYDimID );
    }
    else 
    {
        nLonSize = nXSize;
        nLatSize = nYSize;
        status = nc_def_dim( fpImage, NCDF_DIMNAME_LON, nLonSize, &nLonDimID );
        CPLDebug( "GDAL_netCDF", "status nc_def_dim %s = %d\n", NCDF_DIMNAME_LON, status );   
        status = nc_def_dim( fpImage, NCDF_DIMNAME_LAT, nLatSize, &nLatDimID );
        CPLDebug( "GDAL_netCDF", "status nc_def_dim %s = %d\n", NCDF_DIMNAME_LAT, status );   
        anBandDims[0] = nLatDimID;
        anBandDims[1] = nLonDimID;
    }

    CPLDebug( "GDAL_netCDF", "nLonDimID = %d\n", nXDimID );
    CPLDebug( "GDAL_netCDF", "nLatDimID = %d\n", nYDimID );
    CPLDebug( "GDAL_netCDF", "nXSize = %d\n", nXSize );
    CPLDebug( "GDAL_netCDF", "nYSize = %d\n", nYSize );


    if( oSRS.IsProjected() ) 
    {
        bBottomUp = TRUE;       /* netcdf standard is bottom-up */
        bWriteGeoTransform = TRUE;
        bWriteCFLonLat = CSLFetchBoolean( papszOptions, "WRITELONLAT", FALSE );
        // bWriteGDALTags = CSLFetchBoolean( papszOptions, "WRITEGDALTAGS", FALSE );
        bWriteGDALTags = CSLFetchBoolean( papszOptions, "WRITEGDALTAGS", TRUE );
        eLonLatType = NC_FLOAT;
        pszValue =  CSLFetchNameValue(papszOptions,"TYPELONLAT");
        if ( pszValue && EQUAL(pszValue, "DOUBLE" ) ) 
            eLonLatType = NC_DOUBLE;
    }
    else 
    { 
        bBottomUp = TRUE;       /* netcdf standard is bottom-up */
        /* files without a Datum will not have a grid_mapping variable and geographic information */
        if ( oSRS.IsGeographic() )  bWriteGeoTransform = TRUE;
        else  bWriteGeoTransform = FALSE;
        bWriteCFLonLat = CSLFetchBoolean( papszOptions, "WRITELONLAT", TRUE );
        bWriteGDALTags = CSLFetchBoolean( papszOptions, "WRITEGDALTAGS", FALSE );
        // bWriteGDALTags = CSLFetchBoolean( papszOptions, "WRITEGDALTAGS", TRUE );
        if ( bWriteGDALTags ) bWriteGeoTransform = TRUE;
            
        eLonLatType = NC_DOUBLE;
        pszValue =  CSLFetchNameValue(papszOptions,"TYPELONLAT");
        if ( pszValue && EQUAL(pszValue, "FLOAT" ) ) 
            eLonLatType = NC_FLOAT;
    }
   
    //    poDstDS->SetGeoTransform( adfGeoTransform );

/* -------------------------------------------------------------------- */
/*      Copy global metadata                                            */
/*      Add Conventions, GDAL info and history                          */
/* -------------------------------------------------------------------- */
    CopyMetadata((void *) poSrcDS, fpImage, NC_GLOBAL );

    nc_put_att_text( fpImage, 
                     NC_GLOBAL, 
                     "Conventions", 
                     strlen(NCDF_CONVENTIONS),
                     NCDF_CONVENTIONS ); 
    
    nc_put_att_text( fpImage, 
                     NC_GLOBAL, 
                     "GDAL", 
                     strlen(NCDF_GDAL),
                     NCDF_GDAL ); 
    
    sprintf( szTemp, "GDAL NCDFCreateCopy( %s, ... )",pszFilename );
    NCDFAddHistory( fpImage, 
                    szTemp, 
                    poSrcDS->GetMetadataItem("NC_GLOBAL#history","") );

    /* Variables needed for both projected and geographic */
    int NCDFVarID=0;
 
    double adfGeoTransform[6];
    char   szGeoTransform[ MAX_STR_LEN ];

    double *padLonVal = NULL;
    double *padLatVal = NULL; /* should use float for projected, save space */
    double dfX0=0.0, dfDX=0.0, dfY0=0.0, dfDY=0.0;
    double dfTemp=0.0;
    size_t *startLon = NULL;
    size_t *countLon = NULL;
    size_t *startLat = NULL;
    size_t *countLat = NULL;


/* -------------------------------------------------------------------- */
/*      Copy GeoTransform array from source                             */
/* -------------------------------------------------------------------- */
    poSrcDS->GetGeoTransform( adfGeoTransform );
    *szGeoTransform = '\0';
    for( int i=0; i<6; i++ ) {
        sprintf( szTemp, "%.16g ",
                 adfGeoTransform[i] );
        strcat( szGeoTransform, szTemp );
    }
    CPLDebug( "GDAL_netCDF", "szGeoTranform = %s", szGeoTransform );

/* -------------------------------------------------------------------- */
/*      Get projection values                                           */
/* -------------------------------------------------------------------- */

    if( oSRS.IsProjected() )
    {
        const OGR_SRSNode *poPROJCS = oSRS.GetAttrNode( "PROJCS" );
        const char  *pszProjection;
        OGRSpatialReference *poLatLonCRS = NULL;
        OGRCoordinateTransformation *poTransform = NULL;
        // const char *pszNetCDFSRS;
        // double dfNN, dfSN=0.0, dfEE=0.0, dfWE=0.0;

        double *padYVal = NULL;
        double *padXVal = NULL;
        size_t startX[1];
        size_t countX[1];
        size_t startY[1];
        size_t countY[1];

        pszProjection = oSRS.GetAttrValue( "PROJECTION" );

/* -------------------------------------------------------------------- */
/*      Write projection attributes                                     */
/* -------------------------------------------------------------------- */

        /* Basic Projection info (grid_mapping and datum) */
        for( int i=0; poNetcdf_SRS_PT[i].GDAL_SRS != NULL; i++ ) {
            if( EQUAL( poNetcdf_SRS_PT[i].GDAL_SRS, pszProjection ) ) {
        // for(int i=0; poNetcdfSRS[i].netCDFSRS != NULL; i++ ) {
        //     if( EQUAL( poNetcdfSRS[i].SRS, pszProjection ) ) {
                CPLDebug( "GDAL_netCDF", "GDAL PROJECTION = %s , NCDF PROJECTION = %s", 
                          // poNetcdfSRS[i].netCDFSRS);
                          poNetcdf_SRS_PT[i].GDAL_SRS, 
                          poNetcdf_SRS_PT[i].NCDF_SRS);
                printf( "GDAL_netCDF GDAL PROJECTION = %s , NCDF PROJECTION = %s\n", 
                          poNetcdf_SRS_PT[i].GDAL_SRS, 
                          poNetcdf_SRS_PT[i].NCDF_SRS);
                // strcpy( pszNetcdfProjection, poNetcdfSRS[i].netCDFSRS );
                strcpy( pszNetcdfProjection, poNetcdf_SRS_PT[i].NCDF_SRS );
                status = nc_def_var( fpImage, 
                                     // poNetcdfSRS[i].netCDFSRS, 
                                     poNetcdf_SRS_PT[i].NCDF_SRS,
                                     NC_CHAR, 
                                     0, NULL, &NCDFVarID );
                break;
            }
        }
        nc_put_att_text( fpImage, 
                         NCDFVarID, 
                         GRD_MAPPING_NAME,
                         strlen( pszNetcdfProjection ),
                         pszNetcdfProjection );
        dfTemp = oSRS.GetSemiMajor();
        nc_put_att_double( fpImage,
                           NCDFVarID, 
                           "semi_major_axis",
                           NC_DOUBLE,
                           1,
                           &dfTemp );
        dfTemp = oSRS.GetInvFlattening();
        nc_put_att_double( fpImage,
                           NCDFVarID, 
                           "inverse_flattening",
                           NC_DOUBLE,
                           1,
                           &dfTemp );

        /* Various projection attributes */
        //PDS: poNetcdfSRS, poPROJCS, pszProjection
        //PDS: Write these out: separate function due to complexity, need to
        // keep in synch with SetProjection function
        NCDFWriteProjAttribs(poPROJCS, pszProjection, fpImage, NCDFVarID);
        /////////////////

        /*  Optional GDAL custom projection tags */
        if ( bWriteGDALTags ) {
            // if ( strlen(pszProj4Defn) > 0 ) {
            //     nc_put_att_text( fpImage, NCDFVarID, "proj4",
            //                      strlen( pszProj4Defn ), pszProj4Defn );
            // }
            pszWKT = (char *) poSrcDS->GetProjectionRef() ;
            nc_put_att_text( fpImage, 
                             NCDFVarID, 
                             NCDF_SPATIAL_REF,
                             strlen( pszWKT ),
                             pszWKT );
            nc_put_att_text( fpImage, 
                             NCDFVarID, 
                             "GeoTransform",
                             strlen( szGeoTransform ),
                             szGeoTransform );
        }

/* -------------------------------------------------------------------- */
/*      CF projection X/Y attributes                                    */
/* -------------------------------------------------------------------- */

        padXVal = (double *) CPLMalloc( nXSize * sizeof( double ) );
        padYVal = (double *) CPLMalloc( nYSize * sizeof( double ) );

        /* Get OGR transform */
        if ( bWriteCFLonLat == TRUE ) {
      
            poLatLonCRS = oSRS.CloneGeogCS();
            if ( poLatLonCRS != NULL )
                poTransform = OGRCreateCoordinateTransformation( &oSRS, poLatLonCRS );
            if( poTransform != NULL )
            {
                // printf("TMP ET got transform\n");
                padLatVal = (double *) CPLMalloc( nLatSize * sizeof( double ) );
                padLonVal = (double *) CPLMalloc( nLonSize * sizeof( double ) );
            }
            else /* if no OGR transform, then don't write CF lon/lat */
                bWriteCFLonLat = FALSE;
        }
/* -------------------------------------------------------------------- */
/*      Get Y values                                                    */
/* -------------------------------------------------------------------- */
        if ( ! bBottomUp )
            dfY0 = adfGeoTransform[3];
        else /* invert latitude values */ 
            dfY0 = adfGeoTransform[3] + ( adfGeoTransform[5] * nYSize );
        dfDY = adfGeoTransform[5];
        
        for( int j=0; j<nYSize; j++ ) {
            /* The data point is centered inside the pixel */
            if ( ! bBottomUp )
                padYVal[j] = dfY0 + (j+0.5)*dfDY ;
                //padLatVal[k] = dfY0 + j*dfDY ;
            else /* invert latitude values */ 
                padYVal[j] = dfY0 - (j+0.5)*dfDY ;
                //padLatVal[k] = dfY0 - j*dfDY ;
            
             if ( bWriteCFLonLat == TRUE ) {
                 for( int i=0; i<nXSize; i++ ) {
                     padLatVal[j*nXSize+i] = padYVal[j];
                 }
             }
        }
        startX[0] = 0;
        countX[0] = nXSize;
        if ( bWriteCFLonLat == TRUE ) {
            startLat = (size_t *) CPLMalloc( 2 * sizeof( size_t ) );
            countLat = (size_t *) CPLMalloc( 2 * sizeof( size_t ) );
            startLat[0] = 0;
            startLat[1] = 0;
            countLat[0] = nYSize;
            countLat[1] = nXSize;
        }
/* -------------------------------------------------------------------- */
/*      Get X values                                                    */
/* -------------------------------------------------------------------- */
        dfX0 = adfGeoTransform[0];
        dfDX = adfGeoTransform[1];

        for( int i=0; i<nXSize; i++ ) {
            /* The data point is centered inside the pixel */
            padXVal[i] = dfX0 + (i+0.5)*dfDX ;
            if ( bWriteCFLonLat == TRUE ) {
                for( int j=0; j<nYSize; j++ ) {
                    padLonVal[j*nXSize+i] = padXVal[i];
                    // padLonVal[k] = dfX0 + i*dfDX ;
                }
            }
        }
        startY[0] = 0;
        countY[0] = nYSize;
        if ( bWriteCFLonLat == TRUE ) {
            startLon = (size_t *) CPLMalloc( 2 * sizeof( size_t ) );
            countLon = (size_t *) CPLMalloc( 2 * sizeof( size_t ) );
            startLon[0] = 0;
            startLon[1] = 0;
            countLon[0] = nYSize;
            countLon[1] = nXSize;
        }
/* -------------------------------------------------------------------- */
/*      Transform (X,Y) values to (lon,lat)                             */
/* -------------------------------------------------------------------- */

        // for( int i=0; i<nXSize*nYSize; i++ ) {
        //     printf("%f ",padLonVal[i]);
        // }
        if ( bWriteCFLonLat == TRUE ) {
            if( ! poTransform->Transform( nXSize * nYSize, padLonVal, padLatVal, NULL ) ) {
                CPLError( CE_Failure, CPLE_AppDefined, 
                          "Unable to Transform (X,Y) to (lon,lat).\n" );
            }
        }

        /* Free the srs and transform objects */
        if ( poLatLonCRS != NULL ) CPLFree( poLatLonCRS );
        if ( poTransform != NULL ) CPLFree( poTransform );


/* -------------------------------------------------------------------- */
/*      Write X attributes                                              */
/* -------------------------------------------------------------------- */
        int anXDims[1];
        anXDims[0] = nXDimID;
        status = nc_def_var( fpImage, NCDF_DIMNAME_X, NC_DOUBLE, 1, anXDims, &NCDFVarID );
        nc_put_att_text( fpImage,
                         NCDFVarID,
                         "standard_name",
                         strlen("projection_x_coordinate"),
                         "projection_x_coordinate" );
        nc_put_att_text( fpImage,
                         NCDFVarID,
                         "long_name",
                         strlen("x coordinate of projection"),
                         "x coordinate of projection" );
        nc_put_att_text( fpImage,
                         NCDFVarID,
                         "units",
                         1,
                         "m" ); /*verify this */

/* -------------------------------------------------------------------- */
/*      Write X values                                                  */
/* -------------------------------------------------------------------- */

        /* Temporarily switch to data mode and write data */
        status = nc_enddef( fpImage );
        CPLDebug("GDAL_netCDF", "Writing X values" );
        status = nc_put_vara_double( fpImage, NCDFVarID, startX,
                                     countX, padXVal);
        status = nc_redef( fpImage );
        
        /* free values */
        CPLFree( padXVal );

/* -------------------------------------------------------------------- */
/*      Write Y attributes                                              */
/* -------------------------------------------------------------------- */
        int anYDims[1];
        anYDims[0] = nYDimID;
        status = nc_def_var( fpImage, NCDF_DIMNAME_Y, NC_DOUBLE, 1, anYDims, &NCDFVarID );
        nc_put_att_text( fpImage,
                         NCDFVarID,
                         "standard_name",
                         strlen("projection_y_coordinate"),
                         "projection_y_coordinate" );
        nc_put_att_text( fpImage,
                         NCDFVarID,
                         "long_name",
                         strlen("y coordinate of projection"),
                         "y coordinate of projection" );
        nc_put_att_text( fpImage,
                         NCDFVarID,
                         "units",
                         1,
                         "m" ); /*verify this */

/* -------------------------------------------------------------------- */
/*      Write Y values                                         */
/* -------------------------------------------------------------------- */

        /* Temporarily switch to data mode and write data */
        status = nc_enddef( fpImage );
        CPLDebug("GDAL_netCDF", "Writing Y values" );
        status = nc_put_vara_double( fpImage, NCDFVarID, startY,
                                     countY, padYVal);
        status = nc_redef( fpImage );
        
        /* free values */
        CPLFree( padYVal );

    } // projected

    else  {  /* If not Projected assume Geographic to catch grids without Datum */
        	
/* -------------------------------------------------------------------- */
/*      Write CF-1.x compliant Geographics attributes                   */
/*      Note: WKT information will not be preserved (e.g. WGS84)        */
/*      Don't write custom GDAL values any more as they are not CF-1.x  */
/*      compliant, this is still under discussion                       */
/* -------------------------------------------------------------------- */
        
        if( bWriteGeoTransform == TRUE ) 
 	    {
            strcpy( pszNetcdfProjection, "crs" );
            nc_def_var( fpImage, 
                        pszNetcdfProjection, 
                        NC_CHAR, 
                        0, NULL, &NCDFVarID );
            nc_put_att_text( fpImage, 
                             NCDFVarID, 
                             GRD_MAPPING_NAME,
                             strlen("latitude_longitude"),
                             "latitude_longitude" );
            dfTemp = oSRS.GetPrimeMeridian();
            nc_put_att_double( fpImage,
                               NCDFVarID, 
                               "longitude_of_prime_meridian",
                               NC_DOUBLE,
                               1,
                               &dfTemp );
            dfTemp = oSRS.GetSemiMajor();
            nc_put_att_double( fpImage,
                               NCDFVarID, 
                               "semi_major_axis",
                               NC_DOUBLE,
                               1,
                               &dfTemp );
            dfTemp = oSRS.GetInvFlattening();
            nc_put_att_double( fpImage,
                               NCDFVarID, 
                               "inverse_flattening",
                               NC_DOUBLE,
                               1,
                               &dfTemp );
            /*  Optional GDAL custom projection tags */
            /* ET - this should be written in a common block with projected */
            if ( bWriteGDALTags ) {
            // if ( strlen(pszProj4Defn) > 0 ) {
            //     nc_put_att_text( fpImage, NCDFVarID, "proj4",
            //                      strlen( pszProj4Defn ), pszProj4Defn );
            // }
            pszWKT = (char *) poSrcDS->GetProjectionRef() ;
            nc_put_att_text( fpImage, 
                             NCDFVarID, 
                             NCDF_SPATIAL_REF,
                             strlen( pszWKT ),
                             pszWKT );
            nc_put_att_text( fpImage, 
                             NCDFVarID, 
                             "GeoTransform",
                             strlen( szGeoTransform ),
                             szGeoTransform );
            }
        }


/* -------------------------------------------------------------------- */
/*      Get latitude values                                             */
/* -------------------------------------------------------------------- */
        if ( ! bBottomUp )
            dfY0 = adfGeoTransform[3];
        else /* invert latitude values */ 
            dfY0 = adfGeoTransform[3] + ( adfGeoTransform[5] * nYSize );
        dfDY = adfGeoTransform[5];
        
        padLatVal = (double *) CPLMalloc( nYSize * sizeof( double ) );
        for( int i=0; i<nYSize; i++ ) {
            /* The data point is centered inside the pixel */
            if ( ! bBottomUp )
                padLatVal[i] = dfY0 + (i+0.5)*dfDY ;
            else /* invert latitude values */ 
                padLatVal[i] = dfY0 - (i+0.5)*dfDY ;
        }
        
        startLat = (size_t *) CPLMalloc( sizeof( size_t ) );
        countLat = (size_t *) CPLMalloc( sizeof( size_t ) );
        startLat[0] = 0;
        countLat[0] = nYSize;
                
/* -------------------------------------------------------------------- */
/*      Get longitude values                                            */
/* -------------------------------------------------------------------- */
        dfX0 = adfGeoTransform[0];
        dfDX = adfGeoTransform[1];
        
        padLonVal = (double *) CPLMalloc( nXSize * sizeof( double ) );
        for( int i=0; i<nXSize; i++ ) {
            /* The data point is centered inside the pixel */
            padLonVal[i] = dfX0 + (i+0.5)*dfDX ;
        }
        
        startLon = (size_t *) CPLMalloc( sizeof( size_t ) );
        countLon = (size_t *) CPLMalloc( sizeof( size_t ) );
        startLon[0] = 0;
        countLon[0] = nXSize;
        
    }// not projected
    

/* -------------------------------------------------------------------- */
/*      Write CF projection lat/lon attributes                          */
/* -------------------------------------------------------------------- */
    if ( bWriteCFLonLat ) {

/* -------------------------------------------------------------------- */
/*      Write latitude attributes                                     */
/* -------------------------------------------------------------------- */
        if ( oSRS.IsProjected() ) {
            int anLatDims[2];
            anLatDims[0] = nYDimID;
            anLatDims[1] = nXDimID;
            status = nc_def_var( fpImage, "lat", eLonLatType, 
                                 2, anLatDims, &NCDFVarID );
        }
        else {
            int anLatDims[1];
            anLatDims[0] = nLatDimID;
            status = nc_def_var( fpImage, "lat", eLonLatType, 
                                 1, anLatDims, &NCDFVarID );                  
        }
        status=nc_put_att_text( fpImage,
                         NCDFVarID,
                         "standard_name",
                         8,
                         "latitude" );
        status = nc_put_att_text( fpImage,
                         NCDFVarID,
                         "long_name",
                         8,
                         "latitude" );
        status = nc_put_att_text( fpImage,
                         NCDFVarID,
                         "units",
                         13,
                         "degrees_north" );

/* -------------------------------------------------------------------- */
/*      Write latitude values                                         */
/* -------------------------------------------------------------------- */

        /* Temporarily switch to data mode and write data */
        status = nc_enddef( fpImage );
        CPLDebug("GDAL_netCDF", "Writing lat values" );
        status = nc_put_vara_double( fpImage, NCDFVarID, startLat,
                                     countLat, padLatVal);

        status = nc_redef( fpImage );
        
        /* free values */
        CPLFree( padLatVal );  
        CPLFree( startLat );
        CPLFree( countLat );
        
/* -------------------------------------------------------------------- */
/*      Write longitude attributes                                    */
/* -------------------------------------------------------------------- */
        if ( oSRS.IsProjected() ) {
            int anLonDims[2];
            anLonDims[0] = nYDimID;
            anLonDims[1] = nXDimID;
            status = nc_def_var( fpImage, "lon", eLonLatType, 
                                 2, anLonDims, &NCDFVarID );
        }
        else {
            int anLonDims[1];
            anLonDims[0] = nLonDimID;
            status = nc_def_var( fpImage, "lon", eLonLatType, 
                                 1, anLonDims, &NCDFVarID );
        }
        nc_put_att_text( fpImage,
                         NCDFVarID,
                         "standard_name",
                         9,
                         "longitude" );
        nc_put_att_text( fpImage,
                         NCDFVarID,
                         "long_name",
                         9,
                         "longitude" );
        nc_put_att_text( fpImage,
                         NCDFVarID,
                         "units",
                         12,
                         "degrees_east" );
        
/* -------------------------------------------------------------------- */
/*      Write longitude values                                        */	
/* -------------------------------------------------------------------- */
        
        /* Temporarily switch to data mode and write data */
        status = nc_enddef( fpImage );
        CPLDebug("GDAL_netCDF", "Writing lon values" );
        status = nc_put_vara_double( fpImage, NCDFVarID, startLon,
                                    countLon, padLonVal);
        status = nc_redef( fpImage );
        
        /* free values */
        CPLFree( padLonVal );  
        CPLFree( startLon );
        CPLFree( countLon );
 
    } // bWriteCFLonLat

/* -------------------------------------------------------------------- */
/*      Initialize Band Map                                             */
/* -------------------------------------------------------------------- */

    for(int j=1; j <= nBands; j++ ) {
        anBandMap[j-1]=j;
    }
    
/* -------------------------------------------------------------------- */
/*      Create netCDF variable                                          */
/* -------------------------------------------------------------------- */

    for( int i=1; i <= nBands; i++ ) {

        char      szBandName[ NC_MAX_NAME ];
        GInt16    *pasScanline  = NULL;
        GInt32    *panScanline  = NULL;
        float     *pafScanline  = NULL;
        double    *padScanline  = NULL;
        int       NCDFVarID;
        size_t    start[ GDALNBDIM ];
        size_t    count[ GDALNBDIM ];
        double    dfNoDataValue;
        unsigned char      cNoDataValue;
        float     fNoDataValue;
        int       nlNoDataValue;
        short     nsNoDataValue;
        GDALRasterBandH	hBand;
        const char *tmpMetadata;
        char      szLongName[ NC_MAX_NAME ];

        int nDataType = NC_NAT;
        GDALRasterBand *poSrcBand = poSrcDS->GetRasterBand( i );
        hBand = GDALGetRasterBand( poSrcDS, i );

        /* Get var name from NETCDF_VARNAME */
        tmpMetadata = poSrcBand->GetMetadataItem("NETCDF_VARNAME");
       	if( tmpMetadata != NULL) {
            if(nBands>1) sprintf(szBandName,"%s%d",tmpMetadata,i);
            else strcpy( szBandName, tmpMetadata );
            // poSrcBand->SetMetadataItem("NETCDF_VARNAME","");
        }
        else 
            sprintf( szBandName, "Band%d", i );

        /* Get long_name from <var>#long_name */
        sprintf(szLongName,"%s#long_name",poSrcBand->GetMetadataItem("NETCDF_VARNAME"));
        tmpMetadata = poSrcDS->GetMetadataItem(szLongName);
        if( tmpMetadata != NULL) 
            strcpy( szLongName, tmpMetadata);
        else 
            sprintf( szLongName, "GDAL Band Number %d", i); 

        CPLDebug("GDAL_netCDF", "Writing Band #%d - %s", i, szLongName );

/* -------------------------------------------------------------------- */
/*      Get Data type                                                   */
/* -------------------------------------------------------------------- */
 
        eDT = poSrcDS->GetRasterBand(i)->GetRasterDataType();
        CPLErr      eErr = CE_None;

        dfNoDataValue = poSrcBand->GetNoDataValue(0);
        printf("TMP ET nodata=[%f]\n",dfNoDataValue);

        if( eDT == GDT_Byte ) {
            CPLDebug( "GDAL_netCDF", "%s = GDT_Byte ", szBandName );

            if ( nFileType == NCDF_FILETYPE_NC4 ||
                 nFileType == NCDF_FILETYPE_NC4C )
                nDataType = NC_UBYTE;
            else nDataType = NC_BYTE;
            status = nc_def_var( fpImage, szBandName, nDataType, 
                                 GDALNBDIM, anBandDims, &NCDFVarID );
            // status = nc_def_var( fpImage, szBandName, NC_CHAR, 
            //                      GDALNBDIM, anBandDims, &NCDFVarID );

/* -------------------------------------------------------------------- */
/*      Write Fill Value                                                */
/* -------------------------------------------------------------------- */
                 // cNoDataValue=(unsigned char) dfNoDataValue;
                 cNoDataValue=(unsigned char) dfNoDataValue;
                status=nc_put_att_uchar( fpImage,
                                  NCDFVarID,
                                  _FillValue,
                                  nDataType,
                                  1,
                                  &cNoDataValue );
                // status=nc_put_att_schar( fpImage,
                //                   NCDFVarID,
                //                   _FillValue,
                //                   // NC_CHAR,
                //                   NC_BYTE,
                //                   1,
                //                   &cNoDataValue );
                
                printf("TMP ET nodata=[%d]-[%f]\n",cNoDataValue,dfNoDataValue);
			   if (status != NC_NOERR) 
                   fprintf(stderr, "netcdf error #%d - %s\n", status,nc_strerror(status));
            printf("def_var\n");
            if ( nCompress == NCDF_COMPRESS_DEFLATE ) {
                size_t chunksize[2] = {nXSize,1};
                status = nc_def_var_deflate(fpImage, NCDFVarID,1,1,nZLevel);
                if (status != NC_NOERR) 
                    fprintf(stderr, "%s\n", nc_strerror(status));
else printf("compress ok\n");
                // status = nc_def_var_chunking( fpImage, NCDFVarID, 
                //                               NC_CHUNKED, chunksize );       
                status = nc_def_var_chunking( fpImage, NCDFVarID, 
                                              NC_CHUNKED, NULL );       
         if (status != NC_NOERR) 
                    fprintf(stderr, "%s\n", nc_strerror(status));
else printf("chunking ok\n");
            }
            printf("def_var_deflate\n");

/* -------------------------------------------------------------------- */
/*      Put NetCDF file in data mode.                                   */
/* -------------------------------------------------------------------- */
                status = nc_enddef( fpImage );

/* -------------------------------------------------------------------- */
/*      Write data line per line                                        */
/* -------------------------------------------------------------------- */

             GByte * pabScanline = (GByte *) CPLMalloc( nBands * nXSize );
                // GByte *pabScanline  = (GByte *) CPLMalloc( nBands * nXSize * nYSize );
            // pabScanline2 = (GByte *) CPLMalloc( nBands * nXSize * nYSize );
            for( int iLine = 0; iLine < nYSize && eErr == CE_None; iLine++ )  {
            // for( int iLine = 0; iLine < nYSize && eErr == CE_None; iLine++ )  {
                if (iLine>0 && (iLine % (nYSize/10)) == 0)
                    printf("writing line %d of %d\n",iLine,nYSize);

/* -------------------------------------------------------------------- */
/*      Read data from band i                                           */
/* -------------------------------------------------------------------- */
                eErr = poSrcBand->RasterIO( GF_Read, 0, iLine, nXSize, 1, 
                                            pabScanline, nXSize, 1, GDT_Byte,
                                            0,0);

                // for (int j=0;j<nXSize;j++) if ( pabScanline[j] >126 || pabScanline[j] <0 ) printf("%d ",pabScanline[j]);
           // int startY = iLine;
            // if (bBottomUp)
            //     startY = nYSize - iLine - 1;
            //     eErr = poSrcBand->RasterIO( GF_Read, 0, startY, nXSize, 1, 
            //                                 pabScanline + ( iLine * nXSize * sizeof(GByte ) ),
            //                                 nXSize, 1, GDT_Byte,
            //                                 0,0);
                // eErr = poSrcBand->RasterIO( GF_Read, 0, 0, nXSize, nYSize, 
                //                             pabScanline, nXSize, nYSize, GDT_Byte,
                //                             0,0);

/* -------------------------------------------------------------------- */
/*      Write Data from Band i                                          */
/* -------------------------------------------------------------------- */
                if ( ! bBottomUp )
                    start[0]=iLine;
                else /* invert latitude values */
                    start[0]=nYSize - iLine - 1;
                start[1]=0;
                count[0]=1;
                count[1]=nXSize;

            //     start[0]=0;
            //     // if ( ! bBottomUp )
            //     //     start[0]=iLine;
            //     // else /* invert latitude values */
            //     //     start[0]=nYSize - iLine - 1;
            //     start[1]=0;
            //     count[0]=nYSize;
            //     count[1]=nXSize;
            // printf("start put var\n");
                status = nc_put_vara_uchar (fpImage, NCDFVarID, start,
                                            count, pabScanline);
                if (status != NC_NOERR) 
                    fprintf(stdout, "%s\n", nc_strerror(status));
                    // fprintf(stderr, "%s\n", nc_strerror(status));

            // printf("end put var\n");
                
             }

/* -------------------------------------------------------------------- */
/*      Put NetCDF file back in define mode.                            */
/* -------------------------------------------------------------------- */
                status = nc_redef( fpImage );
		
            CPLFree( pabScanline );
/* -------------------------------------------------------------------- */
/*      Int16                                                           */
/* -------------------------------------------------------------------- */

        } else if( ( eDT == GDT_UInt16 ) || ( eDT == GDT_Int16 ) ) {
            CPLDebug( "GDAL_netCDF", "%s = GDT_Int16 ",szBandName );

            status = nc_def_var( fpImage, szBandName, NC_SHORT, 
                                 GDALNBDIM, anBandDims, &NCDFVarID );

            if ( nCompress == NCDF_COMPRESS_DEFLATE )
                status = nc_def_var_deflate(fpImage, NCDFVarID,1,1,nZLevel);

            pasScanline = (GInt16 *) CPLMalloc( nBands * nXSize *
                                                sizeof( GInt16 ) );
/* -------------------------------------------------------------------- */
/*      Write Fill Value                                                */
/* -------------------------------------------------------------------- */
            nsNoDataValue= (GInt16) dfNoDataValue;
            nc_put_att_short( fpImage,
                              NCDFVarID,
                              _FillValue,
                              NC_SHORT,
                              1,
                              &nsNoDataValue );
            for( int iLine = 0; iLine < nYSize && eErr == CE_None; iLine++ )  {

                eErr = poSrcBand->RasterIO( GF_Read, 0, iLine, nXSize, 1, 
                                            pasScanline, nXSize, 1, GDT_Int16,
                                            0,0);

                if ( ! bBottomUp )
                    start[0]=iLine;
                else /* invert latitude values */
                    start[0]=nYSize - iLine - 1;
                start[1]=0;
                count[0]=1;
                count[1]=nXSize;


                status = nc_enddef( fpImage );
                status = nc_put_vara_short( fpImage, NCDFVarID, start,
                                            count, pasScanline);
                status = nc_redef( fpImage );
            }
            CPLFree( pasScanline );
/* -------------------------------------------------------------------- */
/*      Int32                                                           */
/* -------------------------------------------------------------------- */

        } else if( (eDT == GDT_UInt32) || (eDT == GDT_Int32) ) {
            CPLDebug( "GDAL_netCDF", "%s = GDT_Int32 ",szBandName );

            status = nc_def_var( fpImage, szBandName, NC_INT, 
                                 GDALNBDIM, anBandDims, &NCDFVarID );

            if ( nCompress == NCDF_COMPRESS_DEFLATE )
                status = nc_def_var_deflate(fpImage, NCDFVarID,1,1,nZLevel);

            panScanline = (GInt32 *) CPLMalloc( nBands * nXSize *
                                                sizeof( GInt32 ) );
/* -------------------------------------------------------------------- */
/*      Write Fill Value                                                */
/* -------------------------------------------------------------------- */
            nlNoDataValue= (GInt32) dfNoDataValue;

            nc_put_att_int( fpImage,
                            NCDFVarID,
                            _FillValue,
                            NC_INT,
                            1,
                            &nlNoDataValue );


            for( int iLine = 0; iLine < nYSize && eErr == CE_None; iLine++ )  {

                eErr = poSrcBand->RasterIO( GF_Read, 0, iLine, nXSize, 1, 
                                            panScanline, nXSize, 1, GDT_Int32,
                                            0,0);

                if ( ! bBottomUp )
                    start[0]=iLine;
                else /* invert latitude values */
                    start[0]=nYSize - iLine - 1;
                start[1]=0;
                count[0]=1;
                count[1]=nXSize;


                status = nc_enddef( fpImage );
                status = nc_put_vara_int( fpImage, NCDFVarID, start,
                                          count, panScanline);
                status = nc_redef( fpImage );
            }
            CPLFree( panScanline );
/* -------------------------------------------------------------------- */
/*      float                                                           */
/* -------------------------------------------------------------------- */
        } else if( (eDT == GDT_Float32) ) {
            CPLDebug( "GDAL_netCDF", "%s = GDT_Float32 ",szBandName );

            status = nc_def_var( fpImage, szBandName, NC_FLOAT, 
                                 GDALNBDIM, anBandDims, &NCDFVarID );

            if ( nCompress == NCDF_COMPRESS_DEFLATE )
                status = nc_def_var_deflate(fpImage, NCDFVarID,1,1,nZLevel);

            pafScanline = (float *) CPLMalloc( nBands * nXSize *
                                               sizeof( float ) );

/* -------------------------------------------------------------------- */
/*      Write Fill Value                                                */
/* -------------------------------------------------------------------- */
            fNoDataValue= (float) dfNoDataValue;
            nc_put_att_float( fpImage,
                              NCDFVarID,
                              _FillValue,
                              NC_FLOAT,
                              1,
                              &fNoDataValue );
			   
            for( int iLine = 0; iLine < nYSize && eErr == CE_None; iLine++ )  {

                eErr = poSrcBand->RasterIO( GF_Read, 0, iLine, nXSize, 1, 
                                            pafScanline, nXSize, 1, 
                                            GDT_Float32,
                                            0,0);

                if ( ! bBottomUp )
                    start[0]=iLine;
                else /* invert latitude values */
                    start[0]=nYSize - iLine - 1;
                start[1]=0;
                count[0]=1;
                count[1]=nXSize;


                status = nc_enddef( fpImage );
                status = nc_put_vara_float( fpImage, NCDFVarID, start,
                                            count, pafScanline);
                status = nc_redef( fpImage );
            }
            CPLFree( pafScanline );
/* -------------------------------------------------------------------- */
/*      double                                                          */
/* -------------------------------------------------------------------- */
        } else if( (eDT == GDT_Float64) ) {
            CPLDebug( "GDAL_netCDF", "%s = GDT_Float64 ",szBandName );

            status = nc_def_var( fpImage, szBandName, NC_DOUBLE, 
                                 GDALNBDIM, anBandDims, &NCDFVarID );

            if ( nCompress == NCDF_COMPRESS_DEFLATE )
                status = nc_def_var_deflate(fpImage, NCDFVarID,1,1,nZLevel);

            padScanline = (double *) CPLMalloc( nBands * nXSize *
                                                sizeof( double ) );

/* -------------------------------------------------------------------- */
/*      Write Fill Value                                                */
/* -------------------------------------------------------------------- */
		
            nc_put_att_double( fpImage,
                               NCDFVarID,
                               _FillValue,
                               NC_DOUBLE,
                               1,
                               &dfNoDataValue );

            for( int iLine = 0; iLine < nYSize && eErr == CE_None; iLine++ )  {

                eErr = poSrcBand->RasterIO( GF_Read, 0, iLine, nXSize, 1, 
                                            padScanline, nXSize, 1, 
                                            GDT_Float64,
                                            0,0);

                if ( ! bBottomUp )
                    start[0]=iLine;
                else /* invert latitude values */
                    start[0]=nYSize - iLine - 1;
                start[1]=0;
                count[0]=1;
                count[1]=nXSize;


                status = nc_enddef( fpImage );
                status = nc_put_vara_double( fpImage, NCDFVarID, start,
                                             count, padScanline);
                status = nc_redef( fpImage );
            }
            CPLFree( padScanline );
        }
	
/* -------------------------------------------------------------------- */
/*      Copy Metadata for band                                          */
/* -------------------------------------------------------------------- */

        nc_put_att_text( fpImage, 
                         NCDFVarID, 
                         "long_name", 
                         strlen( szLongName ),
                         szLongName );

        CopyMetadata( (void *) hBand, fpImage, NCDFVarID );

/* -------------------------------------------------------------------- */
/*      Write Projection for band                                       */
/* -------------------------------------------------------------------- */
        if( bWriteGeoTransform == TRUE ) {
            /*	    nc_put_att_text( fpImage, NCDFVarID, 
                    COORDINATES,
                    7,
                    LONLAT );
            */
            // printf("TMP ET writting proj %s %s\n",GRD_MAPPING,pszNetcdfProjection);
            nc_put_att_text( fpImage, NCDFVarID, 
                             GRD_MAPPING,
                             strlen( pszNetcdfProjection ),
                             pszNetcdfProjection );
            if ( bWriteCFLonLat ) {
                nc_put_att_text( fpImage, NCDFVarID, 
                                 COORDINATES,
                                 strlen( LONLAT ),
                                 LONLAT );
            }           
        }

    }




/* -------------------------------------------------------------------- */
/*      Cleanup and close.                                              */
/* -------------------------------------------------------------------- */
//    CPLFree( pabScanline );
nc_close( fpImage );
CPLFree(pszProj4Defn );

/* -------------------------------------------------------------------- */
/*      Re-open dataset, and copy any auxilary pam information.         */
/* -------------------------------------------------------------------- */
    netCDFDataset *poDS = (netCDFDataset *) GDALOpen( pszFilename, GA_ReadOnly );

    if( poDS )
        poDS->CloneInfo( poSrcDS, GCIF_PAM_DEFAULT );

    return poDS;
}

/* PDS */
/* Write any needed projection attributes *
 * poPROJCS: ptr to proj crd system
 * pszProjection: name of projection system in GDAL WKT
 * fpImage: open NetCDF file in writing mode
 * NCDFVarID: NetCDF Var Id of proj system we're writing in to
 *
 * Note: Default behaviour for each projection is just to choose the
 *  appropriate mappings, then call std function that knows how to
 *  write these.
 *
 * However, flexibility exists so that custom code can be written e.g.
 *  if any params need pre-processing in future to transform (i.e.
 *  the straight mappings approach is too simplistic).
 */
void NCDFWriteProjAttribs(const OGR_SRSNode *poPROJCS,
                            const char* pszProjection,
                            const int fpImage, const int NCDFVarID) 
{
    const char *pszParamStr, *pszParamVal;
    double dfStdP[2];
    int bFoundStdP1=FALSE,bFoundStdP2=FALSE;
    double dfTemp=0.0;

    if( EQUAL(pszProjection, SRS_PT_ALBERS_CONIC_EQUAL_AREA) ) {
        NCDFWriteProjAttribsFromMappings(poPROJCS, poAEAMappings, fpImage,
            NCDFVarID);
    }
    else if( EQUAL(pszProjection, SRS_PT_AZIMUTHAL_EQUIDISTANT) ) {
        NCDFWriteProjAttribsFromMappings(poPROJCS, poLAZEQMappings, fpImage,
            NCDFVarID);
    }
    else if( EQUAL(pszProjection, SRS_PT_LAMBERT_AZIMUTHAL_EQUAL_AREA ) ) {
        NCDFWriteProjAttribsFromMappings(poPROJCS, poLAZEQMappings, fpImage,
            NCDFVarID);
    }
    else if( EQUAL(pszProjection, SRS_PT_LAMBERT_CONFORMAL_CONIC_1SP ) ) {
        NCDFWriteProjAttribsFromMappings(poPROJCS, poLC_1SPMappings, fpImage,
            NCDFVarID);
    }
    else if( EQUAL(pszProjection, SRS_PT_LAMBERT_CONFORMAL_CONIC_2SP ) ) {
        NCDFWriteProjAttribsFromMappings(poPROJCS, poLC_2SPMappings, fpImage,
            NCDFVarID);
    }
    else if( EQUAL(pszProjection, SRS_PT_CYLINDRICAL_EQUAL_AREA ) ) {
        NCDFWriteProjAttribsFromMappings(poPROJCS, poLCEAMappings, fpImage,
            NCDFVarID);
    }
    else if( EQUAL(pszProjection, SRS_PT_MERCATOR_1SP) ) {
        NCDFWriteProjAttribsFromMappings(poPROJCS, poM_1SPMappings, fpImage,
            NCDFVarID);
    }        
    else if( EQUAL(pszProjection, SRS_PT_MERCATOR_2SP) ) {
        NCDFWriteProjAttribsFromMappings(poPROJCS, poM_2SPMappings, fpImage,
            NCDFVarID);
    }
    else if( EQUAL(pszProjection, SRS_PT_ORTHOGRAPHIC) ) {
        NCDFWriteProjAttribsFromMappings(poPROJCS, poOrthoMappings, fpImage,
            NCDFVarID);
    }
    else if( EQUAL(pszProjection, SRS_PT_POLAR_STEREOGRAPHIC) ) {
        // TODO: CF-1 says 'standard_parallel' may replace scale factor
        // TODO: In CF-1, LAT_PROJ_ORIGIN is either +90 or -90. Not sure
        //   how this maps to OGC stuff?
        //  eg do we have to check LAT_PROJ_ORIGIN is +90 or -90?
        //  Or in this case, is "std_parallel" in CF-1 really the same as
        //   LAT_OF_ORIGIN in WKT, and we fill in LAT_OF_ORIGIN with +90 or 
        //   -90? (LAT_PROJ_ORIGIN set to 90 if STD_PARALLEL > 0, else -90)
        const oNetcdfSRS_PP mappings[] = {
            {LAT_PROJ_ORIGIN, SRS_PP_LATITUDE_OF_ORIGIN},
            {VERT_LONG_FROM_POLE, SRS_PP_CENTRAL_MERIDIAN},
            {SCALE_FACTOR_ORIGIN, SRS_PP_SCALE_FACTOR},  
            // STD_PARALLEL_1 ?
            {FALSE_EASTING, SRS_PP_FALSE_EASTING },  
            {FALSE_NORTHING, SRS_PP_FALSE_NORTHING },
            {NULL, NULL}
            };
            
        NCDFWriteProjAttribsFromMappings(poPROJCS, mappings, fpImage,
            NCDFVarID);
    }      
    // Rotated Pole: not implemented yet, unsure GDAL supports
    // Currently map Oblique steregraphic to stereographic
    else if( EQUAL(pszProjection, SRS_PT_OBLIQUE_STEREOGRAPHIC) ) {
        NCDFWriteProjAttribsFromMappings(poPROJCS, poStMappings, fpImage,
            NCDFVarID);
    }
    else if( EQUAL(pszProjection, SRS_PT_STEREOGRAPHIC) ) {
        NCDFWriteProjAttribsFromMappings(poPROJCS, poStMappings, fpImage,
            NCDFVarID);
    }
    else if( EQUAL(pszProjection, SRS_PT_TRANSVERSE_MERCATOR ) ) {
        NCDFWriteProjAttribsFromMappings(poPROJCS, poTMMappings, fpImage,
            NCDFVarID);
    }
    else {
        //PDS: PREVIOUS DEFAULT BEHAVIOUR
        for( int iChild = 0; iChild < poPROJCS->GetChildCount(); iChild++ )
        {
            const OGR_SRSNode *poNode;
            
            poNode = poPROJCS->GetChild( iChild );
            if( !EQUAL(poNode->GetValue(),"PARAMETER") 
                || poNode->GetChildCount() != 2 )
                continue;

            pszParamStr = poNode->GetChild(0)->GetValue();
            pszParamVal = poNode->GetChild(1)->GetValue();
        
            // pszNetCDFSRS = NULL;
            /* Find a better way to store and search these values */
            for(int i=0; poNetcdfSRS[i].netCDFSRS != NULL; i++ ) {
                if( EQUAL( poNetcdfSRS[i].SRS, pszParamStr ) ) {
                    CPLDebug( "GDAL_netCDF", "%s = %s", 
                              poNetcdfSRS[i].netCDFSRS, 
                              pszParamVal );
                    // pszNetCDFSRS = poNetcdfSRS[i].netCDFSRS;
                    // sscanf( pszParamVal, "%f", &dfValue );
                    // nc_put_att_float( fpImage, 
                    //                   NCDFVarID, 
                    //                   poNetcdfSRS[i].netCDFSRS,
                    //                   NC_FLOAT,
                    //                   1,
                    //                   &fValue );

                    /* Read STD_PARALLEL attribute */
                    if( EQUAL( poNetcdfSRS[i].SRS, SRS_PP_STANDARD_PARALLEL_1 ) )  {
                        bFoundStdP1 = TRUE;
                        sscanf( pszParamVal, "%lg", &dfStdP[0] );
                    }
                    else if( EQUAL( poNetcdfSRS[i].SRS, SRS_PP_STANDARD_PARALLEL_2 ) )   {
                        bFoundStdP2 = TRUE;
                        sscanf( pszParamVal, "%lg", &dfStdP[1] );
                    } 
                    /* Write attribute */ 
                    else {
                        dfTemp = atof( pszParamVal );
                        nc_put_att_double( fpImage, 
                                           NCDFVarID, 
                                           poNetcdfSRS[i].netCDFSRS,
                                           NC_DOUBLE,
                                           1,
                                           &dfTemp );
                    }
                    break;
                }
            }
        }
    }
    /* Write STD_PARALLEL attribute: special case */
    if ( bFoundStdP1 )  { 
        /* one value or equal values */
        if ( !bFoundStdP2 || dfStdP[0] ==  dfStdP[1] ) {
            nc_put_att_double( fpImage, 
                               NCDFVarID, 
                               STD_PARALLEL, 
                               NC_DOUBLE,
                               1,
                               &dfStdP[0] );
        }
        else { /* two values */
            nc_put_att_double( fpImage, 
                               NCDFVarID, 
                               STD_PARALLEL, 
                               NC_DOUBLE,
                               2,
                               dfStdP );
        }
    }
}

void NCDFWriteProjAttribsFromMappings(const OGR_SRSNode *poPROJCS,
                                        const oNetcdfSRS_PP mappings[],
                                        const int fpImage, const int NCDFVarID) 
{                            
    double dfStdP[2];
    int bFoundStdP1=FALSE,bFoundStdP2=FALSE;
    double dfTemp=0.0;
    const char *pszParamVal;

    for (int iMap = 0; mappings[iMap].GDAL_ATT != NULL; iMap++ ) {
        pszParamVal = GetProjParamVal(poPROJCS, mappings[iMap].GDAL_ATT);
        //Write param to NetCDF
        // include handling of std_parallel special case: this gets 
        // written at very end
        if( EQUAL(mappings[iMap].GDAL_ATT, SRS_PP_STANDARD_PARALLEL_1 ) ) {
            bFoundStdP1 = TRUE;
            sscanf( pszParamVal, "%lg", &dfStdP[0] );
        }
        else if( EQUAL(mappings[iMap].GDAL_ATT, SRS_PP_STANDARD_PARALLEL_2 ) ) { 
            bFoundStdP2 = TRUE;
            sscanf( pszParamVal, "%lg", &dfStdP[1] );
        } 
        else {
            dfTemp = atof( pszParamVal );
            nc_put_att_double( fpImage, 
                               NCDFVarID, 
                               mappings[iMap].NCDF_ATT,
                               NC_FLOAT,
                               1,
                               &dfTemp );
        }
    }
    /* Now handle the STD_PARALLEL attrib special case */
    if ( bFoundStdP1 ) { 
        /* one value or equal values */
        if ( !bFoundStdP2 || dfStdP[0] ==  dfStdP[1] ) {
            nc_put_att_double( fpImage, 
                               NCDFVarID, 
                               STD_PARALLEL, 
                               NC_DOUBLE,
                               1,
                               &dfStdP[0] );
        }
        else { /* two values */
            nc_put_att_double( fpImage, 
                               NCDFVarID, 
                               STD_PARALLEL, 
                               NC_DOUBLE,
                               2,
                               dfStdP );
        }
    }
}

// PDS
/* Helper function to get the value of a parameter of a projection in
 * a GDAL poPROJCS structure */
const char* GetProjParamVal(const OGR_SRSNode *poPROJCS, 
                            const char* findParamStr)
{
    const char *pszParamStr, *pszParamVal;

    for( int iChild = 0; iChild < poPROJCS->GetChildCount(); iChild++ )
            for( int iChild = 0; iChild < poPROJCS->GetChildCount(); iChild++ )
    {
        const OGR_SRSNode *poNode;
        
        poNode = poPROJCS->GetChild( iChild );
        if( !EQUAL(poNode->GetValue(),"PARAMETER") 
            || poNode->GetChildCount() != 2 )
            continue;

        pszParamStr = poNode->GetChild(0)->GetValue();
        pszParamVal = poNode->GetChild(1)->GetValue();
        if( EQUAL( pszParamStr, findParamStr ) ) 
            return pszParamVal;
    }
    return "";    
}

/************************************************************************/
/*                          GDALRegister_netCDF()                       */
/************************************************************************/

void GDALRegister_netCDF()

{
    GDALDriver	*poDriver;

    if (! GDAL_CHECK_VERSION("netCDF driver"))
        return;

    if( GDALGetDriverByName( "netCDF" ) == NULL )
    {
        poDriver = new GDALDriver( );
        
        poDriver->SetDescription( "netCDF" );
        poDriver->SetMetadataItem( GDAL_DMD_LONGNAME, 
                                   "Network Common Data Format" );
        poDriver->SetMetadataItem( GDAL_DMD_HELPTOPIC, 
                                   "frmt_netcdf.html" );
        poDriver->SetMetadataItem( GDAL_DMD_EXTENSION, "nc" );

        poDriver->pfnOpen = netCDFDataset::Open;
        // poDriver->pfnCreateCopy = NCDFCreateCopy;
        poDriver->pfnCreateCopy = NCDFCreateCopy2;
        poDriver->pfnIdentify = netCDFDataset::Identify;

        GetGDALDriverManager( )->RegisterDriver( poDriver );
    }
}

/* code taken from cdo and libcdi, used for writing the history attribute */
//void cdoDefHistory(int fileID, char *histstring)
void NCDFAddHistory(int fpImage, const char *pszAddHist, const char *pszOldHist)
{
    char strtime[32];
    time_t tp;
    struct tm *ltime;

    char *pszNewHist = NULL;
    size_t nNewHistSize = 0;
    int disableHistory = FALSE;

    /* Check pszOldHist - as if there was no previous history, it will be
       a null pointer - if so set as empty. */
    if (NULL == pszOldHist) {
        pszOldHist = "";
    }

    tp = time(NULL);
    if ( tp != -1 )
    {
        ltime = localtime(&tp);
        (void) strftime(strtime, sizeof(strtime), "%a %b %d %H:%M:%S %Y: ", ltime);
    }

    // status = nc_get_att_text( fpImage, NC_GLOBAL, 
    //                           "history", pszOldHist );
    // printf("status: %d pszOldHist: [%s]\n",status,pszOldHist);
    
    nNewHistSize = strlen(pszOldHist)+strlen(strtime)+strlen(pszAddHist)+2;
    pszNewHist = (char *) CPLMalloc(nNewHistSize * sizeof(char));
    
    strcpy(pszNewHist, strtime);
    strcat(pszNewHist, pszAddHist);

    if ( disableHistory == FALSE && pszNewHist )
    {
        strcat(pszNewHist, "\n");
        strcat(pszNewHist, pszOldHist);
    }

    nc_put_att_text( fpImage, NC_GLOBAL, 
                     "history", nNewHistSize,
                     pszNewHist ); 
  
    CPLFree(pszNewHist);
}


/* -------------------------------------------------------------------- */
/*      Set Lambert Conformal Conic Projection                          */
/* -------------------------------------------------------------------- */


    
//Albers equal area
//
//grid_mapping_name = albers_conical_equal_area
//
//Map parameters:
//
//    * standard_parallel - There may be 1 or 2 values.
//    * longitude_of_central_meridian
//    * latitude_of_projection_origin
//    * false_easting
//    * false_northing
//Azimuthal equidistant
//
//grid_mapping_name = azimuthal_equidistant
//
//Map parameters:
//
//    * longitude_of_projection_origin
//    * latitude_of_projection_origin
//    * false_easting
//    * false_northing
//Lambert azimuthal equal area
//
//grid_mapping_name = lambert_azimuthal_equal_area
//
//Map parameters:
//
//    * longitude_of_projection_origin
//    * latitude_of_projection_origin
//    * false_easting
//    * false_northing
//Lambert conformal
//
//grid_mapping_name = lambert_conformal_conic
//
//Map parameters:
//
//    * standard_parallel - There may be 1 or 2 values.
//    * longitude_of_central_meridian
//    * latitude_of_projection_origin
//    * false_easting
//    * false_northing
//Lambert Cylindrical Equal Area
//
//grid_mapping_name = lambert_cylindrical_equal_area
//
//Map parameters:
//
//    * longitude_of_central_meridian
//    * either standard_parallel or scale_factor_at_projection_origin
//    * false_easting
//    * false_northing
//Latitude-Longitude
//
//grid_mapping_name = latitude_longitude
//
//Map parameters:
//
//    * None
//Mercator
//
//grid_mapping_name = mercator
//
//Map parameters:
//
//    * longitude_of_projection_origin
//    * either standard_parallel or scale_factor_at_projection_origin
//    * false_easting
//    * false_northing
//Orthographic//
//grid_mapping_name = orthographic
//
//Map parameters:
//
//    * longitude_of_projection_origin
//    * latitude_of_projection_origin
//    * false_easting
//    * false_northing
//Polar stereographic
//
//grid_mapping_name = polar_stereographic
//
//Map parameters:
//
//    * straight_vertical_longitude_from_pole
//    * latitude_of_projection_origin - Either +90. or -90.
//    * Either standard_parallel or scale_factor_at_projection_origin
//    * false_easting
//    * false_northing
//Rotated pole
//
//grid_mapping_name = rotated_latitude_longitude
//
//Map parameters:
//
//    * grid_north_pole_latitude
//    * grid_north_pole_longitude
//    * north_pole_grid_longitude - This parameter is optional (default is 0.).
//Stereographic
//
//grid_mapping_name = stereographic
//
//Map parameters:
//
//    * longitude_of_projection_origin
//    * latitude_of_projection_origin
//    * scale_factor_at_projection_origin
//    * false_easting
//    * false_northing
//Transverse Mercator
//
//grid_mapping_name = transverse_mercator
//
//Map parameters:
//
//    * scale_factor_at_central_meridian
//    * longitude_of_central_meridian
//    * latitude_of_projection_origin
//    * false_easting
//    * false_northing
//Vertical perspective
//
//grid_mapping_name = vertical_perspective
//
//Map parameters:
//
//    * latitude_of_projection_origin
//    * longitude_of_projection_origin
//    * perspective_point_height
//    * false_easting
//    * false_northing
//
//
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


