////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       LogFilter pipeline I/O routine
//       GRID event report
//       Release: V0.0 -  16/mar/2005
//       Contributors: A.A., A.G., S.V., S.M.
//       Author: Alessio Trois (IASF-Milano)
//				 Alberto Pellizzoni
//
// INPUT
//       TBD
//
// OUTPUT
//       TBD
//
//
// FILE HISTORY
//       25/Apr/2005
//                      First release: V1.0
//                      Authors:       Alessio Trois (IASF-Milano)
//                             		   Alberto Pellizzoni(IASF-Milano)
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

#include <alikeLib.h>

int countsmalibur(AG_params & params) {

	int numcol = 0;
	long nrows = 0; 
	int status = 0;
	double l = 0, b = 0;	
	double x = 0, y = 0;	
	int i = 0, ii = 0;	
	double the = 0;	
	long mxdim= params.mxdim; // dimension (in pixels) of the map
	unsigned short A[mxdim][mxdim];

	for (i = 0; i < mxdim; i++)
	{   for (ii = 0; ii < mxdim; ii++)
		{
			A[i][ii] = 0;
		}
	}	

	double baa = params.ba * D2R;
	double laa = params.la * D2R;
	
	int bitpix   =  USHORT_IMG; /* 16-bit unsigned short pixel values       */
	long naxis    =   2;  /* 2-dimensional image                            */    
	long naxes[2] = { mxdim, mxdim };   /* image is 300 pixels wide by 200 rows */		
	
	fitsfile * evtFits;
	char tempname[FLEN_FILENAME];
	strcpy(tempname, tmpnam(NULL));
	if ( fits_create_file(&evtFits, tempname, &status) != 0 ) {
		printf("Errore in apertura file %s\n",tempname);
		return status;
	}	
	
	char expr[1024];
	strcpy(expr,params.evtexpr().c_str());
	
	std::cout << std::endl << "AG_ctsmapgen....................................adding events files"<< std::endl;
	status = addfile(evtFits, params.evtfile, expr, params.tmin, params.tmax);
	std::cout << "AG_ctsmapgen....................................addfile exiting STATUS : "<< status<< std::endl << std::endl ;	
	

	fitsfile * mapFits;
	if ( fits_create_file(&mapFits, params.outfile, &status) != 0 ) {
		printf("Errore in apertura file '%s'\n",params.outfile);
		return status;
	}	

	
	fits_movabs_hdu(evtFits, 2, NULL, &status);	
	fits_get_num_rows(evtFits, &nrows, &status);
	cout << nrows << endl;

	double ra, dec;
	double dummy;
	switch (params.projection) {
	    case AG_params::ARC:
		for (long k = 0; k<nrows; ++k) {
			fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
			fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);

		    eulerold(ra, dec, &l, &b, 1);
		    l*=D2R;
		    b*=D2R;
		    the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
		    if (the < -1.0) {
			    the = PI;
		    } else if (the > 1.0) {
			    the = 0.0;
		    } else {
			    the = acos(the);
		    }
		    x = R2D/Alikesinaa(the) * cos(b)*sin(l-laa);
		    y = R2D/Alikesinaa(the) * (sin(b)*cos(baa) - cos(b)*sin(baa)*cos(l-laa));

		    i=(int)floor(((-x+(params.mdim/2.))/params.mres));
		    ii=(int)floor(((y+(params.mdim/2.))/params.mres));

		    if (params.inmap(i,ii)) {
			    A[ii][i]+=1;
			}
		}
		break;
		
	    case AG_params::AIT:
		for (long k = 0; k<nrows; ++k) {
			fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
			fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);
		    eulerold(ra, dec, &l, &b, 1);
		    l*=D2R;
		    b*=D2R;
		    the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
		    if (the < -1.0) {
			    the = PI;
		    } else if (the > 1.0) {
			    the = 0.0;
		    } else {
			    the = acos(the);
		    }
		    l=l-laa;

		    if ( l < PI  ) { 
		      l=-l; 
		    }
		    else { 
		      l=2*PI -l; 
		    }

		    x=R2D*(sqrt(2.0)*2.0*cos(b)*sin(l/2.0))/sqrt(1.0 + cos(b)*cos(l/2.0) ) ;         
		    y=R2D*(sqrt(2.0)*sin(b))/sqrt(1.0 + cos(b)*cos(l/2.0) );

		    i=(int)floor(((x+(params.mdim/2.))/params.mres));
		    ii=(int)floor(((y+(params.mdim/2.))/params.mres));

		    if (params.inmap(i,ii)) {
			    A[ii][i]+=1;
			}
		}
		break;
	}
	
		

	long nelement =  naxes[0] * naxes[1];
	std::cout<< "creating Counts Map...................................." << std::endl;	
	fits_create_img(mapFits, bitpix, naxis, naxes, &status);
	std::cout<< "writinig Counts Map...................................." << std::endl;		
	fits_write_img(mapFits, bitpix, 1, nelement, A, &status);	
	std::cout<< "writing header........................................" << std::endl<< std::endl;	

		
	
	params.write_fits_header(mapFits, status);
	
	
	fits_delete_file(evtFits, &status);
	fits_close_file(mapFits, &status);	
	return status;
}




int main(int argc,char **argv)
{

int status = 0;
int numpar=0;


std::cout << " " << std::endl;
std::cout << " " << std::endl;
std::cout << "#################################################################" << std::endl;
std::cout << "########## AG_ctsmapgen.cpp v.0 - 19/12/05 - A.C., A.T. #########" << std::endl;
std::cout << "#################################################################" << std::endl;
std::cout << "#################################################################" << std::endl;
std::cout << " " << std::endl;

AG_params params(AG_params::CTS, argc, argv, numpar, status);
std::cout << "AG_ctsmapgen...............................starting"<< std::endl;	
status = countsmalibur(params);
std::cout << "AG_ctsmapgen............................... exiting"<< std::endl;

		
if (status) {
	if (status != 105) {
		remove(params.outfile);
		}	
	printf("AG_ctsmapgen..................... exiting AG_ctsmapgen ERROR:");		
	fits_report_error(stdout, status);	
	return status;			
	}			
else {
	printf("\n\n\n###################################################################\n");
	printf("#########  Task AG_ctsmapgen........... exiting #################\n");
	printf("#################################################################\n\n\n");					
	}			

return status;

return status;
}

