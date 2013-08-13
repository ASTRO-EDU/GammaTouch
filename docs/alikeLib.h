///////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       alikeLib.h
//       Release: V1.0 -  9/Oct/2010
//       Contributors: 
//       Author: Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano)
//
// INPUT
//       TBD
//
// OUTPUT
//       TBD
//
//
// FILE HISTORY
//       9/Oct/2010
//                      First release: V1.0
//       		Author: Andrew Chen, Alessio Trois (IASF-Milano)
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstring>
#include <stdio.h>
//#include <cstdlib>
#include <fcntl.h>
#include <math.h>
#include <time.h>
#include <TMath.h> //AB1
#include "TH2.h"
#include "TMatrixD.h"
#include "TMatrixF.h"
#include "TVectorD.h"
#include "TVectorF.h"
#include "TMinuit.h"
//#include "TMath.h"
#include <unistd.h>
#include <ctype.h>
//#include <iterator>
#include <cstring>
#include <pil.h>
#include "GridUtilities.h"
#include "alikeQuat.h"
#include "fitsio.h"
#include "wcshdr.h"
#include "wcsmath.h"
#include "wcstrig.h"
#include "sph.h"

void exit_error(const char * error_string);

int excalibur(char * fileList, char *  outfile, char * raeffFileName, double mdim, double mres, double la, double ba, double tmin, double tmax, double ene, double lonpole);


Double_t Alikesinaa(Double_t input);
int alikePolrec(double R, double A, double * X, double * Y, int degrees);
int alikeRecpol(double x, double y, double *r, double *a, int degrees);
double alikeSphdist(double long1, double lat1, double long2, double lat2, int degrees);
double AlikeSphdistDeg(double long1, double lat1, double long2, double lat2);

int eulerold(double ai, double bi, double * ao, double * bo, int select);	
int euler(double ra, double dec, double cc, double * ao, double * bo, double *co);

int alikeReadAeff(const char * raeffFileName, float ****raeffgrid, float **raefftheta, float **raeffphi, float **raeffenergy, long * naxes);	
int alikeReadPsf(const char * psfFileName, float ******psfgrid, float **psfrho, float **psfpsi, float **psftheta, float **psfphi, float **psfenergy, long * naxes);
double alikeAeff(double ene, double theta,  double phi, float ***raeffgrid, float *raefftheta, float *raeffenergy, long * naxes);	
double alikePsf(float ******psfgrid, float **psfrho, float **psfpsi, float **psftheta, float **psfphi, float **psfenergy, double rho, double psi, double theta, double phi, double ene, long * naxes);

int alikeDeletePsf(float ******psfgrid, float **psfrho, float **psfpsi, float **psftheta, float **psfphi, float **psfenergy, long * naxes);
int alikeDeleteAeff(float ****raeffgrid, float **raefftheta, float **raeffphi, float **raeffenergy, long * naxes);

int alikeReadEdp(const char * edpFileName, float *****edpgrid, float **edptrueenergy, float **edpobsenergy, float **edptheta, float **edpphi, long * naxes);
double alikeEpd(float *****edpgrid, float **edptrueenergy, float **edpobsenergy, float **edptheta, float **edpphi, double trueenergy, double obsenergy, double theta, double phi, long * naxes);
int alikeDeleteEdp(float *****edpgrid, float **edptrueenergy, float **edpobsenergy, float **edptheta, float **edpphi, long * naxes);

class AG_params{
public:
	
	char evtfile[FLEN_FILENAME];
	char logfile[FLEN_FILENAME];
	char outfile[FLEN_FILENAME];	
	char raeffFileName[FLEN_FILENAME];
	char edpFileName[FLEN_FILENAME]; //AB1
	enum ProgramType {CTS, EXP, GAS};
	ProgramType program;
	double tmin, tmax, emin, emax, lonpole;
	double mdim, mres, la, ba;
	long mxdim;
	double albrad, fovrad, fovradmin;
	double spectral_index, y_tol, roll_tol, earth_tol;
	int keepmono;
	int phasecode; // 0 = entire orbit, 1 = phase = 0 only
	int filtercode; // 0 = L & G, 1 = G only
	int stepsize;
	int timestep;	
	enum ProjType {ARC, AIT};
	ProjType projection;

	AG_params(ProgramType programin, int &argc, char **&argv, int & numpar, int & status);
	AG_params(const char filename[]);
	AG_params(fitsfile * infile);
	AG_params(const AG_params & in): program(in.program),tmin(in.tmin),tmax(in.tmax),emin(in.emin),emax(in.emax),lonpole(in.lonpole),mdim(in.mdim),mres(in.mres),la(in.la),ba(in.ba),mxdim(in.mxdim),albrad(in.albrad),fovrad(in.fovrad),fovradmin(in.fovradmin),spectral_index(in.spectral_index),y_tol(in.y_tol),roll_tol(in.roll_tol),earth_tol(in.earth_tol),keepmono(in.keepmono),phasecode(in.phasecode),filtercode(in.filtercode),stepsize(in.stepsize),timestep(in.timestep),projection(in.projection){
		strcpy(evtfile,in.evtfile);
		strcpy(logfile,in.logfile);
		strcpy(outfile,in.outfile);
		strcpy(raeffFileName,in.raeffFileName);
	}
//	AG_params(fitsfile *, int * status);
	~AG_params(){}
	
	int ReadFromFITS(fitsfile * infile);
	string evtexpr();
	string logexpr();
	int write_fits_header(fitsfile *mapFits, int & status);
	inline bool inmap(int i, int ii){return (i < mxdim) && (i >= 0) && (ii < mxdim) && (ii >= 0);}
	inline bool etest(double e){return (e >= emin) && (e <= emax);}
	inline bool ttest(double t){return (t >= tmin) && (t <= tmax);}
	inline bool albtest(double ph_earth){return (ph_earth > albrad);}
	inline bool fovtest(double theta){return (theta < fovrad && theta >= fovradmin);}
	inline bool albtest(double ra, double dec, double earth_ra, double earth_dec){return albtest(alikeSphdist(ra, dec, earth_ra, earth_dec, 1));}
	inline bool phasetest(short phase){return (phasecode ? 1 : (phase == 0));}
	bool y_toltest(double & ra_y, double & dec_y, double & ra_y0, double & dec_y0){return alikeSphdist(ra_y, dec_y, ra_y0, dec_y0, 1) > y_tol;}
	bool earth_toltest(double & earth_ra, double & earth_dec, double & earth_ra0, double & earth_dec0){return alikeSphdist(earth_ra, earth_dec, earth_ra0, earth_dec0, 1) > earth_tol;}
	bool roll_toltest(double * psi){return fabs(psi[0] - psi[1]) > roll_tol;}
};

class AlikeAeffGridClass
{

 protected:
	float ****raeffgrid;
	float **raefftheta;
	float **raeffphi;
	float **raeffenergy;
	long naxes[3];

 public:

 	AlikeAeffGridClass(){};
	AlikeAeffGridClass(const char * aeffFileName);
	~AlikeAeffGridClass();

	double Val(double ene, double theta, double phi) ;
	TVectorF energies() ;

};


class AlikePsfGridClass
{

 private:
	float ******psfgrid;
	float **psfrho;
	float **psfpsi;
	float **psftheta;
	float **psfphi;
	float **psfenergy;
	long naxes[5];

 public:

 	AlikePsfGridClass(){};
	AlikePsfGridClass(const char * psfFileName);
	~AlikePsfGridClass();
	
	int getAxes(int n) ;
	double Val(double rho, double psi, double theta, double phi, double energy) ;
	TVectorF energies() ;
	TVectorF rhos() ;
	TVectorF psis() ;
};







class AlikeEdpGridClass
{

 private:
 
	float *****edpgrid;
	float **edptrueenergy;
	float **edpobsenergy;
	float **edptheta;
	float **edpphi;
	long naxes[4];

 public:

 	AlikeEdpGridClass(){}
	AlikeEdpGridClass(const char * edpFileName);
	~AlikeEdpGridClass();
	
	int getAxes(int n) ;

	double Val(double E_true, double E_obs, double theta, double phi) ;
	TVectorF trueenergies() ;
	TVectorF obsenergies() ;
};

