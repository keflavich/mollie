/*-----------------------------------------------------------------------*/
/* REQUIRED INCLUDES and GLOBAL VARIABLES. DON'T CHANGE THESE */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#define MAXADJP 16      /* maximum number of adjustable parameters */

static float maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define MIN(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) < (maxarg2) ?\
        (maxarg1) : (maxarg2))

extern FILE *fpnotes;
extern FILE *fpspectra;
extern FILE *fphist;
extern float crtsph_(float *,float *,float *, float *,float *,float *, int *);
extern int fft2d (int inverse, int nx, int ny, float **beamR, float **beamI,
           float **fftR, float **fftI);
extern int locate_(float *, int *, float *, int *);

float chisq_function( );
int read_data();
int set_model_parameters(int , float *, float **, int, int );
int ulrich_flow(int i, int j, int k, float rnd, float el,
       float* vr, float* vaz, float* vel, float* density);
float flared_disk_density(int i, int j, int k, float x, float z,
                float rstar, float rd);
int smooth ( int nx, int ny, int nz,
                float beamx, float beamy,
                float cellx, float celly,
                float * xcenter, float * ycenter,
                float *** input_cube, float *** smoothed_cube );

/* --------------------------------------------------------------------------- */

struct dcube {
        int nlines;
        int nviews;
        int *nchan;
        int noutx;
        int nouty;
        float celly;
        float cellx;
        float beamy;
        float beamx;
        char **linenames;
	double *restfreq;
        float **chvel;
        float *lng;
        float *lat;
        float *yplane;
        float *xplane;
        float *ycenter;
        float *xcenter;
        float *****cube;
        float ***wt;
};

extern struct dcube dcube;

struct data {
        float vel[65][4];
        float data[65][3][4];
};

struct data data;


struct model {
        float temperature;
        float density;
        float vx,vy,vz;
        float abundance;
        float linewidth;
        int adds;
        int phase;
        float mean_av;
        float dust_temperature;
};
extern struct model ****model;


/*-----------------------------------------------------------------------*/


/*-----------------------------------------------------------------------*/
/* OPTIONAL GLOBAL VARIABLES */

#define PLANCK  6.6252e-27
#define BOLTZ   1.3806e-16
#define C       3.e10
#define AMU     1.660531e-24
#define PC      3.085e18
#define AU      1.5e13
#define GRAV    6.67e-8
#define PI      3.14159265359
#define SQRTPI  1.77245385091
#define MSUN    2.e33
#define YR      3.15360e7
#define RGAS    8.31421e7
#define STEFANB 5.67051E-5
#define RSUN    7.e10
#define LSUN    3.839e33
#define MH      1.67e-24

struct parameters {
        int model_number;
        float env_density;
        float disk_density_factor;
        float point_mass;
	float abundance;
        float gamma0;
        float turbwidth;
        float disk_temperature_factor;
        float env_temperature_factor;
        float planck_exp;
        char  filename[64];
};
struct parameters parameters;


struct hydro{

        int jx0;
        int nsets;
        float shock_smoothing;
        float *sn;
        float *mdn0;
        float *pnt0;
        float *et;
        float *dt;
        float *pext;
        float **rn;
        float **mdn;
        float **tn;
        float **pn;
        float **pntn;
        float **vn;
        float **gn;
        float **qn;
        float **dpn;
        float **dqn;
        float **tdust;
        float **co_abundance;
        float **avgavn;
        float **gpe;
        float **h2o_abundance;

};
struct hydro hydro;

struct lynds{

        float lambda;
        float scale_alpha;
        float scale_beta;
        float bnd_temp;
        double total_mass;
        float pext;
        float boundary;
        float xn;
        int   nsteps;
        int   nt;
        float *x;
        float *s;
        float *z;
        float *th;
        float *tgas;
        float *tdust;

};
struct lynds lynds;



struct avery{

        float x;
        float y;
        float z;
        int   i;
        int   j;
        int   k;
        float density;
        float temperature;
        float ionfraction;
        float neutraldensity;
        float moleculardensity;
        float vx;
        float vy;
        float vz;

};
struct avery ***avery;

int read_avery(char filename[], int si_units);
int read_ulrich(char filename[]);
int read_flash(char filename[]);
int read_lynds(char filename[], int ntp);

/*-----------------------------------------------------------------------*/
/* THE FUNCTION define_model */

/* DON'T CHANGE ARGUMENTS TO THIS FUNCTION */
int define_model(int model_number, int nbox, int nx[],int ny[], int nz[],
                 float **xcenter,  float **ycenter,  float **zcenter,
                 float atoms, float lwmin, float lwmax, int depletion, 
                 int photodissociation, float std_abundance, int molecule) 
{


/* Variables such as these are usually used */

int i,j,k,ibox,m;			/* Loop indices for the model grid */

int idir;			/* parameter for function crtsph that
                                   specifies to convert from cartesian
                                   to spherical (idir=1) or the inverse
                                   (idir = -1) */

float x,y,z;			/* cartesian and spherical coords can */
float r,az,el;			/* be used toegether with the conversion */
				/* function crtsph */
int si_units=1.;		/* For the Avery files, set this to 1 
				   if the input is in SI units. Also converts
				   density from grams to number. 
				   Otherwise 0 for cgs and number density */

/*-----------------------------------------------------------------------*/
/* Any other local variables can be defined here. */

float turbwidth, temperature,r0,t0,t1,d0,abundance;
char  filename[128];
int   jlo,jhi,ip;
float *hsr, *hst, *htd, *hsn, *hsv, *hsa, *hsh;
float vr,vmin,vmax,alpha,rmax;

float sinaz,cosaz;
float sinel,cosel;

FILE  *fpout;
int   ipout;


/*-----------------------------------------------------------------------*/

/* Here is the file to read in */
strcpy(filename,"hydro_std_1E7\0");
printf("Data file is %s\n",filename);
fprintf(fpnotes,"Data file is %s\n",filename);



/* read in the hydro data */
        i = 0;          /* this is a flag that indicates a non-thermal pressure file */
        j = read_lynds(filename,i);
        if (j == 0) {
                printf("read_lynds returned failed\n");
                fprintf(fpnotes,"read_lynds returned failed\n");
                return 0;
        }

/* JLO is the set number corresponding to the evolutionary stage that you want to model. */

        jlo = 1;

/* Copy the hydro data into 1D local arrays. Convert R from CM to PC */
        hsr = (float *)malloc((hydro.jx0+1)*sizeof(float));
        hst = (float *)malloc((hydro.jx0+1)*sizeof(float));
        htd = (float *)malloc((hydro.jx0+1)*sizeof(float));
        hsn = (float *)malloc((hydro.jx0+1)*sizeof(float));
        hsv = (float *)malloc((hydro.jx0+1)*sizeof(float));
        hsa = (float *)malloc((hydro.jx0+1)*sizeof(float));
        hsh = (float *)malloc((hydro.jx0+1)*sizeof(float));
        vmin = 1.e20;
        vmax = -1.e20;
        for (i=0;i<hydro.jx0+1;i++) {
                hsr[i] = hydro.rn[i][jlo]/PC;
                hst[i] = hydro.tn[i][jlo];
                htd[i] = hydro.tdust[i][jlo];
                hsn[i] = hydro.mdn[i][jlo]/(2.33*AMU);
                hsv[i] = hydro.vn[i][jlo];
                hsa[i] = hydro.co_abundance[i][jlo];
                hsh[i] = hydro.h2o_abundance[i][jlo];
                if (hsv[i] < vmin) vmin = hsv[i];
                if (hsv[i] > vmax) vmax = hsv[i];
/* fprintf(fpnotes," %10.3e  %10.3e %10.3e %10.3e\n",hsr[i],hst[i],hydro.tdust[i][jlo],hsn[i]); */
        }

        printf("radius   at center %e and boundary %e\n",hsr[0]*PC,hsr[hydro.jx0]*PC);
        printf("radius   at center %e and boundary %e\n",hsr[0],hsr[hydro.jx0]);
        printf("temp     at center %f and boundary %f\n",hst[0],hst[hydro.jx0]);
        printf("dust T   at center %f and boundary %f\n",htd[0],htd[hydro.jx0]);
        printf("density  at center %e and boundary %e\n",hsn[0],hsn[hydro.jx0]);
        printf("velocity at center %f and boundary %f\n",hsv[0],hsv[hydro.jx0]);
        printf("min max velocity %f %f in cm/s \n",vmin,vmax);
        printf("CO abundance at center %e and boundary %e\n",hsa[0],hsa[hydro.jx0]);
        printf("H2O abundance at center %e and boundary %e\n",hsh[0],hsh[hydro.jx0]);
        printf("now try to fill in the model grid \n");

        fprintf(fpnotes,"radius   at center %e and boundary %e\n",hsr[0]*PC,hsr[hydro.jx0]*PC);
        fprintf(fpnotes,"radius   at center %e and boundary %e\n",hsr[0],hsr[hydro.jx0]);
        fprintf(fpnotes,"temp     at center %f and boundary %f\n",hst[0],hst[hydro.jx0]);
        fprintf(fpnotes,"dust T   at center %f and boundary %f\n",htd[0],htd[hydro.jx0]);
        fprintf(fpnotes,"density  at center %e and boundary %e\n",hsn[0],hsn[hydro.jx0]);
        fprintf(fpnotes,"velocity at center %f and boundary %f\n",hsv[0],hsv[hydro.jx0]);
        fprintf(fpnotes,"min max velocity in cm/s %f %f\n",vmin,vmax);
        fprintf(fpnotes,"CO abundance at center %e and boundary %e\n",hsa[0],hsa[hydro.jx0]);
        fprintf(fpnotes,"H2O abundance at center %e and boundary %e\n",hsh[0],hsh[hydro.jx0]);
        fprintf(fpnotes,"now try to fill in the model grid \n");







/*-----------------------------------------------------------------------*/

/* Fill in the model grid with equations */
/* There are 4 loops below with indices, ibox,i,j,k. The model is defined here, either
   with simple equations or from a file created by another code.

  The minimum information required is
	Temperature
	Density. Molecular hydrogen number density.
	Abundance of the tracer molecule that is being modeled.
	3D velocity vx,vy,vz
	linewidth usually thermal plus some additional width for unresolved turbulence
		within a cell

  Other information that could be entered if the code is set up to use it.
	Phase
		phase = 0 is molecular gas and the default.
		phase = 2 is for ionized gas. A cell with phase = 2 will emit 
			free-free radiation at the gas temperature in the cell.
	Dust temperature
	transmission of starlight. The variable name is mean_av but the quantity required for
		chemistry module is exp(-mean_av) which is the transmission.

These variables are all in the struct model. The one additional variable "adds" is the number of
rays passing through each cell. This is computed by the code to correctly normalize Jbar. 
It is not used to define the model.

*/

        printf("Now fill in the model grid \n");


turbwidth = 0.08e5;

rmax = -xcenter[0][0];

for (ibox=0;ibox<nbox;ibox++) {


 for (i=0;i<nx[ibox];i++) {
 for (j=0;j<ny[ibox];j++) {
 for (k=0;k<nz[ibox];k++) {

    x = xcenter[ibox][i];
    y = ycenter[ibox][j];
    z = zcenter[ibox][k];
    idir = 1;
    crtsph_(&x,&y,&z,&r,&az,&el,&idir);
    sinaz = sin(az);
    sinel = sin(el);
    cosaz = cos(az);
    cosel = cos(el);

// fprintf(fpnotes,"xyz r %e %e %e %e \n",x,y,z,r);

/* ----------------------------------------------------------------------------------- */
/* This is some code for a partially ionized model */

//        model[ibox][i][j][k].phase = 0;
//       if (avery[i][j][k].temperature > 2000.) model[ibox][i][j][k].phase = 2;

/* This is for the ionized gas */
//        model[ibox][i][j][k].temperature = avery[i][j][k].temperature;
//        model[ibox][i][j][k].density = avery[i][j][k].density; /* *avery[i][j][k].ionfraction; */

/* This is for the molecular gas */
//	if (model[ibox][i][j][k].phase == 0) {
//            model[ibox][i][j][k].temperature = avery[i][j][k].temperature
//                * avery[i][j][k].moleculardensity / avery[i][j][k].density;
//            model[ibox][i][j][k].density = avery[i][j][k].moleculardensity;
//	}
/* ----------------------------------------------------------------------------------- */


/* Interpolate the temperature,density, and velocity */
/* The linewidth is set by the temperature */

/* locate the current r in the data array. Here we want to interpolate,
   so we need both JLO and JHI. */

        for (jhi=0;jhi<hydro.jx0;jhi+=10) {
                if (hsr[jhi] > r) break;
        }
        for (jlo=MAX(0,jhi-10);jlo<jhi;jlo++) {
                if (hsr[jlo] >= r) break;
        }
        jlo = jlo - 1;
        jlo = MAX(jlo,0);
        jlo = MIN(jlo,hydro.jx0-2);
        jhi = jlo+1;

        model[ibox][i][j][k].temperature = hst[jlo] +
                (r - hsr[jlo])*(hst[jhi]-hst[jlo])/(hsr[jhi]-hsr[jlo]);

// This keeps the temperature down below 100 K as we go out to larger radii.
// The hydro_rewrite4b models did not have colling by atomic oxygen.
        model[ibox][i][j][k].temperature = model[ibox][i][j][k].temperature  
		* (1. - 1.5*exp(-1./r));

        model[ibox][i][j][k].dust_temperature = htd[jlo] +
                (r - hsr[jlo])*(htd[jhi]-htd[jlo])/(hsr[jhi]-hsr[jlo]);

/* Here the linewidth is the microturbulent linewidth added in
   quadrature to the thermal linewidth */

        model[ibox][i][j][k].linewidth   =
                sqrt(turbwidth*turbwidth
                  + BOLTZ* model[ibox][i][j][k].temperature
                     /(atoms*AMU));

        model[ibox][i][j][k].density = hsn[jlo] +
                (r - hsr[jlo])*(hsn[jhi]-hsn[jlo])/(hsr[jhi]-hsr[jlo]);

// if (r < 0.0045) model[ibox][i][j][k].density = 1.e7;

        if (model[ibox][i][j][k].density < 10.)
                model[ibox][i][j][k].density = 10.;

// The abundance needs to be defined only for the molecule that is modeled.
// This switch is for convenience. I can change setup.c for different molecules
// and not have to change define_model.c. The Cnotes0 output does not keep
// track of which file define_model.c is used. The result, that is the model
// is summarized, but not the file that created it.
switch (molecule) {

	case 7:  /* This is CO abundance */
        model[ibox][i][j][k].abundance = hsa[jlo] +
                (r - hsr[jlo])*(hsa[jhi]-hsa[jlo])/(hsr[jhi]-hsr[jlo]);
	if (r > rmax) model[ibox][i][j][k].abundance = hsa[jlo];
	break;

	case 10:  /* This is 13CO abundance */
        model[ibox][i][j][k].abundance = hsa[jlo] +
                (r - hsr[jlo])*(hsa[jhi]-hsa[jlo])/(hsr[jhi]-hsr[jlo]);
	if (r > rmax) model[ibox][i][j][k].abundance = hsa[jlo];
        model[ibox][i][j][k].abundance /= 77.;
	break;

	case 13:  /* This is C17O abundance */
        model[ibox][i][j][k].abundance = hsa[jlo] +
                (r - hsr[jlo])*(hsa[jhi]-hsa[jlo])/(hsr[jhi]-hsr[jlo]);
	if (r > rmax) model[ibox][i][j][k].abundance = hsa[jlo];
        model[ibox][i][j][k].abundance /=  (77. * 7.3 * 4.0);
	break;

	case 14:  /* This is C18O abundance */
        model[ibox][i][j][k].abundance = hsa[jlo] +
                (r - hsr[jlo])*(hsa[jhi]-hsa[jlo])/(hsr[jhi]-hsr[jlo]);
	if (r > rmax) model[ibox][i][j][k].abundance = hsa[jlo];
        model[ibox][i][j][k].abundance /=  (77. * 7.3);
        model[ibox][i][j][k].abundance = std_abundance;
	break;

	case 6: /* This is H2O abundance */
        model[ibox][i][j][k].abundance = hsh[jlo] +
                (r - hsr[jlo])*(hsh[jhi]-hsh[jlo])/(hsr[jhi]-hsr[jlo]);
/*
		fprintf(fpnotes,"1 ibox ijk jlo jhi %d %d %d %d %d %d\n",ibox,i,j,k,jlo,jhi);
		fprintf(fpnotes,"1 hsrlo hsrhi %e %e %e\n",hsr[jlo],hsr[jhi],r);
		fprintf(fpnotes,"1 hshlo hshhi %e %e %e\n",hsh[jlo],hsh[jhi],model[ibox][i][j][k].abundance);
*/
	if (r > rmax) model[ibox][i][j][k].abundance = hsh[jlo];
//        model[ibox][i][j][k].abundance = 1.e-9;
	break;

        case 4:  /* This is N2H+ abundance */
        model[ibox][i][j][k].abundance = std_abundance;
        break;


	default:
		printf("Bad number for molecule in define_model.c. Molecule = %d\n",molecule);
                fprintf(fpnotes,"Bad number for molecule in define_model.c. Molecule = %d\n",molecule);
		return(0);
	break;
}

	if (model[ibox][i][j][k].abundance <= 0.0){
//	if (k == nz[ibox]/2 && j == ny[ibox]/2) {
		fprintf(fpnotes,"2 ibox ijk jlo jhi %d %d %d %d %d %d\n",ibox,i,j,k,jlo,jhi);
		fprintf(fpnotes,"2 hsrlo hsrhi %e %e %e\n",hsr[jlo],hsr[jhi],r);
		fprintf(fpnotes,"2 hshlo hshhi %e %e %e\n",hsh[jlo],hsh[jhi],model[ibox][i][j][k].abundance);
	}


/*
printf("ijk %d %d %d r hsrLo hsrHi %e %e %e h2o %e %e Abundance %e\n",
	i,j,k,r,hsr[jlo],hsr[jhi],hsh[jlo],hsh[jhi],model[ibox][i][j][k].abundance);
fprintf(fpnotes,"ijk %d %d %d r hsrLo hsrHi %e %e %e h2o %e %e Abundance %e\n",
	i,j,k,r,hsr[jlo],hsr[jhi],hsh[jlo],hsh[jhi],model[ibox][i][j][k].abundance);
*/

        vr = hsv[jlo] +
                (r - hsr[jlo])*(hsv[jhi]-hsv[jlo])/(hsr[jhi]-hsr[jlo]);;

        model[ibox][i][j][k].vx = vr*cos(az)*sin(el);
        model[ibox][i][j][k].vy = vr*sin(az)*sin(el);
        model[ibox][i][j][k].vz = vr*cos(el);

        }}}}

/*-----------------------------------------------------------------------*/
/* Some minimal checks */
	m = 1.;
        for (ibox=0;ibox<nbox;ibox++){
        printf("Model defined over box %d with nx,ny,nz = %d %d %d\n",
		ibox,nx[ibox],ny[ibox],nz[ibox]);
        fprintf(fpnotes,"Model defined over box %d with nx,ny,nz = %d %d %d\n",
		ibox,nx[ibox],ny[ibox],nz[ibox]);
        for (j=0;j<ny[ibox];j++){
        for (k=0;k<nz[ibox];k++){
        for (i=0;i<nx[ibox];i++){
		if (model[ibox][i][j][k].temperature <= 0. ||
		    model[ibox][i][j][k].density     <= 0. ||
		    (model[ibox][i][j][k].abundance  <= 0. && model[ibox][i][j][k].phase == 0) ||
		    model[ibox][i][j][k].linewidth   <= 0. 
		) {

if (m == 1) {
             printf("Part of the model grid is not defined, %d %d %d %d\n",ibox,i,j,k);
             printf("T, n, A, dv %e %e %e %e\n",model[ibox][i][j][k].temperature ,
                    model[ibox][i][j][k].density,
                    model[ibox][i][j][k].abundance,
                    model[ibox][i][j][k].linewidth);
	     m = 0.;
}

             fprintf(fpnotes,"Part of the model grid is not defined, %d %d %d %d\n",ibox,i,j,k);
             fprintf(fpnotes,"T, n, A, dv %e %e %e %e\n",model[ibox][i][j][k].temperature ,
                    model[ibox][i][j][k].density,
                    model[ibox][i][j][k].abundance,
                    model[ibox][i][j][k].linewidth);
                  }

        }}}}

/* This finishes the model inputs. */

/*	printf("Model defined\n");  */
	fflush(stdout);
	fprintf(fpnotes,"Model defined\n");

//  if (m == 0) return 0;
  return 1;
}

/* THIS IS THE FUNCTION TO READ HYDRO OUTPUT */
/* A similar function could be used for other hydro code outputs with
   different formats. */

int read_avery(char filename[], int si_units) 
{

/* Function mfsize determines the size of a file on disk. Sometimes this is
   useful. The function mfsize itself is in cmain.c.  */

extern unsigned long mfsize(FILE *);
 
FILE *fplynds;
int   iplynds;
unsigned long  file_size;

int i,j,k,l,m,n,ii,jj;
long nobj,nbytes,nbytes_remain;
int  nsets,nbytes_per_set,ngrid;
int my_rank;
int ntp = 0;
/* These conversion factors are changed below 
	if the input file is in SI units and mass_density */

float float_var;
double mindensity=1.e20,maxdensity=-1.e20,
	minx=1.e20,maxx=-1.e20,
	miny=1.e20,maxy=-1.e20,
	minz=1.e20,maxz=-1.e20,
	minvx=1.e20,maxvx=-1.e20,
	minvy=1.e20,maxvy=-1.e20,
	minvz=1.e20,maxvz=-1.e20;
long im=0,jm=0,km=0,nm=0;

/*
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
printf("P%d: starting to read data for cloud \n",my_rank);
*/

fplynds = fopen(filename,"r");
if (fplynds == NULL) {
	printf("Could not open data file %s for Avery clouds. Exiting ...\n",
		filename);
	fprintf(fpnotes,"Could not open data file %s for Avery clouds. Exiting ...\n",
		filename);
        return 0;
} else {
	printf("Opened data file %s for Avery clouds\n",
		filename);
	fprintf(fpnotes,"Opened data file %s for Avery clouds\n",
		filename);
}
iplynds = fileno(fplynds);

file_size = mfsize(fplynds);
printf("size of input file is %lu\n",file_size);

/* allocate memory for the model, the size of the model grid will be the first integer
   in the data file */


/*
	nobj = fscanf(fplynds,"%d",&ngrid);
	printf("the number of points in the input grid is %d\n",ngrid);
*/

	ngrid = 60;

    avery = (struct avery***) malloc(ngrid*sizeof(struct avery**));
    for (i=0;i<ngrid;i++){
      avery[i] = (struct avery **) malloc(ngrid*sizeof(struct avery *));
      for (j=0;j<ngrid;j++){
        avery[i][j] = (struct avery *) malloc(ngrid*sizeof(struct avery));
        }}

nbytes = 0;
nobj = 0;

	printf("memory allocated, start reading input file\n");

/* Loop here over ngrid cube */
/* This version reverses the Y axis */

	for (i=0;i<ngrid;i++){
	for (j=0;j<ngrid;j++){
	for (k=0;k<ngrid;k++){
		ii = ngrid - i -1;
		jj = ngrid - j -1;
		n = (long)fscanf(fplynds, "%e %e %e %e %e %e %e ",
			&avery[i][j][k].x,
			&avery[i][j][k].y,
			&avery[i][j][k].z,
			&avery[ii][jj][k].density,
			&avery[ii][jj][k].vx,
			&avery[ii][jj][k].vy,
			&avery[ii][jj][k].vz);
/*
fprintf(fpnotes,"i j k jj kk x %3d %3d %3d     %3d %3d     %f\n",i,j,k,ii,jj,avery[i][j][k].x);
fflush(fpnotes);
*/
			avery[ii][jj][k].vx *= -1.;
			avery[ii][jj][k].vy *= -1.;

		if (si_units){
			avery[i][j][k].x *= 100./PC;
			avery[i][j][k].y *= 100./PC;
			avery[i][j][k].z *= 100./PC;
			avery[ii][jj][k].density /= 1000.*2.33*AMU;
			avery[ii][jj][k].vx *= 100.;
			avery[ii][jj][k].vy *= 100.;
			avery[ii][jj][k].vz *= 100.;
		}
		if (n != EOF) {
			nobj += n;
		} else {
			printf("reached EOF at %d %d %d\n",i,j,k);
			break;
		}
	}}}

	printf("Finished reading data file\n");

	for (i=0;i<ngrid;i++){
	for (j=0;j<ngrid;j++){
	for (k=0;k<ngrid;k++){
		if (avery[i][j][k].x > maxx) maxx = avery[i][j][k].x;
		if (avery[i][j][k].y > maxy) maxy = avery[i][j][k].y;
		if (avery[i][j][k].z > maxz) maxz = avery[i][j][k].z;
		if (avery[i][j][k].vx > maxvx) maxvx = avery[i][j][k].vx;
		if (avery[i][j][k].vy > maxvy) maxvy = avery[i][j][k].vy;
		if (avery[i][j][k].vz > maxvz) maxvz = avery[i][j][k].vz;
		if (avery[i][j][k].density > maxdensity) {
			im = i;
			jm = j;
			km = k;	
			nm = nobj;
			maxdensity = avery[i][j][k].density;
		}
		if (avery[i][j][k].x < minx) minx = avery[i][j][k].x;
		if (avery[i][j][k].y < miny) miny = avery[i][j][k].y;
		if (avery[i][j][k].z < minz) minz = avery[i][j][k].z;
		if (avery[i][j][k].vx < minvx) minvx = avery[i][j][k].vx;
		if (avery[i][j][k].vy < minvy) minvy = avery[i][j][k].vy;
		if (avery[i][j][k].vz < minvz) minvz = avery[i][j][k].vz;
		if (avery[i][j][k].density < mindensity && 
			avery[i][j][k].density > 0.0) 
			mindensity = avery[i][j][k].density;
	}}}

	n = nobj/(ngrid*ngrid*ngrid);
	if (n < 1) {
		printf("Did not finish reading through model grid\n");
		printf("Finished %ld lines, which is %f percent of total\n",
			nobj,(float)nobj/(float)(ngrid*ngrid*ngrid)*100.);
		return 0;
	} else {
		printf("Successfully read Avery data\n");
	}

	printf("range x %f %f (pc)\n",minx,maxx);
	printf("range y %f %f (pc)\n",miny,maxy);
	printf("range z %f %f (pc)\n",minz,maxz);
	printf("range vx %f %f (kms)\n",minvx/1.e5,maxvx/1.e5);
	printf("range vy %f %f (kms)\n",minvy/1.e5,maxvy/1.e5);
	printf("range vz %f %f (kms)\n",minvz/1.e5,maxvz/1.e5);
	printf("range number density %g %g\n",mindensity, maxdensity);
	printf("density max i,j,k,n %ld %ld %ld %ld\n",im,jm,km,nm);

	return 1;

} /* end of read_avery function */

int read_ulrich(char filename[])
{

unsigned long mfsize(FILE *);
 
FILE *fplynds;
int   iplynds;
unsigned long  file_size;

int i,j,k,l,m,n,nx,ny,nz,ngrid;
long nobj,nbytes,nbytes_remain;
int  nsets,nbytes_per_set;
int my_rank;

float cell,rs,vk;

float float_var;
double mindensity=1.e20,maxdensity=-1.e20,
	minx=1.e20,maxx=-1.e20,
	miny=1.e20,maxy=-1.e20,
	minz=1.e20,maxz=-1.e20,
	minvx=1.e20,maxvx=-1.e20,
	minvy=1.e20,maxvy=-1.e20,
	minvz=1.e20,maxvz=-1.e20;
long im=0,jm=0,km=0,nm=0;

/*
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
printf("P%d: starting to read data for cloud \n",my_rank);
*/

fplynds = fopen(filename,"r");
if (fplynds == NULL) {
	printf("Could not open data file %s for cloud. Exiting ...\n",
		filename);
	fprintf(fpnotes,"Could not open data file %s for cloud. Exiting ...\n",
		filename);
        return 0;
} else {
	printf("Opened data file %s for cloud\n",
		filename);
	fprintf(fpnotes,"Opened data file %s for cloud\n",
		filename);
}
iplynds = fileno(fplynds);

file_size = mfsize(fplynds);
printf("size of input file is %lu\n",file_size);


 n = fread(&nx,sizeof(int),1,fplynds);
 n = fread(&ny,sizeof(int),1,fplynds);
 n = fread(&nz,sizeof(int),1,fplynds);
 n = fread(&cell,sizeof(float),1,fplynds);
 n = fread(&rs,sizeof(float),1,fplynds);
 n = fread(&vk,sizeof(float),1,fplynds);

 printf("nx,ny,nz = %d, %d, %d\n",nx,ny,nz);
 printf("cell, rs, vk = %e, %e, %e\n",cell,rs,vk);
 fprintf(fpnotes,"nx,ny,nz = %d, %d, %d\n",nx,ny,nz);
 fprintf(fpnotes,"cell, rs, vk = %e, %e, %e\n",cell,rs,vk);

/* allocate memory for the model, the size of the model grid will be the first integer
   in the data file */

    avery = (struct avery***) malloc(nx*sizeof(struct avery**));
    for (i=0;i<nx;i++){
      avery[i] = (struct avery **) malloc(ny*sizeof(struct avery *));
      for (j=0;j<ny;j++){
        avery[i][j] = (struct avery *) malloc(nz*sizeof(struct avery));
        }}

/* Loop here over cube */

nbytes = 0;
ngrid = 0;

	for (k=0;k<nz;k++){
	for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){

 		n = fread(&avery[i][j][k].vx,sizeof(float),1,fplynds);
 		n = fread(&avery[i][j][k].vy,sizeof(float),1,fplynds);
 		n = fread(&avery[i][j][k].vz,sizeof(float),1,fplynds);
 		n = fread(&avery[i][j][k].density,sizeof(float),1,fplynds);

		if (n != EOF) {
			ngrid += 1;
		} else {
			printf("reached EOF at %d %d %d\n",i,j,k);
			break;
		}

		if (avery[i][j][k].vx > maxvx) maxvx = avery[i][j][k].vx;
		if (avery[i][j][k].vy > maxvy) maxvy = avery[i][j][k].vy;
		if (avery[i][j][k].vz > maxvz) maxvz = avery[i][j][k].vz;
		if (avery[i][j][k].density > maxdensity) {
			im = i;
			jm = j;
			km = k;	
			nm = ngrid;
			maxdensity = avery[i][j][k].density;
		}
		if (avery[i][j][k].vx < minvx) minvx = avery[i][j][k].vx;
		if (avery[i][j][k].vy < minvy) minvy = avery[i][j][k].vy;
		if (avery[i][j][k].vz < minvz) minvz = avery[i][j][k].vz;
		if (avery[i][j][k].density < mindensity && 
			avery[i][j][k].density > 0.0) 
			mindensity = avery[i][j][k].density;
	}}}

	n = ngrid/(nx*ny*nz);
	if (n < 1) {
		printf("Did not finish reading through model grid\n");
		printf("Finished %d lines, which is %f percent of total\n",
			ngrid,(float)ngrid/(float)(nx*ny*nz)*100.);
		fprintf(fpnotes,"Did not finish reading through model grid\n");
		fprintf(fpnotes,"Finished %d lines, which is %f percent of total\n",
			ngrid,(float)ngrid/(float)(nx*ny*nz)*100.);
		return 0;
	} else {
		printf("Appears to have successfully read Ulrich data\n");
		fprintf(fpnotes,"Appears to have successfully read Ulrich data\n");
	}

	printf("range vx %f %f (kms)\n",minvx/100000.,maxvx/100000.);
	printf("range vy %f %f (kms)\n",minvy/100000.,maxvy/100000.);
	printf("range vz %f %f (kms)\n",minvz/100000.,maxvz/100000.);
	printf("range density %g %g\n",mindensity, maxdensity);
	printf("density max i,j,k %ld %ld %ld\n",im,jm,km);
	fprintf(fpnotes,"range vx %f %f (kms)\n",minvx/100000.,maxvx/100000.);
	fprintf(fpnotes,"range vy %f %f (kms)\n",minvy/100000.,maxvy/100000.);
	fprintf(fpnotes,"range vz %f %f (kms)\n",minvz/100000.,maxvz/100000.);
	fprintf(fpnotes,"range density %g %g\n",mindensity, maxdensity);
	fprintf(fpnotes,"density max i,j,k %ld %ld %ld \n",im,jm,km);

	return 1;

} /* end of read_ulrich function */

int read_lynds(char filename[], int ntp) 
{

unsigned long mfsize(FILE *);
 
FILE *fplynds;
int   iplynds;
unsigned long  file_size;

int i,j,k,l,m,n;
int nobj;
long nbytes,nbytes_remain;
int  nsets,nbytes_per_set;
int my_rank;

float float_var;

/*
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
printf("P%d: starting to read data for cloud \n",my_rank);
*/

fplynds = fopen(filename,"r");
if (fplynds == NULL) {
	printf("Could not open data file %s for Lynds clouds. Exiting ...\n",
		filename);
        return 0;
} else {
	printf("Opened data file %s for Lynds clouds\n",
		filename);
}
iplynds = fileno(fplynds);

file_size = mfsize(fplynds);
printf("size of input file is %lu\n",file_size);

nbytes = 0;

nobj = fread(&lynds.lambda,sizeof(float),1,fplynds); 
nobj = fread(&lynds.scale_alpha,sizeof(float),1,fplynds); 
nobj = fread(&lynds.scale_beta,sizeof(float),1,fplynds); 
nobj = fread(&lynds.bnd_temp,sizeof(float),1,fplynds);
nobj = fread(&float_var,sizeof(float),1,fplynds);
lynds.total_mass = float_var*MSUN;
nobj = fread(&lynds.pext,sizeof(float),1,fplynds);
nobj = fread(&lynds.boundary,sizeof(float),1,fplynds);
nobj = fread(&lynds.xn,sizeof(float),1,fplynds);
nobj = fread(&hydro.jx0,sizeof(int),1,fplynds);
nobj = fread(&lynds.nsteps,sizeof(int),1,fplynds);
nobj = fread(&lynds.nt,sizeof(int),1,fplynds);

nbytes += 8*sizeof(float);
nbytes += 3*sizeof(int);

if (hydro.jx0 > 30000) {
        printf("jx0 %d looks wrong.\n",hydro.jx0);
        printf("The file probably has the wrong byte order.\n");
        fprintf(fpnotes,"jx0 %d looks wrong.\n",hydro.jx0);
        fprintf(fpnotes,"The file probably has the wrong byte order.\n");
	return 0;
}
/*
nobj = fread(&lynds.,sizeof(int),1,fplynds);
nobj = fread(&lynds.,sizeof(int),1,fplynds);

nobj = fread(&lynds.,sizeof(float),1,fplynds);
nobj = fread(&lynds.,sizeof(float),1,fplynds);
*/

printf("lambda = %e\n",lynds.lambda/(2.33*AMU));
printf("scale_alpha = %e\n",lynds.scale_alpha);
printf("scale_beta = %e\n",lynds.scale_beta);
printf("bnd_temp = %f\n",lynds.bnd_temp);
printf("total_mass = %e\n",lynds.total_mass/MSUN);
printf("pext = %f\n",lynds.pext/BOLTZ);
printf("boundary = %f\n",lynds.boundary*lynds.scale_alpha/PC);
printf("xn = %f\n",lynds.xn);
printf("jx0 = %d\n",hydro.jx0);
printf("nsteps = %d\n",lynds.nsteps);
printf("nt = %d\n",lynds.nt);

lynds.x 	= (float *)malloc(lynds.nsteps*sizeof(float));
lynds.s 	= (float *)malloc(lynds.nsteps*sizeof(float));
lynds.z 	= (float *)malloc(lynds.nsteps*sizeof(float));
lynds.th	= (float *)malloc(lynds.nsteps*sizeof(float));
lynds.tgas	= (float *)malloc(lynds.nt*sizeof(float));
lynds.tdust	= (float *)malloc(lynds.nt*sizeof(float));
hydro.sn	= (float *)malloc((hydro.jx0+1)*sizeof(float));
hydro.mdn0	= (float *)malloc((hydro.jx0+1)*sizeof(float));
hydro.pnt0	= (float *)malloc((hydro.jx0+1)*sizeof(float));

for (i=0;i<lynds.nsteps;i++) nobj = fread(&lynds.x[i], sizeof(float),1,fplynds);
for (i=0;i<lynds.nsteps;i++) nobj = fread(&lynds.s[i], sizeof(float),1,fplynds);
for (i=0;i<lynds.nsteps;i++) nobj = fread(&lynds.z[i], sizeof(float),1,fplynds);
for (i=0;i<lynds.nsteps;i++) nobj = fread(&lynds.th[i],sizeof(float),1,fplynds);
for (i=0;i<lynds.nt;i++) nobj = fread(&lynds.tgas[i], sizeof(float),1,fplynds);
for (i=0;i<lynds.nt;i++) nobj = fread(&lynds.tdust[i],sizeof(float),1,fplynds);

nbytes += 4*lynds.nsteps*sizeof(float);
nbytes += 2*lynds.nt*sizeof(float);

printf("ND radius at center %f  at boundary %f\n",
	lynds.x[0],lynds.x[lynds.nt-1]);
printf("tgas  at center %f  at boundary %f\n",
	lynds.tgas[0],lynds.tgas[lynds.nt-1]);
printf("tdust at center %f  at boundary %f\n",
	lynds.tdust[0],lynds.tdust[lynds.nt-1]);

nobj = fread(&hydro.shock_smoothing,sizeof(float),1,fplynds);
printf("shock_smoothing = %f\n",hydro.shock_smoothing);

nbytes += sizeof(float);

for (i=0;i<hydro.jx0+1;i++) nobj=fread(&hydro.sn[i],  sizeof(float),1,fplynds);
for (i=0;i<hydro.jx0+1;i++) nobj=fread(&hydro.mdn0[i],sizeof(float),1,fplynds);
nbytes += 2*(hydro.jx0+1)*sizeof(float);
if (ntp){ 
    for (i=0;i<hydro.jx0+1;i++) 
        nobj=fread(&hydro.pnt0[i],sizeof(float),1,fplynds);
    nbytes +=   (hydro.jx0+1)*sizeof(float);
}

printf("density at center %e  at boundary %e\n",
	hydro.mdn0[0],hydro.mdn0[lynds.nt-1]);

printf("finished reading %ld bytes from start of file\n",nbytes);

nbytes_remain = file_size - nbytes;
printf("there are %ld bytes left in the file\n",nbytes_remain);

/* Look below at the list of mallocs. There are 3 scalars, et,dt,pext,
and there are 14 vectors. One of these is the non-thermal pressure, pntn,
that is not always written. So there are 13 always written and 1 more
if the hydro file was written with non-thermal pressure. Now calculate
the number of bytes per set for regular files and add the length
of the extra vector for files with non-thermal pressure */

nbytes_per_set = 3*sizeof(float) + 13*(hydro.jx0+1)*sizeof(float);;
if (ntp) nbytes_per_set += (hydro.jx0+1)*sizeof(float);
printf("there are %d bytes in each set\n",nbytes_per_set);

nsets = nbytes_remain/nbytes_per_set;
if (nsets < 1) nsets = 1;

printf("Expecting %d sets in the evolution\n",nsets);
if (nsets > 1) {
	printf("Drop the last set because might have nan\n");
	nsets -= 1;
}
printf("Will read %d sets of data\n",nsets);
hydro.nsets = nsets;

hydro.et 	=  (float *)malloc((hydro.jx0+1)*sizeof(float));
hydro.dt 	=  (float *)malloc((hydro.jx0+1)*sizeof(float));
hydro.pext 	=  (float *)malloc((hydro.jx0+1)*sizeof(float));
hydro.rn 	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.mdn 	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.tn 	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.pn 	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.pntn 	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.vn 	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.gn 	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.qn 	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.dpn 	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.dqn 	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.tdust 	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.co_abundance	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.avgavn	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.gpe	= (float **)malloc((hydro.jx0+1)*sizeof(float *));
hydro.h2o_abundance	= (float **)malloc((hydro.jx0+1)*sizeof(float *));

for (i=0;i<hydro.jx0+1;i++) {

    hydro.rn[i] 	= (float *)malloc(nsets*sizeof(float));
    hydro.mdn[i]	= (float *)malloc(nsets*sizeof(float));
    hydro.tn[i]		= (float *)malloc(nsets*sizeof(float));
    hydro.pn[i]		= (float *)malloc(nsets*sizeof(float));
    hydro.pntn[i]	= (float *)malloc(nsets*sizeof(float));
    hydro.vn[i]		= (float *)malloc(nsets*sizeof(float));
    hydro.gn[i]		= (float *)malloc(nsets*sizeof(float));
    hydro.qn[i]		= (float *)malloc(nsets*sizeof(float));
    hydro.dpn[i]	= (float *)malloc(nsets*sizeof(float));
    hydro.dqn[i]	= (float *)malloc(nsets*sizeof(float));
    hydro.tdust[i]	= (float *)malloc(nsets*sizeof(float));
    hydro.co_abundance[i]	= (float *)malloc(nsets*sizeof(float));
    hydro.avgavn[i]	= (float *)malloc(nsets*sizeof(float));
    hydro.gpe[i]	= (float *)malloc(nsets*sizeof(float));
    hydro.h2o_abundance[i]	= (float *)malloc(nsets*sizeof(float));

}


for (i=0;i<nsets;i++){

    fread(&hydro.et[i],sizeof(float),1,fplynds);
    printf("hydro.et %g\n",hydro.et[i]); 
    if (i == 0 && hydro.et[i] < 10.) {
        printf("Elapsed time %e looks wrong.\n",hydro.et[i]);
        printf("The nonthermal pressure switch is probably not set correctly.  \n");
        printf("Exiting now ... \n");
	return 0;
    }
    printf("%d et %e \n",i,hydro.et[i]); 
    fread(&hydro.dt[i],sizeof(float),1,fplynds);
    printf("hydro.dt %e \n",hydro.dt[i]); 
    fread(&hydro.pext[i],sizeof(float),1,fplynds);
    printf("hydro.pext %e \n",hydro.pext[i]); 

    for(j=0;j<hydro.jx0+1;j++){
	nobj = fread(&hydro.rn[j][i],sizeof(float),1,fplynds);
/*        printf("hydro.rn %e\n",hydro.rn[j][i]); */
    }
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.mdn[j][i],sizeof(float),1,fplynds);
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.tn[j][i],sizeof(float),1,fplynds);
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.pn[j][i],sizeof(float),1,fplynds);
    if (ntp) {
        for(j=0;j<hydro.jx0+1;j++)
	    nobj = fread(&hydro.pntn[j][i],sizeof(float),1,fplynds);
    }
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.vn[j][i],sizeof(float),1,fplynds);
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.gn[j][i],sizeof(float),1,fplynds);
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.qn[j][i],sizeof(float),1,fplynds);
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.dpn[j][i],sizeof(float),1,fplynds);
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.dqn[j][i],sizeof(float),1,fplynds);
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.tdust[j][i],sizeof(float),1,fplynds);
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.co_abundance[j][i],sizeof(float),1,fplynds);
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.avgavn[j][i],sizeof(float),1,fplynds);
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.gpe[j][i],sizeof(float),1,fplynds);
    for(j=0;j<hydro.jx0+1;j++)
	nobj = fread(&hydro.h2o_abundance[j][i],sizeof(float),1,fplynds);

    if (feof(fplynds) || ferror(fplynds)) {
        printf("EOF or ERROR in readings sets of evolution\n");
        break;
    }


}

printf("radius first set at center %e at boundary %e\n",
	hydro.rn[0][0],hydro.rn[hydro.jx0][0]);
printf("gas temperature first set at center %f at boundary %f\n",
	hydro.tn[0][0],hydro.tn[hydro.jx0][0]);
printf("gas density first set at center %e at boundary %e\n",
	hydro.mdn[0][0],hydro.mdn[hydro.jx0][0]);
printf("gas velocity first set at center %f at boundary %f\n",
	hydro.vn[0][0],hydro.vn[hydro.jx0][0]);

printf("Total elapsed time in data set %e\n",hydro.et[nsets-1]/YR);

return 1;

} /* end of read_lynds function */

int read_data()
{
  FILE *fpdata ;
  int i,j,k,l,m,n,p,line;
  int nobj;
  float fwhm2;
  char fname[64];
  char textline[LINE_MAX];

  strcpy(fname,"./nh3_11line_4arc_new.dat");
  fpdata = fopen(fname,"r");
  printf("Opening spectral line data file %s\n",fname);
  if (fpdata == NULL) {
          printf("problem opening data file\n");
          printf("Program transferring to exit ...\n");
          return 0;
  }

/* Read in the NH3(1,1) line, 62 channels */
    fwhm2 = 4.0*4.0;

/* skip 5 lines of text */
  for (i=0;i<5;i++) fgets(textline,LINE_MAX,fpdata);

  for(i=0;i<62;i++){
  j = 61 - i;

    nobj = fscanf (fpdata, "%f %f %f %f ",&data.vel[j][0],&data.data[j][0][0],&data.data[j][1][0],&data.data[j][2][0]);


    for (m=0;m<3;m++) data.data[j][m][0] = data.data[j][m][0] * pow(30./23.694495,2)
       / (2.0*BOLTZ*fwhm2*2.663e-11)*1.e-23;



  }
  fclose(fpdata);

/* Read in the NH3(2,2) line, 63 channels */
    fwhm2 = 4.0*4.0;

  strcpy(fname,"./nh3_22line_4arc_new.dat");
  fpdata = fopen(fname,"r");
  printf("Opening spectral line data file %s\n",fname);
  if (fpdata == NULL) {
          printf("problem opening data file\n");
          printf("Program transferring to exit ...\n");
          return 0;
  }

/* skip 5 lines of text */
  for (i=0;i<5;i++) fgets(textline,LINE_MAX,fpdata);

  for(i=0;i<62;i++){
  j = 61 - i;

    nobj = fscanf (fpdata, "%f %f %f %f ",&data.vel[j][1],&data.data[j][0][1],&data.data[j][1][1],&data.data[j][2][1]);


    for (m=0;m<3;m++) data.data[j][m][1] = data.data[j][m][1] * pow(30./23.722633000,2)
       / (2.0*BOLTZ*fwhm2*2.663e-11)*1.e-23;




  }
  fclose(fpdata);

/* Read in the NH3(3,3) line, 58 channels*/
    fwhm2 = 4.0*4.0;

  strcpy(fname,"./nh3_33line.dat");
  fpdata = fopen(fname,"r");
  printf("Opening spectral line data file %s\n",fname);
  if (fpdata == NULL) {
          printf("problem opening data file\n");
          printf("Program transferring to exit ...\n");
          return 0;
  }

/* skip 6 lines of text */
  for (i=0;i<6;i++) fgets(textline,LINE_MAX,fpdata);

  for(i=0;i<58;i++){
  j = 57 - i;

    nobj = fscanf (fpdata, "%f %f %f %f ",&data.vel[j][2],&data.data[j][0][2],&data.data[j][1][2],&data.data[j][2][2]);

    for (m=0;m<3;m++) data.data[j][m][2] = data.data[j][m][2] * pow(30./23.870296000,2)
       / (2.0*BOLTZ*fwhm2*2.663e-11)*1.e-23;



  }
  fclose(fpdata);
  data.vel[58][2] = data.vel[57][2] + (data.vel[57][2] - data.vel[56][2]);
  data.vel[59][2] = data.vel[58][2] + 2.*(data.vel[57][2] - data.vel[56][2]);


/* Read in the NH3(4,4) line, 57 channels */
    fwhm2 = 4.0*4.0;

  strcpy(fname,"./nh3_44line.dat");
  fpdata = fopen(fname,"r");
  printf("Opening spectral line data file %s\n",fname);
  if (fpdata == NULL) {
          printf("problem opening data file\n");
          printf("Program transferring to exit ...\n");
          return 0;
  }

/* skip 6 lines of text */
  for (i=0;i<6;i++) fgets(textline,LINE_MAX,fpdata);

  for(i=0;i<57;i++){
  j = 56 - i;

    nobj = fscanf (fpdata, "%f %f %f %f ",&data.vel[j][3],&data.data[j][0][3],&data.data[j][1][3],&data.data[j][2][3]);


    for (m=0;m<3;m++) data.data[j][m][3] = data.data[j][m][3] * pow(30./24.139390000,2)
       / (2.0*BOLTZ*fwhm2*2.663e-11)*1.e-23;


  }
  fclose(fpdata);
  data.vel[57][3] = data.vel[56][3] +    (data.vel[56][3] - data.vel[55][3]);
  data.vel[58][3] = data.vel[56][3] + 2.*(data.vel[56][3] - data.vel[55][3]);
  data.vel[59][3] = data.vel[56][3] + 3.*(data.vel[56][3] - data.vel[55][3]);


  printf("Spectral line data read in OK\n");

    for (line=0;line<4;line++){
    fprintf(fpnotes,"data for line %d\n",line);
    for (i=0;i<60;i++){
        fprintf(fpnotes,"i vel I %4d %7.3f %7.3f %7.3f %7.3f\n",
          i,data.vel[i][line],data.data[i][0][line],data.data[i][1][line],data.data[i][2][line]); 
    }}

/* now we have the data and the rms of each line and position */

  return 1;

}





float chisq_function( )
{

        int i,j,k,n,m,ibest,iloc,ic,ipos[3],jc,iv,line,modl,ipspectra,factor;
        float chisqmin,chisq,voffset,vm,vbest,xpos[3],aposx[3],aposy[3];
	float model_velocities[2000];
	float spectra[2][60][3][4];
	float velocity[60][4];
	char fname[64],  suffix[4], pathname[64];
	float x2,y2;

/* spectra[model][channel][position][line] */

/*
printf("nchan %d \n",dcube.nchan[0]);
printf("starting velocity data  %e\n",data.vel[0]);
printf("starting velocity model %e\n",dcube.chvel[0][0]);
printf("data  channel width %e\n",data.vel[0] - data.vel[1]);
printf("model channel width %e\n",dcube.chvel[0][0] - dcube.chvel[0][1]);
*/

/* there are spectra at 3 positions for each line.
   The data are in the plane of the disk so j = 0.
   Here are the X coordinates of the 3 positions in pc. */

        strcpy(pathname,"./");
        strcpy(fname,pathname);
        strcat(fname,"check_chisq");
        sprintf(suffix,"%d",parameters.model_number);
        strcat(fname,suffix);

       fpspectra = fopen(fname,"w");
             if (fpspectra == NULL) {
                printf("Could not open %s. Exiting ...\n",fname);
                fflush(stdout);
             }
       printf("File %s opened\n",fname);
       ipspectra = fileno(fpspectra);

/* These are the numbers I calculated by hand.
   The calculation below is identical */
	xpos[0] = -0.0085;
	xpos[1] = -0.0364;
	xpos[2] = +0.0221;

/* These are the offsets from the data file */
	aposx[0] = -0.23;
	aposy[0] =  1.00;
	aposx[1] = -3.10;
	aposy[1] =  3.15;
	aposx[2] =  2.50;
	aposy[2] = -0.94;

	for (i=0;i<3;i++) {
		x2 = pow( aposx[i] - aposx[0]  ,2);
		y2 = pow( aposy[i] - aposy[0]  ,2);
		xpos[i] = sqrt(x2 + y2);
	fprintf(fpnotes,"radius arcsec %10.4f\n",xpos[i]);
		xpos[i] = 1700.*tan(xpos[i]/3600.*PI/180.);
	fprintf(fpnotes,"radius pc     %10.4f\n",xpos[i]);
		ipos[i] = xpos[i]/dcube.cellx;
	fprintf(fpnotes,"radius cells  %5d\n",ipos[i]);
	}



/* points 1 and 2 are on the left side */
	xpos[0] = -xpos[0];
	xpos[1] = -xpos[1];
	
	for (i=0;i<3;i++) {
		ipos[i] = xpos[i]/dcube.cellx + dcube.noutx/2 - 1;
		printf("X pos %8.4f I pos %d \n",xpos[i],ipos[i]);
		fprintf(fpnotes,"X pos %8.4f I pos %d \n",xpos[i],ipos[i]);
	}


        /* Slide the model spectra over the data spectra to find the velocity
         *    that fits the best. */
	chisqmin = 1.e20;

  for (iv=0;iv<50;iv++){
		voffset = -4.00 + (iv-25)*0.05;

  chisq = 0.;

/* --------------------------------------------------------------------------------------------- */
/* Loop over all lines, channels, positions */

	for (line=0;line<4;line++){
	factor = 1.0;
	if (line == 3) factor = 4.0;
	n = dcube.nchan[line];
	for (k=0;k<n;k++) 
             model_velocities[k]=dcube.chvel[line][k]/1.e5;

	for (i=0;i<3;i++) {
	ic = ipos[i];
	jc = dcube.nouty/2 - 1;
		for (k=0;k<60;k++) {
		   vm = data.vel[k][line] - voffset;
		   iloc = locate_(model_velocities,&n,&vm,&m); 
/* for model velocity vm, k is the index of the data and m is the index of the model */
                   if (iloc >= 0 && iloc < n)
                   chisq = chisq + factor*pow((dcube.cube[line][0][ic][jc][m] - data.data[k][i][line]) ,2) ;
		}
	}
	}

	if (chisq < chisqmin) {
		chisqmin = chisq;
		ibest = iv;
		vbest = voffset;
	}
  }

	printf("iv best %d voffset %f\n",ibest,vbest);

        for (line=0;line<4;line++){
        n = dcube.nchan[line];
        for (k=0;k<n;k++) 
             model_velocities[k]=dcube.chvel[line][k]/1.e5;

        for (i=0;i<3;i++) {
        ic = ipos[i];
	jc = dcube.nouty/2 - 1;
                for (k=0;k<60;k++) {
                   vm = data.vel[k][line] - vbest;
		   velocity[k][line] = vm;
                   iloc = locate_(model_velocities,&n,&vm,&m);
/* for model velocity vm, k is the index of the data and m is the index of the model */
                   if (iloc >= 0 && iloc < n)
                   spectra[0][k][i][line] = data.data[k][i][line];
                   spectra[1][k][i][line] = dcube.cube[line][0][ic][jc][m];
                }
        }
        }

     modl = 0;
     for (k=0;k<60;k++) {
	fprintf(fpspectra,"%5.1f %8.4f %8.4f %8.4f %5.1f %8.4f %8.4f %8.4f %5.1f %8.4f %8.4f %8.4f %5.1f %8.4f %8.4f %8.4f\n",
	velocity[k][0],spectra[modl][k][0][0],spectra[modl][k][1][0],spectra[modl][k][2][0],
	velocity[k][1],spectra[modl][k][0][1],spectra[modl][k][1][1],spectra[modl][k][2][1],
	velocity[k][2],spectra[modl][k][0][2],spectra[modl][k][1][2],spectra[modl][k][2][2],
	velocity[k][3],spectra[modl][k][0][3],spectra[modl][k][1][3],spectra[modl][k][2][3]);
     }
     modl = 1;
     for (k=0;k<60;k++) {
	fprintf(fpspectra,"%5.1f %8.4f %8.4f %8.4f %5.1f %8.4f %8.4f %8.4f %5.1f %8.4f %8.4f %8.4f %5.1f %8.4f %8.4f %8.4f\n",
	velocity[k][0],spectra[modl][k][0][0],spectra[modl][k][1][0],spectra[modl][k][2][0],
	velocity[k][1],spectra[modl][k][0][1],spectra[modl][k][1][1],spectra[modl][k][2][1],
	velocity[k][2],spectra[modl][k][0][2],spectra[modl][k][1][2],spectra[modl][k][2][2],
	velocity[k][3],spectra[modl][k][0][3],spectra[modl][k][1][3],spectra[modl][k][2][3]);
     }



	fflush(fpspectra);
  	fsync(ipspectra);
        fclose(fpspectra);
        printf("Finished writing plotting data\n");


	return chisqmin;

}

int set_model_parameters(int n, float *p, float **plim, int model_number, int my_rank)
{

int i,j;
float del[MAXADJP],p0[MAXADJP];

parameters.model_number = model_number;

/* Here are the initial guesses for the adjustable model parameters */

/* env density */
p0[0] = log10(pow(10.,8.701));
/* point mass */
p0[1] = 12.7  ;
/* Abundance */
 p0[2] = log10(pow(10.,-7.246));
/* gamma0 = Vr * Rr rotational velocity = 2.e5, rotational radius = 2250.*AU
   is a factor to multiply this standard value */
p0[3] = 7379.;
/* planck mean exponent */
p0[4] = -0.1;
/* disk density factor*/
/* p0[5] = 17.19; */
/* disk temperature factor*/
/* p0[6] = log10(1.533); */
/* width     */
/*p0[7] = 0.500 ; */
/* env  temperature factor*/
/* p0[8] = log10(1.000); */

  parameters.turbwidth  =        1.0;
  parameters.env_temperature_factor = 1.0;
  parameters.disk_temperature_factor = 1.0;
  parameters.disk_density_factor = 1.0;

/* These are the initial changes to be made to the adjustable
   model parameters. */

del[0] = 1.0;
del[1] = 5.0;
del[2] = 0.5;
del[3] = 200.;
del[4] = 0.2;
/*
del[5] = 5.0;
del[6] = 0.2;
del[7] = 0.5;
del[8] = 0.2;
*/

/* These are the limits between which the model parameters are
   allowed to be adjusted. */

plim[0][0] = log10(1.e7);
plim[0][1] = log10(1.e12);
plim[1][0] = 10.;
plim[1][1] = 30.;
plim[2][0] = log10(1.e-10);
plim[2][1] = log10(1.e-7);
plim[3][0] = 1000.  ;
plim[3][1] = 1.e4 ;
plim[4][0] = -1.0;
plim[4][1] = 2.0;
/*
plim[5][0] = 0.01;
plim[5][1] = 1000.0;
plim[6][0] = -1.0;
plim[6][1] = 1.0;
plim[7][0] = 1.0;
plim[7][1] = 5.0;
plim[8][0] = -1.0;
plim[8][1] = 1.0;
*/

/* This scales the limits so that the code can use them */

for (i=0;i<n;i++) {
for (j=0;j<2;j++) {
/* printf("%d Limits: old %f",i,plim[i][j]); */
  fprintf(fpnotes,"%d Limits: old %f",i,plim[i][j]);
  plim[i][j] = (plim[i][j] - p0[i])/del[i];
/*  printf("  new %f\n",plim[i][j]); */
  fprintf(fpnotes,"  new %f\n",plim[i][j]);
}}

/* The model parameters are adjusted here */



  parameters.env_density	= pow(10.,p0[0] + del[0]*p[0]);
  parameters.point_mass	=         p0[1] + del[1]*p[1];
  parameters.abundance  = pow(10.,p0[2] + del[2]*p[2]);
  parameters.gamma0  =            p0[3] + del[3]*p[3];
  parameters.planck_exp = 	p0[4] + del[4]*p[4];
/*
  parameters.disk_temperature_factor = pow(10.,p0[5] + del[5]*p[5]);
  parameters.disk_density_factor= p0[6] + del[6]*p[6] ;
  parameters.turbwidth  =        (p0[7] + del[7]*p[7]);
  parameters.env_temperature_factor = pow(10.,p0[8] + del[8]*p[8]);
*/


  if (my_rank == 0) {
      printf("P%d: Define a new model for simulation %d \n",
                  my_rank,model_number);
      fprintf(fpnotes,"P%d: Define a new model for simulation %d \n",
                  my_rank,model_number);
      printf("Density is now %e\n",parameters.env_density);
      printf("Disk factor is now %f\n",parameters.disk_density_factor);
      printf("Point mass is now %f\n",parameters.point_mass);
      printf("Abundance is now %e\n",parameters.abundance);
      printf("Gamma0    is now %f\n",parameters.gamma0);
      printf("Width     is now %f\n",parameters.turbwidth);
      printf("Disk T fact  now %f\n",parameters.disk_temperature_factor);
      printf("Env  T fact  now %f\n",parameters.env_temperature_factor);
      printf("Planck exp   now %f\n",parameters.planck_exp);
      fflush(stdout);

      if (model_number == 0)
/*                            12345678 12345678 12345678 12345678 12345678 12345678  */
              fprintf(fphist,"  Model Density    DiskD  Mass     Abundanc Gamma      Width    DiskT      EnvT      Planck     X^2\n");
      fprintf(fphist,"%6d %8.3f %8.2f %8.4f %8.3f %8.3f %8.3f %8.3f %8.3f %10.3f ",
                      model_number,
                      log10(parameters.env_density),
                      parameters.disk_density_factor,
                      parameters.point_mass,
                      log10(parameters.abundance),
                      parameters.gamma0, 
                      parameters.turbwidth,
		      parameters.disk_temperature_factor,
		      parameters.env_temperature_factor,
		      parameters.planck_exp);

  }

  fprintf(fpnotes,"Density is now %e\n",parameters.env_density);
  fprintf(fpnotes,"Disk factor is now %f\n",parameters.disk_density_factor);
  fprintf(fpnotes,"Point mass is now %f\n",parameters.point_mass);
  fprintf(fpnotes,"Abundance is now %e\n",parameters.abundance);
  fprintf(fpnotes,"Gamma0    is now %f\n",parameters.gamma0);
  fprintf(fpnotes,"Width     is now %f\n",parameters.turbwidth);
  fprintf(fpnotes,"Disk T fact  now %f\n",parameters.disk_temperature_factor);
  fprintf(fpnotes,"Env  T fact  now %f\n",parameters.env_temperature_factor);
  fprintf(fpnotes,"Planck exp   now %f\n",parameters.planck_exp);



  return 1;
}

/*
int locate(float* xx, int n, float x, int* j)
{

int jl,ju,jm;

       jl = 0;
       ju = n+1;
 10    if (ju-jl > 1) then
           jm = (ju+jl)/2
           if ((xx(n) > xx(1)) == (x > xx(jm))) then
               jl = jm
           else
               ju = jm
           endif
           goto 10
       endif
       j = jl


       return 1;
}
*/



/* Use the functions in the math library
   instead of these

float acosh(float);
float asinh(float);

float acosh( float x )
{
 y = log( fabs(x) + sqrt( x*x - 1.0) );

 return y;

}

float asinh( float x )
{

 y = log( fabs(x) + sqrt( x*x + 1.0) )

 return y;

}
*/


int ulrich_flow(int i, int j, int k, float r, float el,
	float* vr, float* vaz, float* vel, float* den)
{

/* these algorithms were not designed for el > 90 degrees,
*/

double f1,f2,f3,f4,f5,f6,g1,g2,g3,g4,g5,g6,p2;
double cosel,sinel,cosel0,sinel0;

    cosel = cos(el);
    sinel = sin(el);

/* This is the singular point at r = rd */
    if ( fabs(r-1.0) < 1.e-6)  {
	g5 = 0.;
	cosel0 = pow(cosel,(1./3.));
        goto label1;
    }

/* Outside of Rd */
    if (r > 1.0) {

/* This is the special case of theta = 90 */
	if ( fabs(el - PI/2.) < 1.e-6 ) {
		cosel0 = 0.0;
		goto label1;
	}

/* For all other angles and theta not 90 degrees, still outside of Rd */
	g1 = (r - 1.0)/3.0;
	g2 = sqrt(g1);
	g3 = pow(g2,3);
	g4 = r*cosel / (2.*g3);
	g5 = 0.;
	cosel0 = 2.0*g2 * sinh( 1./3.* asinh( g4 ) );
	goto  label1;
    }

/* only get here if r < 1, inside of Rd */

/* This is the special case of theta = 90 degrees, plane of the disk*/
	if (fabs(el - PI/2.) < 1.e-6 ) {
/*		cosel0 = sqrt(1.0-r);  */
/* This doesn't work because vel != 0 in the plane of the disk,
but the line below does. Maybe the limit of theta = 90. */
		cosel0 = 0.0; 
		goto  label1;
	}

        g1 = (1.0 - r)/3.0;
        g2 = sqrt(g1);
        g3 = pow(g2,3);
        g4 = r*cosel / (2.*g3);
	g5 = pow((r/2.0*cosel),2) - pow(g3,2);
	if (g5 > 0.) {
        	cosel0 = 2.0*g2 * cosh( 1./3.* acosh( g4 ) );
	} else {
        	cosel0 = 2.0*g2 * cos( 1./3.* acos( g4 ) );
	}

label1: 
/*	printf("r,el*180./PI,g5,cosel0 %10.3e %10.3e %10.3e %10.3e \n",
		r,el*180./PI,g5,cosel0); 
*/

/* watch out for round-off error when el=0 */
	if (fabs(el) < 1.e-6 ) { 
		sinel0 = 0.0;
	} else {
		sinel0 = sin(acos(cosel0));
	}

/* calculate the density except where r = 0, which it never is. */
    if (r != 0.) {
        f1 = sqrt(1./r);
	f4 = pow(r,(-3./2.));
	f5 = 2./r;
    } else {
	f1 = 0.0;
	f4 = 0.0;
	f5 = 0.0;
    }

/* In the line above, the absolute value function handles
   the case with cosel = cosel0 = 1.0, but because of
   round off error the arg 1. - cosel/cosel0 is negative
*/

    if (fabs(cosel0) > 1.e-6) {
	f2 = sqrt( fabs(1. + cosel/cosel0) );
	p2 = 0.5*(3.*pow(cosel0,2) - 1.);
	f6 = sqrt( fabs(1. - cosel/cosel0) ); 
    } else {
	if (r < 1.) {
	    f2 = 1.;
	    f6 = 1.;
	    p2 = 0.5*(2.-3.*r);
	} else {
	    f2 = sqrt(2. - 1./r);
	    f6 = sqrt(1./r);
	    p2 = -0.5;
	}
    }


	*den = f4 * 1./f2 * 1./( 1. + f5 * p2);

/*
if (i == 10 && j == 10) printf("%3d %3d %3d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",
i,j,k,*den,r,f2,f4,f5,p2,el,cosel,cosel0);
*/


    if (fabs(sinel) > 1.e-6) {
	f3 = ((cosel0 - cosel)/sinel);
	f4 = sinel0/sinel;
/* if (k==12 && j==12) printf("A %10.3e %10.3e %10.3e %10.3e \n",f3,cosel0,cosel,sinel); */
    } else {
        f3 = 0.0;
	f4 = 0.0;
/* if (k==12 && j==12) printf("B %10.3e \n",f3); */
    }

    *vr  = -f1 * f2;
    *vel = f1 * f2 * f3;
    *vaz = f1 * f6 * f4;

/*
if (k==12 && j==12) printf("%d %d r %10.3e f1 %10.3e f2 %10.3e f3 %10.3e vr %10.3e vel %10.3e\n",
i,k,r,f1,f2,f3,*vr,*vel); 
*/


    return 1;

}


float flared_disk_density(int i, int j, int k, float x, float z, float rstar, float rd)
{

/* 

Equation 3 of Whitney 2003 for the density of a flared disk

input coordinates 
	x, z 			cylindrical coords, non-dim in units of rstar 

fixed disk parameters
	hscale0			scale height of the disk at rstar in
					non-dim units of rstar
	disk_alpha		power law exponent for the density
	disk_beta		power law exponent for the flare

output variables:
	disk_density		gas density in non-dim units of density

*/

float disk_alpha, disk_beta, hscale, hscale0, expzh, rpowr, inbnd;

	x = fabs(x);
	z = fabs(z);

	if (x < rstar ) {
/*		printf("x z %e %e return 0\n",x,z); */
		return 0;
	}

/* disk model parameters */
	disk_alpha 	= 2.25;
	disk_beta  	= 1.25;
	hscale0		= 0.01 * rstar;

/* The scale height is a function of cylindrical radius */
	hscale = hscale0 * pow(x/rstar,disk_beta);


/* The density is a function of vertical scale height and
   a radial power law */
	expzh = exp( -0.5*pow((z/hscale),2) ) ;
	rpowr = pow((rstar/x),disk_alpha);
	inbnd = (1. - sqrt(rstar/x));
/*
if (k == 12 && j == 12) printf("x %10.3e hscale %10.3e z %10.3e expzh %10.3e rpowr %10.3e inbnd %10.3e dens %10.3e\n",
	x/rd,hscale/rd,z/rd,expzh,rpowr,inbnd,expzh*rpowr*inbnd); 
*/

/* printf("expzh %e rpowr %e inbnd %e\n",expzh,rpowr,inbnd); */

	return inbnd * rpowr * expzh ;


}


int smooth ( int nx, int ny, int nz, 
		float beamx, float beamy,
		float cellx, float celly,
		float * xcenter, float * ycenter,
		float *** input_cube, float *** smoothed_cube )
{

  float xcent, ycent;
  float s2x, s2y, area, dx, dy;
  int inverse, status;
  int i,j,k;

  float **beamR,    **beamI;
  float **beamRfft, **beamIfft;
  float **mapR,     **mapI;
  float **mapRfft,  **mapIfft;
  float **convR,    **convI;
  float **convRfft, **convIfft;

  beamR = (float **) malloc(nx*sizeof(float *));
  for (i=0;i<nx;i++) beamR[i] = (float *) malloc(ny*sizeof(float));
  beamI = (float **) malloc(nx*sizeof(float *));
  for (i=0;i<nx;i++) beamI[i] = (float *) malloc(ny*sizeof(float));

  beamRfft = (float **) malloc(nx*sizeof(float *));
  for (i=0;i<nx;i++) beamRfft[i] = (float *) malloc(ny*sizeof(float));
  beamIfft = (float **) malloc(nx*sizeof(float *));
  for (i=0;i<nx;i++) beamIfft[i] = (float *) malloc(ny*sizeof(float));

  convR = (float **) malloc(nx*sizeof(float *));
  for (i=0;i<nx;i++) convR[i] = (float *) malloc(ny*sizeof(float));
  convI = (float **) malloc(nx*sizeof(float *));
  for (i=0;i<nx;i++) convI[i] = (float *) malloc(ny*sizeof(float));

  convRfft = (float **) malloc(nx*sizeof(float *));
  for (i=0;i<nx;i++) convRfft[i] = (float *) malloc(ny*sizeof(float));
  convIfft = (float **) malloc(nx*sizeof(float *));
  for (i=0;i<nx;i++) convIfft[i] = (float *) malloc(ny*sizeof(float));

  mapR = (float **) malloc(nx*sizeof(float *));
  for (i=0;i<nx;i++) mapR[i] = (float *) malloc(ny*sizeof(float));
  mapI = (float **) malloc(nx*sizeof(float *));
  for (i=0;i<nx;i++) mapI[i] = (float *) malloc(ny*sizeof(float));

  mapRfft = (float **) malloc(nx*sizeof(float *));
  for (i=0;i<nx;i++) mapRfft[i] = (float *) malloc(ny*sizeof(float));
  mapIfft = (float **) malloc(nx*sizeof(float *));
  for (i=0;i<nx;i++) mapIfft[i] = (float *) malloc(ny*sizeof(float));

/* This is the center of the model grid */
  xcent = xcenter[nx/2] - 0.5*cellx;
  ycent = ycenter[ny/2] - 0.5*celly;

/* sigma squared in the x and y directions */
  s2x = beamx*beamx;
  s2y = beamy*beamy;

/* ratio of smoothed beam to input cell size */
  area = PI * beamx*beamy / (cellx*celly);

/* Build the 2D beam */
     for (i=0;i<nx;i++) {
     for (j=0;j<ny;j++) {
            dx = xcenter[i] - xcent;
            dy = ycenter[j] - ycent;
            dx = dx * dx;
            dy = dy * dy;
            beamR[i][j] = exp( -dx/s2x - dy/s2y );
            beamI[i][j] = 0.0;
     }}

/* FFT of the beam */
  inverse = 0;
  status = fft2d(inverse, nx, ny, beamR, beamI, beamRfft, beamIfft);

  printf("Beam FFT\n");

     for (k=0;k<nz;k++){

/* Copy the 3D cube into 2D array */
     for (i=0;i<nx;i++) {
     for (j=0;j<ny;j++) {
         mapR[i][j] = input_cube[i][j][k];
         mapI[i][j] = 0.0;
     }}

/* FFT the map */
     inverse = 0;
     status = fft2d(inverse, nx, ny, mapR, mapI, mapRfft, mapIfft);

/* Multiply the FFT beam and the FFT map to perform the convolution */
     for (i=0;i<nx;i++) {
     for (j=0;j<ny;j++) {
        convRfft[i][j] = mapRfft[i][j]*beamRfft[i][j];
        convIfft[i][j] = mapIfft[i][j]*beamRfft[i][j];
     }}

/* Inverse FFT of the product to real space */
     inverse = 1;
     status = fft2d(inverse, nx, ny, convR, convI, convRfft, convIfft);

/* Copy 2D map back to 3D cube */
     for (i=0;i<nx;i++) {
     for (j=0;j<ny;j++) {
        smoothed_cube[i][j][k] = convR[i][j]/area;
     }}

     }

     return 1;
}

int read_flash(char filename[])
{

unsigned long mfsize(FILE *);
 
FILE *fplynds;
int   iplynds;
unsigned long  file_size;

int i,j,k,l,m,n,nx,ny,nz,ngrid;
long nobj,nbytes,nbytes_remain;
int  nsets,nbytes_per_set;
int my_rank;

float cell,rs,vk;
float x,y,z;

float float_var;
double mindensity=1.e20,maxdensity=-1.e20,
        mint=1.e20,maxt=-1.e20,
	minx=1.e20,maxx=-1.e20,
	miny=1.e20,maxy=-1.e20,
	minz=1.e20,maxz=-1.e20,
	minvx=1.e20,maxvx=-1.e20,
	minvy=1.e20,maxvy=-1.e20,
	minvz=1.e20,maxvz=-1.e20;
long im=0,jm=0,km=0,nm=0;

/*
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
printf("P%d: starting to read data for cloud \n",my_rank);
*/

fplynds = fopen(filename,"r");
if (fplynds == NULL) {
	printf("Could not open data file %s for cloud. Exiting ...\n",
		filename);
	fprintf(fpnotes,"Could not open data file %s for cloud. Exiting ...\n",
		filename);
        return 0;
} else {
	printf("Opened data file %s for cloud\n",
		filename);
	fprintf(fpnotes,"Opened data file %s for cloud\n",
		filename);
}
iplynds = fileno(fplynds);

file_size = mfsize(fplynds);
printf("size of input file is %lu bytes\n",file_size);
printf("size of input file is %12.6e Mb\n",file_size/1.e6);


 n = fread(&nx,sizeof(int),1,fplynds);
 n = fread(&ny,sizeof(int),1,fplynds);
 n = fread(&nz,sizeof(int),1,fplynds);
 n = fread(&cell,sizeof(float),1,fplynds);

 printf("nx,ny,nz = %d, %d, %d\n",nx,ny,nz);
 printf("cell, = %e, \n",cell);
 printf("grid size, = %e, \n",cell * nx/2 / 3.085e18);
 fprintf(fpnotes,"nx,ny,nz = %d, %d, %d\n",nx,ny,nz);
 fprintf(fpnotes,"cell = %e, \n",cell);
 fprintf(fpnotes,"grid size, = %e, \n",cell * nx/2 / 3.085e18);

/* allocate memory for the model, the size of the model grid will be the first integer
   in the data file */

    avery = (struct avery***) malloc(nx*sizeof(struct avery**));
    for (i=0;i<nx;i++){
      avery[i] = (struct avery **) malloc(ny*sizeof(struct avery *));
      for (j=0;j<ny;j++){
        avery[i][j] = (struct avery *) malloc(nz*sizeof(struct avery));
        }}

/* Loop here over cube */

nbytes = 0;
ngrid = 0;

	for (k=0;k<nx;k++){
	for (j=0;j<ny;j++){
	for (i=0;i<nz;i++){

 		n = fread(&avery[i][j][k].vx,sizeof(float),1,fplynds);
 		n = fread(&avery[i][j][k].vy,sizeof(float),1,fplynds);
 		n = fread(&avery[i][j][k].vz,sizeof(float),1,fplynds);
 		n = fread(&avery[i][j][k].density,sizeof(float),1,fplynds);
 		n = fread(&avery[i][j][k].temperature,sizeof(float),1,fplynds);
 		n = fread(&avery[i][j][k].ionfraction,sizeof(float),1,fplynds);
 		n = fread(&avery[i][j][k].neutraldensity,sizeof(float),1,fplynds);
 		n = fread(&avery[i][j][k].moleculardensity,sizeof(float),1,fplynds);

if (i < 2 && j < 2 && k < 2) printf("ijk %3d %3d %3d vxyz density temp %12.3e %12.3e %12.3e %12.3e %12.3e\n",
i,j,k,avery[i][j][k].vx,avery[i][j][k].vy,avery[i][j][k].vz,avery[i][j][k].density,avery[i][j][k].temperature);

 		avery[i][j][k].density /= (2.33*MH);

		if (n != EOF) {
			ngrid += 1;
		} else {
			printf("reached EOF at %d %d %d\n",i,j,k);
			break;
		}

		if (avery[i][j][k].vx > maxvx) maxvx = avery[i][j][k].vx;
		if (avery[i][j][k].vy > maxvy) maxvy = avery[i][j][k].vy;
		if (avery[i][j][k].vz > maxvz) maxvz = avery[i][j][k].vz;
		if (avery[i][j][k].temperature > maxt) maxt = avery[i][j][k].temperature;
		if (avery[i][j][k].density > maxdensity) {
			im = i;
			jm = j;
			km = k;	
			nm = ngrid;
			maxdensity = avery[i][j][k].density;
		}
		if (avery[i][j][k].vx < minvx) minvx = avery[i][j][k].vx;
		if (avery[i][j][k].vy < minvy) minvy = avery[i][j][k].vy;
		if (avery[i][j][k].vz < minvz) minvz = avery[i][j][k].vz;
		if (avery[i][j][k].temperature < mint) mint = avery[i][j][k].temperature;
		if (avery[i][j][k].density < mindensity && 
			avery[i][j][k].density > 0.0) 
			mindensity = avery[i][j][k].density;
	}}}

	n = ngrid/(nx*ny*nz);
	if (n < 1) {
		printf("Did not finish reading through model grid\n");
		printf("Finished %d lines, which is %f percent of total\n",
			ngrid,(float)ngrid/(float)(nx*ny*nz)*100.);
		fprintf(fpnotes,"Did not finish reading through model grid\n");
		fprintf(fpnotes,"Finished %d lines, which is %f percent of total\n",
			ngrid,(float)ngrid/(float)(nx*ny*nz)*100.);
		return 0;
	} else {
		printf("Appears to have successfully read Ulrich data\n");
		fprintf(fpnotes,"Appears to have successfully read Ulrich data\n");
	}

	printf("range vx %f %f (kms)\n",minvx/100000.,maxvx/100000.);
	printf("range vy %f %f (kms)\n",minvy/100000.,maxvy/100000.);
	printf("range vz %f %f (kms)\n",minvz/100000.,maxvz/100000.);
	printf("range density %12.3e %12.3e\n",mindensity, maxdensity);
	printf("range tempera %12.3e %12.3e\n",mint, maxt);
	printf("density max i,j,k %ld %ld %ld \n",im,jm,km);
	fprintf(fpnotes,"range vx %f %f (kms)\n",minvx/100000.,maxvx/100000.);
	fprintf(fpnotes,"range vy %f %f (kms)\n",minvy/100000.,maxvy/100000.);
	fprintf(fpnotes,"range vz %f %f (kms)\n",minvz/100000.,maxvz/100000.);
	fprintf(fpnotes,"range density %g %g\n",mindensity, maxdensity);
	fprintf(fpnotes,"density max i,j,k %ld %ld %ld \n",im,jm,km);

	return 1;

} /* end of read_flash function */
