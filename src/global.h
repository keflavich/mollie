#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include "definitions.h"
#include "nlines_C.h"
 
#define PLANCK  6.6252e-27
#define BOLTZ   1.3806e-16
#define C       3.e10
#define AMU     1.660531e-24
#define PC      3.085e18
#define AU      1.5e13
#define GRAV    6.67e-8
#define PI      3.14159265359
#define SQRTPI  1.77245385091
#define MSUN	2.e33
#define YR	3.15360e7
#define RGAS	8.31421e7
 
 
#define MAXHYP 50	/* number of hyperfine lines per transition */
#define MAXNW  500	/* number of line widths stored */
#define MAXCH  2000	/* number of channels */
#define MAXBXS 6	/* number of model cubes */
#define MAXVW 50        /* maximum number of angles to view */
#define MAXVOX 1024	/* maximum pathlength through model cube */
#define NL 2            /* parameters for a plotting file. Not normally used */
#define NP 4            /* parameters for a plotting file. Not normally used */
#define MAXTRANS 500    /* number of transitions for non-LTE hyperfine */
#define MAXNAME 20      /* number of characters for linenames */

 
FILE *fpnotes;
FILE *fpnotes1;
FILE *fpconv;
FILE *fphist;
FILE *fpchkpt;
FILE *fpspectra;
int   ipnotes;
int   iphist;
int   ipconv;
int   ipchkpt;
int   ipspectra;

struct radex_data {
        int molecular_weight;
        float brot;
        float dipole;
        int number_levels;
        float *energies;
        float *gweights;
        int nqs;
        int **qnumbers;
        int number_transitions;
        int *upper;
        int *lower;
        double *frequency;
        double *einsteina;
        double *relative_intensity;
        int number_collisions;
        int number_temperatures;
        float *temperatures;
        int *jc_upper;
        int *jc_lower;
        float **jc_rates;
};
struct radex_data radex_data;
 
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

struct model ****model;

struct dcube {
        int nlines;
        int nviews;
        int *nchan;
        int nx;
        int ny;
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
 
struct dcube dcube;
 
/* In these structures, junk is there because there needs to
   be a variable after the [][]. That was the case for one
   of the compilers at one time. No explanation. 
*/
struct apLJbar {
        float apL[NLINES][MAXHYP];
        float Jbar[NLINES][MAXHYP];
	float norm[NLINES][MAXHYP];
	float junk;
};
 
struct apLJbar ****apLJbar;
struct apLJbar ****apLJbar0;
 
struct srcopc {
        float emissivity[NLINES][MAXHYP];
        float opacity[NLINES][MAXHYP];
	float cont_emissivity[NLINES];
	float cont_opacity[NLINES];
	float junk;
};

struct srcopc ****srcopc;

float dust_opacity[NLINES];
 
struct atomic {
        double freq[NLINES];
        double einsta[NLINES];
        double aulhyp[NLINES][MAXHYP];
        double statdg[NSTATE];
        int    upper[MAXTRANS];
        int    lower[MAXTRANS];
	int    jn[NSTATE];
	int    kn[NSTATE];
	int    ln[NSTATE];
	int    indexu[NLINES];
	int    indexl[NLINES];
	int    nhyp[NLINES];
        int    start[NLINES];
};
 
struct atomic atomic;

extern float setup_(int *, int *, int *, int *);
int stateq(int , int[], int[], int[], int, int, int, int, int);
 
int iterations ;
int ngrid,ngridx;
int nltehyp;


float chisq_function( );
int read_data();

int read_radex(char[]);

int define_model(int, int , int nx[MAXBXS],int ny[MAXBXS], int nz[MAXBXS],
                 float **,  float **,  float **, 
                 float, float, float, int, int, float, int);

int close_output_cube(int, char[] );

int replaceHi(int, int, float, float *, float**, float *, float psum[]);

int anneal(int , float, float **, float * , float , int , float);
float xcl(int * ,float );
float pcl(float ,float );
float ran1(int *);
void hunt(float *, int, float, int *);
int move(float , float, int *);
int newStep( int, int, float**, float[], float[], float *);
int reorder( int, int *, int *, int *, float *);
extern int sourceopac_( float *, double[], double[], float *, float[], float[], int *);
extern int acceleratedlambda_( float *, float *, 
                               float [NLINES][MAXHYP], float [NLINES][MAXHYP], 
                               float [NLINES][MAXHYP], 
			       float [NLINES][MAXHYP], float [NLINES][MAXHYP],
		               double [], int*, int*, int*, int*, int*);
int avg_output_cube(int);
extern int newrt_(int *, int *, int *, int *,
                  float *, float *, float *, 
                  float *, float *, float *, 
		  float *, float *, 
                  float *, float *, 
                  float *, float *, 
                  float[],float[],float[],float[]);
int convergence( int, int[], int[], int[], int, int *);
int chkpt( int, int[], int[], int[], int *, int *, int *, int *, int*);
extern int setcmb_( int *, float *);
extern int initrad_( int *, int[], int[], int[], int * ); 
int set_model_parameters(int , float *, float **, int, int );
extern int raypath_(float *, float *, float *,
		float *, float *, float *,
		int *, int[], int[], int[], float[], float[],
		int *, float *, float *, float *,
		int *, int*, float *, int *);
extern int rays_( int *, int *, float *, int *,
		int *, int *, int[], int[], int[],
		float *,float *,float *,float *,float *,float *,
		int *, float *, float *,
		float *,float *,float *,float *,float *,float *,
		float *,float *,float *,int *);
extern int channels_( float *, float *, float *, float *, int[], float[][MAXCH] );
extern int n2hphyp_( float *);
extern int n2dphyp_( float *);
extern int nh3hyp_( float *);
extern int hcnhyp_( float *);
extern int c17ohyp_( float *);
extern float copyatomic_( double[],double[],double[],int[],int[],int[],int[],
                          int[],int[],int[],int[],int[],double[][MAXHYP]);
extern int f77open_();
extern float crtsph_(float *,float *,float *, float *,float *,float *, int *);
extern void getmodel_(int *, int *, int *, int *,
                float *, float *, 
                float *, float *, float *, float *, float *,
                int *); 
extern void putsrcopc_(int *, int *, int *, int *, int *,
                        float *, float *);
extern void getalbj_(int *, int *, int *, int *, int *, int *,
                        float *, float *, float *);
extern void putalbj_(int *, int *, int *, int *, int *, int *,
                        float *, float *, float *);
int pacdata(int *, int *, int *, float *, float *, float *, float *);

extern int ltenh3_(float *, double[] );
extern int ltech3cn_(float *, double[] );
extern int subcritical_(int *, int *, int *, int *, float *, float *, float *, float *, float *, float *, float [NLINES][MAXHYP], double[] );
extern float fraction_para_(float *);

/*
extern int states2lines_(double[], double[], double[], float *, int *, int *, int *, int *);
*/

extern int depcoeffsn2hp_(int *, int *, int *, int *, float *, 
         double[], double[], double[], float[NSTATE]);

extern int locate_(float *, int *, float *, int *);

void fft(int N, double **x, double **X);
void fft_rec(int N, int offset, int delta,
             double **x, double **X, double **XX);
void ifft(int N, double **x, double **X);
int fft2d (int inverse, int nx, int ny, float **beamR, float **beamI,
           float **fftR, float **fftI);

unsigned long mfsize(FILE *);
