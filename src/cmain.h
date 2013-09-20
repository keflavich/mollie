int nbox;
int ibox,i,j,k,m,n,it;
int l,q,p,nh;
 
int nx[MAXBXS], ny[MAXBXS], nz[MAXBXS];
float cellsizex[MAXBXS],cellsizeyz[MAXBXS];
 
float ** xplane, ** yplane, ** zplane;
float ** xcenter, ** ycenter, ** zcenter;
 
float radius;
 
char fname[64],  suffix[4], pathname[64];
 
int idir;
float r,az,el;
float x,y,z;
float fq,qfu,qju,statwt;
float sinlong,sinlat,coslong,coslat;
int rw,status,master_status;
 
float anneal_t0,h20,h21,abundance,vr,turbwidth;
 
float chanwd,velrange,lwmin,lwmax;
 
float atoms;
 
int noutx,nouty;
float cellx,celly;
float beamx,beamy;
 
int molecule, dust_continuum, acceleration;
 
float *x1, *x2, *y1, *y2, *z1, *z2, *dxdu, *dydu, *dzdu;
float dxdu1,dydu1,dzdu1;
 
float *t1, *t2, *t3, *t4, *t5, *t6;
int *t7, *raybox;
float * t8,*t9;
int maxgrid;
int nP;
 
int nlng,nlat;
float *vlng, *vlat;
int nangles,outviews, *outAngleIndex, *outrpv;
float  outlng[MAXVW],  outlat[MAXVW];

int iray,nrays,rayc;
int *outray1, *rayview1, *rpview;
 
int nllist;
int linelist[NLINES];
int nchan[NLINES];
 
float chvel[NLINES][MAXCH];
 
int ivox,maxvox,ic;
int *icount;
int ** iarc;
float * pathlength; 
float * column_density_segment;
float column_density_sum;

float * arcl;
float * arclmin, * arclmax;
int * intersects;
 
int line,ihyp;
char **linenames;
 
int hanning;
 
float sumopac,cellopac,celltau,rymin,rr;
int minray;
 
int converged;
 
int iter;
int model_number=0,number_models,best_model=0;
int l1,l2;

float chisqmin=0.,chisqbest=1.e20;
float vc,vcmin,vcbest;

int accel;
 
void setupOutputCube(int, int*, int *, int, int, int,
  float, float, float, float, float [], float [], char **, float[NLINES][MAXCH]);
void setup_grid(
        int ,int *,int *,int *,
        float **,float **,float **,
        float **,float **,float **,
        float *,float*,
        int , int , int *, float *, float *,
        int, float [], float [], int *);
void  allocate_grid(
        int, int [], int [], int [],int ,
        int, int,
        float ***, float ***, float ***,
        float ***, float ***, float ***,
        int **, int **, float **, float **); 
void set_linenames(int , char **);

int simplex(int, float , float *);
int replaceHi(int , int, float , float *, float **, float *, float []); 

float **position_limits;
float  new_vertex_value;
float * new_vertex_position;

int nParameters;

int add_rays; 


int comparison;

int initial_radiation = 1;
int number_lines = 1;

float chksum;
int * sendcheck, * recvcheck;

int  check_model(int , int , int *, int *, int *,
                 float **,  float **,  float **, 
                 float , float) ;

int n_lines,n_state;

struct timeval starttime,splittime1,splittime2,splittime3,splittime4,stoptime,
	splittimeA,splittimeB,splittimeC;
struct timezone zone;

/* float * vxray, * vyray, * vzray, * pathray, * wdthray; */
/* float * dustemis, * dustopac; */
/* float emisray[MAXVOX][MAXHYP], opacray[MAXVOX][MAXHYP]; */
float apLray[MAXHYP][MAXVOX], barJray[MAXHYP][MAXVOX], normray[MAXHYP][MAXVOX];
float apL1d[MAXVOX*MAXHYP],barJ1d[MAXVOX*MAXHYP],yspec[MAXCH];
float  * ycmb;
float apLsum,barJsum,normsum,opacsum,emissum,srcsum;

int nsets, iset, istart, iend,  kview, iline, mline, iret, sline, lineCount;

float * comm_floats;
float * comm_floats4;
float * comm_floats5;
int * comm_ints; 
int * raylist;

int rays_to_compute, rays_finished, first_proc, last_proc, send_iray, recvd_iray, 
    current_ray, free_processes;

int depletion,photodissociation;
float std_abundance;
void chemistry(int depletion, int photodissociation, int molecule,
	int nbox, int nx[],int ny[], int nz[]);

int iterc, orig, lline;

float xst,yst,zst,xnd,ynd,znd;
