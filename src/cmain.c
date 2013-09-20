/* Version 24 Jan 2013

Corrections:

Mar 8, 2013
Cnotes0 includes the optical depth per cell printed after
every iteration along with source, Jbar, opacity, Lambda.
Corrected to print this data along the ray through the 
center of the cube. This changes cmain.h.

Jan 24, 2013
Corrected a problem in cmain.c which initialized
the molecular opacity and emissivity to zero
before the statistical equilibrium calculation.
This prevented convergence at very high optical depth.
The initialization to zero was meant to apply only
to cells defined as ionized gas.

Jan 11, 2013
This version modifies the ALI for consistency in
case the normalization of the line profile function is not
exactly one. The line profile function is required to be
correct to within 1%. The exact number is now propagated
from the slaves to the master and included in the rate
equations with ALI. This prevents drift of the solution.

May 4, 2011
Corrected gridding for the case where the
model sphere is much larger than the model box.
Renamed alng,alat,nviews to outlng, outlat, outviews.

May 2 2011
Output file now includes rest frequencies. Requires
updated version of load_model.

8 April 2011
Argument list to define_model.c now includes molecule.
read_lynds reads new hydro output files with H2O abundance.

18 March 2011
New chemistry function with improved CO and new H2O ice
model. Fixed an asymmetry in the calculation of Av.

3 March 2011 
Corrected the geometry to fix the problem that 
the number of angles in the C and F77 parts of 
the code were not the same when the view down 
the north pole was included. 

Be sure to turn off the chemistry by setting 
depletion = NO
photodissociation = NO
when using the abundances from the hydro output

*/
#include "definitions.h"
#include "global.h"
#include "mpi.h"
#include <unistd.h>
#include <limits.h>

int main(int argc, char *argv[]) {
#include "cmain.h"

    int         my_rank;       /* rank of process      */
    int         processes;     /* number of processes  */
    int         source;        /* rank of sender       */
    int         dest;          /* rank of receiver     */
    int         tag = 0;       /* tag for messages     */
    char        message[100];  /* storage for message  */
    MPI_Status  mpi_status;    /* return status for    */
                               /* receive              */
    MPI_Request *mpi_request;



    char*	my_host;	/* string for host name */
    int 	result;		/* return result	*/
    long	elapsed;	/* time elapsed 	*/

/* Done with the variable definitions */
/* SETUP.C includes executable statements */
#include "setup.c"

    printf("Starting mpi \n"); 
    fflush(stdout);

    my_host = (char *)malloc(32*sizeof(char));
    result = gethostname(my_host,32);
    gettimeofday(&starttime, &zone);


    /* Start up MPI */
    MPI_Init(&argc, &argv);

    /* Find out process rank  */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

/* Allocate the memory for error handling before there are any arrors */
    comm_ints    = (int *)malloc(5*sizeof(int));

    if (processes < 1) {
	printf("Run this program with at least 2 processes even if the computer has 1 CPU\n");
	printf("for example, mpiexec -n 2 ./c.x\n");
	printf("Going to exit ...\n");
        fflush(stdout);
/* On a detected fatal error the slaves go to their work loop and wait for instructions.
   The master goes to exitnow for orderly shutdown which includes calling the slave 
   processes to shut down. */
        if (my_rank != 0) goto slavecage; 
	goto exitnow;
    }


    printf("Starting process %d on %s \n",my_rank,my_host);
    fflush(stdout);

/* ------------------------------------------------------------------- */

/* Open some input and output files */

/* Slaves cannot write to Cnotes0. No fprintf(fpnotes ... if my_rank != 0 */

    if (my_rank == 0) {

    printf("size of int %d\n",(int)sizeof(int));
    printf("int max %d\n",INT_MAX);
    printf("size of long %d\n",(int)sizeof(long));
    printf("long max %ld\n",LONG_MAX);

    }
    if ('/' != pathname[strlen(pathname)-1]) strcat(pathname,"/");


/* Only process 0 and 1 will open an Fnotes file, but all the processes
   call this function to set the process ID in F77 */
	i = f77open_(&my_rank);
	if (i == 0) {
		printf("Function f77open returned with a failed status\n");
		printf("Going to exit now ...\n");
		fflush(stdout);
/* On a detected fatal error the slaves go to their work loop and wait for instructions.
   The master goes to exitnow for orderly shutdown which includes calling the slave 
   processes to shut down. */
                if (my_rank != 0) goto slavecage;
        	goto exitnow;
	}


    if (my_rank == 1) {
        strcpy(fname,pathname);
        strcat(fname,"Cnotes");
        sprintf(suffix,"%d",my_rank);
        strcat(fname,suffix);
        fpnotes1 = fopen(fname,"w");
             if (fpnotes1 == NULL) {
                printf("Could not open %s. Exiting ...\n",fname);
                goto exitnow;
             }
        printf("File %s opened \n",fname);
        fflush(stdout);
    }

    if (my_rank == 0)  {

        strcpy(fname,pathname);
        strcat(fname,"Cnotes");
	sprintf(suffix,"%d",my_rank);
	strcat(fname,suffix);
        fpnotes = fopen(fname,"w");
             if (fpnotes == NULL) {
                printf("Could not open %s. Exiting ...\n",fname);
                goto exitnow;
             } 
        printf("File %s opened \n",fname);
        ipnotes = fileno(fpnotes);
	fflush(stdout);
    
             strcpy(fname,pathname);
             strcat(fname,"Convergence");
	     strcat(fname,suffix);
             fpconv = fopen(fname,"w");
             if (fpconv == NULL) {
                printf("Could not open %s. Exiting ...\n",fname);
		fflush(stdout);
                goto exitnow;
             } 
             printf("File %s opened\n",fname);
             ipconv = fileno(fpconv);
             strcpy(fname,pathname);
             strcat(fname,"Search_history");
             fphist = fopen(fname,"w");
             if (fphist == NULL) {
                printf("Could not open %s. Exiting ...\n",fname);
		fflush(stdout);
                goto exitnow;
             } 
             printf("File %s opened\n",fname);
             iphist = fileno(fphist);

/*         printf("initial rw = %d \n",rw); */

/* If there is data to be compared read it here */
/* and open a file for plotting the comparison using my IDL program chx2.pro */

   if (comparison != NO) {
       result = read_data();
       if (result == 0) goto exitnow;
/*
The check_chisq file is now opened in define_model.c
       fpspectra = fopen("./check_chisq","w");
             if (fpspectra == NULL) {
                printf("Could not open %s. Exiting ...\n",fname);
		fflush(stdout);
                goto exitnow;
             } 
       printf("File %s opened\n",fname);
       ipspectra = fileno(fpspectra);
*/
   }

	fflush(stdout);
   } /* finished with master process i/o files */

/* Finished opening some input and output files */

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

/* Some checks for common setup errors */
/* NH3 is special in that it needs to have a specific number of lines and states */

/*
    fprintf(fpnotes,"size of int %d\n",(int)sizeof(int));
    fprintf(fpnotes,"int max %ld\n",INT_MAX);
    fprintf(fpnotes,"size of long %d\n",(int)sizeof(long));
    fprintf(fpnotes,"long max %ld\n",LONG_MAX);
*/

	j = 1;
	if (molecule == 3 && (NSTATE != 2*NLINES)) {
        if (my_rank < 2) printf("*** INPUT ERROR ***\n");
	if (my_rank < 2) printf("NH3 requires NSTATE = 2*NLINES \n");
	if (my_rank < 2) fflush(stdout);
	j = 0;
	}

	if (NSTATE > 10 ) {
	  if (molecule > 6 && molecule < 15) {
        if (my_rank < 2) printf("*** INPUT ERROR ***\n");
	if (my_rank < 2) printf("NSTATE > 10 is too many states for a rotor molecule\n");
	if (my_rank < 2) fflush(stdout);
	j = 0;
	}} 

	if (MAXVOX < NLINES ) {
        if (my_rank < 2) printf("*** INPUT ERROR ***\n");
	if (my_rank < 2) printf("MAXVOX must be > or = to NLINES\n");
	if (my_rank < 2) fflush(stdout);
	j = 0;
	}

	for (i=0;i<number_lines;i++) {
		if (linelist[i] > NLINES-1) {
	if (my_rank < 2) printf("*** INPUT ERROR ***\n");
	if (my_rank < 2) printf("The linelist in setup.c has too many lines\n");
	if (my_rank < 2) fflush(stdout);
        j = 0.;
	}}

        if (j == 0) { 
		printf("P%d failed initial input checks %d\n",my_rank,j);
		fflush(stdout);
/* On a detected fatal error the slaves go to their work loop and wait for instructions.
   The master goes to exitnow for orderly shutdown which includes calling the slave 
   processes to shut down. */
		if (my_rank != 0 )  goto slavecage;
		goto exitnow;
	}
	
/*        printf( "P%d finished initial input checks with good status  = %d\n",
        my_rank,j); 
	fflush(stdout);
*/

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

/* The next several blocks use the parameters for initializations including
   molecular constants, grids, memory ... */

/* ------------------------------------------------------------------- */
/* These calls set up the molecular constants. Have to initialize the
   pointer for linenames before passing it through the function set_linenames. */

/* ALL PROCESSES */

        nllist = number_lines;
        i = NLINES;
        j = NSTATE;

/* If we are computing non-LTE hyperfine emission, set the global variable,
   nltehyp, to true. Only N2H+ and HCN can be done with non-LTE hyperfines.
*/
	nltehyp = 0;
        if (
	    molecule == 4 ||
	    molecule == 18 ||
	    molecule == 19 ||
	    molecule == 21
	    ) nltehyp = 1;

/*	printf("Starting initial set up for molecule %d\n",molecule); */
	atoms = setup_(&molecule,&i,&j,&my_rank);
/*	printf("P%d: setup for molecule returned with status = %f non-LTE hyperfine: %d\n",
		my_rank,atoms,nltehyp);  
*/
	fflush(stdout);
	if (atoms < 1.) {
/* On a detected fatal error the slaves go to their work loop and wait for instructions.
   The master goes to exitnow for orderly shutdown which includes calling the slave 
   processes to shut down. */
		if (my_rank != 0) goto slavecage;
		printf("Received error condition from f77 function setup. Exiting\n");
		fprintf(fpnotes,"Received error condition from f77 function setup. Exiting\n");
		goto exitnow;
	}

/* Set up the LTE hyperfine structure if any. */

	i = 1;
	if (molecule == 9) {
		i = n2hphyp_(&lwmin);
	}
	if (molecule == 16) {
		i = n2dphyp_(&lwmin);
	}

        if (molecule == 13) {
                i = c17ohyp_(&lwmin);
        }
        if (molecule == 3 || molecule == 2) {
                nh3hyp_(&lwmin);
        }
        if (molecule == 15) {
           i = hcnhyp_(&lwmin);
        }
        if (i == 0) {
           printf("Look for error in Fnotes%d\n",my_rank);
/* On a detected fatal error the slaves go to their work loop and wait for instructions.
   The master goes to exitnow for orderly shutdown which includes calling the slave 
   processes to shut down. */
           if (my_rank != 0) goto slavecage;
           goto exitnow;
        }

	copyatomic_(atomic.freq,atomic.einsta,atomic.statdg,
			atomic.jn,atomic.kn,atomic.ln,
			atomic.indexu,atomic.indexl,atomic.nhyp,atomic.start,
			atomic.upper,atomic.lower,atomic.aulhyp);
/* printf("P%d: finished copyatomic\n",my_rank); */

if (my_rank == 0) {
        linenames = (char **)malloc(NLINES*sizeof(char*));
        set_linenames(molecule,linenames);
}
/* printf("P%d: finished linenames\n",my_rank); */

/*
 printf("P%d: finished hyperfine structure\n",my_rank); 
 fflush(stdout);
*/


/* This call defines the velocity grid. All processes need the velocity grid  */

/* The minimum and maximum linewidths defined in the model are required
   for the velocity gridding. The hyperfine structure is also required.
   Therefore the model and the hyperfine functions must be called before
   the velocity grid. 
*/

/* Channnels sets up the velocity grid. Returns a float. Return of 0 is bad, 1 is good */

	i = channels_(&chanwd,&velrange,&lwmin,&lwmax,nchan,chvel);
	if (i < 1) {
	printf("Process %d received stop signal %d from function channels. \n", my_rank,i);
	if (my_rank < 2) printf("Look for error in Fnotes%d\n",my_rank);
/* On a detected fatal error the slaves go to their work loop and wait for instructions.
   The master goes to exitnow for orderly shutdown which includes calling the slave 
   processes to shut down. */
	if (my_rank != 0) goto slavecage; 
	goto exitnow;
	} 
/*
	printf("P%d: Velocity grid defined\n",my_rank); 
	printf("nchan %d %d %d %d \n",nchan[0],nchan[1],nchan[2],
			nchan[3]); 
        fflush(stdout);
*/

/* ------------------------------------------------------------------- */

	ycmb = (float *)malloc(NLINES*sizeof(float));

        status = setcmb_(&initial_radiation,ycmb);

	if (my_rank == 0) {
	fprintf(fpnotes,"The CMB brightness at each transition freq\n");
		for (i=0;i<NLINES;i++) 
                fprintf(fpnotes,"Freq (GHz) CMB (K) %e %e\n",
		atomic.freq[i]/1.e9,ycmb[i]);
	}


slavecage:

/* Only process id 1 can write to Cnotes1 with file pointer fpnotes1. 
   No slaves can write to fpnotes. */

/*	printf (        "Slave P%d: allocating memory\n",my_rank); */
/* Allocate enough memory for comm_floats. Make it the larger of the
   sending and receiving data */
	l = 2*MAXVOX*MAXHYP + MAXCH;
	q = 3 + (5+2*MAXHYP)*MAXVOX;
	if (l > q) q = l;
        if (q > INT_MAX && my_rank==0) {
	printf("Model size requires %d variables for message passing\n",q);
	printf("The maximum size of an INT is %d\n",INT_MAX);
	printf("You need to reduce the sizes of \n");
        printf("MAXVOX = %d\n",MAXVOX);
        printf("MAXHYP = %d\n",MAXHYP); 
	goto exitnow;
        }
	comm_floats = (float *)malloc(q*sizeof(float));
	comm_floats4 = (float *)malloc(q*sizeof(float));
	comm_floats5 = (float *)malloc(q*sizeof(float));

/* 	printf ("P%d: allocated comm memory %ld\n",my_rank,k); */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

if (my_rank != 0) {
/* slaves only */


/*	printf (        "Slave P%d: waiting instructions\n",my_rank); */
        fflush(stdout);
	while (1) {

/* The master has solved the stat. eq. equations to calculate the 
   source, opacity, and departure coefficients. Here the slave
   receives this data for one line along one ray. The slave also
   gets the velocity and width info along the ray. The slave will
   solve the RT equation along the ray. */

/* The slaves wait here on this MPI_Recv. This call blocks a slave's further progress 
   until the message is received and the call is satisfied. */

        orig = 0;	/* coming from the master */
	n = 5;
	tag = 0;
        MPI_Recv(comm_ints,n,MPI_LONG,orig,tag,MPI_COMM_WORLD,&mpi_status); 

        if (comm_ints[3] == 0) {
/*		printf (        "Slave P%d: Received final call to exit\n",my_rank); */
		goto exitnow; 
	}

/* still alive?, toil for the master */

/* These were the assignments as sent from the master
*	comm_ints[0] 	 line;
*	comm_ints[1] 	 icount[send_iray];
*	comm_ints[2] 	 send_iray (or rayc here);
*	comm_ints[3] 	 status;
*	comm_ints[4] 	 iterc; iteration
*/


		sline = comm_ints[0];
		ic   = comm_ints[1]; 
		rayc = comm_ints[2];
		iterc= comm_ints[4];

/* In the slave block, we use sline for line so that line is not inadvertently
   changed. */

/* ic is now the number of data along each ray, the same as icount in the master process. 
   Use this to unpack the data. 
*/

/*if (my_rank == 1) fprintf(fpnotes1,"slave P%d got ray %d for line %d\n",
my_rank,rayc,sline); */

	nh = 1; if (nltehyp == 1) nh = atomic.nhyp[sline];
        q = 3 + (7+2*nh)*ic;
//	printf("Slave %d Receiving q = %d line = %d\n",my_rank,q,line); 
        tag = 1;
        n = 3 + 3*ic;
//	printf("Slave %d Receiving q = %d line = %d tag = % d\n",my_rank,n,line,tag); 
        MPI_Recv(&comm_floats[0],n,MPI_FLOAT,orig,tag,MPI_COMM_WORLD,&mpi_status);
        tag = 2;
        n = 2*ic;
//	printf("Slave %d Receiving q = %d line = %d tag = % d\n",my_rank,n,line,tag); 
        MPI_Recv(&comm_floats[3+3*ic],n,MPI_FLOAT,orig,tag,MPI_COMM_WORLD,&mpi_status);
        tag = 3;
        n = 2*ic;
//	printf("Slave %d Receiving q = %d line = %d tag = % d\n",my_rank,n,line,tag); 
        MPI_Recv(&comm_floats[3+5*ic],n,MPI_FLOAT,orig,tag,MPI_COMM_WORLD,&mpi_status);
        tag = 4;
        n = nh*ic;
//	printf("Slave %d Receiving q = %d line = %d tag = % d\n",my_rank,n,line,tag); 
        MPI_Recv(&comm_floats4[0],n,MPI_FLOAT,orig,tag,MPI_COMM_WORLD,&mpi_status);
//        MPI_Recv(&comm_floats[3+7*ic],n,MPI_FLOAT,orig,tag,MPI_COMM_WORLD,&mpi_status);
        tag = 5;
        n = nh*ic;
//	printf("Slave %d Receiving q = %d line = %d tag = % d\n",my_rank,n,line,tag); 
        MPI_Recv(&comm_floats5[0],n,MPI_FLOAT,orig,tag,MPI_COMM_WORLD,&mpi_status);
//        MPI_Recv(&comm_floats[3+(7+nh)*ic],n,MPI_FLOAT,orig,tag,MPI_COMM_WORLD,&mpi_status);

/* augment the LINE index from C to F77. 
   ICOUNT is a total not an index, so do not increment
   We count the iteration, iterc, from 1 in both our C and F77 programs.
   The ray number, rayc, is not used in F77 except for debugging.
   OK to leave rayc in C indexing. */

/* printf("Starting newrt for line %d\n",sline); */

/* Increment the line by 1 before sending to newrt */
		j = sline + 1;

                status = newrt_(&rayc,		/* iray  ray number */
			&j,			/* line  */
			&ic,			/* ic number of cells */
			&iterc,			/* iteration number */
                        &comm_floats[0], 	/* dxdu */
			&comm_floats[1], 	/* dydu */
			&comm_floats[2],	/* dzdu */
                        &comm_floats[3+0*ic],	/* vxray  velocity for cell */
			&comm_floats[3+1*ic],	/* vyray */
			&comm_floats[3+2*ic],	/* vzray */
                        &comm_floats[3+3*ic],	/* wdthray line width for cell */
			&comm_floats[3+4*ic], 	/* pathray length of each segment */
			&comm_floats[3+5*ic],	/* cont_emis continuum emission */
                        &comm_floats[3+6*ic],	/* cont_opac continuum opacity */
			&comm_floats4[0],	/* emis1d line emissivity */
                        &comm_floats5[0],  	/* opac1d line opacity */
                        &comm_floats[0],            /* start of apl1d  */
                        &comm_floats[ic*nh],        /* start of barj1d */
                        &comm_floats[2*ic*nh],      /* start of norm1d */
                        &comm_floats[3*ic*nh]);     /* start of yspec  */

/* The last 4 comm_floats in the arg list are the output apL1d,barJ1d,xnorm1d,yspec. 
   These occupy the same memory space in comm_floats as the inputs. The output
   variables are filled in at the end of NEWRT after the inputs are no longer
   needed.  These output variables are
   sent to the master in the second MPI_Send below. The first MPI_Send sends back the
   ray number and pixel count, followed by some status info and the iteration number.
*/

/*
for (q=0;q<ic*nh;q++){
	printf("rank q bj al %d %d %f %f\n",my_rank,q,comm_floats[ic*nh+q],comm_floats[q]);
}
*/

/* printf("Finished newrt for line %d\n",sline); */

/* There are 3 variables output by newrt to send to the master, 
   apLray, barJray and yspec.
   ic is the number of voxels along the ray. The total number of floats 
   is 2*ICOUNT*NHYP[LINE] + NCHAN.

*/


        comm_ints[3] = 1;			/* This is the status flag */
        if (status == 0) comm_ints[3] = 0;  	/* This is an error condition from newrt */
	comm_ints[4] = iterc;			/* Same as received */
	
        n = 5;
	tag = 0; 				/* 0 for comm_ints, 1 for comm_floats */
	MPI_Send(comm_ints,  n,MPI_LONG, 0, tag,MPI_COMM_WORLD);

/* Send the array of comm_floats */

/* Make Q the total number of floats */
/*	if (q != 2*ic*nh + nchan[sline]) printf("Help, slave cant count\n"); */
	q = 3*ic*nh + nchan[sline];
        tag = 1;

/*	printf("Slave %d sending %d vars to master with final %f\n",
                my_rank,q,comm_floats[q-1]);
*/
        MPI_Send(comm_floats,q,MPI_FLOAT,0,tag,MPI_COMM_WORLD);

/*if (my_rank == 1) fprintf(fpnotes1,"slave P%d sent ray %d for line %d\n",
my_rank,rayc,sline);*/

	} /* slaves remain in this endless loop. can get no further
		unless they are sent to die  */

} /* end of slave if block */
/* ---------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */


/* Everything below is MASTER ONLY */
/* except the exitnow block */

/* More consistency checks */
          i = NLINES;
          j = NSTATE;
          if (i != n_lines || j != n_state) {
            printf("The number of lines and levels must agree in all 3 files,");
            printf("setup.c, nlines_C.h and nlines_f77.h\n");
            printf("P%d: In setup.c lines = %d and levels = %d\n",
                      my_rank,n_lines,n_state);
            printf("P%d: In nlines_C.h lines = %d and levels = %d\n",
                      my_rank,i,j);
            goto exitnow;
          }

          fprintf (fpnotes,"The number of atoms in this molecule is %f\n",
             atoms);
          fprintf(fpnotes,"Transition      Frequency       Einstein A statwt(upper) statwt(lower)\n");
          printf("Transition      Frequency       Einstein A statwt(upper) statwt(lower)\n");
          for (line=0;line<NLINES;line++){
              printf("%s %16.9f %12.6g %8.0f %8.0f\n",linenames[line],
                  atomic.freq[line]/1.e9,atomic.einsta[line],
		  atomic.statdg[atomic.indexu[line]],atomic.statdg[atomic.indexl[line]]);
              fprintf(fpnotes,"%s %16.9f %12.6g %8.0f %8.0f\n",linenames[line],
              atomic.freq[line]/1.e9,atomic.einsta[line],
		  atomic.statdg[atomic.indexu[line]],atomic.statdg[atomic.indexl[line]]);
//printf("line indexl indexu %d %d %d\n",line,atomic.indexu[line],atomic.indexl[line]);
	  }

	gettimeofday(&splittime1, &zone);
	elapsed = (splittime1.tv_sec  - starttime.tv_sec ) * 1000
                + (splittime1.tv_usec - starttime.tv_usec) / 1000;
	fprintf(fpnotes,"P%d: Finished molecule setup at %ld msecs\n",
    		my_rank,elapsed);


/* ------------------------------------------------------------------- */

printf("start allocate grid\n");


/* Allocate memory for the grid */
/* If there are n cells, there will be n+1 planes defining the cells */

        xplane = (float **) malloc(nbox*sizeof(float*));
        yplane = (float **) malloc(nbox*sizeof(float*));
        zplane = (float **) malloc(nbox*sizeof(float*));
        xcenter = (float **) malloc(nbox*sizeof(float*));
        ycenter = (float **) malloc(nbox*sizeof(float*));
        zcenter = (float **) malloc(nbox*sizeof(float*));

  for (ibox=0;ibox<nbox;ibox++){

        xplane[ibox] = (float *) malloc((nx[ibox]+1)*sizeof(float));
        yplane[ibox] = (float *) malloc((ny[ibox]+1)*sizeof(float));
        zplane[ibox] = (float *) malloc((nz[ibox]+1)*sizeof(float));
        xcenter[ibox] = (float *) malloc(nx[ibox]*sizeof(float));
        ycenter[ibox] = (float *) malloc(ny[ibox]*sizeof(float));
        zcenter[ibox] = (float *) malloc(nz[ibox]*sizeof(float));
  }

	vlng = (float *) malloc(nlng*nlat*sizeof(float));
	vlat = (float *) malloc(nlng*nlat*sizeof(float));
	outAngleIndex = (int *) malloc(nlng*nlat*sizeof(int));

        printf("Memory allocation for grid completed\n");
        fprintf(fpnotes,"Memory allocation for grid completed\n");

printf("start grid setup\n");

  setup_grid( 
	nbox,nx,ny,nz,
	xplane,yplane,zplane,
	xcenter,ycenter,zcenter,
	cellsizex,cellsizeyz,
	nlng, nlat, &nangles, vlng, vlat,
	outviews, outlng, outlat, outAngleIndex
  );

printf("finished grid\n");

/* allocate memory for the model */

        model = (struct model ****) malloc(nbox*sizeof(struct model ***));
        for (ibox=0;ibox<nbox;ibox++){
          model[ibox] = (struct model ***) malloc(nx[ibox]*sizeof(struct model **));
          for (i=0;i<nx[ibox];i++){
            model[ibox][i] = (struct model **) malloc(ny[ibox]*sizeof(struct model *));
            for (j=0;j<ny[ibox];j++){
              model[ibox][i][j] = (struct model *) malloc(nz[ibox]*sizeof(struct model));
        }}}

printf("finished allocating memory for the model cube\n");

/* Zero the ray counter and the check counter */
	for (ibox=0;ibox<nbox;ibox++){
	for (i=0;i<nx[ibox];i++){
	for (j=0;j<ny[ibox];j++){
	for (k=0;k<nz[ibox];k++){
            model[ibox][i][j][k].adds = 0.;
        }}}}




/* ------------------------------------------------------------------- */

/* Check the output grid */

/* 
   Calling setupOutputCube here prints out cell centers and planes so that
   the output grid can be checked without waiting for the radiative transfer
   calculations to finish. 

   The velocity grid is used in the output data cubes, so the velocity
   grid must be defined before the output data grid.

   Make a list of the output angles (outlng and outlat) for the output header.
   These differ from vlng and vlat in that the lists outlng and outlat contain
   only the angles of the views that were selected for output. vlng, vlat
   contain the angles of all the views that are computed.
*/
        setupOutputCube(nllist,linelist,nchan,noutx,nouty,outviews,
          beamx,beamy,cellx,celly,outlng,outlat,linenames,chvel);

printf("output cube defined\n");

/* ------------------------------------------------------------------- */

/* Calculate the beginning, end, and direction cosines of each ray. 
   Beginnings and ends are stored in x1,y1,z1 and x2,y2,z2, direction 
   cosines in dxdu,ydu,dzdu */

/* The maximum possible number of rays is the number of angles times
   the number of y,z cells in each box. Allocate this much memory, but the
   actual number of rays will be less because some rays will not intersect
   the model sphere. The number of rays will be calculated in the function
   rays. */


	nrays = 0;
        for (ibox=0;ibox<nbox;ibox++) nrays += nlng*nlat*ny[ibox]*nz[ibox];
	fprintf(fpnotes,"The maximum number of rays is %d\n",nrays);
	printf("The maximum number of rays is %d\n",nrays);

	x1 = (float *) malloc(nrays*sizeof(float));
	x2 = (float *) malloc(nrays*sizeof(float));
	y1 = (float *) malloc(nrays*sizeof(float));
	y2 = (float *) malloc(nrays*sizeof(float));
	z1 = (float *) malloc(nrays*sizeof(float));
	z2 = (float *) malloc(nrays*sizeof(float));
	dxdu = (float *) malloc(nrays*sizeof(float));
	dydu = (float *) malloc(nrays*sizeof(float));
	dzdu = (float *) malloc(nrays*sizeof(float));
        raybox = (int *) malloc(nrays*sizeof(int));
        raylist = (int *) malloc(nrays*sizeof(int));
        sendcheck = (int *) malloc(nrays*sizeof(int));
        recvcheck = (int *) malloc(nrays*sizeof(int));

        rpview = (int *) malloc(nangles*sizeof(int));
        rayview1 = (int *) malloc(nangles*sizeof(int));

/* ------------------------------------------------------------------- */


/* maxgrid and all the t vectors are temporary variables used to
   pass the arrays xcenter...,xplane..., to the function rays. 
   Whereas the x,y,z center and plane variables have variable
   dimensions for each nbox, the t variables have the same 
   dimension maxgrid. This is required for the F77 functions. */

	maxgrid =0;
	for (ibox=0;ibox<nbox;ibox++){
		if (nx[ibox] > maxgrid) maxgrid = nx[ibox];
		if (ny[ibox] > maxgrid) maxgrid = ny[ibox];
		if (nz[ibox] > maxgrid) maxgrid = nz[ibox];
	}

	t1 = (float *) malloc((maxgrid)*nbox*sizeof(float));
	t2 = (float *) malloc((maxgrid)*nbox*sizeof(float));
	t3 = (float *) malloc((maxgrid)*nbox*sizeof(float));
	t4 = (float *) malloc((maxgrid+1)*nbox*sizeof(float));
	t5 = (float *) malloc((maxgrid+1)*nbox*sizeof(float));
	t6 = (float *) malloc((maxgrid+1)*nbox*sizeof(float));

	for (ibox=0;ibox<nbox;ibox++){
	for (i=0;i<maxgrid;i++){
		k = maxgrid*ibox + i;
		if (i < nx[ibox]) t1[k] = xcenter[ibox][i];
		if (i < ny[ibox]) t2[k] = ycenter[ibox][i];
		if (i < nz[ibox]) t3[k] = zcenter[ibox][i];
	}}

	fprintf(fpnotes,"The maximum number of cells across any box is %d\n",maxgrid);
	printf("The maximum number of cells across any box is %d\n",maxgrid);

	for (ibox=0;ibox<nbox;ibox++){
	for (i=0;i<maxgrid+1;i++){
		k = (maxgrid+1)*ibox + i;
		if (i < nx[ibox]+1) t4[k] = xplane[ibox][i];
		if (i < ny[ibox]+1) t5[k] = yplane[ibox][i];
		if (i < nz[ibox]+1) t6[k] = zplane[ibox][i];
	}}
	

/*
	for (i=0;i<nlng;i++){
	for (j=0;j<nlat;j++){
	fprintf(fpnotes,"%d %d %f %f\n",i,j,vlng[i],vlat[j]);
	}}
*/

	    gettimeofday(&splittime1, &zone);
	    elapsed = (splittime1.tv_sec  - starttime.tv_sec ) * 1000
                    + (splittime1.tv_usec - starttime.tv_usec) / 1000;
	    fprintf(fpnotes,"P%d: Finished grids at %ld msecs\n",
    		    my_rank,elapsed);
	    printf (        "P%d: Finished grids at %ld msecs\n",
    		    my_rank,elapsed); 


        result = rays_(&nrays,&maxgrid,&radius,
               &nbox,rayview1,rpview,nx,ny,nz,
               t1,t2,t3,t4,t5,t6,
               &nangles,vlng,vlat,
               x1,y1,z1,x2,y2,z2,
               dxdu,dydu,dzdu,raybox);
	    printf ("P%d: Finished rays_ \n", my_rank); 
        if (result == 0) {
		printf("Received bad status %d from function rays_ \n",status);
		fprintf(fpnotes,"Received bad status %d from function rays_ \n",status);
		goto exitnow;
	}

/* subtract 1 from raybox to get from fortran to C */
	for (i=0;i<nrays;i++) raybox[i] -= 1;

	for (i=0;i<nangles;i++) {
	    printf("The number of rays in direction %d %f %f is %d\n",i,vlng[i],vlat[i],rpview[i]);
	    fprintf(fpnotes,"The number of rays in direction %d %f %f is %d\n",i,vlng[i],vlat[i],rpview[i]);
        }

	    printf("The total number of rays intersecting the model sphere is %d\n",nrays);
	    fprintf(fpnotes,"The total number of rays intersecting the model sphere is %d\n",nrays);


/* For the first view, there is one ray on each grid cell. YZ is the plane
   of the sky, and X is the line of sight. So for the first view, the rays
   run from some -X to +X, points where the ray intersects the sphere set 
   by radius. */




/*
	    fprintf(fpnotes,"Starting and ending points of each ray\n");
	    for (iray=0;iray<nrays;iray++)
		fprintf(fpnotes,"%5d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %d\n",
			iray,
			x1[iray], y1[iray], z1[iray], 
			x2[iray], y2[iray], z2[iray], raybox[iray]);
*/

/* Find the ray that is closest to the center of the cube. Use the first view.
   This is used later to print out the source, Jbar, and opacities after
   every iteration.
*/
        rymin = 1.e20;
        i = 0;
        for (iray=0;iray<rpview[0]-1;iray++) {
                rr = y1[iray]*y1[iray] + z1[iray]*z1[iray];
                if ( rr < rymin ){
                        i++;
                        rymin = rr;
                        minray = iray;
                }
        }
/*
        fprintf(fpnotes,"Searched %d rays. Found min at %d\n",rpview[0],minray);
        fprintf(fpnotes,"y1,z1,rmin %12.3e %12.3e %12.3e\n",y1[minray],z1[minray],rymin);
        printf("Searched %d %d rays. Found min at %d\n",i,rpview[0],minray);
        printf("y1,z1,rmin %12.3e %12.3e %12.3e\n",y1[minray],z1[minray],rymin);
*/






/* Count how many rays pass through each pixel. Store the number in
   model.adds The greatest number of pixels that a ray can cross is
   along a diagonal across the cube or twice the number along one side
   of the cube. Set the dimensions for the F77 arrays at 3*maximum grid
   * the number of boxes, and allocate this memory.*/

	gettimeofday(&splittime1, &zone);

	icount = (int *) malloc(nrays*sizeof(int));

	pathlength = (float *) malloc(MAXVOX*sizeof(float));
	column_density_segment = (float *) malloc(MAXVOX*sizeof(float));

		iarc = (int **) malloc(4*sizeof(int *));
		for (j=0;j<4;j++) {
			iarc[j] = (int *) malloc(MAXVOX*sizeof(int));
                }

	t7 = (int *)malloc(4*MAXVOX*sizeof(int));
	t8 = (float *)malloc(MAXVOX*sizeof(float));

/* add_rays counts up by twos because there is a forward and backward ray */
        add_rays = 2;
        maxvox = MAXVOX;

/* This computation of the ray paths is to load up add_rays into the 
   model cube. add_rays is the number of rays that cross each pixel. 
*/ 
  
	printf("Start computing ray paths\n"); 

	for (iray=0;iray<nrays;iray++){
		result = raypath_(&x1[iray],&y1[iray],&z1[iray],
		         &x2[iray],&y2[iray],&z2[iray],
                         &nbox,nx,ny,nz,cellsizex,cellsizeyz,
                         &maxgrid,t4,t5,t6,
                         &maxvox,t7,t8,&i);
                if (result == 0) {
			printf("Got zero as a result from raypath\n");
			goto exitnow;
		}
		icount[iray] = i-1; /* i is the last arg of raypath */	
		if (icount[iray] == 0) {
			fprintf(fpnotes,"1 Bad icount from raypath %d %d\n",icount[iray],iray);
			printf("1 Bad icount from raypath %d %d\n",icount[iray],iray);
		}

		for (ivox=0;ivox<icount[iray];ivox++){
			pathlength[ivox] = t8[ivox];
		for (i=0;i<4;i++){
			k = 4*ivox + i;
			iarc[i][ivox] = t7[k]-1; 
		}}
/* fprintf(fpnotes,"New ray # %d with %d segments\n",iray,icount[iray]); */
		for (ivox=0;ivox<icount[iray];ivox++){
			i = iarc[0][ivox];
			j = iarc[1][ivox];
			k = iarc[2][ivox];
			ibox = iarc[3][ivox];
			model[ibox][i][j][k].adds += add_rays;
/* fprintf(fpnotes,"ijk box %3d %3d %3d %3d adds %d\n",i,j,k,ibox,model[ibox][i][j][k].adds); */
		}

	} /* end of loop over rays */

	printf("Finished computing ray paths\n"); 

/* If there are nested boxes, then where the boxes overlap, the model cells 
   in the outer box will have adds = 0 because there will be cells in the
   smaller box that fill in this portion of the model. 

   If the nested boxes are concentric on the center, then there will be some
   cells in the outer region of the larger boxes that have a number of rays
   equal to 2 x the number of smaller cells that are seen through the larger
   cell.
*/

/*
		fprintf(fpnotes,"box    i      x       adds\n");
                for (ibox=0;ibox<nbox;ibox++)	{
		j = ny[ibox]/2 ;
		k = nz[ibox]/2 ;
	fprintf(fpnotes,"number of rays through each cell along a line j,k = %d %d\n", j,k);
                for (i=0;i<nx[ibox];i++) 	{
			fprintf(fpnotes,"%3d %4d  %9.6f %6d \n",
                        ibox,i,xcenter[ibox][i], model[ibox][i][j][k].adds);
		}}
*/
	printf("All rays defined\n"); 

/*
        for (i=0;i<nx[ibox];i++){
fprintf(fpnotes,"New layer %d\n",i);
        for (k=0;k<nz[ibox];k++){
        for (j=0;j<ny[ibox];j++){
		l = model[ibox][i][j][k].adds/2;
		if (l > 9) l = l;
		fprintf(fpnotes,"%2d",l);
	} fprintf(fpnotes,"\n");
	}}
*/

/* outray1 is a vector of length NVIEWS, containing the first ray 
   for each of the views to be saved in the output.  outlng and outlat
   are a list of the output angles. 
   rayview1 is the list of first rays of each computed view 
   rayview1 comes from the fsub.f function rays_ and counts in f77
   beginning with 1. Subtract off 1 to count from 0 in C. */
        outray1    = (int   *) malloc(outviews*sizeof(int  ));
        outrpv     = (int   *) malloc(outviews*sizeof(int  ));
	for (i=0;i<outviews;i++){
		fprintf(fpnotes,"i oAI %d %d\n",i,outAngleIndex[i]);
		outray1[i] = rayview1[outAngleIndex[i]] - 1;
		outlng[i] = vlng[outAngleIndex[i]];
		outlat[i] = vlat[outAngleIndex[i]];
		outrpv[i] = rpview[outAngleIndex[i]];
		fprintf(fpnotes,"view %3d:  first and last rays %5d %5d\n",
			i,outray1[i],outray1[i]+outrpv[i]-1);
		j = outray1[i] + outrpv[i]/2 + nx[raybox[outray1[i]]]/2;
		fprintf(fpnotes,"middle ray = %4d\n",j);
		fprintf(fpnotes,"x1 y1 z1 : x2 y2 z2  %7.4f %7.4f %7.4f : %7.4f %7.4f %7.4f\n",
			x1[j],y1[j],z1[j],x2[j],y2[j],z2[j]);
        }
	gettimeofday(&splittime2, &zone);

/* Finished defining all the rays through the model. */

	elapsed = (splittime2.tv_sec  - splittime1.tv_sec ) * 1000
                + (splittime2.tv_usec - splittime1.tv_usec) / 1000;
	fprintf(fpnotes,"P%d: Ray tracing requires %ld msecs\n",
    		my_rank,elapsed);
	printf (        "P%d: Ray tracing requires %ld msecs\n",
    		my_rank,elapsed);

/* ------------------------------------------------------------------------- */
/* Allocate memory for the approximate lambda operator, average intensity,
   source, and opacity */

        apLJbar = (struct apLJbar ****) malloc(nbox*sizeof(struct apLJbar ***));
        for (ibox=0;ibox<nbox;ibox++){
          apLJbar[ibox] = 
            (struct apLJbar ***) malloc(nx[ibox]*sizeof(struct apLJbar **));
          for (i=0;i<nx[ibox];i++){
            apLJbar[ibox][i] = 
              (struct apLJbar **) malloc(ny[ibox]*sizeof(struct apLJbar *));
            for (j=0;j<ny[ibox];j++){
              apLJbar[ibox][i][j] = 
                 (struct apLJbar *) malloc(nz[ibox]*sizeof(struct apLJbar));
        }}}


        srcopc = (struct srcopc ****) malloc(nbox*sizeof(struct srcopc ***));
        for (ibox=0;ibox<nbox;ibox++){
          srcopc[ibox] = 
            (struct srcopc ***) malloc(nx[ibox]*sizeof(struct srcopc **));
          for (i=0;i<nx[ibox];i++){
            srcopc[ibox][i] = 
              (struct srcopc **) malloc(ny[ibox]*sizeof(struct srcopc *));
            for (j=0;j<ny[ibox];j++){
              srcopc[ibox][i][j] = 
                 (struct srcopc *) malloc(nz[ibox]*sizeof(struct srcopc));
        }}}


/* ------------------------------------------------------------------------- */


/* Finished with all the initializations that do not depend on
   the variable model parameters. */
/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */

/* simulation for new model starts here */


/* ------------------------------------------------------------------------- */

/* Allocate some memory for the simplex vertices even if they are not
 * used 
 * There may be some number of parameters used for models even when
 * there is no comparison, but if the number of parameters has not
 * been set, then set it to 1 */

      if (comparison == NO && nParameters == 0) nParameters = 1;

      new_vertex_position = (float *)malloc(nParameters*sizeof(float));

      position_limits = (float **)malloc(nParameters*sizeof(float*));
      for (i=0;i<nParameters+1;i++){
            position_limits[i] = (float *)malloc(2*sizeof(float));
      }

      for (i=0;i<nParameters;i++){
          new_vertex_position[i] = 0.;
          position_limits[i][0] = 0.;
          position_limits[i][1] = 0.;
      }


/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
/* Loop over some number of models, changing parameters. If there is no
   fitting, then the loop below is run once. */

/* This loop runs through a number of models to search for
   a best fit */

  for (model_number=0;model_number<number_models;model_number++){


/* chisqmin is initialized to 0.0 in header, then set to
 * computed value at bottom of model loop */
    new_vertex_value = chisqmin;


    if (comparison == SIMPLEX) {
      simplex(nParameters,
         new_vertex_value,new_vertex_position);
    }

    if (comparison == ANNEAL) {
      status = anneal(nParameters, anneal_t0, position_limits,
	 new_vertex_position, new_vertex_value, model_number, chisqbest);
         if (status == 0){
      printf("anneal returned failed\n");
      fprintf(fpnotes,"anneal returned failed\n");
            goto exitnow;
         }
    }

/* Here the master defines the parameters for the new model using the
   vertex values that were just set by the simplex or simulated annealing
   functions. Those functions process non-dimensional numbers that
   define the simplex polygon, for example. The next function set_model_parameters
   turns those into constants or power law indices to define temperature, density,
   etc.
*/

   if (comparison != NO) 
        set_model_parameters(nParameters,new_vertex_position,position_limits,
	    model_number,my_rank);
/*        printf("Standard abundance in cmain %e \n",std_abundance); 
        fprintf(fpnotes,"Standard abundance in cmain %e \n",std_abundance);
*/


/* Here we actually compute the temperature, density, etc for each i,j,k cell and
   write those values into memory. This is done in two cases.
   First: If we are using number_models to write out intermediate results then
   are using the same physical model throughout the calculation. In this case
   we compute the model structure only once when model_number == 0.
   Second: If we are doing comparisons then we are computing a different model
   structure for each model_number.
*/
  if (model_number == 0 || comparison != NO) {
    j = define_model(model_number, nbox, nx, ny, nz,
                 xcenter,  ycenter,  zcenter, 
                 atoms,lwmin,lwmax,depletion,photodissociation,std_abundance,molecule) ;
    if (j == 0) {
      printf("define_model returned failed\n");
      fprintf(fpnotes,"define_model returned failed\n");
      printf("Going to exit ...\n");
      fprintf(fpnotes,"Going to exit ...\n");
      goto exitnow;
    }

	printf("Finished defining model. Chemistry next.\n"); 

/* Now that the model is defined with its density, compute the mean visual extinction 
   to each cell. Run through all the rays to get the visual extinction along each
   ray and accumulate the mean in each pixel. To do this we have to re-compute the 
   pathlengths along all the rays. First, zero out the mean_av in each pixel.
*/
        for (ibox=0;ibox<nbox;ibox++){
        for (i=0;i<nx[ibox];i++) {
        for (j=0;j<ny[ibox];j++) {
        for (k=0;k<nz[ibox];k++) { 
                model[ibox][i][j][k].mean_av = 0;
        }}}}

// Now loop over all the rays in the usual way to identify ibox,i,j,k for each segment
        for (iray=0;iray<nrays;iray++){
                result = raypath_(&x1[iray],&y1[iray],&z1[iray],
                         &x2[iray],&y2[iray],&z2[iray],
                         &nbox,nx,ny,nz,cellsizex,cellsizeyz,
                         &maxgrid,t4,t5,t6,
                         &maxvox,t7,t8,&i);
                if (result == 0) {
                        printf("Got zero as a result from raypath\n");
                        goto exitnow;
                }
                if (i-1 != icount[iray]) {
                    printf("Loop 1, not getting the same ray segments %d %d\n",i-1,icount[iray]);
                    goto exitnow;
                }
		if (icount[iray] == 0) {
			fprintf(fpnotes,"2 Bad icount from raypath %d %d\n",icount[iray],iray);
			printf("2 Bad icount from raypath %d %d\n",icount[iray],iray);
		}

                for (ivox=0;ivox<icount[iray];ivox++){
                        pathlength[ivox] = t8[ivox];
                for (i=0;i<4;i++){
                        k = 4*ivox + i;
                        iarc[i][ivox] = t7[k]-1;
                }}
/*  Initialize the column density to zero at the start of each ray */
/*  The pathlength is in units of PC. To convert to Av, divide by 9.4e20 */
                column_density_sum = 0.;        
                for (ivox=0;ivox<icount[iray];ivox++){
                        i = iarc[0][ivox];
                        j = iarc[1][ivox];
                        k = iarc[2][ivox];
                        ibox = iarc[3][ivox];
                        column_density_sum += (PC/9.4e20)*pathlength[ivox]*model[ibox][i][j][k].density;
                        column_density_segment[ivox] = column_density_sum;
/* fprintf(fpnotes,"iray %d ivox %d path %e density %e cd sum %e cd seg %e\n",
    iray,ivox,pathlength[ivox], model[ibox][i][j][k].density,column_density_sum,
    column_density_segment[ivox]); */
                }
    
/* We have the column densities for each segment for the forward direction 
   from above and the total sum for that ray.  Now add the extinctions for 
   both the forward and reverse directions along the ray, normalize by 
   the number of rays through each model cell.  The result is the mean visual 
   extinction for each cell. In order to make the Av symmetric, the 
   case ivox=0 has to be done with its own equation. The cases ivox>0 
   use [ivox-1] in the equation.

   The visual extinction is converted into transmission by the exponential
   transmission = exp(-av). The transmission is the variable mean_av
   and the extinction is the variable column_density_sum or _segment is the
   visual extinction, the column density itself is Sum(PC*pathlength*density)
   along the ray. 
*/
		  ivox = 0;
                  i = iarc[0][ivox];
                  j = iarc[1][ivox];
                  k = iarc[2][ivox];
                  ibox = iarc[3][ivox];
                  model[ibox][i][j][k].mean_av += 
                     (exp(-column_density_segment[0])
                     + exp(-column_density_sum ))
                     / model[ibox][i][j][k].adds;

                for (ivox=1;ivox<icount[iray];ivox++){
                  i = iarc[0][ivox];
                  j = iarc[1][ivox];
                  k = iarc[2][ivox];
                  ibox = iarc[3][ivox];
                  model[ibox][i][j][k].mean_av += 
                     (exp(-column_density_segment[ivox])
                     + exp(-column_density_sum + column_density_segment[ivox-1]))
                     / model[ibox][i][j][k].adds;
/* fprintf(fpnotes,"iray ivox %d %d %10.3e %10.3e %10.3e\n",iray,ivox,(-column_density_segment[ivox]),
(-column_density_sum + column_density_segment[ivox]),column_density_sum ); */
                }
        } /* end of loop over rays */

  printf("End of loop over rays. Start chemistry\n"); 

  if (depletion == 1 || photodissociation == 1){
    chemistry(depletion,photodissociation,molecule,nbox,nx,ny,nz);
    printf("Finished chemistry. Now check the model.\n"); 
  } else {
    printf("No chemistry. Now check the model.\n"); 
  }

    i = 1;
    i = check_model(model_number, nbox, nx, ny, nz,
                 xcenter,  ycenter,  zcenter,
                 lwmin,lwmax) ;
   fflush(fpnotes); 
   if (i == 0) {
        printf("Going to exit ...\n");
        goto exitnow;
   }


/* This initializes the radiation field.
   Do this for only the first model if we are using number_models to write
   out intermediate results. Do this for every model if we are doing comparisons.
   Start with the zero, CMB, or LTE according to the value of initial_radiation.
   initial_radiation:  0 = ZERO and 1 = CMB and 2 = LTE 
*/
   if (model_number == 0 || comparison != 0) {
	printf("Finished model check. Initialize the radiation field.\n"); 
	i = initrad_(&nbox,nx,ny,nz,&initial_radiation);
	printf("Finished initializing the radiation field.\n"); 
   }

	gettimeofday(&splittime3, &zone);
	elapsed = (splittime3.tv_sec  - starttime.tv_sec ) * 1000
                + (splittime3.tv_usec - starttime.tv_usec) / 1000;
	fprintf(fpnotes,"P%d: Finished model structure at %ld msecs\n",
    		my_rank,elapsed);
	printf (        "P%d: Finished model structure at %ld msecs\n",
    		my_rank,elapsed);

   } /* This is the end of the if block that initializes the model structure if
        we are at the first model (model_number == 0) or if we are doing 
        comparisons. 
     */


/* --------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------- */

/* Iterate between statistical equilibrium equations and radiative transfer */

/* Before computing a new model make sure to
 * reset the output cube velocity dimensions and channel velocities
 * which might have been averaged down by hanning smoothing at the end of the 
 * last simulation. Re-zero the output data cube brightness and weights
 * to accumulate sums.  */


    for (i=0;i<nllist;i++){
      dcube.nchan[i] = nchan[linelist[i]];
      fprintf(fpnotes,"Reset nchan for line %d to %d\n",i, dcube.nchan[i]);
      for (n=0;n<nchan[linelist[i]];n++) {
        dcube.chvel[i][n] = chvel[linelist[i]][n];
/*        fprintf(fpnotes,"%d %f\n",n,1.e-5*dcube.chvel[i][n]); */
      }
    }

          for (i=0;i<dcube.nviews;i++) {
          for (j=0;j<dcube.nx;j++) {
          for (k=0;k<dcube.ny;k++) {
            dcube.wt[i][j][k] = 0.;
          for (m=0;m<dcube.nlines;m++) {
          for (n=0;n<dcube.nchan[m];n++) {
            dcube.cube[m][i][j][k][n] = 0.;
          }}}}}


    for (i=0;i<dcube.nx;i++) {
    for (j=0;j<dcube.ny;j++) {
    }}

/* Open a new convergence file for the this simulation */

    close(ipconv);
    sprintf(suffix,"%d",model_number);
             strcpy(fname,pathname);
             strcat(fname,"Convergence");
	     strcat(fname,suffix);
             fpconv = fopen(fname,"w");
             if (fpconv == NULL) {
                printf("Could not open %s. Exiting ...\n",fname);
                exit(0);
             } else {
                printf("File %s opened\n",fname);
             }
             ipconv = fileno(fpconv);
/*    printf("Convergence file number is %d\n",ipconv); */


/* Here is the start of the lambda iteration loop */

    iter = 1;
    while (iter<=iterations) {

    printf("P%d: Iteration number %d out of %d\n",my_rank,iter,iterations);
/*		printf("P%d: nchan %d %d %d %d\n",my_rank,nchan[0],nchan[1],nchan[2],
			nchan[3]); */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* The first step is to compute the stat. eq. equations. The actual computation
   is done in fsub.f. The function stateq below calls one of the fortran functions,
   either accleratedlambda or one of the LTE functions as appropriate. The non-LTE
   calculation requires the mean radiation field as input in the stat. eq. equations. 
   For the first iteration, the main radiation field is set by the fortran function 
   initrad called earlier. On subsequent iterations it is carried over from the 
   previous iteration in the variable apLJbar.Jbar.
*/

/* Lambda acceleration only starts with iteration 3 */

	accel = 0;
        if (acceleration) accel = 1;	/* If accel = 0, no acceleration */
	if (iter <= 3 && (model_number == 0 || comparison != 0) ) {
           accel = 0;
        }

	gettimeofday(&splittime1, &zone);
	printf("Starting stat. eq.\n");  

/*	fprintf(fpnotes,"ijk3 %d %d %d line %d ihyp %d xn %e\n",
            16,16,16,1,1,apLJbar[0][16][16][16].norm[1][1]);  */

	  status = stateq(nbox,nx,ny,nz,iter,accel,
                   molecule,dust_continuum,initial_radiation);
printf("Finished stat. eq. with status = %d\n",status);
	  if (status == 0) goto exitnow;


	  convergence(nbox,nx,ny,nz,iter,&converged);

/*	  printf("Converged: %d\n",converged); */
	  fprintf(fpnotes,"Converged: %d\n",converged);
//	printf("finished stat. eq. \n"); 
//	fprintf(fpnotes,"finished stat. eq. \n");


/* If the convergence criteria are satisfied, then set the maximum number
   of iterations to the current number so that the last rt calculation
   will be written to disk.  */
/*
		if (converged == 1) iterations = iter;
*/
/*		printf("iterations %d    iter %d\n",iterations,iter);*/


/* Print out some stuff */
/* The definition of opacity: What is calculated as the opacity is the 
opacity formula without the line profile function. In other words, 
the srcopc.opacity is independent of frequency and should be multiplied 
by the line profile function. Since we don't know this function in 
this part of the code, we can divide by sqrtpi*linewidth(Hz).  
The reason is because this factor is the normalization of the line 
profile function, or the value it would have at the line peak. 
Here is the formula with the normalization times the exponential.
phi(nu) = 1/(sqrtpi delta nu) * exp(((nu - nu_channel)/delta nu)^2)
At any channel, the line profile function is always much less than 1. So 
it is correct to divide. */
    if (iter == iterations || iter >= 1 || 1) {
      fprintf(fpnotes,"NewIteration %d\n",iter);

/* As a representative point, take the middle ray of all the rays in the 
   first view. Then run the raypath function to get the segments along 
   this line of sight.
*/

	iray = minray;
	fprintf(fpnotes,"AtPosition %d %f %f\n",iray,y1[iray],z1[iray]);

                result = raypath_(&x1[iray],&y1[iray],&z1[iray],
                         &x2[iray],&y2[iray],&z2[iray],
                         &nbox,nx,ny,nz,cellsizex,cellsizeyz,
                         &maxgrid,t4,t5,t6,
                         &maxvox,t7,t8,&i);
                if (result == 0) {
                        printf("Error from raypath for ray number %d\n",iray);
                        goto exitnow;
                }
		if (icount[iray] == 0) {
			fprintf(fpnotes,"3 Bad icount from raypath %d %d\n",icount[iray],iray);
			printf("3 Bad icount from raypath %d %d\n",icount[iray],iray);
		}
                for (ivox=0;ivox<icount[iray];ivox++){
                        pathlength[ivox] = t8[ivox];
                for (i=0;i<4;i++){
                        k = 4*ivox + i;
                        iarc[i][ivox] = t7[k]-1;
                }}
 


	for (l=0;l<NLINES;l++){
          fprintf(fpnotes,"LineTransition %d %f\n",
			l,atomic.freq[l]/1.e9);
          fprintf(fpnotes,"box  ix  xcenter    source(K)    opacity(cm-1) cell tau    Jbar(K)     apL\n");

	  sumopac = 0.;

/* Here are the voxels along this ray */
                for (ivox=0;ivox<icount[iray];ivox++){
                  i = iarc[0][ivox];
                  j = iarc[1][ivox];
                  k = iarc[2][ivox];
                  ibox = iarc[3][ivox];

/* Set cellopac to zero so that we can sum over hyperfines, if there are any */
            cellopac = 0.;

/* For N2H+ with non-LTE hyperfines, average the optical depth as for Einstein A */

	  if (molecule == 4) {
            for (ihyp=0;ihyp<atomic.nhyp[l];ihyp++){
           		it = atomic.start[l] + ihyp;
			qfu = atomic.ln[atomic.upper[it]];
			qju = atomic.jn[atomic.upper[it]];
			statwt = (2.0*qfu + 1.0) / (9.0*(2.0*qju + 1.0));
	       cellopac = cellopac + srcopc[ibox][i][j][k].opacity[l][ihyp]
		  /SQRTPI/(atomic.freq[l]/C*model[ibox][i][j][k].linewidth)
		  * statwt * PLANCK*atomic.freq[l]/(4.*PI);
            }
	  } 

/* For HCNhyp with non-LTE hyperfines, average the optical depth as for Einstein A */

	  if (molecule == 21) {
            for (ihyp=0;ihyp<atomic.nhyp[l];ihyp++){
           		it = atomic.start[l] + ihyp;
			qfu = atomic.ln[atomic.upper[it]];
			qju = atomic.jn[atomic.upper[it]];
			statwt = (2.0*qfu + 1.0) / (3.0*(2.0*qju + 1.0));
	       cellopac = cellopac + srcopc[ibox][i][j][k].opacity[l][ihyp]
		  /SQRTPI/(atomic.freq[l]/C*model[ibox][i][j][k].linewidth)
		  * statwt * PLANCK*atomic.freq[l]/(4.*PI);
            }
	  } 


	if (nltehyp == 0) {
/* For molecules without hyperfine structure or with LTE hyperfines*/
	       cellopac = srcopc[ibox][i][j][k].opacity[l][0]
		  /SQRTPI/(atomic.freq[l]/C*model[ibox][i][j][k].linewidth)
		* PLANCK*atomic.freq[l]/(4.*PI);
	  }

	    sumopac = sumopac + cellopac
		* model[ibox][i][j][k].density
		* cellsizex[ibox]*PC;
            celltau = cellopac
                * model[ibox][i][j][k].density
                * cellsizex[ibox]*PC;


/* fprintf(fpnotes,"sumopac cellopac density cell %12.4e %12.4e %12.4e\n",
	sumopac,cellopac,model[ibox][i][j][k].density,cellsizex[ibox]); */
            
            apLsum = 0.;
            barJsum = 0.;
	    emissum = 0.;
	    opacsum = 0.;
	    nh = 1; if (nltehyp == 1) nh = atomic.nhyp[l];
/*	fprintf(fpnotes,"nh for print %d\n",nh);  */
            for (ihyp=0;ihyp<nh;ihyp++){
		apLsum  += apLJbar[ibox][i][j][k].apL[l][ihyp];
		barJsum += apLJbar[ibox][i][j][k].Jbar[l][ihyp];
		emissum += srcopc[ibox][i][j][k].emissivity[l][ihyp];
		opacsum += srcopc[ibox][i][j][k].opacity[l][ihyp];
/*      fprintf(fpnotes,"l ihyp %d %d opac %e\n",
          l,ihyp,srcopc[ibox][i][j][k].opacity[l][ihyp]); */
            } 

	if (opacsum != 0.0) {
            srcsum = emissum/opacsum / (2.*BOLTZ) * pow((C/atomic.freq[l]),2);
	} else {
	    srcsum = 0.0;
	}
            fprintf(fpnotes,"%3d %3d %9.5f %12.4e %12.3e %12.4e %12.4e %12.4e\n",
                ibox,i,xcenter[ibox][i],srcsum,
                cellopac, celltau, barJsum/nh, apLsum/nh);
          } /* End of loop along the ray */
	  fprintf(fpnotes,"Avg. opt.dpth. per cell %g     Total optical depth  %g\n",
		sumopac/nx[ibox],sumopac);
	} /* End of loop on lines */
    } /*  End of IF block: if (iter == iterations ...  */
/* Done printing */

	gettimeofday(&splittime2, &zone);
	elapsed = (splittime2.tv_sec  - splittime1.tv_sec ) * 1000
                + (splittime2.tv_usec - splittime1.tv_usec) / 1000;
	fprintf(fpnotes,"P%d: Stat eq equations required %ld msecs\n",
    		my_rank,elapsed);

	printf (        "P%d: Stat eq equations required %ld msecs\n",
    		my_rank,elapsed);


/* Done with all the S.E. calculations. */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* Compute the radiative transfer */

/* Figure out how many lines to compute */
/* On the last iteration, only compute the nllist lines of 
   interest that are stored in the linelist. On all other
   iterations, compute all the lines.  */
	if (iter == iterations && (model_number == number_models || comparison != 0)) {
		l1 = 0;
		l2 = nllist;
        } else {
		l1 = 0;
		l2 = NLINES;
	}

	printf ("P%d: Starting Radiative Transfer for %d lines\n", 
               my_rank,l2 - l1 ); 
	fprintf (fpnotes,"P%d: Starting Radiative Transfer for %d lines \n",
                my_rank,l2 - l1 ); 

	for (lineCount=l1;lineCount<l2;lineCount++){

	  if (iter == iterations && (model_number == number_models || comparison != 0)) {
            mline = linelist[lineCount];
          } else {
            mline = lineCount;
          }

/* The master is currently working on MLINE */

	gettimeofday(&splittime1, &zone);

/* Compute a new Jbar and approx lambda term for the current line.
   Zero out the variables in preparation for the accumulation of sums
*/
	nh = 1; if (nltehyp == 1) nh = atomic.nhyp[mline];
        for (ibox=0;ibox<nbox;ibox++){
        for (i=0;i<nx[ibox];i++){
        for (j=0;j<ny[ibox];j++){
        for (k=0;k<nz[ibox];k++){
        for (ihyp=0;ihyp<nh;ihyp++){
        apLJbar[ibox][i][j][k].apL[mline][ihyp]  = 0.;
        apLJbar[ibox][i][j][k].Jbar[mline][ihyp] = 0.;
        apLJbar[ibox][i][j][k].norm[mline][ihyp] = 0.;
	}}}}}

/* ------------------------------------------------------------ */
/* Figure out which rays we want to compute.  On the last iteration
   of the last model
   only compute the rays in the views that are requested in the output.
   This section puts the index numbers of the rays to be computed
   into a list, raylist. For all iterations but the last, raylist
   contains the numbers, 0 -> nrays. On the last iteration, the
   list could be shorter if not all the views are to be written
   to output.
*/

	if (iter == iterations && (model_number == number_models || comparison != 0)) {
          nsets = outviews;
        } else {
          nsets = 1;
        }


/*        printf("iteration = %d Computing %d sets of views\n",iter,nsets); */
        fprintf(fpnotes,"iteration = %d Computing %d sets of views\n",iter,nsets);


	j = 0;
        for (iset=0;iset<nsets;iset++){

	  if (iter == iterations && (model_number == number_models || comparison != 0)) {
            istart = outray1[iset] ;
            iend   = outray1[iset] + outrpv[iset]; /* OK because we loop on i < iend */
          } else {
            istart = 0;
            iend   = nrays;
          }
          for (iray=istart;iray<iend;iray++){
	    raylist[j] = iray;
	    j++;
	  }

/*        printf("Computing rays from %d to %d \n",istart,iend-1); */
        fprintf(fpnotes,"Computing rays from %d to %d \n",istart,iend-1);

	} /* end of iset loop */
	rays_to_compute = j-1; /* This looks OK because just below 
				the while loop runs to <= rays_to_compute */
        fprintf(fpnotes,"rays_to_compute %d \n",rays_to_compute);

/* Mark up the send and receive checklists. */

for (i=0;i<nrays;i++) {
	sendcheck[i] = 1;
	recvcheck[i] = 1;
}

/*
	printf("rays to compute % d\n",rays_to_compute);
	fprintf(fpnotes,"rays to compute % d\n",rays_to_compute);
*/

/*	 printf("nchan %d %d %d %d\n",nchan[0],nchan[1],nchan[2],
			nchan[3]); */


/* The list of rays to be computed is now in raylist 
   The number of rays to be computed is rays_to_compute
*/
/* ------------------------------------------------------------ */
/* the start of the parallelization
*/

/* Now we run through the list of rays, sending each ray to a free processor. 
   Start with some initializations
*/
	rays_finished = 0;
	current_ray = 0;
        master_status = 1;

	free_processes = processes - 1;
	first_proc = 1;



/* A ray is finished when the master receives the computed radiation fields
   from the slaves. This while loop is over all rays. Inside this loop
   are sections to send and receive data from the slaves.
*/
        while (rays_finished <= rays_to_compute) { 

/*
	fprintf(fpnotes,"P%d: current ray %d, ray from list %d\n",
		my_rank,current_ray,raylist[current_ray]);
*/

/* This section to send data to slaves. It is executed for each free processor
   as long as there remain rays in the list to send. First we find the 
   minimum of the number of rays and the number of free processors (slaves). 
   This is the number of times to execute the loop before waiting for the slave
   processes to report back. If there are more available processes than 
   rays to compute then reduce the number of "free" processes to the number 
   of rays because we won't need all the extra processes.
*/

/*if (iter==iterations) fprintf(fpnotes,"rays_to_compute %d current_ray %d\n",
      rays_to_compute,current_ray);
*/

	if ((rays_to_compute - current_ray + 1) < free_processes)
		free_processes = (rays_to_compute - current_ray + 1);

	last_proc = first_proc + free_processes;


/*
if (iter==iterations){
	fprintf(fpnotes,"P%d: first proc %d last proc %d free processes %d\n",my_rank,
		first_proc,last_proc,free_processes);
}
*/



/* safety to be sure the following loop does not run once if there is nothing to do */
	if (last_proc > first_proc) {      

	gettimeofday(&splittimeC, &zone);

	        for (p=first_proc; p<last_proc; p++) {

/* Here is the ray from the list to send. Load up all the vectors along this ray. */
		send_iray = raylist[current_ray];
/* if (iter==iterations) fprintf(fpnotes,"send_iray current_ray nrays %d %d %d\n",
	send_iray,current_ray,nrays); 
*/

/* ---------------------------------------------------------------------- */

/* The master now sends the pathlength, emissivity, opacity, 
   velocity, and linewidth for one line and for each voxel along
   one ray to a slave. Also sends the pathlength for each voxel.

   This section recomputes the segments for each ray and loads up the 
   model data corresponding to each segment. 
*/

		result = raypath_(&x1[send_iray],&y1[send_iray],&z1[send_iray],
		         &x2[send_iray],&y2[send_iray],&z2[send_iray],
                         &nbox,nx,ny,nz,cellsizex,cellsizeyz,
                         &maxgrid,t4,t5,t6,
                         &maxvox,t7,t8,&i);

                if (result == 0) {
			printf("Got zero as a result from raypath\n");
			master_status = 0;
			current_ray = rays_to_compute;
			goto end_parallel;
		}
                if (i-1 != icount[send_iray]) {
                    printf("Loop 2, not getting the same ray segments %d %d\n",i-1,icount[send_iray]);
                    goto end_parallel;
                }
		if (icount[send_iray] == 0) {
			fprintf(fpnotes,"4 Bad icount from raypath %d %d\n",icount[send_iray],send_iray);
			printf("4 Bad icount from raypath %d %d\n",icount[send_iray],send_iray);
		}

		for (ivox=0;ivox<icount[send_iray];ivox++){
			pathlength[ivox] = t8[ivox];
		for (i=0;i<4;i++){
			k = 4*ivox + i;
			iarc[i][ivox] = t7[k]-1; 
		}}

/* ---------------------------------------------------------------------- */
	dest = p;	/* Which processor the data will go to  */

/* Pack up the data into 2 arrays for transmission as ints and floats */
/* A short array of 5 ints */
/* Using mline for master:line as an index instead of just line */
	q = 0;
	comm_ints[q++] = mline;			/* 0 */
	comm_ints[q++] = icount[send_iray];	/* 1 */
	comm_ints[q++] = send_iray;		/* 2 */
/* 1 = active status, 0 = ending program */
	comm_ints[q++] = 1;			/* 3 */
	comm_ints[q++] = iter;	 		/* 4 */

/*
printf("Comm ints %d %d %d %d %d\n",mline,icount[send_iray],send_iray,1,iter);
printf("Comm ints %d %d %d %d %d\n",comm_ints[0],comm_ints[1],comm_ints[2],comm_ints[3],comm_ints[4]);
*/

/* All the float data are packed into one array */
/* First pack up the 3 angles of the ray */
        q = 0;
	comm_floats[q++] = dxdu[send_iray];
	comm_floats[q++] = dydu[send_iray];
	comm_floats[q++] = dzdu[send_iray];

/* Next pack up the variables one by one */
	for (ivox=0;ivox<icount[send_iray];ivox++) {
                i = iarc[0][ivox];
                j = iarc[1][ivox];
                k = iarc[2][ivox];
                ibox = iarc[3][ivox];
		comm_floats[q++] = model[ibox][i][j][k].vx;
	}
	for (ivox=0;ivox<icount[send_iray];ivox++) {
                i = iarc[0][ivox];
                j = iarc[1][ivox];
                k = iarc[2][ivox];
                ibox = iarc[3][ivox];
		comm_floats[q++] = model[ibox][i][j][k].vy;
	}
	for (ivox=0;ivox<icount[send_iray];ivox++) {
                i = iarc[0][ivox];
                j = iarc[1][ivox];
                k = iarc[2][ivox];
                ibox = iarc[3][ivox];
		comm_floats[q++] = model[ibox][i][j][k].vz;
	}
	for (ivox=0;ivox<icount[send_iray];ivox++) {
                i = iarc[0][ivox];
                j = iarc[1][ivox];
                k = iarc[2][ivox];
                ibox = iarc[3][ivox];
		comm_floats[q++] = model[ibox][i][j][k].linewidth;
	}
/* Send the column density per cell in cgs units */
	for (ivox=0;ivox<icount[send_iray];ivox++) {
                i = iarc[0][ivox];
                j = iarc[1][ivox];
                k = iarc[2][ivox];
                ibox = iarc[3][ivox];
		comm_floats[q++] = pathlength[ivox]*model[ibox][i][j][k].density*PC;
/* fprintf(fpnotes,"comm_floats colmn d %d %12.3e \n",q-1,comm_floats[q-1]); */
	}
/*	for (ivox=0;ivox<icount[send_iray];ivox++) {
                i = iarc[0][ivox];
                j = iarc[1][ivox];
                k = iarc[2][ivox];
                ibox = iarc[3][ivox];
		comm_floats[q++] = model[ibox][i][j][k].abundance;
	} */
	for (ivox=0;ivox<icount[send_iray];ivox++) {
                i = iarc[0][ivox];
                j = iarc[1][ivox];
                k = iarc[2][ivox];
                ibox = iarc[3][ivox];
		comm_floats[q++] = srcopc[ibox][i][j][k].cont_emissivity[mline];
/* fprintf(fpnotes,"comm_floats dust em %d %12.3e \n",q-1,comm_floats[q-1]); */
	}

        for (ivox=0;ivox<icount[send_iray];ivox++) {
                i = iarc[0][ivox];
                j = iarc[1][ivox];
                k = iarc[2][ivox];
                ibox = iarc[3][ivox];
                comm_floats[q++] = srcopc[ibox][i][j][k].cont_opacity[mline];
/* fprintf(fpnotes,"comm_floats dust op %d %12.3e \n",q-1,comm_floats[q-1]); */
        }

	n = 0;
	for (ivox=0;ivox<icount[send_iray];ivox++) {
                i = iarc[0][ivox];
                j = iarc[1][ivox];
                k = iarc[2][ivox];
                ibox = iarc[3][ivox];
                for (ihyp=0;ihyp<nh;ihyp++){
                        comm_floats4[n++] = srcopc[ibox][i][j][k].emissivity[mline][ihyp];
                        q++;
/* fprintf(fpnotes,"comm_floats emissiv %d %12.3e \n",q-1,comm_floats[q-1]); */
                }
	}
        n = 0;
	for (ivox=0;ivox<icount[send_iray];ivox++) {
                i = iarc[0][ivox];
                j = iarc[1][ivox];
                k = iarc[2][ivox];
                ibox = iarc[3][ivox];
                for (ihyp=0;ihyp<nh;ihyp++){
                        comm_floats5[n++] = srcopc[ibox][i][j][k].opacity[mline][ihyp];
                        q++;
/* fprintf(fpnotes,"comm_floats opacity %d %12.3e \n",q-1,comm_floats[q-1]); */
                }
	}

/* Send the int data. n is the number of variables to send */
/* tag=0 is for the first message with comm_ints, tag=1 for the
   second message with comm_floats */
        n = 5;
        tag = 0;
	MPI_Send(comm_ints,  n,MPI_LONG,  dest,tag,MPI_COMM_WORLD);

/*fprintf(fpnotes,"P%d sent int ray %d for line %d to P%d\n",my_rank,send_iray,line,dest);*/

/* q is the number of variables packed into the float array */
/* In the syntax comm_floats[q++], the increment to q is done after q is used
   as the index. Also the last increment, q++, means that the final q is equal
   to the number of variables which is 1 more than the index of the last variable.
*/
   
       i = 3 + (7+2*nh)*icount[send_iray];
       if (i != q)
          printf("Help, master cant count, counted = %d, calculated = %d\n",q,i);
        tag = 1;
//	MPI_Send(comm_floats,q,MPI_FLOAT,dest,tag,MPI_COMM_WORLD);
	ic = icount[send_iray];	
        n = 3 + 3*ic;
//	printf("Master %d Sending q = %d dest = %d tag = % d\n",my_rank,n,dest,tag); 
        MPI_Send(&comm_floats[0],n,MPI_FLOAT,dest,tag,MPI_COMM_WORLD);
	tag = 2;
        n = 2*ic;
//	printf("Master %d Sending q = %d dest = %d tag = % d\n",my_rank,n,dest,tag); 
        MPI_Send(&comm_floats[3+3*ic],n,MPI_FLOAT,dest,tag,MPI_COMM_WORLD);
	tag = 3;
        n = 2*ic;
//	printf("Master %d Sending q = %d dest = %d tag = % d\n",my_rank,n,dest,tag); 
        MPI_Send(&comm_floats[3+5*ic],n,MPI_FLOAT,dest,tag,MPI_COMM_WORLD);
	tag = 4;
        n = nh*ic;
//	printf("Master %d Sending q = %d dest = %d tag = % d\n",my_rank,n,dest,tag); 
        MPI_Send(&comm_floats4[0],n,MPI_FLOAT,dest,tag,MPI_COMM_WORLD);
//        MPI_Send(&comm_floats[3+7*ic],n,MPI_FLOAT,dest,tag,MPI_COMM_WORLD);
	tag = 5;
        n = nh*ic;
//	printf("Master %d Sending q = %d dest = %d tag = % d\n",my_rank,n,dest,tag); 
        MPI_Send(&comm_floats5[0],n,MPI_FLOAT,dest,tag,MPI_COMM_WORLD);
//        MPI_Send(&comm_floats[3+(7+nh)*ic],n,MPI_FLOAT,dest,tag,MPI_COMM_WORLD);

/*fprintf(fpnotes,"P%d sent flt ray %d for line %d to P%d\n",my_rank,send_iray,line,dest);*/

/* Sent this ray. Mark the checklist */
	sendcheck[send_iray]--;

/* if (mline == 1 ) goto exitnow; */

	current_ray++; 		/* ready to send the next ray */
	free_processes--;	/* used up a process */


/*	fprintf(fpnotes,"P%d: current ray now %d, free_processes %d\n",
		my_rank,current_ray,free_processes);
*/


	} /* end of for loop p=first_proc to last_proc */
	} /* end of safety if block with last_rays > first_ray */

/* Careful, there are two curly brackets above that are not comments, but required */

/* Now there are no free processes */

/* Wait here for data to come back from slave.  */
/* The master receives the lamda factor, the radiation field barJ, and
   the spectral line from the slave. This is data for one line, including
   hyperfines, for each voxel along the ray. */


/*	fprintf(fpnotes,"P%d: waiting for data from a slave \n",my_rank); */


	n = 5;
        tag = 0;			/* tag=0 indicates that we want a comm_ints message 	*/
					/* MPI_ANY_SOURCE receives from any slave 		*/
	MPI_Recv(comm_ints,n,MPI_LONG,MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,&mpi_status); 

        orig = mpi_status.MPI_SOURCE;	/* This is process number of the slave */

        sline = comm_ints[0];		/* the line */
        if (sline != mline) printf("Lost my line number! %d %d \n",sline,mline);
	ic = comm_ints[1];		/* The length of the ray */
        recvd_iray = comm_ints[2];	/* the number of the ray */

        if (ic != icount[recvd_iray]) {
               printf("Some problem with location of received ray, %d %d\n",
                      ic,icount[recvd_iray]);
		master_status = 0;
		current_ray = rays_to_compute;
		goto end_parallel;
        }

/* We received this ray. Mark the checklist */
	recvcheck[recvd_iray]--;

/* fprintf(fpnotes,"P%d got ray %d for line %d from P%d\n",
my_rank,recvd_iray,mline,orig);
*/

	nh = 1; if (nltehyp == 1) nh = atomic.nhyp[mline];
	q = 3*ic*nh + nchan[sline];
					/* q = number of variables to receive in comm_floats */
        tag = 1;			/* tag=1 indicates that we want a comm_floats message 	*/
					/* orig indicates that we want the message from the   	*/
                                        /* slave that just sent us the comm_ints, not from    	*/
					/* another slave 					*/

/*	printf("Master Receiving q = %d from %d\n",q,orig); */
	MPI_Recv(comm_floats,q,MPI_FLOAT,orig,tag,MPI_COMM_WORLD,&mpi_status);
	dest = mpi_status.MPI_SOURCE;

/*	fprintf(fpnotes,"P%d: received data from slave %d\n",my_rank,dest); */

/*for (i=0;i<q;i++)
fprintf(fpnotes,"recv i %d cf %e\n",i,comm_floats[i]);*/

	if (comm_ints[3] == 0)  {
	    printf("Master received an error condition from slave %d working on ray %d\n"
		,dest,recvd_iray);
	    printf("Look for errors from slaves in file Fnotes1 and/or fort.22 \n");
	    fprintf(fpnotes,"Master received an error condition from slave %d working on ray %d\n"
		,dest,recvd_iray);
	    fprintf(fpnotes,"Look for errors from slaves in file Fnotes1 and/or fort.22 \n");
			master_status = 0;
			current_ray = rays_to_compute;
			goto end_parallel;
        }

/* Unpack the comm float array into apLray, barJray, and yspec. For the first 2 arrays,
   the NHYP index varies most rapidly (inner loop). */


	n = ic*nh;
	q = 0;
	for (j=0;j<ic;j++) {
	for (i=0;i<nh;i++) {
		 apLray[i][j] = comm_floats[q];
/*		fprintf(fpnotes,"unpack ivox ihyp aL %d %d %e\n",j,i,apLray[i][j]);  */
                q++;
	}} 
	if (q != n) printf("Master lost count in unpacking apL %d %d\n",q,n);

	for (j=0;j<ic;j++) {
	chksum = 0.;
	for (i=0;i<nh;i++) {
		barJray[i][j] = comm_floats[q];
/*		fprintf(fpnotes,"unpack ivox ihyp bJ %d %d %e\n",j,i,barJray[i][j]);  */
		chksum = chksum + comm_floats[q];
		q++;
	}
	if (chksum == 0. && (initial_radiation != 0 || iter != 1)) {
/*	   printf("chksum initial_radiation iter %f %d %d\n",chksum,initial_radiation,iter); */
           fprintf(fpnotes,"Jbar unpacked to zero for line=%d, iter=%d, i=%d,j=%d,k=%d\n",
                     mline,iter,i,j,k);
           printf("Jbar unpacked to zero for line=%d, iter=%d, i=%d,j=%d,k=%d\n",
                     mline,iter,i,j,k);
	}

	} 
        if (q != 2*n) printf("Master lost count in unpacking bJ  %d %d\n",q,n);

	for (j=0;j<ic;j++) {
	chksum = 0.;
	for (i=0;i<nh;i++) {
		normray[i][j] = comm_floats[q];
/*		fprintf(fpnotes,"unpack ivox ihyp bJ %d %d %e\n",j,i,barJray[i][j]);  */
		chksum = chksum + comm_floats[q];
		q++;
	}
	if (chksum == 0. && (initial_radiation != 0 || iter != 1)) {
/*	   printf("chksum initial_radiation iter %f %d %d\n",chksum,initial_radiation,iter); */
	for (j=0;j<ic;j++) {
	for (i=0;i<nh;i++) {
           fprintf(fpnotes,"Norm unpacked to zero for line=%d, iter=%d, i=%d,j=%d,k=%d %e\n",
                     mline,iter,i,j,k,normray[i][j]);
           printf("Norm unpacked to zero for line=%d, iter=%d, i=%d,j=%d,k=%d %e\n",
                     mline,iter,i,j,k,normray[i][j]);
        }}
	}

	} 
        if (q != 3*n) printf("Master lost count in unpacking bJ  %d %d\n",q,n);


	for (j=0;j<nchan[mline];j++){
		yspec[j] = comm_floats[q];
                q++;
	}


// This block for printing debugging
        if ((recvd_iray == -147 || recvd_iray == -148) && (mline == 0 || mline == 1)) {
      printf("Y test C: %d %d %d %d %d\n",
         recvd_iray,mline,ic,nchan[mline]/2-4,nchan[mline]/2+2);
      fprintf(fpnotes,"Y test C: %d %d %d %d %d\n",
         recvd_iray,mline,ic,nchan[mline]/2-4,nchan[mline]/2+2);
/*
	printf("ray position %8.4f %8.4f \n",y1[recvd_iray],z1[recvd_iray]);
	fprintf(fpnotes,"ray position %8.4f %8.4f \n",y1[recvd_iray],z1[recvd_iray]);
*/
               fprintf(fpnotes,"YY ");
               printf("YY  ");
               for (i=nchan[mline]/2-4;i<nchan[mline]/2+3;i++){
                   fprintf(fpnotes,"%12.3e ",yspec[i]);
                   printf("%12.3e ",yspec[i]);
               }
           fprintf(fpnotes,"\n");
           printf("\n");
        }

/* Just a printing block */
        if (recvd_iray > 165 && recvd_iray < -170) {
      fprintf(fpnotes,"bJ index test C: ray line count nhyp %d %d %d %d\n",
         recvd_iray,mline,ic,nh);
           k=0;
           for (i=0;i<nh;i++){
           for (j=0;j<icount[recvd_iray];j++){
      fprintf(fpnotes,"ih ic k barJray %5d %5d %5d %12.3e \n", i,j,k,barJray[i][j]);
           k++;
           }}
        }

/*fprintf(fpnotes,"P%d packing ray %d for line %d from P%d\n",
my_rank,recvd_iray,mline,orig);
*/

/* Re-run the ray path function for a 2nd time to get the segments of the received ray */
		result = raypath_(&x1[recvd_iray],&y1[recvd_iray],&z1[recvd_iray],
		         &x2[recvd_iray],&y2[recvd_iray],&z2[recvd_iray],
                         &nbox,nx,ny,nz,cellsizex,cellsizeyz,
                         &maxgrid,t4,t5,t6,
                         &maxvox,t7,t8,&i);
                if (result == 0) {
			printf("Bad status from 2nd call to raypath\n");
			master_status = 0;
			current_ray = rays_to_compute;
			goto end_parallel;
		}
                if (i-1 != icount[recvd_iray]) {
			printf("Loop 3, not getting the same ray segments %d %d\n",i-1,icount[iray]);
			master_status = 0;
			current_ray = rays_to_compute;
			goto end_parallel;
                }
		if (icount[recvd_iray] == 0) {
			fprintf(fpnotes,"5 Bad icount from raypath %d %d\n",icount[recvd_iray],recvd_iray);
			printf("5 Bad icount from raypath %d %d\n",icount[recvd_iray],recvd_iray);
		}
                for (ivox=0;ivox<icount[recvd_iray];ivox++){
                        pathlength[ivox] = t8[ivox];
                for (i=0;i<4;i++){
                        k = 4*ivox + i;
                        iarc[i][ivox] = t7[k]-1;
                }}


/* Add the lambda factor and mean radiation into sums indexed on the model cube. */
if (icount[recvd_iray] <= 0 || nh <= 0) printf("lost icount nhyp %d %d\n",
icount[recvd_iray],nh);

                for (ivox=0;ivox<icount[recvd_iray];ivox++){
                        i = iarc[0][ivox];
                        j = iarc[1][ivox];
                        k = iarc[2][ivox];
                        ibox = iarc[3][ivox];    

           	chksum = 0.;
		for (ihyp=0;ihyp<nh;ihyp++){
     			apLJbar[ibox][i][j][k].apL[mline][ihyp]  +=  apLray[ihyp][ivox];
       			apLJbar[ibox][i][j][k].Jbar[mline][ihyp] += barJray[ihyp][ivox];
       			apLJbar[ibox][i][j][k].norm[mline][ihyp] += normray[ihyp][ivox];
           		chksum = chksum + barJray[ihyp][ivox];
/*	fprintf(fpnotes,"ijk1 %d %d %d line %d norm %e %e\n",
            i,j,k,mline,apLJbar[ibox][i][j][k].norm[mline][ihyp],normray[ihyp][ivox]); */
		}
	if (chksum == 0.) {
           fprintf(fpnotes,"Jbar to grid is zero for line=%d, iter=%d, i=%d,j=%d,k=%d\n",
                     mline,iter,i,j,k);
           printf("Jbar to grid is zero for line=%d, iter=%d, i=%d,j=%d,k=%d\n",
                     mline,iter,i,j,k);
	}
                } /* end of ivox loop */

/* ------------------------------------------------------- pacdata --------- */
	gettimeofday(&splittimeA, &zone);
/* pacdata puts the radiation into the output cube on the last iteration. */
/* If we are on the last iteration, then check the line number that we
   are working on, mline, to see if it is in the list of output lines.
   It should be because on the last iteration, we only compute the lines in
   the linelist. Nonetheless, loop through the linelist for a match. 
   Do the same procedure to see if the ray number that we just received 
   is in one of the output viewing angles. The rotation matrix is required
   to make each view look as if it were in the plane of the sky Y,Z looking
   down the X axis. The rotation matrix turns the coordinates to this position
   and pacdata writes them out in traditional imaging format with velocity
   as the third axis instead of X. The rotation matrix is the inverse of the
   one in fsub.f subroutine, ROTATE. More comments on the rotation in there.
*/
        if (iter == iterations) {
        for (iline=0;iline<nllist;iline++) {
            if ( mline == linelist[iline] ) {
/*  fprintf(fpnotes,"recvd iray %d\n",recvd_iray);  */
            for (kview=0;kview<outviews;kview++) {
                if (recvd_iray >= outray1[kview] && 
			recvd_iray < outray1[kview]+rpview[outAngleIndex[kview]]) {

			sinlong = sin(PI/180.*vlng[outAngleIndex[kview]]);
			coslong = cos(PI/180.*vlng[outAngleIndex[kview]]);
			sinlat  = sin(PI/180.*vlat[outAngleIndex[kview]]);
			coslat  = cos(PI/180.*vlat[outAngleIndex[kview]]);

			y = -sinlong*x1[recvd_iray] + coslong*y1[recvd_iray];
			z =  sinlat*coslong*x1[recvd_iray] 
				+ sinlat*sinlong*y1[recvd_iray] + coslat*z1[recvd_iray];

/*  fprintf(fpnotes,"output matched %d %d %d\n",iter,mline,recvd_iray);  */
            iret = pacdata(&kview,&iline,&nchan[mline],
                           &y, &z,
                           &cellsizeyz[raybox[recvd_iray]],yspec);

	       }      /*  if iray > <= block         */
            }         /*  for kview loop             */
	    }         /*  mline = linelist if block   */
        }             /*  for iline loop             */
        } /*  if iter = iterations block */
/* all this above is for pacdata */

/*
if(iter==iterations){
	gettimeofday(&splittimeB, &zone);
	elapsed = (splittimeB.tv_sec  - splittimeA.tv_sec ) * 1000000
                + (splittimeB.tv_usec - splittimeA.tv_usec) ;
	fprintf(fpnotes,"P%d: pacdata for line %d required  %d usecs\n",
    		my_rank,mline,elapsed);
	printf (        "P%d: pacdata for line %d required  %d usecs\n",
    		my_rank,mline,elapsed);
}
*/
/* ------------------------------------------------------- pacdata --------- */

/* That ray is completely finished (sent and received). */
		rays_finished++;
		free_processes++; 
		first_proc = mpi_status.MPI_SOURCE;

/*
if(iter==iterations){
	elapsed = (splittimeB.tv_sec  - splittimeC.tv_sec ) * 1000000
                + (splittimeB.tv_usec - splittimeC.tv_usec) ;
	fprintf(fpnotes,"P%d: slave loop for slave %d required  %d usecs\n",
    		my_rank,first_proc,elapsed);
	printf (        "P%d: slave loop for slave %d required  %d usecs\n",
    		my_rank,first_proc,elapsed);
	gettimeofday(&splittimeC, &zone);
}
*/


/*
		printf("P%d: finished %d rays, free processes %d\n",
			my_rank,rays_finished,free_processes);
*/


	    } /* end of while loop for rays finished*/

/* This prints out the ADDS pattern again. The l=1 prevents printing. */
	l = 1;
	if (l == 0) {
        for (ibox=0;ibox<nbox;ibox++){
fprintf(fpnotes,"New box   %d\n",ibox);
        for (i=0;i<nx[ibox];i++){
fprintf(fpnotes,"New layer %d\n",i);
        for (k=0;k<nz[ibox];k++){
        for (j=0;j<ny[ibox];j++){
		l = model[ibox][i][j][k].adds/2;
		if (l > 9) l = l;
		fprintf(fpnotes,"%2d",l);
	} fprintf(fpnotes,"\n");
	}}}
	}

if (iter != iterations){
j = 1;
chksum = 0;
for (i=0;i<nrays;i++) chksum += sendcheck[i];
if (chksum != 0) {
	printf("The send checklist was not completed \n");
	fprintf(fpnotes,"The send checklist was not completed \n");
	for (i=0;i<nrays;i++) fprintf(fpnotes,"i sendcheck %d %d\n",i,sendcheck[i]);
	j = 0;
}  /* else {
printf("Send ray check list completed OK!\n");
fprintf(fpnotes,"Send ray check list completed OK!\n");
} */
for (i=0;i<nrays;i++) chksum += recvcheck[i];
if (chksum != 0) {
	printf("The receive checklist was not completed \n");
	fprintf(fpnotes,"The receive checklist was not completed \n");
	for (i=0;i<nrays;i++) fprintf(fpnotes,"i recvcheck %d %d\n",i,recvcheck[i]);
	j = 0;
}  /* else {
printf("Recv ray check list completed OK!\n");
fprintf(fpnotes,"Recv ray check list completed OK!\n");
} */
if (j ==0) master_status = 0;
} /* This is the end of the if (iter != iterations) block for the chksum */

end_parallel:
if (master_status == 0) goto exitnow;

/* Finished looping over all the rays. Normalize the brightness in each voxel */

/*
   printf("Finished all rays for iteration %d line %d \n",iter,mline);
   printf(" Normalizing the mean radiation field and the Lambda parameter\n");
   fprintf(fpnotes,"Finished all rays for iteration %d line %d \n",iter,mline);
   fprintf(fpnotes,"Normalizing the mean radiation field and the Lambda parameter\n");
*/

        for (ibox=0;ibox<nbox;ibox++){
        for (i=0;i<nx[ibox];i++){
        for (j=0;j<ny[ibox];j++){
        for (k=0;k<nz[ibox];k++){
          if ( model[ibox][i][j][k].adds !=  0) {
	     for (ihyp=0;ihyp<nh;ihyp++){
		apLJbar[ibox][i][j][k].apL[mline][ihyp] /= model[ibox][i][j][k].adds;
		apLJbar[ibox][i][j][k].Jbar[mline][ihyp] /= model[ibox][i][j][k].adds;
		apLJbar[ibox][i][j][k].norm[mline][ihyp] /= model[ibox][i][j][k].adds;
/*	fprintf(fpnotes,"ijk2 %d %d %d line %d ihyp %d xn %e\n",
            i,j,k,mline,ihyp,apLJbar[ibox][i][j][k].norm[mline][ihyp]); */
	     }
          } else { 
	     for (ihyp=0;ihyp<nh;ihyp++){
		apLJbar[ibox][i][j][k].Jbar[mline][ihyp] = ycmb[mline];
	     }
          }
	}}}}

/*	fprintf(fpnotes,"ijk4 %d %d %d line %d ihyp %d xn %e\n",
            16,16,16,1,1,apLJbar[0][16][16][16].norm[1][1]); */

	gettimeofday(&splittime4, &zone);
	elapsed = (splittime4.tv_sec  - splittime1.tv_sec ) * 1000
                + (splittime4.tv_usec - splittime1.tv_usec) / 1000;
	fprintf(fpnotes,"P%d: RT for line %d required  %ld secs\n",
    		my_rank,mline,elapsed/1000);
	printf (        "P%d: RT for line %d required  %ld secs\n",
    		my_rank,mline,elapsed/1000);

/*
	printf ("P%d: Finished Radiative Transfer for line %d\n",
               my_rank,mline); 
	fprintf (fpnotes,"P%d: Finished Radiative Transfer for line %d\n",
                my_rank,mline); 
*/

		} /* end of for loop on lines, mline */

	gettimeofday(&splittime4, &zone);
	elapsed = (splittime4.tv_sec  - splittime2.tv_sec ) * 1000
                + (splittime4.tv_usec - splittime2.tv_usec) / 1000;
	fprintf(fpnotes,"P%d: Lambda iteration required  %ld sec\n",
    		my_rank,elapsed/1000);
	printf (        "P%d: Lambda iteration required  %ld sec\n",
    		my_rank,elapsed/1000);

/* End of radiative transfer calculation for this iteration. */
/* ---------------------------------------------------------------------- */
/* If this is the last iteration, write the output and close the output grids */

            if (iter == iterations) {
		  printf("Maximum iterations reached. Writing output. \n");
                  avg_output_cube(hanning);
                  close_output_cube(model_number,pathname);
                  printf("Output file written\n");
            }

/* What does this do? If true should exit the while loop anyway */
	if (iter == iterations) break; 
	iter++;
	} /* end of while loop over iterations */


/* Radiation has been stored in dcube and written to disk */
/* Finished with complete simulation for the current model */


/* ---------------------------------------------------------------------- */

/* data and model comparison section */

  if (comparison != NO) {

    chisqmin = chisq_function();

    
    printf("Minimum chisq %f \n",chisqmin);
    fprintf(fpnotes,"Minimum chisq %f \n",chisqmin);
    fprintf(fphist,"%12.6f\n",chisqmin);
    fflush(fphist);

    if (chisqmin < chisqbest) {
      best_model = model_number;
      chisqbest = chisqmin;
    }

    printf("Best chisq for model %d is %f\n",model_number,chisqmin);
    fprintf(fpnotes,"Best chisq for model %d is %f\n",model_number,chisqmin);

    printf("Best chisq so far %f for model number %d\n",chisqbest,best_model);
    fprintf(fpnotes,"Best chisq so far %f for model number %d\n",chisqbest,best_model);

    fflush(stdout);


  }

	gettimeofday(&splittime4, &zone);
	elapsed = (splittime4.tv_sec  - splittime3.tv_sec ) * 1000
                + (splittime4.tv_usec - splittime3.tv_usec) / 1000;
	fprintf(fpnotes,"P%d: Elapsed time to complete model %d was  %ld secs\n",
    		my_rank,model_number,elapsed/1000);
	printf (        "P%d: Elapsed time to complete model %d was  %ld secs\n",
    		my_rank,model_number,elapsed/1000);


/* finished comparing the model and the data */
/* Now go back and define a new model for the next comparison */
} /* end of loop over model numbers */

/* Finished with all the models. Clean up and exit */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

exitnow:

fflush(stdout);

if (my_rank == 0) {
	gettimeofday(&splittime4, &zone);
	elapsed = (splittime4.tv_sec  - starttime.tv_sec ) * 1000
                + (splittime4.tv_usec - starttime.tv_usec) / 1000;
	fprintf(fpnotes,"P%d: Total elapsed time  %f minutes\n",
    		my_rank,elapsed/1000./60.);
	printf (        "P%d: Total elapsed time  %f minutes\n",
    		my_rank,elapsed/1000./60.);
        fflush(fpnotes); 

        printf("The number of processes is %d\n",processes);

/* Call home the slaves for their reward. */
    for (dest = 1; dest < processes; dest++) {
	printf ("P%d: Sending final call for slave %d\n", my_rank,dest); 
      for (i=0;i<5;i++) comm_ints[i] = 0;
/* MPI_Isend is non-blocking. MPI_Send is blocking */
/*      MPI_Isend(comm_ints,n,MPI_LONG,dest,tag,MPI_COMM_WORLD,mpi_request);  */
      n = 5;
      tag = 0; 
      MPI_Send(comm_ints,n,MPI_LONG,dest,tag,MPI_COMM_WORLD); 
    }
	printf ("P%d: Sent final call for slaves\n", my_rank);

}


/* Shut down MPI. All processes must call MPI_Finalize */

    MPI_Finalize();
//	 if (my_rank == 0) printf("P%d: Process complete, exiting now ...\n",my_rank);
	 printf("P%d: Process complete, exiting now ...\n",my_rank);

	exit(0);  

}/* end of main */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

void getmodel_(int *ibox, int *i, int *j, int *k,
        float *temperature, float *density,
        float *vx, float *vy, float *vz, float *abund, float *lwidth,
        int *adds)
{   

/* model is global */
    
        *temperature    = model[*ibox][*i][*j][*k].temperature;
        *density        = model[*ibox][*i][*j][*k].density;
        *vx             = model[*ibox][*i][*j][*k].vx;
        *vy             = model[*ibox][*i][*j][*k].vy;
        *vz             = model[*ibox][*i][*j][*k].vz;
        *abund          = model[*ibox][*i][*j][*k].abundance;
        *lwidth         = model[*ibox][*i][*j][*k].linewidth;
        *adds           = model[*ibox][*i][*j][*k].adds;

}


void getalbj_(int *ibox, int *i, int *j, int *k, int *line, int *ihyp,
        float *apL, float *Jbar, float *norm)
{

/* apLJbar is global */

        *apL = apLJbar[*ibox][*i][*j][*k].apL[*line][*ihyp];
        *Jbar= apLJbar[*ibox][*i][*j][*k].Jbar[*line][*ihyp];
        *norm= apLJbar[*ibox][*i][*j][*k].norm[*line][*ihyp];

}

void putalbj_(int *ibox, int *i, int *j, int *k, int *line, int *ihyp,
        float *apL, float *Jbar, float *norm)
{

/* apLJbar is global */

/*        fprintf(fpnotes,"ijk line hyp %d %d %d %d %d %d a b %e %e\n",
               *ibox,*i,*j,*k,*line,*ihyp,*apL,*Jbar);
*/

/*if (*apL != 0.) printf("Bad apL line ihyp %d %d %e\n",*line,*ihyp,*apL);*/

        apLJbar[*ibox][*i][*j][*k].apL[*line][*ihyp]  = *apL;
        apLJbar[*ibox][*i][*j][*k].Jbar[*line][*ihyp] = *Jbar;
        apLJbar[*ibox][*i][*j][*k].norm[*line][*ihyp] = *norm;

}




int stateq(int nbox,int nx[] ,int ny[], int nz[],int iter,
           int accel, int molecule, int dust_continuum, 
	   int initial_radiation)
{

/* Values for the frac variables are transferred from F77. Hold these
   static to prevent memory problems. The other transferred variables
   are global. */
/*	static double frac_up[NLINES],frac_lo[NLINES], frac[NSTATE]; */
	static double frac[NSTATE]; 
	double db,aulhyp[NLINES][MAXHYP],bulhyp[NLINES][MAXHYP],
                  bluhyp[NLINES][MAXHYP];

	int i,j,k,ibox,line,q,ihyp,ii,jj,kk,it,u,l;
	double qfu,qfl,qju,qjl;
	double gff;
	int m,n,status,nh,select;
        float sumj,suma,sums,chksum,tk,h2,factor;
        float iontemp,a1,a2,plf,src;


/* For this iteration, calculate the level populations in the upper
   and lower states of each line, the source function, and the line center 
   opacity for all the voxels. 

   In this function, we calculate the level populations in each cell,
   but we don't save them. Instead, for each voxel, we go on to calculate 
   the source function and opacity and save those.

   For simple molecules, the upper state of a transition is the same as
   the lower state of the next higher transition. In these cases the
   level populations, frac_up[n] and frac_lo[n+1] are the same.
*/

	status = 1;

/* printf("Starting inside stateq\n"); */

/* Loop over all the cells in all the boxes. */

        for (ibox=0;ibox<nbox;ibox++){
        for (k=0;k<nz[ibox];k++){
/* printf("stateq box %d level %d \n",ibox,k); */
        for (j=0;j<ny[ibox];j++){
        for (i=0;i<nx[ibox];i++){

/* But if there were no rays through the cell, there is not much point. */
	if (model[ibox][i][j][k].adds == 0) goto skip; 

/* If the gas is ionized, do not need the molecular emissivity and
   opacity. Not much point, because the temperature should be 10^4 K. 
   The solution might go bad at very high temperature.
   Skip over the calculation of the molecular level pops
*/   
   
/* If the gas is ionized, zero the molecular emissivity and opacity assuming that
   there are no molecules in the ionized gas. Skip over the calculation of
   molecular level populations*/

        if (model[ibox][i][j][k].phase == 2) {
	    for (line=0;line<NLINES;line++){
	    for (ihyp=0;ihyp<atomic.nhyp[line];ihyp++){
                srcopc[ibox][i][j][k].emissivity[line][0] = 0.0;
                srcopc[ibox][i][j][k].opacity[line][0] = 0.0;
	    }}
	goto continuum; 
        }


/* Use these variables to pass down pointers to F77 functions so 
   that i,j,k do not come back altered by loop indices inside 
   the functions. Should be OK because the F77 functions have
   different names in the arg list, but just in case.
   iter is OK because it is not used inside
   stateq and is not passed back up (not a pointer).*/
        ii = i;
        jj = j;
        kk = k;

/*
c Given the temperature, density, and mean radiation J bar
c compute the level populations in each cell of the model. 
c If you want the optically thin limit, run LEVELPOPS 
c with BARJGRID set to zero. In the optically thin case,
c no need to iterate.
*/

/* Alternatively calculate the population in LTE. This LTE function
   is for linear rotors. */

/*	
        ltepops_(
                model[ibox][i][j][k].temperature,
                frac); 
*/

/*
c ACCELERATEDLAMBDA includes 
c the approximate lambda operator and the source function.
c If approxlambda is zero, then the routine calculates the level
c populations as for normal lambda iteration. 
c Should be ok on the first iteration because the approx lambda
c operator is zero.

*/

/*  printf("Starting grid cell ijk %d %d %d\n",i,j,k);  */

/* NH3      is molecule 2
   CH3CN    is molecule 3
   H2Opara  is molecule 5
   H2Oortho is molecule 6

   NH3 and CH3CN populations are calculated in LTE
   H2Oortho is calculated with the sub-critical approximation
   The select variable is used in the switch statement to
   select the method of calculating the level populations.
*/

	select = 0;
	if (molecule == 2 || molecule == 3) select  = 1;
// If the next line is commented, the code will use ALI for H2O
// active means subcritical approximation
//	if (molecule == 6) select = 2;   

	switch(select) {

/* NH3 and CH3CN are LTE only */
	case 1: {

/* printf("Starting LTE for molecule %d\n"); */

/* The level populations are in FRAC */

	if (molecule ==3) ltenh3_(&model[ibox][i][j][k].temperature,frac);
	if (molecule ==2) ltech3cn_(&model[ibox][i][j][k].temperature,frac);


//        if ( (j == 30 ||j == 31) && (k >= 30 && k <= 32))
//	fprintf(fpnotes,"ibox %3d i,j,k %3d %3d %3d pops %f %f %f %f %f %f %f\n",
//		ibox,i,j,k,frac[0],frac[1],frac[2],frac[3],frac[4],frac[5],
//		frac[0]+frac[1]+frac[2]+frac[3]+frac[4]+frac[5]);


	break;
		} 

/* This is a special case using the subcritical approximation for H2O */
	case 2: {

		status = subcritical_(&iter,&ii,&jj,&kk,
				&model[ibox][i][j][k].temperature,
				&model[ibox][i][j][k].dust_temperature,
				&model[ibox][i][j][k].density,
				&model[ibox][i][j][k].abundance,
				&model[ibox][i][j][k].linewidth,
				&model[ibox][i][j][k].mean_av,
				apLJbar[ibox][i][j][k].Jbar,
				frac);
		status = 1;
//		factor = fraction_para_(&model[ibox][i][j][k].temperature);
                factor = 0.;
        if ( (j == 31) && k == 31 )
	fprintf(fpnotes,"status %3d ibox %3d i,j,k %3d %3d %3d pops %f %f %f %f \n",
		status,ibox,i,j,k,frac[0],frac[1],
		frac[0]+frac[1],factor);

                if (status == 0) {
                    printf("Failed status check in subcritical\n");
                    fprintf(fpnotes,"Failed status check in subcritical\n");
                    return 0;
                }
	break;
		}

/* All the other molecules are non-LTE. The function acceleratedlambda
   handles both LTE and non-LTE hyperfines. */
	default: {

/* Copy the 2D arrays Jbar and apL, that got here though global.h
   to 1D vectors for transfer to F77 */

/*
fprintf(fpnotes,"Transfer bJ ApL emiss opac: ijk iter %3d %3d %3d %3d\n",i,j,k,iter);
*/

	q=0;
        for(line=0;line<NLINES;line++){
        chksum = 0.;
	nh = 1; if (nltehyp == 1) nh = atomic.nhyp[line];
/*	fprintf(fpnotes,"nh for stateq %d\n",nh); */
        for(ihyp=0;ihyp<nh;ihyp++){

/* 
if (i <= 16 && j == 16 && k == 16)
          fprintf(fpnotes,"ijk3 %3d %3d %3d l h bJ al em op %3d %3d %12.4e %12.4e %12.4e %12.4e\n",
            i,j,k,line,ihyp,
		apLJbar[ibox][i][j][k].Jbar[line][ihyp],
                apLJbar[ibox][i][j][k].apL[line][ihyp],
                srcopc[ibox][i][j][k].emissivity[line][ihyp],
		srcopc[ibox][i][j][k].opacity[line][ihyp]);
*/

		chksum = chksum + apLJbar[ibox][i][j][k].Jbar[line][ihyp];
		q++;
        }

/* initial_radiation is set in setup.c. Zero means that the initial radiation is zero
   The chksum is the sum of Jbar in all the cells in the model
*/
	if (chksum == 0. && initial_radiation != 0) {
           fprintf(fpnotes,"Zero Jbar for line=%d, iter=%d, i=%d,j=%d,k=%d\n",
                     line,iter,i,j,k);
           printf("Zero Jbar for line=%d, iter=%d, i=%d,j=%d,k=%d\n",
                     line,iter,i,j,k);
	}

	}

/*   fprintf(fpnotes,"A norm %d %d %d %e \n",
		i,j,k,apLJbar[ibox][i][j][k].norm[0][0]);  */

/*   printf("Starting acceleratedlambda\n");  */
 
	status = acceleratedlambda_(
		&model[ibox][i][j][k].temperature,
		&model[ibox][i][j][k].density,
		apLJbar[ibox][i][j][k].norm,
		apLJbar[ibox][i][j][k].Jbar,
           	apLJbar[ibox][i][j][k].apL,
           	srcopc[ibox][i][j][k].emissivity,
           	srcopc[ibox][i][j][k].opacity,
		frac,
           	&accel,&ii,&jj,&kk,&iter);

/*   fprintf(fpnotes,"B norm %d %d %d %e \n",
		i,j,k,apLJbar[ibox][i][j][k].norm[0][0]); */

/*   printf("Finished acceleratedlambda with status = %d\n",status);     */

	if (status == 0) return 0;

	break;
		} /* This is the end of CASE default block */
		} /* This is the end of the SWITCH statement */

/*
	fprintf(fpnotes,"bijk %d %d %d %d tdalu %f %g %g %f %f\n",
		ibox,i,j,k,
		model[ibox][i][j][k].temperature,
		model[ibox][i][j][k].density,
		model[ibox][i][j][k].abundance,
		frac_up[0],frac_lo[0]);
*/

/* This section checks the populations for negative values.
 */
	for (u=0;u<NSTATE;u++){
	   if (frac[u] < -1.e-10 ) break;
	}
	if (u<NSTATE){
	printf("Negative population %16.9e in line %3d at cell %3d %3d %3d %3d\n",
		frac[u],u,ibox,i,j,k);
	fprintf(fpnotes,"Negative population %16.9e in line %3d at cell %3d %3d %3d %3d\n",
		frac[u],u,ibox,i,j,k);
	} 

/* With the new level populations, calculate the source and opacity for the 
current voxel. Then store them in the model grid. */

/* printf("Starting sourceopac\n");  */
/* if (iter == 2) fprintf(fpnotes,"Starting sourceopac, 2nd iter\n"); */

/* Calculate the Einstein A,s and B,s for non-LTE hyperfine transitions */
/* There are 2 cases, N2H+ molecule=4 and HCN molecule=21. Each is handled
   separately */
        if (nltehyp == 1) {

/* fprintf(fpnotes,"i j k %d %d %d\n",i,j,k); */

        if (molecule == 4) {
	for (line=0;line<NLINES;line++){
	for (ihyp=0;ihyp<atomic.nhyp[line];ihyp++){

           it = atomic.start[line] + ihyp;
	   u = atomic.upper[it];
	   l = atomic.lower[it];
           qfu  = atomic.ln[u];
           qju  = atomic.jn[u];
           qfl  = atomic.ln[l];
           qjl  = atomic.jn[l];

           db = pow((C/atomic.freq[line]),2)
              / 2.0 /(PLANCK*atomic.freq[line]) ;

           aulhyp[line][ihyp] = atomic.aulhyp[line][ihyp] ;
           bulhyp[line][ihyp] = db*aulhyp[line][ihyp] ;
           bluhyp[line][ihyp] = db*aulhyp[line][ihyp]
                   * (2.0*qfu + 1.0) / (2.0*qfl + 1.0);

/* The factor hv/4pi is in comments as a reminder that we don,t
   include this factor in the variables opacity and emissivity.
   It gets factored in later in NEWRT of fsub.f for example. */

          srcopc[ibox][i][j][k].opacity[line][ihyp] =
               (  frac[l]*bluhyp[line][ihyp]
             -    frac[u]*bulhyp[line][ihyp]
               ); /* * PLANCK*atomic.freq[line]/(4.*PI); */

          srcopc[ibox][i][j][k].emissivity[line][ihyp] =
                  frac[u]*aulhyp[line][ihyp] 
		; /* * PLANCK*atomic.freq[line]/(4.*PI); */

	}}
	} /* end of if molecule == 4 */


	if (molecule == 21)
       {
	for (line=0;line<NLINES;line++){
	for (ihyp=0;ihyp<atomic.nhyp[line];ihyp++){

           it = atomic.start[line] + ihyp;
           u = atomic.upper[it];
           l = atomic.lower[it];

           db = pow((C/atomic.freq[line]),2)
              / 2.0 /(PLANCK*atomic.freq[line]) ;

           aulhyp[line][ihyp] = atomic.aulhyp[line][ihyp] ;
           bulhyp[line][ihyp] = db*aulhyp[line][ihyp] ;
           bluhyp[line][ihyp] = db*aulhyp[line][ihyp]
                   * atomic.statdg[u] / atomic.statdg[l];

          srcopc[ibox][i][j][k].opacity[line][ihyp] =
               (  frac[l]*bluhyp[line][ihyp]
             -    frac[u]*bulhyp[line][ihyp]
               ); /* * PLANCK*atomic.freq[line]/(4.*PI); */

          srcopc[ibox][i][j][k].emissivity[line][ihyp] =
                  frac[u]*aulhyp[line][ihyp] 
		; /* * PLANCK*atomic.freq[line]/(4.*PI); */

	}}
       }
     
	} else {

/* This section calculates the Einstein A,s and B,s for molecules with
   no hyperfine structure or for molecules with the LTE-hyperfine approx. 

   This works because for both cases because we are using the values of 
   Einstein A and statistical weights for the total rotational 
   transition (not the hyperfine transitions).
*/
	for (line=0;line<NLINES;line++){

 	    u = atomic.indexu[line];
	    l = atomic.indexl[line];

           db = pow((C/atomic.freq[line]),2)
              / 2.0 /(PLANCK*atomic.freq[line]) ;

           aulhyp[line][0] = atomic.einsta[line];
           bulhyp[line][0] = db*aulhyp[line][0];
           bluhyp[line][0] = db*aulhyp[line][0]
		   * atomic.statdg[u] /atomic.statdg[l] ;


/*
	fprintf(fpnotes,"U L %d %d Gu Gl %e %e db %e bul blu %e %e %e\n",
		atomic.indexu[line],atomic.indexl[line],
		atomic.statdg[u],atomic.statdg[l],db,
                bulhyp[line][0],bluhyp[line][0],bluhyp[line][0]/bulhyp[line][0]);

	printf("U L %d %d Gu Gl %e %e db %e bul blu %e %e %e\n",
		atomic.indexu[line],atomic.indexl[line],
		atomic.statdg[u],atomic.statdg[l],db,
                bulhyp[line][0],bluhyp[line][0],bluhyp[line][0]/bulhyp[line][0]);
*/

          srcopc[ibox][i][j][k].opacity[line][0] =
               (  frac[l]*bluhyp[line][0]
             -    frac[u]*bulhyp[line][0]
               ); /* * PLANCK*atomic.freq[line]/(4.*PI); */

          srcopc[ibox][i][j][k].emissivity[line][0] =
                  frac[u]*aulhyp[line][0] 
		; /* * PLANCK*atomic.freq[line]/(4.*PI); */

plf = 2.*PLANCK*atomic.freq[line]/C*atomic.freq[line]/C*atomic.freq[line] 
/ (exp(PLANCK*atomic.freq[line]/(BOLTZ*model[ibox][i][j][k].temperature)) - 1.);
src = srcopc[ibox][i][j][k].emissivity[line][0]/srcopc[ibox][i][j][k].opacity[line][0];
//if (src > plf) fprintf(fpnotes,"src > planck \n"); 

/*
 if ((srcopc[ibox][i][j][k].opacity[line][0] < 0.0) || (j == 16 && k == 16)) { 
ihyp = 0;
fprintf(fpnotes,"opacities %3d %3d %3d %3d %3d %4.1f %4.1f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
i,j,k,
line,ihyp,atomic.statdg[u] ,atomic.statdg[l],
frac[l],bluhyp[line][ihyp],
frac[u],bulhyp[line][ihyp],
aulhyp[line][ihyp],
srcopc[ibox][i][j][k].emissivity[line][ihyp],
srcopc[ibox][i][j][k].opacity[line][ihyp]);
}
*/

        }
     }


/* In addition to the factor of hnu/4pi above, the emissivities and opacities
   need to take account of the relative abundance of the molecule. Include
   the abundance here. This applies to all molecules with or without hyperfines.
*/

/*
if (j == ny[ibox]/2 && k == nz[ibox]/2)
        fprintf(fpnotes,"bijk %d %d %d %d TDJEOFF %f %g %g %e %e %f %f %e %e\n",
                ibox,i,j,k,
                model[ibox][i][j][k].temperature,
                model[ibox][i][j][k].density,
                apLJbar[ibox][i][j][k].Jbar[0][0],
                srcopc[ibox][i][j][k].emissivity[0][0],
                srcopc[ibox][i][j][k].opacity[0][0] ,
                frac[0],frac[1],
                srcopc[ibox][i][j][k].emissivity[0][0]/srcopc[ibox][i][j][k].opacity[0][0],
                plf);
*/
//                frac[1]*aulhyp[0][0]/(frac[0]*bluhyp[0][0] - frac[1]*bulhyp[0][0]));
//   / (2.*BOLTZ) * pow((C/atomic.freq[0]),2)    This is the factor to convert the source function to K


        for (line=0;line<NLINES;line++){
	nh = 1; if (nltehyp == 1) nh = atomic.nhyp[line];
        for (ihyp=0;ihyp<nh;ihyp++){ 
          srcopc[ibox][i][j][k].emissivity[line][ihyp] *= model[ibox][i][j][k].abundance;
          srcopc[ibox][i][j][k].opacity[line][ihyp]    *= model[ibox][i][j][k].abundance;
        }
            if (srcopc[ibox][i][j][k].opacity[line][0] < 0.0) {
	      fprintf(fpnotes,"line ihyp %d %d emis opac abund %12.4e %12.4e %12.4e\n",line,ihyp,
	      srcopc[ibox][i][j][k].emissivity[line][0],
	      srcopc[ibox][i][j][k].opacity[line][0], model[ibox][i][j][k].abundance); 
            }
	}



/*  printf("Finished calculating emissivity and opacity for the molecular line\n"); */
/*  Now calculate the continuum emissivity and opacity */

continuum:


/* Initialize the continuum opacity and emissivity to zero for summing. 
   The final continuum is the sum of the dust and ionized */
        for (line=0;line<NLINES;line++){
              srcopc[ibox][i][j][k].cont_opacity[line] = 0.0;
              srcopc[ibox][i][j][k].cont_emissivity[line] = 0.0;
        }

/* If there is dust add its continuum radiation */
      if (dust_continuum == 1) {

/*
c determine the continuum opacity for each line frequency,
c unit is cm² per H_2 molecule. This needs to be multiplied by
c the number density of H_2 to get to cm-1 or the column
c density to become non-dimensional optical depth. In NEWRT
c the opacity is multiplied by the column density.
c This parameterization is from Zucconi et al. A&A 376, 650 (2001)

The factor of 4 is from the model of KC08 where we assume 4x the normal
dust opacity to account for fluffy dust.
*/

        for (line=0;line<NLINES;line++){

          if (atomic.freq[line] > 3.e13) {
              srcopc[ibox][i][j][k].cont_opacity[line]
                  += 4.0*3.9e-22 * pow((atomic.freq[line]/3.e14),1.4);
          } else {
              if (atomic.freq[line] > 7.5e11) {
                   srcopc[ibox][i][j][k].cont_opacity[line]
                      += 4.0*1.5e-24 * pow((atomic.freq[line]/2.142857e12),1.6);
              } else {
                   srcopc[ibox][i][j][k].cont_opacity[line]
                      += 4.0*3.3e-26 * pow((atomic.freq[line]/2.83e11),2.0);
              }
          }
/* Finished with the opacity. Find the emissivity from Kirchoff's law
   In radiative equilibrium, the dust emissivity/opacity = black body
*/
          srcopc[ibox][i][j][k].cont_emissivity[line]
            += 2.0*PLANCK * pow((atomic.freq[line]/C),2) * atomic.freq[line]
            / ( exp(
                PLANCK*atomic.freq[line]
                /(BOLTZ*model[ibox][i][j][k].dust_temperature)
                   ) -1.0
              )
            * srcopc[ibox][i][j][k].cont_opacity[line] ;

        } /* End of loop over lines */

/* Now divide out the factor of hv/4pi because the emissivity and
   opacity are multiplied by this factor in NEWRT */
        for (line=0;line<NLINES;line++){
          factor = 1.0;
          factor = PLANCK*atomic.freq[line]/(4.*PI);
	  srcopc[ibox][i][j][k].cont_emissivity[line] /= factor;
	  srcopc[ibox][i][j][k].cont_opacity[line] /= factor;
/*
if (jj == 31 && kk == 31) {
fprintf(fpnotes,"cont line %12.4e %12.4e %12.4e %12.4e \n",
srcopc[ibox][i][j][k].cont_emissivity[line],
srcopc[ibox][i][j][k].cont_opacity[line],
srcopc[ibox][i][j][k].emissivity[line][0],
srcopc[ibox][i][j][k].opacity[line][0]);
}
*/
	}

      } /* End of dust continuum */

/* If there is ionized gas add its continuum radiation */
        if (model[ibox][i][j][k].phase == 2) {

            iontemp = model[ibox][i][j][k].temperature;

            for (line=0;line<NLINES;line++){

/* In NEWRT, the opacities are multiplied by this factor, so 
   divide it out here. Also take it out of the emissivity so
   that when we form the source function as jv/kv we get the
   right answer. */
                factor = PLANCK*atomic.freq[line]/(4.*PI);

                gff = 9.77 * ( 1.0 +
                    0.13*log10( pow(iontemp,1.5)
                        /atomic.freq[line] ) );

/* The dust emissivity and opacity contain ne^2. Here there is a single
   factor of density because in NEWRT we multiply by the column density
   to get the optical depth.
*/
                srcopc[ibox][i][j][k].cont_emissivity[line] +=
                        5.44e-39 /factor * gff / sqrt(iontemp)
                        * exp( -PLANCK*atomic.freq[line]
                        /(BOLTZ*iontemp) )
                        * model[ibox][i][j][k].density;

                srcopc[ibox][i][j][k].cont_opacity[line] +=
                        0.1731*gff/9.77 /factor
                                * pow(iontemp,-1.5)
                                * pow(atomic.freq[line],-2)
                                * model[ibox][i][j][k].density;
/* In this print, multiply the emis and opac by factor and the density to put 
   them in cgs units. Compare opac with the approx formulae below.
   
        fprintf(fpnotes,"continuum %d %d %d  gff emis opac source %12.3e %12.3e %12.3e %12.3e\n",
                i,j,k,gff,
                srcopc[ibox][i][j][k].cont_emissivity[line]*model[ibox][i][j][k].density*factor,
                srcopc[ibox][i][j][k].cont_opacity[line]*model[ibox][i][j][k].density*factor,
                srcopc[ibox][i][j][k].cont_emissivity[line] / srcopc[ibox][i][j][k].cont_opacity[line]
                *pow(3.e10/atomic.freq[line],2)/(2.0*BOLTZ) );

                a1 = 3.014e-2 * pow(iontemp,-1.5) * pow((atomic.freq[line]/1.e9),-2)
                * (
                        log( 4.955e-2 / (atomic.freq[line]/1.e9) )
                      + 1.5*log(iontemp)
                  )
                * model[ibox][i][j][k].density*model[ibox][i][j][k].density /PC ;

                a2 = 8.235e-2 * pow(iontemp,-1.35) * pow((atomic.freq[line]/1.e9),-2.1) 
                * model[ibox][i][j][k].density*model[ibox][i][j][k].density /PC;

                fprintf(fpnotes,"approx opacity %12.3e, %12.3e %9.3f \n",a1,a2,a1/a2);
*/


             } /* End of loop over lines */
        } /* End of if phase == 2 */

/*
        for (line=0;line<NLINES;line++){
	printf("Lijk %3d %3d %3d %3d dust emis opac %12.3e %12.3e\n",
	line,i,j,k,
	srcopc[ibox][i][j][k].cont_emissivity[line],
	srcopc[ibox][i][j][k].cont_opacity[line]);
	}
*/

/* 	printf("Finished sourceopac with status = %d\n",status);  */

/*
	fprintf(fpnotes,"bijk %d %d %d %d dluaso %e %f %f %e %f %e\n",
		ibox,i,j,k,
		model[ibox][i][j][k].density,
		frac_up[0],frac_lo[0], 
		model[ibox][i][j][k].abundance,
           	srcopc[ibox][i][j][k].source[0],
           	srcopc[ibox][i][j][k].opacity[0]);
*/

		
skip:

ii = ii;  /* A do nothing statement to avoid a compiler warning */


	}}}}

	return 1;

} /* end of stateq */

int convergence(int nbox,int nx[] ,int ny[], int nz[], int iter, int *converged)
{

        int i,j,k,ibox,line,ihyp;
        int m,npt1,npt2,nh;
	static double srcold[MAXVOX][NLINES];
	double diff,avgdiff,diffmax,pctmax,avgpct,src[MAXVOX][NLINES];
	int linemax;


	ibox = 0;
//	printf("Convergence for outer box %d only of %d\n",ibox,nbox);

/* Pick a radius */
	j = ny[ibox]/2;
	k = nz[ibox]/2;

/* Highest line to consider for convergence */
	linemax = 3;
	if (linemax > NLINES) linemax = NLINES;

/*	printf("iteration number in convergence %d\n",iter); */

	*converged = 0;

	avgdiff = 0.;
	diffmax = 0.;
	avgpct = 0.;
        pctmax = 0.;
	npt1 = 0;
	npt2 = 0;

/* Define an "average" source function by summing over the hyperfines */
        for (i=0;i<nx[ibox];i++){
	for (line=0;line<linemax;line++){
            src[i][line] = 0.;
	    nh = 1; if (nltehyp == 1) nh = atomic.nhyp[line];
            for (ihyp=0;ihyp<nh;ihyp++) {
		if (srcopc[ibox][i][j][k].opacity[line][ihyp] != 0.0) {
                  src[i][line] = src[i][line] 
                    + srcopc[ibox][i][j][k].emissivity[line][ihyp]
                    / srcopc[ibox][i][j][k].opacity[line][ihyp];
		} 
/*
printf("hyperfine source %d %d %e %e\n",i,j,k,srcopc[ibox][i][j][k].emissivity[line][ihyp],
srcopc[ibox][i][j][k].opacity[line][ihyp]);
*/
            }
	    src[i][line] = src[i][line]/nh
                           / (2.*BOLTZ) * pow((C/atomic.freq[line]),2);
/*printf("average source %d %d %e\n",i,line,src[i][line]);*/
	}}


/* If we have an earlier iteration, then compare the current source */
	if (iter > 1) {

        for (i=0;i<nx[ibox];i++){
	for (line=0;line<linemax;line++){
	    diff = fabs(srcold[i][line] - src[i][line]);
	    avgdiff += diff;
	npt1++;
            if (diff > diffmax) diffmax = diff;
            if (srcold[i][line] != 0.) {
		avgpct  += diff/srcold[i][line];
		npt2++;
               	if (diff/srcold[i][line] > pctmax ) 
		   pctmax  = diff/srcold[i][line];
	    }
	}}
	if (npt1 != 0) avgdiff /= npt1;
	if (npt2 != 0) avgpct /= npt2;

	printf("Convergence: Absolute  :   average %10.6g maximum %10.6g\n",
            avgdiff,diffmax);
	printf("Convergence: Percentage:   average %10.6g maximum %10.6g\n",
            avgpct,pctmax);
	fprintf(fpnotes,"Convergence: Absolute:   average %10.6g maximum %10.6g\n",
		avgdiff,diffmax);
	fprintf(fpnotes,"Convergence: Percentage: average %10.6g maximum %10.6g\n",
		avgpct,pctmax);
	if (pctmax < 0.001 && pctmax > 0.) *converged = 1;
	
	fprintf(fpconv,"%d %g %g %g %g\n",iter,avgdiff,diffmax,avgpct,pctmax);
	fprintf(fpconv,"%d %d\n",nx[ibox],linemax);
        for (i=0;i<nx[ibox];i++){
	for (line=0;line<linemax;line++){
		fprintf(fpconv,"%f ",src[i][line]);
	} fprintf(fpconv,"\n");
	}
	fflush(fpconv);

	} /* This is the end of the if (iter > 1) block */

        for (i=0;i<nx[ibox];i++){
	for (line=0;line<NLINES;line++) {
		srcold[i][line] = src[i][line]; 
	}}


  return 1;

} /* end of function convergence */



void setupOutputCube(int nllist,int *linelist,int *nchan,int noutx,int nouty,int outviews,
float beamx,float  beamy,float cellx,float celly, float lng[], float lat[],
char **linenames,
float chvel[NLINES][MAXCH])

{
  int i,j,k,m,n;
  float sizex,sizey;

/* Set up the output data cubes */

/* dcube is global */
/* Assign values, allocate memory, initialize to zero */


 /*
        printf("nllist,noutx,nouty,beamx,beamy,cellx,celly,outviews %d %d %d %f %f %f %f %d\n",
        nllist,noutx,nouty,beamx,beamy,cellx,celly,outviews);
*/

/* All the axes and coordinates in this routine refer to 
   the output image axes */
	dcube.nlines = nllist;
	dcube.nx = noutx;
	dcube.ny = nouty;
	dcube.beamx = beamx;
	dcube.beamy = beamy;
	dcube.cellx = cellx;
	dcube.celly = celly;
	dcube.nviews= outviews;
	sizex = dcube.cellx*dcube.nx;
	sizey = dcube.celly*dcube.ny;


printf("Start allocating output cube %d %d %d %d\n",
	dcube.nlines,dcube.nviews,dcube.nx,dcube.ny);

          dcube.linenames = (char **) malloc(nllist*sizeof(char *));
          dcube.chvel = (float **)malloc(nllist*sizeof(float *));
          dcube.restfreq = (double *)malloc(nllist*sizeof(double *));
          for (i=0;i<nllist;i++){
             dcube.linenames[i] = (char *) malloc(MAXNAME*sizeof(char));
             dcube.chvel[i] = (float *)malloc(nchan[linelist[i]]*sizeof(float));
          }

printf("OK on linenames and ch vels\n"); 

          dcube.cube = (float *****) malloc(nllist*sizeof(float****));
          for (i=0;i<nllist;i++){
            dcube.cube[i] = (float ****) malloc(dcube.nviews*sizeof(float***));
          for (j=0;j<outviews;j++){
            dcube.cube[i][j] = 
              (float ***) malloc(dcube.nx*sizeof(float**));
            for (k=0;k<dcube.nx;k++){
              dcube.cube[i][j][k] = 
              (float **) malloc(dcube.ny*sizeof(float*));
              for (n=0;n<dcube.ny;n++){
                dcube.cube[i][j][k][n] = 
                 (float *) malloc(nchan[linelist[i]]*sizeof(float));
        }}}}

        
printf("OK on cube\n");
 
          dcube.wt = (float ***) malloc(outviews*sizeof(float**));
          for (j=0;j<outviews;j++){
            dcube.wt[j] = 
              (float **) malloc(dcube.nx*sizeof(float*));
            for (k=0;k<dcube.nx;k++){
              dcube.wt[j][k] = 
              (float *) malloc(dcube.ny*sizeof(float));
        }}

printf("OK on wts\n"); 
 

        dcube.nchan = (int *)malloc(nllist*sizeof(int));
        dcube.xplane = (float *)malloc((dcube.nx+1)*sizeof(float));
        dcube.yplane = (float *)malloc((dcube.ny+1)*sizeof(float));
        dcube.xcenter = (float *)malloc(dcube.nx*sizeof(float));
        dcube.ycenter = (float *)malloc(dcube.ny*sizeof(float));
        dcube.lng = (float *)malloc(dcube.nviews*sizeof(float));
        dcube.lat = (float *)malloc(dcube.nviews*sizeof(float));


	fprintf(fpnotes,"The output cell size is %f %f\n",dcube.cellx,dcube.celly);
	fprintf(fpnotes,"The grid size is        %d %d\n",dcube.nx,dcube.ny);
	fprintf(fpnotes,"The output beam size is %f %f\n",dcube.beamx,dcube.beamy);

	printf("The output cell size is %f %f\n",dcube.cellx,dcube.celly);
	printf("The grid size is        %d %d\n",dcube.nx,dcube.ny);
	printf("The output beam size is %f %f\n",dcube.beamx,dcube.beamy);

        for (i=0;i<nllist;i++){
	  dcube.nchan[i] = nchan[linelist[i]];
        }

	fprintf(fpnotes,"Number of lines in output list %d\n",nllist);

	printf("Number of lines in output list %d\n",nllist);

        for (i=0;i<nllist;i++){
        fprintf(fpnotes,"Line %d with %d channels\n",i,nchan[linelist[i]]);
        for (n=0;n<dcube.nchan[i];n++) {
          dcube.chvel[i][n] = chvel[linelist[i]][n];
/*
          fprintf(fpnotes,"%d %f\n",n,1.e-5*dcube.chvel[i][n]);  
*/
        }}


        for (i=0;i<nllist;i++){
          dcube.linenames[i] = linenames[linelist[i]];
          dcube.restfreq[i] = atomic.freq[linelist[i]];
/*	  fprintf(fpnotes,"output line names: %s\n",dcube.linenames[i]); */
        }

        for (j=0;j<dcube.nviews;j++) {
          dcube.lng[j] = lng[j];
          dcube.lat[j] = lat[j];
        }
	
        for (j=0;j<dcube.nx+1;j++)
          dcube.xplane[j] = -0.5*sizex + j*dcube.cellx;
        for (k=0;k<dcube.ny+1;k++)
          dcube.yplane[k] = -0.5*sizey + k*dcube.celly;

/*
        fprintf(fpnotes,"Here are the output grid planes\n");
        for (j=0;j<dcube.nx+1;j++)
          fprintf(fpnotes,"x: %d %9.6f \n",
            j,dcube.xplane[j]);
        for (k=0;k<dcube.ny+1;k++)
          fprintf(fpnotes,"y: %d %9.6f \n",
            k,dcube.yplane[k]);
*/



        for (j=0;j<dcube.nx;j++)
          dcube.xcenter[j] = -0.5*(sizex - dcube.cellx)
            + j*dcube.cellx;
        for (k=0;k<dcube.ny;k++)
          dcube.ycenter[k] = -0.5*(sizey - dcube.celly)
            + k*dcube.celly;


/*
        fprintf(fpnotes,"Here are the output box centers\n");
        for (j=0;j<dcube.nx;j++)
          fprintf(fpnotes,"x: %d %9.6f  \n",
            j,dcube.xcenter[j]);
        for (k=0;k<dcube.ny;k++)
          fprintf(fpnotes,"y: %d %9.6f  \n",
            k,dcube.ycenter[k]);
*/

	  for (i=0;i<dcube.nviews;i++) {
	  for (j=0;j<dcube.nx;j++) {
	  for (k=0;k<dcube.ny;k++) {
	    dcube.wt[i][j][k] = 0.;
	  for (m=0;m<dcube.nlines;m++) {
	  for (n=0;n<dcube.nchan[m];n++) {
	    dcube.cube[m][i][j][k][n] = 0.;
	  }}}}}
/*
	printf("finished cube setup\n");
*/
}


int avg_output_cube(int hanning)
{

  int i,j,ii,jj,k,n,ih,line,v,view,wt,i1,i2,j1,j2,inverse,status;
  int nx,ny;
  float *dtmp;
  int ntmp[40];
  float **beamR,    **beamI;
  float **beamRfft, **beamIfft;
  float **mapR,     **mapI;
  float **mapRfft,  **mapIfft;
  float **convR,    **convI;
  float **convRfft, **convIfft;
  float dx,dy, sx,sy,wt1,xcent,ycent,sumM,sumC,area;
  float dcubemax;

printf("Start smoothing the output data\n");
fprintf(fpnotes,"Start smoothing the output data\n");

/* Use the stored weights to normalize the spectra */
/* In pacdata, when we summed into dcube.wt, we also summed over every
   line. The cells are the same for each line so
   dcube.wt is now larger by a factor of nlines than it should be 
   for each single line. We need to divide dcube.wt by nlines.
   Then we divide the brightness in the output cell by this
   normalization which now represents the sum of all the weights
   of the rays that went into the cell.
*/

fprintf(fpnotes,"avg_output nx ny %d %d\n",dcube.nx,dcube.ny);
    for (v=0;v<dcube.nviews;v++) {
    for (i=0;i<dcube.nx;i++) {
    for (j=0;j<dcube.ny;j++) {
      dcube.wt[v][i][j] /= dcube.nlines;
/*fprintf(fpnotes,"i j wt %d %d %12.3e \n",i,j,dcube.wt[v][i][j]); */
      if (dcube.wt[v][i][j] > 0.) { 
        for (line=0;line<dcube.nlines;line++){
dcubemax = -1.e20;
        for (n=0;n<dcube.nchan[line];n++){
          dcube.cube[line][v][i][j][n] /= dcube.wt[v][i][j];
	  if (dcube.cube[line][v][i][j][n] > dcubemax) dcubemax = dcube.cube[line][v][i][j][n];
        }
/*
fprintf(fpnotes,"line view %3d %3d  i j  %3d %3d weight %12.3e dcube max %12.3e \n",
    line,v,i,j,dcube.wt[v][i][j],dcubemax);
*/
      }}
    }}}



  for (line=0;line<dcube.nlines;line++){

  dtmp = (float *)malloc(dcube.nchan[line]*sizeof(float));

  if (hanning == 1) {
      printf("Hanning smoothing the spectra\n");
      fprintf(fpnotes,"Hanning smoothing the spectra\n");
  }


/*  
    for (i=0;i<dcube.nx;i++) {
    for (j=0;j<dcube.ny;j++) {
fprintf(fpnotes,"avg %d %d %f \n",i,j,dcube.cube[0][0][i][j][dcube.nchan[0]/2]);
    }}
*/    


  for (ih=0;ih<hanning;ih++) {
    for (i=0;i<dcube.nviews;i++) {
    for (j=0;j<dcube.nx;j++) {
    for (k=0;k<dcube.ny;k++) {
        for (n=1;n<dcube.nchan[line]-1;n++) {
          dtmp[n] = 0.25*dcube.cube[line][i][j][k][n-1] 
                      + 0.50*dcube.cube[line][i][j][k][n] 
                      + 0.25*dcube.cube[line][i][j][k][n+1] ;
        }
        for (n=1;n<dcube.nchan[line]-1;n+=2) {
          dcube.cube[line][i][j][k][n/2] = dtmp[n];
        }
    }}}
    for (n=1;n<dcube.nchan[line]-1;n+=2) {
      dcube.chvel[line][n/2] = dcube.chvel[line][n];
    }
    printf("Reduced number of channels from %d ",dcube.nchan[line]);
    fprintf(fpnotes,"Reduced number of channels from %d ",dcube.nchan[line]);

    dcube.nchan[line] = (dcube.nchan[line] - 3)/2 + 1;

    printf("  to %d\n",dcube.nchan[line]);
    fprintf(fpnotes,"  to %d\n",dcube.nchan[line]);
  }
/*
  fprintf(fpnotes,"New velocity grid\n");
  for (n=0;n<dcube.nchan[line];n++) {
    fprintf(fpnotes,"%d %f\n",n,dcube.chvel[line][n]/1.e5);
  }
*/




/*
  for (j=0;j<dcube.nx;j++) {
  for (k=0;k<dcube.ny;k++) {
      fprintf(fpnotes,"i j wt %d %d %f\n",j,k,dcube.wt[0][j][k]/dcube.nlines);
  }}
*/

  free (dtmp);
  }

if (dcube.beamx > 0. || dcube.beamy > 0.) {

    nx = dcube.nx;
    ny = dcube.ny;

    fprintf(fpnotes,"Convolve the output by the model observing beam\n");
    fprintf(fpnotes,"X: beam %f cell %f\n",dcube.beamx,dcube.cellx);
    fprintf(fpnotes,"Y: beam %f cell %f\n",dcube.beamy,dcube.celly);

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

    printf("Convolve the output by the model observing beam\n");
    printf("X: beam %f cell %f\n",dcube.beamx,dcube.cellx);
    printf("Y: beam %f cell %f\n",dcube.beamy,dcube.celly);

/* This is the center of the model grid */
     xcent = dcube.xcenter[nx/2] - 0.5*dcube.cellx;
     ycent = dcube.ycenter[ny/2] - 0.5*dcube.celly;
/* Square the beam */
            sx = dcube.beamx*dcube.beamx;
            sy = dcube.beamy*dcube.beamy;

	area = PI * dcube.beamx*dcube.beamy / (dcube.cellx*dcube.celly);

     for (i=0;i<dcube.nx;i++) {
     for (j=0;j<dcube.ny;j++) {
	    dx = dcube.xcenter[i] - xcent;
	    dy = dcube.ycenter[j] - ycent;
            dx = dx * dx;
            dy = dy * dy;
            beamR[i][j] = exp( -dx/sx - dy/sy );
            beamI[i][j] = 0.0;
     }}

/*
    j = 16;
    for (i=0;i<nx;i++) printf(" beamR %d %d %f\n",i,j,beamR[i][j]);
    printf("Beam defined\n");
*/
     inverse = 0;
     status = fft2d(inverse, nx, ny, beamR, beamI, beamRfft, beamIfft);

    printf("Beam FFT\n");

/*
    j = 16;
    for (i=0;i<nx;i++) printf(" beamRfft %d %d %f\n",i,j,beamRfft[i][j]);
*/

for (line=0;line<dcube.nlines;line++) {
for (view=0;view<dcube.nviews;view++) {
for (n=0;n<dcube.nchan[line];n++) {

     for (i=0;i<dcube.nx;i++) {
     for (j=0;j<dcube.ny;j++) {
         mapR[i][j] = dcube.cube[line][view][i][j][n];
         mapI[i][j] = 0.0;
     }}

     inverse = 0;
     status = fft2d(inverse, nx, ny, mapR, mapI, mapRfft, mapIfft);

     for (i=0;i<dcube.nx;i++) {
     for (j=0;j<dcube.ny;j++) {
	convRfft[i][j] = mapRfft[i][j]*beamRfft[i][j];
	convIfft[i][j] = mapIfft[i][j]*beamRfft[i][j];
/*
	convRfft[i][j] = mapRfft[i][j]*beamRfft[i][j] - mapIfft[i][j]*beamIfft[i][j];
	convIfft[i][j] = mapRfft[i][j]*beamIfft[i][j] + mapIfft[i][j]*beamRfft[i][j];
*/
     }}

     inverse = 1;
     status = fft2d(inverse, nx, ny, convR, convI, convRfft, convIfft);

     for (i=0;i<dcube.nx;i++) {
     for (j=0;j<dcube.ny;j++) {
	dcube.cube[line][view][i][j][n] = convR[i][j]/area;
     }}


    }}} /* end of for line view channel */
 
  free(beamR);
  free(beamI);
  free(beamRfft);
  free(beamIfft);
  free(convR);
  free(convI);
  free(convRfft);
  free(convIfft);
  free(mapR);
  free(mapI);
  free(mapRfft);
  free(mapIfft);


} else {

    fprintf(fpnotes,"No convolution: because beamsizes < cellsizes\n");
    fprintf(fpnotes,"X: beam %f cell %f\n",dcube.beamx,dcube.cellx);
    fprintf(fpnotes,"Y: beam %f cell %f\n",dcube.beamy,dcube.celly);

}/* end of if block for convolution */

  return 1;
}


int close_output_cube(int model_number, char pathname[])

{

  char fname[64],  suffix[4];
  FILE *fpout;
  int i,j,k,l,m,n;
  long int nobj,nbytes;
  int ipout;

        strcpy(fname,pathname);
        strcat(fname,"ModelCube");
        sprintf(suffix,"%d",model_number);
        strcat(fname,suffix);
        fpout = fopen(fname,"w");
             if (fpout == NULL) {
                printf("Could not open %s. Exiting ...\n",fname);
                exit(0);
             } else {
                printf("File %s opened\n",fname);
             }




/* write dcube structure */

      nobj = nbytes = 0;
      nobj += fwrite(&dcube.nlines,sizeof(dcube.nlines),1,fpout);
      nbytes += sizeof(dcube.nlines);
      nobj += fwrite(&dcube.nviews,sizeof(dcube.nviews),1,fpout);
      nbytes += sizeof(dcube.nviews);
      for (m=0;m<dcube.nlines;m++) {
        nobj += fwrite(&dcube.nchan[m],sizeof(dcube.nchan[m]),1,fpout);
        nbytes += sizeof(dcube.nchan[m]);
      }
      nobj += fwrite(&dcube.nx,sizeof(dcube.nx),1,fpout);
      nbytes += sizeof(dcube.nx );
      nobj += fwrite(&dcube.ny,sizeof(dcube.ny),1,fpout);
      nbytes += sizeof(dcube.ny );
      nobj += fwrite(&dcube.celly,sizeof(dcube.celly),1,fpout);
      nbytes += sizeof(dcube.celly );
      nobj += fwrite(&dcube.cellx,sizeof(dcube.cellx),1,fpout);
      nbytes += sizeof(dcube.cellx );
      nobj += fwrite(&dcube.beamx,sizeof(dcube.beamx),1,fpout);
      nbytes += sizeof(dcube.beamx );
      nobj += fwrite(&dcube.beamy,sizeof(dcube.beamy),1,fpout);
      nbytes += sizeof(dcube.beamy );
      for (m=0;m<dcube.nlines;m++){
      for (i=0;i<MAXNAME;i++){
      nobj += fwrite(&dcube.linenames[m][i],sizeof(char),1,fpout);
      nbytes += sizeof(char);
/*      printf("m i nobj nbytes dcube.linenames[m][i] %d %d %d %d    %c\n",m,i,nobj,nbytes,dcube.linenames[m][i]); */
      }}
      for (m=0;m<dcube.nlines;m++){
      nobj += fwrite(&dcube.restfreq[m],sizeof(dcube.restfreq[m]),1,fpout);
      }
      for (m=0;m<dcube.nlines;m++){
      for (i=0;i<dcube.nchan[m];i++) {
        nobj += fwrite(&dcube.chvel[m][i],sizeof(dcube.chvel[m][i]),1,fpout);
        nbytes += sizeof(dcube.chvel[m][i]);
/*	if (i < 10) fprintf(fpnotes,"out line %d i %d v %f\n",m,i,dcube.chvel[m][i]);
	if (i < 10) printf("out line %d i %d v %f\n",m,i,dcube.chvel[m][i]); */
      }}
      for (i=0;i<dcube.nviews;i++) {
        nobj += fwrite(&dcube.lng[i],sizeof(dcube.lng[i]),1,fpout);
	nbytes += sizeof(dcube.lng[i]);
        nobj += fwrite(&dcube.lat[i],sizeof(dcube.lat[i]),1,fpout);
	nbytes += sizeof(dcube.lat[i]);
/*	printf("view %d long %f lat %f\n",i,dcube.lng[i],dcube.lat[i]);*/
      }
      for (i=0;i<dcube.nx;i++) {
        nobj += fwrite(&dcube.xcenter[i],sizeof(dcube.xcenter[i]),1,fpout);
	nbytes += sizeof(dcube.xcenter[i]);
      }
      for (i=0;i<dcube.ny;i++) {
        nobj += fwrite(&dcube.ycenter[i],sizeof(dcube.ycenter[i]),1,fpout);
	nbytes += sizeof(dcube.ycenter[i]);
      }
      for (m=0;m<dcube.nlines;m++){
/*
      printf("nlines,outviews,noutx,nouty,nchan %d %d %d %d %d \n",
          dcube.nlines,dcube.nviews,dcube.nx,dcube.ny,dcube.nchan[m]);
*/
		      }
      for (m=0;m<dcube.nlines;m++){
      for (i=0;i<dcube.nviews;i++) {
      for (j=0;j<dcube.nx;j++) {
      for (k=0;k<dcube.ny;k++) {
      for (n=0;n<dcube.nchan[m];n++) {
        nobj += fwrite(&dcube.cube[m][i][j][k][n],
          sizeof(dcube.cube[m][i][j][k][n]),1,fpout);
	nbytes += sizeof(dcube.cube[m][i][j][k][n]);
      }}}}}
      printf("wrote %ld data to cube\n",nobj);
      printf("wrote %ld bytes to cube\n",nbytes);

  ipout = fileno(fpout);
  fflush(fpout);
  fsync(ipout);
  fclose(fpout);

  return 1;

}


int pacdata(int *iview, int *line, int *nchan, float *yray, float *zray, 
  float *cellsize, float *yspec)
{

  int i,j,k,l,v,n;
  int i1,i2,j1,j2;
  float x1,x2,y1,y2,cellratiox,cellratioy;
  float xin,yin,wt,xout,yout,weight,yspecmax;

/*
  printf("iview %d\n",*iview);
  printf("line %d\n",*line);
  printf("y %f\n",*y1);
  printf("z %f\n",*z1);
  printf("b %f\n",yspec[250]);
*/


/* Line here is not the line number, but the index of linelist which is what we
   want to store in dcube. */
  l = *line ;
  v = *iview;

/* The model cube is oriented so that at 0 longitude and 0 latitude 
   the rays run parallel to the x axis of the model cube and y,z are
   the east-west and north-south axes. These axes line up with the
   traditional cartesian and spherical axes which have the x,y plane
   horizontal and the z axis vertical.

   On the output modeled spectra we think of the east-west axis as x 
   and the north-south axis as y with the 3rd axis being frequency
   or velocity. This corresponds to the traditional axes in an
   image. 

   To go from the physical model axes to the image model axes we have
   to swap axes for the rays. We don't swap the axes for the output
   model because it was set up with image style axes.

   Here is the axis swap for the endpoints of the ray.
*/

  xin = *yray;
  yin = *zray;

  

/*
  printf("nx %d\n", dcube.nx);
  printf("ny %d\n", dcube.ny);
  printf("nc %d\n", dcube.nchan[l]);
  printf("bx,by %f %f \n", dcube.beamx,dcube.beamy);
*/


/* The cellsize for each ray is the width of the square cell in pc units.
   These cells are always square so that the area associated with a ray
   is its cellsize squared. */

    cellratiox = *cellsize/dcube.cellx;
    cellratioy = *cellsize/dcube.celly;

/* If the ray cellsize is less than the output cell, then there is at
   least one ray inside each cell. Then the ray is assigned to one
   output cell.  I2 is set to I1+1 because there is a loop below
   in which (i=i1;i<i2;i++) and the loop will run once with I=I1.
*/ 
   
    if (cellratiox <= 1.) {
	i1 = xin/(dcube.cellx) + 0.5*dcube.nx ;
	i2 = i1 + 1;
/* If the ray cellsize is larger than the output cellsize, then the 
   spectrum for this ray needs to be assigned to more than one output cell. 
*/
    } else {
        x1 = xin - (*cellsize)/2.;
        x2 = xin + (*cellsize)/2.;
        i1 = x1/(dcube.cellx) + 0.5*dcube.nx ;
        i2 = x2/(dcube.cellx) + 0.5*dcube.nx + 1;
    }

    if (cellratioy <= 1.) {
	j1 = yin/(dcube.celly) + 0.5*dcube.ny ;
	j2 = j1 + 1;
    } else {
        y1 = yin - (*cellsize)/2.;
        y2 = yin + (*cellsize)/2.;
        j1 = y1/(dcube.celly) + 0.5*dcube.ny ;
        j2 = y2/(dcube.celly) + 0.5*dcube.ny + 1;
    }

/*
    printf("ray size %7.4f cellx celly %7.4f %7.4f ratios %7.4f %7.4f xin yin %8.5f %8.5f\n",
	*cellsize,dcube.cellx,dcube.celly,cellratiox,cellratioy,xin,yin);
    fprintf(fpnotes,"ray size %7.4f cellx celly %7.4f %7.4f ratios %7.4f %7.4f xin yin %8.5f %8.5f\n",
	*cellsize,dcube.cellx,dcube.celly,cellratiox,cellratioy,xin,yin);

    printf("pacdata: x1,x2,y1,y2,i1,i2,j1,j2 %7.4f %7.4f %7.4f %7.4f %4d %4d %4d %4d\n",
	x1,x2,y1,y2,i1,i2,j1,j2);
    fprintf(fpnotes,"pacdata: x1,x2,y1,y2,i1,i2,j1,j2 %7.4f %7.4f %7.4f %7.4f %4d %4d %4d %4d\n",
	x1,x2,y1,y2,i1,i2,j1,j2);
*/

    for(i=i1;i<i2;i++) {
    for(j=j1;j<j2;j++) {

/*
   printf("line *y1 *z1 xin,yin i j %d %7.4f %7.4f %2d %2d\n",
	l,xin,yin,i,j);
   fprintf(fpnotes,"line *y1 *z1 xin,yin i j %d %7.4f %7.4f %2d %2d\n",
	l,xin,yin,i,j);
*/

/* safety check to stay inside the box*/
    if (i < 0) goto outsidebox;
    if (j < 0)  goto outsidebox;
    if (i > dcube.nx-1)  goto outsidebox;
    if (j > dcube.ny-1)  goto outsidebox;

/* The area of an output cell is dcube.cellx * dcube.celly. Weight the rays
   in each output cell according to the ratio of the input area to the
   output area.  The normalization of the weights is the sum of the input
   weights. This is stored in dcube.wt. The normalization itself is done
   in the function avg_output_cube.*/

    weight = (*cellsize/dcube.cellx) * (*cellsize/dcube.celly);

/* If the ray is bigger than the cell, then the ray
   fills the cell and the weight is 1. The weight should
   never be larger than 1. */

    if (weight > 1.0) weight = 1.0;

/*    weight = 1.0; */
    dcube.wt[v][i][j] += weight; 

    yspecmax = -1.e20;
    for (n=0;n<dcube.nchan[l];n++){
      dcube.cube[l][v][i][j][n] += yspec[n] * weight;
      if (yspec[n] > yspecmax) yspecmax = yspec[n];
    }

/*
fprintf(fpnotes,"pac lvij %3d %3d %3d %3d weight %12.3e max yspec %12.3e\n",
	l,v,i,j,weight,yspecmax);
*/

outsidebox:
	k = k; /* does nothing but prevent compiler warning */
     }}

/*
  printf("finished pacdata\n");
*/

 return 1;
}

void setup_grid( 
	int nbox,int *nx,int *ny,int *nz,
	float **xplane,float **yplane,float **zplane,
	float **xcenter,float **ycenter,float **zcenter,
	float *cellsizex, float *cellsizeyz,
	int nlng, int nlat, int *nangles, float *vlng, float *vlat,
	int outviews , float outlng[], float outlat[], int *outAngleIndex
)

/*
 outviews is the number of views in the output. This is set in setup.c.
 The vectors outlng and outlat contain the longitude and latitude for each of
 the outviews, also set in setup.c.

 nlng and nlat are specified in setup.c as the number of angles to be
 used in computing the mean radiation field.
 These get converted into a single number nangles which is nlng*nlat
 except when nlat = 4 or a multiple of 4. In this case, there is a view
 down the north pole and all the longitudes for this view are redundant.
 vlng and vlat that contain the longitude and latitude of the angles to 
 be used in the computation. 

 So we have 2 sets of angles.  outviews, outlng, outlat are what the observer 
 wants to see in the output nangles, vlng, vlat are what are used in the 
 computation of the mean radiation field. The vector outAngleIndex specifies 
 how these 2 set match up. outAngleIndex is a vector outviews long that contains 
 the index to the vectors vlng and vlat that are closest to outlng and outlat.
 Suppose nlng=3 and nlat=3. The code computes that vlng is incremented
 by 180/nlng so vlng=0,60,120. Same for vlat. If the observer wants
 to output angles 0,0 and 0,100, then outviews=2, outlng=0,0 and outlat=0,100.
 outAngleIndex = 0,6 since the first (index=0) angle in the vlng,vlat 
 list is 0,0 and the 7th (index=6) is 0,120 which is the closest to 0,100.
*/
 
{

  int i,j,k,l,m,n,ibox;
  float x,y,xs,ys,zs,a,b;

/*  Calculate and print the locations of the planes in the grid. There
    are NX cells (voxels) and NX+1 planes across one dimension of the
    model grid. */

        fprintf(fpnotes,"Here is the main computational grid\n");
	fprintf(fpnotes,"%d ",nbox );
        for (ibox=0;ibox<nbox;ibox++){fprintf(fpnotes," %d ",nx[ibox]);}
	fprintf(fpnotes," \n");
	
        for (ibox=0;ibox<nbox;ibox++){
	  xs = cellsizex[ibox]*nx[ibox];
	  ys = cellsizeyz[ibox]*ny[ibox];
	  zs = cellsizeyz[ibox]*nz[ibox];


                fprintf(fpnotes,"Box number %d \n",ibox);
                fprintf(fpnotes,"x size = %f\n",xs);
                fprintf(fpnotes,"y size = %f\n",ys);
                fprintf(fpnotes,"z size = %f\n",zs);
                fprintf(fpnotes,"voxel size = %f %f\n",cellsizex[ibox],cellsizeyz[ibox]);
                fprintf(fpnotes,"nx,ny,nz = %d, %d, %d\n",nx[ibox],ny[ibox],nz[ibox]);

                for (i=0;i<nx[ibox]+1;i++)
                        xplane[ibox][i] = 0.5*cellsizex[ibox]*(2*i - nx[ibox]);
                for (j=0;j<ny[ibox]+1;j++)
                        yplane[ibox][j] = 0.5*cellsizeyz[ibox]*(2*j - ny[ibox]);
                for (k=0;k<nz[ibox]+1;k++)
                        zplane[ibox][k] = 0.5*cellsizeyz[ibox]*(2*k - nz[ibox]);

/*
                fprintf(fpnotes,"Here are the grid planes for X\n");
                for (i=0;i<nx[ibox]+1;i++)
                        fprintf(fpnotes,"ibox i %4d %4d %9.6f \n",
                                ibox,i,xplane[ibox][i]); 
                fprintf(fpnotes,"Here are the grid planes for Y\n");
                for (j=0;j<ny[ibox]+1;j++)
                        fprintf(fpnotes,"ibox j %4d %4d %9.6f \n",
                                ibox,j,yplane[ibox][j]); 
                fprintf(fpnotes,"Here are the grid planes for Z\n");
                for (k=0;k<nz[ibox]+1;k++)
                        fprintf(fpnotes,"ibox k %4d %4d %9.6f \n",
                                ibox,k,zplane[ibox][k]); 
*/


                for (i=0;i<nx[ibox];i++)
                        xcenter[ibox][i] = 0.5*cellsizex[ibox]*(1 - nx[ibox] + 2*i);
                for (j=0;j<ny[ibox];j++)
                        ycenter[ibox][j] = 0.5*cellsizeyz[ibox]*(1 - ny[ibox] + 2*j);
                for (k=0;k<nz[ibox];k++)
                        zcenter[ibox][k] = 0.5*cellsizeyz[ibox]*(1 - nz[ibox] + 2*k);

/*
                fprintf(fpnotes,"Here are the box centers for X\n");
                for (i=0;i<nx[ibox];i++)
                        fprintf(fpnotes,"%d %9.6f \n",
                                i,xcenter[ibox][i]);
                fprintf(fpnotes,"Here are the box centers for Y\n");
                for (i=0;i<ny[ibox];i++)
                        fprintf(fpnotes,"%d %9.6f \n",
                                i,ycenter[ibox][i]);
                fprintf(fpnotes,"Here are the box centers for Z\n");
                for (i=0;i<nz[ibox];i++)
                        fprintf(fpnotes,"%d %9.6f \n",
                                i,zcenter[ibox][i]);
*/

        } /* end of loop over boxes */

/*        printf("Finished grid set up\n"); */

/* ---------------- Define the viewing angles ------------------------ */

/* Here 0 and 180 are the start and end of the angular range. 
   Use 180 range for longitude because the radiation will be computed
   forward and backward along each ray. */

/* nlng and nlat are the number of angles in the longitude and latitude 
   directions. These are input parameters that come from setup.c.
   Here we assign the angles and store them in a 1 dimensional list.
   Since all views down the z-axis (north pole) are the same, only allow
   one longitude for latitude = 90 degrees. 
*/

	k = 0;
	fprintf(fpnotes,"Here is the angular grid\n");
//	printf("Here is the angular grid\n");
	for (i=0;i<nlng;i++){
	for (j=0;j<nlat;j++){
		a = i*180./nlng;
		b = j*180./nlat;
		if ( (i == 0) || (abs(b-90.) > 4.)) {
			vlng[k] = a;
			vlat[k] = b;
			fprintf(fpnotes,"%d long, lat   %f  %f\n",k,vlng[k],vlat[k]);
//			printf("%d long, lat   %f  %f\n",k,vlng[k],vlat[k]);
			k++;
		}
	}}
	*nangles = k;

/* outviews, outlng and outlat are the number of views to be output and their angles.
   These come from setup.c.
   In the loop below, for each requested angle (outlng and outlat) find the nearest 
   computed angle (vlng and vlat). 
*/

//	printf("locate angles %d %d\n",*nangles, outviews);
	for (k=0;k<outviews;k++) {
//	printf("view %d outlng outlat %f %f\n",k,outlng[k],outlat[k]);  
	x = 1.e6;
	for (i=0;i<*nangles;i++){
		y = (outlat[k] - vlat[i])*(outlat[k] - vlat[i])
                       + (outlng[k] - vlng[i])*(outlng[k] - vlng[i]) ;
//	printf("k i %d %d vlng vlat %f %f x y %f %f\n",k,i,vlng[i],vlat[j],x,y); 
                if (y < x) {
			m = i;
			x = y;
//	printf("m x %d %f\n",m,x); 
		}
        }
//		printf("set outAngleIndex %d to %d\n",k,m);
		outAngleIndex[k] = m;
//		printf("view number %d long, lat   %f  %f\n",k, vlng[m],vlat[m]);
//		fprintf(fpnotes,"view number %d long, lat   %f  %f\n",k, vlng[m],vlat[m]);
	}

/* Replace the requested angles with the angles */
	for (i=0;i<outviews;i++){
		outlng[i] = vlng[outAngleIndex[i]];
		outlat[i] = vlat[outAngleIndex[i]];
		fprintf(fpnotes,"viewing angle %d long, lat   %f  %f\n",i, outlng[i],outlat[i]);
//		printf("viewing angle %d long, lat   %f  %f\n",i, outlng[i],outlat[i]);
	}

//		printf("finished grid setup \n");

}

void set_linenames(int molecule, char ** linenames)
{

int i,j,k,line;
char qn[8];

  for (line=0;line<NLINES;line++){
    linenames[line] = (char *)malloc(MAXNAME*sizeof(char));
  }

  if (molecule == 2) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"CH3CN");
      strcat(linenames[line],"(");
      sprintf(qn,"%d",atomic.jn[atomic.indexu[line]]);
      strcat(linenames[line],qn);
      strcat(linenames[line],",");
      sprintf(qn,"%d",atomic.kn[atomic.indexu[line]]);
      strcat(linenames[line],qn);
      strcat(linenames[line],"-");
      sprintf(qn,"%d",atomic.jn[atomic.indexl[line]]);
      strcat(linenames[line],qn);
      strcat(linenames[line],",");
      sprintf(qn,"%d",atomic.kn[atomic.indexl[line]]);
      strcat(linenames[line],qn);
      strcat(linenames[line],")");
/*      printf("%d %s\n",i,linenames[line]); */
    }
  } 

  if (molecule == 3 ) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"NH3");
      strcat(linenames[line],"(");
      sprintf(qn,"%d",atomic.jn[atomic.indexu[line]]);
      strcat(linenames[line],qn);
      strcat(linenames[line],",");
      sprintf(qn,"%d",atomic.kn[atomic.indexu[line]]);
      strcat(linenames[line],qn);
      strcat(linenames[line],")");
/*      printf("%d %s\n",i,linenames[line]); */
    }
  } 

  if (molecule == 22) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"H2D:");
      sprintf(qn,"%d",line+1);
      strcat(linenames[line],qn);
    }
    strcpy(linenames[0],"H2D+(110->111)");
    if (NLINES > 1) {
      line = 1;
      strcpy(linenames[line],"H2D+(212->111)");
    }
    if (NLINES > 2) {
      line = 2;
      strcpy(linenames[line],"H2D+(211->110)");
    }
    if (NLINES > 3) {
      line = 3;
      strcpy(linenames[line],"H2D+(211->212)");
    }
    if (NLINES > 4) {
      line = 4;
      strcpy(linenames[line],"H2D+(313->212)");
    }
    if (NLINES > 5) {
      line = 5;
      strcpy(linenames[line],"H2D+(331->212)");
    }
    if (NLINES > 6) {
      line = 6;
      strcpy(linenames[line],"H2D+(312->211)");
    }
    if (NLINES > 7) {
      line = 7;
      strcpy(linenames[line],"H2D+(330->211)");
    }
    if (NLINES > 8) {
      line = 8;
      strcpy(linenames[line],"H2D+(312->313)");
    }
    if (NLINES > 9) {
      line = 9;
      strcpy(linenames[line],"H2D+(330->313)");
    }
    if (NLINES > 10) {
      line = 10;
      strcpy(linenames[line],"H2D+(414->313)");
    }
    if (NLINES > 11) {
      line = 11;
      strcpy(linenames[line],"H2D+(331->312)");
    }
    if (NLINES > 12) {
      line = 12;
      strcpy(linenames[line],"H2D+(413->312)");
    }
    if (NLINES > 13) {
      line = 13;
      strcpy(linenames[line],"H2D+(431->312)");
    }
    if (NLINES > 14) {
      line = 14;
      strcpy(linenames[line],"H2D+(413->414)");
    }
    if (NLINES > 15) {
      line = 15;
      strcpy(linenames[line],"H2D+(431->414)");
    }
    if (NLINES > 16) {
      line = 16;
      strcpy(linenames[line],"H2D+(515->414)");
    }
    if (NLINES > 17) {
      line = 17;
      strcpy(linenames[line],"H2D+(330->331)");
    }
    if (NLINES > 18) {
      line = 18;
      strcpy(linenames[line],"H2D+(432->331)");
    }
    if (NLINES > 19) {
      line = 19;
      strcpy(linenames[line],"H2D+(413->330)");
    }
    if (NLINES > 20) {
      line = 20;
      strcpy(linenames[line],"H2D+(431->330)");
    }
    if (NLINES > 21) {
      line = 21;
      strcpy(linenames[line],"H2D+(432->413)");
    }
    if (NLINES > 22) {
      line = 22;
      strcpy(linenames[line],"H2D+(514->413)");
    }
    if (NLINES > 23) {
      line = 23;
      strcpy(linenames[line],"H2D+(514->515)");
    }
    if (NLINES > 24) {
      line = 24;
      strcpy(linenames[line],"H2D+(616->515)");
    }
    if (NLINES > 25) {
      line = 25;
      strcpy(linenames[line],"H2D+(431->432)");
    }
    if (NLINES > 26) {
      line = 26;
      strcpy(linenames[line],"H2D+(533->432)");
    }
    if (NLINES > 27) {
      line = 27;
      strcpy(linenames[line],"H2D+(514->431)");
    }
    if (NLINES > 28) {
      line = 28;
      strcpy(linenames[line],"H2D+(532->431)");
    }
    if (NLINES > 29) {
      line = 29;
      strcpy(linenames[line],"H2D+(533->514)");
    }
    if (NLINES > 30) {
      line = 30;
      strcpy(linenames[line],"H2D+(615->514)");
    }
  } 
  if (molecule == 6) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"H2O:");
      sprintf(qn,"%d",line+1);
      strcat(linenames[line],qn);
    }
    strcpy(linenames[0],"H2O(110->101)");
    if (NLINES > 1) {
      line = 1;
      strcpy(linenames[line],"H2O(212->101)");
    }
    if (NLINES > 2) {
      line = 2;
      strcpy(linenames[line],"H2O(221->110)");
    }
    if (NLINES > 3) {
      line = 3;
      strcpy(linenames[line],"H2O(221->212)");
    }
    if (NLINES > 4) {
      line = 4;
      strcpy(linenames[line],"H2O(303->212)");
    }
    if (NLINES > 5) {
      line = 5;
      strcpy(linenames[line],"H2O(312->303)");
    }
    if (NLINES > 6) {
      line = 6;
      strcpy(linenames[line],"H2O(312->221)");
    }
    if (NLINES > 7) {
      line = 7;
      strcpy(linenames[line],"H2O(321->312)");
    }
    if (NLINES > 8) {
      line = 8;
      strcpy(linenames[line],"H2O(321->212)");
    }
    if (NLINES > 9) {
      line = 9;
      strcpy(linenames[line],"H2O(414->321)");
    }
    if (NLINES > 10) {
      line = 10;
      strcpy(linenames[line],"H2O(414->303)");
    }
    if (NLINES > 11) {
      line = 11;
      strcpy(linenames[line],"H2O(330->321)");
    }
    if (NLINES > 12) {
      line = 12;
      strcpy(linenames[line],"H2O(330->303)");
    }
    if (NLINES > 13) {
      line = 13;
      strcpy(linenames[line],"H2O(330->221)");
    }
    if (NLINES > 14) {
      line = 14;
      strcpy(linenames[line],"H2O(423->330)");
    }
    if (NLINES > 15) {
      line = 15;
      strcpy(linenames[line],"H2O(423->414)");
    }
    if (NLINES > 16) {
      line = 16;
      strcpy(linenames[line],"H2O(423->312)");
    }
  } 
  if (molecule == 7) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"CO");
      strcat(linenames[line],"(");
      sprintf(qn,"%d",line+1);
      strcat(linenames[line],qn);
      strcat(linenames[line],"-");
      sprintf(qn,"%d",line);
      strcat(linenames[line],qn);
      strcat(linenames[line],")");
    }
  } 
  if (molecule == 8) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"CS");
      strcat(linenames[line],"(");
      sprintf(qn,"%d",line+1);
      strcat(linenames[line],qn);
      strcat(linenames[line],"-");
      sprintf(qn,"%d",line);
      strcat(linenames[line],qn);
      strcat(linenames[line],")");
    }
  } 
  if (molecule == 9 || molecule == 4) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"N2H+");
      strcat(linenames[line],"(");
      sprintf(qn,"%d",line+1);
      strcat(linenames[line],qn);
      strcat(linenames[line],"-");
      sprintf(qn,"%d",line);
      strcat(linenames[line],qn);
      strcat(linenames[line],")");
    }
  } 
  if (molecule == 10) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"13CO");
      strcat(linenames[line],"(");
      sprintf(qn,"%d",line+1);
      strcat(linenames[line],qn);
      strcat(linenames[line],"-");
      sprintf(qn,"%d",line);
      strcat(linenames[line],qn);
      strcat(linenames[line],")");
    }
  }
  if (molecule == 12) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"HCO+");
      strcat(linenames[line],"(");
      sprintf(qn,"%d",line+1);
      strcat(linenames[line],qn);
      strcat(linenames[line],"-");
      sprintf(qn,"%d",line);
      strcat(linenames[line],qn);
      strcat(linenames[line],")");
    }
  } 
  if (molecule == 13) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"C17O");
      strcat(linenames[line],"(");
      sprintf(qn,"%d",line+1);
      strcat(linenames[line],qn);
      strcat(linenames[line],"-");
      sprintf(qn,"%d",line);
      strcat(linenames[line],qn);
      strcat(linenames[line],")");
    }
  }
  if (molecule == 14) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"C18O");
      strcat(linenames[line],"(");
      sprintf(qn,"%d",line+1);
      strcat(linenames[line],qn);
      strcat(linenames[line],"-");
      sprintf(qn,"%d",line);
      strcat(linenames[line],qn);
      strcat(linenames[line],")");
    }
  }
  if (molecule == 15 || molecule == 21) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"HCN");
      strcat(linenames[line],"(");
      sprintf(qn,"%d",line+1);
      strcat(linenames[line],qn);
      strcat(linenames[line],"-");
      sprintf(qn,"%d",line);
      strcat(linenames[line],qn);
      strcat(linenames[line],")");
    }
  }
  if (molecule == 17) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"HCN,v2=");
      sprintf(qn,"%d",atomic.kn[line]);
      strcat(linenames[line],qn);
      strcat(linenames[line],"(");
      sprintf(qn,"%d",atomic.upper[line]);
      strcat(linenames[line],qn);
      strcat(linenames[line],"-");
      sprintf(qn,"%d",atomic.lower[line]);
      strcat(linenames[line],qn);
      strcat(linenames[line],")");
    }
  }
  if (molecule == 16) {
    for (line=0;line<NLINES;line++){
      strcpy(linenames[line],"N2D+");
      strcat(linenames[line],"(");
      sprintf(qn,"%d",line+1);
      strcat(linenames[line],qn);
      strcat(linenames[line],"-");
      sprintf(qn,"%d",line);
      strcat(linenames[line],qn);
      strcat(linenames[line],")");
    }
  } 

  if (molecule < 2 || molecule > 22
      ) {
    printf("No line names for this molecule.\n");
    fprintf(fpnotes,"No line names for this molecule.\n");
    for (line=0;line<NLINES;line++){
       strcpy(linenames[line],"                ");
    }
    return;
  }


  for (line=0;line<NLINES;line++){
    j = strlen(linenames[line]);
/*    printf("line %d string length %d %s\n",line,j,linenames[line]);  */
    k = MAXNAME - j ;
    for (i=0;i<k;i++){
      sprintf(qn,"%d",i%10);
      strcat(linenames[line]," ");
    }
    strcat(linenames[line],"\0"); 
    j = strlen(linenames[line]); 
/*    printf("line %d string length %d %s\n",line,j,linenames[line]);  */
  }

}




int read_output_cube(int best_model)

{

  char fname[32], suffix[4];
  FILE *fpmodel;
  int i,j,k,l,m,n;
  int nobj,ipmodel;

        strcpy(fname,"ModelCube");
        sprintf(suffix,"%d",best_model);
        strcat(fname,suffix); 
        fpmodel = fopen(fname,"r");
             if (fpmodel == NULL) {
                printf("Could not open %s. Exiting ...\n",fname);
                exit(0);
             } else {
                printf("File %s opened\n",fname);
             }




/* read dcube structure */

      nobj = fread(&dcube.nlines,sizeof(dcube.nlines),1,fpmodel);
      nobj = fread(&dcube.nviews,sizeof(dcube.nviews),1,fpmodel);

      printf("Start reading output cubes. nlines outviews %d %d\n",
        dcube.nlines,dcube.nviews);

      dcube.nchan = (int *)malloc(dcube.nlines*sizeof(int));
      for (m=0;m<dcube.nlines;m++) {
        nobj = fread(&dcube.nchan[m],sizeof(dcube.nchan[m]),1,fpmodel);
        printf("The number of channels in line %d is %d \n",m,dcube.nchan[m]);
      }

      nobj = fread(&dcube.nx,sizeof(dcube.nx),1,fpmodel);
      nobj = fread(&dcube.ny,sizeof(dcube.ny),1,fpmodel);
      printf("The grid size is %d  x  %d\n",dcube.nx,dcube.ny);

      nobj = fread(&dcube.celly,sizeof(dcube.celly),1,fpmodel);
      nobj = fread(&dcube.cellx,sizeof(dcube.cellx),1,fpmodel);
      nobj = fread(&dcube.beamx,sizeof(dcube.beamx),1,fpmodel);
      nobj = fread(&dcube.beamy,sizeof(dcube.beamy),1,fpmodel);
      printf("The output cell size is %f  x  %f\n",dcube.cellx,dcube.celly);
      printf("The output beam size is %f  x  %f\n",dcube.beamx,dcube.beamy);

      dcube.linenames = (char **) malloc(dcube.nlines*sizeof(char *));
      for (m=0;m<dcube.nlines;m++){
        dcube.linenames[m] = (char *) malloc(MAXNAME*sizeof(char));
      }
      for (m=0;m<dcube.nlines;m++){
      for (i=0;i<MAXNAME;i++){
        nobj = fread(&dcube.linenames[m][i],sizeof(char),1,fpmodel);
      }
        dcube.linenames[m][MAXNAME]= '\0';
        printf("m linenames %d %s\n",m,dcube.linenames[m]);
      }

      dcube.chvel = (float **)malloc(dcube.nlines*sizeof(float *));
      for (i=0;i<dcube.nlines;i++){
         dcube.chvel[i] = (float *)malloc(dcube.nchan[i]*sizeof(float));
      }
      for (m=0;m<dcube.nlines;m++){
      printf("Velocities for line %s\n",dcube.linenames[m]);
      for (i=0;i<dcube.nchan[m];i++) {
        nobj = fread(&dcube.chvel[m][i],sizeof(dcube.chvel[m][i]),1,fpmodel);
/*        printf("%d %f \n",i,dcube.chvel[m][i]/1.e5); */
      }}
      printf("channel velocities read in.\n");

      dcube.lng = (float *)malloc(dcube.nviews*sizeof(float));
      dcube.lat = (float *)malloc(dcube.nviews*sizeof(float));
      for (i=0;i<dcube.nviews;i++) {
        nobj = fread(&dcube.lng[i],sizeof(dcube.lng[i]),1,fpmodel);
        nobj = fread(&dcube.lat[i],sizeof(dcube.lat[i]),1,fpmodel);
        printf("view %d, longitude latitude %f %f\n",i,dcube.lng[i],dcube.lat[i]);
      }

      dcube.xplane = (float *)malloc((dcube.nx+1)*sizeof(float));
      dcube.yplane = (float *)malloc((dcube.ny+1)*sizeof(float));
      dcube.xcenter = (float *)malloc(dcube.nx*sizeof(float));
      dcube.ycenter = (float *)malloc(dcube.ny*sizeof(float));
      printf("x grid centers\n");
      for (i=0;i<dcube.nx;i++) {
        nobj = fread(&dcube.xcenter[i],sizeof(dcube.xcenter[i]),1,fpmodel);
/*        printf("%d %f\n",i,dcube.xcenter[i]); */
      }
      printf("x grid centers\n");
      for (i=0;i<dcube.ny;i++) {
        nobj = fread(&dcube.ycenter[i],sizeof(dcube.ycenter[i]),1,fpmodel);
/*        printf("%d %f\n",i,dcube.ycenter[i]); */
      }
      printf("cell centers read in\n");


      dcube.cube = (float *****) malloc(dcube.nlines*sizeof(float****));
      for (i=0;i<dcube.nlines;i++){
        dcube.cube[i] = (float ****) malloc(dcube.nviews*sizeof(float***));
      for (j=0;j<dcube.nviews;j++){
        dcube.cube[i][j] = 
          (float ***) malloc(dcube.nx*sizeof(float**));
        for (k=0;k<dcube.nx;k++){
          dcube.cube[i][j][k] = 
          (float **) malloc(dcube.ny*sizeof(float*));
          for (n=0;n<dcube.ny;n++){
            dcube.cube[i][j][k][n] = 
             (float *) malloc(dcube.nchan[i]*sizeof(float));
      }}}}

      printf("OK on cube allocation\n");

      for (m=0;m<dcube.nlines;m++){
      for (i=0;i<dcube.nviews;i++) {
      for (j=0;j<dcube.nx;j++) {
      for (k=0;k<dcube.ny;k++) {
      for (n=0;n<dcube.nchan[m];n++) {
        nobj = fread(&dcube.cube[m][i][j][k][n],
          sizeof(dcube.cube[m][i][j][k][n]),1,fpmodel);
      }}}}}

      printf("data cube read in\n");


  ipmodel = fileno(fpmodel);
  close(ipmodel);

  return 1;
}





int simplex(int n, float newy, float *newp)
{

int i,j,k;
static int ihi,inhi,ilo,firstcall=1;
#define MAXN 20
static float psum[MAXN],deltap[MAXN];
static int step=0,count=0;
static float stepsize;
static float ysave;
static int nextensions;
static float *y;
static float **p;



/* firstcall remains 1 until simplex is called n+1 times. 
   These n+1 calls initialize the chisq values at the vertices. */

    if (firstcall == 1) {
        printf("First call\n");

/* First time called ever, allocate memory for vertices and initialize */

        if (count == 0) {
            y = (float *)malloc((n+1)*sizeof(float));
            p = (float **)malloc((n+1)*sizeof(float *));
            for (i=0;i<n+1;i++){
                  p[i] = (float *)malloc(n*sizeof(float));
            } 
            for (i=0;i<n+1;i++){
            for (j=0;j<n;j++){
              p[i][j] = 0.;
            }}
            for (i=0;i<n;i++){
	          p[i+1][i] = 1.;
            }
        }

/* If we have finished with the initial n+1 calls, 
   all the vertices of the simplex have chisq values. 
   Set firstcall to 0 (false) and we can start the 
   minimization. */

        if (count > n) {
            printf("count-1 = %d\n",count-1);
            y[count-1] = newy;
            printf("assigned \n");
            printf("Initial Simplex \n");
            for (i=0;i<n+1;i++){
              printf("i %2d ",i);
            for (j=0;j<n;j++){
              printf(" %f ",p[i][j]);
            }
              printf(" %f \n",y[i]);
            }
            reorder(n, &ilo, &inhi, &ihi, y);
            printf("reordered\n");
            newStep(n, ihi, p, psum, deltap, newp);
            printf("new step\n");

            step = 1;
            nextensions = 0;
            firstcall = 0;
            count = 0;
            return 1;
        }

/* Here we are still in the initialization process, working
   through the n+1 calls to set chisq at the simplex vertices */

    	    for (j=0;j<n;j++) newp[j] = p[count][j];
            if (count > 0) {
                y[count-1] = newy;
                printf("Initial value for vertex %d  y = %f\n",
                    count,y[count-1]);
            }
            printf("Initializing vertex %d\n",count);
            for (j=0;j<n;j++) { 
              printf(" %f ",newp[j]);
            }
            printf("\n");
            count++;
            return 1;
      

    }

  printf("Current Simplex: \n");
            for (i=0;i<n+1;i++){
              printf("i %2d ",i);
            for (j=0;j<n;j++){
              printf(" %f ",p[i][j]);
            }
              printf(" %f \n",y[i]);
            }

  switch (step) {
   case 1:

/* Based on the result of the first step, choose the next stepsize */
  
/* Compute a new step if y is between y lo and y next highest  or
   if we have already done one extension */

        if ( (newy > y[ilo] && newy < y[inhi]) || nextensions >= 1 ) {
          printf("Regular step\n");

/* The if statement below covers the case where we get here by
   the number of extensions */

          if (newy < y[ihi]) replaceHi(n, ihi, newy, y, p, newp, psum) ;
          reorder(n, &ilo, &inhi, &ihi, y);
          newStep(n, ihi, p, psum, deltap, newp);
          step = 1;
          nextensions = 0;
          return 1;
        }

/* In these cases repeat the old step with a different size or direction */
        printf("repeat step\n");
        if (newy <= y[ilo])   stepsize = 1.0;
        if (newy >= y[inhi])  stepsize = -0.5;

/* Replace y hi with the new y for the 2 cases above, but 
   do not replace if new y is higher than y hi */
/* Do not re-order because on the repeat we want to replace
   the first point if the repeat point is better. Otherwise
   we could have 2 vertices on the same line from the
   2 extensions. */

        if (newy > y[ihi])  { stepsize = -1.5;
        } else {
            replaceHi(n, ihi, newy, y, p, newp, psum) ;
        }

        ysave = y[ihi];

	for (j=0;j<n;j++) { 
          newp[j] += stepsize*deltap[j];
        }

        printf("Case 1: stepsize %f  : New point: ",stepsize);
	for (j=0;j<n;j++) { 
          printf(" %f ",newp[j]);
        }
        step = 2;
        if (stepsize > -1.)  step = 1; 
        nextensions++;
        printf("  Step set to %d next = %d \n",step,nextensions);
        return 1;

  case 2:

/* If the repeat move was successful, compute a new step  */

        if (newy < ysave) {
          printf("Repeat was successful. Do a regular step. \n");
          replaceHi(n, ihi, newy, y, p, newp, psum) ;
          reorder(n, &ilo, &inhi, &ihi, y);
          newStep(n, ihi, p, psum, deltap, newp);
          step = 1;
          nextensions = 0;
          return 1;
        }

/* If the new solution is no better than the high point, then
   contract around the low point. */

          count = 0;
          step = 3;
          printf("Step set to %d\n",step);

  case 3:

        if (count == ilo) count++;
        if (count > n) {
           step = 1;
           reorder(n, &ilo, &inhi, &ihi, y);
           newStep(n, ihi, p, psum, deltap, newp);
           return 1;
        }

        for (j=0;j<n;j++) { 
          newp[j] = 0.5*(p[count][j] + p[ilo][j]);
          p[count][j] = newp[j]; 
        }
        printf("Case 3: count = %d contract: New point: ",count);
        for (j=0;j<n;j++) { 
          printf(" %f ",newp[j]);
        }
        printf("\n");
        count++;
        return 1;

  } /* end of switch */

   printf("Can not get here \n");
   exit(0);

} /* end of function simplex */

int replaceHi(int n, int ihi, float newy, float *y, float**p, float *newp, float psum[])

{
int i,j;

      y[ihi] = newy;
        for (j=0;j<n;j++)  {
          p[ihi][j] = newp[j];
          psum[j] += newp[j] - p[ihi][j];
        }
        printf("Replaced highest point\n");
  printf("Current Simplex \n");
    for (i=0;i<n+1;i++){
    printf("i %d vector %f %f %f value %f\n",i,p[i][0],p[i][1],p[i][2],y[i]);
  }
  return 1;
}


int newStep(int n, int ihi, float **p, float psum[], float deltap[], float * newp)
{

int i,j;

    /* Find the centroid of the simplex without the high point */

        for (j=0;j<n;j++) { 
          psum[j] = 0.0;
          for (i=0;i<n+1;i++) {
              if (i != ihi) {
                psum[j] += p[i][j]; 
    /*            printf(" j %d %d pij %f psum %f\n",i,j,p[i][j],psum[j]); */
              }
          }
          psum[j] /= n;
        }

        printf("Centroid: ");
        for (j=0;j<n;j++) { 
           printf(" %f ",psum[j]);
        }
        printf("\n");

	for (j=0;j<n;j++) { 
          deltap[j] = psum[j] - p[ihi][j];
        }
        printf("Difference: ");
	for (j=0;j<n;j++) { 
          printf(" %f ",deltap[j]);
        }
        printf("\n");
 
	for (j=0;j<n;j++) { 
          newp[j] = psum[j] + deltap[j];
        }
        printf("New point: ");
	for (j=0;j<n;j++) { 
          printf(" %f ",newp[j]);
        }
        printf("\n");

  return 1;
}

int reorder(int n, int *ilo, int *inhi, int *ihi, float *y)
{

/* ihi = highest; inhi = next highest; ilo = lowest */

int i;

     *ilo=0;
     for (i=0;i<n+1;i++) if (y[i] < y[*ilo]) *ilo=i;

     *ihi = 0;
     for (i=0;i<n+1;i++) if (y[i] > y[*ihi]) *ihi = i;

     if (*ihi != 0) {
       *inhi = 0;
     } else {
       *inhi = 1;
     }
     for (i=0;i<n+1;i++) if (y[i] > y[*inhi] && i != *ihi) *inhi = i;

     printf("Highest point %d next highest %d lowest %d\n",*ihi,*inhi,*ilo);

  return 1;
}

int anneal(int n, float anneal_t0, float ** plim, float * pn, float yn, int it, float chisqbest)
{

/* yn		input value of new chisq
   y		current value of chisq at position p
   p[n]		current position vector
   pn[n]	new position vector
   t		current annealing temperature
   dy[n]	difference between current and new chisq
   dp[n]	difference between current and new position
   v		current position variable to be changed
   pw[n]	width of probability distribution for new positions
*/

  int i,j,k;
  static int firstcall=1,ic=0;
  static int v,seed;
  static float y,t0,t,fac=1.2,avgdy,new_width;
  static float *dy;
  static float *sdy;
  static float *dp;
  static float *p;
  static int *app;
  static float *yv;
  static float *pw;

static float maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2)) 
#define MIN(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) < (maxarg2) ?\
        (maxarg1) : (maxarg2)) 

  if (firstcall == 1) {

    dy = (float *)malloc(n*sizeof(float));
    sdy = (float *)malloc(n*sizeof(float));
    dp = (float *)malloc(n*sizeof(float));
    pw = (float *)malloc(n*sizeof(float));
    p  = (float *)malloc(n*sizeof(float));
    app  = (int *)malloc(n*sizeof(int));
    yv  = (float *)malloc(n*sizeof(float));

/* Initialize random number generator by calling with negative arg*/
    seed = -2;
    ran1(&seed);

/* Initialize widths for probability distributions */
    for (i=0;i<n;i++) {
      pw[i] = 1.;
      p[i] =  0.;
      pn[i] = 0.;
    }

    v = 0;
    y = yn;
    t0 = anneal_t0; /* This looks arbitrary. Should be about chisq. Maybe the average over the first n rounds */
    t = t0;

    firstcall = 2;
    return 1;

  } /* end of firstcall == 1 */

  if (firstcall == 2) {
/*
    printf("starting 2nd call firstcall=%d\n",firstcall);
    fprintf(fpnotes,"starting 2nd call firstcall=%d\n",firstcall);
*/
    y = yn;
    pn[v] = p[v] + xcl(&seed,pw[v]);
    printf("new vertex=%d  old value %e  new %e\n",v,p[v],pn[v]);
    fprintf(fpnotes,"new vertex=%d  old value %e new %e\n",v,p[v],pn[v]);
    firstcall = 0; 
    return 1;
  }

/* we have a new chisq yn at new position vector pn 
   previous vector is still in p */

   dy[v] = yn - y;
   yv[v] = yn;
   dp[v] = pn[v] - p[v];

   printf("it %d v %d y %f yn %f dy %f dp %f\n",it,v,y,yn,dy[v],dp[v]);
   fprintf(fpnotes,"it %d v %d y %f yn %f dy %f dp %f\n",it,v,y,yn,dy[v],dp[v]);

/* Decide whether to accept the current step */

   app[v] = 1;
   if (move(dy[v],t,&seed) == 1) {
     app[v] = 1;
     p[v] = pn[v];
     y = yn;
     printf("Updating for approval  %f %f\n",p[v],y);
     fprintf(fpnotes,"Updating for approval %f %f\n",p[v],y);
   } else {
     app[v] = 0;
     pn[v] = p[v];
     printf("restoring disapproval  %f %f\n",pn[v],y);
     fprintf(fpnotes,"restoring disapproval  %f %f\n",pn[v],y);
   }

   printf("Current position:\n");
   for (i=0;i<n;i++) printf(" %f ",p[i]);
   printf(" chisq %f\n",y);
   fprintf(fpnotes,"Current position:\n");
   for (i=0;i<n;i++) fprintf(fpnotes," %f ",p[i]);
   fprintf(fpnotes," chisq %f\n",y);

/* Next variable */
  v++;

  if (v == n) {
/* Finished a pass through each of the variables. 
   Chill the solution. */
     v = 0;
     t = t0/pow(it,1./n) ;
     printf("Chill factor %f  New T %f\n",pow(it,1./n),t);
     fprintf(fpnotes,"Chill factor %f  New T %f\n",pow(it,1./n),t);
     avgdy = 0.;
     for (i=0;i<n;i++) {
      sdy[i] = dy[i];
       dy[i] = MAX(fabs(dy[i]),1.e-6);
       avgdy += dy[i];
/* Change jump width if move for this var was approved. Limit magnitude of change to fac */
/* If the solution has drifted far off the global minimum, then keep the jump large             */
/* Try adjusting for both approved and non-approved moves */
       if (app[i] == 1 || app[i] == 0) {			
           new_width = t*fabs(dp[i]/dy[i]);
           if (new_width > pw[i]) {
/* this adjustment was to keep the width not too high following bad guesses
   but I didn't like the result 
               pw[i] = MIN(pw[i]*fac ,t*fabs(dp[i]/dy[i])*yv[i]/chisqbest);
*/
               pw[i] = MIN(pw[i]*fac ,t*fabs(dp[i]/dy[i]));
           } else {
               pw[i] = MAX(pw[i]/fac ,  t * fabs(dp[i]/dy[i]));
           } 
       }
       printf("i app dy nw pw %d %d %f %f %f\n",i,app[i],sdy[i],new_width,pw[i]);
       fprintf(fpnotes,"i app dy nw pw %d %d %f %f %f\n",i,app[i],sdy[i],new_width,pw[i]);
     }
     printf("Average jump %f compared to previous T %f\n",
        avgdy/n,t0/pow(it,1./n));
     fprintf(fpnotes,"Average jump %f compared to previous T %f\n",
        avgdy/n,t0/pow(it,1./n));
  }

/* Generate new vector for next variable */

  for (i=0;i<2000;i++) {
    pn[v] = p[v] + xcl(&seed,pw[v]);
    if (pn[v] > plim[v][0] && pn[v] < plim[v][1]) break;
  }
  if (i == 2000) {
     printf("Failed to find good variable in 1000 tries.\n");
     printf("Check limits for variable %d: %f %f %f %f\n",v,plim[v][0],plim[v][1],p[v],pw[v]);
     fprintf(fpnotes,"Failed to find good variable in 1000 tries.\n");
     fprintf(fpnotes,"Check limits for variable %d: %f %f %f %f\n",
         v,plim[v][0],plim[v][1],p[v],pw[v]);
     return 0;
  }
  printf("v pn %d %f\n",v,pn[v]);
  fprintf(fpnotes,"v pn %d %f\n",v,pn[v]);
   printf("Next position:\n");
   for (i=0;i<n;i++) printf(" %f ",pn[i]);
   printf("\n");
   fprintf(fpnotes,"Next position:\n");
   for (i=0;i<n;i++) fprintf(fpnotes," %f ",pn[i]);
   fprintf(fpnotes,"\n");


/* Check new vector for allowed range */

/* Go back and compute new chisq */

   return 1; 

}

float xcl(int * seed,float t)
{
/* Returns a random variable with a Cauchy-Lorentz distribution
   with zero mean and half-width t.
   Generates a uniform deviate using ran1 and transforms that
   to a C-L distribution. 
   Setting seed to 0 or negative
   re-initializes the random sequence from ran1. */

   return t * tan(ran1(seed) * PI);
}

float pcl(float x, float t)
{

/* Cauchy-Lorentz distribution function. Returns the probability
   of obtaining x from a distribution of width t. 

   return t*t / (t*t +  x*x);
   printf("t %f x %f pcl %f\n",t,x,t*t / (t*t +  x*x));
   fprintf(fpnotes,"t %f x %f pcl %f\n",t,x,t*t / (t*t +  x*x));
*/

/* Exponential probability function. Theoretically, use the exponential. */
   return exp(-x/t);
   
}
 
int move(float dy,float t,int *seed)
{
float r,p;
    r = ran1(seed);
    p = pcl(dy,t);
    if (dy < 0. || r < p) {
       if (dy < 0.) printf("Move approved dy %f r %f p %f\n",dy,r,p);
       if (dy < 0.) fprintf(fpnotes,"Move approved dy %f r %f p %e\n",dy,r,p);
       if (dy >= 0.) printf("MOVE APPROVED ANYWAY dy %f r %f p %f\n",dy,r,p);
       if (dy >= 0.) fprintf(fpnotes,"MOVE APPROVED ANYWAY dy %f r %f p %f\n",dy,r,p);
       return 1;
    } else {
       printf("Move not approved dy %f r %f p %f\n",dy,r,p);
       fprintf(fpnotes,"Move not approved dy %f r %f p %f\n",dy,r,p);
       return 0;
    }
    
}

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

float ran1(idum)
int *idum;
{
	static long ix1,ix2,ix3;
	static float r[98];
	float temp;
	static int iff=0;
	int j;

	if (*idum < 0 || iff == 0) {
		iff=1;
		ix1=(IC1-(*idum)) % M1;
		ix1=(IA1*ix1+IC1) % M1;
		ix2=ix1 % M2;
		ix1=(IA1*ix1+IC1) % M1;
		ix3=ix1 % M3;
		for (j=1;j<=97;j++) {
			ix1=(IA1*ix1+IC1) % M1;
			ix2=(IA2*ix2+IC2) % M2;
			r[j]=(ix1+ix2*RM2)*RM1;
		}
		*idum=1;
	}
	ix1=(IA1*ix1+IC1) % M1;
	ix2=(IA2*ix2+IC2) % M2;
	ix3=(IA3*ix3+IC3) % M3;
	j=1 + ((97*ix3)/M3);
	if (j > 97 || j < 1) printf("RAN1: This cannot happen.");
	temp=r[j];
	r[j]=(ix1+ix2*RM2)*RM1;
	return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3


int  check_model(int model_number, int nbox, int *nx, int *ny, int *nz,
                 float **xcenter,  float **ycenter,  float **zcenter, 
                 float lwmin,float lwmax) 
{

  int i,j,k,ibox;
  float lwmin_model,lwmax_model;

/* This makes a constant, temperature, density, abundance model
regardless of the definition in define_model.c
for (ibox=0;ibox<nbox;ibox++){
*/
/*
                for (i=0;i<nx[ibox];i++){
                for (j=0;j<ny[ibox];j++){
                for (k=0;k<nz[ibox];k++){
                        model[ibox][i][j][k].abundance = 1.e-6;
                        model[ibox][i][j][k].temperature = 2.71000;
                        model[ibox][i][j][k].density = 1.e1;
}}}
*/

        for (ibox=0;ibox<nbox;ibox++){
		i = nx[ibox]/2 ;
		j = ny[ibox]/2 ;
		fprintf(fpnotes,"Box number %d out of %d\n",ibox,nbox);
		fprintf(fpnotes,"Number of cells %d\n",nz[ibox]);
		fprintf(fpnotes,"Check the model at i,j = %d %d x,y = %f %f\n",
			i,j,xcenter[ibox][i],ycenter[ibox][j]);
		fprintf(fpnotes,"  k      z       T     Tdust  density       ");
		fprintf(fpnotes,"vx         vy         vz     linewidth  mean Av    abundance\n");

        	for (k=0;k<nz[ibox];k++){
	fprintf(fpnotes,"%3d %9.6f %6.2f %6.2f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
			k,zcenter[ibox][k],
			model[ibox][i][j][k].temperature,
			model[ibox][i][j][k].dust_temperature,
			model[ibox][i][j][k].density,
			model[ibox][i][j][k].vx,
			model[ibox][i][j][k].vy,
			model[ibox][i][j][k].vz,
			model[ibox][i][j][k].linewidth,
                        model[ibox][i][j][k].mean_av,
			model[ibox][i][j][k].abundance);
		}

	}

       for (ibox=0;ibox<nbox;ibox++){
                i = nx[ibox]/2 ;
                k = nz[ibox]/2 ;
		fprintf(fpnotes,"Box number %d out of %d\n",ibox,nbox);
		fprintf(fpnotes,"Number of cells %d\n",ny[ibox]);
                fprintf(fpnotes,"Check the model at i,k = %d %d x,z = %f %f\n",
                        i,k,xcenter[ibox][j],zcenter[ibox][k]);
		fprintf(fpnotes,"  j      y       T     Tdust  density       ");
                fprintf(fpnotes,"vx         vy         vz     linewidth  mean Av    abundance\n");

                for (j=0;j<ny[ibox];j++){
	fprintf(fpnotes,"%3d %9.6f %6.2f %6.2f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
                        j,ycenter[ibox][j],
                        model[ibox][i][j][k].temperature,
			model[ibox][i][j][k].dust_temperature,
                        model[ibox][i][j][k].density,
                        model[ibox][i][j][k].vx,
                        model[ibox][i][j][k].vy,
                        model[ibox][i][j][k].vz,
                        model[ibox][i][j][k].linewidth,
                        model[ibox][i][j][k].mean_av,
                        model[ibox][i][j][k].abundance);
                }

        }

       for (ibox=0;ibox<nbox;ibox++){
                j = ny[ibox]/2 ;
                k = nz[ibox]/2 ;
		fprintf(fpnotes,"Box number %d out of %d\n",ibox,nbox);
		fprintf(fpnotes,"Number of cells %d\n",nx[ibox]);
                fprintf(fpnotes,"Check the model at j,k = %d %d y,z = %f %f\n",
                        j,k,ycenter[ibox][j],zcenter[ibox][k]);
		fprintf(fpnotes,"  i      x       T     Tdust  density       ");
                fprintf(fpnotes,"vx         vy         vz     linewidth  mean Av    abundance\n");

                for (i=0;i<nx[ibox];i++){
	fprintf(fpnotes,"%3d %9.6f %6.2f %6.2f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
                        i,xcenter[ibox][i],
                        model[ibox][i][j][k].temperature,
			model[ibox][i][j][k].dust_temperature,
                        model[ibox][i][j][k].density,
                        model[ibox][i][j][k].vx,
                        model[ibox][i][j][k].vy,
                        model[ibox][i][j][k].vz,
                        model[ibox][i][j][k].linewidth,
                        model[ibox][i][j][k].mean_av,
                        model[ibox][i][j][k].abundance);
                }

        }
       fflush(fpnotes);
       fflush(stdout);

/* Finished printing the model to the notes file */
/* ------------------------------------------------------------------------------ */

/* Find the minimum and maximum linewidths */
	lwmin_model =  1.e20;
	lwmax_model = -1.e20;
        for (ibox=0;ibox<nbox;ibox++){
        for (j=0;j<ny[ibox];j++){
        for (k=0;k<nz[ibox];k++){
        for (i=0;i<nx[ibox];i++){
	if (model[ibox][i][j][k].linewidth < lwmin_model) 
		lwmin_model = model[ibox][i][j][k].linewidth;
	if (model[ibox][i][j][k].linewidth > lwmax_model) 
		lwmax_model = model[ibox][i][j][k].linewidth;

        }}}}
/*	printf("min max linewidth found in model %f %f\n",lwmin_model,lwmax_model); */

/* Check if the linewidth in the model exceeds the range originally specified. */
  if (lwmin_model < lwmin || lwmax_model > lwmax) {
    printf("linewidths outside originally specified range: %f %f\n",
      lwmin_model/1.e5,lwmax_model/1.e5);
    printf("Model number %d. linewidth min and max %f %f\n",
      model_number,lwmin/1.e5,lwmax/1.e5);
    fprintf(fpnotes,"linewidths outside originally specified range: %f %f\n",
      lwmin_model/1.e5,lwmax_model/1.e5);
    fprintf(fpnotes,"Model number %d. linewidth min and max %f %f\n",
      model_number,lwmin/1.e5,lwmax/1.e5);
    printf("Linewidth interpolation will fail. Must exit now...\n");
    fprintf(fpnotes,"Linewidth interpolation will fail. Must exit now...\n");
       fflush(fpnotes);
       fflush(stdout);
       return 0;
  }

/* Finished checking linewidths of model */
/* ------------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------------ */
	return 1; 
}



unsigned long mfsize(FILE *fp)
{
/*
The LONG_MAX used in function mfsize
//  is a system variable depending on the WORDSIZE
//  of the computers. The function mfsize worked
//  only for LONG_MAX=2147483647L, 
//  the value for 32 bit computers.
//  The LONG_MAX is replaced by LONG_MAX32 which
//  is defined to 2147483647L. 
*/
#define LONG_MAX32 2147483647L

    /* Optimization stuff */
    char temp[BUFSIZ];
    static const long DATALENGTH_MAX=LONG_MAX32%2!=0?LONG_MAX32-1:LONG_MAX32;
    long datalength=DATALENGTH_MAX;

    unsigned long i, counter, fsize;
    fsize = 0;

    if (fp==NULL)
    {
     printf("In mfsize the pointer to the input file is NULL\n");
     return 0;
    }


/* fseek() doesn't signal EOF so i use fread() to detect the end of file */
for (fseek(fp, datalength-1, SEEK_SET); 
    datalength>0 && fread(temp, 1, 1, fp)==0; 
       fseek(fp, datalength-1, SEEK_SET)) datalength/=128;

    fseek(fp, 0, SEEK_SET);
	
    if (datalength==0 && fread(temp, 1, 1, fp)==0)
    {
        return fsize;
    }

    else if (datalength==0)
           datalength=BUFSIZ;

    fseek(fp, datalength-1, SEEK_SET);


    /* fseek() doesn't signal EOF so i use fread() to detect the end of file */
    for(counter=0; fread(temp, 1, 1, fp)!=0; ++counter)
    	fseek(fp, datalength-1, SEEK_CUR);

    fseek(fp, 0, SEEK_SET);

    for( ; counter>0; --counter)
    {
        fseek(fp, datalength, SEEK_CUR);
        fsize += datalength;

    }

    do
    {
        fsize += datalength=fread(temp, 1, BUFSIZ, fp);

    }while(datalength!=0);

    fseek(fp, 0, SEEK_SET);

    return fsize;
}





void chemistry(int depletion, int photodissociation, int molecule, int nbox, int nx[],int ny[], int nz[])
{

float g0 = 0.001;                          /* units of Habing flux */
float cr_ion_rate = 3.e-17;              /* There are several estimates of the cosmic ray ionization
                                           rate: 1.3e-17, 3.0e-17, 1.0e-16 Use the one in the middle. */
float additional_desorption = 3.e-16/cr_ion_rate;

float sticking = 1.0;
float EdCO = 1100.0;                    /* binding energy from Oberg et al 2005 (CO on H2O ice) */
float grain_radius_max = 2.5e-5 ;       /* cm */
float grain_radius_min = 5.0e-7 ;       /* cm */
float grain_cross_section = 3.37e-12 ;  /* cm^2 from Paola's MRN distr. */

float cr_desorb_time, evol_time0, depl_time, v_thermal, x_grain;
float cplus_dest,cplus_creat,co_dest,co_creat,cneut_dest,cneut_creat,cpRcn,coRcp,cnRco,co,cneut,cplus;
double av, av_limit=1.e-6;

int i,j,k,ibox;
double aa,bb,cc,alpha,bbeta,ggama,delta;
double density,tgas,tdust;

/* x_grain is the ratio of the number density of grains to H2 molecules. In KC08, this is
   Rdg = 4 x 10^-10 The number below is 6 x 10^-10 with the grain min and max radii above. */
x_grain = 5.36e-6*pow( (grain_radius_max/1.e-4),(-0.5) ) * pow( (grain_radius_min/1.e-8),(-2.5) );
cr_desorb_time = 3.3e6 / (cr_ion_rate / 1.e-17) * exp( EdCO/ 70.0) / additional_desorption;

double rate[4][4];


/*  Grain cross-section and mass of 2 grams per cm3 is from KC08, just above eqn 10 */
double grainXsection = grain_cross_section;  
double grainRadius = sqrt( grainXsection / PI);
double grainVolume = 4.0*PI/3.0 * pow(grainRadius,3) ;
double grainWeight = 2.0*grainVolume / AMU ;
double grainArea = 4.0*grainXsection;


double fsi = 1.0 ; // fraction of surface covered by ice
double f0 = 1.0e8 ; // photons cm-2 s-1 for G0 = 1

double thermalDesorption[4];

// The vector order is OX , OH , H2O , ICE
double weight[4];
weight[0] = 16.0;
weight[1] = 17.0;
weight[2] = 18.0;
weight[3] = grainWeight;

double adsorptionEnergy[4]; //  from table 1 Hollenbach 2009
adsorptionEnergy[0]  =  800.; 
adsorptionEnergy[1]  = 1300.;
adsorptionEnergy[2] = 4800.;
adsorptionEnergy[3] = 0.0;

double totalDesorption[4];
double oxygen_abundance = 3.2e-4;
double relativeAbundance;
double freezing[4];

double sq,water;

/* ; thermal desorption is equation 2 of H09
; this is the rate per molecule and needs to be multiplied by the number of ice
; molecules. The number of molecules cm-2 on the grain surface 
; is (from H09 below eqn 3)  */
double Nsi = 1.e15 ; // cm-2 

/* ; so the number of molecules released into the gas phase cm-3 is 
; the (rate from eqn2) * Nsi * grain area * number of grains cm-3.
; The number of grains cm-3 is the ratio of grain to gas number density 
; times the gas density. Get this ratio from KC08 just above eqn 10. */
double Rdg = 4.e-10 ;

// ; photo desorption yield is from table 1 of H09
double yield[4];
yield[0] = 0.0;
yield[1] = 2.e-3;
yield[2] = 1.e-3;
yield[3] = 0.0 ;

// ; cosmic ray desorption from H09 table 1
double CRDesorption[4];
double CRDesorption_yield[4];
CRDesorption_yield[0] = 0.0 ; 
CRDesorption_yield[1] = 0.0 ; 
CRDesorption_yield[2] = 4.4e-17 ; 
CRDesorption_yield[2] = 4.4e-17 ; 
CRDesorption_yield[3] = 0.0  ;

double fuv,ggamma,co_undepleted;
double photoDesorption[4],speed[4];
double a,b,c,gg,dd;
int ii,jj;

    for (ibox=0;ibox<nbox;ibox++){
    for (j=0;j<ny[ibox];j++){
    for (k=0;k<nz[ibox];k++){
    for (i=0;i<nx[ibox];i++){

/* The mean_av is not the visual extinction but the transmission
   which is exp(-Av). Therefore in the Hollenbach formulas which
   have factors such as exp(-2.6*Av) we use instead pow(av,2.6)
*/
        av = model[ibox][i][j][k].mean_av;
        if (av < av_limit) av = av_limit;

	density = model[ibox][i][j][k].density;
	tgas    = model[ibox][i][j][k].temperature;
	tdust   = model[ibox][i][j][k].dust_temperature;
	if (tdust <= 0.1) tdust = tgas;

//printf("Starting chemistry ijk %d %d %d\n",i,j,k);
//printf("density tgas tdust av %e %e %e %e \n",density,tgas,tdust,av);

/* The H2O chemistry requires a dust temperature. If the dust temperature is
   not set then it will be zero. In this case, use the gas temperature. 
*/

/* All the molecules except H2O follow the CO model */

if (molecule != 6) {
/* Starting here, all the CO is in the gas phase. We then reduce
   for photodissociation and then for freeze-out if both or either
   are requested */

//printf("molecule photo depletion %d %d %d\n",molecule,photodissociation,depletion);

co = 1.0;

/* This calculation gives us CO in the gas phase in the presence
   of a UV field specified by G0. The units are the fractional
   abundance of total carbon which is a number between 0 and 1. */
      if (photodissociation) {

cplus_dest  =  density * 6.e-16                 ; // cplus
cplus_creat = 2.1e-10 * g0 * pow(av,2.6)       ; // cneut

co_dest  = 1.4e-11 * g0 * pow(av,3.2)          ; // co
co_creat = cplus_dest                           ; // cplus

cneut_dest  = cplus_creat                       ; // cneut
cneut_creat = co_dest                           ; // co


alpha = cplus_dest;
bbeta = cneut_dest;
ggamma = 4.4e-12 * pow((300./tgas),0.61);
delta = co_dest;

aa = ggamma/bbeta * pow((delta/alpha),2);
bb = 1.e0 + delta/alpha + delta/bbeta;
cc = -1.0;

co = (-bb + sqrt(pow(bb,2) - 4.0*aa*cc))/(2.0*aa);
cplus = delta/alpha*co;
cneut = (ggamma*pow((co*delta/alpha),2) + delta*co) / bbeta;

if (av < 1.e-3) {
        cpRcn = cplus_creat / cplus_dest;
        coRcp = co_creat    / co_dest;
        cnRco = cneut_creat / cneut_dest;
        co = coRcp / ( 1. + coRcp + 1./cpRcn);
        cneut = co*cnRco;
        cplus = cneut*cpRcn;
}

// if (j == 32 && k == 32) printf("CO C+ C %d %d %e %e %e %e\n",ibox,i,co,cplus,cneut,av);

      } /* end of IF block for photodissociation */

/* This calculation is for CO freeze out onto grains */
/* Molecules other than H2O follow CO if depletion is set. */
      if (depletion) {
        v_thermal = sqrt( 8.*BOLTZ*tgas /(PI*28.*AMU));
        depl_time = 1. /
            (sticking*x_grain*density*grain_cross_section*v_thermal);
        evol_time0 = depl_time * cr_desorb_time / (depl_time + cr_desorb_time);

//printf("depletion %e\n",evol_time0 / cr_desorb_time);

	co *= evol_time0 / cr_desorb_time;

      }

/* Now we have reduced the fractional abundance of CO by photodissociation and
   depletion if either or both were requested. Apply this fraction to the abundance
   in the model. */

//printf("final CO %e\n",co);

//printf("standard abundance %e\n",model[ibox][i][j][k].abundance);

model[ibox][i][j][k].abundance *= co;

//if (j == 32 && k == 32)
//printf("adjusted abundance %d %d %d %d %e %e %e %e\n",
//ibox,i,j,k,model[ibox][i][j][k].density, model[ibox][i][j][k].temperature,
//model[ibox][i][j][k].abundance, evol_time0/cr_desorb_time);



} else { /* end of IF molecule not equal to H2O */


/* If we got here, then the molecule is water, molecule=6 */

// The velocity of each molecule at the gas temperature
for (ii=0;ii<4;ii++) speed[ii] = sqrt(8.0 * BOLTZ * tgas / (PI * weight[ii]*AMU));

/*
;print,'Total grain cross section per H nucleus ', grainXsection*Rdg

; freeze out from KC08 equation 8 units are molecules cm-3. Need to
; multiply this rate by the total abundance of oxygen and then by
; the fractional abundance of each species. So, totalOxygen*density*fractional abundance
; is the density of oxygen. The fractional abundance will be calculated
; as the solution. Rdg*density is the number density of grains. We end
; up with density^2 and a single power each of Rdg and totalOxygen.
*/
//for (ii=0;ii<4;ii++)
//freezing[ii] = sticking*Rdg*grainXsection*speed[ii]*oxygen_abundance*pow(density,2);

// This is the freezing rate from Hollenbach09 eqn 12, rate=number of dust grains / timescale
for (ii=0;ii<4;ii++)
freezing[ii] = 4.e4 * Rdg * pow(density,2) * pow(tgas/10.0,0.5) * pow(grainRadius,2) * oxygen_abundance;


for (ii=0;ii<4;ii++) {
thermalDesorption[ii] = 1.6e11 * sqrt(adsorptionEnergy[ii]/BOLTZ * AMU/weight[ii])
	* exp(-adsorptionEnergy[ii] / (BOLTZ * tdust)) ;
thermalDesorption[ii] = thermalDesorption[ii] * Nsi * fsi * grainArea * Rdg * density ; //  molecules cm-3
}

/* the units of cosmic ray desorption from table 1 of H09 
   are molecule-1 s-1 so we need to treat this the same
   way as the thermal desorption. 
*/
for (ii=0;ii<4;ii++)
CRDesorption[ii] = CRDesorption_yield[ii] * Nsi * fsi * grainArea * Rdg * density ; // number cm-3

// ; photo desorption is equation 6 of H09
fuv = f0 * 0.001;  // use this for no PD but with CR
fuv = g0 * f0 * (pow(av,1.8) + 0.0005);  // half
fuv = g0 * f0 * (pow(av,1.8) );         // use this for no CR
fuv = g0 * f0 * pow(av,1.8) + f0* 0.0001;  // tenth
fuv = g0 * f0 * pow(av,1.8) + f0* 0.00001;  // hundredth
fuv = g0 * f0 * pow(av,1.8) + f0* 1.e-6  ;  // thousandth

for (ii=0;ii<4;ii++)
photoDesorption[ii] = yield[ii] * fuv * fsi ;


/*
printf("photo  desorption %10.3e %10.3e %10.3e %10.3e\n",
photoDesorption[0],
photoDesorption[1],
photoDesorption[2],
photoDesorption[3]);
*/

/*
; The only units here are the flux f0 in photons cm-2 s-1
; Need to multiply this by the total grain area assuming that photons
; can hit the grain from all directions. Use either the grain area
; or the grain cross-section. Here we use the grain cross-section.
*/

for (ii=0;ii<4;ii++)
photoDesorption[ii] = photoDesorption[ii] * grainXsection * Rdg * density ;

// ; The total desorption is a vector with the rate per cm-3 for each species

for (ii=0;ii<4;ii++)
totalDesorption[ii] = thermalDesorption[ii] + CRDesorption[ii] + photoDesorption[ii] ;

/*
printf("CR desorption %d %10.3e %10.3e %10.3e %10.3e %10.3e\n",
i,av,
CRDesorption[0],
CRDesorption[1],
CRDesorption[2],
CRDesorption[3]);
*/

/*
printf("photo desorption %d %10.3e %10.3e %10.3e %10.3e %10.3e\n",
i,av,
photoDesorption[0],
photoDesorption[1],
photoDesorption[2],
photoDesorption[3]);
*/

/*
printf("i Av %d %e %e %e %e %e\n",i,av,g0,f0,fuv,fsi);
printf("total desorption %d %10.3e %10.3e %10.3e %10.3e %10.3e\n",
i,av,
totalDesorption[0],
totalDesorption[1],
totalDesorption[2],
totalDesorption[3]);
*/

/* ; O + H2 -> OH
; this is rate alpha or K1 */
a = 9.0e-12 ;
b = 1.0e0 ;
c = 4.5e3 ;
aa = oxygen_abundance * pow(density,2) * a * pow((tgas/300.0),b) * exp(-c/tgas) ;

/* ; OH + UV -> O + H
; this is rate beta or P1 */
a = 2.2e-10 ;
b = 2.0e0 ;
bb = oxygen_abundance * density * a * g0 * pow(av,b) ;


/* ; H2O + UV -> OH + H
; this is rate gamma or P2 */
a = 5.1e-10 ;
b = 1.8e0 ;
gg = oxygen_abundance * density * a * g0 * pow(av,b) ;

/* ; OH + H2 -> H2O + H
; this is rate delta or K2 */
a = 3.6e-11 ;
b = 0.0 ;
c = 2.6e3 ;
dd = oxygen_abundance * pow(density,2) * a * pow((tgas/300.0),b) * exp(-c/tgas) ;

/* Now we have all the reactions. Combine these into a rate matrix */

/* ; Creation and destruction of atomic oxygen, OX
; rate[0,0] is oxygen destruction consisting of freezing onto ice 
; and slow molecular reaction with H2 */
rate[0][0] = -(freezing[0] + aa) ;
/* ; rate[1,0] = creation of OX by p.d. of OH */
rate[1][0] = bb ;
rate[2][0] = 0.0 ;
rate[3][0] = 0.0 ;

/* ; Creation and destruction of molecule OH *.
; rate[0,1] is the creation of OH by the molecular reaction H + O -> OH */
rate[0][1] = aa ;
//; rate[1,1] is the destruction of OH by freezing and photodissoc.
rate[1][1] = -(freezing[1] + bb) ;
//; rate[2,1] is the creation of OH by the photodissociation of H2O
rate[2][1] = gg ;
rate[3][1] = 0.0 ;

/* ; Creation and destruction of molecule H2O
; There is no creation of H2O from atomic OX */
rate[0][2] = 0.0 ;
// ; rate[2,1] is the creation of H2O from the slow molecular reaction OH + H2
rate[1][2] = dd ;
// ; rate[2,2] is the destruction of H2O by freezing and photodissociation
rate[2][2] = -(freezing[2] + gg) ;
// ; rate[3,2] is the creation of H2O by desorption off ice
rate[3][2] = totalDesorption[2] ;

/* Creation and destruction of ICE
;rate[0,3] to [2,3] are the creation of ice by freezing OX, OH, and H2O */
rate[0][3] = freezing[0] ;
rate[1][3] = freezing[1] ;
rate[2][3] = freezing[2] ;
// ; rate[3,3] is the destruction ice by photodesorption
rate[3][3] = -(photoDesorption[1] + totalDesorption[2]) ;

// ;replace the last row by the conservation equation
// and multiply all the other rates by 1E12 to get them
// closer to unity.
for (ii=0;ii<4;ii++) {
	rate[ii][3] = 1.00 ;
	for (jj=0;jj<4;jj++) rate[ii][jj] *= 1.00e12 ;
}

/*
; IDL matrix order is (columns,rows)
; A00	A10	A20
; A01	A11	A21
; A02	A12	A22
*/

/* Now I have a 4x4 matrix

|	r00	r10	0	0	|	W		0
|	r01	r11	r21	0	|	X		0
|	0	r12	r22	r32	|	Y	=	0
|	1	1	1	1	|	Z		1

with W = oxygen
     X = OH
     Y = H2O
     Z = ice

write out 4 equations

r00W + r10X = 0
r01W + r11X + r21Y = 0
r12X + r22Y + r32Z = 0
W + X + Y + Z = 1

Solve for Y

W = -r10/r00 X
X = -r21 / sq Y
where sq = (r11 - r01r10/r00)
Z = -1/r32(r22 - r12r21/sq) Y

rewrite as

W = a X
X = b Y
Z = c Y

put these into last (conservation) equation

abY + bY + 1 + cY = 1

y = 1/(ab + b + 1 + c)

*/
aa = -rate[1][0]/rate[0][0] ;
sq = (rate[1][1] - rate[0][1]*rate[1][0]/rate[0][0]) ;
bb = -rate[2][1]/sq ;
cc = -1.0/rate[3][2]*(rate[2][2] - rate[1][2]*rate[2][1]/sq) ;
relativeAbundance = 1.0/(aa*bb + bb + 1.0 + cc) ;

// Experiment to check solution. At very high density the water abundance
// is just this ratio
//if (density > 8.e4) relativeAbundance = totalDesorption[2]/freezing[2];

//printf("relative abund H2O %e \n",relativeAbundance);
//model[ibox][i][j][k].abundance = relativeAbundance ;
model[ibox][i][j][k].abundance = relativeAbundance * oxygen_abundance;


/*
oxy = relativeAbundance[0] ;
oh  = relativeAbundance[1] ;
h2o = relativeAbundance[2] ;
ice = relativeAbundance[3] ;

h2o = h2o * oxygen_abundance ;
ice = ice * oxygen_abundance ;
oh  = oh  * oxygen_abundance ;
oxy = oxy * oxygen_abundance ;
*/

} /* End of the IF block for molecule != 6 */

    }}}} /* end of loop over model cube */



} /* end of chemistry function */


int fft2d (int inverse, int nx, int ny, float **beamR, float **beamI, 
           float **fftR, float **fftI)
{

  int i,j,k;                /* generic index */
  int N;                    /* number of points in FFT */

/*
  double (*normal_var)[2],(*fft_var)[2];
  double (*split_normal_var)[2],(*split_fft_var)[2];
*/
  double **normal_var, **fft_var;
  double **split_normal_var, **split_fft_var;

  N = nx;
  if (ny > nx) N = ny;

/*
  normal_var = malloc(2 * N * sizeof(double));
  fft_var    = malloc(2 * N * sizeof(double));
  split_normal_var = malloc(2 * N * sizeof(double));
  split_fft_var    = malloc(2 * N * sizeof(double));
*/

  normal_var = (double **) malloc(N * sizeof(double *));
  fft_var    = (double **) malloc(N * sizeof(double *));
  split_normal_var = (double **) malloc(N * sizeof(double *));
  split_fft_var    = (double **) malloc(N * sizeof(double *));

  for (i=0;i<N;i++) {
      normal_var[i] = (double *) malloc(2 * sizeof(double));
      fft_var[i]    = (double *) malloc(2 * sizeof(double));
      split_normal_var[i] = (double *) malloc(2 * sizeof(double));
      split_fft_var[i]    = (double *) malloc(2 * sizeof(double));
  }
  

//  printf("Finished allocating memory in FFT2D\n");

/* For each IX, FFT the JY dimension */

     for (i=0;i<nx;i++) {

/* Copy the JY dimension to a 1D vector */
      if (inverse == 0) {
         for (j=0;j<ny;j++) {
                normal_var[j][0] = beamR[i][j];
                normal_var[j][1] = beamI[i][j];
         }
      } else {
         for (j=0;j<ny;j++) {
                normal_var[j][0] = fftR[i][j];
                normal_var[j][1] = fftI[i][j];
         }
      }

//      printf("Finished copying beam on Y axis I = %d\n",i);

/* Split the 1D vector to prepare for FFT */
/* First half of splitting */
         for (j=0;j<ny/2;j++) {
                k = j + ny/2;
//printf("N = %d, j = %d, k = %d\n",N,j,k);
                split_normal_var[j][0] = normal_var[k][0];
                split_normal_var[j][1] = normal_var[k][1];
         }
//      printf("Finished first split Y axis I = %d\n",i);
/* Second half of splitting */
         for (j=ny/2;j<ny;j++) {
                k = j - ny/2;
//printf("N = %d, j = %d, k = %d\n",N,j,k);
                split_normal_var[j][0] = normal_var[k][0];
                split_normal_var[j][1] = normal_var[k][1];
         }
//      printf("Finished second split on Y axis I = %d\n",i);

/*
   if (i == nx/2) {
   printf("split normal \n");
   for (j=0;j<ny;j++) {
       printf("%d %f\n",j,split_normal_var[j][0]);
   } }
*/

//    printf("Starting first FFT\n");
/* The FFT */
	if (inverse == 0) {
           fft(nx,split_normal_var,split_fft_var);
	} else {
           ifft(nx,split_fft_var,split_normal_var);
	}
//    printf("Finished first FFT\n");

/*
   if (i == nx/2) {
   printf("split fft \n");
   for (j=0;j<ny;j++) {
       printf("%d %f\n",j,split_fft_var[j][0]);
   } }
*/

/* Now unsplit the fft variable */
/* First half of splitting */
         for (j=0;j<ny/2;j++) {
		k = j + ny/2;
                fft_var[j][0] = split_fft_var[k][0];
                fft_var[j][1] = split_fft_var[k][1];
         }
/* Second half of splitting */
         for (j=ny/2;j<ny;j++) {
		k = j - ny/2;
                fft_var[j][0] = split_fft_var[k][0];
                fft_var[j][1] = split_fft_var[k][1];
         }

/* Copy 1D vector back to 2D array */
       if (inverse == 0) {
         for (j=0;j<ny;j++) {
		fftR[i][j] = fft_var[j][0];
		fftI[i][j] = fft_var[j][1];
         }
       } else {
         for (j=0;j<ny;j++) {
		beamR[i][j] = fft_var[j][0];
		beamI[i][j] = fft_var[j][1];
         }
       }


/*
   if (i == nx/2) {
   printf("fft \n");
   for (j=0;j<ny;j++) {
       printf("%d %f\n",j,fft_var[j][0]);
   } }
*/

     }

/* Now we have the 2D array FFT,ed on the JY dimension 
   Repeat all of the above for the IX dimension        */


/* For each JY, FFT the IX dimension */

     for (j=0;j<ny;j++) {


/* Copy the IX dimension to a 1D vector */
      if (inverse == 0) {
         for (i=0;i<nx;i++) {
                normal_var[i][0] = fftR[i][j];
                normal_var[i][1] = fftI[i][j];
         }
      } else {
         for (i=0;i<nx;i++) {
                normal_var[i][0] = beamR[i][j];
                normal_var[i][1] = beamI[i][j];
         }
      }

/* Split the 1D vector to prepare for FFT */
/* First half of splitting */
         for (i=0;i<nx/2;i++) {
                k = i + nx/2;
                split_normal_var[i][0] = normal_var[k][0];
                split_normal_var[i][1] = normal_var[k][1];
         }
/* Second half of splitting */
         for (i=nx/2;i<nx;i++) {
                k = i - nx/2;
                split_normal_var[i][0] = normal_var[k][0];
                split_normal_var[i][1] = normal_var[k][1];
         }

/*
   if (j == ny/2) {
   printf("2nd split normal \n");
   for (i=0;i<ny;i++) {
       printf("%d %f\n",i,split_normal_var[i][0]);
   } }
*/

/* The FFT */
	if (inverse == 0) {
           fft(nx,split_normal_var,split_fft_var);
	} else {
           ifft(nx,split_fft_var,split_normal_var);
	}

/*
   if (j == ny/2) {
   printf("2nd split fft \n");
   for (i=0;i<ny;i++) {
       printf("%d %f\n",i,split_fft_var[i][0]);
   } }
*/

/* Now unsplit the fft variable */
/* First half of splitting */
         for (i=0;i<nx/2;i++) {
		k = i + nx/2;
                fft_var[i][0] = split_fft_var[k][0];
                fft_var[i][1] = split_fft_var[k][1];
         }
/* Second half of splitting */
         for (i=nx/2;i<nx;i++) {
		k = i - ny/2;
                fft_var[i][0] = split_fft_var[k][0];
                fft_var[i][1] = split_fft_var[k][1];
         }

/*
   if (j == ny/2) {
   printf("2nd fft \n");
   for (i=0;i<ny;i++) {
       printf("%d %f\n",i,fft_var[i][0]);
   } }
*/

/* Copy 1D vector back to 2D array */
       if (inverse == 0) {
         for (i=0;i<nx;i++) {
		fftR[i][j] = fft_var[i][0];
		fftI[i][j] = fft_var[i][1];
         }
       } else {
         for (i=0;i<nx;i++) {
		beamR[i][j] = fft_var[i][0];
		beamI[i][j] = fft_var[i][1];
         }
       }

     }

  free(normal_var);
  free(   fft_var);
  free(split_normal_var);
  free(   split_fft_var);

  return 1;

}

/*----------------------------------------------------------------------------
   fft.c - fast Fourier transform and its inverse (both recursively)
   Copyright (C) 2004, Jerome R. Breitenbach.  All rights reserved.

   The author gives permission to anyone to freely copy, distribute, and use
   this file, under the following conditions:
      - No changes are made. [Yeah right. If only it worked as written.]
      - No direct commercial advantage is obtained.
      - No liability is attributed to the author for any damages incurred.
  ----------------------------------------------------------------------------*/

/******************************************************************************
 * This file defines a C function fft that, by calling another function       *
 * fft_rec (also defined), calculates an FFT recursively.  Usage:             *
 *   fft(N, x, X);                                                            *
 * Parameters:                                                                *
 *   N: number of points in FFT (must equal 2^n for some integer n >= 1)      *
 *   x: pointer to N time-domain samples given in rectangular form (Re x,     *
 *      Im x)                                                                 *
 *   X: pointer to N frequency-domain samples calculated in rectangular form  *
 *      (Re X, Im X)                                                          *
 * Similarly, a function ifft with the same parameters is defined that        *
 * calculates an inverse FFT (IFFT) recursively.  Usage:                      *
 *   ifft(N, x, X);                                                           *
 * Here, N and X are given, and x is calculated.                              *
 ******************************************************************************/

/* macros */
#define TWO_PI (6.2831853071795864769252867665590057683943L)

/* function prototypes */
/*
void fft(int N, double (*x)[2], double (*X)[2]);
void fft_rec(int N, int offset, int delta,
             double (*x)[2], double (*X)[2], double (*XX)[2]);
void ifft(int N, double (*x)[2], double (*X)[2]);
*/

/* FFT */
void fft(int N, double **x, double **X)
{
  /* Declare a pointer to scratch space. */
  double **XX; 
  int i;
  XX = (double **) malloc(N * sizeof(double *));
  for (i=0;i<N;i++) XX[i] = (double *) malloc(2 * sizeof(double));

  /* Calculate FFT by a recursion. */
  fft_rec(N, 0, 1, x, X, XX);

  /* Free memory. */
  free(XX);
}

/* FFT recursion */
void fft_rec(int N, int offset, int delta,
             double **x, double **X, double **XX)
{
  int N2 = N/2;            /* half the number of points in FFT */
  int k;                   /* generic index */
  double cs, sn;           /* cosine and sine */
  int k00, k01, k10, k11;  /* indices for butterflies */
  double tmp0, tmp1;       /* temporary storage */

  if(N != 2)  /* Perform recursive step. */
    {
      /* Calculate two (N/2)-point DFT's. */
      fft_rec(N2, offset, 2*delta, x, XX, X);
      fft_rec(N2, offset+delta, 2*delta, x, XX, X);

      /* Combine the two (N/2)-point DFT's into one N-point DFT. */
      for(k=0; k<N2; k++)
        {
          k00 = offset + k*delta;    k01 = k00 + N2*delta;
          k10 = offset + 2*k*delta;  k11 = k10 + delta;
          cs = cos(TWO_PI*k/(double)N); sn = sin(TWO_PI*k/(double)N);
          tmp0 = cs * XX[k11][0] + sn * XX[k11][1];
          tmp1 = cs * XX[k11][1] - sn * XX[k11][0];
          X[k01][0] = XX[k10][0] - tmp0;
          X[k01][1] = XX[k10][1] - tmp1;
          X[k00][0] = XX[k10][0] + tmp0;
          X[k00][1] = XX[k10][1] + tmp1;
        }
    }
  else  /* Perform 2-point DFT. */
    {
      k00 = offset; k01 = k00 + delta;
      X[k01][0] = x[k00][0] - x[k01][0];
      X[k01][1] = x[k00][1] - x[k01][1];
      X[k00][0] = x[k00][0] + x[k01][0];
      X[k00][1] = x[k00][1] + x[k01][1];
    }
}

/* IFFT */
void ifft(int N, double **x, double **X)
{
  int N2 = N/2;       /* half the number of points in IFFT */
  int i;              /* generic index */
  double tmp0, tmp1;  /* temporary storage */

  /* Calculate IFFT via reciprocity property of DFT. */
  fft(N, X, x);
  x[0][0] = x[0][0]/N;    x[0][1] = x[0][1]/N;
  x[N2][0] = x[N2][0]/N;  x[N2][1] = x[N2][1]/N;
  for(i=1; i<N2; i++)
    {
      tmp0 = x[i][0]/N;       tmp1 = x[i][1]/N;
      x[i][0] = x[N-i][0]/N;  x[i][1] = x[N-i][1]/N;
      x[N-i][0] = tmp0;       x[N-i][1] = tmp1;
    }
}
