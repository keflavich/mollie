/* Setup file for radiative transfer simulation. This
   file is part of a C program and C syntax is required. */

/* Pathname to the directory where most of the data will be 
 * stored. This directory should be accessible 
 * to the machine where the master is running. 
 * Other files are written by the slaves in the
 * directory where the executable resides. This directory
 * must be cross-mounted to all the machines that are
 * running the slaves. At least "./" is required. */

/*        strcpy(pathname,"/data2/rest_3/keto/L1544/R24/data"); */

        strcpy(pathname,"./");


/* Select how many iterations of lambda iteration. For 
  molecular cloud cores of moderate optical depth 10 - 100
  try 10 to 30 iterations. 

  ITERATIONS and NUMBER_MODELS have been updated (7-10-13)
  to work together to write out intermediate results to check
  the convergence of the Lambda iteration. 
  For example,
	iterations = 10;
	number_models = 10;
  will run a total of 100 iterations and write output
  every 10 iterations. The output files are named
  ModelCube0, ModelCube1, ...
  In this example, the final model, ModelCube10, contains
  the last iteration. If the solution is converged, the
  solution in ModelCube9 should be same as in 10.
  The variable "number_models" is set at the bottom of this
  file (setup.c).
*/

/* NH3 and CH3CN are full LTE only. No point in more than
   one iteration for these 2 molecules. */

	iterations = 4;

/* Set the initial radiation field to either CMB or ZERO or LTE.  
   ZERO means initial Jbar = 0 and I0 = 0
   CMB  means initial Jbar = CMB and I0 = CMB
   LTE  means initial Jbar = Planck function and I0 = CMB */

	initial_radiation = CMB;

/* ----------------------------------------------------*/

/* Select the molecule for the simulation.
Choices:

	H2Oortho	checked OK with accel 1-31-13
	CO		checked OK with accel 1-8-09
	C13O		checked OK with accel 1-7-09
	C17O		checked OK with accel 1-8-09
	C18O		checked OK with accel 1-7-09
	CS		checked OK with accel 1-8-09
	N2Hplus		checked OK with accel 1-7-09
	N2Dplus		checked OK with accel 1-8-09
	HCOplus		checked OK with accel 1-19-09
	H13COplus	checked OK with accel 1-19-09
	SiO		checked OK with accel 1-19-09
	HCN		checked OK with accel 1-19-09

All the molecules above have either no hyperfine structure or 
use the approximation of hyperfine statistical equilibrium (LTE hyperfines)

N2Hhyp below is for N2H+ with non-LTE hyperfines. N2Hplus above
is N2H+ with the LTE hyperfine approximation
	N2Hhyp		checked OK with accel 1-7-09
	HCNhyp
	HCNvib

NH3 and CH3CN are LTE for both the rotational and hyperfine transitions
	NH3		checked OK LTE 1-9-09
	CH3CN  		checked OK LTE 1-8-09

*/

	molecule = N2Hhyp; 

/* Lambda acceleration works most of the time. 
   Tests with HCO+ and H13CO+ failed in that the 
   computed source function was much higher than 
   the gas temperature. Other molecules generally OK. 
   Not relevant for LTE molecules NH3 and CH3CN. */

	acceleration = YES;

	dust_continuum = NO;

/* Photodissociation and depletion are applicable to 
   clouds such as starless cores */

/* Normally we photodissociate all the molecules at the edge
   of the cloud */
	photodissociation = YES ;

/* Normally we deplete all the C molecules, but not the N molecules */
	depletion = YES;
	if (molecule == N2Hplus ||
	    molecule == N2Hhyp  ||
	    molecule == N2Dplus ||
	    molecule == NH3 ) depletion = NO;

/* Some standard abundances for some molecules */
/* These abundances are not automatically used. Depends
   on how define_model.c is written. Also set_model_parameters
   may overwrite these abundances depending on how it is 
   written */
	if (molecule == CO)	std_abundance = 5.625e-5;
	if (molecule == C13O)	std_abundance = 5.625e-5 / 77.;
	if (molecule == C18O)	std_abundance = 5.625e-5 / 77. / 7.3;
        if (molecule == C18O)   std_abundance = 5.625e-5 / 77. / 7.3 * 1.8 ;
        if (molecule == C17O)   std_abundance = 5.625e-5 / 77. / 7.3 / 4.0 * 1.8;
/*	if (molecule == CS)	std_abundance = 3.e-9; */
	if (molecule == CS)	std_abundance = 5.e-10;  /* use this for C34S */
	if (molecule == N2Hplus)std_abundance = 3.e-10;
	if (molecule == N2Hhyp )std_abundance = 2.e-10;
	if (molecule == N2Dplus)std_abundance = 3.e-10;
	if (molecule == NH3)    std_abundance = 0.5*2.4e-8;
        if (molecule == HCN)    std_abundance = 1.5e-8;
        if (molecule == HCNvib) std_abundance = 1.5e-8;
        if (molecule == CH3CN)  std_abundance = 6.3e-10;

/* HCN abundance,
        Irvine, Goldsmith, & Hjalmarson 1987
        Paglione et al 1988
        1.5 x 10-8

   CH3CN abundance from Zhang's paper on AFGL5192
*/


/* Set the number of lines and number of levels. For linear rotors
   this is straight forward, and the number of lines should be 
   less than the number of levels. The lines and levels also have
   to be set in the 2 files nlines_C. and nlines_F77.h

   At the cold temperatures of molecular clouds, the upper states
   have very low populations. In this case, only a few states and
   lines are necessary.

*/

/* There are some special cases that require a specific combination of
   states and lines.

	NH3: can be run with up to 6 lines. The number of states 
	must be twice the number of lines. 
        For example,
        n_state=12	n_lines=6
	produces NH3(1,1) through NH3(6,6)
	nstate=2	n_lines=1
	produces only NH3(1,1)

-------

	N2Hhyp can be run with one of the following pairs

	n_state = 10;
	n_lines = 1;

	n_state = 19;
	n_lines = 2;

	n_state = 28;
	n_lines = 3;

*/
	n_state = 37;
	n_lines = 4;


/*
	n_state = 64;
	n_lines = 7;

*/

//-------

//	CH3CN: requires 	nstate = 120	nlines = 105

/* This is for Rainer's HCN
	n_state = 62;
	n_lines = 125;
*/

//-------

/* This is for CH3CN 
	n_state = 120;
	n_lines = 105;
*/

//-------

/* A nice pair for rotors */
/*
	n_state = 8;
	n_lines = 7;
*/

//-------

/* H2Oortho can be run with 8 states and 11 lines 
   There are only collision rates for 8 states.

These are the combinations that can be used for H2Oortho.
   8 states and 11 lines
   7 states and 9 lines
   6 states and 7 lines
   5 states and 5 lines
   4 states and 4 lines
   3 states and 2 lines
   2 states and 1 lines
The upper states of H2O have high excitation energies and
are probably not populated in cold clouds. 

	n_state = 3;
	n_lines = 2;

*/



/* Set the channel width in cm/s */     
/* Set the velocity range in cm/s. The velocity range has to
   be large enough to cover the linewidth and range of gas
   velocities which are both set in define_model.c.
   If the molecule has hyperfine structure, do not make
   velrange larger to span the hyperfine lines. velrange
   should be set for a single line regardless of the
   hyperfine structure. Mollie knows about hyperfine lines
   and will set the total velocity range according to
   the molecule. */    

	chanwd = 0.010e5;
	velrange = 1.50e5;

/* Choose Hanning smoothing YES or NO. This only affects the
   output. You can always Hanning smooth the output yourself
   with another program. */          
        hanning = YES ;

/* Set the minimum linewidth to decide which hyperfine lines to compute. */
/* Set the maximum linewidth for computing a range of line profiles. */
/* These lwmin,lwmax numbers are the gaussian widths (sigma) of the lines.
   This is a code related parameter that I hope to eliminate from the
   input in the future. For the moment, it has to be set. If the gas
   temperature is 10K, then the sound speed of the molecules in Mollie 
   is larger than 0.1 kms. So lwmin = 3000 cm/s would be a good 
   lower limit. The maximum, lwmax, has to accomodate the maximum
   linewidth set in define_model.c This usually includes both the 
   thermal width and turbulent width. 

   Mollie uses the minimum width here to combine closely spaced hyperfine
   lines in the radiative transfer calculation. for example, N2H+(3-2) 
   has over 40 hyperfine lines, some of which are very close together.
   In the approximation of hyperfine statistical equilibrium, the lines
   that are closer together than lwmin are combined. Another way to
   do this would be to use the channel width as the criterion.

   When calculating non-LTE hyperfine intensities, the hyperfine lines
   are never combined because the overlap is critical to the development
   of the hyperfine anomalies in intensity.

*/
        lwmin =  3000.;
        lwmax =  45000.;

/* ----------------------------------------------------*/

/* Number of nested boxes. All the boxes are centered on the origin. */

	nbox = 3;

/* For each box set the dimensions.  */

/* 
   ny and nz can be as small as 1, but nx must be > 3. 

   If you want to use the beam convolution feature, then
   ny and nz must be a power of 2 because the convolution
   is done by FFT. If beamx and beamy = 0 (below), then no
   convolution is done and the grid does not have to be a 
   power of 2.

   The gridding has been updated (7-10-13) to allow for
   different dimensions in NX and NY,NZ. The line of
   sight is NX and the map or image plane is NY-NZ. 
   So you could make a model that is very long in the
   X-direction and contains just enough Y,Z cells to
   allow for a beam convolution. This would not calculate
   a proper mean intensity, J-bar, but could be useful
   for testing.
*/

        ngrid = 32;

        ibox = 0;
        nx[ibox] = ngrid;
        ny[ibox] = ngrid;
        nz[ibox] = ngrid;


        ibox = 1;
        nx[ibox] = ngrid;
        ny[ibox] = ngrid;
        nz[ibox] = ngrid;

        ibox = 2;
        nx[ibox] = ngrid;
        ny[ibox] = ngrid;
        nz[ibox] = ngrid;



/* Radius of model sphere in pc. The calculations are done inside
  a sphere which can be smaller or larger than the input grid.

  Normally the diameter of the model sphere (2*radius) would 
  be the same as the width of the grid.

  If the model sphere is the same size or smaller than the grid,
  the pathlengths will start and stop on the boundary of the sphere
  rather than the boundary of the grid.  Also the pixels in the 
  corners of the model grid will have no rays through them
  Suppose nx = 100 and the cellsize is 0.002 (pc). Then the
  grid will be 0.200 pc across. If the radius of the sphere
  is 0.100 pc, then it just fits inside the grid. If this is 
  set up as follows then it is easy to change nx, small for
  testing and larger for the final runs.
	radius = 0.100 ;
        cellsize[ibox] = 2.*radius/nx[ibox];

  If the sphere is large enough to contain the entire model grid 
  all the pathlengths (on angles aligned parallel to the grid) 
  will be the same.  All pixels, including those in the corners 
  will have rays through them.
	radius = 0.200 ;
        cellsize[ibox] = radius/nx[ibox];


*/      

	radius = 0.320 ;

/*
  There are now two sizes for each cell. Each box or subgrid
  also has its own sizes. 
  For example, suppose there is one box, no subgrid. Mollie
  can set up a calculation with 8192 cells along X and 24
  cells in Y and Z. The cells can have different sizes along
  X and along the Y,Z. For example,
	cellsizex = radius/nx[0];
	cellsizeyz = radius/ny[0];

  I have not thought of any reason for the cells to have different
  sizes in Y and Z separately. So there are only two cellsizes,
  one for X and another for Y,Z.
*/

	ibox = 0;
        cellsizex[ibox] = 2.0*radius/nx[ibox];

        ibox = 1;
        cellsizex[ibox] = 1.0*radius/nx[ibox];

        ibox = 2;
        cellsizex[ibox] = 0.5*radius/nx[ibox];

/* Make the YZ cell size the same as X. You almost always want to do this */
        for (ibox=0;ibox<nbox;ibox++) cellsizeyz[ibox] = cellsizex[ibox];

        ibox = 0;



/* Set the number of angles through the cloud. These are the
   angles along which radiative transfer will be done. 
*/

/*
   These are rotation angles and the rotation by
   latitude is applied first. The order makes a
   difference. For example, 90 lat followed 0 long
   is a view down the north pole along the z axis
   with x,y as the plane. */

/* The commonly used "6 ray approximation" would be
	nlng = 2;
	nlat = 2;
   The 6-ray approx. puts forward and backward rays 
   along the 3 axes of the Cartesian grid.

   For the LTE molecules, NH3 and CH3CN, use
	nlng = 1;
	nlat = 1;
   because the average intensity J-bar is not used in
   the LTE approximation.
   
*/
       nlng = 2;
       nlat = 2;

/* ----------------------------------------------------*/
/* ----------------------------------------------------*/

/* Output parameters */

/* Set the number of lines to view in number_lines.
   number_lines must be 1 or more.
   Then select the lines. First line to view is
   linelist[0], next line (if number_lines > 1) is
   linelist[1]. For dipole molecules the (1-0)
   transition is number 0, the (2-1) transition is
   number 1, etc. */

        number_lines = 2;	/* output lines to view */

/*
	for (i=0;i<number_lines-1;i++){
        linelist[i] = i;        
	}
*/

	linelist[0] = 0;
	linelist[1] = 2;

 /* For the CH3CN(12,0 -> 11,0) you need the following
    specific lines to be written. */

/*
        number_lines = 12;      
        for (i=0;i<number_lines;i++){
                linelist[i] = 67+i-1;
        }
*/


/* 
   These parameters follow the style of number_lines 
   and linelist above.  The numbers to be entered here 
   are the degrees of the viewing angles rather than
   the array index. There is a section in cmain.c that 
   picks the closest grid angle to each output angle 
   specified here.*/


        outviews = 1;
          outlng[0] =  0. ;
          outlat[0] =  0. ;


/* Set the parameters for the output grid */


/* Output axes are image-like (X and Y are axes on the plane of the sky), 
   and there is no Z axis. This is different from the model grid which
   is physics-like (x and y are the equatorial plane and z points up.)

   Mollie can only write out one grid regardless of the number of boxes
   each of which could have its own cell size. If you choose the finest
   grid, the output intensities in the larger cells will be the same
   within each larger cell. If the output cell is larger than the 
   computational cell, the computed intensities are averaged together.
*/

        noutx = nouty = ngrid;

        cellx = cellsizeyz[0];
        celly = cellsizeyz[0];


/* beamx and beamy are the gaussian widths (sigma) used to weight 
   the output data. The relation between beam and FWHM is
   beam = FWHM/(2.sqrt(ln(2)) The beam size is in pc. */

/* ; N2H+ beam is from 37m Haystack, 27" FWHM
   ; At 150 pc,
   ; sigma = 150.*tan(27./3600*PI/180.)/2.35 pc
   ; sigma = 0.00835530
*/
	beamx = 0.00835530;
        beamy = beamx ;


/* ----------------------------------------------------*/
/* ----------------------------------------------------*/

/* Comparison can be either  NO,  ANNEAL, or  SIMPLEX, 
 * Simplex function is not working. YES is equivalent to
 * ANNEAL */

	comparison = NO;

/* Set the number of models to compute. 

   number_models can be different than 1 in two situations.

   First: The code is set up to read data and compare model spectra against the data. 
          Each output ModelCube corresponds to a new physical model which is defined
          by the search algorithm. comparison = YES

   Second: We want to output intermediate results to check the convergence. See the
           notes at the start of this file (setup.c). If comparison = NO then Mollie
           writes out a new ModelCube after the number of iterations specified at
           the start of this file. The total number of iterations will be
           iterations * number_models. The total number of ModelCubes to produce
           is number_models.
*/

        number_models = 1;

/* The number of adjustable parameters for each model. This is
   irrelevant unless the code is set up for comparisons. Has to be
   set up with the function set_model(). */

	nParameters = 6;
 
