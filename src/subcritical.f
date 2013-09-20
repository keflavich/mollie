        function subcritical(iter, igrid, jgrid, kgrid, 
     &    tk, tdust, h2, abundance, widthline, trans, bJ, arg_frac)

        include 'nlines_f77.h'

        parameter (maxhyp=50)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline,
     &  frac_up,frac_lo, arg_frac,balance

        double precision bracket1,bracket2,accuracy,f2,df
        real *4 tk, h2, abundance
        real *4 fnh3,atoms,brot,dipole
        logical docont,doprint,doprint2,nh3tau
        double precision rtbis,apops

        common /procid/ myid
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /frac/ frac(nstate),xfrac(nstate),gg(nstate)
        common /coll/ cr(nstate,nstate)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /gridcom/ igridcom,jgridcom,kgridcom
        common /cmb/ ycmb(nlines)

        dimension arg_frac(nstate)
        dimension bJ(maxhyp,nlines)
        double precision params(10)

c Get the collision rates. Units of cm^3 s^-1 at temperature Tk
c Multiply by H2 to get s^-1 and by n2 to get cm^-3 s^-1

        doprint = .false.
        call colrateh2o(tk, doprint)

        goto 100

        bracket1 = 0.0
        bracket2 = 0.1
        accuracy = 1.e-5

c These are the variable to pass through rtbis to the function
c apops that calculates the equilibrium between upward and
c downward rates.
        lun = 6
        params(1) = dble(h2)
        params(2) = dble(tdust)
        params(3) = dble(abundance)
        if (igrid .eq. 31 .and. jgrid .eq. 32 .and. kgrid .eq. 32) 
     &      doprint = .true.
        if (doprint) then
            params(5) = 6.d0
        else
            params(5) = 0.d0
        endif
        params(6) = 0.d0
        params(7) = dble(widthline)
        params(8) = dble(bJ(1,1))
        params(9) = dble(trans)

c The next 2 blocks print some data for file 33
        if ( jgrid .eq. 31 .and. kgrid .eq. 32 .and. 
     &       igrid .lt. 33) then

            params(6) = 1.d0
            if (igrid .eq. -1) then
                params(5) = 33.d0
                doprint = .true.
                lun = 33
            endif

            write(33,1010) igrid,jgrid,kgrid,h2,tdust,abundance
 1010       format(3i5,1p,3e12.3)

c This loop does not do anything other than write to file 33
            nguesses = 100
            df = 1.d-2/nguesses
            do 20 i = 1,nguesses
               f2 = df*i
               balance = apops(f2,params)               
               write(33,1009) i,f2,balance
 1009          format(i5,f10.6,1p4e12.3)
 20         continue
            params(6) = 0.d0

        endif
c end of the if block to print to file 33


        if (doprint) write(lun,1001) igrid,jgrid,kgrid,
     &    h2,tk,tdust,abundance
 1001   format('ijk H2 Td Ab dx ',3i5,1p,4e12.4)

        frac(2) = rtbis(bracket1,bracket2,accuracy,params,
     &                  istatus)
        frac(1) = 1.d0 - frac(2)

        if (doprint) write(lun,1004) igrid,jgrid,kgrid,
     &    frac(1),frac(2),istatus
 1004   format('ijk f1 f2 status ',3i5,1p2e12.3,i5)

        arg_frac(1) = frac(1)
        arg_frac(2) = frac(2)

        return
 
 100    continue 

c Below is the two-step approximation. First iteration calculates
c the level pops with the CMB, the next step calculates the level
c pops with the radiation from the first step.


c A21 n2 has units of cm-3 s-1
c Assume that almost all of the molecules are in the ground state.
c This means that to a first approximation n1 = h2*abundance.
c n2 is found from the approximation that every collision results in
c a downward radiative transition because at sub-critical density,
c the collision rates are slower than radiative.
c The number abundance in cm-3 is n2 = cr(1,2)*h2 * h2*abundance / a(2,1)
c Since the total number of molecules is approximately h2*abundance,
c the fractional abundance is f2 = c(1,2)*h2 / a(2,1)
c Finally, set the fractional abundance of the lower level so that
c the population in both levels adds to one.


c This is the CMB. Use this for the first iteration. YCMB is in K.
        if (iter .eq. 1) then

            barJ = ycmb(1) 
            barJ = barJ * boltz/planck/freqline(1)

        else

c On the 2nd iteration, the CMB is already included.
c bJ is indexed for (hyperfine,line)
            barJ = bJ(1,1) * boltz/planck/freqline(1)

        endif

        upJ = statdg(1)/statdg(2)*a(2,1)*barJ
        dnJ = a(2,1)*barJ

        if (jgrid .eq. 32 .and. kgrid .eq. 32)
     &     write(6,1000) igrid,iter,
     &     ycmb(1), barJ,upJ,dnJ, cr(1,2)*h2  , a(2,1)
 1000   format('subcritical ',2i5, 6e12.4)

        frac(2) = ( upJ + cr(1,2)*h2 ) / (a(2,1) + dnJ)
        frac(1) = 1.d0 - frac(2)


        arg_frac(1) = frac(1)
        arg_frac(2) = frac(2)

c The returned variable is the status flag. 

        subcritical = 1

c        if (nstate .ne. 2) subcritical = 0

        return


        end


        function apops(f2,params)

c PARAMS is an array that holds whatever information is necessary
c for the calculation that cannot be obtained from common. Info
c in the C program has to be passed through arguments to get into
c the F77 program

c f2 is the guess of the fractional population of the upper state.
c With this guess and the information in PARAMS, this function
c returns the difference between upward and downward transitions.
c Equilibrium balance is zero.

        include 'nlines_f77.h'

        parameter (maxhyp=50)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline,
     &  frac_up,frac_lo, arg_frac

        real  tk, h2
        real  fnh3,atoms,brot,dipole
        logical docont,doprint,doprint2,nh3tau
c Almost every variable in the calculation is made double prec.
        double precision f1,f2,b12,b21,cont_opac,cont_emis,
     &    line_opac,line_emis,density,tdust,abundance,pathlength,
     &    tau,source,jbar,upJ,dnJ,apops,widthline,coldensity

        common /procid/ myid
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /frac/ frac(nstate),xfrac(nstate),gg(nstate)
        common /coll/ cr(nstate,nstate)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /gridcom/ igridcom,jgridcom,kgridcom
        common /cmb/ ycmb(nlines)

        double precision params(10)

c Information passed through params. PARAMS(4) is not used.
c 5 and 6 are used to control the diagnostic print statements
        density = params(1)
        tdust = params(2)
        abundance = params(3)
        widthline = params(7)
        jbar = params(8)
        transmission = params(9)

        doprint = .false.
        if (params(5) .gt. 0.5) then
            doprint = .true.
            lun = nint(params(5))
        endif

        f1 = 1.d0 - f2

c Dust emission and opacity using the same calculation in cmain.c
        cont_opac = 4.0*3.3e-26 * (freq(2,1)/2.83e11)**2.0
        cont_emis = 2.0*PLANCK * (freq(2,1)/C)**2 * freq(2,1)
     &   / ( exp( PLANCK*freq(2,1) /(BOLTZ*tdust)) -1.0)
     &   * cont_opac

c Einstein coefficients
        b21 = a(2,1) *c /freq(2,1)*c/freq(2,1)/(2.d0*planck*freq(2,1))
        b12 = b21 * statdg(2)/statdg(1)

c Convert the linewidth from cm/s to GHz
        deltanu = widthline * freq(2,1) / c

c The line opacity including the normalization by the linewidth.
c This normalization is an approximation for the integration over
c the line profile function
        line_opac = ( f1*b12 - f2*b21) 
     &    * planck*freq(2,1)/(4.*pi) * abundance / (sqrtpi * deltanu)

        if (doprint) write(lun,1007) f1,f2,b12,b21
 1007   format('f1*b12 f2*b21 diff',1p4e12.3)
        if (doprint) write(lun,1008) f1*b12,f2*b21,f1*b12-f2*b21,
     &    planck*freq(2,1)/(4.*pi) * abundance,
     &    line_opac
 1008   format('* * opac ',1p6e12.3)

c The line emissivity
        line_emis = f2*a(2,1) * planck*freq(2,1) /(4.*pi) * abundance
     &              / (sqrtpi *deltanu)

        if (doprint) write(lun,1005) 
     &    cont_opac,cont_emis,line_opac,line_emis
 1005   format('cont op em line op em ',1p,4e12.3)

c The C program calculates the mean transmission (named mean_av
c in the C program) for each cell. The mean transmission is
c derived from a calculation of the column density. Here we reverse
c this calculation to get the mean column density around a cell and
c finally the mean optical depth surrounding each cell.
        coldensity = -log(transmission) * 9.4d20
        tau = coldensity*(line_opac + cont_opac)
c This includes only the continuum opacity in the optical depth
        tau = coldensity*( cont_opac)

c        tau = 1.d0
c        pathlength = tau / density / (line_opac + cont+opac)
        if (doprint) write(lun,1010) coldensity,tau
 1010   format('columnDensity tau ',1p,2e12.3)

c The units of the source function are degrees K
c        source = (line_emis + cont_emis) / (line_opac + cont_opac) 
c This includes only the continuum emission in the source function
c and Jbar.
        source = ( cont_emis) / ( cont_opac) 
     &           * ((c/freq(2,1))**2) / (2.*boltz)

        if (doprint) write(lun,1002) ,source,ycmb(1)
 1002   format('source ycmb ',1p,2e12.3)

c Approximate Jbar by the source function modified by the
c optical depth according to the radiative transfer eqn.
c The CMB is in units of degrees K
        jbar = source*(1.d0 - dexp(-tau)) + dble(ycmb(1))

c Convert the radiation into non-dimensional units
        jbar = jbar * boltz/planck/freqline(1)

c With the radiation in non-dimensional units we use the
c Einstein A instead of B J(energy)*B21 = J(non-dim)*A21
c Upward and downward Einstein B radiative transitions
        upJ = statdg(1)/statdg(2)*a(2,1)*jbar
        dnJ = a(2,1)*jbar

        if (params(6) .gt. 0.5) 
     &     write(33,1003),upj,dnJ,cr(1,2)*density,a(2,1)
        if (doprint) then 
            write(lun,1003),upj,dnJ,cr(1,2)*density,a(2,1)
            write(lun,1012),b12,b21,jbar/boltz*planck*freqline(1)
        endif
 1012   format('b12 b21 Jbar ',1p,3e12.3)
 1003   format('upJ dnJ cup A ',1p,4e12.3)

c Return the equilibrium condition, transitions up - down. Upward
c includes the Einstein B and the collisions. Downward, the Einstein
c A and B
        apops = ( upJ + cr(1,2)*density )*f1 - (a(2,1) + dnJ)*f2
c        apops = ( upJ + cr(1,2)*density )*f1 - (a(2,1) )*f2
        if (doprint) write(lun,1006) f1,f2,apops
 1006   format('f1 f2 apops ',1p,3e12.3)

        return
        end

      double precision  FUNCTION RTBIS
     &  (X1,X2,XACC,params,istatus)
c parms is a an array that is passed through to
c the function FUNC. It can be any size. The
c declared size of 1 here is enough to pass the
c pointer to the first element.
      implicit double precision (a-h,o-z)
      PARAMETER (JMAX=40)
      dimension params(5)
      istatus = 1
      FMID=apops(X2,params)
      F=apops(X1,params)
      if (params(5) .gt. 0.5) 
     &  write(lun,1000) x1,x2,f,fmid
 1000   format('rtbis x1 x2 f1 f2 ',1p,4e12.3)
      IF(F*FMID.GE.0.) then
         istatus = 0
         return
      endif
      IF(F.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5
        XMID=RTBIS+DX
        FMID=apops(XMID,params)
        IF(FMID.LT.0.)RTBIS=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
11    CONTINUE
      istatus = 2
      return
      END
