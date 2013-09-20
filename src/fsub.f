c Version 24 Jan 2013. 

c Corrections:

c Apr 5, 2013
c The upward rates are now calculated from the downward
c rates after the latter are interpolated to the gas
c temperature. Previously, both the upward and downward
c rates were calculated by interpolation. Small deviations
c from detailed balance can prevent convergence at high
c optical depth

c April 2, 2013
c The LTE ammonia now includes the magnetic hyperfine 
c splitting for the (1,1), (2,2), and (3,3) lines.
c The magnetic splitting in the (3,3) line will already
c be difficult to observe so there is no reason to go
c to higher rotational transitions.

c Feb 15, 2013
c Some of the collision rates were calculated by linearly
c interpolating both the upward and downward rates between
c temperatures. It is possible that the rates could be
c out of detailed balance and cause the solution to drift
c at very high optical depth. Now interpolate the downward
c rates and set the upward by detailed balance.

c Jan 24, 2013
c Simplified the explanation of the units in the rate
c equations. No difference in the results, just easier
c to read the code.

c Jan 11, 2013
c This version carries the normalization of the line profile
c function into the ALI equations. Changes to functions
c acceleratedlambda, and newrt.

c Jul 31, 2012
c In NEWRT, changed diagnostic output file from 22 to lun22
c which is myid+22. Prevents the output printed for one specific
c ray from being overwritten by a different slave.

c in subcritical.f defined rtbis, apops, and balance as *8 for
c compatibility with latest GNU compilers.

c Nov 1, 2011
c Corrected more problems in fsub.f where the code checks
c the rotation of the rays. The check did not allow for
c the fact that the projected size of the model box is
c smaller when viewed from the perspective of a rotated 
c ray. The check code has caused more problems than
c code it is supposed to check. Deleted.

c July 27, 2011 Corrected typo in dmatt, gt corrected to .gt.
c This is a print subroutine, not normally used.

c May 4, 2011
c Corrected gridding for the case where the
c model sphere is much larger than the model box.

c 3 March 2010 
c Corrected the geometry to fix the problem that 
c the number of angles in the C and F77 parts of 
c the code were not the same when the view down 
c the north pole was included. 

c23456789112345678921234567893123456789412345678951234567896123456789712
c NSTATES refers to quantum levels of a molecule. If the molecule has
c non-LTE hyperfine structure, the NSTATES includes the hyperfine 
c levels. If the calculation is set up for LTE hyperfine structure, 
c then NSTATES refers only to the rotational levels.

c LINES always refers to the delta J rotational lines. Because
c of line overlaps between the hyperfine lines, we never calculate
c the hyperfine lines alone, but always group them into lines between
c rotational levels. This applies for all the calculations whether
c the hyperfine line brightness is calculated allowing for
c non-LTE hyperfine effects or the hyperfine brightness is assumed 
c to be LTE

c TRANSITIONS (maxtrans) refers to all the allowed radiative 
c transitions (lines) between individual hyperfine levels. For 
c molecules without hyperfine structure or if we assume LTE 
c hyperfine structure, these are never used. The hyperfine line 
c transitions are only used in the case of non-LTE 
c hyperfine structure and only in the stat. eq. calculations.
c We never calculate the radiative transfer for the hyperfine lines. 
c Instead we derive the brightness for the Jbar of the individual
c hyperfines by dividing up the Jbar of the total rotational 
c transition.

c INDEXU, INDEXL refer to the upper and lower levels of the 
c rotational lines. These are always rotational levels and never 
c hyperfine levels. These indices count in F77 style 1,N.

c JINDX, KINDX, LINDX refer to the quantum numbers of the 
c rotational levels. They are defined differently for different
c molecules. One way to get a quantum numbers of the upper
c level of a line is, JINDX(INDEXU(LINE)).

c UPPER, LOWER are similar to indexu and indexl except that
c UPPER and LOWER refer to the two levels of an allowed radiative
c transition between hyperfine levels. These indices count
c in C/IDL style from 0,N-1 even in the F77 code. In the F77
c code you need to add 1 to use them as array indices.


c23456789112345678921234567893123456789412345678951234567896123456789712
        integer function f77open(idno)

        character *8 fnotes(2)
        character *10 spectrum(100)
        logical there

        common /procid/ myid

        data fnotes / 'Fnotes0','Fnotes1' /

c Copy the process ID into a common variable
c Only process IDs 0 and 1 will be allowed to print to notes files.

        f77open = 1

        myid = idno
c        write(6,*) 'f77open: proc id number = ',myid
        if (myid .gt. 1) then 
            return
        endif
        i = myid 

        inquire(file=fnotes(i+1),exist=there)
        if (there) then
            open (unit=22,file=fnotes(i+1),STATUS='OLD')
c            open (unit=22,name=fnotes(i+1),STATUS='OLD')
            close(unit=22,STATUS='DELETE')
            inquire(file=fnotes(i+1),exist=there)
            if (there) then
                write(6,*) 'old file ',fnotes(i+1),
     &               '  could not be deleted'
                f77open = 0
            else 
                write(6,*) 'old file ',fnotes(i+1),'  was deleted'
                f77open = 1
            endif
        endif

c        open (unit=22,name=fnotes(i+1),form='formatted',
        open (unit=22,file=fnotes(i+1),form='formatted',
     &        STATUS='UNKNOWN')

        write(6,*) 'File opened  ',fnotes(i+1),' as unit 22'
        write(22,*)'File opened  ',fnotes(i+1),' as unit 22'

        inquire(file=fnotes(i+1),exist=there)
        if (there) then f77open = 1
        if (.not. there) then 
           write(6,*) 'File was not found   ',fnotes(i+1),' as unit 22'
           write(22,*)'File was not found   ',fnotes(i+1),' as unit 22'
           f77open = 0
        endif

c        write(6,*) 'f77open exit status is ',f77open
        write(22,*) 'f77open exit status is ',f77open

        return
        end

        integer function channels(achanwd,velrange,avwmin,avwmax,
     &                    nchanup,chvelup)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        logical doprint

        include 'nlines_f77.h'
        parameter(maxch=2000)
        parameter(maxhyp=50)
        parameter(maxnw=500)

        common /procid/ myid
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /vgrid/ chvel(maxch,nlines),chfreq(maxch,nlines),
     &      nchan(nlines),profile0(maxch,maxnw),chanwd,
     &      vwmin,vwmax
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        dimension startv(maxhyp),endv(maxhyp),nchanup(nlines)
        dimension hvstart(maxhyp),hvend(maxhyp)
        dimension chvelup(maxch,nlines)
c        dimension sumprofile0(maxch,maxnw,nlines)


c ------------------------------------------------------------

        doprint = .false.
        doprint = .true. 
        if (myid .ge. 2) doprint = .false.

c transfer from arguments to common
        chanwd = achanwd
        vwmin = avwmin
        vwmax = avwmax
c this is the successful return value, will be set to 0 if there
c is an error
        channels = 1

c          write(6,*) 'starting function channels'
c          write(6,*) 'chanwd ',chanwd/1.e5,' velrange ',velrange/1.e5
c          write(6,*) 'vwmin' ,vwmin/1.e5,' vwmax ',vwmax/1.e5

        if (doprint) then
          write(22,*) 'starting function channels'
          write(22,*) 'chanwd ',chanwd/1.e5,' velrange ',velrange/1.e5
          write(22,*) 'vwmin' ,vwmin/1.e5,' vwmax ',vwmax/1.e5
        endif


c Compute the channel velocities

        do 600 line = 1,nlines

        if (doprint) 
     &      write(22,*) 'Start and end velocities for new line',line

c For each hyperfine belonging to LINE, the STARTV and ENDV are the
c minus and plus the VELRANGE set in SETUP.C
            do 601 ihyp = 1,nhyp(line)
                hvstart(ihyp)   = hypvel(ihyp,line) - velrange/2.
                hvend(ihyp)     = hypvel(ihyp,line) + velrange/2.

        if (doprint) 
     &      write(22,1077) ihyp,hvstart(ihyp)/1.e5,hvend(ihyp)/1.e5
 1077   format('ihyp startv endv ',i5,1p,2e12.4)

 601        continue

c Start with the first 2 hyperfine lines 1 and 2. Because the
c hyperfines are in order of increasing frequency, the lowest
c velocity hyperfine is the last one (nhyp(line) in the list.
            ihyp1 = nhyp(line)
            ihyp2 = ihyp1 - 1

c The first velocity band is initially set to cover the velocity
c range around the first hyperfine.
            nbands = 1
            startv(nbands) = hvstart(ihyp1)
            endv(nbands)   = hvend(ihyp1)

c If there is only one hyperfine line then we are all done
                if (ihyp2 .le. 0) go to 603

c This is a counter in case we get stuck in an infinite loop
            isafe = 1

 602        continue

c Just print out.
        if (doprint) then
        write(22,*) 'band=',nbands,' start end vel ',
     &    startv(nbands)/1.e5,endv(nbands)/1.e5
        write(22,*) 'h1 h2 = ',ihyp1,ihyp2,' and vel ',
     &    hvstart(ihyp1)/1.e5,hvend(ihyp1)/1.e5,
     &    hvstart(ihyp2)/1.e5,hvend(ihyp2)/1.e5
        endif

c If the endv of the first band is > than the startv
c of the next hyperfine, then expand the band to cover this
c hyperfine. Set the endv of the 1st band to be the endv of the 
c 2nd hyperfine. Replace the second hyperfine by the 3rd hyperfine. 
c If the combined velocity range of the 1st and 2nd overlaps with 
c the start of the velocity range of the 3rd line, expand the 
c 1st band again, and keep going (back to 602).
            if (endv(nbands) .ge. hvstart(ihyp2)) then
                endv(nbands) = hvend(ihyp2)
                ihyp2 = ihyp2 - 1

c If there is no next line then we are all done
                if (ihyp2 .le. 0) go to 603
                go to 602
            endif

c We come to this point if the velocity range of the first
c several (or initial single) hyperfine does not overlap with the 
c the next. This means that the velrange we have so far is one
c band and there is a gap until we come to the velrange around
c the next hyperfine. 

c The startv and endv of the first band are the velrange
c identified in the 602 loop above. Done with the first band.

c Move on to the next pair of lines to repeat the process for this
c second band. 
                ihyp1 = ihyp2
                ihyp2 = ihyp1 - 1

c Set the initial velocity range of the next band to be the 
c velrange of the next hyperfine (ihyp1).
                nbands = nbands + 1
                startv(nbands) = hvstart(ihyp1)
                endv(nbands)   = hvend(ihyp1)

c If it turns out that this ihyp1 is the last line (counting backwards
c in the list) then we just set up the band for this last line,
c and we can quit. Use .le. instead of .eq. just to make sure
c we dont get in an absurd loop running down a negative list.
                if (ihyp2 .le. 0) go to 603

c If there is more than one line left, then we are going to go 
c back to the start, 602.

c First make sure we are not in an infinite loop. 
c This cant happen, but just in case.
            isafe = isafe + 1
            if (isafe .gt. 200) then
                if (myid .lt. 2) 
     &            write(22,*) 'stopping because of infinite loop in',
     &            ' setting up velocity grid'
                  write(6,*) 'stopping because of infinite loop in',
     &            ' setting up velocity grid'
                channels = 0
                return
            endif

c We are not at the end of the hyperfine list, we are not in an
c infinite loop.  Go back to the start to work on this next band.
            go to 602

c All done with all the hyperfines.
 603        continue

        if (doprint) write(22,*) 'Nbands ',nbands

        do 604 n = 1,nbands
                ichan = (endv(n)-startv(n))/chanwd + 1
                nchan(line) = nchan(line) + ichan
                if (myid .lt. 2 .and. doprint) 
     &              write(22,*) 'new start and end velocities, nchan',
     &              startv(n)/1.e5,endv(n)/1.e5,
     &              nchan(line)

                if (nchan(line) .gt. maxch) then
                  if (myid .lt.2) 
     &              write( 6,*) 'stopping because the number of',
     &              ' channels required exceeds the maximum of '
                 if (myid .lt.2) 
     &              write( 6,*)'line=',line,'  channels=',nchan(line)
     &              ,maxch
                 if (myid .lt. 2) write(22,*) 
     &              'stopping because the number of',
     &              ' channels required exceeds the maximum of '
     &              ,maxch
                 if (myid .lt. 2) 
     &              write(22,*)'line=',line,'  channels=',nchan(line)
                 channels = 0
                 return
                endif

                do 605 i = 1,ichan
                        ii = nchan(line)-ichan + i
                        chvel(ii,line) = startv(n) + (i-1)*chanwd
                        chvelup(ii,line) = chvel(ii,line) 
 605            continue

 604    continue

 600    continue

        nchanmax = 0
        do 962 line = 1,nlines
                nchanmax = max(nchan(line),nchanmax)
                nchanup(line) = nchan(line)
 962    continue


c For all the lines, compute the channel frequencies from
c the velocities.
        do 961 line = 1,nlines
        do 961 l = 1,nchan(line)
             chfreq(l,line) = freqline(line)
     &                  * (1.d0 - dble(chvel(l,line))/c)
 961    continue

c This section prints out some information about the velocity grid
        if (doprint) then
         do 260 line = 1,nlines
             write(22,1028) freqline(line)/1.e9
             write(22,*) 'The number of channels is ', nchan(line)
             write(22,*) 'ch number, ch velocity, ',
     &        'ch frequency'
             do 63 l = 1,nchan(line)
 63              write(22,1000) 
     &              l,chvel(l,line)/1.e5,chfreq(l,line)/1.e9
 260        continue
        endif
 1000   format(i5,1p,2e18.10)
 1028   format(//,'Here is the velocity grid',
     &    ' for transition frequency ',f16.8)

c Finished with the velocity grid. 
c Set up the grid of Gaussian profiles at different linewidths

        if (myid .lt. 2) 
     &    write(22,*) 'Minimum and maximum linewidths ',vwmin,vwmax

        nwidths = (vwmax - vwmin)/chanwd + 2

        if (nwidths .gt. maxnw ) then
         if (myid .lt. 2) then
          write(22,*) 'dimension MAXNW is too small'
          write(6,*) 'dimension of MAXNW is too small'
          write(22,*) 'The number of template profiles is '
          write(22,*) '(max width - min width)/chan width'
          write(22,*) 'This is ',nwidths, ' and exceeds the dimension '
          write(22,*) 'of the template, MAXNW,  is ',maxnw
         endif
                channels = 0
                return
        endif
        do 166 line = 1,nlines
        if (nchan(line) .gt. maxch) then
           if (myid .lt. 2) then
                write(22,*) 'The number of channels for line ',line
                write(22,*) 'exceeds the dimension, MAXCH ',maxch
                write(6,*) 'The number of channels for line ',line
                write(6,*) 'exceeds the dimension, MAXCH ',maxch
            endif
            channels = 0
            return
        endif
 166    continue


c The profile function is in units of 1/frequency and normalized
c so that the integral over frequency is equal to unity.
c The normalization will be applied in SHIFTPROFILE. Here
c we build the profiles, each with a height of 1.

c tau = (line center opacity)*(path length)*(profile)

c        write(6,*) 'Computing profiles'
        if (myid .lt. 2) write(22,*) 'Computing ',nwidths,' profiles',
     &     ' starting from width of ',vwmin,' step = ',chanwd

c The widths in the grid of profiles increase by 1 channel width
c The channel maxch/2 is the zero velocity channel for both
c cases maxch = even or odd.
        do 494 nw = 1,nwidths
            vwidth = vwmin + (nw-1)*chanwd
c        if (doprint) write(22,*) 
c     &      'New profile, width =',vwidth,maxch
            do 494 l = 1,maxch
                profile0(l,nw) = 
     &            exp( -( (l - maxch/2)*chanwd/vwidth )**2 )
c        if (doprint) write(22,*),l,profile0(l,nw)
 494    continue
 
c This is for debugging.
c        write(6,*) 'Process ',myid,' i
c     &     Sending stop signal at end of channels '
c        channels = 0

        return
        end 

        function setup(m,nl,nst,iid)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline,
     &  energy
        integer molecule,m,u,l
        logical doprint

        include 'nlines_f77.h'

        parameter(maxhyp=50)

        common /procid/ myid
        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /lines/ indexu(nlines),indexl(nlines)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /coll/ cr(nstate,nstate)
        common /cenergy/ energy(nstate)

c This function loads the f77 common blocks with information about the 
c selected molecule. All of the information in the common blocks above
c is filled in except the collision rates, CR. These are here to calculate
c the critical density that is printed in Cnotes0 for information. The
c critical density is not used in the calculations.

        doprint = .true.
        doprint = .false.
        molecule = m
        setup = 1.

        nltehyp = 0
        if (molecule .eq. 4 ) nltehyp = 1
        if (molecule .eq. 21) nltehyp = 1
        if (molecule .eq. 18) nltehyp = 1

c This is a poor if statement >= 18 and < 19 means 18 only
c        if (molecule .eq. 4 
c     &      .or. (molecule .ge. 18 .and. molecule .lt. 19) ) 
c     &        nltehyp = 1

c This function is called by cmain. The process ID is set
c in cmain and passed through the arguments. Does not use 
c the MYID from common /proc/

        if (iid .eq. 0) then
        if (nl .ne. nlines .or. nst .ne. nstate) then
            write(6,*) 'ERROR: The number of lines and levels ',
     &       'in the f77 and C nline headers must be set equal'
            write(6,*) 'in nlines_f77.h nlines = ',nlines,
     &       ' nstate = ',nstate
            write(6,*) 'in nlines_C.h   nlines = ',nl,
     &       ' nstate = ',nst
            setup = 0.
            return
        endif

        i = nstate
        j = nlines
c This could be changed to test more rotors than just CO
        if (m .eq. 7 .and. i .le. j ) then
           write(6,*) 'PROBLEM: For molecules that are rotors ',
     &      'the number of states should be more than the number ',
     &      'of lines'
           write(6,*) 'nlines = ',nlines,' but nstate = ',nstate
           setup = 0.
           return
        endif
        endif

c This calls a function to define physical constants. Call
c this for every molecule.
        call setcon

        istatus = 1

c The rotors are currently molecules 7 through 16
        if (m .ge. 7 .and. m .le. 16) call rotor(istatus)
        if (m .eq. 20) call rotor(istatus)

c 21 is HCN with non-LTE hyperfines
        if (m .eq. 21) call hcn_hyp_init(istatus)

c 6 is H2O
        if (m .eq. 6) call h2o_init(istatus)

c 22 is H2D+
c        write(6,*) 'starting init for H2D+'
        if (m .eq. 22) call h2d_init(istatus)

c 3 is NH3 in LTE
        if (m .eq. 3) call nh3int

c 4 is N2H+ with non-LTE hyperfines
        if (m .eq. 4) call n2h_hyp_init(istatus)

c        do 79 i = 1,nstate
c            write(22,*) 'energies from common ',i,energy(i)
c 79     continue

c 2 is CH3CN in LTE
        if (m .eq. 2) call ch3cnint

c 17,18,19 are HCN vibrational
        if (m .ge. 17) call moldata(istatus)

        if (istatus .eq. 0) then
           if (iid .lt. 2) write(6,*) 
     &       'molecule setup returned with bad status = ',istatus
           if (iid .lt. 2) write(22,*) 
     &       'molecule setup returned with bad status = ',istatus
           setup = 0.
           return
        endif


        setup = atoms
c        write(6,*),'P',iid,':  next load collision rates '



c The only process that needs collision rates is process ZERO,
c the master because it does the stat eq calculations
c But it doesn't much matter


c Load the collision rates for HCO+=12, H13CO+=11,N2H+=9,N2D+=16
        if (m .eq.  9) call loadionrates
        if (m .eq. 11) call loadionrates
        if (m .eq. 12) call loadionrates
        if (m .eq. 16) call loadionrates

c These are the collision rates for H2O
c Just for debugging, use these rates for H2D,
c but for H2D, the H2O rates are not even in the right order.
        if (m .eq. 6) call loadh2orates
        if (m .eq. 22) call loadh2orates

c These are the rates for the rotational levels of
c the ground vibrational state of HCN
        if (m .eq. 15) call loadhcnrates

c Use the delta J rates here and divide into
c hyperfine rates in colrate_hcnhyp2
        if (m .eq. 21) call loadhcnrates


c Load the collision rates for N2H+ with non-LTE hyperfines
c If you edit colrateion, you can run N2H+ using these ion rates
c instead of the LAMDBA rates. Instructions in colrateion.
        if (m .eq. 4) call n2h_deltaJ_rates

c -------------------------------------------------------------
c Calculate the critical density for de-excitation. To do this
c we need the collision rates. Set the temperature to 10 K.

        tk = 10.
        istatus = 1
c Set the collision rates for this temperature
        if (molecule .eq. 6) then
                call colrateh2o(tk,doprint)
        endif

        if (molecule .eq. 22) then
                call colrateh2o(tk,doprint)
        endif

c Collision rates for molecules CO=7, C13O=10, C17O=13, C18O=14
        if (     molecule .eq. 7 
     &      .or. molecule .eq. 10 
     &      .or. molecule .eq. 13
     &      .or. molecule .eq. 14 
     &     ) then
                call colratco(tk)
        endif

c Collision rates for rotational levels of HCN, ground vibrational
        if (molecule .eq. 15) then
                call colrate_hcn_ground(tk)
        endif

c Collision rates for CS = 8 and SiO = 20
        if (molecule .eq. 8 .or. molecule .eq. 20) then
                call colratcs(tk)
        endif

c Collision rates for N2Hplus = 9 and N2Dplus = 16
        if (molecule .eq. 9 .or. molecule .eq. 16) then
                call colrateion(tk)
        endif

c Collision rates for H13COplus = 11 and HCOplus = 12
        if (molecule .eq. 11 .or. molecule .eq. 12) then
c                call colratehco(tk,istatus)
                call colrateion(tk)
        endif

c The subroutine loadionrates should have been run in F77 function
c setup to use colrateion. This is done automatically, and this
c comment is a reminder that colrateion depends on loadionrates

c Collision rates for N2Hhyp = 4 (N2H+ with non-LTE hyperfines)
c     colrate_nhyp uses the rotational rates from the LAMBDA data base
c     and approximates the rates between hyperfine levels using the
c     statistical weights. The subroutine n2h_deltaJ_rates should have
c     been run in F77 function setup to use colrate_nhyp

c If you edit colrateion, you can use the ion rates for N2H+ with 
c non-LTE hyperfines. Also edit F77 setup to run loadionrates

        if (molecule .eq. 4 ) then
                iprint = 0
c                if (igrid .lt. 4 .and. jgrid .eq. 0 
c     &              .and. kgrid .eq. 0) iprint = 1
                call colrate_nhyp(tk,iprint)
        endif



c New colrathcn (HCN has used colratco before)
        if (     molecule .eq. 19
     &      .or. molecule .eq. 17
     &      .or. molecule .eq. 18) call colrathcn(tk)



c Use colrate_hcnhyp2 with loadhcnrates
       if (molecule .eq. 21 ) then
                iprint = 0
c                if (igrid .lt. 4 .and. jgrid .eq. 0 
c     &              .and. kgrid .eq. 0) iprint = 1
                call colrate_hcnhyp2(tk,iprint)
        endif

c The critical density

        if (iid .lt. 2) write(22,*) 
     &     'Critical density: line upper lower  A CR density '
        do 1 i = 1,nlines
             u = indexu(i)
             l = indexl(i)
             critical = a(u,l)/cr(u,l)
             if (iid .lt. 2)write(22,1000) i,u,l,a(u,l),cr(u,l),critical
 1      continue
 1000   format(3i5,1p,3e12.3)

c -------------------------------------------------------------

        if (iid .lt. 2) write(22,*) 'Setup finished with atoms', atoms
c        write(6,*),'P',iid,':  function setup: atoms = ',atoms,setup

        return
        end

        function crtsph(x,y,z,r,az,el,idir)

        real x,y,z,r,az,el,crtsph

c Converts cartesian coords to spherical or back.

        pi =  3.14159265359d0

c eliminate the warning about crtsph used before defined
        crtsph = pi

c IDIR > 1 cart to sph.
        if (idir .le. 0) go to 1
c First find the radius:
        r = sqrt(x*x + y*y + z*z)
c Then the elevation angle:
        if ( r.eq.0.0) then
                el = 0.0
                az = 0.0
                return
        else
                zsr = z/r
                if (abs(zsr) .gt. 1.0) zsr = sign(1.0,zsr)
                el = acos(zsr)
        endif
c Now the azimuth angle:
        sinel = sin(el)
        if (sinel .eq. 0.0) then
                az = 0.0
                return
        endif
        sinaz = y/(r*sinel)
        cosaz = x/(r*sinel)
        if (abs(cosaz) .gt. 1.0) cosaz = sign(1.0,cosaz)
        az = acos(cosaz)
c Those are the principle branches. Now find the correct quadrant.
c Both sin and cos pos ==> first quad.
c sin=+, cos=- ==> 2nd quad
c both - ==> 3rd quad
c sin=-, cos=+ ==> 4th quad.
        if (sinaz .ge. 0.0) then
                        return
        else
                        az = 2.*pi - az
                        if (az .eq. 2*pi) az = 0.
                        return
        endif
 1      continue
c
c Convert sph to cart
c
        sinel = sin(el)
        x = r*sinel*cos(az)
        y = r*sinel*sin(az)
        z = r*cos(el)
c
        return
        end

        subroutine setcon

c Sets universal constants.

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,aline,freqline,guline,glline
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav


c This is all straightforward.

                planck = 6.6252d-27
                boltz = 1.38d-16
                c = 3.d10
                amu = 1.660531d-24
                pc = 3.085d18
                grav = 6.67d-8
                pi = 3.14159265359
                sqrtpi = 1.77245385091


        return
        end

c23456789112345678921234567893123456789412345678951234567896123456789712
        subroutine h2d_init(istatus)

        include 'nlines_f77.h'
        parameter (maxch=2000)
        parameter (maxhyp=50)
        parameter (maxtrans=500)
 
        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,gco,aline,freqline,guline,glline
 
        character *80 msg
        integer upper,lower,istatus 
        logical docont,doprint,doprint2,nh3tau
 
        common /procid/ myid
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /lines/ indexu(nlines),indexl(nlines)
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /nltehyp/ upper(maxtrans),lower(maxtrans),
     &          istartline(nlines)

c This subroutine loads up the molecular line data into freq and a
c that are indexed by state. Depends on copyatomic to transfer this
c data to the variables freqline and aline that are indexed by line.

c Energy levels are in inverse cm from NRAO splatalogue
        dimension energy(nstate)

        write(6,*) 'Start defining',nstate,'  states and',nlines,
     &              '  lines for H2D+'

        atoms = 4.d0
        istatus = 1

        nst = nstate
        nli = nlines
        write(6,*) 'Start defining',nstate,'  states and',nlines,
     &              '  lines for H2D+'
        if (nst .gt. 18 .or. nli .gt. 31) then
          istatus = 0
          write(6,*)'Only have H2D+ data for 18 levels, 31 lines'
          write(6,*) 'nstate = ',nstate,'  nlines = ',nlines
          write(6,*) 'Program will not work, must exit ...'
          if (myid .lt. 2) then 
            write(22,*)'Only have H2D+ data for 18 levels, 31 lines'
            write(22,*)'nstate = ',nstate,'  nlines = ',nlines
            write(22,*)'Program will not work, must exit ...'
          endif
          return
        endif

c Information for each level

        i = 1
        jindx(i) = 1
        kindx(i) = 1
        lindx(i) = 1
        statdg(i) = 9
        energy(i) = 60.0310

        i = 2
        jindx(i) = 1
        kindx(i) = 1
        lindx(i) = 0
        statdg(i) = 9
        energy(i) = 72.4537

        if (nst .ge. 3) then
          i = 3
          jindx(i) = 2
          kindx(i) = 1
          lindx(i) = 2
          statdg(i) = 15
          energy(i) = 138.8631
        endif

        if (nst .ge. 4) then
          i = 4
          jindx(i) = 2
          kindx(i) = 1
          lindx(i) = 1
          statdg(i) = 15
          energy(i) = 175.9389
        endif

        if (nst .ge. 5) then
          i = 5
          jindx(i) = 3
          kindx(i) = 1
          lindx(i) = 3
          statdg(i) = 27
          energy(i) = 254.06630
        endif

        if (nst .ge. 6) then
          i = 6
          jindx(i) = 3
          kindx(i) = 1
          lindx(i) = 2
          statdg(i) = 21
          energy(i) = 326.16520
        endif

        if (nst .ge. 7) then
          i = 7
          jindx(i) = 4
          kindx(i) = 1
          lindx(i) = 4
          statdg(i) = 27
          energy(i) = 403.68130
        endif

        if (nst .ge. 8) then
          i = 8
          jindx(i) = 3
          kindx(i) = 3
          lindx(i) = 1
          statdg(i) = 21
          energy(i) = 458.34840
        endif

        if (nst .ge. 9) then
          i = 9
          jindx(i) = 3
          kindx(i) = 3
          lindx(i) = 0
          statdg(i) = 21
          energy(i) = 459.83370
        endif

        if (nst .ge. 10) then
          i = 10
          jindx(i) = 4
          kindx(i) = 1
          lindx(i) = 3
          statdg(i) = 27
          energy(i) = 516.15200
        endif

        if (nst .ge. 11) then
          i = 11
          jindx(i) = 5
          kindx(i) = 1
          lindx(i) = 5
          statdg(i) = 33
          energy(i) = 586.47600
        endif

        if (nst .ge. 12) then
          i = 12
          jindx(i) = 4
          kindx(i) = 3
          lindx(i) = 2
          statdg(i) = 27
          energy(i) = 645.57140
        endif

        if (nst .ge. 13) then
          i = 13
          jindx(i) = 4
          kindx(i) = 3
          lindx(i) = 1
          statdg(i) = 27
          energy(i) = 654.50710
        endif

        if (nst .ge. 14) then
          i = 14
          jindx(i) = 5
          kindx(i) = 1
          lindx(i) = 4
          statdg(i) = 33
          energy(i) = 738.74580
        endif

        if (nst .ge. 15) then
          i = 15
          jindx(i) = 6
          kindx(i) = 1
          lindx(i) = 6
          statdg(i) = 39
          energy(i) = 801.75063
        endif

        if (nst .ge. 16) then
          i = 16
          jindx(i) = 5
          kindx(i) = 3
          lindx(i) = 3
          statdg(i) = 33
          energy(i) = 876.54889
        endif

        if (nst .ge. 17) then
          i = 17
          jindx(i) = 5
          kindx(i) = 3
          lindx(i) = 2
          statdg(i) = 33
          energy(i) = 904.07052
        endif

        if (nst .ge. 18) then
          i = 18
          jindx(i) = 6
          kindx(i) = 1
          lindx(i) = 5
          statdg(i) = 39
          energy(i) = 991.05674
        endif

        write(6,*) 'Finished defining H2D+ levels'

c Calculate the frequencies of all transitions using the energy levels.
c Convert the energies in cm-1 to GHz by multiplying by C

        do 2 i = 1,nstate
        do 2 j = 1,nstate
          freq(i,j) = c*abs(energy(i) - energy(j))
 2      continue

c Information for each line. Transition frequencies are from the
c JPL catalog. These should be the most accurate and will overwrite
c the frequencies calculated above.

c JINDX(LINE) is the quantum number for line number LINE

c INDEXU(LINE) and INDEXL(LINE) are the Fortran array index numbers for 
c upper and lower states of the line.These are quantum number + 1

        do 1 i = 1,nstate
        do 1 j = 1,nstate
          a(i,j) = 0.
 1        continue

c 110 -> 111
        line = 1
        indexu(line) = 2 
        indexl(line) = 1 
        a(indexu(line),indexl(line)) = 1.082222e-4
        a(indexl(line),indexu(line)) = 
     &        a(indexu(line),indexl(line)) 
        freq(indexu(line),indexl(line)) = 372.421385e9 
        freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
c        write(6,*)'h2d_init: lin,indx_u,indx_l,freqd,frequ ',
c     &        line, indexu(line), indexl(line),
c     &         freq(indexu(line),indexl(line)),
c     &          freq(indexl(line),indexu(line))


c 212 -> 111
        if (nli .ge. 2) then
          line = 2
          indexu(line) = 3
          indexl(line) = 1
          a(indexu(line),indexl(line)) = 1.659336e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 2363.325058e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 211 -> 110
        if (nli .ge. 3) then
          line = 3
          indexu(line) = 4
          indexl(line) = 2
          a(indexu(line),indexl(line)) = 3.753891e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 3102.457981e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 211 -> 212
        if (nli .ge. 4) then
          line = 4
          indexu(line) = 4
          indexl(line) = 3
          a(indexu(line),indexl(line)) = 9.589778e-4
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 1111.505229e9 
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 313 -> 212
        if (nli .ge. 5) then
          line = 5
          indexu(line) = 5
          indexl(line) = 3
          a(indexu(line),indexl(line)) = 6.482111e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 3453.796250e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 331 -> 212
        if (nli .ge. 6) then
          line = 6
          indexu(line) = 8
          indexl(line) = 3
          a(indexu(line),indexl(line)) = 2.002191e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 9577.929144e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 312 -> 211
        if (nli .ge. 7) then
          line = 7
          indexu(line) = 6
          indexl(line) = 4
          a(indexu(line),indexl(line)) = 1.405535e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 4503.672623e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 330 -> 211
        if (nli .ge. 8) then
          line = 8
          indexu(line) = 9
          indexl(line) = 4 
          a(indexu(line),indexl(line)) = 3.554286e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 8510.952209e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 312 -> 313
        if (nli .ge. 9) then
          line = 9
          indexu(line) = 6
          indexl(line) = 5
          a(indexu(line),indexl(line)) = 3.835454e-3
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 2161.471602e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 330 -> 313
        if (nli .ge. 10) then
          line = 10
          indexu(line) = 9
          indexl(line) = 5
          a(indexu(line),indexl(line)) = 2.188888e-3
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 6168.751188e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 414 -> 313
        if (nli .ge. 11) then
          line = 11
          indexu(line) = 7 
          indexl(line) = 5 
          a(indexu(line),indexl(line)) = 1.538892e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 4485.343733e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 331 -> 312
        if (nli .ge. 12) then
          line = 12
          indexu(line) = 8
          indexl(line) = 6
          a(indexu(line),indexl(line)) = 4.373793e-3
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 3962.751292e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 413 -> 312
        if (nli .ge. 13) then
          line = 13
          indexu(line) = 10
          indexl(line) = 6
          a(indexu(line),indexl(line)) = 2.943074e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 5695.661172e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 431 -> 312
        if (nli .ge. 14) then
          line = 14
          indexu(line) = 13
          indexl(line) = 6
          a(indexu(line),indexl(line)) = 1.139447e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 9843.441697e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 413 -> 414
        if (nli .ge. 15) then
          line = 15
          indexu(line) = 10
          indexl(line) = 7
          a(indexu(line),indexl(line)) = 1.047663e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 3371.789041e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 431 -> 414
        if (nli .ge. 16) then
          line = 16
          indexu(line) = 13
          indexl(line) = 7
          a(indexu(line),indexl(line)) = 3.217750e-3
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 7519.569565e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 515 -> 414
        if (nli .ge. 17) then
          line = 17
          indexu(line) = 11
          indexl(line) = 7
          a(indexu(line),indexl(line)) = 2.927744e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 5480.048905e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 330 -> 331
        if (nli .ge. 18) then
          line = 18
          indexu(line) = 9
          indexl(line) = 8
          a(indexu(line),indexl(line)) = 2.67738e-7
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 44.528294e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 432 -> 331
        if (nli .ge. 19) then
          line = 19
          indexu(line) = 12
          indexl(line) = 8
          a(indexu(line),indexl(line)) = 1.462558e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 5612.806087e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 413 -> 330
        if (nli .ge. 20) then
          line = 20
          indexu(line) = 10
          indexl(line) = 9
          a(indexu(line),indexl(line)) = 6.237181e-5
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 1688.381586e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 431 -> 330
        if (nli .ge. 21) then
          line = 21
          indexu(line) = 13
          indexl(line) = 9
          a(indexu(line),indexl(line)) = 1.660701e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 5836.162110e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 432 -> 413
        if (nli .ge. 22) then
          line = 22
          indexu(line) = 12
          indexl(line) = 10
          a(indexu(line),indexl(line)) = 1.100810e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 3879.896207e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 514 -> 413
        if (nli .ge. 23) then
          line = 23
          indexu(line) = 14
          indexl(line) = 10
          a(indexu(line),indexl(line)) = 4.798519e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 6673.194028e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 514 -> 515
        if (nli .ge. 24) then
          line = 24
          indexu(line) = 14
          indexl(line) = 11
          a(indexu(line),indexl(line)) = 2.136213e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 4564.934164e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 616 -> 515
        if (nli .ge. 25) then
          line = 25
          indexu(line) = 15
          indexl(line) = 11
          a(indexu(line),indexl(line)) = 4.913472e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 6453.771060e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 431 -> 432
        if (nli .ge. 26) then
          line = 26
          indexu(line) = 13
          indexl(line) = 12
          a(indexu(line),indexl(line)) = 3.125427e-5
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 267.884317e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 533 -> 432
        if (nli .ge. 27) then
          line = 27
          indexu(line) = 16
          indexl(line) = 12
          a(indexu(line),indexl(line)) = 4.055505e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 6924.531062e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 514 -> 431
        if (nli .ge. 28) then
          line = 28
          indexu(line) = 14
          indexl(line) = 13
          a(indexu(line),indexl(line)) = 3.036282e-4
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 2525.413503e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 532 -> 431
        if (nli .ge. 29) then
          line = 29
          indexu(line) = 17
          indexl(line) = 13
          a(indexu(line),indexl(line)) = 5.311592e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 7481.723256e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 533 -> 514
        if (nli .ge. 30) then
          line = 30
          indexu(line) = 16
          indexl(line) = 14
          a(indexu(line),indexl(line)) = 1.917319e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 4131.233241e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 615 -> 514
        if (nli .ge. 31) then
          line = 31
          indexu(line) = 18
          indexl(line) = 14
          a(indexu(line),indexl(line)) = 7.183288e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 7564.091578e9
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

        line = 1
c        write(6,*)'h2o_init: lin,indx_u,indx_l,freqd,frequ ',
c     &        line, indexu(line), indexl(line),
c     &         freq(indexu(line),indexl(line)),
c     &          freq(indexl(line),indexu(line))

c Here set up the molecule's hyperfine structure for no hyperfine
c splitting. In this case the line will have one "hyperfine" whose relative
c intensity is 1 and whose velocity is not shifted

c These values will be reset by hyperfine subroutine if the
c molecule has hyperfine structure.

        do 411 line = 1,nlines
          nhyp(line) = 1
 411    continue

        do 408 line = 1,nlines
          do 409 ihyp = 2,maxhyp
            relint(ihyp,line) = 0.
            hypvel(ihyp,line) = 0.
 409      continue
          relint(1,line) = 1.
          hypvel(1,line) = 0.
 408    continue


        do 440 line = 1,nlines
           istartline(line) = line - 1
           n = indexu(line)
           m = indexl(line)
           aline(line) = a(n,m)
           freqline(line) = freq(n,m)
           glline(line) = statdg(m)
           guline(line) = statdg(n)
           lower(line)  = line - 1
           upper(line)  = line 
 440    continue

        lun = 22
c        msg = 'matrix of Einstein A coefficients'
c        call dmattmt(a,msg,nstate,lun)
c        msg = 'matrix of frequencies'
c        call dmattmt(freq,msg,nstate,lun)

        end

c23456789112345678921234567893123456789412345678951234567896123456789712
        subroutine h2o_init(istatus)

        include 'nlines_f77.h'
        parameter (maxch=2000)
        parameter (maxhyp=50)
        parameter (maxtrans=500)
 
        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,gco,aline,freqline,guline,glline
 
        character *80 msg
        integer upper,lower 
        logical docont,doprint,doprint2,nh3tau
 
        common /procid/ myid
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /lines/ indexu(nlines),indexl(nlines)
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /nltehyp/ upper(maxtrans),lower(maxtrans),
     &          istartline(nlines)
        dimension energy(nstate)


c This subroutine loads up the molecular line data into freq and a
c that are indexed by state. Depends on copyatomic to transfer this
c data to the variables freqline and aline that are indexed by line.
c Energy levels are in inverse cm from Green 1980 ApJS 42 103

        istatus = 1

        atoms = 18.d0

c There are more 10 states and 17 lines in this subroutine,
c but there are only collision rates for 8 states. Therefore
c cannot run H2O with more than 8 states (11 lines).

        nst = nstate
        nli = nlines
        if (nst .gt. 8 .or. nli .gt. 11) then
          write(6,*)'Only have coll. rates for 8 levels, 11 lines'
          write(6,*) 'nstate = ',nstate,'  nlines = ',nlines
          write(6,*) 'Program will not work, must exit ...'
          if (myid .lt. 2) then 
            write(22,*)'Only have coll. rates for 8 levels, 11 lines'
            write(22,*)'nstate = ',nstate,'  nlines = ',nlines
            write(22,*)'Program will not work, must exit ...'
          endif
          istatus = 0
          return
        endif

c Information for each level

c Levels are in order of energy same as LAMDB
c	1	101	
c	2	110	
c	3	212
c	4	221
c	5	303	limit of low T in LAMDB
c	6	312
c	7	321
c	8	414
c	9	330
c	10	423

        i = 1
        jindx(i) = 1
        kindx(i) = 0
        lindx(i) = 1
        statdg(i) = 9
        energy(i) = 23.79

        i = 2
        jindx(i) = 1
        kindx(i) = 1
        lindx(i) = 0
        statdg(i) = 9
        energy(i) = 42.36

        if (nst .ge. 3) then
          i = 3
          jindx(i) = 2
          kindx(i) = 1
          lindx(i) = 2
          statdg(i) = 15
          energy(i) = 79.50
        endif

        if (nst .ge. 4) then
          i = 4
          jindx(i) = 2
          kindx(i) = 2
          lindx(i) = 1
          statdg(i) = 15
          energy(i) = 134.90
        endif

        if (nst .ge. 5) then
          i = 5
          jindx(i) = 3
          kindx(i) = 0
          lindx(i) = 3
          statdg(i) = 21
          energy(i) = 136.76
        endif

        if (nst .ge. 6) then
          i = 6
          jindx(i) = 3
          kindx(i) = 1
          lindx(i) = 2
          statdg(i) = 21
          energy(i) = 173.37
        endif

        if (nst .ge. 7) then
          i = 7
          jindx(i) = 3
          kindx(i) = 2
          lindx(i) = 1
          statdg(i) = 21
          energy(i) = 212.16
        endif

        if (nst .ge. 8) then
          i = 8
          jindx(i) = 4
          kindx(i) = 1
          lindx(i) = 4
          statdg(i) = 27
          energy(i) = 224.84
        endif

        if (nst .ge. 9) then
          i = 9
          jindx(i) = 3
          kindx(i) = 3
          lindx(i) = 0
          statdg(i) = 21
          energy(i) = 285.42
        endif

        if (nst .ge. 10) then
          i = 10
          jindx(i) = 4
          kindx(i) = 2
          lindx(i) = 3
          statdg(i) = 27
          energy(i) = 300.37
        endif

c Calculate the frequencies of all transitions using the energy levels.

        do 2 i = 1,nstate
        do 2 j = 1,nstate
          freq(i,j) = c*abs(energy(i) - energy(j))
 2      continue

c Information for each line. Transition frequencies are from the
c JPL catalog. These should be the most accurate and will overwrite
c the frequencies calculated above.

c JINDX(LINE) is the quantum number for line number LINE
c INDEXU(LINE) is the Fortran array index number for the line
c that is quantum number + 1

        do 1 i = 1,nstate
        do 1 j = 1,nstate
          a(i,j) = 0.
 1        continue

c 110 -> 101	transition 1 states 2 1
        line = 1
        indexu(line) = 1 + 1
        indexl(line) = 0 + 1
        a(indexu(line),indexl(line)) = 3.4697e-3
        a(indexl(line),indexu(line)) = 
     &        a(indexu(line),indexl(line)) 
        freq(indexu(line),indexl(line)) = 556936.0020e6
        freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
c        write(6,*)'h2o_init: lin,indx_u,indx_l,freqd,frequ ',
c     &        line, indexu(line), indexl(line),
c     &         freq(indexu(line),indexl(line)),
c     &          freq(indexl(line),indexu(line))


c 212 -> 101	transition 2 states 3 1
        if (nli .ge. 2) then
          line = 2
          indexu(line) = 2 + 1
          indexl(line) = 0 + 1
          a(indexu(line),indexl(line)) = 5.6107e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 1669904.7750e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 221 -> 110	transition 3 states 4 2
        if (nli .ge. 3) then
          line = 3
          indexu(line) = 3 + 1
          indexl(line) = 1 + 1
          a(indexu(line),indexl(line)) = 2.5824e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 2773976.5880e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 221 -> 212	transition 4 states 4 3
        if (nli .ge. 4) then
          line = 4
          indexu(line) = 3 + 1
          indexl(line) = 2 + 1
          a(indexu(line),indexl(line)) = 3.0660e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 1661007.6370e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 303 -> 212	transition 5 states 5 3
        if (nli .ge. 5) then
          line = 5
          indexu(line) = 4 + 1
          indexl(line) = 2 + 1
          a(indexu(line),indexl(line)) = 5.0712e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 1716769.6330e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 312 -> 303	transition 7 states 6 5
        if (nli .ge. 7) then
          line = 7
          indexu(line) = 5 + 1
          indexl(line) = 4 + 1
          a(indexu(line),indexl(line)) = 1.6481e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 1097364.7910e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 312 -> 221	transition 6 states 6 4
        if (nli .ge. 6) then
          line = 6
          indexu(line) = 5 + 1
          indexl(line) = 3 + 1
          a(indexu(line),indexl(line)) = 2.6698e-3
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 1153126.8220e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 321 -> 312	transition 9 states 7 6
        if (nli .ge. 9) then
          line = 9
          indexu(line) = 6 + 1
          indexl(line) = 5 + 1
          a(indexu(line),indexl(line)) = 2.2893e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 1162911.5930e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 321 -> 212	transition 8 states 7 3
        if (nli .ge. 8) then
          line = 8
          indexu(line) = 6 + 1
          indexl(line) = 2 + 1
          a(indexu(line),indexl(line)) = 3.3380e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 3977046.4810e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 414 -> 321	transition 11 states 8 7
        if (nli .ge. 11) then
          line = 11
          indexu(line) = 7 + 1
          indexl(line) = 6 + 1
          a(indexu(line),indexl(line)) = 3.0593e-5
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
c this low frequency has been double checked in jpl catalog 
          freq(indexu(line),indexl(line)) = 380197.3720e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 414 -> 303	transition 10 states 8 5
        if (nli .ge. 10) then
          line = 10
          indexu(line) = 7 + 1
          indexl(line) = 4 + 1
          a(indexu(line),indexl(line)) = 2.4780e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 2640473.8360e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 330 -> 321	transition 14 states 9 7
        if (nli .ge. 14) then
          line = 14
          indexu(line) = 8 + 1
          indexl(line) = 6 + 1
          a(indexu(line),indexl(line)) = 6.65e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 2196345.7560e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 330 -> 303	transition 13 states 9 5
        if (nli .ge. 13) then
          line = 13
          indexu(line) = 8 + 1
          indexl(line) = 4 + 1
          a(indexu(line),indexl(line)) = 8.40e-3
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 4456621.9840e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 330 -> 221	transition 12 states 9 4
        if (nli .ge. 12) then
          line = 12
          indexu(line) = 8 + 1
          indexl(line) = 3 + 1
          a(indexu(line),indexl(line)) = 1.25e-0
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 4512384.1210e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 423 -> 330	transition 17 states 10 9
        if (nli .ge. 17) then
          line = 17
          indexu(line) = 9 + 1
          indexl(line) = 8 + 1
          a(indexu(line),indexl(line)) = 5.26e-5
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 448001.0750e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 423 -> 414	transition 16 states 10 8
        if (nli .ge. 16) then
          line = 16
          indexu(line) = 9 + 1
          indexl(line) = 7 + 1
          a(indexu(line),indexl(line)) = 8.09e-2
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 2264149.6500e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

c 423 -> 312	transition 15 states 10 6
        if (nli .ge. 15) then
          line = 15
          indexu(line) = 9 + 1
          indexl(line) = 5 + 1
          a(indexu(line),indexl(line)) = 4.83e-1
          a(indexl(line),indexu(line)) = 
     &          a(indexu(line),indexl(line)) 
          freq(indexu(line),indexl(line)) = 3807258.4120e6
          freq(indexl(line),indexu(line)) = 
     &        freq(indexu(line),indexl(line)) 
        endif

        line = 1
c        write(6,*)'h2o_init: lin,indx_u,indx_l,freqd,frequ ',
c     &        line, indexu(line), indexl(line),
c     &         freq(indexu(line),indexl(line)),
c     &          freq(indexl(line),indexu(line))

c Here set up the molecule's hyperfine structure for no hyperfine
c splitting. In this case the line will have one "hyperfine" whose relative
c intensity is 1 and whose velocity is not shifted

c These values will be reset by hyperfine subroutine if the
c molecule has hyperfine structure.

        do 411 line = 1,nlines
          nhyp(line) = 1
 411    continue

        do 408 line = 1,nlines
          do 409 ihyp = 2,maxhyp
            relint(ihyp,line) = 0.
            hypvel(ihyp,line) = 0.
 409      continue
          relint(1,line) = 1.
          hypvel(1,line) = 0.
 408    continue


        do 440 line = 1,nlines
           istartline(line) = line - 1
           n = indexu(line)
           m = indexl(line)
           aline(line) = a(n,m)
           freqline(line) = freq(n,m)
           glline(line) = statdg(m)
           guline(line) = statdg(n)
           lower(line)  = line - 1
           upper(line)  = line 
 440    continue

        lun = 22
c        msg = 'matrix of Einstein A coefficients'
c        call dmattmt(a,msg,nstate,lun)
c        msg = 'matrix of frequencies'
c        call dmattmt(freq,msg,nstate,lun)

        end

c23456789112345678921234567893123456789412345678951234567896123456789712

        subroutine rotor(istatus)

C This subroutine sets up some constants for transitions of
c dipole molecules, the Frequencies, Einstein A's,
C
C
C Since the linear molecules are so similar
C it suffices to select the proper moment of inertia in the
C subroutines frqrotor or arotor.
c This routine works for transitions of dipole molecules which are not split
c so that there are no hyperfine lines.
c Store the angular momentum and/or the parity of the molecule in JINDX
c KINDX, and LINDX.

c This subroutine loads up the molecular line data into freq and a
c that are indexed by state. Depends on copyatomic to transfer this
c data to the variables freqline and aline that are indexed by line.

        include 'nlines_f77.h'
        parameter (maxch=2000)
        parameter (maxhyp=50)
        parameter (maxtrans=500)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,aline,freqline,guline,glline

        double precision gco
        integer upper,lower


        logical docont,doprint,doprint2,nh3tau

        common /procid/ myid
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /lines/ indexu(nlines),indexl(nlines)
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /nltehyp/ upper(maxtrans),lower(maxtrans),
     &          istartline(nlines)

        istatus = 1

c All states have the same kindx, lindx. The levels will be
c indexed as quantum number + 1. So the ground state level
c which is zero, will be the first state, i=1, jindx(1) = 0.
        do 731 i = 1,nstate
                jindx(i) = i-1
                kindx(i) = 0
                lindx(i) = 0
 731    continue

c JINDX(LINE) is the quantum number for line number LINE
c INDEXU(LINE) is the Fortran array index number for the line

        do 4 line = 1,nlines
                indexu(line) = line+1
                indexl(line) = line
 4      continue

c For CO the dipole is about 0.1 debye or 0.1x10^-18 esu
c Estimates range from 0.10 to 0.128 in "Tables of
c Electric Dipole Moments", L.G. Wesson, 1948, The
c Technology Press, MIT
c The value of 0.1121 yields the Einstein A values
c listed in Martin and Barrett, 1978, ApJS, 36, 1
c The JPL catalogue dipole moment is 0.11011

c For 13CO JPL has 0.11046


c The dipole moment of H12CO+   is 3.3 Debye (JPL catalog)
c The dipole moment of H13CO+   is 3.3 Debye (JPL catalog)
c The dipole moment of 12C32S is 2.0 Debye (5-29-97)
c The dipole moment of N2H+   is 3.4 Debye (JPL catalog)
c The dipole moment of N2D+   is 3.4 Debye (JPL catalog)
c The dipole moment of C17O   is 0.11034 Debye (JPL catalog)
c The dipole moment of C18O   is 0.11079 Debye (JPL catalog)
c The dipole moment of HCN    is 2.984   Debye (JPL catalog)


c B constant for H12CO+   is 44.5944 GHz (JPL catalog)
c B constant for H13CO+   is 43.37732 GHz (JPL catalog)
c B constant for 12CO   is 57.8975 GHz.
c B constant for 13CO   is 55.3449 GHz
c B constant for 12C32S is 24584.352 MHz
c B constant for N2H+   is 46586.867 MHz (JPL catalog)
c B constant for N2D+   is 38554.719 MHz (JPL catalog)
c B constant for C17O   is 56.179990 GHz (JPL catalog)
c B constant for C18O   is 54.891420 GHz (JPL catalog)
c B constant for HCN    is 44.315975 GHz (JPL catalog)


c Checks: B for 12CO from Flower and Launey 1985 MNRAS 214, 271
c 1.92265 cm^-1 = 0.520115 cm = 57679.5 MHz
c 2*57679.5 MHz = 115.359 GHz, 12CO(1-0) = 115.2712 GHz

c B constants from JPL catalogue
c 12CO 57635.968
c 13CO 55101.011

        if (myid .lt. 2) write(22,*) 'Setup for molecule ',molecule

c None of the variables, atoms, brot, or dipole are double precision
c despite the assignments below with doubles on RHS

c 12CO
        if (molecule .eq. 7) then
                atoms = 28.d0
                brot = 57635.968d6
                dipole = 0.11011d-18
        endif

c CS
        if (molecule .eq. 8) then
                atoms = 44.d0
                brot = 24584.352d6
                dipole = 2.0d-18
        endif

c N2H+
        if (molecule .eq. 9) then
                atoms = 29.d0
                brot = 46586.867d6
                dipole = 3.4d-18
        endif

c 13CO
        if (molecule .eq. 10) then
                atoms = 29.d0
                brot = 55101.011d6
                dipole = 0.11046d-18
        endif

c H13CO+
        if (molecule .eq. 11) then
                atoms = 30.d0
                brot = 4.337732d10
                dipole = 3.3d-18
        endif

c H12CO+
        if (molecule .eq. 12) then
                atoms = 29.d0
                brot = 4.4594d10
                dipole = 3.3d-18
        endif

c C17O
        if (molecule .eq. 13) then
                atoms = 29.d0
                brot = 56.179990d9
                dipole = 0.11034d-18
        endif
c C18O
        if (molecule .eq. 14) then
                atoms = 30.d0
                brot = 54.891420d9
                dipole = 0.11079d-18
        endif
c HCN
        if (molecule .eq. 15) then
                atoms = 27.d0
                brot = 44.315975d9
                dipole = 2.984d-18
        endif
c N2D+
        if (molecule .eq. 16) then
                atoms = 30.d0
                brot = 38554.719d6
                dipole = 3.4d-18
        endif
c SiO
        if (molecule .eq. 20) then
                atoms = 44.d0
                brot = 21711.967d6
                dipole = 3.098d-18
        endif



        if (.not.(molecule .eq. 7 .or.
     &      molecule .eq. 8 .or.
     &      molecule .eq. 9 .or.
     &      molecule .eq. 10 .or.
     &      molecule .eq. 11 .or.
     &      molecule .eq. 12 .or.
     &      molecule .eq. 13 .or.
     &      molecule .eq. 14 .or.
     &      molecule .eq. 15 .or.
     &      molecule .eq. 16 .or.
     &      molecule .eq. 20 )) then
            if (myid .lt. 2) then
              write(22,*) 'Subroutine rotor: No values',
     &          ' for this transition. (molecule = ',molecule,')'
              write(22,*) 'No B value',
     &           ' for this transition. (molecule = ',molecule,')'
              write(22,*) 'No dipole value',
     &          ' for this transition. (molecule = ',molecule,')'
            endif
            istatus = 0
            return
        endif

        if (myid .lt. 2) write(22,*) 'molecular weight ',atoms
        if (myid .lt. 2) write(22,*) 'B rotational constant ',brot
        if (myid .lt. 2) write(22,*) 'Dipole moment ',dipole

c Compute the energy levels for all transition, store as frequency
c Compute the Einstein A's for all transitions
c Load the constants used to determine the collision rates

                call frqrotor(istatus)
                if (istatus .eq. 0) return

        if (myid .lt. 2) 
     &    write(22,*) 'Here are the transition frequencies'

                call arotor

       do 812 kk = 1,nstate
               statdg(kk) = gco(kk)
 812   continue

        if (myid .lt. 2) then
        write(22,*) ' '
        if(molecule .eq. 13) write(22,*) 'C17O      frequency GHz  ',
     &          '   Einstein A','        stat. deg. lower state'
        if(molecule .eq. 15) write(22,*) 'HCN       frequency GHz  ',
     &          '   Einstein A','        stat. deg. lower state'
        if(molecule .eq. 14) write(22,*) 'C18O      frequency GHz  ',
     &          '   Einstein A','        stat. deg. lower state'
        if(molecule .eq. 10) write(22,*) '13CO      frequency GHz  ',
     &          '   Einstein A','        stat. deg. lower state'
        if(molecule .eq. 11) write(22,*) 'H13CO+    frequency GHz  ',
     &          '   Einstein A','        stat. deg. lower state'
        if(molecule .eq. 12) write(22,*) 'H12CO+    frequency GHz  ',
     &          '   Einstein A','        stat. deg. lower state'
        if(molecule .eq. 7) write(22,*) '12CO      frequency GHz  ',
     &          '   Einstein A','        stat. deg. lower state',nstate
        if(molecule .eq. 8) write(22,*) 'CS        frequency GHz  ',
     &          '   Einstein A','        stat. deg. lower state'
        if(molecule .eq. 9) write(22,*) 'N2H+      frequency GHz  ',
     &          '   Einstein A','        stat. deg. lower state'
        if(molecule .eq. 16) write(22,*) 'N2D+      frequency GHz  ',
     &          '   Einstein A','        stat. deg. lower state'
        if(molecule .eq. 20) write(22,*) 'SiO+      frequency GHz  ',
     &          '   Einstein A','        stat. deg. lower state'
        do 813 kk = 1,nlines
                n = indexu(kk )
                m = indexl(kk )
                write (22,1002) jindx(n),jindx(m),freq(m,n)/1.d9,
     &                  a(m,n),int(statdg(n)),int(statdg(m))
 1002    format(1x,i2,'-->',i2,2x,f15.9,3x,1pe16.9,3x,i5,3x,i5 )
 813    continue 
        endif


c Here set up the molecule's hyperfine structure for no hyperfine
c splitting. In this case the line will have one "hyperfine" whose relative
c intensity is 1 and whose velocity is not shifted

c These values will be reset by hyperfine subroutine if the
c molecule has hyperfine structure.

       do 411 line = 1,nlines
          nhyp(line) = 1
 411    continue

        do 408 line = 1,nlines
          do 409 ihyp = 2,maxhyp
            relint(ihyp,line) = 0.
            hypvel(ihyp,line) = 0.
 409      continue
          relint(1,line) = 1.
          hypvel(1,line) = 0.
 408    continue

        do 440 line = 1,nlines
           istartline(line) = line - 1
           n = indexu(line)
           m = indexl(line)
           aline(line) = a(n,m)
           freqline(line) = freq(n,m)
           glline(line) = statdg(m)
           guline(line) = statdg(n)
           lower(line)  = line - 1
           upper(line)  = line 
 440    continue



        if (myid .lt. 2) write(22,*) 'Finished loading dipole constants'
        return
        end

        subroutine frqrotor(istatus)

c This routine computes the energy difference between the rotational
c levels in units of frequency. The energy  of the rotational
c states, W, is computed from a formula (eqn. 1-8 of Townes and Schawlow).
c The results are stored in COMMON /AFRC/.


        include 'nlines_f77.h'
        parameter (ndim=nstate-1)

        integer u,l

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        double precision flag

        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)
        common /frqinv/ frqinv(nstate),sgn(nstate)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /lines/ indexu(nlines),indexl(nlines)


        flag = 0.d0
        istatus = 1

c Initialize the array of frequencies
        do 1 u = 1,nstate
        do 1 l = 1,nstate
                freq(u,l) = flag
 1      continue

c Simple calculation for frequencies.
        do 2 u = 1,nstate
        do 2 l = 1,u-1

                freq(u,l) = dble (brot*( jindx(u)*(jindx(u) + 1)
     &                          - jindx(l)*(jindx(l) + 1)  ) )
 2      continue

c 12CO frequencies from JPL spectral line catalogue.
c These are more accurate than the 2B(J+1) values
c These values overwrite the calculated values for the
c allowed dipole transitions. Remember, the indices
c on freq are array indices, not quantum numbers but q.n.+1.

c Use a non-parameter variable to avoid warning messages about
c no path to the following if statements
        nst = nstate

        if (molecule .eq. 7) then
                if (nst .ge.  2) then
                    i = 2
                    freq(i,i-1)   =  115271.2018 d6
                endif
                if (nst .ge.  3) then
                    i = 3
                    freq(i,i-1)   =  230538.0000 d6
                endif
                if (nst .ge.  4) then
                    i = 4
                    freq(i,i-1)   =  345795.9899 d6
                endif
                if (nst .ge.  5) then
                    i = 5
                    freq(i,i-1)   =  461040.7682 d6
                endif
                if (nst .ge.  6) then
                    i = 6
                    freq(i,i-1)   =  576267.9305 d6
                endif
                if (nst .ge.  7) then
                    i = 7
                    freq(i,i-1)   =  691473.0763 d6
                endif
                if (nst .ge.  8) then
                    i = 8
                    freq(i,i-1)   =  806651.8060 d6
                endif
                if (nst .ge.  9) then
                    i = 9
                    freq(i,i-1)   =  921799.7000 d6
                endif
                if (nst .ge. 10) then
                    i = 10
                    freq(i,i-1)  = 1036912.3930 d6
                endif
                if (nst .ge. 11) then
                    i = 11
                    freq(i,i-1) = 1151985.4520 d6
                endif
                if (nst .ge. 12) then
                    i = 12
                    freq(i,i-1) = 1267014.4860 d6
                endif
                if (nst .ge. 13) then
                    i = 13
                    freq(i,i-1) = 1381995.1050 d6
                endif
                if (nst .ge. 14) then
                    i = 14
                    freq(i,i-1) = 1496922.9090 d6
                endif
                if (nst .ge. 15) then
                    i = 15
                    freq(i,i-1) = 1611793.5180 d6
                endif
                if (nst .ge. 16) then
                    i = 16
                    freq(i,i-1) = 1726602.5066 d6
                endif
        endif

c N2H+ frequencies from the JPL catalog

        if (molecule .eq. 9) then
                if (nst .ge.  2) then
                    i = 2
                    freq(i,i-1)   =   93173.7000 d6
                endif
                if (nst .ge.  3) then
                    i = 3
                    freq(i,i-1)   =  186344.7710 d6
                endif
                if (nst .ge.  4) then
                    i = 4
                    freq(i,i-1)   =  279511.7010 d6
                endif
                if (nst .ge.  5) then
                    i = 5
                    freq(i,i-1)   =  372672.5090 d6
                endif
                if (nst .ge.  6) then
                    i = 6
                    freq(i,i-1)   =  465824.9470 d6
                endif
                if (nst .ge.  7) then
                    i = 7
                    freq(i,i-1)   =  558966.6656 d6
                endif
                if (nst .ge.  8) then
                    i = 8
                    freq(i,i-1)   =  652095.8651 d6
                endif
                if (nst .ge.  9) then
                    i = 9
                    freq(i,i-1)   =  745210.3361 d6
                endif
                if (nst .ge. 10) then
                    i = 10
                    freq(i,i-1)  =  838307.9748 d6
                endif
                if (nst .ge. 11) then
                    i = 11
                    freq(i,i-1) =  931386.6769 d6
                endif
                if (nst .ge. 12) then
                    i = 12
                    freq(i,i-1) = 1024444.2000 d6
                endif
                if (nst .ge. 13) then
                    i = 13
                    freq(i,i-1) = 1117478.8554 d6
                endif
                if (nst .ge. 14) then
                    i = 14
                    freq(i,i-1) = 1210488.1238 d6
                endif
                if (nst .ge. 15) then
                    i = 15
                    freq(i,i-1) = 1303470.0394 d6
                endif
                if (nst .ge. 16) then
                    i = 16
                    freq(i,i-1) = 1396422.4982 d6
                endif
        endif

        if (molecule .eq. 10) then
                if (nst .ge.  2) then
                    i = 2
                    freq(i,i-1)   =  110201.3541 d6
                endif
                if (nst .ge.  3) then
                    i = 3
                    freq(i,i-1)   =  220398.6765 d6
                endif
                if (nst .ge.  4) then
                    i = 4
                    freq(i,i-1)   =  330587.9601 d6
                endif
                if (nst .ge.  5) then
                    i = 5
                    freq(i,i-1)   =  440765.1668 d6
                endif
                if (nst .ge.  6) then
                    i = 6
                    freq(i,i-1)   =  550926.3029 d6
                endif
                if (nst .ge.  7) then
                    i = 7
                    freq(i,i-1)   =  661067.2801 d6
                endif
                if (nst .ge.  8) then
                    i = 8
                    freq(i,i-1)   =  771184.1376 d6
                endif
                if (nst .ge.  9) then
                    i = 9
                    freq(i,i-1)   =  881271.8339 d6
                endif
                if (nst .ge. 10) then
                    i = 10
                    freq(i,i-1)   =  991329.3479 d6
                endif
        endif


c H12CO+ frequencies from JPL catalog

        if (molecule .eq. 12) then
                if (nst .ge.  2) then
                    i = 2
                    freq(i,i-1)   =   89188.5230 d6
                endif
                if (nst .ge.  3) then
                    i = 3
                    freq(i,i-1)   =  178375.0650 d6
                endif
                if (nst .ge.  4) then
                    i = 4
                    freq(i,i-1)   =  267557.6190 d6
                endif
                if (nst .ge.  5) then
                    i = 5
                    freq(i,i-1)   =  356734.2880 d6
                endif
                if (nst .ge.  6) then
                    i = 6
                    freq(i,i-1)   =  445902.9960 d6
                endif
                if (nst .ge.  7) then
                    i = 7
                    freq(i,i-1)   =  535061.7755 d6
                endif
                if (nst .ge.  8) then
                    i = 8
                    freq(i,i-1)   =  623882.9970 d6
                endif
                if (nst .ge.  9) then
                    i = 9
                    freq(i,i-1)   =  713342.0900 d6
                endif
                if (nst .ge. 10) then
                    i = 10
                    freq(i,i-1)  =  802458.3290 d6
                endif
                if (nst .ge. 11) then
                    i = 11
                    freq(i,i-1) =  891557.9242 d6
                endif
                if (nst .ge. 12) then
                    i = 12
                    freq(i,i-1) =  983924.5097 d6
                endif
                if (nst .ge. 13) then
                    i = 13
                    freq(i,i-1) = 1069693.8000 d6
                endif
                if (nst .ge. 14) then
                    i = 14
                    freq(i,i-1) = 1158727.6478 d6
                endif
                if (nst .ge. 15) then
                    i = 15
                    freq(i,i-1) = 1247734.8251 d6
                endif
                if (nst .ge. 16) then
                    i = 16
                    freq(i,i-1) = 1336713.8728 d6
                endif
        endif

        if (molecule .eq. 13) then
                if (nst .ge.  2) then 
                    i = 2
                    freq(i,i-1)   =  112359.2837 d6
                endif
                if (nst .ge.  3) then
                    i = 3
                    freq(i,i-1)   =  224714.3850 d6
                endif
                if (nst .ge.  4) then
                    i = 4
                    freq(i,i-1)   =  337061.1298 d6
                endif
                if (nst .ge.  5) then
                    i = 5
                    freq(i,i-1)   =  449395.3412 d6
                endif
                if (nst .ge.  6) then
                    i = 6
                    freq(i,i-1)   =  561712.7845 d6
                endif
                if (nst .ge.  7) then
                    i = 7
                    freq(i,i-1)   =  674009.3443 d6
                endif
                if (nst .ge.  8) then
                    i = 8
                    freq(i,i-1)   =  786280.8166 d6
                endif
                if (nst .ge.  9) then
                    i = 9
                    freq(i,i-1)   =  898523.0217 d6
                endif
        endif


        if (molecule .eq. 14) then
                if (nst .ge.  2) then
                    i = 2
                    freq(i,i-1)   =  109782.1734 d6
                endif
                if (nst .ge.  3) then
                    i = 3
                    freq(i,i-1)   =  219560.3568 d6
                endif
                if (nst .ge.  4) then
                    i = 4
                    freq(i,i-1)   =  329330.5453 d6
                endif
                if (nst .ge.  5) then
                    i = 5
                    freq(i,i-1)   =  439088.7631 d6
                endif
                if (nst .ge.  6) then
                    i = 6
                    freq(i,i-1)   =  548830.9775 d6
                endif
                if (nst .ge.  7) then
                    i = 7
                    freq(i,i-1)   =  658553.2728 d6
                endif
                if (nst .ge.  8) then
                    i = 8
                    freq(i,i-1)   =  768251.5890 d6
                endif
                if (nst .ge.  9) then
                    i = 9
                    freq(i,i-1)   =  877921.9561 d6
                endif
                if (nst .ge. 10) then
                    i = 10
                    freq(i,i-1)   =  987560.3837 d6
                endif
        endif

c HCN  frequencies from the JPL catalog
        if (molecule .eq. 15) then
                if (nst .ge.  2) then
                    i = 2
                    freq(i,i-1)   =   88631.8470 d6
                endif
                if (nst .ge.  3) then
                    i = 3
                    freq(i,i-1)   =  177261.1100 d6
                endif
                if (nst .ge.  4) then
                    i = 4
                    freq(i,i-1)   =  265886.1800 d6
                endif
                if (nst .ge.  5) then
                    i = 5
                    freq(i,i-1)   =  354505.4759 d6
                endif
                if (nst .ge.  6) then
                    i = 6
                    freq(i,i-1)   =  443116.1554 d6
                endif
                if (nst .ge.  7) then
                    i = 7
                    freq(i,i-1)   =  531716.3875 d6
                endif
                if (nst .ge.  8) then
                    i = 8
                    freq(i,i-1)   =  620304.0952 d6
                endif
                if (nst .ge.  9) then
                    i = 9
                    freq(i,i-1)   =  708877.2081 d6
                endif
                if (nst .ge. 10) then
                    i = 10
                    freq(i,i-1)  =  797433.6638 d6
                endif
                if (nst .ge. 11) then
                    i = 11
                    freq(i,i-1) =  885971.4087 d6
                endif
                if (nst .ge. 12) then
                    i = 12
                    freq(i,i-1) =  974488.4000 d6
                endif
        endif

c N2D+ frequencies from the JPL catalog

        if (molecule .eq. 16) then
                if (nst .ge.  2) then
                    i = 2
                    freq(i,i-1)   =   77109.6100 d6
                endif
                if (nst .ge.  3) then
                    i = 3
                    freq(i,i-1)   =  154217.0960 d6
                endif
                if (nst .ge.  4) then
                    i = 4
                    freq(i,i-1)   =  231321.6650 d6
                endif
                if (nst .ge.  5) then
                    i = 5
                    freq(i,i-1)   =  308422.2220 d6
                endif
                if (nst .ge.  6) then
                    i = 6
                    freq(i,i-1)   =  385516.7620 d6
                endif
                if (nst .ge.  7) then
                    i = 7
                    freq(i,i-1)   =  462603.9320 d6
                endif
                if (nst .ge.  8) then
                    i = 8
                    freq(i,i-1)   =  539682.4810 d6
                endif
                if (nst .ge.  9) then
                    i = 9
                    freq(i,i-1)   =  616750.7155 d6
                endif
                if (nst .ge. 10) then
                    i = 10
                    freq(i,i-1)  =  693807.2448 d6
                endif
                if (nst .ge. 11) then
                    i = 11
                    freq(i,i-1) =  770850.6058 d6
                endif
                if (nst .ge. 12) then
                    i = 12
                    freq(i,i-1) =  847879.3354 d6
                endif
                if (nst .ge. 13) then
                    i = 13
                    freq(i,i-1) =  924891.9705 d6
                endif
                if (nst .ge. 14) then
                    i = 14
                    freq(i,i-1) = 1001887.0478 d6
                endif
                if (nst .ge. 15) then
                    i = 15
                    freq(i,i-1) = 1078863.1043 d6
                endif
                if (nst .ge. 16) then
                    i = 16
                    freq(i,i-1) = 1155818.6768 d6
                endif
        endif

c SiO+ frequencies from the JPL catalog

        if (molecule .eq. 20) then
                if (nst .ge.  2) then
                    i = 2
                    freq(i,i-1)   =   43423.7600 d6
                endif
                if (nst .ge.  3) then
                    i = 3
                    freq(i,i-1)   =   86846.9660 d6
                endif
                if (nst .ge.  4) then
                    i = 4
                    freq(i,i-1)   =  130268.6100 d6
                endif
                if (nst .ge.  5) then
                    i = 5
                    freq(i,i-1)   =  173688.3100 d6
                endif
                if (nst .ge.  6) then
                    i = 6
                    freq(i,i-1)   =  217104.9800 d6
                endif
                if (nst .ge.  7) then
                    i = 7
                    freq(i,i-1)   =  260518.0200 d6
                endif
                if (nst .ge.  8) then
                    i = 8
                    freq(i,i-1)   =  303926.9600 d6
                endif
                if (nst .ge.  9) then
                    i = 9
                    freq(i,i-1)   =  347330.6310 d6
                endif
                if (nst .ge. 10) then
                    i = 10
                    freq(i,i-1)  =  390728.4483 d6
                endif
                if (nst .ge. 11) then
                    i = 11
                    freq(i,i-1) =  434119.5521 d6
                endif
                if (nst .ge. 12) then
                    i = 12
                    freq(i,i-1) =  477503.0965 d6
                endif
                if (nst .ge. 13) then
                    i = 13
                    freq(i,i-1) =  520878.2390 d6
                endif
                if (nst .ge. 14) then
                    i = 14
                    freq(i,i-1) =  564203.9620 d6
                endif
                if (nst .ge. 15) then
                    i = 15
                    freq(i,i-1) =  607599.4207 d6
                endif
                if (nst .ge. 16) then
                    i = 16
                    freq(i,i-1) =  650943.5888 d6
                endif
        endif


c                do 103 line = 1,nlines
c                      freqline(line) = freq(indexu(line),indexl(line))
c                      write(6,*) 'frqrotor ',
c     &                    line,indexu(line),indexl(line),freqline(line)
c                      write(22,*) 'frqrotor ',
c     &                    line,indexu(line),indexl(line),freqline(line)
c 103            continue


        iflag = 1
c Since the frequency may be used to compute the energy separation
c fill in the reverse half of the array.
        do 6 u = 1,nstate
        do 6 l = u+1,nstate
                freq(u,l) = freq(l,u)
 6      continue

        do 4 u = 1,nstate
        do 4 l = 1,nstate

        if (( u .ne. l) .and. freq(u,l) .lt. 1.) then
                iflag = 0
        write (6,100) jindx(l),jindx(u) ,freq(u,l),brot
        if (myid .lt. 2) 
     &      write (22,100) jindx(l),jindx(u),freq(u,l),brot
 100            format(1x, 'freq missing: ',i3,'<--',i3,
     &                  2x,1pe16.9,' rot const ',e13.6)
        istatus = 0
        return
        endif
 4      continue
        if (iflag .eq. 0) return


        return
        end

        subroutine arotor

c This subroutine computes the Einstein A coefficients of the allowed
c transitions of dipole molecules.

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        double precision flag,factor,coeff

        integer u,l

        include 'nlines_f77.h'
        parameter (ndim=nstate-1)

        common /procid/ myid
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav

        flag = 0.d0

c Initialize the array of einstein A coefficients to value for
c disallowed transitions. Use zero.
        do 2 u = 1,nstate
        do 2 l = 1,nstate
                a(u,l) = flag
 2      continue


c Downward transitions state M to state N. Only transitions for
c DeltaJ = 1 are allowed. The angular momenta of
c states U and L are stored in integer arrays JINDX(U) and JINDX(L)
c The difference in energy expressed in Hz is stored in FREQ(U,L)

c The equation for A is in Rybicki and Lightman 10-28b, pg 274.
c The factor 1/g S |d|^2 in this equation, which is the
c average of the dipole matrix elements for the rotational
c transitions of a diatomic molecule is in Townes and
c Schawlow 1-76, pg 23. This is called factor below.


        do 1 u = 2,nstate

                l = u - 1

                factor = dble(jindx(l) + 1)/dble(2*jindx(l) + 3)

c       if (myid .lt. 2) write(22,*) 'factor',factor

                coeff = 64.d0*pi**4*freq(u,l)*dble(dipole)*freq(u,l)
     &                  *dble(dipole)*freq(u,l)/(3.d0*c**2*planck*c)

c       write(22,*) pi,freq(u,l),c,planck
c       write(22,*) 'coeff',coeff

                a(l,u) = coeff * factor
                a(u,l) = a(l,u)


 1      continue


        return
        end



        function gco(j)
        double precision gco
        include 'nlines_f77.h'
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)

c Statistical degeneracy of CO. Remember that the index J is
c not the quantum number, which is JINDX(J). For the case
c of CO and other linear molecules, JINDX(J) = J-1

        gco = 2.d0*dble(jindx(j)) + 1.d0

        return
        end



        function n2dphyp(vwmin)

        include 'nlines_f77.h'
        parameter (maxhyp=50)

        common /procid/ myid
        double precision hypfrq(maxhyp,nlines),frq
        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &      grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        dimension tmp1(maxhyp),tmp2(maxhyp),tmp3(maxhyp)
        dimension hypfrqn2d(maxhyp,8),nhypn2d(8),relintn2d(maxhyp,8)

c        write(6,*) 'Working out N2D+ hyperfine structure'
        n2dphyp = 1

c These are in GHz
        hypfrqn2d( 1, 1) =   77.107498
        relintn2d( 1, 1) = 0.0370
        hypfrqn2d( 2, 1) =   77.107798 
        relintn2d( 2, 1) = 0.1852
        hypfrqn2d( 3, 1) =   77.107934 
        relintn2d( 3, 1) = 0.1111
        hypfrqn2d( 4, 1) =   77.109338 
        relintn2d( 4, 1) = 0.1852
        hypfrqn2d( 5, 1) =   77.109632 
        relintn2d( 5, 1) = 0.2593
        hypfrqn2d( 6, 1) =   77.109830 
        relintn2d( 6, 1) = 0.1111
        hypfrqn2d( 7, 1) =   77.112130 
        relintn2d( 7, 1) = 0.1111
        nhypn2d( 1) =  7

c These are from Dore's tables and are in MHz. Divide by 1000 later
       hypfrqn2d(1,2) = 154214.845
       relintn2d(1,2) = 0.002
       hypfrqn2d(2,2) = 154215.040
       relintn2d(2,2) = 0.004
       hypfrqn2d(3,2) = 154215.168
       relintn2d(3,2) = 0.020
       hypfrqn2d(4,2) = 154215.288
       relintn2d(4,2) = 0.039
       hypfrqn2d(5,2) = 154215.331
       relintn2d(5,2) = 0.011
       hypfrqn2d(6,2) = 154215.425
       relintn2d(6,2) = 0.036
       hypfrqn2d(7,2) = 154215.579
       relintn2d(7,2) = 0.005
       hypfrqn2d(8,2) = 154215.628
       relintn2d(8,2) = 0.062
       hypfrqn2d(9,2) = 154215.653
       relintn2d(9,2) = 0.003
       hypfrqn2d(10,2) = 154215.888
       relintn2d(10,2) = 0.013
       hypfrqn2d(11,2) = 154216.754
       relintn2d(11,2) = 0.066
       hypfrqn2d(12,2) = 154216.818
       relintn2d(12,2) = 0.013
       hypfrqn2d(13,2) = 154216.889
       relintn2d(13,2) = 0.027
       hypfrqn2d(14,2) = 154217.076
       relintn2d(14,2) = 0.012
       hypfrqn2d(15,2) = 154217.109
       relintn2d(15,2) = 0.143
       hypfrqn2d(16,2) = 154217.131
       relintn2d(16,2) = 0.088
       hypfrqn2d(17,2) = 154217.137
       relintn2d(17,2) = 0.112
       hypfrqn2d(18,2) = 154217.206
       relintn2d(18,2) = 0.200
       hypfrqn2d(19,2) = 154217.502
       relintn2d(19,2) = 0.031
       hypfrqn2d(20,2) = 154217.617
       relintn2d(20,2) = 0.021
       hypfrqn2d(21,2) = 154218.120
       relintn2d(21,2) = 0.003
       hypfrqn2d(22,2) = 154219.631
       relintn2d(22,2) = 0.010
       hypfrqn2d(23,2) = 154219.766
       relintn2d(23,2) = 0.013
       hypfrqn2d(24,2) = 154219.834
       relintn2d(24,2) = 0.014
       hypfrqn2d(25,2) = 154219.969
       relintn2d(25,2) = 0.031
       hypfrqn2d(26,2) = 154220.058
       relintn2d(26,2) = 0.006
       hypfrqn2d(27,2) = 154220.094
       relintn2d(27,2) = 0.008
       nhypn2d(2) =       27


       hypfrqn2d(1,3) = 231319.945
       relintn2d(1,3) = 0.010
       hypfrqn2d(2,3) = 231319.995
       relintn2d(2,3) = 0.017
       hypfrqn2d(3,3) = 231320.017
       relintn2d(3,3) = 0.007
       hypfrqn2d(4,3) = 231321.253
       relintn2d(4,3) = 0.015
       hypfrqn2d(5,3) = 231321.456
       relintn2d(5,3) = 0.048
       hypfrqn2d(6,3) = 231321.499
       relintn2d(6,3) = 0.022
       hypfrqn2d(7,3) = 231321.530
       relintn2d(7,3) = 0.007
       hypfrqn2d(8,3) = 231321.547
       relintn2d(8,3) = 0.014
       hypfrqn2d(9,3) = 231321.617
       relintn2d(9,3) = 0.094
       hypfrqn2d(10,3) = 231321.795
       relintn2d(10,3) = 0.089
       hypfrqn2d(11,3) = 231321.908
       relintn2d(11,3) = 0.061
       hypfrqn2d(12,3) = 231321.914
       relintn2d(12,3) = 0.125
       hypfrqn2d(13,3) = 231321.918
       relintn2d(13,3) = 0.136
       hypfrqn2d(14,3) = 231321.919
       relintn2d(14,3) = 0.101
       hypfrqn2d(15,3) = 231321.961
       relintn2d(15,3) = 0.019
       hypfrqn2d(16,3) = 231321.966
       relintn2d(16,3) = 0.175
       hypfrqn2d(17,3) = 231322.231
       relintn2d(17,3) = 0.007
       hypfrqn2d(18,3) = 231322.426
       relintn2d(18,3) = 0.010
       hypfrqn2d(19,3) = 231324.012
       relintn2d(19,3) = 0.002
       hypfrqn2d(20,3) = 231324.086
       relintn2d(20,3) = 0.002
       hypfrqn2d(21,3) = 231324.334
       relintn2d(21,3) = 0.011
       hypfrqn2d(22,3) = 231324.450
       relintn2d(22,3) = 0.014
       hypfrqn2d(23,3) = 231324.517
       relintn2d(23,3) = 0.004
       hypfrqn2d(24,3) = 231324.698
       relintn2d(24,3) = 0.002
       hypfrqn2d(25,3) = 231324.839
       relintn2d(25,3) = 0.002
       nhypn2d(3) =       25


       hypfrqn2d(1,4) = 308422.092
       relintn2d(1,4) = 0.069
       hypfrqn2d(2,4) = 308422.094
       relintn2d(2,4) = 0.048
       hypfrqn2d(3,4) = 308422.157
       relintn2d(3,4) = 0.103
       hypfrqn2d(4,4) = 308422.233
       relintn2d(4,4) = 0.098
       hypfrqn2d(5,4) = 308422.263
       relintn2d(5,4) = 0.075
       hypfrqn2d(6,4) = 308422.290
       relintn2d(6,4) = 0.126
       hypfrqn2d(7,4) = 308422.296
       relintn2d(7,4) = 0.105
       hypfrqn2d(8,4) = 308422.301
       relintn2d(8,4) = 0.131
       hypfrqn2d(9,4) = 308422.331
       relintn2d(9,4) = 0.160
       hypfrqn2d(10,4) = 308420.289
       relintn2d(10,4) = 0.006
       hypfrqn2d(11,4) = 308420.320
       relintn2d(11,4) = 0.009
       hypfrqn2d(12,4) = 308420.332
       relintn2d(12,4) = 0.005
       hypfrqn2d(13,4) = 308421.728
       relintn2d(13,4) = 0.008
       hypfrqn2d(14,4) = 308421.865
       relintn2d(14,4) = 0.005
       hypfrqn2d(15,4) = 308421.866
       relintn2d(15,4) = 0.008
       hypfrqn2d(16,4) = 308422.599
       relintn2d(16,4) = 0.010
       hypfrqn2d(17,4) = 308422.699
       relintn2d(17,4) = 0.005
       hypfrqn2d(18,4) = 308422.804
       relintn2d(18,4) = 0.006
       hypfrqn2d(19,4) = 308424.631
       relintn2d(19,4) = 0.007
       hypfrqn2d(20,4) = 308424.693
       relintn2d(20,4) = 0.008
       hypfrqn2d(21,4) = 308424.703
       relintn2d(21,4) = 0.003
       nhypn2d(4) =       21


       hypfrqn2d(1,5) = 385516.647
       relintn2d(1,5) = 0.062
       hypfrqn2d(2,5) = 385516.653
       relintn2d(2,5) = 0.080
       hypfrqn2d(3,5) = 385516.689
       relintn2d(3,5) = 0.106
       hypfrqn2d(4,5) = 385516.733
       relintn2d(4,5) = 0.102
       hypfrqn2d(5,5) = 385516.741
       relintn2d(5,5) = 0.083
       hypfrqn2d(6,5) = 385516.768
       relintn2d(6,5) = 0.126
       hypfrqn2d(7,5) = 385516.772
       relintn2d(7,5) = 0.107
       hypfrqn2d(8,5) = 385516.779
       relintn2d(8,5) = 0.128
       hypfrqn2d(9,5) = 385516.800
       relintn2d(9,5) = 0.152
       hypfrqn2d(10,5) = 385514.735
       relintn2d(10,5) = 0.004
       hypfrqn2d(11,5) = 385514.757
       relintn2d(11,5) = 0.006
       hypfrqn2d(12,5) = 385514.764
       relintn2d(12,5) = 0.003
       hypfrqn2d(13,5) = 385516.224
       relintn2d(13,5) = 0.005
       hypfrqn2d(14,5) = 385516.309
       relintn2d(14,5) = 0.005
       hypfrqn2d(15,5) = 385516.314
       relintn2d(15,5) = 0.003
       hypfrqn2d(16,5) = 385517.154
       relintn2d(16,5) = 0.006
       hypfrqn2d(17,5) = 385517.207
       relintn2d(17,5) = 0.003
       hypfrqn2d(18,5) = 385517.275
       relintn2d(18,5) = 0.004
       hypfrqn2d(19,5) = 385519.050
       relintn2d(19,5) = 0.005
       hypfrqn2d(20,5) = 385519.086
       relintn2d(20,5) = 0.003
       hypfrqn2d(21,5) = 385519.091
       relintn2d(21,5) = 0.005
       nhypn2d(5) =       21


       hypfrqn2d(1,6) = 462603.777
       relintn2d(1,6) = 0.071
       hypfrqn2d(2,6) = 462603.785
       relintn2d(2,6) = 0.087
       hypfrqn2d(3,6) = 462603.809
       relintn2d(3,6) = 0.108
       hypfrqn2d(4,6) = 462603.839
       relintn2d(4,6) = 0.105
       hypfrqn2d(5,6) = 462603.839
       relintn2d(5,6) = 0.089
       hypfrqn2d(6,6) = 462603.863
       relintn2d(6,6) = 0.124
       hypfrqn2d(7,6) = 462603.866
       relintn2d(7,6) = 0.108
       hypfrqn2d(8,6) = 462603.874
       relintn2d(8,6) = 0.126
       hypfrqn2d(9,6) = 462603.890
       relintn2d(9,6) = 0.145
       hypfrqn2d(10,6) = 462601.801
       relintn2d(10,6) = 0.003
       hypfrqn2d(11,6) = 462601.819
       relintn2d(11,6) = 0.004
       hypfrqn2d(12,6) = 462601.823
       relintn2d(12,6) = 0.002
       hypfrqn2d(13,6) = 462603.320
       relintn2d(13,6) = 0.003
       hypfrqn2d(14,6) = 462603.380
       relintn2d(14,6) = 0.003
       hypfrqn2d(15,6) = 462603.387
       relintn2d(15,6) = 0.002
       hypfrqn2d(16,6) = 462604.278
       relintn2d(16,6) = 0.004
       hypfrqn2d(17,6) = 462604.312
       relintn2d(17,6) = 0.003
       hypfrqn2d(18,6) = 462604.361
       relintn2d(18,6) = 0.003
       hypfrqn2d(19,6) = 462606.102
       relintn2d(19,6) = 0.003
       hypfrqn2d(20,6) = 462606.122
       relintn2d(20,6) = 0.002
       hypfrqn2d(21,6) = 462606.132
       relintn2d(21,6) = 0.003
       nhypn2d(6) =       21

       hypfrqn2d(1,7) = 539682.029
       relintn2d(1,7) = 0.077
       hypfrqn2d(2,7) = 539682.038
       relintn2d(2,7) = 0.091
       hypfrqn2d(3,7) = 539682.056
       relintn2d(3,7) = 0.109
       hypfrqn2d(4,7) = 539682.074
       relintn2d(4,7) = 0.092
       hypfrqn2d(5,7) = 539682.078
       relintn2d(5,7) = 0.107
       hypfrqn2d(6,7) = 539682.096
       relintn2d(6,7) = 0.123
       hypfrqn2d(7,7) = 539682.098
       relintn2d(7,7) = 0.109
       hypfrqn2d(8,7) = 539682.106
       relintn2d(8,7) = 0.124
       hypfrqn2d(9,7) = 539682.119
       relintn2d(9,7) = 0.141
       nhypn2d(7) =        9


       hypfrqn2d(1,8) = 616749.934
       relintn2d(1,8) = 0.081
       hypfrqn2d(2,8) = 616749.942
       relintn2d(2,8) = 0.094
       hypfrqn2d(3,8) = 616749.957
       relintn2d(3,8) = 0.109
       hypfrqn2d(4,8) = 616749.969
       relintn2d(4,8) = 0.095
       hypfrqn2d(5,8) = 616749.974
       relintn2d(5,8) = 0.108
       hypfrqn2d(6,8) = 616749.989
       relintn2d(6,8) = 0.122
       hypfrqn2d(7,8) = 616749.991
       relintn2d(7,8) = 0.110
       hypfrqn2d(8,8) = 616749.999
       relintn2d(8,8) = 0.123
       hypfrqn2d(9,8) = 616750.010
       relintn2d(9,8) = 0.137
       nhypn2d(8) =        9


c These numbers below are good, but I only made dimensions for
c 8 hyperfines. Should be easy to change the dimensions.
c       hypfrqn2d(1,9) = 693806.018
c       relintn2d(1,9) = 0.085
c       hypfrqn2d(2,9) = 693806.026
c       relintn2d(2,9) = 0.096
c       hypfrqn2d(3,9) = 693806.038
c       relintn2d(3,9) = 0.110
c       hypfrqn2d(4,9) = 693806.047
c       relintn2d(4,9) = 0.097
c       hypfrqn2d(5,9) = 693806.052
c       relintn2d(5,9) = 0.108
c       hypfrqn2d(6,9) = 693806.065
c       relintn2d(6,9) = 0.121
c       hypfrqn2d(7,9) = 693806.066
c       relintn2d(7,9) = 0.110
c       hypfrqn2d(8,9) = 693806.074
c       relintn2d(8,9) = 0.122
c       hypfrqn2d(9,9) = 693806.084
c       relintn2d(9,9) = 0.134
c       nhypn2d(9) =        9
c
c       hypfrqn2d(1,10) = 770848.804
c       relintn2d(1,10) = 0.088
c       hypfrqn2d(2,10) = 770848.812
c       relintn2d(2,10) = 0.098
c       hypfrqn2d(3,10) = 770848.823
c       relintn2d(3,10) = 0.110
c       hypfrqn2d(4,10) = 770848.828
c       relintn2d(4,10) = 0.098
c       hypfrqn2d(5,10) = 770848.835
c       relintn2d(5,10) = 0.109
c       hypfrqn2d(6,10) = 770848.846
c       relintn2d(6,10) = 0.120
c       hypfrqn2d(7,10) = 770848.846
c       relintn2d(7,10) = 0.110
c       hypfrqn2d(8,10) = 770848.854
c       relintn2d(8,10) = 0.121
c       hypfrqn2d(9,10) = 770848.864
c       relintn2d(9,10) = 0.132
c       nhypn2d(10) =        9
c       hypfrqn2d(1,11) = 847876.816
c       relintn2d(1,11) = 0.090
c       hypfrqn2d(2,11) = 847876.824
c       relintn2d(2,11) = 0.099
c       hypfrqn2d(3,11) = 847876.834
c       relintn2d(3,11) = 0.110
c       hypfrqn2d(4,11) = 847876.838
c       relintn2d(4,11) = 0.100
c       hypfrqn2d(5,11) = 847876.844
c       relintn2d(5,11) = 0.109
c       hypfrqn2d(6,11) = 847876.854
c       relintn2d(6,11) = 0.110
c       hypfrqn2d(7,11) = 847876.855
c       relintn2d(7,11) = 0.120
c       hypfrqn2d(8,11) = 847876.862
c       relintn2d(8,11) = 0.120
c       hypfrqn2d(9,11) = 847876.871
c       relintn2d(9,11) = 0.130
c       nhypn2d(11) =        9
c       hypfrqn2d(1,12) = 924888.578
c       relintn2d(1,12) = 0.092
c       hypfrqn2d(2,12) = 924888.586
c       relintn2d(2,12) = 0.101
c       hypfrqn2d(3,12) = 924888.595
c       relintn2d(3,12) = 0.110
c       hypfrqn2d(4,12) = 924888.597
c       relintn2d(4,12) = 0.101
c       hypfrqn2d(5,12) = 924888.604
c       relintn2d(5,12) = 0.110
c       hypfrqn2d(6,12) = 924888.613
c       relintn2d(6,12) = 0.110
c       hypfrqn2d(7,12) = 924888.614
c       relintn2d(7,12) = 0.119
c       hypfrqn2d(8,12) = 924888.620
c       relintn2d(8,12) = 0.119
c       hypfrqn2d(9,12) = 924888.629
c       relintn2d(9,12) = 0.129
c       nhypn2d(12) =        9
c



        i = nlines
        if (i .gt. 8) then
        write(6,*) 'There is data for only 8 N2H+ lines'
        write(6,*) 'Set nlines to 8 or less for this molecule'
        write(6,*) 'Stopping the program'
        if (myid .lt. 2) then
          write(22,*) 'There is data for only 8 N2H+ lines'
          write(22,*) 'Set nlines to 8 or less for this molecule'
          write(22,*) 'Stopping the program'
        endif
        n2dphyp = 0
        return
        endif


        do 77 line = 1,nlines
                nhyp(line) = nhypn2d(line)
 77        continue

        do 78 line = 1,nlines
        do 78 ihyp = 1,nhyp(line)
                relint(ihyp,line) = relintn2d(ihyp,line)
                hypfrq(ihyp,line) = hypfrqn2d(ihyp,line)
 78        continue
        do 79 line = 2,nlines
        do 79 ihyp = 1,nhyp(line)
                hypfrq(ihyp,line) = hypfrq(ihyp,line)/1000.
 79        continue

c Fix up the hyperfines
        i = isqueeze(vwmin,hypfrq)

        return
        end


        function n2hphyp(vwmin)

        include 'nlines_f77.h'
        parameter (maxhyp=50)

        common /procid/ myid
        double precision hypfrq(maxhyp,nlines),frq
        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &      grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        dimension tmp1(maxhyp),tmp2(maxhyp),tmp3(maxhyp)
        dimension hypfrqn2h(22,8),nhypn2h(8),relintn2h(22,8)

c        write(6,*) 'Working out N2H+ hyperfine structure'
        n2hphyp = 1

        hypfrqn2h( 1, 1) =   93.1716086
        relintn2h( 1, 1) = 0.0370
        hypfrqn2h( 2, 1) =   93.1719054
        relintn2h( 2, 1) = 0.1852
        hypfrqn2h( 3, 1) =   93.1720403
        relintn2h( 3, 1) = 0.1111
        hypfrqn2h( 4, 1) =   93.1734675
        relintn2h( 4, 1) = 0.1852
        hypfrqn2h( 5, 1) =   93.1737643
        relintn2h( 5, 1) = 0.2593
        hypfrqn2h( 6, 1) =   93.1739546
        relintn2h( 6, 1) = 0.1111
        hypfrqn2h( 7, 1) =   93.1762527
        relintn2h( 7, 1) = 0.1111
        nhypn2h( 1) =  7
        hypfrqn2h( 1, 2) =  186.3427964
        relintn2h( 1, 2) = 0.0202
        hypfrqn2h( 2, 2) =  186.3429204
        relintn2h( 2, 2) = 0.0400
        hypfrqn2h( 3, 2) =  186.3429618
        relintn2h( 3, 2) = 0.0113
        hypfrqn2h( 4, 2) =  186.3430541
        relintn2h( 4, 2) = 0.0367
        hypfrqn2h( 5, 2) =  186.3432658
        relintn2h( 5, 2) = 0.0634
        hypfrqn2h( 6, 2) =  186.3435179
        relintn2h( 6, 2) = 0.0134
        hypfrqn2h( 7, 2) =  186.3443890
        relintn2h( 7, 2) = 0.0679
        hypfrqn2h( 8, 2) =  186.3444527
        relintn2h( 8, 2) = 0.0131
        hypfrqn2h( 9, 2) =  186.3445240
        relintn2h( 9, 2) = 0.0283
        hypfrqn2h(10, 2) =  186.3447107
        relintn2h(10, 2) = 0.0118
        hypfrqn2h(11, 2) =  186.3447494
        relintn2h(11, 2) = 0.1460
        hypfrqn2h(12, 2) =  186.3447698
        relintn2h(12, 2) = 0.0902
        hypfrqn2h(13, 2) =  186.3447793
        relintn2h(13, 2) = 0.1142
        hypfrqn2h(14, 2) =  186.3448501
        relintn2h(14, 2) = 0.2048
        hypfrqn2h(15, 2) =  186.3451424
        relintn2h(15, 2) = 0.0319
        hypfrqn2h(16, 2) =  186.3452569
        relintn2h(16, 2) = 0.0218
        hypfrqn2h(17, 2) =  186.3472664
        relintn2h(17, 2) = 0.0107
        hypfrqn2h(18, 2) =  186.3474014
        relintn2h(18, 2) = 0.0137
        hypfrqn2h(19, 2) =  186.3474781
        relintn2h(19, 2) = 0.0142
        hypfrqn2h(20, 2) =  186.3476131
        relintn2h(20, 2) = 0.0320
        hypfrqn2h(21, 2) =  186.3476981
        relintn2h(21, 2) = 0.0060
        hypfrqn2h(22, 2) =  186.3477303
        relintn2h(22, 2) = 0.0085
        nhypn2h( 2) = 22
        hypfrqn2h( 1, 3) =  279.5098245
        relintn2h( 1, 3) = 0.0100
        hypfrqn2h( 2, 3) =  279.5098785
        relintn2h( 2, 3) = 0.0169
        hypfrqn2h( 3, 3) =  279.5098982
        relintn2h( 3, 3) = 0.0072
        hypfrqn2h( 4, 3) =  279.5111328
        relintn2h( 4, 3) = 0.0157
        hypfrqn2h( 5, 3) =  279.5113445
        relintn2h( 5, 3) = 0.0488
        hypfrqn2h( 6, 3) =  279.5113847
        relintn2h( 6, 3) = 0.0226
        hypfrqn2h( 7, 3) =  279.5114140
        relintn2h( 7, 3) = 0.0073
        hypfrqn2h( 8, 3) =  279.5114305
        relintn2h( 8, 3) = 0.0143
        hypfrqn2h( 9, 3) =  279.5115089
        relintn2h( 9, 3) = 0.0963
        hypfrqn2h(10, 3) =  279.5116858
        relintn2h(10, 3) = 0.0901
        hypfrqn2h(11, 3) =  279.5117978
        relintn2h(11, 3) = 0.0621
        hypfrqn2h(12, 3) =  279.5118083
        relintn2h(12, 3) = 0.1272
        hypfrqn2h(13, 3) =  279.5118097
        relintn2h(13, 3) = 0.1024
        hypfrqn2h(14, 3) =  279.5118114
        relintn2h(14, 3) = 0.1381
        hypfrqn2h(15, 3) =  279.5118486
        relintn2h(15, 3) = 0.0190
        hypfrqn2h(16, 3) =  279.5118621
        relintn2h(16, 3) = 0.1779
        hypfrqn2h(17, 3) =  279.5121195
        relintn2h(17, 3) = 0.0076
        hypfrqn2h(18, 3) =  279.5123171
        relintn2h(18, 3) = 0.0104
        hypfrqn2h(19, 3) =  279.5142219
        relintn2h(19, 3) = 0.0116
        hypfrqn2h(20, 3) =  279.5143427
        relintn2h(20, 3) = 0.0143
        nhypn2h( 3) = 20
        hypfrqn2h( 1, 4) =  372.6705393
        relintn2h( 1, 4) = 0.0059
        hypfrqn2h( 2, 4) =  372.6705729
        relintn2h( 2, 4) = 0.0092
        hypfrqn2h( 3, 4) =  372.6719768
        relintn2h( 3, 4) = 0.0083
        hypfrqn2h( 4, 4) =  372.6721184
        relintn2h( 4, 4) = 0.0128
        hypfrqn2h( 5, 4) =  372.6723529
        relintn2h( 5, 4) = 0.1192
        hypfrqn2h( 6, 4) =  372.6724208
        relintn2h( 6, 4) = 0.1045
        hypfrqn2h( 7, 4) =  372.6724961
        relintn2h( 7, 4) = 0.0996
        hypfrqn2h( 8, 4) =  372.6725246
        relintn2h( 8, 4) = 0.0766
        hypfrqn2h( 9, 4) =  372.6725564
        relintn2h( 9, 4) = 0.1286
        hypfrqn2h(10, 4) =  372.6725589
        relintn2h(10, 4) = 0.1069
        hypfrqn2h(11, 4) =  372.6725665
        relintn2h(11, 4) = 0.1336
        hypfrqn2h(12, 4) =  372.6725984
        relintn2h(12, 4) = 0.1634
        hypfrqn2h(13, 4) =  372.6728570
        relintn2h(13, 4) = 0.0101
        hypfrqn2h(14, 4) =  372.6730646
        relintn2h(14, 4) = 0.0061
        hypfrqn2h(15, 4) =  372.6748890
        relintn2h(15, 4) = 0.0074
        hypfrqn2h(16, 4) =  372.6749552
        relintn2h(16, 4) = 0.0079
        nhypn2h( 4) = 16
        hypfrqn2h( 1, 5) =  465.8228727
        relintn2h( 1, 5) = 0.0058
        hypfrqn2h( 2, 5) =  465.8243347
        relintn2h( 2, 5) = 0.0052
        hypfrqn2h( 3, 5) =  465.8244229
        relintn2h( 3, 5) = 0.0053
        hypfrqn2h( 4, 5) =  465.8247706
        relintn2h( 4, 5) = 0.0638
        hypfrqn2h( 5, 5) =  465.8247787
        relintn2h( 5, 5) = 0.0831
        hypfrqn2h( 6, 5) =  465.8248172
        relintn2h( 6, 5) = 0.1095
        hypfrqn2h( 7, 5) =  465.8248609
        relintn2h( 7, 5) = 0.1058
        hypfrqn2h( 8, 5) =  465.8248669
        relintn2h( 8, 5) = 0.0862
        hypfrqn2h( 9, 5) =  465.8248982
        relintn2h( 9, 5) = 0.1297
        hypfrqn2h(10, 5) =  465.8248999
        relintn2h(10, 5) = 0.1107
        hypfrqn2h(11, 5) =  465.8249092
        relintn2h(11, 5) = 0.1324
        hypfrqn2h(12, 5) =  465.8249325
        relintn2h(12, 5) = 0.1565
        hypfrqn2h(13, 5) =  465.8252747
        relintn2h(13, 5) = 0.0062
        nhypn2h( 5) = 13
        hypfrqn2h( 1, 6) =  558.9666312
        relintn2h( 1, 6) = 0.0734
        hypfrqn2h( 2, 6) =  558.9666414
        relintn2h( 2, 6) = 0.0903
        hypfrqn2h( 3, 6) =  558.9666677
        relintn2h( 3, 6) = 0.1118
        hypfrqn2h( 4, 6) =  558.9666949
        relintn2h( 4, 6) = 0.0921
        hypfrqn2h( 5, 6) =  558.9666971
        relintn2h( 5, 6) = 0.1091
        hypfrqn2h( 6, 6) =  558.9667240
        relintn2h( 6, 6) = 0.2416
        hypfrqn2h( 7, 6) =  558.9667342
        relintn2h( 7, 6) = 0.1307
        hypfrqn2h( 8, 6) =  558.9667527
        relintn2h( 8, 6) = 0.1509
        nhypn2h( 6) =  8
        hypfrqn2h( 1, 7) =  652.0958533
        relintn2h( 1, 7) = 0.0791
        hypfrqn2h( 2, 7) =  652.0958640
        relintn2h( 2, 7) = 0.0938
        hypfrqn2h( 3, 7) =  652.0958837
        relintn2h( 3, 7) = 0.1117
        hypfrqn2h( 4, 7) =  652.0959001
        relintn2h( 4, 7) = 0.0949
        hypfrqn2h( 5, 7) =  652.0959058
        relintn2h( 5, 7) = 0.1096
        hypfrqn2h( 6, 7) =  652.0959264
        relintn2h( 6, 7) = 0.2387
        hypfrqn2h( 7, 7) =  652.0959362
        relintn2h( 7, 7) = 0.1276
        hypfrqn2h( 8, 7) =  652.0959520
        relintn2h( 8, 7) = 0.1447
        nhypn2h( 7) =  8
        hypfrqn2h( 1, 8) =  745.2103402
        relintn2h( 1, 8) = 0.0832
        hypfrqn2h( 2, 8) =  745.2103509
        relintn2h( 2, 8) = 0.0962
        hypfrqn2h( 3, 8) =  745.2103676
        relintn2h( 3, 8) = 0.1116
        hypfrqn2h( 4, 8) =  745.2103772
        relintn2h( 4, 8) = 0.0969
        hypfrqn2h( 5, 8) =  745.2103845
        relintn2h( 5, 8) = 0.1100
        hypfrqn2h( 6, 8) =  745.2104007
        relintn2h( 6, 8) = 0.1118
        hypfrqn2h( 7, 8) =  745.2104019
        relintn2h( 7, 8) = 0.1247
        hypfrqn2h( 8, 8) =  745.2104107
        relintn2h( 8, 8) = 0.1254
        hypfrqn2h( 9, 8) =  745.2104247
        relintn2h( 9, 8) = 0.1402
        nhypn2h( 8) =  9

        i = nlines
        if (i .gt. 8) then
          write(6,*) 'There is data for only 8 N2H+ lines'
          write(6,*) 'Set nlines to 8 or less for this molecule'
          write(6,*) 'Stopping the program'
        if (myid .lt. 2) then
            write(22,*) 'There is data for only 8 N2H+ lines'
            write(22,*) 'Set nlines to 8 or less for this molecule'
            write(22,*) 'Stopping the program'
          endif
          n2hphyp = 0
          return
        endif


        do 77 line = 1,nlines
                nhyp(line) = nhypn2h(line)
 77        continue

        do 78 line = 1,nlines
        do 78 ihyp = 1,nhyp(line)
                relint(ihyp,line) = relintn2h(ihyp,line)
                hypfrq(ihyp,line) = hypfrqn2h(ihyp,line)
 78        continue

c Fix up the hyperfines
        i = isqueeze(vwmin,hypfrq)

        return
        end

c23456789112345678921234567893123456789412345678951234567896123456789712
        function isqueeze(vwmin,hypfrq)

        include 'nlines_f77.h'
        parameter (maxhyp=50)

        double precision hypfrq(maxhyp,nlines),frq
        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &      grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        common /procid/ myid
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /lines/ indexu(nlines),indexl(nlines)
        dimension tmp1(maxhyp),tmp2(maxhyp),tmp3(maxhyp)

c        write(6,*) 'Squeeze the hyperfines'
        isqueeze = 1


c For each line
        do 1 line=1,nlines

c Sum the relative intensities.
c Pick out the strongest hyperfine
        ihypmax = 0
        relintmax = 0.
        sum = 0.

        do 2 ihyp = 1,nhyp(line)

        sum = sum + relint(ihyp,line)

        if (relint(ihyp,line) .gt. relintmax) then
                ihypmax = ihyp
                relintmax = relint(ihyp,line)
        endif

 2      continue

c Normalize the relative intensities to unity and convert
c the frequencies to offsets from the main line. The main
c line is defined as in the JPL catalog even if there is
c no hyperfine line at this frequency. Convert to velocity.

c       if (myid .lt. 2) write(22,*) line,indexu(line),indexl(line),
c     &         freq(indexu(line),indexl(line))
        frq = freq(indexu(line),indexl(line))

        do 3 ihyp = 1,nhyp(line)
          hyptmp = hypfrq(ihyp,line)
          relint(ihyp,line) = relint(ihyp,line)/sum
          hypfrq(ihyp,line) = hypfrq(ihyp,line)*1.d9 - frq
          hypvel(ihyp,line) = -hypfrq(ihyp,line)/frq*c
c	if (myid .lt. 2) write(22,*) 
c     &     ihyp,line,hyptmp,hypfrq(ihyp,line)/1.d9,frq/1.d9,
c     &     hypvel(ihyp,line)/1.e5
 3      continue


        if (myid .lt. 2) then
            write(22,*) 'new line # ',line,
     &        freq(indexu(line),indexl(line))/1.e9,' GHz'
            write(22,*) 'number, hyperfine freq shift in MHz, ',
     &        ' shift in kms, relative intensity'
            do 4 ihyp = 1,nhyp(line)
              write(22,*) ihyp,hypfrq(ihyp,line)/1.e6,
     &          hypvel(ihyp,line)/1.e5,
     &          relint(ihyp,line)
 4          continue
        endif


c Normalize the relative intensities to unity.
        sum = 0.
        do 24 ihyp = 1,nhyp(line)
                sum = sum + relint(ihyp,line)
 24     continue
        do 34 ihyp = 1,nhyp(line)
          relint(ihyp,line) = relint(ihyp,line)/sum
 34     continue

        if (myid .lt. 2) then
          write(22,*) 'relative intensities renormalized'
          do 13 ihyp = 1,nhyp(line)
            write(22,*) ihyp,hypfrq(ihyp,line)/1.e6,
     &          hypvel(ihyp,line)/1.e5,
     &          relint(ihyp,line)
 13       continue
        endif

 1      continue

        if (myid .lt. 2) then
          do 26 line = 1,nlines
            write(22,*) 'new frequency ',line,
     &        freq(indexu(line),indexl(line))/1.e9
c              write(6,*) 'new frequency ',line,
c     &        freq(indexu(line),indexl(line))/1.e9
 26     continue
        endif

        return
        end



        integer function rays(nrays,maxgrid,sphere,
     &                nboxes,nray1,nrpview,nx,ny,nz,
     &                xcenter,ycenter,zcenter,
     &          xplane,yplane,zplane,
     &                nangles,vlong,vlat,
     &                x1,y1,z1,x2,y2,z2,
     &          dxdu,dydu,dzdu,raybox)

        logical printgeom,compute

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        common /procid/ myid
        common /orient/ along,alat,coslat,coslong,
     &    sinlat,sinlong,long,lat
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav


        dimension nx(nboxes),ny(nboxes),nz(nboxes)
        dimension vlat(nangles),vlong(nangles)
        dimension nray1(nangles),nrpview(nangles)

        dimension 
     &                xplane(maxgrid+1,nboxes),
     &                yplane(maxgrid+1,nboxes),
     &                zplane(maxgrid+1,nboxes)
        dimension 
     &                xcenter(maxgrid,nboxes),
     &                ycenter(maxgrid,nboxes),
     &                zcenter(maxgrid,nboxes)
        dimension 
     &                x1(nrays),
     &                y1(nrays),
     &                z1(nrays),
     &                x2(nrays),
     &                y2(nrays),
     &                z2(nrays)

        dimension
     &                dxdu(nrays),
     &                dydu(nrays),
     &                dzdu(nrays)

        integer raybox(nrays)

        rays = 1

        printgeom = .true.
        printgeom = .false.
        if (myid .ge. 2) printgeom = .false.


        if (printgeom) then
        write(22,*) 'function rays'
        write(22,*) 'nrays ',nrays
        write(22,*) 'maxgrid ',maxgrid
        write(22,*) 'sphere ',sphere
        write(22,*) 'nboxes ',nboxes
        write(22,*) 'nx ',nx
        write(22,*) 'ny ',ny
        write(22,*) 'nz ',nz

        do ibox = 1,nboxes
            write(22,*) 'planes for box #',ibox
            write(22,*) 'X planes'
            do i = 1,nx(ibox) + 1
                write(22,*) i,xplane(i,ibox)
            enddo
            write(22,*) 'Y planes'
            do j = 1,ny(ibox) + 1
                write(22,*) j,yplane(j,ibox)
            enddo
            write(22,*) 'Z planes'
            do k = 1,nz(ibox) + 1
                write(22,*) k,zplane(k,ibox)
            enddo
        enddo

        do ibox = 1,nboxes
            write(22,*) 'centers for box #',ibox
            write(22,*) 'X centers'
            do i = 1,nx(ibox) 
                write(22,*) i,xcenter(i,ibox)
            enddo
            write(22,*) 'Y centers'
            do j = 1,ny(ibox) 
                write(22,*) j,ycenter(j,ibox)
            enddo
            write(22,*) 'Z centers'
            do k = 1,nz(ibox) 
                write(22,*) k,zcenter(k,ibox)
            enddo
        enddo


        write(22,*) 'angles '
        do ka = 1,nangles
            write(22,*) ka,vlong(ka),vlat(ka)
        enddo

        endif 
c End of if (printgeom) then  ...

        rays = nrays
        maxrays = nrays

        if (myid .lt. 2) write(22,*) 'maxrays = ',maxrays

        nrays = 0
        iiray = 0

c Here is the loop over the viewing angles
        do 93 ka = 1,nangles

c Copy the vector variable for the latitude and longitude
c to scalar variables. There is probably no point to this, but I
c didn't want to retype and find out. ALAT and ALONG are stored
c in COMMON /ORIENT/ and passed through ORIENT to subroutine
c ROTATE

        alat = vlat(ka)
        along = vlong(ka)

        if (printgeom .and. myid .lt. 2) then
          write( 6,1297) ka,long,lat,along,alat,
     &        pi/180.*along,pi/180.*alat
          write(22,1297) ka,long,lat,along,alat, 
     &       pi/180.*along,pi/180.*alat
 1297        format('new angle long,lat ',3i5,4f10.3)
        endif

c Calculate the sines and cosines for these angles

        coslat = cos(pi/180.*alat)
        sinlat = sin(pi/180.*alat)
        coslong = cos(pi/180.*along)
        sinlong = sin(pi/180.*along)

c Notes on the orientation of the axes and the angles are in 
c subroutine ROTATE
c Done with the angular grid 
c ------------------------------------------------------------

c Move on to the image grid in the plane of the sky

c For the current angular direction raster scan the 2 D viewing grid.
c The vertical dimension is NZ and we use NY for the horizontal
c dimension. Remember that we move the viewing grid 
c to all orientations around the model. The important point
c is that the viewing grid be large enough to cover the model.
c There is very little computational penalty paid for grid points
c whose rays do not intersect the model since these are identified
c and skipped.

c Make the viewing grid the same height as Z and same width
c as the Y axis of the model.
c When the viewing grid is parallel to the Y,Z planes of the
c the model there will be one ray through each voxel of each 
c of the boxes of grids of the model. When we loop through the
c voxels of the bigger boxes, we will identify and skip over
c those rays which also pass through the smaller boxes. 
c Those rays will be computed later while looping over the 
c voxels of the smaller box.

        do 33 iibox = 1,nboxes
        do 33 jgrid = 1,ny(iibox)
        do 33 kgrid = 1,nz(iibox)       

        if (printgeom) write(22,*) 'new grid point ',iibox,jgrid,kgrid
        if (printgeom) write(6,*) 'new grid point ',iibox,jgrid,kgrid

c Check if the v,w coordinates of the ray are within the Y,Z coordinates of
c a smaller box. If so, don't
c compute this ray in the current box. Do it when we get to looping
c through the voxels of the smaller box.

        v = ycenter(jgrid,iibox)
        w = zcenter(kgrid,iibox)

        if (printgeom) write(22,*) 'coordinates at this grid point',v,w

        compute = .true.
        if (iibox .lt. nboxes ) then
          if (    v .ge. yplane(1,            iibox+1) 
     &      .and. v .le. yplane(ny(iibox+1)+1,iibox+1)
     &      .and. w .ge. zplane(1,            iibox+1)
     &      .and. w .le. zplane(nz(iibox+1)+1,iibox+1) ) 
     &                compute = .false.

          if(printgeom)write(22,*)'v,w,ypln_1,ypln_nx,zpln_1,zpln_nx',
     &              v,w,yplane(1,iibox+1),yplane(ny(iibox+1)+1,iibox+1),
     &          zplane(1,iibox+1),zplane(nz(iibox+1)+1,iibox+1)
          if (printgeom) write(22,*) 'j,k,ibox,compute = ',jgrid,kgrid,
     &                iibox,compute
        endif

        if (.not. compute) go to 33

c Done with the preliminaries related to the image grid
c -------------------------------------------------------------

c First we define the starting and ending points of
c each ray. These will be x1,y1,z1 and x2,y2,z2.

c In the plane of the viewing grid, the V,W coordinates
c of the rays located at the voxel centers are,   

c Check if the ray intersects the sphere

        if ( (v*v + w*w) .ge. (sphere*sphere)) compute = .false.

        if (printgeom) then
          if ((ka .eq. 1 ) .and. (myid .lt. 2))
     &       write(22,1076) v,w,v*v + w*w,sphere*sphere
 1076        format ('v,w,vv+ww,sphere**2 ',4f10.6)
        endif

        if (.not. compute) then
          if (ka .eq. 1 
     &        .and. printgeom .and. myid .lt. 2) 
     &             write(22,1972)  iibox,jgrid,kgrid, v,w
 1972        format('At grid point bjk',3i5,2f10.5,'   skip')
        else
          if (ka .eq. 1 ) then 
             iiray = iiray + 1
             if (printgeom .and. myid .lt. 2) 
     &           write(22,1973)  iibox,jgrid,kgrid, v,w,iiray
          endif
 1973        format('At grid point bjk',3i5,2f10.5,' ray # ',i5)
          u1 = -sqrt(sphere*sphere - v*v - w*w)
          u2 = -u1
        endif

c u1 and u2 are usually large, about the same as the radius of the 
c sphere except on the edges, even for the smaller boxes. This is
c because the rays run through the whole model sphere, so the starting
c and ending points are at the edge of the sphere. v and w on the other
c hand are small in the smaller boxes. These specify the location of
c the ray away from the center in the rotated coordinates.

c If we are at a grid point whose ray does not intersect the
c SPHERE, then there is no radiation
c along this ray. In this case COMPUTE 
c will be false. Go on to next angle or grid point. 


        if (.not. compute) go to 33

        nrays = nrays + 1
        if (nrays .gt. maxrays) then
           if (myid .lt. 2) 
     &       write(22,*) 'Need to increase maximum number of rays'
           rays = 0
           return
        endif

        raybox(nrays) = iibox

        dxdu(nrays) = coslat*coslong
        dydu(nrays) = coslat*sinlong
        dzdu(nrays) = -sinlat

        if (printgeom) then
c             write(6,1092) nrays,dxdu(nrays),dydu(nrays),dzdu(nrays)
              write(22,1092) nrays,dxdu(nrays),dydu(nrays),dzdu(nrays)
1092        format('dxdu dydu dzdu',i6,3f10.4)
        endif

c Get the coordinates of the beginning and end of the ray
c into x,y,z

        call rotate(u1,v,w,coslong,sinlong,coslat,sinlat,
     &                x1(nrays),y1(nrays),z1(nrays))
        call rotate(u2,v,w,coslong,sinlong,coslat,sinlat,
     &                x2(nrays),y2(nrays),z2(nrays))

        if (printgeom .and. myid .lt. 2) then
              write(22,*) 'starting coords uvw ',u1,v,w
              write(22,*) 'ending   coords uvw ',u2,v,w
              write(22,*) 'starting coords xyz ',
     &                x1(nrays),y1(nrays),z1(nrays)
              write(22,*) 'ending   coords xyz ',
     &                x2(nrays),y2(nrays),z2(nrays)
        endif

c The check section below has caused more problems than it's worth.
c skip over it

       goto 33

c As a check, find the minimum distance on the line to the origin and
c see if this is inside the model grid. The equations here are
c y = dm*x + bb and z = dn*x + cc. Minimize r^2 = x^2 + y^2 + z^2
c w/res/to x after substituting for y and z.
                dm = (y2(nrays)-y1(nrays))/(x2(nrays)-x1(nrays))
                bb = y1(nrays) - x1(nrays)*dm
                dn = (z2(nrays)-z1(nrays))/(x2(nrays)-x1(nrays))
                cc = z1(nrays) - x1(nrays)*dn
                xmn = -(dm*bb + dn*cc)/(1. + dm**2 + dn**2)
                ymn = dm*xmn + bb
                zmn = dn*xmn + cc


c                write(6,*) 'xmn,ymn,zmn ',xmn,ymn,zmn 
c The check is done on the largest box. The inner boxes must be inside the
c largest box. It is possible for a ray to rotate around so that it's x,y, or z
c is larger than the box x,y,z when projected onto the ray's orientation.
c JJBOX is always one, never changes. 

c Checking the largest box is OK if the setup puts the model sphere inside
c the grid, but not otherwise. Still not right.

c The +1 beyond nx is thre because there is one more plane than there 
c are cells.

            jjbox = 1
            iok = 1
          if (    xmn .lt. xplane(1,          jjbox)
     &      .or. xmn .gt. xplane(nx(jjbox)+1,jjbox) )then
                 iok = 0
                 write(22,*) 'Error in ray location X dimension'
                 write(22,*) 'box nx(box) ',iibox, nx(jjbox)
                 write(22,*) 'xmn,xplane(1),xplane(nx) ',
     &               xmn,xplane(1,jjbox),xplane(nx(jjbox)+1,jjbox)
          endif
          if (    ymn .lt. yplane(1,          jjbox)
     &      .or. ymn .gt. yplane(ny(jjbox)+1,jjbox) )then
                 iok = 0
                 write(22,*) 'Error in ray location Y dimension'
                 write(22,*) 'box ny(box) ',iibox, ny(jjbox)
                 write(22,*) 'ymn,yplane(1),yplane(ny) ',
     &               ymn,yplane(1,jjbox),yplane(ny(jjbox)+1,jjbox)
          endif
          if (    zmn .lt. zplane(1,          jjbox)
     &      .or. zmn .gt. zplane(nz(jjbox)+1,jjbox) )then
                 iok = 0
                 write(22,*) 'Error in ray location Z dimension'
                 write(22,*) 'box nz(box) ',iibox, nz(jjbox)
                 write(22,*) 'zmn,zplane(1),zplane(nz) ',
     &               zmn,yplane(1,jjbox),zplane(nz(jjbox)+1,jjbox)
          endif
          if (iok .eq. 0) then
            write(22,*) 'Problem with ray rotation at ray ',nrays
            write(22,*) 'u1,u2,v,w ',u1,u2,v,w
            write(22,*) 'x1,y1,z1 ',x1(nrays),y1(nrays),z1(nrays)
            write(22,*) 'x2,y2,z2 ',x2(nrays),y2(nrays),z2(nrays)
            write(22,*) 'xmn,ymn,zmn ',xmn,ymn,zmn
c            write(25,1025) iok,nrays,u1,u2,v,w,
c     &      x1(nrays),y1(nrays),z1(nrays),x2(nrays),y2(nrays),z2(nrays)
            nrays = nrays - 1
            goto 33
          endif
 1025     format(i3,i8,10f8.4)

c Output file fort.25 above is read by an IDL program used for 
c debugging the ray tracing, plot3d.pro. If you do not intend
c to run this idl program, there is no need for fort.25


c It does not matter if the sphere is larger than the model grid.
c The subroutine arclengths sets up the pathlengths through the
c grid, and starts them at the grid boundary, not the sphere
c boundary.

 33        continue

            if (ka .eq. 1) then
                nray1(ka) = 1
                nrpview(ka) = nrays
            else
                nray1(ka) = nrpview(ka-1) + nray1(ka-1)
                nrpview(ka) = nrays - nray1(ka) + 1
            endif


           if (myid .lt. 2) then
             write(22,*) 'for angle ',ka,' first ray and ',
     &          'number of rays ',nray1(ka),nrpview(ka)
             write( 6,*) 'for angle ',ka,' first ray and ',
     &          'number of rays ',nray1(ka),nrpview(ka)
           endif


           if(printgeom) write(6,*) 
     &           'finished loop 93 for angle number ',ka
 93        continue
c End of the loop over long,lat angles

c ------------------------------------------------------------

        rays = 1
        return
        end

        subroutine rotate(u,v,w,coslong,sinlong,coslat,sinlat,x,y,z)


c u with x, v with y, and w with z.

c The angle called latitude is the angle measured off
c the vertical axis Z between Z and W. In rotating
c the coords, this angle is applied first.

c The angle called longitude is measured off the X
c axis. This angle is applied second.

c With these definitions, the axis W corresponds to
c the radial coordinate of the usual spherical coordinates.

c We will let the model exist in X,Y,Z coordinates with
c X and Y defining the plane of the equator and Z the
c vertical axis. We agree that the viewing angle
c is along the U axis with the radiation propagating from
c -U to + U (observer at +U). The image plane is the V,W
c plane. When the angles are zero 
c we are looking at the model along the X axis and the
c plane of the sky is defined by Y,Z.

c The 3D rotation to go from U,V,W to X,Y,Z is found as
c the matrix multiplication of the 2 2D rotations.

c l = longitude measured off x axis
c c = latitude measured off vertical axis

c First matrix:


c       x      cos(c)    sin(c)   u
c         = 
c       z     -sin(c)    cos(c)   w


c Second matrix

c       x     cos(l)    -sin(l)   u
c         = 
c       y     sin(l)     cos(l)   v


c       second applied to first

c       x     cos(l)   -sin(l)   0    cos(c)  0   sin(c)   u
c       y =   sin(l)    cos(l)   0    0       1   0        v
c       z     0         0        1   -sin(c)  0   cos(c)   w

        x =  u*coslong*coslat - v*sinlong        + w*coslong*sinlat
        y =  u*sinlong*coslat + v*coslong        + w*sinlong*sinlat
        z = -u*sinlat                            + w*coslat

        return
        end


        integer function raypath(
     &                x1,y1,z1,
     &                x2,y2,z2,
     &                nboxes,nx,ny,nz,voxsizex,voxsizeyz,
     &                maxgrid,xplane,yplane,zplane,
     &                maxvox,iarc,pathlength,icount)

        parameter(maxboxes=12)

        logical perpx,perpy,perpz
        logical startinside,endinside,intersects
        logical printgeom
        logical path0

        dimension 
     &          nx(nboxes),ny(nboxes),nz(nboxes),
     &          voxsizex(nboxes),voxsizeyz(nboxes)
        dimension
     &          xplane(maxgrid+1,nboxes),
     &          yplane(maxgrid+1,nboxes),
     &          zplane(maxgrid+1,nboxes)
        dimension
     &          arcl(maxvox),
     &          pathlength(maxvox),
     &          iarc(4,maxvox)

        dimension
     &      arclmin(maxboxes),arclmax(maxboxes),intersects(maxboxes)

        raypath = 1
        printgeom = .true. 
        printgeom = .false.
        if (myid .ge. 2) printgeom = .false.

        pathmin = 1.e-6

        if (printgeom) write(22,*) 'x1,y1,z1',x1,y1,z1
        if (printgeom) write(22,*) 'x2,y2,z2',x2,y2,z2

c For each ray calculate the intersections of the ray with 
c all the planes defined in all the boxes that the ray passes
c through. ICOUNT is set to zero, then increments with each
c box (that has intersections).

        icount = 0
        do 476 ibox = 1,nboxes

          call arclengths(
     &            x1,y1,z1,x2,y2,z2,
     &            xplane(1,ibox),yplane(1,ibox),zplane(1,ibox),
     &            nx(ibox),ny(ibox),nz(ibox),
     &            voxsizex(ibox),voxsizeyz(ibox),
     &            arclmin(ibox),arclmax(ibox),intersects(ibox),
     &            arcl,icount,maxvox,
     &            printgeom,istatus)

        if (istatus .ne. 1) then 
            write(22,*) 'Bad result from subroutine arclengths '
            write(6,*) 'Bad result from subroutine arclengths '
            raypath = 0
            return
        endif

        if (printgeom) write(22,*) 'intersects ?',intersects(ibox)
        if (printgeom) write(22,*) 'ibox,icount ',ibox,icount


 476        continue

        if (icount .eq. 0) then
            write(6,*) 'Ray should intersect the model, but did not'
            write(22,*) 'Ray should intersect the model, but did not'
            raypath = 0
            return
        endif

c Sort the intersections into ascending order.
        if (printgeom) then
            do 789 i = 1,icount 
               write(22,*) 'i,arcl(i) ',i,arcl(i)
 789            continue
        endif

        call sort(icount,arcl)

        if (printgeom) then
            write(22,*) 'sorted'
            do 779 i = 1,icount 
               write(22,*) 'i,arcl(i) ',i,arcl(i)
 779            continue
        endif


c Now that they are in ascending order, differences of 
c sequential pairs of these arc length variables are
c the lengths of the line segments through the voxels
c pierced by the ray. Multiply by the total raylength
c to get back the units.

        raylength = sqrt( 
     &                  (x2-x1)**2 
     &                + (y2-y1)**2 
     &                + (z2-z1)**2 )

        do 76 i = 1,icount-1
                pathlength(i) = raylength*(arcl(i+1) - arcl(i))
                if (printgeom) 
     &               write(22,*) 'i,arcl,path',i,arcl(i),pathlength(i)
 76        continue
                if (printgeom) write(22,*) 'i,arcl     ',i,arcl(i)

c Now we are done with computing the path lengths. 

c Because we have multiple boxes, the planes of some boxes 
c might be identical. These will look like pathlengths of
c zero length. Throw these out. And we can throw out any tiny
c pathlengths.

        path0 = .false.
        newcount = icount
        do 152 ivox = 1,icount-1
 301        continue
          if (pathlength(ivox) .lt. pathmin) then
            if (ivox .ge. newcount) go to 300
            path0 = .true.
            newcount = newcount - 1
            if (printgeom) then
                write(22,*) 'Zero pathlength segment'
                write(22,*) 'ivox,path',ivox,pathlength(ivox)
                write(22,*) 'new count',newcount
            endif
            do 156 jj = ivox,newcount-1
                arcl(jj) = arcl(jj+1)
                pathlength(jj) = pathlength(jj+1)
                if (printgeom) write(22,*) 'jj,arcl(jj),path(jj)',
     &                        jj,arcl(jj),pathlength(jj)
 156            continue
            arcl(newcount) = arcl(newcount+1)
                if (printgeom) write(22,*) 'jj,arcl(jj)         ',
     &                        jj,arcl(jj)
            go to 301
          endif
 152        continue
 300        continue

        icount = newcount

        if (printgeom) then
        if (path0) then
        write(22,*) 'Number of voxels along ray ',icount-1
        write(22,*) 'ivox   pathlength  arcl'
        do 732 ivox = 1,icount-1
            write(22,2009) ivox,pathlength(ivox),arcl(ivox)
 2009                format(i5,2f12.6)
 732    continue
            write(22,2089) ivox,arcl(ivox)
 2089                format(i5,'            ',2f12.6)
        endif
        endif

c The tiny segments have been eliminated.

c Now identify the voxel that each path segment is in. Use
c the midpoint of the path segment.

        delx = x2 - x1
        dely = y2 - y1
        delz = z2 - z1

        do 77 i = 1,icount-1

c The midpoint of the segment is

        arcmid = 0.5*(arcl(i+1) + arcl(i))

c Identify the smallest box containing the midpoint. Use 
c .gt. and .lt. because if there is a midpoint on cell
c boundary, the next algorithm will pick the outside cell
c which would be beyond the box boundary.

        ibox = 1
        do 272 ib = 1,nboxes
            if (intersects(ib) .and. arcmid .gt. arclmin(ib) 
     &        .and. arcmid .lt. arclmax(ib)) 
     &                ibox = ib
        if (printgeom) write(22,1559) 
     &                i,arcmid,arclmin(ib),arclmax(ib),ib,ibox
 1559   format('ivox,arcmid,arclmin(ib),arclmax(ib),ib,ibox',i5,
     &    3f10.5,2i5)
 272        continue

        if (printgeom) write(22,*) 'The working box number for ivox'
     &                ,i,' is ',ibox

c Now that we know the box, find the voxel indices.

        iarc(1,i) = (x1 + arcmid*(delx) - xplane(1,ibox))
     &      /voxsizex(ibox)
        iarc(1,i) = iarc(1,i) + 1
        iarc(2,i) = (y1 + arcmid*(dely) - yplane(1,ibox))
     &      /voxsizeyz(ibox)
        iarc(2,i) = iarc(2,i) + 1
        iarc(3,i) = (z1 + arcmid*(delz) - zplane(1,ibox))
     &      /voxsizeyz(ibox)
        iarc(3,i) = iarc(3,i) + 1

        iarc(4,i) = ibox

        if (printgeom) write(22,1078) i,iarc(1,i),iarc(2,i),iarc(3,i),
     &                iarc(4,i)
 1078   format('i,iarc(1,i),iarc(2,i),iarc(3,i),iarc(4,i)  ',5(i6,2x))


c No indices beyond the box boundaries allowed.
        if (    
     &          iarc(1,i) .lt. 1 .or.
     &                iarc(2,i) .lt. 1 .or.
     &                iarc(3,i) .lt. 1 .or.
     &          iarc(4,i) .lt. 1 .or.
     &                iarc(1,i) .gt. nx(ibox) .or.
     &                iarc(2,i) .gt. ny(ibox) .or.
     &                iarc(3,i) .gt. nz(ibox) .or.
     &                iarc(4,i) .gt. nboxes
     &      )then
                write(6,*) 'Bad index for some voxel, cannot continue '
                if (myid .lt. 2) then
                write(22,*) 'Bad index for some voxel, cannot continue '
                write(22,1078) i,iarc(1,i),iarc(2,i),iarc(3,i),iarc(4,i)
                write(22,2543)  ibox,nboxes,nx(ibox),ny(ibox),nz(ibox)
 2543   format ('ibox,nboxes,nx(ibox),ny(ibox),nz(ibox)  ',5(i6,2x))
                write(22,*) 'x1,y1,z1',x1,y1,z1
                write(22,*) 'x2,y2,z2',x2,y2,z2
                write(22,*) 'i arcl(i+1) arcl(i) ',i,arcl(i+1),arcl(i)
                write(22,*) 'arcmid ',arcmid
                write(22,*) 'delx,dely,delz',delx,dely,delz
                write(22,*)'coordinates of midpoint',x1 + arcmid*(delx),
     &                        y1 + arcmid*(dely),z1 + arcmid*(delz)
                write(22,*) 'xplane,yplane,zplane',
     &                      xplane(1,ibox),yplane(1,ibox),yplane(1,ibox)
                write(22,*) 'voxsize',voxsizex(ibox),voxsizeyz(ibox)
                write(22,*) 'Check for pathlengths of zero'
        do 734 ivox = 1,icount-1
            write(22,2009) ivox,pathlength(ivox),arcl(ivox)
 734    continue
                write(6,*) 'Finished diagnostics, ',
     &              'return with stop signal'
                write(22,*) 'Finished diagnostics, ',
     &              'return with stop signal'
                endif
                raypath = 0
                return
        endif

c Done identifying all the voxels along the ray.
 77        continue

        if (printgeom) then
        write(22,*) 'Here are the voxel indices along the ray'
        do 471 i = 1,icount-1
                write(22,7832) i,iarc(1,i),iarc(2,i),iarc(3,i),iarc(4,i)
 7832        format('i,iarc(1,i),iarc(2,i),iarc(3,i),iarc(4,i)',5i5)
 471        continue
        endif

c Now we have the path of the radiation through the pixels
c computed for this ray

        return
        end


        subroutine arclengths(x1,y1,z1,x2,y2,z2,
     &             xplane,yplane,zplane,nx,ny,nz,voxsizex,voxsizeyz,
     &             arclmin,arclmax,intersects,arcl,icount,
     &             maxvox,
     &             printgeom,istatus)


        logical perpx,perpy,perpz,intersects
        logical startinside,endinside
        logical printgeom

        dimension xplane(nx),yplane(ny),
     &                  zplane(nz)
        dimension arcl(maxvox)

c Inputs

c        x1,y1,z1,x2,y2,z2
c        nx,ny,nz
c        xplane,yplane,zplane

c Outputs

c        arclmin,arclmax arc length parameters for intersection of
c                        ray with the boundaries of the cube
c        intersects        false if ray does not intersect the cube
c        icount                number of intersections
c        arcl(i)                intersections of ray with grid planes of 
c                        this cube




c Find where the ray intersects all the planes
c of the model grid.

c The comments below refer to variables along the X axis, but
c these comments should be thought of as referring to a 3D
c vector X or a list. For example, where x1 is mentioned,
c think of y1, z1 as well. 

c Define arc length variable along ray, normalized so that
c ray has a length of 1 between x1 and x2. 

c If the ray is perpendicular to the x axis, then x2-x1 = 0, and
c the ray never intersects the planes of that axis.

c        write(22,*) 'printgeom ',printgeom 

        istatus = 1

        lun = 22

        delx = x2 - x1
        dely = y2 - y1
        delz = z2 - z1

        perpx = .false.
        perpy = .false.
        perpz = .false.
        if (abs(delx) .lt. 1.e-12) perpx = .true.
        if (abs(dely) .lt. 1.e-12) perpy = .true.
        if (abs(delz) .lt. 1.e-12) perpz = .true.

        if (printgeom) then
        write(lun,*) 'maxvox ',maxvox,icount
        write(lun,*) 'delx,dely,delz',delx,dely,delz
        write(lun,*) 'perpx,perpy,perpz',perpx,perpy,perpz
        endif

c If the ray is perpendicular
c to X, then the starting and ending points of the ray, X1 and X2,
c are the same. If the point X1 is less than XPLANE(1) or greater
c than XPLANE(NX), then the ray is outside and parallel to one 
c face of the model grid.

        if (perpx) then
          if (x1 .lt. xplane(1) .or. x1 .gt. xplane(nx+1)) then
                intersects = .false.
              if (printgeom) write(lun,*) 'x1,xpl(1),xpl(nx+1) ',
     &            x1,xplane(1),xplane(nx+1),intersects
              return
          endif
        endif
        if (perpy) then
          if (y1 .lt. yplane(1) .or. y1 .gt. yplane(ny+1)) then
               intersects = .false.
               if (printgeom) write(lun,*) 'y1,ypl(1),ypl(ny+1) ',
     &             y1,yplane(1),yplane(ny+1),intersects
               return
          endif
        endif
        if (perpz) then
          if (z1 .lt. zplane(1) .or. z1 .gt. zplane(nz+1)) then
                intersects = .false.
              if (printgeom) write(lun,*) 'z1,zpl(1),zpl(nz+1) ',
     &            z1,zplane(1),zplane(nz+1),intersects
              return
          endif
        endif
                
c FInd the arc lengths at the intersections of the ray with the 
c edges of the model grid. 

        if (.not. perpx) then
              arclx1 = (xplane(1) - x1)/delx
              arclnx = (xplane(nx+1) - x1)/delx
              if (printgeom) write(lun,*) 'xp1 xpnx x1 delx',
     &        xplane(1),xplane(nx+1),x1,delx
              if (printgeom) write(lun,*) 'arclx1,arclnx',arclx1,arclnx
        endif
        if (.not. perpy) then
              arcly1 = (yplane(1) - y1)/dely
              arclny = (yplane(ny+1) - y1)/dely
              if (printgeom) write(lun,*) 'arcly1,arclny',arcly1,arclny
        endif
        if (.not. perpz) then
              arclz1 = (zplane(1) - z1)/delz
              arclnz = (zplane(nz+1) - z1)/delz
              if (printgeom) write(lun,*) 'arclz1,arclnz',arclz1,arclnz
        endif


c Find the first point on the ray where it enters the model 
c grid and the last point where it leaves.

        arclmin = -1.0e30
        if (.not. perpx) arclmin = max( arclmin, min(arclx1,arclnx) ) 
        if (.not. perpy) arclmin = max( arclmin, min(arcly1,arclny) )
        if (.not. perpz) arclmin = max( arclmin, min(arclz1,arclnz) )

c Check if the start point is inside the model cube
        startinside = .false.

        if (arclmin .lt. 0.) then
                arclmin = 0.
                startinside = .true.
        endif

        arclmax = +1.0e30
        if (.not. perpx) arclmax = min( arclmax, max(arclx1,arclnx) )
        if (.not. perpy) arclmax = min( arclmax, max(arcly1,arclny) )
        if (.not. perpz) arclmax = min( arclmax, max(arclz1,arclnz) )

c Check if the end point is inside the model cube
        endinside = .false.
        if (arclmax .gt. 1.) then
                arclmax = 1.
                endinside = .true.
        endif
        if (printgeom) then
             write(lun,*)'arclmin,startinside ',arclmin,startinside 
             write(lun,*)'arclmax,endinside ',arclmax,endinside
             write(lun,*) 'arclmin,arclmax ',arclmin,arclmax
        endif



c If ray does not intersect the model, go on to next grid point.

        intersects = .true.
        if (arclmax .le. arclmin) then
                intersects = .false.
                return
        endif

c Find the minimum and maximum indices of all the planes crossed
c by the ray.

          aa = arclmin
          ab = arclmax
          if (delx .lt. 0.) then
            aa = arclmax
            ab = arclmin
          endif

        if (.not. perpx) then
          i1 =  (xplane(nx+1) - (aa*delx + x1) )/voxsizex
          iarcmin = nx - i1
          i2 = -(xplane(1)    - (ab*delx + x1) )/voxsizex
          iarcmax = i2
          if (printgeom) then
              write(lun,1010) i1,i2,nx,xplane(1),xplane(nx+1)
              write(lun,1011) aa,ab,delx,x1,voxsize,iarcmin,iarcmax
 1010     format('i1,i2,nx,xp1,xpnx,aa,ab,delx,vsz ',3i5,1p,2e13.6)
 1011     format(1p,5e13.6,2x,2i5)
          endif
        endif

          aa = arclmin
          ab = arclmax
          if (dely .lt. 0.) then
            aa = arclmax
            ab = arclmin
          endif

        if (.not. perpy) then
          j =  (yplane(ny+1) - (aa*dely + y1) )/voxsizeyz
          jarcmin = ny - j
          j = -(yplane(1)    - (ab*dely + y1) )/voxsizeyz
          jarcmax = j
        endif

          aa = arclmin
          ab = arclmax
          if (delz .lt. 0.) then
            aa = arclmax
            ab = arclmin
          endif

        if (.not. perpz) then
          k =  (zplane(nz+1) - (aa*delz + z1) )/voxsizeyz
          karcmin = nz - k
          k = -(zplane(1)    - (ab*delz + z1) )/voxsizeyz
          karcmax = k
        endif

        if (printgeom .and. .not. perpx) 
     &     write(lun,*) 'iarcmin,iarcmax',iarcmin,iarcmax
        if (printgeom .and. .not. perpy) 
     &     write(lun,*) 'jarcmin,jarcmax',jarcmin,jarcmax
        if (printgeom .and. .not. perpz) 
     &     write(lun,*) 'karcmin,karcmax',karcmin,karcmax

c Now find the arclengths corresponding to all the
c intesections with the planes. 

        if (.not. perpx) then
          if (iarcmax .ge. iarcmin) then
            do 73 i = iarcmin,iarcmax
              arc = (xplane(i) - x1)/delx
              if (arc .ge. arclmin .and. arc .le. arclmax) then
                  icount = icount + 1
                  arcl(icount) = arc
                  if (printgeom) write(lun,1001) 
     &              icount,xplane(i),x1,delx,arcl(icount)
 1001   format('ic,xp,x1,dx,arcl ',i5,1p,4e13.6)
              endif
 73         continue
          endif
        endif

        if (.not. perpy) then
          if (jarcmax .ge. jarcmin) then
            do 74 j = jarcmin,jarcmax
              arc = (yplane(j) - y1)/dely
              if(arc .ge. arclmin .and. arc .le. arclmax) then
                  icount = icount + 1
                  arcl(icount) = arc
                  if (printgeom) write(lun,*) 'ic,yp,y1,dy,arcl ',
     &              icount,yplane(j),y1,dely,arcl(icount)
              endif
 74         continue
          endif
        endif

        if (.not. perpz) then
          if (karcmax .ge. karcmin) then
            do 75 k = karcmin,karcmax
              arc = (zplane(k) - z1)/delz
              if(arc .ge. arclmin .and. arc .le. arclmax) then
                  icount = icount + 1
                  arcl(icount) = arc
                  if (printgeom) write(lun,*) 'ic,zp,z1,dz,arcl ',
     &              icount,zplane(k),z1,delz,arcl(icount)
              endif
 75         continue
          endif
        endif

        if (printgeom) then
        do 675 i = 1,icount
                write(lun,*) 'sub arclengths: i,arcl',i,arcl(i)
 675        continue
        endif

        do 676 i = 1,icount
           if (arcl(i) .ge. 0.) goto 676
              write(lun,*) 'Bad arc length ',i,arcl(i)
              write(lun,*) 'nx,ny,nz ',nx,ny,nz 
              write(lun,*) 'perpx,perpy,perpz ',perpx,perpy,perpz 
              if (.not. perpx) write(lun,*) 'ic,xp,x1,dx,arcl ',
     &              icount,xplane(i),x1,delx,arcl(icount)
              if (.not. perpy) write(lun,*) 'ic,yp,y1,dy,arcl ',
     &              icount,yplane(j),y1,dely,arcl(icount)
              if (.not. perpz) write(lun,*) 'ic,zp,z1,dz,arcl ',
     &              icount,zplane(k),z1,delz,arcl(icount)
              istatus = 0
              return
 676    continue

c If the start is inside the model grid, then add in the
c arc length parameter for the start point. This is arcl = 0.
c If the end is inside, add in arcl = 1.


        if (startinside) then
                icount = icount + 1
                arcl(icount) = 0.
        endif

        if (endinside) then
                icount = icount + 1
                arcl(icount) = 1.
        endif

        if (icount .gt. maxvox) then 
            write(22,*) 'More segments than maximum allowed ',
     &          icount,maxvox
            write(6,*) 'More segments than maximum allowed ',
     &          icount,maxvox
            istatus = 0
            return
        endif


c Now we have all the arc length variables corresponding to
c the intersections with all of the x,y, and z planes, and
c if the ray is inside the cube, the starting and ending 
c points of the ray.

        return
        end

      SUBROUTINE SORT(N,RA)
      DIMENSION RA(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END

        integer function setcmb(initialradiation,yycmb)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        include 'nlines_f77.h'
        common /cmb/ ycmb(nlines)
        common /procid/ myid
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav

        dimension yycmb(nlines)

        setcmb = 1

        if (myid .lt. 2) then
        if (initialradiation .eq. 0)
     &   write(22,*) 'initial radiation field is ZERO  ',
     &      initialradiation
        if (initialradiation .ge. 1)
     &    write(22,*) 'initial radiation field is the CMB  ',
     &      initialradiation
        endif

        if (initialradiation .eq. 0) then
           do 85 line = 1,nlines
              ycmb(line) = 1.e-16
              yycmb(line) = 1.e-16
 85        continue

        else 

c        write(6,*) 'planck,boltz ',planck,boltz

c The YCMB variable loads the F77 COMMON block CMB
c YYCMB goes back to the C program through the function args
        do 84 line = 1,nlines
             ycmb(line) = planck*freqline(line)/boltz
     &         / (dexp(planck*freqline(line)/(boltz*2.7d0))-1.d0)
             yycmb(line) = ycmb(line)
 84     continue

        endif

        if (myid .lt. 2) then
c          write(6 ,*) 'Here is the CMB initial line brightness'
c          write(6 ,*) 'at each transition frequency '
          write(22,*) 'Here is the CMB initial line brightness'
          write(22,*) 'at each transition frequency '
          write(22,*)  'Line  Freq      CMB (K)'
          do 996 line = 1,nlines
c             write(6 ,*) 'Freq CMB (K)',freqline(line)/1.e9,
c     &           ycmb(line)
 1222        format(i5,f12.3,f10.6)
             write(22,1222) line,freqline(line)/1.e9,
     &           ycmb(line)
 996        continue
        endif

        return
        end

        function initrad(nboxes,nx,ny,nz,initial_radiation)

        EXTERNAL putalbj

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        integer u,l

        include 'nlines_f77.h'
        parameter (maxhyp=50)
c        integer putalbj


        common /procid/ myid
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /cmb/ ycmb(nlines)
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)

        dimension nx(nboxes),ny(nboxes),nz(nboxes)

c This is the return value for this function
        initrad = 1


        if (myid .lt. 2) write(22,*) 'nboxes ',nboxes
        if (myid .lt. 2) write(22,*) 'nx,ny,nz ',nx(1),ny(1),nz(1)

c Check whether the CMB is set to ZERO or not and report. The CMB
c was set in SETCMB.
        if (myid .lt. 2) then
         if (initial_radiation .eq. 1) then 
           write(22,*)'Using CMB approx for initial local radiation'
           write( 6,*)'Using CMB approx for initial local radiation'
         endif
         if (initial_radiation .eq. 0) then 
           write(22,*)'Using ZERO approx for initial local radiation'
           write( 6,*)'Using ZERO approx for initial local radiation'
         endif
         if (initial_radiation .eq. 2) then 
           write(22,*)'Using LTE  approx for initial local radiation'
           write( 6,*)'Using LTE  approx for initial local radiation'
         endif
        endif

c Before starting the loop, initialize the grid of the average 
c radiation J bar. This can be the CMB or zero 

c        if (myid .lt. 2) write(22,*) 'Initializing J bar'
c        if (myid .lt. 2) write(6,*) 'Initializing J bar'
 
        do 991 ibox = 1,nboxes
        do 991 i = 1,nx(ibox)
        do 991 j = 1,ny(ibox)
        do 991 k = 1,nz(ibox)

c Decrement indices for the C function PUTSRCOPC
        ib = ibox-1
        ii = i - 1
        jj = j - 1
        kk = k - 1

c Loop over the lines, decrement the line index for the 
c C function PUTSRCOPC
        do 992 line = 1,nlines
          lline = line - 1

         if (initial_radiation .eq. 0) bj = 0.

c This would set the mean radiation field to the planck function. 
c This has never been a good idea in any of the applications so far.

        if (initial_radiation .eq. 2) then
            call getmodel(ib,ii,jj,kk,temperature,density,vx,vy,vz,
     &                abund,widthl,iadds)
            bj = planck*freqline(line)/boltz
     &        / ( dexp(planck*freqline(line)
     &          / boltz/dble(temperature))
     &             -1.d0 )
        endif



c Initialize the approximate lambda operator to zero
c Initialize the normalization of the line profile to 1.0
c Dont forget to subtract 1 from ihyp to go from F77 to C

        al = 0.0
        xn = 1.0
        if (initial_radiation .eq. 1) bj = ycmb(line)

        nh = 1
        if (nltehyp .eq. 1) nh = nhyp(line)
        do 500 ihyp = 1,nh
          iihyp = ihyp - 1
c        write(22,1000) ib,ii,jj,kk,lline,iihyp,al,bj
c 1000   format('send to putalbj ',6i4,1p2e12.4)
          call putalbj(ib,ii,jj,kk,lline,iihyp,al,bj,xn)
 500    continue


 992        continue
 991        continue

c        if (myid .lt. 2) write(22,*) 'Finished initializing J bar'
c        if (myid .lt. 2) write(6,*) 'Finished initializing J bar'
 

        return
        end

c23456789112345678921234567893123456789412345678951234567896123456789712
        integer function acceleratedlambda(
     &    tk,h2,
     &    xnm,bJ,aL,emis,opac,
     &    frac,
     &    iaccel,igrid,jgrid,kgrid,iter
     &    )

c Computes the fractional population given a mean radiation
c density bJ and lambda acceleration parameter aL.

c In this code M is the upper state and N is the lower state.
c The Einstein A values A(N,M) are set equal to A(M,N) for
c convenience. Similarly FREQ(N,M) = FREQ(M,N) and 
c barJmn(N,M) = barJmn(M,N). If M=N all these, A,FREQ,BARJMN
c are zero.

c Downward collision rates are written C(U,L): upper to lower
c Upward collision rates are written C(L,U): lower to upper
c Collision rates to same level, U=L, are set to zero.

c The Einstein A*c^2/(2hv^3) is used where the 
c Einstein B is needed.

c The total downward rates, collision+radiative, are written 
c as P(L <-- U), P(L,U). 
c The upward rates are written P(U <-- L), P(U,L).
c This is the reverse of the left-to-right order of the collision
c rates. The indices are written in this reversed order so that
c the matrices of rates are ordered correctly as for a linear
c system of equations.

c INDEXU(LINE), INDEXL(LINE) are the Fortran array indices that
c point to the upper and lower states of the line. For rotors,
c with rotational transitions, for example N2H+ described
c with LTE hyperfines, these are the quantum numbers 
c of the states+1. We usually use 5 to 9 rotational
c transition. For non-LTE hyperfine N2H+, these point
c to one of the 64 states which is defined as the 
c "main" hyperfine line.
c INDEXU and INDEXL are F77 style, running from 1 to N.

c UPPER and LOWER are indices that point to the upper and lower
c states of all the hyperfines of a molecule defined with
c non-lte hyperfine lines. For example,
c There are 280 hyperfine lines for N2H+
c UPPER and LOWER are defined in C style, c running from 0 to N-1. 

c For molecules with non-LTE hyperfines, INDEXU and INDEXL 
c still refer to the rotational lines, and point to the 
c the states that are defined as the "main" hyperfine
c transition. They are then equal to UPPER+1 and LOWER+1
c for that transition.

c APL1D and BARJ1D are 1D vectors from the C program that
c contain the accleration factor, lambda, and the mean
c radiation field. Unpack the 1D vectors into 2D arrays indexed 
c by MAXHYP,NLINES. 

c Because the output FRAC is a 1D vector, there is no 2D->1D packing
c of the output at the end of this subroutine.

        include 'nlines_f77.h'
        parameter (maxch=2000)
        parameter (maxhyp=50)
        parameter (maxvox=8192)
        parameter (maxtrans=500)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline,
     &  frac_up,frac_lo

        character *80 msg
        double precision p,d,ct,db,sumpops
        double precision barJmn,barJeffmn,alambda,xnorm
        integer upper,lower,dludcmp

        logical doprint

c The double precision variable d is returned by LUDCMP
c but is probably not used.

        dimension bJ(maxhyp,nlines)
        dimension barJmn(nstate,nstate)
        dimension xnm(maxhyp,nlines)
        dimension aL(maxhyp,nlines),source(maxhyp,nlines)
        dimension emis(maxhyp,nlines),opac(maxhyp,nlines)
        dimension barJeffmn(nstate,nstate)
        dimension alambda(nstate,nstate),xnorm(nstate,nstate)
        dimension frac_up(nlines),frac_lo(nlines),frac(nstate)

        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /gridcom/ igridcom,jgridcom,kgridcom
        common /p/ p(nstate,nstate)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /coll/ cr(nstate,nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /lines/ indexu(nlines),indexl(nlines)
        common /nltehyp/ upper(maxtrans),lower(maxtrans),
     &          istartline(nlines)
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
c        common /linelevs/ jloline(nlines),jupline(nlines)
        common /cmb/ ycmb(nlines)

        dimension indx(nstate)

c Debugging variables
c        double precision bt(nstate),at(nstate,nstate),berror
c        double precision ap(nstate,nstate),vp(nstate,nstate),wp(nstate)
c        double precision wpmin,wpmax

c        write(22,*) 'ijk iter ',igrid,jgrid,kgrid,iter
c        write(22,*) 'tk h2 ',tk,h2
c        write(22,*) 'iaccel ',iaccel
c        if (iter .eq. 1) write(22,*) 'first iter al=em=op=0.0'
c        if (iter .eq. 2) write(22,*) '2nd iter'

c        do 901 line = 1,nlines
c        nh = 1
c        if (nltehyp .eq. 1) nh = nhyp(line)
c        write(22,*) 'nltehyp nh ',nltehyp,nh
c        do 901 ihyp = 1,nh
c            write(22,1901) 
c     &         line,ihyp,bJ(ihyp,line),aL(ihyp,line),
c     &         emis(ihyp,line),opac(ihyp,line)
c 1901       format('input l h bJ al em op ',2i5,1p,4e12.4)
c 901    continue

c Return code to signal success of this function. Will be set
c to zero if an error is detected
        acceleratedlambda = 1

c Turns on debug printing for j,k = 16 and i < 16. Good for 32^3 grid
        doprint = .false.
        if (igrid .le. 16 .and. jgrid .eq. 16 .and. kgrid .eq. 16 
     &        .and. iter .le. 100) doprint = .true.
        if (myid .ge. 2) doprint = .false.

c This turns off the debug printing
        doprint = .false.

        igridcom = igrid
        jgridcom = jgrid
        kgridcom = kgrid

        if (doprint .and. myid .lt. 2) then
          write(22,1322) igridcom,jgridcom,kgridcom,iter
          write(22,*) 'doprint is ',doprint,' iaccel = ',iaccel
 1322     format(' accel lambda: i,j,k,iter ',4i5)
        endif

c        write(22,2003) igrid,jgrid,kgrid,xnm(1,1)
c 2003   format('i,j,k xn',3i5,1p,2x,e12.6)

c This provides a return value for the function 
c 1 = good
c 0 = failed
        acceleratedlambda = 1

c Calculate the source function as the ratio of the emissivity and opacity
c The units are degrees K. 
        do 902 line = 1,nlines
        nh = 1
        if (nltehyp .eq. 1) nh = nhyp(line)
        do 902 ihyp = 1,nh
          if (opac(ihyp,line) .ne. 0.0) then
          source(ihyp,line) = emis(ihyp,line)
     &                      / opac(ihyp,line)
     &        * ((c/freqline(line))**2)/(2.d0*boltz)
          else 
            source(ihyp,line) = 0.0
          endif
        if (doprint) write(22,2019) 
     &      ihyp,line,opac(ihyp,line),emis(ihyp,line),source(ihyp,line)
 2019     format('ih l opac emis source ',2i5,1p,3e16.8)
 902    continue

c load up the matrix P with radiative and collision rates
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c First compute the matrix of total transition rates. 

        flag = 0.
        do 2 m = 1,nstate
        do 2 n = 1,nstate
                p(m,n)=dble(flag)
 2      continue

        istatus = 1
c Set the collision rates for this temperature
        if (molecule .eq. 6 .or. molecule .eq. 22) then
                call colrateh2o(tk,doprint)
 76     continue
        endif

c Collision rates for molecules CO=7, C13O=10, C17O=13, C18O=14
        if (     molecule .eq. 7 
     &      .or. molecule .eq. 10 
     &      .or. molecule .eq. 13
     &      .or. molecule .eq. 14 
     &     ) then
                call colratco(tk)
        endif

c Collision rates for rotational levels of HCN, ground vibrational
        if (molecule .eq. 15) then
                call colrate_hcn_ground(tk)
        endif

c Collision rates for CS = 8 and SiO = 20
        if (molecule .eq. 8 .or. molecule .eq. 20) then
                call colratcs(tk)
        endif

c Collision rates for N2Hplus = 9 and N2Dplus = 16
        if (molecule .eq. 9 .or. molecule .eq. 16) then
                call colrateion(tk)
        endif

c Collision rates for H13COplus = 11 and HCOplus = 12
        if (molecule .eq. 11 .or. molecule .eq. 12) then
c                call colratehco(tk,istatus)
                call colrateion(tk)
        endif

c The subroutine loadionrates should have been run in F77 function
c setup to use colrateion. This is done automatically, and this
c comment is a reminder that colrateion depends on loadionrates

c Collision rates for N2Hhyp = 4 (N2H+ with non-LTE hyperfines)
c     colrate_nhyp uses the rotational rates from the LAMBDA data base
c     and approximates the rates between hyperfine levels using the
c     statistical weights. The subroutine n2h_deltaJ_rates should have
c     been run in F77 function setup to use colrate_nhyp

c If you edit colrateion, you can use the ion rates for N2H+ with 
c non-LTE hyperfines. Also edit F77 setup to run loadionrates

        if (molecule .eq. 4 ) then
                iprint = 0
c                if (igrid .lt. 4 .and. jgrid .eq. 0 
c     &              .and. kgrid .eq. 0) iprint = 1
                call colrate_nhyp(tk,iprint)
        endif



c New colrathcn (HCN has used colratco before), 4 Sep 08, R.Rolffs
        if (     molecule .eq. 19
     &      .or. molecule .eq. 17
     &      .or. molecule .eq. 18) call colrathcn(tk)



c Use colrate_hcnhyp2 with loadhcnrates
       if (molecule .eq. 21 ) then
                iprint = 0
c                if (igrid .lt. 4 .and. jgrid .eq. 0 
c     &              .and. kgrid .eq. 0) iprint = 1
                call colrate_hcnhyp2(tk,iprint)
        endif

c        write(*,*) 'Finished colrate_hcnhyp'

        if (istatus .eq. 0) then
              acceleratedlambda = 0
              return
        endif
c        if(doprint)write(6,*) 'collision rates loaded'

        if (doprint) then
          lun = 22
          msg = ' matrix of collision rates'
          call dmattmt(cr,msg,nstate,lun)
        endif

c Now set the radiative rates.

c The next loop (100) handles the case with
c fewer lines than we have transitions between molecular states.
c For the transitions that are not computed,
c we can either set the radiation to the CMB or to
c ZERO. 
c The if (a(m,n .ne. 0.) means
c that we leave barJmn at zero if there is no line between
c states m and n.

c Rather than figure out the Delta J line for each of these
c transitions, determine whether we are using the CMB or
c ZERO approximations as set in SETCMB. Then recompute
c the CMB using FREQ(M,N) or set the radiation to zero.

        cmbyes = 1.
        if (ycmb(1) .lt. 1.e-6) cmbyes = 0.

        if (doprint) then
          if (ycmb(1) .lt. 1.e-6) then
          write(22,*) 'Initialize barJ to ZERO'
        else
          write(22,*) 'Initialize barJ to CMB '
        endif
        endif

        do 100 m = 1,nstate
        do 100 n = 1,nstate

c This is ZERO. This is also the initialization for
c states with no radiative transitions (A(M,N) = 0), so this 
c needs to be first.
          barJmn(m,n) = 0.d0
c This is LTE. Not using this.
c          if (a(m,n) .ne. 0. ) 
c     &      barJmn(m,n) = boltz/planck/freq(m,n) * tk
c This is CMB in units of K.
          if (a(m,n) .ne. 0. ) barJmn(m,n) 
     &      = dble(cmbyes) * planck*freq(m,n)/boltz
     &        / (dexp(planck*freq(m,n)/(boltz*2.7d0))-1.d0)

c 1078     format('Init A barJmn ',2i5,1p2e12.4)
c          if (doprint) write(22,1078) 
c     &        m,n,a(m,n),barJmn(m,n)

c Initial Lambda acceleration is off
                alambda(m,n) = 0.d0
                barJeffmn(m,n) = barJmn(m,n)
 100        continue

c For those lines where we have calculated the radiation, transfer
c the radiation from indexing by lines, bJ(nlines), to indexing by 
c state, barJmn(nstate,nstate). There are 2 different cases depending
c on whether the molecule has hyperfines and how these are calculated.
c
c 1. No hyperfines. NHYP = 1 for all levels.
c
c  or LTE hyperfines. In this case we are not computing the populations
c    in each hyperfine state. The brightness refers to rotational
c    line.
c
c  In this case, copy aL,barJ, and xn with a 1 in the first index.
c  This is equivalent to setting ihyp = 1 and aL(ihyp,line), for example
c
c 2. NLTE hyperfines. Distribute bJ to the all states.
c
c  In this case, loop over the hyperfines, do ihyp = 1,nhyp(line)
c 
c In all cases put the radiation into the right units in barJmn.
c
c Do the same for the accelerated lamda parameter
c
c The barJ equations below have frequencies in the denominator
c This is OK because the loop is over NLINES, and freq and freqline
c are never zero for a line transition. 

c -----------------------------------------------------------------
c This THEN block of the IF THEN ELSE ENDIF is for molecules with 
c no hyperfine structure or LTE hyperfine structure. 

        if (nltehyp .eq. 0) then

c -----------------------------------------------------------------

          do 101 line = 1,nlines
c INDEXU, INDEXL are defined c F77 style, 1 to N.
            m = indexu(line)
            n = indexl(line)
            if (freqline(line) .lt. 10.) then
              write (22,*) 'BAD FREQUENCY A1'
              write ( 6,*) 'BAD FREQUENCY A1'
            endif

            alambda(m,n) = 0.d0
            if (iaccel .eq. 1) alambda(m,n) = dble(aL(1,line))

            barJmn(m,n) = 
     &           dble(bJ(1,line))
            barJmn(n,m) = barJmn(m,n)

              xnorm(m,n) = dble(xnm(1,line))
              xnorm(n,m) = xnorm(m,n)

c            write(22,2001) m,n,igrid,jgrid,kgrid,xnm(1,1),
c     &        xnorm(m,n)
c 2001       format('m,n,i,j,k,xnorm',5i5,3x,1p,2e12.6)

c barJeff is the effective barJ used in ALI
c Normalization should be OK here because barJmn and alambda are
c calculated the same way
              barJeffmn(m,n) = barJmn(m,n) 
     &           - alambda(m,n) * dble(source(1,line))

        if (doprint) then
        write(22,*) 'in loop m,n ',m,n
        write(22,*) 'barJmn barJeffmn ',barJmn(m,n),barJeffmn(m,n)
        write(22,*) 'alambda source   ',alambda(m,n),source(1,line)
        endif

            if (barJeffmn(m,n) .lt. 0) then
              alambda(m,n) = barJmn(m,n) / dble(source(1,line))
              write(22,1200) barJeffmn(m,n),m,n,i,j,k,alambda(m,n)
 1200         format('Case 1 Negative Jb',1pe12.4,' mnijk ',5i5,e12.4)
              barJeffmn(m,n) = 0.
            endif

            alambda(n,m) = alambda(m,n)
            barJeffmn(n,m) = barJeffmn(m,n)
 101      continue

c -----------------------------------------------------------------

c This ELSE block of the IF THEN ELSE ENDIF is for molecules with 
c hyperfine structure
        else 

c -----------------------------------------------------------------

            do 102 line = 1,nlines
c Here we assign each bJ(ihyp,line) to a separate hyperfine
c transition
            do 102 ihyp = 1,nhyp(line)
c This indexing works because istartline begins at 0 (C or QN style)
c  and ihyp begins at 1 (F77 style). 
c UPPER, LOWER are defined C style, 0 to N-1 so need to add 1.
              k = istartline(line) + ihyp
              m = upper(k) + 1
              n = lower(k) + 1

c Check the indexing
              if ( (n .le. 0 .or. n .gt. nstate)
     &        .or.  (m .le. 0 .or. m .gt. nstate) ) then
                  write(6,*) 'BAD INDEXING ',m,n,nstate
                  write(22,*) 'BAD INDEXING ',m,n,nstate
              endif
            if (freqline(line) .lt. 10. .or. freq(m,n) .lt. 10.) then
              write (22,*) 'BAD FREQUENCY A2'
     &          ,line,k,ihyp,istartline(line),m,n,
     &           freqline(line),freq(m,n)
              write ( 6,*) 'BAD FREQUENCY A2'
     &          ,line,m,n,freqline(line),freq(m,n)
            endif

              alambda(m,n) = 0.d0
              if (iaccel .eq. 1) alambda(m,n) = dble(aL(ihyp,line))

              barJmn(m,n) = 
     &            dble(bJ(ihyp,line))
              barJmn(n,m) = barJmn(m,n)

              xnorm(m,n) = dble(xnm(ihyp,line))
              xnorm(n,m) = xnorm(m,n)

c barJeff is the effective barJ used in ALI
c Normalization should be OK here because barJmn and alambda are
c calculated the same way
              barJeffmn(m,n) = barJmn(m,n)
     &         - alambda(m,n) * dble(source(ihyp,line))

        if (doprint) then
        write(22,*) 'in loop m,n ',m,n
        write(22,*) 'barJmn barJeffmn ',barJmn(m,n),barJeffmn(m,n)
        write(22,*) 'alambda source   ',alambda(m,n),source(ihyp,line)
        endif

          if (barJeffmn(m,n) .lt. 0) then
              alambda(m,n) = barJmn(m,n) / dble(source(ihyp,line))
              write(22,1201) barJeffmn(m,n),m,n,i,j,k,alambda(m,n)
 1201         format('Case 2 Negative Jb',1pe12.4,' mnijk ',5i5,e12.4)
              barJeffmn(m,n) = 0.
          endif
              alambda(n,m) = alambda(m,n)
              barJeffmn(n,m) = barJeffmn(m,n)
 102        continue

c -----------------------------------------------------------------

c This is the end of the IF THEN ELSE ENDIF for non-LTE hyperfines
        endif

        if (doprint) then
c          write(6,*) 'finished re-indexing Jbar'
          lun = 22
          msg = ' matrix of barJ'
          call dmattmt(barJmn,msg,nstate,lun)
        endif


c        if (doprint) write(6,*) 'finished re-indexing alambda'


c        if (doprint) write(6,*) 'finished computing Jbar effective'
        if (doprint) then
c          write(6,*) 'finished re-indexing Jbar'
          lun = 22
          msg = ' matrix of barJeffmn'
          call dmattmt(barJeffmn,msg,nstate,lun)
        endif

c more debugging
        if (doprint) then
        write(22,*) 'barJmn barJeffmn(1,2) ',barJmn(1,2),barJeffmn(1,2)
        write(22,*) 'alambda source        ',alambda(1,2),source(1,1)
        endif

c Now compute the transition rates as the sum of collisional and 
c radiative rates.

c        if (doprint)write(6,*) 'start to compute transition rates'
        if (doprint) then
            write(22,1022) ,iter
        endif
 1022   format ('iter=',i3,'    i -> f',t26,'Freq.',t36,
     &    'Einst. A',t46,'X norm',t56,'JbarEff',t70,
     &    'gi',t73,'gf',t77,'Col rate',t86,'H2 cm-3',t96,
     &    'total rate')


        do 1 m = 1,nstate
        do 1 n = 1,nstate

        if (m.eq.n) go to 1

c Store the rates in the format P(n <-- m). This
c will be standard matrix notation where the first index increases
c downward (y axis) and the second increase rightward (x axis).
c This is the transpose of "image notation" where the first index
c increases along the x axis and the second on y.

c The Einstein B's are the A multiplied by c^2/(2hv^3)
c assuming Jbar were in brightness units. Since Jbar is in units of K,
c need a factor 2k*(v/c)^2 to restore Jbar to brightness. Combined
c with the c^2/(2hv^3) factor to convert A to B results in a
c factor k/hv. 

c Downward rates. Einstein A, Einstein B (stimulated emission), 
c and collisions. 
c Here we need to use the numerical normalization (N - lambda) instead
c of (1 - alambda).

c Some of the N2H+ hyperfine states have the same energy, and FREQ 
c which is the energy difference in GHz is zero. FREQ should be
c at least MHz, 1E6. Check so as not to divide by zero.
        db = 0.d0
        if (freq(m,n) .gt. 1.) db = boltz / (planck*freq(m,n))
c        cr(2,1) = 0.d0
c        cr(1,2) = 0.d0

c             if(m .gt. n) p(n,m) = a(m,n)*(1.d0 - alambda(m,n))
             if(m .gt. n) then
                 p(n,m) = a(m,n)*(xnorm(m,n) - alambda(m,n))
     &  + db*a(m,n)*barJeffmn(m,n) + cr(m,n)*dble(h2)
c                 p(n,m) = cr(m,n)*dble(h2)   for LTE only
             endif

c Upward rates. Einstein B (absorption) and collisions

              if (m .lt. n) then
                  p(n,m) =
     &  statdg(n)/statdg(m)*db*a(m,n)*barJeffmn(m,n) + cr(m,n)*dble(h2)
c                  p(n,m) = cr(m,n)*dble(h2)   for LTE only
              endif

        if (doprint) then
              write(22,1000) m,n,
     &              jindx(m),kindx(m),lindx(m),
     &              jindx(n),kindx(n),lindx(n),freq(m,n),
     &              a(m,n),xnorm(m,n),
     &              barJeffmn(m,n)*planck*freq(m,n)/boltz,
     &              int(statdg(m)),int(statdg(n)),
     &              cr(m,n),
     &              h2,
     &              p(n,m)
        endif
 1000   format (i2,i3,2x,3i2,' ->',3i2,2x,1p,3e10.3,e13.6,0p,
     &          2i3,1x,1p,3e10.3)

 1      continue

c save the original rates for debugging, 2 levels only
        p21 = p(2,1)
        p12 = p(1,2)

        do 3 k = 1,nstate
                p(k,k) = 0.
                do 4 n = 1,nstate
                        if(k .eq. n) go to 4
                        p(k,k) = -p(n,k) + p(k,k)
 4              continue
 3      continue

c Drop the last row of the rate matrix P and replace
c it with the normalization: the sum of all level
c populations is unity. When calling ludcmp/lubksb, the
c vector FRAC contains b in Ax=b in the input to these
c routines and on the return contains the solution x.


c The 2 lines  
c        p(nstate,i) = 1.d0
c        frac(nstate) = 1.d0
c Need to be activated (uncommentated) to use ludcmp and lubksb
c Deactivated to use George Rybicki's rate matrix solver

        do 94 i = 1,nstate
                 frac(i) = 0.d0 
c                p(nstate,i) = 1.d0
 94        continue
c        frac(nstate) = 1.d0

c        if(doprint)write(6,*) 'transitions computed'

c DLUDCMP is an integer function
c        acceleratedlambda = dludcmp(p,nstate,nstate,indx,d)
c        if (acceleratedlambda .eq. 0) return
c        call dlubksb(p,nstate,nstate,indx,frac)

c This is George Rybicki's special rate matrix solver.
c There is a maximum matrix size in RATE
c Currently dimensioned to 200 so nstate must be no more than 200
c Change the line below to active for George's solver.
c Comment out for LUDCMP
        call rate(p,nstate,nstate,1.0d0,frac)

c The fraction in each level relative to a total population of
c 1 is in frac.

        sumpops = 0.d0
        if (doprint) then
          write(22,*) 'iter,ijk=',iter,igrid,jgrid,kgrid
          write(22,*) 'level pops: populations in levels'
          do 88 i = 1,nstate
                  write(22,*) i,frac(i)
                  sumpops = sumpops + frac(i)
 88          continue
c Check the solution. Remember that the indexing is P(1 <- 2),
c second index is the initial state
           balance = frac(2)*p12 - frac(1)*p21
        write(22,*) 'Sum of level pops, balance  ',sumpops,balance
        endif

c debugging: 
        if (doprint) then
c compare (1-0) source function to planck function
           db = 2.d0*planck*freq(2,1) * ((freq(2,1)/c)**2) 
           src = a(2,1)*frac(2) * db
     &         / (frac(1)*a(1,2)*statdg(2)/statdg(1) - frac(2)*a(2,1))
           plf   =  db / (dexp(planck*freq(2,1)/(boltz*dble(tk)))-1.d0)
           if (src .gt. plf) then
               write(22,*) 'Bad source > planck ',src,plf
           else
               write(22,*) 'source  planck ', src,plf
           endif
        endif





c This kills the program gracefully at this point. For debugging.
c        if (igrid .eq. 10 .and. jgrid .eq. 10 .and. kgrid .eq. 10)
c     &  acceleratedlambda = 0

        return
        end

c23456789112345678921234567893123456789412345678951234567896123456789712

c23456789112345678921234567893123456789412345678951234567896123456789712

        subroutine colratehco(tk,istatus)

        include 'nlines_f77.h' 
        parameter (nt=5,nr=5)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        double precision dekt
        logical doprint

        integer u,l,dj
        common /procid/ myid
        dimension rate(nt,nr,nr)
        dimension temp(nt)
        dimension ac2(nr,nr-1),bc2(nr,nr-1)
        dimension ac1(15),bc1(15)

        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)

        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /coll/ cr(nstate,nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole

        data temp / 5., 10., 20., 30., 40. /

        data rate/
     &   0.000 ,  0.000 ,  0.000 ,  0.000 ,  0.000 ,
     &   2.584 ,  2.479 ,  2.247 ,  2.091 ,  2.154 ,
     &   1.622 ,  1.499 ,  1.329 ,  1.219 ,  1.217 ,
     &   0.796 ,  0.832 ,  0.725 ,  0.631 ,  0.633 ,
     &   0.804 ,  0.825 ,  0.741 ,  0.667 ,  0.631 ,
     &   3.295 ,  4.849 ,  5.443 ,  5.438 ,  5.805 ,
     &   0.000 ,  0.000 ,  0.000 ,  0.000 ,  0.000 ,
     &   4.955 ,  4.493 ,  3.983 ,  3.705 ,  3.719 ,
     &   2.732 ,  2.865 ,  2.592 ,  2.330 ,  2.350 ,
     &   1.471 ,  1.591 ,  1.538 ,  1.430 ,  1.374 ,
     &   0.623 ,  2.076 ,  3.499 ,  3.974 ,  4.414 ,
     &   1.492 ,  3.183 ,  4.328 ,  4.643 ,  5.004 ,
     &   0.000 ,  0.000 ,  0.000 ,  0.000 ,  0.000 ,
     &   3.561 ,  4.023 ,  4.007 ,  3.868 ,  4.020 ,
     &   2.780 ,  2.892 ,  2.876 ,  2.757 ,  2.655 ,
     &   0.033 ,  0.447 ,  1.406 ,  1.877 ,  2.332 ,
     &   0.088 ,  0.787 ,  2.075 ,  2.664 ,  3.210 ,
     &   0.383 ,  1.561 ,  2.953 ,  3.530 ,  4.082 ,
     &   0.000 ,  0.000 ,  0.000 ,  0.000 ,  0.000 ,
     &   3.737 ,  3.816 ,  3.816 ,  3.714 ,  3.691 ,
     &   0.001 ,  0.103 ,  0.786 ,  1.441 ,  1.948 ,
     &   0.002 ,  0.101 ,  0.673 ,  1.188 ,  1.574 ,
     &   0.013 ,  0.261 ,  1.158 ,  1.828 ,  2.259 ,
     &   0.157 ,  0.886 ,  2.085 ,  2.699 ,  3.093 ,
     &   0.000 ,  0.000 ,  0.000 ,  0.000 ,  0.000 /

        data ac2  /
     &  0.000 ,   6.978 ,   6.810 ,   7.871 ,   6.244 , 
     &  0.000 ,   0.000 ,   6.831 ,   6.488 ,   5.824 , 
     &  0.000 ,   0.000 ,   0.000 ,   4.900 ,   4.659 , 
     &  0.000 ,   0.000 ,   0.000 ,   0.000 ,   5.639 /

        data bc2  /
     &  0.000 ,  -0.546 ,  -0.623 ,  -1.084 ,  -0.964 , 
     &  0.000 ,   0.000 ,  -0.695 ,  -0.816 ,  -0.875 , 
     &  0.000 ,   0.000 ,   0.000 ,  -0.744 ,  -0.812 , 
     &  0.000 ,   0.000 ,   0.000 ,   0.000 ,  -0.696 /

        data ac1  /
     &   8.926571 ,    7.663315 ,    5.105464 ,    5.639232, 
     &   3.321160 ,    2.998259 ,    2.706752 ,    2.443587,     
     &   2.206008 ,    1.991528 ,    1.797901 ,    1.623100,   
     &   1.465293 ,    1.322829 ,    1.194217  /
        data bc1   /
     &   -1.120110 ,   -0.946826 ,   -0.821651 ,   -0.696373,
     &   -0.658220 ,   -0.534964 ,   -0.434789 ,   -0.353372,
     &   -0.287201 ,   -0.233421 ,   -0.189711 ,   -0.154187,
     &   -0.125314 ,   -0.101849 ,   -0.082777 /

        doprint = .false.
        istatus = 1

        if (doprint)then
        write(6,1004) (temp(i),i=1,nt)
 1004        format (9x,5f8.0)
        do 22 u = 1,nr
        do 22 l = 1,nr
                write(6,1005) u-1,l-1,(rate(i,u,l),i=1,nt)
 1005   format(2i5,7f8.4)
 22     continue
 11     continue

        do 24 i = 1,nstate
                write(6,*) 'level =',i,'     stat deg =',statdg(i)
 24        continue

        do 23 u = 2,nr
        do 23 dj = 1,u-1
                l = u - dj
                write(6,1008) u-1,l-1,dj,ac2(u,dj),bc2(u,dj)
 1008        format(i2,'-->',i2,'   dj=',i2,4x,2f8.4)
 23        continue
        endif

c If the temperature is less than 40, interpolate from the tables for J up to NR
c Must use .lt. 40 rather than .le. 40 
c For J > NR use the formula with coefficients which
c are functions of j and dj. 

        if (tk .lt. 40.) then

        it = tk/10. + 1
c        write(6,*) 'tk it ',tk,it
c        write(6,*) 'T < 40 interpolate rates for lower states'


        do 2 u = 1,nr
        do 2 l = 1,nr

c If the temperature, tk, is less than 5 K, do a linear interpolation between
c zero rate at zero K and the rates in the table at 5 K. Can't extrapolate
c downward using the 5 K and 10 K rates because that results in negative
c rates at low temperature.

        if (tk .lt. 5.) then
                r2 = rate(it,u,l)
                r1 = 0.
                t2 = temp(it)
                t1 = 0.
        else
                r2 = rate(it+1,u,l)
                r1 = rate(it,u,l)
                t2 = temp(it+1)
                t1 = temp(it)
        endif

c This is a simple linear interpolation.
        diff = r2 - r1
        if (diff .eq. 0.) then
                cr(u,l) = dble (r2)
        else
                cr(u,l) = dble (r1
     &          + (tk - t1)/(t2 - t1)
     &          * (r2 - r1)  )
        endif

        if (cr(u,l) .lt. 0.) then
           if (myid .lt. 2) then
                write(6,*) 'negative collision rate ',cr(u,l)
                write(22,*) 'negative collision rate ',cr(u,l)
           endif
           istatus = 0
           return
        endif

c        if (myid .lt. 2) write(22,1003) tk,u-1,l-1,r1,r2,cr(u,l)
c        write(6,1003) tk,u-1,l-1,r1,r2,cr(u,l)
 1003        format(f5.0,i5,'-->',i5,3f8.4)

 2      continue

c Now for temperatures > 40 K, use the 2 parameter coefficients,
c functions of J and delta J.
        else

c        write(6,*) 'T > 40, use formula with 2 parameter coefficients'
        do 1 u = 2,nr
        do 1 dj = 1,u-1

        l = u - dj

c        write(6,*) 'l,u',l,u
c        write(6,*) 'planck,boltz,brot',planck,boltz,brot
c        write(6,*) 'tk',tk  

c The formula for the energy difference is dE = h*b(u(u+1) - l(l+1))
c but remember that in fortran, the array indices are not the quantum
c numbers as they are in C. So here the formula is modified because
c U and L are the quantum numbers + 1.
        dekt = planck/boltz * dble(   brot*((u-1)*u - (l-1)*l)   )
     &      / dble(tk)


c Downward rate
        cr(u,l) = ac2(u,dj)*statdg(l)/statdg(u)
     &        * (1.d0 + dekt)
     &        * exp(bc2(u,dj)*sqrt(dekt))

c Upward rate
        cr(l,u) = cr(u,l)*(statdg(u)/statdg(l))*exp(-dekt)

c        if (myid .lt. 2) write(22,1007) tk,u-1,l-1,cr(u,l),cr(l,u)
 1007        format(f5.0,i5,'-->',i5,2f8.4)
        write(6,1010) u-1,l-1,dj,(u-1)*u - (l-1)*l,
     &                planck/boltz*brot,dekt,cr(u,l),cr(l,u)
 1010        format(4i5,4f10.6)




        if (cr(u,l) .lt. 0.) then
         if (myid .lt. 2) write(22,*) 'negative collision rate ',cr(u,l)
         istatus = 0
         return
        endif

 1      continue

        endif
c This is the end of the if tk < 40 block




c Now do the higher transitions using the one parameter a and b
c coefficients, fucntions of DJ alone.

        do 5 u = nr+1,nstate
        do 5 dj = 1,u-1
        l =  u - dj
c Downward rate
        dekt = planck/boltz * dble(   brot*((u-1)*u - (l-1)*l)   )
     &      / dble(tk)
        cr(u,l) = ac1(dj)*statdg(l)/statdg(u)*(1. + dekt)
     &      * exp(bc1(dj)*sqrt(dekt))
c Upward rate
        cr(l,u) = cr(u,l)*(statdg(u)/statdg(l))*exp(-dekt)

c        write(6,1011) u-1,l-1,dj,cr(u,l),cr(l,u)
 1011        format(3i5,2f8.4)
 5        continue

        do 67 u = 1,nstate
        do 67 l = 1,nstate
                cr(u,l) = cr(u,l) * 1.e-10
 67        continue

        if (doprint) then

        write(6,1002) tk
        if (myid .lt. 2) write(22,1002) tk
 1002   format(t1,'  u',t6,'  l',t12,'rate',t24,'T = ',f5.0)

        do 7 l = 1,nstate-1
        write(6,1001)
        if (myid .lt. 2) write(22,1001)
 1001   format(/)
        do 7 u = l+1,nstate
        write(6,1000) u-1,l-1,cr(u,l)
        if (myid .lt. 2) write(22,1000) u-1,l-1,cr(u,l)
 1000   format(t1,i3,t6,i3,1p,t12,d14.6)
 7      continue
        endif


        return
        end




        subroutine colratcs(tk)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,hvkt,ct,aline,freqline,
     &  guline,glline

        integer u,l

        include 'nlines_f77.h'

        common /procid/ myid
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)

        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /coll/ cr(nstate,nstate)


             dimension acs(20),bcs(20),ccs(20)

c These are the coefficients for Turner etal 1992 ApJ 339 114
        acs(1) = 10.35759
        acs(2) = 11.03258
        acs(3) = 8.883451
        acs(4) = 9.641931
        acs(5) = 7.809433
        acs(6) = 9.199441
        acs(7) = 7.425523
        acs(8) = 8.857458
        acs(9) = 7.884091
        acs(10) = 8.701254
        acs(11) = 9.093239
        acs(12) = 9.102847
        acs(13) = 10.15139
        acs(14) = 10.41153
        acs(15) = 10.98320
        acs(16) = 11.65907
        acs(17) = 12.23172
        acs(18) = 12.56786
        acs(19) = 12.98000
        acs(20) = 13.46871

        bcs(1) = -14.04851
        bcs(2) = -13.87741
        bcs(3) = -11.72670
        bcs(4) = -12.36264
        bcs(5) = -10.37700
        bcs(6) = -11.94096
        bcs(7) = -9.617364
        bcs(8) = -11.83214
        bcs(9) = -10.09489
        bcs(10) = -11.66852
        bcs(11) = -12.07070
        bcs(12) = -12.20918
        bcs(13) = -13.60956
        bcs(14) = -14.00692
        bcs(15) = -14.73210
        bcs(16) = -15.58413
        bcs(17) = -16.36966
        bcs(18) = -16.78471
        bcs(19) = -17.32250
        bcs(20) = -17.99528

        ccs(1) = 5.088159
        ccs(2) = 4.820635
        ccs(3) = 3.472995
        ccs(4) = 3.672912
        ccs(5) = 2.615620
        ccs(6) = 3.154852
        ccs(7) = 2.137866
        ccs(8) = 2.902678
        ccs(9) = 2.189489
        ccs(10) = 2.616295
        ccs(11) = 2.699072
        ccs(12) = 2.661824
        ccs(13) = 3.018327
        ccs(14) = 3.059366
        ccs(15) = 3.227232
        ccs(16) = 3.388685
        ccs(17) = 3.575172
        ccs(18) = 3.667707
        ccs(19) = 3.788750
        acs(20) = 3.972845


c Set all rates to zero
        do 3 u = 1,nstate
        do 3 l = 1,nstate
                cr(u,l) = 0.d0
 3        continue

c We only have 16 states, so we will not use all 20 coefficients.
c Transitions to same level are zero rate.
c Remember that the index numbers u,l are not the quantum
c numbers. The ground state quantum number 0 corresponds to 
c index number 1.

c        write(6,1002) tk
c        if (myid .lt. 2) write(22,1002) tk
c 1002        format(t1,'  u',t6,'  l',t12,'rate',t24,'T = ',f5.0)


        do 1 l = 1,nstate-1
c        write(6,1001)
c        if (myid .lt. 2) write(22,1001)
c 1001        format(/)
        do 1 u = l+1,nstate

        hvkt = planck*freq(u,l)/(boltz*dble(tk))

c Downward rate
        cr(u,l) = statdg(l)/statdg(u)/1.d11*hvkt*
     &              exp(dble(acs(u-l))
     &          + dble(bcs(u-l))*(hvkt**.25d0)
     &                + dble(ccs(u-l))*(hvkt**.5d0)
     &           ) 


c Upward rate
        cr(l,u) = cr(u,l)*(statdg(u)/statdg(l))*exp(-hvkt)


c        write(6,1000) u-1,l-1,cr(u,l)
c        if (myid .lt. 2) write(22,1000) u-1,l-1,cr(u,l)
c 1000        format(t1,i3,t6,i3,1p,t12,d12.6)

 1        continue


        return
        end

        subroutine colratco(tk)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,hvkt,ct,cr,aline,freqline,
     &  guline,glline

        integer u,l

        include 'nlines_f77.h'

        common /procid/ myid
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &         statdg(nstate),aline(nlines),freqline(nlines),
     &         guline(nlines),glline(nlines)


        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /coll/ cr(nstate,nstate)


             dimension aparaco(6),bparaco(6),aorthoco(6),borthoco(6)

        irate = 1

c        write(6,*) 'COLRATCO'

        if (irate .eq. 1) then
c These are the coefficients for the Flower and Launey rates

c        if (myid .lt. 2) write(22,*) 'Using Flower and Launey rates'

        aparaco(1) = 1.786d-10
        aparaco(2) = 2.528d-10
        aparaco(3) = 7.531d-11
        aparaco(4) = 8.170d-11
        aparaco(5) = 6.302d-11
        aparaco(6) = 3.521d-11

        bparaco(1) = 1.412d0
        bparaco(2) = 1.199d0
        bparaco(3) = 1.244d0
        bparaco(4) = 1.305d0
        bparaco(5) = 1.684d0
        bparaco(6) = 1.427d0

        aorthoco(1) = 2.403d-10
        aorthoco(2) = 3.553d-10
        aorthoco(3) = 8.814d-11
        aorthoco(4) = 1.053d-10
        aorthoco(5) = 4.422d-11
        aorthoco(6) = 3.983d-11

        borthoco(1) = 1.809d0
        borthoco(2) = 1.076d0
        borthoco(3) = 1.336d0
        borthoco(4) = 0.9683d0
        borthoco(5) = 1.167d0
        borthoco(6) = 0.9850d0

        else

c These are the coefficients for the Green and Thaddeus rates

c        if (myid .lt. 2) write(22,*) 'Using Green and Thaddeus rates'

        aparaco(1) = 1.660d-10
        aparaco(2) = 2.800d-10
        aparaco(3) = 1.190d-10
        aparaco(4) = 1.000d-10
        aparaco(5) = 1.300d-10
        aparaco(6) = 0.000d 00

        bparaco(1) = 1.670d0
        bparaco(2) = 1.470d0
        bparaco(3) = 1.850d0
        bparaco(4) = 1.550d0
        bparaco(5) = 2.240d0
        bparaco(6) = 0.000d0

        aorthoco(1) = 1.660d-10
        aorthoco(2) = 2.800d-10
        aorthoco(3) = 1.190d-10
        aorthoco(4) = 1.000d-10
        aorthoco(5) = 1.300d-10
        aorthoco(6) = 0.000d 00

        borthoco(1) = 1.670d0
        borthoco(2) = 1.470d0
        borthoco(3) = 1.850d0
        borthoco(4) = 1.550d0
        borthoco(5) = 2.240d0
        borthoco(6) = 0.000d0

        endif

c Set all rates to zero
        do 3 u = 1,nstate
        do 3 l = 1,nstate
                cr(u,l) = 0.
 3        continue

c Transitions which are more than 6 levels apart will remain
c at a rate of zero. Transitions to same level also zero rate.
c Remember that the index numbers u,l are not the quantum
c numbers. The ground state quantum number 0 corresponds to 
c index number 1.

c        write(6,*) 'statdg',statdg
c        if (myid .lt. 2) write(22,1002) tk
 1002   format(t1,'  u',t6,'  l',t12,'para',t24,'ortho',t35,'downward',
     &          t51,'upward',t60,'T = ',f5.0)


        do 1 l = 1,nstate-1
c        if (myid .lt. 2) write(22,1001)
 1001        format(/)
        do 1 u = l+1,min(l+6,nstate)

        hvkt = planck*freq(u,l)/(boltz*dble(tk))
c        write(6,*) 'hvkt ',hvkt

        crpara = dble(aparaco(u-l))*statdg(l)/statdg(u) *(1.d0 + hvkt)
     &                         *exp(-dble(bparaco(u-l))*sqrt(hvkt) )
        crortho = dble(aorthoco(u-l))*statdg(l)/statdg(u) *(1.d0 + hvkt)
     &                  *exp(-dble(borthoco(u-l))*sqrt(hvkt) )

        cr(u,l) = 0.25d0*crpara + 0.75d0*crortho

        cr(l,u) = cr(u,l)*(statdg(u)/statdg(l))*exp(-hvkt)

c        if (myid .lt. 2) write(22,1000) u-1,l-1,crpara,crortho,cr(u,l),cr(l,u)
 1000        format(t1,i3,t6,i3,1p,t12,e8.1,t24,e8.1,2(4x,e13.6))

 1        continue

c        write(6,*) 'finished loop'

        return
        end



        integer function DLUDCMP(A,N,NP,INDX,D)
        implicit double precision (a-h,o-z)
        PARAMETER (NMAX=100,TINY=1.0d-20)
        DIMENSION A(NP,NP),INDX(N),VV(NMAX)
        D=1.d0
        DO 12 I=1,N
                AAMAX=0.d0
                DO 11 J=1,N
                        IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
 11             CONTINUE
                IF (AAMAX.EQ.0.d0) then
                   write (22,*) 'Singular matrix in fnct dludcmp.'
                   write (6,*) 'Singular matrix in fnct dludcmp.'
                   dludcmp = 0
                   return
                endif
                VV(I)=1./AAMAX
 12     CONTINUE
        DO 19 J=1,N
                IF (J.GT.1) THEN
                        DO 14 I=1,J-1
                                SUM=A(I,J)
                                IF (I.GT.1)THEN
                                        DO 13 K=1,I-1
                                                SUM=SUM-A(I,K)*A(K,J)
 13                                     CONTINUE
                                        A(I,J)=SUM
                                ENDIF
 14                     CONTINUE
                ENDIF
                AAMAX=0.d0
                DO 16 I=J,N
                        SUM=A(I,J)
                        IF (J.GT.1)THEN
                                DO 15 K=1,J-1
                                        SUM=SUM-A(I,K)*A(K,J)
 15                             CONTINUE
                                A(I,J)=SUM
                        ENDIF
                        DUM=VV(I)*ABS(SUM)
                        IF (DUM.GE.AAMAX) THEN
                                IMAX=I
                                AAMAX=DUM
                        ENDIF
 16             CONTINUE
                IF (J.NE.IMAX)THEN
                        DO 17 K=1,N
                                DUM=A(IMAX,K)
                                A(IMAX,K)=A(J,K)
                                A(J,K)=DUM
 17                     CONTINUE
                        D=-D
                        VV(IMAX)=VV(J)
                ENDIF
                INDX(J)=IMAX
                IF(J.NE.N)THEN
                        IF(A(J,J).EQ.0.)A(J,J)=TINY
                        DUM=1.d0/A(J,J)
                        DO 18 I=J+1,N
                                A(I,J)=A(I,J)*DUM
 18                     CONTINUE
                ENDIF
 19     CONTINUE
        IF(A(N,N).EQ.0.d0)A(N,N)=TINY
        dludcmp = 1
        RETURN
        END


        SUBROUTINE DLUBKSB(A,N,NP,INDX,B)
        implicit double precision (a-h,o-z)
        DIMENSION A(NP,NP),INDX(N),B(N)
        II=0
        DO 12 I=1,N
                LL=INDX(I)
                SUM=B(LL)
                B(LL)=B(I)
                IF (II.NE.0)THEN
                        DO 11 J=II,I-1
                                SUM=SUM-A(I,J)*B(J)
 11                     CONTINUE
                ELSE IF (SUM.NE.0.) THEN
                        II=I
                ENDIF
                B(I)=SUM
 12     CONTINUE
        DO 14 I=N,1,-1
                SUM=B(I)
                IF(I.LT.N)THEN
                        DO 13 J=I+1,N
                                SUM=SUM-A(I,J)*B(J)
 13                     CONTINUE
                ENDIF
                B(I)=SUM/A(I,I)
 14     CONTINUE
        RETURN
        END





        subroutine dmattmt(a,msg,n,lun)
c Routine prints out a matrix of double precision
c values with a message.
c Only the first 16 columns are printed. Columns 1 to 8 are in
c the first block, and 9 to 16 in the second
        character * 80 msg
        double precision  a(n,n)

        write(lun,102)
        write (lun,100) msg
 100    format(a)
        write(lun,*) 
     &          'A(i,j): state (i) down ; state (j) across ',
     &          'first block is columns 1-8, 2nd block 9-16'

        do 1 i = 1,n
        write (lun,101) (a(i,j),j = 1,min(n,8))
 1      continue
        write (lun,102)
 102    format (/)
        if (n .gt. 8) then
        do 2 i = 1,n
 2      write (lun,101) (a(i,j),j = 9,min(n,16))
 101    format(1x,8(1pd10.3))
        write (lun,102)
        endif
        return
        end


        subroutine dmultmt(a,x,ncolsA,nrowsA,ncolsX,b)

c Double precision routine. See comments in single precision
c version MMULTMT

        double precision a,x,b

        dimension a(ncolsA,nrowsA),x(ncolsX,ncolsA),
     &          b(ncolsX,nrowsA)

        do 1 i = 1,nrowsA
            do 3 j = 1,ncolsX
c                write(6,*) 'i,j',i,j
                b(i,j) = 0.
                do 2 k = 1,ncolsA
c                        write(6,*) 'k,a,x,b',k,a(j,k),x(k,j),b(i,j)
                        b(i,j) = b(i,j) + a(i,k)*x(k,j)
 2              continue
 3           continue
 1      continue

        return
        end


      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
c       write(6,*) 'xa',xa
c       write(6,*) 'ya',ya
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

        function newrt(iray,line,icount,iter,
     &                dxdu,dydu,dzdu,
     &                vxray,vyray,vzray,
     &                wdthray,pathray,
     &                cont_emis,cont_opac,
     &                emis1d,opac1d,
     &                apL1d,barJ1d,xnorm1d,yspec)

c line          is the current line (in fortran indexing)
c vxray, etc    vectors with the x,y,z components of velocity 
c                       at each vox on the ray. 
c wdthray       line width at each vox on the ray
c emisray       emissivity at each vox on the ray
c opacray       opacity at each vox on the ray
c dxdu,etc      scalar in this function: angle the ray makes 
c                       with x,y,z axes
c vwmin,vwmax   min and max linewidth in velocity units, cm/s
c chanwd        channel width 

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline
        logical doprint
        double precision tau,exptau,opacity,emissivity
        double precision aulhyp,bulhyp,bluhyp

        include 'nlines_f77.h'
        parameter(maxch=2000)
        parameter(maxhyp=50)
        parameter(maxnw=500)
        parameter(maxvox=8192)
        parameter (maxtrans=500)

        common /procid/ myid
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /vgrid/ chvel(maxch,nlines),chfreq(maxch,nlines),
     &      nchan(nlines),profile0(maxch,maxnw),chanwd,
     &      vwmin,vwmax
c23456789112345678921234567893123456789412345678951234567896123456789712
        common /cmb/ ycmb(nlines)
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &         statdg(nstate),aline(nlines),freqline(nlines),
     &         guline(nlines),glline(nlines)
        common /mn/ molecule,fnh3,atoms,brot,dipole



        dimension profhyp(maxch,maxhyp)
        dimension profile(maxch,maxhyp)
        dimension apLray(maxhyp,maxvox),barJray(maxhyp,maxvox)
        dimension xnormray(maxhyp,maxvox)
        dimension vxray(maxvox),vyray(maxvox),vzray(maxvox)
        dimension emis1d(maxvox*maxhyp),opac1d(maxvox*maxhyp)
        dimension emisray(maxhyp,maxvox),opacray(maxhyp,maxvox)
        dimension pathray(maxvox),wdthray(maxvox),cont_emis(maxvox)
        dimension cont_opac(maxvox)

        dimension apL1d(maxvox*maxhyp),barJ1d(maxvox*maxhyp)
        dimension xnorm1d(maxvox*maxhyp)

        dimension yline(maxch,maxvox),yspec(maxch)


        newrt = 1

        iprint = 0
        iprint6 = 0
        doprint = .false.
        lun22 = 22 + myid
c        write(lun22,*) 'myid iray line = ',myid,iray,line
        if (iray .eq. 333 .and. line .eq. 1 .and. iter .eq. 18)  then
          iprint = iray + 1
          iprint6 = 0
          doprint = .true.
        endif
        doprint = .false.

c This shuts off the debugging print statements
c        doprint = .false.

c Cannot have a specific iray and the if myid statement below.
c NEWRT is run by the slaves, and iray = 1646 (center of
c 64 cubed, 3 angles) could be run by any slave number.
c        if (myid .ge. 2) doprint = .false.

        if (doprint) write(lun22,2201) line,iray,iter
 2201   format('line =',i4,'  ray =',i5,'  iter =',i4)

c Zero out apL and Jbar so that we can accumulute the contributions
c from the forward and back directions.
c Unpack the molecular emissivities and opacities 

        if (doprint) then
            write(lun22,*) 'emissivities and opacities for line ',line
            write(lun22,*) 'icount nhyp ',icount,nhyp(line)
        endif

        nh = 1
        if (nltehyp .eq. 1) nh = nhyp(line)

        k = 1
        do 34 ivox=1,icount
        do 34 ihyp=1,nh
            apLray(ihyp,ivox) = 0.
            barJray(ihyp,ivox) = 0.
            xnormray(ihyp,ivox) = 0.
            emisray(ihyp,ivox) = emis1d(k)
            opacray(ihyp,ivox) = opac1d(k)
            if (doprint .or. opac1d(k) .lt. 0. 
     &                  .or. opacray(ihyp,ivox) .lt. 0.) 
     &             write(lun22,7888) k,iray,ivox,ihyp,
     &                    emisray(ihyp,ivox),
     &                    opacray(ihyp,ivox),pathray(ivox)
 7888      format('# iray ivox ihyp emis opac path ',
     &               i4,i8,2i4,1p,4e12.4)
            k = k + 1
 34     continue

        if (doprint) then
          do 35 icur = 1,icount
           write(lun22,*) 'cont emis opac ',
     &         cont_emis(icur),cont_opac(icur)
 35       continue
        endif

c In this subroutine we are working on one line, passed
c down through the args

c Move along the ray, add the contribution of each voxel to
c the radiation. The contribution of each voxel is given by
c the simple solution of the radiative transfer equation for
c uniform conditions 
c I_i = I_i-1 * exp(-tau_i) + S_i*(1 - exp(-tau_i).

c For this ray, first compute the radiation forward away
c the observer, and then backward toward the observer.
c The radiation is not the same in both directions because 
c the projection of the gas velocity on the line of sight 
c has the opposite sign. But for the reverse direction
c we can use the other values calculated for the forward
c direction.

c Solve the radiative transfer equation.
c Loop 804 does the forward and then the backward directions.

        do 804 idir = 1,2

c Loop 44 is along the ray.

        do 44 i = 1,icount

c Here is some algebra for counting forward and back

          im1 = i-1
          im2 = i-2

          if (idir .eq. 1) then
            icur = i
            iprev = im1
          endif

          if (idir .eq. 2) then
            icur  = icount - im1
            iprev = icount - im2
          endif

c Now we know which voxel along the ray we are in. Build
c the line profiles for this voxel. PROFHYP is a
c 2D array with the profile for each hyperfine.

          call shiftprofile(iprint,line,iter,idir,icur,
     &            dxdu,dydu,dzdu,
     &            vxray(icur),vyray(icur),vzray(icur),wdthray(icur),
     &            profhyp,istatus
     &    )

          newrt = istatus
          if (istatus .eq. 0 .and. myid .lt. 2) then
            write(6,*) 'shiftprofile failed on line, voxel ',line,icur
            write(lun22,*) 'shiftprofile failed on line, voxel ',
     &          line,icur
            return
          endif

          ymax = 0.
          ichanmax = nchan(line)/2

          do 177 l = 1,nchan(line)

c At each channel, loop over all the hyperfine lines to get the
c total line opacity at this frequency

        opac = 0.0
        emis = 0.0
        do 402 ihyp = 1,nh
          opac = opac + opacray(ihyp,icur)*profhyp(l,ihyp)
          emis = emis + emisray(ihyp,icur)*profhyp(l,ihyp)
 402    continue

c PATHRAY is the column density of H2 per cell in cgs units.
c Here we put in the factor of hv/4pi that was left off the
c calculation of opacity. Since it was also left off the
c emissivity, this factor cancels and does not appear in
c the source function.

        tau = dble(( opac + cont_opac(icur) ) * pathray(icur))
     &      * planck*freqline(line)/(4.d0*pi);

c        write(lun22,*) 'tau ', tau

        if ((opac + cont_opac(icur)) .ne. 0.0) then
            srce = dble((emis + cont_emis(icur))
     &         / (opac + cont_opac(icur)))
     &         * ((c/freqline(line))**2)/(2.d0*boltz)
        else
            srce = 0.d0
        endif

c        if (doprint .and. l .eq. ichanmax) then
c            write(lun22,4020) icur,emis,opac,cont_emis(icur),
c     &         cont_opac(icur),srce
c 4020       format('opac emis ',1i5,1p6e12.4)
c        endif

c This keeps the exponential from overflowing. 
c Never happens, but just in case
        if (tau .lt. -20.) then
c                if (myid .lt. 2) then
                    write(lun22,3000) tau,icur,opac,opacray(1,icur)
     &                   ,cont_opac(icur),profhyp(l,1)
 3000   format(' Tau very negative ',1p,e12.4,
     &          ' icur opac opacray cont_opac path',i5,4e12.4)
c                endif
                tau = -20.d0
        endif
        exptau = dexp(-tau)


c If at the first voxel of the ray, start with the background radiation
            if (i .eq. 1) then
                yprev = ycmb(line)
            else
                yprev = yline(l,iprev)
            endif

         yline(l,icur) = yprev*exptau + srce*(1.d0 - exptau)


c Add to the current value, the contribution to the approximate 
c lambda operator. We add both the contributions from the forward 
c and backward rays

            do 101 ihyp = 1,nh

             if (opac .ne. 0.0) then
              apLray(ihyp,icur)  = apLray(ihyp,icur)
     &          + (1.d0 - exptau)
     &          * profhyp(l,ihyp) * chanwd/c*chfreq(l,line)
     &          * (opacray(ihyp,icur)*profhyp(l,ihyp) / opac)
c           write(lun22,2002) apLray(ihyp,icur),(1.0 - exptau),
c     &            profhyp(l,ihyp),opacray(ihyp,icur),opac
c 2002      format('APL ',1p,6e12.4)
             endif

c Add to the average brightness in the voxel. 

              barJray(ihyp,icur) = barJray(ihyp,icur)
     &          + yline(l,icur) 
     &          * profhyp(l,ihyp) * chanwd/c*chfreq(l,line)

              xnormray(ihyp,icur) = xnormray(ihyp,icur) 
     &          + profhyp(l,ihyp) * chanwd/c*chfreq(l,line)

c            write(lun22,1050) barJray(ihyp,icur),yline(l,icur),
c     &                     profhyp(l,ihyp),chanwd,c,chfreq(l,line)
c 1050       format(1p,6e10.2)

 101        continue

            if (doprint .and. l .eq. ichanmax
c     &         .and. line .eq. 1
c     &         .and. idir .eq. 1 .and. (icur .eq. 10 .or. icur .eq. 18)
     &         ) then
 2000       format('SS ',4i4,1p9e10.3)
               if (iprint6 .ne. 0) then
               write(6 ,2000) idir,ichanmax,icur,iprev,yprev,
     &            yline(l,icur),srce,
     &            exptau,tau,
     &            barJray(1,icur)
               endif
               write(lun22,2000) idir,ichanmax,icur,iprev,yprev,
     &            yline(l,icur),srce,
     &            exptau,tau,
     &            barJray(1,icur),apLray(ihyp,icur)
            endif

c Loop 177 is the loop over channels
 177        continue

c        write(6,*) 'Finished loop over channels '
            do 178 l = 1,nchan(line)
            if (yline(l,icur) .gt. ymax) then
               ymax = yline(l,icur)
               ichanmax = l
            endif
 178        continue


c Loop 44 is the loop icur 1 to icount for one ray
 44        continue

c        write(6,*) 'finished loop 44 for one ray '

c End of the loop over the forward and back directions IDIR
 804    continue



c        write(6,*) 'loops 44 and 804: voxels= ',icount-1

c YSPEC is the line brightness with the initial CMB subtracted off.
c This is the outgoing spectrum. Along with apLray and barJray, yspec
c is transferred to cmain.c through the arguments of F77 function NEWRT.
c There it is packed and sent from the slave to the master.
        yspecmin = 1.e20
        yspecmax = -1.e20
        do 649 k = 1,nchan(line)
            yspec(k) = yline(k,1) - ycmb(line)
c            if (myid .lt. 2 .and. doprint)
c     &           write(lun22,*) 'output ',k,yspec(k)
            if (yspec(k) .gt. yspecmax) yspecmax = yspec(k)
            if (yspec(k) .lt. yspecmin) yspecmin = yspec(k)
 649    continue

        if (doprint) write(lun22,*) 'yspec min max ',yspecmin,yspecmax

c apLray and barJray are 2D vectors (maxhyp,maxvox) that
c are output by NEWRT. These contain the acceleration
c factor and the average brightness calculated in NEWRT.
c The 2D arrays are packed into the 1D vectors apL1d and barJ1d 
c for transfer to the C program. Stay with the ordering of
c the NHYP index varying most rapidly. We are working on one
c line only, so NHYP(LINE) is constant even though indexed
c on LINE.
        m = 1
        do 240 j = 1,icount
        do 240 i = 1,nh
            k = i + (j-1)*nh
            if (k .ne. m) write(6,*)'NEWRT cannot count, loop 240: ',k,m
            apL1d(k) = apLray(i,j)
            barJ1d(k) = barJray(i,j)
            xnorm1d(k) = xnormray(i,j)

            if (doprint .and. myid .lt. 2) then
               write(lun22,1028) i,j,k,barJ1d(k),apL1d(k),xnorm1d(k)
            endif
 1028       format ('newrt ih ic k bj al xn',3i5,1p5e12.4)

            m = m + 1
 240    continue

        return
        end


        subroutine shiftprofile(iprint,line,iter,idir,ivox,
     &                dxdu,dydu,dzdu,
     &                vx,vy,vz,vwidth,
     &                profhyp,
     &                istatus)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        logical doprint

        integer upper,lower,istartline
        real maxbrt,maxvel,mb1,mb2,mb3,mb4,mb5,mb6,mb7

        include 'nlines_f77.h'

        parameter(maxch=2000)
        parameter(maxnw=500)
        parameter(maxhyp=50)
        parameter(maxtrans=500)


        common /procid/ myid
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav

        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /vgrid/ chvel(maxch,nlines),chfreq(maxch,nlines),
     &      nchan(nlines),profile0(maxch,maxnw),chanwd,
     &      vwmin,vwmax
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /nltehyp/ upper(maxtrans),lower(maxtrans),
     &          istartline(nlines)

c23456789112345678921234567893123456789412345678951234567896123456789712

c        dimension profile(maxch,maxhyp)
        dimension profhyp(maxch,maxhyp),p1(maxch)
c        dimension sumprofile(maxch)
        dimension maxbrt(maxhyp),maxvel(maxhyp)

        chksumerr = 0.1
        istatus = 1

        doprint = .true.
        doprint = .false.
        if (myid .ge. 2) doprint = .false.

        iprint6 = 0

c The theoretical line profile function is a Gaussian. In this code,
c rather than using the exponential function to calculate the line 
c profile, we interpolate from a set of previously calculated and 
c stored Gaussian profiles. This is supposed to save time. The line 
c profile function has to be calculated twice (forward and reverse 
c directions) for each cellalong each ray.

c Start building the profile. First interpolate by width.
c In units of velocity, the line width is the same for all lines 
c regardless of frequency.

c Interpolate in width.  Find the profiles with the widths greater 
c and less than the width VWIDTH. 

c VWIDTH is the line width in the current cell. This comes from the
c model through the MPI communication. In the model, the line width 
c depends on temperature and the microturbulence specified for
c the current cell.

c VWMIN is set as LWMIN in setup.c. It is the minimum width expected.
c The number of Gaussian profiles, computed in function channels is
c nwidths = (vwmax - vwmin)/chanwd + 2. So IW is the array index for
c the stored profile that is just larger than the width in the current
c cell.
        iw = (vwidth - vwmin)/chanwd + 1
        vwlo = vwmin + (iw-1)*chanwd
        vwhi = vwmin + (iw  )*chanwd

c23456789112345678921234567893123456789412345678951234567896123456789712
        if(doprint)write(22,*) 'vwidth,vwmin,vwmax',vwidth,vwmin,vwmax
        if(doprint)write(22,*) 'nchan,chanwd,iw',nchan(line),chanwd,iw

        if (iw .lt. 0 .and. myid .lt. 2) then
            write(22,*) 'Looks like the line width is less than ',
     &          ' the minimum that was set in setup.c '
            write(22,*) 'Is this going to work? Probably not.'
            write(22,*) 'vwidth = ',wvidth, 'vwmin = ',vwmin,' iw= ',iw
        endif   

c23456789112345678921234567893123456789412345678951234567896123456789712
c Linear interpolation between the grid of profiles of different width

c This first case should not occur, but to avoid dividing by zero
        if (abs(vwmax - vwmin) .lt. 1.e-2*chanwd) then
            do 2 i = 1,nchan(line)
                p1(i) = profile0(i,iw) 
 2          continue
        else

        f1 = (vwidth - vwlo)/(chanwd)
        if (doprint .and. myid .lt. 2) 
     &      write(22,*) 'f1,vwidth,vwlo ',f1,vwidth,vwlo

            do 1 i = 1,maxch
               p1(i) = 
     &           f1 * (profile0(i,iw+1) - profile0(i,iw))
     &           +  profile0(i,iw)
 1          continue
        endif

c This normalization applies to both cases above.

        fwidth = vwidth / c * freqline(line)
        do 604 i = 1,maxch
           p1(i) = p1(i)/(sqrtpi*fwidth)
 604    continue

c The profile p1 is now normalized. The integration involves
c multiplication by the line width in frequency units, and
c summation.

c23456789112345678921234567893123456789412345678951234567896123456789712
c This is some printout
        if (if iprint .ne. 0 .and. doprint) then
            write(22,*) 'Check interpolation by width' 
            write(22,*) 'fwidth,chanwd_freq ',fwidth,
     &        (chanwd/c*freqline(line))
            chksum = 0.
            do 303 i = 1,maxch
               chksum = chksum + p1(i)*(chanwd/c*freqline(line))
               write(22,1200) i,profile0(i,iw),p1(i),profile0(i,iw+1),
     &           chksum
 303        continue
        endif
 1200   format(i5,1p,4e12.3)

c The profile p1 should never be negative anywhere
c Comment this out after debugging
        do 500 i = 1,maxch
           if (p1(i) .lt. 0) then
               write(22,*) 'Negative line profile '
               do 504 l = maxch/2 - maxch/4,maxch/2 + maxch/4
                  write(22,1061) l,p1(l)
 504           continue
               istatus = 0
               return
           endif
 500    continue

c Here is the check on the normalization
        chksum = 0.
        do 304 i = 1,maxch
           chksum = chksum + p1(i)*(chanwd/c*freqline(line))
 304    continue

        if (doprint) then
          write(22,*) 'first check sum on interp. profile ',chksum
c          write(6 ,*) 'first check sum on interp. profile ',chksum
        endif

c This IF statement will catch both NAN or a bad profile
        if (abs(chksum - 1.0) .lt. 0.01) goto 404

            write(6 ,*) 'Bad width-interpolated profile ',chksum,
     &          freqline(line)*chanwd/c
            write(22,*) 'Bad width-interpolated profile ',chksum,
     &          freqline(line)*chanwd/c
            write(6 ,*) 'line ',line
            write(22,*) 'line ',line
            do 405 l = maxch/2 - maxch/4,maxch/2 + maxch/4
               if (myid .lt. 2) write(22,1061) l,p1(l)
 1061       format(i5,1p4e12.4)
 405        continue
            istatus = 0
            return

 404    continue
c --------------------------------------------------------------

c --------------------------------------------------------------

c Now we need the velocity in the voxel projected onto the ray.

c Take the projection along line of sight S.
c Because we started at the back of the cloud, S negative,
c and the observer is looking at the front of the cloud,
c velocities which are positive are blueshifted or negative.
c Switch the sign.

            vls =  -( 
     &        vx*dxdu
     &            + vy*dydu 
     &            + vz*dzdu )

        if (doprint .and. myid .lt. 2) then
            write(22,*) 'vx,vy,vz ',vx,vy,vz
            write(22,*) 'dxdu,dydu,dzdu',dxdu,dydu,dzdu
            write(22,*) 'The line of sight velocity is ',vls
        endif
c ---------------------------------------------------------

c Interpolate in frequency

c23456789112345678921234567893123456789412345678951234567896123456789712
c The output profile will have 3 dimensions which are the channel 
c number, direction, and hyperfine. The input profile0 also had 
c 3 dimensions, but these were the channel number, width, and 
c hyperfine. (Outside of this subroutine, profile has 4 dimensions. 
c See notes at the beginning of this routine.)

c The profile, PROFILE, is no longer used in the calculation and
c can be eliminated. 

c PROFHYP is the profile for each hyperfine, shifted by the
c hyperfine velocity and the gas velocity. Each hyperfine profile has
c the same height and is normalized to 1. In other words, the
c relative intensities are not yet applied.

c SUMPROFILE is used in checking.

c Zero the output profile
        do 4 ihyp = 1,nhyp(line)
        do 4 i = 1,nchan(line)
c                profile(i,ihyp) = 0.  
                profhyp(i,ihyp) = 0.  
c                sumprofile(i) = 0.
 4        continue

c Find the velocity shift including the gas velocity and hyperfine
c velocity.

        do 33 ihyp = 1,nhyp(line)

                vlsdir =  vls 
                if (idir .eq. 2) vlsdir = -vls
                vlsdir = vlsdir + hypvel(ihyp,line)

c In the output profile, find the  channel number of 
c the line of sight velocity. Call this IV. This is where we
c want to put the center of our input template profile.

c We will add the input profile, multiplied by the
c relative intensity and shifted by IV to the output profile

c If the velocity grid were uniform, then we could calculate
c the channel number, IV, by dividing by the channel width.
c Because the velocity grid is made up of bands around the
c hyperfine lines each with its own starting velocity, 
c then we might as well search. LOCATE gives us the
c channel number just below VLSDIR.

c chvel(1,line) is the pointer to the start of the vector
c of velocities for line. The vector has length nchan(line).
c Search for the index corresponding to velocity vlsdir.
c The index is iv

       iloc = locate(chvel(1,line),nchan(line),vlsdir,iv)

c Find the remaining fraction between the desired shift and
c the channel number just below VLSRDIR for linear interpolation
        f2 = (vlsdir - chvel(iv,line) )/chanwd

c f2 has to be between zero and one
        if (f2 .lt. 0.0 .or. f2 .gt. 1.0) then
        write(6,*) 'WARNING: The velocity range velrange is too small'
        write(22,*) 'WARNING: The velocity range velrange is too small'
            write(22,1201)  idir,line,ihyp,nchan(line)
            write(22,1202) vls,hypvel(ihyp,line),vlsdir
            write(22,1203) chvel(iv,line),chvel(iv+1,line),iv,f2
 1201  format('idir line ihyp nchan ', 4i5)
 1202  format('vls, hypvel, vlsdir ',1p3e12.3)
 1203  format('chvel1, chvel2, ichan, fraction ',1p2e12.3,0p,i5,1pe12.3)
         endif

c23456789112345678921234567893123456789412345678951234567896123456789712
c If I is the index of the output profile, then we want to shift
c by IV from the center of the input profile.
        i1 = maxch/2 - iv + 1

        if (doprint) write(22,*) 'iv,i1,f2',iv,i1,f2

c Here we loop over the entire velocity grid. This might include
c more than 1 band with gaps between the bands. Looping over
c the entire grid is only correct if the bands are contiguous
c in velocity. But even if they are not, it should not matter
c because the width + shift of a line should not exceed the
c vrange used to set the band boundaries.

c Dont let IP go below 1 or IP+1 beyond MAXCH
c PROFHYP for each hyperfine is normalized to 1
        do 3 i = 1,nchan(line)
            ip = i1 + (i-1)
            if (ip .gt. 0 .and. ip .lt. maxch) then
c               profile(i,ihyp) = 
c     &          + relint(ihyp,line)
c     &          * (f2*( p1(ip+1) - p1(ip) ) + p1(ip))
                profhyp(i,ihyp) = 
     &            (f2*( p1(ip+1) - p1(ip) ) + p1(ip))
               if (profhyp(i,ihyp) .lt. 0.) then
               write(22,*) 'FATAL ERROR: Negative line profile function'
               write(22,*) 'The velocity range is probably too small '
               write(22,*) 'Process ',myid,' sending STOP signal '
               write(6,*) 'FATAL ERROR: Negative line profile function'
               write(6,*) 'The velocity range is probably too small '
               write(6,*) 'Process ',myid,' sending STOP signal '
                   write(22,3001) i,i1,ip,f2,p1(ip),p1(ip+1),
     &                profhyp(i,ihyp)
 3001     format('Negative line profile: i,i1,ip,f2,p1(ip),p1(ip+1) ',
     &            3i6,1p,4e12.4)
                   istatus = 0
                   return     
               endif
            endif
 3        continue
 33       continue

c        do 98 ihyp = 1,nhyp(line)
c        do 6 i = 1,nchan(line)
c            sumprofile(i) = sumprofile(i) + profhyp(i,ihyp)
c 6      continue
c 98      continue
c          if (line .eq. 1) then
c            do 198 i = 1,nchan(line)
c            write(24,*) i,sumprofile(i)
c 198        continue
c          endif

c Loop 3 is over channels i
c Loop 33 is over hyperfines ihyp


c Now make a profile for each hyperfine with the factor relint
c The sum of these profiles over all hyperfines is normalized to 1
       if (nltehyp .eq. 0 .and. nhyp(line) .gt. 1) then

c This multiplies the profile of the first hyperfine by relint
           do 887 i = 1,nchan(line)
               profhyp(i,1) = profhyp(i,1)*relint(1,line)
 887       continue

c For the 2nd and higher hyp profiles, multiply by relint
c and add to the first profile which is accumulating the total
           do 888 i = 1,nchan(line)
           do 888 ihyp = 2,nhyp(line)
               profhyp(i,1) = profhyp(i,1)
     &           + relint(ihyp,line)*profhyp(i,ihyp)
 888       continue
       endif

       if (nltehyp .eq. 0) then
           if (doprint) then
           do 889 i = 1,nchan(line)
               write(22,1088),line,i,chvel(i,line)/1.e5,
     &            profhyp(i,1)
 889       continue
 1088  format('single profile ',2i5,1p3e12.4)
           endif
       endif

       istatus = 1

        return
        end





        subroutine checkprofile(nchan,line,
     &                chvel,chanwd,chfreq,
     &                dxdu,dydu,dzdu,
     &                vx,vy,vz,vwidth,
     &                profile
     &        )

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        logical printgeom

        include 'nlines_f77.h'
        parameter(maxch=2000)
        parameter(maxvox=8192)
        parameter(maxhyp=50)
        parameter(maxnw=500)

        common /procid/ myid
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)

c This section is used only to check the accuracy of the
c SHIFTPROFILE routine, but is too slow to use normally.

        dimension chvel(maxch,nlines),chfreq(maxch,nlines)
        dimension testprofile(maxch)
        dimension nchan(nlines),profile(maxch,2)

        printgeom = .true. 
        printgeom = .false.
        if (myid .ge. 2) printgeom = .false.

            vls =  (
     &        vx*dxdu
     &      + vy*dydu
     &      + vz*dzdu )

        do 406 l = 1,nchan(line)
            testprofile(l) = 0.
 406    continue

        do 1 idir = 1,2
        is = 1
        if (idir .eq. 2) is = -1

        do 407 ihyp = 1,nhyp(line)
        if (printgeom) write(22,*) 'l,nhyp,hv',line,
     &          nhyp(line),hypvel(ihyp,line)

        do 5 l = 1,nchan(line)
                testprofile(l) = 0.
 5        continue

        sum = 0.
        do 404 l = 1,nchan(line)

          testprofile(l) = testprofile(l)
     &     + relint(ihyp,line)
     &       * exp( -( (is*vls + hypvel(ihyp,line) - chvel(l,line))
     &     / vwidth )**2 )
     &     / sqrtpi/(vwidth / c*freqline(line))

          sum = sum + chfreq(l,line)*chanwd/c 
     &          * abs(testprofile(l) - profile(l,idir))


 404    continue
 407    continue

        checksum = 0.
        do 4 l = 1,nchan(line)
                checksum = checksum + 
     &         chfreq(l,line)*chanwd/c *  testprofile(l) 
 4        continue


        if (printgeom) then
        if (idir .eq. 1)
     &        write(22,*) 'testprofile checksum forward',checksum
        if (idir .eq. 2)
     &        write(22,*) 'testprofile checksum reverse',checksum
        endif

        if (sum .gt. 1.e-2) then
                if (idir .eq. 1)
     &             write(22,*) 'forward profiles do not match'
                if (idir .eq. 2)
     &             write(22,*) 'reverse profiles do not match'
             write(22,*) 'number of channels for line',line,' is ',
     &                nchan(line)
             write(22,*) 'vwidth',vwidth
             do 304 l = 1,nchan(line)
                write(22,*) 
     &                        l,
     &                        testprofile(l),
     &                        profile(l,idir)
 304         continue
             continue
             stop
        else
                if (printgeom) write(22,*) 'profiles checked ok'
        endif

 1        continue
c ---------------------------------------------------------------

        return
        end

c23456789112345678921234567893123456789412345678951234567896123456789712

        function copyatomic(freqlinx,alinx,stdgx,jdx,kdx,ldx,
     &    indux,indlx,nhypx,istartx,upperx,lowerx,aulhyx)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline,
     &  freqlinx,alinx,stdgx,gulineC,gllineC
        double precision aulhyp,bulhyp,bluhyp,aulhyx
        integer u,l,upper,lower,istartline,upperx,lowerx

        logical printgeom

        include 'nlines_f77.h'
        parameter(maxch=2000)
        parameter(maxhyp=50)
        parameter(maxtrans=500)

        common /procid/ myid
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          gulineC(nlines),gllineC(nlines)
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /lines/ indexu(nlines),indexl(nlines)
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /nltehyp/ upper(maxtrans),lower(maxtrans),
     &          istartline(nlines)
        common /einstein/ aulhyp(maxhyp,nlines),
     &                      bulhyp(maxhyp,nlines),bluhyp(maxhyp,nlines)



        dimension freqlinx(nlines),alinx(nlines),stdgx(nstate)
        dimension jdx(nstate),kdx(nstate),ldx(nstate),indux(nlines),
     &      indlx(nlines),nhypx(nlines),istartx(nlines)
        dimension upperx(maxtrans),lowerx(maxtrans)
        dimension aulhyx(maxhyp,nlines)


        copyatomic = 1

c This subroutine does 2 things.

c 1) Copies the frequency and Einstein A of the line
c from the Fortran part of the program to the C part. 
c The information is passed through the function arguments. 

c Load the function arguments. Do this for all lines including N2H+
c with NLTE hyperfines. 
c The variables INDEXU and INDEXL are dimensioned with NLINES.
c In the case of NLTE hyperfines, these refer to the main hyperfine
c transition. 
c UPPER and LOWER are dimensioned with MAXTRANS, the number of hyperfine
c transitions. These refer to the upper and lower states of each
c hyperfine transition. If there is no hyperfine splitting, then
c INDEXU and INDEXL are the same as UPPER and LOWER except that
c INDEXU, INDEXL are indexed from 1 (F77 style)
c and UPPER and LOWER beginning from 0 (C or QN style).

c In copying INDEXL and INDEXU to the function arguments where they
c are passed to the C program, subtract 1 to make them correspond to 
c C indexing.

c If the molecule does not have NLTE structure, then aulhyp(ihyp,line)
c might not be defined at all (= zero), or aulhyp(1,line) = aline(line)
c and zero for ihyp > 1.

          do 4 i = 1,maxtrans
             upperx(i) = upper(i)
             lowerx(i) = lower(i)
 4        continue

        do 3 line = 1,nlines
                indux(line) = indexu(line) - 1
                indlx(line) = indexl(line) - 1
                alinx(line) = aline(line)
                freqlinx(line) = freqline(line)
                nhypx(line) = nhyp(line)
c                write(22,*)'copyatomic: line nhyp nhypx ',
c     &             line,nhyp(line),nhypx(line)
                istartx(line) = istartline(line)
            do 5 ihyp = 1,nhyp(line)
                aulhyx(ihyp,line) = aulhyp(ihyp,line)
c                write(22,*)'copyatomic ',line,ihyp,aulhyx(ihyp,line)
 5          continue
 3      continue

        do 2 n = 1,nstate
                stdgx(n) = statdg(n)
                jdx(n) = jindx(n)
                kdx(n) = kindx(n)
                ldx(n) = lindx(n)
 2        continue


        return
        end

        subroutine colrate_hcn_ground(tk)

c A linear interpolation of the ion rates.
c Using the tables of temperatures and rates in
c common HCNCOLL, interpolates the rates for
c the temperature TK passed as an argument.
c The interpolated rates are stored in the common
c COLL.

c CR                Interpolated collision rates cm^3/sec at 
c                        temperature TK
c NSTATE        number of levels
c TEMPS                temperatures in K
c RATE                collision rates at temperatures in TEMPS
c IT                array index 

c Use this subroutine for HCO+ and N2H+ with no hyperfine
c structure

        integer u,l
       real *8 frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr


        include 'nlines_f77.h'
        common /procid/ myid
        common /coll/ cr(nstate,nstate)
        common /hcncoll/ temps(15),rate(15,30,30)
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav


        it = tk/10. + 1
        it = max(1,min(it,15))

c Interpolate downward rates. Set the upper rates by stat. eq.
        do 2 u = 2,nstate
        do 2 l = 1,u-1

c This is a simple linear interpolation.
                cr(u,l) = 1.d-10 * dble( rate(it,u,l)
     &          + (tk - temps(it))/(temps(it+1) - temps(it))
     &          * (rate(it+1,u,l) - rate(it,u,l)) )

              cr(l,u) = cr(u,l) * statdg(u)/statdg(l) 
     &            * exp( -planck/boltz/tk*freq(u,l) )

        if (cr(u,l) .lt. 0.)
     &          write(6,1000) u,l,it,tk,rate(it,u,l),
     &                        rate(it+1,u,l),cr(u,l)

 1000    format('negative collision rate ',3i5,1p,4e12.4)


 2      continue

        return
        end


        subroutine colrateion(tk)

c A linear interpolation of the ion rates.
c Using the tables of temperatures and rates in
c common IONRATES, interpolates the rates for
c the temperature TK passed as an argument.
c The interpolated rates are stored in the common
c COLL.

c These rates are from Flower 1999

c CR                Interpolated collision rates cm^3/sec at 
c                        temperature TK
c NSTATE        number of levels
c TEMPS                temperatures in K
c RATE                collision rates at temperatures in TEMPS
c IT                array index 

c Use this subroutine for HCO+ and N2H+ with no hyperfine
c structure

        integer u,l
        real *8 frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr


        include 'nlines_f77.h'
        common /procid/ myid
        common /coll/ cr(nstate,nstate)
        common /ionrates/ temps(41),rate(41,21,21)
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav


        it = tk/10. + 1
        it = max(1,min(it,40))

        do 2 u = 2,nstate
        do 2 l = 1,u-1

           ju = jindx(u) + 1
           jl = jindx(l) + 1

c Use this for rates between N2H+ hyperfine states
c           gf = (2.*lindx(l) + 1.)/ (9.*(2.*jindx(l) + 1.))
c Use this for the total rates between rotational transitions
           gf = 1.0

c This is a simple linear interpolation.
                cr(u,l) = dble(gf) * dble (rate(it,ju,jl)
     &          + (tk - temps(it))/(temps(it+1) - temps(it))
     &          * (rate(it+1,ju,jl) - rate(it,ju,jl)) )

         cr(l,u) = cr(u,l) * statdg(u)/statdg(l)
     &            * dexp( -planck/boltz/dble(tk)*freq(u,l) )


        if (cr(u,l) .lt. 0.)
     &          write(6,1000) ju,jl,it,tk,rate(it,ju,jl),
     &                        rate(it+1,ju,jl),cr(ju,jl)

 1000    format('negative collision rate ',3i5,1p,4e12.4)


 2      continue

        return
        end

        subroutine colrate_nhyp(tk,iprint)


c This subroutine takes the collision rates for N2H+ deltaJ
c transitions and divides them up between the hyperfine 
c transitions and interpolates for temperature.

c The rates were loaded into common N2H_HYP in subroutine
c N2H_DELTAJ_RATES


c interpolates the rates for
c the temperature TK passed as an argument.
c The interpolated rates are stored in the common
c COLL.


c CR      Output: Interpolated collision rates cm^3/sec at 
c TK      temperature TK
c NSTATE  number of states
c TEMPS   temperatures in K
c RATE    collision rates at temperatures in TEMPS
c IT      array index 

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline,
     &  energy

        integer u,l

        include 'nlines_f77.h'
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /procid/ myid
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)
        common /coll/ cr(nstate,nstate)
        common /n2h_hyp/ temps(15),rate(15,8,8)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /cenergy/ energy(nstate)

        dimension cbal(nstate,nstate)

        do 1 it=1,15
           if (tk .lt. temps(it)) goto 4
 1      continue
 4      continue
c        if (iprint .eq. 1) write(22,*)'tk,it,temps(it) ',
c     &        tk,it,temps(it)

c tk is between temperature(it-1) and temperature(it) 

c According to my notes and Keto+ 2004, and 
c K&R 2010, to get the hyperfine rate, multiply the 
c total rate by 2F+1 / 9(2J+1) where these are the
c quantum numbers for the final states. Then when
c we sum back using the same factor, but with the
c initial states, we should get the total rate.

c This is a simple linear interpolation. The 4th line
c converts the col rate between J levels into the rate 
c between hyperfine levels. 


c Finally dble the whole thing and multiply by 1.e-10

c ji initial = i
c jf final   = j
        do 2 j = 2,nstate
        do 2 i = 1,j-1
c        do 2 j = 1,nstate
c        do 2 i = 1,nstate

           ji = jindx(j) + 1
           jf = jindx(i) + 1
           gf = (2.*lindx(i) + 1.)/ (9.*(2.*jindx(i) + 1.))

c This should be the downward rate
           cr(j,i) = dble (
     &        rate(it-1,ji,jf)
     &        + (tk - temps(it-1))/(temps(it) - temps(it-1))
     &        * (rate(it,ji,jf) - rate(it-1,ji,jf)) 
     &                      )*dble(gf)*1.d-10

c This should be the upward rate
         if (j .ne. i) then
c         if (energy(j) .gt. energy(i)) then
c            cbal(i,j) = cr(j,i) * statdg(j)/statdg(i)
            cr(i,j) = cr(j,i) * statdg(j)/statdg(i)
     &            * dexp( -planck/boltz/dble(tk)*freq(j,i) )
c         endif
c         if (energy(i) .gt. energy(j)) then
c            cr(j,i) = cr(i,j) * statdg(i)/statdg(j)
c            cbal(i,j) = cr(j,i) * statdg(j)/statdg(i)
c     &            * dexp( planck/boltz/dble(tk)*freq(i,j) )
c         write(6,*) 'Backward rate '
c         endif
c         if (energy(i) .eq. energy(j)) then
c            write(6,*) 'Problem with energies in detailed balance '
c            write(22,*) 'Problem with energies in detailed balance '
c            write(22,*) 'energy ',i,j,ji,jf,energy(i),energy(j)
c         endif
         endif


         
        if (cr(i,j) .lt. 0.)
     &          write(6,*) 'NEG RATE ',tk,rate(it,ji,jf),
     &                        rate(it-1,ji,jf),cr(i,j)

c           if (jindx(j) .eq. 1 .and. jindx(i) .eq. 0
c     &        .and. myid .lt. 2 .and. iprint .eq. 1) 
c     &        write(22,*) 'rate used ',ji,'->',jf,
c     &        rate(it-1,ji,jf)
 2      continue

c        do 32 j = 2,nstate
c        do 32 i = 1,j-1
c        write(22,2020) i,j,ji,jf,cr(i,j),cr(j,i),cbal(j,i)
c 2020   format('check DB ',4i5,1p3e14.6)
c 32     continue

c Make sure that rates between the same levels are zero
        do 3 j=1,nstate
           cr(j,j) = 0.d0
 3      continue

        if (iprint .eq. 1 .and. myid .lt. 2) then
          write(22,*) '------------------------------'
          write(22,*) 'Check collision rates '
          write(22,*) 'Make sure that the hyperfine rates sum to',
     &     ' the rotational rates'
          write(22,*)'Ji   Jf   Tk   sumiHrates  Jrate '
          do 6 k = 1,nlines
          do 6 idir = 1,2
            if (idir .eq. 1) then
              itf = k
              iti = k-1
            else
              iti = k
              itf = k-1
            endif
            sumcr = 0.
            do 5 j = 1,nstate
            do 5 i = 1,nstate
               if (jindx(i) .eq. iti .and. jindx(j) 
     &           .eq. itf) then 
                 gi = (2.*lindx(i)+1.)/(9.*(2.*jindx(i)+1.))
                 sumcr = sumcr + gi*cr(i,j)*1.d10
               endif
 5          continue
 1000       format(i3,' ->',i3,f6.2,1p2e12.4)
            write(22,1000) iti,itf,tk,sumcr, 
     &        rate(it-1,iti+1,itf+1) + (tk-temps(it-1))
     &         * (rate(it,iti+1,itf+1) - rate(it-1,iti+1,itf+1))
     &         / (temps(it) - temps(it-1))
 6      continue
        write(22,*) '------------------------------'

        write(22,*) 'Check the rotational rates: colrate_nhyp',
     &    ' Tk = ',tk,' it = ',it
        write(22,*) 'line   Ju   Jl  freq        CRul       CRlu',
     &              '       CRlu Det.Bal.   Gu    Gl'
        do 10 line = 1,nlines
              ji = line+1
              jf = line
              rif = rate(it-1,ji,jf)
     &        + (tk - temps(it-1))/(temps(it) - temps(it-1))
     &        * (rate(it,ji,jf) - rate(it-1,ji,jf))
              rfi = rate(it-1,jf,ji)
     &        + (tk - temps(it-1))/(temps(it) - temps(it-1))
     &        * (rate(it,jf,ji) - rate(it-1,jf,ji))
 1010   format(3i5,1p4e12.4,2i6)
              write(22,1010) line,ji-1,jf-1,freqline(line),rif,rfi,
     &          rif*float(2*(ji-1)+1)/float(2*(jf-1)+1)
     &          *exp( -planck*freqline(line)/boltz/tk ),
     &          (2*(ji-1)+1),(2*(jf-1)+1)
 10     continue


        endif


        return
        end


        subroutine colrate_hcnhyp2(tk,iprint)



c CR                Interpolated collision rates cm^3/sec at 
c T                 temperature TK
c NSTATE            number of levels
c TEMPS             temperatures in K
c RATE                collision rates at temperatures in TEMPS
c IT                array index 

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        integer u,l

C must define findx here because it is not implicitly an integer array
        integer findx

        include 'nlines_f77.h'
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /procid/ myid
        common /index/ jindx(nstate),findx(nstate),lindx(nstate)
        common /coll/ cr(nstate,nstate)
        common /hcncoll/ temps(15),rates(15,30,30)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)


        do 1 it=1,15
           if (tk .lt. temps(it)) goto 4

c        write(*,*)'tk,it,temps(it) ',
c     &        tk,it,temps(it)

 1      continue
 4      continue

        if (it .lt. 2) it = 2

c tk is between temperature(it-1) and temperature(it) 

c Finally dble the whole thing and multiply by 1.e-10

c ji initial = i
c jf final   = j

        do 2 j = 2,nstate
        do 2 i = 1,j-1
           ji = jindx(j) + 1
           jf = jindx(i) + 1


           gf = (2.*findx(i) + 1.)/ (3.*(2.*jindx(i) + 1.))
c           write(22,1030) j,i,ji,jf,findx(j),jindx(j),gf
c 1030      format('Indices ',6i5,f7.2)

c This should be the downward rate
           cr(j,i) = dble (
     &        rates(it-1,ji,jf)
     &        + (tk - temps(it-1))/(temps(it) - temps(it-1))
     &        * (rates(it,ji,jf) - rates(it-1,ji,jf))
     &                      )*dble(gf)*1.d-10

c This should be the upward rate
            cr(i,j) = cr(j,i) * statdg(j)/statdg(i)
     &            * dexp( -planck/boltz/dble(tk)*freq(j,i) )

        if (cr(i,j) .lt. 0.)
     &          write(6,*) 'NEG RATE ',tk,rates(it,ji,jf),
     &                        rates(it-1,ji,jf),cr(i,j)

c           if (jindx(j) .eq. 1 .and. jindx(i) .eq. 0
c     &        .and. myid .lt. 2 .and. iprint .eq. 1)
c     &        write(22,*) 'rate used ',ji,'->',jf,
c     &        rate(it-1,ji,jf)
 2      continue

c Make sure that rates between the same levels are zero
        do 3 j=1,nstate
           cr(j,j) = 0.d0
 3      continue

c        iprint = 1
        if (iprint .eq. 1 .and. myid .lt. 2) then
          write(22,*) '------------------------------'
          write(22,*) 'Check collision rates '
          write(22,*) 'Make sure that the hyperfine rates sum to',
     &     ' the rotational rates'
          write(6,*)'Ji   Jf   Tk   sumiHrates  Jrate '
          write(6,*) '------------------------------'
          write(6,*) 'Check collision rates '
          write(6,*) 'Make sure that the hyperfine rates sum to',
     &     ' the rotational rates'
          write(6,*)'Ji   Jf   Tk   sumiHrates  Jrate '
          do 6 k = 1,nlines
          do 6 idir = 1,2
            if (idir .eq. 1) then
              itf = k
              iti = k-1
            else
              iti = k
              itf = k-1
            endif
            sumcr = 0.
c            write(6,*) 'iti, itf ',iti,itf
            do 5 j = 1,nstate
            do 5 i = 1,nstate
               if (jindx(i) .eq. iti .and. jindx(j)
     &           .eq. itf) then 
           gi = (2.*findx(i) + 1.)/ (3.*(2.*jindx(i) + 1.))
                 sumcr = sumcr + gi*cr(i,j)*1.d10
               endif
 5          continue

 1000       format(i3,' ->',i3,f6.2,1p2e12.4)
            write(6,*) iti,itf,tk,sumcr,
     &        rates(it-1,iti+1,itf+1) + (tk-temps(it-1))
     &         * (rates(it,iti+1,itf+1) - rates(it-1,iti+1,itf+1))
     &         / (temps(it) - temps(it-1))
            write(22,1000) iti,itf,tk,sumcr,
     &        rates(it-1,iti+1,itf+1) + (tk-temps(it-1))
     &         * (rates(it,iti+1,itf+1) - rates(it-1,iti+1,itf+1))
     &         / (temps(it) - temps(it-1))
 6      continue
        write(22,*) '------------------------------'
        endif


c This is a long print statement for debugging. Compares
c the rates derived from linear interpolation to make sure
c they are in detailed balance.
        if (iprint .ne. 0 .and. j .eq. -1) then 
        write(22,*) 'Check the rotational rates: colrate_hcnhyp2',
     &    ' Tk = ',tk,' it = ',it
        write(22,*) '   line   u  l  freq        CRul       CRlu',
     &              '       CRlu Det.Bal.   Gu    Gl'
        do 20 u = 1,nstate
        do 20 l = 1,nstate

              ji = u
              jf = l
c XXXX sign change in exponential required here from previous version
        crludet = cr(u,l)*exp(planck*freq(ji,jf)/boltz/tk)
     &   *  statdg(u)/statdg(l)
 1010   format('CRDET ',2i5,1p4e12.4)
c        write(22,1010) u,l,freq(ji,jf),cr(u,l)*1.0D10,
c     &       cr(l,u)*1.0D10,crludet*1.0D10
 20     continue
        endif

        return
        end



        function ltepops(tk,frac)

c This computes the lte level populations for CO and other
c dipole molecules.


        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,ric,fsc,a,freq,chfreq,xfrac,statdg,cr,
     &      aline,freqline,guline,glline
        double precision b,q,n,ap1,ap2

        include 'nlines_f77.h'


        common /procid/ myid
        dimension frac(nstate)
        common /flag/ flag
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /mn/ molecule,fnh3,atoms,brot,dipole


c This will add up the first 1000 levels to create the
c partition function
c       q = 0.d0
c       do 1 i = 1,1000
c               n = dble(i-1)
c               q = q + (2.d0*n+1.d0)
c     &                   * exp(-brot*n*(n+1.d0)*planck/(boltz*dble(tk)))
c 1     continue

c Use approximation to the partition function from
c Townes and Schawlow page 20 equations 1-54 and 1-55

        ap1 = boltz*dble(tk)/(planck*brot)
        ap2 = ap1 + 1.d0/3.d0 + 1.d0/15.d0/ap1
     &     + 4.d0/315.d0/(ap1*ap1) + 1.d0/315.d0/(ap1*ap1*ap1)
c       write(6,*) 'partition function ',q,ap1,ap2
c       if (myid .lt. 2) write(22,*) 'partition function ',q,ap1,ap2

        do 2 i = 1,nstate
                n = dble(i-1)
                frac(i) = (2.d0*n+1.d0)
     &            * exp(-brot*n*(n+1.d0)*planck/(boltz*dble(tk)))
     &                  / ap2
 2      continue

        ltepops = 1

        return
        end


        subroutine load_n2h_rates


        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        double precision dekt

        logical doprint

        integer u,l
        include 'nlines_f77.h'

        common /procid/ myid
        common /n2hrates/ temp(5),rate(5,7,7)
        character *80 msg
        double precision cs(7,7)

        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)

        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /coll/ cr(nstate,nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole

        data temp / 5., 10., 20., 30., 40. /

        data rate /
     &  0.0000 , 0.0000 , 0.0000 , 0.0000 , 0.0000 , 
     &  1.0435 , 0.9454 , 0.8660 , 0.8329 , 0.8211 , 
     &  0.5100 , 0.5144 , 0.4657 , 0.4200 , 0.3835 , 
     &  0.2667 , 0.2946 , 0.2835 , 0.2662 , 0.2525 , 
     &  0.1651 , 0.1906 , 0.1989 , 0.1968 , 0.1924 , 
     &  0.0677 , 0.0996 , 0.1162 , 0.1203 , 0.1219 , 
     &  0.0998 , 0.1059 , 0.0958 , 0.0888 , 0.0842 , 
     &  1.2790 , 1.8074 , 2.0639 , 2.1390 , 2.1906 , 
     &  0.0000 , 0.0000 , 0.0000 , 0.0000 , 0.0000 , 
     &  1.8527 , 1.7626 , 1.6240 , 1.5282 , 1.4620 , 
     &  0.7542 , 0.8148 , 0.7721 , 0.7258 , 0.6898 , 
     &  0.3527 , 0.4513 , 0.5044 , 0.5184 , 0.5227 , 
     &  0.2946 , 0.3961 , 0.4258 , 0.4209 , 0.4130 , 
     &  0.3363 , 0.3583 , 0.3320 , 0.3126 , 0.2976 , 
     &  0.1740 , 0.6719 , 1.1900 , 1.3424 , 1.3710 , 
     &  0.5160 , 1.2035 , 1.7359 , 1.8954 , 1.9524 , 
     &  0.0000 , 0.0000 , 0.0000 , 0.0000 , 0.0000 , 
     &  1.3901 , 1.5509 , 1.5515 , 1.5378 , 1.5308 , 
     &  1.0561 , 1.1175 , 1.0348 , 0.9617 , 0.9082 , 
     &  0.4800 , 0.6564 , 0.7070 , 0.6988 , 0.6869 , 
     &  0.6348 , 0.6889 , 0.6400 , 0.5976 , 0.5639 , 
     &  0.0087 , 0.1410 , 0.5189 , 0.7617 , 0.9026 , 
     &  0.0201 , 0.2031 , 0.5892 , 0.8042 , 0.9211 , 
     &  0.1328 , 0.5676 , 1.1108 , 1.3769 , 1.5330 , 
     &  0.0000 , 0.0000 , 0.0000 , 0.0000 , 0.0000 , 
     &  1.5623 , 1.6876 , 1.6612 , 1.6251 , 1.6041 , 
     &  0.6362 , 0.8509 , 0.8886 , 0.8610 , 0.8347 , 
     &  0.6454 , 0.7431 , 0.7287 , 0.7026 , 0.6806 , 
     &  0.0002 , 0.0195 , 0.1908 , 0.3981 , 0.5657 , 
     &  0.0003 , 0.0242 , 0.2023 , 0.4069 , 0.5742 , 
     &  0.0036 , 0.0886 , 0.3918 , 0.6125 , 0.7505 , 
     &  0.0562 , 0.3641 , 0.8765 , 1.1545 , 1.3227 , 
     &  0.0000 , 0.0000 , 0.0000 , 0.0000 , 0.0000 , 
     &  0.8101 , 1.1382 , 1.2776 , 1.3159 , 1.3478 , 
     &  0.8656 , 0.9713 , 0.9294 , 0.8872 , 0.8568 , 
     &  0.0000 , 0.0013 , 0.0446 , 0.1415 , 0.2515 , 
     &  0.0000 , 0.0028 , 0.0680 , 0.1912 , 0.3172 , 
     &  0.0000 , 0.0067 , 0.1061 , 0.2569 , 0.3963 , 
     &  0.0030 , 0.0238 , 0.1865 , 0.3537 , 0.4807 , 
     &  0.0133 , 0.1484 , 0.5099 , 0.7632 , 0.9422 , 
     &  0.0000 , 0.0000 , 0.0000 , 0.0000 , 0.0000 , 
     &  1.0544 , 1.1690 , 1.1494 , 1.1577 , 1.1872 , 
     &  0.0000 , 0.0001 , 0.0113 , 0.0505 , 0.1058 , 
     &  0.0000 , 0.0002 , 0.0164 , 0.0688 , 0.1392 , 
     &  0.0000 , 0.0006 , 0.0297 , 0.1064 , 0.1977 , 
     &  0.0000 , 0.0017 , 0.0472 , 0.1398 , 0.2389 , 
     &  0.0001 , 0.0102 , 0.1146 , 0.2492 , 0.3653 , 
     &  0.0058 , 0.0943 , 0.3548 , 0.5590 , 0.7180 , 
     &  0.0000 , 0.0000 , 0.0000 , 0.0000 , 0.0000  
     &  /

        doprint = .false.


        do 30 i = 1,5
        do 19 u = 1,7
        do 19 l = 1,7
          cs(u,l) = dble(rate(i,u,l))
 19     continue
        lun = 22
        msg = 'original rates'
        if (myid .lt. 2) call dmattmt(cs,msg,7,lun)
 30     continue

        return
        end




        integer function c17ohyp(vwmin)

        include 'nlines_f77.h'
        parameter (maxhyp=50)

        double precision hypfrq(maxhyp,nlines),frq
        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &      grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        common /procid/ myid
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /lines/ indexu(nlines),indexl(nlines)
        dimension tmp1(maxhyp),tmp2(maxhyp),tmp3(maxhyp)
        dimension hypfrqc17o(14,8),nhypc17o(8),relintc17o(14,8)

c        write(6,*) 'Working out C17O hyperfine structure'
c c17ohyp is the return value for the function, but it does nothing
        c17ohyp = 1

        hypfrqc17o( 1, 1) =   112.358774 
        relintc17o( 1, 1) = 0.222 
        hypfrqc17o( 2, 1) =   112.358982
        relintc17o( 2, 1) = 0.444 
        hypfrqc17o( 3, 1) =   112.359999
        relintc17o( 3, 1) = 0.333 
        nhypc17o( 1) =  3
        hypfrqc17o( 1, 2) =  224.713501 
        relintc17o( 1, 2) = 0.040 
        hypfrqc17o( 2, 2) =  224.714045  
        relintc17o( 2, 2) = 0.122
        hypfrqc17o( 3, 2) =  224.714155  
        relintc17o( 3, 2) = 0.171 
        hypfrqc17o( 4, 2) =  224.714199
        relintc17o( 4, 2) = 0.333
        hypfrqc17o( 5, 2) =  224.714214 
        relintc17o( 5, 2) = 0.067 
        hypfrqc17o( 6, 2) =  224.714726 
        relintc17o( 6, 2) = 0.093 
        hypfrqc17o( 7, 2) =  224.715062  
        relintc17o( 7, 2) = 0.016 
        hypfrqc17o( 8, 2) =  224.715172 
        relintc17o( 8, 2) = 0.095 
        hypfrqc17o( 9, 2) =  224.715270 
        relintc17o( 9, 2) = 0.062 
        nhypc17o( 2) = 9
        hypfrqc17o( 1, 3) =  337.060475 
        relintc17o( 1, 3) = 0.011
        hypfrqc17o( 2, 3) =  337.060674 
        relintc17o( 2, 3) = 0.011 
        hypfrqc17o( 3, 3) =  337.060795 
        relintc17o( 3, 3) = 0.008
        hypfrqc17o( 4, 3) =  337.060905  
        relintc17o( 4, 3) = 0.066 
        hypfrqc17o( 5, 3) =  337.060956
        relintc17o( 5, 3) = 0.194  
        hypfrqc17o( 6, 3) =  337.060969 
        relintc17o( 6, 3) = 0.286 
        hypfrqc17o( 7, 3) =  337.061019 
        relintc17o( 7, 3) = 0.054
        hypfrqc17o( 8, 3) =  337.061093 
        relintc17o( 8, 3) = 0.065 
        hypfrqc17o( 9, 3) =  337.061186  
        relintc17o( 9, 3) = 0.037 
        hypfrqc17o(10, 3) =  337.061204 
        relintc17o(10, 3) = 0.122 
        hypfrqc17o(11, 3) =  337.061449 
        relintc17o(11, 3) = 0.069
        hypfrqc17o(12, 3) =  337.061531 
        relintc17o(12, 3) = 0.030 
        hypfrqc17o(13, 3) =  337.061930 
        relintc17o(13, 3) = 0.044 
        hypfrqc17o(14, 3) =  337.062956 
        relintc17o(14, 3) = 0.004 
        nhypc17o( 3) = 14

        if (nlines .gt. 8) then
           write(6,*) 'Cannot use C17O with more than 8 lines '
           write(6,*) 'Program exiting ... '
           c17ohyp = 0
           return
        endif

        if (nlines .gt. 4) then
        line1 = 4
        do 79 line = line1,nlines
                nhypc17o(line) = 1
                relintc17o(1,line) = 1.
                hypfrqc17o(1,line) = freq(indexu(line),indexl(line))
 79     continue
        endif


        do 77 line = 1,nlines
             nhyp(line) = nhypc17o(line)
             if (myid .lt. 2) 
     &         write(22,*) 'line number ',line,' number hyperfines ',
     &               nhyp(line)
             do 78 ihyp = 1,nhyp(line)
                 relint(ihyp,line) = relintc17o(ihyp,line)
                 hypfrq(ihyp,line) = hypfrqc17o(ihyp,line)
                 if (myid .lt. 2) 
     &         write(22,*) 'relint hypfrq ', relint(ihyp,line),
     &                hypfrq(ihyp,line)
 78        continue
 77        continue

c For each of the 3 lowest hyperfines 
c (no hyperfines for higher lines of C17O)
        do 1 line=1,3

c Sum the relative intensities.
c Pick out the strongest hyperfine
        ihypmax = 0
        relintmax = 0.
        sum = 0.

        do 2 ihyp = 1,nhyp(line)

        sum = sum + relint(ihyp,line)

        if (relint(ihyp,line) .gt. relintmax) then
                ihypmax = ihyp
                relintmax = relint(ihyp,line)
        endif

 2      continue

c Normalize the relative intensities to unity and convert
c the frequencies to offsets from the main line. The main
c line is defined as in the JPL catalog even if there is
c no hyperfine line at this frequency. Convert to velocity.

c       if (myid .lt. 2) write(22,*) line,indexu(line),indexl(line),
c     &         freq(indexu(line),indexl(line))

        frq = freq(indexu(line),indexl(line))

        do 3 ihyp = 1,nhyp(line)
          relint(ihyp,line) = relint(ihyp,line)/sum
          hypfrq(ihyp,line) = hypfrq(ihyp,line)*1.d9 - frq
          hypvel(ihyp,line) = -hypfrq(ihyp,line)/frq*c
 3      continue


        if (myid .lt. 2) 
     &    write(22,*) 'new line # ',line,
     &    freq(indexu(line),indexl(line))/1.e9,' GHz'
        if (myid .lt. 2) 
     &    write(22,*) 'number, hyperfine freq shift in MHz, ',
     &     ' shift in kms, relative intensity'
        do 4 ihyp = 1,nhyp(line)
          if (myid .lt. 2) write(22,*) ihyp,hypfrq(ihyp,line)/1.e6,
     &          hypvel(ihyp,line)/1.e5,
     &          relint(ihyp,line)
 4      continue

 1      continue

        return

        end

c23456789112345678921234567893123456789412345678951234567896123456789712
        subroutine n2h_hyp_init(status)

        include 'nlines_f77.h'
        parameter (maxch=2000,maxhyp=50,maxtrans=500)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,gco,aline,freqline,guline,glline,
     &  energy

        integer status,jq,f1q,fq,upper,lower,upper_data,lower_data
        double precision enq,frequency,dpdipole,relative
        double precision aulhyp,bulhyp,bluhyp
        double precision qju,qfu,einstj,einsteina_sum


        common /procid/ myid
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /lines/ indexu(nlines),indexl(nlines)
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /nltehyp/ upper(maxtrans),lower(maxtrans),
     &          istartline(nlines)
        common /linelevs/ jloline(nlines),jupline(nlines)
        common /einstein/ aulhyp(maxhyp,nlines),
     &                      bulhyp(maxhyp,nlines),bluhyp(maxhyp,nlines)
        common /cenergy/ energy(nstate)

c We have data on N2H+ for 64 states and 280 transitions
c The energy levels are from table 1 of 
c Daniel, Cernicharo, Dubernet, 2005, MNRAS, 363, 1083.
c The frequencies and relative intensities are from Luca Dore's
c calculation with Pickett's code

        dimension jq(64),f1q(64),fq(64),enq(64)
        dimension upper_data(280),lower_data(280)
        dimension frequency(280),relative(280),einstein(280)
        dimension istartline_data(7),mainline(7),frequency_line(7),
     &      einsteina_line(7),einsteina_sum(7)
        dimension tmp1(maxhyp),tmp2(maxhyp),tmp3(maxhyp)

c Data statements for the states and transitions. There is an F77
c limit on 19 continuation statements, so break up the data
c statements into groups of less than 19. Use the double precision
c format with a "D" exponent, even if the exponent is zero, so
c that the back half of the double does not get filled with junk.
c RELATIVE only has 6 significant digits, so it can be single
c to be converted to double later.

        dpdipole = 3.4d-18
        atoms = 29.
        brot = 46586.867d6

c This defines the J quantum numbers for the lower and
c upper rotational states of each of the 7 lines
        do 121 j = 1,nlines
          jloline(j) = j-1
          jupline(j) = j
c        write(22,*) 'line jloline jupline '
c     &    ,j,jloline(j),jupline(j)
 121    continue

c Here are the 3 quantum numbers for each of the 64 states

        data jq /
     &   0, 0, 0, 1, 1, 1, 1, 1,
     &   1, 1, 2, 2, 2, 2, 2, 2,
     &   2, 2, 2, 3, 3, 3, 3, 3,
     &   3, 3, 3, 3, 4, 4, 4, 4,
     &   4, 4, 4, 4, 4, 5, 5, 5,
     &   5, 5, 5, 5, 5, 5, 6, 6,
     &   6, 6, 6, 6, 6, 6, 6, 7,
     &   7, 7, 7, 7, 7, 7, 7, 7 /
        data f1q /
     &   1, 1, 1, 1, 1, 1, 2, 2,
     &   2, 0, 2, 2, 2, 3, 3, 3,
     &   1, 1, 1, 3, 3, 3, 4, 4,
     &   4, 2, 2, 2, 4, 4, 4, 5,
     &   3, 5, 5, 3, 3, 5, 5, 5,
     &   6, 4, 6, 6, 4, 4, 6, 6,
     &   6, 7, 5, 7, 7, 5, 5, 7,
     &   7, 7, 8, 6, 8, 8, 6, 6 /
        data fq /
     &   2, 1, 0, 0, 2, 1, 2, 3,
     &   1, 1, 2, 3, 1, 3, 4, 2,
     &   1, 2, 0, 3, 4, 2, 4, 5,
     &   3, 2, 3, 1, 4, 5, 3, 5,
     &   3, 6, 4, 4, 2, 5, 4, 6,
     &   6, 4, 5, 7, 5, 3, 6, 5,
     &   7, 7, 5, 6, 8, 4, 6, 7,
     &   6, 8, 8, 6, 7, 9, 5, 7 /

c Here is the energy in GHz for each of the 64 states

        data enq /
     &      000.0000D0,    000.0000D0,    000.0000D0,  93171.6167D0,
     &    93171.9134D0,  93172.0484D0,  93173.4755D0,  93173.7723D0,
     &    93173.9626D0,  93176.2608D0, 279516.4477D0, 279516.7030D0,
     &   279516.7694D0, 279518.2351D0, 279518.6326D0, 279518.7426D0,
     &   279519.3252D0, 279519.5369D0, 279519.7891D0, 559028.1346D0,
     &   559028.5123D0, 559028.5682D0, 559030.0476D0, 559030.4957D0,
     &   559030.5533D0, 559030.6708D0, 559031.0470D0, 559031.1749D0,
     &   931700.6052D0, 931701.0433D0, 931701.0673D0, 931702.5885D0,
     &   931702.9983D0, 931703.0686D0, 931703.0866D0, 931703.4423D0,
     &   931703.5023D0,1397525.3911D0,1397525.8591D0,1397525.8665D0,
     &  1397527.4227D0,1397527.7020D0,1397527.9114D0,1397527.9260D0,
     &  1397528.1846D0,1397528.1979D0,1956491.9348D0,1956492.4006D0,
     &  1956492.4368D0,1956494.0034D0,1956494.1900D0,1956494.4822D0,
     &  1956494.5253D0,1956494.6758D0,1956494.6989D0,2608587.5742D0,
     &  2608588.0344D0,2608588.0970D0,2608589.6732D0,2608589.7877D0,
     &  2608590.1420D0,2608590.2109D0,2608590.2628D0,2608590.3168D0 /


c Now we have defined the data for each of the 64 states.
c Below we set up the data for the 280 transitions between
c states.

c Here are the pointers to the upper and lower levels of each
c of the 280 transitions. These are defined in quantum number
c (QN) or C style indexing (starting at 0).

        data (upper_data(i),i=       1,      90) /
     &    3,  4,  4,  5,  5,  5,  6,  6,  7,  8,
     &    8,  8,  9,  9,  9, 10, 12, 15, 10, 10,
     &   12, 11, 10, 16, 11, 17, 12, 18, 10, 13,
     &   10, 12, 13, 15, 11, 12, 14, 15, 12, 15,
     &   16, 17, 17, 18, 16, 17, 13, 15, 15, 16,
     &   16, 17, 17, 16, 18, 19, 21, 21, 19, 19,
     &   21, 20, 19, 20, 21, 24, 25, 25, 27, 22,
     &   19, 26, 27, 19, 21, 20, 24, 22, 27, 23,
     &   21, 24, 25, 21, 26, 24, 26, 27, 25, 26 /
        data (upper_data(i),i=      91,     180) /
     &   22, 24, 25, 25, 24, 25, 26, 27, 26, 27,
     &   28, 30, 28, 28, 30, 30, 29, 28, 29, 30,
     &   32, 34, 28, 31, 36, 32, 35, 32, 36, 28,
     &   30, 29, 34, 31, 30, 33, 34, 36, 35, 30,
     &   35, 36, 32, 34, 35, 31, 32, 32, 34, 32,
     &   35, 36, 34, 35, 36, 37, 37, 37, 38, 38,
     &   39, 37, 38, 38, 39, 41, 37, 40, 42, 41,
     &   45, 41, 44, 45, 37, 38, 38, 39, 42, 40,
     &   42, 43, 44, 45, 41, 44, 45, 38, 42, 44 /
        data (upper_data(i),i=     181,     270) /
     &   40, 41, 41, 42, 41, 45, 44, 42, 44, 45,
     &   46, 46, 46, 47, 47, 48, 46, 47, 47, 48,
     &   50, 46, 49, 50, 51, 53, 50, 53, 54, 47,
     &   47, 46, 51, 48, 51, 49, 52, 53, 50, 54,
     &   54, 53, 47, 51, 54, 49, 50, 50, 51, 50,
     &   53, 54, 51, 53, 54, 55, 55, 55, 56, 56,
     &   55, 57, 56, 56, 57, 59, 55, 58, 59, 60,
     &   62, 62, 59, 56, 60, 63, 56, 55, 60, 57,
     &   58, 61, 62, 59, 63, 63, 62, 56, 60, 63 /
        data (upper_data(i),i=     271,     280) /
     &   58, 59, 59, 60, 59, 62, 63, 60, 62, 63 /
        data (lower_data(i),i=       1,      90) /
     &    1,  1,  0,  2,  1,  0,  1,  0,  0,  1,
     &    0,  2,  2,  1,  0,  9,  9,  9,  8,  7,
     &    8,  7,  6,  9,  6,  9,  6,  9,  5,  7,
     &    4,  5,  6,  8,  4,  4,  7,  7,  3,  6,
     &    8,  8,  7,  8,  6,  6,  4,  5,  4,  5,
     &    4,  5,  4,  3,  5, 17, 17, 16, 15, 14,
     &   15, 14, 13, 13, 13, 17, 17, 16, 18, 14,
     &   11, 17, 17, 10, 12, 11, 15, 13, 16, 14,
     &   11, 14, 15, 10, 15, 13, 14, 15, 13, 13 /
        data (lower_data(i),i=      91,     180) /
     &   11, 11, 12, 11, 10, 10, 11, 12, 10, 10,
     &   26, 26, 24, 23, 25, 24, 23, 22, 22, 22,
     &   26, 26, 20, 23, 27, 25, 26, 24, 26, 19,
     &   21, 20, 24, 22, 20, 23, 23, 25, 24, 19,
     &   23, 24, 22, 22, 22, 20, 21, 20, 20, 19,
     &   20, 21, 19, 19, 19, 35, 34, 33, 35, 34,
     &   33, 31, 32, 31, 31, 35, 29, 33, 35, 34,
     &   36, 32, 35, 35, 28, 30, 29, 29, 34, 31,
     &   33, 33, 34, 34, 31, 33, 32, 28, 31, 31 /
        data (lower_data(i),i=     181,     270) /
     &   29, 30, 29, 29, 28, 30, 29, 28, 28, 28,
     &   44, 43, 42, 44, 42, 43, 40, 41, 40, 40,
     &   44, 39, 43, 42, 44, 45, 41, 44, 44, 39,
     &   38, 37, 43, 39, 42, 40, 43, 42, 40, 43,
     &   42, 41, 37, 40, 40, 39, 39, 38, 39, 37,
     &   38, 39, 37, 37, 37, 54, 52, 51, 54, 51,
     &   49, 52, 50, 49, 49, 54, 48, 52, 51, 54,
     &   54, 53, 50, 48, 52, 54, 47, 46, 51, 48,
     &   49, 52, 51, 49, 52, 51, 50, 46, 49, 49 /
        data (lower_data(i),i=     271,     280) /
     &   48, 48, 47, 48, 46, 47, 48, 46, 46, 46 /

c Here is the energy difference in GHz between the 2 states of
c each transition.

        data (frequency(i),i=       1,      40) /
     &   93171.6086D0, 93171.9054D0, 93171.9054D0, 93172.0403D0,
     &   93172.0403D0, 93172.0403D0, 93173.4675D0, 93173.4675D0,
     &   93173.7643D0, 93173.9546D0, 93173.9546D0, 93173.9547D0,
     &   93176.2527D0, 93176.2527D0, 93176.2527D0,186340.1767D0,
     &  186340.4984D0,186342.4717D0,186342.4747D0,186342.6651D0,
     &  186342.7964D0,186342.9204D0,186342.9618D0,186343.0541D0,
     &  186343.2171D0,186343.2658D0,186343.2835D0,186343.5179D0,
     &  186344.3890D0,186344.4527D0,186344.5240D0,186344.7107D0,
     &  186344.7494D0,186344.7698D0,186344.7793D0,186344.8457D0,
     &  186344.8501D0,186344.9601D0,186345.1424D0,186345.2569D0 /
        data (frequency(i),i=      41,      80) /
     &  186345.3521D0,186345.5638D0,186345.7542D0,186345.8160D0,
     &  186345.8392D0,186346.0509D0,186346.3116D0,186346.6841D0,
     &  186346.8190D0,186347.2664D0,186347.4014D0,186347.4781D0,
     &  186347.6131D0,186347.6981D0,186347.7303D0,279508.5967D0,
     &  279509.0304D0,279509.2421D0,279509.3908D0,279509.5008D0,
     &  279509.8245D0,279509.8785D0,279509.8982D0,279510.2760D0,
     &  279510.3319D0,279511.0156D0,279511.1328D0,279511.3445D0,
     &  279511.3847D0,279511.4140D0,279511.4305D0,279511.5089D0,
     &  279511.6369D0,279511.6858D0,279511.7978D0,279511.8083D0,
     &  279511.8097D0,279511.8114D0,279511.8486D0,279511.8621D0 /
        data (frequency(i),i=      81,     120) /
     &  279511.8642D0,279511.9197D0,279511.9269D0,279512.1195D0,
     &  279512.3029D0,279512.3171D0,279512.4130D0,279512.4309D0,
     &  279512.4343D0,279512.8104D0,279513.3437D0,279513.8494D0,
     &  279513.9002D0,279513.9666D0,279514.1047D0,279514.2219D0,
     &  279514.3427D0,279514.4042D0,279514.5980D0,279514.7260D0,
     &  372669.5840D0,372670.0461D0,372670.0772D0,372670.1349D0,
     &  372670.4221D0,372670.5393D0,372670.5729D0,372670.5829D0,
     &  372671.0209D0,372671.0450D0,372671.9768D0,372672.0656D0,
     &  372672.1184D0,372672.1184D0,372672.3529D0,372672.3529D0,
     &  372672.4208D0,372672.4701D0,372672.4809D0,372672.4961D0 /
        data (frequency(i),i=     121,     160) /
     &  372672.5246D0,372672.5564D0,372672.5589D0,372672.5665D0,
     &  372672.5805D0,372672.5984D0,372672.6165D0,372672.8570D0,
     &  372672.9140D0,372672.9582D0,372672.9717D0,372672.9742D0,
     &  372672.9758D0,372673.0646D0,372673.4197D0,372674.1019D0,
     &  372674.4553D0,372674.5112D0,372674.6001D0,372674.8890D0,
     &  372674.9552D0,372674.9594D0,372674.9778D0,372675.3329D0,
     &  372675.3931D0,465822.0241D0,465822.3792D0,465822.3974D0,
     &  465822.4922D0,465822.8473D0,465822.8727D0,465822.8774D0,
     &  465822.9362D0,465823.3455D0,465823.3527D0,465824.3347D0,
     &  465824.4229D0,465824.4293D0,465824.5447D0,465824.6899D0 /
        data (frequency(i),i=     161,     200) /
     &  465824.7706D0,465824.7787D0,465824.8172D0,465824.8307D0,
     &  465824.8609D0,465824.8669D0,465824.8910D0,465824.8982D0,
     &  465824.8999D0,465824.9092D0,465824.9180D0,465824.9325D0,
     &  465825.1724D0,465825.1859D0,465825.1880D0,465825.1905D0,
     &  465825.2747D0,465825.3290D0,465825.3980D0,465825.6705D0,
     &  465826.4548D0,465826.7094D0,465826.7335D0,465826.9436D0,
     &  465827.1716D0,465827.2054D0,465827.2160D0,465827.3816D0,
     &  465827.6541D0,465827.6675D0,558963.9040D0,558964.1620D0,
     &  558964.1765D0,558964.3698D0,558964.6423D0,558964.6639D0,
     &  558964.6653D0,558964.8524D0,558965.1311D0,558965.1672D0 /
        data (frequency(i),i=     201,     240) /
     &  558966.1589D0,558966.2218D0,558966.2309D0,558966.4314D0,
     &  558966.4518D0,558966.6312D0,558966.6414D0,558966.6447D0,
     &  558966.6677D0,558966.6877D0,558966.6949D0,558966.6971D0,
     &  558966.7098D0,558966.7237D0,558966.7243D0,558966.7342D0,
     &  558966.7527D0,558966.9172D0,558966.9202D0,558966.9258D0,
     &  558966.9402D0,558967.1272D0,558967.1630D0,558967.2131D0,
     &  558967.4290D0,558968.2907D0,558968.4767D0,558968.4840D0,
     &  558968.7696D0,558968.9521D0,558968.9698D0,558968.9855D0,
     &  558969.2449D0,558969.4378D0,558969.4608D0,652093.1420D0,
     &  652093.3151D0,652093.3579D0,652093.6022D0,652093.8182D0 /
        data (frequency(i),i=     241,     280) /
     &  652093.8368D0,652093.8378D0,652094.1110D0,652094.2971D0,
     &  652094.3595D0,652095.3552D0,652095.4039D0,652095.4144D0,
     &  652095.5711D0,652095.7103D0,652095.8303D0,652095.8533D0,
     &  652095.8640D0,652095.8641D0,652095.8833D0,652095.8841D0,
     &  652095.9001D0,652095.9058D0,652095.9262D0,652095.9266D0,
     &  652095.9362D0,652095.9520D0,652096.0462D0,652096.0500D0,
     &  652096.0571D0,652096.1000D0,652096.3391D0,652096.3660D0,
     &  652096.4051D0,652096.5789D0,652097.5032D0,652097.6170D0,
     &  652097.6530D0,652097.9721D0,652098.1189D0,652098.1281D0,
     &  652098.1459D0,652098.4740D0,652098.5940D0,652098.6478D0 /

c Here are the relative intensities for each transition.

        data (relative(i),i=       1,      40) /
     &   3.33338D-01, 2.58367D-01, 1.40832D+00, 5.07969D-01,
     &   1.19788D-01, 3.72252D-01, 1.40830D+00, 2.58364D-01,
     &   2.33333D+00, 6.46597D-01, 3.93832D-02, 3.14016D-01,
     &   1.78015D-01, 2.33608D-01, 5.88356D-01, 8.05468D-05,
     &   5.28403D-04, 4.40354D-03, 3.26578D-02, 7.74858D-02,
     &   3.55803D-01, 7.02922D-01, 1.98482D-01, 6.45816D-01,
     &   8.84652D-02, 1.11384D+00, 5.75724D-02, 2.35350D-01,
     &   1.19437D+00, 2.30417D-01, 4.96930D-01, 2.07372D-01,
     &   2.56685D+00, 1.58628D+00, 2.00862D+00, 1.82270D-02,
     &   3.60000D+00, 1.03631D-02, 5.60499D-01, 3.82761D-01 /
        data (relative(i),i=      41,      80) /
     &   2.94595D-03, 6.56565D-03, 4.54872D-02, 1.57538D-02,
     &   1.70997D-02, 2.21027D-02, 2.73193D-03, 1.22744D-02,
     &   3.91777D-03, 1.87543D-01, 2.40431D-01, 2.49537D-01,
     &   5.62461D-01, 1.06163D-01, 1.48896D-01, 2.77579D-04,
     &   1.05740D-03, 1.85249D-03, 1.31613D-02, 2.69119D-02,
     &   2.63957D-01, 4.48685D-01, 1.91912D-01, 3.70579D-02,
     &   2.14350D-02, 1.92448D-03, 4.15735D-01, 1.29370D+00,
     &   6.00001D-01, 1.94175D-01, 3.78766D-01, 2.55177D+00,
     &   2.92409D-02, 2.38897D+00, 1.64609D+00, 3.37140D+00,
     &   2.71424D+00, 3.66138D+00, 5.04451D-01, 4.71428D+00 /
        data (relative(i),i=      81,     120) /
     &   7.39720D-03, 4.17287D-03, 3.73519D-07, 2.01066D-01,
     &   1.97783D-03, 2.75611D-01, 1.17732D-02, 6.66220D-03,
     &   6.89861D-03, 5.70589D-03, 1.58868D-03, 1.64278D-03,
     &   5.83562D-02, 5.95072D-02, 2.40447D-03, 3.08660D-01,
     &   3.79691D-01, 9.55479D-02, 4.90817D-02, 4.98109D-02,
     &   2.41957D-04, 5.94700D-04, 6.43728D-03, 1.22995D-02,
     &   1.22710D-03, 2.09797D-01, 3.23701D-01, 1.70272D-01,
     &   1.84113D-02, 1.03077D-02, 2.93798D-01, 8.06904D-04,
     &   2.88785D-01, 1.65189D-01, 1.71429D+00, 2.50040D+00,
     &   3.69584D+00, 3.79509D-05, 8.72060D-03, 3.52197D+00 /
        data (relative(i),i=     121,     160) /
     &   2.70843D+00, 4.54678D+00, 3.77985D+00, 4.72281D+00,
     &   3.49312D-03, 5.77778D+00, 2.08296D-03, 3.55513D-01,
     &   7.84206D-04, 1.77265D-01, 4.66546D-03, 3.09111D-03,
     &   3.24951D-03, 2.15650D-01, 2.15777D-03, 8.90688D-04,
     &   2.80384D-02, 2.26323D-02, 7.48814D-04, 2.62951D-01,
     &   2.79528D-01, 1.20676D-01, 8.57661D-04, 1.70248D-02,
     &   1.99343D-02, 1.63631D-04, 3.60030D-03, 6.61015D-03,
     &   3.18910D-04, 1.73850D-01, 2.51026D-01, 1.49478D-01,
     &   7.32918D-04, 5.72875D-03, 1.03387D-02, 2.26749D-01,
     &   2.29909D-01, 1.42915D-01, 3.74587D-04, 3.73661D-05 /
        data (relative(i),i=     161,     200) /
     &   2.77778D+00, 3.61972D+00, 4.76868D+00, 3.71102D-03,
     &   4.61024D+00, 3.75354D+00, 1.88606D-03, 5.64773D+00,
     &   4.82050D+00, 5.76565D+00, 1.18550D-03, 6.81818D+00,
     &   3.56052D-04, 1.65695D-03, 1.75232D-03, 2.30552D-03,
     &   2.68440D-01, 1.54857D-01, 1.77168D-01, 1.00110D-03,
     &   5.30600D-04, 1.49036D-02, 1.08350D-02, 3.73945D-04,
     &   2.16915D-01, 1.20448D-01, 2.19848D-01, 3.97702D-04,
     &   7.80498D-03, 9.78356D-03, 1.10272D-04, 3.95185D-03,
     &   2.20999D-03, 1.74504D-04, 1.48297D-01, 2.04066D-01,
     &   1.32007D-01, 4.53466D-04, 3.50481E-03, 6.34298D-03 /
        data (relative(i),i=     201,     240) /
     &   1.84583D-01, 1.89841D-01, 1.25605D-01, 2.29398D-05,
     &   1.79276D-04, 3.81818D+00, 4.69466D+00, 1.91528D-03,
     &   5.81304D+00, 1.12377D-03, 4.78708D+00, 5.67188D+00,
     &   7.37291D-04, 6.71267D+00, 5.84832D+00, 6.79713D+00,
     &   7.84615D+00, 9.84911D-04, 1.04334D-03, 1.30495D-03,
     &   1.67649D-04, 2.13982D-01, 1.36286D-01, 1.50354D-01,
     &   5.22959D-04, 3.36686D-04, 5.98314D-03, 8.71657D-03,
     &   1.96951D-04, 1.81918D-01, 1.13290D-01, 1.80757D-01,
     &   2.14450D-04, 5.49144D-03, 4.20833D-03, 7.67332D-05,
     &   2.54745D-03, 1.45154D-03, 9.61228D-05, 1.29231D-01 /
        data (relative(i),i=     241,     280) /
     &   1.17688D-01, 1.71451D-01, 2.95136D-04, 2.29715D-03,
     &   4.15718D-03, 1.55638D-01, 1.61198D-01, 1.11883D-01,
     &   1.08806D-05, 7.97371D-05, 1.11642D-03, 4.84615D+00,
     &   5.74549D+00, 7.20317D-04, 4.88792D-04, 6.84299D+00,
     &   5.81282D+00, 6.71704D+00, 6.86860D+00, 7.75773D+00,
     &   7.82123D+00, 8.86667D+00, 6.30814D-04, 6.68468D-04,
     &   8.09969D-04, 7.30132D-05, 1.77293D-01, 1.21204D-01,
     &   1.30595D-01, 2.91680D-04, 2.25403D-04, 3.64011D-03,
     &   5.49610D-03, 1.04914D-04, 1.55724D-01, 1.04757D-01,
     &   1.53308D-01, 1.27579D-04, 3.38099D-03, 2.52343D-03 /

c Here are the locations within the 280 length vector where
c the hyperfine transitions for each deltaJ start. The 15
c hyperfine transitions that make up the J = 1->0 transition
c are located at positions 0-14 in C indexing or 1-15 in F77
c indexing. The 40 hyperfine transitions of J=2->1 are
c at the vector locations between 15-54 (C) or 16-55 (F77),
c and so on.

        data istartline_data/ 0, 15, 55, 100, 145, 190, 235 /

c In addition to the 280 transitions between hyperfine levels
c we also need data on the rotational transitions between the 
c J levels. These are averages over the hyperfine transitions
c belonging to each J transition and do not correspond to the 
c values of any single hyperfine transition.

c However, pick out the strongest hyperfine transition
c for each J level, and mark its location in the 280 length vector

        data mainline / 8, 36, 79, 125, 171, 216, 261 /

c Now we collect some info on the total transition rate
c between J levels. These total rates are the sum of all
c the hyperfine rates between each pair of J levels.

c These frequencies are for the average energy difference
c between J levels, where the average is over all hyperfine
c transitions. The actual numbers are from the JPL database
        data frequency_line /  
     &    93173.7000D6, 186344.7710D6, 279511.7010D6, 
     &   372672.5090D6, 465824.9470D6, 558966.6656D6,
     &   652095.8651D6 /

c These A's for radiative transitions between J levels
c are calculated assuming the molecule is a linear rotor
c with no hyperfine structure.
        data einsteina_line /
     &   3.6280E-05, 3.4817E-04, 1.2586E-03, 3.0925E-03,
     &   6.1752E-03, 1.0827E-02, 1.7375E-02 /    

c -----------------------------------------------------------
c Now we are done with all the data statements for the states
c and the lines. Start processing this data. Begin with the
c data for the states.

c This block checks the quantum numbers in the data statements
c above. This just produces a print statement.
        if (myid .lt. 2) then 
        write(22,*) '------------------------------------'
        write(22,*) 'Check quantum numbers. The sum of the ',
     &    'hyperfine ',
     &    'degeneracies should equal to total degeneracy'
        write(22,*) 'The last number should be 1'
        do 78 jj = 1,8
        sumj = 0.
        do 77 i = 1,64
            if (jq(i) .eq. jj-1) then
               sumj = sumj + 2.*fq(i) + 1.
            endif
 77     continue
        write(22,*) 'Check for level J = ',jj-1,
     &    sumj,sumj / (9.*(2.*(jj-1) + 1))
 78     continue
        write(22,*) '------------------------------------'
        endif

c This block writes the energy difference between each pair of
c states into FREQ(m,n)
        do 5 j = 1,nstate
        do 6 i = 1,nstate
             freq(j,i) = dabs(enq(j) - enq(i))*1.d6
             freq(i,j) = freq(j,i)
 6      continue
 5      continue

c This block copies the quantum numbers from DATA to COMMON
        do 4 i = 1,nstate
           jindx(i) = jq(i)
           kindx(i) = f1q(i)
           lindx(i) = fq(i)
           statdg(i)= 2.d0*dble(fq(i)) + 1.d0
           energy(i) = enq(i)
 4      continue

c Print out the states with the energy in inverse cm to
c compare with the Leiden data.

c        write(22,*) 'Compare with Leiden data '
        do 49 m = 1,64
           dg = 2.d0*dble(fq(m)) + 1.d0
           cminv = enq(m)/c*1.e6
c           write(22,1088) m,cminv,dg,jq(m),f1q(m),fq(m)
 49     continue
 1088   format(i4,2x,f11.6,2x,f5.1,2x,i1,'_',i1,'_',i1)

c ---------------------------------------------------------

c Now we start working on the transitions between states.
c Begin with the data for the rotational transitions between
c J states

c Copy the data from the DATA statements above to COMMON
c Calculate the statistical degeneracy of each of the J levels
c This data will be used in the RT calculations.

        do 9 line = 1,nlines
          istartline(line) = istartline_data(line)
          freqline(line) = frequency_line(line)
          aline(line) = einsteina_line(line)
          guline(line) = 2.d0*dble(line) + 1.d0
          glline(line) = 2.d0*dble(line-1) + 1.d0
c          write(22,*) 'line freqline ',line,freqline(line)
 9      continue

c ---------------------------------------------------------

c Now we start working on the hyperfine transitions.

c Put in the number of hyperfines betweeen each of the J levels

        do 411 line = 1,nlines
          nhyp(line) = 45
          if (line .eq. 1) nhyp(line) = 15
          if (line .eq. 2) nhyp(line) = 40
 411    continue

c Here we copy the data for the upper and lower states of
c each hyperfine transition from DATA to COMMON. 
        do 48 it = 1,280
            upper(it) = upper_data(it)
            lower(it) = lower_data(it)
 48     continue


c Copy the relative intensities from the data statements to
c the array variables that are in COMMON. 

c We want the frequency of the mainline to be the JPL frequency
c for the rotational transition. This will make the velocities
c for the NLTE and LTE hyperfine cases the same.

c This method of counting will cover all the 280 transitions.
c IT is the number of the transition, and if NLINES = 7, 
c runs to 280
c The indexing for IT works in F77 because ihyp starts at 1 and
c istartline starts at 0. So the IT starts at 1.

        do 408 line = 1,nlines
          do 409 ihyp = 1,nhyp(line)
            it = istartline(line) + ihyp
            relint(ihyp,line) = relative(it)
            hypvel(ihyp,line) =
     &        -(frequency(it)*1.d6/frequency_line(line)-1.d0)*c
 409      continue
 408    continue

c Fill in the Einstein A's using formula (1) from 
c Daniel et al 2006, ApJ, 648, 461 
c and the line strengths from Luca Dore's data
c UPPER and LOWER data statements were written by IDL
c so, we need to add 1 to get the indexing to start at
c 1 instead of 0

c This definition is appropriate because we use
c the relative strengths in the line profile
c function

c Here we are calculating the Einstein As. UPPER is the
c the index (0,63) of the upper state of transition IT.
c Add 1 for F77. So the QN F = FQ(UPPER(IT)+1) because
c UPPER is the state number and UPPER+1 is the F77 index
c for the array FQ.
        do 2 it = 1,280 
            einstein(it) = 64.d0*(pi**4)
     &          *((frequency(it)/c*1.d6)**3)
     &          /(3.d0*planck)*(dpdipole**2)*relative(it)
     &          /(2.d0*dble(fq(upper(it)+1))+1.d0)
 2      continue

c Check the Einstein As. Sum the As for all hyperfine
c transitions within a delta J. 
c In the sum, multiply the A by the normalized
c statistical weight of the initial hyperfine state.
c We use the upper state as the initial state.
c Still working with frequency ordering of the DATA statements
        if (myid .lt. 2) then 
        write(22,*) '------------------------------------'
        write(22,*) 'Check the Einstein A of the hyperfines',
     &    'against the total A for deltaJ. Ratio should be 1.'
        write(22,*) 'Sum of A_h    A_DJ        Ratio'
        do 79 line = 1,nlines
           einstj = 0.
        do 76 ihyp = 1,nhyp(line)
            it = istartline(line) + ihyp
            aulhyp(ihyp,line) = einstein(it)
            qju = jq(upper(it)+1)
            qfu = fq(upper(it)+1)
            einstj = einstj + (2.d0*qfu + 1.d0)/(9.d0*(2.d0*qju + 1.d0))
     &               *einstein(it)
c        write(22,1066) line,ihyp,fq(upper(it)+1),jq(upper(it)+1),
c     &            (2.*qfu + 1.)/( 9.*(2.*qju + 1.)),
c     &      einstj,einstein(it)
 1066      format(4i4,1p,3e12.5)
 76     continue
        einsteina_sum(line) = einstj
        write(22,999) 
     &      einsteina_sum(line),einsteina_line(line),
     &      einsteina_sum(line)/einsteina_line(line)


c Try this normalization to make the Pickett As sum to
c the known rotational A
c        do 86 k = k1,k2
c            einstein(k) = 
c     &           einstein(k)*einsteina_line(jj)/einsteina_sum(jj)
c 86     continue
c End of trial normalization

 79     continue
        write(22,*) '------------------------------------'
        endif
 999    format(1P,3E12.3)

c Use the Einstein A sums in the Einstein A DeltaJ for consistency
        do 19 line = 1,nlines
          aline(line) = einsteina_sum(line)
c Next line for the trial normalization
c          aline(line) = einsteina_line(line)
 19     continue


c Copy the Einstein A and the exact hyperfine transition
c frequencies into the 2D arrays A and FREQ that are
c used throughout the code. Because NSTATE can be smaller
c than 64, don't copy beyond the array limits

c        write(22,*) 'line    hyp  it  upp  low freq          A'
        do 3 line = 1,nlines
        do 3 ihyp = 1,nhyp(line)
            it = istartline(line) + ihyp 
            m = upper(it)+1
            n = lower(it)+1
c        write(22,2400) line,ihyp,it,m,n,frequency(it),einstein(it)
 2400   format(5i5,1p2e14.6)
            if (m .le. nstate .and. n .le. nstate) then
                a(m,n) = einstein(it)
                a(n,m) = a(m,n)
                freq(m,n) = frequency(it)*1.d6
                freq(n,m) = freq(m,n)
            endif
 3      continue

c Make sure that freq & einstein A between the same levels are 0
        do 700 j = 1,nstate
           freq(j,j) = 0.d0
           a(j,j)    = 0.d0
 700    continue

c The individual Einstein As can be larger than the sum because
c in the sum, the As are multiplied by the ratio of the
c hyperfine to total statistical weight which is less than 1. 
        if (myid .lt. 2) then 
        write(22,*) 'Check indexed by states'
        write(22,*) 'Energy differences, Einstein As and stat dg for ',
     &   'all downward transitions. Only the rad. allowed transitions ',
     &    'have nonzero Einstein A. All numbers C indexing.'
        write(22,*) ' M  N     Q numbers       frequency GHz  ',
     &          '    Einstein A','         g(final)'
        do 813 m = 2,nstate
        do 812 n = 1,m-1
           write (22,1002) m,n,jindx(m),kindx(m),lindx(m),
     &             jindx(n),kindx(n),lindx(n),freq(m,n),
     &                  a(m,n),int(statdg(n))
 1002    format(2i3,1x,3i2,' -->',3i2,2x,1pe18.10,3x,e16.9,3x,i5)
 812    continue
 813    continue
        endif

c Define INDEXU and INDEXL to be the "main line" in F77 indexing (1-n)
c This is OK as a safety, but we dont want anything to depend on
c this choice of hyperfine states for the rotational transitions, delta J. 
c The reason is because the frequency and Einstein A of the
c the rotational transitions are averages over the hyperfine states
c and not the values for any one state.
c Remember to add one to convert from C style indexing 0:N-1 to F77 1:N
        do 1 line = 1,nlines
            indexu(line) = upper(mainline(line)+1) + 1
            indexl(line) = lower(mainline(line)+1) + 1
c            write(22,*) 'main line ',line,' states ',indexu(line),
c     &        ' --> ',indexl(line)
 1      continue

c This is just a print statement to check that all the above
c has been done properly.
c        write(22,*) 'Check indexed by lines'
c        write(22,1501) 
 1501   format('# line ihyp  quantum numbers',t31,'Rel. Int.',
     &          t45,'Hyp. Vel.',t59,'Einstein A',t73,'Frequency')
        do 208 line = 1,nlines
          do 209 ihyp = 1,nhyp(line)
            ih = istartline(line) + ihyp
            m = upper(ih)+1
            n = lower(ih)+1
 2020   format(3i3,2x,3i2,' ->',3i2,2x,1p3e14.6,1x,e16.9)
            if (myid .lt. 2) write(22,2020),ih,line,ihyp,
     &          jq(m),f1q(m),fq(m),
     &          jq(n),f1q(n),fq(n),
     &          relint(ihyp,line),hypvel(ihyp,line)/1.e5,
     &          a(m,n),freq(m,n)
 209      continue
 208      continue

c ------------------------------------------------------------
c Finally, normalize the line strengths into  relative intensities
c We do not use RELINT for non-LTE hyperfines, so we do not need this.

        do 608 line = 1,nlines
          sumrel = 0.
          do 618 ihyp = 1,nhyp(line)
            sumrel =  sumrel + relint(ihyp,line)
 618      continue
          do 619 ihyp = 1,nhyp(line)
            relint(ihyp,line) = relint(ihyp,line)/sumrel
 619      continue
 608    continue

        status = 1

        return
        end


c23456789112345678921234567893123456789412345678951234567896123456789712
        subroutine n2h_deltaJ_rates

        include 'nlines_f77.h'
        parameter (maxch=2000,maxhyp=50)

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,gco,aline,freqline,guline,glline

        integer u,l

        common /procid/ myid
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /hyperfine/ nltehyp,nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav

        common /n2h_hyp/ temperature(15),rate(15,8,8)

c Data statements for the rates between rotational J levels
c including the deltaJ = 0 rates. These rates came from
c the LAMDA database. The dJ=0 rates are estimated from
c deJong's empirical formula, eqn 17, de Jong, Chu, Dalgarno,
c 1975, ApJ, 199, 69
c The rates are multiplied by 1.e10

c There are rates for 8 levels from J=0 to J=7
c There are no rates for higher levels

c The rates include a zero temperature with zero rate
c to prevent the rates from going negative at temperatures
c below 10 K

c Each data statement below has collision rates for the 15 
c temperatures of each deltaJ transition. The transition is given 
c by the the 2nd and 3rd indices of the array RATE. 
c Only the downward rates and the elastic (final=initial)
c rates are given.

        data temperature /
     &     0.0,  10.0,  20.0,  30.0,  50.0,
     &    70.0, 100.0, 150.0, 200.0, 250.0,
     &   300.0, 350.0, 400.0, 500.0,1000.0 /
        data (rate(i,       1,       1),i=1,15) /
     &    0.00, 11.50, 11.54, 10.26, 10.01,
     &    9.76, 10.30, 11.81, 13.30, 14.72,
     &   16.48, 18.67, 18.55, 18.82, 19.64 /
        data (rate(i,       2,       1),i=1,15) /
     &    0.00,  2.60,  2.30,  2.10,  2.00,
     &    1.90,  1.80,  2.00,  2.20,  2.30,
     &    2.50,  2.70,  2.80,  2.90,  3.20 /
        data (rate(i,       2,       2),i=1,15) /
     &    0.00,  9.92,  9.84,  9.04,  9.01,
     &    8.92,  9.13, 10.92, 12.77, 14.12,
     &   15.98, 18.36, 18.55, 18.77, 19.83 /
        data (rate(i,       3,       1),i=1,15) /
     &    0.00,  1.40,  1.20,  1.10,  1.00,
     &    0.92,  0.88,  0.84,  0.82,  0.81,
     &    0.83,  0.81,  0.85,  0.88,  0.97 /
        data (rate(i,       3,       2),i=1,15) /
     &    0.00,  3.80,  3.70,  3.40,  3.40,
     &    3.30,  3.50,  3.60,  3.80,  4.00,
     &    4.20,  4.40,  4.40,  4.60,  5.10 /
        data (rate(i,       3,       3),i=1,15) /
     &    0.00,  6.15,  7.44,  7.20,  7.86,
     &    8.12,  9.46, 10.61, 11.99, 13.41,
     &   14.70, 16.41, 16.02, 16.40, 17.48 /
        data (rate(i,       4,       1),i=1,15) /
     &    0.00,  1.00,  0.91,  0.84,  0.73,
     &    0.68,  0.63,  0.59,  0.56,  0.55,
     &    0.55,  0.55,  0.55,  0.58,  0.66 /
        data (rate(i,       4,       2),i=1,15) /
     &    0.00,  2.50,  2.30,  2.10,  1.90,
     &    1.80,  1.70,  1.70,  1.70,  1.80,
     &    1.80,  1.80,  1.80,  1.90,  2.10 /
        data (rate(i,       4,       3),i=1,15) /
     &    0.00,  4.30,  4.20,  4.00,  3.90,
     &    3.90,  4.00,  4.20,  4.30,  4.50,
     &    4.70,  4.90,  4.90,  5.10,  5.70 /
        data (rate(i,       4,       4),i=1,15) /
     &    0.00,  4.73,  6.14,  6.38,  7.04,
     &    7.63,  8.72, 10.11, 11.16, 12.46,
     &   13.62, 15.16, 14.82, 15.14, 16.33 /
        data (rate(i,       5,       1),i=1,15) /
     &    0.00,  0.59,  0.58,  0.55,  0.52,
     &    0.46,  0.44,  0.41,  0.41,  0.42,
     &    0.42,  0.43,  0.44,  0.47,  0.55 /
        data (rate(i,       5,       2),i=1,15) /
     &    0.00,  1.60,  1.60,  1.50,  1.50,
     &    1.40,  1.30,  1.20,  1.20,  1.20,
     &    1.20,  1.10,  1.10,  1.20,  1.30 /
        data (rate(i,       5,       3),i=1,15) /
     &    0.00,  2.90,  2.70,  2.50,  2.40,
     &    2.30,  2.20,  2.20,  2.20,  2.20,
     &    2.20,  2.20,  2.20,  2.30,  2.60 /
        data (rate(i,       5,       4),i=1,15) /
     &    0.00,  4.00,  4.00,  3.90,  3.90,
     &    4.00,  4.20,  4.40,  4.50,  4.70,
     &    4.90,  5.10,  5.10,  5.30,  6.00 /
        data (rate(i,       5,       5),i=1,15) /
     &    0.00,  3.39,  4.74,  5.18,  6.04,
     &    6.82,  8.09,  9.47, 10.50, 11.75,
     &   12.86, 14.32, 14.01, 14.32, 15.72 /
        data (rate(i,       6,       1),i=1,15) /
     &    0.00,  0.33,  0.34,  0.34,  0.33,
     &    0.33,  0.31,  0.30,  0.29,  0.27,
     &    0.26,  0.26,  0.26,  0.27,  0.32 /
        data (rate(i,       6,       2),i=1,15) /
     &    0.00,  1.20,  1.20,  1.10,  1.10,
     &    1.00,  0.98,  0.96,  0.96,  0.95,
     &    0.93,  0.94,  0.94,  0.99,  1.20 /
        data (rate(i,       6,       3),i=1,15) /
     &    0.00,  2.10,  2.00,  1.90,  1.80,
     &    1.70,  1.60,  1.50,  1.50,  1.50,
     &    1.40,  1.40,  1.40,  1.50,  1.70 /
        data (rate(i,       6,       4),i=1,15) /
     &    0.00,  3.00,  2.90,  2.80,  2.70,
     &    2.60,  2.50,  2.50,  2.40,  2.40,
     &    2.40,  2.50,  2.50,  2.60,  2.90 /
        data (rate(i,       6,       5),i=1,15) /
     &    0.00,  4.20,  4.10,  4.00,  4.00,
     &    4.00,  4.20,  4.40,  4.60,  4.80,
     &    5.00,  5.10,  5.10,  5.30,  6.10 /
        data (rate(i,       6,       6),i=1,15) /
     &    0.00,  2.92,  4.13,  4.62,  5.52,
     &    6.17,  7.41,  8.77, 10.00, 11.22,
     &   12.30, 13.45, 13.18, 13.50, 15.13 /
        data (rate(i,       7,       1),i=1,15) /
     &    0.00,  0.32,  0.30,  0.29,  0.28,
     &    0.28,  0.28,  0.27,  0.25,  0.24,
     &    0.24,  0.22,  0.22,  0.23,  0.28 /
        data (rate(i,       7,       2),i=1,15) /
     &    0.00,  0.96,  0.93,  0.89,  0.85,
     &    0.81,  0.75,  0.69,  0.66,  0.63,
     &    0.61,  0.59,  0.59,  0.62,  0.72 /
        data (rate(i,       7,       3),i=1,15) /
     &    0.00,  1.70,  1.70,  1.60,  1.50,
     &    1.40,  1.30,  1.30,  1.20,  1.20,
     &    1.20,  1.20,  1.20,  1.30,  1.50 /
        data (rate(i,       7,       4),i=1,15) /
     &    0.00,  2.10,  2.10,  2.00,  1.90,
     &    1.80,  1.70,  1.70,  1.60,  1.60,
     &    1.60,  1.60,  1.60,  1.70,  1.90 /
        data (rate(i,       7,       5),i=1,15) /
     &    0.00,  2.90,  2.80,  2.70,  2.60,
     &    2.50,  2.50,  2.50,  2.50,  2.50,
     &    2.50,  2.50,  2.50,  2.60,  2.90 /
        data (rate(i,       7,       6),i=1,15) /
     &    0.00,  3.80,  3.80,  3.80,  3.90,
     &    4.00,  4.20,  4.40,  4.60,  4.80,
     &    5.00,  5.20,  5.20,  5.40,  6.20 /
        data (rate(i,       7,       7),i=1,15) /
     &    0.00,  2.24,  3.35,  3.91,  4.90,
     &    5.69,  6.91,  8.26,  9.48, 10.67,
     &   11.73, 13.10, 12.86, 13.19, 14.80 /
        data (rate(i,       8,       1),i=1,15) /
     &    0.00,  0.31,  0.30,  0.29,  0.26,
     &    0.24,  0.22,  0.19,  0.17,  0.16,
     &    0.15,  0.14,  0.13,  0.14,  0.16 /
        data (rate(i,       8,       2),i=1,15) /
     &    0.00,  0.82,  0.80,  0.76,  0.72,
     &    0.68,  0.62,  0.57,  0.53,  0.51,
     &    0.48,  0.46,  0.46,  0.49,  0.59 /
        data (rate(i,       8,       3),i=1,15) /
     &    0.00,  1.30,  1.30,  1.20,  1.10,
     &    1.10,  0.98,  0.90,  0.86,  0.82,
     &    0.79,  0.77,  0.77,  0.81,  0.94 /
        data (rate(i,       8,       4),i=1,15) /
     &    0.00,  1.60,  1.60,  1.60,  1.50,
     &    1.50,  1.40,  1.40,  1.30,  1.30,
     &    1.30,  1.30,  1.30,  1.40,  1.60 /
        data (rate(i,       8,       5),i=1,15) /
     &    0.00,  2.00,  2.10,  2.00,  1.90,
     &    1.90,  1.80,  1.70,  1.70,  1.70,
     &    1.70,  1.70,  1.70,  1.80,  2.00 /
        data (rate(i,       8,       6),i=1,15) /
     &    0.00,  2.50,  2.50,  2.50,  2.50,
     &    2.50,  2.50,  2.50,  2.50,  2.50,
     &    2.60,  2.60,  2.60,  2.70,  3.10 /
        data (rate(i,       8,       7),i=1,15) /
     &    0.00,  3.30,  3.50,  3.60,  3.80,
     &    4.00,  4.20,  4.40,  4.70,  4.90,
     &    5.00,  5.20,  5.20,  5.40,  6.20 /
        data (rate(i,       8,       8),i=1,15) /
     &    0.00,  1.70,  2.75,  3.35,  4.41,
     &    5.31,  6.52,  7.87,  9.27, 10.47,
     &   11.30, 12.64, 12.42, 12.77, 14.39 /

c Fill in the upward rates from detailed balance
c These rates refer to the rotational transitions
c of N2H+. The J quantum number of the state is
c U-1 or L-1. The statistical degeneracy of the 
c rotational state is 2J+1. The energy difference
c between the rotational states is calculated off
c the quantum J. Same as in FRQROTOR.
        do 1 u = 2,8
        do 1 l = 1,u
        do 1 i = 2,15
           ju = u - 1
           jl = l - 1
           fq = dble (brot*( ju*(ju + 1)
     &                     - jl*(jl + 1)  ) )
           dekt = planck*fq
     &         / (boltz*dble(temperature(i)))
           rate(i,l,u) = rate(i,u,l)
     &         *float(2*ju+1)/float(2*jl+1)*exp(-dekt)
c        write(22,1001) ju,jl,i,temperature(i),
c     &                dekt,2*ju+1,2*jl+1,rate(i,u,l),rate(i,l,u)
c 1001   format( 'uli T de gu gl rul rlu ',3i5,f8.1,1pe12.4,
c     &     2i5,3e12.4)
 1      continue

c Try setting the elastic rates to zero
        do 2 u=1,8
        do 2 i=1,15
           rate(i,u,u) = 0.*rate(i,u,u)
 2      continue

        if (myid .le. 1) then 
           write(22,*) 'Collision rates for N2H+ hyperfines loaded'
           write(6 ,*) 'Collision rates for N2H+ hyperfines loaded'
           write(22,*) 'Elastic rate for N2H+ hyperfines set to zero'
           write(22,*) 'Elastic rate for N2H+ hyperfines set to zero'
        endif

        return
        end
c23456789112345678921234567893123456789412345678951234567896123456789712

c23456789112345678921234567893123456789412345678951234567896123456789712

       integer function locate(xx,n,x,j)
c XX is the array of data
c N is the length
c X is the point to be located
c J is the index just below X in the array XX
c J is returned and also in the arg list
       dimension xx(n)
c       write(22,*) 'n,x ',n,x
c       do 12 i=1,n
c       write(22,*),i,xx(i)
c 12    continue
       jl = 0
       ju = n+1
 10    if (ju-jl .gt. 1) then
           jm = (ju+jl)/2
           if ((xx(n) .gt. xx(1)) .eqv. (x .gt. xx(jm))) then
               jl = jm
           else
               ju = jm
           endif
           goto 10
       endif
       j = jl

       locate = j
       return
       end


        subroutine colrathcn(tk)

c This is a modified version of the colrateion subroutine

c A linear interpolation of the ion rates.
c Using the tables of temperatures and rates in
c common HCNRATES, interpolates the rates for
c the temperature TK passed as an argument.
c The interpolated rates are stored in the common
c COLL.

c CR                Interpolated collision rates cm^3/sec at
c                        temperature TK
c NSTATE        number of levels
c TEMPS                temperatures in K
c RATE                collision rates at temperatures in TEMPS
c IT                array index

        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline

        double precision energy

        integer u,l
        parameter (ntemp=9)
        include 'nlines_f77.h'
        common /procid/ myid
        common /coll/ cr(nstate,nstate)
        common /hcnrates/ energy(nstate),temps(ntemp),
     &                    rate(ntemp,nstate,nstate)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)

c23456789112345678921234567893123456789412345678951234567896123456789712

c initialize cr as 0
        do 9 u = 1,nstate
        do 9 l = 1,nstate
           cr(u,l) = 0
 9         continue

        it = 0
        do 13 i =1,ntemp
           if (tk .gt. temps(i))        it = it + 1
 13     continue

c This is different from the other collision rate routines
c which do not calculate the elastic rates u=l. But they
c are all zero in the data statements, so it comes to the
c same. Do not know why this routine is different. Do not
c change it because it works.
        do 2 u = 1,nstate
        do 2 l = 1,u

           if (it .eq. 0) then
              cr(u,l) = dble (rate(1,u,l))
           elseif (it .eq. ntemp) then
              cr(u,l) = dble (rate(it,u,l))
           else
c This is a simple linear interpolation.
              cr(u,l) = dble (rate(it,u,l)
     &             + (tk - temps(it))/(temps(it+1) - temps(it))
     &             * (rate(it+1,u,l) - rate(it,u,l)) )
           endif

c Upward rate
           cr(l,u) = cr(u,l)*(statdg(u)/statdg(l))*
     &               exp(-(energy(u)-energy(l))/dble(tk))


        if (cr(u,l) .lt. 0.)
     &          write(6,*) 'NEG RATE ',tk,rate(it,u,l),
     &                        rate(it+1,u,l),cr(u,l)

 2      continue

        return
        end

