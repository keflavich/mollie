        function fraction_para(tk)

        real *8 frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg

        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
 
c rotational constant in cm-1 converted to ergs
        b = 60.853 * c * planck

        j0 = 0
        j1 = 1

        g0 = 2.0*j0 + 1.0
        g1 = 2.0*j1 + 1.0

        e0 = b * j0*(j0 + 1.0)
        e1 = b * j1*(j1 + 1.0)

        p0 = (2.0*j0 + 1.0)*g0 * exp(-e0/(boltz*tk))
        p1 = (2.0*j1 + 1.0)*g1 * exp(-e1/(boltz*tk))

        fraction_para = p0/(p0 + p1)
c        fraction_para = 0.25
        fraction_para = 0.999

c        write(22,*) j0,j1,g0,g1,e0,e1,-e0/(boltz*tk),-e1/(boltz*tk)
c        write(22,*) 'T, p0,p1 fpara',tk,p0,p1,fraction_para

        return
        end

        subroutine colrateh2o(tk,doprint)
 
c A linear interpolation of the collision rates.
 
        character *80 msg
        real *8 frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr
        logical doprint
        integer u
 
        include 'nlines_f77.h'

        common /coll/ cr(nstate,nstate)
        common /h2orates/ temps(9),para_rate(9,8,8),ortho_rate(9,8,8)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        dimension para(nstate,nstate),ortho(nstate,nstate)

c        doprint = .false.
 
        do 67 it = 1,9
              if (temps(it) .gt. tk) go to 68
 67     continue
 68     continue
        it = it - 1
        it = max(1,min(it,8))

        if (doprint) write(22,*) 'tk,it', tk,it
 
c Use the rate at 5 K for the rate at any lower T. Keeps the rates
c from going negative in the interpolation (extrapolation).
        tp = max(tk,5.)
        fp = fraction_para(tp)
        tp = (tp - temps(it))/(temps(it+1) - temps(it)) 

        if (doprint) write(22,*) 'temps it it+1 ',temps(it),temps(it+1)
        if (doprint) write(22,*) 'temperature and o/p factors ',tp,fp

        do 2 u = 1,nstate
        do 2 l = 1,nstate
 
 
c This is a simple linear interpolation.
                para(u,l) = para_rate(it,u,l)
     &          + tp * (para_rate(it+1,u,l) - para_rate(it,u,l)) 
                ortho(u,l) = ortho_rate(it,u,l)
     &          + tp * (ortho_rate(it+1,u,l) - ortho_rate(it,u,l)) 
 

 2      continue

        lun = 23

        if (doprint) then
          msg = ' matrix of para rates'
          do 3 u = 1,nstate
          do 3 l = 1,nstate
            cr(u,l) = dble(para(u,l))
 3        continue
          call dmattmt(cr,msg,nstate,lun) 

          msg = ' matrix of ortho rates'
          do 4 u = 1,nstate
          do 4 l = 1,nstate
            cr(u,l) = dble(ortho(u,l))
 4        continue
          call dmattmt(cr,msg,nstate,lun) 
        endif

        do 5 u = 2,nstate
        do 5 l = 1,u-1
          cr(u,l) = dble(fp*para(u,l) + (1.0-fp)*ortho(u,l))
          cr(l,u) = cr(u,l) * statdg(u)/statdg(l)
     &            * exp( -planck/boltz/tk*freq(u,l) )

 5      continue

        msg = ' matrix of collision rates'
        if (doprint) call dmattmt(cr,msg,nstate,lun) 


 
        return
        end

        subroutine loadh2orates
 
c Copies the molecular ion rates from data
c statements to the common block /H2ORATES/
c Add to the common block rates a zero rate at zero
c temperature so that at temperatures less than 10 K, the
c interpolated rates will not go below zero.


        include 'nlines_f77.h'

c We have rates for 10 states, but there are only 8 states
c set in h2o_init

        real *8 frac,planck,boltz,c,pi,amu,pc,sqrtpi,
     &  grav,a,freq,chfreq,statdg,cr


        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate)
 
        integer u,l
 
        common /h2orates/ temps(9),para_rate(9,8,8),ortho_rate(9,8,8)
        dimension t(12),prates(12,10,10),orates(12,10,10)
 
c23456789112345678921234567893123456789412345678951234567896123456789712

c The 3 indices are in order of temperature, upper state, lower state


c23456789112345678921234567893123456789412345678951234567896123456789712
        data t / 5., 8., 12., 16., 20., 100., 200., 400., 
     &         800., 1200., 1600., 2000. / 

        data (prates(i,2,1),i=1,12) /
     &   3.38E-11, 3.42E-11, 3.40E-11, 3.37E-11, 3.40E-11, 5.40E-11,
     &   8.30E-11, 1.10E-10, 1.40E-10, 1.60E-10, 1.70E-10, 1.80E-10/
        data (prates(i,3,1),i=1,12) /
     &   3.27E-11, 3.26E-11, 3.31E-11, 3.36E-11, 3.40E-11, 5.00E-11,
     &   6.00E-11, 7.20E-11, 8.70E-11, 9.70E-11, 1.10E-10, 1.10E-10/
        data (prates(i,3,2),i=1,12) /
     &   1.19E-11, 1.30E-11, 1.44E-11, 1.54E-11, 1.60E-11, 5.20E-11,
     &   7.70E-11, 9.70E-11, 1.10E-10, 1.10E-10, 1.10E-10, 1.10E-10/
        data (prates(i,4,1),i=1,12) /
     &   2.54E-12, 2.74E-12, 2.97E-12, 3.13E-12, 3.30E-12, 1.60E-12,
     &   6.10E-12, 1.40E-11, 2.30E-11, 2.70E-11, 2.80E-11, 2.90E-11/
        data (prates(i,4,2),i=1,12) /
     &   3.11E-11, 3.19E-11, 3.27E-11, 3.30E-11, 3.30E-11, 4.20E-11,
     &   5.70E-11, 7.50E-11, 9.50E-11, 1.10E-10, 1.20E-10, 1.20E-10/
        data (prates(i,4,3),i=1,12) /
     &   2.15E-11, 2.07E-11, 1.99E-11, 1.93E-11, 1.90E-11, 6.50E-12,
     &   1.80E-11, 3.20E-11, 4.20E-11, 4.40E-11, 4.30E-11, 4.20E-11/
        data (prates(i,5,1),i=1,12) /
     &   6.75E-12, 6.90E-12, 7.08E-12, 7.23E-12, 7.30E-12, 9.70E-12,
     &   1.40E-11, 1.70E-11, 2.00E-11, 2.10E-11, 2.20E-11, 2.20E-11/
        data (prates(i,5,2),i=1,12) /
     &   7.06E-12, 6.77E-12, 6.47E-12, 6.24E-12, 6.10E-12, 5.40E-12,
     &   1.20E-11, 1.80E-11, 2.10E-11, 2.20E-11, 2.20E-11, 2.10E-11/
        data (prates(i,5,3),i=1,12) /
     &   4.00E-11, 3.86E-11, 3.70E-11, 3.57E-11, 3.50E-11, 5.40E-11,
     &   8.70E-11, 1.20E-10, 1.40E-10, 1.40E-10, 1.40E-10, 1.40E-10/
        data (prates(i,5,4),i=1,12) /
     &   5.05E-12, 5.12E-12, 4.92E-12, 4.68E-12, 4.40E-12, 1.60E-12,
     &   5.40E-12, 1.20E-11, 1.90E-11, 2.20E-11, 2.30E-11, 2.30E-11/
        data (prates(i,6,1),i=1,12) /
     &   4.25E-13, 6.80E-13, 1.02E-12, 1.36E-12, 1.70E-12, 1.70E-12,
     &   5.60E-12, 1.10E-11, 1.60E-11, 1.70E-11, 1.70E-11, 1.70E-11/
        data (prates(i,6,2),i=1,12) /
     &   5.00E-13, 8.00E-13, 1.20E-12, 1.60E-12, 2.00E-12, 2.00E-12,
     &   7.40E-12, 1.60E-11, 2.30E-11, 2.50E-11, 2.50E-11, 2.40E-11/
        data (prates(i,6,3),i=1,12) /
     &   8.50E-12, 1.36E-11, 2.04E-11, 2.72E-11, 3.40E-11, 3.40E-11,
     &   4.20E-11, 5.50E-11, 7.50E-11, 9.00E-11, 1.00E-10, 1.10E-10/
        data (prates(i,6,4),i=1,12) /
     &   3.25E-12, 5.20E-12, 7.80E-12, 1.04E-11, 1.30E-11, 1.30E-11,
     &   2.20E-11, 3.20E-11, 4.20E-11, 4.80E-11, 5.20E-11, 5.40E-11/
        data (prates(i,6,5),i=1,12) /
     &   2.75E-11, 4.40E-11, 6.60E-11, 8.80E-11, 1.10E-10, 1.10E-10,
     &   1.50E-10, 1.90E-10, 2.40E-10, 2.70E-10, 2.90E-10, 3.00E-10/
        data (prates(i,7,1),i=1,12) /
     &   4.00E-13, 6.40E-13, 9.60E-13, 1.28E-12, 1.60E-12, 1.60E-12,
     &   3.70E-12, 8.60E-12, 1.70E-11, 2.30E-11, 2.70E-11, 3.10E-11/
        data (prates(i,7,2),i=1,12) /
     &   9.00E-13, 1.44E-12, 2.16E-12, 2.88E-12, 3.60E-12, 3.60E-12,
     &   7.50E-12, 1.30E-11, 1.80E-11, 2.20E-11, 2.40E-11, 2.50E-11/
        data (prates(i,7,3),i=1,12) /
     &   1.98E-12, 3.16E-12, 4.74E-12, 6.32E-12, 7.90E-12, 7.90E-12,
     &   1.40E-11, 1.80E-11, 2.40E-11, 2.90E-11, 3.20E-11, 3.50E-11/
        data (prates(i,7,4),i=1,12) /
     &   1.33E-11, 2.12E-11, 3.18E-11, 4.24E-11, 5.30E-11, 5.30E-11,
     &   6.80E-11, 8.60E-11, 1.10E-10, 1.20E-10, 1.30E-10, 1.40E-10/
        data (prates(i,7,5),i=1,12) /
     &   1.50E-12, 2.40E-12, 3.60E-12, 4.80E-12, 6.00E-12, 6.00E-12,
     &   1.20E-11, 1.90E-11, 2.60E-11, 2.90E-11, 3.10E-11, 3.20E-11/
        data (prates(i,7,6),i=1,12) /
     &   2.22E-12, 3.56E-12, 5.34E-12, 7.12E-12, 8.90E-12, 8.90E-12,
     &   2.00E-11, 3.30E-11, 4.40E-11, 4.90E-11, 5.10E-11, 5.30E-11/
        data (prates(i,8,1),i=1,12) /
     &   7.75E-13, 1.24E-12, 1.86E-12, 2.48E-12, 3.10E-12, 3.10E-12,
     &   5.10E-12, 7.60E-12, 1.10E-11, 1.30E-11, 1.40E-11, 1.60E-11/
        data (prates(i,8,2),i=1,12) /
     &   5.25E-13, 8.40E-13, 1.26E-12, 1.68E-12, 2.10E-12, 2.10E-12,
     &   4.70E-12, 5.60E-12, 7.80E-12, 9.20E-12, 1.00E-11, 1.10E-11/
        data (prates(i,8,3),i=1,12) /
     &   2.75E-12, 4.40E-12, 6.60E-12, 8.80E-12, 1.10E-11, 1.10E-11,
     &   1.50E-11, 2.10E-11, 2.80E-11, 3.30E-11, 3.70E-11, 4.00E-11/
        data (prates(i,8,4),i=1,12) /
     &   5.00E-13, 8.00E-13, 1.20E-12, 1.60E-12, 2.00E-12, 2.00E-12,
     &   3.90E-12, 6.50E-12, 9.00E-12, 1.00E-11, 1.10E-11, 1.20E-11/
        data (prates(i,8,5),i=1,12) /
     &   1.65E-11, 2.64E-11, 3.96E-11, 5.28E-11, 6.60E-11, 6.60E-11,
     &   8.70E-11, 1.10E-10, 1.40E-10, 1.50E-10, 1.70E-10, 1.80E-10/
        data (prates(i,8,6),i=1,12) /
     &   3.00E-12, 4.80E-12, 7.20E-12, 9.60E-12, 1.20E-11, 1.20E-11,
     &   2.20E-11, 3.10E-11, 3.60E-11, 3.60E-11, 3.60E-11, 3.50E-11/
        data (prates(i,8,7),i=1,12) /
     &   9.50E-13, 1.52E-12, 2.28E-12, 3.04E-12, 3.80E-12, 3.80E-12,
     &   9.10E-12, 1.60E-11, 2.30E-11, 2.60E-11, 2.80E-11, 2.80E-11/
        data (prates(i,9,1),i=1,12) /
     &   7.50E-15, 1.20E-14, 1.80E-14, 2.40E-14, 3.00E-14, 3.00E-14,
     &   8.20E-14, 2.90E-13, 8.60E-13, 1.30E-12, 1.60E-12, 1.90E-12/
        data (prates(i,9,2),i=1,12) /
     &   1.95E-13, 3.12E-13, 4.68E-13, 6.24E-13, 7.80E-13, 7.80E-13,
     &   2.10E-12, 4.30E-12, 8.90E-12, 1.30E-11, 1.70E-11, 2.10E-11/
        data (prates(i,9,3),i=1,12) /
     &   1.82E-13, 2.92E-13, 4.38E-13, 5.84E-13, 7.30E-13, 7.30E-13,
     &   1.70E-12, 3.10E-12, 6.20E-12, 8.90E-12, 1.10E-11, 1.40E-11/
        data (prates(i,9,4),i=1,12) /
     &   1.02E-11, 1.64E-11, 2.46E-11, 3.28E-11, 4.10E-11, 4.10E-11,
     &   5.70E-11, 7.60E-11, 9.90E-11, 1.10E-10, 1.20E-10, 1.30E-10/
        data (prates(i,9,5),i=1,12) /
     &   2.30E-13, 3.68E-13, 5.52E-13, 7.36E-13, 9.20E-13, 9.20E-13,
     &   1.90E-12, 3.10E-12, 6.60E-12, 1.00E-11, 1.30E-11, 1.60E-11/
        data (prates(i,9,6),i=1,12) /
     &   3.75E-13, 6.00E-13, 9.00E-13, 1.20E-12, 1.50E-12, 1.50E-12,
     &   3.20E-12, 4.80E-12, 8.70E-12, 1.20E-11, 1.50E-11, 1.70E-11/
        data (prates(i,9,7),i=1,12) /
     &   9.75E-12, 1.56E-11, 2.34E-11, 3.12E-11, 3.90E-11, 3.90E-11,
     &   5.60E-11, 7.40E-11, 9.30E-11, 1.00E-10, 1.10E-10, 1.20E-10/
        data (prates(i,9,8),i=1,12) /
     &   4.00E-13, 6.40E-13, 9.60E-13, 1.28E-12, 1.60E-12, 1.60E-12,
     &   3.70E-12, 3.60E-12, 6.20E-12, 8.30E-12, 1.00E-11, 1.20E-11/
        data (prates(i,10,1),i=1,12) /
     &   7.75E-14, 1.24E-13, 1.86E-13, 2.48E-13, 3.10E-13, 3.10E-13,
     &   7.40E-13, 1.60E-12, 2.70E-12, 3.30E-12, 3.70E-12, 4.00E-12/
        data (prates(i,10,2),i=1,12) /
     &   1.70E-13, 2.72E-13, 4.08E-13, 5.44E-13, 6.80E-13, 6.80E-13,
     &   1.40E-12, 2.80E-12, 5.70E-12, 8.00E-12, 9.90E-12, 1.10E-11/
        data (prates(i,10,3),i=1,12) /
     &   7.25E-14, 1.16E-13, 1.74E-13, 2.32E-13, 2.90E-13, 2.90E-13,
     &   5.10E-13, 6.60E-12, 1.10E-11, 1.40E-11, 1.60E-11, 1.80E-11/
        data (prates(i,10,4),i=1,12) /
     &   1.40E-13, 2.24E-13, 3.36E-13, 4.48E-13, 5.60E-13, 5.60E-13,
     &   1.50E-12, 1.00E-11, 1.60E-11, 2.10E-11, 2.40E-11, 2.80E-11/
        data (prates(i,10,5),i=1,12) /
     &   1.02E-12, 1.64E-12, 2.46E-12, 3.28E-12, 4.10E-12, 4.10E-12,
     &   1.20E-11, 2.10E-11, 2.90E-11, 3.10E-11, 3.10E-11, 3.00E-11/
        data (prates(i,10,6),i=1,12) /
     &   4.00E-12, 6.40E-12, 9.60E-12, 1.28E-11, 1.60E-11, 1.60E-11,
     &   2.60E-11, 4.10E-11, 6.30E-11, 8.10E-11, 9.50E-11, 1.10E-10/
        data (prates(i,10,7),i=1,12) /
     &   7.25E-12, 1.16E-11, 1.74E-11, 2.32E-11, 2.90E-11, 2.90E-11,
     &   5.10E-11, 7.80E-11, 1.10E-10, 1.30E-10, 1.40E-10, 1.50E-10/
        data (prates(i,10,8),i=1,12) /
     &   1.77E-12, 2.84E-12, 4.26E-12, 5.68E-12, 7.10E-12, 7.10E-12,
     &   2.60E-11, 4.80E-11, 5.40E-11, 4.80E-11, 4.20E-11, 3.60E-11/
        data (prates(i,10,9),i=1,12) /
     &   1.17E-12, 1.88E-12, 2.82E-12, 3.76E-12, 4.70E-12, 4.70E-12,
     &   1.20E-11, 2.00E-11, 2.60E-11, 2.60E-11, 2.60E-11, 2.50E-11/


        data (orates(i,2,1),i=1,12) /
     &   2.35E-10, 2.51E-10, 2.67E-10, 2.80E-10, 2.90E-10, 1.20E-10,
     &   1.30E-10, 1.40E-10, 1.60E-10, 1.60E-10, 1.60E-10, 1.70E-10/
        data (orates(i,3,1),i=1,12) /
     &   6.21E-11, 7.96E-11, 9.32E-11, 1.02E-10, 1.10E-10, 6.40E-11,
     &   7.50E-11, 8.70E-11, 1.00E-10, 1.10E-10, 1.20E-10, 1.20E-10/
        data (orates(i,3,2),i=1,12) /
     &   6.54E-11, 8.38E-11, 9.71E-11, 1.05E-10, 1.10E-10, 1.10E-10,
     &   1.10E-10, 1.10E-10, 1.10E-10, 1.20E-10, 1.20E-10, 1.20E-10/
        data (orates(i,4,1),i=1,12) /
     &   1.64E-11, 2.14E-11, 2.47E-11, 2.66E-11, 2.80E-11, 1.30E-11,
     &   1.40E-11, 1.70E-11, 2.20E-11, 2.60E-11, 3.00E-11, 3.30E-11/
        data (orates(i,4,2),i=1,12) /
     &   3.84E-11, 5.07E-11, 5.89E-11, 6.36E-11, 6.70E-11, 6.20E-11,
     &   7.10E-11, 8.40E-11, 1.00E-10, 1.20E-10, 1.30E-10, 1.30E-10/
        data (orates(i,4,3),i=1,12) /
     &   4.62E-11, 6.37E-11, 7.72E-11, 8.54E-11, 9.10E-11, 3.20E-11,
     &   3.30E-11, 3.60E-11, 3.90E-11, 4.10E-11, 4.30E-11, 4.40E-11/
        data (orates(i,5,1),i=1,12) /
     &   1.66E-11, 1.87E-11, 2.01E-11, 2.08E-11, 2.10E-11, 1.50E-11,
     &   1.80E-11, 2.10E-11, 2.10E-11, 2.10E-11, 2.10E-11, 2.10E-11/
        data (orates(i,5,2),i=1,12) /
     &   2.15E-11, 2.33E-11, 2.40E-11, 2.41E-11, 2.40E-11, 1.40E-11,
     &   1.50E-11, 1.80E-11, 2.30E-11, 2.70E-11, 3.00E-11, 3.30E-11/
        data (orates(i,5,3),i=1,12) /
     &   7.33E-11, 8.47E-11, 9.35E-11, 9.93E-11, 1.00E-10, 1.00E-10,
     &   1.20E-10, 1.40E-10, 1.40E-10, 1.50E-10, 1.50E-10, 1.50E-10/
        data (orates(i,5,4),i=1,12) /
     &   3.16E-11, 3.67E-11, 4.11E-11, 4.37E-11, 4.50E-11, 2.20E-11,
     &   1.70E-11, 1.70E-11, 2.10E-11, 2.40E-11, 2.70E-11, 3.00E-11/
        data (orates(i,6,1),i=1,12) /
     &   1.50E-12, 2.40E-12, 3.60E-12, 4.80E-12, 6.00E-12, 6.00E-12,
     &   8.00E-12, 1.00E-11, 1.30E-11, 1.40E-11, 1.50E-11, 1.60E-11/
        data (orates(i,6,2),i=1,12) /
     &   3.50E-12, 5.60E-12, 8.40E-12, 1.12E-11, 1.40E-11, 1.40E-11,
     &   1.50E-11, 1.80E-11, 2.30E-11, 2.70E-11, 3.00E-11, 3.30E-11/
        data (orates(i,6,3),i=1,12) /
     &   1.27E-11, 2.04E-11, 3.06E-11, 4.08E-11, 5.10E-11, 5.10E-11,
     &   5.70E-11, 6.90E-11, 8.70E-11, 1.00E-10, 1.10E-10, 1.20E-10/
        data (orates(i,6,4),i=1,12) /
     &   1.33E-11, 2.12E-11, 3.18E-11, 4.24E-11, 5.30E-11, 5.30E-11,
     &   5.10E-11, 5.10E-11, 5.10E-11, 5.20E-11, 5.20E-11, 5.20E-11/
        data (orates(i,6,5),i=1,12) /
     &   3.75E-11, 6.00E-11, 9.00E-11, 1.20E-10, 1.50E-10, 1.50E-10,
     &   1.90E-10, 2.40E-10, 2.90E-10, 3.20E-10, 3.50E-10, 3.60E-10/
        data (orates(i,7,1),i=1,12) /
     &   6.75E-13, 1.08E-12, 1.62E-12, 2.16E-12, 2.70E-12, 2.70E-12,
     &   2.70E-12, 3.60E-12, 5.90E-12, 8.30E-12, 1.10E-11, 1.30E-11/
        data (orates(i,7,2),i=1,12) /
     &   1.50E-12, 2.40E-12, 3.60E-12, 4.80E-12, 6.00E-12, 6.00E-12,
     &   1.10E-11, 1.60E-11, 2.00E-11, 2.20E-11, 2.40E-11, 2.40E-11/
        data (orates(i,7,3),i=1,12) /
     &   3.00E-12, 4.80E-12, 7.20E-12, 9.60E-12, 1.20E-11, 1.20E-11,
     &   1.60E-11, 2.00E-11, 2.50E-11, 2.80E-11, 3.00E-11, 3.20E-11/
        data (orates(i,7,4),i=1,12) /
     &   2.05E-11, 3.28E-11, 4.92E-11, 6.56E-11, 8.20E-11, 8.20E-11,
     &   9.70E-11, 1.10E-10, 1.30E-10, 1.30E-10, 1.40E-10, 1.40E-10/
        data (orates(i,7,5),i=1,12) /
     &   6.25E-12, 1.00E-11, 1.50E-11, 2.00E-11, 2.50E-11, 2.50E-11,
     &   2.20E-11, 2.30E-11, 2.60E-11, 2.90E-11, 3.20E-11, 3.40E-11/
        data (orates(i,7,6),i=1,12) /
     &   9.50E-12, 1.52E-11, 2.28E-11, 3.04E-11, 3.80E-11, 3.80E-11,
     &   4.10E-11, 4.30E-11, 4.50E-11, 4.60E-11, 4.70E-11, 4.80E-11/
        data (orates(i,8,1),i=1,12) /
     &   7.75E-13, 1.24E-12, 1.86E-12, 2.48E-12, 3.10E-12, 3.10E-12,
     &   5.10E-12, 7.60E-12, 1.10E-11, 1.30E-11, 1.40E-11, 1.60E-11/
        data (orates(i,8,2),i=1,12) /
     &   1.00E-12, 1.60E-12, 2.40E-12, 3.20E-12, 4.00E-12, 4.00E-12,
     &   4.00E-12, 5.00E-12, 7.00E-12, 9.00E-12, 1.10E-11, 1.30E-11/
        data (orates(i,8,3),i=1,12) /
     &   3.25E-12, 5.20E-12, 7.80E-12, 1.04E-11, 1.30E-11, 1.30E-11,
     &   1.80E-11, 2.40E-11, 3.00E-11, 3.40E-11, 3.70E-11, 3.90E-11/
        data (orates(i,8,4),i=1,12) /
     &   1.33E-12, 2.12E-12, 3.18E-12, 4.24E-12, 5.30E-12, 5.30E-12,
     &   4.70E-12, 5.70E-12, 8.80E-12, 1.20E-11, 1.60E-11, 2.00E-11/
        data (orates(i,8,5),i=1,12) /
     &   2.23E-11, 3.56E-11, 5.34E-11, 7.12E-11, 8.90E-11, 8.90E-11,
     &   1.20E-10, 1.50E-10, 1.70E-10, 1.80E-10, 1.80E-10, 1.80E-10/
        data (orates(i,8,6),i=1,12) /
     &   8.50E-12, 1.36E-11, 2.04E-11, 2.72E-11, 3.40E-11, 3.40E-11,
     &   3.20E-11, 3.30E-11, 3.80E-11, 4.20E-11, 4.50E-11, 4.80E-11/
        data (orates(i,8,7),i=1,12) /
     &   6.25E-12, 1.00E-11, 1.50E-11, 2.00E-11, 2.50E-11, 2.50E-11,
     &   1.80E-11, 1.80E-11, 2.20E-11, 2.60E-11, 3.10E-11, 3.60E-11/
        data (orates(i,9,1),i=1,12) /
     &   7.50E-15, 1.20E-14, 1.80E-14, 2.40E-14, 3.00E-14, 3.00E-14,
     &   8.20E-14, 2.90E-13, 8.60E-13, 1.30E-12, 1.60E-12, 1.90E-12/
        data (orates(i,9,2),i=1,12) /
     &   9.00E-13, 1.44E-12, 2.16E-12, 2.88E-12, 3.60E-12, 3.60E-12,
     &   8.20E-12, 1.40E-11, 1.90E-11, 2.10E-11, 2.20E-11, 2.30E-11/
        data (orates(i,9,3),i=1,12) /
     &   5.75E-13, 9.20E-13, 1.38E-12, 1.84E-12, 2.30E-12, 2.30E-12,
     &   3.30E-12, 4.90E-12, 7.30E-12, 9.10E-12, 1.10E-11, 1.20E-11/
        data (orates(i,9,4),i=1,12) /
     &   1.75E-11, 2.80E-11, 4.20E-11, 5.60E-11, 7.00E-11, 7.00E-11,
     &   8.50E-11, 1.00E-10, 1.20E-10, 1.30E-10, 1.40E-10, 1.50E-10/
        data (orates(i,9,5),i=1,12) /
     &   1.60E-12, 2.56E-12, 3.84E-12, 5.12E-12, 6.40E-12, 6.40E-12,
     &   3.70E-12, 3.80E-12, 5.80E-12, 8.60E-12, 1.20E-11, 1.60E-11/
        data (orates(i,9,6),i=1,12) /
     &   2.20E-12, 3.52E-12, 5.28E-12, 7.04E-12, 8.80E-12, 8.80E-12,
     &   6.60E-12, 6.90E-12, 9.30E-12, 1.20E-11, 1.50E-11, 1.80E-11/
        data (orates(i,9,7),i=1,12) /
     &   2.42E-11, 3.88E-11, 5.82E-11, 7.76E-11, 9.70E-11, 9.70E-11,
     &   1.10E-10, 1.10E-10, 1.20E-10, 1.30E-10, 1.30E-10, 1.30E-10/
        data (orates(i,9,8),i=1,12) /
     &   8.00E-13, 1.28E-12, 1.92E-12, 2.56E-12, 3.20E-12, 3.20E-12,
     &   3.80E-12, 5.30E-12, 8.00E-12, 1.00E-11, 1.30E-11, 1.50E-11/
        data (orates(i,10,1),i=1,12) /
     &   7.75E-14, 1.24E-13, 1.86E-13, 2.48E-13, 3.10E-13, 3.10E-13,
     &   7.40E-13, 1.60E-12, 2.70E-12, 3.30E-12, 3.70E-12, 4.00E-12/
        data (orates(i,10,2),i=1,12) /
     &   1.70E-13, 2.72E-13, 4.08E-13, 5.44E-13, 6.80E-13, 6.80E-13,
     &   1.40E-12, 1.80E-12, 4.00E-12, 6.30E-12, 8.40E-12, 1.00E-11/
        data (orates(i,10,3),i=1,12) /
     &   1.40E-12, 2.24E-12, 3.36E-12, 4.48E-12, 5.60E-12, 5.60E-12,
     &   5.90E-12, 7.70E-12, 1.20E-11, 1.60E-11, 2.00E-11, 2.40E-11/
        data (orates(i,10,4),i=1,12) /
     &   2.75E-12, 4.40E-12, 6.60E-12, 8.80E-12, 1.10E-11, 1.10E-11,
     &   1.70E-11, 2.10E-11, 2.30E-11, 2.30E-11, 2.30E-11, 2.30E-11/
        data (orates(i,10,5),i=1,12) /
     &   8.75E-12, 1.40E-11, 2.10E-11, 2.80E-11, 3.50E-11, 3.50E-11,
     &   3.80E-11, 4.20E-11, 4.40E-11, 4.60E-11, 4.70E-11, 4.80E-11/
        data (orates(i,10,6),i=1,12) /
     &   9.00E-12, 1.44E-11, 2.16E-11, 2.88E-11, 3.60E-11, 3.60E-11,
     &   4.60E-11, 6.00E-11, 7.90E-11, 9.30E-11, 1.10E-10, 1.20E-10/
        data (orates(i,10,7),i=1,12) /
     &   2.75E-11, 4.40E-11, 6.60E-11, 8.80E-11, 1.10E-10, 1.10E-10,
     &   1.20E-10, 1.30E-10, 1.30E-10, 1.40E-10, 1.40E-10, 1.40E-10/
        data (orates(i,10,8),i=1,12) /
     &   1.62E-11, 2.60E-11, 3.90E-11, 5.20E-11, 6.50E-11, 6.50E-11,
     &   7.30E-11, 7.10E-11, 6.30E-11, 5.60E-11, 5.10E-11, 4.80E-11/
        data (orates(i,10,9),i=1,12) /
     &   7.25E-12, 1.16E-11, 1.74E-11, 2.32E-11, 2.90E-11, 2.90E-11,
     &   3.40E-11, 3.50E-11, 3.40E-11, 3.30E-11, 3.10E-11, 3.00E-11/

c ---------------------------------------------------------------
c In the LAMDB there are low temperature rates for the first 5
c levels, but no rates for T less than 20 for higher levels. In
c creating the tables above, the low temperature rates for the
c higher levels were calculated by interpolation to a zero rate
c at zero temperature. Instead, to use the calculated rate at
c 20 K for lower temperatures (for levels 6,7,8), overwrite the
c interpolated rates at lower T with the rates at 20 K.
       do 4 i = 1,4
       do 4 u = 6,8
       do 4 l = 6,8

           orates(i,u,l) = orates(5,u,l)
           prates(i,u,l) = prates(5,u,l)

 4     continue
c ---------------------------------------------------------------


c We have rates for 12 temperatures and 10 states. The original code
c used 9 temperatures and 8 states. The arrays ORTHO_RATE and PARA_RATE
c are the ones that are passed through COMMON H2ORATES. Leave these
c at the original size and copy the subset of PRATES and ORATES into
c these smaller arrays.

c Compute upward rates using detailed balance eqn.

       do 2 i = 1,9
       do 2 u = 2,nstate
       do 2 l = 1,u-1
          ortho_rate(i,u,l) = orates(i,u,l)
          para_rate (i,u,l) = prates(i,u,l)

          ortho_rate(i,l,u) = orates(i,u,l) *
     &       statdg(u)/statdg(l) *
     &       exp( -planck/boltz/t(i)*freq(u,l)  )

           para_rate(i,l,u) = prates(i,u,l) *
     &       statdg(u)/statdg(l) *
     &       exp( -planck/boltz/t(i)*freq(u,l) )


 2     continue

       do 3 i = 1,9
         temps(i) = t(i)
 3     continue



       return
       end
