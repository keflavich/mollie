        subroutine ch3cnint

        include 'nlines_f77.h'
        parameter (maxtrans=500)
        implicit real *8 (a-h,o-z)
        real *4 fnh3,atoms,brot,dipole
        integer parity,upper,lower
        logical docont,doprint,doprint2,nh3tau
        common /procid/ myid
        common /ri/ ric(nstate,nstate,9),fsc(nstate,nstate,9)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /frac/ frac(nstate),xfrac(nstate),gg(nstate)
        common /index/ jindx(nstate),kindx(nstate),parity(nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /lines/ indexu(nlines),indexl(nlines)
        common /abc / arot11,brot11,crot11,rc,inuc
        common /frqinv/ frqinv(nstate),sgn(nstate)
        common /nltehyp/ upper(maxtrans),lower(maxtrans),
     &          istartline(nlines)


        dimension v0(nstate,nstate),w(nstate)

        atoms = 41.0
        arot = 158.0990d9
        brot = 9198.899d6
        crot = brot
        arot11 = arot/1.d11
        brot11 = brot/1.d11
        crot11 = brot/1.d11
        rc = 4.22534d6
        dipole = 3.92197d-18
        inuc = 1
        dj  = 3.8048e3
        djk = 177.417e3
        hjjj = 0.0140
        hjjk = 1.071
        hjkk = 6.006
        planck = 6.627d-27
        boltz = 1.38d-16

c FRQINV and SGN are not used with the CH3CN molecules. 
c Set these so that we can use the LTE function written for NH3.
        do 16 m = 1,nstate
            frqinv(m) = 0.d0
            sgn(m)    = 0.d0
 16     continue

c Build up the J,K states until there are NSTATES of them.
        if (myid.lt.1)write(22,*)
     &     '   #    J    K    parity  energy level'
        i = 0
        do 7 jj = 1,20
        do 7 kk = 1,jj
                j = jj - 1
c                k = kk - 1 ; Reverse the order so that
c the frequencies increase in the order of the lines?
                k = jj-kk
                i = i + 1
                jindx(i) = j
                kindx(i) = k
                parity(i)= 0 
c Energy levels do not depend on K because C=B
c W = B *J*(J+1) - (C-B)*K^2
                w(i) = brot*j*(j+1)  + (arot-brot)*k*k
        if (myid.lt.2) write(22,1001) i,j,k,parity(i) , w(i)/1.e9
                if (i .eq. nstate) go to 8
 7      continue
 8      continue
 1001   format (4i5,f16.6)

c Selection rules for line transitions. Go through the 
c combinations of states and select those transitions 
c that follow the selection rules
        if (myid.lt.2)write(22,*)'Jinit Jfinal Kfinal     Freq',
     &     '            Energy Difference'
        li = 0
        do 9 m = 1,nstate
        do 9 n = 1,m-1
        ji = jindx(m)
        jf = jindx(n)
        ki = kindx(m)
        kf = kindx(n)
c Fill in the energy difference for each pair of states even if there
c is no line.
       freq(m,n) = w(m) - w(n)
c       write(22,*) 'line jk',li,ji,ki,' -?-> ',jf,kf
         if ( (ji-jf) .eq. 1 
     &          .and. (ki-kf) .eq. 0) then
c li is the running number of the line. With the ordering of
c KINDX above, the lines will increase in frequency with LI
                li = li + 1
                indexu(li) = m
                indexl(li) = n
c        write(22,*)'new li  ',li,':',indexu(li),'->',indexl(li)
c        write(22,*)'transition',li,':',
c     &     jindx(indexu(li)),kindx(indexu(li)),
c     &          '->',jindx(indexl(li)),kindx(indexl(li))
        j = jf
        k = kf
c If there is a line, use the more accurate equation for the line freq.
        freq(m,n) = 2.*brot*(j+1) - 4.*dj*(j+1)**3 - 2.*djk*(J+1)*k**2
     &   +(hjjj*(j+1)**3)*((j+2)**3 - j**3) + (4.*hjjk*(J+1)**3)*(k**2)
     &    + 2.*hjjk*(j+1)*k**4
            if (myid .lt. 2) 
     &  write(22,1002)ji,jf,kf,freq(m,n)/1.e9,(w(m)-w(n))/1.e9
 1002   format(3i6,2f16.6)
            endif
            if (li .eq. nlines) go to 10
 9      continue
 10     continue

            if (myid .lt. 2) then
        write(22,*) 'These are the Energy levels of the ',
     &      ' rotational states of CH3CN '
        write(22,*) '   #    J    K             GHz',
     &      '          cm^-1               K'
        do 11 i = 1,nstate
          write(22,1003) i,jindx(i),kindx(i),w(i)/1.d9,w(i)/3.d10,
     &     w(i)*planck/boltz
 1003   format(3i5,3f16.6)

 11     continue
        write(22,*) 'These are the line transitions. Freq in GHz'
        do 12 li = 1,nlines
        write(22,1004)li,
     &     jindx(indexu(li)),kindx(indexu(li)),
     &          jindx(indexl(li)),kindx(indexl(li)),
     &       freq(indexu(li),indexl(li))/1.d9
 1004   format(i5,5x,i3,',',i3,' ->',i3,',',i3,f16.6)

 12     continue
            endif

c Next step. Hyperfine shifting and relative intensities


        call fsri
        call normal
        call comprs

        do 14 m = 1,nstate
          gg(m) = gch3cn(m)
          statdg(m) = gg(m)
 14     continue
        call einsta

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


c LTE population test. Store the population in W, no longer used

            if (myid .lt. 2) then
        write(22,*) 'Here are the level pops *1000 for 100 K'
        j = ltech3cn(50.,w)
        do 363 i = 1,nstate
            write(22,*) i,' J K =',jindx(i),kindx(i),w(i)*1000.
 363    continue

c        write(22,*) 'Check for population inversion '
c        do 364 line = 1,nlines
c            m = indexu(line)
c            n = indexl(line)
c            write(22,*),jindx(m),kindx(m),jindx(n),kindx(n),
c     &        int(statdg(m)),int(statdg(n)),w(m),w(n),
c     &        w(n)*statdg(m) - w(m)*statdg(n)
c 364    continue

c This is not used
c        do 31 line = 1,nlines
c        m = indexu(line)
c        n = indexl(line)
c           fact = ( float(jindx(m)**2) - float(kindx(m)**2))
c     &        / float(jindx(m))
c        do 32 l = 1,9
c          ric(m,n,l) = ric(m,n,l)*fact
c 32    continue
c 31    continue


        do 30 line = 1,nlines
        m = indexu(line)
        n = indexl(line)

        do 77 lunsel = 2,2
        lun = 6
        if (lunsel .eq. 2) lun = 22
        write(lun,*) ' '
        write(lun,*) '------------------------------'
        write(lun,*) ' '
        write(lun,*) ' CH3CN '
        write(lun,*) 'transition ',jindx(n),kindx(n),' <-- ',
     &     jindx(m),kindx(m)
        write(lun,*) 'frequency ',freq(m,n)
        write(lun,*) 'einstein a ',a(m,n)
        write(lun,*) 'statistical degeneracies ',gg(n),' <-- ',gg(m)
        write(lun,*) 'relative intensities'
        write(lun,1200) (ric(m,n,l),l=1,9)
        write(lun,*) 'frequency shifts'
        write(lun,1200) (fsc(m,n,l),l=1,9)
 1200   format (1x,1pe12.4,2e12.4,/,3e12.4,/,3e12.4)
        write(lun,*) ' '
            fact = ( float(jindx(m)**2) - float(kindx(m)**2))
     &         / float(jindx(m))
        write(22,*) 'factor ',fact
        l = 0
        do 78 li=1,3
        do 78 lf=1,3
           l = l+1
           if (ric(m,n,l) .gt. 0) write(22,*) ,
     &       (freq(m,n)+fsc(m,n,l))/1.e6, ric(m,n,l)*fact
 78     continue
        write(lun,*) ' '
        write(lun,*) '------------------------------'
        write(lun,*) ' '
 77     continue

 30     continue
            endif


        end

        function g2ch3cn(j,k)

c Computes the statistical degeneracy of a J,K level.
c Identical to GCH3CN expect that here there are 2 ARGS J,K
c This is used in the partition function so that we can sum
c more levels than are defined by NSTATE

c The formula for the statistical deg of J levels is  2J+1.
c The Leiden data base has twice these values for both para and ortho
c so I included another factor of 2.

        include 'nlines_f77.h'
        implicit real *8 (a-h,o-z)
        integer parity
        common /index/ jindx(nstate),kindx(nstate),parity(nstate)
        spin=0.5d0

c First compute the degeneracy due to spin. This depends on K.
c The formulae are 3-46 of Townes & Schawlow
        if (mod(k,3) .eq. 0 ) then
                if (k .ne. 0) then
                    s = 2.*(4.d0*spin**2 + 4.d0*spin + 3.d0)
                else
                    s = (4.d0*spin**2 + 4.d0*spin + 3.d0)
                endif
        else
                s = 2.*(4.d0*spin**2 + 4.d0*spin)
        endif
                s = s * (2.d0*spin + 1.d0)/3.d0
c Here S is the spin degeneracy.

c Now include the degeneracy for the J levels.
        g2ch3cn = s*(2.d0*j+1)
c        write(22,*) 'j,k,s,g2ch3cn',j,k,s,g2ch3cn

        return

        end

        function gch3cn(m)

c Computes the statistical degeneracy of a J,K level.
c Identical to G2CH3CN expect that here there is one
c argument for the function which is the line number m 
c instead of two, the state numbers, j,k.
c The formula for the statistical deg of J levels is  2J+1.
c The Leiden data base has twice these values for both para and ortho
c so I included another factor of 2.

        include 'nlines_f77.h'
        implicit real *8 (a-h,o-z)
        integer parity
        common /index/ jindx(nstate),kindx(nstate),parity(nstate)
        spin=0.5d0
        j = jindx(m)
        k = kindx(m)

c First compute the degeneracy due to spin. This depends on K.
c The formulae are 3-46 of Townes & Schawlow
        if (mod(k,3) .eq. 0 ) then
            if (k .ne. 0) then
               s = 2.*(4.d0*spin**2 + 4.d0*spin + 3.d0)
            else
               s = (4.d0*spin**2 + 4.d0*spin + 3.d0)
            endif
        else
                s = 2.*(4.d0*spin**2 + 4.d0*spin)
        endif
        s = s * (2.d0*spin + 1.d0)/3.d0
c Here S is the spin degeneracy. 

c Now include the degeneracy for the J levels.
        gch3cn = s*(2.d0*j+1)
c        write(22,*) 'j,k,s,gch3cn',j,k,s,gch3cn

        return

        end

        integer function ch3cnhyp(vwmin)

c the linewidth vwmin is not used in this function
        include 'nlines_f77.h'
        parameter (maxhyp=50)

        implicit real *8 (a-h,o-z)
        real *4 hypvel,relint

        common /procid/ myid
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /hyperfine/ nhyp(nlines),
     &          hypvel(maxhyp,nlines),relint(maxhyp,nlines)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /lines/ indexu(nlines),indexl(nlines)


c Common index was filled in by nh3int. In nh3int, the parity has
c the same memory space as lindx
        common /index/ jindx(nstate),kindx(nstate),lindx(nstate)

c Common ri contains the relative intensities and freq shifts of
c the NH3 hyperfines filled in by nh3int.
        common /ri/ ric(nstate,nstate,9),fsc(nstate,nstate,9)

        dimension hypfrq(maxhyp,nlines)


        write(6,*) 'Working out CH3CN hyperfine structure'
c nh3hyp is the return value for the function, but it does nothing
        ch3cnhyp = 1



        do 79 line = 1,nlines
            nhyp(line) = 0
            do 80 ihyp = 1,9
                if (ric(indexu(line),indexl(line),ihyp) .gt. 0.)
     &             nhyp(line) = nhyp(line) + 1
 80         continue
                    
            do 81 ihyp = 1,nhyp(line)
                hypfrq(ihyp,line) = fsc(indexu(line),indexl(line),ihyp)
                relint(ihyp,line) = ric(indexu(line),indexl(line),ihyp)
 81         continue
 
 79     continue
 




c For each line
        do 1 line=1,nlines

c  Convert to velocity.

            if (myid .lt. 2) 
     &  write(22,*) line,indexu(line),indexl(line),
     &         freq(indexu(line),indexl(line))
        frq = freq(indexu(line),indexl(line))

        do 3 ihyp = 1,nhyp(line)
          hypfrq(ihyp,line) = hypfrq(ihyp,line) + frq
          hypvel(ihyp,line) = -(hypfrq(ihyp,line)-frq)/frq*c

 3      continue


c Sort the hyperfines into ascending order in frequency

        do 31 ihyp = 1,nhyp(line)

            minvel = 1
            velmin = 1.e20
            do 32 j = ihyp,nhyp(line)
                if (hypvel(j,line) .lt. velmin) then
                        minvel = j
                        velmin = hypvel(j,line)
                endif
 32         continue

            tmp1 = relint(ihyp,line)
            tmp2 = hypfrq(ihyp,line)
            tmp3 = hypvel(ihyp,line)
            relint(ihyp,line) = relint(minvel,line)
            hypfrq(ihyp,line) = hypfrq(minvel,line)
            hypvel(ihyp,line) = hypvel(minvel,line)
            relint(minvel,line) = tmp1
            hypfrq(minvel,line) = tmp2
            hypvel(minvel,line) = tmp3

 31     continue

            if (myid .lt. 2) then
        write(22,*) 'new line # ',line,
     &          freq(indexu(line),indexl(line))/1.e9,' GHz'
        write(22,*) 'number, hyperfine freq shift in MHz, ',
     &          ' shift in kms, relative intensity'
        do 4 ihyp = 1,nhyp(line)
          write(22,1000) ihyp,hypfrq(ihyp,line)/1.e6,
     &          hypvel(ihyp,line)/1.e5,
     &          relint(ihyp,line)
 4      continue
 1000   format(i3,2f14.4,f10.6)
            endif

 1      continue

        return

        end

        function ltech3cn(tk,arg_frac)

        include 'nlines_f77.h'

        implicit real *8 (a-h,o-z)
        real *4 tk

        common /frac/ frac(nstate),xfrac(nstate),gg(nstate)
        integer parity
        common /index/ jindx(nstate),kindx(nstate),parity(nstate)
        common /temp/ ilev

c arg_frac is used to pass the values of common frac through the args
        dimension arg_frac(nstate)

c ltech3cn is just to provide a return value
        ltech3cn = 1

c Compute the first 20 levels only.
        q = 0.d0
        do 3 jj = 1,20
        do 3 kk = 1,jj
                j = jj - 1
                k = kk - 1
                q = q + fch3cn(j,k,tk)

 3      continue

c Loop on m levels to set up pop in NSTATE=16 levels including radio
c frq inv

        do 5 m = 1,nstate
                frac(m) = fch3cn(jindx(m),kindx(m),tk)/q
                arg_frac(m) = frac(m)
c       write(6,*) 'q,frac ',q,frac(m)
 5      continue

        return
        end

        function fch3cn(j,k,tk)

c Computes the fractional population in each J,K level of a
c symmetric top. The formula is 3-44 of Townes and Schawlow.

        include 'nlines_f77.h'
        implicit real *8 (a-h,o-z)
        real *4 tk
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /abc / arot11,brot11,crot11,rc,inuc
         br = brot11*1.d11
         ar = arot11*1.d11

        tk8 = dble(tk)

c This is the exponential argument
        x = ((BR*J*(J+1) + (AR-BR)*K**2) )*
     &          planck/boltz/Tk8

c Use either the CH3CN or NH3 degeneracies
                fch3cn = g2ch3cn(j,k) * EXP(-x)
c       write(22,*) 'x,flte',x,fch3cn

        return
        END

