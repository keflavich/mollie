c These are the oldest subroutines and functions in the code.
c In most of these there is the implicit double precision statement
c        implicit real *8 (a-h,o-z)
c Some of the variables that connect to the subroutines and functions
c in FSUB.F are not double precision. 
c In particular, the variables in common block /MN/ are real *4
c Also common block /HYPERFINE/ is real *4
c Other differences. There is a function with the name C here 
c In FSUB.F C is a variable = speed of light
c The rotational constants here AROT11, BROT11 are in units of 10^11
c so these values are about one.
c In FSUB.F, BROT has a value of a few times 10^11.
c The common block ABC is used here to hold rotational constants
c In FSUB.F, the common block MN is used to hold molecular constants

        subroutine nh3int

        include 'nlines_f77.h'
        parameter (maxtrans=500,maxhyp=50)

        implicit real *8 (a-h,o-z)
        real *4 fnh3,atoms,brot,dipole
        integer parity,upper,lower
        logical docont,doprint,doprint2,nh3tau
        common /procid/ myid
        common /index/ jindx(nstate),kindx(nstate),parity(nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /ri/ ric(nstate,nstate,9),fsc(nstate,nstate,9)
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)
        common /frac/ frac(nstate),xfrac(nstate),gg(nstate)
        common /lines/ indexu(nlines),indexl(nlines)
        common /abc / arot11,brot11,crot11,rc,inuc
        common /nltehyp/ upper(maxtrans),lower(maxtrans),
     &          istartline(nlines)
        common /einstein/ aulhyp(maxhyp,nlines),
     &                      bulhyp(maxhyp,nlines),bluhyp(maxhyp,nlines)


c constants for NH3, both ortho and para. AROT11 and BROT11 are
c the same rotational constants as CROT and BROT, except for the
c units of 10^11. AROT11 and BROT11 are used in the TFF and SFF
c functions that define the hyperfine relative intensities.
c CROT and BROT are used in SETFREQ, SETFRQ33, EINSTA etc where
c frequencies are defined.

c These are for NH3

            atoms = 17.0
            arot11 = 2.9811706d0
            crot11 = 1.8672636d0
            brot11 = arot11
            rc = 4.09d6
            dipole = 1.469711d-18
            inuc = 1

            if (myid .lt. 2) then
                write(22,*) 'Loading NH3 constants  '
                write(6,*) 'Loading NH3 constants  '
            endif

                jindx(1)=1
                kindx(1)=1
                parity(1)=0

                jindx(2)=1
                kindx(2)=1
                parity(2)=1

        if (nstate .ge. 4) then
                jindx(3)=2
                kindx(3)=2
                parity(3)=1

                jindx(4)=2
                kindx(4)=2
                parity(4)=0
        endif

        if (nstate .ge. 6) then
                jindx(5)=3
                kindx(5)=3
                parity(5)=1

                jindx(6)=3
                kindx(6)=3
                parity(6)=0
        endif

        if (nstate .ge. 8) then
                jindx(7)=4
                kindx(7)=4
                parity(7)=1

                jindx(8)=4
                kindx(8)=4
                parity(8)=0
        endif

        if (nstate .ge. 10) then
                jindx(9)=5
                kindx(9)=5
                parity(9)=0

                jindx(10)=5
                kindx(10)=5
                parity(10)=1
        endif

        if (nstate .ge. 12) then
                jindx(11)=6
                kindx(11)=6
                parity(11)=1

                jindx(12)=6
                kindx(12)=6
                parity(12)=0
        endif

                do 101 n = 1,nlines
                    indexl(n) = 2*n-1 
                    indexu(n) = 2*n
 101            continue
 


                call fsri
                call normal
                call comprs
                call setfreq(molecule)

c Statistical weights
c statdg is in common afrc transferred to cmain
c gg is for older f77 code
        do 811 kk = 1,nstate
                gg(kk) = gnh3(kk)
                statdg(kk) = gg(kk)
 811    continue
                call einsta

        do 440 line = 1,nlines
           istartline(line) = line - 1
           n = indexu(line)
           m = indexl(line)
           aline(line) = a(n,m)
           aulhyp(1,line) = aline(line)
           freqline(line) = freq(n,m)
           glline(line) = statdg(m)
           guline(line) = statdg(n)
           lower(line)  = line - 1
           upper(line)  = line
 440    continue


                lun11 = 22
            if (myid .lt. 2) then 
        call printr(lun11)

        do 30 line = 1,nlines
        m = indexu(line)
        n = indexl(line)
 
        do 77 lunsel = 2,2
        lun = 6
        if (lunsel .eq. 2) lun = 22
        write(lun,*) ' '
        write(lun,*) '------------------------------'
        write(lun,*) ' '
        write(lun,*) ' NH3 '
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
        write(lun,*) '------------------------------'
        write(lun,*) ' '
 77     continue

 30     continue
            endif


        return
        end
 
        function c(j,f,i)
c See comments in subroutine FSRI. A quantum mechanical formula:
        implicit real *8 (a-h,o-z)
        integer f
        c = f*(f + 1) - i*(i + 1) - j*(j + 1)
        return
        end


        subroutine tff(j,f,m,n,li,lf)
 
c See comments in subroutine FSRI.
 
        include 'nlines_f77.h'
        implicit real *8 (a-h,o-z)
        integer f
        common /ricomp/ ri(nstate,nstate,3,3),fs(nstate,nstate,3,3)
        common /abc / arot11,brot11,crot11,rc,inuc
        a = arot11
        i = inuc
 
 
        ri(m,n,li,lf) = a*((j*(j+1) + f*(f+1) - i*(i + 1))**2)
     &          * (2*f+1)/(f*(f+1))
 
 
        call delfrq(m,n,li,lf)
 
        return
        end
 
        subroutine tffp1(j,f,m,n,li,lf)
 
c See comments in subroutine FSRI.
 
        include 'nlines_f77.h'
        implicit real *8 (a-h,o-z)
        integer f
        common /ricomp/ ri(nstate,nstate,3,3),fs(nstate,nstate,3,3)
        common /abc / arot11,brot11,crot11,rc,inuc
        a = arot11
        i = inuc
 
 
        ri(m,n,li,lf) = -a*(j + f + i + 2)*(j + f - i + 1)
     &          * (j - f + 1)*(j - f - i - 1)/(f + 1)
 
 
        call delfrq(m,n,li,lf)
 
        return
        end
 
        subroutine tffm1(j,f,m,n,li,lf)
 
c See comments in subroutine FSRI.
 
        include 'nlines_f77.h'
        implicit real *8 (a-h,o-z)
        integer f
        common /ricomp/ ri(nstate,nstate,3,3),fs(nstate,nstate,3,3)
        common /abc / arot11,brot11,crot11,rc,inuc
        a = arot11
        i = inuc
 
        ri(m,n,li,lf) = -a*(j + f + i + 1)*(j + f - i)
     &          * (j - f + i + 1)*(j - f - i)/f
 
        call delfrq(m,n,li,lf)
 
        return
        end


        subroutine einsta
 
c This subroutine computes the Einstein A coefficients of the allowed
c NH3 transitions.
 
        include 'nlines_f77.h'
        parameter (ndim=nstate-1)

        implicit real *8 (a-h,o-z)
        real *4 fnh3,atoms,brot,dipole
        integer parity

        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)

        common /index/ jindx(nstate),kindx(nstate),parity(nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole

        c = 3.d10
        planck = 6.625d-27
        flag = 0.d0
 
        pi = 4.d0*datan(1.d0)
 
        do 2 m = 1,nstate
        do 2 n = 1,nstate
                a(m,n) = flag
 2      continue
 
 
        do 1 m = 1,nstate
        do 1 n = 1,m-1
 
                ji = jindx(m)
                jf = jindx(n)
                ki = kindx(m)
                kf = kindx(n)
 
c Selection rules: transitions across k-ladders and to
c same parity state are forbidden.
                if (ki .ne. kf) go to 1
                if ( (molecule .eq. 3) .and.
     &                  (parity(m) + parity(n) .ne. 1 ) ) go to 1

c Factor expresses the j,k state dependence of the square of the
c dipole moment. There are 2 formulas depending on whether j <-- j
c or j-1 <-- j. See eqn. 3-41,42 of Townes and Schawlow. See also
c 3-40 for another formula jf-ji=1 which we dont use for einstein a
c since we want downward rates.
                if (ji .eq. jf) then
                        factor = ki**2
                        factor = factor/(ji*(ji + 1))
                else
                        if ((ji-jf) .ne. 1) go to 1
                        factor = (ji**2 - ki**2)/(ji*(2.*ji + 1.))
                endif
 
                coeff = 64.d0*pi**4*freq(m,n)*dipole*freq(m,n)
     &                  *dipole*freq(m,n)/(3.d0*c**2*planck*c)
 
                a(n,m) = coeff * factor
                a(m,n) = a(n,m)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c               a(m,n) = 0.0
c               a(n,m) = 0.0
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 1      continue
 
        return
        end
 

        subroutine compar(a,b,eps,itest)
c Compares two numbers and returns ITEST
        implicit real *8 (a-h,o-z)
        itest = 0
c        if((a*b) .lt. 0.d0)return
        if(abs(a-b) .gt. eps)return
        itest = 1
        return
        end
 
        subroutine delfrq(m,n,li,lf)
 
c See comments in subroutine FSRI. This routine computes the frequency
c spacing of the hyperfine lines using equation 6-8 from Townes and
c Schawlow. i is the nuclear spin, rc the quadrupole coupling constant
c in frequency units. j,k,f states are initial or final as indicated
c by i or f. li and lf are passed by the calling routine and are related
c to the quadrupole hyperfine state numbers by the algorithm below.
 
        include 'nlines_f77.h'
        implicit real *8 (a-h,o-z)
        integer fi,ff,parity
        common /ricomp/ ri(nstate,nstate,3,3),fs(nstate,nstate,3,3)
        common /index/ jindx(nstate),kindx(nstate),parity(nstate)
        common /abc / arot11,brot11,crot11,rc,inuc
        i = inuc
        a = arot11
        b = brot11

        ji = jindx(m)
        jf = jindx(n)
        ki = kindx(m)
        kf = kindx(n)
 
        fi = (li - 2) + jindx(m)
        ff = (lf - 2) + jindx(n)
 
        if (kf.eq.0) then
        wf = rc*(- 1.)
     &          /(2.*i*(2.*i - 1.)*(2.*jf - 1.)*(2.*jf + 3.))
     &          *(3./4.*c(jf,ff,i)*(c(jf,ff,i) + 1.)
     &          - i*(i + 1.)*jf*(jf + 1.))
        else
        wf = rc*(3.*(kf**2)/(jf*(jf + 1.)) - 1.)
     &          /(2.*i*(2.*i - 1.)*(2.*jf - 1.)*(2.*jf + 3.))
     &          *(3./4.*c(jf,ff,i)*(c(jf,ff,i) + 1.)
     &          - i*(i + 1.)*jf*(jf + 1.))
        endif
 
        if (ki.eq.0) then
        wi = rc*(- 1.)
     &          /(2.*i*(2.*i - 1.)*(2.*ji - 1.)*(2.*ji + 3.))
     &          *(3./4.*c(ji,fi,i)*(c(ji,fi,i) + 1.)
     &          - i*(i + 1.)*ji*(ji + 1.))
        else
        wi = rc*(3.*(ki**2)/(ji*(ji + 1.)) - 1.)
     &          /(2.*i*(2.*i - 1.)*(2.*ji - 1.)*(2.*ji + 3.))
     &          *(3./4.*c(ji,fi,i)*(c(ji,fi,i) + 1.)
     &          - i*(i + 1.)*ji*(ji + 1.))
        endif
 
        fs(m,n,li,lf) = wf - wi
 
 
        return
        end

 
        subroutine sff(j,f,m,n,li,lf)
 
c See comments in subroutine FSRI.
 
        include 'nlines_f77.h'

        implicit real *8 (a-h,o-z)
        integer f
        common /ricomp/ ri(nstate,nstate,3,3),fs(nstate,nstate,3,3)
        common /abc / arot11,brot11,crot11,rc,inuc
        b = brot11
        i = inuc
 
        ri(m,n,li,lf) = -b*(j + f + i + 1)*(j + f - i)
     &          * (j - f + 1)*(j - f - i - 1)*(2*f + 1)/(f*(f + 1))
 
        call delfrq(m,n,li,lf)
 
        return
        end
 
        subroutine sffm1(j,f,m,n,li,lf)
 
c See comments in subroutine FSRI.
 
        include 'nlines_f77.h'

        implicit real *8 (a-h,o-z)
        integer f
        common /ricomp/ ri(nstate,nstate,3,3),fs(nstate,nstate,3,3)
        common /abc / arot11,brot11,crot11,rc,inuc
        b = brot11
        i = inuc
 
        ri(m,n,li,lf) = b*(j + f + i + 1)*(j + f + i)
     &          * (j + f - 1)*(j + f - i - 1)/f
 
        call delfrq(m,n,li,lf)
 
        return
        end
 
        subroutine sffp1(j,f,m,n,li,lf)
 
c See comments in subroutine FSRI.
 
        include 'nlines_f77.h'

        implicit real *8 (a-h,o-z)
        integer f
        common /ricomp/ ri(nstate,nstate,3,3),fs(nstate,nstate,3,3)
        common /abc / arot11,brot11,crot11,rc,inuc
        b = brot11
        i = inuc
 
        ri(m,n,li,lf) = b*(j - f + i)*(j - f + i - 1)
     &          * (j - f - i - 1)*(j - f - i - 2)/(f + 1)
 
        call delfrq(m,n,li,lf)
 
        return
        end

        subroutine setfreq(molecule)
 
c This routine computes the energy difference between the rotational and
c inversion levels in units of frequency. The energy  of the rotational
c states, W, is computed from a formula (eqn. 3-5 of Townes and Schawlow).
c The energy of the inversion
c transitions is put in by hand in the DATA statement. The program further
c assumes that the energies of the two inversion levels are symmetrically
c displaced around the energy computed for the rotational state. The results
c are stored in COMMON /AFRC/.
 
c For this routine to work properly (not compute negative frequencies)
c the levels have to be in ascending order in energy. This is way the
c J,K states of NH3 were loaded into nh3int
 
        include 'nlines_f77.h'
        parameter (ndim=nstate-1)
        implicit real *8 (a-h,o-z)
        integer parity
        common /afrc/ a(nstate,nstate),freq(nstate,nstate),
     &          statdg(nstate),aline(nlines),freqline(nlines),
     &          guline(nlines),glline(nlines)

        common /index/ jindx(nstate),kindx(nstate),parity(nstate)
        common /frqinv/ frqinv(nstate),sgn(nstate)
        common /abc / arot11,brot11,crot11,rc,inuc
        common /lines/ indexu(nlines),indexl(nlines)


        dimension k(2),j(2),w(2),sgntmp(2)

        planck = 6.625d-27
        boltz = 1.3806d-16
        flag = 0.d0

c Put in the inversion transition frequencies by hand.
c Put in different transitions for para and ortho NH3
c molecule=3 is for para    

 
            frqinv(1) = 23.694495 d9
            if (nstate .ge. 4)  frqinv(3) = 23.722633 d9
            if (nstate .ge. 6)  frqinv(5) = 23.870296 d9
            if (nstate .ge. 8)  frqinv(7)= 24.13939  d9
            if (nstate .ge. 10) frqinv(9)= 24.53299  d9
            if (nstate .ge. 12) frqinv(11) = 25.056025 d9
            cc = crot11 * 1.d11 
            bc = brot11 * 1.d11
 
c Here we set which level of the inversion
c level is upper and which is lower 

            do 7 i = 2,nstate,2
                frqinv(i) = frqinv(i-1)
                sgn(i) = 0.5d0
                sgn(i-1) = -0.5d0
 7          continue

c That finishes the NH3 inversion transitions


c Initialize all frequencies to flag=zero
        do 1 m = 1,nstate
        do 1 n = 1,nstate
                freq(m,n) = flag
 1      continue
 
c        do 99 i = 1,16
c 99        write(6,*),i,frqinv(i) 
 
        do 22 m = 1,nstate
        do 2 n = 1,m-1
 
                j(1) = jindx(m)
                j(2) = jindx(n)
                k(1) = kindx(m)
                k(2) = kindx(n)
                sgntmp(1) = sgn(m)
                sgntmp(2) = sgn(n)
 
                do 3 l = 1,2
c w is the energy of the j,k state not accounting for the inversion
c splitting.
        w(l) = (bc*j(l)*(j(l) + 1) + (cc - bc)*k(l)**2)
c       write(22,*),m,n,j(l),k(l),w(l)
 3              continue
c Add in the inversion splitting assuming it is symmetric about
c the energy level.
                freq(m,n) = w(1) + sgntmp(1)*frqinv(m)
     &                  - (w(2) + sgntmp(2)*frqinv(n))
c            if (myid .lt. 2)
c     &        write(22,*),'energy level ',m,freq(m,n)
 2      continue

 22     continue

c Since the frequency may be used to compute the energy separation
c fill in the reverse half of the array.
        do 6 m = 1,nstate
        do 6 n = m+1,nstate
                freq(m,n) = freq(n,m)
 6      continue
 
 
        return
        end
 



        subroutine comprs
 
c There are 9 transitions between the 3 hyperfine levels of each rotational
c quantum state. Because of symmetry, some of these transitions between any
c two particular levels will have the same frequency and thus form one
c line with a relative intensity which is the sum of the relative intensities
c of the contributing transitions. This routine compares the frequencies of
c the 9 transitions to some specified precision, EPS, which should have some
c relation to the observed linewidth, sums the relative intensities of
c lines at the same frequencies, and eliminates one of the lines from the
c array of 9 possible lines. The resulting fewer lines are stored in the
c arrays  RIC and FSC.

c l is the number of remaining lines
 
        include 'nlines_f77.h'
        implicit real *8 (a-h,o-z)
        logical doprint2
        common /procid/ myid
        common /ricomp/ ri(nstate,nstate,3,3),fs(nstate,nstate,3,3)
        common /ri/ ric(nstate,nstate,9),fsc(nstate,nstate,9)
        common /flag/ flag
        dimension tmpri(9),tmpfs(9)
c        data eps /1.d3/
        eps = 1.d3
        doprint2 = .false.
            if (myid .lt. 2) doprint2 = .false.
 
 
        do 1 m = 1,nstate
        do 1 n = 1,m
c        if (m .eq. 5 .and. n .eq. 4) doprint2 = .true.
 
                do 6 kk = 1,9
                        ric(m,n,kk)=flag
                        fsc(m,n,kk)=flag
 6              continue
 
            l = 0
            do 2 j = 1,3
            do 2 k = 1,3
                l = l + 1
                tmpri(l) = ri(m,n,j,k)
                tmpfs(l) = fs(m,n,j,k)
 2          continue
 
 
                do 7 jj = 1,9
                    if (tmpri(jj) .ne. flag) then
                        ric(m,n,1) = tmpri(jj)
                        fsc(m,n,1) = tmpfs(jj)
       if (doprint2)
     &    write(22,*)'1st loaded,jj,ric,fsc',jj,ric(m,n,1),fsc(m,n,1)
                        go to 8
                    endif
 7              continue
 8           continue
 
            if (jj .gt. 8) go to 1
            l = 1
            do 3 j = jj+1,9
        if (doprint2) write(22,*)'start comparing. Is line valid?'
        if (doprint2) write(22,*)'j,tmpri,tmpfs',j,tmpri(j),tmpfs(j)
                if(tmpri(j) .eq. flag) go to 3
        if (doprint2) write(22,*)'Line is valid: ',
     &     ' Compare to prev lines. l = ',l
                do 4 k = 1,l
                call compar(tmpfs(j),fsc(m,n,k),eps,itest)
        if (doprint2) write(22,*)'is it ',k,'th fsc? :itest ',
     &     fsc(m,n,k),itest
                if (itest .eq. 1) then
        if (doprint2) write(22,*)'here it is'
                        ric(m,n,k) = ric(m,n,k) + tmpri(j)
                        go to 3
                else

        if (doprint2) write(22,*)'not this one'
                endif
 4      continue
        if (doprint2) write(22,*)'couldnt match. must be different'
                        l = l + 1
                        ric(m,n,l) = tmpri(j)
                        fsc(m,n,l) = tmpfs(j)
 3          continue
 1      continue

 
        do 19 m = 1,nstate
        do 19 n = 1,m
        do 19 k = 1,9
                ric(n,m,k) = ric(m,n,k)
                fsc(n,m,k) = fsc(m,n,k)
 19     continue

        return
        end
 

        subroutine normal
 
c This routine is called after subroutine FSRI and normalizes the
c relative intensities of the hyperfine lines.
 
        include 'nlines_f77.h'
        implicit real *8 (a-h,o-z)
        common /ricomp/ ri(nstate,nstate,3,3),fs(nstate,nstate,3,3)
        common /flag/ flag
 
        do 1 m = 1,nstate
        do 1 n = 1,m
 
                sum = 0.
 
                do 2 li = 1,3
                do 2 lf = 1,3
 
                if (ri(m,n,li,lf) .ne. flag) sum = ri(m,n,li,lf) + sum
 
 2              continue
 
                if (sum .le. 0) go to 1
 
                do 3 li = 1,3
                do 3 lf = 1,3
 
                if (ri(m,n,li,lf) .gt. flag)
     &                  ri(m,n,li,lf) = ri(m,n,li,lf)/sum
 
 3              continue
 1      continue
 
        return
        end
 
        subroutine fsri
 
c Calculates frequency spacing and relative intensities of the
c hyperfine lines according to formulas of quantum mechanics. See
c equation 6-6a,b in Townes and Schawlow. The
c routine uses the J,K indices stored in the common block /INDEX/
c and fills the arrays in the common block /RICOMP/. The routine
c loops through nstate number of J,K states and for each J,K state
c through the 9 possible hyperfine transitions. The routine determines
c which of the quantum mechanical formulas applies to a particular
c transition and then calls one of six subroutines, TFF, TFFP1,
c TFFM1, SFF, SFFP1, SFFM1, each of which contains one formula for
c determining the relative intensity and a call to DELFRQ which
c computes the frequency spacing.
c
c The subroutines are named so that
c       ff   = f <-- f
c       ffp1 = f <-- f+1
c       ffm1 = f <-- f-1
c Subroutines starting with t are for j <-- j. Starting with s for
c j <-- j-1 and j-1 <-- j. For the latter, the arrows on the f
c transitions are also reversed.
c
c li and lf are do-loop variables used to set the fi and ff states.
c ji,fi are initial states, jf,ff are final states.
 
        include 'nlines_f77.h'
        implicit real *8 (a-h,o-z)
        integer fi,ff,parity
        real *4 fnh3,atoms,brot,dipole
        common /ricomp/ ri(nstate,nstate,3,3),fs(nstate,nstate,3,3)
        common /index/ jindx(nstate),kindx(nstate),parity(nstate)
        common /mn/ molecule,fnh3,atoms,brot,dipole
        common /flag/ flag
        flag = 0.d0
 
        do 1 m = 1,nstate
        do 1 n = 1,nstate
        do 1 li = 1,3
        do 1 lf = 1,3
                ri(m,n,li,lf) = flag
                fs(m,n,li,lf) = flag
 1      continue
 
 
        do 2 m = 1,nstate
        do 2 n = 1,m
 
                ji = jindx(m)
                jf = jindx(n)
 
c Selection rules: transitions across k-ladders and to 
c same parity are not allowed.

                if (molecule .eq. 2 .and. 
     &                        (jindx(m) - jindx(n)) .ne. 1) go to 2
                if (kindx(m) .ne. kindx(n)) go to 2
                if ((molecule .eq. 3 .or. molecule .eq. 4) 
     &               .and. (parity(m) + parity(n) .ne. 1))  go to 2
 
c now select which formula to use.
 
                if (ji .eq. jf) then
                        do 3 li = 1,3
                        do 3 lf = 1,3
                        fi = (li - 2) + jindx(m)
                        ff = (lf - 2) + jindx(n)
                if (fi .eq. ff .and. ff .ne. 0)
     &                  call tff  (jf,ff,m,n,li,lf)
                if (fi .eq. ff + 1)
     &                  call tffp1(jf,ff,m,n,li,lf)
                if (fi .eq. ff - 1 .and. ff .ne. 0)
     &                  call tffm1(jf,ff,m,n,li,lf)
 
 3                      continue
                endif
 
                if (jf .eq. ji - 1) then
                        do 4 li = 1,3
                        do 4 lf = 1,3
                        fi = (li - 2) + jindx(m)
                        ff = (lf - 2) + jindx(n)
                if (fi .eq. ff .and. fi .ne. 0)
     &                  call sff  (ji,fi,m,n,li,lf)
                if (ff .eq. fi - 1 .and. fi .ne. 0)
     &                  call sffm1(ji,fi,m,n,li,lf)
                if (ff .eq. fi + 1) 
     &                  call sffp1(ji,fi,m,n,li,lf)
 4                      continue
                endif
 
c There are never any of these transitions.
 
                if (ji .eq. jf - 1) then
                        do 5 li = 1,3
                        do 5 lf = 1,3
                        fi = (li - 2) + jindx(m)
                        ff = (lf - 2) + jindx(n)
                if (fi .le. 0 .or. ff .le. 0) go to 5
                if (fi .eq. ff .and. ff .ne. 0)
     &                  call sff(jf,ff,m,n,li,lf)
                if (fi .eq. ff + 1) call sffp1(jf,ff,m,n,li,lf)
                if (fi .eq. ff - 1 .and. ff .ne. 0)
     &                  call sffm1(jf,ff,m,n,li,lf)
 5                      continue
                endif
 
 2      continue
 
 
        return
        end


        function ltenh3(tk,arg_frac)

        include 'nlines_f77.h'

        implicit real *8 (a-h,o-z)
        integer parity
        real *4 tk

        common /frac/ frac(nstate),xfrac(nstate),gg(nstate)
        common /index/ jindx(nstate),kindx(nstate),parity(nstate)
        common /temp/ ilev

c arg_frac is used to pass the values of common frac through the args
        dimension arg_frac(nstate)

c ltenh3 is just to provide a return value
        ltenh3 = 1

c Set m = 0 to calculate partition function. Ignores pop diff
c in radio freq inv transitions. Compute the first 20 levels only.
        m = 0
        q = 0.d0
        do 3 jj = 1,20
        do 3 kk = 1,jj
                j = jj - 1
                k = kk - 1
                levs = 2
                if (k .eq. 0) levs = 1

                q = q + levs*fltenh3(j,k,m,tk)

 3      continue
 
c Loop on m levels to set up pop in NSTATE=16 levels including radio 
c frq inv
 
        do 5 m = 1,nstate
                frac(m) = fltenh3(jindx(m),kindx(m),m,tk)/q
                arg_frac(m) = frac(m)
c        write(6,*) 'q,frac ',q,frac(m)
 5      continue   
  
        return
        end

        function fltenh3(j,k,m,tk)

c Computes the fractional population in each J,K level of a
c symmetric top. The formula is 3-44 of Townes and Schawlow.

        include 'nlines_f77.h'
        implicit real *8 (a-h,o-z)
        real *4 tk
        real *4 fnh3,atoms,brot,dipole
        common /frqinv/ frqinv(nstate),sgn(nstate)
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /abc / arot11,brot11,crot11,rc,inuc
         br = brot11*1.d11
         cr = crot11*1.d11

c If m=0 it means we should ignore the radio freq inversion
c splitting. Use this to calculate the partition function.
c If m is not zero it is the pointer to an array element and
c we should use the splitting to get the population in that
c level. 

        if (m.eq.0) then
                frqtmp = 0.d0
                sgntmp = 0.d0
        else
                frqtmp = frqinv(m)
                sgntmp = sgn(m)
        endif

c        if (tk .lt. 4. .or. tk .gt. 4000.) write(6,*) 'bad tk ',tk
        tk8 = dble(tk)
 
c This is the exponential argument
        x = ((BR*J*(J+1) + (CR-BR)*K**2) + sgntmp*frqtmp)*
     &          planck/boltz/Tk8

                fltenh3 = g2nh3(j,k)  * EXP(-x)
c       write(22,*) 'big,x,flte',big,x,fltenh3
 
        return
        END


        function g2nh3(j,k)
c Computes the statistical degeneracy of a J,K level.
c This function is identical to GNH3 except this takes 2 ARGS J,K
c This is used in the partition function so that we can sum
c more levels than NSTATE

c The formula for the statistical deg of J levels is 2J+1.
c The Leiden data base has twice these values for both para and ortho
c so I included another factor of 2
        implicit real *8 (a-h,o-z)
        if (mod(k,3) .eq. 0) then
                s = 2.d0
        else
                s = 1.d0
        endif
        g2nh3 = s*2.d0*(2.d0*j+1)
        return
        end

        function gnh3(m)
c Computes the statistical degeneracy of a J,K level.
c This function is identical to G2NH3 except this takes 1 ARG M
c The formula for the statistical deg of J levels is  2J+1.

        include 'nlines_f77.h'
        implicit real *8 (a-h,o-z)
        integer parity
        common /index/ jindx(nstate),kindx(nstate),parity(nstate)
        j = jindx(m)
        k = kindx(m)
        if (mod(k,3) .eq. 0) then
                s = 2.d0
        else
                s = 1.d0
        endif
        gnh3 = s*(2.d0*j+1)
c        write(22,*) 'j,k,s,gnh3',j,k,s,gnh3
        return
        end


 
        subroutine printr(lun)
 
c This routine prints out the relative intensities and frequency spacings
c of the hyperfine lines in a tidy format. Writes to logical unit LUN.
 
        implicit real *8 (a-h,o-z)
        include 'nlines_f77.h'
        integer fi,ff,parity
        common /ricomp/ ri(nstate,nstate,3,3),fs(nstate,nstate,3,3)
        common /index/ jindx(nstate),kindx(nstate),parity(nstate)
        character *2 plus,minus,symm(2)
        plus = '+'
        minus = '-'
 
        npage = 0
 
        do 1 m = 1,nstate
        do 1 n = 1,m
 
                if (kindx(m) .ne. kindx(n)) go to 1
                if (iabs(jindx(m)-jindx(n)) .gt. 1) go to 1
                if (parity(m) + parity(n) .ne. 1) go to 1
 
                symm(2) = plus
                if (parity(m) .eq. 0) symm(2) = minus
                symm(1) = plus
                if (parity(n) .eq. 0) symm(1) = minus
                if (mod(npage,4) .eq. 0 .and. npage .ne. 0)
     &          write (lun,115)
 115            format(1h1)
                npage = npage + 1
                write(lun,104)
                write(lun,100)jindx(n),kindx(n),symm(1),
     &                  jindx(m),kindx(m),symm(2)
                write(lun,301)
                write(lun,102)
 
 
 104    format(////)
 100    format(1x,'(',i2,',',i2,')',a1,' <-- (',i2,',',i2,')',a1)
 301    format(1x,60('-'))
 102    format(1x,t15,'freq',t30,'rel. int.')
 
                do 2 li = 1,3
                do 2 lf = 1,3
 
                        if (ri(m,n,li,lf) . gt. 0.) then
        fi = (li - 2) + jindx(m)
        ff = (lf - 2) + jindx(n)
        write(lun,103) ff,fi,fs(m,n,li,lf)/1.e6,ri(m,n,li,lf)
 103    format(1x,i2,'<--',i2,t15,f8.4,t30,f8.4)
 
                        endif
 2              continue
 1      continue
 
        return
        end

c        double precision frac,planck,boltz,c,pi,amu,pc,sqrtpi,
c     &  grav,a,freq,chfreq,statdg,cr,aline,freqline,guline,glline


        integer function nh3hyp(vwmin)

c the linewidth vwmin is not used in this function
        include 'nlines_f77.h'
        parameter (maxhyp=50)

        implicit real *8 (a-h,o-z)
        real *4 hypvel,relint

        common /procid/ myid
        common /consts/ planck,boltz,c,pi,amu,pc,sqrtpi,grav
        common /hyperfine/ nltehyp,nhyp(nlines),
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


c nh3hyp is the return value for the function, but it does nothing
        nh3hyp = 1



        do 79 line = 1,nlines
            nhyp(line) = 0
            do 80 ihyp = 1,9
                if (ric(indexu(line),indexl(line),ihyp) .gt. 0.)
     &             nhyp(line) = nhyp(line) + 1
 80            continue
                        
            do 81 ihyp = 1,nhyp(line)
                hypfrq(ihyp,line) = fsc(indexu(line),indexl(line),ihyp) 
                relint(ihyp,line) = ric(indexu(line),indexl(line),ihyp) 
 81            continue

 79     continue


c For each line
        do 1 line=1,nlines

c  Convert to velocity.

         if (myid .lt. 2) 
     &       write(22,*) line,indexu(line),indexl(line),
     &         freq(indexu(line),indexl(line))

        frq = freq(indexu(line),indexl(line))


c This sets up the electric quadrupole hyperfine splitting
        do 3 ihyp = 1,nhyp(line)
          hypfrq(ihyp,line) = hypfrq(ihyp,line) + frq
          hypvel(ihyp,line) = -(hypfrq(ihyp,line)-frq)/frq*c
 3      continue

c For the (1,1), (2,2), and (3,3) lines we have data on the
c magnetic hyperfine splitting. The next 3 blocks replace
c the electric quad hyperfine splitting above with the
c magnetic splitting below.

c These data are from Shadi Chitsazeddeh, processed into velocities
c from the original data in Kukolich 1967 Phys Rev 156, 83

        if (line .eq. 1) then
        nhyp(1) = 18

        hypvel(18, 1) = -19.4155
        hypvel(17, 1) = -19.2763
        hypvel(16, 1) =  -7.6823
        hypvel(15, 1) =  -7.2406
        hypvel(14, 1) =  -7.1014
        hypvel(13, 1) =  -0.1191
        hypvel(12, 1) =  -0.0801
        hypvel(11, 1) =   0.0  
        hypvel(10, 1) =   0.0591
        hypvel( 9, 1) =   0.3226
        hypvel( 8, 1) =   0.4417
        hypvel( 7, 1) =   0.4562
        hypvel( 6, 1) =   0.5954
        hypvel( 5, 1) =   7.4831
        hypvel( 4, 1) =   7.6023
        hypvel( 3, 1) =   8.0195
        hypvel( 2, 1) =  19.9844
        hypvel( 1, 1) =  19.4527

        relint(18, 1) = 0.44 
        relint(17, 1) = 0.22 
        relint(16, 1) = 0.50 
        relint(15, 1) = 0.06 
        relint(14, 1) = 0.28 
        relint(13, 1) = 0.10 
        relint(12, 1) = 0.28 
        relint(11, 1) = 1.40 
        relint(10, 1) = 0.06 
        relint( 9, 1) = 0.90
        relint( 8, 1) = 0.10
        relint( 7, 1) = 0.06
        relint( 6, 1) = 0.11
        relint( 5, 1) = 0.06
        relint( 4, 1) = 0.50
        relint( 3, 1) = 0.28
        relint( 2, 1) = 0.22
        relint( 1, 1) = 0.44

        rsum = 0.
        do 100 i = 1,18
           hypvel(i,1) = hypvel(i,1)*1.e5
           hypfrq(i,1) = frq - (hypvel(i,1)/c) * frq
           rsum = rsum + relint(i,1)
 100    continue
        do 102 i = 1,18
            relint(i,1) = relint(i,1)/rsum
 102    continue

        endif

        if (line .eq. 2) then
        nhyp(2) = 21

        hypvel(21, 2) = -26.526
        hypvel(20, 2) = -26.011
        hypvel(19, 2) = -25.951
        hypvel(18, 2) = -16.382
        hypvel(17, 2) = -16.370
        hypvel(16, 2) = -15.854
        hypvel(15, 2) =  -0.589
        hypvel(14, 2) =  -0.531
        hypvel(13, 2) =  -0.502
        hypvel(12, 2) =  -0.013
        hypvel(11, 2) =  -0.004
        hypvel(10, 2) =   0.013
        hypvel( 9, 2) =   0.524
        hypvel( 8, 2) =   0.528
        hypvel( 7, 2) =   0.562
        hypvel( 6, 2) =  15.865
        hypvel( 5, 2) =  16.379
        hypvel( 4, 2) =  16.392
        hypvel( 3, 2) =  25.95
        hypvel( 2, 2) =  26.011
        hypvel( 1, 2) =  26.526

        relint(21, 2) = 0.020
        relint(20, 2) = 0.180
        relint(19, 2) = 0.100
        relint(18, 2) = 0.178
        relint(17, 2) = 0.124
        relint(16, 2) = 0.010
        relint(15, 2) = 0.100
        relint(14, 2) = 0.050
        relint(13, 2) = 0.056
        relint(12, 2) = 0.700
        relint(11, 2) = 2.387
        relint(10, 2) = 1.278
        relint( 9, 2) = 0.050
        relint( 8, 2) = 0.056
        relint( 7, 2) = 0.100
        relint( 6, 2) = 0.010
        relint( 5, 2) = 0.124
        relint( 4, 2) = 0.178
        relint( 3, 2) = 0.100
        relint( 2, 2) = 0.180
        relint( 1, 2) = 0.020

        rsum = 0.
        do 101 i = 1,21
           hypvel(i,2) = hypvel(i,2)*1.e5
           hypfrq(i,2) = frq - (hypvel(i,2)/c) * frq
           rsum = rsum + relint(i,2)
 101    continue
        do 103 i = 1,21
            relint(i,2) = relint(i,2)/rsum
 103    continue

        endif

        if (line .eq. 3) then

        nhyp(3) = 19

c This data is from Kukolich 1967, Phys Rev 156, 83 . 
c Missing some frequencies, but the relative intensities of the 
c missing lines are c small (1%). Due to the low optical depth
c and close spacing of the (3,3) magnetic hyperfines, 
c they are very difficult to resolve. Could use just the quadrupole
c splitting.Probably not worth the effort to have typed these in.

c 2 3.5 -> 3 4.5
        hypfrq( 1, 3) = -2324.089
        relint( 1, 3) = 0.11338
c 2 2.5 -> 3 3.5
        hypfrq( 2, 3) = -2312.492
        relint( 2, 3) = 0.07775
c 2 0.5 -> 3 1.5
        hypfrq( 3, 3) = -2304.667
        relint( 3, 3) = 0.03175 
c 2 1.5 -> 3 2.5
        hypfrq( 4, 3) = -2302.375
        relint( 4, 3) = 0.05079
c 4 3.5 -> 3 2.5
        hypfrq( 5, 3) = -1690.939
        relint( 5, 3) = 0.06150
c 4 4.5 -> 3 3.5
        hypfrq( 6, 3) = -1688.839
        relint( 6, 3) = 0.08185
c 4 2.5 -> 3 1.5
        hypfrq( 7, 3) = -1682.922
        relint( 7, 3) = 0.04592
c 4 5.5 -> 3 4.5
        hypfrq( 8, 3) = -1679.057
        relint( 8, 3) = 0.10714
c 3 -> 3
        hypfrq( 9, 3) =  (-64.308 - 1.302 + 59.143)/3.
        relint( 9, 3) = 2.71175
c 4 -> 4
        hypfrq(10, 3) =  (-50.183 + 0.433 + 51.127)/3.
        relint(10, 3) = 4.24474
c 2 -> 2
        hypfrq(11, 2) =  (-80.104 + 1.041 + 82.240)/3.
        relint(11, 2) = 1.63990
c 3 1.5 -> 4 2.5
        hypfrq(12, 3) =  1682.148
        relint(12, 3) = 0.10714
c 3 3.5 -> 4 4.5
        hypfrq(13, 3) =  1687.971
        relint(13, 3) = 0.08185
c 3 2.5 -> 4 3.5
        hypfrq(14, 3) =  1690.070
        relint(14, 3) = 0.06150
c 3 4.5 -> 4 5.5
        hypfrq(15, 3) =  1678.235
        relint(15, 3) = 0.10714
c 3 2.5 -> 2 1.5
        hypfrq(16, 3) =  2302.080
        relint(16, 3) = 0.05079
c 3 1.5 -> 2 0.5
        hypfrq(17, 3) =  2304.227
        relint(17, 3) = 0.03175
c 3 3.5 -> 2 2.5
        hypfrq(18, 3) =  2312.291
        relint(18, 3) = 0.07775
c 3 4.5 -> 2 3.5
        hypfrq(19, 3) =  2323.792
        relint(19, 3) = 0.11338
c missing frequencies	relative intensity
c 2 2.5 -> 3 1.5 	0.01270
c 2 2.5 -> 3 2.5	0.01659
c 2 3.5 -> 3 3.5	0.01296
c 4 2.5 -> 3 2.5	0.00738
c 4 3.5 -> 3 3.5	0.00972
c 4 4.5 -> 3 4.5	0.00744

c Figure 8.5 page 222 T&S reference has no data
c tables page 510

c frq was set above just after the DO 1 loop statement
        rsum = 0.
        do 301 ihyp = 1,nhyp(line)
          hypfrq(ihyp,line) = hypfrq(ihyp,line)*1.e3 + frq
          hypvel(ihyp,line) = -(hypfrq(ihyp,line)-frq)/frq*c
          rsum = rsum + relint(ihyp,3)
 301    continue

        do 203 i = 1,nhyp(line)
            relint(i,3) = relint(i,3)/rsum
 203    continue

        endif




c Sort the hyperfines into descending order in velocity

        do 31 ihyp = 1,nhyp(line)

            maxvel = 1
            velmax = -1.e20
            do 32 j = ihyp,nhyp(line)
                if (hypvel(j,line) .gt. velmax) then
                        maxvel = j
                        velmax = hypvel(j,line)
                endif
 32         continue

            tmp1 = relint(ihyp,line)
            tmp2 = hypfrq(ihyp,line)
            tmp3 = hypvel(ihyp,line)
            relint(ihyp,line) = relint(maxvel,line)
            hypfrq(ihyp,line) = hypfrq(maxvel,line)
            hypvel(ihyp,line) = hypvel(maxvel,line)
            relint(maxvel,line) = tmp1
            hypfrq(maxvel,line) = tmp2
            hypvel(maxvel,line) = tmp3

 31     continue



            if (myid .lt. 2) then
        write(22,*) 'new line # ',line,
     &          freq(indexu(line),indexl(line))/1.e9,' GHz'
        write(22,*) 'number, hyperfine freq in MHz, ',
     &          ' shift in kms, relative intensity'
        do 4 ihyp = 1,nhyp(line)
          write(22,1000) ihyp,hypfrq(ihyp,line)/1.e6,
     &          hypvel(ihyp,line)/1.e5,
     &          relint(ihyp,line)
 4      continue
 1000   format(i3,2f14.4,f10.6)
            endif

c        write(22,*) 'nh3hyp: line nhyp ',line,nhyp(line)

 1      continue

        return

        end
