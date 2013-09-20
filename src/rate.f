	SUBROUTINE rate(gamma,n,ndim,poptot,pop)

c INPUTS
c gamma     The rate matrix
c n         The number of levels of the molecule
c ndim      The allocated size of gamma in memory
c poptot    The total population. For example, 1.0

c OUTPUT
c pop       The level populations

c	solve the statistical rate equations
c
c	   sum(j=1,n) gamma(i,j)*pop(j)=sum(j=1,n) gamma(j,i)*pop(i)
c
c	for the populations pop(i), with auxilliary condition
c
c	   sum(i=1,n) pop(i) = poptot
c
c	The physical dimensions of gamma are ndim x ndim.
c	Note that the order of indices in gamma are opposite to
c	the usual ordering of the rate coefficients.  That is,
c
c		gamma(i,j)=rate(j,i)
c
c       Any inputted diagonal elements gamma(j,j) are ignored, but
c       the appropriate values are generated within the routine,
c       that is,
c
c               gamma(j,j) = - sum(k=1,n, k.ne.j) gamma(k,j)
c
c
c       The solution is done by a LU decomposition of the matrix gamma.
c       The decomposition is applied to the full homogeneous system and only
c       afterwards is the condition of fixed total population poptot imposed.
c       This maintains the validity of the sum rules throughout the  
c	decomposition process.
c
c       A novel feature is that the diagonal elements of L (the pivot
c       elements) are computed using sum rules, which eliminates
c       cancellation problems, maximizes the robustness of the method,
c       and guarantees non-negative solutions.
c
	INTEGER n,ndim
	DOUBLE PRECISION gamma(ndim,ndim),poptot,pop(*)

	INTEGER i,j,k
	DOUBLE PRECISION lambda,sum,const
	DOUBLE PRECISION ZERO,ONE
	PARAMETER (ZERO=0.d0, ONE=1.d0)

	INTEGER NLOC
	PARAMETER (NLOC=200)
	DOUBLE PRECISION lmat(NLOC,NLOC),umat(NLOC,NLOC),x(NLOC)

c       [Note: The local arrays lmat, umat, and x have been introduced for
c	clarity, and to avoid altering the original matrix gamma.  If
c	one does not care about preserving gamma and wants to conserve
c	storage, the whole calculation can be done within gamma and pop, 
c	without defining these local arrays.]

c	   Check if the dimension of local arrays NLOC is large enough

	if(n.gt.NLOC)then
	   stop 'subroutine rate_eqs: NLOC too small'
	endif

c          Initialize lmat with the non-diagonal elements of gamma
c            (lmat will eventually become the lower triangular matrix L.)
c          Initialize umat to the unit matrix.
c            (umat will eventually become the upper triangular matrix U.)
c
	do i=1,n
	   do j=1,n
	      lmat(i,j)=gamma(i,j)
	      umat(i,j)=ZERO
	   enddo
	   umat(i,i)=ONE
	enddo
c
c	   The diagonal elements of lmat should be equal to minus the sum
c	   over the non-diagonal elements in the same column. We only 
c          need the (1,1) element, for initialization, so we compute this
c          here just to be sure.
c
	
        sum=ZERO
	do i=2,n
	      sum=sum+lmat(i,1)
        enddo
	lmat(1,1)=-sum

c
c          Now perform an LU decomposition of the matrix. No "pivoting"
c	   is done, because the problem does not require it.
c
	do k=1,n-1
	   do j=k+1,n
	      lambda=lmat(k,j)/lmat(k,k)
	      lmat(k,j)=ZERO
	      umat(k,j)=lambda
              sum=ZERO
	      do i=k+1,n
		 if(i.ne.j)then   !determine nondiagonal elements
                    lmat(i,j)=lmat(i,j)-lambda*lmat(i,k)
                    sum=sum+lmat(i,j)
                 endif
	      enddo
              lmat(j,j)=-sum   !diagonal element determined by sum rule
	   enddo
	enddo

c
c          The statistical equilibrium equations are homogeneous, so the
c          solution is arbitrary to within a scaling factor.  We first
c          find the solution x for which the last element x(n) is unity.
c
	x(n)=ONE
	do i=n-1,1,-1
	   sum=ZERO
	   do j=i+1,n
	      sum=sum+umat(i,j)*x(j)
	   enddo
	   x(i)=-sum
	enddo
c
c           The desired solution for pop is proportional to x.  We  
c           determine the proportionality constant const with the
c           condition of fixed total population poptot.
c
	sum=ZERO
	do j=1,n
	   sum=sum+x(j)
	enddo

	const=poptot/sum

	do i=1,n
	   pop(i)=const*x(i)
	enddo

	return
	end
