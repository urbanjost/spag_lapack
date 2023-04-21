!*==ssterf.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SSTERF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSTERF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssterf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssterf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssterf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSTERF( N, D, E, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSTERF computes all eigenvalues of a symmetric tridiagonal matrix
!> using the Pal-Walker-Kahan variant of the QL or QR algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          On entry, the n diagonal elements of the tridiagonal matrix.
!>          On exit, if INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          On entry, the (n-1) subdiagonal elements of the tridiagonal
!>          matrix.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  the algorithm failed to find all of the eigenvalues in
!>                a total of 30*N iterations; if INFO = i, then i
!>                elements of E have not converged to zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SSTERF(N,D,E,Info)
      IMPLICIT NONE
!*--SSTERF90
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , N
!     ..
!     .. Array Arguments ..
      REAL D(*) , E(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO , THREE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,TWO=2.0E0,THREE=3.0E0)
      INTEGER MAXIT
      PARAMETER (MAXIT=30)
!     ..
!     .. Local Scalars ..
      INTEGER i , iscale , jtot , l , l1 , lend , lendsv , lsv , m ,    &
     &        nmaxit
      REAL alpha , anorm , bb , c , eps , eps2 , gamma , oldc , oldgam ,&
     &     p , r , rt1 , rt2 , rte , s , safmax , safmin , sigma ,      &
     &     ssfmax , ssfmin
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANST , SLAPY2
      EXTERNAL SLAMCH , SLANST , SLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL SLAE2 , SLASCL , SLASRT , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SIGN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
!     Quick return if possible
!
      IF ( N<0 ) THEN
         Info = -1
         CALL XERBLA('SSTERF',-Info)
         RETURN
      ENDIF
      IF ( N<=1 ) RETURN
!
!     Determine the unit roundoff for this environment.
!
      eps = SLAMCH('E')
      eps2 = eps**2
      safmin = SLAMCH('S')
      safmax = ONE/safmin
      ssfmax = SQRT(safmax)/THREE
      ssfmin = SQRT(safmin)/eps2
!
!     Compute the eigenvalues of the tridiagonal matrix.
!
      nmaxit = N*MAXIT
      sigma = ZERO
      jtot = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      l1 = 1
!
      DO WHILE ( l1<=N )
         IF ( l1>1 ) E(l1-1) = ZERO
         DO m = l1 , N - 1
            IF ( ABS(E(m))<=(SQRT(ABS(D(m)))*SQRT(ABS(D(m+1))))*eps )   &
     &           THEN
               E(m) = ZERO
               GOTO 50
            ENDIF
         ENDDO
         m = N
!
 50      l = l1
         lsv = l
         lend = m
         lendsv = lend
         l1 = m + 1
         IF ( lend/=l ) THEN
!
!     Scale submatrix in rows and columns L to LEND
!
            anorm = SLANST('M',lend-l+1,D(l),E(l))
            iscale = 0
            IF ( anorm/=ZERO ) THEN
               IF ( anorm>ssfmax ) THEN
                  iscale = 1
                  CALL SLASCL('G',0,0,anorm,ssfmax,lend-l+1,1,D(l),N,   &
     &                        Info)
                  CALL SLASCL('G',0,0,anorm,ssfmax,lend-l,1,E(l),N,Info)
               ELSEIF ( anorm<ssfmin ) THEN
                  iscale = 2
                  CALL SLASCL('G',0,0,anorm,ssfmin,lend-l+1,1,D(l),N,   &
     &                        Info)
                  CALL SLASCL('G',0,0,anorm,ssfmin,lend-l,1,E(l),N,Info)
               ENDIF
!
               DO i = l , lend - 1
                  E(i) = E(i)**2
               ENDDO
!
!     Choose between QL and QR iteration
!
               IF ( ABS(D(lend))<ABS(D(l)) ) THEN
                  lend = lsv
                  l = lendsv
               ENDIF
!
               IF ( lend>=l ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
 55               IF ( l/=lend ) THEN
                     DO m = l , lend - 1
                        IF ( ABS(E(m))<=eps2*ABS(D(m)*D(m+1)) ) GOTO 60
                     ENDDO
                  ENDIF
                  m = lend
!
 60               IF ( m<lend ) E(m) = ZERO
                  p = D(l)
                  IF ( m==l ) THEN
!
!        Eigenvalue found.
!
                     D(l) = p
!
                     l = l + 1
                     IF ( l<=lend ) GOTO 55
                  ELSE
!
!        If remaining matrix is 2 by 2, use SLAE2 to compute its
!        eigenvalues.
!
                     IF ( m==l+1 ) THEN
                        rte = SQRT(E(l))
                        CALL SLAE2(D(l),rte,D(l+1),rt1,rt2)
                        D(l) = rt1
                        D(l+1) = rt2
                        E(l) = ZERO
                        l = l + 2
                        IF ( l>lend ) GOTO 80
                        GOTO 55
                     ENDIF
!
                     IF ( jtot/=nmaxit ) THEN
                        jtot = jtot + 1
!
!        Form shift.
!
                        rte = SQRT(E(l))
                        sigma = (D(l+1)-p)/(TWO*rte)
                        r = SLAPY2(sigma,ONE)
                        sigma = p - (rte/(sigma+SIGN(r,sigma)))
!
                        c = ONE
                        s = ZERO
                        gamma = D(m) - sigma
                        p = gamma*gamma
!
!        Inner loop
!
                        DO i = m - 1 , l , -1
                           bb = E(i)
                           r = p + bb
                           IF ( i/=m-1 ) E(i+1) = s*r
                           oldc = c
                           c = p/r
                           s = bb/r
                           oldgam = gamma
                           alpha = D(i)
                           gamma = c*(alpha-sigma) - s*oldgam
                           D(i+1) = oldgam + (alpha-gamma)
                           IF ( c/=ZERO ) THEN
                              p = (gamma*gamma)/c
                           ELSE
                              p = oldc*bb
                           ENDIF
                        ENDDO
!
                        E(l) = s*p
                        D(l) = sigma + gamma
                        GOTO 55
                     ENDIF
                  ENDIF
!
               ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
 65               DO m = l , lend + 1 , -1
                     IF ( ABS(E(m-1))<=eps2*ABS(D(m)*D(m-1)) ) GOTO 70
                  ENDDO
                  m = lend
!
 70               IF ( m>lend ) E(m-1) = ZERO
                  p = D(l)
                  IF ( m==l ) THEN
!
!        Eigenvalue found.
!
                     D(l) = p
!
                     l = l - 1
                     IF ( l>=lend ) GOTO 65
                  ELSE
!
!        If remaining matrix is 2 by 2, use SLAE2 to compute its
!        eigenvalues.
!
                     IF ( m==l-1 ) THEN
                        rte = SQRT(E(l-1))
                        CALL SLAE2(D(l),rte,D(l-1),rt1,rt2)
                        D(l) = rt1
                        D(l-1) = rt2
                        E(l-1) = ZERO
                        l = l - 2
                        IF ( l<lend ) GOTO 80
                        GOTO 65
                     ENDIF
!
                     IF ( jtot/=nmaxit ) THEN
                        jtot = jtot + 1
!
!        Form shift.
!
                        rte = SQRT(E(l-1))
                        sigma = (D(l-1)-p)/(TWO*rte)
                        r = SLAPY2(sigma,ONE)
                        sigma = p - (rte/(sigma+SIGN(r,sigma)))
!
                        c = ONE
                        s = ZERO
                        gamma = D(m) - sigma
                        p = gamma*gamma
!
!        Inner loop
!
                        DO i = m , l - 1
                           bb = E(i)
                           r = p + bb
                           IF ( i/=m ) E(i-1) = s*r
                           oldc = c
                           c = p/r
                           s = bb/r
                           oldgam = gamma
                           alpha = D(i+1)
                           gamma = c*(alpha-sigma) - s*oldgam
                           D(i) = oldgam + (alpha-gamma)
                           IF ( c/=ZERO ) THEN
                              p = (gamma*gamma)/c
                           ELSE
                              p = oldc*bb
                           ENDIF
                        ENDDO
!
                        E(l-1) = s*p
                        D(l) = sigma + gamma
                        GOTO 65
                     ENDIF
                  ENDIF
!
               ENDIF
!
!     Undo scaling if necessary
!
 80            IF ( iscale==1 )                                         &
     &              CALL SLASCL('G',0,0,ssfmax,anorm,lendsv-lsv+1,1,    &
     &              D(lsv),N,Info)
               IF ( iscale==2 )                                         &
     &              CALL SLASCL('G',0,0,ssfmin,anorm,lendsv-lsv+1,1,    &
     &              D(lsv),N,Info)
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
               IF ( jtot>=nmaxit ) THEN
                  DO i = 1 , N - 1
                     IF ( E(i)/=ZERO ) Info = Info + 1
                  ENDDO
                  GOTO 99999
               ENDIF
            ENDIF
         ENDIF
      ENDDO
!
!     Sort eigenvalues in increasing order.
!
      CALL SLASRT('I',N,D,Info)
!
!
!     End of SSTERF
!
99999 END SUBROUTINE SSTERF
