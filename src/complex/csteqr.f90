!*==csteqr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CSTEQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSTEQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csteqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csteqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csteqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPZ
!       INTEGER            INFO, LDZ, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), WORK( * )
!       COMPLEX            Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!> symmetric tridiagonal matrix using the implicit QL or QR method.
!> The eigenvectors of a full or band complex Hermitian matrix can also
!> be found if CHETRD or CHPTRD or CHBTRD has been used to reduce this
!> matrix to tridiagonal form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only.
!>          = 'V':  Compute eigenvalues and eigenvectors of the original
!>                  Hermitian matrix.  On entry, Z must contain the
!>                  unitary matrix used to reduce the original matrix
!>                  to tridiagonal form.
!>          = 'I':  Compute eigenvalues and eigenvectors of the
!>                  tridiagonal matrix.  Z is initialized to the identity
!>                  matrix.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          On entry, the diagonal elements of the tridiagonal matrix.
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
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, N)
!>          On entry, if  COMPZ = 'V', then Z contains the unitary
!>          matrix used in the reduction to tridiagonal form.
!>          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
!>          orthonormal eigenvectors of the original Hermitian matrix,
!>          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!>          of the symmetric tridiagonal matrix.
!>          If COMPZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          eigenvectors are desired, then  LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (max(1,2*N-2))
!>          If COMPZ = 'N', then WORK is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  the algorithm has failed to find all the eigenvalues in
!>                a total of 30*N iterations; if INFO = i, then i
!>                elements of E have not converged to zero; on exit, D
!>                and E contain the elements of a symmetric tridiagonal
!>                matrix which is unitarily similar to the original
!>                matrix.
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
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CSTEQR(Compz,N,D,E,Z,Ldz,Work,Info)
      IMPLICIT NONE
!*--CSTEQR136
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Compz
      INTEGER Info , Ldz , N
!     ..
!     .. Array Arguments ..
      REAL D(*) , E(*) , Work(*)
      COMPLEX Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO , THREE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,TWO=2.0E0,THREE=3.0E0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E0,0.0E0),CONE=(1.0E0,0.0E0))
      INTEGER MAXIT
      PARAMETER (MAXIT=30)
!     ..
!     .. Local Scalars ..
      INTEGER i , icompz , ii , iscale , j , jtot , k , l , l1 , lend , &
     &        lendm1 , lendp1 , lendsv , lm1 , lsv , m , mm , mm1 ,     &
     &        nm1 , nmaxit
      REAL anorm , b , c , eps , eps2 , f , g , p , r , rt1 , rt2 , s , &
     &     safmax , safmin , ssfmax , ssfmin , tst
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLAMCH , SLANST , SLAPY2
      EXTERNAL LSAME , SLAMCH , SLANST , SLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL CLASET , CLASR , CSWAP , SLAE2 , SLAEV2 , SLARTG ,       &
     &         SLASCL , SLASRT , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SIGN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( LSAME(Compz,'N') ) THEN
         icompz = 0
      ELSEIF ( LSAME(Compz,'V') ) THEN
         icompz = 1
      ELSEIF ( LSAME(Compz,'I') ) THEN
         icompz = 2
      ELSE
         icompz = -1
      ENDIF
      IF ( icompz<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( (Ldz<1) .OR. (icompz>0 .AND. Ldz<MAX(1,N)) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CSTEQR',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( N==1 ) THEN
         IF ( icompz==2 ) Z(1,1) = CONE
         RETURN
      ENDIF
!
!     Determine the unit roundoff and over/underflow thresholds.
!
      eps = SLAMCH('E')
      eps2 = eps**2
      safmin = SLAMCH('S')
      safmax = ONE/safmin
      ssfmax = SQRT(safmax)/THREE
      ssfmin = SQRT(safmin)/eps2
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
      IF ( icompz==2 ) CALL CLASET('Full',N,N,CZERO,CONE,Z,Ldz)
!
      nmaxit = N*MAXIT
      jtot = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      l1 = 1
      nm1 = N - 1
!
      DO WHILE ( l1<=N )
         IF ( l1>1 ) E(l1-1) = ZERO
         IF ( l1<=nm1 ) THEN
            DO m = l1 , nm1
               tst = ABS(E(m))
               IF ( tst==ZERO ) GOTO 50
               IF ( tst<=(SQRT(ABS(D(m)))*SQRT(ABS(D(m+1))))*eps ) THEN
                  E(m) = ZERO
                  GOTO 50
               ENDIF
            ENDDO
         ENDIF
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
            anorm = SLANST('I',lend-l+1,D(l),E(l))
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
!     Choose between QL and QR iteration
!
               IF ( ABS(D(lend))<ABS(D(l)) ) THEN
                  lend = lsv
                  l = lendsv
               ENDIF
!
               IF ( lend>l ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
 55               IF ( l/=lend ) THEN
                     lendm1 = lend - 1
                     DO m = l , lendm1
                        tst = ABS(E(m))**2
                        IF ( tst<=(eps2*ABS(D(m)))*ABS(D(m+1))+safmin ) &
     &                       GOTO 60
                     ENDDO
                  ENDIF
!
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
!        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
!        to compute its eigensystem.
!
                     IF ( m==l+1 ) THEN
                        IF ( icompz>0 ) THEN
                           CALL SLAEV2(D(l),E(l),D(l+1),rt1,rt2,c,s)
                           Work(l) = c
                           Work(N-1+l) = s
                           CALL CLASR('R','V','B',N,2,Work(l),          &
     &                                Work(N-1+l),Z(1,l),Ldz)
                        ELSE
                           CALL SLAE2(D(l),E(l),D(l+1),rt1,rt2)
                        ENDIF
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
                        g = (D(l+1)-p)/(TWO*E(l))
                        r = SLAPY2(g,ONE)
                        g = D(m) - p + (E(l)/(g+SIGN(r,g)))
!
                        s = ONE
                        c = ONE
                        p = ZERO
!
!        Inner loop
!
                        mm1 = m - 1
                        DO i = mm1 , l , -1
                           f = s*E(i)
                           b = c*E(i)
                           CALL SLARTG(g,f,c,s,r)
                           IF ( i/=m-1 ) E(i+1) = r
                           g = D(i+1) - p
                           r = (D(i)-g)*s + TWO*c*b
                           p = s*r
                           D(i+1) = g + p
                           g = c*r - b
!
!           If eigenvectors are desired, then save rotations.
!
                           IF ( icompz>0 ) THEN
                              Work(i) = c
                              Work(N-1+i) = -s
                           ENDIF
!
                        ENDDO
!
!        If eigenvectors are desired, then apply saved rotations.
!
                        IF ( icompz>0 ) THEN
                           mm = m - l + 1
                           CALL CLASR('R','V','B',N,mm,Work(l),         &
     &                                Work(N-1+l),Z(1,l),Ldz)
                        ENDIF
!
                        D(l) = D(l) - p
                        E(l) = g
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
 65               IF ( l/=lend ) THEN
                     lendp1 = lend + 1
                     DO m = l , lendp1 , -1
                        tst = ABS(E(m-1))**2
                        IF ( tst<=(eps2*ABS(D(m)))*ABS(D(m-1))+safmin ) &
     &                       GOTO 70
                     ENDDO
                  ENDIF
!
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
!        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
!        to compute its eigensystem.
!
                     IF ( m==l-1 ) THEN
                        IF ( icompz>0 ) THEN
                           CALL SLAEV2(D(l-1),E(l-1),D(l),rt1,rt2,c,s)
                           Work(m) = c
                           Work(N-1+m) = s
                           CALL CLASR('R','V','F',N,2,Work(m),          &
     &                                Work(N-1+m),Z(1,l-1),Ldz)
                        ELSE
                           CALL SLAE2(D(l-1),E(l-1),D(l),rt1,rt2)
                        ENDIF
                        D(l-1) = rt1
                        D(l) = rt2
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
                        g = (D(l-1)-p)/(TWO*E(l-1))
                        r = SLAPY2(g,ONE)
                        g = D(m) - p + (E(l-1)/(g+SIGN(r,g)))
!
                        s = ONE
                        c = ONE
                        p = ZERO
!
!        Inner loop
!
                        lm1 = l - 1
                        DO i = m , lm1
                           f = s*E(i)
                           b = c*E(i)
                           CALL SLARTG(g,f,c,s,r)
                           IF ( i/=m ) E(i-1) = r
                           g = D(i) - p
                           r = (D(i+1)-g)*s + TWO*c*b
                           p = s*r
                           D(i) = g + p
                           g = c*r - b
!
!           If eigenvectors are desired, then save rotations.
!
                           IF ( icompz>0 ) THEN
                              Work(i) = c
                              Work(N-1+i) = s
                           ENDIF
!
                        ENDDO
!
!        If eigenvectors are desired, then apply saved rotations.
!
                        IF ( icompz>0 ) THEN
                           mm = l - m + 1
                           CALL CLASR('R','V','F',N,mm,Work(m),         &
     &                                Work(N-1+m),Z(1,m),Ldz)
                        ENDIF
!
                        D(l) = D(l) - p
                        E(lm1) = g
                        GOTO 65
                     ENDIF
                  ENDIF
!
               ENDIF
!
!     Undo scaling if necessary
!
 80            IF ( iscale==1 ) THEN
                  CALL SLASCL('G',0,0,ssfmax,anorm,lendsv-lsv+1,1,D(lsv)&
     &                        ,N,Info)
                  CALL SLASCL('G',0,0,ssfmax,anorm,lendsv-lsv,1,E(lsv), &
     &                        N,Info)
               ELSEIF ( iscale==2 ) THEN
                  CALL SLASCL('G',0,0,ssfmin,anorm,lendsv-lsv+1,1,D(lsv)&
     &                        ,N,Info)
                  CALL SLASCL('G',0,0,ssfmin,anorm,lendsv-lsv,1,E(lsv), &
     &                        N,Info)
               ENDIF
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
               IF ( jtot==nmaxit ) THEN
                  DO i = 1 , N - 1
                     IF ( E(i)/=ZERO ) Info = Info + 1
                  ENDDO
                  RETURN
               ENDIF
            ENDIF
         ENDIF
      ENDDO
!
!     Order eigenvalues and eigenvectors.
!
      IF ( icompz==0 ) THEN
!
!        Use Quick Sort
!
         CALL SLASRT('I',N,D,Info)
!
      ELSE
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
         DO ii = 2 , N
            i = ii - 1
            k = i
            p = D(i)
            DO j = ii , N
               IF ( D(j)<p ) THEN
                  k = j
                  p = D(j)
               ENDIF
            ENDDO
            IF ( k/=i ) THEN
               D(k) = D(i)
               D(i) = p
               CALL CSWAP(N,Z(1,i),1,Z(1,k),1)
            ENDIF
         ENDDO
      ENDIF
!
!     End of CSTEQR
!
      END SUBROUTINE CSTEQR
