!*==clalsd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLALSD uses the singular value decomposition of A to solve the least squares problem.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLALSD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clalsd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clalsd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clalsd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND,
!                          RANK, WORK, RWORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ
!       REAL               RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               D( * ), E( * ), RWORK( * )
!       COMPLEX            B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLALSD uses the singular value decomposition of A to solve the least
!> squares problem of finding X to minimize the Euclidean norm of each
!> column of A*X-B, where A is N-by-N upper bidiagonal, and X and B
!> are N-by-NRHS. The solution X overwrites B.
!>
!> The singular values of A smaller than RCOND times the largest
!> singular value are treated as zero in solving the least squares
!> problem; in this case a minimum norm solution is returned.
!> The actual singular values are returned in D in ascending order.
!>
!> This code makes very mild assumptions about floating point
!> arithmetic. It will work on machines with a guard digit in
!> add/subtract, or on those binary machines without guard digits
!> which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
!> It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>         = 'U': D and E define an upper bidiagonal matrix.
!>         = 'L': D and E define a  lower bidiagonal matrix.
!> \endverbatim
!>
!> \param[in] SMLSIZ
!> \verbatim
!>          SMLSIZ is INTEGER
!>         The maximum size of the subproblems at the bottom of the
!>         computation tree.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the  bidiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>         The number of columns of B. NRHS must be at least 1.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>         On entry D contains the main diagonal of the bidiagonal
!>         matrix. On exit, if INFO = 0, D contains its singular values.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>         Contains the super-diagonal entries of the bidiagonal matrix.
!>         On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>         On input, B contains the right hand sides of the least
!>         squares problem. On output, B contains the solution X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>         The leading dimension of B in the calling subprogram.
!>         LDB must be at least max(1,N).
!> \endverbatim
!>
!> \param[in] RCOND
!> \verbatim
!>          RCOND is REAL
!>         The singular values of A less than or equal to RCOND times
!>         the largest singular value are treated as zero in solving
!>         the least squares problem. If RCOND is negative,
!>         machine precision is used instead.
!>         For example, if diag(S)*X=B were the least squares problem,
!>         where diag(S) is a diagonal matrix of singular values, the
!>         solution would be X(i) = B(i) / S(i) if S(i) is greater than
!>         RCOND*max(S), and X(i) = 0 if S(i) is less than or equal to
!>         RCOND*max(S).
!> \endverbatim
!>
!> \param[out] RANK
!> \verbatim
!>          RANK is INTEGER
!>         The number of singular values of A greater than RCOND times
!>         the largest singular value.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N * NRHS).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension at least
!>         (9*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS +
!>         MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS ),
!>         where
!>         NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (3*N*NLVL + 11*N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>         = 0:  successful exit.
!>         < 0:  if INFO = -i, the i-th argument had an illegal value.
!>         > 0:  The algorithm failed to compute a singular value while
!>               working on the submatrix lying in rows and columns
!>               INFO/(N+1) through MOD(INFO,N+1).
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
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!>       California at Berkeley, USA \n
!>     Osni Marques, LBNL/NERSC, USA \n
!
!  =====================================================================
      SUBROUTINE CLALSD(Uplo,Smlsiz,N,Nrhs,D,E,B,Ldb,Rcond,Rank,Work,   &
     &                  Rwork,Iwork,Info)
      USE S_CCOPY
      USE S_CLACPY
      USE S_CLALSA
      USE S_CLASCL
      USE S_CLASET
      USE S_CSROT
      USE S_ISAMAX
      USE S_SGEMM
      USE S_SLAMCH
      USE S_SLANST
      USE S_SLARTG
      USE S_SLASCL
      USE S_SLASDA
      USE S_SLASDQ
      USE S_SLASET
      USE S_SLASRT
      USE S_XERBLA
      IMPLICIT NONE
!*--CLALSD207
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER , INTENT(IN) :: Uplo
      INTEGER :: Smlsiz
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      COMPLEX , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: bx , bxst , c , difl , difr , givcol , givnum ,        &
     &           givptr , i , icmpq1 , icmpq2 , irwb , irwib , irwrb ,  &
     &           irwu , irwvt , irwwrk , iwk , j , jcol , jimag ,       &
     &           jreal , jrow , k , nlvl , nm1 , nrwork , nsize , nsub ,&
     &           perm , poles , s , sizei , smlszp , sqre , st , st1 ,  &
     &           u , vt , z
      REAL :: cs , eps , orgnrm , r , rcnd , sn , tol
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Nrhs<1 ) THEN
         Info = -4
      ELSEIF ( (Ldb<1) .OR. (Ldb<N) ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CLALSD',-Info)
         RETURN
      ENDIF
!
      eps = SLAMCH('Epsilon')
!
!     Set up the tolerance.
!
      IF ( (Rcond<=ZERO) .OR. (Rcond>=ONE) ) THEN
         rcnd = eps
      ELSE
         rcnd = Rcond
      ENDIF
!
      Rank = 0
!
!     Quick return if possible.
!
      IF ( N==0 ) THEN
         RETURN
      ELSEIF ( N==1 ) THEN
         IF ( D(1)==ZERO ) THEN
            CALL CLASET('A',1,Nrhs,CZERO,CZERO,B,Ldb)
         ELSE
            Rank = 1
            CALL CLASCL('G',0,0,D(1),ONE,1,Nrhs,B,Ldb,Info)
            D(1) = ABS(D(1))
         ENDIF
         RETURN
      ENDIF
!
!     Rotate the matrix if it is lower bidiagonal.
!
      IF ( Uplo=='L' ) THEN
         DO i = 1 , N - 1
            CALL SLARTG(D(i),E(i),cs,sn,r)
            D(i) = r
            E(i) = sn*D(i+1)
            D(i+1) = cs*D(i+1)
            IF ( Nrhs==1 ) THEN
               CALL CSROT(1,B(i,1),1,B(i+1,1),1,cs,sn)
            ELSE
               Rwork(i*2-1) = cs
               Rwork(i*2) = sn
            ENDIF
         ENDDO
         IF ( Nrhs>1 ) THEN
            DO i = 1 , Nrhs
               DO j = 1 , N - 1
                  cs = Rwork(j*2-1)
                  sn = Rwork(j*2)
                  CALL CSROT(1,B(j,i),1,B(j+1,i),1,cs,sn)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!
!     Scale.
!
      nm1 = N - 1
      orgnrm = SLANST('M',N,D,E)
      IF ( orgnrm==ZERO ) THEN
         CALL CLASET('A',N,Nrhs,CZERO,CZERO,B,Ldb)
         RETURN
      ENDIF
!
      CALL SLASCL('G',0,0,orgnrm,ONE,N,1,D,N,Info)
      CALL SLASCL('G',0,0,orgnrm,ONE,nm1,1,E,nm1,Info)
!
!     If N is smaller than the minimum divide size SMLSIZ, then solve
!     the problem with another solver.
!
      IF ( N<=Smlsiz ) THEN
         irwu = 1
         irwvt = irwu + N*N
         irwwrk = irwvt + N*N
         irwrb = irwwrk
         irwib = irwrb + N*Nrhs
         irwb = irwib + N*Nrhs
         CALL SLASET('A',N,N,ZERO,ONE,Rwork(irwu),N)
         CALL SLASET('A',N,N,ZERO,ONE,Rwork(irwvt),N)
         CALL SLASDQ('U',0,N,N,N,0,D,E,Rwork(irwvt),N,Rwork(irwu),N,    &
     &               Rwork(irwwrk),1,Rwork(irwwrk),Info)
         IF ( Info/=0 ) RETURN
!
!        In the real version, B is passed to SLASDQ and multiplied
!        internally by Q**H. Here B is complex and that product is
!        computed below in two steps (real and imaginary parts).
!
         j = irwb - 1
         DO jcol = 1 , Nrhs
            DO jrow = 1 , N
               j = j + 1
               Rwork(j) = REAL(B(jrow,jcol))
            ENDDO
         ENDDO
         CALL SGEMM('T','N',N,Nrhs,N,ONE,Rwork(irwu),N,Rwork(irwb),N,   &
     &              ZERO,Rwork(irwrb),N)
         j = irwb - 1
         DO jcol = 1 , Nrhs
            DO jrow = 1 , N
               j = j + 1
               Rwork(j) = AIMAG(B(jrow,jcol))
            ENDDO
         ENDDO
         CALL SGEMM('T','N',N,Nrhs,N,ONE,Rwork(irwu),N,Rwork(irwb),N,   &
     &              ZERO,Rwork(irwib),N)
         jreal = irwrb - 1
         jimag = irwib - 1
         DO jcol = 1 , Nrhs
            DO jrow = 1 , N
               jreal = jreal + 1
               jimag = jimag + 1
               B(jrow,jcol) = CMPLX(Rwork(jreal),Rwork(jimag))
            ENDDO
         ENDDO
!
         tol = rcnd*ABS(D(ISAMAX(N,D,1)))
         DO i = 1 , N
            IF ( D(i)<=tol ) THEN
               CALL CLASET('A',1,Nrhs,CZERO,CZERO,B(i,1),Ldb)
            ELSE
               CALL CLASCL('G',0,0,D(i),ONE,1,Nrhs,B(i,1),Ldb,Info)
               Rank = Rank + 1
            ENDIF
         ENDDO
!
!        Since B is complex, the following call to SGEMM is performed
!        in two steps (real and imaginary parts). That is for V * B
!        (in the real version of the code V**H is stored in WORK).
!
!        CALL SGEMM( 'T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO,
!    $               WORK( NWORK ), N )
!
         j = irwb - 1
         DO jcol = 1 , Nrhs
            DO jrow = 1 , N
               j = j + 1
               Rwork(j) = REAL(B(jrow,jcol))
            ENDDO
         ENDDO
         CALL SGEMM('T','N',N,Nrhs,N,ONE,Rwork(irwvt),N,Rwork(irwb),N,  &
     &              ZERO,Rwork(irwrb),N)
         j = irwb - 1
         DO jcol = 1 , Nrhs
            DO jrow = 1 , N
               j = j + 1
               Rwork(j) = AIMAG(B(jrow,jcol))
            ENDDO
         ENDDO
         CALL SGEMM('T','N',N,Nrhs,N,ONE,Rwork(irwvt),N,Rwork(irwb),N,  &
     &              ZERO,Rwork(irwib),N)
         jreal = irwrb - 1
         jimag = irwib - 1
         DO jcol = 1 , Nrhs
            DO jrow = 1 , N
               jreal = jreal + 1
               jimag = jimag + 1
               B(jrow,jcol) = CMPLX(Rwork(jreal),Rwork(jimag))
            ENDDO
         ENDDO
!
!        Unscale.
!
         CALL SLASCL('G',0,0,ONE,orgnrm,N,1,D,N,Info)
         CALL SLASRT('D',N,D,Info)
         CALL CLASCL('G',0,0,orgnrm,ONE,N,Nrhs,B,Ldb,Info)
!
         RETURN
      ENDIF
!
!     Book-keeping and setting up some constants.
!
      nlvl = INT(LOG(REAL(N)/REAL(Smlsiz+1))/LOG(TWO)) + 1
!
      smlszp = Smlsiz + 1
!
      u = 1
      vt = 1 + Smlsiz*N
      difl = vt + smlszp*N
      difr = difl + nlvl*N
      z = difr + nlvl*N*2
      c = z + nlvl*N
      s = c + N
      poles = s + N
      givnum = poles + 2*nlvl*N
      nrwork = givnum + 2*nlvl*N
      bx = 1
!
      irwrb = nrwork
      irwib = irwrb + Smlsiz*Nrhs
      irwb = irwib + Smlsiz*Nrhs
!
      sizei = 1 + N
      k = sizei + N
      givptr = k + N
      perm = givptr + N
      givcol = perm + nlvl*N
      iwk = givcol + nlvl*N*2
!
      st = 1
      sqre = 0
      icmpq1 = 1
      icmpq2 = 0
      nsub = 0
!
      DO i = 1 , N
         IF ( ABS(D(i))<eps ) D(i) = SIGN(eps,D(i))
      ENDDO
!
      DO i = 1 , nm1
         IF ( (ABS(E(i))<eps) .OR. (i==nm1) ) THEN
            nsub = nsub + 1
            Iwork(nsub) = st
!
!           Subproblem found. First determine its size and then
!           apply divide and conquer on it.
!
            IF ( i<nm1 ) THEN
!
!              A subproblem with E(I) small for I < NM1.
!
               nsize = i - st + 1
               Iwork(sizei+nsub-1) = nsize
            ELSEIF ( ABS(E(i))>=eps ) THEN
!
!              A subproblem with E(NM1) not too small but I = NM1.
!
               nsize = N - st + 1
               Iwork(sizei+nsub-1) = nsize
            ELSE
!
!              A subproblem with E(NM1) small. This implies an
!              1-by-1 subproblem at D(N), which is not solved
!              explicitly.
!
               nsize = i - st + 1
               Iwork(sizei+nsub-1) = nsize
               nsub = nsub + 1
               Iwork(nsub) = N
               Iwork(sizei+nsub-1) = 1
               CALL CCOPY(Nrhs,B(N,1),Ldb,Work(bx+nm1),N)
            ENDIF
            st1 = st - 1
            IF ( nsize==1 ) THEN
!
!              This is a 1-by-1 subproblem and is not solved
!              explicitly.
!
               CALL CCOPY(Nrhs,B(st,1),Ldb,Work(bx+st1),N)
            ELSEIF ( nsize<=Smlsiz ) THEN
!
!              This is a small subproblem and is solved by SLASDQ.
!
               CALL SLASET('A',nsize,nsize,ZERO,ONE,Rwork(vt+st1),N)
               CALL SLASET('A',nsize,nsize,ZERO,ONE,Rwork(u+st1),N)
               CALL SLASDQ('U',0,nsize,nsize,nsize,0,D(st),E(st),       &
     &                     Rwork(vt+st1),N,Rwork(u+st1),N,Rwork(nrwork),&
     &                     1,Rwork(nrwork),Info)
               IF ( Info/=0 ) RETURN
!
!              In the real version, B is passed to SLASDQ and multiplied
!              internally by Q**H. Here B is complex and that product is
!              computed below in two steps (real and imaginary parts).
!
               j = irwb - 1
               DO jcol = 1 , Nrhs
                  DO jrow = st , st + nsize - 1
                     j = j + 1
                     Rwork(j) = REAL(B(jrow,jcol))
                  ENDDO
               ENDDO
               CALL SGEMM('T','N',nsize,Nrhs,nsize,ONE,Rwork(u+st1),N,  &
     &                    Rwork(irwb),nsize,ZERO,Rwork(irwrb),nsize)
               j = irwb - 1
               DO jcol = 1 , Nrhs
                  DO jrow = st , st + nsize - 1
                     j = j + 1
                     Rwork(j) = AIMAG(B(jrow,jcol))
                  ENDDO
               ENDDO
               CALL SGEMM('T','N',nsize,Nrhs,nsize,ONE,Rwork(u+st1),N,  &
     &                    Rwork(irwb),nsize,ZERO,Rwork(irwib),nsize)
               jreal = irwrb - 1
               jimag = irwib - 1
               DO jcol = 1 , Nrhs
                  DO jrow = st , st + nsize - 1
                     jreal = jreal + 1
                     jimag = jimag + 1
                     B(jrow,jcol) = CMPLX(Rwork(jreal),Rwork(jimag))
                  ENDDO
               ENDDO
!
               CALL CLACPY('A',nsize,Nrhs,B(st,1),Ldb,Work(bx+st1),N)
            ELSE
!
!              A large problem. Solve it using divide and conquer.
!
               CALL SLASDA(icmpq1,Smlsiz,nsize,sqre,D(st),E(st),        &
     &                     Rwork(u+st1),N,Rwork(vt+st1),Iwork(k+st1),   &
     &                     Rwork(difl+st1),Rwork(difr+st1),Rwork(z+st1),&
     &                     Rwork(poles+st1),Iwork(givptr+st1),          &
     &                     Iwork(givcol+st1),N,Iwork(perm+st1),         &
     &                     Rwork(givnum+st1),Rwork(c+st1),Rwork(s+st1), &
     &                     Rwork(nrwork),Iwork(iwk),Info)
               IF ( Info/=0 ) RETURN
               bxst = bx + st1
               CALL CLALSA(icmpq2,Smlsiz,nsize,Nrhs,B(st,1),Ldb,        &
     &                     Work(bxst),N,Rwork(u+st1),N,Rwork(vt+st1),   &
     &                     Iwork(k+st1),Rwork(difl+st1),Rwork(difr+st1),&
     &                     Rwork(z+st1),Rwork(poles+st1),               &
     &                     Iwork(givptr+st1),Iwork(givcol+st1),N,       &
     &                     Iwork(perm+st1),Rwork(givnum+st1),           &
     &                     Rwork(c+st1),Rwork(s+st1),Rwork(nrwork),     &
     &                     Iwork(iwk),Info)
               IF ( Info/=0 ) RETURN
            ENDIF
            st = i + 1
         ENDIF
      ENDDO
!
!     Apply the singular values and treat the tiny ones as zero.
!
      tol = rcnd*ABS(D(ISAMAX(N,D,1)))
!
      DO i = 1 , N
!
!        Some of the elements in D can be negative because 1-by-1
!        subproblems were not solved explicitly.
!
         IF ( ABS(D(i))<=tol ) THEN
            CALL CLASET('A',1,Nrhs,CZERO,CZERO,Work(bx+i-1),N)
         ELSE
            Rank = Rank + 1
            CALL CLASCL('G',0,0,D(i),ONE,1,Nrhs,Work(bx+i-1),N,Info)
         ENDIF
         D(i) = ABS(D(i))
      ENDDO
!
!     Now apply back the right singular vectors.
!
      icmpq2 = 1
      DO i = 1 , nsub
         st = Iwork(i)
         st1 = st - 1
         nsize = Iwork(sizei+i-1)
         bxst = bx + st1
         IF ( nsize==1 ) THEN
            CALL CCOPY(Nrhs,Work(bxst),N,B(st,1),Ldb)
         ELSEIF ( nsize<=Smlsiz ) THEN
!
!           Since B and BX are complex, the following call to SGEMM
!           is performed in two steps (real and imaginary parts).
!
!           CALL SGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, ONE,
!    $                  RWORK( VT+ST1 ), N, RWORK( BXST ), N, ZERO,
!    $                  B( ST, 1 ), LDB )
!
            j = bxst - N - 1
            jreal = irwb - 1
            DO jcol = 1 , Nrhs
               j = j + N
               DO jrow = 1 , nsize
                  jreal = jreal + 1
                  Rwork(jreal) = REAL(Work(j+jrow))
               ENDDO
            ENDDO
            CALL SGEMM('T','N',nsize,Nrhs,nsize,ONE,Rwork(vt+st1),N,    &
     &                 Rwork(irwb),nsize,ZERO,Rwork(irwrb),nsize)
            j = bxst - N - 1
            jimag = irwb - 1
            DO jcol = 1 , Nrhs
               j = j + N
               DO jrow = 1 , nsize
                  jimag = jimag + 1
                  Rwork(jimag) = AIMAG(Work(j+jrow))
               ENDDO
            ENDDO
            CALL SGEMM('T','N',nsize,Nrhs,nsize,ONE,Rwork(vt+st1),N,    &
     &                 Rwork(irwb),nsize,ZERO,Rwork(irwib),nsize)
            jreal = irwrb - 1
            jimag = irwib - 1
            DO jcol = 1 , Nrhs
               DO jrow = st , st + nsize - 1
                  jreal = jreal + 1
                  jimag = jimag + 1
                  B(jrow,jcol) = CMPLX(Rwork(jreal),Rwork(jimag))
               ENDDO
            ENDDO
         ELSE
            CALL CLALSA(icmpq2,Smlsiz,nsize,Nrhs,Work(bxst),N,B(st,1),  &
     &                  Ldb,Rwork(u+st1),N,Rwork(vt+st1),Iwork(k+st1),  &
     &                  Rwork(difl+st1),Rwork(difr+st1),Rwork(z+st1),   &
     &                  Rwork(poles+st1),Iwork(givptr+st1),             &
     &                  Iwork(givcol+st1),N,Iwork(perm+st1),            &
     &                  Rwork(givnum+st1),Rwork(c+st1),Rwork(s+st1),    &
     &                  Rwork(nrwork),Iwork(iwk),Info)
            IF ( Info/=0 ) RETURN
         ENDIF
      ENDDO
!
!     Unscale and sort the singular values.
!
      CALL SLASCL('G',0,0,ONE,orgnrm,N,1,D,N,Info)
      CALL SLASRT('D',N,D,Info)
      CALL CLASCL('G',0,0,orgnrm,ONE,N,Nrhs,B,Ldb,Info)
!
!
!     End of CLALSD
!
      END SUBROUTINE CLALSD
