!*==dlalsd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLALSD uses the singular value decomposition of A to solve the least squares problem.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLALSD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlalsd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlalsd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlalsd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND,
!                          RANK, WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ
!       DOUBLE PRECISION   RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   B( LDB, * ), D( * ), E( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLALSD uses the singular value decomposition of A to solve the least
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
!>          D is DOUBLE PRECISION array, dimension (N)
!>         On entry D contains the main diagonal of the bidiagonal
!>         matrix. On exit, if INFO = 0, D contains its singular values.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>         Contains the super-diagonal entries of the bidiagonal matrix.
!>         On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
!>          RCOND is DOUBLE PRECISION
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
!>          WORK is DOUBLE PRECISION array, dimension at least
!>         (9*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2),
!>         where NLVL = max(0, INT(log_2 (N/(SMLSIZ+1))) + 1).
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension at least
!>         (3*N*NLVL + 11*N)
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
!> \ingroup doubleOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!>       California at Berkeley, USA \n
!>     Osni Marques, LBNL/NERSC, USA \n
!
!  =====================================================================
      SUBROUTINE DLALSD(Uplo,Smlsiz,N,Nrhs,D,E,B,Ldb,Rcond,Rank,Work,   &
     &                  Iwork,Info)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DGEMM
      USE S_DLACPY
      USE S_DLALSA
      USE S_DLAMCH
      USE S_DLANST
      USE S_DLARTG
      USE S_DLASCL
      USE S_DLASDA
      USE S_DLASDQ
      USE S_DLASET
      USE S_DLASRT
      USE S_DROT
      USE S_IDAMAX
      USE S_XERBLA
      IMPLICIT NONE
!*--DLALSD199
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER , INTENT(IN) :: Uplo
      INTEGER :: Smlsiz
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: bx , bxst , c , difl , difr , givcol , givnum ,        &
     &           givptr , i , icmpq1 , icmpq2 , iwk , j , k , nlvl ,    &
     &           nm1 , nsize , nsub , nwork , perm , poles , s , sizei ,&
     &           smlszp , sqre , st , st1 , u , vt , z
      REAL(R8KIND) :: cs , eps , orgnrm , r , rcnd , sn , tol
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
         CALL XERBLA('DLALSD',-Info)
         RETURN
      ENDIF
!
      eps = DLAMCH('Epsilon')
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
            CALL DLASET('A',1,Nrhs,ZERO,ZERO,B,Ldb)
         ELSE
            Rank = 1
            CALL DLASCL('G',0,0,D(1),ONE,1,Nrhs,B,Ldb,Info)
            D(1) = ABS(D(1))
         ENDIF
         RETURN
      ENDIF
!
!     Rotate the matrix if it is lower bidiagonal.
!
      IF ( Uplo=='L' ) THEN
         DO i = 1 , N - 1
            CALL DLARTG(D(i),E(i),cs,sn,r)
            D(i) = r
            E(i) = sn*D(i+1)
            D(i+1) = cs*D(i+1)
            IF ( Nrhs==1 ) THEN
               CALL DROT(1,B(i,1),1,B(i+1,1),1,cs,sn)
            ELSE
               Work(i*2-1) = cs
               Work(i*2) = sn
            ENDIF
         ENDDO
         IF ( Nrhs>1 ) THEN
            DO i = 1 , Nrhs
               DO j = 1 , N - 1
                  cs = Work(j*2-1)
                  sn = Work(j*2)
                  CALL DROT(1,B(j,i),1,B(j+1,i),1,cs,sn)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!
!     Scale.
!
      nm1 = N - 1
      orgnrm = DLANST('M',N,D,E)
      IF ( orgnrm==ZERO ) THEN
         CALL DLASET('A',N,Nrhs,ZERO,ZERO,B,Ldb)
         RETURN
      ENDIF
!
      CALL DLASCL('G',0,0,orgnrm,ONE,N,1,D,N,Info)
      CALL DLASCL('G',0,0,orgnrm,ONE,nm1,1,E,nm1,Info)
!
!     If N is smaller than the minimum divide size SMLSIZ, then solve
!     the problem with another solver.
!
      IF ( N<=Smlsiz ) THEN
         nwork = 1 + N*N
         CALL DLASET('A',N,N,ZERO,ONE,Work,N)
         CALL DLASDQ('U',0,N,N,0,Nrhs,D,E,Work,N,Work,N,B,Ldb,          &
     &               Work(nwork),Info)
         IF ( Info/=0 ) RETURN
         tol = rcnd*ABS(D(IDAMAX(N,D,1)))
         DO i = 1 , N
            IF ( D(i)<=tol ) THEN
               CALL DLASET('A',1,Nrhs,ZERO,ZERO,B(i,1),Ldb)
            ELSE
               CALL DLASCL('G',0,0,D(i),ONE,1,Nrhs,B(i,1),Ldb,Info)
               Rank = Rank + 1
            ENDIF
         ENDDO
         CALL DGEMM('T','N',N,Nrhs,N,ONE,Work,N,B,Ldb,ZERO,Work(nwork), &
     &              N)
         CALL DLACPY('A',N,Nrhs,Work(nwork),N,B,Ldb)
!
!        Unscale.
!
         CALL DLASCL('G',0,0,ONE,orgnrm,N,1,D,N,Info)
         CALL DLASRT('D',N,D,Info)
         CALL DLASCL('G',0,0,orgnrm,ONE,N,Nrhs,B,Ldb,Info)
!
         RETURN
      ENDIF
!
!     Book-keeping and setting up some constants.
!
      nlvl = INT(LOG(DBLE(N)/DBLE(Smlsiz+1))/LOG(TWO)) + 1
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
      bx = givnum + 2*nlvl*N
      nwork = bx + N*Nrhs
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
               CALL DCOPY(Nrhs,B(N,1),Ldb,Work(bx+nm1),N)
            ENDIF
            st1 = st - 1
            IF ( nsize==1 ) THEN
!
!              This is a 1-by-1 subproblem and is not solved
!              explicitly.
!
               CALL DCOPY(Nrhs,B(st,1),Ldb,Work(bx+st1),N)
            ELSEIF ( nsize<=Smlsiz ) THEN
!
!              This is a small subproblem and is solved by DLASDQ.
!
               CALL DLASET('A',nsize,nsize,ZERO,ONE,Work(vt+st1),N)
               CALL DLASDQ('U',0,nsize,nsize,0,Nrhs,D(st),E(st),        &
     &                     Work(vt+st1),N,Work(nwork),N,B(st,1),Ldb,    &
     &                     Work(nwork),Info)
               IF ( Info/=0 ) RETURN
               CALL DLACPY('A',nsize,Nrhs,B(st,1),Ldb,Work(bx+st1),N)
            ELSE
!
!              A large problem. Solve it using divide and conquer.
!
               CALL DLASDA(icmpq1,Smlsiz,nsize,sqre,D(st),E(st),        &
     &                     Work(u+st1),N,Work(vt+st1),Iwork(k+st1),     &
     &                     Work(difl+st1),Work(difr+st1),Work(z+st1),   &
     &                     Work(poles+st1),Iwork(givptr+st1),           &
     &                     Iwork(givcol+st1),N,Iwork(perm+st1),         &
     &                     Work(givnum+st1),Work(c+st1),Work(s+st1),    &
     &                     Work(nwork),Iwork(iwk),Info)
               IF ( Info/=0 ) RETURN
               bxst = bx + st1
               CALL DLALSA(icmpq2,Smlsiz,nsize,Nrhs,B(st,1),Ldb,        &
     &                     Work(bxst),N,Work(u+st1),N,Work(vt+st1),     &
     &                     Iwork(k+st1),Work(difl+st1),Work(difr+st1),  &
     &                     Work(z+st1),Work(poles+st1),Iwork(givptr+st1)&
     &                     ,Iwork(givcol+st1),N,Iwork(perm+st1),        &
     &                     Work(givnum+st1),Work(c+st1),Work(s+st1),    &
     &                     Work(nwork),Iwork(iwk),Info)
               IF ( Info/=0 ) RETURN
            ENDIF
            st = i + 1
         ENDIF
      ENDDO
!
!     Apply the singular values and treat the tiny ones as zero.
!
      tol = rcnd*ABS(D(IDAMAX(N,D,1)))
!
      DO i = 1 , N
!
!        Some of the elements in D can be negative because 1-by-1
!        subproblems were not solved explicitly.
!
         IF ( ABS(D(i))<=tol ) THEN
            CALL DLASET('A',1,Nrhs,ZERO,ZERO,Work(bx+i-1),N)
         ELSE
            Rank = Rank + 1
            CALL DLASCL('G',0,0,D(i),ONE,1,Nrhs,Work(bx+i-1),N,Info)
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
            CALL DCOPY(Nrhs,Work(bxst),N,B(st,1),Ldb)
         ELSEIF ( nsize<=Smlsiz ) THEN
            CALL DGEMM('T','N',nsize,Nrhs,nsize,ONE,Work(vt+st1),N,     &
     &                 Work(bxst),N,ZERO,B(st,1),Ldb)
         ELSE
            CALL DLALSA(icmpq2,Smlsiz,nsize,Nrhs,Work(bxst),N,B(st,1),  &
     &                  Ldb,Work(u+st1),N,Work(vt+st1),Iwork(k+st1),    &
     &                  Work(difl+st1),Work(difr+st1),Work(z+st1),      &
     &                  Work(poles+st1),Iwork(givptr+st1),              &
     &                  Iwork(givcol+st1),N,Iwork(perm+st1),            &
     &                  Work(givnum+st1),Work(c+st1),Work(s+st1),       &
     &                  Work(nwork),Iwork(iwk),Info)
            IF ( Info/=0 ) RETURN
         ENDIF
      ENDDO
!
!     Unscale and sort the singular values.
!
      CALL DLASCL('G',0,0,ONE,orgnrm,N,1,D,N,Info)
      CALL DLASRT('D',N,D,Info)
      CALL DLASCL('G',0,0,orgnrm,ONE,N,Nrhs,B,Ldb,Info)
!
!
!     End of DLALSD
!
      END SUBROUTINE DLALSD
