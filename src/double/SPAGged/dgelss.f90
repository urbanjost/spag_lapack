!*==dgelss.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> DGELSS solves overdetermined or underdetermined systems for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGELSS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelss.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelss.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelss.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
!       DOUBLE PRECISION   RCOND
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGELSS computes the minimum norm solution to a real linear least
!> squares problem:
!>
!> Minimize 2-norm(| b - A*x |).
!>
!> using the singular value decomposition (SVD) of A. A is an M-by-N
!> matrix which may be rank-deficient.
!>
!> Several right hand side vectors b and solution vectors x can be
!> handled in a single call; they are stored as the columns of the
!> M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
!> X.
!>
!> The effective rank of A is determined by treating as zero those
!> singular values which are less than RCOND times the largest singular
!> value.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A. N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrices B and X. NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the first min(m,n) rows of A are overwritten with
!>          its right singular vectors, stored rowwise.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!>          On entry, the M-by-NRHS right hand side matrix B.
!>          On exit, B is overwritten by the N-by-NRHS solution
!>          matrix X.  If m >= n and RANK = n, the residual
!>          sum-of-squares for the solution in the i-th column is given
!>          by the sum of squares of elements n+1:m in that column.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,max(M,N)).
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (min(M,N))
!>          The singular values of A in decreasing order.
!>          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
!> \endverbatim
!>
!> \param[in] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          RCOND is used to determine the effective rank of A.
!>          Singular values S(i) <= RCOND*S(1) are treated as zero.
!>          If RCOND < 0, machine precision is used instead.
!> \endverbatim
!>
!> \param[out] RANK
!> \verbatim
!>          RANK is INTEGER
!>          The effective rank of A, i.e., the number of singular values
!>          which are greater than RCOND*S(1).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= 1, and also:
!>          LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )
!>          For good performance, LWORK should generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  the algorithm for computing the SVD failed to converge;
!>                if INFO = i, i off-diagonal elements of an intermediate
!>                bidiagonal form did not converge to zero.
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
!> \ingroup doubleGEsolve
!
!  =====================================================================
      SUBROUTINE DGELSS(M,N,Nrhs,A,Lda,B,Ldb,S,Rcond,Rank,Work,Lwork,   &
     &                  Info)
      USE F77KINDS                        
      USE S_DBDSQR
      USE S_DCOPY
      USE S_DGEBRD
      USE S_DGELQF
      USE S_DGEMM
      USE S_DGEMV
      USE S_DGEQRF
      USE S_DLABAD
      USE S_DLACPY
      USE S_DLAMCH
      USE S_DLANGE
      USE S_DLASCL
      USE S_DLASET
      USE S_DORGBR
      USE S_DORMBR
      USE S_DORMLQ
      USE S_DORMQR
      USE S_DRSCL
      USE S_ILAENV
      USE S_XERBLA
      IMPLICIT NONE
!*--DGELSS197
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anrm , bignum , bnrm , eps , sfmin , smlnum , thr
      INTEGER :: bdspac , bl , chunk , i , iascl , ibscl , ie , il ,    &
     &           itau , itaup , itauq , iwork , ldwork , lwork_dgebrd , &
     &           lwork_dgelqf , lwork_dgeqrf , lwork_dorgbr ,           &
     &           lwork_dormbr , lwork_dormlq , lwork_dormqr , maxmn ,   &
     &           maxwrk , minmn , minwrk , mm , mnthr
      REAL(R8KIND) , DIMENSION(1) :: dum
      LOGICAL :: lquery
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
!     .. Local Arrays ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      minmn = MIN(M,N)
      maxmn = MAX(M,N)
      lquery = (Lwork==-1)
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,maxmn) ) THEN
         Info = -7
      ENDIF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.)
!
      IF ( Info==0 ) THEN
         minwrk = 1
         maxwrk = 1
         IF ( minmn>0 ) THEN
            mm = M
            mnthr = ILAENV(6,'DGELSS',' ',M,N,Nrhs,-1)
            IF ( M>=N .AND. M>=mnthr ) THEN
!
!              Path 1a - overdetermined, with many more rows than
!                        columns
!
!              Compute space needed for DGEQRF
               CALL DGEQRF(M,N,A,Lda,dum(1),dum(1),-1,Info)
               lwork_dgeqrf = dum(1)
!              Compute space needed for DORMQR
               CALL DORMQR('L','T',M,Nrhs,N,A,Lda,dum(1),B,Ldb,dum(1),  &
     &                     -1,Info)
               lwork_dormqr = dum(1)
               mm = N
               maxwrk = MAX(maxwrk,N+lwork_dgeqrf)
               maxwrk = MAX(maxwrk,N+lwork_dormqr)
            ENDIF
            IF ( M>=N ) THEN
!
!              Path 1 - overdetermined or exactly determined
!
!              Compute workspace needed for DBDSQR
!
               bdspac = MAX(1,5*N)
!              Compute space needed for DGEBRD
               CALL DGEBRD(mm,N,A,Lda,S,dum(1),dum(1),dum(1),dum(1),-1, &
     &                     Info)
               lwork_dgebrd = dum(1)
!              Compute space needed for DORMBR
               CALL DORMBR('Q','L','T',mm,Nrhs,N,A,Lda,dum(1),B,Ldb,    &
     &                     dum(1),-1,Info)
               lwork_dormbr = dum(1)
!              Compute space needed for DORGBR
               CALL DORGBR('P',N,N,N,A,Lda,dum(1),dum(1),-1,Info)
               lwork_dorgbr = dum(1)
!              Compute total workspace needed
               maxwrk = MAX(maxwrk,3*N+lwork_dgebrd)
               maxwrk = MAX(maxwrk,3*N+lwork_dormbr)
               maxwrk = MAX(maxwrk,3*N+lwork_dorgbr)
               maxwrk = MAX(maxwrk,bdspac)
               maxwrk = MAX(maxwrk,N*Nrhs)
               minwrk = MAX(3*N+mm,3*N+Nrhs,bdspac)
               maxwrk = MAX(minwrk,maxwrk)
            ENDIF
            IF ( N>M ) THEN
!
!              Compute workspace needed for DBDSQR
!
               bdspac = MAX(1,5*M)
               minwrk = MAX(3*M+Nrhs,3*M+N,bdspac)
               IF ( N>=mnthr ) THEN
!
!                 Path 2a - underdetermined, with many more columns
!                 than rows
!
!                 Compute space needed for DGELQF
                  CALL DGELQF(M,N,A,Lda,dum(1),dum(1),-1,Info)
                  lwork_dgelqf = dum(1)
!                 Compute space needed for DGEBRD
                  CALL DGEBRD(M,M,A,Lda,S,dum(1),dum(1),dum(1),dum(1),  &
     &                        -1,Info)
                  lwork_dgebrd = dum(1)
!                 Compute space needed for DORMBR
                  CALL DORMBR('Q','L','T',M,Nrhs,N,A,Lda,dum(1),B,Ldb,  &
     &                        dum(1),-1,Info)
                  lwork_dormbr = dum(1)
!                 Compute space needed for DORGBR
                  CALL DORGBR('P',M,M,M,A,Lda,dum(1),dum(1),-1,Info)
                  lwork_dorgbr = dum(1)
!                 Compute space needed for DORMLQ
                  CALL DORMLQ('L','T',N,Nrhs,M,A,Lda,dum(1),B,Ldb,dum(1)&
     &                        ,-1,Info)
                  lwork_dormlq = dum(1)
!                 Compute total workspace needed
                  maxwrk = M + lwork_dgelqf
                  maxwrk = MAX(maxwrk,M*M+4*M+lwork_dgebrd)
                  maxwrk = MAX(maxwrk,M*M+4*M+lwork_dormbr)
                  maxwrk = MAX(maxwrk,M*M+4*M+lwork_dorgbr)
                  maxwrk = MAX(maxwrk,M*M+M+bdspac)
                  IF ( Nrhs>1 ) THEN
                     maxwrk = MAX(maxwrk,M*M+M+M*Nrhs)
                  ELSE
                     maxwrk = MAX(maxwrk,M*M+2*M)
                  ENDIF
                  maxwrk = MAX(maxwrk,M+lwork_dormlq)
               ELSE
!
!                 Path 2 - underdetermined
!
!                 Compute space needed for DGEBRD
                  CALL DGEBRD(M,N,A,Lda,S,dum(1),dum(1),dum(1),dum(1),  &
     &                        -1,Info)
                  lwork_dgebrd = dum(1)
!                 Compute space needed for DORMBR
                  CALL DORMBR('Q','L','T',M,Nrhs,M,A,Lda,dum(1),B,Ldb,  &
     &                        dum(1),-1,Info)
                  lwork_dormbr = dum(1)
!                 Compute space needed for DORGBR
                  CALL DORGBR('P',M,N,M,A,Lda,dum(1),dum(1),-1,Info)
                  lwork_dorgbr = dum(1)
                  maxwrk = 3*M + lwork_dgebrd
                  maxwrk = MAX(maxwrk,3*M+lwork_dormbr)
                  maxwrk = MAX(maxwrk,3*M+lwork_dorgbr)
                  maxwrk = MAX(maxwrk,bdspac)
                  maxwrk = MAX(maxwrk,N*Nrhs)
               ENDIF
            ENDIF
            maxwrk = MAX(minwrk,maxwrk)
         ENDIF
         Work(1) = maxwrk
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -12
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGELSS',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) THEN
         Rank = 0
         RETURN
      ENDIF
!
!     Get machine parameters
!
      eps = DLAMCH('P')
      sfmin = DLAMCH('S')
      smlnum = sfmin/eps
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      anrm = DLANGE('M',M,N,A,Lda,Work)
      iascl = 0
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL DLASCL('G',0,0,anrm,smlnum,M,N,A,Lda,Info)
         iascl = 1
      ELSEIF ( anrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL DLASCL('G',0,0,anrm,bignum,M,N,A,Lda,Info)
         iascl = 2
      ELSEIF ( anrm==ZERO ) THEN
!
!        Matrix all zero. Return zero solution.
!
         CALL DLASET('F',MAX(M,N),Nrhs,ZERO,ZERO,B,Ldb)
         CALL DLASET('F',minmn,1,ZERO,ZERO,S,minmn)
         Rank = 0
         GOTO 100
      ENDIF
!
!     Scale B if max element outside range [SMLNUM,BIGNUM]
!
      bnrm = DLANGE('M',M,Nrhs,B,Ldb,Work)
      ibscl = 0
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL DLASCL('G',0,0,bnrm,smlnum,M,Nrhs,B,Ldb,Info)
         ibscl = 1
      ELSEIF ( bnrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL DLASCL('G',0,0,bnrm,bignum,M,Nrhs,B,Ldb,Info)
         ibscl = 2
      ENDIF
!
!     Overdetermined case
!
      IF ( M>=N ) THEN
!
!        Path 1 - overdetermined or exactly determined
!
         mm = M
         IF ( M>=mnthr ) THEN
!
!           Path 1a - overdetermined, with many more rows than columns
!
            mm = N
            itau = 1
            iwork = itau + N
!
!           Compute A=Q*R
!           (Workspace: need 2*N, prefer N+N*NB)
!
            CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),Lwork-iwork+1, &
     &                  Info)
!
!           Multiply B by transpose(Q)
!           (Workspace: need N+NRHS, prefer N+NRHS*NB)
!
            CALL DORMQR('L','T',M,Nrhs,N,A,Lda,Work(itau),B,Ldb,        &
     &                  Work(iwork),Lwork-iwork+1,Info)
!
!           Zero out below R
!
            IF ( N>1 ) CALL DLASET('L',N-1,N-1,ZERO,ZERO,A(2,1),Lda)
         ENDIF
!
         ie = 1
         itauq = ie + N
         itaup = itauq + N
         iwork = itaup + N
!
!        Bidiagonalize R in A
!        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)
!
         CALL DGEBRD(mm,N,A,Lda,S,Work(ie),Work(itauq),Work(itaup),     &
     &               Work(iwork),Lwork-iwork+1,Info)
!
!        Multiply B by transpose of left bidiagonalizing vectors of R
!        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)
!
         CALL DORMBR('Q','L','T',mm,Nrhs,N,A,Lda,Work(itauq),B,Ldb,     &
     &               Work(iwork),Lwork-iwork+1,Info)
!
!        Generate right bidiagonalizing vectors of R in A
!        (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
!
         CALL DORGBR('P',N,N,N,A,Lda,Work(itaup),Work(iwork),           &
     &               Lwork-iwork+1,Info)
         iwork = ie + N
!
!        Perform bidiagonal QR iteration
!          multiply B by transpose of left singular vectors
!          compute right singular vectors in A
!        (Workspace: need BDSPAC)
!
         CALL DBDSQR('U',N,N,0,Nrhs,S,Work(ie),A,Lda,dum,1,B,Ldb,       &
     &               Work(iwork),Info)
         IF ( Info/=0 ) GOTO 100
!
!        Multiply B by reciprocals of singular values
!
         thr = MAX(Rcond*S(1),sfmin)
         IF ( Rcond<ZERO ) thr = MAX(eps*S(1),sfmin)
         Rank = 0
         DO i = 1 , N
            IF ( S(i)>thr ) THEN
               CALL DRSCL(Nrhs,S(i),B(i,1),Ldb)
               Rank = Rank + 1
            ELSE
               CALL DLASET('F',1,Nrhs,ZERO,ZERO,B(i,1),Ldb)
            ENDIF
         ENDDO
!
!        Multiply B by right singular vectors
!        (Workspace: need N, prefer N*NRHS)
!
         IF ( Lwork>=Ldb*Nrhs .AND. Nrhs>1 ) THEN
            CALL DGEMM('T','N',N,Nrhs,N,ONE,A,Lda,B,Ldb,ZERO,Work,Ldb)
            CALL DLACPY('G',N,Nrhs,Work,Ldb,B,Ldb)
         ELSEIF ( Nrhs>1 ) THEN
            chunk = Lwork/N
            DO i = 1 , Nrhs , chunk
               bl = MIN(Nrhs-i+1,chunk)
               CALL DGEMM('T','N',N,bl,N,ONE,A,Lda,B(1,i),Ldb,ZERO,Work,&
     &                    N)
               CALL DLACPY('G',N,bl,Work,N,B(1,i),Ldb)
            ENDDO
         ELSE
            CALL DGEMV('T',N,N,ONE,A,Lda,B,1,ZERO,Work,1)
            CALL DCOPY(N,Work,1,B,1)
         ENDIF
!
      ELSEIF ( N>=mnthr .AND. Lwork>=4*M+M*M+MAX(M,2*M-4,Nrhs,N-3*M) )  &
     &         THEN
!
!        Path 2a - underdetermined, with many more columns than rows
!        and sufficient workspace for an efficient algorithm
!
         ldwork = M
         IF ( Lwork>=MAX(4*M+M*Lda+MAX(M,2*M-4,Nrhs,N-3*M),             &
     &        M*Lda+M+M*Nrhs) ) ldwork = Lda
         itau = 1
         iwork = M + 1
!
!        Compute A=L*Q
!        (Workspace: need 2*M, prefer M+M*NB)
!
         CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),Lwork-iwork+1,    &
     &               Info)
         il = iwork
!
!        Copy L to WORK(IL), zeroing out above it
!
         CALL DLACPY('L',M,M,A,Lda,Work(il),ldwork)
         CALL DLASET('U',M-1,M-1,ZERO,ZERO,Work(il+ldwork),ldwork)
         ie = il + ldwork*M
         itauq = ie + M
         itaup = itauq + M
         iwork = itaup + M
!
!        Bidiagonalize L in WORK(IL)
!        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)
!
         CALL DGEBRD(M,M,Work(il),ldwork,S,Work(ie),Work(itauq),        &
     &               Work(itaup),Work(iwork),Lwork-iwork+1,Info)
!
!        Multiply B by transpose of left bidiagonalizing vectors of L
!        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
!
         CALL DORMBR('Q','L','T',M,Nrhs,M,Work(il),ldwork,Work(itauq),B,&
     &               Ldb,Work(iwork),Lwork-iwork+1,Info)
!
!        Generate right bidiagonalizing vectors of R in WORK(IL)
!        (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB)
!
         CALL DORGBR('P',M,M,M,Work(il),ldwork,Work(itaup),Work(iwork), &
     &               Lwork-iwork+1,Info)
         iwork = ie + M
!
!        Perform bidiagonal QR iteration,
!           computing right singular vectors of L in WORK(IL) and
!           multiplying B by transpose of left singular vectors
!        (Workspace: need M*M+M+BDSPAC)
!
         CALL DBDSQR('U',M,M,0,Nrhs,S,Work(ie),Work(il),ldwork,A,Lda,B, &
     &               Ldb,Work(iwork),Info)
         IF ( Info/=0 ) GOTO 100
!
!        Multiply B by reciprocals of singular values
!
         thr = MAX(Rcond*S(1),sfmin)
         IF ( Rcond<ZERO ) thr = MAX(eps*S(1),sfmin)
         Rank = 0
         DO i = 1 , M
            IF ( S(i)>thr ) THEN
               CALL DRSCL(Nrhs,S(i),B(i,1),Ldb)
               Rank = Rank + 1
            ELSE
               CALL DLASET('F',1,Nrhs,ZERO,ZERO,B(i,1),Ldb)
            ENDIF
         ENDDO
         iwork = ie
!
!        Multiply B by right singular vectors of L in WORK(IL)
!        (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS)
!
         IF ( Lwork>=Ldb*Nrhs+iwork-1 .AND. Nrhs>1 ) THEN
            CALL DGEMM('T','N',M,Nrhs,M,ONE,Work(il),ldwork,B,Ldb,ZERO, &
     &                 Work(iwork),Ldb)
            CALL DLACPY('G',M,Nrhs,Work(iwork),Ldb,B,Ldb)
         ELSEIF ( Nrhs>1 ) THEN
            chunk = (Lwork-iwork+1)/M
            DO i = 1 , Nrhs , chunk
               bl = MIN(Nrhs-i+1,chunk)
               CALL DGEMM('T','N',M,bl,M,ONE,Work(il),ldwork,B(1,i),Ldb,&
     &                    ZERO,Work(iwork),M)
               CALL DLACPY('G',M,bl,Work(iwork),M,B(1,i),Ldb)
            ENDDO
         ELSE
            CALL DGEMV('T',M,M,ONE,Work(il),ldwork,B(1,1),1,ZERO,       &
     &                 Work(iwork),1)
            CALL DCOPY(M,Work(iwork),1,B(1,1),1)
         ENDIF
!
!        Zero out below first M rows of B
!
         CALL DLASET('F',N-M,Nrhs,ZERO,ZERO,B(M+1,1),Ldb)
         iwork = itau + M
!
!        Multiply transpose(Q) by B
!        (Workspace: need M+NRHS, prefer M+NRHS*NB)
!
         CALL DORMLQ('L','T',N,Nrhs,M,A,Lda,Work(itau),B,Ldb,Work(iwork)&
     &               ,Lwork-iwork+1,Info)
!
      ELSE
!
!        Path 2 - remaining underdetermined cases
!
         ie = 1
         itauq = ie + M
         itaup = itauq + M
         iwork = itaup + M
!
!        Bidiagonalize A
!        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
!
         CALL DGEBRD(M,N,A,Lda,S,Work(ie),Work(itauq),Work(itaup),      &
     &               Work(iwork),Lwork-iwork+1,Info)
!
!        Multiply B by transpose of left bidiagonalizing vectors
!        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)
!
         CALL DORMBR('Q','L','T',M,Nrhs,N,A,Lda,Work(itauq),B,Ldb,      &
     &               Work(iwork),Lwork-iwork+1,Info)
!
!        Generate right bidiagonalizing vectors in A
!        (Workspace: need 4*M, prefer 3*M+M*NB)
!
         CALL DORGBR('P',M,N,M,A,Lda,Work(itaup),Work(iwork),           &
     &               Lwork-iwork+1,Info)
         iwork = ie + M
!
!        Perform bidiagonal QR iteration,
!           computing right singular vectors of A in A and
!           multiplying B by transpose of left singular vectors
!        (Workspace: need BDSPAC)
!
         CALL DBDSQR('L',M,N,0,Nrhs,S,Work(ie),A,Lda,dum,1,B,Ldb,       &
     &               Work(iwork),Info)
         IF ( Info/=0 ) GOTO 100
!
!        Multiply B by reciprocals of singular values
!
         thr = MAX(Rcond*S(1),sfmin)
         IF ( Rcond<ZERO ) thr = MAX(eps*S(1),sfmin)
         Rank = 0
         DO i = 1 , M
            IF ( S(i)>thr ) THEN
               CALL DRSCL(Nrhs,S(i),B(i,1),Ldb)
               Rank = Rank + 1
            ELSE
               CALL DLASET('F',1,Nrhs,ZERO,ZERO,B(i,1),Ldb)
            ENDIF
         ENDDO
!
!        Multiply B by right singular vectors of A
!        (Workspace: need N, prefer N*NRHS)
!
         IF ( Lwork>=Ldb*Nrhs .AND. Nrhs>1 ) THEN
            CALL DGEMM('T','N',N,Nrhs,M,ONE,A,Lda,B,Ldb,ZERO,Work,Ldb)
            CALL DLACPY('F',N,Nrhs,Work,Ldb,B,Ldb)
         ELSEIF ( Nrhs>1 ) THEN
            chunk = Lwork/N
            DO i = 1 , Nrhs , chunk
               bl = MIN(Nrhs-i+1,chunk)
               CALL DGEMM('T','N',N,bl,M,ONE,A,Lda,B(1,i),Ldb,ZERO,Work,&
     &                    N)
               CALL DLACPY('F',N,bl,Work,N,B(1,i),Ldb)
            ENDDO
         ELSE
            CALL DGEMV('T',M,N,ONE,A,Lda,B,1,ZERO,Work,1)
            CALL DCOPY(N,Work,1,B,1)
         ENDIF
      ENDIF
!
!     Undo scaling
!
      IF ( iascl==1 ) THEN
         CALL DLASCL('G',0,0,anrm,smlnum,N,Nrhs,B,Ldb,Info)
         CALL DLASCL('G',0,0,smlnum,anrm,minmn,1,S,minmn,Info)
      ELSEIF ( iascl==2 ) THEN
         CALL DLASCL('G',0,0,anrm,bignum,N,Nrhs,B,Ldb,Info)
         CALL DLASCL('G',0,0,bignum,anrm,minmn,1,S,minmn,Info)
      ENDIF
      IF ( ibscl==1 ) THEN
         CALL DLASCL('G',0,0,smlnum,bnrm,N,Nrhs,B,Ldb,Info)
      ELSEIF ( ibscl==2 ) THEN
         CALL DLASCL('G',0,0,bignum,bnrm,N,Nrhs,B,Ldb,Info)
      ENDIF
!
 100  Work(1) = maxwrk
!
!     End of DGELSS
!
      END SUBROUTINE DGELSS
