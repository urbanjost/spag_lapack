!*==cgelss.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> CGELSS solves overdetermined or underdetermined systems for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGELSS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgelss.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgelss.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgelss.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
!                          WORK, LWORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
!       REAL               RCOND
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * ), S( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGELSS computes the minimum norm solution to a complex linear
!> least squares problem:
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the first min(m,n) rows of A are overwritten with
!>          its right singular vectors, stored rowwise.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          On entry, the M-by-NRHS right hand side matrix B.
!>          On exit, B is overwritten by the N-by-NRHS solution matrix X.
!>          If m >= n and RANK = n, the residual sum-of-squares for
!>          the solution in the i-th column is given by the sum of
!>          squares of the modulus of elements n+1:m in that column.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M,N).
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (min(M,N))
!>          The singular values of A in decreasing order.
!>          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
!> \endverbatim
!>
!> \param[in] RCOND
!> \verbatim
!>          RCOND is REAL
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
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= 1, and also:
!>          LWORK >=  2*min(M,N) + max(M,N,NRHS)
!>          For good performance, LWORK should generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (5*min(M,N))
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
!> \date June 2016
!
!> \ingroup complexGEsolve
!
!  =====================================================================
      SUBROUTINE CGELSS(M,N,Nrhs,A,Lda,B,Ldb,S,Rcond,Rank,Work,Lwork,   &
     &                  Rwork,Info)
      IMPLICIT NONE
!*--CGELSS182
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldb , Lwork , M , N , Nrhs , Rank
      REAL Rcond
!     ..
!     .. Array Arguments ..
      REAL Rwork(*) , S(*)
      COMPLEX A(Lda,*) , B(Ldb,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL lquery
      INTEGER bl , chunk , i , iascl , ibscl , ie , il , irwork , itau ,&
     &        itaup , itauq , iwork , ldwork , maxmn , maxwrk , minmn , &
     &        minwrk , mm , mnthr
      INTEGER lwork_cgeqrf , lwork_cunmqr , lwork_cgebrd ,              &
     &        lwork_cunmbr , lwork_cungbr , lwork_cunmlq , lwork_cgelqf
      REAL anrm , bignum , bnrm , eps , sfmin , smlnum , thr
!     ..
!     .. Local Arrays ..
      COMPLEX dum(1)
!     ..
!     .. External Subroutines ..
      EXTERNAL CBDSQR , CCOPY , CGEBRD , CGELQF , CGEMM , CGEMV ,       &
     &         CGEQRF , CLACPY , CLASCL , CLASET , CSRSCL , CUNGBR ,    &
     &         CUNMBR , CUNMLQ , CUNMQR , SLABAD , SLASCL , SLASET ,    &
     &         XERBLA
!     ..
!     .. External Functions ..
      INTEGER ILAENV
      REAL CLANGE , SLAMCH
      EXTERNAL ILAENV , CLANGE , SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
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
!       CWorkspace refers to complex workspace, and RWorkspace refers
!       to real workspace. NB refers to the optimal block size for the
!       immediately following subroutine, as returned by ILAENV.)
!
      IF ( Info==0 ) THEN
         minwrk = 1
         maxwrk = 1
         IF ( minmn>0 ) THEN
            mm = M
            mnthr = ILAENV(6,'CGELSS',' ',M,N,Nrhs,-1)
            IF ( M>=N .AND. M>=mnthr ) THEN
!
!              Path 1a - overdetermined, with many more rows than
!                        columns
!
!              Compute space needed for CGEQRF
               CALL CGEQRF(M,N,A,Lda,dum(1),dum(1),-1,Info)
               lwork_cgeqrf = dum(1)
!              Compute space needed for CUNMQR
               CALL CUNMQR('L','C',M,Nrhs,N,A,Lda,dum(1),B,Ldb,dum(1),  &
     &                     -1,Info)
               lwork_cunmqr = dum(1)
               mm = N
               maxwrk = MAX(maxwrk,N+N*ILAENV(1,'CGEQRF',' ',M,N,-1,-1))
               maxwrk = MAX(maxwrk,                                     &
     &                  N+Nrhs*ILAENV(1,'CUNMQR','LC',M,Nrhs,N,-1))
            ENDIF
            IF ( M>=N ) THEN
!
!              Path 1 - overdetermined or exactly determined
!
!              Compute space needed for CGEBRD
               CALL CGEBRD(mm,N,A,Lda,S,S,dum(1),dum(1),dum(1),-1,Info)
               lwork_cgebrd = dum(1)
!              Compute space needed for CUNMBR
               CALL CUNMBR('Q','L','C',mm,Nrhs,N,A,Lda,dum(1),B,Ldb,    &
     &                     dum(1),-1,Info)
               lwork_cunmbr = dum(1)
!              Compute space needed for CUNGBR
               CALL CUNGBR('P',N,N,N,A,Lda,dum(1),dum(1),-1,Info)
               lwork_cungbr = dum(1)
!              Compute total workspace needed
               maxwrk = MAX(maxwrk,2*N+lwork_cgebrd)
               maxwrk = MAX(maxwrk,2*N+lwork_cunmbr)
               maxwrk = MAX(maxwrk,2*N+lwork_cungbr)
               maxwrk = MAX(maxwrk,N*Nrhs)
               minwrk = 2*N + MAX(Nrhs,M)
            ENDIF
            IF ( N>M ) THEN
               minwrk = 2*M + MAX(Nrhs,N)
               IF ( N>=mnthr ) THEN
!
!                 Path 2a - underdetermined, with many more columns
!                 than rows
!
!                 Compute space needed for CGELQF
                  CALL CGELQF(M,N,A,Lda,dum(1),dum(1),-1,Info)
                  lwork_cgelqf = dum(1)
!                 Compute space needed for CGEBRD
                  CALL CGEBRD(M,M,A,Lda,S,S,dum(1),dum(1),dum(1),-1,    &
     &                        Info)
                  lwork_cgebrd = dum(1)
!                 Compute space needed for CUNMBR
                  CALL CUNMBR('Q','L','C',M,Nrhs,N,A,Lda,dum(1),B,Ldb,  &
     &                        dum(1),-1,Info)
                  lwork_cunmbr = dum(1)
!                 Compute space needed for CUNGBR
                  CALL CUNGBR('P',M,M,M,A,Lda,dum(1),dum(1),-1,Info)
                  lwork_cungbr = dum(1)
!                 Compute space needed for CUNMLQ
                  CALL CUNMLQ('L','C',N,Nrhs,M,A,Lda,dum(1),B,Ldb,dum(1)&
     &                        ,-1,Info)
                  lwork_cunmlq = dum(1)
!                 Compute total workspace needed
                  maxwrk = M + lwork_cgelqf
                  maxwrk = MAX(maxwrk,3*M+M*M+lwork_cgebrd)
                  maxwrk = MAX(maxwrk,3*M+M*M+lwork_cunmbr)
                  maxwrk = MAX(maxwrk,3*M+M*M+lwork_cungbr)
                  IF ( Nrhs>1 ) THEN
                     maxwrk = MAX(maxwrk,M*M+M+M*Nrhs)
                  ELSE
                     maxwrk = MAX(maxwrk,M*M+2*M)
                  ENDIF
                  maxwrk = MAX(maxwrk,M+lwork_cunmlq)
               ELSE
!
!                 Path 2 - underdetermined
!
!                 Compute space needed for CGEBRD
                  CALL CGEBRD(M,N,A,Lda,S,S,dum(1),dum(1),dum(1),-1,    &
     &                        Info)
                  lwork_cgebrd = dum(1)
!                 Compute space needed for CUNMBR
                  CALL CUNMBR('Q','L','C',M,Nrhs,M,A,Lda,dum(1),B,Ldb,  &
     &                        dum(1),-1,Info)
                  lwork_cunmbr = dum(1)
!                 Compute space needed for CUNGBR
                  CALL CUNGBR('P',M,N,M,A,Lda,dum(1),dum(1),-1,Info)
                  lwork_cungbr = dum(1)
                  maxwrk = 2*M + lwork_cgebrd
                  maxwrk = MAX(maxwrk,2*M+lwork_cunmbr)
                  maxwrk = MAX(maxwrk,2*M+lwork_cungbr)
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
         CALL XERBLA('CGELSS',-Info)
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
      eps = SLAMCH('P')
      sfmin = SLAMCH('S')
      smlnum = sfmin/eps
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      anrm = CLANGE('M',M,N,A,Lda,Rwork)
      iascl = 0
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL CLASCL('G',0,0,anrm,smlnum,M,N,A,Lda,Info)
         iascl = 1
      ELSEIF ( anrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL CLASCL('G',0,0,anrm,bignum,M,N,A,Lda,Info)
         iascl = 2
      ELSEIF ( anrm==ZERO ) THEN
!
!        Matrix all zero. Return zero solution.
!
         CALL CLASET('F',MAX(M,N),Nrhs,CZERO,CZERO,B,Ldb)
         CALL SLASET('F',minmn,1,ZERO,ZERO,S,minmn)
         Rank = 0
         GOTO 100
      ENDIF
!
!     Scale B if max element outside range [SMLNUM,BIGNUM]
!
      bnrm = CLANGE('M',M,Nrhs,B,Ldb,Rwork)
      ibscl = 0
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL CLASCL('G',0,0,bnrm,smlnum,M,Nrhs,B,Ldb,Info)
         ibscl = 1
      ELSEIF ( bnrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL CLASCL('G',0,0,bnrm,bignum,M,Nrhs,B,Ldb,Info)
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
!           (CWorkspace: need 2*N, prefer N+N*NB)
!           (RWorkspace: none)
!
            CALL CGEQRF(M,N,A,Lda,Work(itau),Work(iwork),Lwork-iwork+1, &
     &                  Info)
!
!           Multiply B by transpose(Q)
!           (CWorkspace: need N+NRHS, prefer N+NRHS*NB)
!           (RWorkspace: none)
!
            CALL CUNMQR('L','C',M,Nrhs,N,A,Lda,Work(itau),B,Ldb,        &
     &                  Work(iwork),Lwork-iwork+1,Info)
!
!           Zero out below R
!
            IF ( N>1 ) CALL CLASET('L',N-1,N-1,CZERO,CZERO,A(2,1),Lda)
         ENDIF
!
         ie = 1
         itauq = 1
         itaup = itauq + N
         iwork = itaup + N
!
!        Bidiagonalize R in A
!        (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)
!        (RWorkspace: need N)
!
         CALL CGEBRD(mm,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),    &
     &               Work(iwork),Lwork-iwork+1,Info)
!
!        Multiply B by transpose of left bidiagonalizing vectors of R
!        (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)
!        (RWorkspace: none)
!
         CALL CUNMBR('Q','L','C',mm,Nrhs,N,A,Lda,Work(itauq),B,Ldb,     &
     &               Work(iwork),Lwork-iwork+1,Info)
!
!        Generate right bidiagonalizing vectors of R in A
!        (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!        (RWorkspace: none)
!
         CALL CUNGBR('P',N,N,N,A,Lda,Work(itaup),Work(iwork),           &
     &               Lwork-iwork+1,Info)
         irwork = ie + N
!
!        Perform bidiagonal QR iteration
!          multiply B by transpose of left singular vectors
!          compute right singular vectors in A
!        (CWorkspace: none)
!        (RWorkspace: need BDSPAC)
!
         CALL CBDSQR('U',N,N,0,Nrhs,S,Rwork(ie),A,Lda,dum,1,B,Ldb,      &
     &               Rwork(irwork),Info)
         IF ( Info/=0 ) GOTO 100
!
!        Multiply B by reciprocals of singular values
!
         thr = MAX(Rcond*S(1),sfmin)
         IF ( Rcond<ZERO ) thr = MAX(eps*S(1),sfmin)
         Rank = 0
         DO i = 1 , N
            IF ( S(i)>thr ) THEN
               CALL CSRSCL(Nrhs,S(i),B(i,1),Ldb)
               Rank = Rank + 1
            ELSE
               CALL CLASET('F',1,Nrhs,CZERO,CZERO,B(i,1),Ldb)
            ENDIF
         ENDDO
!
!        Multiply B by right singular vectors
!        (CWorkspace: need N, prefer N*NRHS)
!        (RWorkspace: none)
!
         IF ( Lwork>=Ldb*Nrhs .AND. Nrhs>1 ) THEN
            CALL CGEMM('C','N',N,Nrhs,N,CONE,A,Lda,B,Ldb,CZERO,Work,Ldb)
            CALL CLACPY('G',N,Nrhs,Work,Ldb,B,Ldb)
         ELSEIF ( Nrhs>1 ) THEN
            chunk = Lwork/N
            DO i = 1 , Nrhs , chunk
               bl = MIN(Nrhs-i+1,chunk)
               CALL CGEMM('C','N',N,bl,N,CONE,A,Lda,B(1,i),Ldb,CZERO,   &
     &                    Work,N)
               CALL CLACPY('G',N,bl,Work,N,B(1,i),Ldb)
            ENDDO
         ELSE
            CALL CGEMV('C',N,N,CONE,A,Lda,B,1,CZERO,Work,1)
            CALL CCOPY(N,Work,1,B,1)
         ENDIF
!
      ELSEIF ( N>=mnthr .AND. Lwork>=3*M+M*M+MAX(M,Nrhs,N-2*M) ) THEN
!
!        Underdetermined case, M much less than N
!
!        Path 2a - underdetermined, with many more columns than rows
!        and sufficient workspace for an efficient algorithm
!
         ldwork = M
         IF ( Lwork>=3*M+M*Lda+MAX(M,Nrhs,N-2*M) ) ldwork = Lda
         itau = 1
         iwork = M + 1
!
!        Compute A=L*Q
!        (CWorkspace: need 2*M, prefer M+M*NB)
!        (RWorkspace: none)
!
         CALL CGELQF(M,N,A,Lda,Work(itau),Work(iwork),Lwork-iwork+1,    &
     &               Info)
         il = iwork
!
!        Copy L to WORK(IL), zeroing out above it
!
         CALL CLACPY('L',M,M,A,Lda,Work(il),ldwork)
         CALL CLASET('U',M-1,M-1,CZERO,CZERO,Work(il+ldwork),ldwork)
         ie = 1
         itauq = il + ldwork*M
         itaup = itauq + M
         iwork = itaup + M
!
!        Bidiagonalize L in WORK(IL)
!        (CWorkspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
!        (RWorkspace: need M)
!
         CALL CGEBRD(M,M,Work(il),ldwork,S,Rwork(ie),Work(itauq),       &
     &               Work(itaup),Work(iwork),Lwork-iwork+1,Info)
!
!        Multiply B by transpose of left bidiagonalizing vectors of L
!        (CWorkspace: need M*M+3*M+NRHS, prefer M*M+3*M+NRHS*NB)
!        (RWorkspace: none)
!
         CALL CUNMBR('Q','L','C',M,Nrhs,M,Work(il),ldwork,Work(itauq),B,&
     &               Ldb,Work(iwork),Lwork-iwork+1,Info)
!
!        Generate right bidiagonalizing vectors of R in WORK(IL)
!        (CWorkspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)
!        (RWorkspace: none)
!
         CALL CUNGBR('P',M,M,M,Work(il),ldwork,Work(itaup),Work(iwork), &
     &               Lwork-iwork+1,Info)
         irwork = ie + M
!
!        Perform bidiagonal QR iteration, computing right singular
!        vectors of L in WORK(IL) and multiplying B by transpose of
!        left singular vectors
!        (CWorkspace: need M*M)
!        (RWorkspace: need BDSPAC)
!
         CALL CBDSQR('U',M,M,0,Nrhs,S,Rwork(ie),Work(il),ldwork,A,Lda,B,&
     &               Ldb,Rwork(irwork),Info)
         IF ( Info/=0 ) GOTO 100
!
!        Multiply B by reciprocals of singular values
!
         thr = MAX(Rcond*S(1),sfmin)
         IF ( Rcond<ZERO ) thr = MAX(eps*S(1),sfmin)
         Rank = 0
         DO i = 1 , M
            IF ( S(i)>thr ) THEN
               CALL CSRSCL(Nrhs,S(i),B(i,1),Ldb)
               Rank = Rank + 1
            ELSE
               CALL CLASET('F',1,Nrhs,CZERO,CZERO,B(i,1),Ldb)
            ENDIF
         ENDDO
         iwork = il + M*ldwork
!
!        Multiply B by right singular vectors of L in WORK(IL)
!        (CWorkspace: need M*M+2*M, prefer M*M+M+M*NRHS)
!        (RWorkspace: none)
!
         IF ( Lwork>=Ldb*Nrhs+iwork-1 .AND. Nrhs>1 ) THEN
            CALL CGEMM('C','N',M,Nrhs,M,CONE,Work(il),ldwork,B,Ldb,     &
     &                 CZERO,Work(iwork),Ldb)
            CALL CLACPY('G',M,Nrhs,Work(iwork),Ldb,B,Ldb)
         ELSEIF ( Nrhs>1 ) THEN
            chunk = (Lwork-iwork+1)/M
            DO i = 1 , Nrhs , chunk
               bl = MIN(Nrhs-i+1,chunk)
               CALL CGEMM('C','N',M,bl,M,CONE,Work(il),ldwork,B(1,i),   &
     &                    Ldb,CZERO,Work(iwork),M)
               CALL CLACPY('G',M,bl,Work(iwork),M,B(1,i),Ldb)
            ENDDO
         ELSE
            CALL CGEMV('C',M,M,CONE,Work(il),ldwork,B(1,1),1,CZERO,     &
     &                 Work(iwork),1)
            CALL CCOPY(M,Work(iwork),1,B(1,1),1)
         ENDIF
!
!        Zero out below first M rows of B
!
         CALL CLASET('F',N-M,Nrhs,CZERO,CZERO,B(M+1,1),Ldb)
         iwork = itau + M
!
!        Multiply transpose(Q) by B
!        (CWorkspace: need M+NRHS, prefer M+NHRS*NB)
!        (RWorkspace: none)
!
         CALL CUNMLQ('L','C',N,Nrhs,M,A,Lda,Work(itau),B,Ldb,Work(iwork)&
     &               ,Lwork-iwork+1,Info)
!
      ELSE
!
!        Path 2 - remaining underdetermined cases
!
         ie = 1
         itauq = 1
         itaup = itauq + M
         iwork = itaup + M
!
!        Bidiagonalize A
!        (CWorkspace: need 3*M, prefer 2*M+(M+N)*NB)
!        (RWorkspace: need N)
!
         CALL CGEBRD(M,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),     &
     &               Work(iwork),Lwork-iwork+1,Info)
!
!        Multiply B by transpose of left bidiagonalizing vectors
!        (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)
!        (RWorkspace: none)
!
         CALL CUNMBR('Q','L','C',M,Nrhs,N,A,Lda,Work(itauq),B,Ldb,      &
     &               Work(iwork),Lwork-iwork+1,Info)
!
!        Generate right bidiagonalizing vectors in A
!        (CWorkspace: need 3*M, prefer 2*M+M*NB)
!        (RWorkspace: none)
!
         CALL CUNGBR('P',M,N,M,A,Lda,Work(itaup),Work(iwork),           &
     &               Lwork-iwork+1,Info)
         irwork = ie + M
!
!        Perform bidiagonal QR iteration,
!           computing right singular vectors of A in A and
!           multiplying B by transpose of left singular vectors
!        (CWorkspace: none)
!        (RWorkspace: need BDSPAC)
!
         CALL CBDSQR('L',M,N,0,Nrhs,S,Rwork(ie),A,Lda,dum,1,B,Ldb,      &
     &               Rwork(irwork),Info)
         IF ( Info/=0 ) GOTO 100
!
!        Multiply B by reciprocals of singular values
!
         thr = MAX(Rcond*S(1),sfmin)
         IF ( Rcond<ZERO ) thr = MAX(eps*S(1),sfmin)
         Rank = 0
         DO i = 1 , M
            IF ( S(i)>thr ) THEN
               CALL CSRSCL(Nrhs,S(i),B(i,1),Ldb)
               Rank = Rank + 1
            ELSE
               CALL CLASET('F',1,Nrhs,CZERO,CZERO,B(i,1),Ldb)
            ENDIF
         ENDDO
!
!        Multiply B by right singular vectors of A
!        (CWorkspace: need N, prefer N*NRHS)
!        (RWorkspace: none)
!
         IF ( Lwork>=Ldb*Nrhs .AND. Nrhs>1 ) THEN
            CALL CGEMM('C','N',N,Nrhs,M,CONE,A,Lda,B,Ldb,CZERO,Work,Ldb)
            CALL CLACPY('G',N,Nrhs,Work,Ldb,B,Ldb)
         ELSEIF ( Nrhs>1 ) THEN
            chunk = Lwork/N
            DO i = 1 , Nrhs , chunk
               bl = MIN(Nrhs-i+1,chunk)
               CALL CGEMM('C','N',N,bl,M,CONE,A,Lda,B(1,i),Ldb,CZERO,   &
     &                    Work,N)
               CALL CLACPY('F',N,bl,Work,N,B(1,i),Ldb)
            ENDDO
         ELSE
            CALL CGEMV('C',M,N,CONE,A,Lda,B,1,CZERO,Work,1)
            CALL CCOPY(N,Work,1,B,1)
         ENDIF
      ENDIF
!
!     Undo scaling
!
      IF ( iascl==1 ) THEN
         CALL CLASCL('G',0,0,anrm,smlnum,N,Nrhs,B,Ldb,Info)
         CALL SLASCL('G',0,0,smlnum,anrm,minmn,1,S,minmn,Info)
      ELSEIF ( iascl==2 ) THEN
         CALL CLASCL('G',0,0,anrm,bignum,N,Nrhs,B,Ldb,Info)
         CALL SLASCL('G',0,0,bignum,anrm,minmn,1,S,minmn,Info)
      ENDIF
      IF ( ibscl==1 ) THEN
         CALL CLASCL('G',0,0,smlnum,bnrm,N,Nrhs,B,Ldb,Info)
      ELSEIF ( ibscl==2 ) THEN
         CALL CLASCL('G',0,0,bignum,bnrm,N,Nrhs,B,Ldb,Info)
      ENDIF
 100  Work(1) = maxwrk
!
!     End of CGELSS
!
      END SUBROUTINE CGELSS
