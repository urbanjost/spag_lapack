!*==cgelsd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> CGELSD computes the minimum-norm solution to a linear least squares problem for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGELSD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgelsd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgelsd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgelsd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
!                          WORK, LWORK, RWORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
!       REAL               RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
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
!> CGELSD computes the minimum-norm solution to a real linear least
!> squares problem:
!>     minimize 2-norm(| b - A*x |)
!> using the singular value decomposition (SVD) of A. A is an M-by-N
!> matrix which may be rank-deficient.
!>
!> Several right hand side vectors b and solution vectors x can be
!> handled in a single call; they are stored as the columns of the
!> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
!> matrix X.
!>
!> The problem is solved in three steps:
!> (1) Reduce the coefficient matrix A to bidiagonal form with
!>     Householder transformations, reducing the original problem
!>     into a "bidiagonal least squares problem" (BLS)
!> (2) Solve the BLS using a divide and conquer approach.
!> (3) Apply back all the Householder transformations to solve
!>     the original least squares problem.
!>
!> The effective rank of A is determined by treating as zero those
!> singular values which are less than RCOND times the largest singular
!> value.
!>
!> The divide and conquer algorithm makes very mild assumptions about
!> floating point arithmetic. It will work on machines with a guard
!> digit in add/subtract, or on those binary machines without guard
!> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!> Cray-2. It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.
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
!>          On exit, A has been destroyed.
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
!>          The dimension of the array WORK. LWORK must be at least 1.
!>          The exact minimum amount of workspace needed depends on M,
!>          N and NRHS. As long as LWORK is at least
!>              2 * N + N * NRHS
!>          if M is greater than or equal to N or
!>              2 * M + M * NRHS
!>          if M is less than N, the code will execute correctly.
!>          For good performance, LWORK should generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the array WORK and the
!>          minimum sizes of the arrays RWORK and IWORK, and returns
!>          these values as the first entries of the WORK, RWORK and
!>          IWORK arrays, and no error message related to LWORK is issued
!>          by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (MAX(1,LRWORK))
!>          LRWORK >=
!>             10*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS +
!>             MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS )
!>          if M is greater than or equal to N or
!>             10*M + 2*M*SMLSIZ + 8*M*NLVL + 3*SMLSIZ*NRHS +
!>             MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS )
!>          if M is less than N, the code will execute correctly.
!>          SMLSIZ is returned by ILAENV and is equal to the maximum
!>          size of the subproblems at the bottom of the computation
!>          tree (usually about 25), and
!>             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
!>          On exit, if INFO = 0, RWORK(1) returns the minimum LRWORK.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          LIWORK >= max(1, 3*MINMN*NLVL + 11*MINMN),
!>          where MINMN = MIN( M,N ).
!>          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup complexGEsolve
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!>       California at Berkeley, USA \n
!>     Osni Marques, LBNL/NERSC, USA \n
!
!  =====================================================================
      SUBROUTINE CGELSD(M,N,Nrhs,A,Lda,B,Ldb,S,Rcond,Rank,Work,Lwork,   &
     &                  Rwork,Iwork,Info)
      USE S_CGEBRD
      USE S_CGELQF
      USE S_CGEQRF
      USE S_CLACPY
      USE S_CLALSD
      USE S_CLANGE
      USE S_CLASCL
      USE S_CLASET
      USE S_CUNMBR
      USE S_CUNMLQ
      USE S_CUNMQR
      USE S_ILAENV
      USE S_SLABAD
      USE S_SLAMCH
      USE S_SLASCL
      USE S_SLASET
      USE S_XERBLA
      IMPLICIT NONE
!*--CGELSD246
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: S
      REAL :: Rcond
      INTEGER :: Rank
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: anrm , bignum , bnrm , eps , sfmin , smlnum
      INTEGER :: iascl , ibscl , ie , il , itau , itaup , itauq ,       &
     &           ldwork , liwork , lrwork , maxmn , maxwrk , minmn ,    &
     &           minwrk , mm , mnthr , nlvl , nrwork , nwork , smlsiz
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
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
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
!     Compute workspace.
!     (Note: Comments in the code beginning "Workspace:" describe the
!     minimal amount of workspace needed at that point in the code,
!     as well as the preferred amount for good performance.
!     NB refers to the optimal block size for the immediately
!     following subroutine, as returned by ILAENV.)
!
      IF ( Info==0 ) THEN
         minwrk = 1
         maxwrk = 1
         liwork = 1
         lrwork = 1
         IF ( minmn>0 ) THEN
            smlsiz = ILAENV(9,'CGELSD',' ',0,0,0,0)
            mnthr = ILAENV(6,'CGELSD',' ',M,N,Nrhs,-1)
            nlvl = MAX(INT(LOG(REAL(minmn)/REAL(smlsiz+1))/LOG(TWO))+1, &
     &             0)
            liwork = 3*minmn*nlvl + 11*minmn
            mm = M
            IF ( M>=N .AND. M>=mnthr ) THEN
!
!              Path 1a - overdetermined, with many more rows than
!                        columns.
!
               mm = N
               maxwrk = MAX(maxwrk,N*ILAENV(1,'CGEQRF',' ',M,N,-1,-1))
               maxwrk = MAX(maxwrk,                                     &
     &                  Nrhs*ILAENV(1,'CUNMQR','LC',M,Nrhs,N,-1))
            ENDIF
            IF ( M>=N ) THEN
!
!              Path 1 - overdetermined or exactly determined.
!
               lrwork = 10*N + 2*N*smlsiz + 8*N*nlvl + 3*smlsiz*Nrhs +  &
     &                  MAX((smlsiz+1)**2,N*(1+Nrhs)+2*Nrhs)
               maxwrk = MAX(maxwrk,2*N+(mm+N)                           &
     &                  *ILAENV(1,'CGEBRD',' ',mm,N,-1,-1))
               maxwrk = MAX(maxwrk,                                     &
     &                  2*N+Nrhs*ILAENV(1,'CUNMBR','QLC',mm,Nrhs,N,-1))
               maxwrk = MAX(maxwrk,2*N+(N-1)                            &
     &                  *ILAENV(1,'CUNMBR','PLN',N,Nrhs,N,-1))
               maxwrk = MAX(maxwrk,2*N+N*Nrhs)
               minwrk = MAX(2*N+mm,2*N+N*Nrhs)
            ENDIF
            IF ( N>M ) THEN
               lrwork = 10*M + 2*M*smlsiz + 8*M*nlvl + 3*smlsiz*Nrhs +  &
     &                  MAX((smlsiz+1)**2,N*(1+Nrhs)+2*Nrhs)
               IF ( N>=mnthr ) THEN
!
!                 Path 2a - underdetermined, with many more columns
!                           than rows.
!
                  maxwrk = M + M*ILAENV(1,'CGELQF',' ',M,N,-1,-1)
                  maxwrk = MAX(maxwrk,                                  &
     &                     M*M+4*M+2*M*ILAENV(1,'CGEBRD',' ',M,M,-1,-1))
                  maxwrk = MAX(maxwrk,                                  &
     &                     M*M+4*M+Nrhs*ILAENV(1,'CUNMBR','QLC',M,Nrhs, &
     &                     M,-1))
                  maxwrk = MAX(maxwrk,M*M+4*M+(M-1)                     &
     &                     *ILAENV(1,'CUNMLQ','LC',N,Nrhs,M,-1))
                  IF ( Nrhs>1 ) THEN
                     maxwrk = MAX(maxwrk,M*M+M+M*Nrhs)
                  ELSE
                     maxwrk = MAX(maxwrk,M*M+2*M)
                  ENDIF
                  maxwrk = MAX(maxwrk,M*M+4*M+M*Nrhs)
!     XXX: Ensure the Path 2a case below is triggered.  The workspace
!     calculation should use queries for all routines eventually.
                  maxwrk = MAX(maxwrk,4*M+M*M+MAX(M,2*M-4,Nrhs,N-3*M))
               ELSE
!
!                 Path 2 - underdetermined.
!
                  maxwrk = 2*M + (N+M)*ILAENV(1,'CGEBRD',' ',M,N,-1,-1)
                  maxwrk = MAX(maxwrk,                                  &
     &                     2*M+Nrhs*ILAENV(1,'CUNMBR','QLC',M,Nrhs,M,-1)&
     &                     )
                  maxwrk = MAX(maxwrk,                                  &
     &                     2*M+M*ILAENV(1,'CUNMBR','PLN',N,Nrhs,M,-1))
                  maxwrk = MAX(maxwrk,2*M+M*Nrhs)
               ENDIF
               minwrk = MAX(2*M+N,2*M+M*Nrhs)
            ENDIF
         ENDIF
         minwrk = MIN(minwrk,maxwrk)
         Work(1) = maxwrk
         Iwork(1) = liwork
         Rwork(1) = lrwork
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -12
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGELSD',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( M==0 .OR. N==0 ) THEN
         Rank = 0
         RETURN
      ENDIF
!
!     Get machine parameters.
!
      eps = SLAMCH('P')
      sfmin = SLAMCH('S')
      smlnum = sfmin/eps
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
!
!     Scale A if max entry outside range [SMLNUM,BIGNUM].
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
!        Scale matrix norm down to BIGNUM.
!
         CALL CLASCL('G',0,0,anrm,bignum,M,N,A,Lda,Info)
         iascl = 2
      ELSEIF ( anrm==ZERO ) THEN
!
!        Matrix all zero. Return zero solution.
!
         CALL CLASET('F',MAX(M,N),Nrhs,CZERO,CZERO,B,Ldb)
         CALL SLASET('F',minmn,1,ZERO,ZERO,S,1)
         Rank = 0
         GOTO 100
      ENDIF
!
!     Scale B if max entry outside range [SMLNUM,BIGNUM].
!
      bnrm = CLANGE('M',M,Nrhs,B,Ldb,Rwork)
      ibscl = 0
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM.
!
         CALL CLASCL('G',0,0,bnrm,smlnum,M,Nrhs,B,Ldb,Info)
         ibscl = 1
      ELSEIF ( bnrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM.
!
         CALL CLASCL('G',0,0,bnrm,bignum,M,Nrhs,B,Ldb,Info)
         ibscl = 2
      ENDIF
!
!     If M < N make sure B(M+1:N,:) = 0
!
      IF ( M<N ) CALL CLASET('F',N-M,Nrhs,CZERO,CZERO,B(M+1,1),Ldb)
!
!     Overdetermined case.
!
      IF ( M>=N ) THEN
!
!        Path 1 - overdetermined or exactly determined.
!
         mm = M
         IF ( M>=mnthr ) THEN
!
!           Path 1a - overdetermined, with many more rows than columns
!
            mm = N
            itau = 1
            nwork = itau + N
!
!           Compute A=Q*R.
!           (RWorkspace: need N)
!           (CWorkspace: need N, prefer N*NB)
!
            CALL CGEQRF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  Info)
!
!           Multiply B by transpose(Q).
!           (RWorkspace: need N)
!           (CWorkspace: need NRHS, prefer NRHS*NB)
!
            CALL CUNMQR('L','C',M,Nrhs,N,A,Lda,Work(itau),B,Ldb,        &
     &                  Work(nwork),Lwork-nwork+1,Info)
!
!           Zero out below R.
!
            IF ( N>1 ) CALL CLASET('L',N-1,N-1,CZERO,CZERO,A(2,1),Lda)
         ENDIF
!
         itauq = 1
         itaup = itauq + N
         nwork = itaup + N
         ie = 1
         nrwork = ie + N
!
!        Bidiagonalize R in A.
!        (RWorkspace: need N)
!        (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)
!
         CALL CGEBRD(mm,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),    &
     &               Work(nwork),Lwork-nwork+1,Info)
!
!        Multiply B by transpose of left bidiagonalizing vectors of R.
!        (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)
!
         CALL CUNMBR('Q','L','C',mm,Nrhs,N,A,Lda,Work(itauq),B,Ldb,     &
     &               Work(nwork),Lwork-nwork+1,Info)
!
!        Solve the bidiagonal least squares problem.
!
         CALL CLALSD('U',smlsiz,N,Nrhs,S,Rwork(ie),B,Ldb,Rcond,Rank,    &
     &               Work(nwork),Rwork(nrwork),Iwork,Info)
         IF ( Info/=0 ) GOTO 100
!
!        Multiply B by right bidiagonalizing vectors of R.
!
         CALL CUNMBR('P','L','N',N,Nrhs,N,A,Lda,Work(itaup),B,Ldb,      &
     &               Work(nwork),Lwork-nwork+1,Info)
!
      ELSEIF ( N>=mnthr .AND. Lwork>=4*M+M*M+MAX(M,2*M-4,Nrhs,N-3*M) )  &
     &         THEN
!
!        Path 2a - underdetermined, with many more columns than rows
!        and sufficient workspace for an efficient algorithm.
!
         ldwork = M
         IF ( Lwork>=MAX(4*M+M*Lda+MAX(M,2*M-4,Nrhs,N-3*M),             &
     &        M*Lda+M+M*Nrhs) ) ldwork = Lda
         itau = 1
         nwork = M + 1
!
!        Compute A=L*Q.
!        (CWorkspace: need 2*M, prefer M+M*NB)
!
         CALL CGELQF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1,    &
     &               Info)
         il = nwork
!
!        Copy L to WORK(IL), zeroing out above its diagonal.
!
         CALL CLACPY('L',M,M,A,Lda,Work(il),ldwork)
         CALL CLASET('U',M-1,M-1,CZERO,CZERO,Work(il+ldwork),ldwork)
         itauq = il + ldwork*M
         itaup = itauq + M
         nwork = itaup + M
         ie = 1
         nrwork = ie + M
!
!        Bidiagonalize L in WORK(IL).
!        (RWorkspace: need M)
!        (CWorkspace: need M*M+4*M, prefer M*M+4*M+2*M*NB)
!
         CALL CGEBRD(M,M,Work(il),ldwork,S,Rwork(ie),Work(itauq),       &
     &               Work(itaup),Work(nwork),Lwork-nwork+1,Info)
!
!        Multiply B by transpose of left bidiagonalizing vectors of L.
!        (CWorkspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
!
         CALL CUNMBR('Q','L','C',M,Nrhs,M,Work(il),ldwork,Work(itauq),B,&
     &               Ldb,Work(nwork),Lwork-nwork+1,Info)
!
!        Solve the bidiagonal least squares problem.
!
         CALL CLALSD('U',smlsiz,M,Nrhs,S,Rwork(ie),B,Ldb,Rcond,Rank,    &
     &               Work(nwork),Rwork(nrwork),Iwork,Info)
         IF ( Info/=0 ) GOTO 100
!
!        Multiply B by right bidiagonalizing vectors of L.
!
         CALL CUNMBR('P','L','N',M,Nrhs,M,Work(il),ldwork,Work(itaup),B,&
     &               Ldb,Work(nwork),Lwork-nwork+1,Info)
!
!        Zero out below first M rows of B.
!
         CALL CLASET('F',N-M,Nrhs,CZERO,CZERO,B(M+1,1),Ldb)
         nwork = itau + M
!
!        Multiply transpose(Q) by B.
!        (CWorkspace: need NRHS, prefer NRHS*NB)
!
         CALL CUNMLQ('L','C',N,Nrhs,M,A,Lda,Work(itau),B,Ldb,Work(nwork)&
     &               ,Lwork-nwork+1,Info)
!
      ELSE
!
!        Path 2 - remaining underdetermined cases.
!
         itauq = 1
         itaup = itauq + M
         nwork = itaup + M
         ie = 1
         nrwork = ie + M
!
!        Bidiagonalize A.
!        (RWorkspace: need M)
!        (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
!
         CALL CGEBRD(M,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),     &
     &               Work(nwork),Lwork-nwork+1,Info)
!
!        Multiply B by transpose of left bidiagonalizing vectors.
!        (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)
!
         CALL CUNMBR('Q','L','C',M,Nrhs,N,A,Lda,Work(itauq),B,Ldb,      &
     &               Work(nwork),Lwork-nwork+1,Info)
!
!        Solve the bidiagonal least squares problem.
!
         CALL CLALSD('L',smlsiz,M,Nrhs,S,Rwork(ie),B,Ldb,Rcond,Rank,    &
     &               Work(nwork),Rwork(nrwork),Iwork,Info)
         IF ( Info/=0 ) GOTO 100
!
!        Multiply B by right bidiagonalizing vectors of A.
!
         CALL CUNMBR('P','L','N',N,Nrhs,M,A,Lda,Work(itaup),B,Ldb,      &
     &               Work(nwork),Lwork-nwork+1,Info)
!
      ENDIF
!
!     Undo scaling.
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
!
 100  Work(1) = maxwrk
      Iwork(1) = liwork
      Rwork(1) = lrwork
!
!     End of CGELSD
!
      END SUBROUTINE CGELSD
