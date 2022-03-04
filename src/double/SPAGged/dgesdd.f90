!*==dgesdd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DGESDD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGESDD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesdd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesdd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesdd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT,
!                          WORK, LWORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ
!       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
!      $                   VT( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGESDD computes the singular value decomposition (SVD) of a real
!> M-by-N matrix A, optionally computing the left and right singular
!> vectors.  If singular vectors are desired, it uses a
!> divide-and-conquer algorithm.
!>
!> The SVD is written
!>
!>      A = U * SIGMA * transpose(V)
!>
!> where SIGMA is an M-by-N matrix which is zero except for its
!> min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
!> V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
!> are the singular values of A; they are real and non-negative, and
!> are returned in descending order.  The first min(m,n) columns of
!> U and V are the left and right singular vectors of A.
!>
!> Note that the routine returns VT = V**T, not V.
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
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER*1
!>          Specifies options for computing all or part of the matrix U:
!>          = 'A':  all M columns of U and all N rows of V**T are
!>                  returned in the arrays U and VT;
!>          = 'S':  the first min(M,N) columns of U and the first
!>                  min(M,N) rows of V**T are returned in the arrays U
!>                  and VT;
!>          = 'O':  If M >= N, the first N columns of U are overwritten
!>                  on the array A and all rows of V**T are returned in
!>                  the array VT;
!>                  otherwise, all columns of U are returned in the
!>                  array U and the first M rows of V**T are overwritten
!>                  in the array A;
!>          = 'N':  no columns of U or rows of V**T are computed.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the input matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the input matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit,
!>          if JOBZ = 'O',  A is overwritten with the first N columns
!>                          of U (the left singular vectors, stored
!>                          columnwise) if M >= N;
!>                          A is overwritten with the first M rows
!>                          of V**T (the right singular vectors, stored
!>                          rowwise) otherwise.
!>          if JOBZ .ne. 'O', the contents of A are destroyed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (min(M,N))
!>          The singular values of A, sorted so that S(i) >= S(i+1).
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU,UCOL)
!>          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
!>          UCOL = min(M,N) if JOBZ = 'S'.
!>          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
!>          orthogonal matrix U;
!>          if JOBZ = 'S', U contains the first min(M,N) columns of U
!>          (the left singular vectors, stored columnwise);
!>          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= 1; if
!>          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is DOUBLE PRECISION array, dimension (LDVT,N)
!>          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
!>          N-by-N orthogonal matrix V**T;
!>          if JOBZ = 'S', VT contains the first min(M,N) rows of
!>          V**T (the right singular vectors, stored rowwise);
!>          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.  LDVT >= 1;
!>          if JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
!>          if JOBZ = 'S', LDVT >= min(M,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= 1.
!>          If LWORK = -1, a workspace query is assumed.  The optimal
!>          size for the WORK array is calculated and stored in WORK(1),
!>          and no other work except argument checking is performed.
!>
!>          Let mx = max(M,N) and mn = min(M,N).
!>          If JOBZ = 'N', LWORK >= 3*mn + max( mx, 7*mn ).
!>          If JOBZ = 'O', LWORK >= 3*mn + max( mx, 5*mn*mn + 4*mn ).
!>          If JOBZ = 'S', LWORK >= 4*mn*mn + 7*mn.
!>          If JOBZ = 'A', LWORK >= 4*mn*mn + 6*mn + mx.
!>          These are not tight minimums in all cases; see comments inside code.
!>          For good performance, LWORK should generally be larger;
!>          a query is recommended.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (8*min(M,N))
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  DBDSDC did not converge, updating process failed.
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
!> \ingroup doubleGEsing
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Huan Ren, Computer Science Division, University of
!>     California at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE DGESDD(Jobz,M,N,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,Iwork,&
     &                  Info)
      USE F77KINDS                        
      USE S_DBDSDC
      USE S_DGEBRD
      USE S_DGELQF
      USE S_DGEMM
      USE S_DGEQRF
      USE S_DISNAN
      USE S_DLACPY
      USE S_DLAMCH
      USE S_DLANGE
      USE S_DLASCL
      USE S_DLASET
      USE S_DORGBR
      USE S_DORGLQ
      USE S_DORGQR
      USE S_DORMBR
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DGESDD240
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobz
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anrm , bignum , eps , smlnum
      INTEGER :: bdspac , blk , chunk , i , ie , ierr , il , ir , iscl ,&
     &           itau , itaup , itauq , iu , ivt , ldwkvt , ldwrkl ,    &
     &           ldwrkr , ldwrku , lwork_dgebrd_mm , lwork_dgebrd_mn ,  &
     &           lwork_dgebrd_nn , lwork_dgelqf_mn , lwork_dgeqrf_mn ,  &
     &           lwork_dorgbr_p_mm , lwork_dorgbr_q_nn ,                &
     &           lwork_dorglq_mn , lwork_dorglq_nn , lwork_dorgqr_mm ,  &
     &           lwork_dorgqr_mn , lwork_dormbr_prt_mm ,                &
     &           lwork_dormbr_prt_mn , lwork_dormbr_prt_nn ,            &
     &           lwork_dormbr_qln_mm , lwork_dormbr_qln_mn ,            &
     &           lwork_dormbr_qln_nn , maxwrk , minmn , minwrk , mnthr ,&
     &           nwork , wrkbl
      REAL(R8KIND) , DIMENSION(1) :: dum
      INTEGER , DIMENSION(1) :: idum
      LOGICAL :: lquery , wntqa , wntqas , wntqn , wntqo , wntqs
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
      wntqa = LSAME(Jobz,'A')
      wntqs = LSAME(Jobz,'S')
      wntqas = wntqa .OR. wntqs
      wntqo = LSAME(Jobz,'O')
      wntqn = LSAME(Jobz,'N')
      lquery = (Lwork==-1)
!
      IF ( .NOT.(wntqa .OR. wntqs .OR. wntqo .OR. wntqn) ) THEN
         Info = -1
      ELSEIF ( M<0 ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ELSEIF ( Ldu<1 .OR. (wntqas .AND. Ldu<M) .OR.                     &
     &         (wntqo .AND. M<N .AND. Ldu<M) ) THEN
         Info = -8
      ELSEIF ( Ldvt<1 .OR. (wntqa .AND. Ldvt<N) .OR.                    &
     &         (wntqs .AND. Ldvt<minmn) .OR.                            &
     &         (wntqo .AND. M>=N .AND. Ldvt<N) ) THEN
         Info = -10
      ENDIF
!
!     Compute workspace
!       Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace allocated at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.
!
      IF ( Info==0 ) THEN
         minwrk = 1
         maxwrk = 1
         bdspac = 0
         mnthr = INT(minmn*11.0D0/6.0D0)
         IF ( M>=N .AND. minmn>0 ) THEN
!
!           Compute space needed for DBDSDC
!
            IF ( wntqn ) THEN
!              dbdsdc needs only 4*N (or 6*N for uplo=L for LAPACK <= 3.6)
!              keep 7*N for backwards compatibility.
               bdspac = 7*N
            ELSE
               bdspac = 3*N*N + 4*N
            ENDIF
!
!           Compute space preferred for each routine
            CALL DGEBRD(M,N,dum(1),M,dum(1),dum(1),dum(1),dum(1),dum(1),&
     &                  -1,ierr)
            lwork_dgebrd_mn = INT(dum(1))
!
            CALL DGEBRD(N,N,dum(1),N,dum(1),dum(1),dum(1),dum(1),dum(1),&
     &                  -1,ierr)
            lwork_dgebrd_nn = INT(dum(1))
!
            CALL DGEQRF(M,N,dum(1),M,dum(1),dum(1),-1,ierr)
            lwork_dgeqrf_mn = INT(dum(1))
!
            CALL DORGBR('Q',N,N,N,dum(1),N,dum(1),dum(1),-1,ierr)
            lwork_dorgbr_q_nn = INT(dum(1))
!
            CALL DORGQR(M,M,N,dum(1),M,dum(1),dum(1),-1,ierr)
            lwork_dorgqr_mm = INT(dum(1))
!
            CALL DORGQR(M,N,N,dum(1),M,dum(1),dum(1),-1,ierr)
            lwork_dorgqr_mn = INT(dum(1))
!
            CALL DORMBR('P','R','T',N,N,N,dum(1),N,dum(1),dum(1),N,     &
     &                  dum(1),-1,ierr)
            lwork_dormbr_prt_nn = INT(dum(1))
!
            CALL DORMBR('Q','L','N',N,N,N,dum(1),N,dum(1),dum(1),N,     &
     &                  dum(1),-1,ierr)
            lwork_dormbr_qln_nn = INT(dum(1))
!
            CALL DORMBR('Q','L','N',M,N,N,dum(1),M,dum(1),dum(1),M,     &
     &                  dum(1),-1,ierr)
            lwork_dormbr_qln_mn = INT(dum(1))
!
            CALL DORMBR('Q','L','N',M,M,N,dum(1),M,dum(1),dum(1),M,     &
     &                  dum(1),-1,ierr)
            lwork_dormbr_qln_mm = INT(dum(1))
!
            IF ( M<mnthr ) THEN
!
!              Path 5 (M >= N, but not much larger)
!
               wrkbl = 3*N + lwork_dgebrd_mn
               IF ( wntqn ) THEN
!                 Path 5n (M >= N, jobz='N')
                  maxwrk = MAX(wrkbl,3*N+bdspac)
                  minwrk = 3*N + MAX(M,bdspac)
               ELSEIF ( wntqo ) THEN
!                 Path 5o (M >= N, jobz='O')
                  wrkbl = MAX(wrkbl,3*N+lwork_dormbr_prt_nn)
                  wrkbl = MAX(wrkbl,3*N+lwork_dormbr_qln_mn)
                  wrkbl = MAX(wrkbl,3*N+bdspac)
                  maxwrk = wrkbl + M*N
                  minwrk = 3*N + MAX(M,N*N+bdspac)
               ELSEIF ( wntqs ) THEN
!                 Path 5s (M >= N, jobz='S')
                  wrkbl = MAX(wrkbl,3*N+lwork_dormbr_qln_mn)
                  wrkbl = MAX(wrkbl,3*N+lwork_dormbr_prt_nn)
                  maxwrk = MAX(wrkbl,3*N+bdspac)
                  minwrk = 3*N + MAX(M,bdspac)
               ELSEIF ( wntqa ) THEN
!                 Path 5a (M >= N, jobz='A')
                  wrkbl = MAX(wrkbl,3*N+lwork_dormbr_qln_mm)
                  wrkbl = MAX(wrkbl,3*N+lwork_dormbr_prt_nn)
                  maxwrk = MAX(wrkbl,3*N+bdspac)
                  minwrk = 3*N + MAX(M,bdspac)
               ENDIF
            ELSEIF ( wntqn ) THEN
!
!                 Path 1 (M >> N, JOBZ='N')
!
               wrkbl = N + lwork_dgeqrf_mn
               wrkbl = MAX(wrkbl,3*N+lwork_dgebrd_nn)
               maxwrk = MAX(wrkbl,bdspac+N)
               minwrk = bdspac + N
            ELSEIF ( wntqo ) THEN
!
!                 Path 2 (M >> N, JOBZ='O')
!
               wrkbl = N + lwork_dgeqrf_mn
               wrkbl = MAX(wrkbl,N+lwork_dorgqr_mn)
               wrkbl = MAX(wrkbl,3*N+lwork_dgebrd_nn)
               wrkbl = MAX(wrkbl,3*N+lwork_dormbr_qln_nn)
               wrkbl = MAX(wrkbl,3*N+lwork_dormbr_prt_nn)
               wrkbl = MAX(wrkbl,3*N+bdspac)
               maxwrk = wrkbl + 2*N*N
               minwrk = bdspac + 2*N*N + 3*N
            ELSEIF ( wntqs ) THEN
!
!                 Path 3 (M >> N, JOBZ='S')
!
               wrkbl = N + lwork_dgeqrf_mn
               wrkbl = MAX(wrkbl,N+lwork_dorgqr_mn)
               wrkbl = MAX(wrkbl,3*N+lwork_dgebrd_nn)
               wrkbl = MAX(wrkbl,3*N+lwork_dormbr_qln_nn)
               wrkbl = MAX(wrkbl,3*N+lwork_dormbr_prt_nn)
               wrkbl = MAX(wrkbl,3*N+bdspac)
               maxwrk = wrkbl + N*N
               minwrk = bdspac + N*N + 3*N
            ELSEIF ( wntqa ) THEN
!
!                 Path 4 (M >> N, JOBZ='A')
!
               wrkbl = N + lwork_dgeqrf_mn
               wrkbl = MAX(wrkbl,N+lwork_dorgqr_mm)
               wrkbl = MAX(wrkbl,3*N+lwork_dgebrd_nn)
               wrkbl = MAX(wrkbl,3*N+lwork_dormbr_qln_nn)
               wrkbl = MAX(wrkbl,3*N+lwork_dormbr_prt_nn)
               wrkbl = MAX(wrkbl,3*N+bdspac)
               maxwrk = wrkbl + N*N
               minwrk = N*N + MAX(3*N+bdspac,N+M)
            ENDIF
         ELSEIF ( minmn>0 ) THEN
!
!           Compute space needed for DBDSDC
!
            IF ( wntqn ) THEN
!              dbdsdc needs only 4*N (or 6*N for uplo=L for LAPACK <= 3.6)
!              keep 7*N for backwards compatibility.
               bdspac = 7*M
            ELSE
               bdspac = 3*M*M + 4*M
            ENDIF
!
!           Compute space preferred for each routine
            CALL DGEBRD(M,N,dum(1),M,dum(1),dum(1),dum(1),dum(1),dum(1),&
     &                  -1,ierr)
            lwork_dgebrd_mn = INT(dum(1))
!
            CALL DGEBRD(M,M,A,M,S,dum(1),dum(1),dum(1),dum(1),-1,ierr)
            lwork_dgebrd_mm = INT(dum(1))
!
            CALL DGELQF(M,N,A,M,dum(1),dum(1),-1,ierr)
            lwork_dgelqf_mn = INT(dum(1))
!
            CALL DORGLQ(N,N,M,dum(1),N,dum(1),dum(1),-1,ierr)
            lwork_dorglq_nn = INT(dum(1))
!
            CALL DORGLQ(M,N,M,A,M,dum(1),dum(1),-1,ierr)
            lwork_dorglq_mn = INT(dum(1))
!
            CALL DORGBR('P',M,M,M,A,N,dum(1),dum(1),-1,ierr)
            lwork_dorgbr_p_mm = INT(dum(1))
!
            CALL DORMBR('P','R','T',M,M,M,dum(1),M,dum(1),dum(1),M,     &
     &                  dum(1),-1,ierr)
            lwork_dormbr_prt_mm = INT(dum(1))
!
            CALL DORMBR('P','R','T',M,N,M,dum(1),M,dum(1),dum(1),M,     &
     &                  dum(1),-1,ierr)
            lwork_dormbr_prt_mn = INT(dum(1))
!
            CALL DORMBR('P','R','T',N,N,M,dum(1),N,dum(1),dum(1),N,     &
     &                  dum(1),-1,ierr)
            lwork_dormbr_prt_nn = INT(dum(1))
!
            CALL DORMBR('Q','L','N',M,M,M,dum(1),M,dum(1),dum(1),M,     &
     &                  dum(1),-1,ierr)
            lwork_dormbr_qln_mm = INT(dum(1))
!
            IF ( N<mnthr ) THEN
!
!              Path 5t (N > M, but not much larger)
!
               wrkbl = 3*M + lwork_dgebrd_mn
               IF ( wntqn ) THEN
!                 Path 5tn (N > M, jobz='N')
                  maxwrk = MAX(wrkbl,3*M+bdspac)
                  minwrk = 3*M + MAX(N,bdspac)
               ELSEIF ( wntqo ) THEN
!                 Path 5to (N > M, jobz='O')
                  wrkbl = MAX(wrkbl,3*M+lwork_dormbr_qln_mm)
                  wrkbl = MAX(wrkbl,3*M+lwork_dormbr_prt_mn)
                  wrkbl = MAX(wrkbl,3*M+bdspac)
                  maxwrk = wrkbl + M*N
                  minwrk = 3*M + MAX(N,M*M+bdspac)
               ELSEIF ( wntqs ) THEN
!                 Path 5ts (N > M, jobz='S')
                  wrkbl = MAX(wrkbl,3*M+lwork_dormbr_qln_mm)
                  wrkbl = MAX(wrkbl,3*M+lwork_dormbr_prt_mn)
                  maxwrk = MAX(wrkbl,3*M+bdspac)
                  minwrk = 3*M + MAX(N,bdspac)
               ELSEIF ( wntqa ) THEN
!                 Path 5ta (N > M, jobz='A')
                  wrkbl = MAX(wrkbl,3*M+lwork_dormbr_qln_mm)
                  wrkbl = MAX(wrkbl,3*M+lwork_dormbr_prt_nn)
                  maxwrk = MAX(wrkbl,3*M+bdspac)
                  minwrk = 3*M + MAX(N,bdspac)
               ENDIF
            ELSEIF ( wntqn ) THEN
!
!                 Path 1t (N >> M, JOBZ='N')
!
               wrkbl = M + lwork_dgelqf_mn
               wrkbl = MAX(wrkbl,3*M+lwork_dgebrd_mm)
               maxwrk = MAX(wrkbl,bdspac+M)
               minwrk = bdspac + M
            ELSEIF ( wntqo ) THEN
!
!                 Path 2t (N >> M, JOBZ='O')
!
               wrkbl = M + lwork_dgelqf_mn
               wrkbl = MAX(wrkbl,M+lwork_dorglq_mn)
               wrkbl = MAX(wrkbl,3*M+lwork_dgebrd_mm)
               wrkbl = MAX(wrkbl,3*M+lwork_dormbr_qln_mm)
               wrkbl = MAX(wrkbl,3*M+lwork_dormbr_prt_mm)
               wrkbl = MAX(wrkbl,3*M+bdspac)
               maxwrk = wrkbl + 2*M*M
               minwrk = bdspac + 2*M*M + 3*M
            ELSEIF ( wntqs ) THEN
!
!                 Path 3t (N >> M, JOBZ='S')
!
               wrkbl = M + lwork_dgelqf_mn
               wrkbl = MAX(wrkbl,M+lwork_dorglq_mn)
               wrkbl = MAX(wrkbl,3*M+lwork_dgebrd_mm)
               wrkbl = MAX(wrkbl,3*M+lwork_dormbr_qln_mm)
               wrkbl = MAX(wrkbl,3*M+lwork_dormbr_prt_mm)
               wrkbl = MAX(wrkbl,3*M+bdspac)
               maxwrk = wrkbl + M*M
               minwrk = bdspac + M*M + 3*M
            ELSEIF ( wntqa ) THEN
!
!                 Path 4t (N >> M, JOBZ='A')
!
               wrkbl = M + lwork_dgelqf_mn
               wrkbl = MAX(wrkbl,M+lwork_dorglq_nn)
               wrkbl = MAX(wrkbl,3*M+lwork_dgebrd_mm)
               wrkbl = MAX(wrkbl,3*M+lwork_dormbr_qln_mm)
               wrkbl = MAX(wrkbl,3*M+lwork_dormbr_prt_mm)
               wrkbl = MAX(wrkbl,3*M+bdspac)
               maxwrk = wrkbl + M*M
               minwrk = M*M + MAX(3*M+bdspac,M+N)
            ENDIF
         ENDIF
 
         maxwrk = MAX(maxwrk,minwrk)
         Work(1) = maxwrk
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -12
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGESDD',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) RETURN
!
!     Get machine constants
!
      eps = DLAMCH('P')
      smlnum = SQRT(DLAMCH('S'))/eps
      bignum = ONE/smlnum
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      anrm = DLANGE('M',M,N,A,Lda,dum)
      IF ( DISNAN(anrm) ) THEN
         Info = -4
         RETURN
      ENDIF
      iscl = 0
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         iscl = 1
         CALL DLASCL('G',0,0,anrm,smlnum,M,N,A,Lda,ierr)
      ELSEIF ( anrm>bignum ) THEN
         iscl = 1
         CALL DLASCL('G',0,0,anrm,bignum,M,N,A,Lda,ierr)
      ENDIF
!
      IF ( M>=N ) THEN
!
!        A has at least as many rows as columns. If A has sufficiently
!        more rows than columns, first reduce using the QR
!        decomposition (if sufficient workspace available)
!
         IF ( M<mnthr ) THEN
!
!           M .LT. MNTHR
!
!           Path 5 (M >= N, but not much larger)
!           Reduce to bidiagonal form without QR decomposition
!
            ie = 1
            itauq = ie + N
            itaup = itauq + N
            nwork = itaup + N
!
!           Bidiagonalize A
!           Workspace: need   3*N [e, tauq, taup] + M        [work]
!           Workspace: prefer 3*N [e, tauq, taup] + (M+N)*NB [work]
!
            CALL DGEBRD(M,N,A,Lda,S,Work(ie),Work(itauq),Work(itaup),   &
     &                  Work(nwork),Lwork-nwork+1,ierr)
            IF ( wntqn ) THEN
!
!              Path 5n (M >= N, JOBZ='N')
!              Perform bidiagonal SVD, only computing singular values
!              Workspace: need   3*N [e, tauq, taup] + BDSPAC
!
               CALL DBDSDC('U','N',N,S,Work(ie),dum,1,dum,1,dum,idum,   &
     &                     Work(nwork),Iwork,Info)
            ELSEIF ( wntqo ) THEN
!              Path 5o (M >= N, JOBZ='O')
               iu = nwork
               IF ( Lwork>=M*N+3*N+bdspac ) THEN
!
!                 WORK( IU ) is M by N
!
                  ldwrku = M
                  nwork = iu + ldwrku*N
                  CALL DLASET('F',M,N,ZERO,ZERO,Work(iu),ldwrku)
!                 IR is unused; silence compile warnings
                  ir = -1
               ELSE
!
!                 WORK( IU ) is N by N
!
                  ldwrku = N
                  nwork = iu + ldwrku*N
!
!                 WORK(IR) is LDWRKR by N
!
                  ir = nwork
                  ldwrkr = (Lwork-N*N-3*N)/N
               ENDIF
               nwork = iu + ldwrku*N
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in WORK(IU) and computing right
!              singular vectors of bidiagonal matrix in VT
!              Workspace: need   3*N [e, tauq, taup] + N*N [U] + BDSPAC
!
               CALL DBDSDC('U','I',N,S,Work(ie),Work(iu),ldwrku,Vt,Ldvt,&
     &                     dum,idum,Work(nwork),Iwork,Info)
!
!              Overwrite VT by right singular vectors of A
!              Workspace: need   3*N [e, tauq, taup] + N*N [U] + N    [work]
!              Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + N*NB [work]
!
               CALL DORMBR('P','R','T',N,N,N,A,Lda,Work(itaup),Vt,Ldvt, &
     &                     Work(nwork),Lwork-nwork+1,ierr)
!
               IF ( Lwork>=M*N+3*N+bdspac ) THEN
!
!                 Path 5o-fast
!                 Overwrite WORK(IU) by left singular vectors of A
!                 Workspace: need   3*N [e, tauq, taup] + M*N [U] + N    [work]
!                 Workspace: prefer 3*N [e, tauq, taup] + M*N [U] + N*NB [work]
!
                  CALL DORMBR('Q','L','N',M,N,N,A,Lda,Work(itauq),      &
     &                        Work(iu),ldwrku,Work(nwork),Lwork-nwork+1,&
     &                        ierr)
!
!                 Copy left singular vectors of A from WORK(IU) to A
!
                  CALL DLACPY('F',M,N,Work(iu),ldwrku,A,Lda)
               ELSE
!
!                 Path 5o-slow
!                 Generate Q in A
!                 Workspace: need   3*N [e, tauq, taup] + N*N [U] + N    [work]
!                 Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + N*NB [work]
!
                  CALL DORGBR('Q',M,N,N,A,Lda,Work(itauq),Work(nwork),  &
     &                        Lwork-nwork+1,ierr)
!
!                 Multiply Q in A by left singular vectors of
!                 bidiagonal matrix in WORK(IU), storing result in
!                 WORK(IR) and copying to A
!                 Workspace: need   3*N [e, tauq, taup] + N*N [U] + NB*N [R]
!                 Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + M*N  [R]
!
                  DO i = 1 , M , ldwrkr
                     chunk = MIN(M-i+1,ldwrkr)
                     CALL DGEMM('N','N',chunk,N,N,ONE,A(i,1),Lda,       &
     &                          Work(iu),ldwrku,ZERO,Work(ir),ldwrkr)
                     CALL DLACPY('F',chunk,N,Work(ir),ldwrkr,A(i,1),Lda)
                  ENDDO
               ENDIF
!
            ELSEIF ( wntqs ) THEN
!
!              Path 5s (M >= N, JOBZ='S')
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in VT
!              Workspace: need   3*N [e, tauq, taup] + BDSPAC
!
               CALL DLASET('F',M,N,ZERO,ZERO,U,Ldu)
               CALL DBDSDC('U','I',N,S,Work(ie),U,Ldu,Vt,Ldvt,dum,idum, &
     &                     Work(nwork),Iwork,Info)
!
!              Overwrite U by left singular vectors of A and VT
!              by right singular vectors of A
!              Workspace: need   3*N [e, tauq, taup] + N    [work]
!              Workspace: prefer 3*N [e, tauq, taup] + N*NB [work]
!
               CALL DORMBR('Q','L','N',M,N,N,A,Lda,Work(itauq),U,Ldu,   &
     &                     Work(nwork),Lwork-nwork+1,ierr)
               CALL DORMBR('P','R','T',N,N,N,A,Lda,Work(itaup),Vt,Ldvt, &
     &                     Work(nwork),Lwork-nwork+1,ierr)
            ELSEIF ( wntqa ) THEN
!
!              Path 5a (M >= N, JOBZ='A')
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in VT
!              Workspace: need   3*N [e, tauq, taup] + BDSPAC
!
               CALL DLASET('F',M,M,ZERO,ZERO,U,Ldu)
               CALL DBDSDC('U','I',N,S,Work(ie),U,Ldu,Vt,Ldvt,dum,idum, &
     &                     Work(nwork),Iwork,Info)
!
!              Set the right corner of U to identity matrix
!
               IF ( M>N ) CALL DLASET('F',M-N,M-N,ZERO,ONE,U(N+1,N+1),  &
     &                                Ldu)
!
!              Overwrite U by left singular vectors of A and VT
!              by right singular vectors of A
!              Workspace: need   3*N [e, tauq, taup] + M    [work]
!              Workspace: prefer 3*N [e, tauq, taup] + M*NB [work]
!
               CALL DORMBR('Q','L','N',M,M,N,A,Lda,Work(itauq),U,Ldu,   &
     &                     Work(nwork),Lwork-nwork+1,ierr)
               CALL DORMBR('P','R','T',N,N,M,A,Lda,Work(itaup),Vt,Ldvt, &
     &                     Work(nwork),Lwork-nwork+1,ierr)
            ENDIF
!
         ELSEIF ( wntqn ) THEN
!
!              Path 1 (M >> N, JOBZ='N')
!              No singular vectors to be computed
!
            itau = 1
            nwork = itau + N
!
!              Compute A=Q*R
!              Workspace: need   N [tau] + N    [work]
!              Workspace: prefer N [tau] + N*NB [work]
!
            CALL DGEQRF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  ierr)
!
!              Zero out below R
!
            CALL DLASET('L',N-1,N-1,ZERO,ZERO,A(2,1),Lda)
            ie = 1
            itauq = ie + N
            itaup = itauq + N
            nwork = itaup + N
!
!              Bidiagonalize R in A
!              Workspace: need   3*N [e, tauq, taup] + N      [work]
!              Workspace: prefer 3*N [e, tauq, taup] + 2*N*NB [work]
!
            CALL DGEBRD(N,N,A,Lda,S,Work(ie),Work(itauq),Work(itaup),   &
     &                  Work(nwork),Lwork-nwork+1,ierr)
            nwork = ie + N
!
!              Perform bidiagonal SVD, computing singular values only
!              Workspace: need   N [e] + BDSPAC
!
            CALL DBDSDC('U','N',N,S,Work(ie),dum,1,dum,1,dum,idum,      &
     &                  Work(nwork),Iwork,Info)
!
         ELSEIF ( wntqo ) THEN
!
!              Path 2 (M >> N, JOBZ = 'O')
!              N left singular vectors to be overwritten on A and
!              N right singular vectors to be computed in VT
!
            ir = 1
!
!              WORK(IR) is LDWRKR by N
!
            IF ( Lwork>=Lda*N+N*N+3*N+bdspac ) THEN
               ldwrkr = Lda
            ELSE
               ldwrkr = (Lwork-N*N-3*N-bdspac)/N
            ENDIF
            itau = ir + ldwrkr*N
            nwork = itau + N
!
!              Compute A=Q*R
!              Workspace: need   N*N [R] + N [tau] + N    [work]
!              Workspace: prefer N*N [R] + N [tau] + N*NB [work]
!
            CALL DGEQRF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  ierr)
!
!              Copy R to WORK(IR), zeroing out below it
!
            CALL DLACPY('U',N,N,A,Lda,Work(ir),ldwrkr)
            CALL DLASET('L',N-1,N-1,ZERO,ZERO,Work(ir+1),ldwrkr)
!
!              Generate Q in A
!              Workspace: need   N*N [R] + N [tau] + N    [work]
!              Workspace: prefer N*N [R] + N [tau] + N*NB [work]
!
            CALL DORGQR(M,N,N,A,Lda,Work(itau),Work(nwork),             &
     &                  Lwork-nwork+1,ierr)
            ie = itau
            itauq = ie + N
            itaup = itauq + N
            nwork = itaup + N
!
!              Bidiagonalize R in WORK(IR)
!              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N      [work]
!              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + 2*N*NB [work]
!
            CALL DGEBRD(N,N,Work(ir),ldwrkr,S,Work(ie),Work(itauq),     &
     &                  Work(itaup),Work(nwork),Lwork-nwork+1,ierr)
!
!              WORK(IU) is N by N
!
            iu = nwork
            nwork = iu + N*N
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in WORK(IU) and computing right
!              singular vectors of bidiagonal matrix in VT
!              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] + BDSPAC
!
            CALL DBDSDC('U','I',N,S,Work(ie),Work(iu),N,Vt,Ldvt,dum,    &
     &                  idum,Work(nwork),Iwork,Info)
!
!              Overwrite WORK(IU) by left singular vectors of R
!              and VT by right singular vectors of R
!              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] + N    [work]
!              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + N*N [U] + N*NB [work]
!
            CALL DORMBR('Q','L','N',N,N,N,Work(ir),ldwrkr,Work(itauq),  &
     &                  Work(iu),N,Work(nwork),Lwork-nwork+1,ierr)
            CALL DORMBR('P','R','T',N,N,N,Work(ir),ldwrkr,Work(itaup),  &
     &                  Vt,Ldvt,Work(nwork),Lwork-nwork+1,ierr)
!
!              Multiply Q in A by left singular vectors of R in
!              WORK(IU), storing result in WORK(IR) and copying to A
!              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U]
!              Workspace: prefer M*N [R] + 3*N [e, tauq, taup] + N*N [U]
!
            DO i = 1 , M , ldwrkr
               chunk = MIN(M-i+1,ldwrkr)
               CALL DGEMM('N','N',chunk,N,N,ONE,A(i,1),Lda,Work(iu),N,  &
     &                    ZERO,Work(ir),ldwrkr)
               CALL DLACPY('F',chunk,N,Work(ir),ldwrkr,A(i,1),Lda)
            ENDDO
!
         ELSEIF ( wntqs ) THEN
!
!              Path 3 (M >> N, JOBZ='S')
!              N left singular vectors to be computed in U and
!              N right singular vectors to be computed in VT
!
            ir = 1
!
!              WORK(IR) is N by N
!
            ldwrkr = N
            itau = ir + ldwrkr*N
            nwork = itau + N
!
!              Compute A=Q*R
!              Workspace: need   N*N [R] + N [tau] + N    [work]
!              Workspace: prefer N*N [R] + N [tau] + N*NB [work]
!
            CALL DGEQRF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  ierr)
!
!              Copy R to WORK(IR), zeroing out below it
!
            CALL DLACPY('U',N,N,A,Lda,Work(ir),ldwrkr)
            CALL DLASET('L',N-1,N-1,ZERO,ZERO,Work(ir+1),ldwrkr)
!
!              Generate Q in A
!              Workspace: need   N*N [R] + N [tau] + N    [work]
!              Workspace: prefer N*N [R] + N [tau] + N*NB [work]
!
            CALL DORGQR(M,N,N,A,Lda,Work(itau),Work(nwork),             &
     &                  Lwork-nwork+1,ierr)
            ie = itau
            itauq = ie + N
            itaup = itauq + N
            nwork = itaup + N
!
!              Bidiagonalize R in WORK(IR)
!              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N      [work]
!              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + 2*N*NB [work]
!
            CALL DGEBRD(N,N,Work(ir),ldwrkr,S,Work(ie),Work(itauq),     &
     &                  Work(itaup),Work(nwork),Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagoal matrix in U and computing right singular
!              vectors of bidiagonal matrix in VT
!              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + BDSPAC
!
            CALL DBDSDC('U','I',N,S,Work(ie),U,Ldu,Vt,Ldvt,dum,idum,    &
     &                  Work(nwork),Iwork,Info)
!
!              Overwrite U by left singular vectors of R and VT
!              by right singular vectors of R
!              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N    [work]
!              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + N*NB [work]
!
            CALL DORMBR('Q','L','N',N,N,N,Work(ir),ldwrkr,Work(itauq),U,&
     &                  Ldu,Work(nwork),Lwork-nwork+1,ierr)
!
            CALL DORMBR('P','R','T',N,N,N,Work(ir),ldwrkr,Work(itaup),  &
     &                  Vt,Ldvt,Work(nwork),Lwork-nwork+1,ierr)
!
!              Multiply Q in A by left singular vectors of R in
!              WORK(IR), storing result in U
!              Workspace: need   N*N [R]
!
            CALL DLACPY('F',N,N,U,Ldu,Work(ir),ldwrkr)
            CALL DGEMM('N','N',M,N,N,ONE,A,Lda,Work(ir),ldwrkr,ZERO,U,  &
     &                 Ldu)
!
         ELSEIF ( wntqa ) THEN
!
!              Path 4 (M >> N, JOBZ='A')
!              M left singular vectors to be computed in U and
!              N right singular vectors to be computed in VT
!
            iu = 1
!
!              WORK(IU) is N by N
!
            ldwrku = N
            itau = iu + ldwrku*N
            nwork = itau + N
!
!              Compute A=Q*R, copying result to U
!              Workspace: need   N*N [U] + N [tau] + N    [work]
!              Workspace: prefer N*N [U] + N [tau] + N*NB [work]
!
            CALL DGEQRF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  ierr)
            CALL DLACPY('L',M,N,A,Lda,U,Ldu)
!
!              Generate Q in U
!              Workspace: need   N*N [U] + N [tau] + M    [work]
!              Workspace: prefer N*N [U] + N [tau] + M*NB [work]
            CALL DORGQR(M,M,N,U,Ldu,Work(itau),Work(nwork),             &
     &                  Lwork-nwork+1,ierr)
!
!              Produce R in A, zeroing out other entries
!
            CALL DLASET('L',N-1,N-1,ZERO,ZERO,A(2,1),Lda)
            ie = itau
            itauq = ie + N
            itaup = itauq + N
            nwork = itaup + N
!
!              Bidiagonalize R in A
!              Workspace: need   N*N [U] + 3*N [e, tauq, taup] + N      [work]
!              Workspace: prefer N*N [U] + 3*N [e, tauq, taup] + 2*N*NB [work]
!
            CALL DGEBRD(N,N,A,Lda,S,Work(ie),Work(itauq),Work(itaup),   &
     &                  Work(nwork),Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in WORK(IU) and computing right
!              singular vectors of bidiagonal matrix in VT
!              Workspace: need   N*N [U] + 3*N [e, tauq, taup] + BDSPAC
!
            CALL DBDSDC('U','I',N,S,Work(ie),Work(iu),N,Vt,Ldvt,dum,    &
     &                  idum,Work(nwork),Iwork,Info)
!
!              Overwrite WORK(IU) by left singular vectors of R and VT
!              by right singular vectors of R
!              Workspace: need   N*N [U] + 3*N [e, tauq, taup] + N    [work]
!              Workspace: prefer N*N [U] + 3*N [e, tauq, taup] + N*NB [work]
!
            CALL DORMBR('Q','L','N',N,N,N,A,Lda,Work(itauq),Work(iu),   &
     &                  ldwrku,Work(nwork),Lwork-nwork+1,ierr)
            CALL DORMBR('P','R','T',N,N,N,A,Lda,Work(itaup),Vt,Ldvt,    &
     &                  Work(nwork),Lwork-nwork+1,ierr)
!
!              Multiply Q in U by left singular vectors of R in
!              WORK(IU), storing result in A
!              Workspace: need   N*N [U]
!
            CALL DGEMM('N','N',M,N,N,ONE,U,Ldu,Work(iu),ldwrku,ZERO,A,  &
     &                 Lda)
!
!              Copy left singular vectors of A from A to U
!
            CALL DLACPY('F',M,N,A,Lda,U,Ldu)
!
         ENDIF
!
!
!        A has more columns than rows. If A has sufficiently more
!        columns than rows, first reduce using the LQ decomposition (if
!        sufficient workspace available)
!
      ELSEIF ( N>=mnthr ) THEN
!
         IF ( wntqn ) THEN
!
!              Path 1t (N >> M, JOBZ='N')
!              No singular vectors to be computed
!
            itau = 1
!
!
            nwork = itau + M
!
!              Compute A=L*Q
!              Workspace: need   M [tau] + M [work]
!              Workspace: prefer M [tau] + M*NB [work]
!
            CALL DGELQF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  ierr)
!
!              Zero out above L
!
            CALL DLASET('U',M-1,M-1,ZERO,ZERO,A(1,2),Lda)
            ie = 1
            itauq = ie + M
            itaup = itauq + M
            nwork = itaup + M
!
!              Bidiagonalize L in A
!              Workspace: need   3*M [e, tauq, taup] + M      [work]
!              Workspace: prefer 3*M [e, tauq, taup] + 2*M*NB [work]
!
            CALL DGEBRD(M,M,A,Lda,S,Work(ie),Work(itauq),Work(itaup),   &
     &                  Work(nwork),Lwork-nwork+1,ierr)
            nwork = ie + M
!
!              Perform bidiagonal SVD, computing singular values only
!              Workspace: need   M [e] + BDSPAC
!
            CALL DBDSDC('U','N',M,S,Work(ie),dum,1,dum,1,dum,idum,      &
     &                  Work(nwork),Iwork,Info)
!
         ELSEIF ( wntqo ) THEN
!
!              Path 2t (N >> M, JOBZ='O')
!              M right singular vectors to be overwritten on A and
!              M left singular vectors to be computed in U
!
            ivt = 1
!
!              WORK(IVT) is M by M
!              WORK(IL)  is M by M; it is later resized to M by chunk for gemm
!
            il = ivt + M*M
            IF ( Lwork>=M*N+M*M+3*M+bdspac ) THEN
               ldwrkl = M
               chunk = N
            ELSE
               ldwrkl = M
               chunk = (Lwork-M*M)/M
            ENDIF
            itau = il + ldwrkl*M
            nwork = itau + M
!
!              Compute A=L*Q
!              Workspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
!              Workspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]
!
            CALL DGELQF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  ierr)
!
!              Copy L to WORK(IL), zeroing about above it
!
            CALL DLACPY('L',M,M,A,Lda,Work(il),ldwrkl)
            CALL DLASET('U',M-1,M-1,ZERO,ZERO,Work(il+ldwrkl),ldwrkl)
!
!              Generate Q in A
!              Workspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
!              Workspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]
!
            CALL DORGLQ(M,N,M,A,Lda,Work(itau),Work(nwork),             &
     &                  Lwork-nwork+1,ierr)
            ie = itau
            itauq = ie + M
            itaup = itauq + M
            nwork = itaup + M
!
!              Bidiagonalize L in WORK(IL)
!              Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M      [work]
!              Workspace: prefer M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + 2*M*NB [work]
!
            CALL DGEBRD(M,M,Work(il),ldwrkl,S,Work(ie),Work(itauq),     &
     &                  Work(itaup),Work(nwork),Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U, and computing right singular
!              vectors of bidiagonal matrix in WORK(IVT)
!              Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + BDSPAC
!
            CALL DBDSDC('U','I',M,S,Work(ie),U,Ldu,Work(ivt),M,dum,idum,&
     &                  Work(nwork),Iwork,Info)
!
!              Overwrite U by left singular vectors of L and WORK(IVT)
!              by right singular vectors of L
!              Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M    [work]
!              Workspace: prefer M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M*NB [work]
!
            CALL DORMBR('Q','L','N',M,M,M,Work(il),ldwrkl,Work(itauq),U,&
     &                  Ldu,Work(nwork),Lwork-nwork+1,ierr)
            CALL DORMBR('P','R','T',M,M,M,Work(il),ldwrkl,Work(itaup),  &
     &                  Work(ivt),M,Work(nwork),Lwork-nwork+1,ierr)
!
!              Multiply right singular vectors of L in WORK(IVT) by Q
!              in A, storing result in WORK(IL) and copying to A
!              Workspace: need   M*M [VT] + M*M [L]
!              Workspace: prefer M*M [VT] + M*N [L]
!              At this point, L is resized as M by chunk.
!
            DO i = 1 , N , chunk
               blk = MIN(N-i+1,chunk)
               CALL DGEMM('N','N',M,blk,M,ONE,Work(ivt),M,A(1,i),Lda,   &
     &                    ZERO,Work(il),ldwrkl)
               CALL DLACPY('F',M,blk,Work(il),ldwrkl,A(1,i),Lda)
            ENDDO
!
         ELSEIF ( wntqs ) THEN
!
!              Path 3t (N >> M, JOBZ='S')
!              M right singular vectors to be computed in VT and
!              M left singular vectors to be computed in U
!
            il = 1
!
!              WORK(IL) is M by M
!
            ldwrkl = M
            itau = il + ldwrkl*M
            nwork = itau + M
!
!              Compute A=L*Q
!              Workspace: need   M*M [L] + M [tau] + M    [work]
!              Workspace: prefer M*M [L] + M [tau] + M*NB [work]
!
            CALL DGELQF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  ierr)
!
!              Copy L to WORK(IL), zeroing out above it
!
            CALL DLACPY('L',M,M,A,Lda,Work(il),ldwrkl)
            CALL DLASET('U',M-1,M-1,ZERO,ZERO,Work(il+ldwrkl),ldwrkl)
!
!              Generate Q in A
!              Workspace: need   M*M [L] + M [tau] + M    [work]
!              Workspace: prefer M*M [L] + M [tau] + M*NB [work]
!
            CALL DORGLQ(M,N,M,A,Lda,Work(itau),Work(nwork),             &
     &                  Lwork-nwork+1,ierr)
            ie = itau
            itauq = ie + M
            itaup = itauq + M
            nwork = itaup + M
!
!              Bidiagonalize L in WORK(IU).
!              Workspace: need   M*M [L] + 3*M [e, tauq, taup] + M      [work]
!              Workspace: prefer M*M [L] + 3*M [e, tauq, taup] + 2*M*NB [work]
!
            CALL DGEBRD(M,M,Work(il),ldwrkl,S,Work(ie),Work(itauq),     &
     &                  Work(itaup),Work(nwork),Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in VT
!              Workspace: need   M*M [L] + 3*M [e, tauq, taup] + BDSPAC
!
            CALL DBDSDC('U','I',M,S,Work(ie),U,Ldu,Vt,Ldvt,dum,idum,    &
     &                  Work(nwork),Iwork,Info)
!
!              Overwrite U by left singular vectors of L and VT
!              by right singular vectors of L
!              Workspace: need   M*M [L] + 3*M [e, tauq, taup] + M    [work]
!              Workspace: prefer M*M [L] + 3*M [e, tauq, taup] + M*NB [work]
!
            CALL DORMBR('Q','L','N',M,M,M,Work(il),ldwrkl,Work(itauq),U,&
     &                  Ldu,Work(nwork),Lwork-nwork+1,ierr)
            CALL DORMBR('P','R','T',M,M,M,Work(il),ldwrkl,Work(itaup),  &
     &                  Vt,Ldvt,Work(nwork),Lwork-nwork+1,ierr)
!
!              Multiply right singular vectors of L in WORK(IL) by
!              Q in A, storing result in VT
!              Workspace: need   M*M [L]
!
            CALL DLACPY('F',M,M,Vt,Ldvt,Work(il),ldwrkl)
            CALL DGEMM('N','N',M,N,M,ONE,Work(il),ldwrkl,A,Lda,ZERO,Vt, &
     &                 Ldvt)
!
         ELSEIF ( wntqa ) THEN
!
!              Path 4t (N >> M, JOBZ='A')
!              N right singular vectors to be computed in VT and
!              M left singular vectors to be computed in U
!
            ivt = 1
!
!              WORK(IVT) is M by M
!
            ldwkvt = M
            itau = ivt + ldwkvt*M
            nwork = itau + M
!
!              Compute A=L*Q, copying result to VT
!              Workspace: need   M*M [VT] + M [tau] + M    [work]
!              Workspace: prefer M*M [VT] + M [tau] + M*NB [work]
!
            CALL DGELQF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  ierr)
            CALL DLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!              Generate Q in VT
!              Workspace: need   M*M [VT] + M [tau] + N    [work]
!              Workspace: prefer M*M [VT] + M [tau] + N*NB [work]
!
            CALL DORGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(nwork),           &
     &                  Lwork-nwork+1,ierr)
!
!              Produce L in A, zeroing out other entries
!
            CALL DLASET('U',M-1,M-1,ZERO,ZERO,A(1,2),Lda)
            ie = itau
            itauq = ie + M
            itaup = itauq + M
            nwork = itaup + M
!
!              Bidiagonalize L in A
!              Workspace: need   M*M [VT] + 3*M [e, tauq, taup] + M      [work]
!              Workspace: prefer M*M [VT] + 3*M [e, tauq, taup] + 2*M*NB [work]
!
            CALL DGEBRD(M,M,A,Lda,S,Work(ie),Work(itauq),Work(itaup),   &
     &                  Work(nwork),Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in WORK(IVT)
!              Workspace: need   M*M [VT] + 3*M [e, tauq, taup] + BDSPAC
!
            CALL DBDSDC('U','I',M,S,Work(ie),U,Ldu,Work(ivt),ldwkvt,dum,&
     &                  idum,Work(nwork),Iwork,Info)
!
!              Overwrite U by left singular vectors of L and WORK(IVT)
!              by right singular vectors of L
!              Workspace: need   M*M [VT] + 3*M [e, tauq, taup]+ M    [work]
!              Workspace: prefer M*M [VT] + 3*M [e, tauq, taup]+ M*NB [work]
!
            CALL DORMBR('Q','L','N',M,M,M,A,Lda,Work(itauq),U,Ldu,      &
     &                  Work(nwork),Lwork-nwork+1,ierr)
            CALL DORMBR('P','R','T',M,M,M,A,Lda,Work(itaup),Work(ivt),  &
     &                  ldwkvt,Work(nwork),Lwork-nwork+1,ierr)
!
!              Multiply right singular vectors of L in WORK(IVT) by
!              Q in VT, storing result in A
!              Workspace: need   M*M [VT]
!
            CALL DGEMM('N','N',M,N,M,ONE,Work(ivt),ldwkvt,Vt,Ldvt,ZERO, &
     &                 A,Lda)
!
!              Copy right singular vectors of A from A to VT
!
            CALL DLACPY('F',M,N,A,Lda,Vt,Ldvt)
!
         ENDIF
!
      ELSE
!
!           N .LT. MNTHR
!
!           Path 5t (N > M, but not much larger)
!           Reduce to bidiagonal form without LQ decomposition
!
         ie = 1
         itauq = ie + M
         itaup = itauq + M
         nwork = itaup + M
!
!           Bidiagonalize A
!           Workspace: need   3*M [e, tauq, taup] + N        [work]
!           Workspace: prefer 3*M [e, tauq, taup] + (M+N)*NB [work]
!
         CALL DGEBRD(M,N,A,Lda,S,Work(ie),Work(itauq),Work(itaup),      &
     &               Work(nwork),Lwork-nwork+1,ierr)
         IF ( wntqn ) THEN
!
!              Path 5tn (N > M, JOBZ='N')
!              Perform bidiagonal SVD, only computing singular values
!              Workspace: need   3*M [e, tauq, taup] + BDSPAC
!
            CALL DBDSDC('L','N',M,S,Work(ie),dum,1,dum,1,dum,idum,      &
     &                  Work(nwork),Iwork,Info)
         ELSEIF ( wntqo ) THEN
!              Path 5to (N > M, JOBZ='O')
            ldwkvt = M
            ivt = nwork
            IF ( Lwork>=M*N+3*M+bdspac ) THEN
!
!                 WORK( IVT ) is M by N
!
               CALL DLASET('F',M,N,ZERO,ZERO,Work(ivt),ldwkvt)
               nwork = ivt + ldwkvt*N
!                 IL is unused; silence compile warnings
               il = -1
            ELSE
!
!                 WORK( IVT ) is M by M
!
               nwork = ivt + ldwkvt*M
               il = nwork
!
!                 WORK(IL) is M by CHUNK
!
               chunk = (Lwork-M*M-3*M)/M
            ENDIF
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in WORK(IVT)
!              Workspace: need   3*M [e, tauq, taup] + M*M [VT] + BDSPAC
!
            CALL DBDSDC('L','I',M,S,Work(ie),U,Ldu,Work(ivt),ldwkvt,dum,&
     &                  idum,Work(nwork),Iwork,Info)
!
!              Overwrite U by left singular vectors of A
!              Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M    [work]
!              Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*NB [work]
!
            CALL DORMBR('Q','L','N',M,M,N,A,Lda,Work(itauq),U,Ldu,      &
     &                  Work(nwork),Lwork-nwork+1,ierr)
!
            IF ( Lwork>=M*N+3*M+bdspac ) THEN
!
!                 Path 5to-fast
!                 Overwrite WORK(IVT) by left singular vectors of A
!                 Workspace: need   3*M [e, tauq, taup] + M*N [VT] + M    [work]
!                 Workspace: prefer 3*M [e, tauq, taup] + M*N [VT] + M*NB [work]
!
               CALL DORMBR('P','R','T',M,N,M,A,Lda,Work(itaup),Work(ivt)&
     &                     ,ldwkvt,Work(nwork),Lwork-nwork+1,ierr)
!
!                 Copy right singular vectors of A from WORK(IVT) to A
!
               CALL DLACPY('F',M,N,Work(ivt),ldwkvt,A,Lda)
            ELSE
!
!                 Path 5to-slow
!                 Generate P**T in A
!                 Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M    [work]
!                 Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*NB [work]
!
               CALL DORGBR('P',M,N,M,A,Lda,Work(itaup),Work(nwork),     &
     &                     Lwork-nwork+1,ierr)
!
!                 Multiply Q in A by right singular vectors of
!                 bidiagonal matrix in WORK(IVT), storing result in
!                 WORK(IL) and copying to A
!                 Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M*NB [L]
!                 Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*N  [L]
!
               DO i = 1 , N , chunk
                  blk = MIN(N-i+1,chunk)
                  CALL DGEMM('N','N',M,blk,M,ONE,Work(ivt),ldwkvt,A(1,i)&
     &                       ,Lda,ZERO,Work(il),M)
                  CALL DLACPY('F',M,blk,Work(il),M,A(1,i),Lda)
               ENDDO
            ENDIF
         ELSEIF ( wntqs ) THEN
!
!              Path 5ts (N > M, JOBZ='S')
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in VT
!              Workspace: need   3*M [e, tauq, taup] + BDSPAC
!
            CALL DLASET('F',M,N,ZERO,ZERO,Vt,Ldvt)
            CALL DBDSDC('L','I',M,S,Work(ie),U,Ldu,Vt,Ldvt,dum,idum,    &
     &                  Work(nwork),Iwork,Info)
!
!              Overwrite U by left singular vectors of A and VT
!              by right singular vectors of A
!              Workspace: need   3*M [e, tauq, taup] + M    [work]
!              Workspace: prefer 3*M [e, tauq, taup] + M*NB [work]
!
            CALL DORMBR('Q','L','N',M,M,N,A,Lda,Work(itauq),U,Ldu,      &
     &                  Work(nwork),Lwork-nwork+1,ierr)
            CALL DORMBR('P','R','T',M,N,M,A,Lda,Work(itaup),Vt,Ldvt,    &
     &                  Work(nwork),Lwork-nwork+1,ierr)
         ELSEIF ( wntqa ) THEN
!
!              Path 5ta (N > M, JOBZ='A')
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in U and computing right singular
!              vectors of bidiagonal matrix in VT
!              Workspace: need   3*M [e, tauq, taup] + BDSPAC
!
            CALL DLASET('F',N,N,ZERO,ZERO,Vt,Ldvt)
            CALL DBDSDC('L','I',M,S,Work(ie),U,Ldu,Vt,Ldvt,dum,idum,    &
     &                  Work(nwork),Iwork,Info)
!
!              Set the right corner of VT to identity matrix
!
            IF ( N>M ) CALL DLASET('F',N-M,N-M,ZERO,ONE,Vt(M+1,M+1),    &
     &                             Ldvt)
!
!              Overwrite U by left singular vectors of A and VT
!              by right singular vectors of A
!              Workspace: need   3*M [e, tauq, taup] + N    [work]
!              Workspace: prefer 3*M [e, tauq, taup] + N*NB [work]
!
            CALL DORMBR('Q','L','N',M,M,N,A,Lda,Work(itauq),U,Ldu,      &
     &                  Work(nwork),Lwork-nwork+1,ierr)
            CALL DORMBR('P','R','T',N,N,M,A,Lda,Work(itaup),Vt,Ldvt,    &
     &                  Work(nwork),Lwork-nwork+1,ierr)
         ENDIF
!
!
      ENDIF
!
!     Undo scaling if necessary
!
      IF ( iscl==1 ) THEN
         IF ( anrm>bignum ) CALL DLASCL('G',0,0,bignum,anrm,minmn,1,S,  &
     &                                  minmn,ierr)
         IF ( anrm<smlnum ) CALL DLASCL('G',0,0,smlnum,anrm,minmn,1,S,  &
     &                                  minmn,ierr)
      ENDIF
!
!     Return optimal workspace in WORK(1)
!
      Work(1) = maxwrk
!
!
!     End of DGESDD
!
      END SUBROUTINE DGESDD
