!*==cgesdd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CGESDD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGESDD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesdd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesdd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesdd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT,
!                          WORK, LWORK, RWORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ
!       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               RWORK( * ), S( * )
!       COMPLEX            A( LDA, * ), U( LDU, * ), VT( LDVT, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGESDD computes the singular value decomposition (SVD) of a complex
!> M-by-N matrix A, optionally computing the left and/or right singular
!> vectors, by using divide-and-conquer method. The SVD is written
!>
!>      A = U * SIGMA * conjugate-transpose(V)
!>
!> where SIGMA is an M-by-N matrix which is zero except for its
!> min(m,n) diagonal elements, U is an M-by-M unitary matrix, and
!> V is an N-by-N unitary matrix.  The diagonal elements of SIGMA
!> are the singular values of A; they are real and non-negative, and
!> are returned in descending order.  The first min(m,n) columns of
!> U and V are the left and right singular vectors of A.
!>
!> Note that the routine returns VT = V**H, not V.
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
!>          = 'A':  all M columns of U and all N rows of V**H are
!>                  returned in the arrays U and VT;
!>          = 'S':  the first min(M,N) columns of U and the first
!>                  min(M,N) rows of V**H are returned in the arrays U
!>                  and VT;
!>          = 'O':  If M >= N, the first N columns of U are overwritten
!>                  in the array A and all rows of V**H are returned in
!>                  the array VT;
!>                  otherwise, all columns of U are returned in the
!>                  array U and the first M rows of V**H are overwritten
!>                  in the array A;
!>          = 'N':  no columns of U or rows of V**H are computed.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit,
!>          if JOBZ = 'O',  A is overwritten with the first N columns
!>                          of U (the left singular vectors, stored
!>                          columnwise) if M >= N;
!>                          A is overwritten with the first M rows
!>                          of V**H (the right singular vectors, stored
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
!>          S is REAL array, dimension (min(M,N))
!>          The singular values of A, sorted so that S(i) >= S(i+1).
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX array, dimension (LDU,UCOL)
!>          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
!>          UCOL = min(M,N) if JOBZ = 'S'.
!>          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
!>          unitary matrix U;
!>          if JOBZ = 'S', U contains the first min(M,N) columns of U
!>          (the left singular vectors, stored columnwise);
!>          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= 1;
!>          if JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is COMPLEX array, dimension (LDVT,N)
!>          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
!>          N-by-N unitary matrix V**H;
!>          if JOBZ = 'S', VT contains the first min(M,N) rows of
!>          V**H (the right singular vectors, stored rowwise);
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
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
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
!>          If JOBZ = 'N', LWORK >= 2*mn + mx.
!>          If JOBZ = 'O', LWORK >= 2*mn*mn + 2*mn + mx.
!>          If JOBZ = 'S', LWORK >=   mn*mn + 3*mn.
!>          If JOBZ = 'A', LWORK >=   mn*mn + 2*mn + mx.
!>          These are not tight minimums in all cases; see comments inside code.
!>          For good performance, LWORK should generally be larger;
!>          a query is recommended.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (MAX(1,LRWORK))
!>          Let mx = max(M,N) and mn = min(M,N).
!>          If JOBZ = 'N',    LRWORK >= 5*mn (LAPACK <= 3.6 needs 7*mn);
!>          else if mx >> mn, LRWORK >= 5*mn*mn + 5*mn;
!>          else              LRWORK >= max( 5*mn*mn + 5*mn,
!>                                           2*mx*mn + 2*mn*mn + mn ).
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
!>          > 0:  The updating process of SBDSDC did not converge.
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
!> \ingroup complexGEsing
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Huan Ren, Computer Science Division, University of
!>     California at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE CGESDD(Jobz,M,N,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,Rwork,&
     &                  Iwork,Info)
      IMPLICIT NONE
!*--CGESDD230
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobz
      INTEGER Info , Lda , Ldu , Ldvt , Lwork , M , N
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      REAL Rwork(*) , S(*)
      COMPLEX A(Lda,*) , U(Ldu,*) , Vt(Ldvt,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , wntqa , wntqas , wntqn , wntqo , wntqs
      INTEGER blk , chunk , i , ie , ierr , il , ir , iru , irvt ,      &
     &        iscl , itau , itaup , itauq , iu , ivt , ldwkvt , ldwrkl ,&
     &        ldwrkr , ldwrku , maxwrk , minmn , minwrk , mnthr1 ,      &
     &        mnthr2 , nrwork , nwork , wrkbl
      INTEGER lwork_cgebrd_mn , lwork_cgebrd_mm , lwork_cgebrd_nn ,     &
     &        lwork_cgelqf_mn , lwork_cgeqrf_mn , lwork_cungbr_p_mn ,   &
     &        lwork_cungbr_p_nn , lwork_cungbr_q_mn ,                   &
     &        lwork_cungbr_q_mm , lwork_cunglq_mn , lwork_cunglq_nn ,   &
     &        lwork_cungqr_mm , lwork_cungqr_mn , lwork_cunmbr_prc_mm , &
     &        lwork_cunmbr_qln_mm , lwork_cunmbr_prc_mn ,               &
     &        lwork_cunmbr_qln_mn , lwork_cunmbr_prc_nn ,               &
     &        lwork_cunmbr_qln_nn
      REAL anrm , bignum , eps , smlnum
!     ..
!     .. Local Arrays ..
      INTEGER idum(1)
      REAL dum(1)
      COMPLEX cdum(1)
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEBRD , CGELQF , CGEMM , CGEQRF , CLACP2 , CLACPY ,     &
     &         CLACRM , CLARCM , CLASCL , CLASET , CUNGBR , CUNGLQ ,    &
     &         CUNGQR , CUNMBR , SBDSDC , SLASCL , XERBLA
!     ..
!     .. External Functions ..
      LOGICAL LSAME , SISNAN
      REAL SLAMCH , CLANGE
      EXTERNAL LSAME , SLAMCH , CLANGE , SISNAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC INT , MAX , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      minmn = MIN(M,N)
      mnthr1 = INT(minmn*17.0E0/9.0E0)
      mnthr2 = INT(minmn*5.0E0/3.0E0)
      wntqa = LSAME(Jobz,'A')
      wntqs = LSAME(Jobz,'S')
      wntqas = wntqa .OR. wntqs
      wntqo = LSAME(Jobz,'O')
      wntqn = LSAME(Jobz,'N')
      lquery = (Lwork==-1)
      minwrk = 1
      maxwrk = 1
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
!       CWorkspace refers to complex workspace, and RWorkspace to
!       real workspace. NB refers to the optimal block size for the
!       immediately following subroutine, as returned by ILAENV.)
!
      IF ( Info==0 ) THEN
         minwrk = 1
         maxwrk = 1
         IF ( M>=N .AND. minmn>0 ) THEN
!
!           There is no complex work space needed for bidiagonal SVD
!           The real work space needed for bidiagonal SVD (sbdsdc) is
!           BDSPAC = 3*N*N + 4*N for singular values and vectors;
!           BDSPAC = 4*N         for singular values only;
!           not including e, RU, and RVT matrices.
!
!           Compute space preferred for each routine
            CALL CGEBRD(M,N,cdum(1),M,dum(1),dum(1),cdum(1),cdum(1),    &
     &                  cdum(1),-1,ierr)
            lwork_cgebrd_mn = INT(cdum(1))
!
            CALL CGEBRD(N,N,cdum(1),N,dum(1),dum(1),cdum(1),cdum(1),    &
     &                  cdum(1),-1,ierr)
            lwork_cgebrd_nn = INT(cdum(1))
!
            CALL CGEQRF(M,N,cdum(1),M,cdum(1),cdum(1),-1,ierr)
            lwork_cgeqrf_mn = INT(cdum(1))
!
            CALL CUNGBR('P',N,N,N,cdum(1),N,cdum(1),cdum(1),-1,ierr)
            lwork_cungbr_p_nn = INT(cdum(1))
!
            CALL CUNGBR('Q',M,M,N,cdum(1),M,cdum(1),cdum(1),-1,ierr)
            lwork_cungbr_q_mm = INT(cdum(1))
!
            CALL CUNGBR('Q',M,N,N,cdum(1),M,cdum(1),cdum(1),-1,ierr)
            lwork_cungbr_q_mn = INT(cdum(1))
!
            CALL CUNGQR(M,M,N,cdum(1),M,cdum(1),cdum(1),-1,ierr)
            lwork_cungqr_mm = INT(cdum(1))
!
            CALL CUNGQR(M,N,N,cdum(1),M,cdum(1),cdum(1),-1,ierr)
            lwork_cungqr_mn = INT(cdum(1))
!
            CALL CUNMBR('P','R','C',N,N,N,cdum(1),N,cdum(1),cdum(1),N,  &
     &                  cdum(1),-1,ierr)
            lwork_cunmbr_prc_nn = INT(cdum(1))
!
            CALL CUNMBR('Q','L','N',M,M,N,cdum(1),M,cdum(1),cdum(1),M,  &
     &                  cdum(1),-1,ierr)
            lwork_cunmbr_qln_mm = INT(cdum(1))
!
            CALL CUNMBR('Q','L','N',M,N,N,cdum(1),M,cdum(1),cdum(1),M,  &
     &                  cdum(1),-1,ierr)
            lwork_cunmbr_qln_mn = INT(cdum(1))
!
            CALL CUNMBR('Q','L','N',N,N,N,cdum(1),N,cdum(1),cdum(1),N,  &
     &                  cdum(1),-1,ierr)
            lwork_cunmbr_qln_nn = INT(cdum(1))
!
            IF ( M>=mnthr1 ) THEN
               IF ( wntqn ) THEN
!
!                 Path 1 (M >> N, JOBZ='N')
!
                  maxwrk = N + lwork_cgeqrf_mn
                  maxwrk = MAX(maxwrk,2*N+lwork_cgebrd_nn)
                  minwrk = 3*N
               ELSEIF ( wntqo ) THEN
!
!                 Path 2 (M >> N, JOBZ='O')
!
                  wrkbl = N + lwork_cgeqrf_mn
                  wrkbl = MAX(wrkbl,N+lwork_cungqr_mn)
                  wrkbl = MAX(wrkbl,2*N+lwork_cgebrd_nn)
                  wrkbl = MAX(wrkbl,2*N+lwork_cunmbr_qln_nn)
                  wrkbl = MAX(wrkbl,2*N+lwork_cunmbr_prc_nn)
                  maxwrk = M*N + N*N + wrkbl
                  minwrk = 2*N*N + 3*N
               ELSEIF ( wntqs ) THEN
!
!                 Path 3 (M >> N, JOBZ='S')
!
                  wrkbl = N + lwork_cgeqrf_mn
                  wrkbl = MAX(wrkbl,N+lwork_cungqr_mn)
                  wrkbl = MAX(wrkbl,2*N+lwork_cgebrd_nn)
                  wrkbl = MAX(wrkbl,2*N+lwork_cunmbr_qln_nn)
                  wrkbl = MAX(wrkbl,2*N+lwork_cunmbr_prc_nn)
                  maxwrk = N*N + wrkbl
                  minwrk = N*N + 3*N
               ELSEIF ( wntqa ) THEN
!
!                 Path 4 (M >> N, JOBZ='A')
!
                  wrkbl = N + lwork_cgeqrf_mn
                  wrkbl = MAX(wrkbl,N+lwork_cungqr_mm)
                  wrkbl = MAX(wrkbl,2*N+lwork_cgebrd_nn)
                  wrkbl = MAX(wrkbl,2*N+lwork_cunmbr_qln_nn)
                  wrkbl = MAX(wrkbl,2*N+lwork_cunmbr_prc_nn)
                  maxwrk = N*N + wrkbl
                  minwrk = N*N + MAX(3*N,N+M)
               ENDIF
            ELSEIF ( M>=mnthr2 ) THEN
!
!              Path 5 (M >> N, but not as much as MNTHR1)
!
               maxwrk = 2*N + lwork_cgebrd_mn
               minwrk = 2*N + M
               IF ( wntqo ) THEN
!                 Path 5o (M >> N, JOBZ='O')
                  maxwrk = MAX(maxwrk,2*N+lwork_cungbr_p_nn)
                  maxwrk = MAX(maxwrk,2*N+lwork_cungbr_q_mn)
                  maxwrk = maxwrk + M*N
                  minwrk = minwrk + N*N
               ELSEIF ( wntqs ) THEN
!                 Path 5s (M >> N, JOBZ='S')
                  maxwrk = MAX(maxwrk,2*N+lwork_cungbr_p_nn)
                  maxwrk = MAX(maxwrk,2*N+lwork_cungbr_q_mn)
               ELSEIF ( wntqa ) THEN
!                 Path 5a (M >> N, JOBZ='A')
                  maxwrk = MAX(maxwrk,2*N+lwork_cungbr_p_nn)
                  maxwrk = MAX(maxwrk,2*N+lwork_cungbr_q_mm)
               ENDIF
            ELSE
!
!              Path 6 (M >= N, but not much larger)
!
               maxwrk = 2*N + lwork_cgebrd_mn
               minwrk = 2*N + M
               IF ( wntqo ) THEN
!                 Path 6o (M >= N, JOBZ='O')
                  maxwrk = MAX(maxwrk,2*N+lwork_cunmbr_prc_nn)
                  maxwrk = MAX(maxwrk,2*N+lwork_cunmbr_qln_mn)
                  maxwrk = maxwrk + M*N
                  minwrk = minwrk + N*N
               ELSEIF ( wntqs ) THEN
!                 Path 6s (M >= N, JOBZ='S')
                  maxwrk = MAX(maxwrk,2*N+lwork_cunmbr_qln_mn)
                  maxwrk = MAX(maxwrk,2*N+lwork_cunmbr_prc_nn)
               ELSEIF ( wntqa ) THEN
!                 Path 6a (M >= N, JOBZ='A')
                  maxwrk = MAX(maxwrk,2*N+lwork_cunmbr_qln_mm)
                  maxwrk = MAX(maxwrk,2*N+lwork_cunmbr_prc_nn)
               ENDIF
            ENDIF
         ELSEIF ( minmn>0 ) THEN
!
!           There is no complex work space needed for bidiagonal SVD
!           The real work space needed for bidiagonal SVD (sbdsdc) is
!           BDSPAC = 3*M*M + 4*M for singular values and vectors;
!           BDSPAC = 4*M         for singular values only;
!           not including e, RU, and RVT matrices.
!
!           Compute space preferred for each routine
            CALL CGEBRD(M,N,cdum(1),M,dum(1),dum(1),cdum(1),cdum(1),    &
     &                  cdum(1),-1,ierr)
            lwork_cgebrd_mn = INT(cdum(1))
!
            CALL CGEBRD(M,M,cdum(1),M,dum(1),dum(1),cdum(1),cdum(1),    &
     &                  cdum(1),-1,ierr)
            lwork_cgebrd_mm = INT(cdum(1))
!
            CALL CGELQF(M,N,cdum(1),M,cdum(1),cdum(1),-1,ierr)
            lwork_cgelqf_mn = INT(cdum(1))
!
            CALL CUNGBR('P',M,N,M,cdum(1),M,cdum(1),cdum(1),-1,ierr)
            lwork_cungbr_p_mn = INT(cdum(1))
!
            CALL CUNGBR('P',N,N,M,cdum(1),N,cdum(1),cdum(1),-1,ierr)
            lwork_cungbr_p_nn = INT(cdum(1))
!
            CALL CUNGBR('Q',M,M,N,cdum(1),M,cdum(1),cdum(1),-1,ierr)
            lwork_cungbr_q_mm = INT(cdum(1))
!
            CALL CUNGLQ(M,N,M,cdum(1),M,cdum(1),cdum(1),-1,ierr)
            lwork_cunglq_mn = INT(cdum(1))
!
            CALL CUNGLQ(N,N,M,cdum(1),N,cdum(1),cdum(1),-1,ierr)
            lwork_cunglq_nn = INT(cdum(1))
!
            CALL CUNMBR('P','R','C',M,M,M,cdum(1),M,cdum(1),cdum(1),M,  &
     &                  cdum(1),-1,ierr)
            lwork_cunmbr_prc_mm = INT(cdum(1))
!
            CALL CUNMBR('P','R','C',M,N,M,cdum(1),M,cdum(1),cdum(1),M,  &
     &                  cdum(1),-1,ierr)
            lwork_cunmbr_prc_mn = INT(cdum(1))
!
            CALL CUNMBR('P','R','C',N,N,M,cdum(1),N,cdum(1),cdum(1),N,  &
     &                  cdum(1),-1,ierr)
            lwork_cunmbr_prc_nn = INT(cdum(1))
!
            CALL CUNMBR('Q','L','N',M,M,M,cdum(1),M,cdum(1),cdum(1),M,  &
     &                  cdum(1),-1,ierr)
            lwork_cunmbr_qln_mm = INT(cdum(1))
!
            IF ( N>=mnthr1 ) THEN
               IF ( wntqn ) THEN
!
!                 Path 1t (N >> M, JOBZ='N')
!
                  maxwrk = M + lwork_cgelqf_mn
                  maxwrk = MAX(maxwrk,2*M+lwork_cgebrd_mm)
                  minwrk = 3*M
               ELSEIF ( wntqo ) THEN
!
!                 Path 2t (N >> M, JOBZ='O')
!
                  wrkbl = M + lwork_cgelqf_mn
                  wrkbl = MAX(wrkbl,M+lwork_cunglq_mn)
                  wrkbl = MAX(wrkbl,2*M+lwork_cgebrd_mm)
                  wrkbl = MAX(wrkbl,2*M+lwork_cunmbr_qln_mm)
                  wrkbl = MAX(wrkbl,2*M+lwork_cunmbr_prc_mm)
                  maxwrk = M*N + M*M + wrkbl
                  minwrk = 2*M*M + 3*M
               ELSEIF ( wntqs ) THEN
!
!                 Path 3t (N >> M, JOBZ='S')
!
                  wrkbl = M + lwork_cgelqf_mn
                  wrkbl = MAX(wrkbl,M+lwork_cunglq_mn)
                  wrkbl = MAX(wrkbl,2*M+lwork_cgebrd_mm)
                  wrkbl = MAX(wrkbl,2*M+lwork_cunmbr_qln_mm)
                  wrkbl = MAX(wrkbl,2*M+lwork_cunmbr_prc_mm)
                  maxwrk = M*M + wrkbl
                  minwrk = M*M + 3*M
               ELSEIF ( wntqa ) THEN
!
!                 Path 4t (N >> M, JOBZ='A')
!
                  wrkbl = M + lwork_cgelqf_mn
                  wrkbl = MAX(wrkbl,M+lwork_cunglq_nn)
                  wrkbl = MAX(wrkbl,2*M+lwork_cgebrd_mm)
                  wrkbl = MAX(wrkbl,2*M+lwork_cunmbr_qln_mm)
                  wrkbl = MAX(wrkbl,2*M+lwork_cunmbr_prc_mm)
                  maxwrk = M*M + wrkbl
                  minwrk = M*M + MAX(3*M,M+N)
               ENDIF
            ELSEIF ( N>=mnthr2 ) THEN
!
!              Path 5t (N >> M, but not as much as MNTHR1)
!
               maxwrk = 2*M + lwork_cgebrd_mn
               minwrk = 2*M + N
               IF ( wntqo ) THEN
!                 Path 5to (N >> M, JOBZ='O')
                  maxwrk = MAX(maxwrk,2*M+lwork_cungbr_q_mm)
                  maxwrk = MAX(maxwrk,2*M+lwork_cungbr_p_mn)
                  maxwrk = maxwrk + M*N
                  minwrk = minwrk + M*M
               ELSEIF ( wntqs ) THEN
!                 Path 5ts (N >> M, JOBZ='S')
                  maxwrk = MAX(maxwrk,2*M+lwork_cungbr_q_mm)
                  maxwrk = MAX(maxwrk,2*M+lwork_cungbr_p_mn)
               ELSEIF ( wntqa ) THEN
!                 Path 5ta (N >> M, JOBZ='A')
                  maxwrk = MAX(maxwrk,2*M+lwork_cungbr_q_mm)
                  maxwrk = MAX(maxwrk,2*M+lwork_cungbr_p_nn)
               ENDIF
            ELSE
!
!              Path 6t (N > M, but not much larger)
!
               maxwrk = 2*M + lwork_cgebrd_mn
               minwrk = 2*M + N
               IF ( wntqo ) THEN
!                 Path 6to (N > M, JOBZ='O')
                  maxwrk = MAX(maxwrk,2*M+lwork_cunmbr_qln_mm)
                  maxwrk = MAX(maxwrk,2*M+lwork_cunmbr_prc_mn)
                  maxwrk = maxwrk + M*N
                  minwrk = minwrk + M*M
               ELSEIF ( wntqs ) THEN
!                 Path 6ts (N > M, JOBZ='S')
                  maxwrk = MAX(maxwrk,2*M+lwork_cunmbr_qln_mm)
                  maxwrk = MAX(maxwrk,2*M+lwork_cunmbr_prc_mn)
               ELSEIF ( wntqa ) THEN
!                 Path 6ta (N > M, JOBZ='A')
                  maxwrk = MAX(maxwrk,2*M+lwork_cunmbr_qln_mm)
                  maxwrk = MAX(maxwrk,2*M+lwork_cunmbr_prc_nn)
               ENDIF
            ENDIF
         ENDIF
         maxwrk = MAX(maxwrk,minwrk)
      ENDIF
      IF ( Info==0 ) THEN
         Work(1) = maxwrk
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -12
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGESDD',-Info)
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
      eps = SLAMCH('P')
      smlnum = SQRT(SLAMCH('S'))/eps
      bignum = ONE/smlnum
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      anrm = CLANGE('M',M,N,A,Lda,dum)
      IF ( SISNAN(anrm) ) THEN
         Info = -4
         RETURN
      ENDIF
      iscl = 0
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         iscl = 1
         CALL CLASCL('G',0,0,anrm,smlnum,M,N,A,Lda,ierr)
      ELSEIF ( anrm>bignum ) THEN
         iscl = 1
         CALL CLASCL('G',0,0,anrm,bignum,M,N,A,Lda,ierr)
      ENDIF
!
      IF ( M>=N ) THEN
!
!        A has at least as many rows as columns. If A has sufficiently
!        more rows than columns, first reduce using the QR
!        decomposition (if sufficient workspace available)
!
         IF ( M>=mnthr1 ) THEN
!
            IF ( wntqn ) THEN
!
!              Path 1 (M >> N, JOBZ='N')
!              No singular vectors to be computed
!
               itau = 1
               nwork = itau + N
!
!              Compute A=Q*R
!              CWorkspace: need   N [tau] + N    [work]
!              CWorkspace: prefer N [tau] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CGEQRF(M,N,A,Lda,Work(itau),Work(nwork),            &
     &                     Lwork-nwork+1,ierr)
!
!              Zero out below R
!
               CALL CLASET('L',N-1,N-1,CZERO,CZERO,A(2,1),Lda)
               ie = 1
               itauq = 1
               itaup = itauq + N
               nwork = itaup + N
!
!              Bidiagonalize R in A
!              CWorkspace: need   2*N [tauq, taup] + N      [work]
!              CWorkspace: prefer 2*N [tauq, taup] + 2*N*NB [work]
!              RWorkspace: need   N [e]
!
               CALL CGEBRD(N,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup)&
     &                     ,Work(nwork),Lwork-nwork+1,ierr)
               nrwork = ie + N
!
!              Perform bidiagonal SVD, compute singular values only
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + BDSPAC
!
               CALL SBDSDC('U','N',N,S,Rwork(ie),dum,1,dum,1,dum,idum,  &
     &                     Rwork(nrwork),Iwork,Info)
!
            ELSEIF ( wntqo ) THEN
!
!              Path 2 (M >> N, JOBZ='O')
!              N left singular vectors to be overwritten on A and
!              N right singular vectors to be computed in VT
!
               iu = 1
!
!              WORK(IU) is N by N
!
               ldwrku = N
               ir = iu + ldwrku*N
               IF ( Lwork>=M*N+N*N+3*N ) THEN
!
!                 WORK(IR) is M by N
!
                  ldwrkr = M
               ELSE
                  ldwrkr = (Lwork-N*N-3*N)/N
               ENDIF
               itau = ir + ldwrkr*N
               nwork = itau + N
!
!              Compute A=Q*R
!              CWorkspace: need   N*N [U] + N*N [R] + N [tau] + N    [work]
!              CWorkspace: prefer N*N [U] + N*N [R] + N [tau] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CGEQRF(M,N,A,Lda,Work(itau),Work(nwork),            &
     &                     Lwork-nwork+1,ierr)
!
!              Copy R to WORK( IR ), zeroing out below it
!
               CALL CLACPY('U',N,N,A,Lda,Work(ir),ldwrkr)
               CALL CLASET('L',N-1,N-1,CZERO,CZERO,Work(ir+1),ldwrkr)
!
!              Generate Q in A
!              CWorkspace: need   N*N [U] + N*N [R] + N [tau] + N    [work]
!              CWorkspace: prefer N*N [U] + N*N [R] + N [tau] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CUNGQR(M,N,N,A,Lda,Work(itau),Work(nwork),          &
     &                     Lwork-nwork+1,ierr)
               ie = 1
               itauq = itau
               itaup = itauq + N
               nwork = itaup + N
!
!              Bidiagonalize R in WORK(IR)
!              CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N      [work]
!              CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + 2*N*NB [work]
!              RWorkspace: need   N [e]
!
               CALL CGEBRD(N,N,Work(ir),ldwrkr,S,Rwork(ie),Work(itauq), &
     &                     Work(itaup),Work(nwork),Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of R in WORK(IRU) and computing right singular vectors
!              of R in WORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
!
               iru = ie + N
               irvt = iru + N*N
               nrwork = irvt + N*N
               CALL SBDSDC('U','I',N,S,Rwork(ie),Rwork(iru),N,          &
     &                     Rwork(irvt),N,dum,idum,Rwork(nrwork),Iwork,  &
     &                     Info)
!
!              Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
!              Overwrite WORK(IU) by the left singular vectors of R
!              CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CLACP2('F',N,N,Rwork(iru),N,Work(iu),ldwrku)
               CALL CUNMBR('Q','L','N',N,N,N,Work(ir),ldwrkr,Work(itauq)&
     &                     ,Work(iu),ldwrku,Work(nwork),Lwork-nwork+1,  &
     &                     ierr)
!
!              Copy real matrix RWORK(IRVT) to complex matrix VT
!              Overwrite VT by the right singular vectors of R
!              CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CLACP2('F',N,N,Rwork(irvt),N,Vt,Ldvt)
               CALL CUNMBR('P','R','C',N,N,N,Work(ir),ldwrkr,Work(itaup)&
     &                     ,Vt,Ldvt,Work(nwork),Lwork-nwork+1,ierr)
!
!              Multiply Q in A by left singular vectors of R in
!              WORK(IU), storing result in WORK(IR) and copying to A
!              CWorkspace: need   N*N [U] + N*N [R]
!              CWorkspace: prefer N*N [U] + M*N [R]
!              RWorkspace: need   0
!
               DO i = 1 , M , ldwrkr
                  chunk = MIN(M-i+1,ldwrkr)
                  CALL CGEMM('N','N',chunk,N,N,CONE,A(i,1),Lda,Work(iu),&
     &                       ldwrku,CZERO,Work(ir),ldwrkr)
                  CALL CLACPY('F',chunk,N,Work(ir),ldwrkr,A(i,1),Lda)
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
!              CWorkspace: need   N*N [R] + N [tau] + N    [work]
!              CWorkspace: prefer N*N [R] + N [tau] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CGEQRF(M,N,A,Lda,Work(itau),Work(nwork),            &
     &                     Lwork-nwork+1,ierr)
!
!              Copy R to WORK(IR), zeroing out below it
!
               CALL CLACPY('U',N,N,A,Lda,Work(ir),ldwrkr)
               CALL CLASET('L',N-1,N-1,CZERO,CZERO,Work(ir+1),ldwrkr)
!
!              Generate Q in A
!              CWorkspace: need   N*N [R] + N [tau] + N    [work]
!              CWorkspace: prefer N*N [R] + N [tau] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CUNGQR(M,N,N,A,Lda,Work(itau),Work(nwork),          &
     &                     Lwork-nwork+1,ierr)
               ie = 1
               itauq = itau
               itaup = itauq + N
               nwork = itaup + N
!
!              Bidiagonalize R in WORK(IR)
!              CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N      [work]
!              CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + 2*N*NB [work]
!              RWorkspace: need   N [e]
!
               CALL CGEBRD(N,N,Work(ir),ldwrkr,S,Rwork(ie),Work(itauq), &
     &                     Work(itaup),Work(nwork),Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
!
               iru = ie + N
               irvt = iru + N*N
               nrwork = irvt + N*N
               CALL SBDSDC('U','I',N,S,Rwork(ie),Rwork(iru),N,          &
     &                     Rwork(irvt),N,dum,idum,Rwork(nrwork),Iwork,  &
     &                     Info)
!
!              Copy real matrix RWORK(IRU) to complex matrix U
!              Overwrite U by left singular vectors of R
!              CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CLACP2('F',N,N,Rwork(iru),N,U,Ldu)
               CALL CUNMBR('Q','L','N',N,N,N,Work(ir),ldwrkr,Work(itauq)&
     &                     ,U,Ldu,Work(nwork),Lwork-nwork+1,ierr)
!
!              Copy real matrix RWORK(IRVT) to complex matrix VT
!              Overwrite VT by right singular vectors of R
!              CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CLACP2('F',N,N,Rwork(irvt),N,Vt,Ldvt)
               CALL CUNMBR('P','R','C',N,N,N,Work(ir),ldwrkr,Work(itaup)&
     &                     ,Vt,Ldvt,Work(nwork),Lwork-nwork+1,ierr)
!
!              Multiply Q in A by left singular vectors of R in
!              WORK(IR), storing result in U
!              CWorkspace: need   N*N [R]
!              RWorkspace: need   0
!
               CALL CLACPY('F',N,N,U,Ldu,Work(ir),ldwrkr)
               CALL CGEMM('N','N',M,N,N,CONE,A,Lda,Work(ir),ldwrkr,     &
     &                    CZERO,U,Ldu)
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
!              CWorkspace: need   N*N [U] + N [tau] + N    [work]
!              CWorkspace: prefer N*N [U] + N [tau] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CGEQRF(M,N,A,Lda,Work(itau),Work(nwork),            &
     &                     Lwork-nwork+1,ierr)
               CALL CLACPY('L',M,N,A,Lda,U,Ldu)
!
!              Generate Q in U
!              CWorkspace: need   N*N [U] + N [tau] + M    [work]
!              CWorkspace: prefer N*N [U] + N [tau] + M*NB [work]
!              RWorkspace: need   0
!
               CALL CUNGQR(M,M,N,U,Ldu,Work(itau),Work(nwork),          &
     &                     Lwork-nwork+1,ierr)
!
!              Produce R in A, zeroing out below it
!
               CALL CLASET('L',N-1,N-1,CZERO,CZERO,A(2,1),Lda)
               ie = 1
               itauq = itau
               itaup = itauq + N
               nwork = itaup + N
!
!              Bidiagonalize R in A
!              CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N      [work]
!              CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + 2*N*NB [work]
!              RWorkspace: need   N [e]
!
               CALL CGEBRD(N,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup)&
     &                     ,Work(nwork),Lwork-nwork+1,ierr)
               iru = ie + N
               irvt = iru + N*N
               nrwork = irvt + N*N
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
!
               CALL SBDSDC('U','I',N,S,Rwork(ie),Rwork(iru),N,          &
     &                     Rwork(irvt),N,dum,idum,Rwork(nrwork),Iwork,  &
     &                     Info)
!
!              Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
!              Overwrite WORK(IU) by left singular vectors of R
!              CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CLACP2('F',N,N,Rwork(iru),N,Work(iu),ldwrku)
               CALL CUNMBR('Q','L','N',N,N,N,A,Lda,Work(itauq),Work(iu),&
     &                     ldwrku,Work(nwork),Lwork-nwork+1,ierr)
!
!              Copy real matrix RWORK(IRVT) to complex matrix VT
!              Overwrite VT by right singular vectors of R
!              CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CLACP2('F',N,N,Rwork(irvt),N,Vt,Ldvt)
               CALL CUNMBR('P','R','C',N,N,N,A,Lda,Work(itaup),Vt,Ldvt, &
     &                     Work(nwork),Lwork-nwork+1,ierr)
!
!              Multiply Q in U by left singular vectors of R in
!              WORK(IU), storing result in A
!              CWorkspace: need   N*N [U]
!              RWorkspace: need   0
!
               CALL CGEMM('N','N',M,N,N,CONE,U,Ldu,Work(iu),ldwrku,     &
     &                    CZERO,A,Lda)
!
!              Copy left singular vectors of A from A to U
!
               CALL CLACPY('F',M,N,A,Lda,U,Ldu)
!
            ENDIF
!
         ELSEIF ( M>=mnthr2 ) THEN
!
!           MNTHR2 <= M < MNTHR1
!
!           Path 5 (M >> N, but not as much as MNTHR1)
!           Reduce to bidiagonal form without QR decomposition, use
!           CUNGBR and matrix multiplication to compute singular vectors
!
            ie = 1
            nrwork = ie + N
            itauq = 1
            itaup = itauq + N
            nwork = itaup + N
!
!           Bidiagonalize A
!           CWorkspace: need   2*N [tauq, taup] + M        [work]
!           CWorkspace: prefer 2*N [tauq, taup] + (M+N)*NB [work]
!           RWorkspace: need   N [e]
!
            CALL CGEBRD(M,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),  &
     &                  Work(nwork),Lwork-nwork+1,ierr)
            IF ( wntqn ) THEN
!
!              Path 5n (M >> N, JOBZ='N')
!              Compute singular values only
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + BDSPAC
!
               CALL SBDSDC('U','N',N,S,Rwork(ie),dum,1,dum,1,dum,idum,  &
     &                     Rwork(nrwork),Iwork,Info)
            ELSEIF ( wntqo ) THEN
               iu = nwork
               iru = nrwork
               irvt = iru + N*N
               nrwork = irvt + N*N
!
!              Path 5o (M >> N, JOBZ='O')
!              Copy A to VT, generate P**H
!              CWorkspace: need   2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CLACPY('U',N,N,A,Lda,Vt,Ldvt)
               CALL CUNGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(nwork),   &
     &                     Lwork-nwork+1,ierr)
!
!              Generate Q in A
!              CWorkspace: need   2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CUNGBR('Q',M,N,N,A,Lda,Work(itauq),Work(nwork),     &
     &                     Lwork-nwork+1,ierr)
!
               IF ( Lwork>=M*N+3*N ) THEN
!
!                 WORK( IU ) is M by N
!
                  ldwrku = M
               ELSE
!
!                 WORK(IU) is LDWRKU by N
!
                  ldwrku = (Lwork-3*N)/N
               ENDIF
               nwork = iu + ldwrku*N
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
!
               CALL SBDSDC('U','I',N,S,Rwork(ie),Rwork(iru),N,          &
     &                     Rwork(irvt),N,dum,idum,Rwork(nrwork),Iwork,  &
     &                     Info)
!
!              Multiply real matrix RWORK(IRVT) by P**H in VT,
!              storing the result in WORK(IU), copying to VT
!              CWorkspace: need   2*N [tauq, taup] + N*N [U]
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork]
!
               CALL CLARCM(N,N,Rwork(irvt),N,Vt,Ldvt,Work(iu),ldwrku,   &
     &                     Rwork(nrwork))
               CALL CLACPY('F',N,N,Work(iu),ldwrku,Vt,Ldvt)
!
!              Multiply Q in A by real matrix RWORK(IRU), storing the
!              result in WORK(IU), copying to A
!              CWorkspace: need   2*N [tauq, taup] + N*N [U]
!              CWorkspace: prefer 2*N [tauq, taup] + M*N [U]
!              RWorkspace: need   N [e] + N*N [RU] + 2*N*N [rwork]
!              RWorkspace: prefer N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here
!
               nrwork = irvt
               DO i = 1 , M , ldwrku
                  chunk = MIN(M-i+1,ldwrku)
                  CALL CLACRM(chunk,N,A(i,1),Lda,Rwork(iru),N,Work(iu), &
     &                        ldwrku,Rwork(nrwork))
                  CALL CLACPY('F',chunk,N,Work(iu),ldwrku,A(i,1),Lda)
               ENDDO
!
            ELSEIF ( wntqs ) THEN
!
!              Path 5s (M >> N, JOBZ='S')
!              Copy A to VT, generate P**H
!              CWorkspace: need   2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CLACPY('U',N,N,A,Lda,Vt,Ldvt)
               CALL CUNGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(nwork),   &
     &                     Lwork-nwork+1,ierr)
!
!              Copy A to U, generate Q
!              CWorkspace: need   2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CLACPY('L',M,N,A,Lda,U,Ldu)
               CALL CUNGBR('Q',M,N,N,U,Ldu,Work(itauq),Work(nwork),     &
     &                     Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
!
               iru = nrwork
               irvt = iru + N*N
               nrwork = irvt + N*N
               CALL SBDSDC('U','I',N,S,Rwork(ie),Rwork(iru),N,          &
     &                     Rwork(irvt),N,dum,idum,Rwork(nrwork),Iwork,  &
     &                     Info)
!
!              Multiply real matrix RWORK(IRVT) by P**H in VT,
!              storing the result in A, copying to VT
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork]
!
               CALL CLARCM(N,N,Rwork(irvt),N,Vt,Ldvt,A,Lda,Rwork(nrwork)&
     &                     )
               CALL CLACPY('F',N,N,A,Lda,Vt,Ldvt)
!
!              Multiply Q in U by real matrix RWORK(IRU), storing the
!              result in A, copying to U
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here
!
               nrwork = irvt
               CALL CLACRM(M,N,U,Ldu,Rwork(iru),N,A,Lda,Rwork(nrwork))
               CALL CLACPY('F',M,N,A,Lda,U,Ldu)
            ELSE
!
!              Path 5a (M >> N, JOBZ='A')
!              Copy A to VT, generate P**H
!              CWorkspace: need   2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   0
!
               CALL CLACPY('U',N,N,A,Lda,Vt,Ldvt)
               CALL CUNGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(nwork),   &
     &                     Lwork-nwork+1,ierr)
!
!              Copy A to U, generate Q
!              CWorkspace: need   2*N [tauq, taup] + M    [work]
!              CWorkspace: prefer 2*N [tauq, taup] + M*NB [work]
!              RWorkspace: need   0
!
               CALL CLACPY('L',M,N,A,Lda,U,Ldu)
               CALL CUNGBR('Q',M,M,N,U,Ldu,Work(itauq),Work(nwork),     &
     &                     Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
!
               iru = nrwork
               irvt = iru + N*N
               nrwork = irvt + N*N
               CALL SBDSDC('U','I',N,S,Rwork(ie),Rwork(iru),N,          &
     &                     Rwork(irvt),N,dum,idum,Rwork(nrwork),Iwork,  &
     &                     Info)
!
!              Multiply real matrix RWORK(IRVT) by P**H in VT,
!              storing the result in A, copying to VT
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork]
!
               CALL CLARCM(N,N,Rwork(irvt),N,Vt,Ldvt,A,Lda,Rwork(nrwork)&
     &                     )
               CALL CLACPY('F',N,N,A,Lda,Vt,Ldvt)
!
!              Multiply Q in U by real matrix RWORK(IRU), storing the
!              result in A, copying to U
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here
!
               nrwork = irvt
               CALL CLACRM(M,N,U,Ldu,Rwork(iru),N,A,Lda,Rwork(nrwork))
               CALL CLACPY('F',M,N,A,Lda,U,Ldu)
            ENDIF
!
         ELSE
!
!           M .LT. MNTHR2
!
!           Path 6 (M >= N, but not much larger)
!           Reduce to bidiagonal form without QR decomposition
!           Use CUNMBR to compute singular vectors
!
            ie = 1
            nrwork = ie + N
            itauq = 1
            itaup = itauq + N
            nwork = itaup + N
!
!           Bidiagonalize A
!           CWorkspace: need   2*N [tauq, taup] + M        [work]
!           CWorkspace: prefer 2*N [tauq, taup] + (M+N)*NB [work]
!           RWorkspace: need   N [e]
!
            CALL CGEBRD(M,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),  &
     &                  Work(nwork),Lwork-nwork+1,ierr)
            IF ( wntqn ) THEN
!
!              Path 6n (M >= N, JOBZ='N')
!              Compute singular values only
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + BDSPAC
!
               CALL SBDSDC('U','N',N,S,Rwork(ie),dum,1,dum,1,dum,idum,  &
     &                     Rwork(nrwork),Iwork,Info)
            ELSEIF ( wntqo ) THEN
               iu = nwork
               iru = nrwork
               irvt = iru + N*N
               nrwork = irvt + N*N
               IF ( Lwork>=M*N+3*N ) THEN
!
!                 WORK( IU ) is M by N
!
                  ldwrku = M
               ELSE
!
!                 WORK( IU ) is LDWRKU by N
!
                  ldwrku = (Lwork-3*N)/N
               ENDIF
               nwork = iu + ldwrku*N
!
!              Path 6o (M >= N, JOBZ='O')
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
!
               CALL SBDSDC('U','I',N,S,Rwork(ie),Rwork(iru),N,          &
     &                     Rwork(irvt),N,dum,idum,Rwork(nrwork),Iwork,  &
     &                     Info)
!
!              Copy real matrix RWORK(IRVT) to complex matrix VT
!              Overwrite VT by right singular vectors of A
!              CWorkspace: need   2*N [tauq, taup] + N*N [U] + N    [work]
!              CWorkspace: prefer 2*N [tauq, taup] + N*N [U] + N*NB [work]
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]
!
               CALL CLACP2('F',N,N,Rwork(irvt),N,Vt,Ldvt)
               CALL CUNMBR('P','R','C',N,N,N,A,Lda,Work(itaup),Vt,Ldvt, &
     &                     Work(nwork),Lwork-nwork+1,ierr)
!
               IF ( Lwork>=M*N+3*N ) THEN
!
!                 Path 6o-fast
!                 Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
!                 Overwrite WORK(IU) by left singular vectors of A, copying
!                 to A
!                 CWorkspace: need   2*N [tauq, taup] + M*N [U] + N    [work]
!                 CWorkspace: prefer 2*N [tauq, taup] + M*N [U] + N*NB [work]
!                 RWorkspace: need   N [e] + N*N [RU]
!
                  CALL CLASET('F',M,N,CZERO,CZERO,Work(iu),ldwrku)
                  CALL CLACP2('F',N,N,Rwork(iru),N,Work(iu),ldwrku)
                  CALL CUNMBR('Q','L','N',M,N,N,A,Lda,Work(itauq),      &
     &                        Work(iu),ldwrku,Work(nwork),Lwork-nwork+1,&
     &                        ierr)
                  CALL CLACPY('F',M,N,Work(iu),ldwrku,A,Lda)
               ELSE
!
!                 Path 6o-slow
!                 Generate Q in A
!                 CWorkspace: need   2*N [tauq, taup] + N*N [U] + N    [work]
!                 CWorkspace: prefer 2*N [tauq, taup] + N*N [U] + N*NB [work]
!                 RWorkspace: need   0
!
                  CALL CUNGBR('Q',M,N,N,A,Lda,Work(itauq),Work(nwork),  &
     &                        Lwork-nwork+1,ierr)
!
!                 Multiply Q in A by real matrix RWORK(IRU), storing the
!                 result in WORK(IU), copying to A
!                 CWorkspace: need   2*N [tauq, taup] + N*N [U]
!                 CWorkspace: prefer 2*N [tauq, taup] + M*N [U]
!                 RWorkspace: need   N [e] + N*N [RU] + 2*N*N [rwork]
!                 RWorkspace: prefer N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here
!
                  nrwork = irvt
                  DO i = 1 , M , ldwrku
                     chunk = MIN(M-i+1,ldwrku)
                     CALL CLACRM(chunk,N,A(i,1),Lda,Rwork(iru),N,       &
     &                           Work(iu),ldwrku,Rwork(nrwork))
                     CALL CLACPY('F',chunk,N,Work(iu),ldwrku,A(i,1),Lda)
                  ENDDO
               ENDIF
!
            ELSEIF ( wntqs ) THEN
!
!              Path 6s (M >= N, JOBZ='S')
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
!
               iru = nrwork
               irvt = iru + N*N
               nrwork = irvt + N*N
               CALL SBDSDC('U','I',N,S,Rwork(ie),Rwork(iru),N,          &
     &                     Rwork(irvt),N,dum,idum,Rwork(nrwork),Iwork,  &
     &                     Info)
!
!              Copy real matrix RWORK(IRU) to complex matrix U
!              Overwrite U by left singular vectors of A
!              CWorkspace: need   2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]
!
               CALL CLASET('F',M,N,CZERO,CZERO,U,Ldu)
               CALL CLACP2('F',N,N,Rwork(iru),N,U,Ldu)
               CALL CUNMBR('Q','L','N',M,N,N,A,Lda,Work(itauq),U,Ldu,   &
     &                     Work(nwork),Lwork-nwork+1,ierr)
!
!              Copy real matrix RWORK(IRVT) to complex matrix VT
!              Overwrite VT by right singular vectors of A
!              CWorkspace: need   2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]
!
               CALL CLACP2('F',N,N,Rwork(irvt),N,Vt,Ldvt)
               CALL CUNMBR('P','R','C',N,N,N,A,Lda,Work(itaup),Vt,Ldvt, &
     &                     Work(nwork),Lwork-nwork+1,ierr)
            ELSE
!
!              Path 6a (M >= N, JOBZ='A')
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
!
               iru = nrwork
               irvt = iru + N*N
               nrwork = irvt + N*N
               CALL SBDSDC('U','I',N,S,Rwork(ie),Rwork(iru),N,          &
     &                     Rwork(irvt),N,dum,idum,Rwork(nrwork),Iwork,  &
     &                     Info)
!
!              Set the right corner of U to identity matrix
!
               CALL CLASET('F',M,M,CZERO,CZERO,U,Ldu)
               IF ( M>N ) CALL CLASET('F',M-N,M-N,CZERO,CONE,U(N+1,N+1),&
     &                                Ldu)
!
!              Copy real matrix RWORK(IRU) to complex matrix U
!              Overwrite U by left singular vectors of A
!              CWorkspace: need   2*N [tauq, taup] + M    [work]
!              CWorkspace: prefer 2*N [tauq, taup] + M*NB [work]
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]
!
               CALL CLACP2('F',N,N,Rwork(iru),N,U,Ldu)
               CALL CUNMBR('Q','L','N',M,M,N,A,Lda,Work(itauq),U,Ldu,   &
     &                     Work(nwork),Lwork-nwork+1,ierr)
!
!              Copy real matrix RWORK(IRVT) to complex matrix VT
!              Overwrite VT by right singular vectors of A
!              CWorkspace: need   2*N [tauq, taup] + N    [work]
!              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
!              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]
!
               CALL CLACP2('F',N,N,Rwork(irvt),N,Vt,Ldvt)
               CALL CUNMBR('P','R','C',N,N,N,A,Lda,Work(itaup),Vt,Ldvt, &
     &                     Work(nwork),Lwork-nwork+1,ierr)
            ENDIF
!
         ENDIF
!
!
!        A has more columns than rows. If A has sufficiently more
!        columns than rows, first reduce using the LQ decomposition (if
!        sufficient workspace available)
!
      ELSEIF ( N>=mnthr1 ) THEN
!
         IF ( wntqn ) THEN
!
!              Path 1t (N >> M, JOBZ='N')
!              No singular vectors to be computed
!
            itau = 1
            nwork = itau + M
!
!              Compute A=L*Q
!              CWorkspace: need   M [tau] + M    [work]
!              CWorkspace: prefer M [tau] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CGELQF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  ierr)
!
!              Zero out above L
!
            CALL CLASET('U',M-1,M-1,CZERO,CZERO,A(1,2),Lda)
            ie = 1
            itauq = 1
            itaup = itauq + M
            nwork = itaup + M
!
!              Bidiagonalize L in A
!              CWorkspace: need   2*M [tauq, taup] + M      [work]
!              CWorkspace: prefer 2*M [tauq, taup] + 2*M*NB [work]
!              RWorkspace: need   M [e]
!
            CALL CGEBRD(M,M,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),  &
     &                  Work(nwork),Lwork-nwork+1,ierr)
            nrwork = ie + M
!
!              Perform bidiagonal SVD, compute singular values only
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + BDSPAC
!
            CALL SBDSDC('U','N',M,S,Rwork(ie),dum,1,dum,1,dum,idum,     &
     &                  Rwork(nrwork),Iwork,Info)
!
         ELSEIF ( wntqo ) THEN
!
!              Path 2t (N >> M, JOBZ='O')
!              M right singular vectors to be overwritten on A and
!              M left singular vectors to be computed in U
!
            ivt = 1
            ldwkvt = M
!
!              WORK(IVT) is M by M
!
            il = ivt + ldwkvt*M
            IF ( Lwork>=M*N+M*M+3*M ) THEN
!
!                 WORK(IL) M by N
!
               ldwrkl = M
               chunk = N
            ELSE
!
!                 WORK(IL) is M by CHUNK
!
               ldwrkl = M
               chunk = (Lwork-M*M-3*M)/M
            ENDIF
            itau = il + ldwrkl*chunk
            nwork = itau + M
!
!              Compute A=L*Q
!              CWorkspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
!              CWorkspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CGELQF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  ierr)
!
!              Copy L to WORK(IL), zeroing about above it
!
            CALL CLACPY('L',M,M,A,Lda,Work(il),ldwrkl)
            CALL CLASET('U',M-1,M-1,CZERO,CZERO,Work(il+ldwrkl),ldwrkl)
!
!              Generate Q in A
!              CWorkspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
!              CWorkspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CUNGLQ(M,N,M,A,Lda,Work(itau),Work(nwork),             &
     &                  Lwork-nwork+1,ierr)
            ie = 1
            itauq = itau
            itaup = itauq + M
            nwork = itaup + M
!
!              Bidiagonalize L in WORK(IL)
!              CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M      [work]
!              CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + 2*M*NB [work]
!              RWorkspace: need   M [e]
!
            CALL CGEBRD(M,M,Work(il),ldwrkl,S,Rwork(ie),Work(itauq),    &
     &                  Work(itaup),Work(nwork),Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC
!
            iru = ie + M
            irvt = iru + M*M
            nrwork = irvt + M*M
            CALL SBDSDC('U','I',M,S,Rwork(ie),Rwork(iru),M,Rwork(irvt), &
     &                  M,dum,idum,Rwork(nrwork),Iwork,Info)
!
!              Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
!              Overwrite WORK(IU) by the left singular vectors of L
!              CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CLACP2('F',M,M,Rwork(iru),M,U,Ldu)
            CALL CUNMBR('Q','L','N',M,M,M,Work(il),ldwrkl,Work(itauq),U,&
     &                  Ldu,Work(nwork),Lwork-nwork+1,ierr)
!
!              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
!              Overwrite WORK(IVT) by the right singular vectors of L
!              CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CLACP2('F',M,M,Rwork(irvt),M,Work(ivt),ldwkvt)
            CALL CUNMBR('P','R','C',M,M,M,Work(il),ldwrkl,Work(itaup),  &
     &                  Work(ivt),ldwkvt,Work(nwork),Lwork-nwork+1,ierr)
!
!              Multiply right singular vectors of L in WORK(IL) by Q
!              in A, storing result in WORK(IL) and copying to A
!              CWorkspace: need   M*M [VT] + M*M [L]
!              CWorkspace: prefer M*M [VT] + M*N [L]
!              RWorkspace: need   0
!
            DO i = 1 , N , chunk
               blk = MIN(N-i+1,chunk)
               CALL CGEMM('N','N',M,blk,M,CONE,Work(ivt),M,A(1,i),Lda,  &
     &                    CZERO,Work(il),ldwrkl)
               CALL CLACPY('F',M,blk,Work(il),ldwrkl,A(1,i),Lda)
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
!              CWorkspace: need   M*M [L] + M [tau] + M    [work]
!              CWorkspace: prefer M*M [L] + M [tau] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CGELQF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  ierr)
!
!              Copy L to WORK(IL), zeroing out above it
!
            CALL CLACPY('L',M,M,A,Lda,Work(il),ldwrkl)
            CALL CLASET('U',M-1,M-1,CZERO,CZERO,Work(il+ldwrkl),ldwrkl)
!
!              Generate Q in A
!              CWorkspace: need   M*M [L] + M [tau] + M    [work]
!              CWorkspace: prefer M*M [L] + M [tau] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CUNGLQ(M,N,M,A,Lda,Work(itau),Work(nwork),             &
     &                  Lwork-nwork+1,ierr)
            ie = 1
            itauq = itau
            itaup = itauq + M
            nwork = itaup + M
!
!              Bidiagonalize L in WORK(IL)
!              CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M      [work]
!              CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + 2*M*NB [work]
!              RWorkspace: need   M [e]
!
            CALL CGEBRD(M,M,Work(il),ldwrkl,S,Rwork(ie),Work(itauq),    &
     &                  Work(itaup),Work(nwork),Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC
!
            iru = ie + M
            irvt = iru + M*M
            nrwork = irvt + M*M
            CALL SBDSDC('U','I',M,S,Rwork(ie),Rwork(iru),M,Rwork(irvt), &
     &                  M,dum,idum,Rwork(nrwork),Iwork,Info)
!
!              Copy real matrix RWORK(IRU) to complex matrix U
!              Overwrite U by left singular vectors of L
!              CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CLACP2('F',M,M,Rwork(iru),M,U,Ldu)
            CALL CUNMBR('Q','L','N',M,M,M,Work(il),ldwrkl,Work(itauq),U,&
     &                  Ldu,Work(nwork),Lwork-nwork+1,ierr)
!
!              Copy real matrix RWORK(IRVT) to complex matrix VT
!              Overwrite VT by left singular vectors of L
!              CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CLACP2('F',M,M,Rwork(irvt),M,Vt,Ldvt)
            CALL CUNMBR('P','R','C',M,M,M,Work(il),ldwrkl,Work(itaup),  &
     &                  Vt,Ldvt,Work(nwork),Lwork-nwork+1,ierr)
!
!              Copy VT to WORK(IL), multiply right singular vectors of L
!              in WORK(IL) by Q in A, storing result in VT
!              CWorkspace: need   M*M [L]
!              RWorkspace: need   0
!
            CALL CLACPY('F',M,M,Vt,Ldvt,Work(il),ldwrkl)
            CALL CGEMM('N','N',M,N,M,CONE,Work(il),ldwrkl,A,Lda,CZERO,  &
     &                 Vt,Ldvt)
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
!              CWorkspace: need   M*M [VT] + M [tau] + M    [work]
!              CWorkspace: prefer M*M [VT] + M [tau] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CGELQF(M,N,A,Lda,Work(itau),Work(nwork),Lwork-nwork+1, &
     &                  ierr)
            CALL CLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!              Generate Q in VT
!              CWorkspace: need   M*M [VT] + M [tau] + N    [work]
!              CWorkspace: prefer M*M [VT] + M [tau] + N*NB [work]
!              RWorkspace: need   0
!
            CALL CUNGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(nwork),           &
     &                  Lwork-nwork+1,ierr)
!
!              Produce L in A, zeroing out above it
!
            CALL CLASET('U',M-1,M-1,CZERO,CZERO,A(1,2),Lda)
            ie = 1
            itauq = itau
            itaup = itauq + M
            nwork = itaup + M
!
!              Bidiagonalize L in A
!              CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M      [work]
!              CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + 2*M*NB [work]
!              RWorkspace: need   M [e]
!
            CALL CGEBRD(M,M,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),  &
     &                  Work(nwork),Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC
!
            iru = ie + M
            irvt = iru + M*M
            nrwork = irvt + M*M
            CALL SBDSDC('U','I',M,S,Rwork(ie),Rwork(iru),M,Rwork(irvt), &
     &                  M,dum,idum,Rwork(nrwork),Iwork,Info)
!
!              Copy real matrix RWORK(IRU) to complex matrix U
!              Overwrite U by left singular vectors of L
!              CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CLACP2('F',M,M,Rwork(iru),M,U,Ldu)
            CALL CUNMBR('Q','L','N',M,M,M,A,Lda,Work(itauq),U,Ldu,      &
     &                  Work(nwork),Lwork-nwork+1,ierr)
!
!              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
!              Overwrite WORK(IVT) by right singular vectors of L
!              CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CLACP2('F',M,M,Rwork(irvt),M,Work(ivt),ldwkvt)
            CALL CUNMBR('P','R','C',M,M,M,A,Lda,Work(itaup),Work(ivt),  &
     &                  ldwkvt,Work(nwork),Lwork-nwork+1,ierr)
!
!              Multiply right singular vectors of L in WORK(IVT) by
!              Q in VT, storing result in A
!              CWorkspace: need   M*M [VT]
!              RWorkspace: need   0
!
            CALL CGEMM('N','N',M,N,M,CONE,Work(ivt),ldwkvt,Vt,Ldvt,     &
     &                 CZERO,A,Lda)
!
!              Copy right singular vectors of A from A to VT
!
            CALL CLACPY('F',M,N,A,Lda,Vt,Ldvt)
!
         ENDIF
!
      ELSEIF ( N>=mnthr2 ) THEN
!
!           MNTHR2 <= N < MNTHR1
!
!           Path 5t (N >> M, but not as much as MNTHR1)
!           Reduce to bidiagonal form without QR decomposition, use
!           CUNGBR and matrix multiplication to compute singular vectors
!
         ie = 1
         nrwork = ie + M
         itauq = 1
         itaup = itauq + M
         nwork = itaup + M
!
!           Bidiagonalize A
!           CWorkspace: need   2*M [tauq, taup] + N        [work]
!           CWorkspace: prefer 2*M [tauq, taup] + (M+N)*NB [work]
!           RWorkspace: need   M [e]
!
         CALL CGEBRD(M,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),     &
     &               Work(nwork),Lwork-nwork+1,ierr)
!
         IF ( wntqn ) THEN
!
!              Path 5tn (N >> M, JOBZ='N')
!              Compute singular values only
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + BDSPAC
!
            CALL SBDSDC('L','N',M,S,Rwork(ie),dum,1,dum,1,dum,idum,     &
     &                  Rwork(nrwork),Iwork,Info)
         ELSEIF ( wntqo ) THEN
            irvt = nrwork
            iru = irvt + M*M
            nrwork = iru + M*M
            ivt = nwork
!
!              Path 5to (N >> M, JOBZ='O')
!              Copy A to U, generate Q
!              CWorkspace: need   2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CLACPY('L',M,M,A,Lda,U,Ldu)
            CALL CUNGBR('Q',M,M,N,U,Ldu,Work(itauq),Work(nwork),        &
     &                  Lwork-nwork+1,ierr)
!
!              Generate P**H in A
!              CWorkspace: need   2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CUNGBR('P',M,N,M,A,Lda,Work(itaup),Work(nwork),        &
     &                  Lwork-nwork+1,ierr)
!
            ldwkvt = M
            IF ( Lwork>=M*N+3*M ) THEN
!
!                 WORK( IVT ) is M by N
!
               nwork = ivt + ldwkvt*N
               chunk = N
            ELSE
!
!                 WORK( IVT ) is M by CHUNK
!
               chunk = (Lwork-3*M)/M
               nwork = ivt + ldwkvt*chunk
            ENDIF
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC
!
            CALL SBDSDC('L','I',M,S,Rwork(ie),Rwork(iru),M,Rwork(irvt), &
     &                  M,dum,idum,Rwork(nrwork),Iwork,Info)
!
!              Multiply Q in U by real matrix RWORK(IRVT)
!              storing the result in WORK(IVT), copying to U
!              CWorkspace: need   2*M [tauq, taup] + M*M [VT]
!              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork]
!
            CALL CLACRM(M,M,U,Ldu,Rwork(iru),M,Work(ivt),ldwkvt,        &
     &                  Rwork(nrwork))
            CALL CLACPY('F',M,M,Work(ivt),ldwkvt,U,Ldu)
!
!              Multiply RWORK(IRVT) by P**H in A, storing the
!              result in WORK(IVT), copying to A
!              CWorkspace: need   2*M [tauq, taup] + M*M [VT]
!              CWorkspace: prefer 2*M [tauq, taup] + M*N [VT]
!              RWorkspace: need   M [e] + M*M [RVT] + 2*M*M [rwork]
!              RWorkspace: prefer M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here
!
            nrwork = iru
            DO i = 1 , N , chunk
               blk = MIN(N-i+1,chunk)
               CALL CLARCM(M,blk,Rwork(irvt),M,A(1,i),Lda,Work(ivt),    &
     &                     ldwkvt,Rwork(nrwork))
               CALL CLACPY('F',M,blk,Work(ivt),ldwkvt,A(1,i),Lda)
            ENDDO
         ELSEIF ( wntqs ) THEN
!
!              Path 5ts (N >> M, JOBZ='S')
!              Copy A to U, generate Q
!              CWorkspace: need   2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CLACPY('L',M,M,A,Lda,U,Ldu)
            CALL CUNGBR('Q',M,M,N,U,Ldu,Work(itauq),Work(nwork),        &
     &                  Lwork-nwork+1,ierr)
!
!              Copy A to VT, generate P**H
!              CWorkspace: need   2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CLACPY('U',M,N,A,Lda,Vt,Ldvt)
            CALL CUNGBR('P',M,N,M,Vt,Ldvt,Work(itaup),Work(nwork),      &
     &                  Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC
!
            irvt = nrwork
            iru = irvt + M*M
            nrwork = iru + M*M
            CALL SBDSDC('L','I',M,S,Rwork(ie),Rwork(iru),M,Rwork(irvt), &
     &                  M,dum,idum,Rwork(nrwork),Iwork,Info)
!
!              Multiply Q in U by real matrix RWORK(IRU), storing the
!              result in A, copying to U
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork]
!
            CALL CLACRM(M,M,U,Ldu,Rwork(iru),M,A,Lda,Rwork(nrwork))
            CALL CLACPY('F',M,M,A,Lda,U,Ldu)
!
!              Multiply real matrix RWORK(IRVT) by P**H in VT,
!              storing the result in A, copying to VT
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here
!
            nrwork = iru
            CALL CLARCM(M,N,Rwork(irvt),M,Vt,Ldvt,A,Lda,Rwork(nrwork))
            CALL CLACPY('F',M,N,A,Lda,Vt,Ldvt)
         ELSE
!
!              Path 5ta (N >> M, JOBZ='A')
!              Copy A to U, generate Q
!              CWorkspace: need   2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   0
!
            CALL CLACPY('L',M,M,A,Lda,U,Ldu)
            CALL CUNGBR('Q',M,M,N,U,Ldu,Work(itauq),Work(nwork),        &
     &                  Lwork-nwork+1,ierr)
!
!              Copy A to VT, generate P**H
!              CWorkspace: need   2*M [tauq, taup] + N    [work]
!              CWorkspace: prefer 2*M [tauq, taup] + N*NB [work]
!              RWorkspace: need   0
!
            CALL CLACPY('U',M,N,A,Lda,Vt,Ldvt)
            CALL CUNGBR('P',N,N,M,Vt,Ldvt,Work(itaup),Work(nwork),      &
     &                  Lwork-nwork+1,ierr)
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC
!
            irvt = nrwork
            iru = irvt + M*M
            nrwork = iru + M*M
            CALL SBDSDC('L','I',M,S,Rwork(ie),Rwork(iru),M,Rwork(irvt), &
     &                  M,dum,idum,Rwork(nrwork),Iwork,Info)
!
!              Multiply Q in U by real matrix RWORK(IRU), storing the
!              result in A, copying to U
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork]
!
            CALL CLACRM(M,M,U,Ldu,Rwork(iru),M,A,Lda,Rwork(nrwork))
            CALL CLACPY('F',M,M,A,Lda,U,Ldu)
!
!              Multiply real matrix RWORK(IRVT) by P**H in VT,
!              storing the result in A, copying to VT
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here
!
            nrwork = iru
            CALL CLARCM(M,N,Rwork(irvt),M,Vt,Ldvt,A,Lda,Rwork(nrwork))
            CALL CLACPY('F',M,N,A,Lda,Vt,Ldvt)
         ENDIF
!
      ELSE
!
!           N .LT. MNTHR2
!
!           Path 6t (N > M, but not much larger)
!           Reduce to bidiagonal form without LQ decomposition
!           Use CUNMBR to compute singular vectors
!
         ie = 1
         nrwork = ie + M
         itauq = 1
         itaup = itauq + M
         nwork = itaup + M
!
!           Bidiagonalize A
!           CWorkspace: need   2*M [tauq, taup] + N        [work]
!           CWorkspace: prefer 2*M [tauq, taup] + (M+N)*NB [work]
!           RWorkspace: need   M [e]
!
         CALL CGEBRD(M,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),     &
     &               Work(nwork),Lwork-nwork+1,ierr)
         IF ( wntqn ) THEN
!
!              Path 6tn (N > M, JOBZ='N')
!              Compute singular values only
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + BDSPAC
!
            CALL SBDSDC('L','N',M,S,Rwork(ie),dum,1,dum,1,dum,idum,     &
     &                  Rwork(nrwork),Iwork,Info)
         ELSEIF ( wntqo ) THEN
!              Path 6to (N > M, JOBZ='O')
            ldwkvt = M
            ivt = nwork
            IF ( Lwork>=M*N+3*M ) THEN
!
!                 WORK( IVT ) is M by N
!
               CALL CLASET('F',M,N,CZERO,CZERO,Work(ivt),ldwkvt)
               nwork = ivt + ldwkvt*N
            ELSE
!
!                 WORK( IVT ) is M by CHUNK
!
               chunk = (Lwork-3*M)/M
               nwork = ivt + ldwkvt*chunk
            ENDIF
!
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC
!
            irvt = nrwork
            iru = irvt + M*M
            nrwork = iru + M*M
            CALL SBDSDC('L','I',M,S,Rwork(ie),Rwork(iru),M,Rwork(irvt), &
     &                  M,dum,idum,Rwork(nrwork),Iwork,Info)
!
!              Copy real matrix RWORK(IRU) to complex matrix U
!              Overwrite U by left singular vectors of A
!              CWorkspace: need   2*M [tauq, taup] + M*M [VT] + M    [work]
!              CWorkspace: prefer 2*M [tauq, taup] + M*M [VT] + M*NB [work]
!              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU]
!
            CALL CLACP2('F',M,M,Rwork(iru),M,U,Ldu)
            CALL CUNMBR('Q','L','N',M,M,N,A,Lda,Work(itauq),U,Ldu,      &
     &                  Work(nwork),Lwork-nwork+1,ierr)
!
            IF ( Lwork>=M*N+3*M ) THEN
!
!                 Path 6to-fast
!                 Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
!                 Overwrite WORK(IVT) by right singular vectors of A,
!                 copying to A
!                 CWorkspace: need   2*M [tauq, taup] + M*N [VT] + M    [work]
!                 CWorkspace: prefer 2*M [tauq, taup] + M*N [VT] + M*NB [work]
!                 RWorkspace: need   M [e] + M*M [RVT]
!
               CALL CLACP2('F',M,M,Rwork(irvt),M,Work(ivt),ldwkvt)
               CALL CUNMBR('P','R','C',M,N,M,A,Lda,Work(itaup),Work(ivt)&
     &                     ,ldwkvt,Work(nwork),Lwork-nwork+1,ierr)
               CALL CLACPY('F',M,N,Work(ivt),ldwkvt,A,Lda)
            ELSE
!
!                 Path 6to-slow
!                 Generate P**H in A
!                 CWorkspace: need   2*M [tauq, taup] + M*M [VT] + M    [work]
!                 CWorkspace: prefer 2*M [tauq, taup] + M*M [VT] + M*NB [work]
!                 RWorkspace: need   0
!
               CALL CUNGBR('P',M,N,M,A,Lda,Work(itaup),Work(nwork),     &
     &                     Lwork-nwork+1,ierr)
!
!                 Multiply Q in A by real matrix RWORK(IRU), storing the
!                 result in WORK(IU), copying to A
!                 CWorkspace: need   2*M [tauq, taup] + M*M [VT]
!                 CWorkspace: prefer 2*M [tauq, taup] + M*N [VT]
!                 RWorkspace: need   M [e] + M*M [RVT] + 2*M*M [rwork]
!                 RWorkspace: prefer M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here
!
               nrwork = iru
               DO i = 1 , N , chunk
                  blk = MIN(N-i+1,chunk)
                  CALL CLARCM(M,blk,Rwork(irvt),M,A(1,i),Lda,Work(ivt), &
     &                        ldwkvt,Rwork(nrwork))
                  CALL CLACPY('F',M,blk,Work(ivt),ldwkvt,A(1,i),Lda)
               ENDDO
            ENDIF
         ELSEIF ( wntqs ) THEN
!
!              Path 6ts (N > M, JOBZ='S')
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC
!
            irvt = nrwork
            iru = irvt + M*M
            nrwork = iru + M*M
            CALL SBDSDC('L','I',M,S,Rwork(ie),Rwork(iru),M,Rwork(irvt), &
     &                  M,dum,idum,Rwork(nrwork),Iwork,Info)
!
!              Copy real matrix RWORK(IRU) to complex matrix U
!              Overwrite U by left singular vectors of A
!              CWorkspace: need   2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU]
!
            CALL CLACP2('F',M,M,Rwork(iru),M,U,Ldu)
            CALL CUNMBR('Q','L','N',M,M,N,A,Lda,Work(itauq),U,Ldu,      &
     &                  Work(nwork),Lwork-nwork+1,ierr)
!
!              Copy real matrix RWORK(IRVT) to complex matrix VT
!              Overwrite VT by right singular vectors of A
!              CWorkspace: need   2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   M [e] + M*M [RVT]
!
            CALL CLASET('F',M,N,CZERO,CZERO,Vt,Ldvt)
            CALL CLACP2('F',M,M,Rwork(irvt),M,Vt,Ldvt)
            CALL CUNMBR('P','R','C',M,N,M,A,Lda,Work(itaup),Vt,Ldvt,    &
     &                  Work(nwork),Lwork-nwork+1,ierr)
         ELSE
!
!              Path 6ta (N > M, JOBZ='A')
!              Perform bidiagonal SVD, computing left singular vectors
!              of bidiagonal matrix in RWORK(IRU) and computing right
!              singular vectors of bidiagonal matrix in RWORK(IRVT)
!              CWorkspace: need   0
!              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC
!
            irvt = nrwork
            iru = irvt + M*M
            nrwork = iru + M*M
!
            CALL SBDSDC('L','I',M,S,Rwork(ie),Rwork(iru),M,Rwork(irvt), &
     &                  M,dum,idum,Rwork(nrwork),Iwork,Info)
!
!              Copy real matrix RWORK(IRU) to complex matrix U
!              Overwrite U by left singular vectors of A
!              CWorkspace: need   2*M [tauq, taup] + M    [work]
!              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
!              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU]
!
            CALL CLACP2('F',M,M,Rwork(iru),M,U,Ldu)
            CALL CUNMBR('Q','L','N',M,M,N,A,Lda,Work(itauq),U,Ldu,      &
     &                  Work(nwork),Lwork-nwork+1,ierr)
!
!              Set all of VT to identity matrix
!
            CALL CLASET('F',N,N,CZERO,CONE,Vt,Ldvt)
!
!              Copy real matrix RWORK(IRVT) to complex matrix VT
!              Overwrite VT by right singular vectors of A
!              CWorkspace: need   2*M [tauq, taup] + N    [work]
!              CWorkspace: prefer 2*M [tauq, taup] + N*NB [work]
!              RWorkspace: need   M [e] + M*M [RVT]
!
            CALL CLACP2('F',M,M,Rwork(irvt),M,Vt,Ldvt)
            CALL CUNMBR('P','R','C',N,N,M,A,Lda,Work(itaup),Vt,Ldvt,    &
     &                  Work(nwork),Lwork-nwork+1,ierr)
         ENDIF
!
!
      ENDIF
!
!     Undo scaling if necessary
!
      IF ( iscl==1 ) THEN
         IF ( anrm>bignum ) CALL SLASCL('G',0,0,bignum,anrm,minmn,1,S,  &
     &                                  minmn,ierr)
         IF ( Info/=0 .AND. anrm>bignum )                               &
     &        CALL SLASCL('G',0,0,bignum,anrm,minmn-1,1,Rwork(ie),minmn,&
     &        ierr)
         IF ( anrm<smlnum ) CALL SLASCL('G',0,0,smlnum,anrm,minmn,1,S,  &
     &                                  minmn,ierr)
         IF ( Info/=0 .AND. anrm<smlnum )                               &
     &        CALL SLASCL('G',0,0,smlnum,anrm,minmn-1,1,Rwork(ie),minmn,&
     &        ierr)
      ENDIF
!
!     Return optimal workspace in WORK(1)
!
      Work(1) = maxwrk
!
!
!     End of CGESDD
!
      END SUBROUTINE CGESDD
