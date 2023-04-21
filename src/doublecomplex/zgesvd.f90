!*==zgesvd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> ZGESVD computes the singular value decomposition (SVD) for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGESVD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgesvd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgesvd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgesvd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
!                          WORK, LWORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBU, JOBVT
!       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * ), S( * )
!       COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGESVD computes the singular value decomposition (SVD) of a complex
!> M-by-N matrix A, optionally computing the left and/or right singular
!> vectors. The SVD is written
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
!> Note that the routine returns V**H, not V.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBU
!> \verbatim
!>          JOBU is CHARACTER*1
!>          Specifies options for computing all or part of the matrix U:
!>          = 'A':  all M columns of U are returned in array U:
!>          = 'S':  the first min(m,n) columns of U (the left singular
!>                  vectors) are returned in the array U;
!>          = 'O':  the first min(m,n) columns of U (the left singular
!>                  vectors) are overwritten on the array A;
!>          = 'N':  no columns of U (no left singular vectors) are
!>                  computed.
!> \endverbatim
!>
!> \param[in] JOBVT
!> \verbatim
!>          JOBVT is CHARACTER*1
!>          Specifies options for computing all or part of the matrix
!>          V**H:
!>          = 'A':  all N rows of V**H are returned in the array VT;
!>          = 'S':  the first min(m,n) rows of V**H (the right singular
!>                  vectors) are returned in the array VT;
!>          = 'O':  the first min(m,n) rows of V**H (the right singular
!>                  vectors) are overwritten on the array A;
!>          = 'N':  no rows of V**H (no right singular vectors) are
!>                  computed.
!>
!>          JOBVT and JOBU cannot both be 'O'.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit,
!>          if JOBU = 'O',  A is overwritten with the first min(m,n)
!>                          columns of U (the left singular vectors,
!>                          stored columnwise);
!>          if JOBVT = 'O', A is overwritten with the first min(m,n)
!>                          rows of V**H (the right singular vectors,
!>                          stored rowwise);
!>          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
!>                          are destroyed.
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
!>          U is COMPLEX*16 array, dimension (LDU,UCOL)
!>          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
!>          If JOBU = 'A', U contains the M-by-M unitary matrix U;
!>          if JOBU = 'S', U contains the first min(m,n) columns of U
!>          (the left singular vectors, stored columnwise);
!>          if JOBU = 'N' or 'O', U is not referenced.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= 1; if
!>          JOBU = 'S' or 'A', LDU >= M.
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is COMPLEX*16 array, dimension (LDVT,N)
!>          If JOBVT = 'A', VT contains the N-by-N unitary matrix
!>          V**H;
!>          if JOBVT = 'S', VT contains the first min(m,n) rows of
!>          V**H (the right singular vectors, stored rowwise);
!>          if JOBVT = 'N' or 'O', VT is not referenced.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.  LDVT >= 1; if
!>          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          LWORK >=  MAX(1,2*MIN(M,N)+MAX(M,N)).
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
!>          RWORK is DOUBLE PRECISION array, dimension (5*min(M,N))
!>          On exit, if INFO > 0, RWORK(1:MIN(M,N)-1) contains the
!>          unconverged superdiagonal elements of an upper bidiagonal
!>          matrix B whose diagonal is in S (not necessarily sorted).
!>          B satisfies A = U * B * VT, so it has the same singular
!>          values as A, and singular vectors related by U and VT.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if ZBDSQR did not converge, INFO specifies how many
!>                superdiagonals of an intermediate bidiagonal form B
!>                did not converge to zero. See the description of RWORK
!>                above for details.
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
!> \date April 2012
!
!> \ingroup complex16GEsing
!
!  =====================================================================
      SUBROUTINE ZGESVD(Jobu,Jobvt,M,N,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,&
     &                  Rwork,Info)
      IMPLICIT NONE
!*--ZGESVD218
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      CHARACTER Jobu , Jobvt
      INTEGER Info , Lda , Ldu , Ldvt , Lwork , M , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Rwork(*) , S(*)
      COMPLEX*16 A(Lda,*) , U(Ldu,*) , Vt(Ldvt,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0))
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , wntua , wntuas , wntun , wntuo , wntus , wntva , &
     &        wntvas , wntvn , wntvo , wntvs
      INTEGER blk , chunk , i , ie , ierr , ir , irwork , iscl , itau , &
     &        itaup , itauq , iu , iwork , ldwrkr , ldwrku , maxwrk ,   &
     &        minmn , minwrk , mnthr , ncu , ncvt , nru , nrvt , wrkbl
      INTEGER lwork_zgeqrf , lwork_zungqr_n , lwork_zungqr_m ,          &
     &        lwork_zgebrd , lwork_zungbr_p , lwork_zungbr_q ,          &
     &        lwork_zgelqf , lwork_zunglq_n , lwork_zunglq_m
      DOUBLE PRECISION anrm , bignum , eps , smlnum
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION dum(1)
      COMPLEX*16 cdum(1)
!     ..
!     .. External Subroutines ..
      EXTERNAL DLASCL , XERBLA , ZBDSQR , ZGEBRD , ZGELQF , ZGEMM ,     &
     &         ZGEQRF , ZLACPY , ZLASCL , ZLASET , ZUNGBR , ZUNGLQ ,    &
     &         ZUNGQR , ZUNMBR
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL LSAME , ILAENV , DLAMCH , ZLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      minmn = MIN(M,N)
      wntua = LSAME(Jobu,'A')
      wntus = LSAME(Jobu,'S')
      wntuas = wntua .OR. wntus
      wntuo = LSAME(Jobu,'O')
      wntun = LSAME(Jobu,'N')
      wntva = LSAME(Jobvt,'A')
      wntvs = LSAME(Jobvt,'S')
      wntvas = wntva .OR. wntvs
      wntvo = LSAME(Jobvt,'O')
      wntvn = LSAME(Jobvt,'N')
      lquery = (Lwork==-1)
!
      IF ( .NOT.(wntua .OR. wntus .OR. wntuo .OR. wntun) ) THEN
         Info = -1
      ELSEIF ( .NOT.(wntva .OR. wntvs .OR. wntvo .OR. wntvn) .OR.       &
     &         (wntvo .AND. wntuo) ) THEN
         Info = -2
      ELSEIF ( M<0 ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -6
      ELSEIF ( Ldu<1 .OR. (wntuas .AND. Ldu<M) ) THEN
         Info = -9
      ELSEIF ( Ldvt<1 .OR. (wntva .AND. Ldvt<N) .OR.                    &
     &         (wntvs .AND. Ldvt<minmn) ) THEN
         Info = -11
      ENDIF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
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
!           Space needed for ZBDSQR is BDSPAC = 5*N
!
            mnthr = ILAENV(6,'ZGESVD',Jobu//Jobvt,M,N,0,0)
!           Compute space needed for ZGEQRF
            CALL ZGEQRF(M,N,A,Lda,cdum(1),cdum(1),-1,ierr)
            lwork_zgeqrf = INT(cdum(1))
!           Compute space needed for ZUNGQR
            CALL ZUNGQR(M,N,N,A,Lda,cdum(1),cdum(1),-1,ierr)
            lwork_zungqr_n = INT(cdum(1))
            CALL ZUNGQR(M,M,N,A,Lda,cdum(1),cdum(1),-1,ierr)
            lwork_zungqr_m = INT(cdum(1))
!           Compute space needed for ZGEBRD
            CALL ZGEBRD(N,N,A,Lda,S,dum(1),cdum(1),cdum(1),cdum(1),-1,  &
     &                  ierr)
            lwork_zgebrd = INT(cdum(1))
!           Compute space needed for ZUNGBR
            CALL ZUNGBR('P',N,N,N,A,Lda,cdum(1),cdum(1),-1,ierr)
            lwork_zungbr_p = INT(cdum(1))
            CALL ZUNGBR('Q',N,N,N,A,Lda,cdum(1),cdum(1),-1,ierr)
            lwork_zungbr_q = INT(cdum(1))
!
            IF ( M<mnthr ) THEN
!
!              Path 10 (M at least N, but not much larger)
!
               CALL ZGEBRD(M,N,A,Lda,S,dum(1),cdum(1),cdum(1),cdum(1),  &
     &                     -1,ierr)
               lwork_zgebrd = INT(cdum(1))
               maxwrk = 2*N + lwork_zgebrd
               IF ( wntus .OR. wntuo ) THEN
                  CALL ZUNGBR('Q',M,N,N,A,Lda,cdum(1),cdum(1),-1,ierr)
                  lwork_zungbr_q = INT(cdum(1))
                  maxwrk = MAX(maxwrk,2*N+lwork_zungbr_q)
               ENDIF
               IF ( wntua ) THEN
                  CALL ZUNGBR('Q',M,M,N,A,Lda,cdum(1),cdum(1),-1,ierr)
                  lwork_zungbr_q = INT(cdum(1))
                  maxwrk = MAX(maxwrk,2*N+lwork_zungbr_q)
               ENDIF
               IF ( .NOT.wntvn ) maxwrk = MAX(maxwrk,2*N+lwork_zungbr_p)
               minwrk = 2*N + M
            ELSEIF ( wntun ) THEN
!
!                 Path 1 (M much larger than N, JOBU='N')
!
               maxwrk = N + lwork_zgeqrf
               maxwrk = MAX(maxwrk,2*N+lwork_zgebrd)
               IF ( wntvo .OR. wntvas )                                 &
     &              maxwrk = MAX(maxwrk,2*N+lwork_zungbr_p)
               minwrk = 3*N
            ELSEIF ( wntuo .AND. wntvn ) THEN
!
!                 Path 2 (M much larger than N, JOBU='O', JOBVT='N')
!
               wrkbl = N + lwork_zgeqrf
               wrkbl = MAX(wrkbl,N+lwork_zungqr_n)
               wrkbl = MAX(wrkbl,2*N+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_q)
               maxwrk = MAX(N*N+wrkbl,N*N+M*N)
               minwrk = 2*N + M
            ELSEIF ( wntuo .AND. wntvas ) THEN
!
!                 Path 3 (M much larger than N, JOBU='O', JOBVT='S' or
!                 'A')
!
               wrkbl = N + lwork_zgeqrf
               wrkbl = MAX(wrkbl,N+lwork_zungqr_n)
               wrkbl = MAX(wrkbl,2*N+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_q)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_p)
               maxwrk = MAX(N*N+wrkbl,N*N+M*N)
               minwrk = 2*N + M
            ELSEIF ( wntus .AND. wntvn ) THEN
!
!                 Path 4 (M much larger than N, JOBU='S', JOBVT='N')
!
               wrkbl = N + lwork_zgeqrf
               wrkbl = MAX(wrkbl,N+lwork_zungqr_n)
               wrkbl = MAX(wrkbl,2*N+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_q)
               maxwrk = N*N + wrkbl
               minwrk = 2*N + M
            ELSEIF ( wntus .AND. wntvo ) THEN
!
!                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
!
               wrkbl = N + lwork_zgeqrf
               wrkbl = MAX(wrkbl,N+lwork_zungqr_n)
               wrkbl = MAX(wrkbl,2*N+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_q)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_p)
               maxwrk = 2*N*N + wrkbl
               minwrk = 2*N + M
            ELSEIF ( wntus .AND. wntvas ) THEN
!
!                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' or
!                 'A')
!
               wrkbl = N + lwork_zgeqrf
               wrkbl = MAX(wrkbl,N+lwork_zungqr_n)
               wrkbl = MAX(wrkbl,2*N+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_q)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_p)
               maxwrk = N*N + wrkbl
               minwrk = 2*N + M
            ELSEIF ( wntua .AND. wntvn ) THEN
!
!                 Path 7 (M much larger than N, JOBU='A', JOBVT='N')
!
               wrkbl = N + lwork_zgeqrf
               wrkbl = MAX(wrkbl,N+lwork_zungqr_m)
               wrkbl = MAX(wrkbl,2*N+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_q)
               maxwrk = N*N + wrkbl
               minwrk = 2*N + M
            ELSEIF ( wntua .AND. wntvo ) THEN
!
!                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
!
               wrkbl = N + lwork_zgeqrf
               wrkbl = MAX(wrkbl,N+lwork_zungqr_m)
               wrkbl = MAX(wrkbl,2*N+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_q)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_p)
               maxwrk = 2*N*N + wrkbl
               minwrk = 2*N + M
            ELSEIF ( wntua .AND. wntvas ) THEN
!
!                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' or
!                 'A')
!
               wrkbl = N + lwork_zgeqrf
               wrkbl = MAX(wrkbl,N+lwork_zungqr_m)
               wrkbl = MAX(wrkbl,2*N+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_q)
               wrkbl = MAX(wrkbl,2*N+lwork_zungbr_p)
               maxwrk = N*N + wrkbl
               minwrk = 2*N + M
            ENDIF
         ELSEIF ( minmn>0 ) THEN
!
!           Space needed for ZBDSQR is BDSPAC = 5*M
!
            mnthr = ILAENV(6,'ZGESVD',Jobu//Jobvt,M,N,0,0)
!           Compute space needed for ZGELQF
            CALL ZGELQF(M,N,A,Lda,cdum(1),cdum(1),-1,ierr)
            lwork_zgelqf = INT(cdum(1))
!           Compute space needed for ZUNGLQ
            CALL ZUNGLQ(N,N,M,cdum(1),N,cdum(1),cdum(1),-1,ierr)
            lwork_zunglq_n = INT(cdum(1))
            CALL ZUNGLQ(M,N,M,A,Lda,cdum(1),cdum(1),-1,ierr)
            lwork_zunglq_m = INT(cdum(1))
!           Compute space needed for ZGEBRD
            CALL ZGEBRD(M,M,A,Lda,S,dum(1),cdum(1),cdum(1),cdum(1),-1,  &
     &                  ierr)
            lwork_zgebrd = INT(cdum(1))
!            Compute space needed for ZUNGBR P
            CALL ZUNGBR('P',M,M,M,A,N,cdum(1),cdum(1),-1,ierr)
            lwork_zungbr_p = INT(cdum(1))
!           Compute space needed for ZUNGBR Q
            CALL ZUNGBR('Q',M,M,M,A,N,cdum(1),cdum(1),-1,ierr)
            lwork_zungbr_q = INT(cdum(1))
            IF ( N<mnthr ) THEN
!
!              Path 10t(N greater than M, but not much larger)
!
               CALL ZGEBRD(M,N,A,Lda,S,dum(1),cdum(1),cdum(1),cdum(1),  &
     &                     -1,ierr)
               lwork_zgebrd = INT(cdum(1))
               maxwrk = 2*M + lwork_zgebrd
               IF ( wntvs .OR. wntvo ) THEN
!                Compute space needed for ZUNGBR P
                  CALL ZUNGBR('P',M,N,M,A,N,cdum(1),cdum(1),-1,ierr)
                  lwork_zungbr_p = INT(cdum(1))
                  maxwrk = MAX(maxwrk,2*M+lwork_zungbr_p)
               ENDIF
               IF ( wntva ) THEN
                  CALL ZUNGBR('P',N,N,M,A,N,cdum(1),cdum(1),-1,ierr)
                  lwork_zungbr_p = INT(cdum(1))
                  maxwrk = MAX(maxwrk,2*M+lwork_zungbr_p)
               ENDIF
               IF ( .NOT.wntun ) maxwrk = MAX(maxwrk,2*M+lwork_zungbr_q)
               minwrk = 2*M + N
            ELSEIF ( wntvn ) THEN
!
!                 Path 1t(N much larger than M, JOBVT='N')
!
               maxwrk = M + lwork_zgelqf
               maxwrk = MAX(maxwrk,2*M+lwork_zgebrd)
               IF ( wntuo .OR. wntuas )                                 &
     &              maxwrk = MAX(maxwrk,2*M+lwork_zungbr_q)
               minwrk = 3*M
            ELSEIF ( wntvo .AND. wntun ) THEN
!
!                 Path 2t(N much larger than M, JOBU='N', JOBVT='O')
!
               wrkbl = M + lwork_zgelqf
               wrkbl = MAX(wrkbl,M+lwork_zunglq_m)
               wrkbl = MAX(wrkbl,2*M+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_p)
               maxwrk = MAX(M*M+wrkbl,M*M+M*N)
               minwrk = 2*M + N
            ELSEIF ( wntvo .AND. wntuas ) THEN
!
!                 Path 3t(N much larger than M, JOBU='S' or 'A',
!                 JOBVT='O')
!
               wrkbl = M + lwork_zgelqf
               wrkbl = MAX(wrkbl,M+lwork_zunglq_m)
               wrkbl = MAX(wrkbl,2*M+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_p)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_q)
               maxwrk = MAX(M*M+wrkbl,M*M+M*N)
               minwrk = 2*M + N
            ELSEIF ( wntvs .AND. wntun ) THEN
!
!                 Path 4t(N much larger than M, JOBU='N', JOBVT='S')
!
               wrkbl = M + lwork_zgelqf
               wrkbl = MAX(wrkbl,M+lwork_zunglq_m)
               wrkbl = MAX(wrkbl,2*M+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_p)
               maxwrk = M*M + wrkbl
               minwrk = 2*M + N
            ELSEIF ( wntvs .AND. wntuo ) THEN
!
!                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
!
               wrkbl = M + lwork_zgelqf
               wrkbl = MAX(wrkbl,M+lwork_zunglq_m)
               wrkbl = MAX(wrkbl,2*M+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_p)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_q)
               maxwrk = 2*M*M + wrkbl
               minwrk = 2*M + N
            ELSEIF ( wntvs .AND. wntuas ) THEN
!
!                 Path 6t(N much larger than M, JOBU='S' or 'A',
!                 JOBVT='S')
!
               wrkbl = M + lwork_zgelqf
               wrkbl = MAX(wrkbl,M+lwork_zunglq_m)
               wrkbl = MAX(wrkbl,2*M+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_p)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_q)
               maxwrk = M*M + wrkbl
               minwrk = 2*M + N
            ELSEIF ( wntva .AND. wntun ) THEN
!
!                 Path 7t(N much larger than M, JOBU='N', JOBVT='A')
!
               wrkbl = M + lwork_zgelqf
               wrkbl = MAX(wrkbl,M+lwork_zunglq_n)
               wrkbl = MAX(wrkbl,2*M+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_p)
               maxwrk = M*M + wrkbl
               minwrk = 2*M + N
            ELSEIF ( wntva .AND. wntuo ) THEN
!
!                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
!
               wrkbl = M + lwork_zgelqf
               wrkbl = MAX(wrkbl,M+lwork_zunglq_n)
               wrkbl = MAX(wrkbl,2*M+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_p)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_q)
               maxwrk = 2*M*M + wrkbl
               minwrk = 2*M + N
            ELSEIF ( wntva .AND. wntuas ) THEN
!
!                 Path 9t(N much larger than M, JOBU='S' or 'A',
!                 JOBVT='A')
!
               wrkbl = M + lwork_zgelqf
               wrkbl = MAX(wrkbl,M+lwork_zunglq_n)
               wrkbl = MAX(wrkbl,2*M+lwork_zgebrd)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_p)
               wrkbl = MAX(wrkbl,2*M+lwork_zungbr_q)
               maxwrk = M*M + wrkbl
               minwrk = 2*M + N
            ENDIF
         ENDIF
         maxwrk = MAX(maxwrk,minwrk)
         Work(1) = maxwrk
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -13
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGESVD',-Info)
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
      anrm = ZLANGE('M',M,N,A,Lda,dum)
      iscl = 0
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         iscl = 1
         CALL ZLASCL('G',0,0,anrm,smlnum,M,N,A,Lda,ierr)
      ELSEIF ( anrm>bignum ) THEN
         iscl = 1
         CALL ZLASCL('G',0,0,anrm,bignum,M,N,A,Lda,ierr)
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
!           Path 10 (M at least N, but not much larger)
!           Reduce to bidiagonal form without QR decomposition
!
            ie = 1
            itauq = 1
            itaup = itauq + N
            iwork = itaup + N
!
!           Bidiagonalize A
!           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB)
!           (RWorkspace: need N)
!
            CALL ZGEBRD(M,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),  &
     &                  Work(iwork),Lwork-iwork+1,ierr)
            IF ( wntuas ) THEN
!
!              If left singular vectors desired in U, copy result to U
!              and generate left bidiagonalizing vectors in U
!              (CWorkspace: need 2*N+NCU, prefer 2*N+NCU*NB)
!              (RWorkspace: 0)
!
               CALL ZLACPY('L',M,N,A,Lda,U,Ldu)
               IF ( wntus ) ncu = N
               IF ( wntua ) ncu = M
               CALL ZUNGBR('Q',M,ncu,N,U,Ldu,Work(itauq),Work(iwork),   &
     &                     Lwork-iwork+1,ierr)
            ENDIF
            IF ( wntvas ) THEN
!
!              If right singular vectors desired in VT, copy result to
!              VT and generate right bidiagonalizing vectors in VT
!              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!              (RWorkspace: 0)
!
               CALL ZLACPY('U',N,N,A,Lda,Vt,Ldvt)
               CALL ZUNGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),   &
     &                     Lwork-iwork+1,ierr)
            ENDIF
!
!              If left singular vectors desired in A, generate left
!              bidiagonalizing vectors in A
!              (CWorkspace: need 3*N, prefer 2*N+N*NB)
!              (RWorkspace: 0)
!
            IF ( wntuo ) CALL ZUNGBR('Q',M,N,N,A,Lda,Work(itauq),       &
     &                               Work(iwork),Lwork-iwork+1,ierr)
!
!              If right singular vectors desired in A, generate right
!              bidiagonalizing vectors in A
!              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!              (RWorkspace: 0)
!
            IF ( wntvo ) CALL ZUNGBR('P',N,N,N,A,Lda,Work(itaup),       &
     &                               Work(iwork),Lwork-iwork+1,ierr)
            irwork = ie + N
            IF ( wntuas .OR. wntuo ) nru = M
            IF ( wntun ) nru = 0
            IF ( wntvas .OR. wntvo ) ncvt = N
            IF ( wntvn ) ncvt = 0
            IF ( (.NOT.wntuo) .AND. (.NOT.wntvo) ) THEN
!
!              Perform bidiagonal QR iteration, if desired, computing
!              left singular vectors in U and computing right singular
!              vectors in VT
!              (CWorkspace: 0)
!              (RWorkspace: need BDSPAC)
!
               CALL ZBDSQR('U',N,ncvt,nru,0,S,Rwork(ie),Vt,Ldvt,U,Ldu,  &
     &                     cdum,1,Rwork(irwork),Info)
            ELSEIF ( (.NOT.wntuo) .AND. wntvo ) THEN
!
!              Perform bidiagonal QR iteration, if desired, computing
!              left singular vectors in U and computing right singular
!              vectors in A
!              (CWorkspace: 0)
!              (RWorkspace: need BDSPAC)
!
               CALL ZBDSQR('U',N,ncvt,nru,0,S,Rwork(ie),A,Lda,U,Ldu,    &
     &                     cdum,1,Rwork(irwork),Info)
            ELSE
!
!              Perform bidiagonal QR iteration, if desired, computing
!              left singular vectors in A and computing right singular
!              vectors in VT
!              (CWorkspace: 0)
!              (RWorkspace: need BDSPAC)
!
               CALL ZBDSQR('U',N,ncvt,nru,0,S,Rwork(ie),Vt,Ldvt,A,Lda,  &
     &                     cdum,1,Rwork(irwork),Info)
            ENDIF
!
         ELSEIF ( wntun ) THEN
!
!              Path 1 (M much larger than N, JOBU='N')
!              No left singular vectors to be computed
!
            itau = 1
            iwork = itau + N
!
!              Compute A=Q*R
!              (CWorkspace: need 2*N, prefer N+N*NB)
!              (RWorkspace: need 0)
!
            CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),Lwork-iwork+1, &
     &                  ierr)
!
!              Zero out below R
!
            IF ( N>1 ) CALL ZLASET('L',N-1,N-1,CZERO,CZERO,A(2,1),Lda)
            ie = 1
            itauq = 1
            itaup = itauq + N
            iwork = itaup + N
!
!              Bidiagonalize R in A
!              (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
!              (RWorkspace: need N)
!
            CALL ZGEBRD(N,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),  &
     &                  Work(iwork),Lwork-iwork+1,ierr)
            ncvt = 0
            IF ( wntvo .OR. wntvas ) THEN
!
!                 If right singular vectors desired, generate P'.
!                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGBR('P',N,N,N,A,Lda,Work(itaup),Work(iwork),     &
     &                     Lwork-iwork+1,ierr)
               ncvt = N
            ENDIF
            irwork = ie + N
!
!              Perform bidiagonal QR iteration, computing right
!              singular vectors of A in A if desired
!              (CWorkspace: 0)
!              (RWorkspace: need BDSPAC)
!
            CALL ZBDSQR('U',N,ncvt,0,0,S,Rwork(ie),A,Lda,cdum,1,cdum,1, &
     &                  Rwork(irwork),Info)
!
!              If right singular vectors desired in VT, copy them there
!
            IF ( wntvas ) CALL ZLACPY('F',N,N,A,Lda,Vt,Ldvt)
!
         ELSEIF ( wntuo .AND. wntvn ) THEN
!
!              Path 2 (M much larger than N, JOBU='O', JOBVT='N')
!              N left singular vectors to be overwritten on A and
!              no right singular vectors to be computed
!
            IF ( Lwork>=N*N+3*N ) THEN
!
!                 Sufficient workspace for a fast algorithm
!
               ir = 1
               IF ( Lwork>=MAX(wrkbl,Lda*N)+Lda*N ) THEN
!
!                    WORK(IU) is LDA by N, WORK(IR) is LDA by N
!
                  ldwrku = Lda
                  ldwrkr = Lda
               ELSEIF ( Lwork>=MAX(wrkbl,Lda*N)+N*N ) THEN
!
!                    WORK(IU) is LDA by N, WORK(IR) is N by N
!
                  ldwrku = Lda
                  ldwrkr = N
               ELSE
!
!                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N
!
                  ldwrku = (Lwork-N*N)/N
                  ldwrkr = N
               ENDIF
               itau = ir + ldwrkr*N
               iwork = itau + N
!
!                 Compute A=Q*R
!                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
!                 (RWorkspace: 0)
!
               CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),            &
     &                     Lwork-iwork+1,ierr)
!
!                 Copy R to WORK(IR) and zero out below it
!
               CALL ZLACPY('U',N,N,A,Lda,Work(ir),ldwrkr)
               CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Work(ir+1),ldwrkr)
!
!                 Generate Q in A
!                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGQR(M,N,N,A,Lda,Work(itau),Work(iwork),          &
     &                     Lwork-iwork+1,ierr)
               ie = 1
               itauq = itau
               itaup = itauq + N
               iwork = itaup + N
!
!                 Bidiagonalize R in WORK(IR)
!                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
!                 (RWorkspace: need N)
!
               CALL ZGEBRD(N,N,Work(ir),ldwrkr,S,Rwork(ie),Work(itauq), &
     &                     Work(itaup),Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate left vectors bidiagonalizing R
!                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
!                 (RWorkspace: need 0)
!
               CALL ZUNGBR('Q',N,N,N,Work(ir),ldwrkr,Work(itauq),       &
     &                     Work(iwork),Lwork-iwork+1,ierr)
               irwork = ie + N
!
!                 Perform bidiagonal QR iteration, computing left
!                 singular vectors of R in WORK(IR)
!                 (CWorkspace: need N*N)
!                 (RWorkspace: need BDSPAC)
!
               CALL ZBDSQR('U',N,0,N,0,S,Rwork(ie),cdum,1,Work(ir),     &
     &                     ldwrkr,cdum,1,Rwork(irwork),Info)
               iu = itauq
!
!                 Multiply Q in A by left singular vectors of R in
!                 WORK(IR), storing result in WORK(IU) and copying to A
!                 (CWorkspace: need N*N+N, prefer N*N+M*N)
!                 (RWorkspace: 0)
!
               DO i = 1 , M , ldwrku
                  chunk = MIN(M-i+1,ldwrku)
                  CALL ZGEMM('N','N',chunk,N,N,CONE,A(i,1),Lda,Work(ir),&
     &                       ldwrkr,CZERO,Work(iu),ldwrku)
                  CALL ZLACPY('F',chunk,N,Work(iu),ldwrku,A(i,1),Lda)
               ENDDO
!
            ELSE
!
!                 Insufficient workspace for a fast algorithm
!
               ie = 1
               itauq = 1
               itaup = itauq + N
               iwork = itaup + N
!
!                 Bidiagonalize A
!                 (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB)
!                 (RWorkspace: N)
!
               CALL ZGEBRD(M,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup)&
     &                     ,Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate left vectors bidiagonalizing A
!                 (CWorkspace: need 3*N, prefer 2*N+N*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGBR('Q',M,N,N,A,Lda,Work(itauq),Work(iwork),     &
     &                     Lwork-iwork+1,ierr)
               irwork = ie + N
!
!                 Perform bidiagonal QR iteration, computing left
!                 singular vectors of A in A
!                 (CWorkspace: need 0)
!                 (RWorkspace: need BDSPAC)
!
               CALL ZBDSQR('U',N,0,M,0,S,Rwork(ie),cdum,1,A,Lda,cdum,1, &
     &                     Rwork(irwork),Info)
!
            ENDIF
!
         ELSEIF ( wntuo .AND. wntvas ) THEN
!
!              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A')
!              N left singular vectors to be overwritten on A and
!              N right singular vectors to be computed in VT
!
            IF ( Lwork>=N*N+3*N ) THEN
!
!                 Sufficient workspace for a fast algorithm
!
               ir = 1
               IF ( Lwork>=MAX(wrkbl,Lda*N)+Lda*N ) THEN
!
!                    WORK(IU) is LDA by N and WORK(IR) is LDA by N
!
                  ldwrku = Lda
                  ldwrkr = Lda
               ELSEIF ( Lwork>=MAX(wrkbl,Lda*N)+N*N ) THEN
!
!                    WORK(IU) is LDA by N and WORK(IR) is N by N
!
                  ldwrku = Lda
                  ldwrkr = N
               ELSE
!
!                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N
!
                  ldwrku = (Lwork-N*N)/N
                  ldwrkr = N
               ENDIF
               itau = ir + ldwrkr*N
               iwork = itau + N
!
!                 Compute A=Q*R
!                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
!                 (RWorkspace: 0)
!
               CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),            &
     &                     Lwork-iwork+1,ierr)
!
!                 Copy R to VT, zeroing out below it
!
               CALL ZLACPY('U',N,N,A,Lda,Vt,Ldvt)
               IF ( N>1 ) CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Vt(2,1),  &
     &                                Ldvt)
!
!                 Generate Q in A
!                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGQR(M,N,N,A,Lda,Work(itau),Work(iwork),          &
     &                     Lwork-iwork+1,ierr)
               ie = 1
               itauq = itau
               itaup = itauq + N
               iwork = itaup + N
!
!                 Bidiagonalize R in VT, copying result to WORK(IR)
!                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
!                 (RWorkspace: need N)
!
               CALL ZGEBRD(N,N,Vt,Ldvt,S,Rwork(ie),Work(itauq),         &
     &                     Work(itaup),Work(iwork),Lwork-iwork+1,ierr)
               CALL ZLACPY('L',N,N,Vt,Ldvt,Work(ir),ldwrkr)
!
!                 Generate left vectors bidiagonalizing R in WORK(IR)
!                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGBR('Q',N,N,N,Work(ir),ldwrkr,Work(itauq),       &
     &                     Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate right vectors bidiagonalizing R in VT
!                 (CWorkspace: need N*N+3*N-1, prefer N*N+2*N+(N-1)*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),   &
     &                     Lwork-iwork+1,ierr)
               irwork = ie + N
!
!                 Perform bidiagonal QR iteration, computing left
!                 singular vectors of R in WORK(IR) and computing right
!                 singular vectors of R in VT
!                 (CWorkspace: need N*N)
!                 (RWorkspace: need BDSPAC)
!
               CALL ZBDSQR('U',N,N,N,0,S,Rwork(ie),Vt,Ldvt,Work(ir),    &
     &                     ldwrkr,cdum,1,Rwork(irwork),Info)
               iu = itauq
!
!                 Multiply Q in A by left singular vectors of R in
!                 WORK(IR), storing result in WORK(IU) and copying to A
!                 (CWorkspace: need N*N+N, prefer N*N+M*N)
!                 (RWorkspace: 0)
!
               DO i = 1 , M , ldwrku
                  chunk = MIN(M-i+1,ldwrku)
                  CALL ZGEMM('N','N',chunk,N,N,CONE,A(i,1),Lda,Work(ir),&
     &                       ldwrkr,CZERO,Work(iu),ldwrku)
                  CALL ZLACPY('F',chunk,N,Work(iu),ldwrku,A(i,1),Lda)
               ENDDO
!
            ELSE
!
!                 Insufficient workspace for a fast algorithm
!
               itau = 1
               iwork = itau + N
!
!                 Compute A=Q*R
!                 (CWorkspace: need 2*N, prefer N+N*NB)
!                 (RWorkspace: 0)
!
               CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),            &
     &                     Lwork-iwork+1,ierr)
!
!                 Copy R to VT, zeroing out below it
!
               CALL ZLACPY('U',N,N,A,Lda,Vt,Ldvt)
               IF ( N>1 ) CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Vt(2,1),  &
     &                                Ldvt)
!
!                 Generate Q in A
!                 (CWorkspace: need 2*N, prefer N+N*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGQR(M,N,N,A,Lda,Work(itau),Work(iwork),          &
     &                     Lwork-iwork+1,ierr)
               ie = 1
               itauq = itau
               itaup = itauq + N
               iwork = itaup + N
!
!                 Bidiagonalize R in VT
!                 (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
!                 (RWorkspace: N)
!
               CALL ZGEBRD(N,N,Vt,Ldvt,S,Rwork(ie),Work(itauq),         &
     &                     Work(itaup),Work(iwork),Lwork-iwork+1,ierr)
!
!                 Multiply Q in A by left vectors bidiagonalizing R
!                 (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNMBR('Q','R','N',M,N,N,Vt,Ldvt,Work(itauq),A,Lda, &
     &                     Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate right vectors bidiagonalizing R in VT
!                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),   &
     &                     Lwork-iwork+1,ierr)
               irwork = ie + N
!
!                 Perform bidiagonal QR iteration, computing left
!                 singular vectors of A in A and computing right
!                 singular vectors of A in VT
!                 (CWorkspace: 0)
!                 (RWorkspace: need BDSPAC)
!
               CALL ZBDSQR('U',N,N,M,0,S,Rwork(ie),Vt,Ldvt,A,Lda,cdum,1,&
     &                     Rwork(irwork),Info)
!
            ENDIF
!
         ELSEIF ( wntus ) THEN
!
            IF ( wntvn ) THEN
!
!                 Path 4 (M much larger than N, JOBU='S', JOBVT='N')
!                 N left singular vectors to be computed in U and
!                 no right singular vectors to be computed
!
               IF ( Lwork>=N*N+3*N ) THEN
!
!                    Sufficient workspace for a fast algorithm
!
                  ir = 1
                  IF ( Lwork>=wrkbl+Lda*N ) THEN
!
!                       WORK(IR) is LDA by N
!
                     ldwrkr = Lda
                  ELSE
!
!                       WORK(IR) is N by N
!
                     ldwrkr = N
                  ENDIF
                  itau = ir + ldwrkr*N
                  iwork = itau + N
!
!                    Compute A=Q*R
!                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R to WORK(IR), zeroing out below it
!
                  CALL ZLACPY('U',N,N,A,Lda,Work(ir),ldwrkr)
                  CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Work(ir+1),ldwrkr)
!
!                    Generate Q in A
!                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGQR(M,N,N,A,Lda,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in WORK(IR)
!                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
!                    (RWorkspace: need N)
!
                  CALL ZGEBRD(N,N,Work(ir),ldwrkr,S,Rwork(ie),          &
     &                        Work(itauq),Work(itaup),Work(iwork),      &
     &                        Lwork-iwork+1,ierr)
!
!                    Generate left vectors bidiagonalizing R in WORK(IR)
!                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',N,N,N,Work(ir),ldwrkr,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  irwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of R in WORK(IR)
!                    (CWorkspace: need N*N)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',N,0,N,0,S,Rwork(ie),cdum,1,Work(ir),  &
     &                        ldwrkr,cdum,1,Rwork(irwork),Info)
!
!                    Multiply Q in A by left singular vectors of R in
!                    WORK(IR), storing result in U
!                    (CWorkspace: need N*N)
!                    (RWorkspace: 0)
!
                  CALL ZGEMM('N','N',M,N,N,CONE,A,Lda,Work(ir),ldwrkr,  &
     &                       CZERO,U,Ldu)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (CWorkspace: need 2*N, prefer N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (CWorkspace: need 2*N, prefer N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGQR(M,N,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Zero out below R in A
!
                  IF ( N>1 ) CALL ZLASET('L',N-1,N-1,CZERO,CZERO,A(2,1),&
     &                 Lda)
!
!                    Bidiagonalize R in A
!                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
!                    (RWorkspace: need N)
!
                  CALL ZGEBRD(N,N,A,Lda,S,Rwork(ie),Work(itauq),        &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply Q in U by left vectors bidiagonalizing R
!                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNMBR('Q','R','N',M,N,N,A,Lda,Work(itauq),U,Ldu,&
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  irwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U
!                    (CWorkspace: 0)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',N,0,M,0,S,Rwork(ie),cdum,1,U,Ldu,cdum,&
     &                        1,Rwork(irwork),Info)
!
               ENDIF
!
            ELSEIF ( wntvo ) THEN
!
!                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
!                 N left singular vectors to be computed in U and
!                 N right singular vectors to be overwritten on A
!
               IF ( Lwork>=2*N*N+3*N ) THEN
!
!                    Sufficient workspace for a fast algorithm
!
                  iu = 1
                  IF ( Lwork>=wrkbl+2*Lda*N ) THEN
!
!                       WORK(IU) is LDA by N and WORK(IR) is LDA by N
!
                     ldwrku = Lda
                     ir = iu + ldwrku*N
                     ldwrkr = Lda
                  ELSEIF ( Lwork>=wrkbl+(Lda+N)*N ) THEN
!
!                       WORK(IU) is LDA by N and WORK(IR) is N by N
!
                     ldwrku = Lda
                     ir = iu + ldwrku*N
                     ldwrkr = N
                  ELSE
!
!                       WORK(IU) is N by N and WORK(IR) is N by N
!
                     ldwrku = N
                     ir = iu + ldwrku*N
                     ldwrkr = N
                  ENDIF
                  itau = ir + ldwrkr*N
                  iwork = itau + N
!
!                    Compute A=Q*R
!                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R to WORK(IU), zeroing out below it
!
                  CALL ZLACPY('U',N,N,A,Lda,Work(iu),ldwrku)
                  CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Work(iu+1),ldwrku)
!
!                    Generate Q in A
!                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGQR(M,N,N,A,Lda,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in WORK(IU), copying result to
!                    WORK(IR)
!                    (CWorkspace: need   2*N*N+3*N,
!                                 prefer 2*N*N+2*N+2*N*NB)
!                    (RWorkspace: need   N)
!
                  CALL ZGEBRD(N,N,Work(iu),ldwrku,S,Rwork(ie),          &
     &                        Work(itauq),Work(itaup),Work(iwork),      &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('U',N,N,Work(iu),ldwrku,Work(ir),ldwrkr)
!
!                    Generate left bidiagonalizing vectors in WORK(IU)
!                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',N,N,N,Work(iu),ldwrku,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in WORK(IR)
!                    (CWorkspace: need   2*N*N+3*N-1,
!                                 prefer 2*N*N+2*N+(N-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',N,N,N,Work(ir),ldwrkr,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  irwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of R in WORK(IU) and computing
!                    right singular vectors of R in WORK(IR)
!                    (CWorkspace: need 2*N*N)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',N,N,N,0,S,Rwork(ie),Work(ir),ldwrkr,  &
     &                        Work(iu),ldwrku,cdum,1,Rwork(irwork),Info)
!
!                    Multiply Q in A by left singular vectors of R in
!                    WORK(IU), storing result in U
!                    (CWorkspace: need N*N)
!                    (RWorkspace: 0)
!
                  CALL ZGEMM('N','N',M,N,N,CONE,A,Lda,Work(iu),ldwrku,  &
     &                       CZERO,U,Ldu)
!
!                    Copy right singular vectors of R to A
!                    (CWorkspace: need N*N)
!                    (RWorkspace: 0)
!
                  CALL ZLACPY('F',N,N,Work(ir),ldwrkr,A,Lda)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (CWorkspace: need 2*N, prefer N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (CWorkspace: need 2*N, prefer N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGQR(M,N,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Zero out below R in A
!
                  IF ( N>1 ) CALL ZLASET('L',N-1,N-1,CZERO,CZERO,A(2,1),&
     &                 Lda)
!
!                    Bidiagonalize R in A
!                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
!                    (RWorkspace: need N)
!
                  CALL ZGEBRD(N,N,A,Lda,S,Rwork(ie),Work(itauq),        &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply Q in U by left vectors bidiagonalizing R
!                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNMBR('Q','R','N',M,N,N,A,Lda,Work(itauq),U,Ldu,&
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right vectors bidiagonalizing R in A
!                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',N,N,N,A,Lda,Work(itaup),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  irwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U and computing right
!                    singular vectors of A in A
!                    (CWorkspace: 0)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',N,N,M,0,S,Rwork(ie),A,Lda,U,Ldu,cdum, &
     &                        1,Rwork(irwork),Info)
!
               ENDIF
!
            ELSEIF ( wntvas ) THEN
!
!                 Path 6 (M much larger than N, JOBU='S', JOBVT='S'
!                         or 'A')
!                 N left singular vectors to be computed in U and
!                 N right singular vectors to be computed in VT
!
               IF ( Lwork>=N*N+3*N ) THEN
!
!                    Sufficient workspace for a fast algorithm
!
                  iu = 1
                  IF ( Lwork>=wrkbl+Lda*N ) THEN
!
!                       WORK(IU) is LDA by N
!
                     ldwrku = Lda
                  ELSE
!
!                       WORK(IU) is N by N
!
                     ldwrku = N
                  ENDIF
                  itau = iu + ldwrku*N
                  iwork = itau + N
!
!                    Compute A=Q*R
!                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R to WORK(IU), zeroing out below it
!
                  CALL ZLACPY('U',N,N,A,Lda,Work(iu),ldwrku)
                  CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Work(iu+1),ldwrku)
!
!                    Generate Q in A
!                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGQR(M,N,N,A,Lda,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in WORK(IU), copying result to VT
!                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
!                    (RWorkspace: need N)
!
                  CALL ZGEBRD(N,N,Work(iu),ldwrku,S,Rwork(ie),          &
     &                        Work(itauq),Work(itaup),Work(iwork),      &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('U',N,N,Work(iu),ldwrku,Vt,Ldvt)
!
!                    Generate left bidiagonalizing vectors in WORK(IU)
!                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',N,N,N,Work(iu),ldwrku,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in VT
!                    (CWorkspace: need   N*N+3*N-1,
!                                 prefer N*N+2*N+(N-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),&
     &                        Lwork-iwork+1,ierr)
                  irwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of R in WORK(IU) and computing
!                    right singular vectors of R in VT
!                    (CWorkspace: need N*N)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',N,N,N,0,S,Rwork(ie),Vt,Ldvt,Work(iu), &
     &                        ldwrku,cdum,1,Rwork(irwork),Info)
!
!                    Multiply Q in A by left singular vectors of R in
!                    WORK(IU), storing result in U
!                    (CWorkspace: need N*N)
!                    (RWorkspace: 0)
!
                  CALL ZGEMM('N','N',M,N,N,CONE,A,Lda,Work(iu),ldwrku,  &
     &                       CZERO,U,Ldu)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (CWorkspace: need 2*N, prefer N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (CWorkspace: need 2*N, prefer N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGQR(M,N,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R to VT, zeroing out below it
!
                  CALL ZLACPY('U',N,N,A,Lda,Vt,Ldvt)
                  IF ( N>1 ) CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Vt(2,1)&
     &                 ,Ldvt)
                  ie = 1
                  itauq = itau
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in VT
!                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
!                    (RWorkspace: need N)
!
                  CALL ZGEBRD(N,N,Vt,Ldvt,S,Rwork(ie),Work(itauq),      &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply Q in U by left bidiagonalizing vectors
!                    in VT
!                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNMBR('Q','R','N',M,N,N,Vt,Ldvt,Work(itauq),U,  &
     &                        Ldu,Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in VT
!                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),&
     &                        Lwork-iwork+1,ierr)
                  irwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U and computing right
!                    singular vectors of A in VT
!                    (CWorkspace: 0)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',N,N,M,0,S,Rwork(ie),Vt,Ldvt,U,Ldu,    &
     &                        cdum,1,Rwork(irwork),Info)
!
               ENDIF
!
            ENDIF
!
         ELSEIF ( wntua ) THEN
!
            IF ( wntvn ) THEN
!
!                 Path 7 (M much larger than N, JOBU='A', JOBVT='N')
!                 M left singular vectors to be computed in U and
!                 no right singular vectors to be computed
!
               IF ( Lwork>=N*N+MAX(N+M,3*N) ) THEN
!
!                    Sufficient workspace for a fast algorithm
!
                  ir = 1
                  IF ( Lwork>=wrkbl+Lda*N ) THEN
!
!                       WORK(IR) is LDA by N
!
                     ldwrkr = Lda
                  ELSE
!
!                       WORK(IR) is N by N
!
                     ldwrkr = N
                  ENDIF
                  itau = ir + ldwrkr*N
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Copy R to WORK(IR), zeroing out below it
!
                  CALL ZLACPY('U',N,N,A,Lda,Work(ir),ldwrkr)
                  CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Work(ir+1),ldwrkr)
!
!                    Generate Q in U
!                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGQR(M,M,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in WORK(IR)
!                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
!                    (RWorkspace: need N)
!
                  CALL ZGEBRD(N,N,Work(ir),ldwrkr,S,Rwork(ie),          &
     &                        Work(itauq),Work(itaup),Work(iwork),      &
     &                        Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in WORK(IR)
!                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',N,N,N,Work(ir),ldwrkr,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  irwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of R in WORK(IR)
!                    (CWorkspace: need N*N)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',N,0,N,0,S,Rwork(ie),cdum,1,Work(ir),  &
     &                        ldwrkr,cdum,1,Rwork(irwork),Info)
!
!                    Multiply Q in U by left singular vectors of R in
!                    WORK(IR), storing result in A
!                    (CWorkspace: need N*N)
!                    (RWorkspace: 0)
!
                  CALL ZGEMM('N','N',M,N,N,CONE,U,Ldu,Work(ir),ldwrkr,  &
     &                       CZERO,A,Lda)
!
!                    Copy left singular vectors of A from A to U
!
                  CALL ZLACPY('F',M,N,A,Lda,U,Ldu)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (CWorkspace: need 2*N, prefer N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (CWorkspace: need N+M, prefer N+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGQR(M,M,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Zero out below R in A
!
                  IF ( N>1 ) CALL ZLASET('L',N-1,N-1,CZERO,CZERO,A(2,1),&
     &                 Lda)
!
!                    Bidiagonalize R in A
!                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
!                    (RWorkspace: need N)
!
                  CALL ZGEBRD(N,N,A,Lda,S,Rwork(ie),Work(itauq),        &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply Q in U by left bidiagonalizing vectors
!                    in A
!                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNMBR('Q','R','N',M,N,N,A,Lda,Work(itauq),U,Ldu,&
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  irwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U
!                    (CWorkspace: 0)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',N,0,M,0,S,Rwork(ie),cdum,1,U,Ldu,cdum,&
     &                        1,Rwork(irwork),Info)
!
               ENDIF
!
            ELSEIF ( wntvo ) THEN
!
!                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
!                 M left singular vectors to be computed in U and
!                 N right singular vectors to be overwritten on A
!
               IF ( Lwork>=2*N*N+MAX(N+M,3*N) ) THEN
!
!                    Sufficient workspace for a fast algorithm
!
                  iu = 1
                  IF ( Lwork>=wrkbl+2*Lda*N ) THEN
!
!                       WORK(IU) is LDA by N and WORK(IR) is LDA by N
!
                     ldwrku = Lda
                     ir = iu + ldwrku*N
                     ldwrkr = Lda
                  ELSEIF ( Lwork>=wrkbl+(Lda+N)*N ) THEN
!
!                       WORK(IU) is LDA by N and WORK(IR) is N by N
!
                     ldwrku = Lda
                     ir = iu + ldwrku*N
                     ldwrkr = N
                  ELSE
!
!                       WORK(IU) is N by N and WORK(IR) is N by N
!
                     ldwrku = N
                     ir = iu + ldwrku*N
                     ldwrkr = N
                  ENDIF
                  itau = ir + ldwrkr*N
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (CWorkspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGQR(M,M,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R to WORK(IU), zeroing out below it
!
                  CALL ZLACPY('U',N,N,A,Lda,Work(iu),ldwrku)
                  CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Work(iu+1),ldwrku)
                  ie = 1
                  itauq = itau
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in WORK(IU), copying result to
!                    WORK(IR)
!                    (CWorkspace: need   2*N*N+3*N,
!                                 prefer 2*N*N+2*N+2*N*NB)
!                    (RWorkspace: need   N)
!
                  CALL ZGEBRD(N,N,Work(iu),ldwrku,S,Rwork(ie),          &
     &                        Work(itauq),Work(itaup),Work(iwork),      &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('U',N,N,Work(iu),ldwrku,Work(ir),ldwrkr)
!
!                    Generate left bidiagonalizing vectors in WORK(IU)
!                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',N,N,N,Work(iu),ldwrku,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in WORK(IR)
!                    (CWorkspace: need   2*N*N+3*N-1,
!                                 prefer 2*N*N+2*N+(N-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',N,N,N,Work(ir),ldwrkr,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  irwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of R in WORK(IU) and computing
!                    right singular vectors of R in WORK(IR)
!                    (CWorkspace: need 2*N*N)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',N,N,N,0,S,Rwork(ie),Work(ir),ldwrkr,  &
     &                        Work(iu),ldwrku,cdum,1,Rwork(irwork),Info)
!
!                    Multiply Q in U by left singular vectors of R in
!                    WORK(IU), storing result in A
!                    (CWorkspace: need N*N)
!                    (RWorkspace: 0)
!
                  CALL ZGEMM('N','N',M,N,N,CONE,U,Ldu,Work(iu),ldwrku,  &
     &                       CZERO,A,Lda)
!
!                    Copy left singular vectors of A from A to U
!
                  CALL ZLACPY('F',M,N,A,Lda,U,Ldu)
!
!                    Copy right singular vectors of R from WORK(IR) to A
!
                  CALL ZLACPY('F',N,N,Work(ir),ldwrkr,A,Lda)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (CWorkspace: need 2*N, prefer N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (CWorkspace: need N+M, prefer N+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGQR(M,M,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Zero out below R in A
!
                  IF ( N>1 ) CALL ZLASET('L',N-1,N-1,CZERO,CZERO,A(2,1),&
     &                 Lda)
!
!                    Bidiagonalize R in A
!                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
!                    (RWorkspace: need N)
!
                  CALL ZGEBRD(N,N,A,Lda,S,Rwork(ie),Work(itauq),        &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply Q in U by left bidiagonalizing vectors
!                    in A
!                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNMBR('Q','R','N',M,N,N,A,Lda,Work(itauq),U,Ldu,&
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in A
!                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',N,N,N,A,Lda,Work(itaup),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  irwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U and computing right
!                    singular vectors of A in A
!                    (CWorkspace: 0)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',N,N,M,0,S,Rwork(ie),A,Lda,U,Ldu,cdum, &
     &                        1,Rwork(irwork),Info)
!
               ENDIF
!
            ELSEIF ( wntvas ) THEN
!
!                 Path 9 (M much larger than N, JOBU='A', JOBVT='S'
!                         or 'A')
!                 M left singular vectors to be computed in U and
!                 N right singular vectors to be computed in VT
!
               IF ( Lwork>=N*N+MAX(N+M,3*N) ) THEN
!
!                    Sufficient workspace for a fast algorithm
!
                  iu = 1
                  IF ( Lwork>=wrkbl+Lda*N ) THEN
!
!                       WORK(IU) is LDA by N
!
                     ldwrku = Lda
                  ELSE
!
!                       WORK(IU) is N by N
!
                     ldwrku = N
                  ENDIF
                  itau = iu + ldwrku*N
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGQR(M,M,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R to WORK(IU), zeroing out below it
!
                  CALL ZLACPY('U',N,N,A,Lda,Work(iu),ldwrku)
                  CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Work(iu+1),ldwrku)
                  ie = 1
                  itauq = itau
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in WORK(IU), copying result to VT
!                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
!                    (RWorkspace: need N)
!
                  CALL ZGEBRD(N,N,Work(iu),ldwrku,S,Rwork(ie),          &
     &                        Work(itauq),Work(itaup),Work(iwork),      &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('U',N,N,Work(iu),ldwrku,Vt,Ldvt)
!
!                    Generate left bidiagonalizing vectors in WORK(IU)
!                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',N,N,N,Work(iu),ldwrku,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in VT
!                    (CWorkspace: need   N*N+3*N-1,
!                                 prefer N*N+2*N+(N-1)*NB)
!                    (RWorkspace: need   0)
!
                  CALL ZUNGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),&
     &                        Lwork-iwork+1,ierr)
                  irwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of R in WORK(IU) and computing
!                    right singular vectors of R in VT
!                    (CWorkspace: need N*N)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',N,N,N,0,S,Rwork(ie),Vt,Ldvt,Work(iu), &
     &                        ldwrku,cdum,1,Rwork(irwork),Info)
!
!                    Multiply Q in U by left singular vectors of R in
!                    WORK(IU), storing result in A
!                    (CWorkspace: need N*N)
!                    (RWorkspace: 0)
!
                  CALL ZGEMM('N','N',M,N,N,CONE,U,Ldu,Work(iu),ldwrku,  &
     &                       CZERO,A,Lda)
!
!                    Copy left singular vectors of A from A to U
!
                  CALL ZLACPY('F',M,N,A,Lda,U,Ldu)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (CWorkspace: need 2*N, prefer N+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (CWorkspace: need N+M, prefer N+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGQR(M,M,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R from A to VT, zeroing out below it
!
                  CALL ZLACPY('U',N,N,A,Lda,Vt,Ldvt)
                  IF ( N>1 ) CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Vt(2,1)&
     &                 ,Ldvt)
                  ie = 1
                  itauq = itau
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in VT
!                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
!                    (RWorkspace: need N)
!
                  CALL ZGEBRD(N,N,Vt,Ldvt,S,Rwork(ie),Work(itauq),      &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply Q in U by left bidiagonalizing vectors
!                    in VT
!                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNMBR('Q','R','N',M,N,N,Vt,Ldvt,Work(itauq),U,  &
     &                        Ldu,Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in VT
!                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),&
     &                        Lwork-iwork+1,ierr)
                  irwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U and computing right
!                    singular vectors of A in VT
!                    (CWorkspace: 0)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',N,N,M,0,S,Rwork(ie),Vt,Ldvt,U,Ldu,    &
     &                        cdum,1,Rwork(irwork),Info)
!
               ENDIF
!
            ENDIF
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
         IF ( wntvn ) THEN
!
!              Path 1t(N much larger than M, JOBVT='N')
!              No right singular vectors to be computed
!
            itau = 1
!
!
            iwork = itau + M
!
!              Compute A=L*Q
!              (CWorkspace: need 2*M, prefer M+M*NB)
!              (RWorkspace: 0)
!
            CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),Lwork-iwork+1, &
     &                  ierr)
!
!              Zero out above L
!
            CALL ZLASET('U',M-1,M-1,CZERO,CZERO,A(1,2),Lda)
            ie = 1
            itauq = 1
            itaup = itauq + M
            iwork = itaup + M
!
!              Bidiagonalize L in A
!              (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
!              (RWorkspace: need M)
!
            CALL ZGEBRD(M,M,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),  &
     &                  Work(iwork),Lwork-iwork+1,ierr)
!
!                 If left singular vectors desired, generate Q
!                 (CWorkspace: need 3*M, prefer 2*M+M*NB)
!                 (RWorkspace: 0)
!
            IF ( wntuo .OR. wntuas )                                    &
     &           CALL ZUNGBR('Q',M,M,M,A,Lda,Work(itauq),Work(iwork),   &
     &           Lwork-iwork+1,ierr)
            irwork = ie + M
            nru = 0
            IF ( wntuo .OR. wntuas ) nru = M
!
!              Perform bidiagonal QR iteration, computing left singular
!              vectors of A in A if desired
!              (CWorkspace: 0)
!              (RWorkspace: need BDSPAC)
!
            CALL ZBDSQR('U',M,0,nru,0,S,Rwork(ie),cdum,1,A,Lda,cdum,1,  &
     &                  Rwork(irwork),Info)
!
!              If left singular vectors desired in U, copy them there
!
            IF ( wntuas ) CALL ZLACPY('F',M,M,A,Lda,U,Ldu)
!
         ELSEIF ( wntvo .AND. wntun ) THEN
!
!              Path 2t(N much larger than M, JOBU='N', JOBVT='O')
!              M right singular vectors to be overwritten on A and
!              no left singular vectors to be computed
!
            IF ( Lwork>=M*M+3*M ) THEN
!
!                 Sufficient workspace for a fast algorithm
!
               ir = 1
               IF ( Lwork>=MAX(wrkbl,Lda*N)+Lda*M ) THEN
!
!                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
!
                  ldwrku = Lda
                  chunk = N
                  ldwrkr = Lda
               ELSEIF ( Lwork>=MAX(wrkbl,Lda*N)+M*M ) THEN
!
!                    WORK(IU) is LDA by N and WORK(IR) is M by M
!
                  ldwrku = Lda
                  chunk = N
                  ldwrkr = M
               ELSE
!
!                    WORK(IU) is M by CHUNK and WORK(IR) is M by M
!
                  ldwrku = M
                  chunk = (Lwork-M*M)/M
                  ldwrkr = M
               ENDIF
               itau = ir + ldwrkr*M
               iwork = itau + M
!
!                 Compute A=L*Q
!                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
!                 (RWorkspace: 0)
!
               CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),            &
     &                     Lwork-iwork+1,ierr)
!
!                 Copy L to WORK(IR) and zero out above it
!
               CALL ZLACPY('L',M,M,A,Lda,Work(ir),ldwrkr)
               CALL ZLASET('U',M-1,M-1,CZERO,CZERO,Work(ir+ldwrkr),     &
     &                     ldwrkr)
!
!                 Generate Q in A
!                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGLQ(M,N,M,A,Lda,Work(itau),Work(iwork),          &
     &                     Lwork-iwork+1,ierr)
               ie = 1
               itauq = itau
               itaup = itauq + M
               iwork = itaup + M
!
!                 Bidiagonalize L in WORK(IR)
!                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
!                 (RWorkspace: need M)
!
               CALL ZGEBRD(M,M,Work(ir),ldwrkr,S,Rwork(ie),Work(itauq), &
     &                     Work(itaup),Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate right vectors bidiagonalizing L
!                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGBR('P',M,M,M,Work(ir),ldwrkr,Work(itaup),       &
     &                     Work(iwork),Lwork-iwork+1,ierr)
               irwork = ie + M
!
!                 Perform bidiagonal QR iteration, computing right
!                 singular vectors of L in WORK(IR)
!                 (CWorkspace: need M*M)
!                 (RWorkspace: need BDSPAC)
!
               CALL ZBDSQR('U',M,M,0,0,S,Rwork(ie),Work(ir),ldwrkr,cdum,&
     &                     1,cdum,1,Rwork(irwork),Info)
               iu = itauq
!
!                 Multiply right singular vectors of L in WORK(IR) by Q
!                 in A, storing result in WORK(IU) and copying to A
!                 (CWorkspace: need M*M+M, prefer M*M+M*N)
!                 (RWorkspace: 0)
!
               DO i = 1 , N , chunk
                  blk = MIN(N-i+1,chunk)
                  CALL ZGEMM('N','N',M,blk,M,CONE,Work(ir),ldwrkr,A(1,i)&
     &                       ,Lda,CZERO,Work(iu),ldwrku)
                  CALL ZLACPY('F',M,blk,Work(iu),ldwrku,A(1,i),Lda)
               ENDDO
!
            ELSE
!
!                 Insufficient workspace for a fast algorithm
!
               ie = 1
               itauq = 1
               itaup = itauq + M
               iwork = itaup + M
!
!                 Bidiagonalize A
!                 (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
!                 (RWorkspace: need M)
!
               CALL ZGEBRD(M,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup)&
     &                     ,Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate right vectors bidiagonalizing A
!                 (CWorkspace: need 3*M, prefer 2*M+M*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGBR('P',M,N,M,A,Lda,Work(itaup),Work(iwork),     &
     &                     Lwork-iwork+1,ierr)
               irwork = ie + M
!
!                 Perform bidiagonal QR iteration, computing right
!                 singular vectors of A in A
!                 (CWorkspace: 0)
!                 (RWorkspace: need BDSPAC)
!
               CALL ZBDSQR('L',M,N,0,0,S,Rwork(ie),A,Lda,cdum,1,cdum,1, &
     &                     Rwork(irwork),Info)
!
            ENDIF
!
         ELSEIF ( wntvo .AND. wntuas ) THEN
!
!              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
!              M right singular vectors to be overwritten on A and
!              M left singular vectors to be computed in U
!
            IF ( Lwork>=M*M+3*M ) THEN
!
!                 Sufficient workspace for a fast algorithm
!
               ir = 1
               IF ( Lwork>=MAX(wrkbl,Lda*N)+Lda*M ) THEN
!
!                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
!
                  ldwrku = Lda
                  chunk = N
                  ldwrkr = Lda
               ELSEIF ( Lwork>=MAX(wrkbl,Lda*N)+M*M ) THEN
!
!                    WORK(IU) is LDA by N and WORK(IR) is M by M
!
                  ldwrku = Lda
                  chunk = N
                  ldwrkr = M
               ELSE
!
!                    WORK(IU) is M by CHUNK and WORK(IR) is M by M
!
                  ldwrku = M
                  chunk = (Lwork-M*M)/M
                  ldwrkr = M
               ENDIF
               itau = ir + ldwrkr*M
               iwork = itau + M
!
!                 Compute A=L*Q
!                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
!                 (RWorkspace: 0)
!
               CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),            &
     &                     Lwork-iwork+1,ierr)
!
!                 Copy L to U, zeroing about above it
!
               CALL ZLACPY('L',M,M,A,Lda,U,Ldu)
               CALL ZLASET('U',M-1,M-1,CZERO,CZERO,U(1,2),Ldu)
!
!                 Generate Q in A
!                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGLQ(M,N,M,A,Lda,Work(itau),Work(iwork),          &
     &                     Lwork-iwork+1,ierr)
               ie = 1
               itauq = itau
               itaup = itauq + M
               iwork = itaup + M
!
!                 Bidiagonalize L in U, copying result to WORK(IR)
!                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
!                 (RWorkspace: need M)
!
               CALL ZGEBRD(M,M,U,Ldu,S,Rwork(ie),Work(itauq),Work(itaup)&
     &                     ,Work(iwork),Lwork-iwork+1,ierr)
               CALL ZLACPY('U',M,M,U,Ldu,Work(ir),ldwrkr)
!
!                 Generate right vectors bidiagonalizing L in WORK(IR)
!                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGBR('P',M,M,M,Work(ir),ldwrkr,Work(itaup),       &
     &                     Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate left vectors bidiagonalizing L in U
!                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGBR('Q',M,M,M,U,Ldu,Work(itauq),Work(iwork),     &
     &                     Lwork-iwork+1,ierr)
               irwork = ie + M
!
!                 Perform bidiagonal QR iteration, computing left
!                 singular vectors of L in U, and computing right
!                 singular vectors of L in WORK(IR)
!                 (CWorkspace: need M*M)
!                 (RWorkspace: need BDSPAC)
!
               CALL ZBDSQR('U',M,M,M,0,S,Rwork(ie),Work(ir),ldwrkr,U,   &
     &                     Ldu,cdum,1,Rwork(irwork),Info)
               iu = itauq
!
!                 Multiply right singular vectors of L in WORK(IR) by Q
!                 in A, storing result in WORK(IU) and copying to A
!                 (CWorkspace: need M*M+M, prefer M*M+M*N))
!                 (RWorkspace: 0)
!
               DO i = 1 , N , chunk
                  blk = MIN(N-i+1,chunk)
                  CALL ZGEMM('N','N',M,blk,M,CONE,Work(ir),ldwrkr,A(1,i)&
     &                       ,Lda,CZERO,Work(iu),ldwrku)
                  CALL ZLACPY('F',M,blk,Work(iu),ldwrku,A(1,i),Lda)
               ENDDO
!
            ELSE
!
!                 Insufficient workspace for a fast algorithm
!
               itau = 1
               iwork = itau + M
!
!                 Compute A=L*Q
!                 (CWorkspace: need 2*M, prefer M+M*NB)
!                 (RWorkspace: 0)
!
               CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),            &
     &                     Lwork-iwork+1,ierr)
!
!                 Copy L to U, zeroing out above it
!
               CALL ZLACPY('L',M,M,A,Lda,U,Ldu)
               CALL ZLASET('U',M-1,M-1,CZERO,CZERO,U(1,2),Ldu)
!
!                 Generate Q in A
!                 (CWorkspace: need 2*M, prefer M+M*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGLQ(M,N,M,A,Lda,Work(itau),Work(iwork),          &
     &                     Lwork-iwork+1,ierr)
               ie = 1
               itauq = itau
               itaup = itauq + M
               iwork = itaup + M
!
!                 Bidiagonalize L in U
!                 (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
!                 (RWorkspace: need M)
!
               CALL ZGEBRD(M,M,U,Ldu,S,Rwork(ie),Work(itauq),Work(itaup)&
     &                     ,Work(iwork),Lwork-iwork+1,ierr)
!
!                 Multiply right vectors bidiagonalizing L by Q in A
!                 (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNMBR('P','L','C',M,N,M,U,Ldu,Work(itaup),A,Lda,   &
     &                     Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate left vectors bidiagonalizing L in U
!                 (CWorkspace: need 3*M, prefer 2*M+M*NB)
!                 (RWorkspace: 0)
!
               CALL ZUNGBR('Q',M,M,M,U,Ldu,Work(itauq),Work(iwork),     &
     &                     Lwork-iwork+1,ierr)
               irwork = ie + M
!
!                 Perform bidiagonal QR iteration, computing left
!                 singular vectors of A in U and computing right
!                 singular vectors of A in A
!                 (CWorkspace: 0)
!                 (RWorkspace: need BDSPAC)
!
               CALL ZBDSQR('U',M,N,M,0,S,Rwork(ie),A,Lda,U,Ldu,cdum,1,  &
     &                     Rwork(irwork),Info)
!
            ENDIF
!
         ELSEIF ( wntvs ) THEN
!
            IF ( wntun ) THEN
!
!                 Path 4t(N much larger than M, JOBU='N', JOBVT='S')
!                 M right singular vectors to be computed in VT and
!                 no left singular vectors to be computed
!
               IF ( Lwork>=M*M+3*M ) THEN
!
!                    Sufficient workspace for a fast algorithm
!
                  ir = 1
                  IF ( Lwork>=wrkbl+Lda*M ) THEN
!
!                       WORK(IR) is LDA by M
!
                     ldwrkr = Lda
                  ELSE
!
!                       WORK(IR) is M by M
!
                     ldwrkr = M
                  ENDIF
                  itau = ir + ldwrkr*M
                  iwork = itau + M
!
!                    Compute A=L*Q
!                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to WORK(IR), zeroing out above it
!
                  CALL ZLACPY('L',M,M,A,Lda,Work(ir),ldwrkr)
                  CALL ZLASET('U',M-1,M-1,CZERO,CZERO,Work(ir+ldwrkr),  &
     &                        ldwrkr)
!
!                    Generate Q in A
!                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGLQ(M,N,M,A,Lda,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in WORK(IR)
!                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
!                    (RWorkspace: need M)
!
                  CALL ZGEBRD(M,M,Work(ir),ldwrkr,S,Rwork(ie),          &
     &                        Work(itauq),Work(itaup),Work(iwork),      &
     &                        Lwork-iwork+1,ierr)
!
!                    Generate right vectors bidiagonalizing L in
!                    WORK(IR)
!                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',M,M,M,Work(ir),ldwrkr,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  irwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing right
!                    singular vectors of L in WORK(IR)
!                    (CWorkspace: need M*M)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',M,M,0,0,S,Rwork(ie),Work(ir),ldwrkr,  &
     &                        cdum,1,cdum,1,Rwork(irwork),Info)
!
!                    Multiply right singular vectors of L in WORK(IR) by
!                    Q in A, storing result in VT
!                    (CWorkspace: need M*M)
!                    (RWorkspace: 0)
!
                  CALL ZGEMM('N','N',M,N,M,CONE,Work(ir),ldwrkr,A,Lda,  &
     &                       CZERO,Vt,Ldvt)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + M
!
!                    Compute A=L*Q
!                    (CWorkspace: need 2*M, prefer M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy result to VT
!
                  CALL ZLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (CWorkspace: need 2*M, prefer M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGLQ(M,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Zero out above L in A
!
                  CALL ZLASET('U',M-1,M-1,CZERO,CZERO,A(1,2),Lda)
!
!                    Bidiagonalize L in A
!                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
!                    (RWorkspace: need M)
!
                  CALL ZGEBRD(M,M,A,Lda,S,Rwork(ie),Work(itauq),        &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply right vectors bidiagonalizing L by Q in VT
!                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNMBR('P','L','C',M,N,M,A,Lda,Work(itaup),Vt,   &
     &                        Ldvt,Work(iwork),Lwork-iwork+1,ierr)
                  irwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing right
!                    singular vectors of A in VT
!                    (CWorkspace: 0)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',M,N,0,0,S,Rwork(ie),Vt,Ldvt,cdum,1,   &
     &                        cdum,1,Rwork(irwork),Info)
!
               ENDIF
!
            ELSEIF ( wntuo ) THEN
!
!                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
!                 M right singular vectors to be computed in VT and
!                 M left singular vectors to be overwritten on A
!
               IF ( Lwork>=2*M*M+3*M ) THEN
!
!                    Sufficient workspace for a fast algorithm
!
                  iu = 1
                  IF ( Lwork>=wrkbl+2*Lda*M ) THEN
!
!                       WORK(IU) is LDA by M and WORK(IR) is LDA by M
!
                     ldwrku = Lda
                     ir = iu + ldwrku*M
                     ldwrkr = Lda
                  ELSEIF ( Lwork>=wrkbl+(Lda+M)*M ) THEN
!
!                       WORK(IU) is LDA by M and WORK(IR) is M by M
!
                     ldwrku = Lda
                     ir = iu + ldwrku*M
                     ldwrkr = M
                  ELSE
!
!                       WORK(IU) is M by M and WORK(IR) is M by M
!
                     ldwrku = M
                     ir = iu + ldwrku*M
                     ldwrkr = M
                  ENDIF
                  itau = ir + ldwrkr*M
                  iwork = itau + M
!
!                    Compute A=L*Q
!                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to WORK(IU), zeroing out below it
!
                  CALL ZLACPY('L',M,M,A,Lda,Work(iu),ldwrku)
                  CALL ZLASET('U',M-1,M-1,CZERO,CZERO,Work(iu+ldwrku),  &
     &                        ldwrku)
!
!                    Generate Q in A
!                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGLQ(M,N,M,A,Lda,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in WORK(IU), copying result to
!                    WORK(IR)
!                    (CWorkspace: need   2*M*M+3*M,
!                                 prefer 2*M*M+2*M+2*M*NB)
!                    (RWorkspace: need   M)
!
                  CALL ZGEBRD(M,M,Work(iu),ldwrku,S,Rwork(ie),          &
     &                        Work(itauq),Work(itaup),Work(iwork),      &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,M,Work(iu),ldwrku,Work(ir),ldwrkr)
!
!                    Generate right bidiagonalizing vectors in WORK(IU)
!                    (CWorkspace: need   2*M*M+3*M-1,
!                                 prefer 2*M*M+2*M+(M-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',M,M,M,Work(iu),ldwrku,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in WORK(IR)
!                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',M,M,M,Work(ir),ldwrkr,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  irwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of L in WORK(IR) and computing
!                    right singular vectors of L in WORK(IU)
!                    (CWorkspace: need 2*M*M)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',M,M,M,0,S,Rwork(ie),Work(iu),ldwrku,  &
     &                        Work(ir),ldwrkr,cdum,1,Rwork(irwork),Info)
!
!                    Multiply right singular vectors of L in WORK(IU) by
!                    Q in A, storing result in VT
!                    (CWorkspace: need M*M)
!                    (RWorkspace: 0)
!
                  CALL ZGEMM('N','N',M,N,M,CONE,Work(iu),ldwrku,A,Lda,  &
     &                       CZERO,Vt,Ldvt)
!
!                    Copy left singular vectors of L to A
!                    (CWorkspace: need M*M)
!                    (RWorkspace: 0)
!
                  CALL ZLACPY('F',M,M,Work(ir),ldwrkr,A,Lda)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (CWorkspace: need 2*M, prefer M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (CWorkspace: need 2*M, prefer M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGLQ(M,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Zero out above L in A
!
                  CALL ZLASET('U',M-1,M-1,CZERO,CZERO,A(1,2),Lda)
!
!                    Bidiagonalize L in A
!                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
!                    (RWorkspace: need M)
!
                  CALL ZGEBRD(M,M,A,Lda,S,Rwork(ie),Work(itauq),        &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply right vectors bidiagonalizing L by Q in VT
!                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNMBR('P','L','C',M,N,M,A,Lda,Work(itaup),Vt,   &
     &                        Ldvt,Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors of L in A
!                    (CWorkspace: need 3*M, prefer 2*M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',M,M,M,A,Lda,Work(itauq),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  irwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in A and computing right
!                    singular vectors of A in VT
!                    (CWorkspace: 0)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',M,N,M,0,S,Rwork(ie),Vt,Ldvt,A,Lda,    &
     &                        cdum,1,Rwork(irwork),Info)
!
               ENDIF
!
            ELSEIF ( wntuas ) THEN
!
!                 Path 6t(N much larger than M, JOBU='S' or 'A',
!                         JOBVT='S')
!                 M right singular vectors to be computed in VT and
!                 M left singular vectors to be computed in U
!
               IF ( Lwork>=M*M+3*M ) THEN
!
!                    Sufficient workspace for a fast algorithm
!
                  iu = 1
                  IF ( Lwork>=wrkbl+Lda*M ) THEN
!
!                       WORK(IU) is LDA by N
!
                     ldwrku = Lda
                  ELSE
!
!                       WORK(IU) is LDA by M
!
                     ldwrku = M
                  ENDIF
                  itau = iu + ldwrku*M
                  iwork = itau + M
!
!                    Compute A=L*Q
!                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to WORK(IU), zeroing out above it
!
                  CALL ZLACPY('L',M,M,A,Lda,Work(iu),ldwrku)
                  CALL ZLASET('U',M-1,M-1,CZERO,CZERO,Work(iu+ldwrku),  &
     &                        ldwrku)
!
!                    Generate Q in A
!                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGLQ(M,N,M,A,Lda,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in WORK(IU), copying result to U
!                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
!                    (RWorkspace: need M)
!
                  CALL ZGEBRD(M,M,Work(iu),ldwrku,S,Rwork(ie),          &
     &                        Work(itauq),Work(itaup),Work(iwork),      &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,M,Work(iu),ldwrku,U,Ldu)
!
!                    Generate right bidiagonalizing vectors in WORK(IU)
!                    (CWorkspace: need   M*M+3*M-1,
!                                 prefer M*M+2*M+(M-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',M,M,M,Work(iu),ldwrku,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in U
!                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',M,M,M,U,Ldu,Work(itauq),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  irwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of L in U and computing right
!                    singular vectors of L in WORK(IU)
!                    (CWorkspace: need M*M)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',M,M,M,0,S,Rwork(ie),Work(iu),ldwrku,U,&
     &                        Ldu,cdum,1,Rwork(irwork),Info)
!
!                    Multiply right singular vectors of L in WORK(IU) by
!                    Q in A, storing result in VT
!                    (CWorkspace: need M*M)
!                    (RWorkspace: 0)
!
                  CALL ZGEMM('N','N',M,N,M,CONE,Work(iu),ldwrku,A,Lda,  &
     &                       CZERO,Vt,Ldvt)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (CWorkspace: need 2*M, prefer M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (CWorkspace: need 2*M, prefer M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGLQ(M,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to U, zeroing out above it
!
                  CALL ZLACPY('L',M,M,A,Lda,U,Ldu)
                  CALL ZLASET('U',M-1,M-1,CZERO,CZERO,U(1,2),Ldu)
                  ie = 1
                  itauq = itau
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in U
!                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
!                    (RWorkspace: need M)
!
                  CALL ZGEBRD(M,M,U,Ldu,S,Rwork(ie),Work(itauq),        &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply right bidiagonalizing vectors in U by Q
!                    in VT
!                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNMBR('P','L','C',M,N,M,U,Ldu,Work(itaup),Vt,   &
     &                        Ldvt,Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in U
!                    (CWorkspace: need 3*M, prefer 2*M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',M,M,M,U,Ldu,Work(itauq),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  irwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U and computing right
!                    singular vectors of A in VT
!                    (CWorkspace: 0)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',M,N,M,0,S,Rwork(ie),Vt,Ldvt,U,Ldu,    &
     &                        cdum,1,Rwork(irwork),Info)
!
               ENDIF
!
            ENDIF
!
         ELSEIF ( wntva ) THEN
!
            IF ( wntun ) THEN
!
!                 Path 7t(N much larger than M, JOBU='N', JOBVT='A')
!                 N right singular vectors to be computed in VT and
!                 no left singular vectors to be computed
!
               IF ( Lwork>=M*M+MAX(N+M,3*M) ) THEN
!
!                    Sufficient workspace for a fast algorithm
!
                  ir = 1
                  IF ( Lwork>=wrkbl+Lda*M ) THEN
!
!                       WORK(IR) is LDA by M
!
                     ldwrkr = Lda
                  ELSE
!
!                       WORK(IR) is M by M
!
                     ldwrkr = M
                  ENDIF
                  itau = ir + ldwrkr*M
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Copy L to WORK(IR), zeroing out above it
!
                  CALL ZLACPY('L',M,M,A,Lda,Work(ir),ldwrkr)
                  CALL ZLASET('U',M-1,M-1,CZERO,CZERO,Work(ir+ldwrkr),  &
     &                        ldwrkr)
!
!                    Generate Q in VT
!                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in WORK(IR)
!                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
!                    (RWorkspace: need M)
!
                  CALL ZGEBRD(M,M,Work(ir),ldwrkr,S,Rwork(ie),          &
     &                        Work(itauq),Work(itaup),Work(iwork),      &
     &                        Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in WORK(IR)
!                    (CWorkspace: need   M*M+3*M-1,
!                                 prefer M*M+2*M+(M-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',M,M,M,Work(ir),ldwrkr,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  irwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing right
!                    singular vectors of L in WORK(IR)
!                    (CWorkspace: need M*M)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',M,M,0,0,S,Rwork(ie),Work(ir),ldwrkr,  &
     &                        cdum,1,cdum,1,Rwork(irwork),Info)
!
!                    Multiply right singular vectors of L in WORK(IR) by
!                    Q in VT, storing result in A
!                    (CWorkspace: need M*M)
!                    (RWorkspace: 0)
!
                  CALL ZGEMM('N','N',M,N,M,CONE,Work(ir),ldwrkr,Vt,Ldvt,&
     &                       CZERO,A,Lda)
!
!                    Copy right singular vectors of A from A to VT
!
                  CALL ZLACPY('F',M,N,A,Lda,Vt,Ldvt)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (CWorkspace: need 2*M, prefer M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (CWorkspace: need M+N, prefer M+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Zero out above L in A
!
                  CALL ZLASET('U',M-1,M-1,CZERO,CZERO,A(1,2),Lda)
!
!                    Bidiagonalize L in A
!                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
!                    (RWorkspace: need M)
!
                  CALL ZGEBRD(M,M,A,Lda,S,Rwork(ie),Work(itauq),        &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply right bidiagonalizing vectors in A by Q
!                    in VT
!                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNMBR('P','L','C',M,N,M,A,Lda,Work(itaup),Vt,   &
     &                        Ldvt,Work(iwork),Lwork-iwork+1,ierr)
                  irwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing right
!                    singular vectors of A in VT
!                    (CWorkspace: 0)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',M,N,0,0,S,Rwork(ie),Vt,Ldvt,cdum,1,   &
     &                        cdum,1,Rwork(irwork),Info)
!
               ENDIF
!
            ELSEIF ( wntuo ) THEN
!
!                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
!                 N right singular vectors to be computed in VT and
!                 M left singular vectors to be overwritten on A
!
               IF ( Lwork>=2*M*M+MAX(N+M,3*M) ) THEN
!
!                    Sufficient workspace for a fast algorithm
!
                  iu = 1
                  IF ( Lwork>=wrkbl+2*Lda*M ) THEN
!
!                       WORK(IU) is LDA by M and WORK(IR) is LDA by M
!
                     ldwrku = Lda
                     ir = iu + ldwrku*M
                     ldwrkr = Lda
                  ELSEIF ( Lwork>=wrkbl+(Lda+M)*M ) THEN
!
!                       WORK(IU) is LDA by M and WORK(IR) is M by M
!
                     ldwrku = Lda
                     ir = iu + ldwrku*M
                     ldwrkr = M
                  ELSE
!
!                       WORK(IU) is M by M and WORK(IR) is M by M
!
                     ldwrku = M
                     ir = iu + ldwrku*M
                     ldwrkr = M
                  ENDIF
                  itau = ir + ldwrkr*M
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (CWorkspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to WORK(IU), zeroing out above it
!
                  CALL ZLACPY('L',M,M,A,Lda,Work(iu),ldwrku)
                  CALL ZLASET('U',M-1,M-1,CZERO,CZERO,Work(iu+ldwrku),  &
     &                        ldwrku)
                  ie = 1
                  itauq = itau
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in WORK(IU), copying result to
!                    WORK(IR)
!                    (CWorkspace: need   2*M*M+3*M,
!                                 prefer 2*M*M+2*M+2*M*NB)
!                    (RWorkspace: need   M)
!
                  CALL ZGEBRD(M,M,Work(iu),ldwrku,S,Rwork(ie),          &
     &                        Work(itauq),Work(itaup),Work(iwork),      &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,M,Work(iu),ldwrku,Work(ir),ldwrkr)
!
!                    Generate right bidiagonalizing vectors in WORK(IU)
!                    (CWorkspace: need   2*M*M+3*M-1,
!                                 prefer 2*M*M+2*M+(M-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',M,M,M,Work(iu),ldwrku,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in WORK(IR)
!                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',M,M,M,Work(ir),ldwrkr,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  irwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of L in WORK(IR) and computing
!                    right singular vectors of L in WORK(IU)
!                    (CWorkspace: need 2*M*M)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',M,M,M,0,S,Rwork(ie),Work(iu),ldwrku,  &
     &                        Work(ir),ldwrkr,cdum,1,Rwork(irwork),Info)
!
!                    Multiply right singular vectors of L in WORK(IU) by
!                    Q in VT, storing result in A
!                    (CWorkspace: need M*M)
!                    (RWorkspace: 0)
!
                  CALL ZGEMM('N','N',M,N,M,CONE,Work(iu),ldwrku,Vt,Ldvt,&
     &                       CZERO,A,Lda)
!
!                    Copy right singular vectors of A from A to VT
!
                  CALL ZLACPY('F',M,N,A,Lda,Vt,Ldvt)
!
!                    Copy left singular vectors of A from WORK(IR) to A
!
                  CALL ZLACPY('F',M,M,Work(ir),ldwrkr,A,Lda)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (CWorkspace: need 2*M, prefer M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (CWorkspace: need M+N, prefer M+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
                  ie = 1
                  itauq = itau
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Zero out above L in A
!
                  CALL ZLASET('U',M-1,M-1,CZERO,CZERO,A(1,2),Lda)
!
!                    Bidiagonalize L in A
!                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
!                    (RWorkspace: need M)
!
                  CALL ZGEBRD(M,M,A,Lda,S,Rwork(ie),Work(itauq),        &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply right bidiagonalizing vectors in A by Q
!                    in VT
!                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNMBR('P','L','C',M,N,M,A,Lda,Work(itaup),Vt,   &
     &                        Ldvt,Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in A
!                    (CWorkspace: need 3*M, prefer 2*M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',M,M,M,A,Lda,Work(itauq),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  irwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in A and computing right
!                    singular vectors of A in VT
!                    (CWorkspace: 0)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',M,N,M,0,S,Rwork(ie),Vt,Ldvt,A,Lda,    &
     &                        cdum,1,Rwork(irwork),Info)
!
               ENDIF
!
            ELSEIF ( wntuas ) THEN
!
!                 Path 9t(N much larger than M, JOBU='S' or 'A',
!                         JOBVT='A')
!                 N right singular vectors to be computed in VT and
!                 M left singular vectors to be computed in U
!
               IF ( Lwork>=M*M+MAX(N+M,3*M) ) THEN
!
!                    Sufficient workspace for a fast algorithm
!
                  iu = 1
                  IF ( Lwork>=wrkbl+Lda*M ) THEN
!
!                       WORK(IU) is LDA by M
!
                     ldwrku = Lda
                  ELSE
!
!                       WORK(IU) is M by M
!
                     ldwrku = M
                  ENDIF
                  itau = iu + ldwrku*M
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to WORK(IU), zeroing out above it
!
                  CALL ZLACPY('L',M,M,A,Lda,Work(iu),ldwrku)
                  CALL ZLASET('U',M-1,M-1,CZERO,CZERO,Work(iu+ldwrku),  &
     &                        ldwrku)
                  ie = 1
                  itauq = itau
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in WORK(IU), copying result to U
!                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
!                    (RWorkspace: need M)
!
                  CALL ZGEBRD(M,M,Work(iu),ldwrku,S,Rwork(ie),          &
     &                        Work(itauq),Work(itaup),Work(iwork),      &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('L',M,M,Work(iu),ldwrku,U,Ldu)
!
!                    Generate right bidiagonalizing vectors in WORK(IU)
!                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('P',M,M,M,Work(iu),ldwrku,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in U
!                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',M,M,M,U,Ldu,Work(itauq),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  irwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of L in U and computing right
!                    singular vectors of L in WORK(IU)
!                    (CWorkspace: need M*M)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',M,M,M,0,S,Rwork(ie),Work(iu),ldwrku,U,&
     &                        Ldu,cdum,1,Rwork(irwork),Info)
!
!                    Multiply right singular vectors of L in WORK(IU) by
!                    Q in VT, storing result in A
!                    (CWorkspace: need M*M)
!                    (RWorkspace: 0)
!
                  CALL ZGEMM('N','N',M,N,M,CONE,Work(iu),ldwrku,Vt,Ldvt,&
     &                       CZERO,A,Lda)
!
!                    Copy right singular vectors of A from A to VT
!
                  CALL ZLACPY('F',M,N,A,Lda,Vt,Ldvt)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (CWorkspace: need 2*M, prefer M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL ZLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (CWorkspace: need M+N, prefer M+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to U, zeroing out above it
!
                  CALL ZLACPY('L',M,M,A,Lda,U,Ldu)
                  CALL ZLASET('U',M-1,M-1,CZERO,CZERO,U(1,2),Ldu)
                  ie = 1
                  itauq = itau
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in U
!                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
!                    (RWorkspace: need M)
!
                  CALL ZGEBRD(M,M,U,Ldu,S,Rwork(ie),Work(itauq),        &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply right bidiagonalizing vectors in U by Q
!                    in VT
!                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNMBR('P','L','C',M,N,M,U,Ldu,Work(itaup),Vt,   &
     &                        Ldvt,Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in U
!                    (CWorkspace: need 3*M, prefer 2*M+M*NB)
!                    (RWorkspace: 0)
!
                  CALL ZUNGBR('Q',M,M,M,U,Ldu,Work(itauq),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  irwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U and computing right
!                    singular vectors of A in VT
!                    (CWorkspace: 0)
!                    (RWorkspace: need BDSPAC)
!
                  CALL ZBDSQR('U',M,N,M,0,S,Rwork(ie),Vt,Ldvt,U,Ldu,    &
     &                        cdum,1,Rwork(irwork),Info)
!
               ENDIF
!
            ENDIF
!
         ENDIF
!
      ELSE
!
!           N .LT. MNTHR
!
!           Path 10t(N greater than M, but not much larger)
!           Reduce to bidiagonal form without LQ decomposition
!
         ie = 1
         itauq = 1
         itaup = itauq + M
         iwork = itaup + M
!
!           Bidiagonalize A
!           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
!           (RWorkspace: M)
!
         CALL ZGEBRD(M,N,A,Lda,S,Rwork(ie),Work(itauq),Work(itaup),     &
     &               Work(iwork),Lwork-iwork+1,ierr)
         IF ( wntuas ) THEN
!
!              If left singular vectors desired in U, copy result to U
!              and generate left bidiagonalizing vectors in U
!              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB)
!              (RWorkspace: 0)
!
            CALL ZLACPY('L',M,M,A,Lda,U,Ldu)
            CALL ZUNGBR('Q',M,M,N,U,Ldu,Work(itauq),Work(iwork),        &
     &                  Lwork-iwork+1,ierr)
         ENDIF
         IF ( wntvas ) THEN
!
!              If right singular vectors desired in VT, copy result to
!              VT and generate right bidiagonalizing vectors in VT
!              (CWorkspace: need 2*M+NRVT, prefer 2*M+NRVT*NB)
!              (RWorkspace: 0)
!
            CALL ZLACPY('U',M,N,A,Lda,Vt,Ldvt)
            IF ( wntva ) nrvt = N
            IF ( wntvs ) nrvt = M
            CALL ZUNGBR('P',nrvt,N,M,Vt,Ldvt,Work(itaup),Work(iwork),   &
     &                  Lwork-iwork+1,ierr)
         ENDIF
!
!              If left singular vectors desired in A, generate left
!              bidiagonalizing vectors in A
!              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB)
!              (RWorkspace: 0)
!
         IF ( wntuo ) CALL ZUNGBR('Q',M,M,N,A,Lda,Work(itauq),          &
     &                            Work(iwork),Lwork-iwork+1,ierr)
!
!              If right singular vectors desired in A, generate right
!              bidiagonalizing vectors in A
!              (CWorkspace: need 3*M, prefer 2*M+M*NB)
!              (RWorkspace: 0)
!
         IF ( wntvo ) CALL ZUNGBR('P',M,N,M,A,Lda,Work(itaup),          &
     &                            Work(iwork),Lwork-iwork+1,ierr)
         irwork = ie + M
         IF ( wntuas .OR. wntuo ) nru = M
         IF ( wntun ) nru = 0
         IF ( wntvas .OR. wntvo ) ncvt = N
         IF ( wntvn ) ncvt = 0
         IF ( (.NOT.wntuo) .AND. (.NOT.wntvo) ) THEN
!
!              Perform bidiagonal QR iteration, if desired, computing
!              left singular vectors in U and computing right singular
!              vectors in VT
!              (CWorkspace: 0)
!              (RWorkspace: need BDSPAC)
!
            CALL ZBDSQR('L',M,ncvt,nru,0,S,Rwork(ie),Vt,Ldvt,U,Ldu,cdum,&
     &                  1,Rwork(irwork),Info)
         ELSEIF ( (.NOT.wntuo) .AND. wntvo ) THEN
!
!              Perform bidiagonal QR iteration, if desired, computing
!              left singular vectors in U and computing right singular
!              vectors in A
!              (CWorkspace: 0)
!              (RWorkspace: need BDSPAC)
!
            CALL ZBDSQR('L',M,ncvt,nru,0,S,Rwork(ie),A,Lda,U,Ldu,cdum,1,&
     &                  Rwork(irwork),Info)
         ELSE
!
!              Perform bidiagonal QR iteration, if desired, computing
!              left singular vectors in A and computing right singular
!              vectors in VT
!              (CWorkspace: 0)
!              (RWorkspace: need BDSPAC)
!
            CALL ZBDSQR('L',M,ncvt,nru,0,S,Rwork(ie),Vt,Ldvt,A,Lda,cdum,&
     &                  1,Rwork(irwork),Info)
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
         IF ( Info/=0 .AND. anrm>bignum )                               &
     &        CALL DLASCL('G',0,0,bignum,anrm,minmn-1,1,Rwork(ie),minmn,&
     &        ierr)
         IF ( anrm<smlnum ) CALL DLASCL('G',0,0,smlnum,anrm,minmn,1,S,  &
     &                                  minmn,ierr)
         IF ( Info/=0 .AND. anrm<smlnum )                               &
     &        CALL DLASCL('G',0,0,smlnum,anrm,minmn-1,1,Rwork(ie),minmn,&
     &        ierr)
      ENDIF
!
!     Return optimal workspace in WORK(1)
!
      Work(1) = maxwrk
!
!
!     End of ZGESVD
!
      END SUBROUTINE ZGESVD
