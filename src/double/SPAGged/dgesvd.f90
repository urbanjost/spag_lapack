!*==dgesvd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> DGESVD computes the singular value decomposition (SVD) for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGESVD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBU, JOBVT
!       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
!       ..
!       .. Array Arguments ..
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
!> DGESVD computes the singular value decomposition (SVD) of a real
!> M-by-N matrix A, optionally computing the left and/or right singular
!> vectors. The SVD is written
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
!> Note that the routine returns V**T, not V.
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
!>          V**T:
!>          = 'A':  all N rows of V**T are returned in the array VT;
!>          = 'S':  the first min(m,n) rows of V**T (the right singular
!>                  vectors) are returned in the array VT;
!>          = 'O':  the first min(m,n) rows of V**T (the right singular
!>                  vectors) are overwritten on the array A;
!>          = 'N':  no rows of V**T (no right singular vectors) are
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit,
!>          if JOBU = 'O',  A is overwritten with the first min(m,n)
!>                          columns of U (the left singular vectors,
!>                          stored columnwise);
!>          if JOBVT = 'O', A is overwritten with the first min(m,n)
!>                          rows of V**T (the right singular vectors,
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
!>          U is DOUBLE PRECISION array, dimension (LDU,UCOL)
!>          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
!>          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
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
!>          VT is DOUBLE PRECISION array, dimension (LDVT,N)
!>          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
!>          V**T;
!>          if JOBVT = 'S', VT contains the first min(m,n) rows of
!>          V**T (the right singular vectors, stored rowwise);
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
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!>          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
!>          superdiagonal elements of an upper bidiagonal matrix B
!>          whose diagonal is in S (not necessarily sorted). B
!>          satisfies A = U * B * VT, so it has the same singular values
!>          as A, and singular vectors related by U and VT.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
!>             - PATH 1  (M much larger than N, JOBU='N')
!>             - PATH 1t (N much larger than M, JOBVT='N')
!>          LWORK >= MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) for the other paths
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
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if DBDSQR did not converge, INFO specifies how many
!>                superdiagonals of an intermediate bidiagonal form B
!>                did not converge to zero. See the description of WORK
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
!> \ingroup doubleGEsing
!
!  =====================================================================
      SUBROUTINE DGESVD(Jobu,Jobvt,M,N,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,&
     &                  Info)
      USE F77KINDS                        
      USE S_DBDSQR
      USE S_DGEBRD
      USE S_DGELQF
      USE S_DGEMM
      USE S_DGEQRF
      USE S_DLACPY
      USE S_DLAMCH
      USE S_DLANGE
      USE S_DLASCL
      USE S_DLASET
      USE S_DORGBR
      USE S_DORGLQ
      USE S_DORGQR
      USE S_DORMBR
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DGESVD233
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobu
      CHARACTER :: Jobvt
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anrm , bignum , eps , smlnum
      INTEGER :: bdspac , blk , chunk , i , ie , ierr , ir , iscl ,     &
     &           itau , itaup , itauq , iu , iwork , ldwrkr , ldwrku ,  &
     &           lwork_dgebrd , lwork_dgelqf , lwork_dgeqrf ,           &
     &           lwork_dorgbr_p , lwork_dorgbr_q , lwork_dorglq_m ,     &
     &           lwork_dorglq_n , lwork_dorgqr_m , lwork_dorgqr_n ,     &
     &           maxwrk , minmn , minwrk , mnthr , ncu , ncvt , nru ,   &
     &           nrvt , wrkbl
      REAL(R8KIND) , DIMENSION(1) :: dum
      LOGICAL :: lquery , wntua , wntuas , wntun , wntuo , wntus ,      &
     &           wntva , wntvas , wntvn , wntvo , wntvs
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
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.)
!
      IF ( Info==0 ) THEN
         minwrk = 1
         maxwrk = 1
         IF ( M>=N .AND. minmn>0 ) THEN
!
!           Compute space needed for DBDSQR
!
            mnthr = ILAENV(6,'DGESVD',Jobu//Jobvt,M,N,0,0)
            bdspac = 5*N
!           Compute space needed for DGEQRF
            CALL DGEQRF(M,N,A,Lda,dum(1),dum(1),-1,ierr)
            lwork_dgeqrf = INT(dum(1))
!           Compute space needed for DORGQR
            CALL DORGQR(M,N,N,A,Lda,dum(1),dum(1),-1,ierr)
            lwork_dorgqr_n = INT(dum(1))
            CALL DORGQR(M,M,N,A,Lda,dum(1),dum(1),-1,ierr)
            lwork_dorgqr_m = INT(dum(1))
!           Compute space needed for DGEBRD
            CALL DGEBRD(N,N,A,Lda,S,dum(1),dum(1),dum(1),dum(1),-1,ierr)
            lwork_dgebrd = INT(dum(1))
!           Compute space needed for DORGBR P
            CALL DORGBR('P',N,N,N,A,Lda,dum(1),dum(1),-1,ierr)
            lwork_dorgbr_p = INT(dum(1))
!           Compute space needed for DORGBR Q
            CALL DORGBR('Q',N,N,N,A,Lda,dum(1),dum(1),-1,ierr)
            lwork_dorgbr_q = INT(dum(1))
!
            IF ( M<mnthr ) THEN
!
!              Path 10 (M at least N, but not much larger)
!
               CALL DGEBRD(M,N,A,Lda,S,dum(1),dum(1),dum(1),dum(1),-1,  &
     &                     ierr)
               lwork_dgebrd = INT(dum(1))
               maxwrk = 3*N + lwork_dgebrd
               IF ( wntus .OR. wntuo ) THEN
                  CALL DORGBR('Q',M,N,N,A,Lda,dum(1),dum(1),-1,ierr)
                  lwork_dorgbr_q = INT(dum(1))
                  maxwrk = MAX(maxwrk,3*N+lwork_dorgbr_q)
               ENDIF
               IF ( wntua ) THEN
                  CALL DORGBR('Q',M,M,N,A,Lda,dum(1),dum(1),-1,ierr)
                  lwork_dorgbr_q = INT(dum(1))
                  maxwrk = MAX(maxwrk,3*N+lwork_dorgbr_q)
               ENDIF
               IF ( .NOT.wntvn ) maxwrk = MAX(maxwrk,3*N+lwork_dorgbr_p)
               maxwrk = MAX(maxwrk,bdspac)
               minwrk = MAX(3*N+M,bdspac)
            ELSEIF ( wntun ) THEN
!
!                 Path 1 (M much larger than N, JOBU='N')
!
               maxwrk = N + lwork_dgeqrf
               maxwrk = MAX(maxwrk,3*N+lwork_dgebrd)
               IF ( wntvo .OR. wntvas )                                 &
     &              maxwrk = MAX(maxwrk,3*N+lwork_dorgbr_p)
               maxwrk = MAX(maxwrk,bdspac)
               minwrk = MAX(4*N,bdspac)
            ELSEIF ( wntuo .AND. wntvn ) THEN
!
!                 Path 2 (M much larger than N, JOBU='O', JOBVT='N')
!
               wrkbl = N + lwork_dgeqrf
               wrkbl = MAX(wrkbl,N+lwork_dorgqr_n)
               wrkbl = MAX(wrkbl,3*N+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = MAX(N*N+wrkbl,N*N+M*N+N)
               minwrk = MAX(3*N+M,bdspac)
            ELSEIF ( wntuo .AND. wntvas ) THEN
!
!                 Path 3 (M much larger than N, JOBU='O', JOBVT='S' or
!                 'A')
!
               wrkbl = N + lwork_dgeqrf
               wrkbl = MAX(wrkbl,N+lwork_dorgqr_n)
               wrkbl = MAX(wrkbl,3*N+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = MAX(N*N+wrkbl,N*N+M*N+N)
               minwrk = MAX(3*N+M,bdspac)
            ELSEIF ( wntus .AND. wntvn ) THEN
!
!                 Path 4 (M much larger than N, JOBU='S', JOBVT='N')
!
               wrkbl = N + lwork_dgeqrf
               wrkbl = MAX(wrkbl,N+lwork_dorgqr_n)
               wrkbl = MAX(wrkbl,3*N+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = N*N + wrkbl
               minwrk = MAX(3*N+M,bdspac)
            ELSEIF ( wntus .AND. wntvo ) THEN
!
!                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
!
               wrkbl = N + lwork_dgeqrf
               wrkbl = MAX(wrkbl,N+lwork_dorgqr_n)
               wrkbl = MAX(wrkbl,3*N+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = 2*N*N + wrkbl
               minwrk = MAX(3*N+M,bdspac)
            ELSEIF ( wntus .AND. wntvas ) THEN
!
!                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' or
!                 'A')
!
               wrkbl = N + lwork_dgeqrf
               wrkbl = MAX(wrkbl,N+lwork_dorgqr_n)
               wrkbl = MAX(wrkbl,3*N+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = N*N + wrkbl
               minwrk = MAX(3*N+M,bdspac)
            ELSEIF ( wntua .AND. wntvn ) THEN
!
!                 Path 7 (M much larger than N, JOBU='A', JOBVT='N')
!
               wrkbl = N + lwork_dgeqrf
               wrkbl = MAX(wrkbl,N+lwork_dorgqr_m)
               wrkbl = MAX(wrkbl,3*N+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = N*N + wrkbl
               minwrk = MAX(3*N+M,bdspac)
            ELSEIF ( wntua .AND. wntvo ) THEN
!
!                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
!
               wrkbl = N + lwork_dgeqrf
               wrkbl = MAX(wrkbl,N+lwork_dorgqr_m)
               wrkbl = MAX(wrkbl,3*N+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = 2*N*N + wrkbl
               minwrk = MAX(3*N+M,bdspac)
            ELSEIF ( wntua .AND. wntvas ) THEN
!
!                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' or
!                 'A')
!
               wrkbl = N + lwork_dgeqrf
               wrkbl = MAX(wrkbl,N+lwork_dorgqr_m)
               wrkbl = MAX(wrkbl,3*N+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,3*N+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = N*N + wrkbl
               minwrk = MAX(3*N+M,bdspac)
            ENDIF
         ELSEIF ( minmn>0 ) THEN
!
!           Compute space needed for DBDSQR
!
            mnthr = ILAENV(6,'DGESVD',Jobu//Jobvt,M,N,0,0)
            bdspac = 5*M
!           Compute space needed for DGELQF
            CALL DGELQF(M,N,A,Lda,dum(1),dum(1),-1,ierr)
            lwork_dgelqf = INT(dum(1))
!           Compute space needed for DORGLQ
            CALL DORGLQ(N,N,M,dum(1),N,dum(1),dum(1),-1,ierr)
            lwork_dorglq_n = INT(dum(1))
            CALL DORGLQ(M,N,M,A,Lda,dum(1),dum(1),-1,ierr)
            lwork_dorglq_m = INT(dum(1))
!           Compute space needed for DGEBRD
            CALL DGEBRD(M,M,A,Lda,S,dum(1),dum(1),dum(1),dum(1),-1,ierr)
            lwork_dgebrd = INT(dum(1))
!            Compute space needed for DORGBR P
            CALL DORGBR('P',M,M,M,A,N,dum(1),dum(1),-1,ierr)
            lwork_dorgbr_p = INT(dum(1))
!           Compute space needed for DORGBR Q
            CALL DORGBR('Q',M,M,M,A,N,dum(1),dum(1),-1,ierr)
            lwork_dorgbr_q = INT(dum(1))
            IF ( N<mnthr ) THEN
!
!              Path 10t(N greater than M, but not much larger)
!
               CALL DGEBRD(M,N,A,Lda,S,dum(1),dum(1),dum(1),dum(1),-1,  &
     &                     ierr)
               lwork_dgebrd = INT(dum(1))
               maxwrk = 3*M + lwork_dgebrd
               IF ( wntvs .OR. wntvo ) THEN
!                Compute space needed for DORGBR P
                  CALL DORGBR('P',M,N,M,A,N,dum(1),dum(1),-1,ierr)
                  lwork_dorgbr_p = INT(dum(1))
                  maxwrk = MAX(maxwrk,3*M+lwork_dorgbr_p)
               ENDIF
               IF ( wntva ) THEN
                  CALL DORGBR('P',N,N,M,A,N,dum(1),dum(1),-1,ierr)
                  lwork_dorgbr_p = INT(dum(1))
                  maxwrk = MAX(maxwrk,3*M+lwork_dorgbr_p)
               ENDIF
               IF ( .NOT.wntun ) maxwrk = MAX(maxwrk,3*M+lwork_dorgbr_q)
               maxwrk = MAX(maxwrk,bdspac)
               minwrk = MAX(3*M+N,bdspac)
            ELSEIF ( wntvn ) THEN
!
!                 Path 1t(N much larger than M, JOBVT='N')
!
               maxwrk = M + lwork_dgelqf
               maxwrk = MAX(maxwrk,3*M+lwork_dgebrd)
               IF ( wntuo .OR. wntuas )                                 &
     &              maxwrk = MAX(maxwrk,3*M+lwork_dorgbr_q)
               maxwrk = MAX(maxwrk,bdspac)
               minwrk = MAX(4*M,bdspac)
            ELSEIF ( wntvo .AND. wntun ) THEN
!
!                 Path 2t(N much larger than M, JOBU='N', JOBVT='O')
!
               wrkbl = M + lwork_dgelqf
               wrkbl = MAX(wrkbl,M+lwork_dorglq_m)
               wrkbl = MAX(wrkbl,3*M+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = MAX(M*M+wrkbl,M*M+M*N+M)
               minwrk = MAX(3*M+N,bdspac)
            ELSEIF ( wntvo .AND. wntuas ) THEN
!
!                 Path 3t(N much larger than M, JOBU='S' or 'A',
!                 JOBVT='O')
!
               wrkbl = M + lwork_dgelqf
               wrkbl = MAX(wrkbl,M+lwork_dorglq_m)
               wrkbl = MAX(wrkbl,3*M+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = MAX(M*M+wrkbl,M*M+M*N+M)
               minwrk = MAX(3*M+N,bdspac)
            ELSEIF ( wntvs .AND. wntun ) THEN
!
!                 Path 4t(N much larger than M, JOBU='N', JOBVT='S')
!
               wrkbl = M + lwork_dgelqf
               wrkbl = MAX(wrkbl,M+lwork_dorglq_m)
               wrkbl = MAX(wrkbl,3*M+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = M*M + wrkbl
               minwrk = MAX(3*M+N,bdspac)
            ELSEIF ( wntvs .AND. wntuo ) THEN
!
!                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
!
               wrkbl = M + lwork_dgelqf
               wrkbl = MAX(wrkbl,M+lwork_dorglq_m)
               wrkbl = MAX(wrkbl,3*M+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = 2*M*M + wrkbl
               minwrk = MAX(3*M+N,bdspac)
            ELSEIF ( wntvs .AND. wntuas ) THEN
!
!                 Path 6t(N much larger than M, JOBU='S' or 'A',
!                 JOBVT='S')
!
               wrkbl = M + lwork_dgelqf
               wrkbl = MAX(wrkbl,M+lwork_dorglq_m)
               wrkbl = MAX(wrkbl,3*M+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = M*M + wrkbl
               minwrk = MAX(3*M+N,bdspac)
            ELSEIF ( wntva .AND. wntun ) THEN
!
!                 Path 7t(N much larger than M, JOBU='N', JOBVT='A')
!
               wrkbl = M + lwork_dgelqf
               wrkbl = MAX(wrkbl,M+lwork_dorglq_n)
               wrkbl = MAX(wrkbl,3*M+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = M*M + wrkbl
               minwrk = MAX(3*M+N,bdspac)
            ELSEIF ( wntva .AND. wntuo ) THEN
!
!                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
!
               wrkbl = M + lwork_dgelqf
               wrkbl = MAX(wrkbl,M+lwork_dorglq_n)
               wrkbl = MAX(wrkbl,3*M+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = 2*M*M + wrkbl
               minwrk = MAX(3*M+N,bdspac)
            ELSEIF ( wntva .AND. wntuas ) THEN
!
!                 Path 9t(N much larger than M, JOBU='S' or 'A',
!                 JOBVT='A')
!
               wrkbl = M + lwork_dgelqf
               wrkbl = MAX(wrkbl,M+lwork_dorglq_n)
               wrkbl = MAX(wrkbl,3*M+lwork_dgebrd)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_p)
               wrkbl = MAX(wrkbl,3*M+lwork_dorgbr_q)
               wrkbl = MAX(wrkbl,bdspac)
               maxwrk = M*M + wrkbl
               minwrk = MAX(3*M+N,bdspac)
            ENDIF
         ENDIF
         maxwrk = MAX(maxwrk,minwrk)
         Work(1) = maxwrk
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -13
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGESVD',-Info)
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
!           Path 10 (M at least N, but not much larger)
!           Reduce to bidiagonal form without QR decomposition
!
            ie = 1
            itauq = ie + N
            itaup = itauq + N
            iwork = itaup + N
!
!           Bidiagonalize A
!           (Workspace: need 3*N + M, prefer 3*N + (M + N)*NB)
!
            CALL DGEBRD(M,N,A,Lda,S,Work(ie),Work(itauq),Work(itaup),   &
     &                  Work(iwork),Lwork-iwork+1,ierr)
            IF ( wntuas ) THEN
!
!              If left singular vectors desired in U, copy result to U
!              and generate left bidiagonalizing vectors in U
!              (Workspace: need 3*N + NCU, prefer 3*N + NCU*NB)
!
               CALL DLACPY('L',M,N,A,Lda,U,Ldu)
               IF ( wntus ) ncu = N
               IF ( wntua ) ncu = M
               CALL DORGBR('Q',M,ncu,N,U,Ldu,Work(itauq),Work(iwork),   &
     &                     Lwork-iwork+1,ierr)
            ENDIF
            IF ( wntvas ) THEN
!
!              If right singular vectors desired in VT, copy result to
!              VT and generate right bidiagonalizing vectors in VT
!              (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
!
               CALL DLACPY('U',N,N,A,Lda,Vt,Ldvt)
               CALL DORGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),   &
     &                     Lwork-iwork+1,ierr)
            ENDIF
!
!              If left singular vectors desired in A, generate left
!              bidiagonalizing vectors in A
!              (Workspace: need 4*N, prefer 3*N + N*NB)
!
            IF ( wntuo ) CALL DORGBR('Q',M,N,N,A,Lda,Work(itauq),       &
     &                               Work(iwork),Lwork-iwork+1,ierr)
!
!              If right singular vectors desired in A, generate right
!              bidiagonalizing vectors in A
!              (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
!
            IF ( wntvo ) CALL DORGBR('P',N,N,N,A,Lda,Work(itaup),       &
     &                               Work(iwork),Lwork-iwork+1,ierr)
            iwork = ie + N
            IF ( wntuas .OR. wntuo ) nru = M
            IF ( wntun ) nru = 0
            IF ( wntvas .OR. wntvo ) ncvt = N
            IF ( wntvn ) ncvt = 0
            IF ( (.NOT.wntuo) .AND. (.NOT.wntvo) ) THEN
!
!              Perform bidiagonal QR iteration, if desired, computing
!              left singular vectors in U and computing right singular
!              vectors in VT
!              (Workspace: need BDSPAC)
!
               CALL DBDSQR('U',N,ncvt,nru,0,S,Work(ie),Vt,Ldvt,U,Ldu,   &
     &                     dum,1,Work(iwork),Info)
            ELSEIF ( (.NOT.wntuo) .AND. wntvo ) THEN
!
!              Perform bidiagonal QR iteration, if desired, computing
!              left singular vectors in U and computing right singular
!              vectors in A
!              (Workspace: need BDSPAC)
!
               CALL DBDSQR('U',N,ncvt,nru,0,S,Work(ie),A,Lda,U,Ldu,dum, &
     &                     1,Work(iwork),Info)
            ELSE
!
!              Perform bidiagonal QR iteration, if desired, computing
!              left singular vectors in A and computing right singular
!              vectors in VT
!              (Workspace: need BDSPAC)
!
               CALL DBDSQR('U',N,ncvt,nru,0,S,Work(ie),Vt,Ldvt,A,Lda,   &
     &                     dum,1,Work(iwork),Info)
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
!              (Workspace: need 2*N, prefer N + N*NB)
!
            CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),Lwork-iwork+1, &
     &                  ierr)
!
!              Zero out below R
!
            IF ( N>1 ) CALL DLASET('L',N-1,N-1,ZERO,ZERO,A(2,1),Lda)
            ie = 1
            itauq = ie + N
            itaup = itauq + N
            iwork = itaup + N
!
!              Bidiagonalize R in A
!              (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!
            CALL DGEBRD(N,N,A,Lda,S,Work(ie),Work(itauq),Work(itaup),   &
     &                  Work(iwork),Lwork-iwork+1,ierr)
            ncvt = 0
            IF ( wntvo .OR. wntvas ) THEN
!
!                 If right singular vectors desired, generate P'.
!                 (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
!
               CALL DORGBR('P',N,N,N,A,Lda,Work(itaup),Work(iwork),     &
     &                     Lwork-iwork+1,ierr)
               ncvt = N
            ENDIF
            iwork = ie + N
!
!              Perform bidiagonal QR iteration, computing right
!              singular vectors of A in A if desired
!              (Workspace: need BDSPAC)
!
            CALL DBDSQR('U',N,ncvt,0,0,S,Work(ie),A,Lda,dum,1,dum,1,    &
     &                  Work(iwork),Info)
!
!              If right singular vectors desired in VT, copy them there
!
            IF ( wntvas ) CALL DLACPY('F',N,N,A,Lda,Vt,Ldvt)
!
         ELSEIF ( wntuo .AND. wntvn ) THEN
!
!              Path 2 (M much larger than N, JOBU='O', JOBVT='N')
!              N left singular vectors to be overwritten on A and
!              no right singular vectors to be computed
!
            IF ( Lwork>=N*N+MAX(4*N,bdspac) ) THEN
!
!                 Sufficient workspace for a fast algorithm
!
               ir = 1
               IF ( Lwork>=MAX(wrkbl,Lda*N+N)+Lda*N ) THEN
!
!                    WORK(IU) is LDA by N, WORK(IR) is LDA by N
!
                  ldwrku = Lda
                  ldwrkr = Lda
               ELSEIF ( Lwork>=MAX(wrkbl,Lda*N+N)+N*N ) THEN
!
!                    WORK(IU) is LDA by N, WORK(IR) is N by N
!
                  ldwrku = Lda
                  ldwrkr = N
               ELSE
!
!                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N
!
                  ldwrku = (Lwork-N*N-N)/N
                  ldwrkr = N
               ENDIF
               itau = ir + ldwrkr*N
               iwork = itau + N
!
!                 Compute A=Q*R
!                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!
               CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),            &
     &                     Lwork-iwork+1,ierr)
!
!                 Copy R to WORK(IR) and zero out below it
!
               CALL DLACPY('U',N,N,A,Lda,Work(ir),ldwrkr)
               CALL DLASET('L',N-1,N-1,ZERO,ZERO,Work(ir+1),ldwrkr)
!
!                 Generate Q in A
!                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!
               CALL DORGQR(M,N,N,A,Lda,Work(itau),Work(iwork),          &
     &                     Lwork-iwork+1,ierr)
               ie = itau
               itauq = ie + N
               itaup = itauq + N
               iwork = itaup + N
!
!                 Bidiagonalize R in WORK(IR)
!                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
!
               CALL DGEBRD(N,N,Work(ir),ldwrkr,S,Work(ie),Work(itauq),  &
     &                     Work(itaup),Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate left vectors bidiagonalizing R
!                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
!
               CALL DORGBR('Q',N,N,N,Work(ir),ldwrkr,Work(itauq),       &
     &                     Work(iwork),Lwork-iwork+1,ierr)
               iwork = ie + N
!
!                 Perform bidiagonal QR iteration, computing left
!                 singular vectors of R in WORK(IR)
!                 (Workspace: need N*N + BDSPAC)
!
               CALL DBDSQR('U',N,0,N,0,S,Work(ie),dum,1,Work(ir),ldwrkr,&
     &                     dum,1,Work(iwork),Info)
               iu = ie + N
!
!                 Multiply Q in A by left singular vectors of R in
!                 WORK(IR), storing result in WORK(IU) and copying to A
!                 (Workspace: need N*N + 2*N, prefer N*N + M*N + N)
!
               DO i = 1 , M , ldwrku
                  chunk = MIN(M-i+1,ldwrku)
                  CALL DGEMM('N','N',chunk,N,N,ONE,A(i,1),Lda,Work(ir), &
     &                       ldwrkr,ZERO,Work(iu),ldwrku)
                  CALL DLACPY('F',chunk,N,Work(iu),ldwrku,A(i,1),Lda)
               ENDDO
!
            ELSE
!
!                 Insufficient workspace for a fast algorithm
!
               ie = 1
               itauq = ie + N
               itaup = itauq + N
               iwork = itaup + N
!
!                 Bidiagonalize A
!                 (Workspace: need 3*N + M, prefer 3*N + (M + N)*NB)
!
               CALL DGEBRD(M,N,A,Lda,S,Work(ie),Work(itauq),Work(itaup),&
     &                     Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate left vectors bidiagonalizing A
!                 (Workspace: need 4*N, prefer 3*N + N*NB)
!
               CALL DORGBR('Q',M,N,N,A,Lda,Work(itauq),Work(iwork),     &
     &                     Lwork-iwork+1,ierr)
               iwork = ie + N
!
!                 Perform bidiagonal QR iteration, computing left
!                 singular vectors of A in A
!                 (Workspace: need BDSPAC)
!
               CALL DBDSQR('U',N,0,M,0,S,Work(ie),dum,1,A,Lda,dum,1,    &
     &                     Work(iwork),Info)
!
            ENDIF
!
         ELSEIF ( wntuo .AND. wntvas ) THEN
!
!              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A')
!              N left singular vectors to be overwritten on A and
!              N right singular vectors to be computed in VT
!
            IF ( Lwork>=N*N+MAX(4*N,bdspac) ) THEN
!
!                 Sufficient workspace for a fast algorithm
!
               ir = 1
               IF ( Lwork>=MAX(wrkbl,Lda*N+N)+Lda*N ) THEN
!
!                    WORK(IU) is LDA by N and WORK(IR) is LDA by N
!
                  ldwrku = Lda
                  ldwrkr = Lda
               ELSEIF ( Lwork>=MAX(wrkbl,Lda*N+N)+N*N ) THEN
!
!                    WORK(IU) is LDA by N and WORK(IR) is N by N
!
                  ldwrku = Lda
                  ldwrkr = N
               ELSE
!
!                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N
!
                  ldwrku = (Lwork-N*N-N)/N
                  ldwrkr = N
               ENDIF
               itau = ir + ldwrkr*N
               iwork = itau + N
!
!                 Compute A=Q*R
!                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!
               CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),            &
     &                     Lwork-iwork+1,ierr)
!
!                 Copy R to VT, zeroing out below it
!
               CALL DLACPY('U',N,N,A,Lda,Vt,Ldvt)
               IF ( N>1 ) CALL DLASET('L',N-1,N-1,ZERO,ZERO,Vt(2,1),    &
     &                                Ldvt)
!
!                 Generate Q in A
!                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!
               CALL DORGQR(M,N,N,A,Lda,Work(itau),Work(iwork),          &
     &                     Lwork-iwork+1,ierr)
               ie = itau
               itauq = ie + N
               itaup = itauq + N
               iwork = itaup + N
!
!                 Bidiagonalize R in VT, copying result to WORK(IR)
!                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
!
               CALL DGEBRD(N,N,Vt,Ldvt,S,Work(ie),Work(itauq),          &
     &                     Work(itaup),Work(iwork),Lwork-iwork+1,ierr)
               CALL DLACPY('L',N,N,Vt,Ldvt,Work(ir),ldwrkr)
!
!                 Generate left vectors bidiagonalizing R in WORK(IR)
!                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
!
               CALL DORGBR('Q',N,N,N,Work(ir),ldwrkr,Work(itauq),       &
     &                     Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate right vectors bidiagonalizing R in VT
!                 (Workspace: need N*N + 4*N-1, prefer N*N + 3*N + (N-1)*NB)
!
               CALL DORGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),   &
     &                     Lwork-iwork+1,ierr)
               iwork = ie + N
!
!                 Perform bidiagonal QR iteration, computing left
!                 singular vectors of R in WORK(IR) and computing right
!                 singular vectors of R in VT
!                 (Workspace: need N*N + BDSPAC)
!
               CALL DBDSQR('U',N,N,N,0,S,Work(ie),Vt,Ldvt,Work(ir),     &
     &                     ldwrkr,dum,1,Work(iwork),Info)
               iu = ie + N
!
!                 Multiply Q in A by left singular vectors of R in
!                 WORK(IR), storing result in WORK(IU) and copying to A
!                 (Workspace: need N*N + 2*N, prefer N*N + M*N + N)
!
               DO i = 1 , M , ldwrku
                  chunk = MIN(M-i+1,ldwrku)
                  CALL DGEMM('N','N',chunk,N,N,ONE,A(i,1),Lda,Work(ir), &
     &                       ldwrkr,ZERO,Work(iu),ldwrku)
                  CALL DLACPY('F',chunk,N,Work(iu),ldwrku,A(i,1),Lda)
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
!                 (Workspace: need 2*N, prefer N + N*NB)
!
               CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),            &
     &                     Lwork-iwork+1,ierr)
!
!                 Copy R to VT, zeroing out below it
!
               CALL DLACPY('U',N,N,A,Lda,Vt,Ldvt)
               IF ( N>1 ) CALL DLASET('L',N-1,N-1,ZERO,ZERO,Vt(2,1),    &
     &                                Ldvt)
!
!                 Generate Q in A
!                 (Workspace: need 2*N, prefer N + N*NB)
!
               CALL DORGQR(M,N,N,A,Lda,Work(itau),Work(iwork),          &
     &                     Lwork-iwork+1,ierr)
               ie = itau
               itauq = ie + N
               itaup = itauq + N
               iwork = itaup + N
!
!                 Bidiagonalize R in VT
!                 (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!
               CALL DGEBRD(N,N,Vt,Ldvt,S,Work(ie),Work(itauq),          &
     &                     Work(itaup),Work(iwork),Lwork-iwork+1,ierr)
!
!                 Multiply Q in A by left vectors bidiagonalizing R
!                 (Workspace: need 3*N + M, prefer 3*N + M*NB)
!
               CALL DORMBR('Q','R','N',M,N,N,Vt,Ldvt,Work(itauq),A,Lda, &
     &                     Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate right vectors bidiagonalizing R in VT
!                 (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
!
               CALL DORGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),   &
     &                     Lwork-iwork+1,ierr)
               iwork = ie + N
!
!                 Perform bidiagonal QR iteration, computing left
!                 singular vectors of A in A and computing right
!                 singular vectors of A in VT
!                 (Workspace: need BDSPAC)
!
               CALL DBDSQR('U',N,N,M,0,S,Work(ie),Vt,Ldvt,A,Lda,dum,1,  &
     &                     Work(iwork),Info)
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
               IF ( Lwork>=N*N+MAX(4*N,bdspac) ) THEN
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
!                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!
                  CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R to WORK(IR), zeroing out below it
!
                  CALL DLACPY('U',N,N,A,Lda,Work(ir),ldwrkr)
                  CALL DLASET('L',N-1,N-1,ZERO,ZERO,Work(ir+1),ldwrkr)
!
!                    Generate Q in A
!                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!
                  CALL DORGQR(M,N,N,A,Lda,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + N
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in WORK(IR)
!                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
!
                  CALL DGEBRD(N,N,Work(ir),ldwrkr,S,Work(ie),Work(itauq)&
     &                        ,Work(itaup),Work(iwork),Lwork-iwork+1,   &
     &                        ierr)
!
!                    Generate left vectors bidiagonalizing R in WORK(IR)
!                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
!
                  CALL DORGBR('Q',N,N,N,Work(ir),ldwrkr,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  iwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of R in WORK(IR)
!                    (Workspace: need N*N + BDSPAC)
!
                  CALL DBDSQR('U',N,0,N,0,S,Work(ie),dum,1,Work(ir),    &
     &                        ldwrkr,dum,1,Work(iwork),Info)
!
!                    Multiply Q in A by left singular vectors of R in
!                    WORK(IR), storing result in U
!                    (Workspace: need N*N)
!
                  CALL DGEMM('N','N',M,N,N,ONE,A,Lda,Work(ir),ldwrkr,   &
     &                       ZERO,U,Ldu)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (Workspace: need 2*N, prefer N + N*NB)
!
                  CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (Workspace: need 2*N, prefer N + N*NB)
!
                  CALL DORGQR(M,N,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + N
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Zero out below R in A
!
                  IF ( N>1 ) CALL DLASET('L',N-1,N-1,ZERO,ZERO,A(2,1),  &
     &                 Lda)
!
!                    Bidiagonalize R in A
!                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!
                  CALL DGEBRD(N,N,A,Lda,S,Work(ie),Work(itauq),         &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply Q in U by left vectors bidiagonalizing R
!                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
!
                  CALL DORMBR('Q','R','N',M,N,N,A,Lda,Work(itauq),U,Ldu,&
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  iwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U
!                    (Workspace: need BDSPAC)
!
                  CALL DBDSQR('U',N,0,M,0,S,Work(ie),dum,1,U,Ldu,dum,1, &
     &                        Work(iwork),Info)
!
               ENDIF
!
            ELSEIF ( wntvo ) THEN
!
!                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
!                 N left singular vectors to be computed in U and
!                 N right singular vectors to be overwritten on A
!
               IF ( Lwork>=2*N*N+MAX(4*N,bdspac) ) THEN
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
!                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB)
!
                  CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R to WORK(IU), zeroing out below it
!
                  CALL DLACPY('U',N,N,A,Lda,Work(iu),ldwrku)
                  CALL DLASET('L',N-1,N-1,ZERO,ZERO,Work(iu+1),ldwrku)
!
!                    Generate Q in A
!                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB)
!
                  CALL DORGQR(M,N,N,A,Lda,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + N
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in WORK(IU), copying result to
!                    WORK(IR)
!                    (Workspace: need 2*N*N + 4*N,
!                                prefer 2*N*N+3*N+2*N*NB)
!
                  CALL DGEBRD(N,N,Work(iu),ldwrku,S,Work(ie),Work(itauq)&
     &                        ,Work(itaup),Work(iwork),Lwork-iwork+1,   &
     &                        ierr)
                  CALL DLACPY('U',N,N,Work(iu),ldwrku,Work(ir),ldwrkr)
!
!                    Generate left bidiagonalizing vectors in WORK(IU)
!                    (Workspace: need 2*N*N + 4*N, prefer 2*N*N + 3*N + N*NB)
!
                  CALL DORGBR('Q',N,N,N,Work(iu),ldwrku,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in WORK(IR)
!                    (Workspace: need 2*N*N + 4*N-1,
!                                prefer 2*N*N+3*N+(N-1)*NB)
!
                  CALL DORGBR('P',N,N,N,Work(ir),ldwrkr,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  iwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of R in WORK(IU) and computing
!                    right singular vectors of R in WORK(IR)
!                    (Workspace: need 2*N*N + BDSPAC)
!
                  CALL DBDSQR('U',N,N,N,0,S,Work(ie),Work(ir),ldwrkr,   &
     &                        Work(iu),ldwrku,dum,1,Work(iwork),Info)
!
!                    Multiply Q in A by left singular vectors of R in
!                    WORK(IU), storing result in U
!                    (Workspace: need N*N)
!
                  CALL DGEMM('N','N',M,N,N,ONE,A,Lda,Work(iu),ldwrku,   &
     &                       ZERO,U,Ldu)
!
!                    Copy right singular vectors of R to A
!                    (Workspace: need N*N)
!
                  CALL DLACPY('F',N,N,Work(ir),ldwrkr,A,Lda)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (Workspace: need 2*N, prefer N + N*NB)
!
                  CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (Workspace: need 2*N, prefer N + N*NB)
!
                  CALL DORGQR(M,N,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + N
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Zero out below R in A
!
                  IF ( N>1 ) CALL DLASET('L',N-1,N-1,ZERO,ZERO,A(2,1),  &
     &                 Lda)
!
!                    Bidiagonalize R in A
!                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!
                  CALL DGEBRD(N,N,A,Lda,S,Work(ie),Work(itauq),         &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply Q in U by left vectors bidiagonalizing R
!                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
!
                  CALL DORMBR('Q','R','N',M,N,N,A,Lda,Work(itauq),U,Ldu,&
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right vectors bidiagonalizing R in A
!                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
!
                  CALL DORGBR('P',N,N,N,A,Lda,Work(itaup),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  iwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U and computing right
!                    singular vectors of A in A
!                    (Workspace: need BDSPAC)
!
                  CALL DBDSQR('U',N,N,M,0,S,Work(ie),A,Lda,U,Ldu,dum,1, &
     &                        Work(iwork),Info)
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
               IF ( Lwork>=N*N+MAX(4*N,bdspac) ) THEN
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
!                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!
                  CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R to WORK(IU), zeroing out below it
!
                  CALL DLACPY('U',N,N,A,Lda,Work(iu),ldwrku)
                  CALL DLASET('L',N-1,N-1,ZERO,ZERO,Work(iu+1),ldwrku)
!
!                    Generate Q in A
!                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!
                  CALL DORGQR(M,N,N,A,Lda,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + N
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in WORK(IU), copying result to VT
!                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
!
                  CALL DGEBRD(N,N,Work(iu),ldwrku,S,Work(ie),Work(itauq)&
     &                        ,Work(itaup),Work(iwork),Lwork-iwork+1,   &
     &                        ierr)
                  CALL DLACPY('U',N,N,Work(iu),ldwrku,Vt,Ldvt)
!
!                    Generate left bidiagonalizing vectors in WORK(IU)
!                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
!
                  CALL DORGBR('Q',N,N,N,Work(iu),ldwrku,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in VT
!                    (Workspace: need N*N + 4*N-1,
!                                prefer N*N+3*N+(N-1)*NB)
!
                  CALL DORGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),&
     &                        Lwork-iwork+1,ierr)
                  iwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of R in WORK(IU) and computing
!                    right singular vectors of R in VT
!                    (Workspace: need N*N + BDSPAC)
!
                  CALL DBDSQR('U',N,N,N,0,S,Work(ie),Vt,Ldvt,Work(iu),  &
     &                        ldwrku,dum,1,Work(iwork),Info)
!
!                    Multiply Q in A by left singular vectors of R in
!                    WORK(IU), storing result in U
!                    (Workspace: need N*N)
!
                  CALL DGEMM('N','N',M,N,N,ONE,A,Lda,Work(iu),ldwrku,   &
     &                       ZERO,U,Ldu)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (Workspace: need 2*N, prefer N + N*NB)
!
                  CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (Workspace: need 2*N, prefer N + N*NB)
!
                  CALL DORGQR(M,N,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R to VT, zeroing out below it
!
                  CALL DLACPY('U',N,N,A,Lda,Vt,Ldvt)
                  IF ( N>1 ) CALL DLASET('L',N-1,N-1,ZERO,ZERO,Vt(2,1), &
     &                 Ldvt)
                  ie = itau
                  itauq = ie + N
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in VT
!                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!
                  CALL DGEBRD(N,N,Vt,Ldvt,S,Work(ie),Work(itauq),       &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply Q in U by left bidiagonalizing vectors
!                    in VT
!                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
!
                  CALL DORMBR('Q','R','N',M,N,N,Vt,Ldvt,Work(itauq),U,  &
     &                        Ldu,Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in VT
!                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
!
                  CALL DORGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),&
     &                        Lwork-iwork+1,ierr)
                  iwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U and computing right
!                    singular vectors of A in VT
!                    (Workspace: need BDSPAC)
!
                  CALL DBDSQR('U',N,N,M,0,S,Work(ie),Vt,Ldvt,U,Ldu,dum, &
     &                        1,Work(iwork),Info)
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
               IF ( Lwork>=N*N+MAX(N+M,4*N,bdspac) ) THEN
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
!                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!
                  CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Copy R to WORK(IR), zeroing out below it
!
                  CALL DLACPY('U',N,N,A,Lda,Work(ir),ldwrkr)
                  CALL DLASET('L',N-1,N-1,ZERO,ZERO,Work(ir+1),ldwrkr)
!
!                    Generate Q in U
!                    (Workspace: need N*N + N + M, prefer N*N + N + M*NB)
!
                  CALL DORGQR(M,M,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + N
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in WORK(IR)
!                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
!
                  CALL DGEBRD(N,N,Work(ir),ldwrkr,S,Work(ie),Work(itauq)&
     &                        ,Work(itaup),Work(iwork),Lwork-iwork+1,   &
     &                        ierr)
!
!                    Generate left bidiagonalizing vectors in WORK(IR)
!                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
!
                  CALL DORGBR('Q',N,N,N,Work(ir),ldwrkr,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  iwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of R in WORK(IR)
!                    (Workspace: need N*N + BDSPAC)
!
                  CALL DBDSQR('U',N,0,N,0,S,Work(ie),dum,1,Work(ir),    &
     &                        ldwrkr,dum,1,Work(iwork),Info)
!
!                    Multiply Q in U by left singular vectors of R in
!                    WORK(IR), storing result in A
!                    (Workspace: need N*N)
!
                  CALL DGEMM('N','N',M,N,N,ONE,U,Ldu,Work(ir),ldwrkr,   &
     &                       ZERO,A,Lda)
!
!                    Copy left singular vectors of A from A to U
!
                  CALL DLACPY('F',M,N,A,Lda,U,Ldu)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (Workspace: need 2*N, prefer N + N*NB)
!
                  CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (Workspace: need N + M, prefer N + M*NB)
!
                  CALL DORGQR(M,M,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + N
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Zero out below R in A
!
                  IF ( N>1 ) CALL DLASET('L',N-1,N-1,ZERO,ZERO,A(2,1),  &
     &                 Lda)
!
!                    Bidiagonalize R in A
!                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!
                  CALL DGEBRD(N,N,A,Lda,S,Work(ie),Work(itauq),         &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply Q in U by left bidiagonalizing vectors
!                    in A
!                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
!
                  CALL DORMBR('Q','R','N',M,N,N,A,Lda,Work(itauq),U,Ldu,&
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  iwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U
!                    (Workspace: need BDSPAC)
!
                  CALL DBDSQR('U',N,0,M,0,S,Work(ie),dum,1,U,Ldu,dum,1, &
     &                        Work(iwork),Info)
!
               ENDIF
!
            ELSEIF ( wntvo ) THEN
!
!                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
!                 M left singular vectors to be computed in U and
!                 N right singular vectors to be overwritten on A
!
               IF ( Lwork>=2*N*N+MAX(N+M,4*N,bdspac) ) THEN
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
!                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB)
!
                  CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (Workspace: need 2*N*N + N + M, prefer 2*N*N + N + M*NB)
!
                  CALL DORGQR(M,M,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R to WORK(IU), zeroing out below it
!
                  CALL DLACPY('U',N,N,A,Lda,Work(iu),ldwrku)
                  CALL DLASET('L',N-1,N-1,ZERO,ZERO,Work(iu+1),ldwrku)
                  ie = itau
                  itauq = ie + N
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in WORK(IU), copying result to
!                    WORK(IR)
!                    (Workspace: need 2*N*N + 4*N,
!                                prefer 2*N*N+3*N+2*N*NB)
!
                  CALL DGEBRD(N,N,Work(iu),ldwrku,S,Work(ie),Work(itauq)&
     &                        ,Work(itaup),Work(iwork),Lwork-iwork+1,   &
     &                        ierr)
                  CALL DLACPY('U',N,N,Work(iu),ldwrku,Work(ir),ldwrkr)
!
!                    Generate left bidiagonalizing vectors in WORK(IU)
!                    (Workspace: need 2*N*N + 4*N, prefer 2*N*N + 3*N + N*NB)
!
                  CALL DORGBR('Q',N,N,N,Work(iu),ldwrku,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in WORK(IR)
!                    (Workspace: need 2*N*N + 4*N-1,
!                                prefer 2*N*N+3*N+(N-1)*NB)
!
                  CALL DORGBR('P',N,N,N,Work(ir),ldwrkr,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  iwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of R in WORK(IU) and computing
!                    right singular vectors of R in WORK(IR)
!                    (Workspace: need 2*N*N + BDSPAC)
!
                  CALL DBDSQR('U',N,N,N,0,S,Work(ie),Work(ir),ldwrkr,   &
     &                        Work(iu),ldwrku,dum,1,Work(iwork),Info)
!
!                    Multiply Q in U by left singular vectors of R in
!                    WORK(IU), storing result in A
!                    (Workspace: need N*N)
!
                  CALL DGEMM('N','N',M,N,N,ONE,U,Ldu,Work(iu),ldwrku,   &
     &                       ZERO,A,Lda)
!
!                    Copy left singular vectors of A from A to U
!
                  CALL DLACPY('F',M,N,A,Lda,U,Ldu)
!
!                    Copy right singular vectors of R from WORK(IR) to A
!
                  CALL DLACPY('F',N,N,Work(ir),ldwrkr,A,Lda)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (Workspace: need 2*N, prefer N + N*NB)
!
                  CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (Workspace: need N + M, prefer N + M*NB)
!
                  CALL DORGQR(M,M,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + N
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Zero out below R in A
!
                  IF ( N>1 ) CALL DLASET('L',N-1,N-1,ZERO,ZERO,A(2,1),  &
     &                 Lda)
!
!                    Bidiagonalize R in A
!                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!
                  CALL DGEBRD(N,N,A,Lda,S,Work(ie),Work(itauq),         &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply Q in U by left bidiagonalizing vectors
!                    in A
!                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
!
                  CALL DORMBR('Q','R','N',M,N,N,A,Lda,Work(itauq),U,Ldu,&
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in A
!                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
!
                  CALL DORGBR('P',N,N,N,A,Lda,Work(itaup),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  iwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U and computing right
!                    singular vectors of A in A
!                    (Workspace: need BDSPAC)
!
                  CALL DBDSQR('U',N,N,M,0,S,Work(ie),A,Lda,U,Ldu,dum,1, &
     &                        Work(iwork),Info)
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
               IF ( Lwork>=N*N+MAX(N+M,4*N,bdspac) ) THEN
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
!                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!
                  CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (Workspace: need N*N + N + M, prefer N*N + N + M*NB)
!
                  CALL DORGQR(M,M,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R to WORK(IU), zeroing out below it
!
                  CALL DLACPY('U',N,N,A,Lda,Work(iu),ldwrku)
                  CALL DLASET('L',N-1,N-1,ZERO,ZERO,Work(iu+1),ldwrku)
                  ie = itau
                  itauq = ie + N
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in WORK(IU), copying result to VT
!                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
!
                  CALL DGEBRD(N,N,Work(iu),ldwrku,S,Work(ie),Work(itauq)&
     &                        ,Work(itaup),Work(iwork),Lwork-iwork+1,   &
     &                        ierr)
                  CALL DLACPY('U',N,N,Work(iu),ldwrku,Vt,Ldvt)
!
!                    Generate left bidiagonalizing vectors in WORK(IU)
!                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
!
                  CALL DORGBR('Q',N,N,N,Work(iu),ldwrku,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in VT
!                    (Workspace: need N*N + 4*N-1,
!                                prefer N*N+3*N+(N-1)*NB)
!
                  CALL DORGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),&
     &                        Lwork-iwork+1,ierr)
                  iwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of R in WORK(IU) and computing
!                    right singular vectors of R in VT
!                    (Workspace: need N*N + BDSPAC)
!
                  CALL DBDSQR('U',N,N,N,0,S,Work(ie),Vt,Ldvt,Work(iu),  &
     &                        ldwrku,dum,1,Work(iwork),Info)
!
!                    Multiply Q in U by left singular vectors of R in
!                    WORK(IU), storing result in A
!                    (Workspace: need N*N)
!
                  CALL DGEMM('N','N',M,N,N,ONE,U,Ldu,Work(iu),ldwrku,   &
     &                       ZERO,A,Lda)
!
!                    Copy left singular vectors of A from A to U
!
                  CALL DLACPY('F',M,N,A,Lda,U,Ldu)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + N
!
!                    Compute A=Q*R, copying result to U
!                    (Workspace: need 2*N, prefer N + N*NB)
!
                  CALL DGEQRF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('L',M,N,A,Lda,U,Ldu)
!
!                    Generate Q in U
!                    (Workspace: need N + M, prefer N + M*NB)
!
                  CALL DORGQR(M,M,N,U,Ldu,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy R from A to VT, zeroing out below it
!
                  CALL DLACPY('U',N,N,A,Lda,Vt,Ldvt)
                  IF ( N>1 ) CALL DLASET('L',N-1,N-1,ZERO,ZERO,Vt(2,1), &
     &                 Ldvt)
                  ie = itau
                  itauq = ie + N
                  itaup = itauq + N
                  iwork = itaup + N
!
!                    Bidiagonalize R in VT
!                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!
                  CALL DGEBRD(N,N,Vt,Ldvt,S,Work(ie),Work(itauq),       &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply Q in U by left bidiagonalizing vectors
!                    in VT
!                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
!
                  CALL DORMBR('Q','R','N',M,N,N,Vt,Ldvt,Work(itauq),U,  &
     &                        Ldu,Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate right bidiagonalizing vectors in VT
!                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
!
                  CALL DORGBR('P',N,N,N,Vt,Ldvt,Work(itaup),Work(iwork),&
     &                        Lwork-iwork+1,ierr)
                  iwork = ie + N
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U and computing right
!                    singular vectors of A in VT
!                    (Workspace: need BDSPAC)
!
                  CALL DBDSQR('U',N,N,M,0,S,Work(ie),Vt,Ldvt,U,Ldu,dum, &
     &                        1,Work(iwork),Info)
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
!              (Workspace: need 2*M, prefer M + M*NB)
!
            CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),Lwork-iwork+1, &
     &                  ierr)
!
!              Zero out above L
!
            CALL DLASET('U',M-1,M-1,ZERO,ZERO,A(1,2),Lda)
            ie = 1
            itauq = ie + M
            itaup = itauq + M
            iwork = itaup + M
!
!              Bidiagonalize L in A
!              (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!
            CALL DGEBRD(M,M,A,Lda,S,Work(ie),Work(itauq),Work(itaup),   &
     &                  Work(iwork),Lwork-iwork+1,ierr)
!
!                 If left singular vectors desired, generate Q
!                 (Workspace: need 4*M, prefer 3*M + M*NB)
!
            IF ( wntuo .OR. wntuas )                                    &
     &           CALL DORGBR('Q',M,M,M,A,Lda,Work(itauq),Work(iwork),   &
     &           Lwork-iwork+1,ierr)
            iwork = ie + M
            nru = 0
            IF ( wntuo .OR. wntuas ) nru = M
!
!              Perform bidiagonal QR iteration, computing left singular
!              vectors of A in A if desired
!              (Workspace: need BDSPAC)
!
            CALL DBDSQR('U',M,0,nru,0,S,Work(ie),dum,1,A,Lda,dum,1,     &
     &                  Work(iwork),Info)
!
!              If left singular vectors desired in U, copy them there
!
            IF ( wntuas ) CALL DLACPY('F',M,M,A,Lda,U,Ldu)
!
         ELSEIF ( wntvo .AND. wntun ) THEN
!
!              Path 2t(N much larger than M, JOBU='N', JOBVT='O')
!              M right singular vectors to be overwritten on A and
!              no left singular vectors to be computed
!
            IF ( Lwork>=M*M+MAX(4*M,bdspac) ) THEN
!
!                 Sufficient workspace for a fast algorithm
!
               ir = 1
               IF ( Lwork>=MAX(wrkbl,Lda*N+M)+Lda*M ) THEN
!
!                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
!
                  ldwrku = Lda
                  chunk = N
                  ldwrkr = Lda
               ELSEIF ( Lwork>=MAX(wrkbl,Lda*N+M)+M*M ) THEN
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
                  chunk = (Lwork-M*M-M)/M
                  ldwrkr = M
               ENDIF
               itau = ir + ldwrkr*M
               iwork = itau + M
!
!                 Compute A=L*Q
!                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!
               CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),            &
     &                     Lwork-iwork+1,ierr)
!
!                 Copy L to WORK(IR) and zero out above it
!
               CALL DLACPY('L',M,M,A,Lda,Work(ir),ldwrkr)
               CALL DLASET('U',M-1,M-1,ZERO,ZERO,Work(ir+ldwrkr),ldwrkr)
!
!                 Generate Q in A
!                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!
               CALL DORGLQ(M,N,M,A,Lda,Work(itau),Work(iwork),          &
     &                     Lwork-iwork+1,ierr)
               ie = itau
               itauq = ie + M
               itaup = itauq + M
               iwork = itaup + M
!
!                 Bidiagonalize L in WORK(IR)
!                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
!
               CALL DGEBRD(M,M,Work(ir),ldwrkr,S,Work(ie),Work(itauq),  &
     &                     Work(itaup),Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate right vectors bidiagonalizing L
!                 (Workspace: need M*M + 4*M-1, prefer M*M + 3*M + (M-1)*NB)
!
               CALL DORGBR('P',M,M,M,Work(ir),ldwrkr,Work(itaup),       &
     &                     Work(iwork),Lwork-iwork+1,ierr)
               iwork = ie + M
!
!                 Perform bidiagonal QR iteration, computing right
!                 singular vectors of L in WORK(IR)
!                 (Workspace: need M*M + BDSPAC)
!
               CALL DBDSQR('U',M,M,0,0,S,Work(ie),Work(ir),ldwrkr,dum,1,&
     &                     dum,1,Work(iwork),Info)
               iu = ie + M
!
!                 Multiply right singular vectors of L in WORK(IR) by Q
!                 in A, storing result in WORK(IU) and copying to A
!                 (Workspace: need M*M + 2*M, prefer M*M + M*N + M)
!
               DO i = 1 , N , chunk
                  blk = MIN(N-i+1,chunk)
                  CALL DGEMM('N','N',M,blk,M,ONE,Work(ir),ldwrkr,A(1,i),&
     &                       Lda,ZERO,Work(iu),ldwrku)
                  CALL DLACPY('F',M,blk,Work(iu),ldwrku,A(1,i),Lda)
               ENDDO
!
            ELSE
!
!                 Insufficient workspace for a fast algorithm
!
               ie = 1
               itauq = ie + M
               itaup = itauq + M
               iwork = itaup + M
!
!                 Bidiagonalize A
!                 (Workspace: need 3*M + N, prefer 3*M + (M + N)*NB)
!
               CALL DGEBRD(M,N,A,Lda,S,Work(ie),Work(itauq),Work(itaup),&
     &                     Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate right vectors bidiagonalizing A
!                 (Workspace: need 4*M, prefer 3*M + M*NB)
!
               CALL DORGBR('P',M,N,M,A,Lda,Work(itaup),Work(iwork),     &
     &                     Lwork-iwork+1,ierr)
               iwork = ie + M
!
!                 Perform bidiagonal QR iteration, computing right
!                 singular vectors of A in A
!                 (Workspace: need BDSPAC)
!
               CALL DBDSQR('L',M,N,0,0,S,Work(ie),A,Lda,dum,1,dum,1,    &
     &                     Work(iwork),Info)
!
            ENDIF
!
         ELSEIF ( wntvo .AND. wntuas ) THEN
!
!              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
!              M right singular vectors to be overwritten on A and
!              M left singular vectors to be computed in U
!
            IF ( Lwork>=M*M+MAX(4*M,bdspac) ) THEN
!
!                 Sufficient workspace for a fast algorithm
!
               ir = 1
               IF ( Lwork>=MAX(wrkbl,Lda*N+M)+Lda*M ) THEN
!
!                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
!
                  ldwrku = Lda
                  chunk = N
                  ldwrkr = Lda
               ELSEIF ( Lwork>=MAX(wrkbl,Lda*N+M)+M*M ) THEN
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
                  chunk = (Lwork-M*M-M)/M
                  ldwrkr = M
               ENDIF
               itau = ir + ldwrkr*M
               iwork = itau + M
!
!                 Compute A=L*Q
!                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!
               CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),            &
     &                     Lwork-iwork+1,ierr)
!
!                 Copy L to U, zeroing about above it
!
               CALL DLACPY('L',M,M,A,Lda,U,Ldu)
               CALL DLASET('U',M-1,M-1,ZERO,ZERO,U(1,2),Ldu)
!
!                 Generate Q in A
!                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!
               CALL DORGLQ(M,N,M,A,Lda,Work(itau),Work(iwork),          &
     &                     Lwork-iwork+1,ierr)
               ie = itau
               itauq = ie + M
               itaup = itauq + M
               iwork = itaup + M
!
!                 Bidiagonalize L in U, copying result to WORK(IR)
!                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
!
               CALL DGEBRD(M,M,U,Ldu,S,Work(ie),Work(itauq),Work(itaup),&
     &                     Work(iwork),Lwork-iwork+1,ierr)
               CALL DLACPY('U',M,M,U,Ldu,Work(ir),ldwrkr)
!
!                 Generate right vectors bidiagonalizing L in WORK(IR)
!                 (Workspace: need M*M + 4*M-1, prefer M*M + 3*M + (M-1)*NB)
!
               CALL DORGBR('P',M,M,M,Work(ir),ldwrkr,Work(itaup),       &
     &                     Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate left vectors bidiagonalizing L in U
!                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB)
!
               CALL DORGBR('Q',M,M,M,U,Ldu,Work(itauq),Work(iwork),     &
     &                     Lwork-iwork+1,ierr)
               iwork = ie + M
!
!                 Perform bidiagonal QR iteration, computing left
!                 singular vectors of L in U, and computing right
!                 singular vectors of L in WORK(IR)
!                 (Workspace: need M*M + BDSPAC)
!
               CALL DBDSQR('U',M,M,M,0,S,Work(ie),Work(ir),ldwrkr,U,Ldu,&
     &                     dum,1,Work(iwork),Info)
               iu = ie + M
!
!                 Multiply right singular vectors of L in WORK(IR) by Q
!                 in A, storing result in WORK(IU) and copying to A
!                 (Workspace: need M*M + 2*M, prefer M*M + M*N + M))
!
               DO i = 1 , N , chunk
                  blk = MIN(N-i+1,chunk)
                  CALL DGEMM('N','N',M,blk,M,ONE,Work(ir),ldwrkr,A(1,i),&
     &                       Lda,ZERO,Work(iu),ldwrku)
                  CALL DLACPY('F',M,blk,Work(iu),ldwrku,A(1,i),Lda)
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
!                 (Workspace: need 2*M, prefer M + M*NB)
!
               CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),            &
     &                     Lwork-iwork+1,ierr)
!
!                 Copy L to U, zeroing out above it
!
               CALL DLACPY('L',M,M,A,Lda,U,Ldu)
               CALL DLASET('U',M-1,M-1,ZERO,ZERO,U(1,2),Ldu)
!
!                 Generate Q in A
!                 (Workspace: need 2*M, prefer M + M*NB)
!
               CALL DORGLQ(M,N,M,A,Lda,Work(itau),Work(iwork),          &
     &                     Lwork-iwork+1,ierr)
               ie = itau
               itauq = ie + M
               itaup = itauq + M
               iwork = itaup + M
!
!                 Bidiagonalize L in U
!                 (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!
               CALL DGEBRD(M,M,U,Ldu,S,Work(ie),Work(itauq),Work(itaup),&
     &                     Work(iwork),Lwork-iwork+1,ierr)
!
!                 Multiply right vectors bidiagonalizing L by Q in A
!                 (Workspace: need 3*M + N, prefer 3*M + N*NB)
!
               CALL DORMBR('P','L','T',M,N,M,U,Ldu,Work(itaup),A,Lda,   &
     &                     Work(iwork),Lwork-iwork+1,ierr)
!
!                 Generate left vectors bidiagonalizing L in U
!                 (Workspace: need 4*M, prefer 3*M + M*NB)
!
               CALL DORGBR('Q',M,M,M,U,Ldu,Work(itauq),Work(iwork),     &
     &                     Lwork-iwork+1,ierr)
               iwork = ie + M
!
!                 Perform bidiagonal QR iteration, computing left
!                 singular vectors of A in U and computing right
!                 singular vectors of A in A
!                 (Workspace: need BDSPAC)
!
               CALL DBDSQR('U',M,N,M,0,S,Work(ie),A,Lda,U,Ldu,dum,1,    &
     &                     Work(iwork),Info)
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
               IF ( Lwork>=M*M+MAX(4*M,bdspac) ) THEN
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
!                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!
                  CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to WORK(IR), zeroing out above it
!
                  CALL DLACPY('L',M,M,A,Lda,Work(ir),ldwrkr)
                  CALL DLASET('U',M-1,M-1,ZERO,ZERO,Work(ir+ldwrkr),    &
     &                        ldwrkr)
!
!                    Generate Q in A
!                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!
                  CALL DORGLQ(M,N,M,A,Lda,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + M
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in WORK(IR)
!                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
!
                  CALL DGEBRD(M,M,Work(ir),ldwrkr,S,Work(ie),Work(itauq)&
     &                        ,Work(itaup),Work(iwork),Lwork-iwork+1,   &
     &                        ierr)
!
!                    Generate right vectors bidiagonalizing L in
!                    WORK(IR)
!                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + (M-1)*NB)
!
                  CALL DORGBR('P',M,M,M,Work(ir),ldwrkr,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  iwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing right
!                    singular vectors of L in WORK(IR)
!                    (Workspace: need M*M + BDSPAC)
!
                  CALL DBDSQR('U',M,M,0,0,S,Work(ie),Work(ir),ldwrkr,   &
     &                        dum,1,dum,1,Work(iwork),Info)
!
!                    Multiply right singular vectors of L in WORK(IR) by
!                    Q in A, storing result in VT
!                    (Workspace: need M*M)
!
                  CALL DGEMM('N','N',M,N,M,ONE,Work(ir),ldwrkr,A,Lda,   &
     &                       ZERO,Vt,Ldvt)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + M
!
!                    Compute A=L*Q
!                    (Workspace: need 2*M, prefer M + M*NB)
!
                  CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy result to VT
!
                  CALL DLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (Workspace: need 2*M, prefer M + M*NB)
!
                  CALL DORGLQ(M,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + M
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Zero out above L in A
!
                  CALL DLASET('U',M-1,M-1,ZERO,ZERO,A(1,2),Lda)
!
!                    Bidiagonalize L in A
!                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!
                  CALL DGEBRD(M,M,A,Lda,S,Work(ie),Work(itauq),         &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply right vectors bidiagonalizing L by Q in VT
!                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
!
                  CALL DORMBR('P','L','T',M,N,M,A,Lda,Work(itaup),Vt,   &
     &                        Ldvt,Work(iwork),Lwork-iwork+1,ierr)
                  iwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing right
!                    singular vectors of A in VT
!                    (Workspace: need BDSPAC)
!
                  CALL DBDSQR('U',M,N,0,0,S,Work(ie),Vt,Ldvt,dum,1,dum, &
     &                        1,Work(iwork),Info)
!
               ENDIF
!
            ELSEIF ( wntuo ) THEN
!
!                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
!                 M right singular vectors to be computed in VT and
!                 M left singular vectors to be overwritten on A
!
               IF ( Lwork>=2*M*M+MAX(4*M,bdspac) ) THEN
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
!                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB)
!
                  CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to WORK(IU), zeroing out below it
!
                  CALL DLACPY('L',M,M,A,Lda,Work(iu),ldwrku)
                  CALL DLASET('U',M-1,M-1,ZERO,ZERO,Work(iu+ldwrku),    &
     &                        ldwrku)
!
!                    Generate Q in A
!                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB)
!
                  CALL DORGLQ(M,N,M,A,Lda,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + M
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in WORK(IU), copying result to
!                    WORK(IR)
!                    (Workspace: need 2*M*M + 4*M,
!                                prefer 2*M*M+3*M+2*M*NB)
!
                  CALL DGEBRD(M,M,Work(iu),ldwrku,S,Work(ie),Work(itauq)&
     &                        ,Work(itaup),Work(iwork),Lwork-iwork+1,   &
     &                        ierr)
                  CALL DLACPY('L',M,M,Work(iu),ldwrku,Work(ir),ldwrkr)
!
!                    Generate right bidiagonalizing vectors in WORK(IU)
!                    (Workspace: need 2*M*M + 4*M-1,
!                                prefer 2*M*M+3*M+(M-1)*NB)
!
                  CALL DORGBR('P',M,M,M,Work(iu),ldwrku,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in WORK(IR)
!                    (Workspace: need 2*M*M + 4*M, prefer 2*M*M + 3*M + M*NB)
!
                  CALL DORGBR('Q',M,M,M,Work(ir),ldwrkr,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  iwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of L in WORK(IR) and computing
!                    right singular vectors of L in WORK(IU)
!                    (Workspace: need 2*M*M + BDSPAC)
!
                  CALL DBDSQR('U',M,M,M,0,S,Work(ie),Work(iu),ldwrku,   &
     &                        Work(ir),ldwrkr,dum,1,Work(iwork),Info)
!
!                    Multiply right singular vectors of L in WORK(IU) by
!                    Q in A, storing result in VT
!                    (Workspace: need M*M)
!
                  CALL DGEMM('N','N',M,N,M,ONE,Work(iu),ldwrku,A,Lda,   &
     &                       ZERO,Vt,Ldvt)
!
!                    Copy left singular vectors of L to A
!                    (Workspace: need M*M)
!
                  CALL DLACPY('F',M,M,Work(ir),ldwrkr,A,Lda)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (Workspace: need 2*M, prefer M + M*NB)
!
                  CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (Workspace: need 2*M, prefer M + M*NB)
!
                  CALL DORGLQ(M,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + M
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Zero out above L in A
!
                  CALL DLASET('U',M-1,M-1,ZERO,ZERO,A(1,2),Lda)
!
!                    Bidiagonalize L in A
!                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!
                  CALL DGEBRD(M,M,A,Lda,S,Work(ie),Work(itauq),         &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply right vectors bidiagonalizing L by Q in VT
!                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
!
                  CALL DORMBR('P','L','T',M,N,M,A,Lda,Work(itaup),Vt,   &
     &                        Ldvt,Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors of L in A
!                    (Workspace: need 4*M, prefer 3*M + M*NB)
!
                  CALL DORGBR('Q',M,M,M,A,Lda,Work(itauq),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  iwork = ie + M
!
!                    Perform bidiagonal QR iteration, compute left
!                    singular vectors of A in A and compute right
!                    singular vectors of A in VT
!                    (Workspace: need BDSPAC)
!
                  CALL DBDSQR('U',M,N,M,0,S,Work(ie),Vt,Ldvt,A,Lda,dum, &
     &                        1,Work(iwork),Info)
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
               IF ( Lwork>=M*M+MAX(4*M,bdspac) ) THEN
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
!                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!
                  CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to WORK(IU), zeroing out above it
!
                  CALL DLACPY('L',M,M,A,Lda,Work(iu),ldwrku)
                  CALL DLASET('U',M-1,M-1,ZERO,ZERO,Work(iu+ldwrku),    &
     &                        ldwrku)
!
!                    Generate Q in A
!                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!
                  CALL DORGLQ(M,N,M,A,Lda,Work(itau),Work(iwork),       &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + M
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in WORK(IU), copying result to U
!                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
!
                  CALL DGEBRD(M,M,Work(iu),ldwrku,S,Work(ie),Work(itauq)&
     &                        ,Work(itaup),Work(iwork),Lwork-iwork+1,   &
     &                        ierr)
                  CALL DLACPY('L',M,M,Work(iu),ldwrku,U,Ldu)
!
!                    Generate right bidiagonalizing vectors in WORK(IU)
!                    (Workspace: need M*M + 4*M-1,
!                                prefer M*M+3*M+(M-1)*NB)
!
                  CALL DORGBR('P',M,M,M,Work(iu),ldwrku,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in U
!                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB)
!
                  CALL DORGBR('Q',M,M,M,U,Ldu,Work(itauq),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  iwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of L in U and computing right
!                    singular vectors of L in WORK(IU)
!                    (Workspace: need M*M + BDSPAC)
!
                  CALL DBDSQR('U',M,M,M,0,S,Work(ie),Work(iu),ldwrku,U, &
     &                        Ldu,dum,1,Work(iwork),Info)
!
!                    Multiply right singular vectors of L in WORK(IU) by
!                    Q in A, storing result in VT
!                    (Workspace: need M*M)
!
                  CALL DGEMM('N','N',M,N,M,ONE,Work(iu),ldwrku,A,Lda,   &
     &                       ZERO,Vt,Ldvt)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (Workspace: need 2*M, prefer M + M*NB)
!
                  CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (Workspace: need 2*M, prefer M + M*NB)
!
                  CALL DORGLQ(M,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to U, zeroing out above it
!
                  CALL DLACPY('L',M,M,A,Lda,U,Ldu)
                  CALL DLASET('U',M-1,M-1,ZERO,ZERO,U(1,2),Ldu)
                  ie = itau
                  itauq = ie + M
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in U
!                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!
                  CALL DGEBRD(M,M,U,Ldu,S,Work(ie),Work(itauq),         &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply right bidiagonalizing vectors in U by Q
!                    in VT
!                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
!
                  CALL DORMBR('P','L','T',M,N,M,U,Ldu,Work(itaup),Vt,   &
     &                        Ldvt,Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in U
!                    (Workspace: need 4*M, prefer 3*M + M*NB)
!
                  CALL DORGBR('Q',M,M,M,U,Ldu,Work(itauq),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  iwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U and computing right
!                    singular vectors of A in VT
!                    (Workspace: need BDSPAC)
!
                  CALL DBDSQR('U',M,N,M,0,S,Work(ie),Vt,Ldvt,U,Ldu,dum, &
     &                        1,Work(iwork),Info)
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
               IF ( Lwork>=M*M+MAX(N+M,4*M,bdspac) ) THEN
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
!                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!
                  CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Copy L to WORK(IR), zeroing out above it
!
                  CALL DLACPY('L',M,M,A,Lda,Work(ir),ldwrkr)
                  CALL DLASET('U',M-1,M-1,ZERO,ZERO,Work(ir+ldwrkr),    &
     &                        ldwrkr)
!
!                    Generate Q in VT
!                    (Workspace: need M*M + M + N, prefer M*M + M + N*NB)
!
                  CALL DORGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + M
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in WORK(IR)
!                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
!
                  CALL DGEBRD(M,M,Work(ir),ldwrkr,S,Work(ie),Work(itauq)&
     &                        ,Work(itaup),Work(iwork),Lwork-iwork+1,   &
     &                        ierr)
!
!                    Generate right bidiagonalizing vectors in WORK(IR)
!                    (Workspace: need M*M + 4*M-1,
!                                prefer M*M+3*M+(M-1)*NB)
!
                  CALL DORGBR('P',M,M,M,Work(ir),ldwrkr,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  iwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing right
!                    singular vectors of L in WORK(IR)
!                    (Workspace: need M*M + BDSPAC)
!
                  CALL DBDSQR('U',M,M,0,0,S,Work(ie),Work(ir),ldwrkr,   &
     &                        dum,1,dum,1,Work(iwork),Info)
!
!                    Multiply right singular vectors of L in WORK(IR) by
!                    Q in VT, storing result in A
!                    (Workspace: need M*M)
!
                  CALL DGEMM('N','N',M,N,M,ONE,Work(ir),ldwrkr,Vt,Ldvt, &
     &                       ZERO,A,Lda)
!
!                    Copy right singular vectors of A from A to VT
!
                  CALL DLACPY('F',M,N,A,Lda,Vt,Ldvt)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (Workspace: need 2*M, prefer M + M*NB)
!
                  CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (Workspace: need M + N, prefer M + N*NB)
!
                  CALL DORGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + M
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Zero out above L in A
!
                  CALL DLASET('U',M-1,M-1,ZERO,ZERO,A(1,2),Lda)
!
!                    Bidiagonalize L in A
!                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!
                  CALL DGEBRD(M,M,A,Lda,S,Work(ie),Work(itauq),         &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply right bidiagonalizing vectors in A by Q
!                    in VT
!                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
!
                  CALL DORMBR('P','L','T',M,N,M,A,Lda,Work(itaup),Vt,   &
     &                        Ldvt,Work(iwork),Lwork-iwork+1,ierr)
                  iwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing right
!                    singular vectors of A in VT
!                    (Workspace: need BDSPAC)
!
                  CALL DBDSQR('U',M,N,0,0,S,Work(ie),Vt,Ldvt,dum,1,dum, &
     &                        1,Work(iwork),Info)
!
               ENDIF
!
            ELSEIF ( wntuo ) THEN
!
!                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
!                 N right singular vectors to be computed in VT and
!                 M left singular vectors to be overwritten on A
!
               IF ( Lwork>=2*M*M+MAX(N+M,4*M,bdspac) ) THEN
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
!                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB)
!
                  CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (Workspace: need 2*M*M + M + N, prefer 2*M*M + M + N*NB)
!
                  CALL DORGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to WORK(IU), zeroing out above it
!
                  CALL DLACPY('L',M,M,A,Lda,Work(iu),ldwrku)
                  CALL DLASET('U',M-1,M-1,ZERO,ZERO,Work(iu+ldwrku),    &
     &                        ldwrku)
                  ie = itau
                  itauq = ie + M
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in WORK(IU), copying result to
!                    WORK(IR)
!                    (Workspace: need 2*M*M + 4*M,
!                                prefer 2*M*M+3*M+2*M*NB)
!
                  CALL DGEBRD(M,M,Work(iu),ldwrku,S,Work(ie),Work(itauq)&
     &                        ,Work(itaup),Work(iwork),Lwork-iwork+1,   &
     &                        ierr)
                  CALL DLACPY('L',M,M,Work(iu),ldwrku,Work(ir),ldwrkr)
!
!                    Generate right bidiagonalizing vectors in WORK(IU)
!                    (Workspace: need 2*M*M + 4*M-1,
!                                prefer 2*M*M+3*M+(M-1)*NB)
!
                  CALL DORGBR('P',M,M,M,Work(iu),ldwrku,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in WORK(IR)
!                    (Workspace: need 2*M*M + 4*M, prefer 2*M*M + 3*M + M*NB)
!
                  CALL DORGBR('Q',M,M,M,Work(ir),ldwrkr,Work(itauq),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
                  iwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of L in WORK(IR) and computing
!                    right singular vectors of L in WORK(IU)
!                    (Workspace: need 2*M*M + BDSPAC)
!
                  CALL DBDSQR('U',M,M,M,0,S,Work(ie),Work(iu),ldwrku,   &
     &                        Work(ir),ldwrkr,dum,1,Work(iwork),Info)
!
!                    Multiply right singular vectors of L in WORK(IU) by
!                    Q in VT, storing result in A
!                    (Workspace: need M*M)
!
                  CALL DGEMM('N','N',M,N,M,ONE,Work(iu),ldwrku,Vt,Ldvt, &
     &                       ZERO,A,Lda)
!
!                    Copy right singular vectors of A from A to VT
!
                  CALL DLACPY('F',M,N,A,Lda,Vt,Ldvt)
!
!                    Copy left singular vectors of A from WORK(IR) to A
!
                  CALL DLACPY('F',M,M,Work(ir),ldwrkr,A,Lda)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (Workspace: need 2*M, prefer M + M*NB)
!
                  CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (Workspace: need M + N, prefer M + N*NB)
!
                  CALL DORGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
                  ie = itau
                  itauq = ie + M
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Zero out above L in A
!
                  CALL DLASET('U',M-1,M-1,ZERO,ZERO,A(1,2),Lda)
!
!                    Bidiagonalize L in A
!                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!
                  CALL DGEBRD(M,M,A,Lda,S,Work(ie),Work(itauq),         &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply right bidiagonalizing vectors in A by Q
!                    in VT
!                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
!
                  CALL DORMBR('P','L','T',M,N,M,A,Lda,Work(itaup),Vt,   &
     &                        Ldvt,Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in A
!                    (Workspace: need 4*M, prefer 3*M + M*NB)
!
                  CALL DORGBR('Q',M,M,M,A,Lda,Work(itauq),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  iwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in A and computing right
!                    singular vectors of A in VT
!                    (Workspace: need BDSPAC)
!
                  CALL DBDSQR('U',M,N,M,0,S,Work(ie),Vt,Ldvt,A,Lda,dum, &
     &                        1,Work(iwork),Info)
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
               IF ( Lwork>=M*M+MAX(N+M,4*M,bdspac) ) THEN
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
!                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!
                  CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (Workspace: need M*M + M + N, prefer M*M + M + N*NB)
!
                  CALL DORGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to WORK(IU), zeroing out above it
!
                  CALL DLACPY('L',M,M,A,Lda,Work(iu),ldwrku)
                  CALL DLASET('U',M-1,M-1,ZERO,ZERO,Work(iu+ldwrku),    &
     &                        ldwrku)
                  ie = itau
                  itauq = ie + M
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in WORK(IU), copying result to U
!                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
!
                  CALL DGEBRD(M,M,Work(iu),ldwrku,S,Work(ie),Work(itauq)&
     &                        ,Work(itaup),Work(iwork),Lwork-iwork+1,   &
     &                        ierr)
                  CALL DLACPY('L',M,M,Work(iu),ldwrku,U,Ldu)
!
!                    Generate right bidiagonalizing vectors in WORK(IU)
!                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + (M-1)*NB)
!
                  CALL DORGBR('P',M,M,M,Work(iu),ldwrku,Work(itaup),    &
     &                        Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in U
!                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB)
!
                  CALL DORGBR('Q',M,M,M,U,Ldu,Work(itauq),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  iwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of L in U and computing right
!                    singular vectors of L in WORK(IU)
!                    (Workspace: need M*M + BDSPAC)
!
                  CALL DBDSQR('U',M,M,M,0,S,Work(ie),Work(iu),ldwrku,U, &
     &                        Ldu,dum,1,Work(iwork),Info)
!
!                    Multiply right singular vectors of L in WORK(IU) by
!                    Q in VT, storing result in A
!                    (Workspace: need M*M)
!
                  CALL DGEMM('N','N',M,N,M,ONE,Work(iu),ldwrku,Vt,Ldvt, &
     &                       ZERO,A,Lda)
!
!                    Copy right singular vectors of A from A to VT
!
                  CALL DLACPY('F',M,N,A,Lda,Vt,Ldvt)
!
               ELSE
!
!                    Insufficient workspace for a fast algorithm
!
                  itau = 1
                  iwork = itau + M
!
!                    Compute A=L*Q, copying result to VT
!                    (Workspace: need 2*M, prefer M + M*NB)
!
                  CALL DGELQF(M,N,A,Lda,Work(itau),Work(iwork),         &
     &                        Lwork-iwork+1,ierr)
                  CALL DLACPY('U',M,N,A,Lda,Vt,Ldvt)
!
!                    Generate Q in VT
!                    (Workspace: need M + N, prefer M + N*NB)
!
                  CALL DORGLQ(N,N,M,Vt,Ldvt,Work(itau),Work(iwork),     &
     &                        Lwork-iwork+1,ierr)
!
!                    Copy L to U, zeroing out above it
!
                  CALL DLACPY('L',M,M,A,Lda,U,Ldu)
                  CALL DLASET('U',M-1,M-1,ZERO,ZERO,U(1,2),Ldu)
                  ie = itau
                  itauq = ie + M
                  itaup = itauq + M
                  iwork = itaup + M
!
!                    Bidiagonalize L in U
!                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!
                  CALL DGEBRD(M,M,U,Ldu,S,Work(ie),Work(itauq),         &
     &                        Work(itaup),Work(iwork),Lwork-iwork+1,    &
     &                        ierr)
!
!                    Multiply right bidiagonalizing vectors in U by Q
!                    in VT
!                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
!
                  CALL DORMBR('P','L','T',M,N,M,U,Ldu,Work(itaup),Vt,   &
     &                        Ldvt,Work(iwork),Lwork-iwork+1,ierr)
!
!                    Generate left bidiagonalizing vectors in U
!                    (Workspace: need 4*M, prefer 3*M + M*NB)
!
                  CALL DORGBR('Q',M,M,M,U,Ldu,Work(itauq),Work(iwork),  &
     &                        Lwork-iwork+1,ierr)
                  iwork = ie + M
!
!                    Perform bidiagonal QR iteration, computing left
!                    singular vectors of A in U and computing right
!                    singular vectors of A in VT
!                    (Workspace: need BDSPAC)
!
                  CALL DBDSQR('U',M,N,M,0,S,Work(ie),Vt,Ldvt,U,Ldu,dum, &
     &                        1,Work(iwork),Info)
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
         itauq = ie + M
         itaup = itauq + M
         iwork = itaup + M
!
!           Bidiagonalize A
!           (Workspace: need 3*M + N, prefer 3*M + (M + N)*NB)
!
         CALL DGEBRD(M,N,A,Lda,S,Work(ie),Work(itauq),Work(itaup),      &
     &               Work(iwork),Lwork-iwork+1,ierr)
         IF ( wntuas ) THEN
!
!              If left singular vectors desired in U, copy result to U
!              and generate left bidiagonalizing vectors in U
!              (Workspace: need 4*M-1, prefer 3*M + (M-1)*NB)
!
            CALL DLACPY('L',M,M,A,Lda,U,Ldu)
            CALL DORGBR('Q',M,M,N,U,Ldu,Work(itauq),Work(iwork),        &
     &                  Lwork-iwork+1,ierr)
         ENDIF
         IF ( wntvas ) THEN
!
!              If right singular vectors desired in VT, copy result to
!              VT and generate right bidiagonalizing vectors in VT
!              (Workspace: need 3*M + NRVT, prefer 3*M + NRVT*NB)
!
            CALL DLACPY('U',M,N,A,Lda,Vt,Ldvt)
            IF ( wntva ) nrvt = N
            IF ( wntvs ) nrvt = M
            CALL DORGBR('P',nrvt,N,M,Vt,Ldvt,Work(itaup),Work(iwork),   &
     &                  Lwork-iwork+1,ierr)
         ENDIF
!
!              If left singular vectors desired in A, generate left
!              bidiagonalizing vectors in A
!              (Workspace: need 4*M-1, prefer 3*M + (M-1)*NB)
!
         IF ( wntuo ) CALL DORGBR('Q',M,M,N,A,Lda,Work(itauq),          &
     &                            Work(iwork),Lwork-iwork+1,ierr)
!
!              If right singular vectors desired in A, generate right
!              bidiagonalizing vectors in A
!              (Workspace: need 4*M, prefer 3*M + M*NB)
!
         IF ( wntvo ) CALL DORGBR('P',M,N,M,A,Lda,Work(itaup),          &
     &                            Work(iwork),Lwork-iwork+1,ierr)
         iwork = ie + M
         IF ( wntuas .OR. wntuo ) nru = M
         IF ( wntun ) nru = 0
         IF ( wntvas .OR. wntvo ) ncvt = N
         IF ( wntvn ) ncvt = 0
         IF ( (.NOT.wntuo) .AND. (.NOT.wntvo) ) THEN
!
!              Perform bidiagonal QR iteration, if desired, computing
!              left singular vectors in U and computing right singular
!              vectors in VT
!              (Workspace: need BDSPAC)
!
            CALL DBDSQR('L',M,ncvt,nru,0,S,Work(ie),Vt,Ldvt,U,Ldu,dum,1,&
     &                  Work(iwork),Info)
         ELSEIF ( (.NOT.wntuo) .AND. wntvo ) THEN
!
!              Perform bidiagonal QR iteration, if desired, computing
!              left singular vectors in U and computing right singular
!              vectors in A
!              (Workspace: need BDSPAC)
!
            CALL DBDSQR('L',M,ncvt,nru,0,S,Work(ie),A,Lda,U,Ldu,dum,1,  &
     &                  Work(iwork),Info)
         ELSE
!
!              Perform bidiagonal QR iteration, if desired, computing
!              left singular vectors in A and computing right singular
!              vectors in VT
!              (Workspace: need BDSPAC)
!
            CALL DBDSQR('L',M,ncvt,nru,0,S,Work(ie),Vt,Ldvt,A,Lda,dum,1,&
     &                  Work(iwork),Info)
         ENDIF
!
!
      ENDIF
!
!     If DBDSQR failed to converge, copy unconverged superdiagonals
!     to WORK( 2:MINMN )
!
      IF ( Info/=0 ) THEN
         IF ( ie>2 ) THEN
            DO i = 1 , minmn - 1
               Work(i+1) = Work(i+ie-1)
            ENDDO
         ENDIF
         IF ( ie<2 ) THEN
            DO i = minmn - 1 , 1 , -1
               Work(i+1) = Work(i+ie-1)
            ENDDO
         ENDIF
      ENDIF
!
!     Undo scaling if necessary
!
      IF ( iscl==1 ) THEN
         IF ( anrm>bignum ) CALL DLASCL('G',0,0,bignum,anrm,minmn,1,S,  &
     &                                  minmn,ierr)
         IF ( Info/=0 .AND. anrm>bignum )                               &
     &        CALL DLASCL('G',0,0,bignum,anrm,minmn-1,1,Work(2),minmn,  &
     &        ierr)
         IF ( anrm<smlnum ) CALL DLASCL('G',0,0,smlnum,anrm,minmn,1,S,  &
     &                                  minmn,ierr)
         IF ( Info/=0 .AND. anrm<smlnum )                               &
     &        CALL DLASCL('G',0,0,smlnum,anrm,minmn-1,1,Work(2),minmn,  &
     &        ierr)
      ENDIF
!
!     Return optimal workspace in WORK(1)
!
      Work(1) = maxwrk
!
!
!     End of DGESVD
!
      END SUBROUTINE DGESVD
