!*==dqrt14.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b dqrt14
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION dqrt14( TRANS, M, N, NRHS, A, LDA, X,
!                        LDX, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            LDA, LDX, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), WORK( LWORK ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dqrt14 checks whether X is in the row space of A or A'.  It does so
!> by scaling both X and A such that their norms are in the range
!> [sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of [A,X]
!> (if TRANS = 'T') or an LQ factorization of [A',X]' (if TRANS = 'N'),
!> and returning the norm of the trailing triangle, scaled by
!> MAX(M,N,NRHS)*eps.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, check for X in the row space of A
!>          = 'T':  Transpose, check for X in the row space of A'.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of X.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
!>          If TRANS = 'N', the N-by-NRHS matrix X.
!>          IF TRANS = 'T', the M-by-NRHS matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          length of workspace array required
!>          If TRANS = 'N', LWORK >= (M+NRHS)*(N+2);
!>          if TRANS = 'T', LWORK >= (N+NRHS)*(M+2).
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
!> \ingroup double_lin
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION dqrt14(Trans,M,N,Nrhs,A,Lda,X,Ldx,Work, &
     &                                 Lwork)
      IMPLICIT NONE
!*--DQRT14120
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Lda , Ldx , Lwork , M , N , Nrhs
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , Work(Lwork) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      LOGICAL tpsd
      INTEGER i , info , j , ldwork
      DOUBLE PRECISION anrm , err , xnrm
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION rwork(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , DLANGE
EXTERNAL    LSAME,   DLAMCH,   DLANGE
!     ..
!     .. External Subroutines ..
EXTERNAL    DGELQ2,   DGEQR2,   DLACPY,   DLASCL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN
!     ..
!     .. Executable Statements ..
!
      dqrt14 = ZERO
      IF ( LSAME(Trans,'N') ) THEN
         ldwork = M + Nrhs
         tpsd = .FALSE.
         IF ( Lwork<(M+Nrhs)*(N+2) ) THEN
            CALL xerbla('dqrt14',10)
            RETURN
         ELSEIF ( N<=0 .OR. Nrhs<=0 ) THEN
            RETURN
         ENDIF
      ELSEIF ( LSAME(Trans,'T') ) THEN
         ldwork = M
         tpsd = .TRUE.
         IF ( Lwork<(N+Nrhs)*(M+2) ) THEN
            CALL xerbla('dqrt14',10)
            RETURN
         ELSEIF ( M<=0 .OR. Nrhs<=0 ) THEN
            RETURN
         ENDIF
      ELSE
         CALL xerbla('dqrt14',1)
         RETURN
      ENDIF
!
!     Copy and scale A
!
      CALL DLACPY('All',M,N,A,Lda,Work,ldwork)
      anrm = DLANGE('M',M,N,Work,ldwork,rwork)
      IF ( anrm/=ZERO ) CALL DLASCL('G',0,0,anrm,ONE,M,N,Work,ldwork,   &
     &                              info)
!
!     Copy X or X' into the right place and scale it
!
      IF ( tpsd ) THEN
!
!        Copy X into columns n+1:n+nrhs of work
!
         CALL DLACPY('All',M,Nrhs,X,Ldx,Work(N*ldwork+1),ldwork)
         xnrm = DLANGE('M',M,Nrhs,Work(N*ldwork+1),ldwork,rwork)
         IF ( xnrm/=ZERO ) CALL DLASCL('G',0,0,xnrm,ONE,M,Nrhs,         &
     &                                 Work(N*ldwork+1),ldwork,info)
         anrm = DLANGE('One-norm',M,N+Nrhs,Work,ldwork,rwork)
!
!        Compute QR factorization of X
!
         CALL DGEQR2(M,N+Nrhs,Work,ldwork,Work(ldwork*(N+Nrhs)+1),      &
     &               Work(ldwork*(N+Nrhs)+MIN(M,N+Nrhs)+1),info)
!
!        Compute largest entry in upper triangle of
!        work(n+1:m,n+1:n+nrhs)
!
         err = ZERO
         DO j = N + 1 , N + Nrhs
            DO i = N + 1 , MIN(M,j)
               err = MAX(err,ABS(Work(i+(j-1)*M)))
            ENDDO
         ENDDO
!
      ELSE
!
!        Copy X' into rows m+1:m+nrhs of work
!
         DO i = 1 , N
            DO j = 1 , Nrhs
               Work(M+j+(i-1)*ldwork) = X(i,j)
            ENDDO
         ENDDO
!
         xnrm = DLANGE('M',Nrhs,N,Work(M+1),ldwork,rwork)
         IF ( xnrm/=ZERO ) CALL DLASCL('G',0,0,xnrm,ONE,Nrhs,N,Work(M+1)&
     &                                 ,ldwork,info)
!
!        Compute LQ factorization of work
!
         CALL DGELQ2(ldwork,N,Work,ldwork,Work(ldwork*N+1),             &
     &               Work(ldwork*(N+1)+1),info)
!
!        Compute largest entry in lower triangle in
!        work(m+1:m+nrhs,m+1:n)
!
         err = ZERO
         DO j = M + 1 , N
            DO i = j , ldwork
               err = MAX(err,ABS(Work(i+(j-1)*ldwork)))
            ENDDO
         ENDDO
!
      ENDIF
!
      dqrt14 = err/(DBLE(MAX(M,N,Nrhs))*DLAMCH('Epsilon'))
!
!
!     End of dqrt14
!
      END FUNCTION dqrt14
