!*==sqrt01p.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SQRT01P
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SQRT01P( M, N, A, AF, Q, R, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), AF( LDA, * ), Q( LDA, * ),
!      $                   R( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SQRT01P tests SGEQRFP, which computes the QR factorization of an m-by-n
!> matrix A, and partially tests SORGQR which forms the m-by-m
!> orthogonal matrix Q.
!>
!> SQRT01P compares R with Q'*A, and checks that Q is orthogonal.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The m-by-n matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is REAL array, dimension (LDA,N)
!>          Details of the QR factorization of A, as returned by SGEQRFP.
!>          See SGEQRFP for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDA,M)
!>          The m-by-m orthogonal matrix Q.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is REAL array, dimension (LDA,max(M,N))
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, Q and R.
!>          LDA >= max(M,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors, as returned
!>          by SGEQRFP.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
!>          The test ratios:
!>          RESULT(1) = norm( R - Q'*A ) / ( M * norm(A) * EPS )
!>          RESULT(2) = norm( I - Q'*Q ) / ( M * EPS )
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE SQRT01P(M,N,A,Af,Q,R,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--SQRT01P129
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Lwork , M , N
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , Af(Lda,*) , Q(Lda,*) , R(Lda,*) , Result(*) ,     &
     &     Rwork(*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      REAL ROGUE
      PARAMETER (ROGUE=-1.0E+10)
!     ..
!     .. Local Scalars ..
      INTEGER info , minmn
      REAL anorm , eps , resid
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE , SLANSY
      EXTERNAL SLAMCH , SLANGE , SLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEMM , SGEQRFP , SLACPY , SLASET , SORGQR , SSYRK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , REAL
!     ..
!     .. Scalars in Common ..
      CHARACTER*32 SRNamt
!     ..
!     .. Common blocks ..
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Executable Statements ..
!
      minmn = MIN(M,N)
      eps = SLAMCH('Epsilon')
!
!     Copy the matrix A to the array AF.
!
      CALL SLACPY('Full',M,N,A,Lda,Af,Lda)
!
!     Factorize the matrix A in the array AF.
!
      SRNamt = 'SGEQRFP'
      CALL SGEQRFP(M,N,Af,Lda,Tau,Work,Lwork,info)
!
!     Copy details of Q
!
      CALL SLASET('Full',M,M,ROGUE,ROGUE,Q,Lda)
      CALL SLACPY('Lower',M-1,N,Af(2,1),Lda,Q(2,1),Lda)
!
!     Generate the m-by-m matrix Q
!
      SRNamt = 'SORGQR'
      CALL SORGQR(M,M,minmn,Q,Lda,Tau,Work,Lwork,info)
!
!     Copy R
!
      CALL SLASET('Full',M,N,ZERO,ZERO,R,Lda)
      CALL SLACPY('Upper',M,N,Af,Lda,R,Lda)
!
!     Compute R - Q'*A
!
      CALL SGEMM('Transpose','No transpose',M,N,M,-ONE,Q,Lda,A,Lda,ONE, &
     &           R,Lda)
!
!     Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .
!
      anorm = SLANGE('1',M,N,A,Lda,Rwork)
      resid = SLANGE('1',M,N,R,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/REAL(MAX(1,M)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q'*Q
!
      CALL SLASET('Full',M,M,ZERO,ONE,R,Lda)
      CALL SSYRK('Upper','Transpose',M,M,-ONE,Q,Lda,ONE,R,Lda)
!
!     Compute norm( I - Q'*Q ) / ( M * EPS ) .
!
      resid = SLANSY('1','Upper',M,R,Lda,Rwork)
!
      Result(2) = (resid/REAL(MAX(1,M)))/eps
!
!
!     End of SQRT01P
!
      END SUBROUTINE SQRT01P
