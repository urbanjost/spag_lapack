!*==srqt03.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SRQT03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SRQT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               AF( LDA, * ), C( LDA, * ), CC( LDA, * ),
!      $                   Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SRQT03 tests SORMRQ, which computes Q*C, Q'*C, C*Q or C*Q'.
!>
!> SRQT03 compares the results of a call to SORMRQ with the results of
!> forming Q explicitly by a call to SORGRQ and then performing matrix
!> multiplication by a call to SGEMM.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows or columns of the matrix C; C is n-by-m if
!>          Q is applied from the left, or m-by-n if Q is applied from
!>          the right.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the orthogonal matrix Q.  N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          orthogonal matrix Q.  N >= K >= 0.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is REAL array, dimension (LDA,N)
!>          Details of the RQ factorization of an m-by-n matrix, as
!>          returned by SGERQF. See SGERQF for further details.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] CC
!> \verbatim
!>          CC is REAL array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays AF, C, CC, and Q.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors corresponding
!>          to the RQ factorization in AF.
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
!>          The length of WORK.  LWORK must be at least M, and should be
!>          M*NB, where NB is the blocksize for this environment.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (4)
!>          The test ratios compare two techniques for multiplying a
!>          random matrix C by an n-by-n orthogonal matrix Q.
!>          RESULT(1) = norm( Q*C - Q*C )  / ( N * norm(C) * EPS )
!>          RESULT(2) = norm( C*Q - C*Q )  / ( N * norm(C) * EPS )
!>          RESULT(3) = norm( Q'*C - Q'*C )/ ( N * norm(C) * EPS )
!>          RESULT(4) = norm( C*Q' - C*Q' )/ ( N * norm(C) * EPS )
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
      SUBROUTINE SRQT03(M,N,K,Af,C,Cc,Q,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--SRQT03139
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER K , Lda , Lwork , M , N
!     ..
!     .. Array Arguments ..
      REAL Af(Lda,*) , C(Lda,*) , Cc(Lda,*) , Q(Lda,*) , Result(*) ,    &
     &     Rwork(*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
      REAL ROGUE
      PARAMETER (ROGUE=-1.0E+10)
!     ..
!     .. Local Scalars ..
      CHARACTER side , trans
      INTEGER info , iside , itrans , j , mc , minmn , nc
      REAL cnorm , eps , resid
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLAMCH , SLANGE
      EXTERNAL LSAME , SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEMM , SLACPY , SLARNV , SLASET , SORGRQ , SORMRQ
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
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
!     .. Data statements ..
      DATA iseed/1988 , 1989 , 1990 , 1991/
!     ..
!     .. Executable Statements ..
!
      eps = SLAMCH('Epsilon')
      minmn = MIN(M,N)
!
!     Quick return if possible
!
      IF ( minmn==0 ) THEN
         Result(1) = ZERO
         Result(2) = ZERO
         Result(3) = ZERO
         Result(4) = ZERO
         RETURN
      ENDIF
!
!     Copy the last k rows of the factorization to the array Q
!
      CALL SLASET('Full',N,N,ROGUE,ROGUE,Q,Lda)
      IF ( K>0 .AND. N>K ) CALL SLACPY('Full',K,N-K,Af(M-K+1,1),Lda,    &
     &                                 Q(N-K+1,1),Lda)
      IF ( K>1 ) CALL SLACPY('Lower',K-1,K-1,Af(M-K+2,N-K+1),Lda,       &
     &                       Q(N-K+2,N-K+1),Lda)
!
!     Generate the n-by-n matrix Q
!
      SRNamt = 'SORGRQ'
      CALL SORGRQ(N,N,K,Q,Lda,Tau(minmn-K+1),Work,Lwork,info)
!
      DO iside = 1 , 2
         IF ( iside==1 ) THEN
            side = 'L'
            mc = N
            nc = M
         ELSE
            side = 'R'
            mc = M
            nc = N
         ENDIF
!
!        Generate MC by NC matrix C
!
         DO j = 1 , nc
            CALL SLARNV(2,iseed,mc,C(1,j))
         ENDDO
         cnorm = SLANGE('1',mc,nc,C,Lda,Rwork)
         IF ( cnorm==0.0 ) cnorm = ONE
!
         DO itrans = 1 , 2
            IF ( itrans==1 ) THEN
               trans = 'N'
            ELSE
               trans = 'T'
            ENDIF
!
!           Copy C
!
            CALL SLACPY('Full',mc,nc,C,Lda,Cc,Lda)
!
!           Apply Q or Q' to C
!
            SRNamt = 'SORMRQ'
            IF ( K>0 ) CALL SORMRQ(side,trans,mc,nc,K,Af(M-K+1,1),Lda,  &
     &                             Tau(minmn-K+1),Cc,Lda,Work,Lwork,    &
     &                             info)
!
!           Form explicit product and subtract
!
            IF ( LSAME(side,'L') ) THEN
               CALL SGEMM(trans,'No transpose',mc,nc,mc,-ONE,Q,Lda,C,   &
     &                    Lda,ONE,Cc,Lda)
            ELSE
               CALL SGEMM('No transpose',trans,mc,nc,nc,-ONE,C,Lda,Q,   &
     &                    Lda,ONE,Cc,Lda)
            ENDIF
!
!           Compute error in the difference
!
            resid = SLANGE('1',mc,nc,Cc,Lda,Rwork)
            Result((iside-1)*2+itrans)                                  &
     &         = resid/(REAL(MAX(1,N))*cnorm*eps)
!
         ENDDO
      ENDDO
!
!
!     End of SRQT03
!
      END SUBROUTINE SRQT03
