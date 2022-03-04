!*==dqlt03.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DQLT03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DQLT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AF( LDA, * ), C( LDA, * ), CC( LDA, * ),
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
!> DQLT03 tests DORMQL, which computes Q*C, Q'*C, C*Q or C*Q'.
!>
!> DQLT03 compares the results of a call to DORMQL with the results of
!> forming Q explicitly by a call to DORGQL and then performing matrix
!> multiplication by a call to DGEMM.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the orthogonal matrix Q.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows or columns of the matrix C; C is m-by-n if
!>          Q is applied from the left, or n-by-m if Q is applied from
!>          the right.  N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          orthogonal matrix Q.  M >= K >= 0.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!>          Details of the QL factorization of an m-by-n matrix, as
!>          returned by DGEQLF. See SGEQLF for further details.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] CC
!> \verbatim
!>          CC is DOUBLE PRECISION array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDA,M)
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
!>          TAU is DOUBLE PRECISION array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors corresponding
!>          to the QL factorization in AF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
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
!>          RWORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (4)
!>          The test ratios compare two techniques for multiplying a
!>          random matrix C by an m-by-m orthogonal matrix Q.
!>          RESULT(1) = norm( Q*C - Q*C )  / ( M * norm(C) * EPS )
!>          RESULT(2) = norm( C*Q - C*Q )  / ( M * norm(C) * EPS )
!>          RESULT(3) = norm( Q'*C - Q'*C )/ ( M * norm(C) * EPS )
!>          RESULT(4) = norm( C*Q' - C*Q' )/ ( M * norm(C) * EPS )
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
      SUBROUTINE DQLT03(M,N,K,Af,C,Cc,Q,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--DQLT03139
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
      DOUBLE PRECISION Af(Lda,*) , C(Lda,*) , Cc(Lda,*) , Q(Lda,*) ,    &
     &                 Result(*) , Rwork(*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION ROGUE
      PARAMETER (ROGUE=-1.0D+10)
!     ..
!     .. Local Scalars ..
      CHARACTER side , trans
      INTEGER info , iside , itrans , j , mc , minmn , nc
      DOUBLE PRECISION cnorm , eps , resid
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , DLANGE
      EXTERNAL LSAME , DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DLACPY , DLARNV , DLASET , DORGQL , DORMQL
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX , MIN
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
      eps = DLAMCH('Epsilon')
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
!     Copy the last k columns of the factorization to the array Q
!
      CALL DLASET('Full',M,M,ROGUE,ROGUE,Q,Lda)
      IF ( K>0 .AND. M>K ) CALL DLACPY('Full',M-K,K,Af(1,N-K+1),Lda,    &
     &                                 Q(1,M-K+1),Lda)
      IF ( K>1 ) CALL DLACPY('Upper',K-1,K-1,Af(M-K+1,N-K+2),Lda,       &
     &                       Q(M-K+1,M-K+2),Lda)
!
!     Generate the m-by-m matrix Q
!
      SRNamt = 'DORGQL'
      CALL DORGQL(M,M,K,Q,Lda,Tau(minmn-K+1),Work,Lwork,info)
!
      DO iside = 1 , 2
         IF ( iside==1 ) THEN
            side = 'L'
            mc = M
            nc = N
         ELSE
            side = 'R'
            mc = N
            nc = M
         ENDIF
!
!        Generate MC by NC matrix C
!
         DO j = 1 , nc
            CALL DLARNV(2,iseed,mc,C(1,j))
         ENDDO
         cnorm = DLANGE('1',mc,nc,C,Lda,Rwork)
         IF ( cnorm==0.0D0 ) cnorm = ONE
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
            CALL DLACPY('Full',mc,nc,C,Lda,Cc,Lda)
!
!           Apply Q or Q' to C
!
            SRNamt = 'DORMQL'
            IF ( K>0 ) CALL DORMQL(side,trans,mc,nc,K,Af(1,N-K+1),Lda,  &
     &                             Tau(minmn-K+1),Cc,Lda,Work,Lwork,    &
     &                             info)
!
!           Form explicit product and subtract
!
            IF ( LSAME(side,'L') ) THEN
               CALL DGEMM(trans,'No transpose',mc,nc,mc,-ONE,Q,Lda,C,   &
     &                    Lda,ONE,Cc,Lda)
            ELSE
               CALL DGEMM('No transpose',trans,mc,nc,nc,-ONE,C,Lda,Q,   &
     &                    Lda,ONE,Cc,Lda)
            ENDIF
!
!           Compute error in the difference
!
            resid = DLANGE('1',mc,nc,Cc,Lda,Rwork)
            Result((iside-1)*2+itrans)                                  &
     &         = resid/(DBLE(MAX(1,M))*cnorm*eps)
!
         ENDDO
      ENDDO
!
!
!     End of DQLT03
!
      END SUBROUTINE DQLT03
