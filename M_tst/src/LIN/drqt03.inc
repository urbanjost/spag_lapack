!*==drqt03.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b drqt03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE drqt03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK,
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
!> drqt03 tests DORMRQ, which computes Q*C, Q'*C, C*Q or C*Q'.
!>
!> drqt03 compares the results of a call to DORMRQ with the results of
!> forming Q explicitly by a call to DORGRQ and then performing matrix
!> multiplication by a call to DGEMM.
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
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!>          Details of the RQ factorization of an m-by-n matrix, as
!>          returned by DGERQF. See SGERQF for further details.
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
!>          Q is DOUBLE PRECISION array, dimension (LDA,N)
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
!>          to the RQ factorization in AF.
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE drqt03(M,N,K,Af,C,Cc,Q,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--DRQT03139
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
EXTERNAL    LSAME,   DLAMCH,   DLANGE
!     ..
!     .. External Subroutines ..
EXTERNAL    DGEMM,   DLACPY,   DLARNV,   DLASET,   DORGRQ,   DORMRQ
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
!     Copy the last k rows of the factorization to the array Q
!
      CALL DLASET('Full',N,N,ROGUE,ROGUE,Q,Lda)
      IF ( K>0 .AND. N>K ) CALL DLACPY('Full',K,N-K,Af(M-K+1,1),Lda,    &
     &                                 Q(N-K+1,1),Lda)
      IF ( K>1 ) CALL DLACPY('Lower',K-1,K-1,Af(M-K+2,N-K+1),Lda,       &
     &                       Q(N-K+2,N-K+1),Lda)
!
!     Generate the n-by-n matrix Q
!
      SRNamt = 'DORGRQ'
      CALL DORGRQ(N,N,K,Q,Lda,Tau(minmn-K+1),Work,Lwork,info)
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
            SRNamt = 'DORMRQ'
            IF ( K>0 ) CALL DORMRQ(side,trans,mc,nc,K,Af(M-K+1,1),Lda,  &
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
     &         = resid/(DBLE(MAX(1,N))*cnorm*eps)
!
         ENDDO
      ENDDO
!
!
!     End of drqt03
!
      END SUBROUTINE drqt03
