!*==zqlt03.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZQLT03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZQLT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RESULT( * ), RWORK( * )
!       COMPLEX*16         AF( LDA, * ), C( LDA, * ), CC( LDA, * ),
!      $                   Q( LDA, * ), TAU( * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZQLT03 tests ZUNMQL, which computes Q*C, Q'*C, C*Q or C*Q'.
!>
!> ZQLT03 compares the results of a call to ZUNMQL with the results of
!> forming Q explicitly by a call to ZUNGQL and then performing matrix
!> multiplication by a call to ZGEMM.
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
!>          AF is COMPLEX*16 array, dimension (LDA,N)
!>          Details of the QL factorization of an m-by-n matrix, as
!>          returned by ZGEQLF. See CGEQLF for further details.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] CC
!> \verbatim
!>          CC is COMPLEX*16 array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDA,M)
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
!>          TAU is COMPLEX*16 array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors corresponding
!>          to the QL factorization in AF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZQLT03(M,N,K,Af,C,Cc,Q,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--ZQLT03139
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
      DOUBLE PRECISION Result(*) , Rwork(*)
      COMPLEX*16 Af(Lda,*) , C(Lda,*) , Cc(Lda,*) , Q(Lda,*) , Tau(*) , &
     &           Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 ROGUE
      PARAMETER (ROGUE=(-1.0D+10,-1.0D+10))
!     ..
!     .. Local Scalars ..
      CHARACTER side , trans
      INTEGER info , iside , itrans , j , mc , minmn , nc
      DOUBLE PRECISION cnorm , eps , resid
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL LSAME , DLAMCH , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMM , ZLACPY , ZLARNV , ZLASET , ZUNGQL , ZUNMQL
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCMPLX , MAX , MIN
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
      CALL ZLASET('Full',M,M,ROGUE,ROGUE,Q,Lda)
      IF ( K>0 .AND. M>K ) CALL ZLACPY('Full',M-K,K,Af(1,N-K+1),Lda,    &
     &                                 Q(1,M-K+1),Lda)
      IF ( K>1 ) CALL ZLACPY('Upper',K-1,K-1,Af(M-K+1,N-K+2),Lda,       &
     &                       Q(M-K+1,M-K+2),Lda)
!
!     Generate the m-by-m matrix Q
!
      SRNamt = 'ZUNGQL'
      CALL ZUNGQL(M,M,K,Q,Lda,Tau(minmn-K+1),Work,Lwork,info)
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
            CALL ZLARNV(2,iseed,mc,C(1,j))
         ENDDO
         cnorm = ZLANGE('1',mc,nc,C,Lda,Rwork)
         IF ( cnorm==ZERO ) cnorm = ONE
!
         DO itrans = 1 , 2
            IF ( itrans==1 ) THEN
               trans = 'N'
            ELSE
               trans = 'C'
            ENDIF
!
!           Copy C
!
            CALL ZLACPY('Full',mc,nc,C,Lda,Cc,Lda)
!
!           Apply Q or Q' to C
!
            SRNamt = 'ZUNMQL'
            IF ( K>0 ) CALL ZUNMQL(side,trans,mc,nc,K,Af(1,N-K+1),Lda,  &
     &                             Tau(minmn-K+1),Cc,Lda,Work,Lwork,    &
     &                             info)
!
!           Form explicit product and subtract
!
            IF ( LSAME(side,'L') ) THEN
               CALL ZGEMM(trans,'No transpose',mc,nc,mc,DCMPLX(-ONE),Q, &
     &                    Lda,C,Lda,DCMPLX(ONE),Cc,Lda)
            ELSE
               CALL ZGEMM('No transpose',trans,mc,nc,nc,DCMPLX(-ONE),C, &
     &                    Lda,Q,Lda,DCMPLX(ONE),Cc,Lda)
            ENDIF
!
!           Compute error in the difference
!
            resid = ZLANGE('1',mc,nc,Cc,Lda,Rwork)
            Result((iside-1)*2+itrans)                                  &
     &         = resid/(DBLE(MAX(1,M))*cnorm*eps)
!
         ENDDO
      ENDDO
!
!
!     End of ZQLT03
!
      END SUBROUTINE ZQLT03
