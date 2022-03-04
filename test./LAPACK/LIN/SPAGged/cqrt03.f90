!*==cqrt03.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CQRT03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CQRT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               RESULT( * ), RWORK( * )
!       COMPLEX            AF( LDA, * ), C( LDA, * ), CC( LDA, * ),
!      $                   Q( LDA, * ), TAU( * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CQRT03 tests CUNMQR, which computes Q*C, Q'*C, C*Q or C*Q'.
!>
!> CQRT03 compares the results of a call to CUNMQR with the results of
!> forming Q explicitly by a call to CUNGQR and then performing matrix
!> multiplication by a call to CGEMM.
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
!>          AF is COMPLEX array, dimension (LDA,N)
!>          Details of the QR factorization of an m-by-n matrix, as
!>          returned by CGEQRF. See CGEQRF for further details.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] CC
!> \verbatim
!>          CC is COMPLEX array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDA,M)
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
!>          TAU is COMPLEX array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors corresponding
!>          to the QR factorization in AF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CQRT03(M,N,K,Af,C,Cc,Q,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--CQRT03139
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
      REAL Result(*) , Rwork(*)
      COMPLEX Af(Lda,*) , C(Lda,*) , Cc(Lda,*) , Q(Lda,*) , Tau(*) ,    &
     &        Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX ROGUE
      PARAMETER (ROGUE=(-1.0E+10,-1.0E+10))
!     ..
!     .. Local Scalars ..
      CHARACTER side , trans
      INTEGER info , iside , itrans , j , mc , nc
      REAL cnorm , eps , resid
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANGE , SLAMCH
      EXTERNAL LSAME , CLANGE , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMM , CLACPY , CLARNV , CLASET , CUNGQR , CUNMQR
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX , REAL
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
!
!     Copy the first k columns of the factorization to the array Q
!
      CALL CLASET('Full',M,M,ROGUE,ROGUE,Q,Lda)
      CALL CLACPY('Lower',M-1,K,Af(2,1),Lda,Q(2,1),Lda)
!
!     Generate the m-by-m matrix Q
!
      SRNamt = 'CUNGQR'
      CALL CUNGQR(M,M,K,Q,Lda,Tau,Work,Lwork,info)
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
            CALL CLARNV(2,iseed,mc,C(1,j))
         ENDDO
         cnorm = CLANGE('1',mc,nc,C,Lda,Rwork)
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
            CALL CLACPY('Full',mc,nc,C,Lda,Cc,Lda)
!
!           Apply Q or Q' to C
!
            SRNamt = 'CUNMQR'
            CALL CUNMQR(side,trans,mc,nc,K,Af,Lda,Tau,Cc,Lda,Work,Lwork,&
     &                  info)
!
!           Form explicit product and subtract
!
            IF ( LSAME(side,'L') ) THEN
               CALL CGEMM(trans,'No transpose',mc,nc,mc,CMPLX(-ONE),Q,  &
     &                    Lda,C,Lda,CMPLX(ONE),Cc,Lda)
            ELSE
               CALL CGEMM('No transpose',trans,mc,nc,nc,CMPLX(-ONE),C,  &
     &                    Lda,Q,Lda,CMPLX(ONE),Cc,Lda)
            ENDIF
!
!           Compute error in the difference
!
            resid = CLANGE('1',mc,nc,Cc,Lda,Rwork)
            Result((iside-1)*2+itrans)                                  &
     &         = resid/(REAL(MAX(1,M))*cnorm*eps)
!
         ENDDO
      ENDDO
!
!
!     End of CQRT03
!
      END SUBROUTINE CQRT03
