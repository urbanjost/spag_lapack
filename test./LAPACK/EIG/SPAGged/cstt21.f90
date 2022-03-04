!*==cstt21.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CSTT21
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSTT21( N, KBAND, AD, AE, SD, SE, U, LDU, WORK, RWORK,
!                          RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            KBAND, LDU, N
!       ..
!       .. Array Arguments ..
!       REAL               AD( * ), AE( * ), RESULT( 2 ), RWORK( * ),
!      $                   SD( * ), SE( * )
!       COMPLEX            U( LDU, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSTT21  checks a decomposition of the form
!>
!>    A = U S U**H
!>
!> where **H means conjugate transpose, A is real symmetric tridiagonal,
!> U is unitary, and S is real and diagonal (if KBAND=0) or symmetric
!> tridiagonal (if KBAND=1).  Two tests are performed:
!>
!>    RESULT(1) = | A - U S U**H | / ( |A| n ulp )
!>
!>    RESULT(2) = | I - U U**H | / ( n ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrix.  If it is zero, CSTT21 does nothing.
!>          It must be at least zero.
!> \endverbatim
!>
!> \param[in] KBAND
!> \verbatim
!>          KBAND is INTEGER
!>          The bandwidth of the matrix S.  It may only be zero or one.
!>          If zero, then S is diagonal, and SE is not referenced.  If
!>          one, then S is symmetric tri-diagonal.
!> \endverbatim
!>
!> \param[in] AD
!> \verbatim
!>          AD is REAL array, dimension (N)
!>          The diagonal of the original (unfactored) matrix A.  A is
!>          assumed to be real symmetric tridiagonal.
!> \endverbatim
!>
!> \param[in] AE
!> \verbatim
!>          AE is REAL array, dimension (N-1)
!>          The off-diagonal of the original (unfactored) matrix A.  A
!>          is assumed to be symmetric tridiagonal.  AE(1) is the (1,2)
!>          and (2,1) element, AE(2) is the (2,3) and (3,2) element, etc.
!> \endverbatim
!>
!> \param[in] SD
!> \verbatim
!>          SD is REAL array, dimension (N)
!>          The diagonal of the real (symmetric tri-) diagonal matrix S.
!> \endverbatim
!>
!> \param[in] SE
!> \verbatim
!>          SE is REAL array, dimension (N-1)
!>          The off-diagonal of the (symmetric tri-) diagonal matrix S.
!>          Not referenced if KBSND=0.  If KBAND=1, then AE(1) is the
!>          (1,2) and (2,1) element, SE(2) is the (2,3) and (3,2)
!>          element, etc.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is COMPLEX array, dimension (LDU, N)
!>          The unitary matrix in the decomposition.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  LDU must be at least N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N**2)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
!>          The values computed by the two tests described above.  The
!>          values are currently limited to 1/ulp, to avoid overflow.
!>          RESULT(1) is always modified.
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE CSTT21(N,Kband,Ad,Ae,Sd,Se,U,Ldu,Work,Rwork,Result)
      IMPLICIT NONE
!*--CSTT21136
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Kband , Ldu , N
!     ..
!     .. Array Arguments ..
      REAL Ad(*) , Ae(*) , Result(2) , Rwork(*) , Sd(*) , Se(*)
      COMPLEX U(Ldu,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER j
      REAL anorm , temp1 , temp2 , ulp , unfl , wnorm
!     ..
!     .. External Functions ..
      REAL CLANGE , CLANHE , SLAMCH
      EXTERNAL CLANGE , CLANHE , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMM , CHER , CHER2 , CLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CMPLX , MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
!     1)      Constants
!
      Result(1) = ZERO
      Result(2) = ZERO
      IF ( N<=0 ) RETURN
!
      unfl = SLAMCH('Safe minimum')
      ulp = SLAMCH('Precision')
!
!     Do Test 1
!
!     Copy A & Compute its 1-Norm:
!
      CALL CLASET('Full',N,N,CZERO,CZERO,Work,N)
!
      anorm = ZERO
      temp1 = ZERO
!
      DO j = 1 , N - 1
         Work((N+1)*(j-1)+1) = Ad(j)
         Work((N+1)*(j-1)+2) = Ae(j)
         temp2 = ABS(Ae(j))
         anorm = MAX(anorm,ABS(Ad(j))+temp1+temp2)
         temp1 = temp2
      ENDDO
!
      Work(N**2) = Ad(N)
      anorm = MAX(anorm,ABS(Ad(N))+temp1,unfl)
!
!     Norm of A - U S U**H
!
      DO j = 1 , N
         CALL CHER('L',N,-Sd(j),U(1,j),1,Work,N)
      ENDDO
!
      IF ( N>1 .AND. Kband==1 ) THEN
         DO j = 1 , N - 1
            CALL CHER2('L',N,-CMPLX(Se(j)),U(1,j),1,U(1,j+1),1,Work,N)
         ENDDO
      ENDIF
!
      wnorm = CLANHE('1','L',N,Work,N,Rwork)
!
      IF ( anorm>wnorm ) THEN
         Result(1) = (wnorm/anorm)/(N*ulp)
      ELSEIF ( anorm<ONE ) THEN
         Result(1) = (MIN(wnorm,N*anorm)/anorm)/(N*ulp)
      ELSE
         Result(1) = MIN(wnorm/anorm,REAL(N))/(N*ulp)
      ENDIF
!
!     Do Test 2
!
!     Compute  U U**H - I
!
      CALL CGEMM('N','C',N,N,N,CONE,U,Ldu,U,Ldu,CZERO,Work,N)
!
      DO j = 1 , N
         Work((N+1)*(j-1)+1) = Work((N+1)*(j-1)+1) - CONE
      ENDDO
!
      Result(2) = MIN(REAL(N),CLANGE('1',N,N,Work,N,Rwork))/(N*ulp)
!
!
!     End of CSTT21
!
      END SUBROUTINE CSTT21
