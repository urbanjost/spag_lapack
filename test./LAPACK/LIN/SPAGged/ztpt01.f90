!*==ztpt01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZTPT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTPT01( UPLO, DIAG, N, AP, AINVP, RCOND, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            N
!       DOUBLE PRECISION   RCOND, RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         AINVP( * ), AP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTPT01 computes the residual for a triangular matrix A times its
!> inverse when A is stored in packed format:
!>    RESID = norm(A*AINV - I) / ( N * norm(A) * norm(AINV) * EPS ),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The original upper or lower triangular matrix A, packed
!>          columnwise in a linear array.  The j-th column of A is stored
!>          in the array AP as follows:
!>          if UPLO = 'U', AP((j-1)*j/2 + i) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L',
!>             AP((j-1)*(n-j) + j*(j+1)/2 + i-j) = A(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[in] AINVP
!> \verbatim
!>          AINVP is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          On entry, the (triangular) inverse of the matrix A, packed
!>          columnwise in a linear array as in AP.
!>          On exit, the contents of AINVP are destroyed.
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The reciprocal condition number of A, computed as
!>          1/(norm(A) * norm(AINV)).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          norm(A*AINV - I) / ( N * norm(A) * norm(AINV) * EPS )
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
      SUBROUTINE ZTPT01(Uplo,Diag,N,Ap,Ainvp,Rcond,Rwork,Resid)
      IMPLICIT NONE
!*--ZTPT01113
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Uplo
      INTEGER N
      DOUBLE PRECISION Rcond , Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 Ainvp(*) , Ap(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL unitd
      INTEGER j , jc
      DOUBLE PRECISION ainvnm , anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , ZLANTP
      EXTERNAL LSAME , DLAMCH , ZLANTP
!     ..
!     .. External Subroutines ..
      EXTERNAL ZTPMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0.
!
      IF ( N<=0 ) THEN
         Rcond = ONE
         Resid = ZERO
         RETURN
      ENDIF
!
!     Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
!
      eps = DLAMCH('Epsilon')
      anorm = ZLANTP('1',Uplo,Diag,N,Ap,Rwork)
      ainvnm = ZLANTP('1',Uplo,Diag,N,Ainvp,Rwork)
      IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
         Rcond = ZERO
         Resid = ONE/eps
         RETURN
      ENDIF
      Rcond = (ONE/anorm)/ainvnm
!
!     Compute A * AINV, overwriting AINV.
!
      unitd = LSAME(Diag,'U')
      IF ( LSAME(Uplo,'U') ) THEN
         jc = 1
         DO j = 1 , N
            IF ( unitd ) Ainvp(jc+j-1) = ONE
!
!           Form the j-th column of A*AINV.
!
            CALL ZTPMV('Upper','No transpose',Diag,j,Ap,Ainvp(jc),1)
!
!           Subtract 1 from the diagonal to form A*AINV - I.
!
            Ainvp(jc+j-1) = Ainvp(jc+j-1) - ONE
            jc = jc + j
         ENDDO
      ELSE
         jc = 1
         DO j = 1 , N
            IF ( unitd ) Ainvp(jc) = ONE
!
!           Form the j-th column of A*AINV.
!
            CALL ZTPMV('Lower','No transpose',Diag,N-j+1,Ap(jc),        &
     &                 Ainvp(jc),1)
!
!           Subtract 1 from the diagonal to form A*AINV - I.
!
            Ainvp(jc) = Ainvp(jc) - ONE
            jc = jc + N - j + 1
         ENDDO
      ENDIF
!
!     Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)
!
      Resid = ZLANTP('1',Uplo,'Non-unit',N,Ainvp,Rwork)
!
      Resid = ((Resid*Rcond)/DBLE(N))/eps
!
!
!     End of ZTPT01
!
      END SUBROUTINE ZTPT01
