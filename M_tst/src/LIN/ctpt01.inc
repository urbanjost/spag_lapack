!*==ctpt01.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b ctpt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ctpt01( UPLO, DIAG, N, AP, AINVP, RCOND, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            N
!       REAL               RCOND, RESID
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            AINVP( * ), AP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ctpt01 computes the residual for a triangular matrix A times its
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
!>          AP is COMPLEX array, dimension (N*(N+1)/2)
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
!>          AINVP is COMPLEX array, dimension (N*(N+1)/2)
!>          On entry, the (triangular) inverse of the matrix A, packed
!>          columnwise in a linear array as in AP.
!>          On exit, the contents of AINVP are destroyed.
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The reciprocal condition number of A, computed as
!>          1/(norm(A) * norm(AINV)).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE ctpt01(Uplo,Diag,N,Ap,Ainvp,Rcond,Rwork,Resid)
      IMPLICIT NONE
!*--CTPT01113
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Uplo
      INTEGER N
      REAL Rcond , Resid
!     ..
!     .. Array Arguments ..
      REAL Rwork(*)
      COMPLEX Ainvp(*) , Ap(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL unitd
      INTEGER j , jc
      REAL ainvnm , anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANTP , SLAMCH
EXTERNAL    LSAME,   CLANTP,   SLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    CTPMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC REAL
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
      eps = SLAMCH('Epsilon')
      anorm = CLANTP('1',Uplo,Diag,N,Ap,Rwork)
      ainvnm = CLANTP('1',Uplo,Diag,N,Ainvp,Rwork)
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
            CALL CTPMV('Upper','No transpose',Diag,j,Ap,Ainvp(jc),1)
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
            CALL CTPMV('Lower','No transpose',Diag,N-j+1,Ap(jc),        &
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
      Resid = CLANTP('1',Uplo,'Non-unit',N,Ainvp,Rwork)
!
      Resid = ((Resid*Rcond)/REAL(N))/eps
!
!
!     End of ctpt01
!
      END SUBROUTINE ctpt01
