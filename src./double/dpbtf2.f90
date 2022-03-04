!*==dpbtf2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DPBTF2 computes the Cholesky factorization of a symmetric/Hermitian positive definite band matrix (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DPBTF2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpbtf2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpbtf2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpbtf2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KD, LDAB, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DPBTF2 computes the Cholesky factorization of a real symmetric
!> positive definite band matrix A.
!>
!> The factorization has the form
!>    A = U**T * U ,  if UPLO = 'U', or
!>    A = L  * L**T,  if UPLO = 'L',
!> where U is an upper triangular matrix, U**T is the transpose of U, and
!> L is lower triangular.
!>
!> This is the unblocked version of the algorithm, calling Level 2 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of super-diagonals of the matrix A if UPLO = 'U',
!>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          On entry, the upper or lower triangle of the symmetric band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>
!>          On exit, if INFO = 0, the triangular factor U or L from the
!>          Cholesky factorization A = U**T*U or A = L*L**T of the band
!>          matrix A, in the same storage format as A.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
!>          > 0: if INFO = k, the leading minor of order k is not
!>               positive definite, and the factorization could not be
!>               completed.
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
!> \ingroup doubleOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The band storage scheme is illustrated by the following example, when
!>  N = 6, KD = 2, and UPLO = 'U':
!>
!>  On entry:                       On exit:
!>
!>      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!>
!>  Similarly, if UPLO = 'L' the format of A is as follows:
!>
!>  On entry:                       On exit:
!>
!>     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!>     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!>     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!>
!>  Array elements marked * are not used by the routine.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DPBTF2(Uplo,N,Kd,Ab,Ldab,Info)
      USE F77KINDS                        
      USE S_DSCAL
      USE S_DSYR
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DPBTF2151
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ajj
      INTEGER :: j , kld , kn
      LOGICAL :: upper
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Kd<0 ) THEN
         Info = -3
      ELSEIF ( Ldab<Kd+1 ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DPBTF2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      kld = MAX(1,Ldab-1)
!
      IF ( upper ) THEN
!
!        Compute the Cholesky factorization A = U**T*U.
!
         DO j = 1 , N
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            ajj = Ab(Kd+1,j)
            IF ( ajj<=ZERO ) GOTO 100
            ajj = SQRT(ajj)
            Ab(Kd+1,j) = ajj
!
!           Compute elements J+1:J+KN of row J and update the
!           trailing submatrix within the band.
!
            kn = MIN(Kd,N-j)
            IF ( kn>0 ) THEN
               CALL DSCAL(kn,ONE/ajj,Ab(Kd,j+1),kld)
               CALL DSYR('Upper',kn,-ONE,Ab(Kd,j+1),kld,Ab(Kd+1,j+1),   &
     &                   kld)
            ENDIF
         ENDDO
      ELSE
!
!        Compute the Cholesky factorization A = L*L**T.
!
         DO j = 1 , N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            ajj = Ab(1,j)
            IF ( ajj<=ZERO ) GOTO 100
            ajj = SQRT(ajj)
            Ab(1,j) = ajj
!
!           Compute elements J+1:J+KN of column J and update the
!           trailing submatrix within the band.
!
            kn = MIN(Kd,N-j)
            IF ( kn>0 ) THEN
               CALL DSCAL(kn,ONE/ajj,Ab(2,j),1)
               CALL DSYR('Lower',kn,-ONE,Ab(2,j),1,Ab(1,j+1),kld)
            ENDIF
         ENDDO
      ENDIF
      RETURN
!
 100  Info = j
!
!     End of DPBTF2
!
      END SUBROUTINE DPBTF2
