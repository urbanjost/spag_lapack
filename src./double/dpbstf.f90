!*==dpbstf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DPBSTF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DPBSTF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpbstf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpbstf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpbstf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DPBSTF( UPLO, N, KD, AB, LDAB, INFO )
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
!> DPBSTF computes a split Cholesky factorization of a real
!> symmetric positive definite band matrix A.
!>
!> This routine is designed to be used in conjunction with DSBGST.
!>
!> The factorization has the form  A = S**T*S  where S is a band matrix
!> of the same bandwidth as A and the following structure:
!>
!>   S = ( U    )
!>       ( M  L )
!>
!> where U is upper triangular of order m = (n+kd)/2, and L is lower
!> triangular of order n-m.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
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
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          On entry, the upper or lower triangle of the symmetric band
!>          matrix A, stored in the first kd+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>
!>          On exit, if INFO = 0, the factor S from the split Cholesky
!>          factorization A = S**T*S. See Further Details.
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
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, the factorization could not be completed,
!>               because the updated element a(i,i) was negative; the
!>               matrix A is not positive definite.
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
!>  N = 7, KD = 2:
!>
!>  S = ( s11  s12  s13                     )
!>      (      s22  s23  s24                )
!>      (           s33  s34                )
!>      (                s44                )
!>      (           s53  s54  s55           )
!>      (                s64  s65  s66      )
!>      (                     s75  s76  s77 )
!>
!>  If UPLO = 'U', the array AB holds:
!>
!>  on entry:                          on exit:
!>
!>   *    *   a13  a24  a35  a46  a57   *    *   s13  s24  s53  s64  s75
!>   *   a12  a23  a34  a45  a56  a67   *   s12  s23  s34  s54  s65  s76
!>  a11  a22  a33  a44  a55  a66  a77  s11  s22  s33  s44  s55  s66  s77
!>
!>  If UPLO = 'L', the array AB holds:
!>
!>  on entry:                          on exit:
!>
!>  a11  a22  a33  a44  a55  a66  a77  s11  s22  s33  s44  s55  s66  s77
!>  a21  a32  a43  a54  a65  a76   *   s12  s23  s34  s54  s65  s76   *
!>  a31  a42  a53  a64  a64   *    *   s13  s24  s53  s64  s75   *    *
!>
!>  Array elements marked * are not used by the routine.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DPBSTF(Uplo,N,Kd,Ab,Ldab,Info)
      USE F77KINDS                        
      USE S_DSCAL
      USE S_DSYR
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DPBSTF161
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
      INTEGER :: j , kld , km , m
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
         CALL XERBLA('DPBSTF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      kld = MAX(1,Ldab-1)
!
!     Set the splitting point m.
!
      m = (N+Kd)/2
!
      IF ( upper ) THEN
!
!        Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m).
!
         DO j = N , m + 1 , -1
!
!           Compute s(j,j) and test for non-positive-definiteness.
!
            ajj = Ab(Kd+1,j)
            IF ( ajj<=ZERO ) GOTO 100
            ajj = SQRT(ajj)
            Ab(Kd+1,j) = ajj
            km = MIN(j-1,Kd)
!
!           Compute elements j-km:j-1 of the j-th column and update the
!           the leading submatrix within the band.
!
            CALL DSCAL(km,ONE/ajj,Ab(Kd+1-km,j),1)
            CALL DSYR('Upper',km,-ONE,Ab(Kd+1-km,j),1,Ab(Kd+1,j-km),kld)
         ENDDO
!
!        Factorize the updated submatrix A(1:m,1:m) as U**T*U.
!
         DO j = 1 , m
!
!           Compute s(j,j) and test for non-positive-definiteness.
!
            ajj = Ab(Kd+1,j)
            IF ( ajj<=ZERO ) GOTO 100
            ajj = SQRT(ajj)
            Ab(Kd+1,j) = ajj
            km = MIN(Kd,m-j)
!
!           Compute elements j+1:j+km of the j-th row and update the
!           trailing submatrix within the band.
!
            IF ( km>0 ) THEN
               CALL DSCAL(km,ONE/ajj,Ab(Kd,j+1),kld)
               CALL DSYR('Upper',km,-ONE,Ab(Kd,j+1),kld,Ab(Kd+1,j+1),   &
     &                   kld)
            ENDIF
         ENDDO
      ELSE
!
!        Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m).
!
         DO j = N , m + 1 , -1
!
!           Compute s(j,j) and test for non-positive-definiteness.
!
            ajj = Ab(1,j)
            IF ( ajj<=ZERO ) GOTO 100
            ajj = SQRT(ajj)
            Ab(1,j) = ajj
            km = MIN(j-1,Kd)
!
!           Compute elements j-km:j-1 of the j-th row and update the
!           trailing submatrix within the band.
!
            CALL DSCAL(km,ONE/ajj,Ab(km+1,j-km),kld)
            CALL DSYR('Lower',km,-ONE,Ab(km+1,j-km),kld,Ab(1,j-km),kld)
         ENDDO
!
!        Factorize the updated submatrix A(1:m,1:m) as U**T*U.
!
         DO j = 1 , m
!
!           Compute s(j,j) and test for non-positive-definiteness.
!
            ajj = Ab(1,j)
            IF ( ajj<=ZERO ) GOTO 100
            ajj = SQRT(ajj)
            Ab(1,j) = ajj
            km = MIN(Kd,m-j)
!
!           Compute elements j+1:j+km of the j-th column and update the
!           trailing submatrix within the band.
!
            IF ( km>0 ) THEN
               CALL DSCAL(km,ONE/ajj,Ab(2,j),1)
               CALL DSYR('Lower',km,-ONE,Ab(2,j),1,Ab(1,j+1),kld)
            ENDIF
         ENDDO
      ENDIF
      RETURN
!
 100  Info = j
!
!     End of DPBSTF
!
      END SUBROUTINE DPBSTF