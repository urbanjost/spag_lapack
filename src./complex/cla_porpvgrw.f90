!*==cla_porpvgrw.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLA_PORPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric or Hermitian positive-definite matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLA_PORPVGRW + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_porpvgrw.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_porpvgrw.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_porpvgrw.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL FUNCTION CLA_PORPVGRW( UPLO, NCOLS, A, LDA, AF, LDAF, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER*1        UPLO
!       INTEGER            NCOLS, LDA, LDAF
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), AF( LDAF, * )
!       REAL               WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>
!> CLA_PORPVGRW computes the reciprocal pivot growth factor
!> norm(A)/norm(U). The "max absolute element" norm is used. If this is
!> much less than 1, the stability of the LU factorization of the
!> (equilibrated) matrix A could be poor. This also means that the
!> solution X, estimated condition numbers, and error bounds could be
!> unreliable.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>       = 'U':  Upper triangle of A is stored;
!>       = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] NCOLS
!> \verbatim
!>          NCOLS is INTEGER
!>     The number of columns of the matrix A. NCOLS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>     On entry, the N-by-N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>     The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is COMPLEX array, dimension (LDAF,N)
!>     The triangular factor U or L from the Cholesky factorization
!>     A = U**T*U or A = L*L**T, as computed by CPOTRF.
!> \endverbatim
!>
!> \param[in] LDAF
!> \verbatim
!>          LDAF is INTEGER
!>     The leading dimension of the array AF.  LDAF >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (2*N)
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
!> \date June 2016
!
!> \ingroup complexPOcomputational
!
!  =====================================================================
      FUNCTION CLA_PORPVGRW(Uplo,Ncols,A,Lda,Af,Ldaf,Work)
      USE S_LSAME
      IMPLICIT NONE
!*--CLA_PORPVGRW110
      REAL :: CLA_PORPVGRW
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER(1) :: Uplo
      INTEGER , INTENT(IN) :: Ncols
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      REAL :: amax , rpvgrw , umax
      REAL :: CABS1
      INTEGER :: i , j
      LOGICAL :: upper
      COMPLEX :: zdum
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Statement Functions ..
!     ..
!     .. Statement Function Definitions ..
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
!     ..
!     .. Executable Statements ..
      upper = LSAME('Upper',Uplo)
!
!     SPOTRF will have factored only the NCOLSxNCOLS leading minor, so
!     we restrict the growth search to that minor and use only the first
!     2*NCOLS workspace entries.
!
      rpvgrw = 1.0
      DO i = 1 , 2*Ncols
         Work(i) = 0.0
      ENDDO
!
!     Find the max magnitude entry of each column.
!
      IF ( upper ) THEN
         DO j = 1 , Ncols
            DO i = 1 , j
               Work(Ncols+j) = MAX(CABS1(A(i,j)),Work(Ncols+j))
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , Ncols
            DO i = j , Ncols
               Work(Ncols+j) = MAX(CABS1(A(i,j)),Work(Ncols+j))
            ENDDO
         ENDDO
      ENDIF
!
!     Now find the max magnitude entry of each column of the factor in
!     AF.  No pivoting, so no permutations.
!
      IF ( LSAME('Upper',Uplo) ) THEN
         DO j = 1 , Ncols
            DO i = 1 , j
               Work(j) = MAX(CABS1(Af(i,j)),Work(j))
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , Ncols
            DO i = j , Ncols
               Work(j) = MAX(CABS1(Af(i,j)),Work(j))
            ENDDO
         ENDDO
      ENDIF
!
!     Compute the *inverse* of the max element growth factor.  Dividing
!     by zero would imply the largest entry of the factor's column is
!     zero.  Than can happen when either the column of A is zero or
!     massive pivots made the factor underflow to zero.  Neither counts
!     as growth in itself, so simply ignore terms with zero
!     denominators.
!
      IF ( LSAME('Upper',Uplo) ) THEN
         DO i = 1 , Ncols
            umax = Work(i)
            amax = Work(Ncols+i)
            IF ( umax/=0.0 ) rpvgrw = MIN(amax/umax,rpvgrw)
         ENDDO
      ELSE
         DO i = 1 , Ncols
            umax = Work(i)
            amax = Work(Ncols+i)
            IF ( umax/=0.0 ) rpvgrw = MIN(amax/umax,rpvgrw)
         ENDDO
      ENDIF
 
      CLA_PORPVGRW = rpvgrw
      END FUNCTION CLA_PORPVGRW
