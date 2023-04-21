!*==dla_porpvgrw.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLA_PORPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric or Hermitian positive-definite matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLA_PORPVGRW + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_porpvgrw.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_porpvgrw.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_porpvgrw.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLA_PORPVGRW( UPLO, NCOLS, A, LDA, AF,
!                                               LDAF, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER*1        UPLO
!       INTEGER            NCOLS, LDA, LDAF
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>
!> DLA_PORPVGRW computes the reciprocal pivot growth factor
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
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
!>          AF is DOUBLE PRECISION array, dimension (LDAF,N)
!>     The triangular factor U or L from the Cholesky factorization
!>     A = U**T*U or A = L*L**T, as computed by DPOTRF.
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
!>          WORK is DOUBLE PRECISION array, dimension (2*N)
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
!> \ingroup doublePOcomputational
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION DLA_PORPVGRW(Uplo,Ncols,A,Lda,Af,Ldaf,  &
     &   Work)
      IMPLICIT NONE
!*--DLA_PORPVGRW110
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*1 Uplo
      INTEGER Ncols , Lda , Ldaf
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , Af(Ldaf,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , j
      DOUBLE PRECISION amax , umax , rpvgrw
      LOGICAL upper
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     ..
!     .. External Functions ..
      EXTERNAL LSAME
      LOGICAL LSAME
!     ..
!     .. Executable Statements ..
!
      upper = LSAME('Upper',Uplo)
!
!     DPOTRF will have factored only the NCOLSxNCOLS leading minor, so
!     we restrict the growth search to that minor and use only the first
!     2*NCOLS workspace entries.
!
      rpvgrw = 1.0D+0
      DO i = 1 , 2*Ncols
         Work(i) = 0.0D+0
      ENDDO
!
!     Find the max magnitude entry of each column.
!
      IF ( upper ) THEN
         DO j = 1 , Ncols
            DO i = 1 , j
               Work(Ncols+j) = MAX(ABS(A(i,j)),Work(Ncols+j))
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , Ncols
            DO i = j , Ncols
               Work(Ncols+j) = MAX(ABS(A(i,j)),Work(Ncols+j))
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
               Work(j) = MAX(ABS(Af(i,j)),Work(j))
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , Ncols
            DO i = j , Ncols
               Work(j) = MAX(ABS(Af(i,j)),Work(j))
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
            IF ( umax/=0.0D+0 ) rpvgrw = MIN(amax/umax,rpvgrw)
         ENDDO
      ELSE
         DO i = 1 , Ncols
            umax = Work(i)
            amax = Work(Ncols+i)
            IF ( umax/=0.0D+0 ) rpvgrw = MIN(amax/umax,rpvgrw)
         ENDDO
      ENDIF
 
      DLA_PORPVGRW = rpvgrw
      END FUNCTION DLA_PORPVGRW
