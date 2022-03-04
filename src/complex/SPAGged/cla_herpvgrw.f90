!*==cla_herpvgrw.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLA_HERPVGRW
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLA_HERPVGRW + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_herpvgrw.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_herpvgrw.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_herpvgrw.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL FUNCTION CLA_HERPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV,
!                                   WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER*1        UPLO
!       INTEGER            N, INFO, LDA, LDAF
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
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
!> CLA_HERPVGRW computes the reciprocal pivot growth factor
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
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>     The number of linear equations, i.e., the order of the
!>     matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] INFO
!> \verbatim
!>          INFO is INTEGER
!>     The value of INFO returned from SSYTRF, .i.e., the pivot in
!>     column INFO is exactly 0.
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
!>     The block diagonal matrix D and the multipliers used to
!>     obtain the factor U or L as computed by CHETRF.
!> \endverbatim
!>
!> \param[in] LDAF
!> \verbatim
!>          LDAF is INTEGER
!>     The leading dimension of the array AF.  LDAF >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>     Details of the interchanges and the block structure of D
!>     as determined by CHETRF.
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
!> \ingroup complexHEcomputational
!
!  =====================================================================
      FUNCTION CLA_HERPVGRW(Uplo,N,Info,A,Lda,Af,Ldaf,Ipiv,Work)
      USE S_LSAME
      IMPLICIT NONE
!*--CLA_HERPVGRW127
      REAL :: CLA_HERPVGRW
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER(1) :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Info
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      REAL :: amax , rpvgrw , tmp , umax
      REAL :: CABS1
      INTEGER :: i , j , k , kp , ncols
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
!
      upper = LSAME('Upper',Uplo)
      IF ( Info/=0 ) THEN
         ncols = Info
      ELSEIF ( upper ) THEN
         ncols = 1
      ELSE
         ncols = N
      ENDIF
 
      rpvgrw = 1.0
      DO i = 1 , 2*N
         Work(i) = 0.0
      ENDDO
!
!     Find the max magnitude entry of each column of A.  Compute the max
!     for all N columns so we can apply the pivot permutation while
!     looping below.  Assume a full factorization is the common case.
!
      IF ( upper ) THEN
         DO j = 1 , N
            DO i = 1 , j
               Work(N+i) = MAX(CABS1(A(i,j)),Work(N+i))
               Work(N+j) = MAX(CABS1(A(i,j)),Work(N+j))
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , N
            DO i = j , N
               Work(N+i) = MAX(CABS1(A(i,j)),Work(N+i))
               Work(N+j) = MAX(CABS1(A(i,j)),Work(N+j))
            ENDDO
         ENDDO
      ENDIF
!
!     Now find the max magnitude entry of each column of U or L.  Also
!     permute the magnitudes of A above so they're in the same order as
!     the factor.
!
!     The iteration orders and permutations were copied from csytrs.
!     Calls to SSWAP would be severe overkill.
!
      IF ( upper ) THEN
         k = N
         DO WHILE ( k<ncols .AND. k>0 )
            IF ( Ipiv(k)>0 ) THEN
!              1x1 pivot
               kp = Ipiv(k)
               IF ( kp/=k ) THEN
                  tmp = Work(N+k)
                  Work(N+k) = Work(N+kp)
                  Work(N+kp) = tmp
               ENDIF
               DO i = 1 , k
                  Work(k) = MAX(CABS1(Af(i,k)),Work(k))
               ENDDO
               k = k - 1
            ELSE
!              2x2 pivot
               kp = -Ipiv(k)
               tmp = Work(N+k-1)
               Work(N+k-1) = Work(N+kp)
               Work(N+kp) = tmp
               DO i = 1 , k - 1
                  Work(k) = MAX(CABS1(Af(i,k)),Work(k))
                  Work(k-1) = MAX(CABS1(Af(i,k-1)),Work(k-1))
               ENDDO
               Work(k) = MAX(CABS1(Af(k,k)),Work(k))
               k = k - 2
            ENDIF
         ENDDO
         k = ncols
         DO WHILE ( k<=N )
            IF ( Ipiv(k)>0 ) THEN
               kp = Ipiv(k)
               IF ( kp/=k ) THEN
                  tmp = Work(N+k)
                  Work(N+k) = Work(N+kp)
                  Work(N+kp) = tmp
               ENDIF
               k = k + 1
            ELSE
               kp = -Ipiv(k)
               tmp = Work(N+k)
               Work(N+k) = Work(N+kp)
               Work(N+kp) = tmp
               k = k + 2
            ENDIF
         ENDDO
      ELSE
         k = 1
         DO WHILE ( k<=ncols )
            IF ( Ipiv(k)>0 ) THEN
!              1x1 pivot
               kp = Ipiv(k)
               IF ( kp/=k ) THEN
                  tmp = Work(N+k)
                  Work(N+k) = Work(N+kp)
                  Work(N+kp) = tmp
               ENDIF
               DO i = k , N
                  Work(k) = MAX(CABS1(Af(i,k)),Work(k))
               ENDDO
               k = k + 1
            ELSE
!              2x2 pivot
               kp = -Ipiv(k)
               tmp = Work(N+k+1)
               Work(N+k+1) = Work(N+kp)
               Work(N+kp) = tmp
               DO i = k + 1 , N
                  Work(k) = MAX(CABS1(Af(i,k)),Work(k))
                  Work(k+1) = MAX(CABS1(Af(i,k+1)),Work(k+1))
               ENDDO
               Work(k) = MAX(CABS1(Af(k,k)),Work(k))
               k = k + 2
            ENDIF
         ENDDO
         k = ncols
         DO WHILE ( k>=1 )
            IF ( Ipiv(k)>0 ) THEN
               kp = Ipiv(k)
               IF ( kp/=k ) THEN
                  tmp = Work(N+k)
                  Work(N+k) = Work(N+kp)
                  Work(N+kp) = tmp
               ENDIF
               k = k - 1
            ELSE
               kp = -Ipiv(k)
               tmp = Work(N+k)
               Work(N+k) = Work(N+kp)
               Work(N+kp) = tmp
               k = k - 2
            ENDIF
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
      IF ( upper ) THEN
         DO i = ncols , N
            umax = Work(i)
            amax = Work(N+i)
            IF ( umax/=0.0 ) rpvgrw = MIN(amax/umax,rpvgrw)
         ENDDO
      ELSE
         DO i = 1 , ncols
            umax = Work(i)
            amax = Work(N+i)
            IF ( umax/=0.0 ) rpvgrw = MIN(amax/umax,rpvgrw)
         ENDDO
      ENDIF
 
      CLA_HERPVGRW = rpvgrw
      END FUNCTION CLA_HERPVGRW
