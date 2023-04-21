!*==zla_syrpvgrw.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLA_SYRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric indefinite matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLA_SYRPVGRW + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_syrpvgrw.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_syrpvgrw.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_syrpvgrw.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF,
!                                               LDAF, IPIV, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER*1        UPLO
!       INTEGER            N, INFO, LDA, LDAF
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), AF( LDAF, * )
!       DOUBLE PRECISION   WORK( * )
!       INTEGER            IPIV( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>
!> ZLA_SYRPVGRW computes the reciprocal pivot growth factor
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
!>     The value of INFO returned from ZSYTRF, .i.e., the pivot in
!>     column INFO is exactly 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
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
!>          AF is COMPLEX*16 array, dimension (LDAF,N)
!>     The block diagonal matrix D and the multipliers used to
!>     obtain the factor U or L as computed by ZSYTRF.
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
!>     as determined by ZSYTRF.
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
!> \ingroup complex16SYcomputational
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION ZLA_SYRPVGRW(Uplo,N,Info,A,Lda,Af,Ldaf, &
     &   Ipiv,Work)
      IMPLICIT NONE
!*--ZLA_SYRPVGRW127
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*1 Uplo
      INTEGER N , Info , Lda , Ldaf
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , Af(Ldaf,*)
      DOUBLE PRECISION Work(*)
      INTEGER Ipiv(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER ncols , i , j , k , kp
      DOUBLE PRECISION amax , umax , rpvgrw , tmp
      LOGICAL upper
      COMPLEX*16 zdum
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , REAL , DIMAG , MAX , MIN
!     ..
!     .. External Subroutines ..
      EXTERNAL LSAME
      LOGICAL LSAME
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION CABS1
!     ..
!     .. Statement Function Definitions ..
      CABS1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))
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
 
      rpvgrw = 1.0D+0
      DO i = 1 , 2*N
         Work(i) = 0.0D+0
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
!     The iteration orders and permutations were copied from zsytrs.
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
            IF ( umax/=0.0D+0 ) rpvgrw = MIN(amax/umax,rpvgrw)
         ENDDO
      ELSE
         DO i = 1 , ncols
            umax = Work(i)
            amax = Work(N+i)
            IF ( umax/=0.0D+0 ) rpvgrw = MIN(amax/umax,rpvgrw)
         ENDDO
      ENDIF
 
      ZLA_SYRPVGRW = rpvgrw
      END FUNCTION ZLA_SYRPVGRW
