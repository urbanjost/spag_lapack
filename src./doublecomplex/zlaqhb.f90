!*==zlaqhb.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLAQHB scales a Hermitian band matrix, using scaling factors computed by cpbequ.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAQHB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqhb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqhb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqhb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAQHB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED )
!
!       .. Scalar Arguments ..
!       CHARACTER          EQUED, UPLO
!       INTEGER            KD, LDAB, N
!       DOUBLE PRECISION   AMAX, SCOND
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   S( * )
!       COMPLEX*16         AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAQHB equilibrates a Hermitian band matrix A
!> using the scaling factors in the vector S.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored.
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
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          On entry, the upper or lower triangle of the symmetric band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>
!>          On exit, if INFO = 0, the triangular factor U or L from the
!>          Cholesky factorization A = U**H *U or A = L*L**H of the band
!>          matrix A, in the same storage format as A.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (N)
!>          The scale factors for A.
!> \endverbatim
!>
!> \param[in] SCOND
!> \verbatim
!>          SCOND is DOUBLE PRECISION
!>          Ratio of the smallest S(i) to the largest S(i).
!> \endverbatim
!>
!> \param[in] AMAX
!> \verbatim
!>          AMAX is DOUBLE PRECISION
!>          Absolute value of largest matrix entry.
!> \endverbatim
!>
!> \param[out] EQUED
!> \verbatim
!>          EQUED is CHARACTER*1
!>          Specifies whether or not equilibration was done.
!>          = 'N':  No equilibration.
!>          = 'Y':  Equilibration was done, i.e., A has been replaced by
!>                  diag(S) * A * diag(S).
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  THRESH is a threshold value used to decide if scaling should be done
!>  based on the ratio of the scaling factors.  If SCOND < THRESH,
!>  scaling is done.
!>
!>  LARGE and SMALL are threshold values used to decide if scaling should
!>  be done based on the absolute size of the largest matrix element.
!>  If AMAX > LARGE or AMAX < SMALL, scaling is done.
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLAQHB(Uplo,N,Kd,Ab,Ldab,S,Scond,Amax,Equed)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_LSAME
      IMPLICIT NONE
!*--ZLAQHB148
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Scond
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: cj , large , small
      INTEGER :: i , j
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
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N<=0 ) THEN
         Equed = 'N'
         RETURN
      ENDIF
!
!     Initialize LARGE and SMALL.
!
      small = DLAMCH('Safe minimum')/DLAMCH('Precision')
      large = ONE/small
!
      IF ( Scond>=THRESH .AND. Amax>=small .AND. Amax<=large ) THEN
!
!        No equilibration
!
         Equed = 'N'
      ELSE
!
!        Replace A by diag(S) * A * diag(S).
!
         IF ( LSAME(Uplo,'U') ) THEN
!
!           Upper triangle of A is stored in band format.
!
            DO j = 1 , N
               cj = S(j)
               DO i = MAX(1,j-Kd) , j - 1
                  Ab(Kd+1+i-j,j) = cj*S(i)*Ab(Kd+1+i-j,j)
               ENDDO
               Ab(Kd+1,j) = cj*cj*DBLE(Ab(Kd+1,j))
            ENDDO
         ELSE
!
!           Lower triangle of A is stored.
!
            DO j = 1 , N
               cj = S(j)
               Ab(1,j) = cj*cj*DBLE(Ab(1,j))
               DO i = j + 1 , MIN(N,j+Kd)
                  Ab(1+i-j,j) = cj*S(i)*Ab(1+i-j,j)
               ENDDO
            ENDDO
         ENDIF
         Equed = 'Y'
      ENDIF
!
!
!     End of ZLAQHB
!
      END SUBROUTINE ZLAQHB
