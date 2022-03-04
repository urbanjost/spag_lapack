!*==claqhp.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLAQHP scales a Hermitian matrix stored in packed form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAQHP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqhp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqhp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqhp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAQHP( UPLO, N, AP, S, SCOND, AMAX, EQUED )
!
!       .. Scalar Arguments ..
!       CHARACTER          EQUED, UPLO
!       INTEGER            N
!       REAL               AMAX, SCOND
!       ..
!       .. Array Arguments ..
!       REAL               S( * )
!       COMPLEX            AP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAQHP equilibrates a Hermitian matrix A using the scaling factors
!> in the vector S.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix A is stored.
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
!> \param[in,out] AP
!> \verbatim
!>          AP is COMPLEX array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the Hermitian matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!>
!>          On exit, the equilibrated matrix:  diag(S) * A * diag(S), in
!>          the same storage format as A.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension (N)
!>          The scale factors for A.
!> \endverbatim
!>
!> \param[in] SCOND
!> \verbatim
!>          SCOND is REAL
!>          Ratio of the smallest S(i) to the largest S(i).
!> \endverbatim
!>
!> \param[in] AMAX
!> \verbatim
!>          AMAX is REAL
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLAQHP(Uplo,N,Ap,S,Scond,Amax,Equed)
      USE S_LSAME
      USE S_SLAMCH
      IMPLICIT NONE
!*--CLAQHP132
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL , INTENT(IN) , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Scond
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
!
! Local variable declarations rewritten by SPAG
!
      REAL :: cj , large , small
      INTEGER :: i , j , jc
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
      small = SLAMCH('Safe minimum')/SLAMCH('Precision')
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
!           Upper triangle of A is stored.
!
            jc = 1
            DO j = 1 , N
               cj = S(j)
               DO i = 1 , j - 1
                  Ap(jc+i-1) = cj*S(i)*Ap(jc+i-1)
               ENDDO
               Ap(jc+j-1) = cj*cj*REAL(Ap(jc+j-1))
               jc = jc + j
            ENDDO
         ELSE
!
!           Lower triangle of A is stored.
!
            jc = 1
            DO j = 1 , N
               cj = S(j)
               Ap(jc) = cj*cj*REAL(Ap(jc))
               DO i = j + 1 , N
                  Ap(jc+i-j) = cj*S(i)*Ap(jc+i-j)
               ENDDO
               jc = jc + N - j + 1
            ENDDO
         ENDIF
         Equed = 'Y'
      ENDIF
!
!
!     End of CLAQHP
!
      END SUBROUTINE CLAQHP
