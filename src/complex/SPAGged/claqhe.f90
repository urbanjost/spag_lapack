!*==claqhe.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLAQHE scales a Hermitian matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAQHE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqhe.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqhe.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqhe.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAQHE( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )
!
!       .. Scalar Arguments ..
!       CHARACTER          EQUED, UPLO
!       INTEGER            LDA, N
!       REAL               AMAX, SCOND
!       ..
!       .. Array Arguments ..
!       REAL               S( * )
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAQHE equilibrates a Hermitian matrix A using the scaling factors
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
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          n by n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n by n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if EQUED = 'Y', the equilibrated matrix:
!>          diag(S) * A * diag(S).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(N,1).
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
!> \ingroup complexHEauxiliary
!
!  =====================================================================
      SUBROUTINE CLAQHE(Uplo,N,A,Lda,S,Scond,Amax,Equed)
      USE S_LSAME
      USE S_SLAMCH
      IMPLICIT NONE
!*--CLAQHE140
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Scond
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
!
! Local variable declarations rewritten by SPAG
!
      REAL :: cj , large , small
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
            DO j = 1 , N
               cj = S(j)
               DO i = 1 , j - 1
                  A(i,j) = cj*S(i)*A(i,j)
               ENDDO
               A(j,j) = cj*cj*REAL(A(j,j))
            ENDDO
         ELSE
!
!           Lower triangle of A is stored.
!
            DO j = 1 , N
               cj = S(j)
               A(j,j) = cj*cj*REAL(A(j,j))
               DO i = j + 1 , N
                  A(i,j) = cj*S(i)*A(i,j)
               ENDDO
            ENDDO
         ENDIF
         Equed = 'Y'
      ENDIF
!
!
!     End of CLAQHE
!
      END SUBROUTINE CLAQHE
