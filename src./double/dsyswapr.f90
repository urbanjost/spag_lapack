!*==dsyswapr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DSYSWAPR applies an elementary permutation on the rows and columns of a symmetric matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSYSWAPR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyswapr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyswapr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyswapr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYSWAPR( UPLO, N, A, LDA, I1, I2)
!
!       .. Scalar Arguments ..
!       CHARACTER        UPLO
!       INTEGER          I1, I2, LDA, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A( LDA, N )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYSWAPR applies an elementary permutation on the rows and the columns of
!> a symmetric matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are stored
!>          as an upper or lower triangular matrix.
!>          = 'U':  Upper triangular, form is A = U*D*U**T;
!>          = 'L':  Lower triangular, form is A = L*D*L**T.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the NB diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by DSYTRF.
!>
!>          On exit, if INFO = 0, the (symmetric) inverse of the original
!>          matrix.  If UPLO = 'U', the upper triangular part of the
!>          inverse is formed and the part of A below the diagonal is not
!>          referenced; if UPLO = 'L' the lower triangular part of the
!>          inverse is formed and the part of A above the diagonal is
!>          not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] I1
!> \verbatim
!>          I1 is INTEGER
!>          Index of the first row to swap
!> \endverbatim
!>
!> \param[in] I2
!> \verbatim
!>          I2 is INTEGER
!>          Index of the second row to swap
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
!> \ingroup doubleSYauxiliary
!
!  =====================================================================
      SUBROUTINE DSYSWAPR(Uplo,N,A,Lda,I1,I2)
      USE F77KINDS                        
      USE S_DSWAP
      USE S_LSAME
      IMPLICIT NONE
!*--DSYSWAPR109
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,N) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) :: I1
      INTEGER , INTENT(IN) :: I2
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i
      REAL(R8KIND) :: tmp
      LOGICAL :: upper
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
!
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
      upper = LSAME(Uplo,'U')
      IF ( upper ) THEN
!
!         UPPER
!         first swap
!          - swap column I1 and I2 from I1 to I1-1
         CALL DSWAP(I1-1,A(1,I1),1,A(1,I2),1)
!
!          second swap :
!          - swap A(I1,I1) and A(I2,I2)
!          - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1
         tmp = A(I1,I1)
         A(I1,I1) = A(I2,I2)
         A(I2,I2) = tmp
!
         DO i = 1 , I2 - I1 - 1
            tmp = A(I1,I1+i)
            A(I1,I1+i) = A(I1+i,I2)
            A(I1+i,I2) = tmp
         ENDDO
!
!          third swap
!          - swap row I1 and I2 from I2+1 to N
         DO i = I2 + 1 , N
            tmp = A(I1,i)
            A(I1,i) = A(I2,i)
            A(I2,i) = tmp
         ENDDO
!
      ELSE
!
!         LOWER
!         first swap
!          - swap row I1 and I2 from I1 to I1-1
         CALL DSWAP(I1-1,A(I1,1),Lda,A(I2,1),Lda)
!
!         second swap :
!          - swap A(I1,I1) and A(I2,I2)
!          - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1
         tmp = A(I1,I1)
         A(I1,I1) = A(I2,I2)
         A(I2,I2) = tmp
!
         DO i = 1 , I2 - I1 - 1
            tmp = A(I1+i,I1)
            A(I1+i,I1) = A(I2,I1+i)
            A(I2,I1+i) = tmp
         ENDDO
!
!         third swap
!          - swap col I1 and I2 from I2+1 to N
         DO i = I2 + 1 , N
            tmp = A(i,I1)
            A(i,I1) = A(i,I2)
            A(i,I2) = tmp
         ENDDO
!
      ENDIF
      END SUBROUTINE DSYSWAPR
