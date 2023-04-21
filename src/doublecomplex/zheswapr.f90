!*==zheswapr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZHESWAPR applies an elementary permutation on the rows and columns of a Hermitian matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHESWAPR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheswapr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheswapr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheswapr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHESWAPR( UPLO, N, A, LDA, I1, I2)
!
!       .. Scalar Arguments ..
!       CHARACTER        UPLO
!       INTEGER          I1, I2, LDA, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16          A( LDA, N )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHESWAPR applies an elementary permutation on the rows and the columns of
!> a hermitian matrix.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the NB diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by CSYTRF.
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
!> \ingroup complex16HEauxiliary
!
!  =====================================================================
      SUBROUTINE ZHESWAPR(Uplo,N,A,Lda,I1,I2)
      IMPLICIT NONE
!*--ZHESWAPR106
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER I1 , I2 , Lda , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,N)
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER i
      COMPLEX*16 tmp
!
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL ZSWAP
!     ..
!     .. Executable Statements ..
!
      upper = LSAME(Uplo,'U')
      IF ( upper ) THEN
!
!         UPPER
!         first swap
!          - swap column I1 and I2 from I1 to I1-1
         CALL ZSWAP(I1-1,A(1,I1),1,A(1,I2),1)
!
!          second swap :
!          - swap A(I1,I1) and A(I2,I2)
!          - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1
!          - swap A(I2,I1) and A(I1,I2)
 
         tmp = A(I1,I1)
         A(I1,I1) = A(I2,I2)
         A(I2,I2) = tmp
!
         DO i = 1 , I2 - I1 - 1
            tmp = A(I1,I1+i)
            A(I1,I1+i) = DCONJG(A(I1+i,I2))
            A(I1+i,I2) = DCONJG(tmp)
         ENDDO
!
         A(I1,I2) = DCONJG(A(I1,I2))
 
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
!          - swap row I1 and I2 from 1 to I1-1
         CALL ZSWAP(I1-1,A(I1,1),Lda,A(I2,1),Lda)
!
!         second swap :
!          - swap A(I1,I1) and A(I2,I2)
!          - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1
!          - swap A(I2,I1) and A(I1,I2)
 
         tmp = A(I1,I1)
         A(I1,I1) = A(I2,I2)
         A(I2,I2) = tmp
!
         DO i = 1 , I2 - I1 - 1
            tmp = A(I1+i,I1)
            A(I1+i,I1) = DCONJG(A(I2,I1+i))
            A(I2,I1+i) = DCONJG(tmp)
         ENDDO
!
         A(I2,I1) = DCONJG(A(I2,I1))
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
 
      END SUBROUTINE ZHESWAPR
