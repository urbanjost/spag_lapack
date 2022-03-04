!*==ztpttr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZTPTTR copies a triangular matrix from the standard packed format (TP) to the standard full format (TR).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTPTTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpttr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpttr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpttr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTPTTR( UPLO, N, AP, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, N, LDA
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), AP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTPTTR copies a triangular matrix A from standard packed format (TP)
!> to standard full format (TR).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  A is upper triangular.
!>          = 'L':  A is lower triangular.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A. N >= 0.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension ( N*(N+1)/2 ),
!>          On entry, the upper or lower triangular matrix A, packed
!>          columnwise in a linear array. The j-th column of A is stored
!>          in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, N )
!>          On exit, the triangular matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZTPTTR(Uplo,N,Ap,A,Lda,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--ZTPTTR111
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , j , k
      LOGICAL :: lower
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
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      lower = LSAME(Uplo,'L')
      IF ( .NOT.lower .AND. .NOT.LSAME(Uplo,'U') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZTPTTR',-Info)
         RETURN
      ENDIF
!
      IF ( lower ) THEN
         k = 0
         DO j = 1 , N
            DO i = j , N
               k = k + 1
               A(i,j) = Ap(k)
            ENDDO
         ENDDO
      ELSE
         k = 0
         DO j = 1 , N
            DO i = 1 , j
               k = k + 1
               A(i,j) = Ap(k)
            ENDDO
         ENDDO
      ENDIF
!
!
!
!     End of ZTPTTR
!
      END SUBROUTINE ZTPTTR
