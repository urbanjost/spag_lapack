!*==dtrtri.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DTRTRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DTRTRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrtri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrtri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrtri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRTRI computes the inverse of a real upper or lower triangular
!> matrix A.
!>
!> This is the Level 3 BLAS version of the algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  A is upper triangular;
!>          = 'L':  A is lower triangular.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          = 'N':  A is non-unit triangular;
!>          = 'U':  A is unit triangular.
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
!>          On entry, the triangular matrix A.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of the array A contains
!>          the upper triangular matrix, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of the array A contains
!>          the lower triangular matrix, and the strictly upper
!>          triangular part of A is not referenced.  If DIAG = 'U', the
!>          diagonal elements of A are also not referenced and are
!>          assumed to be 1.
!>          On exit, the (triangular) inverse of the original matrix, in
!>          the same storage format.
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
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
!>               matrix is singular and its inverse can not be computed.
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
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DTRTRI(Uplo,Diag,N,A,Lda,Info)
      IMPLICIT NONE
!*--DTRTRI113
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Uplo
      INTEGER Info , Lda , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL nounit , upper
      INTEGER j , jb , nb , nn
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      EXTERNAL LSAME , ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL DTRMM , DTRSM , DTRTI2 , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      nounit = LSAME(Diag,'N')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( .NOT.nounit .AND. .NOT.LSAME(Diag,'U') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DTRTRI',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Check for singularity if non-unit.
!
      IF ( nounit ) THEN
         DO Info = 1 , N
            IF ( A(Info,Info)==ZERO ) RETURN
         ENDDO
         Info = 0
      ENDIF
!
!     Determine the block size for this environment.
!
      nb = ILAENV(1,'DTRTRI',Uplo//Diag,N,-1,-1,-1)
      IF ( nb<=1 .OR. nb>=N ) THEN
!
!        Use unblocked code
!
         CALL DTRTI2(Uplo,Diag,N,A,Lda,Info)
!
!        Use blocked code
!
      ELSEIF ( upper ) THEN
!
!           Compute inverse of upper triangular matrix
!
         DO j = 1 , N , nb
            jb = MIN(nb,N-j+1)
!
!              Compute rows 1:j-1 of current block column
!
            CALL DTRMM('Left','Upper','No transpose',Diag,j-1,jb,ONE,A, &
     &                 Lda,A(1,j),Lda)
            CALL DTRSM('Right','Upper','No transpose',Diag,j-1,jb,-ONE, &
     &                 A(j,j),Lda,A(1,j),Lda)
!
!              Compute inverse of current diagonal block
!
            CALL DTRTI2('Upper',Diag,jb,A(j,j),Lda,Info)
         ENDDO
      ELSE
!
!           Compute inverse of lower triangular matrix
!
         nn = ((N-1)/nb)*nb + 1
         DO j = nn , 1 , -nb
            jb = MIN(nb,N-j+1)
            IF ( j+jb<=N ) THEN
!
!                 Compute rows j+jb:n of current block column
!
               CALL DTRMM('Left','Lower','No transpose',Diag,N-j-jb+1,  &
     &                    jb,ONE,A(j+jb,j+jb),Lda,A(j+jb,j),Lda)
               CALL DTRSM('Right','Lower','No transpose',Diag,N-j-jb+1, &
     &                    jb,-ONE,A(j,j),Lda,A(j+jb,j),Lda)
            ENDIF
!
!              Compute inverse of current diagonal block
!
            CALL DTRTI2('Lower',Diag,jb,A(j,j),Lda,Info)
         ENDDO
      ENDIF
!
!
!     End of DTRTRI
!
      END SUBROUTINE DTRTRI
