!*==slauu2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLAUU2 computes the product UUH or LHL, where U and L are upper or lower triangular matrices (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAUU2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slauu2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slauu2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slauu2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAUU2( UPLO, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAUU2 computes the product U * U**T or L**T * L, where the triangular
!> factor U or L is stored in the upper or lower triangular part of
!> the array A.
!>
!> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
!> overwriting the factor U in A.
!> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
!> overwriting the factor L in A.
!>
!> This is the unblocked form of the algorithm, calling Level 2 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the triangular factor stored in the array A
!>          is upper or lower triangular:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the triangular factor U or L.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the triangular factor U or L.
!>          On exit, if UPLO = 'U', the upper triangle of A is
!>          overwritten with the upper triangle of the product U * U**T;
!>          if UPLO = 'L', the lower triangle of A is overwritten with
!>          the lower triangle of the product L**T * L.
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
!>          < 0: if INFO = -k, the k-th argument had an illegal value
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
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLAUU2(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
!*--SLAUU2106
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , N
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER i
      REAL aii
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SDOT
      EXTERNAL LSAME , SDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEMV , SSCAL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLAUU2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Compute the product U * U**T.
!
         DO i = 1 , N
            aii = A(i,i)
            IF ( i<N ) THEN
               A(i,i) = SDOT(N-i+1,A(i,i),Lda,A(i,i),Lda)
               CALL SGEMV('No transpose',i-1,N-i,ONE,A(1,i+1),Lda,      &
     &                    A(i,i+1),Lda,aii,A(1,i),1)
            ELSE
               CALL SSCAL(i,aii,A(1,i),1)
            ENDIF
         ENDDO
!
      ELSE
!
!        Compute the product L**T * L.
!
         DO i = 1 , N
            aii = A(i,i)
            IF ( i<N ) THEN
               A(i,i) = SDOT(N-i+1,A(i,i),1,A(i,i),1)
               CALL SGEMV('Transpose',N-i,i-1,ONE,A(i+1,1),Lda,A(i+1,i),&
     &                    1,aii,A(i,1),Lda)
            ELSE
               CALL SSCAL(i,aii,A(i,1),Lda)
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of SLAUU2
!
      END SUBROUTINE SLAUU2
