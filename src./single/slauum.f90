!*==slauum.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLAUUM computes the product UUH or LHL, where U and L are upper or lower triangular matrices (blocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAUUM + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slauum.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slauum.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slauum.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAUUM( UPLO, N, A, LDA, INFO )
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
!> SLAUUM computes the product U * U**T or L**T * L, where the triangular
!> factor U or L is stored in the upper or lower triangular part of
!> the array A.
!>
!> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
!> overwriting the factor U in A.
!> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
!> overwriting the factor L in A.
!>
!> This is the blocked form of the algorithm, calling Level 3 BLAS.
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
      SUBROUTINE SLAUUM(Uplo,N,A,Lda,Info)
      USE S_ILAENV
      USE S_LSAME
      USE S_SGEMM
      USE S_SLAUU2
      USE S_SSYRK
      USE S_STRMM
      USE S_XERBLA
      IMPLICIT NONE
!*--SLAUUM113
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ib , nb
      LOGICAL :: upper
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
!     .. Intrinsic Functions ..
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
         CALL XERBLA('SLAUUM',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Determine the block size for this environment.
!
      nb = ILAENV(1,'SLAUUM',Uplo,N,-1,-1,-1)
!
      IF ( nb<=1 .OR. nb>=N ) THEN
!
!        Use unblocked code
!
         CALL SLAUU2(Uplo,N,A,Lda,Info)
!
!        Use blocked code
!
      ELSEIF ( upper ) THEN
!
!           Compute the product U * U**T.
!
         DO i = 1 , N , nb
            ib = MIN(nb,N-i+1)
            CALL STRMM('Right','Upper','Transpose','Non-unit',i-1,ib,   &
     &                 ONE,A(i,i),Lda,A(1,i),Lda)
            CALL SLAUU2('Upper',ib,A(i,i),Lda,Info)
            IF ( i+ib<=N ) THEN
               CALL SGEMM('No transpose','Transpose',i-1,ib,N-i-ib+1,   &
     &                    ONE,A(1,i+ib),Lda,A(i,i+ib),Lda,ONE,A(1,i),   &
     &                    Lda)
               CALL SSYRK('Upper','No transpose',ib,N-i-ib+1,ONE,       &
     &                    A(i,i+ib),Lda,ONE,A(i,i),Lda)
            ENDIF
         ENDDO
      ELSE
!
!           Compute the product L**T * L.
!
         DO i = 1 , N , nb
            ib = MIN(nb,N-i+1)
            CALL STRMM('Left','Lower','Transpose','Non-unit',ib,i-1,ONE,&
     &                 A(i,i),Lda,A(i,1),Lda)
            CALL SLAUU2('Lower',ib,A(i,i),Lda,Info)
            IF ( i+ib<=N ) THEN
               CALL SGEMM('Transpose','No transpose',ib,i-1,N-i-ib+1,   &
     &                    ONE,A(i+ib,i),Lda,A(i+ib,1),Lda,ONE,A(i,1),   &
     &                    Lda)
               CALL SSYRK('Lower','Transpose',ib,N-i-ib+1,ONE,A(i+ib,i),&
     &                    Lda,ONE,A(i,i),Lda)
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of SLAUUM
!
      END SUBROUTINE SLAUUM
