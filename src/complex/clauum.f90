!*==clauum.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLAUUM computes the product UUH or LHL, where U and L are upper or lower triangular matrices (blocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAUUM + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clauum.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clauum.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clauum.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAUUM( UPLO, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAUUM computes the product U * U**H or L**H * L, where the triangular
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the triangular factor U or L.
!>          On exit, if UPLO = 'U', the upper triangle of A is
!>          overwritten with the upper triangle of the product U * U**H;
!>          if UPLO = 'L', the lower triangle of A is overwritten with
!>          the lower triangle of the product L**H * L.
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLAUUM(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
!*--CLAUUM106
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
      COMPLEX A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E+0)
      COMPLEX CONE
      PARAMETER (CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER i , ib , nb
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      EXTERNAL LSAME , ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMM , CHERK , CLAUU2 , CTRMM , XERBLA
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
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CLAUUM',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Determine the block size for this environment.
!
      nb = ILAENV(1,'CLAUUM',Uplo,N,-1,-1,-1)
!
      IF ( nb<=1 .OR. nb>=N ) THEN
!
!        Use unblocked code
!
         CALL CLAUU2(Uplo,N,A,Lda,Info)
!
!        Use blocked code
!
      ELSEIF ( upper ) THEN
!
!           Compute the product U * U**H.
!
         DO i = 1 , N , nb
            ib = MIN(nb,N-i+1)
            CALL CTRMM('Right','Upper','Conjugate transpose','Non-unit',&
     &                 i-1,ib,CONE,A(i,i),Lda,A(1,i),Lda)
            CALL CLAUU2('Upper',ib,A(i,i),Lda,Info)
            IF ( i+ib<=N ) THEN
               CALL CGEMM('No transpose','Conjugate transpose',i-1,ib,  &
     &                    N-i-ib+1,CONE,A(1,i+ib),Lda,A(i,i+ib),Lda,    &
     &                    CONE,A(1,i),Lda)
               CALL CHERK('Upper','No transpose',ib,N-i-ib+1,ONE,       &
     &                    A(i,i+ib),Lda,ONE,A(i,i),Lda)
            ENDIF
         ENDDO
      ELSE
!
!           Compute the product L**H * L.
!
         DO i = 1 , N , nb
            ib = MIN(nb,N-i+1)
            CALL CTRMM('Left','Lower','Conjugate transpose','Non-unit', &
     &                 ib,i-1,CONE,A(i,i),Lda,A(i,1),Lda)
            CALL CLAUU2('Lower',ib,A(i,i),Lda,Info)
            IF ( i+ib<=N ) THEN
               CALL CGEMM('Conjugate transpose','No transpose',ib,i-1,  &
     &                    N-i-ib+1,CONE,A(i+ib,i),Lda,A(i+ib,1),Lda,    &
     &                    CONE,A(i,1),Lda)
               CALL CHERK('Lower','Conjugate transpose',ib,N-i-ib+1,ONE,&
     &                    A(i+ib,i),Lda,ONE,A(i,i),Lda)
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of CLAUUM
!
      END SUBROUTINE CLAUUM
