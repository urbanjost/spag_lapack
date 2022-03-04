!*==chetri2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CHETRI2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHETRI2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetri2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetri2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetri2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHETRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHETRI2 computes the inverse of a COMPLEX hermitian indefinite matrix
!> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
!> CHETRF. CHETRI2 set the LEADING DIMENSION of the workspace
!> before calling CHETRI2X that actually computes the inverse.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the block diagonal matrix D and the multipliers
!>          used to obtain the factor U or L as computed by CHETRF.
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
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by CHETRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N+NB+1)*(NB+3)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          WORK is size >= (N+NB+1)*(NB+3)
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>           calculates:
!>              - the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array,
!>              - and no error message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
!>               inverse could not be computed.
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
!> \date November 2017
!
!> \ingroup complexHEcomputational
!
!  =====================================================================
      SUBROUTINE CHETRI2(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE S_CHETRI
      USE S_CHETRI2X
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CHETRI2136
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: lquery , upper
      INTEGER :: minsize , nbmax
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
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
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1)
!     Get blocksize
      nbmax = ILAENV(1,'CHETRF',Uplo,N,-1,-1,-1)
      IF ( nbmax>=N ) THEN
         minsize = N
      ELSE
         minsize = (N+nbmax+1)*(nbmax+3)
      ENDIF
!
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Lwork<minsize .AND. .NOT.lquery ) THEN
         Info = -7
      ENDIF
!
!     Quick return if possible
!
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CHETRI2',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         Work(1) = minsize
         RETURN
      ENDIF
      IF ( N==0 ) RETURN
 
      IF ( nbmax>=N ) THEN
         CALL CHETRI(Uplo,N,A,Lda,Ipiv,Work,Info)
      ELSE
         CALL CHETRI2X(Uplo,N,A,Lda,Ipiv,Work,nbmax,Info)
      ENDIF
!
!     End of CHETRI2
!
      END SUBROUTINE CHETRI2
