!*==dsyev.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> DSYEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSYEV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyev.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyev.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyev.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYEV computes all eigenvalues and, optionally, eigenvectors of a
!> real symmetric matrix A.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only;
!>          = 'V':  Compute eigenvalues and eigenvectors.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
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
!>          A is DOUBLE PRECISION array, dimension (LDA, N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of A contains the
!>          upper triangular part of the matrix A.  If UPLO = 'L',
!>          the leading N-by-N lower triangular part of A contains
!>          the lower triangular part of the matrix A.
!>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!>          orthonormal eigenvectors of the matrix A.
!>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!>          or the upper triangle (if UPLO='U') of A, including the
!>          diagonal, is destroyed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= max(1,3*N-1).
!>          For optimal efficiency, LWORK >= (NB+2)*N,
!>          where NB is the blocksize for DSYTRD returned by ILAENV.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the algorithm failed to converge; i
!>                off-diagonal elements of an intermediate tridiagonal
!>                form did not converge to zero.
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
!> \ingroup doubleSYeigen
!
!  =====================================================================
      SUBROUTINE DSYEV(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_DLANSY
      USE S_DLASCL
      USE S_DORGTR
      USE S_DSCAL
      USE S_DSTEQR
      USE S_DSTERF
      USE S_DSYTRD
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DSYEV148
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anrm , bignum , eps , rmax , rmin , safmin ,      &
     &                sigma , smlnum
      INTEGER :: iinfo , imax , inde , indtau , indwrk , iscale ,       &
     &           llwork , lwkopt , nb
      LOGICAL :: lower , lquery , wantz
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
      wantz = LSAME(Jobz,'V')
      lower = LSAME(Uplo,'L')
      lquery = (Lwork==-1)
!
      Info = 0
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(lower .OR. LSAME(Uplo,'U')) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ENDIF
!
      IF ( Info==0 ) THEN
         nb = ILAENV(1,'DSYTRD',Uplo,N,-1,-1,-1)
         lwkopt = MAX(1,(nb+2)*N)
         Work(1) = lwkopt
!
         IF ( Lwork<MAX(1,3*N-1) .AND. .NOT.lquery ) Info = -8
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSYEV ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( N==1 ) THEN
         W(1) = A(1,1)
         Work(1) = 2
         IF ( wantz ) A(1,1) = ONE
         RETURN
      ENDIF
!
!     Get machine constants.
!
      safmin = DLAMCH('Safe minimum')
      eps = DLAMCH('Precision')
      smlnum = safmin/eps
      bignum = ONE/smlnum
      rmin = SQRT(smlnum)
      rmax = SQRT(bignum)
!
!     Scale matrix to allowable range, if necessary.
!
      anrm = DLANSY('M',Uplo,N,A,Lda,Work)
      iscale = 0
      IF ( anrm>ZERO .AND. anrm<rmin ) THEN
         iscale = 1
         sigma = rmin/anrm
      ELSEIF ( anrm>rmax ) THEN
         iscale = 1
         sigma = rmax/anrm
      ENDIF
      IF ( iscale==1 ) CALL DLASCL(Uplo,0,0,ONE,sigma,N,N,A,Lda,Info)
!
!     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
!
      inde = 1
      indtau = inde + N
      indwrk = indtau + N
      llwork = Lwork - indwrk + 1
      CALL DSYTRD(Uplo,N,A,Lda,W,Work(inde),Work(indtau),Work(indwrk),  &
     &            llwork,iinfo)
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, first call
!     DORGTR to generate the orthogonal matrix, then call DSTEQR.
!
      IF ( .NOT.wantz ) THEN
         CALL DSTERF(N,W,Work(inde),Info)
      ELSE
         CALL DORGTR(Uplo,N,A,Lda,Work(indtau),Work(indwrk),llwork,     &
     &               iinfo)
         CALL DSTEQR(Jobz,N,W,Work(inde),A,Lda,Work(indtau),Info)
      ENDIF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      IF ( iscale==1 ) THEN
         IF ( Info==0 ) THEN
            imax = N
         ELSE
            imax = Info - 1
         ENDIF
         CALL DSCAL(imax,ONE/sigma,W,1)
      ENDIF
!
!     Set WORK(1) to optimal workspace size.
!
      Work(1) = lwkopt
!
!
!     End of DSYEV
!
      END SUBROUTINE DSYEV
