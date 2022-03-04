!*==sspevd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> SSPEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSPEVD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspevd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspevd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspevd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK,
!                          IWORK, LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, LDZ, LIWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               AP( * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSPEVD computes all the eigenvalues and, optionally, eigenvectors
!> of a real symmetric matrix A in packed storage. If eigenvectors are
!> desired, it uses a divide and conquer algorithm.
!>
!> The divide and conquer algorithm makes very mild assumptions about
!> floating point arithmetic. It will work on machines with a guard
!> digit in add/subtract, or on those binary machines without guard
!> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!> Cray-2. It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.
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
!> \param[in,out] AP
!> \verbatim
!>          AP is REAL array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the symmetric matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!>
!>          On exit, AP is overwritten by values generated during the
!>          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
!>          and first superdiagonal of the tridiagonal matrix T overwrite
!>          the corresponding elements of A, and if UPLO = 'L', the
!>          diagonal and first subdiagonal of T overwrite the
!>          corresponding elements of A.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is REAL array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is REAL array, dimension (LDZ, N)
!>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
!>          eigenvectors of the matrix A, with the i-th column of Z
!>          holding the eigenvector associated with W(i).
!>          If JOBZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          JOBZ = 'V', LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the required LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If N <= 1,               LWORK must be at least 1.
!>          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N.
!>          If JOBZ = 'V' and N > 1, LWORK must be at least
!>                                                 1 + 6*N + N**2.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the required sizes of the WORK and IWORK
!>          arrays, returns these values as the first entries of the WORK
!>          and IWORK arrays, and no error message related to LWORK or
!>          LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the required LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.
!>          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the required sizes of the WORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK and IWORK arrays, and no error message related to
!>          LWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup realOTHEReigen
!
!  =====================================================================
      SUBROUTINE SSPEVD(Jobz,Uplo,N,Ap,W,Z,Ldz,Work,Lwork,Iwork,Liwork, &
     &                  Info)
      USE S_LSAME
      USE S_SLAMCH
      USE S_SLANSP
      USE S_SOPMTR
      USE S_SSCAL
      USE S_SSPTRD
      USE S_SSTEDC
      USE S_SSTERF
      USE S_XERBLA
      IMPLICIT NONE
!*--SSPEVD191
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: anrm , bignum , eps , rmax , rmin , safmin , sigma ,      &
     &        smlnum
      INTEGER :: iinfo , inde , indtau , indwrk , iscale , liwmin ,     &
     &           llwork , lwmin
      LOGICAL :: lquery , wantz
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
      lquery = (Lwork==-1 .OR. Liwork==-1)
!
      Info = 0
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(LSAME(Uplo,'U') .OR. LSAME(Uplo,'L')) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) THEN
         Info = -7
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( N<=1 ) THEN
            liwmin = 1
            lwmin = 1
         ELSEIF ( wantz ) THEN
            liwmin = 3 + 5*N
            lwmin = 1 + 6*N + N**2
         ELSE
            liwmin = 1
            lwmin = 2*N
         ENDIF
         Iwork(1) = liwmin
         Work(1) = lwmin
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -9
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -11
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSPEVD',-Info)
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
         W(1) = Ap(1)
         IF ( wantz ) Z(1,1) = ONE
         RETURN
      ENDIF
!
!     Get machine constants.
!
      safmin = SLAMCH('Safe minimum')
      eps = SLAMCH('Precision')
      smlnum = safmin/eps
      bignum = ONE/smlnum
      rmin = SQRT(smlnum)
      rmax = SQRT(bignum)
!
!     Scale matrix to allowable range, if necessary.
!
      anrm = SLANSP('M',Uplo,N,Ap,Work)
      iscale = 0
      IF ( anrm>ZERO .AND. anrm<rmin ) THEN
         iscale = 1
         sigma = rmin/anrm
      ELSEIF ( anrm>rmax ) THEN
         iscale = 1
         sigma = rmax/anrm
      ENDIF
      IF ( iscale==1 ) CALL SSCAL((N*(N+1))/2,sigma,Ap,1)
!
!     Call SSPTRD to reduce symmetric packed matrix to tridiagonal form.
!
      inde = 1
      indtau = inde + N
      CALL SSPTRD(Uplo,N,Ap,W,Work(inde),Work(indtau),iinfo)
!
!     For eigenvalues only, call SSTERF.  For eigenvectors, first call
!     SSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
!     tridiagonal matrix, then call SOPMTR to multiply it by the
!     Householder transformations represented in AP.
!
      IF ( .NOT.wantz ) THEN
         CALL SSTERF(N,W,Work(inde),Info)
      ELSE
         indwrk = indtau + N
         llwork = Lwork - indwrk + 1
         CALL SSTEDC('I',N,W,Work(inde),Z,Ldz,Work(indwrk),llwork,Iwork,&
     &               Liwork,Info)
         CALL SOPMTR('L',Uplo,'N',N,N,Ap,Work(indtau),Z,Ldz,Work(indwrk)&
     &               ,iinfo)
      ENDIF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      IF ( iscale==1 ) CALL SSCAL(N,ONE/sigma,W,1)
!
      Work(1) = lwmin
      Iwork(1) = liwmin
!
!     End of SSPEVD
!
      END SUBROUTINE SSPEVD
