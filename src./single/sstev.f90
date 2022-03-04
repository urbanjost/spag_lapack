!*==sstev.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> SSTEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSTEV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstev.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstev.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstev.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ
!       INTEGER            INFO, LDZ, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSTEV computes all eigenvalues and, optionally, eigenvectors of a
!> real symmetric tridiagonal matrix A.
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
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          On entry, the n diagonal elements of the tridiagonal matrix
!>          A.
!>          On exit, if INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          On entry, the (n-1) subdiagonal elements of the tridiagonal
!>          matrix A, stored in elements 1 to N-1 of E.
!>          On exit, the contents of E are destroyed.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is REAL array, dimension (LDZ, N)
!>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
!>          eigenvectors of the matrix A, with the i-th column of Z
!>          holding the eigenvector associated with D(i).
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
!>          WORK is REAL array, dimension (max(1,2*N-2))
!>          If JOBZ = 'N', WORK is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the algorithm failed to converge; i
!>                off-diagonal elements of E did not converge to zero.
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
      SUBROUTINE SSTEV(Jobz,N,D,E,Z,Ldz,Work,Info)
      USE S_LSAME
      USE S_SLAMCH
      USE S_SLANST
      USE S_SSCAL
      USE S_SSTEQR
      USE S_SSTERF
      USE S_XERBLA
      IMPLICIT NONE
!*--SSTEV127
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobz
      INTEGER :: N
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: bignum , eps , rmax , rmin , safmin , sigma , smlnum ,    &
     &        tnrm
      INTEGER :: imax , iscale
      LOGICAL :: wantz
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
!
      Info = 0
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) THEN
         Info = -6
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSTEV ',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( N==1 ) THEN
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
      iscale = 0
      tnrm = SLANST('M',N,D,E)
      IF ( tnrm>ZERO .AND. tnrm<rmin ) THEN
         iscale = 1
         sigma = rmin/tnrm
      ELSEIF ( tnrm>rmax ) THEN
         iscale = 1
         sigma = rmax/tnrm
      ENDIF
      IF ( iscale==1 ) THEN
         CALL SSCAL(N,sigma,D,1)
         CALL SSCAL(N-1,sigma,E(1),1)
      ENDIF
!
!     For eigenvalues only, call SSTERF.  For eigenvalues and
!     eigenvectors, call SSTEQR.
!
      IF ( .NOT.wantz ) THEN
         CALL SSTERF(N,D,E,Info)
      ELSE
         CALL SSTEQR('I',N,D,E,Z,Ldz,Work,Info)
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
         CALL SSCAL(imax,ONE/sigma,D,1)
      ENDIF
!
!
!     End of SSTEV
!
      END SUBROUTINE SSTEV
