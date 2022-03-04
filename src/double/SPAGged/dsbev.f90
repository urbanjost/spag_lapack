!*==dsbev.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> DSBEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSBEV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbev.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbev.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbev.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,
!                         INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, KD, LDAB, LDZ, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSBEV computes all the eigenvalues and, optionally, eigenvectors of
!> a real symmetric band matrix A.
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
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB, N)
!>          On entry, the upper or lower triangle of the symmetric band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>
!>          On exit, AB is overwritten by values generated during the
!>          reduction to tridiagonal form.  If UPLO = 'U', the first
!>          superdiagonal and the diagonal of the tridiagonal matrix T
!>          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',
!>          the diagonal and first subdiagonal of T are returned in the
!>          first two rows of AB.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD + 1.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
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
!>          WORK is DOUBLE PRECISION array, dimension (max(1,3*N-2))
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
!> \ingroup doubleOTHEReigen
!
!  =====================================================================
      SUBROUTINE DSBEV(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_DLANSB
      USE S_DLASCL
      USE S_DSBTRD
      USE S_DSCAL
      USE S_DSTEQR
      USE S_DSTERF
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DSBEV159
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
      INTEGER :: Kd
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anrm , bignum , eps , rmax , rmin , safmin ,      &
     &                sigma , smlnum
      INTEGER :: iinfo , imax , inde , indwrk , iscale
      LOGICAL :: lower , wantz
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
!
      Info = 0
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(lower .OR. LSAME(Uplo,'U')) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Kd<0 ) THEN
         Info = -4
      ELSEIF ( Ldab<Kd+1 ) THEN
         Info = -6
      ELSEIF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) THEN
         Info = -9
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSBEV ',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( N==1 ) THEN
         IF ( lower ) THEN
            W(1) = Ab(1,1)
         ELSE
            W(1) = Ab(Kd+1,1)
         ENDIF
         IF ( wantz ) Z(1,1) = ONE
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
      anrm = DLANSB('M',Uplo,N,Kd,Ab,Ldab,Work)
      iscale = 0
      IF ( anrm>ZERO .AND. anrm<rmin ) THEN
         iscale = 1
         sigma = rmin/anrm
      ELSEIF ( anrm>rmax ) THEN
         iscale = 1
         sigma = rmax/anrm
      ENDIF
      IF ( iscale==1 ) THEN
         IF ( lower ) THEN
            CALL DLASCL('B',Kd,Kd,ONE,sigma,N,N,Ab,Ldab,Info)
         ELSE
            CALL DLASCL('Q',Kd,Kd,ONE,sigma,N,N,Ab,Ldab,Info)
         ENDIF
      ENDIF
!
!     Call DSBTRD to reduce symmetric band matrix to tridiagonal form.
!
      inde = 1
      indwrk = inde + N
      CALL DSBTRD(Jobz,Uplo,N,Kd,Ab,Ldab,W,Work(inde),Z,Ldz,Work(indwrk)&
     &            ,iinfo)
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR.
!
      IF ( .NOT.wantz ) THEN
         CALL DSTERF(N,W,Work(inde),Info)
      ELSE
         CALL DSTEQR(Jobz,N,W,Work(inde),Z,Ldz,Work(indwrk),Info)
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
!
!     End of DSBEV
!
      END SUBROUTINE DSBEV
