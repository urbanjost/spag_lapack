!*==dstevd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> DSTEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSTEVD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstevd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstevd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstevd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSTEVD( JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK,
!                          LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ
!       INTEGER            INFO, LDZ, LIWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSTEVD computes all eigenvalues and, optionally, eigenvectors of a
!> real symmetric tridiagonal matrix. If eigenvectors are desired, it
!> uses a divide and conquer algorithm.
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
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the n diagonal elements of the tridiagonal matrix
!>          A.
!>          On exit, if INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          On entry, the (n-1) subdiagonal elements of the tridiagonal
!>          matrix A, stored in elements 1 to N-1 of E.
!>          On exit, the contents of E are destroyed.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
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
!>          WORK is DOUBLE PRECISION array,
!>                                         dimension (LWORK)
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If JOBZ  = 'N' or N <= 1 then LWORK must be at least 1.
!>          If JOBZ  = 'V' and N > 1 then LWORK must be at least
!>                         ( 1 + 4*N + N**2 ).
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal sizes of the WORK and IWORK
!>          arrays, returns these values as the first entries of the WORK
!>          and IWORK arrays, and no error message related to LWORK or
!>          LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If JOBZ  = 'N' or N <= 1 then LIWORK must be at least 1.
!>          If JOBZ  = 'V' and N > 1 then LIWORK must be at least 3+5*N.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK and IWORK arrays, and no error message related to
!>          LWORK or LIWORK is issued by XERBLA.
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
!> \ingroup doubleOTHEReigen
!
!  =====================================================================
      SUBROUTINE DSTEVD(Jobz,N,D,E,Z,Ldz,Work,Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!*--DSTEVD166
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobz
      INTEGER Info , Ldz , Liwork , Lwork , N
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      DOUBLE PRECISION D(*) , E(*) , Work(*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , wantz
      INTEGER iscale , liwmin , lwmin
      DOUBLE PRECISION bignum , eps , rmax , rmin , safmin , sigma ,    &
     &                 smlnum , tnrm
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , DLANST
      EXTERNAL LSAME , DLAMCH , DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL DSCAL , DSTEDC , DSTERF , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      wantz = LSAME(Jobz,'V')
      lquery = (Lwork==-1 .OR. Liwork==-1)
!
      Info = 0
      liwmin = 1
      lwmin = 1
      IF ( N>1 .AND. wantz ) THEN
         lwmin = 1 + 4*N + N**2
         liwmin = 3 + 5*N
      ENDIF
!
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) THEN
         Info = -6
      ENDIF
!
      IF ( Info==0 ) THEN
         Work(1) = lwmin
         Iwork(1) = liwmin
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -8
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -10
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSTEVD',-Info)
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
      iscale = 0
      tnrm = DLANST('M',N,D,E)
      IF ( tnrm>ZERO .AND. tnrm<rmin ) THEN
         iscale = 1
         sigma = rmin/tnrm
      ELSEIF ( tnrm>rmax ) THEN
         iscale = 1
         sigma = rmax/tnrm
      ENDIF
      IF ( iscale==1 ) THEN
         CALL DSCAL(N,sigma,D,1)
         CALL DSCAL(N-1,sigma,E(1),1)
      ENDIF
!
!     For eigenvalues only, call DSTERF.  For eigenvalues and
!     eigenvectors, call DSTEDC.
!
      IF ( .NOT.wantz ) THEN
         CALL DSTERF(N,D,E,Info)
      ELSE
         CALL DSTEDC('I',N,D,E,Z,Ldz,Work,Lwork,Iwork,Liwork,Info)
      ENDIF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      IF ( iscale==1 ) CALL DSCAL(N,ONE/sigma,D,1)
!
      Work(1) = lwmin
      Iwork(1) = liwmin
!
!
!     End of DSTEVD
!
      END SUBROUTINE DSTEVD
