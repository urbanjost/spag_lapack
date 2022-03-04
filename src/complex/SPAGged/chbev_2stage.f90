!*==chbev_2stage.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> CHBEV_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  @generated from zhbev_2stage.f, fortran z -> c, Sat Nov  5 23:18:20 2016
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHBEV_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbev_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbev_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbev_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHBEV_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ,
!                                WORK, LWORK, RWORK, INFO )
!
!       IMPLICIT NONE
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, KD, LDAB, LDZ, N, LWORK
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * ), W( * )
!       COMPLEX            AB( LDAB, * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHBEV_2STAGE computes all the eigenvalues and, optionally, eigenvectors of
!> a complex Hermitian band matrix A using the 2stage technique for
!> the reduction to tridiagonal.
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
!>                  Not available in this release.
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
!>          AB is COMPLEX array, dimension (LDAB, N)
!>          On entry, the upper or lower triangle of the Hermitian band
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
!>          W is REAL array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, N)
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
!>          WORK is COMPLEX array, dimension LWORK
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK. LWORK >= 1, when N <= 1;
!>          otherwise
!>          If JOBZ = 'N' and N > 1, LWORK must be queried.
!>                                   LWORK = MAX(1, dimension) where
!>                                   dimension = (2KD+1)*N + KD*NTHREADS
!>                                   where KD is the size of the band.
!>                                   NTHREADS is the number of threads used when
!>                                   openMP compilation is enabled, otherwise =1.
!>          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal sizes of the WORK, RWORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (max(1,3*N-2))
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
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
!> \date November 2017
!
!> \ingroup complexOTHEReigen
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  All details about the 2stage techniques are available in:
!>
!>  Azzam Haidar, Hatem Ltaief, and Jack Dongarra.
!>  Parallel reduction to condensed forms for symmetric eigenvalue problems
!>  using aggregated fine-grained and memory-aware kernels. In Proceedings
!>  of 2011 International Conference for High Performance Computing,
!>  Networking, Storage and Analysis (SC '11), New York, NY, USA,
!>  Article 8 , 11 pages.
!>  http://doi.acm.org/10.1145/2063384.2063394
!>
!>  A. Haidar, J. Kurzak, P. Luszczek, 2013.
!>  An improved parallel singular value algorithm and its implementation
!>  for multicore hardware, In Proceedings of 2013 International Conference
!>  for High Performance Computing, Networking, Storage and Analysis (SC '13).
!>  Denver, Colorado, USA, 2013.
!>  Article 90, 12 pages.
!>  http://doi.acm.org/10.1145/2503210.2503292
!>
!>  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra.
!>  A novel hybrid CPU-GPU generalized eigensolver for electronic structure
!>  calculations based on fine-grained memory aware tasks.
!>  International Journal of High Performance Computing Applications.
!>  Volume 28 Issue 2, Pages 196-209, May 2014.
!>  http://hpc.sagepub.com/content/28/2/196
!>
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE CHBEV_2STAGE(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Lwork,&
     &                        Rwork,Info)
!
      USE S_CHETRD_2STAGE
      USE S_CHETRD_HB2ST
      USE S_CLANHB
      USE S_CLASCL
      USE S_CSTEQR
      USE S_ILAENV2STAGE
      USE S_LSAME
      USE S_SLAMCH
      USE S_SSCAL
      USE S_SSTERF
      USE S_XERBLA
      IMPLICIT NONE
!*--CHBEV_2STAGE227
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: anrm , bignum , eps , rmax , rmin , safmin , sigma ,      &
     &        smlnum
      INTEGER :: ib , iinfo , imax , inde , indhous , indrwk , indwrk , &
     &           iscale , lhtrd , llwork , lwmin , lwtrd
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
      IF ( .NOT.(LSAME(Jobz,'N')) ) THEN
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
      IF ( Info==0 ) THEN
         IF ( N<=1 ) THEN
            lwmin = 1
            Work(1) = lwmin
         ELSE
            ib = ILAENV2STAGE(2,'CHETRD_HB2ST',Jobz,N,Kd,-1,-1)
            lhtrd = ILAENV2STAGE(3,'CHETRD_HB2ST',Jobz,N,Kd,ib,-1)
            lwtrd = ILAENV2STAGE(4,'CHETRD_HB2ST',Jobz,N,Kd,ib,-1)
            lwmin = lhtrd + lwtrd
            Work(1) = lwmin
         ENDIF
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) Info = -11
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CHBEV_2STAGE ',-Info)
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
         IF ( lower ) THEN
            W(1) = REAL(Ab(1,1))
         ELSE
            W(1) = REAL(Ab(Kd+1,1))
         ENDIF
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
      anrm = CLANHB('M',Uplo,N,Kd,Ab,Ldab,Rwork)
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
            CALL CLASCL('B',Kd,Kd,ONE,sigma,N,N,Ab,Ldab,Info)
         ELSE
            CALL CLASCL('Q',Kd,Kd,ONE,sigma,N,N,Ab,Ldab,Info)
         ENDIF
      ENDIF
!
!     Call CHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form.
!
      inde = 1
      indhous = 1
      indwrk = indhous + lhtrd
      llwork = Lwork - indwrk + 1
!
      CALL CHETRD_HB2ST("N",Jobz,Uplo,N,Kd,Ab,Ldab,W,Rwork(inde),       &
     &                  Work(indhous),lhtrd,Work(indwrk),llwork,iinfo)
!
!     For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEQR.
!
      IF ( .NOT.wantz ) THEN
         CALL SSTERF(N,W,Rwork(inde),Info)
      ELSE
         indrwk = inde + N
         CALL CSTEQR(Jobz,N,W,Rwork(inde),Z,Ldz,Rwork(indrwk),Info)
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
         CALL SSCAL(imax,ONE/sigma,W,1)
      ENDIF
!
!     Set WORK(1) to optimal workspace size.
!
      Work(1) = lwmin
!
!
!     End of CHBEV_2STAGE
!
      END SUBROUTINE CHBEV_2STAGE
