!*==dsbevd_2stage.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> DSBEVD_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  @precisions fortran d -> s
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSBEVD_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbevd_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbevd_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbevd_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSBEVD_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ,
!                                 WORK, LWORK, IWORK, LIWORK, INFO )
!
!       IMPLICIT NONE
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, KD, LDAB, LDZ, LIWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSBEVD_2STAGE computes all the eigenvalues and, optionally, eigenvectors of
!> a real symmetric band matrix A using the 2stage technique for
!> the reduction to tridiagonal. If eigenvectors are desired, it uses
!> a divide and conquer algorithm.
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
!>          WORK is DOUBLE PRECISION array, dimension LWORK
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
!>                                   dimension = (2KD+1)*N + KD*NTHREADS + N
!>                                   where KD is the size of the band.
!>                                   NTHREADS is the number of threads used when
!>                                   openMP compilation is enabled, otherwise =1.
!>          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available.
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
!>          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.
!>          If JOBZ  = 'V' and N > 2, LIWORK must be at least 3 + 5*N.
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
!> \ingroup doubleOTHEReigen
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
      SUBROUTINE DSBEVD_2STAGE(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,     &
     &                         Lwork,Iwork,Liwork,Info)
!
      USE F77KINDS                        
      USE S_DGEMM
      USE S_DLACPY
      USE S_DLAMCH
      USE S_DLANSB
      USE S_DLASCL
      USE S_DSCAL
      USE S_DSTEDC
      USE S_DSTERF
      USE S_DSYTRD_SB2ST
      USE S_ILAENV2STAGE
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DSBEVD_2STAGE252
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
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
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anrm , bignum , eps , rmax , rmin , safmin ,      &
     &                sigma , smlnum
      INTEGER :: ib , iinfo , inde , indhous , indwk2 , indwrk ,        &
     &           iscale , lhtrd , liwmin , llwork , llwrk2 , lwmin ,    &
     &           lwtrd
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
      lquery = (Lwork==-1 .OR. Liwork==-1)
!
      Info = 0
      IF ( N<=1 ) THEN
         liwmin = 1
         lwmin = 1
      ELSE
         ib = ILAENV2STAGE(2,'DSYTRD_SB2ST',Jobz,N,Kd,-1,-1)
         lhtrd = ILAENV2STAGE(3,'DSYTRD_SB2ST',Jobz,N,Kd,ib,-1)
         lwtrd = ILAENV2STAGE(4,'DSYTRD_SB2ST',Jobz,N,Kd,ib,-1)
         IF ( wantz ) THEN
            liwmin = 3 + 5*N
            lwmin = 1 + 5*N + 2*N**2
         ELSE
            liwmin = 1
            lwmin = MAX(2*N,N+lhtrd+lwtrd)
         ENDIF
      ENDIF
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
         Work(1) = lwmin
         Iwork(1) = liwmin
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -11
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -13
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSBEVD_2STAGE',-Info)
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
         W(1) = Ab(1,1)
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
!     Call DSYTRD_SB2ST to reduce band symmetric matrix to tridiagonal form.
!
      inde = 1
      indhous = inde + N
      indwrk = indhous + lhtrd
      llwork = Lwork - indwrk + 1
      indwk2 = indwrk + N*N
      llwrk2 = Lwork - indwk2 + 1
!
      CALL DSYTRD_SB2ST("N",Jobz,Uplo,N,Kd,Ab,Ldab,W,Work(inde),        &
     &                  Work(indhous),lhtrd,Work(indwrk),llwork,iinfo)
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEDC.
!
      IF ( .NOT.wantz ) THEN
         CALL DSTERF(N,W,Work(inde),Info)
      ELSE
         CALL DSTEDC('I',N,W,Work(inde),Work(indwrk),N,Work(indwk2),    &
     &               llwrk2,Iwork,Liwork,Info)
         CALL DGEMM('N','N',N,N,N,ONE,Z,Ldz,Work(indwrk),N,ZERO,        &
     &              Work(indwk2),N)
         CALL DLACPY('A',N,N,Work(indwk2),N,Z,Ldz)
      ENDIF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      IF ( iscale==1 ) CALL DSCAL(N,ONE/sigma,W,1)
!
      Work(1) = lwmin
      Iwork(1) = liwmin
!
!     End of DSBEVD_2STAGE
!
      END SUBROUTINE DSBEVD_2STAGE
