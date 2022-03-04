!*==zheev_2stage.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> ZHEEV_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices</b>
!
!  @precisions fortran z -> s d c
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHEEV_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheev_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheev_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheev_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHEEV_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK,
!                                RWORK, INFO )
!
!       IMPLICIT NONE
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * ), W( * )
!       COMPLEX*16         A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHEEV_2STAGE computes all eigenvalues and, optionally, eigenvectors of a
!> complex Hermitian matrix A using the 2stage technique for
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
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA, N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the
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
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
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
!>                                   dimension = max(stage1,stage2) + (KD+1)*N + N
!>                                             = N*KD + N*max(KD+1,FACTOPTNB)
!>                                               + max(2*KD*KD, KD*NTHREADS)
!>                                               + (KD+1)*N + N
!>                                   where KD is the blocking size of the reduction,
!>                                   FACTOPTNB is the blocking used by the QR or LQ
!>                                   algorithm, usually FACTOPTNB=128 is a good choice
!>                                   NTHREADS is the number of threads used when
!>                                   openMP compilation is enabled, otherwise =1.
!>          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (max(1, 3*N-2))
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
!> \ingroup complex16HEeigen
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
      SUBROUTINE ZHEEV_2STAGE(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Rwork,Info)
!
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_DSCAL
      USE S_DSTERF
      USE S_ILAENV2STAGE
      USE S_LSAME
      USE S_XERBLA
      USE S_ZHETRD_2STAGE
      USE S_ZLANHE
      USE S_ZLASCL
      USE S_ZSTEQR
      USE S_ZUNGTR
      IMPLICIT NONE
!*--ZHEEV_2STAGE205
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D0,0.0D0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anrm , bignum , eps , rmax , rmin , safmin ,      &
     &                sigma , smlnum
      INTEGER :: ib , iinfo , imax , inde , indhous , indtau , indwrk , &
     &           iscale , kd , lhtrd , llwork , lwmin , lwtrd
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
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ENDIF
!
      IF ( Info==0 ) THEN
         kd = ILAENV2STAGE(1,'ZHETRD_2STAGE',Jobz,N,-1,-1,-1)
         ib = ILAENV2STAGE(2,'ZHETRD_2STAGE',Jobz,N,kd,-1,-1)
         lhtrd = ILAENV2STAGE(3,'ZHETRD_2STAGE',Jobz,N,kd,ib,-1)
         lwtrd = ILAENV2STAGE(4,'ZHETRD_2STAGE',Jobz,N,kd,ib,-1)
         lwmin = N + lhtrd + lwtrd
         Work(1) = lwmin
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) Info = -8
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHEEV_2STAGE ',-Info)
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
         W(1) = DBLE(A(1,1))
         Work(1) = 1
         IF ( wantz ) A(1,1) = CONE
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
      anrm = ZLANHE('M',Uplo,N,A,Lda,Rwork)
      iscale = 0
      IF ( anrm>ZERO .AND. anrm<rmin ) THEN
         iscale = 1
         sigma = rmin/anrm
      ELSEIF ( anrm>rmax ) THEN
         iscale = 1
         sigma = rmax/anrm
      ENDIF
      IF ( iscale==1 ) CALL ZLASCL(Uplo,0,0,ONE,sigma,N,N,A,Lda,Info)
!
!     Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.
!
      inde = 1
      indtau = 1
      indhous = indtau + N
      indwrk = indhous + lhtrd
      llwork = Lwork - indwrk + 1
!
      CALL ZHETRD_2STAGE(Jobz,Uplo,N,A,Lda,W,Rwork(inde),Work(indtau),  &
     &                   Work(indhous),lhtrd,Work(indwrk),llwork,iinfo)
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, first call
!     ZUNGTR to generate the unitary matrix, then call ZSTEQR.
!
      IF ( .NOT.wantz ) THEN
         CALL DSTERF(N,W,Rwork(inde),Info)
      ELSE
         CALL ZUNGTR(Uplo,N,A,Lda,Work(indtau),Work(indwrk),llwork,     &
     &               iinfo)
         indwrk = inde + N
         CALL ZSTEQR(Jobz,N,W,Rwork(inde),A,Lda,Rwork(indwrk),Info)
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
!     Set WORK(1) to optimal complex workspace size.
!
      Work(1) = lwmin
!
!
!     End of ZHEEV_2STAGE
!
      END SUBROUTINE ZHEEV_2STAGE
