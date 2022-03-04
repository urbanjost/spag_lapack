!*==cheevr_2stage.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> CHEEVR_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices</b>
!
!  @generated from zheevr_2stage.f, fortran z -> c, Sat Nov  5 23:18:11 2016
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHEEVR_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cheevr_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cheevr_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cheevr_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHEEVR_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU,
!                                 IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ,
!                                 WORK, LWORK, RWORK, LRWORK, IWORK,
!                                 LIWORK, INFO )
!
!       IMPLICIT NONE
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, RANGE, UPLO
!       INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LRWORK, LWORK,
!      $                   M, N
!       REAL               ABSTOL, VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            ISUPPZ( * ), IWORK( * )
!       REAL               RWORK( * ), W( * )
!       COMPLEX            A( LDA, * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHEEVR_2STAGE computes selected eigenvalues and, optionally, eigenvectors
!> of a complex Hermitian matrix A using the 2stage technique for
!> the reduction to tridiagonal.  Eigenvalues and eigenvectors can
!> be selected by specifying either a range of values or a range of
!> indices for the desired eigenvalues.
!>
!> CHEEVR_2STAGE first reduces the matrix A to tridiagonal form T with a call
!> to CHETRD.  Then, whenever possible, CHEEVR_2STAGE calls CSTEMR to compute
!> eigenspectrum using Relatively Robust Representations.  CSTEMR
!> computes eigenvalues by the dqds algorithm, while orthogonal
!> eigenvectors are computed from various "good" L D L^T representations
!> (also known as Relatively Robust Representations). Gram-Schmidt
!> orthogonalization is avoided as far as possible. More specifically,
!> the various steps of the algorithm are as follows.
!>
!> For each unreduced block (submatrix) of T,
!>    (a) Compute T - sigma I  = L D L^T, so that L and D
!>        define all the wanted eigenvalues to high relative accuracy.
!>        This means that small relative changes in the entries of D and L
!>        cause only small relative changes in the eigenvalues and
!>        eigenvectors. The standard (unfactored) representation of the
!>        tridiagonal matrix T does not have this property in general.
!>    (b) Compute the eigenvalues to suitable accuracy.
!>        If the eigenvectors are desired, the algorithm attains full
!>        accuracy of the computed eigenvalues only right before
!>        the corresponding vectors have to be computed, see steps c) and d).
!>    (c) For each cluster of close eigenvalues, select a new
!>        shift close to the cluster, find a new factorization, and refine
!>        the shifted eigenvalues to suitable accuracy.
!>    (d) For each eigenvalue with a large enough relative separation compute
!>        the corresponding eigenvector by forming a rank revealing twisted
!>        factorization. Go back to (c) for any clusters that remain.
!>
!> The desired accuracy of the output can be specified by the input
!> parameter ABSTOL.
!>
!> For more details, see DSTEMR's documentation and:
!> - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
!>   to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
!>   Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
!> - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
!>   Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
!>   2004.  Also LAPACK Working Note 154.
!> - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
!>   tridiagonal eigenvalue/eigenvector problem",
!>   Computer Science Division Technical Report No. UCB/CSD-97-971,
!>   UC Berkeley, May 1997.
!>
!>
!> Note 1 : CHEEVR_2STAGE calls CSTEMR when the full spectrum is requested
!> on machines which conform to the ieee-754 floating point standard.
!> CHEEVR_2STAGE calls SSTEBZ and CSTEIN on non-ieee machines and
!> when partial spectrum requests are made.
!>
!> Normal execution of CSTEMR may create NaNs and infinities and
!> hence may abort due to a floating point exception in environments
!> which do not handle NaNs and infinities in the ieee standard default
!> manner.
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
!> \param[in] RANGE
!> \verbatim
!>          RANGE is CHARACTER*1
!>          = 'A': all eigenvalues will be found.
!>          = 'V': all eigenvalues in the half-open interval (VL,VU]
!>                 will be found.
!>          = 'I': the IL-th through IU-th eigenvalues will be found.
!>          For RANGE = 'V' or 'I' and IU - IL < N - 1, SSTEBZ and
!>          CSTEIN are called
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
!>          A is COMPLEX array, dimension (LDA, N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of A contains the
!>          upper triangular part of the matrix A.  If UPLO = 'L',
!>          the leading N-by-N lower triangular part of A contains
!>          the lower triangular part of the matrix A.
!>          On exit, the lower triangle (if UPLO='L') or the upper
!>          triangle (if UPLO='U') of A, including the diagonal, is
!>          destroyed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] VL
!> \verbatim
!>          VL is REAL
!>          If RANGE='V', the lower bound of the interval to
!>          be searched for eigenvalues. VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>          VU is REAL
!>          If RANGE='V', the upper bound of the interval to
!>          be searched for eigenvalues. VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] IL
!> \verbatim
!>          IL is INTEGER
!>          If RANGE='I', the index of the
!>          smallest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] IU
!> \verbatim
!>          IU is INTEGER
!>          If RANGE='I', the index of the
!>          largest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] ABSTOL
!> \verbatim
!>          ABSTOL is REAL
!>          The absolute error tolerance for the eigenvalues.
!>          An approximate eigenvalue is accepted as converged
!>          when it is determined to lie in an interval [a,b]
!>          of width less than or equal to
!>
!>                  ABSTOL + EPS *   max( |a|,|b| ) ,
!>
!>          where EPS is the machine precision.  If ABSTOL is less than
!>          or equal to zero, then  EPS*|T|  will be used in its place,
!>          where |T| is the 1-norm of the tridiagonal matrix obtained
!>          by reducing A to tridiagonal form.
!>
!>          See "Computing Small Singular Values of Bidiagonal Matrices
!>          with Guaranteed High Relative Accuracy," by Demmel and
!>          Kahan, LAPACK Working Note #3.
!>
!>          If high relative accuracy is important, set ABSTOL to
!>          SLAMCH( 'Safe minimum' ).  Doing so will guarantee that
!>          eigenvalues are computed to high relative accuracy when
!>          possible in future releases.  The current code does not
!>          make any guarantees about high relative accuracy, but
!>          future releases will. See J. Barlow and J. Demmel,
!>          "Computing Accurate Eigensystems of Scaled Diagonally
!>          Dominant Matrices", LAPACK Working Note #7, for a discussion
!>          of which matrices define their eigenvalues to high relative
!>          accuracy.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The total number of eigenvalues found.  0 <= M <= N.
!>          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is REAL array, dimension (N)
!>          The first M elements contain the selected eigenvalues in
!>          ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, max(1,M))
!>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!>          contain the orthonormal eigenvectors of the matrix A
!>          corresponding to the selected eigenvalues, with the i-th
!>          column of Z holding the eigenvector associated with W(i).
!>          If JOBZ = 'N', then Z is not referenced.
!>          Note: the user must ensure that at least max(1,M) columns are
!>          supplied in the array Z; if RANGE = 'V', the exact value of M
!>          is not known in advance and an upper bound must be used.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          JOBZ = 'V', LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] ISUPPZ
!> \verbatim
!>          ISUPPZ is INTEGER array, dimension ( 2*max(1,M) )
!>          The support of the eigenvectors in Z, i.e., the indices
!>          indicating the nonzero elements in Z. The i-th eigenvector
!>          is nonzero only in elements ISUPPZ( 2*i-1 ) through
!>          ISUPPZ( 2*i ). This is an output of CSTEMR (tridiagonal
!>          matrix). The support of the eigenvectors of A is typically
!>          1:N because of the unitary transformations applied by CUNMTR.
!>          Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If JOBZ = 'N' and N > 1, LWORK must be queried.
!>                                   LWORK = MAX(1, 26*N, dimension) where
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
!>          only calculates the optimal sizes of the WORK, RWORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (MAX(1,LRWORK))
!>          On exit, if INFO = 0, RWORK(1) returns the optimal
!>          (and minimal) LRWORK.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The length of the array RWORK.  LRWORK >= max(1,24*N).
!>
!>          If LRWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal
!>          (and minimal) LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.  LIWORK >= max(1,10*N).
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  Internal error
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
!> \date June 2016
!
!> \ingroup complexHEeigen
!
!> \par Contributors:
!  ==================
!>
!>     Inderjit Dhillon, IBM Almaden, USA \n
!>     Osni Marques, LBNL/NERSC, USA \n
!>     Ken Stanley, Computer Science Division, University of
!>       California at Berkeley, USA \n
!>     Jason Riedy, Computer Science Division, University of
!>       California at Berkeley, USA \n
!>
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
      SUBROUTINE CHEEVR_2STAGE(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,     &
     &                         Abstol,M,W,Z,Ldz,Isuppz,Work,Lwork,Rwork,&
     &                         Lrwork,Iwork,Liwork,Info)
!
      USE S_CHETRD_2STAGE
      USE S_CLANSY
      USE S_CSSCAL
      USE S_CSTEIN
      USE S_CSTEMR
      USE S_CSWAP
      USE S_CUNMTR
      USE S_ILAENV
      USE S_ILAENV2STAGE
      USE S_LSAME
      USE S_SCOPY
      USE S_SLAMCH
      USE S_SSCAL
      USE S_SSTEBZ
      USE S_SSTERF
      USE S_XERBLA
      IMPLICIT NONE
!*--CHEEVR_2STAGE426
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: abstll , anrm , bignum , eps , rmax , rmin , safmin ,     &
     &        sigma , smlnum , tmp1 , vll , vuu
      LOGICAL :: alleig , indeig , lower , lquery , test , tryrac ,     &
     &           valeig , wantz
      INTEGER :: i , ib , ieeeok , iinfo , imax , indhous , indibl ,    &
     &           indifl , indisp , indiwo , indrd , indrdd , indre ,    &
     &           indree , indrwk , indtau , indwk , indwkn , iscale ,   &
     &           itmp1 , j , jj , kd , lhtrd , liwmin , llrwork ,       &
     &           llwork , llwrkn , lrwmin , lwmin , lwtrd , nsplit
      CHARACTER :: order
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
      ieeeok = ILAENV(10,'CHEEVR','N',1,2,3,4)
!
      lower = LSAME(Uplo,'L')
      wantz = LSAME(Jobz,'V')
      alleig = LSAME(Range,'A')
      valeig = LSAME(Range,'V')
      indeig = LSAME(Range,'I')
!
      lquery = ((Lwork==-1) .OR. (Lrwork==-1) .OR. (Liwork==-1))
!
      kd = ILAENV2STAGE(1,'CHETRD_2STAGE',Jobz,N,-1,-1,-1)
      ib = ILAENV2STAGE(2,'CHETRD_2STAGE',Jobz,N,kd,-1,-1)
      lhtrd = ILAENV2STAGE(3,'CHETRD_2STAGE',Jobz,N,kd,ib,-1)
      lwtrd = ILAENV2STAGE(4,'CHETRD_2STAGE',Jobz,N,kd,ib,-1)
      lwmin = N + lhtrd + lwtrd
      lrwmin = MAX(1,24*N)
      liwmin = MAX(1,10*N)
!
      Info = 0
      IF ( .NOT.(LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(alleig .OR. valeig .OR. indeig) ) THEN
         Info = -2
      ELSEIF ( .NOT.(lower .OR. LSAME(Uplo,'U')) ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( valeig ) THEN
         IF ( N>0 .AND. Vu<=Vl ) Info = -8
      ELSEIF ( indeig ) THEN
         IF ( Il<1 .OR. Il>MAX(1,N) ) THEN
            Info = -9
         ELSEIF ( Iu<MIN(N,Il) .OR. Iu>N ) THEN
            Info = -10
         ENDIF
      ENDIF
      IF ( Info==0 ) THEN
         IF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) Info = -15
      ENDIF
!
      IF ( Info==0 ) THEN
         Work(1) = lwmin
         Rwork(1) = lrwmin
         Iwork(1) = liwmin
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -18
         ELSEIF ( Lrwork<lrwmin .AND. .NOT.lquery ) THEN
            Info = -20
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -22
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CHEEVR_2STAGE',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      M = 0
      IF ( N==0 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
      IF ( N==1 ) THEN
         Work(1) = 2
         IF ( alleig .OR. indeig ) THEN
            M = 1
            W(1) = REAL(A(1,1))
         ELSEIF ( Vl<REAL(A(1,1)) .AND. Vu>=REAL(A(1,1)) ) THEN
            M = 1
            W(1) = REAL(A(1,1))
         ENDIF
         IF ( wantz ) THEN
            Z(1,1) = ONE
            Isuppz(1) = 1
            Isuppz(2) = 1
         ENDIF
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
      rmax = MIN(SQRT(bignum),ONE/SQRT(SQRT(safmin)))
!
!     Scale matrix to allowable range, if necessary.
!
      iscale = 0
      abstll = Abstol
      IF ( valeig ) THEN
         vll = Vl
         vuu = Vu
      ENDIF
      anrm = CLANSY('M',Uplo,N,A,Lda,Rwork)
      IF ( anrm>ZERO .AND. anrm<rmin ) THEN
         iscale = 1
         sigma = rmin/anrm
      ELSEIF ( anrm>rmax ) THEN
         iscale = 1
         sigma = rmax/anrm
      ENDIF
      IF ( iscale==1 ) THEN
         IF ( lower ) THEN
            DO j = 1 , N
               CALL CSSCAL(N-j+1,sigma,A(j,j),1)
            ENDDO
         ELSE
            DO j = 1 , N
               CALL CSSCAL(j,sigma,A(1,j),1)
            ENDDO
         ENDIF
         IF ( Abstol>0 ) abstll = Abstol*sigma
         IF ( valeig ) THEN
            vll = Vl*sigma
            vuu = Vu*sigma
         ENDIF
      ENDIF
 
!     Initialize indices into workspaces.  Note: The IWORK indices are
!     used only if SSTERF or CSTEMR fail.
 
!     WORK(INDTAU:INDTAU+N-1) stores the complex scalar factors of the
!     elementary reflectors used in CHETRD.
      indtau = 1
!     INDWK is the starting offset of the remaining complex workspace,
!     and LLWORK is the remaining complex workspace size.
      indhous = indtau + N
      indwk = indhous + lhtrd
      llwork = Lwork - indwk + 1
 
!     RWORK(INDRD:INDRD+N-1) stores the real tridiagonal's diagonal
!     entries.
      indrd = 1
!     RWORK(INDRE:INDRE+N-1) stores the off-diagonal entries of the
!     tridiagonal matrix from CHETRD.
      indre = indrd + N
!     RWORK(INDRDD:INDRDD+N-1) is a copy of the diagonal entries over
!     -written by CSTEMR (the SSTERF path copies the diagonal to W).
      indrdd = indre + N
!     RWORK(INDREE:INDREE+N-1) is a copy of the off-diagonal entries over
!     -written while computing the eigenvalues in SSTERF and CSTEMR.
      indree = indrdd + N
!     INDRWK is the starting offset of the left-over real workspace, and
!     LLRWORK is the remaining workspace size.
      indrwk = indree + N
      llrwork = Lrwork - indrwk + 1
 
!     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in SSTEBZ and
!     stores the block indices of each of the M<=N eigenvalues.
      indibl = 1
!     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in SSTEBZ and
!     stores the starting and finishing indices of each block.
      indisp = indibl + N
!     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
!     that corresponding to eigenvectors that fail to converge in
!     CSTEIN.  This information is discarded; if any fail, the driver
!     returns INFO > 0.
      indifl = indisp + N
!     INDIWO is the offset of the remaining integer workspace.
      indiwo = indifl + N
 
!
!     Call CHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.
!
      CALL CHETRD_2STAGE(Jobz,Uplo,N,A,Lda,Rwork(indrd),Rwork(indre),   &
     &                   Work(indtau),Work(indhous),lhtrd,Work(indwk),  &
     &                   llwork,iinfo)
!
!     If all eigenvalues are desired
!     then call SSTERF or CSTEMR and CUNMTR.
!
      test = .FALSE.
      IF ( indeig ) THEN
         IF ( Il==1 .AND. Iu==N ) test = .TRUE.
      ENDIF
      IF ( (alleig .OR. test) .AND. (ieeeok==1) ) THEN
         IF ( .NOT.wantz ) THEN
            CALL SCOPY(N,Rwork(indrd),1,W,1)
            CALL SCOPY(N-1,Rwork(indre),1,Rwork(indree),1)
            CALL SSTERF(N,W,Rwork(indree),Info)
         ELSE
            CALL SCOPY(N-1,Rwork(indre),1,Rwork(indree),1)
            CALL SCOPY(N,Rwork(indrd),1,Rwork(indrdd),1)
!
            IF ( Abstol<=TWO*N*eps ) THEN
               tryrac = .TRUE.
            ELSE
               tryrac = .FALSE.
            ENDIF
            CALL CSTEMR(Jobz,'A',N,Rwork(indrdd),Rwork(indree),Vl,Vu,Il,&
     &                  Iu,M,W,Z,Ldz,N,Isuppz,tryrac,Rwork(indrwk),     &
     &                  llrwork,Iwork,Liwork,Info)
!
!           Apply unitary matrix used in reduction to tridiagonal
!           form to eigenvectors returned by CSTEMR.
!
            IF ( wantz .AND. Info==0 ) THEN
               indwkn = indwk
               llwrkn = Lwork - indwkn + 1
               CALL CUNMTR('L',Uplo,'N',N,M,A,Lda,Work(indtau),Z,Ldz,   &
     &                     Work(indwkn),llwrkn,iinfo)
            ENDIF
         ENDIF
!
!
         IF ( Info==0 ) THEN
            M = N
            GOTO 100
         ENDIF
         Info = 0
      ENDIF
!
!     Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN.
!     Also call SSTEBZ and CSTEIN if CSTEMR fails.
!
      IF ( wantz ) THEN
         order = 'B'
      ELSE
         order = 'E'
      ENDIF
 
      CALL SSTEBZ(Range,order,N,vll,vuu,Il,Iu,abstll,Rwork(indrd),      &
     &            Rwork(indre),M,nsplit,W,Iwork(indibl),Iwork(indisp),  &
     &            Rwork(indrwk),Iwork(indiwo),Info)
!
      IF ( wantz ) THEN
         CALL CSTEIN(N,Rwork(indrd),Rwork(indre),M,W,Iwork(indibl),     &
     &               Iwork(indisp),Z,Ldz,Rwork(indrwk),Iwork(indiwo),   &
     &               Iwork(indifl),Info)
!
!        Apply unitary matrix used in reduction to tridiagonal
!        form to eigenvectors returned by CSTEIN.
!
         indwkn = indwk
         llwrkn = Lwork - indwkn + 1
         CALL CUNMTR('L',Uplo,'N',N,M,A,Lda,Work(indtau),Z,Ldz,         &
     &               Work(indwkn),llwrkn,iinfo)
      ENDIF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
 100  IF ( iscale==1 ) THEN
         IF ( Info==0 ) THEN
            imax = M
         ELSE
            imax = Info - 1
         ENDIF
         CALL SSCAL(imax,ONE/sigma,W,1)
      ENDIF
!
!     If eigenvalues are not in order, then sort them, along with
!     eigenvectors.
!
      IF ( wantz ) THEN
         DO j = 1 , M - 1
            i = 0
            tmp1 = W(j)
            DO jj = j + 1 , M
               IF ( W(jj)<tmp1 ) THEN
                  i = jj
                  tmp1 = W(jj)
               ENDIF
            ENDDO
!
            IF ( i/=0 ) THEN
               itmp1 = Iwork(indibl+i-1)
               W(i) = W(j)
               Iwork(indibl+i-1) = Iwork(indibl+j-1)
               W(j) = tmp1
               Iwork(indibl+j-1) = itmp1
               CALL CSWAP(N,Z(1,i),1,Z(1,j),1)
            ENDIF
         ENDDO
      ENDIF
!
!     Set WORK(1) to optimal workspace size.
!
      Work(1) = lwmin
      Rwork(1) = lrwmin
      Iwork(1) = liwmin
!
!
!     End of CHEEVR_2STAGE
!
      END SUBROUTINE CHEEVR_2STAGE
