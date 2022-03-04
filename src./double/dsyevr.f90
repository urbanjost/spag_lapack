!*==dsyevr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> DSYEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSYEVR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyevr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyevr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyevr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
!                          ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,
!                          IWORK, LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, RANGE, UPLO
!       INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N
!       DOUBLE PRECISION   ABSTOL, VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            ISUPPZ( * ), IWORK( * )
!       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYEVR computes selected eigenvalues and, optionally, eigenvectors
!> of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
!> selected by specifying either a range of values or a range of
!> indices for the desired eigenvalues.
!>
!> DSYEVR first reduces the matrix A to tridiagonal form T with a call
!> to DSYTRD.  Then, whenever possible, DSYEVR calls DSTEMR to compute
!> the eigenspectrum using Relatively Robust Representations.  DSTEMR
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
!> Note 1 : DSYEVR calls DSTEMR when the full spectrum is requested
!> on machines which conform to the ieee-754 floating point standard.
!> DSYEVR calls DSTEBZ and DSTEIN on non-ieee machines and
!> when partial spectrum requests are made.
!>
!> Normal execution of DSTEMR may create NaNs and infinities and
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
!> \endverbatim
!>
!> \param[in] RANGE
!> \verbatim
!>          RANGE is CHARACTER*1
!>          = 'A': all eigenvalues will be found.
!>          = 'V': all eigenvalues in the half-open interval (VL,VU]
!>                 will be found.
!>          = 'I': the IL-th through IU-th eigenvalues will be found.
!>          For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and
!>          DSTEIN are called
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
!>          VL is DOUBLE PRECISION
!>          If RANGE='V', the lower bound of the interval to
!>          be searched for eigenvalues. VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>          VU is DOUBLE PRECISION
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
!>          ABSTOL is DOUBLE PRECISION
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
!>          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that
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
!>          W is DOUBLE PRECISION array, dimension (N)
!>          The first M elements contain the selected eigenvalues in
!>          ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M))
!>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!>          contain the orthonormal eigenvectors of the matrix A
!>          corresponding to the selected eigenvalues, with the i-th
!>          column of Z holding the eigenvector associated with W(i).
!>          If JOBZ = 'N', then Z is not referenced.
!>          Note: the user must ensure that at least max(1,M) columns are
!>          supplied in the array Z; if RANGE = 'V', the exact value of M
!>          is not known in advance and an upper bound must be used.
!>          Supplying N columns is always safe.
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
!>          ISUPPZ( 2*i ). This is an output of DSTEMR (tridiagonal
!>          matrix). The support of the eigenvectors of A is typically
!>          1:N because of the orthogonal transformations applied by DORMTR.
!>          Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
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
!>          The dimension of the array WORK.  LWORK >= max(1,26*N).
!>          For optimal efficiency, LWORK >= (NB+6)*N,
!>          where NB is the max of the blocksize for DSYTRD and DORMTR
!>          returned by ILAENV.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.  LIWORK >= max(1,10*N).
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of the IWORK array,
!>          returns this value as the first entry of the IWORK array, and
!>          no error message related to LIWORK is issued by XERBLA.
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
!> \ingroup doubleSYeigen
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
!  =====================================================================
      SUBROUTINE DSYEVR(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,Abstol,M,W, &
     &                  Z,Ldz,Isuppz,Work,Lwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DLAMCH
      USE S_DLANSY
      USE S_DORMTR
      USE S_DSCAL
      USE S_DSTEBZ
      USE S_DSTEIN
      USE S_DSTEMR
      USE S_DSTERF
      USE S_DSWAP
      USE S_DSYTRD
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DSYEVR352
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: abstll , anrm , bignum , eps , rmax , rmin ,      &
     &                safmin , sigma , smlnum , tmp1 , vll , vuu
      LOGICAL :: alleig , indeig , lower , lquery , tryrac , valeig ,   &
     &           wantz
      INTEGER :: i , ieeeok , iinfo , imax , indd , inddd , inde ,      &
     &           indee , indibl , indifl , indisp , indiwo , indtau ,   &
     &           indwk , indwkn , iscale , j , jj , liwmin , llwork ,   &
     &           llwrkn , lwkopt , lwmin , nb , nsplit
      CHARACTER :: order
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
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
      ieeeok = ILAENV(10,'DSYEVR','N',1,2,3,4)
!
      lower = LSAME(Uplo,'L')
      wantz = LSAME(Jobz,'V')
      alleig = LSAME(Range,'A')
      valeig = LSAME(Range,'V')
      indeig = LSAME(Range,'I')
!
      lquery = ((Lwork==-1) .OR. (Liwork==-1))
!
      lwmin = MAX(1,26*N)
      liwmin = MAX(1,10*N)
!
      Info = 0
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
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
         IF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) THEN
            Info = -15
         ELSEIF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -18
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -20
         ENDIF
      ENDIF
!
      IF ( Info==0 ) THEN
         nb = ILAENV(1,'DSYTRD',Uplo,N,-1,-1,-1)
         nb = MAX(nb,ILAENV(1,'DORMTR',Uplo,N,-1,-1,-1))
         lwkopt = MAX((nb+1)*N,lwmin)
         Work(1) = lwkopt
         Iwork(1) = liwmin
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSYEVR',-Info)
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
         Work(1) = 7
         IF ( alleig .OR. indeig ) THEN
            M = 1
            W(1) = A(1,1)
         ELSEIF ( Vl<A(1,1) .AND. Vu>=A(1,1) ) THEN
            M = 1
            W(1) = A(1,1)
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
      safmin = DLAMCH('Safe minimum')
      eps = DLAMCH('Precision')
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
      anrm = DLANSY('M',Uplo,N,A,Lda,Work)
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
               CALL DSCAL(N-j+1,sigma,A(j,j),1)
            ENDDO
         ELSE
            DO j = 1 , N
               CALL DSCAL(j,sigma,A(1,j),1)
            ENDDO
         ENDIF
         IF ( Abstol>0 ) abstll = Abstol*sigma
         IF ( valeig ) THEN
            vll = Vl*sigma
            vuu = Vu*sigma
         ENDIF
      ENDIF
 
!     Initialize indices into workspaces.  Note: The IWORK indices are
!     used only if DSTERF or DSTEMR fail.
 
!     WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the
!     elementary reflectors used in DSYTRD.
      indtau = 1
!     WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries.
      indd = indtau + N
!     WORK(INDE:INDE+N-1) stores the off-diagonal entries of the
!     tridiagonal matrix from DSYTRD.
      inde = indd + N
!     WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over
!     -written by DSTEMR (the DSTERF path copies the diagonal to W).
      inddd = inde + N
!     WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over
!     -written while computing the eigenvalues in DSTERF and DSTEMR.
      indee = inddd + N
!     INDWK is the starting offset of the left-over workspace, and
!     LLWORK is the remaining workspace size.
      indwk = indee + N
      llwork = Lwork - indwk + 1
 
!     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and
!     stores the block indices of each of the M<=N eigenvalues.
      indibl = 1
!     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and
!     stores the starting and finishing indices of each block.
      indisp = indibl + N
!     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
!     that corresponding to eigenvectors that fail to converge in
!     DSTEIN.  This information is discarded; if any fail, the driver
!     returns INFO > 0.
      indifl = indisp + N
!     INDIWO is the offset of the remaining integer workspace.
      indiwo = indifl + N
 
!
!     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
!
      CALL DSYTRD(Uplo,N,A,Lda,Work(indd),Work(inde),Work(indtau),      &
     &            Work(indwk),llwork,iinfo)
!
!     If all eigenvalues are desired
!     then call DSTERF or DSTEMR and DORMTR.
!
      IF ( (alleig .OR. (indeig .AND. Il==1 .AND. Iu==N)) .AND.         &
     &     ieeeok==1 ) THEN
         IF ( .NOT.wantz ) THEN
            CALL DCOPY(N,Work(indd),1,W,1)
            CALL DCOPY(N-1,Work(inde),1,Work(indee),1)
            CALL DSTERF(N,W,Work(indee),Info)
         ELSE
            CALL DCOPY(N-1,Work(inde),1,Work(indee),1)
            CALL DCOPY(N,Work(indd),1,Work(inddd),1)
!
            IF ( Abstol<=TWO*N*eps ) THEN
               tryrac = .TRUE.
            ELSE
               tryrac = .FALSE.
            ENDIF
            CALL DSTEMR(Jobz,'A',N,Work(inddd),Work(indee),Vl,Vu,Il,Iu, &
     &                  M,W,Z,Ldz,N,Isuppz,tryrac,Work(indwk),Lwork,    &
     &                  Iwork,Liwork,Info)
!
!
!
!        Apply orthogonal matrix used in reduction to tridiagonal
!        form to eigenvectors returned by DSTEMR.
!
            IF ( wantz .AND. Info==0 ) THEN
               indwkn = inde
               llwrkn = Lwork - indwkn + 1
               CALL DORMTR('L',Uplo,'N',N,M,A,Lda,Work(indtau),Z,Ldz,   &
     &                     Work(indwkn),llwrkn,iinfo)
            ENDIF
         ENDIF
!
!
         IF ( Info==0 ) THEN
!           Everything worked.  Skip DSTEBZ/DSTEIN.  IWORK(:) are
!           undefined.
            M = N
            GOTO 100
         ENDIF
         Info = 0
      ENDIF
!
!     Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN.
!     Also call DSTEBZ and DSTEIN if DSTEMR fails.
!
      IF ( wantz ) THEN
         order = 'B'
      ELSE
         order = 'E'
      ENDIF
 
      CALL DSTEBZ(Range,order,N,vll,vuu,Il,Iu,abstll,Work(indd),        &
     &            Work(inde),M,nsplit,W,Iwork(indibl),Iwork(indisp),    &
     &            Work(indwk),Iwork(indiwo),Info)
!
      IF ( wantz ) THEN
         CALL DSTEIN(N,Work(indd),Work(inde),M,W,Iwork(indibl),         &
     &               Iwork(indisp),Z,Ldz,Work(indwk),Iwork(indiwo),     &
     &               Iwork(indifl),Info)
!
!        Apply orthogonal matrix used in reduction to tridiagonal
!        form to eigenvectors returned by DSTEIN.
!
         indwkn = inde
         llwrkn = Lwork - indwkn + 1
         CALL DORMTR('L',Uplo,'N',N,M,A,Lda,Work(indtau),Z,Ldz,         &
     &               Work(indwkn),llwrkn,iinfo)
      ENDIF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
!  Jump here if DSTEMR/DSTEIN succeeded.
 100  IF ( iscale==1 ) THEN
         IF ( Info==0 ) THEN
            imax = M
         ELSE
            imax = Info - 1
         ENDIF
         CALL DSCAL(imax,ONE/sigma,W,1)
      ENDIF
!
!     If eigenvalues are not in order, then sort them, along with
!     eigenvectors.  Note: We do not sort the IFAIL portion of IWORK.
!     It may not be initialized (if DSTEMR/DSTEIN succeeded), and we do
!     not return this detailed information to the user.
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
               W(i) = W(j)
               W(j) = tmp1
               CALL DSWAP(N,Z(1,i),1,Z(1,j),1)
            ENDIF
         ENDDO
      ENDIF
!
!     Set WORK(1) to optimal workspace size.
!
      Work(1) = lwkopt
      Iwork(1) = liwmin
!
!
!     End of DSYEVR
!
      END SUBROUTINE DSYEVR
