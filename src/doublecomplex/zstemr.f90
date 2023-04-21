!*==zstemr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZSTEMR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZSTEMR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zstemr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zstemr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zstemr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
!                          M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC, WORK, LWORK,
!                          IWORK, LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, RANGE
!       LOGICAL            TRYRAC
!       INTEGER            IL, INFO, IU, LDZ, NZC, LIWORK, LWORK, M, N
!       DOUBLE PRECISION VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            ISUPPZ( * ), IWORK( * )
!       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
!       COMPLEX*16         Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSTEMR computes selected eigenvalues and, optionally, eigenvectors
!> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
!> a well defined set of pairwise different real eigenvalues, the corresponding
!> real eigenvectors are pairwise orthogonal.
!>
!> The spectrum may be computed either completely or partially by specifying
!> either an interval (VL,VU] or a range of indices IL:IU for the desired
!> eigenvalues.
!>
!> Depending on the number of desired eigenvalues, these are computed either
!> by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are
!> computed by the use of various suitable L D L^T factorizations near clusters
!> of close eigenvalues (referred to as RRRs, Relatively Robust
!> Representations). An informal sketch of the algorithm follows.
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
!> For more details, see:
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
!> Further Details
!> 1.ZSTEMR works only on machines which follow IEEE-754
!> floating-point standard in their handling of infinities and NaNs.
!> This permits the use of efficient inner loops avoiding a check for
!> zero divisors.
!>
!> 2. LAPACK routines can be used to reduce a complex Hermitean matrix to
!> real symmetric tridiagonal form.
!>
!> (Any complex Hermitean tridiagonal matrix has real values on its diagonal
!> and potentially complex numbers on its off-diagonals. By applying a
!> similarity transform with an appropriate diagonal matrix
!> diag(1,e^{i \phy_1}, ... , e^{i \phy_{n-1}}), the complex Hermitean
!> matrix can be transformed into a real symmetric matrix and complex
!> arithmetic can be entirely avoided.)
!>
!> While the eigenvectors of the real symmetric tridiagonal matrix are real,
!> the eigenvectors of original complex Hermitean matrix have complex entries
!> in general.
!> Since LAPACK drivers overwrite the matrix data with the eigenvectors,
!> ZSTEMR accepts complex workspace to facilitate interoperability
!> with ZUNMTR or ZUPMTR.
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
!>          On entry, the N diagonal elements of the tridiagonal matrix
!>          T. On exit, D is overwritten.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N)
!>          On entry, the (N-1) subdiagonal elements of the tridiagonal
!>          matrix T in elements 1 to N-1 of E. E(N) need not be set on
!>          input, but is used internally as workspace.
!>          On exit, E is overwritten.
!> \endverbatim
!>
!> \param[in] VL
!> \verbatim
!>          VL is DOUBLE PRECISION
!>
!>          If RANGE='V', the lower bound of the interval to
!>          be searched for eigenvalues. VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>          VU is DOUBLE PRECISION
!>
!>          If RANGE='V', the upper bound of the interval to
!>          be searched for eigenvalues. VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] IL
!> \verbatim
!>          IL is INTEGER
!>
!>          If RANGE='I', the index of the
!>          smallest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] IU
!> \verbatim
!>          IU is INTEGER
!>
!>          If RANGE='I', the index of the
!>          largest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0.
!>          Not referenced if RANGE = 'A' or 'V'.
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
!>          Z is COMPLEX*16 array, dimension (LDZ, max(1,M) )
!>          If JOBZ = 'V', and if INFO = 0, then the first M columns of Z
!>          contain the orthonormal eigenvectors of the matrix T
!>          corresponding to the selected eigenvalues, with the i-th
!>          column of Z holding the eigenvector associated with W(i).
!>          If JOBZ = 'N', then Z is not referenced.
!>          Note: the user must ensure that at least max(1,M) columns are
!>          supplied in the array Z; if RANGE = 'V', the exact value of M
!>          is not known in advance and can be computed with a workspace
!>          query by setting NZC = -1, see below.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          JOBZ = 'V', then LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[in] NZC
!> \verbatim
!>          NZC is INTEGER
!>          The number of eigenvectors to be held in the array Z.
!>          If RANGE = 'A', then NZC >= max(1,N).
!>          If RANGE = 'V', then NZC >= the number of eigenvalues in (VL,VU].
!>          If RANGE = 'I', then NZC >= IU-IL+1.
!>          If NZC = -1, then a workspace query is assumed; the
!>          routine calculates the number of columns of the array Z that
!>          are needed to hold the eigenvectors.
!>          This value is returned as the first entry of the Z array, and
!>          no error message related to NZC is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] ISUPPZ
!> \verbatim
!>          ISUPPZ is INTEGER array, dimension ( 2*max(1,M) )
!>          The support of the eigenvectors in Z, i.e., the indices
!>          indicating the nonzero elements in Z. The i-th computed eigenvector
!>          is nonzero only in elements ISUPPZ( 2*i-1 ) through
!>          ISUPPZ( 2*i ). This is relevant in the case when the matrix
!>          is split. ISUPPZ is only accessed when JOBZ is 'V' and N > 0.
!> \endverbatim
!>
!> \param[in,out] TRYRAC
!> \verbatim
!>          TRYRAC is LOGICAL
!>          If TRYRAC = .TRUE., indicates that the code should check whether
!>          the tridiagonal matrix defines its eigenvalues to high relative
!>          accuracy.  If so, the code uses relative-accuracy preserving
!>          algorithms that might be (a bit) slower depending on the matrix.
!>          If the matrix does not define its eigenvalues to high relative
!>          accuracy, the code can uses possibly faster algorithms.
!>          If TRYRAC = .FALSE., the code is not required to guarantee
!>          relatively accurate eigenvalues and can use the fastest possible
!>          techniques.
!>          On exit, a .TRUE. TRYRAC will be set to .FALSE. if the matrix
!>          does not define its eigenvalues to high relative accuracy.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!>          On exit, if INFO = 0, WORK(1) returns the optimal
!>          (and minimal) LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= max(1,18*N)
!>          if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'.
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (LIWORK)
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.  LIWORK >= max(1,10*N)
!>          if the eigenvectors are desired, and LIWORK >= max(1,8*N)
!>          if only the eigenvalues are to be computed.
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of the IWORK array,
!>          returns this value as the first entry of the IWORK array, and
!>          no error message related to LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          On exit, INFO
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = 1X, internal error in DLARRE,
!>                if INFO = 2X, internal error in ZLARRV.
!>                Here, the digit X = ABS( IINFO ) < 10, where IINFO is
!>                the nonzero error code returned by DLARRE or
!>                ZLARRV, respectively.
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
!> \ingroup complex16OTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Beresford Parlett, University of California, Berkeley, USA \n
!> Jim Demmel, University of California, Berkeley, USA \n
!> Inderjit Dhillon, University of Texas, Austin, USA \n
!> Osni Marques, LBNL/NERSC, USA \n
!> Christof Voemel, University of California, Berkeley, USA \n
!
!  =====================================================================
      SUBROUTINE ZSTEMR(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,M,W,Z,Ldz,Nzc,     &
     &                  Isuppz,Tryrac,Work,Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!*--ZSTEMR341
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobz , Range
      LOGICAL Tryrac
      INTEGER Il , Info , Iu , Ldz , Nzc , Liwork , Lwork , M , N
      DOUBLE PRECISION Vl , Vu
!     ..
!     .. Array Arguments ..
      INTEGER Isuppz(*) , Iwork(*)
      DOUBLE PRECISION D(*) , E(*) , W(*) , Work(*)
      COMPLEX*16 Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , FOUR , MINRGP
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,FOUR=4.0D0,MINRGP=1.0D-3)
!     ..
!     .. Local Scalars ..
      LOGICAL alleig , indeig , lquery , valeig , wantz , zquery
      INTEGER i , ibegin , iend , ifirst , iil , iindbl , iindw ,       &
     &        iindwk , iinfo , iinspl , iiu , ilast , in , indd ,       &
     &        inde2 , inderr , indgp , indgrs , indwrk , itmp , itmp2 , &
     &        j , jblk , jj , liwmin , lwmin , nsplit , nzcmin ,        &
     &        offset , wbegin , wend
      DOUBLE PRECISION bignum , cs , eps , pivmin , r1 , r2 , rmax ,    &
     &                 rmin , rtol1 , rtol2 , safmin , scale , smlnum , &
     &                 sn , thresh , tmp , tnrm , wl , wu
!     ..
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , DLANST
      EXTERNAL LSAME , DLAMCH , DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DLAE2 , DLAEV2 , DLARRC , DLARRE , DLARRJ ,      &
     &         DLARRR , DLASRT , DSCAL , XERBLA , ZLARRV , ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , SQRT
 
 
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      wantz = LSAME(Jobz,'V')
      alleig = LSAME(Range,'A')
      valeig = LSAME(Range,'V')
      indeig = LSAME(Range,'I')
!
      lquery = ((Lwork==-1) .OR. (Liwork==-1))
      zquery = (Nzc==-1)
 
!     DSTEMR needs WORK of size 6*N, IWORK of size 3*N.
!     In addition, DLARRE needs WORK of size 6*N, IWORK of size 5*N.
!     Furthermore, ZLARRV needs WORK of size 12*N, IWORK of size 7*N.
      IF ( wantz ) THEN
         lwmin = 18*N
         liwmin = 10*N
      ELSE
!        need less workspace if only the eigenvalues are wanted
         lwmin = 12*N
         liwmin = 8*N
      ENDIF
 
      wl = ZERO
      wu = ZERO
      iil = 0
      iiu = 0
      nsplit = 0
 
      IF ( valeig ) THEN
!        We do not reference VL, VU in the cases RANGE = 'I','A'
!        The interval (WL, WU] contains all the wanted eigenvalues.
!        It is either given by the user or computed in DLARRE.
         wl = Vl
         wu = Vu
      ELSEIF ( indeig ) THEN
!        We do not reference IL, IU in the cases RANGE = 'V','A'
         iil = Il
         iiu = Iu
      ENDIF
!
      Info = 0
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(alleig .OR. valeig .OR. indeig) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( valeig .AND. N>0 .AND. wu<=wl ) THEN
         Info = -7
      ELSEIF ( indeig .AND. (iil<1 .OR. iil>N) ) THEN
         Info = -8
      ELSEIF ( indeig .AND. (iiu<iil .OR. iiu>N) ) THEN
         Info = -9
      ELSEIF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) THEN
         Info = -13
      ELSEIF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
         Info = -17
      ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
         Info = -19
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
      IF ( Info==0 ) THEN
         Work(1) = lwmin
         Iwork(1) = liwmin
!
         IF ( wantz .AND. alleig ) THEN
            nzcmin = N
         ELSEIF ( wantz .AND. valeig ) THEN
            CALL DLARRC('T',N,Vl,Vu,D,E,safmin,nzcmin,itmp,itmp2,Info)
         ELSEIF ( wantz .AND. indeig ) THEN
            nzcmin = iiu - iil + 1
         ELSE
!           WANTZ .EQ. FALSE.
            nzcmin = 0
         ENDIF
         IF ( zquery .AND. Info==0 ) THEN
            Z(1,1) = nzcmin
         ELSEIF ( Nzc<nzcmin .AND. .NOT.zquery ) THEN
            Info = -14
         ENDIF
      ENDIF
 
      IF ( Info/=0 ) THEN
!
         CALL XERBLA('ZSTEMR',-Info)
!
         RETURN
      ELSEIF ( lquery .OR. zquery ) THEN
         RETURN
      ENDIF
!
!     Handle N = 0, 1, and 2 cases immediately
!
      M = 0
      IF ( N==0 ) RETURN
!
      IF ( N==1 ) THEN
         IF ( alleig .OR. indeig ) THEN
            M = 1
            W(1) = D(1)
         ELSEIF ( wl<D(1) .AND. wu>=D(1) ) THEN
            M = 1
            W(1) = D(1)
         ENDIF
         IF ( wantz .AND. (.NOT.zquery) ) THEN
            Z(1,1) = ONE
            Isuppz(1) = 1
            Isuppz(2) = 1
         ENDIF
         RETURN
      ENDIF
!
      IF ( N==2 ) THEN
         IF ( .NOT.wantz ) THEN
            CALL DLAE2(D(1),E(1),D(2),r1,r2)
         ELSEIF ( wantz .AND. (.NOT.zquery) ) THEN
            CALL DLAEV2(D(1),E(1),D(2),r1,r2,cs,sn)
         ENDIF
         IF ( alleig .OR. (valeig .AND. (r2>wl) .AND. (r2<=wu)) .OR.    &
     &        (indeig .AND. (iil==1)) ) THEN
            M = M + 1
            W(M) = r2
            IF ( wantz .AND. (.NOT.zquery) ) THEN
               Z(1,M) = -sn
               Z(2,M) = cs
!              Note: At most one of SN and CS can be zero.
               IF ( sn==ZERO ) THEN
                  Isuppz(2*M-1) = 2
                  Isuppz(2*M) = 2
               ELSEIF ( cs/=ZERO ) THEN
                  Isuppz(2*M-1) = 1
                  Isuppz(2*M) = 2
               ELSE
                  Isuppz(2*M-1) = 1
                  Isuppz(2*M) = 1
               ENDIF
            ENDIF
         ENDIF
         IF ( alleig .OR. (valeig .AND. (r1>wl) .AND. (r1<=wu)) .OR.    &
     &        (indeig .AND. (iiu==2)) ) THEN
            M = M + 1
            W(M) = r1
            IF ( wantz .AND. (.NOT.zquery) ) THEN
               Z(1,M) = cs
               Z(2,M) = sn
!              Note: At most one of SN and CS can be zero.
               IF ( sn==ZERO ) THEN
                  Isuppz(2*M-1) = 2
                  Isuppz(2*M) = 2
               ELSEIF ( cs/=ZERO ) THEN
                  Isuppz(2*M-1) = 1
                  Isuppz(2*M) = 2
               ELSE
                  Isuppz(2*M-1) = 1
                  Isuppz(2*M) = 1
               ENDIF
            ENDIF
         ENDIF
      ELSE
 
!        Continue with general N
 
         indgrs = 1
         inderr = 2*N + 1
         indgp = 3*N + 1
         indd = 4*N + 1
         inde2 = 5*N + 1
         indwrk = 6*N + 1
!
         iinspl = 1
         iindbl = N + 1
         iindw = 2*N + 1
         iindwk = 3*N + 1
!
!        Scale matrix to allowable range, if necessary.
!        The allowable range is related to the PIVMIN parameter; see the
!        comments in DLARRD.  The preference for scaling small values
!        up is heuristic; we expect users' matrices not to be close to the
!        RMAX threshold.
!
         scale = ONE
         tnrm = DLANST('M',N,D,E)
         IF ( tnrm>ZERO .AND. tnrm<rmin ) THEN
            scale = rmin/tnrm
         ELSEIF ( tnrm>rmax ) THEN
            scale = rmax/tnrm
         ENDIF
         IF ( scale/=ONE ) THEN
            CALL DSCAL(N,scale,D,1)
            CALL DSCAL(N-1,scale,E,1)
            tnrm = tnrm*scale
            IF ( valeig ) THEN
!              If eigenvalues in interval have to be found,
!              scale (WL, WU] accordingly
               wl = wl*scale
               wu = wu*scale
            ENDIF
         ENDIF
!
!        Compute the desired eigenvalues of the tridiagonal after splitting
!        into smaller subblocks if the corresponding off-diagonal elements
!        are small
!        THRESH is the splitting parameter for DLARRE
!        A negative THRESH forces the old splitting criterion based on the
!        size of the off-diagonal. A positive THRESH switches to splitting
!        which preserves relative accuracy.
!
         IF ( Tryrac ) THEN
!           Test whether the matrix warrants the more expensive relative approach.
            CALL DLARRR(N,D,E,iinfo)
         ELSE
!           The user does not care about relative accurately eigenvalues
            iinfo = -1
         ENDIF
!        Set the splitting criterion
         IF ( iinfo==0 ) THEN
            thresh = eps
         ELSE
            thresh = -eps
!           relative accuracy is desired but T does not guarantee it
            Tryrac = .FALSE.
         ENDIF
!
!           Copy original diagonal, needed to guarantee relative accuracy
         IF ( Tryrac ) CALL DCOPY(N,D,1,Work(indd),1)
!        Store the squares of the offdiagonal values of T
         DO j = 1 , N - 1
            Work(inde2+j-1) = E(j)**2
         ENDDO
 
!        Set the tolerance parameters for bisection
         IF ( .NOT.wantz ) THEN
!           DLARRE computes the eigenvalues to full precision.
            rtol1 = FOUR*eps
            rtol2 = FOUR*eps
         ELSE
!           DLARRE computes the eigenvalues to less than full precision.
!           ZLARRV will refine the eigenvalue approximations, and we only
!           need less accurate initial bisection in DLARRE.
!           Note: these settings do only affect the subset case and DLARRE
            rtol1 = SQRT(eps)
            rtol2 = MAX(SQRT(eps)*5.0D-3,FOUR*eps)
         ENDIF
         CALL DLARRE(Range,N,wl,wu,iil,iiu,D,E,Work(inde2),rtol1,rtol2, &
     &               thresh,nsplit,Iwork(iinspl),M,W,Work(inderr),      &
     &               Work(indgp),Iwork(iindbl),Iwork(iindw),Work(indgrs)&
     &               ,pivmin,Work(indwrk),Iwork(iindwk),iinfo)
         IF ( iinfo/=0 ) THEN
            Info = 10 + ABS(iinfo)
            RETURN
         ENDIF
!        Note that if RANGE .NE. 'V', DLARRE computes bounds on the desired
!        part of the spectrum. All desired eigenvalues are contained in
!        (WL,WU]
 
 
         IF ( wantz ) THEN
!
!           Compute the desired eigenvectors corresponding to the computed
!           eigenvalues
!
            CALL ZLARRV(N,wl,wu,D,E,pivmin,Iwork(iinspl),M,1,M,MINRGP,  &
     &                  rtol1,rtol2,W,Work(inderr),Work(indgp),         &
     &                  Iwork(iindbl),Iwork(iindw),Work(indgrs),Z,Ldz,  &
     &                  Isuppz,Work(indwrk),Iwork(iindwk),iinfo)
            IF ( iinfo/=0 ) THEN
               Info = 20 + ABS(iinfo)
               RETURN
            ENDIF
         ELSE
!           DLARRE computes eigenvalues of the (shifted) root representation
!           ZLARRV returns the eigenvalues of the unshifted matrix.
!           However, if the eigenvectors are not desired by the user, we need
!           to apply the corresponding shifts from DLARRE to obtain the
!           eigenvalues of the original matrix.
            DO j = 1 , M
               itmp = Iwork(iindbl+j-1)
               W(j) = W(j) + E(Iwork(iinspl+itmp-1))
            ENDDO
         ENDIF
!
 
         IF ( Tryrac ) THEN
!           Refine computed eigenvalues so that they are relatively accurate
!           with respect to the original matrix T.
            ibegin = 1
            wbegin = 1
            DO jblk = 1 , Iwork(iindbl+M-1)
               iend = Iwork(iinspl+jblk-1)
               in = iend - ibegin + 1
               wend = wbegin - 1
               DO
!              check if any eigenvalues have to be refined in this block
                  IF ( wend<M ) THEN
                     IF ( Iwork(iindbl+wend)==jblk ) THEN
                        wend = wend + 1
                        CYCLE
                     ENDIF
                  ENDIF
                  IF ( wend<wbegin ) THEN
                     ibegin = iend + 1
                     EXIT
                  ENDIF
 
                  offset = Iwork(iindw+wbegin-1) - 1
                  ifirst = Iwork(iindw+wbegin-1)
                  ilast = Iwork(iindw+wend-1)
                  rtol2 = FOUR*eps
                  CALL DLARRJ(in,Work(indd+ibegin-1),                   &
     &                        Work(inde2+ibegin-1),ifirst,ilast,rtol2,  &
     &                        offset,W(wbegin),Work(inderr+wbegin-1),   &
     &                        Work(indwrk),Iwork(iindwk),pivmin,tnrm,   &
     &                        iinfo)
                  ibegin = iend + 1
                  wbegin = wend + 1
                  EXIT
               ENDDO
            ENDDO
         ENDIF
!
!        If matrix was scaled, then rescale eigenvalues appropriately.
!
         IF ( scale/=ONE ) CALL DSCAL(M,ONE/scale,W,1)
      ENDIF
!
!     If eigenvalues are not in increasing order, then sort them,
!     possibly along with eigenvectors.
!
      IF ( nsplit>1 .OR. N==2 ) THEN
         IF ( .NOT.wantz ) THEN
            CALL DLASRT('I',M,W,iinfo)
            IF ( iinfo/=0 ) THEN
               Info = 3
               RETURN
            ENDIF
         ELSE
            DO j = 1 , M - 1
               i = 0
               tmp = W(j)
               DO jj = j + 1 , M
                  IF ( W(jj)<tmp ) THEN
                     i = jj
                     tmp = W(jj)
                  ENDIF
               ENDDO
               IF ( i/=0 ) THEN
                  W(i) = W(j)
                  W(j) = tmp
                  IF ( wantz ) THEN
                     CALL ZSWAP(N,Z(1,i),1,Z(1,j),1)
                     itmp = Isuppz(2*i-1)
                     Isuppz(2*i-1) = Isuppz(2*j-1)
                     Isuppz(2*j-1) = itmp
                     itmp = Isuppz(2*i)
                     Isuppz(2*i) = Isuppz(2*j)
                     Isuppz(2*j) = itmp
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDIF
!
!
      Work(1) = lwmin
      Iwork(1) = liwmin
!
!     End of ZSTEMR
!
      END SUBROUTINE ZSTEMR
