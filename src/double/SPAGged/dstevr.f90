!*==dstevr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> DSTEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSTEVR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstevr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstevr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstevr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSTEVR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL,
!                          M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
!                          LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, RANGE
!       INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, M, N
!       DOUBLE PRECISION   ABSTOL, VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            ISUPPZ( * ), IWORK( * )
!       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSTEVR computes selected eigenvalues and, optionally, eigenvectors
!> of a real symmetric tridiagonal matrix T.  Eigenvalues and
!> eigenvectors can be selected by specifying either a range of values
!> or a range of indices for the desired eigenvalues.
!>
!> Whenever possible, DSTEVR calls DSTEMR to compute the
!> eigenspectrum using Relatively Robust Representations.  DSTEMR
!> computes eigenvalues by the dqds algorithm, while orthogonal
!> eigenvectors are computed from various "good" L D L^T representations
!> (also known as Relatively Robust Representations). Gram-Schmidt
!> orthogonalization is avoided as far as possible. More specifically,
!> the various steps of the algorithm are as follows. For the i-th
!> unreduced block of T,
!>    (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T
!>         is a relatively robust representation,
!>    (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high
!>        relative accuracy by the dqds algorithm,
!>    (c) If there is a cluster of close eigenvalues, "choose" sigma_i
!>        close to the cluster, and go to step (a),
!>    (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,
!>        compute the corresponding eigenvector by forming a
!>        rank-revealing twisted factorization.
!> The desired accuracy of the output can be specified by the input
!> parameter ABSTOL.
!>
!> For more details, see "A new O(n^2) algorithm for the symmetric
!> tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,
!> Computer Science Division Technical Report No. UCB//CSD-97-971,
!> UC Berkeley, May 1997.
!>
!>
!> Note 1 : DSTEVR calls DSTEMR when the full spectrum is requested
!> on machines which conform to the ieee-754 floating point standard.
!> DSTEVR calls DSTEBZ and DSTEIN on non-ieee machines and
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
!>          On exit, D may be multiplied by a constant factor chosen
!>          to avoid over/underflow in computing the eigenvalues.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (max(1,N-1))
!>          On entry, the (n-1) subdiagonal elements of the tridiagonal
!>          matrix A in elements 1 to N-1 of E.
!>          On exit, E may be multiplied by a constant factor chosen
!>          to avoid over/underflow in computing the eigenvalues.
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
!>          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M) )
!>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!>          contain the orthonormal eigenvectors of the matrix A
!>          corresponding to the selected eigenvalues, with the i-th
!>          column of Z holding the eigenvector associated with W(i).
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
!>          ISUPPZ( 2*i ).
!>          Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal (and
!>          minimal) LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,20*N).
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
!>          On exit, if INFO = 0, IWORK(1) returns the optimal (and
!>          minimal) LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.  LIWORK >= max(1,10*N).
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
!> \ingroup doubleOTHEReigen
!
!> \par Contributors:
!  ==================
!>
!>     Inderjit Dhillon, IBM Almaden, USA \n
!>     Osni Marques, LBNL/NERSC, USA \n
!>     Ken Stanley, Computer Science Division, University of
!>       California at Berkeley, USA \n
!>
!  =====================================================================
      SUBROUTINE DSTEVR(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,  &
     &                  Isuppz,Work,Lwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DLAMCH
      USE S_DLANST
      USE S_DSCAL
      USE S_DSTEBZ
      USE S_DSTEIN
      USE S_DSTEMR
      USE S_DSTERF
      USE S_DSWAP
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DSTEVR320
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
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: alleig , indeig , lquery , test , tryrac , valeig ,    &
     &           wantz
      REAL(R8KIND) :: bignum , eps , rmax , rmin , safmin , sigma ,     &
     &                smlnum , tmp1 , tnrm , vll , vuu
      INTEGER :: i , ieeeok , imax , indibl , indifl , indisp , indiwo ,&
     &           iscale , itmp1 , j , jj , liwmin , lwmin , nsplit
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
!
!     Test the input parameters.
!
      ieeeok = ILAENV(10,'DSTEVR','N',1,2,3,4)
!
      wantz = LSAME(Jobz,'V')
      alleig = LSAME(Range,'A')
      valeig = LSAME(Range,'V')
      indeig = LSAME(Range,'I')
!
      lquery = ((Lwork==-1) .OR. (Liwork==-1))
      lwmin = MAX(1,20*N)
      liwmin = MAX(1,10*N)
!
!
      Info = 0
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(alleig .OR. valeig .OR. indeig) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( valeig ) THEN
         IF ( N>0 .AND. Vu<=Vl ) Info = -7
      ELSEIF ( indeig ) THEN
         IF ( Il<1 .OR. Il>MAX(1,N) ) THEN
            Info = -8
         ELSEIF ( Iu<MIN(N,Il) .OR. Iu>N ) THEN
            Info = -9
         ENDIF
      ENDIF
      IF ( Info==0 ) THEN
         IF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) Info = -14
      ENDIF
!
      IF ( Info==0 ) THEN
         Work(1) = lwmin
         Iwork(1) = liwmin
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -17
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -19
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSTEVR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      M = 0
      IF ( N==0 ) RETURN
!
      IF ( N==1 ) THEN
         IF ( alleig .OR. indeig ) THEN
            M = 1
            W(1) = D(1)
         ELSEIF ( Vl<D(1) .AND. Vu>=D(1) ) THEN
            M = 1
            W(1) = D(1)
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
      rmax = MIN(SQRT(bignum),ONE/SQRT(SQRT(safmin)))
!
!
!     Scale matrix to allowable range, if necessary.
!
      iscale = 0
      IF ( valeig ) THEN
         vll = Vl
         vuu = Vu
      ENDIF
!
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
         IF ( valeig ) THEN
            vll = Vl*sigma
            vuu = Vu*sigma
         ENDIF
      ENDIF
 
!     Initialize indices into workspaces.  Note: These indices are used only
!     if DSTERF or DSTEMR fail.
 
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
      indiwo = indisp + N
!
!     If all eigenvalues are desired, then
!     call DSTERF or DSTEMR.  If this fails for some eigenvalue, then
!     try DSTEBZ.
!
!
      test = .FALSE.
      IF ( indeig ) THEN
         IF ( Il==1 .AND. Iu==N ) test = .TRUE.
      ENDIF
      IF ( (alleig .OR. test) .AND. ieeeok==1 ) THEN
         CALL DCOPY(N-1,E(1),1,Work(1),1)
         IF ( .NOT.wantz ) THEN
            CALL DCOPY(N,D,1,W,1)
            CALL DSTERF(N,W,Work,Info)
         ELSE
            CALL DCOPY(N,D,1,Work(N+1),1)
            IF ( Abstol<=TWO*N*eps ) THEN
               tryrac = .TRUE.
            ELSE
               tryrac = .FALSE.
            ENDIF
            CALL DSTEMR(Jobz,'A',N,Work(N+1),Work,Vl,Vu,Il,Iu,M,W,Z,Ldz,&
     &                  N,Isuppz,tryrac,Work(2*N+1),Lwork-2*N,Iwork,    &
     &                  Liwork,Info)
!
         ENDIF
         IF ( Info==0 ) THEN
            M = N
            GOTO 100
         ENDIF
         Info = 0
      ENDIF
!
!     Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN.
!
      IF ( wantz ) THEN
         order = 'B'
      ELSE
         order = 'E'
      ENDIF
 
      CALL DSTEBZ(Range,order,N,vll,vuu,Il,Iu,Abstol,D,E,M,nsplit,W,    &
     &            Iwork(indibl),Iwork(indisp),Work,Iwork(indiwo),Info)
!
      IF ( wantz ) CALL DSTEIN(N,D,E,M,W,Iwork(indibl),Iwork(indisp),Z, &
     &                         Ldz,Work,Iwork(indiwo),Iwork(indifl),    &
     &                         Info)
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
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
               itmp1 = Iwork(i)
               W(i) = W(j)
               Iwork(i) = Iwork(j)
               W(j) = tmp1
               Iwork(j) = itmp1
               CALL DSWAP(N,Z(1,i),1,Z(1,j),1)
            ENDIF
         ENDDO
      ENDIF
!
!      Causes problems with tests 19 & 20:
!      IF (wantz .and. INDEIG ) Z( 1,1) = Z(1,1) / 1.002 + .002
!
!
      Work(1) = lwmin
      Iwork(1) = liwmin
!
!     End of DSTEVR
!
      END SUBROUTINE DSTEVR
