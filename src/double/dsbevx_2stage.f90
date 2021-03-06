!*==dsbevx_2stage.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> DSBEVX_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  @precisions fortran d -> s
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSBEVX_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbevx_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbevx_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbevx_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSBEVX_2STAGE( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q,
!                                 LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z,
!                                 LDZ, WORK, LWORK, IWORK, IFAIL, INFO )
!
!       IMPLICIT NONE
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, RANGE, UPLO
!       INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N, LWORK
!       DOUBLE PRECISION   ABSTOL, VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            IFAIL( * ), IWORK( * )
!       DOUBLE PRECISION   AB( LDAB, * ), Q( LDQ, * ), W( * ), WORK( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSBEVX_2STAGE computes selected eigenvalues and, optionally, eigenvectors
!> of a real symmetric band matrix A using the 2stage technique for
!> the reduction to tridiagonal. Eigenvalues and eigenvectors can
!> be selected by specifying either a range of values or a range of
!> indices for the desired eigenvalues.
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
!>          = 'A': all eigenvalues will be found;
!>          = 'V': all eigenvalues in the half-open interval (VL,VU]
!>                 will be found;
!>          = 'I': the IL-th through IU-th eigenvalues will be found.
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
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
!>          If JOBZ = 'V', the N-by-N orthogonal matrix used in the
!>                         reduction to tridiagonal form.
!>          If JOBZ = 'N', the array Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  If JOBZ = 'V', then
!>          LDQ >= max(1,N).
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
!>          by reducing AB to tridiagonal form.
!>
!>          Eigenvalues will be computed most accurately when ABSTOL is
!>          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!>          If this routine returns with INFO>0, indicating that some
!>          eigenvectors did not converge, try setting ABSTOL to
!>          2*DLAMCH('S').
!>
!>          See "Computing Small Singular Values of Bidiagonal Matrices
!>          with Guaranteed High Relative Accuracy," by Demmel and
!>          Kahan, LAPACK Working Note #3.
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
!>          If an eigenvector fails to converge, then that column of Z
!>          contains the latest approximation to the eigenvector, and the
!>          index of the eigenvector is returned in IFAIL.
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
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK. LWORK >= 1, when N <= 1;
!>          otherwise
!>          If JOBZ = 'N' and N > 1, LWORK must be queried.
!>                                   LWORK = MAX(1, 7*N, dimension) where
!>                                   dimension = (2KD+1)*N + KD*NTHREADS + 2*N
!>                                   where KD is the size of the band.
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
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (5*N)
!> \endverbatim
!>
!> \param[out] IFAIL
!> \verbatim
!>          IFAIL is INTEGER array, dimension (N)
!>          If JOBZ = 'V', then if INFO = 0, the first M elements of
!>          IFAIL are zero.  If INFO > 0, then IFAIL contains the
!>          indices of the eigenvectors that failed to converge.
!>          If JOBZ = 'N', then IFAIL is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, then i eigenvectors failed to converge.
!>                Their indices are stored in array IFAIL.
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
      SUBROUTINE DSBEVX_2STAGE(Jobz,Range,Uplo,N,Kd,Ab,Ldab,Q,Ldq,Vl,Vu,&
     &                         Il,Iu,Abstol,M,W,Z,Ldz,Work,Lwork,Iwork, &
     &                         Ifail,Info)
!
      IMPLICIT NONE
!*--DSBEVX_2STAGE327
!
!  -- LAPACK driver routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobz , Range , Uplo
      INTEGER Il , Info , Iu , Kd , Ldab , Ldq , Ldz , M , N , Lwork
      DOUBLE PRECISION Abstol , Vl , Vu
!     ..
!     .. Array Arguments ..
      INTEGER Ifail(*) , Iwork(*)
      DOUBLE PRECISION Ab(Ldab,*) , Q(Ldq,*) , W(*) , Work(*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      LOGICAL alleig , indeig , lower , test , valeig , wantz , lquery
      CHARACTER order
      INTEGER i , iinfo , imax , indd , inde , indee , indibl , indisp ,&
     &        indiwo , indwrk , iscale , itmp1 , j , jj , llwork ,      &
     &        lwmin , lhtrd , lwtrd , ib , indhous , nsplit
      DOUBLE PRECISION abstll , anrm , bignum , eps , rmax , rmin ,     &
     &                 safmin , sigma , smlnum , tmp1 , vll , vuu
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV2STAGE
      DOUBLE PRECISION DLAMCH , DLANSB
      EXTERNAL LSAME , DLAMCH , DLANSB , ILAENV2STAGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DGEMV , DLACPY , DLASCL , DSCAL , DSTEBZ ,       &
     &         DSTEIN , DSTEQR , DSTERF , DSWAP , XERBLA , DSYTRD_SB2ST
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
      lower = LSAME(Uplo,'L')
      lquery = (Lwork==-1)
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
      ELSEIF ( Kd<0 ) THEN
         Info = -5
      ELSEIF ( Ldab<Kd+1 ) THEN
         Info = -7
      ELSEIF ( wantz .AND. Ldq<MAX(1,N) ) THEN
         Info = -9
      ELSEIF ( valeig ) THEN
         IF ( N>0 .AND. Vu<=Vl ) Info = -11
      ELSEIF ( indeig ) THEN
         IF ( Il<1 .OR. Il>MAX(1,N) ) THEN
            Info = -12
         ELSEIF ( Iu<MIN(N,Il) .OR. Iu>N ) THEN
            Info = -13
         ENDIF
      ENDIF
      IF ( Info==0 ) THEN
         IF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) Info = -18
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( N<=1 ) THEN
            lwmin = 1
            Work(1) = lwmin
         ELSE
            ib = ILAENV2STAGE(2,'DSYTRD_SB2ST',Jobz,N,Kd,-1,-1)
            lhtrd = ILAENV2STAGE(3,'DSYTRD_SB2ST',Jobz,N,Kd,ib,-1)
            lwtrd = ILAENV2STAGE(4,'DSYTRD_SB2ST',Jobz,N,Kd,ib,-1)
            lwmin = 2*N + lhtrd + lwtrd
            Work(1) = lwmin
         ENDIF
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) Info = -20
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSBEVX_2STAGE ',-Info)
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
         M = 1
         IF ( lower ) THEN
            tmp1 = Ab(1,1)
         ELSE
            tmp1 = Ab(Kd+1,1)
         ENDIF
         IF ( valeig ) THEN
            IF ( Vl>=tmp1 .OR. Vu<tmp1 ) M = 0
         ENDIF
         IF ( M==1 ) THEN
            W(1) = tmp1
            IF ( wantz ) Z(1,1) = ONE
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
      ELSE
         vll = ZERO
         vuu = ZERO
      ENDIF
      anrm = DLANSB('M',Uplo,N,Kd,Ab,Ldab,Work)
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
         IF ( Abstol>0 ) abstll = Abstol*sigma
         IF ( valeig ) THEN
            vll = Vl*sigma
            vuu = Vu*sigma
         ENDIF
      ENDIF
!
!     Call DSYTRD_SB2ST to reduce symmetric band matrix to tridiagonal form.
!
      indd = 1
      inde = indd + N
      indhous = inde + N
      indwrk = indhous + lhtrd
      llwork = Lwork - indwrk + 1
!
      CALL DSYTRD_SB2ST("N",Jobz,Uplo,N,Kd,Ab,Ldab,Work(indd),Work(inde)&
     &                  ,Work(indhous),lhtrd,Work(indwrk),llwork,iinfo)
!
!     If all eigenvalues are desired and ABSTOL is less than or equal
!     to zero, then call DSTERF or SSTEQR.  If this fails for some
!     eigenvalue, then try DSTEBZ.
!
      test = .FALSE.
      IF ( indeig ) THEN
         IF ( Il==1 .AND. Iu==N ) test = .TRUE.
      ENDIF
      IF ( (alleig .OR. test) .AND. (Abstol<=ZERO) ) THEN
         CALL DCOPY(N,Work(indd),1,W,1)
         indee = indwrk + 2*N
         IF ( .NOT.wantz ) THEN
            CALL DCOPY(N-1,Work(inde),1,Work(indee),1)
            CALL DSTERF(N,W,Work(indee),Info)
         ELSE
            CALL DLACPY('A',N,N,Q,Ldq,Z,Ldz)
            CALL DCOPY(N-1,Work(inde),1,Work(indee),1)
            CALL DSTEQR(Jobz,N,W,Work(indee),Z,Ldz,Work(indwrk),Info)
            IF ( Info==0 ) THEN
               DO i = 1 , N
                  Ifail(i) = 0
               ENDDO
            ENDIF
         ENDIF
         IF ( Info==0 ) THEN
            M = N
            GOTO 100
         ENDIF
         Info = 0
      ENDIF
!
!     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN.
!
      IF ( wantz ) THEN
         order = 'B'
      ELSE
         order = 'E'
      ENDIF
      indibl = 1
      indisp = indibl + N
      indiwo = indisp + N
      CALL DSTEBZ(Range,order,N,vll,vuu,Il,Iu,abstll,Work(indd),        &
     &            Work(inde),M,nsplit,W,Iwork(indibl),Iwork(indisp),    &
     &            Work(indwrk),Iwork(indiwo),Info)
!
      IF ( wantz ) THEN
         CALL DSTEIN(N,Work(indd),Work(inde),M,W,Iwork(indibl),         &
     &               Iwork(indisp),Z,Ldz,Work(indwrk),Iwork(indiwo),    &
     &               Ifail,Info)
!
!        Apply orthogonal matrix used in reduction to tridiagonal
!        form to eigenvectors returned by DSTEIN.
!
         DO j = 1 , M
            CALL DCOPY(N,Z(1,j),1,Work(1),1)
            CALL DGEMV('N',N,N,ONE,Q,Ldq,Work,1,ZERO,Z(1,j),1)
         ENDDO
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
               itmp1 = Iwork(indibl+i-1)
               W(i) = W(j)
               Iwork(indibl+i-1) = Iwork(indibl+j-1)
               W(j) = tmp1
               Iwork(indibl+j-1) = itmp1
               CALL DSWAP(N,Z(1,i),1,Z(1,j),1)
               IF ( Info/=0 ) THEN
                  itmp1 = Ifail(i)
                  Ifail(i) = Ifail(j)
                  Ifail(j) = itmp1
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
!     Set WORK(1) to optimal workspace size.
!
      Work(1) = lwmin
!
!
!     End of DSBEVX_2STAGE
!
      END SUBROUTINE DSBEVX_2STAGE
