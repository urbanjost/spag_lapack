!*==chpevx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> CHPEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHPEVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpevx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpevx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpevx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,
!                          ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK,
!                          IFAIL, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, RANGE, UPLO
!       INTEGER            IL, INFO, IU, LDZ, M, N
!       REAL               ABSTOL, VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            IFAIL( * ), IWORK( * )
!       REAL               RWORK( * ), W( * )
!       COMPLEX            AP( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHPEVX computes selected eigenvalues and, optionally, eigenvectors
!> of a complex Hermitian matrix A in packed storage.
!> Eigenvalues/vectors can be selected by specifying either a range of
!> values or a range of indices for the desired eigenvalues.
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
!> \param[in,out] AP
!> \verbatim
!>          AP is COMPLEX array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the Hermitian matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!>
!>          On exit, AP is overwritten by values generated during the
!>          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
!>          and first superdiagonal of the tridiagonal matrix T overwrite
!>          the corresponding elements of A, and if UPLO = 'L', the
!>          diagonal and first subdiagonal of T overwrite the
!>          corresponding elements of A.
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
!>          by reducing AP to tridiagonal form.
!>
!>          Eigenvalues will be computed most accurately when ABSTOL is
!>          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
!>          If this routine returns with INFO>0, indicating that some
!>          eigenvectors did not converge, try setting ABSTOL to
!>          2*SLAMCH('S').
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
!>          W is REAL array, dimension (N)
!>          If INFO = 0, the selected eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, max(1,M))
!>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!>          contain the orthonormal eigenvectors of the matrix A
!>          corresponding to the selected eigenvalues, with the i-th
!>          column of Z holding the eigenvector associated with W(i).
!>          If an eigenvector fails to converge, then that column of Z
!>          contains the latest approximation to the eigenvector, and
!>          the index of the eigenvector is returned in IFAIL.
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
!>          WORK is COMPLEX array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (7*N)
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
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complexOTHEReigen
!
!  =====================================================================
      SUBROUTINE CHPEVX(Jobz,Range,Uplo,N,Ap,Vl,Vu,Il,Iu,Abstol,M,W,Z,  &
     &                  Ldz,Work,Rwork,Iwork,Ifail,Info)
      USE S_CHPTRD
      USE S_CLANHP
      USE S_CSSCAL
      USE S_CSTEIN
      USE S_CSTEQR
      USE S_CSWAP
      USE S_CUPGTR
      USE S_CUPMTR
      USE S_LSAME
      USE S_SCOPY
      USE S_SLAMCH
      USE S_SSCAL
      USE S_SSTEBZ
      USE S_SSTERF
      USE S_XERBLA
      IMPLICIT NONE
!*--CHPEVX258
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CONE = (1.0E0,0.0E0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: abstll , anrm , bignum , eps , rmax , rmin , safmin ,     &
     &        sigma , smlnum , tmp1 , vll , vuu
      LOGICAL :: alleig , indeig , test , valeig , wantz
      INTEGER :: i , iinfo , imax , indd , inde , indee , indibl ,      &
     &           indisp , indiwk , indrwk , indtau , indwrk , iscale ,  &
     &           itmp1 , j , jj , nsplit
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
      wantz = LSAME(Jobz,'V')
      alleig = LSAME(Range,'A')
      valeig = LSAME(Range,'V')
      indeig = LSAME(Range,'I')
!
      Info = 0
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(alleig .OR. valeig .OR. indeig) ) THEN
         Info = -2
      ELSEIF ( .NOT.(LSAME(Uplo,'L') .OR. LSAME(Uplo,'U')) ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
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
      IF ( Info/=0 ) THEN
         CALL XERBLA('CHPEVX',-Info)
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
            W(1) = Ap(1)
         ELSEIF ( Vl<REAL(Ap(1)) .AND. Vu>=REAL(Ap(1)) ) THEN
            M = 1
            W(1) = Ap(1)
         ENDIF
         IF ( wantz ) Z(1,1) = CONE
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
      ELSE
         vll = ZERO
         vuu = ZERO
      ENDIF
      anrm = CLANHP('M',Uplo,N,Ap,Rwork)
      IF ( anrm>ZERO .AND. anrm<rmin ) THEN
         iscale = 1
         sigma = rmin/anrm
      ELSEIF ( anrm>rmax ) THEN
         iscale = 1
         sigma = rmax/anrm
      ENDIF
      IF ( iscale==1 ) THEN
         CALL CSSCAL((N*(N+1))/2,sigma,Ap,1)
         IF ( Abstol>0 ) abstll = Abstol*sigma
         IF ( valeig ) THEN
            vll = Vl*sigma
            vuu = Vu*sigma
         ENDIF
      ENDIF
!
!     Call CHPTRD to reduce Hermitian packed matrix to tridiagonal form.
!
      indd = 1
      inde = indd + N
      indrwk = inde + N
      indtau = 1
      indwrk = indtau + N
      CALL CHPTRD(Uplo,N,Ap,Rwork(indd),Rwork(inde),Work(indtau),iinfo)
!
!     If all eigenvalues are desired and ABSTOL is less than or equal
!     to zero, then call SSTERF or CUPGTR and CSTEQR.  If this fails
!     for some eigenvalue, then try SSTEBZ.
!
      test = .FALSE.
      IF ( indeig ) THEN
         IF ( Il==1 .AND. Iu==N ) test = .TRUE.
      ENDIF
      IF ( (alleig .OR. test) .AND. (Abstol<=ZERO) ) THEN
         CALL SCOPY(N,Rwork(indd),1,W,1)
         indee = indrwk + 2*N
         IF ( .NOT.wantz ) THEN
            CALL SCOPY(N-1,Rwork(inde),1,Rwork(indee),1)
            CALL SSTERF(N,W,Rwork(indee),Info)
         ELSE
            CALL CUPGTR(Uplo,N,Ap,Work(indtau),Z,Ldz,Work(indwrk),iinfo)
            CALL SCOPY(N-1,Rwork(inde),1,Rwork(indee),1)
            CALL CSTEQR(Jobz,N,W,Rwork(indee),Z,Ldz,Rwork(indrwk),Info)
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
!     Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN.
!
      IF ( wantz ) THEN
         order = 'B'
      ELSE
         order = 'E'
      ENDIF
      indibl = 1
      indisp = indibl + N
      indiwk = indisp + N
      CALL SSTEBZ(Range,order,N,vll,vuu,Il,Iu,abstll,Rwork(indd),       &
     &            Rwork(inde),M,nsplit,W,Iwork(indibl),Iwork(indisp),   &
     &            Rwork(indrwk),Iwork(indiwk),Info)
!
      IF ( wantz ) THEN
         CALL CSTEIN(N,Rwork(indd),Rwork(inde),M,W,Iwork(indibl),       &
     &               Iwork(indisp),Z,Ldz,Rwork(indrwk),Iwork(indiwk),   &
     &               Ifail,Info)
!
!        Apply unitary matrix used in reduction to tridiagonal
!        form to eigenvectors returned by CSTEIN.
!
         indwrk = indtau + N
         CALL CUPMTR('L',Uplo,'N',N,M,Ap,Work(indtau),Z,Ldz,Work(indwrk)&
     &               ,iinfo)
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
               IF ( Info/=0 ) THEN
                  itmp1 = Ifail(i)
                  Ifail(i) = Ifail(j)
                  Ifail(j) = itmp1
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of CHPEVX
!
      END SUBROUTINE CHPEVX
