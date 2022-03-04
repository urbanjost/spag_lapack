!*==chpgvx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CHPGVX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHPGVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpgvx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpgvx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpgvx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU,
!                          IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK,
!                          IWORK, IFAIL, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, RANGE, UPLO
!       INTEGER            IL, INFO, ITYPE, IU, LDZ, M, N
!       REAL               ABSTOL, VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            IFAIL( * ), IWORK( * )
!       REAL               RWORK( * ), W( * )
!       COMPLEX            AP( * ), BP( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHPGVX computes selected eigenvalues and, optionally, eigenvectors
!> of a complex generalized Hermitian-definite eigenproblem, of the form
!> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
!> B are assumed to be Hermitian, stored in packed format, and B is also
!> positive definite.  Eigenvalues and eigenvectors can be selected by
!> specifying either a range of values or a range of indices for the
!> desired eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          Specifies the problem type to be solved:
!>          = 1:  A*x = (lambda)*B*x
!>          = 2:  A*B*x = (lambda)*x
!>          = 3:  B*A*x = (lambda)*x
!> \endverbatim
!>
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
!>          = 'U':  Upper triangles of A and B are stored;
!>          = 'L':  Lower triangles of A and B are stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B.  N >= 0.
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
!>          On exit, the contents of AP are destroyed.
!> \endverbatim
!>
!> \param[in,out] BP
!> \verbatim
!>          BP is COMPLEX array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the Hermitian matrix
!>          B, packed columnwise in a linear array.  The j-th column of B
!>          is stored in the array BP as follows:
!>          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;
!>          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.
!>
!>          On exit, the triangular factor U or L from the Cholesky
!>          factorization B = U**H*U or B = L*L**H, in the same storage
!>          format as B.
!> \endverbatim
!>
!> \param[in] VL
!> \verbatim
!>          VL is REAL
!>
!>          If RANGE='V', the lower bound of the interval to
!>          be searched for eigenvalues. VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>          VU is REAL
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
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] IU
!> \verbatim
!>          IU is INTEGER
!>
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
!>          On normal exit, the first M elements contain the selected
!>          eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, N)
!>          If JOBZ = 'N', then Z is not referenced.
!>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!>          contain the orthonormal eigenvectors of the matrix A
!>          corresponding to the selected eigenvalues, with the i-th
!>          column of Z holding the eigenvector associated with W(i).
!>          The eigenvectors are normalized as follows:
!>          if ITYPE = 1 or 2, Z**H*B*Z = I;
!>          if ITYPE = 3, Z**H*inv(B)*Z = I.
!>
!>          If an eigenvector fails to converge, then that column of Z
!>          contains the latest approximation to the eigenvector, and the
!>          index of the eigenvector is returned in IFAIL.
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
!>          > 0:  CPPTRF or CHPEVX returned an error code:
!>             <= N:  if INFO = i, CHPEVX failed to converge;
!>                    i eigenvectors failed to converge.  Their indices
!>                    are stored in array IFAIL.
!>             > N:   if INFO = N + i, for 1 <= i <= n, then the leading
!>                    minor of order i of B is not positive definite.
!>                    The factorization of B could not be completed and
!>                    no eigenvalues or eigenvectors were computed.
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
!> \par Contributors:
!  ==================
!>
!>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
!
!  =====================================================================
      SUBROUTINE CHPGVX(Itype,Jobz,Range,Uplo,N,Ap,Bp,Vl,Vu,Il,Iu,      &
     &                  Abstol,M,W,Z,Ldz,Work,Rwork,Iwork,Ifail,Info)
      USE S_CHPEVX
      USE S_CHPGST
      USE S_CPPTRF
      USE S_CTPMV
      USE S_CTPSV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CHPGVX287
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Bp
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: alleig , indeig , upper , valeig , wantz
      INTEGER :: j
      CHARACTER :: trans
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
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
      upper = LSAME(Uplo,'U')
      alleig = LSAME(Range,'A')
      valeig = LSAME(Range,'V')
      indeig = LSAME(Range,'I')
!
      Info = 0
      IF ( Itype<1 .OR. Itype>3 ) THEN
         Info = -1
      ELSEIF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -2
      ELSEIF ( .NOT.(alleig .OR. valeig .OR. indeig) ) THEN
         Info = -3
      ELSEIF ( .NOT.(upper .OR. LSAME(Uplo,'L')) ) THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( valeig ) THEN
         IF ( N>0 .AND. Vu<=Vl ) Info = -9
      ELSEIF ( indeig ) THEN
         IF ( Il<1 ) THEN
            Info = -10
         ELSEIF ( Iu<MIN(N,Il) .OR. Iu>N ) THEN
            Info = -11
         ENDIF
      ENDIF
      IF ( Info==0 ) THEN
         IF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) Info = -16
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CHPGVX',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Form a Cholesky factorization of B.
!
      CALL CPPTRF(Uplo,N,Bp,Info)
      IF ( Info/=0 ) THEN
         Info = N + Info
         RETURN
      ENDIF
!
!     Transform problem to standard eigenvalue problem and solve.
!
      CALL CHPGST(Itype,Uplo,N,Ap,Bp,Info)
      CALL CHPEVX(Jobz,Range,Uplo,N,Ap,Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,    &
     &            Work,Rwork,Iwork,Ifail,Info)
!
      IF ( wantz ) THEN
!
!        Backtransform eigenvectors to the original problem.
!
         IF ( Info>0 ) M = Info - 1
         IF ( Itype==1 .OR. Itype==2 ) THEN
!
!           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
!           backtransform eigenvectors: x = inv(L)**H*y or inv(U)*y
!
            IF ( upper ) THEN
               trans = 'N'
            ELSE
               trans = 'C'
            ENDIF
!
            DO j = 1 , M
               CALL CTPSV(Uplo,trans,'Non-unit',N,Bp,Z(1,j),1)
            ENDDO
!
         ELSEIF ( Itype==3 ) THEN
!
!           For B*A*x=(lambda)*x;
!           backtransform eigenvectors: x = L*y or U**H*y
!
            IF ( upper ) THEN
               trans = 'C'
            ELSE
               trans = 'N'
            ENDIF
!
            DO j = 1 , M
               CALL CTPMV(Uplo,trans,'Non-unit',N,Bp,Z(1,j),1)
            ENDDO
         ENDIF
      ENDIF
!
!
!     End of CHPGVX
!
      END SUBROUTINE CHPGVX
