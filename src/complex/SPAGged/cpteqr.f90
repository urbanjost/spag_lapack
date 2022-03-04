!*==cpteqr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CPTEQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPTEQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpteqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpteqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpteqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPZ
!       INTEGER            INFO, LDZ, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), WORK( * )
!       COMPLEX            Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPTEQR computes all eigenvalues and, optionally, eigenvectors of a
!> symmetric positive definite tridiagonal matrix by first factoring the
!> matrix using SPTTRF and then calling CBDSQR to compute the singular
!> values of the bidiagonal factor.
!>
!> This routine computes the eigenvalues of the positive definite
!> tridiagonal matrix to high relative accuracy.  This means that if the
!> eigenvalues range over many orders of magnitude in size, then the
!> small eigenvalues and corresponding eigenvectors will be computed
!> more accurately than, for example, with the standard QR method.
!>
!> The eigenvectors of a full or band positive definite Hermitian matrix
!> can also be found if CHETRD, CHPTRD, or CHBTRD has been used to
!> reduce this matrix to tridiagonal form.  (The reduction to
!> tridiagonal form, however, may preclude the possibility of obtaining
!> high relative accuracy in the small eigenvalues of the original
!> matrix, if these eigenvalues range over many orders of magnitude.)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only.
!>          = 'V':  Compute eigenvectors of original Hermitian
!>                  matrix also.  Array Z contains the unitary matrix
!>                  used to reduce the original matrix to tridiagonal
!>                  form.
!>          = 'I':  Compute eigenvectors of tridiagonal matrix also.
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
!>          D is REAL array, dimension (N)
!>          On entry, the n diagonal elements of the tridiagonal matrix.
!>          On normal exit, D contains the eigenvalues, in descending
!>          order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          On entry, the (n-1) subdiagonal elements of the tridiagonal
!>          matrix.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, N)
!>          On entry, if COMPZ = 'V', the unitary matrix used in the
!>          reduction to tridiagonal form.
!>          On exit, if COMPZ = 'V', the orthonormal eigenvectors of the
!>          original Hermitian matrix;
!>          if COMPZ = 'I', the orthonormal eigenvectors of the
!>          tridiagonal matrix.
!>          If INFO > 0 on exit, Z contains the eigenvectors associated
!>          with only the stored eigenvalues.
!>          If  COMPZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          COMPZ = 'V' or 'I', LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (4*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, and i is:
!>                <= N  the Cholesky factorization of the matrix could
!>                      not be performed because the i-th principal minor
!>                      was not positive definite.
!>                > N   the SVD algorithm failed to converge;
!>                      if INFO = N+i, i off-diagonal elements of the
!>                      bidiagonal factor did not converge to zero.
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
!> \date December 2016
!
!> \ingroup complexPTcomputational
!
!  =====================================================================
      SUBROUTINE CPTEQR(Compz,N,D,E,Z,Ldz,Work,Info)
      USE S_CBDSQR
      USE S_CLASET
      USE S_LSAME
      USE S_SPTTRF
      USE S_XERBLA
      IMPLICIT NONE
!*--CPTEQR154
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Compz
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX , DIMENSION(1,1) :: c , vt
      INTEGER :: i , icompz , nru
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  ====================================================================
!
!     .. Parameters ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( LSAME(Compz,'N') ) THEN
         icompz = 0
      ELSEIF ( LSAME(Compz,'V') ) THEN
         icompz = 1
      ELSEIF ( LSAME(Compz,'I') ) THEN
         icompz = 2
      ELSE
         icompz = -1
      ENDIF
      IF ( icompz<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( (Ldz<1) .OR. (icompz>0 .AND. Ldz<MAX(1,N)) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CPTEQR',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( N==1 ) THEN
         IF ( icompz>0 ) Z(1,1) = CONE
         RETURN
      ENDIF
      IF ( icompz==2 ) CALL CLASET('Full',N,N,CZERO,CONE,Z,Ldz)
!
!     Call SPTTRF to factor the matrix.
!
      CALL SPTTRF(N,D,E,Info)
      IF ( Info/=0 ) RETURN
      DO i = 1 , N
         D(i) = SQRT(D(i))
      ENDDO
      DO i = 1 , N - 1
         E(i) = E(i)*D(i)
      ENDDO
!
!     Call CBDSQR to compute the singular values/vectors of the
!     bidiagonal factor.
!
      IF ( icompz>0 ) THEN
         nru = N
      ELSE
         nru = 0
      ENDIF
      CALL CBDSQR('Lower',N,0,nru,0,D,E,vt,1,Z,Ldz,c,1,Work,Info)
!
!     Square the singular values.
!
      IF ( Info==0 ) THEN
         DO i = 1 , N
            D(i) = D(i)*D(i)
         ENDDO
      ELSE
         Info = N + Info
      ENDIF
!
!
!     End of CPTEQR
!
      END SUBROUTINE CPTEQR
