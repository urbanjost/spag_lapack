!*==sspgv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SSPGV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSPGV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspgv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspgv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspgv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,
!                         INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, ITYPE, LDZ, N
!       ..
!       .. Array Arguments ..
!       REAL               AP( * ), BP( * ), W( * ), WORK( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSPGV computes all the eigenvalues and, optionally, the eigenvectors
!> of a real generalized symmetric-definite eigenproblem, of the form
!> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
!> Here A and B are assumed to be symmetric, stored in packed format,
!> and B is also positive definite.
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
!>          AP is REAL array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the symmetric matrix
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
!>          BP is REAL array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the symmetric matrix
!>          B, packed columnwise in a linear array.  The j-th column of B
!>          is stored in the array BP as follows:
!>          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;
!>          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.
!>
!>          On exit, the triangular factor U or L from the Cholesky
!>          factorization B = U**T*U or B = L*L**T, in the same storage
!>          format as B.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is REAL array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is REAL array, dimension (LDZ, N)
!>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
!>          eigenvectors.  The eigenvectors are normalized as follows:
!>          if ITYPE = 1 or 2, Z**T*B*Z = I;
!>          if ITYPE = 3, Z**T*inv(B)*Z = I.
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
!>          WORK is REAL array, dimension (3*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  SPPTRF or SSPEV returned an error code:
!>             <= N:  if INFO = i, SSPEV failed to converge;
!>                    i off-diagonal elements of an intermediate
!>                    tridiagonal form did not converge to zero.
!>             > N:   if INFO = n + i, for 1 <= i <= n, then the leading
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
!> \date June 2017
!
!> \ingroup realOTHEReigen
!
!  =====================================================================
      SUBROUTINE SSPGV(Itype,Jobz,Uplo,N,Ap,Bp,W,Z,Ldz,Work,Info)
      IMPLICIT NONE
!*--SSPGV163
!
!  -- LAPACK driver routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      CHARACTER Jobz , Uplo
      INTEGER Info , Itype , Ldz , N
!     ..
!     .. Array Arguments ..
      REAL Ap(*) , Bp(*) , W(*) , Work(*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL upper , wantz
      CHARACTER trans
      INTEGER j , neig
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL SPPTRF , SSPEV , SSPGST , STPMV , STPSV , XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      wantz = LSAME(Jobz,'V')
      upper = LSAME(Uplo,'U')
!
      Info = 0
      IF ( Itype<1 .OR. Itype>3 ) THEN
         Info = -1
      ELSEIF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -2
      ELSEIF ( .NOT.(upper .OR. LSAME(Uplo,'L')) ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSPGV ',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Form a Cholesky factorization of B.
!
      CALL SPPTRF(Uplo,N,Bp,Info)
      IF ( Info/=0 ) THEN
         Info = N + Info
         RETURN
      ENDIF
!
!     Transform problem to standard eigenvalue problem and solve.
!
      CALL SSPGST(Itype,Uplo,N,Ap,Bp,Info)
      CALL SSPEV(Jobz,Uplo,N,Ap,W,Z,Ldz,Work,Info)
!
      IF ( wantz ) THEN
!
!        Backtransform eigenvectors to the original problem.
!
         neig = N
         IF ( Info>0 ) neig = Info - 1
         IF ( Itype==1 .OR. Itype==2 ) THEN
!
!           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
!           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
!
            IF ( upper ) THEN
               trans = 'N'
            ELSE
               trans = 'T'
            ENDIF
!
            DO j = 1 , neig
               CALL STPSV(Uplo,trans,'Non-unit',N,Bp,Z(1,j),1)
            ENDDO
!
         ELSEIF ( Itype==3 ) THEN
!
!           For B*A*x=(lambda)*x;
!           backtransform eigenvectors: x = L*y or U**T*y
!
            IF ( upper ) THEN
               trans = 'T'
            ELSE
               trans = 'N'
            ENDIF
!
            DO j = 1 , neig
               CALL STPMV(Uplo,trans,'Non-unit',N,Bp,Z(1,j),1)
            ENDDO
         ENDIF
      ENDIF
!
!     End of SSPGV
!
      END SUBROUTINE SSPGV
