!*==cgetrs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGETRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGETRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgetrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgetrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgetrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGETRS solves a system of linear equations
!>    A * X = B,  A**T * X = B,  or  A**H * X = B
!> with a general N-by-N matrix A using the LU factorization computed
!> by CGETRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations:
!>          = 'N':  A * X = B     (No transpose)
!>          = 'T':  A**T * X = B  (Transpose)
!>          = 'C':  A**H * X = B  (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The factors L and U from the factorization A = P*L*U
!>          as computed by CGETRF.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from CGETRF; for 1<=i<=N, row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          On entry, the right hand side matrix B.
!>          On exit, the solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complexGEcomputational
!
!  =====================================================================
      SUBROUTINE CGETRS(Trans,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE S_CLASWP
      USE S_CTRSM
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CGETRS129
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: notran
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
      Info = 0
      notran = LSAME(Trans,'N')
      IF ( .NOT.notran .AND. .NOT.LSAME(Trans,'T') .AND.                &
     &     .NOT.LSAME(Trans,'C') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGETRS',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
      IF ( notran ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL CLASWP(Nrhs,B,Ldb,1,N,Ipiv,1)
!
!        Solve L*X = B, overwriting B with X.
!
         CALL CTRSM('Left','Lower','No transpose','Unit',N,Nrhs,ONE,A,  &
     &              Lda,B,Ldb)
!
!        Solve U*X = B, overwriting B with X.
!
         CALL CTRSM('Left','Upper','No transpose','Non-unit',N,Nrhs,ONE,&
     &              A,Lda,B,Ldb)
      ELSE
!
!        Solve A**T * X = B  or A**H * X = B.
!
!        Solve U**T *X = B or U**H *X = B, overwriting B with X.
!
         CALL CTRSM('Left','Upper',Trans,'Non-unit',N,Nrhs,ONE,A,Lda,B, &
     &              Ldb)
!
!        Solve L**T *X = B, or L**H *X = B overwriting B with X.
!
         CALL CTRSM('Left','Lower',Trans,'Unit',N,Nrhs,ONE,A,Lda,B,Ldb)
!
!        Apply row interchanges to the solution vectors.
!
         CALL CLASWP(Nrhs,B,Ldb,1,N,Ipiv,-1)
      ENDIF
!
!
!     End of CGETRS
!
      END SUBROUTINE CGETRS
