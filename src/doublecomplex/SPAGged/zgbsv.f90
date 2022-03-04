!*==zgbsv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> ZGBSV computes the solution to system of linear equations A * X = B for GB matrices</b> (simple driver)
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGBSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgbsv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgbsv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgbsv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         AB( LDAB, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGBSV computes the solution to a complex system of linear equations
!> A * X = B, where A is a band matrix of order N with KL subdiagonals
!> and KU superdiagonals, and X and B are N-by-NRHS matrices.
!>
!> The LU decomposition with partial pivoting and row interchanges is
!> used to factor A as A = L * U, where L is a product of permutation
!> and unit lower triangular matrices with KL subdiagonals, and U is
!> upper triangular with KL+KU superdiagonals.  The factored form of A
!> is then used to solve the system of equations A * X = B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of linear equations, i.e., the order of the
!>          matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          On entry, the matrix A in band storage, in rows KL+1 to
!>          2*KL+KU+1; rows 1 to KL of the array need not be set.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
!>          On exit, details of the factorization: U is stored as an
!>          upper triangular band matrix with KL+KU superdiagonals in
!>          rows 1 to KL+KU+1, and the multipliers used during the
!>          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!>          See below for further details.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices that define the permutation matrix P;
!>          row i of the matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the N-by-NRHS right hand side matrix B.
!>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
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
!>          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!>                has been completed, but the factor U is exactly
!>                singular, and the solution has not been computed.
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
!> \ingroup complex16GBsolve
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The band storage scheme is illustrated by the following example, when
!>  M = N = 6, KL = 2, KU = 1:
!>
!>  On entry:                       On exit:
!>
!>      *    *    *    +    +    +       *    *    *   u14  u25  u36
!>      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!>     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!>     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!>
!>  Array elements marked * are not used by the routine; elements marked
!>  + need not be set on entry, but are required by the routine to store
!>  elements of U because of fill-in resulting from the row interchanges.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGBSV(N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      USE S_XERBLA
      USE S_ZGBTRF
      USE S_ZGBTRS
      IMPLICIT NONE
!*--ZGBSV170
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Kl<0 ) THEN
         Info = -2
      ELSEIF ( Ku<0 ) THEN
         Info = -3
      ELSEIF ( Nrhs<0 ) THEN
         Info = -4
      ELSEIF ( Ldab<2*Kl+Ku+1 ) THEN
         Info = -6
      ELSEIF ( Ldb<MAX(N,1) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGBSV ',-Info)
         RETURN
      ENDIF
!
!     Compute the LU factorization of the band matrix A.
!
      CALL ZGBTRF(N,N,Kl,Ku,Ab,Ldab,Ipiv,Info)
!
!        Solve the system A*X = B, overwriting B with X.
!
      IF ( Info==0 ) CALL ZGBTRS('No transpose',N,Kl,Ku,Nrhs,Ab,Ldab,   &
     &                           Ipiv,B,Ldb,Info)
!
!     End of ZGBSV
!
      END SUBROUTINE ZGBSV
