!*==dgbtf2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DGBTF2 computes the LU factorization of a general band matrix using the unblocked version of the algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGBTF2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtf2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtf2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtf2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, KL, KU, LDAB, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGBTF2 computes an LU factorization of a real m-by-n band matrix A
!> using partial pivoting with row interchanges.
!>
!> This is the unblocked version of the algorithm, calling Level 2 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
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
!> \param[in,out] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          On entry, the matrix A in band storage, in rows KL+1 to
!>          2*KL+KU+1; rows 1 to KL of the array need not be set.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!>
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
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!>               has been completed, but the factor U is exactly
!>               singular, and division by zero will occur if it is used
!>               to solve a system of equations.
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
!> \ingroup doubleGBcomputational
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
!>  elements of U, because of fill-in resulting from the row
!>  interchanges.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DGBTF2(M,N,Kl,Ku,Ab,Ldab,Ipiv,Info)
      USE F77KINDS                        
      USE S_DGER
      USE S_DSCAL
      USE S_DSWAP
      USE S_IDAMAX
      USE S_XERBLA
      IMPLICIT NONE
!*--DGBTF2155
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , j , jp , ju , km , kv
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
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
      kv = Ku + Kl
!
!     Test the input parameters.
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Kl<0 ) THEN
         Info = -3
      ELSEIF ( Ku<0 ) THEN
         Info = -4
      ELSEIF ( Ldab<Kl+kv+1 ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGBTF2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) RETURN
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
      DO j = Ku + 2 , MIN(kv,N)
         DO i = kv - j + 2 , Kl
            Ab(i,j) = ZERO
         ENDDO
      ENDDO
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
      ju = 1
!
      DO j = 1 , MIN(M,N)
!
!        Set fill-in elements in column J+KV to zero.
!
         IF ( j+kv<=N ) THEN
            DO i = 1 , Kl
               Ab(i,j+kv) = ZERO
            ENDDO
         ENDIF
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
         km = MIN(Kl,M-j)
         jp = IDAMAX(km+1,Ab(kv+1,j),1)
         Ipiv(j) = jp + j - 1
         IF ( Ab(kv+jp,j)/=ZERO ) THEN
            ju = MAX(ju,MIN(j+Ku+jp-1,N))
!
!           Apply interchange to columns J to JU.
!
            IF ( jp/=1 ) CALL DSWAP(ju-j+1,Ab(kv+jp,j),Ldab-1,Ab(kv+1,j)&
     &                              ,Ldab-1)
!
            IF ( km>0 ) THEN
!
!              Compute multipliers.
!
               CALL DSCAL(km,ONE/Ab(kv+1,j),Ab(kv+2,j),1)
!
!              Update trailing submatrix within the band.
!
               IF ( ju>j ) CALL DGER(km,ju-j,-ONE,Ab(kv+2,j),1,         &
     &                               Ab(kv,j+1),Ldab-1,Ab(kv+1,j+1),    &
     &                               Ldab-1)
            ENDIF
         ELSE
!
!           If pivot is zero, set INFO to the index of the pivot
!           unless a zero pivot has already been found.
!
            IF ( Info==0 ) Info = j
         ENDIF
      ENDDO
!
!     End of DGBTF2
!
      END SUBROUTINE DGBTF2
