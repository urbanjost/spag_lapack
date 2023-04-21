!*==zhetri_3.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
 
!> \brief \b ZHETRI_3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETRI_3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri_3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri_3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri_3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK,
!                            INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ), E( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> ZHETRI_3 computes the inverse of a complex Hermitian indefinite
!> matrix A using the factorization computed by ZHETRF_RK or ZHETRF_BK:
!>
!>     A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T),
!>
!> where U (or L) is unit upper (or lower) triangular matrix,
!> U**H (or L**H) is the conjugate of U (or L), P is a permutation
!> matrix, P**T is the transpose of P, and D is Hermitian and block
!> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> ZHETRI_3 sets the leading dimension of the workspace  before calling
!> ZHETRI_3X that actually computes the inverse.  This is the blocked
!> version of the algorithm, calling Level 3 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are
!>          stored as an upper or lower triangular matrix.
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
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, diagonal of the block diagonal matrix D and
!>          factors U or L as computed by ZHETRF_RK and ZHETRF_BK:
!>            a) ONLY diagonal elements of the Hermitian block diagonal
!>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
!>               (superdiagonal (or subdiagonal) elements of D
!>                should be provided on entry in array E), and
!>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
!>               If UPLO = 'L': factor L in the subdiagonal part of A.
!>
!>          On exit, if INFO = 0, the Hermitian inverse of the original
!>          matrix.
!>             If UPLO = 'U': the upper triangular part of the inverse
!>             is formed and the part of A below the diagonal is not
!>             referenced;
!>             If UPLO = 'L': the lower triangular part of the inverse
!>             is formed and the part of A above the diagonal is not
!>             referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX*16 array, dimension (N)
!>          On entry, contains the superdiagonal (or subdiagonal)
!>          elements of the Hermitian block diagonal matrix D
!>          with 1-by-1 or 2-by-2 diagonal blocks, where
!>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
!>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
!>
!>          NOTE: For 1-by-1 diagonal block D(k), where
!>          1 <= k <= N, the element E(k) is not referenced in both
!>          UPLO = 'U' or UPLO = 'L' cases.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by ZHETRF_RK or ZHETRF_BK.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N+NB+1)*(NB+3).
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of WORK. LWORK >= (N+NB+1)*(NB+3).
!>
!>          If LDWORK = -1, then a workspace query is assumed;
!>          the routine only calculates the optimal size of the optimal
!>          size of the WORK array, returns this value as the first
!>          entry of the WORK array, and no error message related to
!>          LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
!>               inverse could not be computed.
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
!> \date November 2017
!
!> \ingroup complex16HEcomputational
!
!> \par Contributors:
!  ==================
!> \verbatim
!>
!>  November 2017,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE ZHETRI_3(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
!*--ZHETRI_3174
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , Lwork , N
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX*16 A(Lda,*) , E(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL upper , lquery
      INTEGER lwkopt , nb
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      EXTERNAL LSAME , ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL ZHETRI_3X , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1)
!
!     Determine the block size
!
      nb = MAX(1,ILAENV(1,'ZHETRI_3',Uplo,N,-1,-1,-1))
      lwkopt = (N+nb+1)*(nb+3)
!
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Lwork<lwkopt .AND. .NOT.lquery ) THEN
         Info = -8
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHETRI_3',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         Work(1) = lwkopt
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      CALL ZHETRI_3X(Uplo,N,A,Lda,E,Ipiv,Work,nb,Info)
!
      Work(1) = lwkopt
!
!
!     End of ZHETRI_3
!
      END SUBROUTINE ZHETRI_3
