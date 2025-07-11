!*==cgeqr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CGEQR
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEQR( M, N, A, LDA, T, TSIZE, WORK, LWORK,
!                         INFO )
!
!       .. Scalar Arguments ..
!       INTEGER           INFO, LDA, M, N, TSIZE, LWORK
!       ..
!       .. Array Arguments ..
!       COMPLEX           A( LDA, * ), T( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEQR computes a QR factorization of a complex M-by-N matrix A:
!>
!>    A = Q * ( R ),
!>            ( 0 )
!>
!> where:
!>
!>    Q is a M-by-M orthogonal matrix;
!>    R is an upper-triangular N-by-N matrix;
!>    0 is a (M-N)-by-N zero matrix, if M > N.
!>
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
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the elements on and above the diagonal of the array
!>          contain the min(M,N)-by-N upper trapezoidal matrix R
!>          (R is upper triangular if M >= N);
!>          the elements below the diagonal are used to store part of the
!>          data structure to represent Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX array, dimension (MAX(5,TSIZE))
!>          On exit, if INFO = 0, T(1) returns optimal (or either minimal
!>          or optimal, if query is assumed) TSIZE. See TSIZE for details.
!>          Remaining T contains part of the data structure used to represent Q.
!>          If one wants to apply or construct Q, then one needs to keep T
!>          (in addition to A) and pass it to further subroutines.
!> \endverbatim
!>
!> \param[in] TSIZE
!> \verbatim
!>          TSIZE is INTEGER
!>          If TSIZE >= 5, the dimension of the array T.
!>          If TSIZE = -1 or -2, then a workspace query is assumed. The routine
!>          only calculates the sizes of the T and WORK arrays, returns these
!>          values as the first entries of the T and WORK arrays, and no error
!>          message related to T or WORK is issued by XERBLA.
!>          If TSIZE = -1, the routine calculates optimal size of T for the
!>          optimum performance and returns this value in T(1).
!>          If TSIZE = -2, the routine calculates minimal size of T and
!>          returns this value in T(1).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          (workspace) COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) contains optimal (or either minimal
!>          or optimal, if query was assumed) LWORK.
!>          See LWORK for details.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If LWORK = -1 or -2, then a workspace query is assumed. The routine
!>          only calculates the sizes of the T and WORK arrays, returns these
!>          values as the first entries of the T and WORK arrays, and no error
!>          message related to T or WORK is issued by XERBLA.
!>          If LWORK = -1, the routine calculates optimal size of WORK for the
!>          optimal performance and returns this value in WORK(1).
!>          If LWORK = -2, the routine calculates minimal size of WORK and
!>          returns this value in WORK(1).
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
!> \par Further Details
!  ====================
!>
!> \verbatim
!>
!> The goal of the interface is to give maximum freedom to the developers for
!> creating any QR factorization algorithm they wish. The triangular
!> (trapezoidal) R has to be stored in the upper part of A. The lower part of A
!> and the array T can be used to store any relevant information for applying or
!> constructing the Q factor. The WORK array can safely be discarded after exit.
!>
!> Caution: One should not expect the sizes of T and WORK to be the same from one
!> LAPACK implementation to the other, or even from one execution to the other.
!> A workspace query (for T and WORK) is needed at each execution. However,
!> for a given execution, the size of T and WORK are fixed and will not change
!> from one query to the next.
!>
!> \endverbatim
!>
!> \par Further Details particular to this LAPACK implementation:
!  ==============================================================
!>
!> \verbatim
!>
!> These details are particular for this LAPACK implementation. Users should not
!> take them for granted. These details may change in the future, and are not likely
!> true for another LAPACK implementation. These details are relevant if one wants
!> to try to understand the code. They are not part of the interface.
!>
!> In this version,
!>
!>          T(2): row block size (MB)
!>          T(3): column block size (NB)
!>          T(6:TSIZE): data structure needed for Q, computed by
!>                           CLATSQR or CGEQRT
!>
!>  Depending on the matrix dimensions M and N, and row and column
!>  block sizes MB and NB returned by ILAENV, CGEQR will use either
!>  CLATSQR (if the matrix is tall-and-skinny) or CGEQRT to compute
!>  the QR factorization.
!>
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CGEQR(M,N,A,Lda,T,Tsize,Work,Lwork,Info)
      IMPLICIT NONE
!*--CGEQR175
!
!  -- LAPACK computational routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
!     November 2019
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , M , N , Tsize , Lwork
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , T(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , lminws , mint , minw
      INTEGER mb , nb , mintsz , nblcks
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL CLATSQR , CGEQRT , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , MOD
!     ..
!     .. External Functions ..
      INTEGER ILAENV
      EXTERNAL ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
!
      lquery = (Tsize==-1 .OR. Tsize==-2 .OR. Lwork==-1 .OR. Lwork==-2)
!
      mint = .FALSE.
      minw = .FALSE.
      IF ( Tsize==-2 .OR. Lwork==-2 ) THEN
         IF ( Tsize/=-1 ) mint = .TRUE.
         IF ( Lwork/=-1 ) minw = .TRUE.
      ENDIF
!
!     Determine the block size
!
      IF ( MIN(M,N)>0 ) THEN
         mb = ILAENV(1,'CGEQR ',' ',M,N,1,-1)
         nb = ILAENV(1,'CGEQR ',' ',M,N,2,-1)
      ELSE
         mb = M
         nb = 1
      ENDIF
      IF ( mb>M .OR. mb<=N ) mb = M
      IF ( nb>MIN(M,N) .OR. nb<1 ) nb = 1
      mintsz = N + 5
      IF ( mb<=N .OR. M<=N ) THEN
         nblcks = 1
      ELSEIF ( MOD(M-N,mb-N)==0 ) THEN
         nblcks = (M-N)/(mb-N)
      ELSE
         nblcks = (M-N)/(mb-N) + 1
      ENDIF
!
!     Determine if the workspace size satisfies minimal size
!
      lminws = .FALSE.
      IF ( (Tsize<MAX(1,nb*N*nblcks+5) .OR. Lwork<nb*N) .AND. (Lwork>=N)&
     &     .AND. (Tsize>=mintsz) .AND. (.NOT.lquery) ) THEN
         IF ( Tsize<MAX(1,nb*N*nblcks+5) ) THEN
            lminws = .TRUE.
            nb = 1
            mb = M
         ENDIF
         IF ( Lwork<nb*N ) THEN
            lminws = .TRUE.
            nb = 1
         ENDIF
      ENDIF
!
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ELSEIF ( Tsize<MAX(1,nb*N*nblcks+5) .AND. (.NOT.lquery) .AND.     &
     &         (.NOT.lminws) ) THEN
         Info = -6
      ELSEIF ( (Lwork<MAX(1,N*nb)) .AND. (.NOT.lquery) .AND.            &
     &         (.NOT.lminws) ) THEN
         Info = -8
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( mint ) THEN
            T(1) = mintsz
         ELSE
            T(1) = nb*N*nblcks + 5
         ENDIF
         T(2) = mb
         T(3) = nb
         IF ( minw ) THEN
            Work(1) = MAX(1,N)
         ELSE
            Work(1) = MAX(1,nb*N)
         ENDIF
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGEQR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N)==0 ) RETURN
!
!     The QR Decomposition
!
      IF ( (M<=N) .OR. (mb<=N) .OR. (mb>=M) ) THEN
         CALL CGEQRT(M,N,nb,A,Lda,T(6),nb,Work,Info)
      ELSE
         CALL CLATSQR(M,N,mb,nb,A,Lda,T(6),nb,Work,Lwork,Info)
      ENDIF
!
      Work(1) = MAX(1,nb*N)
!
!
!     End of CGEQR
!
      END SUBROUTINE CGEQR
