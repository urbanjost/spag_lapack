!*==cgelq.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGELQ
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGELQ( M, N, A, LDA, T, TSIZE, WORK, LWORK,
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
!> CGELQ computes an LQ factorization of a complex M-by-N matrix A:
!>
!>    A = ( L 0 ) *  Q
!>
!> where:
!>
!>    Q is a N-by-N orthogonal matrix;
!>    L is a lower-triangular M-by-M matrix;
!>    0 is a M-by-(N-M) zero matrix, if M < N.
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
!>          On exit, the elements on and below the diagonal of the array
!>          contain the M-by-min(M,N) lower trapezoidal matrix L
!>          (L is lower triangular if M <= N);
!>          the elements above the diagonal are used to store part of the
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
!> creating any LQ factorization algorithm they wish. The triangular
!> (trapezoidal) L has to be stored in the lower part of A. The lower part of A
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
!>                           CLASWLQ or CGELQT
!>
!>  Depending on the matrix dimensions M and N, and row and column
!>  block sizes MB and NB returned by ILAENV, CGELQ will use either
!>  CLASWLQ (if the matrix is short-and-wide) or CGELQT to compute
!>  the LQ factorization.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CGELQ(M,N,A,Lda,T,Tsize,Work,Lwork,Info)
      USE S_CGELQT
      USE S_CLASWLQ
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CGELQ178
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: lminws , lquery , mint , minw
      INTEGER :: lwmin , lwopt , lwreq , mb , mintsz , nb , nblcks
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
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
         mb = ILAENV(1,'CGELQ ',' ',M,N,1,-1)
         nb = ILAENV(1,'CGELQ ',' ',M,N,2,-1)
      ELSE
         mb = 1
         nb = N
      ENDIF
      IF ( mb>MIN(M,N) .OR. mb<1 ) mb = 1
      IF ( nb>N .OR. nb<=M ) nb = N
      mintsz = M + 5
      IF ( nb<=M .OR. N<=M ) THEN
         nblcks = 1
      ELSEIF ( MOD(N-M,nb-M)==0 ) THEN
         nblcks = (N-M)/(nb-M)
      ELSE
         nblcks = (N-M)/(nb-M) + 1
      ENDIF
!
!     Determine if the workspace size satisfies minimal size
!
      IF ( (N<=M) .OR. (nb<=M) .OR. (nb>=N) ) THEN
         lwmin = MAX(1,N)
         lwopt = MAX(1,mb*N)
      ELSE
         lwmin = MAX(1,M)
         lwopt = MAX(1,mb*M)
      ENDIF
      lminws = .FALSE.
      IF ( (Tsize<MAX(1,mb*M*nblcks+5) .OR. Lwork<lwopt) .AND.          &
     &     (Lwork>=lwmin) .AND. (Tsize>=mintsz) .AND. (.NOT.lquery) )   &
     &     THEN
         IF ( Tsize<MAX(1,mb*M*nblcks+5) ) THEN
            lminws = .TRUE.
            mb = 1
            nb = N
         ENDIF
         IF ( Lwork<lwopt ) THEN
            lminws = .TRUE.
            mb = 1
         ENDIF
      ENDIF
      IF ( (N<=M) .OR. (nb<=M) .OR. (nb>=N) ) THEN
         lwreq = MAX(1,mb*N)
      ELSE
         lwreq = MAX(1,mb*M)
      ENDIF
!
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ELSEIF ( Tsize<MAX(1,mb*M*nblcks+5) .AND. (.NOT.lquery) .AND.     &
     &         (.NOT.lminws) ) THEN
         Info = -6
      ELSEIF ( (Lwork<lwreq) .AND. (.NOT.lquery) .AND. (.NOT.lminws) )  &
     &         THEN
         Info = -8
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( mint ) THEN
            T(1) = mintsz
         ELSE
            T(1) = mb*M*nblcks + 5
         ENDIF
         T(2) = mb
         T(3) = nb
         IF ( minw ) THEN
            Work(1) = lwmin
         ELSE
            Work(1) = lwreq
         ENDIF
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGELQ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N)==0 ) RETURN
!
!     The LQ Decomposition
!
      IF ( (N<=M) .OR. (nb<=M) .OR. (nb>=N) ) THEN
         CALL CGELQT(M,N,mb,A,Lda,T(6),mb,Work,Info)
      ELSE
         CALL CLASWLQ(M,N,mb,nb,A,Lda,T(6),mb,Work,Lwork,Info)
      ENDIF
!
      Work(1) = lwreq
!
!
!     End of CGELQ
!
      END SUBROUTINE CGELQ
