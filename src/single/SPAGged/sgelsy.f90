!*==sgelsy.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> SGELSY solves overdetermined or underdetermined systems for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGELSY + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgelsy.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgelsy.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgelsy.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
!       REAL               RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       REAL               A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGELSY computes the minimum-norm solution to a real linear least
!> squares problem:
!>     minimize || A * X - B ||
!> using a complete orthogonal factorization of A.  A is an M-by-N
!> matrix which may be rank-deficient.
!>
!> Several right hand side vectors b and solution vectors x can be
!> handled in a single call; they are stored as the columns of the
!> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
!> matrix X.
!>
!> The routine first computes a QR factorization with column pivoting:
!>     A * P = Q * [ R11 R12 ]
!>                 [  0  R22 ]
!> with R11 defined as the largest leading submatrix whose estimated
!> condition number is less than 1/RCOND.  The order of R11, RANK,
!> is the effective rank of A.
!>
!> Then, R22 is considered to be negligible, and R12 is annihilated
!> by orthogonal transformations from the right, arriving at the
!> complete orthogonal factorization:
!>    A * P = Q * [ T11 0 ] * Z
!>                [  0  0 ]
!> The minimum-norm solution is then
!>    X = P * Z**T [ inv(T11)*Q1**T*B ]
!>                 [        0         ]
!> where Q1 consists of the first RANK columns of Q.
!>
!> This routine is basically identical to the original xGELSX except
!> three differences:
!>   o The call to the subroutine xGEQPF has been substituted by the
!>     the call to the subroutine xGEQP3. This subroutine is a Blas-3
!>     version of the QR factorization with column pivoting.
!>   o Matrix B (the right hand side) is updated with Blas-3.
!>   o The permutation of matrix B (the right hand side) is faster and
!>     more simple.
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
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of
!>          columns of matrices B and X. NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, A has been overwritten by details of its
!>          complete orthogonal factorization.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the M-by-NRHS right hand side matrix B.
!>          On exit, the N-by-NRHS solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,M,N).
!> \endverbatim
!>
!> \param[in,out] JPVT
!> \verbatim
!>          JPVT is INTEGER array, dimension (N)
!>          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
!>          to the front of AP, otherwise column i is a free column.
!>          On exit, if JPVT(i) = k, then the i-th column of AP
!>          was the k-th column of A.
!> \endverbatim
!>
!> \param[in] RCOND
!> \verbatim
!>          RCOND is REAL
!>          RCOND is used to determine the effective rank of A, which
!>          is defined as the order of the largest leading triangular
!>          submatrix R11 in the QR factorization with pivoting of A,
!>          whose estimated condition number < 1/RCOND.
!> \endverbatim
!>
!> \param[out] RANK
!> \verbatim
!>          RANK is INTEGER
!>          The effective rank of A, i.e., the order of the submatrix
!>          R11.  This is the same as the order of the submatrix T11
!>          in the complete orthogonal factorization of A.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          The unblocked strategy requires that:
!>             LWORK >= MAX( MN+3*N+1, 2*MN+NRHS ),
!>          where MN = min( M, N ).
!>          The block algorithm requires that:
!>             LWORK >= MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ),
!>          where NB is an upper bound on the blocksize returned
!>          by ILAENV for the routines SGEQP3, STZRZF, STZRQF, SORMQR,
!>          and SORMRZ.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: If INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup realGEsolve
!
!> \par Contributors:
!  ==================
!>
!>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA \n
!>    E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n
!>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n
!>
!  =====================================================================
      SUBROUTINE SGELSY(M,N,Nrhs,A,Lda,B,Ldb,Jpvt,Rcond,Rank,Work,Lwork,&
     &                  Info)
      USE S_ILAENV
      USE S_SCOPY
      USE S_SGEQP3
      USE S_SLABAD
      USE S_SLAIC1
      USE S_SLAMCH
      USE S_SLANGE
      USE S_SLASCL
      USE S_SLASET
      USE S_SORMQR
      USE S_SORMRZ
      USE S_STRSM
      USE S_STZRZF
      USE S_XERBLA
      IMPLICIT NONE
!*--SGELSY222
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  IMAX = 1 , IMIN = 2
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , DIMENSION(*) :: Jpvt
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: anrm , bignum , bnrm , c1 , c2 , s1 , s2 , smax , smaxpr ,&
     &        smin , sminpr , smlnum , wsize
      INTEGER :: i , iascl , ibscl , ismax , ismin , j , lwkmin ,       &
     &           lwkopt , mn , nb , nb1 , nb2 , nb3 , nb4
      LOGICAL :: lquery
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
      mn = MIN(M,N)
      ismin = mn + 1
      ismax = 2*mn + 1
!
!     Test the input arguments.
!
      Info = 0
      lquery = (Lwork==-1)
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,M,N) ) THEN
         Info = -7
      ENDIF
!
!     Figure out optimal block size
!
      IF ( Info==0 ) THEN
         IF ( mn==0 .OR. Nrhs==0 ) THEN
            lwkmin = 1
            lwkopt = 1
         ELSE
            nb1 = ILAENV(1,'SGEQRF',' ',M,N,-1,-1)
            nb2 = ILAENV(1,'SGERQF',' ',M,N,-1,-1)
            nb3 = ILAENV(1,'SORMQR',' ',M,N,Nrhs,-1)
            nb4 = ILAENV(1,'SORMRQ',' ',M,N,Nrhs,-1)
            nb = MAX(nb1,nb2,nb3,nb4)
            lwkmin = mn + MAX(2*mn,N+1,mn+Nrhs)
            lwkopt = MAX(lwkmin,mn+2*N+nb*(N+1),2*mn+nb*Nrhs)
         ENDIF
         Work(1) = lwkopt
!
         IF ( Lwork<lwkmin .AND. .NOT.lquery ) Info = -12
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGELSY',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( mn==0 .OR. Nrhs==0 ) THEN
         Rank = 0
         RETURN
      ENDIF
!
!     Get machine parameters
!
      smlnum = SLAMCH('S')/SLAMCH('P')
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
!
!     Scale A, B if max entries outside range [SMLNUM,BIGNUM]
!
      anrm = SLANGE('M',M,N,A,Lda,Work)
      iascl = 0
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL SLASCL('G',0,0,anrm,smlnum,M,N,A,Lda,Info)
         iascl = 1
      ELSEIF ( anrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL SLASCL('G',0,0,anrm,bignum,M,N,A,Lda,Info)
         iascl = 2
      ELSEIF ( anrm==ZERO ) THEN
!
!        Matrix all zero. Return zero solution.
!
         CALL SLASET('F',MAX(M,N),Nrhs,ZERO,ZERO,B,Ldb)
         Rank = 0
         GOTO 100
      ENDIF
!
      bnrm = SLANGE('M',M,Nrhs,B,Ldb,Work)
      ibscl = 0
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL SLASCL('G',0,0,bnrm,smlnum,M,Nrhs,B,Ldb,Info)
         ibscl = 1
      ELSEIF ( bnrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL SLASCL('G',0,0,bnrm,bignum,M,Nrhs,B,Ldb,Info)
         ibscl = 2
      ENDIF
!
!     Compute QR factorization with column pivoting of A:
!        A * P = Q * R
!
      CALL SGEQP3(M,N,A,Lda,Jpvt,Work(1),Work(mn+1),Lwork-mn,Info)
      wsize = mn + Work(mn+1)
!
!     workspace: MN+2*N+NB*(N+1).
!     Details of Householder rotations stored in WORK(1:MN).
!
!     Determine RANK using incremental condition estimation
!
      Work(ismin) = ONE
      Work(ismax) = ONE
      smax = ABS(A(1,1))
      smin = smax
      IF ( ABS(A(1,1))==ZERO ) THEN
         Rank = 0
         CALL SLASET('F',MAX(M,N),Nrhs,ZERO,ZERO,B,Ldb)
         GOTO 100
      ELSE
         Rank = 1
      ENDIF
      DO
!
         IF ( Rank<mn ) THEN
            i = Rank + 1
            CALL SLAIC1(IMIN,Rank,Work(ismin),smin,A(1,i),A(i,i),sminpr,&
     &                  s1,c1)
            CALL SLAIC1(IMAX,Rank,Work(ismax),smax,A(1,i),A(i,i),smaxpr,&
     &                  s2,c2)
!
            IF ( smaxpr*Rcond<=sminpr ) THEN
               DO i = 1 , Rank
                  Work(ismin+i-1) = s1*Work(ismin+i-1)
                  Work(ismax+i-1) = s2*Work(ismax+i-1)
               ENDDO
               Work(ismin+Rank) = c1
               Work(ismax+Rank) = c2
               smin = sminpr
               smax = smaxpr
               Rank = Rank + 1
               CYCLE
            ENDIF
         ENDIF
!
!     workspace: 3*MN.
!
!     Logically partition R = [ R11 R12 ]
!                             [  0  R22 ]
!     where R11 = R(1:RANK,1:RANK)
!
!     [R11,R12] = [ T11, 0 ] * Y
!
         IF ( Rank<N ) CALL STZRZF(Rank,N,A,Lda,Work(mn+1),Work(2*mn+1),&
     &                             Lwork-2*mn,Info)
!
!     workspace: 2*MN.
!     Details of Householder rotations stored in WORK(MN+1:2*MN)
!
!     B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
!
         CALL SORMQR('Left','Transpose',M,Nrhs,mn,A,Lda,Work(1),B,Ldb,  &
     &               Work(2*mn+1),Lwork-2*mn,Info)
         wsize = MAX(wsize,2*mn+Work(2*mn+1))
!
!     workspace: 2*MN+NB*NRHS.
!
!     B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)
!
         CALL STRSM('Left','Upper','No transpose','Non-unit',Rank,Nrhs, &
     &              ONE,A,Lda,B,Ldb)
!
         DO j = 1 , Nrhs
            DO i = Rank + 1 , N
               B(i,j) = ZERO
            ENDDO
         ENDDO
!
!     B(1:N,1:NRHS) := Y**T * B(1:N,1:NRHS)
!
         IF ( Rank<N ) CALL SORMRZ('Left','Transpose',N,Nrhs,Rank,      &
     &                             N-Rank,A,Lda,Work(mn+1),B,Ldb,       &
     &                             Work(2*mn+1),Lwork-2*mn,Info)
!
!     workspace: 2*MN+NRHS.
!
!     B(1:N,1:NRHS) := P * B(1:N,1:NRHS)
!
         DO j = 1 , Nrhs
            DO i = 1 , N
               Work(Jpvt(i)) = B(i,j)
            ENDDO
            CALL SCOPY(N,Work(1),1,B(1,j),1)
         ENDDO
!
!     workspace: N.
!
!     Undo scaling
!
         IF ( iascl==1 ) THEN
            CALL SLASCL('G',0,0,anrm,smlnum,N,Nrhs,B,Ldb,Info)
            CALL SLASCL('U',0,0,smlnum,anrm,Rank,Rank,A,Lda,Info)
         ELSEIF ( iascl==2 ) THEN
            CALL SLASCL('G',0,0,anrm,bignum,N,Nrhs,B,Ldb,Info)
            CALL SLASCL('U',0,0,bignum,anrm,Rank,Rank,A,Lda,Info)
         ENDIF
         IF ( ibscl==1 ) THEN
            CALL SLASCL('G',0,0,smlnum,bnrm,N,Nrhs,B,Ldb,Info)
         ELSEIF ( ibscl==2 ) THEN
            CALL SLASCL('G',0,0,bignum,bnrm,N,Nrhs,B,Ldb,Info)
         ENDIF
         EXIT
      ENDDO
!
 100  Work(1) = lwkopt
!
!
!     End of SGELSY
!
      END SUBROUTINE SGELSY
