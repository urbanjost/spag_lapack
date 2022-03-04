!*==zgelsx.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief <b> ZGELSX solves overdetermined or underdetermined systems for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGELSX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgelsx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgelsx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgelsx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,
!                          WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, M, N, NRHS, RANK
!       DOUBLE PRECISION   RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is deprecated and has been replaced by routine ZGELSY.
!>
!> ZGELSX computes the minimum-norm solution to a complex linear least
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
!> by unitary transformations from the right, arriving at the
!> complete orthogonal factorization:
!>    A * P = Q * [ T11 0 ] * Z
!>                [  0  0 ]
!> The minimum-norm solution is then
!>    X = P * Z**H [ inv(T11)*Q1**H*B ]
!>                 [        0         ]
!> where Q1 consists of the first RANK columns of Q.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
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
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the M-by-NRHS right hand side matrix B.
!>          On exit, the N-by-NRHS solution matrix X.
!>          If m >= n and RANK = n, the residual sum-of-squares for
!>          the solution in the i-th column is given by the sum of
!>          squares of elements N+1:M in that column.
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
!>          On entry, if JPVT(i) .ne. 0, the i-th column of A is an
!>          initial column, otherwise it is a free column.  Before
!>          the QR factorization of A, all initial columns are
!>          permuted to the leading positions; only the remaining
!>          free columns are moved as a result of column pivoting
!>          during the factorization.
!>          On exit, if JPVT(i) = k, then the i-th column of A*P
!>          was the k-th column of A.
!> \endverbatim
!>
!> \param[in] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
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
!>          WORK is COMPLEX*16 array, dimension
!>                      (min(M,N) + max( N, 2*min(M,N)+NRHS )),
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (2*N)
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
!> \ingroup complex16GEsolve
!
!  =====================================================================
      SUBROUTINE ZGELSX(M,N,Nrhs,A,Lda,B,Ldb,Jpvt,Rcond,Rank,Work,Rwork,&
     &                  Info)
      IMPLICIT NONE
!*--ZGELSX188
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldb , M , N , Nrhs , Rank
      DOUBLE PRECISION Rcond
!     ..
!     .. Array Arguments ..
      INTEGER Jpvt(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER IMAX , IMIN
      PARAMETER (IMAX=1,IMIN=2)
      DOUBLE PRECISION ZERO , ONE , DONE , NTDONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0,DONE=ZERO,NTDONE=ONE)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , iascl , ibscl , ismax , ismin , j , k , mn
      DOUBLE PRECISION anrm , bignum , bnrm , smax , smaxpr , smin ,    &
     &                 sminpr , smlnum
      COMPLEX*16 c1 , c2 , s1 , s2 , t1 , t2
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZGEQPF , ZLAIC1 , ZLASCL , ZLASET , ZLATZM ,    &
     &         ZTRSM , ZTZRQF , ZUNM2R
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL DLAMCH , ZLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DCONJG , MAX , MIN
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
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGELSX',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N,Nrhs)==0 ) THEN
         Rank = 0
         RETURN
      ENDIF
!
!     Get machine parameters
!
      smlnum = DLAMCH('S')/DLAMCH('P')
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
!
!     Scale A, B if max elements outside range [SMLNUM,BIGNUM]
!
      anrm = ZLANGE('M',M,N,A,Lda,Rwork)
      iascl = 0
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL ZLASCL('G',0,0,anrm,smlnum,M,N,A,Lda,Info)
         iascl = 1
      ELSEIF ( anrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL ZLASCL('G',0,0,anrm,bignum,M,N,A,Lda,Info)
         iascl = 2
      ELSEIF ( anrm==ZERO ) THEN
!
!        Matrix all zero. Return zero solution.
!
         CALL ZLASET('F',MAX(M,N),Nrhs,CZERO,CZERO,B,Ldb)
         Rank = 0
         GOTO 99999
      ENDIF
!
      bnrm = ZLANGE('M',M,Nrhs,B,Ldb,Rwork)
      ibscl = 0
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL ZLASCL('G',0,0,bnrm,smlnum,M,Nrhs,B,Ldb,Info)
         ibscl = 1
      ELSEIF ( bnrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL ZLASCL('G',0,0,bnrm,bignum,M,Nrhs,B,Ldb,Info)
         ibscl = 2
      ENDIF
!
!     Compute QR factorization with column pivoting of A:
!        A * P = Q * R
!
      CALL ZGEQPF(M,N,A,Lda,Jpvt,Work(1),Work(mn+1),Rwork,Info)
!
!     complex workspace MN+N. Real workspace 2*N. Details of Householder
!     rotations stored in WORK(1:MN).
!
!     Determine RANK using incremental condition estimation
!
      Work(ismin) = CONE
      Work(ismax) = CONE
      smax = ABS(A(1,1))
      smin = smax
      IF ( ABS(A(1,1))==ZERO ) THEN
         Rank = 0
         CALL ZLASET('F',MAX(M,N),Nrhs,CZERO,CZERO,B,Ldb)
         GOTO 99999
      ELSE
         Rank = 1
      ENDIF
      DO
!
         IF ( Rank<mn ) THEN
            i = Rank + 1
            CALL ZLAIC1(IMIN,Rank,Work(ismin),smin,A(1,i),A(i,i),sminpr,&
     &                  s1,c1)
            CALL ZLAIC1(IMAX,Rank,Work(ismax),smax,A(1,i),A(i,i),smaxpr,&
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
!     Logically partition R = [ R11 R12 ]
!                             [  0  R22 ]
!     where R11 = R(1:RANK,1:RANK)
!
!     [R11,R12] = [ T11, 0 ] * Y
!
         IF ( Rank<N ) CALL ZTZRQF(Rank,N,A,Lda,Work(mn+1),Info)
!
!     Details of Householder rotations stored in WORK(MN+1:2*MN)
!
!     B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)
!
         CALL ZUNM2R('Left','Conjugate transpose',M,Nrhs,mn,A,Lda,      &
     &               Work(1),B,Ldb,Work(2*mn+1),Info)
!
!     workspace NRHS
!
!      B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)
!
         CALL ZTRSM('Left','Upper','No transpose','Non-unit',Rank,Nrhs, &
     &              CONE,A,Lda,B,Ldb)
!
         DO i = Rank + 1 , N
            DO j = 1 , Nrhs
               B(i,j) = CZERO
            ENDDO
         ENDDO
!
!     B(1:N,1:NRHS) := Y**H * B(1:N,1:NRHS)
!
         IF ( Rank<N ) THEN
            DO i = 1 , Rank
               CALL ZLATZM('Left',N-Rank+1,Nrhs,A(i,Rank+1),Lda,        &
     &                     DCONJG(Work(mn+i)),B(i,1),B(Rank+1,1),Ldb,   &
     &                     Work(2*mn+1))
            ENDDO
         ENDIF
!
!     workspace NRHS
!
!     B(1:N,1:NRHS) := P * B(1:N,1:NRHS)
!
         DO j = 1 , Nrhs
            DO i = 1 , N
               Work(2*mn+i) = NTDONE
            ENDDO
            DO i = 1 , N
               IF ( Work(2*mn+i)==NTDONE ) THEN
                  IF ( Jpvt(i)/=i ) THEN
                     k = i
                     t1 = B(k,j)
                     t2 = B(Jpvt(k),j)
                     DO
                        B(Jpvt(k),j) = t1
                        Work(2*mn+k) = DONE
                        t1 = t2
                        k = Jpvt(k)
                        t2 = B(Jpvt(k),j)
                        IF ( Jpvt(k)==i ) THEN
                           B(i,j) = t1
                           Work(2*mn+k) = DONE
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
!
!     Undo scaling
!
         IF ( iascl==1 ) THEN
            CALL ZLASCL('G',0,0,anrm,smlnum,N,Nrhs,B,Ldb,Info)
            CALL ZLASCL('U',0,0,smlnum,anrm,Rank,Rank,A,Lda,Info)
         ELSEIF ( iascl==2 ) THEN
            CALL ZLASCL('G',0,0,anrm,bignum,N,Nrhs,B,Ldb,Info)
            CALL ZLASCL('U',0,0,bignum,anrm,Rank,Rank,A,Lda,Info)
         ENDIF
         IF ( ibscl==1 ) THEN
            CALL ZLASCL('G',0,0,smlnum,bnrm,N,Nrhs,B,Ldb,Info)
         ELSEIF ( ibscl==2 ) THEN
            CALL ZLASCL('G',0,0,bignum,bnrm,N,Nrhs,B,Ldb,Info)
         ENDIF
         EXIT
      ENDDO
!
!
!
!     End of ZGELSX
!
99999 END SUBROUTINE ZGELSX
