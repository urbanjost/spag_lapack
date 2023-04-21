!*==cgelsy.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> CGELSY solves overdetermined or underdetermined systems for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGELSY + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgelsy.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgelsy.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgelsy.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,
!                          WORK, LWORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
!       REAL               RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGELSY computes the minimum-norm solution to a complex linear least
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
!>
!> This routine is basically identical to the original xGELSX except
!> three differences:
!>   o The permutation of matrix B (the right hand side) is faster and
!>     more simple.
!>   o The call to the subroutine xGEQPF has been substituted by the
!>     the call to the subroutine xGEQP3. This subroutine is a Blas-3
!>     version of the QR factorization with column pivoting.
!>   o Matrix B (the right hand side) is updated with Blas-3.
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          B is COMPLEX array, dimension (LDB,NRHS)
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
!>          On exit, if JPVT(i) = k, then the i-th column of A*P
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
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          The unblocked strategy requires that:
!>            LWORK >= MN + MAX( 2*MN, N+1, MN+NRHS )
!>          where MN = min(M,N).
!>          The block algorithm requires that:
!>            LWORK >= MN + MAX( 2*MN, NB*(N+1), MN+MN*NB, MN+NB*NRHS )
!>          where NB is an upper bound on the blocksize returned
!>          by ILAENV for the routines CGEQP3, CTZRZF, CTZRQF, CUNMQR,
!>          and CUNMRZ.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complexGEsolve
!
!> \par Contributors:
!  ==================
!>
!>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA \n
!>    E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n
!>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n
!>
!  =====================================================================
      SUBROUTINE CGELSY(M,N,Nrhs,A,Lda,B,Ldb,Jpvt,Rcond,Rank,Work,Lwork,&
     &                  Rwork,Info)
      IMPLICIT NONE
!*--CGELSY214
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldb , Lwork , M , N , Nrhs , Rank
      REAL Rcond
!     ..
!     .. Array Arguments ..
      INTEGER Jpvt(*)
      REAL Rwork(*)
      COMPLEX A(Lda,*) , B(Ldb,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER IMAX , IMIN
      PARAMETER (IMAX=1,IMIN=2)
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL lquery
      INTEGER i , iascl , ibscl , ismax , ismin , j , lwkopt , mn , nb ,&
     &        nb1 , nb2 , nb3 , nb4
      REAL anrm , bignum , bnrm , smax , smaxpr , smin , sminpr ,       &
     &     smlnum , wsize
      COMPLEX c1 , c2 , s1 , s2
!     ..
!     .. External Subroutines ..
      EXTERNAL CCOPY , CGEQP3 , CLAIC1 , CLASCL , CLASET , CTRSM ,      &
     &         CTZRZF , CUNMQR , CUNMRZ , SLABAD , XERBLA
!     ..
!     .. External Functions ..
      INTEGER ILAENV
      REAL CLANGE , SLAMCH
      EXTERNAL CLANGE , ILAENV , SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , REAL , CMPLX
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
      nb1 = ILAENV(1,'CGEQRF',' ',M,N,-1,-1)
      nb2 = ILAENV(1,'CGERQF',' ',M,N,-1,-1)
      nb3 = ILAENV(1,'CUNMQR',' ',M,N,Nrhs,-1)
      nb4 = ILAENV(1,'CUNMRQ',' ',M,N,Nrhs,-1)
      nb = MAX(nb1,nb2,nb3,nb4)
      lwkopt = MAX(1,mn+2*N+nb*(N+1),2*mn+nb*Nrhs)
      Work(1) = CMPLX(lwkopt)
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
      ELSEIF ( Lwork<(mn+MAX(2*mn,N+1,mn+Nrhs)) .AND. .NOT.lquery ) THEN
         Info = -12
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGELSY',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
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
      smlnum = SLAMCH('S')/SLAMCH('P')
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
!
!     Scale A, B if max entries outside range [SMLNUM,BIGNUM]
!
      anrm = CLANGE('M',M,N,A,Lda,Rwork)
      iascl = 0
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL CLASCL('G',0,0,anrm,smlnum,M,N,A,Lda,Info)
         iascl = 1
      ELSEIF ( anrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL CLASCL('G',0,0,anrm,bignum,M,N,A,Lda,Info)
         iascl = 2
      ELSEIF ( anrm==ZERO ) THEN
!
!        Matrix all zero. Return zero solution.
!
         CALL CLASET('F',MAX(M,N),Nrhs,CZERO,CZERO,B,Ldb)
         Rank = 0
         GOTO 100
      ENDIF
!
      bnrm = CLANGE('M',M,Nrhs,B,Ldb,Rwork)
      ibscl = 0
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL CLASCL('G',0,0,bnrm,smlnum,M,Nrhs,B,Ldb,Info)
         ibscl = 1
      ELSEIF ( bnrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL CLASCL('G',0,0,bnrm,bignum,M,Nrhs,B,Ldb,Info)
         ibscl = 2
      ENDIF
!
!     Compute QR factorization with column pivoting of A:
!        A * P = Q * R
!
      CALL CGEQP3(M,N,A,Lda,Jpvt,Work(1),Work(mn+1),Lwork-mn,Rwork,Info)
      wsize = mn + REAL(Work(mn+1))
!
!     complex workspace: MN+NB*(N+1). real workspace 2*N.
!     Details of Householder rotations stored in WORK(1:MN).
!
!     Determine RANK using incremental condition estimation
!
      Work(ismin) = CONE
      Work(ismax) = CONE
      smax = ABS(A(1,1))
      smin = smax
      IF ( ABS(A(1,1))==ZERO ) THEN
         Rank = 0
         CALL CLASET('F',MAX(M,N),Nrhs,CZERO,CZERO,B,Ldb)
         GOTO 100
      ELSE
         Rank = 1
      ENDIF
      DO
!
         IF ( Rank<mn ) THEN
            i = Rank + 1
            CALL CLAIC1(IMIN,Rank,Work(ismin),smin,A(1,i),A(i,i),sminpr,&
     &                  s1,c1)
            CALL CLAIC1(IMAX,Rank,Work(ismax),smax,A(1,i),A(i,i),smaxpr,&
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
!     complex workspace: 3*MN.
!
!     Logically partition R = [ R11 R12 ]
!                             [  0  R22 ]
!     where R11 = R(1:RANK,1:RANK)
!
!     [R11,R12] = [ T11, 0 ] * Y
!
         IF ( Rank<N ) CALL CTZRZF(Rank,N,A,Lda,Work(mn+1),Work(2*mn+1),&
     &                             Lwork-2*mn,Info)
!
!     complex workspace: 2*MN.
!     Details of Householder rotations stored in WORK(MN+1:2*MN)
!
!     B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)
!
         CALL CUNMQR('Left','Conjugate transpose',M,Nrhs,mn,A,Lda,      &
     &               Work(1),B,Ldb,Work(2*mn+1),Lwork-2*mn,Info)
         wsize = MAX(wsize,2*mn+REAL(Work(2*mn+1)))
!
!     complex workspace: 2*MN+NB*NRHS.
!
!     B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)
!
         CALL CTRSM('Left','Upper','No transpose','Non-unit',Rank,Nrhs, &
     &              CONE,A,Lda,B,Ldb)
!
         DO j = 1 , Nrhs
            DO i = Rank + 1 , N
               B(i,j) = CZERO
            ENDDO
         ENDDO
!
!     B(1:N,1:NRHS) := Y**H * B(1:N,1:NRHS)
!
         IF ( Rank<N ) CALL CUNMRZ('Left','Conjugate transpose',N,Nrhs, &
     &                             Rank,N-Rank,A,Lda,Work(mn+1),B,Ldb,  &
     &                             Work(2*mn+1),Lwork-2*mn,Info)
!
!     complex workspace: 2*MN+NRHS.
!
!     B(1:N,1:NRHS) := P * B(1:N,1:NRHS)
!
         DO j = 1 , Nrhs
            DO i = 1 , N
               Work(Jpvt(i)) = B(i,j)
            ENDDO
            CALL CCOPY(N,Work(1),1,B(1,j),1)
         ENDDO
!
!     complex workspace: N.
!
!     Undo scaling
!
         IF ( iascl==1 ) THEN
            CALL CLASCL('G',0,0,anrm,smlnum,N,Nrhs,B,Ldb,Info)
            CALL CLASCL('U',0,0,smlnum,anrm,Rank,Rank,A,Lda,Info)
         ELSEIF ( iascl==2 ) THEN
            CALL CLASCL('G',0,0,anrm,bignum,N,Nrhs,B,Ldb,Info)
            CALL CLASCL('U',0,0,bignum,anrm,Rank,Rank,A,Lda,Info)
         ENDIF
         IF ( ibscl==1 ) THEN
            CALL CLASCL('G',0,0,smlnum,bnrm,N,Nrhs,B,Ldb,Info)
         ELSEIF ( ibscl==2 ) THEN
            CALL CLASCL('G',0,0,bignum,bnrm,N,Nrhs,B,Ldb,Info)
         ENDIF
         EXIT
      ENDDO
!
 100  Work(1) = CMPLX(lwkopt)
!
!
!     End of CGELSY
!
      END SUBROUTINE CGELSY
