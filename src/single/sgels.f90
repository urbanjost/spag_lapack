!*==sgels.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> SGELS solves overdetermined or underdetermined systems for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGELS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgels.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgels.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgels.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
!                         INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGELS solves overdetermined or underdetermined real linear systems
!> involving an M-by-N matrix A, or its transpose, using a QR or LQ
!> factorization of A.  It is assumed that A has full rank.
!>
!> The following options are provided:
!>
!> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
!>    an overdetermined system, i.e., solve the least squares problem
!>                 minimize || B - A*X ||.
!>
!> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
!>    an underdetermined system A * X = B.
!>
!> 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
!>    an underdetermined system A**T * X = B.
!>
!> 4. If TRANS = 'T' and m < n:  find the least squares solution of
!>    an overdetermined system, i.e., solve the least squares problem
!>                 minimize || B - A**T * X ||.
!>
!> Several right hand side vectors b and solution vectors x can be
!> handled in a single call; they are stored as the columns of the
!> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
!> matrix X.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': the linear system involves A;
!>          = 'T': the linear system involves A**T.
!> \endverbatim
!>
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
!>          columns of the matrices B and X. NRHS >=0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit,
!>            if M >= N, A is overwritten by details of its QR
!>                       factorization as returned by SGEQRF;
!>            if M <  N, A is overwritten by details of its LQ
!>                       factorization as returned by SGELQF.
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
!>          On entry, the matrix B of right hand side vectors, stored
!>          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
!>          if TRANS = 'T'.
!>          On exit, if INFO = 0, B is overwritten by the solution
!>          vectors, stored columnwise:
!>          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
!>          squares solution vectors; the residual sum of squares for the
!>          solution in each column is given by the sum of squares of
!>          elements N+1 to M in that column;
!>          if TRANS = 'N' and m < n, rows 1 to N of B contain the
!>          minimum norm solution vectors;
!>          if TRANS = 'T' and m >= n, rows 1 to M of B contain the
!>          minimum norm solution vectors;
!>          if TRANS = 'T' and m < n, rows 1 to M of B contain the
!>          least squares solution vectors; the residual sum of squares
!>          for the solution in each column is given by the sum of
!>          squares of elements M+1 to N in that column.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= MAX(1,M,N).
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
!>          LWORK >= max( 1, MN + max( MN, NRHS ) ).
!>          For optimal performance,
!>          LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
!>          where MN = min(M,N) and NB is the optimum block size.
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
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO =  i, the i-th diagonal element of the
!>                triangular factor of A is zero, so that A does not have
!>                full rank; the least squares solution could not be
!>                computed.
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
!  =====================================================================
      SUBROUTINE SGELS(Trans,M,N,Nrhs,A,Lda,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
!*--SGELS186
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Info , Lda , Ldb , Lwork , M , N , Nrhs
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , B(Ldb,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , tpsd
      INTEGER brow , i , iascl , ibscl , j , mn , nb , scllen , wsize
      REAL anrm , bignum , bnrm , smlnum
!     ..
!     .. Local Arrays ..
      REAL rwork(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      REAL SLAMCH , SLANGE
      EXTERNAL LSAME , ILAENV , SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SGELQF , SGEQRF , SLABAD , SLASCL , SLASET , SORMLQ ,    &
     &         SORMQR , STRTRS , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      Info = 0
      mn = MIN(M,N)
      lquery = (Lwork==-1)
      IF ( .NOT.(LSAME(Trans,'N') .OR. LSAME(Trans,'T')) ) THEN
         Info = -1
      ELSEIF ( M<0 ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Nrhs<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -6
      ELSEIF ( Ldb<MAX(1,M,N) ) THEN
         Info = -8
      ELSEIF ( Lwork<MAX(1,mn+MAX(mn,Nrhs)) .AND. .NOT.lquery ) THEN
         Info = -10
      ENDIF
!
!     Figure out optimal block size
!
      IF ( Info==0 .OR. Info==-10 ) THEN
!
         tpsd = .TRUE.
         IF ( LSAME(Trans,'N') ) tpsd = .FALSE.
!
         IF ( M>=N ) THEN
            nb = ILAENV(1,'SGEQRF',' ',M,N,-1,-1)
            IF ( tpsd ) THEN
               nb = MAX(nb,ILAENV(1,'SORMQR','LN',M,Nrhs,N,-1))
            ELSE
               nb = MAX(nb,ILAENV(1,'SORMQR','LT',M,Nrhs,N,-1))
            ENDIF
         ELSE
            nb = ILAENV(1,'SGELQF',' ',M,N,-1,-1)
            IF ( tpsd ) THEN
               nb = MAX(nb,ILAENV(1,'SORMLQ','LT',N,Nrhs,M,-1))
            ELSE
               nb = MAX(nb,ILAENV(1,'SORMLQ','LN',N,Nrhs,M,-1))
            ENDIF
         ENDIF
!
         wsize = MAX(1,mn+MAX(mn,Nrhs)*nb)
         Work(1) = REAL(wsize)
!
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGELS ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N,Nrhs)==0 ) THEN
         CALL SLASET('Full',MAX(M,N),Nrhs,ZERO,ZERO,B,Ldb)
         RETURN
      ENDIF
!
!     Get machine parameters
!
      smlnum = SLAMCH('S')/SLAMCH('P')
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
!
!     Scale A, B if max element outside range [SMLNUM,BIGNUM]
!
      anrm = SLANGE('M',M,N,A,Lda,rwork)
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
         GOTO 100
      ENDIF
!
      brow = M
      IF ( tpsd ) brow = N
      bnrm = SLANGE('M',brow,Nrhs,B,Ldb,rwork)
      ibscl = 0
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL SLASCL('G',0,0,bnrm,smlnum,brow,Nrhs,B,Ldb,Info)
         ibscl = 1
      ELSEIF ( bnrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL SLASCL('G',0,0,bnrm,bignum,brow,Nrhs,B,Ldb,Info)
         ibscl = 2
      ENDIF
!
      IF ( M>=N ) THEN
!
!        compute QR factorization of A
!
         CALL SGEQRF(M,N,A,Lda,Work(1),Work(mn+1),Lwork-mn,Info)
!
!        workspace at least N, optimally N*NB
!
         IF ( .NOT.tpsd ) THEN
!
!           Least-Squares Problem min || A * X - B ||
!
!           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
!
            CALL SORMQR('Left','Transpose',M,Nrhs,N,A,Lda,Work(1),B,Ldb,&
     &                  Work(mn+1),Lwork-mn,Info)
!
!           workspace at least NRHS, optimally NRHS*NB
!
!           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
!
            CALL STRTRS('Upper','No transpose','Non-unit',N,Nrhs,A,Lda, &
     &                  B,Ldb,Info)
!
            IF ( Info>0 ) RETURN
!
            scllen = N
!
         ELSE
!
!           Underdetermined system of equations A**T * X = B
!
!           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
!
            CALL STRTRS('Upper','Transpose','Non-unit',N,Nrhs,A,Lda,B,  &
     &                  Ldb,Info)
!
            IF ( Info>0 ) RETURN
!
!           B(N+1:M,1:NRHS) = ZERO
!
            DO j = 1 , Nrhs
               DO i = N + 1 , M
                  B(i,j) = ZERO
               ENDDO
            ENDDO
!
!           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
!
            CALL SORMQR('Left','No transpose',M,Nrhs,N,A,Lda,Work(1),B, &
     &                  Ldb,Work(mn+1),Lwork-mn,Info)
!
!           workspace at least NRHS, optimally NRHS*NB
!
            scllen = M
!
         ENDIF
!
      ELSE
!
!        Compute LQ factorization of A
!
         CALL SGELQF(M,N,A,Lda,Work(1),Work(mn+1),Lwork-mn,Info)
!
!        workspace at least M, optimally M*NB.
!
         IF ( .NOT.tpsd ) THEN
!
!           underdetermined system of equations A * X = B
!
!           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
!
            CALL STRTRS('Lower','No transpose','Non-unit',M,Nrhs,A,Lda, &
     &                  B,Ldb,Info)
!
            IF ( Info>0 ) RETURN
!
!           B(M+1:N,1:NRHS) = 0
!
            DO j = 1 , Nrhs
               DO i = M + 1 , N
                  B(i,j) = ZERO
               ENDDO
            ENDDO
!
!           B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)
!
            CALL SORMLQ('Left','Transpose',N,Nrhs,M,A,Lda,Work(1),B,Ldb,&
     &                  Work(mn+1),Lwork-mn,Info)
!
!           workspace at least NRHS, optimally NRHS*NB
!
            scllen = N
!
         ELSE
!
!           overdetermined system min || A**T * X - B ||
!
!           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
!
            CALL SORMLQ('Left','No transpose',N,Nrhs,M,A,Lda,Work(1),B, &
     &                  Ldb,Work(mn+1),Lwork-mn,Info)
!
!           workspace at least NRHS, optimally NRHS*NB
!
!           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
!
            CALL STRTRS('Lower','Transpose','Non-unit',M,Nrhs,A,Lda,B,  &
     &                  Ldb,Info)
!
            IF ( Info>0 ) RETURN
!
            scllen = M
!
         ENDIF
!
      ENDIF
!
!     Undo scaling
!
      IF ( iascl==1 ) THEN
         CALL SLASCL('G',0,0,anrm,smlnum,scllen,Nrhs,B,Ldb,Info)
      ELSEIF ( iascl==2 ) THEN
         CALL SLASCL('G',0,0,anrm,bignum,scllen,Nrhs,B,Ldb,Info)
      ENDIF
      IF ( ibscl==1 ) THEN
         CALL SLASCL('G',0,0,smlnum,bnrm,scllen,Nrhs,B,Ldb,Info)
      ELSEIF ( ibscl==2 ) THEN
         CALL SLASCL('G',0,0,bignum,bnrm,scllen,Nrhs,B,Ldb,Info)
      ENDIF
!
 100  Work(1) = REAL(wsize)
!
!
!     End of SGELS
!
      END SUBROUTINE SGELS
