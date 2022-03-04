!*==zgels.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> ZGELS solves overdetermined or underdetermined systems for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGELS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgels.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgels.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgels.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
!                         INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGELS solves overdetermined or underdetermined complex linear systems
!> involving an M-by-N matrix A, or its conjugate-transpose, using a QR
!> or LQ factorization of A.  It is assumed that A has full rank.
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
!> 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
!>    an underdetermined system A**H * X = B.
!>
!> 4. If TRANS = 'C' and m < n:  find the least squares solution of
!>    an overdetermined system, i.e., solve the least squares problem
!>                 minimize || B - A**H * X ||.
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
!>          = 'C': the linear system involves A**H.
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
!>          columns of the matrices B and X. NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>            if M >= N, A is overwritten by details of its QR
!>                       factorization as returned by ZGEQRF;
!>            if M <  N, A is overwritten by details of its LQ
!>                       factorization as returned by ZGELQF.
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
!>          On entry, the matrix B of right hand side vectors, stored
!>          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
!>          if TRANS = 'C'.
!>          On exit, if INFO = 0, B is overwritten by the solution
!>          vectors, stored columnwise:
!>          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
!>          squares solution vectors; the residual sum of squares for the
!>          solution in each column is given by the sum of squares of the
!>          modulus of elements N+1 to M in that column;
!>          if TRANS = 'N' and m < n, rows 1 to N of B contain the
!>          minimum norm solution vectors;
!>          if TRANS = 'C' and m >= n, rows 1 to M of B contain the
!>          minimum norm solution vectors;
!>          if TRANS = 'C' and m < n, rows 1 to M of B contain the
!>          least squares solution vectors; the residual sum of squares
!>          for the solution in each column is given by the sum of
!>          squares of the modulus of elements M+1 to N in that column.
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
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
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
!> \ingroup complex16GEsolve
!
!  =====================================================================
      SUBROUTINE ZGELS(Trans,M,N,Nrhs,A,Lda,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_DLABAD
      USE S_DLAMCH
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      USE S_ZGELQF
      USE S_ZGEQRF
      USE S_ZLANGE
      USE S_ZLASCL
      USE S_ZLASET
      USE S_ZTRTRS
      USE S_ZUNMLQ
      USE S_ZUNMQR
      IMPLICIT NONE
!*--ZGELS199
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anrm , bignum , bnrm , smlnum
      INTEGER :: brow , i , iascl , ibscl , j , mn , nb , scllen , wsize
      LOGICAL :: lquery , tpsd
      REAL(R8KIND) , DIMENSION(1) :: rwork
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
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      Info = 0
      mn = MIN(M,N)
      lquery = (Lwork==-1)
      IF ( .NOT.(LSAME(Trans,'N') .OR. LSAME(Trans,'C')) ) THEN
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
            nb = ILAENV(1,'ZGEQRF',' ',M,N,-1,-1)
            IF ( tpsd ) THEN
               nb = MAX(nb,ILAENV(1,'ZUNMQR','LN',M,Nrhs,N,-1))
            ELSE
               nb = MAX(nb,ILAENV(1,'ZUNMQR','LC',M,Nrhs,N,-1))
            ENDIF
         ELSE
            nb = ILAENV(1,'ZGELQF',' ',M,N,-1,-1)
            IF ( tpsd ) THEN
               nb = MAX(nb,ILAENV(1,'ZUNMLQ','LC',N,Nrhs,M,-1))
            ELSE
               nb = MAX(nb,ILAENV(1,'ZUNMLQ','LN',N,Nrhs,M,-1))
            ENDIF
         ENDIF
!
         wsize = MAX(1,mn+MAX(mn,Nrhs)*nb)
         Work(1) = DBLE(wsize)
!
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGELS ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N,Nrhs)==0 ) THEN
         CALL ZLASET('Full',MAX(M,N),Nrhs,CZERO,CZERO,B,Ldb)
         RETURN
      ENDIF
!
!     Get machine parameters
!
      smlnum = DLAMCH('S')/DLAMCH('P')
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
!
!     Scale A, B if max element outside range [SMLNUM,BIGNUM]
!
      anrm = ZLANGE('M',M,N,A,Lda,rwork)
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
         GOTO 100
      ENDIF
!
      brow = M
      IF ( tpsd ) brow = N
      bnrm = ZLANGE('M',brow,Nrhs,B,Ldb,rwork)
      ibscl = 0
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL ZLASCL('G',0,0,bnrm,smlnum,brow,Nrhs,B,Ldb,Info)
         ibscl = 1
      ELSEIF ( bnrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL ZLASCL('G',0,0,bnrm,bignum,brow,Nrhs,B,Ldb,Info)
         ibscl = 2
      ENDIF
!
      IF ( M>=N ) THEN
!
!        compute QR factorization of A
!
         CALL ZGEQRF(M,N,A,Lda,Work(1),Work(mn+1),Lwork-mn,Info)
!
!        workspace at least N, optimally N*NB
!
         IF ( .NOT.tpsd ) THEN
!
!           Least-Squares Problem min || A * X - B ||
!
!           B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)
!
            CALL ZUNMQR('Left','Conjugate transpose',M,Nrhs,N,A,Lda,    &
     &                  Work(1),B,Ldb,Work(mn+1),Lwork-mn,Info)
!
!           workspace at least NRHS, optimally NRHS*NB
!
!           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
!
            CALL ZTRTRS('Upper','No transpose','Non-unit',N,Nrhs,A,Lda, &
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
!           B(1:N,1:NRHS) := inv(R**H) * B(1:N,1:NRHS)
!
            CALL ZTRTRS('Upper','Conjugate transpose','Non-unit',N,Nrhs,&
     &                  A,Lda,B,Ldb,Info)
!
            IF ( Info>0 ) RETURN
!
!           B(N+1:M,1:NRHS) = ZERO
!
            DO j = 1 , Nrhs
               DO i = N + 1 , M
                  B(i,j) = CZERO
               ENDDO
            ENDDO
!
!           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
!
            CALL ZUNMQR('Left','No transpose',M,Nrhs,N,A,Lda,Work(1),B, &
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
         CALL ZGELQF(M,N,A,Lda,Work(1),Work(mn+1),Lwork-mn,Info)
!
!        workspace at least M, optimally M*NB.
!
         IF ( .NOT.tpsd ) THEN
!
!           underdetermined system of equations A * X = B
!
!           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
!
            CALL ZTRTRS('Lower','No transpose','Non-unit',M,Nrhs,A,Lda, &
     &                  B,Ldb,Info)
!
            IF ( Info>0 ) RETURN
!
!           B(M+1:N,1:NRHS) = 0
!
            DO j = 1 , Nrhs
               DO i = M + 1 , N
                  B(i,j) = CZERO
               ENDDO
            ENDDO
!
!           B(1:N,1:NRHS) := Q(1:N,:)**H * B(1:M,1:NRHS)
!
            CALL ZUNMLQ('Left','Conjugate transpose',N,Nrhs,M,A,Lda,    &
     &                  Work(1),B,Ldb,Work(mn+1),Lwork-mn,Info)
!
!           workspace at least NRHS, optimally NRHS*NB
!
            scllen = N
!
         ELSE
!
!           overdetermined system min || A**H * X - B ||
!
!           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
!
            CALL ZUNMLQ('Left','No transpose',N,Nrhs,M,A,Lda,Work(1),B, &
     &                  Ldb,Work(mn+1),Lwork-mn,Info)
!
!           workspace at least NRHS, optimally NRHS*NB
!
!           B(1:M,1:NRHS) := inv(L**H) * B(1:M,1:NRHS)
!
            CALL ZTRTRS('Lower','Conjugate transpose','Non-unit',M,Nrhs,&
     &                  A,Lda,B,Ldb,Info)
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
         CALL ZLASCL('G',0,0,anrm,smlnum,scllen,Nrhs,B,Ldb,Info)
      ELSEIF ( iascl==2 ) THEN
         CALL ZLASCL('G',0,0,anrm,bignum,scllen,Nrhs,B,Ldb,Info)
      ENDIF
      IF ( ibscl==1 ) THEN
         CALL ZLASCL('G',0,0,smlnum,bnrm,scllen,Nrhs,B,Ldb,Info)
      ELSEIF ( ibscl==2 ) THEN
         CALL ZLASCL('G',0,0,bignum,bnrm,scllen,Nrhs,B,Ldb,Info)
      ENDIF
!
 100  Work(1) = DBLE(wsize)
!
!
!     End of ZGELS
!
      END SUBROUTINE ZGELS
