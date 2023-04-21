!*==sgetsls.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SGETSLS
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB,
!     $                     WORK, LWORK, INFO )
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
!> SGETSLS solves overdetermined or underdetermined real linear systems
!> involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
!> factorization of A.  It is assumed that A has full rank.
!>
!>
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
!>    an undetermined system A**T * X = B.
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
!>          A is overwritten by details of its QR or LQ
!>          factorization as returned by SGEQR or SGELQ.
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
!>          squares solution vectors.
!>          if TRANS = 'N' and m < n, rows 1 to N of B contain the
!>          minimum norm solution vectors;
!>          if TRANS = 'T' and m >= n, rows 1 to M of B contain the
!>          minimum norm solution vectors;
!>          if TRANS = 'T' and m < n, rows 1 to M of B contain the
!>          least squares solution vectors.
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
!>          (workspace) REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) contains optimal (or either minimal
!>          or optimal, if query was assumed) LWORK.
!>          See LWORK for details.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If LWORK = -1 or -2, then a workspace query is assumed.
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
!> \date June 2017
!
!> \ingroup realGEsolve
!
!  =====================================================================
      SUBROUTINE SGETSLS(Trans,M,N,Nrhs,A,Lda,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
!*--SGETSLS165
!
!  -- LAPACK driver routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Info , Lda , Ldb , Lwork , M , N , Nrhs
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , B(Ldb,*) , Work(*)
!
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , tran
      INTEGER i , iascl , ibscl , j , minmn , maxmn , brow , scllen ,   &
     &        mnk , tszo , tszm , lwo , lwm , lw1 , lw2 , wsizeo ,      &
     &        wsizem , info2
      REAL anrm , bignum , bnrm , smlnum , tq(5) , workq(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      REAL SLAMCH , SLANGE
      EXTERNAL LSAME , ILAENV , SLABAD , SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEQR , SGEMQR , SLASCL , SLASET , STRTRS , XERBLA ,     &
     &         SGELQ , SGEMLQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC REAL , MAX , MIN , INT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      Info = 0
      minmn = MIN(M,N)
      maxmn = MAX(M,N)
      mnk = MAX(minmn,Nrhs)
      tran = LSAME(Trans,'T')
!
      lquery = (Lwork==-1 .OR. Lwork==-2)
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
      ENDIF
!
      IF ( Info==0 ) THEN
!
!     Determine the block size and minimum LWORK
!
         IF ( M>=N ) THEN
            CALL SGEQR(M,N,A,Lda,tq,-1,workq,-1,info2)
            tszo = INT(tq(1))
            lwo = INT(workq(1))
            CALL SGEMQR('L',Trans,M,Nrhs,N,A,Lda,tq,tszo,B,Ldb,workq,-1,&
     &                  info2)
            lwo = MAX(lwo,INT(workq(1)))
            CALL SGEQR(M,N,A,Lda,tq,-2,workq,-2,info2)
            tszm = INT(tq(1))
            lwm = INT(workq(1))
            CALL SGEMQR('L',Trans,M,Nrhs,N,A,Lda,tq,tszm,B,Ldb,workq,-1,&
     &                  info2)
            lwm = MAX(lwm,INT(workq(1)))
            wsizeo = tszo + lwo
            wsizem = tszm + lwm
         ELSE
            CALL SGELQ(M,N,A,Lda,tq,-1,workq,-1,info2)
            tszo = INT(tq(1))
            lwo = INT(workq(1))
            CALL SGEMLQ('L',Trans,N,Nrhs,M,A,Lda,tq,tszo,B,Ldb,workq,-1,&
     &                  info2)
            lwo = MAX(lwo,INT(workq(1)))
            CALL SGELQ(M,N,A,Lda,tq,-2,workq,-2,info2)
            tszm = INT(tq(1))
            lwm = INT(workq(1))
            CALL SGEMLQ('L',Trans,N,Nrhs,M,A,Lda,tq,tszm,B,Ldb,workq,-1,&
     &                  info2)
            lwm = MAX(lwm,INT(workq(1)))
            wsizeo = tszo + lwo
            wsizem = tszm + lwm
         ENDIF
!
         IF ( (Lwork<wsizem) .AND. (.NOT.lquery) ) Info = -10
!
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGETSLS',-Info)
         Work(1) = REAL(wsizeo)
         RETURN
      ENDIF
      IF ( lquery ) THEN
         IF ( Lwork==-1 ) Work(1) = REAL(wsizeo)
         IF ( Lwork==-2 ) Work(1) = REAL(wsizem)
         RETURN
      ENDIF
      IF ( Lwork<wsizeo ) THEN
         lw1 = tszm
         lw2 = lwm
      ELSE
         lw1 = tszo
         lw2 = lwo
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N,Nrhs)==0 ) THEN
         CALL SLASET('FULL',MAX(M,N),Nrhs,ZERO,ZERO,B,Ldb)
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
         CALL SLASET('F',maxmn,Nrhs,ZERO,ZERO,B,Ldb)
         GOTO 100
      ENDIF
!
      brow = M
      IF ( tran ) brow = N
      bnrm = SLANGE('M',brow,Nrhs,B,Ldb,Work)
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
         CALL SGEQR(M,N,A,Lda,Work(lw2+1),lw1,Work(1),lw2,Info)
         IF ( .NOT.tran ) THEN
!
!           Least-Squares Problem min || A * X - B ||
!
!           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
!
            CALL SGEMQR('L','T',M,Nrhs,N,A,Lda,Work(lw2+1),lw1,B,Ldb,   &
     &                  Work(1),lw2,Info)
!
!           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
!
            CALL STRTRS('U','N','N',N,Nrhs,A,Lda,B,Ldb,Info)
            IF ( Info>0 ) RETURN
            scllen = N
         ELSE
!
!           Overdetermined system of equations A**T * X = B
!
!           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
!
            CALL STRTRS('U','T','N',N,Nrhs,A,Lda,B,Ldb,Info)
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
            CALL SGEMQR('L','N',M,Nrhs,N,A,Lda,Work(lw2+1),lw1,B,Ldb,   &
     &                  Work(1),lw2,Info)
!
            scllen = M
!
         ENDIF
!
      ELSE
!
!        Compute LQ factorization of A
!
         CALL SGELQ(M,N,A,Lda,Work(lw2+1),lw1,Work(1),lw2,Info)
!
!        workspace at least M, optimally M*NB.
!
         IF ( .NOT.tran ) THEN
!
!           underdetermined system of equations A * X = B
!
!           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
!
            CALL STRTRS('L','N','N',M,Nrhs,A,Lda,B,Ldb,Info)
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
            CALL SGEMLQ('L','T',N,Nrhs,M,A,Lda,Work(lw2+1),lw1,B,Ldb,   &
     &                  Work(1),lw2,Info)
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
            CALL SGEMLQ('L','N',N,Nrhs,M,A,Lda,Work(lw2+1),lw1,B,Ldb,   &
     &                  Work(1),lw2,Info)
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
 100  Work(1) = REAL(tszo+lwo)
!
!     End of SGETSLS
!
      END SUBROUTINE SGETSLS
