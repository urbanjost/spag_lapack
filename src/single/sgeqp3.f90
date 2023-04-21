!*==sgeqp3.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SGEQP3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEQP3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqp3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqp3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqp3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       REAL               A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGEQP3 computes a QR factorization with column pivoting of a
!> matrix A:  A*P = Q*R  using Level 3 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A. M >= 0.
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
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the upper triangle of the array contains the
!>          min(M,N)-by-N upper trapezoidal matrix R; the elements below
!>          the diagonal, together with the array TAU, represent the
!>          orthogonal matrix Q as a product of min(M,N) elementary
!>          reflectors.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] JPVT
!> \verbatim
!>          JPVT is INTEGER array, dimension (N)
!>          On entry, if JPVT(J).ne.0, the J-th column of A is permuted
!>          to the front of A*P (a leading column); if JPVT(J)=0,
!>          the J-th column of A is a free column.
!>          On exit, if JPVT(J)=K, then the J-th column of A*P was the
!>          the K-th column of A.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO=0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= 3*N+1.
!>          For optimal performance LWORK >= 2*N+( N+1 )*NB, where NB
!>          is the optimal blocksize.
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
!>          = 0: successful exit.
!>          < 0: if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup realGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real/complex vector
!>  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
!>  A(i+1:m,i), and tau in TAU(i).
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
!>    X. Sun, Computer Science Dept., Duke University, USA
!>
!  =====================================================================
      SUBROUTINE SGEQP3(M,N,A,Lda,Jpvt,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!*--SGEQP3155
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Lwork , M , N
!     ..
!     .. Array Arguments ..
      INTEGER Jpvt(*)
      REAL A(Lda,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER INB , INBMIN , IXOVER
      PARAMETER (INB=1,INBMIN=2,IXOVER=3)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery
      INTEGER fjb , iws , j , jb , lwkopt , minmn , minws , na , nb ,   &
     &        nbmin , nfxd , nx , sm , sminmn , sn , topbmn
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEQRF , SLAQP2 , SLAQPS , SORMQR , SSWAP , XERBLA
!     ..
!     .. External Functions ..
      INTEGER ILAENV
      REAL SNRM2
      EXTERNAL ILAENV , SNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC INT , MAX , MIN
!     Test input arguments
!  ====================
!
      Info = 0
      lquery = (Lwork==-1)
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ENDIF
!
      IF ( Info==0 ) THEN
         minmn = MIN(M,N)
         IF ( minmn==0 ) THEN
            iws = 1
            lwkopt = 1
         ELSE
            iws = 3*N + 1
            nb = ILAENV(INB,'SGEQRF',' ',M,N,-1,-1)
            lwkopt = 2*N + (N+1)*nb
         ENDIF
         Work(1) = lwkopt
!
         IF ( (Lwork<iws) .AND. .NOT.lquery ) Info = -8
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGEQP3',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Move initial columns up front.
!
      nfxd = 1
      DO j = 1 , N
         IF ( Jpvt(j)/=0 ) THEN
            IF ( j/=nfxd ) THEN
               CALL SSWAP(M,A(1,j),1,A(1,nfxd),1)
               Jpvt(j) = Jpvt(nfxd)
               Jpvt(nfxd) = j
            ELSE
               Jpvt(j) = j
            ENDIF
            nfxd = nfxd + 1
         ELSE
            Jpvt(j) = j
         ENDIF
      ENDDO
      nfxd = nfxd - 1
!
!     Factorize fixed columns
!  =======================
!
!     Compute the QR factorization of fixed columns and update
!     remaining columns.
!
      IF ( nfxd>0 ) THEN
         na = MIN(M,nfxd)
!CC      CALL SGEQR2( M, NA, A, LDA, TAU, WORK, INFO )
         CALL SGEQRF(M,na,A,Lda,Tau,Work,Lwork,Info)
         iws = MAX(iws,INT(Work(1)))
         IF ( na<N ) THEN
!CC         CALL SORM2R( 'Left', 'Transpose', M, N-NA, NA, A, LDA,
!CC  $                   TAU, A( 1, NA+1 ), LDA, WORK, INFO )
            CALL SORMQR('Left','Transpose',M,N-na,na,A,Lda,Tau,A(1,na+1)&
     &                  ,Lda,Work,Lwork,Info)
            iws = MAX(iws,INT(Work(1)))
         ENDIF
      ENDIF
!
!     Factorize free columns
!  ======================
!
      IF ( nfxd<minmn ) THEN
!
         sm = M - nfxd
         sn = N - nfxd
         sminmn = minmn - nfxd
!
!        Determine the block size.
!
         nb = ILAENV(INB,'SGEQRF',' ',sm,sn,-1,-1)
         nbmin = 2
         nx = 0
!
         IF ( (nb>1) .AND. (nb<sminmn) ) THEN
!
!           Determine when to cross over from blocked to unblocked code.
!
            nx = MAX(0,ILAENV(IXOVER,'SGEQRF',' ',sm,sn,-1,-1))
!
!
            IF ( nx<sminmn ) THEN
!
!              Determine if workspace is large enough for blocked code.
!
               minws = 2*sn + (sn+1)*nb
               iws = MAX(iws,minws)
               IF ( Lwork<minws ) THEN
!
!                 Not enough workspace to use optimal NB: Reduce NB and
!                 determine the minimum value of NB.
!
                  nb = (Lwork-2*sn)/(sn+1)
                  nbmin = MAX(2,ILAENV(INBMIN,'SGEQRF',' ',sm,sn,-1,-1))
!
!
               ENDIF
            ENDIF
         ENDIF
!
!        Initialize partial column norms. The first N elements of work
!        store the exact column norms.
!
         DO j = nfxd + 1 , N
            Work(j) = SNRM2(sm,A(nfxd+1,j),1)
            Work(N+j) = Work(j)
         ENDDO
!
         IF ( (nb>=nbmin) .AND. (nb<sminmn) .AND. (nx<sminmn) ) THEN
!
!           Use blocked code initially.
!
            j = nfxd + 1
!
!           Compute factorization: while loop.
!
!
            topbmn = minmn - nx
            DO WHILE ( j<=topbmn )
               jb = MIN(nb,topbmn-j+1)
!
!              Factorize JB columns among columns J:N.
!
               CALL SLAQPS(M,N-j+1,j-1,jb,fjb,A(1,j),Lda,Jpvt(j),Tau(j),&
     &                     Work(j),Work(N+j),Work(2*N+1),Work(2*N+jb+1),&
     &                     N-j+1)
!
               j = j + fjb
            ENDDO
         ELSE
            j = nfxd + 1
         ENDIF
!
!        Use unblocked code to factor the last or only block.
!
!
         IF ( j<=minmn ) CALL SLAQP2(M,N-j+1,j-1,A(1,j),Lda,Jpvt(j),    &
     &                               Tau(j),Work(j),Work(N+j),          &
     &                               Work(2*N+1))
!
      ENDIF
!
      Work(1) = iws
!
!     End of SGEQP3
!
      END SUBROUTINE SGEQP3
