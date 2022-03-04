!*==zunbdb.f90  processed by SPAG 7.51RB at 20:09 on  3 Mar 2022
 
!> \brief \b ZUNBDB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNBDB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunbdb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunbdb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunbdb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12,
!                          X21, LDX21, X22, LDX22, THETA, PHI, TAUP1,
!                          TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIGNS, TRANS
!       INTEGER            INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P,
!      $                   Q
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   PHI( * ), THETA( * )
!       COMPLEX*16         TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ),
!      $                   WORK( * ), X11( LDX11, * ), X12( LDX12, * ),
!      $                   X21( LDX21, * ), X22( LDX22, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNBDB simultaneously bidiagonalizes the blocks of an M-by-M
!> partitioned unitary matrix X:
!>
!>                                 [ B11 | B12 0  0 ]
!>     [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**H
!> X = [-----------] = [---------] [----------------] [---------]   .
!>     [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ]
!>                                 [  0  |  0  0  I ]
!>
!> X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is
!> not the case, then X must be transposed and/or permuted. This can be
!> done in constant time using the TRANS and SIGNS options. See ZUNCSD
!> for details.)
!>
!> The unitary matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by-
!> (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are
!> represented implicitly by Householder vectors.
!>
!> B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented
!> implicitly by angles THETA, PHI.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER
!>          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major
!>                      order;
!>          otherwise:  X, U1, U2, V1T, and V2T are stored in column-
!>                      major order.
!> \endverbatim
!>
!> \param[in] SIGNS
!> \verbatim
!>          SIGNS is CHARACTER
!>          = 'O':      The lower-left block is made nonpositive (the
!>                      "other" convention);
!>          otherwise:  The upper-right block is made nonpositive (the
!>                      "default" convention).
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows and columns in X.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of rows in X11 and X12. 0 <= P <= M.
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is INTEGER
!>          The number of columns in X11 and X21. 0 <= Q <=
!>          MIN(P,M-P,M-Q).
!> \endverbatim
!>
!> \param[in,out] X11
!> \verbatim
!>          X11 is COMPLEX*16 array, dimension (LDX11,Q)
!>          On entry, the top-left block of the unitary matrix to be
!>          reduced. On exit, the form depends on TRANS:
!>          If TRANS = 'N', then
!>             the columns of tril(X11) specify reflectors for P1,
!>             the rows of triu(X11,1) specify reflectors for Q1;
!>          else TRANS = 'T', and
!>             the rows of triu(X11) specify reflectors for P1,
!>             the columns of tril(X11,-1) specify reflectors for Q1.
!> \endverbatim
!>
!> \param[in] LDX11
!> \verbatim
!>          LDX11 is INTEGER
!>          The leading dimension of X11. If TRANS = 'N', then LDX11 >=
!>          P; else LDX11 >= Q.
!> \endverbatim
!>
!> \param[in,out] X12
!> \verbatim
!>          X12 is COMPLEX*16 array, dimension (LDX12,M-Q)
!>          On entry, the top-right block of the unitary matrix to
!>          be reduced. On exit, the form depends on TRANS:
!>          If TRANS = 'N', then
!>             the rows of triu(X12) specify the first P reflectors for
!>             Q2;
!>          else TRANS = 'T', and
!>             the columns of tril(X12) specify the first P reflectors
!>             for Q2.
!> \endverbatim
!>
!> \param[in] LDX12
!> \verbatim
!>          LDX12 is INTEGER
!>          The leading dimension of X12. If TRANS = 'N', then LDX12 >=
!>          P; else LDX11 >= M-Q.
!> \endverbatim
!>
!> \param[in,out] X21
!> \verbatim
!>          X21 is COMPLEX*16 array, dimension (LDX21,Q)
!>          On entry, the bottom-left block of the unitary matrix to
!>          be reduced. On exit, the form depends on TRANS:
!>          If TRANS = 'N', then
!>             the columns of tril(X21) specify reflectors for P2;
!>          else TRANS = 'T', and
!>             the rows of triu(X21) specify reflectors for P2.
!> \endverbatim
!>
!> \param[in] LDX21
!> \verbatim
!>          LDX21 is INTEGER
!>          The leading dimension of X21. If TRANS = 'N', then LDX21 >=
!>          M-P; else LDX21 >= Q.
!> \endverbatim
!>
!> \param[in,out] X22
!> \verbatim
!>          X22 is COMPLEX*16 array, dimension (LDX22,M-Q)
!>          On entry, the bottom-right block of the unitary matrix to
!>          be reduced. On exit, the form depends on TRANS:
!>          If TRANS = 'N', then
!>             the rows of triu(X22(Q+1:M-P,P+1:M-Q)) specify the last
!>             M-P-Q reflectors for Q2,
!>          else TRANS = 'T', and
!>             the columns of tril(X22(P+1:M-Q,Q+1:M-P)) specify the last
!>             M-P-Q reflectors for P2.
!> \endverbatim
!>
!> \param[in] LDX22
!> \verbatim
!>          LDX22 is INTEGER
!>          The leading dimension of X22. If TRANS = 'N', then LDX22 >=
!>          M-P; else LDX22 >= M-Q.
!> \endverbatim
!>
!> \param[out] THETA
!> \verbatim
!>          THETA is DOUBLE PRECISION array, dimension (Q)
!>          The entries of the bidiagonal blocks B11, B12, B21, B22 can
!>          be computed from the angles THETA and PHI. See Further
!>          Details.
!> \endverbatim
!>
!> \param[out] PHI
!> \verbatim
!>          PHI is DOUBLE PRECISION array, dimension (Q-1)
!>          The entries of the bidiagonal blocks B11, B12, B21, B22 can
!>          be computed from the angles THETA and PHI. See Further
!>          Details.
!> \endverbatim
!>
!> \param[out] TAUP1
!> \verbatim
!>          TAUP1 is COMPLEX*16 array, dimension (P)
!>          The scalar factors of the elementary reflectors that define
!>          P1.
!> \endverbatim
!>
!> \param[out] TAUP2
!> \verbatim
!>          TAUP2 is COMPLEX*16 array, dimension (M-P)
!>          The scalar factors of the elementary reflectors that define
!>          P2.
!> \endverbatim
!>
!> \param[out] TAUQ1
!> \verbatim
!>          TAUQ1 is COMPLEX*16 array, dimension (Q)
!>          The scalar factors of the elementary reflectors that define
!>          Q1.
!> \endverbatim
!>
!> \param[out] TAUQ2
!> \verbatim
!>          TAUQ2 is COMPLEX*16 array, dimension (M-Q)
!>          The scalar factors of the elementary reflectors that define
!>          Q2.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= M-Q.
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
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup complex16OTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The bidiagonal blocks B11, B12, B21, and B22 are represented
!>  implicitly by angles THETA(1), ..., THETA(Q) and PHI(1), ...,
!>  PHI(Q-1). B11 and B21 are upper bidiagonal, while B21 and B22 are
!>  lower bidiagonal. Every entry in each bidiagonal band is a product
!>  of a sine or cosine of a THETA with a sine or cosine of a PHI. See
!>  [1] or ZUNCSD for details.
!>
!>  P1, P2, Q1, and Q2 are represented as products of elementary
!>  reflectors. See ZUNCSD for details on generating P1, P2, Q1, and Q2
!>  using ZUNGQR and ZUNGLQ.
!> \endverbatim
!
!> \par References:
!  ================
!>
!>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.
!>      Algorithms, 50(1):33-65, 2009.
!>
!  =====================================================================
      SUBROUTINE ZUNBDB(Trans,Signs,M,P,Q,X11,Ldx11,X12,Ldx12,X21,Ldx21,&
     &                  X22,Ldx22,Theta,Phi,Taup1,Taup2,Tauq1,Tauq2,    &
     &                  Work,Lwork,Info)
      IMPLICIT NONE
!*--ZUNBDB292
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Signs , Trans
      INTEGER Info , Ldx11 , Ldx12 , Ldx21 , Ldx22 , Lwork , M , P , Q
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Phi(*) , Theta(*)
      COMPLEX*16 Taup1(*) , Taup2(*) , Tauq1(*) , Tauq2(*) , Work(*) ,  &
     &           X11(Ldx11,*) , X12(Ldx12,*) , X21(Ldx21,*) ,           &
     &           X22(Ldx22,*)
!     ..
!
!  ====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION REALONE
      PARAMETER (REALONE=1.0D0)
      COMPLEX*16 ONE
      PARAMETER (ONE=(1.0D0,0.0D0))
!     ..
!     .. Local Scalars ..
      LOGICAL colmajor , lquery
      INTEGER i , lworkmin , lworkopt
      DOUBLE PRECISION z1 , z2 , z3 , z4
!     ..
!     .. External Subroutines ..
      EXTERNAL ZAXPY , ZLARF , ZLARFGP , ZSCAL , XERBLA
      EXTERNAL ZLACGV
!
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DZNRM2
      LOGICAL LSAME
      EXTERNAL DZNRM2 , LSAME
!     ..
!     .. Intrinsic Functions
      INTRINSIC ATAN2 , COS , MAX , MIN , SIN
      INTRINSIC DCMPLX , DCONJG
!     ..
!     .. Executable Statements ..
!
!     Test input arguments
!
      Info = 0
      colmajor = .NOT.LSAME(Trans,'T')
      IF ( .NOT.LSAME(Signs,'O') ) THEN
         z1 = REALONE
         z2 = REALONE
         z3 = REALONE
         z4 = REALONE
      ELSE
         z1 = REALONE
         z2 = -REALONE
         z3 = REALONE
         z4 = -REALONE
      ENDIF
      lquery = Lwork== - 1
!
      IF ( M<0 ) THEN
         Info = -3
      ELSEIF ( P<0 .OR. P>M ) THEN
         Info = -4
      ELSEIF ( Q<0 .OR. Q>P .OR. Q>M-P .OR. Q>M-Q ) THEN
         Info = -5
      ELSEIF ( colmajor .AND. Ldx11<MAX(1,P) ) THEN
         Info = -7
      ELSEIF ( .NOT.colmajor .AND. Ldx11<MAX(1,Q) ) THEN
         Info = -7
      ELSEIF ( colmajor .AND. Ldx12<MAX(1,P) ) THEN
         Info = -9
      ELSEIF ( .NOT.colmajor .AND. Ldx12<MAX(1,M-Q) ) THEN
         Info = -9
      ELSEIF ( colmajor .AND. Ldx21<MAX(1,M-P) ) THEN
         Info = -11
      ELSEIF ( .NOT.colmajor .AND. Ldx21<MAX(1,Q) ) THEN
         Info = -11
      ELSEIF ( colmajor .AND. Ldx22<MAX(1,M-P) ) THEN
         Info = -13
      ELSEIF ( .NOT.colmajor .AND. Ldx22<MAX(1,M-Q) ) THEN
         Info = -13
      ENDIF
!
!     Compute workspace
!
      IF ( Info==0 ) THEN
         lworkopt = M - Q
         lworkmin = M - Q
         Work(1) = lworkopt
         IF ( Lwork<lworkmin .AND. .NOT.lquery ) Info = -21
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('xORBDB',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Handle column-major and row-major separately
!
      IF ( colmajor ) THEN
!
!        Reduce columns 1, ..., Q of X11, X12, X21, and X22
!
         DO i = 1 , Q
!
            IF ( i==1 ) THEN
               CALL ZSCAL(P-i+1,DCMPLX(z1,0.0D0),X11(i,i),1)
            ELSE
               CALL ZSCAL(P-i+1,DCMPLX(z1*COS(Phi(i-1)),0.0D0),X11(i,i),&
     &                    1)
               CALL ZAXPY(P-i+1,DCMPLX(-z1*z3*z4*SIN(Phi(i-1)),0.0D0),  &
     &                    X12(i,i-1),1,X11(i,i),1)
            ENDIF
            IF ( i==1 ) THEN
               CALL ZSCAL(M-P-i+1,DCMPLX(z2,0.0D0),X21(i,i),1)
            ELSE
               CALL ZSCAL(M-P-i+1,DCMPLX(z2*COS(Phi(i-1)),0.0D0),       &
     &                    X21(i,i),1)
               CALL ZAXPY(M-P-i+1,DCMPLX(-z2*z3*z4*SIN(Phi(i-1)),0.0D0),&
     &                    X22(i,i-1),1,X21(i,i),1)
            ENDIF
!
            Theta(i) = ATAN2(DZNRM2(M-P-i+1,X21(i,i),1),                &
     &                 DZNRM2(P-i+1,X11(i,i),1))
!
            IF ( P>i ) THEN
               CALL ZLARFGP(P-i+1,X11(i,i),X11(i+1,i),1,Taup1(i))
            ELSEIF ( P==i ) THEN
               CALL ZLARFGP(P-i+1,X11(i,i),X11(i,i),1,Taup1(i))
            ENDIF
            X11(i,i) = ONE
            IF ( M-P>i ) THEN
               CALL ZLARFGP(M-P-i+1,X21(i,i),X21(i+1,i),1,Taup2(i))
            ELSEIF ( M-P==i ) THEN
               CALL ZLARFGP(M-P-i+1,X21(i,i),X21(i,i),1,Taup2(i))
            ENDIF
            X21(i,i) = ONE
!
            IF ( Q>i ) THEN
               CALL ZLARF('L',P-i+1,Q-i,X11(i,i),1,DCONJG(Taup1(i)),    &
     &                    X11(i,i+1),Ldx11,Work)
               CALL ZLARF('L',M-P-i+1,Q-i,X21(i,i),1,DCONJG(Taup2(i)),  &
     &                    X21(i,i+1),Ldx21,Work)
            ENDIF
            IF ( M-Q+1>i ) THEN
               CALL ZLARF('L',P-i+1,M-Q-i+1,X11(i,i),1,DCONJG(Taup1(i)),&
     &                    X12(i,i),Ldx12,Work)
               CALL ZLARF('L',M-P-i+1,M-Q-i+1,X21(i,i),1,               &
     &                    DCONJG(Taup2(i)),X22(i,i),Ldx22,Work)
            ENDIF
!
            IF ( i<Q ) THEN
               CALL ZSCAL(Q-i,DCMPLX(-z1*z3*SIN(Theta(i)),0.0D0),       &
     &                    X11(i,i+1),Ldx11)
               CALL ZAXPY(Q-i,DCMPLX(z2*z3*COS(Theta(i)),0.0D0),        &
     &                    X21(i,i+1),Ldx21,X11(i,i+1),Ldx11)
            ENDIF
            CALL ZSCAL(M-Q-i+1,DCMPLX(-z1*z4*SIN(Theta(i)),0.0D0),      &
     &                 X12(i,i),Ldx12)
            CALL ZAXPY(M-Q-i+1,DCMPLX(z2*z4*COS(Theta(i)),0.0D0),       &
     &                 X22(i,i),Ldx22,X12(i,i),Ldx12)
!
            IF ( i<Q ) Phi(i) = ATAN2(DZNRM2(Q-i,X11(i,i+1),Ldx11),     &
     &                          DZNRM2(M-Q-i+1,X12(i,i),Ldx12))
!
            IF ( i<Q ) THEN
               CALL ZLACGV(Q-i,X11(i,i+1),Ldx11)
               IF ( i==Q-1 ) THEN
                  CALL ZLARFGP(Q-i,X11(i,i+1),X11(i,i+1),Ldx11,Tauq1(i))
               ELSE
                  CALL ZLARFGP(Q-i,X11(i,i+1),X11(i,i+2),Ldx11,Tauq1(i))
               ENDIF
               X11(i,i+1) = ONE
            ENDIF
            IF ( M-Q+1>i ) THEN
               CALL ZLACGV(M-Q-i+1,X12(i,i),Ldx12)
               IF ( M-Q==i ) THEN
                  CALL ZLARFGP(M-Q-i+1,X12(i,i),X12(i,i),Ldx12,Tauq2(i))
               ELSE
                  CALL ZLARFGP(M-Q-i+1,X12(i,i),X12(i,i+1),Ldx12,       &
     &                         Tauq2(i))
               ENDIF
            ENDIF
            X12(i,i) = ONE
!
            IF ( i<Q ) THEN
               CALL ZLARF('R',P-i,Q-i,X11(i,i+1),Ldx11,Tauq1(i),        &
     &                    X11(i+1,i+1),Ldx11,Work)
               CALL ZLARF('R',M-P-i,Q-i,X11(i,i+1),Ldx11,Tauq1(i),      &
     &                    X21(i+1,i+1),Ldx21,Work)
            ENDIF
            IF ( P>i ) CALL ZLARF('R',P-i,M-Q-i+1,X12(i,i),Ldx12,       &
     &                            Tauq2(i),X12(i+1,i),Ldx12,Work)
            IF ( M-P>i ) CALL ZLARF('R',M-P-i,M-Q-i+1,X12(i,i),Ldx12,   &
     &                              Tauq2(i),X22(i+1,i),Ldx22,Work)
!
            IF ( i<Q ) CALL ZLACGV(Q-i,X11(i,i+1),Ldx11)
            CALL ZLACGV(M-Q-i+1,X12(i,i),Ldx12)
!
         ENDDO
!
!        Reduce columns Q + 1, ..., P of X12, X22
!
         DO i = Q + 1 , P
!
            CALL ZSCAL(M-Q-i+1,DCMPLX(-z1*z4,0.0D0),X12(i,i),Ldx12)
            CALL ZLACGV(M-Q-i+1,X12(i,i),Ldx12)
            IF ( i>=M-Q ) THEN
               CALL ZLARFGP(M-Q-i+1,X12(i,i),X12(i,i),Ldx12,Tauq2(i))
            ELSE
               CALL ZLARFGP(M-Q-i+1,X12(i,i),X12(i,i+1),Ldx12,Tauq2(i))
            ENDIF
            X12(i,i) = ONE
!
            IF ( P>i ) CALL ZLARF('R',P-i,M-Q-i+1,X12(i,i),Ldx12,       &
     &                            Tauq2(i),X12(i+1,i),Ldx12,Work)
            IF ( M-P-Q>=1 ) CALL ZLARF('R',M-P-Q,M-Q-i+1,X12(i,i),Ldx12,&
     &                                 Tauq2(i),X22(Q+1,i),Ldx22,Work)
!
            CALL ZLACGV(M-Q-i+1,X12(i,i),Ldx12)
!
         ENDDO
!
!        Reduce columns P + 1, ..., M - Q of X12, X22
!
         DO i = 1 , M - P - Q
!
            CALL ZSCAL(M-P-Q-i+1,DCMPLX(z2*z4,0.0D0),X22(Q+i,P+i),Ldx22)
            CALL ZLACGV(M-P-Q-i+1,X22(Q+i,P+i),Ldx22)
            CALL ZLARFGP(M-P-Q-i+1,X22(Q+i,P+i),X22(Q+i,P+i+1),Ldx22,   &
     &                   Tauq2(P+i))
            X22(Q+i,P+i) = ONE
            CALL ZLARF('R',M-P-Q-i,M-P-Q-i+1,X22(Q+i,P+i),Ldx22,        &
     &                 Tauq2(P+i),X22(Q+i+1,P+i),Ldx22,Work)
!
            CALL ZLACGV(M-P-Q-i+1,X22(Q+i,P+i),Ldx22)
!
         ENDDO
!
      ELSE
!
!        Reduce columns 1, ..., Q of X11, X12, X21, X22
!
         DO i = 1 , Q
!
            IF ( i==1 ) THEN
               CALL ZSCAL(P-i+1,DCMPLX(z1,0.0D0),X11(i,i),Ldx11)
            ELSE
               CALL ZSCAL(P-i+1,DCMPLX(z1*COS(Phi(i-1)),0.0D0),X11(i,i),&
     &                    Ldx11)
               CALL ZAXPY(P-i+1,DCMPLX(-z1*z3*z4*SIN(Phi(i-1)),0.0D0),  &
     &                    X12(i-1,i),Ldx12,X11(i,i),Ldx11)
            ENDIF
            IF ( i==1 ) THEN
               CALL ZSCAL(M-P-i+1,DCMPLX(z2,0.0D0),X21(i,i),Ldx21)
            ELSE
               CALL ZSCAL(M-P-i+1,DCMPLX(z2*COS(Phi(i-1)),0.0D0),       &
     &                    X21(i,i),Ldx21)
               CALL ZAXPY(M-P-i+1,DCMPLX(-z2*z3*z4*SIN(Phi(i-1)),0.0D0),&
     &                    X22(i-1,i),Ldx22,X21(i,i),Ldx21)
            ENDIF
!
            Theta(i) = ATAN2(DZNRM2(M-P-i+1,X21(i,i),Ldx21),            &
     &                 DZNRM2(P-i+1,X11(i,i),Ldx11))
!
            CALL ZLACGV(P-i+1,X11(i,i),Ldx11)
            CALL ZLACGV(M-P-i+1,X21(i,i),Ldx21)
!
            CALL ZLARFGP(P-i+1,X11(i,i),X11(i,i+1),Ldx11,Taup1(i))
            X11(i,i) = ONE
            IF ( i==M-P ) THEN
               CALL ZLARFGP(M-P-i+1,X21(i,i),X21(i,i),Ldx21,Taup2(i))
            ELSE
               CALL ZLARFGP(M-P-i+1,X21(i,i),X21(i,i+1),Ldx21,Taup2(i))
            ENDIF
            X21(i,i) = ONE
!
            CALL ZLARF('R',Q-i,P-i+1,X11(i,i),Ldx11,Taup1(i),X11(i+1,i),&
     &                 Ldx11,Work)
            CALL ZLARF('R',M-Q-i+1,P-i+1,X11(i,i),Ldx11,Taup1(i),       &
     &                 X12(i,i),Ldx12,Work)
            CALL ZLARF('R',Q-i,M-P-i+1,X21(i,i),Ldx21,Taup2(i),         &
     &                 X21(i+1,i),Ldx21,Work)
            CALL ZLARF('R',M-Q-i+1,M-P-i+1,X21(i,i),Ldx21,Taup2(i),     &
     &                 X22(i,i),Ldx22,Work)
!
            CALL ZLACGV(P-i+1,X11(i,i),Ldx11)
            CALL ZLACGV(M-P-i+1,X21(i,i),Ldx21)
!
            IF ( i<Q ) THEN
               CALL ZSCAL(Q-i,DCMPLX(-z1*z3*SIN(Theta(i)),0.0D0),       &
     &                    X11(i+1,i),1)
               CALL ZAXPY(Q-i,DCMPLX(z2*z3*COS(Theta(i)),0.0D0),        &
     &                    X21(i+1,i),1,X11(i+1,i),1)
            ENDIF
            CALL ZSCAL(M-Q-i+1,DCMPLX(-z1*z4*SIN(Theta(i)),0.0D0),      &
     &                 X12(i,i),1)
            CALL ZAXPY(M-Q-i+1,DCMPLX(z2*z4*COS(Theta(i)),0.0D0),       &
     &                 X22(i,i),1,X12(i,i),1)
!
            IF ( i<Q ) Phi(i) = ATAN2(DZNRM2(Q-i,X11(i+1,i),1),DZNRM2(M-&
     &                          Q-i+1,X12(i,i),1))
!
            IF ( i<Q ) THEN
               CALL ZLARFGP(Q-i,X11(i+1,i),X11(i+2,i),1,Tauq1(i))
               X11(i+1,i) = ONE
            ENDIF
            CALL ZLARFGP(M-Q-i+1,X12(i,i),X12(i+1,i),1,Tauq2(i))
            X12(i,i) = ONE
!
            IF ( i<Q ) THEN
               CALL ZLARF('L',Q-i,P-i,X11(i+1,i),1,DCONJG(Tauq1(i)),    &
     &                    X11(i+1,i+1),Ldx11,Work)
               CALL ZLARF('L',Q-i,M-P-i,X11(i+1,i),1,DCONJG(Tauq1(i)),  &
     &                    X21(i+1,i+1),Ldx21,Work)
            ENDIF
            CALL ZLARF('L',M-Q-i+1,P-i,X12(i,i),1,DCONJG(Tauq2(i)),     &
     &                 X12(i,i+1),Ldx12,Work)
            IF ( M-P>i ) CALL ZLARF('L',M-Q-i+1,M-P-i,X12(i,i),1,       &
     &                              DCONJG(Tauq2(i)),X22(i,i+1),Ldx22,  &
     &                              Work)
!
         ENDDO
!
!        Reduce columns Q + 1, ..., P of X12, X22
!
         DO i = Q + 1 , P
!
            CALL ZSCAL(M-Q-i+1,DCMPLX(-z1*z4,0.0D0),X12(i,i),1)
            CALL ZLARFGP(M-Q-i+1,X12(i,i),X12(i+1,i),1,Tauq2(i))
            X12(i,i) = ONE
!
            IF ( P>i ) CALL ZLARF('L',M-Q-i+1,P-i,X12(i,i),1,           &
     &                            DCONJG(Tauq2(i)),X12(i,i+1),Ldx12,    &
     &                            Work)
            IF ( M-P-Q>=1 ) CALL ZLARF('L',M-Q-i+1,M-P-Q,X12(i,i),1,    &
     &                                 DCONJG(Tauq2(i)),X22(i,Q+1),     &
     &                                 Ldx22,Work)
!
         ENDDO
!
!        Reduce columns P + 1, ..., M - Q of X12, X22
!
         DO i = 1 , M - P - Q
!
            CALL ZSCAL(M-P-Q-i+1,DCMPLX(z2*z4,0.0D0),X22(P+i,Q+i),1)
            CALL ZLARFGP(M-P-Q-i+1,X22(P+i,Q+i),X22(P+i+1,Q+i),1,       &
     &                   Tauq2(P+i))
            X22(P+i,Q+i) = ONE
!
            IF ( M-P-Q/=i ) CALL ZLARF('L',M-P-Q-i+1,M-P-Q-i,X22(P+i,Q+i&
     &                                 ),1,DCONJG(Tauq2(P+i)),          &
     &                                 X22(P+i,Q+i+1),Ldx22,Work)
!
         ENDDO
!
      ENDIF
!
!
!     End of ZUNBDB
!
      END SUBROUTINE ZUNBDB
