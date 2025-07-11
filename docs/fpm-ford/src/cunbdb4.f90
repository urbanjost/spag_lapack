!*==cunbdb4.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
 
!> \brief \b CUNBDB4
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CUNBDB4 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunbdb4.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunbdb4.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunbdb4.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNBDB4( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI,
!                           TAUP1, TAUP2, TAUQ1, PHANTOM, WORK, LWORK,
!                           INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21
!       ..
!       .. Array Arguments ..
!       REAL               PHI(*), THETA(*)
!       COMPLEX            PHANTOM(*), TAUP1(*), TAUP2(*), TAUQ1(*),
!      $                   WORK(*), X11(LDX11,*), X21(LDX21,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!>\verbatim
!>
!> CUNBDB4 simultaneously bidiagonalizes the blocks of a tall and skinny
!> matrix X with orthonomal columns:
!>
!>                            [ B11 ]
!>      [ X11 ]   [ P1 |    ] [  0  ]
!>      [-----] = [---------] [-----] Q1**T .
!>      [ X21 ]   [    | P2 ] [ B21 ]
!>                            [  0  ]
!>
!> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P,
!> M-P, or Q. Routines CUNBDB1, CUNBDB2, and CUNBDB3 handle cases in
!> which M-Q is not the minimum dimension.
!>
!> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
!> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
!> Householder vectors.
!>
!> B11 and B12 are (M-Q)-by-(M-Q) bidiagonal matrices represented
!> implicitly by angles THETA, PHI.
!>
!>\endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           The number of rows X11 plus the number of rows in X21.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>           The number of rows in X11. 0 <= P <= M.
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is INTEGER
!>           The number of columns in X11 and X21. 0 <= Q <= M and
!>           M-Q <= min(P,M-P,Q).
!> \endverbatim
!>
!> \param[in,out] X11
!> \verbatim
!>          X11 is COMPLEX array, dimension (LDX11,Q)
!>           On entry, the top block of the matrix X to be reduced. On
!>           exit, the columns of tril(X11) specify reflectors for P1 and
!>           the rows of triu(X11,1) specify reflectors for Q1.
!> \endverbatim
!>
!> \param[in] LDX11
!> \verbatim
!>          LDX11 is INTEGER
!>           The leading dimension of X11. LDX11 >= P.
!> \endverbatim
!>
!> \param[in,out] X21
!> \verbatim
!>          X21 is COMPLEX array, dimension (LDX21,Q)
!>           On entry, the bottom block of the matrix X to be reduced. On
!>           exit, the columns of tril(X21) specify reflectors for P2.
!> \endverbatim
!>
!> \param[in] LDX21
!> \verbatim
!>          LDX21 is INTEGER
!>           The leading dimension of X21. LDX21 >= M-P.
!> \endverbatim
!>
!> \param[out] THETA
!> \verbatim
!>          THETA is REAL array, dimension (Q)
!>           The entries of the bidiagonal blocks B11, B21 are defined by
!>           THETA and PHI. See Further Details.
!> \endverbatim
!>
!> \param[out] PHI
!> \verbatim
!>          PHI is REAL array, dimension (Q-1)
!>           The entries of the bidiagonal blocks B11, B21 are defined by
!>           THETA and PHI. See Further Details.
!> \endverbatim
!>
!> \param[out] TAUP1
!> \verbatim
!>          TAUP1 is COMPLEX array, dimension (P)
!>           The scalar factors of the elementary reflectors that define
!>           P1.
!> \endverbatim
!>
!> \param[out] TAUP2
!> \verbatim
!>          TAUP2 is COMPLEX array, dimension (M-P)
!>           The scalar factors of the elementary reflectors that define
!>           P2.
!> \endverbatim
!>
!> \param[out] TAUQ1
!> \verbatim
!>          TAUQ1 is COMPLEX array, dimension (Q)
!>           The scalar factors of the elementary reflectors that define
!>           Q1.
!> \endverbatim
!>
!> \param[out] PHANTOM
!> \verbatim
!>          PHANTOM is COMPLEX array, dimension (M)
!>           The routine computes an M-by-1 column vector Y that is
!>           orthogonal to the columns of [ X11; X21 ]. PHANTOM(1:P) and
!>           PHANTOM(P+1:M) contain Householder vectors for Y(1:P) and
!>           Y(P+1:M), respectively.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK. LWORK >= M-Q.
!>
!>           If LWORK = -1, then a workspace query is assumed; the routine
!>           only calculates the optimal size of the WORK array, returns
!>           this value as the first entry of the WORK array, and no error
!>           message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           = 0:  successful exit.
!>           < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \date July 2012
!
!> \ingroup complexOTHERcomputational
!
!> \par Further Details:
!  =====================
 
!> \verbatim
!>
!>  The upper-bidiagonal blocks B11, B21 are represented implicitly by
!>  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry
!>  in each bidiagonal band is a product of a sine or cosine of a THETA
!>  with a sine or cosine of a PHI. See [1] or CUNCSD for details.
!>
!>  P1, P2, and Q1 are represented as products of elementary reflectors.
!>  See CUNCSD2BY1 for details on generating P1, P2, and Q1 using CUNGQR
!>  and CUNGLQ.
!> \endverbatim
!
!> \par References:
!  ================
!>
!>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.
!>      Algorithms, 50(1):33-65, 2009.
!>
!  =====================================================================
      SUBROUTINE CUNBDB4(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Phantom,Work,Lwork,Info)
      IMPLICIT NONE
!*--CUNBDB4217
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     July 2012
!
!     .. Scalar Arguments ..
      INTEGER Info , Lwork , M , P , Q , Ldx11 , Ldx21
!     ..
!     .. Array Arguments ..
      REAL Phi(*) , Theta(*)
      COMPLEX Phantom(*) , Taup1(*) , Taup2(*) , Tauq1(*) , Work(*) ,   &
     &        X11(Ldx11,*) , X21(Ldx21,*)
!     ..
!
!  ====================================================================
!
!     .. Parameters ..
      COMPLEX NEGONE , ONE , ZERO
      PARAMETER (NEGONE=(-1.0E0,0.0E0),ONE=(1.0E0,0.0E0),               &
     &           ZERO=(0.0E0,0.0E0))
!     ..
!     .. Local Scalars ..
      REAL c , s
      INTEGER childinfo , i , ilarf , iorbdb5 , j , llarf , lorbdb5 ,   &
     &        lworkmin , lworkopt
      LOGICAL lquery
!     ..
!     .. External Subroutines ..
      EXTERNAL CLARF , CLARFGP , CUNBDB5 , CSROT , CSCAL , CLACGV ,     &
     &         XERBLA
!     ..
!     .. External Functions ..
      REAL SCNRM2
      EXTERNAL SCNRM2
!     ..
!     .. Intrinsic Function ..
      INTRINSIC ATAN2 , COS , MAX , SIN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test input arguments
!
      Info = 0
      lquery = Lwork== - 1
!
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( P<M-Q .OR. M-P<M-Q ) THEN
         Info = -2
      ELSEIF ( Q<M-Q .OR. Q>M ) THEN
         Info = -3
      ELSEIF ( Ldx11<MAX(1,P) ) THEN
         Info = -5
      ELSEIF ( Ldx21<MAX(1,M-P) ) THEN
         Info = -7
      ENDIF
!
!     Compute workspace
!
      IF ( Info==0 ) THEN
         ilarf = 2
         llarf = MAX(Q-1,P-1,M-P-1)
         iorbdb5 = 2
         lorbdb5 = Q
         lworkopt = ilarf + llarf - 1
         lworkopt = MAX(lworkopt,iorbdb5+lorbdb5-1)
         lworkmin = lworkopt
         Work(1) = lworkopt
         IF ( Lwork<lworkmin .AND. .NOT.lquery ) Info = -14
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CUNBDB4',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Reduce columns 1, ..., M-Q of X11 and X21
!
      DO i = 1 , M - Q
!
         IF ( i==1 ) THEN
            DO j = 1 , M
               Phantom(j) = ZERO
            ENDDO
            CALL CUNBDB5(P,M-P,Q,Phantom(1),1,Phantom(P+1),1,X11,Ldx11, &
     &                   X21,Ldx21,Work(iorbdb5),lorbdb5,childinfo)
            CALL CSCAL(P,NEGONE,Phantom(1),1)
            CALL CLARFGP(P,Phantom(1),Phantom(2),1,Taup1(1))
            CALL CLARFGP(M-P,Phantom(P+1),Phantom(P+2),1,Taup2(1))
            Theta(i) = ATAN2(REAL(Phantom(1)),REAL(Phantom(P+1)))
            c = COS(Theta(i))
            s = SIN(Theta(i))
            Phantom(1) = ONE
            Phantom(P+1) = ONE
            CALL CLARF('L',P,Q,Phantom(1),1,CONJG(Taup1(1)),X11,Ldx11,  &
     &                 Work(ilarf))
            CALL CLARF('L',M-P,Q,Phantom(P+1),1,CONJG(Taup2(1)),X21,    &
     &                 Ldx21,Work(ilarf))
         ELSE
            CALL CUNBDB5(P-i+1,M-P-i+1,Q-i+1,X11(i,i-1),1,X21(i,i-1),1, &
     &                   X11(i,i),Ldx11,X21(i,i),Ldx21,Work(iorbdb5),   &
     &                   lorbdb5,childinfo)
            CALL CSCAL(P-i+1,NEGONE,X11(i,i-1),1)
            CALL CLARFGP(P-i+1,X11(i,i-1),X11(i+1,i-1),1,Taup1(i))
            CALL CLARFGP(M-P-i+1,X21(i,i-1),X21(i+1,i-1),1,Taup2(i))
            Theta(i) = ATAN2(REAL(X11(i,i-1)),REAL(X21(i,i-1)))
            c = COS(Theta(i))
            s = SIN(Theta(i))
            X11(i,i-1) = ONE
            X21(i,i-1) = ONE
            CALL CLARF('L',P-i+1,Q-i+1,X11(i,i-1),1,CONJG(Taup1(i)),    &
     &                 X11(i,i),Ldx11,Work(ilarf))
            CALL CLARF('L',M-P-i+1,Q-i+1,X21(i,i-1),1,CONJG(Taup2(i)),  &
     &                 X21(i,i),Ldx21,Work(ilarf))
         ENDIF
!
         CALL CSROT(Q-i+1,X11(i,i),Ldx11,X21(i,i),Ldx21,s,-c)
         CALL CLACGV(Q-i+1,X21(i,i),Ldx21)
         CALL CLARFGP(Q-i+1,X21(i,i),X21(i,i+1),Ldx21,Tauq1(i))
         c = REAL(X21(i,i))
         X21(i,i) = ONE
         CALL CLARF('R',P-i,Q-i+1,X21(i,i),Ldx21,Tauq1(i),X11(i+1,i),   &
     &              Ldx11,Work(ilarf))
         CALL CLARF('R',M-P-i,Q-i+1,X21(i,i),Ldx21,Tauq1(i),X21(i+1,i), &
     &              Ldx21,Work(ilarf))
         CALL CLACGV(Q-i+1,X21(i,i),Ldx21)
         IF ( i<M-Q ) THEN
            s = SQRT(SCNRM2(P-i,X11(i+1,i),1)                           &
     &          **2+SCNRM2(M-P-i,X21(i+1,i),1)**2)
            Phi(i) = ATAN2(s,c)
         ENDIF
!
      ENDDO
!
!     Reduce the bottom-right portion of X11 to [ I 0 ]
!
      DO i = M - Q + 1 , P
         CALL CLACGV(Q-i+1,X11(i,i),Ldx11)
         CALL CLARFGP(Q-i+1,X11(i,i),X11(i,i+1),Ldx11,Tauq1(i))
         X11(i,i) = ONE
         CALL CLARF('R',P-i,Q-i+1,X11(i,i),Ldx11,Tauq1(i),X11(i+1,i),   &
     &              Ldx11,Work(ilarf))
         CALL CLARF('R',Q-P,Q-i+1,X11(i,i),Ldx11,Tauq1(i),X21(M-Q+1,i), &
     &              Ldx21,Work(ilarf))
         CALL CLACGV(Q-i+1,X11(i,i),Ldx11)
      ENDDO
!
!     Reduce the bottom-right portion of X21 to [ 0 I ]
!
      DO i = P + 1 , Q
         CALL CLACGV(Q-i+1,X21(M-Q+i-P,i),Ldx21)
         CALL CLARFGP(Q-i+1,X21(M-Q+i-P,i),X21(M-Q+i-P,i+1),Ldx21,      &
     &                Tauq1(i))
         X21(M-Q+i-P,i) = ONE
         CALL CLARF('R',Q-i,Q-i+1,X21(M-Q+i-P,i),Ldx21,Tauq1(i),        &
     &              X21(M-Q+i-P+1,i),Ldx21,Work(ilarf))
         CALL CLACGV(Q-i+1,X21(M-Q+i-P,i),Ldx21)
      ENDDO
!
!
!     End of CUNBDB4
!
      END SUBROUTINE CUNBDB4
