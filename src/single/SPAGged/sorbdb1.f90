!*==sorbdb1.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SORBDB1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SORBDB1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SORBDB1( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI,
!                           TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21
!       ..
!       .. Array Arguments ..
!       REAL               PHI(*), THETA(*)
!       REAL               TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*),
!      $                   X11(LDX11,*), X21(LDX21,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!>\verbatim
!>
!> SORBDB1 simultaneously bidiagonalizes the blocks of a tall and skinny
!> matrix X with orthonomal columns:
!>
!>                            [ B11 ]
!>      [ X11 ]   [ P1 |    ] [  0  ]
!>      [-----] = [---------] [-----] Q1**T .
!>      [ X21 ]   [    | P2 ] [ B21 ]
!>                            [  0  ]
!>
!> X11 is P-by-Q, and X21 is (M-P)-by-Q. Q must be no larger than P,
!> M-P, or M-Q. Routines SORBDB2, SORBDB3, and SORBDB4 handle cases in
!> which Q is not the minimum dimension.
!>
!> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
!> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
!> Householder vectors.
!>
!> B11 and B12 are Q-by-Q bidiagonal matrices represented implicitly by
!> angles THETA, PHI.
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
!>           The number of columns in X11 and X21. 0 <= Q <=
!>           MIN(P,M-P,M-Q).
!> \endverbatim
!>
!> \param[in,out] X11
!> \verbatim
!>          X11 is REAL array, dimension (LDX11,Q)
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
!>          X21 is REAL array, dimension (LDX21,Q)
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
!>          TAUP1 is REAL array, dimension (P)
!>           The scalar factors of the elementary reflectors that define
!>           P1.
!> \endverbatim
!>
!> \param[out] TAUP2
!> \verbatim
!>          TAUP2 is REAL array, dimension (M-P)
!>           The scalar factors of the elementary reflectors that define
!>           P2.
!> \endverbatim
!>
!> \param[out] TAUQ1
!> \verbatim
!>          TAUQ1 is REAL array, dimension (Q)
!>           The scalar factors of the elementary reflectors that define
!>           Q1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
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
!>
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
!> \ingroup realOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The upper-bidiagonal blocks B11, B21 are represented implicitly by
!>  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry
!>  in each bidiagonal band is a product of a sine or cosine of a THETA
!>  with a sine or cosine of a PHI. See [1] or SORCSD for details.
!>
!>  P1, P2, and Q1 are represented as products of elementary reflectors.
!>  See SORCSD2BY1 for details on generating P1, P2, and Q1 using SORGQR
!>  and SORGLQ.
!> \endverbatim
!
!> \par References:
!  ================
!>
!>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.
!>      Algorithms, 50(1):33-65, 2009.
!>
!  =====================================================================
      SUBROUTINE SORBDB1(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      USE S_SLARF
      USE S_SLARFGP
      USE S_SNRM2
      USE S_SORBDB5
      USE S_SROT
      USE S_XERBLA
      IMPLICIT NONE
!*--SORBDB1213
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      REAL , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL , INTENT(OUT) , DIMENSION(*) :: Phi
      REAL , DIMENSION(*) :: Taup1
      REAL , DIMENSION(*) :: Taup2
      REAL , DIMENSION(*) :: Tauq1
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: c , s
      INTEGER :: childinfo , i , ilarf , iorbdb5 , llarf , lorbdb5 ,    &
     &           lworkmin , lworkopt
      LOGICAL :: lquery
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  ====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Function ..
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
      ELSEIF ( P<Q .OR. M-P<Q ) THEN
         Info = -2
      ELSEIF ( Q<0 .OR. M-Q<Q ) THEN
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
         llarf = MAX(P-1,M-P-1,Q-1)
         iorbdb5 = 2
         lorbdb5 = Q - 2
         lworkopt = MAX(ilarf+llarf-1,iorbdb5+lorbdb5-1)
         lworkmin = lworkopt
         Work(1) = lworkopt
         IF ( Lwork<lworkmin .AND. .NOT.lquery ) Info = -14
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SORBDB1',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Reduce columns 1, ..., Q of X11 and X21
!
      DO i = 1 , Q
!
         CALL SLARFGP(P-i+1,X11(i,i),X11(i+1,i),1,Taup1(i))
         CALL SLARFGP(M-P-i+1,X21(i,i),X21(i+1,i),1,Taup2(i))
         Theta(i) = ATAN2(X21(i,i),X11(i,i))
         c = COS(Theta(i))
         s = SIN(Theta(i))
         X11(i,i) = ONE
         X21(i,i) = ONE
         CALL SLARF('L',P-i+1,Q-i,X11(i,i),1,Taup1(i),X11(i,i+1),Ldx11, &
     &              Work(ilarf))
         CALL SLARF('L',M-P-i+1,Q-i,X21(i,i),1,Taup2(i),X21(i,i+1),     &
     &              Ldx21,Work(ilarf))
!
         IF ( i<Q ) THEN
            CALL SROT(Q-i,X11(i,i+1),Ldx11,X21(i,i+1),Ldx21,c,s)
            CALL SLARFGP(Q-i,X21(i,i+1),X21(i,i+2),Ldx21,Tauq1(i))
            s = X21(i,i+1)
            X21(i,i+1) = ONE
            CALL SLARF('R',P-i,Q-i,X21(i,i+1),Ldx21,Tauq1(i),           &
     &                 X11(i+1,i+1),Ldx11,Work(ilarf))
            CALL SLARF('R',M-P-i,Q-i,X21(i,i+1),Ldx21,Tauq1(i),         &
     &                 X21(i+1,i+1),Ldx21,Work(ilarf))
            c = SQRT(SNRM2(P-i,X11(i+1,i+1),1)                          &
     &          **2+SNRM2(M-P-i,X21(i+1,i+1),1)**2)
            Phi(i) = ATAN2(s,c)
            CALL SORBDB5(P-i,M-P-i,Q-i-1,X11(i+1,i+1),1,X21(i+1,i+1),1, &
     &                   X11(i+1,i+2),Ldx11,X21(i+1,i+2),Ldx21,         &
     &                   Work(iorbdb5),lorbdb5,childinfo)
         ENDIF
!
      ENDDO
!
!
!     End of SORBDB1
!
      END SUBROUTINE SORBDB1
