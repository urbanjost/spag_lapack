!*==dtgex2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DTGEX2 swaps adjacent diagonal blocks in an upper (quasi) triangular matrix pair by an orthogonal equivalence transformation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DTGEX2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgex2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgex2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgex2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
!                          LDZ, J1, N1, N2, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            WANTQ, WANTZ
!       INTEGER            INFO, J1, LDA, LDB, LDQ, LDZ, LWORK, N, N1, N2
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
!      $                   WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTGEX2 swaps adjacent diagonal blocks (A11, B11) and (A22, B22)
!> of size 1-by-1 or 2-by-2 in an upper (quasi) triangular matrix pair
!> (A, B) by an orthogonal equivalence transformation.
!>
!> (A, B) must be in generalized real Schur canonical form (as returned
!> by DGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
!> diagonal blocks. B is upper triangular.
!>
!> Optionally, the matrices Q and Z of generalized Schur vectors are
!> updated.
!>
!>        Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T
!>        Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTQ
!> \verbatim
!>          WANTQ is LOGICAL
!>          .TRUE. : update the left transformation matrix Q;
!>          .FALSE.: do not update Q.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          .TRUE. : update the right transformation matrix Z;
!>          .FALSE.: do not update Z.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B. N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimensions (LDA,N)
!>          On entry, the matrix A in the pair (A, B).
!>          On exit, the updated matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimensions (LDB,N)
!>          On entry, the matrix B in the pair (A, B).
!>          On exit, the updated matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>          On entry, if WANTQ = .TRUE., the orthogonal matrix Q.
!>          On exit, the updated matrix Q.
!>          Not referenced if WANTQ = .FALSE..
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q. LDQ >= 1.
!>          If WANTQ = .TRUE., LDQ >= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ,N)
!>          On entry, if WANTZ =.TRUE., the orthogonal matrix Z.
!>          On exit, the updated matrix Z.
!>          Not referenced if WANTZ = .FALSE..
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= 1.
!>          If WANTZ = .TRUE., LDZ >= N.
!> \endverbatim
!>
!> \param[in] J1
!> \verbatim
!>          J1 is INTEGER
!>          The index to the first block (A11, B11). 1 <= J1 <= N.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>          The order of the first block (A11, B11). N1 = 0, 1 or 2.
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!>          The order of the second block (A22, B22). N2 = 0, 1 or 2.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)).
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          LWORK >=  MAX( 1, N*(N2+N1), (N2+N1)*(N2+N1)*2 )
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>            =0: Successful exit
!>            >0: If INFO = 1, the transformed matrix (A, B) would be
!>                too far from generalized Schur form; the blocks are
!>                not swapped and (A, B) and (Q, Z) are unchanged.
!>                The problem of swapping is too ill-conditioned.
!>            <0: If INFO = -16: LWORK is too small. Appropriate value
!>                for LWORK is returned in WORK(1).
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
!> \ingroup doubleGEauxiliary
!
!> \par Further Details:
!  =====================
!>
!>  In the current code both weak and strong stability tests are
!>  performed. The user can omit the strong stability test by changing
!>  the internal logical parameter WANDS to .FALSE.. See ref. [2] for
!>  details.
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!> \par References:
!  ================
!>
!> \verbatim
!>
!>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
!>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
!>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
!>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
!>
!>  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
!>      Eigenvalues of a Regular Matrix Pair (A, B) and Condition
!>      Estimation: Theory, Algorithms and Software,
!>      Report UMINF - 94.04, Department of Computing Science, Umea
!>      University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working
!>      Note 87. To appear in Numerical Algorithms, 1996.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DTGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,J1,N1,N2, &
     &                  Work,Lwork,Info)
      USE F77KINDS                        
      USE S_DGEMM
      USE S_DGEQR2
      USE S_DGERQ2
      USE S_DLACPY
      USE S_DLAGV2
      USE S_DLAMCH
      USE S_DLARTG
      USE S_DLASET
      USE S_DLASSQ
      USE S_DORG2R
      USE S_DORGR2
      USE S_DORM2R
      USE S_DORMR2
      USE S_DROT
      USE S_DSCAL
      USE S_DTGSY2
      IMPLICIT NONE
!*--DTGEX2242
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWENTY = 2.0D+01
      INTEGER , PARAMETER  ::  LDST = 4
      LOGICAL , PARAMETER  ::  WANDS = .TRUE.
!
! Dummy argument declarations rewritten by SPAG
!
      LOGICAL , INTENT(IN) :: Wantq
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(IN) :: J1
      INTEGER :: N1
      INTEGER :: N2
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) , DIMENSION(2) :: ai , ar , be
      REAL(R8KIND) :: bqra21 , brqa21 , ddum , dnorma , dnormb ,        &
     &                dscale , dsum , eps , f , g , sa , sb , scale ,   &
     &                smlnum , thresha , threshb
      INTEGER :: i , idum , linfo , m
      REAL(R8KIND) , DIMENSION(LDST,LDST) :: ir , ircop , li , licop ,  &
     &            s , scpy , t , tcpy
      INTEGER , DIMENSION(LDST) :: iwork
      LOGICAL :: strong , weak
      REAL(R8KIND) , DIMENSION(LDST) :: taul , taur
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!  Replaced various illegal calls to DCOPY by calls to DLASET, or by DO
!  loops. Sven Hammarling, 1/5/02.
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
      Info = 0
!
!     Quick return if possible
!
      IF ( N<=1 .OR. N1<=0 .OR. N2<=0 ) RETURN
      IF ( N1>N .OR. (J1+N1)>N ) RETURN
      m = N1 + N2
      IF ( Lwork<MAX(1,N*m,m*m*2) ) THEN
         Info = -16
         Work(1) = MAX(1,N*m,m*m*2)
         RETURN
      ENDIF
!
      weak = .FALSE.
      strong = .FALSE.
!
!     Make a local copy of selected block
!
      CALL DLASET('Full',LDST,LDST,ZERO,ZERO,li,LDST)
      CALL DLASET('Full',LDST,LDST,ZERO,ZERO,ir,LDST)
      CALL DLACPY('Full',m,m,A(J1,J1),Lda,s,LDST)
      CALL DLACPY('Full',m,m,B(J1,J1),Ldb,t,LDST)
!
!     Compute threshold for testing acceptance of swapping.
!
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')/eps
      dscale = ZERO
      dsum = ONE
      CALL DLACPY('Full',m,m,s,LDST,Work,m)
      CALL DLASSQ(m*m,Work,1,dscale,dsum)
      dnorma = dscale*SQRT(dsum)
      dscale = ZERO
      dsum = ONE
      CALL DLACPY('Full',m,m,t,LDST,Work,m)
      CALL DLASSQ(m*m,Work,1,dscale,dsum)
      dnormb = dscale*SQRT(dsum)
!
!     THRES has been changed from
!        THRESH = MAX( TEN*EPS*SA, SMLNUM )
!     to
!        THRESH = MAX( TWENTY*EPS*SA, SMLNUM )
!     on 04/01/10.
!     "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by
!     Jim Demmel and Guillaume Revy. See forum post 1783.
!
      thresha = MAX(TWENTY*eps*dnorma,smlnum)
      threshb = MAX(TWENTY*eps*dnormb,smlnum)
!
      IF ( m==2 ) THEN
!
!        CASE 1: Swap 1-by-1 and 1-by-1 blocks.
!
!        Compute orthogonal QL and RQ that swap 1-by-1 and 1-by-1 blocks
!        using Givens rotations and perform the swap tentatively.
!
         f = s(2,2)*t(1,1) - t(2,2)*s(1,1)
         g = s(2,2)*t(1,2) - t(2,2)*s(1,2)
         sa = ABS(s(2,2))*ABS(t(1,1))
         sb = ABS(s(1,1))*ABS(t(2,2))
         CALL DLARTG(f,g,ir(1,2),ir(1,1),ddum)
         ir(2,1) = -ir(1,2)
         ir(2,2) = ir(1,1)
         CALL DROT(2,s(1,1),1,s(1,2),1,ir(1,1),ir(2,1))
         CALL DROT(2,t(1,1),1,t(1,2),1,ir(1,1),ir(2,1))
         IF ( sa>=sb ) THEN
            CALL DLARTG(s(1,1),s(2,1),li(1,1),li(2,1),ddum)
         ELSE
            CALL DLARTG(t(1,1),t(2,1),li(1,1),li(2,1),ddum)
         ENDIF
         CALL DROT(2,s(1,1),LDST,s(2,1),LDST,li(1,1),li(2,1))
         CALL DROT(2,t(1,1),LDST,t(2,1),LDST,li(1,1),li(2,1))
         li(2,2) = li(1,1)
         li(1,2) = -li(2,1)
!
!        Weak stability test: |S21| <= O(EPS F-norm((A)))
!                           and  |T21| <= O(EPS F-norm((B)))
!
         weak = ABS(s(2,1))<=thresha .AND. ABS(t(2,1))<=threshb
         IF ( .NOT.weak ) THEN
!
!     Exit with INFO = 1 if swap was rejected.
!
!
            Info = 1
         ELSE
!
            IF ( WANDS ) THEN
!
!           Strong stability test:
!               F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
!               and
!               F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))
!
               CALL DLACPY('Full',m,m,A(J1,J1),Lda,Work(m*m+1),m)
               CALL DGEMM('N','N',m,m,m,ONE,li,LDST,s,LDST,ZERO,Work,m)
               CALL DGEMM('N','T',m,m,m,-ONE,Work,m,ir,LDST,ONE,        &
     &                    Work(m*m+1),m)
               dscale = ZERO
               dsum = ONE
               CALL DLASSQ(m*m,Work(m*m+1),1,dscale,dsum)
               sa = dscale*SQRT(dsum)
!
               CALL DLACPY('Full',m,m,B(J1,J1),Ldb,Work(m*m+1),m)
               CALL DGEMM('N','N',m,m,m,ONE,li,LDST,t,LDST,ZERO,Work,m)
               CALL DGEMM('N','T',m,m,m,-ONE,Work,m,ir,LDST,ONE,        &
     &                    Work(m*m+1),m)
               dscale = ZERO
               dsum = ONE
               CALL DLASSQ(m*m,Work(m*m+1),1,dscale,dsum)
               sb = dscale*SQRT(dsum)
               strong = sa<=thresha .AND. sb<=threshb
               IF ( .NOT.strong ) THEN
                  Info = 1
                  GOTO 99999
               ENDIF
            ENDIF
!
!        Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and
!               (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)).
!
            CALL DROT(J1+1,A(1,J1),1,A(1,J1+1),1,ir(1,1),ir(2,1))
            CALL DROT(J1+1,B(1,J1),1,B(1,J1+1),1,ir(1,1),ir(2,1))
            CALL DROT(N-J1+1,A(J1,J1),Lda,A(J1+1,J1),Lda,li(1,1),li(2,1)&
     &                )
            CALL DROT(N-J1+1,B(J1,J1),Ldb,B(J1+1,J1),Ldb,li(1,1),li(2,1)&
     &                )
!
!        Set  N1-by-N2 (2,1) - blocks to ZERO.
!
            A(J1+1,J1) = ZERO
            B(J1+1,J1) = ZERO
!
!        Accumulate transformations into Q and Z if requested.
!
            IF ( Wantz ) CALL DROT(N,Z(1,J1),1,Z(1,J1+1),1,ir(1,1),     &
     &                             ir(2,1))
            IF ( Wantq ) CALL DROT(N,Q(1,J1),1,Q(1,J1+1),1,li(1,1),     &
     &                             li(2,1))
!
!        Exit with INFO = 0 if swap was successfully performed.
!
            RETURN
         ENDIF
!
      ELSE
!
!        CASE 2: Swap 1-by-1 and 2-by-2 blocks, or 2-by-2
!                and 2-by-2 blocks.
!
!        Solve the generalized Sylvester equation
!                 S11 * R - L * S22 = SCALE * S12
!                 T11 * R - L * T22 = SCALE * T12
!        for R and L. Solutions in LI and IR.
!
         CALL DLACPY('Full',N1,N2,t(1,N1+1),LDST,li,LDST)
         CALL DLACPY('Full',N1,N2,s(1,N1+1),LDST,ir(N2+1,N1+1),LDST)
         CALL DTGSY2('N',0,N1,N2,s,LDST,s(N1+1,N1+1),LDST,ir(N2+1,N1+1),&
     &               LDST,t,LDST,t(N1+1,N1+1),LDST,li,LDST,scale,dsum,  &
     &               dscale,iwork,idum,linfo)
         IF ( linfo/=0 ) THEN
            Info = 1
         ELSE
!
!        Compute orthogonal matrix QL:
!
!                    QL**T * LI = [ TL ]
!                                 [ 0  ]
!        where
!                    LI =  [      -L              ]
!                          [ SCALE * identity(N2) ]
!
            DO i = 1 , N2
               CALL DSCAL(N1,-ONE,li(1,i),1)
               li(N1+i,i) = scale
            ENDDO
            CALL DGEQR2(m,N2,li,LDST,taul,Work,linfo)
            IF ( linfo/=0 ) THEN
               Info = 1
            ELSE
               CALL DORG2R(m,m,N2,li,LDST,taul,Work,linfo)
               IF ( linfo/=0 ) THEN
                  Info = 1
               ELSE
!
!        Compute orthogonal matrix RQ:
!
!                    IR * RQ**T =   [ 0  TR],
!
!         where IR = [ SCALE * identity(N1), R ]
!
                  DO i = 1 , N1
                     ir(N2+i,i) = scale
                  ENDDO
                  CALL DGERQ2(N1,m,ir(N2+1,1),LDST,taur,Work,linfo)
                  IF ( linfo/=0 ) THEN
                     Info = 1
                  ELSE
                     CALL DORGR2(m,m,N1,ir,LDST,taur,Work,linfo)
                     IF ( linfo/=0 ) THEN
                        Info = 1
                     ELSE
!
!        Perform the swapping tentatively:
!
                        CALL DGEMM('T','N',m,m,m,ONE,li,LDST,s,LDST,    &
     &                             ZERO,Work,m)
                        CALL DGEMM('N','T',m,m,m,ONE,Work,m,ir,LDST,    &
     &                             ZERO,s,LDST)
                        CALL DGEMM('T','N',m,m,m,ONE,li,LDST,t,LDST,    &
     &                             ZERO,Work,m)
                        CALL DGEMM('N','T',m,m,m,ONE,Work,m,ir,LDST,    &
     &                             ZERO,t,LDST)
                        CALL DLACPY('F',m,m,s,LDST,scpy,LDST)
                        CALL DLACPY('F',m,m,t,LDST,tcpy,LDST)
                        CALL DLACPY('F',m,m,ir,LDST,ircop,LDST)
                        CALL DLACPY('F',m,m,li,LDST,licop,LDST)
!
!        Triangularize the B-part by an RQ factorization.
!        Apply transformation (from left) to A-part, giving S.
!
                        CALL DGERQ2(m,m,t,LDST,taur,Work,linfo)
                        IF ( linfo/=0 ) THEN
                           Info = 1
                        ELSE
                           CALL DORMR2('R','T',m,m,m,t,LDST,taur,s,LDST,&
     &                                 Work,linfo)
                           IF ( linfo/=0 ) THEN
                              Info = 1
                           ELSE
                              CALL DORMR2('L','N',m,m,m,t,LDST,taur,ir, &
     &                           LDST,Work,linfo)
                              IF ( linfo/=0 ) THEN
                                 Info = 1
                              ELSE
!
!        Compute F-norm(S21) in BRQA21. (T21 is 0.)
!
                                 dscale = ZERO
                                 dsum = ONE
                                 DO i = 1 , N2
                                    CALL DLASSQ(N1,s(N2+1,i),1,dscale,  &
     &                                 dsum)
                                 ENDDO
                                 brqa21 = dscale*SQRT(dsum)
!
!        Triangularize the B-part by a QR factorization.
!        Apply transformation (from right) to A-part, giving S.
!
                                 CALL DGEQR2(m,m,tcpy,LDST,taul,Work,   &
     &                              linfo)
                                 IF ( linfo/=0 ) THEN
                                    Info = 1
                                 ELSE
                                    CALL DORM2R('L','T',m,m,m,tcpy,LDST,&
     &                                 taul,scpy,LDST,Work,Info)
                                    CALL DORM2R('R','N',m,m,m,tcpy,LDST,&
     &                                 taul,licop,LDST,Work,Info)
                                    IF ( linfo/=0 ) THEN
                                       Info = 1
                                    ELSE
!
!        Compute F-norm(S21) in BQRA21. (T21 is 0.)
!
                                       dscale = ZERO
                                       dsum = ONE
                                       DO i = 1 , N2
                                         CALL DLASSQ(N1,scpy(N2+1,i),1, &
     &                                      dscale,dsum)
                                       ENDDO
                                       bqra21 = dscale*SQRT(dsum)
!
!        Decide which method to use.
!          Weak stability test:
!             F-norm(S21) <= O(EPS * F-norm((S)))
!
                                       IF ( bqra21<=brqa21 .AND.        &
     &                                    bqra21<=thresha ) THEN
                                         CALL DLACPY('F',m,m,scpy,LDST, &
     &                                      s,LDST)
                                         CALL DLACPY('F',m,m,tcpy,LDST, &
     &                                      t,LDST)
                                         CALL DLACPY('F',m,m,ircop,LDST,&
     &                                      ir,LDST)
                                         CALL DLACPY('F',m,m,licop,LDST,&
     &                                      li,LDST)
                                       ELSEIF ( brqa21>=thresha ) THEN
                                         Info = 1
                                         GOTO 99999
                                       ENDIF
!
!        Set lower triangle of B-part to zero
!
                                       CALL DLASET('Lower',m-1,m-1,ZERO,&
     &                                    ZERO,t(2,1),LDST)
!
                                       IF ( WANDS ) THEN
!
!           Strong stability test:
!               F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
!               and
!               F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))
!
                                         CALL DLACPY('Full',m,m,A(J1,J1)&
     &                                      ,Lda,Work(m*m+1),m)
                                         CALL DGEMM('N','N',m,m,m,ONE,  &
     &                                      li,LDST,s,LDST,ZERO,Work,m)
                                         CALL DGEMM('N','N',m,m,m,-ONE, &
     &                                      Work,m,ir,LDST,ONE,         &
     &                                      Work(m*m+1),m)
                                         dscale = ZERO
                                         dsum = ONE
                                         CALL DLASSQ(m*m,Work(m*m+1),1, &
     &                                      dscale,dsum)
                                         sa = dscale*SQRT(dsum)
!
                                         CALL DLACPY('Full',m,m,B(J1,J1)&
     &                                      ,Ldb,Work(m*m+1),m)
                                         CALL DGEMM('N','N',m,m,m,ONE,  &
     &                                      li,LDST,t,LDST,ZERO,Work,m)
                                         CALL DGEMM('N','N',m,m,m,-ONE, &
     &                                      Work,m,ir,LDST,ONE,         &
     &                                      Work(m*m+1),m)
                                         dscale = ZERO
                                         dsum = ONE
                                         CALL DLASSQ(m*m,Work(m*m+1),1, &
     &                                      dscale,dsum)
                                         sb = dscale*SQRT(dsum)
                                         strong = sa<=thresha .AND.     &
     &                                      sb<=threshb
                                         IF ( .NOT.strong ) THEN
                                         Info = 1
                                         GOTO 99999
                                         ENDIF
!
                                       ENDIF
!
!        If the swap is accepted ("weakly" and "strongly"), apply the
!        transformations and set N1-by-N2 (2,1)-block to zero.
!
                                       CALL DLASET('Full',N1,N2,ZERO,   &
     &                                    ZERO,s(N2+1,1),LDST)
!
!        copy back M-by-M diagonal block starting at index J1 of (A, B)
!
                                       CALL DLACPY('F',m,m,s,LDST,      &
     &                                    A(J1,J1),Lda)
                                       CALL DLACPY('F',m,m,t,LDST,      &
     &                                    B(J1,J1),Ldb)
                                       CALL DLASET('Full',LDST,LDST,    &
     &                                    ZERO,ZERO,t,LDST)
!
!        Standardize existing 2-by-2 blocks.
!
                                       CALL DLASET('Full',m,m,ZERO,ZERO,&
     &                                    Work,m)
                                       Work(1) = ONE
                                       t(1,1) = ONE
                                       idum = Lwork - m*m - 2
                                       IF ( N2>1 ) THEN
                                         CALL DLAGV2(A(J1,J1),Lda,      &
     &                                      B(J1,J1),Ldb,ar,ai,be,      &
     &                                      Work(1),Work(2),t(1,1),     &
     &                                      t(2,1))
                                         Work(m+1) = -Work(2)
                                         Work(m+2) = Work(1)
                                         t(N2,N2) = t(1,1)
                                         t(1,2) = -t(2,1)
                                       ENDIF
                                       Work(m*m) = ONE
                                       t(m,m) = ONE
!
                                       IF ( N1>1 ) THEN
                                         CALL DLAGV2(A(J1+N2,J1+N2),Lda,&
     &                                      B(J1+N2,J1+N2),Ldb,taur,    &
     &                                      taul,Work(m*m+1),           &
     &                                      Work(N2*m+N2+1),            &
     &                                      Work(N2*m+N2+2),t(N2+1,N2+1)&
     &                                      ,t(m,m-1))
                                         Work(m*m) = Work(N2*m+N2+1)
                                         Work(m*m-1) = -Work(N2*m+N2+2)
                                         t(m,m) = t(N2+1,N2+1)
                                         t(m-1,m) = -t(m,m-1)
                                       ENDIF
                                       CALL DGEMM('T','N',N2,N1,N2,ONE, &
     &                                    Work,m,A(J1,J1+N2),Lda,ZERO,  &
     &                                    Work(m*m+1),N2)
                                       CALL DLACPY('Full',N2,N1,        &
     &                                    Work(m*m+1),N2,A(J1,J1+N2),   &
     &                                    Lda)
                                       CALL DGEMM('T','N',N2,N1,N2,ONE, &
     &                                    Work,m,B(J1,J1+N2),Ldb,ZERO,  &
     &                                    Work(m*m+1),N2)
                                       CALL DLACPY('Full',N2,N1,        &
     &                                    Work(m*m+1),N2,B(J1,J1+N2),   &
     &                                    Ldb)
                                       CALL DGEMM('N','N',m,m,m,ONE,li, &
     &                                    LDST,Work,m,ZERO,Work(m*m+1), &
     &                                    m)
                                       CALL DLACPY('Full',m,m,          &
     &                                    Work(m*m+1),m,li,LDST)
                                       CALL DGEMM('N','N',N2,N1,N1,ONE, &
     &                                    A(J1,J1+N2),Lda,t(N2+1,N2+1), &
     &                                    LDST,ZERO,Work,N2)
                                       CALL DLACPY('Full',N2,N1,Work,N2,&
     &                                    A(J1,J1+N2),Lda)
                                       CALL DGEMM('N','N',N2,N1,N1,ONE, &
     &                                    B(J1,J1+N2),Ldb,t(N2+1,N2+1), &
     &                                    LDST,ZERO,Work,N2)
                                       CALL DLACPY('Full',N2,N1,Work,N2,&
     &                                    B(J1,J1+N2),Ldb)
                                       CALL DGEMM('T','N',m,m,m,ONE,ir, &
     &                                    LDST,t,LDST,ZERO,Work,m)
                                       CALL DLACPY('Full',m,m,Work,m,ir,&
     &                                    LDST)
!
!        Accumulate transformations into Q and Z if requested.
!
                                       IF ( Wantq ) THEN
                                         CALL DGEMM('N','N',N,m,m,ONE,  &
     &                                      Q(1,J1),Ldq,li,LDST,ZERO,   &
     &                                      Work,N)
                                         CALL DLACPY('Full',N,m,Work,N, &
     &                                      Q(1,J1),Ldq)
!
                                       ENDIF
!
                                       IF ( Wantz ) THEN
                                         CALL DGEMM('N','N',N,m,m,ONE,  &
     &                                      Z(1,J1),Ldz,ir,LDST,ZERO,   &
     &                                      Work,N)
                                         CALL DLACPY('Full',N,m,Work,N, &
     &                                      Z(1,J1),Ldz)
!
                                       ENDIF
!
!        Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and
!                (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)).
!
                                       i = J1 + m
                                       IF ( i<=N ) THEN
                                         CALL DGEMM('T','N',m,N-i+1,m,  &
     &                                      ONE,li,LDST,A(J1,i),Lda,    &
     &                                      ZERO,Work,m)
                                         CALL DLACPY('Full',m,N-i+1,    &
     &                                      Work,m,A(J1,i),Lda)
                                         CALL DGEMM('T','N',m,N-i+1,m,  &
     &                                      ONE,li,LDST,B(J1,i),Ldb,    &
     &                                      ZERO,Work,m)
                                         CALL DLACPY('Full',m,N-i+1,    &
     &                                      Work,m,B(J1,i),Ldb)
                                       ENDIF
                                       i = J1 - 1
                                       IF ( i>0 ) THEN
                                         CALL DGEMM('N','N',i,m,m,ONE,  &
     &                                      A(1,J1),Lda,ir,LDST,ZERO,   &
     &                                      Work,i)
                                         CALL DLACPY('Full',i,m,Work,i, &
     &                                      A(1,J1),Lda)
                                         CALL DGEMM('N','N',i,m,m,ONE,  &
     &                                      B(1,J1),Ldb,ir,LDST,ZERO,   &
     &                                      Work,i)
                                         CALL DLACPY('Full',i,m,Work,i, &
     &                                      B(1,J1),Ldb)
                                       ENDIF
!
!        Exit with INFO = 0 if swap was successfully performed.
!
                                       RETURN
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
!
      ENDIF
!
!     End of DTGEX2
!
99999 END SUBROUTINE DTGEX2
