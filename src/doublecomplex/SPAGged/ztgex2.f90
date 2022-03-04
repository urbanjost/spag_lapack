!*==ztgex2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZTGEX2 swaps adjacent diagonal blocks in an upper (quasi) triangular matrix pair by an unitary equivalence transformation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTGEX2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgex2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgex2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgex2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
!                          LDZ, J1, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            WANTQ, WANTZ
!       INTEGER            INFO, J1, LDA, LDB, LDQ, LDZ, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTGEX2 swaps adjacent diagonal 1 by 1 blocks (A11,B11) and (A22,B22)
!> in an upper triangular matrix pair (A, B) by an unitary equivalence
!> transformation.
!>
!> (A, B) must be in generalized Schur canonical form, that is, A and
!> B are both upper triangular.
!>
!> Optionally, the matrices Q and Z of generalized Schur vectors are
!> updated.
!>
!>        Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
!>        Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H
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
!>          A is COMPLEX*16 array, dimensions (LDA,N)
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
!>          B is COMPLEX*16 array, dimensions (LDB,N)
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
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>          If WANTQ = .TRUE, on entry, the unitary matrix Q. On exit,
!>          the updated matrix Q.
!>          Not referenced if WANTQ = .FALSE..
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q. LDQ >= 1;
!>          If WANTQ = .TRUE., LDQ >= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ,N)
!>          If WANTZ = .TRUE, on entry, the unitary matrix Z. On exit,
!>          the updated matrix Z.
!>          Not referenced if WANTZ = .FALSE..
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= 1;
!>          If WANTZ = .TRUE., LDZ >= N.
!> \endverbatim
!>
!> \param[in] J1
!> \verbatim
!>          J1 is INTEGER
!>          The index to the first block (A11, B11).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           =0:  Successful exit.
!>           =1:  The transformed matrix pair (A, B) would be too far
!>                from generalized Schur form; the problem is ill-
!>                conditioned.
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
!> \ingroup complex16GEauxiliary
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
!>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
!>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
!>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
!>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
!> \n
!>  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
!>      Eigenvalues of a Regular Matrix Pair (A, B) and Condition
!>      Estimation: Theory, Algorithms and Software, Report UMINF-94.04,
!>      Department of Computing Science, Umea University, S-901 87 Umea,
!>      Sweden, 1994. Also as LAPACK Working Note 87. To appear in
!>      Numerical Algorithms, 1996.
!>
!  =====================================================================
      SUBROUTINE ZTGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,J1,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_ZLACPY
      USE S_ZLARTG
      USE S_ZLASSQ
      USE S_ZROT
      IMPLICIT NONE
!*--ZTGEX2199
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  TWENTY = 2.0D+1
      INTEGER , PARAMETER  ::  LDST = 2
      LOGICAL , PARAMETER  ::  WANDS = .TRUE.
!
! Dummy argument declarations rewritten by SPAG
!
      LOGICAL , INTENT(IN) :: Wantq
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER , INTENT(IN) :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      INTEGER , INTENT(IN) :: J1
      INTEGER , INTENT(OUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) :: cdum , f , g , sq , sz
      REAL(R8KIND) :: cq , cz , eps , sa , sb , scale , smlnum , sum ,  &
     &                thresha , threshb
      INTEGER :: i , m
      COMPLEX(CX16KIND) , DIMENSION(LDST,LDST) :: s , t
      LOGICAL :: strong , weak
      COMPLEX(CX16KIND) , DIMENSION(8) :: work
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
      Info = 0
!
!     Quick return if possible
!
      IF ( N<=1 ) RETURN
!
      m = LDST
      weak = .FALSE.
      strong = .FALSE.
!
!     Make a local copy of selected block in (A, B)
!
      CALL ZLACPY('Full',m,m,A(J1,J1),Lda,s,LDST)
      CALL ZLACPY('Full',m,m,B(J1,J1),Ldb,t,LDST)
!
!     Compute the threshold for testing the acceptance of swapping.
!
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')/eps
      scale = DBLE(CZERO)
      sum = DBLE(CONE)
      CALL ZLACPY('Full',m,m,s,LDST,work,m)
      CALL ZLACPY('Full',m,m,t,LDST,work(m*m+1),m)
      CALL ZLASSQ(m*m,work,1,scale,sum)
      sa = scale*SQRT(sum)
      scale = DBLE(CZERO)
      sum = DBLE(CONE)
      CALL ZLASSQ(m*m,work(m*m+1),1,scale,sum)
      sb = scale*SQRT(sum)
!
!     THRES has been changed from
!        THRESH = MAX( TEN*EPS*SA, SMLNUM )
!     to
!        THRESH = MAX( TWENTY*EPS*SA, SMLNUM )
!     on 04/01/10.
!     "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by
!     Jim Demmel and Guillaume Revy. See forum post 1783.
!
      thresha = MAX(TWENTY*eps*sa,smlnum)
      threshb = MAX(TWENTY*eps*sb,smlnum)
!
!     Compute unitary QL and RQ that swap 1-by-1 and 1-by-1 blocks
!     using Givens rotations and perform the swap tentatively.
!
      f = s(2,2)*t(1,1) - t(2,2)*s(1,1)
      g = s(2,2)*t(1,2) - t(2,2)*s(1,2)
      sa = ABS(s(2,2))*ABS(t(1,1))
      sb = ABS(s(1,1))*ABS(t(2,2))
      CALL ZLARTG(g,f,cz,sz,cdum)
      sz = -sz
      CALL ZROT(2,s(1,1),1,s(1,2),1,cz,DCONJG(sz))
      CALL ZROT(2,t(1,1),1,t(1,2),1,cz,DCONJG(sz))
      IF ( sa>=sb ) THEN
         CALL ZLARTG(s(1,1),s(2,1),cq,sq,cdum)
      ELSE
         CALL ZLARTG(t(1,1),t(2,1),cq,sq,cdum)
      ENDIF
      CALL ZROT(2,s(1,1),LDST,s(2,1),LDST,cq,sq)
      CALL ZROT(2,t(1,1),LDST,t(2,1),LDST,cq,sq)
!
!     Weak stability test: |S21| <= O(EPS F-norm((A)))
!                          and  |T21| <= O(EPS F-norm((B)))
!
      weak = ABS(s(2,1))<=thresha .AND. ABS(t(2,1))<=threshb
      IF ( .NOT.weak ) THEN
!
!     Exit with INFO = 1 if swap was rejected.
!
         Info = 1
      ELSE
!
         IF ( WANDS ) THEN
!
!        Strong stability test:
!           F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
!           and
!           F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))
!
            CALL ZLACPY('Full',m,m,s,LDST,work,m)
            CALL ZLACPY('Full',m,m,t,LDST,work(m*m+1),m)
            CALL ZROT(2,work,1,work(3),1,cz,-DCONJG(sz))
            CALL ZROT(2,work(5),1,work(7),1,cz,-DCONJG(sz))
            CALL ZROT(2,work,2,work(2),2,cq,-sq)
            CALL ZROT(2,work(5),2,work(6),2,cq,-sq)
            DO i = 1 , 2
               work(i) = work(i) - A(J1+i-1,J1)
               work(i+2) = work(i+2) - A(J1+i-1,J1+1)
               work(i+4) = work(i+4) - B(J1+i-1,J1)
               work(i+6) = work(i+6) - B(J1+i-1,J1+1)
            ENDDO
            scale = DBLE(CZERO)
            sum = DBLE(CONE)
            CALL ZLASSQ(m*m,work,1,scale,sum)
            sa = scale*SQRT(sum)
            scale = DBLE(CZERO)
            sum = DBLE(CONE)
            CALL ZLASSQ(m*m,work(m*m+1),1,scale,sum)
            sb = scale*SQRT(sum)
            strong = sa<=thresha .AND. sb<=threshb
            IF ( .NOT.strong ) THEN
               Info = 1
               GOTO 99999
            ENDIF
         ENDIF
!
!     If the swap is accepted ("weakly" and "strongly"), apply the
!     equivalence transformations to the original matrix pair (A,B)
!
         CALL ZROT(J1+1,A(1,J1),1,A(1,J1+1),1,cz,DCONJG(sz))
         CALL ZROT(J1+1,B(1,J1),1,B(1,J1+1),1,cz,DCONJG(sz))
         CALL ZROT(N-J1+1,A(J1,J1),Lda,A(J1+1,J1),Lda,cq,sq)
         CALL ZROT(N-J1+1,B(J1,J1),Ldb,B(J1+1,J1),Ldb,cq,sq)
!
!     Set  N1 by N2 (2,1) blocks to 0
!
         A(J1+1,J1) = CZERO
         B(J1+1,J1) = CZERO
!
!     Accumulate transformations into Q and Z if requested.
!
         IF ( Wantz ) CALL ZROT(N,Z(1,J1),1,Z(1,J1+1),1,cz,DCONJG(sz))
         IF ( Wantq ) CALL ZROT(N,Q(1,J1),1,Q(1,J1+1),1,cq,DCONJG(sq))
!
!     Exit with INFO = 0 if swap was successfully performed.
!
         RETURN
      ENDIF
!
!     End of ZTGEX2
!
99999 END SUBROUTINE ZTGEX2
