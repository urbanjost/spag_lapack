!*==ctgexc.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CTGEXC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTGEXC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctgexc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctgexc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctgexc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
!                          LDZ, IFST, ILST, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            WANTQ, WANTZ
!       INTEGER            IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTGEXC reorders the generalized Schur decomposition of a complex
!> matrix pair (A,B), using an unitary equivalence transformation
!> (A, B) := Q * (A, B) * Z**H, so that the diagonal block of (A, B) with
!> row index IFST is moved to row ILST.
!>
!> (A, B) must be in generalized Schur canonical form, that is, A and
!> B are both upper triangular.
!>
!> Optionally, the matrices Q and Z of generalized Schur vectors are
!> updated.
!>
!>        Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
!>        Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the upper triangular matrix A in the pair (A, B).
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
!>          B is COMPLEX array, dimension (LDB,N)
!>          On entry, the upper triangular matrix B in the pair (A, B).
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
!>          Q is COMPLEX array, dimension (LDQ,N)
!>          On entry, if WANTQ = .TRUE., the unitary matrix Q.
!>          On exit, the updated matrix Q.
!>          If WANTQ = .FALSE., Q is not referenced.
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
!>          Z is COMPLEX array, dimension (LDZ,N)
!>          On entry, if WANTZ = .TRUE., the unitary matrix Z.
!>          On exit, the updated matrix Z.
!>          If WANTZ = .FALSE., Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= 1;
!>          If WANTZ = .TRUE., LDZ >= N.
!> \endverbatim
!>
!> \param[in] IFST
!> \verbatim
!>          IFST is INTEGER
!> \endverbatim
!>
!> \param[in,out] ILST
!> \verbatim
!>          ILST is INTEGER
!>          Specify the reordering of the diagonal blocks of (A, B).
!>          The block with row index IFST is moved to row ILST, by a
!>          sequence of swapping between adjacent blocks.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           =0:  Successful exit.
!>           <0:  if INFO = -i, the i-th argument had an illegal value.
!>           =1:  The transformed matrix pair (A, B) would be too far
!>                from generalized Schur form; the problem is ill-
!>                conditioned. (A, B) may have been partially reordered,
!>                and ILST points to the first row of the current
!>                position of the block being moved.
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
!> \ingroup complexGEcomputational
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
!>      Estimation: Theory, Algorithms and Software, Report
!>      UMINF - 94.04, Department of Computing Science, Umea University,
!>      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.
!>      To appear in Numerical Algorithms, 1996.
!> \n
!>  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
!>      for Solving the Generalized Sylvester Equation and Estimating the
!>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,
!>      Department of Computing Science, Umea University, S-901 87 Umea,
!>      Sweden, December 1993, Revised April 1994, Also as LAPACK working
!>      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,
!>      1996.
!>
!  =====================================================================
      SUBROUTINE CTGEXC(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,Ifst,Ilst,&
     &                  Info)
      USE S_CTGEX2
      USE S_XERBLA
      IMPLICIT NONE
!*--CTGEXC206
!
! Dummy argument declarations rewritten by SPAG
!
      LOGICAL :: Wantq
      LOGICAL :: Wantz
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(IN) :: Ifst
      INTEGER , INTENT(INOUT) :: Ilst
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: here
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Decode and test input arguments.
      Info = 0
      IF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldq<1 .OR. Wantq .AND. (Ldq<MAX(1,N)) ) THEN
         Info = -9
      ELSEIF ( Ldz<1 .OR. Wantz .AND. (Ldz<MAX(1,N)) ) THEN
         Info = -11
      ELSEIF ( Ifst<1 .OR. Ifst>N ) THEN
         Info = -12
      ELSEIF ( Ilst<1 .OR. Ilst>N ) THEN
         Info = -13
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CTGEXC',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=1 ) RETURN
      IF ( Ifst==Ilst ) RETURN
!
      IF ( Ifst<Ilst ) THEN
!
         here = Ifst
         DO
!
!
!        Swap with next one below
!
            CALL CTGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,here,Info)
            IF ( Info/=0 ) THEN
               Ilst = here
               RETURN
            ENDIF
            here = here + 1
            IF ( here>=Ilst ) THEN
               here = here - 1
               EXIT
            ENDIF
         ENDDO
      ELSE
         here = Ifst - 1
         DO
!
!
!        Swap with next one above
!
            CALL CTGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,here,Info)
            IF ( Info/=0 ) THEN
               Ilst = here
               RETURN
            ENDIF
            here = here - 1
            IF ( here<Ilst ) THEN
               here = here + 1
               EXIT
            ENDIF
         ENDDO
      ENDIF
      Ilst = here
!
!     End of CTGEXC
!
      END SUBROUTINE CTGEXC
