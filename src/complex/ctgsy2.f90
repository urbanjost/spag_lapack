!*==ctgsy2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CTGSY2 solves the generalized Sylvester equation (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTGSY2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctgsy2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctgsy2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctgsy2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,
!                          LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N
!       REAL               RDSCAL, RDSUM, SCALE
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * ),
!      $                   D( LDD, * ), E( LDE, * ), F( LDF, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTGSY2 solves the generalized Sylvester equation
!>
!>             A * R - L * B = scale *  C               (1)
!>             D * R - L * E = scale * F
!>
!> using Level 1 and 2 BLAS, where R and L are unknown M-by-N matrices,
!> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
!> N-by-N and M-by-N, respectively. A, B, D and E are upper triangular
!> (i.e., (A,D) and (B,E) in generalized Schur form).
!>
!> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
!> scaling factor chosen to avoid overflow.
!>
!> In matrix notation solving equation (1) corresponds to solve
!> Zx = scale * b, where Z is defined as
!>
!>        Z = [ kron(In, A)  -kron(B**H, Im) ]             (2)
!>            [ kron(In, D)  -kron(E**H, Im) ],
!>
!> Ik is the identity matrix of size k and X**H is the transpose of X.
!> kron(X, Y) is the Kronecker product between the matrices X and Y.
!>
!> If TRANS = 'C', y in the conjugate transposed system Z**H*y = scale*b
!> is solved for, which is equivalent to solve for R and L in
!>
!>             A**H * R  + D**H * L   = scale * C           (3)
!>             R  * B**H + L  * E**H  = scale * -F
!>
!> This case is used to compute an estimate of Dif[(A, D), (B, E)] =
!> = sigma_min(Z) using reverse communication with CLACON.
!>
!> CTGSY2 also (IJOB >= 1) contributes to the computation in CTGSYL
!> of an upper bound on the separation between to matrix pairs. Then
!> the input (A, D), (B, E) are sub-pencils of two matrix pairs in
!> CTGSYL.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': solve the generalized Sylvester equation (1).
!>          = 'T': solve the 'transposed' system (3).
!> \endverbatim
!>
!> \param[in] IJOB
!> \verbatim
!>          IJOB is INTEGER
!>          Specifies what kind of functionality to be performed.
!>          = 0: solve (1) only.
!>          = 1: A contribution from this subsystem to a Frobenius
!>               norm-based estimate of the separation between two matrix
!>               pairs is computed. (look ahead strategy is used).
!>          = 2: A contribution from this subsystem to a Frobenius
!>               norm-based estimate of the separation between two matrix
!>               pairs is computed. (SGECON on sub-systems is used.)
!>          Not referenced if TRANS = 'T'.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          On entry, M specifies the order of A and D, and the row
!>          dimension of C, F, R and L.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          On entry, N specifies the order of B and E, and the column
!>          dimension of C, F, R and L.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, M)
!>          On entry, A contains an upper triangular matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the matrix A. LDA >= max(1, M).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB, N)
!>          On entry, B contains an upper triangular matrix.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the matrix B. LDB >= max(1, N).
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC, N)
!>          On entry, C contains the right-hand-side of the first matrix
!>          equation in (1).
!>          On exit, if IJOB = 0, C has been overwritten by the solution
!>          R.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the matrix C. LDC >= max(1, M).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX array, dimension (LDD, M)
!>          On entry, D contains an upper triangular matrix.
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of the matrix D. LDD >= max(1, M).
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX array, dimension (LDE, N)
!>          On entry, E contains an upper triangular matrix.
!> \endverbatim
!>
!> \param[in] LDE
!> \verbatim
!>          LDE is INTEGER
!>          The leading dimension of the matrix E. LDE >= max(1, N).
!> \endverbatim
!>
!> \param[in,out] F
!> \verbatim
!>          F is COMPLEX array, dimension (LDF, N)
!>          On entry, F contains the right-hand-side of the second matrix
!>          equation in (1).
!>          On exit, if IJOB = 0, F has been overwritten by the solution
!>          L.
!> \endverbatim
!>
!> \param[in] LDF
!> \verbatim
!>          LDF is INTEGER
!>          The leading dimension of the matrix F. LDF >= max(1, M).
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is REAL
!>          On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the solutions
!>          R and L (C and F on entry) will hold the solutions to a
!>          slightly perturbed system but the input matrices A, B, D and
!>          E have not been changed. If SCALE = 0, R and L will hold the
!>          solutions to the homogeneous system with C = F = 0.
!>          Normally, SCALE = 1.
!> \endverbatim
!>
!> \param[in,out] RDSUM
!> \verbatim
!>          RDSUM is REAL
!>          On entry, the sum of squares of computed contributions to
!>          the Dif-estimate under computation by CTGSYL, where the
!>          scaling factor RDSCAL (see below) has been factored out.
!>          On exit, the corresponding sum of squares updated with the
!>          contributions from the current sub-system.
!>          If TRANS = 'T' RDSUM is not touched.
!>          NOTE: RDSUM only makes sense when CTGSY2 is called by
!>          CTGSYL.
!> \endverbatim
!>
!> \param[in,out] RDSCAL
!> \verbatim
!>          RDSCAL is REAL
!>          On entry, scaling factor used to prevent overflow in RDSUM.
!>          On exit, RDSCAL is updated w.r.t. the current contributions
!>          in RDSUM.
!>          If TRANS = 'T', RDSCAL is not touched.
!>          NOTE: RDSCAL only makes sense when CTGSY2 is called by
!>          CTGSYL.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          On exit, if INFO is set to
!>            =0: Successful exit
!>            <0: If INFO = -i, input argument number i is illegal.
!>            >0: The matrix pairs (A, D) and (B, E) have common or very
!>                close eigenvalues.
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
!> \ingroup complexSYauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!  =====================================================================
      SUBROUTINE CTGSY2(Trans,Ijob,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F, &
     &                  Ldf,Scale,Rdsum,Rdscal,Info)
      IMPLICIT NONE
!*--CTGSY2262
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Ijob , Info , Lda , Ldb , Ldc , Ldd , Lde , Ldf , M , N
      REAL Rdscal , Rdsum , Scale
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , B(Ldb,*) , C(Ldc,*) , D(Ldd,*) , E(Lde,*) ,    &
     &        F(Ldf,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      INTEGER LDZ
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0,LDZ=2)
!     ..
!     .. Local Scalars ..
      LOGICAL notran
      INTEGER i , ierr , j , k
      REAL scaloc
      COMPLEX alpha
!     ..
!     .. Local Arrays ..
      INTEGER ipiv(LDZ) , jpiv(LDZ)
      COMPLEX rhs(LDZ) , z(LDZ,LDZ)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY , CGESC2 , CGETC2 , CSCAL , CLATDF , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , CONJG , MAX
!     ..
!     .. Executable Statements ..
!
!     Decode and test input parameters
!
      Info = 0
      ierr = 0
      notran = LSAME(Trans,'N')
      IF ( .NOT.notran .AND. .NOT.LSAME(Trans,'C') ) THEN
         Info = -1
      ELSEIF ( notran ) THEN
         IF ( (Ijob<0) .OR. (Ijob>2) ) Info = -2
      ENDIF
      IF ( Info==0 ) THEN
         IF ( M<=0 ) THEN
            Info = -3
         ELSEIF ( N<=0 ) THEN
            Info = -4
         ELSEIF ( Lda<MAX(1,M) ) THEN
            Info = -6
         ELSEIF ( Ldb<MAX(1,N) ) THEN
            Info = -8
         ELSEIF ( Ldc<MAX(1,M) ) THEN
            Info = -10
         ELSEIF ( Ldd<MAX(1,M) ) THEN
            Info = -12
         ELSEIF ( Lde<MAX(1,N) ) THEN
            Info = -14
         ELSEIF ( Ldf<MAX(1,M) ) THEN
            Info = -16
         ENDIF
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CTGSY2',-Info)
         RETURN
      ENDIF
!
      IF ( notran ) THEN
!
!        Solve (I, J) - system
!           A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
!           D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
!        for I = M, M - 1, ..., 1; J = 1, 2, ..., N
!
         Scale = ONE
         scaloc = ONE
         DO j = 1 , N
            DO i = M , 1 , -1
!
!              Build 2 by 2 system
!
               z(1,1) = A(i,i)
               z(2,1) = D(i,i)
               z(1,2) = -B(j,j)
               z(2,2) = -E(j,j)
!
!              Set up right hand side(s)
!
               rhs(1) = C(i,j)
               rhs(2) = F(i,j)
!
!              Solve Z * x = RHS
!
               CALL CGETC2(LDZ,z,LDZ,ipiv,jpiv,ierr)
               IF ( ierr>0 ) Info = ierr
               IF ( Ijob==0 ) THEN
                  CALL CGESC2(LDZ,z,LDZ,rhs,ipiv,jpiv,scaloc)
                  IF ( scaloc/=ONE ) THEN
                     DO k = 1 , N
                        CALL CSCAL(M,CMPLX(scaloc,ZERO),C(1,k),1)
                        CALL CSCAL(M,CMPLX(scaloc,ZERO),F(1,k),1)
                     ENDDO
                     Scale = Scale*scaloc
                  ENDIF
               ELSE
                  CALL CLATDF(Ijob,LDZ,z,LDZ,rhs,Rdsum,Rdscal,ipiv,jpiv)
               ENDIF
!
!              Unpack solution vector(s)
!
               C(i,j) = rhs(1)
               F(i,j) = rhs(2)
!
!              Substitute R(I, J) and L(I, J) into remaining equation.
!
               IF ( i>1 ) THEN
                  alpha = -rhs(1)
                  CALL CAXPY(i-1,alpha,A(1,i),1,C(1,j),1)
                  CALL CAXPY(i-1,alpha,D(1,i),1,F(1,j),1)
               ENDIF
               IF ( j<N ) THEN
                  CALL CAXPY(N-j,rhs(2),B(j,j+1),Ldb,C(i,j+1),Ldc)
                  CALL CAXPY(N-j,rhs(2),E(j,j+1),Lde,F(i,j+1),Ldf)
               ENDIF
!
            ENDDO
         ENDDO
      ELSE
!
!        Solve transposed (I, J) - system:
!           A(I, I)**H * R(I, J) + D(I, I)**H * L(J, J) = C(I, J)
!           R(I, I) * B(J, J) + L(I, J) * E(J, J)   = -F(I, J)
!        for I = 1, 2, ..., M, J = N, N - 1, ..., 1
!
         Scale = ONE
         scaloc = ONE
         DO i = 1 , M
            DO j = N , 1 , -1
!
!              Build 2 by 2 system Z**H
!
               z(1,1) = CONJG(A(i,i))
               z(2,1) = -CONJG(B(j,j))
               z(1,2) = CONJG(D(i,i))
               z(2,2) = -CONJG(E(j,j))
!
!
!              Set up right hand side(s)
!
               rhs(1) = C(i,j)
               rhs(2) = F(i,j)
!
!              Solve Z**H * x = RHS
!
               CALL CGETC2(LDZ,z,LDZ,ipiv,jpiv,ierr)
               IF ( ierr>0 ) Info = ierr
               CALL CGESC2(LDZ,z,LDZ,rhs,ipiv,jpiv,scaloc)
               IF ( scaloc/=ONE ) THEN
                  DO k = 1 , N
                     CALL CSCAL(M,CMPLX(scaloc,ZERO),C(1,k),1)
                     CALL CSCAL(M,CMPLX(scaloc,ZERO),F(1,k),1)
                  ENDDO
                  Scale = Scale*scaloc
               ENDIF
!
!              Unpack solution vector(s)
!
               C(i,j) = rhs(1)
               F(i,j) = rhs(2)
!
!              Substitute R(I, J) and L(I, J) into remaining equation.
!
               DO k = 1 , j - 1
                  F(i,k) = F(i,k) + rhs(1)*CONJG(B(k,j)) + rhs(2)       &
     &                     *CONJG(E(k,j))
               ENDDO
               DO k = i + 1 , M
                  C(k,j) = C(k,j) - CONJG(A(i,k))*rhs(1) - CONJG(D(i,k))&
     &                     *rhs(2)
               ENDDO
!
            ENDDO
         ENDDO
      ENDIF
!
!     End of CTGSY2
!
      END SUBROUTINE CTGSY2
