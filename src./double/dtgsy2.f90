!*==dtgsy2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DTGSY2 solves the generalized Sylvester equation (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DTGSY2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgsy2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgsy2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgsy2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,
!                          LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL,
!                          IWORK, PQ, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N,
!      $                   PQ
!       DOUBLE PRECISION   RDSCAL, RDSUM, SCALE
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
!      $                   D( LDD, * ), E( LDE, * ), F( LDF, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTGSY2 solves the generalized Sylvester equation:
!>
!>             A * R - L * B = scale * C                (1)
!>             D * R - L * E = scale * F,
!>
!> using Level 1 and 2 BLAS. where R and L are unknown M-by-N matrices,
!> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
!> N-by-N and M-by-N, respectively, with real entries. (A, D) and (B, E)
!> must be in generalized Schur canonical form, i.e. A, B are upper
!> quasi triangular and D, E are upper triangular. The solution (R, L)
!> overwrites (C, F). 0 <= SCALE <= 1 is an output scaling factor
!> chosen to avoid overflow.
!>
!> In matrix notation solving equation (1) corresponds to solve
!> Z*x = scale*b, where Z is defined as
!>
!>        Z = [ kron(In, A)  -kron(B**T, Im) ]             (2)
!>            [ kron(In, D)  -kron(E**T, Im) ],
!>
!> Ik is the identity matrix of size k and X**T is the transpose of X.
!> kron(X, Y) is the Kronecker product between the matrices X and Y.
!> In the process of solving (1), we solve a number of such systems
!> where Dim(In), Dim(In) = 1 or 2.
!>
!> If TRANS = 'T', solve the transposed system Z**T*y = scale*b for y,
!> which is equivalent to solve for R and L in
!>
!>             A**T * R  + D**T * L   = scale * C           (3)
!>             R  * B**T + L  * E**T  = scale * -F
!>
!> This case is used to compute an estimate of Dif[(A, D), (B, E)] =
!> sigma_min(Z) using reverse communication with DLACON.
!>
!> DTGSY2 also (IJOB >= 1) contributes to the computation in DTGSYL
!> of an upper bound on the separation between to matrix pairs. Then
!> the input (A, D), (B, E) are sub-pencils of the matrix pair in
!> DTGSYL. See DTGSYL for details.
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
!>               pairs is computed. (DGECON on sub-systems is used.)
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
!>          A is DOUBLE PRECISION array, dimension (LDA, M)
!>          On entry, A contains an upper quasi triangular matrix.
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
!>          B is DOUBLE PRECISION array, dimension (LDB, N)
!>          On entry, B contains an upper quasi triangular matrix.
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
!>          C is DOUBLE PRECISION array, dimension (LDC, N)
!>          On entry, C contains the right-hand-side of the first matrix
!>          equation in (1).
!>          On exit, if IJOB = 0, C has been overwritten by the
!>          solution R.
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
!>          D is DOUBLE PRECISION array, dimension (LDD, M)
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
!>          E is DOUBLE PRECISION array, dimension (LDE, N)
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
!>          F is DOUBLE PRECISION array, dimension (LDF, N)
!>          On entry, F contains the right-hand-side of the second matrix
!>          equation in (1).
!>          On exit, if IJOB = 0, F has been overwritten by the
!>          solution L.
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
!>          SCALE is DOUBLE PRECISION
!>          On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the solutions
!>          R and L (C and F on entry) will hold the solutions to a
!>          slightly perturbed system but the input matrices A, B, D and
!>          E have not been changed. If SCALE = 0, R and L will hold the
!>          solutions to the homogeneous system with C = F = 0. Normally,
!>          SCALE = 1.
!> \endverbatim
!>
!> \param[in,out] RDSUM
!> \verbatim
!>          RDSUM is DOUBLE PRECISION
!>          On entry, the sum of squares of computed contributions to
!>          the Dif-estimate under computation by DTGSYL, where the
!>          scaling factor RDSCAL (see below) has been factored out.
!>          On exit, the corresponding sum of squares updated with the
!>          contributions from the current sub-system.
!>          If TRANS = 'T' RDSUM is not touched.
!>          NOTE: RDSUM only makes sense when DTGSY2 is called by DTGSYL.
!> \endverbatim
!>
!> \param[in,out] RDSCAL
!> \verbatim
!>          RDSCAL is DOUBLE PRECISION
!>          On entry, scaling factor used to prevent overflow in RDSUM.
!>          On exit, RDSCAL is updated w.r.t. the current contributions
!>          in RDSUM.
!>          If TRANS = 'T', RDSCAL is not touched.
!>          NOTE: RDSCAL only makes sense when DTGSY2 is called by
!>                DTGSYL.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (M+N+2)
!> \endverbatim
!>
!> \param[out] PQ
!> \verbatim
!>          PQ is INTEGER
!>          On exit, the number of subsystems (of size 2-by-2, 4-by-4 and
!>          8-by-8) solved by this routine.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          On exit, if INFO is set to
!>            =0: Successful exit
!>            <0: If INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup doubleSYauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!  =====================================================================
      SUBROUTINE DTGSY2(Trans,Ijob,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F, &
     &                  Ldf,Scale,Rdsum,Rdscal,Iwork,Pq,Info)
      USE F77KINDS                        
      USE S_DAXPY
      USE S_DCOPY
      USE S_DGEMM
      USE S_DGEMV
      USE S_DGER
      USE S_DGESC2
      USE S_DGETC2
      USE S_DLASET
      USE S_DLATDF
      USE S_DSCAL
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DTGSY2290
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  LDZ = 8
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trans
      INTEGER :: Ijob
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(Ldd,*) :: D
      INTEGER :: Ldd
      REAL(R8KIND) , DIMENSION(Lde,*) :: E
      INTEGER :: Lde
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) :: Rdsum
      REAL(R8KIND) :: Rdscal
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) :: Pq
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: alpha , scaloc
      INTEGER :: i , ie , ierr , ii , is , isp1 , j , je , jj , js ,    &
     &           jsp1 , k , mb , nb , p , q , zdim
      INTEGER , DIMENSION(LDZ) :: ipiv , jpiv
      LOGICAL :: notran
      REAL(R8KIND) , DIMENSION(LDZ) :: rhs
      REAL(R8KIND) , DIMENSION(LDZ,LDZ) :: z
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!  Replaced various illegal calls to DCOPY by calls to DLASET.
!  Sven Hammarling, 27/5/02.
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
!     Decode and test input parameters
!
      Info = 0
      ierr = 0
      notran = LSAME(Trans,'N')
      IF ( .NOT.notran .AND. .NOT.LSAME(Trans,'T') ) THEN
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
         CALL XERBLA('DTGSY2',-Info)
         RETURN
      ENDIF
!
!     Determine block structure of A
!
      Pq = 0
      p = 0
      i = 1
      DO WHILE ( i<=M )
         p = p + 1
         Iwork(p) = i
         IF ( i==M ) EXIT
         IF ( A(i+1,i)/=ZERO ) THEN
            i = i + 2
         ELSE
            i = i + 1
         ENDIF
      ENDDO
      Iwork(p+1) = M + 1
!
!     Determine block structure of B
!
      q = p + 1
      j = 1
      DO WHILE ( j<=N )
         q = q + 1
         Iwork(q) = j
         IF ( j==N ) EXIT
         IF ( B(j+1,j)/=ZERO ) THEN
            j = j + 2
         ELSE
            j = j + 1
         ENDIF
      ENDDO
      Iwork(q+1) = N + 1
      Pq = p*(q-p-1)
!
      IF ( notran ) THEN
!
!        Solve (I, J) - subsystem
!           A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
!           D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
!        for I = P, P - 1, ..., 1; J = 1, 2, ..., Q
!
         Scale = ONE
         scaloc = ONE
         DO j = p + 2 , q
            js = Iwork(j)
            jsp1 = js + 1
            je = Iwork(j+1) - 1
            nb = je - js + 1
            DO i = p , 1 , -1
!
               is = Iwork(i)
               isp1 = is + 1
               ie = Iwork(i+1) - 1
               mb = ie - is + 1
               zdim = mb*nb*2
!
               IF ( (mb==1) .AND. (nb==1) ) THEN
!
!                 Build a 2-by-2 system Z * x = RHS
!
                  z(1,1) = A(is,is)
                  z(2,1) = D(is,is)
                  z(1,2) = -B(js,js)
                  z(2,2) = -E(js,js)
!
!                 Set up right hand side(s)
!
                  rhs(1) = C(is,js)
                  rhs(2) = F(is,js)
!
!                 Solve Z * x = RHS
!
                  CALL DGETC2(zdim,z,LDZ,ipiv,jpiv,ierr)
                  IF ( ierr>0 ) Info = ierr
!
                  IF ( Ijob==0 ) THEN
                     CALL DGESC2(zdim,z,LDZ,rhs,ipiv,jpiv,scaloc)
                     IF ( scaloc/=ONE ) THEN
                        DO k = 1 , N
                           CALL DSCAL(M,scaloc,C(1,k),1)
                           CALL DSCAL(M,scaloc,F(1,k),1)
                        ENDDO
                        Scale = Scale*scaloc
                     ENDIF
                  ELSE
                     CALL DLATDF(Ijob,zdim,z,LDZ,rhs,Rdsum,Rdscal,ipiv, &
     &                           jpiv)
                  ENDIF
!
!                 Unpack solution vector(s)
!
                  C(is,js) = rhs(1)
                  F(is,js) = rhs(2)
!
!                 Substitute R(I, J) and L(I, J) into remaining
!                 equation.
!
                  IF ( i>1 ) THEN
                     alpha = -rhs(1)
                     CALL DAXPY(is-1,alpha,A(1,is),1,C(1,js),1)
                     CALL DAXPY(is-1,alpha,D(1,is),1,F(1,js),1)
                  ENDIF
                  IF ( j<q ) THEN
                     CALL DAXPY(N-je,rhs(2),B(js,je+1),Ldb,C(is,je+1),  &
     &                          Ldc)
                     CALL DAXPY(N-je,rhs(2),E(js,je+1),Lde,F(is,je+1),  &
     &                          Ldf)
                  ENDIF
!
               ELSEIF ( (mb==1) .AND. (nb==2) ) THEN
!
!                 Build a 4-by-4 system Z * x = RHS
!
                  z(1,1) = A(is,is)
                  z(2,1) = ZERO
                  z(3,1) = D(is,is)
                  z(4,1) = ZERO
!
                  z(1,2) = ZERO
                  z(2,2) = A(is,is)
                  z(3,2) = ZERO
                  z(4,2) = D(is,is)
!
                  z(1,3) = -B(js,js)
                  z(2,3) = -B(js,jsp1)
                  z(3,3) = -E(js,js)
                  z(4,3) = -E(js,jsp1)
!
                  z(1,4) = -B(jsp1,js)
                  z(2,4) = -B(jsp1,jsp1)
                  z(3,4) = ZERO
                  z(4,4) = -E(jsp1,jsp1)
!
!                 Set up right hand side(s)
!
                  rhs(1) = C(is,js)
                  rhs(2) = C(is,jsp1)
                  rhs(3) = F(is,js)
                  rhs(4) = F(is,jsp1)
!
!                 Solve Z * x = RHS
!
                  CALL DGETC2(zdim,z,LDZ,ipiv,jpiv,ierr)
                  IF ( ierr>0 ) Info = ierr
!
                  IF ( Ijob==0 ) THEN
                     CALL DGESC2(zdim,z,LDZ,rhs,ipiv,jpiv,scaloc)
                     IF ( scaloc/=ONE ) THEN
                        DO k = 1 , N
                           CALL DSCAL(M,scaloc,C(1,k),1)
                           CALL DSCAL(M,scaloc,F(1,k),1)
                        ENDDO
                        Scale = Scale*scaloc
                     ENDIF
                  ELSE
                     CALL DLATDF(Ijob,zdim,z,LDZ,rhs,Rdsum,Rdscal,ipiv, &
     &                           jpiv)
                  ENDIF
!
!                 Unpack solution vector(s)
!
                  C(is,js) = rhs(1)
                  C(is,jsp1) = rhs(2)
                  F(is,js) = rhs(3)
                  F(is,jsp1) = rhs(4)
!
!                 Substitute R(I, J) and L(I, J) into remaining
!                 equation.
!
                  IF ( i>1 ) THEN
                     CALL DGER(is-1,nb,-ONE,A(1,is),1,rhs(1),1,C(1,js), &
     &                         Ldc)
                     CALL DGER(is-1,nb,-ONE,D(1,is),1,rhs(1),1,F(1,js), &
     &                         Ldf)
                  ENDIF
                  IF ( j<q ) THEN
                     CALL DAXPY(N-je,rhs(3),B(js,je+1),Ldb,C(is,je+1),  &
     &                          Ldc)
                     CALL DAXPY(N-je,rhs(3),E(js,je+1),Lde,F(is,je+1),  &
     &                          Ldf)
                     CALL DAXPY(N-je,rhs(4),B(jsp1,je+1),Ldb,C(is,je+1),&
     &                          Ldc)
                     CALL DAXPY(N-je,rhs(4),E(jsp1,je+1),Lde,F(is,je+1),&
     &                          Ldf)
                  ENDIF
!
               ELSEIF ( (mb==2) .AND. (nb==1) ) THEN
!
!                 Build a 4-by-4 system Z * x = RHS
!
                  z(1,1) = A(is,is)
                  z(2,1) = A(isp1,is)
                  z(3,1) = D(is,is)
                  z(4,1) = ZERO
!
                  z(1,2) = A(is,isp1)
                  z(2,2) = A(isp1,isp1)
                  z(3,2) = D(is,isp1)
                  z(4,2) = D(isp1,isp1)
!
                  z(1,3) = -B(js,js)
                  z(2,3) = ZERO
                  z(3,3) = -E(js,js)
                  z(4,3) = ZERO
!
                  z(1,4) = ZERO
                  z(2,4) = -B(js,js)
                  z(3,4) = ZERO
                  z(4,4) = -E(js,js)
!
!                 Set up right hand side(s)
!
                  rhs(1) = C(is,js)
                  rhs(2) = C(isp1,js)
                  rhs(3) = F(is,js)
                  rhs(4) = F(isp1,js)
!
!                 Solve Z * x = RHS
!
                  CALL DGETC2(zdim,z,LDZ,ipiv,jpiv,ierr)
                  IF ( ierr>0 ) Info = ierr
                  IF ( Ijob==0 ) THEN
                     CALL DGESC2(zdim,z,LDZ,rhs,ipiv,jpiv,scaloc)
                     IF ( scaloc/=ONE ) THEN
                        DO k = 1 , N
                           CALL DSCAL(M,scaloc,C(1,k),1)
                           CALL DSCAL(M,scaloc,F(1,k),1)
                        ENDDO
                        Scale = Scale*scaloc
                     ENDIF
                  ELSE
                     CALL DLATDF(Ijob,zdim,z,LDZ,rhs,Rdsum,Rdscal,ipiv, &
     &                           jpiv)
                  ENDIF
!
!                 Unpack solution vector(s)
!
                  C(is,js) = rhs(1)
                  C(isp1,js) = rhs(2)
                  F(is,js) = rhs(3)
                  F(isp1,js) = rhs(4)
!
!                 Substitute R(I, J) and L(I, J) into remaining
!                 equation.
!
                  IF ( i>1 ) THEN
                     CALL DGEMV('N',is-1,mb,-ONE,A(1,is),Lda,rhs(1),1,  &
     &                          ONE,C(1,js),1)
                     CALL DGEMV('N',is-1,mb,-ONE,D(1,is),Ldd,rhs(1),1,  &
     &                          ONE,F(1,js),1)
                  ENDIF
                  IF ( j<q ) THEN
                     CALL DGER(mb,N-je,ONE,rhs(3),1,B(js,je+1),Ldb,     &
     &                         C(is,je+1),Ldc)
                     CALL DGER(mb,N-je,ONE,rhs(3),1,E(js,je+1),Lde,     &
     &                         F(is,je+1),Ldf)
                  ENDIF
!
               ELSEIF ( (mb==2) .AND. (nb==2) ) THEN
!
!                 Build an 8-by-8 system Z * x = RHS
!
                  CALL DLASET('F',LDZ,LDZ,ZERO,ZERO,z,LDZ)
!
                  z(1,1) = A(is,is)
                  z(2,1) = A(isp1,is)
                  z(5,1) = D(is,is)
!
                  z(1,2) = A(is,isp1)
                  z(2,2) = A(isp1,isp1)
                  z(5,2) = D(is,isp1)
                  z(6,2) = D(isp1,isp1)
!
                  z(3,3) = A(is,is)
                  z(4,3) = A(isp1,is)
                  z(7,3) = D(is,is)
!
                  z(3,4) = A(is,isp1)
                  z(4,4) = A(isp1,isp1)
                  z(7,4) = D(is,isp1)
                  z(8,4) = D(isp1,isp1)
!
                  z(1,5) = -B(js,js)
                  z(3,5) = -B(js,jsp1)
                  z(5,5) = -E(js,js)
                  z(7,5) = -E(js,jsp1)
!
                  z(2,6) = -B(js,js)
                  z(4,6) = -B(js,jsp1)
                  z(6,6) = -E(js,js)
                  z(8,6) = -E(js,jsp1)
!
                  z(1,7) = -B(jsp1,js)
                  z(3,7) = -B(jsp1,jsp1)
                  z(7,7) = -E(jsp1,jsp1)
!
                  z(2,8) = -B(jsp1,js)
                  z(4,8) = -B(jsp1,jsp1)
                  z(8,8) = -E(jsp1,jsp1)
!
!                 Set up right hand side(s)
!
                  k = 1
                  ii = mb*nb + 1
                  DO jj = 0 , nb - 1
                     CALL DCOPY(mb,C(is,js+jj),1,rhs(k),1)
                     CALL DCOPY(mb,F(is,js+jj),1,rhs(ii),1)
                     k = k + mb
                     ii = ii + mb
                  ENDDO
!
!                 Solve Z * x = RHS
!
                  CALL DGETC2(zdim,z,LDZ,ipiv,jpiv,ierr)
                  IF ( ierr>0 ) Info = ierr
                  IF ( Ijob==0 ) THEN
                     CALL DGESC2(zdim,z,LDZ,rhs,ipiv,jpiv,scaloc)
                     IF ( scaloc/=ONE ) THEN
                        DO k = 1 , N
                           CALL DSCAL(M,scaloc,C(1,k),1)
                           CALL DSCAL(M,scaloc,F(1,k),1)
                        ENDDO
                        Scale = Scale*scaloc
                     ENDIF
                  ELSE
                     CALL DLATDF(Ijob,zdim,z,LDZ,rhs,Rdsum,Rdscal,ipiv, &
     &                           jpiv)
                  ENDIF
!
!                 Unpack solution vector(s)
!
                  k = 1
                  ii = mb*nb + 1
                  DO jj = 0 , nb - 1
                     CALL DCOPY(mb,rhs(k),1,C(is,js+jj),1)
                     CALL DCOPY(mb,rhs(ii),1,F(is,js+jj),1)
                     k = k + mb
                     ii = ii + mb
                  ENDDO
!
!                 Substitute R(I, J) and L(I, J) into remaining
!                 equation.
!
                  IF ( i>1 ) THEN
                     CALL DGEMM('N','N',is-1,nb,mb,-ONE,A(1,is),Lda,    &
     &                          rhs(1),mb,ONE,C(1,js),Ldc)
                     CALL DGEMM('N','N',is-1,nb,mb,-ONE,D(1,is),Ldd,    &
     &                          rhs(1),mb,ONE,F(1,js),Ldf)
                  ENDIF
                  IF ( j<q ) THEN
                     k = mb*nb + 1
                     CALL DGEMM('N','N',mb,N-je,nb,ONE,rhs(k),mb,       &
     &                          B(js,je+1),Ldb,ONE,C(is,je+1),Ldc)
                     CALL DGEMM('N','N',mb,N-je,nb,ONE,rhs(k),mb,       &
     &                          E(js,je+1),Lde,ONE,F(is,je+1),Ldf)
                  ENDIF
!
               ENDIF
!
            ENDDO
         ENDDO
      ELSE
!
!        Solve (I, J) - subsystem
!             A(I, I)**T * R(I, J) + D(I, I)**T * L(J, J)  =  C(I, J)
!             R(I, I)  * B(J, J) + L(I, J)  * E(J, J)  = -F(I, J)
!        for I = 1, 2, ..., P, J = Q, Q - 1, ..., 1
!
         Scale = ONE
         scaloc = ONE
         DO i = 1 , p
!
            is = Iwork(i)
            isp1 = is + 1
            ie = Iwork(i+1) - 1
            mb = ie - is + 1
            DO j = q , p + 2 , -1
!
               js = Iwork(j)
               jsp1 = js + 1
               je = Iwork(j+1) - 1
               nb = je - js + 1
               zdim = mb*nb*2
               IF ( (mb==1) .AND. (nb==1) ) THEN
!
!                 Build a 2-by-2 system Z**T * x = RHS
!
                  z(1,1) = A(is,is)
                  z(2,1) = -B(js,js)
                  z(1,2) = D(is,is)
                  z(2,2) = -E(js,js)
!
!                 Set up right hand side(s)
!
                  rhs(1) = C(is,js)
                  rhs(2) = F(is,js)
!
!                 Solve Z**T * x = RHS
!
                  CALL DGETC2(zdim,z,LDZ,ipiv,jpiv,ierr)
                  IF ( ierr>0 ) Info = ierr
!
                  CALL DGESC2(zdim,z,LDZ,rhs,ipiv,jpiv,scaloc)
                  IF ( scaloc/=ONE ) THEN
                     DO k = 1 , N
                        CALL DSCAL(M,scaloc,C(1,k),1)
                        CALL DSCAL(M,scaloc,F(1,k),1)
                     ENDDO
                     Scale = Scale*scaloc
                  ENDIF
!
!                 Unpack solution vector(s)
!
                  C(is,js) = rhs(1)
                  F(is,js) = rhs(2)
!
!                 Substitute R(I, J) and L(I, J) into remaining
!                 equation.
!
                  IF ( j>p+2 ) THEN
                     alpha = rhs(1)
                     CALL DAXPY(js-1,alpha,B(1,js),1,F(is,1),Ldf)
                     alpha = rhs(2)
                     CALL DAXPY(js-1,alpha,E(1,js),1,F(is,1),Ldf)
                  ENDIF
                  IF ( i<p ) THEN
                     alpha = -rhs(1)
                     CALL DAXPY(M-ie,alpha,A(is,ie+1),Lda,C(ie+1,js),1)
                     alpha = -rhs(2)
                     CALL DAXPY(M-ie,alpha,D(is,ie+1),Ldd,C(ie+1,js),1)
                  ENDIF
!
               ELSEIF ( (mb==1) .AND. (nb==2) ) THEN
!
!                 Build a 4-by-4 system Z**T * x = RHS
!
                  z(1,1) = A(is,is)
                  z(2,1) = ZERO
                  z(3,1) = -B(js,js)
                  z(4,1) = -B(jsp1,js)
!
                  z(1,2) = ZERO
                  z(2,2) = A(is,is)
                  z(3,2) = -B(js,jsp1)
                  z(4,2) = -B(jsp1,jsp1)
!
                  z(1,3) = D(is,is)
                  z(2,3) = ZERO
                  z(3,3) = -E(js,js)
                  z(4,3) = ZERO
!
                  z(1,4) = ZERO
                  z(2,4) = D(is,is)
                  z(3,4) = -E(js,jsp1)
                  z(4,4) = -E(jsp1,jsp1)
!
!                 Set up right hand side(s)
!
                  rhs(1) = C(is,js)
                  rhs(2) = C(is,jsp1)
                  rhs(3) = F(is,js)
                  rhs(4) = F(is,jsp1)
!
!                 Solve Z**T * x = RHS
!
                  CALL DGETC2(zdim,z,LDZ,ipiv,jpiv,ierr)
                  IF ( ierr>0 ) Info = ierr
                  CALL DGESC2(zdim,z,LDZ,rhs,ipiv,jpiv,scaloc)
                  IF ( scaloc/=ONE ) THEN
                     DO k = 1 , N
                        CALL DSCAL(M,scaloc,C(1,k),1)
                        CALL DSCAL(M,scaloc,F(1,k),1)
                     ENDDO
                     Scale = Scale*scaloc
                  ENDIF
!
!                 Unpack solution vector(s)
!
                  C(is,js) = rhs(1)
                  C(is,jsp1) = rhs(2)
                  F(is,js) = rhs(3)
                  F(is,jsp1) = rhs(4)
!
!                 Substitute R(I, J) and L(I, J) into remaining
!                 equation.
!
                  IF ( j>p+2 ) THEN
                     CALL DAXPY(js-1,rhs(1),B(1,js),1,F(is,1),Ldf)
                     CALL DAXPY(js-1,rhs(2),B(1,jsp1),1,F(is,1),Ldf)
                     CALL DAXPY(js-1,rhs(3),E(1,js),1,F(is,1),Ldf)
                     CALL DAXPY(js-1,rhs(4),E(1,jsp1),1,F(is,1),Ldf)
                  ENDIF
                  IF ( i<p ) THEN
                     CALL DGER(M-ie,nb,-ONE,A(is,ie+1),Lda,rhs(1),1,    &
     &                         C(ie+1,js),Ldc)
                     CALL DGER(M-ie,nb,-ONE,D(is,ie+1),Ldd,rhs(3),1,    &
     &                         C(ie+1,js),Ldc)
                  ENDIF
!
               ELSEIF ( (mb==2) .AND. (nb==1) ) THEN
!
!                 Build a 4-by-4 system Z**T * x = RHS
!
                  z(1,1) = A(is,is)
                  z(2,1) = A(is,isp1)
                  z(3,1) = -B(js,js)
                  z(4,1) = ZERO
!
                  z(1,2) = A(isp1,is)
                  z(2,2) = A(isp1,isp1)
                  z(3,2) = ZERO
                  z(4,2) = -B(js,js)
!
                  z(1,3) = D(is,is)
                  z(2,3) = D(is,isp1)
                  z(3,3) = -E(js,js)
                  z(4,3) = ZERO
!
                  z(1,4) = ZERO
                  z(2,4) = D(isp1,isp1)
                  z(3,4) = ZERO
                  z(4,4) = -E(js,js)
!
!                 Set up right hand side(s)
!
                  rhs(1) = C(is,js)
                  rhs(2) = C(isp1,js)
                  rhs(3) = F(is,js)
                  rhs(4) = F(isp1,js)
!
!                 Solve Z**T * x = RHS
!
                  CALL DGETC2(zdim,z,LDZ,ipiv,jpiv,ierr)
                  IF ( ierr>0 ) Info = ierr
!
                  CALL DGESC2(zdim,z,LDZ,rhs,ipiv,jpiv,scaloc)
                  IF ( scaloc/=ONE ) THEN
                     DO k = 1 , N
                        CALL DSCAL(M,scaloc,C(1,k),1)
                        CALL DSCAL(M,scaloc,F(1,k),1)
                     ENDDO
                     Scale = Scale*scaloc
                  ENDIF
!
!                 Unpack solution vector(s)
!
                  C(is,js) = rhs(1)
                  C(isp1,js) = rhs(2)
                  F(is,js) = rhs(3)
                  F(isp1,js) = rhs(4)
!
!                 Substitute R(I, J) and L(I, J) into remaining
!                 equation.
!
                  IF ( j>p+2 ) THEN
                     CALL DGER(mb,js-1,ONE,rhs(1),1,B(1,js),1,F(is,1),  &
     &                         Ldf)
                     CALL DGER(mb,js-1,ONE,rhs(3),1,E(1,js),1,F(is,1),  &
     &                         Ldf)
                  ENDIF
                  IF ( i<p ) THEN
                     CALL DGEMV('T',mb,M-ie,-ONE,A(is,ie+1),Lda,rhs(1), &
     &                          1,ONE,C(ie+1,js),1)
                     CALL DGEMV('T',mb,M-ie,-ONE,D(is,ie+1),Ldd,rhs(3), &
     &                          1,ONE,C(ie+1,js),1)
                  ENDIF
!
               ELSEIF ( (mb==2) .AND. (nb==2) ) THEN
!
!                 Build an 8-by-8 system Z**T * x = RHS
!
                  CALL DLASET('F',LDZ,LDZ,ZERO,ZERO,z,LDZ)
!
                  z(1,1) = A(is,is)
                  z(2,1) = A(is,isp1)
                  z(5,1) = -B(js,js)
                  z(7,1) = -B(jsp1,js)
!
                  z(1,2) = A(isp1,is)
                  z(2,2) = A(isp1,isp1)
                  z(6,2) = -B(js,js)
                  z(8,2) = -B(jsp1,js)
!
                  z(3,3) = A(is,is)
                  z(4,3) = A(is,isp1)
                  z(5,3) = -B(js,jsp1)
                  z(7,3) = -B(jsp1,jsp1)
!
                  z(3,4) = A(isp1,is)
                  z(4,4) = A(isp1,isp1)
                  z(6,4) = -B(js,jsp1)
                  z(8,4) = -B(jsp1,jsp1)
!
                  z(1,5) = D(is,is)
                  z(2,5) = D(is,isp1)
                  z(5,5) = -E(js,js)
!
                  z(2,6) = D(isp1,isp1)
                  z(6,6) = -E(js,js)
!
                  z(3,7) = D(is,is)
                  z(4,7) = D(is,isp1)
                  z(5,7) = -E(js,jsp1)
                  z(7,7) = -E(jsp1,jsp1)
!
                  z(4,8) = D(isp1,isp1)
                  z(6,8) = -E(js,jsp1)
                  z(8,8) = -E(jsp1,jsp1)
!
!                 Set up right hand side(s)
!
                  k = 1
                  ii = mb*nb + 1
                  DO jj = 0 , nb - 1
                     CALL DCOPY(mb,C(is,js+jj),1,rhs(k),1)
                     CALL DCOPY(mb,F(is,js+jj),1,rhs(ii),1)
                     k = k + mb
                     ii = ii + mb
                  ENDDO
!
!
!                 Solve Z**T * x = RHS
!
                  CALL DGETC2(zdim,z,LDZ,ipiv,jpiv,ierr)
                  IF ( ierr>0 ) Info = ierr
!
                  CALL DGESC2(zdim,z,LDZ,rhs,ipiv,jpiv,scaloc)
                  IF ( scaloc/=ONE ) THEN
                     DO k = 1 , N
                        CALL DSCAL(M,scaloc,C(1,k),1)
                        CALL DSCAL(M,scaloc,F(1,k),1)
                     ENDDO
                     Scale = Scale*scaloc
                  ENDIF
!
!                 Unpack solution vector(s)
!
                  k = 1
                  ii = mb*nb + 1
                  DO jj = 0 , nb - 1
                     CALL DCOPY(mb,rhs(k),1,C(is,js+jj),1)
                     CALL DCOPY(mb,rhs(ii),1,F(is,js+jj),1)
                     k = k + mb
                     ii = ii + mb
                  ENDDO
!
!                 Substitute R(I, J) and L(I, J) into remaining
!                 equation.
!
                  IF ( j>p+2 ) THEN
                     CALL DGEMM('N','T',mb,js-1,nb,ONE,C(is,js),Ldc,    &
     &                          B(1,js),Ldb,ONE,F(is,1),Ldf)
                     CALL DGEMM('N','T',mb,js-1,nb,ONE,F(is,js),Ldf,    &
     &                          E(1,js),Lde,ONE,F(is,1),Ldf)
                  ENDIF
                  IF ( i<p ) THEN
                     CALL DGEMM('T','N',M-ie,nb,mb,-ONE,A(is,ie+1),Lda, &
     &                          C(is,js),Ldc,ONE,C(ie+1,js),Ldc)
                     CALL DGEMM('T','N',M-ie,nb,mb,-ONE,D(is,ie+1),Ldd, &
     &                          F(is,js),Ldf,ONE,C(ie+1,js),Ldc)
                  ENDIF
!
               ENDIF
!
            ENDDO
         ENDDO
!
      ENDIF
!
!     End of DTGSY2
!
      END SUBROUTINE DTGSY2
