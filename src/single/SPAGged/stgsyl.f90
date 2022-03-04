!*==stgsyl.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b STGSYL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STGSYL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsyl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsyl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsyl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,
!                          LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK,
!                          IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF,
!      $                   LWORK, M, N
!       REAL               DIF, SCALE
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               A( LDA, * ), B( LDB, * ), C( LDC, * ),
!      $                   D( LDD, * ), E( LDE, * ), F( LDF, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STGSYL solves the generalized Sylvester equation:
!>
!>             A * R - L * B = scale * C                 (1)
!>             D * R - L * E = scale * F
!>
!> where R and L are unknown m-by-n matrices, (A, D), (B, E) and
!> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
!> respectively, with real entries. (A, D) and (B, E) must be in
!> generalized (real) Schur canonical form, i.e. A, B are upper quasi
!> triangular and D, E are upper triangular.
!>
!> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
!> scaling factor chosen to avoid overflow.
!>
!> In matrix notation (1) is equivalent to solve  Zx = scale b, where
!> Z is defined as
!>
!>            Z = [ kron(In, A)  -kron(B**T, Im) ]         (2)
!>                [ kron(In, D)  -kron(E**T, Im) ].
!>
!> Here Ik is the identity matrix of size k and X**T is the transpose of
!> X. kron(X, Y) is the Kronecker product between the matrices X and Y.
!>
!> If TRANS = 'T', STGSYL solves the transposed system Z**T*y = scale*b,
!> which is equivalent to solve for R and L in
!>
!>             A**T * R + D**T * L = scale * C           (3)
!>             R * B**T + L * E**T = scale * -F
!>
!> This case (TRANS = 'T') is used to compute an one-norm-based estimate
!> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
!> and (B,E), using SLACON.
!>
!> If IJOB >= 1, STGSYL computes a Frobenius norm-based estimate
!> of Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
!> reciprocal of the smallest singular value of Z. See [1-2] for more
!> information.
!>
!> This is a level 3 BLAS algorithm.
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
!>          = 1: The functionality of 0 and 3.
!>          = 2: The functionality of 0 and 4.
!>          = 3: Only an estimate of Dif[(A,D), (B,E)] is computed.
!>               (look ahead strategy IJOB  = 1 is used).
!>          = 4: Only an estimate of Dif[(A,D), (B,E)] is computed.
!>               ( SGECON on sub-systems is used ).
!>          Not referenced if TRANS = 'T'.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and D, and the row dimension of
!>          the matrices C, F, R and L.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices B and E, and the column dimension
!>          of the matrices C, F, R and L.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA, M)
!>          The upper quasi triangular matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1, M).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB, N)
!>          The upper quasi triangular matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1, N).
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is REAL array, dimension (LDC, N)
!>          On entry, C contains the right-hand-side of the first matrix
!>          equation in (1) or (3).
!>          On exit, if IJOB = 0, 1 or 2, C has been overwritten by
!>          the solution R. If IJOB = 3 or 4 and TRANS = 'N', C holds R,
!>          the solution achieved during the computation of the
!>          Dif-estimate.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1, M).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (LDD, M)
!>          The upper triangular matrix D.
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of the array D. LDD >= max(1, M).
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (LDE, N)
!>          The upper triangular matrix E.
!> \endverbatim
!>
!> \param[in] LDE
!> \verbatim
!>          LDE is INTEGER
!>          The leading dimension of the array E. LDE >= max(1, N).
!> \endverbatim
!>
!> \param[in,out] F
!> \verbatim
!>          F is REAL array, dimension (LDF, N)
!>          On entry, F contains the right-hand-side of the second matrix
!>          equation in (1) or (3).
!>          On exit, if IJOB = 0, 1 or 2, F has been overwritten by
!>          the solution L. If IJOB = 3 or 4 and TRANS = 'N', F holds L,
!>          the solution achieved during the computation of the
!>          Dif-estimate.
!> \endverbatim
!>
!> \param[in] LDF
!> \verbatim
!>          LDF is INTEGER
!>          The leading dimension of the array F. LDF >= max(1, M).
!> \endverbatim
!>
!> \param[out] DIF
!> \verbatim
!>          DIF is REAL
!>          On exit DIF is the reciprocal of a lower bound of the
!>          reciprocal of the Dif-function, i.e. DIF is an upper bound of
!>          Dif[(A,D), (B,E)] = sigma_min(Z), where Z as in (2).
!>          IF IJOB = 0 or TRANS = 'T', DIF is not touched.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is REAL
!>          On exit SCALE is the scaling factor in (1) or (3).
!>          If 0 < SCALE < 1, C and F hold the solutions R and L, resp.,
!>          to a slightly perturbed system but the input matrices A, B, D
!>          and E have not been changed. If SCALE = 0, C and F hold the
!>          solutions R and L, respectively, to the homogeneous system
!>          with C = F = 0. Normally, SCALE = 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK > = 1.
!>          If IJOB = 1 or 2 and TRANS = 'N', LWORK >= max(1,2*M*N).
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (M+N+6)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>            =0: successful exit
!>            <0: If INFO = -i, the i-th argument had an illegal value.
!>            >0: (A, D) and (B, E) have common or close eigenvalues.
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
!> \ingroup realSYcomputational
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
!>  [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
!>      for Solving the Generalized Sylvester Equation and Estimating the
!>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,
!>      Department of Computing Science, Umea University, S-901 87 Umea,
!>      Sweden, December 1993, Revised April 1994, Also as LAPACK Working
!>      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,
!>      No 1, 1996.
!>
!>  [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester
!>      Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal.
!>      Appl., 15(4):1045-1060, 1994
!>
!>  [3] B. Kagstrom and L. Westin, Generalized Schur Methods with
!>      Condition Estimators for Solving the Generalized Sylvester
!>      Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7,
!>      July 1989, pp 745-751.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE STGSYL(Trans,Ijob,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F, &
     &                  Ldf,Scale,Dif,Work,Lwork,Iwork,Info)
      USE S_ILAENV
      USE S_LSAME
      USE S_SGEMM
      USE S_SLACPY
      USE S_SLASET
      USE S_SSCAL
      USE S_STGSY2
      USE S_XERBLA
      IMPLICIT NONE
!*--STGSYL310
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: Ijob
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(Ldd,*) :: D
      INTEGER :: Ldd
      REAL , DIMENSION(Lde,*) :: E
      INTEGER :: Lde
      REAL , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(OUT) :: Dif
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: dscale , dsum , scale2 , scaloc
      INTEGER :: i , ie , ifunc , iround , is , isolve , j , je , js ,  &
     &           k , linfo , lwmin , mb , nb , p , ppqq , pq , q
      LOGICAL :: lquery , notran
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!  Replaced various illegal calls to SCOPY by calls to SLASET.
!  Sven Hammarling, 1/5/02.
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
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
      notran = LSAME(Trans,'N')
      lquery = (Lwork==-1)
!
      IF ( .NOT.notran .AND. .NOT.LSAME(Trans,'T') ) THEN
         Info = -1
      ELSEIF ( notran ) THEN
         IF ( (Ijob<0) .OR. (Ijob>4) ) Info = -2
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
!
      IF ( Info==0 ) THEN
         IF ( .NOT.(notran) ) THEN
            lwmin = 1
         ELSEIF ( Ijob==1 .OR. Ijob==2 ) THEN
            lwmin = MAX(1,2*M*N)
         ELSE
            lwmin = 1
         ENDIF
         Work(1) = lwmin
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) Info = -20
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('STGSYL',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) THEN
         Scale = 1
         IF ( notran ) THEN
            IF ( Ijob/=0 ) Dif = 0
         ENDIF
         RETURN
      ENDIF
!
!     Determine optimal block sizes MB and NB
!
      mb = ILAENV(2,'STGSYL',Trans,M,N,-1,-1)
      nb = ILAENV(5,'STGSYL',Trans,M,N,-1,-1)
!
      isolve = 1
      ifunc = 0
      IF ( notran ) THEN
         IF ( Ijob>=3 ) THEN
            ifunc = Ijob - 2
            CALL SLASET('F',M,N,ZERO,ZERO,C,Ldc)
            CALL SLASET('F',M,N,ZERO,ZERO,F,Ldf)
         ELSEIF ( Ijob>=1 .AND. notran ) THEN
            isolve = 2
         ENDIF
      ENDIF
!
      IF ( (mb<=1 .AND. nb<=1) .OR. (mb>=M .AND. nb>=N) ) THEN
!
         DO iround = 1 , isolve
!
!           Use unblocked Level 2 solver
!
            dscale = ZERO
            dsum = ONE
            pq = 0
            CALL STGSY2(Trans,ifunc,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F,&
     &                  Ldf,Scale,dsum,dscale,Iwork,pq,Info)
            IF ( dscale/=ZERO ) THEN
               IF ( Ijob==1 .OR. Ijob==3 ) THEN
                  Dif = SQRT(REAL(2*M*N))/(dscale*SQRT(dsum))
               ELSE
                  Dif = SQRT(REAL(pq))/(dscale*SQRT(dsum))
               ENDIF
            ENDIF
!
            IF ( isolve==2 .AND. iround==1 ) THEN
               IF ( notran ) ifunc = Ijob
               scale2 = Scale
               CALL SLACPY('F',M,N,C,Ldc,Work,M)
               CALL SLACPY('F',M,N,F,Ldf,Work(M*N+1),M)
               CALL SLASET('F',M,N,ZERO,ZERO,C,Ldc)
               CALL SLASET('F',M,N,ZERO,ZERO,F,Ldf)
            ELSEIF ( isolve==2 .AND. iround==2 ) THEN
               CALL SLACPY('F',M,N,Work,M,C,Ldc)
               CALL SLACPY('F',M,N,Work(M*N+1),M,F,Ldf)
               Scale = scale2
            ENDIF
         ENDDO
!
         RETURN
      ENDIF
!
!     Determine block structure of A
!
      p = 0
      i = 1
      DO WHILE ( i<=M )
         p = p + 1
         Iwork(p) = i
         i = i + mb
         IF ( i>=M ) EXIT
         IF ( A(i,i-1)/=ZERO ) i = i + 1
      ENDDO
!
      Iwork(p+1) = M + 1
      IF ( Iwork(p)==Iwork(p+1) ) p = p - 1
!
!     Determine block structure of B
!
      q = p + 1
      j = 1
      DO WHILE ( j<=N )
         q = q + 1
         Iwork(q) = j
         j = j + nb
         IF ( j>=N ) EXIT
         IF ( B(j,j-1)/=ZERO ) j = j + 1
      ENDDO
!
      Iwork(q+1) = N + 1
      IF ( Iwork(q)==Iwork(q+1) ) q = q - 1
!
      IF ( notran ) THEN
!
         DO iround = 1 , isolve
!
!           Solve (I, J)-subsystem
!               A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
!               D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
!           for I = P, P - 1,..., 1; J = 1, 2,..., Q
!
            dscale = ZERO
            dsum = ONE
            pq = 0
            Scale = ONE
            DO j = p + 2 , q
               js = Iwork(j)
               je = Iwork(j+1) - 1
               nb = je - js + 1
               DO i = p , 1 , -1
                  is = Iwork(i)
                  ie = Iwork(i+1) - 1
                  mb = ie - is + 1
                  ppqq = 0
                  CALL STGSY2(Trans,ifunc,mb,nb,A(is,is),Lda,B(js,js),  &
     &                        Ldb,C(is,js),Ldc,D(is,is),Ldd,E(js,js),   &
     &                        Lde,F(is,js),Ldf,scaloc,dsum,dscale,      &
     &                        Iwork(q+2),ppqq,linfo)
                  IF ( linfo>0 ) Info = linfo
!
                  pq = pq + ppqq
                  IF ( scaloc/=ONE ) THEN
                     DO k = 1 , js - 1
                        CALL SSCAL(M,scaloc,C(1,k),1)
                        CALL SSCAL(M,scaloc,F(1,k),1)
                     ENDDO
                     DO k = js , je
                        CALL SSCAL(is-1,scaloc,C(1,k),1)
                        CALL SSCAL(is-1,scaloc,F(1,k),1)
                     ENDDO
                     DO k = js , je
                        CALL SSCAL(M-ie,scaloc,C(ie+1,k),1)
                        CALL SSCAL(M-ie,scaloc,F(ie+1,k),1)
                     ENDDO
                     DO k = je + 1 , N
                        CALL SSCAL(M,scaloc,C(1,k),1)
                        CALL SSCAL(M,scaloc,F(1,k),1)
                     ENDDO
                     Scale = Scale*scaloc
                  ENDIF
!
!                 Substitute R(I, J) and L(I, J) into remaining
!                 equation.
!
                  IF ( i>1 ) THEN
                     CALL SGEMM('N','N',is-1,nb,mb,-ONE,A(1,is),Lda,    &
     &                          C(is,js),Ldc,ONE,C(1,js),Ldc)
                     CALL SGEMM('N','N',is-1,nb,mb,-ONE,D(1,is),Ldd,    &
     &                          C(is,js),Ldc,ONE,F(1,js),Ldf)
                  ENDIF
                  IF ( j<q ) THEN
                     CALL SGEMM('N','N',mb,N-je,nb,ONE,F(is,js),Ldf,    &
     &                          B(js,je+1),Ldb,ONE,C(is,je+1),Ldc)
                     CALL SGEMM('N','N',mb,N-je,nb,ONE,F(is,js),Ldf,    &
     &                          E(js,je+1),Lde,ONE,F(is,je+1),Ldf)
                  ENDIF
               ENDDO
            ENDDO
            IF ( dscale/=ZERO ) THEN
               IF ( Ijob==1 .OR. Ijob==3 ) THEN
                  Dif = SQRT(REAL(2*M*N))/(dscale*SQRT(dsum))
               ELSE
                  Dif = SQRT(REAL(pq))/(dscale*SQRT(dsum))
               ENDIF
            ENDIF
            IF ( isolve==2 .AND. iround==1 ) THEN
               IF ( notran ) ifunc = Ijob
               scale2 = Scale
               CALL SLACPY('F',M,N,C,Ldc,Work,M)
               CALL SLACPY('F',M,N,F,Ldf,Work(M*N+1),M)
               CALL SLASET('F',M,N,ZERO,ZERO,C,Ldc)
               CALL SLASET('F',M,N,ZERO,ZERO,F,Ldf)
            ELSEIF ( isolve==2 .AND. iround==2 ) THEN
               CALL SLACPY('F',M,N,Work,M,C,Ldc)
               CALL SLACPY('F',M,N,Work(M*N+1),M,F,Ldf)
               Scale = scale2
            ENDIF
         ENDDO
!
      ELSE
!
!        Solve transposed (I, J)-subsystem
!             A(I, I)**T * R(I, J)  + D(I, I)**T * L(I, J)  =  C(I, J)
!             R(I, J)  * B(J, J)**T + L(I, J)  * E(J, J)**T = -F(I, J)
!        for I = 1,2,..., P; J = Q, Q-1,..., 1
!
         Scale = ONE
         DO i = 1 , p
            is = Iwork(i)
            ie = Iwork(i+1) - 1
            mb = ie - is + 1
            DO j = q , p + 2 , -1
               js = Iwork(j)
               je = Iwork(j+1) - 1
               nb = je - js + 1
               CALL STGSY2(Trans,ifunc,mb,nb,A(is,is),Lda,B(js,js),Ldb, &
     &                     C(is,js),Ldc,D(is,is),Ldd,E(js,js),Lde,      &
     &                     F(is,js),Ldf,scaloc,dsum,dscale,Iwork(q+2),  &
     &                     ppqq,linfo)
               IF ( linfo>0 ) Info = linfo
               IF ( scaloc/=ONE ) THEN
                  DO k = 1 , js - 1
                     CALL SSCAL(M,scaloc,C(1,k),1)
                     CALL SSCAL(M,scaloc,F(1,k),1)
                  ENDDO
                  DO k = js , je
                     CALL SSCAL(is-1,scaloc,C(1,k),1)
                     CALL SSCAL(is-1,scaloc,F(1,k),1)
                  ENDDO
                  DO k = js , je
                     CALL SSCAL(M-ie,scaloc,C(ie+1,k),1)
                     CALL SSCAL(M-ie,scaloc,F(ie+1,k),1)
                  ENDDO
                  DO k = je + 1 , N
                     CALL SSCAL(M,scaloc,C(1,k),1)
                     CALL SSCAL(M,scaloc,F(1,k),1)
                  ENDDO
                  Scale = Scale*scaloc
               ENDIF
!
!              Substitute R(I, J) and L(I, J) into remaining equation.
!
               IF ( j>p+2 ) THEN
                  CALL SGEMM('N','T',mb,js-1,nb,ONE,C(is,js),Ldc,B(1,js)&
     &                       ,Ldb,ONE,F(is,1),Ldf)
                  CALL SGEMM('N','T',mb,js-1,nb,ONE,F(is,js),Ldf,E(1,js)&
     &                       ,Lde,ONE,F(is,1),Ldf)
               ENDIF
               IF ( i<p ) THEN
                  CALL SGEMM('T','N',M-ie,nb,mb,-ONE,A(is,ie+1),Lda,    &
     &                       C(is,js),Ldc,ONE,C(ie+1,js),Ldc)
                  CALL SGEMM('T','N',M-ie,nb,mb,-ONE,D(is,ie+1),Ldd,    &
     &                       F(is,js),Ldf,ONE,C(ie+1,js),Ldc)
               ENDIF
            ENDDO
         ENDDO
!
      ENDIF
!
      Work(1) = lwmin
!
!
!     End of STGSYL
!
      END SUBROUTINE STGSYL
