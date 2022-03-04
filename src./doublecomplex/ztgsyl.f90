!*==ztgsyl.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZTGSYL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTGSYL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsyl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsyl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsyl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,
!                          LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK,
!                          IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF,
!      $                   LWORK, M, N
!       DOUBLE PRECISION   DIF, SCALE
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ),
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
!> ZTGSYL solves the generalized Sylvester equation:
!>
!>             A * R - L * B = scale * C            (1)
!>             D * R - L * E = scale * F
!>
!> where R and L are unknown m-by-n matrices, (A, D), (B, E) and
!> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
!> respectively, with complex entries. A, B, D and E are upper
!> triangular (i.e., (A,D) and (B,E) in generalized Schur form).
!>
!> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1
!> is an output scaling factor chosen to avoid overflow.
!>
!> In matrix notation (1) is equivalent to solve Zx = scale*b, where Z
!> is defined as
!>
!>        Z = [ kron(In, A)  -kron(B**H, Im) ]        (2)
!>            [ kron(In, D)  -kron(E**H, Im) ],
!>
!> Here Ix is the identity matrix of size x and X**H is the conjugate
!> transpose of X. Kron(X, Y) is the Kronecker product between the
!> matrices X and Y.
!>
!> If TRANS = 'C', y in the conjugate transposed system Z**H *y = scale*b
!> is solved for, which is equivalent to solve for R and L in
!>
!>             A**H * R + D**H * L = scale * C           (3)
!>             R * B**H + L * E**H = scale * -F
!>
!> This case (TRANS = 'C') is used to compute an one-norm-based estimate
!> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
!> and (B,E), using ZLACON.
!>
!> If IJOB >= 1, ZTGSYL computes a Frobenius norm-based estimate of
!> Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
!> reciprocal of the smallest singular value of Z.
!>
!> This is a level-3 BLAS algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': solve the generalized sylvester equation (1).
!>          = 'C': solve the "conjugate transposed" system (3).
!> \endverbatim
!>
!> \param[in] IJOB
!> \verbatim
!>          IJOB is INTEGER
!>          Specifies what kind of functionality to be performed.
!>          =0: solve (1) only.
!>          =1: The functionality of 0 and 3.
!>          =2: The functionality of 0 and 4.
!>          =3: Only an estimate of Dif[(A,D), (B,E)] is computed.
!>              (look ahead strategy is used).
!>          =4: Only an estimate of Dif[(A,D), (B,E)] is computed.
!>              (ZGECON on sub-systems is used).
!>          Not referenced if TRANS = 'C'.
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
!>          A is COMPLEX*16 array, dimension (LDA, M)
!>          The upper triangular matrix A.
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
!>          B is COMPLEX*16 array, dimension (LDB, N)
!>          The upper triangular matrix B.
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
!>          C is COMPLEX*16 array, dimension (LDC, N)
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
!>          D is COMPLEX*16 array, dimension (LDD, M)
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
!>          E is COMPLEX*16 array, dimension (LDE, N)
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
!>          F is COMPLEX*16 array, dimension (LDF, N)
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
!>          DIF is DOUBLE PRECISION
!>          On exit DIF is the reciprocal of a lower bound of the
!>          reciprocal of the Dif-function, i.e. DIF is an upper bound of
!>          Dif[(A,D), (B,E)] = sigma-min(Z), where Z as in (2).
!>          IF IJOB = 0 or TRANS = 'C', DIF is not referenced.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          On exit SCALE is the scaling factor in (1) or (3).
!>          If 0 < SCALE < 1, C and F hold the solutions R and L, resp.,
!>          to a slightly perturbed system but the input matrices A, B,
!>          D and E have not been changed. If SCALE = 0, R and L will
!>          hold the solutions to the homogeneous system with C = F = 0.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
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
!>          IWORK is INTEGER array, dimension (M+N+2)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>            =0: successful exit
!>            <0: If INFO = -i, the i-th argument had an illegal value.
!>            >0: (A, D) and (B, E) have common or very close
!>                eigenvalues.
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
!> \ingroup complex16SYcomputational
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
!>  [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
!>      for Solving the Generalized Sylvester Equation and Estimating the
!>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,
!>      Department of Computing Science, Umea University, S-901 87 Umea,
!>      Sweden, December 1993, Revised April 1994, Also as LAPACK Working
!>      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,
!>      No 1, 1996.
!> \n
!>  [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester
!>      Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal.
!>      Appl., 15(4):1045-1060, 1994.
!> \n
!>  [3] B. Kagstrom and L. Westin, Generalized Schur Methods with
!>      Condition Estimators for Solving the Generalized Sylvester
!>      Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7,
!>      July 1989, pp 745-751.
!>
!  =====================================================================
      SUBROUTINE ZTGSYL(Trans,Ijob,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F, &
     &                  Ldf,Scale,Dif,Work,Lwork,Iwork,Info)
      USE F77KINDS                        
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      USE S_ZGEMM
      USE S_ZLACPY
      USE S_ZLASET
      USE S_ZSCAL
      USE S_ZTGSY2
      IMPLICIT NONE
!*--ZTGSYL307
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: Ijob
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(Ldd,*) :: D
      INTEGER :: Ldd
      COMPLEX(CX16KIND) , DIMENSION(Lde,*) :: E
      INTEGER :: Lde
      COMPLEX(CX16KIND) , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(OUT) :: Dif
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: dscale , dsum , scale2 , scaloc
      INTEGER :: i , ie , ifunc , iround , is , isolve , j , je , js ,  &
     &           k , linfo , lwmin , mb , nb , p , pq , q
      LOGICAL :: lquery , notran
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!  Replaced various illegal calls to CCOPY by calls to CLASET.
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
      IF ( .NOT.notran .AND. .NOT.LSAME(Trans,'C') ) THEN
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
         CALL XERBLA('ZTGSYL',-Info)
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
!     Determine  optimal block sizes MB and NB
!
      mb = ILAENV(2,'ZTGSYL',Trans,M,N,-1,-1)
      nb = ILAENV(5,'ZTGSYL',Trans,M,N,-1,-1)
!
      isolve = 1
      ifunc = 0
      IF ( notran ) THEN
         IF ( Ijob>=3 ) THEN
            ifunc = Ijob - 2
            CALL ZLASET('F',M,N,CZERO,CZERO,C,Ldc)
            CALL ZLASET('F',M,N,CZERO,CZERO,F,Ldf)
         ELSEIF ( Ijob>=1 .AND. notran ) THEN
            isolve = 2
         ENDIF
      ENDIF
!
      IF ( (mb<=1 .AND. nb<=1) .OR. (mb>=M .AND. nb>=N) ) THEN
!
!        Use unblocked Level 2 solver
!
         DO iround = 1 , isolve
!
            Scale = ONE
            dscale = ZERO
            dsum = ONE
            pq = M*N
            CALL ZTGSY2(Trans,ifunc,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F,&
     &                  Ldf,Scale,dsum,dscale,Info)
            IF ( dscale/=ZERO ) THEN
               IF ( Ijob==1 .OR. Ijob==3 ) THEN
                  Dif = SQRT(DBLE(2*M*N))/(dscale*SQRT(dsum))
               ELSE
                  Dif = SQRT(DBLE(pq))/(dscale*SQRT(dsum))
               ENDIF
            ENDIF
            IF ( isolve==2 .AND. iround==1 ) THEN
               IF ( notran ) ifunc = Ijob
               scale2 = Scale
               CALL ZLACPY('F',M,N,C,Ldc,Work,M)
               CALL ZLACPY('F',M,N,F,Ldf,Work(M*N+1),M)
               CALL ZLASET('F',M,N,CZERO,CZERO,C,Ldc)
               CALL ZLASET('F',M,N,CZERO,CZERO,F,Ldf)
            ELSEIF ( isolve==2 .AND. iround==2 ) THEN
               CALL ZLACPY('F',M,N,Work,M,C,Ldc)
               CALL ZLACPY('F',M,N,Work(M*N+1),M,F,Ldf)
               Scale = scale2
            ENDIF
         ENDDO
!
         RETURN
!
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
      ENDDO
      Iwork(p+1) = M + 1
      IF ( Iwork(p)==Iwork(p+1) ) p = p - 1
!
!     Determine block structure of B
!
      q = p + 1
      j = 1
      DO WHILE ( j<=N )
!
         q = q + 1
         Iwork(q) = j
         j = j + nb
         IF ( j>=N ) EXIT
      ENDDO
!
      Iwork(q+1) = N + 1
      IF ( Iwork(q)==Iwork(q+1) ) q = q - 1
!
      IF ( notran ) THEN
         DO iround = 1 , isolve
!
!           Solve (I, J) - subsystem
!               A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
!               D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
!           for I = P, P - 1, ..., 1; J = 1, 2, ..., Q
!
            pq = 0
            Scale = ONE
            dscale = ZERO
            dsum = ONE
            DO j = p + 2 , q
               js = Iwork(j)
               je = Iwork(j+1) - 1
               nb = je - js + 1
               DO i = p , 1 , -1
                  is = Iwork(i)
                  ie = Iwork(i+1) - 1
                  mb = ie - is + 1
                  CALL ZTGSY2(Trans,ifunc,mb,nb,A(is,is),Lda,B(js,js),  &
     &                        Ldb,C(is,js),Ldc,D(is,is),Ldd,E(js,js),   &
     &                        Lde,F(is,js),Ldf,scaloc,dsum,dscale,linfo)
                  IF ( linfo>0 ) Info = linfo
                  pq = pq + mb*nb
                  IF ( scaloc/=ONE ) THEN
                     DO k = 1 , js - 1
                        CALL ZSCAL(M,DCMPLX(scaloc,ZERO),C(1,k),1)
                        CALL ZSCAL(M,DCMPLX(scaloc,ZERO),F(1,k),1)
                     ENDDO
                     DO k = js , je
                        CALL ZSCAL(is-1,DCMPLX(scaloc,ZERO),C(1,k),1)
                        CALL ZSCAL(is-1,DCMPLX(scaloc,ZERO),F(1,k),1)
                     ENDDO
                     DO k = js , je
                        CALL ZSCAL(M-ie,DCMPLX(scaloc,ZERO),C(ie+1,k),1)
                        CALL ZSCAL(M-ie,DCMPLX(scaloc,ZERO),F(ie+1,k),1)
                     ENDDO
                     DO k = je + 1 , N
                        CALL ZSCAL(M,DCMPLX(scaloc,ZERO),C(1,k),1)
                        CALL ZSCAL(M,DCMPLX(scaloc,ZERO),F(1,k),1)
                     ENDDO
                     Scale = Scale*scaloc
                  ENDIF
!
!                 Substitute R(I,J) and L(I,J) into remaining equation.
!
                  IF ( i>1 ) THEN
                     CALL ZGEMM('N','N',is-1,nb,mb,DCMPLX(-ONE,ZERO),   &
     &                          A(1,is),Lda,C(is,js),Ldc,               &
     &                          DCMPLX(ONE,ZERO),C(1,js),Ldc)
                     CALL ZGEMM('N','N',is-1,nb,mb,DCMPLX(-ONE,ZERO),   &
     &                          D(1,is),Ldd,C(is,js),Ldc,               &
     &                          DCMPLX(ONE,ZERO),F(1,js),Ldf)
                  ENDIF
                  IF ( j<q ) THEN
                     CALL ZGEMM('N','N',mb,N-je,nb,DCMPLX(ONE,ZERO),    &
     &                          F(is,js),Ldf,B(js,je+1),Ldb,            &
     &                          DCMPLX(ONE,ZERO),C(is,je+1),Ldc)
                     CALL ZGEMM('N','N',mb,N-je,nb,DCMPLX(ONE,ZERO),    &
     &                          F(is,js),Ldf,E(js,je+1),Lde,            &
     &                          DCMPLX(ONE,ZERO),F(is,je+1),Ldf)
                  ENDIF
               ENDDO
            ENDDO
            IF ( dscale/=ZERO ) THEN
               IF ( Ijob==1 .OR. Ijob==3 ) THEN
                  Dif = SQRT(DBLE(2*M*N))/(dscale*SQRT(dsum))
               ELSE
                  Dif = SQRT(DBLE(pq))/(dscale*SQRT(dsum))
               ENDIF
            ENDIF
            IF ( isolve==2 .AND. iround==1 ) THEN
               IF ( notran ) ifunc = Ijob
               scale2 = Scale
               CALL ZLACPY('F',M,N,C,Ldc,Work,M)
               CALL ZLACPY('F',M,N,F,Ldf,Work(M*N+1),M)
               CALL ZLASET('F',M,N,CZERO,CZERO,C,Ldc)
               CALL ZLASET('F',M,N,CZERO,CZERO,F,Ldf)
            ELSEIF ( isolve==2 .AND. iround==2 ) THEN
               CALL ZLACPY('F',M,N,Work,M,C,Ldc)
               CALL ZLACPY('F',M,N,Work(M*N+1),M,F,Ldf)
               Scale = scale2
            ENDIF
         ENDDO
      ELSE
!
!        Solve transposed (I, J)-subsystem
!            A(I, I)**H * R(I, J) + D(I, I)**H * L(I, J) = C(I, J)
!            R(I, J) * B(J, J)  + L(I, J) * E(J, J) = -F(I, J)
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
               CALL ZTGSY2(Trans,ifunc,mb,nb,A(is,is),Lda,B(js,js),Ldb, &
     &                     C(is,js),Ldc,D(is,is),Ldd,E(js,js),Lde,      &
     &                     F(is,js),Ldf,scaloc,dsum,dscale,linfo)
               IF ( linfo>0 ) Info = linfo
               IF ( scaloc/=ONE ) THEN
                  DO k = 1 , js - 1
                     CALL ZSCAL(M,DCMPLX(scaloc,ZERO),C(1,k),1)
                     CALL ZSCAL(M,DCMPLX(scaloc,ZERO),F(1,k),1)
                  ENDDO
                  DO k = js , je
                     CALL ZSCAL(is-1,DCMPLX(scaloc,ZERO),C(1,k),1)
                     CALL ZSCAL(is-1,DCMPLX(scaloc,ZERO),F(1,k),1)
                  ENDDO
                  DO k = js , je
                     CALL ZSCAL(M-ie,DCMPLX(scaloc,ZERO),C(ie+1,k),1)
                     CALL ZSCAL(M-ie,DCMPLX(scaloc,ZERO),F(ie+1,k),1)
                  ENDDO
                  DO k = je + 1 , N
                     CALL ZSCAL(M,DCMPLX(scaloc,ZERO),C(1,k),1)
                     CALL ZSCAL(M,DCMPLX(scaloc,ZERO),F(1,k),1)
                  ENDDO
                  Scale = Scale*scaloc
               ENDIF
!
!              Substitute R(I,J) and L(I,J) into remaining equation.
!
               IF ( j>p+2 ) THEN
                  CALL ZGEMM('N','C',mb,js-1,nb,DCMPLX(ONE,ZERO),       &
     &                       C(is,js),Ldc,B(1,js),Ldb,DCMPLX(ONE,ZERO), &
     &                       F(is,1),Ldf)
                  CALL ZGEMM('N','C',mb,js-1,nb,DCMPLX(ONE,ZERO),       &
     &                       F(is,js),Ldf,E(1,js),Lde,DCMPLX(ONE,ZERO), &
     &                       F(is,1),Ldf)
               ENDIF
               IF ( i<p ) THEN
                  CALL ZGEMM('C','N',M-ie,nb,mb,DCMPLX(-ONE,ZERO),      &
     &                       A(is,ie+1),Lda,C(is,js),Ldc,               &
     &                       DCMPLX(ONE,ZERO),C(ie+1,js),Ldc)
                  CALL ZGEMM('C','N',M-ie,nb,mb,DCMPLX(-ONE,ZERO),      &
     &                       D(is,ie+1),Ldd,F(is,js),Ldf,               &
     &                       DCMPLX(ONE,ZERO),C(ie+1,js),Ldc)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
      Work(1) = lwmin
!
!
!     End of ZTGSYL
!
      END SUBROUTINE ZTGSYL
