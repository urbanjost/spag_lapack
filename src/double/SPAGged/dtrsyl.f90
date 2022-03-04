!*==dtrsyl.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DTRSYL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DTRSYL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrsyl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrsyl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrsyl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
!                          LDC, SCALE, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANA, TRANB
!       INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
!       DOUBLE PRECISION   SCALE
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRSYL solves the real Sylvester matrix equation:
!>
!>    op(A)*X + X*op(B) = scale*C or
!>    op(A)*X - X*op(B) = scale*C,
!>
!> where op(A) = A or A**T, and  A and B are both upper quasi-
!> triangular. A is M-by-M and B is N-by-N; the right hand side C and
!> the solution X are M-by-N; and scale is an output scale factor, set
!> <= 1 to avoid overflow in X.
!>
!> A and B must be in Schur canonical form (as returned by DHSEQR), that
!> is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
!> each 2-by-2 diagonal block has its diagonal elements equal and its
!> off-diagonal elements of opposite sign.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANA
!> \verbatim
!>          TRANA is CHARACTER*1
!>          Specifies the option op(A):
!>          = 'N': op(A) = A    (No transpose)
!>          = 'T': op(A) = A**T (Transpose)
!>          = 'C': op(A) = A**H (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] TRANB
!> \verbatim
!>          TRANB is CHARACTER*1
!>          Specifies the option op(B):
!>          = 'N': op(B) = B    (No transpose)
!>          = 'T': op(B) = B**T (Transpose)
!>          = 'C': op(B) = B**H (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] ISGN
!> \verbatim
!>          ISGN is INTEGER
!>          Specifies the sign in the equation:
!>          = +1: solve op(A)*X + X*op(B) = scale*C
!>          = -1: solve op(A)*X - X*op(B) = scale*C
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrix A, and the number of rows in the
!>          matrices X and C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B, and the number of columns in the
!>          matrices X and C. N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,M)
!>          The upper quasi-triangular matrix A, in Schur canonical form.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          The upper quasi-triangular matrix B, in Schur canonical form.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the M-by-N right hand side matrix C.
!>          On exit, C is overwritten by the solution matrix X.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M)
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          The scale factor, scale, set <= 1 to avoid overflow in X.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          = 1: A and B have common or very close eigenvalues; perturbed
!>               values were used to solve the equation (but the matrices
!>               A and B are unchanged).
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
!> \ingroup doubleSYcomputational
!
!  =====================================================================
      SUBROUTINE DTRSYL(Trana,Tranb,Isgn,M,N,A,Lda,B,Ldb,C,Ldc,Scale,   &
     &                  Info)
      USE F77KINDS                        
      USE S_DDOT
      USE S_DLABAD
      USE S_DLALN2
      USE S_DLAMCH
      USE S_DLANGE
      USE S_DLASY2
      USE S_DSCAL
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DTRSYL178
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trana
      CHARACTER :: Tranb
      INTEGER :: Isgn
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: a11 , bignum , da11 , db , eps , scaloc , sgn ,   &
     &                smin , smlnum , suml , sumr , xnorm
      REAL(R8KIND) , DIMENSION(1) :: dum
      INTEGER :: ierr , j , k , k1 , k2 , knext , l , l1 , l2 , lnext
      LOGICAL :: notrna , notrnb
      REAL(R8KIND) , DIMENSION(2,2) :: vec , x
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
!     Decode and Test input parameters
!
      notrna = LSAME(Trana,'N')
      notrnb = LSAME(Tranb,'N')
!
      Info = 0
      IF ( .NOT.notrna .AND. .NOT.LSAME(Trana,'T') .AND.                &
     &     .NOT.LSAME(Trana,'C') ) THEN
         Info = -1
      ELSEIF ( .NOT.notrnb .AND. .NOT.LSAME(Tranb,'T') .AND.            &
     &         .NOT.LSAME(Tranb,'C') ) THEN
         Info = -2
      ELSEIF ( Isgn/=1 .AND. Isgn/=-1 ) THEN
         Info = -3
      ELSEIF ( M<0 ) THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -7
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -9
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -11
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DTRSYL',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      Scale = ONE
      IF ( M==0 .OR. N==0 ) RETURN
!
!     Set constants to control overflow
!
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
      smlnum = smlnum*DBLE(M*N)/eps
      bignum = ONE/smlnum
!
      smin = MAX(smlnum,eps*DLANGE('M',M,M,A,Lda,dum),                  &
     &       eps*DLANGE('M',N,N,B,Ldb,dum))
!
      sgn = Isgn
!
      IF ( notrna .AND. notrnb ) THEN
!
!        Solve    A*X + ISGN*X*B = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        bottom-left corner column by column by
!
!         A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        Where
!                  M                         L-1
!        R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
!                I=K+1                       J=1
!
!        Start column loop (index = L)
!        L1 (L2) : column index of the first (first) row of X(K,L).
!
         lnext = 1
         DO l = 1 , N
            IF ( l>=lnext ) THEN
               IF ( l==N ) THEN
                  l1 = l
                  l2 = l
               ELSEIF ( B(l+1,l)/=ZERO ) THEN
                  l1 = l
                  l2 = l + 1
                  lnext = l + 2
               ELSE
                  l1 = l
                  l2 = l
                  lnext = l + 1
               ENDIF
!
!           Start row loop (index = K)
!           K1 (K2): row index of the first (last) row of X(K,L).
!
               knext = M
               DO k = M , 1 , -1
                  IF ( k<=knext ) THEN
                     IF ( k==1 ) THEN
                        k1 = k
                        k2 = k
                     ELSEIF ( A(k,k-1)/=ZERO ) THEN
                        k1 = k - 1
                        k2 = k
                        knext = k - 2
                     ELSE
                        k1 = k
                        k2 = k
                        knext = k - 1
                     ENDIF
!
                     IF ( l1==l2 .AND. k1==k2 ) THEN
                        suml = DDOT(M-k1,A(k1,MIN(k1+1,M)),Lda,         &
     &                         C(MIN(k1+1,M),l1),1)
                        sumr = DDOT(l1-1,C(k1,1),Ldc,B(1,l1),1)
                        vec(1,1) = C(k1,l1) - (suml+sgn*sumr)
                        scaloc = ONE
!
                        a11 = A(k1,k1) + sgn*B(l1,l1)
                        da11 = ABS(a11)
                        IF ( da11<=smin ) THEN
                           a11 = smin
                           da11 = smin
                           Info = 1
                        ENDIF
                        db = ABS(vec(1,1))
                        IF ( da11<ONE .AND. db>ONE ) THEN
                           IF ( db>bignum*da11 ) scaloc = ONE/db
                        ENDIF
                        x(1,1) = (vec(1,1)*scaloc)/a11
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
!
                     ELSEIF ( l1==l2 .AND. k1/=k2 ) THEN
!
                        suml = DDOT(M-k2,A(k1,MIN(k2+1,M)),Lda,         &
     &                         C(MIN(k2+1,M),l1),1)
                        sumr = DDOT(l1-1,C(k1,1),Ldc,B(1,l1),1)
                        vec(1,1) = C(k1,l1) - (suml+sgn*sumr)
!
                        suml = DDOT(M-k2,A(k2,MIN(k2+1,M)),Lda,         &
     &                         C(MIN(k2+1,M),l1),1)
                        sumr = DDOT(l1-1,C(k2,1),Ldc,B(1,l1),1)
                        vec(2,1) = C(k2,l1) - (suml+sgn*sumr)
!
                        CALL DLALN2(.FALSE.,2,1,smin,ONE,A(k1,k1),Lda,  &
     &                              ONE,ONE,vec,2,-sgn*B(l1,l1),ZERO,x, &
     &                              2,scaloc,xnorm,ierr)
                        IF ( ierr/=0 ) Info = 1
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
                        C(k2,l1) = x(2,1)
!
                     ELSEIF ( l1/=l2 .AND. k1==k2 ) THEN
!
                        suml = DDOT(M-k1,A(k1,MIN(k1+1,M)),Lda,         &
     &                         C(MIN(k1+1,M),l1),1)
                        sumr = DDOT(l1-1,C(k1,1),Ldc,B(1,l1),1)
                        vec(1,1) = sgn*(C(k1,l1)-(suml+sgn*sumr))
!
                        suml = DDOT(M-k1,A(k1,MIN(k1+1,M)),Lda,         &
     &                         C(MIN(k1+1,M),l2),1)
                        sumr = DDOT(l1-1,C(k1,1),Ldc,B(1,l2),1)
                        vec(2,1) = sgn*(C(k1,l2)-(suml+sgn*sumr))
!
                        CALL DLALN2(.TRUE.,2,1,smin,ONE,B(l1,l1),Ldb,   &
     &                              ONE,ONE,vec,2,-sgn*A(k1,k1),ZERO,x, &
     &                              2,scaloc,xnorm,ierr)
                        IF ( ierr/=0 ) Info = 1
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
                        C(k1,l2) = x(2,1)
!
                     ELSEIF ( l1/=l2 .AND. k1/=k2 ) THEN
!
                        suml = DDOT(M-k2,A(k1,MIN(k2+1,M)),Lda,         &
     &                         C(MIN(k2+1,M),l1),1)
                        sumr = DDOT(l1-1,C(k1,1),Ldc,B(1,l1),1)
                        vec(1,1) = C(k1,l1) - (suml+sgn*sumr)
!
                        suml = DDOT(M-k2,A(k1,MIN(k2+1,M)),Lda,         &
     &                         C(MIN(k2+1,M),l2),1)
                        sumr = DDOT(l1-1,C(k1,1),Ldc,B(1,l2),1)
                        vec(1,2) = C(k1,l2) - (suml+sgn*sumr)
!
                        suml = DDOT(M-k2,A(k2,MIN(k2+1,M)),Lda,         &
     &                         C(MIN(k2+1,M),l1),1)
                        sumr = DDOT(l1-1,C(k2,1),Ldc,B(1,l1),1)
                        vec(2,1) = C(k2,l1) - (suml+sgn*sumr)
!
                        suml = DDOT(M-k2,A(k2,MIN(k2+1,M)),Lda,         &
     &                         C(MIN(k2+1,M),l2),1)
                        sumr = DDOT(l1-1,C(k2,1),Ldc,B(1,l2),1)
                        vec(2,2) = C(k2,l2) - (suml+sgn*sumr)
!
                        CALL DLASY2(.FALSE.,.FALSE.,Isgn,2,2,A(k1,k1),  &
     &                              Lda,B(l1,l1),Ldb,vec,2,scaloc,x,2,  &
     &                              xnorm,ierr)
                        IF ( ierr/=0 ) Info = 1
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
                        C(k1,l2) = x(1,2)
                        C(k2,l1) = x(2,1)
                        C(k2,l2) = x(2,2)
                     ENDIF
                  ENDIF
!
               ENDDO
            ENDIF
!
         ENDDO
!
      ELSEIF ( .NOT.notrna .AND. notrnb ) THEN
!
!        Solve    A**T *X + ISGN*X*B = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        upper-left corner column by column by
!
!          A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        Where
!                   K-1        T                    L-1
!          R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
!                   I=1                          J=1
!
!        Start column loop (index = L)
!        L1 (L2): column index of the first (last) row of X(K,L)
!
         lnext = 1
         DO l = 1 , N
            IF ( l>=lnext ) THEN
               IF ( l==N ) THEN
                  l1 = l
                  l2 = l
               ELSEIF ( B(l+1,l)/=ZERO ) THEN
                  l1 = l
                  l2 = l + 1
                  lnext = l + 2
               ELSE
                  l1 = l
                  l2 = l
                  lnext = l + 1
               ENDIF
!
!           Start row loop (index = K)
!           K1 (K2): row index of the first (last) row of X(K,L)
!
               knext = 1
               DO k = 1 , M
                  IF ( k>=knext ) THEN
                     IF ( k==M ) THEN
                        k1 = k
                        k2 = k
                     ELSEIF ( A(k+1,k)/=ZERO ) THEN
                        k1 = k
                        k2 = k + 1
                        knext = k + 2
                     ELSE
                        k1 = k
                        k2 = k
                        knext = k + 1
                     ENDIF
!
                     IF ( l1==l2 .AND. k1==k2 ) THEN
                        suml = DDOT(k1-1,A(1,k1),1,C(1,l1),1)
                        sumr = DDOT(l1-1,C(k1,1),Ldc,B(1,l1),1)
                        vec(1,1) = C(k1,l1) - (suml+sgn*sumr)
                        scaloc = ONE
!
                        a11 = A(k1,k1) + sgn*B(l1,l1)
                        da11 = ABS(a11)
                        IF ( da11<=smin ) THEN
                           a11 = smin
                           da11 = smin
                           Info = 1
                        ENDIF
                        db = ABS(vec(1,1))
                        IF ( da11<ONE .AND. db>ONE ) THEN
                           IF ( db>bignum*da11 ) scaloc = ONE/db
                        ENDIF
                        x(1,1) = (vec(1,1)*scaloc)/a11
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
!
                     ELSEIF ( l1==l2 .AND. k1/=k2 ) THEN
!
                        suml = DDOT(k1-1,A(1,k1),1,C(1,l1),1)
                        sumr = DDOT(l1-1,C(k1,1),Ldc,B(1,l1),1)
                        vec(1,1) = C(k1,l1) - (suml+sgn*sumr)
!
                        suml = DDOT(k1-1,A(1,k2),1,C(1,l1),1)
                        sumr = DDOT(l1-1,C(k2,1),Ldc,B(1,l1),1)
                        vec(2,1) = C(k2,l1) - (suml+sgn*sumr)
!
                        CALL DLALN2(.TRUE.,2,1,smin,ONE,A(k1,k1),Lda,   &
     &                              ONE,ONE,vec,2,-sgn*B(l1,l1),ZERO,x, &
     &                              2,scaloc,xnorm,ierr)
                        IF ( ierr/=0 ) Info = 1
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
                        C(k2,l1) = x(2,1)
!
                     ELSEIF ( l1/=l2 .AND. k1==k2 ) THEN
!
                        suml = DDOT(k1-1,A(1,k1),1,C(1,l1),1)
                        sumr = DDOT(l1-1,C(k1,1),Ldc,B(1,l1),1)
                        vec(1,1) = sgn*(C(k1,l1)-(suml+sgn*sumr))
!
                        suml = DDOT(k1-1,A(1,k1),1,C(1,l2),1)
                        sumr = DDOT(l1-1,C(k1,1),Ldc,B(1,l2),1)
                        vec(2,1) = sgn*(C(k1,l2)-(suml+sgn*sumr))
!
                        CALL DLALN2(.TRUE.,2,1,smin,ONE,B(l1,l1),Ldb,   &
     &                              ONE,ONE,vec,2,-sgn*A(k1,k1),ZERO,x, &
     &                              2,scaloc,xnorm,ierr)
                        IF ( ierr/=0 ) Info = 1
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
                        C(k1,l2) = x(2,1)
!
                     ELSEIF ( l1/=l2 .AND. k1/=k2 ) THEN
!
                        suml = DDOT(k1-1,A(1,k1),1,C(1,l1),1)
                        sumr = DDOT(l1-1,C(k1,1),Ldc,B(1,l1),1)
                        vec(1,1) = C(k1,l1) - (suml+sgn*sumr)
!
                        suml = DDOT(k1-1,A(1,k1),1,C(1,l2),1)
                        sumr = DDOT(l1-1,C(k1,1),Ldc,B(1,l2),1)
                        vec(1,2) = C(k1,l2) - (suml+sgn*sumr)
!
                        suml = DDOT(k1-1,A(1,k2),1,C(1,l1),1)
                        sumr = DDOT(l1-1,C(k2,1),Ldc,B(1,l1),1)
                        vec(2,1) = C(k2,l1) - (suml+sgn*sumr)
!
                        suml = DDOT(k1-1,A(1,k2),1,C(1,l2),1)
                        sumr = DDOT(l1-1,C(k2,1),Ldc,B(1,l2),1)
                        vec(2,2) = C(k2,l2) - (suml+sgn*sumr)
!
                        CALL DLASY2(.TRUE.,.FALSE.,Isgn,2,2,A(k1,k1),   &
     &                              Lda,B(l1,l1),Ldb,vec,2,scaloc,x,2,  &
     &                              xnorm,ierr)
                        IF ( ierr/=0 ) Info = 1
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
                        C(k1,l2) = x(1,2)
                        C(k2,l1) = x(2,1)
                        C(k2,l2) = x(2,2)
                     ENDIF
                  ENDIF
!
               ENDDO
            ENDIF
         ENDDO
!
      ELSEIF ( .NOT.notrna .AND. .NOT.notrnb ) THEN
!
!        Solve    A**T*X + ISGN*X*B**T = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        top-right corner column by column by
!
!           A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)
!
!        Where
!                     K-1                            N
!            R(K,L) = SUM [A(I,K)**T*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
!                     I=1                          J=L+1
!
!        Start column loop (index = L)
!        L1 (L2): column index of the first (last) row of X(K,L)
!
         lnext = N
         DO l = N , 1 , -1
            IF ( l<=lnext ) THEN
               IF ( l==1 ) THEN
                  l1 = l
                  l2 = l
               ELSEIF ( B(l,l-1)/=ZERO ) THEN
                  l1 = l - 1
                  l2 = l
                  lnext = l - 2
               ELSE
                  l1 = l
                  l2 = l
                  lnext = l - 1
               ENDIF
!
!           Start row loop (index = K)
!           K1 (K2): row index of the first (last) row of X(K,L)
!
               knext = 1
               DO k = 1 , M
                  IF ( k>=knext ) THEN
                     IF ( k==M ) THEN
                        k1 = k
                        k2 = k
                     ELSEIF ( A(k+1,k)/=ZERO ) THEN
                        k1 = k
                        k2 = k + 1
                        knext = k + 2
                     ELSE
                        k1 = k
                        k2 = k
                        knext = k + 1
                     ENDIF
!
                     IF ( l1==l2 .AND. k1==k2 ) THEN
                        suml = DDOT(k1-1,A(1,k1),1,C(1,l1),1)
                        sumr = DDOT(N-l1,C(k1,MIN(l1+1,N)),Ldc,         &
     &                         B(l1,MIN(l1+1,N)),Ldb)
                        vec(1,1) = C(k1,l1) - (suml+sgn*sumr)
                        scaloc = ONE
!
                        a11 = A(k1,k1) + sgn*B(l1,l1)
                        da11 = ABS(a11)
                        IF ( da11<=smin ) THEN
                           a11 = smin
                           da11 = smin
                           Info = 1
                        ENDIF
                        db = ABS(vec(1,1))
                        IF ( da11<ONE .AND. db>ONE ) THEN
                           IF ( db>bignum*da11 ) scaloc = ONE/db
                        ENDIF
                        x(1,1) = (vec(1,1)*scaloc)/a11
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
!
                     ELSEIF ( l1==l2 .AND. k1/=k2 ) THEN
!
                        suml = DDOT(k1-1,A(1,k1),1,C(1,l1),1)
                        sumr = DDOT(N-l2,C(k1,MIN(l2+1,N)),Ldc,         &
     &                         B(l1,MIN(l2+1,N)),Ldb)
                        vec(1,1) = C(k1,l1) - (suml+sgn*sumr)
!
                        suml = DDOT(k1-1,A(1,k2),1,C(1,l1),1)
                        sumr = DDOT(N-l2,C(k2,MIN(l2+1,N)),Ldc,         &
     &                         B(l1,MIN(l2+1,N)),Ldb)
                        vec(2,1) = C(k2,l1) - (suml+sgn*sumr)
!
                        CALL DLALN2(.TRUE.,2,1,smin,ONE,A(k1,k1),Lda,   &
     &                              ONE,ONE,vec,2,-sgn*B(l1,l1),ZERO,x, &
     &                              2,scaloc,xnorm,ierr)
                        IF ( ierr/=0 ) Info = 1
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
                        C(k2,l1) = x(2,1)
!
                     ELSEIF ( l1/=l2 .AND. k1==k2 ) THEN
!
                        suml = DDOT(k1-1,A(1,k1),1,C(1,l1),1)
                        sumr = DDOT(N-l2,C(k1,MIN(l2+1,N)),Ldc,         &
     &                         B(l1,MIN(l2+1,N)),Ldb)
                        vec(1,1) = sgn*(C(k1,l1)-(suml+sgn*sumr))
!
                        suml = DDOT(k1-1,A(1,k1),1,C(1,l2),1)
                        sumr = DDOT(N-l2,C(k1,MIN(l2+1,N)),Ldc,         &
     &                         B(l2,MIN(l2+1,N)),Ldb)
                        vec(2,1) = sgn*(C(k1,l2)-(suml+sgn*sumr))
!
                        CALL DLALN2(.FALSE.,2,1,smin,ONE,B(l1,l1),Ldb,  &
     &                              ONE,ONE,vec,2,-sgn*A(k1,k1),ZERO,x, &
     &                              2,scaloc,xnorm,ierr)
                        IF ( ierr/=0 ) Info = 1
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
                        C(k1,l2) = x(2,1)
!
                     ELSEIF ( l1/=l2 .AND. k1/=k2 ) THEN
!
                        suml = DDOT(k1-1,A(1,k1),1,C(1,l1),1)
                        sumr = DDOT(N-l2,C(k1,MIN(l2+1,N)),Ldc,         &
     &                         B(l1,MIN(l2+1,N)),Ldb)
                        vec(1,1) = C(k1,l1) - (suml+sgn*sumr)
!
                        suml = DDOT(k1-1,A(1,k1),1,C(1,l2),1)
                        sumr = DDOT(N-l2,C(k1,MIN(l2+1,N)),Ldc,         &
     &                         B(l2,MIN(l2+1,N)),Ldb)
                        vec(1,2) = C(k1,l2) - (suml+sgn*sumr)
!
                        suml = DDOT(k1-1,A(1,k2),1,C(1,l1),1)
                        sumr = DDOT(N-l2,C(k2,MIN(l2+1,N)),Ldc,         &
     &                         B(l1,MIN(l2+1,N)),Ldb)
                        vec(2,1) = C(k2,l1) - (suml+sgn*sumr)
!
                        suml = DDOT(k1-1,A(1,k2),1,C(1,l2),1)
                        sumr = DDOT(N-l2,C(k2,MIN(l2+1,N)),Ldc,         &
     &                         B(l2,MIN(l2+1,N)),Ldb)
                        vec(2,2) = C(k2,l2) - (suml+sgn*sumr)
!
                        CALL DLASY2(.TRUE.,.TRUE.,Isgn,2,2,A(k1,k1),Lda,&
     &                              B(l1,l1),Ldb,vec,2,scaloc,x,2,xnorm,&
     &                              ierr)
                        IF ( ierr/=0 ) Info = 1
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
                        C(k1,l2) = x(1,2)
                        C(k2,l1) = x(2,1)
                        C(k2,l2) = x(2,2)
                     ENDIF
                  ENDIF
!
               ENDDO
            ENDIF
         ENDDO
!
      ELSEIF ( notrna .AND. .NOT.notrnb ) THEN
!
!        Solve    A*X + ISGN*X*B**T = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        bottom-right corner column by column by
!
!            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)
!
!        Where
!                      M                          N
!            R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
!                    I=K+1                      J=L+1
!
!        Start column loop (index = L)
!        L1 (L2): column index of the first (last) row of X(K,L)
!
         lnext = N
         DO l = N , 1 , -1
            IF ( l<=lnext ) THEN
               IF ( l==1 ) THEN
                  l1 = l
                  l2 = l
               ELSEIF ( B(l,l-1)/=ZERO ) THEN
                  l1 = l - 1
                  l2 = l
                  lnext = l - 2
               ELSE
                  l1 = l
                  l2 = l
                  lnext = l - 1
               ENDIF
!
!           Start row loop (index = K)
!           K1 (K2): row index of the first (last) row of X(K,L)
!
               knext = M
               DO k = M , 1 , -1
                  IF ( k<=knext ) THEN
                     IF ( k==1 ) THEN
                        k1 = k
                        k2 = k
                     ELSEIF ( A(k,k-1)/=ZERO ) THEN
                        k1 = k - 1
                        k2 = k
                        knext = k - 2
                     ELSE
                        k1 = k
                        k2 = k
                        knext = k - 1
                     ENDIF
!
                     IF ( l1==l2 .AND. k1==k2 ) THEN
                        suml = DDOT(M-k1,A(k1,MIN(k1+1,M)),Lda,         &
     &                         C(MIN(k1+1,M),l1),1)
                        sumr = DDOT(N-l1,C(k1,MIN(l1+1,N)),Ldc,         &
     &                         B(l1,MIN(l1+1,N)),Ldb)
                        vec(1,1) = C(k1,l1) - (suml+sgn*sumr)
                        scaloc = ONE
!
                        a11 = A(k1,k1) + sgn*B(l1,l1)
                        da11 = ABS(a11)
                        IF ( da11<=smin ) THEN
                           a11 = smin
                           da11 = smin
                           Info = 1
                        ENDIF
                        db = ABS(vec(1,1))
                        IF ( da11<ONE .AND. db>ONE ) THEN
                           IF ( db>bignum*da11 ) scaloc = ONE/db
                        ENDIF
                        x(1,1) = (vec(1,1)*scaloc)/a11
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
!
                     ELSEIF ( l1==l2 .AND. k1/=k2 ) THEN
!
                        suml = DDOT(M-k2,A(k1,MIN(k2+1,M)),Lda,         &
     &                         C(MIN(k2+1,M),l1),1)
                        sumr = DDOT(N-l2,C(k1,MIN(l2+1,N)),Ldc,         &
     &                         B(l1,MIN(l2+1,N)),Ldb)
                        vec(1,1) = C(k1,l1) - (suml+sgn*sumr)
!
                        suml = DDOT(M-k2,A(k2,MIN(k2+1,M)),Lda,         &
     &                         C(MIN(k2+1,M),l1),1)
                        sumr = DDOT(N-l2,C(k2,MIN(l2+1,N)),Ldc,         &
     &                         B(l1,MIN(l2+1,N)),Ldb)
                        vec(2,1) = C(k2,l1) - (suml+sgn*sumr)
!
                        CALL DLALN2(.FALSE.,2,1,smin,ONE,A(k1,k1),Lda,  &
     &                              ONE,ONE,vec,2,-sgn*B(l1,l1),ZERO,x, &
     &                              2,scaloc,xnorm,ierr)
                        IF ( ierr/=0 ) Info = 1
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
                        C(k2,l1) = x(2,1)
!
                     ELSEIF ( l1/=l2 .AND. k1==k2 ) THEN
!
                        suml = DDOT(M-k1,A(k1,MIN(k1+1,M)),Lda,         &
     &                         C(MIN(k1+1,M),l1),1)
                        sumr = DDOT(N-l2,C(k1,MIN(l2+1,N)),Ldc,         &
     &                         B(l1,MIN(l2+1,N)),Ldb)
                        vec(1,1) = sgn*(C(k1,l1)-(suml+sgn*sumr))
!
                        suml = DDOT(M-k1,A(k1,MIN(k1+1,M)),Lda,         &
     &                         C(MIN(k1+1,M),l2),1)
                        sumr = DDOT(N-l2,C(k1,MIN(l2+1,N)),Ldc,         &
     &                         B(l2,MIN(l2+1,N)),Ldb)
                        vec(2,1) = sgn*(C(k1,l2)-(suml+sgn*sumr))
!
                        CALL DLALN2(.FALSE.,2,1,smin,ONE,B(l1,l1),Ldb,  &
     &                              ONE,ONE,vec,2,-sgn*A(k1,k1),ZERO,x, &
     &                              2,scaloc,xnorm,ierr)
                        IF ( ierr/=0 ) Info = 1
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
                        C(k1,l2) = x(2,1)
!
                     ELSEIF ( l1/=l2 .AND. k1/=k2 ) THEN
!
                        suml = DDOT(M-k2,A(k1,MIN(k2+1,M)),Lda,         &
     &                         C(MIN(k2+1,M),l1),1)
                        sumr = DDOT(N-l2,C(k1,MIN(l2+1,N)),Ldc,         &
     &                         B(l1,MIN(l2+1,N)),Ldb)
                        vec(1,1) = C(k1,l1) - (suml+sgn*sumr)
!
                        suml = DDOT(M-k2,A(k1,MIN(k2+1,M)),Lda,         &
     &                         C(MIN(k2+1,M),l2),1)
                        sumr = DDOT(N-l2,C(k1,MIN(l2+1,N)),Ldc,         &
     &                         B(l2,MIN(l2+1,N)),Ldb)
                        vec(1,2) = C(k1,l2) - (suml+sgn*sumr)
!
                        suml = DDOT(M-k2,A(k2,MIN(k2+1,M)),Lda,         &
     &                         C(MIN(k2+1,M),l1),1)
                        sumr = DDOT(N-l2,C(k2,MIN(l2+1,N)),Ldc,         &
     &                         B(l1,MIN(l2+1,N)),Ldb)
                        vec(2,1) = C(k2,l1) - (suml+sgn*sumr)
!
                        suml = DDOT(M-k2,A(k2,MIN(k2+1,M)),Lda,         &
     &                         C(MIN(k2+1,M),l2),1)
                        sumr = DDOT(N-l2,C(k2,MIN(l2+1,N)),Ldc,         &
     &                         B(l2,MIN(l2+1,N)),Ldb)
                        vec(2,2) = C(k2,l2) - (suml+sgn*sumr)
!
                        CALL DLASY2(.FALSE.,.TRUE.,Isgn,2,2,A(k1,k1),   &
     &                              Lda,B(l1,l1),Ldb,vec,2,scaloc,x,2,  &
     &                              xnorm,ierr)
                        IF ( ierr/=0 ) Info = 1
!
                        IF ( scaloc/=ONE ) THEN
                           DO j = 1 , N
                              CALL DSCAL(M,scaloc,C(1,j),1)
                           ENDDO
                           Scale = Scale*scaloc
                        ENDIF
                        C(k1,l1) = x(1,1)
                        C(k1,l2) = x(1,2)
                        C(k2,l1) = x(2,1)
                        C(k2,l2) = x(2,2)
                     ENDIF
                  ENDIF
!
               ENDDO
            ENDIF
         ENDDO
!
      ENDIF
!
!
!     End of DTRSYL
!
      END SUBROUTINE DTRSYL
