!*==ctrsyl.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CTRSYL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTRSYL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrsyl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrsyl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrsyl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
!                          LDC, SCALE, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANA, TRANB
!       INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
!       REAL               SCALE
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTRSYL solves the complex Sylvester matrix equation:
!>
!>    op(A)*X + X*op(B) = scale*C or
!>    op(A)*X - X*op(B) = scale*C,
!>
!> where op(A) = A or A**H, and A and B are both upper triangular. A is
!> M-by-M and B is N-by-N; the right hand side C and the solution X are
!> M-by-N; and scale is an output scale factor, set <= 1 to avoid
!> overflow in X.
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
!>          = 'C': op(A) = A**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] TRANB
!> \verbatim
!>          TRANB is CHARACTER*1
!>          Specifies the option op(B):
!>          = 'N': op(B) = B    (No transpose)
!>          = 'C': op(B) = B**H (Conjugate transpose)
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
!>          A is COMPLEX array, dimension (LDA,M)
!>          The upper triangular matrix A.
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
!>          B is COMPLEX array, dimension (LDB,N)
!>          The upper triangular matrix B.
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
!>          C is COMPLEX array, dimension (LDC,N)
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
!>          SCALE is REAL
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
!> \ingroup complexSYcomputational
!
!  =====================================================================
      SUBROUTINE CTRSYL(Trana,Tranb,Isgn,M,N,A,Lda,B,Ldb,C,Ldc,Scale,   &
     &                  Info)
      USE S_CDOTC
      USE S_CDOTU
      USE S_CLADIV
      USE S_CLANGE
      USE S_CSSCAL
      USE S_LSAME
      USE S_SLABAD
      USE S_SLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--CTRSYL170
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trana
      CHARACTER :: Tranb
      INTEGER , INTENT(IN) :: Isgn
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , INTENT(INOUT) :: Scale
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX :: a11 , suml , sumr , vec , x11
      REAL :: bignum , da11 , db , eps , scaloc , sgn , smin , smlnum
      REAL , DIMENSION(1) :: dum
      INTEGER :: j , k , l
      LOGICAL :: notrna , notrnb
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
      IF ( .NOT.notrna .AND. .NOT.LSAME(Trana,'C') ) THEN
         Info = -1
      ELSEIF ( .NOT.notrnb .AND. .NOT.LSAME(Tranb,'C') ) THEN
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
         CALL XERBLA('CTRSYL',-Info)
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
      eps = SLAMCH('P')
      smlnum = SLAMCH('S')
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
      smlnum = smlnum*REAL(M*N)/eps
      bignum = ONE/smlnum
      smin = MAX(smlnum,eps*CLANGE('M',M,M,A,Lda,dum),                  &
     &       eps*CLANGE('M',N,N,B,Ldb,dum))
      sgn = Isgn
!
      IF ( notrna .AND. notrnb ) THEN
!
!        Solve    A*X + ISGN*X*B = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        bottom-left corner column by column by
!
!            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        Where
!                    M                        L-1
!          R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)].
!                  I=K+1                      J=1
!
         DO l = 1 , N
            DO k = M , 1 , -1
!
               suml = CDOTU(M-k,A(k,MIN(k+1,M)),Lda,C(MIN(k+1,M),l),1)
               sumr = CDOTU(l-1,C(k,1),Ldc,B(1,l),1)
               vec = C(k,l) - (suml+sgn*sumr)
!
               scaloc = ONE
               a11 = A(k,k) + sgn*B(l,l)
               da11 = ABS(REAL(a11)) + ABS(AIMAG(a11))
               IF ( da11<=smin ) THEN
                  a11 = smin
                  da11 = smin
                  Info = 1
               ENDIF
               db = ABS(REAL(vec)) + ABS(AIMAG(vec))
               IF ( da11<ONE .AND. db>ONE ) THEN
                  IF ( db>bignum*da11 ) scaloc = ONE/db
               ENDIF
               x11 = CLADIV(vec*CMPLX(scaloc),a11)
!
               IF ( scaloc/=ONE ) THEN
                  DO j = 1 , N
                     CALL CSSCAL(M,scaloc,C(1,j),1)
                  ENDDO
                  Scale = Scale*scaloc
               ENDIF
               C(k,l) = x11
!
            ENDDO
         ENDDO
!
      ELSEIF ( .NOT.notrna .AND. notrnb ) THEN
!
!        Solve    A**H *X + ISGN*X*B = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        upper-left corner column by column by
!
!            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        Where
!                   K-1                           L-1
!          R(K,L) = SUM [A**H(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]
!                   I=1                           J=1
!
         DO l = 1 , N
            DO k = 1 , M
!
               suml = CDOTC(k-1,A(1,k),1,C(1,l),1)
               sumr = CDOTU(l-1,C(k,1),Ldc,B(1,l),1)
               vec = C(k,l) - (suml+sgn*sumr)
!
               scaloc = ONE
               a11 = CONJG(A(k,k)) + sgn*B(l,l)
               da11 = ABS(REAL(a11)) + ABS(AIMAG(a11))
               IF ( da11<=smin ) THEN
                  a11 = smin
                  da11 = smin
                  Info = 1
               ENDIF
               db = ABS(REAL(vec)) + ABS(AIMAG(vec))
               IF ( da11<ONE .AND. db>ONE ) THEN
                  IF ( db>bignum*da11 ) scaloc = ONE/db
               ENDIF
!
               x11 = CLADIV(vec*CMPLX(scaloc),a11)
!
               IF ( scaloc/=ONE ) THEN
                  DO j = 1 , N
                     CALL CSSCAL(M,scaloc,C(1,j),1)
                  ENDDO
                  Scale = Scale*scaloc
               ENDIF
               C(k,l) = x11
!
            ENDDO
         ENDDO
!
      ELSEIF ( .NOT.notrna .AND. .NOT.notrnb ) THEN
!
!        Solve    A**H*X + ISGN*X*B**H = C.
!
!        The (K,L)th block of X is determined starting from
!        upper-right corner column by column by
!
!            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)
!
!        Where
!                    K-1
!           R(K,L) = SUM [A**H(I,K)*X(I,L)] +
!                    I=1
!                           N
!                     ISGN*SUM [X(K,J)*B**H(L,J)].
!                          J=L+1
!
         DO l = N , 1 , -1
            DO k = 1 , M
!
               suml = CDOTC(k-1,A(1,k),1,C(1,l),1)
               sumr = CDOTC(N-l,C(k,MIN(l+1,N)),Ldc,B(l,MIN(l+1,N)),Ldb)
               vec = C(k,l) - (suml+sgn*CONJG(sumr))
!
               scaloc = ONE
               a11 = CONJG(A(k,k)+sgn*B(l,l))
               da11 = ABS(REAL(a11)) + ABS(AIMAG(a11))
               IF ( da11<=smin ) THEN
                  a11 = smin
                  da11 = smin
                  Info = 1
               ENDIF
               db = ABS(REAL(vec)) + ABS(AIMAG(vec))
               IF ( da11<ONE .AND. db>ONE ) THEN
                  IF ( db>bignum*da11 ) scaloc = ONE/db
               ENDIF
!
               x11 = CLADIV(vec*CMPLX(scaloc),a11)
!
               IF ( scaloc/=ONE ) THEN
                  DO j = 1 , N
                     CALL CSSCAL(M,scaloc,C(1,j),1)
                  ENDDO
                  Scale = Scale*scaloc
               ENDIF
               C(k,l) = x11
!
            ENDDO
         ENDDO
!
      ELSEIF ( notrna .AND. .NOT.notrnb ) THEN
!
!        Solve    A*X + ISGN*X*B**H = C.
!
!        The (K,L)th block of X is determined starting from
!        bottom-left corner column by column by
!
!           A(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)
!
!        Where
!                    M                          N
!          R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B**H(L,J)]
!                  I=K+1                      J=L+1
!
         DO l = N , 1 , -1
            DO k = M , 1 , -1
!
               suml = CDOTU(M-k,A(k,MIN(k+1,M)),Lda,C(MIN(k+1,M),l),1)
               sumr = CDOTC(N-l,C(k,MIN(l+1,N)),Ldc,B(l,MIN(l+1,N)),Ldb)
               vec = C(k,l) - (suml+sgn*CONJG(sumr))
!
               scaloc = ONE
               a11 = A(k,k) + sgn*CONJG(B(l,l))
               da11 = ABS(REAL(a11)) + ABS(AIMAG(a11))
               IF ( da11<=smin ) THEN
                  a11 = smin
                  da11 = smin
                  Info = 1
               ENDIF
               db = ABS(REAL(vec)) + ABS(AIMAG(vec))
               IF ( da11<ONE .AND. db>ONE ) THEN
                  IF ( db>bignum*da11 ) scaloc = ONE/db
               ENDIF
!
               x11 = CLADIV(vec*CMPLX(scaloc),a11)
!
               IF ( scaloc/=ONE ) THEN
                  DO j = 1 , N
                     CALL CSSCAL(M,scaloc,C(1,j),1)
                  ENDDO
                  Scale = Scale*scaloc
               ENDIF
               C(k,l) = x11
!
            ENDDO
         ENDDO
!
      ENDIF
!
!
!     End of CTRSYL
!
      END SUBROUTINE CTRSYL
