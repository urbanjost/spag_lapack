!*==ztrsyl.f90  processed by SPAG 7.51RB at 20:09 on  3 Mar 2022
!> \brief \b ZTRSYL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTRSYL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrsyl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrsyl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrsyl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
!                          LDC, SCALE, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANA, TRANB
!       INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
!       DOUBLE PRECISION   SCALE
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTRSYL solves the complex Sylvester matrix equation:
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
!>          A is COMPLEX*16 array, dimension (LDA,M)
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
!>          B is COMPLEX*16 array, dimension (LDB,N)
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
!>          C is COMPLEX*16 array, dimension (LDC,N)
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
!> \ingroup complex16SYcomputational
!
!  =====================================================================
      SUBROUTINE ZTRSYL(Trana,Tranb,Isgn,M,N,A,Lda,B,Ldb,C,Ldc,Scale,   &
     &                  Info)
      IMPLICIT NONE
!*--ZTRSYL161
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trana , Tranb
      INTEGER Info , Isgn , Lda , Ldb , Ldc , M , N
      DOUBLE PRECISION Scale
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , C(Ldc,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL notrna , notrnb
      INTEGER j , k , l
      DOUBLE PRECISION bignum , da11 , db , eps , scaloc , sgn , smin , &
     &                 smlnum
      COMPLEX*16 a11 , suml , sumr , vec , x11
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION dum(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , ZLANGE
      COMPLEX*16 ZDOTC , ZDOTU , ZLADIV
      EXTERNAL LSAME , DLAMCH , ZLANGE , ZDOTC , ZDOTU , ZLADIV
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD , XERBLA , ZDSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DCMPLX , DCONJG , DIMAG , MAX , MIN
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
         CALL XERBLA('ZTRSYL',-Info)
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
      smin = MAX(smlnum,eps*ZLANGE('M',M,M,A,Lda,dum),                  &
     &       eps*ZLANGE('M',N,N,B,Ldb,dum))
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
               suml = ZDOTU(M-k,A(k,MIN(k+1,M)),Lda,C(MIN(k+1,M),l),1)
               sumr = ZDOTU(l-1,C(k,1),Ldc,B(1,l),1)
               vec = C(k,l) - (suml+sgn*sumr)
!
               scaloc = ONE
               a11 = A(k,k) + sgn*B(l,l)
               da11 = ABS(DBLE(a11)) + ABS(DIMAG(a11))
               IF ( da11<=smin ) THEN
                  a11 = smin
                  da11 = smin
                  Info = 1
               ENDIF
               db = ABS(DBLE(vec)) + ABS(DIMAG(vec))
               IF ( da11<ONE .AND. db>ONE ) THEN
                  IF ( db>bignum*da11 ) scaloc = ONE/db
               ENDIF
               x11 = ZLADIV(vec*DCMPLX(scaloc),a11)
!
               IF ( scaloc/=ONE ) THEN
                  DO j = 1 , N
                     CALL ZDSCAL(M,scaloc,C(1,j),1)
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
               suml = ZDOTC(k-1,A(1,k),1,C(1,l),1)
               sumr = ZDOTU(l-1,C(k,1),Ldc,B(1,l),1)
               vec = C(k,l) - (suml+sgn*sumr)
!
               scaloc = ONE
               a11 = DCONJG(A(k,k)) + sgn*B(l,l)
               da11 = ABS(DBLE(a11)) + ABS(DIMAG(a11))
               IF ( da11<=smin ) THEN
                  a11 = smin
                  da11 = smin
                  Info = 1
               ENDIF
               db = ABS(DBLE(vec)) + ABS(DIMAG(vec))
               IF ( da11<ONE .AND. db>ONE ) THEN
                  IF ( db>bignum*da11 ) scaloc = ONE/db
               ENDIF
!
               x11 = ZLADIV(vec*DCMPLX(scaloc),a11)
!
               IF ( scaloc/=ONE ) THEN
                  DO j = 1 , N
                     CALL ZDSCAL(M,scaloc,C(1,j),1)
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
               suml = ZDOTC(k-1,A(1,k),1,C(1,l),1)
               sumr = ZDOTC(N-l,C(k,MIN(l+1,N)),Ldc,B(l,MIN(l+1,N)),Ldb)
               vec = C(k,l) - (suml+sgn*DCONJG(sumr))
!
               scaloc = ONE
               a11 = DCONJG(A(k,k)+sgn*B(l,l))
               da11 = ABS(DBLE(a11)) + ABS(DIMAG(a11))
               IF ( da11<=smin ) THEN
                  a11 = smin
                  da11 = smin
                  Info = 1
               ENDIF
               db = ABS(DBLE(vec)) + ABS(DIMAG(vec))
               IF ( da11<ONE .AND. db>ONE ) THEN
                  IF ( db>bignum*da11 ) scaloc = ONE/db
               ENDIF
!
               x11 = ZLADIV(vec*DCMPLX(scaloc),a11)
!
               IF ( scaloc/=ONE ) THEN
                  DO j = 1 , N
                     CALL ZDSCAL(M,scaloc,C(1,j),1)
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
               suml = ZDOTU(M-k,A(k,MIN(k+1,M)),Lda,C(MIN(k+1,M),l),1)
               sumr = ZDOTC(N-l,C(k,MIN(l+1,N)),Ldc,B(l,MIN(l+1,N)),Ldb)
               vec = C(k,l) - (suml+sgn*DCONJG(sumr))
!
               scaloc = ONE
               a11 = A(k,k) + sgn*DCONJG(B(l,l))
               da11 = ABS(DBLE(a11)) + ABS(DIMAG(a11))
               IF ( da11<=smin ) THEN
                  a11 = smin
                  da11 = smin
                  Info = 1
               ENDIF
               db = ABS(DBLE(vec)) + ABS(DIMAG(vec))
               IF ( da11<ONE .AND. db>ONE ) THEN
                  IF ( db>bignum*da11 ) scaloc = ONE/db
               ENDIF
!
               x11 = ZLADIV(vec*DCMPLX(scaloc),a11)
!
               IF ( scaloc/=ONE ) THEN
                  DO j = 1 , N
                     CALL ZDSCAL(M,scaloc,C(1,j),1)
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
!     End of ZTRSYL
!
      END SUBROUTINE ZTRSYL
