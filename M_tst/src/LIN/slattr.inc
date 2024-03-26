!*==slattr.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b slattr
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE slattr( IMAT, UPLO, TRANS, DIAG, ISEED, N, A, LDA, B,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            IMAT, INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               A( LDA, * ), B( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> slattr generates a triangular test matrix.
!> IMAT and UPLO uniquely specify the properties of the test
!> matrix, which is returned in the array A.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IMAT
!> \verbatim
!>          IMAT is INTEGER
!>          An integer key describing which matrix to generate for this
!>          path.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A will be upper or lower
!>          triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies whether the matrix or its transpose will be used.
!>          = 'N':  No transpose
!>          = 'T':  Transpose
!>          = 'C':  Conjugate transpose (= Transpose)
!> \endverbatim
!>
!> \param[out] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          The seed vector for the random number generator (used in
!>          SLATMS).  Modified on exit.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix to be generated.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The triangular matrix A.  If UPLO = 'U', the leading n by n
!>          upper triangular part of the array A contains the upper
!>          triangular matrix, and the strictly lower triangular part of
!>          A is not referenced.  If UPLO = 'L', the leading n by n lower
!>          triangular part of the array A contains the lower triangular
!>          matrix, and the strictly upper triangular part of A is not
!>          referenced.  If DIAG = 'U', the diagonal elements of A are
!>          set so that A(k,k) = k for 1 <= k <= n.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (N)
!>          The right hand side vector, if IMAT > 10.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (3*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE slattr(Imat,Uplo,Trans,Diag,Iseed,N,A,Lda,B,Work,Info)
      IMPLICIT NONE
!*--SLATTR136
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Trans , Uplo
      INTEGER Imat , Info , Lda , N
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      REAL A(Lda,*) , B(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , TWO , ZERO
      PARAMETER (ONE=1.0E+0,TWO=2.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      CHARACTER dist , type
      CHARACTER*3 path
      INTEGER i , iy , j , jcount , kl , ku , mode
      REAL anorm , bignum , bnorm , bscal , c , cndnum , plus1 , plus2 ,&
     &     ra , rb , rexp , s , sfac , smlnum , star1 , texp , tleft ,  &
     &     tscal , ulp , unfl , x , y , z
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ISAMAX
      REAL SLAMCH 
EXTERNAL    LSAME,   ISAMAX,   SLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    SCOPY,   SLABAD,   SLARNV,   SROT
EXTERNAL    SROTG,   SSCAL,   SSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , REAL , SIGN , SQRT
!     ..
!     .. Executable Statements ..
!
      path(1:1) = 'Single precision'
      path(2:3) = 'TR'
      unfl = SLAMCH('Safe minimum')
      ulp = SLAMCH('Epsilon')*SLAMCH('Base')
      smlnum = unfl
      bignum = (ONE-ulp)/smlnum
      CALL SLABAD(smlnum,bignum)
      IF ( (Imat>=7 .AND. Imat<=10) .OR. Imat==18 ) THEN
         Diag = 'U'
      ELSE
         Diag = 'N'
      ENDIF
      Info = 0
!
!     Quick return if N.LE.0.
!
      IF ( N<=0 ) RETURN
!
!     Call slatb4 to set parameters for SLATMS.
!
      upper = LSAME(Uplo,'U')
      IF ( upper ) THEN
         CALL slatb4(path,Imat,N,N,type,kl,ku,anorm,mode,cndnum,dist)
      ELSE
         CALL slatb4(path,-Imat,N,N,type,kl,ku,anorm,mode,cndnum,dist)
      ENDIF
!
!     IMAT <= 6:  Non-unit triangular matrix
!
      IF ( Imat<=6 ) THEN
         CALL SLATMS(N,N,dist,Iseed,type,B,mode,cndnum,anorm,kl,ku,     &
     &               'No packing',A,Lda,Work,Info)
!
!     IMAT > 6:  Unit triangular matrix
!     The diagonal is deliberately set to something other than 1.
!
!     IMAT = 7:  Matrix is the identity
!
      ELSEIF ( Imat==7 ) THEN
         IF ( upper ) THEN
            DO j = 1 , N
               DO i = 1 , j - 1
                  A(i,j) = ZERO
               ENDDO
               A(j,j) = j
            ENDDO
         ELSE
            DO j = 1 , N
               A(j,j) = j
               DO i = j + 1 , N
                  A(i,j) = ZERO
               ENDDO
            ENDDO
         ENDIF
!
!     IMAT > 7:  Non-trivial unit triangular matrix
!
!     Generate a unit triangular matrix T with condition CNDNUM by
!     forming a triangular matrix with known singular values and
!     filling in the zero entries with Givens rotations.
!
      ELSEIF ( Imat<=10 ) THEN
         IF ( upper ) THEN
            DO j = 1 , N
               DO i = 1 , j - 1
                  A(i,j) = ZERO
               ENDDO
               A(j,j) = j
            ENDDO
         ELSE
            DO j = 1 , N
               A(j,j) = j
               DO i = j + 1 , N
                  A(i,j) = ZERO
               ENDDO
            ENDDO
         ENDIF
!
!        Since the trace of a unit triangular matrix is 1, the product
!        of its singular values must be 1.  Let s = sqrt(CNDNUM),
!        x = sqrt(s) - 1/sqrt(s), y = sqrt(2/(n-2))*x, and z = x**2.
!        The following triangular matrix has singular values s, 1, 1,
!        ..., 1, 1/s:
!
!        1  y  y  y  ...  y  y  z
!           1  0  0  ...  0  0  y
!              1  0  ...  0  0  y
!                 .  ...  .  .  .
!                     .   .  .  .
!                         1  0  y
!                            1  y
!                               1
!
!        To fill in the zeros, we first multiply by a matrix with small
!        condition number of the form
!
!        1  0  0  0  0  ...
!           1  +  *  0  0  ...
!              1  +  0  0  0
!                 1  +  *  0  0
!                    1  +  0  0
!                       ...
!                          1  +  0
!                             1  0
!                                1
!
!        Each element marked with a '*' is formed by taking the product
!        of the adjacent elements marked with '+'.  The '*'s can be
!        chosen freely, and the '+'s are chosen so that the inverse of
!        T will have elements of the same magnitude as T.  If the *'s in
!        both T and inv(T) have small magnitude, T is well conditioned.
!        The two offdiagonals of T are stored in WORK.
!
!        The product of these two matrices has the form
!
!        1  y  y  y  y  y  .  y  y  z
!           1  +  *  0  0  .  0  0  y
!              1  +  0  0  .  0  0  y
!                 1  +  *  .  .  .  .
!                    1  +  .  .  .  .
!                       .  .  .  .  .
!                          .  .  .  .
!                             1  +  y
!                                1  y
!                                   1
!
!        Now we multiply by Givens rotations, using the fact that
!
!              [  c   s ] [  1   w ] [ -c  -s ] =  [  1  -w ]
!              [ -s   c ] [  0   1 ] [  s  -c ]    [  0   1 ]
!        and
!              [ -c  -s ] [  1   0 ] [  c   s ] =  [  1   0 ]
!              [  s  -c ] [  w   1 ] [ -s   c ]    [ -w   1 ]
!
!        where c = w / sqrt(w**2+4) and s = 2 / sqrt(w**2+4).
!
         star1 = 0.25
         sfac = 0.5
         plus1 = sfac
         DO j = 1 , N , 2
            plus2 = star1/plus1
            Work(j) = plus1
            Work(N+j) = star1
            IF ( j+1<=N ) THEN
               Work(j+1) = plus2
               Work(N+j+1) = ZERO
               plus1 = star1/plus2
               rexp = SLARND(2,Iseed)
               star1 = star1*(sfac**rexp)
               IF ( rexp<ZERO ) THEN
                  star1 = -sfac**(ONE-rexp)
               ELSE
                  star1 = sfac**(ONE+rexp)
               ENDIF
            ENDIF
         ENDDO
!
         x = SQRT(cndnum) - 1/SQRT(cndnum)
         IF ( N>2 ) THEN
            y = SQRT(2./(N-2))*x
         ELSE
            y = ZERO
         ENDIF
         z = x*x
!
         IF ( upper ) THEN
            IF ( N>3 ) THEN
               CALL SCOPY(N-3,Work,1,A(2,3),Lda+1)
               IF ( N>4 ) CALL SCOPY(N-4,Work(N+1),1,A(2,4),Lda+1)
            ENDIF
            DO j = 2 , N - 1
               A(1,j) = y
               A(j,N) = y
            ENDDO
            A(1,N) = z
         ELSE
            IF ( N>3 ) THEN
               CALL SCOPY(N-3,Work,1,A(3,2),Lda+1)
               IF ( N>4 ) CALL SCOPY(N-4,Work(N+1),1,A(4,2),Lda+1)
            ENDIF
            DO j = 2 , N - 1
               A(j,1) = y
               A(N,j) = y
            ENDDO
            A(N,1) = z
         ENDIF
!
!        Fill in the zeros using Givens rotations.
!
         IF ( upper ) THEN
            DO j = 1 , N - 1
               ra = A(j,j+1)
               rb = 2.0
               CALL SROTG(ra,rb,c,s)
!
!              Multiply by [ c  s; -s  c] on the left.
!
               IF ( N>j+1 ) CALL SROT(N-j-1,A(j,j+2),Lda,A(j+1,j+2),Lda,&
     &                                c,s)
!
!              Multiply by [-c -s;  s -c] on the right.
!
               IF ( j>1 ) CALL SROT(j-1,A(1,j+1),1,A(1,j),1,-c,-s)
!
!              Negate A(J,J+1).
!
               A(j,j+1) = -A(j,j+1)
            ENDDO
         ELSE
            DO j = 1 , N - 1
               ra = A(j+1,j)
               rb = 2.0
               CALL SROTG(ra,rb,c,s)
!
!              Multiply by [ c -s;  s  c] on the right.
!
               IF ( N>j+1 ) CALL SROT(N-j-1,A(j+2,j+1),1,A(j+2,j),1,c,  &
     &                                -s)
!
!              Multiply by [-c  s; -s -c] on the left.
!
               IF ( j>1 ) CALL SROT(j-1,A(j,1),Lda,A(j+1,1),Lda,-c,s)
!
!              Negate A(J+1,J).
!
               A(j+1,j) = -A(j+1,j)
            ENDDO
         ENDIF
!
!     IMAT > 10:  Pathological test cases.  These triangular matrices
!     are badly scaled or badly conditioned, so when used in solving a
!     triangular system they may cause overflow in the solution vector.
!
      ELSEIF ( Imat==11 ) THEN
!
!        Type 11:  Generate a triangular matrix with elements between
!        -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
!        Make the right hand side large so that it requires scaling.
!
         IF ( upper ) THEN
            DO j = 1 , N
               CALL SLARNV(2,Iseed,j,A(1,j))
               A(j,j) = SIGN(TWO,A(j,j))
            ENDDO
         ELSE
            DO j = 1 , N
               CALL SLARNV(2,Iseed,N-j+1,A(j,j))
               A(j,j) = SIGN(TWO,A(j,j))
            ENDDO
         ENDIF
!
!        Set the right hand side so that the largest value is BIGNUM.
!
         CALL SLARNV(2,Iseed,N,B)
         iy = ISAMAX(N,B,1)
         bnorm = ABS(B(iy))
         bscal = bignum/MAX(ONE,bnorm)
         CALL SSCAL(N,bscal,B,1)
!
      ELSEIF ( Imat==12 ) THEN
!
!        Type 12:  Make the first diagonal element in the solve small to
!        cause immediate overflow when dividing by T(j,j).
!        In type 12, the offdiagonal elements are small (CNORM(j) < 1).
!
         CALL SLARNV(2,Iseed,N,B)
         tscal = ONE/MAX(ONE,REAL(N-1))
         IF ( upper ) THEN
            DO j = 1 , N
               CALL SLARNV(2,Iseed,j,A(1,j))
               CALL SSCAL(j-1,tscal,A(1,j),1)
               A(j,j) = SIGN(ONE,A(j,j))
            ENDDO
            A(N,N) = smlnum*A(N,N)
         ELSE
            DO j = 1 , N
               CALL SLARNV(2,Iseed,N-j+1,A(j,j))
               IF ( N>j ) CALL SSCAL(N-j,tscal,A(j+1,j),1)
               A(j,j) = SIGN(ONE,A(j,j))
            ENDDO
            A(1,1) = smlnum*A(1,1)
         ENDIF
!
      ELSEIF ( Imat==13 ) THEN
!
!        Type 13:  Make the first diagonal element in the solve small to
!        cause immediate overflow when dividing by T(j,j).
!        In type 13, the offdiagonal elements are O(1) (CNORM(j) > 1).
!
         CALL SLARNV(2,Iseed,N,B)
         IF ( upper ) THEN
            DO j = 1 , N
               CALL SLARNV(2,Iseed,j,A(1,j))
               A(j,j) = SIGN(ONE,A(j,j))
            ENDDO
            A(N,N) = smlnum*A(N,N)
         ELSE
            DO j = 1 , N
               CALL SLARNV(2,Iseed,N-j+1,A(j,j))
               A(j,j) = SIGN(ONE,A(j,j))
            ENDDO
            A(1,1) = smlnum*A(1,1)
         ENDIF
!
      ELSEIF ( Imat==14 ) THEN
!
!        Type 14:  T is diagonal with small numbers on the diagonal to
!        make the growth factor underflow, but a small right hand side
!        chosen so that the solution does not overflow.
!
         IF ( upper ) THEN
            jcount = 1
            DO j = N , 1 , -1
               DO i = 1 , j - 1
                  A(i,j) = ZERO
               ENDDO
               IF ( jcount<=2 ) THEN
                  A(j,j) = smlnum
               ELSE
                  A(j,j) = ONE
               ENDIF
               jcount = jcount + 1
               IF ( jcount>4 ) jcount = 1
            ENDDO
         ELSE
            jcount = 1
            DO j = 1 , N
               DO i = j + 1 , N
                  A(i,j) = ZERO
               ENDDO
               IF ( jcount<=2 ) THEN
                  A(j,j) = smlnum
               ELSE
                  A(j,j) = ONE
               ENDIF
               jcount = jcount + 1
               IF ( jcount>4 ) jcount = 1
            ENDDO
         ENDIF
!
!        Set the right hand side alternately zero and small.
!
         IF ( upper ) THEN
            B(1) = ZERO
            DO i = N , 2 , -2
               B(i) = ZERO
               B(i-1) = smlnum
            ENDDO
         ELSE
            B(N) = ZERO
            DO i = 1 , N - 1 , 2
               B(i) = ZERO
               B(i+1) = smlnum
            ENDDO
         ENDIF
!
      ELSEIF ( Imat==15 ) THEN
!
!        Type 15:  Make the diagonal elements small to cause gradual
!        overflow when dividing by T(j,j).  To control the amount of
!        scaling needed, the matrix is bidiagonal.
!
         texp = ONE/MAX(ONE,REAL(N-1))
         tscal = smlnum**texp
         CALL SLARNV(2,Iseed,N,B)
         IF ( upper ) THEN
            DO j = 1 , N
               DO i = 1 , j - 2
                  A(i,j) = 0.
               ENDDO
               IF ( j>1 ) A(j-1,j) = -ONE
               A(j,j) = tscal
            ENDDO
            B(N) = ONE
         ELSE
            DO j = 1 , N
               DO i = j + 2 , N
                  A(i,j) = 0.
               ENDDO
               IF ( j<N ) A(j+1,j) = -ONE
               A(j,j) = tscal
            ENDDO
            B(1) = ONE
         ENDIF
!
      ELSEIF ( Imat==16 ) THEN
!
!        Type 16:  One zero diagonal element.
!
         iy = N/2 + 1
         IF ( upper ) THEN
            DO j = 1 , N
               CALL SLARNV(2,Iseed,j,A(1,j))
               IF ( j/=iy ) THEN
                  A(j,j) = SIGN(TWO,A(j,j))
               ELSE
                  A(j,j) = ZERO
               ENDIF
            ENDDO
         ELSE
            DO j = 1 , N
               CALL SLARNV(2,Iseed,N-j+1,A(j,j))
               IF ( j/=iy ) THEN
                  A(j,j) = SIGN(TWO,A(j,j))
               ELSE
                  A(j,j) = ZERO
               ENDIF
            ENDDO
         ENDIF
         CALL SLARNV(2,Iseed,N,B)
         CALL SSCAL(N,TWO,B,1)
!
      ELSEIF ( Imat==17 ) THEN
!
!        Type 17:  Make the offdiagonal elements large to cause overflow
!        when adding a column of T.  In the non-transposed case, the
!        matrix is constructed to cause overflow when adding a column in
!        every other step.
!
         tscal = unfl/ulp
         tscal = (ONE-ulp)/tscal
         DO j = 1 , N
            DO i = 1 , N
               A(i,j) = 0.
            ENDDO
         ENDDO
         texp = ONE
         IF ( upper ) THEN
            DO j = N , 2 , -2
               A(1,j) = -tscal/REAL(N+1)
               A(j,j) = ONE
               B(j) = texp*(ONE-ulp)
               A(1,j-1) = -(tscal/REAL(N+1))/REAL(N+2)
               A(j-1,j-1) = ONE
               B(j-1) = texp*REAL(N*N+N-1)
               texp = texp*2.
            ENDDO
            B(1) = (REAL(N+1)/REAL(N+2))*tscal
         ELSE
            DO j = 1 , N - 1 , 2
               A(N,j) = -tscal/REAL(N+1)
               A(j,j) = ONE
               B(j) = texp*(ONE-ulp)
               A(N,j+1) = -(tscal/REAL(N+1))/REAL(N+2)
               A(j+1,j+1) = ONE
               B(j+1) = texp*REAL(N*N+N-1)
               texp = texp*2.
            ENDDO
            B(N) = (REAL(N+1)/REAL(N+2))*tscal
         ENDIF
!
      ELSEIF ( Imat==18 ) THEN
!
!        Type 18:  Generate a unit triangular matrix with elements
!        between -1 and 1, and make the right hand side large so that it
!        requires scaling.
!
         IF ( upper ) THEN
            DO j = 1 , N
               CALL SLARNV(2,Iseed,j-1,A(1,j))
               A(j,j) = ZERO
            ENDDO
         ELSE
            DO j = 1 , N
               IF ( j<N ) CALL SLARNV(2,Iseed,N-j,A(j+1,j))
               A(j,j) = ZERO
            ENDDO
         ENDIF
!
!        Set the right hand side so that the largest value is BIGNUM.
!
         CALL SLARNV(2,Iseed,N,B)
         iy = ISAMAX(N,B,1)
         bnorm = ABS(B(iy))
         bscal = bignum/MAX(ONE,bnorm)
         CALL SSCAL(N,bscal,B,1)
!
      ELSEIF ( Imat==19 ) THEN
!
!        Type 19:  Generate a triangular matrix with elements between
!        BIGNUM/(n-1) and BIGNUM so that at least one of the column
!        norms will exceed BIGNUM.
!        1/3/91:  SLATRS no longer can handle this case
!
         tleft = bignum/MAX(ONE,REAL(N-1))
         tscal = bignum*(REAL(N-1)/MAX(ONE,REAL(N)))
         IF ( upper ) THEN
            DO j = 1 , N
               CALL SLARNV(2,Iseed,j,A(1,j))
               DO i = 1 , j
                  A(i,j) = SIGN(tleft,A(i,j)) + tscal*A(i,j)
               ENDDO
            ENDDO
         ELSE
            DO j = 1 , N
               CALL SLARNV(2,Iseed,N-j+1,A(j,j))
               DO i = j , N
                  A(i,j) = SIGN(tleft,A(i,j)) + tscal*A(i,j)
               ENDDO
            ENDDO
         ENDIF
         CALL SLARNV(2,Iseed,N,B)
         CALL SSCAL(N,TWO,B,1)
      ENDIF
!
!     Flip the matrix if the transpose will be used.
!
      IF ( .NOT.LSAME(Trans,'N') ) THEN
         IF ( upper ) THEN
            DO j = 1 , N/2
               CALL SSWAP(N-2*j+1,A(j,j),Lda,A(j+1,N-j+1),-1)
            ENDDO
         ELSE
            DO j = 1 , N/2
               CALL SSWAP(N-2*j+1,A(j,j),1,A(N-j+1,j+1),-Lda)
            ENDDO
         ENDIF
      ENDIF
!
!
!     End of slattr
!
      END SUBROUTINE slattr
