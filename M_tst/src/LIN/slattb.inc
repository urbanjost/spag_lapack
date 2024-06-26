!*==slattb.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b slattb
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE slattb( IMAT, UPLO, TRANS, DIAG, ISEED, N, KD, AB,
!                          LDAB, B, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            IMAT, INFO, KD, LDAB, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               AB( LDAB, * ), B( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> slattb generates a triangular test matrix in 2-dimensional storage.
!> IMAT and UPLO uniquely specify the properties of the test matrix,
!> which is returned in the array A.
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
!>          = 'C':  Conjugate transpose (= transpose)
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
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals or subdiagonals of the banded
!>          triangular matrix A.  KD >= 0.
!> \endverbatim
!>
!> \param[out] AB
!> \verbatim
!>          AB is REAL array, dimension (LDAB,N)
!>          The upper or lower triangular banded matrix A, stored in the
!>          first KD+1 rows of AB.  Let j be a column of A, 1<=j<=n.
!>          If UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j.
!>          If UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (2*N)
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
      SUBROUTINE slattb(Imat,Uplo,Trans,Diag,Iseed,N,Kd,Ab,Ldab,B,Work, &
     &                  Info)
      IMPLICIT NONE
!*--SLATTB139
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Trans , Uplo
      INTEGER Imat , Info , Kd , Ldab , N
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      REAL Ab(Ldab,*) , B(*) , Work(*)
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
      CHARACTER dist , packit , type
      CHARACTER*3 path
      INTEGER i , ioff , iy , j , jcount , kl , ku , lenj , mode
      REAL anorm , bignum , bnorm , bscal , cndnum , plus1 , plus2 ,    &
     &     rexp , sfac , smlnum , star1 , texp , tleft , tnorm , tscal ,&
     &     ulp , unfl
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ISAMAX
      REAL SLAMCH 
EXTERNAL    LSAME,   ISAMAX,   SLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    SCOPY,   SLABAD,   SLARNV,   SSCAL,   SSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , REAL , SIGN , SQRT
!     ..
!     .. Executable Statements ..
!
      path(1:1) = 'Single precision'
      path(2:3) = 'TB'
      unfl = SLAMCH('Safe minimum')
      ulp = SLAMCH('Epsilon')*SLAMCH('Base')
      smlnum = unfl
      bignum = (ONE-ulp)/smlnum
      CALL SLABAD(smlnum,bignum)
      IF ( (Imat>=6 .AND. Imat<=9) .OR. Imat==17 ) THEN
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
         ku = Kd
         ioff = 1 + MAX(0,Kd-N+1)
         kl = 0
         packit = 'Q'
      ELSE
         CALL slatb4(path,-Imat,N,N,type,kl,ku,anorm,mode,cndnum,dist)
         kl = Kd
         ioff = 1
         ku = 0
         packit = 'B'
      ENDIF
!
!     IMAT <= 5:  Non-unit triangular matrix
!
      IF ( Imat<=5 ) THEN
         CALL SLATMS(N,N,dist,Iseed,type,B,mode,cndnum,anorm,kl,ku,     &
     &               packit,Ab(ioff,1),Ldab,Work,Info)
!
!     IMAT > 5:  Unit triangular matrix
!     The diagonal is deliberately set to something other than 1.
!
!     IMAT = 6:  Matrix is the identity
!
      ELSEIF ( Imat==6 ) THEN
         IF ( upper ) THEN
            DO j = 1 , N
               DO i = MAX(1,Kd+2-j) , Kd
                  Ab(i,j) = ZERO
               ENDDO
               Ab(Kd+1,j) = j
            ENDDO
         ELSE
            DO j = 1 , N
               Ab(1,j) = j
               DO i = 2 , MIN(Kd+1,N-j+1)
                  Ab(i,j) = ZERO
               ENDDO
            ENDDO
         ENDIF
!
!     IMAT > 6:  Non-trivial unit triangular matrix
!
!     A unit triangular matrix T with condition CNDNUM is formed.
!     In this version, T only has bandwidth 2, the rest of it is zero.
!
      ELSEIF ( Imat<=9 ) THEN
         tnorm = SQRT(cndnum)
!
!        Initialize AB to zero.
!
         IF ( upper ) THEN
            DO j = 1 , N
               DO i = MAX(1,Kd+2-j) , Kd
                  Ab(i,j) = ZERO
               ENDDO
               Ab(Kd+1,j) = REAL(j)
            ENDDO
         ELSE
            DO j = 1 , N
               DO i = 2 , MIN(Kd+1,N-j+1)
                  Ab(i,j) = ZERO
               ENDDO
               Ab(1,j) = REAL(j)
            ENDDO
         ENDIF
!
!        Special case:  T is tridiagonal.  Set every other offdiagonal
!        so that the matrix has norm TNORM+1.
!
         IF ( Kd==1 ) THEN
            IF ( upper ) THEN
               Ab(1,2) = SIGN(tnorm,SLARND(2,Iseed))
               lenj = (N-3)/2
               CALL SLARNV(2,Iseed,lenj,Work)
               DO j = 1 , lenj
                  Ab(1,2*(j+1)) = tnorm*Work(j)
               ENDDO
            ELSE
               Ab(2,1) = SIGN(tnorm,SLARND(2,Iseed))
               lenj = (N-3)/2
               CALL SLARNV(2,Iseed,lenj,Work)
               DO j = 1 , lenj
                  Ab(2,2*j+1) = tnorm*Work(j)
               ENDDO
            ENDIF
         ELSEIF ( Kd>1 ) THEN
!
!           Form a unit triangular matrix T with condition CNDNUM.  T is
!           given by
!                   | 1   +   *                      |
!                   |     1   +                      |
!               T = |         1   +   *              |
!                   |             1   +              |
!                   |                 1   +   *      |
!                   |                     1   +      |
!                   |                          . . . |
!        Each element marked with a '*' is formed by taking the product
!        of the adjacent elements marked with '+'.  The '*'s can be
!        chosen freely, and the '+'s are chosen so that the inverse of
!        T will have elements of the same magnitude as T.
!
!        The two offdiagonals of T are stored in WORK.
!
            star1 = SIGN(tnorm,SLARND(2,Iseed))
            sfac = SQRT(tnorm)
            plus1 = SIGN(sfac,SLARND(2,Iseed))
            DO j = 1 , N , 2
               plus2 = star1/plus1
               Work(j) = plus1
               Work(N+j) = star1
               IF ( j+1<=N ) THEN
                  Work(j+1) = plus2
                  Work(N+j+1) = ZERO
                  plus1 = star1/plus2
!
!                 Generate a new *-value with norm between sqrt(TNORM)
!                 and TNORM.
!
                  rexp = SLARND(2,Iseed)
                  IF ( rexp<ZERO ) THEN
                     star1 = -sfac**(ONE-rexp)
                  ELSE
                     star1 = sfac**(ONE+rexp)
                  ENDIF
               ENDIF
            ENDDO
!
!           Copy the tridiagonal T to AB.
!
            IF ( upper ) THEN
               CALL SCOPY(N-1,Work,1,Ab(Kd,2),Ldab)
               CALL SCOPY(N-2,Work(N+1),1,Ab(Kd-1,3),Ldab)
            ELSE
               CALL SCOPY(N-1,Work,1,Ab(2,1),Ldab)
               CALL SCOPY(N-2,Work(N+1),1,Ab(3,1),Ldab)
            ENDIF
         ENDIF
!
!     IMAT > 9:  Pathological test cases.  These triangular matrices
!     are badly scaled or badly conditioned, so when used in solving a
!     triangular system they may cause overflow in the solution vector.
!
      ELSEIF ( Imat==10 ) THEN
!
!        Type 10:  Generate a triangular matrix with elements between
!        -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
!        Make the right hand side large so that it requires scaling.
!
         IF ( upper ) THEN
            DO j = 1 , N
               lenj = MIN(j,Kd+1)
               CALL SLARNV(2,Iseed,lenj,Ab(Kd+2-lenj,j))
               Ab(Kd+1,j) = SIGN(TWO,Ab(Kd+1,j))
            ENDDO
         ELSE
            DO j = 1 , N
               lenj = MIN(N-j+1,Kd+1)
               IF ( lenj>0 ) CALL SLARNV(2,Iseed,lenj,Ab(1,j))
               Ab(1,j) = SIGN(TWO,Ab(1,j))
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
      ELSEIF ( Imat==11 ) THEN
!
!        Type 11:  Make the first diagonal element in the solve small to
!        cause immediate overflow when dividing by T(j,j).
!        In type 11, the offdiagonal elements are small (CNORM(j) < 1).
!
         CALL SLARNV(2,Iseed,N,B)
         tscal = ONE/REAL(Kd+1)
         IF ( upper ) THEN
            DO j = 1 , N
               lenj = MIN(j,Kd+1)
               CALL SLARNV(2,Iseed,lenj,Ab(Kd+2-lenj,j))
               CALL SSCAL(lenj-1,tscal,Ab(Kd+2-lenj,j),1)
               Ab(Kd+1,j) = SIGN(ONE,Ab(Kd+1,j))
            ENDDO
            Ab(Kd+1,N) = smlnum*Ab(Kd+1,N)
         ELSE
            DO j = 1 , N
               lenj = MIN(N-j+1,Kd+1)
               CALL SLARNV(2,Iseed,lenj,Ab(1,j))
               IF ( lenj>1 ) CALL SSCAL(lenj-1,tscal,Ab(2,j),1)
               Ab(1,j) = SIGN(ONE,Ab(1,j))
            ENDDO
            Ab(1,1) = smlnum*Ab(1,1)
         ENDIF
!
      ELSEIF ( Imat==12 ) THEN
!
!        Type 12:  Make the first diagonal element in the solve small to
!        cause immediate overflow when dividing by T(j,j).
!        In type 12, the offdiagonal elements are O(1) (CNORM(j) > 1).
!
         CALL SLARNV(2,Iseed,N,B)
         IF ( upper ) THEN
            DO j = 1 , N
               lenj = MIN(j,Kd+1)
               CALL SLARNV(2,Iseed,lenj,Ab(Kd+2-lenj,j))
               Ab(Kd+1,j) = SIGN(ONE,Ab(Kd+1,j))
            ENDDO
            Ab(Kd+1,N) = smlnum*Ab(Kd+1,N)
         ELSE
            DO j = 1 , N
               lenj = MIN(N-j+1,Kd+1)
               CALL SLARNV(2,Iseed,lenj,Ab(1,j))
               Ab(1,j) = SIGN(ONE,Ab(1,j))
            ENDDO
            Ab(1,1) = smlnum*Ab(1,1)
         ENDIF
!
      ELSEIF ( Imat==13 ) THEN
!
!        Type 13:  T is diagonal with small numbers on the diagonal to
!        make the growth factor underflow, but a small right hand side
!        chosen so that the solution does not overflow.
!
         IF ( upper ) THEN
            jcount = 1
            DO j = N , 1 , -1
               DO i = MAX(1,Kd+1-(j-1)) , Kd
                  Ab(i,j) = ZERO
               ENDDO
               IF ( jcount<=2 ) THEN
                  Ab(Kd+1,j) = smlnum
               ELSE
                  Ab(Kd+1,j) = ONE
               ENDIF
               jcount = jcount + 1
               IF ( jcount>4 ) jcount = 1
            ENDDO
         ELSE
            jcount = 1
            DO j = 1 , N
               DO i = 2 , MIN(N-j+1,Kd+1)
                  Ab(i,j) = ZERO
               ENDDO
               IF ( jcount<=2 ) THEN
                  Ab(1,j) = smlnum
               ELSE
                  Ab(1,j) = ONE
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
      ELSEIF ( Imat==14 ) THEN
!
!        Type 14:  Make the diagonal elements small to cause gradual
!        overflow when dividing by T(j,j).  To control the amount of
!        scaling needed, the matrix is bidiagonal.
!
         texp = ONE/REAL(Kd+1)
         tscal = smlnum**texp
         CALL SLARNV(2,Iseed,N,B)
         IF ( upper ) THEN
            DO j = 1 , N
               DO i = MAX(1,Kd+2-j) , Kd
                  Ab(i,j) = ZERO
               ENDDO
               IF ( j>1 .AND. Kd>0 ) Ab(Kd,j) = -ONE
               Ab(Kd+1,j) = tscal
            ENDDO
            B(N) = ONE
         ELSE
            DO j = 1 , N
               DO i = 3 , MIN(N-j+1,Kd+1)
                  Ab(i,j) = ZERO
               ENDDO
               IF ( j<N .AND. Kd>0 ) Ab(2,j) = -ONE
               Ab(1,j) = tscal
            ENDDO
            B(1) = ONE
         ENDIF
!
      ELSEIF ( Imat==15 ) THEN
!
!        Type 15:  One zero diagonal element.
!
         iy = N/2 + 1
         IF ( upper ) THEN
            DO j = 1 , N
               lenj = MIN(j,Kd+1)
               CALL SLARNV(2,Iseed,lenj,Ab(Kd+2-lenj,j))
               IF ( j/=iy ) THEN
                  Ab(Kd+1,j) = SIGN(TWO,Ab(Kd+1,j))
               ELSE
                  Ab(Kd+1,j) = ZERO
               ENDIF
            ENDDO
         ELSE
            DO j = 1 , N
               lenj = MIN(N-j+1,Kd+1)
               CALL SLARNV(2,Iseed,lenj,Ab(1,j))
               IF ( j/=iy ) THEN
                  Ab(1,j) = SIGN(TWO,Ab(1,j))
               ELSE
                  Ab(1,j) = ZERO
               ENDIF
            ENDDO
         ENDIF
         CALL SLARNV(2,Iseed,N,B)
         CALL SSCAL(N,TWO,B,1)
!
      ELSEIF ( Imat==16 ) THEN
!
!        Type 16:  Make the offdiagonal elements large to cause overflow
!        when adding a column of T.  In the non-transposed case, the
!        matrix is constructed to cause overflow when adding a column in
!        every other step.
!
         tscal = unfl/ulp
         tscal = (ONE-ulp)/tscal
         DO j = 1 , N
            DO i = 1 , Kd + 1
               Ab(i,j) = ZERO
            ENDDO
         ENDDO
         texp = ONE
         IF ( Kd<=0 ) THEN
            DO j = 1 , N
               Ab(1,j) = ONE
               B(j) = REAL(j)
            ENDDO
         ELSEIF ( upper ) THEN
            DO j = N , 1 , -Kd
               DO i = j , MAX(1,j-Kd+1) , -2
                  Ab(1+(j-i),i) = -tscal/REAL(Kd+2)
                  Ab(Kd+1,i) = ONE
                  B(i) = texp*(ONE-ulp)
                  IF ( i>MAX(1,j-Kd+1) ) THEN
                     Ab(2+(j-i),i-1) = -(tscal/REAL(Kd+2))/REAL(Kd+3)
                     Ab(Kd+1,i-1) = ONE
                     B(i-1) = texp*REAL((Kd+1)*(Kd+1)+Kd)
                  ENDIF
                  texp = texp*TWO
               ENDDO
               B(MAX(1,j-Kd+1)) = (REAL(Kd+2)/REAL(Kd+3))*tscal
            ENDDO
         ELSE
            DO j = 1 , N , Kd
               texp = ONE
               lenj = MIN(Kd+1,N-j+1)
               DO i = j , MIN(N,j+Kd-1) , 2
                  Ab(lenj-(i-j),j) = -tscal/REAL(Kd+2)
                  Ab(1,j) = ONE
                  B(j) = texp*(ONE-ulp)
                  IF ( i<MIN(N,j+Kd-1) ) THEN
                     Ab(lenj-(i-j+1),i+1) = -(tscal/REAL(Kd+2))         &
     &                  /REAL(Kd+3)
                     Ab(1,i+1) = ONE
                     B(i+1) = texp*REAL((Kd+1)*(Kd+1)+Kd)
                  ENDIF
                  texp = texp*TWO
               ENDDO
               B(MIN(N,j+Kd-1)) = (REAL(Kd+2)/REAL(Kd+3))*tscal
            ENDDO
         ENDIF
!
      ELSEIF ( Imat==17 ) THEN
!
!        Type 17:  Generate a unit triangular matrix with elements
!        between -1 and 1, and make the right hand side large so that it
!        requires scaling.
!
         IF ( upper ) THEN
            DO j = 1 , N
               lenj = MIN(j-1,Kd)
               CALL SLARNV(2,Iseed,lenj,Ab(Kd+1-lenj,j))
               Ab(Kd+1,j) = REAL(j)
            ENDDO
         ELSE
            DO j = 1 , N
               lenj = MIN(N-j,Kd)
               IF ( lenj>0 ) CALL SLARNV(2,Iseed,lenj,Ab(2,j))
               Ab(1,j) = REAL(j)
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
      ELSEIF ( Imat==18 ) THEN
!
!        Type 18:  Generate a triangular matrix with elements between
!        BIGNUM/KD and BIGNUM so that at least one of the column
!        norms will exceed BIGNUM.
!
         tleft = bignum/MAX(ONE,REAL(Kd))
         tscal = bignum*(REAL(Kd)/REAL(Kd+1))
         IF ( upper ) THEN
            DO j = 1 , N
               lenj = MIN(j,Kd+1)
               CALL SLARNV(2,Iseed,lenj,Ab(Kd+2-lenj,j))
               DO i = Kd + 2 - lenj , Kd + 1
                  Ab(i,j) = SIGN(tleft,Ab(i,j)) + tscal*Ab(i,j)
               ENDDO
            ENDDO
         ELSE
            DO j = 1 , N
               lenj = MIN(N-j+1,Kd+1)
               CALL SLARNV(2,Iseed,lenj,Ab(1,j))
               DO i = 1 , lenj
                  Ab(i,j) = SIGN(tleft,Ab(i,j)) + tscal*Ab(i,j)
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
               lenj = MIN(N-2*j+1,Kd+1)
               CALL SSWAP(lenj,Ab(Kd+1,j),Ldab-1,Ab(Kd+2-lenj,N-j+1),-1)
            ENDDO
         ELSE
            DO j = 1 , N/2
               lenj = MIN(N-2*j+1,Kd+1)
               CALL SSWAP(lenj,Ab(1,j),1,Ab(lenj,N-j+2-lenj),-Ldab+1)
            ENDDO
         ENDIF
      ENDIF
!
!
!     End of slattb
!
      END SUBROUTINE slattb
