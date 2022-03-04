!*==zlatbs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLATBS solves a triangular banded system of equations.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLATBS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatbs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatbs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatbs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X,
!                          SCALE, CNORM, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORMIN, TRANS, UPLO
!       INTEGER            INFO, KD, LDAB, N
!       DOUBLE PRECISION   SCALE
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   CNORM( * )
!       COMPLEX*16         AB( LDAB, * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLATBS solves one of the triangular systems
!>
!>    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
!>
!> with scaling to prevent overflow, where A is an upper or lower
!> triangular band matrix.  Here A**T denotes the transpose of A, x and b
!> are n-element vectors, and s is a scaling factor, usually less than
!> or equal to 1, chosen so that the components of x will be less than
!> the overflow threshold.  If the unscaled problem will not cause
!> overflow, the Level 2 BLAS routine ZTBSV is called.  If the matrix A
!> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
!> non-trivial solution to A*x = 0 is returned.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the operation applied to A.
!>          = 'N':  Solve A * x = s*b     (No transpose)
!>          = 'T':  Solve A**T * x = s*b  (Transpose)
!>          = 'C':  Solve A**H * x = s*b  (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] NORMIN
!> \verbatim
!>          NORMIN is CHARACTER*1
!>          Specifies whether CNORM has been set or not.
!>          = 'Y':  CNORM contains the column norms on entry
!>          = 'N':  CNORM is not set on entry.  On exit, the norms will
!>                  be computed and stored in CNORM.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of subdiagonals or superdiagonals in the
!>          triangular matrix A.  KD >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          The upper or lower triangular band matrix A, stored in the
!>          first KD+1 rows of the array. The j-th column of A is stored
!>          in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (N)
!>          On entry, the right hand side b of the triangular system.
!>          On exit, X is overwritten by the solution vector x.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          The scaling factor s for the triangular system
!>             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.
!>          If SCALE = 0, the matrix A is singular or badly scaled, and
!>          the vector x is an exact or approximate solution to A*x = 0.
!> \endverbatim
!>
!> \param[in,out] CNORM
!> \verbatim
!>          CNORM is DOUBLE PRECISION array, dimension (N)
!>
!>          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
!>          contains the norm of the off-diagonal part of the j-th column
!>          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
!>          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
!>          must be greater than or equal to the 1-norm.
!>
!>          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
!>          returns the 1-norm of the offdiagonal part of the j-th column
!>          of A.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -k, the k-th argument had an illegal value
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
!> \date November 2017
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  A rough bound on x is computed; if that is less than overflow, ZTBSV
!>  is called, otherwise, specific code is used which checks for possible
!>  overflow or divide-by-zero at every operation.
!>
!>  A columnwise scheme is used for solving A*x = b.  The basic algorithm
!>  if A is lower triangular is
!>
!>       x[1:n] := b[1:n]
!>       for j = 1, ..., n
!>            x(j) := x(j) / A(j,j)
!>            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
!>       end
!>
!>  Define bounds on the components of x after j iterations of the loop:
!>     M(j) = bound on x[1:j]
!>     G(j) = bound on x[j+1:n]
!>  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
!>
!>  Then for iteration j+1 we have
!>     M(j+1) <= G(j) / | A(j+1,j+1) |
!>     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
!>            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
!>
!>  where CNORM(j+1) is greater than or equal to the infinity-norm of
!>  column j+1 of A, not counting the diagonal.  Hence
!>
!>     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
!>                  1<=i<=j
!>  and
!>
!>     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
!>                                   1<=i< j
!>
!>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine ZTBSV if the
!>  reciprocal of the largest M(j), j=1,..,n, is larger than
!>  max(underflow, 1/overflow).
!>
!>  The bound on x(j) is also used to determine when a step in the
!>  columnwise method can be performed without fear of overflow.  If
!>  the computed bound is greater than a large constant, x is scaled to
!>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
!>  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
!>
!>  Similarly, a row-wise scheme is used to solve A**T *x = b  or
!>  A**H *x = b.  The basic algorithm for A upper triangular is
!>
!>       for j = 1, ..., n
!>            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
!>       end
!>
!>  We simultaneously compute two bounds
!>       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
!>       M(j) = bound on x(i), 1<=i<=j
!>
!>  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
!>  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
!>  Then the bound on x(j) is
!>
!>       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
!>
!>            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
!>                      1<=i<=j
!>
!>  and we can safely call ZTBSV if 1/M(n) and 1/G(n) are both greater
!>  than max(underflow, 1/overflow).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZLATBS(Uplo,Trans,Diag,Normin,N,Kd,Ab,Ldab,X,Scale,    &
     &                  Cnorm,Info)
      USE F77KINDS                        
      USE S_DLABAD
      USE S_DLAMCH
      USE S_DSCAL
      USE S_DZASUM
      USE S_IDAMAX
      USE S_IZAMAX
      USE S_LSAME
      USE S_XERBLA
      USE S_ZAXPY
      USE S_ZDOTC
      USE S_ZDOTU
      USE S_ZDSCAL
      USE S_ZLADIV
      USE S_ZTBSV
      IMPLICIT NONE
!*--ZLATBS262
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , HALF = 0.5D+0 ,     &
     &                              ONE = 1.0D+0 , TWO = 2.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: bignum , grow , rec , smlnum , tjj , tmax ,       &
     &                tscal , xbnd , xj , xmax
      REAL(R8KIND) :: CABS1 , CABS2
      COMPLEX(CX16KIND) :: csumj , tjjs , uscal , zdum
      INTEGER :: i , imax , j , jfirst , jinc , jlast , jlen , maind
      LOGICAL :: notran , nounit , upper
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Statement Functions ..
!     ..
!     .. Statement Function definitions ..
      CABS1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))
      CABS2(zdum) = ABS(DBLE(zdum)/2.D0) + ABS(DIMAG(zdum)/2.D0)
!     ..
!     .. Executable Statements ..
!
      Info = 0
      upper = LSAME(Uplo,'U')
      notran = LSAME(Trans,'N')
      nounit = LSAME(Diag,'N')
!
!     Test the input parameters.
!
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( .NOT.notran .AND. .NOT.LSAME(Trans,'T') .AND.            &
     &         .NOT.LSAME(Trans,'C') ) THEN
         Info = -2
      ELSEIF ( .NOT.nounit .AND. .NOT.LSAME(Diag,'U') ) THEN
         Info = -3
      ELSEIF ( .NOT.LSAME(Normin,'Y') .AND. .NOT.LSAME(Normin,'N') )    &
     &         THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( Kd<0 ) THEN
         Info = -6
      ELSEIF ( Ldab<Kd+1 ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZLATBS',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Determine machine dependent parameters to control overflow.
!
      smlnum = DLAMCH('Safe minimum')
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
      smlnum = smlnum/DLAMCH('Precision')
      bignum = ONE/smlnum
      Scale = ONE
!
      IF ( LSAME(Normin,'N') ) THEN
!
!        Compute the 1-norm of each column, not including the diagonal.
!
         IF ( upper ) THEN
!
!           A is upper triangular.
!
            DO j = 1 , N
               jlen = MIN(Kd,j-1)
               Cnorm(j) = DZASUM(jlen,Ab(Kd+1-jlen,j),1)
            ENDDO
         ELSE
!
!           A is lower triangular.
!
            DO j = 1 , N
               jlen = MIN(Kd,N-j)
               IF ( jlen>0 ) THEN
                  Cnorm(j) = DZASUM(jlen,Ab(2,j),1)
               ELSE
                  Cnorm(j) = ZERO
               ENDIF
            ENDDO
         ENDIF
      ENDIF
!
!     Scale the column norms by TSCAL if the maximum element in CNORM is
!     greater than BIGNUM/2.
!
      imax = IDAMAX(N,Cnorm,1)
      tmax = Cnorm(imax)
      IF ( tmax<=bignum*HALF ) THEN
         tscal = ONE
      ELSE
         tscal = HALF/(smlnum*tmax)
         CALL DSCAL(N,tscal,Cnorm,1)
      ENDIF
!
!     Compute a bound on the computed solution vector to see if the
!     Level 2 BLAS routine ZTBSV can be used.
!
      xmax = ZERO
      DO j = 1 , N
         xmax = MAX(xmax,CABS2(X(j)))
      ENDDO
      xbnd = xmax
      IF ( notran ) THEN
!
!        Compute the growth in A * x = b.
!
         IF ( upper ) THEN
            jfirst = N
            jlast = 1
            jinc = -1
            maind = Kd + 1
         ELSE
            jfirst = 1
            jlast = N
            jinc = 1
            maind = 1
         ENDIF
!
         IF ( tscal/=ONE ) THEN
            grow = ZERO
            GOTO 100
         ENDIF
!
         IF ( nounit ) THEN
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, G(0) = max{x(i), i=1,...,n}.
!
            grow = HALF/MAX(xbnd,smlnum)
            xbnd = grow
            DO j = jfirst , jlast , jinc
!
!              Exit the loop if the growth factor is too small.
!
               IF ( grow<=smlnum ) GOTO 100
!
               tjjs = Ab(maind,j)
               tjj = CABS1(tjjs)
!
               IF ( tjj>=smlnum ) THEN
!
!                 M(j) = G(j-1) / abs(A(j,j))
!
                  xbnd = MIN(xbnd,MIN(ONE,tjj)*grow)
               ELSE
!
!                 M(j) could overflow, set XBND to 0.
!
                  xbnd = ZERO
               ENDIF
!
               IF ( tjj+Cnorm(j)>=smlnum ) THEN
!
!                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
!
                  grow = grow*(tjj/(tjj+Cnorm(j)))
               ELSE
!
!                 G(j) could overflow, set GROW to 0.
!
                  grow = ZERO
               ENDIF
            ENDDO
            grow = xbnd
         ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
            grow = MIN(ONE,HALF/MAX(xbnd,smlnum))
            DO j = jfirst , jlast , jinc
!
!              Exit the loop if the growth factor is too small.
!
               IF ( grow<=smlnum ) EXIT
!
!              G(j) = G(j-1)*( 1 + CNORM(j) )
!
               grow = grow*(ONE/(ONE+Cnorm(j)))
            ENDDO
         ENDIF
!
      ELSE
!
!        Compute the growth in A**T * x = b  or  A**H * x = b.
!
         IF ( upper ) THEN
            jfirst = 1
            jlast = N
            jinc = 1
            maind = Kd + 1
         ELSE
            jfirst = N
            jlast = 1
            jinc = -1
            maind = 1
         ENDIF
!
         IF ( tscal/=ONE ) THEN
            grow = ZERO
            GOTO 100
         ENDIF
!
         IF ( nounit ) THEN
!
!           A is non-unit triangular.
!
!           Compute GROW = 1/G(j) and XBND = 1/M(j).
!           Initially, M(0) = max{x(i), i=1,...,n}.
!
            grow = HALF/MAX(xbnd,smlnum)
            xbnd = grow
            DO j = jfirst , jlast , jinc
!
!              Exit the loop if the growth factor is too small.
!
               IF ( grow<=smlnum ) GOTO 100
!
!              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
!
               xj = ONE + Cnorm(j)
               grow = MIN(grow,xbnd/xj)
!
               tjjs = Ab(maind,j)
               tjj = CABS1(tjjs)
!
               IF ( tjj>=smlnum ) THEN
!
!                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
!
                  IF ( xj>tjj ) xbnd = xbnd*(tjj/xj)
               ELSE
!
!                 M(j) could overflow, set XBND to 0.
!
                  xbnd = ZERO
               ENDIF
            ENDDO
            grow = MIN(grow,xbnd)
         ELSE
!
!           A is unit triangular.
!
!           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
!
            grow = MIN(ONE,HALF/MAX(xbnd,smlnum))
            DO j = jfirst , jlast , jinc
!
!              Exit the loop if the growth factor is too small.
!
               IF ( grow<=smlnum ) EXIT
!
!              G(j) = ( 1 + CNORM(j) )*G(j-1)
!
               xj = ONE + Cnorm(j)
               grow = grow/xj
            ENDDO
         ENDIF
      ENDIF
!
 100  IF ( (grow*tscal)>smlnum ) THEN
!
!        Use the Level 2 BLAS solve if the reciprocal of the bound on
!        elements of X is not too small.
!
         CALL ZTBSV(Uplo,Trans,Diag,N,Kd,Ab,Ldab,X,1)
      ELSE
!
!        Use a Level 1 BLAS solve, scaling intermediate results.
!
         IF ( xmax>bignum*HALF ) THEN
!
!           Scale X so that its components are less than or equal to
!           BIGNUM in absolute value.
!
            Scale = (bignum*HALF)/xmax
            CALL ZDSCAL(N,Scale,X,1)
            xmax = bignum
         ELSE
            xmax = xmax*TWO
         ENDIF
!
         IF ( notran ) THEN
!
!           Solve A * x = b
!
            DO j = jfirst , jlast , jinc
!
!              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
!
               xj = CABS1(X(j))
               IF ( nounit ) THEN
                  tjjs = Ab(maind,j)*tscal
               ELSE
                  tjjs = tscal
                  IF ( tscal==ONE ) GOTO 110
               ENDIF
               tjj = CABS1(tjjs)
               IF ( tjj>smlnum ) THEN
!
!                    abs(A(j,j)) > SMLNUM:
!
                  IF ( tjj<ONE ) THEN
                     IF ( xj>tjj*bignum ) THEN
!
!                          Scale x by 1/b(j).
!
                        rec = ONE/xj
                        CALL ZDSCAL(N,rec,X,1)
                        Scale = Scale*rec
                        xmax = xmax*rec
                     ENDIF
                  ENDIF
                  X(j) = ZLADIV(X(j),tjjs)
                  xj = CABS1(X(j))
               ELSEIF ( tjj>ZERO ) THEN
!
!                    0 < abs(A(j,j)) <= SMLNUM:
!
                  IF ( xj>tjj*bignum ) THEN
!
!                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
!                       to avoid overflow when dividing by A(j,j).
!
                     rec = (tjj*bignum)/xj
!
!                          Scale by 1/CNORM(j) to avoid overflow when
!                          multiplying x(j) times column j.
!
                     IF ( Cnorm(j)>ONE ) rec = rec/Cnorm(j)
                     CALL ZDSCAL(N,rec,X,1)
                     Scale = Scale*rec
                     xmax = xmax*rec
                  ENDIF
                  X(j) = ZLADIV(X(j),tjjs)
                  xj = CABS1(X(j))
               ELSE
!
!                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                    scale = 0, and compute a solution to A*x = 0.
!
                  DO i = 1 , N
                     X(i) = ZERO
                  ENDDO
                  X(j) = ONE
                  xj = ONE
                  Scale = ZERO
                  xmax = ZERO
               ENDIF
!
!              Scale x if necessary to avoid overflow when adding a
!              multiple of column j of A.
!
 110           IF ( xj>ONE ) THEN
                  rec = ONE/xj
                  IF ( Cnorm(j)>(bignum-xmax)*rec ) THEN
!
!                    Scale x by 1/(2*abs(x(j))).
!
                     rec = rec*HALF
                     CALL ZDSCAL(N,rec,X,1)
                     Scale = Scale*rec
                  ENDIF
               ELSEIF ( xj*Cnorm(j)>(bignum-xmax) ) THEN
!
!                 Scale x by 1/2.
!
                  CALL ZDSCAL(N,HALF,X,1)
                  Scale = Scale*HALF
               ENDIF
!
               IF ( upper ) THEN
                  IF ( j>1 ) THEN
!
!                    Compute the update
!                       x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
!                                             x(j)* A(max(1,j-kd):j-1,j)
!
                     jlen = MIN(Kd,j-1)
                     CALL ZAXPY(jlen,-X(j)*tscal,Ab(Kd+1-jlen,j),1,     &
     &                          X(j-jlen),1)
                     i = IZAMAX(j-1,X,1)
                     xmax = CABS1(X(i))
                  ENDIF
               ELSEIF ( j<N ) THEN
!
!                 Compute the update
!                    x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
!                                          x(j) * A(j+1:min(j+kd,n),j)
!
                  jlen = MIN(Kd,N-j)
                  IF ( jlen>0 ) CALL ZAXPY(jlen,-X(j)*tscal,Ab(2,j),1,  &
     &                 X(j+1),1)
                  i = j + IZAMAX(N-j,X(j+1),1)
                  xmax = CABS1(X(i))
               ENDIF
            ENDDO
!
         ELSEIF ( LSAME(Trans,'T') ) THEN
!
!           Solve A**T * x = b
!
            DO j = jfirst , jlast , jinc
!
!              Compute x(j) = b(j) - sum A(k,j)*x(k).
!                                    k<>j
!
               xj = CABS1(X(j))
               uscal = tscal
               rec = ONE/MAX(xmax,ONE)
               IF ( Cnorm(j)>(bignum-xj)*rec ) THEN
!
!                 If x(j) could overflow, scale x by 1/(2*XMAX).
!
                  rec = rec*HALF
                  IF ( nounit ) THEN
                     tjjs = Ab(maind,j)*tscal
                  ELSE
                     tjjs = tscal
                  ENDIF
                  tjj = CABS1(tjjs)
                  IF ( tjj>ONE ) THEN
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     rec = MIN(ONE,rec*tjj)
                     uscal = ZLADIV(uscal,tjjs)
                  ENDIF
                  IF ( rec<ONE ) THEN
                     CALL ZDSCAL(N,rec,X,1)
                     Scale = Scale*rec
                     xmax = xmax*rec
                  ENDIF
               ENDIF
!
               csumj = ZERO
               IF ( uscal==DCMPLX(ONE) ) THEN
!
!                 If the scaling needed for A in the dot product is 1,
!                 call ZDOTU to perform the dot product.
!
                  IF ( upper ) THEN
                     jlen = MIN(Kd,j-1)
                     csumj = ZDOTU(jlen,Ab(Kd+1-jlen,j),1,X(j-jlen),1)
                  ELSE
                     jlen = MIN(Kd,N-j)
                     IF ( jlen>1 ) csumj = ZDOTU(jlen,Ab(2,j),1,X(j+1), &
     &                    1)
                  ENDIF
!
!                 Otherwise, use in-line code for the dot product.
!
               ELSEIF ( upper ) THEN
                  jlen = MIN(Kd,j-1)
                  DO i = 1 , jlen
                     csumj = csumj + (Ab(Kd+i-jlen,j)*uscal)            &
     &                       *X(j-jlen-1+i)
                  ENDDO
               ELSE
                  jlen = MIN(Kd,N-j)
                  DO i = 1 , jlen
                     csumj = csumj + (Ab(i+1,j)*uscal)*X(j+i)
                  ENDDO
               ENDIF
!
               IF ( uscal==DCMPLX(tscal) ) THEN
!
!                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
                  X(j) = X(j) - csumj
                  xj = CABS1(X(j))
                  IF ( nounit ) THEN
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
                     tjjs = Ab(maind,j)*tscal
                  ELSE
                     tjjs = tscal
                     IF ( tscal==ONE ) GOTO 120
                  ENDIF
                  tjj = CABS1(tjjs)
                  IF ( tjj>smlnum ) THEN
!
!                       abs(A(j,j)) > SMLNUM:
!
                     IF ( tjj<ONE ) THEN
                        IF ( xj>tjj*bignum ) THEN
!
!                             Scale X by 1/abs(x(j)).
!
                           rec = ONE/xj
                           CALL ZDSCAL(N,rec,X,1)
                           Scale = Scale*rec
                           xmax = xmax*rec
                        ENDIF
                     ENDIF
                     X(j) = ZLADIV(X(j),tjjs)
                  ELSEIF ( tjj>ZERO ) THEN
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     IF ( xj>tjj*bignum ) THEN
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        rec = (tjj*bignum)/xj
                        CALL ZDSCAL(N,rec,X,1)
                        Scale = Scale*rec
                        xmax = xmax*rec
                     ENDIF
                     X(j) = ZLADIV(X(j),tjjs)
                  ELSE
!
!                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                       scale = 0 and compute a solution to A**T *x = 0.
!
                     DO i = 1 , N
                        X(i) = ZERO
                     ENDDO
                     X(j) = ONE
                     Scale = ZERO
                     xmax = ZERO
                  ENDIF
               ELSE
!
!                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
!                 product has already been divided by 1/A(j,j).
!
                  X(j) = ZLADIV(X(j),tjjs) - csumj
               ENDIF
 120           xmax = MAX(xmax,CABS1(X(j)))
            ENDDO
!
         ELSE
!
!           Solve A**H * x = b
!
            DO j = jfirst , jlast , jinc
!
!              Compute x(j) = b(j) - sum A(k,j)*x(k).
!                                    k<>j
!
               xj = CABS1(X(j))
               uscal = tscal
               rec = ONE/MAX(xmax,ONE)
               IF ( Cnorm(j)>(bignum-xj)*rec ) THEN
!
!                 If x(j) could overflow, scale x by 1/(2*XMAX).
!
                  rec = rec*HALF
                  IF ( nounit ) THEN
                     tjjs = DCONJG(Ab(maind,j))*tscal
                  ELSE
                     tjjs = tscal
                  ENDIF
                  tjj = CABS1(tjjs)
                  IF ( tjj>ONE ) THEN
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     rec = MIN(ONE,rec*tjj)
                     uscal = ZLADIV(uscal,tjjs)
                  ENDIF
                  IF ( rec<ONE ) THEN
                     CALL ZDSCAL(N,rec,X,1)
                     Scale = Scale*rec
                     xmax = xmax*rec
                  ENDIF
               ENDIF
!
               csumj = ZERO
               IF ( uscal==DCMPLX(ONE) ) THEN
!
!                 If the scaling needed for A in the dot product is 1,
!                 call ZDOTC to perform the dot product.
!
                  IF ( upper ) THEN
                     jlen = MIN(Kd,j-1)
                     csumj = ZDOTC(jlen,Ab(Kd+1-jlen,j),1,X(j-jlen),1)
                  ELSE
                     jlen = MIN(Kd,N-j)
                     IF ( jlen>1 ) csumj = ZDOTC(jlen,Ab(2,j),1,X(j+1), &
     &                    1)
                  ENDIF
!
!                 Otherwise, use in-line code for the dot product.
!
               ELSEIF ( upper ) THEN
                  jlen = MIN(Kd,j-1)
                  DO i = 1 , jlen
                     csumj = csumj + (DCONJG(Ab(Kd+i-jlen,j))*uscal)    &
     &                       *X(j-jlen-1+i)
                  ENDDO
               ELSE
                  jlen = MIN(Kd,N-j)
                  DO i = 1 , jlen
                     csumj = csumj + (DCONJG(Ab(i+1,j))*uscal)*X(j+i)
                  ENDDO
               ENDIF
!
               IF ( uscal==DCMPLX(tscal) ) THEN
!
!                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
                  X(j) = X(j) - csumj
                  xj = CABS1(X(j))
                  IF ( nounit ) THEN
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
                     tjjs = DCONJG(Ab(maind,j))*tscal
                  ELSE
                     tjjs = tscal
                     IF ( tscal==ONE ) GOTO 130
                  ENDIF
                  tjj = CABS1(tjjs)
                  IF ( tjj>smlnum ) THEN
!
!                       abs(A(j,j)) > SMLNUM:
!
                     IF ( tjj<ONE ) THEN
                        IF ( xj>tjj*bignum ) THEN
!
!                             Scale X by 1/abs(x(j)).
!
                           rec = ONE/xj
                           CALL ZDSCAL(N,rec,X,1)
                           Scale = Scale*rec
                           xmax = xmax*rec
                        ENDIF
                     ENDIF
                     X(j) = ZLADIV(X(j),tjjs)
                  ELSEIF ( tjj>ZERO ) THEN
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     IF ( xj>tjj*bignum ) THEN
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        rec = (tjj*bignum)/xj
                        CALL ZDSCAL(N,rec,X,1)
                        Scale = Scale*rec
                        xmax = xmax*rec
                     ENDIF
                     X(j) = ZLADIV(X(j),tjjs)
                  ELSE
!
!                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
!                       scale = 0 and compute a solution to A**H *x = 0.
!
                     DO i = 1 , N
                        X(i) = ZERO
                     ENDDO
                     X(j) = ONE
                     Scale = ZERO
                     xmax = ZERO
                  ENDIF
               ELSE
!
!                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
!                 product has already been divided by 1/A(j,j).
!
                  X(j) = ZLADIV(X(j),tjjs) - csumj
               ENDIF
 130           xmax = MAX(xmax,CABS1(X(j)))
            ENDDO
         ENDIF
         Scale = Scale/tscal
      ENDIF
!
!     Scale the column norms by 1/TSCAL for return.
!
      IF ( tscal/=ONE ) CALL DSCAL(N,ONE/tscal,Cnorm,1)
!
!
!     End of ZLATBS
!
      END SUBROUTINE ZLATBS
