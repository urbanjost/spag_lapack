!*==clatrs.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLATRS solves a triangular system of equations with the scale factor set to prevent overflow.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLATRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,
!                          CNORM, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORMIN, TRANS, UPLO
!       INTEGER            INFO, LDA, N
!       REAL               SCALE
!       ..
!       .. Array Arguments ..
!       REAL               CNORM( * )
!       COMPLEX            A( LDA, * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLATRS solves one of the triangular systems
!>
!>    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
!>
!> with scaling to prevent overflow.  Here A is an upper or lower
!> triangular matrix, A**T denotes the transpose of A, A**H denotes the
!> conjugate transpose of A, x and b are n-element vectors, and s is a
!> scaling factor, usually less than or equal to 1, chosen so that the
!> components of x will be less than the overflow threshold.  If the
!> unscaled problem will not cause overflow, the Level 2 BLAS routine
!> CTRSV is called. If the matrix A is singular (A(j,j) = 0 for some j),
!> then s is set to 0 and a non-trivial solution to A*x = 0 is returned.
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
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The triangular matrix A.  If UPLO = 'U', the leading n by n
!>          upper triangular part of the array A contains the upper
!>          triangular matrix, and the strictly lower triangular part of
!>          A is not referenced.  If UPLO = 'L', the leading n by n lower
!>          triangular part of the array A contains the lower triangular
!>          matrix, and the strictly upper triangular part of A is not
!>          referenced.  If DIAG = 'U', the diagonal elements of A are
!>          also not referenced and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max (1,N).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX array, dimension (N)
!>          On entry, the right hand side b of the triangular system.
!>          On exit, X is overwritten by the solution vector x.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is REAL
!>          The scaling factor s for the triangular system
!>             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.
!>          If SCALE = 0, the matrix A is singular or badly scaled, and
!>          the vector x is an exact or approximate solution to A*x = 0.
!> \endverbatim
!>
!> \param[in,out] CNORM
!> \verbatim
!>          CNORM is REAL array, dimension (N)
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
!> \date December 2016
!
!> \ingroup complexOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  A rough bound on x is computed; if that is less than overflow, CTRSV
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
!>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine CTRSV if the
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
!>  and we can safely call CTRSV if 1/M(n) and 1/G(n) are both greater
!>  than max(underflow, 1/overflow).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CLATRS(Uplo,Trans,Diag,Normin,N,A,Lda,X,Scale,Cnorm,   &
     &                  Info)
      IMPLICIT NONE
!*--CLATRS243
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Normin , Trans , Uplo
      INTEGER Info , Lda , N
      REAL Scale
!     ..
!     .. Array Arguments ..
      REAL Cnorm(*)
      COMPLEX A(Lda,*) , X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , HALF , ONE , TWO
      PARAMETER (ZERO=0.0E+0,HALF=0.5E+0,ONE=1.0E+0,TWO=2.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL notran , nounit , upper
      INTEGER i , imax , j , jfirst , jinc , jlast
      REAL bignum , grow , rec , smlnum , tjj , tmax , tscal , xbnd ,   &
     &     xj , xmax
      COMPLEX csumj , tjjs , uscal , zdum
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ICAMAX , ISAMAX
      REAL SCASUM , SLAMCH
      COMPLEX CDOTC , CDOTU , CLADIV
      EXTERNAL LSAME , ICAMAX , ISAMAX , SCASUM , SLAMCH , CDOTC ,      &
     &         CDOTU , CLADIV
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY , CSSCAL , CTRSV , SLABAD , SSCAL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , CMPLX , CONJG , MAX , MIN , REAL
!     ..
!     .. Statement Functions ..
      REAL CABS1 , CABS2
!     ..
!     .. Statement Function definitions ..
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
      CABS2(zdum) = ABS(REAL(zdum)/2.) + ABS(AIMAG(zdum)/2.)
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
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CLATRS',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Determine machine dependent parameters to control overflow.
!
      smlnum = SLAMCH('Safe minimum')
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
      smlnum = smlnum/SLAMCH('Precision')
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
               Cnorm(j) = SCASUM(j-1,A(1,j),1)
            ENDDO
         ELSE
!
!           A is lower triangular.
!
            DO j = 1 , N - 1
               Cnorm(j) = SCASUM(N-j,A(j+1,j),1)
            ENDDO
            Cnorm(N) = ZERO
         ENDIF
      ENDIF
!
!     Scale the column norms by TSCAL if the maximum element in CNORM is
!     greater than BIGNUM/2.
!
      imax = ISAMAX(N,Cnorm,1)
      tmax = Cnorm(imax)
      IF ( tmax<=bignum*HALF ) THEN
         tscal = ONE
      ELSE
         tscal = HALF/(smlnum*tmax)
         CALL SSCAL(N,tscal,Cnorm,1)
      ENDIF
!
!     Compute a bound on the computed solution vector to see if the
!     Level 2 BLAS routine CTRSV can be used.
!
      xmax = ZERO
      DO j = 1 , N
         xmax = MAX(xmax,CABS2(X(j)))
      ENDDO
      xbnd = xmax
!
      IF ( notran ) THEN
!
!        Compute the growth in A * x = b.
!
         IF ( upper ) THEN
            jfirst = N
            jlast = 1
            jinc = -1
         ELSE
            jfirst = 1
            jlast = N
            jinc = 1
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
               tjjs = A(j,j)
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
         ELSE
            jfirst = N
            jlast = 1
            jinc = -1
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
               tjjs = A(j,j)
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
         CALL CTRSV(Uplo,Trans,Diag,N,A,Lda,X,1)
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
            CALL CSSCAL(N,Scale,X,1)
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
                  tjjs = A(j,j)*tscal
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
                        CALL CSSCAL(N,rec,X,1)
                        Scale = Scale*rec
                        xmax = xmax*rec
                     ENDIF
                  ENDIF
                  X(j) = CLADIV(X(j),tjjs)
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
                     CALL CSSCAL(N,rec,X,1)
                     Scale = Scale*rec
                     xmax = xmax*rec
                  ENDIF
                  X(j) = CLADIV(X(j),tjjs)
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
                     CALL CSSCAL(N,rec,X,1)
                     Scale = Scale*rec
                  ENDIF
               ELSEIF ( xj*Cnorm(j)>(bignum-xmax) ) THEN
!
!                 Scale x by 1/2.
!
                  CALL CSSCAL(N,HALF,X,1)
                  Scale = Scale*HALF
               ENDIF
!
               IF ( upper ) THEN
                  IF ( j>1 ) THEN
!
!                    Compute the update
!                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
!
                     CALL CAXPY(j-1,-X(j)*tscal,A(1,j),1,X,1)
                     i = ICAMAX(j-1,X,1)
                     xmax = CABS1(X(i))
                  ENDIF
               ELSEIF ( j<N ) THEN
!
!                    Compute the update
!                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
!
                  CALL CAXPY(N-j,-X(j)*tscal,A(j+1,j),1,X(j+1),1)
                  i = j + ICAMAX(N-j,X(j+1),1)
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
                     tjjs = A(j,j)*tscal
                  ELSE
                     tjjs = tscal
                  ENDIF
                  tjj = CABS1(tjjs)
                  IF ( tjj>ONE ) THEN
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     rec = MIN(ONE,rec*tjj)
                     uscal = CLADIV(uscal,tjjs)
                  ENDIF
                  IF ( rec<ONE ) THEN
                     CALL CSSCAL(N,rec,X,1)
                     Scale = Scale*rec
                     xmax = xmax*rec
                  ENDIF
               ENDIF
!
               csumj = ZERO
               IF ( uscal==CMPLX(ONE) ) THEN
!
!                 If the scaling needed for A in the dot product is 1,
!                 call CDOTU to perform the dot product.
!
                  IF ( upper ) THEN
                     csumj = CDOTU(j-1,A(1,j),1,X,1)
                  ELSEIF ( j<N ) THEN
                     csumj = CDOTU(N-j,A(j+1,j),1,X(j+1),1)
                  ENDIF
!
!                 Otherwise, use in-line code for the dot product.
!
               ELSEIF ( upper ) THEN
                  DO i = 1 , j - 1
                     csumj = csumj + (A(i,j)*uscal)*X(i)
                  ENDDO
               ELSEIF ( j<N ) THEN
                  DO i = j + 1 , N
                     csumj = csumj + (A(i,j)*uscal)*X(i)
                  ENDDO
               ENDIF
!
               IF ( uscal==CMPLX(tscal) ) THEN
!
!                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
                  X(j) = X(j) - csumj
                  xj = CABS1(X(j))
                  IF ( nounit ) THEN
                     tjjs = A(j,j)*tscal
                  ELSE
                     tjjs = tscal
                     IF ( tscal==ONE ) GOTO 120
                  ENDIF
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
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
                           CALL CSSCAL(N,rec,X,1)
                           Scale = Scale*rec
                           xmax = xmax*rec
                        ENDIF
                     ENDIF
                     X(j) = CLADIV(X(j),tjjs)
                  ELSEIF ( tjj>ZERO ) THEN
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     IF ( xj>tjj*bignum ) THEN
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        rec = (tjj*bignum)/xj
                        CALL CSSCAL(N,rec,X,1)
                        Scale = Scale*rec
                        xmax = xmax*rec
                     ENDIF
                     X(j) = CLADIV(X(j),tjjs)
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
                  X(j) = CLADIV(X(j),tjjs) - csumj
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
                     tjjs = CONJG(A(j,j))*tscal
                  ELSE
                     tjjs = tscal
                  ENDIF
                  tjj = CABS1(tjjs)
                  IF ( tjj>ONE ) THEN
!
!                       Divide by A(j,j) when scaling x if A(j,j) > 1.
!
                     rec = MIN(ONE,rec*tjj)
                     uscal = CLADIV(uscal,tjjs)
                  ENDIF
                  IF ( rec<ONE ) THEN
                     CALL CSSCAL(N,rec,X,1)
                     Scale = Scale*rec
                     xmax = xmax*rec
                  ENDIF
               ENDIF
!
               csumj = ZERO
               IF ( uscal==CMPLX(ONE) ) THEN
!
!                 If the scaling needed for A in the dot product is 1,
!                 call CDOTC to perform the dot product.
!
                  IF ( upper ) THEN
                     csumj = CDOTC(j-1,A(1,j),1,X,1)
                  ELSEIF ( j<N ) THEN
                     csumj = CDOTC(N-j,A(j+1,j),1,X(j+1),1)
                  ENDIF
!
!                 Otherwise, use in-line code for the dot product.
!
               ELSEIF ( upper ) THEN
                  DO i = 1 , j - 1
                     csumj = csumj + (CONJG(A(i,j))*uscal)*X(i)
                  ENDDO
               ELSEIF ( j<N ) THEN
                  DO i = j + 1 , N
                     csumj = csumj + (CONJG(A(i,j))*uscal)*X(i)
                  ENDDO
               ENDIF
!
               IF ( uscal==CMPLX(tscal) ) THEN
!
!                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
!                 was not used to scale the dotproduct.
!
                  X(j) = X(j) - csumj
                  xj = CABS1(X(j))
                  IF ( nounit ) THEN
                     tjjs = CONJG(A(j,j))*tscal
                  ELSE
                     tjjs = tscal
                     IF ( tscal==ONE ) GOTO 130
                  ENDIF
!
!                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
!
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
                           CALL CSSCAL(N,rec,X,1)
                           Scale = Scale*rec
                           xmax = xmax*rec
                        ENDIF
                     ENDIF
                     X(j) = CLADIV(X(j),tjjs)
                  ELSEIF ( tjj>ZERO ) THEN
!
!                       0 < abs(A(j,j)) <= SMLNUM:
!
                     IF ( xj>tjj*bignum ) THEN
!
!                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
!
                        rec = (tjj*bignum)/xj
                        CALL CSSCAL(N,rec,X,1)
                        Scale = Scale*rec
                        xmax = xmax*rec
                     ENDIF
                     X(j) = CLADIV(X(j),tjjs)
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
                  X(j) = CLADIV(X(j),tjjs) - csumj
               ENDIF
 130           xmax = MAX(xmax,CABS1(X(j)))
            ENDDO
         ENDIF
         Scale = Scale/tscal
      ENDIF
!
!     Scale the column norms by 1/TSCAL for return.
!
      IF ( tscal/=ONE ) CALL SSCAL(N,ONE/tscal,Cnorm,1)
!
!
!     End of CLATRS
!
      END SUBROUTINE CLATRS
