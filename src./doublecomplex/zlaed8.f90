!*==zlaed8.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLAED8 used by sstedc. Merges eigenvalues and deflates secular equation. Used when the original matrix is dense.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAED8 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaed8.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaed8.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaed8.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, Z, DLAMDA,
!                          Q2, LDQ2, W, INDXP, INDX, INDXQ, PERM, GIVPTR,
!                          GIVCOL, GIVNUM, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            CUTPNT, GIVPTR, INFO, K, LDQ, LDQ2, N, QSIZ
!       DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
!       INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),
!      $                   INDXQ( * ), PERM( * )
!       DOUBLE PRECISION   D( * ), DLAMDA( * ), GIVNUM( 2, * ), W( * ),
!      $                   Z( * )
!       COMPLEX*16         Q( LDQ, * ), Q2( LDQ2, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAED8 merges the two sets of eigenvalues together into a single
!> sorted set.  Then it tries to deflate the size of the problem.
!> There are two ways in which deflation can occur:  when two or more
!> eigenvalues are close together or if there is a tiny element in the
!> Z vector.  For each such occurrence the order of the related secular
!> equation problem is reduced by one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] K
!> \verbatim
!>          K is INTEGER
!>         Contains the number of non-deflated eigenvalues.
!>         This is the order of the related secular equation.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is INTEGER
!>         The dimension of the unitary matrix used to reduce
!>         the dense or band matrix to tridiagonal form.
!>         QSIZ >= N if ICOMPQ = 1.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>         On entry, Q contains the eigenvectors of the partially solved
!>         system which has been previously updated in matrix
!>         multiplies with other partially solved eigensystems.
!>         On exit, Q contains the trailing (N-K) updated eigenvectors
!>         (those which were deflated) in its last N-K columns.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>         The leading dimension of the array Q.  LDQ >= max( 1, N ).
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>         On entry, D contains the eigenvalues of the two submatrices to
!>         be combined.  On exit, D contains the trailing (N-K) updated
!>         eigenvalues (those which were deflated) sorted into increasing
!>         order.
!> \endverbatim
!>
!> \param[in,out] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>         Contains the off diagonal element associated with the rank-1
!>         cut which originally split the two submatrices which are now
!>         being recombined. RHO is modified during the computation to
!>         the value required by DLAED3.
!> \endverbatim
!>
!> \param[in] CUTPNT
!> \verbatim
!>          CUTPNT is INTEGER
!>         Contains the location of the last eigenvalue in the leading
!>         sub-matrix.  MIN(1,N) <= CUTPNT <= N.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (N)
!>         On input this vector contains the updating vector (the last
!>         row of the first sub-eigenvector matrix and the first row of
!>         the second sub-eigenvector matrix).  The contents of Z are
!>         destroyed during the updating process.
!> \endverbatim
!>
!> \param[out] DLAMDA
!> \verbatim
!>          DLAMDA is DOUBLE PRECISION array, dimension (N)
!>         Contains a copy of the first K eigenvalues which will be used
!>         by DLAED3 to form the secular equation.
!> \endverbatim
!>
!> \param[out] Q2
!> \verbatim
!>          Q2 is COMPLEX*16 array, dimension (LDQ2,N)
!>         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,
!>         Contains a copy of the first K eigenvectors which will be used
!>         by DLAED7 in a matrix multiply (DGEMM) to update the new
!>         eigenvectors.
!> \endverbatim
!>
!> \param[in] LDQ2
!> \verbatim
!>          LDQ2 is INTEGER
!>         The leading dimension of the array Q2.  LDQ2 >= max( 1, N ).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>         This will hold the first k values of the final
!>         deflation-altered z-vector and will be passed to DLAED3.
!> \endverbatim
!>
!> \param[out] INDXP
!> \verbatim
!>          INDXP is INTEGER array, dimension (N)
!>         This will contain the permutation used to place deflated
!>         values of D at the end of the array. On output INDXP(1:K)
!>         points to the nondeflated D-values and INDXP(K+1:N)
!>         points to the deflated eigenvalues.
!> \endverbatim
!>
!> \param[out] INDX
!> \verbatim
!>          INDX is INTEGER array, dimension (N)
!>         This will contain the permutation used to sort the contents of
!>         D into ascending order.
!> \endverbatim
!>
!> \param[in] INDXQ
!> \verbatim
!>          INDXQ is INTEGER array, dimension (N)
!>         This contains the permutation which separately sorts the two
!>         sub-problems in D into ascending order.  Note that elements in
!>         the second half of this permutation must first have CUTPNT
!>         added to their values in order to be accurate.
!> \endverbatim
!>
!> \param[out] PERM
!> \verbatim
!>          PERM is INTEGER array, dimension (N)
!>         Contains the permutations (from deflation and sorting) to be
!>         applied to each eigenblock.
!> \endverbatim
!>
!> \param[out] GIVPTR
!> \verbatim
!>          GIVPTR is INTEGER
!>         Contains the number of Givens rotations which took place in
!>         this subproblem.
!> \endverbatim
!>
!> \param[out] GIVCOL
!> \verbatim
!>          GIVCOL is INTEGER array, dimension (2, N)
!>         Each pair of numbers indicates a pair of columns to take place
!>         in a Givens rotation.
!> \endverbatim
!>
!> \param[out] GIVNUM
!> \verbatim
!>          GIVNUM is DOUBLE PRECISION array, dimension (2, N)
!>         Each number indicates the S value to be used in the
!>         corresponding Givens rotation.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZLAED8(K,N,Qsiz,Q,Ldq,D,Rho,Cutpnt,Z,Dlamda,Q2,Ldq2,W, &
     &                  Indxp,Indx,Indxq,Perm,Givptr,Givcol,Givnum,Info)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DLAMCH
      USE S_DLAMRG
      USE S_DLAPY2
      USE S_DSCAL
      USE S_IDAMAX
      USE S_XERBLA
      USE S_ZCOPY
      USE S_ZDROT
      USE S_ZLACPY
      IMPLICIT NONE
!*--ZLAED8242
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  MONE = -1.0D0 , ZERO = 0.0D0 ,      &
     &                              ONE = 1.0D0 , TWO = 2.0D0 ,         &
     &                              EIGHT = 8.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(INOUT) :: K
      INTEGER :: N
      INTEGER :: Qsiz
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) :: Rho
      INTEGER , INTENT(IN) :: Cutpnt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dlamda
      COMPLEX(CX16KIND) , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxp
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxq
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Perm
      INTEGER , INTENT(INOUT) :: Givptr
      INTEGER , INTENT(OUT) , DIMENSION(2,*) :: Givcol
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(2,*) :: Givnum
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: c , eps , s , t , tau , tol
      INTEGER :: i , imax , j , jlam , jmax , jp , k2 , n1 , n1p1 , n2
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
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Qsiz<N ) THEN
         Info = -3
      ELSEIF ( Ldq<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Cutpnt<MIN(1,N) .OR. Cutpnt>N ) THEN
         Info = -8
      ELSEIF ( Ldq2<MAX(1,N) ) THEN
         Info = -12
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZLAED8',-Info)
         RETURN
      ENDIF
!
!     Need to initialize GIVPTR to O here in case of quick exit
!     to prevent an unspecified code behavior (usually sigfault)
!     when IWORK array on entry to *stedc is not zeroed
!     (or at least some IWORK entries which used in *laed7 for GIVPTR).
!
      Givptr = 0
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      n1 = Cutpnt
      n2 = N - n1
      n1p1 = n1 + 1
!
      IF ( Rho<ZERO ) CALL DSCAL(n2,MONE,Z(n1p1),1)
!
!     Normalize z so that norm(z) = 1
!
      t = ONE/SQRT(TWO)
      DO j = 1 , N
         Indx(j) = j
      ENDDO
      CALL DSCAL(N,t,Z,1)
      Rho = ABS(TWO*Rho)
!
!     Sort the eigenvalues into increasing order
!
      DO i = Cutpnt + 1 , N
         Indxq(i) = Indxq(i) + Cutpnt
      ENDDO
      DO i = 1 , N
         Dlamda(i) = D(Indxq(i))
         W(i) = Z(Indxq(i))
      ENDDO
      i = 1
      j = Cutpnt + 1
      CALL DLAMRG(n1,n2,Dlamda,1,1,Indx)
      DO i = 1 , N
         D(i) = Dlamda(Indx(i))
         Z(i) = W(Indx(i))
      ENDDO
!
!     Calculate the allowable deflation tolerance
!
      imax = IDAMAX(N,Z,1)
      jmax = IDAMAX(N,D,1)
      eps = DLAMCH('Epsilon')
      tol = EIGHT*eps*ABS(D(jmax))
!
!     If the rank-1 modifier is small enough, no more needs to be done
!     -- except to reorganize Q so that its columns correspond with the
!     elements in D.
!
      IF ( Rho*ABS(Z(imax))<=tol ) THEN
         K = 0
         DO j = 1 , N
            Perm(j) = Indxq(Indx(j))
            CALL ZCOPY(Qsiz,Q(1,Perm(j)),1,Q2(1,j),1)
         ENDDO
         CALL ZLACPY('A',Qsiz,N,Q2(1,1),Ldq2,Q(1,1),Ldq)
         RETURN
      ENDIF
!
!     If there are multiple eigenvalues then the problem deflates.  Here
!     the number of equal eigenvalues are found.  As each equal
!     eigenvalue is found, an elementary reflector is computed to rotate
!     the corresponding eigensubspace so that the corresponding
!     components of Z are zero in this new basis.
!
      K = 0
      k2 = N + 1
      DO j = 1 , N
         IF ( Rho*ABS(Z(j))<=tol ) THEN
!
!           Deflate due to small z component.
!
            k2 = k2 - 1
            Indxp(k2) = j
            IF ( j==N ) GOTO 100
         ELSE
            jlam = j
            EXIT
         ENDIF
      ENDDO
      DO
         j = j + 1
         IF ( j>N ) THEN
!
!     Record the last eigenvalue.
!
            K = K + 1
            W(K) = Z(jlam)
            Dlamda(K) = D(jlam)
            Indxp(K) = jlam
            EXIT
         ELSEIF ( Rho*ABS(Z(j))<=tol ) THEN
!
!        Deflate due to small z component.
!
            k2 = k2 - 1
            Indxp(k2) = j
         ELSE
!
!        Check if eigenvalues are close enough to allow deflation.
!
            s = Z(jlam)
            c = Z(j)
!
!        Find sqrt(a**2+b**2) without overflow or
!        destructive underflow.
!
            tau = DLAPY2(c,s)
            t = D(j) - D(jlam)
            c = c/tau
            s = -s/tau
            IF ( ABS(t*c*s)<=tol ) THEN
!
!           Deflation is possible.
!
               Z(j) = tau
               Z(jlam) = ZERO
!
!           Record the appropriate Givens rotation
!
               Givptr = Givptr + 1
               Givcol(1,Givptr) = Indxq(Indx(jlam))
               Givcol(2,Givptr) = Indxq(Indx(j))
               Givnum(1,Givptr) = c
               Givnum(2,Givptr) = s
               CALL ZDROT(Qsiz,Q(1,Indxq(Indx(jlam))),1,                &
     &                    Q(1,Indxq(Indx(j))),1,c,s)
               t = D(jlam)*c*c + D(j)*s*s
               D(j) = D(jlam)*s*s + D(j)*c*c
               D(jlam) = t
               k2 = k2 - 1
               i = 1
               DO
                  IF ( k2+i>N ) THEN
                     Indxp(k2+i-1) = jlam
                  ELSEIF ( D(jlam)<D(Indxp(k2+i)) ) THEN
                     Indxp(k2+i-1) = Indxp(k2+i)
                     Indxp(k2+i) = jlam
                     i = i + 1
                     CYCLE
                  ELSE
                     Indxp(k2+i-1) = jlam
                  ENDIF
                  jlam = j
                  EXIT
               ENDDO
            ELSE
               K = K + 1
               W(K) = Z(jlam)
               Dlamda(K) = D(jlam)
               Indxp(K) = jlam
               jlam = j
            ENDIF
         ENDIF
      ENDDO
!
!
!     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
!     and Q2 respectively.  The eigenvalues/vectors which were not
!     deflated go into the first K slots of DLAMDA and Q2 respectively,
!     while those which were deflated go into the last N - K slots.
!
 100  DO j = 1 , N
         jp = Indxp(j)
         Dlamda(j) = D(jp)
         Perm(j) = Indxq(Indx(jp))
         CALL ZCOPY(Qsiz,Q(1,Perm(j)),1,Q2(1,j),1)
      ENDDO
!
!     The deflated eigenvalues and their corresponding vectors go back
!     into the last N - K slots of D and Q respectively.
!
      IF ( K<N ) THEN
         CALL DCOPY(N-K,Dlamda(K+1),1,D(K+1),1)
         CALL ZLACPY('A',Qsiz,N-K,Q2(1,K+1),Ldq2,Q(1,K+1),Ldq)
      ENDIF
!
!
!     End of ZLAED8
!
      END SUBROUTINE ZLAED8
