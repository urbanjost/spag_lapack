!*==slasd2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLASD2 merges the two sets of singular values together into a single sorted set. Used by sbdsdc.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASD2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASD2( NL, NR, SQRE, K, D, Z, ALPHA, BETA, U, LDU, VT,
!                          LDVT, DSIGMA, U2, LDU2, VT2, LDVT2, IDXP, IDX,
!                          IDXC, IDXQ, COLTYP, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDU, LDU2, LDVT, LDVT2, NL, NR, SQRE
!       REAL               ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       INTEGER            COLTYP( * ), IDX( * ), IDXC( * ), IDXP( * ),
!      $                   IDXQ( * )
!       REAL               D( * ), DSIGMA( * ), U( LDU, * ),
!      $                   U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ),
!      $                   Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLASD2 merges the two sets of singular values together into a single
!> sorted set.  Then it tries to deflate the size of the problem.
!> There are two ways in which deflation can occur:  when two or more
!> singular values are close together or if there is a tiny entry in the
!> Z vector.  For each such occurrence the order of the related secular
!> equation problem is reduced by one.
!>
!> SLASD2 is called from SLASD1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NL
!> \verbatim
!>          NL is INTEGER
!>         The row dimension of the upper block.  NL >= 1.
!> \endverbatim
!>
!> \param[in] NR
!> \verbatim
!>          NR is INTEGER
!>         The row dimension of the lower block.  NR >= 1.
!> \endverbatim
!>
!> \param[in] SQRE
!> \verbatim
!>          SQRE is INTEGER
!>         = 0: the lower block is an NR-by-NR square matrix.
!>         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
!>
!>         The bidiagonal matrix has N = NL + NR + 1 rows and
!>         M = N + SQRE >= N columns.
!> \endverbatim
!>
!> \param[out] K
!> \verbatim
!>          K is INTEGER
!>         Contains the dimension of the non-deflated matrix,
!>         This is the order of the related secular equation. 1 <= K <=N.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>         On entry D contains the singular values of the two submatrices
!>         to be combined.  On exit D contains the trailing (N-K) updated
!>         singular values (those which were deflated) sorted into
!>         increasing order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is REAL array, dimension (N)
!>         On exit Z contains the updating row vector in the secular
!>         equation.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is REAL
!>         Contains the diagonal element associated with the added row.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is REAL
!>         Contains the off-diagonal element associated with the added
!>         row.
!> \endverbatim
!>
!> \param[in,out] U
!> \verbatim
!>          U is REAL array, dimension (LDU,N)
!>         On entry U contains the left singular vectors of two
!>         submatrices in the two square blocks with corners at (1,1),
!>         (NL, NL), and (NL+2, NL+2), (N,N).
!>         On exit U contains the trailing (N-K) updated left singular
!>         vectors (those which were deflated) in its last N-K columns.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>         The leading dimension of the array U.  LDU >= N.
!> \endverbatim
!>
!> \param[in,out] VT
!> \verbatim
!>          VT is REAL array, dimension (LDVT,M)
!>         On entry VT**T contains the right singular vectors of two
!>         submatrices in the two square blocks with corners at (1,1),
!>         (NL+1, NL+1), and (NL+2, NL+2), (M,M).
!>         On exit VT**T contains the trailing (N-K) updated right singular
!>         vectors (those which were deflated) in its last N-K columns.
!>         In case SQRE =1, the last row of VT spans the right null
!>         space.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>         The leading dimension of the array VT.  LDVT >= M.
!> \endverbatim
!>
!> \param[out] DSIGMA
!> \verbatim
!>          DSIGMA is REAL array, dimension (N)
!>         Contains a copy of the diagonal elements (K-1 singular values
!>         and one zero) in the secular equation.
!> \endverbatim
!>
!> \param[out] U2
!> \verbatim
!>          U2 is REAL array, dimension (LDU2,N)
!>         Contains a copy of the first K-1 left singular vectors which
!>         will be used by SLASD3 in a matrix multiply (SGEMM) to solve
!>         for the new left singular vectors. U2 is arranged into four
!>         blocks. The first block contains a column with 1 at NL+1 and
!>         zero everywhere else; the second block contains non-zero
!>         entries only at and above NL; the third contains non-zero
!>         entries only below NL+1; and the fourth is dense.
!> \endverbatim
!>
!> \param[in] LDU2
!> \verbatim
!>          LDU2 is INTEGER
!>         The leading dimension of the array U2.  LDU2 >= N.
!> \endverbatim
!>
!> \param[out] VT2
!> \verbatim
!>          VT2 is REAL array, dimension (LDVT2,N)
!>         VT2**T contains a copy of the first K right singular vectors
!>         which will be used by SLASD3 in a matrix multiply (SGEMM) to
!>         solve for the new right singular vectors. VT2 is arranged into
!>         three blocks. The first block contains a row that corresponds
!>         to the special 0 diagonal element in SIGMA; the second block
!>         contains non-zeros only at and before NL +1; the third block
!>         contains non-zeros only at and after  NL +2.
!> \endverbatim
!>
!> \param[in] LDVT2
!> \verbatim
!>          LDVT2 is INTEGER
!>         The leading dimension of the array VT2.  LDVT2 >= M.
!> \endverbatim
!>
!> \param[out] IDXP
!> \verbatim
!>          IDXP is INTEGER array, dimension (N)
!>         This will contain the permutation used to place deflated
!>         values of D at the end of the array. On output IDXP(2:K)
!>         points to the nondeflated D-values and IDXP(K+1:N)
!>         points to the deflated singular values.
!> \endverbatim
!>
!> \param[out] IDX
!> \verbatim
!>          IDX is INTEGER array, dimension (N)
!>         This will contain the permutation used to sort the contents of
!>         D into ascending order.
!> \endverbatim
!>
!> \param[out] IDXC
!> \verbatim
!>          IDXC is INTEGER array, dimension (N)
!>         This will contain the permutation used to arrange the columns
!>         of the deflated U matrix into three groups:  the first group
!>         contains non-zero entries only at and above NL, the second
!>         contains non-zero entries only below NL+2, and the third is
!>         dense.
!> \endverbatim
!>
!> \param[in,out] IDXQ
!> \verbatim
!>          IDXQ is INTEGER array, dimension (N)
!>         This contains the permutation which separately sorts the two
!>         sub-problems in D into ascending order.  Note that entries in
!>         the first hlaf of this permutation must first be moved one
!>         position backward; and entries in the second half
!>         must first have NL+1 added to their values.
!> \endverbatim
!>
!> \param[out] COLTYP
!> \verbatim
!>          COLTYP is INTEGER array, dimension (N)
!>         As workspace, this will contain a label which will indicate
!>         which of the following types a column in the U2 matrix or a
!>         row in the VT2 matrix is:
!>         1 : non-zero in the upper half only
!>         2 : non-zero in the lower half only
!>         3 : dense
!>         4 : deflated
!>
!>         On exit, it is an array of dimension 4, with COLTYP(I) being
!>         the dimension of the I-th type columns.
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
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Huan Ren, Computer Science Division, University of
!>     California at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE SLASD2(Nl,Nr,Sqre,K,D,Z,Alpha,Beta,U,Ldu,Vt,Ldvt,      &
     &                  Dsigma,U2,Ldu2,Vt2,Ldvt2,Idxp,Idx,Idxc,Idxq,    &
     &                  Coltyp,Info)
      USE S_SCOPY
      USE S_SLACPY
      USE S_SLAMCH
      USE S_SLAMRG
      USE S_SLAPY2
      USE S_SLASET
      USE S_SROT
      USE S_XERBLA
      IMPLICIT NONE
!*--SLASD2281
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , EIGHT = 8.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: Nl
      INTEGER :: Nr
      INTEGER , INTENT(IN) :: Sqre
      INTEGER , INTENT(INOUT) :: K
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) :: Beta
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , INTENT(INOUT) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dsigma
      REAL , INTENT(INOUT) , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      REAL , DIMENSION(Ldvt2,*) :: Vt2
      INTEGER :: Ldvt2
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxp
      INTEGER , DIMENSION(*) :: Idx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxc
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxq
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Coltyp
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: c , eps , hlftol , s , tau , tol , z1
      INTEGER :: ct , i , idxi , idxj , idxjp , j , jp , jprev , k2 ,   &
     &           m , n , nlp1 , nlp2
      INTEGER , DIMENSION(4) :: ctot , psm
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
!     .. Local Arrays ..
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
      IF ( Nl<1 ) THEN
         Info = -1
      ELSEIF ( Nr<1 ) THEN
         Info = -2
      ELSEIF ( (Sqre/=1) .AND. (Sqre/=0) ) THEN
         Info = -3
      ENDIF
!
      n = Nl + Nr + 1
      m = n + Sqre
!
      IF ( Ldu<n ) THEN
         Info = -10
      ELSEIF ( Ldvt<m ) THEN
         Info = -12
      ELSEIF ( Ldu2<n ) THEN
         Info = -15
      ELSEIF ( Ldvt2<m ) THEN
         Info = -17
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLASD2',-Info)
         RETURN
      ENDIF
!
      nlp1 = Nl + 1
      nlp2 = Nl + 2
!
!     Generate the first part of the vector Z; and move the singular
!     values in the first part of D one position backward.
!
      z1 = Alpha*Vt(nlp1,nlp1)
      Z(1) = z1
      DO i = Nl , 1 , -1
         Z(i+1) = Alpha*Vt(i,nlp1)
         D(i+1) = D(i)
         Idxq(i+1) = Idxq(i) + 1
      ENDDO
!
!     Generate the second part of the vector Z.
!
      DO i = nlp2 , m
         Z(i) = Beta*Vt(i,nlp2)
      ENDDO
!
!     Initialize some reference arrays.
!
      DO i = 2 , nlp1
         Coltyp(i) = 1
      ENDDO
      DO i = nlp2 , n
         Coltyp(i) = 2
      ENDDO
!
!     Sort the singular values into increasing order
!
      DO i = nlp2 , n
         Idxq(i) = Idxq(i) + nlp1
      ENDDO
!
!     DSIGMA, IDXC, IDXC, and the first column of U2
!     are used as storage space.
!
      DO i = 2 , n
         Dsigma(i) = D(Idxq(i))
         U2(i,1) = Z(Idxq(i))
         Idxc(i) = Coltyp(Idxq(i))
      ENDDO
!
      CALL SLAMRG(Nl,Nr,Dsigma(2),1,1,Idx(2))
!
      DO i = 2 , n
         idxi = 1 + Idx(i)
         D(i) = Dsigma(idxi)
         Z(i) = U2(idxi,1)
         Coltyp(i) = Idxc(idxi)
      ENDDO
!
!     Calculate the allowable deflation tolerance
!
      eps = SLAMCH('Epsilon')
      tol = MAX(ABS(Alpha),ABS(Beta))
      tol = EIGHT*eps*MAX(ABS(D(n)),tol)
!
!     There are 2 kinds of deflation -- first a value in the z-vector
!     is small, second two (or more) singular values are very close
!     together (their difference is small).
!
!     If the value in the z-vector is small, we simply permute the
!     array so that the corresponding singular value is moved to the
!     end.
!
!     If two values in the D-vector are close, we perform a two-sided
!     rotation designed to make one of the corresponding z-vector
!     entries zero, and then permute the array so that the deflated
!     singular value is moved to the end.
!
!     If there are multiple singular values then the problem deflates.
!     Here the number of equal singular values are found.  As each equal
!     singular value is found, an elementary reflector is computed to
!     rotate the corresponding singular subspace so that the
!     corresponding components of Z are zero in this new basis.
!
      K = 1
      k2 = n + 1
      DO j = 2 , n
         IF ( ABS(Z(j))<=tol ) THEN
!
!           Deflate due to small z component.
!
            k2 = k2 - 1
            Idxp(k2) = j
            Coltyp(j) = 4
            IF ( j==n ) GOTO 100
         ELSE
            jprev = j
            EXIT
         ENDIF
      ENDDO
      j = jprev
      DO
         j = j + 1
         IF ( j>n ) THEN
!
!     Record the last singular value.
!
            K = K + 1
            U2(K,1) = Z(jprev)
            Dsigma(K) = D(jprev)
            Idxp(K) = jprev
            EXIT
         ELSEIF ( ABS(Z(j))<=tol ) THEN
!
!        Deflate due to small z component.
!
            k2 = k2 - 1
            Idxp(k2) = j
            Coltyp(j) = 4
!
!        Check if singular values are close enough to allow deflation.
!
         ELSEIF ( ABS(D(j)-D(jprev))<=tol ) THEN
!
!           Deflation is possible.
!
            s = Z(jprev)
            c = Z(j)
!
!           Find sqrt(a**2+b**2) without overflow or
!           destructive underflow.
!
            tau = SLAPY2(c,s)
            c = c/tau
            s = -s/tau
            Z(j) = tau
            Z(jprev) = ZERO
!
!           Apply back the Givens rotation to the left and right
!           singular vector matrices.
!
            idxjp = Idxq(Idx(jprev)+1)
            idxj = Idxq(Idx(j)+1)
            IF ( idxjp<=nlp1 ) idxjp = idxjp - 1
            IF ( idxj<=nlp1 ) idxj = idxj - 1
            CALL SROT(n,U(1,idxjp),1,U(1,idxj),1,c,s)
            CALL SROT(m,Vt(idxjp,1),Ldvt,Vt(idxj,1),Ldvt,c,s)
            IF ( Coltyp(j)/=Coltyp(jprev) ) Coltyp(j) = 3
            Coltyp(jprev) = 4
            k2 = k2 - 1
            Idxp(k2) = jprev
            jprev = j
         ELSE
            K = K + 1
            U2(K,1) = Z(jprev)
            Dsigma(K) = D(jprev)
            Idxp(K) = jprev
            jprev = j
         ENDIF
      ENDDO
!
!
!     Count up the total number of the various types of columns, then
!     form a permutation which positions the four column types into
!     four groups of uniform structure (although one or more of these
!     groups may be empty).
!
 100  DO j = 1 , 4
         ctot(j) = 0
      ENDDO
      DO j = 2 , n
         ct = Coltyp(j)
         ctot(ct) = ctot(ct) + 1
      ENDDO
!
!     PSM(*) = Position in SubMatrix (of types 1 through 4)
!
      psm(1) = 2
      psm(2) = 2 + ctot(1)
      psm(3) = psm(2) + ctot(2)
      psm(4) = psm(3) + ctot(3)
!
!     Fill out the IDXC array so that the permutation which it induces
!     will place all type-1 columns first, all type-2 columns next,
!     then all type-3's, and finally all type-4's, starting from the
!     second column. This applies similarly to the rows of VT.
!
      DO j = 2 , n
         jp = Idxp(j)
         ct = Coltyp(jp)
         Idxc(psm(ct)) = j
         psm(ct) = psm(ct) + 1
      ENDDO
!
!     Sort the singular values and corresponding singular vectors into
!     DSIGMA, U2, and VT2 respectively.  The singular values/vectors
!     which were not deflated go into the first K slots of DSIGMA, U2,
!     and VT2 respectively, while those which were deflated go into the
!     last N - K slots, except that the first column/row will be treated
!     separately.
!
      DO j = 2 , n
         jp = Idxp(j)
         Dsigma(j) = D(jp)
         idxj = Idxq(Idx(Idxp(Idxc(j)))+1)
         IF ( idxj<=nlp1 ) idxj = idxj - 1
         CALL SCOPY(n,U(1,idxj),1,U2(1,j),1)
         CALL SCOPY(m,Vt(idxj,1),Ldvt,Vt2(j,1),Ldvt2)
      ENDDO
!
!     Determine DSIGMA(1), DSIGMA(2) and Z(1)
!
      Dsigma(1) = ZERO
      hlftol = tol/TWO
      IF ( ABS(Dsigma(2))<=hlftol ) Dsigma(2) = hlftol
      IF ( m>n ) THEN
         Z(1) = SLAPY2(z1,Z(m))
         IF ( Z(1)<=tol ) THEN
            c = ONE
            s = ZERO
            Z(1) = tol
         ELSE
            c = z1/Z(1)
            s = Z(m)/Z(1)
         ENDIF
      ELSEIF ( ABS(z1)<=tol ) THEN
         Z(1) = tol
      ELSE
         Z(1) = z1
      ENDIF
!
!     Move the rest of the updating row to Z.
!
      CALL SCOPY(K-1,U2(2,1),1,Z(2),1)
!
!     Determine the first column of U2, the first row of VT2 and the
!     last row of VT.
!
      CALL SLASET('A',n,1,ZERO,ZERO,U2,Ldu2)
      U2(nlp1,1) = ONE
      IF ( m>n ) THEN
         DO i = 1 , nlp1
            Vt(m,i) = -s*Vt(nlp1,i)
            Vt2(1,i) = c*Vt(nlp1,i)
         ENDDO
         DO i = nlp2 , m
            Vt2(1,i) = s*Vt(m,i)
            Vt(m,i) = c*Vt(m,i)
         ENDDO
      ELSE
         CALL SCOPY(m,Vt(nlp1,1),Ldvt,Vt2(1,1),Ldvt2)
      ENDIF
      IF ( m>n ) CALL SCOPY(m,Vt(m,1),Ldvt,Vt2(m,1),Ldvt2)
!
!     The deflated singular values and their corresponding vectors go
!     into the back of D, U, and V respectively.
!
      IF ( n>K ) THEN
         CALL SCOPY(n-K,Dsigma(K+1),1,D(K+1),1)
         CALL SLACPY('A',n,n-K,U2(1,K+1),Ldu2,U(1,K+1),Ldu)
         CALL SLACPY('A',n-K,m,Vt2(K+1,1),Ldvt2,Vt(K+1,1),Ldvt)
      ENDIF
!
!     Copy CTOT into COLTYP for referencing in SLASD3.
!
      DO j = 1 , 4
         Coltyp(j) = ctot(j)
      ENDDO
!
!
!     End of SLASD2
!
      END SUBROUTINE SLASD2
