!*==slasd7.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLASD7 merges the two sets of singular values together into a single sorted set. Then it tries to deflate the size of the problem. Used by sbdsdc.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASD7 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd7.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd7.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd7.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASD7( ICOMPQ, NL, NR, SQRE, K, D, Z, ZW, VF, VFW, VL,
!                          VLW, ALPHA, BETA, DSIGMA, IDX, IDXP, IDXQ,
!                          PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM,
!                          C, S, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL,
!      $                   NR, SQRE
!       REAL               ALPHA, BETA, C, S
!       ..
!       .. Array Arguments ..
!       INTEGER            GIVCOL( LDGCOL, * ), IDX( * ), IDXP( * ),
!      $                   IDXQ( * ), PERM( * )
!       REAL               D( * ), DSIGMA( * ), GIVNUM( LDGNUM, * ),
!      $                   VF( * ), VFW( * ), VL( * ), VLW( * ), Z( * ),
!      $                   ZW( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLASD7 merges the two sets of singular values together into a single
!> sorted set. Then it tries to deflate the size of the problem. There
!> are two ways in which deflation can occur:  when two or more singular
!> values are close together or if there is a tiny entry in the Z
!> vector. For each such occurrence the order of the related
!> secular equation problem is reduced by one.
!>
!> SLASD7 is called from SLASD6.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ICOMPQ
!> \verbatim
!>          ICOMPQ is INTEGER
!>          Specifies whether singular vectors are to be computed
!>          in compact form, as follows:
!>          = 0: Compute singular values only.
!>          = 1: Compute singular vectors of upper
!>               bidiagonal matrix in compact form.
!> \endverbatim
!>
!> \param[in] NL
!> \verbatim
!>          NL is INTEGER
!>         The row dimension of the upper block. NL >= 1.
!> \endverbatim
!>
!> \param[in] NR
!> \verbatim
!>          NR is INTEGER
!>         The row dimension of the lower block. NR >= 1.
!> \endverbatim
!>
!> \param[in] SQRE
!> \verbatim
!>          SQRE is INTEGER
!>         = 0: the lower block is an NR-by-NR square matrix.
!>         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
!>
!>         The bidiagonal matrix has
!>         N = NL + NR + 1 rows and
!>         M = N + SQRE >= N columns.
!> \endverbatim
!>
!> \param[out] K
!> \verbatim
!>          K is INTEGER
!>         Contains the dimension of the non-deflated matrix, this is
!>         the order of the related secular equation. 1 <= K <=N.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension ( N )
!>         On entry D contains the singular values of the two submatrices
!>         to be combined. On exit D contains the trailing (N-K) updated
!>         singular values (those which were deflated) sorted into
!>         increasing order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is REAL array, dimension ( M )
!>         On exit Z contains the updating row vector in the secular
!>         equation.
!> \endverbatim
!>
!> \param[out] ZW
!> \verbatim
!>          ZW is REAL array, dimension ( M )
!>         Workspace for Z.
!> \endverbatim
!>
!> \param[in,out] VF
!> \verbatim
!>          VF is REAL array, dimension ( M )
!>         On entry, VF(1:NL+1) contains the first components of all
!>         right singular vectors of the upper block; and VF(NL+2:M)
!>         contains the first components of all right singular vectors
!>         of the lower block. On exit, VF contains the first components
!>         of all right singular vectors of the bidiagonal matrix.
!> \endverbatim
!>
!> \param[out] VFW
!> \verbatim
!>          VFW is REAL array, dimension ( M )
!>         Workspace for VF.
!> \endverbatim
!>
!> \param[in,out] VL
!> \verbatim
!>          VL is REAL array, dimension ( M )
!>         On entry, VL(1:NL+1) contains the  last components of all
!>         right singular vectors of the upper block; and VL(NL+2:M)
!>         contains the last components of all right singular vectors
!>         of the lower block. On exit, VL contains the last components
!>         of all right singular vectors of the bidiagonal matrix.
!> \endverbatim
!>
!> \param[out] VLW
!> \verbatim
!>          VLW is REAL array, dimension ( M )
!>         Workspace for VL.
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
!> \param[out] DSIGMA
!> \verbatim
!>          DSIGMA is REAL array, dimension ( N )
!>         Contains a copy of the diagonal elements (K-1 singular values
!>         and one zero) in the secular equation.
!> \endverbatim
!>
!> \param[out] IDX
!> \verbatim
!>          IDX is INTEGER array, dimension ( N )
!>         This will contain the permutation used to sort the contents of
!>         D into ascending order.
!> \endverbatim
!>
!> \param[out] IDXP
!> \verbatim
!>          IDXP is INTEGER array, dimension ( N )
!>         This will contain the permutation used to place deflated
!>         values of D at the end of the array. On output IDXP(2:K)
!>         points to the nondeflated D-values and IDXP(K+1:N)
!>         points to the deflated singular values.
!> \endverbatim
!>
!> \param[in] IDXQ
!> \verbatim
!>          IDXQ is INTEGER array, dimension ( N )
!>         This contains the permutation which separately sorts the two
!>         sub-problems in D into ascending order.  Note that entries in
!>         the first half of this permutation must first be moved one
!>         position backward; and entries in the second half
!>         must first have NL+1 added to their values.
!> \endverbatim
!>
!> \param[out] PERM
!> \verbatim
!>          PERM is INTEGER array, dimension ( N )
!>         The permutations (from deflation and sorting) to be applied
!>         to each singular block. Not referenced if ICOMPQ = 0.
!> \endverbatim
!>
!> \param[out] GIVPTR
!> \verbatim
!>          GIVPTR is INTEGER
!>         The number of Givens rotations which took place in this
!>         subproblem. Not referenced if ICOMPQ = 0.
!> \endverbatim
!>
!> \param[out] GIVCOL
!> \verbatim
!>          GIVCOL is INTEGER array, dimension ( LDGCOL, 2 )
!>         Each pair of numbers indicates a pair of columns to take place
!>         in a Givens rotation. Not referenced if ICOMPQ = 0.
!> \endverbatim
!>
!> \param[in] LDGCOL
!> \verbatim
!>          LDGCOL is INTEGER
!>         The leading dimension of GIVCOL, must be at least N.
!> \endverbatim
!>
!> \param[out] GIVNUM
!> \verbatim
!>          GIVNUM is REAL array, dimension ( LDGNUM, 2 )
!>         Each number indicates the C or S value to be used in the
!>         corresponding Givens rotation. Not referenced if ICOMPQ = 0.
!> \endverbatim
!>
!> \param[in] LDGNUM
!> \verbatim
!>          LDGNUM is INTEGER
!>         The leading dimension of GIVNUM, must be at least N.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL
!>         C contains garbage if SQRE =0 and the C-value of a Givens
!>         rotation related to the right null space if SQRE = 1.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL
!>         S contains garbage if SQRE =0 and the S-value of a Givens
!>         rotation related to the right null space if SQRE = 1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>         = 0:  successful exit.
!>         < 0:  if INFO = -i, the i-th argument had an illegal value.
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
      SUBROUTINE SLASD7(Icompq,Nl,Nr,Sqre,K,D,Z,Zw,Vf,Vfw,Vl,Vlw,Alpha, &
     &                  Beta,Dsigma,Idx,Idxp,Idxq,Perm,Givptr,Givcol,   &
     &                  Ldgcol,Givnum,Ldgnum,C,S,Info)
      IMPLICIT NONE
!*--SLASD7283
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Givptr , Icompq , Info , K , Ldgcol , Ldgnum , Nl , Nr ,  &
     &        Sqre
      REAL Alpha , Beta , C , S
!     ..
!     .. Array Arguments ..
      INTEGER Givcol(Ldgcol,*) , Idx(*) , Idxp(*) , Idxq(*) , Perm(*)
      REAL D(*) , Dsigma(*) , Givnum(Ldgnum,*) , Vf(*) , Vfw(*) ,       &
     &     Vl(*) , Vlw(*) , Z(*) , Zw(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO , EIGHT
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0,TWO=2.0E+0,EIGHT=8.0E+0)
!     ..
!     .. Local Scalars ..
!
      INTEGER i , idxi , idxj , idxjp , j , jp , jprev , k2 , m , n ,   &
     &        nlp1 , nlp2
      REAL eps , hlftol , tau , tol , z1
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SLAMRG , SROT , XERBLA
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLAPY2
      EXTERNAL SLAMCH , SLAPY2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      n = Nl + Nr + 1
      m = n + Sqre
!
      IF ( (Icompq<0) .OR. (Icompq>1) ) THEN
         Info = -1
      ELSEIF ( Nl<1 ) THEN
         Info = -2
      ELSEIF ( Nr<1 ) THEN
         Info = -3
      ELSEIF ( (Sqre<0) .OR. (Sqre>1) ) THEN
         Info = -4
      ELSEIF ( Ldgcol<n ) THEN
         Info = -22
      ELSEIF ( Ldgnum<n ) THEN
         Info = -24
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLASD7',-Info)
         RETURN
      ENDIF
!
      nlp1 = Nl + 1
      nlp2 = Nl + 2
      IF ( Icompq==1 ) Givptr = 0
!
!     Generate the first part of the vector Z and move the singular
!     values in the first part of D one position backward.
!
      z1 = Alpha*Vl(nlp1)
      Vl(nlp1) = ZERO
      tau = Vf(nlp1)
      DO i = Nl , 1 , -1
         Z(i+1) = Alpha*Vl(i)
         Vl(i) = ZERO
         Vf(i+1) = Vf(i)
         D(i+1) = D(i)
         Idxq(i+1) = Idxq(i) + 1
      ENDDO
      Vf(1) = tau
!
!     Generate the second part of the vector Z.
!
      DO i = nlp2 , m
         Z(i) = Beta*Vf(i)
         Vf(i) = ZERO
      ENDDO
!
!     Sort the singular values into increasing order
!
      DO i = nlp2 , n
         Idxq(i) = Idxq(i) + nlp1
      ENDDO
!
!     DSIGMA, IDXC, IDXC, and ZW are used as storage space.
!
      DO i = 2 , n
         Dsigma(i) = D(Idxq(i))
         Zw(i) = Z(Idxq(i))
         Vfw(i) = Vf(Idxq(i))
         Vlw(i) = Vl(Idxq(i))
      ENDDO
!
      CALL SLAMRG(Nl,Nr,Dsigma(2),1,1,Idx(2))
!
      DO i = 2 , n
         idxi = 1 + Idx(i)
         D(i) = Dsigma(idxi)
         Z(i) = Zw(idxi)
         Vf(i) = Vfw(idxi)
         Vl(i) = Vlw(idxi)
      ENDDO
!
!     Calculate the allowable deflation tolerance
!
      eps = SLAMCH('Epsilon')
      tol = MAX(ABS(Alpha),ABS(Beta))
      tol = EIGHT*EIGHT*eps*MAX(ABS(D(n)),tol)
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
            Zw(K) = Z(jprev)
            Dsigma(K) = D(jprev)
            Idxp(K) = jprev
            EXIT
         ELSEIF ( ABS(Z(j))<=tol ) THEN
!
!        Deflate due to small z component.
!
            k2 = k2 - 1
            Idxp(k2) = j
!
!        Check if singular values are close enough to allow deflation.
!
         ELSEIF ( ABS(D(j)-D(jprev))<=tol ) THEN
!
!           Deflation is possible.
!
            S = Z(jprev)
            C = Z(j)
!
!           Find sqrt(a**2+b**2) without overflow or
!           destructive underflow.
!
            tau = SLAPY2(C,S)
            Z(j) = tau
            Z(jprev) = ZERO
            C = C/tau
            S = -S/tau
!
!           Record the appropriate Givens rotation
!
            IF ( Icompq==1 ) THEN
               Givptr = Givptr + 1
               idxjp = Idxq(Idx(jprev)+1)
               idxj = Idxq(Idx(j)+1)
               IF ( idxjp<=nlp1 ) idxjp = idxjp - 1
               IF ( idxj<=nlp1 ) idxj = idxj - 1
               Givcol(Givptr,2) = idxjp
               Givcol(Givptr,1) = idxj
               Givnum(Givptr,2) = C
               Givnum(Givptr,1) = S
            ENDIF
            CALL SROT(1,Vf(jprev),1,Vf(j),1,C,S)
            CALL SROT(1,Vl(jprev),1,Vl(j),1,C,S)
            k2 = k2 - 1
            Idxp(k2) = jprev
            jprev = j
         ELSE
            K = K + 1
            Zw(K) = Z(jprev)
            Dsigma(K) = D(jprev)
            Idxp(K) = jprev
            jprev = j
         ENDIF
      ENDDO
!
!
!     Sort the singular values into DSIGMA. The singular values which
!     were not deflated go into the first K slots of DSIGMA, except
!     that DSIGMA(1) is treated separately.
!
 100  DO j = 2 , n
         jp = Idxp(j)
         Dsigma(j) = D(jp)
         Vfw(j) = Vf(jp)
         Vlw(j) = Vl(jp)
      ENDDO
      IF ( Icompq==1 ) THEN
         DO j = 2 , n
            jp = Idxp(j)
            Perm(j) = Idxq(Idx(jp)+1)
            IF ( Perm(j)<=nlp1 ) Perm(j) = Perm(j) - 1
         ENDDO
      ENDIF
!
!     The deflated singular values go back into the last N - K slots of
!     D.
!
      CALL SCOPY(n-K,Dsigma(K+1),1,D(K+1),1)
!
!     Determine DSIGMA(1), DSIGMA(2), Z(1), VF(1), VL(1), VF(M), and
!     VL(M).
!
      Dsigma(1) = ZERO
      hlftol = tol/TWO
      IF ( ABS(Dsigma(2))<=hlftol ) Dsigma(2) = hlftol
      IF ( m>n ) THEN
         Z(1) = SLAPY2(z1,Z(m))
         IF ( Z(1)<=tol ) THEN
            C = ONE
            S = ZERO
            Z(1) = tol
         ELSE
            C = z1/Z(1)
            S = -Z(m)/Z(1)
         ENDIF
         CALL SROT(1,Vf(m),1,Vf(1),1,C,S)
         CALL SROT(1,Vl(m),1,Vl(1),1,C,S)
      ELSEIF ( ABS(z1)<=tol ) THEN
         Z(1) = tol
      ELSE
         Z(1) = z1
      ENDIF
!
!     Restore Z, VF, and VL.
!
      CALL SCOPY(K-1,Zw(2),1,Z(2),1)
      CALL SCOPY(n-1,Vfw(2),1,Vf(2),1)
      CALL SCOPY(n-1,Vlw(2),1,Vl(2),1)
!
!
!     End of SLASD7
!
      END SUBROUTINE SLASD7
