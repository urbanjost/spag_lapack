!*==dlaed9.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAED9 used by sstedc. Finds the roots of the secular equation and updates the eigenvectors. Used when the original matrix is dense.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAED9 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed9.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed9.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed9.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMDA, W,
!                          S, LDS, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, KSTART, KSTOP, LDQ, LDS, N
!       DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), S( LDS, * ),
!      $                   W( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED9 finds the roots of the secular equation, as defined by the
!> values in D, Z, and RHO, between KSTART and KSTOP.  It makes the
!> appropriate calls to DLAED4 and then stores the new matrix of
!> eigenvectors for use in calculating the next level of Z vectors.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of terms in the rational function to be solved by
!>          DLAED4.  K >= 0.
!> \endverbatim
!>
!> \param[in] KSTART
!> \verbatim
!>          KSTART is INTEGER
!> \endverbatim
!>
!> \param[in] KSTOP
!> \verbatim
!>          KSTOP is INTEGER
!>          The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP
!>          are to be computed.  1 <= KSTART <= KSTOP <= K.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns in the Q matrix.
!>          N >= K (delation may result in N > K).
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          D(I) contains the updated eigenvalues
!>          for KSTART <= I <= KSTOP.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max( 1, N ).
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>          The value of the parameter in the rank one update equation.
!>          RHO >= 0 required.
!> \endverbatim
!>
!> \param[in] DLAMDA
!> \verbatim
!>          DLAMDA is DOUBLE PRECISION array, dimension (K)
!>          The first K elements of this array contain the old roots
!>          of the deflated updating problem.  These are the poles
!>          of the secular equation.
!> \endverbatim
!>
!> \param[in] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (K)
!>          The first K elements of this array contain the components
!>          of the deflation-adjusted updating vector.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (LDS, K)
!>          Will contain the eigenvectors of the repaired matrix which
!>          will be stored for subsequent Z vector calculation and
!>          multiplied by the previously accumulated eigenvectors
!>          to update the system.
!> \endverbatim
!>
!> \param[in] LDS
!> \verbatim
!>          LDS is INTEGER
!>          The leading dimension of S.  LDS >= max( 1, K ).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, an eigenvalue did not converge
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
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
      SUBROUTINE DLAED9(K,Kstart,Kstop,N,D,Q,Ldq,Rho,Dlamda,W,S,Lds,    &
     &                  Info)
      IMPLICIT NONE
!*--DLAED9160
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , K , Kstart , Kstop , Ldq , Lds , N
      DOUBLE PRECISION Rho
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION D(*) , Dlamda(*) , Q(Ldq,*) , S(Lds,*) , W(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , j
      DOUBLE PRECISION temp
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMC3 , DNRM2
      EXTERNAL DLAMC3 , DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DLAED4 , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SIGN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( K<0 ) THEN
         Info = -1
      ELSEIF ( Kstart<1 .OR. Kstart>MAX(1,K) ) THEN
         Info = -2
      ELSEIF ( MAX(1,Kstop)<Kstart .OR. Kstop>MAX(1,K) ) THEN
         Info = -3
      ELSEIF ( N<K ) THEN
         Info = -4
      ELSEIF ( Ldq<MAX(1,K) ) THEN
         Info = -7
      ELSEIF ( Lds<MAX(1,K) ) THEN
         Info = -12
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DLAED9',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( K==0 ) RETURN
!
!     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DLAMDA(I) if it is 1; this makes the subsequent
!     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DLAMDA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DLAMDA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO i = 1 , N
         Dlamda(i) = DLAMC3(Dlamda(i),Dlamda(i)) - Dlamda(i)
      ENDDO
!
      DO j = Kstart , Kstop
         CALL DLAED4(K,j,Dlamda,W,Q(1,j),Rho,D(j),Info)
!
!        If the zero finder fails, the computation is terminated.
!
         IF ( Info/=0 ) GOTO 99999
      ENDDO
!
      IF ( K==1 .OR. K==2 ) THEN
         DO i = 1 , K
            DO j = 1 , K
               S(j,i) = Q(j,i)
            ENDDO
         ENDDO
         GOTO 99999
      ENDIF
!
!     Compute updated W.
!
      CALL DCOPY(K,W,1,S,1)
!
!     Initialize W(I) = Q(I,I)
!
      CALL DCOPY(K,Q,Ldq+1,W,1)
      DO j = 1 , K
         DO i = 1 , j - 1
            W(i) = W(i)*(Q(i,j)/(Dlamda(i)-Dlamda(j)))
         ENDDO
         DO i = j + 1 , K
            W(i) = W(i)*(Q(i,j)/(Dlamda(i)-Dlamda(j)))
         ENDDO
      ENDDO
      DO i = 1 , K
         W(i) = SIGN(SQRT(-W(i)),S(i,1))
      ENDDO
!
!     Compute eigenvectors of the modified rank-1 modification.
!
      DO j = 1 , K
         DO i = 1 , K
            Q(i,j) = W(i)/Q(i,j)
         ENDDO
         temp = DNRM2(K,Q(1,j),1)
         DO i = 1 , K
            S(i,j) = Q(i,j)/temp
         ENDDO
      ENDDO
!
!
!     End of DLAED9
!
99999 END SUBROUTINE DLAED9
