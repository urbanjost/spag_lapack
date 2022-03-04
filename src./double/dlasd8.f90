!*==dlasd8.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLASD8 finds the square roots of the roots of the secular equation, and stores, for each element in D, the distance to its two nearest poles. Used by sbdsdc.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASD8 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasd8.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasd8.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasd8.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR,
!                          DSIGMA, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            ICOMPQ, INFO, K, LDDIFR
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), DIFL( * ), DIFR( LDDIFR, * ),
!      $                   DSIGMA( * ), VF( * ), VL( * ), WORK( * ),
!      $                   Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASD8 finds the square roots of the roots of the secular equation,
!> as defined by the values in DSIGMA and Z. It makes the appropriate
!> calls to DLASD4, and stores, for each  element in D, the distance
!> to its two nearest poles (elements in DSIGMA). It also updates
!> the arrays VF and VL, the first and last components of all the
!> right singular vectors of the original bidiagonal matrix.
!>
!> DLASD8 is called from DLASD6.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ICOMPQ
!> \verbatim
!>          ICOMPQ is INTEGER
!>          Specifies whether singular vectors are to be computed in
!>          factored form in the calling routine:
!>          = 0: Compute singular values only.
!>          = 1: Compute singular vectors in factored form as well.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of terms in the rational function to be solved
!>          by DLASD4.  K >= 1.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension ( K )
!>          On output, D contains the updated singular values.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension ( K )
!>          On entry, the first K elements of this array contain the
!>          components of the deflation-adjusted updating row vector.
!>          On exit, Z is updated.
!> \endverbatim
!>
!> \param[in,out] VF
!> \verbatim
!>          VF is DOUBLE PRECISION array, dimension ( K )
!>          On entry, VF contains  information passed through DBEDE8.
!>          On exit, VF contains the first K components of the first
!>          components of all right singular vectors of the bidiagonal
!>          matrix.
!> \endverbatim
!>
!> \param[in,out] VL
!> \verbatim
!>          VL is DOUBLE PRECISION array, dimension ( K )
!>          On entry, VL contains  information passed through DBEDE8.
!>          On exit, VL contains the first K components of the last
!>          components of all right singular vectors of the bidiagonal
!>          matrix.
!> \endverbatim
!>
!> \param[out] DIFL
!> \verbatim
!>          DIFL is DOUBLE PRECISION array, dimension ( K )
!>          On exit, DIFL(I) = D(I) - DSIGMA(I).
!> \endverbatim
!>
!> \param[out] DIFR
!> \verbatim
!>          DIFR is DOUBLE PRECISION array,
!>                   dimension ( LDDIFR, 2 ) if ICOMPQ = 1 and
!>                   dimension ( K ) if ICOMPQ = 0.
!>          On exit, DIFR(I,1) = D(I) - DSIGMA(I+1), DIFR(K,1) is not
!>          defined and will not be referenced.
!>
!>          If ICOMPQ = 1, DIFR(1:K,2) is an array containing the
!>          normalizing factors for the right singular vector matrix.
!> \endverbatim
!>
!> \param[in] LDDIFR
!> \verbatim
!>          LDDIFR is INTEGER
!>          The leading dimension of DIFR, must be at least K.
!> \endverbatim
!>
!> \param[in,out] DSIGMA
!> \verbatim
!>          DSIGMA is DOUBLE PRECISION array, dimension ( K )
!>          On entry, the first K elements of this array contain the old
!>          roots of the deflated updating problem.  These are the poles
!>          of the secular equation.
!>          On exit, the elements of DSIGMA may be very slightly altered
!>          in value.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (3*K)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, a singular value did not converge
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
!> \date June 2017
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
      SUBROUTINE DLASD8(Icompq,K,D,Z,Vf,Vl,Difl,Difr,Lddifr,Dsigma,Work,&
     &                  Info)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DDOT
      USE S_DLAMC3
      USE S_DLASCL
      USE S_DLASD4
      USE S_DLASET
      USE S_DNRM2
      USE S_XERBLA
      IMPLICIT NONE
!*--DLASD8179
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Icompq
      INTEGER :: K
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL(R8KIND) , DIMENSION(*) :: Vf
      REAL(R8KIND) , DIMENSION(*) :: Vl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Difl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lddifr,*) :: Difr
      INTEGER , INTENT(IN) :: Lddifr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dsigma
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: diflj , difrj , dj , dsigj , dsigjp , rho , temp
      INTEGER :: i , iwk1 , iwk2 , iwk2i , iwk3 , iwk3i , j
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
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( (Icompq<0) .OR. (Icompq>1) ) THEN
         Info = -1
      ELSEIF ( K<1 ) THEN
         Info = -2
      ELSEIF ( Lddifr<K ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DLASD8',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( K==1 ) THEN
         D(1) = ABS(Z(1))
         Difl(1) = D(1)
         IF ( Icompq==1 ) THEN
            Difl(2) = ONE
            Difr(1,2) = ONE
         ENDIF
         RETURN
      ENDIF
!
!     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DSIGMA(I) if it is 1; this makes the subsequent
!     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DSIGMA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DSIGMA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO i = 1 , K
         Dsigma(i) = DLAMC3(Dsigma(i),Dsigma(i)) - Dsigma(i)
      ENDDO
!
!     Book keeping.
!
      iwk1 = 1
      iwk2 = iwk1 + K
      iwk3 = iwk2 + K
      iwk2i = iwk2 - 1
      iwk3i = iwk3 - 1
!
!     Normalize Z.
!
      rho = DNRM2(K,Z,1)
      CALL DLASCL('G',0,0,rho,ONE,K,1,Z,K,Info)
      rho = rho*rho
!
!     Initialize WORK(IWK3).
!
      CALL DLASET('A',K,1,ONE,ONE,Work(iwk3),K)
!
!     Compute the updated singular values, the arrays DIFL, DIFR,
!     and the updated Z.
!
      DO j = 1 , K
         CALL DLASD4(K,j,Dsigma,Z,Work(iwk1),rho,D(j),Work(iwk2),Info)
!
!        If the root finder fails, report the convergence failure.
!
         IF ( Info/=0 ) RETURN
         Work(iwk3i+j) = Work(iwk3i+j)*Work(j)*Work(iwk2i+j)
         Difl(j) = -Work(j)
         Difr(j,1) = -Work(j+1)
         DO i = 1 , j - 1
            Work(iwk3i+i) = Work(iwk3i+i)*Work(i)*Work(iwk2i+i)         &
     &                      /(Dsigma(i)-Dsigma(j))/(Dsigma(i)+Dsigma(j))
         ENDDO
         DO i = j + 1 , K
            Work(iwk3i+i) = Work(iwk3i+i)*Work(i)*Work(iwk2i+i)         &
     &                      /(Dsigma(i)-Dsigma(j))/(Dsigma(i)+Dsigma(j))
         ENDDO
      ENDDO
!
!     Compute updated Z.
!
      DO i = 1 , K
         Z(i) = SIGN(SQRT(ABS(Work(iwk3i+i))),Z(i))
      ENDDO
!
!     Update VF and VL.
!
      DO j = 1 , K
         diflj = Difl(j)
         dj = D(j)
         dsigj = -Dsigma(j)
         IF ( j<K ) THEN
            difrj = -Difr(j,1)
            dsigjp = -Dsigma(j+1)
         ENDIF
         Work(j) = -Z(j)/diflj/(Dsigma(j)+dj)
         DO i = 1 , j - 1
            Work(i) = Z(i)/(DLAMC3(Dsigma(i),dsigj)-diflj)              &
     &                /(Dsigma(i)+dj)
         ENDDO
         DO i = j + 1 , K
            Work(i) = Z(i)/(DLAMC3(Dsigma(i),dsigjp)+difrj)             &
     &                /(Dsigma(i)+dj)
         ENDDO
         temp = DNRM2(K,Work,1)
         Work(iwk2i+j) = DDOT(K,Work,1,Vf,1)/temp
         Work(iwk3i+j) = DDOT(K,Work,1,Vl,1)/temp
         IF ( Icompq==1 ) Difr(j,2) = temp
      ENDDO
!
      CALL DCOPY(K,Work(iwk2),1,Vf,1)
      CALL DCOPY(K,Work(iwk3),1,Vl,1)
!
!
!     End of DLASD8
!
      END SUBROUTINE DLASD8