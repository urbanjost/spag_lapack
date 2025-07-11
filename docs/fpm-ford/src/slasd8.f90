!*==slasd8.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLASD8 finds the square roots of the roots of the secular equation, and stores, for each element in D, the distance to its two nearest poles. Used by sbdsdc.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASD8 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd8.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd8.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd8.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR,
!                          DSIGMA, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            ICOMPQ, INFO, K, LDDIFR
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), DIFL( * ), DIFR( LDDIFR, * ),
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
!> SLASD8 finds the square roots of the roots of the secular equation,
!> as defined by the values in DSIGMA and Z. It makes the appropriate
!> calls to SLASD4, and stores, for each  element in D, the distance
!> to its two nearest poles (elements in DSIGMA). It also updates
!> the arrays VF and VL, the first and last components of all the
!> right singular vectors of the original bidiagonal matrix.
!>
!> SLASD8 is called from SLASD6.
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
!>          by SLASD4.  K >= 1.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is REAL array, dimension ( K )
!>          On output, D contains the updated singular values.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is REAL array, dimension ( K )
!>          On entry, the first K elements of this array contain the
!>          components of the deflation-adjusted updating row vector.
!>          On exit, Z is updated.
!> \endverbatim
!>
!> \param[in,out] VF
!> \verbatim
!>          VF is REAL array, dimension ( K )
!>          On entry, VF contains  information passed through DBEDE8.
!>          On exit, VF contains the first K components of the first
!>          components of all right singular vectors of the bidiagonal
!>          matrix.
!> \endverbatim
!>
!> \param[in,out] VL
!> \verbatim
!>          VL is REAL array, dimension ( K )
!>          On entry, VL contains  information passed through DBEDE8.
!>          On exit, VL contains the first K components of the last
!>          components of all right singular vectors of the bidiagonal
!>          matrix.
!> \endverbatim
!>
!> \param[out] DIFL
!> \verbatim
!>          DIFL is REAL array, dimension ( K )
!>          On exit, DIFL(I) = D(I) - DSIGMA(I).
!> \endverbatim
!>
!> \param[out] DIFR
!> \verbatim
!>          DIFR is REAL array,
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
!>          DSIGMA is REAL array, dimension ( K )
!>          On entry, the first K elements of this array contain the old
!>          roots of the deflated updating problem.  These are the poles
!>          of the secular equation.
!>          On exit, the elements of DSIGMA may be very slightly altered
!>          in value.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (3*K)
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
      SUBROUTINE SLASD8(Icompq,K,D,Z,Vf,Vl,Difl,Difr,Lddifr,Dsigma,Work,&
     &                  Info)
      IMPLICIT NONE
!*--SLASD8170
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Icompq , Info , K , Lddifr
!     ..
!     .. Array Arguments ..
      REAL D(*) , Difl(*) , Difr(Lddifr,*) , Dsigma(*) , Vf(*) , Vl(*) ,&
     &     Work(*) , Z(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , iwk1 , iwk2 , iwk2i , iwk3 , iwk3i , j
      REAL diflj , difrj , dj , dsigj , dsigjp , rho , temp
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SLASCL , SLASD4 , SLASET , XERBLA
!     ..
!     .. External Functions ..
      REAL SDOT , SLAMC3 , SNRM2
      EXTERNAL SDOT , SLAMC3 , SNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SIGN , SQRT
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
         CALL XERBLA('SLASD8',-Info)
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
         Dsigma(i) = SLAMC3(Dsigma(i),Dsigma(i)) - Dsigma(i)
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
      rho = SNRM2(K,Z,1)
      CALL SLASCL('G',0,0,rho,ONE,K,1,Z,K,Info)
      rho = rho*rho
!
!     Initialize WORK(IWK3).
!
      CALL SLASET('A',K,1,ONE,ONE,Work(iwk3),K)
!
!     Compute the updated singular values, the arrays DIFL, DIFR,
!     and the updated Z.
!
      DO j = 1 , K
         CALL SLASD4(K,j,Dsigma,Z,Work(iwk1),rho,D(j),Work(iwk2),Info)
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
            Work(i) = Z(i)/(SLAMC3(Dsigma(i),dsigj)-diflj)              &
     &                /(Dsigma(i)+dj)
         ENDDO
         DO i = j + 1 , K
            Work(i) = Z(i)/(SLAMC3(Dsigma(i),dsigjp)+difrj)             &
     &                /(Dsigma(i)+dj)
         ENDDO
         temp = SNRM2(K,Work,1)
         Work(iwk2i+j) = SDOT(K,Work,1,Vf,1)/temp
         Work(iwk3i+j) = SDOT(K,Work,1,Vl,1)/temp
         IF ( Icompq==1 ) Difr(j,2) = temp
      ENDDO
!
      CALL SCOPY(K,Work(iwk2),1,Vf,1)
      CALL SCOPY(K,Work(iwk3),1,Vl,1)
!
!
!     End of SLASD8
!
      END SUBROUTINE SLASD8
