!*==slarrr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLARRR performs tests to decide whether the symmetric tridiagonal matrix T warrants expensive computations which guarantee high relative accuracy in the eigenvalues.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARRR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARRR( N, D, E, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            N, INFO
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * )
!       ..
!
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Perform tests to decide whether the symmetric tridiagonal matrix T
!> warrants expensive computations which guarantee high relative accuracy
!> in the eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix. N > 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The N diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N)
!>          On entry, the first (N-1) entries contain the subdiagonal
!>          elements of the tridiagonal matrix T; E(N) is set to ZERO.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          INFO = 0(default) : the matrix warrants computations preserving
!>                              relative accuracy.
!>          INFO = 1          : the matrix warrants computations guaranteeing
!>                              only absolute accuracy.
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
!> Beresford Parlett, University of California, Berkeley, USA \n
!> Jim Demmel, University of California, Berkeley, USA \n
!> Inderjit Dhillon, University of Texas, Austin, USA \n
!> Osni Marques, LBNL/NERSC, USA \n
!> Christof Voemel, University of California, Berkeley, USA
!
!  =====================================================================
      SUBROUTINE SLARRR(N,D,E,Info)
      USE S_SLAMCH
      IMPLICIT NONE
!*--SLARRR99
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , RELCOND = 0.999E0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: eps , offdig , offdig2 , rmin , safmin , smlnum , tmp ,   &
     &        tmp2
      INTEGER :: i
      LOGICAL :: yesrel
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
 
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N<=0 ) THEN
         Info = 0
         RETURN
      ENDIF
!
!     As a default, do NOT go for relative-accuracy preserving computations.
      Info = 1
 
      safmin = SLAMCH('Safe minimum')
      eps = SLAMCH('Precision')
      smlnum = safmin/eps
      rmin = SQRT(smlnum)
 
!     Tests for relative accuracy
!
!     Test for scaled diagonal dominance
!     Scale the diagonal entries to one and check whether the sum of the
!     off-diagonals is less than one
!
!     The sdd relative error bounds have a 1/(1- 2*x) factor in them,
!     x = max(OFFDIG + OFFDIG2), so when x is close to 1/2, no relative
!     accuracy is promised.  In the notation of the code fragment below,
!     1/(1 - (OFFDIG + OFFDIG2)) is the condition number.
!     We don't think it is worth going into "sdd mode" unless the relative
!     condition number is reasonable, not 1/macheps.
!     The threshold should be compatible with other thresholds used in the
!     code. We set  OFFDIG + OFFDIG2 <= .999 =: RELCOND, it corresponds
!     to losing at most 3 decimal digits: 1 / (1 - (OFFDIG + OFFDIG2)) <= 1000
!     instead of the current OFFDIG + OFFDIG2 < 1
!
      yesrel = .TRUE.
      offdig = ZERO
      tmp = SQRT(ABS(D(1)))
      IF ( tmp<rmin ) yesrel = .FALSE.
      IF ( yesrel ) THEN
         DO i = 2 , N
            tmp2 = SQRT(ABS(D(i)))
            IF ( tmp2<rmin ) yesrel = .FALSE.
            IF ( .NOT.yesrel ) EXIT
            offdig2 = ABS(E(i-1))/(tmp*tmp2)
            IF ( offdig+offdig2>=RELCOND ) yesrel = .FALSE.
            IF ( .NOT.yesrel ) EXIT
            tmp = tmp2
            offdig = offdig2
         ENDDO
      ENDIF
 
      IF ( yesrel ) THEN
         Info = 0
         RETURN
      ENDIF
!
 
!
!     *** MORE TO BE IMPLEMENTED ***
!
 
!
!     Test if the lower bidiagonal matrix L from T = L D L^T
!     (zero shift facto) is well conditioned
!
 
!
!     Test if the upper bidiagonal matrix U from T = U D U^T
!     (zero shift facto) is well conditioned.
!     In this case, the matrix needs to be flipped and, at the end
!     of the eigenvector computation, the flip needs to be applied
!     to the computed eigenvectors (and the support)
!
 
!
!
!     END OF SLARRR
!
      END SUBROUTINE SLARRR