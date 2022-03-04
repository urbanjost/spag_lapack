!*==dlaneg.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLANEG computes the Sturm count.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLANEG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaneg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaneg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaneg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION DLANEG( N, D, LLD, SIGMA, PIVMIN, R )
!
!       .. Scalar Arguments ..
!       INTEGER            N, R
!       DOUBLE PRECISION   PIVMIN, SIGMA
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), LLD( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLANEG computes the Sturm count, the number of negative pivots
!> encountered while factoring tridiagonal T - sigma I = L D L^T.
!> This implementation works directly on the factors without forming
!> the tridiagonal matrix T.  The Sturm count is also the number of
!> eigenvalues of T less than sigma.
!>
!> This routine is called from DLARRB.
!>
!> The current routine does not use the PIVMIN parameter but rather
!> requires IEEE-754 propagation of Infinities and NaNs.  This
!> routine also has no input range restrictions but does require
!> default exception handling such that x/0 produces Inf when x is
!> non-zero, and Inf/Inf produces NaN.  For more information, see:
!>
!>   Marques, Riedy, and Voemel, "Benefits of IEEE-754 Features in
!>   Modern Symmetric Tridiagonal Eigensolvers," SIAM Journal on
!>   Scientific Computing, v28, n5, 2006.  DOI 10.1137/050641624
!>   (Tech report version in LAWN 172 with the same title.)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The N diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[in] LLD
!> \verbatim
!>          LLD is DOUBLE PRECISION array, dimension (N-1)
!>          The (N-1) elements L(i)*L(i)*D(i).
!> \endverbatim
!>
!> \param[in] SIGMA
!> \verbatim
!>          SIGMA is DOUBLE PRECISION
!>          Shift amount in T - sigma I = L D L^T.
!> \endverbatim
!>
!> \param[in] PIVMIN
!> \verbatim
!>          PIVMIN is DOUBLE PRECISION
!>          The minimum pivot in the Sturm sequence.  May be used
!>          when zero pivots are encountered on non-IEEE-754
!>          architectures.
!> \endverbatim
!>
!> \param[in] R
!> \verbatim
!>          R is INTEGER
!>          The twist index for the twisted factorization that is used
!>          for the negcount.
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
!>     Osni Marques, LBNL/NERSC, USA \n
!>     Christof Voemel, University of California, Berkeley, USA \n
!>     Jason Riedy, University of California, Berkeley, USA \n
!>
!  =====================================================================
      FUNCTION DLANEG(N,D,Lld,Sigma,Pivmin,R)
      USE F77KINDS                        
      USE S_DISNAN
      IMPLICIT NONE
!*--DLANEG124
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      INTEGER , PARAMETER  ::  BLKLEN = 128
      INTEGER :: DLANEG
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Lld
      REAL(R8KIND) , INTENT(IN) :: Sigma
      REAL(R8KIND) :: Pivmin
      INTEGER , INTENT(IN) :: R
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: bj , j , neg1 , neg2 , negcnt
      REAL(R8KIND) :: bsav , dminus , dplus , gamma , p , t , tmp
      LOGICAL :: sawnan
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
!     Some architectures propagate Infinities and NaNs very slowly, so
!     the code computes counts in BLKLEN chunks.  Then a NaN can
!     propagate at most BLKLEN columns before being detected.  This is
!     not a general tuning parameter; it needs only to be just large
!     enough that the overhead is tiny in common cases.
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
 
      negcnt = 0
 
!     I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
      t = -Sigma
      DO bj = 1 , R - 1 , BLKLEN
         neg1 = 0
         bsav = t
         DO j = bj , MIN(bj+BLKLEN-1,R-1)
            dplus = D(j) + t
            IF ( dplus<ZERO ) neg1 = neg1 + 1
            tmp = t/dplus
            t = tmp*Lld(j) - Sigma
         ENDDO
         sawnan = DISNAN(t)
!     Run a slower version of the above loop if a NaN is detected.
!     A NaN should occur only with a zero pivot after an infinite
!     pivot.  In that case, substituting 1 for T/DPLUS is the
!     correct limit.
         IF ( sawnan ) THEN
            neg1 = 0
            t = bsav
            DO j = bj , MIN(bj+BLKLEN-1,R-1)
               dplus = D(j) + t
               IF ( dplus<ZERO ) neg1 = neg1 + 1
               tmp = t/dplus
               IF ( DISNAN(tmp) ) tmp = ONE
               t = tmp*Lld(j) - Sigma
            ENDDO
         ENDIF
         negcnt = negcnt + neg1
      ENDDO
!
!     II) lower part: L D L^T - SIGMA I = U- D- U-^T
      p = D(N) - Sigma
      DO bj = N - 1 , R , -BLKLEN
         neg2 = 0
         bsav = p
         DO j = bj , MAX(bj-BLKLEN+1,R) , -1
            dminus = Lld(j) + p
            IF ( dminus<ZERO ) neg2 = neg2 + 1
            tmp = p/dminus
            p = tmp*D(j) - Sigma
         ENDDO
         sawnan = DISNAN(p)
!     As above, run a slower version that substitutes 1 for Inf/Inf.
!
         IF ( sawnan ) THEN
            neg2 = 0
            p = bsav
            DO j = bj , MAX(bj-BLKLEN+1,R) , -1
               dminus = Lld(j) + p
               IF ( dminus<ZERO ) neg2 = neg2 + 1
               tmp = p/dminus
               IF ( DISNAN(tmp) ) tmp = ONE
               p = tmp*D(j) - Sigma
            ENDDO
         ENDIF
         negcnt = negcnt + neg2
      ENDDO
!
!     III) Twist index
!       T was shifted by SIGMA initially.
      gamma = (t+Sigma) + p
      IF ( gamma<ZERO ) negcnt = negcnt + 1
 
      DLANEG = negcnt
      END FUNCTION DLANEG
