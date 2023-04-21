!*==slarrc.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLARRC computes the number of eigenvalues of the symmetric tridiagonal matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARRC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARRC( JOBT, N, VL, VU, D, E, PIVMIN,
!                                   EIGCNT, LCNT, RCNT, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBT
!       INTEGER            EIGCNT, INFO, LCNT, N, RCNT
!       REAL               PIVMIN, VL, VU
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Find the number of eigenvalues of the symmetric tridiagonal matrix T
!> that are in the interval (VL,VU] if JOBT = 'T', and of L D L^T
!> if JOBT = 'L'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBT
!> \verbatim
!>          JOBT is CHARACTER*1
!>          = 'T':  Compute Sturm count for matrix T.
!>          = 'L':  Compute Sturm count for matrix L D L^T.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix. N > 0.
!> \endverbatim
!>
!> \param[in] VL
!> \verbatim
!>          VL is REAL
!>          The lower bound for the eigenvalues.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>          VU is REAL
!>          The upper bound for the eigenvalues.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          JOBT = 'T': The N diagonal elements of the tridiagonal matrix T.
!>          JOBT = 'L': The N diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N)
!>          JOBT = 'T': The N-1 offdiagonal elements of the matrix T.
!>          JOBT = 'L': The N-1 offdiagonal elements of the matrix L.
!> \endverbatim
!>
!> \param[in] PIVMIN
!> \verbatim
!>          PIVMIN is REAL
!>          The minimum pivot in the Sturm sequence for T.
!> \endverbatim
!>
!> \param[out] EIGCNT
!> \verbatim
!>          EIGCNT is INTEGER
!>          The number of eigenvalues of the symmetric tridiagonal matrix T
!>          that are in the interval (VL,VU]
!> \endverbatim
!>
!> \param[out] LCNT
!> \verbatim
!>          LCNT is INTEGER
!> \endverbatim
!>
!> \param[out] RCNT
!> \verbatim
!>          RCNT is INTEGER
!>          The left and right negcounts of the interval.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
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
!> \date June 2016
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
      SUBROUTINE SLARRC(Jobt,N,Vl,Vu,D,E,Pivmin,Eigcnt,Lcnt,Rcnt,Info)
      IMPLICIT NONE
!*--SLARRC140
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobt
      INTEGER Eigcnt , Info , Lcnt , N , Rcnt
      REAL Pivmin , Vl , Vu
!     ..
!     .. Array Arguments ..
      REAL D(*) , E(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      LOGICAL matt
      REAL lpivot , rpivot , sl , su , tmp , tmp2
 
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. Executable Statements ..
!
      Info = 0
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
      Lcnt = 0
      Rcnt = 0
      Eigcnt = 0
      matt = LSAME(Jobt,'T')
 
 
      IF ( matt ) THEN
!        Sturm sequence count on T
         lpivot = D(1) - Vl
         rpivot = D(1) - Vu
         IF ( lpivot<=ZERO ) Lcnt = Lcnt + 1
         IF ( rpivot<=ZERO ) Rcnt = Rcnt + 1
         DO i = 1 , N - 1
            tmp = E(i)**2
            lpivot = (D(i+1)-Vl) - tmp/lpivot
            rpivot = (D(i+1)-Vu) - tmp/rpivot
            IF ( lpivot<=ZERO ) Lcnt = Lcnt + 1
            IF ( rpivot<=ZERO ) Rcnt = Rcnt + 1
         ENDDO
      ELSE
!        Sturm sequence count on L D L^T
         sl = -Vl
         su = -Vu
         DO i = 1 , N - 1
            lpivot = D(i) + sl
            rpivot = D(i) + su
            IF ( lpivot<=ZERO ) Lcnt = Lcnt + 1
            IF ( rpivot<=ZERO ) Rcnt = Rcnt + 1
            tmp = E(i)*D(i)*E(i)
!
            tmp2 = tmp/lpivot
            IF ( tmp2==ZERO ) THEN
               sl = tmp - Vl
            ELSE
               sl = sl*tmp2 - Vl
            ENDIF
!
            tmp2 = tmp/rpivot
            IF ( tmp2==ZERO ) THEN
               su = tmp - Vu
            ELSE
               su = su*tmp2 - Vu
            ENDIF
         ENDDO
         lpivot = D(N) + sl
         rpivot = D(N) + su
         IF ( lpivot<=ZERO ) Lcnt = Lcnt + 1
         IF ( rpivot<=ZERO ) Rcnt = Rcnt + 1
      ENDIF
      Eigcnt = Rcnt - Lcnt
 
!
!     end of SLARRC
!
      END SUBROUTINE SLARRC
