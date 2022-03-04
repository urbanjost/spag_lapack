!*==dlarrf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARRF finds a new relatively robust representation such that at least one of the eigenvalues is relatively isolated.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARRF( N, D, L, LD, CLSTRT, CLEND,
!                          W, WGAP, WERR,
!                          SPDIAM, CLGAPL, CLGAPR, PIVMIN, SIGMA,
!                          DPLUS, LPLUS, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            CLSTRT, CLEND, INFO, N
!       DOUBLE PRECISION   CLGAPL, CLGAPR, PIVMIN, SIGMA, SPDIAM
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), DPLUS( * ), L( * ), LD( * ),
!      $          LPLUS( * ), W( * ), WGAP( * ), WERR( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Given the initial representation L D L^T and its cluster of close
!> eigenvalues (in a relative measure), W( CLSTRT ), W( CLSTRT+1 ), ...
!> W( CLEND ), DLARRF finds a new relatively robust representation
!> L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the
!> eigenvalues of L(+) D(+) L(+)^T is relatively isolated.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix (subblock, if the matrix split).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The N diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[in] L
!> \verbatim
!>          L is DOUBLE PRECISION array, dimension (N-1)
!>          The (N-1) subdiagonal elements of the unit bidiagonal
!>          matrix L.
!> \endverbatim
!>
!> \param[in] LD
!> \verbatim
!>          LD is DOUBLE PRECISION array, dimension (N-1)
!>          The (N-1) elements L(i)*D(i).
!> \endverbatim
!>
!> \param[in] CLSTRT
!> \verbatim
!>          CLSTRT is INTEGER
!>          The index of the first eigenvalue in the cluster.
!> \endverbatim
!>
!> \param[in] CLEND
!> \verbatim
!>          CLEND is INTEGER
!>          The index of the last eigenvalue in the cluster.
!> \endverbatim
!>
!> \param[in] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension
!>          dimension is >=  (CLEND-CLSTRT+1)
!>          The eigenvalue APPROXIMATIONS of L D L^T in ascending order.
!>          W( CLSTRT ) through W( CLEND ) form the cluster of relatively
!>          close eigenalues.
!> \endverbatim
!>
!> \param[in,out] WGAP
!> \verbatim
!>          WGAP is DOUBLE PRECISION array, dimension
!>          dimension is >=  (CLEND-CLSTRT+1)
!>          The separation from the right neighbor eigenvalue in W.
!> \endverbatim
!>
!> \param[in] WERR
!> \verbatim
!>          WERR is DOUBLE PRECISION array, dimension
!>          dimension is  >=  (CLEND-CLSTRT+1)
!>          WERR contain the semiwidth of the uncertainty
!>          interval of the corresponding eigenvalue APPROXIMATION in W
!> \endverbatim
!>
!> \param[in] SPDIAM
!> \verbatim
!>          SPDIAM is DOUBLE PRECISION
!>          estimate of the spectral diameter obtained from the
!>          Gerschgorin intervals
!> \endverbatim
!>
!> \param[in] CLGAPL
!> \verbatim
!>          CLGAPL is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] CLGAPR
!> \verbatim
!>          CLGAPR is DOUBLE PRECISION
!>          absolute gap on each end of the cluster.
!>          Set by the calling routine to protect against shifts too close
!>          to eigenvalues outside the cluster.
!> \endverbatim
!>
!> \param[in] PIVMIN
!> \verbatim
!>          PIVMIN is DOUBLE PRECISION
!>          The minimum pivot allowed in the Sturm sequence.
!> \endverbatim
!>
!> \param[out] SIGMA
!> \verbatim
!>          SIGMA is DOUBLE PRECISION
!>          The shift used to form L(+) D(+) L(+)^T.
!> \endverbatim
!>
!> \param[out] DPLUS
!> \verbatim
!>          DPLUS is DOUBLE PRECISION array, dimension (N)
!>          The N diagonal elements of the diagonal matrix D(+).
!> \endverbatim
!>
!> \param[out] LPLUS
!> \verbatim
!>          LPLUS is DOUBLE PRECISION array, dimension (N-1)
!>          The first (N-1) elements of LPLUS contain the subdiagonal
!>          elements of the unit bidiagonal matrix L(+).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (2*N)
!>          Workspace.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          Signals processing OK (=0) or failure (=1)
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
      SUBROUTINE DLARRF(N,D,L,Ld,Clstrt,Clend,W,Wgap,Werr,Spdiam,Clgapl,&
     &                  Clgapr,Pivmin,Sigma,Dplus,Lplus,Work,Info)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DISNAN
      USE S_DLAMCH
      IMPLICIT NONE
!*--DLARRF199
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , TWO = 2.0D0 ,         &
     &                              FOUR = 4.0D0 , QUART = 0.25D0 ,     &
     &                              MAXGROWTH1 = 8.D0 ,                 &
     &                              MAXGROWTH2 = 8.D0
      INTEGER , PARAMETER  ::  KTRYMAX = 1 , SLEFT = 1 , SRIGHT = 2
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: L
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Ld
      INTEGER , INTENT(IN) :: Clstrt
      INTEGER , INTENT(IN) :: Clend
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Wgap
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Werr
      REAL(R8KIND) , INTENT(IN) :: Spdiam
      REAL(R8KIND) , INTENT(IN) :: Clgapl
      REAL(R8KIND) , INTENT(IN) :: Clgapr
      REAL(R8KIND) , INTENT(IN) :: Pivmin
      REAL(R8KIND) , INTENT(OUT) :: Sigma
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dplus
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Lplus
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: avgap , bestshift , clwdth , eps , fact , fail ,  &
     &                fail2 , growthbound , ldelta , ldmax , lsigma ,   &
     &                max1 , max2 , mingap , oldp , prod , rdelta ,     &
     &                rdmax , rrr1 , rrr2 , rsigma , s , smlgrowth ,    &
     &                tmp , znm2
      LOGICAL :: dorrr1 , forcer , nofail , sawnan1 , sawnan2 , tryrrr1
      INTEGER :: i , indx , ktry , shift
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
      Info = 0
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
      fact = DBLE(2**KTRYMAX)
      eps = DLAMCH('Precision')
      shift = 0
      forcer = .FALSE.
 
 
!     Note that we cannot guarantee that for any of the shifts tried,
!     the factorization has a small or even moderate element growth.
!     There could be Ritz values at both ends of the cluster and despite
!     backing off, there are examples where all factorizations tried
!     (in IEEE mode, allowing zero pivots & infinities) have INFINITE
!     element growth.
!     For this reason, we should use PIVMIN in this subroutine so that at
!     least the L D L^T factorization exists. It can be checked afterwards
!     whether the element growth caused bad residuals/orthogonality.
 
!     Decide whether the code should accept the best among all
!     representations despite large element growth or signal INFO=1
!     Setting NOFAIL to .FALSE. for quick fix for bug 113
      nofail = .FALSE.
!
 
!     Compute the average gap length of the cluster
      clwdth = ABS(W(Clend)-W(Clstrt)) + Werr(Clend) + Werr(Clstrt)
      avgap = clwdth/DBLE(Clend-Clstrt)
      mingap = MIN(Clgapl,Clgapr)
!     Initial values for shifts to both ends of cluster
      lsigma = MIN(W(Clstrt),W(Clend)) - Werr(Clstrt)
      rsigma = MAX(W(Clstrt),W(Clend)) + Werr(Clend)
 
!     Use a small fudge to make sure that we really shift to the outside
      lsigma = lsigma - ABS(lsigma)*FOUR*eps
      rsigma = rsigma + ABS(rsigma)*FOUR*eps
 
!     Compute upper bounds for how much to back off the initial shifts
      ldmax = QUART*mingap + TWO*Pivmin
      rdmax = QUART*mingap + TWO*Pivmin
 
      ldelta = MAX(avgap,Wgap(Clstrt))/fact
      rdelta = MAX(avgap,Wgap(Clend-1))/fact
!
!     Initialize the record of the best representation found
!
      s = DLAMCH('S')
      smlgrowth = ONE/s
      fail = DBLE(N-1)*mingap/(Spdiam*eps)
      fail2 = DBLE(N-1)*mingap/(Spdiam*SQRT(eps))
      bestshift = lsigma
!
!     while (KTRY <= KTRYMAX)
      ktry = 0
      growthbound = MAXGROWTH1*Spdiam
 
 100  sawnan1 = .FALSE.
      sawnan2 = .FALSE.
!     Ensure that we do not back off too much of the initial shifts
      ldelta = MIN(ldmax,ldelta)
      rdelta = MIN(rdmax,rdelta)
 
!     Compute the element growth when shifting to both ends of the cluster
!     accept the shift if there is no element growth at one of the two ends
 
!     Left end
      s = -lsigma
      Dplus(1) = D(1) + s
      IF ( ABS(Dplus(1))<Pivmin ) THEN
         Dplus(1) = -Pivmin
!        Need to set SAWNAN1 because refined RRR test should not be used
!        in this case
         sawnan1 = .TRUE.
      ENDIF
      max1 = ABS(Dplus(1))
      DO i = 1 , N - 1
         Lplus(i) = Ld(i)/Dplus(i)
         s = s*Lplus(i)*L(i) - lsigma
         Dplus(i+1) = D(i+1) + s
         IF ( ABS(Dplus(i+1))<Pivmin ) THEN
            Dplus(i+1) = -Pivmin
!           Need to set SAWNAN1 because refined RRR test should not be used
!           in this case
            sawnan1 = .TRUE.
         ENDIF
         max1 = MAX(max1,ABS(Dplus(i+1)))
      ENDDO
      sawnan1 = sawnan1 .OR. DISNAN(max1)
 
      IF ( forcer .OR. (max1<=growthbound .AND. .NOT.sawnan1) ) THEN
         Sigma = lsigma
         shift = SLEFT
         GOTO 200
      ENDIF
 
!     Right end
      s = -rsigma
      Work(1) = D(1) + s
      IF ( ABS(Work(1))<Pivmin ) THEN
         Work(1) = -Pivmin
!        Need to set SAWNAN2 because refined RRR test should not be used
!        in this case
         sawnan2 = .TRUE.
      ENDIF
      max2 = ABS(Work(1))
      DO i = 1 , N - 1
         Work(N+i) = Ld(i)/Work(i)
         s = s*Work(N+i)*L(i) - rsigma
         Work(i+1) = D(i+1) + s
         IF ( ABS(Work(i+1))<Pivmin ) THEN
            Work(i+1) = -Pivmin
!           Need to set SAWNAN2 because refined RRR test should not be used
!           in this case
            sawnan2 = .TRUE.
         ENDIF
         max2 = MAX(max2,ABS(Work(i+1)))
      ENDDO
      sawnan2 = sawnan2 .OR. DISNAN(max2)
 
      IF ( forcer .OR. (max2<=growthbound .AND. .NOT.sawnan2) ) THEN
         Sigma = rsigma
         shift = SRIGHT
         GOTO 200
      ENDIF
!     If we are at this point, both shifts led to too much element growth
 
!     Record the better of the two shifts (provided it didn't lead to NaN)
!        both MAX1 and MAX2 are NaN
      IF ( .NOT.(sawnan1 .AND. sawnan2) ) THEN
         IF ( .NOT.sawnan1 ) THEN
            indx = 1
            IF ( max1<=smlgrowth ) THEN
               smlgrowth = max1
               bestshift = lsigma
            ENDIF
         ENDIF
         IF ( .NOT.sawnan2 ) THEN
            IF ( sawnan1 .OR. max2<=max1 ) indx = 2
            IF ( max2<=smlgrowth ) THEN
               smlgrowth = max2
               bestshift = rsigma
            ENDIF
         ENDIF
 
!     If we are here, both the left and the right shift led to
!     element growth. If the element growth is moderate, then
!     we may still accept the representation, if it passes a
!     refined test for RRR. This test supposes that no NaN occurred.
!     Moreover, we use the refined RRR test only for isolated clusters.
         IF ( (clwdth<mingap/DBLE(128)) .AND. (MIN(max1,max2)<fail2)    &
     &        .AND. (.NOT.sawnan1) .AND. (.NOT.sawnan2) ) THEN
            dorrr1 = .TRUE.
         ELSE
            dorrr1 = .FALSE.
         ENDIF
         tryrrr1 = .TRUE.
         IF ( tryrrr1 .AND. dorrr1 ) THEN
            IF ( indx==1 ) THEN
               tmp = ABS(Dplus(N))
               znm2 = ONE
               prod = ONE
               oldp = ONE
               DO i = N - 1 , 1 , -1
                  IF ( prod<=eps ) THEN
                     prod = ((Dplus(i+1)*Work(N+i+1))                   &
     &                      /(Dplus(i)*Work(N+i)))*oldp
                  ELSE
                     prod = prod*ABS(Work(N+i))
                  ENDIF
                  oldp = prod
                  znm2 = znm2 + prod**2
                  tmp = MAX(tmp,ABS(Dplus(i)*prod))
               ENDDO
               rrr1 = tmp/(Spdiam*SQRT(znm2))
               IF ( rrr1<=MAXGROWTH2 ) THEN
                  Sigma = lsigma
                  shift = SLEFT
                  GOTO 200
               ENDIF
            ELSEIF ( indx==2 ) THEN
               tmp = ABS(Work(N))
               znm2 = ONE
               prod = ONE
               oldp = ONE
               DO i = N - 1 , 1 , -1
                  IF ( prod<=eps ) THEN
                     prod = ((Work(i+1)*Lplus(i+1))/(Work(i)*Lplus(i))) &
     &                      *oldp
                  ELSE
                     prod = prod*ABS(Lplus(i))
                  ENDIF
                  oldp = prod
                  znm2 = znm2 + prod**2
                  tmp = MAX(tmp,ABS(Work(i)*prod))
               ENDDO
               rrr2 = tmp/(Spdiam*SQRT(znm2))
               IF ( rrr2<=MAXGROWTH2 ) THEN
                  Sigma = rsigma
                  shift = SRIGHT
                  GOTO 200
               ENDIF
            ENDIF
         ENDIF
      ENDIF
 
 
      IF ( ktry<KTRYMAX ) THEN
!        If we are here, both shifts failed also the RRR test.
!        Back off to the outside
         lsigma = MAX(lsigma-ldelta,lsigma-ldmax)
         rsigma = MIN(rsigma+rdelta,rsigma+rdmax)
         ldelta = TWO*ldelta
         rdelta = TWO*rdelta
         ktry = ktry + 1
         GOTO 100
!        None of the representations investigated satisfied our
!        criteria. Take the best one we found.
      ELSEIF ( (smlgrowth<fail) .OR. nofail ) THEN
         lsigma = bestshift
         rsigma = bestshift
         forcer = .TRUE.
         GOTO 100
      ELSE
         Info = 1
         RETURN
      ENDIF
 
 200  IF ( shift==SLEFT ) THEN
      ELSEIF ( shift==SRIGHT ) THEN
!        store new L and D back into DPLUS, LPLUS
         CALL DCOPY(N,Work,1,Dplus,1)
         CALL DCOPY(N-1,Work(N+1),1,Lplus,1)
      ENDIF
 
!
!     End of DLARRF
!
      END SUBROUTINE DLARRF
