!*==zlarrv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLARRV computes the eigenvectors of the tridiagonal matrix T = L D LT given L, D and the eigenvalues of L D LT.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLARRV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarrv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarrv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarrv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARRV( N, VL, VU, D, L, PIVMIN,
!                          ISPLIT, M, DOL, DOU, MINRGP,
!                          RTOL1, RTOL2, W, WERR, WGAP,
!                          IBLOCK, INDEXW, GERS, Z, LDZ, ISUPPZ,
!                          WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            DOL, DOU, INFO, LDZ, M, N
!       DOUBLE PRECISION   MINRGP, PIVMIN, RTOL1, RTOL2, VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            IBLOCK( * ), INDEXW( * ), ISPLIT( * ),
!      $                   ISUPPZ( * ), IWORK( * )
!       DOUBLE PRECISION   D( * ), GERS( * ), L( * ), W( * ), WERR( * ),
!      $                   WGAP( * ), WORK( * )
!       COMPLEX*16        Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARRV computes the eigenvectors of the tridiagonal matrix
!> T = L D L**T given L, D and APPROXIMATIONS to the eigenvalues of L D L**T.
!> The input eigenvalues should have been computed by DLARRE.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] VL
!> \verbatim
!>          VL is DOUBLE PRECISION
!>          Lower bound of the interval that contains the desired
!>          eigenvalues. VL < VU. Needed to compute gaps on the left or right
!>          end of the extremal eigenvalues in the desired RANGE.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>          VU is DOUBLE PRECISION
!>          Upper bound of the interval that contains the desired
!>          eigenvalues. VL < VU. Needed to compute gaps on the left or right
!>          end of the extremal eigenvalues in the desired RANGE.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the N diagonal elements of the diagonal matrix D.
!>          On exit, D may be overwritten.
!> \endverbatim
!>
!> \param[in,out] L
!> \verbatim
!>          L is DOUBLE PRECISION array, dimension (N)
!>          On entry, the (N-1) subdiagonal elements of the unit
!>          bidiagonal matrix L are in elements 1 to N-1 of L
!>          (if the matrix is not split.) At the end of each block
!>          is stored the corresponding shift as given by DLARRE.
!>          On exit, L is overwritten.
!> \endverbatim
!>
!> \param[in] PIVMIN
!> \verbatim
!>          PIVMIN is DOUBLE PRECISION
!>          The minimum pivot allowed in the Sturm sequence.
!> \endverbatim
!>
!> \param[in] ISPLIT
!> \verbatim
!>          ISPLIT is INTEGER array, dimension (N)
!>          The splitting points, at which T breaks up into blocks.
!>          The first block consists of rows/columns 1 to
!>          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
!>          through ISPLIT( 2 ), etc.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The total number of input eigenvalues.  0 <= M <= N.
!> \endverbatim
!>
!> \param[in] DOL
!> \verbatim
!>          DOL is INTEGER
!> \endverbatim
!>
!> \param[in] DOU
!> \verbatim
!>          DOU is INTEGER
!>          If the user wants to compute only selected eigenvectors from all
!>          the eigenvalues supplied, he can specify an index range DOL:DOU.
!>          Or else the setting DOL=1, DOU=M should be applied.
!>          Note that DOL and DOU refer to the order in which the eigenvalues
!>          are stored in W.
!>          If the user wants to compute only selected eigenpairs, then
!>          the columns DOL-1 to DOU+1 of the eigenvector space Z contain the
!>          computed eigenvectors. All other columns of Z are set to zero.
!> \endverbatim
!>
!> \param[in] MINRGP
!> \verbatim
!>          MINRGP is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] RTOL1
!> \verbatim
!>          RTOL1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] RTOL2
!> \verbatim
!>          RTOL2 is DOUBLE PRECISION
!>           Parameters for bisection.
!>           An interval [LEFT,RIGHT] has converged if
!>           RIGHT-LEFT < MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )
!> \endverbatim
!>
!> \param[in,out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          The first M elements of W contain the APPROXIMATE eigenvalues for
!>          which eigenvectors are to be computed.  The eigenvalues
!>          should be grouped by split-off block and ordered from
!>          smallest to largest within the block ( The output array
!>          W from DLARRE is expected here ). Furthermore, they are with
!>          respect to the shift of the corresponding root representation
!>          for their block. On exit, W holds the eigenvalues of the
!>          UNshifted matrix.
!> \endverbatim
!>
!> \param[in,out] WERR
!> \verbatim
!>          WERR is DOUBLE PRECISION array, dimension (N)
!>          The first M elements contain the semiwidth of the uncertainty
!>          interval of the corresponding eigenvalue in W
!> \endverbatim
!>
!> \param[in,out] WGAP
!> \verbatim
!>          WGAP is DOUBLE PRECISION array, dimension (N)
!>          The separation from the right neighbor eigenvalue in W.
!> \endverbatim
!>
!> \param[in] IBLOCK
!> \verbatim
!>          IBLOCK is INTEGER array, dimension (N)
!>          The indices of the blocks (submatrices) associated with the
!>          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue
!>          W(i) belongs to the first block from the top, =2 if W(i)
!>          belongs to the second block, etc.
!> \endverbatim
!>
!> \param[in] INDEXW
!> \verbatim
!>          INDEXW is INTEGER array, dimension (N)
!>          The indices of the eigenvalues within each block (submatrix);
!>          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the
!>          i-th eigenvalue W(i) is the 10-th eigenvalue in the second block.
!> \endverbatim
!>
!> \param[in] GERS
!> \verbatim
!>          GERS is DOUBLE PRECISION array, dimension (2*N)
!>          The N Gerschgorin intervals (the i-th Gerschgorin interval
!>          is (GERS(2*i-1), GERS(2*i)). The Gerschgorin intervals should
!>          be computed from the original UNshifted matrix.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ, max(1,M) )
!>          If INFO = 0, the first M columns of Z contain the
!>          orthonormal eigenvectors of the matrix T
!>          corresponding to the input eigenvalues, with the i-th
!>          column of Z holding the eigenvector associated with W(i).
!>          Note: the user must ensure that at least max(1,M) columns are
!>          supplied in the array Z.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          JOBZ = 'V', LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] ISUPPZ
!> \verbatim
!>          ISUPPZ is INTEGER array, dimension ( 2*max(1,M) )
!>          The support of the eigenvectors in Z, i.e., the indices
!>          indicating the nonzero elements in Z. The I-th eigenvector
!>          is nonzero only in elements ISUPPZ( 2*I-1 ) through
!>          ISUPPZ( 2*I ).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (12*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (7*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>
!>          > 0:  A problem occurred in ZLARRV.
!>          < 0:  One of the called subroutines signaled an internal problem.
!>                Needs inspection of the corresponding parameter IINFO
!>                for further information.
!>
!>          =-1:  Problem in DLARRB when refining a child's eigenvalues.
!>          =-2:  Problem in DLARRF when computing the RRR of a child.
!>                When a child is inside a tight cluster, it can be difficult
!>                to find an RRR. A partial remedy from the user's point of
!>                view is to make the parameter MINRGP smaller and recompile.
!>                However, as the orthogonality of the computed vectors is
!>                proportional to 1/MINRGP, the user should be aware that
!>                he might be trading in precision when he decreases MINRGP.
!>          =-3:  Problem in DLARRB when refining a single eigenvalue
!>                after the Rayleigh correction was rejected.
!>          = 5:  The Rayleigh Quotient Iteration failed to converge to
!>                full accuracy in MAXITR steps.
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
!> \ingroup complex16OTHERauxiliary
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
      SUBROUTINE ZLARRV(N,Vl,Vu,D,L,Pivmin,Isplit,M,Dol,Dou,Minrgp,     &
     &                  Rtol1,Rtol2,W,Werr,Wgap,Iblock,Indexw,Gers,Z,   &
     &                  Ldz,Isuppz,Work,Iwork,Info)
      IMPLICIT NONE
!*--ZLARRV288
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Dol , Dou , Info , Ldz , M , N
      DOUBLE PRECISION Minrgp , Pivmin , Rtol1 , Rtol2 , Vl , Vu
!     ..
!     .. Array Arguments ..
      INTEGER Iblock(*) , Indexw(*) , Isplit(*) , Isuppz(*) , Iwork(*)
      DOUBLE PRECISION D(*) , Gers(*) , L(*) , W(*) , Werr(*) , Wgap(*) &
     &                 , Work(*)
      COMPLEX*16 Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER MAXITR
      PARAMETER (MAXITR=10)
      COMPLEX*16 CZERO
      PARAMETER (CZERO=(0.0D0,0.0D0))
      DOUBLE PRECISION ZERO , ONE , TWO , THREE , FOUR , HALF
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0, &
     &           HALF=0.5D0)
!     ..
!     .. Local Scalars ..
      LOGICAL eskip , needbs , stp2ii , tryrqc , usedbs , usedrq
      INTEGER done , i , ibegin , idone , iend , ii , iindc1 , iindc2 , &
     &        iindr , iindwk , iinfo , im , in , indeig , indld ,       &
     &        indlld , indwrk , isupmn , isupmx , iter , itmp1 , j ,    &
     &        jblk , k , miniwsize , minwsize , nclus , ndepth ,        &
     &        negcnt , newcls , newfst , newftt , newlst , newsiz ,     &
     &        offset , oldcls , oldfst , oldien , oldlst , oldncl , p , &
     &        parity , q , wbegin , wend , windex , windmn , windpl ,   &
     &        zfrom , zto , zusedl , zusedu , zusedw
      INTEGER indin1 , indin2
      DOUBLE PRECISION bstres , bstw , eps , fudge , gap , gaptol , gl ,&
     &                 gu , lambda , left , lgap , mingma , nrminv ,    &
     &                 resid , rgap , right , rqcorr , rqtol , savgap , &
     &                 sgndef , sigma , spdiam , ssigma , tau , tmp ,   &
     &                 tol , ztz
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DLARRB , DLARRF , ZDSCAL , ZLAR1V , ZLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN
      INTRINSIC DCMPLX
!     ..
!     .. Executable Statements ..
!     ..
 
      Info = 0
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
!     The first N entries of WORK are reserved for the eigenvalues
      indld = N + 1
      indlld = 2*N + 1
      indin1 = 3*N + 1
      indin2 = 4*N + 1
      indwrk = 5*N + 1
      minwsize = 12*N
 
      DO i = 1 , minwsize
         Work(i) = ZERO
      ENDDO
 
!     IWORK(IINDR+1:IINDR+N) hold the twist indices R for the
!     factorization used to compute the FP vector
      iindr = 0
!     IWORK(IINDC1+1:IINC2+N) are used to store the clusters of the current
!     layer and the one above.
      iindc1 = N
      iindc2 = 2*N
      iindwk = 3*N + 1
 
      miniwsize = 7*N
      DO i = 1 , miniwsize
         Iwork(i) = 0
      ENDDO
 
      zusedl = 1
!        Set lower bound for use of Z
      IF ( Dol>1 ) zusedl = Dol - 1
      zusedu = M
!        Set lower bound for use of Z
      IF ( Dou<M ) zusedu = Dou + 1
!     The width of the part of Z that is used
      zusedw = zusedu - zusedl + 1
 
 
      CALL ZLASET('Full',N,zusedw,CZERO,CZERO,Z(1,zusedl),Ldz)
 
      eps = DLAMCH('Precision')
      rqtol = TWO*eps
!
!     Set expert flags for standard code.
      tryrqc = .TRUE.
 
      IF ( (Dol/=1) .OR. (Dou/=M) ) THEN
!        Only selected eigenpairs are computed. Since the other evalues
!        are not refined by RQ iteration, bisection has to compute to full
!        accuracy.
         Rtol1 = FOUR*eps
         Rtol2 = FOUR*eps
      ENDIF
 
!     The entries WBEGIN:WEND in W, WERR, WGAP correspond to the
!     desired eigenvalues. The support of the nonzero eigenvector
!     entries is contained in the interval IBEGIN:IEND.
!     Remark that if k eigenpairs are desired, then the eigenvectors
!     are stored in k contiguous columns of Z.
 
!     DONE is the number of eigenvectors already computed
      done = 0
      ibegin = 1
      wbegin = 1
      DO jblk = 1 , Iblock(M)
         iend = Isplit(jblk)
         sigma = L(iend)
!        Find the eigenvectors of the submatrix indexed IBEGIN
!        through IEND.
         wend = wbegin - 1
         DO
            IF ( wend<M ) THEN
               IF ( Iblock(wend+1)==jblk ) THEN
                  wend = wend + 1
                  CYCLE
               ENDIF
            ENDIF
            IF ( wend<wbegin ) THEN
               ibegin = iend + 1
               GOTO 100
            ELSEIF ( (wend<Dol) .OR. (wbegin>Dou) ) THEN
               ibegin = iend + 1
               wbegin = wend + 1
               GOTO 100
            ENDIF
 
!        Find local spectral diameter of the block
            gl = Gers(2*ibegin-1)
            gu = Gers(2*ibegin)
            DO i = ibegin + 1 , iend
               gl = MIN(Gers(2*i-1),gl)
               gu = MAX(Gers(2*i),gu)
            ENDDO
            spdiam = gu - gl
 
!        OLDIEN is the last index of the previous block
            oldien = ibegin - 1
!        Calculate the size of the current block
            in = iend - ibegin + 1
!        The number of eigenvalues in the current block
            im = wend - wbegin + 1
 
!        This is for a 1x1 block
            IF ( ibegin==iend ) THEN
               done = done + 1
               Z(ibegin,wbegin) = DCMPLX(ONE,ZERO)
               Isuppz(2*wbegin-1) = ibegin
               Isuppz(2*wbegin) = ibegin
               W(wbegin) = W(wbegin) + sigma
               Work(wbegin) = W(wbegin)
               ibegin = iend + 1
               wbegin = wbegin + 1
               GOTO 100
            ENDIF
 
!        The desired (shifted) eigenvalues are stored in W(WBEGIN:WEND)
!        Note that these can be approximations, in this case, the corresp.
!        entries of WERR give the size of the uncertainty interval.
!        The eigenvalue approximations will be refined when necessary as
!        high relative accuracy is required for the computation of the
!        corresponding eigenvectors.
            CALL DCOPY(im,W(wbegin),1,Work(wbegin),1)
 
!        We store in W the eigenvalue approximations w.r.t. the original
!        matrix T.
            DO i = 1 , im
               W(wbegin+i-1) = W(wbegin+i-1) + sigma
            ENDDO
 
 
!        NDEPTH is the current depth of the representation tree
            ndepth = 0
!        PARITY is either 1 or 0
            parity = 1
!        NCLUS is the number of clusters for the next level of the
!        representation tree, we start with NCLUS = 1 for the root
            nclus = 1
            Iwork(iindc1+1) = 1
            Iwork(iindc1+2) = im
 
!        IDONE is the number of eigenvectors already computed in the current
!        block
            idone = 0
            EXIT
         ENDDO
         DO
!        loop while( IDONE.LT.IM )
!        generate the representation tree for the current block and
!        compute the eigenvectors
            IF ( idone<im ) THEN
!           This is a crude protection against infinitely deep trees
               IF ( ndepth>M ) THEN
                  Info = -2
                  RETURN
               ENDIF
!           breadth first processing of the current level of the representation
!           tree: OLDNCL = number of clusters on current level
               oldncl = nclus
!           reset NCLUS to count the number of child clusters
               nclus = 0
!
               parity = 1 - parity
               IF ( parity==0 ) THEN
                  oldcls = iindc1
                  newcls = iindc2
               ELSE
                  oldcls = iindc2
                  newcls = iindc1
               ENDIF
!           Process the clusters on the current level
               DO i = 1 , oldncl
                  j = oldcls + 2*i
!              OLDFST, OLDLST = first, last index of current cluster.
!                               cluster indices start with 1 and are relative
!                               to WBEGIN when accessing W, WGAP, WERR, Z
                  oldfst = Iwork(j-1)
                  oldlst = Iwork(j)
                  IF ( ndepth>0 ) THEN
!                 Retrieve relatively robust representation (RRR) of cluster
!                 that has been computed at the previous level
!                 The RRR is stored in Z and overwritten once the eigenvectors
!                 have been computed or when the cluster is refined
 
                     IF ( (Dol==1) .AND. (Dou==M) ) THEN
!                    Get representation from location of the leftmost evalue
!                    of the cluster
                        j = wbegin + oldfst - 1
                     ELSEIF ( wbegin+oldfst-1<Dol ) THEN
!                       Get representation from the left end of Z array
                        j = Dol - 1
                     ELSEIF ( wbegin+oldfst-1>Dou ) THEN
!                       Get representation from the right end of Z array
                        j = Dou
                     ELSE
                        j = wbegin + oldfst - 1
                     ENDIF
                     DO k = 1 , in - 1
                        D(ibegin+k-1) = DBLE(Z(ibegin+k-1,j))
                        L(ibegin+k-1) = DBLE(Z(ibegin+k-1,j+1))
                     ENDDO
                     D(iend) = DBLE(Z(iend,j))
                     sigma = DBLE(Z(iend,j+1))
 
!                 Set the corresponding entries in Z to zero
                     CALL ZLASET('Full',in,2,CZERO,CZERO,Z(ibegin,j),   &
     &                           Ldz)
                  ENDIF
 
!              Compute DL and DLL of current RRR
                  DO j = ibegin , iend - 1
                     tmp = D(j)*L(j)
                     Work(indld-1+j) = tmp
                     Work(indlld-1+j) = tmp*L(j)
                  ENDDO
 
                  IF ( ndepth>0 ) THEN
!                 P and Q are index of the first and last eigenvalue to compute
!                 within the current block
                     p = Indexw(wbegin-1+oldfst)
                     q = Indexw(wbegin-1+oldlst)
!                 Offset for the arrays WORK, WGAP and WERR, i.e., the P-OFFSET
!                 through the Q-OFFSET elements of these arrays are to be used.
!                  OFFSET = P-OLDFST
                     offset = Indexw(wbegin) - 1
!                 perform limited bisection (if necessary) to get approximate
!                 eigenvalues to the precision needed.
                     CALL DLARRB(in,D(ibegin),Work(indlld+ibegin-1),p,q,&
     &                           Rtol1,Rtol2,offset,Work(wbegin),       &
     &                           Wgap(wbegin),Werr(wbegin),Work(indwrk),&
     &                           Iwork(iindwk),Pivmin,spdiam,in,iinfo)
                     IF ( iinfo/=0 ) THEN
                        Info = -1
                        RETURN
                     ENDIF
!                 We also recompute the extremal gaps. W holds all eigenvalues
!                 of the unshifted matrix and must be used for computation
!                 of WGAP, the entries of WORK might stem from RRRs with
!                 different shifts. The gaps from WBEGIN-1+OLDFST to
!                 WBEGIN-1+OLDLST are correctly computed in DLARRB.
!                 However, we only allow the gaps to become greater since
!                 this is what should happen when we decrease WERR
                     IF ( oldfst>1 ) Wgap(wbegin+oldfst-2)              &
     &                    = MAX(Wgap(wbegin+oldfst-2),W(wbegin+oldfst-1)&
     &                    -Werr(wbegin+oldfst-1)-W(wbegin+oldfst-2)     &
     &                    -Werr(wbegin+oldfst-2))
                     IF ( wbegin+oldlst-1<wend ) Wgap(wbegin+oldlst-1)  &
     &                    = MAX(Wgap(wbegin+oldlst-1),W(wbegin+oldlst)  &
     &                    -Werr(wbegin+oldlst)-W(wbegin+oldlst-1)       &
     &                    -Werr(wbegin+oldlst-1))
!                 Each time the eigenvalues in WORK get refined, we store
!                 the newly found approximation with all shifts applied in W
                     DO j = oldfst , oldlst
                        W(wbegin+j-1) = Work(wbegin+j-1) + sigma
                     ENDDO
                  ENDIF
 
!              Process the current node.
                  newfst = oldfst
                  DO j = oldfst , oldlst
                     IF ( j==oldlst ) THEN
!                    we are at the right end of the cluster, this is also the
!                    boundary of the child cluster
                        newlst = j
                     ELSEIF ( Wgap(wbegin+j-1)                          &
     &                        >=Minrgp*ABS(Work(wbegin+j-1)) ) THEN
!                    the right relative gap is big enough, the child cluster
!                    (NEWFST,..,NEWLST) is well separated from the following
                        newlst = j
                     ELSE
!                    inside a child cluster, the relative gap is not
!                    big enough.
                        CYCLE
                     ENDIF
 
!                 Compute size of child cluster found
                     newsiz = newlst - newfst + 1
 
!                 NEWFTT is the place in Z where the new RRR or the computed
!                 eigenvector is to be stored
                     IF ( (Dol==1) .AND. (Dou==M) ) THEN
!                    Store representation at location of the leftmost evalue
!                    of the cluster
                        newftt = wbegin + newfst - 1
                     ELSEIF ( wbegin+newfst-1<Dol ) THEN
!                       Store representation at the left end of Z array
                        newftt = Dol - 1
                     ELSEIF ( wbegin+newfst-1>Dou ) THEN
!                       Store representation at the right end of Z array
                        newftt = Dou
                     ELSE
                        newftt = wbegin + newfst - 1
                     ENDIF
 
                     IF ( newsiz>1 ) THEN
!
!                    Current child is not a singleton but a cluster.
!                    Compute and store new representation of child.
!
!
!                    Compute left and right cluster gap.
!
!                    LGAP and RGAP are not computed from WORK because
!                    the eigenvalue approximations may stem from RRRs
!                    different shifts. However, W hold all eigenvalues
!                    of the unshifted matrix. Still, the entries in WGAP
!                    have to be computed from WORK since the entries
!                    in W might be of the same order so that gaps are not
!                    exhibited correctly for very close eigenvalues.
                        IF ( newfst==1 ) THEN
                           lgap = MAX(ZERO,W(wbegin)-Werr(wbegin)-Vl)
                        ELSE
                           lgap = Wgap(wbegin+newfst-2)
                        ENDIF
                        rgap = Wgap(wbegin+newlst-1)
!
!                    Compute left- and rightmost eigenvalue of child
!                    to high precision in order to shift as close
!                    as possible and obtain as large relative gaps
!                    as possible
!
                        DO k = 1 , 2
                           IF ( k==1 ) THEN
                              p = Indexw(wbegin-1+newfst)
                           ELSE
                              p = Indexw(wbegin-1+newlst)
                           ENDIF
                           offset = Indexw(wbegin) - 1
                           CALL DLARRB(in,D(ibegin),                    &
     &                                 Work(indlld+ibegin-1),p,p,rqtol, &
     &                                 rqtol,offset,Work(wbegin),       &
     &                                 Wgap(wbegin),Werr(wbegin),       &
     &                                 Work(indwrk),Iwork(iindwk),      &
     &                                 Pivmin,spdiam,in,iinfo)
                        ENDDO
!
                        IF ( (wbegin+newlst-1<Dol) .OR.                 &
     &                       (wbegin+newfst-1>Dou) ) THEN
!                       if the cluster contains no desired eigenvalues
!                       skip the computation of that branch of the rep. tree
!
!                       We could skip before the refinement of the extremal
!                       eigenvalues of the child, but then the representation
!                       tree could be different from the one when nothing is
!                       skipped. For this reason we skip at this place.
                           idone = idone + newlst - newfst + 1
                           GOTO 4
                        ENDIF
!
!                    Compute RRR of child cluster.
!                    Note that the new RRR is stored in Z
!
!                    DLARRF needs LWORK = 2*N
                        CALL DLARRF(in,D(ibegin),L(ibegin),             &
     &                              Work(indld+ibegin-1),newfst,newlst, &
     &                              Work(wbegin),Wgap(wbegin),          &
     &                              Werr(wbegin),spdiam,lgap,rgap,      &
     &                              Pivmin,tau,Work(indin1),Work(indin2)&
     &                              ,Work(indwrk),iinfo)
!                    In the complex case, DLARRF cannot write
!                    the new RRR directly into Z and needs an intermediate
!                    workspace
                        DO k = 1 , in - 1
                           Z(ibegin+k-1,newftt)                         &
     &                        = DCMPLX(Work(indin1+k-1),ZERO)
                           Z(ibegin+k-1,newftt+1)                       &
     &                        = DCMPLX(Work(indin2+k-1),ZERO)
                        ENDDO
                        Z(iend,newftt) = DCMPLX(Work(indin1+in-1),ZERO)
                        IF ( iinfo==0 ) THEN
!                       a new RRR for the cluster was found by DLARRF
!                       update shift and store it
                           ssigma = sigma + tau
                           Z(iend,newftt+1) = DCMPLX(ssigma,ZERO)
!                       WORK() are the midpoints and WERR() the semi-width
!                       Note that the entries in W are unchanged.
                           DO k = newfst , newlst
                              fudge = THREE*eps*ABS(Work(wbegin+k-1))
                              Work(wbegin+k-1) = Work(wbegin+k-1) - tau
                              fudge = fudge +                           &
     &                                FOUR*eps*ABS(Work(wbegin+k-1))
!                          Fudge errors
                              Werr(wbegin+k-1) = Werr(wbegin+k-1)       &
     &                           + fudge
!                          Gaps are not fudged. Provided that WERR is small
!                          when eigenvalues are close, a zero gap indicates
!                          that a new representation is needed for resolving
!                          the cluster. A fudge could lead to a wrong decision
!                          of judging eigenvalues 'separated' which in
!                          reality are not. This could have a negative impact
!                          on the orthogonality of the computed eigenvectors.
                           ENDDO
 
                           nclus = nclus + 1
                           k = newcls + 2*nclus
                           Iwork(k-1) = newfst
                           Iwork(k) = newlst
                        ELSE
                           Info = -2
                           RETURN
                        ENDIF
                     ELSE
!
!                    Compute eigenvector of singleton
!
                        iter = 0
!
                        tol = FOUR*LOG(DBLE(in))*eps
!
                        k = newfst
                        windex = wbegin + k - 1
                        windmn = MAX(windex-1,1)
                        windpl = MIN(windex+1,M)
                        lambda = Work(windex)
                        done = done + 1
!                    Check if eigenvector computation is to be skipped
                        IF ( (windex<Dol) .OR. (windex>Dou) ) THEN
                           eskip = .TRUE.
                           GOTO 2
                        ELSE
                           eskip = .FALSE.
                        ENDIF
                        left = Work(windex) - Werr(windex)
                        right = Work(windex) + Werr(windex)
                        indeig = Indexw(windex)
!                    Note that since we compute the eigenpairs for a child,
!                    all eigenvalue approximations are w.r.t the same shift.
!                    In this case, the entries in WORK should be used for
!                    computing the gaps since they exhibit even very small
!                    differences in the eigenvalues, as opposed to the
!                    entries in W which might "look" the same.
 
                        IF ( k==1 ) THEN
!                       In the case RANGE='I' and with not much initial
!                       accuracy in LAMBDA and VL, the formula
!                       LGAP = MAX( ZERO, (SIGMA - VL) + LAMBDA )
!                       can lead to an overestimation of the left gap and
!                       thus to inadequately early RQI 'convergence'.
!                       Prevent this by forcing a small left gap.
                           lgap = eps*MAX(ABS(left),ABS(right))
                        ELSE
                           lgap = Wgap(windmn)
                        ENDIF
                        IF ( k==im ) THEN
!                       In the case RANGE='I' and with not much initial
!                       accuracy in LAMBDA and VU, the formula
!                       can lead to an overestimation of the right gap and
!                       thus to inadequately early RQI 'convergence'.
!                       Prevent this by forcing a small right gap.
                           rgap = eps*MAX(ABS(left),ABS(right))
                        ELSE
                           rgap = Wgap(windex)
                        ENDIF
                        gap = MIN(lgap,rgap)
                        IF ( (k==1) .OR. (k==im) ) THEN
!                       The eigenvector support can become wrong
!                       because significant entries could be cut off due to a
!                       large GAPTOL parameter in LAR1V. Prevent this.
                           gaptol = ZERO
                        ELSE
                           gaptol = gap*eps
                        ENDIF
                        isupmn = in
                        isupmx = 1
!                    Update WGAP so that it holds the minimum gap
!                    to the left or the right. This is crucial in the
!                    case where bisection is used to ensure that the
!                    eigenvalue is refined up to the required precision.
!                    The correct value is restored afterwards.
                        savgap = Wgap(windex)
                        Wgap(windex) = gap
!                    We want to use the Rayleigh Quotient Correction
!                    as often as possible since it converges quadratically
!                    when we are close enough to the desired eigenvalue.
!                    However, the Rayleigh Quotient can have the wrong sign
!                    and lead us away from the desired eigenvalue. In this
!                    case, the best we can do is to use bisection.
                        usedbs = .FALSE.
                        usedrq = .FALSE.
!                    Bisection is initially turned off unless it is forced
                        needbs = .NOT.tryrqc
                        DO
!                    Check if bisection should be used to refine eigenvalue
                           IF ( needbs ) THEN
!                       Take the bisection as new iterate
                              usedbs = .TRUE.
                              itmp1 = Iwork(iindr+windex)
                              offset = Indexw(wbegin) - 1
                              CALL DLARRB(in,D(ibegin),                 &
     &                           Work(indlld+ibegin-1),indeig,indeig,   &
     &                           ZERO,TWO*eps,offset,Work(wbegin),      &
     &                           Wgap(wbegin),Werr(wbegin),Work(indwrk),&
     &                           Iwork(iindwk),Pivmin,spdiam,itmp1,     &
     &                           iinfo)
                              IF ( iinfo/=0 ) THEN
                                 Info = -3
                                 RETURN
                              ENDIF
                              lambda = Work(windex)
!                       Reset twist index from inaccurate LAMBDA to
!                       force computation of true MINGMA
                              Iwork(iindr+windex) = 0
                           ENDIF
!                    Given LAMBDA, compute the eigenvector.
                           CALL ZLAR1V(in,1,in,lambda,D(ibegin),        &
     &                                 L(ibegin),Work(indld+ibegin-1),  &
     &                                 Work(indlld+ibegin-1),Pivmin,    &
     &                                 gaptol,Z(ibegin,windex),         &
     &                                 .NOT.usedbs,negcnt,ztz,mingma,   &
     &                                 Iwork(iindr+windex),             &
     &                                 Isuppz(2*windex-1),nrminv,resid, &
     &                                 rqcorr,Work(indwrk))
                           IF ( iter==0 ) THEN
                              bstres = resid
                              bstw = lambda
                           ELSEIF ( resid<bstres ) THEN
                              bstres = resid
                              bstw = lambda
                           ENDIF
                           isupmn = MIN(isupmn,Isuppz(2*windex-1))
                           isupmx = MAX(isupmx,Isuppz(2*windex))
                           iter = iter + 1
 
!                    sin alpha <= |resid|/gap
!                    Note that both the residual and the gap are
!                    proportional to the matrix, so ||T|| doesn't play
!                    a role in the quotient
 
!
!                    Convergence test for Rayleigh-Quotient iteration
!                    (omitted when Bisection has been used)
!
                           IF ( resid>tol*gap .AND. ABS(rqcorr)         &
     &                          >rqtol*ABS(lambda) .AND. .NOT.usedbs )  &
     &                          THEN
!                       We need to check that the RQCORR update doesn't
!                       move the eigenvalue away from the desired one and
!                       towards a neighbor. -> protection with bisection
                              IF ( indeig<=negcnt ) THEN
!                          The wanted eigenvalue lies to the left
                                 sgndef = -ONE
                              ELSE
!                          The wanted eigenvalue lies to the right
                                 sgndef = ONE
                              ENDIF
!                       We only use the RQCORR if it improves the
!                       the iterate reasonably.
                              IF ( (rqcorr*sgndef>=ZERO) .AND.          &
     &                             (lambda+rqcorr<=right) .AND.         &
     &                             (lambda+rqcorr>=left) ) THEN
                                 usedrq = .TRUE.
!                          Store new midpoint of bisection interval in WORK
                                 IF ( sgndef==ONE ) THEN
!                             The current LAMBDA is on the left of the true
!                             eigenvalue
                                    left = lambda
!                             We prefer to assume that the error estimate
!                             is correct. We could make the interval not
!                             as a bracket but to be modified if the RQCORR
!                             chooses to. In this case, the RIGHT side should
!                             be modified as follows:
!                              RIGHT = MAX(RIGHT, LAMBDA + RQCORR)
                                 ELSE
!                             The current LAMBDA is on the right of the true
!                             eigenvalue
                                    right = lambda
!                             See comment about assuming the error estimate is
!                             correct above.
!                              LEFT = MIN(LEFT, LAMBDA + RQCORR)
                                 ENDIF
                                 Work(windex) = HALF*(right+left)
!                          Take RQCORR since it has the correct sign and
!                          improves the iterate reasonably
                                 lambda = lambda + rqcorr
!                          Update width of error interval
                                 Werr(windex) = HALF*(right-left)
                              ELSE
                                 needbs = .TRUE.
                              ENDIF
                              IF ( right-left<rqtol*ABS(lambda) ) THEN
!                             The eigenvalue is computed to bisection accuracy
!                             compute eigenvector and stop
                                 usedbs = .TRUE.
                                 CYCLE
                              ELSEIF ( iter<MAXITR ) THEN
                                 CYCLE
                              ELSEIF ( iter==MAXITR ) THEN
                                 needbs = .TRUE.
                                 CYCLE
                              ELSE
                                 Info = 5
                                 RETURN
                              ENDIF
                           ELSE
                              stp2ii = .FALSE.
                              IF ( usedrq .AND. usedbs .AND.            &
     &                             bstres<=resid ) THEN
                                 lambda = bstw
                                 stp2ii = .TRUE.
                              ENDIF
!                          improve error angle by second step
                              IF ( stp2ii )                             &
     &                             CALL ZLAR1V(in,1,in,lambda,D(ibegin),&
     &                             L(ibegin),Work(indld+ibegin-1),      &
     &                             Work(indlld+ibegin-1),Pivmin,gaptol, &
     &                             Z(ibegin,windex),.NOT.usedbs,negcnt, &
     &                             ztz,mingma,Iwork(iindr+windex),      &
     &                             Isuppz(2*windex-1),nrminv,resid,     &
     &                             rqcorr,Work(indwrk))
                              Work(windex) = lambda
                           ENDIF
!
!                    Compute FP-vector support w.r.t. whole matrix
!
                           Isuppz(2*windex-1) = Isuppz(2*windex-1)      &
     &                        + oldien
                           Isuppz(2*windex) = Isuppz(2*windex) + oldien
                           zfrom = Isuppz(2*windex-1)
                           zto = Isuppz(2*windex)
                           isupmn = isupmn + oldien
                           isupmx = isupmx + oldien
!                    Ensure vector is ok if support in the RQI has changed
                           IF ( isupmn<zfrom ) THEN
                              DO ii = isupmn , zfrom - 1
                                 Z(ii,windex) = ZERO
                              ENDDO
                           ENDIF
                           IF ( isupmx>zto ) THEN
                              DO ii = zto + 1 , isupmx
                                 Z(ii,windex) = ZERO
                              ENDDO
                           ENDIF
                           CALL ZDSCAL(zto-zfrom+1,nrminv,              &
     &                                 Z(zfrom,windex),1)
                           EXIT
                        ENDDO
!                    Update W
 2                      W(windex) = lambda + sigma
!                    Recompute the gaps on the left and right
!                    But only allow them to become larger and not
!                    smaller (which can only happen through "bad"
!                    cancellation and doesn't reflect the theory
!                    where the initial gaps are underestimated due
!                    to WERR being too crude.)
                        IF ( .NOT.eskip ) THEN
                           IF ( k>1 ) Wgap(windmn)                      &
     &                          = MAX(Wgap(windmn),W(windex)            &
     &                          -Werr(windex)-W(windmn)-Werr(windmn))
                           IF ( windex<wend ) Wgap(windex)              &
     &                          = MAX(savgap,W(windpl)-Werr(windpl)     &
     &                          -W(windex)-Werr(windex))
                        ENDIF
                        idone = idone + 1
                     ENDIF
!                 here ends the code for the current child
!
!                 Proceed to any remaining child nodes
 4                   newfst = j + 1
                  ENDDO
               ENDDO
               ndepth = ndepth + 1
               CYCLE
            ENDIF
            ibegin = iend + 1
            wbegin = wend + 1
            EXIT
         ENDDO
 100  ENDDO
!
 
!
!     End of ZLARRV
!
      END SUBROUTINE ZLARRV
