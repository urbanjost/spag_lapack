!*==cgsvj1.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGSVJ1 pre-processor for the routine cgesvj, applies Jacobi rotations targeting only particular pivots.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGSVJ1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgsvj1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgsvj1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgsvj1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV,
!                          EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       REAL               EPS, SFMIN, TOL
!       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP
!       CHARACTER*1        JOBV
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), D( N ), V( LDV, * ), WORK( LWORK )
!       REAL               SVA( N )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGSVJ1 is called from CGESVJ as a pre-processor and that is its main
!> purpose. It applies Jacobi rotations in the same way as CGESVJ does, but
!> it targets only particular pivots and it does not check convergence
!> (stopping criterion). Few tuning parameters (marked by [TP]) are
!> available for the implementer.
!>
!> Further Details
!> ~~~~~~~~~~~~~~~
!> CGSVJ1 applies few sweeps of Jacobi rotations in the column space of
!> the input M-by-N matrix A. The pivot pairs are taken from the (1,2)
!> off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The
!> block-entries (tiles) of the (1,2) off-diagonal block are marked by the
!> [x]'s in the following scheme:
!>
!>    | *  *  * [x] [x] [x]|
!>    | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
!>    | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
!>    |[x] [x] [x] *  *  * |
!>    |[x] [x] [x] *  *  * |
!>    |[x] [x] [x] *  *  * |
!>
!> In terms of the columns of A, the first N1 columns are rotated 'against'
!> the remaining N-N1 columns, trying to increase the angle between the
!> corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is
!> tiled using quadratic tiles of side KBL. Here, KBL is a tuning parameter.
!> The number of sweeps is given in NSWEEP and the orthogonality threshold
!> is given in TOL.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBV
!> \verbatim
!>          JOBV is CHARACTER*1
!>          Specifies whether the output from this procedure is used
!>          to compute the matrix V:
!>          = 'V': the product of the Jacobi rotations is accumulated
!>                 by postmulyiplying the N-by-N array V.
!>                (See the description of V.)
!>          = 'A': the product of the Jacobi rotations is accumulated
!>                 by postmulyiplying the MV-by-N array V.
!>                (See the descriptions of MV and V.)
!>          = 'N': the Jacobi rotations are not accumulated.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the input matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the input matrix A.
!>          M >= N >= 0.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>          N1 specifies the 2 x 2 block partition, the first N1 columns are
!>          rotated 'against' the remaining N-N1 columns of A.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, M-by-N matrix A, such that A*diag(D) represents
!>          the input matrix.
!>          On exit,
!>          A_onexit * D_onexit represents the input matrix A*diag(D)
!>          post-multiplied by a sequence of Jacobi rotations, where the
!>          rotation threshold and the total number of sweeps are given in
!>          TOL and NSWEEP, respectively.
!>          (See the descriptions of N1, D, TOL and NSWEEP.)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is COMPLEX array, dimension (N)
!>          The array D accumulates the scaling factors from the fast scaled
!>          Jacobi rotations.
!>          On entry, A*diag(D) represents the input matrix.
!>          On exit, A_onexit*diag(D_onexit) represents the input matrix
!>          post-multiplied by a sequence of Jacobi rotations, where the
!>          rotation threshold and the total number of sweeps are given in
!>          TOL and NSWEEP, respectively.
!>          (See the descriptions of N1, A, TOL and NSWEEP.)
!> \endverbatim
!>
!> \param[in,out] SVA
!> \verbatim
!>          SVA is REAL array, dimension (N)
!>          On entry, SVA contains the Euclidean norms of the columns of
!>          the matrix A*diag(D).
!>          On exit, SVA contains the Euclidean norms of the columns of
!>          the matrix onexit*diag(D_onexit).
!> \endverbatim
!>
!> \param[in] MV
!> \verbatim
!>          MV is INTEGER
!>          If JOBV = 'A', then MV rows of V are post-multipled by a
!>                           sequence of Jacobi rotations.
!>          If JOBV = 'N',   then MV is not referenced.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is COMPLEX array, dimension (LDV,N)
!>          If JOBV = 'V' then N rows of V are post-multipled by a
!>                           sequence of Jacobi rotations.
!>          If JOBV = 'A' then MV rows of V are post-multipled by a
!>                           sequence of Jacobi rotations.
!>          If JOBV = 'N',   then V is not referenced.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V,  LDV >= 1.
!>          If JOBV = 'V', LDV >= N.
!>          If JOBV = 'A', LDV >= MV.
!> \endverbatim
!>
!> \param[in] EPS
!> \verbatim
!>          EPS is REAL
!>          EPS = SLAMCH('Epsilon')
!> \endverbatim
!>
!> \param[in] SFMIN
!> \verbatim
!>          SFMIN is REAL
!>          SFMIN = SLAMCH('Safe Minimum')
!> \endverbatim
!>
!> \param[in] TOL
!> \verbatim
!>          TOL is REAL
!>          TOL is the threshold for Jacobi rotations. For a pair
!>          A(:,p), A(:,q) of pivot columns, the Jacobi rotation is
!>          applied only if ABS(COS(angle(A(:,p),A(:,q)))) > TOL.
!> \endverbatim
!>
!> \param[in] NSWEEP
!> \verbatim
!>          NSWEEP is INTEGER
!>          NSWEEP is the number of sweeps of Jacobi rotations to be
!>          performed.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>         WORK is COMPLEX array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          LWORK is the dimension of WORK. LWORK >= M.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, then the i-th argument had an illegal value
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
!> \ingroup complexOTHERcomputational
!
!> \par Contributor:
!  ==================
!>
!> Zlatko Drmac (Zagreb, Croatia)
!
!  =====================================================================
      SUBROUTINE CGSVJ1(Jobv,M,N,N1,A,Lda,D,Sva,Mv,V,Ldv,Eps,Sfmin,Tol, &
     &                  Nsweep,Work,Lwork,Info)
      USE S_CAXPY
      USE S_CCOPY
      USE S_CDOTC
      USE S_CLASCL
      USE S_CLASSQ
      USE S_CROT
      USE S_CSWAP
      USE S_ISAMAX
      USE S_LSAME
      USE S_SCNRM2
      USE S_XERBLA
      IMPLICIT NONE
!*--CGSVJ1251
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , HALF = 0.5E0 , ONE = 1.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER(1) :: Jobv
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: N1
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(N) :: D
      REAL , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      REAL , INTENT(IN) :: Eps
      REAL , INTENT(IN) :: Sfmin
      REAL , INTENT(IN) :: Tol
      INTEGER , INTENT(IN) :: Nsweep
      COMPLEX , DIMENSION(Lwork) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: aapp , aapp0 , aapq1 , aaqq , apoaq , aqoap , big ,       &
     &        bigtheta , cs , mxaapq , mxsinj , rootbig , rooteps ,     &
     &        rootsfmin , roottol , small , sn , t , temp1 , theta ,    &
     &        thsign
      COMPLEX :: aapq , ompq
      LOGICAL :: applv , rotok , rsvec
      INTEGER :: blskip , emptsw , i , ibr , ierr , igl , ijblsk ,      &
     &           iswrot , jbc , jgl , kbl , mvl , nblc , nblr , notrot ,&
     &           p , pskipped , q , rowskip , swband
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     .. from BLAS
!     .. from LAPACK
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      applv = LSAME(Jobv,'A')
      rsvec = LSAME(Jobv,'V')
      IF ( .NOT.(rsvec .OR. applv .OR. LSAME(Jobv,'N')) ) THEN
         Info = -1
      ELSEIF ( M<0 ) THEN
         Info = -2
      ELSEIF ( (N<0) .OR. (N>M) ) THEN
         Info = -3
      ELSEIF ( N1<0 ) THEN
         Info = -4
      ELSEIF ( Lda<M ) THEN
         Info = -6
      ELSEIF ( (rsvec .OR. applv) .AND. (Mv<0) ) THEN
         Info = -9
      ELSEIF ( (rsvec .AND. (Ldv<N)) .OR. (applv .AND. (Ldv<Mv)) ) THEN
         Info = -11
      ELSEIF ( Tol<=Eps ) THEN
         Info = -14
      ELSEIF ( Nsweep<0 ) THEN
         Info = -15
      ELSEIF ( Lwork<M ) THEN
         Info = -17
      ELSE
         Info = 0
      ENDIF
!
!     #:(
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGSVJ1',-Info)
         RETURN
      ENDIF
!
      IF ( rsvec ) THEN
         mvl = N
      ELSEIF ( applv ) THEN
         mvl = Mv
      ENDIF
      rsvec = rsvec .OR. applv
 
      rooteps = SQRT(Eps)
      rootsfmin = SQRT(Sfmin)
      small = Sfmin/Eps
      big = ONE/Sfmin
      rootbig = ONE/rootsfmin
!     LARGE = BIG / SQRT( REAL( M*N ) )
      bigtheta = ONE/rooteps
      roottol = SQRT(Tol)
!
!     .. Initialize the right singular vector matrix ..
!
!     RSVEC = LSAME( JOBV, 'Y' )
!
      emptsw = N1*(N-N1)
      notrot = 0
!
!     .. Row-cyclic pivot strategy with de Rijk's pivoting ..
!
      kbl = MIN(8,N)
      nblr = N1/kbl
      IF ( (nblr*kbl)/=N1 ) nblr = nblr + 1
 
!     .. the tiling is nblr-by-nblc [tiles]
 
      nblc = (N-N1)/kbl
      IF ( (nblc*kbl)/=(N-N1) ) nblc = nblc + 1
      blskip = (kbl**2) + 1
![TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.
 
      rowskip = MIN(5,kbl)
![TP] ROWSKIP is a tuning parameter.
      swband = 0
![TP] SWBAND is a tuning parameter. It is meaningful and effective
!     if CGESVJ is used as a computational routine in the preconditioned
!     Jacobi SVD algorithm CGEJSV.
!
!
!     | *   *   * [x] [x] [x]|
!     | *   *   * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
!     | *   *   * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
!     |[x] [x] [x] *   *   * |
!     |[x] [x] [x] *   *   * |
!     |[x] [x] [x] *   *   * |
!
!
      DO i = 1 , Nsweep
!
!     .. go go go ...
!
         mxaapq = ZERO
         mxsinj = ZERO
         iswrot = 0
!
         notrot = 0
         pskipped = 0
!
!     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs
!     1 <= p < q <= N. This is the first step toward a blocked implementation
!     of the rotations. New implementation, based on block transformations,
!     is under development.
!
         DO ibr = 1 , nblr
!
            igl = (ibr-1)*kbl + 1
!
 
!
! ... go to the off diagonal blocks
!
            igl = (ibr-1)*kbl + 1
!
!            DO 2010 jbc = ibr + 1, NBL
            DO jbc = 1 , nblc
!
               jgl = (jbc-1)*kbl + N1 + 1
!
!        doing the block at ( ibr, jbc )
!
               ijblsk = 0
               DO p = igl , MIN(igl+kbl-1,N1)
!
                  aapp = Sva(p)
                  IF ( aapp>ZERO ) THEN
!
                     pskipped = 0
!
                     DO q = jgl , MIN(jgl+kbl-1,N)
!
                        aaqq = Sva(q)
                        IF ( aaqq>ZERO ) THEN
                           aapp0 = aapp
!
!     .. M x 2 Jacobi SVD ..
!
!        Safe Gram matrix computation
!
                           IF ( aaqq>=ONE ) THEN
                              IF ( aapp>=aaqq ) THEN
                                 rotok = (small*aapp)<=aaqq
                              ELSE
                                 rotok = (small*aaqq)<=aapp
                              ENDIF
                              IF ( aapp<(big/aaqq) ) THEN
                                 aapq = (CDOTC(M,A(1,p),1,A(1,q),1)     &
     &                                  /aaqq)/aapp
                              ELSE
                                 CALL CCOPY(M,A(1,p),1,Work,1)
                                 CALL CLASCL('G',0,0,aapp,ONE,M,1,Work, &
     &                              Lda,ierr)
                                 aapq = CDOTC(M,Work,1,A(1,q),1)/aaqq
                              ENDIF
                           ELSE
                              IF ( aapp>=aaqq ) THEN
                                 rotok = aapp<=(aaqq/small)
                              ELSE
                                 rotok = aaqq<=(aapp/small)
                              ENDIF
                              IF ( aapp>(small/aaqq) ) THEN
                                 aapq = (CDOTC(M,A(1,p),1,A(1,q),1)     &
     &                                  /MAX(aaqq,aapp))/MIN(aaqq,aapp)
                              ELSE
                                 CALL CCOPY(M,A(1,q),1,Work,1)
                                 CALL CLASCL('G',0,0,aaqq,ONE,M,1,Work, &
     &                              Lda,ierr)
                                 aapq = CDOTC(M,A(1,p),1,Work,1)/aapp
                              ENDIF
                           ENDIF
!
!                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q)
                           aapq1 = -ABS(aapq)
                           mxaapq = MAX(mxaapq,-aapq1)
!
!        TO rotate or NOT to rotate, THAT is the question ...
!
                           IF ( ABS(aapq1)>Tol ) THEN
                              ompq = aapq/ABS(aapq)
                              notrot = 0
![RTD]      ROTATED  = ROTATED + 1
                              pskipped = 0
                              iswrot = iswrot + 1
!
                              IF ( rotok ) THEN
!
                                 aqoap = aaqq/aapp
                                 apoaq = aapp/aaqq
                                 theta = -HALF*ABS(aqoap-apoaq)/aapq1
                                 IF ( aaqq>aapp0 ) theta = -theta
!
                                 IF ( ABS(theta)>bigtheta ) THEN
                                    t = HALF/theta
                                    cs = ONE
                                    CALL CROT(M,A(1,p),1,A(1,q),1,cs,   &
     &                                 CONJG(ompq)*t)
                                    IF ( rsvec )                        &
     &                                 CALL CROT(mvl,V(1,p),1,V(1,q),1, &
     &                                 cs,CONJG(ompq)*t)
                                    Sva(q)                              &
     &                                 = aaqq*SQRT(MAX(ZERO,ONE+t*apoaq*&
     &                                 aapq1))
                                    aapp = aapp*SQRT                    &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq1))
                                    mxsinj = MAX(mxsinj,ABS(t))
                                 ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                    thsign = -SIGN(ONE,aapq1)
                                    IF ( aaqq>aapp0 ) thsign = -thsign
                                    t = ONE/                            &
     &                                  (theta+thsign*SQRT(ONE+theta*   &
     &                                  theta))
                                    cs = SQRT(ONE/(ONE+t*t))
                                    sn = t*cs
                                    mxsinj = MAX(mxsinj,ABS(sn))
                                    Sva(q)                              &
     &                                 = aaqq*SQRT(MAX(ZERO,ONE+t*apoaq*&
     &                                 aapq1))
                                    aapp = aapp*SQRT                    &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq1))
!
                                    CALL CROT(M,A(1,p),1,A(1,q),1,cs,   &
     &                                 CONJG(ompq)*sn)
                                    IF ( rsvec )                        &
     &                                 CALL CROT(mvl,V(1,p),1,V(1,q),1, &
     &                                 cs,CONJG(ompq)*sn)
                                 ENDIF
                                 D(p) = -D(q)*ompq
!
!              .. have to use modified Gram-Schmidt like transformation
                              ELSEIF ( aapp>aaqq ) THEN
                                 CALL CCOPY(M,A(1,p),1,Work,1)
                                 CALL CLASCL('G',0,0,aapp,ONE,M,1,Work, &
     &                              Lda,ierr)
                                 CALL CLASCL('G',0,0,aaqq,ONE,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 CALL CAXPY(M,-aapq,Work,1,A(1,q),1)
                                 CALL CLASCL('G',0,0,ONE,aaqq,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 Sva(q) = aaqq*SQRT(MAX(ZERO,ONE-aapq1* &
     &                              aapq1))
                                 mxsinj = MAX(mxsinj,Sfmin)
                              ELSE
                                 CALL CCOPY(M,A(1,q),1,Work,1)
                                 CALL CLASCL('G',0,0,aaqq,ONE,M,1,Work, &
     &                              Lda,ierr)
                                 CALL CLASCL('G',0,0,aapp,ONE,M,1,A(1,p)&
     &                              ,Lda,ierr)
                                 CALL CAXPY(M,-CONJG(aapq),Work,1,A(1,p)&
     &                              ,1)
                                 CALL CLASCL('G',0,0,ONE,aapp,M,1,A(1,p)&
     &                              ,Lda,ierr)
                                 Sva(p) = aapp*SQRT(MAX(ZERO,ONE-aapq1* &
     &                              aapq1))
                                 mxsinj = MAX(mxsinj,Sfmin)
                              ENDIF
!           END IF ROTOK THEN ... ELSE
!
!           In the case of cancellation in updating SVA(q), SVA(p)
!           .. recompute SVA(q), SVA(p)
                              IF ( (Sva(q)/aaqq)**2<=rooteps ) THEN
                                 IF ( (aaqq<rootbig) .AND.              &
     &                                (aaqq>rootsfmin) ) THEN
                                    Sva(q) = SCNRM2(M,A(1,q),1)
                                 ELSE
                                    t = ZERO
                                    aaqq = ONE
                                    CALL CLASSQ(M,A(1,q),1,t,aaqq)
                                    Sva(q) = t*SQRT(aaqq)
                                 ENDIF
                              ENDIF
                              IF ( (aapp/aapp0)**2<=rooteps ) THEN
                                 IF ( (aapp<rootbig) .AND.              &
     &                                (aapp>rootsfmin) ) THEN
                                    aapp = SCNRM2(M,A(1,p),1)
                                 ELSE
                                    t = ZERO
                                    aapp = ONE
                                    CALL CLASSQ(M,A(1,p),1,t,aapp)
                                    aapp = t*SQRT(aapp)
                                 ENDIF
                                 Sva(p) = aapp
                              ENDIF
!              end of OK rotation
                           ELSE
                              notrot = notrot + 1
![RTD]      SKIPPED  = SKIPPED  + 1
                              pskipped = pskipped + 1
                              ijblsk = ijblsk + 1
                           ENDIF
                        ELSE
                           notrot = notrot + 1
                           pskipped = pskipped + 1
                           ijblsk = ijblsk + 1
                        ENDIF
!
                        IF ( (i<=swband) .AND. (ijblsk>=blskip) ) THEN
                           Sva(p) = aapp
                           notrot = 0
                           GOTO 20
                        ENDIF
                        IF ( (i<=swband) .AND. (pskipped>rowskip) ) THEN
                           aapp = -aapp
                           notrot = 0
                           EXIT
                        ENDIF
!
                     ENDDO
!        end of the q-loop
!
                     Sva(p) = aapp
!
                  ELSE
!
                     IF ( aapp==ZERO ) notrot = notrot +                &
     &                    MIN(jgl+kbl-1,N) - jgl + 1
                     IF ( aapp<ZERO ) notrot = 0
!
                  ENDIF
!
               ENDDO
!     end of the p-loop
            ENDDO
!     end of the jbc-loop
!2011 bailed out of the jbc-loop
 20         DO p = igl , MIN(igl+kbl-1,N)
               Sva(p) = ABS(Sva(p))
            ENDDO
!**
         ENDDO
!2000 :: end of the ibr-loop
!
!     .. update SVA(N)
         IF ( (Sva(N)<rootbig) .AND. (Sva(N)>rootsfmin) ) THEN
            Sva(N) = SCNRM2(M,A(1,N),1)
         ELSE
            t = ZERO
            aapp = ONE
            CALL CLASSQ(M,A(1,N),1,t,aapp)
            Sva(N) = t*SQRT(aapp)
         ENDIF
!
!     Additional steering devices
!
         IF ( (i<swband) .AND. ((mxaapq<=roottol) .OR. (iswrot<=N)) )   &
     &        swband = i
!
         IF ( (i>swband+1) .AND. (mxaapq<SQRT(REAL(N))*Tol) .AND.       &
     &        (REAL(N)*mxaapq*mxsinj<Tol) ) GOTO 100
!
         IF ( notrot>=emptsw ) GOTO 100
!
      ENDDO
!     end i=1:NSWEEP loop
!
! #:( Reaching this point means that the procedure has not converged.
      Info = Nsweep - 1
      GOTO 200
!
! #:) Reaching this point means numerical convergence after the i-th
!     sweep.
!
 100  Info = 0
! #:) INFO = 0 confirms successful iterations.
!
!     Sort the vector SVA() of column norms.
 200  DO p = 1 , N - 1
         q = ISAMAX(N-p+1,Sva(p),1) + p - 1
         IF ( p/=q ) THEN
            temp1 = Sva(p)
            Sva(p) = Sva(q)
            Sva(q) = temp1
            aapq = D(p)
            D(p) = D(q)
            D(q) = aapq
            CALL CSWAP(M,A(1,p),1,A(1,q),1)
            IF ( rsvec ) CALL CSWAP(mvl,V(1,p),1,V(1,q),1)
         ENDIF
      ENDDO
!
!
!     ..
!     .. END OF CGSVJ1
!     ..
      END SUBROUTINE CGSVJ1
