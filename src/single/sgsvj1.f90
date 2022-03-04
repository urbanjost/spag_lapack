!*==sgsvj1.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SGSVJ1 pre-processor for the routine sgesvj, applies Jacobi rotations targeting only particular pivots.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGSVJ1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgsvj1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgsvj1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgsvj1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV,
!                          EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       REAL               EPS, SFMIN, TOL
!       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP
!       CHARACTER*1        JOBV
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), D( N ), SVA( N ), V( LDV, * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGSVJ1 is called from SGESVJ as a pre-processor and that is its main
!> purpose. It applies Jacobi rotations in the same way as SGESVJ does, but
!> it targets only particular pivots and it does not check convergence
!> (stopping criterion). Few tuning parameters (marked by [TP]) are
!> available for the implementer.
!>
!> Further Details
!> ~~~~~~~~~~~~~~~
!> SGSVJ1 applies few sweeps of Jacobi rotations in the column space of
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
!>          A is REAL array, dimension (LDA,N)
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
!>          D is REAL array, dimension (N)
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
!>          V is REAL array, dimension (LDV,N)
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
!>          WORK is REAL array, dimension (LWORK)
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
!> \date November 2017
!
!> \ingroup realOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)
!
!  =====================================================================
      SUBROUTINE SGSVJ1(Jobv,M,N,N1,A,Lda,D,Sva,Mv,V,Ldv,Eps,Sfmin,Tol, &
     &                  Nsweep,Work,Lwork,Info)
      IMPLICIT NONE
!*--SGSVJ1240
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      REAL Eps , Sfmin , Tol
      INTEGER Info , Lda , Ldv , Lwork , M , Mv , N , N1 , Nsweep
      CHARACTER*1 Jobv
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , D(N) , Sva(N) , V(Ldv,*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Local Parameters ..
      REAL ZERO , HALF , ONE
      PARAMETER (ZERO=0.0E0,HALF=0.5E0,ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      REAL aapp , aapp0 , aapq , aaqq , apoaq , aqoap , big , bigtheta ,&
     &     cs , large , mxaapq , mxsinj , rootbig , rooteps ,           &
     &     rootsfmin , roottol , small , sn , t , temp1 , theta , thsign
      INTEGER blskip , emptsw , i , ibr , igl , ierr , ijblsk , iswrot ,&
     &        jbc , jgl , kbl , mvl , notrot , nblc , nblr , p ,        &
     &        pskipped , q , rowskip , swband
      LOGICAL applv , rotok , rsvec
!     ..
!     .. Local Arrays ..
      REAL fastr(5)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , FLOAT , MIN , SIGN , SQRT
!     ..
!     .. External Functions ..
      REAL SDOT , SNRM2
      INTEGER ISAMAX
      LOGICAL LSAME
      EXTERNAL ISAMAX , LSAME , SDOT , SNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL SAXPY , SCOPY , SLASCL , SLASSQ , SROTM , SSWAP , XERBLA
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
         CALL XERBLA('SGSVJ1',-Info)
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
      large = big/SQRT(FLOAT(M*N))
      bigtheta = ONE/rooteps
      roottol = SQRT(Tol)
!
!     .. Initialize the right singular vector matrix ..
!
!     RSVEC = LSAME( JOBV, 'Y' )
!
      emptsw = N1*(N-N1)
      notrot = 0
      fastr(1) = ZERO
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
!     if SGESVJ is used as a computational routine in the preconditioned
!     Jacobi SVD algorithm SGESVJ.
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
!     .. go go go ...
!
         mxaapq = ZERO
         mxsinj = ZERO
         iswrot = 0
!
         notrot = 0
         pskipped = 0
!
         DO ibr = 1 , nblr
 
            igl = (ibr-1)*kbl + 1
!
!
!........................................................
! ... go to the off diagonal blocks
 
            igl = (ibr-1)*kbl + 1
 
            DO jbc = 1 , nblc
 
               jgl = N1 + (jbc-1)*kbl + 1
 
!        doing the block at ( ibr, jbc )
 
               ijblsk = 0
               DO p = igl , MIN(igl+kbl-1,N1)
 
                  aapp = Sva(p)
 
                  IF ( aapp>ZERO ) THEN
 
                     pskipped = 0
 
                     DO q = jgl , MIN(jgl+kbl-1,N)
!
                        aaqq = Sva(q)
 
                        IF ( aaqq>ZERO ) THEN
                           aapp0 = aapp
!
!     .. M x 2 Jacobi SVD ..
!
!        .. Safe Gram matrix computation ..
!
                           IF ( aaqq>=ONE ) THEN
                              IF ( aapp>=aaqq ) THEN
                                 rotok = (small*aapp)<=aaqq
                              ELSE
                                 rotok = (small*aaqq)<=aapp
                              ENDIF
                              IF ( aapp<(big/aaqq) ) THEN
                                 aapq = (SDOT(M,A(1,p),1,A(1,q),1)*D(p) &
     &                                  *D(q)/aaqq)/aapp
                              ELSE
                                 CALL SCOPY(M,A(1,p),1,Work,1)
                                 CALL SLASCL('G',0,0,aapp,D(p),M,1,Work,&
     &                              Lda,ierr)
                                 aapq = SDOT(M,Work,1,A(1,q),1)*D(q)    &
     &                                  /aaqq
                              ENDIF
                           ELSE
                              IF ( aapp>=aaqq ) THEN
                                 rotok = aapp<=(aaqq/small)
                              ELSE
                                 rotok = aaqq<=(aapp/small)
                              ENDIF
                              IF ( aapp>(small/aaqq) ) THEN
                                 aapq = (SDOT(M,A(1,p),1,A(1,q),1)*D(p) &
     &                                  *D(q)/aaqq)/aapp
                              ELSE
                                 CALL SCOPY(M,A(1,q),1,Work,1)
                                 CALL SLASCL('G',0,0,aaqq,D(q),M,1,Work,&
     &                              Lda,ierr)
                                 aapq = SDOT(M,Work,1,A(1,p),1)*D(p)    &
     &                                  /aapp
                              ENDIF
                           ENDIF
 
                           mxaapq = MAX(mxaapq,ABS(aapq))
 
!        TO rotate or NOT to rotate, THAT is the question ...
!
                           IF ( ABS(aapq)>Tol ) THEN
                              notrot = 0
!           ROTATED  = ROTATED + 1
                              pskipped = 0
                              iswrot = iswrot + 1
!
                              IF ( rotok ) THEN
!
                                 aqoap = aaqq/aapp
                                 apoaq = aapp/aaqq
                                 theta = -HALF*ABS(aqoap-apoaq)/aapq
                                 IF ( aaqq>aapp0 ) theta = -theta
 
                                 IF ( ABS(theta)>bigtheta ) THEN
                                    t = HALF/theta
                                    fastr(3) = t*D(p)/D(q)
                                    fastr(4) = -t*D(q)/D(p)
                                    CALL SROTM(M,A(1,p),1,A(1,q),1,     &
     &                                 fastr)
                                    IF ( rsvec )                        &
     &                                 CALL SROTM(mvl,V(1,p),1,V(1,q),1,&
     &                                 fastr)
                                    Sva(q)                              &
     &                                 = aaqq*SQRT(MAX(ZERO,ONE+t*apoaq*&
     &                                 aapq))
                                    aapp = aapp*SQRT                    &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq))
                                    mxsinj = MAX(mxsinj,ABS(t))
                                 ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                    thsign = -SIGN(ONE,aapq)
                                    IF ( aaqq>aapp0 ) thsign = -thsign
                                    t = ONE/                            &
     &                                  (theta+thsign*SQRT(ONE+theta*   &
     &                                  theta))
                                    cs = SQRT(ONE/(ONE+t*t))
                                    sn = t*cs
                                    mxsinj = MAX(mxsinj,ABS(sn))
                                    Sva(q)                              &
     &                                 = aaqq*SQRT(MAX(ZERO,ONE+t*apoaq*&
     &                                 aapq))
                                    aapp = aapp*SQRT                    &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq))
 
                                    apoaq = D(p)/D(q)
                                    aqoap = D(q)/D(p)
                                    IF ( D(p)>=ONE ) THEN
!
                                       IF ( D(q)>=ONE ) THEN
                                         fastr(3) = t*apoaq
                                         fastr(4) = -t*aqoap
                                         D(p) = D(p)*cs
                                         D(q) = D(q)*cs
                                         CALL SROTM(M,A(1,p),1,A(1,q),1,&
     &                                      fastr)
                                         IF ( rsvec )                   &
     &                                      CALL SROTM(mvl,V(1,p),1,    &
     &                                      V(1,q),1,fastr)
                                       ELSE
                                         CALL SAXPY(M,-t*aqoap,A(1,q),1,&
     &                                      A(1,p),1)
                                         CALL SAXPY(M,cs*sn*apoaq,A(1,p)&
     &                                      ,1,A(1,q),1)
                                         IF ( rsvec ) THEN
                                         CALL SAXPY(mvl,-t*aqoap,V(1,q),&
     &                                      1,V(1,p),1)
                                         CALL SAXPY(mvl,cs*sn*apoaq,    &
     &                                      V(1,p),1,V(1,q),1)
                                         ENDIF
                                         D(p) = D(p)*cs
                                         D(q) = D(q)/cs
                                       ENDIF
                                    ELSEIF ( D(q)>=ONE ) THEN
                                       CALL SAXPY(M,t*apoaq,A(1,p),1,   &
     &                                    A(1,q),1)
                                       CALL SAXPY(M,-cs*sn*aqoap,A(1,q),&
     &                                    1,A(1,p),1)
                                       IF ( rsvec ) THEN
                                         CALL SAXPY(mvl,t*apoaq,V(1,p), &
     &                                      1,V(1,q),1)
                                         CALL SAXPY(mvl,-cs*sn*aqoap,   &
     &                                      V(1,q),1,V(1,p),1)
                                       ENDIF
                                       D(p) = D(p)/cs
                                       D(q) = D(q)*cs
                                    ELSEIF ( D(p)>=D(q) ) THEN
                                       CALL SAXPY(M,-t*aqoap,A(1,q),1,  &
     &                                    A(1,p),1)
                                       CALL SAXPY(M,cs*sn*apoaq,A(1,p), &
     &                                    1,A(1,q),1)
                                       D(p) = D(p)*cs
                                       D(q) = D(q)/cs
                                       IF ( rsvec ) THEN
                                         CALL SAXPY(mvl,-t*aqoap,V(1,q),&
     &                                      1,V(1,p),1)
                                         CALL SAXPY(mvl,cs*sn*apoaq,    &
     &                                      V(1,p),1,V(1,q),1)
                                       ENDIF
                                    ELSE
                                       CALL SAXPY(M,t*apoaq,A(1,p),1,   &
     &                                    A(1,q),1)
                                       CALL SAXPY(M,-cs*sn*aqoap,A(1,q),&
     &                                    1,A(1,p),1)
                                       D(p) = D(p)/cs
                                       D(q) = D(q)*cs
                                       IF ( rsvec ) THEN
                                         CALL SAXPY(mvl,t*apoaq,V(1,p), &
     &                                      1,V(1,q),1)
                                         CALL SAXPY(mvl,-cs*sn*aqoap,   &
     &                                      V(1,q),1,V(1,p),1)
                                       ENDIF
                                    ENDIF
                                 ENDIF
 
                              ELSEIF ( aapp>aaqq ) THEN
                                 CALL SCOPY(M,A(1,p),1,Work,1)
                                 CALL SLASCL('G',0,0,aapp,ONE,M,1,Work, &
     &                              Lda,ierr)
                                 CALL SLASCL('G',0,0,aaqq,ONE,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 temp1 = -aapq*D(p)/D(q)
                                 CALL SAXPY(M,temp1,Work,1,A(1,q),1)
                                 CALL SLASCL('G',0,0,ONE,aaqq,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 Sva(q)                                 &
     &                              = aaqq*SQRT(MAX(ZERO,ONE-aapq*aapq))
                                 mxsinj = MAX(mxsinj,Sfmin)
                              ELSE
                                 CALL SCOPY(M,A(1,q),1,Work,1)
                                 CALL SLASCL('G',0,0,aaqq,ONE,M,1,Work, &
     &                              Lda,ierr)
                                 CALL SLASCL('G',0,0,aapp,ONE,M,1,A(1,p)&
     &                              ,Lda,ierr)
                                 temp1 = -aapq*D(q)/D(p)
                                 CALL SAXPY(M,temp1,Work,1,A(1,p),1)
                                 CALL SLASCL('G',0,0,ONE,aapp,M,1,A(1,p)&
     &                              ,Lda,ierr)
                                 Sva(p)                                 &
     &                              = aapp*SQRT(MAX(ZERO,ONE-aapq*aapq))
                                 mxsinj = MAX(mxsinj,Sfmin)
                              ENDIF
!           END IF ROTOK THEN ... ELSE
!
!           In the case of cancellation in updating SVA(q)
!           .. recompute SVA(q)
                              IF ( (Sva(q)/aaqq)**2<=rooteps ) THEN
                                 IF ( (aaqq<rootbig) .AND.              &
     &                                (aaqq>rootsfmin) ) THEN
                                    Sva(q) = SNRM2(M,A(1,q),1)*D(q)
                                 ELSE
                                    t = ZERO
                                    aaqq = ONE
                                    CALL SLASSQ(M,A(1,q),1,t,aaqq)
                                    Sva(q) = t*SQRT(aaqq)*D(q)
                                 ENDIF
                              ENDIF
                              IF ( (aapp/aapp0)**2<=rooteps ) THEN
                                 IF ( (aapp<rootbig) .AND.              &
     &                                (aapp>rootsfmin) ) THEN
                                    aapp = SNRM2(M,A(1,p),1)*D(p)
                                 ELSE
                                    t = ZERO
                                    aapp = ONE
                                    CALL SLASSQ(M,A(1,p),1,t,aapp)
                                    aapp = t*SQRT(aapp)*D(p)
                                 ENDIF
                                 Sva(p) = aapp
                              ENDIF
!              end of OK rotation
                           ELSE
                              notrot = notrot + 1
!           SKIPPED  = SKIPPED  + 1
                              pskipped = pskipped + 1
                              ijblsk = ijblsk + 1
                           ENDIF
                        ELSE
                           notrot = notrot + 1
                           pskipped = pskipped + 1
                           ijblsk = ijblsk + 1
                        ENDIF
 
!      IF ( NOTROT .GE. EMPTSW )  GO TO 2011
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
 
                     Sva(p) = aapp
!
                  ELSE
                     IF ( aapp==ZERO ) notrot = notrot +                &
     &                    MIN(jgl+kbl-1,N) - jgl + 1
                     IF ( aapp<ZERO ) notrot = 0
!**      IF ( NOTROT .GE. EMPTSW )  GO TO 2011
                  ENDIF
 
               ENDDO
!     end of the p-loop
            ENDDO
!     end of the jbc-loop
!2011 bailed out of the jbc-loop
 20         DO p = igl , MIN(igl+kbl-1,N)
               Sva(p) = ABS(Sva(p))
            ENDDO
!**   IF ( NOTROT .GE. EMPTSW ) GO TO 1994
         ENDDO
!2000 :: end of the ibr-loop
!
!     .. update SVA(N)
         IF ( (Sva(N)<rootbig) .AND. (Sva(N)>rootsfmin) ) THEN
            Sva(N) = SNRM2(M,A(1,N),1)*D(N)
         ELSE
            t = ZERO
            aapp = ONE
            CALL SLASSQ(M,A(1,N),1,t,aapp)
            Sva(N) = t*SQRT(aapp)*D(N)
         ENDIF
!
!     Additional steering devices
!
         IF ( (i<swband) .AND. ((mxaapq<=roottol) .OR. (iswrot<=N)) )   &
     &        swband = i
 
         IF ( (i>swband+1) .AND. (mxaapq<FLOAT(N)*Tol) .AND.            &
     &        (FLOAT(N)*mxaapq*mxsinj<Tol) ) GOTO 100
 
!
         IF ( notrot>=emptsw ) GOTO 100
 
      ENDDO
!     end i=1:NSWEEP loop
! #:) Reaching this point means that the procedure has completed the given
!     number of sweeps.
      Info = Nsweep - 1
      GOTO 200
! #:) Reaching this point means that during the i-th sweep all pivots were
!     below the given threshold, causing early exit.
 
 100  Info = 0
! #:) INFO = 0 confirms successful iterations.
!
!     Sort the vector D
!
 200  DO p = 1 , N - 1
         q = ISAMAX(N-p+1,Sva(p),1) + p - 1
         IF ( p/=q ) THEN
            temp1 = Sva(p)
            Sva(p) = Sva(q)
            Sva(q) = temp1
            temp1 = D(p)
            D(p) = D(q)
            D(q) = temp1
            CALL SSWAP(M,A(1,p),1,A(1,q),1)
            IF ( rsvec ) CALL SSWAP(mvl,V(1,p),1,V(1,q),1)
         ENDIF
      ENDDO
!
!     ..
!     .. END OF SGSVJ1
!     ..
      END SUBROUTINE SGSVJ1
