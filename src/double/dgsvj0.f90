!*==dgsvj0.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DGSVJ0 pre-processor for the routine dgesvj.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGSVJ0 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgsvj0.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgsvj0.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgsvj0.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS,
!                          SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, NSWEEP
!       DOUBLE PRECISION   EPS, SFMIN, TOL
!       CHARACTER*1        JOBV
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), SVA( N ), D( N ), V( LDV, * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGSVJ0 is called from DGESVJ as a pre-processor and that is its main
!> purpose. It applies Jacobi rotations in the same way as DGESVJ does, but
!> it does not check convergence (stopping criterion). Few tuning
!> parameters (marked by [TP]) are available for the implementer.
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
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, M-by-N matrix A, such that A*diag(D) represents
!>          the input matrix.
!>          On exit,
!>          A_onexit * D_onexit represents the input matrix A*diag(D)
!>          post-multiplied by a sequence of Jacobi rotations, where the
!>          rotation threshold and the total number of sweeps are given in
!>          TOL and NSWEEP, respectively.
!>          (See the descriptions of D, TOL and NSWEEP.)
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
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The array D accumulates the scaling factors from the fast scaled
!>          Jacobi rotations.
!>          On entry, A*diag(D) represents the input matrix.
!>          On exit, A_onexit*diag(D_onexit) represents the input matrix
!>          post-multiplied by a sequence of Jacobi rotations, where the
!>          rotation threshold and the total number of sweeps are given in
!>          TOL and NSWEEP, respectively.
!>          (See the descriptions of A, TOL and NSWEEP.)
!> \endverbatim
!>
!> \param[in,out] SVA
!> \verbatim
!>          SVA is DOUBLE PRECISION array, dimension (N)
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
!>          V is DOUBLE PRECISION array, dimension (LDV,N)
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
!>          EPS is DOUBLE PRECISION
!>          EPS = DLAMCH('Epsilon')
!> \endverbatim
!>
!> \param[in] SFMIN
!> \verbatim
!>          SFMIN is DOUBLE PRECISION
!>          SFMIN = DLAMCH('Safe Minimum')
!> \endverbatim
!>
!> \param[in] TOL
!> \verbatim
!>          TOL is DOUBLE PRECISION
!>          TOL is the threshold for Jacobi rotations. For a pair
!>          A(:,p), A(:,q) of pivot columns, the Jacobi rotation is
!>          applied only if DABS(COS(angle(A(:,p),A(:,q)))) > TOL.
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
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
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
!> \ingroup doubleOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> DGSVJ0 is used just to enable DGESVJ to call a simplified version of
!> itself to work on a submatrix of the original matrix.
!>
!> \par Contributors:
!  ==================
!>
!> Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)
!>
!> \par Bugs, Examples and Comments:
!  =================================
!>
!> Please report all bugs and send interesting test examples and comments to
!> drmac@math.hr. Thank you.
!
!  =====================================================================
      SUBROUTINE DGSVJ0(Jobv,M,N,A,Lda,D,Sva,Mv,V,Ldv,Eps,Sfmin,Tol,    &
     &                  Nsweep,Work,Lwork,Info)
      IMPLICIT NONE
!*--DGSVJ0222
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldv , Lwork , M , Mv , N , Nsweep
      DOUBLE PRECISION Eps , Sfmin , Tol
      CHARACTER*1 Jobv
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , Sva(N) , D(N) , V(Ldv,*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Local Parameters ..
      DOUBLE PRECISION ZERO , HALF , ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION aapp , aapp0 , aapq , aaqq , apoaq , aqoap ,     &
     &                 big , bigtheta , cs , mxaapq , mxsinj , rootbig ,&
     &                 rooteps , rootsfmin , roottol , small , sn , t , &
     &                 temp1 , theta , thsign
      INTEGER blskip , emptsw , i , ibr , ierr , igl , ijblsk , ir1 ,   &
     &        iswrot , jbc , jgl , kbl , lkahead , mvl , nbl , notrot , &
     &        p , pskipped , q , rowskip , swband
      LOGICAL applv , rotok , rsvec
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION fastr(5)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DABS , MAX , DBLE , MIN , DSIGN , DSQRT
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DDOT , DNRM2
      INTEGER IDAMAX
      LOGICAL LSAME
      EXTERNAL IDAMAX , LSAME , DDOT , DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL DAXPY , DCOPY , DLASCL , DLASSQ , DROTM , DSWAP , XERBLA
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
      ELSEIF ( Lda<M ) THEN
         Info = -5
      ELSEIF ( (rsvec .OR. applv) .AND. (Mv<0) ) THEN
         Info = -8
      ELSEIF ( (rsvec .AND. (Ldv<N)) .OR. (applv .AND. (Ldv<Mv)) ) THEN
         Info = -10
      ELSEIF ( Tol<=Eps ) THEN
         Info = -13
      ELSEIF ( Nsweep<0 ) THEN
         Info = -14
      ELSEIF ( Lwork<M ) THEN
         Info = -16
      ELSE
         Info = 0
      ENDIF
!
!     #:(
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGSVJ0',-Info)
         RETURN
      ENDIF
!
      IF ( rsvec ) THEN
         mvl = N
      ELSEIF ( applv ) THEN
         mvl = Mv
      ENDIF
      rsvec = rsvec .OR. applv
 
      rooteps = DSQRT(Eps)
      rootsfmin = DSQRT(Sfmin)
      small = Sfmin/Eps
      big = ONE/Sfmin
      rootbig = ONE/rootsfmin
      bigtheta = ONE/rooteps
      roottol = DSQRT(Tol)
!
!     -#- Row-cyclic Jacobi SVD algorithm with column pivoting -#-
!
      emptsw = (N*(N-1))/2
      notrot = 0
      fastr(1) = ZERO
!
!     -#- Row-cyclic pivot strategy with de Rijk's pivoting -#-
!
 
      swband = 0
![TP] SWBAND is a tuning parameter. It is meaningful and effective
!     if SGESVJ is used as a computational routine in the preconditioned
!     Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure
!     ......
 
      kbl = MIN(8,N)
![TP] KBL is a tuning parameter that defines the tile size in the
!     tiling of the p-q loops of pivot pairs. In general, an optimal
!     value of KBL depends on the matrix dimensions and on the
!     parameters of the computer's memory.
!
      nbl = N/kbl
      IF ( (nbl*kbl)/=N ) nbl = nbl + 1
 
      blskip = (kbl**2) + 1
![TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.
 
      rowskip = MIN(5,kbl)
![TP] ROWSKIP is a tuning parameter.
 
      lkahead = 1
![TP] LKAHEAD is a tuning parameter.
      swband = 0
      pskipped = 0
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
         DO ibr = 1 , nbl
 
            igl = (ibr-1)*kbl + 1
!
            DO ir1 = 0 , MIN(lkahead,nbl-ibr)
!
               igl = igl + ir1*kbl
!
               DO p = igl , MIN(igl+kbl-1,N-1)
 
!     .. de Rijk's pivoting
                  q = IDAMAX(N-p+1,Sva(p),1) + p - 1
                  IF ( p/=q ) THEN
                     CALL DSWAP(M,A(1,p),1,A(1,q),1)
                     IF ( rsvec ) CALL DSWAP(mvl,V(1,p),1,V(1,q),1)
                     temp1 = Sva(p)
                     Sva(p) = Sva(q)
                     Sva(q) = temp1
                     temp1 = D(p)
                     D(p) = D(q)
                     D(q) = temp1
                  ENDIF
!
                  IF ( ir1==0 ) THEN
!
!        Column norms are periodically updated by explicit
!        norm computation.
!        Caveat:
!        Some BLAS implementations compute DNRM2(M,A(1,p),1)
!        as DSQRT(DDOT(M,A(1,p),1,A(1,p),1)), which may result in
!        overflow for ||A(:,p)||_2 > DSQRT(overflow_threshold), and
!        underflow for ||A(:,p)||_2 < DSQRT(underflow_threshold).
!        Hence, DNRM2 cannot be trusted, not even in the case when
!        the true norm is far from the under(over)flow boundaries.
!        If properly implemented DNRM2 is available, the IF-THEN-ELSE
!        below should read "AAPP = DNRM2( M, A(1,p), 1 ) * D(p)".
!
                     IF ( (Sva(p)<rootbig) .AND. (Sva(p)>rootsfmin) )   &
     &                    THEN
                        Sva(p) = DNRM2(M,A(1,p),1)*D(p)
                     ELSE
                        temp1 = ZERO
                        aapp = ONE
                        CALL DLASSQ(M,A(1,p),1,temp1,aapp)
                        Sva(p) = temp1*DSQRT(aapp)*D(p)
                     ENDIF
                     aapp = Sva(p)
                  ELSE
                     aapp = Sva(p)
                  ENDIF
 
!
                  IF ( aapp>ZERO ) THEN
!
                     pskipped = 0
!
                     DO q = p + 1 , MIN(igl+kbl-1,N)
!
                        aaqq = Sva(q)
 
                        IF ( aaqq>ZERO ) THEN
!
                           aapp0 = aapp
                           IF ( aaqq>=ONE ) THEN
                              rotok = (small*aapp)<=aaqq
                              IF ( aapp<(big/aaqq) ) THEN
                                 aapq = (DDOT(M,A(1,p),1,A(1,q),1)*D(p) &
     &                                  *D(q)/aaqq)/aapp
                              ELSE
                                 CALL DCOPY(M,A(1,p),1,Work,1)
                                 CALL DLASCL('G',0,0,aapp,D(p),M,1,Work,&
     &                              Lda,ierr)
                                 aapq = DDOT(M,Work,1,A(1,q),1)*D(q)    &
     &                                  /aaqq
                              ENDIF
                           ELSE
                              rotok = aapp<=(aaqq/small)
                              IF ( aapp>(small/aaqq) ) THEN
                                 aapq = (DDOT(M,A(1,p),1,A(1,q),1)*D(p) &
     &                                  *D(q)/aaqq)/aapp
                              ELSE
                                 CALL DCOPY(M,A(1,q),1,Work,1)
                                 CALL DLASCL('G',0,0,aaqq,D(q),M,1,Work,&
     &                              Lda,ierr)
                                 aapq = DDOT(M,Work,1,A(1,p),1)*D(p)    &
     &                                  /aapp
                              ENDIF
                           ENDIF
!
                           mxaapq = MAX(mxaapq,DABS(aapq))
!
!        TO rotate or NOT to rotate, THAT is the question ...
!
                           IF ( DABS(aapq)>Tol ) THEN
!
!           .. rotate
!           ROTATED = ROTATED + ONE
!
                              IF ( ir1==0 ) THEN
                                 notrot = 0
                                 pskipped = 0
                                 iswrot = iswrot + 1
                              ENDIF
!
                              IF ( rotok ) THEN
!
                                 aqoap = aaqq/aapp
                                 apoaq = aapp/aaqq
                                 theta = -HALF*DABS(aqoap-apoaq)/aapq
!
                                 IF ( DABS(theta)>bigtheta ) THEN
!
                                    t = HALF/theta
                                    fastr(3) = t*D(p)/D(q)
                                    fastr(4) = -t*D(q)/D(p)
                                    CALL DROTM(M,A(1,p),1,A(1,q),1,     &
     &                                 fastr)
                                    IF ( rsvec )                        &
     &                                 CALL DROTM(mvl,V(1,p),1,V(1,q),1,&
     &                                 fastr)
                                    Sva(q) = aaqq*DSQRT(MAX(ZERO,ONE+t* &
     &                                 apoaq*aapq))
                                    aapp = aapp*DSQRT                   &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq))
                                    mxsinj = MAX(mxsinj,DABS(t))
!
                                 ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                    thsign = -DSIGN(ONE,aapq)
                                    t = ONE/                            &
     &                                  (theta+thsign*DSQRT(ONE+theta*  &
     &                                  theta))
                                    cs = DSQRT(ONE/(ONE+t*t))
                                    sn = t*cs
!
                                    mxsinj = MAX(mxsinj,DABS(sn))
                                    Sva(q) = aaqq*DSQRT(MAX(ZERO,ONE+t* &
     &                                 apoaq*aapq))
                                    aapp = aapp*DSQRT                   &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq))
!
                                    apoaq = D(p)/D(q)
                                    aqoap = D(q)/D(p)
                                    IF ( D(p)>=ONE ) THEN
                                       IF ( D(q)>=ONE ) THEN
                                         fastr(3) = t*apoaq
                                         fastr(4) = -t*aqoap
                                         D(p) = D(p)*cs
                                         D(q) = D(q)*cs
                                         CALL DROTM(M,A(1,p),1,A(1,q),1,&
     &                                      fastr)
                                         IF ( rsvec )                   &
     &                                      CALL DROTM(mvl,V(1,p),1,    &
     &                                      V(1,q),1,fastr)
                                       ELSE
                                         CALL DAXPY(M,-t*aqoap,A(1,q),1,&
     &                                      A(1,p),1)
                                         CALL DAXPY(M,cs*sn*apoaq,A(1,p)&
     &                                      ,1,A(1,q),1)
                                         D(p) = D(p)*cs
                                         D(q) = D(q)/cs
                                         IF ( rsvec ) THEN
                                         CALL DAXPY(mvl,-t*aqoap,V(1,q),&
     &                                      1,V(1,p),1)
                                         CALL DAXPY(mvl,cs*sn*apoaq,    &
     &                                      V(1,p),1,V(1,q),1)
                                         ENDIF
                                       ENDIF
                                    ELSEIF ( D(q)>=ONE ) THEN
                                       CALL DAXPY(M,t*apoaq,A(1,p),1,   &
     &                                    A(1,q),1)
                                       CALL DAXPY(M,-cs*sn*aqoap,A(1,q),&
     &                                    1,A(1,p),1)
                                       D(p) = D(p)/cs
                                       D(q) = D(q)*cs
                                       IF ( rsvec ) THEN
                                         CALL DAXPY(mvl,t*apoaq,V(1,p), &
     &                                      1,V(1,q),1)
                                         CALL DAXPY(mvl,-cs*sn*aqoap,   &
     &                                      V(1,q),1,V(1,p),1)
                                       ENDIF
                                    ELSEIF ( D(p)>=D(q) ) THEN
                                       CALL DAXPY(M,-t*aqoap,A(1,q),1,  &
     &                                    A(1,p),1)
                                       CALL DAXPY(M,cs*sn*apoaq,A(1,p), &
     &                                    1,A(1,q),1)
                                       D(p) = D(p)*cs
                                       D(q) = D(q)/cs
                                       IF ( rsvec ) THEN
                                         CALL DAXPY(mvl,-t*aqoap,V(1,q),&
     &                                      1,V(1,p),1)
                                         CALL DAXPY(mvl,cs*sn*apoaq,    &
     &                                      V(1,p),1,V(1,q),1)
                                       ENDIF
                                    ELSE
                                       CALL DAXPY(M,t*apoaq,A(1,p),1,   &
     &                                    A(1,q),1)
                                       CALL DAXPY(M,-cs*sn*aqoap,A(1,q),&
     &                                    1,A(1,p),1)
                                       D(p) = D(p)/cs
                                       D(q) = D(q)*cs
                                       IF ( rsvec ) THEN
                                         CALL DAXPY(mvl,t*apoaq,V(1,p), &
     &                                      1,V(1,q),1)
                                         CALL DAXPY(mvl,-cs*sn*aqoap,   &
     &                                      V(1,q),1,V(1,p),1)
                                       ENDIF
                                    ENDIF
                                 ENDIF
!
                              ELSE
!              .. have to use modified Gram-Schmidt like transformation
                                 CALL DCOPY(M,A(1,p),1,Work,1)
                                 CALL DLASCL('G',0,0,aapp,ONE,M,1,Work, &
     &                              Lda,ierr)
                                 CALL DLASCL('G',0,0,aaqq,ONE,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 temp1 = -aapq*D(p)/D(q)
                                 CALL DAXPY(M,temp1,Work,1,A(1,q),1)
                                 CALL DLASCL('G',0,0,ONE,aaqq,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 Sva(q)                                 &
     &                              = aaqq*DSQRT(MAX(ZERO,ONE-aapq*aapq)&
     &                              )
                                 mxsinj = MAX(mxsinj,Sfmin)
                              ENDIF
!           END IF ROTOK THEN ... ELSE
!
!           In the case of cancellation in updating SVA(q), SVA(p)
!           recompute SVA(q), SVA(p).
                              IF ( (Sva(q)/aaqq)**2<=rooteps ) THEN
                                 IF ( (aaqq<rootbig) .AND.              &
     &                                (aaqq>rootsfmin) ) THEN
                                    Sva(q) = DNRM2(M,A(1,q),1)*D(q)
                                 ELSE
                                    t = ZERO
                                    aaqq = ONE
                                    CALL DLASSQ(M,A(1,q),1,t,aaqq)
                                    Sva(q) = t*DSQRT(aaqq)*D(q)
                                 ENDIF
                              ENDIF
                              IF ( (aapp/aapp0)<=rooteps ) THEN
                                 IF ( (aapp<rootbig) .AND.              &
     &                                (aapp>rootsfmin) ) THEN
                                    aapp = DNRM2(M,A(1,p),1)*D(p)
                                 ELSE
                                    t = ZERO
                                    aapp = ONE
                                    CALL DLASSQ(M,A(1,p),1,t,aapp)
                                    aapp = t*DSQRT(aapp)*D(p)
                                 ENDIF
                                 Sva(p) = aapp
                              ENDIF
!
                           ELSE
!        A(:,p) and A(:,q) already numerically orthogonal
                              IF ( ir1==0 ) notrot = notrot + 1
                              pskipped = pskipped + 1
                           ENDIF
                        ELSE
!        A(:,q) is zero column
                           IF ( ir1==0 ) notrot = notrot + 1
                           pskipped = pskipped + 1
                        ENDIF
!
                        IF ( (i<=swband) .AND. (pskipped>rowskip) ) THEN
                           IF ( ir1==0 ) aapp = -aapp
                           notrot = 0
                           EXIT
                        ENDIF
!
                     ENDDO
!     END q-LOOP
!
!     bailed out of q-loop
 
                     Sva(p) = aapp
 
                  ELSE
                     Sva(p) = aapp
                     IF ( (ir1==0) .AND. (aapp==ZERO) )                 &
     &                    notrot = notrot + MIN(igl+kbl-1,N) - p
                  ENDIF
!
               ENDDO
!     end of the p-loop
!     end of doing the block ( ibr, ibr )
            ENDDO
!     end of ir1-loop
!
!........................................................
! ... go to the off diagonal blocks
!
            igl = (ibr-1)*kbl + 1
!
            DO jbc = ibr + 1 , nbl
!
               jgl = (jbc-1)*kbl + 1
!
!        doing the block at ( ibr, jbc )
!
               ijblsk = 0
               DO p = igl , MIN(igl+kbl-1,N)
!
                  aapp = Sva(p)
!
                  IF ( aapp>ZERO ) THEN
!
                     pskipped = 0
!
                     DO q = jgl , MIN(jgl+kbl-1,N)
!
                        aaqq = Sva(q)
!
                        IF ( aaqq>ZERO ) THEN
                           aapp0 = aapp
!
!     -#- M x 2 Jacobi SVD -#-
!
!        -#- Safe Gram matrix computation -#-
!
                           IF ( aaqq>=ONE ) THEN
                              IF ( aapp>=aaqq ) THEN
                                 rotok = (small*aapp)<=aaqq
                              ELSE
                                 rotok = (small*aaqq)<=aapp
                              ENDIF
                              IF ( aapp<(big/aaqq) ) THEN
                                 aapq = (DDOT(M,A(1,p),1,A(1,q),1)*D(p) &
     &                                  *D(q)/aaqq)/aapp
                              ELSE
                                 CALL DCOPY(M,A(1,p),1,Work,1)
                                 CALL DLASCL('G',0,0,aapp,D(p),M,1,Work,&
     &                              Lda,ierr)
                                 aapq = DDOT(M,Work,1,A(1,q),1)*D(q)    &
     &                                  /aaqq
                              ENDIF
                           ELSE
                              IF ( aapp>=aaqq ) THEN
                                 rotok = aapp<=(aaqq/small)
                              ELSE
                                 rotok = aaqq<=(aapp/small)
                              ENDIF
                              IF ( aapp>(small/aaqq) ) THEN
                                 aapq = (DDOT(M,A(1,p),1,A(1,q),1)*D(p) &
     &                                  *D(q)/aaqq)/aapp
                              ELSE
                                 CALL DCOPY(M,A(1,q),1,Work,1)
                                 CALL DLASCL('G',0,0,aaqq,D(q),M,1,Work,&
     &                              Lda,ierr)
                                 aapq = DDOT(M,Work,1,A(1,p),1)*D(p)    &
     &                                  /aapp
                              ENDIF
                           ENDIF
!
                           mxaapq = MAX(mxaapq,DABS(aapq))
!
!        TO rotate or NOT to rotate, THAT is the question ...
!
                           IF ( DABS(aapq)>Tol ) THEN
                              notrot = 0
!           ROTATED  = ROTATED + 1
                              pskipped = 0
                              iswrot = iswrot + 1
!
                              IF ( rotok ) THEN
!
                                 aqoap = aaqq/aapp
                                 apoaq = aapp/aaqq
                                 theta = -HALF*DABS(aqoap-apoaq)/aapq
                                 IF ( aaqq>aapp0 ) theta = -theta
!
                                 IF ( DABS(theta)>bigtheta ) THEN
                                    t = HALF/theta
                                    fastr(3) = t*D(p)/D(q)
                                    fastr(4) = -t*D(q)/D(p)
                                    CALL DROTM(M,A(1,p),1,A(1,q),1,     &
     &                                 fastr)
                                    IF ( rsvec )                        &
     &                                 CALL DROTM(mvl,V(1,p),1,V(1,q),1,&
     &                                 fastr)
                                    Sva(q) = aaqq*DSQRT(MAX(ZERO,ONE+t* &
     &                                 apoaq*aapq))
                                    aapp = aapp*DSQRT                   &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq))
                                    mxsinj = MAX(mxsinj,DABS(t))
                                 ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                    thsign = -DSIGN(ONE,aapq)
                                    IF ( aaqq>aapp0 ) thsign = -thsign
                                    t = ONE/                            &
     &                                  (theta+thsign*DSQRT(ONE+theta*  &
     &                                  theta))
                                    cs = DSQRT(ONE/(ONE+t*t))
                                    sn = t*cs
                                    mxsinj = MAX(mxsinj,DABS(sn))
                                    Sva(q) = aaqq*DSQRT(MAX(ZERO,ONE+t* &
     &                                 apoaq*aapq))
                                    aapp = aapp*DSQRT                   &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq))
!
                                    apoaq = D(p)/D(q)
                                    aqoap = D(q)/D(p)
                                    IF ( D(p)>=ONE ) THEN
!
                                       IF ( D(q)>=ONE ) THEN
                                         fastr(3) = t*apoaq
                                         fastr(4) = -t*aqoap
                                         D(p) = D(p)*cs
                                         D(q) = D(q)*cs
                                         CALL DROTM(M,A(1,p),1,A(1,q),1,&
     &                                      fastr)
                                         IF ( rsvec )                   &
     &                                      CALL DROTM(mvl,V(1,p),1,    &
     &                                      V(1,q),1,fastr)
                                       ELSE
                                         CALL DAXPY(M,-t*aqoap,A(1,q),1,&
     &                                      A(1,p),1)
                                         CALL DAXPY(M,cs*sn*apoaq,A(1,p)&
     &                                      ,1,A(1,q),1)
                                         IF ( rsvec ) THEN
                                         CALL DAXPY(mvl,-t*aqoap,V(1,q),&
     &                                      1,V(1,p),1)
                                         CALL DAXPY(mvl,cs*sn*apoaq,    &
     &                                      V(1,p),1,V(1,q),1)
                                         ENDIF
                                         D(p) = D(p)*cs
                                         D(q) = D(q)/cs
                                       ENDIF
                                    ELSEIF ( D(q)>=ONE ) THEN
                                       CALL DAXPY(M,t*apoaq,A(1,p),1,   &
     &                                    A(1,q),1)
                                       CALL DAXPY(M,-cs*sn*aqoap,A(1,q),&
     &                                    1,A(1,p),1)
                                       IF ( rsvec ) THEN
                                         CALL DAXPY(mvl,t*apoaq,V(1,p), &
     &                                      1,V(1,q),1)
                                         CALL DAXPY(mvl,-cs*sn*aqoap,   &
     &                                      V(1,q),1,V(1,p),1)
                                       ENDIF
                                       D(p) = D(p)/cs
                                       D(q) = D(q)*cs
                                    ELSEIF ( D(p)>=D(q) ) THEN
                                       CALL DAXPY(M,-t*aqoap,A(1,q),1,  &
     &                                    A(1,p),1)
                                       CALL DAXPY(M,cs*sn*apoaq,A(1,p), &
     &                                    1,A(1,q),1)
                                       D(p) = D(p)*cs
                                       D(q) = D(q)/cs
                                       IF ( rsvec ) THEN
                                         CALL DAXPY(mvl,-t*aqoap,V(1,q),&
     &                                      1,V(1,p),1)
                                         CALL DAXPY(mvl,cs*sn*apoaq,    &
     &                                      V(1,p),1,V(1,q),1)
                                       ENDIF
                                    ELSE
                                       CALL DAXPY(M,t*apoaq,A(1,p),1,   &
     &                                    A(1,q),1)
                                       CALL DAXPY(M,-cs*sn*aqoap,A(1,q),&
     &                                    1,A(1,p),1)
                                       D(p) = D(p)/cs
                                       D(q) = D(q)*cs
                                       IF ( rsvec ) THEN
                                         CALL DAXPY(mvl,t*apoaq,V(1,p), &
     &                                      1,V(1,q),1)
                                         CALL DAXPY(mvl,-cs*sn*aqoap,   &
     &                                      V(1,q),1,V(1,p),1)
                                       ENDIF
                                    ENDIF
                                 ENDIF
!
                              ELSEIF ( aapp>aaqq ) THEN
                                 CALL DCOPY(M,A(1,p),1,Work,1)
                                 CALL DLASCL('G',0,0,aapp,ONE,M,1,Work, &
     &                              Lda,ierr)
                                 CALL DLASCL('G',0,0,aaqq,ONE,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 temp1 = -aapq*D(p)/D(q)
                                 CALL DAXPY(M,temp1,Work,1,A(1,q),1)
                                 CALL DLASCL('G',0,0,ONE,aaqq,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 Sva(q)                                 &
     &                              = aaqq*DSQRT(MAX(ZERO,ONE-aapq*aapq)&
     &                              )
                                 mxsinj = MAX(mxsinj,Sfmin)
                              ELSE
                                 CALL DCOPY(M,A(1,q),1,Work,1)
                                 CALL DLASCL('G',0,0,aaqq,ONE,M,1,Work, &
     &                              Lda,ierr)
                                 CALL DLASCL('G',0,0,aapp,ONE,M,1,A(1,p)&
     &                              ,Lda,ierr)
                                 temp1 = -aapq*D(q)/D(p)
                                 CALL DAXPY(M,temp1,Work,1,A(1,p),1)
                                 CALL DLASCL('G',0,0,ONE,aapp,M,1,A(1,p)&
     &                              ,Lda,ierr)
                                 Sva(p)                                 &
     &                              = aapp*DSQRT(MAX(ZERO,ONE-aapq*aapq)&
     &                              )
                                 mxsinj = MAX(mxsinj,Sfmin)
                              ENDIF
!           END IF ROTOK THEN ... ELSE
!
!           In the case of cancellation in updating SVA(q)
!           .. recompute SVA(q)
                              IF ( (Sva(q)/aaqq)**2<=rooteps ) THEN
                                 IF ( (aaqq<rootbig) .AND.              &
     &                                (aaqq>rootsfmin) ) THEN
                                    Sva(q) = DNRM2(M,A(1,q),1)*D(q)
                                 ELSE
                                    t = ZERO
                                    aaqq = ONE
                                    CALL DLASSQ(M,A(1,q),1,t,aaqq)
                                    Sva(q) = t*DSQRT(aaqq)*D(q)
                                 ENDIF
                              ENDIF
                              IF ( (aapp/aapp0)**2<=rooteps ) THEN
                                 IF ( (aapp<rootbig) .AND.              &
     &                                (aapp>rootsfmin) ) THEN
                                    aapp = DNRM2(M,A(1,p),1)*D(p)
                                 ELSE
                                    t = ZERO
                                    aapp = ONE
                                    CALL DLASSQ(M,A(1,p),1,t,aapp)
                                    aapp = t*DSQRT(aapp)*D(p)
                                 ENDIF
                                 Sva(p) = aapp
                              ENDIF
!              end of OK rotation
                           ELSE
                              notrot = notrot + 1
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
                     IF ( aapp==ZERO ) notrot = notrot +                &
     &                    MIN(jgl+kbl-1,N) - jgl + 1
                     IF ( aapp<ZERO ) notrot = 0
                  ENDIF
 
               ENDDO
!     end of the p-loop
            ENDDO
!     end of the jbc-loop
!2011 bailed out of the jbc-loop
 20         DO p = igl , MIN(igl+kbl-1,N)
               Sva(p) = DABS(Sva(p))
            ENDDO
!
         ENDDO
!2000 :: end of the ibr-loop
!
!     .. update SVA(N)
         IF ( (Sva(N)<rootbig) .AND. (Sva(N)>rootsfmin) ) THEN
            Sva(N) = DNRM2(M,A(1,N),1)*D(N)
         ELSE
            t = ZERO
            aapp = ONE
            CALL DLASSQ(M,A(1,N),1,t,aapp)
            Sva(N) = t*DSQRT(aapp)*D(N)
         ENDIF
!
!     Additional steering devices
!
         IF ( (i<swband) .AND. ((mxaapq<=roottol) .OR. (iswrot<=N)) )   &
     &        swband = i
!
         IF ( (i>swband+1) .AND. (mxaapq<DBLE(N)*Tol) .AND.             &
     &        (DBLE(N)*mxaapq*mxsinj<Tol) ) GOTO 100
!
         IF ( notrot>=emptsw ) GOTO 100
 
      ENDDO
!     end i=1:NSWEEP loop
! #:) Reaching this point means that the procedure has completed the given
!     number of iterations.
      Info = Nsweep - 1
      GOTO 200
! #:) Reaching this point means that during the i-th sweep all pivots were
!     below the given tolerance, causing early exit.
!
 100  Info = 0
! #:) INFO = 0 confirms successful iterations.
!
!     Sort the vector D.
 200  DO p = 1 , N - 1
         q = IDAMAX(N-p+1,Sva(p),1) + p - 1
         IF ( p/=q ) THEN
            temp1 = Sva(p)
            Sva(p) = Sva(q)
            Sva(q) = temp1
            temp1 = D(p)
            D(p) = D(q)
            D(q) = temp1
            CALL DSWAP(M,A(1,p),1,A(1,q),1)
            IF ( rsvec ) CALL DSWAP(mvl,V(1,p),1,V(1,q),1)
         ENDIF
      ENDDO
!
!     ..
!     .. END OF DGSVJ0
!     ..
      END SUBROUTINE DGSVJ0
