!*==ddrvls.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b ddrvls
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ddrvls( DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB,
!                          NBVAL, NXVAL, THRESH, TSTERR, A, COPYA, B,
!                          COPYB, C, S, COPYS, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NM, NN, NNB, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            MVAL( * ), NBVAL( * ), NSVAL( * ),
!      $                   NVAL( * ), NXVAL( * )
!       DOUBLE PRECISION   A( * ), B( * ), C( * ), COPYA( * ), COPYB( * ),
!      $                   COPYS( * ), S( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ddrvls tests the least squares driver routines DGELS, DGETSLS, DGELSS, DGELSY,
!> and DGELSD.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          The matrix types to be used for testing.  Matrices of type j
!>          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
!>          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
!>          The matrix of type j is generated as follows:
!>          j=1: A = U*D*V where U and V are random orthogonal matrices
!>               and D has random entries (> 0.1) taken from a uniform
!>               distribution (0,1). A is full rank.
!>          j=2: The same of 1, but A is scaled up.
!>          j=3: The same of 1, but A is scaled down.
!>          j=4: A = U*D*V where U and V are random orthogonal matrices
!>               and D has 3*min(M,N)/4 random entries (> 0.1) taken
!>               from a uniform distribution (0,1) and the remaining
!>               entries set to 0. A is rank-deficient.
!>          j=5: The same of 4, but A is scaled up.
!>          j=6: The same of 5, but A is scaled down.
!> \endverbatim
!>
!> \param[in] NM
!> \verbatim
!>          NM is INTEGER
!>          The number of values of M contained in the vector MVAL.
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NM)
!>          The values of the matrix row dimension M.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER
!>          The number of values of N contained in the vector NVAL.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NN)
!>          The values of the matrix column dimension N.
!> \endverbatim
!>
!> \param[in] NNS
!> \verbatim
!>          NNS is INTEGER
!>          The number of values of NRHS contained in the vector NSVAL.
!> \endverbatim
!>
!> \param[in] NSVAL
!> \verbatim
!>          NSVAL is INTEGER array, dimension (NNS)
!>          The values of the number of right hand sides NRHS.
!> \endverbatim
!>
!> \param[in] NNB
!> \verbatim
!>          NNB is INTEGER
!>          The number of values of NB and NX contained in the
!>          vectors NBVAL and NXVAL.  The blocking parameters are used
!>          in pairs (NB,NX).
!> \endverbatim
!>
!> \param[in] NBVAL
!> \verbatim
!>          NBVAL is INTEGER array, dimension (NNB)
!>          The values of the blocksize NB.
!> \endverbatim
!>
!> \param[in] NXVAL
!> \verbatim
!>          NXVAL is INTEGER array, dimension (NNB)
!>          The values of the crossover point NX.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  To have
!>          every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[in] TSTERR
!> \verbatim
!>          TSTERR is LOGICAL
!>          Flag that indicates whether error exits are to be tested.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (MMAX*NMAX)
!>          where MMAX is the maximum value of M in MVAL and NMAX is the
!>          maximum value of N in NVAL.
!> \endverbatim
!>
!> \param[out] COPYA
!> \verbatim
!>          COPYA is DOUBLE PRECISION array, dimension (MMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (MMAX*NSMAX)
!>          where MMAX is the maximum value of M in MVAL and NSMAX is the
!>          maximum value of NRHS in NSVAL.
!> \endverbatim
!>
!> \param[out] COPYB
!> \verbatim
!>          COPYB is DOUBLE PRECISION array, dimension (MMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (MMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension
!>                      (min(MMAX,NMAX))
!> \endverbatim
!>
!> \param[out] COPYS
!> \verbatim
!>          COPYS is DOUBLE PRECISION array, dimension
!>                      (min(MMAX,NMAX))
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number for output.
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE ddrvls(Dotype,Nm,Mval,Nn,Nval,Nns,Nsval,Nnb,Nbval,     &
     &                  Nxval,Thresh,Tsterr,A,Copya,B,Copyb,C,S,Copys,  &
     &                  Nout)
      IMPLICIT NONE
!*--DDRVLS196
!
!  -- LAPACK test routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nm , Nn , Nnb , Nns , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Mval(*) , Nbval(*) , Nsval(*) , Nval(*) , Nxval(*)
      DOUBLE PRECISION A(*) , B(*) , C(*) , Copya(*) , Copyb(*) ,       &
     &                 Copys(*) , S(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NTESTS
      PARAMETER (NTESTS=16)
      INTEGER SMLSIZ
      PARAMETER (SMLSIZ=25)
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      CHARACTER trans
      CHARACTER*3 path
      INTEGER crank , i , im , imb , in , inb , info , ins , irank ,    &
     &        iscale , itran , itype , j , k , lda , ldb , ldwork ,     &
     &        lwlsy , lwork , m , mnmin , n , nb , ncols , nerrs ,      &
     &        nfail , nrhs , nrows , nrun , rank , mb , mmax , nmax ,   &
     &        nsmax , liwork , lwork_dgels , lwork_dgetsls ,            &
     &        lwork_dgelss , lwork_dgelsy , lwork_dgelsd
      DOUBLE PRECISION eps , norma , normb , rcond
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4) , iseedy(4) , iwq(1)
      DOUBLE PRECISION result(NTESTS) , wq(1)
!     ..
!     .. Allocatable Arrays ..
      DOUBLE PRECISION , ALLOCATABLE  ::  work(:)
      INTEGER , ALLOCATABLE  ::  iwork(:)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DASUM , DLAMCH 
EXTERNAL    DASUM,   DLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    DAXPY,   DGELS
EXTERNAL    DGELSD,   DGELSS,   DGELSY,   DGEMM,   DLACPY,   DLARNV
EXTERNAL    DLASRT,   DSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , INT , LOG , MAX , MIN , SQRT
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , IOUnit
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , IOUnit , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Double precision'
      path(2:3) = 'LS'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
      eps = DLAMCH('Epsilon')
!
!     Threshold for rank estimation
!
      rcond = SQRT(eps) - (SQRT(eps)-eps)/2
!
!     Test the error exits
!
      CALL xlaenv(2,2)
      CALL xlaenv(9,SMLSIZ)
      IF ( Tsterr ) CALL derrls(path,Nout)
!
!     Print the header if NM = 0 or NN = 0 and THRESH = 0.
!
      IF ( (Nm==0 .OR. Nn==0) .AND. Thresh==ZERO ) CALL alahd(Nout,path)
      INFot = 0
      CALL xlaenv(2,2)
      CALL xlaenv(9,SMLSIZ)
!
!     Compute maximal workspace needed for all routines
!
      nmax = 0
      mmax = 0
      nsmax = 0
      DO i = 1 , Nm
         IF ( Mval(i)>mmax ) mmax = Mval(i)
      ENDDO
      DO i = 1 , Nn
         IF ( Nval(i)>nmax ) nmax = Nval(i)
      ENDDO
      DO i = 1 , Nns
         IF ( Nsval(i)>nsmax ) nsmax = Nsval(i)
      ENDDO
      m = mmax
      n = nmax
      nrhs = nsmax
      mnmin = MAX(MIN(m,n),1)
!
!     Compute workspace needed for routines
!     dqrt14, dqrt17 (two side cases), dqrt15 and dqrt12
!
      lwork = MAX(1,(m+n)*nrhs,(n+nrhs)*(m+2),(m+nrhs)*(n+2),           &
     &        MAX(m+mnmin,nrhs*mnmin,2*n+m),                            &
     &        MAX(m*n+4*mnmin+MAX(m,n),m*n+2*mnmin+4*n))
      liwork = 1
!
!     Iterate through all test cases and compute necessary workspace
!     sizes for ?GELS, ?GETSLS, ?GELSY, ?GELSS and ?GELSD routines.
!
      DO im = 1 , Nm
         m = Mval(im)
         lda = MAX(1,m)
         DO in = 1 , Nn
            n = Nval(in)
            mnmin = MAX(MIN(m,n),1)
            ldb = MAX(1,m,n)
            DO ins = 1 , Nns
               nrhs = Nsval(ins)
               DO irank = 1 , 2
                  DO iscale = 1 , 3
                     itype = (irank-1)*3 + iscale
                     IF ( Dotype(itype) ) THEN
                        IF ( irank==1 ) THEN
                           DO itran = 1 , 2
                              IF ( itran==1 ) THEN
                                 trans = 'N'
                              ELSE
                                 trans = 'T'
                              ENDIF
!
!                             Compute workspace needed for DGELS
                              CALL DGELS(trans,m,n,nrhs,A,lda,B,ldb,wq, &
     &                           -1,info)
                              lwork_dgels = INT(wq(1))
!                             Compute workspace needed for DGETSLS
                              CALL DGETSLS(trans,m,n,nrhs,A,lda,B,ldb,  &
     &                           wq,-1,info)
                              lwork_dgetsls = INT(wq(1))
                           ENDDO
                        ENDIF
!                       Compute workspace needed for DGELSY
                        CALL DGELSY(m,n,nrhs,A,lda,B,ldb,iwq,rcond,     &
     &                              crank,wq,-1,info)
                        lwork_dgelsy = INT(wq(1))
!                       Compute workspace needed for DGELSS
                        CALL DGELSS(m,n,nrhs,A,lda,B,ldb,S,rcond,crank, &
     &                              wq,-1,info)
                        lwork_dgelss = INT(wq(1))
!                       Compute workspace needed for DGELSD
                        CALL DGELSD(m,n,nrhs,A,lda,B,ldb,S,rcond,crank, &
     &                              wq,-1,iwq,info)
                        lwork_dgelsd = INT(wq(1))
!                       Compute LIWORK workspace needed for DGELSY and DGELSD
                        liwork = MAX(liwork,n,iwq(1))
!                       Compute LWORK workspace needed for all functions
                        lwork = MAX(lwork,lwork_dgels,lwork_dgetsls,    &
     &                          lwork_dgelsy,lwork_dgelss,lwork_dgelsd)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
      lwlsy = lwork
!
      ALLOCATE (work(lwork))
      ALLOCATE (iwork(liwork))
!
      DO im = 1 , Nm
         m = Mval(im)
         lda = MAX(1,m)
!
         DO in = 1 , Nn
            n = Nval(in)
            mnmin = MAX(MIN(m,n),1)
            ldb = MAX(1,m,n)
            mb = (mnmin+1)
!
            DO ins = 1 , Nns
               nrhs = Nsval(ins)
!
               DO irank = 1 , 2
                  DO iscale = 1 , 3
                     itype = (irank-1)*3 + iscale
                     IF ( Dotype(itype) ) THEN
!
                        IF ( irank==1 ) THEN
!
!                       Test DGELS
!
!                       Generate a matrix of scaling type ISCALE
!
                           CALL dqrt13(iscale,m,n,Copya,lda,norma,iseed)
                           DO inb = 1 , Nnb
                              nb = Nbval(inb)
                              CALL xlaenv(1,nb)
                              CALL xlaenv(3,Nxval(inb))
!
                              DO itran = 1 , 2
                                 IF ( itran==1 ) THEN
                                    trans = 'N'
                                    nrows = m
                                    ncols = n
                                 ELSE
                                    trans = 'T'
                                    nrows = n
                                    ncols = m
                                 ENDIF
                                 ldwork = MAX(1,ncols)
!
!                             Set up a consistent rhs
!
                                 IF ( ncols>0 ) THEN
                                    CALL DLARNV(2,iseed,ncols*nrhs,work)
                                    CALL DSCAL(ncols*nrhs,              &
     &                                 ONE/DBLE(ncols),work,1)
                                 ENDIF
                                 CALL DGEMM(trans,'No transpose',nrows, &
     &                              nrhs,ncols,ONE,Copya,lda,work,      &
     &                              ldwork,ZERO,B,ldb)
                                 CALL DLACPY('Full',nrows,nrhs,B,ldb,   &
     &                              Copyb,ldb)
!
!                             Solve LS or overdetermined system
!
                                 IF ( m>0 .AND. n>0 ) THEN
                                    CALL DLACPY('Full',m,n,Copya,lda,A, &
     &                                 lda)
                                    CALL DLACPY('Full',nrows,nrhs,Copyb,&
     &                                 ldb,B,ldb)
                                 ENDIF
                                 SRNamt = 'DGELS '
                                 CALL DGELS(trans,m,n,nrhs,A,lda,B,ldb, &
     &                              work,lwork,info)
                                 IF ( info/=0 )                         &
     &                                 CALL alaerh(path,'DGELS ',info,0,&
     &                                trans,m,n,nrhs,-1,nb,itype,nfail, &
     &                                nerrs,Nout)
!
!                             Check correctness of results
!
                                 ldwork = MAX(1,nrows)
                                 IF ( nrows>0 .AND. nrhs>0 )            &
     &                                CALL DLACPY('Full',nrows,nrhs,    &
     &                                Copyb,ldb,C,ldb)
                                 CALL dqrt16(trans,m,n,nrhs,Copya,lda,B,&
     &                              ldb,C,ldb,work,result(1))
!
                                 IF ( (itran==1 .AND. m>=n) .OR.        &
     &                                (itran==2 .AND. m<n) ) THEN
!
!                                Solving LS system
!
                                    result(2)                           &
     &                                 = dqrt17(trans,1,m,n,nrhs,Copya, &
     &                                 lda,B,ldb,Copyb,ldb,C,work,lwork)
                                 ELSE
!
!                                Solving overdetermined system
!
                                    result(2)                           &
     &                                 = dqrt14(trans,m,n,nrhs,Copya,   &
     &                                 lda,B,ldb,work,lwork)
                                 ENDIF
!
!                             Print information about the tests that
!                             did not pass the threshold.
!
                                 DO k = 1 , 2
                                    IF ( result(k)>=Thresh ) THEN
                                       IF ( nfail==0 .AND. nerrs==0 )   &
     &                                    CALL alahd(Nout,path)
                                       WRITE (Nout,FMT=99001) trans ,   &
     &                                    m , n , nrhs , nb , itype ,   &
     &                                    k , result(k)
                                       nfail = nfail + 1
                                    ENDIF
                                 ENDDO
                                 nrun = nrun + 2
                              ENDDO
                           ENDDO
!
!
!                       Test DGETSLS
!
!                       Generate a matrix of scaling type ISCALE
!
                           CALL dqrt13(iscale,m,n,Copya,lda,norma,iseed)
                           DO inb = 1 , Nnb
                              mb = Nbval(inb)
                              CALL xlaenv(1,mb)
                              DO imb = 1 , Nnb
                                 nb = Nbval(imb)
                                 CALL xlaenv(2,nb)
!
                                 DO itran = 1 , 2
                                    IF ( itran==1 ) THEN
                                       trans = 'N'
                                       nrows = m
                                       ncols = n
                                    ELSE
                                       trans = 'T'
                                       nrows = n
                                       ncols = m
                                    ENDIF
                                    ldwork = MAX(1,ncols)
!
!                             Set up a consistent rhs
!
                                    IF ( ncols>0 ) THEN
                                       CALL DLARNV(2,iseed,ncols*nrhs,  &
     &                                    work)
                                       CALL DSCAL(ncols*nrhs,           &
     &                                    ONE/DBLE(ncols),work,1)
                                    ENDIF
                                    CALL DGEMM(trans,'No transpose',    &
     &                                 nrows,nrhs,ncols,ONE,Copya,lda,  &
     &                                 work,ldwork,ZERO,B,ldb)
                                    CALL DLACPY('Full',nrows,nrhs,B,ldb,&
     &                                 Copyb,ldb)
!
!                             Solve LS or overdetermined system
!
                                    IF ( m>0 .AND. n>0 ) THEN
                                       CALL DLACPY('Full',m,n,Copya,lda,&
     &                                    A,lda)
                                       CALL DLACPY('Full',nrows,nrhs,   &
     &                                    Copyb,ldb,B,ldb)
                                    ENDIF
                                    SRNamt = 'DGETSLS '
                                    CALL DGETSLS(trans,m,n,nrhs,A,lda,B,&
     &                                 ldb,work,lwork,info)
                                    IF ( info/=0 )                      &
     &                                  CALL alaerh(path,'DGETSLS ',    &
     &                                 info,0,trans,m,n,nrhs,-1,nb,     &
     &                                 itype,nfail,nerrs,Nout)
!
!                             Check correctness of results
!
                                    ldwork = MAX(1,nrows)
                                    IF ( nrows>0 .AND. nrhs>0 )         &
     &                                 CALL DLACPY('Full',nrows,nrhs,   &
     &                                 Copyb,ldb,C,ldb)
                                    CALL dqrt16(trans,m,n,nrhs,Copya,   &
     &                                 lda,B,ldb,C,ldb,work,result(15))
!
                                    IF ( (itran==1 .AND. m>=n) .OR.     &
     &                                 (itran==2 .AND. m<n) ) THEN
!
!                                Solving LS system
!
                                       result(16)                       &
     &                                    = dqrt17(trans,1,m,n,nrhs,    &
     &                                    Copya,lda,B,ldb,Copyb,ldb,C,  &
     &                                    work,lwork)
                                    ELSE
!
!                                Solving overdetermined system
!
                                       result(16)                       &
     &                                    = dqrt14(trans,m,n,nrhs,Copya,&
     &                                    lda,B,ldb,work,lwork)
                                    ENDIF
!
!                             Print information about the tests that
!                             did not pass the threshold.
!
                                    DO k = 15 , 16
                                       IF ( result(k)>=Thresh ) THEN
                                         IF ( nfail==0 .AND. nerrs==0 ) &
     &                                      CALL alahd(Nout,path)
                                         WRITE (Nout,FMT=99003) trans , &
     &                                      m , n , nrhs , mb , nb ,    &
     &                                      itype , k , result(k)
                                         nfail = nfail + 1
                                       ENDIF
                                    ENDDO
                                    nrun = nrun + 2
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDIF
!
!                    Generate a matrix of scaling type ISCALE and rank
!                    type IRANK.
!
                        CALL dqrt15(iscale,irank,m,n,nrhs,Copya,lda,    &
     &                              Copyb,ldb,Copys,rank,norma,normb,   &
     &                              iseed,work,lwork)
!
!                    workspace used: MAX(M+MIN(M,N),NRHS*MIN(M,N),2*N+M)
!
                        ldwork = MAX(1,m)
!
!                    Loop for testing different block sizes.
!
                        DO inb = 1 , Nnb
                           nb = Nbval(inb)
                           CALL xlaenv(1,nb)
                           CALL xlaenv(3,Nxval(inb))
!
!                       Test DGELSY
!
!                       DGELSY:  Compute the minimum-norm solution X
!                       to min( norm( A * X - B ) )
!                       using the rank-revealing orthogonal
!                       factorization.
!
!                       Initialize vector IWORK.
!
                           DO j = 1 , n
                              iwork(j) = 0
                           ENDDO
!
                           CALL DLACPY('Full',m,n,Copya,lda,A,lda)
                           CALL DLACPY('Full',m,nrhs,Copyb,ldb,B,ldb)
!
                           SRNamt = 'DGELSY'
                           CALL DGELSY(m,n,nrhs,A,lda,B,ldb,iwork,rcond,&
     &                                 crank,work,lwlsy,info)
                           IF ( info/=0 )                               &
     &                          CALL alaerh(path,'DGELSY',info,0,' ',m, &
     &                          n,nrhs,-1,nb,itype,nfail,nerrs,Nout)
!
!                       Test 3:  Compute relative error in svd
!                                workspace: M*N + 4*MIN(M,N) + MAX(M,N)
!
                           result(3) = dqrt12(crank,crank,A,lda,Copys,  &
     &                                 work,lwork)
!
!                       Test 4:  Compute error in solution
!                                workspace:  M*NRHS + M
!
                           CALL DLACPY('Full',m,nrhs,Copyb,ldb,work,    &
     &                                 ldwork)
                           CALL dqrt16('No transpose',m,n,nrhs,Copya,   &
     &                                 lda,B,ldb,work,ldwork,           &
     &                                 work(m*nrhs+1),result(4))
!
!                       Test 5:  Check norm of r'*A
!                                workspace: NRHS*(M+N)
!
                           result(5) = ZERO
                           IF ( m>crank ) result(5)                     &
     &                           = dqrt17('No transpose',1,m,n,nrhs,    &
     &                          Copya,lda,B,ldb,Copyb,ldb,C,work,lwork)
!
!                       Test 6:  Check if x is in the rowspace of A
!                                workspace: (M+NRHS)*(N+2)
!
                           result(6) = ZERO
!
                           IF ( n>crank ) result(6)                     &
     &                           = dqrt14('No transpose',m,n,nrhs,Copya,&
     &                          lda,B,ldb,work,lwork)
!
!                       Test DGELSS
!
!                       DGELSS:  Compute the minimum-norm solution X
!                       to min( norm( A * X - B ) )
!                       using the SVD.
!
                           CALL DLACPY('Full',m,n,Copya,lda,A,lda)
                           CALL DLACPY('Full',m,nrhs,Copyb,ldb,B,ldb)
                           SRNamt = 'DGELSS'
                           CALL DGELSS(m,n,nrhs,A,lda,B,ldb,S,rcond,    &
     &                                 crank,work,lwork,info)
                           IF ( info/=0 )                               &
     &                          CALL alaerh(path,'DGELSS',info,0,' ',m, &
     &                          n,nrhs,-1,nb,itype,nfail,nerrs,Nout)
!
!                       workspace used: 3*min(m,n) +
!                                       max(2*min(m,n),nrhs,max(m,n))
!
!                       Test 7:  Compute relative error in svd
!
                           IF ( rank>0 ) THEN
                              CALL DAXPY(mnmin,-ONE,Copys,1,S,1)
                              result(7) = DASUM(mnmin,S,1)              &
     &                           /DASUM(mnmin,Copys,1)/(eps*DBLE(mnmin))
                           ELSE
                              result(7) = ZERO
                           ENDIF
!
!                       Test 8:  Compute error in solution
!
                           CALL DLACPY('Full',m,nrhs,Copyb,ldb,work,    &
     &                                 ldwork)
                           CALL dqrt16('No transpose',m,n,nrhs,Copya,   &
     &                                 lda,B,ldb,work,ldwork,           &
     &                                 work(m*nrhs+1),result(8))
!
!                       Test 9:  Check norm of r'*A
!
                           result(9) = ZERO
                           IF ( m>crank ) result(9)                     &
     &                           = dqrt17('No transpose',1,m,n,nrhs,    &
     &                          Copya,lda,B,ldb,Copyb,ldb,C,work,lwork)
!
!                       Test 10:  Check if x is in the rowspace of A
!
                           result(10) = ZERO
                           IF ( n>crank ) result(10)                    &
     &                           = dqrt14('No transpose',m,n,nrhs,Copya,&
     &                          lda,B,ldb,work,lwork)
!
!                       Test DGELSD
!
!                       DGELSD:  Compute the minimum-norm solution X
!                       to min( norm( A * X - B ) ) using a
!                       divide and conquer SVD.
!
!                       Initialize vector IWORK.
!
                           DO j = 1 , n
                              iwork(j) = 0
                           ENDDO
!
                           CALL DLACPY('Full',m,n,Copya,lda,A,lda)
                           CALL DLACPY('Full',m,nrhs,Copyb,ldb,B,ldb)
!
                           SRNamt = 'DGELSD'
                           CALL DGELSD(m,n,nrhs,A,lda,B,ldb,S,rcond,    &
     &                                 crank,work,lwork,iwork,info)
                           IF ( info/=0 )                               &
     &                          CALL alaerh(path,'DGELSD',info,0,' ',m, &
     &                          n,nrhs,-1,nb,itype,nfail,nerrs,Nout)
!
!                       Test 11:  Compute relative error in svd
!
                           IF ( rank>0 ) THEN
                              CALL DAXPY(mnmin,-ONE,Copys,1,S,1)
                              result(11) = DASUM(mnmin,S,1)             &
     &                           /DASUM(mnmin,Copys,1)/(eps*DBLE(mnmin))
                           ELSE
                              result(11) = ZERO
                           ENDIF
!
!                       Test 12:  Compute error in solution
!
                           CALL DLACPY('Full',m,nrhs,Copyb,ldb,work,    &
     &                                 ldwork)
                           CALL dqrt16('No transpose',m,n,nrhs,Copya,   &
     &                                 lda,B,ldb,work,ldwork,           &
     &                                 work(m*nrhs+1),result(12))
!
!                       Test 13:  Check norm of r'*A
!
                           result(13) = ZERO
                           IF ( m>crank ) result(13)                    &
     &                           = dqrt17('No transpose',1,m,n,nrhs,    &
     &                          Copya,lda,B,ldb,Copyb,ldb,C,work,lwork)
!
!                       Test 14:  Check if x is in the rowspace of A
!
                           result(14) = ZERO
                           IF ( n>crank ) result(14)                    &
     &                           = dqrt14('No transpose',m,n,nrhs,Copya,&
     &                          lda,B,ldb,work,lwork)
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                           DO k = 3 , 14
                              IF ( result(k)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL alahd(Nout,path)
                                 WRITE (Nout,FMT=99002) m , n , nrhs ,  &
     &                                  nb , itype , k , result(k)
                                 nfail = nfail + 1
                              ENDIF
                           ENDDO
                           nrun = nrun + 12
!
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasvm(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' TRANS=''',A1,''', M=',I5,', N=',I5,', NRHS=',I4,', NB=',&
     &        I4,', type',I2,', test(',I2,')=',G12.5)
99002 FORMAT (' M=',I5,', N=',I5,', NRHS=',I4,', NB=',I4,', type',I2,   &
     &        ', test(',I2,')=',G12.5)
99003 FORMAT (' TRANS=''',A1,' M=',I5,', N=',I5,', NRHS=',I4,', MB=',I4,&
     &        ', NB=',I4,', type',I2,', test(',I2,')=',G12.5)
!
      DEALLOCATE (work)
      DEALLOCATE (iwork)
!
!     End of ddrvls
!
      END SUBROUTINE ddrvls
