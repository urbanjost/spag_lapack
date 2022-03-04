!*==zdrvls.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZDRVLS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDRVLS( DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB,
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
!       DOUBLE PRECISION   COPYS( * ), S( * )
!       COMPLEX*16         A( * ), B( * ), C( * ), COPYA( * ), COPYB( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZDRVLS tests the least squares driver routines ZGELS, ZGETSLS, ZGELSS, ZGELSY
!> and ZGELSD.
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
!>          j=1: A = U*D*V where U and V are random unitary matrices
!>               and D has random entries (> 0.1) taken from a uniform
!>               distribution (0,1). A is full rank.
!>          j=2: The same of 1, but A is scaled up.
!>          j=3: The same of 1, but A is scaled down.
!>          j=4: A = U*D*V where U and V are random unitary matrices
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
!>          A is COMPLEX*16 array, dimension (MMAX*NMAX)
!>          where MMAX is the maximum value of M in MVAL and NMAX is the
!>          maximum value of N in NVAL.
!> \endverbatim
!>
!> \param[out] COPYA
!> \verbatim
!>          COPYA is COMPLEX*16 array, dimension (MMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (MMAX*NSMAX)
!>          where MMAX is the maximum value of M in MVAL and NSMAX is the
!>          maximum value of NRHS in NSVAL.
!> \endverbatim
!>
!> \param[out] COPYB
!> \verbatim
!>          COPYB is COMPLEX*16 array, dimension (MMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (MMAX*NSMAX)
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZDRVLS(Dotype,Nm,Mval,Nn,Nval,Nns,Nsval,Nnb,Nbval,     &
     &                  Nxval,Thresh,Tsterr,A,Copya,B,Copyb,C,S,Copys,  &
     &                  Nout)
      IMPLICIT NONE
!*--ZDRVLS196
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
      DOUBLE PRECISION Copys(*) , S(*)
      COMPLEX*16 A(*) , B(*) , C(*) , Copya(*) , Copyb(*)
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
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      COMPLEX*16 CONE , CZERO
      PARAMETER (CONE=(1.0D+0,0.0D+0),CZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      CHARACTER trans
      CHARACTER*3 path
      INTEGER crank , i , im , imb , in , inb , info , ins , irank ,    &
     &        iscale , itran , itype , j , k , lda , ldb , ldwork ,     &
     &        lwlsy , lwork , m , mnmin , n , nb , ncols , nerrs ,      &
     &        nfail , nrhs , nrows , nrun , rank , mb , mmax , nmax ,   &
     &        nsmax , liwork , lrwork , lwork_zgels , lwork_zgetsls ,   &
     &        lwork_zgelss , lwork_zgelsy , lwork_zgelsd ,              &
     &        lrwork_zgelsy , lrwork_zgelss , lrwork_zgelsd
      DOUBLE PRECISION eps , norma , normb , rcond
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4) , iseedy(4) , iwq(1)
      DOUBLE PRECISION result(NTESTS) , rwq(1)
      COMPLEX*16 wq(1)
!     ..
!     .. Allocatable Arrays ..
      COMPLEX*16 , ALLOCATABLE  ::  work(:)
      DOUBLE PRECISION , ALLOCATABLE  ::  rwork(:) , work2(:)
      INTEGER , ALLOCATABLE  ::  iwork(:)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DASUM , DLAMCH , ZQRT12 , ZQRT14 , ZQRT17
      EXTERNAL DASUM , DLAMCH , ZQRT12 , ZQRT14 , ZQRT17
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAERH , ALAHD , ALASVM , DAXPY , DLASRT , XLAENV ,      &
     &         ZDSCAL , ZERRLS , ZGELS , ZGELSD , ZGELSS , ZGELSY ,     &
     &         ZGEMM , ZLACPY , ZLARNV , ZQRT13 , ZQRT15 , ZQRT16 ,     &
     &         ZGETSLS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX , MIN , INT , SQRT
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
      path(1:1) = 'Zomplex precision'
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
      CALL XLAENV(9,SMLSIZ)
      IF ( Tsterr ) CALL ZERRLS(path,Nout)
!
!     Print the header if NM = 0 or NN = 0 and THRESH = 0.
!
      IF ( (Nm==0 .OR. Nn==0) .AND. Thresh==ZERO ) CALL ALAHD(Nout,path)
      INFot = 0
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
!     ZQRT14, ZQRT17 (two side cases), ZQRT15 and ZQRT12
!
      lwork = MAX(1,(m+n)*nrhs,(n+nrhs)*(m+2),(m+nrhs)*(n+2),           &
     &        MAX(m+mnmin,nrhs*mnmin,2*n+m),                            &
     &        MAX(m*n+4*mnmin+MAX(m,n),m*n+2*mnmin+4*n))
      lrwork = 1
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
                                 trans = 'C'
                              ENDIF
!
!                             Compute workspace needed for ZGELS
                              CALL ZGELS(trans,m,n,nrhs,A,lda,B,ldb,wq, &
     &                           -1,info)
                              lwork_zgels = INT(wq(1))
!                             Compute workspace needed for ZGETSLS
                              CALL ZGETSLS(trans,m,n,nrhs,A,lda,B,ldb,  &
     &                           wq,-1,info)
                              lwork_zgetsls = INT(wq(1))
                           ENDDO
                        ENDIF
!                       Compute workspace needed for ZGELSY
                        CALL ZGELSY(m,n,nrhs,A,lda,B,ldb,iwq,rcond,     &
     &                              crank,wq,-1,rwq,info)
                        lwork_zgelsy = INT(wq(1))
                        lrwork_zgelsy = 2*n
!                       Compute workspace needed for ZGELSS
                        CALL ZGELSS(m,n,nrhs,A,lda,B,ldb,S,rcond,crank, &
     &                              wq,-1,rwq,info)
                        lwork_zgelss = INT(wq(1))
                        lrwork_zgelss = 5*mnmin
!                       Compute workspace needed for ZGELSD
                        CALL ZGELSD(m,n,nrhs,A,lda,B,ldb,S,rcond,crank, &
     &                              wq,-1,rwq,iwq,info)
                        lwork_zgelsd = INT(wq(1))
                        lrwork_zgelsd = INT(rwq(1))
!                       Compute LIWORK workspace needed for ZGELSY and ZGELSD
                        liwork = MAX(liwork,n,iwq(1))
!                       Compute LRWORK workspace needed for ZGELSY, ZGELSS and ZGELSD
                        lrwork = MAX(lrwork,lrwork_zgelsy,lrwork_zgelss,&
     &                           lrwork_zgelsd)
!                       Compute LWORK workspace needed for all functions
                        lwork = MAX(lwork,lwork_zgels,lwork_zgetsls,    &
     &                          lwork_zgelsy,lwork_zgelss,lwork_zgelsd)
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
      ALLOCATE (work2(2*lwork))
      ALLOCATE (iwork(liwork))
      ALLOCATE (rwork(lrwork))
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
!                       Test ZGELS
!
!                       Generate a matrix of scaling type ISCALE
!
                           CALL ZQRT13(iscale,m,n,Copya,lda,norma,iseed)
                           DO inb = 1 , Nnb
                              nb = Nbval(inb)
                              CALL XLAENV(1,nb)
                              CALL XLAENV(3,Nxval(inb))
!
                              DO itran = 1 , 2
                                 IF ( itran==1 ) THEN
                                    trans = 'N'
                                    nrows = m
                                    ncols = n
                                 ELSE
                                    trans = 'C'
                                    nrows = n
                                    ncols = m
                                 ENDIF
                                 ldwork = MAX(1,ncols)
!
!                             Set up a consistent rhs
!
                                 IF ( ncols>0 ) THEN
                                    CALL ZLARNV(2,iseed,ncols*nrhs,work)
                                    CALL ZDSCAL(ncols*nrhs,             &
     &                                 ONE/DBLE(ncols),work,1)
                                 ENDIF
                                 CALL ZGEMM(trans,'No transpose',nrows, &
     &                              nrhs,ncols,CONE,Copya,lda,work,     &
     &                              ldwork,CZERO,B,ldb)
                                 CALL ZLACPY('Full',nrows,nrhs,B,ldb,   &
     &                              Copyb,ldb)
!
!                             Solve LS or overdetermined system
!
                                 IF ( m>0 .AND. n>0 ) THEN
                                    CALL ZLACPY('Full',m,n,Copya,lda,A, &
     &                                 lda)
                                    CALL ZLACPY('Full',nrows,nrhs,Copyb,&
     &                                 ldb,B,ldb)
                                 ENDIF
                                 SRNamt = 'ZGELS '
                                 CALL ZGELS(trans,m,n,nrhs,A,lda,B,ldb, &
     &                              work,lwork,info)
!
                                 IF ( info/=0 )                         &
     &                                 CALL ALAERH(path,'ZGELS ',info,0,&
     &                                trans,m,n,nrhs,-1,nb,itype,nfail, &
     &                                nerrs,Nout)
!
!                             Check correctness of results
!
                                 ldwork = MAX(1,nrows)
                                 IF ( nrows>0 .AND. nrhs>0 )            &
     &                                CALL ZLACPY('Full',nrows,nrhs,    &
     &                                Copyb,ldb,C,ldb)
                                 CALL ZQRT16(trans,m,n,nrhs,Copya,lda,B,&
     &                              ldb,C,ldb,rwork,result(1))
!
                                 IF ( (itran==1 .AND. m>=n) .OR.        &
     &                                (itran==2 .AND. m<n) ) THEN
!
!                                Solving LS system
!
                                    result(2)                           &
     &                                 = ZQRT17(trans,1,m,n,nrhs,Copya, &
     &                                 lda,B,ldb,Copyb,ldb,C,work,lwork)
                                 ELSE
!
!                                Solving overdetermined system
!
                                    result(2)                           &
     &                                 = ZQRT14(trans,m,n,nrhs,Copya,   &
     &                                 lda,B,ldb,work,lwork)
                                 ENDIF
!
!                             Print information about the tests that
!                             did not pass the threshold.
!
                                 DO k = 1 , 2
                                    IF ( result(k)>=Thresh ) THEN
                                       IF ( nfail==0 .AND. nerrs==0 )   &
     &                                    CALL ALAHD(Nout,path)
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
!                       Test ZGETSLS
!
!                       Generate a matrix of scaling type ISCALE
!
                           CALL ZQRT13(iscale,m,n,Copya,lda,norma,iseed)
                           DO inb = 1 , Nnb
                              mb = Nbval(inb)
                              CALL XLAENV(1,mb)
                              DO imb = 1 , Nnb
                                 nb = Nbval(imb)
                                 CALL XLAENV(2,nb)
!
                                 DO itran = 1 , 2
                                    IF ( itran==1 ) THEN
                                       trans = 'N'
                                       nrows = m
                                       ncols = n
                                    ELSE
                                       trans = 'C'
                                       nrows = n
                                       ncols = m
                                    ENDIF
                                    ldwork = MAX(1,ncols)
!
!                             Set up a consistent rhs
!
                                    IF ( ncols>0 ) THEN
                                       CALL ZLARNV(2,iseed,ncols*nrhs,  &
     &                                    work)
                                       CALL ZSCAL(ncols*nrhs,           &
     &                                    CONE/DBLE(ncols),work,1)
                                    ENDIF
                                    CALL ZGEMM(trans,'No transpose',    &
     &                                 nrows,nrhs,ncols,CONE,Copya,lda, &
     &                                 work,ldwork,CZERO,B,ldb)
                                    CALL ZLACPY('Full',nrows,nrhs,B,ldb,&
     &                                 Copyb,ldb)
!
!                             Solve LS or overdetermined system
!
                                    IF ( m>0 .AND. n>0 ) THEN
                                       CALL ZLACPY('Full',m,n,Copya,lda,&
     &                                    A,lda)
                                       CALL ZLACPY('Full',nrows,nrhs,   &
     &                                    Copyb,ldb,B,ldb)
                                    ENDIF
                                    SRNamt = 'ZGETSLS '
                                    CALL ZGETSLS(trans,m,n,nrhs,A,lda,B,&
     &                                 ldb,work,lwork,info)
                                    IF ( info/=0 )                      &
     &                                  CALL ALAERH(path,'ZGETSLS ',    &
     &                                 info,0,trans,m,n,nrhs,-1,nb,     &
     &                                 itype,nfail,nerrs,Nout)
!
!                             Check correctness of results
!
                                    ldwork = MAX(1,nrows)
                                    IF ( nrows>0 .AND. nrhs>0 )         &
     &                                 CALL ZLACPY('Full',nrows,nrhs,   &
     &                                 Copyb,ldb,C,ldb)
                                    CALL ZQRT16(trans,m,n,nrhs,Copya,   &
     &                                 lda,B,ldb,C,ldb,work2,result(15))
!
                                    IF ( (itran==1 .AND. m>=n) .OR.     &
     &                                 (itran==2 .AND. m<n) ) THEN
!
!                                Solving LS system
!
                                       result(16)                       &
     &                                    = ZQRT17(trans,1,m,n,nrhs,    &
     &                                    Copya,lda,B,ldb,Copyb,ldb,C,  &
     &                                    work,lwork)
                                    ELSE
!
!                                Solving overdetermined system
!
                                       result(16)                       &
     &                                    = ZQRT14(trans,m,n,nrhs,Copya,&
     &                                    lda,B,ldb,work,lwork)
                                    ENDIF
!
!                             Print information about the tests that
!                             did not pass the threshold.
!
                                    DO k = 15 , 16
                                       IF ( result(k)>=Thresh ) THEN
                                         IF ( nfail==0 .AND. nerrs==0 ) &
     &                                      CALL ALAHD(Nout,path)
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
                        CALL ZQRT15(iscale,irank,m,n,nrhs,Copya,lda,    &
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
                           CALL XLAENV(1,nb)
                           CALL XLAENV(3,Nxval(inb))
!
!                       Test ZGELSY
!
!                       ZGELSY:  Compute the minimum-norm solution
!                       X to min( norm( A * X - B ) )
!                       using the rank-revealing orthogonal
!                       factorization.
!
                           CALL ZLACPY('Full',m,n,Copya,lda,A,lda)
                           CALL ZLACPY('Full',m,nrhs,Copyb,ldb,B,ldb)
!
!                       Initialize vector IWORK.
!
                           DO j = 1 , n
                              iwork(j) = 0
                           ENDDO
!
                           SRNamt = 'ZGELSY'
                           CALL ZGELSY(m,n,nrhs,A,lda,B,ldb,iwork,rcond,&
     &                                 crank,work,lwlsy,rwork,info)
                           IF ( info/=0 )                               &
     &                          CALL ALAERH(path,'ZGELSY',info,0,' ',m, &
     &                          n,nrhs,-1,nb,itype,nfail,nerrs,Nout)
!
!                       workspace used: 2*MNMIN+NB*NB+NB*MAX(N,NRHS)
!
!                       Test 3:  Compute relative error in svd
!                                workspace: M*N + 4*MIN(M,N) + MAX(M,N)
!
                           result(3) = ZQRT12(crank,crank,A,lda,Copys,  &
     &                                 work,lwork,rwork)
!
!                       Test 4:  Compute error in solution
!                                workspace:  M*NRHS + M
!
                           CALL ZLACPY('Full',m,nrhs,Copyb,ldb,work,    &
     &                                 ldwork)
                           CALL ZQRT16('No transpose',m,n,nrhs,Copya,   &
     &                                 lda,B,ldb,work,ldwork,rwork,     &
     &                                 result(4))
!
!                       Test 5:  Check norm of r'*A
!                                workspace: NRHS*(M+N)
!
                           result(5) = ZERO
                           IF ( m>crank ) result(5)                     &
     &                           = ZQRT17('No transpose',1,m,n,nrhs,    &
     &                          Copya,lda,B,ldb,Copyb,ldb,C,work,lwork)
!
!                       Test 6:  Check if x is in the rowspace of A
!                                workspace: (M+NRHS)*(N+2)
!
                           result(6) = ZERO
!
                           IF ( n>crank ) result(6)                     &
     &                           = ZQRT14('No transpose',m,n,nrhs,Copya,&
     &                          lda,B,ldb,work,lwork)
!
!                       Test ZGELSS
!
!                       ZGELSS:  Compute the minimum-norm solution
!                       X to min( norm( A * X - B ) )
!                       using the SVD.
!
                           CALL ZLACPY('Full',m,n,Copya,lda,A,lda)
                           CALL ZLACPY('Full',m,nrhs,Copyb,ldb,B,ldb)
                           SRNamt = 'ZGELSS'
                           CALL ZGELSS(m,n,nrhs,A,lda,B,ldb,S,rcond,    &
     &                                 crank,work,lwork,rwork,info)
!
                           IF ( info/=0 )                               &
     &                          CALL ALAERH(path,'ZGELSS',info,0,' ',m, &
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
                           CALL ZLACPY('Full',m,nrhs,Copyb,ldb,work,    &
     &                                 ldwork)
                           CALL ZQRT16('No transpose',m,n,nrhs,Copya,   &
     &                                 lda,B,ldb,work,ldwork,rwork,     &
     &                                 result(8))
!
!                       Test 9:  Check norm of r'*A
!
                           result(9) = ZERO
                           IF ( m>crank ) result(9)                     &
     &                           = ZQRT17('No transpose',1,m,n,nrhs,    &
     &                          Copya,lda,B,ldb,Copyb,ldb,C,work,lwork)
!
!                       Test 10:  Check if x is in the rowspace of A
!
                           result(10) = ZERO
                           IF ( n>crank ) result(10)                    &
     &                           = ZQRT14('No transpose',m,n,nrhs,Copya,&
     &                          lda,B,ldb,work,lwork)
!
!                       Test ZGELSD
!
!                       ZGELSD:  Compute the minimum-norm solution X
!                       to min( norm( A * X - B ) ) using a
!                       divide and conquer SVD.
!
                           CALL XLAENV(9,25)
!
                           CALL ZLACPY('Full',m,n,Copya,lda,A,lda)
                           CALL ZLACPY('Full',m,nrhs,Copyb,ldb,B,ldb)
!
                           SRNamt = 'ZGELSD'
                           CALL ZGELSD(m,n,nrhs,A,lda,B,ldb,S,rcond,    &
     &                                 crank,work,lwork,rwork,iwork,    &
     &                                 info)
                           IF ( info/=0 )                               &
     &                          CALL ALAERH(path,'ZGELSD',info,0,' ',m, &
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
                           CALL ZLACPY('Full',m,nrhs,Copyb,ldb,work,    &
     &                                 ldwork)
                           CALL ZQRT16('No transpose',m,n,nrhs,Copya,   &
     &                                 lda,B,ldb,work,ldwork,rwork,     &
     &                                 result(12))
!
!                       Test 13:  Check norm of r'*A
!
                           result(13) = ZERO
                           IF ( m>crank ) result(13)                    &
     &                           = ZQRT17('No transpose',1,m,n,nrhs,    &
     &                          Copya,lda,B,ldb,Copyb,ldb,C,work,lwork)
!
!                       Test 14:  Check if x is in the rowspace of A
!
                           result(14) = ZERO
                           IF ( n>crank ) result(14)                    &
     &                           = ZQRT14('No transpose',m,n,nrhs,Copya,&
     &                          lda,B,ldb,work,lwork)
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                           DO k = 3 , 14
                              IF ( result(k)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL ALAHD(Nout,path)
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
      CALL ALASVM(path,Nout,nfail,nrun,nerrs)
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
      DEALLOCATE (rwork)
!
!     End of ZDRVLS
!
      END SUBROUTINE ZDRVLS
