!*==cchkq3.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cchkq3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cchkq3( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL,
!                          THRESH, A, COPYA, S, TAU, WORK, RWORK,
!                          IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       INTEGER            NM, NN, NNB, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), MVAL( * ), NBVAL( * ), NVAL( * ),
!      $                   NXVAL( * )
!       REAL               S( * ), RWORK( * )
!       COMPLEX            A( * ), COPYA( * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cchkq3 tests CGEQP3.
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
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  To have
!>          every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (MMAX*NMAX)
!>          where MMAX is the maximum value of M in MVAL and NMAX is the
!>          maximum value of N in NVAL.
!> \endverbatim
!>
!> \param[out] COPYA
!> \verbatim
!>          COPYA is COMPLEX array, dimension (MMAX*NMAX)
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension
!>                      (min(MMAX,NMAX))
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (MMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension
!>                      (max(M*max(M,N) + 4*min(M,N) + max(M,N)))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (4*NMAX)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (2*NMAX)
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
!> \date December 2016
!
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cchkq3(Dotype,Nm,Mval,Nn,Nval,Nnb,Nbval,Nxval,Thresh,A,&
     &                  Copya,S,Tau,Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--CCHKQ3161
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Nm , Nn , Nnb , Nout
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Mval(*) , Nbval(*) , Nval(*) , Nxval(*)
      REAL S(*) , Rwork(*)
      COMPLEX A(*) , Copya(*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NTYPES
      PARAMETER (NTYPES=6)
      INTEGER NTESTS
      PARAMETER (NTESTS=3)
      REAL ONE , ZERO
      COMPLEX CZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0,CZERO=(0.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      CHARACTER*3 path
      INTEGER i , ihigh , ilow , im , imode , in , inb , info , istep , &
     &        k , lda , lw , lwork , m , mnmin , mode , n , nb , nerrs ,&
     &        nfail , nrun , nx
      REAL eps
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4) , iseedy(4)
      REAL result(NTESTS)
!     ..
!     .. External Functions ..
      REAL SLAMCH
EXTERNAL    SLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    CGEQP3,   CLACPY,   CLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
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
      path(1:1) = 'Complex precision'
      path(2:3) = 'Q3'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
      eps = SLAMCH('Epsilon')
      INFot = 0
!
      DO im = 1 , Nm
!
!        Do for each value of M in MVAL.
!
         m = Mval(im)
         lda = MAX(1,m)
!
         DO in = 1 , Nn
!
!           Do for each value of N in NVAL.
!
            n = Nval(in)
            mnmin = MIN(m,n)
            lwork = MAX(1,m*MAX(m,n)+4*mnmin+MAX(m,n))
!
            DO imode = 1 , NTYPES
               IF ( Dotype(imode) ) THEN
!
!              Do for each type of matrix
!                 1:  zero matrix
!                 2:  one small singular value
!                 3:  geometric distribution of singular values
!                 4:  first n/2 columns fixed
!                 5:  last n/2 columns fixed
!                 6:  every second column fixed
!
                  mode = imode
                  IF ( imode>3 ) mode = 1
!
!              Generate test matrix of size m by n using
!              singular value distribution indicated by `mode'.
!
                  DO i = 1 , n
                     Iwork(i) = 0
                  ENDDO
                  IF ( imode==1 ) THEN
                     CALL CLASET('Full',m,n,CZERO,CZERO,Copya,lda)
                     DO i = 1 , mnmin
                        S(i) = ZERO
                     ENDDO
                  ELSE
                     CALL CLATMS(m,n,'Uniform',iseed,'Nonsymm',S,mode,  &
     &                           ONE/eps,ONE,m,n,'No packing',Copya,lda,&
     &                           Work,info)
                     IF ( imode>=4 ) THEN
                        IF ( imode==4 ) THEN
                           ilow = 1
                           istep = 1
                           ihigh = MAX(1,n/2)
                        ELSEIF ( imode==5 ) THEN
                           ilow = MAX(1,n/2)
                           istep = 1
                           ihigh = n
                        ELSEIF ( imode==6 ) THEN
                           ilow = 1
                           istep = 2
                           ihigh = n
                        ENDIF
                        DO i = ilow , ihigh , istep
                           Iwork(i) = 1
                        ENDDO
                     ENDIF
                     CALL slaord('Decreasing',mnmin,S,1)
                  ENDIF
!
                  DO inb = 1 , Nnb
!
!                 Do for each pair of values (NB,NX) in NBVAL and NXVAL.
!
                     nb = Nbval(inb)
                     CALL xlaenv(1,nb)
                     nx = Nxval(inb)
                     CALL xlaenv(3,nx)
!
!                 Save A and its singular values and a copy of
!                 vector IWORK.
!
                     CALL CLACPY('All',m,n,Copya,lda,A,lda)
                     CALL icopy(n,Iwork(1),1,Iwork(n+1),1)
!
!                 Workspace needed.
!
                     lw = nb*(n+1)
!
                     SRNamt = 'CGEQP3'
                     CALL CGEQP3(m,n,A,lda,Iwork(n+1),Tau,Work,lw,Rwork,&
     &                           info)
!
!                 Compute norm(svd(a) - svd(r))
!
                     result(1) = cqrt12(m,n,A,lda,S,Work,lwork,Rwork)
!
!                 Compute norm( A*P - Q*R )
!
                     result(2) = cqpt01(m,n,mnmin,Copya,A,lda,Tau,      &
     &                           Iwork(n+1),Work,lwork)
!
!                 Compute Q'*Q
!
                     result(3) = cqrt11(m,mnmin,A,lda,Tau,Work,lwork)
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                     DO k = 1 , NTESTS
                        IF ( result(k)>=Thresh ) THEN
                           IF ( nfail==0 .AND. nerrs==0 )               &
     &                          CALL alahd(Nout,path)
                           WRITE (Nout,FMT=99001) 'CGEQP3' , m , n ,    &
     &                            nb , imode , k , result(k)
                           nfail = nfail + 1
                        ENDIF
                     ENDDO
                     nrun = nrun + NTESTS
!
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasum(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (1X,A,' M =',I5,', N =',I5,', NB =',I4,', type ',I2,       &
     &        ', test ',I2,', ratio =',G12.5)
!
!     End of cchkq3
!
      END SUBROUTINE cchkq3
