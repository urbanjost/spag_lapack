!*==cchklq.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CCHKLQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKLQ( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL,
!                          NRHS, THRESH, TSTERR, NMAX, A, AF, AQ, AL, AC,
!                          B, X, XACT, TAU, WORK, RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NM, NMAX, NN, NNB, NOUT, NRHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            MVAL( * ), NBVAL( * ), NVAL( * ),
!      $                   NXVAL( * )
!       REAL               RWORK( * )
!       COMPLEX            A( * ), AC( * ), AF( * ), AL( * ), AQ( * ),
!      $                   B( * ), TAU( * ), WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCHKLQ tests CGELQF, CUNGLQ and CUNMLQ.
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
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand side vectors to be generated for
!>          each linear system.
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
!> \param[in] TSTERR
!> \verbatim
!>          TSTERR is LOGICAL
!>          Flag that indicates whether error exits are to be tested.
!> \endverbatim
!>
!> \param[in] NMAX
!> \verbatim
!>          NMAX is INTEGER
!>          The maximum value permitted for M or N, used in dimensioning
!>          the work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AQ
!> \verbatim
!>          AQ is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AL
!> \verbatim
!>          AL is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AC
!> \verbatim
!>          AC is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (NMAX)
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
      SUBROUTINE CCHKLQ(Dotype,Nm,Mval,Nn,Nval,Nnb,Nbval,Nxval,Nrhs,    &
     &                  Thresh,Tsterr,Nmax,A,Af,Aq,Al,Ac,B,X,Xact,Tau,  &
     &                  Work,Rwork,Nout)
      IMPLICIT NONE
!*--CCHKLQ200
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nm , Nmax , Nn , Nnb , Nout , Nrhs
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Mval(*) , Nbval(*) , Nval(*) , Nxval(*)
      REAL Rwork(*)
      COMPLEX A(*) , Ac(*) , Af(*) , Al(*) , Aq(*) , B(*) , Tau(*) ,    &
     &        Work(*) , X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NTESTS
      PARAMETER (NTESTS=7)
      INTEGER NTYPES
      PARAMETER (NTYPES=8)
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      CHARACTER dist , type
      CHARACTER*3 path
      INTEGER i , ik , im , imat , in , inb , info , k , kl , ku , lda ,&
     &        lwork , m , minmn , mode , n , nb , nerrs , nfail , nk ,  &
     &        nrun , nt , nx
      REAL anorm , cndnum
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4) , iseedy(4) , kval(4)
      REAL result(NTESTS)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAERH , ALAHD , ALASUM , CERRLQ , CGELQS , CGET02 ,     &
     &         CLACPY , CLARHS , CLATB4 , CLATMS , CLQT01 , CLQT02 ,    &
     &         CLQT03 , XLAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NUNit
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NUNit , OK , LERr
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
      path(2:3) = 'LQ'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL CERRLQ(path,Nout)
      INFot = 0
      CALL XLAENV(2,2)
!
      lda = Nmax
      lwork = Nmax*MAX(Nmax,Nrhs)
!
!     Do for each value of M in MVAL.
!
      DO im = 1 , Nm
         m = Mval(im)
!
!        Do for each value of N in NVAL.
!
         DO in = 1 , Nn
            n = Nval(in)
            minmn = MIN(m,n)
            DO imat = 1 , NTYPES
!
!              Do the tests only if DOTYPE( IMAT ) is true.
!
               IF ( Dotype(imat) ) THEN
!
!              Set up parameters with CLATB4 and generate a test matrix
!              with CLATMS.
!
                  CALL CLATB4(path,imat,m,n,type,kl,ku,anorm,mode,      &
     &                        cndnum,dist)
!
                  SRNamt = 'CLATMS'
                  CALL CLATMS(m,n,dist,iseed,type,Rwork,mode,cndnum,    &
     &                        anorm,kl,ku,'No packing',A,lda,Work,info)
!
!              Check error code from CLATMS.
!
                  IF ( info/=0 ) THEN
                     CALL ALAERH(path,'CLATMS',info,0,' ',m,n,-1,-1,-1, &
     &                           imat,nfail,nerrs,Nout)
                     CYCLE
                  ENDIF
!
!              Set some values for K: the first value must be MINMN,
!              corresponding to the call of CLQT01; other values are
!              used in the calls of CLQT02, and must not exceed MINMN.
!
                  kval(1) = minmn
                  kval(2) = 0
                  kval(3) = 1
                  kval(4) = minmn/2
                  IF ( minmn==0 ) THEN
                     nk = 1
                  ELSEIF ( minmn==1 ) THEN
                     nk = 2
                  ELSEIF ( minmn<=3 ) THEN
                     nk = 3
                  ELSE
                     nk = 4
                  ENDIF
!
!              Do for each value of K in KVAL
!
                  DO ik = 1 , nk
                     k = kval(ik)
!
!                 Do for each pair of values (NB,NX) in NBVAL and NXVAL.
!
                     DO inb = 1 , Nnb
                        nb = Nbval(inb)
                        CALL XLAENV(1,nb)
                        nx = Nxval(inb)
                        CALL XLAENV(3,nx)
                        DO i = 1 , NTESTS
                           result(i) = ZERO
                        ENDDO
                        nt = 2
                        IF ( ik==1 ) THEN
!
!                       Test CGELQF
!
                           CALL CLQT01(m,n,A,Af,Aq,Al,lda,Tau,Work,     &
     &                                 lwork,Rwork,result(1))
                        ELSEIF ( m<=n ) THEN
!
!                       Test CUNGLQ, using factorization
!                       returned by CLQT01
!
                           CALL CLQT02(m,n,k,A,Af,Aq,Al,lda,Tau,Work,   &
     &                                 lwork,Rwork,result(1))
                        ENDIF
                        IF ( m>=k ) THEN
!
!                       Test CUNMLQ, using factorization returned
!                       by CLQT01
!
                           CALL CLQT03(m,n,k,Af,Ac,Al,Aq,lda,Tau,Work,  &
     &                                 lwork,Rwork,result(3))
                           nt = nt + 4
!
!                       If M>=N and K=N, call CGELQS to solve a system
!                       with NRHS right hand sides and compute the
!                       residual.
!
                           IF ( k==m .AND. inb==1 ) THEN
!
!                          Generate a solution and set the right
!                          hand side.
!
                              SRNamt = 'CLARHS'
                              CALL CLARHS(path,'New','Full',            &
     &                           'No transpose',m,n,0,0,Nrhs,A,lda,Xact,&
     &                           lda,B,lda,iseed,info)
!
                              CALL CLACPY('Full',m,Nrhs,B,lda,X,lda)
                              SRNamt = 'CGELQS'
                              CALL CGELQS(m,n,Nrhs,Af,lda,Tau,X,lda,    &
     &                           Work,lwork,info)
!
!                          Check error code from CGELQS.
!
                              IF ( info/=0 )                            &
     &                             CALL ALAERH(path,'CGELQS',info,0,' ',&
     &                             m,n,Nrhs,-1,nb,imat,nfail,nerrs,Nout)
!
                              CALL CGET02('No transpose',m,n,Nrhs,A,lda,&
     &                           X,lda,B,lda,Rwork,result(7))
                              nt = nt + 1
                           ENDIF
                        ENDIF
!
!                    Print information about the tests that did not
!                    pass the threshold.
!
                        DO i = 1 , nt
                           IF ( result(i)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL ALAHD(Nout,path)
                              WRITE (Nout,FMT=99001) m , n , k , nb ,   &
     &                               nx , imat , i , result(i)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + nt
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL ALASUM(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' M=',I5,', N=',I5,', K=',I5,', NB=',I4,', NX=',I5,       &
     &        ', type ',I2,', test(',I2,')=',G12.5)
!
!     End of CCHKLQ
!
      END SUBROUTINE CCHKLQ
