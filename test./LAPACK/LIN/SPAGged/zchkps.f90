!*==zchkps.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZCHKPS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKPS( DOTYPE, NN, NVAL, NNB, NBVAL, NRANK, RANKVAL,
!                          THRESH, TSTERR, NMAX, A, AFAC, PERM, PIV, WORK,
!                          RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   THRESH
!       INTEGER            NMAX, NN, NNB, NOUT, NRANK
!       LOGICAL            TSTERR
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( * ), AFAC( * ), PERM( * ), WORK( * )
!       DOUBLE PRECISION   RWORK( * )
!       INTEGER            NBVAL( * ), NVAL( * ), PIV( * ), RANKVAL( * )
!       LOGICAL            DOTYPE( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCHKPS tests ZPSTRF.
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
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER
!>          The number of values of N contained in the vector NVAL.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NN)
!>          The values of the matrix dimension N.
!> \endverbatim
!>
!> \param[in] NNB
!> \verbatim
!>          NNB is INTEGER
!>          The number of values of NB contained in the vector NBVAL.
!> \endverbatim
!>
!> \param[in] NBVAL
!> \verbatim
!>          NBVAL is INTEGER array, dimension (NBVAL)
!>          The values of the block size NB.
!> \endverbatim
!>
!> \param[in] NRANK
!> \verbatim
!>          NRANK is INTEGER
!>          The number of values of RANK contained in the vector RANKVAL.
!> \endverbatim
!>
!> \param[in] RANKVAL
!> \verbatim
!>          RANKVAL is INTEGER array, dimension (NBVAL)
!>          The values of the block size NB.
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
!> \param[in] NMAX
!> \verbatim
!>          NMAX is INTEGER
!>          The maximum value permitted for N, used in dimensioning the
!>          work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] PERM
!> \verbatim
!>          PERM is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] PIV
!> \verbatim
!>          PIV is INTEGER array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (NMAX*3)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (NMAX)
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZCHKPS(Dotype,Nn,Nval,Nnb,Nbval,Nrank,Rankval,Thresh,  &
     &                  Tsterr,Nmax,A,Afac,Perm,Piv,Work,Rwork,Nout)
      IMPLICIT NONE
!*--ZCHKPS157
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Thresh
      INTEGER Nmax , Nn , Nnb , Nout , Nrank
      LOGICAL Tsterr
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(*) , Afac(*) , Perm(*) , Work(*)
      DOUBLE PRECISION Rwork(*)
      INTEGER Nbval(*) , Nval(*) , Piv(*) , Rankval(*)
      LOGICAL Dotype(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0E+0)
      INTEGER NTYPES
      PARAMETER (NTYPES=9)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION anorm , cndnum , result , tol
      INTEGER comprank , i , imat , in , inb , info , irank , iuplo ,   &
     &        izero , kl , ku , lda , mode , n , nb , nerrs , nfail ,   &
     &        nimat , nrun , rank , rankdiff
      CHARACTER dist , type , uplo
      CHARACTER*3 path
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4) , iseedy(4)
      CHARACTER uplos(2)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAERH , ALAHD , ALASUM , XLAENV , ZERRPS , ZLACPY ,     &
     &         ZLATB5 , ZLATMT , ZPST01 , ZPSTRF
!     ..
!     .. Scalars in Common ..
      INTEGER INFot , NUNit
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NUNit , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX , CEILING
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
      DATA uplos/'U' , 'L'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Zomplex Precision'
      path(2:3) = 'PS'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL ZERRPS(path,Nout)
      INFot = 0
!
!     Do for each value of N in NVAL
!
      DO in = 1 , Nn
         n = Nval(in)
         lda = MAX(n,1)
         nimat = NTYPES
         IF ( n<=0 ) nimat = 1
!
         izero = 0
         DO imat = 1 , nimat
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
            IF ( Dotype(imat) ) THEN
!
!              Do for each value of RANK in RANKVAL
!
               DO irank = 1 , Nrank
!
!              Only repeat test 3 to 5 for different ranks
!              Other tests use full rank
!
                  IF ( .NOT.((imat<3 .OR. imat>5) .AND. irank>1) ) THEN
!
                     rank = CEILING((n*DBLE(Rankval(irank)))/100.E+0)
!
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
                     DO iuplo = 1 , 2
                        uplo = uplos(iuplo)
!
!              Set up parameters with ZLATB5 and generate a test matrix
!              with ZLATMT.
!
                        CALL ZLATB5(path,imat,n,type,kl,ku,anorm,mode,  &
     &                              cndnum,dist)
!
                        SRNamt = 'ZLATMT'
                        CALL ZLATMT(n,n,dist,iseed,type,Rwork,mode,     &
     &                              cndnum,anorm,rank,kl,ku,uplo,A,lda, &
     &                              Work,info)
!
!              Check error code from ZLATMT.
!
                        IF ( info/=0 ) THEN
                           CALL ALAERH(path,'ZLATMT',info,0,uplo,n,n,-1,&
     &                                 -1,-1,imat,nfail,nerrs,Nout)
                           CYCLE
                        ENDIF
!
!              Do for each value of NB in NBVAL
!
                        DO inb = 1 , Nnb
                           nb = Nbval(inb)
                           CALL XLAENV(1,nb)
!
!                 Compute the pivoted L*L' or U'*U factorization
!                 of the matrix.
!
                           CALL ZLACPY(uplo,n,n,A,lda,Afac,lda)
                           SRNamt = 'ZPSTRF'
!
!                 Use default tolerance
!
                           tol = -ONE
                           CALL ZPSTRF(uplo,n,Afac,lda,Piv,comprank,tol,&
     &                                 Rwork,info)
!
!                 Check error code from ZPSTRF.
!
                           IF ( (info<izero) .OR.                       &
     &                          (info/=izero .AND. rank==n) .OR.        &
     &                          (info<=izero .AND. rank<n) ) THEN
                              CALL ALAERH(path,'ZPSTRF',info,izero,uplo,&
     &                           n,n,-1,-1,nb,imat,nfail,nerrs,Nout)
                              CYCLE
                           ENDIF
!
!                 Skip the test if INFO is not 0.
!
                           IF ( info==0 ) THEN
!
!                 Reconstruct matrix from factors and compute residual.
!
!                 PERM holds permuted L*L^T or U^T*U
!
                              CALL ZPST01(uplo,n,A,lda,Afac,lda,Perm,   &
     &                           lda,Piv,Rwork,result,comprank)
!
!                 Print information about the tests that did not pass
!                 the threshold or where computed rank was not RANK.
!
                              IF ( n==0 ) comprank = 0
                              rankdiff = rank - comprank
                              IF ( result>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL ALAHD(Nout,path)
                                 WRITE (Nout,FMT=99001) uplo , n ,      &
     &                                  rank , rankdiff , nb , imat ,   &
     &                                  result
                                 nfail = nfail + 1
                              ENDIF
                              nrun = nrun + 1
                           ENDIF
                        ENDDO
!
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL ALASUM(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' UPLO = ''',A1,''', N =',I5,', RANK =',I3,', Diff =',I5, &
     &        ', NB =',I4,', type ',I2,', Ratio =',G12.5)
!
!     End of ZCHKPS
!
      END SUBROUTINE ZCHKPS
