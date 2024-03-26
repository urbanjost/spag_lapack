!*==cdrvsy_aa.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cdrvsy_aa
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cdrvsy_aa( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
!                             A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK,
!                             NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
!       REAL               RWORK( * )
!       COMPLEX            A( * ), AFAC( * ), AINV( * ), B( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cdrvsy_aa tests the driver routine CSYSV_AA.
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
!>          The maximum value permitted for N, used in dimensioning the
!>          work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (NMAX*max(2,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (NMAX+2*NRHS)
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
!> \date November 2017
!
!  @generated from LIN/ddrvsy_aa.f, fortran d -> c, Thu Nov 17 12:14:51 2016
!
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cdrvsy_aa(Dotype,Nn,Nval,Nrhs,Thresh,Tsterr,Nmax,A,    &
     &                     Afac,Ainv,B,X,Xact,Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--CDRVSY_AA158
!
!  -- LAPACK test routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nmax , Nn , Nout , Nrhs
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nval(*)
      REAL Rwork(*)
      COMPLEX A(*) , Afac(*) , Ainv(*) , B(*) , Work(*) , X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX CZERO
      PARAMETER (CZERO=0.0E+0)
      INTEGER NTYPES , NTESTS
      PARAMETER (NTYPES=10,NTESTS=3)
      INTEGER NFACT
      PARAMETER (NFACT=2)
!     ..
!     .. Local Scalars ..
      LOGICAL zerot
      CHARACTER dist , fact , type , uplo , xtype
      CHARACTER*3 matpath , path
      INTEGER i , i1 , i2 , ifact , imat , in , info , ioff , iuplo ,   &
     &        izero , j , k , kl , ku , lda , lwork , mode , n , nb ,   &
     &        nbmin , nerrs , nfail , nimat , nrun , nt
      REAL anorm , cndnum
!     ..
!     .. Local Arrays ..
      CHARACTER facts(NFACT) , uplos(2)
      INTEGER iseed(4) , iseedy(4)
      REAL result(NTESTS)
!     ..
!     .. External Functions ..
      REAL CLANSY
EXTERNAL    CLANSY
!     ..
!     .. External Subroutines ..
EXTERNAL    CLACPY
EXTERNAL    CLASET,   CSYSV_AA
EXTERNAL    CSYTRF_AA
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
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
      DATA uplos/'U' , 'L'/ , facts/'F' , 'N'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
!     Test path
!
      path(1:1) = 'Complex precision'
      path(2:3) = 'SA'
!
!     Path to generate matrices
!
      matpath(1:1) = 'Complex precision'
      matpath(2:3) = 'SY'
!
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL cerrvx(path,Nout)
      INFot = 0
!
!     Set the block size and minimum block size for testing.
!
      nb = 1
      nbmin = 2
      CALL xlaenv(1,nb)
      CALL xlaenv(2,nbmin)
!
!     Do for each value of N in NVAL
!
      DO in = 1 , Nn
         n = Nval(in)
         lwork = MAX(3*n-2,n*(1+nb))
         lwork = MAX(lwork,1)
         lda = MAX(n,1)
         xtype = 'N'
         nimat = NTYPES
         IF ( n<=0 ) nimat = 1
!
         DO imat = 1 , nimat
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
            IF ( Dotype(imat) ) THEN
!
!           Skip types 3, 4, 5, or 6 if the matrix size is too small.
!
               zerot = imat>=3 .AND. imat<=6
               IF ( .NOT.(zerot .AND. n<imat-2) ) THEN
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
                  DO iuplo = 1 , 2
                     uplo = uplos(iuplo)
!
!              Set up parameters with clatb4 and generate a test matrix
!              with CLATMS.
!
                     CALL clatb4(matpath,imat,n,n,type,kl,ku,anorm,mode,&
     &                           cndnum,dist)
!
                     SRNamt = 'CLATMS'
                     CALL CLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,uplo,A,lda,Work,info)
!
!              Check error code from CLATMS.
!
                     IF ( info/=0 ) THEN
                        CALL alaerh(path,'CLATMS',info,0,uplo,n,n,-1,-1,&
     &                              -1,imat,nfail,nerrs,Nout)
                        CYCLE
                     ENDIF
!
!              For types 3-6, zero one or more rows and columns of the
!              matrix to test that INFO is returned correctly.
!
                     IF ( zerot ) THEN
                        IF ( imat==3 ) THEN
                           izero = 1
                        ELSEIF ( imat==4 ) THEN
                           izero = n
                        ELSE
                           izero = n/2 + 1
                        ENDIF
!
                        IF ( imat>=6 ) THEN
                           ioff = 0
                           IF ( iuplo==1 ) THEN
!
!                       Set the first IZERO rows and columns to zero.
!
                              DO j = 1 , n
                                 i2 = MIN(j,izero)
                                 DO i = 1 , i2
                                    A(ioff+i) = CZERO
                                 ENDDO
                                 ioff = ioff + lda
                              ENDDO
                              izero = 1
                           ELSE
!
!                       Set the last IZERO rows and columns to zero.
!
                              DO j = 1 , n
                                 i1 = MAX(j,izero)
                                 DO i = i1 , n
                                    A(ioff+i) = CZERO
                                 ENDDO
                                 ioff = ioff + lda
                              ENDDO
                           ENDIF
!
!                    Set row and column IZERO to zero.
!
                        ELSEIF ( iuplo==1 ) THEN
                           ioff = (izero-1)*lda
                           DO i = 1 , izero - 1
                              A(ioff+i) = CZERO
                           ENDDO
                           ioff = ioff + izero
                           DO i = izero , n
                              A(ioff) = CZERO
                              ioff = ioff + lda
                           ENDDO
                        ELSE
                           ioff = izero
                           DO i = 1 , izero - 1
                              A(ioff) = CZERO
                              ioff = ioff + lda
                           ENDDO
                           ioff = ioff - izero
                           DO i = izero , n
                              A(ioff+i) = CZERO
                           ENDDO
                        ENDIF
                     ELSE
                        izero = 0
                     ENDIF
!
                     DO ifact = 1 , NFACT
!
!                 Do first for FACT = 'F', then for other values.
!
                        fact = facts(ifact)
!
!                 Form an exact solution and set the right hand side.
!
                        SRNamt = 'clarhs'
                        CALL clarhs(matpath,xtype,uplo,' ',n,n,kl,ku,   &
     &                              Nrhs,A,lda,Xact,lda,B,lda,iseed,    &
     &                              info)
                        xtype = 'C'
!
!                 --- Test CSYSV_AA  ---
!
                        IF ( ifact==2 ) THEN
                           CALL CLACPY(uplo,n,n,A,lda,Afac,lda)
                           CALL CLACPY('Full',n,Nrhs,B,lda,X,lda)
!
!                    Factor the matrix and solve the system using CSYSV_AA.
!
                           SRNamt = 'CSYSV_AA'
                           CALL CSYSV_AA(uplo,n,Nrhs,Afac,lda,Iwork,X,  &
     &                        lda,Work,lwork,info)
!
!                    Adjust the expected value of INFO to account for
!                    pivoting.
!
                           IF ( izero>0 ) THEN
                              j = 1
                              k = izero
                              DO
                                 IF ( j==k ) THEN
                                    k = Iwork(j)
                                 ELSEIF ( Iwork(j)==k ) THEN
                                    k = j
                                 ENDIF
                                 IF ( j<k ) THEN
                                    j = j + 1
                                    CYCLE
                                 ENDIF
                                 EXIT
                              ENDDO
                           ELSE
                              k = 0
                           ENDIF
!
!                    Check error code from CSYSV_AA .
!
                           IF ( info/=k ) THEN
                              CALL alaerh(path,'CSYSV_AA ',info,k,uplo, &
     &                           n,n,-1,-1,Nrhs,imat,nfail,nerrs,Nout)
                              CYCLE
                           ELSEIF ( info/=0 ) THEN
                              CYCLE
                           ENDIF
!
!                    Reconstruct matrix from factors and compute
!                    residual.
!
                           CALL csyt01_aa(uplo,n,A,lda,Afac,lda,Iwork,  &
     &                        Ainv,lda,Rwork,result(1))
!
!                    Compute residual of the computed solution.
!
                           CALL CLACPY('Full',n,Nrhs,B,lda,Work,lda)
                           CALL csyt02(uplo,n,Nrhs,A,lda,X,lda,Work,lda,&
     &                                 Rwork,result(2))
                           nt = 2
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                           DO k = 1 , nt
                              IF ( result(k)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL aladhd(Nout,path)
                                 WRITE (Nout,FMT=99001) 'CSYSV_AA ' ,   &
     &                                  uplo , n , imat , k , result(k)
                                 nfail = nfail + 1
                              ENDIF
                           ENDDO
                           nrun = nrun + nt
                        ENDIF
!
                     ENDDO
!
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasvm(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (1X,A,', UPLO=''',A1,''', N =',I5,', type ',I2,', test ',  &
     &        I2,', ratio =',G12.5)
!
!     End of cdrvsy_aa
!
      END SUBROUTINE cdrvsy_aa
