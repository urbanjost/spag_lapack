!*==sckgqr.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b sckgqr
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sckgqr( NM, MVAL, NP, PVAL, NN, NVAL, NMATS, ISEED,
!                          THRESH, NMAX, A, AF, AQ, AR, TAUA, B, BF, BZ,
!                          BT, BWK, TAUB, WORK, RWORK, NIN, NOUT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, NIN, NM, NMATS, NMAX, NN, NOUT, NP
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * )
!       REAL               A( * ), AF( * ), AQ( * ), AR( * ), B( * ),
!      $                   BF( * ), BT( * ), BWK( * ), BZ( * ),
!      $                   RWORK( * ), TAUA( * ), TAUB( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sckgqr tests
!> SGGQRF: GQR factorization for N-by-M matrix A and N-by-P matrix B,
!> SGGRQF: GRQ factorization for M-by-N matrix A and P-by-N matrix B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NM
!> \verbatim
!>          NM is INTEGER
!>          The number of values of M contained in the vector MVAL.
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NM)
!>          The values of the matrix row(column) dimension M.
!> \endverbatim
!>
!> \param[in] NP
!> \verbatim
!>          NP is INTEGER
!>          The number of values of P contained in the vector PVAL.
!> \endverbatim
!>
!> \param[in] PVAL
!> \verbatim
!>          PVAL is INTEGER array, dimension (NP)
!>          The values of the matrix row(column) dimension P.
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
!>          The values of the matrix column(row) dimension N.
!> \endverbatim
!>
!> \param[in] NMATS
!> \verbatim
!>          NMATS is INTEGER
!>          The number of matrix types to be tested for each combination
!>          of matrix dimensions.  If NMATS >= NTYPES (the maximum
!>          number of matrix types), then all the different types are
!>          generated for testing.  If NMATS < NTYPES, another input line
!>          is read to get the numbers of the matrix types to be used.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry, the seed of the random number generator.  The array
!>          elements should be between 0 and 4095, otherwise they will be
!>          reduced mod 4096, and ISEED(4) must be odd.
!>          On exit, the next seed in the random number sequence after
!>          all the test matrices have been generated.
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
!> \param[in] NMAX
!> \verbatim
!>          NMAX is INTEGER
!>          The maximum value permitted for M or N, used in dimensioning
!>          the work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AQ
!> \verbatim
!>          AQ is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AR
!> \verbatim
!>          AR is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] TAUA
!> \verbatim
!>          TAUA is REAL array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] BZ
!> \verbatim
!>          BZ is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] BT
!> \verbatim
!>          BT is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] BWK
!> \verbatim
!>          BWK is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] TAUB
!> \verbatim
!>          TAUB is REAL array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (NMAX)
!> \endverbatim
!>
!> \param[in] NIN
!> \verbatim
!>          NIN is INTEGER
!>          The unit number for input.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number for output.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0 :  successful exit
!>          > 0 :  If slatms returns an error code, the absolute value
!>                 of it is returned.
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE sckgqr(Nm,Mval,Np,Pval,Nn,Nval,Nmats,Iseed,Thresh,Nmax,&
     &                  A,Af,Aq,Ar,Taua,B,Bf,Bz,Bt,Bwk,Taub,Work,Rwork, &
     &                  Nin,Nout,Info)
      IMPLICIT NONE
!*--SCKGQR214
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Nin , Nm , Nmats , Nmax , Nn , Nout , Np
      REAL Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4) , Mval(*) , Nval(*) , Pval(*)
      REAL A(*) , Af(*) , Aq(*) , Ar(*) , B(*) , Bf(*) , Bt(*) , Bwk(*) &
     &     , Bz(*) , Rwork(*) , Taua(*) , Taub(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NTESTS
      PARAMETER (NTESTS=7)
      INTEGER NTYPES
      PARAMETER (NTYPES=8)
!     ..
!     .. Local Scalars ..
      LOGICAL firstt
      CHARACTER dista , distb , type
      CHARACTER*3 path
      INTEGER i , iinfo , im , imat , in , ip , kla , klb , kua , kub , &
     &        lda , ldb , lwork , m , modea , modeb , n , nfail , nrun ,&
     &        nt , p
      REAL anorm , bnorm , cndnma , cndnmb
!     ..
!     .. Local Arrays ..
      LOGICAL dotype(NTYPES)
      REAL result(NTESTS)
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
!     .. Executable Statements ..
!
!     Initialize constants.
!
      path(1:3) = 'GQR'
      Info = 0
      nrun = 0
      nfail = 0
      firstt = .TRUE.
      CALL alareq(path,Nmats,dotype,NTYPES,Nin,Nout)
      lda = Nmax
      ldb = Nmax
      lwork = Nmax*Nmax
!
!     Do for each value of M in MVAL.
!
      DO im = 1 , Nm
         m = Mval(im)
!
!        Do for each value of P in PVAL.
!
         DO ip = 1 , Np
            p = Pval(ip)
!
!           Do for each value of N in NVAL.
!
            DO in = 1 , Nn
               n = Nval(in)
!
               DO imat = 1 , NTYPES
!
!                 Do the tests only if DOTYPE( IMAT ) is true.
!
                  IF ( dotype(imat) ) THEN
!
!                 Test SGGRQF
!
!                 Set up parameters with slatb9 and generate test
!                 matrices A and B with slatms.
!
                     CALL slatb9('GRQ',imat,m,p,n,type,kla,kua,klb,kub, &
     &                           anorm,bnorm,modea,modeb,cndnma,cndnmb, &
     &                           dista,distb)
!
!                 Generate M by N matrix A
!
                     CALL slatms(m,n,dista,Iseed,type,Rwork,modea,      &
     &                           cndnma,anorm,kla,kua,'No packing',A,   &
     &                           lda,Work,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nout,FMT=99001) iinfo
                        Info = ABS(iinfo)
                        CYCLE
                     ENDIF
!
!                 Generate P by N matrix B
!
                     CALL slatms(p,n,distb,Iseed,type,Rwork,modeb,      &
     &                           cndnmb,bnorm,klb,kub,'No packing',B,   &
     &                           ldb,Work,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nout,FMT=99001) iinfo
                        Info = ABS(iinfo)
                        CYCLE
                     ENDIF
!
                     nt = 4
!
                     CALL sgrqts(m,p,n,A,Af,Aq,Ar,lda,Taua,B,Bf,Bz,Bt,  &
     &                           Bwk,ldb,Taub,Work,lwork,Rwork,result)
!
!                 Print information about the tests that did not
!                 pass the threshold.
!
                     DO i = 1 , nt
                        IF ( result(i)>=Thresh ) THEN
                           IF ( nfail==0 .AND. firstt ) THEN
                              firstt = .FALSE.
                              CALL alahdg(Nout,'GRQ')
                           ENDIF
                           WRITE (Nout,FMT=99002) m , p , n , imat , i ,&
     &                            result(i)
                           nfail = nfail + 1
                        ENDIF
                     ENDDO
                     nrun = nrun + nt
!
!                 Test SGGQRF
!
!                 Set up parameters with slatb9 and generate test
!                 matrices A and B with slatms.
!
                     CALL slatb9('GQR',imat,m,p,n,type,kla,kua,klb,kub, &
     &                           anorm,bnorm,modea,modeb,cndnma,cndnmb, &
     &                           dista,distb)
!
!                 Generate N-by-M matrix  A
!
                     CALL slatms(n,m,dista,Iseed,type,Rwork,modea,      &
     &                           cndnma,anorm,kla,kua,'No packing',A,   &
     &                           lda,Work,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nout,FMT=99001) iinfo
                        Info = ABS(iinfo)
                        CYCLE
                     ENDIF
!
!                 Generate N-by-P matrix  B
!
                     CALL slatms(n,p,distb,Iseed,type,Rwork,modea,      &
     &                           cndnma,bnorm,klb,kub,'No packing',B,   &
     &                           ldb,Work,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nout,FMT=99001) iinfo
                        Info = ABS(iinfo)
                        CYCLE
                     ENDIF
!
                     nt = 4
!
                     CALL sgqrts(n,m,p,A,Af,Aq,Ar,lda,Taua,B,Bf,Bz,Bt,  &
     &                           Bwk,ldb,Taub,Work,lwork,Rwork,result)
!
!                 Print information about the tests that did not
!                 pass the threshold.
!
                     DO i = 1 , nt
                        IF ( result(i)>=Thresh ) THEN
                           IF ( nfail==0 .AND. firstt ) THEN
                              firstt = .FALSE.
                              CALL alahdg(Nout,path)
                           ENDIF
                           WRITE (Nout,FMT=99003) n , m , p , imat , i ,&
     &                            result(i)
                           nfail = nfail + 1
                        ENDIF
                     ENDDO
                     nrun = nrun + nt
                  ENDIF
!
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasum(path,Nout,nfail,nrun,0)
!
99001 FORMAT (' slatms in sckgqr:    INFO = ',I5)
99002 FORMAT (' M=',I4,' P=',I4,', N=',I4,', type ',I2,', test ',I2,    &
     &        ', ratio=',G13.6)
99003 FORMAT (' N=',I4,' M=',I4,', P=',I4,', type ',I2,', test ',I2,    &
     &        ', ratio=',G13.6)
!
!     End of sckgqr
!
      END SUBROUTINE sckgqr
