!*==ccklse.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CCKLSE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCKLSE( NN, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH,
!                          NMAX, A, AF, B, BF, X, WORK, RWORK, NIN, NOUT,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, NIN, NMATS, NMAX, NN, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * )
!       REAL               RWORK( * )
!       COMPLEX            A( * ), AF( * ), B( * ), BF( * ), WORK( * ),
!      $                   X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCKLSE tests CGGLSE - a subroutine for solving linear equality
!> constrained least square problem (LSE).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER
!>          The number of values of (M,P,N) contained in the vectors
!>          (MVAL, PVAL, NVAL).
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NN)
!>          The values of the matrix row(column) dimension M.
!> \endverbatim
!>
!> \param[in] PVAL
!> \verbatim
!>          PVAL is INTEGER array, dimension (NN)
!>          The values of the matrix row(column) dimension P.
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
!>          A is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension (5*NMAX)
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
!>          > 0 :  If CLATMS returns an error code, the absolute value
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE CCKLSE(Nn,Mval,Pval,Nval,Nmats,Iseed,Thresh,Nmax,A,Af, &
     &                  B,Bf,X,Work,Rwork,Nin,Nout,Info)
      IMPLICIT NONE
!*--CCKLSE171
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Nin , Nmats , Nmax , Nn , Nout
      REAL Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4) , Mval(*) , Nval(*) , Pval(*)
      REAL Rwork(*)
      COMPLEX A(*) , Af(*) , B(*) , Bf(*) , Work(*) , X(*)
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
      INTEGER i , iinfo , ik , imat , kla , klb , kua , kub , lda ,     &
     &        ldb , lwork , m , modea , modeb , n , nfail , nrun , nt , &
     &        p
      REAL anorm , bnorm , cndnma , cndnmb
!     ..
!     .. Local Arrays ..
      LOGICAL dotype(NTYPES)
      REAL result(NTESTS)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAHDG , ALAREQ , ALASUM , CLARHS , CLATMS , CLSETS ,    &
     &         SLATB9
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:3) = 'LSE'
      Info = 0
      nrun = 0
      nfail = 0
      firstt = .TRUE.
      CALL ALAREQ(path,Nmats,dotype,NTYPES,Nin,Nout)
      lda = Nmax
      ldb = Nmax
      lwork = Nmax*Nmax
!
!     Check for valid input values.
!
      DO ik = 1 , Nn
         m = Mval(ik)
         p = Pval(ik)
         n = Nval(ik)
         IF ( p>n .OR. n>m+p ) THEN
            IF ( firstt ) THEN
               WRITE (Nout,FMT=*)
               firstt = .FALSE.
            ENDIF
            WRITE (Nout,FMT=99003) m , p , n
         ENDIF
      ENDDO
      firstt = .TRUE.
!
!     Do for each value of M in MVAL.
!
      DO ik = 1 , Nn
         m = Mval(ik)
         p = Pval(ik)
         n = Nval(ik)
         IF ( p<=n .AND. n<=m+p ) THEN
!
            DO imat = 1 , NTYPES
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
               IF ( dotype(imat) ) THEN
!
!           Set up parameters with SLATB9 and generate test
!           matrices A and B with CLATMS.
!
                  CALL SLATB9(path,imat,m,p,n,type,kla,kua,klb,kub,     &
     &                        anorm,bnorm,modea,modeb,cndnma,cndnmb,    &
     &                        dista,distb)
!
                  CALL CLATMS(m,n,dista,Iseed,type,Rwork,modea,cndnma,  &
     &                        anorm,kla,kua,'No packing',A,lda,Work,    &
     &                        iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nout,FMT=99001) iinfo
                     Info = ABS(iinfo)
                     CYCLE
                  ENDIF
!
                  CALL CLATMS(p,n,distb,Iseed,type,Rwork,modeb,cndnmb,  &
     &                        bnorm,klb,kub,'No packing',B,ldb,Work,    &
     &                        iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nout,FMT=99001) iinfo
                     Info = ABS(iinfo)
                     CYCLE
                  ENDIF
!
!           Generate the right-hand sides C and D for the LSE.
!
                  CALL CLARHS('CGE','New solution','Upper','N',m,n,     &
     &                        MAX(m-1,0),MAX(n-1,0),1,A,lda,X(4*Nmax+1),&
     &                        MAX(n,1),X,MAX(m,1),Iseed,iinfo)
!
                  CALL CLARHS('CGE','Computed','Upper','N',p,n,         &
     &                        MAX(p-1,0),MAX(n-1,0),1,B,ldb,X(4*Nmax+1),&
     &                        MAX(n,1),X(2*Nmax+1),MAX(p,1),Iseed,iinfo)
!
                  nt = 2
!
                  CALL CLSETS(m,p,n,A,Af,lda,B,Bf,ldb,X,X(Nmax+1),      &
     &                        X(2*Nmax+1),X(3*Nmax+1),X(4*Nmax+1),Work, &
     &                        lwork,Rwork,result(1))
!
!           Print information about the tests that did not
!           pass the threshold.
!
                  DO i = 1 , nt
                     IF ( result(i)>=Thresh ) THEN
                        IF ( nfail==0 .AND. firstt ) THEN
                           firstt = .FALSE.
                           CALL ALAHDG(Nout,path)
                        ENDIF
                        WRITE (Nout,FMT=99002) m , p , n , imat , i ,   &
     &                         result(i)
                        nfail = nfail + 1
                     ENDIF
                  ENDDO
                  nrun = nrun + nt
               ENDIF
!
            ENDDO
         ENDIF
      ENDDO
!
!     Print a summary of the results.
!
      CALL ALASUM(path,Nout,nfail,nrun,0)
!
99001 FORMAT (' CLATMS in CCKLSE   INFO = ',I5)
99002 FORMAT (' M=',I4,' P=',I4,', N=',I4,', type ',I2,', test ',I2,    &
     &        ', ratio=',G13.6)
99003 FORMAT (' *** Invalid input  for LSE:  M = ',I6,', P = ',I6,      &
     &        ', N = ',I6,';',/'     must satisfy P <= N <= P+M  ',     &
     &        '(this set of values will be skipped)')
!
!     End of CCKLSE
!
      END SUBROUTINE CCKLSE
