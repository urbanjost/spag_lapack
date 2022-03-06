!*==dckgsv.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DCKGSV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCKGSV( NM, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH,
!                          NMAX, A, AF, B, BF, U, V, Q, ALPHA, BETA, R,
!                          IWORK, WORK, RWORK, NIN, NOUT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, NIN, NM, NMATS, NMAX, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * ),
!      $                   PVAL( * )
!       DOUBLE PRECISION   A( * ), AF( * ), ALPHA( * ), B( * ), BETA( * ),
!      $                   BF( * ), Q( * ), R( * ), RWORK( * ), U( * ),
!      $                   V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCKGSV tests DGGSVD:
!>        the GSVD for M-by-N matrix A and P-by-N matrix B.
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
!>          The values of the matrix row dimension M.
!> \endverbatim
!>
!> \param[in] PVAL
!> \verbatim
!>          PVAL is INTEGER array, dimension (NP)
!>          The values of the matrix row dimension P.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NN)
!>          The values of the matrix column dimension N.
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
!>          THRESH is DOUBLE PRECISION
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
!>          A is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (NMAX)
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
!>          > 0 :  If DLATMS returns an error code, the absolute value
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DCKGSV(Nm,Mval,Pval,Nval,Nmats,Iseed,Thresh,Nmax,A,Af, &
     &                  B,Bf,U,V,Q,Alpha,Beta,R,Iwork,Work,Rwork,Nin,   &
     &                  Nout,Info)
      IMPLICIT NONE
!*--DCKGSV202
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Nin , Nm , Nmats , Nmax , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4) , Iwork(*) , Mval(*) , Nval(*) , Pval(*)
      DOUBLE PRECISION A(*) , Af(*) , Alpha(*) , B(*) , Beta(*) ,       &
     &                 Bf(*) , Q(*) , R(*) , Rwork(*) , U(*) , V(*) ,   &
     &                 Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NTESTS
      PARAMETER (NTESTS=12)
      INTEGER NTYPES
      PARAMETER (NTYPES=8)
!     ..
!     .. Local Scalars ..
      LOGICAL firstt
      CHARACTER dista , distb , type
      CHARACTER*3 path
      INTEGER i , iinfo , im , imat , kla , klb , kua , kub , lda ,     &
     &        ldb , ldq , ldr , ldu , ldv , lwork , m , modea , modeb , &
     &        n , nfail , nrun , nt , p
      DOUBLE PRECISION anorm , bnorm , cndnma , cndnmb
!     ..
!     .. Local Arrays ..
      LOGICAL dotype(NTYPES)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAHDG , ALAREQ , ALASUM , DGSVTS3 , DLATB9 , DLATMS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:3) = 'GSV'
      Info = 0
      nrun = 0
      nfail = 0
      firstt = .TRUE.
      CALL ALAREQ(path,Nmats,dotype,NTYPES,Nin,Nout)
      lda = Nmax
      ldb = Nmax
      ldu = Nmax
      ldv = Nmax
      ldq = Nmax
      ldr = Nmax
      lwork = Nmax*Nmax
!
!     Do for each value of M in MVAL.
!
      DO im = 1 , Nm
         m = Mval(im)
         p = Pval(im)
         n = Nval(im)
!
         DO imat = 1 , NTYPES
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
            IF ( dotype(imat) ) THEN
!
!           Set up parameters with DLATB9 and generate test
!           matrices A and B with DLATMS.
!
               CALL DLATB9(path,imat,m,p,n,type,kla,kua,klb,kub,anorm,  &
     &                     bnorm,modea,modeb,cndnma,cndnmb,dista,distb)
!
!           Generate M by N matrix A
!
               CALL DLATMS(m,n,dista,Iseed,type,Rwork,modea,cndnma,     &
     &                     anorm,kla,kua,'No packing',A,lda,Work,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nout,FMT=99001) iinfo
                  Info = ABS(iinfo)
                  CYCLE
               ENDIF
!
               CALL DLATMS(p,n,distb,Iseed,type,Rwork,modeb,cndnmb,     &
     &                     bnorm,klb,kub,'No packing',B,ldb,Work,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nout,FMT=99001) iinfo
                  Info = ABS(iinfo)
                  CYCLE
               ENDIF
!
               nt = 6
!
               CALL DGSVTS3(m,p,n,A,Af,lda,B,Bf,ldb,U,ldu,V,ldv,Q,ldq,  &
     &                      Alpha,Beta,R,ldr,Iwork,Work,lwork,Rwork,    &
     &                      result)
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
                     WRITE (Nout,FMT=99002) m , p , n , imat , i ,      &
     &                      result(i)
                     nfail = nfail + 1
                  ENDIF
               ENDDO
               nrun = nrun + nt
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL ALASUM(path,Nout,nfail,nrun,0)
!
99001 FORMAT (' DLATMS in DCKGSV   INFO = ',I5)
99002 FORMAT (' M=',I4,' P=',I4,', N=',I4,', type ',I2,', test ',I2,    &
     &        ', ratio=',G13.6)
!
!     End of DCKGSV
!
      END SUBROUTINE DCKGSV
