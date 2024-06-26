!*==zckcsd.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b zckcsd
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zckcsd( NM, MVAL, PVAL, QVAL, NMATS, ISEED, THRESH,
!                          MMAX, X, XF, U1, U2, V1T, V2T, THETA, IWORK,
!                          WORK, RWORK, NIN, NOUT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, NIN, NM, NMATS, MMAX, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 ), IWORK( * ), MVAL( * ), PVAL( * ),
!      $                   QVAL( * )
!       DOUBLE PRECISION   RWORK( * ), THETA( * )
!       COMPLEX*16         U1( * ), U2( * ), V1T( * ), V2T( * ),
!      $                   WORK( * ), X( * ), XF( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zckcsd tests ZUNCSD:
!>        the CSD for an M-by-M unitary matrix X partitioned as
!>        [ X11 X12; X21 X22 ]. X11 is P-by-Q.
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
!>          PVAL is INTEGER array, dimension (NM)
!>          The values of the matrix row dimension P.
!> \endverbatim
!>
!> \param[in] QVAL
!> \verbatim
!>          QVAL is INTEGER array, dimension (NM)
!>          The values of the matrix column dimension Q.
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
!> \param[in] MMAX
!> \verbatim
!>          MMAX is INTEGER
!>          The maximum value permitted for M, used in dimensioning the
!>          work arrays.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (MMAX*MMAX)
!> \endverbatim
!>
!> \param[out] XF
!> \verbatim
!>          XF is COMPLEX*16 array, dimension (MMAX*MMAX)
!> \endverbatim
!>
!> \param[out] U1
!> \verbatim
!>          U1 is COMPLEX*16 array, dimension (MMAX*MMAX)
!> \endverbatim
!>
!> \param[out] U2
!> \verbatim
!>          U2 is COMPLEX*16 array, dimension (MMAX*MMAX)
!> \endverbatim
!>
!> \param[out] V1T
!> \verbatim
!>          V1T is COMPLEX*16 array, dimension (MMAX*MMAX)
!> \endverbatim
!>
!> \param[out] V2T
!> \verbatim
!>          V2T is COMPLEX*16 array, dimension (MMAX*MMAX)
!> \endverbatim
!>
!> \param[out] THETA
!> \verbatim
!>          THETA is DOUBLE PRECISION array, dimension (MMAX)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array
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
!>          > 0 :  If zlaror returns an error code, the absolute value
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE zckcsd(Nm,Mval,Pval,Qval,Nmats,Iseed,Thresh,Mmax,X,Xf, &
     &                  U1,U2,V1t,V2t,Theta,Iwork,Work,Rwork,Nin,Nout,  &
     &                  Info)
      IMPLICIT NONE
!*--ZCKCSD188
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Nin , Nm , Nmats , Mmax , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4) , Iwork(*) , Mval(*) , Pval(*) , Qval(*)
      DOUBLE PRECISION Rwork(*) , Theta(*)
      COMPLEX*16 U1(*) , U2(*) , V1t(*) , V2t(*) , Work(*) , X(*) ,     &
     &           Xf(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NTESTS
      PARAMETER (NTESTS=15)
      INTEGER NTYPES
      PARAMETER (NTYPES=4)
      DOUBLE PRECISION GAPDIGIT , ORTH , REALONE , REALZERO , TEN
      PARAMETER (GAPDIGIT=18.0D0,ORTH=1.0D-12,REALONE=1.0D0,            &
     &           REALZERO=0.0D0,TEN=10.0D0)
      COMPLEX*16 ONE , ZERO
      PARAMETER (ONE=(1.0D0,0.0D0),ZERO=(0.0D0,0.0D0))
      DOUBLE PRECISION PIOVER2
      PARAMETER (PIOVER2=1.57079632679489661923132169163975144210D0)
!     ..
!     .. Local Scalars ..
      LOGICAL firstt
      CHARACTER*3 path
      INTEGER i , iinfo , im , imat , j , ldu1 , ldu2 , ldv1t , ldv2t , &
     &        ldx , lwork , m , nfail , nrun , nt , p , q , r
!     ..
!     .. Local Arrays ..
      LOGICAL dotype(NTYPES)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Subroutines ..
      EXTERNAL ZLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MIN
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:3) = 'CSD'
      Info = 0
      nrun = 0
      nfail = 0
      firstt = .TRUE.
      CALL alareq(path,Nmats,dotype,NTYPES,Nin,Nout)
      ldx = Mmax
      ldu1 = Mmax
      ldu2 = Mmax
      ldv1t = Mmax
      ldv2t = Mmax
      lwork = Mmax*Mmax
!
!     Do for each value of M in MVAL.
!
      DO im = 1 , Nm
         m = Mval(im)
         p = Pval(im)
         q = Qval(im)
!
         DO imat = 1 , NTYPES
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
            IF ( dotype(imat) ) THEN
!
!           Generate X
!
               IF ( imat==1 ) THEN
                  CALL zlaror('L','I',m,m,X,ldx,Iseed,Work,iinfo)
                  IF ( m/=0 .AND. iinfo/=0 ) THEN
                     WRITE (Nout,FMT=99001) m , iinfo
                     Info = ABS(iinfo)
                     CYCLE
                  ENDIF
               ELSEIF ( imat==2 ) THEN
                  r = MIN(p,m-p,q,m-q)
                  DO i = 1 , r
                     Theta(i) = PIOVER2*dlarnd(1,Iseed)
                  ENDDO
                  CALL ZLACSG(m,p,q,Theta,Iseed,X,ldx,Work)
                  DO i = 1 , m
                     DO j = 1 , m
                        X(i+(j-1)*ldx) = X(i+(j-1)*ldx)                 &
     &                     + ORTH*dlarnd(2,Iseed)
                     ENDDO
                  ENDDO
               ELSEIF ( imat==3 ) THEN
                  r = MIN(p,m-p,q,m-q)
                  DO i = 1 , r + 1
                     Theta(i) = TEN**(-dlarnd(1,Iseed)*GAPDIGIT)
                  ENDDO
                  DO i = 2 , r + 1
                     Theta(i) = Theta(i-1) + Theta(i)
                  ENDDO
                  DO i = 1 , r
                     Theta(i) = PIOVER2*Theta(i)/Theta(r+1)
                  ENDDO
                  CALL ZLACSG(m,p,q,Theta,Iseed,X,ldx,Work)
               ELSE
                  CALL ZLASET('F',m,m,ZERO,ONE,X,ldx)
                  DO i = 1 , m
                     j = INT(dlaran(Iseed)*m) + 1
                     IF ( j/=i )                                        &
     &                    CALL ZDROT(m,X(1+(i-1)*ldx),1,X(1+(j-1)*ldx), &
     &                    1,REALZERO,REALONE)
                  ENDDO
               ENDIF
!
               nt = 15
!
               CALL zcsdts(m,p,q,X,Xf,ldx,U1,ldu1,U2,ldu2,V1t,ldv1t,V2t,&
     &                     ldv2t,Theta,Iwork,Work,lwork,Rwork,result)
!
!           Print information about the tests that did not
!           pass the threshold.
!
               DO i = 1 , nt
                  IF ( result(i)>=Thresh ) THEN
                     IF ( nfail==0 .AND. firstt ) THEN
                        firstt = .FALSE.
                        CALL alahdg(Nout,path)
                     ENDIF
                     WRITE (Nout,FMT=99002) m , p , q , imat , i ,      &
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
      CALL alasum(path,Nout,nfail,nrun,0)
!
99001 FORMAT (' zlaror in zckcsd: M = ',I5,', INFO = ',I15)
99002 FORMAT (' M=',I4,' P=',I4,', Q=',I4,', type ',I2,', test ',I2,    &
     &        ', ratio=',G13.6)
!
!     End of zckcsd
!
      END SUBROUTINE zckcsd
