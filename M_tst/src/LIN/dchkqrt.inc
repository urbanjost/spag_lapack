!*==dchkqrt.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b dchkqrt
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dchkqrt( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB,
!                           NBVAL, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NM, NN, NNB, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            MVAL( * ), NBVAL( * ), NVAL( * )
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dchkqrt tests DGEQRT and DGEMQRT.
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!>          The number of values of NB contained in the vector NBVAL.
!> \endverbatim
!>
!> \param[in] NBVAL
!> \verbatim
!>          NBVAL is INTEGER array, dimension (NBVAL)
!>          The values of the blocksize NB.
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE dchkqrt(Thresh,Tsterr,Nm,Mval,Nn,Nval,Nnb,Nbval,Nout)
      IMPLICIT NONE
!*--DCHKQRT105
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nm , Nn , Nnb , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Mval(*) , Nbval(*) , Nval(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NTESTS
      PARAMETER (NTESTS=6)
!     ..
!     .. Local Scalars ..
      CHARACTER*3 path
      INTEGER i , j , k , t , m , n , nb , nfail , nerrs , nrun , minmn
!
!     .. Local Arrays ..
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Subroutines ..
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
!     .. Executable Statements ..
!
!     Initialize constants
!
      path(1:1) = 'D'
      path(2:3) = 'QT'
      nrun = 0
      nfail = 0
      nerrs = 0
!
!     Test the error exits
!
      IF ( Tsterr ) CALL derrqrt(path,Nout)
      INFot = 0
!
!     Do for each value of M in MVAL.
!
      DO i = 1 , Nm
         m = Mval(i)
!
!        Do for each value of N in NVAL.
!
         DO j = 1 , Nn
            n = Nval(j)
!
!        Do for each possible value of NB
!
            minmn = MIN(m,n)
            DO k = 1 , Nnb
               nb = Nbval(k)
!
!              Test DGEQRT and DGEMQRT
!
               IF ( (nb<=minmn) .AND. (nb>0) ) THEN
                  CALL dqrt04(m,n,nb,result)
!
!                 Print information about the tests that did not
!                 pass the threshold.
!
                  DO t = 1 , NTESTS
                     IF ( result(t)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL alahd(Nout,path)
                        WRITE (Nout,FMT=99001) m , n , nb , t ,         &
     &                         result(t)
                        nfail = nfail + 1
                     ENDIF
                  ENDDO
                  nrun = nrun + NTESTS
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasum(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' M=',I5,', N=',I5,', NB=',I4,' test(',I2,')=',G12.5)
!
!     End of dchkqrt
!
      END SUBROUTINE dchkqrt
