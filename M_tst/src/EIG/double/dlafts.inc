!*==dlafts.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b dlafts
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dlafts( TYPE, M, N, IMAT, NTESTS, RESULT, ISEED,
!                          THRESH, IOUNIT, IE )
!
!       .. Scalar Arguments ..
!       CHARACTER*3        TYPE
!       INTEGER            IE, IMAT, IOUNIT, M, N, NTESTS
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   RESULT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    dlafts tests the result vector against the threshold value to
!>    see which tests for this matrix type failed to pass the threshold.
!>    Output is to the file given by unit IOUNIT.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  TYPE   - CHARACTER*3
!>           On entry, TYPE specifies the matrix type to be used in the
!>           printed messages.
!>           Not modified.
!>
!>  N      - INTEGER
!>           On entry, N specifies the order of the test matrix.
!>           Not modified.
!>
!>  IMAT   - INTEGER
!>           On entry, IMAT specifies the type of the test matrix.
!>           A listing of the different types is printed by dlahd2
!>           to the output file if a test fails to pass the threshold.
!>           Not modified.
!>
!>  NTESTS - INTEGER
!>           On entry, NTESTS is the number of tests performed on the
!>           subroutines in the path given by TYPE.
!>           Not modified.
!>
!>  RESULT - DOUBLE PRECISION               array of dimension( NTESTS )
!>           On entry, RESULT contains the test ratios from the tests
!>           performed in the calling program.
!>           Not modified.
!>
!>  ISEED  - INTEGER            array of dimension( 4 )
!>           Contains the random seed that generated the matrix used
!>           for the tests whose ratios are in RESULT.
!>           Not modified.
!>
!>  THRESH - DOUBLE PRECISION
!>           On entry, THRESH specifies the acceptable threshold of the
!>           test ratios.  If RESULT( K ) > THRESH, then the K-th test
!>           did not pass the threshold and a message will be printed.
!>           Not modified.
!>
!>  IOUNIT - INTEGER
!>           On entry, IOUNIT specifies the unit number of the file
!>           to which the messages are printed.
!>           Not modified.
!>
!>  IE     - INTEGER
!>           On entry, IE contains the number of tests which have
!>           failed to pass the threshold so far.
!>           Updated on exit if any of the ratios in RESULT also fail.
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
      SUBROUTINE dlafts(Type,M,N,Imat,Ntests,Result,Iseed,Thresh,Iounit,&
     &                  Ie)
      IMPLICIT NONE
!*--DLAFTS103
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*3 Type
      INTEGER Ie , Imat , Iounit , M , N , Ntests
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      DOUBLE PRECISION Result(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER k
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
      IF ( M==N ) THEN
!
!     Output for square matrices:
!
         DO k = 1 , Ntests
            IF ( Result(k)>=Thresh ) THEN
!
!           If this is the first test to fail, call dlahd2
!           to print a header to the data file.
!
               IF ( Ie==0 ) CALL dlahd2(Iounit,Type)
               Ie = Ie + 1
               IF ( Result(k)<10000.0D0 ) THEN
                  WRITE (Iounit,FMT=99001) N , Imat , Iseed , k ,       &
     &                   Result(k)
99001             FORMAT (' Matrix order=',I5,', type=',I2,', seed=',   &
     &                    4(I4,','),' result ',I3,' is',0P,F8.2)
               ELSE
                  WRITE (Iounit,FMT=99002) N , Imat , Iseed , k ,       &
     &                   Result(k)
99002             FORMAT (' Matrix order=',I5,', type=',I2,', seed=',   &
     &                    4(I4,','),' result ',I3,' is',1P,D10.3)
               ENDIF
            ENDIF
         ENDDO
      ELSE
!
!     Output for rectangular matrices
!
         DO k = 1 , Ntests
            IF ( Result(k)>=Thresh ) THEN
!
!              If this is the first test to fail, call dlahd2
!              to print a header to the data file.
!
               IF ( Ie==0 ) CALL dlahd2(Iounit,Type)
               Ie = Ie + 1
               IF ( Result(k)<10000.0D0 ) THEN
                  WRITE (Iounit,FMT=99003) M , N , Imat , Iseed , k ,   &
     &                   Result(k)
99003             FORMAT (1X,I5,' x',I5,' matrix, type=',I2,', s',      &
     &                    'eed=',3(I4,','),I4,': result ',I3,' is',0P,  &
     &                    F8.2)
               ELSE
                  WRITE (Iounit,FMT=99004) M , N , Imat , Iseed , k ,   &
     &                   Result(k)
99004             FORMAT (1X,I5,' x',I5,' matrix, type=',I2,', s',      &
     &                    'eed=',3(I4,','),I4,': result ',I3,' is',1P,  &
     &                    D10.3)
               ENDIF
            ENDIF
         ENDDO
!
      ENDIF
!
!     End of dlafts
!
      END SUBROUTINE dlafts
