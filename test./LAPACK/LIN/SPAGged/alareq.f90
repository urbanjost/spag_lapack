!*==alareq.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ALAREQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
!
!       .. Scalar Arguments ..
!       CHARACTER*3        PATH
!       INTEGER            NIN, NMATS, NOUT, NTYPES
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ALAREQ handles input for the LAPACK test program.  It is called
!> to evaluate the input line which requested NMATS matrix types for
!> PATH.  The flow of control is as follows:
!>
!> If NMATS = NTYPES then
!>    DOTYPE(1:NTYPES) = .TRUE.
!> else
!>    Read the next input line for NMATS matrix types
!>    Set DOTYPE(I) = .TRUE. for each valid type I
!> endif
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3
!>          An LAPACK path name for testing.
!> \endverbatim
!>
!> \param[in] NMATS
!> \verbatim
!>          NMATS is INTEGER
!>          The number of matrix types to be used in testing this path.
!> \endverbatim
!>
!> \param[out] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          The vector of flags indicating if each type will be tested.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The maximum number of matrix types for this path.
!> \endverbatim
!>
!> \param[in] NIN
!> \verbatim
!>          NIN is INTEGER
!>          The unit number for input.  NIN >= 1.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number for output.  NOUT >= 1.
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
!> \ingroup aux_lin
!
!  =====================================================================
      SUBROUTINE ALAREQ(Path,Nmats,Dotype,Ntypes,Nin,Nout)
      IMPLICIT NONE
!*--ALAREQ94
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*3 Path
      INTEGER Nin , Nmats , Nout , Ntypes
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL firstt
      CHARACTER c1
      CHARACTER*10 intstr
      CHARACTER*80 line
      INTEGER i , i1 , ic , j , k , lenp , nt
!     ..
!     .. Local Arrays ..
      INTEGER nreq(100)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC LEN
!     ..
!     .. Data statements ..
      DATA intstr/'0123456789'/
!     ..
!     .. Executable Statements ..
!
      IF ( Nmats>=Ntypes ) THEN
!
!        Test everything if NMATS >= NTYPES.
!
         DO i = 1 , Ntypes
            Dotype(i) = .TRUE.
         ENDDO
      ELSE
         DO i = 1 , Ntypes
            Dotype(i) = .FALSE.
         ENDDO
         firstt = .TRUE.
!
!        Read a line of matrix types if 0 < NMATS < NTYPES.
!
         IF ( Nmats>0 ) THEN
            READ (Nin,FMT='(A80)',END=200) line
            lenp = LEN(line)
            i = 0
            DO j = 1 , Nmats
               nreq(j) = 0
               i1 = 0
               DO
                  i = i + 1
                  IF ( i>lenp ) THEN
                     IF ( j==Nmats .AND. i1>0 ) EXIT
                     WRITE (Nout,FMT=99005) line
                     WRITE (Nout,FMT=99006) Nmats
                     GOTO 100
                  ENDIF
                  IF ( line(i:i)/=' ' .AND. line(i:i)/=',' ) THEN
                     i1 = i
                     c1 = line(i1:i1)
!
!              Check that a valid integer was read
!
                     DO k = 1 , 10
                        IF ( c1==intstr(k:k) ) THEN
                           ic = k - 1
                           GOTO 2
                        ENDIF
                     ENDDO
                     WRITE (Nout,FMT=99004) i , line
                     WRITE (Nout,FMT=99006) Nmats
                     GOTO 100
 2                   nreq(j) = 10*nreq(j) + ic
                  ELSEIF ( i1>0 ) THEN
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         DO i = 1 , Nmats
            nt = nreq(i)
            IF ( nt>0 .AND. nt<=Ntypes ) THEN
               IF ( Dotype(nt) ) THEN
                  IF ( firstt ) WRITE (Nout,FMT=*)
                  firstt = .FALSE.
                  WRITE (Nout,FMT=99003) nt , Path
               ENDIF
               Dotype(nt) = .TRUE.
            ELSE
               WRITE (Nout,FMT=99001) Path , nt , Ntypes
99001          FORMAT (' *** Invalid type request for ',A3,', type  ',  &
     &                 I4,': must satisfy  1 <= type <= ',I2)
            ENDIF
         ENDDO
      ENDIF
 100  RETURN
!
 200  WRITE (Nout,FMT=99002) Path
99002 FORMAT (/' *** End of file reached when trying to read matrix ',  &
     &        'types for ',A3,/' *** Check that you are requesting the',&
     &        ' right number of types for each path',/)
99003 FORMAT (' *** Warning:  duplicate request of matrix type ',I2,    &
     &        ' for ',A3)
99004 FORMAT (//' *** Invalid integer value in column ',I2,' of input', &
     &        ' line:',/A79)
99005 FORMAT (//' *** Not enough matrix types on input line',/A79)
      WRITE (Nout,FMT=*)
      STOP
99006 FORMAT (' ==> Specify ',I4,' matrix types on this line or ',      &
     &        'adjust NTYPES on previous line')
!
!     End of ALAREQ
!
      END SUBROUTINE ALAREQ
