!*==ilaenv2stage.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
      INTEGER FUNCTION ILAENV2STAGE(Ispec,Name,Opts,N1,N2,N3,N4)
      IMPLICIT NONE
!*--ILAENV2STAGE254
!     .. Scalar Arguments ..
      CHARACTER*(*) Name , Opts
      INTEGER Ispec , N1 , N2 , N3 , N4
!     ..
!
!  =====================================================================
!
!     .. Local variables ..
      INTEGER iispec
!     .. External Functions ..
      INTEGER IPARAM2STAGE
      EXTERNAL IPARAM2STAGE
!     ..
!     .. Arrays in Common ..
      INTEGER IPArms(100)
!     ..
!     .. Common blocks ..
      COMMON /CLAENV/ IPArms
!     ..
!     .. Save statement ..
      SAVE /CLAENV/ 
!     ..
!     .. Executable Statements ..
!
      IF ( (Ispec<1) .OR. (Ispec>5) ) THEN
!
!        Invalid value for ISPEC
!
         ILAENV2STAGE = -1
!
!     1 <= ISPEC <= 5: 2stage eigenvalues SVD routines.
!
      ELSEIF ( Ispec==1 ) THEN
         ILAENV2STAGE = IPArms(1)
      ELSE
         iispec = 16 + Ispec
         ILAENV2STAGE = IPARAM2STAGE(iispec,Name,Opts,N1,N2,N3,N4)
!
      ENDIF
!
      END FUNCTION ILAENV2STAGE
