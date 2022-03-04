!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPOTRF
   INTERFACE
      SUBROUTINE CPOTRF(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPOTRF
   END INTERFACE
END MODULE S_CPOTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPOTRF
   INTERFACE
      SUBROUTINE DPOTRF(Uplo,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPOTRF
   END INTERFACE
END MODULE S_DPOTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPOTRF
   INTERFACE
      SUBROUTINE SPOTRF(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPOTRF
   END INTERFACE
END MODULE S_SPOTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPOTRF
   INTERFACE
      SUBROUTINE ZPOTRF(Uplo,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPOTRF
   END INTERFACE
END MODULE S_ZPOTRF
