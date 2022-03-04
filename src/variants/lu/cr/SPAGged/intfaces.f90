!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGETRF
   INTERFACE
      SUBROUTINE CGETRF(M,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGETRF
   END INTERFACE
END MODULE S_CGETRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGETRF
   INTERFACE
      SUBROUTINE DGETRF(M,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGETRF
   END INTERFACE
END MODULE S_DGETRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGETRF
   INTERFACE
      SUBROUTINE SGETRF(M,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGETRF
   END INTERFACE
END MODULE S_SGETRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGETRF
   INTERFACE
      SUBROUTINE ZGETRF(M,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGETRF
   END INTERFACE
END MODULE S_ZGETRF
