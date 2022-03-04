!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEQRF
   INTERFACE
      SUBROUTINE CGEQRF(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEQRF
   END INTERFACE
END MODULE S_CGEQRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEQRF
   INTERFACE
      SUBROUTINE DGEQRF(M,N,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEQRF
   END INTERFACE
END MODULE S_DGEQRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SCEIL
   INTERFACE
      FUNCTION SCEIL(A)
      IMPLICIT NONE
      REAL :: SCEIL
      REAL , INTENT(IN) :: A
      END FUNCTION SCEIL
   END INTERFACE
END MODULE S_SCEIL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEQRF
   INTERFACE
      SUBROUTINE SGEQRF(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEQRF
   END INTERFACE
END MODULE S_SGEQRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEQRF
   INTERFACE
      SUBROUTINE ZGEQRF(M,N,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEQRF
   END INTERFACE
END MODULE S_ZGEQRF
