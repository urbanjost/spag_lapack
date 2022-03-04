!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEGS
   INTERFACE
      SUBROUTINE CGEGS(Jobvsl,Jobvsr,N,A,Lda,B,Ldb,Alpha,Beta,Vsl,Ldvsl,&
     &                 Vsr,Ldvsr,Work,Lwork,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Alpha
      COMPLEX , DIMENSION(*) :: Beta
      COMPLEX , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      COMPLEX , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEGS
   END INTERFACE
END MODULE S_CGEGS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEGV
   INTERFACE
      SUBROUTINE CGEGV(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alpha,Beta,Vl,Ldvl,Vr, &
     &                 Ldvr,Work,Lwork,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Alpha
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEGV
   END INTERFACE
END MODULE S_CGEGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGELSX
   INTERFACE
      SUBROUTINE CGELSX(M,N,Nrhs,A,Lda,B,Ldb,Jpvt,Rcond,Rank,Work,Rwork,&
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  IMAX = 1 , IMIN = 2
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , DONE = ZERO ,&
     &                      NTDONE = ONE
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , DIMENSION(*) :: Jpvt
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGELSX
   END INTERFACE
END MODULE S_CGELSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEQPF
   INTERFACE
      SUBROUTINE CGEQPF(M,N,A,Lda,Jpvt,Tau,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEQPF
   END INTERFACE
END MODULE S_CGEQPF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGSVD
   INTERFACE
      SUBROUTINE CGGSVD(Jobu,Jobv,Jobq,M,N,P,K,L,A,Lda,B,Ldb,Alpha,Beta,&
     &                  U,Ldu,V,Ldv,Q,Ldq,Work,Rwork,Iwork,Info)
      IMPLICIT NONE
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      CHARACTER :: Jobq
      INTEGER :: M
      INTEGER :: N
      INTEGER :: P
      INTEGER :: K
      INTEGER :: L
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Alpha
      REAL , DIMENSION(*) :: Beta
      COMPLEX , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGSVD
   END INTERFACE
END MODULE S_CGGSVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGSVP
   INTERFACE
      SUBROUTINE CGGSVP(Jobu,Jobv,Jobq,M,P,N,A,Lda,B,Ldb,Tola,Tolb,K,L, &
     &                  U,Ldu,V,Ldv,Q,Ldq,Iwork,Rwork,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      CHARACTER :: Jobq
      INTEGER :: M
      INTEGER :: P
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(IN) :: Tola
      REAL , INTENT(IN) :: Tolb
      INTEGER , INTENT(INOUT) :: K
      INTEGER , INTENT(INOUT) :: L
      COMPLEX , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , DIMENSION(*) :: Iwork
      REAL , DIMENSION(*) :: Rwork
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGSVP
   END INTERFACE
END MODULE S_CGGSVP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAHRD
   INTERFACE
      SUBROUTINE CLAHRD(N,K,Nb,A,Lda,Tau,T,Ldt,Y,Ldy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Nb) :: Tau
      COMPLEX , DIMENSION(Ldt,Nb) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(Ldy,Nb) :: Y
      INTEGER :: Ldy
      END SUBROUTINE CLAHRD
   END INTERFACE
END MODULE S_CLAHRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLATZM
   INTERFACE
      SUBROUTINE CLATZM(Side,M,N,V,Incv,Tau,C1,C2,Ldc,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: V
      INTEGER :: Incv
      COMPLEX :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C1
      COMPLEX , DIMENSION(Ldc,*) :: C2
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      END SUBROUTINE CLATZM
   END INTERFACE
END MODULE S_CLATZM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTZRQF
   INTERFACE
      SUBROUTINE CTZRQF(M,N,A,Lda,Tau,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTZRQF
   END INTERFACE
END MODULE S_CTZRQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEGS
   INTERFACE
      SUBROUTINE DGEGS(Jobvsl,Jobvsr,N,A,Lda,B,Ldb,Alphar,Alphai,Beta,  &
     &                 Vsl,Ldvsl,Vsr,Ldvsr,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , DIMENSION(*) :: Beta
      REAL(R8KIND) , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      REAL(R8KIND) , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEGS
   END INTERFACE
END MODULE S_DGEGS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEGV
   INTERFACE
      SUBROUTINE DGEGV(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alphar,Alphai,Beta,Vl, &
     &                 Ldvl,Vr,Ldvr,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEGV
   END INTERFACE
END MODULE S_DGEGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGELSX
   INTERFACE
      SUBROUTINE DGELSX(M,N,Nrhs,A,Lda,B,Ldb,Jpvt,Rcond,Rank,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  IMAX = 1 , IMIN = 2
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              DONE = ZERO , NTDONE = ONE
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , DIMENSION(*) :: Jpvt
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGELSX
   END INTERFACE
END MODULE S_DGELSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEQPF
   INTERFACE
      SUBROUTINE DGEQPF(M,N,A,Lda,Jpvt,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEQPF
   END INTERFACE
END MODULE S_DGEQPF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGSVD
   INTERFACE
      SUBROUTINE DGGSVD(Jobu,Jobv,Jobq,M,N,P,K,L,A,Lda,B,Ldb,Alpha,Beta,&
     &                  U,Ldu,V,Ldv,Q,Ldq,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      CHARACTER :: Jobq
      INTEGER :: M
      INTEGER :: N
      INTEGER :: P
      INTEGER :: K
      INTEGER :: L
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Alpha
      REAL(R8KIND) , DIMENSION(*) :: Beta
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGSVD
   END INTERFACE
END MODULE S_DGGSVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGSVP
   INTERFACE
      SUBROUTINE DGGSVP(Jobu,Jobv,Jobq,M,P,N,A,Lda,B,Ldb,Tola,Tolb,K,L, &
     &                  U,Ldu,V,Ldv,Q,Ldq,Iwork,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      CHARACTER :: Jobq
      INTEGER :: M
      INTEGER :: P
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Tola
      REAL(R8KIND) , INTENT(IN) :: Tolb
      INTEGER , INTENT(INOUT) :: K
      INTEGER , INTENT(INOUT) :: L
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , DIMENSION(*) :: Iwork
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGSVP
   END INTERFACE
END MODULE S_DGGSVP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAHRD
   INTERFACE
      SUBROUTINE DLAHRD(N,K,Nb,A,Lda,Tau,T,Ldt,Y,Ldy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Nb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Nb) :: Tau
      REAL(R8KIND) , DIMENSION(Ldt,Nb) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Ldy,Nb) :: Y
      INTEGER :: Ldy
      END SUBROUTINE DLAHRD
   END INTERFACE
END MODULE S_DLAHRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLATZM
   INTERFACE
      SUBROUTINE DLATZM(Side,M,N,V,Incv,Tau,C1,C2,Ldc,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: V
      INTEGER :: Incv
      REAL(R8KIND) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C1
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C2
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      END SUBROUTINE DLATZM
   END INTERFACE
END MODULE S_DLATZM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTZRQF
   INTERFACE
      SUBROUTINE DTZRQF(M,N,A,Lda,Tau,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTZRQF
   END INTERFACE
END MODULE S_DTZRQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEGS
   INTERFACE
      SUBROUTINE SGEGS(Jobvsl,Jobvsr,N,A,Lda,B,Ldb,Alphar,Alphai,Beta,  &
     &                 Vsl,Ldvsl,Vsr,Ldvsr,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Alphar
      REAL , DIMENSION(*) :: Alphai
      REAL , DIMENSION(*) :: Beta
      REAL , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      REAL , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEGS
   END INTERFACE
END MODULE S_SGEGS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEGV
   INTERFACE
      SUBROUTINE SGEGV(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alphar,Alphai,Beta,Vl, &
     &                 Ldvl,Vr,Ldvr,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alphar
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEGV
   END INTERFACE
END MODULE S_SGEGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGELSX
   INTERFACE
      SUBROUTINE SGELSX(M,N,Nrhs,A,Lda,B,Ldb,Jpvt,Rcond,Rank,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  IMAX = 1 , IMIN = 2
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , DONE = ZERO ,  &
     &                      NTDONE = ONE
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , DIMENSION(*) :: Jpvt
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGELSX
   END INTERFACE
END MODULE S_SGELSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEQPF
   INTERFACE
      SUBROUTINE SGEQPF(M,N,A,Lda,Jpvt,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      REAL , DIMENSION(*) :: Tau
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEQPF
   END INTERFACE
END MODULE S_SGEQPF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGSVD
   INTERFACE
      SUBROUTINE SGGSVD(Jobu,Jobv,Jobq,M,N,P,K,L,A,Lda,B,Ldb,Alpha,Beta,&
     &                  U,Ldu,V,Ldv,Q,Ldq,Work,Iwork,Info)
      IMPLICIT NONE
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      CHARACTER :: Jobq
      INTEGER :: M
      INTEGER :: N
      INTEGER :: P
      INTEGER :: K
      INTEGER :: L
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Alpha
      REAL , DIMENSION(*) :: Beta
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGSVD
   END INTERFACE
END MODULE S_SGGSVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGSVP
   INTERFACE
      SUBROUTINE SGGSVP(Jobu,Jobv,Jobq,M,P,N,A,Lda,B,Ldb,Tola,Tolb,K,L, &
     &                  U,Ldu,V,Ldv,Q,Ldq,Iwork,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      CHARACTER :: Jobq
      INTEGER :: M
      INTEGER :: P
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(IN) :: Tola
      REAL , INTENT(IN) :: Tolb
      INTEGER , INTENT(INOUT) :: K
      INTEGER , INTENT(INOUT) :: L
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , DIMENSION(*) :: Iwork
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGSVP
   END INTERFACE
END MODULE S_SGGSVP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAHRD
   INTERFACE
      SUBROUTINE SLAHRD(N,K,Nb,A,Lda,Tau,T,Ldt,Y,Ldy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Nb
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Nb) :: Tau
      REAL , DIMENSION(Ldt,Nb) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Ldy,Nb) :: Y
      INTEGER :: Ldy
      END SUBROUTINE SLAHRD
   END INTERFACE
END MODULE S_SLAHRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLATZM
   INTERFACE
      SUBROUTINE SLATZM(Side,M,N,V,Incv,Tau,C1,C2,Ldc,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(*) :: V
      INTEGER :: Incv
      REAL :: Tau
      REAL , DIMENSION(Ldc,*) :: C1
      REAL , DIMENSION(Ldc,*) :: C2
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      END SUBROUTINE SLATZM
   END INTERFACE
END MODULE S_SLATZM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STZRQF
   INTERFACE
      SUBROUTINE STZRQF(M,N,A,Lda,Tau,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STZRQF
   END INTERFACE
END MODULE S_STZRQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEGS
   INTERFACE
      SUBROUTINE ZGEGS(Jobvsl,Jobvsr,N,A,Lda,B,Ldb,Alpha,Beta,Vsl,Ldvsl,&
     &                 Vsr,Ldvsr,Work,Lwork,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Alpha
      COMPLEX(CX16KIND) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      COMPLEX(CX16KIND) , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEGS
   END INTERFACE
END MODULE S_ZGEGS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEGV
   INTERFACE
      SUBROUTINE ZGEGV(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alpha,Beta,Vl,Ldvl,Vr, &
     &                 Ldvr,Work,Lwork,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Alpha
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEGV
   END INTERFACE
END MODULE S_ZGEGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGELSX
   INTERFACE
      SUBROUTINE ZGELSX(M,N,Nrhs,A,Lda,B,Ldb,Jpvt,Rcond,Rank,Work,Rwork,&
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  IMAX = 1 , IMIN = 2
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              DONE = ZERO , NTDONE = ONE
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , DIMENSION(*) :: Jpvt
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGELSX
   END INTERFACE
END MODULE S_ZGELSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEQPF
   INTERFACE
      SUBROUTINE ZGEQPF(M,N,A,Lda,Jpvt,Tau,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEQPF
   END INTERFACE
END MODULE S_ZGEQPF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGSVD
   INTERFACE
      SUBROUTINE ZGGSVD(Jobu,Jobv,Jobq,M,N,P,K,L,A,Lda,B,Ldb,Alpha,Beta,&
     &                  U,Ldu,V,Ldv,Q,Ldq,Work,Rwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      CHARACTER :: Jobq
      INTEGER :: M
      INTEGER :: N
      INTEGER :: P
      INTEGER :: K
      INTEGER :: L
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Alpha
      REAL(R8KIND) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGSVD
   END INTERFACE
END MODULE S_ZGGSVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGSVP
   INTERFACE
      SUBROUTINE ZGGSVP(Jobu,Jobv,Jobq,M,P,N,A,Lda,B,Ldb,Tola,Tolb,K,L, &
     &                  U,Ldu,V,Ldv,Q,Ldq,Iwork,Rwork,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      CHARACTER :: Jobq
      INTEGER :: M
      INTEGER :: P
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Tola
      REAL(R8KIND) , INTENT(IN) :: Tolb
      INTEGER , INTENT(INOUT) :: K
      INTEGER , INTENT(INOUT) :: L
      COMPLEX(CX16KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , DIMENSION(*) :: Iwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGSVP
   END INTERFACE
END MODULE S_ZGGSVP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAHRD
   INTERFACE
      SUBROUTINE ZLAHRD(N,K,Nb,A,Lda,Tau,T,Ldt,Y,Ldy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Nb) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldt,Nb) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Ldy,Nb) :: Y
      INTEGER :: Ldy
      END SUBROUTINE ZLAHRD
   END INTERFACE
END MODULE S_ZLAHRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLATZM
   INTERFACE
      SUBROUTINE ZLATZM(Side,M,N,V,Incv,Tau,C1,C2,Ldc,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: V
      INTEGER :: Incv
      COMPLEX(CX16KIND) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C1
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C2
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      END SUBROUTINE ZLATZM
   END INTERFACE
END MODULE S_ZLATZM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTZRQF
   INTERFACE
      SUBROUTINE ZTZRQF(M,N,A,Lda,Tau,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTZRQF
   END INTERFACE
END MODULE S_ZTZRQF
