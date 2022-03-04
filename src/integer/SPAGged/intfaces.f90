!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ICMAX1
   INTERFACE
      FUNCTION ICMAX1(N,Cx,Incx)
      IMPLICIT NONE
      INTEGER :: ICMAX1
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION ICMAX1
   END INTERFACE
END MODULE S_ICMAX1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_IEEECK
   INTERFACE
      FUNCTION IEEECK(Ispec,Zero,One)
      IMPLICIT NONE
      INTEGER :: IEEECK
      INTEGER , INTENT(IN) :: Ispec
      REAL , INTENT(IN) :: Zero
      REAL , INTENT(IN) :: One
      END FUNCTION IEEECK
   END INTERFACE
END MODULE S_IEEECK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILACLC
   INTERFACE
      FUNCTION ILACLC(M,N,A,Lda)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER :: ILACLC
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END FUNCTION ILACLC
   END INTERFACE
END MODULE S_ILACLC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILACLR
   INTERFACE
      FUNCTION ILACLR(M,N,A,Lda)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER :: ILACLR
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END FUNCTION ILACLR
   END INTERFACE
END MODULE S_ILACLR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILADIAG
   INTERFACE
      FUNCTION ILADIAG(Diag)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  BLAS_NON_UNIT_DIAG = 131 ,               &
     &                         BLAS_UNIT_DIAG = 132
      INTEGER :: ILADIAG
      CHARACTER :: Diag
      END FUNCTION ILADIAG
   END INTERFACE
END MODULE S_ILADIAG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILADLC
   INTERFACE
      FUNCTION ILADLC(M,N,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER :: ILADLC
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END FUNCTION ILADLC
   END INTERFACE
END MODULE S_ILADLC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILADLR
   INTERFACE
      FUNCTION ILADLR(M,N,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER :: ILADLR
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END FUNCTION ILADLR
   END INTERFACE
END MODULE S_ILADLR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILAENV2STAGE
   INTERFACE
      FUNCTION ILAENV2STAGE(Ispec,Name,Opts,N1,N2,N3,N4)
      IMPLICIT NONE
      INTEGER :: ILAENV2STAGE
      INTEGER , INTENT(IN) :: Ispec
      CHARACTER(*) :: Name
      CHARACTER(*) :: Opts
      INTEGER :: N1
      INTEGER :: N2
      INTEGER :: N3
      INTEGER :: N4
      END FUNCTION ILAENV2STAGE
   END INTERFACE
END MODULE S_ILAENV2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILAENV
   INTERFACE
      FUNCTION ILAENV(Ispec,Name,Opts,N1,N2,N3,N4)
      IMPLICIT NONE
      INTEGER :: ILAENV
      INTEGER :: Ispec
      CHARACTER(*) :: Name
      CHARACTER(*) :: Opts
      INTEGER :: N1
      INTEGER :: N2
      INTEGER :: N3
      INTEGER :: N4
      END FUNCTION ILAENV
   END INTERFACE
END MODULE S_ILAENV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILAPREC
   INTERFACE
      FUNCTION ILAPREC(Prec)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  BLAS_PREC_SINGLE = 211 ,                 &
     &                         BLAS_PREC_DOUBLE = 212 ,                 &
     &                         BLAS_PREC_INDIGENOUS = 213 ,             &
     &                         BLAS_PREC_EXTRA = 214
      INTEGER :: ILAPREC
      CHARACTER :: Prec
      END FUNCTION ILAPREC
   END INTERFACE
END MODULE S_ILAPREC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILASLC
   INTERFACE
      FUNCTION ILASLC(M,N,A,Lda)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER :: ILASLC
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END FUNCTION ILASLC
   END INTERFACE
END MODULE S_ILASLC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILASLR
   INTERFACE
      FUNCTION ILASLR(M,N,A,Lda)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER :: ILASLR
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END FUNCTION ILASLR
   END INTERFACE
END MODULE S_ILASLR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILATRANS
   INTERFACE
      FUNCTION ILATRANS(Trans)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  BLAS_NO_TRANS = 111 , BLAS_TRANS = 112 , &
     &                         BLAS_CONJ_TRANS = 113
      INTEGER :: ILATRANS
      CHARACTER :: Trans
      END FUNCTION ILATRANS
   END INTERFACE
END MODULE S_ILATRANS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILAUPLO
   INTERFACE
      FUNCTION ILAUPLO(Uplo)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  BLAS_UPPER = 121 , BLAS_LOWER = 122
      INTEGER :: ILAUPLO
      CHARACTER :: Uplo
      END FUNCTION ILAUPLO
   END INTERFACE
END MODULE S_ILAUPLO
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILAZLC
   INTERFACE
      FUNCTION ILAZLC(M,N,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER :: ILAZLC
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END FUNCTION ILAZLC
   END INTERFACE
END MODULE S_ILAZLC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ILAZLR
   INTERFACE
      FUNCTION ILAZLR(M,N,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER :: ILAZLR
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END FUNCTION ILAZLR
   END INTERFACE
END MODULE S_ILAZLR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_IPARMQ
   INTERFACE
      FUNCTION IPARMQ(Ispec,Name,Opts,N,Ilo,Ihi,Lwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  INMIN = 12 , INWIN = 13 , INIBL = 14 ,   &
     &                         ISHFTS = 15 , IACC22 = 16 , NMIN = 75 ,  &
     &                         K22MIN = 14 , KACMIN = 14 , NIBBLE = 14 ,&
     &                         KNWSWP = 500
      REAL , PARAMETER  ::  TWO = 2.0
      INTEGER :: IPARMQ
      INTEGER , INTENT(IN) :: Ispec
      CHARACTER(*) , INTENT(IN) :: Name
      CHARACTER(*) :: Opts
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      INTEGER :: Lwork
      END FUNCTION IPARMQ
   END INTERFACE
END MODULE S_IPARMQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_IZMAX1
   INTERFACE
      FUNCTION IZMAX1(N,Zx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: IZMAX1
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION IZMAX1
   END INTERFACE
END MODULE S_IZMAX1
