!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CAXPY
   INTERFACE
      SUBROUTINE CAXPY(N,Ca,Cx,Incx,Cy,Incy)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX :: Ca
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE CAXPY
   END INTERFACE
END MODULE S_CAXPY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CCOPY
   INTERFACE
      SUBROUTINE CCOPY(N,Cx,Incx,Cy,Incy)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE CCOPY
   END INTERFACE
END MODULE S_CCOPY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CDOTC
   INTERFACE
      FUNCTION CDOTC(N,Cx,Incx,Cy,Incy)
      IMPLICIT NONE
      COMPLEX :: CDOTC
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
      END FUNCTION CDOTC
   END INTERFACE
END MODULE S_CDOTC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CDOTU
   INTERFACE
      FUNCTION CDOTU(N,Cx,Incx,Cy,Incy)
      IMPLICIT NONE
      COMPLEX :: CDOTU
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
      END FUNCTION CDOTU
   END INTERFACE
END MODULE S_CDOTU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBMV
   INTERFACE
      SUBROUTINE CGBMV(Trans,M,N,Kl,Ku,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE CGBMV
   END INTERFACE
END MODULE S_CGBMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEMM
   INTERFACE
      SUBROUTINE CGEMM(Transa,Transb,M,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Transa
      CHARACTER :: Transb
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE CGEMM
   END INTERFACE
END MODULE S_CGEMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEMV
   INTERFACE
      SUBROUTINE CGEMV(Trans,M,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE CGEMV
   END INTERFACE
END MODULE S_CGEMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGERC
   INTERFACE
      SUBROUTINE CGERC(M,N,Alpha,X,Incx,Y,Incy,A,Lda)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE CGERC
   END INTERFACE
END MODULE S_CGERC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGERU
   INTERFACE
      SUBROUTINE CGERU(M,N,Alpha,X,Incx,Y,Incy,A,Lda)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE CGERU
   END INTERFACE
END MODULE S_CGERU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHBMV
   INTERFACE
      SUBROUTINE CHBMV(Uplo,N,K,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE CHBMV
   END INTERFACE
END MODULE S_CHBMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEMM
   INTERFACE
      SUBROUTINE CHEMM(Side,Uplo,M,N,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE CHEMM
   END INTERFACE
END MODULE S_CHEMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEMV
   INTERFACE
      SUBROUTINE CHEMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE CHEMV
   END INTERFACE
END MODULE S_CHEMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHER2
   INTERFACE
      SUBROUTINE CHER2(Uplo,N,Alpha,X,Incx,Y,Incy,A,Lda)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE CHER2
   END INTERFACE
END MODULE S_CHER2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHER2K
   INTERFACE
      SUBROUTINE CHER2K(Uplo,Trans,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE CHER2K
   END INTERFACE
END MODULE S_CHER2K
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHER
   INTERFACE
      SUBROUTINE CHER(Uplo,N,Alpha,X,Incx,A,Lda)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE CHER
   END INTERFACE
END MODULE S_CHER
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHERK
   INTERFACE
      SUBROUTINE CHERK(Uplo,Trans,N,K,Alpha,A,Lda,Beta,C,Ldc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE CHERK
   END INTERFACE
END MODULE S_CHERK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPMV
   INTERFACE
      SUBROUTINE CHPMV(Uplo,N,Alpha,Ap,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE CHPMV
   END INTERFACE
END MODULE S_CHPMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPR2
   INTERFACE
      SUBROUTINE CHPR2(Uplo,N,Alpha,X,Incx,Y,Incy,Ap)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      END SUBROUTINE CHPR2
   END INTERFACE
END MODULE S_CHPR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPR
   INTERFACE
      SUBROUTINE CHPR(Uplo,N,Alpha,X,Incx,Ap)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      END SUBROUTINE CHPR
   END INTERFACE
END MODULE S_CHPR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CROTG
   INTERFACE
      SUBROUTINE CROTG(Ca,Cb,C,S)
      IMPLICIT NONE
      COMPLEX , INTENT(INOUT) :: Ca
      COMPLEX , INTENT(IN) :: Cb
      REAL , INTENT(OUT) :: C
      COMPLEX , INTENT(OUT) :: S
      END SUBROUTINE CROTG
   END INTERFACE
END MODULE S_CROTG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSCAL
   INTERFACE
      SUBROUTINE CSCAL(N,Ca,Cx,Incx)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Ca
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE CSCAL
   END INTERFACE
END MODULE S_CSCAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSROT
   INTERFACE
      SUBROUTINE CSROT(N,Cx,Incx,Cy,Incy,C,S)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: S
      END SUBROUTINE CSROT
   END INTERFACE
END MODULE S_CSROT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSSCAL
   INTERFACE
      SUBROUTINE CSSCAL(N,Sa,Cx,Incx)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Sa
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE CSSCAL
   END INTERFACE
END MODULE S_CSSCAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSWAP
   INTERFACE
      SUBROUTINE CSWAP(N,Cx,Incx,Cy,Incy)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE CSWAP
   END INTERFACE
END MODULE S_CSWAP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYMM
   INTERFACE
      SUBROUTINE CSYMM(Side,Uplo,M,N,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE CSYMM
   END INTERFACE
END MODULE S_CSYMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYR2K
   INTERFACE
      SUBROUTINE CSYR2K(Uplo,Trans,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE CSYR2K
   END INTERFACE
END MODULE S_CSYR2K
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYRK
   INTERFACE
      SUBROUTINE CSYRK(Uplo,Trans,N,K,Alpha,A,Lda,Beta,C,Ldc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE CSYRK
   END INTERFACE
END MODULE S_CSYRK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTBMV
   INTERFACE
      SUBROUTINE CTBMV(Uplo,Trans,Diag,N,K,A,Lda,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE CTBMV
   END INTERFACE
END MODULE S_CTBMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTBSV
   INTERFACE
      SUBROUTINE CTBSV(Uplo,Trans,Diag,N,K,A,Lda,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE CTBSV
   END INTERFACE
END MODULE S_CTBSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPMV
   INTERFACE
      SUBROUTINE CTPMV(Uplo,Trans,Diag,N,Ap,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE CTPMV
   END INTERFACE
END MODULE S_CTPMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPSV
   INTERFACE
      SUBROUTINE CTPSV(Uplo,Trans,Diag,N,Ap,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE CTPSV
   END INTERFACE
END MODULE S_CTPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRMM
   INTERFACE
      SUBROUTINE CTRMM(Side,Uplo,Transa,Diag,M,N,Alpha,A,Lda,B,Ldb)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Transa
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE CTRMM
   END INTERFACE
END MODULE S_CTRMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRMV
   INTERFACE
      SUBROUTINE CTRMV(Uplo,Trans,Diag,N,A,Lda,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE CTRMV
   END INTERFACE
END MODULE S_CTRMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRSM
   INTERFACE
      SUBROUTINE CTRSM(Side,Uplo,Transa,Diag,M,N,Alpha,A,Lda,B,Ldb)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Transa
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE CTRSM
   END INTERFACE
END MODULE S_CTRSM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRSV
   INTERFACE
      SUBROUTINE CTRSV(Uplo,Trans,Diag,N,A,Lda,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE CTRSV
   END INTERFACE
END MODULE S_CTRSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DASUM
   INTERFACE
      FUNCTION DASUM(N,Dx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DASUM
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION DASUM
   END INTERFACE
END MODULE S_DASUM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DAXPY
   INTERFACE
      SUBROUTINE DAXPY(N,Da,Dx,Incx,Dy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Da
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dy
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE DAXPY
   END INTERFACE
END MODULE S_DAXPY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DCABS1
   INTERFACE
      FUNCTION DCABS1(Z)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DCABS1
      COMPLEX(CX16KIND) , INTENT(IN) :: Z
      END FUNCTION DCABS1
   END INTERFACE
END MODULE S_DCABS1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DCOPY
   INTERFACE
      SUBROUTINE DCOPY(N,Dx,Incx,Dy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Dy
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE DCOPY
   END INTERFACE
END MODULE S_DCOPY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DDOT
   INTERFACE
      FUNCTION DDOT(N,Dx,Incx,Dy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DDOT
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Dy
      INTEGER , INTENT(IN) :: Incy
      END FUNCTION DDOT
   END INTERFACE
END MODULE S_DDOT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBMV
   INTERFACE
      SUBROUTINE DGBMV(Trans,M,N,Kl,Ku,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE DGBMV
   END INTERFACE
END MODULE S_DGBMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEMM
   INTERFACE
      SUBROUTINE DGEMM(Transa,Transb,M,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Transa
      CHARACTER :: Transb
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE DGEMM
   END INTERFACE
END MODULE S_DGEMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEMV
   INTERFACE
      SUBROUTINE DGEMV(Trans,M,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE DGEMV
   END INTERFACE
END MODULE S_DGEMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGER
   INTERFACE
      SUBROUTINE DGER(M,N,Alpha,X,Incx,Y,Incy,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE DGER
   END INTERFACE
END MODULE S_DGER
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DNRM2
   INTERFACE
      FUNCTION DNRM2(N,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DNRM2
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION DNRM2
   END INTERFACE
END MODULE S_DNRM2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DROT
   INTERFACE
      SUBROUTINE DROT(N,Dx,Incx,Dy,Incy,C,S)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dy
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(IN) :: C
      REAL(R8KIND) , INTENT(IN) :: S
      END SUBROUTINE DROT
   END INTERFACE
END MODULE S_DROT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DROTG
   INTERFACE
      SUBROUTINE DROTG(Da,Db,C,S)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) , INTENT(INOUT) :: Da
      REAL(R8KIND) , INTENT(INOUT) :: Db
      REAL(R8KIND) , INTENT(INOUT) :: C
      REAL(R8KIND) , INTENT(INOUT) :: S
      END SUBROUTINE DROTG
   END INTERFACE
END MODULE S_DROTG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DROTM
   INTERFACE
      SUBROUTINE DROTM(N,Dx,Incx,Dy,Incy,Dparam)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dy
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(IN) , DIMENSION(5) :: Dparam
      END SUBROUTINE DROTM
   END INTERFACE
END MODULE S_DROTM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DROTMG
   INTERFACE
      SUBROUTINE DROTMG(Dd1,Dd2,Dx1,Dy1,Dparam)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) , INTENT(INOUT) :: Dd1
      REAL(R8KIND) , INTENT(INOUT) :: Dd2
      REAL(R8KIND) , INTENT(INOUT) :: Dx1
      REAL(R8KIND) , INTENT(IN) :: Dy1
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(5) :: Dparam
      END SUBROUTINE DROTMG
   END INTERFACE
END MODULE S_DROTMG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSBMV
   INTERFACE
      SUBROUTINE DSBMV(Uplo,N,K,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE DSBMV
   END INTERFACE
END MODULE S_DSBMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSCAL
   INTERFACE
      SUBROUTINE DSCAL(N,Da,Dx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Da
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE DSCAL
   END INTERFACE
END MODULE S_DSCAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSDOT
   INTERFACE
      FUNCTION DSDOT(N,Sx,Incx,Sy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DSDOT
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
      END FUNCTION DSDOT
   END INTERFACE
END MODULE S_DSDOT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPMV
   INTERFACE
      SUBROUTINE DSPMV(Uplo,N,Alpha,Ap,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE DSPMV
   END INTERFACE
END MODULE S_DSPMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPR2
   INTERFACE
      SUBROUTINE DSPR2(Uplo,N,Alpha,X,Incx,Y,Incy,Ap)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      END SUBROUTINE DSPR2
   END INTERFACE
END MODULE S_DSPR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPR
   INTERFACE
      SUBROUTINE DSPR(Uplo,N,Alpha,X,Incx,Ap)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      END SUBROUTINE DSPR
   END INTERFACE
END MODULE S_DSPR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSWAP
   INTERFACE
      SUBROUTINE DSWAP(N,Dx,Incx,Dy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dy
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE DSWAP
   END INTERFACE
END MODULE S_DSWAP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYMM
   INTERFACE
      SUBROUTINE DSYMM(Side,Uplo,M,N,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Side
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE DSYMM
   END INTERFACE
END MODULE S_DSYMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYMV
   INTERFACE
      SUBROUTINE DSYMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE DSYMV
   END INTERFACE
END MODULE S_DSYMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYR2
   INTERFACE
      SUBROUTINE DSYR2(Uplo,N,Alpha,X,Incx,Y,Incy,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE DSYR2
   END INTERFACE
END MODULE S_DSYR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYR2K
   INTERFACE
      SUBROUTINE DSYR2K(Uplo,Trans,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE DSYR2K
   END INTERFACE
END MODULE S_DSYR2K
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYR
   INTERFACE
      SUBROUTINE DSYR(Uplo,N,Alpha,X,Incx,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE DSYR
   END INTERFACE
END MODULE S_DSYR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYRK
   INTERFACE
      SUBROUTINE DSYRK(Uplo,Trans,N,K,Alpha,A,Lda,Beta,C,Ldc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE DSYRK
   END INTERFACE
END MODULE S_DSYRK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTBMV
   INTERFACE
      SUBROUTINE DTBMV(Uplo,Trans,Diag,N,K,A,Lda,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE DTBMV
   END INTERFACE
END MODULE S_DTBMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTBSV
   INTERFACE
      SUBROUTINE DTBSV(Uplo,Trans,Diag,N,K,A,Lda,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE DTBSV
   END INTERFACE
END MODULE S_DTBSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPMV
   INTERFACE
      SUBROUTINE DTPMV(Uplo,Trans,Diag,N,Ap,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE DTPMV
   END INTERFACE
END MODULE S_DTPMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPSV
   INTERFACE
      SUBROUTINE DTPSV(Uplo,Trans,Diag,N,Ap,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE DTPSV
   END INTERFACE
END MODULE S_DTPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRMM
   INTERFACE
      SUBROUTINE DTRMM(Side,Uplo,Transa,Diag,M,N,Alpha,A,Lda,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Transa
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE DTRMM
   END INTERFACE
END MODULE S_DTRMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRMV
   INTERFACE
      SUBROUTINE DTRMV(Uplo,Trans,Diag,N,A,Lda,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE DTRMV
   END INTERFACE
END MODULE S_DTRMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRSM
   INTERFACE
      SUBROUTINE DTRSM(Side,Uplo,Transa,Diag,M,N,Alpha,A,Lda,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Transa
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE DTRSM
   END INTERFACE
END MODULE S_DTRSM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRSV
   INTERFACE
      SUBROUTINE DTRSV(Uplo,Trans,Diag,N,A,Lda,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE DTRSV
   END INTERFACE
END MODULE S_DTRSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DZASUM
   INTERFACE
      FUNCTION DZASUM(N,Zx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DZASUM
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION DZASUM
   END INTERFACE
END MODULE S_DZASUM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DZNRM2
   INTERFACE
      FUNCTION DZNRM2(N,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DZNRM2
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION DZNRM2
   END INTERFACE
END MODULE S_DZNRM2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ICAMAX
   INTERFACE
      FUNCTION ICAMAX(N,Cx,Incx)
      IMPLICIT NONE
      INTEGER :: ICAMAX
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION ICAMAX
   END INTERFACE
END MODULE S_ICAMAX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_IDAMAX
   INTERFACE
      FUNCTION IDAMAX(N,Dx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: IDAMAX
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION IDAMAX
   END INTERFACE
END MODULE S_IDAMAX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ISAMAX
   INTERFACE
      FUNCTION ISAMAX(N,Sx,Incx)
      IMPLICIT NONE
      INTEGER :: ISAMAX
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION ISAMAX
   END INTERFACE
END MODULE S_ISAMAX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_IZAMAX
   INTERFACE
      FUNCTION IZAMAX(N,Zx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: IZAMAX
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION IZAMAX
   END INTERFACE
END MODULE S_IZAMAX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_LSAME
   INTERFACE
      FUNCTION LSAME(Ca,Cb)
      IMPLICIT NONE
      LOGICAL :: LSAME
      CHARACTER , INTENT(IN) :: Ca
      CHARACTER , INTENT(IN) :: Cb
      END FUNCTION LSAME
   END INTERFACE
END MODULE S_LSAME
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SASUM
   INTERFACE
      FUNCTION SASUM(N,Sx,Incx)
      IMPLICIT NONE
      REAL :: SASUM
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION SASUM
   END INTERFACE
END MODULE S_SASUM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SAXPY
   INTERFACE
      SUBROUTINE SAXPY(N,Sa,Sx,Incx,Sy,Incy)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Sa
      REAL , INTENT(IN) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE SAXPY
   END INTERFACE
END MODULE S_SAXPY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SCABS1
   INTERFACE
      FUNCTION SCABS1(Z)
      IMPLICIT NONE
      REAL :: SCABS1
      COMPLEX , INTENT(IN) :: Z
      END FUNCTION SCABS1
   END INTERFACE
END MODULE S_SCABS1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SCASUM
   INTERFACE
      FUNCTION SCASUM(N,Cx,Incx)
      IMPLICIT NONE
      REAL :: SCASUM
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION SCASUM
   END INTERFACE
END MODULE S_SCASUM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SCNRM2
   INTERFACE
      FUNCTION SCNRM2(N,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SCNRM2
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION SCNRM2
   END INTERFACE
END MODULE S_SCNRM2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SCOPY
   INTERFACE
      SUBROUTINE SCOPY(N,Sx,Incx,Sy,Incy)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(OUT) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE SCOPY
   END INTERFACE
END MODULE S_SCOPY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SDOT
   INTERFACE
      FUNCTION SDOT(N,Sx,Incx,Sy,Incy)
      IMPLICIT NONE
      REAL :: SDOT
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
      END FUNCTION SDOT
   END INTERFACE
END MODULE S_SDOT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SDSDOT
   INTERFACE
      FUNCTION SDSDOT(N,Sb,Sx,Incx,Sy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL :: SDSDOT
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Sb
      REAL , INTENT(IN) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
      END FUNCTION SDSDOT
   END INTERFACE
END MODULE S_SDSDOT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBMV
   INTERFACE
      SUBROUTINE SGBMV(Trans,M,N,Kl,Ku,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE SGBMV
   END INTERFACE
END MODULE S_SGBMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEMM
   INTERFACE
      SUBROUTINE SGEMM(Transa,Transb,M,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Transa
      CHARACTER :: Transb
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE SGEMM
   END INTERFACE
END MODULE S_SGEMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEMV
   INTERFACE
      SUBROUTINE SGEMV(Trans,M,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE SGEMV
   END INTERFACE
END MODULE S_SGEMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGER
   INTERFACE
      SUBROUTINE SGER(M,N,Alpha,X,Incx,Y,Incy,A,Lda)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE SGER
   END INTERFACE
END MODULE S_SGER
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SNRM2
   INTERFACE
      FUNCTION SNRM2(N,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SNRM2
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION SNRM2
   END INTERFACE
END MODULE S_SNRM2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SROT
   INTERFACE
      SUBROUTINE SROT(N,Sx,Incx,Sy,Incy,C,S)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: S
      END SUBROUTINE SROT
   END INTERFACE
END MODULE S_SROT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SROTG
   INTERFACE
      SUBROUTINE SROTG(Sa,Sb,C,S)
      IMPLICIT NONE
      REAL , INTENT(INOUT) :: Sa
      REAL , INTENT(INOUT) :: Sb
      REAL , INTENT(INOUT) :: C
      REAL , INTENT(INOUT) :: S
      END SUBROUTINE SROTG
   END INTERFACE
END MODULE S_SROTG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SROTM
   INTERFACE
      SUBROUTINE SROTM(N,Sx,Incx,Sy,Incy,Sparam)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(IN) , DIMENSION(5) :: Sparam
      END SUBROUTINE SROTM
   END INTERFACE
END MODULE S_SROTM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SROTMG
   INTERFACE
      SUBROUTINE SROTMG(Sd1,Sd2,Sx1,Sy1,Sparam)
      IMPLICIT NONE
      REAL , INTENT(INOUT) :: Sd1
      REAL , INTENT(INOUT) :: Sd2
      REAL , INTENT(INOUT) :: Sx1
      REAL , INTENT(IN) :: Sy1
      REAL , INTENT(OUT) , DIMENSION(5) :: Sparam
      END SUBROUTINE SROTMG
   END INTERFACE
END MODULE S_SROTMG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSBMV
   INTERFACE
      SUBROUTINE SSBMV(Uplo,N,K,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE SSBMV
   END INTERFACE
END MODULE S_SSBMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSCAL
   INTERFACE
      SUBROUTINE SSCAL(N,Sa,Sx,Incx)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Sa
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE SSCAL
   END INTERFACE
END MODULE S_SSCAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPMV
   INTERFACE
      SUBROUTINE SSPMV(Uplo,N,Alpha,Ap,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(*) :: Ap
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE SSPMV
   END INTERFACE
END MODULE S_SSPMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPR2
   INTERFACE
      SUBROUTINE SSPR2(Uplo,N,Alpha,X,Incx,Y,Incy,Ap)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ap
      END SUBROUTINE SSPR2
   END INTERFACE
END MODULE S_SSPR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPR
   INTERFACE
      SUBROUTINE SSPR(Uplo,N,Alpha,X,Incx,Ap)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ap
      END SUBROUTINE SSPR
   END INTERFACE
END MODULE S_SSPR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSWAP
   INTERFACE
      SUBROUTINE SSWAP(N,Sx,Incx,Sy,Incy)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE SSWAP
   END INTERFACE
END MODULE S_SSWAP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYMM
   INTERFACE
      SUBROUTINE SSYMM(Side,Uplo,M,N,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Side
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE SSYMM
   END INTERFACE
END MODULE S_SSYMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYMV
   INTERFACE
      SUBROUTINE SSYMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE SSYMV
   END INTERFACE
END MODULE S_SSYMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYR2
   INTERFACE
      SUBROUTINE SSYR2(Uplo,N,Alpha,X,Incx,Y,Incy,A,Lda)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE SSYR2
   END INTERFACE
END MODULE S_SSYR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYR2K
   INTERFACE
      SUBROUTINE SSYR2K(Uplo,Trans,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE SSYR2K
   END INTERFACE
END MODULE S_SSYR2K
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYR
   INTERFACE
      SUBROUTINE SSYR(Uplo,N,Alpha,X,Incx,A,Lda)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE SSYR
   END INTERFACE
END MODULE S_SSYR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYRK
   INTERFACE
      SUBROUTINE SSYRK(Uplo,Trans,N,K,Alpha,A,Lda,Beta,C,Ldc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE SSYRK
   END INTERFACE
END MODULE S_SSYRK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STBMV
   INTERFACE
      SUBROUTINE STBMV(Uplo,Trans,Diag,N,K,A,Lda,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE STBMV
   END INTERFACE
END MODULE S_STBMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STBSV
   INTERFACE
      SUBROUTINE STBSV(Uplo,Trans,Diag,N,K,A,Lda,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE STBSV
   END INTERFACE
END MODULE S_STBSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPMV
   INTERFACE
      SUBROUTINE STPMV(Uplo,Trans,Diag,N,Ap,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: Ap
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE STPMV
   END INTERFACE
END MODULE S_STPMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPSV
   INTERFACE
      SUBROUTINE STPSV(Uplo,Trans,Diag,N,Ap,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: Ap
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE STPSV
   END INTERFACE
END MODULE S_STPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRMM
   INTERFACE
      SUBROUTINE STRMM(Side,Uplo,Transa,Diag,M,N,Alpha,A,Lda,B,Ldb)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Transa
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE STRMM
   END INTERFACE
END MODULE S_STRMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRMV
   INTERFACE
      SUBROUTINE STRMV(Uplo,Trans,Diag,N,A,Lda,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE STRMV
   END INTERFACE
END MODULE S_STRMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRSM
   INTERFACE
      SUBROUTINE STRSM(Side,Uplo,Transa,Diag,M,N,Alpha,A,Lda,B,Ldb)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Transa
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE STRSM
   END INTERFACE
END MODULE S_STRSM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRSV
   INTERFACE
      SUBROUTINE STRSV(Uplo,Trans,Diag,N,A,Lda,X,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE STRSV
   END INTERFACE
END MODULE S_STRSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_XERBLA_ARRAY
   INTERFACE
      SUBROUTINE XERBLA_ARRAY(Srname_array,Srname_len,Info)
      IMPLICIT NONE
      CHARACTER(1) , INTENT(IN) , DIMENSION(Srname_len) :: Srname_array
      INTEGER , INTENT(IN) :: Srname_len
      INTEGER :: Info
      END SUBROUTINE XERBLA_ARRAY
   END INTERFACE
END MODULE S_XERBLA_ARRAY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_XERBLA
   INTERFACE
      SUBROUTINE XERBLA(Srname,Info)
      IMPLICIT NONE
      CHARACTER(*) , INTENT(IN) :: Srname
      INTEGER , INTENT(IN) :: Info
      END SUBROUTINE XERBLA
   END INTERFACE
END MODULE S_XERBLA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZAXPY
   INTERFACE
      SUBROUTINE ZAXPY(N,Za,Zx,Incx,Zy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) :: Za
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Zy
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE ZAXPY
   END INTERFACE
END MODULE S_ZAXPY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZCOPY
   INTERFACE
      SUBROUTINE ZCOPY(N,Zx,Incx,Zy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Zy
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE ZCOPY
   END INTERFACE
END MODULE S_ZCOPY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZDOTC
   INTERFACE
      FUNCTION ZDOTC(N,Zx,Incx,Zy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
      COMPLEX(CX16KIND) :: ZDOTC
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Zy
      INTEGER , INTENT(IN) :: Incy
      END FUNCTION ZDOTC
   END INTERFACE
END MODULE S_ZDOTC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZDOTU
   INTERFACE
      FUNCTION ZDOTU(N,Zx,Incx,Zy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
      COMPLEX(CX16KIND) :: ZDOTU
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Zy
      INTEGER , INTENT(IN) :: Incy
      END FUNCTION ZDOTU
   END INTERFACE
END MODULE S_ZDOTU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZDROT
   INTERFACE
      SUBROUTINE ZDROT(N,Zx,Incx,Zy,Incy,C,S)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Zy
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(IN) :: C
      REAL(R8KIND) , INTENT(IN) :: S
      END SUBROUTINE ZDROT
   END INTERFACE
END MODULE S_ZDROT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZDSCAL
   INTERFACE
      SUBROUTINE ZDSCAL(N,Da,Zx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Da
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE ZDSCAL
   END INTERFACE
END MODULE S_ZDSCAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBMV
   INTERFACE
      SUBROUTINE ZGBMV(Trans,M,N,Kl,Ku,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE ZGBMV
   END INTERFACE
END MODULE S_ZGBMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEMM
   INTERFACE
      SUBROUTINE ZGEMM(Transa,Transb,M,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Transa
      CHARACTER :: Transb
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE ZGEMM
   END INTERFACE
END MODULE S_ZGEMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEMV
   INTERFACE
      SUBROUTINE ZGEMV(Trans,M,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE ZGEMV
   END INTERFACE
END MODULE S_ZGEMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGERC
   INTERFACE
      SUBROUTINE ZGERC(M,N,Alpha,X,Incx,Y,Incy,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE ZGERC
   END INTERFACE
END MODULE S_ZGERC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGERU
   INTERFACE
      SUBROUTINE ZGERU(M,N,Alpha,X,Incx,Y,Incy,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE ZGERU
   END INTERFACE
END MODULE S_ZGERU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHBMV
   INTERFACE
      SUBROUTINE ZHBMV(Uplo,N,K,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE ZHBMV
   END INTERFACE
END MODULE S_ZHBMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEMM
   INTERFACE
      SUBROUTINE ZHEMM(Side,Uplo,M,N,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE ZHEMM
   END INTERFACE
END MODULE S_ZHEMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEMV
   INTERFACE
      SUBROUTINE ZHEMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE ZHEMV
   END INTERFACE
END MODULE S_ZHEMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHER2
   INTERFACE
      SUBROUTINE ZHER2(Uplo,N,Alpha,X,Incx,Y,Incy,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE ZHER2
   END INTERFACE
END MODULE S_ZHER2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHER2K
   INTERFACE
      SUBROUTINE ZHER2K(Uplo,Trans,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE ZHER2K
   END INTERFACE
END MODULE S_ZHER2K
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHER
   INTERFACE
      SUBROUTINE ZHER(Uplo,N,Alpha,X,Incx,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE ZHER
   END INTERFACE
END MODULE S_ZHER
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHERK
   INTERFACE
      SUBROUTINE ZHERK(Uplo,Trans,N,K,Alpha,A,Lda,Beta,C,Ldc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE ZHERK
   END INTERFACE
END MODULE S_ZHERK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPMV
   INTERFACE
      SUBROUTINE ZHPMV(Uplo,N,Alpha,Ap,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE ZHPMV
   END INTERFACE
END MODULE S_ZHPMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPR2
   INTERFACE
      SUBROUTINE ZHPR2(Uplo,N,Alpha,X,Incx,Y,Incy,Ap)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      END SUBROUTINE ZHPR2
   END INTERFACE
END MODULE S_ZHPR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPR
   INTERFACE
      SUBROUTINE ZHPR(Uplo,N,Alpha,X,Incx,Ap)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      END SUBROUTINE ZHPR
   END INTERFACE
END MODULE S_ZHPR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZROTG
   INTERFACE
      SUBROUTINE ZROTG(Ca,Cb,C,S)
      USE F77KINDS                        
      IMPLICIT NONE
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Ca
      COMPLEX(CX16KIND) , INTENT(IN) :: Cb
      REAL(R8KIND) , INTENT(OUT) :: C
      COMPLEX(CX16KIND) , INTENT(OUT) :: S
      END SUBROUTINE ZROTG
   END INTERFACE
END MODULE S_ZROTG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSCAL
   INTERFACE
      SUBROUTINE ZSCAL(N,Za,Zx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Za
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE ZSCAL
   END INTERFACE
END MODULE S_ZSCAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSWAP
   INTERFACE
      SUBROUTINE ZSWAP(N,Zx,Incx,Zy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Zy
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE ZSWAP
   END INTERFACE
END MODULE S_ZSWAP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYMM
   INTERFACE
      SUBROUTINE ZSYMM(Side,Uplo,M,N,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE ZSYMM
   END INTERFACE
END MODULE S_ZSYMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYR2K
   INTERFACE
      SUBROUTINE ZSYR2K(Uplo,Trans,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE ZSYR2K
   END INTERFACE
END MODULE S_ZSYR2K
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYRK
   INTERFACE
      SUBROUTINE ZSYRK(Uplo,Trans,N,K,Alpha,A,Lda,Beta,C,Ldc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      END SUBROUTINE ZSYRK
   END INTERFACE
END MODULE S_ZSYRK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTBMV
   INTERFACE
      SUBROUTINE ZTBMV(Uplo,Trans,Diag,N,K,A,Lda,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE ZTBMV
   END INTERFACE
END MODULE S_ZTBMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTBSV
   INTERFACE
      SUBROUTINE ZTBSV(Uplo,Trans,Diag,N,K,A,Lda,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE ZTBSV
   END INTERFACE
END MODULE S_ZTBSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPMV
   INTERFACE
      SUBROUTINE ZTPMV(Uplo,Trans,Diag,N,Ap,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE ZTPMV
   END INTERFACE
END MODULE S_ZTPMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPSV
   INTERFACE
      SUBROUTINE ZTPSV(Uplo,Trans,Diag,N,Ap,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE ZTPSV
   END INTERFACE
END MODULE S_ZTPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRMM
   INTERFACE
      SUBROUTINE ZTRMM(Side,Uplo,Transa,Diag,M,N,Alpha,A,Lda,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Transa
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE ZTRMM
   END INTERFACE
END MODULE S_ZTRMM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRMV
   INTERFACE
      SUBROUTINE ZTRMV(Uplo,Trans,Diag,N,A,Lda,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE ZTRMV
   END INTERFACE
END MODULE S_ZTRMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRSM
   INTERFACE
      SUBROUTINE ZTRSM(Side,Uplo,Transa,Diag,M,N,Alpha,A,Lda,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Transa
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE ZTRSM
   END INTERFACE
END MODULE S_ZTRSM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRSV
   INTERFACE
      SUBROUTINE ZTRSV(Uplo,Trans,Diag,N,A,Lda,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE ZTRSV
   END INTERFACE
END MODULE S_ZTRSV
