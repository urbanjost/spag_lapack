!*==s_caxpy.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CAXPY
      IMPLICIT NONE
!*--S_CAXPY4
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY)
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      COMPLEX  ::  CA
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  CX
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  CY
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE CAXPY
      END INTERFACE
      END MODULE S_CAXPY
!*==s_ccopy.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CCOPY
      IMPLICIT NONE
!*--S_CCOPY23
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CCOPY(N,CX,INCX,CY,INCY)
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  CX
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(OUT) , DIMENSION(*)  ::  CY
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE CCOPY
      END INTERFACE
      END MODULE S_CCOPY
!*==s_cdotc.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CDOTC
      IMPLICIT NONE
!*--S_CDOTC41
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION CDOTC(N,CX,INCX,CY,INCY)
      IMPLICIT NONE
      COMPLEX  ::  CDOTC
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  CX
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  CY
      INTEGER , INTENT(IN)  ::  INCY
      END FUNCTION CDOTC
      END INTERFACE
      END MODULE S_CDOTC
!*==s_cdotu.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CDOTU
      IMPLICIT NONE
!*--S_CDOTU60
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION CDOTU(N,CX,INCX,CY,INCY)
      IMPLICIT NONE
      COMPLEX  ::  CDOTU
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  CX
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  CY
      INTEGER , INTENT(IN)  ::  INCY
      END FUNCTION CDOTU
      END INTERFACE
      END MODULE S_CDOTU
!*==s_cgbmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CGBMV
      IMPLICIT NONE
!*--S_CGBMV79
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  KL
      INTEGER , INTENT(IN)  ::  KU
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(IN)  ::  BETA
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE CGBMV
      END INTERFACE
      END MODULE S_CGBMV
!*==s_cgemm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CGEMM
      IMPLICIT NONE
!*--S_CGEMM115
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  TRANSA
      CHARACTER  ::  TRANSB
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      COMPLEX , INTENT(IN)  ::  BETA
      COMPLEX , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE CGEMM
      END INTERFACE
      END MODULE S_CGEMM
!*==s_cgemv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CGEMV
      IMPLICIT NONE
!*--S_CGEMV151
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(IN)  ::  BETA
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE CGEMV
      END INTERFACE
      END MODULE S_CGEMV
!*==s_cgerc.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CGERC
      IMPLICIT NONE
!*--S_CGERC185
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      COMPLEX , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE CGERC
      END INTERFACE
      END MODULE S_CGERC
!*==s_cgeru.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CGERU
      IMPLICIT NONE
!*--S_CGERU215
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      COMPLEX , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE CGERU
      END INTERFACE
      END MODULE S_CGERU
!*==s_chbmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CHBMV
      IMPLICIT NONE
!*--S_CHBMV245
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(IN)  ::  BETA
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE CHBMV
      END INTERFACE
      END MODULE S_CHBMV
!*==s_chemm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CHEMM
      IMPLICIT NONE
!*--S_CHEMM279
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      COMPLEX , INTENT(IN)  ::  BETA
      COMPLEX , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE CHEMM
      END INTERFACE
      END MODULE S_CHEMM
!*==s_chemv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CHEMV
      IMPLICIT NONE
!*--S_CHEMV314
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(IN)  ::  BETA
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE CHEMV
      END INTERFACE
      END MODULE S_CHEMV
!*==s_cher2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CHER2
      IMPLICIT NONE
!*--S_CHER2347
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      COMPLEX , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE CHER2
      END INTERFACE
      END MODULE S_CHER2
!*==s_cher2k.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CHER2K
      IMPLICIT NONE
!*--S_CHER2K377
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      REAL , INTENT(IN)  ::  BETA
      COMPLEX , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE CHER2K
      END INTERFACE
      END MODULE S_CHER2K
!*==s_cher.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CHER
      IMPLICIT NONE
!*--S_CHER412
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CHER(UPLO,N,ALPHA,X,INCX,A,LDA)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE CHER
      END INTERFACE
      END MODULE S_CHER
!*==s_cherk.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CHERK
      IMPLICIT NONE
!*--S_CHERK440
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(IN)  ::  BETA
      COMPLEX , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE CHERK
      END INTERFACE
      END MODULE S_CHERK
!*==s_chpmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CHPMV
      IMPLICIT NONE
!*--S_CHPMV471
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CHPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  AP
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(IN)  ::  BETA
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE CHPMV
      END INTERFACE
      END MODULE S_CHPMV
!*==s_chpr2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CHPR2
      IMPLICIT NONE
!*--S_CHPR2503
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  AP
      END SUBROUTINE CHPR2
      END INTERFACE
      END MODULE S_CHPR2
!*==s_chpr.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CHPR
      IMPLICIT NONE
!*--S_CHPR532
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CHPR(UPLO,N,ALPHA,X,INCX,AP)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  AP
      END SUBROUTINE CHPR
      END INTERFACE
      END MODULE S_CHPR
!*==s_crotg.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CROTG
      IMPLICIT NONE
!*--S_CROTG559
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CROTG(CA,CB,C,S)
      IMPLICIT NONE
      COMPLEX , INTENT(INOUT)  ::  CA
      COMPLEX , INTENT(IN)  ::  CB
      REAL , INTENT(OUT)  ::  C
      COMPLEX , INTENT(OUT)  ::  S
      END SUBROUTINE CROTG
      END INTERFACE
      END MODULE S_CROTG
!*==s_cscal.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CSCAL
      IMPLICIT NONE
!*--S_CSCAL576
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CSCAL(N,CA,CX,INCX)
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN)  ::  CA
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  CX
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE CSCAL
      END INTERFACE
      END MODULE S_CSCAL
!*==s_csrot.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CSROT
      IMPLICIT NONE
!*--S_CSROT593
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CSROT(N,CX,INCX,CY,INCY,C,S)
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  CX
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  CY
      INTEGER , INTENT(IN)  ::  INCY
      REAL , INTENT(IN)  ::  C
      REAL , INTENT(IN)  ::  S
      END SUBROUTINE CSROT
      END INTERFACE
      END MODULE S_CSROT
!*==s_csscal.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CSSCAL
      IMPLICIT NONE
!*--S_CSSCAL613
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CSSCAL(N,SA,CX,INCX)
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  SA
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  CX
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE CSSCAL
      END INTERFACE
      END MODULE S_CSSCAL
!*==s_cswap.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CSWAP
      IMPLICIT NONE
!*--S_CSWAP630
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CSWAP(N,CX,INCX,CY,INCY)
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  CX
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  CY
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE CSWAP
      END INTERFACE
      END MODULE S_CSWAP
!*==s_csymm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CSYMM
      IMPLICIT NONE
!*--S_CSYMM648
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      COMPLEX , INTENT(IN)  ::  BETA
      COMPLEX , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE CSYMM
      END INTERFACE
      END MODULE S_CSYMM
!*==s_csyr2k.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CSYR2K
      IMPLICIT NONE
!*--S_CSYR2K683
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      COMPLEX , INTENT(IN)  ::  BETA
      COMPLEX , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE CSYR2K
      END INTERFACE
      END MODULE S_CSYR2K
!*==s_csyrk.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CSYRK
      IMPLICIT NONE
!*--S_CSYRK718
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(IN)  ::  BETA
      COMPLEX , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE CSYRK
      END INTERFACE
      END MODULE S_CSYRK
!*==s_ctbmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CTBMV
      IMPLICIT NONE
!*--S_CTBMV751
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE CTBMV
      END INTERFACE
      END MODULE S_CTBMV
!*==s_ctbsv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CTBSV
      IMPLICIT NONE
!*--S_CTBSV781
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE CTBSV
      END INTERFACE
      END MODULE S_CTBSV
!*==s_ctpmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CTPMV
      IMPLICIT NONE
!*--S_CTPMV811
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  AP
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE CTPMV
      END INTERFACE
      END MODULE S_CTPMV
!*==s_ctpsv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CTPSV
      IMPLICIT NONE
!*--S_CTPSV839
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  AP
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE CTPSV
      END INTERFACE
      END MODULE S_CTPSV
!*==s_ctrmm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CTRMM
      IMPLICIT NONE
!*--S_CTRMM867
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANSA
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(INOUT) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      END SUBROUTINE CTRMM
      END INTERFACE
      END MODULE S_CTRMM
!*==s_ctrmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CTRMV
      IMPLICIT NONE
!*--S_CTRMV901
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE CTRMV
      END INTERFACE
      END MODULE S_CTRMV
!*==s_ctrsm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CTRSM
      IMPLICIT NONE
!*--S_CTRSM930
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANSA
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN)  ::  ALPHA
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(INOUT) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      END SUBROUTINE CTRSM
      END INTERFACE
      END MODULE S_CTRSM
!*==s_ctrsv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_CTRSV
      IMPLICIT NONE
!*--S_CTRSV964
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE CTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE CTRSV
      END INTERFACE
      END MODULE S_CTRSV
!*==s_dasum.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DASUM
      IMPLICIT NONE
!*--S_DASUM993
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION DASUM(N,DX,INCX)
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND)  ::  DASUM
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  DX
      INTEGER , INTENT(IN)  ::  INCX
      END FUNCTION DASUM
      END INTERFACE
      END MODULE S_DASUM
!*==s_daxpy.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DAXPY
      IMPLICIT NONE
!*--S_DAXPY1011
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  DA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  DX
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  DY
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE DAXPY
      END INTERFACE
      END MODULE S_DAXPY
!*==s_dcabs1.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DCABS1
      IMPLICIT NONE
!*--S_DCABS11031
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION DCABS1(Z)
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND)  ::  DCABS1
      COMPLEX(CX16KIND) , INTENT(IN)  ::  Z
      END FUNCTION DCABS1
      END INTERFACE
      END MODULE S_DCABS1
!*==s_dcopy.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DCOPY
      IMPLICIT NONE
!*--S_DCOPY1047
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  DX
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*)  ::  DY
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE DCOPY
      END INTERFACE
      END MODULE S_DCOPY
!*==s_ddot.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DDOT
      IMPLICIT NONE
!*--S_DDOT1066
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION DDOT(N,DX,INCX,DY,INCY)
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND)  ::  DDOT
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  DX
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  DY
      INTEGER , INTENT(IN)  ::  INCY
      END FUNCTION DDOT
      END INTERFACE
      END MODULE S_DDOT
!*==s_dgbmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DGBMV
      IMPLICIT NONE
!*--S_DGBMV1086
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  KL
      INTEGER , INTENT(IN)  ::  KU
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(IN)  ::  BETA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE DGBMV
      END INTERFACE
      END MODULE S_DGBMV
!*==s_dgemm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DGEMM
      IMPLICIT NONE
!*--S_DGEMM1121
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER  ::  TRANSA
      CHARACTER  ::  TRANSB
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      REAL(R8KIND) , INTENT(IN)  ::  BETA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE DGEMM
      END INTERFACE
      END MODULE S_DGEMM
!*==s_dgemv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DGEMV
      IMPLICIT NONE
!*--S_DGEMV1156
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(IN)  ::  BETA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE DGEMV
      END INTERFACE
      END MODULE S_DGEMV
!*==s_dger.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DGER
      IMPLICIT NONE
!*--S_DGER1189
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE DGER
      END INTERFACE
      END MODULE S_DGER
!*==s_dnrm2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DNRM2
      IMPLICIT NONE
!*--S_DNRM21220
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION DNRM2(N,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND)  ::  DNRM2
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END FUNCTION DNRM2
      END INTERFACE
      END MODULE S_DNRM2
!*==s_drot.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DROT
      IMPLICIT NONE
!*--S_DROT1246
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  DX
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  DY
      INTEGER , INTENT(IN)  ::  INCY
      REAL(R8KIND) , INTENT(IN)  ::  C
      REAL(R8KIND) , INTENT(IN)  ::  S
      END SUBROUTINE DROT
      END INTERFACE
      END MODULE S_DROT
!*==s_drotg.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DROTG
      IMPLICIT NONE
!*--S_DROTG1267
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DROTG(DA,DB,C,S)
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) , INTENT(INOUT)  ::  DA
      REAL(R8KIND) , INTENT(INOUT)  ::  DB
      REAL(R8KIND) , INTENT(INOUT)  ::  C
      REAL(R8KIND) , INTENT(INOUT)  ::  S
      END SUBROUTINE DROTG
      END INTERFACE
      END MODULE S_DROTG
!*==s_drotm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DROTM
      IMPLICIT NONE
!*--S_DROTM1285
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DROTM(N,DX,INCX,DY,INCY,DPARAM)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  DX
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  DY
      INTEGER , INTENT(IN)  ::  INCY
      REAL(R8KIND) , INTENT(IN) , DIMENSION(5)  ::  DPARAM
      END SUBROUTINE DROTM
      END INTERFACE
      END MODULE S_DROTM
!*==s_drotmg.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DROTMG
      IMPLICIT NONE
!*--S_DROTMG1305
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND) , INTENT(INOUT)  ::  DD1
      REAL(R8KIND) , INTENT(INOUT)  ::  DD2
      REAL(R8KIND) , INTENT(INOUT)  ::  DX1
      REAL(R8KIND) , INTENT(IN)  ::  DY1
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(5)  ::  DPARAM
      END SUBROUTINE DROTMG
      END INTERFACE
      END MODULE S_DROTMG
!*==s_dsbmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSBMV
      IMPLICIT NONE
!*--S_DSBMV1324
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(IN)  ::  BETA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE DSBMV
      END INTERFACE
      END MODULE S_DSBMV
!*==s_dscal.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSCAL
      IMPLICIT NONE
!*--S_DSCAL1357
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DSCAL(N,DA,DX,INCX)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  DA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  DX
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE DSCAL
      END INTERFACE
      END MODULE S_DSCAL
!*==s_dsdot.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSDOT
      IMPLICIT NONE
!*--S_DSDOT1375
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION DSDOT(N,SX,INCX,SY,INCY)
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND)  ::  DSDOT
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN) , DIMENSION(*)  ::  SX
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(IN) , DIMENSION(*)  ::  SY
      INTEGER , INTENT(IN)  ::  INCY
      END FUNCTION DSDOT
      END INTERFACE
      END MODULE S_DSDOT
!*==s_dspmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSPMV
      IMPLICIT NONE
!*--S_DSPMV1395
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  AP
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(IN)  ::  BETA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE DSPMV
      END INTERFACE
      END MODULE S_DSPMV
!*==s_dspr2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSPR2
      IMPLICIT NONE
!*--S_DSPR21426
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  AP
      END SUBROUTINE DSPR2
      END INTERFACE
      END MODULE S_DSPR2
!*==s_dspr.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSPR
      IMPLICIT NONE
!*--S_DSPR1456
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DSPR(UPLO,N,ALPHA,X,INCX,AP)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  AP
      END SUBROUTINE DSPR
      END INTERFACE
      END MODULE S_DSPR
!*==s_dswap.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSWAP
      IMPLICIT NONE
!*--S_DSWAP1484
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  DX
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  DY
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE DSWAP
      END INTERFACE
      END MODULE S_DSWAP
!*==s_dsymm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSYMM
      IMPLICIT NONE
!*--S_DSYMM1503
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      REAL(R8KIND) , INTENT(IN)  ::  BETA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE DSYMM
      END INTERFACE
      END MODULE S_DSYMM
!*==s_dsymv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSYMV
      IMPLICIT NONE
!*--S_DSYMV1537
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(IN)  ::  BETA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE DSYMV
      END INTERFACE
      END MODULE S_DSYMV
!*==s_dsyr2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSYR2
      IMPLICIT NONE
!*--S_DSYR21569
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE DSYR2
      END INTERFACE
      END MODULE S_DSYR2
!*==s_dsyr2k.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSYR2K
      IMPLICIT NONE
!*--S_DSYR2K1600
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      REAL(R8KIND) , INTENT(IN)  ::  BETA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE DSYR2K
      END INTERFACE
      END MODULE S_DSYR2K
!*==s_dsyr.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSYR
      IMPLICIT NONE
!*--S_DSYR1634
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE DSYR
      END INTERFACE
      END MODULE S_DSYR
!*==s_dsyrk.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DSYRK
      IMPLICIT NONE
!*--S_DSYRK1663
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(IN)  ::  BETA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE DSYRK
      END INTERFACE
      END MODULE S_DSYRK
!*==s_dtbmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DTBMV
      IMPLICIT NONE
!*--S_DTBMV1695
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE DTBMV
      END INTERFACE
      END MODULE S_DTBMV
!*==s_dtbsv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DTBSV
      IMPLICIT NONE
!*--S_DTBSV1726
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE DTBSV
      END INTERFACE
      END MODULE S_DTBSV
!*==s_dtpmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DTPMV
      IMPLICIT NONE
!*--S_DTPMV1757
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  AP
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE DTPMV
      END INTERFACE
      END MODULE S_DTPMV
!*==s_dtpsv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DTPSV
      IMPLICIT NONE
!*--S_DTPSV1786
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  AP
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE DTPSV
      END INTERFACE
      END MODULE S_DTPSV
!*==s_dtrmm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DTRMM
      IMPLICIT NONE
!*--S_DTRMM1815
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANSA
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      END SUBROUTINE DTRMM
      END INTERFACE
      END MODULE S_DTRMM
!*==s_dtrmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DTRMV
      IMPLICIT NONE
!*--S_DTRMV1848
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE DTRMV
      END INTERFACE
      END MODULE S_DTRMV
!*==s_dtrsm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DTRSM
      IMPLICIT NONE
!*--S_DTRSM1878
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANSA
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      END SUBROUTINE DTRSM
      END INTERFACE
      END MODULE S_DTRSM
!*==s_dtrsv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DTRSV
      IMPLICIT NONE
!*--S_DTRSV1911
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE DTRSV
      END INTERFACE
      END MODULE S_DTRSV
!*==s_dzasum.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DZASUM
      IMPLICIT NONE
!*--S_DZASUM1941
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION DZASUM(N,ZX,INCX)
      USE F77KINDS
      IMPLICIT NONE
      REAL(R8KIND)  ::  DZASUM
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , DIMENSION(*)  ::  ZX
      INTEGER , INTENT(IN)  ::  INCX
      END FUNCTION DZASUM
      END INTERFACE
      END MODULE S_DZASUM
!*==s_dznrm2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_DZNRM2
      IMPLICIT NONE
!*--S_DZNRM21959
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION DZNRM2(N,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND)  ::  DZNRM2
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END FUNCTION DZNRM2
      END INTERFACE
      END MODULE S_DZNRM2
!*==s_icamax.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ICAMAX
      IMPLICIT NONE
!*--S_ICAMAX1985
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION ICAMAX(N,CX,INCX)
      IMPLICIT NONE
      INTEGER  ::  ICAMAX
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , DIMENSION(*)  ::  CX
      INTEGER , INTENT(IN)  ::  INCX
      END FUNCTION ICAMAX
      END INTERFACE
      END MODULE S_ICAMAX
!*==s_idamax.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_IDAMAX
      IMPLICIT NONE
!*--S_IDAMAX2002
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION IDAMAX(N,DX,INCX)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER  ::  IDAMAX
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*)  ::  DX
      INTEGER , INTENT(IN)  ::  INCX
      END FUNCTION IDAMAX
      END INTERFACE
      END MODULE S_IDAMAX
!*==s_isamax.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ISAMAX
      IMPLICIT NONE
!*--S_ISAMAX2020
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION ISAMAX(N,SX,INCX)
      IMPLICIT NONE
      INTEGER  ::  ISAMAX
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN) , DIMENSION(*)  ::  SX
      INTEGER , INTENT(IN)  ::  INCX
      END FUNCTION ISAMAX
      END INTERFACE
      END MODULE S_ISAMAX
!*==s_izamax.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_IZAMAX
      IMPLICIT NONE
!*--S_IZAMAX2037
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION IZAMAX(N,ZX,INCX)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER  ::  IZAMAX
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , DIMENSION(*)  ::  ZX
      INTEGER , INTENT(IN)  ::  INCX
      END FUNCTION IZAMAX
      END INTERFACE
      END MODULE S_IZAMAX
!*==s_lsame.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_LSAME
      IMPLICIT NONE
!*--S_LSAME2055
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION LSAME(CA,CB)
      IMPLICIT NONE
      LOGICAL  ::  LSAME
      CHARACTER , INTENT(IN)  ::  CA
      CHARACTER , INTENT(IN)  ::  CB
      END FUNCTION LSAME
      END INTERFACE
      END MODULE S_LSAME
!*==s_sasum.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SASUM
      IMPLICIT NONE
!*--S_SASUM2071
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION SASUM(N,SX,INCX)
      IMPLICIT NONE
      REAL  ::  SASUM
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN) , DIMENSION(*)  ::  SX
      INTEGER , INTENT(IN)  ::  INCX
      END FUNCTION SASUM
      END INTERFACE
      END MODULE S_SASUM
!*==s_saxpy.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SAXPY
      IMPLICIT NONE
!*--S_SAXPY2088
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  SA
      REAL , INTENT(IN) , DIMENSION(*)  ::  SX
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  SY
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE SAXPY
      END INTERFACE
      END MODULE S_SAXPY
!*==s_scabs1.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SCABS1
      IMPLICIT NONE
!*--S_SCABS12107
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION SCABS1(Z)
      IMPLICIT NONE
      REAL  ::  SCABS1
      COMPLEX , INTENT(IN)  ::  Z
      END FUNCTION SCABS1
      END INTERFACE
      END MODULE S_SCABS1
!*==s_scasum.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SCASUM
      IMPLICIT NONE
!*--S_SCASUM2122
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION SCASUM(N,CX,INCX)
      IMPLICIT NONE
      REAL  ::  SCASUM
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  CX
      INTEGER , INTENT(IN)  ::  INCX
      END FUNCTION SCASUM
      END INTERFACE
      END MODULE S_SCASUM
!*==s_scnrm2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SCNRM2
      IMPLICIT NONE
!*--S_SCNRM22139
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION SCNRM2(N,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL  ::  SCNRM2
      INTEGER , INTENT(IN)  ::  N
      COMPLEX , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END FUNCTION SCNRM2
      END INTERFACE
      END MODULE S_SCNRM2
!*==s_scopy.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SCOPY
      IMPLICIT NONE
!*--S_SCOPY2164
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN) , DIMENSION(*)  ::  SX
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(OUT) , DIMENSION(*)  ::  SY
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE SCOPY
      END INTERFACE
      END MODULE S_SCOPY
!*==s_sdot.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SDOT
      IMPLICIT NONE
!*--S_SDOT2182
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION SDOT(N,SX,INCX,SY,INCY)
      IMPLICIT NONE
      REAL  ::  SDOT
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN) , DIMENSION(*)  ::  SX
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(IN) , DIMENSION(*)  ::  SY
      INTEGER , INTENT(IN)  ::  INCY
      END FUNCTION SDOT
      END INTERFACE
      END MODULE S_SDOT
!*==s_sdsdot.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SDSDOT
      IMPLICIT NONE
!*--S_SDSDOT2201
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION SDSDOT(N,SB,SX,INCX,SY,INCY)
      USE F77KINDS
      IMPLICIT NONE
      REAL  ::  SDSDOT
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  SB
      REAL , INTENT(IN) , DIMENSION(*)  ::  SX
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(IN) , DIMENSION(*)  ::  SY
      INTEGER , INTENT(IN)  ::  INCY
      END FUNCTION SDSDOT
      END INTERFACE
      END MODULE S_SDSDOT
!*==s_sgbmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SGBMV
      IMPLICIT NONE
!*--S_SGBMV2222
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  KL
      INTEGER , INTENT(IN)  ::  KU
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(IN)  ::  BETA
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE SGBMV
      END INTERFACE
      END MODULE S_SGBMV
!*==s_sgemm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SGEMM
      IMPLICIT NONE
!*--S_SGEMM2256
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER  ::  TRANSA
      CHARACTER  ::  TRANSB
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      REAL , INTENT(IN)  ::  BETA
      REAL , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE SGEMM
      END INTERFACE
      END MODULE S_SGEMM
!*==s_sgemv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SGEMV
      IMPLICIT NONE
!*--S_SGEMV2290
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(IN)  ::  BETA
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE SGEMV
      END INTERFACE
      END MODULE S_SGEMV
!*==s_sger.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SGER
      IMPLICIT NONE
!*--S_SGER2322
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      REAL , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE SGER
      END INTERFACE
      END MODULE S_SGER
!*==s_snrm2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SNRM2
      IMPLICIT NONE
!*--S_SNRM22352
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION SNRM2(N,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL  ::  SNRM2
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END FUNCTION SNRM2
      END INTERFACE
      END MODULE S_SNRM2
!*==s_srot.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SROT
      IMPLICIT NONE
!*--S_SROT2377
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SROT(N,SX,INCX,SY,INCY,C,S)
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  SX
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  SY
      INTEGER , INTENT(IN)  ::  INCY
      REAL , INTENT(IN)  ::  C
      REAL , INTENT(IN)  ::  S
      END SUBROUTINE SROT
      END INTERFACE
      END MODULE S_SROT
!*==s_srotg.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SROTG
      IMPLICIT NONE
!*--S_SROTG2397
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SROTG(SA,SB,C,S)
      IMPLICIT NONE
      REAL , INTENT(INOUT)  ::  SA
      REAL , INTENT(INOUT)  ::  SB
      REAL , INTENT(INOUT)  ::  C
      REAL , INTENT(INOUT)  ::  S
      END SUBROUTINE SROTG
      END INTERFACE
      END MODULE S_SROTG
!*==s_srotm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SROTM
      IMPLICIT NONE
!*--S_SROTM2414
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SROTM(N,SX,INCX,SY,INCY,SPARAM)
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  SX
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  SY
      INTEGER , INTENT(IN)  ::  INCY
      REAL , INTENT(IN) , DIMENSION(5)  ::  SPARAM
      END SUBROUTINE SROTM
      END INTERFACE
      END MODULE S_SROTM
!*==s_srotmg.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SROTMG
      IMPLICIT NONE
!*--S_SROTMG2433
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SROTMG(SD1,SD2,SX1,SY1,SPARAM)
      IMPLICIT NONE
      REAL , INTENT(INOUT)  ::  SD1
      REAL , INTENT(INOUT)  ::  SD2
      REAL , INTENT(INOUT)  ::  SX1
      REAL , INTENT(IN)  ::  SY1
      REAL , INTENT(OUT) , DIMENSION(5)  ::  SPARAM
      END SUBROUTINE SROTMG
      END INTERFACE
      END MODULE S_SROTMG
!*==s_ssbmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SSBMV
      IMPLICIT NONE
!*--S_SSBMV2451
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(IN)  ::  BETA
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE SSBMV
      END INTERFACE
      END MODULE S_SSBMV
!*==s_sscal.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SSCAL
      IMPLICIT NONE
!*--S_SSCAL2483
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SSCAL(N,SA,SX,INCX)
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  SA
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  SX
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE SSCAL
      END INTERFACE
      END MODULE S_SSCAL
!*==s_sspmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SSPMV
      IMPLICIT NONE
!*--S_SSPMV2500
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(*)  ::  AP
      REAL , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(IN)  ::  BETA
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE SSPMV
      END INTERFACE
      END MODULE S_SSPMV
!*==s_sspr2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SSPR2
      IMPLICIT NONE
!*--S_SSPR22530
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  AP
      END SUBROUTINE SSPR2
      END INTERFACE
      END MODULE S_SSPR2
!*==s_sspr.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SSPR
      IMPLICIT NONE
!*--S_SSPR2559
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SSPR(UPLO,N,ALPHA,X,INCX,AP)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  AP
      END SUBROUTINE SSPR
      END INTERFACE
      END MODULE S_SSPR
!*==s_sswap.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SSWAP
      IMPLICIT NONE
!*--S_SSWAP2586
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  SX
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  SY
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE SSWAP
      END INTERFACE
      END MODULE S_SSWAP
!*==s_ssymm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SSYMM
      IMPLICIT NONE
!*--S_SSYMM2604
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      REAL , INTENT(IN)  ::  BETA
      REAL , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE SSYMM
      END INTERFACE
      END MODULE S_SSYMM
!*==s_ssymv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SSYMV
      IMPLICIT NONE
!*--S_SSYMV2637
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(IN)  ::  BETA
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE SSYMV
      END INTERFACE
      END MODULE S_SSYMV
!*==s_ssyr2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SSYR2
      IMPLICIT NONE
!*--S_SSYR22668
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      REAL , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE SSYR2
      END INTERFACE
      END MODULE S_SSYR2
!*==s_ssyr2k.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SSYR2K
      IMPLICIT NONE
!*--S_SSYR2K2698
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      REAL , INTENT(IN)  ::  BETA
      REAL , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE SSYR2K
      END INTERFACE
      END MODULE S_SSYR2K
!*==s_ssyr.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SSYR
      IMPLICIT NONE
!*--S_SSYR2731
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      REAL , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE SSYR
      END INTERFACE
      END MODULE S_SSYR
!*==s_ssyrk.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_SSYRK
      IMPLICIT NONE
!*--S_SSYRK2759
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE SSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(IN)  ::  BETA
      REAL , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE SSYRK
      END INTERFACE
      END MODULE S_SSYRK
!*==s_stbmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_STBMV
      IMPLICIT NONE
!*--S_STBMV2790
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE STBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE STBMV
      END INTERFACE
      END MODULE S_STBMV
!*==s_stbsv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_STBSV
      IMPLICIT NONE
!*--S_STBSV2820
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE STBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE STBSV
      END INTERFACE
      END MODULE S_STBSV
!*==s_stpmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_STPMV
      IMPLICIT NONE
!*--S_STPMV2850
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE STPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN) , DIMENSION(*)  ::  AP
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE STPMV
      END INTERFACE
      END MODULE S_STPMV
!*==s_stpsv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_STPSV
      IMPLICIT NONE
!*--S_STPSV2878
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE STPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN) , DIMENSION(*)  ::  AP
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE STPSV
      END INTERFACE
      END MODULE S_STPSV
!*==s_strmm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_STRMM
      IMPLICIT NONE
!*--S_STRMM2906
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANSA
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(INOUT) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      END SUBROUTINE STRMM
      END INTERFACE
      END MODULE S_STRMM
!*==s_strmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_STRMV
      IMPLICIT NONE
!*--S_STRMV2938
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE STRMV
      END INTERFACE
      END MODULE S_STRMV
!*==s_strsm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_STRSM
      IMPLICIT NONE
!*--S_STRSM2967
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANSA
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN)  ::  ALPHA
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(INOUT) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      END SUBROUTINE STRSM
      END INTERFACE
      END MODULE S_STRSM
!*==s_strsv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_STRSV
      IMPLICIT NONE
!*--S_STRSV2999
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE STRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      REAL , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE STRSV
      END INTERFACE
      END MODULE S_STRSV
!*==s_xerbla_array.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_XERBLA_ARRAY
      IMPLICIT NONE
!*--S_XERBLA_ARRAY3028
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE XERBLA_ARRAY(SRNAME_ARRAY,SRNAME_LEN,INFO)
      IMPLICIT NONE
      CHARACTER(1) , INTENT(IN) , DIMENSION(SRNAME_LEN)                 &
     &                       ::  SRNAME_ARRAY
      INTEGER , INTENT(IN)  ::  SRNAME_LEN
      INTEGER  ::  INFO
      END SUBROUTINE XERBLA_ARRAY
      END INTERFACE
      END MODULE S_XERBLA_ARRAY
!*==s_xerbla.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_XERBLA
      IMPLICIT NONE
!*--S_XERBLA3045
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE XERBLA(SRNAME,INFO)
      IMPLICIT NONE
      CHARACTER(*) , INTENT(IN)  ::  SRNAME
      INTEGER , INTENT(IN)  ::  INFO
      END SUBROUTINE XERBLA
      END INTERFACE
      END MODULE S_XERBLA
!*==s_zaxpy.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZAXPY
      IMPLICIT NONE
!*--S_ZAXPY3060
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND)  ::  ZA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  ZX
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  ZY
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE ZAXPY
      END INTERFACE
      END MODULE S_ZAXPY
!*==s_zcopy.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZCOPY
      IMPLICIT NONE
!*--S_ZCOPY3080
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  ZX
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*)  ::  ZY
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE ZCOPY
      END INTERFACE
      END MODULE S_ZCOPY
!*==s_zdotc.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZDOTC
      IMPLICIT NONE
!*--S_ZDOTC3099
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY)
      USE F77KINDS
      IMPLICIT NONE
      COMPLEX(CX16KIND)  ::  ZDOTC
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  ZX
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  ZY
      INTEGER , INTENT(IN)  ::  INCY
      END FUNCTION ZDOTC
      END INTERFACE
      END MODULE S_ZDOTC
!*==s_zdotu.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZDOTU
      IMPLICIT NONE
!*--S_ZDOTU3119
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)
      USE F77KINDS
      IMPLICIT NONE
      COMPLEX(CX16KIND)  ::  ZDOTU
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  ZX
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  ZY
      INTEGER , INTENT(IN)  ::  INCY
      END FUNCTION ZDOTU
      END INTERFACE
      END MODULE S_ZDOTU
!*==s_zdrot.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZDROT
      IMPLICIT NONE
!*--S_ZDROT3139
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZDROT(N,ZX,INCX,ZY,INCY,C,S)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  ZX
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  ZY
      INTEGER , INTENT(IN)  ::  INCY
      REAL(R8KIND) , INTENT(IN)  ::  C
      REAL(R8KIND) , INTENT(IN)  ::  S
      END SUBROUTINE ZDROT
      END INTERFACE
      END MODULE S_ZDROT
!*==s_zdscal.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZDSCAL
      IMPLICIT NONE
!*--S_ZDSCAL3160
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZDSCAL(N,DA,ZX,INCX)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  DA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  ZX
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE ZDSCAL
      END INTERFACE
      END MODULE S_ZDSCAL
!*==s_zgbmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZGBMV
      IMPLICIT NONE
!*--S_ZGBMV3178
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &        ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  KL
      INTEGER , INTENT(IN)  ::  KU
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(IN)  ::  BETA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE ZGBMV
      END INTERFACE
      END MODULE S_ZGBMV
!*==s_zgemm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZGEMM
      IMPLICIT NONE
!*--S_ZGEMM3215
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &        ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  TRANSA
      CHARACTER  ::  TRANSB
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      COMPLEX(CX16KIND) , INTENT(IN)  ::  BETA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE ZGEMM
      END INTERFACE
      END MODULE S_ZGEMM
!*==s_zgemv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZGEMV
      IMPLICIT NONE
!*--S_ZGEMV3252
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &        ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(IN)  ::  BETA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE ZGEMV
      END INTERFACE
      END MODULE S_ZGEMV
!*==s_zgerc.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZGERC
      IMPLICIT NONE
!*--S_ZGERC3287
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE ZGERC
      END INTERFACE
      END MODULE S_ZGERC
!*==s_zgeru.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZGERU
      IMPLICIT NONE
!*--S_ZGERU3318
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE ZGERU
      END INTERFACE
      END MODULE S_ZGERU
!*==s_zhbmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZHBMV
      IMPLICIT NONE
!*--S_ZHBMV3349
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &        ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(IN)  ::  BETA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE ZHBMV
      END INTERFACE
      END MODULE S_ZHBMV
!*==s_zhemm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZHEMM
      IMPLICIT NONE
!*--S_ZHEMM3384
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &        ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      COMPLEX(CX16KIND) , INTENT(IN)  ::  BETA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE ZHEMM
      END INTERFACE
      END MODULE S_ZHEMM
!*==s_zhemv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZHEMV
      IMPLICIT NONE
!*--S_ZHEMV3420
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &        ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(IN)  ::  BETA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE ZHEMV
      END INTERFACE
      END MODULE S_ZHEMV
!*==s_zher2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZHER2
      IMPLICIT NONE
!*--S_ZHER23454
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE ZHER2
      END INTERFACE
      END MODULE S_ZHER2
!*==s_zher2k.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZHER2K
      IMPLICIT NONE
!*--S_ZHER2K3485
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      REAL(R8KIND) , INTENT(IN)  ::  BETA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE ZHER2K
      END INTERFACE
      END MODULE S_ZHER2K
!*==s_zher.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZHER
      IMPLICIT NONE
!*--S_ZHER3521
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      END SUBROUTINE ZHER
      END INTERFACE
      END MODULE S_ZHER
!*==s_zherk.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZHERK
      IMPLICIT NONE
!*--S_ZHERK3550
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      REAL(R8KIND) , INTENT(IN)  ::  BETA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE ZHERK
      END INTERFACE
      END MODULE S_ZHERK
!*==s_zhpmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZHPMV
      IMPLICIT NONE
!*--S_ZHPMV3582
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZHPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &        ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  AP
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(IN)  ::  BETA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE ZHPMV
      END INTERFACE
      END MODULE S_ZHPMV
!*==s_zhpr2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZHPR2
      IMPLICIT NONE
!*--S_ZHPR23615
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  Y
      INTEGER , INTENT(IN)  ::  INCY
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  AP
      END SUBROUTINE ZHPR2
      END INTERFACE
      END MODULE S_ZHPR2
!*==s_zhpr.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZHPR
      IMPLICIT NONE
!*--S_ZHPR3645
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZHPR(UPLO,N,ALPHA,X,INCX,AP)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  N
      REAL(R8KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  AP
      END SUBROUTINE ZHPR
      END INTERFACE
      END MODULE S_ZHPR
!*==s_zrotg.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZROTG
      IMPLICIT NONE
!*--S_ZROTG3673
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZROTG(CA,CB,C,S)
      USE F77KINDS
      IMPLICIT NONE
      COMPLEX(CX16KIND) , INTENT(INOUT)  ::  CA
      COMPLEX(CX16KIND) , INTENT(IN)  ::  CB
      REAL(R8KIND) , INTENT(OUT)  ::  C
      COMPLEX(CX16KIND) , INTENT(OUT)  ::  S
      END SUBROUTINE ZROTG
      END INTERFACE
      END MODULE S_ZROTG
!*==s_zscal.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZSCAL
      IMPLICIT NONE
!*--S_ZSCAL3691
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZSCAL(N,ZA,ZX,INCX)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ZA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  ZX
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE ZSCAL
      END INTERFACE
      END MODULE S_ZSCAL
!*==s_zswap.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZSWAP
      IMPLICIT NONE
!*--S_ZSWAP3709
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY)
      USE F77KINDS
      IMPLICIT NONE
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  ZX
      INTEGER , INTENT(IN)  ::  INCX
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  ZY
      INTEGER , INTENT(IN)  ::  INCY
      END SUBROUTINE ZSWAP
      END INTERFACE
      END MODULE S_ZSWAP
!*==s_zsymm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZSYMM
      IMPLICIT NONE
!*--S_ZSYMM3728
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &        ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      COMPLEX(CX16KIND) , INTENT(IN)  ::  BETA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE ZSYMM
      END INTERFACE
      END MODULE S_ZSYMM
!*==s_zsyr2k.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZSYR2K
      IMPLICIT NONE
!*--S_ZSYR2K3764
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &        ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      COMPLEX(CX16KIND) , INTENT(IN)  ::  BETA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE ZSYR2K
      END INTERFACE
      END MODULE S_ZSYR2K
!*==s_zsyrk.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZSYRK
      IMPLICIT NONE
!*--S_ZSYRK3800
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &        ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(IN)  ::  BETA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDC,*)  ::  C
      INTEGER , INTENT(IN)  ::  LDC
      END SUBROUTINE ZSYRK
      END INTERFACE
      END MODULE S_ZSYRK
!*==s_ztbmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZTBMV
      IMPLICIT NONE
!*--S_ZTBMV3834
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE ZTBMV
      END INTERFACE
      END MODULE S_ZTBMV
!*==s_ztbsv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZTBSV
      IMPLICIT NONE
!*--S_ZTBSV3865
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      INTEGER , INTENT(IN)  ::  K
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE ZTBSV
      END INTERFACE
      END MODULE S_ZTBSV
!*==s_ztpmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZTPMV
      IMPLICIT NONE
!*--S_ZTPMV3896
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  AP
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE ZTPMV
      END INTERFACE
      END MODULE S_ZTPMV
!*==s_ztpsv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZTPSV
      IMPLICIT NONE
!*--S_ZTPSV3925
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*)  ::  AP
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE ZTPSV
      END INTERFACE
      END MODULE S_ZTPSV
!*==s_ztrmm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZTRMM
      IMPLICIT NONE
!*--S_ZTRMM3954
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &        ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANSA
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      END SUBROUTINE ZTRMM
      END INTERFACE
      END MODULE S_ZTRMM
!*==s_ztrmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZTRMV
      IMPLICIT NONE
!*--S_ZTRMV3989
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE ZTRMV
      END INTERFACE
      END MODULE S_ZTRMV
!*==s_ztrsm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZTRSM
      IMPLICIT NONE
!*--S_ZTRSM4019
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &        ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  SIDE
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANSA
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  M
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN)  ::  ALPHA
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(LDB,*)  ::  B
      INTEGER , INTENT(IN)  ::  LDB
      END SUBROUTINE ZTRSM
      END INTERFACE
      END MODULE S_ZTRSM
!*==s_ztrsv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
      MODULE S_ZTRSV
      IMPLICIT NONE
!*--S_ZTRSV4054
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! End of declarations rewritten by SPAG
!
      INTERFACE
      SUBROUTINE ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      USE F77KINDS
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER  ::  UPLO
      CHARACTER  ::  TRANS
      CHARACTER  ::  DIAG
      INTEGER , INTENT(IN)  ::  N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(LDA,*)  ::  A
      INTEGER , INTENT(IN)  ::  LDA
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*)  ::  X
      INTEGER , INTENT(IN)  ::  INCX
      END SUBROUTINE ZTRSV
      END INTERFACE
      END MODULE S_ZTRSV
