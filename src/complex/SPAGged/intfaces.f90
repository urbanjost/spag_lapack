!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CBBCSD
   INTERFACE
      SUBROUTINE CBBCSD(Jobu1,Jobu2,Jobv1t,Jobv2t,Trans,M,P,Q,Theta,Phi,&
     &                  U1,Ldu1,U2,Ldu2,V1t,Ldv1t,V2t,Ldv2t,B11d,B11e,  &
     &                  B12d,B12e,B21d,B21e,B22d,B22e,Rwork,Lrwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXITR = 6
      REAL , PARAMETER  ::  HUNDRED = 100.0E0 , MEIGHTH = -0.125E0 ,    &
     &                      ONE = 1.0E0 , TEN = 10.0E0 , ZERO = 0.0E0
      COMPLEX , PARAMETER  ::  NEGONECOMPLEX = (-1.0E0,0.0E0)
      REAL , PARAMETER  ::  PIOVER2 =                                   &
     &                      1.57079632679489661923132169163975144210E0
      CHARACTER :: Jobu1
      CHARACTER :: Jobu2
      CHARACTER :: Jobv1t
      CHARACTER :: Jobv2t
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER :: P
      INTEGER :: Q
      REAL , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Phi
      COMPLEX , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      COMPLEX , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      COMPLEX , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      COMPLEX , DIMENSION(Ldv2t,*) :: V2t
      INTEGER :: Ldv2t
      REAL , INTENT(INOUT) , DIMENSION(*) :: B11d
      REAL , INTENT(INOUT) , DIMENSION(*) :: B11e
      REAL , INTENT(INOUT) , DIMENSION(*) :: B12d
      REAL , INTENT(INOUT) , DIMENSION(*) :: B12e
      REAL , INTENT(INOUT) , DIMENSION(*) :: B21d
      REAL , INTENT(INOUT) , DIMENSION(*) :: B21e
      REAL , INTENT(INOUT) , DIMENSION(*) :: B22d
      REAL , INTENT(INOUT) , DIMENSION(*) :: B22e
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CBBCSD
   END INTERFACE
END MODULE S_CBBCSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CBDSQR
   INTERFACE
      SUBROUTINE CBDSQR(Uplo,N,Ncvt,Nru,Ncc,D,E,Vt,Ldvt,U,Ldu,C,Ldc,    &
     &                  Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 ,                &
     &                      NEGONE = -1.0E0 , HNDRTH = 0.01E0 ,         &
     &                      TEN = 10.0E0 , HNDRD = 100.0E0 ,            &
     &                      MEIGTH = -0.125E0
      INTEGER , PARAMETER  ::  MAXITR = 6
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ncvt
      INTEGER :: Nru
      INTEGER :: Ncc
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      COMPLEX , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CBDSQR
   END INTERFACE
END MODULE S_CBDSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBBRD
   INTERFACE
      SUBROUTINE CGBBRD(Vect,M,N,Ncc,Kl,Ku,Ab,Ldab,D,E,Q,Ldq,Pt,Ldpt,C, &
     &                  Ldc,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Vect
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Ncc
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(OUT) , DIMENSION(*) :: D
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(Ldpt,*) :: Pt
      INTEGER :: Ldpt
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGBBRD
   END INTERFACE
END MODULE S_CGBBRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBCON
   INTERFACE
      SUBROUTINE CGBCON(Norm,N,Kl,Ku,Ab,Ldab,Ipiv,Anorm,Rcond,Work,     &
     &                  Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Norm
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGBCON
   END INTERFACE
END MODULE S_CGBCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBEQUB
   INTERFACE
      SUBROUTINE CGBEQUB(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: R
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      REAL , INTENT(OUT) :: Rowcnd
      REAL , INTENT(OUT) :: Colcnd
      REAL , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGBEQUB
   END INTERFACE
END MODULE S_CGBEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBEQU
   INTERFACE
      SUBROUTINE CGBEQU(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: R
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      REAL , INTENT(OUT) :: Rowcnd
      REAL , INTENT(OUT) :: Colcnd
      REAL , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGBEQU
   END INTERFACE
END MODULE S_CGBEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBRFS
   INTERFACE
      SUBROUTINE CGBRFS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv,B,Ldb,&
     &                  X,Ldx,Ferr,Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGBRFS
   END INTERFACE
END MODULE S_CGBRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBRFSX
   INTERFACE
      SUBROUTINE CGBRFSX(Trans,Equed,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,    &
     &                   Ipiv,R,C,B,Ldb,X,Ldx,Rcond,Berr,N_err_bnds,    &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ITREF_DEFAULT = 1.0 ,       &
     &                      ITHRESH_DEFAULT = 10.0 ,                    &
     &                      COMPONENTWISE_DEFAULT = 1.0 ,               &
     &                      RTHRESH_DEFAULT = 0.5 ,                     &
     &                      DZTHRESH_DEFAULT = 0.25
      INTEGER , PARAMETER  ::  LA_LINRX_ITREF_I = 1 ,                   &
     &                         LA_LINRX_ITHRESH_I = 2 ,                 &
     &                         LA_LINRX_CWISE_I = 3 ,                   &
     &                         LA_LINRX_TRUST_I = 1 ,                   &
     &                         LA_LINRX_ERR_I = 2 , LA_LINRX_RCOND_I = 3
      CHARACTER :: Trans
      CHARACTER :: Equed
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: R
      REAL , DIMENSION(*) :: C
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL , INTENT(INOUT) , DIMENSION(*) :: Params
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGBRFSX
   END INTERFACE
END MODULE S_CGBRFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBSV
   INTERFACE
      SUBROUTINE CGBSV(N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGBSV
   END INTERFACE
END MODULE S_CGBSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBSVX
   INTERFACE
      SUBROUTINE CGBSVX(Fact,Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv, &
     &                  Equed,R,C,B,Ldb,X,Ldx,Rcond,Ferr,Berr,Work,     &
     &                  Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: R
      REAL , DIMENSION(*) :: C
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGBSVX
   END INTERFACE
END MODULE S_CGBSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBSVXX
   INTERFACE
      SUBROUTINE CGBSVXX(Fact,Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv,&
     &                   Equed,R,C,B,Ldb,X,Ldx,Rcond,Rpvgrw,Berr,       &
     &                   N_err_bnds,Err_bnds_norm,Err_bnds_comp,Nparams,&
     &                   Params,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL , INTENT(INOUT) , DIMENSION(*) :: R
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , INTENT(OUT) :: Rpvgrw
      REAL , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL , DIMENSION(*) :: Params
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGBSVXX
   END INTERFACE
END MODULE S_CGBSVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBTF2
   INTERFACE
      SUBROUTINE CGBTF2(M,N,Kl,Ku,Ab,Ldab,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGBTF2
   END INTERFACE
END MODULE S_CGBTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBTRF
   INTERFACE
      SUBROUTINE CGBTRF(M,N,Kl,Ku,Ab,Ldab,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDWORK = NBMAX + 1
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGBTRF
   END INTERFACE
END MODULE S_CGBTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGBTRS
   INTERFACE
      SUBROUTINE CGBTRS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGBTRS
   END INTERFACE
END MODULE S_CGBTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEBAK
   INTERFACE
      SUBROUTINE CGEBAK(Job,Side,N,Ilo,Ihi,Scale,M,V,Ldv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Job
      CHARACTER :: Side
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL , INTENT(IN) , DIMENSION(*) :: Scale
      INTEGER :: M
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEBAK
   END INTERFACE
END MODULE S_CGEBAK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEBAL
   INTERFACE
      SUBROUTINE CGEBAL(Job,N,A,Lda,Ilo,Ihi,Scale,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      SCLFAC = 2.0E+0 , FACTOR = 0.95E+0
      CHARACTER :: Job
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) :: Ilo
      INTEGER , INTENT(OUT) :: Ihi
      REAL , INTENT(INOUT) , DIMENSION(*) :: Scale
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEBAL
   END INTERFACE
END MODULE S_CGEBAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEBD2
   INTERFACE
      SUBROUTINE CGEBD2(M,N,A,Lda,D,E,Tauq,Taup,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Tauq
      COMPLEX , DIMENSION(*) :: Taup
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEBD2
   END INTERFACE
END MODULE S_CGEBD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEBRD
   INTERFACE
      SUBROUTINE CGEBRD(M,N,A,Lda,D,E,Tauq,Taup,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      COMPLEX , DIMENSION(*) :: Tauq
      COMPLEX , DIMENSION(*) :: Taup
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEBRD
   END INTERFACE
END MODULE S_CGEBRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGECON
   INTERFACE
      SUBROUTINE CGECON(Norm,N,A,Lda,Anorm,Rcond,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Norm
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGECON
   END INTERFACE
END MODULE S_CGECON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEEQUB
   INTERFACE
      SUBROUTINE CGEEQUB(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: R
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      REAL , INTENT(OUT) :: Rowcnd
      REAL , INTENT(OUT) :: Colcnd
      REAL , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEEQUB
   END INTERFACE
END MODULE S_CGEEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEEQU
   INTERFACE
      SUBROUTINE CGEEQU(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: R
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      REAL , INTENT(OUT) :: Rowcnd
      REAL , INTENT(OUT) :: Colcnd
      REAL , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEEQU
   END INTERFACE
END MODULE S_CGEEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEES
   INTERFACE
      SUBROUTINE CGEES(Jobvs,Sort,SELECT,N,A,Lda,Sdim,W,Vs,Ldvs,Work,   &
     &                 Lwork,Rwork,Bwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobvs
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELECT
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER :: Sdim
      COMPLEX , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldvs,*) :: Vs
      INTEGER :: Ldvs
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEES
   END INTERFACE
END MODULE S_CGEES
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEESX
   INTERFACE
      SUBROUTINE CGEESX(Jobvs,Sort,SELECT,Sense,N,A,Lda,Sdim,W,Vs,Ldvs, &
     &                  Rconde,Rcondv,Work,Lwork,Rwork,Bwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobvs
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELECT
      CHARACTER :: Sense
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Sdim
      COMPLEX , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldvs,*) :: Vs
      INTEGER :: Ldvs
      REAL :: Rconde
      REAL , INTENT(INOUT) :: Rcondv
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEESX
   END INTERFACE
END MODULE S_CGEESX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEEV
   INTERFACE
      SUBROUTINE CGEEV(Jobvl,Jobvr,N,A,Lda,W,Vl,Ldvl,Vr,Ldvr,Work,Lwork,&
     &                 Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: W
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEEV
   END INTERFACE
END MODULE S_CGEEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEEVX
   INTERFACE
      SUBROUTINE CGEEVX(Balanc,Jobvl,Jobvr,Sense,N,A,Lda,W,Vl,Ldvl,Vr,  &
     &                  Ldvr,Ilo,Ihi,Scale,Abnrm,Rconde,Rcondv,Work,    &
     &                  Lwork,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Balanc
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      CHARACTER :: Sense
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: W
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL , DIMENSION(*) :: Scale
      REAL , INTENT(INOUT) :: Abnrm
      REAL , DIMENSION(*) :: Rconde
      REAL , DIMENSION(*) :: Rcondv
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEEVX
   END INTERFACE
END MODULE S_CGEEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEHD2
   INTERFACE
      SUBROUTINE CGEHD2(N,Ilo,Ihi,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER :: Ihi
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEHD2
   END INTERFACE
END MODULE S_CGEHD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEHRD
   INTERFACE
      SUBROUTINE CGEHRD(N,Ilo,Ihi,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDT = NBMAX + 1 ,           &
     &                         TSIZE = LDT*NBMAX
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEHRD
   END INTERFACE
END MODULE S_CGEHRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEJSV
   INTERFACE
      SUBROUTINE CGEJSV(Joba,Jobu,Jobv,Jobr,Jobt,Jobp,M,N,A,Lda,Sva,U,  &
     &                  Ldu,V,Ldv,Cwork,Lwork,Rwork,Lrwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      CHARACTER(1) :: Joba
      CHARACTER(1) :: Jobu
      CHARACTER(1) :: Jobv
      CHARACTER(1) :: Jobr
      CHARACTER(1) :: Jobt
      CHARACTER(1) :: Jobp
      INTEGER :: M
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(N) :: Sva
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , INTENT(INOUT) , DIMENSION(Lwork) :: Cwork
      INTEGER :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(Lrwork) :: Rwork
      INTEGER :: Lrwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEJSV
   END INTERFACE
END MODULE S_CGEJSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGELQ2
   INTERFACE
      SUBROUTINE CGELQ2(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGELQ2
   END INTERFACE
END MODULE S_CGELQ2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGELQ
   INTERFACE
      SUBROUTINE CGELQ(M,N,A,Lda,T,Tsize,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGELQ
   END INTERFACE
END MODULE S_CGELQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGELQF
   INTERFACE
      SUBROUTINE CGELQF(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGELQF
   END INTERFACE
END MODULE S_CGELQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGELQT3
   INTERFACE
      RECURSIVE SUBROUTINE CGELQT3(M,N,A,Lda,T,Ldt,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+00,0.0E+00) ,                &
     &                         ZERO = (0.0E+00,0.0E+00)
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGELQT3
   END INTERFACE
END MODULE S_CGELQT3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGELQT
   INTERFACE
      SUBROUTINE CGELQT(M,N,Mb,A,Lda,T,Ldt,Work,Info)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Mb
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGELQT
   END INTERFACE
END MODULE S_CGELQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGELSD
   INTERFACE
      SUBROUTINE CGELSD(M,N,Nrhs,A,Lda,B,Ldb,S,Rcond,Rank,Work,Lwork,   &
     &                  Rwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: S
      REAL :: Rcond
      INTEGER :: Rank
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGELSD
   END INTERFACE
END MODULE S_CGELSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGELS
   INTERFACE
      SUBROUTINE CGELS(Trans,M,N,Nrhs,A,Lda,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGELS
   END INTERFACE
END MODULE S_CGELS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGELSS
   INTERFACE
      SUBROUTINE CGELSS(M,N,Nrhs,A,Lda,B,Ldb,S,Rcond,Rank,Work,Lwork,   &
     &                  Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGELSS
   END INTERFACE
END MODULE S_CGELSS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGELSY
   INTERFACE
      SUBROUTINE CGELSY(M,N,Nrhs,A,Lda,B,Ldb,Jpvt,Rcond,Rank,Work,Lwork,&
     &                  Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  IMAX = 1 , IMIN = 2
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
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
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGELSY
   END INTERFACE
END MODULE S_CGELSY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEMLQ
   INTERFACE
      SUBROUTINE CGEMLQ(Side,Trans,M,N,K,A,Lda,T,Tsize,C,Ldc,Work,Lwork,&
     &                  Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEMLQ
   END INTERFACE
END MODULE S_CGEMLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEMLQT
   INTERFACE
      SUBROUTINE CGEMLQT(Side,Trans,M,N,K,Mb,V,Ldv,T,Ldt,C,Ldc,Work,    &
     &                   Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Mb
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEMLQT
   END INTERFACE
END MODULE S_CGEMLQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEMQR
   INTERFACE
      SUBROUTINE CGEMQR(Side,Trans,M,N,K,A,Lda,T,Tsize,C,Ldc,Work,Lwork,&
     &                  Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEMQR
   END INTERFACE
END MODULE S_CGEMQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEMQRT
   INTERFACE
      SUBROUTINE CGEMQRT(Side,Trans,M,N,K,Nb,V,Ldv,T,Ldt,C,Ldc,Work,    &
     &                   Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEMQRT
   END INTERFACE
END MODULE S_CGEMQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEQL2
   INTERFACE
      SUBROUTINE CGEQL2(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEQL2
   END INTERFACE
END MODULE S_CGEQL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEQLF
   INTERFACE
      SUBROUTINE CGEQLF(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEQLF
   END INTERFACE
END MODULE S_CGEQLF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEQP3
   INTERFACE
      SUBROUTINE CGEQP3(M,N,A,Lda,Jpvt,Tau,Work,Lwork,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  INB = 1 , INBMIN = 2 , IXOVER = 3
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEQP3
   END INTERFACE
END MODULE S_CGEQP3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEQR2
   INTERFACE
      SUBROUTINE CGEQR2(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEQR2
   END INTERFACE
END MODULE S_CGEQR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEQR2P
   INTERFACE
      SUBROUTINE CGEQR2P(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEQR2P
   END INTERFACE
END MODULE S_CGEQR2P
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEQR
   INTERFACE
      SUBROUTINE CGEQR(M,N,A,Lda,T,Tsize,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEQR
   END INTERFACE
END MODULE S_CGEQR
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
MODULE S_CGEQRFP
   INTERFACE
      SUBROUTINE CGEQRFP(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEQRFP
   END INTERFACE
END MODULE S_CGEQRFP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEQRT2
   INTERFACE
      SUBROUTINE CGEQRT2(M,N,A,Lda,T,Ldt,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0,0.0) , ZERO = (0.0,0.0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEQRT2
   END INTERFACE
END MODULE S_CGEQRT2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEQRT3
   INTERFACE
      RECURSIVE SUBROUTINE CGEQRT3(M,N,A,Lda,T,Ldt,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0,0.0)
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEQRT3
   END INTERFACE
END MODULE S_CGEQRT3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGEQRT
   INTERFACE
      SUBROUTINE CGEQRT(M,N,Nb,A,Lda,T,Ldt,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      LOGICAL , PARAMETER  ::  USE_RECURSIVE_QR = .TRUE.
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGEQRT
   END INTERFACE
END MODULE S_CGEQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGERFS
   INTERFACE
      SUBROUTINE CGERFS(Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,    &
     &                  Ferr,Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGERFS
   END INTERFACE
END MODULE S_CGERFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGERFSX
   INTERFACE
      SUBROUTINE CGERFSX(Trans,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,R,C,B,   &
     &                   Ldb,X,Ldx,Rcond,Berr,N_err_bnds,Err_bnds_norm, &
     &                   Err_bnds_comp,Nparams,Params,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ITREF_DEFAULT = 1.0 ,       &
     &                      ITHRESH_DEFAULT = 10.0 ,                    &
     &                      COMPONENTWISE_DEFAULT = 1.0 ,               &
     &                      RTHRESH_DEFAULT = 0.5 ,                     &
     &                      DZTHRESH_DEFAULT = 0.25
      INTEGER , PARAMETER  ::  LA_LINRX_ITREF_I = 1 ,                   &
     &                         LA_LINRX_ITHRESH_I = 2 ,                 &
     &                         LA_LINRX_CWISE_I = 3 ,                   &
     &                         LA_LINRX_TRUST_I = 1 ,                   &
     &                         LA_LINRX_ERR_I = 2 , LA_LINRX_RCOND_I = 3
      CHARACTER :: Trans
      CHARACTER :: Equed
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: R
      REAL , DIMENSION(*) :: C
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL , INTENT(INOUT) , DIMENSION(*) :: Params
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGERFSX
   END INTERFACE
END MODULE S_CGERFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGERQ2
   INTERFACE
      SUBROUTINE CGERQ2(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGERQ2
   END INTERFACE
END MODULE S_CGERQ2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGERQF
   INTERFACE
      SUBROUTINE CGERQF(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGERQF
   END INTERFACE
END MODULE S_CGERQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGESC2
   INTERFACE
      SUBROUTINE CGESC2(N,A,Lda,Rhs,Ipiv,Jpiv,Scale)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Rhs
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Jpiv
      REAL , INTENT(INOUT) :: Scale
      END SUBROUTINE CGESC2
   END INTERFACE
END MODULE S_CGESC2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGESDD
   INTERFACE
      SUBROUTINE CGESDD(Jobz,M,N,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,Rwork,&
     &                  Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobz
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: S
      COMPLEX , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGESDD
   END INTERFACE
END MODULE S_CGESDD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGESVD
   INTERFACE
      SUBROUTINE CGESVD(Jobu,Jobvt,M,N,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,&
     &                  Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobu
      CHARACTER :: Jobvt
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: S
      COMPLEX , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGESVD
   END INTERFACE
END MODULE S_CGESVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGESVDQ
   INTERFACE
      SUBROUTINE CGESVDQ(Joba,Jobp,Jobr,Jobu,Jobv,M,N,A,Lda,S,U,Ldu,V,  &
     &                   Ldv,Numrank,Iwork,Liwork,Cwork,Lcwork,Rwork,   &
     &                   Lrwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      CHARACTER :: Joba
      CHARACTER :: Jobp
      CHARACTER :: Jobr
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      INTEGER :: M
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: S
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(OUT) :: Numrank
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      COMPLEX , DIMENSION(*) :: Cwork
      INTEGER :: Lcwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGESVDQ
   END INTERFACE
END MODULE S_CGESVDQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGESVDX
   INTERFACE
      SUBROUTINE CGESVDX(Jobu,Jobvt,Range,M,N,A,Lda,Vl,Vu,Il,Iu,Ns,S,U, &
     &                   Ldu,Vt,Ldvt,Work,Lwork,Rwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0)
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobu
      CHARACTER :: Jobvt
      CHARACTER :: Range
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL :: Vl
      REAL :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      INTEGER , INTENT(INOUT) :: Ns
      REAL , DIMENSION(*) :: S
      COMPLEX , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGESVDX
   END INTERFACE
END MODULE S_CGESVDX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGESV
   INTERFACE
      SUBROUTINE CGESV(N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGESV
   END INTERFACE
END MODULE S_CGESV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGESVJ
   INTERFACE
      SUBROUTINE CGESVJ(Joba,Jobu,Jobv,M,N,A,Lda,Sva,Mv,V,Ldv,Cwork,    &
     &                  Lwork,Rwork,Lrwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , HALF = 0.5E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      INTEGER , PARAMETER  ::  NSWEEP = 30
      CHARACTER(1) :: Joba
      CHARACTER(1) :: Jobu
      CHARACTER(1) :: Jobv
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , INTENT(INOUT) , DIMENSION(Lwork) :: Cwork
      INTEGER , INTENT(IN) :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(Lrwork) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGESVJ
   END INTERFACE
END MODULE S_CGESVJ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGESVX
   INTERFACE
      SUBROUTINE CGESVX(Fact,Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,R,C, &
     &                  B,Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: R
      REAL , DIMENSION(*) :: C
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGESVX
   END INTERFACE
END MODULE S_CGESVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGESVXX
   INTERFACE
      SUBROUTINE CGESVXX(Fact,Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,R,C,&
     &                   B,Ldb,X,Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,      &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL , INTENT(INOUT) , DIMENSION(*) :: R
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , INTENT(OUT) :: Rpvgrw
      REAL , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL , DIMENSION(*) :: Params
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGESVXX
   END INTERFACE
END MODULE S_CGESVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGETC2
   INTERFACE
      SUBROUTINE CGETC2(N,A,Lda,Ipiv,Jpiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Jpiv
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE CGETC2
   END INTERFACE
END MODULE S_CGETC2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGETF2
   INTERFACE
      SUBROUTINE CGETF2(M,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGETF2
   END INTERFACE
END MODULE S_CGETF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGETRF2
   INTERFACE
      RECURSIVE SUBROUTINE CGETRF2(M,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGETRF2
   END INTERFACE
END MODULE S_CGETRF2
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
MODULE S_CGETRI
   INTERFACE
      SUBROUTINE CGETRI(N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGETRI
   END INTERFACE
END MODULE S_CGETRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGETRS
   INTERFACE
      SUBROUTINE CGETRS(Trans,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGETRS
   END INTERFACE
END MODULE S_CGETRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGETSLS
   INTERFACE
      SUBROUTINE CGETSLS(Trans,M,N,Nrhs,A,Lda,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGETSLS
   END INTERFACE
END MODULE S_CGETSLS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGETSQRHRT
   INTERFACE
      SUBROUTINE CGETSQRHRT(M,N,Mb1,Nb1,Nb2,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb1
      INTEGER , INTENT(IN) :: Nb1
      INTEGER , INTENT(IN) :: Nb2
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGETSQRHRT
   END INTERFACE
END MODULE S_CGETSQRHRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGBAK
   INTERFACE
      SUBROUTINE CGGBAK(Job,Side,N,Ilo,Ihi,Lscale,Rscale,M,V,Ldv,Info)
      IMPLICIT NONE
      CHARACTER :: Job
      CHARACTER :: Side
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL , DIMENSION(*) :: Lscale
      REAL , DIMENSION(*) :: Rscale
      INTEGER :: M
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGBAK
   END INTERFACE
END MODULE S_CGGBAK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGBAL
   INTERFACE
      SUBROUTINE CGGBAL(Job,N,A,Lda,B,Ldb,Ilo,Ihi,Lscale,Rscale,Work,   &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , HALF = 0.5E+0 ,             &
     &                      ONE = 1.0E+0 , THREE = 3.0E+0 ,             &
     &                      SCLFAC = 1.0E+1
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Job
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Ilo
      INTEGER , INTENT(INOUT) :: Ihi
      REAL , INTENT(INOUT) , DIMENSION(*) :: Lscale
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rscale
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGBAL
   END INTERFACE
END MODULE S_CGGBAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGES3
   INTERFACE
      SUBROUTINE CGGES3(Jobvsl,Jobvsr,Sort,SELCTG,N,A,Lda,B,Ldb,Sdim,   &
     &                  Alpha,Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work,Lwork,Rwork,&
     &                  Bwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      COMPLEX , DIMENSION(*) :: Alpha
      COMPLEX , DIMENSION(*) :: Beta
      COMPLEX , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      COMPLEX , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGES3
   END INTERFACE
END MODULE S_CGGES3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGES
   INTERFACE
      SUBROUTINE CGGES(Jobvsl,Jobvsr,Sort,SELCTG,N,A,Lda,B,Ldb,Sdim,    &
     &                 Alpha,Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work,Lwork,Rwork, &
     &                 Bwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      COMPLEX , DIMENSION(*) :: Alpha
      COMPLEX , DIMENSION(*) :: Beta
      COMPLEX , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      COMPLEX , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGES
   END INTERFACE
END MODULE S_CGGES
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGESX
   INTERFACE
      SUBROUTINE CGGESX(Jobvsl,Jobvsr,Sort,SELCTG,Sense,N,A,Lda,B,Ldb,  &
     &                  Sdim,Alpha,Beta,Vsl,Ldvsl,Vsr,Ldvsr,Rconde,     &
     &                  Rcondv,Work,Lwork,Rwork,Iwork,Liwork,Bwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      CHARACTER :: Sense
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      COMPLEX , DIMENSION(*) :: Alpha
      COMPLEX , DIMENSION(*) :: Beta
      COMPLEX , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      COMPLEX , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      REAL , INTENT(OUT) , DIMENSION(2) :: Rconde
      REAL , INTENT(OUT) , DIMENSION(2) :: Rcondv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGESX
   END INTERFACE
END MODULE S_CGGESX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGEV3
   INTERFACE
      SUBROUTINE CGGEV3(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alpha,Beta,Vl,Ldvl,Vr,&
     &                  Ldvr,Work,Lwork,Rwork,Info)
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
      COMPLEX , DIMENSION(*) :: Alpha
      COMPLEX , DIMENSION(*) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGEV3
   END INTERFACE
END MODULE S_CGGEV3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGEV
   INTERFACE
      SUBROUTINE CGGEV(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alpha,Beta,Vl,Ldvl,Vr, &
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
      COMPLEX , DIMENSION(*) :: Alpha
      COMPLEX , DIMENSION(*) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGEV
   END INTERFACE
END MODULE S_CGGEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGEVX
   INTERFACE
      SUBROUTINE CGGEVX(Balanc,Jobvl,Jobvr,Sense,N,A,Lda,B,Ldb,Alpha,   &
     &                  Beta,Vl,Ldvl,Vr,Ldvr,Ilo,Ihi,Lscale,Rscale,     &
     &                  Abnrm,Bbnrm,Rconde,Rcondv,Work,Lwork,Rwork,     &
     &                  Iwork,Bwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Balanc
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      CHARACTER :: Sense
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Alpha
      COMPLEX , DIMENSION(*) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL , DIMENSION(*) :: Lscale
      REAL , DIMENSION(*) :: Rscale
      REAL , INTENT(INOUT) :: Abnrm
      REAL , INTENT(INOUT) :: Bbnrm
      REAL , DIMENSION(*) :: Rconde
      REAL , DIMENSION(*) :: Rcondv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGEVX
   END INTERFACE
END MODULE S_CGGEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGGLM
   INTERFACE
      SUBROUTINE CGGGLM(N,M,P,A,Lda,B,Ldb,D,X,Y,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      INTEGER :: N
      INTEGER :: M
      INTEGER :: P
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: X
      COMPLEX , DIMENSION(*) :: Y
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGGLM
   END INTERFACE
END MODULE S_CGGGLM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGHD3
   INTERFACE
      SUBROUTINE CGGHD3(Compq,Compz,N,Ilo,Ihi,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,  &
     &                  Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGHD3
   END INTERFACE
END MODULE S_CGGHD3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGHRD
   INTERFACE
      SUBROUTINE CGGHRD(Compq,Compz,N,Ilo,Ihi,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,  &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER :: Ihi
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGHRD
   END INTERFACE
END MODULE S_CGGHRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGLSE
   INTERFACE
      SUBROUTINE CGGLSE(M,N,P,A,Lda,B,Ldb,C,D,X,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: P
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: C
      COMPLEX , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: X
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGLSE
   END INTERFACE
END MODULE S_CGGLSE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGQRF
   INTERFACE
      SUBROUTINE CGGQRF(N,M,P,A,Lda,Taua,B,Ldb,Taub,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: M
      INTEGER :: P
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Taua
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Taub
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGQRF
   END INTERFACE
END MODULE S_CGGQRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGRQF
   INTERFACE
      SUBROUTINE CGGRQF(M,P,N,A,Lda,Taua,B,Ldb,Taub,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: P
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Taua
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Taub
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGRQF
   END INTERFACE
END MODULE S_CGGRQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGSVD3
   INTERFACE
      SUBROUTINE CGGSVD3(Jobu,Jobv,Jobq,M,N,P,K,L,A,Lda,B,Ldb,Alpha,    &
     &                   Beta,U,Ldu,V,Ldv,Q,Ldq,Work,Lwork,Rwork,Iwork, &
     &                   Info)
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
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGSVD3
   END INTERFACE
END MODULE S_CGGSVD3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGGSVP3
   INTERFACE
      SUBROUTINE CGGSVP3(Jobu,Jobv,Jobq,M,P,N,A,Lda,B,Ldb,Tola,Tolb,K,L,&
     &                   U,Ldu,V,Ldv,Q,Ldq,Iwork,Rwork,Tau,Work,Lwork,  &
     &                   Info)
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
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGGSVP3
   END INTERFACE
END MODULE S_CGGSVP3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGSVJ0
   INTERFACE
      SUBROUTINE CGSVJ0(Jobv,M,N,A,Lda,D,Sva,Mv,V,Ldv,Eps,Sfmin,Tol,    &
     &                  Nsweep,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , HALF = 0.5E0 , ONE = 1.0E0
      CHARACTER(1) :: Jobv
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(N) :: D
      REAL , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      REAL , INTENT(IN) :: Eps
      REAL , INTENT(IN) :: Sfmin
      REAL , INTENT(IN) :: Tol
      INTEGER , INTENT(IN) :: Nsweep
      COMPLEX , DIMENSION(Lwork) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGSVJ0
   END INTERFACE
END MODULE S_CGSVJ0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGSVJ1
   INTERFACE
      SUBROUTINE CGSVJ1(Jobv,M,N,N1,A,Lda,D,Sva,Mv,V,Ldv,Eps,Sfmin,Tol, &
     &                  Nsweep,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , HALF = 0.5E0 , ONE = 1.0E0
      CHARACTER(1) :: Jobv
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: N1
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(N) :: D
      REAL , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      REAL , INTENT(IN) :: Eps
      REAL , INTENT(IN) :: Sfmin
      REAL , INTENT(IN) :: Tol
      INTEGER , INTENT(IN) :: Nsweep
      COMPLEX , DIMENSION(Lwork) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGSVJ1
   END INTERFACE
END MODULE S_CGSVJ1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGTCON
   INTERFACE
      SUBROUTINE CGTCON(Norm,N,Dl,D,Du,Du2,Ipiv,Anorm,Rcond,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Norm
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Dl
      COMPLEX , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: Du
      COMPLEX , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGTCON
   END INTERFACE
END MODULE S_CGTCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGTRFS
   INTERFACE
      SUBROUTINE CGTRFS(Trans,N,Nrhs,Dl,D,Du,Dlf,Df,Duf,Du2,Ipiv,B,Ldb, &
     &                  X,Ldx,Ferr,Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(*) :: Dl
      COMPLEX , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: Du
      COMPLEX , DIMENSION(*) :: Dlf
      COMPLEX , DIMENSION(*) :: Df
      COMPLEX , DIMENSION(*) :: Duf
      COMPLEX , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGTRFS
   END INTERFACE
END MODULE S_CGTRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGTSV
   INTERFACE
      SUBROUTINE CGTSV(N,Nrhs,Dl,D,Du,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Dl
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: D
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Du
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGTSV
   END INTERFACE
END MODULE S_CGTSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGTSVX
   INTERFACE
      SUBROUTINE CGTSVX(Fact,Trans,N,Nrhs,Dl,D,Du,Dlf,Df,Duf,Du2,Ipiv,B,&
     &                  Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(*) :: Dl
      COMPLEX , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: Du
      COMPLEX , DIMENSION(*) :: Dlf
      COMPLEX , DIMENSION(*) :: Df
      COMPLEX , DIMENSION(*) :: Duf
      COMPLEX , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGTSVX
   END INTERFACE
END MODULE S_CGTSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGTTRF
   INTERFACE
      SUBROUTINE CGTTRF(N,Dl,D,Du,Du2,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Dl
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: D
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Du
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: Du2
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER :: Info
      END SUBROUTINE CGTTRF
   END INTERFACE
END MODULE S_CGTTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGTTRS
   INTERFACE
      SUBROUTINE CGTTRS(Trans,N,Nrhs,Dl,D,Du,Du2,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(*) :: Dl
      COMPLEX , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: Du
      COMPLEX , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CGTTRS
   END INTERFACE
END MODULE S_CGTTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CGTTS2
   INTERFACE
      SUBROUTINE CGTTS2(Itrans,N,Nrhs,Dl,D,Du,Du2,Ipiv,B,Ldb)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: Itrans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Dl
      COMPLEX , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Du
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Du2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE CGTTS2
   END INTERFACE
END MODULE S_CGTTS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHB2ST_KERNELS
   INTERFACE
      SUBROUTINE CHB2ST_KERNELS(Uplo,Wantz,Ttype,St,Ed,Sweep,N,Nb,Ib,A, &
     &                          Lda,V,Tau,Ldvt,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: Ttype
      INTEGER , INTENT(IN) :: St
      INTEGER , INTENT(IN) :: Ed
      INTEGER , INTENT(IN) :: Sweep
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(IN) :: Ib
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , DIMENSION(*) :: V
      COMPLEX , DIMENSION(*) :: Tau
      INTEGER , INTENT(IN) :: Ldvt
      COMPLEX , DIMENSION(*) :: Work
      END SUBROUTINE CHB2ST_KERNELS
   END INTERFACE
END MODULE S_CHB2ST_KERNELS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHBEV_2STAGE
   INTERFACE
      SUBROUTINE CHBEV_2STAGE(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Lwork,&
     &                        Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHBEV_2STAGE
   END INTERFACE
END MODULE S_CHBEV_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHBEVD_2STAGE
   INTERFACE
      SUBROUTINE CHBEVD_2STAGE(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,     &
     &                         Lwork,Rwork,Lrwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHBEVD_2STAGE
   END INTERFACE
END MODULE S_CHBEVD_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHBEVD
   INTERFACE
      SUBROUTINE CHBEVD(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Lwork,Rwork,&
     &                  Lrwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHBEVD
   END INTERFACE
END MODULE S_CHBEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHBEV
   INTERFACE
      SUBROUTINE CHBEV(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHBEV
   END INTERFACE
END MODULE S_CHBEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHBEVX_2STAGE
   INTERFACE
      SUBROUTINE CHBEVX_2STAGE(Jobz,Range,Uplo,N,Kd,Ab,Ldab,Q,Ldq,Vl,Vu,&
     &                         Il,Iu,Abstol,M,W,Z,Ldz,Work,Lwork,Rwork, &
     &                         Iwork,Ifail,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHBEVX_2STAGE
   END INTERFACE
END MODULE S_CHBEVX_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHBEVX
   INTERFACE
      SUBROUTINE CHBEVX(Jobz,Range,Uplo,N,Kd,Ab,Ldab,Q,Ldq,Vl,Vu,Il,Iu, &
     &                  Abstol,M,W,Z,Ldz,Work,Rwork,Iwork,Ifail,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHBEVX
   END INTERFACE
END MODULE S_CHBEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHBGST
   INTERFACE
      SUBROUTINE CHBGST(Vect,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,X,Ldx,Work,   &
     &                  Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ka
      INTEGER , INTENT(IN) :: Kb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX , DIMENSION(Ldbb,*) :: Bb
      INTEGER , INTENT(IN) :: Ldbb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHBGST
   END INTERFACE
END MODULE S_CHBGST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHBGVD
   INTERFACE
      SUBROUTINE CHBGVD(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,W,Z,Ldz,Work, &
     &                  Lwork,Rwork,Lrwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ka
      INTEGER :: Kb
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldbb,*) :: Bb
      INTEGER :: Ldbb
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHBGVD
   END INTERFACE
END MODULE S_CHBGVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHBGV
   INTERFACE
      SUBROUTINE CHBGV(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,W,Z,Ldz,Work,  &
     &                 Rwork,Info)
      IMPLICIT NONE
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ka
      INTEGER :: Kb
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldbb,*) :: Bb
      INTEGER :: Ldbb
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHBGV
   END INTERFACE
END MODULE S_CHBGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHBGVX
   INTERFACE
      SUBROUTINE CHBGVX(Jobz,Range,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,Q,Ldq,  &
     &                  Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,Work,Rwork,Iwork,  &
     &                  Ifail,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ka
      INTEGER :: Kb
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldbb,*) :: Bb
      INTEGER :: Ldbb
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHBGVX
   END INTERFACE
END MODULE S_CHBGVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHBTRD
   INTERFACE
      SUBROUTINE CHBTRD(Vect,Uplo,N,Kd,Ab,Ldab,D,E,Q,Ldq,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHBTRD
   END INTERFACE
END MODULE S_CHBTRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHECON_3
   INTERFACE
      SUBROUTINE CHECON_3(Uplo,N,A,Lda,E,Ipiv,Anorm,Rcond,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHECON_3
   END INTERFACE
END MODULE S_CHECON_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHECON
   INTERFACE
      SUBROUTINE CHECON(Uplo,N,A,Lda,Ipiv,Anorm,Rcond,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHECON
   END INTERFACE
END MODULE S_CHECON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHECON_ROOK
   INTERFACE
      SUBROUTINE CHECON_ROOK(Uplo,N,A,Lda,Ipiv,Anorm,Rcond,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHECON_ROOK
   END INTERFACE
END MODULE S_CHECON_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEEQUB
   INTERFACE
      SUBROUTINE CHEEQUB(Uplo,N,A,Lda,S,Scond,Amax,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , ZERO = 0.0E0
      INTEGER , PARAMETER  ::  MAX_ITER = 100
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) :: Scond
      REAL , INTENT(INOUT) :: Amax
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEEQUB
   END INTERFACE
END MODULE S_CHEEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEEV_2STAGE
   INTERFACE
      SUBROUTINE CHEEV_2STAGE(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEEV_2STAGE
   END INTERFACE
END MODULE S_CHEEV_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEEVD_2STAGE
   INTERFACE
      SUBROUTINE CHEEVD_2STAGE(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Rwork,    &
     &                         Lrwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEEVD_2STAGE
   END INTERFACE
END MODULE S_CHEEVD_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEEVD
   INTERFACE
      SUBROUTINE CHEEVD(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Rwork,Lrwork,    &
     &                  Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEEVD
   END INTERFACE
END MODULE S_CHEEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEEV
   INTERFACE
      SUBROUTINE CHEEV(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEEV
   END INTERFACE
END MODULE S_CHEEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEEVR_2STAGE
   INTERFACE
      SUBROUTINE CHEEVR_2STAGE(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,     &
     &                         Abstol,M,W,Z,Ldz,Isuppz,Work,Lwork,Rwork,&
     &                         Lrwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEEVR_2STAGE
   END INTERFACE
END MODULE S_CHEEVR_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEEVR
   INTERFACE
      SUBROUTINE CHEEVR(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,Abstol,M,W, &
     &                  Z,Ldz,Isuppz,Work,Lwork,Rwork,Lrwork,Iwork,     &
     &                  Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEEVR
   END INTERFACE
END MODULE S_CHEEVR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEEVX_2STAGE
   INTERFACE
      SUBROUTINE CHEEVX_2STAGE(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,     &
     &                         Abstol,M,W,Z,Ldz,Work,Lwork,Rwork,Iwork, &
     &                         Ifail,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEEVX_2STAGE
   END INTERFACE
END MODULE S_CHEEVX_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEEVX
   INTERFACE
      SUBROUTINE CHEEVX(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,Abstol,M,W, &
     &                  Z,Ldz,Work,Lwork,Rwork,Iwork,Ifail,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEEVX
   END INTERFACE
END MODULE S_CHEEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEGS2
   INTERFACE
      SUBROUTINE CHEGS2(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , HALF = 0.5E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: Itype
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEGS2
   END INTERFACE
END MODULE S_CHEGS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEGST
   INTERFACE
      SUBROUTINE CHEGST(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         HALF = (0.5E+0,0.0E+0)
      INTEGER :: Itype
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEGST
   END INTERFACE
END MODULE S_CHEGST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEGV_2STAGE
   INTERFACE
      SUBROUTINE CHEGV_2STAGE(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,     &
     &                        Lwork,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEGV_2STAGE
   END INTERFACE
END MODULE S_CHEGV_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEGVD
   INTERFACE
      SUBROUTINE CHEGVD(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,Lwork,     &
     &                  Rwork,Lrwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: W
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER :: Lrwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEGVD
   END INTERFACE
END MODULE S_CHEGVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEGV
   INTERFACE
      SUBROUTINE CHEGV(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,Lwork,Rwork,&
     &                 Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEGV
   END INTERFACE
END MODULE S_CHEGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHEGVX
   INTERFACE
      SUBROUTINE CHEGVX(Itype,Jobz,Range,Uplo,N,A,Lda,B,Ldb,Vl,Vu,Il,Iu,&
     &                  Abstol,M,W,Z,Ldz,Work,Lwork,Rwork,Iwork,Ifail,  &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL :: Abstol
      INTEGER :: M
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHEGVX
   END INTERFACE
END MODULE S_CHEGVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHERFS
   INTERFACE
      SUBROUTINE CHERFS(Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,Ferr,&
     &                  Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHERFS
   END INTERFACE
END MODULE S_CHERFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHERFSX
   INTERFACE
      SUBROUTINE CHERFSX(Uplo,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,S,B,Ldb,X,&
     &                   Ldx,Rcond,Berr,N_err_bnds,Err_bnds_norm,       &
     &                   Err_bnds_comp,Nparams,Params,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ITREF_DEFAULT = 1.0 ,       &
     &                      ITHRESH_DEFAULT = 10.0 ,                    &
     &                      COMPONENTWISE_DEFAULT = 1.0 ,               &
     &                      RTHRESH_DEFAULT = 0.5 ,                     &
     &                      DZTHRESH_DEFAULT = 0.25
      INTEGER , PARAMETER  ::  LA_LINRX_ITREF_I = 1 ,                   &
     &                         LA_LINRX_ITHRESH_I = 2 ,                 &
     &                         LA_LINRX_CWISE_I = 3 ,                   &
     &                         LA_LINRX_TRUST_I = 1 ,                   &
     &                         LA_LINRX_ERR_I = 2 , LA_LINRX_RCOND_I = 3
      CHARACTER :: Uplo
      CHARACTER :: Equed
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: S
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL , INTENT(INOUT) , DIMENSION(*) :: Params
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHERFSX
   END INTERFACE
END MODULE S_CHERFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHESV_AA_2STAGE
   INTERFACE
      SUBROUTINE CHESV_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B, &
     &                           Ldb,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tb
      INTEGER :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHESV_AA_2STAGE
   END INTERFACE
END MODULE S_CHESV_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHESV_AA
   INTERFACE
      SUBROUTINE CHESV_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHESV_AA
   END INTERFACE
END MODULE S_CHESV_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHESV
   INTERFACE
      SUBROUTINE CHESV(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHESV
   END INTERFACE
END MODULE S_CHESV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHESV_RK
   INTERFACE
      SUBROUTINE CHESV_RK(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Work,Lwork,    &
     &                    Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHESV_RK
   END INTERFACE
END MODULE S_CHESV_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHESV_ROOK
   INTERFACE
      SUBROUTINE CHESV_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,    &
     &                      Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHESV_ROOK
   END INTERFACE
END MODULE S_CHESV_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHESVX
   INTERFACE
      SUBROUTINE CHESVX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,&
     &                  Rcond,Ferr,Berr,Work,Lwork,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHESVX
   END INTERFACE
END MODULE S_CHESVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHESVXX
   INTERFACE
      SUBROUTINE CHESVXX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,S,B, &
     &                   Ldb,X,Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,        &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: S
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , INTENT(OUT) :: Rpvgrw
      REAL , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL , DIMENSION(*) :: Params
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHESVXX
   END INTERFACE
END MODULE S_CHESVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHESWAPR
   INTERFACE
      SUBROUTINE CHESWAPR(Uplo,N,A,Lda,I1,I2)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,N) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) :: I1
      INTEGER , INTENT(IN) :: I2
      END SUBROUTINE CHESWAPR
   END INTERFACE
END MODULE S_CHESWAPR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETD2
   INTERFACE
      SUBROUTINE CHETD2(Uplo,N,A,Lda,D,E,Tau,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         HALF = (0.5E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(OUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETD2
   END INTERFACE
END MODULE S_CHETD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETF2
   INTERFACE
      SUBROUTINE CHETF2(Uplo,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETF2
   END INTERFACE
END MODULE S_CHETF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETF2_RK
   INTERFACE
      SUBROUTINE CHETF2_RK(Uplo,N,A,Lda,E,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETF2_RK
   END INTERFACE
END MODULE S_CHETF2_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETF2_ROOK
   INTERFACE
      SUBROUTINE CHETF2_ROOK(Uplo,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETF2_ROOK
   END INTERFACE
END MODULE S_CHETF2_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRD_2STAGE
   INTERFACE
      SUBROUTINE CHETRD_2STAGE(Vect,Uplo,N,A,Lda,D,E,Tau,Hous2,Lhous2,  &
     &                         Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Hous2
      INTEGER :: Lhous2
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRD_2STAGE
   END INTERFACE
END MODULE S_CHETRD_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRD
   INTERFACE
      SUBROUTINE CHETRD(Uplo,N,A,Lda,D,E,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRD
   END INTERFACE
END MODULE S_CHETRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRD_HE2HB
   INTERFACE
      SUBROUTINE CHETRD_HE2HB(Uplo,N,Kd,A,Lda,Ab,Ldab,Tau,Work,Lwork,   &
     &                        Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  RONE = 1.0E+0
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0) ,                  &
     &                         HALF = (0.5E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRD_HE2HB
   END INTERFACE
END MODULE S_CHETRD_HE2HB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRF_AA_2STAGE
   INTERFACE
      SUBROUTINE CHETRF_AA_2STAGE(Uplo,N,A,Lda,Tb,Ltb,Ipiv,Ipiv2,Work,  &
     &                            Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRF_AA_2STAGE
   END INTERFACE
END MODULE S_CHETRF_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRF_AA
   INTERFACE
      SUBROUTINE CHETRF_AA(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRF_AA
   END INTERFACE
END MODULE S_CHETRF_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRF
   INTERFACE
      SUBROUTINE CHETRF(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRF
   END INTERFACE
END MODULE S_CHETRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRF_RK
   INTERFACE
      SUBROUTINE CHETRF_RK(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRF_RK
   END INTERFACE
END MODULE S_CHETRF_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRF_ROOK
   INTERFACE
      SUBROUTINE CHETRF_ROOK(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRF_ROOK
   END INTERFACE
END MODULE S_CHETRF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRI2
   INTERFACE
      SUBROUTINE CHETRI2(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRI2
   END INTERFACE
END MODULE S_CHETRI2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRI2X
   INTERFACE
      SUBROUTINE CHETRI2X(Uplo,N,A,Lda,Ipiv,Work,Nb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRI2X
   END INTERFACE
END MODULE S_CHETRI2X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRI_3
   INTERFACE
      SUBROUTINE CHETRI_3(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRI_3
   END INTERFACE
END MODULE S_CHETRI_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRI_3X
   INTERFACE
      SUBROUTINE CHETRI_3X(Uplo,N,A,Lda,E,Ipiv,Work,Nb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRI_3X
   END INTERFACE
END MODULE S_CHETRI_3X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRI
   INTERFACE
      SUBROUTINE CHETRI(Uplo,N,A,Lda,Ipiv,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRI
   END INTERFACE
END MODULE S_CHETRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRI_ROOK
   INTERFACE
      SUBROUTINE CHETRI_ROOK(Uplo,N,A,Lda,Ipiv,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRI_ROOK
   END INTERFACE
END MODULE S_CHETRI_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRS2
   INTERFACE
      SUBROUTINE CHETRS2(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRS2
   END INTERFACE
END MODULE S_CHETRS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRS_3
   INTERFACE
      SUBROUTINE CHETRS_3(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRS_3
   END INTERFACE
END MODULE S_CHETRS_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRS_AA_2STAGE
   INTERFACE
      SUBROUTINE CHETRS_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B,&
     &                            Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRS_AA_2STAGE
   END INTERFACE
END MODULE S_CHETRS_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRS_AA
   INTERFACE
      SUBROUTINE CHETRS_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRS_AA
   END INTERFACE
END MODULE S_CHETRS_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRS
   INTERFACE
      SUBROUTINE CHETRS(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRS
   END INTERFACE
END MODULE S_CHETRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHETRS_ROOK
   INTERFACE
      SUBROUTINE CHETRS_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHETRS_ROOK
   END INTERFACE
END MODULE S_CHETRS_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHFRK
   INTERFACE
      SUBROUTINE CHFRK(Transr,Uplo,Trans,N,K,Alpha,A,Lda,Beta,C)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Transr
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: K
      REAL :: Alpha
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL :: Beta
      COMPLEX , DIMENSION(*) :: C
      END SUBROUTINE CHFRK
   END INTERFACE
END MODULE S_CHFRK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHGEQZ
   INTERFACE
      SUBROUTINE CHGEQZ(Job,Compq,Compz,N,Ilo,Ihi,H,Ldh,T,Ldt,Alpha,    &
     &                  Beta,Q,Ldq,Z,Ldz,Work,Lwork,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , HALF = 0.5E+0
      CHARACTER :: Job
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: Alpha
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHGEQZ
   END INTERFACE
END MODULE S_CHGEQZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHLA_TRANSTYPE
   INTERFACE
      FUNCTION CHLA_TRANSTYPE(Trans)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  BLAS_NO_TRANS = 111 , BLAS_TRANS = 112 , &
     &                         BLAS_CONJ_TRANS = 113
      CHARACTER(1) :: CHLA_TRANSTYPE
      INTEGER , INTENT(IN) :: Trans
      END FUNCTION CHLA_TRANSTYPE
   END INTERFACE
END MODULE S_CHLA_TRANSTYPE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPCON
   INTERFACE
      SUBROUTINE CHPCON(Uplo,N,Ap,Ipiv,Anorm,Rcond,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPCON
   END INTERFACE
END MODULE S_CHPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPEVD
   INTERFACE
      SUBROUTINE CHPEVD(Jobz,Uplo,N,Ap,W,Z,Ldz,Work,Lwork,Rwork,Lrwork, &
     &                  Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPEVD
   END INTERFACE
END MODULE S_CHPEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPEV
   INTERFACE
      SUBROUTINE CHPEV(Jobz,Uplo,N,Ap,W,Z,Ldz,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPEV
   END INTERFACE
END MODULE S_CHPEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPEVX
   INTERFACE
      SUBROUTINE CHPEVX(Jobz,Range,Uplo,N,Ap,Vl,Vu,Il,Iu,Abstol,M,W,Z,  &
     &                  Ldz,Work,Rwork,Iwork,Ifail,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CONE = (1.0E0,0.0E0)
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPEVX
   END INTERFACE
END MODULE S_CHPEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPGST
   INTERFACE
      SUBROUTINE CHPGST(Itype,Uplo,N,Ap,Bp,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , HALF = 0.5E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: Itype
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Bp
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPGST
   END INTERFACE
END MODULE S_CHPGST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPGVD
   INTERFACE
      SUBROUTINE CHPGVD(Itype,Jobz,Uplo,N,Ap,Bp,W,Z,Ldz,Work,Lwork,     &
     &                  Rwork,Lrwork,Iwork,Liwork,Info)
      IMPLICIT NONE
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Bp
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER :: Lrwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPGVD
   END INTERFACE
END MODULE S_CHPGVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPGV
   INTERFACE
      SUBROUTINE CHPGV(Itype,Jobz,Uplo,N,Ap,Bp,W,Z,Ldz,Work,Rwork,Info)
      IMPLICIT NONE
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Bp
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPGV
   END INTERFACE
END MODULE S_CHPGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPGVX
   INTERFACE
      SUBROUTINE CHPGVX(Itype,Jobz,Range,Uplo,N,Ap,Bp,Vl,Vu,Il,Iu,      &
     &                  Abstol,M,W,Z,Ldz,Work,Rwork,Iwork,Ifail,Info)
      IMPLICIT NONE
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Bp
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPGVX
   END INTERFACE
END MODULE S_CHPGVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPRFS
   INTERFACE
      SUBROUTINE CHPRFS(Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Ferr,Berr,  &
     &                  Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPRFS
   END INTERFACE
END MODULE S_CHPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPSV
   INTERFACE
      SUBROUTINE CHPSV(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPSV
   END INTERFACE
END MODULE S_CHPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPSVX
   INTERFACE
      SUBROUTINE CHPSVX(Fact,Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Rcond, &
     &                  Ferr,Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPSVX
   END INTERFACE
END MODULE S_CHPSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPTRD
   INTERFACE
      SUBROUTINE CHPTRD(Uplo,N,Ap,D,E,Tau,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         HALF = (0.5E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL , INTENT(OUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPTRD
   END INTERFACE
END MODULE S_CHPTRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPTRF
   INTERFACE
      SUBROUTINE CHPTRF(Uplo,N,Ap,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPTRF
   END INTERFACE
END MODULE S_CHPTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPTRI
   INTERFACE
      SUBROUTINE CHPTRI(Uplo,N,Ap,Ipiv,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPTRI
   END INTERFACE
END MODULE S_CHPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHPTRS
   INTERFACE
      SUBROUTINE CHPTRS(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHPTRS
   END INTERFACE
END MODULE S_CHPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHSEIN
   INTERFACE
      SUBROUTINE CHSEIN(Side,Eigsrc,Initv,Select,N,H,Ldh,W,Vl,Ldvl,Vr,  &
     &                  Ldvr,Mm,M,Work,Rwork,Ifaill,Ifailr,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      REAL , PARAMETER  ::  RZERO = 0.0E+0
      CHARACTER :: Side
      CHARACTER :: Eigsrc
      CHARACTER :: Initv
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      COMPLEX , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifaill
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifailr
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHSEIN
   END INTERFACE
END MODULE S_CHSEIN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CHSEQR
   INTERFACE
      SUBROUTINE CHSEQR(Job,Compz,N,Ilo,Ihi,H,Ldh,W,Z,Ldz,Work,Lwork,   &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NTINY = 15 , NL = 49
      COMPLEX , PARAMETER  ::  ZERO = (0.0E0,0.0E0) ,                   &
     &                         ONE = (1.0E0,0.0E0)
      REAL , PARAMETER  ::  RZERO = 0.0E0
      CHARACTER :: Job
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CHSEQR
   END INTERFACE
END MODULE S_CHSEQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLABRD
   INTERFACE
      SUBROUTINE CLABRD(M,N,Nb,A,Lda,D,E,Tauq,Taup,X,Ldx,Y,Ldy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(OUT) , DIMENSION(*) :: D
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      COMPLEX , DIMENSION(*) :: Tauq
      COMPLEX , DIMENSION(*) :: Taup
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      COMPLEX , DIMENSION(Ldy,*) :: Y
      INTEGER :: Ldy
      END SUBROUTINE CLABRD
   END INTERFACE
END MODULE S_CLABRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLACGV
   INTERFACE
      SUBROUTINE CLACGV(N,X,Incx)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE CLACGV
   END INTERFACE
END MODULE S_CLACGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLACN2
   INTERFACE
      SUBROUTINE CLACN2(N,V,X,Est,Kase,Isave)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ONE = 1.0E0 , TWO = 2.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: V
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      REAL , INTENT(INOUT) :: Est
      INTEGER , INTENT(INOUT) :: Kase
      INTEGER , INTENT(INOUT) , DIMENSION(3) :: Isave
      END SUBROUTINE CLACN2
   END INTERFACE
END MODULE S_CLACN2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLACON
   INTERFACE
      SUBROUTINE CLACON(N,V,X,Est,Kase)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ONE = 1.0E0 , TWO = 2.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(N) :: V
      COMPLEX , INTENT(INOUT) , DIMENSION(N) :: X
      REAL , INTENT(INOUT) :: Est
      INTEGER , INTENT(INOUT) :: Kase
      END SUBROUTINE CLACON
   END INTERFACE
END MODULE S_CLACON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLACP2
   INTERFACE
      SUBROUTINE CLACP2(Uplo,M,N,A,Lda,B,Ldb)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(OUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE CLACP2
   END INTERFACE
END MODULE S_CLACP2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLACPY
   INTERFACE
      SUBROUTINE CLACPY(Uplo,M,N,A,Lda,B,Ldb)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(OUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE CLACPY
   END INTERFACE
END MODULE S_CLACPY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLACRM
   INTERFACE
      SUBROUTINE CLACRM(M,N,A,Lda,B,Ldb,C,Ldc,Rwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , ZERO = 0.0E0
      INTEGER :: M
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END SUBROUTINE CLACRM
   END INTERFACE
END MODULE S_CLACRM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLACRT
   INTERFACE
      SUBROUTINE CLACRT(N,Cx,Incx,Cy,Incy,C,S)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
      COMPLEX , INTENT(IN) :: C
      COMPLEX , INTENT(IN) :: S
      END SUBROUTINE CLACRT
   END INTERFACE
END MODULE S_CLACRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLADIV
   INTERFACE
      FUNCTION CLADIV(X,Y)
      IMPLICIT NONE
      COMPLEX :: CLADIV
      COMPLEX , INTENT(IN) :: X
      COMPLEX , INTENT(IN) :: Y
      END FUNCTION CLADIV
   END INTERFACE
END MODULE S_CLADIV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAED0
   INTERFACE
      SUBROUTINE CLAED0(Qsiz,N,D,E,Q,Ldq,Qstore,Ldqs,Rwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  TWO = 2.E+0
      INTEGER :: Qsiz
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(Ldqs,*) :: Qstore
      INTEGER :: Ldqs
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLAED0
   END INTERFACE
END MODULE S_CLAED0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAED7
   INTERFACE
      SUBROUTINE CLAED7(N,Cutpnt,Qsiz,Tlvls,Curlvl,Curpbm,D,Q,Ldq,Rho,  &
     &                  Indxq,Qstore,Qptr,Prmptr,Perm,Givptr,Givcol,    &
     &                  Givnum,Work,Rwork,Iwork,Info)
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Cutpnt
      INTEGER :: Qsiz
      INTEGER :: Tlvls
      INTEGER :: Curlvl
      INTEGER :: Curpbm
      REAL , DIMENSION(*) :: D
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL :: Rho
      INTEGER , DIMENSION(*) :: Indxq
      REAL , DIMENSION(*) :: Qstore
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Qptr
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Prmptr
      INTEGER , DIMENSION(*) :: Perm
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Givptr
      INTEGER , DIMENSION(2,*) :: Givcol
      REAL , DIMENSION(2,*) :: Givnum
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLAED7
   END INTERFACE
END MODULE S_CLAED7
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAED8
   INTERFACE
      SUBROUTINE CLAED8(K,N,Qsiz,Q,Ldq,D,Rho,Cutpnt,Z,Dlamda,Q2,Ldq2,W, &
     &                  Indxp,Indx,Indxq,Perm,Givptr,Givcol,Givnum,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  MONE = -1.0E0 , ZERO = 0.0E0 , ONE = 1.0E0 ,&
     &                      TWO = 2.0E0 , EIGHT = 8.0E0
      INTEGER , INTENT(INOUT) :: K
      INTEGER :: N
      INTEGER :: Qsiz
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) :: Rho
      INTEGER , INTENT(IN) :: Cutpnt
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dlamda
      COMPLEX , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxp
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxq
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Perm
      INTEGER , INTENT(INOUT) :: Givptr
      INTEGER , INTENT(OUT) , DIMENSION(2,*) :: Givcol
      REAL , INTENT(OUT) , DIMENSION(2,*) :: Givnum
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLAED8
   END INTERFACE
END MODULE S_CLAED8
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAEIN
   INTERFACE
      SUBROUTINE CLAEIN(Rightv,Noinit,N,H,Ldh,W,V,B,Ldb,Rwork,Eps3,     &
     &                  Smlnum,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , TENTH = 1.0E-1
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      LOGICAL , INTENT(IN) :: Rightv
      LOGICAL , INTENT(IN) :: Noinit
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Ldh,*) :: H
      INTEGER , INTENT(IN) :: Ldh
      COMPLEX , INTENT(IN) :: W
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: V
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Rwork
      REAL , INTENT(IN) :: Eps3
      REAL , INTENT(IN) :: Smlnum
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE CLAEIN
   END INTERFACE
END MODULE S_CLAEIN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAESY
   INTERFACE
      SUBROUTINE CLAESY(A,B,C,Rt1,Rt2,Evscal,Cs1,Sn1)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CONE = (1.0E0,0.0E0)
      REAL , PARAMETER  ::  HALF = 0.5E0 , THRESH = 0.1E0
      COMPLEX , INTENT(IN) :: A
      COMPLEX , INTENT(IN) :: B
      COMPLEX , INTENT(IN) :: C
      COMPLEX , INTENT(INOUT) :: Rt1
      COMPLEX , INTENT(INOUT) :: Rt2
      COMPLEX , INTENT(INOUT) :: Evscal
      COMPLEX , INTENT(OUT) :: Cs1
      COMPLEX , INTENT(INOUT) :: Sn1
      END SUBROUTINE CLAESY
   END INTERFACE
END MODULE S_CLAESY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAEV2
   INTERFACE
      SUBROUTINE CLAEV2(A,B,C,Rt1,Rt2,Cs1,Sn1)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , INTENT(IN) :: A
      COMPLEX , INTENT(IN) :: B
      COMPLEX , INTENT(IN) :: C
      REAL :: Rt1
      REAL :: Rt2
      REAL :: Cs1
      COMPLEX , INTENT(OUT) :: Sn1
      END SUBROUTINE CLAEV2
   END INTERFACE
END MODULE S_CLAEV2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAG2Z
   INTERFACE
      SUBROUTINE CLAG2Z(M,N,Sa,Ldsa,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Ldsa,*) :: Sa
      INTEGER , INTENT(IN) :: Ldsa
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE CLAG2Z
   END INTERFACE
END MODULE S_CLAG2Z
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_GBAMV
   INTERFACE
      SUBROUTINE CLA_GBAMV(Trans,M,N,Kl,Ku,Alpha,Ab,Ldab,X,Incx,Beta,Y, &
     &                     Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE CLA_GBAMV
   END INTERFACE
END MODULE S_CLA_GBAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_GBRCOND_C
   INTERFACE
      FUNCTION CLA_GBRCOND_C(Trans,N,Kl,Ku,Ab,Ldab,Afb,Ldafb,Ipiv,C,    &
     &                       Capply,Info,Work,Rwork)
      IMPLICIT NONE
      REAL :: CLA_GBRCOND_C
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      COMPLEX , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) , DIMENSION(*) :: C
      LOGICAL , INTENT(IN) :: Capply
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION CLA_GBRCOND_C
   END INTERFACE
END MODULE S_CLA_GBRCOND_C
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_GBRCOND_X
   INTERFACE
      FUNCTION CLA_GBRCOND_X(Trans,N,Kl,Ku,Ab,Ldab,Afb,Ldafb,Ipiv,X,    &
     &                       Info,Work,Rwork)
      IMPLICIT NONE
      REAL :: CLA_GBRCOND_X
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      COMPLEX , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION CLA_GBRCOND_X
   END INTERFACE
END MODULE S_CLA_GBRCOND_X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_GBRFSX_EXTENDED
   INTERFACE
      SUBROUTINE CLA_GBRFSX_EXTENDED(Prec_type,Trans_type,N,Kl,Ku,Nrhs, &
     &                               Ab,Ldab,Afb,Ldafb,Ipiv,Colequ,C,B, &
     &                               Ldb,Y,Ldy,Berr_out,N_norms,        &
     &                               Err_bnds_norm,Err_bnds_comp,Res,   &
     &                               Ayb,Dy,Y_tail,Rcond,Ithresh,       &
     &                               Rthresh,Dz_ub,Ignore_cwise,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  UNSTABLE_STATE = 0 , WORKING_STATE = 1 , &
     &                         CONV_STATE = 2 , NOPROG_STATE = 3 ,      &
     &                         BASE_RESIDUAL = 0 , EXTRA_RESIDUAL = 1 , &
     &                         EXTRA_Y = 2 , LA_LINRX_ERR_I = 2
      INTEGER :: Prec_type
      INTEGER :: Trans_type
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      COMPLEX , DIMENSION(*) :: Res
      REAL , DIMENSION(*) :: Ayb
      COMPLEX , DIMENSION(*) :: Dy
      COMPLEX , DIMENSION(*) :: Y_tail
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL , INTENT(IN) :: Rthresh
      REAL , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
      END SUBROUTINE CLA_GBRFSX_EXTENDED
   END INTERFACE
END MODULE S_CLA_GBRFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_GBRPVGRW
   INTERFACE
      FUNCTION CLA_GBRPVGRW(N,Kl,Ku,Ncols,Ab,Ldab,Afb,Ldafb)
      IMPLICIT NONE
      REAL :: CLA_GBRPVGRW
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      INTEGER , INTENT(IN) :: Ncols
      COMPLEX , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX , INTENT(IN) , DIMENSION(Ldafb,*) :: Afb
      INTEGER , INTENT(IN) :: Ldafb
      END FUNCTION CLA_GBRPVGRW
   END INTERFACE
END MODULE S_CLA_GBRPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_GEAMV
   INTERFACE
      SUBROUTINE CLA_GEAMV(Trans,M,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE CLA_GEAMV
   END INTERFACE
END MODULE S_CLA_GEAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_GERCOND_C
   INTERFACE
      FUNCTION CLA_GERCOND_C(Trans,N,A,Lda,Af,Ldaf,Ipiv,C,Capply,Info,  &
     &                       Work,Rwork)
      IMPLICIT NONE
      REAL :: CLA_GERCOND_C
      CHARACTER :: Trans
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) , DIMENSION(*) :: C
      LOGICAL , INTENT(IN) :: Capply
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION CLA_GERCOND_C
   END INTERFACE
END MODULE S_CLA_GERCOND_C
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_GERCOND_X
   INTERFACE
      FUNCTION CLA_GERCOND_X(Trans,N,A,Lda,Af,Ldaf,Ipiv,X,Info,Work,    &
     &                       Rwork)
      IMPLICIT NONE
      REAL :: CLA_GERCOND_X
      CHARACTER :: Trans
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION CLA_GERCOND_X
   END INTERFACE
END MODULE S_CLA_GERCOND_X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_GERFSX_EXTENDED
   INTERFACE
      SUBROUTINE CLA_GERFSX_EXTENDED(Prec_type,Trans_type,N,Nrhs,A,Lda, &
     &                               Af,Ldaf,Ipiv,Colequ,C,B,Ldb,Y,Ldy, &
     &                               Berr_out,N_norms,Errs_n,Errs_c,Res,&
     &                               Ayb,Dy,Y_tail,Rcond,Ithresh,       &
     &                               Rthresh,Dz_ub,Ignore_cwise,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  UNSTABLE_STATE = 0 , WORKING_STATE = 1 , &
     &                         CONV_STATE = 2 , NOPROG_STATE = 3 ,      &
     &                         BASE_RESIDUAL = 0 , EXTRA_RESIDUAL = 1 , &
     &                         EXTRA_Y = 2 , LA_LINRX_ERR_I = 2
      INTEGER :: Prec_type
      INTEGER :: Trans_type
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Errs_n
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Errs_c
      COMPLEX , DIMENSION(*) :: Res
      REAL , DIMENSION(*) :: Ayb
      COMPLEX , DIMENSION(*) :: Dy
      COMPLEX , DIMENSION(*) :: Y_tail
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL , INTENT(IN) :: Rthresh
      REAL , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
      END SUBROUTINE CLA_GERFSX_EXTENDED
   END INTERFACE
END MODULE S_CLA_GERFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_GERPVGRW
   INTERFACE
      FUNCTION CLA_GERPVGRW(N,Ncols,A,Lda,Af,Ldaf)
      IMPLICIT NONE
      REAL :: CLA_GERPVGRW
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ncols
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      END FUNCTION CLA_GERPVGRW
   END INTERFACE
END MODULE S_CLA_GERPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAGS2
   INTERFACE
      SUBROUTINE CLAGS2(Upper,A1,A2,A3,B1,B2,B3,Csu,Snu,Csv,Snv,Csq,Snq)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      LOGICAL , INTENT(IN) :: Upper
      REAL , INTENT(IN) :: A1
      COMPLEX , INTENT(IN) :: A2
      REAL , INTENT(IN) :: A3
      REAL , INTENT(IN) :: B1
      COMPLEX , INTENT(IN) :: B2
      REAL , INTENT(IN) :: B3
      REAL , INTENT(OUT) :: Csu
      COMPLEX , INTENT(OUT) :: Snu
      REAL , INTENT(OUT) :: Csv
      COMPLEX , INTENT(OUT) :: Snv
      REAL :: Csq
      COMPLEX :: Snq
      END SUBROUTINE CLAGS2
   END INTERFACE
END MODULE S_CLAGS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAGTM
   INTERFACE
      SUBROUTINE CLAGTM(Trans,N,Nrhs,Alpha,Dl,D,Du,X,Ldx,Beta,B,Ldb)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Dl
      COMPLEX , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Du
      COMPLEX , INTENT(IN) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE CLAGTM
   END INTERFACE
END MODULE S_CLAGTM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_HEAMV
   INTERFACE
      SUBROUTINE CLA_HEAMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE CLA_HEAMV
   END INTERFACE
END MODULE S_CLA_HEAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAHEF_AA
   INTERFACE
      SUBROUTINE CLAHEF_AA(Uplo,J1,M,Nb,A,Lda,Ipiv,H,Ldh,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: J1
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      END SUBROUTINE CLAHEF_AA
   END INTERFACE
END MODULE S_CLAHEF_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAHEF
   INTERFACE
      SUBROUTINE CLAHEF(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLAHEF
   END INTERFACE
END MODULE S_CLAHEF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAHEF_RK
   INTERFACE
      SUBROUTINE CLAHEF_RK(Uplo,N,Nb,Kb,A,Lda,E,Ipiv,W,Ldw,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLAHEF_RK
   END INTERFACE
END MODULE S_CLAHEF_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAHEF_ROOK
   INTERFACE
      SUBROUTINE CLAHEF_ROOK(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLAHEF_ROOK
   END INTERFACE
END MODULE S_CLAHEF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_HERCOND_C
   INTERFACE
      FUNCTION CLA_HERCOND_C(Uplo,N,A,Lda,Af,Ldaf,Ipiv,C,Capply,Info,   &
     &                       Work,Rwork)
      IMPLICIT NONE
      REAL :: CLA_HERCOND_C
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) , DIMENSION(*) :: C
      LOGICAL , INTENT(IN) :: Capply
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION CLA_HERCOND_C
   END INTERFACE
END MODULE S_CLA_HERCOND_C
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_HERCOND_X
   INTERFACE
      FUNCTION CLA_HERCOND_X(Uplo,N,A,Lda,Af,Ldaf,Ipiv,X,Info,Work,     &
     &                       Rwork)
      IMPLICIT NONE
      REAL :: CLA_HERCOND_X
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION CLA_HERCOND_X
   END INTERFACE
END MODULE S_CLA_HERCOND_X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_HERFSX_EXTENDED
   INTERFACE
      SUBROUTINE CLA_HERFSX_EXTENDED(Prec_type,Uplo,N,Nrhs,A,Lda,Af,    &
     &                               Ldaf,Ipiv,Colequ,C,B,Ldb,Y,Ldy,    &
     &                               Berr_out,N_norms,Err_bnds_norm,    &
     &                               Err_bnds_comp,Res,Ayb,Dy,Y_tail,   &
     &                               Rcond,Ithresh,Rthresh,Dz_ub,       &
     &                               Ignore_cwise,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  UNSTABLE_STATE = 0 , WORKING_STATE = 1 , &
     &                         CONV_STATE = 2 , NOPROG_STATE = 3 ,      &
     &                         BASE_RESIDUAL = 0 , EXTRA_RESIDUAL = 1 , &
     &                         EXTRA_Y = 2 , LA_LINRX_ERR_I = 2
      INTEGER :: Prec_type
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      COMPLEX , DIMENSION(*) :: Res
      REAL , DIMENSION(*) :: Ayb
      COMPLEX , DIMENSION(*) :: Dy
      COMPLEX , DIMENSION(*) :: Y_tail
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL , INTENT(IN) :: Rthresh
      REAL , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLA_HERFSX_EXTENDED
   END INTERFACE
END MODULE S_CLA_HERFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_HERPVGRW
   INTERFACE
      FUNCTION CLA_HERPVGRW(Uplo,N,Info,A,Lda,Af,Ldaf,Ipiv,Work)
      IMPLICIT NONE
      REAL :: CLA_HERPVGRW
      CHARACTER(1) :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Info
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLA_HERPVGRW
   END INTERFACE
END MODULE S_CLA_HERPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAHQR
   INTERFACE
      SUBROUTINE CLAHQR(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,W,Iloz,Ihiz,Z,Ldz,  &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E0,0.0E0) ,                   &
     &                         ONE = (1.0E0,0.0E0)
      REAL , PARAMETER  ::  RZERO = 0.0E0 , RONE = 1.0E0 ,              &
     &                      HALF = 0.5E0 , DAT1 = 3.0E0/4.0E0
      INTEGER , PARAMETER  ::  KEXSH = 10
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: W
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE CLAHQR
   END INTERFACE
END MODULE S_CLAHQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAHR2
   INTERFACE
      SUBROUTINE CLAHR2(N,K,Nb,A,Lda,Tau,T,Ldt,Y,Ldy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: N
      INTEGER :: K
      INTEGER :: Nb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Nb) :: Tau
      COMPLEX , DIMENSION(Ldt,Nb) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(Ldy,Nb) :: Y
      INTEGER :: Ldy
      END SUBROUTINE CLAHR2
   END INTERFACE
END MODULE S_CLAHR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAIC1
   INTERFACE
      SUBROUTINE CLAIC1(Job,J,X,Sest,W,Gamma,Sestpr,S,C)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      HALF = 0.5E0 , FOUR = 4.0E0
      INTEGER , INTENT(IN) :: Job
      INTEGER :: J
      COMPLEX , DIMENSION(J) :: X
      REAL , INTENT(IN) :: Sest
      COMPLEX , DIMENSION(J) :: W
      COMPLEX , INTENT(IN) :: Gamma
      REAL , INTENT(OUT) :: Sestpr
      COMPLEX , INTENT(INOUT) :: S
      COMPLEX , INTENT(INOUT) :: C
      END SUBROUTINE CLAIC1
   END INTERFACE
END MODULE S_CLAIC1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_LIN_BERR
   INTERFACE
      SUBROUTINE CLA_LIN_BERR(N,Nz,Nrhs,Res,Ayb,Berr)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nz
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , INTENT(IN) , DIMENSION(N,Nrhs) :: Res
      REAL , INTENT(IN) , DIMENSION(N,Nrhs) :: Ayb
      REAL , INTENT(INOUT) , DIMENSION(Nrhs) :: Berr
      END SUBROUTINE CLA_LIN_BERR
   END INTERFACE
END MODULE S_CLA_LIN_BERR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLALS0
   INTERFACE
      SUBROUTINE CLALS0(Icompq,Nl,Nr,Sqre,Nrhs,B,Ldb,Bx,Ldbx,Perm,      &
     &                  Givptr,Givcol,Ldgcol,Givnum,Ldgnum,Poles,Difl,  &
     &                  Difr,Z,K,C,S,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , ZERO = 0.0E0 , NEGONE = -1.0E0
      INTEGER , INTENT(IN) :: Icompq
      INTEGER , INTENT(IN) :: Nl
      INTEGER , INTENT(IN) :: Nr
      INTEGER , INTENT(IN) :: Sqre
      INTEGER :: Nrhs
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldbx,*) :: Bx
      INTEGER :: Ldbx
      INTEGER , INTENT(IN) , DIMENSION(*) :: Perm
      INTEGER , INTENT(IN) :: Givptr
      INTEGER , INTENT(IN) , DIMENSION(Ldgcol,*) :: Givcol
      INTEGER , INTENT(IN) :: Ldgcol
      REAL , DIMENSION(Ldgnum,*) :: Givnum
      INTEGER , INTENT(IN) :: Ldgnum
      REAL , DIMENSION(Ldgnum,*) :: Poles
      REAL , INTENT(IN) , DIMENSION(*) :: Difl
      REAL , INTENT(IN) , DIMENSION(Ldgnum,*) :: Difr
      REAL , INTENT(IN) , DIMENSION(*) :: Z
      INTEGER :: K
      REAL :: C
      REAL :: S
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLALS0
   END INTERFACE
END MODULE S_CLALS0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLALSA
   INTERFACE
      SUBROUTINE CLALSA(Icompq,Smlsiz,N,Nrhs,B,Ldb,Bx,Ldbx,U,Ldu,Vt,K,  &
     &                  Difl,Difr,Z,Poles,Givptr,Givcol,Ldgcol,Perm,    &
     &                  Givnum,C,S,Rwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      INTEGER :: Icompq
      INTEGER :: Smlsiz
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldbx,*) :: Bx
      INTEGER :: Ldbx
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldu,*) :: Vt
      INTEGER , DIMENSION(*) :: K
      REAL , DIMENSION(Ldu,*) :: Difl
      REAL , DIMENSION(Ldu,*) :: Difr
      REAL , DIMENSION(Ldu,*) :: Z
      REAL , DIMENSION(Ldu,*) :: Poles
      INTEGER , DIMENSION(*) :: Givptr
      INTEGER , DIMENSION(Ldgcol,*) :: Givcol
      INTEGER :: Ldgcol
      INTEGER , DIMENSION(Ldgcol,*) :: Perm
      REAL , DIMENSION(Ldu,*) :: Givnum
      REAL , DIMENSION(*) :: C
      REAL , DIMENSION(*) :: S
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLALSA
   END INTERFACE
END MODULE S_CLALSA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLALSD
   INTERFACE
      SUBROUTINE CLALSD(Uplo,Smlsiz,N,Nrhs,D,E,B,Ldb,Rcond,Rank,Work,   &
     &                  Rwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0)
      CHARACTER , INTENT(IN) :: Uplo
      INTEGER :: Smlsiz
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      COMPLEX , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLALSD
   END INTERFACE
END MODULE S_CLALSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAMSWLQ
   INTERFACE
      SUBROUTINE CLAMSWLQ(Side,Trans,M,N,K,Mb,Nb,A,Lda,T,Ldt,C,Ldc,Work,&
     &                    Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: Mb
      INTEGER :: Nb
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLAMSWLQ
   END INTERFACE
END MODULE S_CLAMSWLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAMTSQR
   INTERFACE
      SUBROUTINE CLAMTSQR(Side,Trans,M,N,K,Mb,Nb,A,Lda,T,Ldt,C,Ldc,Work,&
     &                    Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: Mb
      INTEGER :: Nb
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLAMTSQR
   END INTERFACE
END MODULE S_CLAMTSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANGB
   INTERFACE
      FUNCTION CLANGB(Norm,N,Kl,Ku,Ab,Ldab,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANGB
      CHARACTER :: Norm
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLANGB
   END INTERFACE
END MODULE S_CLANGB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANGE
   INTERFACE
      FUNCTION CLANGE(Norm,M,N,A,Lda,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANGE
      CHARACTER :: Norm
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLANGE
   END INTERFACE
END MODULE S_CLANGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANGT
   INTERFACE
      FUNCTION CLANGT(Norm,N,Dl,D,Du)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANGT
      CHARACTER :: Norm
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Dl
      COMPLEX , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: Du
      END FUNCTION CLANGT
   END INTERFACE
END MODULE S_CLANGT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANHB
   INTERFACE
      FUNCTION CLANHB(Norm,Uplo,N,K,Ab,Ldab,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANHB
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLANHB
   END INTERFACE
END MODULE S_CLANHB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANHE
   INTERFACE
      FUNCTION CLANHE(Norm,Uplo,N,A,Lda,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANHE
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLANHE
   END INTERFACE
END MODULE S_CLANHE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANHF
   INTERFACE
      FUNCTION CLANHF(Norm,Transr,Uplo,N,A,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANHF
      CHARACTER :: Norm
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(0:*) :: A
      REAL , INTENT(INOUT) , DIMENSION(0:*) :: Work
      END FUNCTION CLANHF
   END INTERFACE
END MODULE S_CLANHF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANHP
   INTERFACE
      FUNCTION CLANHP(Norm,Uplo,N,Ap,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANHP
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(*) :: Ap
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLANHP
   END INTERFACE
END MODULE S_CLANHP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANHS
   INTERFACE
      FUNCTION CLANHS(Norm,N,A,Lda,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANHS
      CHARACTER :: Norm
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLANHS
   END INTERFACE
END MODULE S_CLANHS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANHT
   INTERFACE
      FUNCTION CLANHT(Norm,N,D,E)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANHT
      CHARACTER :: Norm
      INTEGER :: N
      REAL , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: E
      END FUNCTION CLANHT
   END INTERFACE
END MODULE S_CLANHT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANSB
   INTERFACE
      FUNCTION CLANSB(Norm,Uplo,N,K,Ab,Ldab,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANSB
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLANSB
   END INTERFACE
END MODULE S_CLANSB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANSP
   INTERFACE
      FUNCTION CLANSP(Norm,Uplo,N,Ap,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANSP
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(*) :: Ap
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLANSP
   END INTERFACE
END MODULE S_CLANSP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANSY
   INTERFACE
      FUNCTION CLANSY(Norm,Uplo,N,A,Lda,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANSY
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLANSY
   END INTERFACE
END MODULE S_CLANSY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANTB
   INTERFACE
      FUNCTION CLANTB(Norm,Uplo,Diag,N,K,Ab,Ldab,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANTB
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLANTB
   END INTERFACE
END MODULE S_CLANTB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANTP
   INTERFACE
      FUNCTION CLANTP(Norm,Uplo,Diag,N,Ap,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANTP
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(*) :: Ap
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLANTP
   END INTERFACE
END MODULE S_CLANTP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLANTR
   INTERFACE
      FUNCTION CLANTR(Norm,Uplo,Diag,M,N,A,Lda,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: CLANTR
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLANTR
   END INTERFACE
END MODULE S_CLANTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAPLL
   INTERFACE
      SUBROUTINE CLAPLL(N,X,Incx,Y,Incy,Ssmin)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER :: Incx
      COMPLEX , DIMENSION(*) :: Y
      INTEGER :: Incy
      REAL :: Ssmin
      END SUBROUTINE CLAPLL
   END INTERFACE
END MODULE S_CLAPLL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAPMR
   INTERFACE
      SUBROUTINE CLAPMR(Forwrd,M,N,X,Ldx,K)
      IMPLICIT NONE
      LOGICAL , INTENT(IN) :: Forwrd
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: K
      END SUBROUTINE CLAPMR
   END INTERFACE
END MODULE S_CLAPMR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAPMT
   INTERFACE
      SUBROUTINE CLAPMT(Forwrd,M,N,X,Ldx,K)
      IMPLICIT NONE
      LOGICAL , INTENT(IN) :: Forwrd
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: K
      END SUBROUTINE CLAPMT
   END INTERFACE
END MODULE S_CLAPMT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_PORCOND_C
   INTERFACE
      FUNCTION CLA_PORCOND_C(Uplo,N,A,Lda,Af,Ldaf,C,Capply,Info,Work,   &
     &                       Rwork)
      IMPLICIT NONE
      REAL :: CLA_PORCOND_C
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      REAL , INTENT(IN) , DIMENSION(*) :: C
      LOGICAL , INTENT(IN) :: Capply
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION CLA_PORCOND_C
   END INTERFACE
END MODULE S_CLA_PORCOND_C
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_PORCOND_X
   INTERFACE
      FUNCTION CLA_PORCOND_X(Uplo,N,A,Lda,Af,Ldaf,X,Info,Work,Rwork)
      IMPLICIT NONE
      REAL :: CLA_PORCOND_X
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION CLA_PORCOND_X
   END INTERFACE
END MODULE S_CLA_PORCOND_X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_PORFSX_EXTENDED
   INTERFACE
      SUBROUTINE CLA_PORFSX_EXTENDED(Prec_type,Uplo,N,Nrhs,A,Lda,Af,    &
     &                               Ldaf,Colequ,C,B,Ldb,Y,Ldy,Berr_out,&
     &                               N_norms,Err_bnds_norm,             &
     &                               Err_bnds_comp,Res,Ayb,Dy,Y_tail,   &
     &                               Rcond,Ithresh,Rthresh,Dz_ub,       &
     &                               Ignore_cwise,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  UNSTABLE_STATE = 0 , WORKING_STATE = 1 , &
     &                         CONV_STATE = 2 , NOPROG_STATE = 3 ,      &
     &                         BASE_RESIDUAL = 0 , EXTRA_RESIDUAL = 1 , &
     &                         EXTRA_Y = 2 , LA_LINRX_ERR_I = 2
      INTEGER :: Prec_type
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      LOGICAL , INTENT(IN) :: Colequ
      REAL , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      COMPLEX , DIMENSION(*) :: Res
      REAL , DIMENSION(*) :: Ayb
      COMPLEX , DIMENSION(*) :: Dy
      COMPLEX , DIMENSION(*) :: Y_tail
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL , INTENT(IN) :: Rthresh
      REAL , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
      END SUBROUTINE CLA_PORFSX_EXTENDED
   END INTERFACE
END MODULE S_CLA_PORFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_PORPVGRW
   INTERFACE
      FUNCTION CLA_PORPVGRW(Uplo,Ncols,A,Lda,Af,Ldaf,Work)
      IMPLICIT NONE
      REAL :: CLA_PORPVGRW
      CHARACTER(1) :: Uplo
      INTEGER , INTENT(IN) :: Ncols
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLA_PORPVGRW
   END INTERFACE
END MODULE S_CLA_PORPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQGB
   INTERFACE
      SUBROUTINE CLAQGB(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(IN) , DIMENSION(*) :: R
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , INTENT(IN) :: Rowcnd
      REAL , INTENT(IN) :: Colcnd
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE CLAQGB
   END INTERFACE
END MODULE S_CLAQGB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQGE
   INTERFACE
      SUBROUTINE CLAQGE(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: R
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , INTENT(IN) :: Rowcnd
      REAL , INTENT(IN) :: Colcnd
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE CLAQGE
   END INTERFACE
END MODULE S_CLAQGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQHB
   INTERFACE
      SUBROUTINE CLAQHB(Uplo,N,Kd,Ab,Ldab,S,Scond,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(IN) , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Scond
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE CLAQHB
   END INTERFACE
END MODULE S_CLAQHB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQHE
   INTERFACE
      SUBROUTINE CLAQHE(Uplo,N,A,Lda,S,Scond,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Scond
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE CLAQHE
   END INTERFACE
END MODULE S_CLAQHE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQHP
   INTERFACE
      SUBROUTINE CLAQHP(Uplo,N,Ap,S,Scond,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL , INTENT(IN) , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Scond
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE CLAQHP
   END INTERFACE
END MODULE S_CLAQHP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQP2
   INTERFACE
      SUBROUTINE CLAQP2(M,N,Offset,A,Lda,Jpvt,Tau,Vn1,Vn2,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Offset
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      COMPLEX , DIMENSION(*) :: Tau
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vn1
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vn2
      COMPLEX , DIMENSION(*) :: Work
      END SUBROUTINE CLAQP2
   END INTERFACE
END MODULE S_CLAQP2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQPS
   INTERFACE
      SUBROUTINE CLAQPS(M,N,Offset,Nb,Kb,A,Lda,Jpvt,Tau,Vn1,Vn2,Auxv,F, &
     &                  Ldf)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Offset
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Kb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      COMPLEX , DIMENSION(*) :: Tau
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vn1
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vn2
      COMPLEX , DIMENSION(*) :: Auxv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      END SUBROUTINE CLAQPS
   END INTERFACE
END MODULE S_CLAQPS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQR0
   INTERFACE
      SUBROUTINE CLAQR0(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,W,Iloz,Ihiz,Z,Ldz,  &
     &                  Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NTINY = 15 , KEXNW = 5 , KEXSH = 6
      REAL , PARAMETER  ::  WILK1 = 0.75E0
      COMPLEX , PARAMETER  ::  ZERO = (0.0E0,0.0E0) ,                   &
     &                         ONE = (1.0E0,0.0E0)
      REAL , PARAMETER  ::  TWO = 2.0E0
      LOGICAL :: Wantt
      LOGICAL :: Wantz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER :: Iloz
      INTEGER :: Ihiz
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER :: Info
      END SUBROUTINE CLAQR0
   END INTERFACE
END MODULE S_CLAQR0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQR1
   INTERFACE
      SUBROUTINE CLAQR1(N,H,Ldh,S1,S2,V)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E0,0.0E0)
      REAL , PARAMETER  ::  RZERO = 0.0E0
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Ldh,*) :: H
      INTEGER , INTENT(IN) :: Ldh
      COMPLEX , INTENT(IN) :: S1
      COMPLEX , INTENT(IN) :: S2
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: V
      END SUBROUTINE CLAQR1
   END INTERFACE
END MODULE S_CLAQR1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQR2
   INTERFACE
      SUBROUTINE CLAQR2(Wantt,Wantz,N,Ktop,Kbot,Nw,H,Ldh,Iloz,Ihiz,Z,   &
     &                  Ldz,Ns,Nd,Sh,V,Ldv,Nh,T,Ldt,Nv,Wv,Ldwv,Work,    &
     &                  Lwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E0,0.0E0) ,                   &
     &                         ONE = (1.0E0,0.0E0)
      REAL , PARAMETER  ::  RZERO = 0.0E0 , RONE = 1.0E0
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nw
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ns
      INTEGER , INTENT(OUT) :: Nd
      COMPLEX , DIMENSION(*) :: Sh
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(IN) :: Nh
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(IN) :: Nv
      COMPLEX , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      END SUBROUTINE CLAQR2
   END INTERFACE
END MODULE S_CLAQR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQR3
   INTERFACE
      SUBROUTINE CLAQR3(Wantt,Wantz,N,Ktop,Kbot,Nw,H,Ldh,Iloz,Ihiz,Z,   &
     &                  Ldz,Ns,Nd,Sh,V,Ldv,Nh,T,Ldt,Nv,Wv,Ldwv,Work,    &
     &                  Lwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E0,0.0E0) ,                   &
     &                         ONE = (1.0E0,0.0E0)
      REAL , PARAMETER  ::  RZERO = 0.0E0 , RONE = 1.0E0
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nw
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ns
      INTEGER , INTENT(OUT) :: Nd
      COMPLEX , DIMENSION(*) :: Sh
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(IN) :: Nh
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(IN) :: Nv
      COMPLEX , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      END SUBROUTINE CLAQR3
   END INTERFACE
END MODULE S_CLAQR3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQR4
   INTERFACE
      SUBROUTINE CLAQR4(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,W,Iloz,Ihiz,Z,Ldz,  &
     &                  Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NTINY = 15 , KEXNW = 5 , KEXSH = 6
      REAL , PARAMETER  ::  WILK1 = 0.75E0
      COMPLEX , PARAMETER  ::  ZERO = (0.0E0,0.0E0) ,                   &
     &                         ONE = (1.0E0,0.0E0)
      REAL , PARAMETER  ::  TWO = 2.0E0
      LOGICAL :: Wantt
      LOGICAL :: Wantz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER :: Iloz
      INTEGER :: Ihiz
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER :: Info
      END SUBROUTINE CLAQR4
   END INTERFACE
END MODULE S_CLAQR4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQR5
   INTERFACE
      SUBROUTINE CLAQR5(Wantt,Wantz,Kacc22,N,Ktop,Kbot,Nshfts,S,H,Ldh,  &
     &                  Iloz,Ihiz,Z,Ldz,V,Ldv,U,Ldu,Nv,Wv,Ldwv,Nh,Wh,   &
     &                  Ldwh)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E0,0.0E0) ,                   &
     &                         ONE = (1.0E0,0.0E0)
      REAL , PARAMETER  ::  RZERO = 0.0E0 , RONE = 1.0E0
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: Kacc22
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nshfts
      COMPLEX , DIMENSION(*) :: S
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      INTEGER , INTENT(IN) :: Nv
      COMPLEX , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      INTEGER , INTENT(IN) :: Nh
      COMPLEX , DIMENSION(Ldwh,*) :: Wh
      INTEGER :: Ldwh
      END SUBROUTINE CLAQR5
   END INTERFACE
END MODULE S_CLAQR5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQSB
   INTERFACE
      SUBROUTINE CLAQSB(Uplo,N,Kd,Ab,Ldab,S,Scond,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(IN) , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Scond
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE CLAQSB
   END INTERFACE
END MODULE S_CLAQSB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQSP
   INTERFACE
      SUBROUTINE CLAQSP(Uplo,N,Ap,S,Scond,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL , INTENT(IN) , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Scond
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE CLAQSP
   END INTERFACE
END MODULE S_CLAQSP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAQSY
   INTERFACE
      SUBROUTINE CLAQSY(Uplo,N,A,Lda,S,Scond,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Scond
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE CLAQSY
   END INTERFACE
END MODULE S_CLAQSY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAR1V
   INTERFACE
      SUBROUTINE CLAR1V(N,B1,Bn,Lambda,D,L,Ld,Lld,Pivmin,Gaptol,Z,      &
     &                  Wantnc,Negcnt,Ztz,Mingma,R,Isuppz,Nrminv,Resid, &
     &                  Rqcorr,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CONE = (1.0E0,0.0E0)
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: B1
      INTEGER , INTENT(IN) :: Bn
      REAL , INTENT(IN) :: Lambda
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: L
      REAL , INTENT(IN) , DIMENSION(*) :: Ld
      REAL , INTENT(IN) , DIMENSION(*) :: Lld
      REAL , INTENT(IN) :: Pivmin
      REAL , INTENT(IN) :: Gaptol
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Z
      LOGICAL , INTENT(IN) :: Wantnc
      INTEGER , INTENT(OUT) :: Negcnt
      REAL , INTENT(INOUT) :: Ztz
      REAL , INTENT(INOUT) :: Mingma
      INTEGER , INTENT(INOUT) :: R
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Isuppz
      REAL , INTENT(INOUT) :: Nrminv
      REAL , INTENT(OUT) :: Resid
      REAL , INTENT(OUT) :: Rqcorr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END SUBROUTINE CLAR1V
   END INTERFACE
END MODULE S_CLAR1V
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAR2V
   INTERFACE
      SUBROUTINE CLAR2V(N,X,Y,Z,Incx,C,S,Incc)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Y
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX , INTENT(IN) , DIMENSION(*) :: S
      INTEGER , INTENT(IN) :: Incc
      END SUBROUTINE CLAR2V
   END INTERFACE
END MODULE S_CLAR2V
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARCM
   INTERFACE
      SUBROUTINE CLARCM(M,N,A,Lda,B,Ldb,C,Ldc,Rwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , ZERO = 0.0E0
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END SUBROUTINE CLARCM
   END INTERFACE
END MODULE S_CLARCM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARFB
   INTERFACE
      SUBROUTINE CLARFB(Side,Trans,Direct,Storev,M,N,K,V,Ldv,T,Ldt,C,   &
     &                  Ldc,Work,Ldwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE CLARFB
   END INTERFACE
END MODULE S_CLARFB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARFB_GETT
   INTERFACE
      SUBROUTINE CLARFB_GETT(Ident,M,N,K,T,Ldt,A,Lda,B,Ldb,Work,Ldwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Ident
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE CLARFB_GETT
   END INTERFACE
END MODULE S_CLARFB_GETT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARF
   INTERFACE
      SUBROUTINE CLARF(Side,M,N,V,Incv,Tau,C,Ldc,Work)
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
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      END SUBROUTINE CLARF
   END INTERFACE
END MODULE S_CLARF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARFG
   INTERFACE
      SUBROUTINE CLARFG(N,Alpha,X,Incx,Tau)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) :: Alpha
      COMPLEX , DIMENSION(*) :: X
      INTEGER :: Incx
      COMPLEX , INTENT(OUT) :: Tau
      END SUBROUTINE CLARFG
   END INTERFACE
END MODULE S_CLARFG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARFGP
   INTERFACE
      SUBROUTINE CLARFGP(N,Alpha,X,Incx,Tau)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  TWO = 2.0E+0 , ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) :: Alpha
      COMPLEX , DIMENSION(*) :: X
      INTEGER :: Incx
      COMPLEX , INTENT(INOUT) :: Tau
      END SUBROUTINE CLARFGP
   END INTERFACE
END MODULE S_CLARFGP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARFT
   INTERFACE
      SUBROUTINE CLARFT(Direct,Storev,N,K,V,Ldv,Tau,T,Ldt)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      END SUBROUTINE CLARFT
   END INTERFACE
END MODULE S_CLARFT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARFX
   INTERFACE
      SUBROUTINE CLARFX(Side,M,N,V,Tau,C,Ldc,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: V
      COMPLEX :: Tau
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      END SUBROUTINE CLARFX
   END INTERFACE
END MODULE S_CLARFX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARFY
   INTERFACE
      SUBROUTINE CLARFY(Uplo,N,V,Incv,Tau,C,Ldc,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         HALF = (0.5E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: V
      INTEGER :: Incv
      COMPLEX :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      END SUBROUTINE CLARFY
   END INTERFACE
END MODULE S_CLARFY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARGV
   INTERFACE
      SUBROUTINE CLARGV(N,X,Incx,Y,Incy,C,Incc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  TWO = 2.0E+0 , ONE = 1.0E+0 , ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(OUT) , DIMENSION(*) :: C
      INTEGER , INTENT(IN) :: Incc
      END SUBROUTINE CLARGV
   END INTERFACE
END MODULE S_CLARGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARNV
   INTERFACE
      SUBROUTINE CLARNV(Idist,Iseed,N,X)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      INTEGER , PARAMETER  ::  LV = 128
      REAL , PARAMETER  ::  TWOPI =                                     &
     &                      6.28318530717958647692528676655900576839E+0
      INTEGER , INTENT(IN) :: Idist
      INTEGER , DIMENSION(4) :: Iseed
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: X
      END SUBROUTINE CLARNV
   END INTERFACE
END MODULE S_CLARNV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARRV
   INTERFACE
      SUBROUTINE CLARRV(N,Vl,Vu,D,L,Pivmin,Isplit,M,Dol,Dou,Minrgp,     &
     &                  Rtol1,Rtol2,W,Werr,Wgap,Iblock,Indexw,Gers,Z,   &
     &                  Ldz,Isuppz,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXITR = 10
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0)
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      THREE = 3.0E0 , FOUR = 4.0E0 , HALF = 0.5E0
      INTEGER :: N
      REAL , INTENT(IN) :: Vl
      REAL :: Vu
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: L
      REAL :: Pivmin
      INTEGER , INTENT(IN) , DIMENSION(*) :: Isplit
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: Dol
      INTEGER , INTENT(IN) :: Dou
      REAL , INTENT(IN) :: Minrgp
      REAL :: Rtol1
      REAL :: Rtol2
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Wgap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(IN) , DIMENSION(*) :: Indexw
      REAL , INTENT(IN) , DIMENSION(*) :: Gers
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isuppz
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE CLARRV
   END INTERFACE
END MODULE S_CLARRV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARSCL2
   INTERFACE
      SUBROUTINE CLARSCL2(M,N,D,X,Ldx)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      END SUBROUTINE CLARSCL2
   END INTERFACE
END MODULE S_CLARSCL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARTG
   INTERFACE
      SUBROUTINE CLARTG(F,G,Cs,Sn,R)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  TWO = 2.0E+0 , ONE = 1.0E+0 , ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
      COMPLEX , INTENT(IN) :: F
      COMPLEX , INTENT(IN) :: G
      REAL , INTENT(INOUT) :: Cs
      COMPLEX , INTENT(INOUT) :: Sn
      COMPLEX , INTENT(INOUT) :: R
      END SUBROUTINE CLARTG
   END INTERFACE
END MODULE S_CLARTG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARTV
   INTERFACE
      SUBROUTINE CLARTV(N,X,Incx,Y,Incy,C,S,Incc)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX , INTENT(IN) , DIMENSION(*) :: S
      INTEGER , INTENT(IN) :: Incc
      END SUBROUTINE CLARTV
   END INTERFACE
END MODULE S_CLARTV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARZB
   INTERFACE
      SUBROUTINE CLARZB(Side,Trans,Direct,Storev,M,N,K,L,V,Ldv,T,Ldt,C, &
     &                  Ldc,Work,Ldwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: L
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE CLARZB
   END INTERFACE
END MODULE S_CLARZB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARZ
   INTERFACE
      SUBROUTINE CLARZ(Side,M,N,L,V,Incv,Tau,C,Ldc,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      INTEGER :: L
      COMPLEX , DIMENSION(*) :: V
      INTEGER :: Incv
      COMPLEX :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      END SUBROUTINE CLARZ
   END INTERFACE
END MODULE S_CLARZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLARZT
   INTERFACE
      SUBROUTINE CLARZT(Direct,Storev,N,K,V,Ldv,Tau,T,Ldt)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      END SUBROUTINE CLARZT
   END INTERFACE
END MODULE S_CLARZT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLASCL2
   INTERFACE
      SUBROUTINE CLASCL2(M,N,D,X,Ldx)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      END SUBROUTINE CLASCL2
   END INTERFACE
END MODULE S_CLASCL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLASCL
   INTERFACE
      SUBROUTINE CLASCL(Type,Kl,Ku,Cfrom,Cto,M,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Type
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL :: Cfrom
      REAL :: Cto
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLASCL
   END INTERFACE
END MODULE S_CLASCL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLASET
   INTERFACE
      SUBROUTINE CLASET(Uplo,M,N,Alpha,Beta,A,Lda)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(OUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE CLASET
   END INTERFACE
END MODULE S_CLASET
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLASR
   INTERFACE
      SUBROUTINE CLASR(Side,Pivot,Direct,M,N,C,S,A,Lda)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Side
      CHARACTER :: Pivot
      CHARACTER :: Direct
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , INTENT(IN) , DIMENSION(*) :: S
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE CLASR
   END INTERFACE
END MODULE S_CLASR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLASSQ
   INTERFACE
      SUBROUTINE CLASSQ(N,X,Incx,Scale,Sumsq)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(INOUT) :: Sumsq
      END SUBROUTINE CLASSQ
   END INTERFACE
END MODULE S_CLASSQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLASWLQ
   INTERFACE
      SUBROUTINE CLASWLQ(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER :: Nb
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLASWLQ
   END INTERFACE
END MODULE S_CLASWLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLASWP
   INTERFACE
      SUBROUTINE CLASWP(N,A,Lda,K1,K2,Ipiv,Incx)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) :: K1
      INTEGER , INTENT(IN) :: K2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE CLASWP
   END INTERFACE
END MODULE S_CLASWP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_SYAMV
   INTERFACE
      SUBROUTINE CLA_SYAMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE CLA_SYAMV
   END INTERFACE
END MODULE S_CLA_SYAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLASYF_AA
   INTERFACE
      SUBROUTINE CLASYF_AA(Uplo,J1,M,Nb,A,Lda,Ipiv,H,Ldh,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: J1
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      END SUBROUTINE CLASYF_AA
   END INTERFACE
END MODULE S_CLASYF_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLASYF
   INTERFACE
      SUBROUTINE CLASYF(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLASYF
   END INTERFACE
END MODULE S_CLASYF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLASYF_RK
   INTERFACE
      SUBROUTINE CLASYF_RK(Uplo,N,Nb,Kb,A,Lda,E,Ipiv,W,Ldw,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLASYF_RK
   END INTERFACE
END MODULE S_CLASYF_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLASYF_ROOK
   INTERFACE
      SUBROUTINE CLASYF_ROOK(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLASYF_ROOK
   END INTERFACE
END MODULE S_CLASYF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_SYRCOND_C
   INTERFACE
      FUNCTION CLA_SYRCOND_C(Uplo,N,A,Lda,Af,Ldaf,Ipiv,C,Capply,Info,   &
     &                       Work,Rwork)
      IMPLICIT NONE
      REAL :: CLA_SYRCOND_C
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) , DIMENSION(*) :: C
      LOGICAL , INTENT(IN) :: Capply
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION CLA_SYRCOND_C
   END INTERFACE
END MODULE S_CLA_SYRCOND_C
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_SYRCOND_X
   INTERFACE
      FUNCTION CLA_SYRCOND_X(Uplo,N,A,Lda,Af,Ldaf,Ipiv,X,Info,Work,     &
     &                       Rwork)
      IMPLICIT NONE
      REAL :: CLA_SYRCOND_X
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION CLA_SYRCOND_X
   END INTERFACE
END MODULE S_CLA_SYRCOND_X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_SYRFSX_EXTENDED
   INTERFACE
      SUBROUTINE CLA_SYRFSX_EXTENDED(Prec_type,Uplo,N,Nrhs,A,Lda,Af,    &
     &                               Ldaf,Ipiv,Colequ,C,B,Ldb,Y,Ldy,    &
     &                               Berr_out,N_norms,Err_bnds_norm,    &
     &                               Err_bnds_comp,Res,Ayb,Dy,Y_tail,   &
     &                               Rcond,Ithresh,Rthresh,Dz_ub,       &
     &                               Ignore_cwise,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  UNSTABLE_STATE = 0 , WORKING_STATE = 1 , &
     &                         CONV_STATE = 2 , NOPROG_STATE = 3 ,      &
     &                         BASE_RESIDUAL = 0 , EXTRA_RESIDUAL = 1 , &
     &                         EXTRA_Y = 2 , LA_LINRX_ERR_I = 2
      INTEGER :: Prec_type
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      COMPLEX , DIMENSION(*) :: Res
      REAL , DIMENSION(*) :: Ayb
      COMPLEX , DIMENSION(*) :: Dy
      COMPLEX , DIMENSION(*) :: Y_tail
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL , INTENT(IN) :: Rthresh
      REAL , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLA_SYRFSX_EXTENDED
   END INTERFACE
END MODULE S_CLA_SYRFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_SYRPVGRW
   INTERFACE
      FUNCTION CLA_SYRPVGRW(Uplo,N,Info,A,Lda,Af,Ldaf,Ipiv,Work)
      IMPLICIT NONE
      REAL :: CLA_SYRPVGRW
      CHARACTER(1) :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Info
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION CLA_SYRPVGRW
   END INTERFACE
END MODULE S_CLA_SYRPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLATBS
   INTERFACE
      SUBROUTINE CLATBS(Uplo,Trans,Diag,Normin,N,Kd,Ab,Ldab,X,Scale,    &
     &                  Cnorm,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , HALF = 0.5E+0 ,             &
     &                      ONE = 1.0E+0 , TWO = 2.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLATBS
   END INTERFACE
END MODULE S_CLATBS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLATDF
   INTERFACE
      SUBROUTINE CLATDF(Ijob,N,Z,Ldz,Rhs,Rdsum,Rdscal,Ipiv,Jpiv)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXDIM = 2
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: Ijob
      INTEGER :: N
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Rhs
      REAL :: Rdsum
      REAL :: Rdscal
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Jpiv
      END SUBROUTINE CLATDF
   END INTERFACE
END MODULE S_CLATDF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLATPS
   INTERFACE
      SUBROUTINE CLATPS(Uplo,Trans,Diag,Normin,N,Ap,X,Scale,Cnorm,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , HALF = 0.5E+0 ,             &
     &                      ONE = 1.0E+0 , TWO = 2.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLATPS
   END INTERFACE
END MODULE S_CLATPS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLATRD
   INTERFACE
      SUBROUTINE CLATRD(Uplo,N,Nb,A,Lda,E,Tau,W,Ldw)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0) ,                  &
     &                         HALF = (0.5E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      END SUBROUTINE CLATRD
   END INTERFACE
END MODULE S_CLATRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLATRS
   INTERFACE
      SUBROUTINE CLATRS(Uplo,Trans,Diag,Normin,N,A,Lda,X,Scale,Cnorm,   &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , HALF = 0.5E+0 ,             &
     &                      ONE = 1.0E+0 , TWO = 2.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLATRS
   END INTERFACE
END MODULE S_CLATRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLATRZ
   INTERFACE
      SUBROUTINE CLATRZ(M,N,L,A,Lda,Tau,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: L
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      END SUBROUTINE CLATRZ
   END INTERFACE
END MODULE S_CLATRZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLATSQR
   INTERFACE
      SUBROUTINE CLATSQR(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER :: Nb
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLATSQR
   END INTERFACE
END MODULE S_CLATSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAUNHR_COL_GETRFNP2
   INTERFACE
      RECURSIVE SUBROUTINE CLAUNHR_COL_GETRFNP2(M,N,A,Lda,D,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLAUNHR_COL_GETRFNP2
   END INTERFACE
END MODULE S_CLAUNHR_COL_GETRFNP2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAUNHR_COL_GETRFNP
   INTERFACE
      SUBROUTINE CLAUNHR_COL_GETRFNP(M,N,A,Lda,D,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLAUNHR_COL_GETRFNP
   END INTERFACE
END MODULE S_CLAUNHR_COL_GETRFNP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAUU2
   INTERFACE
      SUBROUTINE CLAUU2(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CLAUU2
   END INTERFACE
END MODULE S_CLAUU2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLAUUM
   INTERFACE
      SUBROUTINE CLAUUM(Uplo,N,A,Lda,Info)
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
      END SUBROUTINE CLAUUM
   END INTERFACE
END MODULE S_CLAUUM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CLA_WWADDW
   INTERFACE
      SUBROUTINE CLA_WWADDW(N,X,Y,W)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: X
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Y
      COMPLEX , INTENT(IN) , DIMENSION(*) :: W
      END SUBROUTINE CLA_WWADDW
   END INTERFACE
END MODULE S_CLA_WWADDW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPBCON
   INTERFACE
      SUBROUTINE CPBCON(Uplo,N,Kd,Ab,Ldab,Anorm,Rcond,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPBCON
   END INTERFACE
END MODULE S_CPBCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPBEQU
   INTERFACE
      SUBROUTINE CPBEQU(Uplo,N,Kd,Ab,Ldab,S,Scond,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) :: Scond
      REAL , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPBEQU
   END INTERFACE
END MODULE S_CPBEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPBRFS
   INTERFACE
      SUBROUTINE CPBRFS(Uplo,N,Kd,Nrhs,Ab,Ldab,Afb,Ldafb,B,Ldb,X,Ldx,   &
     &                  Ferr,Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPBRFS
   END INTERFACE
END MODULE S_CPBRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPBSTF
   INTERFACE
      SUBROUTINE CPBSTF(Uplo,N,Kd,Ab,Ldab,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPBSTF
   END INTERFACE
END MODULE S_CPBSTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPBSV
   INTERFACE
      SUBROUTINE CPBSV(Uplo,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPBSV
   END INTERFACE
END MODULE S_CPBSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPBSVX
   INTERFACE
      SUBROUTINE CPBSVX(Fact,Uplo,N,Kd,Nrhs,Ab,Ldab,Afb,Ldafb,Equed,S,B,&
     &                  Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: S
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPBSVX
   END INTERFACE
END MODULE S_CPBSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPBTF2
   INTERFACE
      SUBROUTINE CPBTF2(Uplo,N,Kd,Ab,Ldab,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPBTF2
   END INTERFACE
END MODULE S_CPBTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPBTRF
   INTERFACE
      SUBROUTINE CPBTRF(Uplo,N,Kd,Ab,Ldab,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      INTEGER , PARAMETER  ::  NBMAX = 32 , LDWORK = NBMAX + 1
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPBTRF
   END INTERFACE
END MODULE S_CPBTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPBTRS
   INTERFACE
      SUBROUTINE CPBTRS(Uplo,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPBTRS
   END INTERFACE
END MODULE S_CPBTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPFTRF
   INTERFACE
      SUBROUTINE CPFTRF(Transr,Uplo,N,A,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPFTRF
   END INTERFACE
END MODULE S_CPFTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPFTRI
   INTERFACE
      SUBROUTINE CPFTRI(Transr,Uplo,N,A,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPFTRI
   END INTERFACE
END MODULE S_CPFTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPFTRS
   INTERFACE
      SUBROUTINE CPFTRS(Transr,Uplo,N,Nrhs,A,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(0:*) :: A
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPFTRS
   END INTERFACE
END MODULE S_CPFTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPOCON
   INTERFACE
      SUBROUTINE CPOCON(Uplo,N,A,Lda,Anorm,Rcond,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPOCON
   END INTERFACE
END MODULE S_CPOCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPOEQUB
   INTERFACE
      SUBROUTINE CPOEQUB(N,A,Lda,S,Scond,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) :: Scond
      REAL , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPOEQUB
   END INTERFACE
END MODULE S_CPOEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPOEQU
   INTERFACE
      SUBROUTINE CPOEQU(N,A,Lda,S,Scond,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) :: Scond
      REAL , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPOEQU
   END INTERFACE
END MODULE S_CPOEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPORFS
   INTERFACE
      SUBROUTINE CPORFS(Uplo,N,Nrhs,A,Lda,Af,Ldaf,B,Ldb,X,Ldx,Ferr,Berr,&
     &                  Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPORFS
   END INTERFACE
END MODULE S_CPORFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPORFSX
   INTERFACE
      SUBROUTINE CPORFSX(Uplo,Equed,N,Nrhs,A,Lda,Af,Ldaf,S,B,Ldb,X,Ldx, &
     &                   Rcond,Berr,N_err_bnds,Err_bnds_norm,           &
     &                   Err_bnds_comp,Nparams,Params,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ITREF_DEFAULT = 1.0 ,       &
     &                      ITHRESH_DEFAULT = 10.0 ,                    &
     &                      COMPONENTWISE_DEFAULT = 1.0 ,               &
     &                      RTHRESH_DEFAULT = 0.5 ,                     &
     &                      DZTHRESH_DEFAULT = 0.25
      INTEGER , PARAMETER  ::  LA_LINRX_ITREF_I = 1 ,                   &
     &                         LA_LINRX_ITHRESH_I = 2 ,                 &
     &                         LA_LINRX_CWISE_I = 3 ,                   &
     &                         LA_LINRX_TRUST_I = 1 ,                   &
     &                         LA_LINRX_ERR_I = 2 , LA_LINRX_RCOND_I = 3
      CHARACTER :: Uplo
      CHARACTER :: Equed
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      REAL , DIMENSION(*) :: S
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL , INTENT(INOUT) , DIMENSION(*) :: Params
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPORFSX
   END INTERFACE
END MODULE S_CPORFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPOSV
   INTERFACE
      SUBROUTINE CPOSV(Uplo,N,Nrhs,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPOSV
   END INTERFACE
END MODULE S_CPOSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPOSVX
   INTERFACE
      SUBROUTINE CPOSVX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Equed,S,B,Ldb,X, &
     &                  Ldx,Rcond,Ferr,Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: S
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPOSVX
   END INTERFACE
END MODULE S_CPOSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPOSVXX
   INTERFACE
      SUBROUTINE CPOSVXX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Equed,S,B,Ldb,X,&
     &                   Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,Err_bnds_norm,&
     &                   Err_bnds_comp,Nparams,Params,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: S
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , INTENT(OUT) :: Rpvgrw
      REAL , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL , DIMENSION(*) :: Params
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPOSVXX
   END INTERFACE
END MODULE S_CPOSVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPOTF2
   INTERFACE
      SUBROUTINE CPOTF2(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPOTF2
   END INTERFACE
END MODULE S_CPOTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPOTRF2
   INTERFACE
      RECURSIVE SUBROUTINE CPOTRF2(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPOTRF2
   END INTERFACE
END MODULE S_CPOTRF2
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
MODULE S_CPOTRI
   INTERFACE
      SUBROUTINE CPOTRI(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPOTRI
   END INTERFACE
END MODULE S_CPOTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPOTRS
   INTERFACE
      SUBROUTINE CPOTRS(Uplo,N,Nrhs,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPOTRS
   END INTERFACE
END MODULE S_CPOTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPPCON
   INTERFACE
      SUBROUTINE CPPCON(Uplo,N,Ap,Anorm,Rcond,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPPCON
   END INTERFACE
END MODULE S_CPPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPPEQU
   INTERFACE
      SUBROUTINE CPPEQU(Uplo,N,Ap,S,Scond,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Ap
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) :: Scond
      REAL , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPPEQU
   END INTERFACE
END MODULE S_CPPEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPPRFS
   INTERFACE
      SUBROUTINE CPPRFS(Uplo,N,Nrhs,Ap,Afp,B,Ldb,X,Ldx,Ferr,Berr,Work,  &
     &                  Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Afp
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPPRFS
   END INTERFACE
END MODULE S_CPPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPPSV
   INTERFACE
      SUBROUTINE CPPSV(Uplo,N,Nrhs,Ap,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPPSV
   END INTERFACE
END MODULE S_CPPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPPSVX
   INTERFACE
      SUBROUTINE CPPSVX(Fact,Uplo,N,Nrhs,Ap,Afp,Equed,S,B,Ldb,X,Ldx,    &
     &                  Rcond,Ferr,Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Afp
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: S
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPPSVX
   END INTERFACE
END MODULE S_CPPSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPPTRF
   INTERFACE
      SUBROUTINE CPPTRF(Uplo,N,Ap,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPPTRF
   END INTERFACE
END MODULE S_CPPTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPPTRI
   INTERFACE
      SUBROUTINE CPPTRI(Uplo,N,Ap,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPPTRI
   END INTERFACE
END MODULE S_CPPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPPTRS
   INTERFACE
      SUBROUTINE CPPTRS(Uplo,N,Nrhs,Ap,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPPTRS
   END INTERFACE
END MODULE S_CPPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPSTF2
   INTERFACE
      SUBROUTINE CPSTF2(Uplo,N,A,Lda,Piv,Rank,Tol,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(N) :: Piv
      INTEGER , INTENT(OUT) :: Rank
      REAL , INTENT(IN) :: Tol
      REAL , INTENT(INOUT) , DIMENSION(2*N) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPSTF2
   END INTERFACE
END MODULE S_CPSTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPSTRF
   INTERFACE
      SUBROUTINE CPSTRF(Uplo,N,A,Lda,Piv,Rank,Tol,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(N) :: Piv
      INTEGER :: Rank
      REAL :: Tol
      REAL , INTENT(INOUT) , DIMENSION(2*N) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPSTRF
   END INTERFACE
END MODULE S_CPSTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPTCON
   INTERFACE
      SUBROUTINE CPTCON(N,D,E,Anorm,Rcond,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX , INTENT(IN) , DIMENSION(*) :: E
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPTCON
   END INTERFACE
END MODULE S_CPTCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPTEQR
   INTERFACE
      SUBROUTINE CPTEQR(Compz,N,D,E,Z,Ldz,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Compz
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPTEQR
   END INTERFACE
END MODULE S_CPTEQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPTRFS
   INTERFACE
      SUBROUTINE CPTRFS(Uplo,N,Nrhs,D,E,Df,Ef,B,Ldb,X,Ldx,Ferr,Berr,    &
     &                  Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX , INTENT(IN) , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Df
      COMPLEX , DIMENSION(*) :: Ef
      COMPLEX , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPTRFS
   END INTERFACE
END MODULE S_CPTRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPTSV
   INTERFACE
      SUBROUTINE CPTSV(N,Nrhs,D,E,B,Ldb,Info)
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: E
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPTSV
   END INTERFACE
END MODULE S_CPTSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPTSVX
   INTERFACE
      SUBROUTINE CPTSVX(Fact,N,Nrhs,D,E,Df,Ef,B,Ldb,X,Ldx,Rcond,Ferr,   &
     &                  Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Fact
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Df
      COMPLEX , DIMENSION(*) :: Ef
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPTSVX
   END INTERFACE
END MODULE S_CPTSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPTTRF
   INTERFACE
      SUBROUTINE CPTTRF(N,D,E,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER :: Info
      END SUBROUTINE CPTTRF
   END INTERFACE
END MODULE S_CPTTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPTTRS
   INTERFACE
      SUBROUTINE CPTTRS(Uplo,N,Nrhs,D,E,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: E
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CPTTRS
   END INTERFACE
END MODULE S_CPTTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CPTTS2
   INTERFACE
      SUBROUTINE CPTTS2(Iuplo,N,Nrhs,D,E,B,Ldb)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: Iuplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: D
      COMPLEX , INTENT(IN) , DIMENSION(*) :: E
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      END SUBROUTINE CPTTS2
   END INTERFACE
END MODULE S_CPTTS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CROT
   INTERFACE
      SUBROUTINE CROT(N,Cx,Incx,Cy,Incy,C,S)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(IN) :: C
      COMPLEX , INTENT(IN) :: S
      END SUBROUTINE CROT
   END INTERFACE
END MODULE S_CROT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSPCON
   INTERFACE
      SUBROUTINE CSPCON(Uplo,N,Ap,Ipiv,Anorm,Rcond,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSPCON
   END INTERFACE
END MODULE S_CSPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSPMV
   INTERFACE
      SUBROUTINE CSPMV(Uplo,N,Alpha,Ap,X,Incx,Beta,Y,Incy)
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
      END SUBROUTINE CSPMV
   END INTERFACE
END MODULE S_CSPMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSPR
   INTERFACE
      SUBROUTINE CSPR(Uplo,N,Alpha,X,Incx,Ap)
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
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      END SUBROUTINE CSPR
   END INTERFACE
END MODULE S_CSPR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSPRFS
   INTERFACE
      SUBROUTINE CSPRFS(Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Ferr,Berr,  &
     &                  Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSPRFS
   END INTERFACE
END MODULE S_CSPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSPSV
   INTERFACE
      SUBROUTINE CSPSV(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSPSV
   END INTERFACE
END MODULE S_CSPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSPSVX
   INTERFACE
      SUBROUTINE CSPSVX(Fact,Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Rcond, &
     &                  Ferr,Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSPSVX
   END INTERFACE
END MODULE S_CSPSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSPTRF
   INTERFACE
      SUBROUTINE CSPTRF(Uplo,N,Ap,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSPTRF
   END INTERFACE
END MODULE S_CSPTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSPTRI
   INTERFACE
      SUBROUTINE CSPTRI(Uplo,N,Ap,Ipiv,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSPTRI
   END INTERFACE
END MODULE S_CSPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSPTRS
   INTERFACE
      SUBROUTINE CSPTRS(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSPTRS
   END INTERFACE
END MODULE S_CSPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSRSCL
   INTERFACE
      SUBROUTINE CSRSCL(N,Sa,Sx,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: N
      REAL , INTENT(IN) :: Sa
      COMPLEX , DIMENSION(*) :: Sx
      INTEGER :: Incx
      END SUBROUTINE CSRSCL
   END INTERFACE
END MODULE S_CSRSCL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSTEDC
   INTERFACE
      SUBROUTINE CSTEDC(Compz,N,D,E,Z,Ldz,Work,Lwork,Rwork,Lrwork,Iwork,&
     &                  Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
      CHARACTER :: Compz
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSTEDC
   END INTERFACE
END MODULE S_CSTEDC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSTEGR
   INTERFACE
      SUBROUTINE CSTEGR(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,  &
     &                  Isuppz,Work,Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
      CHARACTER :: Jobz
      CHARACTER :: Range
      INTEGER :: N
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL :: Abstol
      INTEGER :: M
      REAL , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER :: Info
      END SUBROUTINE CSTEGR
   END INTERFACE
END MODULE S_CSTEGR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSTEIN
   INTERFACE
      SUBROUTINE CSTEIN(N,D,E,M,W,Iblock,Isplit,Z,Ldz,Work,Iwork,Ifail, &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TEN = 1.0E+1 , ODM3 = 1.0E-3 , ODM1 = 1.0E-1
      INTEGER , PARAMETER  ::  MAXITS = 5 , EXTRA = 2
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) , DIMENSION(*) :: W
      INTEGER , INTENT(IN) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(IN) , DIMENSION(*) :: Isplit
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSTEIN
   END INTERFACE
END MODULE S_CSTEIN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSTEMR
   INTERFACE
      SUBROUTINE CSTEMR(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,M,W,Z,Ldz,Nzc,     &
     &                  Isuppz,Tryrac,Work,Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , FOUR = 4.0E0 , &
     &                      MINRGP = 3.0E-3
      CHARACTER :: Jobz
      CHARACTER :: Range
      INTEGER :: N
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL :: Vl
      REAL :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(IN) :: Nzc
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isuppz
      LOGICAL , INTENT(INOUT) :: Tryrac
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSTEMR
   END INTERFACE
END MODULE S_CSTEMR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSTEQR
   INTERFACE
      SUBROUTINE CSTEQR(Compz,N,D,E,Z,Ldz,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      THREE = 3.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
      INTEGER , PARAMETER  ::  MAXIT = 30
      CHARACTER :: Compz
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSTEQR
   END INTERFACE
END MODULE S_CSTEQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYCON_3
   INTERFACE
      SUBROUTINE CSYCON_3(Uplo,N,A,Lda,E,Ipiv,Anorm,Rcond,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYCON_3
   END INTERFACE
END MODULE S_CSYCON_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYCON
   INTERFACE
      SUBROUTINE CSYCON(Uplo,N,A,Lda,Ipiv,Anorm,Rcond,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYCON
   END INTERFACE
END MODULE S_CSYCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYCON_ROOK
   INTERFACE
      SUBROUTINE CSYCON_ROOK(Uplo,N,A,Lda,Ipiv,Anorm,Rcond,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYCON_ROOK
   END INTERFACE
END MODULE S_CSYCON_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYCONV
   INTERFACE
      SUBROUTINE CSYCONV(Uplo,Way,N,A,Lda,Ipiv,E,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYCONV
   END INTERFACE
END MODULE S_CSYCONV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYCONVF
   INTERFACE
      SUBROUTINE CSYCONVF(Uplo,Way,N,A,Lda,E,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYCONVF
   END INTERFACE
END MODULE S_CSYCONVF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYCONVF_ROOK
   INTERFACE
      SUBROUTINE CSYCONVF_ROOK(Uplo,Way,N,A,Lda,E,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYCONVF_ROOK
   END INTERFACE
END MODULE S_CSYCONVF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYEQUB
   INTERFACE
      SUBROUTINE CSYEQUB(Uplo,N,A,Lda,S,Scond,Amax,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , ZERO = 0.0E0
      INTEGER , PARAMETER  ::  MAX_ITER = 100
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) :: Scond
      REAL , INTENT(INOUT) :: Amax
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYEQUB
   END INTERFACE
END MODULE S_CSYEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYMV
   INTERFACE
      SUBROUTINE CSYMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
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
      END SUBROUTINE CSYMV
   END INTERFACE
END MODULE S_CSYMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYR
   INTERFACE
      SUBROUTINE CSYR(Uplo,N,Alpha,X,Incx,A,Lda)
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
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE CSYR
   END INTERFACE
END MODULE S_CSYR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYRFS
   INTERFACE
      SUBROUTINE CSYRFS(Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,Ferr,&
     &                  Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYRFS
   END INTERFACE
END MODULE S_CSYRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYRFSX
   INTERFACE
      SUBROUTINE CSYRFSX(Uplo,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,S,B,Ldb,X,&
     &                   Ldx,Rcond,Berr,N_err_bnds,Err_bnds_norm,       &
     &                   Err_bnds_comp,Nparams,Params,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ITREF_DEFAULT = 1.0 ,       &
     &                      ITHRESH_DEFAULT = 10.0 ,                    &
     &                      COMPONENTWISE_DEFAULT = 1.0 ,               &
     &                      RTHRESH_DEFAULT = 0.5 ,                     &
     &                      DZTHRESH_DEFAULT = 0.25
      INTEGER , PARAMETER  ::  LA_LINRX_ITREF_I = 1 ,                   &
     &                         LA_LINRX_ITHRESH_I = 2 ,                 &
     &                         LA_LINRX_CWISE_I = 3 ,                   &
     &                         LA_LINRX_TRUST_I = 1 ,                   &
     &                         LA_LINRX_ERR_I = 2 , LA_LINRX_RCOND_I = 3
      CHARACTER :: Uplo
      CHARACTER :: Equed
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: S
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL , INTENT(INOUT) , DIMENSION(*) :: Params
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYRFSX
   END INTERFACE
END MODULE S_CSYRFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYSV_AA_2STAGE
   INTERFACE
      SUBROUTINE CSYSV_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B, &
     &                           Ldb,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tb
      INTEGER :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYSV_AA_2STAGE
   END INTERFACE
END MODULE S_CSYSV_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYSV_AA
   INTERFACE
      SUBROUTINE CSYSV_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYSV_AA
   END INTERFACE
END MODULE S_CSYSV_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYSV
   INTERFACE
      SUBROUTINE CSYSV(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYSV
   END INTERFACE
END MODULE S_CSYSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYSV_RK
   INTERFACE
      SUBROUTINE CSYSV_RK(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Work,Lwork,    &
     &                    Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYSV_RK
   END INTERFACE
END MODULE S_CSYSV_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYSV_ROOK
   INTERFACE
      SUBROUTINE CSYSV_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,    &
     &                      Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYSV_ROOK
   END INTERFACE
END MODULE S_CSYSV_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYSVX
   INTERFACE
      SUBROUTINE CSYSVX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,&
     &                  Rcond,Ferr,Berr,Work,Lwork,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYSVX
   END INTERFACE
END MODULE S_CSYSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYSVXX
   INTERFACE
      SUBROUTINE CSYSVXX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,S,B, &
     &                   Ldb,X,Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,        &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: S
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , INTENT(OUT) :: Rpvgrw
      REAL , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL , DIMENSION(*) :: Params
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYSVXX
   END INTERFACE
END MODULE S_CSYSVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYSWAPR
   INTERFACE
      SUBROUTINE CSYSWAPR(Uplo,N,A,Lda,I1,I2)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,N) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) :: I1
      INTEGER , INTENT(IN) :: I2
      END SUBROUTINE CSYSWAPR
   END INTERFACE
END MODULE S_CSYSWAPR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTF2
   INTERFACE
      SUBROUTINE CSYTF2(Uplo,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTF2
   END INTERFACE
END MODULE S_CSYTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTF2_RK
   INTERFACE
      SUBROUTINE CSYTF2_RK(Uplo,N,A,Lda,E,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTF2_RK
   END INTERFACE
END MODULE S_CSYTF2_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTF2_ROOK
   INTERFACE
      SUBROUTINE CSYTF2_ROOK(Uplo,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTF2_ROOK
   END INTERFACE
END MODULE S_CSYTF2_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRF_AA_2STAGE
   INTERFACE
      SUBROUTINE CSYTRF_AA_2STAGE(Uplo,N,A,Lda,Tb,Ltb,Ipiv,Ipiv2,Work,  &
     &                            Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRF_AA_2STAGE
   END INTERFACE
END MODULE S_CSYTRF_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRF_AA
   INTERFACE
      SUBROUTINE CSYTRF_AA(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRF_AA
   END INTERFACE
END MODULE S_CSYTRF_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRF
   INTERFACE
      SUBROUTINE CSYTRF(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRF
   END INTERFACE
END MODULE S_CSYTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRF_RK
   INTERFACE
      SUBROUTINE CSYTRF_RK(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRF_RK
   END INTERFACE
END MODULE S_CSYTRF_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRF_ROOK
   INTERFACE
      SUBROUTINE CSYTRF_ROOK(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRF_ROOK
   END INTERFACE
END MODULE S_CSYTRF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRI2
   INTERFACE
      SUBROUTINE CSYTRI2(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRI2
   END INTERFACE
END MODULE S_CSYTRI2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRI2X
   INTERFACE
      SUBROUTINE CSYTRI2X(Uplo,N,A,Lda,Ipiv,Work,Nb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRI2X
   END INTERFACE
END MODULE S_CSYTRI2X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRI_3
   INTERFACE
      SUBROUTINE CSYTRI_3(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRI_3
   END INTERFACE
END MODULE S_CSYTRI_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRI_3X
   INTERFACE
      SUBROUTINE CSYTRI_3X(Uplo,N,A,Lda,E,Ipiv,Work,Nb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRI_3X
   END INTERFACE
END MODULE S_CSYTRI_3X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRI
   INTERFACE
      SUBROUTINE CSYTRI(Uplo,N,A,Lda,Ipiv,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRI
   END INTERFACE
END MODULE S_CSYTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRI_ROOK
   INTERFACE
      SUBROUTINE CSYTRI_ROOK(Uplo,N,A,Lda,Ipiv,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRI_ROOK
   END INTERFACE
END MODULE S_CSYTRI_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRS2
   INTERFACE
      SUBROUTINE CSYTRS2(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRS2
   END INTERFACE
END MODULE S_CSYTRS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRS_3
   INTERFACE
      SUBROUTINE CSYTRS_3(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRS_3
   END INTERFACE
END MODULE S_CSYTRS_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRS_AA_2STAGE
   INTERFACE
      SUBROUTINE CSYTRS_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B,&
     &                            Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRS_AA_2STAGE
   END INTERFACE
END MODULE S_CSYTRS_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRS_AA
   INTERFACE
      SUBROUTINE CSYTRS_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRS_AA
   END INTERFACE
END MODULE S_CSYTRS_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRS
   INTERFACE
      SUBROUTINE CSYTRS(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRS
   END INTERFACE
END MODULE S_CSYTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CSYTRS_ROOK
   INTERFACE
      SUBROUTINE CSYTRS_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CSYTRS_ROOK
   END INTERFACE
END MODULE S_CSYTRS_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTBCON
   INTERFACE
      SUBROUTINE CTBCON(Norm,Uplo,Diag,N,Kd,Ab,Ldab,Rcond,Work,Rwork,   &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTBCON
   END INTERFACE
END MODULE S_CTBCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTBRFS
   INTERFACE
      SUBROUTINE CTBRFS(Uplo,Trans,Diag,N,Kd,Nrhs,Ab,Ldab,B,Ldb,X,Ldx,  &
     &                  Ferr,Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(OUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTBRFS
   END INTERFACE
END MODULE S_CTBRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTBTRS
   INTERFACE
      SUBROUTINE CTBTRS(Uplo,Trans,Diag,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTBTRS
   END INTERFACE
END MODULE S_CTBTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTFSM
   INTERFACE
      SUBROUTINE CTFSM(Transr,Side,Uplo,Trans,Diag,M,N,Alpha,A,B,Ldb)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Transr
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: M
      INTEGER :: N
      COMPLEX :: Alpha
      COMPLEX , DIMENSION(0:*) :: A
      COMPLEX , DIMENSION(0:Ldb-1,0:*) :: B
      INTEGER :: Ldb
      END SUBROUTINE CTFSM
   END INTERFACE
END MODULE S_CTFSM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTFTRI
   INTERFACE
      SUBROUTINE CTFTRI(Transr,Uplo,Diag,N,A,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Transr
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      COMPLEX , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTFTRI
   END INTERFACE
END MODULE S_CTFTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTFTTP
   INTERFACE
      SUBROUTINE CTFTTP(Transr,Uplo,N,Arf,Ap,Info)
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(0:*) :: Arf
      COMPLEX , INTENT(OUT) , DIMENSION(0:*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTFTTP
   END INTERFACE
END MODULE S_CTFTTP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTFTTR
   INTERFACE
      SUBROUTINE CTFTTR(Transr,Uplo,N,Arf,A,Lda,Info)
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(0:*) :: Arf
      COMPLEX , INTENT(OUT) , DIMENSION(0:Lda-1,0:*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTFTTR
   END INTERFACE
END MODULE S_CTFTTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTGEVC
   INTERFACE
      SUBROUTINE CTGEVC(Side,Howmny,Select,N,S,Lds,P,Ldp,Vl,Ldvl,Vr,    &
     &                  Ldvr,Mm,M,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lds,*) :: S
      INTEGER , INTENT(IN) :: Lds
      COMPLEX , INTENT(IN) , DIMENSION(Ldp,*) :: P
      INTEGER , INTENT(IN) :: Ldp
      COMPLEX , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(OUT) :: M
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTGEVC
   END INTERFACE
END MODULE S_CTGEVC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTGEX2
   INTERFACE
      SUBROUTINE CTGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,J1,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  TWENTY = 2.0E+1
      INTEGER , PARAMETER  ::  LDST = 2
      LOGICAL , PARAMETER  ::  WANDS = .TRUE.
      LOGICAL , INTENT(IN) :: Wantq
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER , INTENT(IN) :: Ldq
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      INTEGER , INTENT(IN) :: J1
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE CTGEX2
   END INTERFACE
END MODULE S_CTGEX2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTGEXC
   INTERFACE
      SUBROUTINE CTGEXC(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,Ifst,Ilst,&
     &                  Info)
      IMPLICIT NONE
      LOGICAL :: Wantq
      LOGICAL :: Wantz
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(IN) :: Ifst
      INTEGER , INTENT(INOUT) :: Ilst
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTGEXC
   END INTERFACE
END MODULE S_CTGEXC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTGSEN
   INTERFACE
      SUBROUTINE CTGSEN(Ijob,Wantq,Wantz,Select,N,A,Lda,B,Ldb,Alpha,    &
     &                  Beta,Q,Ldq,Z,Ldz,M,Pl,Pr,Dif,Work,Lwork,Iwork,  &
     &                  Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  IDIFJB = 3
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , INTENT(IN) :: Ijob
      LOGICAL :: Wantq
      LOGICAL :: Wantz
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: Alpha
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: Beta
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) :: Pl
      REAL , INTENT(INOUT) :: Pr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dif
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTGSEN
   END INTERFACE
END MODULE S_CTGSEN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTGSJA
   INTERFACE
      SUBROUTINE CTGSJA(Jobu,Jobv,Jobq,M,P,N,K,L,A,Lda,B,Ldb,Tola,Tolb, &
     &                  Alpha,Beta,U,Ldu,V,Ldv,Q,Ldq,Work,Ncycle,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXIT = 40
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      CHARACTER :: Jobq
      INTEGER :: M
      INTEGER :: P
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER :: L
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(IN) :: Tola
      REAL , INTENT(IN) :: Tolb
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alpha
      REAL , INTENT(INOUT) , DIMENSION(*) :: Beta
      COMPLEX , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) :: Ncycle
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTGSJA
   END INTERFACE
END MODULE S_CTGSJA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTGSNA
   INTERFACE
      SUBROUTINE CTGSNA(Job,Howmny,Select,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,&
     &                  S,Dif,Mm,M,Work,Lwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , PARAMETER  ::  IDIFJB = 3
      CHARACTER :: Job
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      COMPLEX , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      REAL , INTENT(OUT) , DIMENSION(*) :: S
      REAL , DIMENSION(*) :: Dif
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTGSNA
   END INTERFACE
END MODULE S_CTGSNA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTGSY2
   INTERFACE
      SUBROUTINE CTGSY2(Trans,Ijob,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F, &
     &                  Ldf,Scale,Rdsum,Rdscal,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , PARAMETER  ::  LDZ = 2
      CHARACTER :: Trans
      INTEGER :: Ijob
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(Ldd,*) :: D
      INTEGER , INTENT(IN) :: Ldd
      COMPLEX , DIMENSION(Lde,*) :: E
      INTEGER :: Lde
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      REAL , INTENT(INOUT) :: Scale
      REAL :: Rdsum
      REAL :: Rdscal
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTGSY2
   END INTERFACE
END MODULE S_CTGSY2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTGSYL
   INTERFACE
      SUBROUTINE CTGSYL(Trans,Ijob,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F, &
     &                  Ldf,Scale,Dif,Work,Lwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: Ijob
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(Ldd,*) :: D
      INTEGER :: Ldd
      COMPLEX , DIMENSION(Lde,*) :: E
      INTEGER :: Lde
      COMPLEX , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(OUT) :: Dif
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTGSYL
   END INTERFACE
END MODULE S_CTGSYL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPCON
   INTERFACE
      SUBROUTINE CTPCON(Norm,Uplo,Diag,N,Ap,Rcond,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: Ap
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTPCON
   END INTERFACE
END MODULE S_CTPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPLQT2
   INTERFACE
      SUBROUTINE CTPLQT2(M,N,L,A,Lda,B,Ldb,T,Ldt,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: L
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTPLQT2
   END INTERFACE
END MODULE S_CTPLQT2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPLQT
   INTERFACE
      SUBROUTINE CTPLQT(M,N,L,Mb,A,Lda,B,Ldb,T,Ldt,Work,Info)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Mb
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTPLQT
   END INTERFACE
END MODULE S_CTPLQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPMLQT
   INTERFACE
      SUBROUTINE CTPMLQT(Side,Trans,M,N,K,L,Mb,V,Ldv,T,Ldt,A,Lda,B,Ldb, &
     &                   Work,Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Mb
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTPMLQT
   END INTERFACE
END MODULE S_CTPMLQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPMQRT
   INTERFACE
      SUBROUTINE CTPMQRT(Side,Trans,M,N,K,L,Nb,V,Ldv,T,Ldt,A,Lda,B,Ldb, &
     &                   Work,Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTPMQRT
   END INTERFACE
END MODULE S_CTPMQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPQRT2
   INTERFACE
      SUBROUTINE CTPQRT2(M,N,L,A,Lda,B,Ldb,T,Ldt,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0,0.0) , ZERO = (0.0,0.0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: L
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTPQRT2
   END INTERFACE
END MODULE S_CTPQRT2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPQRT
   INTERFACE
      SUBROUTINE CTPQRT(M,N,L,Nb,A,Lda,B,Ldb,T,Ldt,Work,Info)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTPQRT
   END INTERFACE
END MODULE S_CTPQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPRFB
   INTERFACE
      SUBROUTINE CTPRFB(Side,Trans,Direct,Storev,M,N,K,L,V,Ldv,T,Ldt,A, &
     &                  Lda,B,Ldb,Work,Ldwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0,0.0) , ZERO = (0.0,0.0)
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: L
      COMPLEX , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE CTPRFB
   END INTERFACE
END MODULE S_CTPRFB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPRFS
   INTERFACE
      SUBROUTINE CTPRFS(Uplo,Trans,Diag,N,Nrhs,Ap,B,Ldb,X,Ldx,Ferr,Berr,&
     &                  Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(OUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTPRFS
   END INTERFACE
END MODULE S_CTPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPTRI
   INTERFACE
      SUBROUTINE CTPTRI(Uplo,Diag,N,Ap,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTPTRI
   END INTERFACE
END MODULE S_CTPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPTRS
   INTERFACE
      SUBROUTINE CTPTRS(Uplo,Trans,Diag,N,Nrhs,Ap,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTPTRS
   END INTERFACE
END MODULE S_CTPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPTTF
   INTERFACE
      SUBROUTINE CTPTTF(Transr,Uplo,N,Ap,Arf,Info)
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(0:*) :: Ap
      COMPLEX , INTENT(OUT) , DIMENSION(0:*) :: Arf
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTPTTF
   END INTERFACE
END MODULE S_CTPTTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTPTTR
   INTERFACE
      SUBROUTINE CTPTTR(Uplo,N,Ap,A,Lda,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX , INTENT(OUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTPTTR
   END INTERFACE
END MODULE S_CTPTTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRCON
   INTERFACE
      SUBROUTINE CTRCON(Norm,Uplo,Diag,N,A,Lda,Rcond,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTRCON
   END INTERFACE
END MODULE S_CTRCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTREVC3
   INTERFACE
      SUBROUTINE CTREVC3(Side,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,Mm, &
     &                   M,Work,Lwork,Rwork,Lrwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      INTEGER , PARAMETER  ::  NBMIN = 8 , NBMAX = 128
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTREVC3
   END INTERFACE
END MODULE S_CTREVC3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTREVC
   INTERFACE
      SUBROUTINE CTREVC(Side,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,Mm,M,&
     &                  Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CMZERO = (0.0E+0,0.0E+0) ,               &
     &                         CMONE = (1.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTREVC
   END INTERFACE
END MODULE S_CTREVC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTREXC
   INTERFACE
      SUBROUTINE CTREXC(Compq,N,T,Ldt,Q,Ldq,Ifst,Ilst,Info)
      IMPLICIT NONE
      CHARACTER :: Compq
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER , INTENT(IN) :: Ldq
      INTEGER , INTENT(IN) :: Ifst
      INTEGER , INTENT(IN) :: Ilst
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTREXC
   END INTERFACE
END MODULE S_CTREXC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRRFS
   INTERFACE
      SUBROUTINE CTRRFS(Uplo,Trans,Diag,N,Nrhs,A,Lda,B,Ldb,X,Ldx,Ferr,  &
     &                  Berr,Work,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(OUT) , DIMENSION(*) :: Berr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTRRFS
   END INTERFACE
END MODULE S_CTRRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRSEN
   INTERFACE
      SUBROUTINE CTRSEN(Job,Compq,Select,N,T,Ldt,Q,Ldq,W,M,S,Sep,Work,  &
     &                  Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Job
      CHARACTER :: Compq
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: W
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(OUT) :: S
      REAL , INTENT(OUT) :: Sep
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTRSEN
   END INTERFACE
END MODULE S_CTRSEN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRSNA
   INTERFACE
      SUBROUTINE CTRSNA(Job,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,S,Sep,&
     &                  Mm,M,Work,Ldwork,Rwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0 + 0
      CHARACTER :: Job
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      COMPLEX , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      REAL , INTENT(OUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) , DIMENSION(*) :: Sep
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTRSNA
   END INTERFACE
END MODULE S_CTRSNA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRSYL
   INTERFACE
      SUBROUTINE CTRSYL(Trana,Tranb,Isgn,M,N,A,Lda,B,Ldb,C,Ldc,Scale,   &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Trana
      CHARACTER :: Tranb
      INTEGER , INTENT(IN) :: Isgn
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , INTENT(INOUT) :: Scale
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTRSYL
   END INTERFACE
END MODULE S_CTRSYL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRTI2
   INTERFACE
      SUBROUTINE CTRTI2(Uplo,Diag,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTRTI2
   END INTERFACE
END MODULE S_CTRTI2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRTRI
   INTERFACE
      SUBROUTINE CTRTRI(Uplo,Diag,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTRTRI
   END INTERFACE
END MODULE S_CTRTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRTRS
   INTERFACE
      SUBROUTINE CTRTRS(Uplo,Trans,Diag,N,Nrhs,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTRTRS
   END INTERFACE
END MODULE S_CTRTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRTTF
   INTERFACE
      SUBROUTINE CTRTTF(Transr,Uplo,N,A,Lda,Arf,Info)
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(0:Lda-1,0:*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(OUT) , DIMENSION(0:*) :: Arf
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTRTTF
   END INTERFACE
END MODULE S_CTRTTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTRTTP
   INTERFACE
      SUBROUTINE CTRTTP(Uplo,N,A,Lda,Ap,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTRTTP
   END INTERFACE
END MODULE S_CTRTTP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CTZRZF
   INTERFACE
      SUBROUTINE CTZRZF(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CTZRZF
   END INTERFACE
END MODULE S_CTZRZF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNBDB1
   INTERFACE
      SUBROUTINE CUNBDB1(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E0,0.0E0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL , INTENT(OUT) , DIMENSION(*) :: Phi
      COMPLEX , DIMENSION(*) :: Taup1
      COMPLEX , DIMENSION(*) :: Taup2
      COMPLEX , DIMENSION(*) :: Tauq1
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNBDB1
   END INTERFACE
END MODULE S_CUNBDB1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNBDB2
   INTERFACE
      SUBROUTINE CUNBDB2(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  NEGONE = (-1.0E0,0.0E0) ,                &
     &                         ONE = (1.0E0,0.0E0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , INTENT(OUT) , DIMENSION(*) :: Theta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Phi
      COMPLEX , DIMENSION(*) :: Taup1
      COMPLEX , DIMENSION(*) :: Taup2
      COMPLEX , DIMENSION(*) :: Tauq1
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNBDB2
   END INTERFACE
END MODULE S_CUNBDB2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNBDB3
   INTERFACE
      SUBROUTINE CUNBDB3(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E0,0.0E0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , INTENT(OUT) , DIMENSION(*) :: Theta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Phi
      COMPLEX , DIMENSION(*) :: Taup1
      COMPLEX , DIMENSION(*) :: Taup2
      COMPLEX , DIMENSION(*) :: Tauq1
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNBDB3
   END INTERFACE
END MODULE S_CUNBDB3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNBDB4
   INTERFACE
      SUBROUTINE CUNBDB4(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Phantom,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  NEGONE = (-1.0E0,0.0E0) ,                &
     &                         ONE = (1.0E0,0.0E0) ,                    &
     &                         ZERO = (0.0E0,0.0E0)
      INTEGER , INTENT(IN) :: M
      INTEGER :: P
      INTEGER :: Q
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL , INTENT(OUT) , DIMENSION(*) :: Phi
      COMPLEX , DIMENSION(*) :: Taup1
      COMPLEX , DIMENSION(*) :: Taup2
      COMPLEX , DIMENSION(*) :: Tauq1
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Phantom
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNBDB4
   END INTERFACE
END MODULE S_CUNBDB4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNBDB5
   INTERFACE
      SUBROUTINE CUNBDB5(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,&
     &                   Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E0,0.0E0) ,                    &
     &                         ZERO = (0.0E0,0.0E0)
      INTEGER :: M1
      INTEGER :: M2
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: X1
      INTEGER :: Incx1
      COMPLEX , DIMENSION(*) :: X2
      INTEGER :: Incx2
      COMPLEX , DIMENSION(Ldq1,*) :: Q1
      INTEGER :: Ldq1
      COMPLEX , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNBDB5
   END INTERFACE
END MODULE S_CUNBDB5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNBDB6
   INTERFACE
      SUBROUTINE CUNBDB6(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,&
     &                   Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ALPHASQ = 0.01E0 , REALONE = 1.0E0 ,        &
     &                      REALZERO = 0.0E0
      COMPLEX , PARAMETER  ::  NEGONE = (-1.0E0,0.0E0) ,                &
     &                         ONE = (1.0E0,0.0E0) ,                    &
     &                         ZERO = (0.0E0,0.0E0)
      INTEGER :: M1
      INTEGER :: M2
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: X1
      INTEGER :: Incx1
      COMPLEX , DIMENSION(*) :: X2
      INTEGER :: Incx2
      COMPLEX , DIMENSION(Ldq1,*) :: Q1
      INTEGER :: Ldq1
      COMPLEX , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNBDB6
   END INTERFACE
END MODULE S_CUNBDB6
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNBDB
   INTERFACE
      SUBROUTINE CUNBDB(Trans,Signs,M,P,Q,X11,Ldx11,X12,Ldx12,X21,Ldx21,&
     &                  X22,Ldx22,Theta,Phi,Taup1,Taup2,Tauq1,Tauq2,    &
     &                  Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  REALONE = 1.0E0
      COMPLEX , PARAMETER  ::  ONE = (1.0E0,0.0E0)
      CHARACTER :: Trans
      CHARACTER :: Signs
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      COMPLEX , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX , DIMENSION(Ldx12,*) :: X12
      INTEGER :: Ldx12
      COMPLEX , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      COMPLEX , DIMENSION(Ldx22,*) :: X22
      INTEGER :: Ldx22
      REAL , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Phi
      COMPLEX , DIMENSION(*) :: Taup1
      COMPLEX , DIMENSION(*) :: Taup2
      COMPLEX , DIMENSION(*) :: Tauq1
      COMPLEX , DIMENSION(*) :: Tauq2
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNBDB
   END INTERFACE
END MODULE S_CUNBDB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNCSD2BY1
   INTERFACE
      SUBROUTINE CUNCSD2BY1(Jobu1,Jobu2,Jobv1t,M,P,Q,X11,Ldx11,X21,     &
     &                      Ldx21,Theta,U1,Ldu1,U2,Ldu2,V1t,Ldv1t,Work, &
     &                      Lwork,Rwork,Lrwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E0,0.0E0) ,                    &
     &                         ZERO = (0.0E0,0.0E0)
      CHARACTER :: Jobu1
      CHARACTER :: Jobu2
      CHARACTER :: Jobv1t
      INTEGER :: M
      INTEGER :: P
      INTEGER :: Q
      COMPLEX , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , DIMENSION(*) :: Theta
      COMPLEX , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      COMPLEX , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      COMPLEX , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNCSD2BY1
   END INTERFACE
END MODULE S_CUNCSD2BY1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNCSD
   INTERFACE
      RECURSIVE SUBROUTINE CUNCSD(Jobu1,Jobu2,Jobv1t,Jobv2t,Trans,Signs,&
     &                            M,P,Q,X11,Ldx11,X12,Ldx12,X21,Ldx21,  &
     &                            X22,Ldx22,Theta,U1,Ldu1,U2,Ldu2,V1t,  &
     &                            Ldv1t,V2t,Ldv2t,Work,Lwork,Rwork,     &
     &                            Lrwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E0,0.0E0) ,                    &
     &                         ZERO = (0.0E0,0.0E0)
      CHARACTER :: Jobu1
      CHARACTER :: Jobu2
      CHARACTER :: Jobv1t
      CHARACTER :: Jobv2t
      CHARACTER :: Trans
      CHARACTER :: Signs
      INTEGER :: M
      INTEGER :: P
      INTEGER :: Q
      COMPLEX , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX , DIMENSION(Ldx12,*) :: X12
      INTEGER :: Ldx12
      COMPLEX , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      COMPLEX , DIMENSION(Ldx22,*) :: X22
      INTEGER :: Ldx22
      REAL , DIMENSION(*) :: Theta
      COMPLEX , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      COMPLEX , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      COMPLEX , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      COMPLEX , DIMENSION(Ldv2t,*) :: V2t
      INTEGER :: Ldv2t
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNCSD
   END INTERFACE
END MODULE S_CUNCSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNG2L
   INTERFACE
      SUBROUTINE CUNG2L(M,N,K,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNG2L
   END INTERFACE
END MODULE S_CUNG2L
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNG2R
   INTERFACE
      SUBROUTINE CUNG2R(M,N,K,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNG2R
   END INTERFACE
END MODULE S_CUNG2R
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNGBR
   INTERFACE
      SUBROUTINE CUNGBR(Vect,M,N,K,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Vect
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNGBR
   END INTERFACE
END MODULE S_CUNGBR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNGHR
   INTERFACE
      SUBROUTINE CUNGHR(N,Ilo,Ihi,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNGHR
   END INTERFACE
END MODULE S_CUNGHR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNGL2
   INTERFACE
      SUBROUTINE CUNGL2(M,N,K,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNGL2
   END INTERFACE
END MODULE S_CUNGL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNGLQ
   INTERFACE
      SUBROUTINE CUNGLQ(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNGLQ
   END INTERFACE
END MODULE S_CUNGLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNGQL
   INTERFACE
      SUBROUTINE CUNGQL(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNGQL
   END INTERFACE
END MODULE S_CUNGQL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNGQR
   INTERFACE
      SUBROUTINE CUNGQR(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNGQR
   END INTERFACE
END MODULE S_CUNGQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNGR2
   INTERFACE
      SUBROUTINE CUNGR2(M,N,K,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNGR2
   END INTERFACE
END MODULE S_CUNGR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNGRQ
   INTERFACE
      SUBROUTINE CUNGRQ(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNGRQ
   END INTERFACE
END MODULE S_CUNGRQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNGTR
   INTERFACE
      SUBROUTINE CUNGTR(Uplo,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNGTR
   END INTERFACE
END MODULE S_CUNGTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNGTSQR
   INTERFACE
      SUBROUTINE CUNGTSQR(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNGTSQR
   END INTERFACE
END MODULE S_CUNGTSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNGTSQR_ROW
   INTERFACE
      SUBROUTINE CUNGTSQR_ROW(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Mb
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNGTSQR_ROW
   END INTERFACE
END MODULE S_CUNGTSQR_ROW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNHR_COL
   INTERFACE
      SUBROUTINE CUNHR_COL(M,N,Nb,A,Lda,T,Ldt,D,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0) ,                 &
     &                         CZERO = (0.0E+0,0.0E+0)
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nb
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNHR_COL
   END INTERFACE
END MODULE S_CUNHR_COL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNM22
   INTERFACE
      SUBROUTINE CUNM22(Side,Trans,M,N,N1,N2,Q,Ldq,C,Ldc,Work,Lwork,    &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: N1
      INTEGER :: N2
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNM22
   END INTERFACE
END MODULE S_CUNM22
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNM2L
   INTERFACE
      SUBROUTINE CUNM2L(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNM2L
   END INTERFACE
END MODULE S_CUNM2L
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNM2R
   INTERFACE
      SUBROUTINE CUNM2R(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNM2R
   END INTERFACE
END MODULE S_CUNM2R
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNMBR
   INTERFACE
      SUBROUTINE CUNMBR(Vect,Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,     &
     &                  Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Vect
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNMBR
   END INTERFACE
END MODULE S_CUNMBR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNMHR
   INTERFACE
      SUBROUTINE CUNMHR(Side,Trans,M,N,Ilo,Ihi,A,Lda,Tau,C,Ldc,Work,    &
     &                  Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNMHR
   END INTERFACE
END MODULE S_CUNMHR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNML2
   INTERFACE
      SUBROUTINE CUNML2(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNML2
   END INTERFACE
END MODULE S_CUNML2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNMLQ
   INTERFACE
      SUBROUTINE CUNMLQ(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDT = NBMAX + 1 ,           &
     &                         TSIZE = LDT*NBMAX
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNMLQ
   END INTERFACE
END MODULE S_CUNMLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNMQL
   INTERFACE
      SUBROUTINE CUNMQL(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDT = NBMAX + 1 ,           &
     &                         TSIZE = LDT*NBMAX
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNMQL
   END INTERFACE
END MODULE S_CUNMQL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNMQR
   INTERFACE
      SUBROUTINE CUNMQR(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDT = NBMAX + 1 ,           &
     &                         TSIZE = LDT*NBMAX
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNMQR
   END INTERFACE
END MODULE S_CUNMQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNMR2
   INTERFACE
      SUBROUTINE CUNMR2(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNMR2
   END INTERFACE
END MODULE S_CUNMR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNMR3
   INTERFACE
      SUBROUTINE CUNMR3(Side,Trans,M,N,K,L,A,Lda,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      INTEGER :: L
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNMR3
   END INTERFACE
END MODULE S_CUNMR3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNMRQ
   INTERFACE
      SUBROUTINE CUNMRQ(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDT = NBMAX + 1 ,           &
     &                         TSIZE = LDT*NBMAX
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNMRQ
   END INTERFACE
END MODULE S_CUNMRQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNMRZ
   INTERFACE
      SUBROUTINE CUNMRZ(Side,Trans,M,N,K,L,A,Lda,Tau,C,Ldc,Work,Lwork,  &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDT = NBMAX + 1 ,           &
     &                         TSIZE = LDT*NBMAX
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: L
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNMRZ
   END INTERFACE
END MODULE S_CUNMRZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUNMTR
   INTERFACE
      SUBROUTINE CUNMTR(Side,Uplo,Trans,M,N,A,Lda,Tau,C,Ldc,Work,Lwork, &
     &                  Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUNMTR
   END INTERFACE
END MODULE S_CUNMTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUPGTR
   INTERFACE
      SUBROUTINE CUPGTR(Uplo,N,Ap,Tau,Q,Ldq,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUPGTR
   END INTERFACE
END MODULE S_CUPGTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_CUPMTR
   INTERFACE
      SUBROUTINE CUPMTR(Side,Uplo,Trans,M,N,Ap,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Ap
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE CUPMTR
   END INTERFACE
END MODULE S_CUPMTR
