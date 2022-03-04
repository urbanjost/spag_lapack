!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SBBCSD
   INTERFACE
      SUBROUTINE SBBCSD(Jobu1,Jobu2,Jobv1t,Jobv2t,Trans,M,P,Q,Theta,Phi,&
     &                  U1,Ldu1,U2,Ldu2,V1t,Ldv1t,V2t,Ldv2t,B11d,B11e,  &
     &                  B12d,B12e,B21d,B21e,B22d,B22e,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXITR = 6
      REAL , PARAMETER  ::  HUNDRED = 100.0E0 , MEIGHTH = -0.125E0 ,    &
     &                      ONE = 1.0E0 , TEN = 10.0E0 , ZERO = 0.0E0 , &
     &                      NEGONE = -1.0E0 , PIOVER2 =                 &
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
      REAL , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      REAL , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      REAL , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      REAL , DIMENSION(Ldv2t,*) :: V2t
      INTEGER :: Ldv2t
      REAL , INTENT(INOUT) , DIMENSION(*) :: B11d
      REAL , INTENT(INOUT) , DIMENSION(*) :: B11e
      REAL , INTENT(INOUT) , DIMENSION(*) :: B12d
      REAL , INTENT(INOUT) , DIMENSION(*) :: B12e
      REAL , INTENT(INOUT) , DIMENSION(*) :: B21d
      REAL , INTENT(INOUT) , DIMENSION(*) :: B21e
      REAL , INTENT(INOUT) , DIMENSION(*) :: B22d
      REAL , INTENT(INOUT) , DIMENSION(*) :: B22e
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SBBCSD
   END INTERFACE
END MODULE S_SBBCSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SBDSDC
   INTERFACE
      SUBROUTINE SBDSDC(Uplo,Compq,N,D,E,U,Ldu,Vt,Ldvt,Q,Iq,Work,Iwork, &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Compq
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL , DIMENSION(*) :: Q
      INTEGER , DIMENSION(*) :: Iq
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SBDSDC
   END INTERFACE
END MODULE S_SBDSDC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SBDSQR
   INTERFACE
      SUBROUTINE SBDSQR(Uplo,N,Ncvt,Nru,Ncc,D,E,Vt,Ldvt,U,Ldu,C,Ldc,    &
     &                  Work,Info)
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
      REAL , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SBDSQR
   END INTERFACE
END MODULE S_SBDSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SBDSVDX
   INTERFACE
      SUBROUTINE SBDSVDX(Uplo,Jobz,Range,N,D,E,Vl,Vu,Il,Iu,Ns,S,Z,Ldz,  &
     &                   Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TEN = 10.0E0 , &
     &                      HNDRD = 100.0E0 , MEIGTH = -0.1250E0 ,      &
     &                      FUDGE = 2.0E0
      CHARACTER :: Uplo
      CHARACTER :: Jobz
      CHARACTER :: Range
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      INTEGER , INTENT(INOUT) :: Ns
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SBDSVDX
   END INTERFACE
END MODULE S_SBDSVDX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SCOMBSSQ
   INTERFACE
      SUBROUTINE SCOMBSSQ(V1,V2)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0D+0
      REAL , INTENT(INOUT) , DIMENSION(2) :: V1
      REAL , INTENT(IN) , DIMENSION(2) :: V2
      END SUBROUTINE SCOMBSSQ
   END INTERFACE
END MODULE S_SCOMBSSQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SCSUM1
   INTERFACE
      FUNCTION SCSUM1(N,Cx,Incx)
      IMPLICIT NONE
      REAL :: SCSUM1
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION SCSUM1
   END INTERFACE
END MODULE S_SCSUM1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SDISNA
   INTERFACE
      SUBROUTINE SDISNA(Job,M,N,D,Sep,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Job
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sep
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SDISNA
   END INTERFACE
END MODULE S_SDISNA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBBRD
   INTERFACE
      SUBROUTINE SGBBRD(Vect,M,N,Ncc,Kl,Ku,Ab,Ldab,D,E,Q,Ldq,Pt,Ldpt,C, &
     &                  Ldc,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Vect
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Ncc
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(OUT) , DIMENSION(*) :: D
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , DIMENSION(Ldpt,*) :: Pt
      INTEGER :: Ldpt
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGBBRD
   END INTERFACE
END MODULE S_SGBBRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBCON
   INTERFACE
      SUBROUTINE SGBCON(Norm,N,Kl,Ku,Ab,Ldab,Ipiv,Anorm,Rcond,Work,     &
     &                  Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Norm
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGBCON
   END INTERFACE
END MODULE S_SGBCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBEQUB
   INTERFACE
      SUBROUTINE SGBEQUB(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: R
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      REAL , INTENT(OUT) :: Rowcnd
      REAL , INTENT(OUT) :: Colcnd
      REAL , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGBEQUB
   END INTERFACE
END MODULE S_SGBEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBEQU
   INTERFACE
      SUBROUTINE SGBEQU(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: R
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      REAL , INTENT(OUT) :: Rowcnd
      REAL , INTENT(OUT) :: Colcnd
      REAL , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGBEQU
   END INTERFACE
END MODULE S_SGBEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBRFS
   INTERFACE
      SUBROUTINE SGBRFS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv,B,Ldb,&
     &                  X,Ldx,Ferr,Berr,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER , INTENT(IN) :: Nrhs
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGBRFS
   END INTERFACE
END MODULE S_SGBRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBRFSX
   INTERFACE
      SUBROUTINE SGBRFSX(Trans,Equed,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,    &
     &                   Ipiv,R,C,B,Ldb,X,Ldx,Rcond,Berr,N_err_bnds,    &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ITREF_DEFAULT = 1.0 ,                       &
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
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: R
      REAL , DIMENSION(*) :: C
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL , INTENT(INOUT) , DIMENSION(*) :: Params
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGBRFSX
   END INTERFACE
END MODULE S_SGBRFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBSV
   INTERFACE
      SUBROUTINE SGBSV(N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER :: Nrhs
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGBSV
   END INTERFACE
END MODULE S_SGBSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBSVX
   INTERFACE
      SUBROUTINE SGBSVX(Fact,Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv, &
     &                  Equed,R,C,B,Ldb,X,Ldx,Rcond,Ferr,Berr,Work,     &
     &                  Iwork,Info)
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
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: R
      REAL , DIMENSION(*) :: C
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGBSVX
   END INTERFACE
END MODULE S_SGBSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBSVXX
   INTERFACE
      SUBROUTINE SGBSVXX(Fact,Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv,&
     &                   Equed,R,C,B,Ldb,X,Ldx,Rcond,Rpvgrw,Berr,       &
     &                   N_err_bnds,Err_bnds_norm,Err_bnds_comp,Nparams,&
     &                   Params,Work,Iwork,Info)
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
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL , INTENT(INOUT) , DIMENSION(*) :: R
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , INTENT(OUT) :: Rpvgrw
      REAL , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL , DIMENSION(*) :: Params
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGBSVXX
   END INTERFACE
END MODULE S_SGBSVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBTF2
   INTERFACE
      SUBROUTINE SGBTF2(M,N,Kl,Ku,Ab,Ldab,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGBTF2
   END INTERFACE
END MODULE S_SGBTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBTRF
   INTERFACE
      SUBROUTINE SGBTRF(M,N,Kl,Ku,Ab,Ldab,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDWORK = NBMAX + 1
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      REAL , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGBTRF
   END INTERFACE
END MODULE S_SGBTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGBTRS
   INTERFACE
      SUBROUTINE SGBTRS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      INTEGER :: Nrhs
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGBTRS
   END INTERFACE
END MODULE S_SGBTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEBAK
   INTERFACE
      SUBROUTINE SGEBAK(Job,Side,N,Ilo,Ihi,Scale,M,V,Ldv,Info)
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
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEBAK
   END INTERFACE
END MODULE S_SGEBAK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEBAL
   INTERFACE
      SUBROUTINE SGEBAL(Job,N,A,Lda,Ilo,Ihi,Scale,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      SCLFAC = 2.0E+0 , FACTOR = 0.95E+0
      CHARACTER :: Job
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) :: Ilo
      INTEGER , INTENT(OUT) :: Ihi
      REAL , INTENT(INOUT) , DIMENSION(*) :: Scale
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEBAL
   END INTERFACE
END MODULE S_SGEBAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEBD2
   INTERFACE
      SUBROUTINE SGEBD2(M,N,A,Lda,D,E,Tauq,Taup,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Tauq
      REAL , DIMENSION(*) :: Taup
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEBD2
   END INTERFACE
END MODULE S_SGEBD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEBRD
   INTERFACE
      SUBROUTINE SGEBRD(M,N,A,Lda,D,E,Tauq,Taup,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Tauq
      REAL , DIMENSION(*) :: Taup
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEBRD
   END INTERFACE
END MODULE S_SGEBRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGECON
   INTERFACE
      SUBROUTINE SGECON(Norm,N,A,Lda,Anorm,Rcond,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Norm
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGECON
   END INTERFACE
END MODULE S_SGECON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEEQUB
   INTERFACE
      SUBROUTINE SGEEQUB(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: R
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      REAL , INTENT(OUT) :: Rowcnd
      REAL , INTENT(OUT) :: Colcnd
      REAL , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEEQUB
   END INTERFACE
END MODULE S_SGEEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEEQU
   INTERFACE
      SUBROUTINE SGEEQU(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: R
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      REAL , INTENT(OUT) :: Rowcnd
      REAL , INTENT(OUT) :: Colcnd
      REAL , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEEQU
   END INTERFACE
END MODULE S_SGEEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEES
   INTERFACE
      SUBROUTINE SGEES(Jobvs,Sort,SELECT,N,A,Lda,Sdim,Wr,Wi,Vs,Ldvs,    &
     &                 Work,Lwork,Bwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobvs
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELECT
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Sdim
      REAL , DIMENSION(*) :: Wr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Wi
      REAL , DIMENSION(Ldvs,*) :: Vs
      INTEGER :: Ldvs
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEES
   END INTERFACE
END MODULE S_SGEES
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEESX
   INTERFACE
      SUBROUTINE SGEESX(Jobvs,Sort,SELECT,Sense,N,A,Lda,Sdim,Wr,Wi,Vs,  &
     &                  Ldvs,Rconde,Rcondv,Work,Lwork,Iwork,Liwork,     &
     &                  Bwork,Info)
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
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Sdim
      REAL , DIMENSION(*) :: Wr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Wi
      REAL , DIMENSION(Ldvs,*) :: Vs
      INTEGER :: Ldvs
      REAL :: Rconde
      REAL , INTENT(INOUT) :: Rcondv
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEESX
   END INTERFACE
END MODULE S_SGEESX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEEV
   INTERFACE
      SUBROUTINE SGEEV(Jobvl,Jobvr,N,A,Lda,Wr,Wi,Vl,Ldvl,Vr,Ldvr,Work,  &
     &                 Lwork,Info)
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
      REAL , DIMENSION(*) :: Wr
      REAL , DIMENSION(*) :: Wi
      REAL , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEEV
   END INTERFACE
END MODULE S_SGEEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEEVX
   INTERFACE
      SUBROUTINE SGEEVX(Balanc,Jobvl,Jobvr,Sense,N,A,Lda,Wr,Wi,Vl,Ldvl, &
     &                  Vr,Ldvr,Ilo,Ihi,Scale,Abnrm,Rconde,Rcondv,Work, &
     &                  Lwork,Iwork,Info)
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Wr
      REAL , DIMENSION(*) :: Wi
      REAL , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL , DIMENSION(*) :: Scale
      REAL , INTENT(INOUT) :: Abnrm
      REAL , DIMENSION(*) :: Rconde
      REAL , DIMENSION(*) :: Rcondv
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEEVX
   END INTERFACE
END MODULE S_SGEEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEHD2
   INTERFACE
      SUBROUTINE SGEHD2(N,Ilo,Ihi,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER :: Ihi
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEHD2
   END INTERFACE
END MODULE S_SGEHD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEHRD
   INTERFACE
      SUBROUTINE SGEHRD(N,Ilo,Ihi,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDT = NBMAX + 1 ,           &
     &                         TSIZE = LDT*NBMAX
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEHRD
   END INTERFACE
END MODULE S_SGEHRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEJSV
   INTERFACE
      SUBROUTINE SGEJSV(Joba,Jobu,Jobv,Jobr,Jobt,Jobp,M,N,A,Lda,Sva,U,  &
     &                  Ldu,V,Ldv,Work,Lwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER(1) :: Joba
      CHARACTER(1) :: Jobu
      CHARACTER(1) :: Jobv
      CHARACTER(1) :: Jobr
      CHARACTER(1) :: Jobt
      CHARACTER(1) :: Jobp
      INTEGER :: M
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(N) :: Sva
      REAL , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , INTENT(INOUT) , DIMENSION(Lwork) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEJSV
   END INTERFACE
END MODULE S_SGEJSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGELQ2
   INTERFACE
      SUBROUTINE SGELQ2(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGELQ2
   END INTERFACE
END MODULE S_SGELQ2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGELQ
   INTERFACE
      SUBROUTINE SGELQ(M,N,A,Lda,T,Tsize,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGELQ
   END INTERFACE
END MODULE S_SGELQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGELQF
   INTERFACE
      SUBROUTINE SGELQF(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGELQF
   END INTERFACE
END MODULE S_SGELQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGELQT3
   INTERFACE
      RECURSIVE SUBROUTINE SGELQT3(M,N,A,Lda,T,Ldt,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+00
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGELQT3
   END INTERFACE
END MODULE S_SGELQT3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGELQT
   INTERFACE
      SUBROUTINE SGELQT(M,N,Mb,A,Lda,T,Ldt,Work,Info)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Mb
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGELQT
   END INTERFACE
END MODULE S_SGELQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGELSD
   INTERFACE
      SUBROUTINE SGELSD(M,N,Nrhs,A,Lda,B,Ldb,S,Rcond,Rank,Work,Lwork,   &
     &                  Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: S
      REAL :: Rcond
      INTEGER :: Rank
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGELSD
   END INTERFACE
END MODULE S_SGELSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGELS
   INTERFACE
      SUBROUTINE SGELS(Trans,M,N,Nrhs,A,Lda,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGELS
   END INTERFACE
END MODULE S_SGELS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGELSS
   INTERFACE
      SUBROUTINE SGELSS(M,N,Nrhs,A,Lda,B,Ldb,S,Rcond,Rank,Work,Lwork,   &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGELSS
   END INTERFACE
END MODULE S_SGELSS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGELSY
   INTERFACE
      SUBROUTINE SGELSY(M,N,Nrhs,A,Lda,B,Ldb,Jpvt,Rcond,Rank,Work,Lwork,&
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  IMAX = 1 , IMIN = 2
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
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
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGELSY
   END INTERFACE
END MODULE S_SGELSY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEMLQ
   INTERFACE
      SUBROUTINE SGEMLQ(Side,Trans,M,N,K,A,Lda,T,Tsize,C,Ldc,Work,Lwork,&
     &                  Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEMLQ
   END INTERFACE
END MODULE S_SGEMLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEMLQT
   INTERFACE
      SUBROUTINE SGEMLQT(Side,Trans,M,N,K,Mb,V,Ldv,T,Ldt,C,Ldc,Work,    &
     &                   Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Mb
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEMLQT
   END INTERFACE
END MODULE S_SGEMLQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEMQR
   INTERFACE
      SUBROUTINE SGEMQR(Side,Trans,M,N,K,A,Lda,T,Tsize,C,Ldc,Work,Lwork,&
     &                  Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEMQR
   END INTERFACE
END MODULE S_SGEMQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEMQRT
   INTERFACE
      SUBROUTINE SGEMQRT(Side,Trans,M,N,K,Nb,V,Ldv,T,Ldt,C,Ldc,Work,    &
     &                   Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Nb
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEMQRT
   END INTERFACE
END MODULE S_SGEMQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEQL2
   INTERFACE
      SUBROUTINE SGEQL2(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEQL2
   END INTERFACE
END MODULE S_SGEQL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEQLF
   INTERFACE
      SUBROUTINE SGEQLF(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEQLF
   END INTERFACE
END MODULE S_SGEQLF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEQP3
   INTERFACE
      SUBROUTINE SGEQP3(M,N,A,Lda,Jpvt,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  INB = 1 , INBMIN = 2 , IXOVER = 3
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      REAL , DIMENSION(*) :: Tau
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEQP3
   END INTERFACE
END MODULE S_SGEQP3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEQR2
   INTERFACE
      SUBROUTINE SGEQR2(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEQR2
   END INTERFACE
END MODULE S_SGEQR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEQR2P
   INTERFACE
      SUBROUTINE SGEQR2P(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEQR2P
   END INTERFACE
END MODULE S_SGEQR2P
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEQR
   INTERFACE
      SUBROUTINE SGEQR(M,N,A,Lda,T,Tsize,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEQR
   END INTERFACE
END MODULE S_SGEQR
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
MODULE S_SGEQRFP
   INTERFACE
      SUBROUTINE SGEQRFP(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEQRFP
   END INTERFACE
END MODULE S_SGEQRFP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEQRT2
   INTERFACE
      SUBROUTINE SGEQRT2(M,N,A,Lda,T,Ldt,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0 , ZERO = 0.0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEQRT2
   END INTERFACE
END MODULE S_SGEQRT2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEQRT3
   INTERFACE
      RECURSIVE SUBROUTINE SGEQRT3(M,N,A,Lda,T,Ldt,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEQRT3
   END INTERFACE
END MODULE S_SGEQRT3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGEQRT
   INTERFACE
      SUBROUTINE SGEQRT(M,N,Nb,A,Lda,T,Ldt,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      LOGICAL , PARAMETER  ::  USE_RECURSIVE_QR = .TRUE.
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGEQRT
   END INTERFACE
END MODULE S_SGEQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGERFS
   INTERFACE
      SUBROUTINE SGERFS(Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,    &
     &                  Ferr,Berr,Work,Iwork,Info)
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGERFS
   END INTERFACE
END MODULE S_SGERFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGERFSX
   INTERFACE
      SUBROUTINE SGERFSX(Trans,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,R,C,B,   &
     &                   Ldb,X,Ldx,Rcond,Berr,N_err_bnds,Err_bnds_norm, &
     &                   Err_bnds_comp,Nparams,Params,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ITREF_DEFAULT = 1.0 ,                       &
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: R
      REAL , DIMENSION(*) :: C
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL , INTENT(INOUT) , DIMENSION(*) :: Params
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGERFSX
   END INTERFACE
END MODULE S_SGERFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGERQ2
   INTERFACE
      SUBROUTINE SGERQ2(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGERQ2
   END INTERFACE
END MODULE S_SGERQ2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGERQF
   INTERFACE
      SUBROUTINE SGERQF(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGERQF
   END INTERFACE
END MODULE S_SGERQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGESC2
   INTERFACE
      SUBROUTINE SGESC2(N,A,Lda,Rhs,Ipiv,Jpiv,Scale)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , TWO = 2.0E+0
      INTEGER :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rhs
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Jpiv
      REAL , INTENT(INOUT) :: Scale
      END SUBROUTINE SGESC2
   END INTERFACE
END MODULE S_SGESC2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGESDD
   INTERFACE
      SUBROUTINE SGESDD(Jobz,M,N,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,Iwork,&
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: S
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGESDD
   END INTERFACE
END MODULE S_SGESDD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGESVD
   INTERFACE
      SUBROUTINE SGESVD(Jobu,Jobvt,M,N,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,&
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobu
      CHARACTER :: Jobvt
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: S
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGESVD
   END INTERFACE
END MODULE S_SGESVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGESVDQ
   INTERFACE
      SUBROUTINE SGESVDQ(Joba,Jobp,Jobr,Jobu,Jobv,M,N,A,Lda,S,U,Ldu,V,  &
     &                   Ldv,Numrank,Iwork,Liwork,Work,Lwork,Rwork,     &
     &                   Lrwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Joba
      CHARACTER :: Jobp
      CHARACTER :: Jobr
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      INTEGER :: M
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: S
      REAL , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(OUT) :: Numrank
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGESVDQ
   END INTERFACE
END MODULE S_SGESVDQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGESVDX
   INTERFACE
      SUBROUTINE SGESVDX(Jobu,Jobvt,Range,M,N,A,Lda,Vl,Vu,Il,Iu,Ns,S,U, &
     &                   Ldu,Vt,Ldvt,Work,Lwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobu
      CHARACTER :: Jobvt
      CHARACTER :: Range
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL :: Vl
      REAL :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      INTEGER , INTENT(INOUT) :: Ns
      REAL , DIMENSION(*) :: S
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGESVDX
   END INTERFACE
END MODULE S_SGESVDX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGESV
   INTERFACE
      SUBROUTINE SGESV(N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGESV
   END INTERFACE
END MODULE S_SGESV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGESVJ
   INTERFACE
      SUBROUTINE SGESVJ(Joba,Jobu,Jobv,M,N,A,Lda,Sva,Mv,V,Ldv,Work,     &
     &                  Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , HALF = 0.5E0 , ONE = 1.0E0
      INTEGER , PARAMETER  ::  NSWEEP = 30
      CHARACTER(1) :: Joba
      CHARACTER(1) :: Jobu
      CHARACTER(1) :: Jobv
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , INTENT(INOUT) , DIMENSION(Lwork) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGESVJ
   END INTERFACE
END MODULE S_SGESVJ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGESVX
   INTERFACE
      SUBROUTINE SGESVX(Fact,Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,R,C, &
     &                  B,Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: R
      REAL , DIMENSION(*) :: C
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGESVX
   END INTERFACE
END MODULE S_SGESVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGESVXX
   INTERFACE
      SUBROUTINE SGESVXX(Fact,Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,R,C,&
     &                   B,Ldb,X,Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,      &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL , INTENT(INOUT) , DIMENSION(*) :: R
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , INTENT(OUT) :: Rpvgrw
      REAL , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL , DIMENSION(*) :: Params
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGESVXX
   END INTERFACE
END MODULE S_SGESVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGETC2
   INTERFACE
      SUBROUTINE SGETC2(N,A,Lda,Ipiv,Jpiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Jpiv
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SGETC2
   END INTERFACE
END MODULE S_SGETC2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGETF2
   INTERFACE
      SUBROUTINE SGETF2(M,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGETF2
   END INTERFACE
END MODULE S_SGETF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGETRF2
   INTERFACE
      RECURSIVE SUBROUTINE SGETRF2(M,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGETRF2
   END INTERFACE
END MODULE S_SGETRF2
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
MODULE S_SGETRI
   INTERFACE
      SUBROUTINE SGETRI(N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGETRI
   END INTERFACE
END MODULE S_SGETRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGETRS
   INTERFACE
      SUBROUTINE SGETRS(Trans,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGETRS
   END INTERFACE
END MODULE S_SGETRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGETSLS
   INTERFACE
      SUBROUTINE SGETSLS(Trans,M,N,Nrhs,A,Lda,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGETSLS
   END INTERFACE
END MODULE S_SGETSLS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGETSQRHRT
   INTERFACE
      SUBROUTINE SGETSQRHRT(M,N,Mb1,Nb1,Nb2,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb1
      INTEGER , INTENT(IN) :: Nb1
      INTEGER , INTENT(IN) :: Nb2
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGETSQRHRT
   END INTERFACE
END MODULE S_SGETSQRHRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGBAK
   INTERFACE
      SUBROUTINE SGGBAK(Job,Side,N,Ilo,Ihi,Lscale,Rscale,M,V,Ldv,Info)
      IMPLICIT NONE
      CHARACTER :: Job
      CHARACTER :: Side
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL , DIMENSION(*) :: Lscale
      REAL , DIMENSION(*) :: Rscale
      INTEGER :: M
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGBAK
   END INTERFACE
END MODULE S_SGGBAK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGBAL
   INTERFACE
      SUBROUTINE SGGBAL(Job,N,A,Lda,B,Ldb,Ilo,Ihi,Lscale,Rscale,Work,   &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , HALF = 0.5E+0 ,             &
     &                      ONE = 1.0E+0 , THREE = 3.0E+0 ,             &
     &                      SCLFAC = 1.0E+1
      CHARACTER :: Job
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Ilo
      INTEGER , INTENT(INOUT) :: Ihi
      REAL , INTENT(INOUT) , DIMENSION(*) :: Lscale
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rscale
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGBAL
   END INTERFACE
END MODULE S_SGGBAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGES3
   INTERFACE
      SUBROUTINE SGGES3(Jobvsl,Jobvsr,Sort,SELCTG,N,A,Lda,B,Ldb,Sdim,   &
     &                  Alphar,Alphai,Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work,    &
     &                  Lwork,Bwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alphar
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      REAL , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGES3
   END INTERFACE
END MODULE S_SGGES3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGES
   INTERFACE
      SUBROUTINE SGGES(Jobvsl,Jobvsr,Sort,SELCTG,N,A,Lda,B,Ldb,Sdim,    &
     &                 Alphar,Alphai,Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work,     &
     &                 Lwork,Bwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alphar
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      REAL , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGES
   END INTERFACE
END MODULE S_SGGES
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGESX
   INTERFACE
      SUBROUTINE SGGESX(Jobvsl,Jobvsr,Sort,SELCTG,Sense,N,A,Lda,B,Ldb,  &
     &                  Sdim,Alphar,Alphai,Beta,Vsl,Ldvsl,Vsr,Ldvsr,    &
     &                  Rconde,Rcondv,Work,Lwork,Iwork,Liwork,Bwork,    &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      CHARACTER :: Sense
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alphar
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      REAL , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      REAL , INTENT(OUT) , DIMENSION(2) :: Rconde
      REAL , INTENT(OUT) , DIMENSION(2) :: Rcondv
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGESX
   END INTERFACE
END MODULE S_SGGESX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGEV3
   INTERFACE
      SUBROUTINE SGGEV3(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alphar,Alphai,Beta,Vl,&
     &                  Ldvl,Vr,Ldvr,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Alphar
      REAL , DIMENSION(*) :: Alphai
      REAL , DIMENSION(*) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGEV3
   END INTERFACE
END MODULE S_SGGEV3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGEV
   INTERFACE
      SUBROUTINE SGGEV(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alphar,Alphai,Beta,Vl, &
     &                 Ldvl,Vr,Ldvr,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Alphar
      REAL , DIMENSION(*) :: Alphai
      REAL , DIMENSION(*) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGEV
   END INTERFACE
END MODULE S_SGGEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGEVX
   INTERFACE
      SUBROUTINE SGGEVX(Balanc,Jobvl,Jobvr,Sense,N,A,Lda,B,Ldb,Alphar,  &
     &                  Alphai,Beta,Vl,Ldvl,Vr,Ldvr,Ilo,Ihi,Lscale,     &
     &                  Rscale,Abnrm,Bbnrm,Rconde,Rcondv,Work,Lwork,    &
     &                  Iwork,Bwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Balanc
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      CHARACTER :: Sense
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Alphar
      REAL , DIMENSION(*) :: Alphai
      REAL , DIMENSION(*) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL , DIMENSION(*) :: Lscale
      REAL , DIMENSION(*) :: Rscale
      REAL , INTENT(INOUT) :: Abnrm
      REAL , INTENT(INOUT) :: Bbnrm
      REAL , DIMENSION(*) :: Rconde
      REAL , DIMENSION(*) :: Rcondv
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGEVX
   END INTERFACE
END MODULE S_SGGEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGGLM
   INTERFACE
      SUBROUTINE SGGGLM(N,M,P,A,Lda,B,Ldb,D,X,Y,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: N
      INTEGER :: M
      INTEGER :: P
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: X
      REAL , DIMENSION(*) :: Y
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGGLM
   END INTERFACE
END MODULE S_SGGGLM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGHD3
   INTERFACE
      SUBROUTINE SGGHD3(Compq,Compz,N,Ilo,Ihi,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,  &
     &                  Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGHD3
   END INTERFACE
END MODULE S_SGGHD3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGHRD
   INTERFACE
      SUBROUTINE SGGHRD(Compq,Compz,N,Ilo,Ihi,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,  &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER :: Ihi
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGHRD
   END INTERFACE
END MODULE S_SGGHRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGLSE
   INTERFACE
      SUBROUTINE SGGLSE(M,N,P,A,Lda,B,Ldb,C,D,X,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: P
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: C
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: X
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGLSE
   END INTERFACE
END MODULE S_SGGLSE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGQRF
   INTERFACE
      SUBROUTINE SGGQRF(N,M,P,A,Lda,Taua,B,Ldb,Taub,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: M
      INTEGER :: P
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Taua
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Taub
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGQRF
   END INTERFACE
END MODULE S_SGGQRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGRQF
   INTERFACE
      SUBROUTINE SGGRQF(M,P,N,A,Lda,Taua,B,Ldb,Taub,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: P
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Taua
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Taub
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGRQF
   END INTERFACE
END MODULE S_SGGRQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGSVD3
   INTERFACE
      SUBROUTINE SGGSVD3(Jobu,Jobv,Jobq,M,N,P,K,L,A,Lda,B,Ldb,Alpha,    &
     &                   Beta,U,Ldu,V,Ldv,Q,Ldq,Work,Lwork,Iwork,Info)
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
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGSVD3
   END INTERFACE
END MODULE S_SGGSVD3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGGSVP3
   INTERFACE
      SUBROUTINE SGGSVP3(Jobu,Jobv,Jobq,M,P,N,A,Lda,B,Ldb,Tola,Tolb,K,L,&
     &                   U,Ldu,V,Ldv,Q,Ldq,Iwork,Tau,Work,Lwork,Info)
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
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGGSVP3
   END INTERFACE
END MODULE S_SGGSVP3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGSVJ0
   INTERFACE
      SUBROUTINE SGSVJ0(Jobv,M,N,A,Lda,D,Sva,Mv,V,Ldv,Eps,Sfmin,Tol,    &
     &                  Nsweep,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , HALF = 0.5E0 , ONE = 1.0E0
      CHARACTER(1) :: Jobv
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(N) :: D
      REAL , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      REAL , INTENT(IN) :: Eps
      REAL , INTENT(IN) :: Sfmin
      REAL , INTENT(IN) :: Tol
      INTEGER , INTENT(IN) :: Nsweep
      REAL , DIMENSION(Lwork) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGSVJ0
   END INTERFACE
END MODULE S_SGSVJ0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGSVJ1
   INTERFACE
      SUBROUTINE SGSVJ1(Jobv,M,N,N1,A,Lda,D,Sva,Mv,V,Ldv,Eps,Sfmin,Tol, &
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(N) :: D
      REAL , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      REAL , INTENT(IN) :: Eps
      REAL , INTENT(IN) :: Sfmin
      REAL , INTENT(IN) :: Tol
      INTEGER , INTENT(IN) :: Nsweep
      REAL , DIMENSION(Lwork) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGSVJ1
   END INTERFACE
END MODULE S_SGSVJ1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGTCON
   INTERFACE
      SUBROUTINE SGTCON(Norm,N,Dl,D,Du,Du2,Ipiv,Anorm,Rcond,Work,Iwork, &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Norm
      INTEGER :: N
      REAL , DIMENSION(*) :: Dl
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: Du
      REAL , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGTCON
   END INTERFACE
END MODULE S_SGTCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGTRFS
   INTERFACE
      SUBROUTINE SGTRFS(Trans,N,Nrhs,Dl,D,Du,Dlf,Df,Duf,Du2,Ipiv,B,Ldb, &
     &                  X,Ldx,Ferr,Berr,Work,Iwork,Info)
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
      REAL , DIMENSION(*) :: Dl
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: Du
      REAL , DIMENSION(*) :: Dlf
      REAL , DIMENSION(*) :: Df
      REAL , DIMENSION(*) :: Duf
      REAL , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGTRFS
   END INTERFACE
END MODULE S_SGTRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGTSV
   INTERFACE
      SUBROUTINE SGTSV(N,Nrhs,Dl,D,Du,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dl
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: Du
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGTSV
   END INTERFACE
END MODULE S_SGTSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGTSVX
   INTERFACE
      SUBROUTINE SGTSVX(Fact,Trans,N,Nrhs,Dl,D,Du,Dlf,Df,Duf,Du2,Ipiv,B,&
     &                  Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: Dl
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: Du
      REAL , DIMENSION(*) :: Dlf
      REAL , DIMENSION(*) :: Df
      REAL , DIMENSION(*) :: Duf
      REAL , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGTSVX
   END INTERFACE
END MODULE S_SGTSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGTTRF
   INTERFACE
      SUBROUTINE SGTTRF(N,Dl,D,Du,Du2,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dl
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: Du
      REAL , INTENT(OUT) , DIMENSION(*) :: Du2
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER :: Info
      END SUBROUTINE SGTTRF
   END INTERFACE
END MODULE S_SGTTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGTTRS
   INTERFACE
      SUBROUTINE SGTTRS(Trans,N,Nrhs,Dl,D,Du,Du2,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: Dl
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: Du
      REAL , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SGTTRS
   END INTERFACE
END MODULE S_SGTTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SGTTS2
   INTERFACE
      SUBROUTINE SGTTS2(Itrans,N,Nrhs,Dl,D,Du,Du2,Ipiv,B,Ldb)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: Itrans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL , INTENT(IN) , DIMENSION(*) :: Dl
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: Du
      REAL , INTENT(IN) , DIMENSION(*) :: Du2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE SGTTS2
   END INTERFACE
END MODULE S_SGTTS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SHGEQZ
   INTERFACE
      SUBROUTINE SHGEQZ(Job,Compq,Compz,N,Ilo,Ihi,H,Ldh,T,Ldt,Alphar,   &
     &                  Alphai,Beta,Q,Ldq,Z,Ldz,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  HALF = 0.5E+0 , ZERO = 0.0E+0 ,             &
     &                      ONE = 1.0E+0 , SAFETY = 1.0E+2
      CHARACTER :: Job
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , INTENT(OUT) , DIMENSION(*) :: Alphar
      REAL , INTENT(OUT) , DIMENSION(*) :: Alphai
      REAL , INTENT(OUT) , DIMENSION(*) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SHGEQZ
   END INTERFACE
END MODULE S_SHGEQZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SHSEIN
   INTERFACE
      SUBROUTINE SHSEIN(Side,Eigsrc,Initv,Select,N,H,Ldh,Wr,Wi,Vl,Ldvl, &
     &                  Vr,Ldvr,Mm,M,Work,Ifaill,Ifailr,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Side
      CHARACTER :: Eigsrc
      CHARACTER :: Initv
      LOGICAL , INTENT(INOUT) , DIMENSION(*) :: Select
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL , INTENT(INOUT) , DIMENSION(*) :: Wr
      REAL , INTENT(IN) , DIMENSION(*) :: Wi
      REAL , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      REAL , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifaill
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifailr
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SHSEIN
   END INTERFACE
END MODULE S_SHSEIN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SHSEQR
   INTERFACE
      SUBROUTINE SHSEQR(Job,Compz,N,Ilo,Ihi,H,Ldh,Wr,Wi,Z,Ldz,Work,     &
     &                  Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NTINY = 15 , NL = 49
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Job
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL , DIMENSION(*) :: Wr
      REAL , DIMENSION(*) :: Wi
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SHSEQR
   END INTERFACE
END MODULE S_SHSEQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SISNAN
   INTERFACE
      FUNCTION SISNAN(Sin)
      IMPLICIT NONE
      LOGICAL :: SISNAN
      REAL , INTENT(IN) :: Sin
      END FUNCTION SISNAN
   END INTERFACE
END MODULE S_SISNAN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLABAD
   INTERFACE
      SUBROUTINE SLABAD(Small,Large)
      IMPLICIT NONE
      REAL , INTENT(INOUT) :: Small
      REAL , INTENT(INOUT) :: Large
      END SUBROUTINE SLABAD
   END INTERFACE
END MODULE S_SLABAD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLABRD
   INTERFACE
      SUBROUTINE SLABRD(M,N,Nb,A,Lda,D,E,Tauq,Taup,X,Ldx,Y,Ldy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(OUT) , DIMENSION(*) :: D
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Tauq
      REAL , DIMENSION(*) :: Taup
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , DIMENSION(Ldy,*) :: Y
      INTEGER :: Ldy
      END SUBROUTINE SLABRD
   END INTERFACE
END MODULE S_SLABRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLACN2
   INTERFACE
      SUBROUTINE SLACN2(N,V,X,Isgn,Est,Kase,Isave)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: V
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isgn
      REAL , INTENT(INOUT) :: Est
      INTEGER , INTENT(INOUT) :: Kase
      INTEGER , INTENT(INOUT) , DIMENSION(3) :: Isave
      END SUBROUTINE SLACN2
   END INTERFACE
END MODULE S_SLACN2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLACON
   INTERFACE
      SUBROUTINE SLACON(N,V,X,Isgn,Est,Kase)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: V
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isgn
      REAL , INTENT(INOUT) :: Est
      INTEGER , INTENT(INOUT) :: Kase
      END SUBROUTINE SLACON
   END INTERFACE
END MODULE S_SLACON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLACPY
   INTERFACE
      SUBROUTINE SLACPY(Uplo,M,N,A,Lda,B,Ldb)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(OUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE SLACPY
   END INTERFACE
END MODULE S_SLACPY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLADIV
   INTERFACE
      SUBROUTINE SLADIV(A,B,C,D,P,Q)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  BS = 2.0E0 , HALF = 0.5E0 , TWO = 2.0E0
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(INOUT) :: P
      REAL , INTENT(INOUT) :: Q
      END SUBROUTINE SLADIV
   END INTERFACE
END MODULE S_SLADIV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLADIV1
   INTERFACE
      SUBROUTINE SLADIV1(A,B,C,D,P,Q)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0
      REAL , INTENT(INOUT) :: A
      REAL :: B
      REAL :: C
      REAL :: D
      REAL , INTENT(OUT) :: P
      REAL , INTENT(OUT) :: Q
      END SUBROUTINE SLADIV1
   END INTERFACE
END MODULE S_SLADIV1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLADIV2
   INTERFACE
      FUNCTION SLADIV2(A,B,C,D,R,T)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0
      REAL :: SLADIV2
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(IN) :: R
      REAL , INTENT(IN) :: T
      END FUNCTION SLADIV2
   END INTERFACE
END MODULE S_SLADIV2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAE2
   INTERFACE
      SUBROUTINE SLAE2(A,B,C,Rt1,Rt2)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , TWO = 2.0E0 , ZERO = 0.0E0 ,  &
     &                      HALF = 0.5E0
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(INOUT) :: Rt1
      REAL , INTENT(OUT) :: Rt2
      END SUBROUTINE SLAE2
   END INTERFACE
END MODULE S_SLAE2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAEBZ
   INTERFACE
      SUBROUTINE SLAEBZ(Ijob,Nitmax,N,Mmax,Minp,Nbmin,Abstol,Reltol,    &
     &                  Pivmin,D,E,E2,Nval,Ab,C,Mout,Nab,Work,Iwork,    &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , TWO = 2.0E0 ,                &
     &                      HALF = 1.0E0/TWO
      INTEGER , INTENT(IN) :: Ijob
      INTEGER , INTENT(IN) :: Nitmax
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Mmax
      INTEGER , INTENT(IN) :: Minp
      INTEGER , INTENT(IN) :: Nbmin
      REAL , INTENT(IN) :: Abstol
      REAL , INTENT(IN) :: Reltol
      REAL , INTENT(IN) :: Pivmin
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , INTENT(IN) , DIMENSION(*) :: E2
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Nval
      REAL , INTENT(INOUT) , DIMENSION(Mmax,*) :: Ab
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      INTEGER , INTENT(INOUT) :: Mout
      INTEGER , INTENT(INOUT) , DIMENSION(Mmax,*) :: Nab
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAEBZ
   END INTERFACE
END MODULE S_SLAEBZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAED0
   INTERFACE
      SUBROUTINE SLAED0(Icompq,Qsiz,N,D,E,Q,Ldq,Qstore,Ldqs,Work,Iwork, &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.E0 , ONE = 1.E0 , TWO = 2.E0
      INTEGER :: Icompq
      INTEGER :: Qsiz
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , DIMENSION(Ldqs,*) :: Qstore
      INTEGER :: Ldqs
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAED0
   END INTERFACE
END MODULE S_SLAED0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAED1
   INTERFACE
      SUBROUTINE SLAED1(N,D,Q,Ldq,Indxq,Rho,Cutpnt,Work,Iwork,Info)
      IMPLICIT NONE
      INTEGER :: N
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , DIMENSION(*) :: Indxq
      REAL :: Rho
      INTEGER :: Cutpnt
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAED1
   END INTERFACE
END MODULE S_SLAED1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAED2
   INTERFACE
      SUBROUTINE SLAED2(K,N,N1,D,Q,Ldq,Indxq,Rho,Z,Dlamda,W,Q2,Indx,    &
     &                  Indxc,Indxp,Coltyp,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  MONE = -1.0E0 , ZERO = 0.0E0 , ONE = 1.0E0 ,&
     &                      TWO = 2.0E0 , EIGHT = 8.0E0
      INTEGER , INTENT(INOUT) :: K
      INTEGER :: N
      INTEGER :: N1
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxq
      REAL , INTENT(INOUT) :: Rho
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL , DIMENSION(*) :: Dlamda
      REAL , INTENT(OUT) , DIMENSION(*) :: W
      REAL , DIMENSION(*) :: Q2
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxc
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxp
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Coltyp
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAED2
   END INTERFACE
END MODULE S_SLAED2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAED3
   INTERFACE
      SUBROUTINE SLAED3(K,N,N1,D,Q,Ldq,Rho,Dlamda,Q2,Indx,Ctot,W,S,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , ZERO = 0.0E0
      INTEGER :: K
      INTEGER , INTENT(IN) :: N
      INTEGER :: N1
      REAL , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL :: Rho
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dlamda
      REAL , DIMENSION(*) :: Q2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Indx
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ctot
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAED3
   END INTERFACE
END MODULE S_SLAED3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAED4
   INTERFACE
      SUBROUTINE SLAED4(N,I,D,Z,Delta,Rho,Dlam,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXIT = 30
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      THREE = 3.0E0 , FOUR = 4.0E0 ,              &
     &                      EIGHT = 8.0E0 , TEN = 10.0E0
      INTEGER , INTENT(IN) :: N
      INTEGER :: I
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: Z
      REAL , INTENT(INOUT) , DIMENSION(*) :: Delta
      REAL :: Rho
      REAL :: Dlam
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAED4
   END INTERFACE
END MODULE S_SLAED4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAED5
   INTERFACE
      SUBROUTINE SLAED5(I,D,Z,Delta,Rho,Dlam)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      FOUR = 4.0E0
      INTEGER , INTENT(IN) :: I
      REAL , INTENT(IN) , DIMENSION(2) :: D
      REAL , INTENT(IN) , DIMENSION(2) :: Z
      REAL , INTENT(INOUT) , DIMENSION(2) :: Delta
      REAL , INTENT(IN) :: Rho
      REAL , INTENT(OUT) :: Dlam
      END SUBROUTINE SLAED5
   END INTERFACE
END MODULE S_SLAED5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAED6
   INTERFACE
      SUBROUTINE SLAED6(Kniter,Orgati,Rho,D,Z,Finit,Tau,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXIT = 40
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      THREE = 3.0E0 , FOUR = 4.0E0 , EIGHT = 8.0E0
      INTEGER , INTENT(IN) :: Kniter
      LOGICAL , INTENT(IN) :: Orgati
      REAL , INTENT(IN) :: Rho
      REAL , INTENT(IN) , DIMENSION(3) :: D
      REAL , INTENT(IN) , DIMENSION(3) :: Z
      REAL , INTENT(IN) :: Finit
      REAL , INTENT(INOUT) :: Tau
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLAED6
   END INTERFACE
END MODULE S_SLAED6
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAED7
   INTERFACE
      SUBROUTINE SLAED7(Icompq,N,Qsiz,Tlvls,Curlvl,Curpbm,D,Q,Ldq,Indxq,&
     &                  Rho,Cutpnt,Qstore,Qptr,Prmptr,Perm,Givptr,      &
     &                  Givcol,Givnum,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , ZERO = 0.0E0
      INTEGER :: Icompq
      INTEGER :: N
      INTEGER :: Qsiz
      INTEGER :: Tlvls
      INTEGER :: Curlvl
      INTEGER :: Curpbm
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , DIMENSION(*) :: Indxq
      REAL :: Rho
      INTEGER :: Cutpnt
      REAL , DIMENSION(*) :: Qstore
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Qptr
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Prmptr
      INTEGER , DIMENSION(*) :: Perm
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Givptr
      INTEGER , DIMENSION(2,*) :: Givcol
      REAL , DIMENSION(2,*) :: Givnum
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAED7
   END INTERFACE
END MODULE S_SLAED7
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAED8
   INTERFACE
      SUBROUTINE SLAED8(Icompq,K,N,Qsiz,D,Q,Ldq,Indxq,Rho,Cutpnt,Z,     &
     &                  Dlamda,Q2,Ldq2,W,Perm,Givptr,Givcol,Givnum,     &
     &                  Indxp,Indx,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  MONE = -1.0E0 , ZERO = 0.0E0 , ONE = 1.0E0 ,&
     &                      TWO = 2.0E0 , EIGHT = 8.0E0
      INTEGER , INTENT(IN) :: Icompq
      INTEGER , INTENT(INOUT) :: K
      INTEGER :: N
      INTEGER :: Qsiz
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxq
      REAL , INTENT(INOUT) :: Rho
      INTEGER , INTENT(IN) :: Cutpnt
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dlamda
      REAL , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Perm
      INTEGER , INTENT(INOUT) :: Givptr
      INTEGER , INTENT(OUT) , DIMENSION(2,*) :: Givcol
      REAL , INTENT(OUT) , DIMENSION(2,*) :: Givnum
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxp
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indx
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAED8
   END INTERFACE
END MODULE S_SLAED8
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAED9
   INTERFACE
      SUBROUTINE SLAED9(K,Kstart,Kstop,N,D,Q,Ldq,Rho,Dlamda,W,S,Lds,    &
     &                  Info)
      IMPLICIT NONE
      INTEGER :: K
      INTEGER , INTENT(IN) :: Kstart
      INTEGER , INTENT(IN) :: Kstop
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER , INTENT(IN) :: Ldq
      REAL :: Rho
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dlamda
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , INTENT(INOUT) , DIMENSION(Lds,*) :: S
      INTEGER , INTENT(IN) :: Lds
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAED9
   END INTERFACE
END MODULE S_SLAED9
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAEDA
   INTERFACE
      SUBROUTINE SLAEDA(N,Tlvls,Curlvl,Curpbm,Prmptr,Perm,Givptr,Givcol,&
     &                  Givnum,Q,Qptr,Z,Ztemp,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , HALF = 0.5E0 , ONE = 1.0E0
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Tlvls
      INTEGER , INTENT(IN) :: Curlvl
      INTEGER , INTENT(IN) :: Curpbm
      INTEGER , INTENT(IN) , DIMENSION(*) :: Prmptr
      INTEGER , INTENT(IN) , DIMENSION(*) :: Perm
      INTEGER , INTENT(IN) , DIMENSION(*) :: Givptr
      INTEGER , INTENT(IN) , DIMENSION(2,*) :: Givcol
      REAL , DIMENSION(2,*) :: Givnum
      REAL , DIMENSION(*) :: Q
      INTEGER , INTENT(IN) , DIMENSION(*) :: Qptr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL , DIMENSION(*) :: Ztemp
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAEDA
   END INTERFACE
END MODULE S_SLAEDA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAEIN
   INTERFACE
      SUBROUTINE SLAEIN(Rightv,Noinit,N,H,Ldh,Wr,Wi,Vr,Vi,B,Ldb,Work,   &
     &                  Eps3,Smlnum,Bignum,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TENTH = 1.0E-1
      LOGICAL , INTENT(IN) :: Rightv
      LOGICAL , INTENT(IN) :: Noinit
      INTEGER :: N
      REAL , INTENT(IN) , DIMENSION(Ldh,*) :: H
      INTEGER , INTENT(IN) :: Ldh
      REAL , INTENT(IN) :: Wr
      REAL , INTENT(IN) :: Wi
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vi
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(IN) :: Eps3
      REAL , INTENT(IN) :: Smlnum
      REAL , INTENT(IN) :: Bignum
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLAEIN
   END INTERFACE
END MODULE S_SLAEIN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAEV2
   INTERFACE
      SUBROUTINE SLAEV2(A,B,C,Rt1,Rt2,Cs1,Sn1)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , TWO = 2.0E0 , ZERO = 0.0E0 ,  &
     &                      HALF = 0.5E0
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(INOUT) :: Rt1
      REAL , INTENT(OUT) :: Rt2
      REAL , INTENT(INOUT) :: Cs1
      REAL , INTENT(INOUT) :: Sn1
      END SUBROUTINE SLAEV2
   END INTERFACE
END MODULE S_SLAEV2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAEXC
   INTERFACE
      SUBROUTINE SLAEXC(Wantq,N,T,Ldt,Q,Ldq,J1,N1,N2,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TEN = 1.0E+1
      INTEGER , PARAMETER  ::  LDD = 4 , LDX = 2
      LOGICAL , INTENT(IN) :: Wantq
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , INTENT(IN) :: J1
      INTEGER :: N1
      INTEGER :: N2
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLAEXC
   END INTERFACE
END MODULE S_SLAEXC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAG2D
   INTERFACE
      SUBROUTINE SLAG2D(M,N,Sa,Ldsa,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Ldsa,*) :: Sa
      INTEGER , INTENT(IN) :: Ldsa
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLAG2D
   END INTERFACE
END MODULE S_SLAG2D
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAG2
   INTERFACE
      SUBROUTINE SLAG2(A,Lda,B,Ldb,Safmin,Scale1,Scale2,Wr1,Wr2,Wi)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , HALF = ONE/TWO ,             &
     &                      FUZZY1 = ONE + 1.0E-5
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , INTENT(IN) :: Safmin
      REAL , INTENT(INOUT) :: Scale1
      REAL , INTENT(OUT) :: Scale2
      REAL , INTENT(INOUT) :: Wr1
      REAL , INTENT(INOUT) :: Wr2
      REAL , INTENT(INOUT) :: Wi
      END SUBROUTINE SLAG2
   END INTERFACE
END MODULE S_SLAG2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_GBAMV
   INTERFACE
      SUBROUTINE SLA_GBAMV(Trans,M,N,Kl,Ku,Alpha,Ab,Ldab,X,Incx,Beta,Y, &
     &                     Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE SLA_GBAMV
   END INTERFACE
END MODULE S_SLA_GBAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_GBRCOND
   INTERFACE
      FUNCTION SLA_GBRCOND(Trans,N,Kl,Ku,Ab,Ldab,Afb,Ldafb,Ipiv,Cmode,C,&
     &                     Info,Work,Iwork)
      IMPLICIT NONE
      REAL :: SLA_GBRCOND
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      REAL , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(IN) :: Cmode
      REAL , INTENT(IN) , DIMENSION(*) :: C
      INTEGER , INTENT(INOUT) :: Info
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      END FUNCTION SLA_GBRCOND
   END INTERFACE
END MODULE S_SLA_GBRCOND
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_GBRFSX_EXTENDED
   INTERFACE
      SUBROUTINE SLA_GBRFSX_EXTENDED(Prec_type,Trans_type,N,Kl,Ku,Nrhs, &
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
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      REAL , DIMENSION(*) :: Res
      REAL , DIMENSION(*) :: Ayb
      REAL , DIMENSION(*) :: Dy
      REAL , DIMENSION(*) :: Y_tail
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL , INTENT(IN) :: Rthresh
      REAL , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
      END SUBROUTINE SLA_GBRFSX_EXTENDED
   END INTERFACE
END MODULE S_SLA_GBRFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_GBRPVGRW
   INTERFACE
      FUNCTION SLA_GBRPVGRW(N,Kl,Ku,Ncols,Ab,Ldab,Afb,Ldafb)
      IMPLICIT NONE
      REAL :: SLA_GBRPVGRW
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      INTEGER , INTENT(IN) :: Ncols
      REAL , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(IN) , DIMENSION(Ldafb,*) :: Afb
      INTEGER , INTENT(IN) :: Ldafb
      END FUNCTION SLA_GBRPVGRW
   END INTERFACE
END MODULE S_SLA_GBRPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_GEAMV
   INTERFACE
      SUBROUTINE SLA_GEAMV(Trans,M,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: Trans
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
      END SUBROUTINE SLA_GEAMV
   END INTERFACE
END MODULE S_SLA_GEAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_GERCOND
   INTERFACE
      FUNCTION SLA_GERCOND(Trans,N,A,Lda,Af,Ldaf,Ipiv,Cmode,C,Info,Work,&
     &                     Iwork)
      IMPLICIT NONE
      REAL :: SLA_GERCOND
      CHARACTER :: Trans
      INTEGER :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(IN) :: Cmode
      REAL , INTENT(IN) , DIMENSION(*) :: C
      INTEGER , INTENT(INOUT) :: Info
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      END FUNCTION SLA_GERCOND
   END INTERFACE
END MODULE S_SLA_GERCOND
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_GERFSX_EXTENDED
   INTERFACE
      SUBROUTINE SLA_GERFSX_EXTENDED(Prec_type,Trans_type,N,Nrhs,A,Lda, &
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Errs_n
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Errs_c
      REAL , DIMENSION(*) :: Res
      REAL , DIMENSION(*) :: Ayb
      REAL , DIMENSION(*) :: Dy
      REAL , DIMENSION(*) :: Y_tail
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL , INTENT(IN) :: Rthresh
      REAL , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
      END SUBROUTINE SLA_GERFSX_EXTENDED
   END INTERFACE
END MODULE S_SLA_GERFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_GERPVGRW
   INTERFACE
      FUNCTION SLA_GERPVGRW(N,Ncols,A,Lda,Af,Ldaf)
      IMPLICIT NONE
      REAL :: SLA_GERPVGRW
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ncols
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      END FUNCTION SLA_GERPVGRW
   END INTERFACE
END MODULE S_SLA_GERPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAGS2
   INTERFACE
      SUBROUTINE SLAGS2(Upper,A1,A2,A3,B1,B2,B3,Csu,Snu,Csv,Snv,Csq,Snq)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      LOGICAL , INTENT(IN) :: Upper
      REAL , INTENT(IN) :: A1
      REAL , INTENT(IN) :: A2
      REAL , INTENT(IN) :: A3
      REAL , INTENT(IN) :: B1
      REAL , INTENT(IN) :: B2
      REAL , INTENT(IN) :: B3
      REAL , INTENT(OUT) :: Csu
      REAL , INTENT(OUT) :: Snu
      REAL , INTENT(OUT) :: Csv
      REAL , INTENT(OUT) :: Snv
      REAL :: Csq
      REAL :: Snq
      END SUBROUTINE SLAGS2
   END INTERFACE
END MODULE S_SLAGS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAGTF
   INTERFACE
      SUBROUTINE SLAGTF(N,A,Lambda,B,C,Tol,D,In,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: A
      REAL , INTENT(IN) :: Lambda
      REAL , INTENT(INOUT) , DIMENSION(*) :: B
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      REAL , INTENT(IN) :: Tol
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: In
      INTEGER :: Info
      END SUBROUTINE SLAGTF
   END INTERFACE
END MODULE S_SLAGTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAGTM
   INTERFACE
      SUBROUTINE SLAGTM(Trans,N,Nrhs,Alpha,Dl,D,Du,X,Ldx,Beta,B,Ldb)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(*) :: Dl
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: Du
      REAL , INTENT(IN) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE SLAGTM
   END INTERFACE
END MODULE S_SLAGTM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAGTS
   INTERFACE
      SUBROUTINE SLAGTS(Job,N,A,B,C,D,In,Y,Tol,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: Job
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: A
      REAL , INTENT(IN) , DIMENSION(*) :: B
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , INTENT(IN) , DIMENSION(*) :: D
      INTEGER , INTENT(IN) , DIMENSION(*) :: In
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      REAL , INTENT(INOUT) :: Tol
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAGTS
   END INTERFACE
END MODULE S_SLAGTS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAGV2
   INTERFACE
      SUBROUTINE SLAGV2(A,Lda,B,Ldb,Alphar,Alphai,Beta,Csl,Snl,Csr,Snr)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(2) :: Alphar
      REAL , INTENT(INOUT) , DIMENSION(2) :: Alphai
      REAL , INTENT(OUT) , DIMENSION(2) :: Beta
      REAL :: Csl
      REAL :: Snl
      REAL :: Csr
      REAL , INTENT(INOUT) :: Snr
      END SUBROUTINE SLAGV2
   END INTERFACE
END MODULE S_SLAGV2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAHQR
   INTERFACE
      SUBROUTINE SLAHQR(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,Wr,Wi,Iloz,Ihiz,Z,  &
     &                  Ldz,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      DAT1 = 3.0E0/4.0E0 , DAT2 = -0.4375E0
      INTEGER , PARAMETER  ::  KEXSH = 10
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL , DIMENSION(*) :: Wr
      REAL , DIMENSION(*) :: Wi
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      REAL , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLAHQR
   END INTERFACE
END MODULE S_SLAHQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAHR2
   INTERFACE
      SUBROUTINE SLAHR2(N,K,Nb,A,Lda,Tau,T,Ldt,Y,Ldy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , INTENT(IN) :: N
      INTEGER :: K
      INTEGER :: Nb
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Nb) :: Tau
      REAL , DIMENSION(Ldt,Nb) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Ldy,Nb) :: Y
      INTEGER :: Ldy
      END SUBROUTINE SLAHR2
   END INTERFACE
END MODULE S_SLAHR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAIC1
   INTERFACE
      SUBROUTINE SLAIC1(Job,J,X,Sest,W,Gamma,Sestpr,S,C)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      HALF = 0.5E0 , FOUR = 4.0E0
      INTEGER , INTENT(IN) :: Job
      INTEGER :: J
      REAL , DIMENSION(J) :: X
      REAL , INTENT(IN) :: Sest
      REAL , DIMENSION(J) :: W
      REAL , INTENT(IN) :: Gamma
      REAL , INTENT(OUT) :: Sestpr
      REAL , INTENT(INOUT) :: S
      REAL , INTENT(INOUT) :: C
      END SUBROUTINE SLAIC1
   END INTERFACE
END MODULE S_SLAIC1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAISNAN
   INTERFACE
      FUNCTION SLAISNAN(Sin1,Sin2)
      IMPLICIT NONE
      LOGICAL :: SLAISNAN
      REAL , INTENT(IN) :: Sin1
      REAL , INTENT(IN) :: Sin2
      END FUNCTION SLAISNAN
   END INTERFACE
END MODULE S_SLAISNAN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_LIN_BERR
   INTERFACE
      SUBROUTINE SLA_LIN_BERR(N,Nz,Nrhs,Res,Ayb,Berr)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nz
      INTEGER , INTENT(IN) :: Nrhs
      REAL , INTENT(IN) , DIMENSION(N,Nrhs) :: Res
      REAL , INTENT(IN) , DIMENSION(N,Nrhs) :: Ayb
      REAL , INTENT(INOUT) , DIMENSION(Nrhs) :: Berr
      END SUBROUTINE SLA_LIN_BERR
   END INTERFACE
END MODULE S_SLA_LIN_BERR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLALN2
   INTERFACE
      SUBROUTINE SLALN2(Ltrans,Na,Nw,Smin,Ca,A,Lda,D1,D2,B,Ldb,Wr,Wi,X, &
     &                  Ldx,Scale,Xnorm,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
      LOGICAL , INTENT(IN) :: Ltrans
      INTEGER , INTENT(IN) :: Na
      INTEGER , INTENT(IN) :: Nw
      REAL , INTENT(IN) :: Smin
      REAL , INTENT(IN) :: Ca
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) :: D1
      REAL , INTENT(IN) :: D2
      REAL , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , INTENT(IN) :: Wr
      REAL , INTENT(IN) :: Wi
      REAL , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(INOUT) :: Xnorm
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLALN2
   END INTERFACE
END MODULE S_SLALN2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLALS0
   INTERFACE
      SUBROUTINE SLALS0(Icompq,Nl,Nr,Sqre,Nrhs,B,Ldb,Bx,Ldbx,Perm,      &
     &                  Givptr,Givcol,Ldgcol,Givnum,Ldgnum,Poles,Difl,  &
     &                  Difr,Z,K,C,S,Work,Info)
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
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldbx,*) :: Bx
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
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLALS0
   END INTERFACE
END MODULE S_SLALS0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLALSA
   INTERFACE
      SUBROUTINE SLALSA(Icompq,Smlsiz,N,Nrhs,B,Ldb,Bx,Ldbx,U,Ldu,Vt,K,  &
     &                  Difl,Difr,Z,Poles,Givptr,Givcol,Ldgcol,Perm,    &
     &                  Givnum,C,S,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      INTEGER :: Icompq
      INTEGER :: Smlsiz
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldbx,*) :: Bx
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
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLALSA
   END INTERFACE
END MODULE S_SLALSA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLALSD
   INTERFACE
      SUBROUTINE SLALSD(Uplo,Smlsiz,N,Nrhs,D,E,B,Ldb,Rcond,Rank,Work,   &
     &                  Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
      CHARACTER , INTENT(IN) :: Uplo
      INTEGER :: Smlsiz
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLALSD
   END INTERFACE
END MODULE S_SLALSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAMRG
   INTERFACE
      SUBROUTINE SLAMRG(N1,N2,A,Strd1,Strd2,Index)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N1
      INTEGER , INTENT(IN) :: N2
      REAL , INTENT(IN) , DIMENSION(*) :: A
      INTEGER , INTENT(IN) :: Strd1
      INTEGER , INTENT(IN) :: Strd2
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Index
      END SUBROUTINE SLAMRG
   END INTERFACE
END MODULE S_SLAMRG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAMSWLQ
   INTERFACE
      SUBROUTINE SLAMSWLQ(Side,Trans,M,N,K,Mb,Nb,A,Lda,T,Ldt,C,Ldc,Work,&
     &                    Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: Mb
      INTEGER :: Nb
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAMSWLQ
   END INTERFACE
END MODULE S_SLAMSWLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAMTSQR
   INTERFACE
      SUBROUTINE SLAMTSQR(Side,Trans,M,N,K,Mb,Nb,A,Lda,T,Ldt,C,Ldc,Work,&
     &                    Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: Mb
      INTEGER :: Nb
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAMTSQR
   END INTERFACE
END MODULE S_SLAMTSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANEG
   INTERFACE
      FUNCTION SLANEG(N,D,Lld,Sigma,Pivmin,R)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      INTEGER , PARAMETER  ::  BLKLEN = 128
      INTEGER :: SLANEG
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: Lld
      REAL , INTENT(IN) :: Sigma
      REAL :: Pivmin
      INTEGER , INTENT(IN) :: R
      END FUNCTION SLANEG
   END INTERFACE
END MODULE S_SLANEG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANGB
   INTERFACE
      FUNCTION SLANGB(Norm,N,Kl,Ku,Ab,Ldab,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANGB
      CHARACTER :: Norm
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION SLANGB
   END INTERFACE
END MODULE S_SLANGB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANGE
   INTERFACE
      FUNCTION SLANGE(Norm,M,N,A,Lda,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANGE
      CHARACTER :: Norm
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION SLANGE
   END INTERFACE
END MODULE S_SLANGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANGT
   INTERFACE
      FUNCTION SLANGT(Norm,N,Dl,D,Du)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANGT
      CHARACTER :: Norm
      INTEGER :: N
      REAL , DIMENSION(*) :: Dl
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: Du
      END FUNCTION SLANGT
   END INTERFACE
END MODULE S_SLANGT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANHS
   INTERFACE
      FUNCTION SLANHS(Norm,N,A,Lda,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANHS
      CHARACTER :: Norm
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION SLANHS
   END INTERFACE
END MODULE S_SLANHS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANSB
   INTERFACE
      FUNCTION SLANSB(Norm,Uplo,N,K,Ab,Ldab,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANSB
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION SLANSB
   END INTERFACE
END MODULE S_SLANSB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANSF
   INTERFACE
      FUNCTION SLANSF(Norm,Transr,Uplo,N,A,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANSF
      CHARACTER :: Norm
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(0:*) :: A
      REAL , INTENT(INOUT) , DIMENSION(0:*) :: Work
      END FUNCTION SLANSF
   END INTERFACE
END MODULE S_SLANSF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANSP
   INTERFACE
      FUNCTION SLANSP(Norm,Uplo,N,Ap,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANSP
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(*) :: Ap
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION SLANSP
   END INTERFACE
END MODULE S_SLANSP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANST
   INTERFACE
      FUNCTION SLANST(Norm,N,D,E)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANST
      CHARACTER :: Norm
      INTEGER :: N
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      END FUNCTION SLANST
   END INTERFACE
END MODULE S_SLANST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANSY
   INTERFACE
      FUNCTION SLANSY(Norm,Uplo,N,A,Lda,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANSY
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION SLANSY
   END INTERFACE
END MODULE S_SLANSY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANTB
   INTERFACE
      FUNCTION SLANTB(Norm,Uplo,Diag,N,K,Ab,Ldab,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANTB
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION SLANTB
   END INTERFACE
END MODULE S_SLANTB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANTP
   INTERFACE
      FUNCTION SLANTP(Norm,Uplo,Diag,N,Ap,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANTP
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(*) :: Ap
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION SLANTP
   END INTERFACE
END MODULE S_SLANTP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANTR
   INTERFACE
      FUNCTION SLANTR(Norm,Uplo,Diag,M,N,A,Lda,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SLANTR
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION SLANTR
   END INTERFACE
END MODULE S_SLANTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLANV2
   INTERFACE
      SUBROUTINE SLANV2(A,B,C,D,Rt1r,Rt1i,Rt2r,Rt2i,Cs,Sn)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , HALF = 0.5E+0 ,             &
     &                      ONE = 1.0E+0 , TWO = 2.0E+0 ,               &
     &                      MULTPL = 4.0E+0
      REAL , INTENT(INOUT) :: A
      REAL , INTENT(INOUT) :: B
      REAL , INTENT(INOUT) :: C
      REAL , INTENT(INOUT) :: D
      REAL , INTENT(OUT) :: Rt1r
      REAL , INTENT(INOUT) :: Rt1i
      REAL , INTENT(OUT) :: Rt2r
      REAL , INTENT(OUT) :: Rt2i
      REAL , INTENT(INOUT) :: Cs
      REAL , INTENT(INOUT) :: Sn
      END SUBROUTINE SLANV2
   END INTERFACE
END MODULE S_SLANV2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAORHR_COL_GETRFNP2
   INTERFACE
      RECURSIVE SUBROUTINE SLAORHR_COL_GETRFNP2(M,N,A,Lda,D,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAORHR_COL_GETRFNP2
   END INTERFACE
END MODULE S_SLAORHR_COL_GETRFNP2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAORHR_COL_GETRFNP
   INTERFACE
      SUBROUTINE SLAORHR_COL_GETRFNP(M,N,A,Lda,D,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAORHR_COL_GETRFNP
   END INTERFACE
END MODULE S_SLAORHR_COL_GETRFNP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAPLL
   INTERFACE
      SUBROUTINE SLAPLL(N,X,Incx,Y,Incy,Ssmin)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER :: Incx
      REAL , DIMENSION(*) :: Y
      INTEGER :: Incy
      REAL :: Ssmin
      END SUBROUTINE SLAPLL
   END INTERFACE
END MODULE S_SLAPLL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAPMR
   INTERFACE
      SUBROUTINE SLAPMR(Forwrd,M,N,X,Ldx,K)
      IMPLICIT NONE
      LOGICAL , INTENT(IN) :: Forwrd
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: K
      END SUBROUTINE SLAPMR
   END INTERFACE
END MODULE S_SLAPMR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAPMT
   INTERFACE
      SUBROUTINE SLAPMT(Forwrd,M,N,X,Ldx,K)
      IMPLICIT NONE
      LOGICAL , INTENT(IN) :: Forwrd
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: K
      END SUBROUTINE SLAPMT
   END INTERFACE
END MODULE S_SLAPMT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_PORCOND
   INTERFACE
      FUNCTION SLA_PORCOND(Uplo,N,A,Lda,Af,Ldaf,Cmode,C,Info,Work,Iwork)
      IMPLICIT NONE
      REAL :: SLA_PORCOND
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , INTENT(IN) :: Cmode
      REAL , INTENT(IN) , DIMENSION(*) :: C
      INTEGER , INTENT(INOUT) :: Info
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      END FUNCTION SLA_PORCOND
   END INTERFACE
END MODULE S_SLA_PORCOND
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_PORFSX_EXTENDED
   INTERFACE
      SUBROUTINE SLA_PORFSX_EXTENDED(Prec_type,Uplo,N,Nrhs,A,Lda,Af,    &
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      LOGICAL , INTENT(IN) :: Colequ
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      REAL , DIMENSION(*) :: Res
      REAL , DIMENSION(*) :: Ayb
      REAL , DIMENSION(*) :: Dy
      REAL , DIMENSION(*) :: Y_tail
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL , INTENT(IN) :: Rthresh
      REAL , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
      END SUBROUTINE SLA_PORFSX_EXTENDED
   END INTERFACE
END MODULE S_SLA_PORFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_PORPVGRW
   INTERFACE
      FUNCTION SLA_PORPVGRW(Uplo,Ncols,A,Lda,Af,Ldaf,Work)
      IMPLICIT NONE
      REAL :: SLA_PORPVGRW
      CHARACTER(1) :: Uplo
      INTEGER , INTENT(IN) :: Ncols
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION SLA_PORPVGRW
   END INTERFACE
END MODULE S_SLA_PORPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAPY2
   INTERFACE
      FUNCTION SLAPY2(X,Y)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      REAL :: SLAPY2
      REAL :: X
      REAL :: Y
      END FUNCTION SLAPY2
   END INTERFACE
END MODULE S_SLAPY2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAPY3
   INTERFACE
      FUNCTION SLAPY3(X,Y,Z)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0
      REAL :: SLAPY3
      REAL , INTENT(IN) :: X
      REAL , INTENT(IN) :: Y
      REAL , INTENT(IN) :: Z
      END FUNCTION SLAPY3
   END INTERFACE
END MODULE S_SLAPY3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQGB
   INTERFACE
      SUBROUTINE SLAQGB(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(IN) , DIMENSION(*) :: R
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , INTENT(IN) :: Rowcnd
      REAL , INTENT(IN) :: Colcnd
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE SLAQGB
   END INTERFACE
END MODULE S_SLAQGB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQGE
   INTERFACE
      SUBROUTINE SLAQGE(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: R
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , INTENT(IN) :: Rowcnd
      REAL , INTENT(IN) :: Colcnd
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE SLAQGE
   END INTERFACE
END MODULE S_SLAQGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQP2
   INTERFACE
      SUBROUTINE SLAQP2(M,N,Offset,A,Lda,Jpvt,Tau,Vn1,Vn2,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Offset
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      REAL , DIMENSION(*) :: Tau
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vn1
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vn2
      REAL , DIMENSION(*) :: Work
      END SUBROUTINE SLAQP2
   END INTERFACE
END MODULE S_SLAQP2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQPS
   INTERFACE
      SUBROUTINE SLAQPS(M,N,Offset,Nb,Kb,A,Lda,Jpvt,Tau,Vn1,Vn2,Auxv,F, &
     &                  Ldf)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Offset
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Kb
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      REAL , DIMENSION(*) :: Tau
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vn1
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vn2
      REAL , DIMENSION(*) :: Auxv
      REAL , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      END SUBROUTINE SLAQPS
   END INTERFACE
END MODULE S_SLAQPS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQR0
   INTERFACE
      SUBROUTINE SLAQR0(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,Wr,Wi,Iloz,Ihiz,Z,  &
     &                  Ldz,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NTINY = 15 , KEXNW = 5 , KEXSH = 6
      REAL , PARAMETER  ::  WILK1 = 0.75E0 , WILK2 = -0.4375E0 ,        &
     &                      ZERO = 0.0E0 , ONE = 1.0E0
      LOGICAL :: Wantt
      LOGICAL :: Wantz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL , INTENT(INOUT) , DIMENSION(*) :: Wr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Wi
      INTEGER :: Iloz
      INTEGER :: Ihiz
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER :: Info
      END SUBROUTINE SLAQR0
   END INTERFACE
END MODULE S_SLAQR0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQR1
   INTERFACE
      SUBROUTINE SLAQR1(N,H,Ldh,Sr1,Si1,Sr2,Si2,V)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Ldh,*) :: H
      INTEGER , INTENT(IN) :: Ldh
      REAL , INTENT(IN) :: Sr1
      REAL , INTENT(IN) :: Si1
      REAL , INTENT(IN) :: Sr2
      REAL , INTENT(IN) :: Si2
      REAL , INTENT(OUT) , DIMENSION(*) :: V
      END SUBROUTINE SLAQR1
   END INTERFACE
END MODULE S_SLAQR1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQR2
   INTERFACE
      SUBROUTINE SLAQR2(Wantt,Wantz,N,Ktop,Kbot,Nw,H,Ldh,Iloz,Ihiz,Z,   &
     &                  Ldz,Ns,Nd,Sr,Si,V,Ldv,Nh,T,Ldt,Nv,Wv,Ldwv,Work, &
     &                  Lwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nw
      REAL , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ns
      INTEGER , INTENT(OUT) :: Nd
      REAL , DIMENSION(*) :: Sr
      REAL , DIMENSION(*) :: Si
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(IN) :: Nh
      REAL , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(IN) :: Nv
      REAL , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      END SUBROUTINE SLAQR2
   END INTERFACE
END MODULE S_SLAQR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQR3
   INTERFACE
      SUBROUTINE SLAQR3(Wantt,Wantz,N,Ktop,Kbot,Nw,H,Ldh,Iloz,Ihiz,Z,   &
     &                  Ldz,Ns,Nd,Sr,Si,V,Ldv,Nh,T,Ldt,Nv,Wv,Ldwv,Work, &
     &                  Lwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nw
      REAL , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ns
      INTEGER , INTENT(OUT) :: Nd
      REAL , DIMENSION(*) :: Sr
      REAL , DIMENSION(*) :: Si
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(IN) :: Nh
      REAL , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(IN) :: Nv
      REAL , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      END SUBROUTINE SLAQR3
   END INTERFACE
END MODULE S_SLAQR3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQR4
   INTERFACE
      SUBROUTINE SLAQR4(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,Wr,Wi,Iloz,Ihiz,Z,  &
     &                  Ldz,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NTINY = 15 , KEXNW = 5 , KEXSH = 6
      REAL , PARAMETER  ::  WILK1 = 0.75E0 , WILK2 = -0.4375E0 ,        &
     &                      ZERO = 0.0E0 , ONE = 1.0E0
      LOGICAL :: Wantt
      LOGICAL :: Wantz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL , INTENT(INOUT) , DIMENSION(*) :: Wr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Wi
      INTEGER :: Iloz
      INTEGER :: Ihiz
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER :: Info
      END SUBROUTINE SLAQR4
   END INTERFACE
END MODULE S_SLAQR4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQR5
   INTERFACE
      SUBROUTINE SLAQR5(Wantt,Wantz,Kacc22,N,Ktop,Kbot,Nshfts,Sr,Si,H,  &
     &                  Ldh,Iloz,Ihiz,Z,Ldz,V,Ldv,U,Ldu,Nv,Wv,Ldwv,Nh,  &
     &                  Wh,Ldwh)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: Kacc22
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nshfts
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Si
      REAL , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      REAL , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      REAL , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      INTEGER , INTENT(IN) :: Nv
      REAL , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      INTEGER , INTENT(IN) :: Nh
      REAL , DIMENSION(Ldwh,*) :: Wh
      INTEGER :: Ldwh
      END SUBROUTINE SLAQR5
   END INTERFACE
END MODULE S_SLAQR5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQSB
   INTERFACE
      SUBROUTINE SLAQSB(Uplo,N,Kd,Ab,Ldab,S,Scond,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      REAL , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(IN) , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Scond
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE SLAQSB
   END INTERFACE
END MODULE S_SLAQSB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQSP
   INTERFACE
      SUBROUTINE SLAQSP(Uplo,N,Ap,S,Scond,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL , INTENT(IN) , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Scond
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE SLAQSP
   END INTERFACE
END MODULE S_SLAQSP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQSY
   INTERFACE
      SUBROUTINE SLAQSY(Uplo,N,A,Lda,S,Scond,Amax,Equed)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , THRESH = 0.1E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: S
      REAL , INTENT(IN) :: Scond
      REAL , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE SLAQSY
   END INTERFACE
END MODULE S_SLAQSY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAQTR
   INTERFACE
      SUBROUTINE SLAQTR(Ltran,Lreal,N,T,Ldt,B,W,Scale,X,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      LOGICAL , INTENT(IN) :: Ltran
      LOGICAL , INTENT(IN) :: Lreal
      INTEGER :: N
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(*) :: B
      REAL :: W
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLAQTR
   END INTERFACE
END MODULE S_SLAQTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAR1V
   INTERFACE
      SUBROUTINE SLAR1V(N,B1,Bn,Lambda,D,L,Ld,Lld,Pivmin,Gaptol,Z,      &
     &                  Wantnc,Negcnt,Ztz,Mingma,R,Isuppz,Nrminv,Resid, &
     &                  Rqcorr,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
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
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
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
      END SUBROUTINE SLAR1V
   END INTERFACE
END MODULE S_SLAR1V
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAR2V
   INTERFACE
      SUBROUTINE SLAR2V(N,X,Y,Z,Incx,C,S,Incc)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , INTENT(IN) , DIMENSION(*) :: S
      INTEGER , INTENT(IN) :: Incc
      END SUBROUTINE SLAR2V
   END INTERFACE
END MODULE S_SLAR2V
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARFB
   INTERFACE
      SUBROUTINE SLARFB(Side,Trans,Direct,Storev,M,N,K,V,Ldv,T,Ldt,C,   &
     &                  Ldc,Work,Ldwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE SLARFB
   END INTERFACE
END MODULE S_SLARFB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARFB_GETT
   INTERFACE
      SUBROUTINE SLARFB_GETT(Ident,M,N,K,T,Ldt,A,Lda,B,Ldb,Work,Ldwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Ident
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: K
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE SLARFB_GETT
   END INTERFACE
END MODULE S_SLARFB_GETT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARF
   INTERFACE
      SUBROUTINE SLARF(Side,M,N,V,Incv,Tau,C,Ldc,Work)
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
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      END SUBROUTINE SLARF
   END INTERFACE
END MODULE S_SLARF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARFG
   INTERFACE
      SUBROUTINE SLARFG(N,Alpha,X,Incx,Tau)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) :: Alpha
      REAL , DIMENSION(*) :: X
      INTEGER :: Incx
      REAL , INTENT(OUT) :: Tau
      END SUBROUTINE SLARFG
   END INTERFACE
END MODULE S_SLARFG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARFGP
   INTERFACE
      SUBROUTINE SLARFGP(N,Alpha,X,Incx,Tau)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  TWO = 2.0E+0 , ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) :: Alpha
      REAL , DIMENSION(*) :: X
      INTEGER :: Incx
      REAL , INTENT(INOUT) :: Tau
      END SUBROUTINE SLARFGP
   END INTERFACE
END MODULE S_SLARFGP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARFT
   INTERFACE
      SUBROUTINE SLARFT(Direct,Storev,N,K,V,Ldv,Tau,T,Ldt)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      END SUBROUTINE SLARFT
   END INTERFACE
END MODULE S_SLARFT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARFX
   INTERFACE
      SUBROUTINE SLARFX(Side,M,N,V,Tau,C,Ldc,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(*) :: V
      REAL :: Tau
      REAL , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      END SUBROUTINE SLARFX
   END INTERFACE
END MODULE S_SLARFX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARFY
   INTERFACE
      SUBROUTINE SLARFY(Uplo,N,V,Incv,Tau,C,Ldc,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0 , HALF = 0.5E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(*) :: V
      INTEGER :: Incv
      REAL :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      END SUBROUTINE SLARFY
   END INTERFACE
END MODULE S_SLARFY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARGV
   INTERFACE
      SUBROUTINE SLARGV(N,X,Incx,Y,Incy,C,Incc)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(INOUT) , DIMENSION(*) :: C
      INTEGER , INTENT(IN) :: Incc
      END SUBROUTINE SLARGV
   END INTERFACE
END MODULE S_SLARGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARNV
   INTERFACE
      SUBROUTINE SLARNV(Idist,Iseed,N,X)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , TWO = 2.0E+0
      INTEGER , PARAMETER  ::  LV = 128
      REAL , PARAMETER  ::  TWOPI =                                     &
     &                      6.28318530717958647692528676655900576839E+0
      INTEGER , INTENT(IN) :: Idist
      INTEGER , DIMENSION(4) :: Iseed
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(OUT) , DIMENSION(*) :: X
      END SUBROUTINE SLARNV
   END INTERFACE
END MODULE S_SLARNV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARRA
   INTERFACE
      SUBROUTINE SLARRA(N,D,E,E2,Spltol,Tnrm,Nsplit,Isplit,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , INTENT(OUT) , DIMENSION(*) :: E2
      REAL , INTENT(IN) :: Spltol
      REAL , INTENT(IN) :: Tnrm
      INTEGER , INTENT(INOUT) :: Nsplit
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Isplit
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLARRA
   END INTERFACE
END MODULE S_SLARRA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARRB
   INTERFACE
      SUBROUTINE SLARRB(N,D,Lld,Ifirst,Ilast,Rtol1,Rtol2,Offset,W,Wgap, &
     &                  Werr,Work,Iwork,Pivmin,Spdiam,Twist,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , TWO = 2.0E0 , HALF = 0.5E0
      INTEGER :: N
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: Lld
      INTEGER , INTENT(IN) :: Ifirst
      INTEGER , INTENT(IN) :: Ilast
      REAL , INTENT(IN) :: Rtol1
      REAL , INTENT(IN) :: Rtol2
      INTEGER , INTENT(IN) :: Offset
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , INTENT(INOUT) , DIMENSION(*) :: Wgap
      REAL , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      REAL :: Pivmin
      REAL , INTENT(IN) :: Spdiam
      INTEGER , INTENT(IN) :: Twist
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLARRB
   END INTERFACE
END MODULE S_SLARRB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARRC
   INTERFACE
      SUBROUTINE SLARRC(Jobt,N,Vl,Vu,D,E,Pivmin,Eigcnt,Lcnt,Rcnt,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0
      CHARACTER :: Jobt
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: E
      REAL :: Pivmin
      INTEGER , INTENT(OUT) :: Eigcnt
      INTEGER , INTENT(INOUT) :: Lcnt
      INTEGER , INTENT(INOUT) :: Rcnt
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLARRC
   END INTERFACE
END MODULE S_SLARRC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARRD
   INTERFACE
      SUBROUTINE SLARRD(Range,Order,N,Vl,Vu,Il,Iu,Gers,Reltol,D,E,E2,   &
     &                  Pivmin,Nsplit,Isplit,M,W,Werr,Wl,Wu,Iblock,     &
     &                  Indexw,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      HALF = ONE/TWO , FUDGE = TWO
      INTEGER , PARAMETER  ::  ALLRNG = 1 , VALRNG = 2 , INDRNG = 3
      CHARACTER :: Range
      CHARACTER :: Order
      INTEGER :: N
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      REAL , INTENT(IN) , DIMENSION(*) :: Gers
      REAL , INTENT(IN) :: Reltol
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: E2
      REAL :: Pivmin
      INTEGER , INTENT(IN) :: Nsplit
      INTEGER , INTENT(IN) , DIMENSION(*) :: Isplit
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL , INTENT(INOUT) :: Wl
      REAL , INTENT(INOUT) :: Wu
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indexw
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLARRD
   END INTERFACE
END MODULE S_SLARRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARRE
   INTERFACE
      SUBROUTINE SLARRE(Range,N,Vl,Vu,Il,Iu,D,E,E2,Rtol1,Rtol2,Spltol,  &
     &                  Nsplit,Isplit,M,W,Werr,Wgap,Iblock,Indexw,Gers, &
     &                  Pivmin,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      FOUR = 4.0E0 , HNDRD = 100.0E0 ,            &
     &                      PERT = 4.0E0 , HALF = ONE/TWO ,             &
     &                      FOURTH = ONE/FOUR , FAC = HALF ,            &
     &                      MAXGROWTH = 64.0E0 , FUDGE = 2.0E0
      INTEGER , PARAMETER  ::  MAXTRY = 6 , ALLRNG = 1 , INDRNG = 2 ,   &
     &                         VALRNG = 3
      CHARACTER :: Range
      INTEGER :: N
      REAL , INTENT(INOUT) :: Vl
      REAL , INTENT(INOUT) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: E2
      REAL :: Rtol1
      REAL :: Rtol2
      REAL :: Spltol
      INTEGER :: Nsplit
      INTEGER , DIMENSION(*) :: Isplit
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Wgap
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indexw
      REAL , INTENT(INOUT) , DIMENSION(*) :: Gers
      REAL , INTENT(INOUT) :: Pivmin
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLARRE
   END INTERFACE
END MODULE S_SLARRE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARRF
   INTERFACE
      SUBROUTINE SLARRF(N,D,L,Ld,Clstrt,Clend,W,Wgap,Werr,Spdiam,Clgapl,&
     &                  Clgapr,Pivmin,Sigma,Dplus,Lplus,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , TWO = 2.0E0 , QUART = 0.25E0 ,&
     &                      MAXGROWTH1 = 8.E0 , MAXGROWTH2 = 8.E0
      INTEGER , PARAMETER  ::  KTRYMAX = 1 , SLEFT = 1 , SRIGHT = 2
      INTEGER :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: L
      REAL , INTENT(IN) , DIMENSION(*) :: Ld
      INTEGER , INTENT(IN) :: Clstrt
      INTEGER , INTENT(IN) :: Clend
      REAL , INTENT(IN) , DIMENSION(*) :: W
      REAL , INTENT(IN) , DIMENSION(*) :: Wgap
      REAL , INTENT(IN) , DIMENSION(*) :: Werr
      REAL , INTENT(IN) :: Spdiam
      REAL , INTENT(IN) :: Clgapl
      REAL , INTENT(IN) :: Clgapr
      REAL , INTENT(IN) :: Pivmin
      REAL , INTENT(OUT) :: Sigma
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dplus
      REAL , INTENT(INOUT) , DIMENSION(*) :: Lplus
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLARRF
   END INTERFACE
END MODULE S_SLARRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARRJ
   INTERFACE
      SUBROUTINE SLARRJ(N,D,E2,Ifirst,Ilast,Rtol,Offset,W,Werr,Work,    &
     &                  Iwork,Pivmin,Spdiam,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      HALF = 0.5E0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: E2
      INTEGER , INTENT(IN) :: Ifirst
      INTEGER , INTENT(IN) :: Ilast
      REAL , INTENT(IN) :: Rtol
      INTEGER , INTENT(IN) :: Offset
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      REAL , INTENT(IN) :: Pivmin
      REAL , INTENT(IN) :: Spdiam
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLARRJ
   END INTERFACE
END MODULE S_SLARRJ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARRK
   INTERFACE
      SUBROUTINE SLARRK(N,Iw,Gl,Gu,D,E2,Pivmin,Reltol,W,Werr,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  HALF = 0.5E0 , TWO = 2.0E0 , FUDGE = TWO ,  &
     &                      ZERO = 0.0E0
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Iw
      REAL , INTENT(IN) :: Gl
      REAL , INTENT(IN) :: Gu
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: E2
      REAL , INTENT(IN) :: Pivmin
      REAL , INTENT(IN) :: Reltol
      REAL , INTENT(OUT) :: W
      REAL , INTENT(OUT) :: Werr
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLARRK
   END INTERFACE
END MODULE S_SLARRK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARRR
   INTERFACE
      SUBROUTINE SLARRR(N,D,E,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , RELCOND = 0.999E0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLARRR
   END INTERFACE
END MODULE S_SLARRR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARRV
   INTERFACE
      SUBROUTINE SLARRV(N,Vl,Vu,D,L,Pivmin,Isplit,M,Dol,Dou,Minrgp,     &
     &                  Rtol1,Rtol2,W,Werr,Wgap,Iblock,Indexw,Gers,Z,   &
     &                  Ldz,Isuppz,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXITR = 10
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      THREE = 3.0E0 , FOUR = 4.0E0 , HALF = 0.5E0
      INTEGER :: N
      REAL , INTENT(IN) :: Vl
      REAL :: Vu
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: L
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
      REAL , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isuppz
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLARRV
   END INTERFACE
END MODULE S_SLARRV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARSCL2
   INTERFACE
      SUBROUTINE SLARSCL2(M,N,D,X,Ldx)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      END SUBROUTINE SLARSCL2
   END INTERFACE
END MODULE S_SLARSCL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARTG
   INTERFACE
      SUBROUTINE SLARTG(F,G,Cs,Sn,R)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
      REAL , INTENT(IN) :: F
      REAL , INTENT(IN) :: G
      REAL , INTENT(INOUT) :: Cs
      REAL , INTENT(INOUT) :: Sn
      REAL , INTENT(INOUT) :: R
      END SUBROUTINE SLARTG
   END INTERFACE
END MODULE S_SLARTG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARTGP
   INTERFACE
      SUBROUTINE SLARTGP(F,G,Cs,Sn,R)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
      REAL , INTENT(IN) :: F
      REAL , INTENT(IN) :: G
      REAL , INTENT(INOUT) :: Cs
      REAL , INTENT(INOUT) :: Sn
      REAL , INTENT(INOUT) :: R
      END SUBROUTINE SLARTGP
   END INTERFACE
END MODULE S_SLARTGP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARTGS
   INTERFACE
      SUBROUTINE SLARTGS(X,Y,Sigma,Cs,Sn)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  NEGONE = -1.0E0 , ONE = 1.0E0 , ZERO = 0.0E0
      REAL , INTENT(IN) :: X
      REAL , INTENT(IN) :: Y
      REAL , INTENT(IN) :: Sigma
      REAL :: Cs
      REAL :: Sn
      END SUBROUTINE SLARTGS
   END INTERFACE
END MODULE S_SLARTGS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARTV
   INTERFACE
      SUBROUTINE SLARTV(N,X,Incx,Y,Incy,C,S,Incc)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , INTENT(IN) , DIMENSION(*) :: S
      INTEGER , INTENT(IN) :: Incc
      END SUBROUTINE SLARTV
   END INTERFACE
END MODULE S_SLARTV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARUV
   INTERFACE
      SUBROUTINE SLARUV(Iseed,N,X)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0
      INTEGER , PARAMETER  ::  LV = 128 , IPW2 = 4096
      REAL , PARAMETER  ::  R = ONE/IPW2
      INTEGER , INTENT(INOUT) , DIMENSION(4) :: Iseed
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(N) :: X
      END SUBROUTINE SLARUV
   END INTERFACE
END MODULE S_SLARUV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARZB
   INTERFACE
      SUBROUTINE SLARZB(Side,Trans,Direct,Storev,M,N,K,L,V,Ldv,T,Ldt,C, &
     &                  Ldc,Work,Ldwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: L
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE SLARZB
   END INTERFACE
END MODULE S_SLARZB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARZ
   INTERFACE
      SUBROUTINE SLARZ(Side,M,N,L,V,Incv,Tau,C,Ldc,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      INTEGER :: L
      REAL , DIMENSION(*) :: V
      INTEGER :: Incv
      REAL :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      END SUBROUTINE SLARZ
   END INTERFACE
END MODULE S_SLARZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLARZT
   INTERFACE
      SUBROUTINE SLARZT(Direct,Storev,N,K,V,Ldv,Tau,T,Ldt)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      END SUBROUTINE SLARZT
   END INTERFACE
END MODULE S_SLARZT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAS2
   INTERFACE
      SUBROUTINE SLAS2(F,G,H,Ssmin,Ssmax)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
      REAL , INTENT(IN) :: F
      REAL , INTENT(IN) :: G
      REAL , INTENT(IN) :: H
      REAL , INTENT(INOUT) :: Ssmin
      REAL , INTENT(OUT) :: Ssmax
      END SUBROUTINE SLAS2
   END INTERFACE
END MODULE S_SLAS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASCL2
   INTERFACE
      SUBROUTINE SLASCL2(M,N,D,X,Ldx)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      END SUBROUTINE SLASCL2
   END INTERFACE
END MODULE S_SLASCL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASCL
   INTERFACE
      SUBROUTINE SLASCL(Type,Kl,Ku,Cfrom,Cto,M,N,A,Lda,Info)
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
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASCL
   END INTERFACE
END MODULE S_SLASCL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASD0
   INTERFACE
      SUBROUTINE SLASD0(N,Sqre,D,E,U,Ldu,Vt,Ldvt,Smlsiz,Iwork,Work,Info)
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Sqre
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      INTEGER :: Smlsiz
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASD0
   END INTERFACE
END MODULE S_SLASD0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASD1
   INTERFACE
      SUBROUTINE SLASD1(Nl,Nr,Sqre,D,Alpha,Beta,U,Ldu,Vt,Ldvt,Idxq,     &
     &                  Iwork,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER :: Nl
      INTEGER :: Nr
      INTEGER :: Sqre
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) :: Alpha
      REAL , INTENT(INOUT) :: Beta
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      INTEGER , DIMENSION(*) :: Idxq
      INTEGER , DIMENSION(*) :: Iwork
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASD1
   END INTERFACE
END MODULE S_SLASD1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASD2
   INTERFACE
      SUBROUTINE SLASD2(Nl,Nr,Sqre,K,D,Z,Alpha,Beta,U,Ldu,Vt,Ldvt,      &
     &                  Dsigma,U2,Ldu2,Vt2,Ldvt2,Idxp,Idx,Idxc,Idxq,    &
     &                  Coltyp,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , EIGHT = 8.0E+0
      INTEGER :: Nl
      INTEGER :: Nr
      INTEGER , INTENT(IN) :: Sqre
      INTEGER , INTENT(INOUT) :: K
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) :: Beta
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , INTENT(INOUT) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dsigma
      REAL , INTENT(INOUT) , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      REAL , DIMENSION(Ldvt2,*) :: Vt2
      INTEGER :: Ldvt2
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxp
      INTEGER , DIMENSION(*) :: Idx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxc
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxq
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Coltyp
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASD2
   END INTERFACE
END MODULE S_SLASD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASD3
   INTERFACE
      SUBROUTINE SLASD3(Nl,Nr,Sqre,K,D,Q,Ldq,Dsigma,U,Ldu,U2,Ldu2,Vt,   &
     &                  Ldvt,Vt2,Ldvt2,Idxc,Ctot,Z,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0 ,              &
     &                      NEGONE = -1.0E+0
      INTEGER :: Nl
      INTEGER :: Nr
      INTEGER , INTENT(IN) :: Sqre
      INTEGER :: K
      REAL , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dsigma
      REAL , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      REAL , INTENT(INOUT) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL , INTENT(INOUT) , DIMENSION(Ldvt2,*) :: Vt2
      INTEGER :: Ldvt2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Idxc
      INTEGER , DIMENSION(*) :: Ctot
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASD3
   END INTERFACE
END MODULE S_SLASD3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASD4
   INTERFACE
      SUBROUTINE SLASD4(N,I,D,Z,Delta,Rho,Sigma,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXIT = 400
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , THREE = 3.0E+0 ,             &
     &                      FOUR = 4.0E+0 , EIGHT = 8.0E+0 ,            &
     &                      TEN = 10.0E+0
      INTEGER , INTENT(IN) :: N
      INTEGER :: I
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: Z
      REAL , INTENT(INOUT) , DIMENSION(*) :: Delta
      REAL :: Rho
      REAL , INTENT(INOUT) :: Sigma
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASD4
   END INTERFACE
END MODULE S_SLASD4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASD5
   INTERFACE
      SUBROUTINE SLASD5(I,D,Z,Delta,Rho,Dsigma,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , THREE = 3.0E+0 ,             &
     &                      FOUR = 4.0E+0
      INTEGER , INTENT(IN) :: I
      REAL , INTENT(IN) , DIMENSION(2) :: D
      REAL , INTENT(IN) , DIMENSION(2) :: Z
      REAL , INTENT(OUT) , DIMENSION(2) :: Delta
      REAL , INTENT(IN) :: Rho
      REAL , INTENT(OUT) :: Dsigma
      REAL , INTENT(OUT) , DIMENSION(2) :: Work
      END SUBROUTINE SLASD5
   END INTERFACE
END MODULE S_SLASD5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASD6
   INTERFACE
      SUBROUTINE SLASD6(Icompq,Nl,Nr,Sqre,D,Vf,Vl,Alpha,Beta,Idxq,Perm, &
     &                  Givptr,Givcol,Ldgcol,Givnum,Ldgnum,Poles,Difl,  &
     &                  Difr,Z,K,C,S,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER :: Icompq
      INTEGER :: Nl
      INTEGER :: Nr
      INTEGER :: Sqre
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: Vf
      REAL , DIMENSION(*) :: Vl
      REAL , INTENT(INOUT) :: Alpha
      REAL , INTENT(INOUT) :: Beta
      INTEGER , DIMENSION(*) :: Idxq
      INTEGER , DIMENSION(*) :: Perm
      INTEGER :: Givptr
      INTEGER , DIMENSION(Ldgcol,*) :: Givcol
      INTEGER :: Ldgcol
      REAL , DIMENSION(Ldgnum,*) :: Givnum
      INTEGER :: Ldgnum
      REAL , DIMENSION(Ldgnum,*) :: Poles
      REAL , DIMENSION(*) :: Difl
      REAL , DIMENSION(*) :: Difr
      REAL , DIMENSION(*) :: Z
      INTEGER :: K
      REAL :: C
      REAL :: S
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASD6
   END INTERFACE
END MODULE S_SLASD6
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASD7
   INTERFACE
      SUBROUTINE SLASD7(Icompq,Nl,Nr,Sqre,K,D,Z,Zw,Vf,Vfw,Vl,Vlw,Alpha, &
     &                  Beta,Dsigma,Idx,Idxp,Idxq,Perm,Givptr,Givcol,   &
     &                  Ldgcol,Givnum,Ldgnum,C,S,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , EIGHT = 8.0E+0
      INTEGER , INTENT(IN) :: Icompq
      INTEGER :: Nl
      INTEGER :: Nr
      INTEGER , INTENT(IN) :: Sqre
      INTEGER , INTENT(INOUT) :: K
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL , INTENT(INOUT) , DIMENSION(*) :: Zw
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vf
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vfw
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vl
      REAL , INTENT(INOUT) , DIMENSION(*) :: Vlw
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dsigma
      INTEGER , DIMENSION(*) :: Idx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxp
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxq
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Perm
      INTEGER , INTENT(INOUT) :: Givptr
      INTEGER , INTENT(OUT) , DIMENSION(Ldgcol,*) :: Givcol
      INTEGER , INTENT(IN) :: Ldgcol
      REAL , INTENT(OUT) , DIMENSION(Ldgnum,*) :: Givnum
      INTEGER , INTENT(IN) :: Ldgnum
      REAL , INTENT(INOUT) :: C
      REAL , INTENT(INOUT) :: S
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASD7
   END INTERFACE
END MODULE S_SLASD7
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASD8
   INTERFACE
      SUBROUTINE SLASD8(Icompq,K,D,Z,Vf,Vl,Difl,Difr,Lddifr,Dsigma,Work,&
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER , INTENT(IN) :: Icompq
      INTEGER :: K
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL , DIMENSION(*) :: Vf
      REAL , DIMENSION(*) :: Vl
      REAL , INTENT(INOUT) , DIMENSION(*) :: Difl
      REAL , INTENT(INOUT) , DIMENSION(Lddifr,*) :: Difr
      INTEGER , INTENT(IN) :: Lddifr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dsigma
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASD8
   END INTERFACE
END MODULE S_SLASD8
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASDA
   INTERFACE
      SUBROUTINE SLASDA(Icompq,Smlsiz,N,Sqre,D,E,U,Ldu,Vt,K,Difl,Difr,Z,&
     &                  Poles,Givptr,Givcol,Ldgcol,Perm,Givnum,C,S,Work,&
     &                  Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER :: Icompq
      INTEGER :: Smlsiz
      INTEGER :: N
      INTEGER :: Sqre
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
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
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASDA
   END INTERFACE
END MODULE S_SLASDA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASDQ
   INTERFACE
      SUBROUTINE SLASDQ(Uplo,Sqre,N,Ncvt,Nru,Ncc,D,E,Vt,Ldvt,U,Ldu,C,   &
     &                  Ldc,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: Sqre
      INTEGER :: N
      INTEGER :: Ncvt
      INTEGER :: Nru
      INTEGER :: Ncc
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASDQ
   END INTERFACE
END MODULE S_SLASDQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASDT
   INTERFACE
      SUBROUTINE SLASDT(N,Lvl,Nd,Inode,Ndiml,Ndimr,Msub)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  TWO = 2.0E+0
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(INOUT) :: Lvl
      INTEGER , INTENT(OUT) :: Nd
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Inode
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ndiml
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ndimr
      INTEGER , INTENT(IN) :: Msub
      END SUBROUTINE SLASDT
   END INTERFACE
END MODULE S_SLASDT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASET
   INTERFACE
      SUBROUTINE SLASET(Uplo,M,N,Alpha,Beta,A,Lda)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(OUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE SLASET
   END INTERFACE
END MODULE S_SLASET
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASQ1
   INTERFACE
      SUBROUTINE SLASQ1(N,D,E,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASQ1
   END INTERFACE
END MODULE S_SLASQ1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASQ2
   INTERFACE
      SUBROUTINE SLASQ2(N,Z,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  CBIAS = 1.50E0 , ZERO = 0.0E0 ,             &
     &                      HALF = 0.5E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      FOUR = 4.0E0 , HUNDRD = 100.0E0
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLASQ2
   END INTERFACE
END MODULE S_SLASQ2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASQ3
   INTERFACE
      SUBROUTINE SLASQ3(I0,N0,Z,Pp,Dmin,Sigma,Desig,Qmax,Nfail,Iter,    &
     &                  Ndiv,Ieee,Ttype,Dmin1,Dmin2,Dn,Dn1,Dn2,G,Tau)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  CBIAS = 1.50E0 , ZERO = 0.0E0 ,             &
     &                      QURTR = 0.250E0 , HALF = 0.5E0 ,            &
     &                      ONE = 1.0E0 , TWO = 2.0E0 , HUNDRD = 100.0E0
      INTEGER :: I0
      INTEGER , INTENT(INOUT) :: N0
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(INOUT) :: Pp
      REAL , INTENT(INOUT) :: Dmin
      REAL , INTENT(INOUT) :: Sigma
      REAL , INTENT(INOUT) :: Desig
      REAL , INTENT(INOUT) :: Qmax
      INTEGER , INTENT(INOUT) :: Nfail
      INTEGER , INTENT(INOUT) :: Iter
      INTEGER , INTENT(INOUT) :: Ndiv
      LOGICAL :: Ieee
      INTEGER , INTENT(INOUT) :: Ttype
      REAL :: Dmin1
      REAL , INTENT(INOUT) :: Dmin2
      REAL :: Dn
      REAL :: Dn1
      REAL :: Dn2
      REAL :: G
      REAL , INTENT(INOUT) :: Tau
      END SUBROUTINE SLASQ3
   END INTERFACE
END MODULE S_SLASQ3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASQ4
   INTERFACE
      SUBROUTINE SLASQ4(I0,N0,Z,Pp,N0in,Dmin,Dmin1,Dmin2,Dn,Dn1,Dn2,Tau,&
     &                  Ttype,G)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  CNST1 = 0.5630E0 , CNST2 = 1.010E0 ,        &
     &                      CNST3 = 1.050E0 , QURTR = 0.250E0 ,         &
     &                      THIRD = 0.3330E0 , HALF = 0.50E0 ,          &
     &                      ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      HUNDRD = 100.0E0
      INTEGER , INTENT(IN) :: I0
      INTEGER , INTENT(IN) :: N0
      REAL , INTENT(IN) , DIMENSION(*) :: Z
      INTEGER , INTENT(IN) :: Pp
      INTEGER , INTENT(IN) :: N0in
      REAL , INTENT(IN) :: Dmin
      REAL , INTENT(IN) :: Dmin1
      REAL , INTENT(IN) :: Dmin2
      REAL , INTENT(IN) :: Dn
      REAL , INTENT(IN) :: Dn1
      REAL , INTENT(IN) :: Dn2
      REAL , INTENT(OUT) :: Tau
      INTEGER , INTENT(INOUT) :: Ttype
      REAL , INTENT(INOUT) :: G
      END SUBROUTINE SLASQ4
   END INTERFACE
END MODULE S_SLASQ4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASQ5
   INTERFACE
      SUBROUTINE SLASQ5(I0,N0,Z,Pp,Tau,Sigma,Dmin,Dmin1,Dmin2,Dn,Dnm1,  &
     &                  Dnm2,Ieee,Eps)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , HALF = 0.5
      INTEGER , INTENT(IN) :: I0
      INTEGER , INTENT(IN) :: N0
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(IN) :: Pp
      REAL , INTENT(INOUT) :: Tau
      REAL , INTENT(IN) :: Sigma
      REAL , INTENT(INOUT) :: Dmin
      REAL , INTENT(OUT) :: Dmin1
      REAL , INTENT(OUT) :: Dmin2
      REAL , INTENT(INOUT) :: Dn
      REAL , INTENT(INOUT) :: Dnm1
      REAL , INTENT(INOUT) :: Dnm2
      LOGICAL , INTENT(IN) :: Ieee
      REAL , INTENT(IN) :: Eps
      END SUBROUTINE SLASQ5
   END INTERFACE
END MODULE S_SLASQ5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASQ6
   INTERFACE
      SUBROUTINE SLASQ6(I0,N0,Z,Pp,Dmin,Dmin1,Dmin2,Dn,Dnm1,Dnm2)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0
      INTEGER , INTENT(IN) :: I0
      INTEGER , INTENT(IN) :: N0
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(IN) :: Pp
      REAL , INTENT(INOUT) :: Dmin
      REAL , INTENT(OUT) :: Dmin1
      REAL , INTENT(OUT) :: Dmin2
      REAL , INTENT(INOUT) :: Dn
      REAL , INTENT(INOUT) :: Dnm1
      REAL , INTENT(INOUT) :: Dnm2
      END SUBROUTINE SLASQ6
   END INTERFACE
END MODULE S_SLASQ6
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASR
   INTERFACE
      SUBROUTINE SLASR(Side,Pivot,Direct,M,N,C,S,A,Lda)
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
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE SLASR
   END INTERFACE
END MODULE S_SLASR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASRT
   INTERFACE
      SUBROUTINE SLASRT(Id,N,D,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  SELECT = 20
      CHARACTER :: Id
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASRT
   END INTERFACE
END MODULE S_SLASRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASSQ
   INTERFACE
      SUBROUTINE SLASSQ(N,X,Incx,Scale,Sumsq)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(INOUT) :: Sumsq
      END SUBROUTINE SLASSQ
   END INTERFACE
END MODULE S_SLASSQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASV2
   INTERFACE
      SUBROUTINE SLASV2(F,G,H,Ssmin,Ssmax,Snr,Csr,Snl,Csl)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , HALF = 0.5E0 , ONE = 1.0E0 , &
     &                      TWO = 2.0E0 , FOUR = 4.0E0
      REAL , INTENT(IN) :: F
      REAL , INTENT(IN) :: G
      REAL , INTENT(IN) :: H
      REAL , INTENT(INOUT) :: Ssmin
      REAL , INTENT(INOUT) :: Ssmax
      REAL , INTENT(INOUT) :: Snr
      REAL , INTENT(INOUT) :: Csr
      REAL , INTENT(INOUT) :: Snl
      REAL , INTENT(INOUT) :: Csl
      END SUBROUTINE SLASV2
   END INTERFACE
END MODULE S_SLASV2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASWLQ
   INTERFACE
      SUBROUTINE SLASWLQ(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER :: Nb
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASWLQ
   END INTERFACE
END MODULE S_SLASWLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASWP
   INTERFACE
      SUBROUTINE SLASWP(N,A,Lda,K1,K2,Ipiv,Incx)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) :: K1
      INTEGER , INTENT(IN) :: K2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE SLASWP
   END INTERFACE
END MODULE S_SLASWP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASY2
   INTERFACE
      SUBROUTINE SLASY2(Ltranl,Ltranr,Isgn,N1,N2,Tl,Ldtl,Tr,Ldtr,B,Ldb, &
     &                  Scale,X,Ldx,Xnorm,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , HALF = 0.5E+0 ,              &
     &                      EIGHT = 8.0E+0
      LOGICAL , INTENT(IN) :: Ltranl
      LOGICAL , INTENT(IN) :: Ltranr
      INTEGER , INTENT(IN) :: Isgn
      INTEGER , INTENT(IN) :: N1
      INTEGER , INTENT(IN) :: N2
      REAL , INTENT(IN) , DIMENSION(Ldtl,*) :: Tl
      INTEGER , INTENT(IN) :: Ldtl
      REAL , INTENT(IN) , DIMENSION(Ldtr,*) :: Tr
      INTEGER , INTENT(IN) :: Ldtr
      REAL , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(OUT) :: Xnorm
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE SLASY2
   END INTERFACE
END MODULE S_SLASY2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_SYAMV
   INTERFACE
      SUBROUTINE SLA_SYAMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE SLA_SYAMV
   END INTERFACE
END MODULE S_SLA_SYAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASYF_AA
   INTERFACE
      SUBROUTINE SLASYF_AA(Uplo,J1,M,Nb,A,Lda,Ipiv,H,Ldh,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: J1
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: Nb
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END SUBROUTINE SLASYF_AA
   END INTERFACE
END MODULE S_SLASYF_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASYF
   INTERFACE
      SUBROUTINE SLASYF(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASYF
   END INTERFACE
END MODULE S_SLASYF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASYF_RK
   INTERFACE
      SUBROUTINE SLASYF_RK(Uplo,N,Nb,Kb,A,Lda,E,Ipiv,W,Ldw,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASYF_RK
   END INTERFACE
END MODULE S_SLASYF_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLASYF_ROOK
   INTERFACE
      SUBROUTINE SLASYF_ROOK(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLASYF_ROOK
   END INTERFACE
END MODULE S_SLASYF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_SYRCOND
   INTERFACE
      FUNCTION SLA_SYRCOND(Uplo,N,A,Lda,Af,Ldaf,Ipiv,Cmode,C,Info,Work, &
     &                     Iwork)
      IMPLICIT NONE
      REAL :: SLA_SYRCOND
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(IN) :: Cmode
      REAL , INTENT(IN) , DIMENSION(*) :: C
      INTEGER , INTENT(INOUT) :: Info
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      END FUNCTION SLA_SYRCOND
   END INTERFACE
END MODULE S_SLA_SYRCOND
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_SYRFSX_EXTENDED
   INTERFACE
      SUBROUTINE SLA_SYRFSX_EXTENDED(Prec_type,Uplo,N,Nrhs,A,Lda,Af,    &
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL , INTENT(IN) , DIMENSION(*) :: C
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      REAL , DIMENSION(*) :: Res
      REAL , DIMENSION(*) :: Ayb
      REAL , DIMENSION(*) :: Dy
      REAL , DIMENSION(*) :: Y_tail
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL , INTENT(IN) :: Rthresh
      REAL , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLA_SYRFSX_EXTENDED
   END INTERFACE
END MODULE S_SLA_SYRFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_SYRPVGRW
   INTERFACE
      FUNCTION SLA_SYRPVGRW(Uplo,N,Info,A,Lda,Af,Ldaf,Ipiv,Work)
      IMPLICIT NONE
      REAL :: SLA_SYRPVGRW
      CHARACTER(1) :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Info
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION SLA_SYRPVGRW
   END INTERFACE
END MODULE S_SLA_SYRPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLATBS
   INTERFACE
      SUBROUTINE SLATBS(Uplo,Trans,Diag,Normin,N,Kd,Ab,Ldab,X,Scale,    &
     &                  Cnorm,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , HALF = 0.5E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      INTEGER :: Kd
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLATBS
   END INTERFACE
END MODULE S_SLATBS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLATDF
   INTERFACE
      SUBROUTINE SLATDF(Ijob,N,Z,Ldz,Rhs,Rdsum,Rdscal,Ipiv,Jpiv)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXDIM = 8
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , INTENT(IN) :: Ijob
      INTEGER :: N
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rhs
      REAL :: Rdsum
      REAL :: Rdscal
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Jpiv
      END SUBROUTINE SLATDF
   END INTERFACE
END MODULE S_SLATDF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLATPS
   INTERFACE
      SUBROUTINE SLATPS(Uplo,Trans,Diag,Normin,N,Ap,X,Scale,Cnorm,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , HALF = 0.5E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      REAL , DIMENSION(*) :: Ap
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLATPS
   END INTERFACE
END MODULE S_SLATPS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLATRD
   INTERFACE
      SUBROUTINE SLATRD(Uplo,N,Nb,A,Lda,E,Tau,W,Ldw)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , HALF = 0.5E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      END SUBROUTINE SLATRD
   END INTERFACE
END MODULE S_SLATRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLATRS
   INTERFACE
      SUBROUTINE SLATRS(Uplo,Trans,Diag,Normin,N,A,Lda,X,Scale,Cnorm,   &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , HALF = 0.5E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLATRS
   END INTERFACE
END MODULE S_SLATRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLATRZ
   INTERFACE
      SUBROUTINE SLATRZ(M,N,L,A,Lda,Tau,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: L
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      END SUBROUTINE SLATRZ
   END INTERFACE
END MODULE S_SLATRZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLATSQR
   INTERFACE
      SUBROUTINE SLATSQR(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER :: Nb
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLATSQR
   END INTERFACE
END MODULE S_SLATSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAUU2
   INTERFACE
      SUBROUTINE SLAUU2(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SLAUU2
   END INTERFACE
END MODULE S_SLAUU2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLAUUM
   INTERFACE
      SUBROUTINE SLAUUM(Uplo,N,A,Lda,Info)
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
      END SUBROUTINE SLAUUM
   END INTERFACE
END MODULE S_SLAUUM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SLA_WWADDW
   INTERFACE
      SUBROUTINE SLA_WWADDW(N,X,Y,W)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: X
      REAL , INTENT(INOUT) , DIMENSION(*) :: Y
      REAL , INTENT(IN) , DIMENSION(*) :: W
      END SUBROUTINE SLA_WWADDW
   END INTERFACE
END MODULE S_SLA_WWADDW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SOPGTR
   INTERFACE
      SUBROUTINE SOPGTR(Uplo,N,Ap,Tau,Q,Ldq,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SOPGTR
   END INTERFACE
END MODULE S_SOPGTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SOPMTR
   INTERFACE
      SUBROUTINE SOPMTR(Side,Uplo,Trans,M,N,Ap,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SOPMTR
   END INTERFACE
END MODULE S_SOPMTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORBDB1
   INTERFACE
      SUBROUTINE SORBDB1(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      REAL , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL , INTENT(OUT) , DIMENSION(*) :: Phi
      REAL , DIMENSION(*) :: Taup1
      REAL , DIMENSION(*) :: Taup2
      REAL , DIMENSION(*) :: Tauq1
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORBDB1
   END INTERFACE
END MODULE S_SORBDB1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORBDB2
   INTERFACE
      SUBROUTINE SORBDB2(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  NEGONE = -1.0E0 , ONE = 1.0E0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      REAL , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , INTENT(OUT) , DIMENSION(*) :: Theta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Phi
      REAL , DIMENSION(*) :: Taup1
      REAL , DIMENSION(*) :: Taup2
      REAL , DIMENSION(*) :: Tauq1
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORBDB2
   END INTERFACE
END MODULE S_SORBDB2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORBDB3
   INTERFACE
      SUBROUTINE SORBDB3(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      REAL , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , INTENT(OUT) , DIMENSION(*) :: Theta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Phi
      REAL , DIMENSION(*) :: Taup1
      REAL , DIMENSION(*) :: Taup2
      REAL , DIMENSION(*) :: Tauq1
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORBDB3
   END INTERFACE
END MODULE S_SORBDB3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORBDB4
   INTERFACE
      SUBROUTINE SORBDB4(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Phantom,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  NEGONE = -1.0E0 , ONE = 1.0E0 , ZERO = 0.0E0
      INTEGER , INTENT(IN) :: M
      INTEGER :: P
      INTEGER :: Q
      REAL , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL , INTENT(OUT) , DIMENSION(*) :: Phi
      REAL , DIMENSION(*) :: Taup1
      REAL , DIMENSION(*) :: Taup2
      REAL , DIMENSION(*) :: Tauq1
      REAL , INTENT(INOUT) , DIMENSION(*) :: Phantom
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORBDB4
   END INTERFACE
END MODULE S_SORBDB4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORBDB5
   INTERFACE
      SUBROUTINE SORBDB5(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,&
     &                   Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , ZERO = 0.0E0
      INTEGER :: M1
      INTEGER :: M2
      INTEGER :: N
      REAL , DIMENSION(*) :: X1
      INTEGER :: Incx1
      REAL , DIMENSION(*) :: X2
      INTEGER :: Incx2
      REAL , DIMENSION(Ldq1,*) :: Q1
      INTEGER :: Ldq1
      REAL , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORBDB5
   END INTERFACE
END MODULE S_SORBDB5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORBDB6
   INTERFACE
      SUBROUTINE SORBDB6(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,&
     &                   Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ALPHASQ = 0.01E0 , REALONE = 1.0E0 ,        &
     &                      REALZERO = 0.0E0 , NEGONE = -1.0E0 ,        &
     &                      ONE = 1.0E0 , ZERO = 0.0E0
      INTEGER :: M1
      INTEGER :: M2
      INTEGER :: N
      REAL , DIMENSION(*) :: X1
      INTEGER :: Incx1
      REAL , DIMENSION(*) :: X2
      INTEGER :: Incx2
      REAL , DIMENSION(Ldq1,*) :: Q1
      INTEGER :: Ldq1
      REAL , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORBDB6
   END INTERFACE
END MODULE S_SORBDB6
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORBDB
   INTERFACE
      SUBROUTINE SORBDB(Trans,Signs,M,P,Q,X11,Ldx11,X12,Ldx12,X21,Ldx21,&
     &                  X22,Ldx22,Theta,Phi,Taup1,Taup2,Tauq1,Tauq2,    &
     &                  Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  REALONE = 1.0E0 , ONE = 1.0E0
      CHARACTER :: Trans
      CHARACTER :: Signs
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      REAL , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL , DIMENSION(Ldx12,*) :: X12
      INTEGER :: Ldx12
      REAL , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , DIMENSION(Ldx22,*) :: X22
      INTEGER :: Ldx22
      REAL , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL , INTENT(INOUT) , DIMENSION(*) :: Phi
      REAL , DIMENSION(*) :: Taup1
      REAL , DIMENSION(*) :: Taup2
      REAL , DIMENSION(*) :: Tauq1
      REAL , DIMENSION(*) :: Tauq2
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORBDB
   END INTERFACE
END MODULE S_SORBDB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORCSD2BY1
   INTERFACE
      SUBROUTINE SORCSD2BY1(Jobu1,Jobu2,Jobv1t,M,P,Q,X11,Ldx11,X21,     &
     &                      Ldx21,Theta,U1,Ldu1,U2,Ldu2,V1t,Ldv1t,Work, &
     &                      Lwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , ZERO = 0.0E0
      CHARACTER :: Jobu1
      CHARACTER :: Jobu2
      CHARACTER :: Jobv1t
      INTEGER :: M
      INTEGER :: P
      INTEGER :: Q
      REAL , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , DIMENSION(*) :: Theta
      REAL , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      REAL , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      REAL , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORCSD2BY1
   END INTERFACE
END MODULE S_SORCSD2BY1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORCSD
   INTERFACE
      RECURSIVE SUBROUTINE SORCSD(Jobu1,Jobu2,Jobv1t,Jobv2t,Trans,Signs,&
     &                            M,P,Q,X11,Ldx11,X12,Ldx12,X21,Ldx21,  &
     &                            X22,Ldx22,Theta,U1,Ldu1,U2,Ldu2,V1t,  &
     &                            Ldv1t,V2t,Ldv2t,Work,Lwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Jobu1
      CHARACTER :: Jobu2
      CHARACTER :: Jobv1t
      CHARACTER :: Jobv2t
      CHARACTER :: Trans
      CHARACTER :: Signs
      INTEGER :: M
      INTEGER :: P
      INTEGER :: Q
      REAL , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL , DIMENSION(Ldx12,*) :: X12
      INTEGER :: Ldx12
      REAL , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , DIMENSION(Ldx22,*) :: X22
      INTEGER :: Ldx22
      REAL , DIMENSION(*) :: Theta
      REAL , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      REAL , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      REAL , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      REAL , DIMENSION(Ldv2t,*) :: V2t
      INTEGER :: Ldv2t
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORCSD
   END INTERFACE
END MODULE S_SORCSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORG2L
   INTERFACE
      SUBROUTINE SORG2L(M,N,K,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORG2L
   END INTERFACE
END MODULE S_SORG2L
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORG2R
   INTERFACE
      SUBROUTINE SORG2R(M,N,K,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORG2R
   END INTERFACE
END MODULE S_SORG2R
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORGBR
   INTERFACE
      SUBROUTINE SORGBR(Vect,M,N,K,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Vect
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORGBR
   END INTERFACE
END MODULE S_SORGBR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORGHR
   INTERFACE
      SUBROUTINE SORGHR(N,Ilo,Ihi,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORGHR
   END INTERFACE
END MODULE S_SORGHR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORGL2
   INTERFACE
      SUBROUTINE SORGL2(M,N,K,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORGL2
   END INTERFACE
END MODULE S_SORGL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORGLQ
   INTERFACE
      SUBROUTINE SORGLQ(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORGLQ
   END INTERFACE
END MODULE S_SORGLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORGQL
   INTERFACE
      SUBROUTINE SORGQL(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORGQL
   END INTERFACE
END MODULE S_SORGQL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORGQR
   INTERFACE
      SUBROUTINE SORGQR(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORGQR
   END INTERFACE
END MODULE S_SORGQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORGR2
   INTERFACE
      SUBROUTINE SORGR2(M,N,K,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORGR2
   END INTERFACE
END MODULE S_SORGR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORGRQ
   INTERFACE
      SUBROUTINE SORGRQ(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORGRQ
   END INTERFACE
END MODULE S_SORGRQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORGTR
   INTERFACE
      SUBROUTINE SORGTR(Uplo,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORGTR
   END INTERFACE
END MODULE S_SORGTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORGTSQR
   INTERFACE
      SUBROUTINE SORGTSQR(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER , INTENT(IN) :: Nb
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORGTSQR
   END INTERFACE
END MODULE S_SORGTSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORGTSQR_ROW
   INTERFACE
      SUBROUTINE SORGTSQR_ROW(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Mb
      INTEGER , INTENT(IN) :: Nb
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORGTSQR_ROW
   END INTERFACE
END MODULE S_SORGTSQR_ROW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORHR_COL
   INTERFACE
      SUBROUTINE SORHR_COL(M,N,Nb,A,Lda,T,Ldt,D,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nb
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORHR_COL
   END INTERFACE
END MODULE S_SORHR_COL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORM22
   INTERFACE
      SUBROUTINE SORM22(Side,Trans,M,N,N1,N2,Q,Ldq,C,Ldc,Work,Lwork,    &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: N1
      INTEGER :: N2
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORM22
   END INTERFACE
END MODULE S_SORM22
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORM2L
   INTERFACE
      SUBROUTINE SORM2L(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORM2L
   END INTERFACE
END MODULE S_SORM2L
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORM2R
   INTERFACE
      SUBROUTINE SORM2R(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORM2R
   END INTERFACE
END MODULE S_SORM2R
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORMBR
   INTERFACE
      SUBROUTINE SORMBR(Vect,Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,     &
     &                  Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Vect
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORMBR
   END INTERFACE
END MODULE S_SORMBR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORMHR
   INTERFACE
      SUBROUTINE SORMHR(Side,Trans,M,N,Ilo,Ihi,A,Lda,Tau,C,Ldc,Work,    &
     &                  Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORMHR
   END INTERFACE
END MODULE S_SORMHR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORML2
   INTERFACE
      SUBROUTINE SORML2(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORML2
   END INTERFACE
END MODULE S_SORML2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORMLQ
   INTERFACE
      SUBROUTINE SORMLQ(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORMLQ
   END INTERFACE
END MODULE S_SORMLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORMQL
   INTERFACE
      SUBROUTINE SORMQL(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORMQL
   END INTERFACE
END MODULE S_SORMQL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORMQR
   INTERFACE
      SUBROUTINE SORMQR(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORMQR
   END INTERFACE
END MODULE S_SORMQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORMR2
   INTERFACE
      SUBROUTINE SORMR2(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORMR2
   END INTERFACE
END MODULE S_SORMR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORMR3
   INTERFACE
      SUBROUTINE SORMR3(Side,Trans,M,N,K,L,A,Lda,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      INTEGER :: L
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORMR3
   END INTERFACE
END MODULE S_SORMR3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORMRQ
   INTERFACE
      SUBROUTINE SORMRQ(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORMRQ
   END INTERFACE
END MODULE S_SORMRQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORMRZ
   INTERFACE
      SUBROUTINE SORMRZ(Side,Trans,M,N,K,L,A,Lda,Tau,C,Ldc,Work,Lwork,  &
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORMRZ
   END INTERFACE
END MODULE S_SORMRZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SORMTR
   INTERFACE
      SUBROUTINE SORMTR(Side,Uplo,Trans,M,N,A,Lda,Tau,C,Ldc,Work,Lwork, &
     &                  Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SORMTR
   END INTERFACE
END MODULE S_SORMTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPBCON
   INTERFACE
      SUBROUTINE SPBCON(Uplo,N,Kd,Ab,Ldab,Anorm,Rcond,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPBCON
   END INTERFACE
END MODULE S_SPBCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPBEQU
   INTERFACE
      SUBROUTINE SPBEQU(Uplo,N,Kd,Ab,Ldab,S,Scond,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      REAL , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) :: Scond
      REAL , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPBEQU
   END INTERFACE
END MODULE S_SPBEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPBRFS
   INTERFACE
      SUBROUTINE SPBRFS(Uplo,N,Kd,Nrhs,Ab,Ldab,Afb,Ldafb,B,Ldb,X,Ldx,   &
     &                  Ferr,Berr,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , THREE = 3.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPBRFS
   END INTERFACE
END MODULE S_SPBRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPBSTF
   INTERFACE
      SUBROUTINE SPBSTF(Uplo,N,Kd,Ab,Ldab,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      REAL , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPBSTF
   END INTERFACE
END MODULE S_SPBSTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPBSV
   INTERFACE
      SUBROUTINE SPBSV(Uplo,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER :: Nrhs
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPBSV
   END INTERFACE
END MODULE S_SPBSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPBSVX
   INTERFACE
      SUBROUTINE SPBSVX(Fact,Uplo,N,Kd,Nrhs,Ab,Ldab,Afb,Ldafb,Equed,S,B,&
     &                  Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Iwork,Info)
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
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: S
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPBSVX
   END INTERFACE
END MODULE S_SPBSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPBTF2
   INTERFACE
      SUBROUTINE SPBTF2(Uplo,N,Kd,Ab,Ldab,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      REAL , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPBTF2
   END INTERFACE
END MODULE S_SPBTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPBTRF
   INTERFACE
      SUBROUTINE SPBTRF(Uplo,N,Kd,Ab,Ldab,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER , PARAMETER  ::  NBMAX = 32 , LDWORK = NBMAX + 1
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPBTRF
   END INTERFACE
END MODULE S_SPBTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPBTRS
   INTERFACE
      SUBROUTINE SPBTRS(Uplo,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPBTRS
   END INTERFACE
END MODULE S_SPBTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPFTRF
   INTERFACE
      SUBROUTINE SPFTRF(Transr,Uplo,N,A,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPFTRF
   END INTERFACE
END MODULE S_SPFTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPFTRI
   INTERFACE
      SUBROUTINE SPFTRI(Transr,Uplo,N,A,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPFTRI
   END INTERFACE
END MODULE S_SPFTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPFTRS
   INTERFACE
      SUBROUTINE SPFTRS(Transr,Uplo,N,Nrhs,A,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(0:*) :: A
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPFTRS
   END INTERFACE
END MODULE S_SPFTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPOCON
   INTERFACE
      SUBROUTINE SPOCON(Uplo,N,A,Lda,Anorm,Rcond,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPOCON
   END INTERFACE
END MODULE S_SPOCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPOEQUB
   INTERFACE
      SUBROUTINE SPOEQUB(N,A,Lda,S,Scond,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) :: Scond
      REAL , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPOEQUB
   END INTERFACE
END MODULE S_SPOEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPOEQU
   INTERFACE
      SUBROUTINE SPOEQU(N,A,Lda,S,Scond,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) :: Scond
      REAL , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPOEQU
   END INTERFACE
END MODULE S_SPOEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPORFS
   INTERFACE
      SUBROUTINE SPORFS(Uplo,N,Nrhs,A,Lda,Af,Ldaf,B,Ldb,X,Ldx,Ferr,Berr,&
     &                  Work,Iwork,Info)
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPORFS
   END INTERFACE
END MODULE S_SPORFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPORFSX
   INTERFACE
      SUBROUTINE SPORFSX(Uplo,Equed,N,Nrhs,A,Lda,Af,Ldaf,S,B,Ldb,X,Ldx, &
     &                   Rcond,Berr,N_err_bnds,Err_bnds_norm,           &
     &                   Err_bnds_comp,Nparams,Params,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ITREF_DEFAULT = 1.0 ,                       &
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      REAL , DIMENSION(*) :: S
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL , INTENT(INOUT) , DIMENSION(*) :: Params
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPORFSX
   END INTERFACE
END MODULE S_SPORFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPOSV
   INTERFACE
      SUBROUTINE SPOSV(Uplo,N,Nrhs,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPOSV
   END INTERFACE
END MODULE S_SPOSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPOSVX
   INTERFACE
      SUBROUTINE SPOSVX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Equed,S,B,Ldb,X, &
     &                  Ldx,Rcond,Ferr,Berr,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: S
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPOSVX
   END INTERFACE
END MODULE S_SPOSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPOSVXX
   INTERFACE
      SUBROUTINE SPOSVXX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Equed,S,B,Ldb,X,&
     &                   Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,Err_bnds_norm,&
     &                   Err_bnds_comp,Nparams,Params,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: S
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , INTENT(OUT) :: Rpvgrw
      REAL , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL , DIMENSION(*) :: Params
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPOSVXX
   END INTERFACE
END MODULE S_SPOSVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPOTF2
   INTERFACE
      SUBROUTINE SPOTF2(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPOTF2
   END INTERFACE
END MODULE S_SPOTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPOTRF2
   INTERFACE
      RECURSIVE SUBROUTINE SPOTRF2(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPOTRF2
   END INTERFACE
END MODULE S_SPOTRF2
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
MODULE S_SPOTRI
   INTERFACE
      SUBROUTINE SPOTRI(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPOTRI
   END INTERFACE
END MODULE S_SPOTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPOTRS
   INTERFACE
      SUBROUTINE SPOTRS(Uplo,N,Nrhs,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPOTRS
   END INTERFACE
END MODULE S_SPOTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPPCON
   INTERFACE
      SUBROUTINE SPPCON(Uplo,N,Ap,Anorm,Rcond,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(*) :: Ap
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPPCON
   END INTERFACE
END MODULE S_SPPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPPEQU
   INTERFACE
      SUBROUTINE SPPEQU(Uplo,N,Ap,S,Scond,Amax,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: Ap
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) :: Scond
      REAL , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPPEQU
   END INTERFACE
END MODULE S_SPPEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPPRFS
   INTERFACE
      SUBROUTINE SPPRFS(Uplo,N,Nrhs,Ap,Afp,B,Ldb,X,Ldx,Ferr,Berr,Work,  &
     &                  Iwork,Info)
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
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: Afp
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPPRFS
   END INTERFACE
END MODULE S_SPPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPPSV
   INTERFACE
      SUBROUTINE SPPSV(Uplo,N,Nrhs,Ap,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPPSV
   END INTERFACE
END MODULE S_SPPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPPSVX
   INTERFACE
      SUBROUTINE SPPSVX(Fact,Uplo,N,Nrhs,Ap,Afp,Equed,S,B,Ldb,X,Ldx,    &
     &                  Rcond,Ferr,Berr,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: Afp
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: S
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPPSVX
   END INTERFACE
END MODULE S_SPPSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPPTRF
   INTERFACE
      SUBROUTINE SPPTRF(Uplo,N,Ap,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPPTRF
   END INTERFACE
END MODULE S_SPPTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPPTRI
   INTERFACE
      SUBROUTINE SPPTRI(Uplo,N,Ap,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPPTRI
   END INTERFACE
END MODULE S_SPPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPPTRS
   INTERFACE
      SUBROUTINE SPPTRS(Uplo,N,Nrhs,Ap,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPPTRS
   END INTERFACE
END MODULE S_SPPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPSTF2
   INTERFACE
      SUBROUTINE SPSTF2(Uplo,N,A,Lda,Piv,Rank,Tol,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(N) :: Piv
      INTEGER , INTENT(OUT) :: Rank
      REAL , INTENT(IN) :: Tol
      REAL , INTENT(INOUT) , DIMENSION(2*N) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPSTF2
   END INTERFACE
END MODULE S_SPSTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPSTRF
   INTERFACE
      SUBROUTINE SPSTRF(Uplo,N,A,Lda,Piv,Rank,Tol,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(N) :: Piv
      INTEGER :: Rank
      REAL :: Tol
      REAL , INTENT(INOUT) , DIMENSION(2*N) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPSTRF
   END INTERFACE
END MODULE S_SPSTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPTCON
   INTERFACE
      SUBROUTINE SPTCON(N,D,E,Anorm,Rcond,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER :: N
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: E
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPTCON
   END INTERFACE
END MODULE S_SPTCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPTEQR
   INTERFACE
      SUBROUTINE SPTEQR(Compz,N,D,E,Z,Ldz,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Compz
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPTEQR
   END INTERFACE
END MODULE S_SPTEQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPTRFS
   INTERFACE
      SUBROUTINE SPTRFS(N,Nrhs,D,E,Df,Ef,B,Ldb,X,Ldx,Ferr,Berr,Work,    &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , THREE = 3.0E+0
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL , INTENT(IN) , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Df
      REAL , DIMENSION(*) :: Ef
      REAL , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPTRFS
   END INTERFACE
END MODULE S_SPTRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPTSV
   INTERFACE
      SUBROUTINE SPTSV(N,Nrhs,D,E,B,Ldb,Info)
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPTSV
   END INTERFACE
END MODULE S_SPTSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPTSVX
   INTERFACE
      SUBROUTINE SPTSVX(Fact,N,Nrhs,D,E,Df,Ef,B,Ldb,X,Ldx,Rcond,Ferr,   &
     &                  Berr,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Fact
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Df
      REAL , DIMENSION(*) :: Ef
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPTSVX
   END INTERFACE
END MODULE S_SPTSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPTTRF
   INTERFACE
      SUBROUTINE SPTTRF(N,D,E,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER :: Info
      END SUBROUTINE SPTTRF
   END INTERFACE
END MODULE S_SPTTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPTTRS
   INTERFACE
      SUBROUTINE SPTTRS(N,Nrhs,D,E,B,Ldb,Info)
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SPTTRS
   END INTERFACE
END MODULE S_SPTTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SPTTS2
   INTERFACE
      SUBROUTINE SPTTS2(N,Nrhs,D,E,B,Ldb)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: D
      REAL , INTENT(IN) , DIMENSION(*) :: E
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      END SUBROUTINE SPTTS2
   END INTERFACE
END MODULE S_SPTTS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SRSCL
   INTERFACE
      SUBROUTINE SRSCL(N,Sa,Sx,Incx)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      INTEGER :: N
      REAL , INTENT(IN) :: Sa
      REAL , DIMENSION(*) :: Sx
      INTEGER :: Incx
      END SUBROUTINE SRSCL
   END INTERFACE
END MODULE S_SRSCL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSB2ST_KERNELS
   INTERFACE
      SUBROUTINE SSB2ST_KERNELS(Uplo,Wantz,Ttype,St,Ed,Sweep,N,Nb,Ib,A, &
     &                          Lda,V,Tau,Ldvt,Work)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: Ttype
      INTEGER , INTENT(IN) :: St
      INTEGER , INTENT(IN) :: Ed
      INTEGER , INTENT(IN) :: Sweep
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(IN) :: Ib
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , DIMENSION(*) :: V
      REAL , DIMENSION(*) :: Tau
      INTEGER , INTENT(IN) :: Ldvt
      REAL , DIMENSION(*) :: Work
      END SUBROUTINE SSB2ST_KERNELS
   END INTERFACE
END MODULE S_SSB2ST_KERNELS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSBEV_2STAGE
   INTERFACE
      SUBROUTINE SSBEV_2STAGE(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Lwork,&
     &                        Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSBEV_2STAGE
   END INTERFACE
END MODULE S_SSBEV_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSBEVD_2STAGE
   INTERFACE
      SUBROUTINE SSBEVD_2STAGE(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,     &
     &                         Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSBEVD_2STAGE
   END INTERFACE
END MODULE S_SSBEVD_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSBEVD
   INTERFACE
      SUBROUTINE SSBEVD(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Lwork,Iwork,&
     &                  Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSBEVD
   END INTERFACE
END MODULE S_SSBEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSBEV
   INTERFACE
      SUBROUTINE SSBEV(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSBEV
   END INTERFACE
END MODULE S_SSBEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSBEVX_2STAGE
   INTERFACE
      SUBROUTINE SSBEVX_2STAGE(Jobz,Range,Uplo,N,Kd,Ab,Ldab,Q,Ldq,Vl,Vu,&
     &                         Il,Iu,Abstol,M,W,Z,Ldz,Work,Lwork,Iwork, &
     &                         Ifail,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSBEVX_2STAGE
   END INTERFACE
END MODULE S_SSBEVX_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSBEVX
   INTERFACE
      SUBROUTINE SSBEVX(Jobz,Range,Uplo,N,Kd,Ab,Ldab,Q,Ldq,Vl,Vu,Il,Iu, &
     &                  Abstol,M,W,Z,Ldz,Work,Iwork,Ifail,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSBEVX
   END INTERFACE
END MODULE S_SSBEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSBGST
   INTERFACE
      SUBROUTINE SSBGST(Vect,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,X,Ldx,Work,   &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ka
      INTEGER , INTENT(IN) :: Kb
      REAL , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , DIMENSION(Ldbb,*) :: Bb
      INTEGER , INTENT(IN) :: Ldbb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSBGST
   END INTERFACE
END MODULE S_SSBGST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSBGVD
   INTERFACE
      SUBROUTINE SSBGVD(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,W,Z,Ldz,Work, &
     &                  Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ka
      INTEGER :: Kb
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldbb,*) :: Bb
      INTEGER :: Ldbb
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSBGVD
   END INTERFACE
END MODULE S_SSBGVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSBGV
   INTERFACE
      SUBROUTINE SSBGV(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,W,Z,Ldz,Work,  &
     &                 Info)
      IMPLICIT NONE
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ka
      INTEGER :: Kb
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldbb,*) :: Bb
      INTEGER :: Ldbb
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSBGV
   END INTERFACE
END MODULE S_SSBGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSBGVX
   INTERFACE
      SUBROUTINE SSBGVX(Jobz,Range,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,Q,Ldq,  &
     &                  Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,Work,Iwork,Ifail,  &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ka
      INTEGER :: Kb
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldbb,*) :: Bb
      INTEGER :: Ldbb
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSBGVX
   END INTERFACE
END MODULE S_SSBGVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSBTRD
   INTERFACE
      SUBROUTINE SSBTRD(Vect,Uplo,N,Kd,Ab,Ldab,D,E,Q,Ldq,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kd
      REAL , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSBTRD
   END INTERFACE
END MODULE S_SSBTRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSFRK
   INTERFACE
      SUBROUTINE SSFRK(Transr,Uplo,Trans,N,K,Alpha,A,Lda,Beta,C)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Transr
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: K
      REAL :: Alpha
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL :: Beta
      REAL , DIMENSION(*) :: C
      END SUBROUTINE SSFRK
   END INTERFACE
END MODULE S_SSFRK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPCON
   INTERFACE
      SUBROUTINE SSPCON(Uplo,N,Ap,Ipiv,Anorm,Rcond,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(*) :: Ap
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPCON
   END INTERFACE
END MODULE S_SSPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPEVD
   INTERFACE
      SUBROUTINE SSPEVD(Jobz,Uplo,N,Ap,W,Z,Ldz,Work,Lwork,Iwork,Liwork, &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPEVD
   END INTERFACE
END MODULE S_SSPEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPEV
   INTERFACE
      SUBROUTINE SSPEV(Jobz,Uplo,N,Ap,W,Z,Ldz,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPEV
   END INTERFACE
END MODULE S_SSPEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPEVX
   INTERFACE
      SUBROUTINE SSPEVX(Jobz,Range,Uplo,N,Ap,Vl,Vu,Il,Iu,Abstol,M,W,Z,  &
     &                  Ldz,Work,Iwork,Ifail,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(*) :: Ap
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPEVX
   END INTERFACE
END MODULE S_SSPEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPGST
   INTERFACE
      SUBROUTINE SSPGST(Itype,Uplo,N,Ap,Bp,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0 , HALF = 0.5
      INTEGER , INTENT(IN) :: Itype
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: Bp
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPGST
   END INTERFACE
END MODULE S_SSPGST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPGVD
   INTERFACE
      SUBROUTINE SSPGVD(Itype,Jobz,Uplo,N,Ap,Bp,W,Z,Ldz,Work,Lwork,     &
     &                  Iwork,Liwork,Info)
      IMPLICIT NONE
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: Bp
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPGVD
   END INTERFACE
END MODULE S_SSPGVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPGV
   INTERFACE
      SUBROUTINE SSPGV(Itype,Jobz,Uplo,N,Ap,Bp,W,Z,Ldz,Work,Info)
      IMPLICIT NONE
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: Bp
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPGV
   END INTERFACE
END MODULE S_SSPGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPGVX
   INTERFACE
      SUBROUTINE SSPGVX(Itype,Jobz,Range,Uplo,N,Ap,Bp,Vl,Vu,Il,Iu,      &
     &                  Abstol,M,W,Z,Ldz,Work,Iwork,Ifail,Info)
      IMPLICIT NONE
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: Bp
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPGVX
   END INTERFACE
END MODULE S_SSPGVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPRFS
   INTERFACE
      SUBROUTINE SSPRFS(Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Ferr,Berr,  &
     &                  Work,Iwork,Info)
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
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPRFS
   END INTERFACE
END MODULE S_SSPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPSV
   INTERFACE
      SUBROUTINE SSPSV(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: Ap
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPSV
   END INTERFACE
END MODULE S_SSPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPSVX
   INTERFACE
      SUBROUTINE SSPSVX(Fact,Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Rcond, &
     &                  Ferr,Berr,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPSVX
   END INTERFACE
END MODULE S_SSPSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPTRD
   INTERFACE
      SUBROUTINE SSPTRD(Uplo,N,Ap,D,E,Tau,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0 , ZERO = 0.0 , HALF = 1.0/2.0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL , INTENT(OUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPTRD
   END INTERFACE
END MODULE S_SSPTRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPTRF
   INTERFACE
      SUBROUTINE SSPTRF(Uplo,N,Ap,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPTRF
   END INTERFACE
END MODULE S_SSPTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPTRI
   INTERFACE
      SUBROUTINE SSPTRI(Uplo,N,Ap,Ipiv,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPTRI
   END INTERFACE
END MODULE S_SSPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSPTRS
   INTERFACE
      SUBROUTINE SSPTRS(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSPTRS
   END INTERFACE
END MODULE S_SSPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSTEBZ
   INTERFACE
      SUBROUTINE SSTEBZ(Range,Order,N,Vl,Vu,Il,Iu,Abstol,D,E,M,Nsplit,W,&
     &                  Iblock,Isplit,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      HALF = 1.0E0/TWO , FUDGE = 2.1E0 ,          &
     &                      RELFAC = 2.0E0
      CHARACTER :: Range
      CHARACTER :: Order
      INTEGER :: N
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      REAL , INTENT(IN) :: Abstol
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) :: M
      INTEGER , INTENT(INOUT) :: Nsplit
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isplit
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSTEBZ
   END INTERFACE
END MODULE S_SSTEBZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSTEDC
   INTERFACE
      SUBROUTINE SSTEDC(Compz,N,D,E,Z,Ldz,Work,Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
      CHARACTER :: Compz
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSTEDC
   END INTERFACE
END MODULE S_SSTEDC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSTEGR
   INTERFACE
      SUBROUTINE SSTEGR(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,  &
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
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER :: Info
      END SUBROUTINE SSTEGR
   END INTERFACE
END MODULE S_SSTEGR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSTEIN
   INTERFACE
      SUBROUTINE SSTEIN(N,D,E,M,W,Iblock,Isplit,Z,Ldz,Work,Iwork,Ifail, &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
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
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSTEIN
   END INTERFACE
END MODULE S_SSTEIN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSTEMR
   INTERFACE
      SUBROUTINE SSTEMR(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,M,W,Z,Ldz,Nzc,     &
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
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(IN) :: Nzc
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isuppz
      LOGICAL , INTENT(INOUT) :: Tryrac
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSTEMR
   END INTERFACE
END MODULE S_SSTEMR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSTEQR
   INTERFACE
      SUBROUTINE SSTEQR(Compz,N,D,E,Z,Ldz,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      THREE = 3.0E0
      INTEGER , PARAMETER  ::  MAXIT = 30
      CHARACTER :: Compz
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSTEQR
   END INTERFACE
END MODULE S_SSTEQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSTERF
   INTERFACE
      SUBROUTINE SSTERF(N,D,E,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      THREE = 3.0E0
      INTEGER , PARAMETER  ::  MAXIT = 30
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSTERF
   END INTERFACE
END MODULE S_SSTERF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSTEVD
   INTERFACE
      SUBROUTINE SSTEVD(Jobz,N,D,E,Z,Ldz,Work,Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      INTEGER :: N
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSTEVD
   END INTERFACE
END MODULE S_SSTEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSTEV
   INTERFACE
      SUBROUTINE SSTEV(Jobz,N,D,E,Z,Ldz,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      INTEGER :: N
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSTEV
   END INTERFACE
END MODULE S_SSTEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSTEVR
   INTERFACE
      SUBROUTINE SSTEVR(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,  &
     &                  Isuppz,Work,Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
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
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSTEVR
   END INTERFACE
END MODULE S_SSTEVR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSTEVX
   INTERFACE
      SUBROUTINE SSTEVX(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,  &
     &                  Work,Iwork,Ifail,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      CHARACTER :: Range
      INTEGER :: N
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSTEVX
   END INTERFACE
END MODULE S_SSTEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYCON_3
   INTERFACE
      SUBROUTINE SSYCON_3(Uplo,N,A,Lda,E,Ipiv,Anorm,Rcond,Work,Iwork,   &
     &                    Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYCON_3
   END INTERFACE
END MODULE S_SSYCON_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYCON
   INTERFACE
      SUBROUTINE SSYCON(Uplo,N,A,Lda,Ipiv,Anorm,Rcond,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYCON
   END INTERFACE
END MODULE S_SSYCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYCON_ROOK
   INTERFACE
      SUBROUTINE SSYCON_ROOK(Uplo,N,A,Lda,Ipiv,Anorm,Rcond,Work,Iwork,  &
     &                       Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYCON_ROOK
   END INTERFACE
END MODULE S_SSYCON_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYCONV
   INTERFACE
      SUBROUTINE SSYCONV(Uplo,Way,N,A,Lda,Ipiv,E,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYCONV
   END INTERFACE
END MODULE S_SSYCONV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYCONVF
   INTERFACE
      SUBROUTINE SSYCONVF(Uplo,Way,N,A,Lda,E,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYCONVF
   END INTERFACE
END MODULE S_SSYCONVF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYCONVF_ROOK
   INTERFACE
      SUBROUTINE SSYCONVF_ROOK(Uplo,Way,N,A,Lda,E,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYCONVF_ROOK
   END INTERFACE
END MODULE S_SSYCONVF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYEQUB
   INTERFACE
      SUBROUTINE SSYEQUB(Uplo,N,A,Lda,S,Scond,Amax,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , ZERO = 0.0E0
      INTEGER , PARAMETER  ::  MAX_ITER = 100
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(OUT) :: Scond
      REAL , INTENT(INOUT) :: Amax
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYEQUB
   END INTERFACE
END MODULE S_SSYEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYEV_2STAGE
   INTERFACE
      SUBROUTINE SSYEV_2STAGE(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYEV_2STAGE
   END INTERFACE
END MODULE S_SSYEV_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYEVD_2STAGE
   INTERFACE
      SUBROUTINE SSYEVD_2STAGE(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Iwork,    &
     &                         Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYEVD_2STAGE
   END INTERFACE
END MODULE S_SSYEVD_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYEVD
   INTERFACE
      SUBROUTINE SSYEVD(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Iwork,Liwork,    &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYEVD
   END INTERFACE
END MODULE S_SSYEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYEV
   INTERFACE
      SUBROUTINE SSYEV(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYEV
   END INTERFACE
END MODULE S_SSYEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYEVR_2STAGE
   INTERFACE
      SUBROUTINE SSYEVR_2STAGE(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,     &
     &                         Abstol,M,W,Z,Ldz,Isuppz,Work,Lwork,Iwork,&
     &                         Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYEVR_2STAGE
   END INTERFACE
END MODULE S_SSYEVR_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYEVR
   INTERFACE
      SUBROUTINE SSYEVR(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,Abstol,M,W, &
     &                  Z,Ldz,Isuppz,Work,Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYEVR
   END INTERFACE
END MODULE S_SSYEVR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYEVX_2STAGE
   INTERFACE
      SUBROUTINE SSYEVX_2STAGE(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,     &
     &                         Abstol,M,W,Z,Ldz,Work,Lwork,Iwork,Ifail, &
     &                         Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYEVX_2STAGE
   END INTERFACE
END MODULE S_SSYEVX_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYEVX
   INTERFACE
      SUBROUTINE SSYEVX(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,Abstol,M,W, &
     &                  Z,Ldz,Work,Lwork,Iwork,Ifail,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(IN) :: Vl
      REAL , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYEVX
   END INTERFACE
END MODULE S_SSYEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYGS2
   INTERFACE
      SUBROUTINE SSYGS2(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0 , HALF = 0.5
      INTEGER , INTENT(IN) :: Itype
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYGS2
   END INTERFACE
END MODULE S_SSYGS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYGST
   INTERFACE
      SUBROUTINE SSYGST(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0 , HALF = 0.5
      INTEGER :: Itype
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYGST
   END INTERFACE
END MODULE S_SSYGST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYGV_2STAGE
   INTERFACE
      SUBROUTINE SSYGV_2STAGE(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,     &
     &                        Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYGV_2STAGE
   END INTERFACE
END MODULE S_SSYGV_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYGVD
   INTERFACE
      SUBROUTINE SSYGVD(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,Lwork,     &
     &                  Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: W
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYGVD
   END INTERFACE
END MODULE S_SSYGVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYGV
   INTERFACE
      SUBROUTINE SSYGV(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYGV
   END INTERFACE
END MODULE S_SSYGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYGVX
   INTERFACE
      SUBROUTINE SSYGVX(Itype,Jobz,Range,Uplo,N,A,Lda,B,Ldb,Vl,Vu,Il,Iu,&
     &                  Abstol,M,W,Z,Ldz,Work,Lwork,Iwork,Ifail,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL :: Abstol
      INTEGER :: M
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYGVX
   END INTERFACE
END MODULE S_SSYGVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYRFS
   INTERFACE
      SUBROUTINE SSYRFS(Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,Ferr,&
     &                  Berr,Work,Iwork,Info)
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYRFS
   END INTERFACE
END MODULE S_SSYRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYRFSX
   INTERFACE
      SUBROUTINE SSYRFSX(Uplo,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,S,B,Ldb,X,&
     &                   Ldx,Rcond,Berr,N_err_bnds,Err_bnds_norm,       &
     &                   Err_bnds_comp,Nparams,Params,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ITREF_DEFAULT = 1.0 ,                       &
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
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: S
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL , INTENT(INOUT) , DIMENSION(*) :: Params
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYRFSX
   END INTERFACE
END MODULE S_SSYRFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYSV_AA_2STAGE
   INTERFACE
      SUBROUTINE SSYSV_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B, &
     &                           Ldb,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tb
      INTEGER :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYSV_AA_2STAGE
   END INTERFACE
END MODULE S_SSYSV_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYSV_AA
   INTERFACE
      SUBROUTINE SSYSV_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYSV_AA
   END INTERFACE
END MODULE S_SSYSV_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYSV
   INTERFACE
      SUBROUTINE SSYSV(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYSV
   END INTERFACE
END MODULE S_SSYSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYSV_RK
   INTERFACE
      SUBROUTINE SSYSV_RK(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Work,Lwork,    &
     &                    Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYSV_RK
   END INTERFACE
END MODULE S_SSYSV_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYSV_ROOK
   INTERFACE
      SUBROUTINE SSYSV_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,    &
     &                      Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYSV_ROOK
   END INTERFACE
END MODULE S_SSYSV_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYSVX
   INTERFACE
      SUBROUTINE SSYSVX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,&
     &                  Rcond,Ferr,Berr,Work,Lwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYSVX
   END INTERFACE
END MODULE S_SSYSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYSVXX
   INTERFACE
      SUBROUTINE SSYSVXX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,S,B, &
     &                   Ldb,X,Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,        &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL , DIMENSION(*) :: S
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , INTENT(OUT) :: Rpvgrw
      REAL , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL , DIMENSION(*) :: Params
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYSVXX
   END INTERFACE
END MODULE S_SSYSVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYSWAPR
   INTERFACE
      SUBROUTINE SSYSWAPR(Uplo,N,A,Lda,I1,I2)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,N) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) :: I1
      INTEGER , INTENT(IN) :: I2
      END SUBROUTINE SSYSWAPR
   END INTERFACE
END MODULE S_SSYSWAPR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTD2
   INTERFACE
      SUBROUTINE SSYTD2(Uplo,N,A,Lda,D,E,Tau,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0 , ZERO = 0.0 , HALF = 1.0/2.0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(OUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTD2
   END INTERFACE
END MODULE S_SSYTD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTF2
   INTERFACE
      SUBROUTINE SSYTF2(Uplo,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTF2
   END INTERFACE
END MODULE S_SSYTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTF2_RK
   INTERFACE
      SUBROUTINE SSYTF2_RK(Uplo,N,A,Lda,E,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTF2_RK
   END INTERFACE
END MODULE S_SSYTF2_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTF2_ROOK
   INTERFACE
      SUBROUTINE SSYTF2_ROOK(Uplo,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      EIGHT = 8.0E+0 , SEVTEN = 17.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTF2_ROOK
   END INTERFACE
END MODULE S_SSYTF2_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRD_2STAGE
   INTERFACE
      SUBROUTINE SSYTRD_2STAGE(Vect,Uplo,N,A,Lda,D,E,Tau,Hous2,Lhous2,  &
     &                         Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Hous2
      INTEGER :: Lhous2
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRD_2STAGE
   END INTERFACE
END MODULE S_SSYTRD_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRD
   INTERFACE
      SUBROUTINE SSYTRD(Uplo,N,A,Lda,D,E,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRD
   END INTERFACE
END MODULE S_SSYTRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRD_SB2ST
   INTERFACE
      SUBROUTINE SSYTRD_SB2ST(Stage1,Vect,Uplo,N,Kd,Ab,Ldab,D,E,Hous,   &
     &                        Lhous,Work,Lwork,Info)
      USE OMP_LIB
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  RZERO = 0.0E+0 , ZERO = 0.0E+0 ,            &
     &                      ONE = 1.0E+0
      CHARACTER :: Stage1
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , INTENT(OUT) , DIMENSION(*) :: D
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Hous
      INTEGER , INTENT(IN) :: Lhous
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRD_SB2ST
   END INTERFACE
END MODULE S_SSYTRD_SB2ST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRD_SY2SB
   INTERFACE
      SUBROUTINE SSYTRD_SY2SB(Uplo,N,Kd,A,Lda,Ab,Ldab,Tau,Work,Lwork,   &
     &                        Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  RONE = 1.0E+0 , ZERO = 0.0E+0 ,             &
     &                      ONE = 1.0E+0 , HALF = 0.5E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRD_SY2SB
   END INTERFACE
END MODULE S_SSYTRD_SY2SB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRF_AA_2STAGE
   INTERFACE
      SUBROUTINE SSYTRF_AA_2STAGE(Uplo,N,A,Lda,Tb,Ltb,Ipiv,Ipiv2,Work,  &
     &                            Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRF_AA_2STAGE
   END INTERFACE
END MODULE S_SSYTRF_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRF_AA
   INTERFACE
      SUBROUTINE SSYTRF_AA(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRF_AA
   END INTERFACE
END MODULE S_SSYTRF_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRF
   INTERFACE
      SUBROUTINE SSYTRF(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRF
   END INTERFACE
END MODULE S_SSYTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRF_RK
   INTERFACE
      SUBROUTINE SSYTRF_RK(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRF_RK
   END INTERFACE
END MODULE S_SSYTRF_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRF_ROOK
   INTERFACE
      SUBROUTINE SSYTRF_ROOK(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRF_ROOK
   END INTERFACE
END MODULE S_SSYTRF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRI2
   INTERFACE
      SUBROUTINE SSYTRI2(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRI2
   END INTERFACE
END MODULE S_SSYTRI2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRI2X
   INTERFACE
      SUBROUTINE SSYTRI2X(Uplo,N,A,Lda,Ipiv,Work,Nb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRI2X
   END INTERFACE
END MODULE S_SSYTRI2X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRI_3
   INTERFACE
      SUBROUTINE SSYTRI_3(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRI_3
   END INTERFACE
END MODULE S_SSYTRI_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRI_3X
   INTERFACE
      SUBROUTINE SSYTRI_3X(Uplo,N,A,Lda,E,Ipiv,Work,Nb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRI_3X
   END INTERFACE
END MODULE S_SSYTRI_3X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRI
   INTERFACE
      SUBROUTINE SSYTRI(Uplo,N,A,Lda,Ipiv,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRI
   END INTERFACE
END MODULE S_SSYTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRI_ROOK
   INTERFACE
      SUBROUTINE SSYTRI_ROOK(Uplo,N,A,Lda,Ipiv,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRI_ROOK
   END INTERFACE
END MODULE S_SSYTRI_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRS2
   INTERFACE
      SUBROUTINE SSYTRS2(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRS2
   END INTERFACE
END MODULE S_SSYTRS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRS_3
   INTERFACE
      SUBROUTINE SSYTRS_3(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRS_3
   END INTERFACE
END MODULE S_SSYTRS_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRS_AA_2STAGE
   INTERFACE
      SUBROUTINE SSYTRS_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B,&
     &                            Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRS_AA_2STAGE
   END INTERFACE
END MODULE S_SSYTRS_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRS_AA
   INTERFACE
      SUBROUTINE SSYTRS_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRS_AA
   END INTERFACE
END MODULE S_SSYTRS_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRS
   INTERFACE
      SUBROUTINE SSYTRS(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRS
   END INTERFACE
END MODULE S_SSYTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_SSYTRS_ROOK
   INTERFACE
      SUBROUTINE SSYTRS_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE SSYTRS_ROOK
   END INTERFACE
END MODULE S_SSYTRS_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STBCON
   INTERFACE
      SUBROUTINE STBCON(Norm,Uplo,Diag,N,Kd,Ab,Ldab,Rcond,Work,Iwork,   &
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
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STBCON
   END INTERFACE
END MODULE S_STBCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STBRFS
   INTERFACE
      SUBROUTINE STBRFS(Uplo,Trans,Diag,N,Kd,Nrhs,Ab,Ldab,B,Ldb,X,Ldx,  &
     &                  Ferr,Berr,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(OUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STBRFS
   END INTERFACE
END MODULE S_STBRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STBTRS
   INTERFACE
      SUBROUTINE STBTRS(Uplo,Trans,Diag,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      REAL , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STBTRS
   END INTERFACE
END MODULE S_STBTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STFSM
   INTERFACE
      SUBROUTINE STFSM(Transr,Side,Uplo,Trans,Diag,M,N,Alpha,A,B,Ldb)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Transr
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: M
      INTEGER :: N
      REAL :: Alpha
      REAL , DIMENSION(0:*) :: A
      REAL , DIMENSION(0:Ldb-1,0:*) :: B
      INTEGER :: Ldb
      END SUBROUTINE STFSM
   END INTERFACE
END MODULE S_STFSM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STFTRI
   INTERFACE
      SUBROUTINE STFTRI(Transr,Uplo,Diag,N,A,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Transr
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      REAL , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STFTRI
   END INTERFACE
END MODULE S_STFTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STFTTP
   INTERFACE
      SUBROUTINE STFTTP(Transr,Uplo,N,Arf,Ap,Info)
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(0:*) :: Arf
      REAL , INTENT(OUT) , DIMENSION(0:*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STFTTP
   END INTERFACE
END MODULE S_STFTTP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STFTTR
   INTERFACE
      SUBROUTINE STFTTR(Transr,Uplo,N,Arf,A,Lda,Info)
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(0:*) :: Arf
      REAL , INTENT(OUT) , DIMENSION(0:Lda-1,0:*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STFTTR
   END INTERFACE
END MODULE S_STFTTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STGEVC
   INTERFACE
      SUBROUTINE STGEVC(Side,Howmny,Select,N,S,Lds,P,Ldp,Vl,Ldvl,Vr,    &
     &                  Ldvr,Mm,M,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      SAFETY = 1.0E+2
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL , DIMENSION(Lds,*) :: S
      INTEGER :: Lds
      REAL , DIMENSION(Ldp,*) :: P
      INTEGER :: Ldp
      REAL , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(OUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STGEVC
   END INTERFACE
END MODULE S_STGEVC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STGEX2
   INTERFACE
      SUBROUTINE STGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,J1,N1,N2, &
     &                  Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWENTY = 2.0E+01
      INTEGER , PARAMETER  ::  LDST = 4
      LOGICAL , PARAMETER  ::  WANDS = .TRUE.
      LOGICAL , INTENT(IN) :: Wantq
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(IN) :: J1
      INTEGER :: N1
      INTEGER :: N2
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER :: Info
      END SUBROUTINE STGEX2
   END INTERFACE
END MODULE S_STGEX2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STGEXC
   INTERFACE
      SUBROUTINE STGEXC(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,Ifst,Ilst,&
     &                  Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      LOGICAL :: Wantq
      LOGICAL :: Wantz
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ifst
      INTEGER , INTENT(INOUT) :: Ilst
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STGEXC
   END INTERFACE
END MODULE S_STGEXC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STGSEN
   INTERFACE
      SUBROUTINE STGSEN(Ijob,Wantq,Wantz,Select,N,A,Lda,B,Ldb,Alphar,   &
     &                  Alphai,Beta,Q,Ldq,Z,Ldz,M,Pl,Pr,Dif,Work,Lwork, &
     &                  Iwork,Liwork,Info)
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
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Alphar
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL , DIMENSION(*) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) :: Pl
      REAL , INTENT(INOUT) :: Pr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dif
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STGSEN
   END INTERFACE
END MODULE S_STGSEN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STGSJA
   INTERFACE
      SUBROUTINE STGSJA(Jobu,Jobv,Jobq,M,P,N,K,L,A,Lda,B,Ldb,Tola,Tolb, &
     &                  Alpha,Beta,U,Ldu,V,Ldv,Q,Ldq,Work,Ncycle,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXIT = 40
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      CHARACTER :: Jobq
      INTEGER :: M
      INTEGER :: P
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER :: L
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(IN) :: Tola
      REAL , INTENT(IN) :: Tolb
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alpha
      REAL , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) :: Ncycle
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STGSJA
   END INTERFACE
END MODULE S_STGSJA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STGSNA
   INTERFACE
      SUBROUTINE STGSNA(Job,Howmny,Select,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,&
     &                  S,Dif,Mm,M,Work,Lwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  DIFDRI = 3
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , FOUR = 4.0E+0
      CHARACTER :: Job
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      REAL , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dif
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STGSNA
   END INTERFACE
END MODULE S_STGSNA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STGSY2
   INTERFACE
      SUBROUTINE STGSY2(Trans,Ijob,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F, &
     &                  Ldf,Scale,Rdsum,Rdscal,Iwork,Pq,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  LDZ = 8
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Trans
      INTEGER :: Ijob
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(Ldd,*) :: D
      INTEGER :: Ldd
      REAL , DIMENSION(Lde,*) :: E
      INTEGER :: Lde
      REAL , INTENT(INOUT) , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      REAL , INTENT(INOUT) :: Scale
      REAL :: Rdsum
      REAL :: Rdscal
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) :: Pq
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STGSY2
   END INTERFACE
END MODULE S_STGSY2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STGSYL
   INTERFACE
      SUBROUTINE STGSYL(Trans,Ijob,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F, &
     &                  Ldf,Scale,Dif,Work,Lwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: Ijob
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , DIMENSION(Ldd,*) :: D
      INTEGER :: Ldd
      REAL , DIMENSION(Lde,*) :: E
      INTEGER :: Lde
      REAL , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      REAL , INTENT(INOUT) :: Scale
      REAL , INTENT(OUT) :: Dif
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STGSYL
   END INTERFACE
END MODULE S_STGSYL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPCON
   INTERFACE
      SUBROUTINE STPCON(Norm,Uplo,Diag,N,Ap,Rcond,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      REAL , DIMENSION(*) :: Ap
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STPCON
   END INTERFACE
END MODULE S_STPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPLQT2
   INTERFACE
      SUBROUTINE STPLQT2(M,N,L,A,Lda,B,Ldb,T,Ldt,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0 , ZERO = 0.0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: L
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STPLQT2
   END INTERFACE
END MODULE S_STPLQT2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPLQT
   INTERFACE
      SUBROUTINE STPLQT(M,N,L,Mb,A,Lda,B,Ldb,T,Ldt,Work,Info)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Mb
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STPLQT
   END INTERFACE
END MODULE S_STPLQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPMLQT
   INTERFACE
      SUBROUTINE STPMLQT(Side,Trans,M,N,K,L,Mb,V,Ldv,T,Ldt,A,Lda,B,Ldb, &
     &                   Work,Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Mb
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STPMLQT
   END INTERFACE
END MODULE S_STPMLQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPMQRT
   INTERFACE
      SUBROUTINE STPMQRT(Side,Trans,M,N,K,L,Nb,V,Ldv,T,Ldt,A,Lda,B,Ldb, &
     &                   Work,Info)
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Nb
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STPMQRT
   END INTERFACE
END MODULE S_STPMQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPQRT2
   INTERFACE
      SUBROUTINE STPQRT2(M,N,L,A,Lda,B,Ldb,T,Ldt,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0 , ZERO = 0.0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: L
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STPQRT2
   END INTERFACE
END MODULE S_STPQRT2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPQRT
   INTERFACE
      SUBROUTINE STPQRT(M,N,L,Nb,A,Lda,B,Ldb,T,Ldt,Work,Info)
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Nb
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STPQRT
   END INTERFACE
END MODULE S_STPQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPRFB
   INTERFACE
      SUBROUTINE STPRFB(Side,Trans,Direct,Storev,M,N,K,L,V,Ldv,T,Ldt,A, &
     &                  Lda,B,Ldb,Work,Ldwork)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0 , ZERO = 0.0
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: L
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE STPRFB
   END INTERFACE
END MODULE S_STPRFB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPRFS
   INTERFACE
      SUBROUTINE STPRFS(Uplo,Trans,Diag,N,Nrhs,Ap,B,Ldb,X,Ldx,Ferr,Berr,&
     &                  Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(OUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STPRFS
   END INTERFACE
END MODULE S_STPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPTRI
   INTERFACE
      SUBROUTINE STPTRI(Uplo,Diag,N,Ap,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STPTRI
   END INTERFACE
END MODULE S_STPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPTRS
   INTERFACE
      SUBROUTINE STPTRS(Uplo,Trans,Diag,N,Nrhs,Ap,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL , DIMENSION(*) :: Ap
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STPTRS
   END INTERFACE
END MODULE S_STPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPTTF
   INTERFACE
      SUBROUTINE STPTTF(Transr,Uplo,N,Ap,Arf,Info)
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(0:*) :: Ap
      REAL , INTENT(OUT) , DIMENSION(0:*) :: Arf
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STPTTF
   END INTERFACE
END MODULE S_STPTTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STPTTR
   INTERFACE
      SUBROUTINE STPTTR(Uplo,N,Ap,A,Lda,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: Ap
      REAL , INTENT(OUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STPTTR
   END INTERFACE
END MODULE S_STPTTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRCON
   INTERFACE
      SUBROUTINE STRCON(Norm,Uplo,Diag,N,A,Lda,Rcond,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STRCON
   END INTERFACE
END MODULE S_STRCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STREVC3
   INTERFACE
      SUBROUTINE STREVC3(Side,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,Mm, &
     &                   M,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      INTEGER , PARAMETER  ::  NBMIN = 8 , NBMAX = 128
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(INOUT) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STREVC3
   END INTERFACE
END MODULE S_STREVC3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STREVC
   INTERFACE
      SUBROUTINE STREVC(Side,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,Mm,M,&
     &                  Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(INOUT) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STREVC
   END INTERFACE
END MODULE S_STREVC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STREXC
   INTERFACE
      SUBROUTINE STREXC(Compq,N,T,Ldt,Q,Ldq,Ifst,Ilst,Work,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      CHARACTER :: Compq
      INTEGER :: N
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , INTENT(INOUT) :: Ifst
      INTEGER , INTENT(INOUT) :: Ilst
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STREXC
   END INTERFACE
END MODULE S_STREXC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRRFS
   INTERFACE
      SUBROUTINE STRRFS(Uplo,Trans,Diag,N,Nrhs,A,Lda,B,Ldb,X,Ldx,Ferr,  &
     &                  Berr,Work,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL , INTENT(OUT) , DIMENSION(*) :: Berr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STRRFS
   END INTERFACE
END MODULE S_STRRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRSEN
   INTERFACE
      SUBROUTINE STRSEN(Job,Compq,Select,N,T,Ldt,Q,Ldq,Wr,Wi,M,S,Sep,   &
     &                  Work,Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Job
      CHARACTER :: Compq
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , INTENT(OUT) , DIMENSION(*) :: Wr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Wi
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(OUT) :: S
      REAL , INTENT(OUT) :: Sep
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STRSEN
   END INTERFACE
END MODULE S_STRSEN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRSNA
   INTERFACE
      SUBROUTINE STRSNA(Job,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,S,Sep,&
     &                  Mm,M,Work,Ldwork,Iwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      CHARACTER :: Job
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      REAL , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      REAL , INTENT(OUT) , DIMENSION(*) :: S
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sep
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STRSNA
   END INTERFACE
END MODULE S_STRSNA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRSYL
   INTERFACE
      SUBROUTINE STRSYL(Trana,Tranb,Isgn,M,N,A,Lda,B,Ldb,C,Ldc,Scale,   &
     &                  Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Trana
      CHARACTER :: Tranb
      INTEGER :: Isgn
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL , INTENT(INOUT) :: Scale
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STRSYL
   END INTERFACE
END MODULE S_STRSYL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRTI2
   INTERFACE
      SUBROUTINE STRTI2(Uplo,Diag,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STRTI2
   END INTERFACE
END MODULE S_STRTI2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRTRI
   INTERFACE
      SUBROUTINE STRTRI(Uplo,Diag,N,A,Lda,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STRTRI
   END INTERFACE
END MODULE S_STRTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRTRS
   INTERFACE
      SUBROUTINE STRTRS(Uplo,Trans,Diag,N,Nrhs,A,Lda,B,Ldb,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STRTRS
   END INTERFACE
END MODULE S_STRTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRTTF
   INTERFACE
      SUBROUTINE STRTTF(Transr,Uplo,N,A,Lda,Arf,Info)
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(0:Lda-1,0:*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(OUT) , DIMENSION(0:*) :: Arf
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STRTTF
   END INTERFACE
END MODULE S_STRTTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STRTTP
   INTERFACE
      SUBROUTINE STRTTP(Uplo,N,A,Lda,Ap,Info)
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(OUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STRTTP
   END INTERFACE
END MODULE S_STRTTP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_STZRZF
   INTERFACE
      SUBROUTINE STZRZF(M,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE STZRZF
   END INTERFACE
END MODULE S_STZRZF
