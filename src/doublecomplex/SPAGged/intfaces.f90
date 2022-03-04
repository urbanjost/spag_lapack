!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZBBCSD
   INTERFACE
      SUBROUTINE ZBBCSD(Jobu1,Jobu2,Jobv1t,Jobv2t,Trans,M,P,Q,Theta,Phi,&
     &                  U1,Ldu1,U2,Ldu2,V1t,Ldv1t,V2t,Ldv2t,B11d,B11e,  &
     &                  B12d,B12e,B21d,B21e,B22d,B22e,Rwork,Lrwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXITR = 6
      REAL(R8KIND) , PARAMETER  ::  HUNDRED = 100.0D0 ,                 &
     &                              MEIGHTH = -0.125D0 , ONE = 1.0D0 ,  &
     &                              TEN = 10.0D0 , ZERO = 0.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  NEGONECOMPLEX = (-1.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  PIOVER2 =                           &
     &                        1.57079632679489661923132169163975144210D0
      CHARACTER :: Jobu1
      CHARACTER :: Jobu2
      CHARACTER :: Jobv1t
      CHARACTER :: Jobv2t
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER :: P
      INTEGER :: Q
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Phi
      COMPLEX(CX16KIND) , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      COMPLEX(CX16KIND) , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      COMPLEX(CX16KIND) , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      COMPLEX(CX16KIND) , DIMENSION(Ldv2t,*) :: V2t
      INTEGER :: Ldv2t
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B11d
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B11e
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B12d
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B12e
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B21d
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B21e
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B22d
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B22e
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZBBCSD
   END INTERFACE
END MODULE S_ZBBCSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZBDSQR
   INTERFACE
      SUBROUTINE ZBDSQR(Uplo,N,Ncvt,Nru,Ncc,D,E,Vt,Ldvt,U,Ldu,C,Ldc,    &
     &                  Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              NEGONE = -1.0D0 , HNDRTH = 0.01D0 , &
     &                              TEN = 10.0D0 , HNDRD = 100.0D0 ,    &
     &                              MEIGTH = -0.125D0
      INTEGER , PARAMETER  ::  MAXITR = 6
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ncvt
      INTEGER :: Nru
      INTEGER :: Ncc
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      COMPLEX(CX16KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZBDSQR
   END INTERFACE
END MODULE S_ZBDSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZCGESV
   INTERFACE
      SUBROUTINE ZCGESV(N,Nrhs,A,Lda,Ipiv,B,Ldb,X,Ldx,Work,Swork,Rwork, &
     &                  Iter,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      LOGICAL , PARAMETER  ::  DOITREF = .TRUE.
      INTEGER , PARAMETER  ::  ITERMAX = 30
      REAL(R8KIND) , PARAMETER  ::  BWDMAX = 1.0E+00
      COMPLEX(CX16KIND) , PARAMETER  ::  NEGONE = (-1.0D+00,0.0D+00) ,  &
     &                 ONE = (1.0D+00,0.0D+00)
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      COMPLEX(CX16KIND) , DIMENSION(N,*) :: Work
      COMPLEX , DIMENSION(*) :: Swork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(OUT) :: Iter
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZCGESV
   END INTERFACE
END MODULE S_ZCGESV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZCPOSV
   INTERFACE
      SUBROUTINE ZCPOSV(Uplo,N,Nrhs,A,Lda,B,Ldb,X,Ldx,Work,Swork,Rwork, &
     &                  Iter,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      LOGICAL , PARAMETER  ::  DOITREF = .TRUE.
      INTEGER , PARAMETER  ::  ITERMAX = 30
      REAL(R8KIND) , PARAMETER  ::  BWDMAX = 1.0E+00
      COMPLEX(CX16KIND) , PARAMETER  ::  NEGONE = (-1.0D+00,0.0D+00) ,  &
     &                 ONE = (1.0D+00,0.0D+00)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      COMPLEX(CX16KIND) , DIMENSION(N,*) :: Work
      COMPLEX , DIMENSION(*) :: Swork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(OUT) :: Iter
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZCPOSV
   END INTERFACE
END MODULE S_ZCPOSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZDRSCL
   INTERFACE
      SUBROUTINE ZDRSCL(N,Sa,Sx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) :: Sa
      COMPLEX(CX16KIND) , DIMENSION(*) :: Sx
      INTEGER :: Incx
      END SUBROUTINE ZDRSCL
   END INTERFACE
END MODULE S_ZDRSCL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBBRD
   INTERFACE
      SUBROUTINE ZGBBRD(Vect,M,N,Ncc,Kl,Ku,Ab,Ldab,D,E,Q,Ldq,Pt,Ldpt,C, &
     &                  Ldc,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Vect
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Ncc
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(Ldpt,*) :: Pt
      INTEGER :: Ldpt
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGBBRD
   END INTERFACE
END MODULE S_ZGBBRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBCON
   INTERFACE
      SUBROUTINE ZGBCON(Norm,N,Kl,Ku,Ab,Ldab,Ipiv,Anorm,Rcond,Work,     &
     &                  Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Norm
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGBCON
   END INTERFACE
END MODULE S_ZGBCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBEQUB
   INTERFACE
      SUBROUTINE ZGBEQUB(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(OUT) :: Rowcnd
      REAL(R8KIND) , INTENT(OUT) :: Colcnd
      REAL(R8KIND) , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGBEQUB
   END INTERFACE
END MODULE S_ZGBEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBEQU
   INTERFACE
      SUBROUTINE ZGBEQU(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(OUT) :: Rowcnd
      REAL(R8KIND) , INTENT(OUT) :: Colcnd
      REAL(R8KIND) , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGBEQU
   END INTERFACE
END MODULE S_ZGBEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBRFS
   INTERFACE
      SUBROUTINE ZGBRFS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv,B,Ldb,&
     &                  X,Ldx,Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , THREE = 3.0D+0
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGBRFS
   END INTERFACE
END MODULE S_ZGBRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBRFSX
   INTERFACE
      SUBROUTINE ZGBRFSX(Trans,Equed,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,    &
     &                   Ipiv,R,C,B,Ldb,X,Ldx,Rcond,Berr,N_err_bnds,    &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 ,                     &
     &                              ITREF_DEFAULT = 1.0D+0 ,            &
     &                              ITHRESH_DEFAULT = 10.0D+0 ,         &
     &                              COMPONENTWISE_DEFAULT = 1.0D+0 ,    &
     &                              RTHRESH_DEFAULT = 0.5D+0 ,          &
     &                              DZTHRESH_DEFAULT = 0.25D+0
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
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: R
      REAL(R8KIND) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Params
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGBRFSX
   END INTERFACE
END MODULE S_ZGBRFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBSV
   INTERFACE
      SUBROUTINE ZGBSV(N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGBSV
   END INTERFACE
END MODULE S_ZGBSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBSVX
   INTERFACE
      SUBROUTINE ZGBSVX(Fact,Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv, &
     &                  Equed,R,C,B,Ldb,X,Ldx,Rcond,Ferr,Berr,Work,     &
     &                  Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: R
      REAL(R8KIND) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGBSVX
   END INTERFACE
END MODULE S_ZGBSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBSVXX
   INTERFACE
      SUBROUTINE ZGBSVXX(Fact,Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv,&
     &                   Equed,R,C,B,Ldb,X,Ldx,Rcond,Rpvgrw,Berr,       &
     &                   N_err_bnds,Err_bnds_norm,Err_bnds_comp,Nparams,&
     &                   Params,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , INTENT(OUT) :: Rpvgrw
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL(R8KIND) , DIMENSION(*) :: Params
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGBSVXX
   END INTERFACE
END MODULE S_ZGBSVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBTF2
   INTERFACE
      SUBROUTINE ZGBTF2(M,N,Kl,Ku,Ab,Ldab,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGBTF2
   END INTERFACE
END MODULE S_ZGBTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBTRF
   INTERFACE
      SUBROUTINE ZGBTRF(M,N,Kl,Ku,Ab,Ldab,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDWORK = NBMAX + 1
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGBTRF
   END INTERFACE
END MODULE S_ZGBTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGBTRS
   INTERFACE
      SUBROUTINE ZGBTRS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGBTRS
   END INTERFACE
END MODULE S_ZGBTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEBAK
   INTERFACE
      SUBROUTINE ZGEBAK(Job,Side,N,Ilo,Ihi,Scale,M,V,Ldv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Job
      CHARACTER :: Side
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Scale
      INTEGER :: M
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEBAK
   END INTERFACE
END MODULE S_ZGEBAK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEBAL
   INTERFACE
      SUBROUTINE ZGEBAL(Job,N,A,Lda,Ilo,Ihi,Scale,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              SCLFAC = 2.0D+0 , FACTOR = 0.95D+0
      CHARACTER :: Job
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) :: Ilo
      INTEGER , INTENT(OUT) :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Scale
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEBAL
   END INTERFACE
END MODULE S_ZGEBAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEBD2
   INTERFACE
      SUBROUTINE ZGEBD2(M,N,A,Lda,D,E,Tauq,Taup,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Tauq
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEBD2
   END INTERFACE
END MODULE S_ZGEBD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEBRD
   INTERFACE
      SUBROUTINE ZGEBRD(M,N,A,Lda,D,E,Tauq,Taup,Work,Lwork,Info)
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
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tauq
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEBRD
   END INTERFACE
END MODULE S_ZGEBRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGECON
   INTERFACE
      SUBROUTINE ZGECON(Norm,N,A,Lda,Anorm,Rcond,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Norm
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGECON
   END INTERFACE
END MODULE S_ZGECON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEEQUB
   INTERFACE
      SUBROUTINE ZGEEQUB(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(OUT) :: Rowcnd
      REAL(R8KIND) , INTENT(OUT) :: Colcnd
      REAL(R8KIND) , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEEQUB
   END INTERFACE
END MODULE S_ZGEEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEEQU
   INTERFACE
      SUBROUTINE ZGEEQU(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(OUT) :: Rowcnd
      REAL(R8KIND) , INTENT(OUT) :: Colcnd
      REAL(R8KIND) , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEEQU
   END INTERFACE
END MODULE S_ZGEEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEES
   INTERFACE
      SUBROUTINE ZGEES(Jobvs,Sort,SELECT,N,A,Lda,Sdim,W,Vs,Ldvs,Work,   &
     &                 Lwork,Rwork,Bwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobvs
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELECT
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER :: Sdim
      COMPLEX(CX16KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldvs,*) :: Vs
      INTEGER :: Ldvs
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEES
   END INTERFACE
END MODULE S_ZGEES
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEESX
   INTERFACE
      SUBROUTINE ZGEESX(Jobvs,Sort,SELECT,Sense,N,A,Lda,Sdim,W,Vs,Ldvs, &
     &                  Rconde,Rcondv,Work,Lwork,Rwork,Bwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobvs
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELECT
      CHARACTER :: Sense
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Sdim
      COMPLEX(CX16KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldvs,*) :: Vs
      INTEGER :: Ldvs
      REAL(R8KIND) :: Rconde
      REAL(R8KIND) , INTENT(INOUT) :: Rcondv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEESX
   END INTERFACE
END MODULE S_ZGEESX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEEV
   INTERFACE
      SUBROUTINE ZGEEV(Jobvl,Jobvr,N,A,Lda,W,Vl,Ldvl,Vr,Ldvr,Work,Lwork,&
     &                 Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEEV
   END INTERFACE
END MODULE S_ZGEEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEEVX
   INTERFACE
      SUBROUTINE ZGEEVX(Balanc,Jobvl,Jobvr,Sense,N,A,Lda,W,Vl,Ldvl,Vr,  &
     &                  Ldvr,Ilo,Ihi,Scale,Abnrm,Rconde,Rcondv,Work,    &
     &                  Lwork,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Balanc
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      CHARACTER :: Sense
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , DIMENSION(*) :: Scale
      REAL(R8KIND) , INTENT(INOUT) :: Abnrm
      REAL(R8KIND) , DIMENSION(*) :: Rconde
      REAL(R8KIND) , DIMENSION(*) :: Rcondv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEEVX
   END INTERFACE
END MODULE S_ZGEEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEHD2
   INTERFACE
      SUBROUTINE ZGEHD2(N,Ilo,Ihi,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER :: Ihi
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEHD2
   END INTERFACE
END MODULE S_ZGEHD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEHRD
   INTERFACE
      SUBROUTINE ZGEHRD(N,Ilo,Ihi,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDT = NBMAX + 1 ,           &
     &                         TSIZE = LDT*NBMAX
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEHRD
   END INTERFACE
END MODULE S_ZGEHRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEJSV
   INTERFACE
      SUBROUTINE ZGEJSV(Joba,Jobu,Jobv,Jobr,Jobt,Jobp,M,N,A,Lda,Sva,U,  &
     &                  Ldu,V,Ldv,Cwork,Lwork,Rwork,Lrwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      CHARACTER(1) :: Joba
      CHARACTER(1) :: Jobu
      CHARACTER(1) :: Jobv
      CHARACTER(1) :: Jobr
      CHARACTER(1) :: Jobt
      CHARACTER(1) :: Jobp
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N) :: Sva
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lwork) :: Cwork
      INTEGER :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lrwork) :: Rwork
      INTEGER :: Lrwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEJSV
   END INTERFACE
END MODULE S_ZGEJSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGELQ2
   INTERFACE
      SUBROUTINE ZGELQ2(M,N,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGELQ2
   END INTERFACE
END MODULE S_ZGELQ2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGELQ
   INTERFACE
      SUBROUTINE ZGELQ(M,N,A,Lda,T,Tsize,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGELQ
   END INTERFACE
END MODULE S_ZGELQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGELQF
   INTERFACE
      SUBROUTINE ZGELQF(M,N,A,Lda,Tau,Work,Lwork,Info)
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
      END SUBROUTINE ZGELQF
   END INTERFACE
END MODULE S_ZGELQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGELQT3
   INTERFACE
      RECURSIVE SUBROUTINE ZGELQT3(M,N,A,Lda,T,Ldt,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+00,0.0D+00) ,      &
     &                 ZERO = (0.0D+00,0.0D+00)
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGELQT3
   END INTERFACE
END MODULE S_ZGELQT3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGELQT
   INTERFACE
      SUBROUTINE ZGELQT(M,N,Mb,A,Lda,T,Ldt,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Mb
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGELQT
   END INTERFACE
END MODULE S_ZGELQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGELSD
   INTERFACE
      SUBROUTINE ZGELSD(M,N,Nrhs,A,Lda,B,Ldb,S,Rcond,Rank,Work,Lwork,   &
     &                  Rwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) :: Rcond
      INTEGER :: Rank
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGELSD
   END INTERFACE
END MODULE S_ZGELSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGELS
   INTERFACE
      SUBROUTINE ZGELS(Trans,M,N,Nrhs,A,Lda,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGELS
   END INTERFACE
END MODULE S_ZGELS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGELSS
   INTERFACE
      SUBROUTINE ZGELSS(M,N,Nrhs,A,Lda,B,Ldb,S,Rcond,Rank,Work,Lwork,   &
     &                  Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGELSS
   END INTERFACE
END MODULE S_ZGELSS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGELSY
   INTERFACE
      SUBROUTINE ZGELSY(M,N,Nrhs,A,Lda,B,Ldb,Jpvt,Rcond,Rank,Work,Lwork,&
     &                  Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  IMAX = 1 , IMIN = 2
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
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
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGELSY
   END INTERFACE
END MODULE S_ZGELSY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEMLQ
   INTERFACE
      SUBROUTINE ZGEMLQ(Side,Trans,M,N,K,A,Lda,T,Tsize,C,Ldc,Work,Lwork,&
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEMLQ
   END INTERFACE
END MODULE S_ZGEMLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEMLQT
   INTERFACE
      SUBROUTINE ZGEMLQT(Side,Trans,M,N,K,Mb,V,Ldv,T,Ldt,C,Ldc,Work,    &
     &                   Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Mb
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEMLQT
   END INTERFACE
END MODULE S_ZGEMLQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEMQR
   INTERFACE
      SUBROUTINE ZGEMQR(Side,Trans,M,N,K,A,Lda,T,Tsize,C,Ldc,Work,Lwork,&
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEMQR
   END INTERFACE
END MODULE S_ZGEMQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEMQRT
   INTERFACE
      SUBROUTINE ZGEMQRT(Side,Trans,M,N,K,Nb,V,Ldv,T,Ldt,C,Ldc,Work,    &
     &                   Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEMQRT
   END INTERFACE
END MODULE S_ZGEMQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEQL2
   INTERFACE
      SUBROUTINE ZGEQL2(M,N,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEQL2
   END INTERFACE
END MODULE S_ZGEQL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEQLF
   INTERFACE
      SUBROUTINE ZGEQLF(M,N,A,Lda,Tau,Work,Lwork,Info)
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
      END SUBROUTINE ZGEQLF
   END INTERFACE
END MODULE S_ZGEQLF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEQP3
   INTERFACE
      SUBROUTINE ZGEQP3(M,N,A,Lda,Jpvt,Tau,Work,Lwork,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  INB = 1 , INBMIN = 2 , IXOVER = 3
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEQP3
   END INTERFACE
END MODULE S_ZGEQP3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEQR2
   INTERFACE
      SUBROUTINE ZGEQR2(M,N,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEQR2
   END INTERFACE
END MODULE S_ZGEQR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEQR2P
   INTERFACE
      SUBROUTINE ZGEQR2P(M,N,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEQR2P
   END INTERFACE
END MODULE S_ZGEQR2P
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEQR
   INTERFACE
      SUBROUTINE ZGEQR(M,N,A,Lda,T,Tsize,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEQR
   END INTERFACE
END MODULE S_ZGEQR
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
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEQRFP
   INTERFACE
      SUBROUTINE ZGEQRFP(M,N,A,Lda,Tau,Work,Lwork,Info)
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
      END SUBROUTINE ZGEQRFP
   END INTERFACE
END MODULE S_ZGEQRFP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEQRT2
   INTERFACE
      SUBROUTINE ZGEQRT2(M,N,A,Lda,T,Ldt,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+00,0.0D+00) ,      &
     &                 ZERO = (0.0D+00,0.0D+00)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEQRT2
   END INTERFACE
END MODULE S_ZGEQRT2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEQRT3
   INTERFACE
      RECURSIVE SUBROUTINE ZGEQRT3(M,N,A,Lda,T,Ldt,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+00,0.0D+00)
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEQRT3
   END INTERFACE
END MODULE S_ZGEQRT3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGEQRT
   INTERFACE
      SUBROUTINE ZGEQRT(M,N,Nb,A,Lda,T,Ldt,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      LOGICAL , PARAMETER  ::  USE_RECURSIVE_QR = .TRUE.
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGEQRT
   END INTERFACE
END MODULE S_ZGEQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGERFS
   INTERFACE
      SUBROUTINE ZGERFS(Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,    &
     &                  Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , THREE = 3.0D+0
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGERFS
   END INTERFACE
END MODULE S_ZGERFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGERFSX
   INTERFACE
      SUBROUTINE ZGERFSX(Trans,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,R,C,B,   &
     &                   Ldb,X,Ldx,Rcond,Berr,N_err_bnds,Err_bnds_norm, &
     &                   Err_bnds_comp,Nparams,Params,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 ,                     &
     &                              ITREF_DEFAULT = 1.0D+0 ,            &
     &                              ITHRESH_DEFAULT = 10.0D+0 ,         &
     &                              COMPONENTWISE_DEFAULT = 1.0D+0 ,    &
     &                              RTHRESH_DEFAULT = 0.5D+0 ,          &
     &                              DZTHRESH_DEFAULT = 0.25D+0
      INTEGER , PARAMETER  ::  LA_LINRX_ITREF_I = 1 ,                   &
     &                         LA_LINRX_ITHRESH_I = 2 ,                 &
     &                         LA_LINRX_CWISE_I = 3 ,                   &
     &                         LA_LINRX_TRUST_I = 1 ,                   &
     &                         LA_LINRX_ERR_I = 2 , LA_LINRX_RCOND_I = 3
      CHARACTER :: Trans
      CHARACTER :: Equed
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: R
      REAL(R8KIND) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Params
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGERFSX
   END INTERFACE
END MODULE S_ZGERFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGERQ2
   INTERFACE
      SUBROUTINE ZGERQ2(M,N,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGERQ2
   END INTERFACE
END MODULE S_ZGERQ2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGERQF
   INTERFACE
      SUBROUTINE ZGERQF(M,N,A,Lda,Tau,Work,Lwork,Info)
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
      END SUBROUTINE ZGERQF
   END INTERFACE
END MODULE S_ZGERQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGESC2
   INTERFACE
      SUBROUTINE ZGESC2(N,A,Lda,Rhs,Ipiv,Jpiv,Scale)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Rhs
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Jpiv
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      END SUBROUTINE ZGESC2
   END INTERFACE
END MODULE S_ZGESC2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGESDD
   INTERFACE
      SUBROUTINE ZGESDD(Jobz,M,N,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,Rwork,&
     &                  Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobz
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX(CX16KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGESDD
   END INTERFACE
END MODULE S_ZGESDD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGESVD
   INTERFACE
      SUBROUTINE ZGESVD(Jobu,Jobvt,M,N,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,&
     &                  Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobu
      CHARACTER :: Jobvt
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX(CX16KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGESVD
   END INTERFACE
END MODULE S_ZGESVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGESVDQ
   INTERFACE
      SUBROUTINE ZGESVDQ(Joba,Jobp,Jobr,Jobu,Jobv,M,N,A,Lda,S,U,Ldu,V,  &
     &                   Ldv,Numrank,Iwork,Liwork,Cwork,Lcwork,Rwork,   &
     &                   Lrwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      CHARACTER :: Joba
      CHARACTER :: Jobp
      CHARACTER :: Jobr
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(OUT) :: Numrank
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      COMPLEX(CX16KIND) , DIMENSION(*) :: Cwork
      INTEGER :: Lcwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGESVDQ
   END INTERFACE
END MODULE S_ZGESVDQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGESVDX
   INTERFACE
      SUBROUTINE ZGESVDX(Jobu,Jobvt,Range,M,N,A,Lda,Vl,Vu,Il,Iu,Ns,S,U, &
     &                   Ldu,Vt,Ldvt,Work,Lwork,Rwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobu
      CHARACTER :: Jobvt
      CHARACTER :: Range
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      INTEGER , INTENT(INOUT) :: Ns
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX(CX16KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGESVDX
   END INTERFACE
END MODULE S_ZGESVDX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGESV
   INTERFACE
      SUBROUTINE ZGESV(N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGESV
   END INTERFACE
END MODULE S_ZGESV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGESVJ
   INTERFACE
      SUBROUTINE ZGESVJ(Joba,Jobu,Jobv,M,N,A,Lda,Sva,Mv,V,Ldv,Cwork,    &
     &                  Lwork,Rwork,Lrwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , HALF = 0.5D0 ,       &
     &                              ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      INTEGER , PARAMETER  ::  NSWEEP = 30
      CHARACTER(1) :: Joba
      CHARACTER(1) :: Jobu
      CHARACTER(1) :: Jobv
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lwork) :: Cwork
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lrwork) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGESVJ
   END INTERFACE
END MODULE S_ZGESVJ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGESVX
   INTERFACE
      SUBROUTINE ZGESVX(Fact,Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,R,C, &
     &                  B,Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: R
      REAL(R8KIND) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGESVX
   END INTERFACE
END MODULE S_ZGESVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGESVXX
   INTERFACE
      SUBROUTINE ZGESVXX(Fact,Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,R,C,&
     &                   B,Ldb,X,Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,      &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , INTENT(OUT) :: Rpvgrw
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL(R8KIND) , DIMENSION(*) :: Params
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGESVXX
   END INTERFACE
END MODULE S_ZGESVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGETC2
   INTERFACE
      SUBROUTINE ZGETC2(N,A,Lda,Ipiv,Jpiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Jpiv
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE ZGETC2
   END INTERFACE
END MODULE S_ZGETC2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGETF2
   INTERFACE
      SUBROUTINE ZGETF2(M,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGETF2
   END INTERFACE
END MODULE S_ZGETF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGETRF2
   INTERFACE
      RECURSIVE SUBROUTINE ZGETRF2(M,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGETRF2
   END INTERFACE
END MODULE S_ZGETRF2
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
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGETRI
   INTERFACE
      SUBROUTINE ZGETRI(N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGETRI
   END INTERFACE
END MODULE S_ZGETRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGETRS
   INTERFACE
      SUBROUTINE ZGETRS(Trans,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGETRS
   END INTERFACE
END MODULE S_ZGETRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGETSLS
   INTERFACE
      SUBROUTINE ZGETSLS(Trans,M,N,Nrhs,A,Lda,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGETSLS
   END INTERFACE
END MODULE S_ZGETSLS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGETSQRHRT
   INTERFACE
      SUBROUTINE ZGETSQRHRT(M,N,Mb1,Nb1,Nb2,A,Lda,T,Ldt,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb1
      INTEGER , INTENT(IN) :: Nb1
      INTEGER , INTENT(IN) :: Nb2
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGETSQRHRT
   END INTERFACE
END MODULE S_ZGETSQRHRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGBAK
   INTERFACE
      SUBROUTINE ZGGBAK(Job,Side,N,Ilo,Ihi,Lscale,Rscale,M,V,Ldv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Job
      CHARACTER :: Side
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL(R8KIND) , DIMENSION(*) :: Lscale
      REAL(R8KIND) , DIMENSION(*) :: Rscale
      INTEGER :: M
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGBAK
   END INTERFACE
END MODULE S_ZGGBAK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGBAL
   INTERFACE
      SUBROUTINE ZGGBAL(Job,N,A,Lda,B,Ldb,Ilo,Ihi,Lscale,Rscale,Work,   &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , HALF = 0.5D+0 ,     &
     &                              ONE = 1.0D+0 , THREE = 3.0D+0 ,     &
     &                              SCLFAC = 1.0D+1
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Job
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Ilo
      INTEGER , INTENT(INOUT) :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Lscale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rscale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGBAL
   END INTERFACE
END MODULE S_ZGGBAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGES3
   INTERFACE
      SUBROUTINE ZGGES3(Jobvsl,Jobvsr,Sort,SELCTG,N,A,Lda,B,Ldb,Sdim,   &
     &                  Alpha,Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work,Lwork,Rwork,&
     &                  Bwork,Info)
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
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      COMPLEX(CX16KIND) , DIMENSION(*) :: Alpha
      COMPLEX(CX16KIND) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      COMPLEX(CX16KIND) , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGES3
   END INTERFACE
END MODULE S_ZGGES3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGES
   INTERFACE
      SUBROUTINE ZGGES(Jobvsl,Jobvsr,Sort,SELCTG,N,A,Lda,B,Ldb,Sdim,    &
     &                 Alpha,Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work,Lwork,Rwork, &
     &                 Bwork,Info)
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
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      COMPLEX(CX16KIND) , DIMENSION(*) :: Alpha
      COMPLEX(CX16KIND) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      COMPLEX(CX16KIND) , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGES
   END INTERFACE
END MODULE S_ZGGES
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGESX
   INTERFACE
      SUBROUTINE ZGGESX(Jobvsl,Jobvsr,Sort,SELCTG,Sense,N,A,Lda,B,Ldb,  &
     &                  Sdim,Alpha,Beta,Vsl,Ldvsl,Vsr,Ldvsr,Rconde,     &
     &                  Rcondv,Work,Lwork,Rwork,Iwork,Liwork,Bwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      CHARACTER :: Sense
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      COMPLEX(CX16KIND) , DIMENSION(*) :: Alpha
      COMPLEX(CX16KIND) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      COMPLEX(CX16KIND) , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(2) :: Rconde
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(2) :: Rcondv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGESX
   END INTERFACE
END MODULE S_ZGGESX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGEV3
   INTERFACE
      SUBROUTINE ZGGEV3(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alpha,Beta,Vl,Ldvl,Vr,&
     &                  Ldvr,Work,Lwork,Rwork,Info)
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
      COMPLEX(CX16KIND) , DIMENSION(*) :: Alpha
      COMPLEX(CX16KIND) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGEV3
   END INTERFACE
END MODULE S_ZGGEV3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGEV
   INTERFACE
      SUBROUTINE ZGGEV(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alpha,Beta,Vl,Ldvl,Vr, &
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
      COMPLEX(CX16KIND) , DIMENSION(*) :: Alpha
      COMPLEX(CX16KIND) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGEV
   END INTERFACE
END MODULE S_ZGGEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGEVX
   INTERFACE
      SUBROUTINE ZGGEVX(Balanc,Jobvl,Jobvr,Sense,N,A,Lda,B,Ldb,Alpha,   &
     &                  Beta,Vl,Ldvl,Vr,Ldvr,Ilo,Ihi,Lscale,Rscale,     &
     &                  Abnrm,Bbnrm,Rconde,Rcondv,Work,Lwork,Rwork,     &
     &                  Iwork,Bwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Balanc
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      CHARACTER :: Sense
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Alpha
      COMPLEX(CX16KIND) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , DIMENSION(*) :: Lscale
      REAL(R8KIND) , DIMENSION(*) :: Rscale
      REAL(R8KIND) , INTENT(INOUT) :: Abnrm
      REAL(R8KIND) , INTENT(INOUT) :: Bbnrm
      REAL(R8KIND) , DIMENSION(*) :: Rconde
      REAL(R8KIND) , DIMENSION(*) :: Rcondv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGEVX
   END INTERFACE
END MODULE S_ZGGEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGGLM
   INTERFACE
      SUBROUTINE ZGGGLM(N,M,P,A,Lda,B,Ldb,D,X,Y,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      INTEGER :: N
      INTEGER :: M
      INTEGER :: P
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: X
      COMPLEX(CX16KIND) , DIMENSION(*) :: Y
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGGLM
   END INTERFACE
END MODULE S_ZGGGLM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGHD3
   INTERFACE
      SUBROUTINE ZGGHD3(Compq,Compz,N,Ilo,Ihi,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,  &
     &                  Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGHD3
   END INTERFACE
END MODULE S_ZGGHD3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGHRD
   INTERFACE
      SUBROUTINE ZGGHRD(Compq,Compz,N,Ilo,Ihi,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,  &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER :: Ihi
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGHRD
   END INTERFACE
END MODULE S_ZGGHRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGLSE
   INTERFACE
      SUBROUTINE ZGGLSE(M,N,P,A,Lda,B,Ldb,C,D,X,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: P
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: X
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGLSE
   END INTERFACE
END MODULE S_ZGGLSE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGQRF
   INTERFACE
      SUBROUTINE ZGGQRF(N,M,P,A,Lda,Taua,B,Ldb,Taub,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: M
      INTEGER :: P
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taua
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taub
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGQRF
   END INTERFACE
END MODULE S_ZGGQRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGRQF
   INTERFACE
      SUBROUTINE ZGGRQF(M,P,N,A,Lda,Taua,B,Ldb,Taub,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: P
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taua
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taub
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGRQF
   END INTERFACE
END MODULE S_ZGGRQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGSVD3
   INTERFACE
      SUBROUTINE ZGGSVD3(Jobu,Jobv,Jobq,M,N,P,K,L,A,Lda,B,Ldb,Alpha,    &
     &                   Beta,U,Ldu,V,Ldv,Q,Ldq,Work,Lwork,Rwork,Iwork, &
     &                   Info)
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
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGSVD3
   END INTERFACE
END MODULE S_ZGGSVD3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGGSVP3
   INTERFACE
      SUBROUTINE ZGGSVP3(Jobu,Jobv,Jobq,M,P,N,A,Lda,B,Ldb,Tola,Tolb,K,L,&
     &                   U,Ldu,V,Ldv,Q,Ldq,Iwork,Rwork,Tau,Work,Lwork,  &
     &                   Info)
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
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGGSVP3
   END INTERFACE
END MODULE S_ZGGSVP3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGSVJ0
   INTERFACE
      SUBROUTINE ZGSVJ0(Jobv,M,N,A,Lda,D,Sva,Mv,V,Ldv,Eps,Sfmin,Tol,    &
     &                  Nsweep,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , HALF = 0.5D0 ,       &
     &                              ONE = 1.0D0
      CHARACTER(1) :: Jobv
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(N) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      REAL(R8KIND) , INTENT(IN) :: Eps
      REAL(R8KIND) , INTENT(IN) :: Sfmin
      REAL(R8KIND) , INTENT(IN) :: Tol
      INTEGER , INTENT(IN) :: Nsweep
      COMPLEX(CX16KIND) , DIMENSION(Lwork) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGSVJ0
   END INTERFACE
END MODULE S_ZGSVJ0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGSVJ1
   INTERFACE
      SUBROUTINE ZGSVJ1(Jobv,M,N,N1,A,Lda,D,Sva,Mv,V,Ldv,Eps,Sfmin,Tol, &
     &                  Nsweep,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , HALF = 0.5D0 ,       &
     &                              ONE = 1.0D0
      CHARACTER(1) :: Jobv
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: N1
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(N) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      REAL(R8KIND) , INTENT(IN) :: Eps
      REAL(R8KIND) , INTENT(IN) :: Sfmin
      REAL(R8KIND) , INTENT(IN) :: Tol
      INTEGER , INTENT(IN) :: Nsweep
      COMPLEX(CX16KIND) , DIMENSION(Lwork) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGSVJ1
   END INTERFACE
END MODULE S_ZGSVJ1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGTCON
   INTERFACE
      SUBROUTINE ZGTCON(Norm,N,Dl,D,Du,Du2,Ipiv,Anorm,Rcond,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Norm
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dl
      COMPLEX(CX16KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: Du
      COMPLEX(CX16KIND) , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGTCON
   END INTERFACE
END MODULE S_ZGTCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGTRFS
   INTERFACE
      SUBROUTINE ZGTRFS(Trans,N,Nrhs,Dl,D,Du,Dlf,Df,Duf,Du2,Ipiv,B,Ldb, &
     &                  X,Ldx,Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0 , THREE = 3.0D+0
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dl
      COMPLEX(CX16KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: Du
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dlf
      COMPLEX(CX16KIND) , DIMENSION(*) :: Df
      COMPLEX(CX16KIND) , DIMENSION(*) :: Duf
      COMPLEX(CX16KIND) , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGTRFS
   END INTERFACE
END MODULE S_ZGTRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGTSV
   INTERFACE
      SUBROUTINE ZGTSV(N,Nrhs,Dl,D,Du,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Dl
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Du
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGTSV
   END INTERFACE
END MODULE S_ZGTSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGTSVX
   INTERFACE
      SUBROUTINE ZGTSVX(Fact,Trans,N,Nrhs,Dl,D,Du,Dlf,Df,Duf,Du2,Ipiv,B,&
     &                  Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dl
      COMPLEX(CX16KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: Du
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dlf
      COMPLEX(CX16KIND) , DIMENSION(*) :: Df
      COMPLEX(CX16KIND) , DIMENSION(*) :: Duf
      COMPLEX(CX16KIND) , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGTSVX
   END INTERFACE
END MODULE S_ZGTSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGTTRF
   INTERFACE
      SUBROUTINE ZGTTRF(N,Dl,D,Du,Du2,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Dl
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Du
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Du2
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER :: Info
      END SUBROUTINE ZGTTRF
   END INTERFACE
END MODULE S_ZGTTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGTTRS
   INTERFACE
      SUBROUTINE ZGTTRS(Trans,N,Nrhs,Dl,D,Du,Du2,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dl
      COMPLEX(CX16KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: Du
      COMPLEX(CX16KIND) , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZGTTRS
   END INTERFACE
END MODULE S_ZGTTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZGTTS2
   INTERFACE
      SUBROUTINE ZGTTS2(Itrans,N,Nrhs,Dl,D,Du,Du2,Ipiv,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: Itrans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Dl
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Du
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Du2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE ZGTTS2
   END INTERFACE
END MODULE S_ZGTTS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHB2ST_KERNELS
   INTERFACE
      SUBROUTINE ZHB2ST_KERNELS(Uplo,Wantz,Ttype,St,Ed,Sweep,N,Nb,Ib,A, &
     &                          Lda,V,Tau,Ldvt,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: Ttype
      INTEGER , INTENT(IN) :: St
      INTEGER , INTENT(IN) :: Ed
      INTEGER , INTENT(IN) :: Sweep
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(IN) :: Ib
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: V
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      INTEGER , INTENT(IN) :: Ldvt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      END SUBROUTINE ZHB2ST_KERNELS
   END INTERFACE
END MODULE S_ZHB2ST_KERNELS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHBEV_2STAGE
   INTERFACE
      SUBROUTINE ZHBEV_2STAGE(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Lwork,&
     &                        Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHBEV_2STAGE
   END INTERFACE
END MODULE S_ZHBEV_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHBEVD_2STAGE
   INTERFACE
      SUBROUTINE ZHBEVD_2STAGE(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,     &
     &                         Lwork,Rwork,Lrwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHBEVD_2STAGE
   END INTERFACE
END MODULE S_ZHBEVD_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHBEVD
   INTERFACE
      SUBROUTINE ZHBEVD(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Lwork,Rwork,&
     &                  Lrwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHBEVD
   END INTERFACE
END MODULE S_ZHBEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHBEV
   INTERFACE
      SUBROUTINE ZHBEV(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHBEV
   END INTERFACE
END MODULE S_ZHBEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHBEVX_2STAGE
   INTERFACE
      SUBROUTINE ZHBEVX_2STAGE(Jobz,Range,Uplo,N,Kd,Ab,Ldab,Q,Ldq,Vl,Vu,&
     &                         Il,Iu,Abstol,M,W,Z,Ldz,Work,Lwork,Rwork, &
     &                         Iwork,Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHBEVX_2STAGE
   END INTERFACE
END MODULE S_ZHBEVX_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHBEVX
   INTERFACE
      SUBROUTINE ZHBEVX(Jobz,Range,Uplo,N,Kd,Ab,Ldab,Q,Ldq,Vl,Vu,Il,Iu, &
     &                  Abstol,M,W,Z,Ldz,Work,Rwork,Iwork,Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHBEVX
   END INTERFACE
END MODULE S_ZHBEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHBGST
   INTERFACE
      SUBROUTINE ZHBGST(Vect,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,X,Ldx,Work,   &
     &                  Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ka
      INTEGER , INTENT(IN) :: Kb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldbb,*) :: Bb
      INTEGER , INTENT(IN) :: Ldbb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHBGST
   END INTERFACE
END MODULE S_ZHBGST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHBGVD
   INTERFACE
      SUBROUTINE ZHBGVD(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,W,Z,Ldz,Work, &
     &                  Lwork,Rwork,Lrwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ka
      INTEGER :: Kb
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldbb,*) :: Bb
      INTEGER :: Ldbb
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHBGVD
   END INTERFACE
END MODULE S_ZHBGVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHBGV
   INTERFACE
      SUBROUTINE ZHBGV(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,W,Z,Ldz,Work,  &
     &                 Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ka
      INTEGER :: Kb
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldbb,*) :: Bb
      INTEGER :: Ldbb
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHBGV
   END INTERFACE
END MODULE S_ZHBGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHBGVX
   INTERFACE
      SUBROUTINE ZHBGVX(Jobz,Range,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,Q,Ldq,  &
     &                  Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,Work,Rwork,Iwork,  &
     &                  Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ka
      INTEGER :: Kb
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldbb,*) :: Bb
      INTEGER :: Ldbb
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHBGVX
   END INTERFACE
END MODULE S_ZHBGVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHBTRD
   INTERFACE
      SUBROUTINE ZHBTRD(Vect,Uplo,N,Kd,Ab,Ldab,D,E,Q,Ldq,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHBTRD
   END INTERFACE
END MODULE S_ZHBTRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHECON_3
   INTERFACE
      SUBROUTINE ZHECON_3(Uplo,N,A,Lda,E,Ipiv,Anorm,Rcond,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHECON_3
   END INTERFACE
END MODULE S_ZHECON_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHECON
   INTERFACE
      SUBROUTINE ZHECON(Uplo,N,A,Lda,Ipiv,Anorm,Rcond,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHECON
   END INTERFACE
END MODULE S_ZHECON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHECON_ROOK
   INTERFACE
      SUBROUTINE ZHECON_ROOK(Uplo,N,A,Lda,Ipiv,Anorm,Rcond,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHECON_ROOK
   END INTERFACE
END MODULE S_ZHECON_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEEQUB
   INTERFACE
      SUBROUTINE ZHEEQUB(Uplo,N,A,Lda,S,Scond,Amax,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0
      INTEGER , PARAMETER  ::  MAX_ITER = 100
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) :: Scond
      REAL(R8KIND) , INTENT(INOUT) :: Amax
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEEQUB
   END INTERFACE
END MODULE S_ZHEEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEEV_2STAGE
   INTERFACE
      SUBROUTINE ZHEEV_2STAGE(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D0,0.0D0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEEV_2STAGE
   END INTERFACE
END MODULE S_ZHEEV_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEEVD_2STAGE
   INTERFACE
      SUBROUTINE ZHEEVD_2STAGE(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Rwork,    &
     &                         Lrwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D0,0.0D0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEEVD_2STAGE
   END INTERFACE
END MODULE S_ZHEEVD_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEEVD
   INTERFACE
      SUBROUTINE ZHEEVD(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Rwork,Lrwork,    &
     &                  Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D0,0.0D0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEEVD
   END INTERFACE
END MODULE S_ZHEEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEEV
   INTERFACE
      SUBROUTINE ZHEEV(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D0,0.0D0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEEV
   END INTERFACE
END MODULE S_ZHEEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEEVR_2STAGE
   INTERFACE
      SUBROUTINE ZHEEVR_2STAGE(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,     &
     &                         Abstol,M,W,Z,Ldz,Isuppz,Work,Lwork,Rwork,&
     &                         Lrwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEEVR_2STAGE
   END INTERFACE
END MODULE S_ZHEEVR_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEEVR
   INTERFACE
      SUBROUTINE ZHEEVR(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,Abstol,M,W, &
     &                  Z,Ldz,Isuppz,Work,Lwork,Rwork,Lrwork,Iwork,     &
     &                  Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEEVR
   END INTERFACE
END MODULE S_ZHEEVR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEEVX_2STAGE
   INTERFACE
      SUBROUTINE ZHEEVX_2STAGE(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,     &
     &                         Abstol,M,W,Z,Ldz,Work,Lwork,Rwork,Iwork, &
     &                         Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEEVX_2STAGE
   END INTERFACE
END MODULE S_ZHEEVX_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEEVX
   INTERFACE
      SUBROUTINE ZHEEVX(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,Abstol,M,W, &
     &                  Z,Ldz,Work,Lwork,Rwork,Iwork,Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEEVX
   END INTERFACE
END MODULE S_ZHEEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEGS2
   INTERFACE
      SUBROUTINE ZHEGS2(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , HALF = 0.5D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: Itype
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEGS2
   END INTERFACE
END MODULE S_ZHEGS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEGST
   INTERFACE
      SUBROUTINE ZHEGST(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 HALF = (0.5D+0,0.0D+0)
      INTEGER :: Itype
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEGST
   END INTERFACE
END MODULE S_ZHEGST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEGV_2STAGE
   INTERFACE
      SUBROUTINE ZHEGV_2STAGE(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,     &
     &                        Lwork,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEGV_2STAGE
   END INTERFACE
END MODULE S_ZHEGV_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEGVD
   INTERFACE
      SUBROUTINE ZHEGVD(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,Lwork,     &
     &                  Rwork,Lrwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER :: Lrwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEGVD
   END INTERFACE
END MODULE S_ZHEGVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEGV
   INTERFACE
      SUBROUTINE ZHEGV(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,Lwork,Rwork,&
     &                 Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEGV
   END INTERFACE
END MODULE S_ZHEGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHEGVX
   INTERFACE
      SUBROUTINE ZHEGVX(Itype,Jobz,Range,Uplo,N,A,Lda,B,Ldb,Vl,Vu,Il,Iu,&
     &                  Abstol,M,W,Z,Ldz,Work,Lwork,Rwork,Iwork,Ifail,  &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) :: Abstol
      INTEGER :: M
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHEGVX
   END INTERFACE
END MODULE S_ZHEGVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHERFS
   INTERFACE
      SUBROUTINE ZHERFS(Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,Ferr,&
     &                  Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , THREE = 3.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHERFS
   END INTERFACE
END MODULE S_ZHERFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHERFSX
   INTERFACE
      SUBROUTINE ZHERFSX(Uplo,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,S,B,Ldb,X,&
     &                   Ldx,Rcond,Berr,N_err_bnds,Err_bnds_norm,       &
     &                   Err_bnds_comp,Nparams,Params,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 ,                     &
     &                              ITREF_DEFAULT = 1.0D+0 ,            &
     &                              ITHRESH_DEFAULT = 10.0D+0 ,         &
     &                              COMPONENTWISE_DEFAULT = 1.0D+0 ,    &
     &                              RTHRESH_DEFAULT = 0.5D+0 ,          &
     &                              DZTHRESH_DEFAULT = 0.25D+0
      INTEGER , PARAMETER  ::  LA_LINRX_ITREF_I = 1 ,                   &
     &                         LA_LINRX_ITHRESH_I = 2 ,                 &
     &                         LA_LINRX_CWISE_I = 3 ,                   &
     &                         LA_LINRX_TRUST_I = 1 ,                   &
     &                         LA_LINRX_ERR_I = 2 , LA_LINRX_RCOND_I = 3
      CHARACTER :: Uplo
      CHARACTER :: Equed
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Params
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHERFSX
   END INTERFACE
END MODULE S_ZHERFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHESV_AA_2STAGE
   INTERFACE
      SUBROUTINE ZHESV_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B, &
     &                           Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tb
      INTEGER :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHESV_AA_2STAGE
   END INTERFACE
END MODULE S_ZHESV_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHESV_AA
   INTERFACE
      SUBROUTINE ZHESV_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHESV_AA
   END INTERFACE
END MODULE S_ZHESV_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHESV
   INTERFACE
      SUBROUTINE ZHESV(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHESV
   END INTERFACE
END MODULE S_ZHESV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHESV_RK
   INTERFACE
      SUBROUTINE ZHESV_RK(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Work,Lwork,    &
     &                    Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHESV_RK
   END INTERFACE
END MODULE S_ZHESV_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHESV_ROOK
   INTERFACE
      SUBROUTINE ZHESV_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,    &
     &                      Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHESV_ROOK
   END INTERFACE
END MODULE S_ZHESV_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHESVX
   INTERFACE
      SUBROUTINE ZHESVX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,&
     &                  Rcond,Ferr,Berr,Work,Lwork,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHESVX
   END INTERFACE
END MODULE S_ZHESVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHESVXX
   INTERFACE
      SUBROUTINE ZHESVXX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,S,B, &
     &                   Ldb,X,Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,        &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , INTENT(OUT) :: Rpvgrw
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL(R8KIND) , DIMENSION(*) :: Params
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHESVXX
   END INTERFACE
END MODULE S_ZHESVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHESWAPR
   INTERFACE
      SUBROUTINE ZHESWAPR(Uplo,N,A,Lda,I1,I2)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,N) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) :: I1
      INTEGER , INTENT(IN) :: I2
      END SUBROUTINE ZHESWAPR
   END INTERFACE
END MODULE S_ZHESWAPR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETD2
   INTERFACE
      SUBROUTINE ZHETD2(Uplo,N,A,Lda,D,E,Tau,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0) , HALF = (0.5D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETD2
   END INTERFACE
END MODULE S_ZHETD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETF2
   INTERFACE
      SUBROUTINE ZHETF2(Uplo,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETF2
   END INTERFACE
END MODULE S_ZHETF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETF2_RK
   INTERFACE
      SUBROUTINE ZHETF2_RK(Uplo,N,A,Lda,E,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETF2_RK
   END INTERFACE
END MODULE S_ZHETF2_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETF2_ROOK
   INTERFACE
      SUBROUTINE ZHETF2_ROOK(Uplo,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETF2_ROOK
   END INTERFACE
END MODULE S_ZHETF2_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRD_2STAGE
   INTERFACE
      SUBROUTINE ZHETRD_2STAGE(Vect,Uplo,N,A,Lda,D,E,Tau,Hous2,Lhous2,  &
     &                         Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Hous2
      INTEGER :: Lhous2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRD_2STAGE
   END INTERFACE
END MODULE S_ZHETRD_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRD
   INTERFACE
      SUBROUTINE ZHETRD(Uplo,N,A,Lda,D,E,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRD
   END INTERFACE
END MODULE S_ZHETRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRD_HE2HB
   INTERFACE
      SUBROUTINE ZHETRD_HE2HB(Uplo,N,Kd,A,Lda,Ab,Ldab,Tau,Work,Lwork,   &
     &                        Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  RONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0) , HALF = (0.5D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRD_HE2HB
   END INTERFACE
END MODULE S_ZHETRD_HE2HB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRF_AA_2STAGE
   INTERFACE
      SUBROUTINE ZHETRF_AA_2STAGE(Uplo,N,A,Lda,Tb,Ltb,Ipiv,Ipiv2,Work,  &
     &                            Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0) ,       &
     &                 ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRF_AA_2STAGE
   END INTERFACE
END MODULE S_ZHETRF_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRF_AA
   INTERFACE
      SUBROUTINE ZHETRF_AA(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRF_AA
   END INTERFACE
END MODULE S_ZHETRF_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRF
   INTERFACE
      SUBROUTINE ZHETRF(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRF
   END INTERFACE
END MODULE S_ZHETRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRF_RK
   INTERFACE
      SUBROUTINE ZHETRF_RK(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRF_RK
   END INTERFACE
END MODULE S_ZHETRF_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRF_ROOK
   INTERFACE
      SUBROUTINE ZHETRF_ROOK(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRF_ROOK
   END INTERFACE
END MODULE S_ZHETRF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRI2
   INTERFACE
      SUBROUTINE ZHETRI2(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRI2
   END INTERFACE
END MODULE S_ZHETRI2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRI2X
   INTERFACE
      SUBROUTINE ZHETRI2X(Uplo,N,A,Lda,Ipiv,Work,Nb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRI2X
   END INTERFACE
END MODULE S_ZHETRI2X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRI_3
   INTERFACE
      SUBROUTINE ZHETRI_3(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRI_3
   END INTERFACE
END MODULE S_ZHETRI_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRI_3X
   INTERFACE
      SUBROUTINE ZHETRI_3X(Uplo,N,A,Lda,E,Ipiv,Work,Nb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRI_3X
   END INTERFACE
END MODULE S_ZHETRI_3X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRI
   INTERFACE
      SUBROUTINE ZHETRI(Uplo,N,A,Lda,Ipiv,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRI
   END INTERFACE
END MODULE S_ZHETRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRI_ROOK
   INTERFACE
      SUBROUTINE ZHETRI_ROOK(Uplo,N,A,Lda,Ipiv,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRI_ROOK
   END INTERFACE
END MODULE S_ZHETRI_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRS2
   INTERFACE
      SUBROUTINE ZHETRS2(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRS2
   END INTERFACE
END MODULE S_ZHETRS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRS_3
   INTERFACE
      SUBROUTINE ZHETRS_3(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRS_3
   END INTERFACE
END MODULE S_ZHETRS_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRS_AA_2STAGE
   INTERFACE
      SUBROUTINE ZHETRS_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B,&
     &                            Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRS_AA_2STAGE
   END INTERFACE
END MODULE S_ZHETRS_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRS_AA
   INTERFACE
      SUBROUTINE ZHETRS_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRS_AA
   END INTERFACE
END MODULE S_ZHETRS_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRS
   INTERFACE
      SUBROUTINE ZHETRS(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRS
   END INTERFACE
END MODULE S_ZHETRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHETRS_ROOK
   INTERFACE
      SUBROUTINE ZHETRS_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHETRS_ROOK
   END INTERFACE
END MODULE S_ZHETRS_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHFRK
   INTERFACE
      SUBROUTINE ZHFRK(Transr,Uplo,Trans,N,K,Alpha,A,Lda,Beta,C)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Transr
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) :: Alpha
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) :: Beta
      COMPLEX(CX16KIND) , DIMENSION(*) :: C
      END SUBROUTINE ZHFRK
   END INTERFACE
END MODULE S_ZHFRK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHGEQZ
   INTERFACE
      SUBROUTINE ZHGEQZ(Job,Compq,Compz,N,Ilo,Ihi,H,Ldh,T,Ldt,Alpha,    &
     &                  Beta,Q,Ldq,Z,Ldz,Work,Lwork,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              HALF = 0.5D+0
      CHARACTER :: Job
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Alpha
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHGEQZ
   END INTERFACE
END MODULE S_ZHGEQZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPCON
   INTERFACE
      SUBROUTINE ZHPCON(Uplo,N,Ap,Ipiv,Anorm,Rcond,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPCON
   END INTERFACE
END MODULE S_ZHPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPEVD
   INTERFACE
      SUBROUTINE ZHPEVD(Jobz,Uplo,N,Ap,W,Z,Ldz,Work,Lwork,Rwork,Lrwork, &
     &                  Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPEVD
   END INTERFACE
END MODULE S_ZHPEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPEV
   INTERFACE
      SUBROUTINE ZHPEV(Jobz,Uplo,N,Ap,W,Z,Ldz,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPEV
   END INTERFACE
END MODULE S_ZHPEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPEVX
   INTERFACE
      SUBROUTINE ZHPEVX(Jobz,Range,Uplo,N,Ap,Vl,Vu,Il,Iu,Abstol,M,W,Z,  &
     &                  Ldz,Work,Rwork,Iwork,Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D0,0.0D0)
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPEVX
   END INTERFACE
END MODULE S_ZHPEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPGST
   INTERFACE
      SUBROUTINE ZHPGST(Itype,Uplo,N,Ap,Bp,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , HALF = 0.5D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: Itype
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Bp
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPGST
   END INTERFACE
END MODULE S_ZHPGST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPGVD
   INTERFACE
      SUBROUTINE ZHPGVD(Itype,Jobz,Uplo,N,Ap,Bp,W,Z,Ldz,Work,Lwork,     &
     &                  Rwork,Lrwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Bp
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER :: Lrwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPGVD
   END INTERFACE
END MODULE S_ZHPGVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPGV
   INTERFACE
      SUBROUTINE ZHPGV(Itype,Jobz,Uplo,N,Ap,Bp,W,Z,Ldz,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Bp
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPGV
   END INTERFACE
END MODULE S_ZHPGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPGVX
   INTERFACE
      SUBROUTINE ZHPGVX(Itype,Jobz,Range,Uplo,N,Ap,Bp,Vl,Vu,Il,Iu,      &
     &                  Abstol,M,W,Z,Ldz,Work,Rwork,Iwork,Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Bp
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPGVX
   END INTERFACE
END MODULE S_ZHPGVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPRFS
   INTERFACE
      SUBROUTINE ZHPRFS(Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Ferr,Berr,  &
     &                  Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , THREE = 3.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPRFS
   END INTERFACE
END MODULE S_ZHPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPSV
   INTERFACE
      SUBROUTINE ZHPSV(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPSV
   END INTERFACE
END MODULE S_ZHPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPSVX
   INTERFACE
      SUBROUTINE ZHPSVX(Fact,Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Rcond, &
     &                  Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPSVX
   END INTERFACE
END MODULE S_ZHPSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPTRD
   INTERFACE
      SUBROUTINE ZHPTRD(Uplo,N,Ap,D,E,Tau,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0) , HALF = (0.5D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPTRD
   END INTERFACE
END MODULE S_ZHPTRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPTRF
   INTERFACE
      SUBROUTINE ZHPTRF(Uplo,N,Ap,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPTRF
   END INTERFACE
END MODULE S_ZHPTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPTRI
   INTERFACE
      SUBROUTINE ZHPTRI(Uplo,N,Ap,Ipiv,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPTRI
   END INTERFACE
END MODULE S_ZHPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHPTRS
   INTERFACE
      SUBROUTINE ZHPTRS(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHPTRS
   END INTERFACE
END MODULE S_ZHPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHSEIN
   INTERFACE
      SUBROUTINE ZHSEIN(Side,Eigsrc,Initv,Select,N,H,Ldh,W,Vl,Ldvl,Vr,  &
     &                  Ldvr,Mm,M,Work,Rwork,Ifaill,Ifailr,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  RZERO = 0.0D+0
      CHARACTER :: Side
      CHARACTER :: Eigsrc
      CHARACTER :: Initv
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      COMPLEX(CX16KIND) , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifaill
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifailr
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHSEIN
   END INTERFACE
END MODULE S_ZHSEIN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZHSEQR
   INTERFACE
      SUBROUTINE ZHSEQR(Job,Compz,N,Ilo,Ihi,H,Ldh,W,Z,Ldz,Work,Lwork,   &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NTINY = 15 , NL = 49
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D0,0.0D0) ,         &
     &                 ONE = (1.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  RZERO = 0.0D0
      CHARACTER :: Job
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX(CX16KIND) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX(CX16KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZHSEQR
   END INTERFACE
END MODULE S_ZHSEQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLABRD
   INTERFACE
      SUBROUTINE ZLABRD(M,N,Nb,A,Lda,D,E,Tauq,Taup,X,Ldx,Y,Ldy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tauq
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      COMPLEX(CX16KIND) , DIMENSION(Ldy,*) :: Y
      INTEGER :: Ldy
      END SUBROUTINE ZLABRD
   END INTERFACE
END MODULE S_ZLABRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLACGV
   INTERFACE
      SUBROUTINE ZLACGV(N,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE ZLACGV
   END INTERFACE
END MODULE S_ZLACGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLACN2
   INTERFACE
      SUBROUTINE ZLACN2(N,V,X,Est,Kase,Isave)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , TWO = 2.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: V
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      REAL(R8KIND) , INTENT(INOUT) :: Est
      INTEGER , INTENT(INOUT) :: Kase
      INTEGER , INTENT(INOUT) , DIMENSION(3) :: Isave
      END SUBROUTINE ZLACN2
   END INTERFACE
END MODULE S_ZLACN2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLACON
   INTERFACE
      SUBROUTINE ZLACON(N,V,X,Est,Kase)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , TWO = 2.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(N) :: V
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(N) :: X
      REAL(R8KIND) , INTENT(INOUT) :: Est
      INTEGER , INTENT(INOUT) :: Kase
      END SUBROUTINE ZLACON
   END INTERFACE
END MODULE S_ZLACON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLACP2
   INTERFACE
      SUBROUTINE ZLACP2(Uplo,M,N,A,Lda,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE ZLACP2
   END INTERFACE
END MODULE S_ZLACP2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLACPY
   INTERFACE
      SUBROUTINE ZLACPY(Uplo,M,N,A,Lda,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE ZLACPY
   END INTERFACE
END MODULE S_ZLACPY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLACRM
   INTERFACE
      SUBROUTINE ZLACRM(M,N,A,Lda,B,Ldb,C,Ldc,Rwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END SUBROUTINE ZLACRM
   END INTERFACE
END MODULE S_ZLACRM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLACRT
   INTERFACE
      SUBROUTINE ZLACRT(N,Cx,Incx,Cy,Incy,C,S)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
      COMPLEX(CX16KIND) , INTENT(IN) :: C
      COMPLEX(CX16KIND) , INTENT(IN) :: S
      END SUBROUTINE ZLACRT
   END INTERFACE
END MODULE S_ZLACRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLADIV
   INTERFACE
      FUNCTION ZLADIV(X,Y)
      USE F77KINDS                        
      IMPLICIT NONE
      COMPLEX(CX16KIND) :: ZLADIV
      COMPLEX(CX16KIND) , INTENT(IN) :: X
      COMPLEX(CX16KIND) , INTENT(IN) :: Y
      END FUNCTION ZLADIV
   END INTERFACE
END MODULE S_ZLADIV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAED0
   INTERFACE
      SUBROUTINE ZLAED0(Qsiz,N,D,E,Q,Ldq,Qstore,Ldqs,Rwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.D+0
      INTEGER :: Qsiz
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(Ldqs,*) :: Qstore
      INTEGER :: Ldqs
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLAED0
   END INTERFACE
END MODULE S_ZLAED0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAED7
   INTERFACE
      SUBROUTINE ZLAED7(N,Cutpnt,Qsiz,Tlvls,Curlvl,Curpbm,D,Q,Ldq,Rho,  &
     &                  Indxq,Qstore,Qptr,Prmptr,Perm,Givptr,Givcol,    &
     &                  Givnum,Work,Rwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Cutpnt
      INTEGER :: Qsiz
      INTEGER :: Tlvls
      INTEGER :: Curlvl
      INTEGER :: Curpbm
      REAL(R8KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) :: Rho
      INTEGER , DIMENSION(*) :: Indxq
      REAL(R8KIND) , DIMENSION(*) :: Qstore
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Qptr
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Prmptr
      INTEGER , DIMENSION(*) :: Perm
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Givptr
      INTEGER , DIMENSION(2,*) :: Givcol
      REAL(R8KIND) , DIMENSION(2,*) :: Givnum
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLAED7
   END INTERFACE
END MODULE S_ZLAED7
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAED8
   INTERFACE
      SUBROUTINE ZLAED8(K,N,Qsiz,Q,Ldq,D,Rho,Cutpnt,Z,Dlamda,Q2,Ldq2,W, &
     &                  Indxp,Indx,Indxq,Perm,Givptr,Givcol,Givnum,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  MONE = -1.0D0 , ZERO = 0.0D0 ,      &
     &                              ONE = 1.0D0 , TWO = 2.0D0 ,         &
     &                              EIGHT = 8.0D0
      INTEGER , INTENT(INOUT) :: K
      INTEGER :: N
      INTEGER :: Qsiz
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) :: Rho
      INTEGER , INTENT(IN) :: Cutpnt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dlamda
      COMPLEX(CX16KIND) , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxp
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxq
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Perm
      INTEGER , INTENT(INOUT) :: Givptr
      INTEGER , INTENT(OUT) , DIMENSION(2,*) :: Givcol
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(2,*) :: Givnum
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLAED8
   END INTERFACE
END MODULE S_ZLAED8
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAEIN
   INTERFACE
      SUBROUTINE ZLAEIN(Rightv,Noinit,N,H,Ldh,W,V,B,Ldb,Rwork,Eps3,     &
     &                  Smlnum,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , TENTH = 1.0D-1
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      LOGICAL , INTENT(IN) :: Rightv
      LOGICAL , INTENT(IN) :: Noinit
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldh,*) :: H
      INTEGER , INTENT(IN) :: Ldh
      COMPLEX(CX16KIND) , INTENT(IN) :: W
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: V
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      REAL(R8KIND) , INTENT(IN) :: Eps3
      REAL(R8KIND) , INTENT(IN) :: Smlnum
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE ZLAEIN
   END INTERFACE
END MODULE S_ZLAEIN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAESY
   INTERFACE
      SUBROUTINE ZLAESY(A,B,C,Rt1,Rt2,Evscal,Cs1,Sn1)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  HALF = 0.5D0 , THRESH = 0.1D0
      COMPLEX(CX16KIND) , INTENT(IN) :: A
      COMPLEX(CX16KIND) , INTENT(IN) :: B
      COMPLEX(CX16KIND) , INTENT(IN) :: C
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Rt1
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Rt2
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Evscal
      COMPLEX(CX16KIND) , INTENT(OUT) :: Cs1
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Sn1
      END SUBROUTINE ZLAESY
   END INTERFACE
END MODULE S_ZLAESY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAEV2
   INTERFACE
      SUBROUTINE ZLAEV2(A,B,C,Rt1,Rt2,Cs1,Sn1)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , INTENT(IN) :: A
      COMPLEX(CX16KIND) , INTENT(IN) :: B
      COMPLEX(CX16KIND) , INTENT(IN) :: C
      REAL(R8KIND) :: Rt1
      REAL(R8KIND) :: Rt2
      REAL(R8KIND) :: Cs1
      COMPLEX(CX16KIND) , INTENT(OUT) :: Sn1
      END SUBROUTINE ZLAEV2
   END INTERFACE
END MODULE S_ZLAEV2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAG2C
   INTERFACE
      SUBROUTINE ZLAG2C(M,N,A,Lda,Sa,Ldsa,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(OUT) , DIMENSION(Ldsa,*) :: Sa
      INTEGER , INTENT(IN) :: Ldsa
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE ZLAG2C
   END INTERFACE
END MODULE S_ZLAG2C
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_GBAMV
   INTERFACE
      SUBROUTINE ZLA_GBAMV(Trans,M,N,Kl,Ku,Alpha,Ab,Ldab,X,Incx,Beta,Y, &
     &                     Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL(R8KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE ZLA_GBAMV
   END INTERFACE
END MODULE S_ZLA_GBAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_GBRCOND_C
   INTERFACE
      FUNCTION ZLA_GBRCOND_C(Trans,N,Kl,Ku,Ab,Ldab,Afb,Ldafb,Ipiv,C,    &
     &                       Capply,Info,Work,Rwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_GBRCOND_C
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      LOGICAL , INTENT(IN) :: Capply
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION ZLA_GBRCOND_C
   END INTERFACE
END MODULE S_ZLA_GBRCOND_C
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_GBRCOND_X
   INTERFACE
      FUNCTION ZLA_GBRCOND_X(Trans,N,Kl,Ku,Ab,Ldab,Afb,Ldafb,Ipiv,X,    &
     &                       Info,Work,Rwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_GBRCOND_X
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION ZLA_GBRCOND_X
   END INTERFACE
END MODULE S_ZLA_GBRCOND_X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_GBRFSX_EXTENDED
   INTERFACE
      SUBROUTINE ZLA_GBRFSX_EXTENDED(Prec_type,Trans_type,N,Kl,Ku,Nrhs, &
     &                               Ab,Ldab,Afb,Ldafb,Ipiv,Colequ,C,B, &
     &                               Ldb,Y,Ldy,Berr_out,N_norms,        &
     &                               Err_bnds_norm,Err_bnds_comp,Res,   &
     &                               Ayb,Dy,Y_tail,Rcond,Ithresh,       &
     &                               Rthresh,Dz_ub,Ignore_cwise,Info)
      USE F77KINDS                        
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
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL(R8KIND) , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      COMPLEX(CX16KIND) , DIMENSION(*) :: Res
      REAL(R8KIND) , DIMENSION(*) :: Ayb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dy
      COMPLEX(CX16KIND) , DIMENSION(*) :: Y_tail
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL(R8KIND) , INTENT(IN) :: Rthresh
      REAL(R8KIND) , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
      END SUBROUTINE ZLA_GBRFSX_EXTENDED
   END INTERFACE
END MODULE S_ZLA_GBRFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_GBRPVGRW
   INTERFACE
      FUNCTION ZLA_GBRPVGRW(N,Kl,Ku,Ncols,Ab,Ldab,Afb,Ldafb)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_GBRPVGRW
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      INTEGER , INTENT(IN) :: Ncols
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldafb,*) :: Afb
      INTEGER , INTENT(IN) :: Ldafb
      END FUNCTION ZLA_GBRPVGRW
   END INTERFACE
END MODULE S_ZLA_GBRPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_GEAMV
   INTERFACE
      SUBROUTINE ZLA_GEAMV(Trans,M,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE ZLA_GEAMV
   END INTERFACE
END MODULE S_ZLA_GEAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_GERCOND_C
   INTERFACE
      FUNCTION ZLA_GERCOND_C(Trans,N,A,Lda,Af,Ldaf,Ipiv,C,Capply,Info,  &
     &                       Work,Rwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_GERCOND_C
      CHARACTER :: Trans
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      LOGICAL , INTENT(IN) :: Capply
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION ZLA_GERCOND_C
   END INTERFACE
END MODULE S_ZLA_GERCOND_C
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_GERCOND_X
   INTERFACE
      FUNCTION ZLA_GERCOND_X(Trans,N,A,Lda,Af,Ldaf,Ipiv,X,Info,Work,    &
     &                       Rwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_GERCOND_X
      CHARACTER :: Trans
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION ZLA_GERCOND_X
   END INTERFACE
END MODULE S_ZLA_GERCOND_X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_GERFSX_EXTENDED
   INTERFACE
      SUBROUTINE ZLA_GERFSX_EXTENDED(Prec_type,Trans_type,N,Nrhs,A,Lda, &
     &                               Af,Ldaf,Ipiv,Colequ,C,B,Ldb,Y,Ldy, &
     &                               Berr_out,N_norms,Errs_n,Errs_c,Res,&
     &                               Ayb,Dy,Y_tail,Rcond,Ithresh,       &
     &                               Rthresh,Dz_ub,Ignore_cwise,Info)
      USE F77KINDS                        
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
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL(R8KIND) , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Errs_n
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Errs_c
      COMPLEX(CX16KIND) , DIMENSION(*) :: Res
      REAL(R8KIND) , DIMENSION(*) :: Ayb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dy
      COMPLEX(CX16KIND) , DIMENSION(*) :: Y_tail
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL(R8KIND) , INTENT(IN) :: Rthresh
      REAL(R8KIND) , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
      END SUBROUTINE ZLA_GERFSX_EXTENDED
   END INTERFACE
END MODULE S_ZLA_GERFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_GERPVGRW
   INTERFACE
      FUNCTION ZLA_GERPVGRW(N,Ncols,A,Lda,Af,Ldaf)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_GERPVGRW
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ncols
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      END FUNCTION ZLA_GERPVGRW
   END INTERFACE
END MODULE S_ZLA_GERPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAGS2
   INTERFACE
      SUBROUTINE ZLAGS2(Upper,A1,A2,A3,B1,B2,B3,Csu,Snu,Csv,Snv,Csq,Snq)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      LOGICAL , INTENT(IN) :: Upper
      REAL(R8KIND) , INTENT(IN) :: A1
      COMPLEX(CX16KIND) , INTENT(IN) :: A2
      REAL(R8KIND) , INTENT(IN) :: A3
      REAL(R8KIND) , INTENT(IN) :: B1
      COMPLEX(CX16KIND) , INTENT(IN) :: B2
      REAL(R8KIND) , INTENT(IN) :: B3
      REAL(R8KIND) , INTENT(OUT) :: Csu
      COMPLEX(CX16KIND) , INTENT(OUT) :: Snu
      REAL(R8KIND) , INTENT(OUT) :: Csv
      COMPLEX(CX16KIND) , INTENT(OUT) :: Snv
      REAL(R8KIND) :: Csq
      COMPLEX(CX16KIND) :: Snq
      END SUBROUTINE ZLAGS2
   END INTERFACE
END MODULE S_ZLAGS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAGTM
   INTERFACE
      SUBROUTINE ZLAGTM(Trans,N,Nrhs,Alpha,Dl,D,Du,X,Ldx,Beta,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Dl
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Du
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE ZLAGTM
   END INTERFACE
END MODULE S_ZLAGTM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_HEAMV
   INTERFACE
      SUBROUTINE ZLA_HEAMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE ZLA_HEAMV
   END INTERFACE
END MODULE S_ZLA_HEAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAHEF_AA
   INTERFACE
      SUBROUTINE ZLAHEF_AA(Uplo,J1,M,Nb,A,Lda,Ipiv,H,Ldh,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: J1
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END SUBROUTINE ZLAHEF_AA
   END INTERFACE
END MODULE S_ZLAHEF_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAHEF
   INTERFACE
      SUBROUTINE ZLAHEF(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLAHEF
   END INTERFACE
END MODULE S_ZLAHEF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAHEF_RK
   INTERFACE
      SUBROUTINE ZLAHEF_RK(Uplo,N,Nb,Kb,A,Lda,E,Ipiv,W,Ldw,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLAHEF_RK
   END INTERFACE
END MODULE S_ZLAHEF_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAHEF_ROOK
   INTERFACE
      SUBROUTINE ZLAHEF_ROOK(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLAHEF_ROOK
   END INTERFACE
END MODULE S_ZLAHEF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_HERCOND_C
   INTERFACE
      FUNCTION ZLA_HERCOND_C(Uplo,N,A,Lda,Af,Ldaf,Ipiv,C,Capply,Info,   &
     &                       Work,Rwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_HERCOND_C
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      LOGICAL , INTENT(IN) :: Capply
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION ZLA_HERCOND_C
   END INTERFACE
END MODULE S_ZLA_HERCOND_C
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_HERCOND_X
   INTERFACE
      FUNCTION ZLA_HERCOND_X(Uplo,N,A,Lda,Af,Ldaf,Ipiv,X,Info,Work,     &
     &                       Rwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_HERCOND_X
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION ZLA_HERCOND_X
   END INTERFACE
END MODULE S_ZLA_HERCOND_X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_HERFSX_EXTENDED
   INTERFACE
      SUBROUTINE ZLA_HERFSX_EXTENDED(Prec_type,Uplo,N,Nrhs,A,Lda,Af,    &
     &                               Ldaf,Ipiv,Colequ,C,B,Ldb,Y,Ldy,    &
     &                               Berr_out,N_norms,Err_bnds_norm,    &
     &                               Err_bnds_comp,Res,Ayb,Dy,Y_tail,   &
     &                               Rcond,Ithresh,Rthresh,Dz_ub,       &
     &                               Ignore_cwise,Info)
      USE F77KINDS                        
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
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL(R8KIND) , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      COMPLEX(CX16KIND) , DIMENSION(*) :: Res
      REAL(R8KIND) , DIMENSION(*) :: Ayb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dy
      COMPLEX(CX16KIND) , DIMENSION(*) :: Y_tail
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL(R8KIND) , INTENT(IN) :: Rthresh
      REAL(R8KIND) , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLA_HERFSX_EXTENDED
   END INTERFACE
END MODULE S_ZLA_HERFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_HERPVGRW
   INTERFACE
      FUNCTION ZLA_HERPVGRW(Uplo,N,Info,A,Lda,Af,Ldaf,Ipiv,Work)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_HERPVGRW
      CHARACTER(1) :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Info
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLA_HERPVGRW
   END INTERFACE
END MODULE S_ZLA_HERPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAHQR
   INTERFACE
      SUBROUTINE ZLAHQR(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,W,Iloz,Ihiz,Z,Ldz,  &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D0,0.0D0) ,         &
     &                 ONE = (1.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  RZERO = 0.0D0 , RONE = 1.0D0 ,      &
     &                              HALF = 0.5D0 , DAT1 = 3.0D0/4.0D0
      INTEGER , PARAMETER  ::  KEXSH = 10
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: W
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE ZLAHQR
   END INTERFACE
END MODULE S_ZLAHQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAHR2
   INTERFACE
      SUBROUTINE ZLAHR2(N,K,Nb,A,Lda,Tau,T,Ldt,Y,Ldy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: N
      INTEGER :: K
      INTEGER :: Nb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Nb) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldt,Nb) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Ldy,Nb) :: Y
      INTEGER :: Ldy
      END SUBROUTINE ZLAHR2
   END INTERFACE
END MODULE S_ZLAHR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAIC1
   INTERFACE
      SUBROUTINE ZLAIC1(Job,J,X,Sest,W,Gamma,Sestpr,S,C)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , HALF = 0.5D0 ,        &
     &                              FOUR = 4.0D0
      INTEGER , INTENT(IN) :: Job
      INTEGER :: J
      COMPLEX(CX16KIND) , DIMENSION(J) :: X
      REAL(R8KIND) , INTENT(IN) :: Sest
      COMPLEX(CX16KIND) , DIMENSION(J) :: W
      COMPLEX(CX16KIND) , INTENT(IN) :: Gamma
      REAL(R8KIND) , INTENT(OUT) :: Sestpr
      COMPLEX(CX16KIND) , INTENT(INOUT) :: S
      COMPLEX(CX16KIND) , INTENT(INOUT) :: C
      END SUBROUTINE ZLAIC1
   END INTERFACE
END MODULE S_ZLAIC1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_LIN_BERR
   INTERFACE
      SUBROUTINE ZLA_LIN_BERR(N,Nz,Nrhs,Res,Ayb,Berr)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nz
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(N,Nrhs) :: Res
      REAL(R8KIND) , INTENT(IN) , DIMENSION(N,Nrhs) :: Ayb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs) :: Berr
      END SUBROUTINE ZLA_LIN_BERR
   END INTERFACE
END MODULE S_ZLA_LIN_BERR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLALS0
   INTERFACE
      SUBROUTINE ZLALS0(Icompq,Nl,Nr,Sqre,Nrhs,B,Ldb,Bx,Ldbx,Perm,      &
     &                  Givptr,Givcol,Ldgcol,Givnum,Ldgnum,Poles,Difl,  &
     &                  Difr,Z,K,C,S,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0 ,        &
     &                              NEGONE = -1.0D0
      INTEGER , INTENT(IN) :: Icompq
      INTEGER , INTENT(IN) :: Nl
      INTEGER , INTENT(IN) :: Nr
      INTEGER , INTENT(IN) :: Sqre
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldbx,*) :: Bx
      INTEGER :: Ldbx
      INTEGER , INTENT(IN) , DIMENSION(*) :: Perm
      INTEGER , INTENT(IN) :: Givptr
      INTEGER , INTENT(IN) , DIMENSION(Ldgcol,*) :: Givcol
      INTEGER , INTENT(IN) :: Ldgcol
      REAL(R8KIND) , DIMENSION(Ldgnum,*) :: Givnum
      INTEGER , INTENT(IN) :: Ldgnum
      REAL(R8KIND) , DIMENSION(Ldgnum,*) :: Poles
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Difl
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldgnum,*) :: Difr
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Z
      INTEGER :: K
      REAL(R8KIND) :: C
      REAL(R8KIND) :: S
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLALS0
   END INTERFACE
END MODULE S_ZLALS0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLALSA
   INTERFACE
      SUBROUTINE ZLALSA(Icompq,Smlsiz,N,Nrhs,B,Ldb,Bx,Ldbx,U,Ldu,Vt,K,  &
     &                  Difl,Difr,Z,Poles,Givptr,Givcol,Ldgcol,Perm,    &
     &                  Givnum,C,S,Rwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      INTEGER :: Icompq
      INTEGER :: Smlsiz
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldbx,*) :: Bx
      INTEGER :: Ldbx
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldu,*) :: Vt
      INTEGER , DIMENSION(*) :: K
      REAL(R8KIND) , DIMENSION(Ldu,*) :: Difl
      REAL(R8KIND) , DIMENSION(Ldu,*) :: Difr
      REAL(R8KIND) , DIMENSION(Ldu,*) :: Z
      REAL(R8KIND) , DIMENSION(Ldu,*) :: Poles
      INTEGER , DIMENSION(*) :: Givptr
      INTEGER , DIMENSION(Ldgcol,*) :: Givcol
      INTEGER :: Ldgcol
      INTEGER , DIMENSION(Ldgcol,*) :: Perm
      REAL(R8KIND) , DIMENSION(Ldu,*) :: Givnum
      REAL(R8KIND) , DIMENSION(*) :: C
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLALSA
   END INTERFACE
END MODULE S_ZLALSA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLALSD
   INTERFACE
      SUBROUTINE ZLALSD(Uplo,Smlsiz,N,Nrhs,D,E,B,Ldb,Rcond,Rank,Work,   &
     &                  Rwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0)
      CHARACTER , INTENT(IN) :: Uplo
      INTEGER :: Smlsiz
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLALSD
   END INTERFACE
END MODULE S_ZLALSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAMSWLQ
   INTERFACE
      SUBROUTINE ZLAMSWLQ(Side,Trans,M,N,K,Mb,Nb,A,Lda,T,Ldt,C,Ldc,Work,&
     &                    Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: Mb
      INTEGER :: Nb
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLAMSWLQ
   END INTERFACE
END MODULE S_ZLAMSWLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAMTSQR
   INTERFACE
      SUBROUTINE ZLAMTSQR(Side,Trans,M,N,K,Mb,Nb,A,Lda,T,Ldt,C,Ldc,Work,&
     &                    Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: Mb
      INTEGER :: Nb
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLAMTSQR
   END INTERFACE
END MODULE S_ZLAMTSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANGB
   INTERFACE
      FUNCTION ZLANGB(Norm,N,Kl,Ku,Ab,Ldab,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANGB
      CHARACTER :: Norm
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLANGB
   END INTERFACE
END MODULE S_ZLANGB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANGE
   INTERFACE
      FUNCTION ZLANGE(Norm,M,N,A,Lda,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANGE
      CHARACTER :: Norm
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLANGE
   END INTERFACE
END MODULE S_ZLANGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANGT
   INTERFACE
      FUNCTION ZLANGT(Norm,N,Dl,D,Du)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANGT
      CHARACTER :: Norm
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dl
      COMPLEX(CX16KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: Du
      END FUNCTION ZLANGT
   END INTERFACE
END MODULE S_ZLANGT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANHB
   INTERFACE
      FUNCTION ZLANHB(Norm,Uplo,N,K,Ab,Ldab,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANHB
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLANHB
   END INTERFACE
END MODULE S_ZLANHB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANHE
   INTERFACE
      FUNCTION ZLANHE(Norm,Uplo,N,A,Lda,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANHE
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLANHE
   END INTERFACE
END MODULE S_ZLANHE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANHF
   INTERFACE
      FUNCTION ZLANHF(Norm,Transr,Uplo,N,A,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANHF
      CHARACTER :: Norm
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(0:*) :: A
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(0:*) :: Work
      END FUNCTION ZLANHF
   END INTERFACE
END MODULE S_ZLANHF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANHP
   INTERFACE
      FUNCTION ZLANHP(Norm,Uplo,N,Ap,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANHP
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLANHP
   END INTERFACE
END MODULE S_ZLANHP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANHS
   INTERFACE
      FUNCTION ZLANHS(Norm,N,A,Lda,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANHS
      CHARACTER :: Norm
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLANHS
   END INTERFACE
END MODULE S_ZLANHS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANHT
   INTERFACE
      FUNCTION ZLANHT(Norm,N,D,E)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANHT
      CHARACTER :: Norm
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      END FUNCTION ZLANHT
   END INTERFACE
END MODULE S_ZLANHT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANSB
   INTERFACE
      FUNCTION ZLANSB(Norm,Uplo,N,K,Ab,Ldab,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANSB
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLANSB
   END INTERFACE
END MODULE S_ZLANSB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANSP
   INTERFACE
      FUNCTION ZLANSP(Norm,Uplo,N,Ap,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANSP
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLANSP
   END INTERFACE
END MODULE S_ZLANSP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANSY
   INTERFACE
      FUNCTION ZLANSY(Norm,Uplo,N,A,Lda,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANSY
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLANSY
   END INTERFACE
END MODULE S_ZLANSY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANTB
   INTERFACE
      FUNCTION ZLANTB(Norm,Uplo,Diag,N,K,Ab,Ldab,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANTB
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLANTB
   END INTERFACE
END MODULE S_ZLANTB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANTP
   INTERFACE
      FUNCTION ZLANTP(Norm,Uplo,Diag,N,Ap,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANTP
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLANTP
   END INTERFACE
END MODULE S_ZLANTP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLANTR
   INTERFACE
      FUNCTION ZLANTR(Norm,Uplo,Diag,M,N,A,Lda,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANTR
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLANTR
   END INTERFACE
END MODULE S_ZLANTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAPLL
   INTERFACE
      SUBROUTINE ZLAPLL(N,X,Incx,Y,Incy,Ssmin)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER :: Incx
      COMPLEX(CX16KIND) , DIMENSION(*) :: Y
      INTEGER :: Incy
      REAL(R8KIND) :: Ssmin
      END SUBROUTINE ZLAPLL
   END INTERFACE
END MODULE S_ZLAPLL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAPMR
   INTERFACE
      SUBROUTINE ZLAPMR(Forwrd,M,N,X,Ldx,K)
      USE F77KINDS                        
      IMPLICIT NONE
      LOGICAL , INTENT(IN) :: Forwrd
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: K
      END SUBROUTINE ZLAPMR
   END INTERFACE
END MODULE S_ZLAPMR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAPMT
   INTERFACE
      SUBROUTINE ZLAPMT(Forwrd,M,N,X,Ldx,K)
      USE F77KINDS                        
      IMPLICIT NONE
      LOGICAL , INTENT(IN) :: Forwrd
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: K
      END SUBROUTINE ZLAPMT
   END INTERFACE
END MODULE S_ZLAPMT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_PORCOND_C
   INTERFACE
      FUNCTION ZLA_PORCOND_C(Uplo,N,A,Lda,Af,Ldaf,C,Capply,Info,Work,   &
     &                       Rwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_PORCOND_C
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      LOGICAL , INTENT(IN) :: Capply
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION ZLA_PORCOND_C
   END INTERFACE
END MODULE S_ZLA_PORCOND_C
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_PORCOND_X
   INTERFACE
      FUNCTION ZLA_PORCOND_X(Uplo,N,A,Lda,Af,Ldaf,X,Info,Work,Rwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_PORCOND_X
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION ZLA_PORCOND_X
   END INTERFACE
END MODULE S_ZLA_PORCOND_X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_PORFSX_EXTENDED
   INTERFACE
      SUBROUTINE ZLA_PORFSX_EXTENDED(Prec_type,Uplo,N,Nrhs,A,Lda,Af,    &
     &                               Ldaf,Colequ,C,B,Ldb,Y,Ldy,Berr_out,&
     &                               N_norms,Err_bnds_norm,             &
     &                               Err_bnds_comp,Res,Ayb,Dy,Y_tail,   &
     &                               Rcond,Ithresh,Rthresh,Dz_ub,       &
     &                               Ignore_cwise,Info)
      USE F77KINDS                        
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
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      LOGICAL , INTENT(IN) :: Colequ
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL(R8KIND) , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      COMPLEX(CX16KIND) , DIMENSION(*) :: Res
      REAL(R8KIND) , DIMENSION(*) :: Ayb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dy
      COMPLEX(CX16KIND) , DIMENSION(*) :: Y_tail
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL(R8KIND) , INTENT(IN) :: Rthresh
      REAL(R8KIND) , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
      END SUBROUTINE ZLA_PORFSX_EXTENDED
   END INTERFACE
END MODULE S_ZLA_PORFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_PORPVGRW
   INTERFACE
      FUNCTION ZLA_PORPVGRW(Uplo,Ncols,A,Lda,Af,Ldaf,Work)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_PORPVGRW
      CHARACTER(1) :: Uplo
      INTEGER , INTENT(IN) :: Ncols
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLA_PORPVGRW
   END INTERFACE
END MODULE S_ZLA_PORPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQGB
   INTERFACE
      SUBROUTINE ZLAQGB(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Equed)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(IN) :: Rowcnd
      REAL(R8KIND) , INTENT(IN) :: Colcnd
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE ZLAQGB
   END INTERFACE
END MODULE S_ZLAQGB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQGE
   INTERFACE
      SUBROUTINE ZLAQGE(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Equed)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(IN) :: Rowcnd
      REAL(R8KIND) , INTENT(IN) :: Colcnd
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE ZLAQGE
   END INTERFACE
END MODULE S_ZLAQGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQHB
   INTERFACE
      SUBROUTINE ZLAQHB(Uplo,N,Kd,Ab,Ldab,S,Scond,Amax,Equed)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Scond
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE ZLAQHB
   END INTERFACE
END MODULE S_ZLAQHB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQHE
   INTERFACE
      SUBROUTINE ZLAQHE(Uplo,N,A,Lda,S,Scond,Amax,Equed)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Scond
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE ZLAQHE
   END INTERFACE
END MODULE S_ZLAQHE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQHP
   INTERFACE
      SUBROUTINE ZLAQHP(Uplo,N,Ap,S,Scond,Amax,Equed)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Scond
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE ZLAQHP
   END INTERFACE
END MODULE S_ZLAQHP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQP2
   INTERFACE
      SUBROUTINE ZLAQP2(M,N,Offset,A,Lda,Jpvt,Tau,Vn1,Vn2,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Offset
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vn1
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vn2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      END SUBROUTINE ZLAQP2
   END INTERFACE
END MODULE S_ZLAQP2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQPS
   INTERFACE
      SUBROUTINE ZLAQPS(M,N,Offset,Nb,Kb,A,Lda,Jpvt,Tau,Vn1,Vn2,Auxv,F, &
     &                  Ldf)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Offset
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Kb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vn1
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vn2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Auxv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      END SUBROUTINE ZLAQPS
   END INTERFACE
END MODULE S_ZLAQPS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQR0
   INTERFACE
      SUBROUTINE ZLAQR0(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,W,Iloz,Ihiz,Z,Ldz,  &
     &                  Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NTINY = 15 , KEXNW = 5 , KEXSH = 6
      REAL(R8KIND) , PARAMETER  ::  WILK1 = 0.75D0
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D0,0.0D0) ,         &
     &                 ONE = (1.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D0
      LOGICAL :: Wantt
      LOGICAL :: Wantz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX(CX16KIND) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER :: Iloz
      INTEGER :: Ihiz
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER :: Info
      END SUBROUTINE ZLAQR0
   END INTERFACE
END MODULE S_ZLAQR0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQR1
   INTERFACE
      SUBROUTINE ZLAQR1(N,H,Ldh,S1,S2,V)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  RZERO = 0.0D0
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldh,*) :: H
      INTEGER , INTENT(IN) :: Ldh
      COMPLEX(CX16KIND) , INTENT(IN) :: S1
      COMPLEX(CX16KIND) , INTENT(IN) :: S2
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: V
      END SUBROUTINE ZLAQR1
   END INTERFACE
END MODULE S_ZLAQR1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQR2
   INTERFACE
      SUBROUTINE ZLAQR2(Wantt,Wantz,N,Ktop,Kbot,Nw,H,Ldh,Iloz,Ihiz,Z,   &
     &                  Ldz,Ns,Nd,Sh,V,Ldv,Nh,T,Ldt,Nv,Wv,Ldwv,Work,    &
     &                  Lwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D0,0.0D0) ,         &
     &                 ONE = (1.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  RZERO = 0.0D0 , RONE = 1.0D0
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nw
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ns
      INTEGER , INTENT(OUT) :: Nd
      COMPLEX(CX16KIND) , DIMENSION(*) :: Sh
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(IN) :: Nh
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(IN) :: Nv
      COMPLEX(CX16KIND) , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      END SUBROUTINE ZLAQR2
   END INTERFACE
END MODULE S_ZLAQR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQR3
   INTERFACE
      SUBROUTINE ZLAQR3(Wantt,Wantz,N,Ktop,Kbot,Nw,H,Ldh,Iloz,Ihiz,Z,   &
     &                  Ldz,Ns,Nd,Sh,V,Ldv,Nh,T,Ldt,Nv,Wv,Ldwv,Work,    &
     &                  Lwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D0,0.0D0) ,         &
     &                 ONE = (1.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  RZERO = 0.0D0 , RONE = 1.0D0
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nw
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ns
      INTEGER , INTENT(OUT) :: Nd
      COMPLEX(CX16KIND) , DIMENSION(*) :: Sh
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(IN) :: Nh
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(IN) :: Nv
      COMPLEX(CX16KIND) , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      END SUBROUTINE ZLAQR3
   END INTERFACE
END MODULE S_ZLAQR3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQR4
   INTERFACE
      SUBROUTINE ZLAQR4(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,W,Iloz,Ihiz,Z,Ldz,  &
     &                  Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NTINY = 15 , KEXNW = 5 , KEXSH = 6
      REAL(R8KIND) , PARAMETER  ::  WILK1 = 0.75D0
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D0,0.0D0) ,         &
     &                 ONE = (1.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D0
      LOGICAL :: Wantt
      LOGICAL :: Wantz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX(CX16KIND) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER :: Iloz
      INTEGER :: Ihiz
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER :: Info
      END SUBROUTINE ZLAQR4
   END INTERFACE
END MODULE S_ZLAQR4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQR5
   INTERFACE
      SUBROUTINE ZLAQR5(Wantt,Wantz,Kacc22,N,Ktop,Kbot,Nshfts,S,H,Ldh,  &
     &                  Iloz,Ihiz,Z,Ldz,V,Ldv,U,Ldu,Nv,Wv,Ldwv,Nh,Wh,   &
     &                  Ldwh)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D0,0.0D0) ,         &
     &                 ONE = (1.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  RZERO = 0.0D0 , RONE = 1.0D0
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: Kacc22
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nshfts
      COMPLEX(CX16KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      INTEGER , INTENT(IN) :: Nv
      COMPLEX(CX16KIND) , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      INTEGER , INTENT(IN) :: Nh
      COMPLEX(CX16KIND) , DIMENSION(Ldwh,*) :: Wh
      INTEGER :: Ldwh
      END SUBROUTINE ZLAQR5
   END INTERFACE
END MODULE S_ZLAQR5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQSB
   INTERFACE
      SUBROUTINE ZLAQSB(Uplo,N,Kd,Ab,Ldab,S,Scond,Amax,Equed)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Scond
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE ZLAQSB
   END INTERFACE
END MODULE S_ZLAQSB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQSP
   INTERFACE
      SUBROUTINE ZLAQSP(Uplo,N,Ap,S,Scond,Amax,Equed)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Scond
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE ZLAQSP
   END INTERFACE
END MODULE S_ZLAQSP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAQSY
   INTERFACE
      SUBROUTINE ZLAQSY(Uplo,N,A,Lda,S,Scond,Amax,Equed)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Scond
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE ZLAQSY
   END INTERFACE
END MODULE S_ZLAQSY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAR1V
   INTERFACE
      SUBROUTINE ZLAR1V(N,B1,Bn,Lambda,D,L,Ld,Lld,Pivmin,Gaptol,Z,      &
     &                  Wantnc,Negcnt,Ztz,Mingma,R,Isuppz,Nrminv,Resid, &
     &                  Rqcorr,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D0,0.0D0)
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: B1
      INTEGER , INTENT(IN) :: Bn
      REAL(R8KIND) , INTENT(IN) :: Lambda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: L
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Ld
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Lld
      REAL(R8KIND) , INTENT(IN) :: Pivmin
      REAL(R8KIND) , INTENT(IN) :: Gaptol
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      LOGICAL , INTENT(IN) :: Wantnc
      INTEGER , INTENT(OUT) :: Negcnt
      REAL(R8KIND) , INTENT(INOUT) :: Ztz
      REAL(R8KIND) , INTENT(INOUT) :: Mingma
      INTEGER , INTENT(INOUT) :: R
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Isuppz
      REAL(R8KIND) , INTENT(INOUT) :: Nrminv
      REAL(R8KIND) , INTENT(OUT) :: Resid
      REAL(R8KIND) , INTENT(OUT) :: Rqcorr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END SUBROUTINE ZLAR1V
   END INTERFACE
END MODULE S_ZLAR1V
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAR2V
   INTERFACE
      SUBROUTINE ZLAR2V(N,X,Y,Z,Incx,C,S,Incc)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: S
      INTEGER , INTENT(IN) :: Incc
      END SUBROUTINE ZLAR2V
   END INTERFACE
END MODULE S_ZLAR2V
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARCM
   INTERFACE
      SUBROUTINE ZLARCM(M,N,A,Lda,B,Ldb,C,Ldc,Rwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END SUBROUTINE ZLARCM
   END INTERFACE
END MODULE S_ZLARCM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARFB
   INTERFACE
      SUBROUTINE ZLARFB(Side,Trans,Direct,Storev,M,N,K,V,Ldv,T,Ldt,C,   &
     &                  Ldc,Work,Ldwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE ZLARFB
   END INTERFACE
END MODULE S_ZLARFB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARFB_GETT
   INTERFACE
      SUBROUTINE ZLARFB_GETT(Ident,M,N,K,T,Ldt,A,Lda,B,Ldb,Work,Ldwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Ident
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: K
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE ZLARFB_GETT
   END INTERFACE
END MODULE S_ZLARFB_GETT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARF
   INTERFACE
      SUBROUTINE ZLARF(Side,M,N,V,Incv,Tau,C,Ldc,Work)
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
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      END SUBROUTINE ZLARF
   END INTERFACE
END MODULE S_ZLARF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARFG
   INTERFACE
      SUBROUTINE ZLARFG(N,Alpha,X,Incx,Tau)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Alpha
      COMPLEX(CX16KIND) , DIMENSION(*) :: X
      INTEGER :: Incx
      COMPLEX(CX16KIND) , INTENT(OUT) :: Tau
      END SUBROUTINE ZLARFG
   END INTERFACE
END MODULE S_ZLARFG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARFGP
   INTERFACE
      SUBROUTINE ZLARFGP(N,Alpha,X,Incx,Tau)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , ONE = 1.0D+0 ,       &
     &                              ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Alpha
      COMPLEX(CX16KIND) , DIMENSION(*) :: X
      INTEGER :: Incx
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Tau
      END SUBROUTINE ZLARFGP
   END INTERFACE
END MODULE S_ZLARFGP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARFT
   INTERFACE
      SUBROUTINE ZLARFT(Direct,Storev,N,K,V,Ldv,Tau,T,Ldt)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      END SUBROUTINE ZLARFT
   END INTERFACE
END MODULE S_ZLARFT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARFX
   INTERFACE
      SUBROUTINE ZLARFX(Side,M,N,V,Tau,C,Ldc,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: V
      COMPLEX(CX16KIND) :: Tau
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      END SUBROUTINE ZLARFX
   END INTERFACE
END MODULE S_ZLARFX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARFY
   INTERFACE
      SUBROUTINE ZLARFY(Uplo,N,V,Incv,Tau,C,Ldc,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0) , HALF = (0.5D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: V
      INTEGER :: Incv
      COMPLEX(CX16KIND) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      END SUBROUTINE ZLARFY
   END INTERFACE
END MODULE S_ZLARFY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARGV
   INTERFACE
      SUBROUTINE ZLARGV(N,X,Incx,Y,Incy,C,Incc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , ONE = 1.0D+0 ,       &
     &                              ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: C
      INTEGER , INTENT(IN) :: Incc
      END SUBROUTINE ZLARGV
   END INTERFACE
END MODULE S_ZLARGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARNV
   INTERFACE
      SUBROUTINE ZLARNV(Idist,Iseed,N,X)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
      INTEGER , PARAMETER  ::  LV = 128
      REAL(R8KIND) , PARAMETER  ::  TWOPI =                             &
     &                       6.28318530717958647692528676655900576839D+0
      INTEGER , INTENT(IN) :: Idist
      INTEGER , DIMENSION(4) :: Iseed
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: X
      END SUBROUTINE ZLARNV
   END INTERFACE
END MODULE S_ZLARNV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARRV
   INTERFACE
      SUBROUTINE ZLARRV(N,Vl,Vu,D,L,Pivmin,Isplit,M,Dol,Dou,Minrgp,     &
     &                  Rtol1,Rtol2,W,Werr,Wgap,Iblock,Indexw,Gers,Z,   &
     &                  Ldz,Isuppz,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXITR = 10
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , THREE = 3.0D0 ,       &
     &                              FOUR = 4.0D0 , HALF = 0.5D0
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) :: Vu
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: L
      REAL(R8KIND) :: Pivmin
      INTEGER , INTENT(IN) , DIMENSION(*) :: Isplit
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: Dol
      INTEGER , INTENT(IN) :: Dou
      REAL(R8KIND) , INTENT(IN) :: Minrgp
      REAL(R8KIND) :: Rtol1
      REAL(R8KIND) :: Rtol2
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wgap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(IN) , DIMENSION(*) :: Indexw
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Gers
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isuppz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE ZLARRV
   END INTERFACE
END MODULE S_ZLARRV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARSCL2
   INTERFACE
      SUBROUTINE ZLARSCL2(M,N,D,X,Ldx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      END SUBROUTINE ZLARSCL2
   END INTERFACE
END MODULE S_ZLARSCL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARTG
   INTERFACE
      SUBROUTINE ZLARTG(F,G,Cs,Sn,R)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , ONE = 1.0D+0 ,       &
     &                              ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
      COMPLEX(CX16KIND) , INTENT(IN) :: F
      COMPLEX(CX16KIND) , INTENT(IN) :: G
      REAL(R8KIND) , INTENT(INOUT) :: Cs
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Sn
      COMPLEX(CX16KIND) , INTENT(INOUT) :: R
      END SUBROUTINE ZLARTG
   END INTERFACE
END MODULE S_ZLARTG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARTV
   INTERFACE
      SUBROUTINE ZLARTV(N,X,Incx,Y,Incy,C,S,Incc)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: S
      INTEGER , INTENT(IN) :: Incc
      END SUBROUTINE ZLARTV
   END INTERFACE
END MODULE S_ZLARTV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARZB
   INTERFACE
      SUBROUTINE ZLARZB(Side,Trans,Direct,Storev,M,N,K,L,V,Ldv,T,Ldt,C, &
     &                  Ldc,Work,Ldwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: L
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE ZLARZB
   END INTERFACE
END MODULE S_ZLARZB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARZ
   INTERFACE
      SUBROUTINE ZLARZ(Side,M,N,L,V,Incv,Tau,C,Ldc,Work)
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
      INTEGER :: L
      COMPLEX(CX16KIND) , DIMENSION(*) :: V
      INTEGER :: Incv
      COMPLEX(CX16KIND) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      END SUBROUTINE ZLARZ
   END INTERFACE
END MODULE S_ZLARZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLARZT
   INTERFACE
      SUBROUTINE ZLARZT(Direct,Storev,N,K,V,Ldv,Tau,T,Ldt)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      END SUBROUTINE ZLARZT
   END INTERFACE
END MODULE S_ZLARZT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLASCL2
   INTERFACE
      SUBROUTINE ZLASCL2(M,N,D,X,Ldx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      END SUBROUTINE ZLASCL2
   END INTERFACE
END MODULE S_ZLASCL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLASCL
   INTERFACE
      SUBROUTINE ZLASCL(Type,Kl,Ku,Cfrom,Cto,M,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Type
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL(R8KIND) :: Cfrom
      REAL(R8KIND) :: Cto
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLASCL
   END INTERFACE
END MODULE S_ZLASCL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLASET
   INTERFACE
      SUBROUTINE ZLASET(Uplo,M,N,Alpha,Beta,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE ZLASET
   END INTERFACE
END MODULE S_ZLASET
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLASR
   INTERFACE
      SUBROUTINE ZLASR(Side,Pivot,Direct,M,N,C,S,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Side
      CHARACTER :: Pivot
      CHARACTER :: Direct
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE ZLASR
   END INTERFACE
END MODULE S_ZLASR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLASSQ
   INTERFACE
      SUBROUTINE ZLASSQ(N,X,Incx,Scale,Sumsq)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) :: Sumsq
      END SUBROUTINE ZLASSQ
   END INTERFACE
END MODULE S_ZLASSQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLASWLQ
   INTERFACE
      SUBROUTINE ZLASWLQ(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER :: Nb
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLASWLQ
   END INTERFACE
END MODULE S_ZLASWLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLASWP
   INTERFACE
      SUBROUTINE ZLASWP(N,A,Lda,K1,K2,Ipiv,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) :: K1
      INTEGER , INTENT(IN) :: K2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE ZLASWP
   END INTERFACE
END MODULE S_ZLASWP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_SYAMV
   INTERFACE
      SUBROUTINE ZLA_SYAMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE ZLA_SYAMV
   END INTERFACE
END MODULE S_ZLA_SYAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLASYF_AA
   INTERFACE
      SUBROUTINE ZLASYF_AA(Uplo,J1,M,Nb,A,Lda,Ipiv,H,Ldh,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: J1
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END SUBROUTINE ZLASYF_AA
   END INTERFACE
END MODULE S_ZLASYF_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLASYF
   INTERFACE
      SUBROUTINE ZLASYF(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLASYF
   END INTERFACE
END MODULE S_ZLASYF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLASYF_RK
   INTERFACE
      SUBROUTINE ZLASYF_RK(Uplo,N,Nb,Kb,A,Lda,E,Ipiv,W,Ldw,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLASYF_RK
   END INTERFACE
END MODULE S_ZLASYF_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLASYF_ROOK
   INTERFACE
      SUBROUTINE ZLASYF_ROOK(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLASYF_ROOK
   END INTERFACE
END MODULE S_ZLASYF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_SYRCOND_C
   INTERFACE
      FUNCTION ZLA_SYRCOND_C(Uplo,N,A,Lda,Af,Ldaf,Ipiv,C,Capply,Info,   &
     &                       Work,Rwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_SYRCOND_C
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      LOGICAL , INTENT(IN) :: Capply
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION ZLA_SYRCOND_C
   END INTERFACE
END MODULE S_ZLA_SYRCOND_C
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_SYRCOND_X
   INTERFACE
      FUNCTION ZLA_SYRCOND_X(Uplo,N,A,Lda,Af,Ldaf,Ipiv,X,Info,Work,     &
     &                       Rwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_SYRCOND_X
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      END FUNCTION ZLA_SYRCOND_X
   END INTERFACE
END MODULE S_ZLA_SYRCOND_X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_SYRFSX_EXTENDED
   INTERFACE
      SUBROUTINE ZLA_SYRFSX_EXTENDED(Prec_type,Uplo,N,Nrhs,A,Lda,Af,    &
     &                               Ldaf,Ipiv,Colequ,C,B,Ldb,Y,Ldy,    &
     &                               Berr_out,N_norms,Err_bnds_norm,    &
     &                               Err_bnds_comp,Res,Ayb,Dy,Y_tail,   &
     &                               Rcond,Ithresh,Rthresh,Dz_ub,       &
     &                               Ignore_cwise,Info)
      USE F77KINDS                        
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
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL(R8KIND) , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      COMPLEX(CX16KIND) , DIMENSION(*) :: Res
      REAL(R8KIND) , DIMENSION(*) :: Ayb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dy
      COMPLEX(CX16KIND) , DIMENSION(*) :: Y_tail
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL(R8KIND) , INTENT(IN) :: Rthresh
      REAL(R8KIND) , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLA_SYRFSX_EXTENDED
   END INTERFACE
END MODULE S_ZLA_SYRFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_SYRPVGRW
   INTERFACE
      FUNCTION ZLA_SYRPVGRW(Uplo,N,Info,A,Lda,Af,Ldaf,Ipiv,Work)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: ZLA_SYRPVGRW
      CHARACTER(1) :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Info
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION ZLA_SYRPVGRW
   END INTERFACE
END MODULE S_ZLA_SYRPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAT2C
   INTERFACE
      SUBROUTINE ZLAT2C(Uplo,N,A,Lda,Sa,Ldsa,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(OUT) , DIMENSION(Ldsa,*) :: Sa
      INTEGER , INTENT(IN) :: Ldsa
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE ZLAT2C
   END INTERFACE
END MODULE S_ZLAT2C
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLATBS
   INTERFACE
      SUBROUTINE ZLATBS(Uplo,Trans,Diag,Normin,N,Kd,Ab,Ldab,X,Scale,    &
     &                  Cnorm,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , HALF = 0.5D+0 ,     &
     &                              ONE = 1.0D+0 , TWO = 2.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLATBS
   END INTERFACE
END MODULE S_ZLATBS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLATDF
   INTERFACE
      SUBROUTINE ZLATDF(Ijob,N,Z,Ldz,Rhs,Rdsum,Rdscal,Ipiv,Jpiv)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXDIM = 2
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: Ijob
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Rhs
      REAL(R8KIND) :: Rdsum
      REAL(R8KIND) :: Rdscal
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Jpiv
      END SUBROUTINE ZLATDF
   END INTERFACE
END MODULE S_ZLATDF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLATPS
   INTERFACE
      SUBROUTINE ZLATPS(Uplo,Trans,Diag,Normin,N,Ap,X,Scale,Cnorm,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , HALF = 0.5D+0 ,     &
     &                              ONE = 1.0D+0 , TWO = 2.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLATPS
   END INTERFACE
END MODULE S_ZLATPS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLATRD
   INTERFACE
      SUBROUTINE ZLATRD(Uplo,N,Nb,A,Lda,E,Tau,W,Ldw)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0) , HALF = (0.5D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      END SUBROUTINE ZLATRD
   END INTERFACE
END MODULE S_ZLATRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLATRS
   INTERFACE
      SUBROUTINE ZLATRS(Uplo,Trans,Diag,Normin,N,A,Lda,X,Scale,Cnorm,   &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , HALF = 0.5D+0 ,     &
     &                              ONE = 1.0D+0 , TWO = 2.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLATRS
   END INTERFACE
END MODULE S_ZLATRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLATRZ
   INTERFACE
      SUBROUTINE ZLATRZ(M,N,L,A,Lda,Tau,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: L
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      END SUBROUTINE ZLATRZ
   END INTERFACE
END MODULE S_ZLATRZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLATSQR
   INTERFACE
      SUBROUTINE ZLATSQR(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER :: Nb
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLATSQR
   END INTERFACE
END MODULE S_ZLATSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAUNHR_COL_GETRFNP2
   INTERFACE
      RECURSIVE SUBROUTINE ZLAUNHR_COL_GETRFNP2(M,N,A,Lda,D,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLAUNHR_COL_GETRFNP2
   END INTERFACE
END MODULE S_ZLAUNHR_COL_GETRFNP2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAUNHR_COL_GETRFNP
   INTERFACE
      SUBROUTINE ZLAUNHR_COL_GETRFNP(M,N,A,Lda,D,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLAUNHR_COL_GETRFNP
   END INTERFACE
END MODULE S_ZLAUNHR_COL_GETRFNP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAUU2
   INTERFACE
      SUBROUTINE ZLAUU2(Uplo,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZLAUU2
   END INTERFACE
END MODULE S_ZLAUU2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLAUUM
   INTERFACE
      SUBROUTINE ZLAUUM(Uplo,N,A,Lda,Info)
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
      END SUBROUTINE ZLAUUM
   END INTERFACE
END MODULE S_ZLAUUM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZLA_WWADDW
   INTERFACE
      SUBROUTINE ZLA_WWADDW(N,X,Y,W)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: W
      END SUBROUTINE ZLA_WWADDW
   END INTERFACE
END MODULE S_ZLA_WWADDW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPBCON
   INTERFACE
      SUBROUTINE ZPBCON(Uplo,N,Kd,Ab,Ldab,Anorm,Rcond,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPBCON
   END INTERFACE
END MODULE S_ZPBCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPBEQU
   INTERFACE
      SUBROUTINE ZPBEQU(Uplo,N,Kd,Ab,Ldab,S,Scond,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) :: Scond
      REAL(R8KIND) , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPBEQU
   END INTERFACE
END MODULE S_ZPBEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPBRFS
   INTERFACE
      SUBROUTINE ZPBRFS(Uplo,N,Kd,Nrhs,Ab,Ldab,Afb,Ldafb,B,Ldb,X,Ldx,   &
     &                  Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , THREE = 3.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPBRFS
   END INTERFACE
END MODULE S_ZPBRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPBSTF
   INTERFACE
      SUBROUTINE ZPBSTF(Uplo,N,Kd,Ab,Ldab,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPBSTF
   END INTERFACE
END MODULE S_ZPBSTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPBSV
   INTERFACE
      SUBROUTINE ZPBSV(Uplo,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPBSV
   END INTERFACE
END MODULE S_ZPBSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPBSVX
   INTERFACE
      SUBROUTINE ZPBSVX(Fact,Uplo,N,Kd,Nrhs,Ab,Ldab,Afb,Ldafb,Equed,S,B,&
     &                  Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPBSVX
   END INTERFACE
END MODULE S_ZPBSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPBTF2
   INTERFACE
      SUBROUTINE ZPBTF2(Uplo,N,Kd,Ab,Ldab,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPBTF2
   END INTERFACE
END MODULE S_ZPBTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPBTRF
   INTERFACE
      SUBROUTINE ZPBTRF(Uplo,N,Kd,Ab,Ldab,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER , PARAMETER  ::  NBMAX = 32 , LDWORK = NBMAX + 1
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPBTRF
   END INTERFACE
END MODULE S_ZPBTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPBTRS
   INTERFACE
      SUBROUTINE ZPBTRS(Uplo,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPBTRS
   END INTERFACE
END MODULE S_ZPBTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPFTRF
   INTERFACE
      SUBROUTINE ZPFTRF(Transr,Uplo,N,A,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPFTRF
   END INTERFACE
END MODULE S_ZPFTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPFTRI
   INTERFACE
      SUBROUTINE ZPFTRI(Transr,Uplo,N,A,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.D0,0.D0)
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPFTRI
   END INTERFACE
END MODULE S_ZPFTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPFTRS
   INTERFACE
      SUBROUTINE ZPFTRS(Transr,Uplo,N,Nrhs,A,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(0:*) :: A
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPFTRS
   END INTERFACE
END MODULE S_ZPFTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPOCON
   INTERFACE
      SUBROUTINE ZPOCON(Uplo,N,A,Lda,Anorm,Rcond,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPOCON
   END INTERFACE
END MODULE S_ZPOCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPOEQUB
   INTERFACE
      SUBROUTINE ZPOEQUB(N,A,Lda,S,Scond,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) :: Scond
      REAL(R8KIND) , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPOEQUB
   END INTERFACE
END MODULE S_ZPOEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPOEQU
   INTERFACE
      SUBROUTINE ZPOEQU(N,A,Lda,S,Scond,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) :: Scond
      REAL(R8KIND) , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPOEQU
   END INTERFACE
END MODULE S_ZPOEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPORFS
   INTERFACE
      SUBROUTINE ZPORFS(Uplo,N,Nrhs,A,Lda,Af,Ldaf,B,Ldb,X,Ldx,Ferr,Berr,&
     &                  Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , THREE = 3.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPORFS
   END INTERFACE
END MODULE S_ZPORFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPORFSX
   INTERFACE
      SUBROUTINE ZPORFSX(Uplo,Equed,N,Nrhs,A,Lda,Af,Ldaf,S,B,Ldb,X,Ldx, &
     &                   Rcond,Berr,N_err_bnds,Err_bnds_norm,           &
     &                   Err_bnds_comp,Nparams,Params,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 ,                     &
     &                              ITREF_DEFAULT = 1.0D+0 ,            &
     &                              ITHRESH_DEFAULT = 10.0D+0 ,         &
     &                              COMPONENTWISE_DEFAULT = 1.0D+0 ,    &
     &                              RTHRESH_DEFAULT = 0.5D+0 ,          &
     &                              DZTHRESH_DEFAULT = 0.25D+0
      INTEGER , PARAMETER  ::  LA_LINRX_ITREF_I = 1 ,                   &
     &                         LA_LINRX_ITHRESH_I = 2 ,                 &
     &                         LA_LINRX_CWISE_I = 3 ,                   &
     &                         LA_LINRX_TRUST_I = 1 ,                   &
     &                         LA_LINRX_ERR_I = 2 , LA_LINRX_RCOND_I = 3
      CHARACTER :: Uplo
      CHARACTER :: Equed
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Params
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPORFSX
   END INTERFACE
END MODULE S_ZPORFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPOSV
   INTERFACE
      SUBROUTINE ZPOSV(Uplo,N,Nrhs,A,Lda,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPOSV
   END INTERFACE
END MODULE S_ZPOSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPOSVX
   INTERFACE
      SUBROUTINE ZPOSVX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Equed,S,B,Ldb,X, &
     &                  Ldx,Rcond,Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPOSVX
   END INTERFACE
END MODULE S_ZPOSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPOSVXX
   INTERFACE
      SUBROUTINE ZPOSVXX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Equed,S,B,Ldb,X,&
     &                   Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,Err_bnds_norm,&
     &                   Err_bnds_comp,Nparams,Params,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , INTENT(OUT) :: Rpvgrw
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL(R8KIND) , DIMENSION(*) :: Params
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPOSVXX
   END INTERFACE
END MODULE S_ZPOSVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPOTF2
   INTERFACE
      SUBROUTINE ZPOTF2(Uplo,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPOTF2
   END INTERFACE
END MODULE S_ZPOTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPOTRF2
   INTERFACE
      RECURSIVE SUBROUTINE ZPOTRF2(Uplo,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPOTRF2
   END INTERFACE
END MODULE S_ZPOTRF2
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
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPOTRI
   INTERFACE
      SUBROUTINE ZPOTRI(Uplo,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPOTRI
   END INTERFACE
END MODULE S_ZPOTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPOTRS
   INTERFACE
      SUBROUTINE ZPOTRS(Uplo,N,Nrhs,A,Lda,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPOTRS
   END INTERFACE
END MODULE S_ZPOTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPPCON
   INTERFACE
      SUBROUTINE ZPPCON(Uplo,N,Ap,Anorm,Rcond,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPPCON
   END INTERFACE
END MODULE S_ZPPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPPEQU
   INTERFACE
      SUBROUTINE ZPPEQU(Uplo,N,Ap,S,Scond,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) :: Scond
      REAL(R8KIND) , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPPEQU
   END INTERFACE
END MODULE S_ZPPEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPPRFS
   INTERFACE
      SUBROUTINE ZPPRFS(Uplo,N,Nrhs,Ap,Afp,B,Ldb,X,Ldx,Ferr,Berr,Work,  &
     &                  Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , THREE = 3.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Afp
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPPRFS
   END INTERFACE
END MODULE S_ZPPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPPSV
   INTERFACE
      SUBROUTINE ZPPSV(Uplo,N,Nrhs,Ap,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPPSV
   END INTERFACE
END MODULE S_ZPPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPPSVX
   INTERFACE
      SUBROUTINE ZPPSVX(Fact,Uplo,N,Nrhs,Ap,Afp,Equed,S,B,Ldb,X,Ldx,    &
     &                  Rcond,Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Afp
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPPSVX
   END INTERFACE
END MODULE S_ZPPSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPPTRF
   INTERFACE
      SUBROUTINE ZPPTRF(Uplo,N,Ap,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPPTRF
   END INTERFACE
END MODULE S_ZPPTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPPTRI
   INTERFACE
      SUBROUTINE ZPPTRI(Uplo,N,Ap,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPPTRI
   END INTERFACE
END MODULE S_ZPPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPPTRS
   INTERFACE
      SUBROUTINE ZPPTRS(Uplo,N,Nrhs,Ap,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPPTRS
   END INTERFACE
END MODULE S_ZPPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPSTF2
   INTERFACE
      SUBROUTINE ZPSTF2(Uplo,N,A,Lda,Piv,Rank,Tol,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(N) :: Piv
      INTEGER , INTENT(OUT) :: Rank
      REAL(R8KIND) , INTENT(IN) :: Tol
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(2*N) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPSTF2
   END INTERFACE
END MODULE S_ZPSTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPSTRF
   INTERFACE
      SUBROUTINE ZPSTRF(Uplo,N,A,Lda,Piv,Rank,Tol,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(N) :: Piv
      INTEGER :: Rank
      REAL(R8KIND) :: Tol
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(2*N) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPSTRF
   END INTERFACE
END MODULE S_ZPSTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPTCON
   INTERFACE
      SUBROUTINE ZPTCON(N,D,E,Anorm,Rcond,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: E
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPTCON
   END INTERFACE
END MODULE S_ZPTCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPTEQR
   INTERFACE
      SUBROUTINE ZPTEQR(Compz,N,D,E,Z,Ldz,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Compz
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPTEQR
   END INTERFACE
END MODULE S_ZPTEQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPTRFS
   INTERFACE
      SUBROUTINE ZPTRFS(Uplo,N,Nrhs,D,E,Df,Ef,B,Ldb,X,Ldx,Ferr,Berr,    &
     &                  Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0 , THREE = 3.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Df
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ef
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPTRFS
   END INTERFACE
END MODULE S_ZPTRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPTSV
   INTERFACE
      SUBROUTINE ZPTSV(N,Nrhs,D,E,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPTSV
   END INTERFACE
END MODULE S_ZPTSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPTSVX
   INTERFACE
      SUBROUTINE ZPTSVX(Fact,N,Nrhs,D,E,Df,Ef,B,Ldb,X,Ldx,Rcond,Ferr,   &
     &                  Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Fact
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Df
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ef
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPTSVX
   END INTERFACE
END MODULE S_ZPTSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPTTRF
   INTERFACE
      SUBROUTINE ZPTTRF(N,D,E,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER :: Info
      END SUBROUTINE ZPTTRF
   END INTERFACE
END MODULE S_ZPTTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPTTRS
   INTERFACE
      SUBROUTINE ZPTTRS(Uplo,N,Nrhs,D,E,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZPTTRS
   END INTERFACE
END MODULE S_ZPTTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZPTTS2
   INTERFACE
      SUBROUTINE ZPTTS2(Iuplo,N,Nrhs,D,E,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: Iuplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      END SUBROUTINE ZPTTS2
   END INTERFACE
END MODULE S_ZPTTS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZROT
   INTERFACE
      SUBROUTINE ZROT(N,Cx,Incx,Cy,Incy,C,S)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(IN) :: C
      COMPLEX(CX16KIND) , INTENT(IN) :: S
      END SUBROUTINE ZROT
   END INTERFACE
END MODULE S_ZROT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSPCON
   INTERFACE
      SUBROUTINE ZSPCON(Uplo,N,Ap,Ipiv,Anorm,Rcond,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSPCON
   END INTERFACE
END MODULE S_ZSPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSPMV
   INTERFACE
      SUBROUTINE ZSPMV(Uplo,N,Alpha,Ap,X,Incx,Beta,Y,Incy)
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
      END SUBROUTINE ZSPMV
   END INTERFACE
END MODULE S_ZSPMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSPR
   INTERFACE
      SUBROUTINE ZSPR(Uplo,N,Alpha,X,Incx,Ap)
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
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      END SUBROUTINE ZSPR
   END INTERFACE
END MODULE S_ZSPR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSPRFS
   INTERFACE
      SUBROUTINE ZSPRFS(Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Ferr,Berr,  &
     &                  Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , THREE = 3.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSPRFS
   END INTERFACE
END MODULE S_ZSPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSPSV
   INTERFACE
      SUBROUTINE ZSPSV(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSPSV
   END INTERFACE
END MODULE S_ZSPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSPSVX
   INTERFACE
      SUBROUTINE ZSPSVX(Fact,Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Rcond, &
     &                  Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSPSVX
   END INTERFACE
END MODULE S_ZSPSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSPTRF
   INTERFACE
      SUBROUTINE ZSPTRF(Uplo,N,Ap,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSPTRF
   END INTERFACE
END MODULE S_ZSPTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSPTRI
   INTERFACE
      SUBROUTINE ZSPTRI(Uplo,N,Ap,Ipiv,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSPTRI
   END INTERFACE
END MODULE S_ZSPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSPTRS
   INTERFACE
      SUBROUTINE ZSPTRS(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSPTRS
   END INTERFACE
END MODULE S_ZSPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSTEDC
   INTERFACE
      SUBROUTINE ZSTEDC(Compz,N,D,E,Z,Ldz,Work,Lwork,Rwork,Lrwork,Iwork,&
     &                  Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0
      CHARACTER :: Compz
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSTEDC
   END INTERFACE
END MODULE S_ZSTEDC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSTEGR
   INTERFACE
      SUBROUTINE ZSTEGR(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,  &
     &                  Isuppz,Work,Lwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Jobz
      CHARACTER :: Range
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) :: Abstol
      INTEGER :: M
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER :: Info
      END SUBROUTINE ZSTEGR
   END INTERFACE
END MODULE S_ZSTEGR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSTEIN
   INTERFACE
      SUBROUTINE ZSTEIN(N,D,E,M,W,Iblock,Isplit,Z,Ldz,Work,Iwork,Ifail, &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TEN = 1.0D+1 , ODM3 = 1.0D-3 ,      &
     &                              ODM1 = 1.0D-1
      INTEGER , PARAMETER  ::  MAXITS = 5 , EXTRA = 2
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) :: M
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: W
      INTEGER , INTENT(IN) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(IN) , DIMENSION(*) :: Isplit
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSTEIN
   END INTERFACE
END MODULE S_ZSTEIN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSTEMR
   INTERFACE
      SUBROUTINE ZSTEMR(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,M,W,Z,Ldz,Nzc,     &
     &                  Isuppz,Tryrac,Work,Lwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              FOUR = 4.0D0 , MINRGP = 1.0D-3
      CHARACTER :: Jobz
      CHARACTER :: Range
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(IN) :: Nzc
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isuppz
      LOGICAL , INTENT(INOUT) :: Tryrac
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSTEMR
   END INTERFACE
END MODULE S_ZSTEMR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSTEQR
   INTERFACE
      SUBROUTINE ZSTEQR(Compz,N,D,E,Z,Ldz,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , THREE = 3.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D0,0.0D0) ,        &
     &                 CONE = (1.0D0,0.0D0)
      INTEGER , PARAMETER  ::  MAXIT = 30
      CHARACTER :: Compz
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSTEQR
   END INTERFACE
END MODULE S_ZSTEQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYCON_3
   INTERFACE
      SUBROUTINE ZSYCON_3(Uplo,N,A,Lda,E,Ipiv,Anorm,Rcond,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYCON_3
   END INTERFACE
END MODULE S_ZSYCON_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYCON
   INTERFACE
      SUBROUTINE ZSYCON(Uplo,N,A,Lda,Ipiv,Anorm,Rcond,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYCON
   END INTERFACE
END MODULE S_ZSYCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYCON_ROOK
   INTERFACE
      SUBROUTINE ZSYCON_ROOK(Uplo,N,A,Lda,Ipiv,Anorm,Rcond,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYCON_ROOK
   END INTERFACE
END MODULE S_ZSYCON_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYCONV
   INTERFACE
      SUBROUTINE ZSYCONV(Uplo,Way,N,A,Lda,Ipiv,E,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYCONV
   END INTERFACE
END MODULE S_ZSYCONV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYCONVF
   INTERFACE
      SUBROUTINE ZSYCONVF(Uplo,Way,N,A,Lda,E,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYCONVF
   END INTERFACE
END MODULE S_ZSYCONVF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYCONVF_ROOK
   INTERFACE
      SUBROUTINE ZSYCONVF_ROOK(Uplo,Way,N,A,Lda,E,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYCONVF_ROOK
   END INTERFACE
END MODULE S_ZSYCONVF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYEQUB
   INTERFACE
      SUBROUTINE ZSYEQUB(Uplo,N,A,Lda,S,Scond,Amax,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0
      INTEGER , PARAMETER  ::  MAX_ITER = 100
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) :: Scond
      REAL(R8KIND) , INTENT(INOUT) :: Amax
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYEQUB
   END INTERFACE
END MODULE S_ZSYEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYMV
   INTERFACE
      SUBROUTINE ZSYMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
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
      END SUBROUTINE ZSYMV
   END INTERFACE
END MODULE S_ZSYMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYR
   INTERFACE
      SUBROUTINE ZSYR(Uplo,N,Alpha,X,Incx,A,Lda)
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
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE ZSYR
   END INTERFACE
END MODULE S_ZSYR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYRFS
   INTERFACE
      SUBROUTINE ZSYRFS(Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,Ferr,&
     &                  Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , THREE = 3.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYRFS
   END INTERFACE
END MODULE S_ZSYRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYRFSX
   INTERFACE
      SUBROUTINE ZSYRFSX(Uplo,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,S,B,Ldb,X,&
     &                   Ldx,Rcond,Berr,N_err_bnds,Err_bnds_norm,       &
     &                   Err_bnds_comp,Nparams,Params,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 ,                     &
     &                              ITREF_DEFAULT = 1.0D+0 ,            &
     &                              ITHRESH_DEFAULT = 10.0D+0 ,         &
     &                              COMPONENTWISE_DEFAULT = 1.0D+0 ,    &
     &                              RTHRESH_DEFAULT = 0.5D+0 ,          &
     &                              DZTHRESH_DEFAULT = 0.25D+0
      INTEGER , PARAMETER  ::  LA_LINRX_ITREF_I = 1 ,                   &
     &                         LA_LINRX_ITHRESH_I = 2 ,                 &
     &                         LA_LINRX_CWISE_I = 3 ,                   &
     &                         LA_LINRX_TRUST_I = 1 ,                   &
     &                         LA_LINRX_ERR_I = 2 , LA_LINRX_RCOND_I = 3
      CHARACTER :: Uplo
      CHARACTER :: Equed
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Params
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYRFSX
   END INTERFACE
END MODULE S_ZSYRFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYSV_AA_2STAGE
   INTERFACE
      SUBROUTINE ZSYSV_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B, &
     &                           Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tb
      INTEGER :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYSV_AA_2STAGE
   END INTERFACE
END MODULE S_ZSYSV_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYSV_AA
   INTERFACE
      SUBROUTINE ZSYSV_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYSV_AA
   END INTERFACE
END MODULE S_ZSYSV_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYSV
   INTERFACE
      SUBROUTINE ZSYSV(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYSV
   END INTERFACE
END MODULE S_ZSYSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYSV_RK
   INTERFACE
      SUBROUTINE ZSYSV_RK(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Work,Lwork,    &
     &                    Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYSV_RK
   END INTERFACE
END MODULE S_ZSYSV_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYSV_ROOK
   INTERFACE
      SUBROUTINE ZSYSV_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,    &
     &                      Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYSV_ROOK
   END INTERFACE
END MODULE S_ZSYSV_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYSVX
   INTERFACE
      SUBROUTINE ZSYSVX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,&
     &                  Rcond,Ferr,Berr,Work,Lwork,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYSVX
   END INTERFACE
END MODULE S_ZSYSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYSVXX
   INTERFACE
      SUBROUTINE ZSYSVXX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,S,B, &
     &                   Ldb,X,Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,        &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Fact
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: S
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , INTENT(OUT) :: Rpvgrw
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL(R8KIND) , DIMENSION(*) :: Params
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYSVXX
   END INTERFACE
END MODULE S_ZSYSVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYSWAPR
   INTERFACE
      SUBROUTINE ZSYSWAPR(Uplo,N,A,Lda,I1,I2)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,N) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) :: I1
      INTEGER , INTENT(IN) :: I2
      END SUBROUTINE ZSYSWAPR
   END INTERFACE
END MODULE S_ZSYSWAPR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTF2
   INTERFACE
      SUBROUTINE ZSYTF2(Uplo,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTF2
   END INTERFACE
END MODULE S_ZSYTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTF2_RK
   INTERFACE
      SUBROUTINE ZSYTF2_RK(Uplo,N,A,Lda,E,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTF2_RK
   END INTERFACE
END MODULE S_ZSYTF2_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTF2_ROOK
   INTERFACE
      SUBROUTINE ZSYTF2_ROOK(Uplo,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTF2_ROOK
   END INTERFACE
END MODULE S_ZSYTF2_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRF_AA_2STAGE
   INTERFACE
      SUBROUTINE ZSYTRF_AA_2STAGE(Uplo,N,A,Lda,Tb,Ltb,Ipiv,Ipiv2,Work,  &
     &                            Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRF_AA_2STAGE
   END INTERFACE
END MODULE S_ZSYTRF_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRF_AA
   INTERFACE
      SUBROUTINE ZSYTRF_AA(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRF_AA
   END INTERFACE
END MODULE S_ZSYTRF_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRF
   INTERFACE
      SUBROUTINE ZSYTRF(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRF
   END INTERFACE
END MODULE S_ZSYTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRF_RK
   INTERFACE
      SUBROUTINE ZSYTRF_RK(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRF_RK
   END INTERFACE
END MODULE S_ZSYTRF_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRF_ROOK
   INTERFACE
      SUBROUTINE ZSYTRF_ROOK(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRF_ROOK
   END INTERFACE
END MODULE S_ZSYTRF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRI2
   INTERFACE
      SUBROUTINE ZSYTRI2(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRI2
   END INTERFACE
END MODULE S_ZSYTRI2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRI2X
   INTERFACE
      SUBROUTINE ZSYTRI2X(Uplo,N,A,Lda,Ipiv,Work,Nb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRI2X
   END INTERFACE
END MODULE S_ZSYTRI2X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRI_3
   INTERFACE
      SUBROUTINE ZSYTRI_3(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRI_3
   END INTERFACE
END MODULE S_ZSYTRI_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRI_3X
   INTERFACE
      SUBROUTINE ZSYTRI_3X(Uplo,N,A,Lda,E,Ipiv,Work,Nb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRI_3X
   END INTERFACE
END MODULE S_ZSYTRI_3X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRI
   INTERFACE
      SUBROUTINE ZSYTRI(Uplo,N,A,Lda,Ipiv,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRI
   END INTERFACE
END MODULE S_ZSYTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRI_ROOK
   INTERFACE
      SUBROUTINE ZSYTRI_ROOK(Uplo,N,A,Lda,Ipiv,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRI_ROOK
   END INTERFACE
END MODULE S_ZSYTRI_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRS2
   INTERFACE
      SUBROUTINE ZSYTRS2(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRS2
   END INTERFACE
END MODULE S_ZSYTRS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRS_3
   INTERFACE
      SUBROUTINE ZSYTRS_3(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRS_3
   END INTERFACE
END MODULE S_ZSYTRS_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRS_AA_2STAGE
   INTERFACE
      SUBROUTINE ZSYTRS_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B,&
     &                            Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRS_AA_2STAGE
   END INTERFACE
END MODULE S_ZSYTRS_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRS_AA
   INTERFACE
      SUBROUTINE ZSYTRS_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRS_AA
   END INTERFACE
END MODULE S_ZSYTRS_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRS
   INTERFACE
      SUBROUTINE ZSYTRS(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRS
   END INTERFACE
END MODULE S_ZSYTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZSYTRS_ROOK
   INTERFACE
      SUBROUTINE ZSYTRS_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZSYTRS_ROOK
   END INTERFACE
END MODULE S_ZSYTRS_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTBCON
   INTERFACE
      SUBROUTINE ZTBCON(Norm,Uplo,Diag,N,Kd,Ab,Ldab,Rcond,Work,Rwork,   &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTBCON
   END INTERFACE
END MODULE S_ZTBCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTBRFS
   INTERFACE
      SUBROUTINE ZTBRFS(Uplo,Trans,Diag,N,Kd,Nrhs,Ab,Ldab,B,Ldb,X,Ldx,  &
     &                  Ferr,Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTBRFS
   END INTERFACE
END MODULE S_ZTBRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTBTRS
   INTERFACE
      SUBROUTINE ZTBTRS(Uplo,Trans,Diag,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTBTRS
   END INTERFACE
END MODULE S_ZTBTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTFSM
   INTERFACE
      SUBROUTINE ZTFSM(Transr,Side,Uplo,Trans,Diag,M,N,Alpha,A,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Transr
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) :: Alpha
      COMPLEX(CX16KIND) , DIMENSION(0:*) :: A
      COMPLEX(CX16KIND) , DIMENSION(0:Ldb-1,0:*) :: B
      INTEGER :: Ldb
      END SUBROUTINE ZTFSM
   END INTERFACE
END MODULE S_ZTFSM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTFTRI
   INTERFACE
      SUBROUTINE ZTFTRI(Transr,Uplo,Diag,N,A,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Transr
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTFTRI
   END INTERFACE
END MODULE S_ZTFTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTFTTP
   INTERFACE
      SUBROUTINE ZTFTTP(Transr,Uplo,N,Arf,Ap,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(0:*) :: Arf
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(0:*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTFTTP
   END INTERFACE
END MODULE S_ZTFTTP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTFTTR
   INTERFACE
      SUBROUTINE ZTFTTR(Transr,Uplo,N,Arf,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(0:*) :: Arf
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(0:Lda-1,0:*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTFTTR
   END INTERFACE
END MODULE S_ZTFTTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTGEVC
   INTERFACE
      SUBROUTINE ZTGEVC(Side,Howmny,Select,N,S,Lds,P,Ldp,Vl,Ldvl,Vr,    &
     &                  Ldvr,Mm,M,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lds,*) :: S
      INTEGER , INTENT(IN) :: Lds
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldp,*) :: P
      INTEGER , INTENT(IN) :: Ldp
      COMPLEX(CX16KIND) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX(CX16KIND) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(OUT) :: M
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTGEVC
   END INTERFACE
END MODULE S_ZTGEVC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTGEX2
   INTERFACE
      SUBROUTINE ZTGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,J1,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  TWENTY = 2.0D+1
      INTEGER , PARAMETER  ::  LDST = 2
      LOGICAL , PARAMETER  ::  WANDS = .TRUE.
      LOGICAL , INTENT(IN) :: Wantq
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER , INTENT(IN) :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      INTEGER , INTENT(IN) :: J1
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE ZTGEX2
   END INTERFACE
END MODULE S_ZTGEX2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTGEXC
   INTERFACE
      SUBROUTINE ZTGEXC(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,Ifst,Ilst,&
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
      LOGICAL :: Wantq
      LOGICAL :: Wantz
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(IN) :: Ifst
      INTEGER , INTENT(INOUT) :: Ilst
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTGEXC
   END INTERFACE
END MODULE S_ZTGEXC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTGSEN
   INTERFACE
      SUBROUTINE ZTGSEN(Ijob,Wantq,Wantz,Select,N,A,Lda,B,Ldb,Alpha,    &
     &                  Beta,Q,Ldq,Z,Ldz,M,Pl,Pr,Dif,Work,Lwork,Iwork,  &
     &                  Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  IDIFJB = 3
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , INTENT(IN) :: Ijob
      LOGICAL :: Wantq
      LOGICAL :: Wantz
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Alpha
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) :: Pl
      REAL(R8KIND) , INTENT(INOUT) :: Pr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dif
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTGSEN
   END INTERFACE
END MODULE S_ZTGSEN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTGSJA
   INTERFACE
      SUBROUTINE ZTGSJA(Jobu,Jobv,Jobq,M,P,N,K,L,A,Lda,B,Ldb,Tola,Tolb, &
     &                  Alpha,Beta,U,Ldu,V,Ldv,Q,Ldq,Work,Ncycle,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXIT = 40
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      CHARACTER :: Jobq
      INTEGER :: M
      INTEGER :: P
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER :: L
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Tola
      REAL(R8KIND) , INTENT(IN) :: Tolb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alpha
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) :: Ncycle
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTGSJA
   END INTERFACE
END MODULE S_ZTGSJA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTGSNA
   INTERFACE
      SUBROUTINE ZTGSNA(Job,Howmny,Select,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,&
     &                  S,Dif,Mm,M,Work,Lwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , PARAMETER  ::  IDIFJB = 3
      CHARACTER :: Job
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      COMPLEX(CX16KIND) , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: S
      REAL(R8KIND) , DIMENSION(*) :: Dif
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTGSNA
   END INTERFACE
END MODULE S_ZTGSNA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTGSY2
   INTERFACE
      SUBROUTINE ZTGSY2(Trans,Ijob,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F, &
     &                  Ldf,Scale,Rdsum,Rdscal,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , PARAMETER  ::  LDZ = 2
      CHARACTER :: Trans
      INTEGER :: Ijob
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(Ldd,*) :: D
      INTEGER , INTENT(IN) :: Ldd
      COMPLEX(CX16KIND) , DIMENSION(Lde,*) :: E
      INTEGER :: Lde
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) :: Rdsum
      REAL(R8KIND) :: Rdscal
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTGSY2
   END INTERFACE
END MODULE S_ZTGSY2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTGSYL
   INTERFACE
      SUBROUTINE ZTGSYL(Trans,Ijob,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F, &
     &                  Ldf,Scale,Dif,Work,Lwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: Ijob
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(Ldd,*) :: D
      INTEGER :: Ldd
      COMPLEX(CX16KIND) , DIMENSION(Lde,*) :: E
      INTEGER :: Lde
      COMPLEX(CX16KIND) , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(OUT) :: Dif
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTGSYL
   END INTERFACE
END MODULE S_ZTGSYL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPCON
   INTERFACE
      SUBROUTINE ZTPCON(Norm,Uplo,Diag,N,Ap,Rcond,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTPCON
   END INTERFACE
END MODULE S_ZTPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPLQT2
   INTERFACE
      SUBROUTINE ZTPLQT2(M,N,L,A,Lda,B,Ldb,T,Ldt,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: L
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTPLQT2
   END INTERFACE
END MODULE S_ZTPLQT2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPLQT
   INTERFACE
      SUBROUTINE ZTPLQT(M,N,L,Mb,A,Lda,B,Ldb,T,Ldt,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Mb
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTPLQT
   END INTERFACE
END MODULE S_ZTPLQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPMLQT
   INTERFACE
      SUBROUTINE ZTPMLQT(Side,Trans,M,N,K,L,Mb,V,Ldv,T,Ldt,A,Lda,B,Ldb, &
     &                   Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Mb
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTPMLQT
   END INTERFACE
END MODULE S_ZTPMLQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPMQRT
   INTERFACE
      SUBROUTINE ZTPMQRT(Side,Trans,M,N,K,L,Nb,V,Ldv,T,Ldt,A,Lda,B,Ldb, &
     &                   Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTPMQRT
   END INTERFACE
END MODULE S_ZTPMQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPQRT2
   INTERFACE
      SUBROUTINE ZTPQRT2(M,N,L,A,Lda,B,Ldb,T,Ldt,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0,0.0) ,              &
     &                 ZERO = (0.0,0.0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: L
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTPQRT2
   END INTERFACE
END MODULE S_ZTPQRT2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPQRT
   INTERFACE
      SUBROUTINE ZTPQRT(M,N,L,Nb,A,Lda,B,Ldb,T,Ldt,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTPQRT
   END INTERFACE
END MODULE S_ZTPQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPRFB
   INTERFACE
      SUBROUTINE ZTPRFB(Side,Trans,Direct,Storev,M,N,K,L,V,Ldv,T,Ldt,A, &
     &                  Lda,B,Ldb,Work,Ldwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0,0.0) ,              &
     &                 ZERO = (0.0,0.0)
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: L
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE ZTPRFB
   END INTERFACE
END MODULE S_ZTPRFB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPRFS
   INTERFACE
      SUBROUTINE ZTPRFS(Uplo,Trans,Diag,N,Nrhs,Ap,B,Ldb,X,Ldx,Ferr,Berr,&
     &                  Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTPRFS
   END INTERFACE
END MODULE S_ZTPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPTRI
   INTERFACE
      SUBROUTINE ZTPTRI(Uplo,Diag,N,Ap,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTPTRI
   END INTERFACE
END MODULE S_ZTPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPTRS
   INTERFACE
      SUBROUTINE ZTPTRS(Uplo,Trans,Diag,N,Nrhs,Ap,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTPTRS
   END INTERFACE
END MODULE S_ZTPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPTTF
   INTERFACE
      SUBROUTINE ZTPTTF(Transr,Uplo,N,Ap,Arf,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(0:*) :: Ap
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(0:*) :: Arf
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTPTTF
   END INTERFACE
END MODULE S_ZTPTTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTPTTR
   INTERFACE
      SUBROUTINE ZTPTTR(Uplo,N,Ap,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTPTTR
   END INTERFACE
END MODULE S_ZTPTTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRCON
   INTERFACE
      SUBROUTINE ZTRCON(Norm,Uplo,Diag,N,A,Lda,Rcond,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTRCON
   END INTERFACE
END MODULE S_ZTRCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTREVC3
   INTERFACE
      SUBROUTINE ZTREVC3(Side,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,Mm, &
     &                   M,Work,Lwork,Rwork,Lrwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      INTEGER , PARAMETER  ::  NBMIN = 8 , NBMAX = 128
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTREVC3
   END INTERFACE
END MODULE S_ZTREVC3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTREVC
   INTERFACE
      SUBROUTINE ZTREVC(Side,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,Mm,M,&
     &                  Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CMZERO = (0.0D+0,0.0D+0) ,     &
     &                 CMONE = (1.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTREVC
   END INTERFACE
END MODULE S_ZTREVC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTREXC
   INTERFACE
      SUBROUTINE ZTREXC(Compq,N,T,Ldt,Q,Ldq,Ifst,Ilst,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Compq
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER , INTENT(IN) :: Ldq
      INTEGER , INTENT(IN) :: Ifst
      INTEGER , INTENT(IN) :: Ilst
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTREXC
   END INTERFACE
END MODULE S_ZTREXC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRRFS
   INTERFACE
      SUBROUTINE ZTRRFS(Uplo,Trans,Diag,N,Nrhs,A,Lda,B,Ldb,X,Ldx,Ferr,  &
     &                  Berr,Work,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Berr
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTRRFS
   END INTERFACE
END MODULE S_ZTRRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRSEN
   INTERFACE
      SUBROUTINE ZTRSEN(Job,Compq,Select,N,T,Ldt,Q,Ldq,W,M,S,Sep,Work,  &
     &                  Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Job
      CHARACTER :: Compq
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: W
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(OUT) :: S
      REAL(R8KIND) , INTENT(OUT) :: Sep
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTRSEN
   END INTERFACE
END MODULE S_ZTRSEN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRSNA
   INTERFACE
      SUBROUTINE ZTRSNA(Job,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,S,Sep,&
     &                  Mm,M,Work,Ldwork,Rwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D0 + 0
      CHARACTER :: Job
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      COMPLEX(CX16KIND) , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Sep
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTRSNA
   END INTERFACE
END MODULE S_ZTRSNA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRSYL
   INTERFACE
      SUBROUTINE ZTRSYL(Trana,Tranb,Isgn,M,N,A,Lda,B,Ldb,C,Ldc,Scale,   &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Trana
      CHARACTER :: Tranb
      INTEGER , INTENT(IN) :: Isgn
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTRSYL
   END INTERFACE
END MODULE S_ZTRSYL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRTI2
   INTERFACE
      SUBROUTINE ZTRTI2(Uplo,Diag,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTRTI2
   END INTERFACE
END MODULE S_ZTRTI2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRTRI
   INTERFACE
      SUBROUTINE ZTRTRI(Uplo,Diag,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTRTRI
   END INTERFACE
END MODULE S_ZTRTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRTRS
   INTERFACE
      SUBROUTINE ZTRTRS(Uplo,Trans,Diag,N,Nrhs,A,Lda,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTRTRS
   END INTERFACE
END MODULE S_ZTRTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRTTF
   INTERFACE
      SUBROUTINE ZTRTTF(Transr,Uplo,N,A,Lda,Arf,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(0:Lda-1,0:*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(0:*) :: Arf
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTRTTF
   END INTERFACE
END MODULE S_ZTRTTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTRTTP
   INTERFACE
      SUBROUTINE ZTRTTP(Uplo,N,A,Lda,Ap,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTRTTP
   END INTERFACE
END MODULE S_ZTRTTP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZTZRZF
   INTERFACE
      SUBROUTINE ZTZRZF(M,N,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZTZRZF
   END INTERFACE
END MODULE S_ZTZRZF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNBDB1
   INTERFACE
      SUBROUTINE ZUNBDB1(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D0,0.0D0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Phi
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup1
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tauq1
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNBDB1
   END INTERFACE
END MODULE S_ZUNBDB1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNBDB2
   INTERFACE
      SUBROUTINE ZUNBDB2(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  NEGONE = (-1.0D0,0.0D0) ,      &
     &                 ONE = (1.0D0,0.0D0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Theta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Phi
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup1
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tauq1
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNBDB2
   END INTERFACE
END MODULE S_ZUNBDB2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNBDB3
   INTERFACE
      SUBROUTINE ZUNBDB3(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D0,0.0D0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Theta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Phi
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup1
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tauq1
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNBDB3
   END INTERFACE
END MODULE S_ZUNBDB3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNBDB4
   INTERFACE
      SUBROUTINE ZUNBDB4(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Phantom,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  NEGONE = (-1.0D0,0.0D0) ,      &
     &                 ONE = (1.0D0,0.0D0) , ZERO = (0.0D0,0.0D0)
      INTEGER , INTENT(IN) :: M
      INTEGER :: P
      INTEGER :: Q
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Phi
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup1
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tauq1
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Phantom
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNBDB4
   END INTERFACE
END MODULE S_ZUNBDB4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNBDB5
   INTERFACE
      SUBROUTINE ZUNBDB5(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,&
     &                   Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D0,0.0D0) ,          &
     &                 ZERO = (0.0D0,0.0D0)
      INTEGER :: M1
      INTEGER :: M2
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: X1
      INTEGER :: Incx1
      COMPLEX(CX16KIND) , DIMENSION(*) :: X2
      INTEGER :: Incx2
      COMPLEX(CX16KIND) , DIMENSION(Ldq1,*) :: Q1
      INTEGER :: Ldq1
      COMPLEX(CX16KIND) , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNBDB5
   END INTERFACE
END MODULE S_ZUNBDB5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNBDB6
   INTERFACE
      SUBROUTINE ZUNBDB6(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,&
     &                   Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ALPHASQ = 0.01D0 , REALONE = 1.0D0 ,&
     &                              REALZERO = 0.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  NEGONE = (-1.0D0,0.0D0) ,      &
     &                 ONE = (1.0D0,0.0D0) , ZERO = (0.0D0,0.0D0)
      INTEGER :: M1
      INTEGER :: M2
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: X1
      INTEGER :: Incx1
      COMPLEX(CX16KIND) , DIMENSION(*) :: X2
      INTEGER :: Incx2
      COMPLEX(CX16KIND) , DIMENSION(Ldq1,*) :: Q1
      INTEGER :: Ldq1
      COMPLEX(CX16KIND) , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNBDB6
   END INTERFACE
END MODULE S_ZUNBDB6
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNBDB
   INTERFACE
      SUBROUTINE ZUNBDB(Trans,Signs,M,P,Q,X11,Ldx11,X12,Ldx12,X21,Ldx21,&
     &                  X22,Ldx22,Theta,Phi,Taup1,Taup2,Tauq1,Tauq2,    &
     &                  Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  REALONE = 1.0D0
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D0,0.0D0)
      CHARACTER :: Trans
      CHARACTER :: Signs
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      COMPLEX(CX16KIND) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX(CX16KIND) , DIMENSION(Ldx12,*) :: X12
      INTEGER :: Ldx12
      COMPLEX(CX16KIND) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      COMPLEX(CX16KIND) , DIMENSION(Ldx22,*) :: X22
      INTEGER :: Ldx22
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Phi
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup1
      COMPLEX(CX16KIND) , DIMENSION(*) :: Taup2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tauq1
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tauq2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNBDB
   END INTERFACE
END MODULE S_ZUNBDB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNCSD2BY1
   INTERFACE
      SUBROUTINE ZUNCSD2BY1(Jobu1,Jobu2,Jobv1t,M,P,Q,X11,Ldx11,X21,     &
     &                      Ldx21,Theta,U1,Ldu1,U2,Ldu2,V1t,Ldv1t,Work, &
     &                      Lwork,Rwork,Lrwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D0,0.0D0) ,          &
     &                 ZERO = (0.0D0,0.0D0)
      CHARACTER :: Jobu1
      CHARACTER :: Jobu2
      CHARACTER :: Jobv1t
      INTEGER :: M
      INTEGER :: P
      INTEGER :: Q
      COMPLEX(CX16KIND) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX(CX16KIND) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL(R8KIND) , DIMENSION(*) :: Theta
      COMPLEX(CX16KIND) , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      COMPLEX(CX16KIND) , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      COMPLEX(CX16KIND) , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNCSD2BY1
   END INTERFACE
END MODULE S_ZUNCSD2BY1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNCSD
   INTERFACE
      RECURSIVE SUBROUTINE ZUNCSD(Jobu1,Jobu2,Jobv1t,Jobv2t,Trans,Signs,&
     &                            M,P,Q,X11,Ldx11,X12,Ldx12,X21,Ldx21,  &
     &                            X22,Ldx22,Theta,U1,Ldu1,U2,Ldu2,V1t,  &
     &                            Ldv1t,V2t,Ldv2t,Work,Lwork,Rwork,     &
     &                            Lrwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D0,0.0D0) ,          &
     &                 ZERO = (0.0D0,0.0D0)
      CHARACTER :: Jobu1
      CHARACTER :: Jobu2
      CHARACTER :: Jobv1t
      CHARACTER :: Jobv2t
      CHARACTER :: Trans
      CHARACTER :: Signs
      INTEGER :: M
      INTEGER :: P
      INTEGER :: Q
      COMPLEX(CX16KIND) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX(CX16KIND) , DIMENSION(Ldx12,*) :: X12
      INTEGER :: Ldx12
      COMPLEX(CX16KIND) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      COMPLEX(CX16KIND) , DIMENSION(Ldx22,*) :: X22
      INTEGER :: Ldx22
      REAL(R8KIND) , DIMENSION(*) :: Theta
      COMPLEX(CX16KIND) , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      COMPLEX(CX16KIND) , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      COMPLEX(CX16KIND) , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      COMPLEX(CX16KIND) , DIMENSION(Ldv2t,*) :: V2t
      INTEGER :: Ldv2t
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNCSD
   END INTERFACE
END MODULE S_ZUNCSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNG2L
   INTERFACE
      SUBROUTINE ZUNG2L(M,N,K,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNG2L
   END INTERFACE
END MODULE S_ZUNG2L
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNG2R
   INTERFACE
      SUBROUTINE ZUNG2R(M,N,K,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNG2R
   END INTERFACE
END MODULE S_ZUNG2R
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNGBR
   INTERFACE
      SUBROUTINE ZUNGBR(Vect,M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Vect
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNGBR
   END INTERFACE
END MODULE S_ZUNGBR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNGHR
   INTERFACE
      SUBROUTINE ZUNGHR(N,Ilo,Ihi,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNGHR
   END INTERFACE
END MODULE S_ZUNGHR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNGL2
   INTERFACE
      SUBROUTINE ZUNGL2(M,N,K,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNGL2
   END INTERFACE
END MODULE S_ZUNGL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNGLQ
   INTERFACE
      SUBROUTINE ZUNGLQ(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNGLQ
   END INTERFACE
END MODULE S_ZUNGLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNGQL
   INTERFACE
      SUBROUTINE ZUNGQL(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNGQL
   END INTERFACE
END MODULE S_ZUNGQL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNGQR
   INTERFACE
      SUBROUTINE ZUNGQR(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNGQR
   END INTERFACE
END MODULE S_ZUNGQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNGR2
   INTERFACE
      SUBROUTINE ZUNGR2(M,N,K,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNGR2
   END INTERFACE
END MODULE S_ZUNGR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNGRQ
   INTERFACE
      SUBROUTINE ZUNGRQ(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNGRQ
   END INTERFACE
END MODULE S_ZUNGRQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNGTR
   INTERFACE
      SUBROUTINE ZUNGTR(Uplo,N,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNGTR
   END INTERFACE
END MODULE S_ZUNGTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNGTSQR
   INTERFACE
      SUBROUTINE ZUNGTSQR(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNGTSQR
   END INTERFACE
END MODULE S_ZUNGTSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNGTSQR_ROW
   INTERFACE
      SUBROUTINE ZUNGTSQR_ROW(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Mb
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNGTSQR_ROW
   END INTERFACE
END MODULE S_ZUNGTSQR_ROW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNHR_COL
   INTERFACE
      SUBROUTINE ZUNHR_COL(M,N,Nb,A,Lda,T,Ldt,D,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0) ,       &
     &                 CZERO = (0.0D+0,0.0D+0)
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNHR_COL
   END INTERFACE
END MODULE S_ZUNHR_COL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNM22
   INTERFACE
      SUBROUTINE ZUNM22(Side,Trans,M,N,N1,N2,Q,Ldq,C,Ldc,Work,Lwork,    &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: N1
      INTEGER :: N2
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNM22
   END INTERFACE
END MODULE S_ZUNM22
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNM2L
   INTERFACE
      SUBROUTINE ZUNM2L(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNM2L
   END INTERFACE
END MODULE S_ZUNM2L
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNM2R
   INTERFACE
      SUBROUTINE ZUNM2R(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNM2R
   END INTERFACE
END MODULE S_ZUNM2R
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNMBR
   INTERFACE
      SUBROUTINE ZUNMBR(Vect,Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,     &
     &                  Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Vect
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNMBR
   END INTERFACE
END MODULE S_ZUNMBR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNMHR
   INTERFACE
      SUBROUTINE ZUNMHR(Side,Trans,M,N,Ilo,Ihi,A,Lda,Tau,C,Ldc,Work,    &
     &                  Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNMHR
   END INTERFACE
END MODULE S_ZUNMHR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNML2
   INTERFACE
      SUBROUTINE ZUNML2(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNML2
   END INTERFACE
END MODULE S_ZUNML2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNMLQ
   INTERFACE
      SUBROUTINE ZUNMLQ(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
     &                  Info)
      USE F77KINDS                        
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
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNMLQ
   END INTERFACE
END MODULE S_ZUNMLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNMQL
   INTERFACE
      SUBROUTINE ZUNMQL(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
     &                  Info)
      USE F77KINDS                        
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
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNMQL
   END INTERFACE
END MODULE S_ZUNMQL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNMQR
   INTERFACE
      SUBROUTINE ZUNMQR(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
     &                  Info)
      USE F77KINDS                        
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
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNMQR
   END INTERFACE
END MODULE S_ZUNMQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNMR2
   INTERFACE
      SUBROUTINE ZUNMR2(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNMR2
   END INTERFACE
END MODULE S_ZUNMR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNMR3
   INTERFACE
      SUBROUTINE ZUNMR3(Side,Trans,M,N,K,L,A,Lda,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      INTEGER :: L
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNMR3
   END INTERFACE
END MODULE S_ZUNMR3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNMRQ
   INTERFACE
      SUBROUTINE ZUNMRQ(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
     &                  Info)
      USE F77KINDS                        
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
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNMRQ
   END INTERFACE
END MODULE S_ZUNMRQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNMRZ
   INTERFACE
      SUBROUTINE ZUNMRZ(Side,Trans,M,N,K,L,A,Lda,Tau,C,Ldc,Work,Lwork,  &
     &                  Info)
      USE F77KINDS                        
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
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNMRZ
   END INTERFACE
END MODULE S_ZUNMRZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUNMTR
   INTERFACE
      SUBROUTINE ZUNMTR(Side,Uplo,Trans,M,N,A,Lda,Tau,C,Ldc,Work,Lwork, &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUNMTR
   END INTERFACE
END MODULE S_ZUNMTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUPGTR
   INTERFACE
      SUBROUTINE ZUPGTR(Uplo,N,Ap,Tau,Q,Ldq,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUPGTR
   END INTERFACE
END MODULE S_ZUPGTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_ZUPMTR
   INTERFACE
      SUBROUTINE ZUPMTR(Side,Uplo,Trans,M,N,Ap,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE ZUPMTR
   END INTERFACE
END MODULE S_ZUPMTR
