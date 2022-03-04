!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DBBCSD
   INTERFACE
      SUBROUTINE DBBCSD(Jobu1,Jobu2,Jobv1t,Jobv2t,Trans,M,P,Q,Theta,Phi,&
     &                  U1,Ldu1,U2,Ldu2,V1t,Ldv1t,V2t,Ldv2t,B11d,B11e,  &
     &                  B12d,B12e,B21d,B21e,B22d,B22e,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXITR = 6
      REAL(R8KIND) , PARAMETER  ::  HUNDRED = 100.0D0 ,                 &
     &                              MEIGHTH = -0.125D0 , ONE = 1.0D0 ,  &
     &                              TEN = 10.0D0 , ZERO = 0.0D0 ,       &
     &                              NEGONE = -1.0D0 , PIOVER2 =         &
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
      REAL(R8KIND) , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      REAL(R8KIND) , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      REAL(R8KIND) , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      REAL(R8KIND) , DIMENSION(Ldv2t,*) :: V2t
      INTEGER :: Ldv2t
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B11d
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B11e
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B12d
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B12e
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B21d
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B21e
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B22d
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B22e
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DBBCSD
   END INTERFACE
END MODULE S_DBBCSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DBDSDC
   INTERFACE
      SUBROUTINE DBDSDC(Uplo,Compq,N,D,E,U,Ldu,Vt,Ldvt,Q,Iq,Work,Iwork, &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Compq
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL(R8KIND) , DIMENSION(*) :: Q
      INTEGER , DIMENSION(*) :: Iq
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DBDSDC
   END INTERFACE
END MODULE S_DBDSDC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DBDSQR
   INTERFACE
      SUBROUTINE DBDSQR(Uplo,N,Ncvt,Nru,Ncc,D,E,Vt,Ldvt,U,Ldu,C,Ldc,    &
     &                  Work,Info)
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
      REAL(R8KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DBDSQR
   END INTERFACE
END MODULE S_DBDSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DBDSVDX
   INTERFACE
      SUBROUTINE DBDSVDX(Uplo,Jobz,Range,N,D,E,Vl,Vu,Il,Iu,Ns,S,Z,Ldz,  &
     &                   Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TEN = 10.0D0 , HNDRD = 100.0D0 ,    &
     &                              MEIGTH = -0.1250D0 , FUDGE = 2.0D0
      CHARACTER :: Uplo
      CHARACTER :: Jobz
      CHARACTER :: Range
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      INTEGER , INTENT(INOUT) :: Ns
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DBDSVDX
   END INTERFACE
END MODULE S_DBDSVDX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DCOMBSSQ
   INTERFACE
      SUBROUTINE DCOMBSSQ(V1,V2)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(2) :: V1
      REAL(R8KIND) , INTENT(IN) , DIMENSION(2) :: V2
      END SUBROUTINE DCOMBSSQ
   END INTERFACE
END MODULE S_DCOMBSSQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DDISNA
   INTERFACE
      SUBROUTINE DDISNA(Job,M,N,D,Sep,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Job
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Sep
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DDISNA
   END INTERFACE
END MODULE S_DDISNA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBBRD
   INTERFACE
      SUBROUTINE DGBBRD(Vect,M,N,Ncc,Kl,Ku,Ab,Ldab,D,E,Q,Ldq,Pt,Ldpt,C, &
     &                  Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Vect
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Ncc
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , DIMENSION(Ldpt,*) :: Pt
      INTEGER :: Ldpt
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGBBRD
   END INTERFACE
END MODULE S_DGBBRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBCON
   INTERFACE
      SUBROUTINE DGBCON(Norm,N,Kl,Ku,Ab,Ldab,Ipiv,Anorm,Rcond,Work,     &
     &                  Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGBCON
   END INTERFACE
END MODULE S_DGBCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBEQUB
   INTERFACE
      SUBROUTINE DGBEQUB(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Info)
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
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(OUT) :: Rowcnd
      REAL(R8KIND) , INTENT(OUT) :: Colcnd
      REAL(R8KIND) , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGBEQUB
   END INTERFACE
END MODULE S_DGBEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBEQU
   INTERFACE
      SUBROUTINE DGBEQU(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Info)
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
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(OUT) :: Rowcnd
      REAL(R8KIND) , INTENT(OUT) :: Colcnd
      REAL(R8KIND) , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGBEQU
   END INTERFACE
END MODULE S_DGBEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBRFS
   INTERFACE
      SUBROUTINE DGBRFS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv,B,Ldb,&
     &                  X,Ldx,Ferr,Berr,Work,Iwork,Info)
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
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGBRFS
   END INTERFACE
END MODULE S_DGBRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBRFSX
   INTERFACE
      SUBROUTINE DGBRFSX(Trans,Equed,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,    &
     &                   Ipiv,R,C,B,Ldb,X,Ldx,Rcond,Berr,N_err_bnds,    &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ITREF_DEFAULT = 1.0D+0 ,            &
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
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: R
      REAL(R8KIND) , DIMENSION(*) :: C
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Params
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGBRFSX
   END INTERFACE
END MODULE S_DGBRFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBSV
   INTERFACE
      SUBROUTINE DGBSV(N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGBSV
   END INTERFACE
END MODULE S_DGBSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBSVX
   INTERFACE
      SUBROUTINE DGBSVX(Fact,Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv, &
     &                  Equed,R,C,B,Ldb,X,Ldx,Rcond,Ferr,Berr,Work,     &
     &                  Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: R
      REAL(R8KIND) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGBSVX
   END INTERFACE
END MODULE S_DGBSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBSVXX
   INTERFACE
      SUBROUTINE DGBSVXX(Fact,Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,Ipiv,&
     &                   Equed,R,C,B,Ldb,X,Ldx,Rcond,Rpvgrw,Berr,       &
     &                   N_err_bnds,Err_bnds_norm,Err_bnds_comp,Nparams,&
     &                   Params,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , INTENT(OUT) :: Rpvgrw
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL(R8KIND) , DIMENSION(*) :: Params
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGBSVXX
   END INTERFACE
END MODULE S_DGBSVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBTF2
   INTERFACE
      SUBROUTINE DGBTF2(M,N,Kl,Ku,Ab,Ldab,Ipiv,Info)
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
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGBTF2
   END INTERFACE
END MODULE S_DGBTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBTRF
   INTERFACE
      SUBROUTINE DGBTRF(M,N,Kl,Ku,Ab,Ldab,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDWORK = NBMAX + 1
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGBTRF
   END INTERFACE
END MODULE S_DGBTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGBTRS
   INTERFACE
      SUBROUTINE DGBTRS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGBTRS
   END INTERFACE
END MODULE S_DGBTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEBAK
   INTERFACE
      SUBROUTINE DGEBAK(Job,Side,N,Ilo,Ihi,Scale,M,V,Ldv,Info)
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
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEBAK
   END INTERFACE
END MODULE S_DGEBAK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEBAL
   INTERFACE
      SUBROUTINE DGEBAL(Job,N,A,Lda,Ilo,Ihi,Scale,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              SCLFAC = 2.0D+0 , FACTOR = 0.95D+0
      CHARACTER :: Job
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) :: Ilo
      INTEGER , INTENT(OUT) :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Scale
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEBAL
   END INTERFACE
END MODULE S_DGEBAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEBD2
   INTERFACE
      SUBROUTINE DGEBD2(M,N,A,Lda,D,E,Tauq,Taup,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Tauq
      REAL(R8KIND) , DIMENSION(*) :: Taup
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEBD2
   END INTERFACE
END MODULE S_DGEBD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEBRD
   INTERFACE
      SUBROUTINE DGEBRD(M,N,A,Lda,D,E,Tauq,Taup,Work,Lwork,Info)
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
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Tauq
      REAL(R8KIND) , DIMENSION(*) :: Taup
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEBRD
   END INTERFACE
END MODULE S_DGEBRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGECON
   INTERFACE
      SUBROUTINE DGECON(Norm,N,A,Lda,Anorm,Rcond,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Norm
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGECON
   END INTERFACE
END MODULE S_DGECON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEEQUB
   INTERFACE
      SUBROUTINE DGEEQUB(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(OUT) :: Rowcnd
      REAL(R8KIND) , INTENT(OUT) :: Colcnd
      REAL(R8KIND) , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEEQUB
   END INTERFACE
END MODULE S_DGEEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEEQU
   INTERFACE
      SUBROUTINE DGEEQU(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(OUT) :: Rowcnd
      REAL(R8KIND) , INTENT(OUT) :: Colcnd
      REAL(R8KIND) , INTENT(OUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEEQU
   END INTERFACE
END MODULE S_DGEEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEES
   INTERFACE
      SUBROUTINE DGEES(Jobvs,Sort,SELECT,N,A,Lda,Sdim,Wr,Wi,Vs,Ldvs,    &
     &                 Work,Lwork,Bwork,Info)
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
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Sdim
      REAL(R8KIND) , DIMENSION(*) :: Wr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wi
      REAL(R8KIND) , DIMENSION(Ldvs,*) :: Vs
      INTEGER :: Ldvs
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEES
   END INTERFACE
END MODULE S_DGEES
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEESX
   INTERFACE
      SUBROUTINE DGEESX(Jobvs,Sort,SELECT,Sense,N,A,Lda,Sdim,Wr,Wi,Vs,  &
     &                  Ldvs,Rconde,Rcondv,Work,Lwork,Iwork,Liwork,     &
     &                  Bwork,Info)
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
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Sdim
      REAL(R8KIND) , DIMENSION(*) :: Wr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wi
      REAL(R8KIND) , DIMENSION(Ldvs,*) :: Vs
      INTEGER :: Ldvs
      REAL(R8KIND) :: Rconde
      REAL(R8KIND) , INTENT(INOUT) :: Rcondv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEESX
   END INTERFACE
END MODULE S_DGEESX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEEV
   INTERFACE
      SUBROUTINE DGEEV(Jobvl,Jobvr,N,A,Lda,Wr,Wi,Vl,Ldvl,Vr,Ldvr,Work,  &
     &                 Lwork,Info)
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
      REAL(R8KIND) , DIMENSION(*) :: Wr
      REAL(R8KIND) , DIMENSION(*) :: Wi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEEV
   END INTERFACE
END MODULE S_DGEEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEEVX
   INTERFACE
      SUBROUTINE DGEEVX(Balanc,Jobvl,Jobvr,Sense,N,A,Lda,Wr,Wi,Vl,Ldvl, &
     &                  Vr,Ldvr,Ilo,Ihi,Scale,Abnrm,Rconde,Rcondv,Work, &
     &                  Lwork,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Wr
      REAL(R8KIND) , DIMENSION(*) :: Wi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , DIMENSION(*) :: Scale
      REAL(R8KIND) , INTENT(INOUT) :: Abnrm
      REAL(R8KIND) , DIMENSION(*) :: Rconde
      REAL(R8KIND) , DIMENSION(*) :: Rcondv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEEVX
   END INTERFACE
END MODULE S_DGEEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEHD2
   INTERFACE
      SUBROUTINE DGEHD2(N,Ilo,Ihi,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEHD2
   END INTERFACE
END MODULE S_DGEHD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEHRD
   INTERFACE
      SUBROUTINE DGEHRD(N,Ilo,Ihi,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDT = NBMAX + 1 ,           &
     &                         TSIZE = LDT*NBMAX
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEHRD
   END INTERFACE
END MODULE S_DGEHRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEJSV
   INTERFACE
      SUBROUTINE DGEJSV(Joba,Jobu,Jobv,Jobr,Jobt,Jobp,M,N,A,Lda,Sva,U,  &
     &                  Ldu,V,Ldv,Work,Lwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER(1) :: Joba
      CHARACTER(1) :: Jobu
      CHARACTER(1) :: Jobv
      CHARACTER(1) :: Jobr
      CHARACTER(1) :: Jobt
      CHARACTER(1) :: Jobp
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N) :: Sva
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lwork) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEJSV
   END INTERFACE
END MODULE S_DGEJSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGELQ2
   INTERFACE
      SUBROUTINE DGELQ2(M,N,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGELQ2
   END INTERFACE
END MODULE S_DGELQ2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGELQ
   INTERFACE
      SUBROUTINE DGELQ(M,N,A,Lda,T,Tsize,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGELQ
   END INTERFACE
END MODULE S_DGELQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGELQF
   INTERFACE
      SUBROUTINE DGELQF(M,N,A,Lda,Tau,Work,Lwork,Info)
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
      END SUBROUTINE DGELQF
   END INTERFACE
END MODULE S_DGELQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGELQT3
   INTERFACE
      RECURSIVE SUBROUTINE DGELQT3(M,N,A,Lda,T,Ldt,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+00
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGELQT3
   END INTERFACE
END MODULE S_DGELQT3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGELQT
   INTERFACE
      SUBROUTINE DGELQT(M,N,Mb,A,Lda,T,Ldt,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Mb
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGELQT
   END INTERFACE
END MODULE S_DGELQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGELSD
   INTERFACE
      SUBROUTINE DGELSD(M,N,Nrhs,A,Lda,B,Ldb,S,Rcond,Rank,Work,Lwork,   &
     &                  Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) :: Rcond
      INTEGER :: Rank
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGELSD
   END INTERFACE
END MODULE S_DGELSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGELS
   INTERFACE
      SUBROUTINE DGELS(Trans,M,N,Nrhs,A,Lda,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGELS
   END INTERFACE
END MODULE S_DGELS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGELSS
   INTERFACE
      SUBROUTINE DGELSS(M,N,Nrhs,A,Lda,B,Ldb,S,Rcond,Rank,Work,Lwork,   &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGELSS
   END INTERFACE
END MODULE S_DGELSS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGELSY
   INTERFACE
      SUBROUTINE DGELSY(M,N,Nrhs,A,Lda,B,Ldb,Jpvt,Rcond,Rank,Work,Lwork,&
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  IMAX = 1 , IMIN = 2
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
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
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGELSY
   END INTERFACE
END MODULE S_DGELSY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEMLQ
   INTERFACE
      SUBROUTINE DGEMLQ(Side,Trans,M,N,K,A,Lda,T,Tsize,C,Ldc,Work,Lwork,&
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEMLQ
   END INTERFACE
END MODULE S_DGEMLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEMLQT
   INTERFACE
      SUBROUTINE DGEMLQT(Side,Trans,M,N,K,Mb,V,Ldv,T,Ldt,C,Ldc,Work,    &
     &                   Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Mb
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEMLQT
   END INTERFACE
END MODULE S_DGEMLQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEMQR
   INTERFACE
      SUBROUTINE DGEMQR(Side,Trans,M,N,K,A,Lda,T,Tsize,C,Ldc,Work,Lwork,&
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEMQR
   END INTERFACE
END MODULE S_DGEMQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEMQRT
   INTERFACE
      SUBROUTINE DGEMQRT(Side,Trans,M,N,K,Nb,V,Ldv,T,Ldt,C,Ldc,Work,    &
     &                   Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Nb
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEMQRT
   END INTERFACE
END MODULE S_DGEMQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEQL2
   INTERFACE
      SUBROUTINE DGEQL2(M,N,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEQL2
   END INTERFACE
END MODULE S_DGEQL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEQLF
   INTERFACE
      SUBROUTINE DGEQLF(M,N,A,Lda,Tau,Work,Lwork,Info)
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
      END SUBROUTINE DGEQLF
   END INTERFACE
END MODULE S_DGEQLF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEQP3
   INTERFACE
      SUBROUTINE DGEQP3(M,N,A,Lda,Jpvt,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  INB = 1 , INBMIN = 2 , IXOVER = 3
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEQP3
   END INTERFACE
END MODULE S_DGEQP3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEQR2
   INTERFACE
      SUBROUTINE DGEQR2(M,N,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEQR2
   END INTERFACE
END MODULE S_DGEQR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEQR2P
   INTERFACE
      SUBROUTINE DGEQR2P(M,N,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEQR2P
   END INTERFACE
END MODULE S_DGEQR2P
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEQR
   INTERFACE
      SUBROUTINE DGEQR(M,N,A,Lda,T,Tsize,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEQR
   END INTERFACE
END MODULE S_DGEQR
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
MODULE S_DGEQRFP
   INTERFACE
      SUBROUTINE DGEQRFP(M,N,A,Lda,Tau,Work,Lwork,Info)
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
      END SUBROUTINE DGEQRFP
   END INTERFACE
END MODULE S_DGEQRFP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEQRT2
   INTERFACE
      SUBROUTINE DGEQRT2(M,N,A,Lda,T,Ldt,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+00 , ZERO = 0.0D+00
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEQRT2
   END INTERFACE
END MODULE S_DGEQRT2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEQRT3
   INTERFACE
      RECURSIVE SUBROUTINE DGEQRT3(M,N,A,Lda,T,Ldt,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+00
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEQRT3
   END INTERFACE
END MODULE S_DGEQRT3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGEQRT
   INTERFACE
      SUBROUTINE DGEQRT(M,N,Nb,A,Lda,T,Ldt,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      LOGICAL , PARAMETER  ::  USE_RECURSIVE_QR = .TRUE.
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGEQRT
   END INTERFACE
END MODULE S_DGEQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGERFS
   INTERFACE
      SUBROUTINE DGERFS(Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,    &
     &                  Ferr,Berr,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGERFS
   END INTERFACE
END MODULE S_DGERFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGERFSX
   INTERFACE
      SUBROUTINE DGERFSX(Trans,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,R,C,B,   &
     &                   Ldb,X,Ldx,Rcond,Berr,N_err_bnds,Err_bnds_norm, &
     &                   Err_bnds_comp,Nparams,Params,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ITREF_DEFAULT = 1.0D+0 ,            &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: R
      REAL(R8KIND) , DIMENSION(*) :: C
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Params
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGERFSX
   END INTERFACE
END MODULE S_DGERFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGERQ2
   INTERFACE
      SUBROUTINE DGERQ2(M,N,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGERQ2
   END INTERFACE
END MODULE S_DGERQ2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGERQF
   INTERFACE
      SUBROUTINE DGERQF(M,N,A,Lda,Tau,Work,Lwork,Info)
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
      END SUBROUTINE DGERQF
   END INTERFACE
END MODULE S_DGERQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGESC2
   INTERFACE
      SUBROUTINE DGESC2(N,A,Lda,Rhs,Ipiv,Jpiv,Scale)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , TWO = 2.0D+0
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rhs
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Jpiv
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      END SUBROUTINE DGESC2
   END INTERFACE
END MODULE S_DGESC2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGESDD
   INTERFACE
      SUBROUTINE DGESDD(Jobz,M,N,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,Iwork,&
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGESDD
   END INTERFACE
END MODULE S_DGESDD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGESVD
   INTERFACE
      SUBROUTINE DGESVD(Jobu,Jobvt,M,N,A,Lda,S,U,Ldu,Vt,Ldvt,Work,Lwork,&
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobu
      CHARACTER :: Jobvt
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGESVD
   END INTERFACE
END MODULE S_DGESVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGESVDQ
   INTERFACE
      SUBROUTINE DGESVDQ(Joba,Jobp,Jobr,Jobu,Jobv,M,N,A,Lda,S,U,Ldu,V,  &
     &                   Ldv,Numrank,Iwork,Liwork,Work,Lwork,Rwork,     &
     &                   Lrwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Joba
      CHARACTER :: Jobp
      CHARACTER :: Jobr
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(OUT) :: Numrank
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGESVDQ
   END INTERFACE
END MODULE S_DGESVDQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGESVDX
   INTERFACE
      SUBROUTINE DGESVDX(Jobu,Jobvt,Range,M,N,A,Lda,Vl,Vu,Il,Iu,Ns,S,U, &
     &                   Ldu,Vt,Ldvt,Work,Lwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobu
      CHARACTER :: Jobvt
      CHARACTER :: Range
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      INTEGER , INTENT(INOUT) :: Ns
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGESVDX
   END INTERFACE
END MODULE S_DGESVDX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGESV
   INTERFACE
      SUBROUTINE DGESV(N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGESV
   END INTERFACE
END MODULE S_DGESV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGESVJ
   INTERFACE
      SUBROUTINE DGESVJ(Joba,Jobu,Jobv,M,N,A,Lda,Sva,Mv,V,Ldv,Work,     &
     &                  Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , HALF = 0.5D0 ,       &
     &                              ONE = 1.0D0
      INTEGER , PARAMETER  ::  NSWEEP = 30
      CHARACTER(1) :: Joba
      CHARACTER(1) :: Jobu
      CHARACTER(1) :: Jobv
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lwork) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGESVJ
   END INTERFACE
END MODULE S_DGESVJ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGESVX
   INTERFACE
      SUBROUTINE DGESVX(Fact,Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,R,C, &
     &                  B,Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: R
      REAL(R8KIND) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGESVX
   END INTERFACE
END MODULE S_DGESVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGESVXX
   INTERFACE
      SUBROUTINE DGESVXX(Fact,Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,R,C,&
     &                   B,Ldb,X,Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,      &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , INTENT(OUT) :: Rpvgrw
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL(R8KIND) , DIMENSION(*) :: Params
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGESVXX
   END INTERFACE
END MODULE S_DGESVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGETC2
   INTERFACE
      SUBROUTINE DGETC2(N,A,Lda,Ipiv,Jpiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Jpiv
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DGETC2
   END INTERFACE
END MODULE S_DGETC2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGETF2
   INTERFACE
      SUBROUTINE DGETF2(M,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGETF2
   END INTERFACE
END MODULE S_DGETF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGETRF2
   INTERFACE
      RECURSIVE SUBROUTINE DGETRF2(M,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGETRF2
   END INTERFACE
END MODULE S_DGETRF2
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
MODULE S_DGETRI
   INTERFACE
      SUBROUTINE DGETRI(N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGETRI
   END INTERFACE
END MODULE S_DGETRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGETRS
   INTERFACE
      SUBROUTINE DGETRS(Trans,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGETRS
   END INTERFACE
END MODULE S_DGETRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGETSLS
   INTERFACE
      SUBROUTINE DGETSLS(Trans,M,N,Nrhs,A,Lda,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGETSLS
   END INTERFACE
END MODULE S_DGETSLS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGETSQRHRT
   INTERFACE
      SUBROUTINE DGETSQRHRT(M,N,Mb1,Nb1,Nb2,A,Lda,T,Ldt,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb1
      INTEGER , INTENT(IN) :: Nb1
      INTEGER , INTENT(IN) :: Nb2
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGETSQRHRT
   END INTERFACE
END MODULE S_DGETSQRHRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGBAK
   INTERFACE
      SUBROUTINE DGGBAK(Job,Side,N,Ilo,Ihi,Lscale,Rscale,M,V,Ldv,Info)
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
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGBAK
   END INTERFACE
END MODULE S_DGGBAK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGBAL
   INTERFACE
      SUBROUTINE DGGBAL(Job,N,A,Lda,B,Ldb,Ilo,Ihi,Lscale,Rscale,Work,   &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , HALF = 0.5D+0 ,     &
     &                              ONE = 1.0D+0 , THREE = 3.0D+0 ,     &
     &                              SCLFAC = 1.0D+1
      CHARACTER :: Job
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Ilo
      INTEGER , INTENT(INOUT) :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Lscale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rscale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGBAL
   END INTERFACE
END MODULE S_DGGBAL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGES3
   INTERFACE
      SUBROUTINE DGGES3(Jobvsl,Jobvsr,Sort,SELCTG,N,A,Lda,B,Ldb,Sdim,   &
     &                  Alphar,Alphai,Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work,    &
     &                  Lwork,Bwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL(R8KIND) , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      REAL(R8KIND) , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGES3
   END INTERFACE
END MODULE S_DGGES3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGES
   INTERFACE
      SUBROUTINE DGGES(Jobvsl,Jobvsr,Sort,SELCTG,N,A,Lda,B,Ldb,Sdim,    &
     &                 Alphar,Alphai,Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work,     &
     &                 Lwork,Bwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL(R8KIND) , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      REAL(R8KIND) , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGES
   END INTERFACE
END MODULE S_DGGES
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGESX
   INTERFACE
      SUBROUTINE DGGESX(Jobvsl,Jobvsr,Sort,SELCTG,Sense,N,A,Lda,B,Ldb,  &
     &                  Sdim,Alphar,Alphai,Beta,Vsl,Ldvsl,Vsr,Ldvsr,    &
     &                  Rconde,Rcondv,Work,Lwork,Iwork,Liwork,Bwork,    &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      CHARACTER :: Sense
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL(R8KIND) , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      REAL(R8KIND) , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(2) :: Rconde
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(2) :: Rcondv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGESX
   END INTERFACE
END MODULE S_DGGESX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGEV3
   INTERFACE
      SUBROUTINE DGGEV3(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alphar,Alphai,Beta,Vl,&
     &                  Ldvl,Vr,Ldvr,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , DIMENSION(*) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGEV3
   END INTERFACE
END MODULE S_DGGEV3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGEV
   INTERFACE
      SUBROUTINE DGGEV(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alphar,Alphai,Beta,Vl, &
     &                 Ldvl,Vr,Ldvr,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , DIMENSION(*) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGEV
   END INTERFACE
END MODULE S_DGGEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGEVX
   INTERFACE
      SUBROUTINE DGGEVX(Balanc,Jobvl,Jobvr,Sense,N,A,Lda,B,Ldb,Alphar,  &
     &                  Alphai,Beta,Vl,Ldvl,Vr,Ldvr,Ilo,Ihi,Lscale,     &
     &                  Rscale,Abnrm,Bbnrm,Rconde,Rcondv,Work,Lwork,    &
     &                  Iwork,Bwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Balanc
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      CHARACTER :: Sense
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , DIMENSION(*) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , DIMENSION(*) :: Lscale
      REAL(R8KIND) , DIMENSION(*) :: Rscale
      REAL(R8KIND) , INTENT(INOUT) :: Abnrm
      REAL(R8KIND) , INTENT(INOUT) :: Bbnrm
      REAL(R8KIND) , DIMENSION(*) :: Rconde
      REAL(R8KIND) , DIMENSION(*) :: Rcondv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGEVX
   END INTERFACE
END MODULE S_DGGEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGGLM
   INTERFACE
      SUBROUTINE DGGGLM(N,M,P,A,Lda,B,Ldb,D,X,Y,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: N
      INTEGER :: M
      INTEGER :: P
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: X
      REAL(R8KIND) , DIMENSION(*) :: Y
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGGLM
   END INTERFACE
END MODULE S_DGGGLM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGHD3
   INTERFACE
      SUBROUTINE DGGHD3(Compq,Compz,N,Ilo,Ihi,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,  &
     &                  Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGHD3
   END INTERFACE
END MODULE S_DGGHD3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGHRD
   INTERFACE
      SUBROUTINE DGGHRD(Compq,Compz,N,Ilo,Ihi,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,  &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGHRD
   END INTERFACE
END MODULE S_DGGHRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGLSE
   INTERFACE
      SUBROUTINE DGGLSE(M,N,P,A,Lda,B,Ldb,C,D,X,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: P
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: C
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: X
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGLSE
   END INTERFACE
END MODULE S_DGGLSE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGQRF
   INTERFACE
      SUBROUTINE DGGQRF(N,M,P,A,Lda,Taua,B,Ldb,Taub,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: M
      INTEGER :: P
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Taua
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Taub
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGQRF
   END INTERFACE
END MODULE S_DGGQRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGRQF
   INTERFACE
      SUBROUTINE DGGRQF(M,P,N,A,Lda,Taua,B,Ldb,Taub,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: P
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Taua
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Taub
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGRQF
   END INTERFACE
END MODULE S_DGGRQF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGSVD3
   INTERFACE
      SUBROUTINE DGGSVD3(Jobu,Jobv,Jobq,M,N,P,K,L,A,Lda,B,Ldb,Alpha,    &
     &                   Beta,U,Ldu,V,Ldv,Q,Ldq,Work,Lwork,Iwork,Info)
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
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGSVD3
   END INTERFACE
END MODULE S_DGGSVD3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGGSVP3
   INTERFACE
      SUBROUTINE DGGSVP3(Jobu,Jobv,Jobq,M,P,N,A,Lda,B,Ldb,Tola,Tolb,K,L,&
     &                   U,Ldu,V,Ldv,Q,Ldq,Iwork,Tau,Work,Lwork,Info)
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
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGGSVP3
   END INTERFACE
END MODULE S_DGGSVP3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGSVJ0
   INTERFACE
      SUBROUTINE DGSVJ0(Jobv,M,N,A,Lda,D,Sva,Mv,V,Ldv,Eps,Sfmin,Tol,    &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      REAL(R8KIND) , INTENT(IN) :: Eps
      REAL(R8KIND) , INTENT(IN) :: Sfmin
      REAL(R8KIND) , INTENT(IN) :: Tol
      INTEGER , INTENT(IN) :: Nsweep
      REAL(R8KIND) , DIMENSION(Lwork) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGSVJ0
   END INTERFACE
END MODULE S_DGSVJ0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGSVJ1
   INTERFACE
      SUBROUTINE DGSVJ1(Jobv,M,N,N1,A,Lda,D,Sva,Mv,V,Ldv,Eps,Sfmin,Tol, &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      REAL(R8KIND) , INTENT(IN) :: Eps
      REAL(R8KIND) , INTENT(IN) :: Sfmin
      REAL(R8KIND) , INTENT(IN) :: Tol
      INTEGER , INTENT(IN) :: Nsweep
      REAL(R8KIND) , DIMENSION(Lwork) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGSVJ1
   END INTERFACE
END MODULE S_DGSVJ1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGTCON
   INTERFACE
      SUBROUTINE DGTCON(Norm,N,Dl,D,Du,Du2,Ipiv,Anorm,Rcond,Work,Iwork, &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Norm
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: Dl
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: Du
      REAL(R8KIND) , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGTCON
   END INTERFACE
END MODULE S_DGTCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGTRFS
   INTERFACE
      SUBROUTINE DGTRFS(Trans,N,Nrhs,Dl,D,Du,Dlf,Df,Duf,Du2,Ipiv,B,Ldb, &
     &                  X,Ldx,Ferr,Berr,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(*) :: Dl
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: Du
      REAL(R8KIND) , DIMENSION(*) :: Dlf
      REAL(R8KIND) , DIMENSION(*) :: Df
      REAL(R8KIND) , DIMENSION(*) :: Duf
      REAL(R8KIND) , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGTRFS
   END INTERFACE
END MODULE S_DGTRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGTSV
   INTERFACE
      SUBROUTINE DGTSV(N,Nrhs,Dl,D,Du,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Du
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGTSV
   END INTERFACE
END MODULE S_DGTSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGTSVX
   INTERFACE
      SUBROUTINE DGTSVX(Fact,Trans,N,Nrhs,Dl,D,Du,Dlf,Df,Duf,Du2,Ipiv,B,&
     &                  Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(*) :: Dl
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: Du
      REAL(R8KIND) , DIMENSION(*) :: Dlf
      REAL(R8KIND) , DIMENSION(*) :: Df
      REAL(R8KIND) , DIMENSION(*) :: Duf
      REAL(R8KIND) , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGTSVX
   END INTERFACE
END MODULE S_DGTSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGTTRF
   INTERFACE
      SUBROUTINE DGTTRF(N,Dl,D,Du,Du2,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Du
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Du2
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER :: Info
      END SUBROUTINE DGTTRF
   END INTERFACE
END MODULE S_DGTTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGTTRS
   INTERFACE
      SUBROUTINE DGTTRS(Trans,N,Nrhs,Dl,D,Du,Du2,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: Dl
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: Du
      REAL(R8KIND) , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DGTTRS
   END INTERFACE
END MODULE S_DGTTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DGTTS2
   INTERFACE
      SUBROUTINE DGTTS2(Itrans,N,Nrhs,Dl,D,Du,Du2,Ipiv,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: Itrans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Dl
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Du
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Du2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE DGTTS2
   END INTERFACE
END MODULE S_DGTTS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DHGEQZ
   INTERFACE
      SUBROUTINE DHGEQZ(Job,Compq,Compz,N,Ilo,Ihi,H,Ldh,T,Ldt,Alphar,   &
     &                  Alphai,Beta,Q,Ldq,Z,Ldz,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  HALF = 0.5D+0 , ZERO = 0.0D+0 ,     &
     &                              ONE = 1.0D+0 , SAFETY = 1.0D+2
      CHARACTER :: Job
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DHGEQZ
   END INTERFACE
END MODULE S_DHGEQZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DHSEIN
   INTERFACE
      SUBROUTINE DHSEIN(Side,Eigsrc,Initv,Select,N,H,Ldh,Wr,Wi,Vl,Ldvl, &
     &                  Vr,Ldvr,Mm,M,Work,Ifaill,Ifailr,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Side
      CHARACTER :: Eigsrc
      CHARACTER :: Initv
      LOGICAL , INTENT(INOUT) , DIMENSION(*) :: Select
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wr
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Wi
      REAL(R8KIND) , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      REAL(R8KIND) , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifaill
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifailr
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DHSEIN
   END INTERFACE
END MODULE S_DHSEIN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DHSEQR
   INTERFACE
      SUBROUTINE DHSEQR(Job,Compz,N,Ilo,Ihi,H,Ldh,Wr,Wi,Z,Ldz,Work,     &
     &                  Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NTINY = 15 , NL = 49
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Job
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL(R8KIND) , DIMENSION(*) :: Wr
      REAL(R8KIND) , DIMENSION(*) :: Wi
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DHSEQR
   END INTERFACE
END MODULE S_DHSEQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DISNAN
   INTERFACE
      FUNCTION DISNAN(Din)
      USE F77KINDS                        
      IMPLICIT NONE
      LOGICAL :: DISNAN
      REAL(R8KIND) , INTENT(IN) :: Din
      END FUNCTION DISNAN
   END INTERFACE
END MODULE S_DISNAN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLABAD
   INTERFACE
      SUBROUTINE DLABAD(Small,Large)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) , INTENT(INOUT) :: Small
      REAL(R8KIND) , INTENT(INOUT) :: Large
      END SUBROUTINE DLABAD
   END INTERFACE
END MODULE S_DLABAD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLABRD
   INTERFACE
      SUBROUTINE DLABRD(M,N,Nb,A,Lda,D,E,Tauq,Taup,X,Ldx,Y,Ldy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Tauq
      REAL(R8KIND) , DIMENSION(*) :: Taup
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , DIMENSION(Ldy,*) :: Y
      INTEGER :: Ldy
      END SUBROUTINE DLABRD
   END INTERFACE
END MODULE S_DLABRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLACN2
   INTERFACE
      SUBROUTINE DLACN2(N,V,X,Isgn,Est,Kase,Isave)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: V
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isgn
      REAL(R8KIND) , INTENT(INOUT) :: Est
      INTEGER , INTENT(INOUT) :: Kase
      INTEGER , INTENT(INOUT) , DIMENSION(3) :: Isave
      END SUBROUTINE DLACN2
   END INTERFACE
END MODULE S_DLACN2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLACON
   INTERFACE
      SUBROUTINE DLACON(N,V,X,Isgn,Est,Kase)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: V
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isgn
      REAL(R8KIND) , INTENT(INOUT) :: Est
      INTEGER , INTENT(INOUT) :: Kase
      END SUBROUTINE DLACON
   END INTERFACE
END MODULE S_DLACON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLACPY
   INTERFACE
      SUBROUTINE DLACPY(Uplo,M,N,A,Lda,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE DLACPY
   END INTERFACE
END MODULE S_DLACPY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLADIV
   INTERFACE
      SUBROUTINE DLADIV(A,B,C,D,P,Q)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  BS = 2.0D0 , HALF = 0.5D0 ,         &
     &                              TWO = 2.0D0
      REAL(R8KIND) , INTENT(IN) :: A
      REAL(R8KIND) , INTENT(IN) :: B
      REAL(R8KIND) , INTENT(IN) :: C
      REAL(R8KIND) , INTENT(IN) :: D
      REAL(R8KIND) , INTENT(INOUT) :: P
      REAL(R8KIND) , INTENT(INOUT) :: Q
      END SUBROUTINE DLADIV
   END INTERFACE
END MODULE S_DLADIV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLADIV1
   INTERFACE
      SUBROUTINE DLADIV1(A,B,C,D,P,Q)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0
      REAL(R8KIND) , INTENT(INOUT) :: A
      REAL(R8KIND) :: B
      REAL(R8KIND) :: C
      REAL(R8KIND) :: D
      REAL(R8KIND) , INTENT(OUT) :: P
      REAL(R8KIND) , INTENT(OUT) :: Q
      END SUBROUTINE DLADIV1
   END INTERFACE
END MODULE S_DLADIV1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLADIV2
   INTERFACE
      FUNCTION DLADIV2(A,B,C,D,R,T)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0
      REAL(R8KIND) :: DLADIV2
      REAL(R8KIND) , INTENT(IN) :: A
      REAL(R8KIND) , INTENT(IN) :: B
      REAL(R8KIND) , INTENT(IN) :: C
      REAL(R8KIND) , INTENT(IN) :: D
      REAL(R8KIND) , INTENT(IN) :: R
      REAL(R8KIND) , INTENT(IN) :: T
      END FUNCTION DLADIV2
   END INTERFACE
END MODULE S_DLADIV2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAE2
   INTERFACE
      SUBROUTINE DLAE2(A,B,C,Rt1,Rt2)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , TWO = 2.0D0 ,         &
     &                              ZERO = 0.0D0 , HALF = 0.5D0
      REAL(R8KIND) , INTENT(IN) :: A
      REAL(R8KIND) , INTENT(IN) :: B
      REAL(R8KIND) , INTENT(IN) :: C
      REAL(R8KIND) , INTENT(INOUT) :: Rt1
      REAL(R8KIND) , INTENT(OUT) :: Rt2
      END SUBROUTINE DLAE2
   END INTERFACE
END MODULE S_DLAE2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAEBZ
   INTERFACE
      SUBROUTINE DLAEBZ(Ijob,Nitmax,N,Mmax,Minp,Nbmin,Abstol,Reltol,    &
     &                  Pivmin,D,E,E2,Nval,Ab,C,Mout,Nab,Work,Iwork,    &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , TWO = 2.0D0 ,        &
     &                              HALF = 1.0D0/TWO
      INTEGER , INTENT(IN) :: Ijob
      INTEGER , INTENT(IN) :: Nitmax
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Mmax
      INTEGER , INTENT(IN) :: Minp
      INTEGER , INTENT(IN) :: Nbmin
      REAL(R8KIND) , INTENT(IN) :: Abstol
      REAL(R8KIND) , INTENT(IN) :: Reltol
      REAL(R8KIND) , INTENT(IN) :: Pivmin
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: E2
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Nval
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Mmax,*) :: Ab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      INTEGER , INTENT(INOUT) :: Mout
      INTEGER , INTENT(INOUT) , DIMENSION(Mmax,*) :: Nab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAEBZ
   END INTERFACE
END MODULE S_DLAEBZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAED0
   INTERFACE
      SUBROUTINE DLAED0(Icompq,Qsiz,N,D,E,Q,Ldq,Qstore,Ldqs,Work,Iwork, &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.D0 , ONE = 1.D0 ,          &
     &                              TWO = 2.D0
      INTEGER :: Icompq
      INTEGER :: Qsiz
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , DIMENSION(Ldqs,*) :: Qstore
      INTEGER :: Ldqs
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAED0
   END INTERFACE
END MODULE S_DLAED0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAED1
   INTERFACE
      SUBROUTINE DLAED1(N,D,Q,Ldq,Indxq,Rho,Cutpnt,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , DIMENSION(*) :: Indxq
      REAL(R8KIND) :: Rho
      INTEGER :: Cutpnt
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAED1
   END INTERFACE
END MODULE S_DLAED1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAED2
   INTERFACE
      SUBROUTINE DLAED2(K,N,N1,D,Q,Ldq,Indxq,Rho,Z,Dlamda,W,Q2,Indx,    &
     &                  Indxc,Indxp,Coltyp,Info)
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
      INTEGER :: N1
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxq
      REAL(R8KIND) , INTENT(INOUT) :: Rho
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL(R8KIND) , DIMENSION(*) :: Dlamda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(*) :: Q2
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxc
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxp
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Coltyp
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAED2
   END INTERFACE
END MODULE S_DLAED2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAED3
   INTERFACE
      SUBROUTINE DLAED3(K,N,N1,D,Q,Ldq,Rho,Dlamda,Q2,Indx,Ctot,W,S,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0
      INTEGER :: K
      INTEGER , INTENT(IN) :: N
      INTEGER :: N1
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) :: Rho
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dlamda
      REAL(R8KIND) , DIMENSION(*) :: Q2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Indx
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ctot
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAED3
   END INTERFACE
END MODULE S_DLAED3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAED4
   INTERFACE
      SUBROUTINE DLAED4(N,I,D,Z,Delta,Rho,Dlam,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXIT = 30
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , THREE = 3.0D0 ,       &
     &                              FOUR = 4.0D0 , EIGHT = 8.0D0 ,      &
     &                              TEN = 10.0D0
      INTEGER , INTENT(IN) :: N
      INTEGER :: I
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: Z
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Delta
      REAL(R8KIND) :: Rho
      REAL(R8KIND) :: Dlam
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAED4
   END INTERFACE
END MODULE S_DLAED4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAED5
   INTERFACE
      SUBROUTINE DLAED5(I,D,Z,Delta,Rho,Dlam)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , FOUR = 4.0D0
      INTEGER , INTENT(IN) :: I
      REAL(R8KIND) , INTENT(IN) , DIMENSION(2) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(2) :: Z
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(2) :: Delta
      REAL(R8KIND) , INTENT(IN) :: Rho
      REAL(R8KIND) , INTENT(OUT) :: Dlam
      END SUBROUTINE DLAED5
   END INTERFACE
END MODULE S_DLAED5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAED6
   INTERFACE
      SUBROUTINE DLAED6(Kniter,Orgati,Rho,D,Z,Finit,Tau,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXIT = 40
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , THREE = 3.0D0 ,       &
     &                              FOUR = 4.0D0 , EIGHT = 8.0D0
      INTEGER , INTENT(IN) :: Kniter
      LOGICAL , INTENT(IN) :: Orgati
      REAL(R8KIND) , INTENT(IN) :: Rho
      REAL(R8KIND) , INTENT(IN) , DIMENSION(3) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(3) :: Z
      REAL(R8KIND) , INTENT(IN) :: Finit
      REAL(R8KIND) , INTENT(INOUT) :: Tau
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLAED6
   END INTERFACE
END MODULE S_DLAED6
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAED7
   INTERFACE
      SUBROUTINE DLAED7(Icompq,N,Qsiz,Tlvls,Curlvl,Curpbm,D,Q,Ldq,Indxq,&
     &                  Rho,Cutpnt,Qstore,Qptr,Prmptr,Perm,Givptr,      &
     &                  Givcol,Givnum,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0
      INTEGER :: Icompq
      INTEGER :: N
      INTEGER :: Qsiz
      INTEGER :: Tlvls
      INTEGER :: Curlvl
      INTEGER :: Curpbm
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , DIMENSION(*) :: Indxq
      REAL(R8KIND) :: Rho
      INTEGER :: Cutpnt
      REAL(R8KIND) , DIMENSION(*) :: Qstore
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Qptr
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Prmptr
      INTEGER , DIMENSION(*) :: Perm
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Givptr
      INTEGER , DIMENSION(2,*) :: Givcol
      REAL(R8KIND) , DIMENSION(2,*) :: Givnum
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAED7
   END INTERFACE
END MODULE S_DLAED7
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAED8
   INTERFACE
      SUBROUTINE DLAED8(Icompq,K,N,Qsiz,D,Q,Ldq,Indxq,Rho,Cutpnt,Z,     &
     &                  Dlamda,Q2,Ldq2,W,Perm,Givptr,Givcol,Givnum,     &
     &                  Indxp,Indx,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  MONE = -1.0D0 , ZERO = 0.0D0 ,      &
     &                              ONE = 1.0D0 , TWO = 2.0D0 ,         &
     &                              EIGHT = 8.0D0
      INTEGER , INTENT(IN) :: Icompq
      INTEGER , INTENT(INOUT) :: K
      INTEGER :: N
      INTEGER :: Qsiz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxq
      REAL(R8KIND) , INTENT(INOUT) :: Rho
      INTEGER , INTENT(IN) :: Cutpnt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dlamda
      REAL(R8KIND) , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Perm
      INTEGER , INTENT(INOUT) :: Givptr
      INTEGER , INTENT(OUT) , DIMENSION(2,*) :: Givcol
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(2,*) :: Givnum
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indxp
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indx
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAED8
   END INTERFACE
END MODULE S_DLAED8
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAED9
   INTERFACE
      SUBROUTINE DLAED9(K,Kstart,Kstop,N,D,Q,Ldq,Rho,Dlamda,W,S,Lds,    &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: K
      INTEGER , INTENT(IN) :: Kstart
      INTEGER , INTENT(IN) :: Kstop
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER , INTENT(IN) :: Ldq
      REAL(R8KIND) :: Rho
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dlamda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lds,*) :: S
      INTEGER , INTENT(IN) :: Lds
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAED9
   END INTERFACE
END MODULE S_DLAED9
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAEDA
   INTERFACE
      SUBROUTINE DLAEDA(N,Tlvls,Curlvl,Curpbm,Prmptr,Perm,Givptr,Givcol,&
     &                  Givnum,Q,Qptr,Z,Ztemp,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , HALF = 0.5D0 ,       &
     &                              ONE = 1.0D0
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Tlvls
      INTEGER , INTENT(IN) :: Curlvl
      INTEGER , INTENT(IN) :: Curpbm
      INTEGER , INTENT(IN) , DIMENSION(*) :: Prmptr
      INTEGER , INTENT(IN) , DIMENSION(*) :: Perm
      INTEGER , INTENT(IN) , DIMENSION(*) :: Givptr
      INTEGER , INTENT(IN) , DIMENSION(2,*) :: Givcol
      REAL(R8KIND) , DIMENSION(2,*) :: Givnum
      REAL(R8KIND) , DIMENSION(*) :: Q
      INTEGER , INTENT(IN) , DIMENSION(*) :: Qptr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL(R8KIND) , DIMENSION(*) :: Ztemp
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAEDA
   END INTERFACE
END MODULE S_DLAEDA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAEIN
   INTERFACE
      SUBROUTINE DLAEIN(Rightv,Noinit,N,H,Ldh,Wr,Wi,Vr,Vi,B,Ldb,Work,   &
     &                  Eps3,Smlnum,Bignum,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TENTH = 1.0D-1
      LOGICAL , INTENT(IN) :: Rightv
      LOGICAL , INTENT(IN) :: Noinit
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldh,*) :: H
      INTEGER , INTENT(IN) :: Ldh
      REAL(R8KIND) , INTENT(IN) :: Wr
      REAL(R8KIND) , INTENT(IN) :: Wi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(IN) :: Eps3
      REAL(R8KIND) , INTENT(IN) :: Smlnum
      REAL(R8KIND) , INTENT(IN) :: Bignum
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLAEIN
   END INTERFACE
END MODULE S_DLAEIN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAEV2
   INTERFACE
      SUBROUTINE DLAEV2(A,B,C,Rt1,Rt2,Cs1,Sn1)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , TWO = 2.0D0 ,         &
     &                              ZERO = 0.0D0 , HALF = 0.5D0
      REAL(R8KIND) , INTENT(IN) :: A
      REAL(R8KIND) , INTENT(IN) :: B
      REAL(R8KIND) , INTENT(IN) :: C
      REAL(R8KIND) , INTENT(INOUT) :: Rt1
      REAL(R8KIND) , INTENT(OUT) :: Rt2
      REAL(R8KIND) , INTENT(INOUT) :: Cs1
      REAL(R8KIND) , INTENT(INOUT) :: Sn1
      END SUBROUTINE DLAEV2
   END INTERFACE
END MODULE S_DLAEV2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAEXC
   INTERFACE
      SUBROUTINE DLAEXC(Wantq,N,T,Ldt,Q,Ldq,J1,N1,N2,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TEN = 1.0D+1
      INTEGER , PARAMETER  ::  LDD = 4 , LDX = 2
      LOGICAL , INTENT(IN) :: Wantq
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , INTENT(IN) :: J1
      INTEGER :: N1
      INTEGER :: N2
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLAEXC
   END INTERFACE
END MODULE S_DLAEXC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAG2
   INTERFACE
      SUBROUTINE DLAG2(A,Lda,B,Ldb,Safmin,Scale1,Scale2,Wr1,Wr2,Wi)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0 , HALF = ONE/TWO ,     &
     &                              FUZZY1 = ONE + 1.0D-5
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Safmin
      REAL(R8KIND) , INTENT(INOUT) :: Scale1
      REAL(R8KIND) , INTENT(OUT) :: Scale2
      REAL(R8KIND) , INTENT(INOUT) :: Wr1
      REAL(R8KIND) , INTENT(INOUT) :: Wr2
      REAL(R8KIND) , INTENT(INOUT) :: Wi
      END SUBROUTINE DLAG2
   END INTERFACE
END MODULE S_DLAG2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAG2S
   INTERFACE
      SUBROUTINE DLAG2S(M,N,A,Lda,Sa,Ldsa,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(OUT) , DIMENSION(Ldsa,*) :: Sa
      INTEGER , INTENT(IN) :: Ldsa
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLAG2S
   END INTERFACE
END MODULE S_DLAG2S
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_GBAMV
   INTERFACE
      SUBROUTINE DLA_GBAMV(Trans,M,N,Kl,Ku,Alpha,Ab,Ldab,X,Incx,Beta,Y, &
     &                     Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE DLA_GBAMV
   END INTERFACE
END MODULE S_DLA_GBAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_GBRCOND
   INTERFACE
      FUNCTION DLA_GBRCOND(Trans,N,Kl,Ku,Ab,Ldab,Afb,Ldafb,Ipiv,Cmode,C,&
     &                     Info,Work,Iwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DLA_GBRCOND
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(IN) :: Cmode
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      INTEGER , INTENT(INOUT) :: Info
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      END FUNCTION DLA_GBRCOND
   END INTERFACE
END MODULE S_DLA_GBRCOND
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_GBRFSX_EXTENDED
   INTERFACE
      SUBROUTINE DLA_GBRFSX_EXTENDED(Prec_type,Trans_type,N,Kl,Ku,Nrhs, &
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
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL(R8KIND) , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      REAL(R8KIND) , DIMENSION(*) :: Res
      REAL(R8KIND) , DIMENSION(*) :: Ayb
      REAL(R8KIND) , DIMENSION(*) :: Dy
      REAL(R8KIND) , DIMENSION(*) :: Y_tail
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL(R8KIND) , INTENT(IN) :: Rthresh
      REAL(R8KIND) , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
      END SUBROUTINE DLA_GBRFSX_EXTENDED
   END INTERFACE
END MODULE S_DLA_GBRFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_GBRPVGRW
   INTERFACE
      FUNCTION DLA_GBRPVGRW(N,Kl,Ku,Ncols,Ab,Ldab,Afb,Ldafb)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DLA_GBRPVGRW
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      INTEGER , INTENT(IN) :: Ncols
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldafb,*) :: Afb
      INTEGER , INTENT(IN) :: Ldafb
      END FUNCTION DLA_GBRPVGRW
   END INTERFACE
END MODULE S_DLA_GBRPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_GEAMV
   INTERFACE
      SUBROUTINE DLA_GEAMV(Trans,M,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: Trans
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
      END SUBROUTINE DLA_GEAMV
   END INTERFACE
END MODULE S_DLA_GEAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_GERCOND
   INTERFACE
      FUNCTION DLA_GERCOND(Trans,N,A,Lda,Af,Ldaf,Ipiv,Cmode,C,Info,Work,&
     &                     Iwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DLA_GERCOND
      CHARACTER :: Trans
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(IN) :: Cmode
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      INTEGER , INTENT(INOUT) :: Info
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      END FUNCTION DLA_GERCOND
   END INTERFACE
END MODULE S_DLA_GERCOND
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_GERFSX_EXTENDED
   INTERFACE
      SUBROUTINE DLA_GERFSX_EXTENDED(Prec_type,Trans_type,N,Nrhs,A,Lda, &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL(R8KIND) , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Errs_n
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Errs_c
      REAL(R8KIND) , DIMENSION(*) :: Res
      REAL(R8KIND) , DIMENSION(*) :: Ayb
      REAL(R8KIND) , DIMENSION(*) :: Dy
      REAL(R8KIND) , DIMENSION(*) :: Y_tail
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL(R8KIND) , INTENT(IN) :: Rthresh
      REAL(R8KIND) , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
      END SUBROUTINE DLA_GERFSX_EXTENDED
   END INTERFACE
END MODULE S_DLA_GERFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_GERPVGRW
   INTERFACE
      FUNCTION DLA_GERPVGRW(N,Ncols,A,Lda,Af,Ldaf)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DLA_GERPVGRW
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ncols
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      END FUNCTION DLA_GERPVGRW
   END INTERFACE
END MODULE S_DLA_GERPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAGS2
   INTERFACE
      SUBROUTINE DLAGS2(Upper,A1,A2,A3,B1,B2,B3,Csu,Snu,Csv,Snv,Csq,Snq)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      LOGICAL , INTENT(IN) :: Upper
      REAL(R8KIND) , INTENT(IN) :: A1
      REAL(R8KIND) , INTENT(IN) :: A2
      REAL(R8KIND) , INTENT(IN) :: A3
      REAL(R8KIND) , INTENT(IN) :: B1
      REAL(R8KIND) , INTENT(IN) :: B2
      REAL(R8KIND) , INTENT(IN) :: B3
      REAL(R8KIND) , INTENT(OUT) :: Csu
      REAL(R8KIND) , INTENT(OUT) :: Snu
      REAL(R8KIND) , INTENT(OUT) :: Csv
      REAL(R8KIND) , INTENT(OUT) :: Snv
      REAL(R8KIND) :: Csq
      REAL(R8KIND) :: Snq
      END SUBROUTINE DLAGS2
   END INTERFACE
END MODULE S_DLAGS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAGTF
   INTERFACE
      SUBROUTINE DLAGTF(N,A,Lambda,B,C,Tol,D,In,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: A
      REAL(R8KIND) , INTENT(IN) :: Lambda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(IN) :: Tol
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: In
      INTEGER :: Info
      END SUBROUTINE DLAGTF
   END INTERFACE
END MODULE S_DLAGTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAGTM
   INTERFACE
      SUBROUTINE DLAGTM(Trans,N,Nrhs,Alpha,Dl,D,Du,X,Ldx,Beta,B,Ldb)
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
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Dl
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Du
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      END SUBROUTINE DLAGTM
   END INTERFACE
END MODULE S_DLAGTM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAGTS
   INTERFACE
      SUBROUTINE DLAGTS(Job,N,A,B,C,D,In,Y,Tol,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: Job
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: A
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: B
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      INTEGER , INTENT(IN) , DIMENSION(*) :: In
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      REAL(R8KIND) , INTENT(INOUT) :: Tol
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAGTS
   END INTERFACE
END MODULE S_DLAGTS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAGV2
   INTERFACE
      SUBROUTINE DLAGV2(A,Lda,B,Ldb,Alphar,Alphai,Beta,Csl,Snl,Csr,Snr)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(2) :: Alphar
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(2) :: Alphai
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(2) :: Beta
      REAL(R8KIND) :: Csl
      REAL(R8KIND) :: Snl
      REAL(R8KIND) :: Csr
      REAL(R8KIND) , INTENT(INOUT) :: Snr
      END SUBROUTINE DLAGV2
   END INTERFACE
END MODULE S_DLAGV2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAHQR
   INTERFACE
      SUBROUTINE DLAHQR(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,Wr,Wi,Iloz,Ihiz,Z,  &
     &                  Ldz,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , DAT1 = 3.0D0/4.0D0 ,  &
     &                              DAT2 = -0.4375D0
      INTEGER , PARAMETER  ::  KEXSH = 10
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL(R8KIND) , DIMENSION(*) :: Wr
      REAL(R8KIND) , DIMENSION(*) :: Wi
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLAHQR
   END INTERFACE
END MODULE S_DLAHQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAHR2
   INTERFACE
      SUBROUTINE DLAHR2(N,K,Nb,A,Lda,Tau,T,Ldt,Y,Ldy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , INTENT(IN) :: N
      INTEGER :: K
      INTEGER :: Nb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Nb) :: Tau
      REAL(R8KIND) , DIMENSION(Ldt,Nb) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Ldy,Nb) :: Y
      INTEGER :: Ldy
      END SUBROUTINE DLAHR2
   END INTERFACE
END MODULE S_DLAHR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAIC1
   INTERFACE
      SUBROUTINE DLAIC1(Job,J,X,Sest,W,Gamma,Sestpr,S,C)
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
      REAL(R8KIND) , DIMENSION(J) :: X
      REAL(R8KIND) , INTENT(IN) :: Sest
      REAL(R8KIND) , DIMENSION(J) :: W
      REAL(R8KIND) , INTENT(IN) :: Gamma
      REAL(R8KIND) , INTENT(OUT) :: Sestpr
      REAL(R8KIND) , INTENT(INOUT) :: S
      REAL(R8KIND) , INTENT(INOUT) :: C
      END SUBROUTINE DLAIC1
   END INTERFACE
END MODULE S_DLAIC1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAISNAN
   INTERFACE
      FUNCTION DLAISNAN(Din1,Din2)
      USE F77KINDS                        
      IMPLICIT NONE
      LOGICAL :: DLAISNAN
      REAL(R8KIND) , INTENT(IN) :: Din1
      REAL(R8KIND) , INTENT(IN) :: Din2
      END FUNCTION DLAISNAN
   END INTERFACE
END MODULE S_DLAISNAN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_LIN_BERR
   INTERFACE
      SUBROUTINE DLA_LIN_BERR(N,Nz,Nrhs,Res,Ayb,Berr)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nz
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , INTENT(IN) , DIMENSION(N,Nrhs) :: Res
      REAL(R8KIND) , INTENT(IN) , DIMENSION(N,Nrhs) :: Ayb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs) :: Berr
      END SUBROUTINE DLA_LIN_BERR
   END INTERFACE
END MODULE S_DLA_LIN_BERR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLALN2
   INTERFACE
      SUBROUTINE DLALN2(Ltrans,Na,Nw,Smin,Ca,A,Lda,D1,D2,B,Ldb,Wr,Wi,X, &
     &                  Ldx,Scale,Xnorm,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0
      LOGICAL , INTENT(IN) :: Ltrans
      INTEGER , INTENT(IN) :: Na
      INTEGER , INTENT(IN) :: Nw
      REAL(R8KIND) , INTENT(IN) :: Smin
      REAL(R8KIND) , INTENT(IN) :: Ca
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) :: D1
      REAL(R8KIND) , INTENT(IN) :: D2
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Wr
      REAL(R8KIND) , INTENT(IN) :: Wi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) :: Xnorm
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLALN2
   END INTERFACE
END MODULE S_DLALN2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLALS0
   INTERFACE
      SUBROUTINE DLALS0(Icompq,Nl,Nr,Sqre,Nrhs,B,Ldb,Bx,Ldbx,Perm,      &
     &                  Givptr,Givcol,Ldgcol,Givnum,Ldgnum,Poles,Difl,  &
     &                  Difr,Z,K,C,S,Work,Info)
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
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldbx,*) :: Bx
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
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLALS0
   END INTERFACE
END MODULE S_DLALS0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLALSA
   INTERFACE
      SUBROUTINE DLALSA(Icompq,Smlsiz,N,Nrhs,B,Ldb,Bx,Ldbx,U,Ldu,Vt,K,  &
     &                  Difl,Difr,Z,Poles,Givptr,Givcol,Ldgcol,Perm,    &
     &                  Givnum,C,S,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldbx,*) :: Bx
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
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLALSA
   END INTERFACE
END MODULE S_DLALSA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLALSD
   INTERFACE
      SUBROUTINE DLALSD(Uplo,Smlsiz,N,Nrhs,D,E,B,Ldb,Rcond,Rank,Work,   &
     &                  Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0
      CHARACTER , INTENT(IN) :: Uplo
      INTEGER :: Smlsiz
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(INOUT) :: Rank
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLALSD
   END INTERFACE
END MODULE S_DLALSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAMRG
   INTERFACE
      SUBROUTINE DLAMRG(N1,N2,A,Dtrd1,Dtrd2,Index)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N1
      INTEGER , INTENT(IN) :: N2
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: A
      INTEGER , INTENT(IN) :: Dtrd1
      INTEGER , INTENT(IN) :: Dtrd2
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Index
      END SUBROUTINE DLAMRG
   END INTERFACE
END MODULE S_DLAMRG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAMSWLQ
   INTERFACE
      SUBROUTINE DLAMSWLQ(Side,Trans,M,N,K,Mb,Nb,A,Lda,T,Ldt,C,Ldc,Work,&
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAMSWLQ
   END INTERFACE
END MODULE S_DLAMSWLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAMTSQR
   INTERFACE
      SUBROUTINE DLAMTSQR(Side,Trans,M,N,K,Mb,Nb,A,Lda,T,Ldt,C,Ldc,Work,&
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAMTSQR
   END INTERFACE
END MODULE S_DLAMTSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANEG
   INTERFACE
      FUNCTION DLANEG(N,D,Lld,Sigma,Pivmin,R)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      INTEGER , PARAMETER  ::  BLKLEN = 128
      INTEGER :: DLANEG
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Lld
      REAL(R8KIND) , INTENT(IN) :: Sigma
      REAL(R8KIND) :: Pivmin
      INTEGER , INTENT(IN) :: R
      END FUNCTION DLANEG
   END INTERFACE
END MODULE S_DLANEG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANGB
   INTERFACE
      FUNCTION DLANGB(Norm,N,Kl,Ku,Ab,Ldab,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANGB
      CHARACTER :: Norm
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION DLANGB
   END INTERFACE
END MODULE S_DLANGB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANGE
   INTERFACE
      FUNCTION DLANGE(Norm,M,N,A,Lda,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANGE
      CHARACTER :: Norm
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION DLANGE
   END INTERFACE
END MODULE S_DLANGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANGT
   INTERFACE
      FUNCTION DLANGT(Norm,N,Dl,D,Du)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANGT
      CHARACTER :: Norm
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: Dl
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: Du
      END FUNCTION DLANGT
   END INTERFACE
END MODULE S_DLANGT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANHS
   INTERFACE
      FUNCTION DLANHS(Norm,N,A,Lda,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANHS
      CHARACTER :: Norm
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION DLANHS
   END INTERFACE
END MODULE S_DLANHS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANSB
   INTERFACE
      FUNCTION DLANSB(Norm,Uplo,N,K,Ab,Ldab,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANSB
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION DLANSB
   END INTERFACE
END MODULE S_DLANSB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANSF
   INTERFACE
      FUNCTION DLANSF(Norm,Transr,Uplo,N,A,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANSF
      CHARACTER :: Norm
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(0:*) :: A
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(0:*) :: Work
      END FUNCTION DLANSF
   END INTERFACE
END MODULE S_DLANSF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANSP
   INTERFACE
      FUNCTION DLANSP(Norm,Uplo,N,Ap,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANSP
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION DLANSP
   END INTERFACE
END MODULE S_DLANSP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANST
   INTERFACE
      FUNCTION DLANST(Norm,N,D,E)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANST
      CHARACTER :: Norm
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      END FUNCTION DLANST
   END INTERFACE
END MODULE S_DLANST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANSY
   INTERFACE
      FUNCTION DLANSY(Norm,Uplo,N,A,Lda,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANSY
      CHARACTER :: Norm
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION DLANSY
   END INTERFACE
END MODULE S_DLANSY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANTB
   INTERFACE
      FUNCTION DLANTB(Norm,Uplo,Diag,N,K,Ab,Ldab,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANTB
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION DLANTB
   END INTERFACE
END MODULE S_DLANTB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANTP
   INTERFACE
      FUNCTION DLANTP(Norm,Uplo,Diag,N,Ap,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANTP
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION DLANTP
   END INTERFACE
END MODULE S_DLANTP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANTR
   INTERFACE
      FUNCTION DLANTR(Norm,Uplo,Diag,M,N,A,Lda,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DLANTR
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION DLANTR
   END INTERFACE
END MODULE S_DLANTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLANV2
   INTERFACE
      SUBROUTINE DLANV2(A,B,C,D,Rt1r,Rt1i,Rt2r,Rt2i,Cs,Sn)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , HALF = 0.5D+0 ,     &
     &                              ONE = 1.0D+0 , TWO = 2.0D0 ,        &
     &                              MULTPL = 4.0D+0
      REAL(R8KIND) , INTENT(INOUT) :: A
      REAL(R8KIND) , INTENT(INOUT) :: B
      REAL(R8KIND) , INTENT(INOUT) :: C
      REAL(R8KIND) , INTENT(INOUT) :: D
      REAL(R8KIND) , INTENT(OUT) :: Rt1r
      REAL(R8KIND) , INTENT(INOUT) :: Rt1i
      REAL(R8KIND) , INTENT(OUT) :: Rt2r
      REAL(R8KIND) , INTENT(OUT) :: Rt2i
      REAL(R8KIND) , INTENT(INOUT) :: Cs
      REAL(R8KIND) , INTENT(INOUT) :: Sn
      END SUBROUTINE DLANV2
   END INTERFACE
END MODULE S_DLANV2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAORHR_COL_GETRFNP2
   INTERFACE
      RECURSIVE SUBROUTINE DLAORHR_COL_GETRFNP2(M,N,A,Lda,D,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAORHR_COL_GETRFNP2
   END INTERFACE
END MODULE S_DLAORHR_COL_GETRFNP2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAORHR_COL_GETRFNP
   INTERFACE
      SUBROUTINE DLAORHR_COL_GETRFNP(M,N,A,Lda,D,Info)
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
      REAL(R8KIND) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAORHR_COL_GETRFNP
   END INTERFACE
END MODULE S_DLAORHR_COL_GETRFNP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAPLL
   INTERFACE
      SUBROUTINE DLAPLL(N,X,Incx,Y,Incy,Ssmin)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER :: Incx
      REAL(R8KIND) , DIMENSION(*) :: Y
      INTEGER :: Incy
      REAL(R8KIND) :: Ssmin
      END SUBROUTINE DLAPLL
   END INTERFACE
END MODULE S_DLAPLL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAPMR
   INTERFACE
      SUBROUTINE DLAPMR(Forwrd,M,N,X,Ldx,K)
      USE F77KINDS                        
      IMPLICIT NONE
      LOGICAL , INTENT(IN) :: Forwrd
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: K
      END SUBROUTINE DLAPMR
   END INTERFACE
END MODULE S_DLAPMR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAPMT
   INTERFACE
      SUBROUTINE DLAPMT(Forwrd,M,N,X,Ldx,K)
      USE F77KINDS                        
      IMPLICIT NONE
      LOGICAL , INTENT(IN) :: Forwrd
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: K
      END SUBROUTINE DLAPMT
   END INTERFACE
END MODULE S_DLAPMT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_PORCOND
   INTERFACE
      FUNCTION DLA_PORCOND(Uplo,N,A,Lda,Af,Ldaf,Cmode,C,Info,Work,Iwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DLA_PORCOND
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , INTENT(IN) :: Cmode
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      INTEGER , INTENT(INOUT) :: Info
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      END FUNCTION DLA_PORCOND
   END INTERFACE
END MODULE S_DLA_PORCOND
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_PORFSX_EXTENDED
   INTERFACE
      SUBROUTINE DLA_PORFSX_EXTENDED(Prec_type,Uplo,N,Nrhs,A,Lda,Af,    &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      LOGICAL , INTENT(IN) :: Colequ
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL(R8KIND) , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      REAL(R8KIND) , DIMENSION(*) :: Res
      REAL(R8KIND) , DIMENSION(*) :: Ayb
      REAL(R8KIND) , DIMENSION(*) :: Dy
      REAL(R8KIND) , DIMENSION(*) :: Y_tail
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL(R8KIND) , INTENT(IN) :: Rthresh
      REAL(R8KIND) , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
      END SUBROUTINE DLA_PORFSX_EXTENDED
   END INTERFACE
END MODULE S_DLA_PORFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_PORPVGRW
   INTERFACE
      FUNCTION DLA_PORPVGRW(Uplo,Ncols,A,Lda,Af,Ldaf,Work)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DLA_PORPVGRW
      CHARACTER(1) :: Uplo
      INTEGER , INTENT(IN) :: Ncols
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION DLA_PORPVGRW
   END INTERFACE
END MODULE S_DLA_PORPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAPY2
   INTERFACE
      FUNCTION DLAPY2(X,Y)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      REAL(R8KIND) :: DLAPY2
      REAL(R8KIND) :: X
      REAL(R8KIND) :: Y
      END FUNCTION DLAPY2
   END INTERFACE
END MODULE S_DLAPY2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAPY3
   INTERFACE
      FUNCTION DLAPY3(X,Y,Z)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0
      REAL(R8KIND) :: DLAPY3
      REAL(R8KIND) , INTENT(IN) :: X
      REAL(R8KIND) , INTENT(IN) :: Y
      REAL(R8KIND) , INTENT(IN) :: Z
      END FUNCTION DLAPY3
   END INTERFACE
END MODULE S_DLAPY3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQGB
   INTERFACE
      SUBROUTINE DLAQGB(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Equed)
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
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(IN) :: Rowcnd
      REAL(R8KIND) , INTENT(IN) :: Colcnd
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE DLAQGB
   END INTERFACE
END MODULE S_DLAQGB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQGE
   INTERFACE
      SUBROUTINE DLAQGE(M,N,A,Lda,R,C,Rowcnd,Colcnd,Amax,Equed)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(IN) :: Rowcnd
      REAL(R8KIND) , INTENT(IN) :: Colcnd
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE DLAQGE
   END INTERFACE
END MODULE S_DLAQGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQP2
   INTERFACE
      SUBROUTINE DLAQP2(M,N,Offset,A,Lda,Jpvt,Tau,Vn1,Vn2,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Offset
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vn1
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vn2
      REAL(R8KIND) , DIMENSION(*) :: Work
      END SUBROUTINE DLAQP2
   END INTERFACE
END MODULE S_DLAQP2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQPS
   INTERFACE
      SUBROUTINE DLAQPS(M,N,Offset,Nb,Kb,A,Lda,Jpvt,Tau,Vn1,Vn2,Auxv,F, &
     &                  Ldf)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Offset
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Kb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vn1
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vn2
      REAL(R8KIND) , DIMENSION(*) :: Auxv
      REAL(R8KIND) , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      END SUBROUTINE DLAQPS
   END INTERFACE
END MODULE S_DLAQPS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQR0
   INTERFACE
      SUBROUTINE DLAQR0(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,Wr,Wi,Iloz,Ihiz,Z,  &
     &                  Ldz,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NTINY = 15 , KEXNW = 5 , KEXSH = 6
      REAL(R8KIND) , PARAMETER  ::  WILK1 = 0.75D0 , WILK2 = -0.4375D0 ,&
     &                              ZERO = 0.0D0 , ONE = 1.0D0
      LOGICAL :: Wantt
      LOGICAL :: Wantz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wi
      INTEGER :: Iloz
      INTEGER :: Ihiz
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER :: Info
      END SUBROUTINE DLAQR0
   END INTERFACE
END MODULE S_DLAQR0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQR1
   INTERFACE
      SUBROUTINE DLAQR1(N,H,Ldh,Sr1,Si1,Sr2,Si2,V)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldh,*) :: H
      INTEGER , INTENT(IN) :: Ldh
      REAL(R8KIND) , INTENT(IN) :: Sr1
      REAL(R8KIND) , INTENT(IN) :: Si1
      REAL(R8KIND) , INTENT(IN) :: Sr2
      REAL(R8KIND) , INTENT(IN) :: Si2
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: V
      END SUBROUTINE DLAQR1
   END INTERFACE
END MODULE S_DLAQR1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQR2
   INTERFACE
      SUBROUTINE DLAQR2(Wantt,Wantz,N,Ktop,Kbot,Nw,H,Ldh,Iloz,Ihiz,Z,   &
     &                  Ldz,Ns,Nd,Sr,Si,V,Ldv,Nh,T,Ldt,Nv,Wv,Ldwv,Work, &
     &                  Lwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nw
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ns
      INTEGER , INTENT(OUT) :: Nd
      REAL(R8KIND) , DIMENSION(*) :: Sr
      REAL(R8KIND) , DIMENSION(*) :: Si
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(IN) :: Nh
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(IN) :: Nv
      REAL(R8KIND) , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      END SUBROUTINE DLAQR2
   END INTERFACE
END MODULE S_DLAQR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQR3
   INTERFACE
      SUBROUTINE DLAQR3(Wantt,Wantz,N,Ktop,Kbot,Nw,H,Ldh,Iloz,Ihiz,Z,   &
     &                  Ldz,Ns,Nd,Sr,Si,V,Ldv,Nh,T,Ldt,Nv,Wv,Ldwv,Work, &
     &                  Lwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nw
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ns
      INTEGER , INTENT(OUT) :: Nd
      REAL(R8KIND) , DIMENSION(*) :: Sr
      REAL(R8KIND) , DIMENSION(*) :: Si
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(IN) :: Nh
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(IN) :: Nv
      REAL(R8KIND) , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      END SUBROUTINE DLAQR3
   END INTERFACE
END MODULE S_DLAQR3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQR4
   INTERFACE
      SUBROUTINE DLAQR4(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,Wr,Wi,Iloz,Ihiz,Z,  &
     &                  Ldz,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  NTINY = 15 , KEXNW = 5 , KEXSH = 6
      REAL(R8KIND) , PARAMETER  ::  WILK1 = 0.75D0 , WILK2 = -0.4375D0 ,&
     &                              ZERO = 0.0D0 , ONE = 1.0D0
      LOGICAL :: Wantt
      LOGICAL :: Wantz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wi
      INTEGER :: Iloz
      INTEGER :: Ihiz
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER :: Info
      END SUBROUTINE DLAQR4
   END INTERFACE
END MODULE S_DLAQR4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQR5
   INTERFACE
      SUBROUTINE DLAQR5(Wantt,Wantz,Kacc22,N,Ktop,Kbot,Nshfts,Sr,Si,H,  &
     &                  Ldh,Iloz,Ihiz,Z,Ldz,V,Ldv,U,Ldu,Nv,Wv,Ldwv,Nh,  &
     &                  Wh,Ldwh)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: Kacc22
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nshfts
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Sr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Si
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER , INTENT(IN) :: Ldv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      INTEGER , INTENT(IN) :: Nv
      REAL(R8KIND) , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      INTEGER , INTENT(IN) :: Nh
      REAL(R8KIND) , DIMENSION(Ldwh,*) :: Wh
      INTEGER :: Ldwh
      END SUBROUTINE DLAQR5
   END INTERFACE
END MODULE S_DLAQR5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQSB
   INTERFACE
      SUBROUTINE DLAQSB(Uplo,N,Kd,Ab,Ldab,S,Scond,Amax,Equed)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Scond
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE DLAQSB
   END INTERFACE
END MODULE S_DLAQSB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQSP
   INTERFACE
      SUBROUTINE DLAQSP(Uplo,N,Ap,S,Scond,Amax,Equed)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Scond
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE DLAQSP
   END INTERFACE
END MODULE S_DLAQSP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQSY
   INTERFACE
      SUBROUTINE DLAQSY(Uplo,N,A,Lda,S,Scond,Amax,Equed)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(IN) :: Scond
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
      END SUBROUTINE DLAQSY
   END INTERFACE
END MODULE S_DLAQSY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAQTR
   INTERFACE
      SUBROUTINE DLAQTR(Ltran,Lreal,N,T,Ldt,B,W,Scale,X,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      LOGICAL , INTENT(IN) :: Ltran
      LOGICAL , INTENT(IN) :: Lreal
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(*) :: B
      REAL(R8KIND) :: W
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLAQTR
   END INTERFACE
END MODULE S_DLAQTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAR1V
   INTERFACE
      SUBROUTINE DLAR1V(N,B1,Bn,Lambda,D,L,Ld,Lld,Pivmin,Gaptol,Z,      &
     &                  Wantnc,Negcnt,Ztz,Mingma,R,Isuppz,Nrminv,Resid, &
     &                  Rqcorr,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
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
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
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
      END SUBROUTINE DLAR1V
   END INTERFACE
END MODULE S_DLAR1V
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAR2V
   INTERFACE
      SUBROUTINE DLAR2V(N,X,Y,Z,Incx,C,S,Incc)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      INTEGER , INTENT(IN) :: Incc
      END SUBROUTINE DLAR2V
   END INTERFACE
END MODULE S_DLAR2V
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARFB
   INTERFACE
      SUBROUTINE DLARFB(Side,Trans,Direct,Storev,M,N,K,V,Ldv,T,Ldt,C,   &
     &                  Ldc,Work,Ldwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE DLARFB
   END INTERFACE
END MODULE S_DLARFB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARFB_GETT
   INTERFACE
      SUBROUTINE DLARFB_GETT(Ident,M,N,K,T,Ldt,A,Lda,B,Ldb,Work,Ldwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Ident
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: K
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE DLARFB_GETT
   END INTERFACE
END MODULE S_DLARFB_GETT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARF
   INTERFACE
      SUBROUTINE DLARF(Side,M,N,V,Incv,Tau,C,Ldc,Work)
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
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      END SUBROUTINE DLARF
   END INTERFACE
END MODULE S_DLARF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARFG
   INTERFACE
      SUBROUTINE DLARFG(N,Alpha,X,Incx,Tau)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) :: Alpha
      REAL(R8KIND) , DIMENSION(*) :: X
      INTEGER :: Incx
      REAL(R8KIND) , INTENT(OUT) :: Tau
      END SUBROUTINE DLARFG
   END INTERFACE
END MODULE S_DLARFG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARFGP
   INTERFACE
      SUBROUTINE DLARFGP(N,Alpha,X,Incx,Tau)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0 , ONE = 1.0D+0 ,       &
     &                              ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) :: Alpha
      REAL(R8KIND) , DIMENSION(*) :: X
      INTEGER :: Incx
      REAL(R8KIND) , INTENT(INOUT) :: Tau
      END SUBROUTINE DLARFGP
   END INTERFACE
END MODULE S_DLARFGP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARFT
   INTERFACE
      SUBROUTINE DLARFT(Direct,Storev,N,K,V,Ldv,Tau,T,Ldt)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      END SUBROUTINE DLARFT
   END INTERFACE
END MODULE S_DLARFT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARFX
   INTERFACE
      SUBROUTINE DLARFX(Side,M,N,V,Tau,C,Ldc,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: V
      REAL(R8KIND) :: Tau
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      END SUBROUTINE DLARFX
   END INTERFACE
END MODULE S_DLARFX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARFY
   INTERFACE
      SUBROUTINE DLARFY(Uplo,N,V,Incv,Tau,C,Ldc,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0 ,      &
     &                              HALF = 0.5D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: V
      INTEGER :: Incv
      REAL(R8KIND) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      END SUBROUTINE DLARFY
   END INTERFACE
END MODULE S_DLARFY
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARGV
   INTERFACE
      SUBROUTINE DLARGV(N,X,Incx,Y,Incy,C,Incc)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      INTEGER , INTENT(IN) :: Incc
      END SUBROUTINE DLARGV
   END INTERFACE
END MODULE S_DLARGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARNV
   INTERFACE
      SUBROUTINE DLARNV(Idist,Iseed,N,X)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , TWO = 2.0D+0
      INTEGER , PARAMETER  ::  LV = 128
      REAL(R8KIND) , PARAMETER  ::  TWOPI =                             &
     &                       6.28318530717958647692528676655900576839D+0
      INTEGER , INTENT(IN) :: Idist
      INTEGER , DIMENSION(4) :: Iseed
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: X
      END SUBROUTINE DLARNV
   END INTERFACE
END MODULE S_DLARNV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARRA
   INTERFACE
      SUBROUTINE DLARRA(N,D,E,E2,Spltol,Tnrm,Nsplit,Isplit,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E2
      REAL(R8KIND) , INTENT(IN) :: Spltol
      REAL(R8KIND) , INTENT(IN) :: Tnrm
      INTEGER , INTENT(INOUT) :: Nsplit
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Isplit
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLARRA
   END INTERFACE
END MODULE S_DLARRA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARRB
   INTERFACE
      SUBROUTINE DLARRB(N,D,Lld,Ifirst,Ilast,Rtol1,Rtol2,Offset,W,Wgap, &
     &                  Werr,Work,Iwork,Pivmin,Spdiam,Twist,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , TWO = 2.0D0 ,        &
     &                              HALF = 0.5D0
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: Lld
      INTEGER , INTENT(IN) :: Ifirst
      INTEGER , INTENT(IN) :: Ilast
      REAL(R8KIND) , INTENT(IN) :: Rtol1
      REAL(R8KIND) , INTENT(IN) :: Rtol2
      INTEGER , INTENT(IN) :: Offset
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wgap
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      REAL(R8KIND) :: Pivmin
      REAL(R8KIND) , INTENT(IN) :: Spdiam
      INTEGER , INTENT(IN) :: Twist
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLARRB
   END INTERFACE
END MODULE S_DLARRB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARRC
   INTERFACE
      SUBROUTINE DLARRC(Jobt,N,Vl,Vu,D,E,Pivmin,Eigcnt,Lcnt,Rcnt,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0
      CHARACTER :: Jobt
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: E
      REAL(R8KIND) :: Pivmin
      INTEGER , INTENT(OUT) :: Eigcnt
      INTEGER , INTENT(INOUT) :: Lcnt
      INTEGER , INTENT(INOUT) :: Rcnt
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLARRC
   END INTERFACE
END MODULE S_DLARRC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARRD
   INTERFACE
      SUBROUTINE DLARRD(Range,Order,N,Vl,Vu,Il,Iu,Gers,Reltol,D,E,E2,   &
     &                  Pivmin,Nsplit,Isplit,M,W,Werr,Wl,Wu,Iblock,     &
     &                  Indexw,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , HALF = ONE/TWO ,      &
     &                              FUDGE = TWO
      INTEGER , PARAMETER  ::  ALLRNG = 1 , VALRNG = 2 , INDRNG = 3
      CHARACTER :: Range
      CHARACTER :: Order
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Gers
      REAL(R8KIND) , INTENT(IN) :: Reltol
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: E2
      REAL(R8KIND) :: Pivmin
      INTEGER , INTENT(IN) :: Nsplit
      INTEGER , INTENT(IN) , DIMENSION(*) :: Isplit
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL(R8KIND) , INTENT(INOUT) :: Wl
      REAL(R8KIND) , INTENT(INOUT) :: Wu
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indexw
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLARRD
   END INTERFACE
END MODULE S_DLARRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARRE
   INTERFACE
      SUBROUTINE DLARRE(Range,N,Vl,Vu,Il,Iu,D,E,E2,Rtol1,Rtol2,Spltol,  &
     &                  Nsplit,Isplit,M,W,Werr,Wgap,Iblock,Indexw,Gers, &
     &                  Pivmin,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , FOUR = 4.0D0 ,        &
     &                              HNDRD = 100.0D0 , PERT = 8.0D0 ,    &
     &                              HALF = ONE/TWO , FOURTH = ONE/FOUR ,&
     &                              FAC = HALF , MAXGROWTH = 64.0D0 ,   &
     &                              FUDGE = 2.0D0
      INTEGER , PARAMETER  ::  MAXTRY = 6 , ALLRNG = 1 , INDRNG = 2 ,   &
     &                         VALRNG = 3
      CHARACTER :: Range
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) :: Vl
      REAL(R8KIND) , INTENT(INOUT) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: E2
      REAL(R8KIND) :: Rtol1
      REAL(R8KIND) :: Rtol2
      REAL(R8KIND) :: Spltol
      INTEGER :: Nsplit
      INTEGER , DIMENSION(*) :: Isplit
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wgap
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indexw
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Gers
      REAL(R8KIND) , INTENT(INOUT) :: Pivmin
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLARRE
   END INTERFACE
END MODULE S_DLARRE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARRF
   INTERFACE
      SUBROUTINE DLARRF(N,D,L,Ld,Clstrt,Clend,W,Wgap,Werr,Spdiam,Clgapl,&
     &                  Clgapr,Pivmin,Sigma,Dplus,Lplus,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , TWO = 2.0D0 ,         &
     &                              FOUR = 4.0D0 , QUART = 0.25D0 ,     &
     &                              MAXGROWTH1 = 8.D0 ,                 &
     &                              MAXGROWTH2 = 8.D0
      INTEGER , PARAMETER  ::  KTRYMAX = 1 , SLEFT = 1 , SRIGHT = 2
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: L
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Ld
      INTEGER , INTENT(IN) :: Clstrt
      INTEGER , INTENT(IN) :: Clend
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Wgap
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Werr
      REAL(R8KIND) , INTENT(IN) :: Spdiam
      REAL(R8KIND) , INTENT(IN) :: Clgapl
      REAL(R8KIND) , INTENT(IN) :: Clgapr
      REAL(R8KIND) , INTENT(IN) :: Pivmin
      REAL(R8KIND) , INTENT(OUT) :: Sigma
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dplus
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Lplus
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLARRF
   END INTERFACE
END MODULE S_DLARRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARRJ
   INTERFACE
      SUBROUTINE DLARRJ(N,D,E2,Ifirst,Ilast,Rtol,Offset,W,Werr,Work,    &
     &                  Iwork,Pivmin,Spdiam,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , HALF = 0.5D0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: E2
      INTEGER , INTENT(IN) :: Ifirst
      INTEGER , INTENT(IN) :: Ilast
      REAL(R8KIND) , INTENT(IN) :: Rtol
      INTEGER , INTENT(IN) :: Offset
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      REAL(R8KIND) , INTENT(IN) :: Pivmin
      REAL(R8KIND) , INTENT(IN) :: Spdiam
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLARRJ
   END INTERFACE
END MODULE S_DLARRJ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARRK
   INTERFACE
      SUBROUTINE DLARRK(N,Iw,Gl,Gu,D,E2,Pivmin,Reltol,W,Werr,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  HALF = 0.5D0 , TWO = 2.0D0 ,        &
     &                              FUDGE = TWO , ZERO = 0.0D0
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Iw
      REAL(R8KIND) , INTENT(IN) :: Gl
      REAL(R8KIND) , INTENT(IN) :: Gu
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: E2
      REAL(R8KIND) , INTENT(IN) :: Pivmin
      REAL(R8KIND) , INTENT(IN) :: Reltol
      REAL(R8KIND) , INTENT(OUT) :: W
      REAL(R8KIND) , INTENT(OUT) :: Werr
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLARRK
   END INTERFACE
END MODULE S_DLARRK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARRR
   INTERFACE
      SUBROUTINE DLARRR(N,D,E,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , RELCOND = 0.999D0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLARRR
   END INTERFACE
END MODULE S_DLARRR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARRV
   INTERFACE
      SUBROUTINE DLARRV(N,Vl,Vu,D,L,Pivmin,Isplit,M,Dol,Dou,Minrgp,     &
     &                  Rtol1,Rtol2,W,Werr,Wgap,Iblock,Indexw,Gers,Z,   &
     &                  Ldz,Isuppz,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXITR = 10
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , THREE = 3.0D0 ,       &
     &                              FOUR = 4.0D0 , HALF = 0.5D0
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) :: Vu
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: L
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
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isuppz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLARRV
   END INTERFACE
END MODULE S_DLARRV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARSCL2
   INTERFACE
      SUBROUTINE DLARSCL2(M,N,D,X,Ldx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      END SUBROUTINE DLARSCL2
   END INTERFACE
END MODULE S_DLARSCL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARTG
   INTERFACE
      SUBROUTINE DLARTG(F,G,Cs,Sn,R)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0
      REAL(R8KIND) , INTENT(IN) :: F
      REAL(R8KIND) , INTENT(IN) :: G
      REAL(R8KIND) , INTENT(INOUT) :: Cs
      REAL(R8KIND) , INTENT(INOUT) :: Sn
      REAL(R8KIND) , INTENT(INOUT) :: R
      END SUBROUTINE DLARTG
   END INTERFACE
END MODULE S_DLARTG
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARTGP
   INTERFACE
      SUBROUTINE DLARTGP(F,G,Cs,Sn,R)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0
      REAL(R8KIND) , INTENT(IN) :: F
      REAL(R8KIND) , INTENT(IN) :: G
      REAL(R8KIND) , INTENT(INOUT) :: Cs
      REAL(R8KIND) , INTENT(INOUT) :: Sn
      REAL(R8KIND) , INTENT(INOUT) :: R
      END SUBROUTINE DLARTGP
   END INTERFACE
END MODULE S_DLARTGP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARTGS
   INTERFACE
      SUBROUTINE DLARTGS(X,Y,Sigma,Cs,Sn)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  NEGONE = -1.0D0 , ONE = 1.0D0 ,     &
     &                              ZERO = 0.0D0
      REAL(R8KIND) , INTENT(IN) :: X
      REAL(R8KIND) , INTENT(IN) :: Y
      REAL(R8KIND) , INTENT(IN) :: Sigma
      REAL(R8KIND) :: Cs
      REAL(R8KIND) :: Sn
      END SUBROUTINE DLARTGS
   END INTERFACE
END MODULE S_DLARTGS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARTV
   INTERFACE
      SUBROUTINE DLARTV(N,X,Incx,Y,Incy,C,S,Incc)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      INTEGER , INTENT(IN) :: Incc
      END SUBROUTINE DLARTV
   END INTERFACE
END MODULE S_DLARTV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARUV
   INTERFACE
      SUBROUTINE DLARUV(Iseed,N,X)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0
      INTEGER , PARAMETER  ::  LV = 128 , IPW2 = 4096
      REAL(R8KIND) , PARAMETER  ::  R = ONE/IPW2
      INTEGER , INTENT(INOUT) , DIMENSION(4) :: Iseed
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N) :: X
      END SUBROUTINE DLARUV
   END INTERFACE
END MODULE S_DLARUV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARZB
   INTERFACE
      SUBROUTINE DLARZB(Side,Trans,Direct,Storev,M,N,K,L,V,Ldv,T,Ldt,C, &
     &                  Ldc,Work,Ldwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: L
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE DLARZB
   END INTERFACE
END MODULE S_DLARZB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARZ
   INTERFACE
      SUBROUTINE DLARZ(Side,M,N,L,V,Incv,Tau,C,Ldc,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Side
      INTEGER :: M
      INTEGER :: N
      INTEGER :: L
      REAL(R8KIND) , DIMENSION(*) :: V
      INTEGER :: Incv
      REAL(R8KIND) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      END SUBROUTINE DLARZ
   END INTERFACE
END MODULE S_DLARZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLARZT
   INTERFACE
      SUBROUTINE DLARZT(Direct,Storev,N,K,V,Ldv,Tau,T,Ldt)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      END SUBROUTINE DLARZT
   END INTERFACE
END MODULE S_DLARZT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAS2
   INTERFACE
      SUBROUTINE DLAS2(F,G,H,Ssmin,Ssmax)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0
      REAL(R8KIND) , INTENT(IN) :: F
      REAL(R8KIND) , INTENT(IN) :: G
      REAL(R8KIND) , INTENT(IN) :: H
      REAL(R8KIND) , INTENT(INOUT) :: Ssmin
      REAL(R8KIND) , INTENT(OUT) :: Ssmax
      END SUBROUTINE DLAS2
   END INTERFACE
END MODULE S_DLAS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASCL2
   INTERFACE
      SUBROUTINE DLASCL2(M,N,D,X,Ldx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      END SUBROUTINE DLASCL2
   END INTERFACE
END MODULE S_DLASCL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASCL
   INTERFACE
      SUBROUTINE DLASCL(Type,Kl,Ku,Cfrom,Cto,M,N,A,Lda,Info)
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
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASCL
   END INTERFACE
END MODULE S_DLASCL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASD0
   INTERFACE
      SUBROUTINE DLASD0(N,Sqre,D,E,U,Ldu,Vt,Ldvt,Smlsiz,Iwork,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Sqre
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      INTEGER :: Smlsiz
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASD0
   END INTERFACE
END MODULE S_DLASD0
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASD1
   INTERFACE
      SUBROUTINE DLASD1(Nl,Nr,Sqre,D,Alpha,Beta,U,Ldu,Vt,Ldvt,Idxq,     &
     &                  Iwork,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER :: Nl
      INTEGER :: Nr
      INTEGER :: Sqre
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) :: Alpha
      REAL(R8KIND) , INTENT(INOUT) :: Beta
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      INTEGER , DIMENSION(*) :: Idxq
      INTEGER , DIMENSION(*) :: Iwork
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASD1
   END INTERFACE
END MODULE S_DLASD1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASD2
   INTERFACE
      SUBROUTINE DLASD2(Nl,Nr,Sqre,K,D,Z,Alpha,Beta,U,Ldu,Vt,Ldvt,      &
     &                  Dsigma,U2,Ldu2,Vt2,Ldvt2,Idxp,Idx,Idxc,Idxq,    &
     &                  Coltyp,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0 , EIGHT = 8.0D+0
      INTEGER :: Nl
      INTEGER :: Nr
      INTEGER , INTENT(IN) :: Sqre
      INTEGER , INTENT(INOUT) :: K
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dsigma
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      REAL(R8KIND) , DIMENSION(Ldvt2,*) :: Vt2
      INTEGER :: Ldvt2
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxp
      INTEGER , DIMENSION(*) :: Idx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxc
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxq
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Coltyp
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASD2
   END INTERFACE
END MODULE S_DLASD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASD3
   INTERFACE
      SUBROUTINE DLASD3(Nl,Nr,Sqre,K,D,Q,Ldq,Dsigma,U,Ldu,U2,Ldu2,Vt,   &
     &                  Ldvt,Vt2,Ldvt2,Idxc,Ctot,Z,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0 ,      &
     &                              NEGONE = -1.0D+0
      INTEGER :: Nl
      INTEGER :: Nr
      INTEGER , INTENT(IN) :: Sqre
      INTEGER :: K
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dsigma
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvt2,*) :: Vt2
      INTEGER :: Ldvt2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Idxc
      INTEGER , DIMENSION(*) :: Ctot
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASD3
   END INTERFACE
END MODULE S_DLASD3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASD4
   INTERFACE
      SUBROUTINE DLASD4(N,I,D,Z,Delta,Rho,Sigma,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXIT = 400
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0 , THREE = 3.0D+0 ,     &
     &                              FOUR = 4.0D+0 , EIGHT = 8.0D+0 ,    &
     &                              TEN = 10.0D+0
      INTEGER , INTENT(IN) :: N
      INTEGER :: I
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: Z
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Delta
      REAL(R8KIND) :: Rho
      REAL(R8KIND) , INTENT(INOUT) :: Sigma
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASD4
   END INTERFACE
END MODULE S_DLASD4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASD5
   INTERFACE
      SUBROUTINE DLASD5(I,D,Z,Delta,Rho,Dsigma,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0 , THREE = 3.0D+0 ,     &
     &                              FOUR = 4.0D+0
      INTEGER , INTENT(IN) :: I
      REAL(R8KIND) , INTENT(IN) , DIMENSION(2) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(2) :: Z
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(2) :: Delta
      REAL(R8KIND) , INTENT(IN) :: Rho
      REAL(R8KIND) , INTENT(OUT) :: Dsigma
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(2) :: Work
      END SUBROUTINE DLASD5
   END INTERFACE
END MODULE S_DLASD5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASD6
   INTERFACE
      SUBROUTINE DLASD6(Icompq,Nl,Nr,Sqre,D,Vf,Vl,Alpha,Beta,Idxq,Perm, &
     &                  Givptr,Givcol,Ldgcol,Givnum,Ldgnum,Poles,Difl,  &
     &                  Difr,Z,K,C,S,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER :: Icompq
      INTEGER :: Nl
      INTEGER :: Nr
      INTEGER :: Sqre
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: Vf
      REAL(R8KIND) , DIMENSION(*) :: Vl
      REAL(R8KIND) , INTENT(INOUT) :: Alpha
      REAL(R8KIND) , INTENT(INOUT) :: Beta
      INTEGER , DIMENSION(*) :: Idxq
      INTEGER , DIMENSION(*) :: Perm
      INTEGER :: Givptr
      INTEGER , DIMENSION(Ldgcol,*) :: Givcol
      INTEGER :: Ldgcol
      REAL(R8KIND) , DIMENSION(Ldgnum,*) :: Givnum
      INTEGER :: Ldgnum
      REAL(R8KIND) , DIMENSION(Ldgnum,*) :: Poles
      REAL(R8KIND) , DIMENSION(*) :: Difl
      REAL(R8KIND) , DIMENSION(*) :: Difr
      REAL(R8KIND) , DIMENSION(*) :: Z
      INTEGER :: K
      REAL(R8KIND) :: C
      REAL(R8KIND) :: S
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASD6
   END INTERFACE
END MODULE S_DLASD6
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASD7
   INTERFACE
      SUBROUTINE DLASD7(Icompq,Nl,Nr,Sqre,K,D,Z,Zw,Vf,Vfw,Vl,Vlw,Alpha, &
     &                  Beta,Dsigma,Idx,Idxp,Idxq,Perm,Givptr,Givcol,   &
     &                  Ldgcol,Givnum,Ldgnum,C,S,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0 , EIGHT = 8.0D+0
      INTEGER , INTENT(IN) :: Icompq
      INTEGER :: Nl
      INTEGER :: Nr
      INTEGER , INTENT(IN) :: Sqre
      INTEGER , INTENT(INOUT) :: K
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Zw
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vf
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vfw
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vlw
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dsigma
      INTEGER , DIMENSION(*) :: Idx
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxp
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idxq
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Perm
      INTEGER , INTENT(INOUT) :: Givptr
      INTEGER , INTENT(OUT) , DIMENSION(Ldgcol,*) :: Givcol
      INTEGER , INTENT(IN) :: Ldgcol
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Ldgnum,*) :: Givnum
      INTEGER , INTENT(IN) :: Ldgnum
      REAL(R8KIND) , INTENT(INOUT) :: C
      REAL(R8KIND) , INTENT(INOUT) :: S
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASD7
   END INTERFACE
END MODULE S_DLASD7
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASD8
   INTERFACE
      SUBROUTINE DLASD8(Icompq,K,D,Z,Vf,Vl,Difl,Difr,Lddifr,Dsigma,Work,&
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER , INTENT(IN) :: Icompq
      INTEGER :: K
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL(R8KIND) , DIMENSION(*) :: Vf
      REAL(R8KIND) , DIMENSION(*) :: Vl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Difl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lddifr,*) :: Difr
      INTEGER , INTENT(IN) :: Lddifr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dsigma
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASD8
   END INTERFACE
END MODULE S_DLASD8
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASDA
   INTERFACE
      SUBROUTINE DLASDA(Icompq,Smlsiz,N,Sqre,D,E,U,Ldu,Vt,K,Difl,Difr,Z,&
     &                  Poles,Givptr,Givcol,Ldgcol,Perm,Givnum,C,S,Work,&
     &                  Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER :: Icompq
      INTEGER :: Smlsiz
      INTEGER :: N
      INTEGER :: Sqre
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
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
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASDA
   END INTERFACE
END MODULE S_DLASDA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASDQ
   INTERFACE
      SUBROUTINE DLASDQ(Uplo,Sqre,N,Ncvt,Nru,Ncc,D,E,Vt,Ldvt,U,Ldu,C,   &
     &                  Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: Sqre
      INTEGER :: N
      INTEGER :: Ncvt
      INTEGER :: Nru
      INTEGER :: Ncc
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldvt,*) :: Vt
      INTEGER :: Ldvt
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASDQ
   END INTERFACE
END MODULE S_DLASDQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASDT
   INTERFACE
      SUBROUTINE DLASDT(N,Lvl,Nd,Inode,Ndiml,Ndimr,Msub)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(INOUT) :: Lvl
      INTEGER , INTENT(OUT) :: Nd
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Inode
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ndiml
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ndimr
      INTEGER , INTENT(IN) :: Msub
      END SUBROUTINE DLASDT
   END INTERFACE
END MODULE S_DLASDT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASET
   INTERFACE
      SUBROUTINE DLASET(Uplo,M,N,Alpha,Beta,A,Lda)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE DLASET
   END INTERFACE
END MODULE S_DLASET
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASQ1
   INTERFACE
      SUBROUTINE DLASQ1(N,D,E,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASQ1
   END INTERFACE
END MODULE S_DLASQ1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASQ2
   INTERFACE
      SUBROUTINE DLASQ2(N,Z,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  CBIAS = 1.50D0 , ZERO = 0.0D0 ,     &
     &                              HALF = 0.5D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , FOUR = 4.0D0 ,        &
     &                              HUNDRD = 100.0D0
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLASQ2
   END INTERFACE
END MODULE S_DLASQ2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASQ3
   INTERFACE
      SUBROUTINE DLASQ3(I0,N0,Z,Pp,Dmin,Sigma,Desig,Qmax,Nfail,Iter,    &
     &                  Ndiv,Ieee,Ttype,Dmin1,Dmin2,Dn,Dn1,Dn2,G,Tau)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  CBIAS = 1.50D0 , ZERO = 0.0D0 ,     &
     &                              QURTR = 0.250D0 , HALF = 0.5D0 ,    &
     &                              ONE = 1.0D0 , TWO = 2.0D0 ,         &
     &                              HUNDRD = 100.0D0
      INTEGER :: I0
      INTEGER , INTENT(INOUT) :: N0
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(INOUT) :: Pp
      REAL(R8KIND) , INTENT(INOUT) :: Dmin
      REAL(R8KIND) , INTENT(INOUT) :: Sigma
      REAL(R8KIND) , INTENT(INOUT) :: Desig
      REAL(R8KIND) , INTENT(INOUT) :: Qmax
      INTEGER , INTENT(INOUT) :: Nfail
      INTEGER , INTENT(INOUT) :: Iter
      INTEGER , INTENT(INOUT) :: Ndiv
      LOGICAL :: Ieee
      INTEGER , INTENT(INOUT) :: Ttype
      REAL(R8KIND) :: Dmin1
      REAL(R8KIND) , INTENT(INOUT) :: Dmin2
      REAL(R8KIND) :: Dn
      REAL(R8KIND) :: Dn1
      REAL(R8KIND) :: Dn2
      REAL(R8KIND) :: G
      REAL(R8KIND) , INTENT(INOUT) :: Tau
      END SUBROUTINE DLASQ3
   END INTERFACE
END MODULE S_DLASQ3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASQ4
   INTERFACE
      SUBROUTINE DLASQ4(I0,N0,Z,Pp,N0in,Dmin,Dmin1,Dmin2,Dn,Dn1,Dn2,Tau,&
     &                  Ttype,G)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  CNST1 = 0.5630D0 , CNST2 = 1.010D0 ,&
     &                              CNST3 = 1.050D0 , QURTR = 0.250D0 , &
     &                              THIRD = 0.3330D0 , HALF = 0.50D0 ,  &
     &                              ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , HUNDRD = 100.0D0
      INTEGER , INTENT(IN) :: I0
      INTEGER , INTENT(IN) :: N0
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Z
      INTEGER , INTENT(IN) :: Pp
      INTEGER , INTENT(IN) :: N0in
      REAL(R8KIND) , INTENT(IN) :: Dmin
      REAL(R8KIND) , INTENT(IN) :: Dmin1
      REAL(R8KIND) , INTENT(IN) :: Dmin2
      REAL(R8KIND) , INTENT(IN) :: Dn
      REAL(R8KIND) , INTENT(IN) :: Dn1
      REAL(R8KIND) , INTENT(IN) :: Dn2
      REAL(R8KIND) , INTENT(OUT) :: Tau
      INTEGER , INTENT(INOUT) :: Ttype
      REAL(R8KIND) , INTENT(INOUT) :: G
      END SUBROUTINE DLASQ4
   END INTERFACE
END MODULE S_DLASQ4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASQ5
   INTERFACE
      SUBROUTINE DLASQ5(I0,N0,Z,Pp,Tau,Sigma,Dmin,Dmin1,Dmin2,Dn,Dnm1,  &
     &                  Dnm2,Ieee,Eps)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , HALF = 0.5
      INTEGER , INTENT(IN) :: I0
      INTEGER , INTENT(IN) :: N0
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(IN) :: Pp
      REAL(R8KIND) , INTENT(INOUT) :: Tau
      REAL(R8KIND) , INTENT(IN) :: Sigma
      REAL(R8KIND) , INTENT(INOUT) :: Dmin
      REAL(R8KIND) , INTENT(OUT) :: Dmin1
      REAL(R8KIND) , INTENT(OUT) :: Dmin2
      REAL(R8KIND) , INTENT(INOUT) :: Dn
      REAL(R8KIND) , INTENT(INOUT) :: Dnm1
      REAL(R8KIND) , INTENT(INOUT) :: Dnm2
      LOGICAL , INTENT(IN) :: Ieee
      REAL(R8KIND) , INTENT(IN) :: Eps
      END SUBROUTINE DLASQ5
   END INTERFACE
END MODULE S_DLASQ5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASQ6
   INTERFACE
      SUBROUTINE DLASQ6(I0,N0,Z,Pp,Dmin,Dmin1,Dmin2,Dn,Dnm1,Dnm2)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0
      INTEGER , INTENT(IN) :: I0
      INTEGER , INTENT(IN) :: N0
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(IN) :: Pp
      REAL(R8KIND) , INTENT(INOUT) :: Dmin
      REAL(R8KIND) , INTENT(OUT) :: Dmin1
      REAL(R8KIND) , INTENT(OUT) :: Dmin2
      REAL(R8KIND) , INTENT(INOUT) :: Dn
      REAL(R8KIND) , INTENT(INOUT) :: Dnm1
      REAL(R8KIND) , INTENT(INOUT) :: Dnm2
      END SUBROUTINE DLASQ6
   END INTERFACE
END MODULE S_DLASQ6
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASR
   INTERFACE
      SUBROUTINE DLASR(Side,Pivot,Direct,M,N,C,S,A,Lda)
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
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      END SUBROUTINE DLASR
   END INTERFACE
END MODULE S_DLASR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASRT
   INTERFACE
      SUBROUTINE DLASRT(Id,N,D,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  SELECT = 20
      CHARACTER :: Id
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASRT
   END INTERFACE
END MODULE S_DLASRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASSQ
   INTERFACE
      SUBROUTINE DLASSQ(N,X,Incx,Scale,Sumsq)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) :: Sumsq
      END SUBROUTINE DLASSQ
   END INTERFACE
END MODULE S_DLASSQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASV2
   INTERFACE
      SUBROUTINE DLASV2(F,G,H,Ssmin,Ssmax,Snr,Csr,Snl,Csl)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , HALF = 0.5D0 ,       &
     &                              ONE = 1.0D0 , TWO = 2.0D0 ,         &
     &                              FOUR = 4.0D0
      REAL(R8KIND) , INTENT(IN) :: F
      REAL(R8KIND) , INTENT(IN) :: G
      REAL(R8KIND) , INTENT(IN) :: H
      REAL(R8KIND) , INTENT(INOUT) :: Ssmin
      REAL(R8KIND) , INTENT(INOUT) :: Ssmax
      REAL(R8KIND) , INTENT(INOUT) :: Snr
      REAL(R8KIND) , INTENT(INOUT) :: Csr
      REAL(R8KIND) , INTENT(INOUT) :: Snl
      REAL(R8KIND) , INTENT(INOUT) :: Csl
      END SUBROUTINE DLASV2
   END INTERFACE
END MODULE S_DLASV2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASWLQ
   INTERFACE
      SUBROUTINE DLASWLQ(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER :: Nb
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASWLQ
   END INTERFACE
END MODULE S_DLASWLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASWP
   INTERFACE
      SUBROUTINE DLASWP(N,A,Lda,K1,K2,Ipiv,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) :: K1
      INTEGER , INTENT(IN) :: K2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(IN) :: Incx
      END SUBROUTINE DLASWP
   END INTERFACE
END MODULE S_DLASWP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASY2
   INTERFACE
      SUBROUTINE DLASY2(Ltranl,Ltranr,Isgn,N1,N2,Tl,Ldtl,Tr,Ldtr,B,Ldb, &
     &                  Scale,X,Ldx,Xnorm,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0 , HALF = 0.5D+0 ,      &
     &                              EIGHT = 8.0D+0
      LOGICAL , INTENT(IN) :: Ltranl
      LOGICAL , INTENT(IN) :: Ltranr
      INTEGER , INTENT(IN) :: Isgn
      INTEGER , INTENT(IN) :: N1
      INTEGER , INTENT(IN) :: N2
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldtl,*) :: Tl
      INTEGER , INTENT(IN) :: Ldtl
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldtr,*) :: Tr
      INTEGER , INTENT(IN) :: Ldtr
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(OUT) :: Xnorm
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLASY2
   END INTERFACE
END MODULE S_DLASY2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_SYAMV
   INTERFACE
      SUBROUTINE DLA_SYAMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      END SUBROUTINE DLA_SYAMV
   END INTERFACE
END MODULE S_DLA_SYAMV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASYF_AA
   INTERFACE
      SUBROUTINE DLASYF_AA(Uplo,J1,M,Nb,A,Lda,Ipiv,H,Ldh,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: J1
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: Nb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END SUBROUTINE DLASYF_AA
   END INTERFACE
END MODULE S_DLASYF_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASYF
   INTERFACE
      SUBROUTINE DLASYF(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASYF
   END INTERFACE
END MODULE S_DLASYF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASYF_RK
   INTERFACE
      SUBROUTINE DLASYF_RK(Uplo,N,Nb,Kb,A,Lda,E,Ipiv,W,Ldw,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASYF_RK
   END INTERFACE
END MODULE S_DLASYF_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLASYF_ROOK
   INTERFACE
      SUBROUTINE DLASYF_ROOK(Uplo,N,Nb,Kb,A,Lda,Ipiv,W,Ldw,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(OUT) :: Kb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLASYF_ROOK
   END INTERFACE
END MODULE S_DLASYF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_SYRCOND
   INTERFACE
      FUNCTION DLA_SYRCOND(Uplo,N,A,Lda,Af,Ldaf,Ipiv,Cmode,C,Info,Work, &
     &                     Iwork)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DLA_SYRCOND
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(IN) :: Cmode
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      INTEGER , INTENT(INOUT) :: Info
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      END FUNCTION DLA_SYRCOND
   END INTERFACE
END MODULE S_DLA_SYRCOND
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_SYRFSX_EXTENDED
   INTERFACE
      SUBROUTINE DLA_SYRFSX_EXTENDED(Prec_type,Uplo,N,Nrhs,A,Lda,Af,    &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL(R8KIND) , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      REAL(R8KIND) , DIMENSION(*) :: Res
      REAL(R8KIND) , DIMENSION(*) :: Ayb
      REAL(R8KIND) , DIMENSION(*) :: Dy
      REAL(R8KIND) , DIMENSION(*) :: Y_tail
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL(R8KIND) , INTENT(IN) :: Rthresh
      REAL(R8KIND) , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLA_SYRFSX_EXTENDED
   END INTERFACE
END MODULE S_DLA_SYRFSX_EXTENDED
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_SYRPVGRW
   INTERFACE
      FUNCTION DLA_SYRPVGRW(Uplo,N,Info,A,Lda,Af,Ldaf,Ipiv,Work)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DLA_SYRPVGRW
      CHARACTER(1) :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Info
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldaf,*) :: Af
      INTEGER , INTENT(IN) :: Ldaf
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      END FUNCTION DLA_SYRPVGRW
   END INTERFACE
END MODULE S_DLA_SYRPVGRW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAT2S
   INTERFACE
      SUBROUTINE DLAT2S(Uplo,N,A,Lda,Sa,Ldsa,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , INTENT(OUT) , DIMENSION(Ldsa,*) :: Sa
      INTEGER , INTENT(IN) :: Ldsa
      INTEGER , INTENT(OUT) :: Info
      END SUBROUTINE DLAT2S
   END INTERFACE
END MODULE S_DLAT2S
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLATBS
   INTERFACE
      SUBROUTINE DLATBS(Uplo,Trans,Diag,Normin,N,Kd,Ab,Ldab,X,Scale,    &
     &                  Cnorm,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , HALF = 0.5D+0 ,     &
     &                              ONE = 1.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      INTEGER :: Kd
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLATBS
   END INTERFACE
END MODULE S_DLATBS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLATDF
   INTERFACE
      SUBROUTINE DLATDF(Ijob,N,Z,Ldz,Rhs,Rdsum,Rdscal,Ipiv,Jpiv)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXDIM = 8
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , INTENT(IN) :: Ijob
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rhs
      REAL(R8KIND) :: Rdsum
      REAL(R8KIND) :: Rdscal
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Jpiv
      END SUBROUTINE DLATDF
   END INTERFACE
END MODULE S_DLATDF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLATPS
   INTERFACE
      SUBROUTINE DLATPS(Uplo,Trans,Diag,Normin,N,Ap,X,Scale,Cnorm,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , HALF = 0.5D+0 ,     &
     &                              ONE = 1.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLATPS
   END INTERFACE
END MODULE S_DLATPS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLATRD
   INTERFACE
      SUBROUTINE DLATRD(Uplo,N,Nb,A,Lda,E,Tau,W,Ldw)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              HALF = 0.5D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldw,*) :: W
      INTEGER :: Ldw
      END SUBROUTINE DLATRD
   END INTERFACE
END MODULE S_DLATRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLATRS
   INTERFACE
      SUBROUTINE DLATRS(Uplo,Trans,Diag,Normin,N,A,Lda,X,Scale,Cnorm,   &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , HALF = 0.5D+0 ,     &
     &                              ONE = 1.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      CHARACTER :: Normin
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Cnorm
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLATRS
   END INTERFACE
END MODULE S_DLATRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLATRZ
   INTERFACE
      SUBROUTINE DLATRZ(M,N,L,A,Lda,Tau,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: L
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      END SUBROUTINE DLATRZ
   END INTERFACE
END MODULE S_DLATRZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLATSQR
   INTERFACE
      SUBROUTINE DLATSQR(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER :: Nb
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLATSQR
   END INTERFACE
END MODULE S_DLATSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAUU2
   INTERFACE
      SUBROUTINE DLAUU2(Uplo,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DLAUU2
   END INTERFACE
END MODULE S_DLAUU2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLAUUM
   INTERFACE
      SUBROUTINE DLAUUM(Uplo,N,A,Lda,Info)
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
      END SUBROUTINE DLAUUM
   END INTERFACE
END MODULE S_DLAUUM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DLA_WWADDW
   INTERFACE
      SUBROUTINE DLA_WWADDW(N,X,Y,W)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: W
      END SUBROUTINE DLA_WWADDW
   END INTERFACE
END MODULE S_DLA_WWADDW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DOPGTR
   INTERFACE
      SUBROUTINE DOPGTR(Uplo,N,Ap,Tau,Q,Ldq,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DOPGTR
   END INTERFACE
END MODULE S_DOPGTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DOPMTR
   INTERFACE
      SUBROUTINE DOPMTR(Side,Uplo,Trans,M,N,Ap,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DOPMTR
   END INTERFACE
END MODULE S_DOPMTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORBDB1
   INTERFACE
      SUBROUTINE DORBDB1(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Phi
      REAL(R8KIND) , DIMENSION(*) :: Taup1
      REAL(R8KIND) , DIMENSION(*) :: Taup2
      REAL(R8KIND) , DIMENSION(*) :: Tauq1
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORBDB1
   END INTERFACE
END MODULE S_DORBDB1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORBDB2
   INTERFACE
      SUBROUTINE DORBDB2(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  NEGONE = -1.0D0 , ONE = 1.0D0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Theta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Phi
      REAL(R8KIND) , DIMENSION(*) :: Taup1
      REAL(R8KIND) , DIMENSION(*) :: Taup2
      REAL(R8KIND) , DIMENSION(*) :: Tauq1
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORBDB2
   END INTERFACE
END MODULE S_DORBDB2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORBDB3
   INTERFACE
      SUBROUTINE DORBDB3(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Theta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Phi
      REAL(R8KIND) , DIMENSION(*) :: Taup1
      REAL(R8KIND) , DIMENSION(*) :: Taup2
      REAL(R8KIND) , DIMENSION(*) :: Tauq1
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORBDB3
   END INTERFACE
END MODULE S_DORBDB3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORBDB4
   INTERFACE
      SUBROUTINE DORBDB4(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Phi,Taup1,     &
     &                   Taup2,Tauq1,Phantom,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  NEGONE = -1.0D0 , ONE = 1.0D0 ,     &
     &                              ZERO = 0.0D0
      INTEGER , INTENT(IN) :: M
      INTEGER :: P
      INTEGER :: Q
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Phi
      REAL(R8KIND) , DIMENSION(*) :: Taup1
      REAL(R8KIND) , DIMENSION(*) :: Taup2
      REAL(R8KIND) , DIMENSION(*) :: Tauq1
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Phantom
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORBDB4
   END INTERFACE
END MODULE S_DORBDB4
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORBDB5
   INTERFACE
      SUBROUTINE DORBDB5(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,&
     &                   Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0
      INTEGER :: M1
      INTEGER :: M2
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: X1
      INTEGER :: Incx1
      REAL(R8KIND) , DIMENSION(*) :: X2
      INTEGER :: Incx2
      REAL(R8KIND) , DIMENSION(Ldq1,*) :: Q1
      INTEGER :: Ldq1
      REAL(R8KIND) , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORBDB5
   END INTERFACE
END MODULE S_DORBDB5
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORBDB6
   INTERFACE
      SUBROUTINE DORBDB6(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,&
     &                   Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ALPHASQ = 0.01D0 , REALONE = 1.0D0 ,&
     &                              REALZERO = 0.0D0 , NEGONE = -1.0D0 ,&
     &                              ONE = 1.0D0 , ZERO = 0.0D0
      INTEGER :: M1
      INTEGER :: M2
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: X1
      INTEGER :: Incx1
      REAL(R8KIND) , DIMENSION(*) :: X2
      INTEGER :: Incx2
      REAL(R8KIND) , DIMENSION(Ldq1,*) :: Q1
      INTEGER :: Ldq1
      REAL(R8KIND) , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORBDB6
   END INTERFACE
END MODULE S_DORBDB6
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORBDB
   INTERFACE
      SUBROUTINE DORBDB(Trans,Signs,M,P,Q,X11,Ldx11,X12,Ldx12,X21,Ldx21,&
     &                  X22,Ldx22,Theta,Phi,Taup1,Taup2,Tauq1,Tauq2,    &
     &                  Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  REALONE = 1.0D0 , ONE = 1.0D0
      CHARACTER :: Trans
      CHARACTER :: Signs
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: P
      INTEGER , INTENT(IN) :: Q
      REAL(R8KIND) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL(R8KIND) , DIMENSION(Ldx12,*) :: X12
      INTEGER :: Ldx12
      REAL(R8KIND) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL(R8KIND) , DIMENSION(Ldx22,*) :: X22
      INTEGER :: Ldx22
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Theta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Phi
      REAL(R8KIND) , DIMENSION(*) :: Taup1
      REAL(R8KIND) , DIMENSION(*) :: Taup2
      REAL(R8KIND) , DIMENSION(*) :: Tauq1
      REAL(R8KIND) , DIMENSION(*) :: Tauq2
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORBDB
   END INTERFACE
END MODULE S_DORBDB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORCSD2BY1
   INTERFACE
      SUBROUTINE DORCSD2BY1(Jobu1,Jobu2,Jobv1t,M,P,Q,X11,Ldx11,X21,     &
     &                      Ldx21,Theta,U1,Ldu1,U2,Ldu2,V1t,Ldv1t,Work, &
     &                      Lwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0
      CHARACTER :: Jobu1
      CHARACTER :: Jobu2
      CHARACTER :: Jobv1t
      INTEGER :: M
      INTEGER :: P
      INTEGER :: Q
      REAL(R8KIND) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL(R8KIND) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL(R8KIND) , DIMENSION(*) :: Theta
      REAL(R8KIND) , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      REAL(R8KIND) , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      REAL(R8KIND) , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORCSD2BY1
   END INTERFACE
END MODULE S_DORCSD2BY1
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORCSD
   INTERFACE
      RECURSIVE SUBROUTINE DORCSD(Jobu1,Jobu2,Jobv1t,Jobv2t,Trans,Signs,&
     &                            M,P,Q,X11,Ldx11,X12,Ldx12,X21,Ldx21,  &
     &                            X22,Ldx22,Theta,U1,Ldu1,U2,Ldu2,V1t,  &
     &                            Ldv1t,V2t,Ldv2t,Work,Lwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0
      CHARACTER :: Jobu1
      CHARACTER :: Jobu2
      CHARACTER :: Jobv1t
      CHARACTER :: Jobv2t
      CHARACTER :: Trans
      CHARACTER :: Signs
      INTEGER :: M
      INTEGER :: P
      INTEGER :: Q
      REAL(R8KIND) , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      REAL(R8KIND) , DIMENSION(Ldx12,*) :: X12
      INTEGER :: Ldx12
      REAL(R8KIND) , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL(R8KIND) , DIMENSION(Ldx22,*) :: X22
      INTEGER :: Ldx22
      REAL(R8KIND) , DIMENSION(*) :: Theta
      REAL(R8KIND) , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      REAL(R8KIND) , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      REAL(R8KIND) , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      REAL(R8KIND) , DIMENSION(Ldv2t,*) :: V2t
      INTEGER :: Ldv2t
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORCSD
   END INTERFACE
END MODULE S_DORCSD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORG2L
   INTERFACE
      SUBROUTINE DORG2L(M,N,K,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORG2L
   END INTERFACE
END MODULE S_DORG2L
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORG2R
   INTERFACE
      SUBROUTINE DORG2R(M,N,K,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORG2R
   END INTERFACE
END MODULE S_DORG2R
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORGBR
   INTERFACE
      SUBROUTINE DORGBR(Vect,M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Vect
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORGBR
   END INTERFACE
END MODULE S_DORGBR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORGHR
   INTERFACE
      SUBROUTINE DORGHR(N,Ilo,Ihi,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORGHR
   END INTERFACE
END MODULE S_DORGHR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORGL2
   INTERFACE
      SUBROUTINE DORGL2(M,N,K,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORGL2
   END INTERFACE
END MODULE S_DORGL2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORGLQ
   INTERFACE
      SUBROUTINE DORGLQ(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORGLQ
   END INTERFACE
END MODULE S_DORGLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORGQL
   INTERFACE
      SUBROUTINE DORGQL(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORGQL
   END INTERFACE
END MODULE S_DORGQL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORGQR
   INTERFACE
      SUBROUTINE DORGQR(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORGQR
   END INTERFACE
END MODULE S_DORGQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORGR2
   INTERFACE
      SUBROUTINE DORGR2(M,N,K,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORGR2
   END INTERFACE
END MODULE S_DORGR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORGRQ
   INTERFACE
      SUBROUTINE DORGRQ(M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORGRQ
   END INTERFACE
END MODULE S_DORGRQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORGTR
   INTERFACE
      SUBROUTINE DORGTR(Uplo,N,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORGTR
   END INTERFACE
END MODULE S_DORGTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORGTSQR
   INTERFACE
      SUBROUTINE DORGTSQR(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER , INTENT(IN) :: Nb
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORGTSQR
   END INTERFACE
END MODULE S_DORGTSQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORGTSQR_ROW
   INTERFACE
      SUBROUTINE DORGTSQR_ROW(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Mb
      INTEGER , INTENT(IN) :: Nb
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORGTSQR_ROW
   END INTERFACE
END MODULE S_DORGTSQR_ROW
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORHR_COL
   INTERFACE
      SUBROUTINE DORHR_COL(M,N,Nb,A,Lda,T,Ldt,D,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nb
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORHR_COL
   END INTERFACE
END MODULE S_DORHR_COL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORM22
   INTERFACE
      SUBROUTINE DORM22(Side,Trans,M,N,N1,N2,Q,Ldq,C,Ldc,Work,Lwork,    &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: N1
      INTEGER :: N2
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORM22
   END INTERFACE
END MODULE S_DORM22
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORM2L
   INTERFACE
      SUBROUTINE DORM2L(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORM2L
   END INTERFACE
END MODULE S_DORM2L
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORM2R
   INTERFACE
      SUBROUTINE DORM2R(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORM2R
   END INTERFACE
END MODULE S_DORM2R
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORMBR
   INTERFACE
      SUBROUTINE DORMBR(Vect,Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,     &
     &                  Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Vect
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORMBR
   END INTERFACE
END MODULE S_DORMBR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORMHR
   INTERFACE
      SUBROUTINE DORMHR(Side,Trans,M,N,Ilo,Ihi,A,Lda,Tau,C,Ldc,Work,    &
     &                  Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORMHR
   END INTERFACE
END MODULE S_DORMHR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORML2
   INTERFACE
      SUBROUTINE DORML2(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORML2
   END INTERFACE
END MODULE S_DORML2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORMLQ
   INTERFACE
      SUBROUTINE DORMLQ(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORMLQ
   END INTERFACE
END MODULE S_DORMLQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORMQL
   INTERFACE
      SUBROUTINE DORMQL(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORMQL
   END INTERFACE
END MODULE S_DORMQL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORMQR
   INTERFACE
      SUBROUTINE DORMQR(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORMQR
   END INTERFACE
END MODULE S_DORMQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORMR2
   INTERFACE
      SUBROUTINE DORMR2(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORMR2
   END INTERFACE
END MODULE S_DORMR2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORMR3
   INTERFACE
      SUBROUTINE DORMR3(Side,Trans,M,N,K,L,A,Lda,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      INTEGER :: L
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORMR3
   END INTERFACE
END MODULE S_DORMR3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORMRQ
   INTERFACE
      SUBROUTINE DORMRQ(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORMRQ
   END INTERFACE
END MODULE S_DORMRQ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORMRZ
   INTERFACE
      SUBROUTINE DORMRZ(Side,Trans,M,N,K,L,A,Lda,Tau,C,Ldc,Work,Lwork,  &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORMRZ
   END INTERFACE
END MODULE S_DORMRZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DORMTR
   INTERFACE
      SUBROUTINE DORMTR(Side,Uplo,Trans,M,N,A,Lda,Tau,C,Ldc,Work,Lwork, &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DORMTR
   END INTERFACE
END MODULE S_DORMTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPBCON
   INTERFACE
      SUBROUTINE DPBCON(Uplo,N,Kd,Ab,Ldab,Anorm,Rcond,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPBCON
   END INTERFACE
END MODULE S_DPBCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPBEQU
   INTERFACE
      SUBROUTINE DPBEQU(Uplo,N,Kd,Ab,Ldab,S,Scond,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) :: Scond
      REAL(R8KIND) , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPBEQU
   END INTERFACE
END MODULE S_DPBEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPBRFS
   INTERFACE
      SUBROUTINE DPBRFS(Uplo,N,Kd,Nrhs,Ab,Ldab,Afb,Ldafb,B,Ldb,X,Ldx,   &
     &                  Ferr,Berr,Work,Iwork,Info)
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
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPBRFS
   END INTERFACE
END MODULE S_DPBRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPBSTF
   INTERFACE
      SUBROUTINE DPBSTF(Uplo,N,Kd,Ab,Ldab,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPBSTF
   END INTERFACE
END MODULE S_DPBSTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPBSV
   INTERFACE
      SUBROUTINE DPBSV(Uplo,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPBSV
   END INTERFACE
END MODULE S_DPBSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPBSVX
   INTERFACE
      SUBROUTINE DPBSVX(Fact,Uplo,N,Kd,Nrhs,Ab,Ldab,Afb,Ldafb,Equed,S,B,&
     &                  Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPBSVX
   END INTERFACE
END MODULE S_DPBSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPBTF2
   INTERFACE
      SUBROUTINE DPBTF2(Uplo,N,Kd,Ab,Ldab,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kd
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPBTF2
   END INTERFACE
END MODULE S_DPBTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPBTRF
   INTERFACE
      SUBROUTINE DPBTRF(Uplo,N,Kd,Ab,Ldab,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER , PARAMETER  ::  NBMAX = 32 , LDWORK = NBMAX + 1
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPBTRF
   END INTERFACE
END MODULE S_DPBTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPBTRS
   INTERFACE
      SUBROUTINE DPBTRS(Uplo,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPBTRS
   END INTERFACE
END MODULE S_DPBTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPFTRF
   INTERFACE
      SUBROUTINE DPFTRF(Transr,Uplo,N,A,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPFTRF
   END INTERFACE
END MODULE S_DPFTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPFTRI
   INTERFACE
      SUBROUTINE DPFTRI(Transr,Uplo,N,A,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPFTRI
   END INTERFACE
END MODULE S_DPFTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPFTRS
   INTERFACE
      SUBROUTINE DPFTRS(Transr,Uplo,N,Nrhs,A,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(0:*) :: A
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPFTRS
   END INTERFACE
END MODULE S_DPFTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPOCON
   INTERFACE
      SUBROUTINE DPOCON(Uplo,N,A,Lda,Anorm,Rcond,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPOCON
   END INTERFACE
END MODULE S_DPOCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPOEQUB
   INTERFACE
      SUBROUTINE DPOEQUB(N,A,Lda,S,Scond,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) :: Scond
      REAL(R8KIND) , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPOEQUB
   END INTERFACE
END MODULE S_DPOEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPOEQU
   INTERFACE
      SUBROUTINE DPOEQU(N,A,Lda,S,Scond,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) :: Scond
      REAL(R8KIND) , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPOEQU
   END INTERFACE
END MODULE S_DPOEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPORFS
   INTERFACE
      SUBROUTINE DPORFS(Uplo,N,Nrhs,A,Lda,Af,Ldaf,B,Ldb,X,Ldx,Ferr,Berr,&
     &                  Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPORFS
   END INTERFACE
END MODULE S_DPORFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPORFSX
   INTERFACE
      SUBROUTINE DPORFSX(Uplo,Equed,N,Nrhs,A,Lda,Af,Ldaf,S,B,Ldb,X,Ldx, &
     &                   Rcond,Berr,N_err_bnds,Err_bnds_norm,           &
     &                   Err_bnds_comp,Nparams,Params,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ITREF_DEFAULT = 1.0D+0 ,            &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Params
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPORFSX
   END INTERFACE
END MODULE S_DPORFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPOSV
   INTERFACE
      SUBROUTINE DPOSV(Uplo,N,Nrhs,A,Lda,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPOSV
   END INTERFACE
END MODULE S_DPOSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPOSVX
   INTERFACE
      SUBROUTINE DPOSVX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Equed,S,B,Ldb,X, &
     &                  Ldx,Rcond,Ferr,Berr,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPOSVX
   END INTERFACE
END MODULE S_DPOSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPOSVXX
   INTERFACE
      SUBROUTINE DPOSVXX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Equed,S,B,Ldb,X,&
     &                   Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,Err_bnds_norm,&
     &                   Err_bnds_comp,Nparams,Params,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , INTENT(OUT) :: Rpvgrw
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL(R8KIND) , DIMENSION(*) :: Params
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPOSVXX
   END INTERFACE
END MODULE S_DPOSVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPOTF2
   INTERFACE
      SUBROUTINE DPOTF2(Uplo,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPOTF2
   END INTERFACE
END MODULE S_DPOTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPOTRF2
   INTERFACE
      RECURSIVE SUBROUTINE DPOTRF2(Uplo,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPOTRF2
   END INTERFACE
END MODULE S_DPOTRF2
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
MODULE S_DPOTRI
   INTERFACE
      SUBROUTINE DPOTRI(Uplo,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPOTRI
   END INTERFACE
END MODULE S_DPOTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPOTRS
   INTERFACE
      SUBROUTINE DPOTRS(Uplo,N,Nrhs,A,Lda,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPOTRS
   END INTERFACE
END MODULE S_DPOTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPPCON
   INTERFACE
      SUBROUTINE DPPCON(Uplo,N,Ap,Anorm,Rcond,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPPCON
   END INTERFACE
END MODULE S_DPPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPPEQU
   INTERFACE
      SUBROUTINE DPPEQU(Uplo,N,Ap,S,Scond,Amax,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) :: Scond
      REAL(R8KIND) , INTENT(INOUT) :: Amax
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPPEQU
   END INTERFACE
END MODULE S_DPPEQU
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPPRFS
   INTERFACE
      SUBROUTINE DPPRFS(Uplo,N,Nrhs,Ap,Afp,B,Ldb,X,Ldx,Ferr,Berr,Work,  &
     &                  Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: Afp
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPPRFS
   END INTERFACE
END MODULE S_DPPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPPSV
   INTERFACE
      SUBROUTINE DPPSV(Uplo,N,Nrhs,Ap,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPPSV
   END INTERFACE
END MODULE S_DPPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPPSVX
   INTERFACE
      SUBROUTINE DPPSVX(Fact,Uplo,N,Nrhs,Ap,Afp,Equed,S,B,Ldb,X,Ldx,    &
     &                  Rcond,Ferr,Berr,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: Afp
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPPSVX
   END INTERFACE
END MODULE S_DPPSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPPTRF
   INTERFACE
      SUBROUTINE DPPTRF(Uplo,N,Ap,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPPTRF
   END INTERFACE
END MODULE S_DPPTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPPTRI
   INTERFACE
      SUBROUTINE DPPTRI(Uplo,N,Ap,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPPTRI
   END INTERFACE
END MODULE S_DPPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPPTRS
   INTERFACE
      SUBROUTINE DPPTRS(Uplo,N,Nrhs,Ap,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPPTRS
   END INTERFACE
END MODULE S_DPPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPSTF2
   INTERFACE
      SUBROUTINE DPSTF2(Uplo,N,A,Lda,Piv,Rank,Tol,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(N) :: Piv
      INTEGER , INTENT(OUT) :: Rank
      REAL(R8KIND) , INTENT(IN) :: Tol
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(2*N) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPSTF2
   END INTERFACE
END MODULE S_DPSTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPSTRF
   INTERFACE
      SUBROUTINE DPSTRF(Uplo,N,A,Lda,Piv,Rank,Tol,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(N) :: Piv
      INTEGER :: Rank
      REAL(R8KIND) :: Tol
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(2*N) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPSTRF
   END INTERFACE
END MODULE S_DPSTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPTCON
   INTERFACE
      SUBROUTINE DPTCON(N,D,E,Anorm,Rcond,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: E
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPTCON
   END INTERFACE
END MODULE S_DPTCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPTEQR
   INTERFACE
      SUBROUTINE DPTEQR(Compz,N,D,E,Z,Ldz,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Compz
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPTEQR
   END INTERFACE
END MODULE S_DPTEQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPTRFS
   INTERFACE
      SUBROUTINE DPTRFS(N,Nrhs,D,E,Df,Ef,B,Ldb,X,Ldx,Ferr,Berr,Work,    &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  ITMAX = 5
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0 , THREE = 3.0D+0
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Df
      REAL(R8KIND) , DIMENSION(*) :: Ef
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPTRFS
   END INTERFACE
END MODULE S_DPTRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPTSV
   INTERFACE
      SUBROUTINE DPTSV(N,Nrhs,D,E,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPTSV
   END INTERFACE
END MODULE S_DPTSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPTSVX
   INTERFACE
      SUBROUTINE DPTSVX(Fact,N,Nrhs,D,E,Df,Ef,B,Ldb,X,Ldx,Rcond,Ferr,   &
     &                  Berr,Work,Info)
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
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Df
      REAL(R8KIND) , DIMENSION(*) :: Ef
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPTSVX
   END INTERFACE
END MODULE S_DPTSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPTTRF
   INTERFACE
      SUBROUTINE DPTTRF(N,D,E,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER :: Info
      END SUBROUTINE DPTTRF
   END INTERFACE
END MODULE S_DPTTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPTTRS
   INTERFACE
      SUBROUTINE DPTTRS(N,Nrhs,D,E,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DPTTRS
   END INTERFACE
END MODULE S_DPTTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DPTTS2
   INTERFACE
      SUBROUTINE DPTTS2(N,Nrhs,D,E,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: E
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      END SUBROUTINE DPTTS2
   END INTERFACE
END MODULE S_DPTTS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DRSCL
   INTERFACE
      SUBROUTINE DRSCL(N,Sa,Sx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) :: Sa
      REAL(R8KIND) , DIMENSION(*) :: Sx
      INTEGER :: Incx
      END SUBROUTINE DRSCL
   END INTERFACE
END MODULE S_DRSCL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSB2ST_KERNELS
   INTERFACE
      SUBROUTINE DSB2ST_KERNELS(Uplo,Wantz,Ttype,St,Ed,Sweep,N,Nb,Ib,A, &
     &                          Lda,V,Tau,Ldvt,Work)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: Ttype
      INTEGER , INTENT(IN) :: St
      INTEGER , INTENT(IN) :: Ed
      INTEGER , INTENT(IN) :: Sweep
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(IN) :: Ib
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , DIMENSION(*) :: V
      REAL(R8KIND) , DIMENSION(*) :: Tau
      INTEGER , INTENT(IN) :: Ldvt
      REAL(R8KIND) , DIMENSION(*) :: Work
      END SUBROUTINE DSB2ST_KERNELS
   END INTERFACE
END MODULE S_DSB2ST_KERNELS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSBEV_2STAGE
   INTERFACE
      SUBROUTINE DSBEV_2STAGE(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Lwork,&
     &                        Info)
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
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSBEV_2STAGE
   END INTERFACE
END MODULE S_DSBEV_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSBEVD_2STAGE
   INTERFACE
      SUBROUTINE DSBEVD_2STAGE(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,     &
     &                         Lwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSBEVD_2STAGE
   END INTERFACE
END MODULE S_DSBEVD_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSBEVD
   INTERFACE
      SUBROUTINE DSBEVD(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Lwork,Iwork,&
     &                  Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSBEVD
   END INTERFACE
END MODULE S_DSBEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSBEV
   INTERFACE
      SUBROUTINE DSBEV(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Info)
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
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSBEV
   END INTERFACE
END MODULE S_DSBEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSBEVX_2STAGE
   INTERFACE
      SUBROUTINE DSBEVX_2STAGE(Jobz,Range,Uplo,N,Kd,Ab,Ldab,Q,Ldq,Vl,Vu,&
     &                         Il,Iu,Abstol,M,W,Z,Ldz,Work,Lwork,Iwork, &
     &                         Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSBEVX_2STAGE
   END INTERFACE
END MODULE S_DSBEVX_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSBEVX
   INTERFACE
      SUBROUTINE DSBEVX(Jobz,Range,Uplo,N,Kd,Ab,Ldab,Q,Ldq,Vl,Vu,Il,Iu, &
     &                  Abstol,M,W,Z,Ldz,Work,Iwork,Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSBEVX
   END INTERFACE
END MODULE S_DSBEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSBGST
   INTERFACE
      SUBROUTINE DSBGST(Vect,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,X,Ldx,Work,   &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ka
      INTEGER , INTENT(IN) :: Kb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , DIMENSION(Ldbb,*) :: Bb
      INTEGER , INTENT(IN) :: Ldbb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSBGST
   END INTERFACE
END MODULE S_DSBGST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSBGVD
   INTERFACE
      SUBROUTINE DSBGVD(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,W,Z,Ldz,Work, &
     &                  Lwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ka
      INTEGER :: Kb
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldbb,*) :: Bb
      INTEGER :: Ldbb
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSBGVD
   END INTERFACE
END MODULE S_DSBGVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSBGV
   INTERFACE
      SUBROUTINE DSBGV(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,W,Z,Ldz,Work,  &
     &                 Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ka
      INTEGER :: Kb
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldbb,*) :: Bb
      INTEGER :: Ldbb
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSBGV
   END INTERFACE
END MODULE S_DSBGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSBGVX
   INTERFACE
      SUBROUTINE DSBGVX(Jobz,Range,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,Q,Ldq,  &
     &                  Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,Work,Iwork,Ifail,  &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Ka
      INTEGER :: Kb
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldbb,*) :: Bb
      INTEGER :: Ldbb
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSBGVX
   END INTERFACE
END MODULE S_DSBGVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSBTRD
   INTERFACE
      SUBROUTINE DSBTRD(Vect,Uplo,N,Kd,Ab,Ldab,D,E,Q,Ldq,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kd
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSBTRD
   END INTERFACE
END MODULE S_DSBTRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSFRK
   INTERFACE
      SUBROUTINE DSFRK(Transr,Uplo,Trans,N,K,Alpha,A,Lda,Beta,C)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Transr
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) :: Alpha
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) :: Beta
      REAL(R8KIND) , DIMENSION(*) :: C
      END SUBROUTINE DSFRK
   END INTERFACE
END MODULE S_DSFRK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSGESV
   INTERFACE
      SUBROUTINE DSGESV(N,Nrhs,A,Lda,Ipiv,B,Ldb,X,Ldx,Work,Swork,Iter,  &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      LOGICAL , PARAMETER  ::  DOITREF = .TRUE.
      INTEGER , PARAMETER  ::  ITERMAX = 30
      REAL(R8KIND) , PARAMETER  ::  BWDMAX = 1.0E+00 ,                  &
     &                              NEGONE = -1.0D+0 , ONE = 1.0D+0
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , DIMENSION(N,*) :: Work
      REAL , DIMENSION(*) :: Swork
      INTEGER , INTENT(OUT) :: Iter
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSGESV
   END INTERFACE
END MODULE S_DSGESV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPCON
   INTERFACE
      SUBROUTINE DSPCON(Uplo,N,Ap,Ipiv,Anorm,Rcond,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: Ap
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPCON
   END INTERFACE
END MODULE S_DSPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPEVD
   INTERFACE
      SUBROUTINE DSPEVD(Jobz,Uplo,N,Ap,W,Z,Ldz,Work,Lwork,Iwork,Liwork, &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPEVD
   END INTERFACE
END MODULE S_DSPEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPEV
   INTERFACE
      SUBROUTINE DSPEV(Jobz,Uplo,N,Ap,W,Z,Ldz,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPEV
   END INTERFACE
END MODULE S_DSPEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPEVX
   INTERFACE
      SUBROUTINE DSPEVX(Jobz,Range,Uplo,N,Ap,Vl,Vu,Il,Iu,Abstol,M,W,Z,  &
     &                  Ldz,Work,Iwork,Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPEVX
   END INTERFACE
END MODULE S_DSPEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPGST
   INTERFACE
      SUBROUTINE DSPGST(Itype,Uplo,N,Ap,Bp,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , HALF = 0.5D0
      INTEGER , INTENT(IN) :: Itype
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: Bp
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPGST
   END INTERFACE
END MODULE S_DSPGST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPGVD
   INTERFACE
      SUBROUTINE DSPGVD(Itype,Jobz,Uplo,N,Ap,Bp,W,Z,Ldz,Work,Lwork,     &
     &                  Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: Bp
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPGVD
   END INTERFACE
END MODULE S_DSPGVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPGV
   INTERFACE
      SUBROUTINE DSPGV(Itype,Jobz,Uplo,N,Ap,Bp,W,Z,Ldz,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: Bp
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPGV
   END INTERFACE
END MODULE S_DSPGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPGVX
   INTERFACE
      SUBROUTINE DSPGVX(Itype,Jobz,Range,Uplo,N,Ap,Bp,Vl,Vu,Il,Iu,      &
     &                  Abstol,M,W,Z,Ldz,Work,Iwork,Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: Bp
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPGVX
   END INTERFACE
END MODULE S_DSPGVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPOSV
   INTERFACE
      SUBROUTINE DSPOSV(Uplo,N,Nrhs,A,Lda,B,Ldb,X,Ldx,Work,Swork,Iter,  &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      LOGICAL , PARAMETER  ::  DOITREF = .TRUE.
      INTEGER , PARAMETER  ::  ITERMAX = 30
      REAL(R8KIND) , PARAMETER  ::  BWDMAX = 1.0E+00 ,                  &
     &                              NEGONE = -1.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , DIMENSION(N,*) :: Work
      REAL , DIMENSION(*) :: Swork
      INTEGER , INTENT(OUT) :: Iter
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPOSV
   END INTERFACE
END MODULE S_DSPOSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPRFS
   INTERFACE
      SUBROUTINE DSPRFS(Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Ferr,Berr,  &
     &                  Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPRFS
   END INTERFACE
END MODULE S_DSPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPSV
   INTERFACE
      SUBROUTINE DSPSV(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: Ap
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPSV
   END INTERFACE
END MODULE S_DSPSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPSVX
   INTERFACE
      SUBROUTINE DSPSVX(Fact,Uplo,N,Nrhs,Ap,Afp,Ipiv,B,Ldb,X,Ldx,Rcond, &
     &                  Ferr,Berr,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: Afp
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPSVX
   END INTERFACE
END MODULE S_DSPSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPTRD
   INTERFACE
      SUBROUTINE DSPTRD(Uplo,N,Ap,D,E,Tau,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0 ,        &
     &                              HALF = 1.0D0/2.0D0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPTRD
   END INTERFACE
END MODULE S_DSPTRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPTRF
   INTERFACE
      SUBROUTINE DSPTRF(Uplo,N,Ap,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPTRF
   END INTERFACE
END MODULE S_DSPTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPTRI
   INTERFACE
      SUBROUTINE DSPTRI(Uplo,N,Ap,Ipiv,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPTRI
   END INTERFACE
END MODULE S_DSPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSPTRS
   INTERFACE
      SUBROUTINE DSPTRS(Uplo,N,Nrhs,Ap,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: Ap
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSPTRS
   END INTERFACE
END MODULE S_DSPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSTEBZ
   INTERFACE
      SUBROUTINE DSTEBZ(Range,Order,N,Vl,Vu,Il,Iu,Abstol,D,E,M,Nsplit,W,&
     &                  Iblock,Isplit,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , HALF = 1.0D0/TWO ,    &
     &                              FUDGE = 2.1D0 , RELFAC = 2.0D0
      CHARACTER :: Range
      CHARACTER :: Order
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) :: M
      INTEGER , INTENT(INOUT) :: Nsplit
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isplit
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSTEBZ
   END INTERFACE
END MODULE S_DSTEBZ
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSTEDC
   INTERFACE
      SUBROUTINE DSTEDC(Compz,N,D,E,Z,Ldz,Work,Lwork,Iwork,Liwork,Info)
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
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSTEDC
   END INTERFACE
END MODULE S_DSTEDC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSTEGR
   INTERFACE
      SUBROUTINE DSTEGR(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,  &
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
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER :: Info
      END SUBROUTINE DSTEGR
   END INTERFACE
END MODULE S_DSTEGR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSTEIN
   INTERFACE
      SUBROUTINE DSTEIN(N,D,E,M,W,Iblock,Isplit,Z,Ldz,Work,Iwork,Ifail, &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
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
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER , INTENT(IN) :: Ldz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSTEIN
   END INTERFACE
END MODULE S_DSTEIN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSTEMR
   INTERFACE
      SUBROUTINE DSTEMR(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,M,W,Z,Ldz,Nzc,     &
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
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(IN) :: Nzc
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Isuppz
      LOGICAL , INTENT(INOUT) :: Tryrac
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSTEMR
   END INTERFACE
END MODULE S_DSTEMR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSTEQR
   INTERFACE
      SUBROUTINE DSTEQR(Compz,N,D,E,Z,Ldz,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , THREE = 3.0D0
      INTEGER , PARAMETER  ::  MAXIT = 30
      CHARACTER :: Compz
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSTEQR
   END INTERFACE
END MODULE S_DSTEQR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSTERF
   INTERFACE
      SUBROUTINE DSTERF(N,D,E,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , THREE = 3.0D0
      INTEGER , PARAMETER  ::  MAXIT = 30
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSTERF
   END INTERFACE
END MODULE S_DSTERF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSTEVD
   INTERFACE
      SUBROUTINE DSTEVD(Jobz,N,D,E,Z,Ldz,Work,Lwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSTEVD
   END INTERFACE
END MODULE S_DSTEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSTEV
   INTERFACE
      SUBROUTINE DSTEV(Jobz,N,D,E,Z,Ldz,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSTEV
   END INTERFACE
END MODULE S_DSTEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSTEVR
   INTERFACE
      SUBROUTINE DSTEVR(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,  &
     &                  Isuppz,Work,Lwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
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
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSTEVR
   END INTERFACE
END MODULE S_DSTEVR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSTEVX
   INTERFACE
      SUBROUTINE DSTEVX(Jobz,Range,N,D,E,Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,  &
     &                  Work,Iwork,Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      CHARACTER :: Range
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSTEVX
   END INTERFACE
END MODULE S_DSTEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYCON_3
   INTERFACE
      SUBROUTINE DSYCON_3(Uplo,N,A,Lda,E,Ipiv,Anorm,Rcond,Work,Iwork,   &
     &                    Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYCON_3
   END INTERFACE
END MODULE S_DSYCON_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYCON
   INTERFACE
      SUBROUTINE DSYCON(Uplo,N,A,Lda,Ipiv,Anorm,Rcond,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYCON
   END INTERFACE
END MODULE S_DSYCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYCON_ROOK
   INTERFACE
      SUBROUTINE DSYCON_ROOK(Uplo,N,A,Lda,Ipiv,Anorm,Rcond,Work,Iwork,  &
     &                       Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYCON_ROOK
   END INTERFACE
END MODULE S_DSYCON_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYCONV
   INTERFACE
      SUBROUTINE DSYCONV(Uplo,Way,N,A,Lda,Ipiv,E,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYCONV
   END INTERFACE
END MODULE S_DSYCONV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYCONVF
   INTERFACE
      SUBROUTINE DSYCONVF(Uplo,Way,N,A,Lda,E,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYCONVF
   END INTERFACE
END MODULE S_DSYCONVF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYCONVF_ROOK
   INTERFACE
      SUBROUTINE DSYCONVF_ROOK(Uplo,Way,N,A,Lda,E,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYCONVF_ROOK
   END INTERFACE
END MODULE S_DSYCONVF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYEQUB
   INTERFACE
      SUBROUTINE DSYEQUB(Uplo,N,A,Lda,S,Scond,Amax,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0
      INTEGER , PARAMETER  ::  MAX_ITER = 100
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) :: Scond
      REAL(R8KIND) , INTENT(INOUT) :: Amax
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYEQUB
   END INTERFACE
END MODULE S_DSYEQUB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYEV_2STAGE
   INTERFACE
      SUBROUTINE DSYEV_2STAGE(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYEV_2STAGE
   END INTERFACE
END MODULE S_DSYEV_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYEVD_2STAGE
   INTERFACE
      SUBROUTINE DSYEVD_2STAGE(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Iwork,    &
     &                         Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYEVD_2STAGE
   END INTERFACE
END MODULE S_DSYEVD_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYEVD
   INTERFACE
      SUBROUTINE DSYEVD(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Iwork,Liwork,    &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYEVD
   END INTERFACE
END MODULE S_DSYEVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYEV
   INTERFACE
      SUBROUTINE DSYEV(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYEV
   END INTERFACE
END MODULE S_DSYEV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYEVR_2STAGE
   INTERFACE
      SUBROUTINE DSYEVR_2STAGE(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,     &
     &                         Abstol,M,W,Z,Ldz,Isuppz,Work,Lwork,Iwork,&
     &                         Liwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYEVR_2STAGE
   END INTERFACE
END MODULE S_DSYEVR_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYEVR
   INTERFACE
      SUBROUTINE DSYEVR(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,Abstol,M,W, &
     &                  Z,Ldz,Isuppz,Work,Lwork,Iwork,Liwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , DIMENSION(*) :: Isuppz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYEVR
   END INTERFACE
END MODULE S_DSYEVR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYEVX_2STAGE
   INTERFACE
      SUBROUTINE DSYEVX_2STAGE(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,     &
     &                         Abstol,M,W,Z,Ldz,Work,Lwork,Iwork,Ifail, &
     &                         Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYEVX_2STAGE
   END INTERFACE
END MODULE S_DSYEVX_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYEVX
   INTERFACE
      SUBROUTINE DSYEVX(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,Abstol,M,W, &
     &                  Z,Ldz,Work,Lwork,Iwork,Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(IN) :: Abstol
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYEVX
   END INTERFACE
END MODULE S_DSYEVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYGS2
   INTERFACE
      SUBROUTINE DSYGS2(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , HALF = 0.5D0
      INTEGER , INTENT(IN) :: Itype
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYGS2
   END INTERFACE
END MODULE S_DSYGS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYGST
   INTERFACE
      SUBROUTINE DSYGST(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , HALF = 0.5D0
      INTEGER :: Itype
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYGST
   END INTERFACE
END MODULE S_DSYGST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYGV_2STAGE
   INTERFACE
      SUBROUTINE DSYGV_2STAGE(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,     &
     &                        Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYGV_2STAGE
   END INTERFACE
END MODULE S_DSYGV_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYGVD
   INTERFACE
      SUBROUTINE DSYGVD(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,Lwork,     &
     &                  Iwork,Liwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYGVD
   END INTERFACE
END MODULE S_DSYGVD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYGV
   INTERFACE
      SUBROUTINE DSYGV(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYGV
   END INTERFACE
END MODULE S_DSYGV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYGVX
   INTERFACE
      SUBROUTINE DSYGVX(Itype,Jobz,Range,Uplo,N,A,Lda,B,Ldb,Vl,Vu,Il,Iu,&
     &                  Abstol,M,W,Z,Ldz,Work,Lwork,Iwork,Ifail,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) :: Vl
      REAL(R8KIND) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) :: Abstol
      INTEGER :: M
      REAL(R8KIND) , DIMENSION(*) :: W
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYGVX
   END INTERFACE
END MODULE S_DSYGVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYRFS
   INTERFACE
      SUBROUTINE DSYRFS(Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,Ferr,&
     &                  Berr,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Berr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYRFS
   END INTERFACE
END MODULE S_DSYRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYRFSX
   INTERFACE
      SUBROUTINE DSYRFSX(Uplo,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,S,B,Ldb,X,&
     &                   Ldx,Rcond,Berr,N_err_bnds,Err_bnds_norm,       &
     &                   Err_bnds_comp,Nparams,Params,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ITREF_DEFAULT = 1.0D+0 ,            &
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Params
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYRFSX
   END INTERFACE
END MODULE S_DSYRFSX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYSV_AA_2STAGE
   INTERFACE
      SUBROUTINE DSYSV_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B, &
     &                           Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tb
      INTEGER :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYSV_AA_2STAGE
   END INTERFACE
END MODULE S_DSYSV_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYSV_AA
   INTERFACE
      SUBROUTINE DSYSV_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYSV_AA
   END INTERFACE
END MODULE S_DSYSV_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYSV
   INTERFACE
      SUBROUTINE DSYSV(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYSV
   END INTERFACE
END MODULE S_DSYSV
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYSV_RK
   INTERFACE
      SUBROUTINE DSYSV_RK(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Work,Lwork,    &
     &                    Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYSV_RK
   END INTERFACE
END MODULE S_DSYSV_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYSV_ROOK
   INTERFACE
      SUBROUTINE DSYSV_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,    &
     &                      Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYSV_ROOK
   END INTERFACE
END MODULE S_DSYSV_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYSVX
   INTERFACE
      SUBROUTINE DSYSVX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,&
     &                  Rcond,Ferr,Berr,Work,Lwork,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , DIMENSION(*) :: Berr
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYSVX
   END INTERFACE
END MODULE S_DSYSVX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYSVXX
   INTERFACE
      SUBROUTINE DSYSVXX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,S,B, &
     &                   Ldb,X,Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,        &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , DIMENSION(*) :: S
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , INTENT(OUT) :: Rpvgrw
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL(R8KIND) , DIMENSION(*) :: Params
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYSVXX
   END INTERFACE
END MODULE S_DSYSVXX
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYSWAPR
   INTERFACE
      SUBROUTINE DSYSWAPR(Uplo,N,A,Lda,I1,I2)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,N) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) :: I1
      INTEGER , INTENT(IN) :: I2
      END SUBROUTINE DSYSWAPR
   END INTERFACE
END MODULE S_DSYSWAPR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTD2
   INTERFACE
      SUBROUTINE DSYTD2(Uplo,N,A,Lda,D,E,Tau,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , ZERO = 0.0D0 ,        &
     &                              HALF = 1.0D0/2.0D0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Tau
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTD2
   END INTERFACE
END MODULE S_DSYTD2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTF2
   INTERFACE
      SUBROUTINE DSYTF2(Uplo,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTF2
   END INTERFACE
END MODULE S_DSYTF2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTF2_RK
   INTERFACE
      SUBROUTINE DSYTF2_RK(Uplo,N,A,Lda,E,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTF2_RK
   END INTERFACE
END MODULE S_DSYTF2_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTF2_ROOK
   INTERFACE
      SUBROUTINE DSYTF2_ROOK(Uplo,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              EIGHT = 8.0D+0 , SEVTEN = 17.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTF2_ROOK
   END INTERFACE
END MODULE S_DSYTF2_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRD_2STAGE
   INTERFACE
      SUBROUTINE DSYTRD_2STAGE(Vect,Uplo,N,A,Lda,D,E,Tau,Hous2,Lhous2,  &
     &                         Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Hous2
      INTEGER :: Lhous2
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRD_2STAGE
   END INTERFACE
END MODULE S_DSYTRD_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRD
   INTERFACE
      SUBROUTINE DSYTRD(Uplo,N,A,Lda,D,E,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRD
   END INTERFACE
END MODULE S_DSYTRD
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRD_SB2ST
   INTERFACE
      SUBROUTINE DSYTRD_SB2ST(Stage1,Vect,Uplo,N,Kd,Ab,Ldab,D,E,Hous,   &
     &                        Lhous,Work,Lwork,Info)
      USE OMP_LIB
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  RZERO = 0.0D+0 , ZERO = 0.0D+0 ,    &
     &                              ONE = 1.0D+0
      CHARACTER :: Stage1
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Hous
      INTEGER , INTENT(IN) :: Lhous
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRD_SB2ST
   END INTERFACE
END MODULE S_DSYTRD_SB2ST
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRD_SY2SB
   INTERFACE
      SUBROUTINE DSYTRD_SY2SB(Uplo,N,Kd,A,Lda,Ab,Ldab,Tau,Work,Lwork,   &
     &                        Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  RONE = 1.0D+0 , ZERO = 0.0D+0 ,     &
     &                              ONE = 1.0D+0 , HALF = 0.5D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRD_SY2SB
   END INTERFACE
END MODULE S_DSYTRD_SY2SB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRF_AA_2STAGE
   INTERFACE
      SUBROUTINE DSYTRF_AA_2STAGE(Uplo,N,A,Lda,Tb,Ltb,Ipiv,Ipiv2,Work,  &
     &                            Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRF_AA_2STAGE
   END INTERFACE
END MODULE S_DSYTRF_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRF_AA
   INTERFACE
      SUBROUTINE DSYTRF_AA(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRF_AA
   END INTERFACE
END MODULE S_DSYTRF_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRF
   INTERFACE
      SUBROUTINE DSYTRF(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRF
   END INTERFACE
END MODULE S_DSYTRF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRF_RK
   INTERFACE
      SUBROUTINE DSYTRF_RK(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRF_RK
   END INTERFACE
END MODULE S_DSYTRF_RK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRF_ROOK
   INTERFACE
      SUBROUTINE DSYTRF_ROOK(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRF_ROOK
   END INTERFACE
END MODULE S_DSYTRF_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRI2
   INTERFACE
      SUBROUTINE DSYTRI2(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRI2
   END INTERFACE
END MODULE S_DSYTRI2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRI2X
   INTERFACE
      SUBROUTINE DSYTRI2X(Uplo,N,A,Lda,Ipiv,Work,Nb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRI2X
   END INTERFACE
END MODULE S_DSYTRI2X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRI_3
   INTERFACE
      SUBROUTINE DSYTRI_3(Uplo,N,A,Lda,E,Ipiv,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: E
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRI_3
   END INTERFACE
END MODULE S_DSYTRI_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRI_3X
   INTERFACE
      SUBROUTINE DSYTRI_3X(Uplo,N,A,Lda,E,Ipiv,Work,Nb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(N+Nb+1,*) :: Work
      INTEGER , INTENT(IN) :: Nb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRI_3X
   END INTERFACE
END MODULE S_DSYTRI_3X
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRI
   INTERFACE
      SUBROUTINE DSYTRI(Uplo,N,A,Lda,Ipiv,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRI
   END INTERFACE
END MODULE S_DSYTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRI_ROOK
   INTERFACE
      SUBROUTINE DSYTRI_ROOK(Uplo,N,A,Lda,Ipiv,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRI_ROOK
   END INTERFACE
END MODULE S_DSYTRI_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRS2
   INTERFACE
      SUBROUTINE DSYTRS2(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRS2
   END INTERFACE
END MODULE S_DSYTRS2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRS_3
   INTERFACE
      SUBROUTINE DSYTRS_3(Uplo,N,Nrhs,A,Lda,E,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: E
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRS_3
   END INTERFACE
END MODULE S_DSYTRS_3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRS_AA_2STAGE
   INTERFACE
      SUBROUTINE DSYTRS_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B,&
     &                            Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRS_AA_2STAGE
   END INTERFACE
END MODULE S_DSYTRS_AA_2STAGE
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRS_AA
   INTERFACE
      SUBROUTINE DSYTRS_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRS_AA
   END INTERFACE
END MODULE S_DSYTRS_AA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRS
   INTERFACE
      SUBROUTINE DSYTRS(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRS
   END INTERFACE
END MODULE S_DSYTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DSYTRS_ROOK
   INTERFACE
      SUBROUTINE DSYTRS_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DSYTRS_ROOK
   END INTERFACE
END MODULE S_DSYTRS_ROOK
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTBCON
   INTERFACE
      SUBROUTINE DTBCON(Norm,Uplo,Diag,N,Kd,Ab,Ldab,Rcond,Work,Iwork,   &
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
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTBCON
   END INTERFACE
END MODULE S_DTBCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTBRFS
   INTERFACE
      SUBROUTINE DTBRFS(Uplo,Trans,Diag,N,Kd,Nrhs,Ab,Ldab,B,Ldb,X,Ldx,  &
     &                  Ferr,Berr,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Berr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTBRFS
   END INTERFACE
END MODULE S_DTBRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTBTRS
   INTERFACE
      SUBROUTINE DTBTRS(Uplo,Trans,Diag,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTBTRS
   END INTERFACE
END MODULE S_DTBTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTFSM
   INTERFACE
      SUBROUTINE DTFSM(Transr,Side,Uplo,Trans,Diag,M,N,Alpha,A,B,Ldb)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Transr
      CHARACTER :: Side
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) :: Alpha
      REAL(R8KIND) , DIMENSION(0:*) :: A
      REAL(R8KIND) , DIMENSION(0:Ldb-1,0:*) :: B
      INTEGER :: Ldb
      END SUBROUTINE DTFSM
   END INTERFACE
END MODULE S_DTFSM
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTFTRI
   INTERFACE
      SUBROUTINE DTFTRI(Transr,Uplo,Diag,N,A,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Transr
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTFTRI
   END INTERFACE
END MODULE S_DTFTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTFTTP
   INTERFACE
      SUBROUTINE DTFTTP(Transr,Uplo,N,Arf,Ap,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(0:*) :: Arf
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(0:*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTFTTP
   END INTERFACE
END MODULE S_DTFTTP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTFTTR
   INTERFACE
      SUBROUTINE DTFTTR(Transr,Uplo,N,Arf,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(0:*) :: Arf
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(0:Lda-1,0:*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTFTTR
   END INTERFACE
END MODULE S_DTFTTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTGEVC
   INTERFACE
      SUBROUTINE DTGEVC(Side,Howmny,Select,N,S,Lds,P,Ldp,Vl,Ldvl,Vr,    &
     &                  Ldvr,Mm,M,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              SAFETY = 1.0D+2
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lds,*) :: S
      INTEGER :: Lds
      REAL(R8KIND) , DIMENSION(Ldp,*) :: P
      INTEGER :: Ldp
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(OUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTGEVC
   END INTERFACE
END MODULE S_DTGEVC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTGEX2
   INTERFACE
      SUBROUTINE DTGEX2(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,J1,N1,N2, &
     &                  Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWENTY = 2.0D+01
      INTEGER , PARAMETER  ::  LDST = 4
      LOGICAL , PARAMETER  ::  WANDS = .TRUE.
      LOGICAL , INTENT(IN) :: Wantq
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(IN) :: J1
      INTEGER :: N1
      INTEGER :: N2
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER :: Info
      END SUBROUTINE DTGEX2
   END INTERFACE
END MODULE S_DTGEX2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTGEXC
   INTERFACE
      SUBROUTINE DTGEXC(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,Ifst,Ilst,&
     &                  Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      LOGICAL :: Wantq
      LOGICAL :: Wantz
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ifst
      INTEGER , INTENT(INOUT) :: Ilst
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTGEXC
   END INTERFACE
END MODULE S_DTGEXC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTGSEN
   INTERFACE
      SUBROUTINE DTGSEN(Ijob,Wantq,Wantz,Select,N,A,Lda,B,Ldb,Alphar,   &
     &                  Alphai,Beta,Q,Ldq,Z,Ldz,M,Pl,Pr,Dif,Work,Lwork, &
     &                  Iwork,Liwork,Info)
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
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , DIMENSION(*) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) :: Pl
      REAL(R8KIND) , INTENT(INOUT) :: Pr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dif
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTGSEN
   END INTERFACE
END MODULE S_DTGSEN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTGSJA
   INTERFACE
      SUBROUTINE DTGSJA(Jobu,Jobv,Jobq,M,P,N,K,L,A,Lda,B,Ldb,Tola,Tolb, &
     &                  Alpha,Beta,U,Ldu,V,Ldv,Q,Ldq,Work,Ncycle,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  MAXIT = 40
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      CHARACTER :: Jobq
      INTEGER :: M
      INTEGER :: P
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER :: L
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(IN) :: Tola
      REAL(R8KIND) , INTENT(IN) :: Tolb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alpha
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL(R8KIND) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) :: Ncycle
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTGSJA
   END INTERFACE
END MODULE S_DTGSJA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTGSNA
   INTERFACE
      SUBROUTINE DTGSNA(Job,Howmny,Select,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,&
     &                  S,Dif,Mm,M,Work,Lwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  DIFDRI = 3
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0 , FOUR = 4.0D+0
      CHARACTER :: Job
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      REAL(R8KIND) , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dif
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTGSNA
   END INTERFACE
END MODULE S_DTGSNA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTGSY2
   INTERFACE
      SUBROUTINE DTGSY2(Trans,Ijob,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F, &
     &                  Ldf,Scale,Rdsum,Rdscal,Iwork,Pq,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER  ::  LDZ = 8
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Trans
      INTEGER :: Ijob
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(Ldd,*) :: D
      INTEGER :: Ldd
      REAL(R8KIND) , DIMENSION(Lde,*) :: E
      INTEGER :: Lde
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) :: Rdsum
      REAL(R8KIND) :: Rdscal
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) :: Pq
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTGSY2
   END INTERFACE
END MODULE S_DTGSY2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTGSYL
   INTERFACE
      SUBROUTINE DTGSYL(Trans,Ijob,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F, &
     &                  Ldf,Scale,Dif,Work,Lwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: Ijob
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(Ldd,*) :: D
      INTEGER :: Ldd
      REAL(R8KIND) , DIMENSION(Lde,*) :: E
      INTEGER :: Lde
      REAL(R8KIND) , DIMENSION(Ldf,*) :: F
      INTEGER :: Ldf
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(OUT) :: Dif
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTGSYL
   END INTERFACE
END MODULE S_DTGSYL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPCON
   INTERFACE
      SUBROUTINE DTPCON(Norm,Uplo,Diag,N,Ap,Rcond,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTPCON
   END INTERFACE
END MODULE S_DTPCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPLQT2
   INTERFACE
      SUBROUTINE DTPLQT2(M,N,L,A,Lda,B,Ldb,T,Ldt,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0 , ZERO = 0.0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: L
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTPLQT2
   END INTERFACE
END MODULE S_DTPLQT2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPLQT
   INTERFACE
      SUBROUTINE DTPLQT(M,N,L,Mb,A,Lda,B,Ldb,T,Ldt,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Mb
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTPLQT
   END INTERFACE
END MODULE S_DTPLQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPMLQT
   INTERFACE
      SUBROUTINE DTPMLQT(Side,Trans,M,N,K,L,Mb,V,Ldv,T,Ldt,A,Lda,B,Ldb, &
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
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTPMLQT
   END INTERFACE
END MODULE S_DTPMLQT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPMQRT
   INTERFACE
      SUBROUTINE DTPMQRT(Side,Trans,M,N,K,L,Nb,V,Ldv,T,Ldt,A,Lda,B,Ldb, &
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
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTPMQRT
   END INTERFACE
END MODULE S_DTPMQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPQRT2
   INTERFACE
      SUBROUTINE DTPQRT2(M,N,L,A,Lda,B,Ldb,T,Ldt,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0 , ZERO = 0.0
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER :: L
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTPQRT2
   END INTERFACE
END MODULE S_DTPQRT2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPQRT
   INTERFACE
      SUBROUTINE DTPQRT(M,N,L,Nb,A,Lda,B,Ldb,T,Ldt,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Nb
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTPQRT
   END INTERFACE
END MODULE S_DTPQRT
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPRFB
   INTERFACE
      SUBROUTINE DTPRFB(Side,Trans,Direct,Storev,M,N,K,L,V,Ldv,T,Ldt,A, &
     &                  Lda,B,Ldb,Work,Ldwork)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0 , ZERO = 0.0
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: L
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      END SUBROUTINE DTPRFB
   END INTERFACE
END MODULE S_DTPRFB
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPRFS
   INTERFACE
      SUBROUTINE DTPRFS(Uplo,Trans,Diag,N,Nrhs,Ap,B,Ldb,X,Ldx,Ferr,Berr,&
     &                  Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Berr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTPRFS
   END INTERFACE
END MODULE S_DTPRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPTRI
   INTERFACE
      SUBROUTINE DTPTRI(Uplo,Diag,N,Ap,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTPTRI
   END INTERFACE
END MODULE S_DTPTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPTRS
   INTERFACE
      SUBROUTINE DTPTRS(Uplo,Trans,Diag,N,Nrhs,Ap,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTPTRS
   END INTERFACE
END MODULE S_DTPTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPTTF
   INTERFACE
      SUBROUTINE DTPTTF(Transr,Uplo,N,Ap,Arf,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(0:*) :: Ap
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(0:*) :: Arf
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTPTTF
   END INTERFACE
END MODULE S_DTPTTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTPTTR
   INTERFACE
      SUBROUTINE DTPTTR(Uplo,N,Ap,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTPTTR
   END INTERFACE
END MODULE S_DTPTTR
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRCON
   INTERFACE
      SUBROUTINE DTRCON(Norm,Uplo,Diag,N,A,Lda,Rcond,Work,Iwork,Info)
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
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTRCON
   END INTERFACE
END MODULE S_DTRCON
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTREVC3
   INTERFACE
      SUBROUTINE DTREVC3(Side,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,Mm, &
     &                   M,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , PARAMETER  ::  NBMIN = 8 , NBMAX = 128
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(INOUT) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTREVC3
   END INTERFACE
END MODULE S_DTREVC3
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTREVC
   INTERFACE
      SUBROUTINE DTREVC(Side,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,Mm,M,&
     &                  Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(INOUT) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTREVC
   END INTERFACE
END MODULE S_DTREVC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTREXC
   INTERFACE
      SUBROUTINE DTREXC(Compq,N,T,Ldt,Q,Ldq,Ifst,Ilst,Work,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      CHARACTER :: Compq
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      INTEGER , INTENT(INOUT) :: Ifst
      INTEGER , INTENT(INOUT) :: Ilst
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTREXC
   END INTERFACE
END MODULE S_DTREXC
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRRFS
   INTERFACE
      SUBROUTINE DTRRFS(Uplo,Trans,Diag,N,Nrhs,A,Lda,B,Ldb,X,Ldx,Ferr,  &
     &                  Berr,Work,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER , INTENT(IN) :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ferr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Berr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTRRFS
   END INTERFACE
END MODULE S_DTRRFS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRSEN
   INTERFACE
      SUBROUTINE DTRSEN(Job,Compq,Select,N,T,Ldt,Q,Ldq,Wr,Wi,M,S,Sep,   &
     &                  Work,Lwork,Iwork,Liwork,Info)
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
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Wr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wi
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(OUT) :: S
      REAL(R8KIND) , INTENT(OUT) :: Sep
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTRSEN
   END INTERFACE
END MODULE S_DTRSEN
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRSNA
   INTERFACE
      SUBROUTINE DTRSNA(Job,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,S,Sep,&
     &                  Mm,M,Work,Ldwork,Iwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
      CHARACTER :: Job
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      REAL(R8KIND) , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Sep
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTRSNA
   END INTERFACE
END MODULE S_DTRSNA
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRSYL
   INTERFACE
      SUBROUTINE DTRSYL(Trana,Tranb,Isgn,M,N,A,Lda,B,Ldb,C,Ldc,Scale,   &
     &                  Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Trana
      CHARACTER :: Tranb
      INTEGER :: Isgn
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTRSYL
   END INTERFACE
END MODULE S_DTRSYL
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRTI2
   INTERFACE
      SUBROUTINE DTRTI2(Uplo,Diag,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTRTI2
   END INTERFACE
END MODULE S_DTRTI2
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRTRI
   INTERFACE
      SUBROUTINE DTRTRI(Uplo,Diag,N,A,Lda,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTRTRI
   END INTERFACE
END MODULE S_DTRTRI
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRTRS
   INTERFACE
      SUBROUTINE DTRTRS(Uplo,Trans,Diag,N,Nrhs,A,Lda,B,Ldb,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      CHARACTER :: Uplo
      CHARACTER :: Trans
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTRTRS
   END INTERFACE
END MODULE S_DTRTRS
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRTTF
   INTERFACE
      SUBROUTINE DTRTTF(Transr,Uplo,N,A,Lda,Arf,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(0:Lda-1,0:*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(0:*) :: Arf
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTRTTF
   END INTERFACE
END MODULE S_DTRTTF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTRTTP
   INTERFACE
      SUBROUTINE DTRTTP(Uplo,N,A,Lda,Ap,Info)
      USE F77KINDS                        
      IMPLICIT NONE
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTRTTP
   END INTERFACE
END MODULE S_DTRTTP
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DTZRZF
   INTERFACE
      SUBROUTINE DTZRZF(M,N,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!
! PARAMETER definitions
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
      INTEGER :: M
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
      END SUBROUTINE DTZRZF
   END INTERFACE
END MODULE S_DTZRZF
!*==intfaces.f90  created by SPAG 7.51RB at 21:36 on  3 Mar 2022
MODULE S_DZSUM1
   INTERFACE
      FUNCTION DZSUM1(N,Cx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
      REAL(R8KIND) :: DZSUM1
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      END FUNCTION DZSUM1
   END INTERFACE
END MODULE S_DZSUM1
