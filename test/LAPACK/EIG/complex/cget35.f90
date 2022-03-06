!*==cget35.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b CGET35
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGET35( RMAX, LMAX, NINFO, KNT, NIN )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NIN, NINFO
!       REAL               RMAX
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGET35 tests CTRSYL, a routine for solving the Sylvester matrix
!> equation
!>
!>    op(A)*X + ISGN*X*op(B) = scale*C,
!>
!> A and B are assumed to be in Schur canonical form, op() represents an
!> optional transpose, and ISGN can be -1 or +1.  Scale is an output
!> less than or equal to 1, chosen to avoid overflow in X.
!>
!> The test code verifies that the following residual is order 1:
!>
!>    norm(op(A)*X + ISGN*X*op(B) - scale*C) /
!>        (EPS*max(norm(A),norm(B))*norm(X))
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] RMAX
!> \verbatim
!>          RMAX is REAL
!>          Value of the largest test ratio.
!> \endverbatim
!>
!> \param[out] LMAX
!> \verbatim
!>          LMAX is INTEGER
!>          Example number where largest test ratio achieved.
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER
!>          Number of examples where INFO is nonzero.
!> \endverbatim
!>
!> \param[out] KNT
!> \verbatim
!>          KNT is INTEGER
!>          Total number of examples tested.
!> \endverbatim
!>
!> \param[in] NIN
!> \verbatim
!>          NIN is INTEGER
!>          Input logical unit number.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE CGET35(Rmax,Lmax,Ninfo,Knt,Nin)
      IMPLICIT NONE
!*--CGET3588
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Knt , Lmax , Nin , Ninfo
      REAL Rmax
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER LDT
      PARAMETER (LDT=10)
      REAL ZERO , ONE , TWO
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,TWO=2.0E0)
      REAL LARGE
      PARAMETER (LARGE=1.0E6)
      COMPLEX CONE
      PARAMETER (CONE=1.0E0)
!     ..
!     .. Local Scalars ..
      CHARACTER trana , tranb
      INTEGER i , imla , imlad , imlb , imlc , info , isgn , itrana ,   &
     &        itranb , j , m , n
      REAL bignum , eps , res , res1 , scale , smlnum , tnrm , xnrm
      COMPLEX rmul
!     ..
!     .. Local Arrays ..
      REAL dum(1) , vm1(3) , vm2(3)
      COMPLEX a(LDT,LDT) , atmp(LDT,LDT) , b(LDT,LDT) , btmp(LDT,LDT) , &
     &        c(LDT,LDT) , csav(LDT,LDT) , ctmp(LDT,LDT)
!     ..
!     .. External Functions ..
      REAL CLANGE , SLAMCH
      EXTERNAL CLANGE , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMM , CTRSYL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , REAL , SQRT
!     ..
!     .. Executable Statements ..
!
!     Get machine parameters
!
      eps = SLAMCH('P')
      smlnum = SLAMCH('S')/eps
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
!
!     Set up test case parameters
!
      vm1(1) = SQRT(smlnum)
      vm1(2) = ONE
      vm1(3) = LARGE
      vm2(1) = ONE
      vm2(2) = ONE + TWO*eps
      vm2(3) = TWO
!
      Knt = 0
      Ninfo = 0
      Lmax = 0
      Rmax = ZERO
      DO
!
!     Begin test loop
!
         READ (Nin,FMT=*) m , n
         IF ( n==0 ) RETURN
         DO i = 1 , m
            READ (Nin,FMT=*) (atmp(i,j),j=1,m)
         ENDDO
         DO i = 1 , n
            READ (Nin,FMT=*) (btmp(i,j),j=1,n)
         ENDDO
         DO i = 1 , m
            READ (Nin,FMT=*) (ctmp(i,j),j=1,n)
         ENDDO
         DO imla = 1 , 3
            DO imlad = 1 , 3
               DO imlb = 1 , 3
                  DO imlc = 1 , 3
                     DO itrana = 1 , 2
                        DO itranb = 1 , 2
                           DO isgn = -1 , 1 , 2
                              IF ( itrana==1 ) trana = 'N'
                              IF ( itrana==2 ) trana = 'C'
                              IF ( itranb==1 ) tranb = 'N'
                              IF ( itranb==2 ) tranb = 'C'
                              tnrm = ZERO
                              DO i = 1 , m
                                 DO j = 1 , m
                                    a(i,j) = atmp(i,j)*vm1(imla)
                                    tnrm = MAX(tnrm,ABS(a(i,j)))
                                 ENDDO
                                 a(i,i) = a(i,i)*vm2(imlad)
                                 tnrm = MAX(tnrm,ABS(a(i,i)))
                              ENDDO
                              DO i = 1 , n
                                 DO j = 1 , n
                                    b(i,j) = btmp(i,j)*vm1(imlb)
                                    tnrm = MAX(tnrm,ABS(b(i,j)))
                                 ENDDO
                              ENDDO
                              IF ( tnrm==ZERO ) tnrm = ONE
                              DO i = 1 , m
                                 DO j = 1 , n
                                    c(i,j) = ctmp(i,j)*vm1(imlc)
                                    csav(i,j) = c(i,j)
                                 ENDDO
                              ENDDO
                              Knt = Knt + 1
                              CALL CTRSYL(trana,tranb,isgn,m,n,a,LDT,b, &
     &                           LDT,c,LDT,scale,info)
                              IF ( info/=0 ) Ninfo = Ninfo + 1
                              xnrm = CLANGE('M',m,n,c,LDT,dum)
                              rmul = CONE
                              IF ( xnrm>ONE .AND. tnrm>ONE ) THEN
                                 IF ( xnrm>bignum/tnrm ) THEN
                                    rmul = MAX(xnrm,tnrm)
                                    rmul = CONE/rmul
                                 ENDIF
                              ENDIF
                              CALL CGEMM(trana,'N',m,n,m,rmul,a,LDT,c,  &
     &                           LDT,-scale*rmul,csav,LDT)
                              CALL CGEMM('N',tranb,m,n,n,REAL(isgn)     &
     &                           *rmul,c,LDT,b,LDT,CONE,csav,LDT)
                              res1 = CLANGE('M',m,n,csav,LDT,dum)
                              res = res1/MAX(smlnum,smlnum*xnrm,        &
     &                              ((ABS(rmul)*tnrm)*eps)*xnrm)
                              IF ( res>Rmax ) THEN
                                 Lmax = Knt
                                 Rmax = res
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     End of CGET35
!
      END SUBROUTINE CGET35
