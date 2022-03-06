!*==dget36.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DGET36
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET36( RMAX, LMAX, NINFO, KNT, NIN )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NIN
!       DOUBLE PRECISION   RMAX
!       ..
!       .. Array Arguments ..
!       INTEGER            NINFO( 3 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGET36 tests DTREXC, a routine for moving blocks (either 1 by 1 or
!> 2 by 2) on the diagonal of a matrix in real Schur form.  Thus, DLAEXC
!> computes an orthogonal matrix Q such that
!>
!>    Q' * T1 * Q  = T2
!>
!> and where one of the diagonal blocks of T1 (the one at row IFST) has
!> been moved to position ILST.
!>
!> The test code verifies that the residual Q'*T1*Q-T2 is small, that T2
!> is in Schur form, and that the final position of the IFST block is
!> ILST (within +-1).
!>
!> The test matrices are read from a file with logical unit number NIN.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] RMAX
!> \verbatim
!>          RMAX is DOUBLE PRECISION
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
!>          NINFO is INTEGER array, dimension (3)
!>          NINFO(J) is the number of examples where INFO=J.
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DGET36(Rmax,Lmax,Ninfo,Knt,Nin)
      IMPLICIT NONE
!*--DGET3692
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Knt , Lmax , Nin
      DOUBLE PRECISION Rmax
!     ..
!     .. Array Arguments ..
      INTEGER Ninfo(3)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      INTEGER LDT , LWORK
      PARAMETER (LDT=10,LWORK=2*LDT*LDT)
!     ..
!     .. Local Scalars ..
      INTEGER i , ifst , ifst1 , ifst2 , ifstsv , ilst , ilst1 , ilst2 ,&
     &        ilstsv , info1 , info2 , j , loc , n
      DOUBLE PRECISION eps , res
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION q(LDT,LDT) , result(2) , t1(LDT,LDT) ,           &
     &                 t2(LDT,LDT) , tmp(LDT,LDT) , work(LWORK)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DHST01 , DLACPY , DLASET , DTREXC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SIGN
!     ..
!     .. Executable Statements ..
!
      eps = DLAMCH('P')
      Rmax = ZERO
      Lmax = 0
      Knt = 0
      Ninfo(1) = 0
      Ninfo(2) = 0
      Ninfo(3) = 0
!
!     Read input data until N=0
!
 100  READ (Nin,FMT=*) n , ifst , ilst
      IF ( n==0 ) RETURN
      Knt = Knt + 1
      DO i = 1 , n
         READ (Nin,FMT=*) (tmp(i,j),j=1,n)
      ENDDO
      CALL DLACPY('F',n,n,tmp,LDT,t1,LDT)
      CALL DLACPY('F',n,n,tmp,LDT,t2,LDT)
      ifstsv = ifst
      ilstsv = ilst
      ifst1 = ifst
      ilst1 = ilst
      ifst2 = ifst
      ilst2 = ilst
      res = ZERO
!
!     Test without accumulating Q
!
      CALL DLASET('Full',n,n,ZERO,ONE,q,LDT)
      CALL DTREXC('N',n,t1,LDT,q,LDT,ifst1,ilst1,work,info1)
      DO i = 1 , n
         DO j = 1 , n
            IF ( i==j .AND. q(i,j)/=ONE ) res = res + ONE/eps
            IF ( i/=j .AND. q(i,j)/=ZERO ) res = res + ONE/eps
         ENDDO
      ENDDO
!
!     Test with accumulating Q
!
      CALL DLASET('Full',n,n,ZERO,ONE,q,LDT)
      CALL DTREXC('V',n,t2,LDT,q,LDT,ifst2,ilst2,work,info2)
!
!     Compare T1 with T2
!
      DO i = 1 , n
         DO j = 1 , n
            IF ( t1(i,j)/=t2(i,j) ) res = res + ONE/eps
         ENDDO
      ENDDO
      IF ( ifst1/=ifst2 ) res = res + ONE/eps
      IF ( ilst1/=ilst2 ) res = res + ONE/eps
      IF ( info1/=info2 ) res = res + ONE/eps
!
!     Test for successful reordering of T2
!
      IF ( info2/=0 ) THEN
         Ninfo(info2) = Ninfo(info2) + 1
      ELSE
         IF ( ABS(ifst2-ifstsv)>1 ) res = res + ONE/eps
         IF ( ABS(ilst2-ilstsv)>1 ) res = res + ONE/eps
      ENDIF
!
!     Test for small residual, and orthogonality of Q
!
      CALL DHST01(n,1,n,tmp,LDT,t2,LDT,q,LDT,work,LWORK,result)
      res = res + result(1) + result(2)
!
!     Test for T2 being in Schur form
!
      loc = 1
      DO
         IF ( t2(loc+1,loc)/=ZERO ) THEN
!
!        2 by 2 block
!
            IF ( t2(loc,loc+1)==ZERO .OR. t2(loc,loc)/=t2(loc+1,loc+1)  &
     &           .OR. SIGN(ONE,t2(loc,loc+1))==SIGN(ONE,t2(loc+1,loc)) )&
     &           res = res + ONE/eps
            DO i = loc + 2 , n
               IF ( t2(i,loc)/=ZERO ) res = res + ONE/res
               IF ( t2(i,loc+1)/=ZERO ) res = res + ONE/res
            ENDDO
            loc = loc + 2
         ELSE
!
!        1 by 1 block
!
            DO i = loc + 1 , n
               IF ( t2(i,loc)/=ZERO ) res = res + ONE/res
            ENDDO
            loc = loc + 1
         ENDIF
         IF ( loc>=n ) THEN
            IF ( res>Rmax ) THEN
               Rmax = res
               Lmax = Knt
            ENDIF
            GOTO 100
         ENDIF
      ENDDO
!
!     End of DGET36
!
      END SUBROUTINE DGET36
