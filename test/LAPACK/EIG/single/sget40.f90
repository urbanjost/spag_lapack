!*==sget40.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b SGET40
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGET40( RMAX, LMAX, NINFO, KNT, NIN )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NIN
!       REAL   RMAX
!       ..
!       .. Array Arguments ..
!       INTEGER            NINFO( 3 )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGET40 tests STGEXC, a routine for swapping adjacent blocks (either
!> 1 by 1 or 2 by 2) on the diagonal of a pencil in real generalized Schur form.
!> Thus, STGEXC computes an orthogonal matrices Q and Z such that
!>
!>     Q' * ( [ A B ], [ D E ] ) * Z  = ( [ C1 B1 ], [ F1 E1 ] )
!>          ( [ 0 C ]  [   F ] )        ( [ 0  A1 ]  [    D1]  )
!>
!> where (C1,F1) is similar to (C,F) and (A1,D1) is similar to (A,D).
!> Both (A,D) and (C,F) are assumed to be in standard form
!> and (A1,D1) and (C1,F1) are returned with the
!> same properties.
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
!> \param[out] NIN
!> \verbatim
!>          NINFO is INTEGER
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
      SUBROUTINE SGET40(Rmax,Lmax,Ninfo,Knt,Nin)
      IMPLICIT NONE
!*--SGET4087
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Knt , Lmax , Nin
      REAL Rmax
!     ..
!     .. Array Arguments ..
      INTEGER Ninfo(3)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0,ONE=1.0)
      INTEGER LDT , LWORK
      PARAMETER (LDT=10,LWORK=100+4*LDT+16)
!     ..
!     .. Local Scalars ..
      INTEGER i , ifst , ifst1 , ifst2 , ifstsv , ilst , ilst1 , ilst2 ,&
     &        ilstsv , info1 , info2 , j , n
      REAL eps , res
!     ..
!     .. Local Arrays ..
      REAL q(LDT,LDT) , z(LDT,LDT) , result(4) , t(LDT,LDT) ,           &
     &     t1(LDT,LDT) , t2(LDT,LDT) , s(LDT,LDT) , s1(LDT,LDT) ,       &
     &     s2(LDT,LDT) , tmp(LDT,LDT) , work(LWORK)
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SGET51 , SLACPY , SLASET , STGEXC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SIGN
!     ..
!     .. Executable Statements ..
!
      eps = SLAMCH('P')
      Rmax = ZERO
      Lmax = 0
      Knt = 0
      Ninfo(1) = 0
      Ninfo(2) = 0
      Ninfo(3) = 0
      DO
!
!     Read input data until N=0
!
         READ (Nin,FMT=*) n , ifst , ilst
         IF ( n==0 ) RETURN
         Knt = Knt + 1
         DO i = 1 , n
            READ (Nin,FMT=*) (tmp(i,j),j=1,n)
         ENDDO
         CALL SLACPY('F',n,n,tmp,LDT,t,LDT)
         CALL SLACPY('F',n,n,tmp,LDT,t1,LDT)
         CALL SLACPY('F',n,n,tmp,LDT,t2,LDT)
         DO i = 1 , n
            READ (Nin,FMT=*) (tmp(i,j),j=1,n)
         ENDDO
         CALL SLACPY('F',n,n,tmp,LDT,s,LDT)
         CALL SLACPY('F',n,n,tmp,LDT,s1,LDT)
         CALL SLACPY('F',n,n,tmp,LDT,s2,LDT)
         ifstsv = ifst
         ilstsv = ilst
         ifst1 = ifst
         ilst1 = ilst
         ifst2 = ifst
         ilst2 = ilst
         res = ZERO
!
!     Test without accumulating Q and Z
!
         CALL SLASET('Full',n,n,ZERO,ONE,q,LDT)
         CALL SLASET('Full',n,n,ZERO,ONE,z,LDT)
         CALL STGEXC(.FALSE.,.FALSE.,n,t1,LDT,s1,LDT,q,LDT,z,LDT,ifst1, &
     &               ilst1,work,LWORK,info1)
         DO i = 1 , n
            DO j = 1 , n
               IF ( i==j .AND. q(i,j)/=ONE ) res = res + ONE/eps
               IF ( i/=j .AND. q(i,j)/=ZERO ) res = res + ONE/eps
               IF ( i==j .AND. z(i,j)/=ONE ) res = res + ONE/eps
               IF ( i/=j .AND. z(i,j)/=ZERO ) res = res + ONE/eps
            ENDDO
         ENDDO
!
!     Test with accumulating Q
!
         CALL SLASET('Full',n,n,ZERO,ONE,q,LDT)
         CALL SLASET('Full',n,n,ZERO,ONE,z,LDT)
         CALL STGEXC(.TRUE.,.TRUE.,n,t2,LDT,s2,LDT,q,LDT,z,LDT,ifst2,   &
     &               ilst2,work,LWORK,info2)
!
!     Compare T1 with T2 and S1 with S2
!
         DO i = 1 , n
            DO j = 1 , n
               IF ( t1(i,j)/=t2(i,j) ) res = res + ONE/eps
               IF ( s1(i,j)/=s2(i,j) ) res = res + ONE/eps
            ENDDO
         ENDDO
         IF ( ifst1/=ifst2 ) res = res + ONE/eps
         IF ( ilst1/=ilst2 ) res = res + ONE/eps
         IF ( info1/=info2 ) res = res + ONE/eps
!
!     Test orthogonality of Q and Z and backward error on T2 and S2
!
         CALL SGET51(1,n,t,LDT,t2,LDT,q,LDT,z,LDT,work,result(1))
         CALL SGET51(1,n,s,LDT,s2,LDT,q,LDT,z,LDT,work,result(2))
         CALL SGET51(3,n,t,LDT,t2,LDT,q,LDT,q,LDT,work,result(3))
         CALL SGET51(3,n,t,LDT,t2,LDT,z,LDT,z,LDT,work,result(4))
!
!     Read next matrix pair
!
         res = res + result(1) + result(2) + result(3) + result(4)
      ENDDO
!
!     End of SGET40
!
      END SUBROUTINE SGET40
