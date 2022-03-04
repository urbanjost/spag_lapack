!*==zget36.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZGET36
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGET36( RMAX, LMAX, NINFO, KNT, NIN )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, LMAX, NIN, NINFO
!       DOUBLE PRECISION   RMAX
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGET36 tests ZTREXC, a routine for reordering diagonal entries of a
!> matrix in complex Schur form. Thus, ZLAEXC computes a unitary matrix
!> Q such that
!>
!>    Q' * T1 * Q  = T2
!>
!> and where one of the diagonal blocks of T1 (the one at row IFST) has
!> been moved to position ILST.
!>
!> The test code verifies that the residual Q'*T1*Q-T2 is small, that T2
!> is in Schur form, and that the final position of the IFST block is
!> ILST.
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZGET36(Rmax,Lmax,Ninfo,Knt,Nin)
      IMPLICIT NONE
!*--ZGET3689
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Knt , Lmax , Nin , Ninfo
      DOUBLE PRECISION Rmax
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
      INTEGER LDT , LWORK
      PARAMETER (LDT=10,LWORK=2*LDT*LDT)
!     ..
!     .. Local Scalars ..
      INTEGER i , ifst , ilst , info1 , info2 , j , n
      DOUBLE PRECISION eps , res
      COMPLEX*16 ctemp
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION result(2) , rwork(LDT)
      COMPLEX*16 diag(LDT) , q(LDT,LDT) , t1(LDT,LDT) , t2(LDT,LDT) ,   &
     &           tmp(LDT,LDT) , work(LWORK)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL ZCOPY , ZHST01 , ZLACPY , ZLASET , ZTREXC
!     ..
!     .. Executable Statements ..
!
      eps = DLAMCH('P')
      Rmax = ZERO
      Lmax = 0
      Knt = 0
      Ninfo = 0
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
         CALL ZLACPY('F',n,n,tmp,LDT,t1,LDT)
         CALL ZLACPY('F',n,n,tmp,LDT,t2,LDT)
         res = ZERO
!
!     Test without accumulating Q
!
         CALL ZLASET('Full',n,n,CZERO,CONE,q,LDT)
         CALL ZTREXC('N',n,t1,LDT,q,LDT,ifst,ilst,info1)
         DO i = 1 , n
            DO j = 1 , n
               IF ( i==j .AND. q(i,j)/=CONE ) res = res + ONE/eps
               IF ( i/=j .AND. q(i,j)/=CZERO ) res = res + ONE/eps
            ENDDO
         ENDDO
!
!     Test with accumulating Q
!
         CALL ZLASET('Full',n,n,CZERO,CONE,q,LDT)
         CALL ZTREXC('V',n,t2,LDT,q,LDT,ifst,ilst,info2)
!
!     Compare T1 with T2
!
         DO i = 1 , n
            DO j = 1 , n
               IF ( t1(i,j)/=t2(i,j) ) res = res + ONE/eps
            ENDDO
         ENDDO
         IF ( info1/=0 .OR. info2/=0 ) Ninfo = Ninfo + 1
         IF ( info1/=info2 ) res = res + ONE/eps
!
!     Test for successful reordering of T2
!
         CALL ZCOPY(n,tmp,LDT+1,diag,1)
         IF ( ifst<ilst ) THEN
            DO i = ifst + 1 , ilst
               ctemp = diag(i)
               diag(i) = diag(i-1)
               diag(i-1) = ctemp
            ENDDO
         ELSEIF ( ifst>ilst ) THEN
            DO i = ifst - 1 , ilst , -1
               ctemp = diag(i+1)
               diag(i+1) = diag(i)
               diag(i) = ctemp
            ENDDO
         ENDIF
         DO i = 1 , n
            IF ( t2(i,i)/=diag(i) ) res = res + ONE/eps
         ENDDO
!
!     Test for small residual, and orthogonality of Q
!
         CALL ZHST01(n,1,n,tmp,LDT,t2,LDT,q,LDT,work,LWORK,rwork,result)
         res = res + result(1) + result(2)
!
!     Test for T2 being in Schur form
!
         DO j = 1 , n - 1
            DO i = j + 1 , n
               IF ( t2(i,j)/=CZERO ) res = res + ONE/eps
            ENDDO
         ENDDO
         IF ( res>Rmax ) THEN
            Rmax = res
            Lmax = Knt
         ENDIF
      ENDDO
!
!     End of ZGET36
!
      END SUBROUTINE ZGET36
