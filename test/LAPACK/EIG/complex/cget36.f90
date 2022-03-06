!*==cget36.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b CGET36
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGET36( RMAX, LMAX, NINFO, KNT, NIN )
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
!> CGET36 tests CTREXC, a routine for reordering diagonal entries of a
!> matrix in complex Schur form. Thus, CLAEXC computes a unitary matrix
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
      SUBROUTINE CGET36(Rmax,Lmax,Ninfo,Knt,Nin)
      IMPLICIT NONE
!*--CGET3689
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
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
      INTEGER LDT , LWORK
      PARAMETER (LDT=10,LWORK=2*LDT*LDT)
!     ..
!     .. Local Scalars ..
      INTEGER i , ifst , ilst , info1 , info2 , j , n
      REAL eps , res
      COMPLEX ctemp
!     ..
!     .. Local Arrays ..
      REAL result(2) , rwork(LDT)
      COMPLEX diag(LDT) , q(LDT,LDT) , t1(LDT,LDT) , t2(LDT,LDT) ,      &
     &        tmp(LDT,LDT) , work(LWORK)
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CCOPY , CHST01 , CLACPY , CLASET , CTREXC
!     ..
!     .. Executable Statements ..
!
      eps = SLAMCH('P')
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
         CALL CLACPY('F',n,n,tmp,LDT,t1,LDT)
         CALL CLACPY('F',n,n,tmp,LDT,t2,LDT)
         res = ZERO
!
!     Test without accumulating Q
!
         CALL CLASET('Full',n,n,CZERO,CONE,q,LDT)
         CALL CTREXC('N',n,t1,LDT,q,LDT,ifst,ilst,info1)
         DO i = 1 , n
            DO j = 1 , n
               IF ( i==j .AND. q(i,j)/=CONE ) res = res + ONE/eps
               IF ( i/=j .AND. q(i,j)/=CZERO ) res = res + ONE/eps
            ENDDO
         ENDDO
!
!     Test with accumulating Q
!
         CALL CLASET('Full',n,n,CZERO,CONE,q,LDT)
         CALL CTREXC('V',n,t2,LDT,q,LDT,ifst,ilst,info2)
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
         CALL CCOPY(n,tmp,LDT+1,diag,1)
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
         CALL CHST01(n,1,n,tmp,LDT,t2,LDT,q,LDT,work,LWORK,rwork,result)
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
!     End of CGET36
!
      END SUBROUTINE CGET36
