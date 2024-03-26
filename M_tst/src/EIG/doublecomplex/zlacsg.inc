!*==zlacsg.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!
!
!
      SUBROUTINE ZLACSG(M,P,Q,Theta,Iseed,X,Ldx,Work)
      IMPLICIT NONE
!*--ZLACSG356
!
      INTEGER Ldx , M , P , Q
      INTEGER Iseed(4)
      DOUBLE PRECISION Theta(*)
      COMPLEX*16 Work(*) , X(Ldx,*)
!
      COMPLEX*16 ONE , ZERO
      PARAMETER (ONE=(1.0D0,0.0D0),ZERO=(0.0D0,0.0D0))
!
      INTEGER i , info , r
!
      r = MIN(P,M-P,Q,M-Q)
!
      CALL ZLASET('Full',M,M,ZERO,ZERO,X,Ldx)
!
      DO i = 1 , MIN(P,Q) - r
         X(i,i) = ONE
      ENDDO
      DO i = 1 , r
         X(MIN(P,Q)-r+i,MIN(P,Q)-r+i) = DCMPLX(COS(Theta(i)),0.0D0)
      ENDDO
      DO i = 1 , MIN(P,M-Q) - r
         X(P-i+1,M-i+1) = -ONE
      ENDDO
      DO i = 1 , r
         X(P-(MIN(P,M-Q)-r)+1-i,M-(MIN(P,M-Q)-r)+1-i)                   &
     &      = DCMPLX(-SIN(Theta(r-i+1)),0.0D0)
      ENDDO
      DO i = 1 , MIN(M-P,Q) - r
         X(M-i+1,Q-i+1) = ONE
      ENDDO
      DO i = 1 , r
         X(M-(MIN(M-P,Q)-r)+1-i,Q-(MIN(M-P,Q)-r)+1-i)                   &
     &      = DCMPLX(SIN(Theta(r-i+1)),0.0D0)
      ENDDO
      DO i = 1 , MIN(M-P,M-Q) - r
         X(P+i,Q+i) = ONE
      ENDDO
      DO i = 1 , r
         X(P+(MIN(M-P,M-Q)-r)+i,Q+(MIN(M-P,M-Q)-r)+i)                   &
     &      = DCMPLX(COS(Theta(i)),0.0D0)
      ENDDO
      CALL zlaror('Left','No init',P,M,X,Ldx,Iseed,Work,info)
      CALL zlaror('Left','No init',M-P,M,X(P+1,1),Ldx,Iseed,Work,info)
      CALL zlaror('Right','No init',M,Q,X,Ldx,Iseed,Work,info)
      CALL zlaror('Right','No init',M,M-Q,X(1,Q+1),Ldx,Iseed,Work,info)
!
      END SUBROUTINE ZLACSG
