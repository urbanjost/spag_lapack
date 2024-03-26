!*==clacsg.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!
!
!
      SUBROUTINE CLACSG(M,P,Q,Theta,Iseed,X,Ldx,Work)
      IMPLICIT NONE
!*--CLACSG355
!
      INTEGER Ldx , M , P , Q
      INTEGER Iseed(4)
      REAL Theta(*)
      COMPLEX Work(*) , X(Ldx,*)
!
      COMPLEX ONE , ZERO
      PARAMETER (ONE=(1.0E0,0.0E0),ZERO=(0.0E0,0.0E0))
!
      INTEGER i , info , r
!
      r = MIN(P,M-P,Q,M-Q)
!
      CALL CLASET('Full',M,M,ZERO,ZERO,X,Ldx)
!
      DO i = 1 , MIN(P,Q) - r
         X(i,i) = ONE
      ENDDO
      DO i = 1 , r
         X(MIN(P,Q)-r+i,MIN(P,Q)-r+i) = CMPLX(COS(Theta(i)),0.0E0)
      ENDDO
      DO i = 1 , MIN(P,M-Q) - r
         X(P-i+1,M-i+1) = -ONE
      ENDDO
      DO i = 1 , r
         X(P-(MIN(P,M-Q)-r)+1-i,M-(MIN(P,M-Q)-r)+1-i)                   &
     &      = CMPLX(-SIN(Theta(r-i+1)),0.0E0)
      ENDDO
      DO i = 1 , MIN(M-P,Q) - r
         X(M-i+1,Q-i+1) = ONE
      ENDDO
      DO i = 1 , r
         X(M-(MIN(M-P,Q)-r)+1-i,Q-(MIN(M-P,Q)-r)+1-i)                   &
     &      = CMPLX(SIN(Theta(r-i+1)),0.0E0)
      ENDDO
      DO i = 1 , MIN(M-P,M-Q) - r
         X(P+i,Q+i) = ONE
      ENDDO
      DO i = 1 , r
         X(P+(MIN(M-P,M-Q)-r)+i,Q+(MIN(M-P,M-Q)-r)+i)                   &
     &      = CMPLX(COS(Theta(i)),0.0E0)
      ENDDO
      CALL claror('Left','No init',P,M,X,Ldx,Iseed,Work,info)
      CALL claror('Left','No init',M-P,M,X(P+1,1),Ldx,Iseed,Work,info)
      CALL claror('Right','No init',M,Q,X,Ldx,Iseed,Work,info)
      CALL claror('Right','No init',M,M-Q,X(1,Q+1),Ldx,Iseed,Work,info)
!
      END SUBROUTINE CLACSG
