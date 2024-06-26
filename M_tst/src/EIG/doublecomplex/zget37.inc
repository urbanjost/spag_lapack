!*==zget37.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b zget37
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zget37( RMAX, LMAX, NINFO, KNT, NIN )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, NIN
!       ..
!       .. Array Arguments ..
!       INTEGER            LMAX( 3 ), NINFO( 3 )
!       DOUBLE PRECISION   RMAX( 3 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zget37 tests ZTRSNA, a routine for estimating condition numbers of
!> eigenvalues and/or right eigenvectors of a matrix.
!>
!> The test matrices are read from a file with logical unit number NIN.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] RMAX
!> \verbatim
!>          RMAX is DOUBLE PRECISION array, dimension (3)
!>          Value of the largest test ratio.
!>          RMAX(1) = largest ratio comparing different calls to ZTRSNA
!>          RMAX(2) = largest error in reciprocal condition
!>                    numbers taking their conditioning into account
!>          RMAX(3) = largest error in reciprocal condition
!>                    numbers not taking their conditioning into
!>                    account (may be larger than RMAX(2))
!> \endverbatim
!>
!> \param[out] LMAX
!> \verbatim
!>          LMAX is INTEGER array, dimension (3)
!>          LMAX(i) is example number where largest test ratio
!>          RMAX(i) is achieved. Also:
!>          If ZGEHRD returns INFO nonzero on example i, LMAX(1)=i
!>          If ZHSEQR returns INFO nonzero on example i, LMAX(2)=i
!>          If ZTRSNA returns INFO nonzero on example i, LMAX(3)=i
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER array, dimension (3)
!>          NINFO(1) = No. of times ZGEHRD returned INFO nonzero
!>          NINFO(2) = No. of times ZHSEQR returned INFO nonzero
!>          NINFO(3) = No. of times ZTRSNA returned INFO nonzero
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
!>          Input logical unit number
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
      SUBROUTINE zget37(Rmax,Lmax,Ninfo,Knt,Nin)
      IMPLICIT NONE
!*--ZGET3794
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Knt , Nin
!     ..
!     .. Array Arguments ..
      INTEGER Lmax(3) , Ninfo(3)
      DOUBLE PRECISION Rmax(3)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      DOUBLE PRECISION EPSIN
      PARAMETER (EPSIN=5.9605D-8)
      INTEGER LDT , LWORK
      PARAMETER (LDT=20,LWORK=2*LDT*(10+LDT))
!     ..
!     .. Local Scalars ..
      INTEGER i , icmp , info , iscl , isrt , j , kmin , m , n
      DOUBLE PRECISION bignum , eps , smlnum , tnrm , tol , tolin , v , &
     &                 vcmin , vmax , vmin , vmul
!     ..
!     .. Local Arrays ..
      LOGICAL select(LDT)
      INTEGER lcmp(3)
      DOUBLE PRECISION dum(1) , rwork(2*LDT) , s(LDT) , sep(LDT) ,      &
     &                 sepin(LDT) , septmp(LDT) , sin(LDT) , stmp(LDT) ,&
     &                 val(3) , wiin(LDT) , wrin(LDT) , wsrt(LDT)
      COMPLEX*16 cdum(1) , le(LDT,LDT) , re(LDT,LDT) , t(LDT,LDT) ,     &
     &           tmp(LDT,LDT) , w(LDT) , work(LWORK) , wtmp(LDT)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL DLAMCH , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DLABAD , DSCAL , ZCOPY , ZDSCAL , ZGEHRD ,       ZHSEQR , ZLACPY , ZTREVC , ZTRSNA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DIMAG , MAX , SQRT
!     ..
!     .. Executable Statements ..
!
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')/eps
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
!
!     EPSIN = 2**(-24) = precision to which input data computed
!
      eps = MAX(eps,EPSIN)
      Rmax(1) = ZERO
      Rmax(2) = ZERO
      Rmax(3) = ZERO
      Lmax(1) = 0
      Lmax(2) = 0
      Lmax(3) = 0
      Knt = 0
      Ninfo(1) = 0
      Ninfo(2) = 0
      Ninfo(3) = 0
      val(1) = SQRT(smlnum)
      val(2) = ONE
      val(3) = SQRT(bignum)
      DO
!
!     Read input data until N=0.  Assume input eigenvalues are sorted
!     lexicographically (increasing by real part if ISRT = 0,
!     increasing by imaginary part if ISRT = 1)
!
         READ (Nin,FMT=*) n , isrt
         IF ( n==0 ) RETURN
         DO i = 1 , n
            READ (Nin,FMT=*) (tmp(i,j),j=1,n)
         ENDDO
         DO i = 1 , n
            READ (Nin,FMT=*) wrin(i) , wiin(i) , sin(i) , sepin(i)
         ENDDO
         tnrm = ZLANGE('M',n,n,tmp,LDT,rwork)
         DO iscl = 1 , 3
!
!        Scale input matrix
!
            Knt = Knt + 1
            CALL ZLACPY('F',n,n,tmp,LDT,t,LDT)
            vmul = val(iscl)
            DO i = 1 , n
               CALL ZDSCAL(n,vmul,t(1,i),1)
            ENDDO
            IF ( tnrm==ZERO ) vmul = ONE
!
!        Compute eigenvalues and eigenvectors
!
            CALL ZGEHRD(n,1,n,t,LDT,work(1),work(n+1),LWORK-n,info)
            IF ( info/=0 ) THEN
               Lmax(1) = Knt
               Ninfo(1) = Ninfo(1) + 1
               CYCLE
            ENDIF
            DO j = 1 , n - 2
               DO i = j + 2 , n
                  t(i,j) = ZERO
               ENDDO
            ENDDO
!
!        Compute Schur form
!
            CALL ZHSEQR('S','N',n,1,n,t,LDT,w,cdum,1,work,LWORK,info)
            IF ( info/=0 ) THEN
               Lmax(2) = Knt
               Ninfo(2) = Ninfo(2) + 1
               CYCLE
            ENDIF
!
!        Compute eigenvectors
!
            DO i = 1 , n
               select(i) = .TRUE.
            ENDDO
            CALL ZTREVC('B','A',select,n,t,LDT,le,LDT,re,LDT,n,m,work,  &
     &                  rwork,info)
!
!        Compute condition numbers
!
            CALL ZTRSNA('B','A',select,n,t,LDT,le,LDT,re,LDT,s,sep,n,m, &
     &                  work,n,rwork,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
!
!        Sort eigenvalues and condition numbers lexicographically
!        to compare with inputs
!
            CALL ZCOPY(n,w,1,wtmp,1)
            IF ( isrt==0 ) THEN
!
!           Sort by increasing real part
!
               DO i = 1 , n
                  wsrt(i) = DBLE(w(i))
               ENDDO
            ELSE
!
!           Sort by increasing imaginary part
!
               DO i = 1 , n
                  wsrt(i) = DIMAG(w(i))
               ENDDO
            ENDIF
            CALL DCOPY(n,s,1,stmp,1)
            CALL DCOPY(n,sep,1,septmp,1)
            CALL DSCAL(n,ONE/vmul,septmp,1)
            DO i = 1 , n - 1
               kmin = i
               vmin = wsrt(i)
               DO j = i + 1 , n
                  IF ( wsrt(j)<vmin ) THEN
                     kmin = j
                     vmin = wsrt(j)
                  ENDIF
               ENDDO
               wsrt(kmin) = wsrt(i)
               wsrt(i) = vmin
               vcmin = wtmp(i)
               wtmp(i) = w(kmin)
               wtmp(kmin) = vcmin
               vmin = stmp(kmin)
               stmp(kmin) = stmp(i)
               stmp(i) = vmin
               vmin = septmp(kmin)
               septmp(kmin) = septmp(i)
               septmp(i) = vmin
            ENDDO
!
!        Compare condition numbers for eigenvalues
!        taking their condition numbers into account
!
            v = MAX(TWO*DBLE(n)*eps*tnrm,smlnum)
            IF ( tnrm==ZERO ) v = ONE
            DO i = 1 , n
               IF ( v>septmp(i) ) THEN
                  tol = ONE
               ELSE
                  tol = v/septmp(i)
               ENDIF
               IF ( v>sepin(i) ) THEN
                  tolin = ONE
               ELSE
                  tolin = v/sepin(i)
               ENDIF
               tol = MAX(tol,smlnum/eps)
               tolin = MAX(tolin,smlnum/eps)
               IF ( eps*(sin(i)-tolin)>stmp(i)+tol ) THEN
                  vmax = ONE/eps
               ELSEIF ( sin(i)-tolin>stmp(i)+tol ) THEN
                  vmax = (sin(i)-tolin)/(stmp(i)+tol)
               ELSEIF ( sin(i)+tolin<eps*(stmp(i)-tol) ) THEN
                  vmax = ONE/eps
               ELSEIF ( sin(i)+tolin<stmp(i)-tol ) THEN
                  vmax = (stmp(i)-tol)/(sin(i)+tolin)
               ELSE
                  vmax = ONE
               ENDIF
               IF ( vmax>Rmax(2) ) THEN
                  Rmax(2) = vmax
                  IF ( Ninfo(2)==0 ) Lmax(2) = Knt
               ENDIF
            ENDDO
!
!        Compare condition numbers for eigenvectors
!        taking their condition numbers into account
!
            DO i = 1 , n
               IF ( v>septmp(i)*stmp(i) ) THEN
                  tol = septmp(i)
               ELSE
                  tol = v/stmp(i)
               ENDIF
               IF ( v>sepin(i)*sin(i) ) THEN
                  tolin = sepin(i)
               ELSE
                  tolin = v/sin(i)
               ENDIF
               tol = MAX(tol,smlnum/eps)
               tolin = MAX(tolin,smlnum/eps)
               IF ( eps*(sepin(i)-tolin)>septmp(i)+tol ) THEN
                  vmax = ONE/eps
               ELSEIF ( sepin(i)-tolin>septmp(i)+tol ) THEN
                  vmax = (sepin(i)-tolin)/(septmp(i)+tol)
               ELSEIF ( sepin(i)+tolin<eps*(septmp(i)-tol) ) THEN
                  vmax = ONE/eps
               ELSEIF ( sepin(i)+tolin<septmp(i)-tol ) THEN
                  vmax = (septmp(i)-tol)/(sepin(i)+tolin)
               ELSE
                  vmax = ONE
               ENDIF
               IF ( vmax>Rmax(2) ) THEN
                  Rmax(2) = vmax
                  IF ( Ninfo(2)==0 ) Lmax(2) = Knt
               ENDIF
            ENDDO
!
!        Compare condition numbers for eigenvalues
!        without taking their condition numbers into account
!
            DO i = 1 , n
               IF ( sin(i)<=DBLE(2*n)*eps .AND. stmp(i)<=DBLE(2*n)*eps )&
     &              THEN
                  vmax = ONE
               ELSEIF ( eps*sin(i)>stmp(i) ) THEN
                  vmax = ONE/eps
               ELSEIF ( sin(i)>stmp(i) ) THEN
                  vmax = sin(i)/stmp(i)
               ELSEIF ( sin(i)<eps*stmp(i) ) THEN
                  vmax = ONE/eps
               ELSEIF ( sin(i)<stmp(i) ) THEN
                  vmax = stmp(i)/sin(i)
               ELSE
                  vmax = ONE
               ENDIF
               IF ( vmax>Rmax(3) ) THEN
                  Rmax(3) = vmax
                  IF ( Ninfo(3)==0 ) Lmax(3) = Knt
               ENDIF
            ENDDO
!
!        Compare condition numbers for eigenvectors
!        without taking their condition numbers into account
!
            DO i = 1 , n
               IF ( sepin(i)<=v .AND. septmp(i)<=v ) THEN
                  vmax = ONE
               ELSEIF ( eps*sepin(i)>septmp(i) ) THEN
                  vmax = ONE/eps
               ELSEIF ( sepin(i)>septmp(i) ) THEN
                  vmax = sepin(i)/septmp(i)
               ELSEIF ( sepin(i)<eps*septmp(i) ) THEN
                  vmax = ONE/eps
               ELSEIF ( sepin(i)<septmp(i) ) THEN
                  vmax = septmp(i)/sepin(i)
               ELSE
                  vmax = ONE
               ENDIF
               IF ( vmax>Rmax(3) ) THEN
                  Rmax(3) = vmax
                  IF ( Ninfo(3)==0 ) Lmax(3) = Knt
               ENDIF
            ENDDO
!
!        Compute eigenvalue condition numbers only and compare
!
            vmax = ZERO
            dum(1) = -ONE
            CALL DCOPY(n,dum,0,stmp,1)
            CALL DCOPY(n,dum,0,septmp,1)
            CALL ZTRSNA('E','A',select,n,t,LDT,le,LDT,re,LDT,stmp,      &
     &                  septmp,n,m,work,n,rwork,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            DO i = 1 , n
               IF ( stmp(i)/=s(i) ) vmax = ONE/eps
               IF ( septmp(i)/=dum(1) ) vmax = ONE/eps
            ENDDO
!
!        Compute eigenvector condition numbers only and compare
!
            CALL DCOPY(n,dum,0,stmp,1)
            CALL DCOPY(n,dum,0,septmp,1)
            CALL ZTRSNA('V','A',select,n,t,LDT,le,LDT,re,LDT,stmp,      &
     &                  septmp,n,m,work,n,rwork,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            DO i = 1 , n
               IF ( stmp(i)/=dum(1) ) vmax = ONE/eps
               IF ( septmp(i)/=sep(i) ) vmax = ONE/eps
            ENDDO
!
!        Compute all condition numbers using SELECT and compare
!
            DO i = 1 , n
               select(i) = .TRUE.
            ENDDO
            CALL DCOPY(n,dum,0,stmp,1)
            CALL DCOPY(n,dum,0,septmp,1)
            CALL ZTRSNA('B','S',select,n,t,LDT,le,LDT,re,LDT,stmp,      &
     &                  septmp,n,m,work,n,rwork,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            DO i = 1 , n
               IF ( septmp(i)/=sep(i) ) vmax = ONE/eps
               IF ( stmp(i)/=s(i) ) vmax = ONE/eps
            ENDDO
!
!        Compute eigenvalue condition numbers using SELECT and compare
!
            CALL DCOPY(n,dum,0,stmp,1)
            CALL DCOPY(n,dum,0,septmp,1)
            CALL ZTRSNA('E','S',select,n,t,LDT,le,LDT,re,LDT,stmp,      &
     &                  septmp,n,m,work,n,rwork,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            DO i = 1 , n
               IF ( stmp(i)/=s(i) ) vmax = ONE/eps
               IF ( septmp(i)/=dum(1) ) vmax = ONE/eps
            ENDDO
!
!        Compute eigenvector condition numbers using SELECT and compare
!
            CALL DCOPY(n,dum,0,stmp,1)
            CALL DCOPY(n,dum,0,septmp,1)
            CALL ZTRSNA('V','S',select,n,t,LDT,le,LDT,re,LDT,stmp,      &
     &                  septmp,n,m,work,n,rwork,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            DO i = 1 , n
               IF ( stmp(i)/=dum(1) ) vmax = ONE/eps
               IF ( septmp(i)/=sep(i) ) vmax = ONE/eps
            ENDDO
            IF ( vmax>Rmax(1) ) THEN
               Rmax(1) = vmax
               IF ( Ninfo(1)==0 ) Lmax(1) = Knt
            ENDIF
!
!        Select second and next to last eigenvalues
!
            DO i = 1 , n
               select(i) = .FALSE.
            ENDDO
            icmp = 0
            IF ( n>1 ) THEN
               icmp = 1
               lcmp(1) = 2
               select(2) = .TRUE.
               CALL ZCOPY(n,re(1,2),1,re(1,1),1)
               CALL ZCOPY(n,le(1,2),1,le(1,1),1)
            ENDIF
            IF ( n>3 ) THEN
               icmp = 2
               lcmp(2) = n - 1
               select(n-1) = .TRUE.
               CALL ZCOPY(n,re(1,n-1),1,re(1,2),1)
               CALL ZCOPY(n,le(1,n-1),1,le(1,2),1)
            ENDIF
!
!        Compute all selected condition numbers
!
            CALL DCOPY(icmp,dum,0,stmp,1)
            CALL DCOPY(icmp,dum,0,septmp,1)
            CALL ZTRSNA('B','S',select,n,t,LDT,le,LDT,re,LDT,stmp,      &
     &                  septmp,n,m,work,n,rwork,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            DO i = 1 , icmp
               j = lcmp(i)
               IF ( septmp(i)/=sep(j) ) vmax = ONE/eps
               IF ( stmp(i)/=s(j) ) vmax = ONE/eps
            ENDDO
!
!        Compute selected eigenvalue condition numbers
!
            CALL DCOPY(icmp,dum,0,stmp,1)
            CALL DCOPY(icmp,dum,0,septmp,1)
            CALL ZTRSNA('E','S',select,n,t,LDT,le,LDT,re,LDT,stmp,      &
     &                  septmp,n,m,work,n,rwork,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            DO i = 1 , icmp
               j = lcmp(i)
               IF ( stmp(i)/=s(j) ) vmax = ONE/eps
               IF ( septmp(i)/=dum(1) ) vmax = ONE/eps
            ENDDO
!
!        Compute selected eigenvector condition numbers
!
            CALL DCOPY(icmp,dum,0,stmp,1)
            CALL DCOPY(icmp,dum,0,septmp,1)
            CALL ZTRSNA('V','S',select,n,t,LDT,le,LDT,re,LDT,stmp,      &
     &                  septmp,n,m,work,n,rwork,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            DO i = 1 , icmp
               j = lcmp(i)
               IF ( stmp(i)/=dum(1) ) vmax = ONE/eps
               IF ( septmp(i)/=sep(j) ) vmax = ONE/eps
            ENDDO
            IF ( vmax>Rmax(1) ) THEN
               Rmax(1) = vmax
               IF ( Ninfo(1)==0 ) Lmax(1) = Knt
            ENDIF
         ENDDO
      ENDDO
!
!     End of zget37
!
      END SUBROUTINE zget37
