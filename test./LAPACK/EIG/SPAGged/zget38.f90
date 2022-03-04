!*==zget38.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZGET38
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGET38( RMAX, LMAX, NINFO, KNT, NIN )
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
!> ZGET38 tests ZTRSEN, a routine for estimating condition numbers of a
!> cluster of eigenvalues and/or its associated right invariant subspace
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
!>          Values of the largest test ratios.
!>          RMAX(1) = largest residuals from ZHST01 or comparing
!>                    different calls to ZTRSEN
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
!>          If ZTRSEN returns INFO nonzero on example i, LMAX(3)=i
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER array, dimension (3)
!>          NINFO(1) = No. of times ZGEHRD returned INFO nonzero
!>          NINFO(2) = No. of times ZHSEQR returned INFO nonzero
!>          NINFO(3) = No. of times ZTRSEN returned INFO nonzero
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
      SUBROUTINE ZGET38(Rmax,Lmax,Ninfo,Knt,Nin)
      IMPLICIT NONE
!*--ZGET3895
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
      INTEGER LDT , LWORK
      PARAMETER (LDT=20,LWORK=2*LDT*(10+LDT))
      DOUBLE PRECISION ZERO , ONE , TWO
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0,TWO=2.0D+0)
      DOUBLE PRECISION EPSIN
      PARAMETER (EPSIN=5.9605D-8)
      COMPLEX*16 CZERO
      PARAMETER (CZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , info , iscl , isrt , itmp , j , kmin , m , n , ndim
      DOUBLE PRECISION bignum , eps , s , sep , sepin , septmp , sin ,  &
     &                 smlnum , stmp , tnrm , tol , tolin , v , vmax ,  &
     &                 vmin , vmul
!     ..
!     .. Local Arrays ..
      LOGICAL select(LDT)
      INTEGER ipnt(LDT) , iselec(LDT)
      DOUBLE PRECISION result(2) , rwork(LDT) , val(3) , wsrt(LDT)
      COMPLEX*16 q(LDT,LDT) , qsav(LDT,LDT) , qtmp(LDT,LDT) , t(LDT,LDT)&
     &           , tmp(LDT,LDT) , tsav(LDT,LDT) , tsav1(LDT,LDT) ,      &
     &           ttmp(LDT,LDT) , w(LDT) , work(LWORK) , wtmp(LDT)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL DLAMCH , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD , ZDSCAL , ZGEHRD , ZHSEQR , ZHST01 , ZLACPY ,    &
     &         ZTRSEN , ZUNGHR
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
      val(3) = SQRT(SQRT(bignum))
      DO
!
!     Read input data until N=0.  Assume input eigenvalues are sorted
!     lexicographically (increasing by real part, then decreasing by
!     imaginary part)
!
         READ (Nin,FMT=*) n , ndim , isrt
         IF ( n==0 ) RETURN
         READ (Nin,FMT=*) (iselec(i),i=1,ndim)
         DO i = 1 , n
            READ (Nin,FMT=*) (tmp(i,j),j=1,n)
         ENDDO
         READ (Nin,FMT=*) sin , sepin
!
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
            CALL ZLACPY('F',n,n,t,LDT,tsav,LDT)
!
!        Compute Schur form
!
            CALL ZGEHRD(n,1,n,t,LDT,work(1),work(n+1),LWORK-n,info)
            IF ( info/=0 ) THEN
               Lmax(1) = Knt
               Ninfo(1) = Ninfo(1) + 1
               CYCLE
            ENDIF
!
!        Generate unitary matrix
!
            CALL ZLACPY('L',n,n,t,LDT,q,LDT)
            CALL ZUNGHR(n,1,n,q,LDT,work(1),work(n+1),LWORK-n,info)
!
!        Compute Schur form
!
            DO j = 1 , n - 2
               DO i = j + 2 , n
                  t(i,j) = CZERO
               ENDDO
            ENDDO
            CALL ZHSEQR('S','V',n,1,n,t,LDT,w,q,LDT,work,LWORK,info)
            IF ( info/=0 ) THEN
               Lmax(2) = Knt
               Ninfo(2) = Ninfo(2) + 1
               CYCLE
            ENDIF
!
!        Sort, select eigenvalues
!
            DO i = 1 , n
               ipnt(i) = i
               select(i) = .FALSE.
            ENDDO
            IF ( isrt==0 ) THEN
               DO i = 1 , n
                  wsrt(i) = DBLE(w(i))
               ENDDO
            ELSE
               DO i = 1 , n
                  wsrt(i) = DIMAG(w(i))
               ENDDO
            ENDIF
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
               itmp = ipnt(i)
               ipnt(i) = ipnt(kmin)
               ipnt(kmin) = itmp
            ENDDO
            DO i = 1 , ndim
               select(ipnt(iselec(i))) = .TRUE.
            ENDDO
!
!        Compute condition numbers
!
            CALL ZLACPY('F',n,n,q,LDT,qsav,LDT)
            CALL ZLACPY('F',n,n,t,LDT,tsav1,LDT)
            CALL ZTRSEN('B','V',select,n,t,LDT,q,LDT,wtmp,m,s,sep,work, &
     &                  LWORK,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            septmp = sep/vmul
            stmp = s
!
!        Compute residuals
!
            CALL ZHST01(n,1,n,tsav,LDT,t,LDT,q,LDT,work,LWORK,rwork,    &
     &                  result)
            vmax = MAX(result(1),result(2))
            IF ( vmax>Rmax(1) ) THEN
               Rmax(1) = vmax
               IF ( Ninfo(1)==0 ) Lmax(1) = Knt
            ENDIF
!
!        Compare condition number for eigenvalue cluster
!        taking its condition number into account
!
            v = MAX(TWO*DBLE(n)*eps*tnrm,smlnum)
            IF ( tnrm==ZERO ) v = ONE
            IF ( v>septmp ) THEN
               tol = ONE
            ELSE
               tol = v/septmp
            ENDIF
            IF ( v>sepin ) THEN
               tolin = ONE
            ELSE
               tolin = v/sepin
            ENDIF
            tol = MAX(tol,smlnum/eps)
            tolin = MAX(tolin,smlnum/eps)
            IF ( eps*(sin-tolin)>stmp+tol ) THEN
               vmax = ONE/eps
            ELSEIF ( sin-tolin>stmp+tol ) THEN
               vmax = (sin-tolin)/(stmp+tol)
            ELSEIF ( sin+tolin<eps*(stmp-tol) ) THEN
               vmax = ONE/eps
            ELSEIF ( sin+tolin<stmp-tol ) THEN
               vmax = (stmp-tol)/(sin+tolin)
            ELSE
               vmax = ONE
            ENDIF
            IF ( vmax>Rmax(2) ) THEN
               Rmax(2) = vmax
               IF ( Ninfo(2)==0 ) Lmax(2) = Knt
            ENDIF
!
!        Compare condition numbers for invariant subspace
!        taking its condition number into account
!
            IF ( v>septmp*stmp ) THEN
               tol = septmp
            ELSE
               tol = v/stmp
            ENDIF
            IF ( v>sepin*sin ) THEN
               tolin = sepin
            ELSE
               tolin = v/sin
            ENDIF
            tol = MAX(tol,smlnum/eps)
            tolin = MAX(tolin,smlnum/eps)
            IF ( eps*(sepin-tolin)>septmp+tol ) THEN
               vmax = ONE/eps
            ELSEIF ( sepin-tolin>septmp+tol ) THEN
               vmax = (sepin-tolin)/(septmp+tol)
            ELSEIF ( sepin+tolin<eps*(septmp-tol) ) THEN
               vmax = ONE/eps
            ELSEIF ( sepin+tolin<septmp-tol ) THEN
               vmax = (septmp-tol)/(sepin+tolin)
            ELSE
               vmax = ONE
            ENDIF
            IF ( vmax>Rmax(2) ) THEN
               Rmax(2) = vmax
               IF ( Ninfo(2)==0 ) Lmax(2) = Knt
            ENDIF
!
!        Compare condition number for eigenvalue cluster
!        without taking its condition number into account
!
            IF ( sin<=DBLE(2*n)*eps .AND. stmp<=DBLE(2*n)*eps ) THEN
               vmax = ONE
            ELSEIF ( eps*sin>stmp ) THEN
               vmax = ONE/eps
            ELSEIF ( sin>stmp ) THEN
               vmax = sin/stmp
            ELSEIF ( sin<eps*stmp ) THEN
               vmax = ONE/eps
            ELSEIF ( sin<stmp ) THEN
               vmax = stmp/sin
            ELSE
               vmax = ONE
            ENDIF
            IF ( vmax>Rmax(3) ) THEN
               Rmax(3) = vmax
               IF ( Ninfo(3)==0 ) Lmax(3) = Knt
            ENDIF
!
!        Compare condition numbers for invariant subspace
!        without taking its condition number into account
!
            IF ( sepin<=v .AND. septmp<=v ) THEN
               vmax = ONE
            ELSEIF ( eps*sepin>septmp ) THEN
               vmax = ONE/eps
            ELSEIF ( sepin>septmp ) THEN
               vmax = sepin/septmp
            ELSEIF ( sepin<eps*septmp ) THEN
               vmax = ONE/eps
            ELSEIF ( sepin<septmp ) THEN
               vmax = septmp/sepin
            ELSE
               vmax = ONE
            ENDIF
            IF ( vmax>Rmax(3) ) THEN
               Rmax(3) = vmax
               IF ( Ninfo(3)==0 ) Lmax(3) = Knt
            ENDIF
!
!        Compute eigenvalue condition number only and compare
!        Update Q
!
            vmax = ZERO
            CALL ZLACPY('F',n,n,tsav1,LDT,ttmp,LDT)
            CALL ZLACPY('F',n,n,qsav,LDT,qtmp,LDT)
            septmp = -ONE
            stmp = -ONE
            CALL ZTRSEN('E','V',select,n,ttmp,LDT,qtmp,LDT,wtmp,m,stmp, &
     &                  septmp,work,LWORK,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            IF ( s/=stmp ) vmax = ONE/eps
            IF ( -ONE/=septmp ) vmax = ONE/eps
            DO i = 1 , n
               DO j = 1 , n
                  IF ( ttmp(i,j)/=t(i,j) ) vmax = ONE/eps
                  IF ( qtmp(i,j)/=q(i,j) ) vmax = ONE/eps
               ENDDO
            ENDDO
!
!        Compute invariant subspace condition number only and compare
!        Update Q
!
            CALL ZLACPY('F',n,n,tsav1,LDT,ttmp,LDT)
            CALL ZLACPY('F',n,n,qsav,LDT,qtmp,LDT)
            septmp = -ONE
            stmp = -ONE
            CALL ZTRSEN('V','V',select,n,ttmp,LDT,qtmp,LDT,wtmp,m,stmp, &
     &                  septmp,work,LWORK,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            IF ( -ONE/=stmp ) vmax = ONE/eps
            IF ( sep/=septmp ) vmax = ONE/eps
            DO i = 1 , n
               DO j = 1 , n
                  IF ( ttmp(i,j)/=t(i,j) ) vmax = ONE/eps
                  IF ( qtmp(i,j)/=q(i,j) ) vmax = ONE/eps
               ENDDO
            ENDDO
!
!        Compute eigenvalue condition number only and compare
!        Do not update Q
!
            CALL ZLACPY('F',n,n,tsav1,LDT,ttmp,LDT)
            CALL ZLACPY('F',n,n,qsav,LDT,qtmp,LDT)
            septmp = -ONE
            stmp = -ONE
            CALL ZTRSEN('E','N',select,n,ttmp,LDT,qtmp,LDT,wtmp,m,stmp, &
     &                  septmp,work,LWORK,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            IF ( s/=stmp ) vmax = ONE/eps
            IF ( -ONE/=septmp ) vmax = ONE/eps
            DO i = 1 , n
               DO j = 1 , n
                  IF ( ttmp(i,j)/=t(i,j) ) vmax = ONE/eps
                  IF ( qtmp(i,j)/=qsav(i,j) ) vmax = ONE/eps
               ENDDO
            ENDDO
!
!        Compute invariant subspace condition number only and compare
!        Do not update Q
!
            CALL ZLACPY('F',n,n,tsav1,LDT,ttmp,LDT)
            CALL ZLACPY('F',n,n,qsav,LDT,qtmp,LDT)
            septmp = -ONE
            stmp = -ONE
            CALL ZTRSEN('V','N',select,n,ttmp,LDT,qtmp,LDT,wtmp,m,stmp, &
     &                  septmp,work,LWORK,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
            IF ( -ONE/=stmp ) vmax = ONE/eps
            IF ( sep/=septmp ) vmax = ONE/eps
            DO i = 1 , n
               DO j = 1 , n
                  IF ( ttmp(i,j)/=t(i,j) ) vmax = ONE/eps
                  IF ( qtmp(i,j)/=qsav(i,j) ) vmax = ONE/eps
               ENDDO
            ENDDO
            IF ( vmax>Rmax(1) ) THEN
               Rmax(1) = vmax
               IF ( Ninfo(1)==0 ) Lmax(1) = Knt
            ENDIF
         ENDDO
      ENDDO
!
!     End of ZGET38
!
      END SUBROUTINE ZGET38
