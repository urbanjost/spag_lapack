!*==sget37.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b sget37
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sget37( RMAX, LMAX, NINFO, KNT, NIN )
!
!       .. Scalar Arguments ..
!       INTEGER            KNT, NIN
!       ..
!       .. Array Arguments ..
!       INTEGER            LMAX( 3 ), NINFO( 3 )
!       REAL               RMAX( 3 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sget37 tests STRSNA, a routine for estimating condition numbers of
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
!>          RMAX is REAL array, dimension (3)
!>          Value of the largest test ratio.
!>          RMAX(1) = largest ratio comparing different calls to STRSNA
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
!>          If SGEHRD returns INFO nonzero on example i, LMAX(1)=i
!>          If SHSEQR returns INFO nonzero on example i, LMAX(2)=i
!>          If STRSNA returns INFO nonzero on example i, LMAX(3)=i
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER array, dimension (3)
!>          NINFO(1) = No. of times SGEHRD returned INFO nonzero
!>          NINFO(2) = No. of times SHSEQR returned INFO nonzero
!>          NINFO(3) = No. of times STRSNA returned INFO nonzero
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE sget37(Rmax,Lmax,Ninfo,Knt,Nin)
      IMPLICIT NONE
!*--SGET3794
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
      REAL Rmax(3)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,TWO=2.0E0)
      REAL EPSIN
      PARAMETER (EPSIN=5.9605E-8)
      INTEGER LDT , LWORK
      PARAMETER (LDT=20,LWORK=2*LDT*(10+LDT))
!     ..
!     .. Local Scalars ..
      INTEGER i , icmp , ifnd , info , iscl , j , kmin , m , n
      REAL bignum , eps , smlnum , tnrm , tol , tolin , v , vimin ,     &
     &     vmax , vmul , vrmin
!     ..
!     .. Local Arrays ..
      LOGICAL select(LDT)
      INTEGER iwork(2*LDT) , lcmp(3)
      REAL dum(1) , le(LDT,LDT) , re(LDT,LDT) , s(LDT) , sep(LDT) ,     &
     &     sepin(LDT) , septmp(LDT) , sin(LDT) , stmp(LDT) , t(LDT,LDT) &
     &     , tmp(LDT,LDT) , val(3) , wi(LDT) , wiin(LDT) , witmp(LDT) , &
     &     work(LWORK) , wr(LDT) , wrin(LDT) , wrtmp(LDT)
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE
      EXTERNAL SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SGEHRD , SHSEQR , SLABAD , SLACPY , SSCAL ,  STREVC , STRSNA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , REAL , SQRT
!     ..
!     .. Executable Statements ..
!
      eps = SLAMCH('P')
      smlnum = SLAMCH('S')/eps
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
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
!
      val(1) = SQRT(smlnum)
      val(2) = ONE
      val(3) = SQRT(bignum)
      DO
!
!     Read input data until N=0.  Assume input eigenvalues are sorted
!     lexicographically (increasing by real part, then decreasing by
!     imaginary part)
!
         READ (Nin,FMT=*) n
         IF ( n==0 ) RETURN
         DO i = 1 , n
            READ (Nin,FMT=*) (tmp(i,j),j=1,n)
         ENDDO
         DO i = 1 , n
            READ (Nin,FMT=*) wrin(i) , wiin(i) , sin(i) , sepin(i)
         ENDDO
         tnrm = SLANGE('M',n,n,tmp,LDT,work)
!
!     Begin test
!
         DO iscl = 1 , 3
!
!        Scale input matrix
!
            Knt = Knt + 1
            CALL SLACPY('F',n,n,tmp,LDT,t,LDT)
            vmul = val(iscl)
            DO i = 1 , n
               CALL SSCAL(n,vmul,t(1,i),1)
            ENDDO
            IF ( tnrm==ZERO ) vmul = ONE
!
!        Compute eigenvalues and eigenvectors
!
            CALL SGEHRD(n,1,n,t,LDT,work(1),work(n+1),LWORK-n,info)
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
            CALL SHSEQR('S','N',n,1,n,t,LDT,wr,wi,dum,1,work,LWORK,info)
            IF ( info/=0 ) THEN
               Lmax(2) = Knt
               Ninfo(2) = Ninfo(2) + 1
               CYCLE
            ENDIF
!
!        Compute eigenvectors
!
            CALL STREVC('Both','All',select,n,t,LDT,le,LDT,re,LDT,n,m,  &
     &                  work,info)
!
!        Compute condition numbers
!
            CALL STRSNA('Both','All',select,n,t,LDT,le,LDT,re,LDT,s,sep,&
     &                  n,m,work,n,iwork,info)
            IF ( info/=0 ) THEN
               Lmax(3) = Knt
               Ninfo(3) = Ninfo(3) + 1
               CYCLE
            ENDIF
!
!        Sort eigenvalues and condition numbers lexicographically
!        to compare with inputs
!
            CALL SCOPY(n,wr,1,wrtmp,1)
            CALL SCOPY(n,wi,1,witmp,1)
            CALL SCOPY(n,s,1,stmp,1)
            CALL SCOPY(n,sep,1,septmp,1)
            CALL SSCAL(n,ONE/vmul,septmp,1)
            DO i = 1 , n - 1
               kmin = i
               vrmin = wrtmp(i)
               vimin = witmp(i)
               DO j = i + 1 , n
                  IF ( wrtmp(j)<vrmin ) THEN
                     kmin = j
                     vrmin = wrtmp(j)
                     vimin = witmp(j)
                  ENDIF
               ENDDO
               wrtmp(kmin) = wrtmp(i)
               witmp(kmin) = witmp(i)
               wrtmp(i) = vrmin
               witmp(i) = vimin
               vrmin = stmp(kmin)
               stmp(kmin) = stmp(i)
               stmp(i) = vrmin
               vrmin = septmp(kmin)
               septmp(kmin) = septmp(i)
               septmp(i) = vrmin
            ENDDO
!
!        Compare condition numbers for eigenvalues
!        taking their condition numbers into account
!
            v = MAX(TWO*REAL(n)*eps*tnrm,smlnum)
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
               IF ( sin(i)<=REAL(2*n)*eps .AND. stmp(i)<=REAL(2*n)*eps )&
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
            CALL SCOPY(n,dum,0,stmp,1)
            CALL SCOPY(n,dum,0,septmp,1)
            CALL STRSNA('Eigcond','All',select,n,t,LDT,le,LDT,re,LDT,   &
     &                  stmp,septmp,n,m,work,n,iwork,info)
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
            CALL SCOPY(n,dum,0,stmp,1)
            CALL SCOPY(n,dum,0,septmp,1)
            CALL STRSNA('Veccond','All',select,n,t,LDT,le,LDT,re,LDT,   &
     &                  stmp,septmp,n,m,work,n,iwork,info)
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
            CALL SCOPY(n,dum,0,stmp,1)
            CALL SCOPY(n,dum,0,septmp,1)
            CALL STRSNA('Bothcond','Some',select,n,t,LDT,le,LDT,re,LDT, &
     &                  stmp,septmp,n,m,work,n,iwork,info)
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
            CALL SCOPY(n,dum,0,stmp,1)
            CALL SCOPY(n,dum,0,septmp,1)
            CALL STRSNA('Eigcond','Some',select,n,t,LDT,le,LDT,re,LDT,  &
     &                  stmp,septmp,n,m,work,n,iwork,info)
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
            CALL SCOPY(n,dum,0,stmp,1)
            CALL SCOPY(n,dum,0,septmp,1)
            CALL STRSNA('Veccond','Some',select,n,t,LDT,le,LDT,re,LDT,  &
     &                  stmp,septmp,n,m,work,n,iwork,info)
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
!        Select first real and first complex eigenvalue
!
            IF ( wi(1)==ZERO ) THEN
               lcmp(1) = 1
               ifnd = 0
               DO i = 2 , n
                  IF ( ifnd==1 .OR. wi(i)==ZERO ) THEN
                     select(i) = .FALSE.
                  ELSE
                     ifnd = 1
                     lcmp(2) = i
                     lcmp(3) = i + 1
                     CALL SCOPY(n,re(1,i),1,re(1,2),1)
                     CALL SCOPY(n,re(1,i+1),1,re(1,3),1)
                     CALL SCOPY(n,le(1,i),1,le(1,2),1)
                     CALL SCOPY(n,le(1,i+1),1,le(1,3),1)
                  ENDIF
               ENDDO
               IF ( ifnd==0 ) THEN
                  icmp = 1
               ELSE
                  icmp = 3
               ENDIF
            ELSE
               lcmp(1) = 1
               lcmp(2) = 2
               ifnd = 0
               DO i = 3 , n
                  IF ( ifnd==1 .OR. wi(i)/=ZERO ) THEN
                     select(i) = .FALSE.
                  ELSE
                     lcmp(3) = i
                     ifnd = 1
                     CALL SCOPY(n,re(1,i),1,re(1,3),1)
                     CALL SCOPY(n,le(1,i),1,le(1,3),1)
                  ENDIF
               ENDDO
               IF ( ifnd==0 ) THEN
                  icmp = 2
               ELSE
                  icmp = 3
               ENDIF
            ENDIF
!
!        Compute all selected condition numbers
!
            CALL SCOPY(icmp,dum,0,stmp,1)
            CALL SCOPY(icmp,dum,0,septmp,1)
            CALL STRSNA('Bothcond','Some',select,n,t,LDT,le,LDT,re,LDT, &
     &                  stmp,septmp,n,m,work,n,iwork,info)
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
            CALL SCOPY(icmp,dum,0,stmp,1)
            CALL SCOPY(icmp,dum,0,septmp,1)
            CALL STRSNA('Eigcond','Some',select,n,t,LDT,le,LDT,re,LDT,  &
     &                  stmp,septmp,n,m,work,n,iwork,info)
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
            CALL SCOPY(icmp,dum,0,stmp,1)
            CALL SCOPY(icmp,dum,0,septmp,1)
            CALL STRSNA('Veccond','Some',select,n,t,LDT,le,LDT,re,LDT,  &
     &                  stmp,septmp,n,m,work,n,iwork,info)
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
!     End of sget37
!
      END SUBROUTINE sget37
