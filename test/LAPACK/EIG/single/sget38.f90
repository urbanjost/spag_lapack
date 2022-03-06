!*==sget38.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b SGET38
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGET38( RMAX, LMAX, NINFO, KNT, NIN )
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
!> SGET38 tests STRSEN, a routine for estimating condition numbers of a
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
!>          RMAX is REAL array, dimension (3)
!>          Values of the largest test ratios.
!>          RMAX(1) = largest residuals from SHST01 or comparing
!>                    different calls to STRSEN
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
!>          If STRSEN returns INFO nonzero on example i, LMAX(3)=i
!> \endverbatim
!>
!> \param[out] NINFO
!> \verbatim
!>          NINFO is INTEGER array, dimension (3)
!>          NINFO(1) = No. of times SGEHRD returned INFO nonzero
!>          NINFO(2) = No. of times SHSEQR returned INFO nonzero
!>          NINFO(3) = No. of times STRSEN returned INFO nonzero
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE SGET38(Rmax,Lmax,Ninfo,Knt,Nin)
      IMPLICIT NONE
!*--SGET3895
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
      INTEGER LIWORK
      PARAMETER (LIWORK=LDT*LDT)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , iscl , itmp , j , kmin , m , n , ndim
      REAL bignum , eps , s , sep , sepin , septmp , sin , smlnum ,     &
     &     stmp , tnrm , tol , tolin , v , vimin , vmax , vmul , vrmin
!     ..
!     .. Local Arrays ..
      LOGICAL select(LDT)
      INTEGER ipnt(LDT) , iselec(LDT) , iwork(LIWORK)
      REAL q(LDT,LDT) , qsav(LDT,LDT) , qtmp(LDT,LDT) , result(2) ,     &
     &     t(LDT,LDT) , tmp(LDT,LDT) , tsav(LDT,LDT) , tsav1(LDT,LDT) , &
     &     ttmp(LDT,LDT) , val(3) , wi(LDT) , witmp(LDT) , work(LWORK) ,&
     &     wr(LDT) , wrtmp(LDT)
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE
      EXTERNAL SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SGEHRD , SHSEQR , SHST01 , SLABAD , SLACPY ,     &
     &         SORGHR , SSCAL , STRSEN
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
      val(3) = SQRT(SQRT(bignum))
      DO
!
!     Read input data until N=0.  Assume input eigenvalues are sorted
!     lexicographically (increasing by real part, then decreasing by
!     imaginary part)
!
         READ (Nin,FMT=*) n , ndim
         IF ( n==0 ) RETURN
         READ (Nin,FMT=*) (iselec(i),i=1,ndim)
         DO i = 1 , n
            READ (Nin,FMT=*) (tmp(i,j),j=1,n)
         ENDDO
         READ (Nin,FMT=*) sin , sepin
!
         tnrm = SLANGE('M',n,n,tmp,LDT,work)
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
            CALL SLACPY('F',n,n,t,LDT,tsav,LDT)
!
!        Compute Schur form
!
            CALL SGEHRD(n,1,n,t,LDT,work(1),work(n+1),LWORK-n,info)
            IF ( info/=0 ) THEN
               Lmax(1) = Knt
               Ninfo(1) = Ninfo(1) + 1
               CYCLE
            ENDIF
!
!        Generate orthogonal matrix
!
            CALL SLACPY('L',n,n,t,LDT,q,LDT)
            CALL SORGHR(n,1,n,q,LDT,work(1),work(n+1),LWORK-n,info)
!
!        Compute Schur form
!
            CALL SHSEQR('S','V',n,1,n,t,LDT,wr,wi,q,LDT,work,LWORK,info)
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
            CALL SCOPY(n,wr,1,wrtmp,1)
            CALL SCOPY(n,wi,1,witmp,1)
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
            CALL SLACPY('F',n,n,q,LDT,qsav,LDT)
            CALL SLACPY('F',n,n,t,LDT,tsav1,LDT)
            CALL STRSEN('B','V',select,n,t,LDT,q,LDT,wrtmp,witmp,m,s,   &
     &                  sep,work,LWORK,iwork,LIWORK,info)
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
            CALL SHST01(n,1,n,tsav,LDT,t,LDT,q,LDT,work,LWORK,result)
            vmax = MAX(result(1),result(2))
            IF ( vmax>Rmax(1) ) THEN
               Rmax(1) = vmax
               IF ( Ninfo(1)==0 ) Lmax(1) = Knt
            ENDIF
!
!        Compare condition number for eigenvalue cluster
!        taking its condition number into account
!
            v = MAX(TWO*REAL(n)*eps*tnrm,smlnum)
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
            IF ( sin<=REAL(2*n)*eps .AND. stmp<=REAL(2*n)*eps ) THEN
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
            CALL SLACPY('F',n,n,tsav1,LDT,ttmp,LDT)
            CALL SLACPY('F',n,n,qsav,LDT,qtmp,LDT)
            septmp = -ONE
            stmp = -ONE
            CALL STRSEN('E','V',select,n,ttmp,LDT,qtmp,LDT,wrtmp,witmp, &
     &                  m,stmp,septmp,work,LWORK,iwork,LIWORK,info)
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
            CALL SLACPY('F',n,n,tsav1,LDT,ttmp,LDT)
            CALL SLACPY('F',n,n,qsav,LDT,qtmp,LDT)
            septmp = -ONE
            stmp = -ONE
            CALL STRSEN('V','V',select,n,ttmp,LDT,qtmp,LDT,wrtmp,witmp, &
     &                  m,stmp,septmp,work,LWORK,iwork,LIWORK,info)
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
            CALL SLACPY('F',n,n,tsav1,LDT,ttmp,LDT)
            CALL SLACPY('F',n,n,qsav,LDT,qtmp,LDT)
            septmp = -ONE
            stmp = -ONE
            CALL STRSEN('E','N',select,n,ttmp,LDT,qtmp,LDT,wrtmp,witmp, &
     &                  m,stmp,septmp,work,LWORK,iwork,LIWORK,info)
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
            CALL SLACPY('F',n,n,tsav1,LDT,ttmp,LDT)
            CALL SLACPY('F',n,n,qsav,LDT,qtmp,LDT)
            septmp = -ONE
            stmp = -ONE
            CALL STRSEN('V','N',select,n,ttmp,LDT,qtmp,LDT,wrtmp,witmp, &
     &                  m,stmp,septmp,work,LWORK,iwork,LIWORK,info)
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
!     End of SGET38
!
      END SUBROUTINE SGET38
