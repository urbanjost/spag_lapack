!*==dchkgk.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DCHKGK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCHKGK( NIN, NOUT )
!
!       .. Scalar Arguments ..
!       INTEGER            NIN, NOUT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCHKGK tests DGGBAK, a routine for backward balancing  of
!> a matrix pair (A, B).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NIN
!> \verbatim
!>          NIN is INTEGER
!>          The logical unit number for input.  NIN > 0.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The logical unit number for output.  NOUT > 0.
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
      SUBROUTINE DCHKGK(Nin,Nout)
      IMPLICIT NONE
!*--DCHKGK58
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Nin , Nout
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER LDA , LDB , LDVL , LDVR
      PARAMETER (LDA=50,LDB=50,LDVL=50,LDVR=50)
      INTEGER LDE , LDF , LDWORK
      PARAMETER (LDE=50,LDF=50,LDWORK=50)
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ihi , ilo , info , j , knt , m , n , ninfo
      DOUBLE PRECISION anorm , bnorm , eps , rmax , vmax
!     ..
!     .. Local Arrays ..
      INTEGER lmax(4)
      DOUBLE PRECISION a(LDA,LDA) , af(LDA,LDA) , b(LDB,LDB) ,          &
     &                 bf(LDB,LDB) , e(LDE,LDE) , f(LDF,LDF) ,          &
     &                 lscale(LDA) , rscale(LDA) , vl(LDVL,LDVL) ,      &
     &                 vlf(LDVL,LDVL) , vr(LDVR,LDVR) , vrf(LDVR,LDVR) ,&
     &                 work(LDWORK,LDWORK)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE
      EXTERNAL DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DGGBAK , DGGBAL , DLACPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
!     Initialization
!
      lmax(1) = 0
      lmax(2) = 0
      lmax(3) = 0
      lmax(4) = 0
      ninfo = 0
      knt = 0
      rmax = ZERO
!
      eps = DLAMCH('Precision')
      DO
!
         READ (Nin,FMT=*) n , m
         IF ( n==0 ) THEN
!
!
            WRITE (Nout,FMT=99001)
99001       FORMAT (1X,'.. test output of DGGBAK .. ')
!
            WRITE (Nout,FMT=99002) rmax
99002       FORMAT (' value of largest test error                  =',  &
     &              D12.3)
            WRITE (Nout,FMT=99003) lmax(1)
99003       FORMAT (' example number where DGGBAL info is not 0    =',  &
     &              I4)
            WRITE (Nout,FMT=99004) lmax(2)
99004       FORMAT (' example number where DGGBAK(L) info is not 0 =',  &
     &              I4)
            WRITE (Nout,FMT=99005) lmax(3)
99005       FORMAT (' example number where DGGBAK(R) info is not 0 =',  &
     &              I4)
            WRITE (Nout,FMT=99006) lmax(4)
99006       FORMAT (' example number having largest error          =',  &
     &              I4)
            WRITE (Nout,FMT=99007) ninfo
99007       FORMAT (' number of examples where info is not 0       =',  &
     &              I4)
            WRITE (Nout,FMT=99008) knt
99008       FORMAT (' total number of examples tested              =',  &
     &              I4)
            EXIT
         ELSE
!
            DO i = 1 , n
               READ (Nin,FMT=*) (a(i,j),j=1,n)
            ENDDO
!
            DO i = 1 , n
               READ (Nin,FMT=*) (b(i,j),j=1,n)
            ENDDO
!
            DO i = 1 , n
               READ (Nin,FMT=*) (vl(i,j),j=1,m)
            ENDDO
!
            DO i = 1 , n
               READ (Nin,FMT=*) (vr(i,j),j=1,m)
            ENDDO
!
            knt = knt + 1
!
            anorm = DLANGE('M',n,n,a,LDA,work)
            bnorm = DLANGE('M',n,n,b,LDB,work)
!
            CALL DLACPY('FULL',n,n,a,LDA,af,LDA)
            CALL DLACPY('FULL',n,n,b,LDB,bf,LDB)
!
            CALL DGGBAL('B',n,a,LDA,b,LDB,ilo,ihi,lscale,rscale,work,   &
     &                  info)
            IF ( info/=0 ) THEN
               ninfo = ninfo + 1
               lmax(1) = knt
            ENDIF
!
            CALL DLACPY('FULL',n,m,vl,LDVL,vlf,LDVL)
            CALL DLACPY('FULL',n,m,vr,LDVR,vrf,LDVR)
!
            CALL DGGBAK('B','L',n,ilo,ihi,lscale,rscale,m,vl,LDVL,info)
            IF ( info/=0 ) THEN
               ninfo = ninfo + 1
               lmax(2) = knt
            ENDIF
!
            CALL DGGBAK('B','R',n,ilo,ihi,lscale,rscale,m,vr,LDVR,info)
            IF ( info/=0 ) THEN
               ninfo = ninfo + 1
               lmax(3) = knt
            ENDIF
!
!     Test of DGGBAK
!
!     Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
!     where tilde(A) denotes the transformed matrix.
!
            CALL DGEMM('N','N',n,m,n,ONE,af,LDA,vr,LDVR,ZERO,work,      &
     &                 LDWORK)
            CALL DGEMM('T','N',m,m,n,ONE,vl,LDVL,work,LDWORK,ZERO,e,LDE)
!
            CALL DGEMM('N','N',n,m,n,ONE,a,LDA,vrf,LDVR,ZERO,work,      &
     &                 LDWORK)
            CALL DGEMM('T','N',m,m,n,ONE,vlf,LDVL,work,LDWORK,ZERO,f,   &
     &                 LDF)
!
            vmax = ZERO
            DO j = 1 , m
               DO i = 1 , m
                  vmax = MAX(vmax,ABS(e(i,j)-f(i,j)))
               ENDDO
            ENDDO
            vmax = vmax/(eps*MAX(anorm,bnorm))
            IF ( vmax>rmax ) THEN
               lmax(4) = knt
               rmax = vmax
            ENDIF
!
!     Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR
!
            CALL DGEMM('N','N',n,m,n,ONE,bf,LDB,vr,LDVR,ZERO,work,      &
     &                 LDWORK)
            CALL DGEMM('T','N',m,m,n,ONE,vl,LDVL,work,LDWORK,ZERO,e,LDE)
!
            CALL DGEMM('N','N',n,m,n,ONE,b,LDB,vrf,LDVR,ZERO,work,      &
     &                 LDWORK)
            CALL DGEMM('T','N',m,m,n,ONE,vlf,LDVL,work,LDWORK,ZERO,f,   &
     &                 LDF)
!
            vmax = ZERO
            DO j = 1 , m
               DO i = 1 , m
                  vmax = MAX(vmax,ABS(e(i,j)-f(i,j)))
               ENDDO
            ENDDO
            vmax = vmax/(eps*MAX(anorm,bnorm))
            IF ( vmax>rmax ) THEN
               lmax(4) = knt
               rmax = vmax
!
            ENDIF
         ENDIF
      ENDDO
!
!
!     End of DCHKGK
!
      END SUBROUTINE DCHKGK
