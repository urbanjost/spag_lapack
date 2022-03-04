!*==zchkgk.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZCHKGK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKGK( NIN, NOUT )
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
!> ZCHKGK tests ZGGBAK, a routine for backward balancing  of
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZCHKGK(Nin,Nout)
      IMPLICIT NONE
!*--ZCHKGK58
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
      INTEGER LDE , LDF , LDWORK , LRWORK
      PARAMETER (LDE=50,LDF=50,LDWORK=50,LRWORK=6*50)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , ihi , ilo , info , j , knt , m , n , ninfo
      DOUBLE PRECISION anorm , bnorm , eps , rmax , vmax
      COMPLEX*16 cdum
!     ..
!     .. Local Arrays ..
      INTEGER lmax(4)
      DOUBLE PRECISION lscale(LDA) , rscale(LDA) , rwork(LRWORK)
      COMPLEX*16 a(LDA,LDA) , af(LDA,LDA) , b(LDB,LDB) , bf(LDB,LDB) ,  &
     &           e(LDE,LDE) , f(LDF,LDF) , vl(LDVL,LDVL) ,              &
     &           vlf(LDVL,LDVL) , vr(LDVR,LDVR) , vrf(LDVR,LDVR) ,      &
     &           work(LDWORK,LDWORK)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL DLAMCH , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMM , ZGGBAK , ZGGBAL , ZLACPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DIMAG , MAX
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1(cdum) = ABS(DBLE(cdum)) + ABS(DIMAG(cdum))
!     ..
!     .. Executable Statements ..
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
99001       FORMAT (1X,'.. test output of ZGGBAK .. ')
!
            WRITE (Nout,FMT=99002) rmax
99002       FORMAT (' value of largest test error                  =',  &
     &              D12.3)
            WRITE (Nout,FMT=99003) lmax(1)
99003       FORMAT (' example number where ZGGBAL info is not 0    =',  &
     &              I4)
            WRITE (Nout,FMT=99004) lmax(2)
99004       FORMAT (' example number where ZGGBAK(L) info is not 0 =',  &
     &              I4)
            WRITE (Nout,FMT=99005) lmax(3)
99005       FORMAT (' example number where ZGGBAK(R) info is not 0 =',  &
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
            anorm = ZLANGE('M',n,n,a,LDA,rwork)
            bnorm = ZLANGE('M',n,n,b,LDB,rwork)
!
            CALL ZLACPY('FULL',n,n,a,LDA,af,LDA)
            CALL ZLACPY('FULL',n,n,b,LDB,bf,LDB)
!
            CALL ZGGBAL('B',n,a,LDA,b,LDB,ilo,ihi,lscale,rscale,rwork,  &
     &                  info)
            IF ( info/=0 ) THEN
               ninfo = ninfo + 1
               lmax(1) = knt
            ENDIF
!
            CALL ZLACPY('FULL',n,m,vl,LDVL,vlf,LDVL)
            CALL ZLACPY('FULL',n,m,vr,LDVR,vrf,LDVR)
!
            CALL ZGGBAK('B','L',n,ilo,ihi,lscale,rscale,m,vl,LDVL,info)
            IF ( info/=0 ) THEN
               ninfo = ninfo + 1
               lmax(2) = knt
            ENDIF
!
            CALL ZGGBAK('B','R',n,ilo,ihi,lscale,rscale,m,vr,LDVR,info)
            IF ( info/=0 ) THEN
               ninfo = ninfo + 1
               lmax(3) = knt
            ENDIF
!
!     Test of ZGGBAK
!
!     Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
!     where tilde(A) denotes the transformed matrix.
!
            CALL ZGEMM('N','N',n,m,n,CONE,af,LDA,vr,LDVR,CZERO,work,    &
     &                 LDWORK)
            CALL ZGEMM('C','N',m,m,n,CONE,vl,LDVL,work,LDWORK,CZERO,e,  &
     &                 LDE)
!
            CALL ZGEMM('N','N',n,m,n,CONE,a,LDA,vrf,LDVR,CZERO,work,    &
     &                 LDWORK)
            CALL ZGEMM('C','N',m,m,n,CONE,vlf,LDVL,work,LDWORK,CZERO,f, &
     &                 LDF)
!
            vmax = ZERO
            DO j = 1 , m
               DO i = 1 , m
                  vmax = MAX(vmax,CABS1(e(i,j)-f(i,j)))
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
            CALL ZGEMM('N','N',n,m,n,CONE,bf,LDB,vr,LDVR,CZERO,work,    &
     &                 LDWORK)
            CALL ZGEMM('C','N',m,m,n,CONE,vl,LDVL,work,LDWORK,CZERO,e,  &
     &                 LDE)
!
            CALL ZGEMM('n','n',n,m,n,CONE,b,LDB,vrf,LDVR,CZERO,work,    &
     &                 LDWORK)
            CALL ZGEMM('C','N',m,m,n,CONE,vlf,LDVL,work,LDWORK,CZERO,f, &
     &                 LDF)
!
            vmax = ZERO
            DO j = 1 , m
               DO i = 1 , m
                  vmax = MAX(vmax,CABS1(e(i,j)-f(i,j)))
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
!     End of ZCHKGK
!
      END SUBROUTINE ZCHKGK
