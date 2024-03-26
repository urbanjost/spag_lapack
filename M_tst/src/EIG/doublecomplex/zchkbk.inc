!*==zchkbk.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b zchkbk
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zchkbk( NIN, NOUT )
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
!> zchkbk tests ZGEBAK, a routine for backward transformation of
!> the computed right or left eigenvectors if the original matrix
!> was preprocessed by balance subroutine ZGEBAL.
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
      SUBROUTINE zchkbk(Nin,Nout)
      IMPLICIT NONE
!*--ZCHKBK59
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
! ======================================================================
!
!     .. Parameters ..
      INTEGER LDE
      PARAMETER (LDE=20)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ihi , ilo , info , j , knt , n , ninfo
      DOUBLE PRECISION eps , rmax , safmin , vmax , x
      COMPLEX*16 cdum
!     ..
!     .. Local Arrays ..
      INTEGER lmax(2)
      DOUBLE PRECISION scale(LDE)
      COMPLEX*16 e(LDE,LDE) , ein(LDE,LDE)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEBAK
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
      ninfo = 0
      knt = 0
      rmax = ZERO
      eps = DLAMCH('E')
      safmin = DLAMCH('S')
      DO
!
!
         READ (Nin,FMT=*) n , ilo , ihi
         IF ( n==0 ) THEN
!
!
            WRITE (Nout,FMT=99001)
99001       FORMAT (1X,'.. test output of ZGEBAK .. ')
!
            WRITE (Nout,FMT=99002) rmax
99002       FORMAT (1X,'value of largest test error             = ',    &
     &              D12.3)
            WRITE (Nout,FMT=99003) lmax(1)
99003       FORMAT (1X,'example number where info is not zero   = ',I4)
            WRITE (Nout,FMT=99004) lmax(2)
99004       FORMAT (1X,'example number having largest error     = ',I4)
            WRITE (Nout,FMT=99005) ninfo
99005       FORMAT (1X,'number of examples where info is not 0  = ',I4)
            WRITE (Nout,FMT=99006) knt
99006       FORMAT (1X,'total number of examples tested         = ',I4)
            EXIT
         ELSE
!
            READ (Nin,FMT=*) (scale(i),i=1,n)
            DO i = 1 , n
               READ (Nin,FMT=*) (e(i,j),j=1,n)
            ENDDO
!
            DO i = 1 , n
               READ (Nin,FMT=*) (ein(i,j),j=1,n)
            ENDDO
!
            knt = knt + 1
            CALL ZGEBAK('B','R',n,ilo,ihi,scale,n,e,LDE,info)
!
            IF ( info/=0 ) THEN
               ninfo = ninfo + 1
               lmax(1) = knt
            ENDIF
!
            vmax = ZERO
            DO i = 1 , n
               DO j = 1 , n
                  x = CABS1(e(i,j)-ein(i,j))/eps
                  IF ( CABS1(e(i,j))>safmin ) x = x/CABS1(e(i,j))
                  vmax = MAX(vmax,x)
               ENDDO
            ENDDO
!
            IF ( vmax>rmax ) THEN
               lmax(2) = knt
               rmax = vmax
!
            ENDIF
         ENDIF
      ENDDO
!
!
!     End of zchkbk
!
      END SUBROUTINE zchkbk
