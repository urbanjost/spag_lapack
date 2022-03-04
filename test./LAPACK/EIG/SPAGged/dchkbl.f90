!*==dchkbl.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DCHKBL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCHKBL( NIN, NOUT )
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
!> DCHKBL tests DGEBAL, a routine for balancing a general real
!> matrix and isolating some of its eigenvalues.
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
      SUBROUTINE DCHKBL(Nin,Nout)
      IMPLICIT NONE
!*--DCHKBL58
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
      INTEGER LDA
      PARAMETER (LDA=20)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ihi , ihiin , ilo , iloin , info , j , knt , n , ninfo
      DOUBLE PRECISION anorm , meps , rmax , sfmin , temp , vmax
!     ..
!     .. Local Arrays ..
      INTEGER lmax(3)
      DOUBLE PRECISION a(LDA,LDA) , ain(LDA,LDA) , dummy(1) , scale(LDA)&
     &                 , scalin(LDA)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE
      EXTERNAL DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEBAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
      lmax(1) = 0
      lmax(2) = 0
      lmax(3) = 0
      ninfo = 0
      knt = 0
      rmax = ZERO
      vmax = ZERO
      sfmin = DLAMCH('S')
      meps = DLAMCH('E')
      DO
!
!
         READ (Nin,FMT=*) n
         IF ( n==0 ) THEN
!
!
            WRITE (Nout,FMT=99001)
99001       FORMAT (1X,'.. test output of DGEBAL .. ')
!
            WRITE (Nout,FMT=99002) rmax
99002       FORMAT (1X,'value of largest test error            = ',     &
     &              D12.3)
            WRITE (Nout,FMT=99003) lmax(1)
99003       FORMAT (1X,'example number where info is not zero  = ',I4)
            WRITE (Nout,FMT=99004) lmax(2)
99004       FORMAT (1X,'example number where ILO or IHI wrong  = ',I4)
            WRITE (Nout,FMT=99005) lmax(3)
99005       FORMAT (1X,'example number having largest error    = ',I4)
            WRITE (Nout,FMT=99006) ninfo
99006       FORMAT (1X,'number of examples where info is not 0 = ',I4)
            WRITE (Nout,FMT=99007) knt
99007       FORMAT (1X,'total number of examples tested        = ',I4)
            EXIT
         ELSE
            DO i = 1 , n
               READ (Nin,FMT=*) (a(i,j),j=1,n)
            ENDDO
!
            READ (Nin,FMT=*) iloin , ihiin
            DO i = 1 , n
               READ (Nin,FMT=*) (ain(i,j),j=1,n)
            ENDDO
            READ (Nin,FMT=*) (scalin(i),i=1,n)
!
            anorm = DLANGE('M',n,n,a,LDA,dummy)
            knt = knt + 1
!
            CALL DGEBAL('B',n,a,LDA,ilo,ihi,scale,info)
!
            IF ( info/=0 ) THEN
               ninfo = ninfo + 1
               lmax(1) = knt
            ENDIF
!
            IF ( ilo/=iloin .OR. ihi/=ihiin ) THEN
               ninfo = ninfo + 1
               lmax(2) = knt
            ENDIF
!
            DO i = 1 , n
               DO j = 1 , n
                  temp = MAX(a(i,j),ain(i,j))
                  temp = MAX(temp,sfmin)
                  vmax = MAX(vmax,ABS(a(i,j)-ain(i,j))/temp)
               ENDDO
            ENDDO
!
            DO i = 1 , n
               temp = MAX(scale(i),scalin(i))
               temp = MAX(temp,sfmin)
               vmax = MAX(vmax,ABS(scale(i)-scalin(i))/temp)
            ENDDO
!
!
            IF ( vmax>rmax ) THEN
               lmax(3) = knt
               rmax = vmax
!
            ENDIF
         ENDIF
      ENDDO
!
!
!     End of DCHKBL
!
      END SUBROUTINE DCHKBL
