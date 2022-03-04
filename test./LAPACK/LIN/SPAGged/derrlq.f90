!*==derrlq.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DERRLQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DERRLQ( PATH, NUNIT )
!
!       .. Scalar Arguments ..
!       CHARACTER*3        PATH
!       INTEGER            NUNIT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DERRLQ tests the error exits for the DOUBLE PRECISION routines
!> that use the LQ decomposition of a general matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3
!>          The LAPACK path name for the routines to be tested.
!> \endverbatim
!>
!> \param[in] NUNIT
!> \verbatim
!>          NUNIT is INTEGER
!>          The unit number for output.
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DERRLQ(Path,Nunit)
      IMPLICIT NONE
!*--DERRLQ59
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*3 Path
      INTEGER Nunit
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=2)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(NMAX) &
     &                 , x(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , DGELQ2 , DGELQF , DGELQS , DORGL2 ,    &
     &         DORGLQ , DORML2 , DORMLQ
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NOUt
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUt , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE
!     ..
!     .. Executable Statements ..
!
      NOUt = Nunit
      WRITE (NOUt,FMT=*)
!
!     Set the variables to innocuous values.
!
      DO j = 1 , NMAX
         DO i = 1 , NMAX
            a(i,j) = 1.D0/DBLE(i+j)
            af(i,j) = 1.D0/DBLE(i+j)
         ENDDO
         b(j) = 0.D0
         w(j) = 0.D0
         x(j) = 0.D0
      ENDDO
      OK = .TRUE.
!
!     Error exits for LQ factorization
!
!     DGELQF
!
      SRNamt = 'DGELQF'
      INFot = 1
      CALL DGELQF(-1,0,a,1,b,w,1,info)
      CALL CHKXER('DGELQF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGELQF(0,-1,a,1,b,w,1,info)
      CALL CHKXER('DGELQF',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGELQF(2,1,a,1,b,w,2,info)
      CALL CHKXER('DGELQF',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DGELQF(2,1,a,2,b,w,1,info)
      CALL CHKXER('DGELQF',INFot,NOUt,LERr,OK)
!
!     DGELQ2
!
      SRNamt = 'DGELQ2'
      INFot = 1
      CALL DGELQ2(-1,0,a,1,b,w,info)
      CALL CHKXER('DGELQ2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGELQ2(0,-1,a,1,b,w,info)
      CALL CHKXER('DGELQ2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DGELQ2(2,1,a,1,b,w,info)
      CALL CHKXER('DGELQ2',INFot,NOUt,LERr,OK)
!
!     DGELQS
!
      SRNamt = 'DGELQS'
      INFot = 1
      CALL DGELQS(-1,0,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('DGELQS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGELQS(0,-1,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('DGELQS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DGELQS(2,1,0,a,2,x,b,1,w,1,info)
      CALL CHKXER('DGELQS',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DGELQS(0,0,-1,a,1,x,b,1,w,1,info)
      CALL CHKXER('DGELQS',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DGELQS(2,2,0,a,1,x,b,2,w,1,info)
      CALL CHKXER('DGELQS',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DGELQS(1,2,0,a,1,x,b,1,w,1,info)
      CALL CHKXER('DGELQS',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DGELQS(1,1,2,a,1,x,b,1,w,1,info)
      CALL CHKXER('DGELQS',INFot,NOUt,LERr,OK)
!
!     DORGLQ
!
      SRNamt = 'DORGLQ'
      INFot = 1
      CALL DORGLQ(-1,0,0,a,1,x,w,1,info)
      CALL CHKXER('DORGLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORGLQ(0,-1,0,a,1,x,w,1,info)
      CALL CHKXER('DORGLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORGLQ(2,1,0,a,2,x,w,2,info)
      CALL CHKXER('DORGLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORGLQ(0,0,-1,a,1,x,w,1,info)
      CALL CHKXER('DORGLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORGLQ(1,1,2,a,1,x,w,1,info)
      CALL CHKXER('DORGLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORGLQ(2,2,0,a,1,x,w,2,info)
      CALL CHKXER('DORGLQ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DORGLQ(2,2,0,a,2,x,w,1,info)
      CALL CHKXER('DORGLQ',INFot,NOUt,LERr,OK)
!
!     DORGL2
!
      SRNamt = 'DORGL2'
      INFot = 1
      CALL DORGL2(-1,0,0,a,1,x,w,info)
      CALL CHKXER('DORGL2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORGL2(0,-1,0,a,1,x,w,info)
      CALL CHKXER('DORGL2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORGL2(2,1,0,a,2,x,w,info)
      CALL CHKXER('DORGL2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORGL2(0,0,-1,a,1,x,w,info)
      CALL CHKXER('DORGL2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORGL2(1,1,2,a,1,x,w,info)
      CALL CHKXER('DORGL2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORGL2(2,2,0,a,1,x,w,info)
      CALL CHKXER('DORGL2',INFot,NOUt,LERr,OK)
!
!     DORMLQ
!
      SRNamt = 'DORMLQ'
      INFot = 1
      CALL DORMLQ('/','N',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMLQ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORMLQ('L','/',0,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMLQ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORMLQ('L','N',-1,0,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMLQ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DORMLQ('L','N',0,-1,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMLQ('L','N',0,0,-1,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMLQ('L','N',0,1,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMLQ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORMLQ('R','N',1,0,1,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMLQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORMLQ('L','N',2,0,2,a,1,x,af,2,w,1,info)
      CALL CHKXER('DORMLQ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORMLQ('R','N',0,2,2,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMLQ',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DORMLQ('L','N',2,1,0,a,2,x,af,1,w,1,info)
      CALL CHKXER('DORMLQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL DORMLQ('L','N',1,2,0,a,1,x,af,1,w,1,info)
      CALL CHKXER('DORMLQ',INFot,NOUt,LERr,OK)
      INFot = 12
      CALL DORMLQ('R','N',2,1,0,a,1,x,af,2,w,1,info)
      CALL CHKXER('DORMLQ',INFot,NOUt,LERr,OK)
!
!     DORML2
!
      SRNamt = 'DORML2'
      INFot = 1
      CALL DORML2('/','N',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORML2',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DORML2('L','/',0,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORML2',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DORML2('L','N',-1,0,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORML2',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DORML2('L','N',0,-1,0,a,1,x,af,1,w,info)
      CALL CHKXER('DORML2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORML2('L','N',0,0,-1,a,1,x,af,1,w,info)
      CALL CHKXER('DORML2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORML2('L','N',0,1,1,a,1,x,af,1,w,info)
      CALL CHKXER('DORML2',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DORML2('R','N',1,0,1,a,1,x,af,1,w,info)
      CALL CHKXER('DORML2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORML2('L','N',2,1,2,a,1,x,af,2,w,info)
      CALL CHKXER('DORML2',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DORML2('R','N',1,2,2,a,1,x,af,1,w,info)
      CALL CHKXER('DORML2',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DORML2('L','N',2,1,0,a,2,x,af,1,w,info)
      CALL CHKXER('DORML2',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of DERRLQ
!
      END SUBROUTINE DERRLQ
