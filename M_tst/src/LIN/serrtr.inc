!*==serrtr.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b serrtr
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE serrtr( PATH, NUNIT )
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
!> serrtr tests the error exits for the REAL triangular
!> routines.
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE serrtr(Path,Nunit)
      IMPLICIT NONE
!*--SERRTR59
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
      CHARACTER*2 c2
      INTEGER info
      REAL rcond , scale
!     ..
!     .. Local Arrays ..
      INTEGER iw(NMAX)
      REAL a(NMAX,NMAX) , b(NMAX) , r1(NMAX) , r2(NMAX) , w(NMAX) ,     &
     &     x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    SLATBS,   SLATPS,   SLATRS,   STBCON
EXTERNAL    STBRFS,   STBTRS,   STPCON,   STPRFS,   STPTRI,   STPTRS
EXTERNAL    STRCON,   STRRFS,   STRTI2,   STRTRI,   STRTRS
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
!     .. Executable Statements ..
!
      NOUt = Nunit
      WRITE (NOUt,FMT=*)
      c2 = Path(2:3)
      a(1,1) = 1.
      a(1,2) = 2.
      a(2,2) = 3.
      a(2,1) = 4.
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'TR') ) THEN
!
!        Test error exits for the general triangular routines.
!
!        STRTRI
!
         SRNamt = 'STRTRI'
         INFot = 1
         CALL STRTRI('/','N',0,a,1,info)
         CALL CHKXER('STRTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STRTRI('U','/',0,a,1,info)
         CALL CHKXER('STRTRI',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STRTRI('U','N',-1,a,1,info)
         CALL CHKXER('STRTRI',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL STRTRI('U','N',2,a,1,info)
         CALL CHKXER('STRTRI',INFot,NOUt,LERr,OK)
!
!        STRTI2
!
         SRNamt = 'STRTI2'
         INFot = 1
         CALL STRTI2('/','N',0,a,1,info)
         CALL CHKXER('STRTI2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STRTI2('U','/',0,a,1,info)
         CALL CHKXER('STRTI2',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STRTI2('U','N',-1,a,1,info)
         CALL CHKXER('STRTI2',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL STRTI2('U','N',2,a,1,info)
         CALL CHKXER('STRTI2',INFot,NOUt,LERr,OK)
!
!        STRTRS
!
         SRNamt = 'STRTRS'
         INFot = 1
         CALL STRTRS('/','N','N',0,0,a,1,x,1,info)
         CALL CHKXER('STRTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STRTRS('U','/','N',0,0,a,1,x,1,info)
         CALL CHKXER('STRTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STRTRS('U','N','/',0,0,a,1,x,1,info)
         CALL CHKXER('STRTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STRTRS('U','N','N',-1,0,a,1,x,1,info)
         CALL CHKXER('STRTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL STRTRS('U','N','N',0,-1,a,1,x,1,info)
         CALL CHKXER('STRTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL STRTRS('U','N','N',2,1,a,1,x,2,info)
         CALL CHKXER('STRTRS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL STRTRS('U','N','N',2,1,a,2,x,1,info)
         CALL CHKXER('STRTRS',INFot,NOUt,LERr,OK)
!
!        STRRFS
!
         SRNamt = 'STRRFS'
         INFot = 1
         CALL STRRFS('/','N','N',0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STRRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STRRFS('U','/','N',0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STRRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STRRFS('U','N','/',0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STRRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STRRFS('U','N','N',-1,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STRRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL STRRFS('U','N','N',0,-1,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STRRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL STRRFS('U','N','N',2,1,a,1,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('STRRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL STRRFS('U','N','N',2,1,a,2,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('STRRFS',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL STRRFS('U','N','N',2,1,a,2,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STRRFS',INFot,NOUt,LERr,OK)
!
!        STRCON
!
         SRNamt = 'STRCON'
         INFot = 1
         CALL STRCON('/','U','N',0,a,1,rcond,w,iw,info)
         CALL CHKXER('STRCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STRCON('1','/','N',0,a,1,rcond,w,iw,info)
         CALL CHKXER('STRCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STRCON('1','U','/',0,a,1,rcond,w,iw,info)
         CALL CHKXER('STRCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STRCON('1','U','N',-1,a,1,rcond,w,iw,info)
         CALL CHKXER('STRCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL STRCON('1','U','N',2,a,1,rcond,w,iw,info)
         CALL CHKXER('STRCON',INFot,NOUt,LERr,OK)
!
!        SLATRS
!
         SRNamt = 'SLATRS'
         INFot = 1
         CALL SLATRS('/','N','N','N',0,a,1,x,scale,w,info)
         CALL CHKXER('SLATRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SLATRS('U','/','N','N',0,a,1,x,scale,w,info)
         CALL CHKXER('SLATRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SLATRS('U','N','/','N',0,a,1,x,scale,w,info)
         CALL CHKXER('SLATRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SLATRS('U','N','N','/',0,a,1,x,scale,w,info)
         CALL CHKXER('SLATRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SLATRS('U','N','N','N',-1,a,1,x,scale,w,info)
         CALL CHKXER('SLATRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SLATRS('U','N','N','N',2,a,1,x,scale,w,info)
         CALL CHKXER('SLATRS',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'TP') ) THEN
!
!        Test error exits for the packed triangular routines.
!
!        STPTRI
!
         SRNamt = 'STPTRI'
         INFot = 1
         CALL STPTRI('/','N',0,a,info)
         CALL CHKXER('STPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STPTRI('U','/',0,a,info)
         CALL CHKXER('STPTRI',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STPTRI('U','N',-1,a,info)
         CALL CHKXER('STPTRI',INFot,NOUt,LERr,OK)
!
!        STPTRS
!
         SRNamt = 'STPTRS'
         INFot = 1
         CALL STPTRS('/','N','N',0,0,a,x,1,info)
         CALL CHKXER('STPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STPTRS('U','/','N',0,0,a,x,1,info)
         CALL CHKXER('STPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STPTRS('U','N','/',0,0,a,x,1,info)
         CALL CHKXER('STPTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STPTRS('U','N','N',-1,0,a,x,1,info)
         CALL CHKXER('STPTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL STPTRS('U','N','N',0,-1,a,x,1,info)
         CALL CHKXER('STPTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL STPTRS('U','N','N',2,1,a,x,1,info)
         CALL CHKXER('STPTRS',INFot,NOUt,LERr,OK)
!
!        STPRFS
!
         SRNamt = 'STPRFS'
         INFot = 1
         CALL STPRFS('/','N','N',0,0,a,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STPRFS('U','/','N',0,0,a,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STPRFS('U','N','/',0,0,a,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STPRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STPRFS('U','N','N',-1,0,a,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STPRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL STPRFS('U','N','N',0,-1,a,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STPRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL STPRFS('U','N','N',2,1,a,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('STPRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL STPRFS('U','N','N',2,1,a,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STPRFS',INFot,NOUt,LERr,OK)
!
!        STPCON
!
         SRNamt = 'STPCON'
         INFot = 1
         CALL STPCON('/','U','N',0,a,rcond,w,iw,info)
         CALL CHKXER('STPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STPCON('1','/','N',0,a,rcond,w,iw,info)
         CALL CHKXER('STPCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STPCON('1','U','/',0,a,rcond,w,iw,info)
         CALL CHKXER('STPCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STPCON('1','U','N',-1,a,rcond,w,iw,info)
         CALL CHKXER('STPCON',INFot,NOUt,LERr,OK)
!
!        SLATPS
!
         SRNamt = 'SLATPS'
         INFot = 1
         CALL SLATPS('/','N','N','N',0,a,x,scale,w,info)
         CALL CHKXER('SLATPS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SLATPS('U','/','N','N',0,a,x,scale,w,info)
         CALL CHKXER('SLATPS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SLATPS('U','N','/','N',0,a,x,scale,w,info)
         CALL CHKXER('SLATPS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SLATPS('U','N','N','/',0,a,x,scale,w,info)
         CALL CHKXER('SLATPS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SLATPS('U','N','N','N',-1,a,x,scale,w,info)
         CALL CHKXER('SLATPS',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'TB') ) THEN
!
!        Test error exits for the banded triangular routines.
!
!        STBTRS
!
         SRNamt = 'STBTRS'
         INFot = 1
         CALL STBTRS('/','N','N',0,0,0,a,1,x,1,info)
         CALL CHKXER('STBTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STBTRS('U','/','N',0,0,0,a,1,x,1,info)
         CALL CHKXER('STBTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STBTRS('U','N','/',0,0,0,a,1,x,1,info)
         CALL CHKXER('STBTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STBTRS('U','N','N',-1,0,0,a,1,x,1,info)
         CALL CHKXER('STBTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL STBTRS('U','N','N',0,-1,0,a,1,x,1,info)
         CALL CHKXER('STBTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL STBTRS('U','N','N',0,0,-1,a,1,x,1,info)
         CALL CHKXER('STBTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL STBTRS('U','N','N',2,1,1,a,1,x,2,info)
         CALL CHKXER('STBTRS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL STBTRS('U','N','N',2,0,1,a,1,x,1,info)
         CALL CHKXER('STBTRS',INFot,NOUt,LERr,OK)
!
!        STBRFS
!
         SRNamt = 'STBRFS'
         INFot = 1
         CALL STBRFS('/','N','N',0,0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STBRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STBRFS('U','/','N',0,0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STBRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STBRFS('U','N','/',0,0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STBRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STBRFS('U','N','N',-1,0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STBRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL STBRFS('U','N','N',0,-1,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STBRFS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL STBRFS('U','N','N',0,0,-1,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STBRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL STBRFS('U','N','N',2,1,1,a,1,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('STBRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL STBRFS('U','N','N',2,1,1,a,2,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('STBRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL STBRFS('U','N','N',2,1,1,a,2,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('STBRFS',INFot,NOUt,LERr,OK)
!
!        STBCON
!
         SRNamt = 'STBCON'
         INFot = 1
         CALL STBCON('/','U','N',0,0,a,1,rcond,w,iw,info)
         CALL CHKXER('STBCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STBCON('1','/','N',0,0,a,1,rcond,w,iw,info)
         CALL CHKXER('STBCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STBCON('1','U','/',0,0,a,1,rcond,w,iw,info)
         CALL CHKXER('STBCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STBCON('1','U','N',-1,0,a,1,rcond,w,iw,info)
         CALL CHKXER('STBCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL STBCON('1','U','N',0,-1,a,1,rcond,w,iw,info)
         CALL CHKXER('STBCON',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL STBCON('1','U','N',2,1,a,1,rcond,w,iw,info)
         CALL CHKXER('STBCON',INFot,NOUt,LERr,OK)
!
!        SLATBS
!
         SRNamt = 'SLATBS'
         INFot = 1
         CALL SLATBS('/','N','N','N',0,0,a,1,x,scale,w,info)
         CALL CHKXER('SLATBS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SLATBS('U','/','N','N',0,0,a,1,x,scale,w,info)
         CALL CHKXER('SLATBS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SLATBS('U','N','/','N',0,0,a,1,x,scale,w,info)
         CALL CHKXER('SLATBS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SLATBS('U','N','N','/',0,0,a,1,x,scale,w,info)
         CALL CHKXER('SLATBS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SLATBS('U','N','N','N',-1,0,a,1,x,scale,w,info)
         CALL CHKXER('SLATBS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SLATBS('U','N','N','N',1,-1,a,1,x,scale,w,info)
         CALL CHKXER('SLATBS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SLATBS('U','N','N','N',2,1,a,1,x,scale,w,info)
         CALL CHKXER('SLATBS',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of serrtr
!
      END SUBROUTINE serrtr
