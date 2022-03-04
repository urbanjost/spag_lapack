!*==dlasdt.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLASDT creates a tree of subproblems for bidiagonal divide and conquer. Used by sbdsdc.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASDT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasdt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasdt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasdt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB )
!
!       .. Scalar Arguments ..
!       INTEGER            LVL, MSUB, N, ND
!       ..
!       .. Array Arguments ..
!       INTEGER            INODE( * ), NDIML( * ), NDIMR( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASDT creates a tree of subproblems for bidiagonal divide and
!> conquer.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          On entry, the number of diagonal elements of the
!>          bidiagonal matrix.
!> \endverbatim
!>
!> \param[out] LVL
!> \verbatim
!>          LVL is INTEGER
!>          On exit, the number of levels on the computation tree.
!> \endverbatim
!>
!> \param[out] ND
!> \verbatim
!>          ND is INTEGER
!>          On exit, the number of nodes on the tree.
!> \endverbatim
!>
!> \param[out] INODE
!> \verbatim
!>          INODE is INTEGER array, dimension ( N )
!>          On exit, centers of subproblems.
!> \endverbatim
!>
!> \param[out] NDIML
!> \verbatim
!>          NDIML is INTEGER array, dimension ( N )
!>          On exit, row dimensions of left children.
!> \endverbatim
!>
!> \param[out] NDIMR
!> \verbatim
!>          NDIMR is INTEGER array, dimension ( N )
!>          On exit, row dimensions of right children.
!> \endverbatim
!>
!> \param[in] MSUB
!> \verbatim
!>          MSUB is INTEGER
!>          On entry, the maximum row dimension each subproblem at the
!>          bottom of the tree can be of.
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
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Huan Ren, Computer Science Division, University of
!>     California at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE DLASDT(N,Lvl,Nd,Inode,Ndiml,Ndimr,Msub)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DLASDT110
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  TWO = 2.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(INOUT) :: Lvl
      INTEGER , INTENT(OUT) :: Nd
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Inode
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ndiml
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ndimr
      INTEGER , INTENT(IN) :: Msub
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , il , ir , llst , maxn , ncrnt , nlvl
      REAL(R8KIND) :: temp
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Find the number of levels on the tree.
!
      maxn = MAX(1,N)
      temp = LOG(DBLE(maxn)/DBLE(Msub+1))/LOG(TWO)
      Lvl = INT(temp) + 1
!
      i = N/2
      Inode(1) = i + 1
      Ndiml(1) = i
      Ndimr(1) = N - i - 1
      il = 0
      ir = 1
      llst = 1
      DO nlvl = 1 , Lvl - 1
!
!        Constructing the tree at (NLVL+1)-st level. The number of
!        nodes created on this level is LLST * 2.
!
         DO i = 0 , llst - 1
            il = il + 2
            ir = ir + 2
            ncrnt = llst + i
            Ndiml(il) = Ndiml(ncrnt)/2
            Ndimr(il) = Ndiml(ncrnt) - Ndiml(il) - 1
            Inode(il) = Inode(ncrnt) - Ndimr(il) - 1
            Ndiml(ir) = Ndimr(ncrnt)/2
            Ndimr(ir) = Ndimr(ncrnt) - Ndiml(ir) - 1
            Inode(ir) = Inode(ncrnt) + Ndiml(ir) + 1
         ENDDO
         llst = llst*2
      ENDDO
      Nd = llst*2 - 1
!
!
!     End of DLASDT
!
      END SUBROUTINE DLASDT
