!*==slasdt.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLASDT creates a tree of subproblems for bidiagonal divide and conquer. Used by sbdsdc.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASDT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasdt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasdt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasdt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB )
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
!> SLASDT creates a tree of subproblems for bidiagonal divide and
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
      SUBROUTINE SLASDT(N,Lvl,Nd,Inode,Ndiml,Ndimr,Msub)
      IMPLICIT NONE
!*--SLASDT109
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lvl , Msub , N , Nd
!     ..
!     .. Array Arguments ..
      INTEGER Inode(*) , Ndiml(*) , Ndimr(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL TWO
      PARAMETER (TWO=2.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , il , ir , llst , maxn , ncrnt , nlvl
      REAL temp
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC INT , LOG , MAX , REAL
!     ..
!     .. Executable Statements ..
!
!     Find the number of levels on the tree.
!
      maxn = MAX(1,N)
      temp = LOG(REAL(maxn)/REAL(Msub+1))/LOG(TWO)
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
!     End of SLASDT
!
      END SUBROUTINE SLASDT
