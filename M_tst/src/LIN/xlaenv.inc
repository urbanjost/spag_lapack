!*==xlaenv.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b xlaenv
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE xlaenv( ISPEC, NVALUE )
!
!       .. Scalar Arguments ..
!       INTEGER            ISPEC, NVALUE
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> xlaenv sets certain machine- and problem-dependent quantities
!> which will later be retrieved by ilaenv.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies the parameter to be set in the COMMON array IPARMS.
!>          = 1: the optimal blocksize; if this value is 1, an unblocked
!>               algorithm will give the best performance.
!>          = 2: the minimum block size for which the block routine
!>               should be used; if the usable block size is less than
!>               this value, an unblocked routine should be used.
!>          = 3: the crossover point (in a block routine, for N less
!>               than this value, an unblocked routine should be used)
!>          = 4: the number of shifts, used in the nonsymmetric
!>               eigenvalue routines
!>          = 5: the minimum column dimension for blocking to be used;
!>               rectangular blocks must have dimension at least k by m,
!>               where k is given by ilaenv(2,...) and m by ilaenv(5,...)
!>          = 6: the crossover point for the SVD (when reducing an m by n
!>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!>               this value, a QR factorization is used first to reduce
!>               the matrix to a triangular form)
!>          = 7: the number of processors
!>          = 8: another crossover point, for the multishift QR and QZ
!>               methods for nonsymmetric eigenvalue problems.
!>          = 9: maximum size of the subproblems at the bottom of the
!>               computation tree in the divide-and-conquer algorithm
!>               (used by xGELSD and xGESDD)
!>          =10: ieee NaN arithmetic can be trusted not to trap
!>          =11: infinity arithmetic can be trusted not to trap
!> \endverbatim
!>
!> \param[in] NVALUE
!> \verbatim
!>          NVALUE is INTEGER
!>          The value of the parameter specified by ISPEC.
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
!> \ingroup aux_lin
!
!  =====================================================================
      SUBROUTINE xlaenv(Ispec,Nvalue)
      IMPLICIT NONE
!*--XLAENV82
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ispec , Nvalue
!     ..
!
!  =====================================================================
!
!     .. Arrays in Common ..
      INTEGER IPArms(100)
!     ..
!     .. Common blocks ..
      COMMON /CLAENV/ IPArms
!     ..
!     .. Save statement ..
      SAVE /CLAENV/ 
!     ..
!     .. Executable Statements ..
!
      IF ( Ispec>=1 .AND. Ispec<=9 ) IPArms(Ispec) = Nvalue
!
!
!     End of xlaenv
!
      END SUBROUTINE xlaenv
