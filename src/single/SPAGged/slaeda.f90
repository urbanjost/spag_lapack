!*==slaeda.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLAEDA used by sstedc. Computes the Z vector determining the rank-one modification of the diagonal matrix. Used when the original matrix is dense.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAEDA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaeda.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaeda.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaeda.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR,
!                          GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            CURLVL, CURPBM, INFO, N, TLVLS
!       ..
!       .. Array Arguments ..
!       INTEGER            GIVCOL( 2, * ), GIVPTR( * ), PERM( * ),
!      $                   PRMPTR( * ), QPTR( * )
!       REAL               GIVNUM( 2, * ), Q( * ), Z( * ), ZTEMP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAEDA computes the Z vector corresponding to the merge step in the
!> CURLVLth step of the merge process with TLVLS steps for the CURPBMth
!> problem.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] TLVLS
!> \verbatim
!>          TLVLS is INTEGER
!>         The total number of merging levels in the overall divide and
!>         conquer tree.
!> \endverbatim
!>
!> \param[in] CURLVL
!> \verbatim
!>          CURLVL is INTEGER
!>         The current level in the overall merge routine,
!>         0 <= curlvl <= tlvls.
!> \endverbatim
!>
!> \param[in] CURPBM
!> \verbatim
!>          CURPBM is INTEGER
!>         The current problem in the current level in the overall
!>         merge routine (counting from upper left to lower right).
!> \endverbatim
!>
!> \param[in] PRMPTR
!> \verbatim
!>          PRMPTR is INTEGER array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in PERM a
!>         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
!>         indicates the size of the permutation and incidentally the
!>         size of the full, non-deflated problem.
!> \endverbatim
!>
!> \param[in] PERM
!> \verbatim
!>          PERM is INTEGER array, dimension (N lg N)
!>         Contains the permutations (from deflation and sorting) to be
!>         applied to each eigenblock.
!> \endverbatim
!>
!> \param[in] GIVPTR
!> \verbatim
!>          GIVPTR is INTEGER array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in GIVCOL a
!>         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
!>         indicates the number of Givens rotations.
!> \endverbatim
!>
!> \param[in] GIVCOL
!> \verbatim
!>          GIVCOL is INTEGER array, dimension (2, N lg N)
!>         Each pair of numbers indicates a pair of columns to take place
!>         in a Givens rotation.
!> \endverbatim
!>
!> \param[in] GIVNUM
!> \verbatim
!>          GIVNUM is REAL array, dimension (2, N lg N)
!>         Each number indicates the S value to be used in the
!>         corresponding Givens rotation.
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is REAL array, dimension (N**2)
!>         Contains the square eigenblocks from previous levels, the
!>         starting positions for blocks are given by QPTR.
!> \endverbatim
!>
!> \param[in] QPTR
!> \verbatim
!>          QPTR is INTEGER array, dimension (N+2)
!>         Contains a list of pointers which indicate where in Q an
!>         eigenblock is stored.  SQRT( QPTR(i+1) - QPTR(i) ) indicates
!>         the size of the block.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is REAL array, dimension (N)
!>         On output this vector contains the updating vector (the last
!>         row of the first sub-eigenvector matrix and the first row of
!>         the second sub-eigenvector matrix).
!> \endverbatim
!>
!> \param[out] ZTEMP
!> \verbatim
!>          ZTEMP is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
      SUBROUTINE SLAEDA(N,Tlvls,Curlvl,Curpbm,Prmptr,Perm,Givptr,Givcol,&
     &                  Givnum,Q,Qptr,Z,Ztemp,Info)
      USE S_SCOPY
      USE S_SGEMV
      USE S_SROT
      USE S_XERBLA
      IMPLICIT NONE
!*--SLAEDA174
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , HALF = 0.5E0 , ONE = 1.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Tlvls
      INTEGER , INTENT(IN) :: Curlvl
      INTEGER , INTENT(IN) :: Curpbm
      INTEGER , INTENT(IN) , DIMENSION(*) :: Prmptr
      INTEGER , INTENT(IN) , DIMENSION(*) :: Perm
      INTEGER , INTENT(IN) , DIMENSION(*) :: Givptr
      INTEGER , INTENT(IN) , DIMENSION(2,*) :: Givcol
      REAL , DIMENSION(2,*) :: Givnum
      REAL , DIMENSION(*) :: Q
      INTEGER , INTENT(IN) , DIMENSION(*) :: Qptr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL , DIMENSION(*) :: Ztemp
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: bsiz1 , bsiz2 , curr , i , k , mid , psiz1 , psiz2 ,   &
     &           ptr , zptr1
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
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( N<0 ) Info = -1
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLAEDA',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Determine location of first number in second half.
!
      mid = N/2 + 1
!
!     Gather last/first rows of appropriate eigenblocks into center of Z
!
      ptr = 1
!
!     Determine location of lowest level subproblem in the full storage
!     scheme
!
      curr = ptr + Curpbm*2**Curlvl + 2**(Curlvl-1) - 1
!
!     Determine size of these matrices.  We add HALF to the value of
!     the SQRT in case the machine underestimates one of these square
!     roots.
!
      bsiz1 = INT(HALF+SQRT(REAL(Qptr(curr+1)-Qptr(curr))))
      bsiz2 = INT(HALF+SQRT(REAL(Qptr(curr+2)-Qptr(curr+1))))
      DO k = 1 , mid - bsiz1 - 1
         Z(k) = ZERO
      ENDDO
      CALL SCOPY(bsiz1,Q(Qptr(curr)+bsiz1-1),bsiz1,Z(mid-bsiz1),1)
      CALL SCOPY(bsiz2,Q(Qptr(curr+1)),bsiz2,Z(mid),1)
      DO k = mid + bsiz2 , N
         Z(k) = ZERO
      ENDDO
!
!     Loop through remaining levels 1 -> CURLVL applying the Givens
!     rotations and permutation and then multiplying the center matrices
!     against the current Z.
!
      ptr = 2**Tlvls + 1
      DO k = 1 , Curlvl - 1
         curr = ptr + Curpbm*2**(Curlvl-k) + 2**(Curlvl-k-1) - 1
         psiz1 = Prmptr(curr+1) - Prmptr(curr)
         psiz2 = Prmptr(curr+2) - Prmptr(curr+1)
         zptr1 = mid - psiz1
!
!       Apply Givens at CURR and CURR+1
!
         DO i = Givptr(curr) , Givptr(curr+1) - 1
            CALL SROT(1,Z(zptr1+Givcol(1,i)-1),1,Z(zptr1+Givcol(2,i)-1),&
     &                1,Givnum(1,i),Givnum(2,i))
         ENDDO
         DO i = Givptr(curr+1) , Givptr(curr+2) - 1
            CALL SROT(1,Z(mid-1+Givcol(1,i)),1,Z(mid-1+Givcol(2,i)),1,  &
     &                Givnum(1,i),Givnum(2,i))
         ENDDO
         psiz1 = Prmptr(curr+1) - Prmptr(curr)
         psiz2 = Prmptr(curr+2) - Prmptr(curr+1)
         DO i = 0 , psiz1 - 1
            Ztemp(i+1) = Z(zptr1+Perm(Prmptr(curr)+i)-1)
         ENDDO
         DO i = 0 , psiz2 - 1
            Ztemp(psiz1+i+1) = Z(mid+Perm(Prmptr(curr+1)+i)-1)
         ENDDO
!
!        Multiply Blocks at CURR and CURR+1
!
!        Determine size of these matrices.  We add HALF to the value of
!        the SQRT in case the machine underestimates one of these
!        square roots.
!
         bsiz1 = INT(HALF+SQRT(REAL(Qptr(curr+1)-Qptr(curr))))
         bsiz2 = INT(HALF+SQRT(REAL(Qptr(curr+2)-Qptr(curr+1))))
         IF ( bsiz1>0 ) CALL SGEMV('T',bsiz1,bsiz1,ONE,Q(Qptr(curr)),   &
     &                             bsiz1,Ztemp(1),1,ZERO,Z(zptr1),1)
         CALL SCOPY(psiz1-bsiz1,Ztemp(bsiz1+1),1,Z(zptr1+bsiz1),1)
         IF ( bsiz2>0 ) CALL SGEMV('T',bsiz2,bsiz2,ONE,Q(Qptr(curr+1)), &
     &                             bsiz2,Ztemp(psiz1+1),1,ZERO,Z(mid),1)
         CALL SCOPY(psiz2-bsiz2,Ztemp(psiz1+bsiz2+1),1,Z(mid+bsiz2),1)
!
         ptr = ptr + 2**(Tlvls-k)
      ENDDO
!
!
!     End of SLAEDA
!
      END SUBROUTINE SLAEDA
