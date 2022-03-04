!*==claed7.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLAED7 used by sstedc. Computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. Used when the original matrix is dense.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAED7 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claed7.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claed7.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claed7.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAED7( N, CUTPNT, QSIZ, TLVLS, CURLVL, CURPBM, D, Q,
!                          LDQ, RHO, INDXQ, QSTORE, QPTR, PRMPTR, PERM,
!                          GIVPTR, GIVCOL, GIVNUM, WORK, RWORK, IWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            CURLVL, CURPBM, CUTPNT, INFO, LDQ, N, QSIZ,
!      $                   TLVLS
!       REAL               RHO
!       ..
!       .. Array Arguments ..
!       INTEGER            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ),
!      $                   IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )
!       REAL               D( * ), GIVNUM( 2, * ), QSTORE( * ), RWORK( * )
!       COMPLEX            Q( LDQ, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAED7 computes the updated eigensystem of a diagonal
!> matrix after modification by a rank-one symmetric matrix. This
!> routine is used only for the eigenproblem which requires all
!> eigenvalues and optionally eigenvectors of a dense or banded
!> Hermitian matrix that has been reduced to tridiagonal form.
!>
!>   T = Q(in) ( D(in) + RHO * Z*Z**H ) Q**H(in) = Q(out) * D(out) * Q**H(out)
!>
!>   where Z = Q**Hu, u is a vector of length N with ones in the
!>   CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
!>
!>    The eigenvectors of the original matrix are stored in Q, and the
!>    eigenvalues are in D.  The algorithm consists of three stages:
!>
!>       The first stage consists of deflating the size of the problem
!>       when there are multiple eigenvalues or if there is a zero in
!>       the Z vector.  For each such occurrence the dimension of the
!>       secular equation problem is reduced by one.  This stage is
!>       performed by the routine SLAED2.
!>
!>       The second stage consists of calculating the updated
!>       eigenvalues. This is done by finding the roots of the secular
!>       equation via the routine SLAED4 (as called by SLAED3).
!>       This routine also calculates the eigenvectors of the current
!>       problem.
!>
!>       The final stage consists of computing the updated eigenvectors
!>       directly using the updated eigenvalues.  The eigenvectors for
!>       the current problem are multiplied with the eigenvectors from
!>       the overall problem.
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
!> \param[in] CUTPNT
!> \verbatim
!>          CUTPNT is INTEGER
!>         Contains the location of the last eigenvalue in the leading
!>         sub-matrix.  min(1,N) <= CUTPNT <= N.
!> \endverbatim
!>
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is INTEGER
!>         The dimension of the unitary matrix used to reduce
!>         the full matrix to tridiagonal form.  QSIZ >= N.
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
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>         On entry, the eigenvalues of the rank-1-perturbed matrix.
!>         On exit, the eigenvalues of the repaired matrix.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDQ,N)
!>         On entry, the eigenvectors of the rank-1-perturbed matrix.
!>         On exit, the eigenvectors of the repaired tridiagonal matrix.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>         The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is REAL
!>         Contains the subdiagonal element used to create the rank-1
!>         modification.
!> \endverbatim
!>
!> \param[out] INDXQ
!> \verbatim
!>          INDXQ is INTEGER array, dimension (N)
!>         This contains the permutation which will reintegrate the
!>         subproblem just solved back into sorted order,
!>         ie. D( INDXQ( I = 1, N ) ) will be in ascending order.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (4*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array,
!>                                 dimension (3*N+2*QSIZ*N)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (QSIZ*N)
!> \endverbatim
!>
!> \param[in,out] QSTORE
!> \verbatim
!>          QSTORE is REAL array, dimension (N**2+1)
!>         Stores eigenvectors of submatrices encountered during
!>         divide and conquer, packed together. QPTR points to
!>         beginning of the submatrices.
!> \endverbatim
!>
!> \param[in,out] QPTR
!> \verbatim
!>          QPTR is INTEGER array, dimension (N+2)
!>         List of indices pointing to beginning of submatrices stored
!>         in QSTORE. The submatrices are numbered starting at the
!>         bottom left of the divide and conquer tree, from left to
!>         right and bottom to top.
!> \endverbatim
!>
!> \param[in] PRMPTR
!> \verbatim
!>          PRMPTR is INTEGER array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in PERM a
!>         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
!>         indicates the size of the permutation and also the size of
!>         the full, non-deflated problem.
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
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, an eigenvalue did not converge
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
!> \date June 2016
!
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CLAED7(N,Cutpnt,Qsiz,Tlvls,Curlvl,Curpbm,D,Q,Ldq,Rho,  &
     &                  Indxq,Qstore,Qptr,Prmptr,Perm,Givptr,Givcol,    &
     &                  Givnum,Work,Rwork,Iwork,Info)
      USE S_CLACRM
      USE S_CLAED8
      USE S_SLAED9
      USE S_SLAEDA
      USE S_SLAMRG
      USE S_XERBLA
      IMPLICIT NONE
!*--CLAED7258
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      INTEGER :: Cutpnt
      INTEGER :: Qsiz
      INTEGER :: Tlvls
      INTEGER :: Curlvl
      INTEGER :: Curpbm
      REAL , DIMENSION(*) :: D
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL :: Rho
      INTEGER , DIMENSION(*) :: Indxq
      REAL , DIMENSION(*) :: Qstore
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Qptr
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Prmptr
      INTEGER , DIMENSION(*) :: Perm
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Givptr
      INTEGER , DIMENSION(2,*) :: Givcol
      REAL , DIMENSION(2,*) :: Givnum
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: coltyp , curr , i , idlmda , indx , indxc , indxp ,    &
     &           iq , iw , iz , k , n1 , n2 , ptr
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
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
!     IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
!        INFO = -1
!     ELSE IF( N.LT.0 ) THEN
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( MIN(1,N)>Cutpnt .OR. N<Cutpnt ) THEN
         Info = -2
      ELSEIF ( Qsiz<N ) THEN
         Info = -3
      ELSEIF ( Ldq<MAX(1,N) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CLAED7',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     The following values are for bookkeeping purposes only.  They are
!     integer pointers which indicate the portion of the workspace
!     used by a particular array in SLAED2 and SLAED3.
!
      iz = 1
      idlmda = iz + N
      iw = idlmda + N
      iq = iw + N
!
      indx = 1
      indxc = indx + N
      coltyp = indxc + N
      indxp = coltyp + N
!
!     Form the z-vector which consists of the last row of Q_1 and the
!     first row of Q_2.
!
      ptr = 1 + 2**Tlvls
      DO i = 1 , Curlvl - 1
         ptr = ptr + 2**(Tlvls-i)
      ENDDO
      curr = ptr + Curpbm
      CALL SLAEDA(N,Tlvls,Curlvl,Curpbm,Prmptr,Perm,Givptr,Givcol,      &
     &            Givnum,Qstore,Qptr,Rwork(iz),Rwork(iz+N),Info)
!
!     When solving the final problem, we no longer need the stored data,
!     so we will overwrite the data from this level onto the previously
!     used storage space.
!
      IF ( Curlvl==Tlvls ) THEN
         Qptr(curr) = 1
         Prmptr(curr) = 1
         Givptr(curr) = 1
      ENDIF
!
!     Sort and Deflate eigenvalues.
!
      CALL CLAED8(k,N,Qsiz,Q,Ldq,D,Rho,Cutpnt,Rwork(iz),Rwork(idlmda),  &
     &            Work,Qsiz,Rwork(iw),Iwork(indxp),Iwork(indx),Indxq,   &
     &            Perm(Prmptr(curr)),Givptr(curr+1),                    &
     &            Givcol(1,Givptr(curr)),Givnum(1,Givptr(curr)),Info)
      Prmptr(curr+1) = Prmptr(curr) + N
      Givptr(curr+1) = Givptr(curr+1) + Givptr(curr)
!
!     Solve Secular Equation.
!
      IF ( k/=0 ) THEN
         CALL SLAED9(k,1,k,N,D,Rwork(iq),k,Rho,Rwork(idlmda),Rwork(iw), &
     &               Qstore(Qptr(curr)),k,Info)
         CALL CLACRM(Qsiz,k,Work,Qsiz,Qstore(Qptr(curr)),k,Q,Ldq,       &
     &               Rwork(iq))
         Qptr(curr+1) = Qptr(curr) + k**2
         IF ( Info/=0 ) RETURN
!
!     Prepare the INDXQ sorting premutation.
!
         n1 = k
         n2 = N - k
         CALL SLAMRG(n1,n2,D,1,-1,Indxq)
      ELSE
         Qptr(curr+1) = Qptr(curr)
         DO i = 1 , N
            Indxq(i) = i
         ENDDO
      ENDIF
!
!
!     End of CLAED7
!
      END SUBROUTINE CLAED7
