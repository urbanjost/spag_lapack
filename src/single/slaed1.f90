!*==slaed1.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLAED1 used by sstedc. Computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. Used when the original matrix is tridiagonal.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAED1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAED1( N, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK, IWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            CUTPNT, INFO, LDQ, N
!       REAL               RHO
!       ..
!       .. Array Arguments ..
!       INTEGER            INDXQ( * ), IWORK( * )
!       REAL               D( * ), Q( LDQ, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAED1 computes the updated eigensystem of a diagonal
!> matrix after modification by a rank-one symmetric matrix.  This
!> routine is used only for the eigenproblem which requires all
!> eigenvalues and eigenvectors of a tridiagonal matrix.  SLAED7 handles
!> the case in which eigenvalues only or eigenvalues and eigenvectors
!> of a full symmetric matrix (which was reduced to tridiagonal form)
!> are desired.
!>
!>   T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
!>
!>    where Z = Q**T*u, u is a vector of length N with ones in the
!>    CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
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
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>         On entry, the eigenvalues of the rank-1-perturbed matrix.
!>         On exit, the eigenvalues of the repaired matrix.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ,N)
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
!> \param[in,out] INDXQ
!> \verbatim
!>          INDXQ is INTEGER array, dimension (N)
!>         On entry, the permutation which separately sorts the two
!>         subproblems in D into ascending order.
!>         On exit, the permutation which will reintegrate the
!>         subproblems back into sorted order,
!>         i.e. D( INDXQ( I = 1, N ) ) will be in ascending order.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is REAL
!>         The subdiagonal entry used to create the rank-1 modification.
!> \endverbatim
!>
!> \param[in] CUTPNT
!> \verbatim
!>          CUTPNT is INTEGER
!>         The location of the last eigenvalue in the leading sub-matrix.
!>         min(1,N) <= CUTPNT <= N/2.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (4*N + N**2)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (4*N)
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
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA \n
!>  Modified by Francoise Tisseur, University of Tennessee
!>
!  =====================================================================
      SUBROUTINE SLAED1(N,D,Q,Ldq,Indxq,Rho,Cutpnt,Work,Iwork,Info)
      IMPLICIT NONE
!*--SLAED1166
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Cutpnt , Info , Ldq , N
      REAL Rho
!     ..
!     .. Array Arguments ..
      INTEGER Indxq(*) , Iwork(*)
      REAL D(*) , Q(Ldq,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER coltyp , cpp1 , i , idlmda , indx , indxc , indxp , iq2 , &
     &        is , iw , iz , k , n1 , n2
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SLAED2 , SLAED3 , SLAMRG , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Ldq<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( MIN(1,N/2)>Cutpnt .OR. (N/2)<Cutpnt ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLAED1',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     The following values are integer pointers which indicate
!     the portion of the workspace
!     used by a particular array in SLAED2 and SLAED3.
!
      iz = 1
      idlmda = iz + N
      iw = idlmda + N
      iq2 = iw + N
!
      indx = 1
      indxc = indx + N
      coltyp = indxc + N
      indxp = coltyp + N
!
!
!     Form the z-vector which consists of the last row of Q_1 and the
!     first row of Q_2.
!
      CALL SCOPY(Cutpnt,Q(Cutpnt,1),Ldq,Work(iz),1)
      cpp1 = Cutpnt + 1
      CALL SCOPY(N-Cutpnt,Q(cpp1,cpp1),Ldq,Work(iz+Cutpnt),1)
!
!     Deflate eigenvalues.
!
      CALL SLAED2(k,N,Cutpnt,D,Q,Ldq,Indxq,Rho,Work(iz),Work(idlmda),   &
     &            Work(iw),Work(iq2),Iwork(indx),Iwork(indxc),          &
     &            Iwork(indxp),Iwork(coltyp),Info)
!
      IF ( Info==0 ) THEN
!
!     Solve Secular Equation.
!
         IF ( k/=0 ) THEN
            is = (Iwork(coltyp)+Iwork(coltyp+1))                        &
     &           *Cutpnt + (Iwork(coltyp+1)+Iwork(coltyp+2))*(N-Cutpnt) &
     &           + iq2
            CALL SLAED3(k,N,Cutpnt,D,Q,Ldq,Rho,Work(idlmda),Work(iq2),  &
     &                  Iwork(indxc),Iwork(coltyp),Work(iw),Work(is),   &
     &                  Info)
            IF ( Info==0 ) THEN
!
!     Prepare the INDXQ sorting permutation.
!
               n1 = k
               n2 = N - k
               CALL SLAMRG(n1,n2,D,1,-1,Indxq)
            ENDIF
         ELSE
            DO i = 1 , N
               Indxq(i) = i
            ENDDO
         ENDIF
      ENDIF
!
!
!     End of SLAED1
!
      END SUBROUTINE SLAED1
