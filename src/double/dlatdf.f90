!*==dlatdf.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLATDF uses the LU factorization of the n-by-n matrix computed by sgetc2 and computes a contribution to the reciprocal Dif-estimate.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLATDF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlatdf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlatdf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlatdf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV,
!                          JPIV )
!
!       .. Scalar Arguments ..
!       INTEGER            IJOB, LDZ, N
!       DOUBLE PRECISION   RDSCAL, RDSUM
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), JPIV( * )
!       DOUBLE PRECISION   RHS( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLATDF uses the LU factorization of the n-by-n matrix Z computed by
!> DGETC2 and computes a contribution to the reciprocal Dif-estimate
!> by solving Z * x = b for x, and choosing the r.h.s. b such that
!> the norm of x is as large as possible. On entry RHS = b holds the
!> contribution from earlier solved sub-systems, and on return RHS = x.
!>
!> The factorization of Z returned by DGETC2 has the form Z = P*L*U*Q,
!> where P and Q are permutation matrices. L is lower triangular with
!> unit diagonal elements and U is upper triangular.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IJOB
!> \verbatim
!>          IJOB is INTEGER
!>          IJOB = 2: First compute an approximative null-vector e
!>              of Z using DGECON, e is normalized and solve for
!>              Zx = +-e - f with the sign giving the greater value
!>              of 2-norm(x). About 5 times as expensive as Default.
!>          IJOB .ne. 2: Local look ahead strategy where all entries of
!>              the r.h.s. b is chosen as either +1 or -1 (Default).
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Z.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
!>          On entry, the LU part of the factorization of the n-by-n
!>          matrix Z computed by DGETC2:  Z = P * L * U * Q
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDA >= max(1, N).
!> \endverbatim
!>
!> \param[in,out] RHS
!> \verbatim
!>          RHS is DOUBLE PRECISION array, dimension (N)
!>          On entry, RHS contains contributions from other subsystems.
!>          On exit, RHS contains the solution of the subsystem with
!>          entries according to the value of IJOB (see above).
!> \endverbatim
!>
!> \param[in,out] RDSUM
!> \verbatim
!>          RDSUM is DOUBLE PRECISION
!>          On entry, the sum of squares of computed contributions to
!>          the Dif-estimate under computation by DTGSYL, where the
!>          scaling factor RDSCAL (see below) has been factored out.
!>          On exit, the corresponding sum of squares updated with the
!>          contributions from the current sub-system.
!>          If TRANS = 'T' RDSUM is not touched.
!>          NOTE: RDSUM only makes sense when DTGSY2 is called by STGSYL.
!> \endverbatim
!>
!> \param[in,out] RDSCAL
!> \verbatim
!>          RDSCAL is DOUBLE PRECISION
!>          On entry, scaling factor used to prevent overflow in RDSUM.
!>          On exit, RDSCAL is updated w.r.t. the current contributions
!>          in RDSUM.
!>          If TRANS = 'T', RDSCAL is not touched.
!>          NOTE: RDSCAL only makes sense when DTGSY2 is called by
!>                DTGSYL.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N).
!>          The pivot indices; for 1 <= i <= N, row i of the
!>          matrix has been interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in] JPIV
!> \verbatim
!>          JPIV is INTEGER array, dimension (N).
!>          The pivot indices; for 1 <= j <= N, column j of the
!>          matrix has been interchanged with column JPIV(j).
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
!> \ingroup doubleOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!>  This routine is a further developed implementation of algorithm
!>  BSOLVE in [1] using complete pivoting in the LU factorization.
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!> \par References:
!  ================
!>
!> \verbatim
!>
!>
!>  [1] Bo Kagstrom and Lars Westin,
!>      Generalized Schur Methods with Condition Estimators for
!>      Solving the Generalized Sylvester Equation, IEEE Transactions
!>      on Automatic Control, Vol. 34, No. 7, July 1989, pp 745-751.
!>
!>  [2] Peter Poromaa,
!>      On Efficient and Robust Estimators for the Separation
!>      between two Regular Matrix Pairs with Applications in
!>      Condition Estimation. Report IMINF-95.05, Departement of
!>      Computing Science, Umea University, S-901 87 Umea, Sweden, 1995.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLATDF(Ijob,N,Z,Ldz,Rhs,Rdsum,Rdscal,Ipiv,Jpiv)
      IMPLICIT NONE
!*--DLATDF174
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Ijob , Ldz , N
      DOUBLE PRECISION Rdscal , Rdsum
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*) , Jpiv(*)
      DOUBLE PRECISION Rhs(*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER MAXDIM
      PARAMETER (MAXDIM=8)
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j , k
      DOUBLE PRECISION bm , bp , pmone , sminu , splus , temp
!     ..
!     .. Local Arrays ..
      INTEGER iwork(MAXDIM)
      DOUBLE PRECISION work(4*MAXDIM) , xm(MAXDIM) , xp(MAXDIM)
!     ..
!     .. External Subroutines ..
      EXTERNAL DAXPY , DCOPY , DGECON , DGESC2 , DLASSQ , DLASWP , DSCAL
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DASUM , DDOT
      EXTERNAL DASUM , DDOT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SQRT
!     ..
!     .. Executable Statements ..
!
      IF ( Ijob/=2 ) THEN
!
!        Apply permutations IPIV to RHS
!
         CALL DLASWP(1,Rhs,Ldz,1,N-1,Ipiv,1)
!
!        Solve for L-part choosing RHS either to +1 or -1.
!
         pmone = -ONE
!
         DO j = 1 , N - 1
            bp = Rhs(j) + ONE
            bm = Rhs(j) - ONE
            splus = ONE
!
!           Look-ahead for L-part RHS(1:N-1) = + or -1, SPLUS and
!           SMIN computed more efficiently than in BSOLVE [1].
!
            splus = splus + DDOT(N-j,Z(j+1,j),1,Z(j+1,j),1)
            sminu = DDOT(N-j,Z(j+1,j),1,Rhs(j+1),1)
            splus = splus*Rhs(j)
            IF ( splus>sminu ) THEN
               Rhs(j) = bp
            ELSEIF ( sminu>splus ) THEN
               Rhs(j) = bm
            ELSE
!
!              In this case the updating sums are equal and we can
!              choose RHS(J) +1 or -1. The first time this happens
!              we choose -1, thereafter +1. This is a simple way to
!              get good estimates of matrices like Byers well-known
!              example (see [1]). (Not done in BSOLVE.)
!
               Rhs(j) = Rhs(j) + pmone
               pmone = ONE
            ENDIF
!
!           Compute the remaining r.h.s.
!
            temp = -Rhs(j)
            CALL DAXPY(N-j,temp,Z(j+1,j),1,Rhs(j+1),1)
!
         ENDDO
!
!        Solve for U-part, look-ahead for RHS(N) = +-1. This is not done
!        in BSOLVE and will hopefully give us a better estimate because
!        any ill-conditioning of the original matrix is transferred to U
!        and not to L. U(N, N) is an approximation to sigma_min(LU).
!
         CALL DCOPY(N-1,Rhs,1,xp,1)
         xp(N) = Rhs(N) + ONE
         Rhs(N) = Rhs(N) - ONE
         splus = ZERO
         sminu = ZERO
         DO i = N , 1 , -1
            temp = ONE/Z(i,i)
            xp(i) = xp(i)*temp
            Rhs(i) = Rhs(i)*temp
            DO k = i + 1 , N
               xp(i) = xp(i) - xp(k)*(Z(i,k)*temp)
               Rhs(i) = Rhs(i) - Rhs(k)*(Z(i,k)*temp)
            ENDDO
            splus = splus + ABS(xp(i))
            sminu = sminu + ABS(Rhs(i))
         ENDDO
         IF ( splus>sminu ) CALL DCOPY(N,xp,1,Rhs,1)
!
!        Apply the permutations JPIV to the computed solution (RHS)
!
         CALL DLASWP(1,Rhs,Ldz,1,N-1,Jpiv,-1)
!
!        Compute the sum of squares
!
         CALL DLASSQ(N,Rhs,1,Rdscal,Rdsum)
!
      ELSE
!
!        IJOB = 2, Compute approximate nullvector XM of Z
!
         CALL DGECON('I',N,Z,Ldz,ONE,temp,work,iwork,info)
         CALL DCOPY(N,work(N+1),1,xm,1)
!
!        Compute RHS
!
         CALL DLASWP(1,xm,Ldz,1,N-1,Ipiv,-1)
         temp = ONE/SQRT(DDOT(N,xm,1,xm,1))
         CALL DSCAL(N,temp,xm,1)
         CALL DCOPY(N,xm,1,xp,1)
         CALL DAXPY(N,ONE,Rhs,1,xp,1)
         CALL DAXPY(N,-ONE,xm,1,Rhs,1)
         CALL DGESC2(N,Z,Ldz,Rhs,Ipiv,Jpiv,temp)
         CALL DGESC2(N,Z,Ldz,xp,Ipiv,Jpiv,temp)
         IF ( DASUM(N,xp,1)>DASUM(N,Rhs,1) ) CALL DCOPY(N,xp,1,Rhs,1)
!
!        Compute the sum of squares
!
         CALL DLASSQ(N,Rhs,1,Rdscal,Rdsum)
!
      ENDIF
!
!
!     End of DLATDF
!
      END SUBROUTINE DLATDF
