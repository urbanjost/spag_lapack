!*==zlatdf.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLATDF uses the LU factorization of the n-by-n matrix computed by sgetc2 and computes a contribution to the reciprocal Dif-estimate.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLATDF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatdf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatdf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatdf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV,
!                          JPIV )
!
!       .. Scalar Arguments ..
!       INTEGER            IJOB, LDZ, N
!       DOUBLE PRECISION   RDSCAL, RDSUM
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), JPIV( * )
!       COMPLEX*16         RHS( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLATDF computes the contribution to the reciprocal Dif-estimate
!> by solving for x in Z * x = b, where b is chosen such that the norm
!> of x is as large as possible. It is assumed that LU decomposition
!> of Z has been computed by ZGETC2. On entry RHS = f holds the
!> contribution from earlier solved sub-systems, and on return RHS = x.
!>
!> The factorization of Z returned by ZGETC2 has the form
!> Z = P * L * U * Q, where P and Q are permutation matrices. L is lower
!> triangular with unit diagonal elements and U is upper triangular.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IJOB
!> \verbatim
!>          IJOB is INTEGER
!>          IJOB = 2: First compute an approximative null-vector e
!>              of Z using ZGECON, e is normalized and solve for
!>              Zx = +-e - f with the sign giving the greater value of
!>              2-norm(x).  About 5 times as expensive as Default.
!>          IJOB .ne. 2: Local look ahead strategy where
!>              all entries of the r.h.s. b is chosen as either +1 or
!>              -1.  Default.
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
!>          Z is COMPLEX*16 array, dimension (LDZ, N)
!>          On entry, the LU part of the factorization of the n-by-n
!>          matrix Z computed by ZGETC2:  Z = P * L * U * Q
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
!>          RHS is COMPLEX*16 array, dimension (N).
!>          On entry, RHS contains contributions from other subsystems.
!>          On exit, RHS contains the solution of the subsystem with
!>          entries according to the value of IJOB (see above).
!> \endverbatim
!>
!> \param[in,out] RDSUM
!> \verbatim
!>          RDSUM is DOUBLE PRECISION
!>          On entry, the sum of squares of computed contributions to
!>          the Dif-estimate under computation by ZTGSYL, where the
!>          scaling factor RDSCAL (see below) has been factored out.
!>          On exit, the corresponding sum of squares updated with the
!>          contributions from the current sub-system.
!>          If TRANS = 'T' RDSUM is not touched.
!>          NOTE: RDSUM only makes sense when ZTGSY2 is called by CTGSYL.
!> \endverbatim
!>
!> \param[in,out] RDSCAL
!> \verbatim
!>          RDSCAL is DOUBLE PRECISION
!>          On entry, scaling factor used to prevent overflow in RDSUM.
!>          On exit, RDSCAL is updated w.r.t. the current contributions
!>          in RDSUM.
!>          If TRANS = 'T', RDSCAL is not touched.
!>          NOTE: RDSCAL only makes sense when ZTGSY2 is called by
!>          ZTGSYL.
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
!> \ingroup complex16OTHERauxiliary
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
!>   [1]   Bo Kagstrom and Lars Westin,
!>         Generalized Schur Methods with Condition Estimators for
!>         Solving the Generalized Sylvester Equation, IEEE Transactions
!>         on Automatic Control, Vol. 34, No. 7, July 1989, pp 745-751.
!>\n
!>   [2]   Peter Poromaa,
!>         On Efficient and Robust Estimators for the Separation
!>         between two Regular Matrix Pairs with Applications in
!>         Condition Estimation. Report UMINF-95.05, Department of
!>         Computing Science, Umea University, S-901 87 Umea, Sweden,
!>         1995.
!
!  =====================================================================
      SUBROUTINE ZLATDF(Ijob,N,Z,Ldz,Rhs,Rdsum,Rdscal,Ipiv,Jpiv)
      IMPLICIT NONE
!*--ZLATDF172
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
      COMPLEX*16 Rhs(*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER MAXDIM
      PARAMETER (MAXDIM=2)
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CONE
      PARAMETER (CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j , k
      DOUBLE PRECISION rtemp , scale , sminu , splus
      COMPLEX*16 bm , bp , pmone , temp
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION rwork(MAXDIM)
      COMPLEX*16 work(4*MAXDIM) , xm(MAXDIM) , xp(MAXDIM)
!     ..
!     .. External Subroutines ..
      EXTERNAL ZAXPY , ZCOPY , ZGECON , ZGESC2 , ZLASSQ , ZLASWP , ZSCAL
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DZASUM
      COMPLEX*16 ZDOTC
      EXTERNAL DZASUM , ZDOTC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , SQRT
!     ..
!     .. Executable Statements ..
!
      IF ( Ijob/=2 ) THEN
!
!        Apply permutations IPIV to RHS
!
         CALL ZLASWP(1,Rhs,Ldz,1,N-1,Ipiv,1)
!
!        Solve for L-part choosing RHS either to +1 or -1.
!
         pmone = -CONE
         DO j = 1 , N - 1
            bp = Rhs(j) + CONE
            bm = Rhs(j) - CONE
            splus = ONE
!
!           Lockahead for L- part RHS(1:N-1) = +-1
!           SPLUS and SMIN computed more efficiently than in BSOLVE[1].
!
            splus = splus + DBLE(ZDOTC(N-j,Z(j+1,j),1,Z(j+1,j),1))
            sminu = DBLE(ZDOTC(N-j,Z(j+1,j),1,Rhs(j+1),1))
            splus = splus*DBLE(Rhs(j))
            IF ( splus>sminu ) THEN
               Rhs(j) = bp
            ELSEIF ( sminu>splus ) THEN
               Rhs(j) = bm
            ELSE
!
!              In this case the updating sums are equal and we can
!              choose RHS(J) +1 or -1. The first time this happens we
!              choose -1, thereafter +1. This is a simple way to get
!              good estimates of matrices like Byers well-known example
!              (see [1]). (Not done in BSOLVE.)
!
               Rhs(j) = Rhs(j) + pmone
               pmone = CONE
            ENDIF
!
!           Compute the remaining r.h.s.
!
            temp = -Rhs(j)
            CALL ZAXPY(N-j,temp,Z(j+1,j),1,Rhs(j+1),1)
         ENDDO
!
!        Solve for U- part, lockahead for RHS(N) = +-1. This is not done
!        In BSOLVE and will hopefully give us a better estimate because
!        any ill-conditioning of the original matrix is transferred to U
!        and not to L. U(N, N) is an approximation to sigma_min(LU).
!
         CALL ZCOPY(N-1,Rhs,1,work,1)
         work(N) = Rhs(N) + CONE
         Rhs(N) = Rhs(N) - CONE
         splus = ZERO
         sminu = ZERO
         DO i = N , 1 , -1
            temp = CONE/Z(i,i)
            work(i) = work(i)*temp
            Rhs(i) = Rhs(i)*temp
            DO k = i + 1 , N
               work(i) = work(i) - work(k)*(Z(i,k)*temp)
               Rhs(i) = Rhs(i) - Rhs(k)*(Z(i,k)*temp)
            ENDDO
            splus = splus + ABS(work(i))
            sminu = sminu + ABS(Rhs(i))
         ENDDO
         IF ( splus>sminu ) CALL ZCOPY(N,work,1,Rhs,1)
!
!        Apply the permutations JPIV to the computed solution (RHS)
!
         CALL ZLASWP(1,Rhs,Ldz,1,N-1,Jpiv,-1)
!
!        Compute the sum of squares
!
         CALL ZLASSQ(N,Rhs,1,Rdscal,Rdsum)
         RETURN
      ENDIF
!
!     ENTRY IJOB = 2
!
!     Compute approximate nullvector XM of Z
!
      CALL ZGECON('I',N,Z,Ldz,ONE,rtemp,work,rwork,info)
      CALL ZCOPY(N,work(N+1),1,xm,1)
!
!     Compute RHS
!
      CALL ZLASWP(1,xm,Ldz,1,N-1,Ipiv,-1)
      temp = CONE/SQRT(ZDOTC(N,xm,1,xm,1))
      CALL ZSCAL(N,temp,xm,1)
      CALL ZCOPY(N,xm,1,xp,1)
      CALL ZAXPY(N,CONE,Rhs,1,xp,1)
      CALL ZAXPY(N,-CONE,xm,1,Rhs,1)
      CALL ZGESC2(N,Z,Ldz,Rhs,Ipiv,Jpiv,scale)
      CALL ZGESC2(N,Z,Ldz,xp,Ipiv,Jpiv,scale)
      IF ( DZASUM(N,xp,1)>DZASUM(N,Rhs,1) ) CALL ZCOPY(N,xp,1,Rhs,1)
!
!     Compute the sum of squares
!
      CALL ZLASSQ(N,Rhs,1,Rdscal,Rdsum)
!
!     End of ZLATDF
!
      END SUBROUTINE ZLATDF
