!*==dlaed2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAED2 used by sstedc. Merges eigenvalues and deflates secular equation. Used when the original matrix is tridiagonal.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAED2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED2( K, N, N1, D, Q, LDQ, INDXQ, RHO, Z, DLAMDA, W,
!                          Q2, INDX, INDXC, INDXP, COLTYP, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDQ, N, N1
!       DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
!       INTEGER            COLTYP( * ), INDX( * ), INDXC( * ), INDXP( * ),
!      $                   INDXQ( * )
!       DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),
!      $                   W( * ), Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED2 merges the two sets of eigenvalues together into a single
!> sorted set.  Then it tries to deflate the size of the problem.
!> There are two ways in which deflation can occur:  when two or more
!> eigenvalues are close together or if there is a tiny entry in the
!> Z vector.  For each such occurrence the order of the related secular
!> equation problem is reduced by one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] K
!> \verbatim
!>          K is INTEGER
!>         The number of non-deflated eigenvalues, and the order of the
!>         related secular equation. 0 <= K <=N.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>         The location of the last eigenvalue in the leading sub-matrix.
!>         min(1,N) <= N1 <= N/2.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>         On entry, D contains the eigenvalues of the two submatrices to
!>         be combined.
!>         On exit, D contains the trailing (N-K) updated eigenvalues
!>         (those which were deflated) sorted into increasing order.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
!>         On entry, Q contains the eigenvectors of two submatrices in
!>         the two square blocks with corners at (1,1), (N1,N1)
!>         and (N1+1, N1+1), (N,N).
!>         On exit, Q contains the trailing (N-K) updated eigenvectors
!>         (those which were deflated) in its last N-K columns.
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
!>         The permutation which separately sorts the two sub-problems
!>         in D into ascending order.  Note that elements in the second
!>         half of this permutation must first have N1 added to their
!>         values. Destroyed on exit.
!> \endverbatim
!>
!> \param[in,out] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>         On entry, the off-diagonal element associated with the rank-1
!>         cut which originally split the two submatrices which are now
!>         being recombined.
!>         On exit, RHO has been modified to the value required by
!>         DLAED3.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (N)
!>         On entry, Z contains the updating vector (the last
!>         row of the first sub-eigenvector matrix and the first row of
!>         the second sub-eigenvector matrix).
!>         On exit, the contents of Z have been destroyed by the updating
!>         process.
!> \endverbatim
!>
!> \param[out] DLAMDA
!> \verbatim
!>          DLAMDA is DOUBLE PRECISION array, dimension (N)
!>         A copy of the first K eigenvalues which will be used by
!>         DLAED3 to form the secular equation.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>         The first k values of the final deflation-altered z-vector
!>         which will be passed to DLAED3.
!> \endverbatim
!>
!> \param[out] Q2
!> \verbatim
!>          Q2 is DOUBLE PRECISION array, dimension (N1**2+(N-N1)**2)
!>         A copy of the first K eigenvectors which will be used by
!>         DLAED3 in a matrix multiply (DGEMM) to solve for the new
!>         eigenvectors.
!> \endverbatim
!>
!> \param[out] INDX
!> \verbatim
!>          INDX is INTEGER array, dimension (N)
!>         The permutation used to sort the contents of DLAMDA into
!>         ascending order.
!> \endverbatim
!>
!> \param[out] INDXC
!> \verbatim
!>          INDXC is INTEGER array, dimension (N)
!>         The permutation used to arrange the columns of the deflated
!>         Q matrix into three groups:  the first group contains non-zero
!>         elements only at and above N1, the second contains
!>         non-zero elements only below N1, and the third is dense.
!> \endverbatim
!>
!> \param[out] INDXP
!> \verbatim
!>          INDXP is INTEGER array, dimension (N)
!>         The permutation used to place deflated values of D at the end
!>         of the array.  INDXP(1:K) points to the nondeflated D-values
!>         and INDXP(K+1:N) points to the deflated eigenvalues.
!> \endverbatim
!>
!> \param[out] COLTYP
!> \verbatim
!>          COLTYP is INTEGER array, dimension (N)
!>         During execution, a label which will indicate which of the
!>         following types a column in the Q2 matrix is:
!>         1 : non-zero in the upper half only;
!>         2 : dense;
!>         3 : non-zero in the lower half only;
!>         4 : deflated.
!>         On exit, COLTYP(i) is the number of columns of type i,
!>         for i=1 to 4 only.
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
!> at Berkeley, USA \n
!>  Modified by Francoise Tisseur, University of Tennessee
!>
!  =====================================================================
      SUBROUTINE DLAED2(K,N,N1,D,Q,Ldq,Indxq,Rho,Z,Dlamda,W,Q2,Indx,    &
     &                  Indxc,Indxp,Coltyp,Info)
      IMPLICIT NONE
!*--DLAED2216
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , K , Ldq , N , N1
      DOUBLE PRECISION Rho
!     ..
!     .. Array Arguments ..
      INTEGER Coltyp(*) , Indx(*) , Indxc(*) , Indxp(*) , Indxq(*)
      DOUBLE PRECISION D(*) , Dlamda(*) , Q(Ldq,*) , Q2(*) , W(*) , Z(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION MONE , ZERO , ONE , TWO , EIGHT
      PARAMETER (MONE=-1.0D0,ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,EIGHT=8.0D0)
!     ..
!     .. Local Arrays ..
      INTEGER ctot(4) , psm(4)
!     ..
!     .. Local Scalars ..
      INTEGER ct , i , imax , iq1 , iq2 , j , jmax , js , k2 , n1p1 ,   &
     &        n2 , nj , pj
      DOUBLE PRECISION c , eps , s , t , tau , tol
!     ..
!     .. External Functions ..
      INTEGER IDAMAX
      DOUBLE PRECISION DLAMCH , DLAPY2
      EXTERNAL IDAMAX , DLAMCH , DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DLACPY , DLAMRG , DROT , DSCAL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Ldq<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( MIN(1,(N/2))>N1 .OR. (N/2)<N1 ) THEN
         Info = -3
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DLAED2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      n2 = N - N1
      n1p1 = N1 + 1
!
      IF ( Rho<ZERO ) CALL DSCAL(n2,MONE,Z(n1p1),1)
!
!     Normalize z so that norm(z) = 1.  Since z is the concatenation of
!     two normalized vectors, norm2(z) = sqrt(2).
!
      t = ONE/SQRT(TWO)
      CALL DSCAL(N,t,Z,1)
!
!     RHO = ABS( norm(z)**2 * RHO )
!
      Rho = ABS(TWO*Rho)
!
!     Sort the eigenvalues into increasing order
!
      DO i = n1p1 , N
         Indxq(i) = Indxq(i) + N1
      ENDDO
!
!     re-integrate the deflated parts from the last pass
!
      DO i = 1 , N
         Dlamda(i) = D(Indxq(i))
      ENDDO
      CALL DLAMRG(N1,n2,Dlamda,1,1,Indxc)
      DO i = 1 , N
         Indx(i) = Indxq(Indxc(i))
      ENDDO
!
!     Calculate the allowable deflation tolerance
!
      imax = IDAMAX(N,Z,1)
      jmax = IDAMAX(N,D,1)
      eps = DLAMCH('Epsilon')
      tol = EIGHT*eps*MAX(ABS(D(jmax)),ABS(Z(imax)))
!
!     If the rank-1 modifier is small enough, no more needs to be done
!     except to reorganize Q so that its columns correspond with the
!     elements in D.
!
      IF ( Rho*ABS(Z(imax))<=tol ) THEN
         K = 0
         iq2 = 1
         DO j = 1 , N
            i = Indx(j)
            CALL DCOPY(N,Q(1,i),1,Q2(iq2),1)
            Dlamda(j) = D(i)
            iq2 = iq2 + N
         ENDDO
         CALL DLACPY('A',N,N,Q2,N,Q,Ldq)
         CALL DCOPY(N,Dlamda,1,D,1)
         GOTO 99999
      ENDIF
!
!     If there are multiple eigenvalues then the problem deflates.  Here
!     the number of equal eigenvalues are found.  As each equal
!     eigenvalue is found, an elementary reflector is computed to rotate
!     the corresponding eigensubspace so that the corresponding
!     components of Z are zero in this new basis.
!
      DO i = 1 , N1
         Coltyp(i) = 1
      ENDDO
      DO i = n1p1 , N
         Coltyp(i) = 3
      ENDDO
!
!
      K = 0
      k2 = N + 1
      DO j = 1 , N
         nj = Indx(j)
         IF ( Rho*ABS(Z(nj))<=tol ) THEN
!
!           Deflate due to small z component.
!
            k2 = k2 - 1
            Coltyp(nj) = 4
            Indxp(k2) = nj
            IF ( j==N ) GOTO 100
         ELSE
            pj = nj
            EXIT
         ENDIF
      ENDDO
      DO
         j = j + 1
         nj = Indx(j)
         IF ( j>N ) EXIT
         IF ( Rho*ABS(Z(nj))<=tol ) THEN
!
!        Deflate due to small z component.
!
            k2 = k2 - 1
            Coltyp(nj) = 4
            Indxp(k2) = nj
         ELSE
!
!        Check if eigenvalues are close enough to allow deflation.
!
            s = Z(pj)
            c = Z(nj)
!
!        Find sqrt(a**2+b**2) without overflow or
!        destructive underflow.
!
            tau = DLAPY2(c,s)
            t = D(nj) - D(pj)
            c = c/tau
            s = -s/tau
            IF ( ABS(t*c*s)<=tol ) THEN
!
!           Deflation is possible.
!
               Z(nj) = tau
               Z(pj) = ZERO
               IF ( Coltyp(nj)/=Coltyp(pj) ) Coltyp(nj) = 2
               Coltyp(pj) = 4
               CALL DROT(N,Q(1,pj),1,Q(1,nj),1,c,s)
               t = D(pj)*c**2 + D(nj)*s**2
               D(nj) = D(pj)*s**2 + D(nj)*c**2
               D(pj) = t
               k2 = k2 - 1
               i = 1
               DO
                  IF ( k2+i>N ) THEN
                     Indxp(k2+i-1) = pj
                  ELSEIF ( D(pj)<D(Indxp(k2+i)) ) THEN
                     Indxp(k2+i-1) = Indxp(k2+i)
                     Indxp(k2+i) = pj
                     i = i + 1
                     CYCLE
                  ELSE
                     Indxp(k2+i-1) = pj
                  ENDIF
                  pj = nj
                  EXIT
               ENDDO
            ELSE
               K = K + 1
               Dlamda(K) = D(pj)
               W(K) = Z(pj)
               Indxp(K) = pj
               pj = nj
            ENDIF
         ENDIF
      ENDDO
!
!     Record the last eigenvalue.
!
 100  K = K + 1
      Dlamda(K) = D(pj)
      W(K) = Z(pj)
      Indxp(K) = pj
!
!     Count up the total number of the various types of columns, then
!     form a permutation which positions the four column types into
!     four uniform groups (although one or more of these groups may be
!     empty).
!
      DO j = 1 , 4
         ctot(j) = 0
      ENDDO
      DO j = 1 , N
         ct = Coltyp(j)
         ctot(ct) = ctot(ct) + 1
      ENDDO
!
!     PSM(*) = Position in SubMatrix (of types 1 through 4)
!
      psm(1) = 1
      psm(2) = 1 + ctot(1)
      psm(3) = psm(2) + ctot(2)
      psm(4) = psm(3) + ctot(3)
      K = N - ctot(4)
!
!     Fill out the INDXC array so that the permutation which it induces
!     will place all type-1 columns first, all type-2 columns next,
!     then all type-3's, and finally all type-4's.
!
      DO j = 1 , N
         js = Indxp(j)
         ct = Coltyp(js)
         Indx(psm(ct)) = js
         Indxc(psm(ct)) = j
         psm(ct) = psm(ct) + 1
      ENDDO
!
!     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
!     and Q2 respectively.  The eigenvalues/vectors which were not
!     deflated go into the first K slots of DLAMDA and Q2 respectively,
!     while those which were deflated go into the last N - K slots.
!
      i = 1
      iq1 = 1
      iq2 = 1 + (ctot(1)+ctot(2))*N1
      DO j = 1 , ctot(1)
         js = Indx(i)
         CALL DCOPY(N1,Q(1,js),1,Q2(iq1),1)
         Z(i) = D(js)
         i = i + 1
         iq1 = iq1 + N1
      ENDDO
!
      DO j = 1 , ctot(2)
         js = Indx(i)
         CALL DCOPY(N1,Q(1,js),1,Q2(iq1),1)
         CALL DCOPY(n2,Q(N1+1,js),1,Q2(iq2),1)
         Z(i) = D(js)
         i = i + 1
         iq1 = iq1 + N1
         iq2 = iq2 + n2
      ENDDO
!
      DO j = 1 , ctot(3)
         js = Indx(i)
         CALL DCOPY(n2,Q(N1+1,js),1,Q2(iq2),1)
         Z(i) = D(js)
         i = i + 1
         iq2 = iq2 + n2
      ENDDO
!
      iq1 = iq2
      DO j = 1 , ctot(4)
         js = Indx(i)
         CALL DCOPY(N,Q(1,js),1,Q2(iq2),1)
         iq2 = iq2 + N
         Z(i) = D(js)
         i = i + 1
      ENDDO
!
!     The deflated eigenvalues and their corresponding vectors go back
!     into the last N - K slots of D and Q respectively.
!
      IF ( K<N ) THEN
         CALL DLACPY('A',N,ctot(4),Q2(iq1),N,Q(1,K+1),Ldq)
         CALL DCOPY(N-K,Z(K+1),1,D(K+1),1)
      ENDIF
!
!     Copy CTOT into COLTYP for referencing in DLAED3.
!
      DO j = 1 , 4
         Coltyp(j) = ctot(j)
      ENDDO
!
!
!     End of DLAED2
!
99999 END SUBROUTINE DLAED2
