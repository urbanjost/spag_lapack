!*==dlaed3.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAED3 used by sstedc. Finds the roots of the secular equation and updates the eigenvectors. Used when the original matrix is tridiagonal.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAED3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMDA, Q2, INDX,
!                          CTOT, W, S, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDQ, N, N1
!       DOUBLE PRECISION   RHO
!       ..
!       .. Array Arguments ..
!       INTEGER            CTOT( * ), INDX( * )
!       DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),
!      $                   S( * ), W( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED3 finds the roots of the secular equation, as defined by the
!> values in D, W, and RHO, between 1 and K.  It makes the
!> appropriate calls to DLAED4 and then updates the eigenvectors by
!> multiplying the matrix of eigenvectors of the pair of eigensystems
!> being combined by the matrix of eigenvectors of the K-by-K system
!> which is solved here.
!>
!> This code makes very mild assumptions about floating point
!> arithmetic. It will work on machines with a guard digit in
!> add/subtract, or on those binary machines without guard digits
!> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!> It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of terms in the rational function to be solved by
!>          DLAED4.  K >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns in the Q matrix.
!>          N >= K (deflation may result in N>K).
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>          The location of the last eigenvalue in the leading submatrix.
!>          min(1,N) <= N1 <= N/2.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          D(I) contains the updated eigenvalues for
!>          1 <= I <= K.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>          Initially the first K columns are used as workspace.
!>          On output the columns 1 to K contain
!>          the updated eigenvectors.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is DOUBLE PRECISION
!>          The value of the parameter in the rank one update equation.
!>          RHO >= 0 required.
!> \endverbatim
!>
!> \param[in,out] DLAMDA
!> \verbatim
!>          DLAMDA is DOUBLE PRECISION array, dimension (K)
!>          The first K elements of this array contain the old roots
!>          of the deflated updating problem.  These are the poles
!>          of the secular equation. May be changed on output by
!>          having lowest order bit set to zero on Cray X-MP, Cray Y-MP,
!>          Cray-2, or Cray C-90, as described above.
!> \endverbatim
!>
!> \param[in] Q2
!> \verbatim
!>          Q2 is DOUBLE PRECISION array, dimension (LDQ2*N)
!>          The first K columns of this matrix contain the non-deflated
!>          eigenvectors for the split problem.
!> \endverbatim
!>
!> \param[in] INDX
!> \verbatim
!>          INDX is INTEGER array, dimension (N)
!>          The permutation used to arrange the columns of the deflated
!>          Q matrix into three groups (see DLAED2).
!>          The rows of the eigenvectors found by DLAED4 must be likewise
!>          permuted before the matrix multiply can take place.
!> \endverbatim
!>
!> \param[in] CTOT
!> \verbatim
!>          CTOT is INTEGER array, dimension (4)
!>          A count of the total number of the various types of columns
!>          in Q, as described in INDX.  The fourth column type is any
!>          column which has been deflated.
!> \endverbatim
!>
!> \param[in,out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (K)
!>          The first K elements of this array contain the components
!>          of the deflation-adjusted updating vector. Destroyed on
!>          output.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (N1 + 1)*K
!>          Will contain the eigenvectors of the repaired matrix which
!>          will be multiplied by the previously accumulated eigenvectors
!>          to update the system.
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
!> \date June 2017
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
      SUBROUTINE DLAED3(K,N,N1,D,Q,Ldq,Rho,Dlamda,Q2,Indx,Ctot,W,S,Info)
      IMPLICIT NONE
!*--DLAED3188
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Info , K , Ldq , N , N1
      DOUBLE PRECISION Rho
!     ..
!     .. Array Arguments ..
      INTEGER Ctot(*) , Indx(*)
      DOUBLE PRECISION D(*) , Dlamda(*) , Q(Ldq,*) , Q2(*) , S(*) , W(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ii , iq2 , j , n12 , n2 , n23
      DOUBLE PRECISION temp
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMC3 , DNRM2
      EXTERNAL DLAMC3 , DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DGEMM , DLACPY , DLAED4 , DLASET , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SIGN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( K<0 ) THEN
         Info = -1
      ELSEIF ( N<K ) THEN
         Info = -2
      ELSEIF ( Ldq<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DLAED3',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( K==0 ) RETURN
!
!     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DLAMDA(I) if it is 1; this makes the subsequent
!     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DLAMDA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DLAMDA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO i = 1 , K
         Dlamda(i) = DLAMC3(Dlamda(i),Dlamda(i)) - Dlamda(i)
      ENDDO
!
      DO j = 1 , K
         CALL DLAED4(K,j,Dlamda,W,Q(1,j),Rho,D(j),Info)
!
!        If the zero finder fails, the computation is terminated.
!
         IF ( Info/=0 ) GOTO 99999
      ENDDO
!
      IF ( K/=1 ) THEN
         IF ( K==2 ) THEN
            DO j = 1 , K
               W(1) = Q(1,j)
               W(2) = Q(2,j)
               ii = Indx(1)
               Q(1,j) = W(ii)
               ii = Indx(2)
               Q(2,j) = W(ii)
            ENDDO
            GOTO 100
         ENDIF
!
!     Compute updated W.
!
         CALL DCOPY(K,W,1,S,1)
!
!     Initialize W(I) = Q(I,I)
!
         CALL DCOPY(K,Q,Ldq+1,W,1)
         DO j = 1 , K
            DO i = 1 , j - 1
               W(i) = W(i)*(Q(i,j)/(Dlamda(i)-Dlamda(j)))
            ENDDO
            DO i = j + 1 , K
               W(i) = W(i)*(Q(i,j)/(Dlamda(i)-Dlamda(j)))
            ENDDO
         ENDDO
         DO i = 1 , K
            W(i) = SIGN(SQRT(-W(i)),S(i))
         ENDDO
!
!     Compute eigenvectors of the modified rank-1 modification.
!
         DO j = 1 , K
            DO i = 1 , K
               S(i) = W(i)/Q(i,j)
            ENDDO
            temp = DNRM2(K,S,1)
            DO i = 1 , K
               ii = Indx(i)
               Q(i,j) = S(ii)/temp
            ENDDO
         ENDDO
      ENDIF
!
!     Compute the updated eigenvectors.
!
!
 100  n2 = N - N1
      n12 = Ctot(1) + Ctot(2)
      n23 = Ctot(2) + Ctot(3)
!
      CALL DLACPY('A',n23,K,Q(Ctot(1)+1,1),Ldq,S,n23)
      iq2 = N1*n12 + 1
      IF ( n23/=0 ) THEN
         CALL DGEMM('N','N',n2,K,n23,ONE,Q2(iq2),n2,S,n23,ZERO,Q(N1+1,1)&
     &              ,Ldq)
      ELSE
         CALL DLASET('A',n2,K,ZERO,ZERO,Q(N1+1,1),Ldq)
      ENDIF
!
      CALL DLACPY('A',n12,K,Q,Ldq,S,n12)
      IF ( n12/=0 ) THEN
         CALL DGEMM('N','N',N1,K,n12,ONE,Q2,N1,S,n12,ZERO,Q,Ldq)
      ELSE
         CALL DLASET('A',N1,K,ZERO,ZERO,Q(1,1),Ldq)
      ENDIF
!
!
!
!     End of DLAED3
!
99999 END SUBROUTINE DLAED3
