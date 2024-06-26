!*==cunt03.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b cunt03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cunt03( RC, MU, MV, N, K, U, LDU, V, LDV, WORK, LWORK,
!                          RWORK, RESULT, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER*( * )    RC
!       INTEGER            INFO, K, LDU, LDV, LWORK, MU, MV, N
!       REAL               RESULT
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            U( LDU, * ), V( LDV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cunt03 compares two unitary matrices U and V to see if their
!> corresponding rows or columns span the same spaces.  The rows are
!> checked if RC = 'R', and the columns are checked if RC = 'C'.
!>
!> RESULT is the maximum of
!>
!>    | V*V' - I | / ( MV ulp ), if RC = 'R', or
!>
!>    | V'*V - I | / ( MV ulp ), if RC = 'C',
!>
!> and the maximum over rows (or columns) 1 to K of
!>
!>    | U(i) - S*V(i) |/ ( N ulp )
!>
!> where abs(S) = 1 (chosen to minimize the expression), U(i) is the
!> i-th row (column) of U, and V(i) is the i-th row (column) of V.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] RC
!> \verbatim
!>          RC is CHARACTER*1
!>          If RC = 'R' the rows of U and V are to be compared.
!>          If RC = 'C' the columns of U and V are to be compared.
!> \endverbatim
!>
!> \param[in] MU
!> \verbatim
!>          MU is INTEGER
!>          The number of rows of U if RC = 'R', and the number of
!>          columns if RC = 'C'.  If MU = 0 cunt03 does nothing.
!>          MU must be at least zero.
!> \endverbatim
!>
!> \param[in] MV
!> \verbatim
!>          MV is INTEGER
!>          The number of rows of V if RC = 'R', and the number of
!>          columns if RC = 'C'.  If MV = 0 cunt03 does nothing.
!>          MV must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          If RC = 'R', the number of columns in the matrices U and V,
!>          and if RC = 'C', the number of rows in U and V.  If N = 0
!>          cunt03 does nothing.  N must be at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of rows or columns of U and V to compare.
!>          0 <= K <= max(MU,MV).
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is COMPLEX array, dimension (LDU,N)
!>          The first matrix to compare.  If RC = 'R', U is MU by N, and
!>          if RC = 'C', U is N by MU.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  If RC = 'R', LDU >= max(1,MU),
!>          and if RC = 'C', LDU >= max(1,N).
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX array, dimension (LDV,N)
!>          The second matrix to compare.  If RC = 'R', V is MV by N, and
!>          if RC = 'C', V is N by MV.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of V.  If RC = 'R', LDV >= max(1,MV),
!>          and if RC = 'C', LDV >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  For best performance, LWORK
!>          should be at least N*N if RC = 'C' or M*M if RC = 'R', but
!>          the tests will be done even if LWORK is 0.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (max(MV,N))
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL
!>          The value computed by the test described above.  RESULT is
!>          limited to 1/ulp to avoid overflow.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          0  indicates a successful exit
!>          -k indicates the k-th parameter had an illegal value
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE cunt03(Rc,Mu,Mv,N,K,U,Ldu,V,Ldv,Work,Lwork,Rwork,      &
     &                  Result,Info)
      IMPLICIT NONE
!*--CUNT03166
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*(*) Rc
      INTEGER Info , K , Ldu , Ldv , Lwork , Mu , Mv , N
      REAL Result
!     ..
!     .. Array Arguments ..
      REAL Rwork(*)
      COMPLEX U(Ldu,*) , V(Ldv,*) , Work(*)
!     ..
!
!  =====================================================================
!
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER i , irc , j , lmx
      REAL res1 , res2 , ulp
      COMPLEX s , su , sv
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ICAMAX
      REAL SLAMCH
      EXTERNAL LSAME , ICAMAX , SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CMPLX , MAX , MIN , REAL
!     ..
!     .. External Subroutines ..
     
!     ..
!     .. Executable Statements ..
!
!     Check inputs
!
      Info = 0
      IF ( LSAME(Rc,'R') ) THEN
         irc = 0
      ELSEIF ( LSAME(Rc,'C') ) THEN
         irc = 1
      ELSE
         irc = -1
      ENDIF
      IF ( irc==-1 ) THEN
         Info = -1
      ELSEIF ( Mu<0 ) THEN
         Info = -2
      ELSEIF ( Mv<0 ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( K<0 .OR. K>MAX(Mu,Mv) ) THEN
         Info = -5
      ELSEIF ( (irc==0 .AND. Ldu<MAX(1,Mu)) .OR.                        &
     &         (irc==1 .AND. Ldu<MAX(1,N)) ) THEN
         Info = -7
      ELSEIF ( (irc==0 .AND. Ldv<MAX(1,Mv)) .OR.                        &
     &         (irc==1 .AND. Ldv<MAX(1,N)) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL xerbla('cunt03',-Info)
         RETURN
      ENDIF
!
!     Initialize result
!
      Result = ZERO
      IF ( Mu==0 .OR. Mv==0 .OR. N==0 ) RETURN
!
!     Machine constants
!
      ulp = SLAMCH('Precision')
!
      IF ( irc==0 ) THEN
!
!        Compare rows
!
         res1 = ZERO
         DO i = 1 , K
            lmx = ICAMAX(N,U(i,1),Ldu)
            IF ( V(i,lmx)==CMPLX(ZERO) ) THEN
               sv = ONE
            ELSE
               sv = ABS(V(i,lmx))/V(i,lmx)
            ENDIF
            IF ( U(i,lmx)==CMPLX(ZERO) ) THEN
               su = ONE
            ELSE
               su = ABS(U(i,lmx))/U(i,lmx)
            ENDIF
            s = sv/su
            DO j = 1 , N
               res1 = MAX(res1,ABS(U(i,j)-s*V(i,j)))
            ENDDO
         ENDDO
         res1 = res1/(REAL(N)*ulp)
!
!        Compute orthogonality of rows of V.
!
         CALL cunt01('Rows',Mv,N,V,Ldv,Work,Lwork,Rwork,res2)
!
      ELSE
!
!        Compare columns
!
         res1 = ZERO
         DO i = 1 , K
            lmx = ICAMAX(N,U(1,i),1)
            IF ( V(lmx,i)==CMPLX(ZERO) ) THEN
               sv = ONE
            ELSE
               sv = ABS(V(lmx,i))/V(lmx,i)
            ENDIF
            IF ( U(lmx,i)==CMPLX(ZERO) ) THEN
               su = ONE
            ELSE
               su = ABS(U(lmx,i))/U(lmx,i)
            ENDIF
            s = sv/su
            DO j = 1 , N
               res1 = MAX(res1,ABS(U(j,i)-s*V(j,i)))
            ENDDO
         ENDDO
         res1 = res1/(REAL(N)*ulp)
!
!        Compute orthogonality of columns of V.
!
         CALL cunt01('Columns',N,Mv,V,Ldv,Work,Lwork,Rwork,res2)
      ENDIF
!
      Result = MIN(MAX(res1,res2),ONE/ulp)
!
!     End of cunt03
!
      END SUBROUTINE cunt03
