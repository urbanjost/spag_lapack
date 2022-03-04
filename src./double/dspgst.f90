!*==dspgst.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DSPGST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSPGST + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dspgst.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dspgst.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dspgst.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSPGST( ITYPE, UPLO, N, AP, BP, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, ITYPE, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( * ), BP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSPGST reduces a real symmetric-definite generalized eigenproblem
!> to standard form, using packed storage.
!>
!> If ITYPE = 1, the problem is A*x = lambda*B*x,
!> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
!>
!> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
!>
!> B must have been previously factorized as U**T*U or L*L**T by DPPTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
!>          = 2 or 3: compute U*A*U**T or L**T*A*L.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored and B is factored as
!>                  U**T*U;
!>          = 'L':  Lower triangle of A is stored and B is factored as
!>                  L*L**T.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the symmetric matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!>
!>          On exit, if INFO = 0, the transformed matrix, stored in the
!>          same format as A.
!> \endverbatim
!>
!> \param[in] BP
!> \verbatim
!>          BP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          The triangular factor from the Cholesky factorization of B,
!>          stored in the same format as A, as returned by DPPTRF.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DSPGST(Itype,Uplo,N,Ap,Bp,Info)
      USE F77KINDS                        
      USE S_DAXPY
      USE S_DDOT
      USE S_DSCAL
      USE S_DSPMV
      USE S_DSPR2
      USE S_DTPMV
      USE S_DTPSV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DSPGST127
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D0 , HALF = 0.5D0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Itype
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: Bp
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ajj , akk , bjj , bkk , ct
      INTEGER :: j , j1 , j1j1 , jj , k , k1 , k1k1 , kk
      LOGICAL :: upper
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
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( Itype<1 .OR. Itype>3 ) THEN
         Info = -1
      ELSEIF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSPGST',-Info)
         RETURN
      ENDIF
!
      IF ( Itype==1 ) THEN
         IF ( upper ) THEN
!
!           Compute inv(U**T)*A*inv(U)
!
!           J1 and JJ are the indices of A(1,j) and A(j,j)
!
            jj = 0
            DO j = 1 , N
               j1 = jj + 1
               jj = jj + j
!
!              Compute the j-th column of the upper triangle of A
!
               bjj = Bp(jj)
               CALL DTPSV(Uplo,'Transpose','Nonunit',j,Bp,Ap(j1),1)
               CALL DSPMV(Uplo,j-1,-ONE,Ap,Bp(j1),1,ONE,Ap(j1),1)
               CALL DSCAL(j-1,ONE/bjj,Ap(j1),1)
               Ap(jj) = (Ap(jj)-DDOT(j-1,Ap(j1),1,Bp(j1),1))/bjj
            ENDDO
         ELSE
!
!           Compute inv(L)*A*inv(L**T)
!
!           KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)
!
            kk = 1
            DO k = 1 , N
               k1k1 = kk + N - k + 1
!
!              Update the lower triangle of A(k:n,k:n)
!
               akk = Ap(kk)
               bkk = Bp(kk)
               akk = akk/bkk**2
               Ap(kk) = akk
               IF ( k<N ) THEN
                  CALL DSCAL(N-k,ONE/bkk,Ap(kk+1),1)
                  ct = -HALF*akk
                  CALL DAXPY(N-k,ct,Bp(kk+1),1,Ap(kk+1),1)
                  CALL DSPR2(Uplo,N-k,-ONE,Ap(kk+1),1,Bp(kk+1),1,       &
     &                       Ap(k1k1))
                  CALL DAXPY(N-k,ct,Bp(kk+1),1,Ap(kk+1),1)
                  CALL DTPSV(Uplo,'No transpose','Non-unit',N-k,Bp(k1k1)&
     &                       ,Ap(kk+1),1)
               ENDIF
               kk = k1k1
            ENDDO
         ENDIF
      ELSEIF ( upper ) THEN
!
!           Compute U*A*U**T
!
!           K1 and KK are the indices of A(1,k) and A(k,k)
!
         kk = 0
         DO k = 1 , N
            k1 = kk + 1
            kk = kk + k
!
!              Update the upper triangle of A(1:k,1:k)
!
            akk = Ap(kk)
            bkk = Bp(kk)
            CALL DTPMV(Uplo,'No transpose','Non-unit',k-1,Bp,Ap(k1),1)
            ct = HALF*akk
            CALL DAXPY(k-1,ct,Bp(k1),1,Ap(k1),1)
            CALL DSPR2(Uplo,k-1,ONE,Ap(k1),1,Bp(k1),1,Ap)
            CALL DAXPY(k-1,ct,Bp(k1),1,Ap(k1),1)
            CALL DSCAL(k-1,bkk,Ap(k1),1)
            Ap(kk) = akk*bkk**2
         ENDDO
      ELSE
!
!           Compute L**T *A*L
!
!           JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)
!
         jj = 1
         DO j = 1 , N
            j1j1 = jj + N - j + 1
!
!              Compute the j-th column of the lower triangle of A
!
            ajj = Ap(jj)
            bjj = Bp(jj)
            Ap(jj) = ajj*bjj + DDOT(N-j,Ap(jj+1),1,Bp(jj+1),1)
            CALL DSCAL(N-j,bjj,Ap(jj+1),1)
            CALL DSPMV(Uplo,N-j,ONE,Ap(j1j1),Bp(jj+1),1,ONE,Ap(jj+1),1)
            CALL DTPMV(Uplo,'Transpose','Non-unit',N-j+1,Bp(jj),Ap(jj), &
     &                 1)
            jj = j1j1
         ENDDO
      ENDIF
!
!     End of DSPGST
!
      END SUBROUTINE DSPGST
