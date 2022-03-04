!*==zget01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZGET01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGET01( M, N, A, LDA, AFAC, LDAFAC, IPIV, RWORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDAFAC, M, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGET01 reconstructs a matrix A from its L*U factorization and
!> computes the residual
!>    norm(L*U - A) / ( N * norm(A) * EPS ),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The original M x N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension (LDAFAC,N)
!>          The factored form of the matrix A.  AFAC contains the factors
!>          L and U from the L*U factorization as computed by ZGETRF.
!>          Overwritten with the reconstructed matrix, and then with the
!>          difference L*U - A.
!> \endverbatim
!>
!> \param[in] LDAFAC
!> \verbatim
!>          LDAFAC is INTEGER
!>          The leading dimension of the array AFAC.  LDAFAC >= max(1,M).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from ZGETRF.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          norm(L*U - A) / ( N * norm(A) * EPS )
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZGET01(M,N,A,Lda,Afac,Ldafac,Ipiv,Rwork,Resid)
      IMPLICIT NONE
!*--ZGET01111
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldafac , M , N
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(Lda,*) , Afac(Ldafac,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CONE
      PARAMETER (CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , j , k
      DOUBLE PRECISION anorm , eps
      COMPLEX*16 t
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE
      COMPLEX*16 ZDOTU
      EXTERNAL DLAMCH , ZLANGE , ZDOTU
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMV , ZLASWP , ZSCAL , ZTRMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MIN
!     ..
!     .. Executable Statements ..
!
!     Quick exit if M = 0 or N = 0.
!
      IF ( M<=0 .OR. N<=0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
!
!     Determine EPS and the norm of A.
!
      eps = DLAMCH('Epsilon')
      anorm = ZLANGE('1',M,N,A,Lda,Rwork)
!
!     Compute the product L*U and overwrite AFAC with the result.
!     A column at a time of the product is obtained, starting with
!     column N.
!
      DO k = N , 1 , -1
         IF ( k>M ) THEN
            CALL ZTRMV('Lower','No transpose','Unit',M,Afac,Ldafac,     &
     &                 Afac(1,k),1)
         ELSE
!
!           Compute elements (K+1:M,K)
!
            t = Afac(k,k)
            IF ( k+1<=M ) THEN
               CALL ZSCAL(M-k,t,Afac(k+1,k),1)
               CALL ZGEMV('No transpose',M-k,k-1,CONE,Afac(k+1,1),      &
     &                    Ldafac,Afac(1,k),1,CONE,Afac(k+1,k),1)
            ENDIF
!
!           Compute the (K,K) element
!
            Afac(k,k) = t + ZDOTU(k-1,Afac(k,1),Ldafac,Afac(1,k),1)
!
!           Compute elements (1:K-1,K)
!
            CALL ZTRMV('Lower','No transpose','Unit',k-1,Afac,Ldafac,   &
     &                 Afac(1,k),1)
         ENDIF
      ENDDO
      CALL ZLASWP(N,Afac,Ldafac,1,MIN(M,N),Ipiv,-1)
!
!     Compute the difference  L*U - A  and store in AFAC.
!
      DO j = 1 , N
         DO i = 1 , M
            Afac(i,j) = Afac(i,j) - A(i,j)
         ENDDO
      ENDDO
!
!     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
!
      Resid = ZLANGE('1',M,N,Afac,Ldafac,Rwork)
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSE
         Resid = ((Resid/DBLE(N))/anorm)/eps
      ENDIF
!
!
!     End of ZGET01
!
      END SUBROUTINE ZGET01
