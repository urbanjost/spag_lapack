!*==cgbt01.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cgbt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cgbt01( M, N, KL, KU, A, LDA, AFAC, LDAFAC, IPIV, WORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            KL, KU, LDA, LDAFAC, M, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cgbt01 reconstructs a band matrix  A  from its L*U factorization and
!> computes the residual:
!>    norm(L*U - A) / ( N * norm(A) * EPS ),
!> where EPS is the machine epsilon.
!>
!> The expression L*U - A is computed one column at a time, so A and
!> AFAC are not modified.
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
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The original matrix A in band storage, stored in rows 1 to
!>          KL+KU+1.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER.
!>          The leading dimension of the array A.  LDA >= max(1,KL+KU+1).
!> \endverbatim
!>
!> \param[in] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (LDAFAC,N)
!>          The factored form of the matrix A.  AFAC contains the banded
!>          factors L and U from the L*U factorization, as computed by
!>          CGBTRF.  U is stored as an upper triangular band matrix with
!>          KL+KU superdiagonals in rows 1 to KL+KU+1, and the
!>          multipliers used during the factorization are stored in rows
!>          KL+KU+2 to 2*KL+KU+1.  See CGBTRF for further details.
!> \endverbatim
!>
!> \param[in] LDAFAC
!> \verbatim
!>          LDAFAC is INTEGER
!>          The leading dimension of the array AFAC.
!>          LDAFAC >= max(1,2*KL*KU+1).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices from CGBTRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*KL+KU+1)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cgbt01(M,N,Kl,Ku,A,Lda,Afac,Ldafac,Ipiv,Work,Resid)
      IMPLICIT NONE
!*--CGBT01129
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Kl , Ku , Lda , Ldafac , M , N
      REAL Resid
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX A(Lda,*) , Afac(Ldafac,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , i1 , i2 , il , ip , iw , j , jl , ju , jua , kd , lenj
      REAL anorm , eps
      COMPLEX t
!     ..
!     .. External Functions ..
      REAL SCASUM , SLAMCH
EXTERNAL    SCASUM,   SLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    CAXPY,   CCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
!     Quick exit if M = 0 or N = 0.
!
      Resid = ZERO
      IF ( M<=0 .OR. N<=0 ) RETURN
!
!     Determine EPS and the norm of A.
!
      eps = SLAMCH('Epsilon')
      kd = Ku + 1
      anorm = ZERO
      DO j = 1 , N
         i1 = MAX(kd+1-j,1)
         i2 = MIN(kd+M-j,Kl+kd)
         IF ( i2>=i1 ) anorm = MAX(anorm,SCASUM(i2-i1+1,A(i1,j),1))
      ENDDO
!
!     Compute one column at a time of L*U - A.
!
      kd = Kl + Ku + 1
      DO j = 1 , N
!
!        Copy the J-th column of U to WORK.
!
         ju = MIN(Kl+Ku,j-1)
         jl = MIN(Kl,M-j)
         lenj = MIN(M,j) - j + ju + 1
         IF ( lenj>0 ) THEN
            CALL CCOPY(lenj,Afac(kd-ju,j),1,Work,1)
            DO i = lenj + 1 , ju + jl + 1
               Work(i) = ZERO
            ENDDO
!
!           Multiply by the unit lower triangular matrix L.  Note that L
!           is stored as a product of transformations and permutations.
!
            DO i = MIN(M-1,j) , j - ju , -1
               il = MIN(Kl,M-i)
               IF ( il>0 ) THEN
                  iw = i - j + ju + 1
                  t = Work(iw)
                  CALL CAXPY(il,t,Afac(kd+1,i),1,Work(iw+1),1)
                  ip = Ipiv(i)
                  IF ( i/=ip ) THEN
                     ip = ip - j + ju + 1
                     Work(iw) = Work(ip)
                     Work(ip) = t
                  ENDIF
               ENDIF
            ENDDO
!
!           Subtract the corresponding column of A.
!
            jua = MIN(ju,Ku)
            IF ( jua+jl+1>0 )                                           &
     &           CALL CAXPY(jua+jl+1,-CMPLX(ONE),A(Ku+1-jua,j),1,       &
     &           Work(ju+1-jua),1)
!
!           Compute the 1-norm of the column.
!
            Resid = MAX(Resid,SCASUM(ju+jl+1,Work,1))
         ENDIF
      ENDDO
!
!     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSE
         Resid = ((Resid/REAL(N))/anorm)/eps
      ENDIF
!
!
!     End of cgbt01
!
      END SUBROUTINE cgbt01
