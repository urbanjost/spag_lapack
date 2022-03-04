!*==zgbtrs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZGBTRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGBTRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgbtrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgbtrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgbtrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         AB( LDAB, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGBTRS solves a system of linear equations
!>    A * X = B,  A**T * X = B,  or  A**H * X = B
!> with a general band matrix A using the LU factorization computed
!> by ZGBTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations.
!>          = 'N':  A * X = B     (No transpose)
!>          = 'T':  A**T * X = B  (Transpose)
!>          = 'C':  A**H * X = B  (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
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
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          Details of the LU factorization of the band matrix A, as
!>          computed by ZGBTRF.  U is stored as an upper triangular band
!>          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!>          the multipliers used during the factorization are stored in
!>          rows KL+KU+2 to 2*KL+KU+1.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices; for 1 <= i <= N, row i of the matrix was
!>          interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the right hand side matrix B.
!>          On exit, the solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
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
!> \ingroup complex16GBcomputational
!
!  =====================================================================
      SUBROUTINE ZGBTRS(Trans,N,Kl,Ku,Nrhs,Ab,Ldab,Ipiv,B,Ldb,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZGEMV
      USE S_ZGERU
      USE S_ZLACGV
      USE S_ZSWAP
      USE S_ZTBSV
      IMPLICIT NONE
!*--ZGBTRS149
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , j , kd , l , lm
      LOGICAL :: lnoti , notran
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
!     .. External Functions ..
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
      notran = LSAME(Trans,'N')
      IF ( .NOT.notran .AND. .NOT.LSAME(Trans,'T') .AND.                &
     &     .NOT.LSAME(Trans,'C') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Kl<0 ) THEN
         Info = -3
      ELSEIF ( Ku<0 ) THEN
         Info = -4
      ELSEIF ( Nrhs<0 ) THEN
         Info = -5
      ELSEIF ( Ldab<(2*Kl+Ku+1) ) THEN
         Info = -7
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -10
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGBTRS',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
      kd = Ku + Kl + 1
      lnoti = Kl>0
!
      IF ( notran ) THEN
!
!        Solve  A*X = B.
!
!        Solve L*X = B, overwriting B with X.
!
!        L is represented as a product of permutations and unit lower
!        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!        where each transformation L(i) is a rank-one modification of
!        the identity matrix.
!
         IF ( lnoti ) THEN
            DO j = 1 , N - 1
               lm = MIN(Kl,N-j)
               l = Ipiv(j)
               IF ( l/=j ) CALL ZSWAP(Nrhs,B(l,1),Ldb,B(j,1),Ldb)
               CALL ZGERU(lm,Nrhs,-ONE,Ab(kd+1,j),1,B(j,1),Ldb,B(j+1,1),&
     &                    Ldb)
            ENDDO
         ENDIF
!
         DO i = 1 , Nrhs
!
!           Solve U*X = B, overwriting B with X.
!
            CALL ZTBSV('Upper','No transpose','Non-unit',N,Kl+Ku,Ab,    &
     &                 Ldab,B(1,i),1)
         ENDDO
!
      ELSEIF ( LSAME(Trans,'T') ) THEN
!
!        Solve A**T * X = B.
!
         DO i = 1 , Nrhs
!
!           Solve U**T * X = B, overwriting B with X.
!
            CALL ZTBSV('Upper','Transpose','Non-unit',N,Kl+Ku,Ab,Ldab,  &
     &                 B(1,i),1)
         ENDDO
!
!        Solve L**T * X = B, overwriting B with X.
!
         IF ( lnoti ) THEN
            DO j = N - 1 , 1 , -1
               lm = MIN(Kl,N-j)
               CALL ZGEMV('Transpose',lm,Nrhs,-ONE,B(j+1,1),Ldb,        &
     &                    Ab(kd+1,j),1,ONE,B(j,1),Ldb)
               l = Ipiv(j)
               IF ( l/=j ) CALL ZSWAP(Nrhs,B(l,1),Ldb,B(j,1),Ldb)
            ENDDO
         ENDIF
!
      ELSE
!
!        Solve A**H * X = B.
!
         DO i = 1 , Nrhs
!
!           Solve U**H * X = B, overwriting B with X.
!
            CALL ZTBSV('Upper','Conjugate transpose','Non-unit',N,Kl+Ku,&
     &                 Ab,Ldab,B(1,i),1)
         ENDDO
!
!        Solve L**H * X = B, overwriting B with X.
!
         IF ( lnoti ) THEN
            DO j = N - 1 , 1 , -1
               lm = MIN(Kl,N-j)
               CALL ZLACGV(Nrhs,B(j,1),Ldb)
               CALL ZGEMV('Conjugate transpose',lm,Nrhs,-ONE,B(j+1,1),  &
     &                    Ldb,Ab(kd+1,j),1,ONE,B(j,1),Ldb)
               CALL ZLACGV(Nrhs,B(j,1),Ldb)
               l = Ipiv(j)
               IF ( l/=j ) CALL ZSWAP(Nrhs,B(l,1),Ldb,B(j,1),Ldb)
            ENDDO
         ENDIF
      ENDIF
!
!     End of ZGBTRS
!
      END SUBROUTINE ZGBTRS
