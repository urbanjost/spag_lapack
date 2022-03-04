!*==zpbtrs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZPBTRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZPBTRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpbtrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpbtrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpbtrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KD, LDAB, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         AB( LDAB, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPBTRS solves a system of linear equations A*X = B with a Hermitian
!> positive definite band matrix A using the Cholesky factorization
!> A = U**H *U or A = L*L**H computed by ZPBTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangular factor stored in AB;
!>          = 'L':  Lower triangular factor stored in AB.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
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
!>          The triangular factor U or L from the Cholesky factorization
!>          A = U**H *U or A = L*L**H of the band matrix A, stored in the
!>          first KD+1 rows of the array.  The j-th column of U or L is
!>          stored in the j-th column of the array AB as follows:
!>          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZPBTRS(Uplo,N,Kd,Nrhs,Ab,Ldab,B,Ldb,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZTBSV
      IMPLICIT NONE
!*--ZPBTRS129
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: j
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
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Kd<0 ) THEN
         Info = -3
      ELSEIF ( Nrhs<0 ) THEN
         Info = -4
      ELSEIF ( Ldab<Kd+1 ) THEN
         Info = -6
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZPBTRS',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Solve A*X = B where A = U**H *U.
!
         DO j = 1 , Nrhs
!
!           Solve U**H *X = B, overwriting B with X.
!
            CALL ZTBSV('Upper','Conjugate transpose','Non-unit',N,Kd,Ab,&
     &                 Ldab,B(1,j),1)
!
!           Solve U*X = B, overwriting B with X.
!
            CALL ZTBSV('Upper','No transpose','Non-unit',N,Kd,Ab,Ldab,  &
     &                 B(1,j),1)
         ENDDO
      ELSE
!
!        Solve A*X = B where A = L*L**H.
!
         DO j = 1 , Nrhs
!
!           Solve L*X = B, overwriting B with X.
!
            CALL ZTBSV('Lower','No transpose','Non-unit',N,Kd,Ab,Ldab,  &
     &                 B(1,j),1)
!
!           Solve L**H *X = B, overwriting B with X.
!
            CALL ZTBSV('Lower','Conjugate transpose','Non-unit',N,Kd,Ab,&
     &                 Ldab,B(1,j),1)
         ENDDO
      ENDIF
!
!
!     End of ZPBTRS
!
      END SUBROUTINE ZPBTRS
