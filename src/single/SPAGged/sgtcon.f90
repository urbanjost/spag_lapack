!*==sgtcon.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SGTCON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGTCON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgtcon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgtcon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgtcon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND,
!                          WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            INFO, N
!       REAL               ANORM, RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IWORK( * )
!       REAL               D( * ), DL( * ), DU( * ), DU2( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGTCON estimates the reciprocal of the condition number of a real
!> tridiagonal matrix A using the LU factorization as computed by
!> SGTTRF.
!>
!> An estimate is obtained for norm(inv(A)), and the reciprocal of the
!> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies whether the 1-norm condition number or the
!>          infinity-norm condition number is required:
!>          = '1' or 'O':  1-norm;
!>          = 'I':         Infinity-norm.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] DL
!> \verbatim
!>          DL is REAL array, dimension (N-1)
!>          The (n-1) multipliers that define the matrix L from the
!>          LU factorization of A as computed by SGTTRF.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the upper triangular matrix U from
!>          the LU factorization of A.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is REAL array, dimension (N-1)
!>          The (n-1) elements of the first superdiagonal of U.
!> \endverbatim
!>
!> \param[in] DU2
!> \verbatim
!>          DU2 is REAL array, dimension (N-2)
!>          The (n-2) elements of the second superdiagonal of U.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices; for 1 <= i <= n, row i of the matrix was
!>          interchanged with row IPIV(i).  IPIV(i) will always be either
!>          i or i+1; IPIV(i) = i indicates a row interchange was not
!>          required.
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is REAL
!>          If NORM = '1' or 'O', the 1-norm of the original matrix A.
!>          If NORM = 'I', the infinity-norm of the original matrix A.
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The reciprocal of the condition number of the matrix A,
!>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
!>          estimate of the 1-norm of inv(A) computed in this routine.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N)
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
!> \ingroup realGTcomputational
!
!  =====================================================================
      SUBROUTINE SGTCON(Norm,N,Dl,D,Du,Du2,Ipiv,Anorm,Rcond,Work,Iwork, &
     &                  Info)
      USE S_LSAME
      USE S_SGTTRS
      USE S_SLACN2
      USE S_XERBLA
      IMPLICIT NONE
!*--SGTCON154
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Norm
      INTEGER :: N
      REAL , DIMENSION(*) :: Dl
      REAL , DIMENSION(*) :: D
      REAL , DIMENSION(*) :: Du
      REAL , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ainvnm
      INTEGER :: i , kase , kase1
      INTEGER , DIMENSION(3) :: isave
      LOGICAL :: onenrm
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
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      Info = 0
      onenrm = Norm=='1' .OR. LSAME(Norm,'O')
      IF ( .NOT.onenrm .AND. .NOT.LSAME(Norm,'I') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Anorm<ZERO ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGTCON',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      Rcond = ZERO
      IF ( N==0 ) THEN
         Rcond = ONE
         RETURN
      ELSEIF ( Anorm==ZERO ) THEN
         RETURN
      ENDIF
!
!     Check that D(1:N) is non-zero.
!
      DO i = 1 , N
         IF ( D(i)==ZERO ) RETURN
      ENDDO
!
      ainvnm = ZERO
      IF ( onenrm ) THEN
         kase1 = 1
      ELSE
         kase1 = 2
      ENDIF
      kase = 0
      DO
         CALL SLACN2(N,Work(N+1),Work,Iwork,ainvnm,kase,isave)
         IF ( kase/=0 ) THEN
            IF ( kase==kase1 ) THEN
!
!           Multiply by inv(U)*inv(L).
!
               CALL SGTTRS('No transpose',N,1,Dl,D,Du,Du2,Ipiv,Work,N,  &
     &                     Info)
            ELSE
!
!           Multiply by inv(L**T)*inv(U**T).
!
               CALL SGTTRS('Transpose',N,1,Dl,D,Du,Du2,Ipiv,Work,N,Info)
            ENDIF
            CYCLE
         ENDIF
!
!     Compute the estimate of the reciprocal condition number.
!
         IF ( ainvnm/=ZERO ) Rcond = (ONE/ainvnm)/Anorm
         EXIT
      ENDDO
!
!
!     End of SGTCON
!
      END SUBROUTINE SGTCON
