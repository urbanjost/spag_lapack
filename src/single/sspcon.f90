!*==sspcon.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SSPCON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSPCON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspcon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspcon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspcon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, IWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, N
!       REAL               ANORM, RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IWORK( * )
!       REAL               AP( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSPCON estimates the reciprocal of the condition number (in the
!> 1-norm) of a real symmetric packed matrix A using the factorization
!> A = U*D*U**T or A = L*D*L**T computed by SSPTRF.
!>
!> An estimate is obtained for norm(inv(A)), and the reciprocal of the
!> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are stored
!>          as an upper or lower triangular matrix.
!>          = 'U':  Upper triangular, form is A = U*D*U**T;
!>          = 'L':  Lower triangular, form is A = L*D*L**T.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is REAL array, dimension (N*(N+1)/2)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by SSPTRF, stored as a
!>          packed triangular matrix.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by SSPTRF.
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is REAL
!>          The 1-norm of the original matrix A.
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
!> \ingroup realOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SSPCON(Uplo,N,Ap,Ipiv,Anorm,Rcond,Work,Iwork,Info)
      IMPLICIT NONE
!*--SSPCON128
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , N
      REAL Anorm , Rcond
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*) , Iwork(*)
      REAL Ap(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER i , ip , kase
      REAL ainvnm
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL SLACN2 , SSPTRS , XERBLA
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
      ELSEIF ( Anorm<ZERO ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSPCON',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      Rcond = ZERO
      IF ( N==0 ) THEN
         Rcond = ONE
         RETURN
      ELSEIF ( Anorm<=ZERO ) THEN
         RETURN
      ENDIF
!
!     Check that the diagonal matrix D is nonsingular.
!
      IF ( upper ) THEN
!
!        Upper triangular storage: examine D from bottom to top
!
         ip = N*(N+1)/2
         DO i = N , 1 , -1
            IF ( Ipiv(i)>0 .AND. Ap(ip)==ZERO ) RETURN
            ip = ip - i
         ENDDO
      ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
         ip = 1
         DO i = 1 , N
            IF ( Ipiv(i)>0 .AND. Ap(ip)==ZERO ) RETURN
            ip = ip + N - i + 1
         ENDDO
      ENDIF
!
!     Estimate the 1-norm of the inverse.
!
      kase = 0
      DO
         CALL SLACN2(N,Work(N+1),Work,Iwork,ainvnm,kase,isave)
         IF ( kase/=0 ) THEN
!
!        Multiply by inv(L*D*L**T) or inv(U*D*U**T).
!
            CALL SSPTRS(Uplo,N,1,Ap,Ipiv,Work,N,Info)
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
!     End of SSPCON
!
      END SUBROUTINE SSPCON
