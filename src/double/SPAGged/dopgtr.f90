!*==dopgtr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DOPGTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DOPGTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dopgtr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dopgtr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dopgtr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDQ, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DOPGTR generates a real orthogonal matrix Q which is defined as the
!> product of n-1 elementary reflectors H(i) of order n, as returned by
!> DSPTRD using packed storage:
!>
!> if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!>
!> if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U': Upper triangular packed storage used in previous
!>                 call to DSPTRD;
!>          = 'L': Lower triangular packed storage used in previous
!>                 call to DSPTRD.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix Q. N >= 0.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          The vectors which define the elementary reflectors, as
!>          returned by DSPTRD.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (N-1)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DSPTRD.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>          The N-by-N orthogonal matrix Q.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q. LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N-1)
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
      SUBROUTINE DOPGTR(Uplo,N,Ap,Tau,Q,Ldq,Work,Info)
      USE F77KINDS                        
      USE S_DORG2L
      USE S_DORG2R
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DOPGTR123
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Ap
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , iinfo , ij , j
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Ldq<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DOPGTR',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Q was determined by a call to DSPTRD with UPLO = 'U'
!
!        Unpack the vectors which define the elementary reflectors and
!        set the last row and column of Q equal to those of the unit
!        matrix
!
         ij = 2
         DO j = 1 , N - 1
            DO i = 1 , j - 1
               Q(i,j) = Ap(ij)
               ij = ij + 1
            ENDDO
            ij = ij + 2
            Q(N,j) = ZERO
         ENDDO
         DO i = 1 , N - 1
            Q(i,N) = ZERO
         ENDDO
         Q(N,N) = ONE
!
!        Generate Q(1:n-1,1:n-1)
!
         CALL DORG2L(N-1,N-1,N-1,Q,Ldq,Tau,Work,iinfo)
!
      ELSE
!
!        Q was determined by a call to DSPTRD with UPLO = 'L'.
!
!        Unpack the vectors which define the elementary reflectors and
!        set the first row and column of Q equal to those of the unit
!        matrix
!
         Q(1,1) = ONE
         DO i = 2 , N
            Q(i,1) = ZERO
         ENDDO
         ij = 3
         DO j = 2 , N
            Q(1,j) = ZERO
            DO i = j + 1 , N
               Q(i,j) = Ap(ij)
               ij = ij + 1
            ENDDO
            ij = ij + 2
         ENDDO
!
!           Generate Q(2:n,2:n)
!
         IF ( N>1 ) CALL DORG2R(N-1,N-1,N-1,Q(2,2),Ldq,Tau,Work,iinfo)
      ENDIF
!
!     End of DOPGTR
!
      END SUBROUTINE DOPGTR
