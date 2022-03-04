!*==cupgtr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CUPGTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CUPGTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cupgtr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cupgtr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cupgtr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDQ, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CUPGTR generates a complex unitary matrix Q which is defined as the
!> product of n-1 elementary reflectors H(i) of order n, as returned by
!> CHPTRD using packed storage:
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
!>                 call to CHPTRD;
!>          = 'L': Lower triangular packed storage used in previous
!>                 call to CHPTRD.
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
!>          AP is COMPLEX array, dimension (N*(N+1)/2)
!>          The vectors which define the elementary reflectors, as
!>          returned by CHPTRD.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (N-1)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by CHPTRD.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDQ,N)
!>          The N-by-N unitary matrix Q.
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
!>          WORK is COMPLEX array, dimension (N-1)
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
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CUPGTR(Uplo,N,Ap,Tau,Q,Ldq,Work,Info)
      USE S_CUNG2L
      USE S_CUNG2R
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CUPGTR122
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX , DIMENSION(*) :: Work
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
         CALL XERBLA('CUPGTR',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Q was determined by a call to CHPTRD with UPLO = 'U'
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
            Q(N,j) = CZERO
         ENDDO
         DO i = 1 , N - 1
            Q(i,N) = CZERO
         ENDDO
         Q(N,N) = CONE
!
!        Generate Q(1:n-1,1:n-1)
!
         CALL CUNG2L(N-1,N-1,N-1,Q,Ldq,Tau,Work,iinfo)
!
      ELSE
!
!        Q was determined by a call to CHPTRD with UPLO = 'L'.
!
!        Unpack the vectors which define the elementary reflectors and
!        set the first row and column of Q equal to those of the unit
!        matrix
!
         Q(1,1) = CONE
         DO i = 2 , N
            Q(i,1) = CZERO
         ENDDO
         ij = 3
         DO j = 2 , N
            Q(1,j) = CZERO
            DO i = j + 1 , N
               Q(i,j) = Ap(ij)
               ij = ij + 1
            ENDDO
            ij = ij + 2
         ENDDO
!
!           Generate Q(2:n,2:n)
!
         IF ( N>1 ) CALL CUNG2R(N-1,N-1,N-1,Q(2,2),Ldq,Tau,Work,iinfo)
      ENDIF
!
!     End of CUPGTR
!
      END SUBROUTINE CUPGTR
