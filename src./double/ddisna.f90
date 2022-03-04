!*==ddisna.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DDISNA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DDISNA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ddisna.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ddisna.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ddisna.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DDISNA( JOB, M, N, D, SEP, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB
!       INTEGER            INFO, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), SEP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DDISNA computes the reciprocal condition numbers for the eigenvectors
!> of a real symmetric or complex Hermitian matrix or for the left or
!> right singular vectors of a general m-by-n matrix. The reciprocal
!> condition number is the 'gap' between the corresponding eigenvalue or
!> singular value and the nearest other one.
!>
!> The bound on the error, measured by angle in radians, in the I-th
!> computed vector is given by
!>
!>        DLAMCH( 'E' ) * ( ANORM / SEP( I ) )
!>
!> where ANORM = 2-norm(A) = max( abs( D(j) ) ).  SEP(I) is not allowed
!> to be smaller than DLAMCH( 'E' )*ANORM in order to limit the size of
!> the error bound.
!>
!> DDISNA may also be used to compute error bounds for eigenvectors of
!> the generalized symmetric definite eigenproblem.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies for which problem the reciprocal condition numbers
!>          should be computed:
!>          = 'E':  the eigenvectors of a symmetric/Hermitian matrix;
!>          = 'L':  the left singular vectors of a general matrix;
!>          = 'R':  the right singular vectors of a general matrix.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          If JOB = 'L' or 'R', the number of columns of the matrix,
!>          in which case N >= 0. Ignored if JOB = 'E'.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (M) if JOB = 'E'
!>                              dimension (min(M,N)) if JOB = 'L' or 'R'
!>          The eigenvalues (if JOB = 'E') or singular values (if JOB =
!>          'L' or 'R') of the matrix, in either increasing or decreasing
!>          order. If singular values, they must be non-negative.
!> \endverbatim
!>
!> \param[out] SEP
!> \verbatim
!>          SEP is DOUBLE PRECISION array, dimension (M) if JOB = 'E'
!>                               dimension (min(M,N)) if JOB = 'L' or 'R'
!>          The reciprocal condition numbers of the vectors.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DDISNA(Job,M,N,D,Sep,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DDISNA125
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Job
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Sep
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anorm , eps , newgap , oldgap , safmin , thresh
      LOGICAL :: decr , eigen , incr , left , right , sing
      INTEGER :: i , k
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
!     .. Intrinsic Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      eigen = LSAME(Job,'E')
      left = LSAME(Job,'L')
      right = LSAME(Job,'R')
      sing = left .OR. right
      IF ( eigen ) THEN
         k = M
      ELSEIF ( sing ) THEN
         k = MIN(M,N)
      ENDIF
      IF ( .NOT.eigen .AND. .NOT.sing ) THEN
         Info = -1
      ELSEIF ( M<0 ) THEN
         Info = -2
      ELSEIF ( k<0 ) THEN
         Info = -3
      ELSE
         incr = .TRUE.
         decr = .TRUE.
         DO i = 1 , k - 1
            IF ( incr ) incr = incr .AND. D(i)<=D(i+1)
            IF ( decr ) decr = decr .AND. D(i)>=D(i+1)
         ENDDO
         IF ( sing .AND. k>0 ) THEN
            IF ( incr ) incr = incr .AND. ZERO<=D(1)
            IF ( decr ) decr = decr .AND. D(k)>=ZERO
         ENDIF
         IF ( .NOT.(incr .OR. decr) ) Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DDISNA',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( k==0 ) RETURN
!
!     Compute reciprocal condition numbers
!
      IF ( k==1 ) THEN
         Sep(1) = DLAMCH('O')
      ELSE
         oldgap = ABS(D(2)-D(1))
         Sep(1) = oldgap
         DO i = 2 , k - 1
            newgap = ABS(D(i+1)-D(i))
            Sep(i) = MIN(oldgap,newgap)
            oldgap = newgap
         ENDDO
         Sep(k) = oldgap
      ENDIF
      IF ( sing ) THEN
         IF ( (left .AND. M>N) .OR. (right .AND. M<N) ) THEN
            IF ( incr ) Sep(1) = MIN(Sep(1),D(1))
            IF ( decr ) Sep(k) = MIN(Sep(k),D(k))
         ENDIF
      ENDIF
!
!     Ensure that reciprocal condition numbers are not less than
!     threshold, in order to limit the size of the error bound
!
      eps = DLAMCH('E')
      safmin = DLAMCH('S')
      anorm = MAX(ABS(D(1)),ABS(D(k)))
      IF ( anorm==ZERO ) THEN
         thresh = eps
      ELSE
         thresh = MAX(eps*anorm,safmin)
      ENDIF
      DO i = 1 , k
         Sep(i) = MAX(Sep(i),thresh)
      ENDDO
!
!
!     End of DDISNA
!
      END SUBROUTINE DDISNA
