!*==slasq1.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLASQ1 computes the singular values of a real square bidiagonal matrix. Used by sbdsqr.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASQ1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASQ1( N, D, E, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLASQ1 computes the singular values of a real N-by-N bidiagonal
!> matrix with diagonal D and off-diagonal E. The singular values
!> are computed to high relative accuracy, in the absence of
!> denormalization, underflow and overflow. The algorithm was first
!> presented in
!>
!> "Accurate singular values and differential qd algorithms" by K. V.
!> Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
!> 1994,
!>
!> and the present implementation is described in "An implementation of
!> the dqds Algorithm (Positive Case)", LAPACK Working Note.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>        The number of rows and columns in the matrix. N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>        On entry, D contains the diagonal elements of the
!>        bidiagonal matrix whose SVD is desired. On normal exit,
!>        D contains the singular values in decreasing order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N)
!>        On entry, elements E(1:N-1) contain the off-diagonal elements
!>        of the bidiagonal matrix whose SVD is desired.
!>        On exit, E is overwritten.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (4*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>        = 0: successful exit
!>        < 0: if INFO = -i, the i-th argument had an illegal value
!>        > 0: the algorithm failed
!>             = 1, a split was marked by a positive value in E
!>             = 2, current block of Z not diagonalized after 100*N
!>                  iterations (in inner while loop)  On exit D and E
!>                  represent a matrix with the same singular values
!>                  which the calling subroutine could use to finish the
!>                  computation, or even feed back into SLASQ1
!>             = 3, termination criterion of outer while loop not met
!>                  (program created more than N unreduced blocks)
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
      SUBROUTINE SLASQ1(N,D,E,Work,Info)
      USE S_SCOPY
      USE S_SLAMCH
      USE S_SLAS2
      USE S_SLASCL
      USE S_SLASQ2
      USE S_SLASRT
      USE S_XERBLA
      IMPLICIT NONE
!*--SLASQ1119
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(INOUT) , DIMENSION(*) :: E
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: eps , safmin , scale , sigmn , sigmx
      INTEGER :: i , iinfo
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
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
         CALL XERBLA('SLASQ1',-Info)
         RETURN
      ELSEIF ( N==0 ) THEN
         RETURN
      ELSEIF ( N==1 ) THEN
         D(1) = ABS(D(1))
         RETURN
      ELSEIF ( N==2 ) THEN
         CALL SLAS2(D(1),E(1),D(2),sigmn,sigmx)
         D(1) = sigmx
         D(2) = sigmn
         RETURN
      ENDIF
!
!     Estimate the largest singular value.
!
      sigmx = ZERO
      DO i = 1 , N - 1
         D(i) = ABS(D(i))
         sigmx = MAX(sigmx,ABS(E(i)))
      ENDDO
      D(N) = ABS(D(N))
!
!     Early return if SIGMX is zero (matrix is already diagonal).
!
      IF ( sigmx==ZERO ) THEN
         CALL SLASRT('D',N,D,iinfo)
         RETURN
      ENDIF
!
      DO i = 1 , N
         sigmx = MAX(sigmx,D(i))
      ENDDO
!
!     Copy D and E into WORK (in the Z format) and scale (squaring the
!     input data makes scaling by a power of the radix pointless).
!
      eps = SLAMCH('Precision')
      safmin = SLAMCH('Safe minimum')
      scale = SQRT(eps/safmin)
      CALL SCOPY(N,D,1,Work(1),2)
      CALL SCOPY(N-1,E,1,Work(2),2)
      CALL SLASCL('G',0,0,sigmx,scale,2*N-1,1,Work,2*N-1,iinfo)
!
!     Compute the q's and e's.
!
      DO i = 1 , 2*N - 1
         Work(i) = Work(i)**2
      ENDDO
      Work(2*N) = ZERO
!
      CALL SLASQ2(N,Work,Info)
!
      IF ( Info==0 ) THEN
         DO i = 1 , N
            D(i) = SQRT(Work(i))
         ENDDO
         CALL SLASCL('G',0,0,scale,sigmx,N,1,D,N,iinfo)
      ELSEIF ( Info==2 ) THEN
!
!     Maximum number of iterations exceeded.  Move data from WORK
!     into D and E so the calling subroutine can try to finish
!
         DO i = 1 , N
            D(i) = SQRT(Work(2*i-1))
            E(i) = SQRT(Work(2*i))
         ENDDO
         CALL SLASCL('G',0,0,scale,sigmx,N,1,D,N,iinfo)
         CALL SLASCL('G',0,0,scale,sigmx,N,1,E,N,iinfo)
      ENDIF
!
!
!     End of SLASQ1
!
      END SUBROUTINE SLASQ1
