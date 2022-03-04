!*==cunbdb6.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
 
!> \brief \b CUNBDB6
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CUNBDB6 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunbdb6.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunbdb6.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunbdb6.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2,
!                           LDQ2, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2,
!      $                   N
!       ..
!       .. Array Arguments ..
!       COMPLEX            Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!>\verbatim
!>
!> CUNBDB6 orthogonalizes the column vector
!>      X = [ X1 ]
!>          [ X2 ]
!> with respect to the columns of
!>      Q = [ Q1 ] .
!>          [ Q2 ]
!> The columns of Q must be orthonormal.
!>
!> If the projection is zero according to Kahan's "twice is enough"
!> criterion, then the zero vector is returned.
!>
!>\endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M1
!> \verbatim
!>          M1 is INTEGER
!>           The dimension of X1 and the number of rows in Q1. 0 <= M1.
!> \endverbatim
!>
!> \param[in] M2
!> \verbatim
!>          M2 is INTEGER
!>           The dimension of X2 and the number of rows in Q2. 0 <= M2.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The number of columns in Q1 and Q2. 0 <= N.
!> \endverbatim
!>
!> \param[in,out] X1
!> \verbatim
!>          X1 is COMPLEX array, dimension (M1)
!>           On entry, the top part of the vector to be orthogonalized.
!>           On exit, the top part of the projected vector.
!> \endverbatim
!>
!> \param[in] INCX1
!> \verbatim
!>          INCX1 is INTEGER
!>           Increment for entries of X1.
!> \endverbatim
!>
!> \param[in,out] X2
!> \verbatim
!>          X2 is COMPLEX array, dimension (M2)
!>           On entry, the bottom part of the vector to be
!>           orthogonalized. On exit, the bottom part of the projected
!>           vector.
!> \endverbatim
!>
!> \param[in] INCX2
!> \verbatim
!>          INCX2 is INTEGER
!>           Increment for entries of X2.
!> \endverbatim
!>
!> \param[in] Q1
!> \verbatim
!>          Q1 is COMPLEX array, dimension (LDQ1, N)
!>           The top part of the orthonormal basis matrix.
!> \endverbatim
!>
!> \param[in] LDQ1
!> \verbatim
!>          LDQ1 is INTEGER
!>           The leading dimension of Q1. LDQ1 >= M1.
!> \endverbatim
!>
!> \param[in] Q2
!> \verbatim
!>          Q2 is COMPLEX array, dimension (LDQ2, N)
!>           The bottom part of the orthonormal basis matrix.
!> \endverbatim
!>
!> \param[in] LDQ2
!> \verbatim
!>          LDQ2 is INTEGER
!>           The leading dimension of Q2. LDQ2 >= M2.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK. LWORK >= N.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           = 0:  successful exit.
!>           < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \date July 2012
!
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CUNBDB6(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,&
     &                   Lwork,Info)
      USE S_CGEMV
      USE S_CLASSQ
      USE S_XERBLA
      IMPLICIT NONE
!*--CUNBDB6162
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ALPHASQ = 0.01E0 , REALONE = 1.0E0 ,        &
     &                      REALZERO = 0.0E0
      COMPLEX , PARAMETER  ::  NEGONE = (-1.0E0,0.0E0) ,                &
     &                         ONE = (1.0E0,0.0E0) ,                    &
     &                         ZERO = (0.0E0,0.0E0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M1
      INTEGER :: M2
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: X1
      INTEGER :: Incx1
      COMPLEX , DIMENSION(*) :: X2
      INTEGER :: Incx2
      COMPLEX , DIMENSION(Ldq1,*) :: Q1
      INTEGER :: Ldq1
      COMPLEX , DIMENSION(Ldq2,*) :: Q2
      INTEGER :: Ldq2
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i
      REAL :: normsq1 , normsq2 , scl1 , scl2 , ssq1 , ssq2
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
!     .. Intrinsic Function ..
!     ..
!     .. Executable Statements ..
!
!     Test input arguments
!
      Info = 0
      IF ( M1<0 ) THEN
         Info = -1
      ELSEIF ( M2<0 ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Incx1<1 ) THEN
         Info = -5
      ELSEIF ( Incx2<1 ) THEN
         Info = -7
      ELSEIF ( Ldq1<MAX(1,M1) ) THEN
         Info = -9
      ELSEIF ( Ldq2<MAX(1,M2) ) THEN
         Info = -11
      ELSEIF ( Lwork<N ) THEN
         Info = -13
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CUNBDB6',-Info)
         RETURN
      ENDIF
!
!     First, project X onto the orthogonal complement of Q's column
!     space
!
      scl1 = REALZERO
      ssq1 = REALONE
      CALL CLASSQ(M1,X1,Incx1,scl1,ssq1)
      scl2 = REALZERO
      ssq2 = REALONE
      CALL CLASSQ(M2,X2,Incx2,scl2,ssq2)
      normsq1 = scl1**2*ssq1 + scl2**2*ssq2
!
      IF ( M1==0 ) THEN
         DO i = 1 , N
            Work(i) = ZERO
         ENDDO
      ELSE
         CALL CGEMV('C',M1,N,ONE,Q1,Ldq1,X1,Incx1,ZERO,Work,1)
      ENDIF
!
      CALL CGEMV('C',M2,N,ONE,Q2,Ldq2,X2,Incx2,ONE,Work,1)
!
      CALL CGEMV('N',M1,N,NEGONE,Q1,Ldq1,Work,1,ONE,X1,Incx1)
      CALL CGEMV('N',M2,N,NEGONE,Q2,Ldq2,Work,1,ONE,X2,Incx2)
!
      scl1 = REALZERO
      ssq1 = REALONE
      CALL CLASSQ(M1,X1,Incx1,scl1,ssq1)
      scl2 = REALZERO
      ssq2 = REALONE
      CALL CLASSQ(M2,X2,Incx2,scl2,ssq2)
      normsq2 = scl1**2*ssq1 + scl2**2*ssq2
!
!     If projection is sufficiently large in norm, then stop.
!     If projection is zero, then stop.
!     Otherwise, project again.
!
      IF ( normsq2>=ALPHASQ*normsq1 ) RETURN
!
      IF ( normsq2==ZERO ) RETURN
!
      normsq1 = normsq2
!
      DO i = 1 , N
         Work(i) = ZERO
      ENDDO
!
      IF ( M1==0 ) THEN
         DO i = 1 , N
            Work(i) = ZERO
         ENDDO
      ELSE
         CALL CGEMV('C',M1,N,ONE,Q1,Ldq1,X1,Incx1,ZERO,Work,1)
      ENDIF
!
      CALL CGEMV('C',M2,N,ONE,Q2,Ldq2,X2,Incx2,ONE,Work,1)
!
      CALL CGEMV('N',M1,N,NEGONE,Q1,Ldq1,Work,1,ONE,X1,Incx1)
      CALL CGEMV('N',M2,N,NEGONE,Q2,Ldq2,Work,1,ONE,X2,Incx2)
!
      scl1 = REALZERO
      ssq1 = REALONE
      CALL CLASSQ(M1,X1,Incx1,scl1,ssq1)
      scl2 = REALZERO
      ssq2 = REALONE
      CALL CLASSQ(M1,X1,Incx1,scl1,ssq1)
      normsq2 = scl1**2*ssq1 + scl2**2*ssq2
!
!     If second projection is sufficiently large in norm, then do
!     nothing more. Alternatively, if it shrunk significantly, then
!     truncate it to zero.
!
      IF ( normsq2<ALPHASQ*normsq1 ) THEN
         DO i = 1 , M1
            X1(i) = ZERO
         ENDDO
         DO i = 1 , M2
            X2(i) = ZERO
         ENDDO
      ENDIF
!
!
!     End of CUNBDB6
!
      END SUBROUTINE CUNBDB6
