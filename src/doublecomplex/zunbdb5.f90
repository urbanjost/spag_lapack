!*==zunbdb5.f90  processed by SPAG 7.51RB at 20:09 on  3 Mar 2022
 
!> \brief \b ZUNBDB5
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNBDB5 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunbdb5.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunbdb5.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunbdb5.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNBDB5( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2,
!                           LDQ2, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2,
!      $                   N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!>\verbatim
!>
!> ZUNBDB5 orthogonalizes the column vector
!>      X = [ X1 ]
!>          [ X2 ]
!> with respect to the columns of
!>      Q = [ Q1 ] .
!>          [ Q2 ]
!> The columns of Q must be orthonormal.
!>
!> If the projection is zero according to Kahan's "twice is enough"
!> criterion, then some other vector from the orthogonal complement
!> is returned. This vector is chosen in an arbitrary but deterministic
!> way.
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
!>          X1 is COMPLEX*16 array, dimension (M1)
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
!>          X2 is COMPLEX*16 array, dimension (M2)
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
!>          Q1 is COMPLEX*16 array, dimension (LDQ1, N)
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
!>          Q2 is COMPLEX*16 array, dimension (LDQ2, N)
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
!>          WORK is COMPLEX*16 array, dimension (LWORK)
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNBDB5(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,&
     &                   Lwork,Info)
      IMPLICIT NONE
!*--ZUNBDB5161
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     July 2012
!
!     .. Scalar Arguments ..
      INTEGER Incx1 , Incx2 , Info , Ldq1 , Ldq2 , Lwork , M1 , M2 , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 Q1(Ldq1,*) , Q2(Ldq2,*) , Work(*) , X1(*) , X2(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE , ZERO
      PARAMETER (ONE=(1.0D0,0.0D0),ZERO=(0.0D0,0.0D0))
!     ..
!     .. Local Scalars ..
      INTEGER childinfo , i , j
!     ..
!     .. External Subroutines ..
      EXTERNAL ZUNBDB6 , XERBLA
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DZNRM2
      EXTERNAL DZNRM2
!     ..
!     .. Intrinsic Function ..
      INTRINSIC MAX
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
         CALL XERBLA('ZUNBDB5',-Info)
         RETURN
      ENDIF
!
!     Project X onto the orthogonal complement of Q
!
      CALL ZUNBDB6(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,Lwork,&
     &             childinfo)
!
!     If the projection is nonzero, then return
!
      IF ( DZNRM2(M1,X1,Incx1)/=ZERO .OR. DZNRM2(M2,X2,Incx2)/=ZERO )   &
     &     RETURN
!
!     Project each standard basis vector e_1,...,e_M1 in turn, stopping
!     when a nonzero projection is found
!
      DO i = 1 , M1
         DO j = 1 , M1
            X1(j) = ZERO
         ENDDO
         X1(i) = ONE
         DO j = 1 , M2
            X2(j) = ZERO
         ENDDO
         CALL ZUNBDB6(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,   &
     &                Lwork,childinfo)
         IF ( DZNRM2(M1,X1,Incx1)/=ZERO .OR. DZNRM2(M2,X2,Incx2)/=ZERO )&
     &        RETURN
      ENDDO
!
!     Project each standard basis vector e_(M1+1),...,e_(M1+M2) in turn,
!     stopping when a nonzero projection is found
!
      DO i = 1 , M2
         DO j = 1 , M1
            X1(j) = ZERO
         ENDDO
         DO j = 1 , M2
            X2(j) = ZERO
         ENDDO
         X2(i) = ONE
         CALL ZUNBDB6(M1,M2,N,X1,Incx1,X2,Incx2,Q1,Ldq1,Q2,Ldq2,Work,   &
     &                Lwork,childinfo)
         IF ( DZNRM2(M1,X1,Incx1)/=ZERO .OR. DZNRM2(M2,X2,Incx2)/=ZERO )&
     &        RETURN
      ENDDO
!
!
!     End of ZUNBDB5
!
      END SUBROUTINE ZUNBDB5
