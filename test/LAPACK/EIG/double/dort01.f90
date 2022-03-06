!*==dort01.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DORT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORT01( ROWCOL, M, N, U, LDU, WORK, LWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          ROWCOL
!       INTEGER            LDU, LWORK, M, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   U( LDU, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORT01 checks that the matrix U is orthogonal by computing the ratio
!>
!>    RESID = norm( I - U*U' ) / ( n * EPS ), if ROWCOL = 'R',
!> or
!>    RESID = norm( I - U'*U ) / ( m * EPS ), if ROWCOL = 'C'.
!>
!> Alternatively, if there isn't sufficient workspace to form
!> I - U*U' or I - U'*U, the ratio is computed as
!>
!>    RESID = abs( I - U*U' ) / ( n * EPS ), if ROWCOL = 'R',
!> or
!>    RESID = abs( I - U'*U ) / ( m * EPS ), if ROWCOL = 'C'.
!>
!> where EPS is the machine precision.  ROWCOL is used only if m = n;
!> if m > n, ROWCOL is assumed to be 'C', and if m < n, ROWCOL is
!> assumed to be 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ROWCOL
!> \verbatim
!>          ROWCOL is CHARACTER
!>          Specifies whether the rows or columns of U should be checked
!>          for orthogonality.  Used only if M = N.
!>          = 'R':  Check for orthogonal rows of U
!>          = 'C':  Check for orthogonal columns of U
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix U.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix U.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU,N)
!>          The orthogonal matrix U.  U is checked for orthogonal columns
!>          if m > n or if m = n and ROWCOL = 'C'.  U is checked for
!>          orthogonal rows if m < n or if m = n and ROWCOL = 'R'.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  For best performance, LWORK
!>          should be at least N*(N+1) if ROWCOL = 'C' or M*(M+1) if
!>          ROWCOL = 'R', but the test will be done even if LWORK is 0.
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          RESID = norm( I - U * U' ) / ( n * EPS ), if ROWCOL = 'R', or
!>          RESID = norm( I - U' * U ) / ( m * EPS ), if ROWCOL = 'C'.
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DORT01(Rowcol,M,N,U,Ldu,Work,Lwork,Resid)
      IMPLICIT NONE
!*--DORT01120
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Rowcol
      INTEGER Ldu , Lwork , M , N
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION U(Ldu,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      CHARACTER transu
      INTEGER i , j , k , ldwork , mnmin
      DOUBLE PRECISION eps , tmp
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DDOT , DLAMCH , DLANSY
      EXTERNAL LSAME , DDOT , DLAMCH , DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL DLASET , DSYRK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN
!     ..
!     .. Executable Statements ..
!
      Resid = ZERO
!
!     Quick return if possible
!
      IF ( M<=0 .OR. N<=0 ) RETURN
!
      eps = DLAMCH('Precision')
      IF ( M<N .OR. (M==N .AND. LSAME(Rowcol,'R')) ) THEN
         transu = 'N'
         k = N
      ELSE
         transu = 'T'
         k = M
      ENDIF
      mnmin = MIN(M,N)
!
      IF ( (mnmin+1)*mnmin<=Lwork ) THEN
         ldwork = mnmin
      ELSE
         ldwork = 0
      ENDIF
      IF ( ldwork>0 ) THEN
!
!        Compute I - U*U' or I - U'*U.
!
         CALL DLASET('Upper',mnmin,mnmin,ZERO,ONE,Work,ldwork)
         CALL DSYRK('Upper',transu,mnmin,k,-ONE,U,Ldu,ONE,Work,ldwork)
!
!        Compute norm( I - U*U' ) / ( K * EPS ) .
!
         Resid = DLANSY('1','Upper',mnmin,Work,ldwork,                  &
     &           Work(ldwork*mnmin+1))
         Resid = (Resid/DBLE(k))/eps
      ELSEIF ( transu=='T' ) THEN
!
!        Find the maximum element in abs( I - U'*U ) / ( m * EPS )
!
         DO j = 1 , N
            DO i = 1 , j
               IF ( i/=j ) THEN
                  tmp = ZERO
               ELSE
                  tmp = ONE
               ENDIF
               tmp = tmp - DDOT(M,U(1,i),1,U(1,j),1)
               Resid = MAX(Resid,ABS(tmp))
            ENDDO
         ENDDO
         Resid = (Resid/DBLE(M))/eps
      ELSE
!
!        Find the maximum element in abs( I - U*U' ) / ( n * EPS )
!
         DO j = 1 , M
            DO i = 1 , j
               IF ( i/=j ) THEN
                  tmp = ZERO
               ELSE
                  tmp = ONE
               ENDIF
               tmp = tmp - DDOT(N,U(j,1),Ldu,U(i,1),Ldu)
               Resid = MAX(Resid,ABS(tmp))
            ENDDO
         ENDDO
         Resid = (Resid/DBLE(N))/eps
      ENDIF
!
!     End of DORT01
!
      END SUBROUTINE DORT01
