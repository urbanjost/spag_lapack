!*==zspt03.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zspt03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zspt03( UPLO, N, A, AINV, WORK, LDW, RWORK, RCOND,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDW, N
!       DOUBLE PRECISION   RCOND, RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( * ), AINV( * ), WORK( LDW, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zspt03 computes the residual for a complex symmetric packed matrix
!> times its inverse:
!>    norm( I - A*AINV ) / ( N * norm(A) * norm(AINV) * EPS ),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          complex symmetric matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The original complex symmetric matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in] AINV
!> \verbatim
!>          AINV is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The (symmetric) inverse of the matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LDW,N)
!> \endverbatim
!>
!> \param[in] LDW
!> \verbatim
!>          LDW is INTEGER
!>          The leading dimension of the array WORK.  LDW >= max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The reciprocal of the condition number of A, computed as
!>          ( 1/norm(A) ) / norm(AINV).
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          norm(I - A*AINV) / ( N * norm(A) * norm(AINV) * EPS )
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zspt03(Uplo,N,A,Ainv,Work,Ldw,Rwork,Rcond,Resid)
      IMPLICIT NONE
!*--ZSPT03113
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Ldw , N
      DOUBLE PRECISION Rcond , Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(*) , Ainv(*) , Work(Ldw,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , icol , j , jcol , k , kcol , nall
      DOUBLE PRECISION ainvnm , anorm , eps
      COMPLEX*16 t
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , ZLANGE , ZLANSP
      COMPLEX*16 ZDOTU
EXTERNAL    LSAME,   DLAMCH,   ZLANGE,   ZLANSP,   ZDOTU
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0.
!
      IF ( N<=0 ) THEN
         Rcond = ONE
         Resid = ZERO
         RETURN
      ENDIF
!
!     Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
!
      eps = DLAMCH('Epsilon')
      anorm = ZLANSP('1',Uplo,N,A,Rwork)
      ainvnm = ZLANSP('1',Uplo,N,Ainv,Rwork)
      IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
         Rcond = ZERO
         Resid = ONE/eps
         RETURN
      ENDIF
      Rcond = (ONE/anorm)/ainvnm
!
!     Case where both A and AINV are upper triangular:
!     Each element of - A * AINV is computed by taking the dot product
!     of a row of A with a column of AINV.
!
      IF ( LSAME(Uplo,'U') ) THEN
         DO i = 1 , N
            icol = ((i-1)*i)/2 + 1
!
!           Code when J <= I
!
            DO j = 1 , i
               jcol = ((j-1)*j)/2 + 1
               t = ZDOTU(j,A(icol),1,Ainv(jcol),1)
               jcol = jcol + 2*j - 1
               kcol = icol - 1
               DO k = j + 1 , i
                  t = t + A(kcol+k)*Ainv(jcol)
                  jcol = jcol + k
               ENDDO
               kcol = kcol + 2*i
               DO k = i + 1 , N
                  t = t + A(kcol)*Ainv(jcol)
                  kcol = kcol + k
                  jcol = jcol + k
               ENDDO
               Work(i,j) = -t
            ENDDO
!
!           Code when J > I
!
            DO j = i + 1 , N
               jcol = ((j-1)*j)/2 + 1
               t = ZDOTU(i,A(icol),1,Ainv(jcol),1)
               jcol = jcol - 1
               kcol = icol + 2*i - 1
               DO k = i + 1 , j
                  t = t + A(kcol)*Ainv(jcol+k)
                  kcol = kcol + k
               ENDDO
               jcol = jcol + 2*j
               DO k = j + 1 , N
                  t = t + A(kcol)*Ainv(jcol)
                  kcol = kcol + k
                  jcol = jcol + k
               ENDDO
               Work(i,j) = -t
            ENDDO
         ENDDO
      ELSE
!
!        Case where both A and AINV are lower triangular
!
         nall = (N*(N+1))/2
         DO i = 1 , N
!
!           Code when J <= I
!
            icol = nall - ((N-i+1)*(N-i+2))/2 + 1
            DO j = 1 , i
               jcol = nall - ((N-j)*(N-j+1))/2 - (N-i)
               t = ZDOTU(N-i+1,A(icol),1,Ainv(jcol),1)
               kcol = i
               jcol = j
               DO k = 1 , j - 1
                  t = t + A(kcol)*Ainv(jcol)
                  jcol = jcol + N - k
                  kcol = kcol + N - k
               ENDDO
               jcol = jcol - j
               DO k = j , i - 1
                  t = t + A(kcol)*Ainv(jcol+k)
                  kcol = kcol + N - k
               ENDDO
               Work(i,j) = -t
            ENDDO
!
!           Code when J > I
!
            icol = nall - ((N-i)*(N-i+1))/2
            DO j = i + 1 , N
               jcol = nall - ((N-j+1)*(N-j+2))/2 + 1
               t = ZDOTU(N-j+1,A(icol-N+j),1,Ainv(jcol),1)
               kcol = i
               jcol = j
               DO k = 1 , i - 1
                  t = t + A(kcol)*Ainv(jcol)
                  jcol = jcol + N - k
                  kcol = kcol + N - k
               ENDDO
               kcol = kcol - i
               DO k = i , j - 1
                  t = t + A(kcol+k)*Ainv(jcol)
                  jcol = jcol + N - k
               ENDDO
               Work(i,j) = -t
            ENDDO
         ENDDO
      ENDIF
!
!     Add the identity matrix to WORK .
!
      DO i = 1 , N
         Work(i,i) = Work(i,i) + ONE
      ENDDO
!
!     Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)
!
      Resid = ZLANGE('1',N,N,Work,Ldw,Rwork)
!
      Resid = ((Resid*Rcond)/eps)/DBLE(N)
!
!
!     End of zspt03
!
      END SUBROUTINE zspt03
