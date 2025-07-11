!*==stptri.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b STPTRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STPTRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stptri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stptri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stptri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STPTRI( UPLO, DIAG, N, AP, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       REAL               AP( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STPTRI computes the inverse of a real upper or lower triangular
!> matrix A stored in packed format.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  A is upper triangular;
!>          = 'L':  A is lower triangular.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          = 'N':  A is non-unit triangular;
!>          = 'U':  A is unit triangular.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] AP
!> \verbatim
!>          AP is REAL array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangular matrix A, stored
!>          columnwise in a linear array.  The j-th column of A is stored
!>          in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*((2*n-j)/2) = A(i,j) for j<=i<=n.
!>          See below for further details.
!>          On exit, the (triangular) inverse of the original matrix, in
!>          the same packed storage format.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular
!>                matrix is singular and its inverse can not be computed.
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
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  A triangular matrix A can be transferred to packed storage using one
!>  of the following program segments:
!>
!>  UPLO = 'U':                      UPLO = 'L':
!>
!>        JC = 1                           JC = 1
!>        DO 2 J = 1, N                    DO 2 J = 1, N
!>           DO 1 I = 1, J                    DO 1 I = J, N
!>              AP(JC+I-1) = A(I,J)              AP(JC+I-J) = A(I,J)
!>      1    CONTINUE                    1    CONTINUE
!>           JC = JC + J                      JC = JC + N - J + 1
!>      2 CONTINUE                       2 CONTINUE
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE STPTRI(Uplo,Diag,N,Ap,Info)
      IMPLICIT NONE
!*--STPTRI121
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Uplo
      INTEGER Info , N
!     ..
!     .. Array Arguments ..
      REAL Ap(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL nounit , upper
      INTEGER j , jc , jclast , jj
      REAL ajj
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL SSCAL , STPMV , XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      nounit = LSAME(Diag,'N')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( .NOT.nounit .AND. .NOT.LSAME(Diag,'U') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('STPTRI',-Info)
         RETURN
      ENDIF
!
!     Check for singularity if non-unit.
!
      IF ( nounit ) THEN
         IF ( upper ) THEN
            jj = 0
            DO Info = 1 , N
               jj = jj + Info
               IF ( Ap(jj)==ZERO ) RETURN
            ENDDO
         ELSE
            jj = 1
            DO Info = 1 , N
               IF ( Ap(jj)==ZERO ) RETURN
               jj = jj + N - Info + 1
            ENDDO
         ENDIF
         Info = 0
      ENDIF
!
      IF ( upper ) THEN
!
!        Compute inverse of upper triangular matrix.
!
         jc = 1
         DO j = 1 , N
            IF ( nounit ) THEN
               Ap(jc+j-1) = ONE/Ap(jc+j-1)
               ajj = -Ap(jc+j-1)
            ELSE
               ajj = -ONE
            ENDIF
!
!           Compute elements 1:j-1 of j-th column.
!
            CALL STPMV('Upper','No transpose',Diag,j-1,Ap,Ap(jc),1)
            CALL SSCAL(j-1,ajj,Ap(jc),1)
            jc = jc + j
         ENDDO
!
      ELSE
!
!        Compute inverse of lower triangular matrix.
!
         jc = N*(N+1)/2
         DO j = N , 1 , -1
            IF ( nounit ) THEN
               Ap(jc) = ONE/Ap(jc)
               ajj = -Ap(jc)
            ELSE
               ajj = -ONE
            ENDIF
            IF ( j<N ) THEN
!
!              Compute elements j+1:n of j-th column.
!
               CALL STPMV('Lower','No transpose',Diag,N-j,Ap(jclast),   &
     &                    Ap(jc+1),1)
               CALL SSCAL(N-j,ajj,Ap(jc+1),1)
            ENDIF
            jclast = jc
            jc = jc - N + j - 2
         ENDDO
      ENDIF
!
!
!     End of STPTRI
!
      END SUBROUTINE STPTRI
