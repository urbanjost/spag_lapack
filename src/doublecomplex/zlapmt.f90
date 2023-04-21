!*==zlapmt.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
 
!> \brief \b ZLAPMT performs a forward or backward permutation of the columns of a matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAPMT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlapmt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlapmt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlapmt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAPMT( FORWRD, M, N, X, LDX, K )
!
!       .. Scalar Arguments ..
!       LOGICAL            FORWRD
!       INTEGER            LDX, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            K( * )
!       COMPLEX*16         X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAPMT rearranges the columns of the M by N matrix X as specified
!> by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.
!> If FORWRD = .TRUE.,  forward permutation:
!>
!>      X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.
!>
!> If FORWRD = .FALSE., backward permutation:
!>
!>      X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] FORWRD
!> \verbatim
!>          FORWRD is LOGICAL
!>          = .TRUE., forward permutation
!>          = .FALSE., backward permutation
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix X. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix X. N >= 0.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (LDX,N)
!>          On entry, the M by N matrix X.
!>          On exit, X contains the permuted matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X, LDX >= MAX(1,M).
!> \endverbatim
!>
!> \param[in,out] K
!> \verbatim
!>          K is INTEGER array, dimension (N)
!>          On entry, K contains the permutation vector. K is used as
!>          internal workspace, but reset to its original value on
!>          output.
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLAPMT(Forwrd,M,N,X,Ldx,K)
      IMPLICIT NONE
!*--ZLAPMT109
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Forwrd
      INTEGER Ldx , M , N
!     ..
!     .. Array Arguments ..
      INTEGER K(*)
      COMPLEX*16 X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , ii , in , j
      COMPLEX*16 temp
!     ..
!     .. Executable Statements ..
!
      IF ( N<=1 ) RETURN
!
      DO i = 1 , N
         K(i) = -K(i)
      ENDDO
!
      IF ( Forwrd ) THEN
!
!        Forward permutation
!
         DO i = 1 , N
!
            IF ( K(i)<=0 ) THEN
!
               j = i
               K(j) = -K(j)
               in = K(j)
!
               DO WHILE ( K(in)<=0 )
!
                  DO ii = 1 , M
                     temp = X(ii,j)
                     X(ii,j) = X(ii,in)
                     X(ii,in) = temp
                  ENDDO
!
                  K(in) = -K(in)
                  j = in
                  in = K(in)
               ENDDO
            ENDIF
!
!
         ENDDO
!
      ELSE
!
!        Backward permutation
!
         DO i = 1 , N
!
            IF ( K(i)<=0 ) THEN
!
               K(i) = -K(i)
               j = K(i)
               DO WHILE ( j/=i )
!
                  DO ii = 1 , M
                     temp = X(ii,i)
                     X(ii,i) = X(ii,j)
                     X(ii,j) = temp
                  ENDDO
!
                  K(j) = -K(j)
                  j = K(j)
               ENDDO
            ENDIF
!
!
         ENDDO
!
      ENDIF
!
!
!     End of ZLAPMT
!
      END SUBROUTINE ZLAPMT
