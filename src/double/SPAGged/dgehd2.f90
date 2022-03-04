!*==dgehd2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DGEHD2 reduces a general square matrix to upper Hessenberg form using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGEHD2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgehd2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgehd2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgehd2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEHD2 reduces a real general matrix A to upper Hessenberg form H by
!> an orthogonal similarity transformation:  Q**T * A * Q = H .
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          It is assumed that A is already upper triangular in rows
!>          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!>          set by a previous call to DGEBAL; otherwise they should be
!>          set to 1 and N respectively. See Further Details.
!>          1 <= ILO <= IHI <= max(1,N).
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the n by n general matrix to be reduced.
!>          On exit, the upper triangle and the first subdiagonal of A
!>          are overwritten with the upper Hessenberg matrix H, and the
!>          elements below the first subdiagonal, with the array TAU,
!>          represent the orthogonal matrix Q as a product of elementary
!>          reflectors. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
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
!> \ingroup doubleGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of (ihi-ilo) elementary
!>  reflectors
!>
!>     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
!>  exit in A(i+2:ihi,i), and tau in TAU(i).
!>
!>  The contents of A are illustrated by the following example, with
!>  n = 7, ilo = 2 and ihi = 6:
!>
!>  on entry,                        on exit,
!>
!>  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
!>  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
!>  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
!>  (                         a )    (                          a )
!>
!>  where a denotes an element of the original matrix A, h denotes a
!>  modified element of the upper Hessenberg matrix H, and vi denotes an
!>  element of the vector defining H(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DGEHD2(N,Ilo,Ihi,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      USE S_DLARF
      USE S_DLARFG
      USE S_XERBLA
      IMPLICIT NONE
!*--DGEHD2157
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: aii
      INTEGER :: i
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
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Ilo<1 .OR. Ilo>MAX(1,N) ) THEN
         Info = -2
      ELSEIF ( Ihi<MIN(Ilo,N) .OR. Ihi>N ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGEHD2',-Info)
         RETURN
      ENDIF
!
      DO i = Ilo , Ihi - 1
!
!        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
!
         CALL DLARFG(Ihi-i,A(i+1,i),A(MIN(i+2,N),i),1,Tau(i))
         aii = A(i+1,i)
         A(i+1,i) = ONE
!
!        Apply H(i) to A(1:ihi,i+1:ihi) from the right
!
         CALL DLARF('Right',Ihi,Ihi-i,A(i+1,i),1,Tau(i),A(1,i+1),Lda,   &
     &              Work)
!
!        Apply H(i) to A(i+1:ihi,i+1:n) from the left
!
         CALL DLARF('Left',Ihi-i,N-i,A(i+1,i),1,Tau(i),A(i+1,i+1),Lda,  &
     &              Work)
!
         A(i+1,i) = aii
      ENDDO
!
!
!     End of DGEHD2
!
      END SUBROUTINE DGEHD2
