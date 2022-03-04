!*==dorghr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DORGHR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORGHR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorghr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorghr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorghr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
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
!> DORGHR generates a real orthogonal matrix Q which is defined as the
!> product of IHI-ILO elementary reflectors of order N, as returned by
!> DGEHRD:
!>
!> Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix Q. N >= 0.
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
!>          ILO and IHI must have the same values as in the previous call
!>          of DGEHRD. Q is equal to the unit matrix except in the
!>          submatrix Q(ilo+1:ihi,ilo+1:ihi).
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the vectors which define the elementary reflectors,
!>          as returned by DGEHRD.
!>          On exit, the N-by-N orthogonal matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (N-1)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEHRD.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= IHI-ILO.
!>          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
!>          the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
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
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DORGHR(N,Ilo,Ihi,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_DORGQR
      USE S_ILAENV
      USE S_XERBLA
      IMPLICIT NONE
!*--DORGHR134
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , iinfo , j , lwkopt , nb , nh
      LOGICAL :: lquery
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
!     Test the input arguments
!
      Info = 0
      nh = Ihi - Ilo
      lquery = (Lwork==-1)
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Ilo<1 .OR. Ilo>MAX(1,N) ) THEN
         Info = -2
      ELSEIF ( Ihi<MIN(Ilo,N) .OR. Ihi>N ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Lwork<MAX(1,nh) .AND. .NOT.lquery ) THEN
         Info = -8
      ENDIF
!
      IF ( Info==0 ) THEN
         nb = ILAENV(1,'DORGQR',' ',nh,nh,nh,-1)
         lwkopt = MAX(1,nh)*nb
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DORGHR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
!     Shift the vectors which define the elementary reflectors one
!     column to the right, and set the first ilo and the last n-ihi
!     rows and columns to those of the unit matrix
!
      DO j = Ihi , Ilo + 1 , -1
         DO i = 1 , j - 1
            A(i,j) = ZERO
         ENDDO
         DO i = j + 1 , Ihi
            A(i,j) = A(i,j-1)
         ENDDO
         DO i = Ihi + 1 , N
            A(i,j) = ZERO
         ENDDO
      ENDDO
      DO j = 1 , Ilo
         DO i = 1 , N
            A(i,j) = ZERO
         ENDDO
         A(j,j) = ONE
      ENDDO
      DO j = Ihi + 1 , N
         DO i = 1 , N
            A(i,j) = ZERO
         ENDDO
         A(j,j) = ONE
      ENDDO
!
!
!        Generate Q(ilo+1:ihi,ilo+1:ihi)
!
      IF ( nh>0 ) CALL DORGQR(nh,nh,nh,A(Ilo+1,Ilo+1),Lda,Tau(Ilo),Work,&
     &                        Lwork,iinfo)
      Work(1) = lwkopt
!
!     End of DORGHR
!
      END SUBROUTINE DORGHR
