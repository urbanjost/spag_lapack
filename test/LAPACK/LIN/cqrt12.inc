!*==cqrt12.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cqrt12
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION cqrt12( M, N, A, LDA, S, WORK, LWORK,
!                        RWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * ), S( * )
!       COMPLEX            A( LDA, * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cqrt12 computes the singular values `svlues' of the upper trapezoid
!> of A(1:M,1:N) and returns the ratio
!>
!>      || s - svlues||/(||svlues||*eps*max(M,N))
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The M-by-N matrix A. Only the upper trapezoid is referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension (min(M,N))
!>          The singular values of the matrix A.
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
!>          The length of the array WORK. LWORK >= M*N + 2*min(M,N) +
!>          max(M,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (4*min(M,N))
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
!> \ingroup complex_lin
!
!  =====================================================================
      REAL FUNCTION cqrt12(M,N,A,Lda,S,Work,Lwork,Rwork)
      IMPLICIT NONE
!*--CQRT12100
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Lwork , M , N
!     ..
!     .. Array Arguments ..
      REAL Rwork(*) , S(*)
      COMPLEX A(Lda,*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , iscl , j , mn
      REAL anrm , bignum , nrmsvl , smlnum
!     ..
!     .. Local Arrays ..
      REAL dummy(1)
!     ..
!     .. External Functions ..
      REAL CLANGE , SASUM , SLAMCH , SNRM2
EXTERNAL    CLANGE,   SASUM,   SLAMCH,   SNRM2
!     ..
!     .. External Subroutines ..
EXTERNAL    CGEBD2,   CLASCL,   CLASET,   SAXPY,   SBDSQR,   SLABAD,   SLASCL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
      cqrt12 = ZERO
!
!     Test that enough workspace is supplied
!
      IF ( Lwork<M*N+2*MIN(M,N)+MAX(M,N) ) THEN
         CALL xerbla('cqrt12',7)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      mn = MIN(M,N)
      IF ( mn<=ZERO ) RETURN
!
      nrmsvl = SNRM2(mn,S,1)
!
!     Copy upper triangle of A into work
!
      CALL CLASET('Full',M,N,CMPLX(ZERO),CMPLX(ZERO),Work,M)
      DO j = 1 , N
         DO i = 1 , MIN(j,M)
            Work((j-1)*M+i) = A(i,j)
         ENDDO
      ENDDO
!
!     Get machine parameters
!
      smlnum = SLAMCH('S')/SLAMCH('P')
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
!
!     Scale work if max entry outside range [SMLNUM,BIGNUM]
!
      anrm = CLANGE('M',M,N,Work,M,dummy)
      iscl = 0
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL CLASCL('G',0,0,anrm,smlnum,M,N,Work,M,info)
         iscl = 1
      ELSEIF ( anrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL CLASCL('G',0,0,anrm,bignum,M,N,Work,M,info)
         iscl = 1
      ENDIF
!
      IF ( anrm/=ZERO ) THEN
!
!        Compute SVD of work
!
         CALL CGEBD2(M,N,Work,M,Rwork(1),Rwork(mn+1),Work(M*N+1),       &
     &               Work(M*N+mn+1),Work(M*N+2*mn+1),info)
         CALL SBDSQR('Upper',mn,0,0,0,Rwork(1),Rwork(mn+1),dummy,mn,    &
     &               dummy,1,dummy,mn,Rwork(2*mn+1),info)
!
         IF ( iscl==1 ) THEN
            IF ( anrm>bignum ) CALL SLASCL('G',0,0,bignum,anrm,mn,1,    &
     &           Rwork(1),mn,info)
            IF ( anrm<smlnum ) CALL SLASCL('G',0,0,smlnum,anrm,mn,1,    &
     &           Rwork(1),mn,info)
         ENDIF
!
      ELSE
!
         DO i = 1 , mn
            Rwork(i) = ZERO
         ENDDO
      ENDIF
!
!     Compare s and singular values of work
!
      CALL SAXPY(mn,-ONE,S,1,Rwork(1),1)
      cqrt12 = SASUM(mn,Rwork(1),1)/(SLAMCH('Epsilon')*REAL(MAX(M,N)))
      IF ( nrmsvl/=ZERO ) cqrt12 = cqrt12/nrmsvl
!
!
!     End of cqrt12
!
      END FUNCTION cqrt12
