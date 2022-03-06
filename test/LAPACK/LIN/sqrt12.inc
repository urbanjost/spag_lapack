!*==sqrt12.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b sqrt12
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION sqrt12( M, N, A, LDA, S, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), S( * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sqrt12 computes the singular values `svlues' of the upper trapezoid
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
!>          A is REAL array, dimension (LDA,N)
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
!>          WORK is REAL array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK. LWORK >= max(M*N + 4*min(M,N) +
!>          max(M,N), M*N+2*MIN( M, N )+4*N).
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
!> \ingroup single_lin
!
!  =====================================================================
      REAL FUNCTION sqrt12(M,N,A,Lda,S,Work,Lwork)
      IMPLICIT NONE
!*--SQRT1293
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
      REAL A(Lda,*) , S(*) , Work(Lwork)
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
!     .. External Functions ..
      REAL SASUM , SLAMCH , SLANGE , SNRM2
EXTERNAL    SASUM,   SLAMCH,   SLANGE,   SNRM2
!     ..
!     .. External Subroutines ..
EXTERNAL    SAXPY,   SBDSQR,   SGEBD2,   SLABAD,   SLASCL,   SLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , REAL
!     ..
!     .. Local Arrays ..
      REAL dummy(1)
!     ..
!     .. Executable Statements ..
!
      sqrt12 = ZERO
!
!     Test that enough workspace is supplied
!
      IF ( Lwork<MAX(M*N+4*MIN(M,N)+MAX(M,N),M*N+2*MIN(M,N)+4*N) ) THEN
         CALL xerbla('sqrt12',7)
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
      CALL SLASET('Full',M,N,ZERO,ZERO,Work,M)
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
      anrm = SLANGE('M',M,N,Work,M,dummy)
      iscl = 0
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL SLASCL('G',0,0,anrm,smlnum,M,N,Work,M,info)
         iscl = 1
      ELSEIF ( anrm>bignum ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL SLASCL('G',0,0,anrm,bignum,M,N,Work,M,info)
         iscl = 1
      ENDIF
!
      IF ( anrm/=ZERO ) THEN
!
!        Compute SVD of work
!
         CALL SGEBD2(M,N,Work,M,Work(M*N+1),Work(M*N+mn+1),             &
     &               Work(M*N+2*mn+1),Work(M*N+3*mn+1),Work(M*N+4*mn+1),&
     &               info)
         CALL SBDSQR('Upper',mn,0,0,0,Work(M*N+1),Work(M*N+mn+1),dummy, &
     &               mn,dummy,1,dummy,mn,Work(M*N+2*mn+1),info)
!
         IF ( iscl==1 ) THEN
            IF ( anrm>bignum ) CALL SLASCL('G',0,0,bignum,anrm,mn,1,    &
     &           Work(M*N+1),mn,info)
            IF ( anrm<smlnum ) CALL SLASCL('G',0,0,smlnum,anrm,mn,1,    &
     &           Work(M*N+1),mn,info)
         ENDIF
!
      ELSE
!
         DO i = 1 , mn
            Work(M*N+i) = ZERO
         ENDDO
      ENDIF
!
!     Compare s and singular values of work
!
      CALL SAXPY(mn,-ONE,S,1,Work(M*N+1),1)
      sqrt12 = SASUM(mn,Work(M*N+1),1)                                  &
     &         /(SLAMCH('Epsilon')*REAL(MAX(M,N)))
      IF ( nrmsvl/=ZERO ) sqrt12 = sqrt12/nrmsvl
!
!
!     End of sqrt12
!
      END FUNCTION sqrt12
