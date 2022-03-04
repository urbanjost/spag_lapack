!*==cgecon.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CGECON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGECON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgecon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgecon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgecon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            INFO, LDA, N
!       REAL               ANORM, RCOND
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGECON estimates the reciprocal of the condition number of a general
!> complex matrix A, in either the 1-norm or the infinity-norm, using
!> the LU factorization computed by CGETRF.
!>
!> An estimate is obtained for norm(inv(A)), and the reciprocal of the
!> condition number is computed as
!>    RCOND = 1 / ( norm(A) * norm(inv(A)) ).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies whether the 1-norm condition number or the
!>          infinity-norm condition number is required:
!>          = '1' or 'O':  1-norm;
!>          = 'I':         Infinity-norm.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The factors L and U from the factorization A = P*L*U
!>          as computed by CGETRF.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is REAL
!>          If NORM = '1' or 'O', the 1-norm of the original matrix A.
!>          If NORM = 'I', the infinity-norm of the original matrix A.
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The reciprocal of the condition number of the matrix A,
!>          computed as RCOND = 1/(norm(A) * norm(inv(A))).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (2*N)
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
!> \ingroup complexGEcomputational
!
!  =====================================================================
      SUBROUTINE CGECON(Norm,N,A,Lda,Anorm,Rcond,Work,Rwork,Info)
      IMPLICIT NONE
!*--CGECON127
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Norm
      INTEGER Info , Lda , N
      REAL Anorm , Rcond
!     ..
!     .. Array Arguments ..
      REAL Rwork(*)
      COMPLEX A(Lda,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL onenrm
      CHARACTER normin
      INTEGER ix , kase , kase1
      REAL ainvnm , scale , sl , smlnum , su
      COMPLEX zdum
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ICAMAX
      REAL SLAMCH
      EXTERNAL LSAME , ICAMAX , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CLACN2 , CLATRS , CSRSCL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , MAX , REAL
!     ..
!     .. Statement Functions ..
      REAL CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      onenrm = Norm=='1' .OR. LSAME(Norm,'O')
      IF ( .NOT.onenrm .AND. .NOT.LSAME(Norm,'I') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Anorm<ZERO ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGECON',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      Rcond = ZERO
      IF ( N==0 ) THEN
         Rcond = ONE
         RETURN
      ELSEIF ( Anorm==ZERO ) THEN
         RETURN
      ENDIF
!
      smlnum = SLAMCH('Safe minimum')
!
!     Estimate the norm of inv(A).
!
      ainvnm = ZERO
      normin = 'N'
      IF ( onenrm ) THEN
         kase1 = 1
      ELSE
         kase1 = 2
      ENDIF
      kase = 0
      DO
         CALL CLACN2(N,Work(N+1),Work,ainvnm,kase,isave)
         IF ( kase/=0 ) THEN
            IF ( kase==kase1 ) THEN
!
!           Multiply by inv(L).
!
               CALL CLATRS('Lower','No transpose','Unit',normin,N,A,Lda,&
     &                     Work,sl,Rwork,Info)
!
!           Multiply by inv(U).
!
               CALL CLATRS('Upper','No transpose','Non-unit',normin,N,A,&
     &                     Lda,Work,su,Rwork(N+1),Info)
            ELSE
!
!           Multiply by inv(U**H).
!
               CALL CLATRS('Upper','Conjugate transpose','Non-unit',    &
     &                     normin,N,A,Lda,Work,su,Rwork(N+1),Info)
!
!           Multiply by inv(L**H).
!
               CALL CLATRS('Lower','Conjugate transpose','Unit',normin, &
     &                     N,A,Lda,Work,sl,Rwork,Info)
            ENDIF
!
!        Divide X by 1/(SL*SU) if doing so will not cause overflow.
!
            scale = sl*su
            normin = 'Y'
            IF ( scale/=ONE ) THEN
               ix = ICAMAX(N,Work,1)
               IF ( scale<CABS1(Work(ix))*smlnum .OR. scale==ZERO ) EXIT
               CALL CSRSCL(N,scale,Work,1)
            ENDIF
            CYCLE
         ENDIF
!
!     Compute the estimate of the reciprocal condition number.
!
         IF ( ainvnm/=ZERO ) Rcond = (ONE/ainvnm)/Anorm
         EXIT
      ENDDO
!
!
!     End of CGECON
!
      END SUBROUTINE CGECON
