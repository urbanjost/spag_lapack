!*==dtrcon.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DTRCON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DTRCON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrcon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrcon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrcon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK,
!                          IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORM, UPLO
!       INTEGER            INFO, LDA, N
!       DOUBLE PRECISION   RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRCON estimates the reciprocal of the condition number of a
!> triangular matrix A, in either the 1-norm or the infinity-norm.
!>
!> The norm of A is computed and an estimate is obtained for
!> norm(inv(A)), then the reciprocal of the condition number is
!> computed as
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
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
!>          upper triangular part of the array A contains the upper
!>          triangular matrix, and the strictly lower triangular part of
!>          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
!>          triangular part of the array A contains the lower triangular
!>          matrix, and the strictly upper triangular part of A is not
!>          referenced.  If DIAG = 'U', the diagonal elements of A are
!>          also not referenced and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The reciprocal of the condition number of the matrix A,
!>          computed as RCOND = 1/(norm(A) * norm(inv(A))).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (3*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N)
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
      SUBROUTINE DTRCON(Norm,Uplo,Diag,N,A,Lda,Rcond,Work,Iwork,Info)
      IMPLICIT NONE
!*--DTRCON140
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Norm , Uplo
      INTEGER Info , Lda , N
      DOUBLE PRECISION Rcond
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      DOUBLE PRECISION A(Lda,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL nounit , onenrm , upper
      CHARACTER normin
      INTEGER ix , kase , kase1
      DOUBLE PRECISION ainvnm , anorm , scale , smlnum , xnorm
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER IDAMAX
      DOUBLE PRECISION DLAMCH , DLANTR
      EXTERNAL LSAME , IDAMAX , DLAMCH , DLANTR
!     ..
!     .. External Subroutines ..
      EXTERNAL DLACN2 , DLATRS , DRSCL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      onenrm = Norm=='1' .OR. LSAME(Norm,'O')
      nounit = LSAME(Diag,'N')
!
      IF ( .NOT.onenrm .AND. .NOT.LSAME(Norm,'I') ) THEN
         Info = -1
      ELSEIF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( .NOT.nounit .AND. .NOT.LSAME(Diag,'U') ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DTRCON',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) THEN
         Rcond = ONE
         RETURN
      ENDIF
!
      Rcond = ZERO
      smlnum = DLAMCH('Safe minimum')*DBLE(MAX(1,N))
!
!     Compute the norm of the triangular matrix A.
!
      anorm = DLANTR(Norm,Uplo,Diag,N,N,A,Lda,Work)
!
!     Continue only if ANORM > 0.
!
      IF ( anorm>ZERO ) THEN
!
!        Estimate the norm of the inverse of A.
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
            CALL DLACN2(N,Work(N+1),Work,Iwork,ainvnm,kase,isave)
            IF ( kase/=0 ) THEN
               IF ( kase==kase1 ) THEN
!
!              Multiply by inv(A).
!
                  CALL DLATRS(Uplo,'No transpose',Diag,normin,N,A,Lda,  &
     &                        Work,scale,Work(2*N+1),Info)
               ELSE
!
!              Multiply by inv(A**T).
!
                  CALL DLATRS(Uplo,'Transpose',Diag,normin,N,A,Lda,Work,&
     &                        scale,Work(2*N+1),Info)
               ENDIF
               normin = 'Y'
!
!           Multiply by 1/SCALE if doing so will not cause overflow.
!
               IF ( scale/=ONE ) THEN
                  ix = IDAMAX(N,Work,1)
                  xnorm = ABS(Work(ix))
                  IF ( scale<xnorm*smlnum .OR. scale==ZERO ) EXIT
                  CALL DRSCL(N,scale,Work,1)
               ENDIF
               CYCLE
            ENDIF
!
!        Compute the estimate of the reciprocal condition number.
!
            IF ( ainvnm/=ZERO ) Rcond = (ONE/anorm)/ainvnm
            EXIT
         ENDDO
      ENDIF
!
!
!     End of DTRCON
!
      END SUBROUTINE DTRCON
