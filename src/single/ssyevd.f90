!*==ssyevd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> SSYEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYEVD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyevd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyevd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyevd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK,
!                          LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, LDA, LIWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               A( LDA, * ), W( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYEVD computes all eigenvalues and, optionally, eigenvectors of a
!> real symmetric matrix A. If eigenvectors are desired, it uses a
!> divide and conquer algorithm.
!>
!> The divide and conquer algorithm makes very mild assumptions about
!> floating point arithmetic. It will work on machines with a guard
!> digit in add/subtract, or on those binary machines without guard
!> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!> Cray-2. It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.
!>
!> Because of large use of BLAS of level 3, SSYEVD needs N**2 more
!> workspace than SSYEVX.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only;
!>          = 'V':  Compute eigenvalues and eigenvectors.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA, N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of A contains the
!>          upper triangular part of the matrix A.  If UPLO = 'L',
!>          the leading N-by-N lower triangular part of A contains
!>          the lower triangular part of the matrix A.
!>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!>          orthonormal eigenvectors of the matrix A.
!>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!>          or the upper triangle (if UPLO='U') of A, including the
!>          diagonal, is destroyed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is REAL array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array,
!>                                         dimension (LWORK)
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If N <= 1,               LWORK must be at least 1.
!>          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.
!>          If JOBZ = 'V' and N > 1, LWORK must be at least
!>                                                1 + 6*N + 2*N**2.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal sizes of the WORK and IWORK
!>          arrays, returns these values as the first entries of the WORK
!>          and IWORK arrays, and no error message related to LWORK or
!>          LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If N <= 1,                LIWORK must be at least 1.
!>          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
!>          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK and IWORK arrays, and no error message related to
!>          LWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed
!>                to converge; i off-diagonal elements of an intermediate
!>                tridiagonal form did not converge to zero;
!>                if INFO = i and JOBZ = 'V', then the algorithm failed
!>                to compute an eigenvalue while working on the submatrix
!>                lying in rows and columns INFO/(N+1) through
!>                mod(INFO,N+1).
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
!> \ingroup realSYeigen
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA \n
!>  Modified by Francoise Tisseur, University of Tennessee \n
!>  Modified description of INFO. Sven, 16 Feb 05. \n
!>
!  =====================================================================
      SUBROUTINE SSYEVD(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Iwork,Liwork,    &
     &                  Info)
      IMPLICIT NONE
!*--SSYEVD187
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobz , Uplo
      INTEGER Info , Lda , Liwork , Lwork , N
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      REAL A(Lda,*) , W(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
!
      LOGICAL lower , lquery , wantz
      INTEGER iinfo , inde , indtau , indwk2 , indwrk , iscale , liopt ,&
     &        liwmin , llwork , llwrk2 , lopt , lwmin
      REAL anrm , bignum , eps , rmax , rmin , safmin , sigma , smlnum
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      REAL SLAMCH , SLANSY
      EXTERNAL ILAENV , LSAME , SLAMCH , SLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL SLACPY , SLASCL , SORMTR , SSCAL , SSTEDC , SSTERF ,     &
     &         SSYTRD , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      wantz = LSAME(Jobz,'V')
      lower = LSAME(Uplo,'L')
      lquery = (Lwork==-1 .OR. Liwork==-1)
!
      Info = 0
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(lower .OR. LSAME(Uplo,'U')) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( N<=1 ) THEN
            liwmin = 1
            lwmin = 1
            lopt = lwmin
            liopt = liwmin
         ELSE
            IF ( wantz ) THEN
               liwmin = 3 + 5*N
               lwmin = 1 + 6*N + 2*N**2
            ELSE
               liwmin = 1
               lwmin = 2*N + 1
            ENDIF
            lopt = MAX(lwmin,2*N+ILAENV(1,'SSYTRD',Uplo,N,-1,-1,-1))
            liopt = liwmin
         ENDIF
         Work(1) = lopt
         Iwork(1) = liopt
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -8
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -10
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSYEVD',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( N==1 ) THEN
         W(1) = A(1,1)
         IF ( wantz ) A(1,1) = ONE
         RETURN
      ENDIF
!
!     Get machine constants.
!
      safmin = SLAMCH('Safe minimum')
      eps = SLAMCH('Precision')
      smlnum = safmin/eps
      bignum = ONE/smlnum
      rmin = SQRT(smlnum)
      rmax = SQRT(bignum)
!
!     Scale matrix to allowable range, if necessary.
!
      anrm = SLANSY('M',Uplo,N,A,Lda,Work)
      iscale = 0
      IF ( anrm>ZERO .AND. anrm<rmin ) THEN
         iscale = 1
         sigma = rmin/anrm
      ELSEIF ( anrm>rmax ) THEN
         iscale = 1
         sigma = rmax/anrm
      ENDIF
      IF ( iscale==1 ) CALL SLASCL(Uplo,0,0,ONE,sigma,N,N,A,Lda,Info)
!
!     Call SSYTRD to reduce symmetric matrix to tridiagonal form.
!
      inde = 1
      indtau = inde + N
      indwrk = indtau + N
      llwork = Lwork - indwrk + 1
      indwk2 = indwrk + N*N
      llwrk2 = Lwork - indwk2 + 1
!
      CALL SSYTRD(Uplo,N,A,Lda,W,Work(inde),Work(indtau),Work(indwrk),  &
     &            llwork,iinfo)
!
!     For eigenvalues only, call SSTERF.  For eigenvectors, first call
!     SSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
!     tridiagonal matrix, then call SORMTR to multiply it by the
!     Householder transformations stored in A.
!
      IF ( .NOT.wantz ) THEN
         CALL SSTERF(N,W,Work(inde),Info)
      ELSE
         CALL SSTEDC('I',N,W,Work(inde),Work(indwrk),N,Work(indwk2),    &
     &               llwrk2,Iwork,Liwork,Info)
         CALL SORMTR('L',Uplo,'N',N,N,A,Lda,Work(indtau),Work(indwrk),N,&
     &               Work(indwk2),llwrk2,iinfo)
         CALL SLACPY('A',N,N,Work(indwrk),N,A,Lda)
      ENDIF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      IF ( iscale==1 ) CALL SSCAL(N,ONE/sigma,W,1)
!
      Work(1) = lopt
      Iwork(1) = liopt
!
!
!     End of SSYEVD
!
      END SUBROUTINE SSYEVD
