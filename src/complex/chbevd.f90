!*==chbevd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> CHBEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHBEVD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbevd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbevd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbevd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,
!                          LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, KD, LDAB, LDZ, LIWORK, LRWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               RWORK( * ), W( * )
!       COMPLEX            AB( LDAB, * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHBEVD computes all the eigenvalues and, optionally, eigenvectors of
!> a complex Hermitian band matrix A.  If eigenvectors are desired, it
!> uses a divide and conquer algorithm.
!>
!> The divide and conquer algorithm makes very mild assumptions about
!> floating point arithmetic. It will work on machines with a guard
!> digit in add/subtract, or on those binary machines without guard
!> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!> Cray-2. It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.
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
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB, N)
!>          On entry, the upper or lower triangle of the Hermitian band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>
!>          On exit, AB is overwritten by values generated during the
!>          reduction to tridiagonal form.  If UPLO = 'U', the first
!>          superdiagonal and the diagonal of the tridiagonal matrix T
!>          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',
!>          the diagonal and first subdiagonal of T are returned in the
!>          first two rows of AB.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD + 1.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is REAL array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ, N)
!>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
!>          eigenvectors of the matrix A, with the i-th column of Z
!>          holding the eigenvector associated with W(i).
!>          If JOBZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          JOBZ = 'V', LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If N <= 1,               LWORK must be at least 1.
!>          If JOBZ = 'N' and N > 1, LWORK must be at least N.
!>          If JOBZ = 'V' and N > 1, LWORK must be at least 2*N**2.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal sizes of the WORK, RWORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array,
!>                                         dimension (LRWORK)
!>          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The dimension of array RWORK.
!>          If N <= 1,               LRWORK must be at least 1.
!>          If JOBZ = 'N' and N > 1, LRWORK must be at least N.
!>          If JOBZ = 'V' and N > 1, LRWORK must be at least
!>                        1 + 5*N + 2*N**2.
!>
!>          If LRWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
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
!>          The dimension of array IWORK.
!>          If JOBZ = 'N' or N <= 1, LIWORK must be at least 1.
!>          If JOBZ = 'V' and N > 1, LIWORK must be at least 3 + 5*N .
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, the algorithm failed to converge; i
!>                off-diagonal elements of an intermediate tridiagonal
!>                form did not converge to zero.
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
!> \ingroup complexOTHEReigen
!
!  =====================================================================
      SUBROUTINE CHBEVD(Jobz,Uplo,N,Kd,Ab,Ldab,W,Z,Ldz,Work,Lwork,Rwork,&
     &                  Lrwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!*--CHBEVD219
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobz , Uplo
      INTEGER Info , Kd , Ldab , Ldz , Liwork , Lrwork , Lwork , N
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      REAL Rwork(*) , W(*)
      COMPLEX Ab(Ldab,*) , Work(*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E0,0.0E0),CONE=(1.0E0,0.0E0))
!     ..
!     .. Local Scalars ..
      LOGICAL lower , lquery , wantz
      INTEGER iinfo , imax , inde , indwk2 , indwrk , iscale , liwmin , &
     &        llrwk , llwk2 , lrwmin , lwmin
      REAL anrm , bignum , eps , rmax , rmin , safmin , sigma , smlnum
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANHB , SLAMCH
      EXTERNAL LSAME , CLANHB , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMM , CHBTRD , CLACPY , CLASCL , CSTEDC , SSCAL ,      &
     &         SSTERF , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      wantz = LSAME(Jobz,'V')
      lower = LSAME(Uplo,'L')
      lquery = (Lwork==-1 .OR. Liwork==-1 .OR. Lrwork==-1)
!
      Info = 0
      IF ( N<=1 ) THEN
         lwmin = 1
         lrwmin = 1
         liwmin = 1
      ELSEIF ( wantz ) THEN
         lwmin = 2*N**2
         lrwmin = 1 + 5*N + 2*N**2
         liwmin = 3 + 5*N
      ELSE
         lwmin = N
         lrwmin = N
         liwmin = 1
      ENDIF
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(lower .OR. LSAME(Uplo,'U')) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Kd<0 ) THEN
         Info = -4
      ELSEIF ( Ldab<Kd+1 ) THEN
         Info = -6
      ELSEIF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) THEN
         Info = -9
      ENDIF
!
      IF ( Info==0 ) THEN
         Work(1) = lwmin
         Rwork(1) = lrwmin
         Iwork(1) = liwmin
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -11
         ELSEIF ( Lrwork<lrwmin .AND. .NOT.lquery ) THEN
            Info = -13
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -15
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CHBEVD',-Info)
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
         W(1) = Ab(1,1)
         IF ( wantz ) Z(1,1) = CONE
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
      anrm = CLANHB('M',Uplo,N,Kd,Ab,Ldab,Rwork)
      iscale = 0
      IF ( anrm>ZERO .AND. anrm<rmin ) THEN
         iscale = 1
         sigma = rmin/anrm
      ELSEIF ( anrm>rmax ) THEN
         iscale = 1
         sigma = rmax/anrm
      ENDIF
      IF ( iscale==1 ) THEN
         IF ( lower ) THEN
            CALL CLASCL('B',Kd,Kd,ONE,sigma,N,N,Ab,Ldab,Info)
         ELSE
            CALL CLASCL('Q',Kd,Kd,ONE,sigma,N,N,Ab,Ldab,Info)
         ENDIF
      ENDIF
!
!     Call CHBTRD to reduce Hermitian band matrix to tridiagonal form.
!
      inde = 1
      indwrk = inde + N
      indwk2 = 1 + N*N
      llwk2 = Lwork - indwk2 + 1
      llrwk = Lrwork - indwrk + 1
      CALL CHBTRD(Jobz,Uplo,N,Kd,Ab,Ldab,W,Rwork(inde),Z,Ldz,Work,iinfo)
!
!     For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEDC.
!
      IF ( .NOT.wantz ) THEN
         CALL SSTERF(N,W,Rwork(inde),Info)
      ELSE
         CALL CSTEDC('I',N,W,Rwork(inde),Work,N,Work(indwk2),llwk2,     &
     &               Rwork(indwrk),llrwk,Iwork,Liwork,Info)
         CALL CGEMM('N','N',N,N,N,CONE,Z,Ldz,Work,N,CZERO,Work(indwk2), &
     &              N)
         CALL CLACPY('A',N,N,Work(indwk2),N,Z,Ldz)
      ENDIF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      IF ( iscale==1 ) THEN
         IF ( Info==0 ) THEN
            imax = N
         ELSE
            imax = Info - 1
         ENDIF
         CALL SSCAL(imax,ONE/sigma,W,1)
      ENDIF
!
      Work(1) = lwmin
      Rwork(1) = lrwmin
      Iwork(1) = liwmin
!
!     End of CHBEVD
!
      END SUBROUTINE CHBEVD
