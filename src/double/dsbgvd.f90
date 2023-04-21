!*==dsbgvd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DSBGVD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSBGVD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbgvd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbgvd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbgvd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSBGVD( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W,
!                          Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   AB( LDAB, * ), BB( LDBB, * ), W( * ),
!      $                   WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSBGVD computes all the eigenvalues, and optionally, the eigenvectors
!> of a real generalized symmetric-definite banded eigenproblem, of the
!> form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric and
!> banded, and B is also positive definite.  If eigenvectors are
!> desired, it uses a divide and conquer algorithm.
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
!>          = 'U':  Upper triangles of A and B are stored;
!>          = 'L':  Lower triangles of A and B are stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in] KA
!> \verbatim
!>          KA is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.
!> \endverbatim
!>
!> \param[in] KB
!> \verbatim
!>          KB is INTEGER
!>          The number of superdiagonals of the matrix B if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KB >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB, N)
!>          On entry, the upper or lower triangle of the symmetric band
!>          matrix A, stored in the first ka+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).
!>
!>          On exit, the contents of AB are destroyed.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KA+1.
!> \endverbatim
!>
!> \param[in,out] BB
!> \verbatim
!>          BB is DOUBLE PRECISION array, dimension (LDBB, N)
!>          On entry, the upper or lower triangle of the symmetric band
!>          matrix B, stored in the first kb+1 rows of the array.  The
!>          j-th column of B is stored in the j-th column of the array BB
!>          as follows:
!>          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;
!>          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).
!>
!>          On exit, the factor S from the split Cholesky factorization
!>          B = S**T*S, as returned by DPBSTF.
!> \endverbatim
!>
!> \param[in] LDBB
!> \verbatim
!>          LDBB is INTEGER
!>          The leading dimension of the array BB.  LDBB >= KB+1.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
!>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
!>          eigenvectors, with the i-th column of Z holding the
!>          eigenvector associated with W(i).  The eigenvectors are
!>          normalized so Z**T*B*Z = I.
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
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If N <= 1,               LWORK >= 1.
!>          If JOBZ = 'N' and N > 1, LWORK >= 2*N.
!>          If JOBZ = 'V' and N > 1, LWORK >= 1 + 5*N + 2*N**2.
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
!>          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If JOBZ  = 'N' or N <= 1, LIWORK >= 1.
!>          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.
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
!>          > 0:  if INFO = i, and i is:
!>             <= N:  the algorithm failed to converge:
!>                    i off-diagonal elements of an intermediate
!>                    tridiagonal form did not converge to zero;
!>             > N:   if INFO = N + i, for 1 <= i <= N, then DPBSTF
!>                    returned INFO = i: B is not positive definite.
!>                    The factorization of B could not be completed and
!>                    no eigenvalues or eigenvectors were computed.
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
!> \date June 2016
!
!> \ingroup doubleOTHEReigen
!
!> \par Contributors:
!  ==================
!>
!>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
!
!  =====================================================================
      SUBROUTINE DSBGVD(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,W,Z,Ldz,Work, &
     &                  Lwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!*--DSBGVD231
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobz , Uplo
      INTEGER Info , Ka , Kb , Ldab , Ldbb , Ldz , Liwork , Lwork , N
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      DOUBLE PRECISION Ab(Ldab,*) , Bb(Ldbb,*) , W(*) , Work(*) ,       &
     &                 Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , upper , wantz
      CHARACTER vect
      INTEGER iinfo , inde , indwk2 , indwrk , liwmin , llwrk2 , lwmin
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DLACPY , DPBSTF , DSBGST , DSBTRD , DSTEDC ,     &
     &         DSTERF , XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      wantz = LSAME(Jobz,'V')
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1 .OR. Liwork==-1)
!
      Info = 0
      IF ( N<=1 ) THEN
         liwmin = 1
         lwmin = 1
      ELSEIF ( wantz ) THEN
         liwmin = 3 + 5*N
         lwmin = 1 + 5*N + 2*N**2
      ELSE
         liwmin = 1
         lwmin = 2*N
      ENDIF
!
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(upper .OR. LSAME(Uplo,'L')) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Ka<0 ) THEN
         Info = -4
      ELSEIF ( Kb<0 .OR. Kb>Ka ) THEN
         Info = -5
      ELSEIF ( Ldab<Ka+1 ) THEN
         Info = -7
      ELSEIF ( Ldbb<Kb+1 ) THEN
         Info = -9
      ELSEIF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) THEN
         Info = -12
      ENDIF
!
      IF ( Info==0 ) THEN
         Work(1) = lwmin
         Iwork(1) = liwmin
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -14
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -16
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSBGVD',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Form a split Cholesky factorization of B.
!
      CALL DPBSTF(Uplo,N,Kb,Bb,Ldbb,Info)
      IF ( Info/=0 ) THEN
         Info = N + Info
         RETURN
      ENDIF
!
!     Transform problem to standard eigenvalue problem.
!
      inde = 1
      indwrk = inde + N
      indwk2 = indwrk + N*N
      llwrk2 = Lwork - indwk2 + 1
      CALL DSBGST(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,Z,Ldz,Work,iinfo)
!
!     Reduce to tridiagonal form.
!
      IF ( wantz ) THEN
         vect = 'U'
      ELSE
         vect = 'N'
      ENDIF
      CALL DSBTRD(vect,Uplo,N,Ka,Ab,Ldab,W,Work(inde),Z,Ldz,Work(indwrk)&
     &            ,iinfo)
!
!     For eigenvalues only, call DSTERF. For eigenvectors, call SSTEDC.
!
      IF ( .NOT.wantz ) THEN
         CALL DSTERF(N,W,Work(inde),Info)
      ELSE
         CALL DSTEDC('I',N,W,Work(inde),Work(indwrk),N,Work(indwk2),    &
     &               llwrk2,Iwork,Liwork,Info)
         CALL DGEMM('N','N',N,N,N,ONE,Z,Ldz,Work(indwrk),N,ZERO,        &
     &              Work(indwk2),N)
         CALL DLACPY('A',N,N,Work(indwk2),N,Z,Ldz)
      ENDIF
!
      Work(1) = lwmin
      Iwork(1) = liwmin
!
!
!     End of DSBGVD
!
      END SUBROUTINE DSBGVD
