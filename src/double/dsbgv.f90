!*==dsbgv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DSBGV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSBGV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbgv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbgv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbgv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z,
!                         LDZ, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, N
!       ..
!       .. Array Arguments ..
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
!> DSBGV computes all the eigenvalues, and optionally, the eigenvectors
!> of a real generalized symmetric-definite banded eigenproblem, of
!> the form A*x=(lambda)*B*x. Here A and B are assumed to be symmetric
!> and banded, and B is also positive definite.
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
!>          or the number of subdiagonals if UPLO = 'L'. KA >= 0.
!> \endverbatim
!>
!> \param[in] KB
!> \verbatim
!>          KB is INTEGER
!>          The number of superdiagonals of the matrix B if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'. KB >= 0.
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
!>          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;
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
!>          eigenvector associated with W(i). The eigenvectors are
!>          normalized so that Z**T*B*Z = I.
!>          If JOBZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          JOBZ = 'V', LDZ >= N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (3*N)
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
!> \date December 2016
!
!> \ingroup doubleOTHEReigen
!
!  =====================================================================
      SUBROUTINE DSBGV(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,W,Z,Ldz,Work,  &
     &                 Info)
      IMPLICIT NONE
!*--DSBGV181
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobz , Uplo
      INTEGER Info , Ka , Kb , Ldab , Ldbb , Ldz , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Ab(Ldab,*) , Bb(Ldbb,*) , W(*) , Work(*) ,       &
     &                 Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL upper , wantz
      CHARACTER vect
      INTEGER iinfo , inde , indwrk
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL DPBSTF , DSBGST , DSBTRD , DSTEQR , DSTERF , XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      wantz = LSAME(Jobz,'V')
      upper = LSAME(Uplo,'U')
!
      Info = 0
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
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSBGV ',-Info)
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
      CALL DSBGST(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,Z,Ldz,Work(indwrk), &
     &            iinfo)
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
!     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR.
!
      IF ( .NOT.wantz ) THEN
         CALL DSTERF(N,W,Work(inde),Info)
      ELSE
         CALL DSTEQR(Jobz,N,W,Work(inde),Z,Ldz,Work(indwrk),Info)
      ENDIF
!
!     End of DSBGV
!
      END SUBROUTINE DSBGV
