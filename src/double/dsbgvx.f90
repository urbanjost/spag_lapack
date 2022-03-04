!*==dsbgvx.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DSBGVX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSBGVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbgvx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbgvx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbgvx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB,
!                          LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z,
!                          LDZ, WORK, IWORK, IFAIL, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, RANGE, UPLO
!       INTEGER            IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ, M,
!      $                   N
!       DOUBLE PRECISION   ABSTOL, VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            IFAIL( * ), IWORK( * )
!       DOUBLE PRECISION   AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ),
!      $                   W( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSBGVX computes selected eigenvalues, and optionally, eigenvectors
!> of a real generalized symmetric-definite banded eigenproblem, of
!> the form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric
!> and banded, and B is also positive definite.  Eigenvalues and
!> eigenvectors can be selected by specifying either all eigenvalues,
!> a range of values or a range of indices for the desired eigenvalues.
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
!> \param[in] RANGE
!> \verbatim
!>          RANGE is CHARACTER*1
!>          = 'A': all eigenvalues will be found.
!>          = 'V': all eigenvalues in the half-open interval (VL,VU]
!>                 will be found.
!>          = 'I': the IL-th through IU-th eigenvalues will be found.
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
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
!>          If JOBZ = 'V', the n-by-n matrix used in the reduction of
!>          A*x = (lambda)*B*x to standard form, i.e. C*x = (lambda)*x,
!>          and consequently C to tridiagonal form.
!>          If JOBZ = 'N', the array Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  If JOBZ = 'N',
!>          LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in] VL
!> \verbatim
!>          VL is DOUBLE PRECISION
!>
!>          If RANGE='V', the lower bound of the interval to
!>          be searched for eigenvalues. VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>          VU is DOUBLE PRECISION
!>
!>          If RANGE='V', the upper bound of the interval to
!>          be searched for eigenvalues. VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] IL
!> \verbatim
!>          IL is INTEGER
!>
!>          If RANGE='I', the index of the
!>          smallest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] IU
!> \verbatim
!>          IU is INTEGER
!>
!>          If RANGE='I', the index of the
!>          largest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] ABSTOL
!> \verbatim
!>          ABSTOL is DOUBLE PRECISION
!>          The absolute error tolerance for the eigenvalues.
!>          An approximate eigenvalue is accepted as converged
!>          when it is determined to lie in an interval [a,b]
!>          of width less than or equal to
!>
!>                  ABSTOL + EPS *   max( |a|,|b| ) ,
!>
!>          where EPS is the machine precision.  If ABSTOL is less than
!>          or equal to zero, then  EPS*|T|  will be used in its place,
!>          where |T| is the 1-norm of the tridiagonal matrix obtained
!>          by reducing A to tridiagonal form.
!>
!>          Eigenvalues will be computed most accurately when ABSTOL is
!>          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!>          If this routine returns with INFO>0, indicating that some
!>          eigenvectors did not converge, try setting ABSTOL to
!>          2*DLAMCH('S').
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The total number of eigenvalues found.  0 <= M <= N.
!>          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
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
!>          WORK is DOUBLE PRECISION array, dimension (7*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (5*N)
!> \endverbatim
!>
!> \param[out] IFAIL
!> \verbatim
!>          IFAIL is INTEGER array, dimension (M)
!>          If JOBZ = 'V', then if INFO = 0, the first M elements of
!>          IFAIL are zero.  If INFO > 0, then IFAIL contains the
!>          indices of the eigenvalues that failed to converge.
!>          If JOBZ = 'N', then IFAIL is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          <= N: if INFO = i, then i eigenvectors failed to converge.
!>                  Their indices are stored in IFAIL.
!>          > N:  DPBSTF returned an error code; i.e.,
!>                if INFO = N + i, for 1 <= i <= N, then the leading
!>                minor of order i of B is not positive definite.
!>                The factorization of B could not be completed and
!>                no eigenvalues or eigenvectors were computed.
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
      SUBROUTINE DSBGVX(Jobz,Range,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,Q,Ldq,  &
     &                  Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz,Work,Iwork,Ifail,  &
     &                  Info)
      IMPLICIT NONE
!*--DSBGVX298
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobz , Range , Uplo
      INTEGER Il , Info , Iu , Ka , Kb , Ldab , Ldbb , Ldq , Ldz , M , N
      DOUBLE PRECISION Abstol , Vl , Vu
!     ..
!     .. Array Arguments ..
      INTEGER Ifail(*) , Iwork(*)
      DOUBLE PRECISION Ab(Ldab,*) , Bb(Ldbb,*) , Q(Ldq,*) , W(*) ,      &
     &                 Work(*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL alleig , indeig , test , upper , valeig , wantz
      CHARACTER order , vect
      INTEGER i , iinfo , indd , inde , indee , indibl , indisp ,       &
     &        indiwo , indwrk , itmp1 , j , jj , nsplit
      DOUBLE PRECISION tmp1
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DGEMV , DLACPY , DPBSTF , DSBGST , DSBTRD ,      &
     &         DSTEBZ , DSTEIN , DSTEQR , DSTERF , DSWAP , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      wantz = LSAME(Jobz,'V')
      upper = LSAME(Uplo,'U')
      alleig = LSAME(Range,'A')
      valeig = LSAME(Range,'V')
      indeig = LSAME(Range,'I')
!
      Info = 0
      IF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(alleig .OR. valeig .OR. indeig) ) THEN
         Info = -2
      ELSEIF ( .NOT.(upper .OR. LSAME(Uplo,'L')) ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Ka<0 ) THEN
         Info = -5
      ELSEIF ( Kb<0 .OR. Kb>Ka ) THEN
         Info = -6
      ELSEIF ( Ldab<Ka+1 ) THEN
         Info = -8
      ELSEIF ( Ldbb<Kb+1 ) THEN
         Info = -10
      ELSEIF ( Ldq<1 .OR. (wantz .AND. Ldq<N) ) THEN
         Info = -12
      ELSEIF ( valeig ) THEN
         IF ( N>0 .AND. Vu<=Vl ) Info = -14
      ELSEIF ( indeig ) THEN
         IF ( Il<1 .OR. Il>MAX(1,N) ) THEN
            Info = -15
         ELSEIF ( Iu<MIN(N,Il) .OR. Iu>N ) THEN
            Info = -16
         ENDIF
      ENDIF
      IF ( Info==0 ) THEN
         IF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) Info = -21
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSBGVX',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      M = 0
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
      CALL DSBGST(Jobz,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,Q,Ldq,Work,iinfo)
!
!     Reduce symmetric band matrix to tridiagonal form.
!
      indd = 1
      inde = indd + N
      indwrk = inde + N
      IF ( wantz ) THEN
         vect = 'U'
      ELSE
         vect = 'N'
      ENDIF
      CALL DSBTRD(vect,Uplo,N,Ka,Ab,Ldab,Work(indd),Work(inde),Q,Ldq,   &
     &            Work(indwrk),iinfo)
!
!     If all eigenvalues are desired and ABSTOL is less than or equal
!     to zero, then call DSTERF or SSTEQR.  If this fails for some
!     eigenvalue, then try DSTEBZ.
!
      test = .FALSE.
      IF ( indeig ) THEN
         IF ( Il==1 .AND. Iu==N ) test = .TRUE.
      ENDIF
      IF ( (alleig .OR. test) .AND. (Abstol<=ZERO) ) THEN
         CALL DCOPY(N,Work(indd),1,W,1)
         indee = indwrk + 2*N
         CALL DCOPY(N-1,Work(inde),1,Work(indee),1)
         IF ( .NOT.wantz ) THEN
            CALL DSTERF(N,W,Work(indee),Info)
         ELSE
            CALL DLACPY('A',N,N,Q,Ldq,Z,Ldz)
            CALL DSTEQR(Jobz,N,W,Work(indee),Z,Ldz,Work(indwrk),Info)
            IF ( Info==0 ) THEN
               DO i = 1 , N
                  Ifail(i) = 0
               ENDDO
            ENDIF
         ENDIF
         IF ( Info==0 ) THEN
            M = N
            GOTO 100
         ENDIF
         Info = 0
      ENDIF
!
!     Otherwise, call DSTEBZ and, if eigenvectors are desired,
!     call DSTEIN.
!
      IF ( wantz ) THEN
         order = 'B'
      ELSE
         order = 'E'
      ENDIF
      indibl = 1
      indisp = indibl + N
      indiwo = indisp + N
      CALL DSTEBZ(Range,order,N,Vl,Vu,Il,Iu,Abstol,Work(indd),Work(inde)&
     &            ,M,nsplit,W,Iwork(indibl),Iwork(indisp),Work(indwrk), &
     &            Iwork(indiwo),Info)
!
      IF ( wantz ) THEN
         CALL DSTEIN(N,Work(indd),Work(inde),M,W,Iwork(indibl),         &
     &               Iwork(indisp),Z,Ldz,Work(indwrk),Iwork(indiwo),    &
     &               Ifail,Info)
!
!        Apply transformation matrix used in reduction to tridiagonal
!        form to eigenvectors returned by DSTEIN.
!
         DO j = 1 , M
            CALL DCOPY(N,Z(1,j),1,Work(1),1)
            CALL DGEMV('N',N,N,ONE,Q,Ldq,Work,1,ZERO,Z(1,j),1)
         ENDDO
      ENDIF
!
!
!     If eigenvalues are not in order, then sort them, along with
!     eigenvectors.
!
 100  IF ( wantz ) THEN
         DO j = 1 , M - 1
            i = 0
            tmp1 = W(j)
            DO jj = j + 1 , M
               IF ( W(jj)<tmp1 ) THEN
                  i = jj
                  tmp1 = W(jj)
               ENDIF
            ENDDO
!
            IF ( i/=0 ) THEN
               itmp1 = Iwork(indibl+i-1)
               W(i) = W(j)
               Iwork(indibl+i-1) = Iwork(indibl+j-1)
               W(j) = tmp1
               Iwork(indibl+j-1) = itmp1
               CALL DSWAP(N,Z(1,i),1,Z(1,j),1)
               IF ( Info/=0 ) THEN
                  itmp1 = Ifail(i)
                  Ifail(i) = Ifail(j)
                  Ifail(j) = itmp1
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of DSBGVX
!
      END SUBROUTINE DSBGVX
