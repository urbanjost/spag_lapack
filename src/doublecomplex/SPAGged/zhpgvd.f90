!*==zhpgvd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZHPGVD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHPGVD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpgvd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpgvd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpgvd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,
!                          LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, ITYPE, LDZ, LIWORK, LRWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   RWORK( * ), W( * )
!       COMPLEX*16         AP( * ), BP( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHPGVD computes all the eigenvalues and, optionally, the eigenvectors
!> of a complex generalized Hermitian-definite eigenproblem, of the form
!> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
!> B are assumed to be Hermitian, stored in packed format, and B is also
!> positive definite.
!> If eigenvectors are desired, it uses a divide and conquer algorithm.
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
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          Specifies the problem type to be solved:
!>          = 1:  A*x = (lambda)*B*x
!>          = 2:  A*B*x = (lambda)*x
!>          = 3:  B*A*x = (lambda)*x
!> \endverbatim
!>
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
!> \param[in,out] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the Hermitian matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!>
!>          On exit, the contents of AP are destroyed.
!> \endverbatim
!>
!> \param[in,out] BP
!> \verbatim
!>          BP is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the Hermitian matrix
!>          B, packed columnwise in a linear array.  The j-th column of B
!>          is stored in the array BP as follows:
!>          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;
!>          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.
!>
!>          On exit, the triangular factor U or L from the Cholesky
!>          factorization B = U**H*U or B = L*L**H, in the same storage
!>          format as B.
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
!>          Z is COMPLEX*16 array, dimension (LDZ, N)
!>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
!>          eigenvectors.  The eigenvectors are normalized as follows:
!>          if ITYPE = 1 or 2, Z**H*B*Z = I;
!>          if ITYPE = 3, Z**H*inv(B)*Z = I.
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
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the required LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If N <= 1,               LWORK >= 1.
!>          If JOBZ = 'N' and N > 1, LWORK >= N.
!>          If JOBZ = 'V' and N > 1, LWORK >= 2*N.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the required sizes of the WORK, RWORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK))
!>          On exit, if INFO = 0, RWORK(1) returns the required LRWORK.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The dimension of array RWORK.
!>          If N <= 1,               LRWORK >= 1.
!>          If JOBZ = 'N' and N > 1, LRWORK >= N.
!>          If JOBZ = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2.
!>
!>          If LRWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the required sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the required LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of array IWORK.
!>          If JOBZ  = 'N' or N <= 1, LIWORK >= 1.
!>          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the required sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  ZPPTRF or ZHPEVD returned an error code:
!>             <= N:  if INFO = i, ZHPEVD failed to converge;
!>                    i off-diagonal elements of an intermediate
!>                    tridiagonal form did not convergeto zero;
!>             > N:   if INFO = N + i, for 1 <= i <= n, then the leading
!>                    minor of order i of B is not positive definite.
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
!> \ingroup complex16OTHEReigen
!
!> \par Contributors:
!  ==================
!>
!>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
!
!  =====================================================================
      SUBROUTINE ZHPGVD(Itype,Jobz,Uplo,N,Ap,Bp,W,Z,Ldz,Work,Lwork,     &
     &                  Rwork,Lrwork,Iwork,Liwork,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZHPEVD
      USE S_ZHPGST
      USE S_ZPPTRF
      USE S_ZTPMV
      USE S_ZTPSV
      IMPLICIT NONE
!*--ZHPGVD243
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(*) :: Bp
      REAL(R8KIND) , DIMENSION(*) :: W
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER :: Lrwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: j , liwmin , lrwmin , lwmin , neig
      LOGICAL :: lquery , upper , wantz
      CHARACTER :: trans
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      wantz = LSAME(Jobz,'V')
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1 .OR. Lrwork==-1 .OR. Liwork==-1)
!
      Info = 0
      IF ( Itype<1 .OR. Itype>3 ) THEN
         Info = -1
      ELSEIF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -2
      ELSEIF ( .NOT.(upper .OR. LSAME(Uplo,'L')) ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) THEN
         Info = -9
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( N<=1 ) THEN
            lwmin = 1
            liwmin = 1
            lrwmin = 1
         ELSEIF ( wantz ) THEN
            lwmin = 2*N
            lrwmin = 1 + 5*N + 2*N**2
            liwmin = 3 + 5*N
         ELSE
            lwmin = N
            lrwmin = N
            liwmin = 1
         ENDIF
!
         Work(1) = lwmin
         Rwork(1) = lrwmin
         Iwork(1) = liwmin
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
         CALL XERBLA('ZHPGVD',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Form a Cholesky factorization of B.
!
      CALL ZPPTRF(Uplo,N,Bp,Info)
      IF ( Info/=0 ) THEN
         Info = N + Info
         RETURN
      ENDIF
!
!     Transform problem to standard eigenvalue problem and solve.
!
      CALL ZHPGST(Itype,Uplo,N,Ap,Bp,Info)
      CALL ZHPEVD(Jobz,Uplo,N,Ap,W,Z,Ldz,Work,Lwork,Rwork,Lrwork,Iwork, &
     &            Liwork,Info)
      lwmin = MAX(DBLE(lwmin),DBLE(Work(1)))
      lrwmin = MAX(DBLE(lrwmin),DBLE(Rwork(1)))
      liwmin = MAX(DBLE(liwmin),DBLE(Iwork(1)))
!
      IF ( wantz ) THEN
!
!        Backtransform eigenvectors to the original problem.
!
         neig = N
         IF ( Info>0 ) neig = Info - 1
         IF ( Itype==1 .OR. Itype==2 ) THEN
!
!           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
!           backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y
!
            IF ( upper ) THEN
               trans = 'N'
            ELSE
               trans = 'C'
            ENDIF
!
            DO j = 1 , neig
               CALL ZTPSV(Uplo,trans,'Non-unit',N,Bp,Z(1,j),1)
            ENDDO
!
         ELSEIF ( Itype==3 ) THEN
!
!           For B*A*x=(lambda)*x;
!           backtransform eigenvectors: x = L*y or U**H *y
!
            IF ( upper ) THEN
               trans = 'C'
            ELSE
               trans = 'N'
            ENDIF
!
            DO j = 1 , neig
               CALL ZTPMV(Uplo,trans,'Non-unit',N,Bp,Z(1,j),1)
            ENDDO
         ENDIF
      ENDIF
!
      Work(1) = lwmin
      Rwork(1) = lrwmin
      Iwork(1) = liwmin
!
!     End of ZHPGVD
!
      END SUBROUTINE ZHPGVD
