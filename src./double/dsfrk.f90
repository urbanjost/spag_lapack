!*==dsfrk.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DSFRK performs a symmetric rank-k operation for matrix in RFP format.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSFRK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsfrk.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsfrk.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsfrk.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA,
!                         C )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   ALPHA, BETA
!       INTEGER            K, LDA, N
!       CHARACTER          TRANS, TRANSR, UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), C( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Level 3 BLAS like routine for C in RFP Format.
!>
!> DSFRK performs one of the symmetric rank--k operations
!>
!>    C := alpha*A*A**T + beta*C,
!>
!> or
!>
!>    C := alpha*A**T*A + beta*C,
!>
!> where alpha and beta are real scalars, C is an n--by--n symmetric
!> matrix and A is an n--by--k matrix in the first case and a k--by--n
!> matrix in the second case.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSR
!> \verbatim
!>          TRANSR is CHARACTER*1
!>          = 'N':  The Normal Form of RFP A is stored;
!>          = 'T':  The Transpose Form of RFP A is stored.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On  entry, UPLO specifies whether the upper or lower
!>           triangular part of the array C is to be referenced as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   Only the upper triangular part of C
!>                                  is to be referenced.
!>
!>              UPLO = 'L' or 'l'   Only the lower triangular part of C
!>                                  is to be referenced.
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.
!>
!>              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix C. N must be
!>           at least zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry with TRANS = 'N' or 'n', K specifies the number
!>           of  columns of the matrix A, and on entry with TRANS = 'T'
!>           or 't', K specifies the number of rows of the matrix A. K
!>           must be at least zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION
!>           On entry, ALPHA specifies the scalar alpha.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,ka)
!>           where KA
!>           is K  when TRANS = 'N' or 'n', and is N otherwise. Before
!>           entry with TRANS = 'N' or 'n', the leading N--by--K part of
!>           the array A must contain the matrix A, otherwise the leading
!>           K--by--N part of the array A must contain the matrix A.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!>           then  LDA must be at least  max( 1, n ), otherwise  LDA must
!>           be at least  max( 1, k ).
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION
!>           On entry, BETA specifies the scalar beta.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (NT)
!>           NT = N*(N+1)/2. On entry, the symmetric matrix C in RFP
!>           Format. RFP Format is described by TRANSR, UPLO and N.
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
      SUBROUTINE DSFRK(Transr,Uplo,Trans,N,K,Alpha,A,Lda,Beta,C)
      USE F77KINDS                        
      USE S_DGEMM
      USE S_DSYRK
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DSFRK174
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Transr
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) :: Alpha
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) :: Beta
      REAL(R8KIND) , DIMENSION(*) :: C
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: info , j , n1 , n2 , nk , nrowa
      LOGICAL :: lower , nisodd , normaltransr , notrans
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     ..
!     .. Parameters ..
!     ..
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
      info = 0
      normaltransr = LSAME(Transr,'N')
      lower = LSAME(Uplo,'L')
      notrans = LSAME(Trans,'N')
!
      IF ( notrans ) THEN
         nrowa = N
      ELSE
         nrowa = K
      ENDIF
!
      IF ( .NOT.normaltransr .AND. .NOT.LSAME(Transr,'T') ) THEN
         info = -1
      ELSEIF ( .NOT.lower .AND. .NOT.LSAME(Uplo,'U') ) THEN
         info = -2
      ELSEIF ( .NOT.notrans .AND. .NOT.LSAME(Trans,'T') ) THEN
         info = -3
      ELSEIF ( N<0 ) THEN
         info = -4
      ELSEIF ( K<0 ) THEN
         info = -5
      ELSEIF ( Lda<MAX(1,nrowa) ) THEN
         info = -8
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('DSFRK ',-info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
!     The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not
!     done (it is in DSYRK for example) and left in the general case.
!
      IF ( (N==0) .OR. (((Alpha==ZERO) .OR. (K==0)) .AND. (Beta==ONE)) )&
     &     RETURN
!
      IF ( (Alpha==ZERO) .AND. (Beta==ZERO) ) THEN
         DO j = 1 , ((N*(N+1))/2)
            C(j) = ZERO
         ENDDO
         RETURN
      ENDIF
!
!     C is N-by-N.
!     If N is odd, set NISODD = .TRUE., and N1 and N2.
!     If N is even, NISODD = .FALSE., and NK.
!
      IF ( MOD(N,2)==0 ) THEN
         nisodd = .FALSE.
         nk = N/2
      ELSE
         nisodd = .TRUE.
         IF ( lower ) THEN
            n2 = N/2
            n1 = N - n2
         ELSE
            n1 = N/2
            n2 = N - n1
         ENDIF
      ENDIF
!
      IF ( nisodd ) THEN
!
!        N is odd
!
         IF ( normaltransr ) THEN
!
!           N is odd and TRANSR = 'N'
!
            IF ( lower ) THEN
!
!              N is odd, TRANSR = 'N', and UPLO = 'L'
!
               IF ( notrans ) THEN
!
!                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'N'
!
                  CALL DSYRK('L','N',n1,K,Alpha,A(1,1),Lda,Beta,C(1),N)
                  CALL DSYRK('U','N',n2,K,Alpha,A(n1+1,1),Lda,Beta,     &
     &                       C(N+1),N)
                  CALL DGEMM('N','T',n2,n1,K,Alpha,A(n1+1,1),Lda,A(1,1),&
     &                       Lda,Beta,C(n1+1),N)
!
               ELSE
!
!                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'T'
!
                  CALL DSYRK('L','T',n1,K,Alpha,A(1,1),Lda,Beta,C(1),N)
                  CALL DSYRK('U','T',n2,K,Alpha,A(1,n1+1),Lda,Beta,     &
     &                       C(N+1),N)
                  CALL DGEMM('T','N',n2,n1,K,Alpha,A(1,n1+1),Lda,A(1,1),&
     &                       Lda,Beta,C(n1+1),N)
!
               ENDIF
!
!
!              N is odd, TRANSR = 'N', and UPLO = 'U'
!
            ELSEIF ( notrans ) THEN
!
!                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'
!
               CALL DSYRK('L','N',n1,K,Alpha,A(1,1),Lda,Beta,C(n2+1),N)
               CALL DSYRK('U','N',n2,K,Alpha,A(n2,1),Lda,Beta,C(n1+1),N)
               CALL DGEMM('N','T',n1,n2,K,Alpha,A(1,1),Lda,A(n2,1),Lda, &
     &                    Beta,C(1),N)
!
            ELSE
!
!                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'T'
!
               CALL DSYRK('L','T',n1,K,Alpha,A(1,1),Lda,Beta,C(n2+1),N)
               CALL DSYRK('U','T',n2,K,Alpha,A(1,n2),Lda,Beta,C(n1+1),N)
               CALL DGEMM('T','N',n1,n2,K,Alpha,A(1,1),Lda,A(1,n2),Lda, &
     &                    Beta,C(1),N)
!
!
            ENDIF
!
!
!           N is odd, and TRANSR = 'T'
!
         ELSEIF ( lower ) THEN
!
!              N is odd, TRANSR = 'T', and UPLO = 'L'
!
            IF ( notrans ) THEN
!
!                 N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'N'
!
               CALL DSYRK('U','N',n1,K,Alpha,A(1,1),Lda,Beta,C(1),n1)
               CALL DSYRK('L','N',n2,K,Alpha,A(n1+1,1),Lda,Beta,C(2),n1)
               CALL DGEMM('N','T',n1,n2,K,Alpha,A(1,1),Lda,A(n1+1,1),   &
     &                    Lda,Beta,C(n1*n1+1),n1)
!
            ELSE
!
!                 N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'T'
!
               CALL DSYRK('U','T',n1,K,Alpha,A(1,1),Lda,Beta,C(1),n1)
               CALL DSYRK('L','T',n2,K,Alpha,A(1,n1+1),Lda,Beta,C(2),n1)
               CALL DGEMM('T','N',n1,n2,K,Alpha,A(1,1),Lda,A(1,n1+1),   &
     &                    Lda,Beta,C(n1*n1+1),n1)
!
            ENDIF
!
!
!              N is odd, TRANSR = 'T', and UPLO = 'U'
!
         ELSEIF ( notrans ) THEN
!
!                 N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'N'
!
            CALL DSYRK('U','N',n1,K,Alpha,A(1,1),Lda,Beta,C(n2*n2+1),n2)
            CALL DSYRK('L','N',n2,K,Alpha,A(n1+1,1),Lda,Beta,C(n1*n2+1),&
     &                 n2)
            CALL DGEMM('N','T',n2,n1,K,Alpha,A(n1+1,1),Lda,A(1,1),Lda,  &
     &                 Beta,C(1),n2)
!
         ELSE
!
!                 N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'T'
!
            CALL DSYRK('U','T',n1,K,Alpha,A(1,1),Lda,Beta,C(n2*n2+1),n2)
            CALL DSYRK('L','T',n2,K,Alpha,A(1,n1+1),Lda,Beta,C(n1*n2+1),&
     &                 n2)
            CALL DGEMM('T','N',n2,n1,K,Alpha,A(1,n1+1),Lda,A(1,1),Lda,  &
     &                 Beta,C(1),n2)
!
!
!
         ENDIF
!
!
!        N is even
!
      ELSEIF ( normaltransr ) THEN
!
!           N is even and TRANSR = 'N'
!
         IF ( lower ) THEN
!
!              N is even, TRANSR = 'N', and UPLO = 'L'
!
            IF ( notrans ) THEN
!
!                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'N'
!
               CALL DSYRK('L','N',nk,K,Alpha,A(1,1),Lda,Beta,C(2),N+1)
               CALL DSYRK('U','N',nk,K,Alpha,A(nk+1,1),Lda,Beta,C(1),   &
     &                    N+1)
               CALL DGEMM('N','T',nk,nk,K,Alpha,A(nk+1,1),Lda,A(1,1),   &
     &                    Lda,Beta,C(nk+2),N+1)
!
            ELSE
!
!                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'T'
!
               CALL DSYRK('L','T',nk,K,Alpha,A(1,1),Lda,Beta,C(2),N+1)
               CALL DSYRK('U','T',nk,K,Alpha,A(1,nk+1),Lda,Beta,C(1),   &
     &                    N+1)
               CALL DGEMM('T','N',nk,nk,K,Alpha,A(1,nk+1),Lda,A(1,1),   &
     &                    Lda,Beta,C(nk+2),N+1)
!
            ENDIF
!
!
!              N is even, TRANSR = 'N', and UPLO = 'U'
!
         ELSEIF ( notrans ) THEN
!
!                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'
!
            CALL DSYRK('L','N',nk,K,Alpha,A(1,1),Lda,Beta,C(nk+2),N+1)
            CALL DSYRK('U','N',nk,K,Alpha,A(nk+1,1),Lda,Beta,C(nk+1),   &
     &                 N+1)
            CALL DGEMM('N','T',nk,nk,K,Alpha,A(1,1),Lda,A(nk+1,1),Lda,  &
     &                 Beta,C(1),N+1)
!
         ELSE
!
!                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'T'
!
            CALL DSYRK('L','T',nk,K,Alpha,A(1,1),Lda,Beta,C(nk+2),N+1)
            CALL DSYRK('U','T',nk,K,Alpha,A(1,nk+1),Lda,Beta,C(nk+1),   &
     &                 N+1)
            CALL DGEMM('T','N',nk,nk,K,Alpha,A(1,1),Lda,A(1,nk+1),Lda,  &
     &                 Beta,C(1),N+1)
!
!
         ENDIF
!
!
!           N is even, and TRANSR = 'T'
!
      ELSEIF ( lower ) THEN
!
!              N is even, TRANSR = 'T', and UPLO = 'L'
!
         IF ( notrans ) THEN
!
!                 N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'N'
!
            CALL DSYRK('U','N',nk,K,Alpha,A(1,1),Lda,Beta,C(nk+1),nk)
            CALL DSYRK('L','N',nk,K,Alpha,A(nk+1,1),Lda,Beta,C(1),nk)
            CALL DGEMM('N','T',nk,nk,K,Alpha,A(1,1),Lda,A(nk+1,1),Lda,  &
     &                 Beta,C(((nk+1)*nk)+1),nk)
!
         ELSE
!
!                 N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'T'
!
            CALL DSYRK('U','T',nk,K,Alpha,A(1,1),Lda,Beta,C(nk+1),nk)
            CALL DSYRK('L','T',nk,K,Alpha,A(1,nk+1),Lda,Beta,C(1),nk)
            CALL DGEMM('T','N',nk,nk,K,Alpha,A(1,1),Lda,A(1,nk+1),Lda,  &
     &                 Beta,C(((nk+1)*nk)+1),nk)
!
         ENDIF
!
!
!              N is even, TRANSR = 'T', and UPLO = 'U'
!
      ELSEIF ( notrans ) THEN
!
!                 N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'N'
!
         CALL DSYRK('U','N',nk,K,Alpha,A(1,1),Lda,Beta,C(nk*(nk+1)+1),  &
     &              nk)
         CALL DSYRK('L','N',nk,K,Alpha,A(nk+1,1),Lda,Beta,C(nk*nk+1),nk)
         CALL DGEMM('N','T',nk,nk,K,Alpha,A(nk+1,1),Lda,A(1,1),Lda,Beta,&
     &              C(1),nk)
!
      ELSE
!
!                 N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'T'
!
         CALL DSYRK('U','T',nk,K,Alpha,A(1,1),Lda,Beta,C(nk*(nk+1)+1),  &
     &              nk)
         CALL DSYRK('L','T',nk,K,Alpha,A(1,nk+1),Lda,Beta,C(nk*nk+1),nk)
         CALL DGEMM('T','N',nk,nk,K,Alpha,A(1,nk+1),Lda,A(1,1),Lda,Beta,&
     &              C(1),nk)
!
!
!
!
      ENDIF
!
!
!     End of DSFRK
!
      END SUBROUTINE DSFRK
