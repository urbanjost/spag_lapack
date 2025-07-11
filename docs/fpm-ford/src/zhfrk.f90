!*==zhfrk.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZHFRK performs a Hermitian rank-k operation for matrix in RFP format.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHFRK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhfrk.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhfrk.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhfrk.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA,
!                         C )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   ALPHA, BETA
!       INTEGER            K, LDA, N
!       CHARACTER          TRANS, TRANSR, UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), C( * )
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
!> ZHFRK performs one of the Hermitian rank--k operations
!>
!>    C := alpha*A*A**H + beta*C,
!>
!> or
!>
!>    C := alpha*A**H*A + beta*C,
!>
!> where alpha and beta are real scalars, C is an n--by--n Hermitian
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
!>          = 'C':  The Conjugate-transpose Form of RFP A is stored.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!>           triangular  part  of the  array  C  is to be  referenced  as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
!>                                  is to be referenced.
!>
!>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
!>                                  is to be referenced.
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry,  TRANS  specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   C := alpha*A*A**H + beta*C.
!>
!>              TRANS = 'C' or 'c'   C := alpha*A**H*A + beta*C.
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry,  N specifies the order of the matrix C.  N must be
!>           at least zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
!>           of  columns   of  the   matrix   A,   and  on   entry   with
!>           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
!>           matrix A.  K must be at least zero.
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
!>          A is COMPLEX*16 array, dimension (LDA,ka)
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
!>          C is COMPLEX*16 array, dimension (N*(N+1)/2)
!>           On entry, the matrix A in RFP Format. RFP Format is
!>           described by TRANSR, UPLO and N. Note that the imaginary
!>           parts of the diagonal elements need not be set, they are
!>           assumed to be zero, and on exit they are set to zero.
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
!> \date June 2017
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZHFRK(Transr,Uplo,Trans,N,K,Alpha,A,Lda,Beta,C)
      IMPLICIT NONE
!*--ZHFRK171
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Alpha , Beta
      INTEGER K , Lda , N
      CHARACTER Trans , Transr , Uplo
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , C(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      COMPLEX*16 CZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      PARAMETER (CZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL lower , normaltransr , nisodd , notrans
      INTEGER info , nrowa , j , nk , n1 , n2
      COMPLEX*16 calpha , cbeta
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZGEMM , ZHERK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , DCMPLX
!     ..
!     .. Executable Statements ..
!
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
      IF ( .NOT.normaltransr .AND. .NOT.LSAME(Transr,'C') ) THEN
         info = -1
      ELSEIF ( .NOT.lower .AND. .NOT.LSAME(Uplo,'U') ) THEN
         info = -2
      ELSEIF ( .NOT.notrans .AND. .NOT.LSAME(Trans,'C') ) THEN
         info = -3
      ELSEIF ( N<0 ) THEN
         info = -4
      ELSEIF ( K<0 ) THEN
         info = -5
      ELSEIF ( Lda<MAX(1,nrowa) ) THEN
         info = -8
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('ZHFRK ',-info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
!     The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not
!     done (it is in ZHERK for example) and left in the general case.
!
      IF ( (N==0) .OR. (((Alpha==ZERO) .OR. (K==0)) .AND. (Beta==ONE)) )&
     &     RETURN
!
      IF ( (Alpha==ZERO) .AND. (Beta==ZERO) ) THEN
         DO j = 1 , ((N*(N+1))/2)
            C(j) = CZERO
         ENDDO
         RETURN
      ENDIF
!
      calpha = DCMPLX(Alpha,ZERO)
      cbeta = DCMPLX(Beta,ZERO)
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
                  CALL ZHERK('L','N',n1,K,Alpha,A(1,1),Lda,Beta,C(1),N)
                  CALL ZHERK('U','N',n2,K,Alpha,A(n1+1,1),Lda,Beta,     &
     &                       C(N+1),N)
                  CALL ZGEMM('N','C',n2,n1,K,calpha,A(n1+1,1),Lda,A(1,1)&
     &                       ,Lda,cbeta,C(n1+1),N)
!
               ELSE
!
!                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'C'
!
                  CALL ZHERK('L','C',n1,K,Alpha,A(1,1),Lda,Beta,C(1),N)
                  CALL ZHERK('U','C',n2,K,Alpha,A(1,n1+1),Lda,Beta,     &
     &                       C(N+1),N)
                  CALL ZGEMM('C','N',n2,n1,K,calpha,A(1,n1+1),Lda,A(1,1)&
     &                       ,Lda,cbeta,C(n1+1),N)
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
               CALL ZHERK('L','N',n1,K,Alpha,A(1,1),Lda,Beta,C(n2+1),N)
               CALL ZHERK('U','N',n2,K,Alpha,A(n2,1),Lda,Beta,C(n1+1),N)
               CALL ZGEMM('N','C',n1,n2,K,calpha,A(1,1),Lda,A(n2,1),Lda,&
     &                    cbeta,C(1),N)
!
            ELSE
!
!                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'C'
!
               CALL ZHERK('L','C',n1,K,Alpha,A(1,1),Lda,Beta,C(n2+1),N)
               CALL ZHERK('U','C',n2,K,Alpha,A(1,n2),Lda,Beta,C(n1+1),N)
               CALL ZGEMM('C','N',n1,n2,K,calpha,A(1,1),Lda,A(1,n2),Lda,&
     &                    cbeta,C(1),N)
!
!
            ENDIF
!
!
!           N is odd, and TRANSR = 'C'
!
         ELSEIF ( lower ) THEN
!
!              N is odd, TRANSR = 'C', and UPLO = 'L'
!
            IF ( notrans ) THEN
!
!                 N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'N'
!
               CALL ZHERK('U','N',n1,K,Alpha,A(1,1),Lda,Beta,C(1),n1)
               CALL ZHERK('L','N',n2,K,Alpha,A(n1+1,1),Lda,Beta,C(2),n1)
               CALL ZGEMM('N','C',n1,n2,K,calpha,A(1,1),Lda,A(n1+1,1),  &
     &                    Lda,cbeta,C(n1*n1+1),n1)
!
            ELSE
!
!                 N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'C'
!
               CALL ZHERK('U','C',n1,K,Alpha,A(1,1),Lda,Beta,C(1),n1)
               CALL ZHERK('L','C',n2,K,Alpha,A(1,n1+1),Lda,Beta,C(2),n1)
               CALL ZGEMM('C','N',n1,n2,K,calpha,A(1,1),Lda,A(1,n1+1),  &
     &                    Lda,cbeta,C(n1*n1+1),n1)
!
            ENDIF
!
!
!              N is odd, TRANSR = 'C', and UPLO = 'U'
!
         ELSEIF ( notrans ) THEN
!
!                 N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'N'
!
            CALL ZHERK('U','N',n1,K,Alpha,A(1,1),Lda,Beta,C(n2*n2+1),n2)
            CALL ZHERK('L','N',n2,K,Alpha,A(n1+1,1),Lda,Beta,C(n1*n2+1),&
     &                 n2)
            CALL ZGEMM('N','C',n2,n1,K,calpha,A(n1+1,1),Lda,A(1,1),Lda, &
     &                 cbeta,C(1),n2)
!
         ELSE
!
!                 N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'C'
!
            CALL ZHERK('U','C',n1,K,Alpha,A(1,1),Lda,Beta,C(n2*n2+1),n2)
            CALL ZHERK('L','C',n2,K,Alpha,A(1,n1+1),Lda,Beta,C(n1*n2+1),&
     &                 n2)
            CALL ZGEMM('C','N',n2,n1,K,calpha,A(1,n1+1),Lda,A(1,1),Lda, &
     &                 cbeta,C(1),n2)
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
               CALL ZHERK('L','N',nk,K,Alpha,A(1,1),Lda,Beta,C(2),N+1)
               CALL ZHERK('U','N',nk,K,Alpha,A(nk+1,1),Lda,Beta,C(1),   &
     &                    N+1)
               CALL ZGEMM('N','C',nk,nk,K,calpha,A(nk+1,1),Lda,A(1,1),  &
     &                    Lda,cbeta,C(nk+2),N+1)
!
            ELSE
!
!                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'C'
!
               CALL ZHERK('L','C',nk,K,Alpha,A(1,1),Lda,Beta,C(2),N+1)
               CALL ZHERK('U','C',nk,K,Alpha,A(1,nk+1),Lda,Beta,C(1),   &
     &                    N+1)
               CALL ZGEMM('C','N',nk,nk,K,calpha,A(1,nk+1),Lda,A(1,1),  &
     &                    Lda,cbeta,C(nk+2),N+1)
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
            CALL ZHERK('L','N',nk,K,Alpha,A(1,1),Lda,Beta,C(nk+2),N+1)
            CALL ZHERK('U','N',nk,K,Alpha,A(nk+1,1),Lda,Beta,C(nk+1),   &
     &                 N+1)
            CALL ZGEMM('N','C',nk,nk,K,calpha,A(1,1),Lda,A(nk+1,1),Lda, &
     &                 cbeta,C(1),N+1)
!
         ELSE
!
!                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'C'
!
            CALL ZHERK('L','C',nk,K,Alpha,A(1,1),Lda,Beta,C(nk+2),N+1)
            CALL ZHERK('U','C',nk,K,Alpha,A(1,nk+1),Lda,Beta,C(nk+1),   &
     &                 N+1)
            CALL ZGEMM('C','N',nk,nk,K,calpha,A(1,1),Lda,A(1,nk+1),Lda, &
     &                 cbeta,C(1),N+1)
!
!
         ENDIF
!
!
!           N is even, and TRANSR = 'C'
!
      ELSEIF ( lower ) THEN
!
!              N is even, TRANSR = 'C', and UPLO = 'L'
!
         IF ( notrans ) THEN
!
!                 N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'N'
!
            CALL ZHERK('U','N',nk,K,Alpha,A(1,1),Lda,Beta,C(nk+1),nk)
            CALL ZHERK('L','N',nk,K,Alpha,A(nk+1,1),Lda,Beta,C(1),nk)
            CALL ZGEMM('N','C',nk,nk,K,calpha,A(1,1),Lda,A(nk+1,1),Lda, &
     &                 cbeta,C(((nk+1)*nk)+1),nk)
!
         ELSE
!
!                 N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'C'
!
            CALL ZHERK('U','C',nk,K,Alpha,A(1,1),Lda,Beta,C(nk+1),nk)
            CALL ZHERK('L','C',nk,K,Alpha,A(1,nk+1),Lda,Beta,C(1),nk)
            CALL ZGEMM('C','N',nk,nk,K,calpha,A(1,1),Lda,A(1,nk+1),Lda, &
     &                 cbeta,C(((nk+1)*nk)+1),nk)
!
         ENDIF
!
!
!              N is even, TRANSR = 'C', and UPLO = 'U'
!
      ELSEIF ( notrans ) THEN
!
!                 N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'N'
!
         CALL ZHERK('U','N',nk,K,Alpha,A(1,1),Lda,Beta,C(nk*(nk+1)+1),  &
     &              nk)
         CALL ZHERK('L','N',nk,K,Alpha,A(nk+1,1),Lda,Beta,C(nk*nk+1),nk)
         CALL ZGEMM('N','C',nk,nk,K,calpha,A(nk+1,1),Lda,A(1,1),Lda,    &
     &              cbeta,C(1),nk)
!
      ELSE
!
!                 N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'C'
!
         CALL ZHERK('U','C',nk,K,Alpha,A(1,1),Lda,Beta,C(nk*(nk+1)+1),  &
     &              nk)
         CALL ZHERK('L','C',nk,K,Alpha,A(1,nk+1),Lda,Beta,C(nk*nk+1),nk)
         CALL ZGEMM('C','N',nk,nk,K,calpha,A(1,nk+1),Lda,A(1,1),Lda,    &
     &              cbeta,C(1),nk)
!
!
!
!
      ENDIF
!
!
!     End of ZHFRK
!
      END SUBROUTINE ZHFRK
