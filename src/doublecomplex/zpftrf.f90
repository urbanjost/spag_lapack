!*==zpftrf.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZPFTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZPFTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpftrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpftrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpftrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPFTRF( TRANSR, UPLO, N, A, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANSR, UPLO
!       INTEGER            N, INFO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( 0: * )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPFTRF computes the Cholesky factorization of a complex Hermitian
!> positive definite matrix A.
!>
!> The factorization has the form
!>    A = U**H * U,  if UPLO = 'U', or
!>    A = L  * L**H,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular.
!>
!> This is the block version of the algorithm, calling Level 3 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSR
!> \verbatim
!>          TRANSR is CHARACTER*1
!>          = 'N':  The Normal TRANSR of RFP A is stored;
!>          = 'C':  The Conjugate-transpose TRANSR of RFP A is stored.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of RFP A is stored;
!>          = 'L':  Lower triangle of RFP A is stored.
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
!>          A is COMPLEX*16 array, dimension ( N*(N+1)/2 );
!>          On entry, the Hermitian matrix A in RFP format. RFP format is
!>          described by TRANSR, UPLO, and N as follows: If TRANSR = 'N'
!>          then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is
!>          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'C' then RFP is
!>          the Conjugate-transpose of RFP A as defined when
!>          TRANSR = 'N'. The contents of RFP A are defined by UPLO as
!>          follows: If UPLO = 'U' the RFP A contains the nt elements of
!>          upper packed A. If UPLO = 'L' the RFP A contains the elements
!>          of lower packed A. The LDA of RFP A is (N+1)/2 when TRANSR =
!>          'C'. When TRANSR is 'N' the LDA is N+1 when N is even and N
!>          is odd. See the Note below for more details.
!>
!>          On exit, if INFO = 0, the factor U or L from the Cholesky
!>          factorization RFP A = U**H*U or RFP A = L*L**H.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the leading minor of order i is not
!>                positive definite, and the factorization could not be
!>                completed.
!>
!>  Further Notes on RFP Format:
!>  ============================
!>
!>  We first consider Standard Packed Format when N is even.
!>  We give an example where N = 6.
!>
!>     AP is Upper             AP is Lower
!>
!>   00 01 02 03 04 05       00
!>      11 12 13 14 15       10 11
!>         22 23 24 25       20 21 22
!>            33 34 35       30 31 32 33
!>               44 45       40 41 42 43 44
!>                  55       50 51 52 53 54 55
!>
!>  Let TRANSR = 'N'. RFP holds AP as follows:
!>  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last
!>  three columns of AP upper. The lower triangle A(4:6,0:2) consists of
!>  conjugate-transpose of the first three columns of AP upper.
!>  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first
!>  three columns of AP lower. The upper triangle A(0:2,0:2) consists of
!>  conjugate-transpose of the last three columns of AP lower.
!>  To denote conjugate we place -- above the element. This covers the
!>  case N even and TRANSR = 'N'.
!>
!>         RFP A                   RFP A
!>
!>                                -- -- --
!>        03 04 05                33 43 53
!>                                   -- --
!>        13 14 15                00 44 54
!>                                      --
!>        23 24 25                10 11 55
!>
!>        33 34 35                20 21 22
!>        --
!>        00 44 45                30 31 32
!>        -- --
!>        01 11 55                40 41 42
!>        -- -- --
!>        02 12 22                50 51 52
!>
!>  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-
!>  transpose of RFP A above. One therefore gets:
!>
!>           RFP A                   RFP A
!>
!>     -- -- -- --                -- -- -- -- -- --
!>     03 13 23 33 00 01 02    33 00 10 20 30 40 50
!>     -- -- -- -- --                -- -- -- -- --
!>     04 14 24 34 44 11 12    43 44 11 21 31 41 51
!>     -- -- -- -- -- --                -- -- -- --
!>     05 15 25 35 45 55 22    53 54 55 22 32 42 52
!>
!>  We next  consider Standard Packed Format when N is odd.
!>  We give an example where N = 5.
!>
!>     AP is Upper                 AP is Lower
!>
!>   00 01 02 03 04              00
!>      11 12 13 14              10 11
!>         22 23 24              20 21 22
!>            33 34              30 31 32 33
!>               44              40 41 42 43 44
!>
!>  Let TRANSR = 'N'. RFP holds AP as follows:
!>  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last
!>  three columns of AP upper. The lower triangle A(3:4,0:1) consists of
!>  conjugate-transpose of the first two   columns of AP upper.
!>  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first
!>  three columns of AP lower. The upper triangle A(0:1,1:2) consists of
!>  conjugate-transpose of the last two   columns of AP lower.
!>  To denote conjugate we place -- above the element. This covers the
!>  case N odd  and TRANSR = 'N'.
!>
!>         RFP A                   RFP A
!>
!>                                   -- --
!>        02 03 04                00 33 43
!>                                      --
!>        12 13 14                10 11 44
!>
!>        22 23 24                20 21 22
!>        --
!>        00 33 34                30 31 32
!>        -- --
!>        01 11 44                40 41 42
!>
!>  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-
!>  transpose of RFP A above. One therefore gets:
!>
!>           RFP A                   RFP A
!>
!>     -- -- --                   -- -- -- -- -- --
!>     02 12 22 00 01             00 10 20 30 40 50
!>     -- -- -- --                   -- -- -- -- --
!>     03 13 23 33 11             33 11 21 31 41 51
!>     -- -- -- -- --                   -- -- -- --
!>     04 14 24 34 44             43 44 22 32 42 52
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZPFTRF(Transr,Uplo,N,A,Info)
      IMPLICIT NONE
!*--ZPFTRF215
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Transr , Uplo
      INTEGER N , Info
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(0:*)
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      COMPLEX*16 CONE
      PARAMETER (ONE=1.0D+0,CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL lower , nisodd , normaltransr
      INTEGER n1 , n2 , k
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZHERK , ZPOTRF , ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      normaltransr = LSAME(Transr,'N')
      lower = LSAME(Uplo,'L')
      IF ( .NOT.normaltransr .AND. .NOT.LSAME(Transr,'C') ) THEN
         Info = -1
      ELSEIF ( .NOT.lower .AND. .NOT.LSAME(Uplo,'U') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZPFTRF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     If N is odd, set NISODD = .TRUE.
!     If N is even, set K = N/2 and NISODD = .FALSE.
!
      IF ( MOD(N,2)==0 ) THEN
         k = N/2
         nisodd = .FALSE.
      ELSE
         nisodd = .TRUE.
      ENDIF
!
!     Set N1 and N2 depending on LOWER
!
      IF ( lower ) THEN
         n2 = N/2
         n1 = N - n2
      ELSE
         n1 = N/2
         n2 = N - n1
      ENDIF
!
!     start execution: there are eight cases
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
!             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
!             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
!             T1 -> a(0), T2 -> a(n), S -> a(n1)
!
               CALL ZPOTRF('L',n1,A(0),N,Info)
               IF ( Info>0 ) RETURN
               CALL ZTRSM('R','L','C','N',n2,n1,CONE,A(0),N,A(n1),N)
               CALL ZHERK('U','N',n2,n1,-ONE,A(n1),N,ONE,A(N),N)
               CALL ZPOTRF('U',n2,A(N),N,Info)
               IF ( Info>0 ) Info = Info + n1
!
            ELSE
!
!             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
!             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
!             T1 -> a(n2), T2 -> a(n1), S -> a(0)
!
               CALL ZPOTRF('L',n1,A(n2),N,Info)
               IF ( Info>0 ) RETURN
               CALL ZTRSM('L','L','N','N',n1,n2,CONE,A(n2),N,A(0),N)
               CALL ZHERK('U','C',n2,n1,-ONE,A(0),N,ONE,A(n1),N)
               CALL ZPOTRF('U',n2,A(n1),N,Info)
               IF ( Info>0 ) Info = Info + n1
!
            ENDIF
!
!
!           N is odd and TRANSR = 'C'
!
         ELSEIF ( lower ) THEN
!
!              SRPA for LOWER, TRANSPOSE and N is odd
!              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
!              T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1
!
            CALL ZPOTRF('U',n1,A(0),n1,Info)
            IF ( Info>0 ) RETURN
            CALL ZTRSM('L','U','C','N',n1,n2,CONE,A(0),n1,A(n1*n1),n1)
            CALL ZHERK('L','C',n2,n1,-ONE,A(n1*n1),n1,ONE,A(1),n1)
            CALL ZPOTRF('L',n2,A(1),n1,Info)
            IF ( Info>0 ) Info = Info + n1
!
         ELSE
!
!              SRPA for UPPER, TRANSPOSE and N is odd
!              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
!              T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2
!
            CALL ZPOTRF('U',n1,A(n2*n2),n2,Info)
            IF ( Info>0 ) RETURN
            CALL ZTRSM('R','U','N','N',n2,n1,CONE,A(n2*n2),n2,A(0),n2)
            CALL ZHERK('L','N',n2,n1,-ONE,A(0),n2,ONE,A(n1*n2),n2)
            CALL ZPOTRF('L',n2,A(n1*n2),n2,Info)
            IF ( Info>0 ) Info = Info + n1
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
!              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
!              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
!              T1 -> a(1), T2 -> a(0), S -> a(k+1)
!
            CALL ZPOTRF('L',k,A(1),N+1,Info)
            IF ( Info>0 ) RETURN
            CALL ZTRSM('R','L','C','N',k,k,CONE,A(1),N+1,A(k+1),N+1)
            CALL ZHERK('U','N',k,k,-ONE,A(k+1),N+1,ONE,A(0),N+1)
            CALL ZPOTRF('U',k,A(0),N+1,Info)
            IF ( Info>0 ) Info = Info + k
!
         ELSE
!
!              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
!              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
!              T1 -> a(k+1), T2 -> a(k), S -> a(0)
!
            CALL ZPOTRF('L',k,A(k+1),N+1,Info)
            IF ( Info>0 ) RETURN
            CALL ZTRSM('L','L','N','N',k,k,CONE,A(k+1),N+1,A(0),N+1)
            CALL ZHERK('U','C',k,k,-ONE,A(0),N+1,ONE,A(k),N+1)
            CALL ZPOTRF('U',k,A(k),N+1,Info)
            IF ( Info>0 ) Info = Info + k
!
         ENDIF
!
!
!           N is even and TRANSR = 'C'
!
      ELSEIF ( lower ) THEN
!
!              SRPA for LOWER, TRANSPOSE and N is even (see paper)
!              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
!              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k
!
         CALL ZPOTRF('U',k,A(0+k),k,Info)
         IF ( Info>0 ) RETURN
         CALL ZTRSM('L','U','C','N',k,k,CONE,A(k),n1,A(k*(k+1)),k)
         CALL ZHERK('L','C',k,k,-ONE,A(k*(k+1)),k,ONE,A(0),k)
         CALL ZPOTRF('L',k,A(0),k,Info)
         IF ( Info>0 ) Info = Info + k
!
      ELSE
!
!              SRPA for UPPER, TRANSPOSE and N is even (see paper)
!              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
!              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k
!
         CALL ZPOTRF('U',k,A(k*(k+1)),k,Info)
         IF ( Info>0 ) RETURN
         CALL ZTRSM('R','U','N','N',k,k,CONE,A(k*(k+1)),k,A(0),k)
         CALL ZHERK('L','N',k,k,-ONE,A(0),k,ONE,A(k*k),k)
         CALL ZPOTRF('L',k,A(k*k),k,Info)
         IF ( Info>0 ) Info = Info + k
!
!
!
      ENDIF
!
!
!     End of ZPFTRF
!
      END SUBROUTINE ZPFTRF