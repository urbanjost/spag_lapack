!*==ctftri.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CTFTRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTFTRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctftri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctftri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctftri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTFTRI( TRANSR, UPLO, DIAG, N, A, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANSR, UPLO, DIAG
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( 0: * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTFTRI computes the inverse of a triangular matrix A stored in RFP
!> format.
!>
!> This is a Level 3 BLAS version of the algorithm.
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
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension ( N*(N+1)/2 );
!>          On entry, the triangular matrix A in RFP format. RFP format
!>          is described by TRANSR, UPLO, and N as follows: If TRANSR =
!>          'N' then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is
!>          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'C' then RFP is
!>          the Conjugate-transpose of RFP A as defined when
!>          TRANSR = 'N'. The contents of RFP A are defined by UPLO as
!>          follows: If UPLO = 'U' the RFP A contains the nt elements of
!>          upper packed A; If UPLO = 'L' the RFP A contains the nt
!>          elements of lower packed A. The LDA of RFP A is (N+1)/2 when
!>          TRANSR = 'C'. When TRANSR is 'N' the LDA is N+1 when N is
!>          even and N is odd. See the Note below for more details.
!>
!>          On exit, the (triangular) inverse of the original matrix, in
!>          the same storage format.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
!>               matrix is singular and its inverse can not be computed.
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
!> \ingroup complexOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  We first consider Standard Packed Format when N is even.
!>  We give an example where N = 6.
!>
!>      AP is Upper             AP is Lower
!>
!>   00 01 02 03 04 05       00
!>      11 12 13 14 15       10 11
!>         22 23 24 25       20 21 22
!>            33 34 35       30 31 32 33
!>               44 45       40 41 42 43 44
!>                  55       50 51 52 53 54 55
!>
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
!>
!  =====================================================================
      SUBROUTINE CTFTRI(Transr,Uplo,Diag,N,A,Info)
      USE S_CTRMM
      USE S_CTRTRI
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CTFTRI229
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Transr
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      COMPLEX , DIMENSION(0:*) :: A
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: k , n1 , n2
      LOGICAL :: lower , nisodd , normaltransr
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
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
      Info = 0
      normaltransr = LSAME(Transr,'N')
      lower = LSAME(Uplo,'L')
      IF ( .NOT.normaltransr .AND. .NOT.LSAME(Transr,'C') ) THEN
         Info = -1
      ELSEIF ( .NOT.lower .AND. .NOT.LSAME(Uplo,'U') ) THEN
         Info = -2
      ELSEIF ( .NOT.LSAME(Diag,'N') .AND. .NOT.LSAME(Diag,'U') ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CTFTRI',-Info)
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
               CALL CTRTRI('L',Diag,n1,A(0),N,Info)
               IF ( Info>0 ) RETURN
               CALL CTRMM('R','L','N',Diag,n2,n1,-CONE,A(0),N,A(n1),N)
               CALL CTRTRI('U',Diag,n2,A(N),N,Info)
               IF ( Info>0 ) Info = Info + n1
               IF ( Info>0 ) RETURN
               CALL CTRMM('L','U','C',Diag,n2,n1,CONE,A(N),N,A(n1),N)
!
            ELSE
!
!             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
!             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
!             T1 -> a(n2), T2 -> a(n1), S -> a(0)
!
               CALL CTRTRI('L',Diag,n1,A(n2),N,Info)
               IF ( Info>0 ) RETURN
               CALL CTRMM('L','L','C',Diag,n1,n2,-CONE,A(n2),N,A(0),N)
               CALL CTRTRI('U',Diag,n2,A(n1),N,Info)
               IF ( Info>0 ) Info = Info + n1
               IF ( Info>0 ) RETURN
               CALL CTRMM('R','U','N',Diag,n1,n2,CONE,A(n1),N,A(0),N)
!
            ENDIF
!
!
!           N is odd and TRANSR = 'C'
!
         ELSEIF ( lower ) THEN
!
!              SRPA for LOWER, TRANSPOSE and N is odd
!              T1 -> a(0), T2 -> a(1), S -> a(0+n1*n1)
!
            CALL CTRTRI('U',Diag,n1,A(0),n1,Info)
            IF ( Info>0 ) RETURN
            CALL CTRMM('L','U','N',Diag,n1,n2,-CONE,A(0),n1,A(n1*n1),n1)
            CALL CTRTRI('L',Diag,n2,A(1),n1,Info)
            IF ( Info>0 ) Info = Info + n1
            IF ( Info>0 ) RETURN
            CALL CTRMM('R','L','C',Diag,n1,n2,CONE,A(1),n1,A(n1*n1),n1)
!
         ELSE
!
!              SRPA for UPPER, TRANSPOSE and N is odd
!              T1 -> a(0+n2*n2), T2 -> a(0+n1*n2), S -> a(0)
!
            CALL CTRTRI('U',Diag,n1,A(n2*n2),n2,Info)
            IF ( Info>0 ) RETURN
            CALL CTRMM('R','U','C',Diag,n2,n1,-CONE,A(n2*n2),n2,A(0),n2)
            CALL CTRTRI('L',Diag,n2,A(n1*n2),n2,Info)
            IF ( Info>0 ) Info = Info + n1
            IF ( Info>0 ) RETURN
            CALL CTRMM('L','L','N',Diag,n2,n1,CONE,A(n1*n2),n2,A(0),n2)
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
            CALL CTRTRI('L',Diag,k,A(1),N+1,Info)
            IF ( Info>0 ) RETURN
            CALL CTRMM('R','L','N',Diag,k,k,-CONE,A(1),N+1,A(k+1),N+1)
            CALL CTRTRI('U',Diag,k,A(0),N+1,Info)
            IF ( Info>0 ) Info = Info + k
            IF ( Info>0 ) RETURN
            CALL CTRMM('L','U','C',Diag,k,k,CONE,A(0),N+1,A(k+1),N+1)
!
         ELSE
!
!              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
!              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
!              T1 -> a(k+1), T2 -> a(k), S -> a(0)
!
            CALL CTRTRI('L',Diag,k,A(k+1),N+1,Info)
            IF ( Info>0 ) RETURN
            CALL CTRMM('L','L','C',Diag,k,k,-CONE,A(k+1),N+1,A(0),N+1)
            CALL CTRTRI('U',Diag,k,A(k),N+1,Info)
            IF ( Info>0 ) Info = Info + k
            IF ( Info>0 ) RETURN
            CALL CTRMM('R','U','N',Diag,k,k,CONE,A(k),N+1,A(0),N+1)
         ENDIF
!
!           N is even and TRANSR = 'C'
!
      ELSEIF ( lower ) THEN
!
!              SRPA for LOWER, TRANSPOSE and N is even (see paper)
!              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
!              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k
!
         CALL CTRTRI('U',Diag,k,A(k),k,Info)
         IF ( Info>0 ) RETURN
         CALL CTRMM('L','U','N',Diag,k,k,-CONE,A(k),k,A(k*(k+1)),k)
         CALL CTRTRI('L',Diag,k,A(0),k,Info)
         IF ( Info>0 ) Info = Info + k
         IF ( Info>0 ) RETURN
         CALL CTRMM('R','L','C',Diag,k,k,CONE,A(0),k,A(k*(k+1)),k)
      ELSE
!
!              SRPA for UPPER, TRANSPOSE and N is even (see paper)
!              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
!              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k
!
         CALL CTRTRI('U',Diag,k,A(k*(k+1)),k,Info)
         IF ( Info>0 ) RETURN
         CALL CTRMM('R','U','C',Diag,k,k,-CONE,A(k*(k+1)),k,A(0),k)
         CALL CTRTRI('L',Diag,k,A(k*k),k,Info)
         IF ( Info>0 ) Info = Info + k
         IF ( Info>0 ) RETURN
         CALL CTRMM('L','L','N',Diag,k,k,CONE,A(k*k),k,A(0),k)
      ENDIF
!
!
!     End of CTFTRI
!
      END SUBROUTINE CTFTRI
