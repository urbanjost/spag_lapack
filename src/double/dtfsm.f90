!*==dtfsm.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DTFSM solves a matrix equation (one operand is a triangular matrix in RFP format).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DTFSM + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtfsm.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtfsm.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtfsm.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A,
!                         B, LDB )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANSR, DIAG, SIDE, TRANS, UPLO
!       INTEGER            LDB, M, N
!       DOUBLE PRECISION   ALPHA
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( 0: * ), B( 0: LDB-1, 0: * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Level 3 BLAS like routine for A in RFP Format.
!>
!> DTFSM  solves the matrix equation
!>
!>    op( A )*X = alpha*B  or  X*op( A ) = alpha*B
!>
!> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!>
!>    op( A ) = A   or   op( A ) = A**T.
!>
!> A is in Rectangular Full Packed (RFP) Format.
!>
!> The matrix X is overwritten on B.
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
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           On entry, SIDE specifies whether op( A ) appears on the left
!>           or right of X as follows:
!>
!>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!>
!>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the RFP matrix A came from
!>           an upper or lower triangular matrix as follows:
!>           UPLO = 'U' or 'u' RFP A came from an upper triangular matrix
!>           UPLO = 'L' or 'l' RFP A came from a  lower triangular matrix
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS  specifies the form of op( A ) to be used
!>           in the matrix multiplication as follows:
!>
!>              TRANS  = 'N' or 'n'   op( A ) = A.
!>
!>              TRANS  = 'T' or 't'   op( A ) = A'.
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not RFP A is unit
!>           triangular as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of B. M must be at
!>           least zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of B.  N must be
!>           at least zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION
!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!>           zero then  A is not referenced and  B need not be set before
!>           entry.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (NT)
!>           NT = N*(N+1)/2. On entry, the matrix A in RFP Format.
!>           RFP Format is described by TRANSR, UPLO and N as follows:
!>           If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even;
!>           K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If
!>           TRANSR = 'T' then RFP is the transpose of RFP A as
!>           defined when TRANSR = 'N'. The contents of RFP A are defined
!>           by UPLO as follows: If UPLO = 'U' the RFP A contains the NT
!>           elements of upper packed A either in normal or
!>           transpose Format. If UPLO = 'L' the RFP A contains
!>           the NT elements of lower packed A either in normal or
!>           transpose Format. The LDA of RFP A is (N+1)/2 when
!>           TRANSR = 'T'. When TRANSR is 'N' the LDA is N+1 when N is
!>           even and is N when is odd.
!>           See the Note below for more details. Unchanged on exit.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>           Before entry,  the leading  m by n part of the array  B must
!>           contain  the  right-hand  side  matrix  B,  and  on exit  is
!>           overwritten by the solution matrix  X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
!>           max( 1, m ).
!>           Unchanged on exit.
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
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  We first consider Rectangular Full Packed (RFP) Format when N is
!>  even. We give an example where N = 6.
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
!>  the transpose of the first three columns of AP upper.
!>  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first
!>  three columns of AP lower. The upper triangle A(0:2,0:2) consists of
!>  the transpose of the last three columns of AP lower.
!>  This covers the case N even and TRANSR = 'N'.
!>
!>         RFP A                   RFP A
!>
!>        03 04 05                33 43 53
!>        13 14 15                00 44 54
!>        23 24 25                10 11 55
!>        33 34 35                20 21 22
!>        00 44 45                30 31 32
!>        01 11 55                40 41 42
!>        02 12 22                50 51 52
!>
!>  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the
!>  transpose of RFP A above. One therefore gets:
!>
!>
!>           RFP A                   RFP A
!>
!>     03 13 23 33 00 01 02    33 00 10 20 30 40 50
!>     04 14 24 34 44 11 12    43 44 11 21 31 41 51
!>     05 15 25 35 45 55 22    53 54 55 22 32 42 52
!>
!>
!>  We then consider Rectangular Full Packed (RFP) Format when N is
!>  odd. We give an example where N = 5.
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
!>  the transpose of the first two columns of AP upper.
!>  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first
!>  three columns of AP lower. The upper triangle A(0:1,1:2) consists of
!>  the transpose of the last two columns of AP lower.
!>  This covers the case N odd and TRANSR = 'N'.
!>
!>         RFP A                   RFP A
!>
!>        02 03 04                00 33 43
!>        12 13 14                10 11 44
!>        22 23 24                20 21 22
!>        00 33 34                30 31 32
!>        01 11 44                40 41 42
!>
!>  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the
!>  transpose of RFP A above. One therefore gets:
!>
!>           RFP A                   RFP A
!>
!>     02 12 22 00 01             00 10 20 30 40 50
!>     03 13 23 33 11             33 11 21 31 41 51
!>     04 14 24 34 44             43 44 22 32 42 52
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE DTFSM(Transr,Side,Uplo,Trans,Diag,M,N,Alpha,A,B,Ldb)
      IMPLICIT NONE
!*--DTFSM280
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Transr , Diag , Side , Trans , Uplo
      INTEGER Ldb , M , N
      DOUBLE PRECISION Alpha
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(0:*) , B(0:Ldb-1,0:*)
!     ..
!
!  =====================================================================
!
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lower , lside , misodd , nisodd , normaltransr , notrans
      INTEGER m1 , m2 , n1 , n2 , k , info , i , j
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , DGEMM , DTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MOD
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0
      normaltransr = LSAME(Transr,'N')
      lside = LSAME(Side,'L')
      lower = LSAME(Uplo,'L')
      notrans = LSAME(Trans,'N')
      IF ( .NOT.normaltransr .AND. .NOT.LSAME(Transr,'T') ) THEN
         info = -1
      ELSEIF ( .NOT.lside .AND. .NOT.LSAME(Side,'R') ) THEN
         info = -2
      ELSEIF ( .NOT.lower .AND. .NOT.LSAME(Uplo,'U') ) THEN
         info = -3
      ELSEIF ( .NOT.notrans .AND. .NOT.LSAME(Trans,'T') ) THEN
         info = -4
      ELSEIF ( .NOT.LSAME(Diag,'N') .AND. .NOT.LSAME(Diag,'U') ) THEN
         info = -5
      ELSEIF ( M<0 ) THEN
         info = -6
      ELSEIF ( N<0 ) THEN
         info = -7
      ELSEIF ( Ldb<MAX(1,M) ) THEN
         info = -11
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('DTFSM ',-info)
         RETURN
      ENDIF
!
!     Quick return when ( (N.EQ.0).OR.(M.EQ.0) )
!
      IF ( (M==0) .OR. (N==0) ) RETURN
!
!     Quick return when ALPHA.EQ.(0D+0)
!
      IF ( Alpha==ZERO ) THEN
         DO j = 0 , N - 1
            DO i = 0 , M - 1
               B(i,j) = ZERO
            ENDDO
         ENDDO
         RETURN
      ENDIF
!
      IF ( lside ) THEN
!
!        SIDE = 'L'
!
!        A is M-by-M.
!        If M is odd, set NISODD = .TRUE., and M1 and M2.
!        If M is even, NISODD = .FALSE., and M.
!
         IF ( MOD(M,2)==0 ) THEN
            misodd = .FALSE.
            k = M/2
         ELSE
            misodd = .TRUE.
            IF ( lower ) THEN
               m2 = M/2
               m1 = M - m2
            ELSE
               m1 = M/2
               m2 = M - m1
            ENDIF
         ENDIF
!
!
         IF ( misodd ) THEN
!
!           SIDE = 'L' and N is odd
!
            IF ( normaltransr ) THEN
!
!              SIDE = 'L', N is odd, and TRANSR = 'N'
!
               IF ( lower ) THEN
!
!                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'L'
!
                  IF ( notrans ) THEN
!
!                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
!                    TRANS = 'N'
!
                     IF ( M==1 ) THEN
                        CALL DTRSM('L','L','N',Diag,m1,N,Alpha,A,M,B,   &
     &                             Ldb)
                     ELSE
                        CALL DTRSM('L','L','N',Diag,m1,N,Alpha,A(0),M,B,&
     &                             Ldb)
                        CALL DGEMM('N','N',m2,N,m1,-ONE,A(m1),M,B,Ldb,  &
     &                             Alpha,B(m1,0),Ldb)
                        CALL DTRSM('L','U','T',Diag,m2,N,ONE,A(M),M,    &
     &                             B(m1,0),Ldb)
                     ENDIF
!
!
!                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
!                    TRANS = 'T'
!
                  ELSEIF ( M==1 ) THEN
                     CALL DTRSM('L','L','T',Diag,m1,N,Alpha,A(0),M,B,   &
     &                          Ldb)
                  ELSE
                     CALL DTRSM('L','U','N',Diag,m2,N,Alpha,A(M),M,     &
     &                          B(m1,0),Ldb)
                     CALL DGEMM('T','N',m1,N,m2,-ONE,A(m1),M,B(m1,0),   &
     &                          Ldb,Alpha,B,Ldb)
                     CALL DTRSM('L','L','T',Diag,m1,N,ONE,A(0),M,B,Ldb)
!
                  ENDIF
!
!
!                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U'
!
               ELSEIF ( .NOT.notrans ) THEN
!
!                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
!                    TRANS = 'N'
!
                  CALL DTRSM('L','L','N',Diag,m1,N,Alpha,A(m2),M,B,Ldb)
                  CALL DGEMM('T','N',m2,N,m1,-ONE,A(0),M,B,Ldb,Alpha,   &
     &                       B(m1,0),Ldb)
                  CALL DTRSM('L','U','T',Diag,m2,N,ONE,A(m1),M,B(m1,0), &
     &                       Ldb)
!
               ELSE
!
!                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
!                    TRANS = 'T'
!
                  CALL DTRSM('L','U','N',Diag,m2,N,Alpha,A(m1),M,B(m1,0)&
     &                       ,Ldb)
                  CALL DGEMM('N','N',m1,N,m2,-ONE,A(0),M,B(m1,0),Ldb,   &
     &                       Alpha,B,Ldb)
                  CALL DTRSM('L','L','T',Diag,m1,N,ONE,A(m2),M,B,Ldb)
!
!
               ENDIF
!
!
!              SIDE = 'L', N is odd, and TRANSR = 'T'
!
            ELSEIF ( lower ) THEN
!
!                 SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'L'
!
               IF ( notrans ) THEN
!
!                    SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and
!                    TRANS = 'N'
!
                  IF ( M==1 ) THEN
                     CALL DTRSM('L','U','T',Diag,m1,N,Alpha,A(0),m1,B,  &
     &                          Ldb)
                  ELSE
                     CALL DTRSM('L','U','T',Diag,m1,N,Alpha,A(0),m1,B,  &
     &                          Ldb)
                     CALL DGEMM('T','N',m2,N,m1,-ONE,A(m1*m1),m1,B,Ldb, &
     &                          Alpha,B(m1,0),Ldb)
                     CALL DTRSM('L','L','N',Diag,m2,N,ONE,A(1),m1,      &
     &                          B(m1,0),Ldb)
                  ENDIF
!
!
!                    SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and
!                    TRANS = 'T'
!
               ELSEIF ( M==1 ) THEN
                  CALL DTRSM('L','U','N',Diag,m1,N,Alpha,A(0),m1,B,Ldb)
               ELSE
                  CALL DTRSM('L','L','T',Diag,m2,N,Alpha,A(1),m1,B(m1,0)&
     &                       ,Ldb)
                  CALL DGEMM('N','N',m1,N,m2,-ONE,A(m1*m1),m1,B(m1,0),  &
     &                       Ldb,Alpha,B,Ldb)
                  CALL DTRSM('L','U','N',Diag,m1,N,ONE,A(0),m1,B,Ldb)
!
               ENDIF
!
!
!                 SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'U'
!
            ELSEIF ( .NOT.notrans ) THEN
!
!                    SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and
!                    TRANS = 'N'
!
               CALL DTRSM('L','U','T',Diag,m1,N,Alpha,A(m2*m2),m2,B,Ldb)
               CALL DGEMM('N','N',m2,N,m1,-ONE,A(0),m2,B,Ldb,Alpha,     &
     &                    B(m1,0),Ldb)
               CALL DTRSM('L','L','N',Diag,m2,N,ONE,A(m1*m2),m2,B(m1,0),&
     &                    Ldb)
!
            ELSE
!
!                    SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and
!                    TRANS = 'T'
!
               CALL DTRSM('L','L','T',Diag,m2,N,Alpha,A(m1*m2),m2,      &
     &                    B(m1,0),Ldb)
               CALL DGEMM('T','N',m1,N,m2,-ONE,A(0),m2,B(m1,0),Ldb,     &
     &                    Alpha,B,Ldb)
               CALL DTRSM('L','U','N',Diag,m1,N,ONE,A(m2*m2),m2,B,Ldb)
!
!
!
            ENDIF
!
!
!           SIDE = 'L' and N is even
!
         ELSEIF ( normaltransr ) THEN
!
!              SIDE = 'L', N is even, and TRANSR = 'N'
!
            IF ( lower ) THEN
!
!                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'L'
!
               IF ( notrans ) THEN
!
!                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
!                    and TRANS = 'N'
!
                  CALL DTRSM('L','L','N',Diag,k,N,Alpha,A(1),M+1,B,Ldb)
                  CALL DGEMM('N','N',k,N,k,-ONE,A(k+1),M+1,B,Ldb,Alpha, &
     &                       B(k,0),Ldb)
                  CALL DTRSM('L','U','T',Diag,k,N,ONE,A(0),M+1,B(k,0),  &
     &                       Ldb)
!
               ELSE
!
!                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
!                    and TRANS = 'T'
!
                  CALL DTRSM('L','U','N',Diag,k,N,Alpha,A(0),M+1,B(k,0),&
     &                       Ldb)
                  CALL DGEMM('T','N',k,N,k,-ONE,A(k+1),M+1,B(k,0),Ldb,  &
     &                       Alpha,B,Ldb)
                  CALL DTRSM('L','L','T',Diag,k,N,ONE,A(1),M+1,B,Ldb)
!
               ENDIF
!
!
!                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U'
!
            ELSEIF ( .NOT.notrans ) THEN
!
!                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
!                    and TRANS = 'N'
!
               CALL DTRSM('L','L','N',Diag,k,N,Alpha,A(k+1),M+1,B,Ldb)
               CALL DGEMM('T','N',k,N,k,-ONE,A(0),M+1,B,Ldb,Alpha,B(k,0)&
     &                    ,Ldb)
               CALL DTRSM('L','U','T',Diag,k,N,ONE,A(k),M+1,B(k,0),Ldb)
!
            ELSE
!
!                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
!                    and TRANS = 'T'
               CALL DTRSM('L','U','N',Diag,k,N,Alpha,A(k),M+1,B(k,0),   &
     &                    Ldb)
               CALL DGEMM('N','N',k,N,k,-ONE,A(0),M+1,B(k,0),Ldb,Alpha, &
     &                    B,Ldb)
               CALL DTRSM('L','L','T',Diag,k,N,ONE,A(k+1),M+1,B,Ldb)
!
!
            ENDIF
!
!
!              SIDE = 'L', N is even, and TRANSR = 'T'
!
         ELSEIF ( lower ) THEN
!
!                 SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'L'
!
            IF ( notrans ) THEN
!
!                    SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L',
!                    and TRANS = 'N'
!
               CALL DTRSM('L','U','T',Diag,k,N,Alpha,A(k),k,B,Ldb)
               CALL DGEMM('T','N',k,N,k,-ONE,A(k*(k+1)),k,B,Ldb,Alpha,  &
     &                    B(k,0),Ldb)
               CALL DTRSM('L','L','N',Diag,k,N,ONE,A(0),k,B(k,0),Ldb)
!
            ELSE
!
!                    SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L',
!                    and TRANS = 'T'
!
               CALL DTRSM('L','L','T',Diag,k,N,Alpha,A(0),k,B(k,0),Ldb)
               CALL DGEMM('N','N',k,N,k,-ONE,A(k*(k+1)),k,B(k,0),Ldb,   &
     &                    Alpha,B,Ldb)
               CALL DTRSM('L','U','N',Diag,k,N,ONE,A(k),k,B,Ldb)
!
            ENDIF
!
!
!                 SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'U'
!
         ELSEIF ( .NOT.notrans ) THEN
!
!                    SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U',
!                    and TRANS = 'N'
!
            CALL DTRSM('L','U','T',Diag,k,N,Alpha,A(k*(k+1)),k,B,Ldb)
            CALL DGEMM('N','N',k,N,k,-ONE,A(0),k,B,Ldb,Alpha,B(k,0),Ldb)
            CALL DTRSM('L','L','N',Diag,k,N,ONE,A(k*k),k,B(k,0),Ldb)
!
         ELSE
!
!                    SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U',
!                    and TRANS = 'T'
!
            CALL DTRSM('L','L','T',Diag,k,N,Alpha,A(k*k),k,B(k,0),Ldb)
            CALL DGEMM('T','N',k,N,k,-ONE,A(0),k,B(k,0),Ldb,Alpha,B,Ldb)
            CALL DTRSM('L','U','N',Diag,k,N,ONE,A(k*(k+1)),k,B,Ldb)
!
!
!
!
         ENDIF
!
      ELSE
!
!        SIDE = 'R'
!
!        A is N-by-N.
!        If N is odd, set NISODD = .TRUE., and N1 and N2.
!        If N is even, NISODD = .FALSE., and K.
!
         IF ( MOD(N,2)==0 ) THEN
            nisodd = .FALSE.
            k = N/2
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
!           SIDE = 'R' and N is odd
!
            IF ( normaltransr ) THEN
!
!              SIDE = 'R', N is odd, and TRANSR = 'N'
!
               IF ( lower ) THEN
!
!                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'L'
!
                  IF ( notrans ) THEN
!
!                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
!                    TRANS = 'N'
!
                     CALL DTRSM('R','U','T',Diag,M,n2,Alpha,A(N),N,     &
     &                          B(0,n1),Ldb)
                     CALL DGEMM('N','N',M,n1,n2,-ONE,B(0,n1),Ldb,A(n1), &
     &                          N,Alpha,B(0,0),Ldb)
                     CALL DTRSM('R','L','N',Diag,M,n1,ONE,A(0),N,B(0,0),&
     &                          Ldb)
!
                  ELSE
!
!                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
!                    TRANS = 'T'
!
                     CALL DTRSM('R','L','T',Diag,M,n1,Alpha,A(0),N,     &
     &                          B(0,0),Ldb)
                     CALL DGEMM('N','T',M,n2,n1,-ONE,B(0,0),Ldb,A(n1),N,&
     &                          Alpha,B(0,n1),Ldb)
                     CALL DTRSM('R','U','N',Diag,M,n2,ONE,A(N),N,B(0,n1)&
     &                          ,Ldb)
!
                  ENDIF
!
!
!                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U'
!
               ELSEIF ( notrans ) THEN
!
!                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
!                    TRANS = 'N'
!
                  CALL DTRSM('R','L','T',Diag,M,n1,Alpha,A(n2),N,B(0,0),&
     &                       Ldb)
                  CALL DGEMM('N','N',M,n2,n1,-ONE,B(0,0),Ldb,A(0),N,    &
     &                       Alpha,B(0,n1),Ldb)
                  CALL DTRSM('R','U','N',Diag,M,n2,ONE,A(n1),N,B(0,n1), &
     &                       Ldb)
!
               ELSE
!
!                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
!                    TRANS = 'T'
!
                  CALL DTRSM('R','U','T',Diag,M,n2,Alpha,A(n1),N,B(0,n1)&
     &                       ,Ldb)
                  CALL DGEMM('N','T',M,n1,n2,-ONE,B(0,n1),Ldb,A(0),N,   &
     &                       Alpha,B(0,0),Ldb)
                  CALL DTRSM('R','L','N',Diag,M,n1,ONE,A(n2),N,B(0,0),  &
     &                       Ldb)
!
!
               ENDIF
!
!
!              SIDE = 'R', N is odd, and TRANSR = 'T'
!
            ELSEIF ( lower ) THEN
!
!                 SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'L'
!
               IF ( notrans ) THEN
!
!                    SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and
!                    TRANS = 'N'
!
                  CALL DTRSM('R','L','N',Diag,M,n2,Alpha,A(1),n1,B(0,n1)&
     &                       ,Ldb)
                  CALL DGEMM('N','T',M,n1,n2,-ONE,B(0,n1),Ldb,A(n1*n1), &
     &                       n1,Alpha,B(0,0),Ldb)
                  CALL DTRSM('R','U','T',Diag,M,n1,ONE,A(0),n1,B(0,0),  &
     &                       Ldb)
!
               ELSE
!
!                    SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and
!                    TRANS = 'T'
!
                  CALL DTRSM('R','U','N',Diag,M,n1,Alpha,A(0),n1,B(0,0),&
     &                       Ldb)
                  CALL DGEMM('N','N',M,n2,n1,-ONE,B(0,0),Ldb,A(n1*n1),  &
     &                       n1,Alpha,B(0,n1),Ldb)
                  CALL DTRSM('R','L','T',Diag,M,n2,ONE,A(1),n1,B(0,n1), &
     &                       Ldb)
!
               ENDIF
!
!
!                 SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'U'
!
            ELSEIF ( notrans ) THEN
!
!                    SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and
!                    TRANS = 'N'
!
               CALL DTRSM('R','U','N',Diag,M,n1,Alpha,A(n2*n2),n2,B(0,0)&
     &                    ,Ldb)
               CALL DGEMM('N','T',M,n2,n1,-ONE,B(0,0),Ldb,A(0),n2,Alpha,&
     &                    B(0,n1),Ldb)
               CALL DTRSM('R','L','T',Diag,M,n2,ONE,A(n1*n2),n2,B(0,n1),&
     &                    Ldb)
!
            ELSE
!
!                    SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and
!                    TRANS = 'T'
!
               CALL DTRSM('R','L','N',Diag,M,n2,Alpha,A(n1*n2),n2,      &
     &                    B(0,n1),Ldb)
               CALL DGEMM('N','N',M,n1,n2,-ONE,B(0,n1),Ldb,A(0),n2,     &
     &                    Alpha,B(0,0),Ldb)
               CALL DTRSM('R','U','T',Diag,M,n1,ONE,A(n2*n2),n2,B(0,0), &
     &                    Ldb)
!
!
!
            ENDIF
!
!
!           SIDE = 'R' and N is even
!
         ELSEIF ( normaltransr ) THEN
!
!              SIDE = 'R', N is even, and TRANSR = 'N'
!
            IF ( lower ) THEN
!
!                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'L'
!
               IF ( notrans ) THEN
!
!                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
!                    and TRANS = 'N'
!
                  CALL DTRSM('R','U','T',Diag,M,k,Alpha,A(0),N+1,B(0,k),&
     &                       Ldb)
                  CALL DGEMM('N','N',M,k,k,-ONE,B(0,k),Ldb,A(k+1),N+1,  &
     &                       Alpha,B(0,0),Ldb)
                  CALL DTRSM('R','L','N',Diag,M,k,ONE,A(1),N+1,B(0,0),  &
     &                       Ldb)
!
               ELSE
!
!                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
!                    and TRANS = 'T'
!
                  CALL DTRSM('R','L','T',Diag,M,k,Alpha,A(1),N+1,B(0,0),&
     &                       Ldb)
                  CALL DGEMM('N','T',M,k,k,-ONE,B(0,0),Ldb,A(k+1),N+1,  &
     &                       Alpha,B(0,k),Ldb)
                  CALL DTRSM('R','U','N',Diag,M,k,ONE,A(0),N+1,B(0,k),  &
     &                       Ldb)
!
               ENDIF
!
!
!                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U'
!
            ELSEIF ( notrans ) THEN
!
!                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
!                    and TRANS = 'N'
!
               CALL DTRSM('R','L','T',Diag,M,k,Alpha,A(k+1),N+1,B(0,0), &
     &                    Ldb)
               CALL DGEMM('N','N',M,k,k,-ONE,B(0,0),Ldb,A(0),N+1,Alpha, &
     &                    B(0,k),Ldb)
               CALL DTRSM('R','U','N',Diag,M,k,ONE,A(k),N+1,B(0,k),Ldb)
!
            ELSE
!
!                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
!                    and TRANS = 'T'
!
               CALL DTRSM('R','U','T',Diag,M,k,Alpha,A(k),N+1,B(0,k),   &
     &                    Ldb)
               CALL DGEMM('N','T',M,k,k,-ONE,B(0,k),Ldb,A(0),N+1,Alpha, &
     &                    B(0,0),Ldb)
               CALL DTRSM('R','L','N',Diag,M,k,ONE,A(k+1),N+1,B(0,0),   &
     &                    Ldb)
!
!
            ENDIF
!
!
!              SIDE = 'R', N is even, and TRANSR = 'T'
!
         ELSEIF ( lower ) THEN
!
!                 SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'L'
!
            IF ( notrans ) THEN
!
!                    SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L',
!                    and TRANS = 'N'
!
               CALL DTRSM('R','L','N',Diag,M,k,Alpha,A(0),k,B(0,k),Ldb)
               CALL DGEMM('N','T',M,k,k,-ONE,B(0,k),Ldb,A((k+1)*k),k,   &
     &                    Alpha,B(0,0),Ldb)
               CALL DTRSM('R','U','T',Diag,M,k,ONE,A(k),k,B(0,0),Ldb)
!
            ELSE
!
!                    SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L',
!                    and TRANS = 'T'
!
               CALL DTRSM('R','U','N',Diag,M,k,Alpha,A(k),k,B(0,0),Ldb)
               CALL DGEMM('N','N',M,k,k,-ONE,B(0,0),Ldb,A((k+1)*k),k,   &
     &                    Alpha,B(0,k),Ldb)
               CALL DTRSM('R','L','T',Diag,M,k,ONE,A(0),k,B(0,k),Ldb)
!
            ENDIF
!
!
!                 SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'U'
!
         ELSEIF ( notrans ) THEN
!
!                    SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U',
!                    and TRANS = 'N'
!
            CALL DTRSM('R','U','N',Diag,M,k,Alpha,A((k+1)*k),k,B(0,0),  &
     &                 Ldb)
            CALL DGEMM('N','T',M,k,k,-ONE,B(0,0),Ldb,A(0),k,Alpha,B(0,k)&
     &                 ,Ldb)
            CALL DTRSM('R','L','T',Diag,M,k,ONE,A(k*k),k,B(0,k),Ldb)
!
         ELSE
!
!                    SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U',
!                    and TRANS = 'T'
!
            CALL DTRSM('R','L','N',Diag,M,k,Alpha,A(k*k),k,B(0,k),Ldb)
            CALL DGEMM('N','N',M,k,k,-ONE,B(0,k),Ldb,A(0),k,Alpha,B(0,0)&
     &                 ,Ldb)
            CALL DTRSM('R','U','T',Diag,M,k,ONE,A((k+1)*k),k,B(0,0),Ldb)
!
!
!
!
         ENDIF
      ENDIF
!
!
!     End of DTFSM
!
      END SUBROUTINE DTFSM
