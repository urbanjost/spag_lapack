!*==dtfttp.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DTFTTP copies a triangular matrix from the rectangular full packed format (TF) to the standard packed format (TP).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DTFTTP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtfttp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtfttp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtfttp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTFTTP( TRANSR, UPLO, N, ARF, AP, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANSR, UPLO
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( 0: * ), ARF( 0: * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTFTTP copies a triangular matrix A from rectangular full packed
!> format (TF) to standard packed format (TP).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSR
!> \verbatim
!>          TRANSR is CHARACTER*1
!>          = 'N':  ARF is in Normal format;
!>          = 'T':  ARF is in Transpose format;
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  A is upper triangular;
!>          = 'L':  A is lower triangular.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A. N >= 0.
!> \endverbatim
!>
!> \param[in] ARF
!> \verbatim
!>          ARF is DOUBLE PRECISION array, dimension ( N*(N+1)/2 ),
!>          On entry, the upper or lower triangular matrix A stored in
!>          RFP format. For a further discussion see Notes below.
!> \endverbatim
!>
!> \param[out] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension ( N*(N+1)/2 ),
!>          On exit, the upper or lower triangular matrix A, packed
!>          columnwise in a linear array. The j-th column of A is stored
!>          in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!>
!  =====================================================================
      SUBROUTINE DTFTTP(Transr,Uplo,N,Arf,Ap,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DTFTTP194
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(0:*) :: Arf
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(0:*) :: Ap
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ij , ijp , j , jp , js , k , lda , n1 , n2 , nt
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
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      normaltransr = LSAME(Transr,'N')
      lower = LSAME(Uplo,'L')
      IF ( .NOT.normaltransr .AND. .NOT.LSAME(Transr,'T') ) THEN
         Info = -1
      ELSEIF ( .NOT.lower .AND. .NOT.LSAME(Uplo,'U') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DTFTTP',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( N==1 ) THEN
         IF ( normaltransr ) THEN
            Ap(0) = Arf(0)
         ELSE
            Ap(0) = Arf(0)
         ENDIF
         RETURN
      ENDIF
!
!     Size of array ARF(0:NT-1)
!
      nt = N*(N+1)/2
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
!     If N is odd, set NISODD = .TRUE.
!     If N is even, set K = N/2 and NISODD = .FALSE.
!
!     set lda of ARF^C; ARF^C is (0:(N+1)/2-1,0:N-noe)
!     where noe = 0 if n is even, noe = 1 if n is odd
!
      IF ( MOD(N,2)==0 ) THEN
         k = N/2
         nisodd = .FALSE.
         lda = N + 1
      ELSE
         nisodd = .TRUE.
         lda = N
      ENDIF
!
!     ARF^C has lda rows and n+1-noe cols
!
      IF ( .NOT.normaltransr ) lda = (N+1)/2
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
!             T1 -> a(0), T2 -> a(n), S -> a(n1); lda = n
!
               ijp = 0
               jp = 0
               DO j = 0 , n2
                  DO i = j , N - 1
                     ij = i + jp
                     Ap(ijp) = Arf(ij)
                     ijp = ijp + 1
                  ENDDO
                  jp = jp + lda
               ENDDO
               DO i = 0 , n2 - 1
                  DO j = 1 + i , n2
                     ij = i + j*lda
                     Ap(ijp) = Arf(ij)
                     ijp = ijp + 1
                  ENDDO
               ENDDO
!
            ELSE
!
!             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
!             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
!             T1 -> a(n2), T2 -> a(n1), S -> a(0)
!
               ijp = 0
               DO j = 0 , n1 - 1
                  ij = n2 + j
                  DO i = 0 , j
                     Ap(ijp) = Arf(ij)
                     ijp = ijp + 1
                     ij = ij + lda
                  ENDDO
               ENDDO
               js = 0
               DO j = n1 , N - 1
                  ij = js
                  DO ij = js , js + j
                     Ap(ijp) = Arf(ij)
                     ijp = ijp + 1
                  ENDDO
                  js = js + lda
               ENDDO
!
            ENDIF
!
!
!           N is odd and TRANSR = 'T'
!
         ELSEIF ( lower ) THEN
!
!              SRPA for LOWER, TRANSPOSE and N is odd
!              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
!              T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1
!
            ijp = 0
            DO i = 0 , n2
               DO ij = i*(lda+1) , N*lda - 1 , lda
                  Ap(ijp) = Arf(ij)
                  ijp = ijp + 1
               ENDDO
            ENDDO
            js = 1
            DO j = 0 , n2 - 1
               DO ij = js , js + n2 - j - 1
                  Ap(ijp) = Arf(ij)
                  ijp = ijp + 1
               ENDDO
               js = js + lda + 1
            ENDDO
!
         ELSE
!
!              SRPA for UPPER, TRANSPOSE and N is odd
!              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
!              T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2
!
            ijp = 0
            js = n2*lda
            DO j = 0 , n1 - 1
               DO ij = js , js + j
                  Ap(ijp) = Arf(ij)
                  ijp = ijp + 1
               ENDDO
               js = js + lda
            ENDDO
            DO i = 0 , n1
               DO ij = i , i + (n1+i)*lda , lda
                  Ap(ijp) = Arf(ij)
                  ijp = ijp + 1
               ENDDO
            ENDDO
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
            ijp = 0
            jp = 0
            DO j = 0 , k - 1
               DO i = j , N - 1
                  ij = 1 + i + jp
                  Ap(ijp) = Arf(ij)
                  ijp = ijp + 1
               ENDDO
               jp = jp + lda
            ENDDO
            DO i = 0 , k - 1
               DO j = i , k - 1
                  ij = i + j*lda
                  Ap(ijp) = Arf(ij)
                  ijp = ijp + 1
               ENDDO
            ENDDO
!
         ELSE
!
!              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
!              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
!              T1 -> a(k+1), T2 -> a(k), S -> a(0)
!
            ijp = 0
            DO j = 0 , k - 1
               ij = k + 1 + j
               DO i = 0 , j
                  Ap(ijp) = Arf(ij)
                  ijp = ijp + 1
                  ij = ij + lda
               ENDDO
            ENDDO
            js = 0
            DO j = k , N - 1
               ij = js
               DO ij = js , js + j
                  Ap(ijp) = Arf(ij)
                  ijp = ijp + 1
               ENDDO
               js = js + lda
            ENDDO
!
         ENDIF
!
!
!           N is even and TRANSR = 'T'
!
      ELSEIF ( lower ) THEN
!
!              SRPA for LOWER, TRANSPOSE and N is even (see paper)
!              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
!              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k
!
         ijp = 0
         DO i = 0 , k - 1
            DO ij = i + (i+1)*lda , (N+1)*lda - 1 , lda
               Ap(ijp) = Arf(ij)
               ijp = ijp + 1
            ENDDO
         ENDDO
         js = 0
         DO j = 0 , k - 1
            DO ij = js , js + k - j - 1
               Ap(ijp) = Arf(ij)
               ijp = ijp + 1
            ENDDO
            js = js + lda + 1
         ENDDO
!
      ELSE
!
!              SRPA for UPPER, TRANSPOSE and N is even (see paper)
!              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
!              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k
!
         ijp = 0
         js = (k+1)*lda
         DO j = 0 , k - 1
            DO ij = js , js + j
               Ap(ijp) = Arf(ij)
               ijp = ijp + 1
            ENDDO
            js = js + lda
         ENDDO
         DO i = 0 , k - 1
            DO ij = i , i + (k+i)*lda , lda
               Ap(ijp) = Arf(ij)
               ijp = ijp + 1
            ENDDO
         ENDDO
!
!
!
      ENDIF
!
!
!     End of DTFTTP
!
      END SUBROUTINE DTFTTP
