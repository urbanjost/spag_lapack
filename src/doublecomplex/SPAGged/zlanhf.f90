!*==zlanhf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLANHF returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a Hermitian matrix in RFP format.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLANHF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlanhf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlanhf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlanhf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZLANHF( NORM, TRANSR, UPLO, N, A, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM, TRANSR, UPLO
!       INTEGER            N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   WORK( 0: * )
!       COMPLEX*16         A( 0: * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLANHF  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> complex Hermitian matrix A in RFP format.
!> \endverbatim
!>
!> \return ZLANHF
!> \verbatim
!>
!>    ZLANHF = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!>             (
!>             ( norm1(A),         NORM = '1', 'O' or 'o'
!>             (
!>             ( normI(A),         NORM = 'I' or 'i'
!>             (
!>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normF  denotes the  Frobenius norm of a matrix (square root of sum of
!> squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER
!>            Specifies the value to be returned in ZLANHF as described
!>            above.
!> \endverbatim
!>
!> \param[in] TRANSR
!> \verbatim
!>          TRANSR is CHARACTER
!>            Specifies whether the RFP format of A is normal or
!>            conjugate-transposed format.
!>            = 'N':  RFP format is Normal
!>            = 'C':  RFP format is Conjugate-transposed
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER
!>            On entry, UPLO specifies whether the RFP matrix A came from
!>            an upper or lower triangular matrix as follows:
!>
!>            UPLO = 'U' or 'u' RFP A came from an upper triangular
!>            matrix
!>
!>            UPLO = 'L' or 'l' RFP A came from a  lower triangular
!>            matrix
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>            The order of the matrix A.  N >= 0.  When N = 0, ZLANHF is
!>            set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( N*(N+1)/2 );
!>            On entry, the matrix A in RFP Format.
!>            RFP Format is described by TRANSR, UPLO and N as follows:
!>            If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even;
!>            K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If
!>            TRANSR = 'C' then RFP is the Conjugate-transpose of RFP A
!>            as defined when TRANSR = 'N'. The contents of RFP A are
!>            defined by UPLO as follows: If UPLO = 'U' the RFP A
!>            contains the ( N*(N+1)/2 ) elements of upper packed A
!>            either in normal or conjugate-transpose Format. If
!>            UPLO = 'L' the RFP A contains the ( N*(N+1) /2 ) elements
!>            of lower packed A either in normal or conjugate-transpose
!>            Format. The LDA of RFP A is (N+1)/2 when TRANSR = 'C'. When
!>            TRANSR is 'N' the LDA is N+1 when N is even and is N when
!>            is odd. See the Note below for more details.
!>            Unchanged on exit.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK),
!>            where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!>            WORK is not referenced.
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
!> \ingroup complex16OTHERcomputational
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
      FUNCTION ZLANHF(Norm,Transr,Uplo,N,A,Work)
      USE F77KINDS                        
      USE S_DISNAN
      USE S_LSAME
      USE S_ZLASSQ
      IMPLICIT NONE
!*--ZLANHF254
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: ZLANHF
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Norm
      CHARACTER :: Transr
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , DIMENSION(0:*) :: A
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(0:*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: aa , s , scale , temp , value
      INTEGER :: i , ifm , ilu , j , k , l , lda , n1 , noe
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
      IF ( N==0 ) THEN
         ZLANHF = ZERO
         RETURN
      ELSEIF ( N==1 ) THEN
         ZLANHF = ABS(DBLE(A(0)))
         RETURN
      ENDIF
!
!     set noe = 1 if n is odd. if n is even set noe=0
!
      noe = 1
      IF ( MOD(N,2)==0 ) noe = 0
!
!     set ifm = 0 when form='C' or 'c' and 1 otherwise
!
      ifm = 1
      IF ( LSAME(Transr,'C') ) ifm = 0
!
!     set ilu = 0 when uplo='U or 'u' and 1 otherwise
!
      ilu = 1
      IF ( LSAME(Uplo,'U') ) ilu = 0
!
!     set lda = (n+1)/2 when ifm = 0
!     set lda = n when ifm = 1 and noe = 1
!     set lda = n+1 when ifm = 1 and noe = 0
!
      IF ( ifm/=1 ) THEN
!        ifm=0
         lda = (N+1)/2
      ELSEIF ( noe==1 ) THEN
         lda = N
      ELSE
!           noe=0
         lda = N + 1
      ENDIF
!
      IF ( LSAME(Norm,'M') ) THEN
!
!       Find max(abs(A(i,j))).
!
         k = (N+1)/2
         value = ZERO
         IF ( noe==1 ) THEN
!           n is odd & n = k + k - 1
            IF ( ifm==1 ) THEN
!              A is n by k
               IF ( ilu==1 ) THEN
!                 uplo ='L'
                  j = 0
!                 -> L(0,0)
                  temp = ABS(DBLE(A(j+j*lda)))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  DO i = 1 , N - 1
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
                  DO j = 1 , k - 1
                     DO i = 0 , j - 2
                        temp = ABS(A(i+j*lda))
                        IF ( value<temp .OR. DISNAN(temp) ) value = temp
                     ENDDO
                     i = j - 1
!                    L(k+j,k+j)
                     temp = ABS(DBLE(A(i+j*lda)))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                     i = j
!                    -> L(j,j)
                     temp = ABS(DBLE(A(i+j*lda)))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                     DO i = j + 1 , N - 1
                        temp = ABS(A(i+j*lda))
                        IF ( value<temp .OR. DISNAN(temp) ) value = temp
                     ENDDO
                  ENDDO
               ELSE
!                 uplo = 'U'
                  DO j = 0 , k - 2
                     DO i = 0 , k + j - 2
                        temp = ABS(A(i+j*lda))
                        IF ( value<temp .OR. DISNAN(temp) ) value = temp
                     ENDDO
                     i = k + j - 1
!                    -> U(i,i)
                     temp = ABS(DBLE(A(i+j*lda)))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                     i = i + 1
!                    =k+j; i -> U(j,j)
                     temp = ABS(DBLE(A(i+j*lda)))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                     DO i = k + j + 1 , N - 1
                        temp = ABS(A(i+j*lda))
                        IF ( value<temp .OR. DISNAN(temp) ) value = temp
                     ENDDO
                  ENDDO
                  DO i = 0 , N - 2
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
!                    j=k-1
                  ENDDO
!                 i=n-1 -> U(n-1,n-1)
                  temp = ABS(DBLE(A(i+j*lda)))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDIF
!              xpose case; A is k by n
            ELSEIF ( ilu==1 ) THEN
!                 uplo ='L'
               DO j = 0 , k - 2
                  DO i = 0 , j - 1
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
                  i = j
!                    L(i,i)
                  temp = ABS(DBLE(A(i+j*lda)))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  i = j + 1
!                    L(j+k,j+k)
                  temp = ABS(DBLE(A(i+j*lda)))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  DO i = j + 2 , k - 1
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
               ENDDO
               j = k - 1
               DO i = 0 , k - 2
                  temp = ABS(A(i+j*lda))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
               i = k - 1
!                 -> L(i,i) is at A(i,j)
               temp = ABS(DBLE(A(i+j*lda)))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
               DO j = k , N - 1
                  DO i = 0 , k - 1
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
               ENDDO
            ELSE
!                 uplo = 'U'
               DO j = 0 , k - 2
                  DO i = 0 , k - 1
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
               ENDDO
               j = k - 1
!                 -> U(j,j) is at A(0,j)
               temp = ABS(DBLE(A(0+j*lda)))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
               DO i = 1 , k - 1
                  temp = ABS(A(i+j*lda))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
               DO j = k , N - 1
                  DO i = 0 , j - k - 1
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
                  i = j - k
!                    -> U(i,i) at A(i,j)
                  temp = ABS(DBLE(A(i+j*lda)))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  i = j - k + 1
!                    U(j,j)
                  temp = ABS(DBLE(A(i+j*lda)))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  DO i = j - k + 2 , k - 1
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
               ENDDO
            ENDIF
!           n is even & k = n/2
         ELSEIF ( ifm==1 ) THEN
!              A is n+1 by k
            IF ( ilu==1 ) THEN
!                 uplo ='L'
               j = 0
!                 -> L(k,k) & j=1 -> L(0,0)
               temp = ABS(DBLE(A(j+j*lda)))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
               temp = ABS(DBLE(A(j+1+j*lda)))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
               DO i = 2 , N
                  temp = ABS(A(i+j*lda))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
               DO j = 1 , k - 1
                  DO i = 0 , j - 1
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
                  i = j
!                    L(k+j,k+j)
                  temp = ABS(DBLE(A(i+j*lda)))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  i = j + 1
!                    -> L(j,j)
                  temp = ABS(DBLE(A(i+j*lda)))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  DO i = j + 2 , N
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
               ENDDO
            ELSE
!                 uplo = 'U'
               DO j = 0 , k - 2
                  DO i = 0 , k + j - 1
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
                  i = k + j
!                    -> U(i,i)
                  temp = ABS(DBLE(A(i+j*lda)))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  i = i + 1
!                    =k+j+1; i -> U(j,j)
                  temp = ABS(DBLE(A(i+j*lda)))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  DO i = k + j + 2 , N
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
               ENDDO
               DO i = 0 , N - 2
                  temp = ABS(A(i+j*lda))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
!                    j=k-1
               ENDDO
!                 i=n-1 -> U(n-1,n-1)
               temp = ABS(DBLE(A(i+j*lda)))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
               i = N
!                 -> U(k-1,k-1)
               temp = ABS(DBLE(A(i+j*lda)))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
            ENDIF
!              xpose case; A is k by n+1
         ELSEIF ( ilu==1 ) THEN
!                 uplo ='L'
            j = 0
!                 -> L(k,k) at A(0,0)
            temp = ABS(DBLE(A(j+j*lda)))
            IF ( value<temp .OR. DISNAN(temp) ) value = temp
            DO i = 1 , k - 1
               temp = ABS(A(i+j*lda))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
            ENDDO
            DO j = 1 , k - 1
               DO i = 0 , j - 2
                  temp = ABS(A(i+j*lda))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
               i = j - 1
!                    L(i,i)
               temp = ABS(DBLE(A(i+j*lda)))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
               i = j
!                    L(j+k,j+k)
               temp = ABS(DBLE(A(i+j*lda)))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
               DO i = j + 1 , k - 1
                  temp = ABS(A(i+j*lda))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
            ENDDO
            j = k
            DO i = 0 , k - 2
               temp = ABS(A(i+j*lda))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
            ENDDO
            i = k - 1
!                 -> L(i,i) is at A(i,j)
            temp = ABS(DBLE(A(i+j*lda)))
            IF ( value<temp .OR. DISNAN(temp) ) value = temp
            DO j = k + 1 , N
               DO i = 0 , k - 1
                  temp = ABS(A(i+j*lda))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
            ENDDO
         ELSE
!                 uplo = 'U'
            DO j = 0 , k - 1
               DO i = 0 , k - 1
                  temp = ABS(A(i+j*lda))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
            ENDDO
            j = k
!                 -> U(j,j) is at A(0,j)
            temp = ABS(DBLE(A(0+j*lda)))
            IF ( value<temp .OR. DISNAN(temp) ) value = temp
            DO i = 1 , k - 1
               temp = ABS(A(i+j*lda))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
            ENDDO
            DO j = k + 1 , N - 1
               DO i = 0 , j - k - 2
                  temp = ABS(A(i+j*lda))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
               i = j - k - 1
!                    -> U(i,i) at A(i,j)
               temp = ABS(DBLE(A(i+j*lda)))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
               i = j - k
!                    U(j,j)
               temp = ABS(DBLE(A(i+j*lda)))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
               DO i = j - k + 1 , k - 1
                  temp = ABS(A(i+j*lda))
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
            ENDDO
            j = N
            DO i = 0 , k - 2
               temp = ABS(A(i+j*lda))
               IF ( value<temp .OR. DISNAN(temp) ) value = temp
            ENDDO
            i = k - 1
!                 U(k,k) at A(i,j)
            temp = ABS(DBLE(A(i+j*lda)))
            IF ( value<temp .OR. DISNAN(temp) ) value = temp
         ENDIF
      ELSEIF ( (LSAME(Norm,'I')) .OR. (LSAME(Norm,'O')) .OR. (Norm=='1')&
     &         ) THEN
!
!       Find normI(A) ( = norm1(A), since A is Hermitian).
!
         IF ( ifm==1 ) THEN
!           A is 'N'
            k = N/2
            IF ( noe==1 ) THEN
!              n is odd & A is n by (n+1)/2
               IF ( ilu==0 ) THEN
!                 uplo = 'U'
                  DO i = 0 , k - 1
                     Work(i) = ZERO
                  ENDDO
                  DO j = 0 , k
                     s = ZERO
                     DO i = 0 , k + j - 1
                        aa = ABS(A(i+j*lda))
!                       -> A(i,j+k)
                        s = s + aa
                        Work(i) = Work(i) + aa
                     ENDDO
                     aa = ABS(DBLE(A(i+j*lda)))
!                    -> A(j+k,j+k)
                     Work(j+k) = s + aa
                     IF ( i==k+k ) EXIT
                     i = i + 1
                     aa = ABS(DBLE(A(i+j*lda)))
!                    -> A(j,j)
                     Work(j) = Work(j) + aa
                     s = ZERO
                     DO l = j + 1 , k - 1
                        i = i + 1
                        aa = ABS(A(i+j*lda))
!                       -> A(l,j)
                        s = s + aa
                        Work(l) = Work(l) + aa
                     ENDDO
                     Work(j) = Work(j) + s
                  ENDDO
                  value = Work(0)
                  DO i = 1 , N - 1
                     temp = Work(i)
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
               ELSE
!                 ilu = 1 & uplo = 'L'
                  k = k + 1
!                 k=(n+1)/2 for n odd and ilu=1
                  DO i = k , N - 1
                     Work(i) = ZERO
                  ENDDO
                  DO j = k - 1 , 0 , -1
                     s = ZERO
                     DO i = 0 , j - 2
                        aa = ABS(A(i+j*lda))
!                       -> A(j+k,i+k)
                        s = s + aa
                        Work(i+k) = Work(i+k) + aa
                     ENDDO
                     IF ( j>0 ) THEN
                        aa = ABS(DBLE(A(i+j*lda)))
!                       -> A(j+k,j+k)
                        s = s + aa
                        Work(i+k) = Work(i+k) + s
!                       i=j
                        i = i + 1
                     ENDIF
                     aa = ABS(DBLE(A(i+j*lda)))
!                    -> A(j,j)
                     Work(j) = aa
                     s = ZERO
                     DO l = j + 1 , N - 1
                        i = i + 1
                        aa = ABS(A(i+j*lda))
!                       -> A(l,j)
                        s = s + aa
                        Work(l) = Work(l) + aa
                     ENDDO
                     Work(j) = Work(j) + s
                  ENDDO
                  value = Work(0)
                  DO i = 1 , N - 1
                     temp = Work(i)
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
               ENDIF
!              n is even & A is n+1 by k = n/2
            ELSEIF ( ilu==0 ) THEN
!                 uplo = 'U'
               DO i = 0 , k - 1
                  Work(i) = ZERO
               ENDDO
               DO j = 0 , k - 1
                  s = ZERO
                  DO i = 0 , k + j - 1
                     aa = ABS(A(i+j*lda))
!                       -> A(i,j+k)
                     s = s + aa
                     Work(i) = Work(i) + aa
                  ENDDO
                  aa = ABS(DBLE(A(i+j*lda)))
!                    -> A(j+k,j+k)
                  Work(j+k) = s + aa
                  i = i + 1
                  aa = ABS(DBLE(A(i+j*lda)))
!                    -> A(j,j)
                  Work(j) = Work(j) + aa
                  s = ZERO
                  DO l = j + 1 , k - 1
                     i = i + 1
                     aa = ABS(A(i+j*lda))
!                       -> A(l,j)
                     s = s + aa
                     Work(l) = Work(l) + aa
                  ENDDO
                  Work(j) = Work(j) + s
               ENDDO
               value = Work(0)
               DO i = 1 , N - 1
                  temp = Work(i)
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
            ELSE
!                 ilu = 1 & uplo = 'L'
               DO i = k , N - 1
                  Work(i) = ZERO
               ENDDO
               DO j = k - 1 , 0 , -1
                  s = ZERO
                  DO i = 0 , j - 1
                     aa = ABS(A(i+j*lda))
!                       -> A(j+k,i+k)
                     s = s + aa
                     Work(i+k) = Work(i+k) + aa
                  ENDDO
                  aa = ABS(DBLE(A(i+j*lda)))
!                    -> A(j+k,j+k)
                  s = s + aa
                  Work(i+k) = Work(i+k) + s
!                    i=j
                  i = i + 1
                  aa = ABS(DBLE(A(i+j*lda)))
!                    -> A(j,j)
                  Work(j) = aa
                  s = ZERO
                  DO l = j + 1 , N - 1
                     i = i + 1
                     aa = ABS(A(i+j*lda))
!                       -> A(l,j)
                     s = s + aa
                     Work(l) = Work(l) + aa
                  ENDDO
                  Work(j) = Work(j) + s
               ENDDO
               value = Work(0)
               DO i = 1 , N - 1
                  temp = Work(i)
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
            ENDIF
         ELSE
!           ifm=0
            k = N/2
            IF ( noe==1 ) THEN
!              n is odd & A is (n+1)/2 by n
               IF ( ilu==0 ) THEN
!                 uplo = 'U'
                  n1 = k
!                 n/2
                  k = k + 1
!                 k is the row size and lda
                  DO i = n1 , N - 1
                     Work(i) = ZERO
                  ENDDO
                  DO j = 0 , n1 - 1
                     s = ZERO
                     DO i = 0 , k - 1
                        aa = ABS(A(i+j*lda))
!                       A(j,n1+i)
                        Work(i+n1) = Work(i+n1) + aa
                        s = s + aa
                     ENDDO
                     Work(j) = s
                  ENDDO
!                 j=n1=k-1 is special
                  s = ABS(DBLE(A(0+j*lda)))
!                 A(k-1,k-1)
                  DO i = 1 , k - 1
                     aa = ABS(A(i+j*lda))
!                    A(k-1,i+n1)
                     Work(i+n1) = Work(i+n1) + aa
                     s = s + aa
                  ENDDO
                  Work(j) = Work(j) + s
                  DO j = k , N - 1
                     s = ZERO
                     DO i = 0 , j - k - 1
                        aa = ABS(A(i+j*lda))
!                       A(i,j-k)
                        Work(i) = Work(i) + aa
                        s = s + aa
                     ENDDO
!                    i=j-k
                     aa = ABS(DBLE(A(i+j*lda)))
!                    A(j-k,j-k)
                     s = s + aa
                     Work(j-k) = Work(j-k) + s
                     i = i + 1
                     s = ABS(DBLE(A(i+j*lda)))
!                    A(j,j)
                     DO l = j + 1 , N - 1
                        i = i + 1
                        aa = ABS(A(i+j*lda))
!                       A(j,l)
                        Work(l) = Work(l) + aa
                        s = s + aa
                     ENDDO
                     Work(j) = Work(j) + s
                  ENDDO
                  value = Work(0)
                  DO i = 1 , N - 1
                     temp = Work(i)
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
               ELSE
!                 ilu=1 & uplo = 'L'
                  k = k + 1
!                 k=(n+1)/2 for n odd and ilu=1
                  DO i = k , N - 1
                     Work(i) = ZERO
                  ENDDO
                  DO j = 0 , k - 2
!                    process
                     s = ZERO
                     DO i = 0 , j - 1
                        aa = ABS(A(i+j*lda))
!                       A(j,i)
                        Work(i) = Work(i) + aa
                        s = s + aa
                     ENDDO
                     aa = ABS(DBLE(A(i+j*lda)))
!                    i=j so process of A(j,j)
                     s = s + aa
                     Work(j) = s
!                    is initialised here
                     i = i + 1
!                    i=j process A(j+k,j+k)
                     aa = ABS(DBLE(A(i+j*lda)))
                     s = aa
                     DO l = k + j + 1 , N - 1
                        i = i + 1
                        aa = ABS(A(i+j*lda))
!                       A(l,k+j)
                        s = s + aa
                        Work(l) = Work(l) + aa
                     ENDDO
                     Work(k+j) = Work(k+j) + s
                  ENDDO
!                 j=k-1 is special :process col A(k-1,0:k-1)
                  s = ZERO
                  DO i = 0 , k - 2
                     aa = ABS(A(i+j*lda))
!                    A(k,i)
                     Work(i) = Work(i) + aa
                     s = s + aa
                  ENDDO
!                 i=k-1
                  aa = ABS(DBLE(A(i+j*lda)))
!                 A(k-1,k-1)
                  s = s + aa
                  Work(i) = s
!                 done with col j=k+1
                  DO j = k , N - 1
!                    process col j of A = A(j,0:k-1)
                     s = ZERO
                     DO i = 0 , k - 1
                        aa = ABS(A(i+j*lda))
!                       A(j,i)
                        Work(i) = Work(i) + aa
                        s = s + aa
                     ENDDO
                     Work(j) = Work(j) + s
                  ENDDO
                  value = Work(0)
                  DO i = 1 , N - 1
                     temp = Work(i)
                     IF ( value<temp .OR. DISNAN(temp) ) value = temp
                  ENDDO
               ENDIF
!              n is even & A is k=n/2 by n+1
            ELSEIF ( ilu==0 ) THEN
!                 uplo = 'U'
               DO i = k , N - 1
                  Work(i) = ZERO
               ENDDO
               DO j = 0 , k - 1
                  s = ZERO
                  DO i = 0 , k - 1
                     aa = ABS(A(i+j*lda))
!                       A(j,i+k)
                     Work(i+k) = Work(i+k) + aa
                     s = s + aa
                  ENDDO
                  Work(j) = s
               ENDDO
!                 j=k
               aa = ABS(DBLE(A(0+j*lda)))
!                 A(k,k)
               s = aa
               DO i = 1 , k - 1
                  aa = ABS(A(i+j*lda))
!                    A(k,k+i)
                  Work(i+k) = Work(i+k) + aa
                  s = s + aa
               ENDDO
               Work(j) = Work(j) + s
               DO j = k + 1 , N - 1
                  s = ZERO
                  DO i = 0 , j - 2 - k
                     aa = ABS(A(i+j*lda))
!                       A(i,j-k-1)
                     Work(i) = Work(i) + aa
                     s = s + aa
                  ENDDO
!                    i=j-1-k
                  aa = ABS(DBLE(A(i+j*lda)))
!                    A(j-k-1,j-k-1)
                  s = s + aa
                  Work(j-k-1) = Work(j-k-1) + s
                  i = i + 1
                  aa = ABS(DBLE(A(i+j*lda)))
!                    A(j,j)
                  s = aa
                  DO l = j + 1 , N - 1
                     i = i + 1
                     aa = ABS(A(i+j*lda))
!                       A(j,l)
                     Work(l) = Work(l) + aa
                     s = s + aa
                  ENDDO
                  Work(j) = Work(j) + s
               ENDDO
!                 j=n
               s = ZERO
               DO i = 0 , k - 2
                  aa = ABS(A(i+j*lda))
!                    A(i,k-1)
                  Work(i) = Work(i) + aa
                  s = s + aa
               ENDDO
!                 i=k-1
               aa = ABS(DBLE(A(i+j*lda)))
!                 A(k-1,k-1)
               s = s + aa
               Work(i) = Work(i) + s
               value = Work(0)
               DO i = 1 , N - 1
                  temp = Work(i)
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
            ELSE
!                 ilu=1 & uplo = 'L'
               DO i = k , N - 1
                  Work(i) = ZERO
               ENDDO
!                 j=0 is special :process col A(k:n-1,k)
               s = ABS(DBLE(A(0)))
!                 A(k,k)
               DO i = 1 , k - 1
                  aa = ABS(A(i))
!                    A(k+i,k)
                  Work(i+k) = Work(i+k) + aa
                  s = s + aa
               ENDDO
               Work(k) = Work(k) + s
               DO j = 1 , k - 1
!                    process
                  s = ZERO
                  DO i = 0 , j - 2
                     aa = ABS(A(i+j*lda))
!                       A(j-1,i)
                     Work(i) = Work(i) + aa
                     s = s + aa
                  ENDDO
                  aa = ABS(DBLE(A(i+j*lda)))
!                    i=j-1 so process of A(j-1,j-1)
                  s = s + aa
                  Work(j-1) = s
!                    is initialised here
                  i = i + 1
!                    i=j process A(j+k,j+k)
                  aa = ABS(DBLE(A(i+j*lda)))
                  s = aa
                  DO l = k + j + 1 , N - 1
                     i = i + 1
                     aa = ABS(A(i+j*lda))
!                       A(l,k+j)
                     s = s + aa
                     Work(l) = Work(l) + aa
                  ENDDO
                  Work(k+j) = Work(k+j) + s
               ENDDO
!                 j=k is special :process col A(k,0:k-1)
               s = ZERO
               DO i = 0 , k - 2
                  aa = ABS(A(i+j*lda))
!                    A(k,i)
                  Work(i) = Work(i) + aa
                  s = s + aa
               ENDDO
!
!                 i=k-1
               aa = ABS(DBLE(A(i+j*lda)))
!                 A(k-1,k-1)
               s = s + aa
               Work(i) = s
!                 done with col j=k+1
               DO j = k + 1 , N
!
!                    process col j-1 of A = A(j-1,0:k-1)
                  s = ZERO
                  DO i = 0 , k - 1
                     aa = ABS(A(i+j*lda))
!                       A(j-1,i)
                     Work(i) = Work(i) + aa
                     s = s + aa
                  ENDDO
                  Work(j-1) = Work(j-1) + s
               ENDDO
               value = Work(0)
               DO i = 1 , N - 1
                  temp = Work(i)
                  IF ( value<temp .OR. DISNAN(temp) ) value = temp
               ENDDO
            ENDIF
         ENDIF
      ELSEIF ( (LSAME(Norm,'F')) .OR. (LSAME(Norm,'E')) ) THEN
!
!       Find normF(A).
!
         k = (N+1)/2
         scale = ZERO
         s = ONE
         IF ( noe==1 ) THEN
!           n is odd
            IF ( ifm==1 ) THEN
!              A is normal & A is n by k
               IF ( ilu==0 ) THEN
!                 A is upper
                  DO j = 0 , k - 3
                     CALL ZLASSQ(k-j-2,A(k+j+1+j*lda),1,scale,s)
!                    L at A(k,0)
                  ENDDO
                  DO j = 0 , k - 1
                     CALL ZLASSQ(k+j-1,A(0+j*lda),1,scale,s)
!                    trap U at A(0,0)
                  ENDDO
                  s = s + s
!                 double s for the off diagonal elements
                  l = k - 1
!                 -> U(k,k) at A(k-1,0)
                  DO i = 0 , k - 2
                     aa = DBLE(A(l))
!                    U(k+i,k+i)
                     IF ( aa/=ZERO ) THEN
                        IF ( scale<aa ) THEN
                           s = ONE + s*(scale/aa)**2
                           scale = aa
                        ELSE
                           s = s + (aa/scale)**2
                        ENDIF
                     ENDIF
                     aa = DBLE(A(l+1))
!                    U(i,i)
                     IF ( aa/=ZERO ) THEN
                        IF ( scale<aa ) THEN
                           s = ONE + s*(scale/aa)**2
                           scale = aa
                        ELSE
                           s = s + (aa/scale)**2
                        ENDIF
                     ENDIF
                     l = l + lda + 1
                  ENDDO
                  aa = DBLE(A(l))
!                 U(n-1,n-1)
                  IF ( aa/=ZERO ) THEN
                     IF ( scale<aa ) THEN
                        s = ONE + s*(scale/aa)**2
                        scale = aa
                     ELSE
                        s = s + (aa/scale)**2
                     ENDIF
                  ENDIF
               ELSE
!                 ilu=1 & A is lower
                  DO j = 0 , k - 1
                     CALL ZLASSQ(N-j-1,A(j+1+j*lda),1,scale,s)
!                    trap L at A(0,0)
                  ENDDO
                  DO j = 1 , k - 2
                     CALL ZLASSQ(j,A(0+(1+j)*lda),1,scale,s)
!                    U at A(0,1)
                  ENDDO
                  s = s + s
!                 double s for the off diagonal elements
                  aa = DBLE(A(0))
!                 L(0,0) at A(0,0)
                  IF ( aa/=ZERO ) THEN
                     IF ( scale<aa ) THEN
                        s = ONE + s*(scale/aa)**2
                        scale = aa
                     ELSE
                        s = s + (aa/scale)**2
                     ENDIF
                  ENDIF
                  l = lda
!                 -> L(k,k) at A(0,1)
                  DO i = 1 , k - 1
                     aa = DBLE(A(l))
!                    L(k-1+i,k-1+i)
                     IF ( aa/=ZERO ) THEN
                        IF ( scale<aa ) THEN
                           s = ONE + s*(scale/aa)**2
                           scale = aa
                        ELSE
                           s = s + (aa/scale)**2
                        ENDIF
                     ENDIF
                     aa = DBLE(A(l+1))
!                    L(i,i)
                     IF ( aa/=ZERO ) THEN
                        IF ( scale<aa ) THEN
                           s = ONE + s*(scale/aa)**2
                           scale = aa
                        ELSE
                           s = s + (aa/scale)**2
                        ENDIF
                     ENDIF
                     l = l + lda + 1
                  ENDDO
               ENDIF
!              A is xpose & A is k by n
            ELSEIF ( ilu==0 ) THEN
!                 A**H is upper
               DO j = 1 , k - 2
                  CALL ZLASSQ(j,A(0+(k+j)*lda),1,scale,s)
!                    U at A(0,k)
               ENDDO
               DO j = 0 , k - 2
                  CALL ZLASSQ(k,A(0+j*lda),1,scale,s)
!                    k by k-1 rect. at A(0,0)
               ENDDO
               DO j = 0 , k - 2
                  CALL ZLASSQ(k-j-1,A(j+1+(j+k-1)*lda),1,scale,s)
!                    L at A(0,k-1)
               ENDDO
               s = s + s
!                 double s for the off diagonal elements
               l = 0 + k*lda - lda
!                 -> U(k-1,k-1) at A(0,k-1)
               aa = DBLE(A(l))
!                 U(k-1,k-1)
               IF ( aa/=ZERO ) THEN
                  IF ( scale<aa ) THEN
                     s = ONE + s*(scale/aa)**2
                     scale = aa
                  ELSE
                     s = s + (aa/scale)**2
                  ENDIF
               ENDIF
               l = l + lda
!                 -> U(0,0) at A(0,k)
               DO j = k , N - 1
                  aa = DBLE(A(l))
!                    -> U(j-k,j-k)
                  IF ( aa/=ZERO ) THEN
                     IF ( scale<aa ) THEN
                        s = ONE + s*(scale/aa)**2
                        scale = aa
                     ELSE
                        s = s + (aa/scale)**2
                     ENDIF
                  ENDIF
                  aa = DBLE(A(l+1))
!                    -> U(j,j)
                  IF ( aa/=ZERO ) THEN
                     IF ( scale<aa ) THEN
                        s = ONE + s*(scale/aa)**2
                        scale = aa
                     ELSE
                        s = s + (aa/scale)**2
                     ENDIF
                  ENDIF
                  l = l + lda + 1
               ENDDO
            ELSE
!                 A**H is lower
               DO j = 1 , k - 1
                  CALL ZLASSQ(j,A(0+j*lda),1,scale,s)
!                    U at A(0,0)
               ENDDO
               DO j = k , N - 1
                  CALL ZLASSQ(k,A(0+j*lda),1,scale,s)
!                    k by k-1 rect. at A(0,k)
               ENDDO
               DO j = 0 , k - 3
                  CALL ZLASSQ(k-j-2,A(j+2+j*lda),1,scale,s)
!                    L at A(1,0)
               ENDDO
               s = s + s
!                 double s for the off diagonal elements
               l = 0
!                 -> L(0,0) at A(0,0)
               DO i = 0 , k - 2
                  aa = DBLE(A(l))
!                    L(i,i)
                  IF ( aa/=ZERO ) THEN
                     IF ( scale<aa ) THEN
                        s = ONE + s*(scale/aa)**2
                        scale = aa
                     ELSE
                        s = s + (aa/scale)**2
                     ENDIF
                  ENDIF
                  aa = DBLE(A(l+1))
!                    L(k+i,k+i)
                  IF ( aa/=ZERO ) THEN
                     IF ( scale<aa ) THEN
                        s = ONE + s*(scale/aa)**2
                        scale = aa
                     ELSE
                        s = s + (aa/scale)**2
                     ENDIF
                  ENDIF
                  l = l + lda + 1
               ENDDO
!                 L-> k-1 + (k-1)*lda or L(k-1,k-1) at A(k-1,k-1)
               aa = DBLE(A(l))
!                 L(k-1,k-1) at A(k-1,k-1)
               IF ( aa/=ZERO ) THEN
                  IF ( scale<aa ) THEN
                     s = ONE + s*(scale/aa)**2
                     scale = aa
                  ELSE
                     s = s + (aa/scale)**2
                  ENDIF
               ENDIF
            ENDIF
!           n is even
         ELSEIF ( ifm==1 ) THEN
!              A is normal
            IF ( ilu==0 ) THEN
!                 A is upper
               DO j = 0 , k - 2
                  CALL ZLASSQ(k-j-1,A(k+j+2+j*lda),1,scale,s)
!                 L at A(k+1,0)
               ENDDO
               DO j = 0 , k - 1
                  CALL ZLASSQ(k+j,A(0+j*lda),1,scale,s)
!                 trap U at A(0,0)
               ENDDO
               s = s + s
!                 double s for the off diagonal elements
               l = k
!                 -> U(k,k) at A(k,0)
               DO i = 0 , k - 1
                  aa = DBLE(A(l))
!                    U(k+i,k+i)
                  IF ( aa/=ZERO ) THEN
                     IF ( scale<aa ) THEN
                        s = ONE + s*(scale/aa)**2
                        scale = aa
                     ELSE
                        s = s + (aa/scale)**2
                     ENDIF
                  ENDIF
                  aa = DBLE(A(l+1))
!                    U(i,i)
                  IF ( aa/=ZERO ) THEN
                     IF ( scale<aa ) THEN
                        s = ONE + s*(scale/aa)**2
                        scale = aa
                     ELSE
                        s = s + (aa/scale)**2
                     ENDIF
                  ENDIF
                  l = l + lda + 1
               ENDDO
            ELSE
!                 ilu=1 & A is lower
               DO j = 0 , k - 1
                  CALL ZLASSQ(N-j-1,A(j+2+j*lda),1,scale,s)
!                    trap L at A(1,0)
               ENDDO
               DO j = 1 , k - 1
                  CALL ZLASSQ(j,A(0+j*lda),1,scale,s)
!                    U at A(0,0)
               ENDDO
               s = s + s
!                 double s for the off diagonal elements
               l = 0
!                 -> L(k,k) at A(0,0)
               DO i = 0 , k - 1
                  aa = DBLE(A(l))
!                    L(k-1+i,k-1+i)
                  IF ( aa/=ZERO ) THEN
                     IF ( scale<aa ) THEN
                        s = ONE + s*(scale/aa)**2
                        scale = aa
                     ELSE
                        s = s + (aa/scale)**2
                     ENDIF
                  ENDIF
                  aa = DBLE(A(l+1))
!                    L(i,i)
                  IF ( aa/=ZERO ) THEN
                     IF ( scale<aa ) THEN
                        s = ONE + s*(scale/aa)**2
                        scale = aa
                     ELSE
                        s = s + (aa/scale)**2
                     ENDIF
                  ENDIF
                  l = l + lda + 1
               ENDDO
            ENDIF
!              A is xpose
         ELSEIF ( ilu==0 ) THEN
!                 A**H is upper
            DO j = 1 , k - 1
               CALL ZLASSQ(j,A(0+(k+1+j)*lda),1,scale,s)
!                 U at A(0,k+1)
            ENDDO
            DO j = 0 , k - 1
               CALL ZLASSQ(k,A(0+j*lda),1,scale,s)
!                 k by k rect. at A(0,0)
            ENDDO
            DO j = 0 , k - 2
               CALL ZLASSQ(k-j-1,A(j+1+(j+k)*lda),1,scale,s)
!                 L at A(0,k)
            ENDDO
            s = s + s
!                 double s for the off diagonal elements
            l = 0 + k*lda
!                 -> U(k,k) at A(0,k)
            aa = DBLE(A(l))
!                 U(k,k)
            IF ( aa/=ZERO ) THEN
               IF ( scale<aa ) THEN
                  s = ONE + s*(scale/aa)**2
                  scale = aa
               ELSE
                  s = s + (aa/scale)**2
               ENDIF
            ENDIF
            l = l + lda
!                 -> U(0,0) at A(0,k+1)
            DO j = k + 1 , N - 1
               aa = DBLE(A(l))
!                    -> U(j-k-1,j-k-1)
               IF ( aa/=ZERO ) THEN
                  IF ( scale<aa ) THEN
                     s = ONE + s*(scale/aa)**2
                     scale = aa
                  ELSE
                     s = s + (aa/scale)**2
                  ENDIF
               ENDIF
               aa = DBLE(A(l+1))
!                    -> U(j,j)
               IF ( aa/=ZERO ) THEN
                  IF ( scale<aa ) THEN
                     s = ONE + s*(scale/aa)**2
                     scale = aa
                  ELSE
                     s = s + (aa/scale)**2
                  ENDIF
               ENDIF
               l = l + lda + 1
            ENDDO
!                 L=k-1+n*lda
!                 -> U(k-1,k-1) at A(k-1,n)
            aa = DBLE(A(l))
!                 U(k,k)
            IF ( aa/=ZERO ) THEN
               IF ( scale<aa ) THEN
                  s = ONE + s*(scale/aa)**2
                  scale = aa
               ELSE
                  s = s + (aa/scale)**2
               ENDIF
            ENDIF
         ELSE
!                 A**H is lower
            DO j = 1 , k - 1
               CALL ZLASSQ(j,A(0+(j+1)*lda),1,scale,s)
!                 U at A(0,1)
            ENDDO
            DO j = k + 1 , N
               CALL ZLASSQ(k,A(0+j*lda),1,scale,s)
!                 k by k rect. at A(0,k+1)
            ENDDO
            DO j = 0 , k - 2
               CALL ZLASSQ(k-j-1,A(j+1+j*lda),1,scale,s)
!                 L at A(0,0)
            ENDDO
            s = s + s
!                 double s for the off diagonal elements
            l = 0
!                 -> L(k,k) at A(0,0)
            aa = DBLE(A(l))
!                 L(k,k) at A(0,0)
            IF ( aa/=ZERO ) THEN
               IF ( scale<aa ) THEN
                  s = ONE + s*(scale/aa)**2
                  scale = aa
               ELSE
                  s = s + (aa/scale)**2
               ENDIF
            ENDIF
            l = lda
!                 -> L(0,0) at A(0,1)
            DO i = 0 , k - 2
               aa = DBLE(A(l))
!                    L(i,i)
               IF ( aa/=ZERO ) THEN
                  IF ( scale<aa ) THEN
                     s = ONE + s*(scale/aa)**2
                     scale = aa
                  ELSE
                     s = s + (aa/scale)**2
                  ENDIF
               ENDIF
               aa = DBLE(A(l+1))
!                    L(k+i+1,k+i+1)
               IF ( aa/=ZERO ) THEN
                  IF ( scale<aa ) THEN
                     s = ONE + s*(scale/aa)**2
                     scale = aa
                  ELSE
                     s = s + (aa/scale)**2
                  ENDIF
               ENDIF
               l = l + lda + 1
            ENDDO
!                 L-> k - 1 + k*lda or L(k-1,k-1) at A(k-1,k)
            aa = DBLE(A(l))
!                 L(k-1,k-1) at A(k-1,k)
            IF ( aa/=ZERO ) THEN
               IF ( scale<aa ) THEN
                  s = ONE + s*(scale/aa)**2
                  scale = aa
               ELSE
                  s = s + (aa/scale)**2
               ENDIF
            ENDIF
         ENDIF
         value = scale*SQRT(s)
      ENDIF
!
      ZLANHF = value
!
!     End of ZLANHF
!
      END FUNCTION ZLANHF
