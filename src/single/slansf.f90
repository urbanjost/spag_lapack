!*==slansf.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLANSF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLANSF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slansf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slansf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slansf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SLANSF( NORM, TRANSR, UPLO, N, A, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM, TRANSR, UPLO
!       INTEGER            N
!       ..
!       .. Array Arguments ..
!       REAL               A( 0: * ), WORK( 0: * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLANSF returns the value of the one norm, or the Frobenius norm, or
!> the infinity norm, or the element of largest absolute value of a
!> real symmetric matrix A in RFP format.
!> \endverbatim
!>
!> \return SLANSF
!> \verbatim
!>
!>    SLANSF = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          NORM is CHARACTER*1
!>          Specifies the value to be returned in SLANSF as described
!>          above.
!> \endverbatim
!>
!> \param[in] TRANSR
!> \verbatim
!>          TRANSR is CHARACTER*1
!>          Specifies whether the RFP format of A is normal or
!>          transposed format.
!>          = 'N':  RFP format is Normal;
!>          = 'T':  RFP format is Transpose.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the RFP matrix A came from
!>           an upper or lower triangular matrix as follows:
!>           = 'U': RFP A came from an upper triangular matrix;
!>           = 'L': RFP A came from a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A. N >= 0. When N = 0, SLANSF is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension ( N*(N+1)/2 );
!>          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L')
!>          part of the symmetric matrix A stored in RFP format. See the
!>          "Notes" below for more details.
!>          Unchanged on exit.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK)),
!>          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!>          WORK is not referenced.
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
!> \ingroup realOTHERcomputational
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
      REAL FUNCTION SLANSF(Norm,Transr,Uplo,N,A,Work)
      IMPLICIT NONE
!*--SLANSF213
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Norm , Transr , Uplo
      INTEGER N
!     ..
!     .. Array Arguments ..
      REAL A(0:*) , Work(0:*)
!     ..
!
!  =====================================================================
!
!     ..
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j , ifm , ilu , noe , n1 , k , l , lda
      REAL scale , s , value , aa , temp
!     ..
!     .. External Functions ..
      LOGICAL LSAME , SISNAN
      EXTERNAL LSAME , SISNAN
!     ..
!     .. External Subroutines ..
      EXTERNAL SLASSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SQRT
!     ..
!     .. Executable Statements ..
!
      IF ( N==0 ) THEN
         SLANSF = ZERO
         RETURN
      ELSEIF ( N==1 ) THEN
         SLANSF = ABS(A(0))
         RETURN
      ENDIF
!
!     set noe = 1 if n is odd. if n is even set noe=0
!
      noe = 1
      IF ( MOD(N,2)==0 ) noe = 0
!
!     set ifm = 0 when form='T or 't' and 1 otherwise
!
      ifm = 1
      IF ( LSAME(Transr,'T') ) ifm = 0
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
!           n is odd
            IF ( ifm==1 ) THEN
!           A is n by k
               DO j = 0 , k - 1
                  DO i = 0 , N - 1
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. SISNAN(temp) ) value = temp
                  ENDDO
               ENDDO
            ELSE
!              xpose case; A is k by n
               DO j = 0 , N - 1
                  DO i = 0 , k - 1
                     temp = ABS(A(i+j*lda))
                     IF ( value<temp .OR. SISNAN(temp) ) value = temp
                  ENDDO
               ENDDO
            ENDIF
!           n is even
         ELSEIF ( ifm==1 ) THEN
!              A is n+1 by k
            DO j = 0 , k - 1
               DO i = 0 , N
                  temp = ABS(A(i+j*lda))
                  IF ( value<temp .OR. SISNAN(temp) ) value = temp
               ENDDO
            ENDDO
         ELSE
!              xpose case; A is k by n+1
            DO j = 0 , N
               DO i = 0 , k - 1
                  temp = ABS(A(i+j*lda))
                  IF ( value<temp .OR. SISNAN(temp) ) value = temp
               ENDDO
            ENDDO
         ENDIF
      ELSEIF ( (LSAME(Norm,'I')) .OR. (LSAME(Norm,'O')) .OR. (Norm=='1')&
     &         ) THEN
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
         IF ( ifm==1 ) THEN
            k = N/2
            IF ( noe==1 ) THEN
!              n is odd
               IF ( ilu==0 ) THEN
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
                     aa = ABS(A(i+j*lda))
!                    -> A(j+k,j+k)
                     Work(j+k) = s + aa
                     IF ( i==k+k ) EXIT
                     i = i + 1
                     aa = ABS(A(i+j*lda))
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
                     IF ( value<temp .OR. SISNAN(temp) ) value = temp
                  ENDDO
               ELSE
!                 ilu = 1
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
                        aa = ABS(A(i+j*lda))
!                       -> A(j+k,j+k)
                        s = s + aa
                        Work(i+k) = Work(i+k) + s
!                       i=j
                        i = i + 1
                     ENDIF
                     aa = ABS(A(i+j*lda))
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
                     IF ( value<temp .OR. SISNAN(temp) ) value = temp
                  ENDDO
               ENDIF
!              n is even
            ELSEIF ( ilu==0 ) THEN
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
                  aa = ABS(A(i+j*lda))
!                    -> A(j+k,j+k)
                  Work(j+k) = s + aa
                  i = i + 1
                  aa = ABS(A(i+j*lda))
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
                  IF ( value<temp .OR. SISNAN(temp) ) value = temp
               ENDDO
            ELSE
!                 ilu = 1
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
                  aa = ABS(A(i+j*lda))
!                    -> A(j+k,j+k)
                  s = s + aa
                  Work(i+k) = Work(i+k) + s
!                    i=j
                  i = i + 1
                  aa = ABS(A(i+j*lda))
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
                  IF ( value<temp .OR. SISNAN(temp) ) value = temp
               ENDDO
            ENDIF
         ELSE
!           ifm=0
            k = N/2
            IF ( noe==1 ) THEN
!              n is odd
               IF ( ilu==0 ) THEN
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
                  s = ABS(A(0+j*lda))
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
                     aa = ABS(A(i+j*lda))
!                    A(j-k,j-k)
                     s = s + aa
                     Work(j-k) = Work(j-k) + s
                     i = i + 1
                     s = ABS(A(i+j*lda))
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
                     IF ( value<temp .OR. SISNAN(temp) ) value = temp
                  ENDDO
               ELSE
!                 ilu=1
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
                     aa = ABS(A(i+j*lda))
!                    i=j so process of A(j,j)
                     s = s + aa
                     Work(j) = s
!                    is initialised here
                     i = i + 1
!                    i=j process A(j+k,j+k)
                     aa = ABS(A(i+j*lda))
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
                  aa = ABS(A(i+j*lda))
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
                     IF ( value<temp .OR. SISNAN(temp) ) value = temp
                  ENDDO
               ENDIF
!              n is even
            ELSEIF ( ilu==0 ) THEN
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
               aa = ABS(A(0+j*lda))
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
!                     i=j-1-k
                  aa = ABS(A(i+j*lda))
!                    A(j-k-1,j-k-1)
                  s = s + aa
                  Work(j-k-1) = Work(j-k-1) + s
                  i = i + 1
                  aa = ABS(A(i+j*lda))
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
               aa = ABS(A(i+j*lda))
!                 A(k-1,k-1)
               s = s + aa
               Work(i) = Work(i) + s
               value = Work(0)
               DO i = 1 , N - 1
                  temp = Work(i)
                  IF ( value<temp .OR. SISNAN(temp) ) value = temp
               ENDDO
            ELSE
!                 ilu=1
               DO i = k , N - 1
                  Work(i) = ZERO
               ENDDO
!                 j=0 is special :process col A(k:n-1,k)
               s = ABS(A(0))
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
                  aa = ABS(A(i+j*lda))
!                    i=j-1 so process of A(j-1,j-1)
                  s = s + aa
                  Work(j-1) = s
!                    is initialised here
                  i = i + 1
!                    i=j process A(j+k,j+k)
                  aa = ABS(A(i+j*lda))
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
!                 i=k-1
               aa = ABS(A(i+j*lda))
!                 A(k-1,k-1)
               s = s + aa
               Work(i) = s
!                 done with col j=k+1
               DO j = k + 1 , N
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
                  IF ( value<temp .OR. SISNAN(temp) ) value = temp
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
!              A is normal
               IF ( ilu==0 ) THEN
!                 A is upper
                  DO j = 0 , k - 3
                     CALL SLASSQ(k-j-2,A(k+j+1+j*lda),1,scale,s)
!                    L at A(k,0)
                  ENDDO
                  DO j = 0 , k - 1
                     CALL SLASSQ(k+j-1,A(0+j*lda),1,scale,s)
!                    trap U at A(0,0)
                  ENDDO
                  s = s + s
!                 double s for the off diagonal elements
                  CALL SLASSQ(k-1,A(k),lda+1,scale,s)
!                 tri L at A(k,0)
                  CALL SLASSQ(k,A(k-1),lda+1,scale,s)
!                 tri U at A(k-1,0)
               ELSE
!                 ilu=1 & A is lower
                  DO j = 0 , k - 1
                     CALL SLASSQ(N-j-1,A(j+1+j*lda),1,scale,s)
!                    trap L at A(0,0)
                  ENDDO
                  DO j = 0 , k - 2
                     CALL SLASSQ(j,A(0+(1+j)*lda),1,scale,s)
!                    U at A(0,1)
                  ENDDO
                  s = s + s
!                 double s for the off diagonal elements
                  CALL SLASSQ(k,A(0),lda+1,scale,s)
!                 tri L at A(0,0)
                  CALL SLASSQ(k-1,A(0+lda),lda+1,scale,s)
!                 tri U at A(0,1)
               ENDIF
!              A is xpose
            ELSEIF ( ilu==0 ) THEN
!                 A**T is upper
               DO j = 1 , k - 2
                  CALL SLASSQ(j,A(0+(k+j)*lda),1,scale,s)
!                    U at A(0,k)
               ENDDO
               DO j = 0 , k - 2
                  CALL SLASSQ(k,A(0+j*lda),1,scale,s)
!                    k by k-1 rect. at A(0,0)
               ENDDO
               DO j = 0 , k - 2
                  CALL SLASSQ(k-j-1,A(j+1+(j+k-1)*lda),1,scale,s)
!                    L at A(0,k-1)
               ENDDO
               s = s + s
!                 double s for the off diagonal elements
               CALL SLASSQ(k-1,A(0+k*lda),lda+1,scale,s)
!                 tri U at A(0,k)
               CALL SLASSQ(k,A(0+(k-1)*lda),lda+1,scale,s)
!                 tri L at A(0,k-1)
            ELSE
!                 A**T is lower
               DO j = 1 , k - 1
                  CALL SLASSQ(j,A(0+j*lda),1,scale,s)
!                    U at A(0,0)
               ENDDO
               DO j = k , N - 1
                  CALL SLASSQ(k,A(0+j*lda),1,scale,s)
!                    k by k-1 rect. at A(0,k)
               ENDDO
               DO j = 0 , k - 3
                  CALL SLASSQ(k-j-2,A(j+2+j*lda),1,scale,s)
!                    L at A(1,0)
               ENDDO
               s = s + s
!                 double s for the off diagonal elements
               CALL SLASSQ(k,A(0),lda+1,scale,s)
!                 tri U at A(0,0)
               CALL SLASSQ(k-1,A(1),lda+1,scale,s)
!                 tri L at A(1,0)
            ENDIF
!           n is even
         ELSEIF ( ifm==1 ) THEN
!              A is normal
            IF ( ilu==0 ) THEN
!                 A is upper
               DO j = 0 , k - 2
                  CALL SLASSQ(k-j-1,A(k+j+2+j*lda),1,scale,s)
!                    L at A(k+1,0)
               ENDDO
               DO j = 0 , k - 1
                  CALL SLASSQ(k+j,A(0+j*lda),1,scale,s)
!                    trap U at A(0,0)
               ENDDO
               s = s + s
!                 double s for the off diagonal elements
               CALL SLASSQ(k,A(k+1),lda+1,scale,s)
!                 tri L at A(k+1,0)
               CALL SLASSQ(k,A(k),lda+1,scale,s)
!                 tri U at A(k,0)
            ELSE
!                 ilu=1 & A is lower
               DO j = 0 , k - 1
                  CALL SLASSQ(N-j-1,A(j+2+j*lda),1,scale,s)
!                    trap L at A(1,0)
               ENDDO
               DO j = 1 , k - 1
                  CALL SLASSQ(j,A(0+j*lda),1,scale,s)
!                    U at A(0,0)
               ENDDO
               s = s + s
!                 double s for the off diagonal elements
               CALL SLASSQ(k,A(1),lda+1,scale,s)
!                 tri L at A(1,0)
               CALL SLASSQ(k,A(0),lda+1,scale,s)
!                 tri U at A(0,0)
            ENDIF
!              A is xpose
         ELSEIF ( ilu==0 ) THEN
!                 A**T is upper
            DO j = 1 , k - 1
               CALL SLASSQ(j,A(0+(k+1+j)*lda),1,scale,s)
!                    U at A(0,k+1)
            ENDDO
            DO j = 0 , k - 1
               CALL SLASSQ(k,A(0+j*lda),1,scale,s)
!                    k by k rect. at A(0,0)
            ENDDO
            DO j = 0 , k - 2
               CALL SLASSQ(k-j-1,A(j+1+(j+k)*lda),1,scale,s)
!                    L at A(0,k)
            ENDDO
            s = s + s
!                 double s for the off diagonal elements
            CALL SLASSQ(k,A(0+(k+1)*lda),lda+1,scale,s)
!                 tri U at A(0,k+1)
            CALL SLASSQ(k,A(0+k*lda),lda+1,scale,s)
!                 tri L at A(0,k)
         ELSE
!                 A**T is lower
            DO j = 1 , k - 1
               CALL SLASSQ(j,A(0+(j+1)*lda),1,scale,s)
!                    U at A(0,1)
            ENDDO
            DO j = k + 1 , N
               CALL SLASSQ(k,A(0+j*lda),1,scale,s)
!                    k by k rect. at A(0,k+1)
            ENDDO
            DO j = 0 , k - 2
               CALL SLASSQ(k-j-1,A(j+1+j*lda),1,scale,s)
!                    L at A(0,0)
            ENDDO
            s = s + s
!                 double s for the off diagonal elements
            CALL SLASSQ(k,A(lda),lda+1,scale,s)
!                 tri L at A(0,1)
            CALL SLASSQ(k,A(0),lda+1,scale,s)
!                 tri U at A(0,0)
         ENDIF
         value = scale*SQRT(s)
      ENDIF
!
      SLANSF = value
!
!     End of SLANSF
!
      END FUNCTION SLANSF
