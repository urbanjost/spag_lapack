!*==stfttr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b STFTTR copies a triangular matrix from the rectangular full packed format (TF) to the standard full format (TR).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STFTTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stfttr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stfttr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stfttr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STFTTR( TRANSR, UPLO, N, ARF, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANSR, UPLO
!       INTEGER            INFO, N, LDA
!       ..
!       .. Array Arguments ..
!       REAL               A( 0: LDA-1, 0: * ), ARF( 0: * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STFTTR copies a triangular matrix A from rectangular full packed
!> format (TF) to standard full format (TR).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSR
!> \verbatim
!>          TRANSR is CHARACTER*1
!>          = 'N':  ARF is in Normal format;
!>          = 'T':  ARF is in Transpose format.
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
!>          The order of the matrices ARF and A. N >= 0.
!> \endverbatim
!>
!> \param[in] ARF
!> \verbatim
!>          ARF is REAL array, dimension (N*(N+1)/2).
!>          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L')
!>          matrix A in RFP format. See the "Notes" below for more
!>          details.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On exit, the triangular matrix A.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of the array A contains
!>          the upper triangular matrix, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of the array A contains
!>          the lower triangular matrix, and the strictly upper
!>          triangular part of A is not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
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
      SUBROUTINE STFTTR(Transr,Uplo,N,Arf,A,Lda,Info)
      IMPLICIT NONE
!*--STFTTR200
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Transr , Uplo
      INTEGER Info , N , Lda
!     ..
!     .. Array Arguments ..
      REAL A(0:Lda-1,0:*) , Arf(0:*)
!     ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
      LOGICAL lower , nisodd , normaltransr
      INTEGER n1 , n2 , k , nt , nx2 , np1x2
      INTEGER i , j , l , ij
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MOD
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
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('STFTTR',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=1 ) THEN
         IF ( N==1 ) A(0,0) = Arf(0)
         RETURN
      ENDIF
!
!     Size of array ARF(0:nt-1)
!
      nt = N*(N+1)/2
!
!     set N1 and N2 depending on LOWER: for N even N1=N2=K
!
      IF ( lower ) THEN
         n2 = N/2
         n1 = N - n2
      ELSE
         n1 = N/2
         n2 = N - n1
      ENDIF
!
!     If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2.
!     If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is
!     N--by--(N+1)/2.
!
      IF ( MOD(N,2)==0 ) THEN
         k = N/2
         nisodd = .FALSE.
         IF ( .NOT.lower ) np1x2 = N + N + 2
      ELSE
         nisodd = .TRUE.
         IF ( .NOT.lower ) nx2 = N + N
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
               ij = 0
               DO j = 0 , n2
                  DO i = n1 , n2 + j
                     A(n2+j,i) = Arf(ij)
                     ij = ij + 1
                  ENDDO
                  DO i = j , N - 1
                     A(i,j) = Arf(ij)
                     ij = ij + 1
                  ENDDO
               ENDDO
!
            ELSE
!
!              N is odd, TRANSR = 'N', and UPLO = 'U'
!
               ij = nt - N
               DO j = N - 1 , n1 , -1
                  DO i = 0 , j
                     A(i,j) = Arf(ij)
                     ij = ij + 1
                  ENDDO
                  DO l = j - n1 , n1 - 1
                     A(j-n1,l) = Arf(ij)
                     ij = ij + 1
                  ENDDO
                  ij = ij - nx2
               ENDDO
!
            ENDIF
!
!
!           N is odd and TRANSR = 'T'
!
         ELSEIF ( lower ) THEN
!
!              N is odd, TRANSR = 'T', and UPLO = 'L'
!
            ij = 0
            DO j = 0 , n2 - 1
               DO i = 0 , j
                  A(j,i) = Arf(ij)
                  ij = ij + 1
               ENDDO
               DO i = n1 + j , N - 1
                  A(i,n1+j) = Arf(ij)
                  ij = ij + 1
               ENDDO
            ENDDO
            DO j = n2 , N - 1
               DO i = 0 , n1 - 1
                  A(j,i) = Arf(ij)
                  ij = ij + 1
               ENDDO
            ENDDO
!
         ELSE
!
!              N is odd, TRANSR = 'T', and UPLO = 'U'
!
            ij = 0
            DO j = 0 , n1
               DO i = n1 , N - 1
                  A(j,i) = Arf(ij)
                  ij = ij + 1
               ENDDO
            ENDDO
            DO j = 0 , n1 - 1
               DO i = 0 , j
                  A(i,j) = Arf(ij)
                  ij = ij + 1
               ENDDO
               DO l = n2 + j , N - 1
                  A(n2+j,l) = Arf(ij)
                  ij = ij + 1
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
!              N is even, TRANSR = 'N', and UPLO = 'L'
!
            ij = 0
            DO j = 0 , k - 1
               DO i = k , k + j
                  A(k+j,i) = Arf(ij)
                  ij = ij + 1
               ENDDO
               DO i = j , N - 1
                  A(i,j) = Arf(ij)
                  ij = ij + 1
               ENDDO
            ENDDO
!
         ELSE
!
!              N is even, TRANSR = 'N', and UPLO = 'U'
!
            ij = nt - N - 1
            DO j = N - 1 , k , -1
               DO i = 0 , j
                  A(i,j) = Arf(ij)
                  ij = ij + 1
               ENDDO
               DO l = j - k , k - 1
                  A(j-k,l) = Arf(ij)
                  ij = ij + 1
               ENDDO
               ij = ij - np1x2
            ENDDO
!
         ENDIF
!
!
!           N is even and TRANSR = 'T'
!
      ELSEIF ( lower ) THEN
!
!              N is even, TRANSR = 'T', and UPLO = 'L'
!
         ij = 0
         j = k
         DO i = k , N - 1
            A(i,j) = Arf(ij)
            ij = ij + 1
         ENDDO
         DO j = 0 , k - 2
            DO i = 0 , j
               A(j,i) = Arf(ij)
               ij = ij + 1
            ENDDO
            DO i = k + 1 + j , N - 1
               A(i,k+1+j) = Arf(ij)
               ij = ij + 1
            ENDDO
         ENDDO
         DO j = k - 1 , N - 1
            DO i = 0 , k - 1
               A(j,i) = Arf(ij)
               ij = ij + 1
            ENDDO
         ENDDO
!
      ELSE
!
!              N is even, TRANSR = 'T', and UPLO = 'U'
!
         ij = 0
         DO j = 0 , k
            DO i = k , N - 1
               A(j,i) = Arf(ij)
               ij = ij + 1
            ENDDO
         ENDDO
         DO j = 0 , k - 2
            DO i = 0 , j
               A(i,j) = Arf(ij)
               ij = ij + 1
            ENDDO
            DO l = k + 1 + j , N - 1
               A(k+1+j,l) = Arf(ij)
               ij = ij + 1
            ENDDO
         ENDDO
!              Note that here, on exit of the loop, J = K-1
         DO i = 0 , j
            A(i,j) = Arf(ij)
            ij = ij + 1
         ENDDO
!
!
!
      ENDIF
!
!
!     End of STFTTR
!
      END SUBROUTINE STFTTR
