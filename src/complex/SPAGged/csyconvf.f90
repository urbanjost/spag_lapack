!*==csyconvf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CSYCONVF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSYCONVF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyconvf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyconvf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyconvf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO, WAY
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> If parameter WAY = 'C':
!> CSYCONVF converts the factorization output format used in
!> CSYTRF provided on entry in parameter A into the factorization
!> output format used in CSYTRF_RK (or CSYTRF_BK) that is stored
!> on exit in parameters A and E. It also converts in place details of
!> the intechanges stored in IPIV from the format used in CSYTRF into
!> the format used in CSYTRF_RK (or CSYTRF_BK).
!>
!> If parameter WAY = 'R':
!> CSYCONVF performs the conversion in reverse direction, i.e.
!> converts the factorization output format used in CSYTRF_RK
!> (or CSYTRF_BK) provided on entry in parameters A and E into
!> the factorization output format used in CSYTRF that is stored
!> on exit in parameter A. It also converts in place details of
!> the intechanges stored in IPIV from the format used in CSYTRF_RK
!> (or CSYTRF_BK) into the format used in CSYTRF.
!>
!> CSYCONVF can also convert in Hermitian matrix case, i.e. between
!> formats used in CHETRF and CHETRF_RK (or CHETRF_BK).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are
!>          stored as an upper or lower triangular matrix A.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] WAY
!> \verbatim
!>          WAY is CHARACTER*1
!>          = 'C': Convert
!>          = 'R': Revert
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
!>          A is COMPLEX array, dimension (LDA,N)
!>
!>          1) If WAY ='C':
!>
!>          On entry, contains factorization details in format used in
!>          CSYTRF:
!>            a) all elements of the symmetric block diagonal
!>               matrix D on the diagonal of A and on superdiagonal
!>               (or subdiagonal) of A, and
!>            b) If UPLO = 'U': multipliers used to obtain factor U
!>               in the superdiagonal part of A.
!>               If UPLO = 'L': multipliers used to obtain factor L
!>               in the superdiagonal part of A.
!>
!>          On exit, contains factorization details in format used in
!>          CSYTRF_RK or CSYTRF_BK:
!>            a) ONLY diagonal elements of the symmetric block diagonal
!>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
!>               (superdiagonal (or subdiagonal) elements of D
!>                are stored on exit in array E), and
!>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
!>               If UPLO = 'L': factor L in the subdiagonal part of A.
!>
!>          2) If WAY = 'R':
!>
!>          On entry, contains factorization details in format used in
!>          CSYTRF_RK or CSYTRF_BK:
!>            a) ONLY diagonal elements of the symmetric block diagonal
!>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
!>               (superdiagonal (or subdiagonal) elements of D
!>                are stored on exit in array E), and
!>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
!>               If UPLO = 'L': factor L in the subdiagonal part of A.
!>
!>          On exit, contains factorization details in format used in
!>          CSYTRF:
!>            a) all elements of the symmetric block diagonal
!>               matrix D on the diagonal of A and on superdiagonal
!>               (or subdiagonal) of A, and
!>            b) If UPLO = 'U': multipliers used to obtain factor U
!>               in the superdiagonal part of A.
!>               If UPLO = 'L': multipliers used to obtain factor L
!>               in the superdiagonal part of A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is COMPLEX array, dimension (N)
!>
!>          1) If WAY ='C':
!>
!>          On entry, just a workspace.
!>
!>          On exit, contains the superdiagonal (or subdiagonal)
!>          elements of the symmetric block diagonal matrix D
!>          with 1-by-1 or 2-by-2 diagonal blocks, where
!>          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
!>          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
!>
!>          2) If WAY = 'R':
!>
!>          On entry, contains the superdiagonal (or subdiagonal)
!>          elements of the symmetric block diagonal matrix D
!>          with 1-by-1 or 2-by-2 diagonal blocks, where
!>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
!>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
!>
!>          On exit, is not changed
!> \endverbatim
!.
!> \param[in,out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>
!>          1) If WAY ='C':
!>          On entry, details of the interchanges and the block
!>          structure of D in the format used in CSYTRF.
!>          On exit, details of the interchanges and the block
!>          structure of D in the format used in CSYTRF_RK
!>          ( or CSYTRF_BK).
!>
!>          1) If WAY ='R':
!>          On entry, details of the interchanges and the block
!>          structure of D in the format used in CSYTRF_RK
!>          ( or CSYTRF_BK).
!>          On exit, details of the interchanges and the block
!>          structure of D in the format used in CSYTRF.
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
!> \date November 2017
!
!> \ingroup complexSYcomputational
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  November 2017,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!> \endverbatim
!  =====================================================================
      SUBROUTINE CSYCONVF(Uplo,Way,N,A,Lda,E,Ipiv,Info)
      USE S_CSWAP
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CSYCONVF216
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      CHARACTER :: Way
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: E
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: convert , upper
      INTEGER :: i , ip
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
!     .. External Functions ..
!
!     .. External Subroutines ..
!     .. Local Scalars ..
!     ..
!     .. Executable Statements ..
!
      Info = 0
      upper = LSAME(Uplo,'U')
      convert = LSAME(Way,'C')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( .NOT.convert .AND. .NOT.LSAME(Way,'R') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
 
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CSYCONVF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Begin A is UPPER
!
         IF ( convert ) THEN
!
!           Convert A (A is upper)
!
!
!           Convert VALUE
!
!           Assign superdiagonal entries of D to array E and zero out
!           corresponding entries in input storage A
!
            i = N
            E(1) = ZERO
            DO WHILE ( i>1 )
               IF ( Ipiv(i)<0 ) THEN
                  E(i) = A(i-1,i)
                  E(i-1) = ZERO
                  A(i-1,i) = ZERO
                  i = i - 1
               ELSE
                  E(i) = ZERO
               ENDIF
               i = i - 1
            ENDDO
!
!           Convert PERMUTATIONS and IPIV
!
!           Apply permutations to submatrices of upper part of A
!           in factorization order where i decreases from N to 1
!
            i = N
            DO WHILE ( i>=1 )
               IF ( Ipiv(i)>0 ) THEN
!
!                 1-by-1 pivot interchange
!
!                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
!
                  ip = Ipiv(i)
                  IF ( i<N ) THEN
                     IF ( ip/=i ) CALL CSWAP(N-i,A(i,i+1),Lda,A(ip,i+1),&
     &                    Lda)
                  ENDIF
!
               ELSE
!
!                 2-by-2 pivot interchange
!
!                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)
!
                  ip = -Ipiv(i)
                  IF ( i<N ) THEN
                     IF ( ip/=(i-1) )                                   &
     &                    CALL CSWAP(N-i,A(i-1,i+1),Lda,A(ip,i+1),Lda)
                  ENDIF
!
!                 Convert IPIV
!                 There is no interchnge of rows i and and IPIV(i),
!                 so this should be reflected in IPIV format for
!                 *SYTRF_RK ( or *SYTRF_BK)
!
                  Ipiv(i) = i
!
                  i = i - 1
!
               ENDIF
               i = i - 1
            ENDDO
!
         ELSE
!
!           Revert A (A is upper)
!
!
!           Revert PERMUTATIONS and IPIV
!
!           Apply permutations to submatrices of upper part of A
!           in reverse factorization order where i increases from 1 to N
!
            i = 1
            DO WHILE ( i<=N )
               IF ( Ipiv(i)>0 ) THEN
!
!                 1-by-1 pivot interchange
!
!                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
!
                  ip = Ipiv(i)
                  IF ( i<N ) THEN
                     IF ( ip/=i ) CALL CSWAP(N-i,A(ip,i+1),Lda,A(i,i+1),&
     &                    Lda)
                  ENDIF
!
               ELSE
!
!                 2-by-2 pivot interchange
!
!                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N)
!
                  i = i + 1
                  ip = -Ipiv(i)
                  IF ( i<N ) THEN
                     IF ( ip/=(i-1) )                                   &
     &                    CALL CSWAP(N-i,A(ip,i+1),Lda,A(i-1,i+1),Lda)
                  ENDIF
!
!                 Convert IPIV
!                 There is one interchange of rows i-1 and IPIV(i-1),
!                 so this should be recorded in two consecutive entries
!                 in IPIV format for *SYTRF
!
                  Ipiv(i) = Ipiv(i-1)
!
               ENDIF
               i = i + 1
            ENDDO
!
!           Revert VALUE
!           Assign superdiagonal entries of D from array E to
!           superdiagonal entries of A.
!
            i = N
            DO WHILE ( i>1 )
               IF ( Ipiv(i)<0 ) THEN
                  A(i-1,i) = E(i)
                  i = i - 1
               ENDIF
               i = i - 1
            ENDDO
!
!        End A is UPPER
!
         ENDIF
!
!
!        Begin A is LOWER
!
      ELSEIF ( convert ) THEN
!
!           Convert A (A is lower)
!
!
!           Convert VALUE
!           Assign subdiagonal entries of D to array E and zero out
!           corresponding entries in input storage A
!
         i = 1
         E(N) = ZERO
         DO WHILE ( i<=N )
            IF ( i<N .AND. Ipiv(i)<0 ) THEN
               E(i) = A(i+1,i)
               E(i+1) = ZERO
               A(i+1,i) = ZERO
               i = i + 1
            ELSE
               E(i) = ZERO
            ENDIF
            i = i + 1
         ENDDO
!
!           Convert PERMUTATIONS and IPIV
!
!           Apply permutations to submatrices of lower part of A
!           in factorization order where k increases from 1 to N
!
         i = 1
         DO WHILE ( i<=N )
            IF ( Ipiv(i)>0 ) THEN
!
!                 1-by-1 pivot interchange
!
!                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
!
               ip = Ipiv(i)
               IF ( i>1 ) THEN
                  IF ( ip/=i ) CALL CSWAP(i-1,A(i,1),Lda,A(ip,1),Lda)
               ENDIF
!
            ELSE
!
!                 2-by-2 pivot interchange
!
!                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)
!
               ip = -Ipiv(i)
               IF ( i>1 ) THEN
                  IF ( ip/=(i+1) ) CALL CSWAP(i-1,A(i+1,1),Lda,A(ip,1), &
     &                 Lda)
               ENDIF
!
!                 Convert IPIV
!                 There is no interchnge of rows i and and IPIV(i),
!                 so this should be reflected in IPIV format for
!                 *SYTRF_RK ( or *SYTRF_BK)
!
               Ipiv(i) = i
!
               i = i + 1
!
            ENDIF
            i = i + 1
         ENDDO
!
      ELSE
!
!           Revert A (A is lower)
!
!
!           Revert PERMUTATIONS and IPIV
!
!           Apply permutations to submatrices of lower part of A
!           in reverse factorization order where i decreases from N to 1
!
         i = N
         DO WHILE ( i>=1 )
            IF ( Ipiv(i)>0 ) THEN
!
!                 1-by-1 pivot interchange
!
!                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
!
               ip = Ipiv(i)
               IF ( i>1 ) THEN
                  IF ( ip/=i ) CALL CSWAP(i-1,A(ip,1),Lda,A(i,1),Lda)
               ENDIF
!
            ELSE
!
!                 2-by-2 pivot interchange
!
!                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1)
!
               i = i - 1
               ip = -Ipiv(i)
               IF ( i>1 ) THEN
                  IF ( ip/=(i+1) ) CALL CSWAP(i-1,A(ip,1),Lda,A(i+1,1), &
     &                 Lda)
               ENDIF
!
!                 Convert IPIV
!                 There is one interchange of rows i+1 and IPIV(i+1),
!                 so this should be recorded in consecutive entries
!                 in IPIV format for *SYTRF
!
               Ipiv(i) = Ipiv(i+1)
!
            ENDIF
            i = i - 1
         ENDDO
!
!           Revert VALUE
!           Assign subdiagonal entries of D from array E to
!           subgiagonal entries of A.
!
         i = 1
         DO WHILE ( i<=N-1 )
            IF ( Ipiv(i)<0 ) THEN
               A(i+1,i) = E(i)
               i = i + 1
            ENDIF
            i = i + 1
         ENDDO
!
!
!        End A is LOWER
!
      ENDIF
 
!
!     End of CSYCONVF
!
      END SUBROUTINE CSYCONVF
