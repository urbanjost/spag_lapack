!*==csyconv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CSYCONV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CSYCONV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyconv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyconv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyconv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
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
!>
!> CSYCONV convert A given by TRF into L and D and vice-versa.
!> Get Non-diag elements of D (returned in workspace) and
!> apply or reverse permutation done in TRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are stored
!>          as an upper or lower triangular matrix.
!>          = 'U':  Upper triangular, form is A = U*D*U**T;
!>          = 'L':  Lower triangular, form is A = L*D*L**T.
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
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by CSYTRF.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by CSYTRF.
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is COMPLEX array, dimension (N)
!>          E stores the supdiagonal/subdiagonal of the symmetric 1-by-1
!>          or 2-by-2 block diagonal matrix D in LDLT.
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
!> \ingroup complexSYcomputational
!
!  =====================================================================
      SUBROUTINE CSYCONV(Uplo,Way,N,A,Lda,Ipiv,E,Info)
      IMPLICIT NONE
!*--CSYCONV118
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo , Way
      INTEGER Info , Lda , N
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX A(Lda,*) , E(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0E+0,0.0E+0))
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     .. Local Scalars ..
      LOGICAL upper , convert
      INTEGER i , ip , j
      COMPLEX temp
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
         CALL XERBLA('CSYCONV',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( upper ) THEN
!
!      A is UPPER
!
!      Convert A (A is upper)
!
!        Convert VALUE
!
         IF ( convert ) THEN
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
!        Convert PERMUTATIONS
!
            i = N
            DO WHILE ( i>=1 )
               IF ( Ipiv(i)>0 ) THEN
                  ip = Ipiv(i)
                  IF ( i<N ) THEN
                     DO j = i + 1 , N
                        temp = A(ip,j)
                        A(ip,j) = A(i,j)
                        A(i,j) = temp
                     ENDDO
                  ENDIF
               ELSE
                  ip = -Ipiv(i)
                  IF ( i<N ) THEN
                     DO j = i + 1 , N
                        temp = A(ip,j)
                        A(ip,j) = A(i-1,j)
                        A(i-1,j) = temp
                     ENDDO
                  ENDIF
                  i = i - 1
               ENDIF
               i = i - 1
            ENDDO
 
         ELSE
!
!      Revert A (A is upper)
!
!
!        Revert PERMUTATIONS
!
            i = 1
            DO WHILE ( i<=N )
               IF ( Ipiv(i)>0 ) THEN
                  ip = Ipiv(i)
                  IF ( i<N ) THEN
                     DO j = i + 1 , N
                        temp = A(ip,j)
                        A(ip,j) = A(i,j)
                        A(i,j) = temp
                     ENDDO
                  ENDIF
               ELSE
                  ip = -Ipiv(i)
                  i = i + 1
                  IF ( i<N ) THEN
                     DO j = i + 1 , N
                        temp = A(ip,j)
                        A(ip,j) = A(i-1,j)
                        A(i-1,j) = temp
                     ENDDO
                  ENDIF
               ENDIF
               i = i + 1
            ENDDO
!
!        Revert VALUE
!
            i = N
            DO WHILE ( i>1 )
               IF ( Ipiv(i)<0 ) THEN
                  A(i-1,i) = E(i)
                  i = i - 1
               ENDIF
               i = i - 1
            ENDDO
         ENDIF
!
!      A is LOWER
!
      ELSEIF ( convert ) THEN
!
!      Convert A (A is lower)
!
!
!        Convert VALUE
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
!        Convert PERMUTATIONS
!
         i = 1
         DO WHILE ( i<=N )
            IF ( Ipiv(i)>0 ) THEN
               ip = Ipiv(i)
               IF ( i>1 ) THEN
                  DO j = 1 , i - 1
                     temp = A(ip,j)
                     A(ip,j) = A(i,j)
                     A(i,j) = temp
                  ENDDO
               ENDIF
            ELSE
               ip = -Ipiv(i)
               IF ( i>1 ) THEN
                  DO j = 1 , i - 1
                     temp = A(ip,j)
                     A(ip,j) = A(i+1,j)
                     A(i+1,j) = temp
                  ENDDO
               ENDIF
               i = i + 1
            ENDIF
            i = i + 1
         ENDDO
      ELSE
!
!      Revert A (A is lower)
!
!
!        Revert PERMUTATIONS
!
         i = N
         DO WHILE ( i>=1 )
            IF ( Ipiv(i)>0 ) THEN
               ip = Ipiv(i)
               IF ( i>1 ) THEN
                  DO j = 1 , i - 1
                     temp = A(i,j)
                     A(i,j) = A(ip,j)
                     A(ip,j) = temp
                  ENDDO
               ENDIF
            ELSE
               ip = -Ipiv(i)
               i = i - 1
               IF ( i>1 ) THEN
                  DO j = 1 , i - 1
                     temp = A(i+1,j)
                     A(i+1,j) = A(ip,j)
                     A(ip,j) = temp
                  ENDDO
               ENDIF
            ENDIF
            i = i - 1
         ENDDO
!
!        Revert VALUE
!
         i = 1
         DO WHILE ( i<=N-1 )
            IF ( Ipiv(i)<0 ) THEN
               A(i+1,i) = E(i)
               i = i + 1
            ENDIF
            i = i + 1
         ENDDO
      ENDIF
 
!
!     End of CSYCONV
!
      END SUBROUTINE CSYCONV
