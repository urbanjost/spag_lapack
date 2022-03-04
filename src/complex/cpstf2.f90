!*==cpstf2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CPSTF2 computes the Cholesky factorization with complete pivoting of complex Hermitian positive semidefinite matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPSTF2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpstf2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpstf2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpstf2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPSTF2( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
!
!       .. Scalar Arguments ..
!       REAL               TOL
!       INTEGER            INFO, LDA, N, RANK
!       CHARACTER          UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * )
!       REAL               WORK( 2*N )
!       INTEGER            PIV( N )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPSTF2 computes the Cholesky factorization with complete
!> pivoting of a complex Hermitian positive semidefinite matrix A.
!>
!> The factorization has the form
!>    P**T * A * P = U**H * U ,  if UPLO = 'U',
!>    P**T * A * P = L  * L**H,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular, and
!> P is stored as vector PIV.
!>
!> This algorithm does not attempt to check that A is positive
!> semidefinite. This version of the algorithm calls level 2 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
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
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          n by n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n by n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the factor U or L from the Cholesky
!>          factorization as above.
!> \endverbatim
!>
!> \param[out] PIV
!> \verbatim
!>          PIV is INTEGER array, dimension (N)
!>          PIV is such that the nonzero entries are P( PIV(K), K ) = 1.
!> \endverbatim
!>
!> \param[out] RANK
!> \verbatim
!>          RANK is INTEGER
!>          The rank of A given by the number of steps the algorithm
!>          completed.
!> \endverbatim
!>
!> \param[in] TOL
!> \verbatim
!>          TOL is REAL
!>          User defined tolerance. If TOL < 0, then N*U*MAX( A( K,K ) )
!>          will be used. The algorithm terminates at the (K-1)st step
!>          if the pivot <= TOL.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (2*N)
!>          Work space.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          < 0: If INFO = -K, the K-th argument had an illegal value,
!>          = 0: algorithm completed successfully, and
!>          > 0: the matrix A is either rank deficient with computed rank
!>               as returned in RANK, or is not positive semidefinite. See
!>               Section 7 of LAPACK Working Note #161 for further
!>               information.
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
!  =====================================================================
      SUBROUTINE CPSTF2(Uplo,N,A,Lda,Piv,Rank,Tol,Work,Info)
      IMPLICIT NONE
!*--CPSTF2146
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      REAL Tol
      INTEGER Info , Lda , N , Rank
      CHARACTER Uplo
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*)
      REAL Work(2*N)
      INTEGER Piv(N)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      COMPLEX CONE
      PARAMETER (CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      COMPLEX ctemp
      REAL ajj , sstop , stemp
      INTEGER i , itemp , j , pvt
      LOGICAL upper
!     ..
!     .. External Functions ..
      REAL SLAMCH
      LOGICAL LSAME , SISNAN
      EXTERNAL SLAMCH , LSAME , SISNAN
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMV , CLACGV , CSSCAL , CSWAP , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CONJG , MAX , REAL , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CPSTF2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Initialize PIV
!
      DO i = 1 , N
         Piv(i) = i
      ENDDO
!
!     Compute stopping value
!
      DO i = 1 , N
         Work(i) = REAL(A(i,i))
      ENDDO
      pvt = MAXLOC(Work(1:N),1)
      ajj = REAL(A(pvt,pvt))
      IF ( ajj<=ZERO .OR. SISNAN(ajj) ) THEN
         Rank = 0
         Info = 1
         GOTO 99999
      ENDIF
!
!     Compute stopping value if not supplied
!
      IF ( Tol<ZERO ) THEN
         sstop = N*SLAMCH('Epsilon')*ajj
      ELSE
         sstop = Tol
      ENDIF
!
!     Set first half of WORK to zero, holds dot products
!
      DO i = 1 , N
         Work(i) = 0
      ENDDO
!
      IF ( upper ) THEN
!
!        Compute the Cholesky factorization P**T * A * P = U**H * U
!
         DO j = 1 , N
!
!        Find pivot, test for exit, else swap rows and columns
!        Update dot products, compute possible pivots which are
!        stored in the second half of WORK
!
            DO i = j , N
!
               IF ( j>1 ) Work(i) = Work(i)                             &
     &                              + REAL(CONJG(A(j-1,i))*A(j-1,i))
               Work(N+i) = REAL(A(i,i)) - Work(i)
!
            ENDDO
!
            IF ( j>1 ) THEN
               itemp = MAXLOC(Work((N+j):(2*N)),1)
               pvt = itemp + j - 1
               ajj = Work(N+pvt)
               IF ( ajj<=sstop .OR. SISNAN(ajj) ) THEN
                  A(j,j) = ajj
                  GOTO 100
               ENDIF
            ENDIF
!
            IF ( j/=pvt ) THEN
!
!              Pivot OK, so can now swap pivot rows and columns
!
               A(pvt,pvt) = A(j,j)
               CALL CSWAP(j-1,A(1,j),1,A(1,pvt),1)
               IF ( pvt<N ) CALL CSWAP(N-pvt,A(j,pvt+1),Lda,A(pvt,pvt+1)&
     &                                 ,Lda)
               DO i = j + 1 , pvt - 1
                  ctemp = CONJG(A(j,i))
                  A(j,i) = CONJG(A(i,pvt))
                  A(i,pvt) = ctemp
               ENDDO
               A(j,pvt) = CONJG(A(j,pvt))
!
!              Swap dot products and PIV
!
               stemp = Work(j)
               Work(j) = Work(pvt)
               Work(pvt) = stemp
               itemp = Piv(pvt)
               Piv(pvt) = Piv(j)
               Piv(j) = itemp
            ENDIF
!
            ajj = SQRT(ajj)
            A(j,j) = ajj
!
!           Compute elements J+1:N of row J
!
            IF ( j<N ) THEN
               CALL CLACGV(j-1,A(1,j),1)
               CALL CGEMV('Trans',j-1,N-j,-CONE,A(1,j+1),Lda,A(1,j),1,  &
     &                    CONE,A(j,j+1),Lda)
               CALL CLACGV(j-1,A(1,j),1)
               CALL CSSCAL(N-j,ONE/ajj,A(j,j+1),Lda)
            ENDIF
!
         ENDDO
!
      ELSE
!
!        Compute the Cholesky factorization P**T * A * P = L * L**H
!
         DO j = 1 , N
!
!        Find pivot, test for exit, else swap rows and columns
!        Update dot products, compute possible pivots which are
!        stored in the second half of WORK
!
            DO i = j , N
!
               IF ( j>1 ) Work(i) = Work(i)                             &
     &                              + REAL(CONJG(A(i,j-1))*A(i,j-1))
               Work(N+i) = REAL(A(i,i)) - Work(i)
!
            ENDDO
!
            IF ( j>1 ) THEN
               itemp = MAXLOC(Work((N+j):(2*N)),1)
               pvt = itemp + j - 1
               ajj = Work(N+pvt)
               IF ( ajj<=sstop .OR. SISNAN(ajj) ) THEN
                  A(j,j) = ajj
                  GOTO 100
               ENDIF
            ENDIF
!
            IF ( j/=pvt ) THEN
!
!              Pivot OK, so can now swap pivot rows and columns
!
               A(pvt,pvt) = A(j,j)
               CALL CSWAP(j-1,A(j,1),Lda,A(pvt,1),Lda)
               IF ( pvt<N ) CALL CSWAP(N-pvt,A(pvt+1,j),1,A(pvt+1,pvt), &
     &                                 1)
               DO i = j + 1 , pvt - 1
                  ctemp = CONJG(A(i,j))
                  A(i,j) = CONJG(A(pvt,i))
                  A(pvt,i) = ctemp
               ENDDO
               A(pvt,j) = CONJG(A(pvt,j))
!
!              Swap dot products and PIV
!
               stemp = Work(j)
               Work(j) = Work(pvt)
               Work(pvt) = stemp
               itemp = Piv(pvt)
               Piv(pvt) = Piv(j)
               Piv(j) = itemp
            ENDIF
!
            ajj = SQRT(ajj)
            A(j,j) = ajj
!
!           Compute elements J+1:N of column J
!
            IF ( j<N ) THEN
               CALL CLACGV(j-1,A(j,1),Lda)
               CALL CGEMV('No Trans',N-j,j-1,-CONE,A(j+1,1),Lda,A(j,1), &
     &                    Lda,CONE,A(j+1,j),1)
               CALL CLACGV(j-1,A(j,1),Lda)
               CALL CSSCAL(N-j,ONE/ajj,A(j+1,j),1)
            ENDIF
!
         ENDDO
!
      ENDIF
!
!     Ran to completion, A has full rank
!
      Rank = N
!
      GOTO 99999
!
!     Rank is number of steps completed.  Set INFO = 1 to signal
!     that the factorization cannot be used to solve a system.
!
 100  Rank = j - 1
      Info = 1
!
!
!     End of CPSTF2
!
99999 END SUBROUTINE CPSTF2
