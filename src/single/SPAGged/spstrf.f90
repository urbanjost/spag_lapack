!*==spstrf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SPSTRF computes the Cholesky factorization with complete pivoting of a real symmetric positive semidefinite matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SPSTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spstrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spstrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spstrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
!
!       .. Scalar Arguments ..
!       REAL               TOL
!       INTEGER            INFO, LDA, N, RANK
!       CHARACTER          UPLO
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), WORK( 2*N )
!       INTEGER            PIV( N )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SPSTRF computes the Cholesky factorization with complete
!> pivoting of a real symmetric positive semidefinite matrix A.
!>
!> The factorization has the form
!>    P**T * A * P = U**T * U ,  if UPLO = 'U',
!>    P**T * A * P = L  * L**T,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular, and
!> P is stored as vector PIV.
!>
!> This algorithm does not attempt to check that A is positive
!> semidefinite. This version of the algorithm calls level 3 BLAS.
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
!>          A is REAL array, dimension (LDA,N)
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
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
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
!>          User defined tolerance. If TOL < 0, then N*U*MAX( A(K,K) )
!>          will be used. The algorithm terminates at the (K-1)st step
!>          if the pivot <= TOL.
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
!> \ingroup realOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SPSTRF(Uplo,N,A,Lda,Piv,Rank,Tol,Work,Info)
      USE S_ILAENV
      USE S_LSAME
      USE S_SGEMV
      USE S_SISNAN
      USE S_SLAMCH
      USE S_SPSTF2
      USE S_SSCAL
      USE S_SSWAP
      USE S_SSYRK
      USE S_XERBLA
      IMPLICIT NONE
!*--SPSTRF155
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(N) :: Piv
      INTEGER :: Rank
      REAL :: Tol
      REAL , INTENT(INOUT) , DIMENSION(2*N) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ajj , sstop , stemp
      INTEGER :: i , itemp , j , jb , k , nb , pvt
      LOGICAL :: upper
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
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SPSTRF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Get block size
!
      nb = ILAENV(1,'SPOTRF',Uplo,N,-1,-1,-1)
      IF ( nb<=1 .OR. nb>=N ) THEN
!
!        Use unblocked code
!
         CALL SPSTF2(Uplo,N,A(1,1),Lda,Piv,Rank,Tol,Work,Info)
         GOTO 99999
!
      ELSE
!
!     Initialize PIV
!
         DO i = 1 , N
            Piv(i) = i
         ENDDO
!
!     Compute stopping value
!
         pvt = 1
         ajj = A(pvt,pvt)
         DO i = 2 , N
            IF ( A(i,i)>ajj ) THEN
               pvt = i
               ajj = A(pvt,pvt)
            ENDIF
         ENDDO
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
!
         IF ( upper ) THEN
!
!           Compute the Cholesky factorization P**T * A * P = U**T * U
!
            DO k = 1 , N , nb
!
!              Account for last block not being NB wide
!
               jb = MIN(nb,N-k+1)
!
!              Set relevant part of first half of WORK to zero,
!              holds dot products
!
               DO i = k , N
                  Work(i) = 0
               ENDDO
!
               DO j = k , k + jb - 1
!
!              Find pivot, test for exit, else swap rows and columns
!              Update dot products, compute possible pivots which are
!              stored in the second half of WORK
!
                  DO i = j , N
!
                     IF ( j>k ) Work(i) = Work(i) + A(j-1,i)**2
                     Work(N+i) = A(i,i) - Work(i)
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
!                    Pivot OK, so can now swap pivot rows and columns
!
                     A(pvt,pvt) = A(j,j)
                     CALL SSWAP(j-1,A(1,j),1,A(1,pvt),1)
                     IF ( pvt<N )                                       &
     &                    CALL SSWAP(N-pvt,A(j,pvt+1),Lda,A(pvt,pvt+1), &
     &                    Lda)
                     CALL SSWAP(pvt-j-1,A(j,j+1),Lda,A(j+1,pvt),1)
!
!                    Swap dot products and PIV
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
!                 Compute elements J+1:N of row J.
!
                  IF ( j<N ) THEN
                     CALL SGEMV('Trans',j-k,N-j,-ONE,A(k,j+1),Lda,A(k,j)&
     &                          ,1,ONE,A(j,j+1),Lda)
                     CALL SSCAL(N-j,ONE/ajj,A(j,j+1),Lda)
                  ENDIF
!
               ENDDO
!
!              Update trailing matrix, J already incremented
!
               IF ( k+jb<=N ) CALL SSYRK('Upper','Trans',N-j+1,jb,-ONE, &
     &              A(k,j),Lda,ONE,A(j,j),Lda)
!
            ENDDO
!
         ELSE
!
!        Compute the Cholesky factorization P**T * A * P = L * L**T
!
            DO k = 1 , N , nb
!
!              Account for last block not being NB wide
!
               jb = MIN(nb,N-k+1)
!
!              Set relevant part of first half of WORK to zero,
!              holds dot products
!
               DO i = k , N
                  Work(i) = 0
               ENDDO
!
               DO j = k , k + jb - 1
!
!              Find pivot, test for exit, else swap rows and columns
!              Update dot products, compute possible pivots which are
!              stored in the second half of WORK
!
                  DO i = j , N
!
                     IF ( j>k ) Work(i) = Work(i) + A(i,j-1)**2
                     Work(N+i) = A(i,i) - Work(i)
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
!                    Pivot OK, so can now swap pivot rows and columns
!
                     A(pvt,pvt) = A(j,j)
                     CALL SSWAP(j-1,A(j,1),Lda,A(pvt,1),Lda)
                     IF ( pvt<N )                                       &
     &                    CALL SSWAP(N-pvt,A(pvt+1,j),1,A(pvt+1,pvt),1)
                     CALL SSWAP(pvt-j-1,A(j+1,j),1,A(pvt,j+1),Lda)
!
!                    Swap dot products and PIV
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
!                 Compute elements J+1:N of column J.
!
                  IF ( j<N ) THEN
                     CALL SGEMV('No Trans',N-j,j-k,-ONE,A(j+1,k),Lda,   &
     &                          A(j,k),Lda,ONE,A(j+1,j),1)
                     CALL SSCAL(N-j,ONE/ajj,A(j+1,j),1)
                  ENDIF
!
               ENDDO
!
!              Update trailing matrix, J already incremented
!
               IF ( k+jb<=N ) CALL SSYRK('Lower','No Trans',N-j+1,jb,   &
     &              -ONE,A(j,k),Lda,ONE,A(j,j),Lda)
!
            ENDDO
!
         ENDIF
      ENDIF
!
!     Ran to completion, A has full rank
!
      Rank = N
!
      GOTO 99999
!
!     Rank is the number of steps completed.  Set INFO = 1 to signal
!     that the factorization cannot be used to solve a system.
!
 100  Rank = j - 1
      Info = 1
!
!
!     End of SPSTRF
!
99999 END SUBROUTINE SPSTRF
