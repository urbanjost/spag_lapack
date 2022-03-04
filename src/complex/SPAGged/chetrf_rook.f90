!*==chetrf_rook.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CHETRF_ROOK computes the factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method (blocked algorithm, calling Level 3 BLAS).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHETRF_ROOK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrf_rook.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrf_rook.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrf_rook.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHETRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHETRF_ROOK computes the factorization of a complex Hermitian matrix A
!> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
!> The form of the factorization is
!>
!>    A = U*D*U**T  or  A = L*D*L**T
!>
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, and D is Hermitian and block diagonal with
!> 1-by-1 and 2-by-2 diagonal blocks.
!>
!> This is the blocked version of the algorithm, calling Level 3 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
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
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, the block diagonal matrix D and the multipliers used
!>          to obtain the factor U or L (see below for further details).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D.
!>
!>          If UPLO = 'U':
!>             Only the last KB elements of IPIV are set.
!>
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>             interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
!>             columns k and -IPIV(k) were interchanged and rows and
!>             columns k-1 and -IPIV(k-1) were inerchaged,
!>             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>             Only the first KB elements of IPIV are set.
!>
!>             If IPIV(k) > 0, then rows and columns k and IPIV(k)
!>             were interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
!>             columns k and -IPIV(k) were interchanged and rows and
!>             columns k+1 and -IPIV(k+1) were inerchaged,
!>             D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK)).
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of WORK.  LWORK >=1.  For best performance
!>          LWORK >= N*NB, where NB is the block size returned by ILAENV.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
!>                has been completed, but the block diagonal matrix D is
!>                exactly singular, and division by zero will occur if it
!>                is used to solve a system of equations.
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
!> \date June 2016
!
!> \ingroup complexHEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', then A = U*D*U**T, where
!>     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!>
!>             (   I    v    0   )   k-s
!>     U(k) =  (   0    I    0   )   s
!>             (   0    0    I   )   n-k
!>                k-s   s   n-k
!>
!>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!>  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!>
!>  If UPLO = 'L', then A = L*D*L**T, where
!>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!>
!>             (   I    0     0   )  k-1
!>     L(k) =  (   0    I     0   )  s
!>             (   0    v     I   )  n-k-s+1
!>                k-1   s  n-k-s+1
!>
!>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>  June 2016,  Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
!>                  School of Mathematics,
!>                  University of Manchester
!>
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE CHETRF_ROOK(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
      USE S_CHETF2_ROOK
      USE S_CLAHEF_ROOK
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CHETRF_ROOK221
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: iinfo , iws , j , k , kb , ldwork , lwkopt , nb , nbmin
      LOGICAL :: lquery , upper
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
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
      lquery = (Lwork==-1)
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Lwork<1 .AND. .NOT.lquery ) THEN
         Info = -7
      ENDIF
!
      IF ( Info==0 ) THEN
!
!        Determine the block size
!
         nb = ILAENV(1,'CHETRF_ROOK',Uplo,N,-1,-1,-1)
         lwkopt = MAX(1,N*nb)
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CHETRF_ROOK',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
      nbmin = 2
      ldwork = N
      IF ( nb>1 .AND. nb<N ) THEN
         iws = ldwork*nb
         IF ( Lwork<iws ) THEN
            nb = MAX(Lwork/ldwork,1)
            nbmin = MAX(2,ILAENV(2,'CHETRF_ROOK',Uplo,N,-1,-1,-1))
         ENDIF
      ELSE
         iws = 1
      ENDIF
      IF ( nb<nbmin ) nb = N
!
      IF ( upper ) THEN
!
!        Factorize A as U*D*U**T using the upper triangle of A
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        KB, where KB is the number of columns factorized by CLAHEF_ROOK;
!        KB is either NB or NB-1, or K for the last block
!
         k = N
!
!        If K < 1, exit from loop
!
         DO WHILE ( k>=1 )
!
            IF ( k>nb ) THEN
!
!           Factorize columns k-kb+1:k of A and use blocked code to
!           update columns 1:k-kb
!
               CALL CLAHEF_ROOK(Uplo,k,nb,kb,A,Lda,Ipiv,Work,ldwork,    &
     &                          iinfo)
            ELSE
!
!           Use unblocked code to factorize columns 1:k of A
!
               CALL CHETF2_ROOK(Uplo,k,A,Lda,Ipiv,iinfo)
               kb = k
            ENDIF
!
!        Set INFO on the first occurrence of a zero pivot
!
            IF ( Info==0 .AND. iinfo>0 ) Info = iinfo
!
!        No need to adjust IPIV
!
!        Decrease K and return to the start of the main loop
!
            k = k - kb
         ENDDO
!
         Work(1) = lwkopt
!
      ELSE
!
!        Factorize A as L*D*L**T using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        KB, where KB is the number of columns factorized by CLAHEF_ROOK;
!        KB is either NB or NB-1, or N-K+1 for the last block
!
         k = 1
!
!        If K > N, exit from loop
!
         DO WHILE ( k<=N )
!
            IF ( k<=N-nb ) THEN
!
!           Factorize columns k:k+kb-1 of A and use blocked code to
!           update columns k+kb:n
!
               CALL CLAHEF_ROOK(Uplo,N-k+1,nb,kb,A(k,k),Lda,Ipiv(k),    &
     &                          Work,ldwork,iinfo)
            ELSE
!
!           Use unblocked code to factorize columns k:n of A
!
               CALL CHETF2_ROOK(Uplo,N-k+1,A(k,k),Lda,Ipiv(k),iinfo)
               kb = N - k + 1
            ENDIF
!
!        Set INFO on the first occurrence of a zero pivot
!
            IF ( Info==0 .AND. iinfo>0 ) Info = iinfo + k - 1
!
!        Adjust IPIV
!
            DO j = k , k + kb - 1
               IF ( Ipiv(j)>0 ) THEN
                  Ipiv(j) = Ipiv(j) + k - 1
               ELSE
                  Ipiv(j) = Ipiv(j) - k + 1
               ENDIF
            ENDDO
!
!        Increase K and return to the start of the main loop
!
            k = k + kb
         ENDDO
         Work(1) = lwkopt
!
      ENDIF
!
!     End of CHETRF_ROOK
!
      END SUBROUTINE CHETRF_ROOK
