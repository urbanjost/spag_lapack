!*==dsytrd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DSYTRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSYTRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYTRD reduces a real symmetric matrix A to real symmetric
!> tridiagonal form T by an orthogonal similarity transformation:
!> Q**T * A * Q = T.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!>          of A are overwritten by the corresponding elements of the
!>          tridiagonal matrix T, and the elements above the first
!>          superdiagonal, with the array TAU, represent the orthogonal
!>          matrix Q as a product of elementary reflectors; if UPLO
!>          = 'L', the diagonal and first subdiagonal of A are over-
!>          written by the corresponding elements of the tridiagonal
!>          matrix T, and the elements below the first subdiagonal, with
!>          the array TAU, represent the orthogonal matrix Q as a product
!>          of elementary reflectors. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal elements of the tridiagonal matrix T:
!>          D(i) = A(i,i).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The off-diagonal elements of the tridiagonal matrix T:
!>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= 1.
!>          For optimum performance LWORK >= N*NB, where NB is the
!>          optimal blocksize.
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
!> \ingroup doubleSYcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(n-1) . . . H(2) H(1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!>  A(1:i-1,i+1), and tau in TAU(i).
!>
!>  If UPLO = 'L', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(1) H(2) . . . H(n-1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!>  and tau in TAU(i).
!>
!>  The contents of A on exit are illustrated by the following examples
!>  with n = 5:
!>
!>  if UPLO = 'U':                       if UPLO = 'L':
!>
!>    (  d   e   v2  v3  v4 )              (  d                  )
!>    (      d   e   v3  v4 )              (  e   d              )
!>    (          d   e   v4 )              (  v1  e   d          )
!>    (              d   e  )              (  v1  v2  e   d      )
!>    (                  d  )              (  v1  v2  v3  e   d  )
!>
!>  where d and e denote diagonal and off-diagonal elements of T, and vi
!>  denotes an element of the vector defining H(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DSYTRD(Uplo,N,A,Lda,D,E,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_DLATRD
      USE S_DSYR2K
      USE S_DSYTD2
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DSYTRD203
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , iinfo , iws , j , kk , ldwork , lwkopt , nb ,      &
     &           nbmin , nx
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
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
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
         Info = -9
      ENDIF
!
      IF ( Info==0 ) THEN
!
!        Determine the block size.
!
         nb = ILAENV(1,'DSYTRD',Uplo,N,-1,-1,-1)
         lwkopt = N*nb
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSYTRD',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
      nx = N
      iws = 1
      IF ( nb>1 .AND. nb<N ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code).
!
         nx = MAX(nb,ILAENV(3,'DSYTRD',Uplo,N,-1,-1,-1))
         IF ( nx<N ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            ldwork = N
            iws = ldwork*nb
            IF ( Lwork<iws ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code by setting NX = N.
!
               nb = MAX(Lwork/ldwork,1)
               nbmin = ILAENV(2,'DSYTRD',Uplo,N,-1,-1,-1)
               IF ( nb<nbmin ) nx = N
            ENDIF
         ELSE
            nx = N
         ENDIF
      ELSE
         nb = 1
      ENDIF
!
      IF ( upper ) THEN
!
!        Reduce the upper triangle of A.
!        Columns 1:kk are handled by the unblocked method.
!
         kk = N - ((N-nx+nb-1)/nb)*nb
         DO i = N - nb + 1 , kk + 1 , -nb
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
            CALL DLATRD(Uplo,i+nb-1,nb,A,Lda,E,Tau,Work,ldwork)
!
!           Update the unreduced submatrix A(1:i-1,1:i-1), using an
!           update of the form:  A := A - V*W**T - W*V**T
!
            CALL DSYR2K(Uplo,'No transpose',i-1,nb,-ONE,A(1,i),Lda,Work,&
     &                  ldwork,ONE,A,Lda)
!
!           Copy superdiagonal elements back into A, and diagonal
!           elements into D
!
            DO j = i , i + nb - 1
               A(j-1,j) = E(j-1)
               D(j) = A(j,j)
            ENDDO
         ENDDO
!
!        Use unblocked code to reduce the last or only block
!
         CALL DSYTD2(Uplo,kk,A,Lda,D,E,Tau,iinfo)
      ELSE
!
!        Reduce the lower triangle of A
!
         DO i = 1 , N - nx , nb
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
            CALL DLATRD(Uplo,N-i+1,nb,A(i,i),Lda,E(i),Tau(i),Work,      &
     &                  ldwork)
!
!           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
!           an update of the form:  A := A - V*W**T - W*V**T
!
            CALL DSYR2K(Uplo,'No transpose',N-i-nb+1,nb,-ONE,A(i+nb,i), &
     &                  Lda,Work(nb+1),ldwork,ONE,A(i+nb,i+nb),Lda)
!
!           Copy subdiagonal elements back into A, and diagonal
!           elements into D
!
            DO j = i , i + nb - 1
               A(j+1,j) = E(j)
               D(j) = A(j,j)
            ENDDO
         ENDDO
!
!        Use unblocked code to reduce the last or only block
!
         CALL DSYTD2(Uplo,N-i+1,A(i,i),Lda,D(i),E(i),Tau(i),iinfo)
      ENDIF
!
      Work(1) = lwkopt
!
!     End of DSYTRD
!
      END SUBROUTINE DSYTRD
