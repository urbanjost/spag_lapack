!*==cgeqrf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGEQRF VARIANT: left-looking Level 3 BLAS version of the algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEQRF ( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!  Purpose
!  =======
!
!>\details \b Purpose:
!>\verbatim
!>
!> CGEQRF computes a QR factorization of a real M-by-N matrix A:
!> A = Q * R.
!>
!> This is the left-looking Level 3 BLAS version of the algorithm.
!>
!>\endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the elements on and above the diagonal of the array
!>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
!>          upper triangular if m >= n); the elements below the diagonal,
!>          with the array TAU, represent the orthogonal matrix Q as a
!>          product of min(m,n) elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!> \endverbatim
!> \verbatim
!>          The dimension of the array WORK. The dimension can be divided into three parts.
!> \endverbatim
!> \verbatim
!>          1) The part for the triangular factor T. If the very last T is not bigger
!>             than any of the rest, then this part is NB x ceiling(K/NB), otherwise,
!>             NB x (K-NT), where K = min(M,N) and NT is the dimension of the very last T
!> \endverbatim
!> \verbatim
!>          2) The part for the very last T when T is bigger than any of the rest T.
!>             The size of this part is NT x NT, where NT = K - ceiling ((K-NX)/NB) x NB,
!>             where K = min(M,N), NX is calculated by
!>                   NX = MAX( 0, ILAENV( 3, 'CGEQRF', ' ', M, N, -1, -1 ) )
!> \endverbatim
!> \verbatim
!>          3) The part for dlarfb is of size max((N-M)*K, (N-M)*NB, K*NB, NB*NB)
!> \endverbatim
!> \verbatim
!>          So LWORK = part1 + part2 + part3
!> \endverbatim
!> \verbatim
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
!>
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
!> \ingroup variantsGEcomputational
!
!  Further Details
!  ===============
!>\details \b Further \b Details
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v'
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!>  and tau in TAU(i).
!>
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CGEQRF(M,N,A,Lda,Tau,Work,Lwork,Info)
      USE S_CGEQR2
      USE S_CLARFB
      USE S_CLARFT
      USE S_ILAENV
      USE S_SCEIL
      USE S_XERBLA
      IMPLICIT NONE
!*--CGEQRF159
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ib , iinfo , iws , j , k , lbwork , llwork ,       &
     &           lwkopt , nb , nbmin , nt , nx
      LOGICAL :: lquery
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
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
 
      Info = 0
      nbmin = 2
      nx = 0
      iws = N
      k = MIN(M,N)
      nb = ILAENV(1,'CGEQRF',' ',M,N,-1,-1)
 
!
!        Determine when to cross over from blocked to unblocked code.
!
      IF ( nb>1 .AND. nb<k ) nx = MAX(0,ILAENV(3,'CGEQRF',' ',M,N,-1,-1)&
     &                            )
!
!     Get NT, the size of the very last T, which is the left-over from in-between K-NX and K to K, eg.:
!
!            NB=3     2NB=6       K=10
!            |        |           |
!      1--2--3--4--5--6--7--8--9--10
!                  |     \________/
!               K-NX=5      NT=4
!
!     So here 4 x 4 is the last T stored in the workspace
!
      nt = k - SCEIL(REAL(k-nx)/REAL(nb))*nb
 
!
!     optimal workspace = space for dlarfb + space for normal T's + space for the last T
!
      llwork = MAX(MAX((N-M)*k,(N-M)*nb),MAX(k*nb,nb*nb))
      llwork = SCEIL(REAL(llwork)/REAL(nb))
 
      IF ( nt>nb ) THEN
 
         lbwork = k - nt
!
!         Optimal workspace for dlarfb = MAX(1,N)*NT
!
         lwkopt = (lbwork+llwork)*nb
         Work(1) = (lwkopt+nt*nt)
 
      ELSE
 
         lbwork = SCEIL(REAL(k)/REAL(nb))*nb
         lwkopt = (lbwork+llwork-nb)*nb
         Work(1) = lwkopt
 
      ENDIF
 
!
!     Test the input arguments
!
      lquery = (Lwork==-1)
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ELSEIF ( Lwork<MAX(1,N) .AND. .NOT.lquery ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGEQRF',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( k==0 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
      IF ( nb>1 .AND. nb<k ) THEN
 
         IF ( nx<k ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            IF ( nt<=nb ) THEN
               iws = (lbwork+llwork-nb)*nb
            ELSE
               iws = (lbwork+llwork)*nb + nt*nt
            ENDIF
 
            IF ( Lwork<iws ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               IF ( nt<=nb ) THEN
                  nb = Lwork/(llwork+(lbwork-nb))
               ELSE
                  nb = (Lwork-nt*nt)/(lbwork+llwork)
               ENDIF
 
               nbmin = MAX(2,ILAENV(2,'CGEQRF',' ',M,N,-1,-1))
            ENDIF
         ENDIF
      ENDIF
!
      IF ( nb>=nbmin .AND. nb<k .AND. nx<k ) THEN
!
!        Use blocked code initially
!
         DO i = 1 , k - nx , nb
            ib = MIN(k-i+1,nb)
!
!           Update the current column using old T's
!
            DO j = 1 , i - nb , nb
!
!              Apply H' to A(J:M,I:I+IB-1) from the left
!
               CALL CLARFB('Left','Transpose','Forward','Columnwise',   &
     &                     M-j+1,ib,nb,A(j,j),Lda,Work(j),lbwork,A(j,i),&
     &                     Lda,Work(lbwork*nb+nt*nt+1),ib)
 
            ENDDO
!
!           Compute the QR factorization of the current block
!           A(I:M,I:I+IB-1)
!
            CALL CGEQR2(M-i+1,ib,A(i,i),Lda,Tau(i),                     &
     &                  Work(lbwork*nb+nt*nt+1),iinfo)
 
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
!
            IF ( i+ib<=N ) CALL CLARFT('Forward','Columnwise',M-i+1,ib, &
     &                                 A(i,i),Lda,Tau(i),Work(i),lbwork)
         ENDDO
      ELSE
         i = 1
      ENDIF
!
!     Use unblocked code to factor the last or only block.
!
      IF ( i<=k ) THEN
 
         IF ( i/=1 ) THEN
 
            DO j = 1 , i - nb , nb
!
!                Apply H' to A(J:M,I:K) from the left
!
               CALL CLARFB('Left','Transpose','Forward','Columnwise',   &
     &                     M-j+1,k-i+1,nb,A(j,j),Lda,Work(j),lbwork,    &
     &                     A(j,i),Lda,Work(lbwork*nb+nt*nt+1),k-i+1)
            ENDDO
 
            CALL CGEQR2(M-i+1,k-i+1,A(i,i),Lda,Tau(i),                  &
     &                  Work(lbwork*nb+nt*nt+1),iinfo)
 
         ELSE
!
!        Use unblocked code to factor the last or only block.
!
            CALL CGEQR2(M-i+1,N-i+1,A(i,i),Lda,Tau(i),Work,iinfo)
 
         ENDIF
      ENDIF
 
 
!
!     Apply update to the column M+1:N when N > M
!
      IF ( M<N .AND. i/=1 ) THEN
!
!         Form the last triangular factor of the block reflector
!         H = H(i) H(i+1) . . . H(i+ib-1)
!
         IF ( nt<=nb ) THEN
            CALL CLARFT('Forward','Columnwise',M-i+1,k-i+1,A(i,i),Lda,  &
     &                  Tau(i),Work(i),lbwork)
         ELSE
            CALL CLARFT('Forward','Columnwise',M-i+1,k-i+1,A(i,i),Lda,  &
     &                  Tau(i),Work(lbwork*nb+1),nt)
         ENDIF
 
!
!         Apply H' to A(1:M,M+1:N) from the left
!
         DO j = 1 , k - nx , nb
 
            ib = MIN(k-j+1,nb)
 
            CALL CLARFB('Left','Transpose','Forward','Columnwise',M-j+1,&
     &                  N-M,ib,A(j,j),Lda,Work(j),lbwork,A(j,M+1),Lda,  &
     &                  Work(lbwork*nb+nt*nt+1),N-M)
 
         ENDDO
 
         IF ( nt<=nb ) THEN
            CALL CLARFB('Left','Transpose','Forward','Columnwise',M-j+1,&
     &                  N-M,k-j+1,A(j,j),Lda,Work(j),lbwork,A(j,M+1),   &
     &                  Lda,Work(lbwork*nb+nt*nt+1),N-M)
         ELSE
            CALL CLARFB('Left','Transpose','Forward','Columnwise',M-j+1,&
     &                  N-M,k-j+1,A(j,j),Lda,Work(lbwork*nb+1),nt,      &
     &                  A(j,M+1),Lda,Work(lbwork*nb+nt*nt+1),N-M)
         ENDIF
 
      ENDIF
 
      Work(1) = iws
!
!     End of CGEQRF
!
      END SUBROUTINE CGEQRF
