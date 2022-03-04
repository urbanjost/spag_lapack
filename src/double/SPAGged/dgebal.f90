!*==dgebal.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DGEBAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGEBAL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebal.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebal.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebal.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB
!       INTEGER            IHI, ILO, INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), SCALE( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEBAL balances a general real matrix A.  This involves, first,
!> permuting A by a similarity transformation to isolate eigenvalues
!> in the first 1 to ILO-1 and last IHI+1 to N elements on the
!> diagonal; and second, applying a diagonal similarity transformation
!> to rows and columns ILO to IHI to make the rows and columns as
!> close in norm as possible.  Both steps are optional.
!>
!> Balancing may reduce the 1-norm of the matrix, and improve the
!> accuracy of the computed eigenvalues and/or eigenvectors.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies the operations to be performed on A:
!>          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0
!>                  for i = 1,...,N;
!>          = 'P':  permute only;
!>          = 'S':  scale only;
!>          = 'B':  both permute and scale.
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
!>          On entry, the input matrix A.
!>          On exit,  A is overwritten by the balanced matrix.
!>          If JOB = 'N', A is not referenced.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!> \param[out] IHI
!> \verbatim
!>          IHI is INTEGER
!>          ILO and IHI are set to integers such that on exit
!>          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
!>          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION array, dimension (N)
!>          Details of the permutations and scaling factors applied to
!>          A.  If P(j) is the index of the row and column interchanged
!>          with row and column j and D(j) is the scaling factor
!>          applied to row and column j, then
!>          SCALE(j) = P(j)    for j = 1,...,ILO-1
!>                   = D(j)    for j = ILO,...,IHI
!>                   = P(j)    for j = IHI+1,...,N.
!>          The order in which the interchanges are made is N to IHI+1,
!>          then 1 to ILO-1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \date June 2017
!
!> \ingroup doubleGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The permutations consist of row and column interchanges which put
!>  the matrix in the form
!>
!>             ( T1   X   Y  )
!>     P A P = (  0   B   Z  )
!>             (  0   0   T2 )
!>
!>  where T1 and T2 are upper triangular matrices whose eigenvalues lie
!>  along the diagonal.  The column indices ILO and IHI mark the starting
!>  and ending columns of the submatrix B. Balancing consists of applying
!>  a diagonal similarity transformation inv(D) * B * D to make the
!>  1-norms of each row of B and its corresponding column nearly equal.
!>  The output matrix is
!>
!>     ( T1     X*D          Y    )
!>     (  0  inv(D)*B*D  inv(D)*Z ).
!>     (  0      0           T2   )
!>
!>  Information about the permutations P and the diagonal matrix D is
!>  returned in the vector SCALE.
!>
!>  This subroutine is based on the EISPACK routine BALANC.
!>
!>  Modified by Tzu-Yi Chen, Computer Science Division, University of
!>    California at Berkeley, USA
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DGEBAL(Job,N,A,Lda,Ilo,Ihi,Scale,Info)
      USE F77KINDS                        
      USE S_DISNAN
      USE S_DLAMCH
      USE S_DNRM2
      USE S_DSCAL
      USE S_DSWAP
      USE S_IDAMAX
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DGEBAL173
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              SCLFAC = 2.0D+0 , FACTOR = 0.95D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Job
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) :: Ilo
      INTEGER , INTENT(OUT) :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Scale
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: c , ca , f , g , r , ra , s , sfmax1 , sfmax2 ,   &
     &                sfmin1 , sfmin2
      INTEGER :: i , ica , iexc , ira , j , k , l , m
      LOGICAL :: noconv
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
!     Test the input parameters
!
      Info = 0
      IF ( .NOT.LSAME(Job,'N') .AND. .NOT.LSAME(Job,'P') .AND.          &
     &     .NOT.LSAME(Job,'S') .AND. .NOT.LSAME(Job,'B') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGEBAL',-Info)
         RETURN
      ENDIF
!
      k = 1
      l = N
!
      IF ( N==0 ) GOTO 700
!
      IF ( LSAME(Job,'N') ) THEN
         DO i = 1 , N
            Scale(i) = ONE
         ENDDO
         GOTO 700
      ENDIF
!
!
!     Permutation to isolate eigenvalues if possible
!
      IF ( LSAME(Job,'S') ) GOTO 600
      GOTO 200
!
!     Row and column exchange.
!
 100  Scale(m) = j
      IF ( j/=m ) THEN
!
         CALL DSWAP(l,A(1,j),1,A(1,m),1)
         CALL DSWAP(N-k+1,A(j,k),Lda,A(m,k),Lda)
      ENDIF
!
      IF ( iexc==2 ) THEN
!
!     Search for columns isolating an eigenvalue and push them left.
!
         k = k + 1
         GOTO 400
      ELSE
!
!     Search for rows isolating an eigenvalue and push them down.
!
         IF ( l==1 ) GOTO 700
         l = l - 1
      ENDIF
!
 200  DO j = l , 1 , -1
!
         DO i = 1 , l
            IF ( i/=j ) THEN
               IF ( A(j,i)/=ZERO ) GOTO 300
            ENDIF
         ENDDO
!
         m = l
         iexc = 1
         GOTO 100
!
 300  ENDDO
!
 400  DO j = k , l
!
         DO i = k , l
            IF ( i/=j ) THEN
               IF ( A(i,j)/=ZERO ) GOTO 500
            ENDIF
         ENDDO
!
         m = k
         iexc = 2
         GOTO 100
 500  ENDDO
!
 600  DO i = k , l
         Scale(i) = ONE
      ENDDO
!
      IF ( .NOT.(LSAME(Job,'P')) ) THEN
!
!     Balance the submatrix in rows K to L.
!
!     Iterative loop for norm reduction
!
         sfmin1 = DLAMCH('S')/DLAMCH('P')
         sfmax1 = ONE/sfmin1
         sfmin2 = sfmin1*SCLFAC
         sfmax2 = ONE/sfmin2
         DO
!
            noconv = .FALSE.
!
            DO i = k , l
!
               c = DNRM2(l-k+1,A(k,i),1)
               r = DNRM2(l-k+1,A(i,k),Lda)
               ica = IDAMAX(l,A(1,i),1)
               ca = ABS(A(ica,i))
               ira = IDAMAX(N-k+1,A(i,k),Lda)
               ra = ABS(A(i,ira+k-1))
!
!        Guard against zero C or R due to underflow.
!
               IF ( c/=ZERO .AND. r/=ZERO ) THEN
                  g = r/SCLFAC
                  f = ONE
                  s = c + r
                  DO WHILE ( c<g .AND. MAX(f,c,ca)<sfmax2 .AND.         &
     &                       MIN(r,g,ra)>sfmin2 )
                     IF ( DISNAN(c+f+ca+r+g+ra) ) THEN
!
!           Exit if NaN to avoid infinite loop
!
                        Info = -3
                        CALL XERBLA('DGEBAL',-Info)
                        RETURN
                     ENDIF
                     f = f*SCLFAC
                     c = c*SCLFAC
                     ca = ca*SCLFAC
                     r = r/SCLFAC
                     g = g/SCLFAC
                     ra = ra/SCLFAC
                  ENDDO
!
                  g = c/SCLFAC
                  DO WHILE ( g>=r .AND. MAX(r,ra)<sfmax2 .AND.          &
     &                       MIN(f,c,g,ca)>sfmin2 )
                     f = f/SCLFAC
                     c = c/SCLFAC
                     g = g/SCLFAC
                     ca = ca/SCLFAC
                     r = r*SCLFAC
                     ra = ra*SCLFAC
                  ENDDO
!
!        Now balance.
!
                  IF ( (c+r)<FACTOR*s ) THEN
                     IF ( f<ONE .AND. Scale(i)<ONE ) THEN
                        IF ( f*Scale(i)<=sfmin1 ) CYCLE
                     ENDIF
                     IF ( f>ONE .AND. Scale(i)>ONE ) THEN
                        IF ( Scale(i)>=sfmax1/f ) CYCLE
                     ENDIF
                     g = ONE/f
                     Scale(i) = Scale(i)*f
                     noconv = .TRUE.
!
                     CALL DSCAL(N-k+1,g,A(i,k),Lda)
                     CALL DSCAL(l,f,A(1,i),1)
                  ENDIF
               ENDIF
!
            ENDDO
!
            IF ( .NOT.(noconv) ) EXIT
         ENDDO
      ENDIF
!
 700  Ilo = k
      Ihi = l
!
!
!     End of DGEBAL
!
      END SUBROUTINE DGEBAL
