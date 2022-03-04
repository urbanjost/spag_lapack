!*==dsbgst.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DSBGST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSBGST + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbgst.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbgst.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbgst.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X,
!                          LDX, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO, VECT
!       INTEGER            INFO, KA, KB, LDAB, LDBB, LDX, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AB( LDAB, * ), BB( LDBB, * ), WORK( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSBGST reduces a real symmetric-definite banded generalized
!> eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y,
!> such that C has the same bandwidth as A.
!>
!> B must have been previously factorized as S**T*S by DPBSTF, using a
!> split Cholesky factorization. A is overwritten by C = X**T*A*X, where
!> X = S**(-1)*Q and Q is an orthogonal matrix chosen to preserve the
!> bandwidth of A.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] VECT
!> \verbatim
!>          VECT is CHARACTER*1
!>          = 'N':  do not form the transformation matrix X;
!>          = 'V':  form X.
!> \endverbatim
!>
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
!>          The order of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in] KA
!> \verbatim
!>          KA is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.
!> \endverbatim
!>
!> \param[in] KB
!> \verbatim
!>          KB is INTEGER
!>          The number of superdiagonals of the matrix B if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KA >= KB >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          On entry, the upper or lower triangle of the symmetric band
!>          matrix A, stored in the first ka+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).
!>
!>          On exit, the transformed matrix X**T*A*X, stored in the same
!>          format as A.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KA+1.
!> \endverbatim
!>
!> \param[in] BB
!> \verbatim
!>          BB is DOUBLE PRECISION array, dimension (LDBB,N)
!>          The banded factor S from the split Cholesky factorization of
!>          B, as returned by DPBSTF, stored in the first KB+1 rows of
!>          the array.
!> \endverbatim
!>
!> \param[in] LDBB
!> \verbatim
!>          LDBB is INTEGER
!>          The leading dimension of the array BB.  LDBB >= KB+1.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,N)
!>          If VECT = 'V', the n-by-n matrix X.
!>          If VECT = 'N', the array X is not referenced.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.
!>          LDX >= max(1,N) if VECT = 'V'; LDX >= 1 otherwise.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
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
!> \date December 2016
!
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DSBGST(Vect,Uplo,N,Ka,Kb,Ab,Ldab,Bb,Ldbb,X,Ldx,Work,   &
     &                  Info)
      USE F77KINDS                        
      USE S_DGER
      USE S_DLAR2V
      USE S_DLARGV
      USE S_DLARTG
      USE S_DLARTV
      USE S_DLASET
      USE S_DROT
      USE S_DSCAL
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DSBGST174
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ka
      INTEGER , INTENT(IN) :: Kb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , DIMENSION(Ldbb,*) :: Bb
      INTEGER , INTENT(IN) :: Ldbb
      REAL(R8KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: bii , ra , ra1 , t
      INTEGER :: i , i0 , i1 , i2 , inca , j , j1 , j1t , j2 , j2t , k ,&
     &           ka1 , kb1 , kbt , l , m , nr , nrt , nx
      LOGICAL :: update , upper , wantx
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
!     Test the input parameters
!
      wantx = LSAME(Vect,'V')
      upper = LSAME(Uplo,'U')
      ka1 = Ka + 1
      kb1 = Kb + 1
      Info = 0
      IF ( .NOT.wantx .AND. .NOT.LSAME(Vect,'N') ) THEN
         Info = -1
      ELSEIF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Ka<0 ) THEN
         Info = -4
      ELSEIF ( Kb<0 .OR. Kb>Ka ) THEN
         Info = -5
      ELSEIF ( Ldab<Ka+1 ) THEN
         Info = -7
      ELSEIF ( Ldbb<Kb+1 ) THEN
         Info = -9
      ELSEIF ( Ldx<1 .OR. wantx .AND. Ldx<MAX(1,N) ) THEN
         Info = -11
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSBGST',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      inca = Ldab*ka1
!
!     Initialize X to the unit matrix, if needed
!
      IF ( wantx ) CALL DLASET('Full',N,N,ZERO,ONE,X,Ldx)
!
!     Set M to the splitting point m. It must be the same value as is
!     used in DPBSTF. The chosen value allows the arrays WORK and RWORK
!     to be of dimension (N).
!
      m = (N+Kb)/2
!
!     The routine works in two phases, corresponding to the two halves
!     of the split Cholesky factorization of B as S**T*S where
!
!     S = ( U    )
!         ( M  L )
!
!     with U upper triangular of order m, and L lower triangular of
!     order n-m. S has the same bandwidth as B.
!
!     S is treated as a product of elementary matrices:
!
!     S = S(m)*S(m-1)*...*S(2)*S(1)*S(m+1)*S(m+2)*...*S(n-1)*S(n)
!
!     where S(i) is determined by the i-th row of S.
!
!     In phase 1, the index i takes the values n, n-1, ... , m+1;
!     in phase 2, it takes the values 1, 2, ... , m.
!
!     For each value of i, the current matrix A is updated by forming
!     inv(S(i))**T*A*inv(S(i)). This creates a triangular bulge outside
!     the band of A. The bulge is then pushed down toward the bottom of
!     A in phase 1, and up toward the top of A in phase 2, by applying
!     plane rotations.
!
!     There are kb*(kb+1)/2 elements in the bulge, but at most 2*kb-1
!     of them are linearly independent, so annihilating a bulge requires
!     only 2*kb-1 plane rotations. The rotations are divided into a 1st
!     set of kb-1 rotations, and a 2nd set of kb rotations.
!
!     Wherever possible, rotations are generated and applied in vector
!     operations of length NR between the indices J1 and J2 (sometimes
!     replaced by modified values NRT, J1T or J2T).
!
!     The cosines and sines of the rotations are stored in the array
!     WORK. The cosines of the 1st set of rotations are stored in
!     elements n+2:n+m-kb-1 and the sines of the 1st set in elements
!     2:m-kb-1; the cosines of the 2nd set are stored in elements
!     n+m-kb+1:2*n and the sines of the second set in elements m-kb+1:n.
!
!     The bulges are not formed explicitly; nonzero elements outside the
!     band are created only when they are required for generating new
!     rotations; they are stored in the array WORK, in positions where
!     they are later overwritten by the sines of the rotations which
!     annihilate them.
!
!     **************************** Phase 1 *****************************
!
!     The logical structure of this phase is:
!
!     UPDATE = .TRUE.
!     DO I = N, M + 1, -1
!        use S(i) to update A and create a new bulge
!        apply rotations to push all bulges KA positions downward
!     END DO
!     UPDATE = .FALSE.
!     DO I = M + KA + 1, N - 1
!        apply rotations to push all bulges KA positions downward
!     END DO
!
!     To avoid duplicating code, the two loops are merged.
!
      update = .TRUE.
      i = N + 1
 100  DO WHILE ( update )
         i = i - 1
         kbt = MIN(Kb,i-1)
         i0 = i - 1
         i1 = MIN(N,i+Ka)
         i2 = i - kbt + ka1
         IF ( i<m+1 ) THEN
            update = .FALSE.
            i = i + 1
            i0 = m
            IF ( Ka/=0 ) CYCLE
            GOTO 300
         ENDIF
         GOTO 200
      ENDDO
      i = i + Ka
      IF ( i>N-1 ) GOTO 300
!
 200  IF ( upper ) THEN
!
!        Transform A, working with the upper triangle
!
         IF ( update ) THEN
!
!           Form  inv(S(i))**T * A * inv(S(i))
!
            bii = Bb(kb1,i)
            DO j = i , i1
               Ab(i-j+ka1,j) = Ab(i-j+ka1,j)/bii
            ENDDO
            DO j = MAX(1,i-Ka) , i
               Ab(j-i+ka1,i) = Ab(j-i+ka1,i)/bii
            ENDDO
            DO k = i - kbt , i - 1
               DO j = i - kbt , k
                  Ab(j-k+ka1,k) = Ab(j-k+ka1,k) - Bb(j-i+kb1,i)         &
     &                            *Ab(k-i+ka1,i) - Bb(k-i+kb1,i)        &
     &                            *Ab(j-i+ka1,i) + Ab(ka1,i)            &
     &                            *Bb(j-i+kb1,i)*Bb(k-i+kb1,i)
               ENDDO
               DO j = MAX(1,i-Ka) , i - kbt - 1
                  Ab(j-k+ka1,k) = Ab(j-k+ka1,k) - Bb(k-i+kb1,i)         &
     &                            *Ab(j-i+ka1,i)
               ENDDO
            ENDDO
            DO j = i , i1
               DO k = MAX(j-Ka,i-kbt) , i - 1
                  Ab(k-j+ka1,j) = Ab(k-j+ka1,j) - Bb(k-i+kb1,i)         &
     &                            *Ab(i-j+ka1,j)
               ENDDO
            ENDDO
!
            IF ( wantx ) THEN
!
!              post-multiply X by inv(S(i))
!
               CALL DSCAL(N-m,ONE/bii,X(m+1,i),1)
               IF ( kbt>0 ) CALL DGER(N-m,kbt,-ONE,X(m+1,i),1,Bb(kb1-kbt&
     &                                ,i),1,X(m+1,i-kbt),Ldx)
            ENDIF
!
!           store a(i,i1) in RA1 for use in next loop over K
!
            ra1 = Ab(i-i1+ka1,i1)
         ENDIF
!
!        Generate and apply vectors of rotations to chase all the
!        existing bulges KA positions down toward the bottom of the
!        band
!
         DO k = 1 , Kb - 1
            IF ( update ) THEN
!
!              Determine the rotations which would annihilate the bulge
!              which has in theory just been created
!
               IF ( i-k+Ka<N .AND. i-k>1 ) THEN
!
!                 generate rotation to annihilate a(i,i-k+ka+1)
!
                  CALL DLARTG(Ab(k+1,i-k+Ka),ra1,Work(N+i-k+Ka-m),      &
     &                        Work(i-k+Ka-m),ra)
!
!                 create nonzero element a(i-k,i-k+ka+1) outside the
!                 band and store it in WORK(i-k)
!
                  t = -Bb(kb1-k,i)*ra1
                  Work(i-k) = Work(N+i-k+Ka-m)*t - Work(i-k+Ka-m)       &
     &                        *Ab(1,i-k+Ka)
                  Ab(1,i-k+Ka) = Work(i-k+Ka-m)*t + Work(N+i-k+Ka-m)    &
     &                           *Ab(1,i-k+Ka)
                  ra1 = ra
               ENDIF
            ENDIF
            j2 = i - k - 1 + MAX(1,k-i0+2)*ka1
            nr = (N-j2+Ka)/ka1
            j1 = j2 + (nr-1)*ka1
            IF ( update ) THEN
               j2t = MAX(j2,i+2*Ka-k+1)
            ELSE
               j2t = j2
            ENDIF
            nrt = (N-j2t+Ka)/ka1
            DO j = j2t , j1 , ka1
!
!              create nonzero element a(j-ka,j+1) outside the band
!              and store it in WORK(j-m)
!
               Work(j-m) = Work(j-m)*Ab(1,j+1)
               Ab(1,j+1) = Work(N+j-m)*Ab(1,j+1)
            ENDDO
!
!           generate rotations in 1st set to annihilate elements which
!           have been created outside the band
!
            IF ( nrt>0 ) CALL DLARGV(nrt,Ab(1,j2t),inca,Work(j2t-m),ka1,&
     &                               Work(N+j2t-m),ka1)
            IF ( nr>0 ) THEN
!
!              apply rotations in 1st set from the right
!
               DO l = 1 , Ka - 1
                  CALL DLARTV(nr,Ab(ka1-l,j2),inca,Ab(Ka-l,j2+1),inca,  &
     &                        Work(N+j2-m),Work(j2-m),ka1)
               ENDDO
!
!              apply rotations in 1st set from both sides to diagonal
!              blocks
!
               CALL DLAR2V(nr,Ab(ka1,j2),Ab(ka1,j2+1),Ab(Ka,j2+1),inca, &
     &                     Work(N+j2-m),Work(j2-m),ka1)
!
            ENDIF
!
!           start applying rotations in 1st set from the left
!
            DO l = Ka - 1 , Kb - k + 1 , -1
               nrt = (N-j2+l)/ka1
               IF ( nrt>0 ) CALL DLARTV(nrt,Ab(l,j2+ka1-l),inca,        &
     &                                  Ab(l+1,j2+ka1-l),inca,          &
     &                                  Work(N+j2-m),Work(j2-m),ka1)
            ENDDO
!
            IF ( wantx ) THEN
!
!              post-multiply X by product of rotations in 1st set
!
               DO j = j2 , j1 , ka1
                  CALL DROT(N-m,X(m+1,j),1,X(m+1,j+1),1,Work(N+j-m),    &
     &                      Work(j-m))
               ENDDO
            ENDIF
         ENDDO
!
         IF ( update ) THEN
!
!              create nonzero element a(i-kbt,i-kbt+ka+1) outside the
!              band and store it in WORK(i-kbt)
!
            IF ( i2<=N .AND. kbt>0 ) Work(i-kbt) = -Bb(kb1-kbt,i)*ra1
         ENDIF
!
         DO k = Kb , 1 , -1
            IF ( update ) THEN
               j2 = i - k - 1 + MAX(2,k-i0+1)*ka1
            ELSE
               j2 = i - k - 1 + MAX(1,k-i0+1)*ka1
            ENDIF
!
!           finish applying rotations in 2nd set from the left
!
            DO l = Kb - k , 1 , -1
               nrt = (N-j2+Ka+l)/ka1
               IF ( nrt>0 ) CALL DLARTV(nrt,Ab(l,j2-l+1),inca,Ab(l+1,j2-&
     &                                  l+1),inca,Work(N+j2-Ka),        &
     &                                  Work(j2-Ka),ka1)
            ENDDO
            nr = (N-j2+Ka)/ka1
            j1 = j2 + (nr-1)*ka1
            DO j = j1 , j2 , -ka1
               Work(j) = Work(j-Ka)
               Work(N+j) = Work(N+j-Ka)
            ENDDO
            DO j = j2 , j1 , ka1
!
!              create nonzero element a(j-ka,j+1) outside the band
!              and store it in WORK(j)
!
               Work(j) = Work(j)*Ab(1,j+1)
               Ab(1,j+1) = Work(N+j)*Ab(1,j+1)
            ENDDO
            IF ( update ) THEN
               IF ( i-k<N-Ka .AND. k<=kbt ) Work(i-k+Ka) = Work(i-k)
            ENDIF
         ENDDO
!
         DO k = Kb , 1 , -1
            j2 = i - k - 1 + MAX(1,k-i0+1)*ka1
            nr = (N-j2+Ka)/ka1
            j1 = j2 + (nr-1)*ka1
            IF ( nr>0 ) THEN
!
!              generate rotations in 2nd set to annihilate elements
!              which have been created outside the band
!
               CALL DLARGV(nr,Ab(1,j2),inca,Work(j2),ka1,Work(N+j2),ka1)
!
!              apply rotations in 2nd set from the right
!
               DO l = 1 , Ka - 1
                  CALL DLARTV(nr,Ab(ka1-l,j2),inca,Ab(Ka-l,j2+1),inca,  &
     &                        Work(N+j2),Work(j2),ka1)
               ENDDO
!
!              apply rotations in 2nd set from both sides to diagonal
!              blocks
!
               CALL DLAR2V(nr,Ab(ka1,j2),Ab(ka1,j2+1),Ab(Ka,j2+1),inca, &
     &                     Work(N+j2),Work(j2),ka1)
!
            ENDIF
!
!           start applying rotations in 2nd set from the left
!
            DO l = Ka - 1 , Kb - k + 1 , -1
               nrt = (N-j2+l)/ka1
               IF ( nrt>0 ) CALL DLARTV(nrt,Ab(l,j2+ka1-l),inca,        &
     &                                  Ab(l+1,j2+ka1-l),inca,Work(N+j2)&
     &                                  ,Work(j2),ka1)
            ENDDO
!
            IF ( wantx ) THEN
!
!              post-multiply X by product of rotations in 2nd set
!
               DO j = j2 , j1 , ka1
                  CALL DROT(N-m,X(m+1,j),1,X(m+1,j+1),1,Work(N+j),      &
     &                      Work(j))
               ENDDO
            ENDIF
         ENDDO
!
         DO k = 1 , Kb - 1
            j2 = i - k - 1 + MAX(1,k-i0+2)*ka1
!
!           finish applying rotations in 1st set from the left
!
            DO l = Kb - k , 1 , -1
               nrt = (N-j2+l)/ka1
               IF ( nrt>0 ) CALL DLARTV(nrt,Ab(l,j2+ka1-l),inca,        &
     &                                  Ab(l+1,j2+ka1-l),inca,          &
     &                                  Work(N+j2-m),Work(j2-m),ka1)
            ENDDO
         ENDDO
!
         IF ( Kb>1 ) THEN
            DO j = N - 1 , i - Kb + 2*Ka + 1 , -1
               Work(N+j-m) = Work(N+j-Ka-m)
               Work(j-m) = Work(j-Ka-m)
            ENDDO
         ENDIF
!
      ELSE
!
!        Transform A, working with the lower triangle
!
         IF ( update ) THEN
!
!           Form  inv(S(i))**T * A * inv(S(i))
!
            bii = Bb(1,i)
            DO j = i , i1
               Ab(j-i+1,i) = Ab(j-i+1,i)/bii
            ENDDO
            DO j = MAX(1,i-Ka) , i
               Ab(i-j+1,j) = Ab(i-j+1,j)/bii
            ENDDO
            DO k = i - kbt , i - 1
               DO j = i - kbt , k
                  Ab(k-j+1,j) = Ab(k-j+1,j) - Bb(i-j+1,j)*Ab(i-k+1,k)   &
     &                          - Bb(i-k+1,k)*Ab(i-j+1,j) + Ab(1,i)     &
     &                          *Bb(i-j+1,j)*Bb(i-k+1,k)
               ENDDO
               DO j = MAX(1,i-Ka) , i - kbt - 1
                  Ab(k-j+1,j) = Ab(k-j+1,j) - Bb(i-k+1,k)*Ab(i-j+1,j)
               ENDDO
            ENDDO
            DO j = i , i1
               DO k = MAX(j-Ka,i-kbt) , i - 1
                  Ab(j-k+1,k) = Ab(j-k+1,k) - Bb(i-k+1,k)*Ab(j-i+1,i)
               ENDDO
            ENDDO
!
            IF ( wantx ) THEN
!
!              post-multiply X by inv(S(i))
!
               CALL DSCAL(N-m,ONE/bii,X(m+1,i),1)
               IF ( kbt>0 ) CALL DGER(N-m,kbt,-ONE,X(m+1,i),1,Bb(kbt+1,i&
     &                                -kbt),Ldbb-1,X(m+1,i-kbt),Ldx)
            ENDIF
!
!           store a(i1,i) in RA1 for use in next loop over K
!
            ra1 = Ab(i1-i+1,i)
         ENDIF
!
!        Generate and apply vectors of rotations to chase all the
!        existing bulges KA positions down toward the bottom of the
!        band
!
         DO k = 1 , Kb - 1
            IF ( update ) THEN
!
!              Determine the rotations which would annihilate the bulge
!              which has in theory just been created
!
               IF ( i-k+Ka<N .AND. i-k>1 ) THEN
!
!                 generate rotation to annihilate a(i-k+ka+1,i)
!
                  CALL DLARTG(Ab(ka1-k,i),ra1,Work(N+i-k+Ka-m),         &
     &                        Work(i-k+Ka-m),ra)
!
!                 create nonzero element a(i-k+ka+1,i-k) outside the
!                 band and store it in WORK(i-k)
!
                  t = -Bb(k+1,i-k)*ra1
                  Work(i-k) = Work(N+i-k+Ka-m)*t - Work(i-k+Ka-m)       &
     &                        *Ab(ka1,i-k)
                  Ab(ka1,i-k) = Work(i-k+Ka-m)*t + Work(N+i-k+Ka-m)     &
     &                          *Ab(ka1,i-k)
                  ra1 = ra
               ENDIF
            ENDIF
            j2 = i - k - 1 + MAX(1,k-i0+2)*ka1
            nr = (N-j2+Ka)/ka1
            j1 = j2 + (nr-1)*ka1
            IF ( update ) THEN
               j2t = MAX(j2,i+2*Ka-k+1)
            ELSE
               j2t = j2
            ENDIF
            nrt = (N-j2t+Ka)/ka1
            DO j = j2t , j1 , ka1
!
!              create nonzero element a(j+1,j-ka) outside the band
!              and store it in WORK(j-m)
!
               Work(j-m) = Work(j-m)*Ab(ka1,j-Ka+1)
               Ab(ka1,j-Ka+1) = Work(N+j-m)*Ab(ka1,j-Ka+1)
            ENDDO
!
!           generate rotations in 1st set to annihilate elements which
!           have been created outside the band
!
            IF ( nrt>0 ) CALL DLARGV(nrt,Ab(ka1,j2t-Ka),inca,Work(j2t-m)&
     &                               ,ka1,Work(N+j2t-m),ka1)
            IF ( nr>0 ) THEN
!
!              apply rotations in 1st set from the left
!
               DO l = 1 , Ka - 1
                  CALL DLARTV(nr,Ab(l+1,j2-l),inca,Ab(l+2,j2-l),inca,   &
     &                        Work(N+j2-m),Work(j2-m),ka1)
               ENDDO
!
!              apply rotations in 1st set from both sides to diagonal
!              blocks
!
               CALL DLAR2V(nr,Ab(1,j2),Ab(1,j2+1),Ab(2,j2),inca,        &
     &                     Work(N+j2-m),Work(j2-m),ka1)
!
            ENDIF
!
!           start applying rotations in 1st set from the right
!
            DO l = Ka - 1 , Kb - k + 1 , -1
               nrt = (N-j2+l)/ka1
               IF ( nrt>0 ) CALL DLARTV(nrt,Ab(ka1-l+1,j2),inca,        &
     &                                  Ab(ka1-l,j2+1),inca,Work(N+j2-m)&
     &                                  ,Work(j2-m),ka1)
            ENDDO
!
            IF ( wantx ) THEN
!
!              post-multiply X by product of rotations in 1st set
!
               DO j = j2 , j1 , ka1
                  CALL DROT(N-m,X(m+1,j),1,X(m+1,j+1),1,Work(N+j-m),    &
     &                      Work(j-m))
               ENDDO
            ENDIF
         ENDDO
!
         IF ( update ) THEN
!
!              create nonzero element a(i-kbt+ka+1,i-kbt) outside the
!              band and store it in WORK(i-kbt)
!
            IF ( i2<=N .AND. kbt>0 ) Work(i-kbt) = -Bb(kbt+1,i-kbt)*ra1
         ENDIF
!
         DO k = Kb , 1 , -1
            IF ( update ) THEN
               j2 = i - k - 1 + MAX(2,k-i0+1)*ka1
            ELSE
               j2 = i - k - 1 + MAX(1,k-i0+1)*ka1
            ENDIF
!
!           finish applying rotations in 2nd set from the right
!
            DO l = Kb - k , 1 , -1
               nrt = (N-j2+Ka+l)/ka1
               IF ( nrt>0 ) CALL DLARTV(nrt,Ab(ka1-l+1,j2-Ka),inca,     &
     &                                  Ab(ka1-l,j2-Ka+1),inca,         &
     &                                  Work(N+j2-Ka),Work(j2-Ka),ka1)
            ENDDO
            nr = (N-j2+Ka)/ka1
            j1 = j2 + (nr-1)*ka1
            DO j = j1 , j2 , -ka1
               Work(j) = Work(j-Ka)
               Work(N+j) = Work(N+j-Ka)
            ENDDO
            DO j = j2 , j1 , ka1
!
!              create nonzero element a(j+1,j-ka) outside the band
!              and store it in WORK(j)
!
               Work(j) = Work(j)*Ab(ka1,j-Ka+1)
               Ab(ka1,j-Ka+1) = Work(N+j)*Ab(ka1,j-Ka+1)
            ENDDO
            IF ( update ) THEN
               IF ( i-k<N-Ka .AND. k<=kbt ) Work(i-k+Ka) = Work(i-k)
            ENDIF
         ENDDO
!
         DO k = Kb , 1 , -1
            j2 = i - k - 1 + MAX(1,k-i0+1)*ka1
            nr = (N-j2+Ka)/ka1
            j1 = j2 + (nr-1)*ka1
            IF ( nr>0 ) THEN
!
!              generate rotations in 2nd set to annihilate elements
!              which have been created outside the band
!
               CALL DLARGV(nr,Ab(ka1,j2-Ka),inca,Work(j2),ka1,Work(N+j2)&
     &                     ,ka1)
!
!              apply rotations in 2nd set from the left
!
               DO l = 1 , Ka - 1
                  CALL DLARTV(nr,Ab(l+1,j2-l),inca,Ab(l+2,j2-l),inca,   &
     &                        Work(N+j2),Work(j2),ka1)
               ENDDO
!
!              apply rotations in 2nd set from both sides to diagonal
!              blocks
!
               CALL DLAR2V(nr,Ab(1,j2),Ab(1,j2+1),Ab(2,j2),inca,        &
     &                     Work(N+j2),Work(j2),ka1)
!
            ENDIF
!
!           start applying rotations in 2nd set from the right
!
            DO l = Ka - 1 , Kb - k + 1 , -1
               nrt = (N-j2+l)/ka1
               IF ( nrt>0 ) CALL DLARTV(nrt,Ab(ka1-l+1,j2),inca,        &
     &                                  Ab(ka1-l,j2+1),inca,Work(N+j2), &
     &                                  Work(j2),ka1)
            ENDDO
!
            IF ( wantx ) THEN
!
!              post-multiply X by product of rotations in 2nd set
!
               DO j = j2 , j1 , ka1
                  CALL DROT(N-m,X(m+1,j),1,X(m+1,j+1),1,Work(N+j),      &
     &                      Work(j))
               ENDDO
            ENDIF
         ENDDO
!
         DO k = 1 , Kb - 1
            j2 = i - k - 1 + MAX(1,k-i0+2)*ka1
!
!           finish applying rotations in 1st set from the right
!
            DO l = Kb - k , 1 , -1
               nrt = (N-j2+l)/ka1
               IF ( nrt>0 ) CALL DLARTV(nrt,Ab(ka1-l+1,j2),inca,        &
     &                                  Ab(ka1-l,j2+1),inca,Work(N+j2-m)&
     &                                  ,Work(j2-m),ka1)
            ENDDO
         ENDDO
!
         IF ( Kb>1 ) THEN
            DO j = N - 1 , i - Kb + 2*Ka + 1 , -1
               Work(N+j-m) = Work(N+j-Ka-m)
               Work(j-m) = Work(j-Ka-m)
            ENDDO
         ENDIF
!
      ENDIF
!
      GOTO 100
!
!
!     **************************** Phase 2 *****************************
!
!     The logical structure of this phase is:
!
!     UPDATE = .TRUE.
!     DO I = 1, M
!        use S(i) to update A and create a new bulge
!        apply rotations to push all bulges KA positions upward
!     END DO
!     UPDATE = .FALSE.
!     DO I = M - KA - 1, 2, -1
!        apply rotations to push all bulges KA positions upward
!     END DO
!
!     To avoid duplicating code, the two loops are merged.
!
 300  update = .TRUE.
      i = 0
      DO
         IF ( update ) THEN
            i = i + 1
            kbt = MIN(Kb,m-i)
            i0 = i + 1
            i1 = MAX(1,i-Ka)
            i2 = i + kbt - ka1
            IF ( i>m ) THEN
               update = .FALSE.
               i = i - 1
               i0 = m + 1
               IF ( Ka==0 ) RETURN
               CYCLE
            ENDIF
         ELSE
            i = i - Ka
            IF ( i<2 ) RETURN
         ENDIF
!
         IF ( i<m-kbt ) THEN
            nx = m
         ELSE
            nx = N
         ENDIF
!
         IF ( upper ) THEN
!
!        Transform A, working with the upper triangle
!
            IF ( update ) THEN
!
!           Form  inv(S(i))**T * A * inv(S(i))
!
               bii = Bb(kb1,i)
               DO j = i1 , i
                  Ab(j-i+ka1,i) = Ab(j-i+ka1,i)/bii
               ENDDO
               DO j = i , MIN(N,i+Ka)
                  Ab(i-j+ka1,j) = Ab(i-j+ka1,j)/bii
               ENDDO
               DO k = i + 1 , i + kbt
                  DO j = k , i + kbt
                     Ab(k-j+ka1,j) = Ab(k-j+ka1,j) - Bb(i-j+kb1,j)      &
     &                               *Ab(i-k+ka1,k) - Bb(i-k+kb1,k)     &
     &                               *Ab(i-j+ka1,j) + Ab(ka1,i)         &
     &                               *Bb(i-j+kb1,j)*Bb(i-k+kb1,k)
                  ENDDO
                  DO j = i + kbt + 1 , MIN(N,i+Ka)
                     Ab(k-j+ka1,j) = Ab(k-j+ka1,j) - Bb(i-k+kb1,k)      &
     &                               *Ab(i-j+ka1,j)
                  ENDDO
               ENDDO
               DO j = i1 , i
                  DO k = i + 1 , MIN(j+Ka,i+kbt)
                     Ab(j-k+ka1,k) = Ab(j-k+ka1,k) - Bb(i-k+kb1,k)      &
     &                               *Ab(j-i+ka1,i)
                  ENDDO
               ENDDO
!
               IF ( wantx ) THEN
!
!              post-multiply X by inv(S(i))
!
                  CALL DSCAL(nx,ONE/bii,X(1,i),1)
                  IF ( kbt>0 ) CALL DGER(nx,kbt,-ONE,X(1,i),1,Bb(Kb,i+1)&
     &                 ,Ldbb-1,X(1,i+1),Ldx)
               ENDIF
!
!           store a(i1,i) in RA1 for use in next loop over K
!
               ra1 = Ab(i1-i+ka1,i)
            ENDIF
!
!        Generate and apply vectors of rotations to chase all the
!        existing bulges KA positions up toward the top of the band
!
            DO k = 1 , Kb - 1
               IF ( update ) THEN
!
!              Determine the rotations which would annihilate the bulge
!              which has in theory just been created
!
                  IF ( i+k-ka1>0 .AND. i+k<m ) THEN
!
!                 generate rotation to annihilate a(i+k-ka-1,i)
!
                     CALL DLARTG(Ab(k+1,i),ra1,Work(N+i+k-Ka),          &
     &                           Work(i+k-Ka),ra)
!
!                 create nonzero element a(i+k-ka-1,i+k) outside the
!                 band and store it in WORK(m-kb+i+k)
!
                     t = -Bb(kb1-k,i+k)*ra1
                     Work(m-Kb+i+k) = Work(N+i+k-Ka)*t - Work(i+k-Ka)   &
     &                                *Ab(1,i+k)
                     Ab(1,i+k) = Work(i+k-Ka)*t + Work(N+i+k-Ka)        &
     &                           *Ab(1,i+k)
                     ra1 = ra
                  ENDIF
               ENDIF
               j2 = i + k + 1 - MAX(1,k+i0-m+1)*ka1
               nr = (j2+Ka-1)/ka1
               j1 = j2 - (nr-1)*ka1
               IF ( update ) THEN
                  j2t = MIN(j2,i-2*Ka+k-1)
               ELSE
                  j2t = j2
               ENDIF
               nrt = (j2t+Ka-1)/ka1
               DO j = j1 , j2t , ka1
!
!              create nonzero element a(j-1,j+ka) outside the band
!              and store it in WORK(j)
!
                  Work(j) = Work(j)*Ab(1,j+Ka-1)
                  Ab(1,j+Ka-1) = Work(N+j)*Ab(1,j+Ka-1)
               ENDDO
!
!           generate rotations in 1st set to annihilate elements which
!           have been created outside the band
!
               IF ( nrt>0 ) CALL DLARGV(nrt,Ab(1,j1+Ka),inca,Work(j1),  &
     &                                  ka1,Work(N+j1),ka1)
               IF ( nr>0 ) THEN
!
!              apply rotations in 1st set from the left
!
                  DO l = 1 , Ka - 1
                     CALL DLARTV(nr,Ab(ka1-l,j1+l),inca,Ab(Ka-l,j1+l),  &
     &                           inca,Work(N+j1),Work(j1),ka1)
                  ENDDO
!
!              apply rotations in 1st set from both sides to diagonal
!              blocks
!
                  CALL DLAR2V(nr,Ab(ka1,j1),Ab(ka1,j1-1),Ab(Ka,j1),inca,&
     &                        Work(N+j1),Work(j1),ka1)
!
               ENDIF
!
!           start applying rotations in 1st set from the right
!
               DO l = Ka - 1 , Kb - k + 1 , -1
                  nrt = (j2+l-1)/ka1
                  j1t = j2 - (nrt-1)*ka1
                  IF ( nrt>0 )                                          &
     &                 CALL DLARTV(nrt,Ab(l,j1t),inca,Ab(l+1,j1t-1),    &
     &                 inca,Work(N+j1t),Work(j1t),ka1)
               ENDDO
!
               IF ( wantx ) THEN
!
!              post-multiply X by product of rotations in 1st set
!
                  DO j = j1 , j2 , ka1
                     CALL DROT(nx,X(1,j),1,X(1,j-1),1,Work(N+j),Work(j))
                  ENDDO
               ENDIF
            ENDDO
!
            IF ( update ) THEN
!
!              create nonzero element a(i+kbt-ka-1,i+kbt) outside the
!              band and store it in WORK(m-kb+i+kbt)
!
               IF ( i2>0 .AND. kbt>0 ) Work(m-Kb+i+kbt)                 &
     &              = -Bb(kb1-kbt,i+kbt)*ra1
            ENDIF
!
            DO k = Kb , 1 , -1
               IF ( update ) THEN
                  j2 = i + k + 1 - MAX(2,k+i0-m)*ka1
               ELSE
                  j2 = i + k + 1 - MAX(1,k+i0-m)*ka1
               ENDIF
!
!           finish applying rotations in 2nd set from the right
!
               DO l = Kb - k , 1 , -1
                  nrt = (j2+Ka+l-1)/ka1
                  j1t = j2 - (nrt-1)*ka1
                  IF ( nrt>0 )                                          &
     &                 CALL DLARTV(nrt,Ab(l,j1t+Ka),inca,Ab(l+1,j1t+Ka- &
     &                 1),inca,Work(N+m-Kb+j1t+Ka),Work(m-Kb+j1t+Ka),   &
     &                 ka1)
               ENDDO
               nr = (j2+Ka-1)/ka1
               j1 = j2 - (nr-1)*ka1
               DO j = j1 , j2 , ka1
                  Work(m-Kb+j) = Work(m-Kb+j+Ka)
                  Work(N+m-Kb+j) = Work(N+m-Kb+j+Ka)
               ENDDO
               DO j = j1 , j2 , ka1
!
!              create nonzero element a(j-1,j+ka) outside the band
!              and store it in WORK(m-kb+j)
!
                  Work(m-Kb+j) = Work(m-Kb+j)*Ab(1,j+Ka-1)
                  Ab(1,j+Ka-1) = Work(N+m-Kb+j)*Ab(1,j+Ka-1)
               ENDDO
               IF ( update ) THEN
                  IF ( i+k>ka1 .AND. k<=kbt ) Work(m-Kb+i+k-Ka)         &
     &                 = Work(m-Kb+i+k)
               ENDIF
            ENDDO
!
            DO k = Kb , 1 , -1
               j2 = i + k + 1 - MAX(1,k+i0-m)*ka1
               nr = (j2+Ka-1)/ka1
               j1 = j2 - (nr-1)*ka1
               IF ( nr>0 ) THEN
!
!              generate rotations in 2nd set to annihilate elements
!              which have been created outside the band
!
                  CALL DLARGV(nr,Ab(1,j1+Ka),inca,Work(m-Kb+j1),ka1,    &
     &                        Work(N+m-Kb+j1),ka1)
!
!              apply rotations in 2nd set from the left
!
                  DO l = 1 , Ka - 1
                     CALL DLARTV(nr,Ab(ka1-l,j1+l),inca,Ab(Ka-l,j1+l),  &
     &                           inca,Work(N+m-Kb+j1),Work(m-Kb+j1),ka1)
                  ENDDO
!
!              apply rotations in 2nd set from both sides to diagonal
!              blocks
!
                  CALL DLAR2V(nr,Ab(ka1,j1),Ab(ka1,j1-1),Ab(Ka,j1),inca,&
     &                        Work(N+m-Kb+j1),Work(m-Kb+j1),ka1)
!
               ENDIF
!
!           start applying rotations in 2nd set from the right
!
               DO l = Ka - 1 , Kb - k + 1 , -1
                  nrt = (j2+l-1)/ka1
                  j1t = j2 - (nrt-1)*ka1
                  IF ( nrt>0 )                                          &
     &                 CALL DLARTV(nrt,Ab(l,j1t),inca,Ab(l+1,j1t-1),    &
     &                 inca,Work(N+m-Kb+j1t),Work(m-Kb+j1t),ka1)
               ENDDO
!
               IF ( wantx ) THEN
!
!              post-multiply X by product of rotations in 2nd set
!
                  DO j = j1 , j2 , ka1
                     CALL DROT(nx,X(1,j),1,X(1,j-1),1,Work(N+m-Kb+j),   &
     &                         Work(m-Kb+j))
                  ENDDO
               ENDIF
            ENDDO
!
            DO k = 1 , Kb - 1
               j2 = i + k + 1 - MAX(1,k+i0-m+1)*ka1
!
!           finish applying rotations in 1st set from the right
!
               DO l = Kb - k , 1 , -1
                  nrt = (j2+l-1)/ka1
                  j1t = j2 - (nrt-1)*ka1
                  IF ( nrt>0 )                                          &
     &                 CALL DLARTV(nrt,Ab(l,j1t),inca,Ab(l+1,j1t-1),    &
     &                 inca,Work(N+j1t),Work(j1t),ka1)
               ENDDO
            ENDDO
!
            IF ( Kb>1 ) THEN
               DO j = 2 , MIN(i+Kb,m) - 2*Ka - 1
                  Work(N+j) = Work(N+j+Ka)
                  Work(j) = Work(j+Ka)
               ENDDO
            ENDIF
!
         ELSE
!
!        Transform A, working with the lower triangle
!
            IF ( update ) THEN
!
!           Form  inv(S(i))**T * A * inv(S(i))
!
               bii = Bb(1,i)
               DO j = i1 , i
                  Ab(i-j+1,j) = Ab(i-j+1,j)/bii
               ENDDO
               DO j = i , MIN(N,i+Ka)
                  Ab(j-i+1,i) = Ab(j-i+1,i)/bii
               ENDDO
               DO k = i + 1 , i + kbt
                  DO j = k , i + kbt
                     Ab(j-k+1,k) = Ab(j-k+1,k) - Bb(j-i+1,i)*Ab(k-i+1,i)&
     &                             - Bb(k-i+1,i)*Ab(j-i+1,i) + Ab(1,i)  &
     &                             *Bb(j-i+1,i)*Bb(k-i+1,i)
                  ENDDO
                  DO j = i + kbt + 1 , MIN(N,i+Ka)
                     Ab(j-k+1,k) = Ab(j-k+1,k) - Bb(k-i+1,i)*Ab(j-i+1,i)
                  ENDDO
               ENDDO
               DO j = i1 , i
                  DO k = i + 1 , MIN(j+Ka,i+kbt)
                     Ab(k-j+1,j) = Ab(k-j+1,j) - Bb(k-i+1,i)*Ab(i-j+1,j)
                  ENDDO
               ENDDO
!
               IF ( wantx ) THEN
!
!              post-multiply X by inv(S(i))
!
                  CALL DSCAL(nx,ONE/bii,X(1,i),1)
                  IF ( kbt>0 ) CALL DGER(nx,kbt,-ONE,X(1,i),1,Bb(2,i),1,&
     &                 X(1,i+1),Ldx)
               ENDIF
!
!           store a(i,i1) in RA1 for use in next loop over K
!
               ra1 = Ab(i-i1+1,i1)
            ENDIF
!
!        Generate and apply vectors of rotations to chase all the
!        existing bulges KA positions up toward the top of the band
!
            DO k = 1 , Kb - 1
               IF ( update ) THEN
!
!              Determine the rotations which would annihilate the bulge
!              which has in theory just been created
!
                  IF ( i+k-ka1>0 .AND. i+k<m ) THEN
!
!                 generate rotation to annihilate a(i,i+k-ka-1)
!
                     CALL DLARTG(Ab(ka1-k,i+k-Ka),ra1,Work(N+i+k-Ka),   &
     &                           Work(i+k-Ka),ra)
!
!                 create nonzero element a(i+k,i+k-ka-1) outside the
!                 band and store it in WORK(m-kb+i+k)
!
                     t = -Bb(k+1,i)*ra1
                     Work(m-Kb+i+k) = Work(N+i+k-Ka)*t - Work(i+k-Ka)   &
     &                                *Ab(ka1,i+k-Ka)
                     Ab(ka1,i+k-Ka) = Work(i+k-Ka)*t + Work(N+i+k-Ka)   &
     &                                *Ab(ka1,i+k-Ka)
                     ra1 = ra
                  ENDIF
               ENDIF
               j2 = i + k + 1 - MAX(1,k+i0-m+1)*ka1
               nr = (j2+Ka-1)/ka1
               j1 = j2 - (nr-1)*ka1
               IF ( update ) THEN
                  j2t = MIN(j2,i-2*Ka+k-1)
               ELSE
                  j2t = j2
               ENDIF
               nrt = (j2t+Ka-1)/ka1
               DO j = j1 , j2t , ka1
!
!              create nonzero element a(j+ka,j-1) outside the band
!              and store it in WORK(j)
!
                  Work(j) = Work(j)*Ab(ka1,j-1)
                  Ab(ka1,j-1) = Work(N+j)*Ab(ka1,j-1)
               ENDDO
!
!           generate rotations in 1st set to annihilate elements which
!           have been created outside the band
!
               IF ( nrt>0 ) CALL DLARGV(nrt,Ab(ka1,j1),inca,Work(j1),   &
     &                                  ka1,Work(N+j1),ka1)
               IF ( nr>0 ) THEN
!
!              apply rotations in 1st set from the right
!
                  DO l = 1 , Ka - 1
                     CALL DLARTV(nr,Ab(l+1,j1),inca,Ab(l+2,j1-1),inca,  &
     &                           Work(N+j1),Work(j1),ka1)
                  ENDDO
!
!              apply rotations in 1st set from both sides to diagonal
!              blocks
!
                  CALL DLAR2V(nr,Ab(1,j1),Ab(1,j1-1),Ab(2,j1-1),inca,   &
     &                        Work(N+j1),Work(j1),ka1)
!
               ENDIF
!
!           start applying rotations in 1st set from the left
!
               DO l = Ka - 1 , Kb - k + 1 , -1
                  nrt = (j2+l-1)/ka1
                  j1t = j2 - (nrt-1)*ka1
                  IF ( nrt>0 ) CALL DLARTV(nrt,Ab(ka1-l+1,j1t-ka1+l),   &
     &                 inca,Ab(ka1-l,j1t-ka1+l),inca,Work(N+j1t),       &
     &                 Work(j1t),ka1)
               ENDDO
!
               IF ( wantx ) THEN
!
!              post-multiply X by product of rotations in 1st set
!
                  DO j = j1 , j2 , ka1
                     CALL DROT(nx,X(1,j),1,X(1,j-1),1,Work(N+j),Work(j))
                  ENDDO
               ENDIF
            ENDDO
!
            IF ( update ) THEN
!
!              create nonzero element a(i+kbt,i+kbt-ka-1) outside the
!              band and store it in WORK(m-kb+i+kbt)
!
               IF ( i2>0 .AND. kbt>0 ) Work(m-Kb+i+kbt) = -Bb(kbt+1,i)  &
     &              *ra1
            ENDIF
!
            DO k = Kb , 1 , -1
               IF ( update ) THEN
                  j2 = i + k + 1 - MAX(2,k+i0-m)*ka1
               ELSE
                  j2 = i + k + 1 - MAX(1,k+i0-m)*ka1
               ENDIF
!
!           finish applying rotations in 2nd set from the left
!
               DO l = Kb - k , 1 , -1
                  nrt = (j2+Ka+l-1)/ka1
                  j1t = j2 - (nrt-1)*ka1
                  IF ( nrt>0 ) CALL DLARTV(nrt,Ab(ka1-l+1,j1t+l-1),inca,&
     &                 Ab(ka1-l,j1t+l-1),inca,Work(N+m-Kb+j1t+Ka),      &
     &                 Work(m-Kb+j1t+Ka),ka1)
               ENDDO
               nr = (j2+Ka-1)/ka1
               j1 = j2 - (nr-1)*ka1
               DO j = j1 , j2 , ka1
                  Work(m-Kb+j) = Work(m-Kb+j+Ka)
                  Work(N+m-Kb+j) = Work(N+m-Kb+j+Ka)
               ENDDO
               DO j = j1 , j2 , ka1
!
!              create nonzero element a(j+ka,j-1) outside the band
!              and store it in WORK(m-kb+j)
!
                  Work(m-Kb+j) = Work(m-Kb+j)*Ab(ka1,j-1)
                  Ab(ka1,j-1) = Work(N+m-Kb+j)*Ab(ka1,j-1)
               ENDDO
               IF ( update ) THEN
                  IF ( i+k>ka1 .AND. k<=kbt ) Work(m-Kb+i+k-Ka)         &
     &                 = Work(m-Kb+i+k)
               ENDIF
            ENDDO
!
            DO k = Kb , 1 , -1
               j2 = i + k + 1 - MAX(1,k+i0-m)*ka1
               nr = (j2+Ka-1)/ka1
               j1 = j2 - (nr-1)*ka1
               IF ( nr>0 ) THEN
!
!              generate rotations in 2nd set to annihilate elements
!              which have been created outside the band
!
                  CALL DLARGV(nr,Ab(ka1,j1),inca,Work(m-Kb+j1),ka1,     &
     &                        Work(N+m-Kb+j1),ka1)
!
!              apply rotations in 2nd set from the right
!
                  DO l = 1 , Ka - 1
                     CALL DLARTV(nr,Ab(l+1,j1),inca,Ab(l+2,j1-1),inca,  &
     &                           Work(N+m-Kb+j1),Work(m-Kb+j1),ka1)
                  ENDDO
!
!              apply rotations in 2nd set from both sides to diagonal
!              blocks
!
                  CALL DLAR2V(nr,Ab(1,j1),Ab(1,j1-1),Ab(2,j1-1),inca,   &
     &                        Work(N+m-Kb+j1),Work(m-Kb+j1),ka1)
!
               ENDIF
!
!           start applying rotations in 2nd set from the left
!
               DO l = Ka - 1 , Kb - k + 1 , -1
                  nrt = (j2+l-1)/ka1
                  j1t = j2 - (nrt-1)*ka1
                  IF ( nrt>0 ) CALL DLARTV(nrt,Ab(ka1-l+1,j1t-ka1+l),   &
     &                 inca,Ab(ka1-l,j1t-ka1+l),inca,Work(N+m-Kb+j1t),  &
     &                 Work(m-Kb+j1t),ka1)
               ENDDO
!
               IF ( wantx ) THEN
!
!              post-multiply X by product of rotations in 2nd set
!
                  DO j = j1 , j2 , ka1
                     CALL DROT(nx,X(1,j),1,X(1,j-1),1,Work(N+m-Kb+j),   &
     &                         Work(m-Kb+j))
                  ENDDO
               ENDIF
            ENDDO
!
            DO k = 1 , Kb - 1
               j2 = i + k + 1 - MAX(1,k+i0-m+1)*ka1
!
!           finish applying rotations in 1st set from the left
!
               DO l = Kb - k , 1 , -1
                  nrt = (j2+l-1)/ka1
                  j1t = j2 - (nrt-1)*ka1
                  IF ( nrt>0 ) CALL DLARTV(nrt,Ab(ka1-l+1,j1t-ka1+l),   &
     &                 inca,Ab(ka1-l,j1t-ka1+l),inca,Work(N+j1t),       &
     &                 Work(j1t),ka1)
               ENDDO
            ENDDO
!
            IF ( Kb>1 ) THEN
               DO j = 2 , MIN(i+Kb,m) - 2*Ka - 1
                  Work(N+j) = Work(N+j+Ka)
                  Work(j) = Work(j+Ka)
               ENDDO
            ENDIF
!
!
         ENDIF
      ENDDO
!
!     End of DSBGST
!
      END SUBROUTINE DSBGST
