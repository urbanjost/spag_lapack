!*==ssbtrd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SSBTRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSBTRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbtrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbtrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbtrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO, VECT
!       INTEGER            INFO, KD, LDAB, LDQ, N
!       ..
!       .. Array Arguments ..
!       REAL               AB( LDAB, * ), D( * ), E( * ), Q( LDQ, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSBTRD reduces a real symmetric band matrix A to symmetric
!> tridiagonal form T by an orthogonal similarity transformation:
!> Q**T * A * Q = T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] VECT
!> \verbatim
!>          VECT is CHARACTER*1
!>          = 'N':  do not form Q;
!>          = 'V':  form Q;
!>          = 'U':  update a matrix X, by forming X*Q.
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
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is REAL array, dimension (LDAB,N)
!>          On entry, the upper or lower triangle of the symmetric band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>          On exit, the diagonal elements of AB are overwritten by the
!>          diagonal elements of the tridiagonal matrix T; if KD > 0, the
!>          elements on the first superdiagonal (if UPLO = 'U') or the
!>          first subdiagonal (if UPLO = 'L') are overwritten by the
!>          off-diagonal elements of T; the rest of AB is overwritten by
!>          values generated during the reduction.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The off-diagonal elements of the tridiagonal matrix T:
!>          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ,N)
!>          On entry, if VECT = 'U', then Q must contain an N-by-N
!>          matrix X; if VECT = 'N' or 'V', then Q need not be set.
!>
!>          On exit:
!>          if VECT = 'V', Q contains the N-by-N orthogonal matrix Q;
!>          if VECT = 'U', Q contains the product X*Q;
!>          if VECT = 'N', the array Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.
!>          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N)
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
!>  Modified by Linda Kaufman, Bell Labs.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SSBTRD(Vect,Uplo,N,Kd,Ab,Ldab,D,E,Q,Ldq,Work,Info)
      USE S_LSAME
      USE S_SLAR2V
      USE S_SLARGV
      USE S_SLARTG
      USE S_SLARTV
      USE S_SLASET
      USE S_SROT
      USE S_XERBLA
      IMPLICIT NONE
!*--SSBTRD174
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kd
      REAL , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL , INTENT(INOUT) , DIMENSION(*) :: D
      REAL , INTENT(OUT) , DIMENSION(*) :: E
      REAL , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , i2 , ibl , inca , incx , iqaend , iqb , iqend , j ,&
     &           j1 , j1end , j1inc , j2 , jend , jin , jinc , k , kd1 ,&
     &           kdm1 , kdn , l , last , lend , nq , nr , nrt
      LOGICAL :: initq , upper , wantq
      REAL :: temp
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
      initq = LSAME(Vect,'V')
      wantq = initq .OR. LSAME(Vect,'U')
      upper = LSAME(Uplo,'U')
      kd1 = Kd + 1
      kdm1 = Kd - 1
      incx = Ldab - 1
      iqend = 1
!
      Info = 0
      IF ( .NOT.wantq .AND. .NOT.LSAME(Vect,'N') ) THEN
         Info = -1
      ELSEIF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Kd<0 ) THEN
         Info = -4
      ELSEIF ( Ldab<kd1 ) THEN
         Info = -6
      ELSEIF ( Ldq<MAX(1,N) .AND. wantq ) THEN
         Info = -10
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSBTRD',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Initialize Q to the unit matrix, if needed
!
      IF ( initq ) CALL SLASET('Full',N,N,ZERO,ONE,Q,Ldq)
!
!     Wherever possible, plane rotations are generated and applied in
!     vector operations of length NR over the index set J1:J2:KD1.
!
!     The cosines and sines of the plane rotations are stored in the
!     arrays D and WORK.
!
      inca = kd1*Ldab
      kdn = MIN(N-1,Kd)
      IF ( upper ) THEN
!
         IF ( Kd>1 ) THEN
!
!           Reduce to tridiagonal form, working with upper triangle
!
            nr = 0
            j1 = kdn + 2
            j2 = 1
!
            DO i = 1 , N - 2
!
!              Reduce i-th row of matrix to tridiagonal form
!
               DO k = kdn + 1 , 2 , -1
                  j1 = j1 + kdn
                  j2 = j2 + kdn
!
                  IF ( nr>0 ) THEN
!
!                    generate plane rotations to annihilate nonzero
!                    elements which have been created outside the band
!
                     CALL SLARGV(nr,Ab(1,j1-1),inca,Work(j1),kd1,D(j1), &
     &                           kd1)
!
!                    apply rotations from the right
!
!
!                    Dependent on the the number of diagonals either
!                    SLARTV or SROT is used
!
                     IF ( nr>=2*Kd-1 ) THEN
                        DO l = 1 , Kd - 1
                           CALL SLARTV(nr,Ab(l+1,j1-1),inca,Ab(l,j1),   &
     &                                 inca,D(j1),Work(j1),kd1)
                        ENDDO
!
                     ELSE
                        jend = j1 + (nr-1)*kd1
                        DO jinc = j1 , jend , kd1
                           CALL SROT(kdm1,Ab(2,jinc-1),1,Ab(1,jinc),1,  &
     &                               D(jinc),Work(jinc))
                        ENDDO
                     ENDIF
                  ENDIF
!
!
                  IF ( k>2 ) THEN
                     IF ( k<=N-i+1 ) THEN
!
!                       generate plane rotation to annihilate a(i,i+k-1)
!                       within the band
!
                        CALL SLARTG(Ab(Kd-k+3,i+k-2),Ab(Kd-k+2,i+k-1),  &
     &                              D(i+k-1),Work(i+k-1),temp)
                        Ab(Kd-k+3,i+k-2) = temp
!
!                       apply rotation from the right
!
                        CALL SROT(k-3,Ab(Kd-k+4,i+k-2),1,               &
     &                            Ab(Kd-k+3,i+k-1),1,D(i+k-1),          &
     &                            Work(i+k-1))
                     ENDIF
                     nr = nr + 1
                     j1 = j1 - kdn - 1
                  ENDIF
!
!                 apply plane rotations from both sides to diagonal
!                 blocks
!
                  IF ( nr>0 ) CALL SLAR2V(nr,Ab(kd1,j1-1),Ab(kd1,j1),   &
     &                 Ab(Kd,j1),inca,D(j1),Work(j1),kd1)
!
!                 apply plane rotations from the left
!
                  IF ( nr>0 ) THEN
                     IF ( 2*Kd-1<nr ) THEN
!
!                    Dependent on the the number of diagonals either
!                    SLARTV or SROT is used
!
                        DO l = 1 , Kd - 1
                           IF ( j2+l>N ) THEN
                              nrt = nr - 1
                           ELSE
                              nrt = nr
                           ENDIF
                           IF ( nrt>0 )                                 &
     &                          CALL SLARTV(nrt,Ab(Kd-l,j1+l),inca,     &
     &                          Ab(Kd-l+1,j1+l),inca,D(j1),Work(j1),kd1)
                        ENDDO
                     ELSE
                        j1end = j1 + kd1*(nr-2)
                        IF ( j1end>=j1 ) THEN
                           DO jin = j1 , j1end , kd1
                              CALL SROT(Kd-1,Ab(Kd-1,jin+1),incx,       &
     &                                  Ab(Kd,jin+1),incx,D(jin),       &
     &                                  Work(jin))
                           ENDDO
                        ENDIF
                        lend = MIN(kdm1,N-j2)
                        last = j1end + kd1
                        IF ( lend>0 )                                   &
     &                       CALL SROT(lend,Ab(Kd-1,last+1),incx,       &
     &                       Ab(Kd,last+1),incx,D(last),Work(last))
                     ENDIF
                  ENDIF
!
                  IF ( wantq ) THEN
!
!                    accumulate product of plane rotations in Q
!
                     IF ( initq ) THEN
!
!                 take advantage of the fact that Q was
!                 initially the Identity matrix
!
                        iqend = MAX(iqend,j2)
                        i2 = MAX(0,k-3)
                        iqaend = 1 + i*Kd
                        IF ( k==2 ) iqaend = iqaend + Kd
                        iqaend = MIN(iqaend,iqend)
                        DO j = j1 , j2 , kd1
                           ibl = i - i2/kdm1
                           i2 = i2 + 1
                           iqb = MAX(1,j-ibl)
                           nq = 1 + iqaend - iqb
                           iqaend = MIN(iqaend+Kd,iqend)
                           CALL SROT(nq,Q(iqb,j-1),1,Q(iqb,j),1,D(j),   &
     &                               Work(j))
                        ENDDO
                     ELSE
!
                        DO j = j1 , j2 , kd1
                           CALL SROT(N,Q(1,j-1),1,Q(1,j),1,D(j),Work(j))
                        ENDDO
                     ENDIF
!
                  ENDIF
!
                  IF ( j2+kdn>N ) THEN
!
!                    adjust J2 to keep within the bounds of the matrix
!
                     nr = nr - 1
                     j2 = j2 - kdn - 1
                  ENDIF
!
                  DO j = j1 , j2 , kd1
!
!                    create nonzero element a(j-1,j+kd) outside the band
!                    and store it in WORK
!
                     Work(j+Kd) = Work(j)*Ab(1,j+Kd)
                     Ab(1,j+Kd) = D(j)*Ab(1,j+Kd)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
!
         IF ( Kd>0 ) THEN
!
!           copy off-diagonal elements to E
!
            DO i = 1 , N - 1
               E(i) = Ab(Kd,i+1)
            ENDDO
         ELSE
!
!           set E to zero if original matrix was diagonal
!
            DO i = 1 , N - 1
               E(i) = ZERO
            ENDDO
         ENDIF
!
!        copy diagonal elements to D
!
         DO i = 1 , N
            D(i) = Ab(kd1,i)
         ENDDO
!
      ELSE
!
         IF ( Kd>1 ) THEN
!
!           Reduce to tridiagonal form, working with lower triangle
!
            nr = 0
            j1 = kdn + 2
            j2 = 1
!
            DO i = 1 , N - 2
!
!              Reduce i-th column of matrix to tridiagonal form
!
               DO k = kdn + 1 , 2 , -1
                  j1 = j1 + kdn
                  j2 = j2 + kdn
!
                  IF ( nr>0 ) THEN
!
!                    generate plane rotations to annihilate nonzero
!                    elements which have been created outside the band
!
                     CALL SLARGV(nr,Ab(kd1,j1-kd1),inca,Work(j1),kd1,   &
     &                           D(j1),kd1)
!
!                    apply plane rotations from one side
!
!
!                    Dependent on the the number of diagonals either
!                    SLARTV or SROT is used
!
                     IF ( nr>2*Kd-1 ) THEN
                        DO l = 1 , Kd - 1
                           CALL SLARTV(nr,Ab(kd1-l,j1-kd1+l),inca,      &
     &                                 Ab(kd1-l+1,j1-kd1+l),inca,D(j1), &
     &                                 Work(j1),kd1)
                        ENDDO
                     ELSE
                        jend = j1 + kd1*(nr-1)
                        DO jinc = j1 , jend , kd1
                           CALL SROT(kdm1,Ab(Kd,jinc-Kd),incx,          &
     &                               Ab(kd1,jinc-Kd),incx,D(jinc),      &
     &                               Work(jinc))
                        ENDDO
                     ENDIF
!
                  ENDIF
!
                  IF ( k>2 ) THEN
                     IF ( k<=N-i+1 ) THEN
!
!                       generate plane rotation to annihilate a(i+k-1,i)
!                       within the band
!
                        CALL SLARTG(Ab(k-1,i),Ab(k,i),D(i+k-1),         &
     &                              Work(i+k-1),temp)
                        Ab(k-1,i) = temp
!
!                       apply rotation from the left
!
                        CALL SROT(k-3,Ab(k-2,i+1),Ldab-1,Ab(k-1,i+1),   &
     &                            Ldab-1,D(i+k-1),Work(i+k-1))
                     ENDIF
                     nr = nr + 1
                     j1 = j1 - kdn - 1
                  ENDIF
!
!                 apply plane rotations from both sides to diagonal
!                 blocks
!
                  IF ( nr>0 ) CALL SLAR2V(nr,Ab(1,j1-1),Ab(1,j1),       &
     &                 Ab(2,j1-1),inca,D(j1),Work(j1),kd1)
!
!                 apply plane rotations from the right
!
!
!                    Dependent on the the number of diagonals either
!                    SLARTV or SROT is used
!
                  IF ( nr>0 ) THEN
                     IF ( nr>2*Kd-1 ) THEN
                        DO l = 1 , Kd - 1
                           IF ( j2+l>N ) THEN
                              nrt = nr - 1
                           ELSE
                              nrt = nr
                           ENDIF
                           IF ( nrt>0 )                                 &
     &                          CALL SLARTV(nrt,Ab(l+2,j1-1),inca,      &
     &                          Ab(l+1,j1),inca,D(j1),Work(j1),kd1)
                        ENDDO
                     ELSE
                        j1end = j1 + kd1*(nr-2)
                        IF ( j1end>=j1 ) THEN
                           DO j1inc = j1 , j1end , kd1
                              CALL SROT(kdm1,Ab(3,j1inc-1),1,Ab(2,j1inc)&
     &                                  ,1,D(j1inc),Work(j1inc))
                           ENDDO
                        ENDIF
                        lend = MIN(kdm1,N-j2)
                        last = j1end + kd1
                        IF ( lend>0 )                                   &
     &                       CALL SROT(lend,Ab(3,last-1),1,Ab(2,last),1,&
     &                       D(last),Work(last))
                     ENDIF
                  ENDIF
!
!
!
                  IF ( wantq ) THEN
!
!                    accumulate product of plane rotations in Q
!
                     IF ( initq ) THEN
!
!                 take advantage of the fact that Q was
!                 initially the Identity matrix
!
                        iqend = MAX(iqend,j2)
                        i2 = MAX(0,k-3)
                        iqaend = 1 + i*Kd
                        IF ( k==2 ) iqaend = iqaend + Kd
                        iqaend = MIN(iqaend,iqend)
                        DO j = j1 , j2 , kd1
                           ibl = i - i2/kdm1
                           i2 = i2 + 1
                           iqb = MAX(1,j-ibl)
                           nq = 1 + iqaend - iqb
                           iqaend = MIN(iqaend+Kd,iqend)
                           CALL SROT(nq,Q(iqb,j-1),1,Q(iqb,j),1,D(j),   &
     &                               Work(j))
                        ENDDO
                     ELSE
!
                        DO j = j1 , j2 , kd1
                           CALL SROT(N,Q(1,j-1),1,Q(1,j),1,D(j),Work(j))
                        ENDDO
                     ENDIF
                  ENDIF
!
                  IF ( j2+kdn>N ) THEN
!
!                    adjust J2 to keep within the bounds of the matrix
!
                     nr = nr - 1
                     j2 = j2 - kdn - 1
                  ENDIF
!
                  DO j = j1 , j2 , kd1
!
!                    create nonzero element a(j+kd,j-1) outside the
!                    band and store it in WORK
!
                     Work(j+Kd) = Work(j)*Ab(kd1,j)
                     Ab(kd1,j) = D(j)*Ab(kd1,j)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
!
         IF ( Kd>0 ) THEN
!
!           copy off-diagonal elements to E
!
            DO i = 1 , N - 1
               E(i) = Ab(2,i)
            ENDDO
         ELSE
!
!           set E to zero if original matrix was diagonal
!
            DO i = 1 , N - 1
               E(i) = ZERO
            ENDDO
         ENDIF
!
!        copy diagonal elements to D
!
         DO i = 1 , N
            D(i) = Ab(1,i)
         ENDDO
      ENDIF
!
!
!     End of SSBTRD
!
      END SUBROUTINE SSBTRD
