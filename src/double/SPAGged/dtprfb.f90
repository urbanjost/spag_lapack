!*==dtprfb.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DTPRFB applies a real or complex "triangular-pentagonal" blocked reflector to a real or complex matrix, which is composed of two blocks.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DTPRFB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtprfb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtprfb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtprfb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTPRFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L,
!                          V, LDV, T, LDT, A, LDA, B, LDB, WORK, LDWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER DIRECT, SIDE, STOREV, TRANS
!       INTEGER   K, L, LDA, LDB, LDT, LDV, LDWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), T( LDT, * ),
!      $          V( LDV, * ), WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTPRFB applies a real "triangular-pentagonal" block reflector H or its
!> transpose H**T to a real matrix C, which is composed of two
!> blocks A and B, either from the left or right.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply H or H**T from the Left
!>          = 'R': apply H or H**T from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply H (No transpose)
!>          = 'T': apply H**T (Transpose)
!> \endverbatim
!>
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Indicates how H is formed from a product of elementary
!>          reflectors
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Indicates how the vectors which define the elementary
!>          reflectors are stored:
!>          = 'C': Columns
!>          = 'R': Rows
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix B.
!>          M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix B.
!>          N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The order of the matrix T, i.e. the number of elementary
!>          reflectors whose product defines the block reflector.
!>          K >= 0.
!> \endverbatim
!>
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>          The order of the trapezoidal part of V.
!>          K >= L >= 0.  See Further Details.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension
!>                                (LDV,K) if STOREV = 'C'
!>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!>          The pentagonal matrix V, which contains the elementary reflectors
!>          H(1), H(2), ..., H(K).  See Further Details.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!>          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!>          if STOREV = 'R', LDV >= K.
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,K)
!>          The triangular K-by-K matrix T in the representation of the
!>          block reflector.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.
!>          LDT >= K.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension
!>          (LDA,N) if SIDE = 'L' or (LDA,K) if SIDE = 'R'
!>          On entry, the K-by-N or M-by-K matrix A.
!>          On exit, A is overwritten by the corresponding block of
!>          H*C or H**T*C or C*H or C*H**T.  See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDA >= max(1,K);
!>          If SIDE = 'R', LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          On entry, the M-by-N matrix B.
!>          On exit, B is overwritten by the corresponding block of
!>          H*C or H**T*C or C*H or C*H**T.  See Further Details.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.
!>          LDB >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension
!>          (LDWORK,N) if SIDE = 'L',
!>          (LDWORK,K) if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of the array WORK.
!>          If SIDE = 'L', LDWORK >= K;
!>          if SIDE = 'R', LDWORK >= M.
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
!> \ingroup doubleOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix C is a composite matrix formed from blocks A and B.
!>  The block B is of size M-by-N; if SIDE = 'R', A is of size M-by-K,
!>  and if SIDE = 'L', A is of size K-by-N.
!>
!>  If SIDE = 'R' and DIRECT = 'F', C = [A B].
!>
!>  If SIDE = 'L' and DIRECT = 'F', C = [A]
!>                                      [B].
!>
!>  If SIDE = 'R' and DIRECT = 'B', C = [B A].
!>
!>  If SIDE = 'L' and DIRECT = 'B', C = [B]
!>                                      [A].
!>
!>  The pentagonal matrix V is composed of a rectangular block V1 and a
!>  trapezoidal block V2.  The size of the trapezoidal block is determined by
!>  the parameter L, where 0<=L<=K.  If L=K, the V2 block of V is triangular;
!>  if L=0, there is no trapezoidal block, thus V = V1 is rectangular.
!>
!>  If DIRECT = 'F' and STOREV = 'C':  V = [V1]
!>                                         [V2]
!>     - V2 is upper trapezoidal (first L rows of K-by-K upper triangular)
!>
!>  If DIRECT = 'F' and STOREV = 'R':  V = [V1 V2]
!>
!>     - V2 is lower trapezoidal (first L columns of K-by-K lower triangular)
!>
!>  If DIRECT = 'B' and STOREV = 'C':  V = [V2]
!>                                         [V1]
!>     - V2 is lower trapezoidal (last L rows of K-by-K lower triangular)
!>
!>  If DIRECT = 'B' and STOREV = 'R':  V = [V2 V1]
!>
!>     - V2 is upper trapezoidal (last L columns of K-by-K upper triangular)
!>
!>  If STOREV = 'C' and SIDE = 'L', V is M-by-K with V2 L-by-K.
!>
!>  If STOREV = 'C' and SIDE = 'R', V is N-by-K with V2 L-by-K.
!>
!>  If STOREV = 'R' and SIDE = 'L', V is K-by-M with V2 K-by-L.
!>
!>  If STOREV = 'R' and SIDE = 'R', V is K-by-N with V2 K-by-L.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DTPRFB(Side,Trans,Direct,Storev,M,N,K,L,V,Ldv,T,Ldt,A, &
     &                  Lda,B,Ldb,Work,Ldwork)
      USE F77KINDS                        
      USE S_DGEMM
      USE S_DTRMM
      USE S_LSAME
      IMPLICIT NONE
!*--DTPRFB259
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0 , ZERO = 0.0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Trans
      CHARACTER :: Direct
      CHARACTER :: Storev
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: L
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: backward , column , forward , left , right , row
      INTEGER :: i , j , kp , mp , np
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  ==========================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( M<=0 .OR. N<=0 .OR. K<=0 .OR. L<0 ) RETURN
!
      IF ( LSAME(Storev,'C') ) THEN
         column = .TRUE.
         row = .FALSE.
      ELSEIF ( LSAME(Storev,'R') ) THEN
         column = .FALSE.
         row = .TRUE.
      ELSE
         column = .FALSE.
         row = .FALSE.
      ENDIF
!
      IF ( LSAME(Side,'L') ) THEN
         left = .TRUE.
         right = .FALSE.
      ELSEIF ( LSAME(Side,'R') ) THEN
         left = .FALSE.
         right = .TRUE.
      ELSE
         left = .FALSE.
         right = .FALSE.
      ENDIF
!
      IF ( LSAME(Direct,'F') ) THEN
         forward = .TRUE.
         backward = .FALSE.
      ELSEIF ( LSAME(Direct,'B') ) THEN
         forward = .FALSE.
         backward = .TRUE.
      ELSE
         forward = .FALSE.
         backward = .FALSE.
      ENDIF
!
! ---------------------------------------------------------------------------
!
      IF ( column .AND. forward .AND. left ) THEN
!
! ---------------------------------------------------------------------------
!
!        Let  W =  [ I ]    (K-by-K)
!                  [ V ]    (M-by-K)
!
!        Form  H C  or  H**T C  where  C = [ A ]  (K-by-N)
!                                          [ B ]  (M-by-N)
!
!        H = I - W T W**T          or  H**T = I - W T**T W**T
!
!        A = A -   T (A + V**T B)  or  A = A -   T**T (A + V**T B)
!        B = B - V T (A + V**T B)  or  B = B - V T**T (A + V**T B)
!
! ---------------------------------------------------------------------------
!
         mp = MIN(M-L+1,M)
         kp = MIN(L+1,K)
!
         DO j = 1 , N
            DO i = 1 , L
               Work(i,j) = B(M-L+i,j)
            ENDDO
         ENDDO
         CALL DTRMM('L','U','T','N',L,N,ONE,V(mp,1),Ldv,Work,Ldwork)
         CALL DGEMM('T','N',L,N,M-L,ONE,V,Ldv,B,Ldb,ONE,Work,Ldwork)
         CALL DGEMM('T','N',K-L,N,M,ONE,V(1,kp),Ldv,B,Ldb,ZERO,         &
     &              Work(kp,1),Ldwork)
!
         DO j = 1 , N
            DO i = 1 , K
               Work(i,j) = Work(i,j) + A(i,j)
            ENDDO
         ENDDO
!
         CALL DTRMM('L','U',Trans,'N',K,N,ONE,T,Ldt,Work,Ldwork)
!
         DO j = 1 , N
            DO i = 1 , K
               A(i,j) = A(i,j) - Work(i,j)
            ENDDO
         ENDDO
!
         CALL DGEMM('N','N',M-L,N,K,-ONE,V,Ldv,Work,Ldwork,ONE,B,Ldb)
         CALL DGEMM('N','N',L,N,K-L,-ONE,V(mp,kp),Ldv,Work(kp,1),Ldwork,&
     &              ONE,B(mp,1),Ldb)
         CALL DTRMM('L','U','N','N',L,N,ONE,V(mp,1),Ldv,Work,Ldwork)
         DO j = 1 , N
            DO i = 1 , L
               B(M-L+i,j) = B(M-L+i,j) - Work(i,j)
            ENDDO
         ENDDO
!
! ---------------------------------------------------------------------------
!
      ELSEIF ( column .AND. forward .AND. right ) THEN
!
! ---------------------------------------------------------------------------
!
!        Let  W =  [ I ]    (K-by-K)
!                  [ V ]    (N-by-K)
!
!        Form  C H or  C H**T  where  C = [ A B ] (A is M-by-K, B is M-by-N)
!
!        H = I - W T W**T          or  H**T = I - W T**T W**T
!
!        A = A - (A + B V) T      or  A = A - (A + B V) T**T
!        B = B - (A + B V) T V**T  or  B = B - (A + B V) T**T V**T
!
! ---------------------------------------------------------------------------
!
         np = MIN(N-L+1,N)
         kp = MIN(L+1,K)
!
         DO j = 1 , L
            DO i = 1 , M
               Work(i,j) = B(i,N-L+j)
            ENDDO
         ENDDO
         CALL DTRMM('R','U','N','N',M,L,ONE,V(np,1),Ldv,Work,Ldwork)
         CALL DGEMM('N','N',M,L,N-L,ONE,B,Ldb,V,Ldv,ONE,Work,Ldwork)
         CALL DGEMM('N','N',M,K-L,N,ONE,B,Ldb,V(1,kp),Ldv,ZERO,         &
     &              Work(1,kp),Ldwork)
!
         DO j = 1 , K
            DO i = 1 , M
               Work(i,j) = Work(i,j) + A(i,j)
            ENDDO
         ENDDO
!
         CALL DTRMM('R','U',Trans,'N',M,K,ONE,T,Ldt,Work,Ldwork)
!
         DO j = 1 , K
            DO i = 1 , M
               A(i,j) = A(i,j) - Work(i,j)
            ENDDO
         ENDDO
!
         CALL DGEMM('N','T',M,N-L,K,-ONE,Work,Ldwork,V,Ldv,ONE,B,Ldb)
         CALL DGEMM('N','T',M,L,K-L,-ONE,Work(1,kp),Ldwork,V(np,kp),Ldv,&
     &              ONE,B(1,np),Ldb)
         CALL DTRMM('R','U','T','N',M,L,ONE,V(np,1),Ldv,Work,Ldwork)
         DO j = 1 , L
            DO i = 1 , M
               B(i,N-L+j) = B(i,N-L+j) - Work(i,j)
            ENDDO
         ENDDO
!
! ---------------------------------------------------------------------------
!
      ELSEIF ( column .AND. backward .AND. left ) THEN
!
! ---------------------------------------------------------------------------
!
!        Let  W =  [ V ]    (M-by-K)
!                  [ I ]    (K-by-K)
!
!        Form  H C  or  H**T C  where  C = [ B ]  (M-by-N)
!                                          [ A ]  (K-by-N)
!
!        H = I - W T W**T          or  H**T = I - W T**T W**T
!
!        A = A -   T (A + V**T B)  or  A = A -   T**T (A + V**T B)
!        B = B - V T (A + V**T B)  or  B = B - V T**T (A + V**T B)
!
! ---------------------------------------------------------------------------
!
         mp = MIN(L+1,M)
         kp = MIN(K-L+1,K)
!
         DO j = 1 , N
            DO i = 1 , L
               Work(K-L+i,j) = B(i,j)
            ENDDO
         ENDDO
!
         CALL DTRMM('L','L','T','N',L,N,ONE,V(1,kp),Ldv,Work(kp,1),     &
     &              Ldwork)
         CALL DGEMM('T','N',L,N,M-L,ONE,V(mp,kp),Ldv,B(mp,1),Ldb,ONE,   &
     &              Work(kp,1),Ldwork)
         CALL DGEMM('T','N',K-L,N,M,ONE,V,Ldv,B,Ldb,ZERO,Work,Ldwork)
!
         DO j = 1 , N
            DO i = 1 , K
               Work(i,j) = Work(i,j) + A(i,j)
            ENDDO
         ENDDO
!
         CALL DTRMM('L','L',Trans,'N',K,N,ONE,T,Ldt,Work,Ldwork)
!
         DO j = 1 , N
            DO i = 1 , K
               A(i,j) = A(i,j) - Work(i,j)
            ENDDO
         ENDDO
!
         CALL DGEMM('N','N',M-L,N,K,-ONE,V(mp,1),Ldv,Work,Ldwork,ONE,   &
     &              B(mp,1),Ldb)
         CALL DGEMM('N','N',L,N,K-L,-ONE,V,Ldv,Work,Ldwork,ONE,B,Ldb)
         CALL DTRMM('L','L','N','N',L,N,ONE,V(1,kp),Ldv,Work(kp,1),     &
     &              Ldwork)
         DO j = 1 , N
            DO i = 1 , L
               B(i,j) = B(i,j) - Work(K-L+i,j)
            ENDDO
         ENDDO
!
! ---------------------------------------------------------------------------
!
      ELSEIF ( column .AND. backward .AND. right ) THEN
!
! ---------------------------------------------------------------------------
!
!        Let  W =  [ V ]    (N-by-K)
!                  [ I ]    (K-by-K)
!
!        Form  C H  or  C H**T  where  C = [ B A ] (B is M-by-N, A is M-by-K)
!
!        H = I - W T W**T          or  H**T = I - W T**T W**T
!
!        A = A - (A + B V) T      or  A = A - (A + B V) T**T
!        B = B - (A + B V) T V**T  or  B = B - (A + B V) T**T V**T
!
! ---------------------------------------------------------------------------
!
         np = MIN(L+1,N)
         kp = MIN(K-L+1,K)
!
         DO j = 1 , L
            DO i = 1 , M
               Work(i,K-L+j) = B(i,j)
            ENDDO
         ENDDO
         CALL DTRMM('R','L','N','N',M,L,ONE,V(1,kp),Ldv,Work(1,kp),     &
     &              Ldwork)
         CALL DGEMM('N','N',M,L,N-L,ONE,B(1,np),Ldb,V(np,kp),Ldv,ONE,   &
     &              Work(1,kp),Ldwork)
         CALL DGEMM('N','N',M,K-L,N,ONE,B,Ldb,V,Ldv,ZERO,Work,Ldwork)
!
         DO j = 1 , K
            DO i = 1 , M
               Work(i,j) = Work(i,j) + A(i,j)
            ENDDO
         ENDDO
!
         CALL DTRMM('R','L',Trans,'N',M,K,ONE,T,Ldt,Work,Ldwork)
!
         DO j = 1 , K
            DO i = 1 , M
               A(i,j) = A(i,j) - Work(i,j)
            ENDDO
         ENDDO
!
         CALL DGEMM('N','T',M,N-L,K,-ONE,Work,Ldwork,V(np,1),Ldv,ONE,   &
     &              B(1,np),Ldb)
         CALL DGEMM('N','T',M,L,K-L,-ONE,Work,Ldwork,V,Ldv,ONE,B,Ldb)
         CALL DTRMM('R','L','T','N',M,L,ONE,V(1,kp),Ldv,Work(1,kp),     &
     &              Ldwork)
         DO j = 1 , L
            DO i = 1 , M
               B(i,j) = B(i,j) - Work(i,K-L+j)
            ENDDO
         ENDDO
!
! ---------------------------------------------------------------------------
!
      ELSEIF ( row .AND. forward .AND. left ) THEN
!
! ---------------------------------------------------------------------------
!
!        Let  W =  [ I V ] ( I is K-by-K, V is K-by-M )
!
!        Form  H C  or  H**T C  where  C = [ A ]  (K-by-N)
!                                          [ B ]  (M-by-N)
!
!        H = I - W**T T W          or  H**T = I - W**T T**T W
!
!        A = A -     T (A + V B)  or  A = A -     T**T (A + V B)
!        B = B - V**T T (A + V B)  or  B = B - V**T T**T (A + V B)
!
! ---------------------------------------------------------------------------
!
         mp = MIN(M-L+1,M)
         kp = MIN(L+1,K)
!
         DO j = 1 , N
            DO i = 1 , L
               Work(i,j) = B(M-L+i,j)
            ENDDO
         ENDDO
         CALL DTRMM('L','L','N','N',L,N,ONE,V(1,mp),Ldv,Work,Ldb)
         CALL DGEMM('N','N',L,N,M-L,ONE,V,Ldv,B,Ldb,ONE,Work,Ldwork)
         CALL DGEMM('N','N',K-L,N,M,ONE,V(kp,1),Ldv,B,Ldb,ZERO,         &
     &              Work(kp,1),Ldwork)
!
         DO j = 1 , N
            DO i = 1 , K
               Work(i,j) = Work(i,j) + A(i,j)
            ENDDO
         ENDDO
!
         CALL DTRMM('L','U',Trans,'N',K,N,ONE,T,Ldt,Work,Ldwork)
!
         DO j = 1 , N
            DO i = 1 , K
               A(i,j) = A(i,j) - Work(i,j)
            ENDDO
         ENDDO
!
         CALL DGEMM('T','N',M-L,N,K,-ONE,V,Ldv,Work,Ldwork,ONE,B,Ldb)
         CALL DGEMM('T','N',L,N,K-L,-ONE,V(kp,mp),Ldv,Work(kp,1),Ldwork,&
     &              ONE,B(mp,1),Ldb)
         CALL DTRMM('L','L','T','N',L,N,ONE,V(1,mp),Ldv,Work,Ldwork)
         DO j = 1 , N
            DO i = 1 , L
               B(M-L+i,j) = B(M-L+i,j) - Work(i,j)
            ENDDO
         ENDDO
!
! ---------------------------------------------------------------------------
!
      ELSEIF ( row .AND. forward .AND. right ) THEN
!
! ---------------------------------------------------------------------------
!
!        Let  W =  [ I V ] ( I is K-by-K, V is K-by-N )
!
!        Form  C H  or  C H**T  where  C = [ A B ] (A is M-by-K, B is M-by-N)
!
!        H = I - W**T T W            or  H**T = I - W**T T**T W
!
!        A = A - (A + B V**T) T      or  A = A - (A + B V**T) T**T
!        B = B - (A + B V**T) T V    or  B = B - (A + B V**T) T**T V
!
! ---------------------------------------------------------------------------
!
         np = MIN(N-L+1,N)
         kp = MIN(L+1,K)
!
         DO j = 1 , L
            DO i = 1 , M
               Work(i,j) = B(i,N-L+j)
            ENDDO
         ENDDO
         CALL DTRMM('R','L','T','N',M,L,ONE,V(1,np),Ldv,Work,Ldwork)
         CALL DGEMM('N','T',M,L,N-L,ONE,B,Ldb,V,Ldv,ONE,Work,Ldwork)
         CALL DGEMM('N','T',M,K-L,N,ONE,B,Ldb,V(kp,1),Ldv,ZERO,         &
     &              Work(1,kp),Ldwork)
!
         DO j = 1 , K
            DO i = 1 , M
               Work(i,j) = Work(i,j) + A(i,j)
            ENDDO
         ENDDO
!
         CALL DTRMM('R','U',Trans,'N',M,K,ONE,T,Ldt,Work,Ldwork)
!
         DO j = 1 , K
            DO i = 1 , M
               A(i,j) = A(i,j) - Work(i,j)
            ENDDO
         ENDDO
!
         CALL DGEMM('N','N',M,N-L,K,-ONE,Work,Ldwork,V,Ldv,ONE,B,Ldb)
         CALL DGEMM('N','N',M,L,K-L,-ONE,Work(1,kp),Ldwork,V(kp,np),Ldv,&
     &              ONE,B(1,np),Ldb)
         CALL DTRMM('R','L','N','N',M,L,ONE,V(1,np),Ldv,Work,Ldwork)
         DO j = 1 , L
            DO i = 1 , M
               B(i,N-L+j) = B(i,N-L+j) - Work(i,j)
            ENDDO
         ENDDO
!
! ---------------------------------------------------------------------------
!
      ELSEIF ( row .AND. backward .AND. left ) THEN
!
! ---------------------------------------------------------------------------
!
!        Let  W =  [ V I ] ( I is K-by-K, V is K-by-M )
!
!        Form  H C  or  H**T C  where  C = [ B ]  (M-by-N)
!                                          [ A ]  (K-by-N)
!
!        H = I - W**T T W          or  H**T = I - W**T T**T W
!
!        A = A -     T (A + V B)  or  A = A -     T**T (A + V B)
!        B = B - V**T T (A + V B)  or  B = B - V**T T**T (A + V B)
!
! ---------------------------------------------------------------------------
!
         mp = MIN(L+1,M)
         kp = MIN(K-L+1,K)
!
         DO j = 1 , N
            DO i = 1 , L
               Work(K-L+i,j) = B(i,j)
            ENDDO
         ENDDO
         CALL DTRMM('L','U','N','N',L,N,ONE,V(kp,1),Ldv,Work(kp,1),     &
     &              Ldwork)
         CALL DGEMM('N','N',L,N,M-L,ONE,V(kp,mp),Ldv,B(mp,1),Ldb,ONE,   &
     &              Work(kp,1),Ldwork)
         CALL DGEMM('N','N',K-L,N,M,ONE,V,Ldv,B,Ldb,ZERO,Work,Ldwork)
!
         DO j = 1 , N
            DO i = 1 , K
               Work(i,j) = Work(i,j) + A(i,j)
            ENDDO
         ENDDO
!
         CALL DTRMM('L','L ',Trans,'N',K,N,ONE,T,Ldt,Work,Ldwork)
!
         DO j = 1 , N
            DO i = 1 , K
               A(i,j) = A(i,j) - Work(i,j)
            ENDDO
         ENDDO
!
         CALL DGEMM('T','N',M-L,N,K,-ONE,V(1,mp),Ldv,Work,Ldwork,ONE,   &
     &              B(mp,1),Ldb)
         CALL DGEMM('T','N',L,N,K-L,-ONE,V,Ldv,Work,Ldwork,ONE,B,Ldb)
         CALL DTRMM('L','U','T','N',L,N,ONE,V(kp,1),Ldv,Work(kp,1),     &
     &              Ldwork)
         DO j = 1 , N
            DO i = 1 , L
               B(i,j) = B(i,j) - Work(K-L+i,j)
            ENDDO
         ENDDO
!
! ---------------------------------------------------------------------------
!
      ELSEIF ( row .AND. backward .AND. right ) THEN
!
! ---------------------------------------------------------------------------
!
!        Let  W =  [ V I ] ( I is K-by-K, V is K-by-N )
!
!        Form  C H  or  C H**T  where  C = [ B A ] (A is M-by-K, B is M-by-N)
!
!        H = I - W**T T W            or  H**T = I - W**T T**T W
!
!        A = A - (A + B V**T) T      or  A = A - (A + B V**T) T**T
!        B = B - (A + B V**T) T V    or  B = B - (A + B V**T) T**T V
!
! ---------------------------------------------------------------------------
!
         np = MIN(L+1,N)
         kp = MIN(K-L+1,K)
!
         DO j = 1 , L
            DO i = 1 , M
               Work(i,K-L+j) = B(i,j)
            ENDDO
         ENDDO
         CALL DTRMM('R','U','T','N',M,L,ONE,V(kp,1),Ldv,Work(1,kp),     &
     &              Ldwork)
         CALL DGEMM('N','T',M,L,N-L,ONE,B(1,np),Ldb,V(kp,np),Ldv,ONE,   &
     &              Work(1,kp),Ldwork)
         CALL DGEMM('N','T',M,K-L,N,ONE,B,Ldb,V,Ldv,ZERO,Work,Ldwork)
!
         DO j = 1 , K
            DO i = 1 , M
               Work(i,j) = Work(i,j) + A(i,j)
            ENDDO
         ENDDO
!
         CALL DTRMM('R','L',Trans,'N',M,K,ONE,T,Ldt,Work,Ldwork)
!
         DO j = 1 , K
            DO i = 1 , M
               A(i,j) = A(i,j) - Work(i,j)
            ENDDO
         ENDDO
!
         CALL DGEMM('N','N',M,N-L,K,-ONE,Work,Ldwork,V(1,np),Ldv,ONE,   &
     &              B(1,np),Ldb)
         CALL DGEMM('N','N',M,L,K-L,-ONE,Work,Ldwork,V,Ldv,ONE,B,Ldb)
         CALL DTRMM('R','U','N','N',M,L,ONE,V(kp,1),Ldv,Work(1,kp),     &
     &              Ldwork)
         DO j = 1 , L
            DO i = 1 , M
               B(i,j) = B(i,j) - Work(i,K-L+j)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!     End of DTPRFB
!
      END SUBROUTINE DTPRFB
