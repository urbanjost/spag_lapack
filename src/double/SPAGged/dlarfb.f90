!*==dlarfb.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARFB applies a block reflector or its transpose to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARFB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
!                          T, LDT, C, LDC, WORK, LDWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, SIDE, STOREV, TRANS
!       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
!      $                   WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARFB applies a real block reflector H or its transpose H**T to a
!> real m by n matrix C, from either the left or the right.
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
!>          = 'C': Columnwise
!>          = 'R': Rowwise
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The order of the matrix T (= the number of elementary
!>          reflectors whose product defines the block reflector).
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension
!>                                (LDV,K) if STOREV = 'C'
!>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!>          The matrix V. See Further Details.
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
!>          The triangular k by k matrix T in the representation of the
!>          block reflector.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= K.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LDWORK,K)
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of the array WORK.
!>          If SIDE = 'L', LDWORK >= max(1,N);
!>          if SIDE = 'R', LDWORK >= max(1,M).
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
!> \date June 2013
!
!> \ingroup doubleOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The shape of the matrix V and the storage of the vectors which define
!>  the H(i) is best illustrated by the following example with n = 5 and
!>  k = 3. The elements equal to 1 are not stored; the corresponding
!>  array elements are modified but restored on exit. The rest of the
!>  array is not used.
!>
!>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!>
!>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!>                   ( v1  1    )                     (     1 v2 v2 v2 )
!>                   ( v1 v2  1 )                     (        1 v3 v3 )
!>                   ( v1 v2 v3 )
!>                   ( v1 v2 v3 )
!>
!>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!>
!>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!>                   (     1 v3 )
!>                   (        1 )
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLARFB(Side,Trans,Direct,Storev,M,N,K,V,Ldv,T,Ldt,C,   &
     &                  Ldc,Work,Ldwork)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DGEMM
      USE S_DTRMM
      USE S_LSAME
      IMPLICIT NONE
!*--DLARFB206
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
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
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , j
      CHARACTER :: transt
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
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( M<=0 .OR. N<=0 ) RETURN
!
      IF ( LSAME(Trans,'N') ) THEN
         transt = 'T'
      ELSE
         transt = 'N'
      ENDIF
!
      IF ( LSAME(Storev,'C') ) THEN
!
         IF ( LSAME(Direct,'F') ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
            IF ( LSAME(Side,'L') ) THEN
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
!
!              W := C1**T
!
               DO j = 1 , K
                  CALL DCOPY(N,C(j,1),Ldc,Work(1,j),1)
               ENDDO
!
!              W := W * V1
!
               CALL DTRMM('Right','Lower','No transpose','Unit',N,K,ONE,&
     &                    V,Ldv,Work,Ldwork)
!
!                 W := W + C2**T * V2
!
               IF ( M>K ) CALL DGEMM('Transpose','No transpose',N,K,M-K,&
     &                               ONE,C(K+1,1),Ldc,V(K+1,1),Ldv,ONE, &
     &                               Work,Ldwork)
!
!              W := W * T**T  or  W * T
!
               CALL DTRMM('Right','Upper',transt,'Non-unit',N,K,ONE,T,  &
     &                    Ldt,Work,Ldwork)
!
!              C := C - V * W**T
!
!
!                 C2 := C2 - V2 * W**T
!
               IF ( M>K ) CALL DGEMM('No transpose','Transpose',M-K,N,K,&
     &                               -ONE,V(K+1,1),Ldv,Work,Ldwork,ONE, &
     &                               C(K+1,1),Ldc)
!
!              W := W * V1**T
!
               CALL DTRMM('Right','Lower','Transpose','Unit',N,K,ONE,V, &
     &                    Ldv,Work,Ldwork)
!
!              C1 := C1 - W**T
!
               DO j = 1 , K
                  DO i = 1 , N
                     C(j,i) = C(j,i) - Work(i,j)
                  ENDDO
               ENDDO
!
            ELSEIF ( LSAME(Side,'R') ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO j = 1 , K
                  CALL DCOPY(M,C(1,j),1,Work(1,j),1)
               ENDDO
!
!              W := W * V1
!
               CALL DTRMM('Right','Lower','No transpose','Unit',M,K,ONE,&
     &                    V,Ldv,Work,Ldwork)
!
!                 W := W + C2 * V2
!
               IF ( N>K ) CALL DGEMM('No transpose','No transpose',M,K, &
     &                               N-K,ONE,C(1,K+1),Ldc,V(K+1,1),Ldv, &
     &                               ONE,Work,Ldwork)
!
!              W := W * T  or  W * T**T
!
               CALL DTRMM('Right','Upper',Trans,'Non-unit',M,K,ONE,T,   &
     &                    Ldt,Work,Ldwork)
!
!              C := C - W * V**T
!
!
!                 C2 := C2 - W * V2**T
!
               IF ( N>K ) CALL DGEMM('No transpose','Transpose',M,N-K,K,&
     &                               -ONE,Work,Ldwork,V(K+1,1),Ldv,ONE, &
     &                               C(1,K+1),Ldc)
!
!              W := W * V1**T
!
               CALL DTRMM('Right','Lower','Transpose','Unit',M,K,ONE,V, &
     &                    Ldv,Work,Ldwork)
!
!              C1 := C1 - W
!
               DO j = 1 , K
                  DO i = 1 , M
                     C(i,j) = C(i,j) - Work(i,j)
                  ENDDO
               ENDDO
            ENDIF
!
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
         ELSEIF ( LSAME(Side,'L') ) THEN
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
!
!              W := C2**T
!
            DO j = 1 , K
               CALL DCOPY(N,C(M-K+j,1),Ldc,Work(1,j),1)
            ENDDO
!
!              W := W * V2
!
            CALL DTRMM('Right','Upper','No transpose','Unit',N,K,ONE,   &
     &                 V(M-K+1,1),Ldv,Work,Ldwork)
!
!                 W := W + C1**T * V1
!
            IF ( M>K ) CALL DGEMM('Transpose','No transpose',N,K,M-K,   &
     &                            ONE,C,Ldc,V,Ldv,ONE,Work,Ldwork)
!
!              W := W * T**T  or  W * T
!
            CALL DTRMM('Right','Lower',transt,'Non-unit',N,K,ONE,T,Ldt, &
     &                 Work,Ldwork)
!
!              C := C - V * W**T
!
!
!                 C1 := C1 - V1 * W**T
!
            IF ( M>K ) CALL DGEMM('No transpose','Transpose',M-K,N,K,   &
     &                            -ONE,V,Ldv,Work,Ldwork,ONE,C,Ldc)
!
!              W := W * V2**T
!
            CALL DTRMM('Right','Upper','Transpose','Unit',N,K,ONE,      &
     &                 V(M-K+1,1),Ldv,Work,Ldwork)
!
!              C2 := C2 - W**T
!
            DO j = 1 , K
               DO i = 1 , N
                  C(M-K+j,i) = C(M-K+j,i) - Work(i,j)
               ENDDO
            ENDDO
!
         ELSEIF ( LSAME(Side,'R') ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
            DO j = 1 , K
               CALL DCOPY(M,C(1,N-K+j),1,Work(1,j),1)
            ENDDO
!
!              W := W * V2
!
            CALL DTRMM('Right','Upper','No transpose','Unit',M,K,ONE,   &
     &                 V(N-K+1,1),Ldv,Work,Ldwork)
!
!                 W := W + C1 * V1
!
            IF ( N>K ) CALL DGEMM('No transpose','No transpose',M,K,N-K,&
     &                            ONE,C,Ldc,V,Ldv,ONE,Work,Ldwork)
!
!              W := W * T  or  W * T**T
!
            CALL DTRMM('Right','Lower',Trans,'Non-unit',M,K,ONE,T,Ldt,  &
     &                 Work,Ldwork)
!
!              C := C - W * V**T
!
!
!                 C1 := C1 - W * V1**T
!
            IF ( N>K ) CALL DGEMM('No transpose','Transpose',M,N-K,K,   &
     &                            -ONE,Work,Ldwork,V,Ldv,ONE,C,Ldc)
!
!              W := W * V2**T
!
            CALL DTRMM('Right','Upper','Transpose','Unit',M,K,ONE,      &
     &                 V(N-K+1,1),Ldv,Work,Ldwork)
!
!              C2 := C2 - W
!
            DO j = 1 , K
               DO i = 1 , M
                  C(i,N-K+j) = C(i,N-K+j) - Work(i,j)
               ENDDO
            ENDDO
         ENDIF
!
      ELSEIF ( LSAME(Storev,'R') ) THEN
!
         IF ( LSAME(Direct,'F') ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
            IF ( LSAME(Side,'L') ) THEN
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
!
!              W := C1**T
!
               DO j = 1 , K
                  CALL DCOPY(N,C(j,1),Ldc,Work(1,j),1)
               ENDDO
!
!              W := W * V1**T
!
               CALL DTRMM('Right','Upper','Transpose','Unit',N,K,ONE,V, &
     &                    Ldv,Work,Ldwork)
!
!                 W := W + C2**T * V2**T
!
               IF ( M>K ) CALL DGEMM('Transpose','Transpose',N,K,M-K,   &
     &                               ONE,C(K+1,1),Ldc,V(1,K+1),Ldv,ONE, &
     &                               Work,Ldwork)
!
!              W := W * T**T  or  W * T
!
               CALL DTRMM('Right','Upper',transt,'Non-unit',N,K,ONE,T,  &
     &                    Ldt,Work,Ldwork)
!
!              C := C - V**T * W**T
!
!
!                 C2 := C2 - V2**T * W**T
!
               IF ( M>K ) CALL DGEMM('Transpose','Transpose',M-K,N,K,   &
     &                               -ONE,V(1,K+1),Ldv,Work,Ldwork,ONE, &
     &                               C(K+1,1),Ldc)
!
!              W := W * V1
!
               CALL DTRMM('Right','Upper','No transpose','Unit',N,K,ONE,&
     &                    V,Ldv,Work,Ldwork)
!
!              C1 := C1 - W**T
!
               DO j = 1 , K
                  DO i = 1 , N
                     C(j,i) = C(j,i) - Work(i,j)
                  ENDDO
               ENDDO
!
            ELSEIF ( LSAME(Side,'R') ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!
!              W := C1
!
               DO j = 1 , K
                  CALL DCOPY(M,C(1,j),1,Work(1,j),1)
               ENDDO
!
!              W := W * V1**T
!
               CALL DTRMM('Right','Upper','Transpose','Unit',M,K,ONE,V, &
     &                    Ldv,Work,Ldwork)
!
!                 W := W + C2 * V2**T
!
               IF ( N>K ) CALL DGEMM('No transpose','Transpose',M,K,N-K,&
     &                               ONE,C(1,K+1),Ldc,V(1,K+1),Ldv,ONE, &
     &                               Work,Ldwork)
!
!              W := W * T  or  W * T**T
!
               CALL DTRMM('Right','Upper',Trans,'Non-unit',M,K,ONE,T,   &
     &                    Ldt,Work,Ldwork)
!
!              C := C - W * V
!
!
!                 C2 := C2 - W * V2
!
               IF ( N>K ) CALL DGEMM('No transpose','No transpose',M,   &
     &                               N-K,K,-ONE,Work,Ldwork,V(1,K+1),   &
     &                               Ldv,ONE,C(1,K+1),Ldc)
!
!              W := W * V1
!
               CALL DTRMM('Right','Upper','No transpose','Unit',M,K,ONE,&
     &                    V,Ldv,Work,Ldwork)
!
!              C1 := C1 - W
!
               DO j = 1 , K
                  DO i = 1 , M
                     C(i,j) = C(i,j) - Work(i,j)
                  ENDDO
               ENDDO
!
            ENDIF
!
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
         ELSEIF ( LSAME(Side,'L') ) THEN
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
!
!              W := C2**T
!
            DO j = 1 , K
               CALL DCOPY(N,C(M-K+j,1),Ldc,Work(1,j),1)
            ENDDO
!
!              W := W * V2**T
!
            CALL DTRMM('Right','Lower','Transpose','Unit',N,K,ONE,      &
     &                 V(1,M-K+1),Ldv,Work,Ldwork)
!
!                 W := W + C1**T * V1**T
!
            IF ( M>K ) CALL DGEMM('Transpose','Transpose',N,K,M-K,ONE,C,&
     &                            Ldc,V,Ldv,ONE,Work,Ldwork)
!
!              W := W * T**T  or  W * T
!
            CALL DTRMM('Right','Lower',transt,'Non-unit',N,K,ONE,T,Ldt, &
     &                 Work,Ldwork)
!
!              C := C - V**T * W**T
!
!
!                 C1 := C1 - V1**T * W**T
!
            IF ( M>K ) CALL DGEMM('Transpose','Transpose',M-K,N,K,-ONE, &
     &                            V,Ldv,Work,Ldwork,ONE,C,Ldc)
!
!              W := W * V2
!
            CALL DTRMM('Right','Lower','No transpose','Unit',N,K,ONE,   &
     &                 V(1,M-K+1),Ldv,Work,Ldwork)
!
!              C2 := C2 - W**T
!
            DO j = 1 , K
               DO i = 1 , N
                  C(M-K+j,i) = C(M-K+j,i) - Work(i,j)
               ENDDO
            ENDDO
!
         ELSEIF ( LSAME(Side,'R') ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!
!              W := C2
!
            DO j = 1 , K
               CALL DCOPY(M,C(1,N-K+j),1,Work(1,j),1)
            ENDDO
!
!              W := W * V2**T
!
            CALL DTRMM('Right','Lower','Transpose','Unit',M,K,ONE,      &
     &                 V(1,N-K+1),Ldv,Work,Ldwork)
!
!                 W := W + C1 * V1**T
!
            IF ( N>K ) CALL DGEMM('No transpose','Transpose',M,K,N-K,   &
     &                            ONE,C,Ldc,V,Ldv,ONE,Work,Ldwork)
!
!              W := W * T  or  W * T**T
!
            CALL DTRMM('Right','Lower',Trans,'Non-unit',M,K,ONE,T,Ldt,  &
     &                 Work,Ldwork)
!
!              C := C - W * V
!
!
!                 C1 := C1 - W * V1
!
            IF ( N>K ) CALL DGEMM('No transpose','No transpose',M,N-K,K,&
     &                            -ONE,Work,Ldwork,V,Ldv,ONE,C,Ldc)
!
!              W := W * V2
!
            CALL DTRMM('Right','Lower','No transpose','Unit',M,K,ONE,   &
     &                 V(1,N-K+1),Ldv,Work,Ldwork)
!
!              C1 := C1 - W
!
            DO j = 1 , K
               DO i = 1 , M
                  C(i,N-K+j) = C(i,N-K+j) - Work(i,j)
               ENDDO
            ENDDO
!
!
         ENDIF
      ENDIF
!
!
!     End of DLARFB
!
      END SUBROUTINE DLARFB
