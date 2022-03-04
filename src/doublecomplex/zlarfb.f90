!*==zlarfb.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLARFB applies a block reflector or its conjugate-transpose to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLARFB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
!                          T, LDT, C, LDC, WORK, LDWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, SIDE, STOREV, TRANS
!       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ),
!      $                   WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARFB applies a complex block reflector H or its transpose H**H to a
!> complex M-by-N matrix C, from either the left or the right.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply H or H**H from the Left
!>          = 'R': apply H or H**H from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply H (No transpose)
!>          = 'C': apply H**H (Conjugate transpose)
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
!>          V is COMPLEX*16 array, dimension
!>                                (LDV,K) if STOREV = 'C'
!>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!>          See Further Details.
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
!>          T is COMPLEX*16 array, dimension (LDT,K)
!>          The triangular K-by-K matrix T in the representation of the
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
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by H*C or H**H*C or C*H or C*H**H.
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
!>          WORK is COMPLEX*16 array, dimension (LDWORK,K)
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
!> \ingroup complex16OTHERauxiliary
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
      SUBROUTINE ZLARFB(Side,Trans,Direct,Storev,M,N,K,V,Ldv,T,Ldt,C,   &
     &                  Ldc,Work,Ldwork)
      IMPLICIT NONE
!*--ZLARFB201
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2013
!
!     .. Scalar Arguments ..
      CHARACTER Direct , Side , Storev , Trans
      INTEGER K , Ldc , Ldt , Ldv , Ldwork , M , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 C(Ldc,*) , T(Ldt,*) , V(Ldv,*) , Work(Ldwork,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      CHARACTER transt
      INTEGER i , j
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL ZCOPY , ZGEMM , ZLACGV , ZTRMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( M<=0 .OR. N<=0 ) RETURN
!
      IF ( LSAME(Trans,'N') ) THEN
         transt = 'C'
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
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
!
!              W := C1**H
!
               DO j = 1 , K
                  CALL ZCOPY(N,C(j,1),Ldc,Work(1,j),1)
                  CALL ZLACGV(N,Work(1,j),1)
               ENDDO
!
!              W := W * V1
!
               CALL ZTRMM('Right','Lower','No transpose','Unit',N,K,ONE,&
     &                    V,Ldv,Work,Ldwork)
!
!                 W := W + C2**H * V2
!
               IF ( M>K ) CALL ZGEMM('Conjugate transpose',             &
     &                               'No transpose',N,K,M-K,ONE,C(K+1,1)&
     &                               ,Ldc,V(K+1,1),Ldv,ONE,Work,Ldwork)
!
!              W := W * T**H  or  W * T
!
               CALL ZTRMM('Right','Upper',transt,'Non-unit',N,K,ONE,T,  &
     &                    Ldt,Work,Ldwork)
!
!              C := C - V * W**H
!
!
!                 C2 := C2 - V2 * W**H
!
               IF ( M>K ) CALL ZGEMM('No transpose',                    &
     &                               'Conjugate transpose',M-K,N,K,-ONE,&
     &                               V(K+1,1),Ldv,Work,Ldwork,ONE,      &
     &                               C(K+1,1),Ldc)
!
!              W := W * V1**H
!
               CALL ZTRMM('Right','Lower','Conjugate transpose','Unit', &
     &                    N,K,ONE,V,Ldv,Work,Ldwork)
!
!              C1 := C1 - W**H
!
               DO j = 1 , K
                  DO i = 1 , N
                     C(j,i) = C(j,i) - DCONJG(Work(i,j))
                  ENDDO
               ENDDO
!
            ELSEIF ( LSAME(Side,'R') ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO j = 1 , K
                  CALL ZCOPY(M,C(1,j),1,Work(1,j),1)
               ENDDO
!
!              W := W * V1
!
               CALL ZTRMM('Right','Lower','No transpose','Unit',M,K,ONE,&
     &                    V,Ldv,Work,Ldwork)
!
!                 W := W + C2 * V2
!
               IF ( N>K ) CALL ZGEMM('No transpose','No transpose',M,K, &
     &                               N-K,ONE,C(1,K+1),Ldc,V(K+1,1),Ldv, &
     &                               ONE,Work,Ldwork)
!
!              W := W * T  or  W * T**H
!
               CALL ZTRMM('Right','Upper',Trans,'Non-unit',M,K,ONE,T,   &
     &                    Ldt,Work,Ldwork)
!
!              C := C - W * V**H
!
!
!                 C2 := C2 - W * V2**H
!
               IF ( N>K ) CALL ZGEMM('No transpose',                    &
     &                               'Conjugate transpose',M,N-K,K,-ONE,&
     &                               Work,Ldwork,V(K+1,1),Ldv,ONE,      &
     &                               C(1,K+1),Ldc)
!
!              W := W * V1**H
!
               CALL ZTRMM('Right','Lower','Conjugate transpose','Unit', &
     &                    M,K,ONE,V,Ldv,Work,Ldwork)
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
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
!
!              W := C2**H
!
            DO j = 1 , K
               CALL ZCOPY(N,C(M-K+j,1),Ldc,Work(1,j),1)
               CALL ZLACGV(N,Work(1,j),1)
            ENDDO
!
!              W := W * V2
!
            CALL ZTRMM('Right','Upper','No transpose','Unit',N,K,ONE,   &
     &                 V(M-K+1,1),Ldv,Work,Ldwork)
!
!                 W := W + C1**H * V1
!
            IF ( M>K ) CALL ZGEMM('Conjugate transpose','No transpose', &
     &                            N,K,M-K,ONE,C,Ldc,V,Ldv,ONE,Work,     &
     &                            Ldwork)
!
!              W := W * T**H  or  W * T
!
            CALL ZTRMM('Right','Lower',transt,'Non-unit',N,K,ONE,T,Ldt, &
     &                 Work,Ldwork)
!
!              C := C - V * W**H
!
!
!                 C1 := C1 - V1 * W**H
!
            IF ( M>K ) CALL ZGEMM('No transpose','Conjugate transpose', &
     &                            M-K,N,K,-ONE,V,Ldv,Work,Ldwork,ONE,C, &
     &                            Ldc)
!
!              W := W * V2**H
!
            CALL ZTRMM('Right','Upper','Conjugate transpose','Unit',N,K,&
     &                 ONE,V(M-K+1,1),Ldv,Work,Ldwork)
!
!              C2 := C2 - W**H
!
            DO j = 1 , K
               DO i = 1 , N
                  C(M-K+j,i) = C(M-K+j,i) - DCONJG(Work(i,j))
               ENDDO
            ENDDO
!
         ELSEIF ( LSAME(Side,'R') ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
            DO j = 1 , K
               CALL ZCOPY(M,C(1,N-K+j),1,Work(1,j),1)
            ENDDO
!
!              W := W * V2
!
            CALL ZTRMM('Right','Upper','No transpose','Unit',M,K,ONE,   &
     &                 V(N-K+1,1),Ldv,Work,Ldwork)
!
!                 W := W + C1 * V1
!
            IF ( N>K ) CALL ZGEMM('No transpose','No transpose',M,K,N-K,&
     &                            ONE,C,Ldc,V,Ldv,ONE,Work,Ldwork)
!
!              W := W * T  or  W * T**H
!
            CALL ZTRMM('Right','Lower',Trans,'Non-unit',M,K,ONE,T,Ldt,  &
     &                 Work,Ldwork)
!
!              C := C - W * V**H
!
!
!                 C1 := C1 - W * V1**H
!
            IF ( N>K ) CALL ZGEMM('No transpose','Conjugate transpose', &
     &                            M,N-K,K,-ONE,Work,Ldwork,V,Ldv,ONE,C, &
     &                            Ldc)
!
!              W := W * V2**H
!
            CALL ZTRMM('Right','Upper','Conjugate transpose','Unit',M,K,&
     &                 ONE,V(N-K+1,1),Ldv,Work,Ldwork)
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
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
!
!              W := C1**H
!
               DO j = 1 , K
                  CALL ZCOPY(N,C(j,1),Ldc,Work(1,j),1)
                  CALL ZLACGV(N,Work(1,j),1)
               ENDDO
!
!              W := W * V1**H
!
               CALL ZTRMM('Right','Upper','Conjugate transpose','Unit', &
     &                    N,K,ONE,V,Ldv,Work,Ldwork)
!
!                 W := W + C2**H * V2**H
!
               IF ( M>K ) CALL ZGEMM('Conjugate transpose',             &
     &                               'Conjugate transpose',N,K,M-K,ONE, &
     &                               C(K+1,1),Ldc,V(1,K+1),Ldv,ONE,Work,&
     &                               Ldwork)
!
!              W := W * T**H  or  W * T
!
               CALL ZTRMM('Right','Upper',transt,'Non-unit',N,K,ONE,T,  &
     &                    Ldt,Work,Ldwork)
!
!              C := C - V**H * W**H
!
!
!                 C2 := C2 - V2**H * W**H
!
               IF ( M>K ) CALL ZGEMM('Conjugate transpose',             &
     &                               'Conjugate transpose',M-K,N,K,-ONE,&
     &                               V(1,K+1),Ldv,Work,Ldwork,ONE,      &
     &                               C(K+1,1),Ldc)
!
!              W := W * V1
!
               CALL ZTRMM('Right','Upper','No transpose','Unit',N,K,ONE,&
     &                    V,Ldv,Work,Ldwork)
!
!              C1 := C1 - W**H
!
               DO j = 1 , K
                  DO i = 1 , N
                     C(j,i) = C(j,i) - DCONJG(Work(i,j))
                  ENDDO
               ENDDO
!
            ELSEIF ( LSAME(Side,'R') ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
!              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
!
!              W := C1
!
               DO j = 1 , K
                  CALL ZCOPY(M,C(1,j),1,Work(1,j),1)
               ENDDO
!
!              W := W * V1**H
!
               CALL ZTRMM('Right','Upper','Conjugate transpose','Unit', &
     &                    M,K,ONE,V,Ldv,Work,Ldwork)
!
!                 W := W + C2 * V2**H
!
               IF ( N>K ) CALL ZGEMM('No transpose',                    &
     &                               'Conjugate transpose',M,K,N-K,ONE, &
     &                               C(1,K+1),Ldc,V(1,K+1),Ldv,ONE,Work,&
     &                               Ldwork)
!
!              W := W * T  or  W * T**H
!
               CALL ZTRMM('Right','Upper',Trans,'Non-unit',M,K,ONE,T,   &
     &                    Ldt,Work,Ldwork)
!
!              C := C - W * V
!
!
!                 C2 := C2 - W * V2
!
               IF ( N>K ) CALL ZGEMM('No transpose','No transpose',M,   &
     &                               N-K,K,-ONE,Work,Ldwork,V(1,K+1),   &
     &                               Ldv,ONE,C(1,K+1),Ldc)
!
!              W := W * V1
!
               CALL ZTRMM('Right','Upper','No transpose','Unit',M,K,ONE,&
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
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
!
!              W := C2**H
!
            DO j = 1 , K
               CALL ZCOPY(N,C(M-K+j,1),Ldc,Work(1,j),1)
               CALL ZLACGV(N,Work(1,j),1)
            ENDDO
!
!              W := W * V2**H
!
            CALL ZTRMM('Right','Lower','Conjugate transpose','Unit',N,K,&
     &                 ONE,V(1,M-K+1),Ldv,Work,Ldwork)
!
!                 W := W + C1**H * V1**H
!
            IF ( M>K ) CALL ZGEMM('Conjugate transpose',                &
     &                            'Conjugate transpose',N,K,M-K,ONE,C,  &
     &                            Ldc,V,Ldv,ONE,Work,Ldwork)
!
!              W := W * T**H  or  W * T
!
            CALL ZTRMM('Right','Lower',transt,'Non-unit',N,K,ONE,T,Ldt, &
     &                 Work,Ldwork)
!
!              C := C - V**H * W**H
!
!
!                 C1 := C1 - V1**H * W**H
!
            IF ( M>K ) CALL ZGEMM('Conjugate transpose',                &
     &                            'Conjugate transpose',M-K,N,K,-ONE,V, &
     &                            Ldv,Work,Ldwork,ONE,C,Ldc)
!
!              W := W * V2
!
            CALL ZTRMM('Right','Lower','No transpose','Unit',N,K,ONE,   &
     &                 V(1,M-K+1),Ldv,Work,Ldwork)
!
!              C2 := C2 - W**H
!
            DO j = 1 , K
               DO i = 1 , N
                  C(M-K+j,i) = C(M-K+j,i) - DCONJG(Work(i,j))
               ENDDO
            ENDDO
!
         ELSEIF ( LSAME(Side,'R') ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
!              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
!
!              W := C2
!
            DO j = 1 , K
               CALL ZCOPY(M,C(1,N-K+j),1,Work(1,j),1)
            ENDDO
!
!              W := W * V2**H
!
            CALL ZTRMM('Right','Lower','Conjugate transpose','Unit',M,K,&
     &                 ONE,V(1,N-K+1),Ldv,Work,Ldwork)
!
!                 W := W + C1 * V1**H
!
            IF ( N>K ) CALL ZGEMM('No transpose','Conjugate transpose', &
     &                            M,K,N-K,ONE,C,Ldc,V,Ldv,ONE,Work,     &
     &                            Ldwork)
!
!              W := W * T  or  W * T**H
!
            CALL ZTRMM('Right','Lower',Trans,'Non-unit',M,K,ONE,T,Ldt,  &
     &                 Work,Ldwork)
!
!              C := C - W * V
!
!
!                 C1 := C1 - W * V1
!
            IF ( N>K ) CALL ZGEMM('No transpose','No transpose',M,N-K,K,&
     &                            -ONE,Work,Ldwork,V,Ldv,ONE,C,Ldc)
!
!              W := W * V2
!
            CALL ZTRMM('Right','Lower','No transpose','Unit',M,K,ONE,   &
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
!     End of ZLARFB
!
      END SUBROUTINE ZLARFB
