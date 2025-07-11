!*==zlarft.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLARFT forms the triangular factor T of a block reflector H = I - vtvH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLARFT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarft.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarft.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarft.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, STOREV
!       INTEGER            K, LDT, LDV, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARFT forms the triangular factor T of a complex block reflector H
!> of order n, which is defined as a product of k elementary reflectors.
!>
!> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!>
!> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!>
!> If STOREV = 'C', the vector which defines the elementary reflector
!> H(i) is stored in the i-th column of the array V, and
!>
!>    H  =  I - V * T * V**H
!>
!> If STOREV = 'R', the vector which defines the elementary reflector
!> H(i) is stored in the i-th row of the array V, and
!>
!>    H  =  I - V**H * T * V
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Specifies the order in which the elementary reflectors are
!>          multiplied to form the block reflector:
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Specifies how the vectors which define the elementary
!>          reflectors are stored (see also Further Details):
!>          = 'C': columnwise
!>          = 'R': rowwise
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the block reflector H. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The order of the triangular factor T (= the number of
!>          elementary reflectors). K >= 1.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension
!>                               (LDV,K) if STOREV = 'C'
!>                               (LDV,N) if STOREV = 'R'
!>          The matrix V. See further details.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,K)
!>          The k by k triangular factor T of the block reflector.
!>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!>          lower triangular. The rest of the array is not used.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= K.
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
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The shape of the matrix V and the storage of the vectors which define
!>  the H(i) is best illustrated by the following example with n = 5 and
!>  k = 3. The elements equal to 1 are not stored.
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
      SUBROUTINE ZLARFT(Direct,Storev,N,K,V,Ldv,Tau,T,Ldt)
      IMPLICIT NONE
!*--ZLARFT167
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Direct , Storev
      INTEGER K , Ldt , Ldv , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 T(Ldt,*) , Tau(*) , V(Ldv,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE , ZERO
      PARAMETER (ONE=(1.0D+0,0.0D+0),ZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , j , prevlastv , lastv
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMV , ZTRMV , ZGEMM
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( LSAME(Direct,'F') ) THEN
         prevlastv = N
         DO i = 1 , K
            prevlastv = MAX(prevlastv,i)
            IF ( Tau(i)==ZERO ) THEN
!
!              H(i)  =  I
!
               DO j = 1 , i
                  T(j,i) = ZERO
               ENDDO
            ELSE
!
!              general case
!
               IF ( LSAME(Storev,'C') ) THEN
!                 Skip any trailing zeros.
                  DO lastv = N , i + 1 , -1
                     IF ( V(lastv,i)/=ZERO ) EXIT
                  ENDDO
                  DO j = 1 , i - 1
                     T(j,i) = -Tau(i)*CONJG(V(i,j))
                  ENDDO
                  j = MIN(lastv,prevlastv)
!
!                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**H * V(i:j,i)
!
                  CALL ZGEMV('Conjugate transpose',j-i,i-1,-Tau(i),     &
     &                       V(i+1,1),Ldv,V(i+1,i),1,ONE,T(1,i),1)
               ELSE
!                 Skip any trailing zeros.
                  DO lastv = N , i + 1 , -1
                     IF ( V(i,lastv)/=ZERO ) EXIT
                  ENDDO
                  DO j = 1 , i - 1
                     T(j,i) = -Tau(i)*V(j,i)
                  ENDDO
                  j = MIN(lastv,prevlastv)
!
!                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**H
!
                  CALL ZGEMM('N','C',i-1,1,j-i,-Tau(i),V(1,i+1),Ldv,    &
     &                       V(i,i+1),Ldv,ONE,T(1,i),Ldt)
               ENDIF
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
               CALL ZTRMV('Upper','No transpose','Non-unit',i-1,T,Ldt,  &
     &                    T(1,i),1)
               T(i,i) = Tau(i)
               IF ( i>1 ) THEN
                  prevlastv = MAX(prevlastv,lastv)
               ELSE
                  prevlastv = lastv
               ENDIF
            ENDIF
         ENDDO
      ELSE
         prevlastv = 1
         DO i = K , 1 , -1
            IF ( Tau(i)==ZERO ) THEN
!
!              H(i)  =  I
!
               DO j = i , K
                  T(j,i) = ZERO
               ENDDO
            ELSE
!
!              general case
!
               IF ( i<K ) THEN
                  IF ( LSAME(Storev,'C') ) THEN
!                    Skip any leading zeros.
                     DO lastv = 1 , i - 1
                        IF ( V(lastv,i)/=ZERO ) EXIT
                     ENDDO
                     DO j = i + 1 , K
                        T(j,i) = -Tau(i)*CONJG(V(N-K+i,j))
                     ENDDO
                     j = MAX(lastv,prevlastv)
!
!                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**H * V(j:n-k+i,i)
!
                     CALL ZGEMV('Conjugate transpose',N-K+i-j,K-i,      &
     &                          -Tau(i),V(j,i+1),Ldv,V(j,i),1,ONE,      &
     &                          T(i+1,i),1)
                  ELSE
!                    Skip any leading zeros.
                     DO lastv = 1 , i - 1
                        IF ( V(i,lastv)/=ZERO ) EXIT
                     ENDDO
                     DO j = i + 1 , K
                        T(j,i) = -Tau(i)*V(j,N-K+i)
                     ENDDO
                     j = MAX(lastv,prevlastv)
!
!                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**H
!
                     CALL ZGEMM('N','C',K-i,1,N-K+i-j,-Tau(i),V(i+1,j), &
     &                          Ldv,V(i,j),Ldv,ONE,T(i+1,i),Ldt)
                  ENDIF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
                  CALL ZTRMV('Lower','No transpose','Non-unit',K-i,     &
     &                       T(i+1,i+1),Ldt,T(i+1,i),1)
                  IF ( i>1 ) THEN
                     prevlastv = MIN(prevlastv,lastv)
                  ELSE
                     prevlastv = lastv
                  ENDIF
               ENDIF
               T(i,i) = Tau(i)
            ENDIF
         ENDDO
      ENDIF
!
!     End of ZLARFT
!
      END SUBROUTINE ZLARFT
