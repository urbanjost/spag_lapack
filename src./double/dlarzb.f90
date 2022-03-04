!*==dlarzb.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARZB applies a block reflector or its transpose to a general matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARZB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarzb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarzb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarzb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V,
!                          LDV, T, LDT, C, LDC, WORK, LDWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, SIDE, STOREV, TRANS
!       INTEGER            K, L, LDC, LDT, LDV, LDWORK, M, N
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
!> DLARZB applies a real block reflector H or its transpose H**T to
!> a real distributed M-by-N  C from the left or the right.
!>
!> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.
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
!>          = 'C': apply H**T (Transpose)
!> \endverbatim
!>
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Indicates how H is formed from a product of elementary
!>          reflectors
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Indicates how the vectors which define the elementary
!>          reflectors are stored:
!>          = 'C': Columnwise                        (not supported yet)
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
!> \endverbatim
!>
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>          The number of columns of the matrix V containing the
!>          meaningful part of the Householder reflectors.
!>          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (LDV,NV).
!>          If STOREV = 'C', NV = K; if STOREV = 'R', NV = L.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If STOREV = 'C', LDV >= L; if STOREV = 'R', LDV >= K.
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
!>          The leading dimension of the array T. LDT >= K.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
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
!> \date December 2016
!
!> \ingroup doubleOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLARZB(Side,Trans,Direct,Storev,M,N,K,L,V,Ldv,T,Ldt,C, &
     &                  Ldc,Work,Ldwork)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DGEMM
      USE S_DTRMM
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DLARZB193
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
      INTEGER :: L
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
      INTEGER :: i , info , j
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
!     Check for currently supported options
!
      info = 0
      IF ( .NOT.LSAME(Direct,'B') ) THEN
         info = -3
      ELSEIF ( .NOT.LSAME(Storev,'R') ) THEN
         info = -4
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('DLARZB',-info)
         RETURN
      ENDIF
!
      IF ( LSAME(Trans,'N') ) THEN
         transt = 'T'
      ELSE
         transt = 'N'
      ENDIF
!
      IF ( LSAME(Side,'L') ) THEN
!
!        Form  H * C  or  H**T * C
!
!        W( 1:n, 1:k ) = C( 1:k, 1:n )**T
!
         DO j = 1 , K
            CALL DCOPY(N,C(j,1),Ldc,Work(1,j),1)
         ENDDO
!
!        W( 1:n, 1:k ) = W( 1:n, 1:k ) + ...
!                        C( m-l+1:m, 1:n )**T * V( 1:k, 1:l )**T
!
         IF ( L>0 ) CALL DGEMM('Transpose','Transpose',N,K,L,ONE,       &
     &                         C(M-L+1,1),Ldc,V,Ldv,ONE,Work,Ldwork)
!
!        W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T
!
         CALL DTRMM('Right','Lower',transt,'Non-unit',N,K,ONE,T,Ldt,    &
     &              Work,Ldwork)
!
!        C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**T
!
         DO j = 1 , N
            DO i = 1 , K
               C(i,j) = C(i,j) - Work(j,i)
            ENDDO
         ENDDO
!
!        C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
!                            V( 1:k, 1:l )**T * W( 1:n, 1:k )**T
!
         IF ( L>0 ) CALL DGEMM('Transpose','Transpose',L,N,K,-ONE,V,Ldv,&
     &                         Work,Ldwork,ONE,C(M-L+1,1),Ldc)
!
      ELSEIF ( LSAME(Side,'R') ) THEN
!
!        Form  C * H  or  C * H**T
!
!        W( 1:m, 1:k ) = C( 1:m, 1:k )
!
         DO j = 1 , K
            CALL DCOPY(M,C(1,j),1,Work(1,j),1)
         ENDDO
!
!        W( 1:m, 1:k ) = W( 1:m, 1:k ) + ...
!                        C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**T
!
         IF ( L>0 ) CALL DGEMM('No transpose','Transpose',M,K,L,ONE,    &
     &                         C(1,N-L+1),Ldc,V,Ldv,ONE,Work,Ldwork)
!
!        W( 1:m, 1:k ) = W( 1:m, 1:k ) * T  or  W( 1:m, 1:k ) * T**T
!
         CALL DTRMM('Right','Lower',Trans,'Non-unit',M,K,ONE,T,Ldt,Work,&
     &              Ldwork)
!
!        C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k )
!
         DO j = 1 , K
            DO i = 1 , M
               C(i,j) = C(i,j) - Work(i,j)
            ENDDO
         ENDDO
!
!        C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
!                            W( 1:m, 1:k ) * V( 1:k, 1:l )
!
         IF ( L>0 ) CALL DGEMM('No transpose','No transpose',M,L,K,-ONE,&
     &                         Work,Ldwork,V,Ldv,ONE,C(1,N-L+1),Ldc)
!
      ENDIF
!
!
!     End of DLARZB
!
      END SUBROUTINE DLARZB
