!*==dlascl.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASCL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlascl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlascl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlascl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TYPE
!       INTEGER            INFO, KL, KU, LDA, M, N
!       DOUBLE PRECISION   CFROM, CTO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASCL multiplies the M by N real matrix A by the real scalar
!> CTO/CFROM.  This is done without over/underflow as long as the final
!> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!> A may be full, upper triangular, lower triangular, upper Hessenberg,
!> or banded.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TYPE
!> \verbatim
!>          TYPE is CHARACTER*1
!>          TYPE indices the storage type of the input matrix.
!>          = 'G':  A is a full matrix.
!>          = 'L':  A is a lower triangular matrix.
!>          = 'U':  A is an upper triangular matrix.
!>          = 'H':  A is an upper Hessenberg matrix.
!>          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the lower
!>                  half stored.
!>          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the upper
!>                  half stored.
!>          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!>                  bandwidth KU. See DGBTRF for storage details.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] CFROM
!> \verbatim
!>          CFROM is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] CTO
!> \verbatim
!>          CTO is DOUBLE PRECISION
!>
!>          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!>          without over/underflow if the final result CTO*A(I,J)/CFROM
!>          can be represented without over/underflow.  CFROM must be
!>          nonzero.
!> \endverbatim
!>
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!>          storage type.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If TYPE = 'G', 'L', 'U', 'H', LDA >= max(1,M);
!>             TYPE = 'B', LDA >= KL+1;
!>             TYPE = 'Q', LDA >= KU+1;
!>             TYPE = 'Z', LDA >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          0  - successful exit
!>          <0 - if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLASCL(Type,Kl,Ku,Cfrom,Cto,M,N,A,Lda,Info)
      USE F77KINDS                        
      USE S_DISNAN
      USE S_DLAMCH
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DLASCL152
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Type
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL(R8KIND) :: Cfrom
      REAL(R8KIND) :: Cto
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: bignum , cfrom1 , cfromc , cto1 , ctoc , mul ,    &
     &                smlnum
      LOGICAL :: done
      INTEGER :: i , itype , j , k1 , k2 , k3 , k4
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
!     .. Intrinsic Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
!
      IF ( LSAME(Type,'G') ) THEN
         itype = 0
      ELSEIF ( LSAME(Type,'L') ) THEN
         itype = 1
      ELSEIF ( LSAME(Type,'U') ) THEN
         itype = 2
      ELSEIF ( LSAME(Type,'H') ) THEN
         itype = 3
      ELSEIF ( LSAME(Type,'B') ) THEN
         itype = 4
      ELSEIF ( LSAME(Type,'Q') ) THEN
         itype = 5
      ELSEIF ( LSAME(Type,'Z') ) THEN
         itype = 6
      ELSE
         itype = -1
      ENDIF
!
      IF ( itype==-1 ) THEN
         Info = -1
      ELSEIF ( Cfrom==ZERO .OR. DISNAN(Cfrom) ) THEN
         Info = -4
      ELSEIF ( DISNAN(Cto) ) THEN
         Info = -5
      ELSEIF ( M<0 ) THEN
         Info = -6
      ELSEIF ( N<0 .OR. (itype==4 .AND. N/=M) .OR. (itype==5 .AND. N/=M)&
     &         ) THEN
         Info = -7
      ELSEIF ( itype<=3 .AND. Lda<MAX(1,M) ) THEN
         Info = -9
      ELSEIF ( itype>=4 ) THEN
         IF ( Kl<0 .OR. Kl>MAX(M-1,0) ) THEN
            Info = -2
         ELSEIF ( Ku<0 .OR. Ku>MAX(N-1,0) .OR.                          &
     &            ((itype==4 .OR. itype==5) .AND. Kl/=Ku) ) THEN
            Info = -3
         ELSEIF ( (itype==4 .AND. Lda<Kl+1) .OR.                        &
     &            (itype==5 .AND. Lda<Ku+1) .OR.                        &
     &            (itype==6 .AND. Lda<2*Kl+Ku+1) ) THEN
            Info = -9
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DLASCL',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. M==0 ) RETURN
!
!     Get machine parameters
!
      smlnum = DLAMCH('S')
      bignum = ONE/smlnum
!
      cfromc = Cfrom
      ctoc = Cto
      DO
!
         cfrom1 = cfromc*smlnum
         IF ( cfrom1==cfromc ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
            mul = ctoc/cfromc
            done = .TRUE.
            cto1 = ctoc
         ELSE
            cto1 = ctoc/bignum
            IF ( cto1==ctoc ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
               mul = ctoc
               done = .TRUE.
               cfromc = ONE
            ELSEIF ( ABS(cfrom1)>ABS(ctoc) .AND. ctoc/=ZERO ) THEN
               mul = smlnum
               done = .FALSE.
               cfromc = cfrom1
            ELSEIF ( ABS(cto1)>ABS(cfromc) ) THEN
               mul = bignum
               done = .FALSE.
               ctoc = cto1
            ELSE
               mul = ctoc/cfromc
               done = .TRUE.
            ENDIF
         ENDIF
!
         IF ( itype==0 ) THEN
!
!        Full matrix
!
            DO j = 1 , N
               DO i = 1 , M
                  A(i,j) = A(i,j)*mul
               ENDDO
            ENDDO
!
         ELSEIF ( itype==1 ) THEN
!
!        Lower triangular matrix
!
            DO j = 1 , N
               DO i = j , M
                  A(i,j) = A(i,j)*mul
               ENDDO
            ENDDO
!
         ELSEIF ( itype==2 ) THEN
!
!        Upper triangular matrix
!
            DO j = 1 , N
               DO i = 1 , MIN(j,M)
                  A(i,j) = A(i,j)*mul
               ENDDO
            ENDDO
!
         ELSEIF ( itype==3 ) THEN
!
!        Upper Hessenberg matrix
!
            DO j = 1 , N
               DO i = 1 , MIN(j+1,M)
                  A(i,j) = A(i,j)*mul
               ENDDO
            ENDDO
!
         ELSEIF ( itype==4 ) THEN
!
!        Lower half of a symmetric band matrix
!
            k3 = Kl + 1
            k4 = N + 1
            DO j = 1 , N
               DO i = 1 , MIN(k3,k4-j)
                  A(i,j) = A(i,j)*mul
               ENDDO
            ENDDO
!
         ELSEIF ( itype==5 ) THEN
!
!        Upper half of a symmetric band matrix
!
            k1 = Ku + 2
            k3 = Ku + 1
            DO j = 1 , N
               DO i = MAX(k1-j,1) , k3
                  A(i,j) = A(i,j)*mul
               ENDDO
            ENDDO
!
         ELSEIF ( itype==6 ) THEN
!
!        Band matrix
!
            k1 = Kl + Ku + 2
            k2 = Kl + 1
            k3 = 2*Kl + Ku + 1
            k4 = Kl + Ku + 1 + M
            DO j = 1 , N
               DO i = MAX(k1-j,k2) , MIN(k3,k4-j)
                  A(i,j) = A(i,j)*mul
               ENDDO
            ENDDO
!
         ENDIF
!
         IF ( done ) EXIT
      ENDDO
!
!
!     End of DLASCL
!
      END SUBROUTINE DLASCL
