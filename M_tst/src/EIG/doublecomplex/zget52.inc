!*==zget52.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b zget52
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zget52( LEFT, N, A, LDA, B, LDB, E, LDE, ALPHA, BETA,
!                          WORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       LOGICAL            LEFT
!       INTEGER            LDA, LDB, LDE, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RESULT( 2 ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
!      $                   BETA( * ), E( LDE, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zget52  does an eigenvector check for the generalized eigenvalue
!> problem.
!>
!> The basic test for right eigenvectors is:
!>
!>                           | b(i) A E(i) -  a(i) B E(i) |
!>         RESULT(1) = max   -------------------------------
!>                      i    n ulp max( |b(i) A|, |a(i) B| )
!>
!> using the 1-norm.  Here, a(i)/b(i) = w is the i-th generalized
!> eigenvalue of A - w B, or, equivalently, b(i)/a(i) = m is the i-th
!> generalized eigenvalue of m A - B.
!>
!>                         H   H  _      _
!> For left eigenvectors, A , B , a, and b  are used.
!>
!> zget52 also tests the normalization of E.  Each eigenvector is
!> supposed to be normalized so that the maximum "absolute value"
!> of its elements is 1, where in this case, "absolute value"
!> of a complex value x is  |Re(x)| + |Im(x)| ; let us call this
!> maximum "absolute value" norm of a vector v  M(v).
!> If a(i)=b(i)=0, then the eigenvector is set to be the jth coordinate
!> vector. The normalization test is:
!>
!>         RESULT(2) =      max       | M(v(i)) - 1 | / ( n ulp )
!>                    eigenvectors v(i)
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] LEFT
!> \verbatim
!>          LEFT is LOGICAL
!>          =.TRUE.:  The eigenvectors in the columns of E are assumed
!>                    to be *left* eigenvectors.
!>          =.FALSE.: The eigenvectors in the columns of E are assumed
!>                    to be *right* eigenvectors.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrices.  If it is zero, zget52 does
!>          nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA, N)
!>          The matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 1
!>          and at least N.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB, N)
!>          The matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of B.  It must be at least 1
!>          and at least N.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX*16 array, dimension (LDE, N)
!>          The matrix of eigenvectors.  It must be O( 1 ).
!> \endverbatim
!>
!> \param[in] LDE
!> \verbatim
!>          LDE is INTEGER
!>          The leading dimension of E.  It must be at least 1 and at
!>          least N.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16 array, dimension (N)
!>          The values a(i) as described above, which, along with b(i),
!>          define the generalized eigenvalues.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX*16 array, dimension (N)
!>          The values b(i) as described above, which, along with a(i),
!>          define the generalized eigenvalues.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N**2)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (2)
!>          The values computed by the test described above.  If A E or
!>          B E is likely to overflow, then RESULT(1:2) is set to
!>          10 / ulp.
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE zget52(Left,N,A,Lda,B,Ldb,E,Lde,Alpha,Beta,Work,Rwork, &
     &                  Result)
      IMPLICIT NONE
!*--ZGET52166
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Left
      INTEGER Lda , Ldb , Lde , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Result(2) , Rwork(*)
      COMPLEX*16 A(Lda,*) , Alpha(*) , B(Ldb,*) , Beta(*) , E(Lde,*) ,  &
     &           Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      CHARACTER normab , trans
      INTEGER j , jvec
      DOUBLE PRECISION abmax , alfmax , anorm , betmax , bnorm , enorm ,&
     &                 enrmer , errnrm , safmax , safmin , scale ,      &
     &                 temp1 , ulp
      COMPLEX*16 acoeff , alphai , bcoeff , betai , x
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL DLAMCH , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DCONJG , DIMAG , MAX
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION ABS1
!     ..
!     .. Statement Function definitions ..
      ABS1(x) = ABS(DBLE(x)) + ABS(DIMAG(x))
!     ..
!     .. Executable Statements ..
!
      Result(1) = ZERO
      Result(2) = ZERO
      IF ( N<=0 ) RETURN
!
      safmin = DLAMCH('Safe minimum')
      safmax = ONE/safmin
      ulp = DLAMCH('Epsilon')*DLAMCH('Base')
!
      IF ( Left ) THEN
         trans = 'C'
         normab = 'I'
      ELSE
         trans = 'N'
         normab = 'O'
      ENDIF
!
!     Norm of A, B, and E:
!
      anorm = MAX(ZLANGE(normab,N,N,A,Lda,Rwork),safmin)
      bnorm = MAX(ZLANGE(normab,N,N,B,Ldb,Rwork),safmin)
      enorm = MAX(ZLANGE('O',N,N,E,Lde,Rwork),ulp)
      alfmax = safmax/MAX(ONE,bnorm)
      betmax = safmax/MAX(ONE,anorm)
!
!     Compute error matrix.
!     Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B| |b(i) A| )
!
      DO jvec = 1 , N
         alphai = Alpha(jvec)
         betai = Beta(jvec)
         abmax = MAX(ABS1(alphai),ABS1(betai))
         IF ( ABS1(alphai)>alfmax .OR. ABS1(betai)>betmax .OR.          &
     &        abmax<ONE ) THEN
            scale = ONE/MAX(abmax,safmin)
            alphai = scale*alphai
            betai = scale*betai
         ENDIF
         scale = ONE/MAX(ABS1(alphai)*bnorm,ABS1(betai)*anorm,safmin)
         acoeff = scale*betai
         bcoeff = scale*alphai
         IF ( Left ) THEN
            acoeff = DCONJG(acoeff)
            bcoeff = DCONJG(bcoeff)
         ENDIF
         CALL ZGEMV(trans,N,N,acoeff,A,Lda,E(1,jvec),1,CZERO,           &
     &              Work(N*(jvec-1)+1),1)
         CALL ZGEMV(trans,N,N,-bcoeff,B,Lda,E(1,jvec),1,CONE,           &
     &              Work(N*(jvec-1)+1),1)
      ENDDO
!
      errnrm = ZLANGE('One',N,N,Work,N,Rwork)/enorm
!
!     Compute RESULT(1)
!
      Result(1) = errnrm/ulp
!
!     Normalization of E:
!
      enrmer = ZERO
      DO jvec = 1 , N
         temp1 = ZERO
         DO j = 1 , N
            temp1 = MAX(temp1,ABS1(E(j,jvec)))
         ENDDO
         enrmer = MAX(enrmer,temp1-ONE)
      ENDDO
!
!     Compute RESULT(2) : the normalization error in E.
!
      Result(2) = enrmer/(DBLE(N)*ulp)
!
!
!     End of zget52
!
      END SUBROUTINE zget52
