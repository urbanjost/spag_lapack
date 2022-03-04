!*==claein.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLAEIN computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse iteration.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAEIN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claein.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claein.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claein.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAEIN( RIGHTV, NOINIT, N, H, LDH, W, V, B, LDB, RWORK,
!                          EPS3, SMLNUM, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            NOINIT, RIGHTV
!       INTEGER            INFO, LDB, LDH, N
!       REAL               EPS3, SMLNUM
!       COMPLEX            W
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            B( LDB, * ), H( LDH, * ), V( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAEIN uses inverse iteration to find a right or left eigenvector
!> corresponding to the eigenvalue W of a complex upper Hessenberg
!> matrix H.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] RIGHTV
!> \verbatim
!>          RIGHTV is LOGICAL
!>          = .TRUE. : compute right eigenvector;
!>          = .FALSE.: compute left eigenvector.
!> \endverbatim
!>
!> \param[in] NOINIT
!> \verbatim
!>          NOINIT is LOGICAL
!>          = .TRUE. : no initial vector supplied in V
!>          = .FALSE.: initial vector supplied in V.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] H
!> \verbatim
!>          H is COMPLEX array, dimension (LDH,N)
!>          The upper Hessenberg matrix H.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H.  LDH >= max(1,N).
!> \endverbatim
!>
!> \param[in] W
!> \verbatim
!>          W is COMPLEX
!>          The eigenvalue of H whose corresponding right or left
!>          eigenvector is to be computed.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is COMPLEX array, dimension (N)
!>          On entry, if NOINIT = .FALSE., V must contain a starting
!>          vector for inverse iteration; otherwise V need not be set.
!>          On exit, V contains the computed eigenvector, normalized so
!>          that the component of largest magnitude has magnitude 1; here
!>          the magnitude of a complex number (x,y) is taken to be
!>          |x| + |y|.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[in] EPS3
!> \verbatim
!>          EPS3 is REAL
!>          A small machine-dependent value which is used to perturb
!>          close eigenvalues, and to replace zero pivots.
!> \endverbatim
!>
!> \param[in] SMLNUM
!> \verbatim
!>          SMLNUM is REAL
!>          A machine-dependent value close to the underflow threshold.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          = 1:  inverse iteration did not converge; V is set to the
!>                last iterate.
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLAEIN(Rightv,Noinit,N,H,Ldh,W,V,B,Ldb,Rwork,Eps3,     &
     &                  Smlnum,Info)
      USE S_CLADIV
      USE S_CLATRS
      USE S_CSSCAL
      USE S_ICAMAX
      USE S_SCASUM
      USE S_SCNRM2
      IMPLICIT NONE
!*--CLAEIN159
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , TENTH = 1.0E-1
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      LOGICAL , INTENT(IN) :: Rightv
      LOGICAL , INTENT(IN) :: Noinit
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Ldh,*) :: H
      INTEGER , INTENT(IN) :: Ldh
      COMPLEX , INTENT(IN) :: W
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: V
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Rwork
      REAL , INTENT(IN) :: Eps3
      REAL , INTENT(IN) :: Smlnum
      INTEGER , INTENT(OUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: CABS1
      COMPLEX :: cdum , ei , ej , temp , x
      REAL :: growto , nrmsml , rootn , rtemp , scale , vnorm
      INTEGER :: i , ierr , its , j
      CHARACTER :: normin , trans
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
!     .. Statement Functions ..
!     ..
!     .. Statement Function definitions ..
      CABS1(cdum) = ABS(REAL(cdum)) + ABS(AIMAG(cdum))
!     ..
!     .. Executable Statements ..
!
      Info = 0
!
!     GROWTO is the threshold used in the acceptance test for an
!     eigenvector.
!
      rootn = SQRT(REAL(N))
      growto = TENTH/rootn
      nrmsml = MAX(ONE,Eps3*rootn)*Smlnum
!
!     Form B = H - W*I (except that the subdiagonal elements are not
!     stored).
!
      DO j = 1 , N
         DO i = 1 , j - 1
            B(i,j) = H(i,j)
         ENDDO
         B(j,j) = H(j,j) - W
      ENDDO
!
      IF ( Noinit ) THEN
!
!        Initialize V.
!
         DO i = 1 , N
            V(i) = Eps3
         ENDDO
      ELSE
!
!        Scale supplied initial vector.
!
         vnorm = SCNRM2(N,V,1)
         CALL CSSCAL(N,(Eps3*rootn)/MAX(vnorm,nrmsml),V,1)
      ENDIF
!
      IF ( Rightv ) THEN
!
!        LU decomposition with partial pivoting of B, replacing zero
!        pivots by EPS3.
!
         DO i = 1 , N - 1
            ei = H(i+1,i)
            IF ( CABS1(B(i,i))<CABS1(ei) ) THEN
!
!              Interchange rows and eliminate.
!
               x = CLADIV(B(i,i),ei)
               B(i,i) = ei
               DO j = i + 1 , N
                  temp = B(i+1,j)
                  B(i+1,j) = B(i,j) - x*temp
                  B(i,j) = temp
               ENDDO
            ELSE
!
!              Eliminate without interchange.
!
               IF ( B(i,i)==ZERO ) B(i,i) = Eps3
               x = CLADIV(ei,B(i,i))
               IF ( x/=ZERO ) THEN
                  DO j = i + 1 , N
                     B(i+1,j) = B(i+1,j) - x*B(i,j)
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
         IF ( B(N,N)==ZERO ) B(N,N) = Eps3
!
         trans = 'N'
!
      ELSE
!
!        UL decomposition with partial pivoting of B, replacing zero
!        pivots by EPS3.
!
         DO j = N , 2 , -1
            ej = H(j,j-1)
            IF ( CABS1(B(j,j))<CABS1(ej) ) THEN
!
!              Interchange columns and eliminate.
!
               x = CLADIV(B(j,j),ej)
               B(j,j) = ej
               DO i = 1 , j - 1
                  temp = B(i,j-1)
                  B(i,j-1) = B(i,j) - x*temp
                  B(i,j) = temp
               ENDDO
            ELSE
!
!              Eliminate without interchange.
!
               IF ( B(j,j)==ZERO ) B(j,j) = Eps3
               x = CLADIV(ej,B(j,j))
               IF ( x/=ZERO ) THEN
                  DO i = 1 , j - 1
                     B(i,j-1) = B(i,j-1) - x*B(i,j)
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
         IF ( B(1,1)==ZERO ) B(1,1) = Eps3
!
         trans = 'C'
!
      ENDIF
!
      normin = 'N'
      DO its = 1 , N
!
!        Solve U*x = scale*v for a right eigenvector
!          or U**H *x = scale*v for a left eigenvector,
!        overwriting x on v.
!
         CALL CLATRS('Upper',trans,'Nonunit',normin,N,B,Ldb,V,scale,    &
     &               Rwork,ierr)
         normin = 'Y'
!
!        Test for sufficient growth in the norm of v.
!
         vnorm = SCASUM(N,V,1)
         IF ( vnorm>=growto*scale ) GOTO 100
!
!        Choose new orthogonal starting vector and try again.
!
         rtemp = Eps3/(rootn+ONE)
         V(1) = Eps3
         DO i = 2 , N
            V(i) = rtemp
         ENDDO
         V(N-its+1) = V(N-its+1) - Eps3*rootn
      ENDDO
!
!     Failure to find eigenvector in N iterations.
!
      Info = 1
!
!
!     Normalize eigenvector.
!
 100  i = ICAMAX(N,V,1)
      CALL CSSCAL(N,ONE/CABS1(V(i)),V,1)
!
!
!     End of CLAEIN
!
      END SUBROUTINE CLAEIN
