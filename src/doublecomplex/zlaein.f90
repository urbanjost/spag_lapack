!*==zlaein.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLAEIN computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse iteration.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAEIN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaein.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaein.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaein.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAEIN( RIGHTV, NOINIT, N, H, LDH, W, V, B, LDB, RWORK,
!                          EPS3, SMLNUM, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            NOINIT, RIGHTV
!       INTEGER            INFO, LDB, LDH, N
!       DOUBLE PRECISION   EPS3, SMLNUM
!       COMPLEX*16         W
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         B( LDB, * ), H( LDH, * ), V( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAEIN uses inverse iteration to find a right or left eigenvector
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
!>          H is COMPLEX*16 array, dimension (LDH,N)
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
!>          W is COMPLEX*16
!>          The eigenvalue of H whose corresponding right or left
!>          eigenvector is to be computed.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (N)
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
!>          B is COMPLEX*16 array, dimension (LDB,N)
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
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[in] EPS3
!> \verbatim
!>          EPS3 is DOUBLE PRECISION
!>          A small machine-dependent value which is used to perturb
!>          close eigenvalues, and to replace zero pivots.
!> \endverbatim
!>
!> \param[in] SMLNUM
!> \verbatim
!>          SMLNUM is DOUBLE PRECISION
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLAEIN(Rightv,Noinit,N,H,Ldh,W,V,B,Ldb,Rwork,Eps3,     &
     &                  Smlnum,Info)
      IMPLICIT NONE
!*--ZLAEIN153
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Noinit , Rightv
      INTEGER Info , Ldb , Ldh , N
      DOUBLE PRECISION Eps3 , Smlnum
      COMPLEX*16 W
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 B(Ldb,*) , H(Ldh,*) , V(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , TENTH
      PARAMETER (ONE=1.0D+0,TENTH=1.0D-1)
      COMPLEX*16 ZERO
      PARAMETER (ZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      CHARACTER normin , trans
      INTEGER i , ierr , its , j
      DOUBLE PRECISION growto , nrmsml , rootn , rtemp , scale , vnorm
      COMPLEX*16 cdum , ei , ej , temp , x
!     ..
!     .. External Functions ..
      INTEGER IZAMAX
      DOUBLE PRECISION DZASUM , DZNRM2
      COMPLEX*16 ZLADIV
      EXTERNAL IZAMAX , DZASUM , DZNRM2 , ZLADIV
!     ..
!     .. External Subroutines ..
      EXTERNAL ZDSCAL , ZLATRS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DIMAG , MAX , SQRT
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1(cdum) = ABS(DBLE(cdum)) + ABS(DIMAG(cdum))
!     ..
!     .. Executable Statements ..
!
      Info = 0
!
!     GROWTO is the threshold used in the acceptance test for an
!     eigenvector.
!
      rootn = SQRT(DBLE(N))
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
         vnorm = DZNRM2(N,V,1)
         CALL ZDSCAL(N,(Eps3*rootn)/MAX(vnorm,nrmsml),V,1)
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
               x = ZLADIV(B(i,i),ei)
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
               x = ZLADIV(ei,B(i,i))
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
               x = ZLADIV(B(j,j),ej)
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
               x = ZLADIV(ej,B(j,j))
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
         CALL ZLATRS('Upper',trans,'Nonunit',normin,N,B,Ldb,V,scale,    &
     &               Rwork,ierr)
         normin = 'Y'
!
!        Test for sufficient growth in the norm of v.
!
         vnorm = DZASUM(N,V,1)
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
 100  i = IZAMAX(N,V,1)
      CALL ZDSCAL(N,ONE/CABS1(V(i)),V,1)
!
!
!     End of ZLAEIN
!
      END SUBROUTINE ZLAEIN
