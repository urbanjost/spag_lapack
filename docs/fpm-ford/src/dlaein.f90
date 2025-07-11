!*==dlaein.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAEIN computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse iteration.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAEIN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaein.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaein.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaein.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAEIN( RIGHTV, NOINIT, N, H, LDH, WR, WI, VR, VI, B,
!                          LDB, WORK, EPS3, SMLNUM, BIGNUM, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            NOINIT, RIGHTV
!       INTEGER            INFO, LDB, LDH, N
!       DOUBLE PRECISION   BIGNUM, EPS3, SMLNUM, WI, WR
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   B( LDB, * ), H( LDH, * ), VI( * ), VR( * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAEIN uses inverse iteration to find a right or left eigenvector
!> corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg
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
!>          = .TRUE. : no initial vector supplied in (VR,VI).
!>          = .FALSE.: initial vector supplied in (VR,VI).
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
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>          The upper Hessenberg matrix H.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H.  LDH >= max(1,N).
!> \endverbatim
!>
!> \param[in] WR
!> \verbatim
!>          WR is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] WI
!> \verbatim
!>          WI is DOUBLE PRECISION
!>          The real and imaginary parts of the eigenvalue of H whose
!>          corresponding right or left eigenvector is to be computed.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[in,out] VI
!> \verbatim
!>          VI is DOUBLE PRECISION array, dimension (N)
!>          On entry, if NOINIT = .FALSE. and WI = 0.0, VR must contain
!>          a real starting vector for inverse iteration using the real
!>          eigenvalue WR; if NOINIT = .FALSE. and WI.ne.0.0, VR and VI
!>          must contain the real and imaginary parts of a complex
!>          starting vector for inverse iteration using the complex
!>          eigenvalue (WR,WI); otherwise VR and VI need not be set.
!>          On exit, if WI = 0.0 (real eigenvalue), VR contains the
!>          computed real eigenvector; if WI.ne.0.0 (complex eigenvalue),
!>          VR and VI contain the real and imaginary parts of the
!>          computed complex eigenvector. The eigenvector is normalized
!>          so that the component of largest magnitude has magnitude 1;
!>          here the magnitude of a complex number (x,y) is taken to be
!>          |x| + |y|.
!>          VI is not referenced if WI = 0.0.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= N+1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
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
!> \param[in] BIGNUM
!> \verbatim
!>          BIGNUM is DOUBLE PRECISION
!>          A machine-dependent value close to the overflow threshold.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          = 1:  inverse iteration did not converge; VR is set to the
!>                last iterate, and so is VI if WI.ne.0.0.
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
!  =====================================================================
      SUBROUTINE DLAEIN(Rightv,Noinit,N,H,Ldh,Wr,Wi,Vr,Vi,B,Ldb,Work,   &
     &                  Eps3,Smlnum,Bignum,Info)
      IMPLICIT NONE
!*--DLAEIN176
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Noinit , Rightv
      INTEGER Info , Ldb , Ldh , N
      DOUBLE PRECISION Bignum , Eps3 , Smlnum , Wi , Wr
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION B(Ldb,*) , H(Ldh,*) , Vi(*) , Vr(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TENTH
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0,TENTH=1.0D-1)
!     ..
!     .. Local Scalars ..
      CHARACTER normin , trans
      INTEGER i , i1 , i2 , i3 , ierr , its , j
      DOUBLE PRECISION absbii , absbjj , ei , ej , growto , norm ,      &
     &                 nrmsml , rec , rootn , scale , temp , vcrit ,    &
     &                 vmax , vnorm , w , w1 , x , xi , xr , y
!     ..
!     .. External Functions ..
      INTEGER IDAMAX
      DOUBLE PRECISION DASUM , DLAPY2 , DNRM2
      EXTERNAL IDAMAX , DASUM , DLAPY2 , DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL DLADIV , DLATRS , DSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , SQRT
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
!     Form B = H - (WR,WI)*I (except that the subdiagonal elements and
!     the imaginary parts of the diagonal elements are not stored).
!
      DO j = 1 , N
         DO i = 1 , j - 1
            B(i,j) = H(i,j)
         ENDDO
         B(j,j) = H(j,j) - Wr
      ENDDO
!
      IF ( Wi==ZERO ) THEN
!
!        Real eigenvalue.
!
         IF ( Noinit ) THEN
!
!           Set initial vector.
!
            DO i = 1 , N
               Vr(i) = Eps3
            ENDDO
         ELSE
!
!           Scale supplied initial vector.
!
            vnorm = DNRM2(N,Vr,1)
            CALL DSCAL(N,(Eps3*rootn)/MAX(vnorm,nrmsml),Vr,1)
         ENDIF
!
         IF ( Rightv ) THEN
!
!           LU decomposition with partial pivoting of B, replacing zero
!           pivots by EPS3.
!
            DO i = 1 , N - 1
               ei = H(i+1,i)
               IF ( ABS(B(i,i))<ABS(ei) ) THEN
!
!                 Interchange rows and eliminate.
!
                  x = B(i,i)/ei
                  B(i,i) = ei
                  DO j = i + 1 , N
                     temp = B(i+1,j)
                     B(i+1,j) = B(i,j) - x*temp
                     B(i,j) = temp
                  ENDDO
               ELSE
!
!                 Eliminate without interchange.
!
                  IF ( B(i,i)==ZERO ) B(i,i) = Eps3
                  x = ei/B(i,i)
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
!           UL decomposition with partial pivoting of B, replacing zero
!           pivots by EPS3.
!
            DO j = N , 2 , -1
               ej = H(j,j-1)
               IF ( ABS(B(j,j))<ABS(ej) ) THEN
!
!                 Interchange columns and eliminate.
!
                  x = B(j,j)/ej
                  B(j,j) = ej
                  DO i = 1 , j - 1
                     temp = B(i,j-1)
                     B(i,j-1) = B(i,j) - x*temp
                     B(i,j) = temp
                  ENDDO
               ELSE
!
!                 Eliminate without interchange.
!
                  IF ( B(j,j)==ZERO ) B(j,j) = Eps3
                  x = ej/B(j,j)
                  IF ( x/=ZERO ) THEN
                     DO i = 1 , j - 1
                        B(i,j-1) = B(i,j-1) - x*B(i,j)
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            IF ( B(1,1)==ZERO ) B(1,1) = Eps3
!
            trans = 'T'
!
         ENDIF
!
         normin = 'N'
         DO its = 1 , N
!
!           Solve U*x = scale*v for a right eigenvector
!             or U**T*x = scale*v for a left eigenvector,
!           overwriting x on v.
!
            CALL DLATRS('Upper',trans,'Nonunit',normin,N,B,Ldb,Vr,scale,&
     &                  Work,ierr)
            normin = 'Y'
!
!           Test for sufficient growth in the norm of v.
!
            vnorm = DASUM(N,Vr,1)
            IF ( vnorm>=growto*scale ) GOTO 50
!
!           Choose new orthogonal starting vector and try again.
!
            temp = Eps3/(rootn+ONE)
            Vr(1) = Eps3
            DO i = 2 , N
               Vr(i) = temp
            ENDDO
            Vr(N-its+1) = Vr(N-its+1) - Eps3*rootn
         ENDDO
!
!        Failure to find eigenvector in N iterations.
!
         Info = 1
!
!
!        Normalize eigenvector.
!
 50      i = IDAMAX(N,Vr,1)
         CALL DSCAL(N,ONE/ABS(Vr(i)),Vr,1)
      ELSE
!
!        Complex eigenvalue.
!
         IF ( Noinit ) THEN
!
!           Set initial vector.
!
            DO i = 1 , N
               Vr(i) = Eps3
               Vi(i) = ZERO
            ENDDO
         ELSE
!
!           Scale supplied initial vector.
!
            norm = DLAPY2(DNRM2(N,Vr,1),DNRM2(N,Vi,1))
            rec = (Eps3*rootn)/MAX(norm,nrmsml)
            CALL DSCAL(N,rec,Vr,1)
            CALL DSCAL(N,rec,Vi,1)
         ENDIF
!
         IF ( Rightv ) THEN
!
!           LU decomposition with partial pivoting of B, replacing zero
!           pivots by EPS3.
!
!           The imaginary part of the (i,j)-th element of U is stored in
!           B(j+1,i).
!
            B(2,1) = -Wi
            DO i = 2 , N
               B(i+1,1) = ZERO
            ENDDO
!
            DO i = 1 , N - 1
               absbii = DLAPY2(B(i,i),B(i+1,i))
               ei = H(i+1,i)
               IF ( absbii<ABS(ei) ) THEN
!
!                 Interchange rows and eliminate.
!
                  xr = B(i,i)/ei
                  xi = B(i+1,i)/ei
                  B(i,i) = ei
                  B(i+1,i) = ZERO
                  DO j = i + 1 , N
                     temp = B(i+1,j)
                     B(i+1,j) = B(i,j) - xr*temp
                     B(j+1,i+1) = B(j+1,i) - xi*temp
                     B(i,j) = temp
                     B(j+1,i) = ZERO
                  ENDDO
                  B(i+2,i) = -Wi
                  B(i+1,i+1) = B(i+1,i+1) - xi*Wi
                  B(i+2,i+1) = B(i+2,i+1) + xr*Wi
               ELSE
!
!                 Eliminate without interchanging rows.
!
                  IF ( absbii==ZERO ) THEN
                     B(i,i) = Eps3
                     B(i+1,i) = ZERO
                     absbii = Eps3
                  ENDIF
                  ei = (ei/absbii)/absbii
                  xr = B(i,i)*ei
                  xi = -B(i+1,i)*ei
                  DO j = i + 1 , N
                     B(i+1,j) = B(i+1,j) - xr*B(i,j) + xi*B(j+1,i)
                     B(j+1,i+1) = -xr*B(j+1,i) - xi*B(i,j)
                  ENDDO
                  B(i+2,i+1) = B(i+2,i+1) - Wi
               ENDIF
!
!              Compute 1-norm of offdiagonal elements of i-th row.
!
               Work(i) = DASUM(N-i,B(i,i+1),Ldb) + DASUM(N-i,B(i+2,i),1)
            ENDDO
            IF ( B(N,N)==ZERO .AND. B(N+1,N)==ZERO ) B(N,N) = Eps3
            Work(N) = ZERO
!
            i1 = N
            i2 = 1
            i3 = -1
         ELSE
!
!           UL decomposition with partial pivoting of conjg(B),
!           replacing zero pivots by EPS3.
!
!           The imaginary part of the (i,j)-th element of U is stored in
!           B(j+1,i).
!
            B(N+1,N) = Wi
            DO j = 1 , N - 1
               B(N+1,j) = ZERO
            ENDDO
!
            DO j = N , 2 , -1
               ej = H(j,j-1)
               absbjj = DLAPY2(B(j,j),B(j+1,j))
               IF ( absbjj<ABS(ej) ) THEN
!
!                 Interchange columns and eliminate
!
                  xr = B(j,j)/ej
                  xi = B(j+1,j)/ej
                  B(j,j) = ej
                  B(j+1,j) = ZERO
                  DO i = 1 , j - 1
                     temp = B(i,j-1)
                     B(i,j-1) = B(i,j) - xr*temp
                     B(j,i) = B(j+1,i) - xi*temp
                     B(i,j) = temp
                     B(j+1,i) = ZERO
                  ENDDO
                  B(j+1,j-1) = Wi
                  B(j-1,j-1) = B(j-1,j-1) + xi*Wi
                  B(j,j-1) = B(j,j-1) - xr*Wi
               ELSE
!
!                 Eliminate without interchange.
!
                  IF ( absbjj==ZERO ) THEN
                     B(j,j) = Eps3
                     B(j+1,j) = ZERO
                     absbjj = Eps3
                  ENDIF
                  ej = (ej/absbjj)/absbjj
                  xr = B(j,j)*ej
                  xi = -B(j+1,j)*ej
                  DO i = 1 , j - 1
                     B(i,j-1) = B(i,j-1) - xr*B(i,j) + xi*B(j+1,i)
                     B(j,i) = -xr*B(j+1,i) - xi*B(i,j)
                  ENDDO
                  B(j,j-1) = B(j,j-1) + Wi
               ENDIF
!
!              Compute 1-norm of offdiagonal elements of j-th column.
!
               Work(j) = DASUM(j-1,B(1,j),1) + DASUM(j-1,B(j+1,1),Ldb)
            ENDDO
            IF ( B(1,1)==ZERO .AND. B(2,1)==ZERO ) B(1,1) = Eps3
            Work(1) = ZERO
!
            i1 = 1
            i2 = N
            i3 = 1
         ENDIF
!
         DO its = 1 , N
            scale = ONE
            vmax = ONE
            vcrit = Bignum
!
!           Solve U*(xr,xi) = scale*(vr,vi) for a right eigenvector,
!             or U**T*(xr,xi) = scale*(vr,vi) for a left eigenvector,
!           overwriting (xr,xi) on (vr,vi).
!
            DO i = i1 , i2 , i3
!
               IF ( Work(i)>vcrit ) THEN
                  rec = ONE/vmax
                  CALL DSCAL(N,rec,Vr,1)
                  CALL DSCAL(N,rec,Vi,1)
                  scale = scale*rec
                  vmax = ONE
                  vcrit = Bignum
               ENDIF
!
               xr = Vr(i)
               xi = Vi(i)
               IF ( Rightv ) THEN
                  DO j = i + 1 , N
                     xr = xr - B(i,j)*Vr(j) + B(j+1,i)*Vi(j)
                     xi = xi - B(i,j)*Vi(j) - B(j+1,i)*Vr(j)
                  ENDDO
               ELSE
                  DO j = 1 , i - 1
                     xr = xr - B(j,i)*Vr(j) + B(i+1,j)*Vi(j)
                     xi = xi - B(j,i)*Vi(j) - B(i+1,j)*Vr(j)
                  ENDDO
               ENDIF
!
               w = ABS(B(i,i)) + ABS(B(i+1,i))
               IF ( w>Smlnum ) THEN
                  IF ( w<ONE ) THEN
                     w1 = ABS(xr) + ABS(xi)
                     IF ( w1>w*Bignum ) THEN
                        rec = ONE/w1
                        CALL DSCAL(N,rec,Vr,1)
                        CALL DSCAL(N,rec,Vi,1)
                        xr = Vr(i)
                        xi = Vi(i)
                        scale = scale*rec
                        vmax = vmax*rec
                     ENDIF
                  ENDIF
!
!                 Divide by diagonal element of B.
!
                  CALL DLADIV(xr,xi,B(i,i),B(i+1,i),Vr(i),Vi(i))
                  vmax = MAX(ABS(Vr(i))+ABS(Vi(i)),vmax)
                  vcrit = Bignum/vmax
               ELSE
                  DO j = 1 , N
                     Vr(j) = ZERO
                     Vi(j) = ZERO
                  ENDDO
                  Vr(i) = ONE
                  Vi(i) = ONE
                  scale = ZERO
                  vmax = ONE
                  vcrit = Bignum
               ENDIF
            ENDDO
!
!           Test for sufficient growth in the norm of (VR,VI).
!
            vnorm = DASUM(N,Vr,1) + DASUM(N,Vi,1)
            IF ( vnorm>=growto*scale ) GOTO 100
!
!           Choose a new orthogonal starting vector and try again.
!
            y = Eps3/(rootn+ONE)
            Vr(1) = Eps3
            Vi(1) = ZERO
!
            DO i = 2 , N
               Vr(i) = y
               Vi(i) = ZERO
            ENDDO
            Vr(N-its+1) = Vr(N-its+1) - Eps3*rootn
         ENDDO
!
!        Failure to find eigenvector in N iterations
!
         Info = 1
!
!
!        Normalize eigenvector.
!
 100     vnorm = ZERO
         DO i = 1 , N
            vnorm = MAX(vnorm,ABS(Vr(i))+ABS(Vi(i)))
         ENDDO
         CALL DSCAL(N,ONE/vnorm,Vr,1)
         CALL DSCAL(N,ONE/vnorm,Vi,1)
!
      ENDIF
!
!
!     End of DLAEIN
!
      END SUBROUTINE DLAEIN
