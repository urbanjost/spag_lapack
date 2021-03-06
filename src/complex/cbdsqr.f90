!*==cbdsqr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
 
!> \brief \b CBDSQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CBDSQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cbdsqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cbdsqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cbdsqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,
!                          LDU, C, LDC, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), RWORK( * )
!       COMPLEX            C( LDC, * ), U( LDU, * ), VT( LDVT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CBDSQR computes the singular values and, optionally, the right and/or
!> left singular vectors from the singular value decomposition (SVD) of
!> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
!> zero-shift QR algorithm.  The SVD of B has the form
!>
!>    B = Q * S * P**H
!>
!> where S is the diagonal matrix of singular values, Q is an orthogonal
!> matrix of left singular vectors, and P is an orthogonal matrix of
!> right singular vectors.  If left singular vectors are requested, this
!> subroutine actually returns U*Q instead of Q, and, if right singular
!> vectors are requested, this subroutine returns P**H*VT instead of
!> P**H, for given complex input matrices U and VT.  When U and VT are
!> the unitary matrices that reduce a general matrix A to bidiagonal
!> form: A = U*B*VT, as computed by CGEBRD, then
!>
!>    A = (U*Q) * S * (P**H*VT)
!>
!> is the SVD of A.  Optionally, the subroutine may also compute Q**H*C
!> for a given complex input matrix C.
!>
!> See "Computing  Small Singular Values of Bidiagonal Matrices With
!> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
!> LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
!> no. 5, pp. 873-912, Sept 1990) and
!> "Accurate singular values and differential qd algorithms," by
!> B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
!> Department, University of California at Berkeley, July 1992
!> for a detailed description of the algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  B is upper bidiagonal;
!>          = 'L':  B is lower bidiagonal.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B.  N >= 0.
!> \endverbatim
!>
!> \param[in] NCVT
!> \verbatim
!>          NCVT is INTEGER
!>          The number of columns of the matrix VT. NCVT >= 0.
!> \endverbatim
!>
!> \param[in] NRU
!> \verbatim
!>          NRU is INTEGER
!>          The number of rows of the matrix U. NRU >= 0.
!> \endverbatim
!>
!> \param[in] NCC
!> \verbatim
!>          NCC is INTEGER
!>          The number of columns of the matrix C. NCC >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          On entry, the n diagonal elements of the bidiagonal matrix B.
!>          On exit, if INFO=0, the singular values of B in decreasing
!>          order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          On entry, the N-1 offdiagonal elements of the bidiagonal
!>          matrix B.
!>          On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E
!>          will contain the diagonal and superdiagonal elements of a
!>          bidiagonal matrix orthogonally equivalent to the one given
!>          as input.
!> \endverbatim
!>
!> \param[in,out] VT
!> \verbatim
!>          VT is COMPLEX array, dimension (LDVT, NCVT)
!>          On entry, an N-by-NCVT matrix VT.
!>          On exit, VT is overwritten by P**H * VT.
!>          Not referenced if NCVT = 0.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.
!>          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.
!> \endverbatim
!>
!> \param[in,out] U
!> \verbatim
!>          U is COMPLEX array, dimension (LDU, N)
!>          On entry, an NRU-by-N matrix U.
!>          On exit, U is overwritten by U * Q.
!>          Not referenced if NRU = 0.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,NRU).
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC, NCC)
!>          On entry, an N-by-NCC matrix C.
!>          On exit, C is overwritten by Q**H * C.
!>          Not referenced if NCC = 0.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.
!>          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (4*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  If INFO = -i, the i-th argument had an illegal value
!>          > 0:  the algorithm did not converge; D and E contain the
!>                elements of a bidiagonal matrix which is orthogonally
!>                similar to the input matrix B;  if INFO = i, i
!>                elements of E have not converged to zero.
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  TOLMUL  REAL, default = max(10,min(100,EPS**(-1/8)))
!>          TOLMUL controls the convergence criterion of the QR loop.
!>          If it is positive, TOLMUL*EPS is the desired relative
!>             precision in the computed singular values.
!>          If it is negative, abs(TOLMUL*EPS*sigma_max) is the
!>             desired absolute accuracy in the computed singular
!>             values (corresponds to relative accuracy
!>             abs(TOLMUL*EPS) in the largest singular value.
!>          abs(TOLMUL) should be between 1 and 1/EPS, and preferably
!>             between 10 (for fast convergence) and .1/EPS
!>             (for there to be some accuracy in the results).
!>          Default is to lose at either one eighth or 2 of the
!>             available decimal digits in each computed singular value
!>             (whichever is smaller).
!>
!>  MAXITR  INTEGER, default = 6
!>          MAXITR controls the maximum number of passes of the
!>          algorithm through its inner loop. The algorithms stops
!>          (and so fails to converge) if the number of passes
!>          through the inner loop exceeds MAXITR*N**2.
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
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CBDSQR(Uplo,N,Ncvt,Nru,Ncc,D,E,Vt,Ldvt,U,Ldu,C,Ldc,    &
     &                  Rwork,Info)
      IMPLICIT NONE
!*--CBDSQR227
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Ldc , Ldu , Ldvt , N , Ncc , Ncvt , Nru
!     ..
!     .. Array Arguments ..
      REAL D(*) , E(*) , Rwork(*)
      COMPLEX C(Ldc,*) , U(Ldu,*) , Vt(Ldvt,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
      REAL ONE
      PARAMETER (ONE=1.0E0)
      REAL NEGONE
      PARAMETER (NEGONE=-1.0E0)
      REAL HNDRTH
      PARAMETER (HNDRTH=0.01E0)
      REAL TEN
      PARAMETER (TEN=10.0E0)
      REAL HNDRD
      PARAMETER (HNDRD=100.0E0)
      REAL MEIGTH
      PARAMETER (MEIGTH=-0.125E0)
      INTEGER MAXITR
      PARAMETER (MAXITR=6)
!     ..
!     .. Local Scalars ..
      LOGICAL lower , rotate
      INTEGER i , idir , isub , iter , j , ll , lll , m , maxit , nm1 , &
     &        nm12 , nm13 , oldll , oldm
      REAL abse , abss , cosl , cosr , cs , eps , f , g , h , mu ,      &
     &     oldcs , oldsn , r , shift , sigmn , sigmx , sinl , sinr ,    &
     &     sll , smax , smin , sminl , sminoa , sn , thresh , tol ,     &
     &     tolmul , unfl
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLAMCH
      EXTERNAL LSAME , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CLASR , CSROT , CSSCAL , CSWAP , SLARTG , SLAS2 ,        &
     &         SLASQ1 , SLASV2 , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , REAL , SIGN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      lower = LSAME(Uplo,'L')
      IF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.lower ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Ncvt<0 ) THEN
         Info = -3
      ELSEIF ( Nru<0 ) THEN
         Info = -4
      ELSEIF ( Ncc<0 ) THEN
         Info = -5
      ELSEIF ( (Ncvt==0 .AND. Ldvt<1) .OR. (Ncvt>0 .AND. Ldvt<MAX(1,N)) &
     &         ) THEN
         Info = -9
      ELSEIF ( Ldu<MAX(1,Nru) ) THEN
         Info = -11
      ELSEIF ( (Ncc==0 .AND. Ldc<1) .OR. (Ncc>0 .AND. Ldc<MAX(1,N)) )   &
     &         THEN
         Info = -13
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CBDSQR',-Info)
         RETURN
      ENDIF
      IF ( N==0 ) RETURN
      IF ( N==1 ) GOTO 300
!
!     ROTATE is true if any singular vectors desired, false otherwise
!
      rotate = (Ncvt>0) .OR. (Nru>0) .OR. (Ncc>0)
!
!     If no singular vectors desired, use qd algorithm
!
      IF ( .NOT.rotate ) THEN
         CALL SLASQ1(N,D,E,Rwork,Info)
!
!     If INFO equals 2, dqds didn't finish, try to finish
!
         IF ( Info/=2 ) RETURN
         Info = 0
      ENDIF
!
      nm1 = N - 1
      nm12 = nm1 + nm1
      nm13 = nm12 + nm1
      idir = 0
!
!     Get machine constants
!
      eps = SLAMCH('Epsilon')
      unfl = SLAMCH('Safe minimum')
!
!     If matrix lower bidiagonal, rotate to be upper bidiagonal
!     by applying Givens rotations on the left
!
      IF ( lower ) THEN
         DO i = 1 , N - 1
            CALL SLARTG(D(i),E(i),cs,sn,r)
            D(i) = r
            E(i) = sn*D(i+1)
            D(i+1) = cs*D(i+1)
            Rwork(i) = cs
            Rwork(nm1+i) = sn
         ENDDO
!
!        Update singular vectors if desired
!
         IF ( Nru>0 ) CALL CLASR('R','V','F',Nru,N,Rwork(1),Rwork(N),U, &
     &                           Ldu)
         IF ( Ncc>0 ) CALL CLASR('L','V','F',N,Ncc,Rwork(1),Rwork(N),C, &
     &                           Ldc)
      ENDIF
!
!     Compute singular values to relative accuracy TOL
!     (By setting TOL to be negative, algorithm will compute
!     singular values to absolute accuracy ABS(TOL)*norm(input matrix))
!
      tolmul = MAX(TEN,MIN(HNDRD,eps**MEIGTH))
      tol = tolmul*eps
!
!     Compute approximate maximum, minimum singular values
!
      smax = ZERO
      DO i = 1 , N
         smax = MAX(smax,ABS(D(i)))
      ENDDO
      DO i = 1 , N - 1
         smax = MAX(smax,ABS(E(i)))
      ENDDO
      sminl = ZERO
      IF ( tol>=ZERO ) THEN
!
!        Relative accuracy desired
!
         sminoa = ABS(D(1))
         IF ( sminoa/=ZERO ) THEN
            mu = sminoa
            DO i = 2 , N
               mu = ABS(D(i))*(mu/(mu+ABS(E(i-1))))
               sminoa = MIN(sminoa,mu)
               IF ( sminoa==ZERO ) EXIT
            ENDDO
         ENDIF
         sminoa = sminoa/SQRT(REAL(N))
         thresh = MAX(tol*sminoa,MAXITR*N*N*unfl)
      ELSE
!
!        Absolute accuracy desired
!
         thresh = MAX(ABS(tol)*smax,MAXITR*N*N*unfl)
      ENDIF
!
!     Prepare for main iteration loop for the singular values
!     (MAXIT is the maximum number of passes through the inner
!     loop permitted before nonconvergence signalled.)
!
      maxit = MAXITR*N*N
      iter = 0
      oldll = -1
      oldm = -1
!
!     M points to last element of unconverged part of matrix
!
      m = N
!
!     Begin main iteration loop
!
!
!     Check for convergence or exceeding iteration count
!
 100  DO WHILE ( m>1 )
         IF ( iter>maxit ) THEN
!
!     Maximum number of iterations exceeded, failure to converge
!
            Info = 0
            DO i = 1 , N - 1
               IF ( E(i)/=ZERO ) Info = Info + 1
            ENDDO
            GOTO 99999
         ELSE
!
!     Find diagonal block of matrix to work on
!
            IF ( tol<ZERO .AND. ABS(D(m))<=thresh ) D(m) = ZERO
            smax = ABS(D(m))
            smin = smax
            DO lll = 1 , m - 1
               ll = m - lll
               abss = ABS(D(ll))
               abse = ABS(E(ll))
               IF ( tol<ZERO .AND. abss<=thresh ) D(ll) = ZERO
               IF ( abse<=thresh ) GOTO 150
               smin = MIN(smin,abss)
               smax = MAX(smax,abss,abse)
            ENDDO
            ll = 0
            GOTO 200
         ENDIF
 150     E(ll) = ZERO
!
!     Matrix splits since E(LL) = 0
!
         IF ( ll==m-1 ) THEN
!
!        Convergence of bottom singular value, return to top of loop
!
            m = m - 1
            CYCLE
         ENDIF
         GOTO 200
      ENDDO
      GOTO 300
 200  ll = ll + 1
!
!     E(LL) through E(M-1) are nonzero, E(LL-1) is zero
!
      IF ( ll==m-1 ) THEN
!
!        2 by 2 block, handle separately
!
         CALL SLASV2(D(m-1),E(m-1),D(m),sigmn,sigmx,sinr,cosr,sinl,cosl)
         D(m-1) = sigmx
         E(m-1) = ZERO
         D(m) = sigmn
!
!        Compute singular vectors, if desired
!
         IF ( Ncvt>0 ) CALL CSROT(Ncvt,Vt(m-1,1),Ldvt,Vt(m,1),Ldvt,cosr,&
     &                            sinr)
         IF ( Nru>0 ) CALL CSROT(Nru,U(1,m-1),1,U(1,m),1,cosl,sinl)
         IF ( Ncc>0 ) CALL CSROT(Ncc,C(m-1,1),Ldc,C(m,1),Ldc,cosl,sinl)
         m = m - 2
         GOTO 100
      ENDIF
!
!     If working on new submatrix, choose shift direction
!     (from larger end diagonal element towards smaller)
!
      IF ( ll>oldm .OR. m<oldll ) THEN
         IF ( ABS(D(ll))>=ABS(D(m)) ) THEN
!
!           Chase bulge from top (big end) to bottom (small end)
!
            idir = 1
         ELSE
!
!           Chase bulge from bottom (big end) to top (small end)
!
            idir = 2
         ENDIF
      ENDIF
!
!     Apply convergence tests
!
      IF ( idir==1 ) THEN
!
!        Run convergence test in forward direction
!        First apply standard test to bottom of matrix
!
         IF ( ABS(E(m-1))<=ABS(tol)*ABS(D(m)) .OR.                      &
     &        (tol<ZERO .AND. ABS(E(m-1))<=thresh) ) THEN
            E(m-1) = ZERO
            GOTO 100
         ENDIF
!
         IF ( tol>=ZERO ) THEN
!
!           If relative accuracy desired,
!           apply convergence criterion forward
!
            mu = ABS(D(ll))
            sminl = mu
            DO lll = ll , m - 1
               IF ( ABS(E(lll))<=tol*mu ) THEN
                  E(lll) = ZERO
                  GOTO 100
               ENDIF
               mu = ABS(D(lll+1))*(mu/(mu+ABS(E(lll))))
               sminl = MIN(sminl,mu)
            ENDDO
         ENDIF
!
      ELSE
!
!        Run convergence test in backward direction
!        First apply standard test to top of matrix
!
         IF ( ABS(E(ll))<=ABS(tol)*ABS(D(ll)) .OR.                      &
     &        (tol<ZERO .AND. ABS(E(ll))<=thresh) ) THEN
            E(ll) = ZERO
            GOTO 100
         ENDIF
!
         IF ( tol>=ZERO ) THEN
!
!           If relative accuracy desired,
!           apply convergence criterion backward
!
            mu = ABS(D(m))
            sminl = mu
            DO lll = m - 1 , ll , -1
               IF ( ABS(E(lll))<=tol*mu ) THEN
                  E(lll) = ZERO
                  GOTO 100
               ENDIF
               mu = ABS(D(lll))*(mu/(mu+ABS(E(lll))))
               sminl = MIN(sminl,mu)
            ENDDO
         ENDIF
      ENDIF
      oldll = ll
      oldm = m
!
!     Compute shift.  First, test if shifting would ruin relative
!     accuracy, and if so set the shift to zero.
!
      IF ( tol>=ZERO .AND. N*tol*(sminl/smax)<=MAX(eps,HNDRTH*tol) )    &
     &     THEN
!
!        Use a zero shift to avoid loss of relative accuracy
!
         shift = ZERO
      ELSE
!
!        Compute the shift from 2-by-2 block at end of matrix
!
         IF ( idir==1 ) THEN
            sll = ABS(D(ll))
            CALL SLAS2(D(m-1),E(m-1),D(m),shift,r)
         ELSE
            sll = ABS(D(m))
            CALL SLAS2(D(ll),E(ll),D(ll+1),shift,r)
         ENDIF
!
!        Test if shift negligible, and if so set to zero
!
         IF ( sll>ZERO ) THEN
            IF ( (shift/sll)**2<eps ) shift = ZERO
         ENDIF
      ENDIF
!
!     Increment iteration count
!
      iter = iter + m - ll
!
!     If SHIFT = 0, do simplified QR iteration
!
      IF ( shift==ZERO ) THEN
         IF ( idir==1 ) THEN
!
!           Chase bulge from top to bottom
!           Save cosines and sines for later singular vector updates
!
            cs = ONE
            oldcs = ONE
            DO i = ll , m - 1
               CALL SLARTG(D(i)*cs,E(i),cs,sn,r)
               IF ( i>ll ) E(i-1) = oldsn*r
               CALL SLARTG(oldcs*r,D(i+1)*sn,oldcs,oldsn,D(i))
               Rwork(i-ll+1) = cs
               Rwork(i-ll+1+nm1) = sn
               Rwork(i-ll+1+nm12) = oldcs
               Rwork(i-ll+1+nm13) = oldsn
            ENDDO
            h = D(m)*cs
            D(m) = h*oldcs
            E(m-1) = h*oldsn
!
!           Update singular vectors
!
            IF ( Ncvt>0 ) CALL CLASR('L','V','F',m-ll+1,Ncvt,Rwork(1),  &
     &                               Rwork(N),Vt(ll,1),Ldvt)
            IF ( Nru>0 ) CALL CLASR('R','V','F',Nru,m-ll+1,Rwork(nm12+1)&
     &                              ,Rwork(nm13+1),U(1,ll),Ldu)
            IF ( Ncc>0 ) CALL CLASR('L','V','F',m-ll+1,Ncc,Rwork(nm12+1)&
     &                              ,Rwork(nm13+1),C(ll,1),Ldc)
!
!           Test convergence
!
            IF ( ABS(E(m-1))<=thresh ) E(m-1) = ZERO
!
         ELSE
!
!           Chase bulge from bottom to top
!           Save cosines and sines for later singular vector updates
!
            cs = ONE
            oldcs = ONE
            DO i = m , ll + 1 , -1
               CALL SLARTG(D(i)*cs,E(i-1),cs,sn,r)
               IF ( i<m ) E(i) = oldsn*r
               CALL SLARTG(oldcs*r,D(i-1)*sn,oldcs,oldsn,D(i))
               Rwork(i-ll) = cs
               Rwork(i-ll+nm1) = -sn
               Rwork(i-ll+nm12) = oldcs
               Rwork(i-ll+nm13) = -oldsn
            ENDDO
            h = D(ll)*cs
            D(ll) = h*oldcs
            E(ll) = h*oldsn
!
!           Update singular vectors
!
            IF ( Ncvt>0 ) CALL CLASR('L','V','B',m-ll+1,Ncvt,           &
     &                               Rwork(nm12+1),Rwork(nm13+1),       &
     &                               Vt(ll,1),Ldvt)
            IF ( Nru>0 ) CALL CLASR('R','V','B',Nru,m-ll+1,Rwork(1),    &
     &                              Rwork(N),U(1,ll),Ldu)
            IF ( Ncc>0 ) CALL CLASR('L','V','B',m-ll+1,Ncc,Rwork(1),    &
     &                              Rwork(N),C(ll,1),Ldc)
!
!           Test convergence
!
            IF ( ABS(E(ll))<=thresh ) E(ll) = ZERO
         ENDIF
!
!        Use nonzero shift
!
      ELSEIF ( idir==1 ) THEN
!
!           Chase bulge from top to bottom
!           Save cosines and sines for later singular vector updates
!
         f = (ABS(D(ll))-shift)*(SIGN(ONE,D(ll))+shift/D(ll))
         g = E(ll)
         DO i = ll , m - 1
            CALL SLARTG(f,g,cosr,sinr,r)
            IF ( i>ll ) E(i-1) = r
            f = cosr*D(i) + sinr*E(i)
            E(i) = cosr*E(i) - sinr*D(i)
            g = sinr*D(i+1)
            D(i+1) = cosr*D(i+1)
            CALL SLARTG(f,g,cosl,sinl,r)
            D(i) = r
            f = cosl*E(i) + sinl*D(i+1)
            D(i+1) = cosl*D(i+1) - sinl*E(i)
            IF ( i<m-1 ) THEN
               g = sinl*E(i+1)
               E(i+1) = cosl*E(i+1)
            ENDIF
            Rwork(i-ll+1) = cosr
            Rwork(i-ll+1+nm1) = sinr
            Rwork(i-ll+1+nm12) = cosl
            Rwork(i-ll+1+nm13) = sinl
         ENDDO
         E(m-1) = f
!
!           Update singular vectors
!
         IF ( Ncvt>0 ) CALL CLASR('L','V','F',m-ll+1,Ncvt,Rwork(1),     &
     &                            Rwork(N),Vt(ll,1),Ldvt)
         IF ( Nru>0 ) CALL CLASR('R','V','F',Nru,m-ll+1,Rwork(nm12+1),  &
     &                           Rwork(nm13+1),U(1,ll),Ldu)
         IF ( Ncc>0 ) CALL CLASR('L','V','F',m-ll+1,Ncc,Rwork(nm12+1),  &
     &                           Rwork(nm13+1),C(ll,1),Ldc)
!
!           Test convergence
!
         IF ( ABS(E(m-1))<=thresh ) E(m-1) = ZERO
!
      ELSE
!
!           Chase bulge from bottom to top
!           Save cosines and sines for later singular vector updates
!
         f = (ABS(D(m))-shift)*(SIGN(ONE,D(m))+shift/D(m))
         g = E(m-1)
         DO i = m , ll + 1 , -1
            CALL SLARTG(f,g,cosr,sinr,r)
            IF ( i<m ) E(i) = r
            f = cosr*D(i) + sinr*E(i-1)
            E(i-1) = cosr*E(i-1) - sinr*D(i)
            g = sinr*D(i-1)
            D(i-1) = cosr*D(i-1)
            CALL SLARTG(f,g,cosl,sinl,r)
            D(i) = r
            f = cosl*E(i-1) + sinl*D(i-1)
            D(i-1) = cosl*D(i-1) - sinl*E(i-1)
            IF ( i>ll+1 ) THEN
               g = sinl*E(i-2)
               E(i-2) = cosl*E(i-2)
            ENDIF
            Rwork(i-ll) = cosr
            Rwork(i-ll+nm1) = -sinr
            Rwork(i-ll+nm12) = cosl
            Rwork(i-ll+nm13) = -sinl
         ENDDO
         E(ll) = f
!
!           Test convergence
!
         IF ( ABS(E(ll))<=thresh ) E(ll) = ZERO
!
!           Update singular vectors if desired
!
         IF ( Ncvt>0 ) CALL CLASR('L','V','B',m-ll+1,Ncvt,Rwork(nm12+1),&
     &                            Rwork(nm13+1),Vt(ll,1),Ldvt)
         IF ( Nru>0 ) CALL CLASR('R','V','B',Nru,m-ll+1,Rwork(1),       &
     &                           Rwork(N),U(1,ll),Ldu)
         IF ( Ncc>0 ) CALL CLASR('L','V','B',m-ll+1,Ncc,Rwork(1),       &
     &                           Rwork(N),C(ll,1),Ldc)
      ENDIF
!
!     QR iteration finished, go back and check convergence
!
      GOTO 100
!
!     All singular values converged, so make them positive
!
 300  DO i = 1 , N
         IF ( D(i)<ZERO ) THEN
            D(i) = -D(i)
!
!           Change sign of singular vectors, if desired
!
            IF ( Ncvt>0 ) CALL CSSCAL(Ncvt,NEGONE,Vt(i,1),Ldvt)
         ENDIF
      ENDDO
!
!     Sort the singular values into decreasing order (insertion sort on
!     singular values, but only one transposition per singular vector)
!
      DO i = 1 , N - 1
!
!        Scan for smallest D(I)
!
         isub = 1
         smin = D(1)
         DO j = 2 , N + 1 - i
            IF ( D(j)<=smin ) THEN
               isub = j
               smin = D(j)
            ENDIF
         ENDDO
         IF ( isub/=N+1-i ) THEN
!
!           Swap singular values and vectors
!
            D(isub) = D(N+1-i)
            D(N+1-i) = smin
            IF ( Ncvt>0 ) CALL CSWAP(Ncvt,Vt(isub,1),Ldvt,Vt(N+1-i,1),  &
     &                               Ldvt)
            IF ( Nru>0 ) CALL CSWAP(Nru,U(1,isub),1,U(1,N+1-i),1)
            IF ( Ncc>0 ) CALL CSWAP(Ncc,C(isub,1),Ldc,C(N+1-i,1),Ldc)
         ENDIF
      ENDDO
!
!     End of CBDSQR
!
99999 END SUBROUTINE CBDSQR
