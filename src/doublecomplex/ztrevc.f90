!*==ztrevc.f90  processed by SPAG 7.51RB at 20:09 on  3 Mar 2022
!> \brief \b ZTREVC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTREVC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrevc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrevc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrevc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
!                          LDVR, MM, M, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, SIDE
!       INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTREVC computes some or all of the right and/or left eigenvectors of
!> a complex upper triangular matrix T.
!> Matrices of this type are produced by the Schur factorization of
!> a complex general matrix:  A = Q*T*Q**H, as computed by ZHSEQR.
!>
!> The right eigenvector x and the left eigenvector y of T corresponding
!> to an eigenvalue w are defined by:
!>
!>              T*x = w*x,     (y**H)*T = w*(y**H)
!>
!> where y**H denotes the conjugate transpose of the vector y.
!> The eigenvalues are not input to this routine, but are read directly
!> from the diagonal of T.
!>
!> This routine returns the matrices X and/or Y of right and left
!> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
!> input matrix.  If Q is the unitary factor that reduces a matrix A to
!> Schur form T, then Q*X and Q*Y are the matrices of right and left
!> eigenvectors of A.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'R':  compute right eigenvectors only;
!>          = 'L':  compute left eigenvectors only;
!>          = 'B':  compute both right and left eigenvectors.
!> \endverbatim
!>
!> \param[in] HOWMNY
!> \verbatim
!>          HOWMNY is CHARACTER*1
!>          = 'A':  compute all right and/or left eigenvectors;
!>          = 'B':  compute all right and/or left eigenvectors,
!>                  backtransformed using the matrices supplied in
!>                  VR and/or VL;
!>          = 'S':  compute selected right and/or left eigenvectors,
!>                  as indicated by the logical array SELECT.
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
!>          computed.
!>          The eigenvector corresponding to the j-th eigenvalue is
!>          computed if SELECT(j) = .TRUE..
!>          Not referenced if HOWMNY = 'A' or 'B'.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,N)
!>          The upper triangular matrix T.  T is modified, but restored
!>          on exit.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] VL
!> \verbatim
!>          VL is COMPLEX*16 array, dimension (LDVL,MM)
!>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!>          contain an N-by-N matrix Q (usually the unitary matrix Q of
!>          Schur vectors returned by ZHSEQR).
!>          On exit, if SIDE = 'L' or 'B', VL contains:
!>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
!>          if HOWMNY = 'B', the matrix Q*Y;
!>          if HOWMNY = 'S', the left eigenvectors of T specified by
!>                           SELECT, stored consecutively in the columns
!>                           of VL, in the same order as their
!>                           eigenvalues.
!>          Not referenced if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL.  LDVL >= 1, and if
!>          SIDE = 'L' or 'B', LDVL >= N.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is COMPLEX*16 array, dimension (LDVR,MM)
!>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!>          contain an N-by-N matrix Q (usually the unitary matrix Q of
!>          Schur vectors returned by ZHSEQR).
!>          On exit, if SIDE = 'R' or 'B', VR contains:
!>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
!>          if HOWMNY = 'B', the matrix Q*X;
!>          if HOWMNY = 'S', the right eigenvectors of T specified by
!>                           SELECT, stored consecutively in the columns
!>                           of VR, in the same order as their
!>                           eigenvalues.
!>          Not referenced if SIDE = 'L'.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.  LDVR >= 1, and if
!>          SIDE = 'R' or 'B'; LDVR >= N.
!> \endverbatim
!>
!> \param[in] MM
!> \verbatim
!>          MM is INTEGER
!>          The number of columns in the arrays VL and/or VR. MM >= M.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The number of columns in the arrays VL and/or VR actually
!>          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M
!>          is set to N.  Each selected eigenvector occupies one
!>          column.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \date November 2017
!
!> \ingroup complex16OTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The algorithm used in this program is basically backward (forward)
!>  substitution, with scaling to make the the code robust against
!>  possible overflow.
!>
!>  Each eigenvector is normalized so that the element of largest
!>  magnitude has magnitude 1; here the magnitude of a complex number
!>  (x,y) is taken to be |x| + |y|.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZTREVC(Side,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,Mm,M,&
     &                  Work,Rwork,Info)
      IMPLICIT NONE
!*--ZTREVC222
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      CHARACTER Howmny , Side
      INTEGER Info , Ldt , Ldvl , Ldvr , M , Mm , N
!     ..
!     .. Array Arguments ..
      LOGICAL Select(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 T(Ldt,*) , Vl(Ldvl,*) , Vr(Ldvr,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CMZERO , CMONE
      PARAMETER (CMZERO=(0.0D+0,0.0D+0),CMONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL allv , bothv , leftv , over , rightv , somev
      INTEGER i , ii , is , j , k , ki
      DOUBLE PRECISION ovfl , remax , scale , smin , smlnum , ulp , unfl
      COMPLEX*16 cdum
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER IZAMAX
      DOUBLE PRECISION DLAMCH , DZASUM
      EXTERNAL LSAME , IZAMAX , DLAMCH , DZASUM
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZCOPY , ZDSCAL , ZGEMV , ZLATRS , DLABAD
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DCMPLX , DCONJG , DIMAG , MAX
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1(cdum) = ABS(DBLE(cdum)) + ABS(DIMAG(cdum))
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      bothv = LSAME(Side,'B')
      rightv = LSAME(Side,'R') .OR. bothv
      leftv = LSAME(Side,'L') .OR. bothv
!
      allv = LSAME(Howmny,'A')
      over = LSAME(Howmny,'B')
      somev = LSAME(Howmny,'S')
!
!     Set M to the number of columns required to store the selected
!     eigenvectors.
!
      IF ( somev ) THEN
         M = 0
         DO j = 1 , N
            IF ( Select(j) ) M = M + 1
         ENDDO
      ELSE
         M = N
      ENDIF
!
      Info = 0
      IF ( .NOT.rightv .AND. .NOT.leftv ) THEN
         Info = -1
      ELSEIF ( .NOT.allv .AND. .NOT.over .AND. .NOT.somev ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Ldt<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldvl<1 .OR. (leftv .AND. Ldvl<N) ) THEN
         Info = -8
      ELSEIF ( Ldvr<1 .OR. (rightv .AND. Ldvr<N) ) THEN
         Info = -10
      ELSEIF ( Mm<M ) THEN
         Info = -11
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZTREVC',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( N==0 ) RETURN
!
!     Set the constants to control overflow.
!
      unfl = DLAMCH('Safe minimum')
      ovfl = ONE/unfl
      CALL DLABAD(unfl,ovfl)
      ulp = DLAMCH('Precision')
      smlnum = unfl*(N/ulp)
!
!     Store the diagonal elements of T in working array WORK.
!
      DO i = 1 , N
         Work(i+N) = T(i,i)
      ENDDO
!
!     Compute 1-norm of each column of strictly upper triangular
!     part of T to control overflow in triangular solver.
!
      Rwork(1) = ZERO
      DO j = 2 , N
         Rwork(j) = DZASUM(j-1,T(1,j),1)
      ENDDO
!
      IF ( rightv ) THEN
!
!        Compute right eigenvectors.
!
         is = M
         DO ki = N , 1 , -1
!
            IF ( somev ) THEN
               IF ( .NOT.Select(ki) ) CYCLE
            ENDIF
            smin = MAX(ulp*(CABS1(T(ki,ki))),smlnum)
!
            Work(1) = CMONE
!
!           Form right-hand side.
!
            DO k = 1 , ki - 1
               Work(k) = -T(k,ki)
            ENDDO
!
!           Solve the triangular system:
!              (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.
!
            DO k = 1 , ki - 1
               T(k,k) = T(k,k) - T(ki,ki)
               IF ( CABS1(T(k,k))<smin ) T(k,k) = smin
            ENDDO
!
            IF ( ki>1 ) THEN
               CALL ZLATRS('Upper','No transpose','Non-unit','Y',ki-1,T,&
     &                     Ldt,Work(1),scale,Rwork,Info)
               Work(ki) = scale
            ENDIF
!
!           Copy the vector x or Q*x to VR and normalize.
!
            IF ( .NOT.over ) THEN
               CALL ZCOPY(ki,Work(1),1,Vr(1,is),1)
!
               ii = IZAMAX(ki,Vr(1,is),1)
               remax = ONE/CABS1(Vr(ii,is))
               CALL ZDSCAL(ki,remax,Vr(1,is),1)
!
               DO k = ki + 1 , N
                  Vr(k,is) = CMZERO
               ENDDO
            ELSE
               IF ( ki>1 ) CALL ZGEMV('N',N,ki-1,CMONE,Vr,Ldvr,Work(1), &
     &                                1,DCMPLX(scale),Vr(1,ki),1)
!
               ii = IZAMAX(N,Vr(1,ki),1)
               remax = ONE/CABS1(Vr(ii,ki))
               CALL ZDSCAL(N,remax,Vr(1,ki),1)
            ENDIF
!
!           Set back the original diagonal elements of T.
!
            DO k = 1 , ki - 1
               T(k,k) = Work(k+N)
            ENDDO
!
            is = is - 1
         ENDDO
      ENDIF
!
      IF ( leftv ) THEN
!
!        Compute left eigenvectors.
!
         is = 1
         DO ki = 1 , N
!
            IF ( somev ) THEN
               IF ( .NOT.Select(ki) ) CYCLE
            ENDIF
            smin = MAX(ulp*(CABS1(T(ki,ki))),smlnum)
!
            Work(N) = CMONE
!
!           Form right-hand side.
!
            DO k = ki + 1 , N
               Work(k) = -DCONJG(T(ki,k))
            ENDDO
!
!           Solve the triangular system:
!              (T(KI+1:N,KI+1:N) - T(KI,KI))**H * X = SCALE*WORK.
!
            DO k = ki + 1 , N
               T(k,k) = T(k,k) - T(ki,ki)
               IF ( CABS1(T(k,k))<smin ) T(k,k) = smin
            ENDDO
!
            IF ( ki<N ) THEN
               CALL ZLATRS('Upper','Conjugate transpose','Non-unit','Y',&
     &                     N-ki,T(ki+1,ki+1),Ldt,Work(ki+1),scale,Rwork,&
     &                     Info)
               Work(ki) = scale
            ENDIF
!
!           Copy the vector x or Q*x to VL and normalize.
!
            IF ( .NOT.over ) THEN
               CALL ZCOPY(N-ki+1,Work(ki),1,Vl(ki,is),1)
!
               ii = IZAMAX(N-ki+1,Vl(ki,is),1) + ki - 1
               remax = ONE/CABS1(Vl(ii,is))
               CALL ZDSCAL(N-ki+1,remax,Vl(ki,is),1)
!
               DO k = 1 , ki - 1
                  Vl(k,is) = CMZERO
               ENDDO
            ELSE
               IF ( ki<N ) CALL ZGEMV('N',N,N-ki,CMONE,Vl(1,ki+1),Ldvl, &
     &                                Work(ki+1),1,DCMPLX(scale),       &
     &                                Vl(1,ki),1)
!
               ii = IZAMAX(N,Vl(1,ki),1)
               remax = ONE/CABS1(Vl(ii,ki))
               CALL ZDSCAL(N,remax,Vl(1,ki),1)
            ENDIF
!
!           Set back the original diagonal elements of T.
!
            DO k = ki + 1 , N
               T(k,k) = Work(k+N)
            ENDDO
!
            is = is + 1
         ENDDO
      ENDIF
!
!
!     End of ZTREVC
!
      END SUBROUTINE ZTREVC
