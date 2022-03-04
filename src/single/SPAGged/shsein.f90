!*==shsein.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SHSEIN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SHSEIN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/shsein.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/shsein.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/shsein.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, WR, WI,
!                          VL, LDVL, VR, LDVR, MM, M, WORK, IFAILL,
!                          IFAILR, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          EIGSRC, INITV, SIDE
!       INTEGER            INFO, LDH, LDVL, LDVR, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       INTEGER            IFAILL( * ), IFAILR( * )
!       REAL               H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WI( * ), WORK( * ), WR( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SHSEIN uses inverse iteration to find specified right and/or left
!> eigenvectors of a real upper Hessenberg matrix H.
!>
!> The right eigenvector x and the left eigenvector y of the matrix H
!> corresponding to an eigenvalue w are defined by:
!>
!>              H * x = w * x,     y**h * H = w * y**h
!>
!> where y**h denotes the conjugate transpose of the vector y.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'R': compute right eigenvectors only;
!>          = 'L': compute left eigenvectors only;
!>          = 'B': compute both right and left eigenvectors.
!> \endverbatim
!>
!> \param[in] EIGSRC
!> \verbatim
!>          EIGSRC is CHARACTER*1
!>          Specifies the source of eigenvalues supplied in (WR,WI):
!>          = 'Q': the eigenvalues were found using SHSEQR; thus, if
!>                 H has zero subdiagonal elements, and so is
!>                 block-triangular, then the j-th eigenvalue can be
!>                 assumed to be an eigenvalue of the block containing
!>                 the j-th row/column.  This property allows SHSEIN to
!>                 perform inverse iteration on just one diagonal block.
!>          = 'N': no assumptions are made on the correspondence
!>                 between eigenvalues and diagonal blocks.  In this
!>                 case, SHSEIN must always perform inverse iteration
!>                 using the whole matrix H.
!> \endverbatim
!>
!> \param[in] INITV
!> \verbatim
!>          INITV is CHARACTER*1
!>          = 'N': no initial vectors are supplied;
!>          = 'U': user-supplied initial vectors are stored in the arrays
!>                 VL and/or VR.
!> \endverbatim
!>
!> \param[in,out] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          Specifies the eigenvectors to be computed. To select the
!>          real eigenvector corresponding to a real eigenvalue WR(j),
!>          SELECT(j) must be set to .TRUE.. To select the complex
!>          eigenvector corresponding to a complex eigenvalue
!>          (WR(j),WI(j)), with complex conjugate (WR(j+1),WI(j+1)),
!>          either SELECT(j) or SELECT(j+1) or both must be set to
!>          .TRUE.; then on exit SELECT(j) is .TRUE. and SELECT(j+1) is
!>          .FALSE..
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
!>          H is REAL array, dimension (LDH,N)
!>          The upper Hessenberg matrix H.
!>          If a NaN is detected in H, the routine will return with INFO=-6.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H.  LDH >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] WR
!> \verbatim
!>          WR is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[in] WI
!> \verbatim
!>          WI is REAL array, dimension (N)
!>
!>          On entry, the real and imaginary parts of the eigenvalues of
!>          H; a complex conjugate pair of eigenvalues must be stored in
!>          consecutive elements of WR and WI.
!>          On exit, WR may have been altered since close eigenvalues
!>          are perturbed slightly in searching for independent
!>          eigenvectors.
!> \endverbatim
!>
!> \param[in,out] VL
!> \verbatim
!>          VL is REAL array, dimension (LDVL,MM)
!>          On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must
!>          contain starting vectors for the inverse iteration for the
!>          left eigenvectors; the starting vector for each eigenvector
!>          must be in the same column(s) in which the eigenvector will
!>          be stored.
!>          On exit, if SIDE = 'L' or 'B', the left eigenvectors
!>          specified by SELECT will be stored consecutively in the
!>          columns of VL, in the same order as their eigenvalues. A
!>          complex eigenvector corresponding to a complex eigenvalue is
!>          stored in two consecutive columns, the first holding the real
!>          part and the second the imaginary part.
!>          If SIDE = 'R', VL is not referenced.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL.
!>          LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is REAL array, dimension (LDVR,MM)
!>          On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must
!>          contain starting vectors for the inverse iteration for the
!>          right eigenvectors; the starting vector for each eigenvector
!>          must be in the same column(s) in which the eigenvector will
!>          be stored.
!>          On exit, if SIDE = 'R' or 'B', the right eigenvectors
!>          specified by SELECT will be stored consecutively in the
!>          columns of VR, in the same order as their eigenvalues. A
!>          complex eigenvector corresponding to a complex eigenvalue is
!>          stored in two consecutive columns, the first holding the real
!>          part and the second the imaginary part.
!>          If SIDE = 'L', VR is not referenced.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.
!>          LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.
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
!>          The number of columns in the arrays VL and/or VR required to
!>          store the eigenvectors; each selected real eigenvector
!>          occupies one column and each selected complex eigenvector
!>          occupies two columns.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension ((N+2)*N)
!> \endverbatim
!>
!> \param[out] IFAILL
!> \verbatim
!>          IFAILL is INTEGER array, dimension (MM)
!>          If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left
!>          eigenvector in the i-th column of VL (corresponding to the
!>          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the
!>          eigenvector converged satisfactorily. If the i-th and (i+1)th
!>          columns of VL hold a complex eigenvector, then IFAILL(i) and
!>          IFAILL(i+1) are set to the same value.
!>          If SIDE = 'R', IFAILL is not referenced.
!> \endverbatim
!>
!> \param[out] IFAILR
!> \verbatim
!>          IFAILR is INTEGER array, dimension (MM)
!>          If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right
!>          eigenvector in the i-th column of VR (corresponding to the
!>          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the
!>          eigenvector converged satisfactorily. If the i-th and (i+1)th
!>          columns of VR hold a complex eigenvector, then IFAILR(i) and
!>          IFAILR(i+1) are set to the same value.
!>          If SIDE = 'L', IFAILR is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, i is the number of eigenvectors which
!>                failed to converge; see IFAILL and IFAILR for further
!>                details.
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
!> \ingroup realOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Each eigenvector is normalized so that the element of largest
!>  magnitude has magnitude 1; here the magnitude of a complex number
!>  (x,y) is taken to be |x|+|y|.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SHSEIN(Side,Eigsrc,Initv,Select,N,H,Ldh,Wr,Wi,Vl,Ldvl, &
     &                  Vr,Ldvr,Mm,M,Work,Ifaill,Ifailr,Info)
      USE S_LSAME
      USE S_SISNAN
      USE S_SLAEIN
      USE S_SLAMCH
      USE S_SLANHS
      USE S_XERBLA
      IMPLICIT NONE
!*--SHSEIN272
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Eigsrc
      CHARACTER :: Initv
      LOGICAL , INTENT(INOUT) , DIMENSION(*) :: Select
      INTEGER , INTENT(IN) :: N
      REAL , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL , INTENT(INOUT) , DIMENSION(*) :: Wr
      REAL , INTENT(IN) , DIMENSION(*) :: Wi
      REAL , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      REAL , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifaill
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ifailr
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: bignum , eps3 , hnorm , smlnum , ulp , unfl , wki , wkr
      LOGICAL :: bothv , fromqr , leftv , noinit , pair , rightv
      INTEGER :: i , iinfo , k , kl , kln , kr , ksi , ksr , ldwork
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
!     .. Executable Statements ..
!
!     Decode and test the input parameters.
!
      bothv = LSAME(Side,'B')
      rightv = LSAME(Side,'R') .OR. bothv
      leftv = LSAME(Side,'L') .OR. bothv
!
      fromqr = LSAME(Eigsrc,'Q')
!
      noinit = LSAME(Initv,'N')
!
!     Set M to the number of columns required to store the selected
!     eigenvectors, and standardize the array SELECT.
!
      M = 0
      pair = .FALSE.
      DO k = 1 , N
         IF ( pair ) THEN
            pair = .FALSE.
            Select(k) = .FALSE.
         ELSEIF ( Wi(k)==ZERO ) THEN
            IF ( Select(k) ) M = M + 1
         ELSE
            pair = .TRUE.
            IF ( Select(k) .OR. Select(k+1) ) THEN
               Select(k) = .TRUE.
               M = M + 2
            ENDIF
         ENDIF
      ENDDO
!
      Info = 0
      IF ( .NOT.rightv .AND. .NOT.leftv ) THEN
         Info = -1
      ELSEIF ( .NOT.fromqr .AND. .NOT.LSAME(Eigsrc,'N') ) THEN
         Info = -2
      ELSEIF ( .NOT.noinit .AND. .NOT.LSAME(Initv,'U') ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( Ldh<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldvl<1 .OR. (leftv .AND. Ldvl<N) ) THEN
         Info = -11
      ELSEIF ( Ldvr<1 .OR. (rightv .AND. Ldvr<N) ) THEN
         Info = -13
      ELSEIF ( Mm<M ) THEN
         Info = -14
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SHSEIN',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( N==0 ) RETURN
!
!     Set machine-dependent constants.
!
      unfl = SLAMCH('Safe minimum')
      ulp = SLAMCH('Precision')
      smlnum = unfl*(N/ulp)
      bignum = (ONE-ulp)/smlnum
!
      ldwork = N + 1
!
      kl = 1
      kln = 0
      IF ( fromqr ) THEN
         kr = 0
      ELSE
         kr = N
      ENDIF
      ksr = 1
!
      DO k = 1 , N
         IF ( Select(k) ) THEN
!
!           Compute eigenvector(s) corresponding to W(K).
!
            IF ( fromqr ) THEN
!
!              If affiliation of eigenvalues is known, check whether
!              the matrix splits.
!
!              Determine KL and KR such that 1 <= KL <= K <= KR <= N
!              and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or
!              KR = N).
!
!              Then inverse iteration can be performed with the
!              submatrix H(KL:N,KL:N) for a left eigenvector, and with
!              the submatrix H(1:KR,1:KR) for a right eigenvector.
!
               DO i = k , kl + 1 , -1
                  IF ( H(i,i-1)==ZERO ) EXIT
               ENDDO
               kl = i
               IF ( k>kr ) THEN
                  DO i = k , N - 1
                     IF ( H(i+1,i)==ZERO ) EXIT
                  ENDDO
                  kr = i
               ENDIF
            ENDIF
!
            IF ( kl/=kln ) THEN
               kln = kl
!
!              Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it
!              has not ben computed before.
!
               hnorm = SLANHS('I',kr-kl+1,H(kl,kl),Ldh,Work)
               IF ( SISNAN(hnorm) ) THEN
                  Info = -6
                  RETURN
               ELSEIF ( hnorm>ZERO ) THEN
                  eps3 = hnorm*ulp
               ELSE
                  eps3 = smlnum
               ENDIF
            ENDIF
!
!           Perturb eigenvalue if it is close to any previous
!           selected eigenvalues affiliated to the submatrix
!           H(KL:KR,KL:KR). Close roots are modified by EPS3.
!
            wkr = Wr(k)
            wki = Wi(k)
            DO
               DO i = k - 1 , kl , -1
                  IF ( Select(i) .AND. ABS(Wr(i)-wkr)+ABS(Wi(i)-wki)    &
     &                 <eps3 ) THEN
                     wkr = wkr + eps3
                     GOTO 20
                  ENDIF
               ENDDO
               Wr(k) = wkr
!
               pair = wki/=ZERO
               IF ( pair ) THEN
                  ksi = ksr + 1
               ELSE
                  ksi = ksr
               ENDIF
               IF ( leftv ) THEN
!
!              Compute left eigenvector.
!
                  CALL SLAEIN(.FALSE.,noinit,N-kl+1,H(kl,kl),Ldh,wkr,   &
     &                        wki,Vl(kl,ksr),Vl(kl,ksi),Work,ldwork,    &
     &                        Work(N*N+N+1),eps3,smlnum,bignum,iinfo)
                  IF ( iinfo>0 ) THEN
                     IF ( pair ) THEN
                        Info = Info + 2
                     ELSE
                        Info = Info + 1
                     ENDIF
                     Ifaill(ksr) = k
                     Ifaill(ksi) = k
                  ELSE
                     Ifaill(ksr) = 0
                     Ifaill(ksi) = 0
                  ENDIF
                  DO i = 1 , kl - 1
                     Vl(i,ksr) = ZERO
                  ENDDO
                  IF ( pair ) THEN
                     DO i = 1 , kl - 1
                        Vl(i,ksi) = ZERO
                     ENDDO
                  ENDIF
               ENDIF
               IF ( rightv ) THEN
!
!              Compute right eigenvector.
!
                  CALL SLAEIN(.TRUE.,noinit,kr,H,Ldh,wkr,wki,Vr(1,ksr), &
     &                        Vr(1,ksi),Work,ldwork,Work(N*N+N+1),eps3, &
     &                        smlnum,bignum,iinfo)
                  IF ( iinfo>0 ) THEN
                     IF ( pair ) THEN
                        Info = Info + 2
                     ELSE
                        Info = Info + 1
                     ENDIF
                     Ifailr(ksr) = k
                     Ifailr(ksi) = k
                  ELSE
                     Ifailr(ksr) = 0
                     Ifailr(ksi) = 0
                  ENDIF
                  DO i = kr + 1 , N
                     Vr(i,ksr) = ZERO
                  ENDDO
                  IF ( pair ) THEN
                     DO i = kr + 1 , N
                        Vr(i,ksi) = ZERO
                     ENDDO
                  ENDIF
               ENDIF
!
               IF ( pair ) THEN
                  ksr = ksr + 2
               ELSE
                  ksr = ksr + 1
               ENDIF
               EXIT
 20         ENDDO
         ENDIF
      ENDDO
!
!
!     End of SHSEIN
!
      END SUBROUTINE SHSEIN
