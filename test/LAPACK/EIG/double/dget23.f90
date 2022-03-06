!*==dget23.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DGET23
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET23( COMP, BALANC, JTYPE, THRESH, ISEED, NOUNIT, N,
!                          A, LDA, H, WR, WI, WR1, WI1, VL, LDVL, VR,
!                          LDVR, LRE, LDLRE, RCONDV, RCNDV1, RCDVIN,
!                          RCONDE, RCNDE1, RCDEIN, SCALE, SCALE1, RESULT,
!                          WORK, LWORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            COMP
!       CHARACTER          BALANC
!       INTEGER            INFO, JTYPE, LDA, LDLRE, LDVL, LDVR, LWORK, N,
!      $                   NOUNIT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 ), IWORK( * )
!       DOUBLE PRECISION   A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ),
!      $                   RCDEIN( * ), RCDVIN( * ), RCNDE1( * ),
!      $                   RCNDV1( * ), RCONDE( * ), RCONDV( * ),
!      $                   RESULT( 11 ), SCALE( * ), SCALE1( * ),
!      $                   VL( LDVL, * ), VR( LDVR, * ), WI( * ),
!      $                   WI1( * ), WORK( * ), WR( * ), WR1( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DGET23  checks the nonsymmetric eigenvalue problem driver SGEEVX.
!>    If COMP = .FALSE., the first 8 of the following tests will be
!>    performed on the input matrix A, and also test 9 if LWORK is
!>    sufficiently large.
!>    if COMP is .TRUE. all 11 tests will be performed.
!>
!>    (1)     | A * VR - VR * W | / ( n |A| ulp )
!>
!>      Here VR is the matrix of unit right eigenvectors.
!>      W is a block diagonal matrix, with a 1x1 block for each
!>      real eigenvalue and a 2x2 block for each complex conjugate
!>      pair.  If eigenvalues j and j+1 are a complex conjugate pair,
!>      so WR(j) = WR(j+1) = wr and WI(j) = - WI(j+1) = wi, then the
!>      2 x 2 block corresponding to the pair will be:
!>
!>              (  wr  wi  )
!>              ( -wi  wr  )
!>
!>      Such a block multiplying an n x 2 matrix  ( ur ui ) on the
!>      right will be the same as multiplying  ur + i*ui  by  wr + i*wi.
!>
!>    (2)     | A**H * VL - VL * W**H | / ( n |A| ulp )
!>
!>      Here VL is the matrix of unit left eigenvectors, A**H is the
!>      conjugate transpose of A, and W is as above.
!>
!>    (3)     | |VR(i)| - 1 | / ulp and largest component real
!>
!>      VR(i) denotes the i-th column of VR.
!>
!>    (4)     | |VL(i)| - 1 | / ulp and largest component real
!>
!>      VL(i) denotes the i-th column of VL.
!>
!>    (5)     0 if W(full) = W(partial), 1/ulp otherwise
!>
!>      W(full) denotes the eigenvalues computed when VR, VL, RCONDV
!>      and RCONDE are also computed, and W(partial) denotes the
!>      eigenvalues computed when only some of VR, VL, RCONDV, and
!>      RCONDE are computed.
!>
!>    (6)     0 if VR(full) = VR(partial), 1/ulp otherwise
!>
!>      VR(full) denotes the right eigenvectors computed when VL, RCONDV
!>      and RCONDE are computed, and VR(partial) denotes the result
!>      when only some of VL and RCONDV are computed.
!>
!>    (7)     0 if VL(full) = VL(partial), 1/ulp otherwise
!>
!>      VL(full) denotes the left eigenvectors computed when VR, RCONDV
!>      and RCONDE are computed, and VL(partial) denotes the result
!>      when only some of VR and RCONDV are computed.
!>
!>    (8)     0 if SCALE, ILO, IHI, ABNRM (full) =
!>                 SCALE, ILO, IHI, ABNRM (partial)
!>            1/ulp otherwise
!>
!>      SCALE, ILO, IHI and ABNRM describe how the matrix is balanced.
!>      (full) is when VR, VL, RCONDE and RCONDV are also computed, and
!>      (partial) is when some are not computed.
!>
!>    (9)     0 if RCONDV(full) = RCONDV(partial), 1/ulp otherwise
!>
!>      RCONDV(full) denotes the reciprocal condition numbers of the
!>      right eigenvectors computed when VR, VL and RCONDE are also
!>      computed. RCONDV(partial) denotes the reciprocal condition
!>      numbers when only some of VR, VL and RCONDE are computed.
!>
!>   (10)     |RCONDV - RCDVIN| / cond(RCONDV)
!>
!>      RCONDV is the reciprocal right eigenvector condition number
!>      computed by DGEEVX and RCDVIN (the precomputed true value)
!>      is supplied as input. cond(RCONDV) is the condition number of
!>      RCONDV, and takes errors in computing RCONDV into account, so
!>      that the resulting quantity should be O(ULP). cond(RCONDV) is
!>      essentially given by norm(A)/RCONDE.
!>
!>   (11)     |RCONDE - RCDEIN| / cond(RCONDE)
!>
!>      RCONDE is the reciprocal eigenvalue condition number
!>      computed by DGEEVX and RCDEIN (the precomputed true value)
!>      is supplied as input.  cond(RCONDE) is the condition number
!>      of RCONDE, and takes errors in computing RCONDE into account,
!>      so that the resulting quantity should be O(ULP). cond(RCONDE)
!>      is essentially given by norm(A)/RCONDV.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMP
!> \verbatim
!>          COMP is LOGICAL
!>          COMP describes which input tests to perform:
!>            = .FALSE. if the computed condition numbers are not to
!>                      be tested against RCDVIN and RCDEIN
!>            = .TRUE.  if they are to be compared
!> \endverbatim
!>
!> \param[in] BALANC
!> \verbatim
!>          BALANC is CHARACTER
!>          Describes the balancing option to be tested.
!>            = 'N' for no permuting or diagonal scaling
!>            = 'P' for permuting but no diagonal scaling
!>            = 'S' for no permuting but diagonal scaling
!>            = 'B' for permuting and diagonal scaling
!> \endverbatim
!>
!> \param[in] JTYPE
!> \verbatim
!>          JTYPE is INTEGER
!>          Type of input matrix. Used to label output if error occurs.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error
!>          is scaled to be O(1), so THRESH should be a reasonably
!>          small multiple of 1, e.g., 10 or 100.  In particular,
!>          it should not depend on the precision (single vs. double)
!>          or the size of the matrix.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          If COMP = .FALSE., the random number generator seed
!>          used to produce matrix.
!>          If COMP = .TRUE., ISEED(1) = the number of the example.
!>          Used to label output if error occurs.
!> \endverbatim
!>
!> \param[in] NOUNIT
!> \verbatim
!>          NOUNIT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns INFO not equal to 0.)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of A. N must be at least 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          Used to hold the matrix whose eigenvalues are to be
!>          computed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A, and H. LDA must be at
!>          least 1 and at least N.
!> \endverbatim
!>
!> \param[out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDA,N)
!>          Another copy of the test matrix A, modified by DGEEVX.
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (N)
!>
!>          The real and imaginary parts of the eigenvalues of A.
!>          On exit, WR + WI*i are the eigenvalues of the matrix in A.
!> \endverbatim
!>
!> \param[out] WR1
!> \verbatim
!>          WR1 is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI1
!> \verbatim
!>          WI1 is DOUBLE PRECISION array, dimension (N)
!>
!>          Like WR, WI, these arrays contain the eigenvalues of A,
!>          but those computed when DGEEVX only computes a partial
!>          eigendecomposition, i.e. not the eigenvalues and left
!>          and right eigenvectors.
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is DOUBLE PRECISION array, dimension (LDVL,N)
!>          VL holds the computed left eigenvectors.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          Leading dimension of VL. Must be at least max(1,N).
!> \endverbatim
!>
!> \param[out] VR
!> \verbatim
!>          VR is DOUBLE PRECISION array, dimension (LDVR,N)
!>          VR holds the computed right eigenvectors.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          Leading dimension of VR. Must be at least max(1,N).
!> \endverbatim
!>
!> \param[out] LRE
!> \verbatim
!>          LRE is DOUBLE PRECISION array, dimension (LDLRE,N)
!>          LRE holds the computed right or left eigenvectors.
!> \endverbatim
!>
!> \param[in] LDLRE
!> \verbatim
!>          LDLRE is INTEGER
!>          Leading dimension of LRE. Must be at least max(1,N).
!> \endverbatim
!>
!> \param[out] RCONDV
!> \verbatim
!>          RCONDV is DOUBLE PRECISION array, dimension (N)
!>          RCONDV holds the computed reciprocal condition numbers
!>          for eigenvectors.
!> \endverbatim
!>
!> \param[out] RCNDV1
!> \verbatim
!>          RCNDV1 is DOUBLE PRECISION array, dimension (N)
!>          RCNDV1 holds more computed reciprocal condition numbers
!>          for eigenvectors.
!> \endverbatim
!>
!> \param[in] RCDVIN
!> \verbatim
!>          RCDVIN is DOUBLE PRECISION array, dimension (N)
!>          When COMP = .TRUE. RCDVIN holds the precomputed reciprocal
!>          condition numbers for eigenvectors to be compared with
!>          RCONDV.
!> \endverbatim
!>
!> \param[out] RCONDE
!> \verbatim
!>          RCONDE is DOUBLE PRECISION array, dimension (N)
!>          RCONDE holds the computed reciprocal condition numbers
!>          for eigenvalues.
!> \endverbatim
!>
!> \param[out] RCNDE1
!> \verbatim
!>          RCNDE1 is DOUBLE PRECISION array, dimension (N)
!>          RCNDE1 holds more computed reciprocal condition numbers
!>          for eigenvalues.
!> \endverbatim
!>
!> \param[in] RCDEIN
!> \verbatim
!>          RCDEIN is DOUBLE PRECISION array, dimension (N)
!>          When COMP = .TRUE. RCDEIN holds the precomputed reciprocal
!>          condition numbers for eigenvalues to be compared with
!>          RCONDE.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION array, dimension (N)
!>          Holds information describing balancing of matrix.
!> \endverbatim
!>
!> \param[out] SCALE1
!> \verbatim
!>          SCALE1 is DOUBLE PRECISION array, dimension (N)
!>          Holds information describing balancing of matrix.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (11)
!>          The values computed by the 11 tests described above.
!>          The values are currently limited to 1/ulp, to avoid
!>          overflow.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          3*N, and 6*N+N**2 if tests 9, 10 or 11 are to be performed.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0,  successful exit.
!>          If <0, input parameter -INFO had an incorrect value.
!>          If >0, DGEEVX returned an error code, the absolute
!>                 value of which is returned.
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DGET23(Comp,Balanc,Jtype,Thresh,Iseed,Nounit,N,A,Lda,H,&
     &                  Wr,Wi,Wr1,Wi1,Vl,Ldvl,Vr,Ldvr,Lre,Ldlre,Rcondv, &
     &                  Rcndv1,Rcdvin,Rconde,Rcnde1,Rcdein,Scale,Scale1,&
     &                  Result,Work,Lwork,Iwork,Info)
      IMPLICIT NONE
!*--DGET23381
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Comp
      CHARACTER Balanc
      INTEGER Info , Jtype , Lda , Ldlre , Ldvl , Ldvr , Lwork , N ,    &
     &        Nounit
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4) , Iwork(*)
      DOUBLE PRECISION A(Lda,*) , H(Lda,*) , Lre(Ldlre,*) , Rcdein(*) , &
     &                 Rcdvin(*) , Rcnde1(*) , Rcndv1(*) , Rconde(*) ,  &
     &                 Rcondv(*) , Result(11) , Scale(*) , Scale1(*) ,  &
     &                 Vl(Ldvl,*) , Vr(Ldvr,*) , Wi(*) , Wi1(*) ,       &
     &                 Work(*) , Wr(*) , Wr1(*)
!     ..
!
!  =====================================================================
!
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      DOUBLE PRECISION EPSIN
      PARAMETER (EPSIN=5.9605D-8)
!     ..
!     .. Local Scalars ..
      LOGICAL balok , nobal
      CHARACTER sense
      INTEGER i , ihi , ihi1 , iinfo , ilo , ilo1 , isens , isensm , j ,&
     &        jj , kmin
      DOUBLE PRECISION abnrm , abnrm1 , eps , smlnum , tnrm , tol ,     &
     &                 tolin , ulp , ulpinv , v , vimin , vmax , vmx ,  &
     &                 vrmin , vrmx , vtst
!     ..
!     .. Local Arrays ..
      CHARACTER sens(2)
      DOUBLE PRECISION dum(1) , res(2)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , DLAPY2 , DNRM2
      EXTERNAL LSAME , DLAMCH , DLAPY2 , DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEEVX , DGET22 , DLACPY , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN
!     ..
!     .. Data statements ..
      DATA sens/'N' , 'V'/
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
      nobal = LSAME(Balanc,'N')
      balok = nobal .OR. LSAME(Balanc,'P') .OR. LSAME(Balanc,'S') .OR.  &
     &        LSAME(Balanc,'B')
      Info = 0
      IF ( .NOT.balok ) THEN
         Info = -2
      ELSEIF ( Thresh<ZERO ) THEN
         Info = -4
      ELSEIF ( Nounit<=0 ) THEN
         Info = -6
      ELSEIF ( N<0 ) THEN
         Info = -7
      ELSEIF ( Lda<1 .OR. Lda<N ) THEN
         Info = -9
      ELSEIF ( Ldvl<1 .OR. Ldvl<N ) THEN
         Info = -16
      ELSEIF ( Ldvr<1 .OR. Ldvr<N ) THEN
         Info = -18
      ELSEIF ( Ldlre<1 .OR. Ldlre<N ) THEN
         Info = -20
      ELSEIF ( Lwork<3*N .OR. (Comp .AND. Lwork<6*N+N*N) ) THEN
         Info = -31
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGET23',-Info)
         RETURN
      ENDIF
!
!     Quick return if nothing to do
!
      DO i = 1 , 11
         Result(i) = -ONE
      ENDDO
!
      IF ( N==0 ) RETURN
!
!     More Important constants
!
      ulp = DLAMCH('Precision')
      smlnum = DLAMCH('S')
      ulpinv = ONE/ulp
!
!     Compute eigenvalues and eigenvectors, and test them
!
      IF ( Lwork>=6*N+N*N ) THEN
         sense = 'B'
         isensm = 2
      ELSE
         sense = 'E'
         isensm = 1
      ENDIF
      CALL DLACPY('F',N,N,A,Lda,H,Lda)
      CALL DGEEVX(Balanc,'V','V',sense,N,H,Lda,Wr,Wi,Vl,Ldvl,Vr,Ldvr,   &
     &            ilo,ihi,Scale,abnrm,Rconde,Rcondv,Work,Lwork,Iwork,   &
     &            iinfo)
      IF ( iinfo/=0 ) THEN
         Result(1) = ulpinv
         IF ( Jtype/=22 ) THEN
            WRITE (Nounit,FMT=99002) 'DGEEVX1' , iinfo , N , Jtype ,    &
     &                               Balanc , Iseed
         ELSE
            WRITE (Nounit,FMT=99001) 'DGEEVX1' , iinfo , N , Iseed(1)
         ENDIF
         Info = ABS(iinfo)
         RETURN
      ENDIF
!
!     Do Test (1)
!
      CALL DGET22('N','N','N',N,A,Lda,Vr,Ldvr,Wr,Wi,Work,res)
      Result(1) = res(1)
!
!     Do Test (2)
!
      CALL DGET22('T','N','T',N,A,Lda,Vl,Ldvl,Wr,Wi,Work,res)
      Result(2) = res(1)
!
!     Do Test (3)
!
      DO j = 1 , N
         tnrm = ONE
         IF ( Wi(j)==ZERO ) THEN
            tnrm = DNRM2(N,Vr(1,j),1)
         ELSEIF ( Wi(j)>ZERO ) THEN
            tnrm = DLAPY2(DNRM2(N,Vr(1,j),1),DNRM2(N,Vr(1,j+1),1))
         ENDIF
         Result(3) = MAX(Result(3),MIN(ulpinv,ABS(tnrm-ONE)/ulp))
         IF ( Wi(j)>ZERO ) THEN
            vmx = ZERO
            vrmx = ZERO
            DO jj = 1 , N
               vtst = DLAPY2(Vr(jj,j),Vr(jj,j+1))
               IF ( vtst>vmx ) vmx = vtst
               IF ( Vr(jj,j+1)==ZERO .AND. ABS(Vr(jj,j))>vrmx )         &
     &              vrmx = ABS(Vr(jj,j))
            ENDDO
            IF ( vrmx/vmx<ONE-TWO*ulp ) Result(3) = ulpinv
         ENDIF
      ENDDO
!
!     Do Test (4)
!
      DO j = 1 , N
         tnrm = ONE
         IF ( Wi(j)==ZERO ) THEN
            tnrm = DNRM2(N,Vl(1,j),1)
         ELSEIF ( Wi(j)>ZERO ) THEN
            tnrm = DLAPY2(DNRM2(N,Vl(1,j),1),DNRM2(N,Vl(1,j+1),1))
         ENDIF
         Result(4) = MAX(Result(4),MIN(ulpinv,ABS(tnrm-ONE)/ulp))
         IF ( Wi(j)>ZERO ) THEN
            vmx = ZERO
            vrmx = ZERO
            DO jj = 1 , N
               vtst = DLAPY2(Vl(jj,j),Vl(jj,j+1))
               IF ( vtst>vmx ) vmx = vtst
               IF ( Vl(jj,j+1)==ZERO .AND. ABS(Vl(jj,j))>vrmx )         &
     &              vrmx = ABS(Vl(jj,j))
            ENDDO
            IF ( vrmx/vmx<ONE-TWO*ulp ) Result(4) = ulpinv
         ENDIF
      ENDDO
!
!     Test for all options of computing condition numbers
!
      DO isens = 1 , isensm
!
         sense = sens(isens)
!
!        Compute eigenvalues only, and test them
!
         CALL DLACPY('F',N,N,A,Lda,H,Lda)
         CALL DGEEVX(Balanc,'N','N',sense,N,H,Lda,Wr1,Wi1,dum,1,dum,1,  &
     &               ilo1,ihi1,Scale1,abnrm1,Rcnde1,Rcndv1,Work,Lwork,  &
     &               Iwork,iinfo)
         IF ( iinfo/=0 ) THEN
            Result(1) = ulpinv
            IF ( Jtype/=22 ) THEN
               WRITE (Nounit,FMT=99002) 'DGEEVX2' , iinfo , N , Jtype , &
     &                                  Balanc , Iseed
            ELSE
               WRITE (Nounit,FMT=99001) 'DGEEVX2' , iinfo , N , Iseed(1)
            ENDIF
            Info = ABS(iinfo)
            CYCLE
         ENDIF
!
!        Do Test (5)
!
         DO j = 1 , N
            IF ( Wr(j)/=Wr1(j) .OR. Wi(j)/=Wi1(j) ) Result(5) = ulpinv
         ENDDO
!
!        Do Test (8)
!
         IF ( .NOT.nobal ) THEN
            DO j = 1 , N
               IF ( Scale(j)/=Scale1(j) ) Result(8) = ulpinv
            ENDDO
            IF ( ilo/=ilo1 ) Result(8) = ulpinv
            IF ( ihi/=ihi1 ) Result(8) = ulpinv
            IF ( abnrm/=abnrm1 ) Result(8) = ulpinv
         ENDIF
!
!        Do Test (9)
!
         IF ( isens==2 .AND. N>1 ) THEN
            DO j = 1 , N
               IF ( Rcondv(j)/=Rcndv1(j) ) Result(9) = ulpinv
            ENDDO
         ENDIF
!
!        Compute eigenvalues and right eigenvectors, and test them
!
         CALL DLACPY('F',N,N,A,Lda,H,Lda)
         CALL DGEEVX(Balanc,'N','V',sense,N,H,Lda,Wr1,Wi1,dum,1,Lre,    &
     &               Ldlre,ilo1,ihi1,Scale1,abnrm1,Rcnde1,Rcndv1,Work,  &
     &               Lwork,Iwork,iinfo)
         IF ( iinfo/=0 ) THEN
            Result(1) = ulpinv
            IF ( Jtype/=22 ) THEN
               WRITE (Nounit,FMT=99002) 'DGEEVX3' , iinfo , N , Jtype , &
     &                                  Balanc , Iseed
            ELSE
               WRITE (Nounit,FMT=99001) 'DGEEVX3' , iinfo , N , Iseed(1)
            ENDIF
            Info = ABS(iinfo)
            CYCLE
         ENDIF
!
!        Do Test (5) again
!
         DO j = 1 , N
            IF ( Wr(j)/=Wr1(j) .OR. Wi(j)/=Wi1(j) ) Result(5) = ulpinv
         ENDDO
!
!        Do Test (6)
!
         DO j = 1 , N
            DO jj = 1 , N
               IF ( Vr(j,jj)/=Lre(j,jj) ) Result(6) = ulpinv
            ENDDO
         ENDDO
!
!        Do Test (8) again
!
         IF ( .NOT.nobal ) THEN
            DO j = 1 , N
               IF ( Scale(j)/=Scale1(j) ) Result(8) = ulpinv
            ENDDO
            IF ( ilo/=ilo1 ) Result(8) = ulpinv
            IF ( ihi/=ihi1 ) Result(8) = ulpinv
            IF ( abnrm/=abnrm1 ) Result(8) = ulpinv
         ENDIF
!
!        Do Test (9) again
!
         IF ( isens==2 .AND. N>1 ) THEN
            DO j = 1 , N
               IF ( Rcondv(j)/=Rcndv1(j) ) Result(9) = ulpinv
            ENDDO
         ENDIF
!
!        Compute eigenvalues and left eigenvectors, and test them
!
         CALL DLACPY('F',N,N,A,Lda,H,Lda)
         CALL DGEEVX(Balanc,'V','N',sense,N,H,Lda,Wr1,Wi1,Lre,Ldlre,dum,&
     &               1,ilo1,ihi1,Scale1,abnrm1,Rcnde1,Rcndv1,Work,Lwork,&
     &               Iwork,iinfo)
         IF ( iinfo/=0 ) THEN
            Result(1) = ulpinv
            IF ( Jtype/=22 ) THEN
               WRITE (Nounit,FMT=99002) 'DGEEVX4' , iinfo , N , Jtype , &
     &                                  Balanc , Iseed
            ELSE
               WRITE (Nounit,FMT=99001) 'DGEEVX4' , iinfo , N , Iseed(1)
            ENDIF
            Info = ABS(iinfo)
            CYCLE
         ENDIF
!
!        Do Test (5) again
!
         DO j = 1 , N
            IF ( Wr(j)/=Wr1(j) .OR. Wi(j)/=Wi1(j) ) Result(5) = ulpinv
         ENDDO
!
!        Do Test (7)
!
         DO j = 1 , N
            DO jj = 1 , N
               IF ( Vl(j,jj)/=Lre(j,jj) ) Result(7) = ulpinv
            ENDDO
         ENDDO
!
!        Do Test (8) again
!
         IF ( .NOT.nobal ) THEN
            DO j = 1 , N
               IF ( Scale(j)/=Scale1(j) ) Result(8) = ulpinv
            ENDDO
            IF ( ilo/=ilo1 ) Result(8) = ulpinv
            IF ( ihi/=ihi1 ) Result(8) = ulpinv
            IF ( abnrm/=abnrm1 ) Result(8) = ulpinv
         ENDIF
!
!        Do Test (9) again
!
         IF ( isens==2 .AND. N>1 ) THEN
            DO j = 1 , N
               IF ( Rcondv(j)/=Rcndv1(j) ) Result(9) = ulpinv
            ENDDO
         ENDIF
!
!
      ENDDO
!
!     If COMP, compare condition numbers to precomputed ones
!
      IF ( Comp ) THEN
         CALL DLACPY('F',N,N,A,Lda,H,Lda)
         CALL DGEEVX('N','V','V','B',N,H,Lda,Wr,Wi,Vl,Ldvl,Vr,Ldvr,ilo, &
     &               ihi,Scale,abnrm,Rconde,Rcondv,Work,Lwork,Iwork,    &
     &               iinfo)
         IF ( iinfo/=0 ) THEN
            Result(1) = ulpinv
            WRITE (Nounit,FMT=99001) 'DGEEVX5' , iinfo , N , Iseed(1)
            Info = ABS(iinfo)
            GOTO 99999
         ENDIF
!
!        Sort eigenvalues and condition numbers lexicographically
!        to compare with inputs
!
         DO i = 1 , N - 1
            kmin = i
            vrmin = Wr(i)
            vimin = Wi(i)
            DO j = i + 1 , N
               IF ( Wr(j)<vrmin ) THEN
                  kmin = j
                  vrmin = Wr(j)
                  vimin = Wi(j)
               ENDIF
            ENDDO
            Wr(kmin) = Wr(i)
            Wi(kmin) = Wi(i)
            Wr(i) = vrmin
            Wi(i) = vimin
            vrmin = Rconde(kmin)
            Rconde(kmin) = Rconde(i)
            Rconde(i) = vrmin
            vrmin = Rcondv(kmin)
            Rcondv(kmin) = Rcondv(i)
            Rcondv(i) = vrmin
         ENDDO
!
!        Compare condition numbers for eigenvectors
!        taking their condition numbers into account
!
         Result(10) = ZERO
         eps = MAX(EPSIN,ulp)
         v = MAX(DBLE(N)*eps*abnrm,smlnum)
         IF ( abnrm==ZERO ) v = ONE
         DO i = 1 , N
            IF ( v>Rcondv(i)*Rconde(i) ) THEN
               tol = Rcondv(i)
            ELSE
               tol = v/Rconde(i)
            ENDIF
            IF ( v>Rcdvin(i)*Rcdein(i) ) THEN
               tolin = Rcdvin(i)
            ELSE
               tolin = v/Rcdein(i)
            ENDIF
            tol = MAX(tol,smlnum/eps)
            tolin = MAX(tolin,smlnum/eps)
            IF ( eps*(Rcdvin(i)-tolin)>Rcondv(i)+tol ) THEN
               vmax = ONE/eps
            ELSEIF ( Rcdvin(i)-tolin>Rcondv(i)+tol ) THEN
               vmax = (Rcdvin(i)-tolin)/(Rcondv(i)+tol)
            ELSEIF ( Rcdvin(i)+tolin<eps*(Rcondv(i)-tol) ) THEN
               vmax = ONE/eps
            ELSEIF ( Rcdvin(i)+tolin<Rcondv(i)-tol ) THEN
               vmax = (Rcondv(i)-tol)/(Rcdvin(i)+tolin)
            ELSE
               vmax = ONE
            ENDIF
            Result(10) = MAX(Result(10),vmax)
         ENDDO
!
!        Compare condition numbers for eigenvalues
!        taking their condition numbers into account
!
         Result(11) = ZERO
         DO i = 1 , N
            IF ( v>Rcondv(i) ) THEN
               tol = ONE
            ELSE
               tol = v/Rcondv(i)
            ENDIF
            IF ( v>Rcdvin(i) ) THEN
               tolin = ONE
            ELSE
               tolin = v/Rcdvin(i)
            ENDIF
            tol = MAX(tol,smlnum/eps)
            tolin = MAX(tolin,smlnum/eps)
            IF ( eps*(Rcdein(i)-tolin)>Rconde(i)+tol ) THEN
               vmax = ONE/eps
            ELSEIF ( Rcdein(i)-tolin>Rconde(i)+tol ) THEN
               vmax = (Rcdein(i)-tolin)/(Rconde(i)+tol)
            ELSEIF ( Rcdein(i)+tolin<eps*(Rconde(i)-tol) ) THEN
               vmax = ONE/eps
            ELSEIF ( Rcdein(i)+tolin<Rconde(i)-tol ) THEN
               vmax = (Rconde(i)-tol)/(Rcdein(i)+tolin)
            ELSE
               vmax = ONE
            ENDIF
            Result(11) = MAX(Result(11),vmax)
         ENDDO
!
      ENDIF
!
99001 FORMAT (' DGET23: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', INPUT EXAMPLE NUMBER = ',I4)
99002 FORMAT (' DGET23: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,', BALANC = ',A,', ISEED=(',3(I5,','),I5,   &
     &        ')')
!
!
!     End of DGET23
!
99999 END SUBROUTINE DGET23
