!*==zgghrd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZGGHRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGGHRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgghrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgghrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgghrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,
!                          LDQ, Z, LDZ, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ, COMPZ
!       INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGGHRD reduces a pair of complex matrices (A,B) to generalized upper
!> Hessenberg form using unitary transformations, where A is a
!> general matrix and B is upper triangular.  The form of the
!> generalized eigenvalue problem is
!>    A*x = lambda*B*x,
!> and B is typically made upper triangular by computing its QR
!> factorization and moving the unitary matrix Q to the left side
!> of the equation.
!>
!> This subroutine simultaneously reduces A to a Hessenberg matrix H:
!>    Q**H*A*Z = H
!> and transforms B to another upper triangular matrix T:
!>    Q**H*B*Z = T
!> in order to reduce the problem to its standard form
!>    H*y = lambda*T*y
!> where y = Z**H*x.
!>
!> The unitary matrices Q and Z are determined as products of Givens
!> rotations.  They may either be formed explicitly, or they may be
!> postmultiplied into input matrices Q1 and Z1, so that
!>      Q1 * A * Z1**H = (Q1*Q) * H * (Z1*Z)**H
!>      Q1 * B * Z1**H = (Q1*Q) * T * (Z1*Z)**H
!> If Q1 is the unitary matrix from the QR factorization of B in the
!> original equation A*x = lambda*B*x, then ZGGHRD reduces the original
!> problem to generalized Hessenberg form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPQ
!> \verbatim
!>          COMPQ is CHARACTER*1
!>          = 'N': do not compute Q;
!>          = 'I': Q is initialized to the unit matrix, and the
!>                 unitary matrix Q is returned;
!>          = 'V': Q must contain a unitary matrix Q1 on entry,
!>                 and the product Q1*Q is returned.
!> \endverbatim
!>
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N': do not compute Z;
!>          = 'I': Z is initialized to the unit matrix, and the
!>                 unitary matrix Z is returned;
!>          = 'V': Z must contain a unitary matrix Z1 on entry,
!>                 and the product Z1*Z is returned.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          ILO and IHI mark the rows and columns of A which are to be
!>          reduced.  It is assumed that A is already upper triangular
!>          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are
!>          normally set by a previous call to ZGGBAL; otherwise they
!>          should be set to 1 and N respectively.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA, N)
!>          On entry, the N-by-N general matrix to be reduced.
!>          On exit, the upper triangle and the first subdiagonal of A
!>          are overwritten with the upper Hessenberg matrix H, and the
!>          rest is set to zero.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB, N)
!>          On entry, the N-by-N upper triangular matrix B.
!>          On exit, the upper triangular matrix T = Q**H B Z.  The
!>          elements below the diagonal are set to zero.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ, N)
!>          On entry, if COMPQ = 'V', the unitary matrix Q1, typically
!>          from the QR factorization of B.
!>          On exit, if COMPQ='I', the unitary matrix Q, and if
!>          COMPQ = 'V', the product Q1*Q.
!>          Not referenced if COMPQ='N'.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.
!>          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ, N)
!>          On entry, if COMPZ = 'V', the unitary matrix Z1.
!>          On exit, if COMPZ='I', the unitary matrix Z, and if
!>          COMPZ = 'V', the product Z1*Z.
!>          Not referenced if COMPZ='N'.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.
!>          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup complex16OTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  This routine reduces A to Hessenberg and B to triangular form by
!>  an unblocked reduction, as described in _Matrix_Computations_,
!>  by Golub and van Loan (Johns Hopkins Press).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGGHRD(Compq,Compz,N,Ilo,Ihi,A,Lda,B,Ldb,Q,Ldq,Z,Ldz,  &
     &                  Info)
      IMPLICIT NONE
!*--ZGGHRD208
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Compq , Compz
      INTEGER Ihi , Ilo , Info , Lda , Ldb , Ldq , Ldz , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , Q(Ldq,*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 CONE , CZERO
      PARAMETER (CONE=(1.0D+0,0.0D+0),CZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL ilq , ilz
      INTEGER icompq , icompz , jcol , jrow
      DOUBLE PRECISION c
      COMPLEX*16 ctemp , s
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZLARTG , ZLASET , ZROT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG , MAX
!     ..
!     .. Executable Statements ..
!
!     Decode COMPQ
!
      IF ( LSAME(Compq,'N') ) THEN
         ilq = .FALSE.
         icompq = 1
      ELSEIF ( LSAME(Compq,'V') ) THEN
         ilq = .TRUE.
         icompq = 2
      ELSEIF ( LSAME(Compq,'I') ) THEN
         ilq = .TRUE.
         icompq = 3
      ELSE
         icompq = 0
      ENDIF
!
!     Decode COMPZ
!
      IF ( LSAME(Compz,'N') ) THEN
         ilz = .FALSE.
         icompz = 1
      ELSEIF ( LSAME(Compz,'V') ) THEN
         ilz = .TRUE.
         icompz = 2
      ELSEIF ( LSAME(Compz,'I') ) THEN
         ilz = .TRUE.
         icompz = 3
      ELSE
         icompz = 0
      ENDIF
!
!     Test the input parameters.
!
      Info = 0
      IF ( icompq<=0 ) THEN
         Info = -1
      ELSEIF ( icompz<=0 ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Ilo<1 ) THEN
         Info = -4
      ELSEIF ( Ihi>N .OR. Ihi<Ilo-1 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -9
      ELSEIF ( (ilq .AND. Ldq<N) .OR. Ldq<1 ) THEN
         Info = -11
      ELSEIF ( (ilz .AND. Ldz<N) .OR. Ldz<1 ) THEN
         Info = -13
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGGHRD',-Info)
         RETURN
      ENDIF
!
!     Initialize Q and Z if desired.
!
      IF ( icompq==3 ) CALL ZLASET('Full',N,N,CZERO,CONE,Q,Ldq)
      IF ( icompz==3 ) CALL ZLASET('Full',N,N,CZERO,CONE,Z,Ldz)
!
!     Quick return if possible
!
      IF ( N<=1 ) RETURN
!
!     Zero out lower triangle of B
!
      DO jcol = 1 , N - 1
         DO jrow = jcol + 1 , N
            B(jrow,jcol) = CZERO
         ENDDO
      ENDDO
!
!     Reduce A and B
!
      DO jcol = Ilo , Ihi - 2
!
         DO jrow = Ihi , jcol + 2 , -1
!
!           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)
!
            ctemp = A(jrow-1,jcol)
            CALL ZLARTG(ctemp,A(jrow,jcol),c,s,A(jrow-1,jcol))
            A(jrow,jcol) = CZERO
            CALL ZROT(N-jcol,A(jrow-1,jcol+1),Lda,A(jrow,jcol+1),Lda,c, &
     &                s)
            CALL ZROT(N+2-jrow,B(jrow-1,jrow-1),Ldb,B(jrow,jrow-1),Ldb, &
     &                c,s)
            IF ( ilq ) CALL ZROT(N,Q(1,jrow-1),1,Q(1,jrow),1,c,DCONJG(s)&
     &                           )
!
!           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)
!
            ctemp = B(jrow,jrow)
            CALL ZLARTG(ctemp,B(jrow,jrow-1),c,s,B(jrow,jrow))
            B(jrow,jrow-1) = CZERO
            CALL ZROT(Ihi,A(1,jrow),1,A(1,jrow-1),1,c,s)
            CALL ZROT(jrow-1,B(1,jrow),1,B(1,jrow-1),1,c,s)
            IF ( ilz ) CALL ZROT(N,Z(1,jrow),1,Z(1,jrow-1),1,c,s)
         ENDDO
      ENDDO
!
!
!     End of ZGGHRD
!
      END SUBROUTINE ZGGHRD
