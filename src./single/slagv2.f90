!*==slagv2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLAGV2 computes the Generalized Schur factorization of a real 2-by-2 matrix pencil (A,B) where B is upper triangular.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAGV2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slagv2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slagv2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slagv2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAGV2( A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, CSL, SNL,
!                          CSR, SNR )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB
!       REAL               CSL, CSR, SNL, SNR
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), ALPHAI( 2 ), ALPHAR( 2 ),
!      $                   B( LDB, * ), BETA( 2 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAGV2 computes the Generalized Schur factorization of a real 2-by-2
!> matrix pencil (A,B) where B is upper triangular. This routine
!> computes orthogonal (rotation) matrices given by CSL, SNL and CSR,
!> SNR such that
!>
!> 1) if the pencil (A,B) has two real eigenvalues (include 0/0 or 1/0
!>    types), then
!>
!>    [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
!>    [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
!>
!>    [ b11 b12 ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
!>    [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ],
!>
!> 2) if the pencil (A,B) has a pair of complex conjugate eigenvalues,
!>    then
!>
!>    [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
!>    [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
!>
!>    [ b11  0  ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
!>    [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ]
!>
!>    where b11 >= b22 > 0.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA, 2)
!>          On entry, the 2 x 2 matrix A.
!>          On exit, A is overwritten by the ``A-part'' of the
!>          generalized Schur form.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          THe leading dimension of the array A.  LDA >= 2.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB, 2)
!>          On entry, the upper triangular 2 x 2 matrix B.
!>          On exit, B is overwritten by the ``B-part'' of the
!>          generalized Schur form.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          THe leading dimension of the array B.  LDB >= 2.
!> \endverbatim
!>
!> \param[out] ALPHAR
!> \verbatim
!>          ALPHAR is REAL array, dimension (2)
!> \endverbatim
!>
!> \param[out] ALPHAI
!> \verbatim
!>          ALPHAI is REAL array, dimension (2)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is REAL array, dimension (2)
!>          (ALPHAR(k)+i*ALPHAI(k))/BETA(k) are the eigenvalues of the
!>          pencil (A,B), k=1,2, i = sqrt(-1).  Note that BETA(k) may
!>          be zero.
!> \endverbatim
!>
!> \param[out] CSL
!> \verbatim
!>          CSL is REAL
!>          The cosine of the left rotation matrix.
!> \endverbatim
!>
!> \param[out] SNL
!> \verbatim
!>          SNL is REAL
!>          The sine of the left rotation matrix.
!> \endverbatim
!>
!> \param[out] CSR
!> \verbatim
!>          CSR is REAL
!>          The cosine of the right rotation matrix.
!> \endverbatim
!>
!> \param[out] SNR
!> \verbatim
!>          SNR is REAL
!>          The sine of the right rotation matrix.
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
!> \ingroup realOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
!
!  =====================================================================
      SUBROUTINE SLAGV2(A,Lda,B,Ldb,Alphar,Alphai,Beta,Csl,Snl,Csr,Snr)
      USE S_SLAG2
      USE S_SLAMCH
      USE S_SLAPY2
      USE S_SLARTG
      USE S_SLASV2
      USE S_SROT
      IMPLICIT NONE
!*--SLAGV2166
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(2) :: Alphar
      REAL , INTENT(INOUT) , DIMENSION(2) :: Alphai
      REAL , INTENT(OUT) , DIMENSION(2) :: Beta
      REAL :: Csl
      REAL :: Snl
      REAL :: Csr
      REAL , INTENT(INOUT) :: Snr
!
! Local variable declarations rewritten by SPAG
!
      REAL :: anorm , ascale , bnorm , bscale , h1 , h2 , h3 , qq , r , &
     &        rr , safmin , scale1 , scale2 , t , ulp , wi , wr1 , wr2
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
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      safmin = SLAMCH('S')
      ulp = SLAMCH('P')
!
!     Scale A
!
      anorm = MAX(ABS(A(1,1))+ABS(A(2,1)),ABS(A(1,2))+ABS(A(2,2)),      &
     &        safmin)
      ascale = ONE/anorm
      A(1,1) = ascale*A(1,1)
      A(1,2) = ascale*A(1,2)
      A(2,1) = ascale*A(2,1)
      A(2,2) = ascale*A(2,2)
!
!     Scale B
!
      bnorm = MAX(ABS(B(1,1)),ABS(B(1,2))+ABS(B(2,2)),safmin)
      bscale = ONE/bnorm
      B(1,1) = bscale*B(1,1)
      B(1,2) = bscale*B(1,2)
      B(2,2) = bscale*B(2,2)
!
!     Check if A can be deflated
!
      IF ( ABS(A(2,1))<=ulp ) THEN
         Csl = ONE
         Snl = ZERO
         Csr = ONE
         Snr = ZERO
         A(2,1) = ZERO
         B(2,1) = ZERO
         wi = ZERO
!
!     Check if B is singular
!
      ELSEIF ( ABS(B(1,1))<=ulp ) THEN
         CALL SLARTG(A(1,1),A(2,1),Csl,Snl,r)
         Csr = ONE
         Snr = ZERO
         CALL SROT(2,A(1,1),Lda,A(2,1),Lda,Csl,Snl)
         CALL SROT(2,B(1,1),Ldb,B(2,1),Ldb,Csl,Snl)
         A(2,1) = ZERO
         B(1,1) = ZERO
         B(2,1) = ZERO
         wi = ZERO
!
      ELSEIF ( ABS(B(2,2))<=ulp ) THEN
         CALL SLARTG(A(2,2),A(2,1),Csr,Snr,t)
         Snr = -Snr
         CALL SROT(2,A(1,1),1,A(1,2),1,Csr,Snr)
         CALL SROT(2,B(1,1),1,B(1,2),1,Csr,Snr)
         Csl = ONE
         Snl = ZERO
         A(2,1) = ZERO
         B(2,1) = ZERO
         B(2,2) = ZERO
         wi = ZERO
!
      ELSE
!
!        B is nonsingular, first compute the eigenvalues of (A,B)
!
         CALL SLAG2(A,Lda,B,Ldb,safmin,scale1,scale2,wr1,wr2,wi)
!
         IF ( wi==ZERO ) THEN
!
!           two real eigenvalues, compute s*A-w*B
!
            h1 = scale1*A(1,1) - wr1*B(1,1)
            h2 = scale1*A(1,2) - wr1*B(1,2)
            h3 = scale1*A(2,2) - wr1*B(2,2)
!
            rr = SLAPY2(h1,h2)
            qq = SLAPY2(scale1*A(2,1),h3)
!
            IF ( rr>qq ) THEN
!
!              find right rotation matrix to zero 1,1 element of
!              (sA - wB)
!
               CALL SLARTG(h2,h1,Csr,Snr,t)
!
            ELSE
!
!              find right rotation matrix to zero 2,1 element of
!              (sA - wB)
!
               CALL SLARTG(h3,scale1*A(2,1),Csr,Snr,t)
!
            ENDIF
!
            Snr = -Snr
            CALL SROT(2,A(1,1),1,A(1,2),1,Csr,Snr)
            CALL SROT(2,B(1,1),1,B(1,2),1,Csr,Snr)
!
!           compute inf norms of A and B
!
            h1 = MAX(ABS(A(1,1))+ABS(A(1,2)),ABS(A(2,1))+ABS(A(2,2)))
            h2 = MAX(ABS(B(1,1))+ABS(B(1,2)),ABS(B(2,1))+ABS(B(2,2)))
!
            IF ( (scale1*h1)>=ABS(wr1)*h2 ) THEN
!
!              find left rotation matrix Q to zero out B(2,1)
!
               CALL SLARTG(B(1,1),B(2,1),Csl,Snl,r)
!
            ELSE
!
!              find left rotation matrix Q to zero out A(2,1)
!
               CALL SLARTG(A(1,1),A(2,1),Csl,Snl,r)
!
            ENDIF
!
            CALL SROT(2,A(1,1),Lda,A(2,1),Lda,Csl,Snl)
            CALL SROT(2,B(1,1),Ldb,B(2,1),Ldb,Csl,Snl)
!
            A(2,1) = ZERO
            B(2,1) = ZERO
!
         ELSE
!
!           a pair of complex conjugate eigenvalues
!           first compute the SVD of the matrix B
!
            CALL SLASV2(B(1,1),B(1,2),B(2,2),r,t,Snr,Csr,Snl,Csl)
!
!           Form (A,B) := Q(A,B)Z**T where Q is left rotation matrix and
!           Z is right rotation matrix computed from SLASV2
!
            CALL SROT(2,A(1,1),Lda,A(2,1),Lda,Csl,Snl)
            CALL SROT(2,B(1,1),Ldb,B(2,1),Ldb,Csl,Snl)
            CALL SROT(2,A(1,1),1,A(1,2),1,Csr,Snr)
            CALL SROT(2,B(1,1),1,B(1,2),1,Csr,Snr)
!
            B(2,1) = ZERO
            B(1,2) = ZERO
!
         ENDIF
!
      ENDIF
!
!     Unscaling
!
      A(1,1) = anorm*A(1,1)
      A(2,1) = anorm*A(2,1)
      A(1,2) = anorm*A(1,2)
      A(2,2) = anorm*A(2,2)
      B(1,1) = bnorm*B(1,1)
      B(2,1) = bnorm*B(2,1)
      B(1,2) = bnorm*B(1,2)
      B(2,2) = bnorm*B(2,2)
!
      IF ( wi==ZERO ) THEN
         Alphar(1) = A(1,1)
         Alphar(2) = A(2,2)
         Alphai(1) = ZERO
         Alphai(2) = ZERO
         Beta(1) = B(1,1)
         Beta(2) = B(2,2)
      ELSE
         Alphar(1) = anorm*wr1/scale1/bnorm
         Alphai(1) = anorm*wi/scale1/bnorm
         Alphar(2) = Alphar(1)
         Alphai(2) = -Alphai(1)
         Beta(1) = ONE
         Beta(2) = ONE
      ENDIF
!
!
!     End of SLAGV2
!
      END SUBROUTINE SLAGV2
