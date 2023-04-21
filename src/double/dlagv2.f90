!*==dlagv2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAGV2 computes the Generalized Schur factorization of a real 2-by-2 matrix pencil (A,B) where B is upper triangular.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAGV2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlagv2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlagv2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlagv2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAGV2( A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, CSL, SNL,
!                          CSR, SNR )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB
!       DOUBLE PRECISION   CSL, CSR, SNL, SNR
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), ALPHAI( 2 ), ALPHAR( 2 ),
!      $                   B( LDB, * ), BETA( 2 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAGV2 computes the Generalized Schur factorization of a real 2-by-2
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
!>          A is DOUBLE PRECISION array, dimension (LDA, 2)
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
!>          B is DOUBLE PRECISION array, dimension (LDB, 2)
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
!>          ALPHAR is DOUBLE PRECISION array, dimension (2)
!> \endverbatim
!>
!> \param[out] ALPHAI
!> \verbatim
!>          ALPHAI is DOUBLE PRECISION array, dimension (2)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION array, dimension (2)
!>          (ALPHAR(k)+i*ALPHAI(k))/BETA(k) are the eigenvalues of the
!>          pencil (A,B), k=1,2, i = sqrt(-1).  Note that BETA(k) may
!>          be zero.
!> \endverbatim
!>
!> \param[out] CSL
!> \verbatim
!>          CSL is DOUBLE PRECISION
!>          The cosine of the left rotation matrix.
!> \endverbatim
!>
!> \param[out] SNL
!> \verbatim
!>          SNL is DOUBLE PRECISION
!>          The sine of the left rotation matrix.
!> \endverbatim
!>
!> \param[out] CSR
!> \verbatim
!>          CSR is DOUBLE PRECISION
!>          The cosine of the right rotation matrix.
!> \endverbatim
!>
!> \param[out] SNR
!> \verbatim
!>          SNR is DOUBLE PRECISION
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
!> \ingroup doubleOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
!
!  =====================================================================
      SUBROUTINE DLAGV2(A,Lda,B,Ldb,Alphar,Alphai,Beta,Csl,Snl,Csr,Snr)
      IMPLICIT NONE
!*--DLAGV2160
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldb
      DOUBLE PRECISION Csl , Csr , Snl , Snr
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , Alphai(2) , Alphar(2) , B(Ldb,*) ,    &
     &                 Beta(2)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION anorm , ascale , bnorm , bscale , h1 , h2 , h3 , &
     &                 qq , r , rr , safmin , scale1 , scale2 , t ,     &
     &                 ulp , wi , wr1 , wr2
!     ..
!     .. External Subroutines ..
      EXTERNAL DLAG2 , DLARTG , DLASV2 , DROT
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLAPY2
      EXTERNAL DLAMCH , DLAPY2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
      safmin = DLAMCH('S')
      ulp = DLAMCH('P')
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
         CALL DLARTG(A(1,1),A(2,1),Csl,Snl,r)
         Csr = ONE
         Snr = ZERO
         CALL DROT(2,A(1,1),Lda,A(2,1),Lda,Csl,Snl)
         CALL DROT(2,B(1,1),Ldb,B(2,1),Ldb,Csl,Snl)
         A(2,1) = ZERO
         B(1,1) = ZERO
         B(2,1) = ZERO
         wi = ZERO
!
      ELSEIF ( ABS(B(2,2))<=ulp ) THEN
         CALL DLARTG(A(2,2),A(2,1),Csr,Snr,t)
         Snr = -Snr
         CALL DROT(2,A(1,1),1,A(1,2),1,Csr,Snr)
         CALL DROT(2,B(1,1),1,B(1,2),1,Csr,Snr)
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
         CALL DLAG2(A,Lda,B,Ldb,safmin,scale1,scale2,wr1,wr2,wi)
!
         IF ( wi==ZERO ) THEN
!
!           two real eigenvalues, compute s*A-w*B
!
            h1 = scale1*A(1,1) - wr1*B(1,1)
            h2 = scale1*A(1,2) - wr1*B(1,2)
            h3 = scale1*A(2,2) - wr1*B(2,2)
!
            rr = DLAPY2(h1,h2)
            qq = DLAPY2(scale1*A(2,1),h3)
!
            IF ( rr>qq ) THEN
!
!              find right rotation matrix to zero 1,1 element of
!              (sA - wB)
!
               CALL DLARTG(h2,h1,Csr,Snr,t)
!
            ELSE
!
!              find right rotation matrix to zero 2,1 element of
!              (sA - wB)
!
               CALL DLARTG(h3,scale1*A(2,1),Csr,Snr,t)
!
            ENDIF
!
            Snr = -Snr
            CALL DROT(2,A(1,1),1,A(1,2),1,Csr,Snr)
            CALL DROT(2,B(1,1),1,B(1,2),1,Csr,Snr)
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
               CALL DLARTG(B(1,1),B(2,1),Csl,Snl,r)
!
            ELSE
!
!              find left rotation matrix Q to zero out A(2,1)
!
               CALL DLARTG(A(1,1),A(2,1),Csl,Snl,r)
!
            ENDIF
!
            CALL DROT(2,A(1,1),Lda,A(2,1),Lda,Csl,Snl)
            CALL DROT(2,B(1,1),Ldb,B(2,1),Ldb,Csl,Snl)
!
            A(2,1) = ZERO
            B(2,1) = ZERO
!
         ELSE
!
!           a pair of complex conjugate eigenvalues
!           first compute the SVD of the matrix B
!
            CALL DLASV2(B(1,1),B(1,2),B(2,2),r,t,Snr,Csr,Snl,Csl)
!
!           Form (A,B) := Q(A,B)Z**T where Q is left rotation matrix and
!           Z is right rotation matrix computed from DLASV2
!
            CALL DROT(2,A(1,1),Lda,A(2,1),Lda,Csl,Snl)
            CALL DROT(2,B(1,1),Ldb,B(2,1),Ldb,Csl,Snl)
            CALL DROT(2,A(1,1),1,A(1,2),1,Csr,Snr)
            CALL DROT(2,B(1,1),1,B(1,2),1,Csr,Snr)
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
!     End of DLAGV2
!
      END SUBROUTINE DLAGV2
