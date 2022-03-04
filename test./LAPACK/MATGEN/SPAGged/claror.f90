!*==claror.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CLAROR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAROR( SIDE, INIT, M, N, A, LDA, ISEED, X, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          INIT, SIDE
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       COMPLEX            A( LDA, * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CLAROR pre- or post-multiplies an M by N matrix A by a random
!>    unitary matrix U, overwriting A. A may optionally be
!>    initialized to the identity matrix before multiplying by U.
!>    U is generated using the method of G.W. Stewart
!>    ( SIAM J. Numer. Anal. 17, 1980, pp. 403-409 ).
!>    (BLAS-2 version)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           SIDE specifies whether A is multiplied on the left or right
!>           by U.
!>       SIDE = 'L'   Multiply A on the left (premultiply) by U
!>       SIDE = 'R'   Multiply A on the right (postmultiply) by UC>       SIDE = 'C'   Multiply A on the left by U and the right by UC>       SIDE = 'T'   Multiply A on the left by U and the right by U'
!>           Not modified.
!> \endverbatim
!>
!> \param[in] INIT
!> \verbatim
!>          INIT is CHARACTER*1
!>           INIT specifies whether or not A should be initialized to
!>           the identity matrix.
!>              INIT = 'I'   Initialize A to (a section of) the
!>                           identity matrix before applying U.
!>              INIT = 'N'   No initialization.  Apply U to the
!>                           input matrix A.
!>
!>           INIT = 'I' may be used to generate square (i.e., unitary)
!>           or rectangular orthogonal matrices (orthogonality being
!>           in the sense of CDOTC):
!>
!>           For square matrices, M=N, and SIDE many be either 'L' or
!>           'R'; the rows will be orthogonal to each other, as will the
!>           columns.
!>           For rectangular matrices where M < N, SIDE = 'R' will
!>           produce a dense matrix whose rows will be orthogonal and
!>           whose columns will not, while SIDE = 'L' will produce a
!>           matrix whose rows will be orthogonal, and whose first M
!>           columns will be orthogonal, the remaining columns being
!>           zero.
!>           For matrices where M > N, just use the previous
!>           explanation, interchanging 'L' and 'R' and "rows" and
!>           "columns".
!>
!>           Not modified.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           Number of rows of A. Not modified.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           Number of columns of A. Not modified.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension ( LDA, N )
!>           Input and output array. Overwritten by U A ( if SIDE = 'L' )
!>           or by A U ( if SIDE = 'R' )
!>           or by U A U* ( if SIDE = 'C')
!>           or by U A U' ( if SIDE = 'T') on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           Leading dimension of A. Must be at least MAX ( 1, M ).
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension ( 4 )
!>           On entry ISEED specifies the seed of the random number
!>           generator. The array elements should be between 0 and 4095;
!>           if not they will be reduced mod 4096.  Also, ISEED(4) must
!>           be odd.  The random number generator uses a linear
!>           congruential sequence limited to small integers, and so
!>           should produce machine independent random numbers. The
!>           values of ISEED are changed on exit, and can be used in the
!>           next call to CLAROR to continue the same random number
!>           sequence.
!>           Modified.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension ( 3*MAX( M, N ) )
!>           Workspace. Of length:
!>               2*M + N if SIDE = 'L',
!>               2*N + M if SIDE = 'R',
!>               3*N     if SIDE = 'C' or 'T'.
!>           Modified.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           An error flag.  It is set to:
!>            0  if no error.
!>            1  if CLARND returned a bad random number (installation
!>               problem)
!>           -1  if SIDE is not L, R, C, or T.
!>           -3  if M is negative.
!>           -4  if N is negative or if SIDE is C or T and N is not equal
!>               to M.
!>           -6  if LDA is less than M.
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
!> \ingroup complex_matgen
!
!  =====================================================================
      SUBROUTINE CLAROR(Side,Init,M,N,A,Lda,Iseed,X,Info)
      IMPLICIT NONE
!*--CLAROR162
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Init , Side
      INTEGER Info , Lda , M , N
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      COMPLEX A(Lda,*) , X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TOOSML
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0,TOOSML=1.0E-20)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER irow , itype , ixfrm , j , jcol , kbeg , nxfrm
      REAL factor , xabs , xnorm
      COMPLEX csign , xnorms
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SCNRM2
      COMPLEX CLARND
      EXTERNAL LSAME , SCNRM2 , CLARND
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMV , CGERC , CLACGV , CLASET , CSCAL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CMPLX , CONJG
!     ..
!     .. Executable Statements ..
!
      Info = 0
      IF ( N==0 .OR. M==0 ) RETURN
!
      itype = 0
      IF ( LSAME(Side,'L') ) THEN
         itype = 1
      ELSEIF ( LSAME(Side,'R') ) THEN
         itype = 2
      ELSEIF ( LSAME(Side,'C') ) THEN
         itype = 3
      ELSEIF ( LSAME(Side,'T') ) THEN
         itype = 4
      ENDIF
!
!     Check for argument errors.
!
      IF ( itype==0 ) THEN
         Info = -1
      ELSEIF ( M<0 ) THEN
         Info = -3
      ELSEIF ( N<0 .OR. (itype==3 .AND. N/=M) ) THEN
         Info = -4
      ELSEIF ( Lda<M ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CLAROR',-Info)
         RETURN
      ENDIF
!
      IF ( itype==1 ) THEN
         nxfrm = M
      ELSE
         nxfrm = N
      ENDIF
!
!     Initialize A to the identity matrix if desired
!
      IF ( LSAME(Init,'I') ) CALL CLASET('Full',M,N,CZERO,CONE,A,Lda)
!
!     If no rotation possible, still multiply by
!     a random complex number from the circle |x| = 1
!
!      2)      Compute Rotation by computing Householder
!              Transformations H(2), H(3), ..., H(n).  Note that the
!              order in which they are computed is irrelevant.
!
      DO j = 1 , nxfrm
         X(j) = CZERO
      ENDDO
!
      DO ixfrm = 2 , nxfrm
         kbeg = nxfrm - ixfrm + 1
!
!        Generate independent normal( 0, 1 ) random numbers
!
         DO j = kbeg , nxfrm
            X(j) = CLARND(3,Iseed)
         ENDDO
!
!        Generate a Householder transformation from the random vector X
!
         xnorm = SCNRM2(ixfrm,X(kbeg),1)
         xabs = ABS(X(kbeg))
         IF ( xabs/=CZERO ) THEN
            csign = X(kbeg)/xabs
         ELSE
            csign = CONE
         ENDIF
         xnorms = csign*xnorm
         X(nxfrm+kbeg) = -csign
         factor = xnorm*(xnorm+xabs)
         IF ( ABS(factor)<TOOSML ) THEN
            Info = 1
            CALL XERBLA('CLAROR',-Info)
            RETURN
         ELSE
            factor = ONE/factor
         ENDIF
         X(kbeg) = X(kbeg) + xnorms
!
!        Apply Householder transformation to A
!
         IF ( itype==1 .OR. itype==3 .OR. itype==4 ) THEN
!
!           Apply H(k) on the left of A
!
            CALL CGEMV('C',ixfrm,N,CONE,A(kbeg,1),Lda,X(kbeg),1,CZERO,  &
     &                 X(2*nxfrm+1),1)
            CALL CGERC(ixfrm,N,-CMPLX(factor),X(kbeg),1,X(2*nxfrm+1),1, &
     &                 A(kbeg,1),Lda)
!
         ENDIF
!
         IF ( itype>=2 .AND. itype<=4 ) THEN
!
!           Apply H(k)* (or H(k)') on the right of A
!
            IF ( itype==4 ) CALL CLACGV(ixfrm,X(kbeg),1)
!
            CALL CGEMV('N',M,ixfrm,CONE,A(1,kbeg),Lda,X(kbeg),1,CZERO,  &
     &                 X(2*nxfrm+1),1)
            CALL CGERC(M,ixfrm,-CMPLX(factor),X(2*nxfrm+1),1,X(kbeg),1, &
     &                 A(1,kbeg),Lda)
!
         ENDIF
      ENDDO
!
      X(1) = CLARND(3,Iseed)
      xabs = ABS(X(1))
      IF ( xabs/=ZERO ) THEN
         csign = X(1)/xabs
      ELSE
         csign = CONE
      ENDIF
      X(2*nxfrm) = csign
!
!     Scale the matrix A by D.
!
      IF ( itype==1 .OR. itype==3 .OR. itype==4 ) THEN
         DO irow = 1 , M
            CALL CSCAL(N,CONJG(X(nxfrm+irow)),A(irow,1),Lda)
         ENDDO
      ENDIF
!
      IF ( itype==2 .OR. itype==3 ) THEN
         DO jcol = 1 , N
            CALL CSCAL(M,X(nxfrm+jcol),A(1,jcol),1)
         ENDDO
      ENDIF
!
      IF ( itype==4 ) THEN
         DO jcol = 1 , N
            CALL CSCAL(M,CONJG(X(nxfrm+jcol)),A(1,jcol),1)
         ENDDO
      ENDIF
!
!     End of CLAROR
!
      END SUBROUTINE CLAROR
