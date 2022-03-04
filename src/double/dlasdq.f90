!*==dlasdq.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLASDQ computes the SVD of a real bidiagonal matrix with diagonal d and off-diagonal e. Used by sbdsdc.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASDQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasdq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasdq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasdq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASDQ( UPLO, SQRE, N, NCVT, NRU, NCC, D, E, VT, LDVT,
!                          U, LDU, C, LDC, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU, SQRE
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ),
!      $                   VT( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASDQ computes the singular value decomposition (SVD) of a real
!> (upper or lower) bidiagonal matrix with diagonal D and offdiagonal
!> E, accumulating the transformations if desired. Letting B denote
!> the input bidiagonal matrix, the algorithm computes orthogonal
!> matrices Q and P such that B = Q * S * P**T (P**T denotes the transpose
!> of P). The singular values S are overwritten on D.
!>
!> The input matrix U  is changed to U  * Q  if desired.
!> The input matrix VT is changed to P**T * VT if desired.
!> The input matrix C  is changed to Q**T * C  if desired.
!>
!> See "Computing  Small Singular Values of Bidiagonal Matrices With
!> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
!> LAPACK Working Note #3, for a detailed description of the algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>        On entry, UPLO specifies whether the input bidiagonal matrix
!>        is upper or lower bidiagonal, and whether it is square are
!>        not.
!>           UPLO = 'U' or 'u'   B is upper bidiagonal.
!>           UPLO = 'L' or 'l'   B is lower bidiagonal.
!> \endverbatim
!>
!> \param[in] SQRE
!> \verbatim
!>          SQRE is INTEGER
!>        = 0: then the input matrix is N-by-N.
!>        = 1: then the input matrix is N-by-(N+1) if UPLU = 'U' and
!>             (N+1)-by-N if UPLU = 'L'.
!>
!>        The bidiagonal matrix has
!>        N = NL + NR + 1 rows and
!>        M = N + SQRE >= N columns.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>        On entry, N specifies the number of rows and columns
!>        in the matrix. N must be at least 0.
!> \endverbatim
!>
!> \param[in] NCVT
!> \verbatim
!>          NCVT is INTEGER
!>        On entry, NCVT specifies the number of columns of
!>        the matrix VT. NCVT must be at least 0.
!> \endverbatim
!>
!> \param[in] NRU
!> \verbatim
!>          NRU is INTEGER
!>        On entry, NRU specifies the number of rows of
!>        the matrix U. NRU must be at least 0.
!> \endverbatim
!>
!> \param[in] NCC
!> \verbatim
!>          NCC is INTEGER
!>        On entry, NCC specifies the number of columns of
!>        the matrix C. NCC must be at least 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>        On entry, D contains the diagonal entries of the
!>        bidiagonal matrix whose SVD is desired. On normal exit,
!>        D contains the singular values in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array.
!>        dimension is (N-1) if SQRE = 0 and N if SQRE = 1.
!>        On entry, the entries of E contain the offdiagonal entries
!>        of the bidiagonal matrix whose SVD is desired. On normal
!>        exit, E will contain 0. If the algorithm does not converge,
!>        D and E will contain the diagonal and superdiagonal entries
!>        of a bidiagonal matrix orthogonally equivalent to the one
!>        given as input.
!> \endverbatim
!>
!> \param[in,out] VT
!> \verbatim
!>          VT is DOUBLE PRECISION array, dimension (LDVT, NCVT)
!>        On entry, contains a matrix which on exit has been
!>        premultiplied by P**T, dimension N-by-NCVT if SQRE = 0
!>        and (N+1)-by-NCVT if SQRE = 1 (not referenced if NCVT=0).
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>        On entry, LDVT specifies the leading dimension of VT as
!>        declared in the calling (sub) program. LDVT must be at
!>        least 1. If NCVT is nonzero LDVT must also be at least N.
!> \endverbatim
!>
!> \param[in,out] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU, N)
!>        On entry, contains a  matrix which on exit has been
!>        postmultiplied by Q, dimension NRU-by-N if SQRE = 0
!>        and NRU-by-(N+1) if SQRE = 1 (not referenced if NRU=0).
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>        On entry, LDU  specifies the leading dimension of U as
!>        declared in the calling (sub) program. LDU must be at
!>        least max( 1, NRU ) .
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC, NCC)
!>        On entry, contains an N-by-NCC matrix which on exit
!>        has been premultiplied by Q**T  dimension N-by-NCC if SQRE = 0
!>        and (N+1)-by-NCC if SQRE = 1 (not referenced if NCC=0).
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>        On entry, LDC  specifies the leading dimension of C as
!>        declared in the calling (sub) program. LDC must be at
!>        least 1. If NCC is nonzero, LDC must also be at least N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (4*N)
!>        Workspace. Only referenced if one of NCVT, NRU, or NCC is
!>        nonzero, and if N is at least 2.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>        On exit, a value of 0 indicates a successful exit.
!>        If INFO < 0, argument number -INFO is illegal.
!>        If INFO > 0, the algorithm did not converge, and INFO
!>        specifies how many superdiagonals did not converge.
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
!> \date June 2016
!
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Huan Ren, Computer Science Division, University of
!>     California at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE DLASDQ(Uplo,Sqre,N,Ncvt,Nru,Ncc,D,E,Vt,Ldvt,U,Ldu,C,   &
     &                  Ldc,Work,Info)
      IMPLICIT NONE
!*--DLASDQ215
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Ldc , Ldu , Ldvt , N , Ncc , Ncvt , Nru , Sqre
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION C(Ldc,*) , D(*) , E(*) , U(Ldu,*) , Vt(Ldvt,*) , &
     &                 Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL rotate
      INTEGER i , isub , iuplo , j , np1 , sqre1
      DOUBLE PRECISION cs , r , smin , sn
!     ..
!     .. External Subroutines ..
      EXTERNAL DBDSQR , DLARTG , DLASR , DSWAP , XERBLA
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      iuplo = 0
      IF ( LSAME(Uplo,'U') ) iuplo = 1
      IF ( LSAME(Uplo,'L') ) iuplo = 2
      IF ( iuplo==0 ) THEN
         Info = -1
      ELSEIF ( (Sqre<0) .OR. (Sqre>1) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Ncvt<0 ) THEN
         Info = -4
      ELSEIF ( Nru<0 ) THEN
         Info = -5
      ELSEIF ( Ncc<0 ) THEN
         Info = -6
      ELSEIF ( (Ncvt==0 .AND. Ldvt<1) .OR. (Ncvt>0 .AND. Ldvt<MAX(1,N)) &
     &         ) THEN
         Info = -10
      ELSEIF ( Ldu<MAX(1,Nru) ) THEN
         Info = -12
      ELSEIF ( (Ncc==0 .AND. Ldc<1) .OR. (Ncc>0 .AND. Ldc<MAX(1,N)) )   &
     &         THEN
         Info = -14
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DLASDQ',-Info)
         RETURN
      ENDIF
      IF ( N==0 ) RETURN
!
!     ROTATE is true if any singular vectors desired, false otherwise
!
      rotate = (Ncvt>0) .OR. (Nru>0) .OR. (Ncc>0)
      np1 = N + 1
      sqre1 = Sqre
!
!     If matrix non-square upper bidiagonal, rotate to be lower
!     bidiagonal.  The rotations are on the right.
!
      IF ( (iuplo==1) .AND. (sqre1==1) ) THEN
         DO i = 1 , N - 1
            CALL DLARTG(D(i),E(i),cs,sn,r)
            D(i) = r
            E(i) = sn*D(i+1)
            D(i+1) = cs*D(i+1)
            IF ( rotate ) THEN
               Work(i) = cs
               Work(N+i) = sn
            ENDIF
         ENDDO
         CALL DLARTG(D(N),E(N),cs,sn,r)
         D(N) = r
         E(N) = ZERO
         IF ( rotate ) THEN
            Work(N) = cs
            Work(N+N) = sn
         ENDIF
         iuplo = 2
         sqre1 = 0
!
!        Update singular vectors if desired.
!
         IF ( Ncvt>0 ) CALL DLASR('L','V','F',np1,Ncvt,Work(1),Work(np1)&
     &                            ,Vt,Ldvt)
      ENDIF
!
!     If matrix lower bidiagonal, rotate to be upper bidiagonal
!     by applying Givens rotations on the left.
!
      IF ( iuplo==2 ) THEN
         DO i = 1 , N - 1
            CALL DLARTG(D(i),E(i),cs,sn,r)
            D(i) = r
            E(i) = sn*D(i+1)
            D(i+1) = cs*D(i+1)
            IF ( rotate ) THEN
               Work(i) = cs
               Work(N+i) = sn
            ENDIF
         ENDDO
!
!        If matrix (N+1)-by-N lower bidiagonal, one additional
!        rotation is needed.
!
         IF ( sqre1==1 ) THEN
            CALL DLARTG(D(N),E(N),cs,sn,r)
            D(N) = r
            IF ( rotate ) THEN
               Work(N) = cs
               Work(N+N) = sn
            ENDIF
         ENDIF
!
!        Update singular vectors if desired.
!
         IF ( Nru>0 ) THEN
            IF ( sqre1==0 ) THEN
               CALL DLASR('R','V','F',Nru,N,Work(1),Work(np1),U,Ldu)
            ELSE
               CALL DLASR('R','V','F',Nru,np1,Work(1),Work(np1),U,Ldu)
            ENDIF
         ENDIF
         IF ( Ncc>0 ) THEN
            IF ( sqre1==0 ) THEN
               CALL DLASR('L','V','F',N,Ncc,Work(1),Work(np1),C,Ldc)
            ELSE
               CALL DLASR('L','V','F',np1,Ncc,Work(1),Work(np1),C,Ldc)
            ENDIF
         ENDIF
      ENDIF
!
!     Call DBDSQR to compute the SVD of the reduced real
!     N-by-N upper bidiagonal matrix.
!
      CALL DBDSQR('U',N,Ncvt,Nru,Ncc,D,E,Vt,Ldvt,U,Ldu,C,Ldc,Work,Info)
!
!     Sort the singular values into ascending order (insertion sort on
!     singular values, but only one transposition per singular vector)
!
      DO i = 1 , N
!
!        Scan for smallest D(I).
!
         isub = i
         smin = D(i)
         DO j = i + 1 , N
            IF ( D(j)<smin ) THEN
               isub = j
               smin = D(j)
            ENDIF
         ENDDO
         IF ( isub/=i ) THEN
!
!           Swap singular values and vectors.
!
            D(isub) = D(i)
            D(i) = smin
            IF ( Ncvt>0 ) CALL DSWAP(Ncvt,Vt(isub,1),Ldvt,Vt(i,1),Ldvt)
            IF ( Nru>0 ) CALL DSWAP(Nru,U(1,isub),1,U(1,i),1)
            IF ( Ncc>0 ) CALL DSWAP(Ncc,C(isub,1),Ldc,C(i,1),Ldc)
         ENDIF
      ENDDO
!
!
!     End of DLASDQ
!
      END SUBROUTINE DLASDQ
