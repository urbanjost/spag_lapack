!*==alahdg.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b alahdg
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE alahdg( IOUNIT, PATH )
!
!       .. Scalar Arguments ..
!       CHARACTER*3       PATH
!       INTEGER           IOUNIT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> alahdg prints header information for the different test paths.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IOUNIT
!> \verbatim
!>          IOUNIT is INTEGER
!>          The unit number to which the header information should be
!>          printed.
!> \endverbatim
!>
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3
!>          The name of the path for which the header information is to
!>          be printed.  Current paths are
!>             GQR:  GQR (general matrices)
!>             GRQ:  GRQ (general matrices)
!>             LSE:  LSE Problem
!>             GLM:  GLM Problem
!>             GSV:  Generalized Singular Value Decomposition
!>             CSD:  CS Decomposition
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
!> \ingroup aux_eig
!
!  =====================================================================
      SUBROUTINE alahdg(Iounit,Path)
      IMPLICIT NONE
!*--ALAHDG66
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*3 Path
      INTEGER Iounit
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      CHARACTER*3 c2
      INTEGER itype
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. Executable Statements ..
!
      IF ( Iounit<=0 ) RETURN
      c2 = Path(1:3)
!
!     First line describing matrices in this path
!
      IF ( LSAMEN(3,c2,'GQR') ) THEN
         itype = 1
         WRITE (Iounit,FMT=99002) Path
      ELSEIF ( LSAMEN(3,c2,'GRQ') ) THEN
         itype = 2
         WRITE (Iounit,FMT=99003) Path
      ELSEIF ( LSAMEN(3,c2,'LSE') ) THEN
         itype = 3
         WRITE (Iounit,FMT=99004) Path
      ELSEIF ( LSAMEN(3,c2,'GLM') ) THEN
         itype = 4
         WRITE (Iounit,FMT=99005) Path
      ELSEIF ( LSAMEN(3,c2,'GSV') ) THEN
         itype = 5
         WRITE (Iounit,FMT=99006) Path
      ELSEIF ( LSAMEN(3,c2,'CSD') ) THEN
         itype = 6
         WRITE (Iounit,FMT=99007) Path
      ENDIF
!
!     Matrix types
!
      WRITE (Iounit,FMT=99001) 'Matrix types: '
!
      IF ( itype==1 ) THEN
         WRITE (Iounit,FMT=99008) 1
         WRITE (Iounit,FMT=99010) 2
         WRITE (Iounit,FMT=99012) 3
         WRITE (Iounit,FMT=99013) 4
         WRITE (Iounit,FMT=99014) 5
         WRITE (Iounit,FMT=99015) 6
         WRITE (Iounit,FMT=99018) 7
         WRITE (Iounit,FMT=99019) 8
      ELSEIF ( itype==2 ) THEN
         WRITE (Iounit,FMT=99009) 1
         WRITE (Iounit,FMT=99011) 2
         WRITE (Iounit,FMT=99012) 3
         WRITE (Iounit,FMT=99013) 4
         WRITE (Iounit,FMT=99014) 5
         WRITE (Iounit,FMT=99015) 6
         WRITE (Iounit,FMT=99018) 7
         WRITE (Iounit,FMT=99019) 8
      ELSEIF ( itype==3 ) THEN
         WRITE (Iounit,FMT=99008) 1
         WRITE (Iounit,FMT=99010) 2
         WRITE (Iounit,FMT=99012) 3
         WRITE (Iounit,FMT=99013) 4
         WRITE (Iounit,FMT=99013) 5
         WRITE (Iounit,FMT=99013) 6
         WRITE (Iounit,FMT=99013) 7
         WRITE (Iounit,FMT=99013) 8
      ELSEIF ( itype==4 ) THEN
         WRITE (Iounit,FMT=99009) 1
         WRITE (Iounit,FMT=99011) 2
         WRITE (Iounit,FMT=99012) 3
         WRITE (Iounit,FMT=99013) 4
         WRITE (Iounit,FMT=99013) 5
         WRITE (Iounit,FMT=99013) 6
         WRITE (Iounit,FMT=99013) 7
         WRITE (Iounit,FMT=99013) 8
      ELSEIF ( itype==5 ) THEN
         WRITE (Iounit,FMT=99008) 1
         WRITE (Iounit,FMT=99010) 2
         WRITE (Iounit,FMT=99012) 3
         WRITE (Iounit,FMT=99013) 4
         WRITE (Iounit,FMT=99014) 5
         WRITE (Iounit,FMT=99015) 6
         WRITE (Iounit,FMT=99016) 7
         WRITE (Iounit,FMT=99017) 8
      ELSEIF ( itype==6 ) THEN
         WRITE (Iounit,FMT=99020) 1
         WRITE (Iounit,FMT=99021) 2
         WRITE (Iounit,FMT=99022) 3
      ENDIF
!
!     Tests performed
!
      WRITE (Iounit,FMT=99001) 'Test ratios: '
!
      IF ( itype==1 ) THEN
!
!        GQR decomposition of rectangular matrices
!
         WRITE (Iounit,FMT=99023) 1
         WRITE (Iounit,FMT=99024) 2
         WRITE (Iounit,FMT=99025) 3
         WRITE (Iounit,FMT=99026) 4
      ELSEIF ( itype==2 ) THEN
!
!        GRQ decomposition of rectangular matrices
!
         WRITE (Iounit,FMT=99027) 1
         WRITE (Iounit,FMT=99028) 2
         WRITE (Iounit,FMT=99025) 3
         WRITE (Iounit,FMT=99026) 4
      ELSEIF ( itype==3 ) THEN
!
!        LSE Problem
!
         WRITE (Iounit,FMT=99029) 1
         WRITE (Iounit,FMT=99030) 2
      ELSEIF ( itype==4 ) THEN
!
!        GLM Problem
!
         WRITE (Iounit,FMT=99031) 1
      ELSEIF ( itype==5 ) THEN
!
!        GSVD
!
         WRITE (Iounit,FMT=99032) 1
         WRITE (Iounit,FMT=99033) 2
         WRITE (Iounit,FMT=99034) 3
         WRITE (Iounit,FMT=99035) 4
         WRITE (Iounit,FMT=99036) 5
      ELSEIF ( itype==6 ) THEN
!
!        CSD
!
         WRITE (Iounit,FMT=99037)
         WRITE (Iounit,FMT=99038) 1
         WRITE (Iounit,FMT=99039) 2
         WRITE (Iounit,FMT=99040) 3
         WRITE (Iounit,FMT=99041) 4
         WRITE (Iounit,FMT=99042) 5
         WRITE (Iounit,FMT=99043) 6
         WRITE (Iounit,FMT=99044) 7
         WRITE (Iounit,FMT=99045) 8
         WRITE (Iounit,FMT=99046) 9
         WRITE (Iounit,FMT=99047)
         WRITE (Iounit,FMT=99048) 10
         WRITE (Iounit,FMT=99049) 11
         WRITE (Iounit,FMT=99050) 12
         WRITE (Iounit,FMT=99051) 13
         WRITE (Iounit,FMT=99052) 14
         WRITE (Iounit,FMT=99053) 15
      ENDIF
!
99001 FORMAT (1X,A)
99002 FORMAT (/1X,A3,': GQR factorization of general matrices')
99003 FORMAT (/1X,A3,': GRQ factorization of general matrices')
99004 FORMAT (/1X,A3,': LSE Problem')
99005 FORMAT (/1X,A3,': GLM Problem')
99006 FORMAT (/1X,A3,': Generalized Singular Value Decomposition')
99007 FORMAT (/1X,A3,': CS Decomposition')
!
99008 FORMAT (3X,I2,': A-diagonal matrix  B-upper triangular')
99009 FORMAT (3X,I2,': A-diagonal matrix  B-lower triangular')
99010 FORMAT (3X,I2,': A-upper triangular B-upper triangular')
99011 FORMAT (3X,I2,': A-lower triangular B-diagonal triangular')
99012 FORMAT (3X,I2,': A-lower triangular B-upper triangular')
!
99013 FORMAT (3X,I2,': Random matrices cond(A)=100, cond(B)=10,')
!
99014 FORMAT (3X,I2,': Random matrices cond(A)= sqrt( 0.1/EPS ) ',      &
     &        'cond(B)= sqrt( 0.1/EPS )')
99015 FORMAT (3X,I2,': Random matrices cond(A)= 0.1/EPS ',              &
     &        'cond(B)= 0.1/EPS')
99016 FORMAT (3X,I2,': Random matrices cond(A)= sqrt( 0.1/EPS ) ',      &
     &        'cond(B)=  0.1/EPS ')
99017 FORMAT (3X,I2,': Random matrices cond(A)= 0.1/EPS ',              &
     &        'cond(B)=  sqrt( 0.1/EPS )')
!
99018 FORMAT (3X,I2,': Matrix scaled near underflow limit')
99019 FORMAT (3X,I2,': Matrix scaled near overflow limit')
99020 FORMAT (3X,I2,': Random orthogonal matrix (Haar measure)')
99021 FORMAT (3X,I2,': Nearly orthogonal matrix with uniformly ',       &
     &        'distributed angles atan2( S, C ) in CS decomposition')
99022 FORMAT (3X,I2,': Random orthogonal matrix with clustered ',       &
     &        'angles atan2( S, C ) in CS decomposition')
!
!
!     GQR test ratio
!
99023 FORMAT (3X,I2,': norm( R - Q'' * A ) / ( min( N, M )*norm( A )',  &
     &        '* EPS )')
99024 FORMAT (3X,I2,': norm( T * Z - Q'' * B )  / ( min(P,N)*norm(B)',  &
     &        '* EPS )')
99025 FORMAT (3X,I2,': norm( I - Q''*Q )   / ( N * EPS )')
99026 FORMAT (3X,I2,': norm( I - Z''*Z )   / ( P * EPS )')
!
!     GRQ test ratio
!
99027 FORMAT (3X,I2,': norm( R - A * Q'' ) / ( min( N,M )*norm(A) * ',  &
     &        'EPS )')
99028 FORMAT (3X,I2,': norm( T * Q - Z'' * B )  / ( min( P,N ) * nor',  &
     &        'm(B)*EPS )')
!
!     LSE test ratio
!
99029 FORMAT (3X,I2,': norm( A*x - c )  / ( norm(A)*norm(x) * EPS )')
99030 FORMAT (3X,I2,': norm( B*x - d )  / ( norm(B)*norm(x) * EPS )')
!
!     GLM test ratio
!
99031 FORMAT (3X,I2,': norm( d - A*x - B*y ) / ( (norm(A)+norm(B) )*',  &
     &        '(norm(x)+norm(y))*EPS )')
!
!     GSVD test ratio
!
99032 FORMAT (3X,I2,': norm( U'' * A * Q - D1 * R ) / ( min( M, N )*',  &
     &        'norm( A ) * EPS )')
99033 FORMAT (3X,I2,': norm( V'' * B * Q - D2 * R ) / ( min( P, N )*',  &
     &        'norm( B ) * EPS )')
99034 FORMAT (3X,I2,': norm( I - U''*U )   / ( M * EPS )')
99035 FORMAT (3X,I2,': norm( I - V''*V )   / ( P * EPS )')
99036 FORMAT (3X,I2,': norm( I - Q''*Q )   / ( N * EPS )')
!
!     CSD test ratio
!
99037 FORMAT (3X,'2-by-2 CSD')
99038 FORMAT (3X,I2,': norm( U1'' * X11 * V1 - C ) / ( max(  P,  Q)',   &
     &        ' * max(norm(I-X''*X),EPS) )')
99039 FORMAT (3X,I2,': norm( U1'' * X12 * V2-(-S)) / ( max(  P,',       &
     &        'M-Q) * max(norm(I-X''*X),EPS) )')
99040 FORMAT (3X,I2,': norm( U2'' * X21 * V1 - S ) / ( max(M-P,',       &
     &        '  Q) * max(norm(I-X''*X),EPS) )')
99041 FORMAT (3X,I2,': norm( U2'' * X22 * V2 - C ) / ( max(M-P,',       &
     &        'M-Q) * max(norm(I-X''*X),EPS) )')
99042 FORMAT (3X,I2,': norm( I - U1''*U1 ) / (   P   * EPS )')
99043 FORMAT (3X,I2,': norm( I - U2''*U2 ) / ( (M-P) * EPS )')
99044 FORMAT (3X,I2,': norm( I - V1''*V1 ) / (   Q   * EPS )')
99045 FORMAT (3X,I2,': norm( I - V2''*V2 ) / ( (M-Q) * EPS )')
99046 FORMAT (3X,I2,': principal angle ordering ( 0 or ULP )')
99047 FORMAT (3X,'2-by-1 CSD')
99048 FORMAT (3X,I2,': norm( U1'' * X11 * V1 - C ) / ( max(  P,  Q)',   &
     &        ' * max(norm(I-X''*X),EPS) )')
99049 FORMAT (3X,I2,': norm( U2'' * X21 * V1 - S ) / ( max(  M-P,',     &
     &        'Q) * max(norm(I-X''*X),EPS) )')
99050 FORMAT (3X,I2,': norm( I - U1''*U1 ) / (   P   * EPS )')
99051 FORMAT (3X,I2,': norm( I - U2''*U2 ) / ( (M-P) * EPS )')
99052 FORMAT (3X,I2,': norm( I - V1''*V1 ) / (   Q   * EPS )')
99053 FORMAT (3X,I2,': principal angle ordering ( 0 or ULP )')
!
!     End of alahdg
!
      END SUBROUTINE alahdg
