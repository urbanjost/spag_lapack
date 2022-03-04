!*==alahd.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ALAHD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ALAHD( IOUNIT, PATH )
!
!       .. Scalar Arguments ..
!       CHARACTER*3        PATH
!       INTEGER            IOUNIT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ALAHD prints header information for the different test paths.
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
!>             _GE:  General matrices
!>             _GB:  General band
!>             _GT:  General Tridiagonal
!>             _PO:  Symmetric or Hermitian positive definite
!>             _PS:  Symmetric or Hermitian positive semi-definite
!>             _PP:  Symmetric or Hermitian positive definite packed
!>             _PB:  Symmetric or Hermitian positive definite band
!>             _PT:  Symmetric or Hermitian positive definite tridiagonal
!>             _SY:  Symmetric indefinite,
!>                     with partial (Bunch-Kaufman) pivoting
!>             _SR:  Symmetric indefinite,
!>                     with rook (bounded Bunch-Kaufman) pivoting
!>             _SK:  Symmetric indefinite,
!>                     with rook (bounded Bunch-Kaufman) pivoting
!>                     ( new storage format for factors:
!>                       L and diagonal of D is stored in A,
!>                       subdiagonal of D is stored in E )
!>             _SP:  Symmetric indefinite packed,
!>                     with partial (Bunch-Kaufman) pivoting
!>             _HA:  (complex) Hermitian ,
!>                     with Aasen Algorithm
!>             _HE:  (complex) Hermitian indefinite,
!>                     with partial (Bunch-Kaufman) pivoting
!>             _HR:  (complex) Hermitian indefinite,
!>                     with rook (bounded Bunch-Kaufman) pivoting
!>             _HK:  (complex) Hermitian indefinite,
!>                     with rook (bounded Bunch-Kaufman) pivoting
!>                     ( new storage format for factors:
!>                       L and diagonal of D is stored in A,
!>                       subdiagonal of D is stored in E )
!>             _HP:  (complex) Hermitian indefinite packed,
!>                     with partial (Bunch-Kaufman) pivoting
!>             _TR:  Triangular
!>             _TP:  Triangular packed
!>             _TB:  Triangular band
!>             _QR:  QR (general matrices)
!>             _LQ:  LQ (general matrices)
!>             _QL:  QL (general matrices)
!>             _RQ:  RQ (general matrices)
!>             _QP:  QR with column pivoting
!>             _TZ:  Trapezoidal
!>             _LS:  Least Squares driver routines
!>             _LU:  LU variants
!>             _CH:  Cholesky variants
!>             _QS:  QR variants
!>             _QT:  QRT (general matrices)
!>             _QX:  QRT (triangular-pentagonal matrices)
!>             _TS:  QR routines for tall-skinny and short-wide matrices
!>             _HH:  Householder reconstruction for tall-skinny matrices
!>          The first character must be one of S, D, C, or Z (C or Z only
!>          if complex).
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
!> \date June 2019
!
!> \ingroup aux_lin
!
!  =====================================================================
      SUBROUTINE ALAHD(Iounit,Path)
      IMPLICIT NONE
!*--ALAHD111
!
!  -- LAPACK test routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2019
!
!     .. Scalar Arguments ..
      CHARACTER*3 Path
      INTEGER Iounit
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL corz , sord
      CHARACTER c1 , c3
      CHARACTER*2 p2
      CHARACTER*4 eigcnm
      CHARACTER*32 subnam
      CHARACTER*9 sym
!     ..
!     .. External Functions ..
      LOGICAL LSAME , LSAMEN
      EXTERNAL LSAME , LSAMEN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC LEN_TRIM
!     ..
!     .. Executable Statements ..
!
      IF ( Iounit<=0 ) RETURN
      c1 = Path(1:1)
      c3 = Path(3:3)
      p2 = Path(2:3)
      sord = LSAME(c1,'S') .OR. LSAME(c1,'D')
      corz = LSAME(c1,'C') .OR. LSAME(c1,'Z')
      IF ( .NOT.(sord .OR. corz) ) RETURN
!
      IF ( LSAMEN(2,p2,'GE') ) THEN
!
!        GE: General dense
!
         WRITE (Iounit,FMT=99001) Path
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99029)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99048) 1
         WRITE (Iounit,FMT=99049) 2
         WRITE (Iounit,FMT=99050) 3
         WRITE (Iounit,FMT=99052) 4
         WRITE (Iounit,FMT=99053) 5
         WRITE (Iounit,FMT=99054) 6
         WRITE (Iounit,FMT=99055) 7
         WRITE (Iounit,FMT=99056) 8
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'GB') ) THEN
!
!        GB: General band
!
         WRITE (Iounit,FMT=99002) Path
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99030)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99048) 1
         WRITE (Iounit,FMT=99050) 2
         WRITE (Iounit,FMT=99052) 3
         WRITE (Iounit,FMT=99053) 4
         WRITE (Iounit,FMT=99054) 5
         WRITE (Iounit,FMT=99055) 6
         WRITE (Iounit,FMT=99056) 7
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'GT') ) THEN
!
!        GT: General tridiagonal
!
         WRITE (Iounit,FMT=99003) Path
         WRITE (Iounit,FMT=99031)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99048) 1
         WRITE (Iounit,FMT=99050) 2
         WRITE (Iounit,FMT=99052) 3
         WRITE (Iounit,FMT=99053) 4
         WRITE (Iounit,FMT=99054) 5
         WRITE (Iounit,FMT=99055) 6
         WRITE (Iounit,FMT=99056) 7
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'PO') .OR. LSAMEN(2,p2,'PP') ) THEN
!
!        PO: Positive definite full
!        PP: Positive definite packed
!
         IF ( sord ) THEN
            sym = 'Symmetric'
         ELSE
            sym = 'Hermitian'
         ENDIF
         IF ( LSAME(c3,'O') ) THEN
            WRITE (Iounit,FMT=99004) Path , sym
         ELSE
            WRITE (Iounit,FMT=99005) Path , sym
         ENDIF
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99033) Path
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99057) 1
         WRITE (Iounit,FMT=99049) 2
         WRITE (Iounit,FMT=99050) 3
         WRITE (Iounit,FMT=99052) 4
         WRITE (Iounit,FMT=99053) 5
         WRITE (Iounit,FMT=99054) 6
         WRITE (Iounit,FMT=99055) 7
         WRITE (Iounit,FMT=99056) 8
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'PS') ) THEN
!
!        PS: Positive semi-definite full
!
         IF ( sord ) THEN
            sym = 'Symmetric'
         ELSE
            sym = 'Hermitian'
         ENDIF
         IF ( LSAME(c1,'S') .OR. LSAME(c1,'C') ) THEN
            eigcnm = '1E04'
         ELSE
            eigcnm = '1D12'
         ENDIF
         WRITE (Iounit,FMT=99005) Path , sym
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99035) eigcnm , eigcnm , eigcnm
         WRITE (Iounit,FMT='( '' Difference:'' )')
         WRITE (Iounit,FMT=99036) c1
         WRITE (Iounit,FMT='( '' Test ratio:'' )')
         WRITE (Iounit,FMT=99058)
         WRITE (Iounit,FMT='( '' Messages:'' )')
      ELSEIF ( LSAMEN(2,p2,'PB') ) THEN
!
!        PB: Positive definite band
!
         IF ( sord ) THEN
            WRITE (Iounit,FMT=99006) Path , 'Symmetric'
         ELSE
            WRITE (Iounit,FMT=99006) Path , 'Hermitian'
         ENDIF
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99037) Path
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99057) 1
         WRITE (Iounit,FMT=99050) 2
         WRITE (Iounit,FMT=99052) 3
         WRITE (Iounit,FMT=99053) 4
         WRITE (Iounit,FMT=99054) 5
         WRITE (Iounit,FMT=99055) 6
         WRITE (Iounit,FMT=99056) 7
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'PT') ) THEN
!
!        PT: Positive definite tridiagonal
!
         IF ( sord ) THEN
            WRITE (Iounit,FMT=99007) Path , 'Symmetric'
         ELSE
            WRITE (Iounit,FMT=99007) Path , 'Hermitian'
         ENDIF
         WRITE (Iounit,FMT=99032)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99060) 1
         WRITE (Iounit,FMT=99050) 2
         WRITE (Iounit,FMT=99052) 3
         WRITE (Iounit,FMT=99053) 4
         WRITE (Iounit,FMT=99054) 5
         WRITE (Iounit,FMT=99055) 6
         WRITE (Iounit,FMT=99056) 7
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'SY') ) THEN
!
!        SY: Symmetric indefinite full,
!            with partial (Bunch-Kaufman) pivoting algorithm
!
         IF ( LSAME(c3,'Y') ) THEN
            WRITE (Iounit,FMT=99008) Path , 'Symmetric'
         ELSE
            WRITE (Iounit,FMT=99009) Path , 'Symmetric'
         ENDIF
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         IF ( sord ) THEN
            WRITE (Iounit,FMT=99038)
         ELSE
            WRITE (Iounit,FMT=99039)
         ENDIF
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99059) 1
         WRITE (Iounit,FMT=99049) 2
         WRITE (Iounit,FMT=99050) 3
         WRITE (Iounit,FMT=99050) 4
         WRITE (Iounit,FMT=99052) 5
         WRITE (Iounit,FMT=99053) 6
         WRITE (Iounit,FMT=99055) 7
         WRITE (Iounit,FMT=99054) 8
         WRITE (Iounit,FMT=99056) 9
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'SR') .OR. LSAMEN(2,p2,'SK') ) THEN
!
!        SR: Symmetric indefinite full,
!            with rook (bounded Bunch-Kaufman) pivoting algorithm
!
!        SK: Symmetric indefinite full,
!            with rook (bounded Bunch-Kaufman) pivoting algorithm,
!            ( new storage format for factors:
!              L and diagonal of D is stored in A,
!              subdiagonal of D is stored in E )
!
         WRITE (Iounit,FMT=99010) Path , 'Symmetric'
!
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         IF ( sord ) THEN
            WRITE (Iounit,FMT=99038)
         ELSE
            WRITE (Iounit,FMT=99039)
         ENDIF
!
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99059) 1
         WRITE (Iounit,FMT=99049) 2
         WRITE (Iounit,FMT=99086) 3
         WRITE (Iounit,FMT=99085)
         WRITE (Iounit,FMT=99087) 4
         WRITE (Iounit,FMT=99085)
         WRITE (Iounit,FMT=99050) 5
         WRITE (Iounit,FMT=99052) 6
         WRITE (Iounit,FMT=99056) 7
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'SP') ) THEN
!
!        SP: Symmetric indefinite packed,
!            with partial (Bunch-Kaufman) pivoting algorithm
!
         IF ( LSAME(c3,'Y') ) THEN
            WRITE (Iounit,FMT=99008) Path , 'Symmetric'
         ELSE
            WRITE (Iounit,FMT=99009) Path , 'Symmetric'
         ENDIF
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         IF ( sord ) THEN
            WRITE (Iounit,FMT=99038)
         ELSE
            WRITE (Iounit,FMT=99039)
         ENDIF
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99059) 1
         WRITE (Iounit,FMT=99049) 2
         WRITE (Iounit,FMT=99050) 3
         WRITE (Iounit,FMT=99052) 4
         WRITE (Iounit,FMT=99053) 5
         WRITE (Iounit,FMT=99055) 6
         WRITE (Iounit,FMT=99054) 7
         WRITE (Iounit,FMT=99056) 8
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'HA') ) THEN
!
!        HA: Hermitian,
!            with Assen Algorithm
!
         WRITE (Iounit,FMT=99008) Path , 'Hermitian'
!
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99038)
!
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99059) 1
         WRITE (Iounit,FMT=99049) 2
         WRITE (Iounit,FMT=99050) 3
         WRITE (Iounit,FMT=99050) 4
         WRITE (Iounit,FMT=99052) 5
         WRITE (Iounit,FMT=99053) 6
         WRITE (Iounit,FMT=99055) 7
         WRITE (Iounit,FMT=99054) 8
         WRITE (Iounit,FMT=99056) 9
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'HE') ) THEN
!
!        HE: Hermitian indefinite full,
!            with partial (Bunch-Kaufman) pivoting algorithm
!
         WRITE (Iounit,FMT=99008) Path , 'Hermitian'
!
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99038)
!
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99059) 1
         WRITE (Iounit,FMT=99049) 2
         WRITE (Iounit,FMT=99050) 3
         WRITE (Iounit,FMT=99050) 4
         WRITE (Iounit,FMT=99052) 5
         WRITE (Iounit,FMT=99053) 6
         WRITE (Iounit,FMT=99055) 7
         WRITE (Iounit,FMT=99054) 8
         WRITE (Iounit,FMT=99056) 9
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'HR') .OR. LSAMEN(2,p2,'HR') ) THEN
!
!        HR: Hermitian indefinite full,
!            with rook (bounded Bunch-Kaufman) pivoting algorithm
!
!        HK: Hermitian indefinite full,
!            with rook (bounded Bunch-Kaufman) pivoting algorithm,
!            ( new storage format for factors:
!              L and diagonal of D is stored in A,
!              subdiagonal of D is stored in E )
!
         WRITE (Iounit,FMT=99010) Path , 'Hermitian'
!
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99038)
!
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99059) 1
         WRITE (Iounit,FMT=99049) 2
         WRITE (Iounit,FMT=99086) 3
         WRITE (Iounit,FMT=99085)
         WRITE (Iounit,FMT=99087) 4
         WRITE (Iounit,FMT=99085)
         WRITE (Iounit,FMT=99050) 5
         WRITE (Iounit,FMT=99052) 6
         WRITE (Iounit,FMT=99056) 7
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'HP') ) THEN
!
!        HP: Hermitian indefinite packed,
!            with partial (Bunch-Kaufman) pivoting algorithm
!
         IF ( LSAME(c3,'E') ) THEN
            WRITE (Iounit,FMT=99008) Path , 'Hermitian'
         ELSE
            WRITE (Iounit,FMT=99009) Path , 'Hermitian'
         ENDIF
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99038)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99059) 1
         WRITE (Iounit,FMT=99049) 2
         WRITE (Iounit,FMT=99050) 3
         WRITE (Iounit,FMT=99052) 4
         WRITE (Iounit,FMT=99053) 5
         WRITE (Iounit,FMT=99055) 6
         WRITE (Iounit,FMT=99054) 7
         WRITE (Iounit,FMT=99056) 8
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'TR') .OR. LSAMEN(2,p2,'TP') ) THEN
!
!        TR: Triangular full
!        TP: Triangular packed
!
         IF ( LSAME(c3,'R') ) THEN
            WRITE (Iounit,FMT=99012) Path
            subnam = Path(1:1)//'LATRS'
         ELSE
            WRITE (Iounit,FMT=99013) Path
            subnam = Path(1:1)//'LATPS'
         ENDIF
         WRITE (Iounit,FMT=99044) Path
         WRITE (Iounit,FMT=99045) subnam(1:LEN_TRIM(subnam))
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99049) 1
         WRITE (Iounit,FMT=99050) 2
         WRITE (Iounit,FMT=99052) 3
         WRITE (Iounit,FMT=99053) 4
         WRITE (Iounit,FMT=99054) 5
         WRITE (Iounit,FMT=99055) 6
         WRITE (Iounit,FMT=99056) 7
         WRITE (Iounit,FMT=99061) subnam(1:LEN_TRIM(subnam)) , 8
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'TB') ) THEN
!
!        TB: Triangular band
!
         WRITE (Iounit,FMT=99014) Path
         subnam = Path(1:1)//'LATBS'
         WRITE (Iounit,FMT=99046) Path
         WRITE (Iounit,FMT=99047) subnam(1:LEN_TRIM(subnam))
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99050) 1
         WRITE (Iounit,FMT=99052) 2
         WRITE (Iounit,FMT=99053) 3
         WRITE (Iounit,FMT=99054) 4
         WRITE (Iounit,FMT=99055) 5
         WRITE (Iounit,FMT=99056) 6
         WRITE (Iounit,FMT=99061) subnam(1:LEN_TRIM(subnam)) , 7
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'QR') ) THEN
!
!        QR decomposition of rectangular matrices
!
         WRITE (Iounit,FMT=99015) Path , 'QR'
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99040)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99062) 1
         WRITE (Iounit,FMT=99063) 8
         WRITE (Iounit,FMT=99067) 2
         WRITE (Iounit,FMT=99069) 3 , 'M'
         WRITE (Iounit,FMT=99070) 4 , 'M'
         WRITE (Iounit,FMT=99071) 5 , 'M'
         WRITE (Iounit,FMT=99072) 6 , 'M'
         WRITE (Iounit,FMT=99050) 7
         WRITE (Iounit,FMT=99051) 9
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'LQ') ) THEN
!
!        LQ decomposition of rectangular matrices
!
         WRITE (Iounit,FMT=99015) Path , 'LQ'
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99040)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99064) 1
         WRITE (Iounit,FMT=99068) 2
         WRITE (Iounit,FMT=99069) 3 , 'N'
         WRITE (Iounit,FMT=99070) 4 , 'N'
         WRITE (Iounit,FMT=99071) 5 , 'N'
         WRITE (Iounit,FMT=99072) 6 , 'N'
         WRITE (Iounit,FMT=99050) 7
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'QL') ) THEN
!
!        QL decomposition of rectangular matrices
!
         WRITE (Iounit,FMT=99015) Path , 'QL'
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99040)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99065) 1
         WRITE (Iounit,FMT=99067) 2
         WRITE (Iounit,FMT=99069) 3 , 'M'
         WRITE (Iounit,FMT=99070) 4 , 'M'
         WRITE (Iounit,FMT=99071) 5 , 'M'
         WRITE (Iounit,FMT=99072) 6 , 'M'
         WRITE (Iounit,FMT=99050) 7
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'RQ') ) THEN
!
!        RQ decomposition of rectangular matrices
!
         WRITE (Iounit,FMT=99015) Path , 'RQ'
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99040)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99066) 1
         WRITE (Iounit,FMT=99068) 2
         WRITE (Iounit,FMT=99069) 3 , 'N'
         WRITE (Iounit,FMT=99070) 4 , 'N'
         WRITE (Iounit,FMT=99071) 5 , 'N'
         WRITE (Iounit,FMT=99072) 6 , 'N'
         WRITE (Iounit,FMT=99050) 7
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'QP') ) THEN
!
!        QR decomposition with column pivoting
!
         WRITE (Iounit,FMT=99016) Path
         WRITE (Iounit,FMT=99041)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99073) 1
         WRITE (Iounit,FMT=99074) 2
         WRITE (Iounit,FMT=99075) 3
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'TZ') ) THEN
!
!        TZ:  Trapezoidal
!
         WRITE (Iounit,FMT=99017) Path
         WRITE (Iounit,FMT=99042)
         WRITE (Iounit,FMT=99082) c1
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99073) 1
         WRITE (Iounit,FMT=99076) 2
         WRITE (Iounit,FMT=99075) 3
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'LS') ) THEN
!
!        LS:  Least Squares driver routines for
!             LS, LSD, LSS, LSX and LSY.
!
         WRITE (Iounit,FMT=99018) Path
         WRITE (Iounit,FMT=99043)
         WRITE (Iounit,FMT=99084) c1 , c1 , c1 , c1
         WRITE (Iounit,FMT=99077) 1
         WRITE (Iounit,FMT=99081) 2
         WRITE (Iounit,FMT=99079) 3
         WRITE (Iounit,FMT=99077) 4
         WRITE (Iounit,FMT=99078) 5
         WRITE (Iounit,FMT=99080) 6
         WRITE (Iounit,FMT=99083)
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'LU') ) THEN
!
!        LU factorization variants
!
         WRITE (Iounit,FMT=99019) Path
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99029)
         WRITE (Iounit,FMT='( '' Test ratio:'' )')
         WRITE (Iounit,FMT=99048) 1
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'CH') ) THEN
!
!        Cholesky factorization variants
!
         WRITE (Iounit,FMT=99020) Path
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99034)
         WRITE (Iounit,FMT='( '' Test ratio:'' )')
         WRITE (Iounit,FMT=99057) 1
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'QS') ) THEN
!
!        QR factorization variants
!
         WRITE (Iounit,FMT=99021) Path
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99040)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'QT') ) THEN
!
!        QRT (general matrices)
!
         WRITE (Iounit,FMT=99023) Path
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99088) 1
         WRITE (Iounit,FMT=99089) 2
         WRITE (Iounit,FMT=99090) 3
         WRITE (Iounit,FMT=99091) 4
         WRITE (Iounit,FMT=99092) 5
         WRITE (Iounit,FMT=99093) 6
!
      ELSEIF ( LSAMEN(2,p2,'QX') ) THEN
!
!        QRT (triangular-pentagonal)
!
         WRITE (Iounit,FMT=99024) Path
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99094) 1
         WRITE (Iounit,FMT=99095) 2
         WRITE (Iounit,FMT=99096) 3
         WRITE (Iounit,FMT=99097) 4
         WRITE (Iounit,FMT=99098) 5
         WRITE (Iounit,FMT=99099) 6
!
      ELSEIF ( LSAMEN(2,p2,'TQ') ) THEN
!
!        QRT (triangular-pentagonal)
!
         WRITE (Iounit,FMT=99025) Path
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99100) 1
         WRITE (Iounit,FMT=99101) 2
         WRITE (Iounit,FMT=99102) 3
         WRITE (Iounit,FMT=99103) 4
         WRITE (Iounit,FMT=99104) 5
         WRITE (Iounit,FMT=99105) 6
!
      ELSEIF ( LSAMEN(2,p2,'XQ') ) THEN
!
!        QRT (triangular-pentagonal)
!
         WRITE (Iounit,FMT=99026) Path
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99106) 1
         WRITE (Iounit,FMT=99107) 2
         WRITE (Iounit,FMT=99108) 3
         WRITE (Iounit,FMT=99109) 4
         WRITE (Iounit,FMT=99110) 5
         WRITE (Iounit,FMT=99111) 6
!
      ELSEIF ( LSAMEN(2,p2,'TS') ) THEN
!
!        TS:  QR routines for tall-skinny and short-wide matrices
!
         WRITE (Iounit,FMT=99027) Path
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99112) 1
         WRITE (Iounit,FMT=99113) 2
         WRITE (Iounit,FMT=99114) 3
         WRITE (Iounit,FMT=99115) 4
         WRITE (Iounit,FMT=99116) 5
         WRITE (Iounit,FMT=99117) 6
!
      ELSEIF ( LSAMEN(2,p2,'HH') ) THEN
!
!        HH:  Householder reconstruction for tall-skinny matrices
!
         WRITE (Iounit,FMT=99028) Path
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99118) 1
         WRITE (Iounit,FMT=99119) 2
         WRITE (Iounit,FMT=99120) 3
         WRITE (Iounit,FMT=99121) 4
         WRITE (Iounit,FMT=99122) 5
         WRITE (Iounit,FMT=99123) 6
!
      ELSE
!
!        Print error message if no header is available.
!
         WRITE (Iounit,FMT=99022) Path
      ENDIF
!
!     First line of header
!
99001 FORMAT (/1X,A3,':  General dense matrices')
99002 FORMAT (/1X,A3,':  General band matrices')
99003 FORMAT (/1X,A3,':  General tridiagonal')
99004 FORMAT (/1X,A3,':  ',A9,' positive definite matrices')
99005 FORMAT (/1X,A3,':  ',A9,' positive definite packed matrices')
99006 FORMAT (/1X,A3,':  ',A9,' positive definite band matrices')
99007 FORMAT (/1X,A3,':  ',A9,' positive definite tridiagonal')
99008 FORMAT (/1X,A3,':  ',A9,' indefinite matrices',                   &
     &        ', partial (Bunch-Kaufman) pivoting')
99009 FORMAT (/1X,A3,':  ',A9,' indefinite packed matrices',            &
     &        ', partial (Bunch-Kaufman) pivoting')
99010 FORMAT (/1X,A3,':  ',A9,' indefinite matrices',                   &
     &        ', "rook" (bounded Bunch-Kaufman) pivoting')
99011 FORMAT (/1X,A3,':  ',A9,' indefinite packed matrices',            &
     &        ', "rook" (bounded Bunch-Kaufman) pivoting')
99012 FORMAT (/1X,A3,':  Triangular matrices')
99013 FORMAT (/1X,A3,':  Triangular packed matrices')
99014 FORMAT (/1X,A3,':  Triangular band matrices')
99015 FORMAT (/1X,A3,':  ',A2,' factorization of general matrices')
99016 FORMAT (/1X,A3,':  QR factorization with column pivoting')
99017 FORMAT (/1X,A3,':  RQ factorization of trapezoidal matrix')
99018 FORMAT (/1X,A3,':  Least squares driver routines')
99019 FORMAT (/1X,A3,':  LU factorization variants')
99020 FORMAT (/1X,A3,':  Cholesky factorization variants')
99021 FORMAT (/1X,A3,':  QR factorization variants')
99022 FORMAT (/1X,A3,':  No header available')
99023 FORMAT (/1X,A3,':  QRT factorization for general matrices')
99024 FORMAT (/1X,A3,':  QRT factorization for ',                       &
     &        'triangular-pentagonal matrices')
99025 FORMAT (/1X,A3,':  LQT factorization for general matrices')
99026 FORMAT (/1X,A3,':  LQT factorization for ',                       &
     &        'triangular-pentagonal matrices')
99027 FORMAT (/1X,A3,':  TS factorization for ',                        &
     &        'tall-skinny or short-wide matrices')
99028 FORMAT (/1X,A3,':  Householder recostruction from TSQR',          &
     &        ' factorization output ',/,' for tall-skinny matrices.')
!
!     GE matrix types
!
99029 FORMAT (4X,'1. Diagonal',24X,'7. Last n/2 columns zero',/4X,      &
     &        '2. Upper triangular',16X,                                &
     &        '8. Random, CNDNUM = sqrt(0.1/EPS)',/4X,                  &
     &        '3. Lower triangular',16X,'9. Random, CNDNUM = 0.1/EPS',  &
     &        /4X,'4. Random, CNDNUM = 2',13X,                          &
     &        '10. Scaled near underflow',/4X,'5. First column zero',   &
     &        14X,'11. Scaled near overflow',/4X,'6. Last column zero')
!
!     GB matrix types
!
99030 FORMAT (4X,'1. Random, CNDNUM = 2',14X,                           &
     &        '5. Random, CNDNUM = sqrt(0.1/EPS)',/4X,                  &
     &        '2. First column zero',15X,'6. Random, CNDNUM = .01/EPS', &
     &        /4X,'3. Last column zero',16X,'7. Scaled near underflow', &
     &        /4X,'4. Last n/2 columns zero',11X,                       &
     &        '8. Scaled near overflow')
!
!     GT matrix types
!
99031 FORMAT (' Matrix types (1-6 have specified condition numbers):',  &
     &        /4X,'1. Diagonal',24X,'7. Random, unspecified CNDNUM',/4X,&
     &        '2. Random, CNDNUM = 2',14X,'8. First column zero',/4X,   &
     &        '3. Random, CNDNUM = sqrt(0.1/EPS)',2X,                   &
     &        '9. Last column zero',/4X,'4. Random, CNDNUM = 0.1/EPS',  &
     &        7X,'10. Last n/2 columns zero',/4X,                       &
     &        '5. Scaled near underflow',10X,                           &
     &        '11. Scaled near underflow',/4X,'6. Scaled near overflow',&
     &        11X,'12. Scaled near overflow')
!
!     PT matrix types
!
99032 FORMAT (' Matrix types (1-6 have specified condition numbers):',  &
     &        /4X,'1. Diagonal',24X,'7. Random, unspecified CNDNUM',/4X,&
     &        '2. Random, CNDNUM = 2',14X,                              &
     &        '8. First row and column zero',/4X,                       &
     &        '3. Random, CNDNUM = sqrt(0.1/EPS)',2X,                   &
     &        '9. Last row and column zero',/4X,                        &
     &        '4. Random, CNDNUM = 0.1/EPS',7X,                         &
     &        '10. Middle row and column zero',/4X,                     &
     &        '5. Scaled near underflow',10X,                           &
     &        '11. Scaled near underflow',/4X,'6. Scaled near overflow',&
     &        11X,'12. Scaled near overflow')
!
!     PO, PP matrix types
!
99033 FORMAT (4X,'1. Diagonal',24X,'6. Random, CNDNUM = sqrt(0.1/EPS)', &
     &        /4X,'2. Random, CNDNUM = 2',14X,                          &
     &        '7. Random, CNDNUM = 0.1/EPS',/3X,                        &
     &        '*3. First row and column zero',7X,                       &
     &        '8. Scaled near underflow',/3X,                           &
     &        '*4. Last row and column zero',8X,                        &
     &        '9. Scaled near overflow',/3X,                            &
     &        '*5. Middle row and column zero',/3X,                     &
     &        '(* - tests error exits from ',A3,                        &
     &        'TRF, no test ratios are computed)')
!
!     CH matrix types
!
99034 FORMAT (4X,'1. Diagonal',24X,'6. Random, CNDNUM = sqrt(0.1/EPS)', &
     &        /4X,'2. Random, CNDNUM = 2',14X,                          &
     &        '7. Random, CNDNUM = 0.1/EPS',/3X,                        &
     &        '*3. First row and column zero',7X,                       &
     &        '8. Scaled near underflow',/3X,                           &
     &        '*4. Last row and column zero',8X,                        &
     &        '9. Scaled near overflow',/3X,                            &
     &        '*5. Middle row and column zero',/3X,                     &
     &        '(* - tests error exits, no test ratios are computed)')
!
!     PS matrix types
!
99035 FORMAT (4X,'1. Diagonal',/4X,'2. Random, CNDNUM = 2',14X,/3X,     &
     &        '*3. Nonzero eigenvalues of: D(1:RANK-1)=1 and ',         &
     &        'D(RANK) = 1.0/',A4,/3X,                                  &
     &        '*4. Nonzero eigenvalues of: D(1)=1 and ',                &
     &        ' D(2:RANK) = 1.0/',A4,/3X,                               &
     &        '*5. Nonzero eigenvalues of: D(I) = ',A4,                 &
     &        '**(-(I-1)/(RANK-1)) ',' I=1:RANK',/4X,                   &
     &        '6. Random, CNDNUM = sqrt(0.1/EPS)',/4X,                  &
     &        '7. Random, CNDNUM = 0.1/EPS',/4X,                        &
     &        '8. Scaled near underflow',/4X,'9. Scaled near overflow', &
     &        /3X,'(* - Semi-definite tests )')
99036 FORMAT (3X,'RANK minus computed rank, returned by ',A,'PSTRF')
!
!     PB matrix types
!
99037 FORMAT (4X,'1. Random, CNDNUM = 2',14X,                           &
     &        '5. Random, CNDNUM = sqrt(0.1/EPS)',/3X,                  &
     &        '*2. First row and column zero',7X,                       &
     &        '6. Random, CNDNUM = 0.1/EPS',/3X,                        &
     &        '*3. Last row and column zero',8X,                        &
     &        '7. Scaled near underflow',/3X,                           &
     &        '*4. Middle row and column zero',6X,                      &
     &        '8. Scaled near overflow',/3X,                            &
     &        '(* - tests error exits from ',A3,                        &
     &        'TRF, no test ratios are computed)')
!
!     SSY, SSR, SSP, CHE, CHR, CHP matrix types
!
99038 FORMAT (4X,'1. Diagonal',24X,'6. Last n/2 rows and columns zero', &
     &        /4X,'2. Random, CNDNUM = 2',14X,                          &
     &        '7. Random, CNDNUM = sqrt(0.1/EPS)',/4X,                  &
     &        '3. First row and column zero',7X,                        &
     &        '8. Random, CNDNUM = 0.1/EPS',/4X,                        &
     &        '4. Last row and column zero',8X,                         &
     &        '9. Scaled near underflow',/4X,                           &
     &        '5. Middle row and column zero',5X,                       &
     &        '10. Scaled near overflow')
!
!     CSY, CSR, CSP matrix types
!
99039 FORMAT (4X,'1. Diagonal',24X,'7. Random, CNDNUM = sqrt(0.1/EPS)', &
     &        /4X,'2. Random, CNDNUM = 2',14X,                          &
     &        '8. Random, CNDNUM = 0.1/EPS',/4X,                        &
     &        '3. First row and column zero',7X,                        &
     &        '9. Scaled near underflow',/4X,                           &
     &        '4. Last row and column zero',7X,                         &
     &        '10. Scaled near overflow',/4X,                           &
     &        '5. Middle row and column zero',5X,                       &
     &        '11. Block diagonal matrix',/4X,                          &
     &        '6. Last n/2 rows and columns zero')
!
!     QR matrix types
!
99040 FORMAT (4X,'1. Diagonal',24X,'5. Random, CNDNUM = sqrt(0.1/EPS)', &
     &        /4X,'2. Upper triangular',16X,                            &
     &        '6. Random, CNDNUM = 0.1/EPS',/4X,'3. Lower triangular',  &
     &        16X,'7. Scaled near underflow',/4X,                       &
     &        '4. Random, CNDNUM = 2',14X,'8. Scaled near overflow')
!
!     QP matrix types
!
99041 FORMAT (' Matrix types (2-6 have condition 1/EPS):',/4X,          &
     &        '1. Zero matrix',21X,'4. First n/2 columns fixed',/4X,    &
     &        '2. One small eigenvalue',12X,'5. Last n/2 columns fixed',&
     &        /4X,'3. Geometric distribution',10X,                      &
     &        '6. Every second column fixed')
!
!     TZ matrix types
!
99042 FORMAT (' Matrix types (2-3 have condition 1/EPS):',/4X,          &
     &        '1. Zero matrix',/4X,'2. One small eigenvalue',/4X,       &
     &        '3. Geometric distribution')
!
!     LS matrix types
!
99043 FORMAT (' Matrix types (1-3: full rank, 4-6: rank deficient):',   &
     &        /4X,'1 and 4. Normal scaling',/4X,                        &
     &        '2 and 5. Scaled near overflow',/4X,                      &
     &        '3 and 6. Scaled near underflow')
!
!     TR, TP matrix types
!
99044 FORMAT (' Matrix types for ',A3,' routines:',/4X,'1. Diagonal',   &
     &        24X,'6. Scaled near overflow',/4X,'2. Random, CNDNUM = 2',&
     &        14X,'7. Identity',/4X,                                    &
     &        '3. Random, CNDNUM = sqrt(0.1/EPS)  ',                    &
     &        '8. Unit triangular, CNDNUM = 2',/4X,                     &
     &        '4. Random, CNDNUM = 0.1/EPS',8X,                         &
     &        '9. Unit, CNDNUM = sqrt(0.1/EPS)',/4X,                    &
     &        '5. Scaled near underflow',10X,                           &
     &        '10. Unit, CNDNUM = 0.1/EPS')
99045 FORMAT (' Special types for testing ',A,':',/3X,                  &
     &        '11. Matrix elements are O(1), large right hand side',/3X,&
     &        '12. First diagonal causes overflow,',                    &
     &        ' offdiagonal column norms < 1',/3X,                      &
     &        '13. First diagonal causes overflow,',                    &
     &        ' offdiagonal column norms > 1',/3X,                      &
     &        '14. Growth factor underflows, solution does not overflow'&
     &        ,/3X,'15. Small diagonal causes gradual overflow',/3X,    &
     &        '16. One zero diagonal element',/3X,                      &
     &      '17. Large offdiagonals cause overflow when adding a column'&
     &      ,/3X,'18. Unit triangular with large right hand side')
!
!     TB matrix types
!
99046 FORMAT (' Matrix types for ',A3,' routines:',/4X,                 &
     &        '1. Random, CNDNUM = 2',14X,'6. Identity',/4X,            &
     &        '2. Random, CNDNUM = sqrt(0.1/EPS)  ',                    &
     &        '7. Unit triangular, CNDNUM = 2',/4X,                     &
     &        '3. Random, CNDNUM = 0.1/EPS',8X,                         &
     &        '8. Unit, CNDNUM = sqrt(0.1/EPS)',/4X,                    &
     &        '4. Scaled near underflow',11X,                           &
     &        '9. Unit, CNDNUM = 0.1/EPS',/4X,'5. Scaled near overflow')
99047 FORMAT (' Special types for testing ',A,':',/3X,                  &
     &        '10. Matrix elements are O(1), large right hand side',/3X,&
     &        '11. First diagonal causes overflow,',                    &
     &        ' offdiagonal column norms < 1',/3X,                      &
     &        '12. First diagonal causes overflow,',                    &
     &        ' offdiagonal column norms > 1',/3X,                      &
     &        '13. Growth factor underflows, solution does not overflow'&
     &        ,/3X,'14. Small diagonal causes gradual overflow',/3X,    &
     &        '15. One zero diagonal element',/3X,                      &
     &      '16. Large offdiagonals cause overflow when adding a column'&
     &      ,/3X,'17. Unit triangular with large right hand side')
!
!     Test ratios
!
99048 FORMAT (3X,I2,': norm( L * U - A )  / ( N * norm(A) * EPS )')
99049 FORMAT (3X,I2,': norm( I - A*AINV ) / ',                          &
     &        '( N * norm(A) * norm(AINV) * EPS )')
99050 FORMAT (3X,I2,': norm( B - A * X )  / ',                          &
     &        '( norm(A) * norm(X) * EPS )')
99051 FORMAT (3X,I2,': diagonal is not non-negative')
99052 FORMAT (3X,I2,': norm( X - XACT )   / ',                          &
     &        '( norm(XACT) * CNDNUM * EPS )')
99053 FORMAT (3X,I2,': norm( X - XACT )   / ',                          &
     &        '( norm(XACT) * CNDNUM * EPS ), refined')
99054 FORMAT (3X,I2,': norm( X - XACT )   / ',                          &
     &        '( norm(XACT) * (error bound) )')
99055 FORMAT (3X,I2,': (backward error)   / EPS')
99056 FORMAT (3X,I2,': RCOND * CNDNUM - 1.0')
99057 FORMAT (3X,I2,': norm( U'' * U - A ) / ( N * norm(A) * EPS )',    &
     &        ', or',/7X,'norm( L * L'' - A ) / ( N * norm(A) * EPS )')
99058 FORMAT (3X,'norm( P * U'' * U * P'' - A ) / ( N * norm(A) * EPS )'&
     &        ,', or',/3X,                                              &
     &        'norm( P * L * L'' * P'' - A ) / ( N * norm(A) * EPS )')
99059 FORMAT (3X,I2,': norm( U*D*U'' - A ) / ( N * norm(A) * EPS )',    &
     &        ', or',/7X,'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )')
99060 FORMAT (3X,I2,': norm( U''*D*U - A ) / ( N * norm(A) * EPS )',    &
     &        ', or',/7X,'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )')
99061 FORMAT (' Test ratio for ',A,':',/3X,I2,                          &
     &        ': norm( s*b - A*x )  / ( norm(A) * norm(x) * EPS )')
99062 FORMAT (3X,I2,': norm( R - Q'' * A ) / ( M * norm(A) * EPS )')
99063 FORMAT (3X,I2,                                                    &
     &': norm( R - Q'' * A ) / ( M * norm(A) * EPS )                    &
     &                                                                  &
     &                                                                  &
     &                                     [RFPG]')
99064 FORMAT (3X,I2,': norm( L - A * Q'' ) / ( N * norm(A) * EPS )')
99065 FORMAT (3X,I2,': norm( L - Q'' * A ) / ( M * norm(A) * EPS )')
99066 FORMAT (3X,I2,': norm( R - A * Q'' ) / ( N * norm(A) * EPS )')
99067 FORMAT (3X,I2,': norm( I - Q''*Q )   / ( M * EPS )')
99068 FORMAT (3X,I2,': norm( I - Q*Q'' )   / ( N * EPS )')
99069 FORMAT (3X,I2,': norm( Q*C - Q*C )  / ','( ',A1,                  &
     &        ' * norm(C) * EPS )')
99070 FORMAT (3X,I2,': norm( C*Q - C*Q )  / ','( ',A1,                  &
     &        ' * norm(C) * EPS )')
99071 FORMAT (3X,I2,': norm( Q''*C - Q''*C )/ ','( ',A1,                &
     &        ' * norm(C) * EPS )')
99072 FORMAT (3X,I2,': norm( C*Q'' - C*Q'' )/ ','( ',A1,                &
     &        ' * norm(C) * EPS )')
99073 FORMAT (3X,I2,': norm(svd(A) - svd(R)) / ',                       &
     &        '( M * norm(svd(R)) * EPS )')
99074 FORMAT (3X,I2,': norm( A*P - Q*R )     / ( M * norm(A) * EPS )')
99075 FORMAT (3X,I2,': norm( I - Q''*Q )      / ( M * EPS )')
99076 FORMAT (3X,I2,': norm( A - R*Q )       / ( M * norm(A) * EPS )')
99077 FORMAT (3X,I2,': norm( B - A * X )   / ',                         &
     &        '( max(M,N) * norm(A) * norm(X) * EPS )')
99078 FORMAT (3X,I2,': norm( (A*X-B)'' *A ) / ',                        &
     &        '( max(M,N,NRHS) * norm(A) * norm(B) * EPS )')
99079 FORMAT (3X,I2,': norm(svd(A)-svd(R)) / ',                         &
     &        '( min(M,N) * norm(svd(R)) * EPS )')
99080 FORMAT (3X,I2,': Check if X is in the row space of A or A''')
99081 FORMAT (3X,I2,': norm( (A*X-B)'' *A ) / ',                        &
     &        '( max(M,N,NRHS) * norm(A) * norm(B) * EPS )',/7X,        &
     &        'if TRANS=''N'' and M.GE.N or TRANS=''T'' and M.LT.N, ',  &
     &        'otherwise',/7X,                                          &
     &        'check if X is in the row space of A or A'' ',            &
     &        '(overdetermined case)')
99082 FORMAT (' Test ratios (1-3: ',A1,'TZRZF):')
99083 FORMAT (3X,' 7-10: same as 3-6',3X,' 11-14: same as 3-6')
99084 FORMAT (' Test ratios:',/'    (1-2: ',A1,'GELS, 3-6: ',A1,        &
     &        'GELSY, 7-10: ',A1,'GELSS, 11-14: ',A1,'GELSD, 15-16: ',  &
     &        A1,'GETSLS)')
99085 FORMAT (7X,'where ALPHA = ( 1 + SQRT( 17 ) ) / 8')
99086 FORMAT (3X,I2,': ABS( Largest element in L )',/12X,               &
     &        ' - ( 1 / ( 1 - ALPHA ) ) + THRESH')
99087 FORMAT (3X,I2,': Largest 2-Norm of 2-by-2 pivots',/12X,           &
     &        ' - ( ( 1 + ALPHA ) / ( 1 - ALPHA ) ) + THRESH')
99088 FORMAT (3X,I2,': norm( R - Q''*A ) / ( M * norm(A) * EPS )')
99089 FORMAT (3X,I2,': norm( I - Q''*Q ) / ( M * EPS )')
99090 FORMAT (3X,I2,': norm( Q*C - Q*C ) / ( M * norm(C) * EPS )')
99091 FORMAT (3X,I2,': norm( Q''*C - Q''*C ) / ( M * norm(C) * EPS )')
99092 FORMAT (3X,I2,': norm( C*Q - C*Q ) / ( M * norm(C) * EPS )')
99093 FORMAT (3X,I2,': norm( C*Q'' - C*Q'' ) / ( M * norm(C) * EPS )')
99094 FORMAT (3X,I2,': norm( R - Q''*A ) / ( (M+N) * norm(A) * EPS )')
99095 FORMAT (3X,I2,': norm( I - Q''*Q ) / ( (M+N) * EPS )')
99096 FORMAT (3X,I2,': norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )')
99097 FORMAT (3X,I2,                                                    &
     &        ': norm( Q''*C - Q''*C ) / ( (M+N) * norm(C) * EPS )')
99098 FORMAT (3X,I2,': norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )')
99099 FORMAT (3X,I2,                                                    &
     &        ': norm( C*Q'' - C*Q'' ) / ( (M+N) * norm(C) * EPS )')
99100 FORMAT (3X,I2,': norm( L - A*Q'' ) / ( (M+N) * norm(A) * EPS )')
99101 FORMAT (3X,I2,': norm( I - Q*Q'' ) / ( (M+N) * EPS )')
99102 FORMAT (3X,I2,': norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )')
99103 FORMAT (3X,I2,                                                    &
     &        ': norm( Q''*C - Q''*C ) / ( (M+N) * norm(C) * EPS )')
99104 FORMAT (3X,I2,': norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )')
99105 FORMAT (3X,I2,                                                    &
     &        ': norm( C*Q'' - C*Q'' ) / ( (M+N) * norm(C) * EPS )')
99106 FORMAT (3X,I2,': norm( L - A*Q'' ) / ( (M+N) * norm(A) * EPS )')
99107 FORMAT (3X,I2,': norm( I - Q*Q'' ) / ( (M+N) * EPS )')
99108 FORMAT (3X,I2,': norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )')
99109 FORMAT (3X,I2,                                                    &
     &        ': norm( Q''*C - Q''*C ) / ( (M+N) * norm(C) * EPS )')
99110 FORMAT (3X,I2,': norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )')
99111 FORMAT (3X,I2,                                                    &
     &        ': norm( C*Q'' - C*Q'' ) / ( (M+N) * norm(C) * EPS )')
99112 FORMAT (3X,I2,': norm( R - Q''*A ) / ( (M+N) * norm(A) * EPS )')
99113 FORMAT (3X,I2,': norm( I - Q''*Q ) / ( (M+N) * EPS )')
99114 FORMAT (3X,I2,': norm( Q*C - Q*C ) / ( (M+N) * norm(C) * EPS )')
99115 FORMAT (3X,I2,                                                    &
     &        ': norm( Q''*C - Q''*C ) / ( (M+N) * norm(C) * EPS )')
99116 FORMAT (3X,I2,': norm( C*Q - C*Q ) / ( (M+N) * norm(C) * EPS )')
99117 FORMAT (3X,I2,                                                    &
     &        ': norm( C*Q'' - C*Q'' ) / ( (M+N) * norm(C) * EPS )')
!
99118 FORMAT (3X,I2,': norm( R - Q''*A ) / ( M * norm(A) * EPS )')
99119 FORMAT (3X,I2,': norm( I - Q''*Q ) / ( M * EPS )')
99120 FORMAT (3X,I2,': norm( Q*C - Q*C ) / ( M * norm(C) * EPS )')
99121 FORMAT (3X,I2,': norm( Q''*C - Q''*C ) / ( M * norm(C) * EPS )')
99122 FORMAT (3X,I2,': norm( C*Q - C*Q ) / ( M * norm(C) * EPS )')
99123 FORMAT (3X,I2,': norm( C*Q'' - C*Q'' ) / ( M * norm(C) * EPS )')
 
!
!
!     End of ALAHD
!
      END SUBROUTINE ALAHD
