!*==aladhd.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b aladhd
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE aladhd( IOUNIT, PATH )
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
!> aladhd prints header information for the driver routines test paths.
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
!>                     Assen Algorithm
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
!> \date December 2016
!
!> \ingroup aux_lin
!
!  =====================================================================
      SUBROUTINE aladhd(Iounit,Path)
      IMPLICIT NONE
!*--ALADHD94
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
      LOGICAL corz , sord
      CHARACTER c1 , c3
      CHARACTER*2 p2
      CHARACTER*9 sym
!     ..
!     .. External Functions ..
      LOGICAL LSAME , LSAMEN
EXTERNAL    LSAME,   LSAMEN
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
         WRITE (Iounit,FMT=99013)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99021) 1
         WRITE (Iounit,FMT=99022) 2
         WRITE (Iounit,FMT=99023) 3
         WRITE (Iounit,FMT=99024) 4
         WRITE (Iounit,FMT=99025) 5
         WRITE (Iounit,FMT=99026) 6
         WRITE (Iounit,FMT=99030) 7
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'GB') ) THEN
!
!        GB: General band
!
         WRITE (Iounit,FMT=99002) Path
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99014)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99021) 1
         WRITE (Iounit,FMT=99022) 2
         WRITE (Iounit,FMT=99023) 3
         WRITE (Iounit,FMT=99024) 4
         WRITE (Iounit,FMT=99025) 5
         WRITE (Iounit,FMT=99026) 6
         WRITE (Iounit,FMT=99030) 7
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'GT') ) THEN
!
!        GT: General tridiagonal
!
         WRITE (Iounit,FMT=99003) Path
         WRITE (Iounit,FMT=99015)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99021) 1
         WRITE (Iounit,FMT=99022) 2
         WRITE (Iounit,FMT=99023) 3
         WRITE (Iounit,FMT=99024) 4
         WRITE (Iounit,FMT=99025) 5
         WRITE (Iounit,FMT=99026) 6
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'PO') .OR. LSAMEN(2,p2,'PP') .OR.            &
     &         LSAMEN(2,p2,'PS') ) THEN
!
!        PO: Positive definite full
!        PS: Positive definite full
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
         WRITE (Iounit,FMT=99017) Path
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99027) 1
         WRITE (Iounit,FMT=99022) 2
         WRITE (Iounit,FMT=99023) 3
         WRITE (Iounit,FMT=99024) 4
         WRITE (Iounit,FMT=99025) 5
         WRITE (Iounit,FMT=99026) 6
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
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
         WRITE (Iounit,FMT=99018) Path
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99027) 1
         WRITE (Iounit,FMT=99022) 2
         WRITE (Iounit,FMT=99023) 3
         WRITE (Iounit,FMT=99024) 4
         WRITE (Iounit,FMT=99025) 5
         WRITE (Iounit,FMT=99026) 6
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
         WRITE (Iounit,FMT=99016)
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99029) 1
         WRITE (Iounit,FMT=99022) 2
         WRITE (Iounit,FMT=99023) 3
         WRITE (Iounit,FMT=99024) 4
         WRITE (Iounit,FMT=99025) 5
         WRITE (Iounit,FMT=99026) 6
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'SY') .OR. LSAMEN(2,p2,'SP') ) THEN
!
!        SY: Symmetric indefinite full
!            with partial (Bunch-Kaufman) pivoting algorithm
!        SP: Symmetric indefinite packed
!            with partial (Bunch-Kaufman) pivoting algorithm
!
         IF ( LSAME(c3,'Y') ) THEN
            WRITE (Iounit,FMT=99009) Path , 'Symmetric'
         ELSE
            WRITE (Iounit,FMT=99010) Path , 'Symmetric'
         ENDIF
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         IF ( sord ) THEN
            WRITE (Iounit,FMT=99019)
         ELSE
            WRITE (Iounit,FMT=99020)
         ENDIF
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99028) 1
         WRITE (Iounit,FMT=99022) 2
         WRITE (Iounit,FMT=99023) 3
         WRITE (Iounit,FMT=99025) 4
         WRITE (Iounit,FMT=99024) 5
         WRITE (Iounit,FMT=99026) 6
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
         WRITE (Iounit,FMT=99009) Path , 'Symmetric'
!
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         IF ( sord ) THEN
            WRITE (Iounit,FMT=99019)
         ELSE
            WRITE (Iounit,FMT=99020)
         ENDIF
!
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99028) 1
         WRITE (Iounit,FMT=99022) 2
         WRITE (Iounit,FMT=99023) 3
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'HA') ) THEN
!
!        HA: Hermitian
!            Aasen algorithm
         WRITE (Iounit,FMT=99008) Path , 'Hermitian'
!
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99019)
!
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99028) 1
         WRITE (Iounit,FMT=99022) 2
         WRITE (Iounit,FMT=99023) 3
         WRITE (Iounit,FMT=99025) 4
         WRITE (Iounit,FMT=99024) 5
         WRITE (Iounit,FMT=99026) 6
         WRITE (Iounit,FMT='( '' Messages:'' )')
 
 
      ELSEIF ( LSAMEN(2,p2,'HE') .OR. LSAMEN(2,p2,'HP') ) THEN
!
!        HE: Hermitian indefinite full
!            with partial (Bunch-Kaufman) pivoting algorithm
!        HP: Hermitian indefinite packed
!            with partial (Bunch-Kaufman) pivoting algorithm
!
         IF ( LSAME(c3,'E') ) THEN
            WRITE (Iounit,FMT=99009) Path , 'Hermitian'
         ELSE
            WRITE (Iounit,FMT=99010) Path , 'Hermitian'
         ENDIF
!
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99019)
!
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99028) 1
         WRITE (Iounit,FMT=99022) 2
         WRITE (Iounit,FMT=99023) 3
         WRITE (Iounit,FMT=99025) 4
         WRITE (Iounit,FMT=99024) 5
         WRITE (Iounit,FMT=99026) 6
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSEIF ( LSAMEN(2,p2,'HR') .OR. LSAMEN(2,p2,'HK') ) THEN
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
         WRITE (Iounit,FMT=99009) Path , 'Hermitian'
!
         WRITE (Iounit,FMT='( '' Matrix types:'' )')
         WRITE (Iounit,FMT=99019)
!
         WRITE (Iounit,FMT='( '' Test ratios:'' )')
         WRITE (Iounit,FMT=99028) 1
         WRITE (Iounit,FMT=99022) 2
         WRITE (Iounit,FMT=99023) 3
         WRITE (Iounit,FMT='( '' Messages:'' )')
!
      ELSE
!
!        Print error message if no header is available.
!
         WRITE (Iounit,FMT=99012) Path
      ENDIF
!
!     First line of header
!
99001 FORMAT (/1X,A3,' drivers:  General dense matrices')
99002 FORMAT (/1X,A3,' drivers:  General band matrices')
99003 FORMAT (/1X,A3,' drivers:  General tridiagonal')
99004 FORMAT (/1X,A3,' drivers:  ',A9,' positive definite matrices')
99005 FORMAT (/1X,A3,' drivers:  ',A9,                                  &
     &        ' positive definite packed matrices')
99006 FORMAT (/1X,A3,' drivers:  ',A9,                                  &
     &        ' positive definite band matrices')
99007 FORMAT (/1X,A3,' drivers:  ',A9,' positive definite tridiagonal')
99008 FORMAT (/1X,A3,' drivers:  ',A9,' indefinite matrices',           &
     &        ', "Aasen" Algorithm')
99009 FORMAT (/1X,A3,' drivers:  ',A9,' indefinite matrices',           &
     &        ', "rook" (bounded Bunch-Kaufman) pivoting')
99010 FORMAT (/1X,A3,' drivers:  ',A9,' indefinite packed matrices',    &
     &        ', partial (Bunch-Kaufman) pivoting')
99011 FORMAT (/1X,A3,' drivers:  ',A9,' indefinite packed matrices',    &
     &        ', "rook" (bounded Bunch-Kaufman) pivoting')
99012 FORMAT (/1X,A3,':  No header available')
!
!     GE matrix types
!
99013 FORMAT (4X,'1. Diagonal',24X,'7. Last n/2 columns zero',/4X,      &
     &        '2. Upper triangular',16X,                                &
     &        '8. Random, CNDNUM = sqrt(0.1/EPS)',/4X,                  &
     &        '3. Lower triangular',16X,'9. Random, CNDNUM = 0.1/EPS',  &
     &        /4X,'4. Random, CNDNUM = 2',13X,                          &
     &        '10. Scaled near underflow',/4X,'5. First column zero',   &
     &        14X,'11. Scaled near overflow',/4X,'6. Last column zero')
!
!     GB matrix types
!
99014 FORMAT (4X,'1. Random, CNDNUM = 2',14X,                           &
     &        '5. Random, CNDNUM = sqrt(0.1/EPS)',/4X,                  &
     &        '2. First column zero',15X,'6. Random, CNDNUM = 0.1/EPS', &
     &        /4X,'3. Last column zero',16X,'7. Scaled near underflow', &
     &        /4X,'4. Last n/2 columns zero',11X,                       &
     &        '8. Scaled near overflow')
!
!     GT matrix types
!
99015 FORMAT (' Matrix types (1-6 have specified condition numbers):',  &
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
99016 FORMAT (' Matrix types (1-6 have specified condition numbers):',  &
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
99017 FORMAT (4X,'1. Diagonal',24X,'6. Random, CNDNUM = sqrt(0.1/EPS)', &
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
!     PB matrix types
!
99018 FORMAT (4X,'1. Random, CNDNUM = 2',14X,                           &
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
!     SSY, SSP, CHE, CHP matrix types
!
99019 FORMAT (4X,'1. Diagonal',24X,'6. Last n/2 rows and columns zero', &
     &        /4X,'2. Random, CNDNUM = 2',14X,                          &
     &        '7. Random, CNDNUM = sqrt(0.1/EPS)',/4X,                  &
     &        '3. First row and column zero',7X,                        &
     &        '8. Random, CNDNUM = 0.1/EPS',/4X,                        &
     &        '4. Last row and column zero',8X,                         &
     &        '9. Scaled near underflow',/4X,                           &
     &        '5. Middle row and column zero',5X,                       &
     &        '10. Scaled near overflow')
!
!     CSY, CSP matrix types
!
99020 FORMAT (4X,'1. Diagonal',24X,'7. Random, CNDNUM = sqrt(0.1/EPS)', &
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
!     Test ratios
!
99021 FORMAT (3X,I2,': norm( L * U - A )  / ( N * norm(A) * EPS )')
99022 FORMAT (3X,I2,': norm( B - A * X )  / ',                          &
     &        '( norm(A) * norm(X) * EPS )')
99023 FORMAT (3X,I2,': norm( X - XACT )   / ',                          &
     &        '( norm(XACT) * CNDNUM * EPS )')
99024 FORMAT (3X,I2,': norm( X - XACT )   / ',                          &
     &        '( norm(XACT) * (error bound) )')
99025 FORMAT (3X,I2,': (backward error)   / EPS')
99026 FORMAT (3X,I2,': RCOND * CNDNUM - 1.0')
99027 FORMAT (3X,I2,': norm( U'' * U - A ) / ( N * norm(A) * EPS )',    &
     &        ', or',/7X,'norm( L * L'' - A ) / ( N * norm(A) * EPS )')
99028 FORMAT (3X,I2,': norm( U*D*U'' - A ) / ( N * norm(A) * EPS )',    &
     &        ', or',/7X,'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )')
99029 FORMAT (3X,I2,': norm( U''*D*U - A ) / ( N * norm(A) * EPS )',    &
     &        ', or',/7X,'norm( L*D*L'' - A ) / ( N * norm(A) * EPS )')
99030 FORMAT (3X,I2,': abs( WORK(1) - RPVGRW ) /',                      &
     &        ' ( max( WORK(1), RPVGRW ) * EPS )')
!
!
!     End of aladhd
!
      END SUBROUTINE aladhd
