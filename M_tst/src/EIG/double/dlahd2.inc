!*==dlahd2.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b dlahd2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dlahd2( IOUNIT, PATH )
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
!> dlahd2 prints header information for the different test paths.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IOUNIT
!> \verbatim
!>          IOUNIT is INTEGER.
!>          On entry, IOUNIT specifies the unit number to which the
!>          header information should be printed.
!> \endverbatim
!>
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3.
!>          On entry, PATH contains the name of the path for which the
!>          header information is to be printed.  Current paths are
!>
!>             DHS, ZHS:  Non-symmetric eigenproblem.
!>             DST, ZST:  Symmetric eigenproblem.
!>             DSG, ZSG:  Symmetric Generalized eigenproblem.
!>             DBD, ZBD:  Singular Value Decomposition (SVD)
!>             DBB, ZBB:  General Banded reduction to bidiagonal form
!>
!>          These paths also are supplied in double precision (replace
!>          leading S by D and leading C by Z in path names).
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
      SUBROUTINE dlahd2(Iounit,Path)
      IMPLICIT NONE
!*--DLAHD269
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
      CHARACTER*2 c2
      INTEGER j
!     ..
!     .. External Functions ..
      LOGICAL LSAME , LSAMEN
      EXTERNAL LSAME , LSAMEN
!     ..
!     .. Executable Statements ..
!
      IF ( Iounit<=0 ) RETURN
      sord = LSAME(Path,'S') .OR. LSAME(Path,'D')
      corz = LSAME(Path,'C') .OR. LSAME(Path,'Z')
      IF ( .NOT.sord .AND. .NOT.corz ) WRITE (Iounit,FMT=99001) Path
      c2 = Path(2:3)
!
      IF ( LSAMEN(2,c2,'HS') ) THEN
         IF ( sord ) THEN
!
!           Real Non-symmetric Eigenvalue Problem:
!
            WRITE (Iounit,FMT=99002) Path
!
!           Matrix types
!
            WRITE (Iounit,FMT=99012)
            WRITE (Iounit,FMT=99013)
            WRITE (Iounit,FMT=99014) 'pairs ' , 'pairs ' , 'prs.' ,     &
     &                               'prs.'
            WRITE (Iounit,FMT=99015)
!
!           Tests performed
!
            WRITE (Iounit,FMT=99016) 'orthogonal' , '''=transpose' ,    &
     &                               ('''',j=1,6)
!
         ELSE
!
!           Complex Non-symmetric Eigenvalue Problem:
!
            WRITE (Iounit,FMT=99003) Path
!
!           Matrix types
!
            WRITE (Iounit,FMT=99012)
            WRITE (Iounit,FMT=99013)
            WRITE (Iounit,FMT=99014) 'e.vals' , 'e.vals' , 'e.vs' ,     &
     &                               'e.vs'
            WRITE (Iounit,FMT=99015)
!
!           Tests performed
!
            WRITE (Iounit,FMT=99016) 'unitary' , '*=conj.transp.' ,     &
     &                               ('*',j=1,6)
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'ST') ) THEN
!
         IF ( sord ) THEN
!
!           Real Symmetric Eigenvalue Problem:
!
            WRITE (Iounit,FMT=99004) Path
!
!           Matrix types
!
            WRITE (Iounit,FMT=99017)
            WRITE (Iounit,FMT=99018)
            WRITE (Iounit,FMT=99019) 'Symmetric'
!
!           Tests performed
!
            WRITE (Iounit,FMT=99032)
!
         ELSE
!
!           Complex Hermitian Eigenvalue Problem:
!
            WRITE (Iounit,FMT=99005) Path
!
!           Matrix types
!
            WRITE (Iounit,FMT=99017)
            WRITE (Iounit,FMT=99018)
            WRITE (Iounit,FMT=99019) 'Hermitian'
!
!           Tests performed
!
            WRITE (Iounit,FMT=99033)
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'SG') ) THEN
!
         IF ( sord ) THEN
!
!           Real Symmetric Generalized Eigenvalue Problem:
!
            WRITE (Iounit,FMT=99008) Path
!
!           Matrix types
!
            WRITE (Iounit,FMT=99020)
            WRITE (Iounit,FMT=99021)
            WRITE (Iounit,FMT=99022) 'Symmetric'
!
!           Tests performed
!
            WRITE (Iounit,FMT=99023)
            WRITE (Iounit,FMT=99024)
!
         ELSE
!
!           Complex Hermitian Generalized Eigenvalue Problem:
!
            WRITE (Iounit,FMT=99009) Path
!
!           Matrix types
!
            WRITE (Iounit,FMT=99020)
            WRITE (Iounit,FMT=99021)
            WRITE (Iounit,FMT=99022) 'Hermitian'
!
!           Tests performed
!
            WRITE (Iounit,FMT=99025)
            WRITE (Iounit,FMT=99026)
!
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'BD') ) THEN
!
         IF ( sord ) THEN
!
!           Real Singular Value Decomposition:
!
            WRITE (Iounit,FMT=99006) Path
!
!           Matrix types
!
            WRITE (Iounit,FMT=99027)
!
!           Tests performed
!
            WRITE (Iounit,FMT=99028) 'orthogonal'
            WRITE (Iounit,FMT=99029)
         ELSE
!
!           Complex Singular Value Decomposition:
!
            WRITE (Iounit,FMT=99007) Path
!
!           Matrix types
!
            WRITE (Iounit,FMT=99027)
!
!           Tests performed
!
            WRITE (Iounit,FMT=99028) 'unitary   '
            WRITE (Iounit,FMT=99029)
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'BB') ) THEN
!
         IF ( sord ) THEN
!
!           Real General Band reduction to bidiagonal form:
!
            WRITE (Iounit,FMT=99010) Path
!
!           Matrix types
!
            WRITE (Iounit,FMT=99030)
!
!           Tests performed
!
            WRITE (Iounit,FMT=99031) 'orthogonal'
         ELSE
!
!           Complex Band reduction to bidiagonal form:
!
            WRITE (Iounit,FMT=99011) Path
!
!           Matrix types
!
            WRITE (Iounit,FMT=99030)
!
!           Tests performed
!
            WRITE (Iounit,FMT=99031) 'unitary   '
         ENDIF
!
      ELSE
!
         WRITE (Iounit,FMT=99001) Path
         RETURN
      ENDIF
!
      RETURN
!
99001 FORMAT (1X,A3,':  no header available')
99002 FORMAT (/1X,A3,' -- Real Non-symmetric eigenvalue problem')
99003 FORMAT (/1X,A3,' -- Complex Non-symmetric eigenvalue problem')
99004 FORMAT (/1X,A3,' -- Real Symmetric eigenvalue problem')
99005 FORMAT (/1X,A3,' -- Complex Hermitian eigenvalue problem')
99006 FORMAT (/1X,A3,' -- Real Singular Value Decomposition')
99007 FORMAT (/1X,A3,' -- Complex Singular Value Decomposition')
99008 FORMAT (/1X,A3,' -- Real Symmetric Generalized eigenvalue ',      &
     &        'problem')
99009 FORMAT (/1X,A3,' -- Complex Hermitian Generalized eigenvalue ',   &
     &        'problem')
99010 FORMAT (/1X,A3,' -- Real Band reduc. to bidiagonal form')
99011 FORMAT (/1X,A3,' -- Complex Band reduc. to bidiagonal form')
!
99012 FORMAT (' Matrix types (see xCHKHS for details): ')
!
99013 FORMAT (/' Special Matrices:',/'  1=Zero matrix.             ',   &
     &        '           ','  5=Diagonal: geometr. spaced entries.',   &
     &        /'  2=Identity matrix.                    ','  6=Diagona',&
     &        'l: clustered entries.',/'  3=Transposed Jordan block.  ',&
     &        '          ','  7=Diagonal: large, evenly spaced.',/'  ', &
     &        '4=Diagonal: evenly spaced entries.    ',                 &
     &        '  8=Diagonal: s','mall, evenly spaced.')
99014 FORMAT (' Dense, Non-Symmetric Matrices:',/'  9=Well-cond., ev',  &
     &        'enly spaced eigenvals.',                                 &
     &        ' 14=Ill-cond., geomet. spaced e','igenals.',             &
     &        /' 10=Well-cond., geom. spaced eigenvals. ',              &
     &        ' 15=Ill-conditioned, clustered e.vals.',/' 11=Well-cond',&
     &        'itioned, clustered e.vals. ',                            &
     &        ' 16=Ill-cond., random comp','lex ',A6,                   &
     &        /' 12=Well-cond., random complex ',A6,'   ',              &
     &        ' 17=Ill-cond., large rand. complx ',A4,/' 13=Ill-condi', &
     &        'tioned, evenly spaced.     ',                            &
     &        ' 18=Ill-cond., small rand.',' complx ',A4)
99015 FORMAT (' 19=Matrix with random O(1) entries.    ',' 21=Matrix ', &
     &        'with small random entries.',/' 20=Matrix with large ran',&
     &        'dom entries.   ')
99016 FORMAT (/' Tests performed:   ','(H is Hessenberg, T is Schur,',  &
     &        ' U and Z are ',A,',',/20X,A,', W is a diagonal matr',    &
     &        'ix of eigenvalues,',/20X,'L and R are the left and rig', &
     &        'ht eigenvector matrices)',/'  1 = | A - U H U',A1,' |',  &
     &        ' / ( |A| n ulp )         ','  2 = | I - U U',A1,' | / ', &
     &        '( n ulp )',/'  3 = | H - Z T Z',A1,' | / ( |H| n ulp ',  &
     &        ')         ','  4 = | I - Z Z',A1,' | / ( n ulp )',       &
     &        /'  5 = | A - UZ T (UZ)',A1,' | / ( |A| n ulp )     ',    &
     &        '  6 = | I - UZ (UZ)',A1,' | / ( n ulp )',/'  7 = | T(',  &
     &        'e.vects.) - T(no e.vects.) | / ( |T| ulp )',/'  8 = | W',&
     &        '(e.vects.) - W(no e.vects.) | / ( |W| ulp )',/'  9 = | ',&
     &        'TR - RW | / ( |T| |R| ulp )     ',                       &
     &        ' 10 = | LT - WL | / (',' |T| |L| ulp )',                 &
     &        /' 11= |HX - XW| / (|H| |X| ulp)  (inv.','it)',           &
     &        ' 12= |YH - WY| / (|H| |Y| ulp)  (inv.it)')
!
!     Symmetric/Hermitian eigenproblem
!
99017 FORMAT (' Matrix types (see xDRVST for details): ')
!
99018 FORMAT (/' Special Matrices:',/'  1=Zero matrix.             ',   &
     &        '           ','  5=Diagonal: clustered entries.',/'  2=', &
     &        'Identity matrix.                    ',                   &
     &        '  6=Diagonal: lar','ge, evenly spaced.',                 &
     &        /'  3=Diagonal: evenly spaced entri','es.    ',           &
     &        '  7=Diagonal: small, evenly spaced.',/'  4=D',           &
     &        'iagonal: geometr. spaced entries.')
99019 FORMAT (' Dense ',A,' Matrices:',/'  8=Evenly spaced eigen',      &
     &        'vals.            ',' 12=Small, evenly spaced eigenvals.',&
     &        /'  9=Geometrically spaced eigenvals.     ',' 13=Matrix ',&
     &        'with random O(1) entries.',/' 10=Clustered eigenvalues.',&
     &        '              ',' 14=Matrix with large random entries.', &
     &        /' 11=Large, evenly spaced eigenvals.     ',' 15=Matrix ',&
     &        'with small random entries.')
!
!     Symmetric/Hermitian Generalized eigenproblem
!
99020 FORMAT (' Matrix types (see xDRVSG for details): ')
!
99021 FORMAT (/' Special Matrices:',/'  1=Zero matrix.             ',   &
     &        '           ','  5=Diagonal: clustered entries.',/'  2=', &
     &        'Identity matrix.                    ',                   &
     &        '  6=Diagonal: lar','ge, evenly spaced.',                 &
     &        /'  3=Diagonal: evenly spaced entri','es.    ',           &
     &        '  7=Diagonal: small, evenly spaced.',/'  4=D',           &
     &        'iagonal: geometr. spaced entries.')
99022 FORMAT (' Dense or Banded ',A,' Matrices: ',                      &
     &        /'  8=Evenly spaced eigenvals.         ',                 &
     &        ' 15=Matrix with small random entries.',                  &
     &        /'  9=Geometrically spaced eigenvals.  ',                 &
     &        ' 16=Evenly spaced eigenvals, KA=1, KB=1.',               &
     &        /' 10=Clustered eigenvalues.           ',                 &
     &        ' 17=Evenly spaced eigenvals, KA=2, KB=1.',               &
     &        /' 11=Large, evenly spaced eigenvals.  ',                 &
     &        ' 18=Evenly spaced eigenvals, KA=2, KB=2.',               &
     &        /' 12=Small, evenly spaced eigenvals.  ',                 &
     &        ' 19=Evenly spaced eigenvals, KA=3, KB=1.',               &
     &        /' 13=Matrix with random O(1) entries. ',                 &
     &        ' 20=Evenly spaced eigenvals, KA=3, KB=2.',               &
     &        /' 14=Matrix with large random entries.',                 &
     &        ' 21=Evenly spaced eigenvals, KA=3, KB=3.')
99023 FORMAT (/' Tests performed:   ',                                  &
     &        /'( For each pair (A,B), where A is of the given type ',  &
     &        /' and B is a random well-conditioned matrix. D is ',     &
     &        /' diagonal, and Z is orthogonal. )',                     &
     &        /' 1 = DSYGV, with ITYPE=1 and UPLO=''U'':',              &
     &        '  | A Z - B Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 2 = DSPGV, with ITYPE=1 and UPLO=''U'':',              &
     &        '  | A Z - B Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 3 = DSBGV, with ITYPE=1 and UPLO=''U'':',              &
     &        '  | A Z - B Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 4 = DSYGV, with ITYPE=1 and UPLO=''L'':',              &
     &        '  | A Z - B Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 5 = DSPGV, with ITYPE=1 and UPLO=''L'':',              &
     &        '  | A Z - B Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 6 = DSBGV, with ITYPE=1 and UPLO=''L'':',              &
     &        '  | A Z - B Z D | / ( |A| |Z| n ulp )     ')
99024 FORMAT (' 7 = DSYGV, with ITYPE=2 and UPLO=''U'':',               &
     &        '  | A B Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 8 = DSPGV, with ITYPE=2 and UPLO=''U'':',              &
     &        '  | A B Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 9 = DSPGV, with ITYPE=2 and UPLO=''L'':',              &
     &        '  | A B Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /'10 = DSPGV, with ITYPE=2 and UPLO=''L'':',              &
     &        '  | A B Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /'11 = DSYGV, with ITYPE=3 and UPLO=''U'':',              &
     &        '  | B A Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /'12 = DSPGV, with ITYPE=3 and UPLO=''U'':',              &
     &        '  | B A Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /'13 = DSYGV, with ITYPE=3 and UPLO=''L'':',              &
     &        '  | B A Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /'14 = DSPGV, with ITYPE=3 and UPLO=''L'':',              &
     &        '  | B A Z - Z D | / ( |A| |Z| n ulp )     ')
99025 FORMAT (/' Tests performed:   ',                                  &
     &        /'( For each pair (A,B), where A is of the given type ',  &
     &        /' and B is a random well-conditioned matrix. D is ',     &
     &        /' diagonal, and Z is unitary. )',                        &
     &        /' 1 = ZHEGV, with ITYPE=1 and UPLO=''U'':',              &
     &        '  | A Z - B Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 2 = ZHPGV, with ITYPE=1 and UPLO=''U'':',              &
     &        '  | A Z - B Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 3 = ZHBGV, with ITYPE=1 and UPLO=''U'':',              &
     &        '  | A Z - B Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 4 = ZHEGV, with ITYPE=1 and UPLO=''L'':',              &
     &        '  | A Z - B Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 5 = ZHPGV, with ITYPE=1 and UPLO=''L'':',              &
     &        '  | A Z - B Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 6 = ZHBGV, with ITYPE=1 and UPLO=''L'':',              &
     &        '  | A Z - B Z D | / ( |A| |Z| n ulp )     ')
99026 FORMAT (' 7 = ZHEGV, with ITYPE=2 and UPLO=''U'':',               &
     &        '  | A B Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 8 = ZHPGV, with ITYPE=2 and UPLO=''U'':',              &
     &        '  | A B Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /' 9 = ZHPGV, with ITYPE=2 and UPLO=''L'':',              &
     &        '  | A B Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /'10 = ZHPGV, with ITYPE=2 and UPLO=''L'':',              &
     &        '  | A B Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /'11 = ZHEGV, with ITYPE=3 and UPLO=''U'':',              &
     &        '  | B A Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /'12 = ZHPGV, with ITYPE=3 and UPLO=''U'':',              &
     &        '  | B A Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /'13 = ZHEGV, with ITYPE=3 and UPLO=''L'':',              &
     &        '  | B A Z - Z D | / ( |A| |Z| n ulp )     ',             &
     &        /'14 = ZHPGV, with ITYPE=3 and UPLO=''L'':',              &
     &        '  | B A Z - Z D | / ( |A| |Z| n ulp )     ')
!
!     Singular Value Decomposition
!
99027 FORMAT (' Matrix types (see xCHKBD for details):',                &
     &        /' Diagonal matrices:',/'   1: Zero',28X,                 &
     &        ' 5: Clustered entries',/'   2: Identity',24X,            &
     &        ' 6: Large, evenly spaced entries',                       &
     &        /'   3: Evenly spaced entries',11X,                       &
     &        ' 7: Small, evenly spaced entries',                       &
     &        /'   4: Geometrically spaced entries',                    &
     &        /' General matrices:',/'   8: Evenly spaced sing. vals.', &
     &        7X,'12: Small, evenly spaced sing vals',                  &
     &        /'   9: Geometrically spaced sing vals  ',                &
     &        '13: Random, O(1) entries',/'  10: Clustered sing. vals.',&
     &        11X,'14: Random, scaled near overflow',                   &
     &        /'  11: Large, evenly spaced sing vals  ',                &
     &        '15: Random, scaled near underflow')
!
99028 FORMAT (/' Test ratios:  ',                                       &
     &        '(B: bidiagonal, S: diagonal, Q, P, U, and V: ',A10,/16X, &
     &        'X: m x nrhs, Y = Q'' X, and Z = U'' Y)')
99029 FORMAT ('   1: norm( A - Q B P'' ) / ( norm(A) max(m,n) ulp )',   &
     &        /'   2: norm( I - Q'' Q )   / ( m ulp )',                 &
     &        /'   3: norm( I - P'' P )   / ( n ulp )',                 &
     &        /'   4: norm( B - U S V'' ) / ( norm(B) min(m,n) ulp )',  &
     &        /'   5: norm( Y - U Z )    / ',                           &
     &        '( norm(Z) max(min(m,n),k) ulp )',                        &
     &        /'   6: norm( I - U'' U )   / ( min(m,n) ulp )',          &
     &        /'   7: norm( I - V'' V )   / ( min(m,n) ulp )',          &
     &        /'   8: Test ordering of S  (0 if nondecreasing, 1/ulp ', &
     &        ' otherwise)',                                            &
     &        /'   9: norm( S - S1 )     / ( norm(S) ulp ),',           &
     &        ' where S1 is computed',/43X,                             &
     &        ' without computing U and V''',                           &
     &        /'  10: Sturm sequence test ',                            &
     &        '(0 if sing. vals of B within THRESH of S)',              &
     &        /'  11: norm( A - (QU) S (V'' P'') ) / ',                 &
     &        '( norm(A) max(m,n) ulp )',                               &
     &        /'  12: norm( X - (QU) Z )         / ( |X| max(M,k) ulp )'&
     &        ,/'  13: norm( I - (QU)''(QU) )      / ( M ulp )',        &
     &        /'  14: norm( I - (V'' P'') (P V) )  / ( N ulp )',        &
     &        /'  15: norm( B - U S V'' ) / ( norm(B) min(m,n) ulp )',  &
     &        /'  16: norm( I - U'' U )   / ( min(m,n) ulp )',          &
     &        /'  17: norm( I - V'' V )   / ( min(m,n) ulp )',          &
     &        /'  18: Test ordering of S  (0 if nondecreasing, 1/ulp ', &
     &        ' otherwise)',                                            &
     &        /'  19: norm( S - S1 )     / ( norm(S) ulp ),',           &
     &        ' where S1 is computed',/43X,                             &
     &        ' without computing U and V''',                           &
     &        /'  20: norm( B - U S V'' )  / ( norm(B) min(m,n) ulp )', &
     &        '  DBDSVX(V,A)',                                          &
     &        /'  21: norm( I - U'' U )    / ( min(m,n) ulp )',         &
     &        /'  22: norm( I - V'' V )    / ( min(m,n) ulp )',         &
     &        /'  23: Test ordering of S  (0 if nondecreasing, 1/ulp ', &
     &        ' otherwise)',                                            &
     &        /'  24: norm( S - S1 )      / ( norm(S) ulp ),',          &
     &        ' where S1 is computed',/44X,                             &
     &        ' without computing U and V''',                           &
     &        /'  25: norm( S - U'' B V ) / ( norm(B) n ulp )',         &
     &        '  DBDSVX(V,I)',                                          &
     &        /'  26: norm( I - U'' U )    / ( min(m,n) ulp )',         &
     &        /'  27: norm( I - V'' V )    / ( min(m,n) ulp )',         &
     &        /'  28: Test ordering of S  (0 if nondecreasing, 1/ulp ', &
     &        ' otherwise)',                                            &
     &        /'  29: norm( S - S1 )      / ( norm(S) ulp ),',          &
     &        ' where S1 is computed',/44X,                             &
     &        ' without computing U and V''',                           &
     &        /'  30: norm( S - U'' B V ) / ( norm(B) n ulp )',         &
     &        '  DBDSVX(V,V)',                                          &
     &        /'  31: norm( I - U'' U )    / ( min(m,n) ulp )',         &
     &        /'  32: norm( I - V'' V )    / ( min(m,n) ulp )',         &
     &        /'  33: Test ordering of S  (0 if nondecreasing, 1/ulp ', &
     &        ' otherwise)',                                            &
     &        /'  34: norm( S - S1 )      / ( norm(S) ulp ),',          &
     &        ' where S1 is computed',/44X,                             &
     &        ' without computing U and V''')
!
!     Band reduction to bidiagonal form
!
99030 FORMAT (' Matrix types (see xCHKBB for details):',                &
     &        /' Diagonal matrices:',/'   1: Zero',28X,                 &
     &        ' 5: Clustered entries',/'   2: Identity',24X,            &
     &        ' 6: Large, evenly spaced entries',                       &
     &        /'   3: Evenly spaced entries',11X,                       &
     &        ' 7: Small, evenly spaced entries',                       &
     &        /'   4: Geometrically spaced entries',                    &
     &        /' General matrices:',/'   8: Evenly spaced sing. vals.', &
     &        7X,'12: Small, evenly spaced sing vals',                  &
     &        /'   9: Geometrically spaced sing vals  ',                &
     &        '13: Random, O(1) entries',/'  10: Clustered sing. vals.',&
     &        11X,'14: Random, scaled near overflow',                   &
     &        /'  11: Large, evenly spaced sing vals  ',                &
     &        '15: Random, scaled near underflow')
!
99031 FORMAT (/' Test ratios:  ','(B: upper bidiagonal, Q and P: ',A10, &
     &        /16X,'C: m x nrhs, PT = P'', Y = Q'' C)',                 &
     &        /' 1: norm( A - Q B PT ) / ( norm(A) max(m,n) ulp )',     &
     &        /' 2: norm( I - Q'' Q )   / ( m ulp )',                   &
     &        /' 3: norm( I - PT PT'' )   / ( n ulp )',                 &
     &        /' 4: norm( Y - Q'' C )   / ( norm(Y) max(m,nrhs) ulp )')
99032 FORMAT (/' Tests performed:  See sdrvst.f')
99033 FORMAT (/' Tests performed:  See cdrvst.f')
!
!     End of dlahd2
!
      END SUBROUTINE dlahd2
