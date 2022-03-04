!*==sdrves.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SDRVES
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SDRVES( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          NOUNIT, A, LDA, H, HT, WR, WI, WRT, WIT, VS,
!                          LDVS, RESULT, WORK, NWORK, IWORK, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDVS, NOUNIT, NSIZES, NTYPES, NWORK
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * ), DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
!       REAL               A( LDA, * ), H( LDA, * ), HT( LDA, * ),
!      $                   RESULT( 13 ), VS( LDVS, * ), WI( * ), WIT( * ),
!      $                   WORK( * ), WR( * ), WRT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SDRVES checks the nonsymmetric eigenvalue (Schur form) problem
!>    driver SGEES.
!>
!>    When SDRVES is called, a number of matrix "sizes" ("n's") and a
!>    number of matrix "types" are specified.  For each size ("n")
!>    and each type of matrix, one matrix will be generated and used
!>    to test the nonsymmetric eigenroutines.  For each matrix, 13
!>    tests will be performed:
!>
!>    (1)     0 if T is in Schur form, 1/ulp otherwise
!>           (no sorting of eigenvalues)
!>
!>    (2)     | A - VS T VS' | / ( n |A| ulp )
!>
!>      Here VS is the matrix of Schur eigenvectors, and T is in Schur
!>      form  (no sorting of eigenvalues).
!>
!>    (3)     | I - VS VS' | / ( n ulp ) (no sorting of eigenvalues).
!>
!>    (4)     0     if WR+sqrt(-1)*WI are eigenvalues of T
!>            1/ulp otherwise
!>            (no sorting of eigenvalues)
!>
!>    (5)     0     if T(with VS) = T(without VS),
!>            1/ulp otherwise
!>            (no sorting of eigenvalues)
!>
!>    (6)     0     if eigenvalues(with VS) = eigenvalues(without VS),
!>            1/ulp otherwise
!>            (no sorting of eigenvalues)
!>
!>    (7)     0 if T is in Schur form, 1/ulp otherwise
!>            (with sorting of eigenvalues)
!>
!>    (8)     | A - VS T VS' | / ( n |A| ulp )
!>
!>      Here VS is the matrix of Schur eigenvectors, and T is in Schur
!>      form  (with sorting of eigenvalues).
!>
!>    (9)     | I - VS VS' | / ( n ulp ) (with sorting of eigenvalues).
!>
!>    (10)    0     if WR+sqrt(-1)*WI are eigenvalues of T
!>            1/ulp otherwise
!>            (with sorting of eigenvalues)
!>
!>    (11)    0     if T(with VS) = T(without VS),
!>            1/ulp otherwise
!>            (with sorting of eigenvalues)
!>
!>    (12)    0     if eigenvalues(with VS) = eigenvalues(without VS),
!>            1/ulp otherwise
!>            (with sorting of eigenvalues)
!>
!>    (13)    if sorting worked and SDIM is the number of
!>            eigenvalues which were SELECTed
!>
!>    The "sizes" are specified by an array NN(1:NSIZES); the value of
!>    each element NN(j) specifies one size.
!>    The "types" are specified by a logical array DOTYPE( 1:NTYPES );
!>    if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!>    Currently, the list of possible types is:
!>
!>    (1)  The zero matrix.
!>    (2)  The identity matrix.
!>    (3)  A (transposed) Jordan block, with 1's on the diagonal.
!>
!>    (4)  A diagonal matrix with evenly spaced entries
!>         1, ..., ULP  and random signs.
!>         (ULP = (first number larger than 1) - 1 )
!>    (5)  A diagonal matrix with geometrically spaced entries
!>         1, ..., ULP  and random signs.
!>    (6)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
!>         and random signs.
!>
!>    (7)  Same as (4), but multiplied by a constant near
!>         the overflow threshold
!>    (8)  Same as (4), but multiplied by a constant near
!>         the underflow threshold
!>
!>    (9)  A matrix of the form  U' T U, where U is orthogonal and
!>         T has evenly spaced entries 1, ..., ULP with random signs
!>         on the diagonal and random O(1) entries in the upper
!>         triangle.
!>
!>    (10) A matrix of the form  U' T U, where U is orthogonal and
!>         T has geometrically spaced entries 1, ..., ULP with random
!>         signs on the diagonal and random O(1) entries in the upper
!>         triangle.
!>
!>    (11) A matrix of the form  U' T U, where U is orthogonal and
!>         T has "clustered" entries 1, ULP,..., ULP with random
!>         signs on the diagonal and random O(1) entries in the upper
!>         triangle.
!>
!>    (12) A matrix of the form  U' T U, where U is orthogonal and
!>         T has real or complex conjugate paired eigenvalues randomly
!>         chosen from ( ULP, 1 ) and random O(1) entries in the upper
!>         triangle.
!>
!>    (13) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has evenly spaced entries 1, ..., ULP
!>         with random signs on the diagonal and random O(1) entries
!>         in the upper triangle.
!>
!>    (14) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has geometrically spaced entries
!>         1, ..., ULP with random signs on the diagonal and random
!>         O(1) entries in the upper triangle.
!>
!>    (15) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has "clustered" entries 1, ULP,..., ULP
!>         with random signs on the diagonal and random O(1) entries
!>         in the upper triangle.
!>
!>    (16) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has real or complex conjugate paired
!>         eigenvalues randomly chosen from ( ULP, 1 ) and random
!>         O(1) entries in the upper triangle.
!>
!>    (17) Same as (16), but multiplied by a constant
!>         near the overflow threshold
!>    (18) Same as (16), but multiplied by a constant
!>         near the underflow threshold
!>
!>    (19) Nonsymmetric matrix with random entries chosen from (-1,1).
!>         If N is at least 4, all entries in first two rows and last
!>         row, and first column and last two columns are zero.
!>    (20) Same as (19), but multiplied by a constant
!>         near the overflow threshold
!>    (21) Same as (19), but multiplied by a constant
!>         near the underflow threshold
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          SDRVES does nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER array, dimension (NSIZES)
!>          An array containing the sizes to be used for the matrices.
!>          Zero values will be skipped.  The values must be at least
!>          zero.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, SDRVES
!>          does nothing.  It must be at least zero.  If it is MAXTYP+1
!>          and NSIZES is 1, then an additional type, MAXTYP+1 is
!>          defined, which is to use whatever matrix is in A.  This
!>          is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>          DOTYPE(MAXTYP+1) is .TRUE. .
!> \endverbatim
!>
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          If DOTYPE(j) is .TRUE., then for each size in NN a
!>          matrix of that size and of type j will be generated.
!>          If NTYPES is smaller than the maximum number of types
!>          defined (PARAMETER MAXTYP), then types NTYPES+1 through
!>          MAXTYP will not be generated.  If NTYPES is larger
!>          than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
!>          will be ignored.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator. The array elements should be between 0 and 4095;
!>          if not they will be reduced mod 4096.  Also, ISEED(4) must
!>          be odd.  The random number generator uses a linear
!>          congruential sequence limited to small integers, and so
!>          should produce machine independent random numbers. The
!>          values of ISEED are changed on exit, and can be used in the
!>          next call to SDRVES to continue the same random number
!>          sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error
!>          is scaled to be O(1), so THRESH should be a reasonably
!>          small multiple of 1, e.g., 10 or 100.  In particular,
!>          it should not depend on the precision (single vs. double)
!>          or the size of the matrix.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] NOUNIT
!> \verbatim
!>          NOUNIT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns INFO not equal to 0.)
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (LDA, max(NN))
!>          Used to hold the matrix whose eigenvalues are to be
!>          computed.  On exit, A contains the last matrix actually used.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A, and H. LDA must be at
!>          least 1 and at least max(NN).
!> \endverbatim
!>
!> \param[out] H
!> \verbatim
!>          H is REAL array, dimension (LDA, max(NN))
!>          Another copy of the test matrix A, modified by SGEES.
!> \endverbatim
!>
!> \param[out] HT
!> \verbatim
!>          HT is REAL array, dimension (LDA, max(NN))
!>          Yet another copy of the test matrix A, modified by SGEES.
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is REAL array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is REAL array, dimension (max(NN))
!>
!>          The real and imaginary parts of the eigenvalues of A.
!>          On exit, WR + WI*i are the eigenvalues of the matrix in A.
!> \endverbatim
!>
!> \param[out] WRT
!> \verbatim
!>          WRT is REAL array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] WIT
!> \verbatim
!>          WIT is REAL array, dimension (max(NN))
!>
!>          Like WR, WI, these arrays contain the eigenvalues of A,
!>          but those computed when SGEES only computes a partial
!>          eigendecomposition, i.e. not Schur vectors
!> \endverbatim
!>
!> \param[out] VS
!> \verbatim
!>          VS is REAL array, dimension (LDVS, max(NN))
!>          VS holds the computed Schur vectors.
!> \endverbatim
!>
!> \param[in] LDVS
!> \verbatim
!>          LDVS is INTEGER
!>          Leading dimension of VS. Must be at least max(1,max(NN)).
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (13)
!>          The values computed by the 13 tests described above.
!>          The values are currently limited to 1/ulp, to avoid overflow.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (NWORK)
!> \endverbatim
!>
!> \param[in] NWORK
!> \verbatim
!>          NWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          5*NN(j)+2*NN(j)**2 for all j.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] BWORK
!> \verbatim
!>          BWORK is LOGICAL array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0, then everything ran OK.
!>           -1: NSIZES < 0
!>           -2: Some NN(j) < 0
!>           -3: NTYPES < 0
!>           -6: THRESH < 0
!>           -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
!>          -17: LDVS < 1 or LDVS < NMAX, where NMAX is max( NN(j) ).
!>          -20: NWORK too small.
!>          If  SLATMR, SLATMS, SLATME or SGEES returns an error code,
!>              the absolute value of it is returned.
!>
!>-----------------------------------------------------------------------
!>
!>     Some Local Variables and Parameters:
!>     ---- ----- --------- --- ----------
!>
!>     ZERO, ONE       Real 0 and 1.
!>     MAXTYP          The number of types defined.
!>     NMAX            Largest value in NN.
!>     NERRS           The number of tests which have exceeded THRESH
!>     COND, CONDS,
!>     IMODE           Values to be passed to the matrix generators.
!>     ANORM           Norm of A; passed to matrix generators.
!>
!>     OVFL, UNFL      Overflow and underflow thresholds.
!>     ULP, ULPINV     Finest relative precision and its inverse.
!>     RTULP, RTULPI   Square roots of the previous 4 values.
!>
!>             The following four arrays decode JTYPE:
!>     KTYPE(j)        The general type (1-10) for type "j".
!>     KMODE(j)        The MODE value to be passed to the matrix
!>                     generator for type "j".
!>     KMAGN(j)        The order of magnitude ( O(1),
!>                     O(overflow^(1/2) ), O(underflow^(1/2) )
!>     KCONDS(j)       Selectw whether CONDS is to be 1 or
!>                     1/sqrt(ulp).  (0 means irrelevant.)
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE SDRVES(Nsizes,Nn,Ntypes,Dotype,Iseed,Thresh,Nounit,A,  &
     &                  Lda,H,Ht,Wr,Wi,Wrt,Wit,Vs,Ldvs,Result,Work,     &
     &                  Nwork,Iwork,Bwork,Info)
      IMPLICIT NONE
!*--SDRVES392
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldvs , Nounit , Nsizes , Ntypes , Nwork
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Bwork(*) , Dotype(*)
      INTEGER Iseed(4) , Iwork(*) , Nn(*)
      REAL A(Lda,*) , H(Lda,*) , Ht(Lda,*) , Result(13) , Vs(Ldvs,*) ,  &
     &     Wi(*) , Wit(*) , Work(*) , Wr(*) , Wrt(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
      INTEGER MAXTYP
      PARAMETER (MAXTYP=21)
!     ..
!     .. Local Scalars ..
      LOGICAL badnn
      CHARACTER sort
      CHARACTER*3 path
      INTEGER i , iinfo , imode , isort , itype , iwk , j , jcol ,      &
     &        jsize , jtype , knteig , lwork , mtypes , n , nerrs ,     &
     &        nfail , nmax , nnwork , ntest , ntestf , ntestt , rsub ,  &
     &        sdim
      REAL anorm , cond , conds , ovfl , rtulp , rtulpi , tmp , ulp ,   &
     &     ulpinv , unfl
!     ..
!     .. Local Arrays ..
      CHARACTER adumma(1)
      INTEGER idumma(1) , ioldsd(4) , kconds(MAXTYP) , kmagn(MAXTYP) ,  &
     &        kmode(MAXTYP) , ktype(MAXTYP)
      REAL res(2)
!     ..
!     .. Arrays in Common ..
      LOGICAL SELval(20)
      REAL SELwi(20) , SELwr(20)
!     ..
!     .. Scalars in Common ..
      INTEGER SELdim , SELopt
!     ..
!     .. Common blocks ..
      COMMON /SSLCT / SELopt , SELdim , SELval , SELwr , SELwi
!     ..
!     .. External Functions ..
      LOGICAL SSLECT
      REAL SLAMCH
      EXTERNAL SSLECT , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEES , SHST01 , SLABAD , SLACPY , SLASUM , SLATME ,     &
     &         SLATMR , SLATMS , SLASET , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SIGN , SQRT
!     ..
!     .. Data statements ..
      DATA ktype/1 , 2 , 3 , 5*4 , 4*6 , 6*6 , 3*9/
      DATA kmagn/3*1 , 1 , 1 , 1 , 2 , 3 , 4*1 , 1 , 1 , 1 , 1 , 2 , 3 ,&
     &     1 , 2 , 3/
      DATA kmode/3*0 , 4 , 3 , 1 , 4 , 4 , 4 , 3 , 1 , 5 , 4 , 3 , 1 ,  &
     &     5 , 5 , 5 , 4 , 3 , 1/
      DATA kconds/3*0 , 5*0 , 4*1 , 6*2 , 3*0/
!     ..
!     .. Executable Statements ..
!
      path(1:1) = 'Single precision'
      path(2:3) = 'ES'
!
!     Check for errors
!
      ntestt = 0
      ntestf = 0
      Info = 0
      SELopt = 0
!
!     Important constants
!
      badnn = .FALSE.
      nmax = 0
      DO j = 1 , Nsizes
         nmax = MAX(nmax,Nn(j))
         IF ( Nn(j)<0 ) badnn = .TRUE.
      ENDDO
!
!     Check for errors
!
      IF ( Nsizes<0 ) THEN
         Info = -1
      ELSEIF ( badnn ) THEN
         Info = -2
      ELSEIF ( Ntypes<0 ) THEN
         Info = -3
      ELSEIF ( Thresh<ZERO ) THEN
         Info = -6
      ELSEIF ( Nounit<=0 ) THEN
         Info = -7
      ELSEIF ( Lda<1 .OR. Lda<nmax ) THEN
         Info = -9
      ELSEIF ( Ldvs<1 .OR. Ldvs<nmax ) THEN
         Info = -17
      ELSEIF ( 5*nmax+2*nmax**2>Nwork ) THEN
         Info = -20
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SDRVES',-Info)
         RETURN
      ENDIF
!
!     Quick return if nothing to do
!
      IF ( Nsizes==0 .OR. Ntypes==0 ) RETURN
!
!     More Important constants
!
      unfl = SLAMCH('Safe minimum')
      ovfl = ONE/unfl
      CALL SLABAD(unfl,ovfl)
      ulp = SLAMCH('Precision')
      ulpinv = ONE/ulp
      rtulp = SQRT(ulp)
      rtulpi = ONE/rtulp
!
!     Loop over sizes, types
!
      nerrs = 0
!
      DO jsize = 1 , Nsizes
         n = Nn(jsize)
         mtypes = MAXTYP
         IF ( Nsizes==1 .AND. Ntypes==MAXTYP+1 ) mtypes = mtypes + 1
!
         DO jtype = 1 , mtypes
            IF ( Dotype(jtype) ) THEN
!
!           Save ISEED in case of an error.
!
               DO j = 1 , 4
                  ioldsd(j) = Iseed(j)
               ENDDO
!
!           Compute "A"
!
!           Control parameters:
!
!           KMAGN  KCONDS  KMODE        KTYPE
!       =1  O(1)   1       clustered 1  zero
!       =2  large  large   clustered 2  identity
!       =3  small          exponential  Jordan
!       =4                 arithmetic   diagonal, (w/ eigenvalues)
!       =5                 random log   symmetric, w/ eigenvalues
!       =6                 random       general, w/ eigenvalues
!       =7                              random diagonal
!       =8                              random symmetric
!       =9                              random general
!       =10                             random triangular
!
               IF ( mtypes<=MAXTYP ) THEN
!
                  itype = ktype(jtype)
                  imode = kmode(jtype)
!
!           Compute norm
!
                  IF ( kmagn(jtype)==2 ) THEN
!
                     anorm = ovfl*ulp
                  ELSEIF ( kmagn(jtype)==3 ) THEN
!
                     anorm = unfl*ulpinv
                  ELSE
!
                     anorm = ONE
                  ENDIF
!
!
                  CALL SLASET('Full',Lda,n,ZERO,ZERO,A,Lda)
                  iinfo = 0
                  cond = ulpinv
!
!           Special Matrices -- Identity & Jordan block
!
!              Zero
!
                  IF ( itype==1 ) THEN
                     iinfo = 0
!
                  ELSEIF ( itype==2 ) THEN
!
!              Identity
!
                     DO jcol = 1 , n
                        A(jcol,jcol) = anorm
                     ENDDO
!
                  ELSEIF ( itype==3 ) THEN
!
!              Jordan Block
!
                     DO jcol = 1 , n
                        A(jcol,jcol) = anorm
                        IF ( jcol>1 ) A(jcol,jcol-1) = ONE
                     ENDDO
!
                  ELSEIF ( itype==4 ) THEN
!
!              Diagonal Matrix, [Eigen]values Specified
!
                     CALL SLATMS(n,n,'S',Iseed,'S',Work,imode,cond,     &
     &                           anorm,0,0,'N',A,Lda,Work(n+1),iinfo)
!
                  ELSEIF ( itype==5 ) THEN
!
!              Symmetric, eigenvalues specified
!
                     CALL SLATMS(n,n,'S',Iseed,'S',Work,imode,cond,     &
     &                           anorm,n,n,'N',A,Lda,Work(n+1),iinfo)
!
                  ELSEIF ( itype==6 ) THEN
!
!              General, eigenvalues specified
!
                     IF ( kconds(jtype)==1 ) THEN
                        conds = ONE
                     ELSEIF ( kconds(jtype)==2 ) THEN
                        conds = rtulpi
                     ELSE
                        conds = ZERO
                     ENDIF
!
                     adumma(1) = ' '
                     CALL SLATME(n,'S',Iseed,Work,imode,cond,ONE,adumma,&
     &                           'T','T','T',Work(n+1),4,conds,n,n,     &
     &                           anorm,A,Lda,Work(2*n+1),iinfo)
!
                  ELSEIF ( itype==7 ) THEN
!
!              Diagonal, random eigenvalues
!
                     CALL SLATMR(n,n,'S',Iseed,'S',Work,6,ONE,ONE,'T',  &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,0,0,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==8 ) THEN
!
!              Symmetric, random eigenvalues
!
                     CALL SLATMR(n,n,'S',Iseed,'S',Work,6,ONE,ONE,'T',  &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,n,n,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==9 ) THEN
!
!              General, random eigenvalues
!
                     CALL SLATMR(n,n,'S',Iseed,'N',Work,6,ONE,ONE,'T',  &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,n,n,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
                     IF ( n>=4 ) THEN
                        CALL SLASET('Full',2,n,ZERO,ZERO,A,Lda)
                        CALL SLASET('Full',n-3,1,ZERO,ZERO,A(3,1),Lda)
                        CALL SLASET('Full',n-3,2,ZERO,ZERO,A(3,n-1),Lda)
                        CALL SLASET('Full',1,n,ZERO,ZERO,A(n,1),Lda)
                     ENDIF
!
                  ELSEIF ( itype==10 ) THEN
!
!              Triangular, random eigenvalues
!
                     CALL SLATMR(n,n,'S',Iseed,'N',Work,6,ONE,ONE,'T',  &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,n,0,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSE
!
                     iinfo = 1
                  ENDIF
!
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99008) 'Generator' , iinfo , n , &
     &                      jtype , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
               ENDIF
!
!
!           Test for minimal and generous workspace
!
               DO iwk = 1 , 2
                  IF ( iwk==1 ) THEN
                     nnwork = 3*n
                  ELSE
                     nnwork = 5*n + 2*n**2
                  ENDIF
                  nnwork = MAX(nnwork,1)
!
!              Initialize RESULT
!
                  DO j = 1 , 13
                     Result(j) = -ONE
                  ENDDO
!
!              Test with and without sorting of eigenvalues
!
                  DO isort = 0 , 1
                     IF ( isort==0 ) THEN
                        sort = 'N'
                        rsub = 0
                     ELSE
                        sort = 'S'
                        rsub = 6
                     ENDIF
!
!                 Compute Schur form and Schur vectors, and test them
!
                     CALL SLACPY('F',n,n,A,Lda,H,Lda)
                     CALL SGEES('V',sort,SSLECT,n,H,Lda,sdim,Wr,Wi,Vs,  &
     &                          Ldvs,Work,nnwork,Bwork,iinfo)
                     IF ( iinfo/=0 .AND. iinfo/=n+2 ) THEN
                        Result(1+rsub) = ulpinv
                        WRITE (Nounit,FMT=99008) 'SGEES1' , iinfo , n , &
     &                         jtype , ioldsd
                        Info = ABS(iinfo)
                        EXIT
                     ENDIF
!
!                 Do Test (1) or Test (7)
!
                     Result(1+rsub) = ZERO
                     DO j = 1 , n - 2
                        DO i = j + 2 , n
                           IF ( H(i,j)/=ZERO ) Result(1+rsub) = ulpinv
                        ENDDO
                     ENDDO
                     DO i = 1 , n - 2
                        IF ( H(i+1,i)/=ZERO .AND. H(i+2,i+1)/=ZERO )    &
     &                       Result(1+rsub) = ulpinv
                     ENDDO
                     DO i = 1 , n - 1
                        IF ( H(i+1,i)/=ZERO ) THEN
                           IF ( H(i,i)/=H(i+1,i+1) .OR. H(i,i+1)        &
     &                          ==ZERO .OR. SIGN(ONE,H(i+1,i))          &
     &                          ==SIGN(ONE,H(i,i+1)) ) Result(1+rsub)   &
     &                          = ulpinv
                        ENDIF
                     ENDDO
!
!                 Do Tests (2) and (3) or Tests (8) and (9)
!
                     lwork = MAX(1,2*n*n)
                     CALL SHST01(n,1,n,A,Lda,H,Lda,Vs,Ldvs,Work,lwork,  &
     &                           res)
                     Result(2+rsub) = res(1)
                     Result(3+rsub) = res(2)
!
!                 Do Test (4) or Test (10)
!
                     Result(4+rsub) = ZERO
                     DO i = 1 , n
                        IF ( H(i,i)/=Wr(i) ) Result(4+rsub) = ulpinv
                     ENDDO
                     IF ( n>1 ) THEN
                        IF ( H(2,1)==ZERO .AND. Wi(1)/=ZERO )           &
     &                       Result(4+rsub) = ulpinv
                        IF ( H(n,n-1)==ZERO .AND. Wi(n)/=ZERO )         &
     &                       Result(4+rsub) = ulpinv
                     ENDIF
                     DO i = 1 , n - 1
                        IF ( H(i+1,i)/=ZERO ) THEN
                           tmp = SQRT(ABS(H(i+1,i)))*SQRT(ABS(H(i,i+1)))
                           Result(4+rsub)                               &
     &                        = MAX(Result(4+rsub),ABS(Wi(i)-tmp)       &
     &                        /MAX(ulp*tmp,unfl))
                           Result(4+rsub)                               &
     &                        = MAX(Result(4+rsub),ABS(Wi(i+1)+tmp)     &
     &                        /MAX(ulp*tmp,unfl))
                        ELSEIF ( i>1 ) THEN
                           IF ( H(i+1,i)==ZERO .AND. H(i,i-1)           &
     &                          ==ZERO .AND. Wi(i)/=ZERO )              &
     &                          Result(4+rsub) = ulpinv
                        ENDIF
                     ENDDO
!
!                 Do Test (5) or Test (11)
!
                     CALL SLACPY('F',n,n,A,Lda,Ht,Lda)
                     CALL SGEES('N',sort,SSLECT,n,Ht,Lda,sdim,Wrt,Wit,  &
     &                          Vs,Ldvs,Work,nnwork,Bwork,iinfo)
                     IF ( iinfo/=0 .AND. iinfo/=n+2 ) THEN
                        Result(5+rsub) = ulpinv
                        WRITE (Nounit,FMT=99008) 'SGEES2' , iinfo , n , &
     &                         jtype , ioldsd
                        Info = ABS(iinfo)
                        EXIT
                     ENDIF
!
                     Result(5+rsub) = ZERO
                     DO j = 1 , n
                        DO i = 1 , n
                           IF ( H(i,j)/=Ht(i,j) ) Result(5+rsub)        &
     &                          = ulpinv
                        ENDDO
                     ENDDO
!
!                 Do Test (6) or Test (12)
!
                     Result(6+rsub) = ZERO
                     DO i = 1 , n
                        IF ( Wr(i)/=Wrt(i) .OR. Wi(i)/=Wit(i) )         &
     &                       Result(6+rsub) = ulpinv
                     ENDDO
!
!                 Do Test (13)
!
                     IF ( isort==1 ) THEN
                        Result(13) = ZERO
                        knteig = 0
                        DO i = 1 , n
                           IF ( SSLECT(Wr(i),Wi(i)) .OR.                &
     &                          SSLECT(Wr(i),-Wi(i)) ) knteig = knteig +&
     &                          1
                           IF ( i<n ) THEN
                              IF ( (SSLECT(Wr(i+1),Wi(i+1)) .OR. SSLECT(&
     &                             Wr(i+1),-Wi(i+1))) .AND.             &
     &                             (.NOT.(SSLECT(Wr(i),Wi(i)) .OR.      &
     &                             SSLECT(Wr(i),-Wi(i)))) .AND.         &
     &                             iinfo/=n+2 ) Result(13) = ulpinv
                           ENDIF
                        ENDDO
                        IF ( sdim/=knteig ) Result(13) = ulpinv
                     ENDIF
!
                  ENDDO
!
!              End of Loop -- Check for RESULT(j) > THRESH
!
!
                  ntest = 0
                  nfail = 0
                  DO j = 1 , 13
                     IF ( Result(j)>=ZERO ) ntest = ntest + 1
                     IF ( Result(j)>=Thresh ) nfail = nfail + 1
                  ENDDO
!
                  IF ( nfail>0 ) ntestf = ntestf + 1
                  IF ( ntestf==1 ) THEN
                     WRITE (Nounit,FMT=99001) path
                     WRITE (Nounit,FMT=99002)
                     WRITE (Nounit,FMT=99003)
                     WRITE (Nounit,FMT=99004)
                     WRITE (Nounit,FMT=99005) Thresh
                     WRITE (Nounit,FMT=99006)
                     ntestf = 2
                  ENDIF
!
                  DO j = 1 , 13
                     IF ( Result(j)>=Thresh ) WRITE (Nounit,FMT=99007)  &
     &                    n , iwk , ioldsd , jtype , j , Result(j)
                  ENDDO
!
                  nerrs = nerrs + nfail
                  ntestt = ntestt + ntest
!
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!     Summary
!
      CALL SLASUM(path,Nounit,nerrs,ntestt)
!
99001 FORMAT (/1X,A3,' -- Real Schur Form Decomposition Driver',        &
     &        /' Matrix types (see SDRVES for details): ')
!
99002 FORMAT (/' Special Matrices:',/'  1=Zero matrix.             ',   &
     &        '           ','  5=Diagonal: geometr. spaced entries.',   &
     &        /'  2=Identity matrix.                    ','  6=Diagona',&
     &        'l: clustered entries.',/'  3=Transposed Jordan block.  ',&
     &        '          ','  7=Diagonal: large, evenly spaced.',/'  ', &
     &        '4=Diagonal: evenly spaced entries.    ',                 &
     &        '  8=Diagonal: s','mall, evenly spaced.')
99003 FORMAT (' Dense, Non-Symmetric Matrices:',/'  9=Well-cond., ev',  &
     &        'enly spaced eigenvals.',                                 &
     &        ' 14=Ill-cond., geomet. spaced e','igenals.',             &
     &        /' 10=Well-cond., geom. spaced eigenvals. ',              &
     &        ' 15=Ill-conditioned, clustered e.vals.',/' 11=Well-cond',&
     &        'itioned, clustered e.vals. ',                            &
     &        ' 16=Ill-cond., random comp','lex ',                      &
     &        /' 12=Well-cond., random complex ',6X,'   ',              &
     &        ' 17=Ill-cond., large rand. complx ',/' 13=Ill-condi',    &
     &        'tioned, evenly spaced.     ',                            &
     &        ' 18=Ill-cond., small rand.',' complx ')
99004 FORMAT (' 19=Matrix with random O(1) entries.    ',' 21=Matrix ', &
     &        'with small random entries.',/' 20=Matrix with large ran',&
     &        'dom entries.   ',/)
99005 FORMAT (' Tests performed with test threshold =',F8.2,            &
     &        /' ( A denotes A on input and T denotes A on output)',    &
     &        //' 1 = 0 if T in Schur form (no sort), ',                &
     &        '  1/ulp otherwise',/                                     &
     &       ' 2 = | A - VS T transpose(VS) | / ( n |A| ulp ) (no sort)'&
     &       ,/' 3 = | I - VS transpose(VS) | / ( n ulp ) (no sort) ',  &
     &       /' 4 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (no sort),'&
     &       ,'  1/ulp otherwise',                                      &
     &       /' 5 = 0 if T same no matter if VS computed (no sort),',   &
     &       '  1/ulp otherwise',                                       &
     &       /' 6 = 0 if WR, WI same no matter if VS computed (no sort)'&
     &       ,',  1/ulp otherwise')
99006 FORMAT (' 7 = 0 if T in Schur form (sort), ','  1/ulp otherwise', &
     &        /' 8 = | A - VS T transpose(VS) | / ( n |A| ulp ) (sort)',&
     &        /' 9 = | I - VS transpose(VS) | / ( n ulp ) (sort) ',     &
     &        /' 10 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (sort),',&
     &        '  1/ulp otherwise',                                      &
     &        /' 11 = 0 if T same no matter if VS computed (sort),',    &
     &        '  1/ulp otherwise',                                      &
     &        /' 12 = 0 if WR, WI same no matter if VS computed (sort),'&
     &        ,'  1/ulp otherwise',                                     &
     &        /' 13 = 0 if sorting successful, 1/ulp otherwise',/)
99007 FORMAT (' N=',I5,', IWK=',I2,', seed=',4(I4,','),' type ',I2,     &
     &        ', test(',I2,')=',G10.3)
99008 FORMAT (' SDRVES: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
!
!     End of SDRVES
!
      END SUBROUTINE SDRVES
