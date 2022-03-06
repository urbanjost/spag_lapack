!*==cdrvst.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b CDRVST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CDRVST( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          NOUNIT, A, LDA, D1, D2, D3, WA1, WA2, WA3, U,
!                          LDU, V, TAU, Z, WORK, LWORK, RWORK, LRWORK,
!                          IWORK, LIWORK, RESULT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDU, LIWORK, LRWORK, LWORK, NOUNIT,
!      $                   NSIZES, NTYPES
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
!       REAL               D1( * ), D2( * ), D3( * ), RESULT( * ),
!      $                   RWORK( * ), WA1( * ), WA2( * ), WA3( * )
!       COMPLEX            A( LDA, * ), TAU( * ), U( LDU, * ),
!      $                   V( LDU, * ), WORK( * ), Z( LDU, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      CDRVST  checks the Hermitian eigenvalue problem drivers.
!>
!>              CHEEVD computes all eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian matrix,
!>              using a divide-and-conquer algorithm.
!>
!>              CHEEVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian matrix.
!>
!>              CHEEVR computes selected eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian matrix
!>              using the Relatively Robust Representation where it can.
!>
!>              CHPEVD computes all eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian matrix in packed
!>              storage, using a divide-and-conquer algorithm.
!>
!>              CHPEVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian matrix in packed
!>              storage.
!>
!>              CHBEVD computes all eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian band matrix,
!>              using a divide-and-conquer algorithm.
!>
!>              CHBEVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian band matrix.
!>
!>              CHEEV computes all eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian matrix.
!>
!>              CHPEV computes all eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian matrix in packed
!>              storage.
!>
!>              CHBEV computes all eigenvalues and, optionally,
!>              eigenvectors of a complex Hermitian band matrix.
!>
!>      When CDRVST is called, a number of matrix "sizes" ("n's") and a
!>      number of matrix "types" are specified.  For each size ("n")
!>      and each type of matrix, one matrix will be generated and used
!>      to test the appropriate drivers.  For each matrix and each
!>      driver routine called, the following tests will be performed:
!>
!>      (1)     | A - Z D Z' | / ( |A| n ulp )
!>
!>      (2)     | I - Z Z' | / ( n ulp )
!>
!>      (3)     | D1 - D2 | / ( |D1| ulp )
!>
!>      where Z is the matrix of eigenvectors returned when the
!>      eigenvector option is given and D1 and D2 are the eigenvalues
!>      returned with and without the eigenvector option.
!>
!>      The "sizes" are specified by an array NN(1:NSIZES); the value of
!>      each element NN(j) specifies one size.
!>      The "types" are specified by a logical array DOTYPE( 1:NTYPES );
!>      if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!>      Currently, the list of possible types is:
!>
!>      (1)  The zero matrix.
!>      (2)  The identity matrix.
!>
!>      (3)  A diagonal matrix with evenly spaced entries
!>           1, ..., ULP  and random signs.
!>           (ULP = (first number larger than 1) - 1 )
!>      (4)  A diagonal matrix with geometrically spaced entries
!>           1, ..., ULP  and random signs.
!>      (5)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
!>           and random signs.
!>
!>      (6)  Same as (4), but multiplied by SQRT( overflow threshold )
!>      (7)  Same as (4), but multiplied by SQRT( underflow threshold )
!>
!>      (8)  A matrix of the form  U* D U, where U is unitary and
!>           D has evenly spaced entries 1, ..., ULP with random signs
!>           on the diagonal.
!>
!>      (9)  A matrix of the form  U* D U, where U is unitary and
!>           D has geometrically spaced entries 1, ..., ULP with random
!>           signs on the diagonal.
!>
!>      (10) A matrix of the form  U* D U, where U is unitary and
!>           D has "clustered" entries 1, ULP,..., ULP with random
!>           signs on the diagonal.
!>
!>      (11) Same as (8), but multiplied by SQRT( overflow threshold )
!>      (12) Same as (8), but multiplied by SQRT( underflow threshold )
!>
!>      (13) Symmetric matrix with random entries chosen from (-1,1).
!>      (14) Same as (13), but multiplied by SQRT( overflow threshold )
!>      (15) Same as (13), but multiplied by SQRT( underflow threshold )
!>      (16) A band matrix with half bandwidth randomly chosen between
!>           0 and N-1, with evenly spaced eigenvalues 1, ..., ULP
!>           with random signs.
!>      (17) Same as (16), but multiplied by SQRT( overflow threshold )
!>      (18) Same as (16), but multiplied by SQRT( underflow threshold )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  NSIZES  INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          CDRVST does nothing.  It must be at least zero.
!>          Not modified.
!>
!>  NN      INTEGER array, dimension (NSIZES)
!>          An array containing the sizes to be used for the matrices.
!>          Zero values will be skipped.  The values must be at least
!>          zero.
!>          Not modified.
!>
!>  NTYPES  INTEGER
!>          The number of elements in DOTYPE.   If it is zero, CDRVST
!>          does nothing.  It must be at least zero.  If it is MAXTYP+1
!>          and NSIZES is 1, then an additional type, MAXTYP+1 is
!>          defined, which is to use whatever matrix is in A.  This
!>          is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>          DOTYPE(MAXTYP+1) is .TRUE. .
!>          Not modified.
!>
!>  DOTYPE  LOGICAL array, dimension (NTYPES)
!>          If DOTYPE(j) is .TRUE., then for each size in NN a
!>          matrix of that size and of type j will be generated.
!>          If NTYPES is smaller than the maximum number of types
!>          defined (PARAMETER MAXTYP), then types NTYPES+1 through
!>          MAXTYP will not be generated.  If NTYPES is larger
!>          than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
!>          will be ignored.
!>          Not modified.
!>
!>  ISEED   INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator. The array elements should be between 0 and 4095;
!>          if not they will be reduced mod 4096.  Also, ISEED(4) must
!>          be odd.  The random number generator uses a linear
!>          congruential sequence limited to small integers, and so
!>          should produce machine independent random numbers. The
!>          values of ISEED are changed on exit, and can be used in the
!>          next call to CDRVST to continue the same random number
!>          sequence.
!>          Modified.
!>
!>  THRESH  REAL
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error
!>          is scaled to be O(1), so THRESH should be a reasonably
!>          small multiple of 1, e.g., 10 or 100.  In particular,
!>          it should not depend on the precision (single vs. double)
!>          or the size of the matrix.  It must be at least zero.
!>          Not modified.
!>
!>  NOUNIT  INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns IINFO not equal to 0.)
!>          Not modified.
!>
!>  A       COMPLEX array, dimension (LDA , max(NN))
!>          Used to hold the matrix whose eigenvalues are to be
!>          computed.  On exit, A contains the last matrix actually
!>          used.
!>          Modified.
!>
!>  LDA     INTEGER
!>          The leading dimension of A.  It must be at
!>          least 1 and at least max( NN ).
!>          Not modified.
!>
!>  D1      REAL array, dimension (max(NN))
!>          The eigenvalues of A, as computed by CSTEQR simlutaneously
!>          with Z.  On exit, the eigenvalues in D1 correspond with the
!>          matrix in A.
!>          Modified.
!>
!>  D2      REAL array, dimension (max(NN))
!>          The eigenvalues of A, as computed by CSTEQR if Z is not
!>          computed.  On exit, the eigenvalues in D2 correspond with
!>          the matrix in A.
!>          Modified.
!>
!>  D3      REAL array, dimension (max(NN))
!>          The eigenvalues of A, as computed by SSTERF.  On exit, the
!>          eigenvalues in D3 correspond with the matrix in A.
!>          Modified.
!>
!>  WA1     REAL array, dimension
!>
!>  WA2     REAL array, dimension
!>
!>  WA3     REAL array, dimension
!>
!>  U       COMPLEX array, dimension (LDU, max(NN))
!>          The unitary matrix computed by CHETRD + CUNGC3.
!>          Modified.
!>
!>  LDU     INTEGER
!>          The leading dimension of U, Z, and V.  It must be at
!>          least 1 and at least max( NN ).
!>          Not modified.
!>
!>  V       COMPLEX array, dimension (LDU, max(NN))
!>          The Housholder vectors computed by CHETRD in reducing A to
!>          tridiagonal form.
!>          Modified.
!>
!>  TAU     COMPLEX array, dimension (max(NN))
!>          The Householder factors computed by CHETRD in reducing A
!>          to tridiagonal form.
!>          Modified.
!>
!>  Z       COMPLEX array, dimension (LDU, max(NN))
!>          The unitary matrix of eigenvectors computed by CHEEVD,
!>          CHEEVX, CHPEVD, CHPEVX, CHBEVD, and CHBEVX.
!>          Modified.
!>
!>  WORK  - COMPLEX array of dimension ( LWORK )
!>           Workspace.
!>           Modified.
!>
!>  LWORK - INTEGER
!>           The number of entries in WORK.  This must be at least
!>           2*max( NN(j), 2 )**2.
!>           Not modified.
!>
!>  RWORK   REAL array, dimension (3*max(NN))
!>           Workspace.
!>           Modified.
!>
!>  LRWORK - INTEGER
!>           The number of entries in RWORK.
!>
!>  IWORK   INTEGER array, dimension (6*max(NN))
!>          Workspace.
!>          Modified.
!>
!>  LIWORK - INTEGER
!>           The number of entries in IWORK.
!>
!>  RESULT  REAL array, dimension (??)
!>          The values computed by the tests described above.
!>          The values are currently limited to 1/ulp, to avoid
!>          overflow.
!>          Modified.
!>
!>  INFO    INTEGER
!>          If 0, then everything ran OK.
!>           -1: NSIZES < 0
!>           -2: Some NN(j) < 0
!>           -3: NTYPES < 0
!>           -5: THRESH < 0
!>           -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
!>          -16: LDU < 1 or LDU < NMAX.
!>          -21: LWORK too small.
!>          If  SLATMR, SLATMS, CHETRD, SORGC3, CSTEQR, SSTERF,
!>              or SORMC2 returns an error code, the
!>              absolute value of it is returned.
!>          Modified.
!>
!>-----------------------------------------------------------------------
!>
!>       Some Local Variables and Parameters:
!>       ---- ----- --------- --- ----------
!>       ZERO, ONE       Real 0 and 1.
!>       MAXTYP          The number of types defined.
!>       NTEST           The number of tests performed, or which can
!>                       be performed so far, for the current matrix.
!>       NTESTT          The total number of tests performed so far.
!>       NMAX            Largest value in NN.
!>       NMATS           The number of matrices generated so far.
!>       NERRS           The number of tests which have exceeded THRESH
!>                       so far (computed by SLAFTS).
!>       COND, IMODE     Values to be passed to the matrix generators.
!>       ANORM           Norm of A; passed to matrix generators.
!>
!>       OVFL, UNFL      Overflow and underflow thresholds.
!>       ULP, ULPINV     Finest relative precision and its inverse.
!>       RTOVFL, RTUNFL  Square roots of the previous 2 values.
!>               The following four arrays decode JTYPE:
!>       KTYPE(j)        The general type (1-10) for type "j".
!>       KMODE(j)        The MODE value to be passed to the matrix
!>                       generator for type "j".
!>       KMAGN(j)        The order of magnitude ( O(1),
!>                       O(overflow^(1/2) ), O(underflow^(1/2) )
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE CDRVST(Nsizes,Nn,Ntypes,Dotype,Iseed,Thresh,Nounit,A,  &
     &                  Lda,D1,D2,D3,Wa1,Wa2,Wa3,U,Ldu,V,Tau,Z,Work,    &
     &                  Lwork,Rwork,Lrwork,Iwork,Liwork,Result,Info)
      IMPLICIT NONE
!*--CDRVST341
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldu , Liwork , Lrwork , Lwork , Nounit ,     &
     &        Nsizes , Ntypes
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iseed(4) , Iwork(*) , Nn(*)
      REAL D1(*) , D2(*) , D3(*) , Result(*) , Rwork(*) , Wa1(*) ,      &
     &     Wa2(*) , Wa3(*)
      COMPLEX A(Lda,*) , Tau(*) , U(Ldu,*) , V(Ldu,*) , Work(*) ,       &
     &        Z(Ldu,*)
!     ..
!
!  =====================================================================
!
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO , TEN
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0,TWO=2.0E+0,TEN=10.0E+0)
      REAL HALF
      PARAMETER (HALF=ONE/TWO)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
      INTEGER MAXTYP
      PARAMETER (MAXTYP=18)
!     ..
!     .. Local Scalars ..
      LOGICAL badnn
      CHARACTER uplo
      INTEGER i , idiag , ihbw , iinfo , il , imode , indwrk , indx ,   &
     &        irow , itemp , itype , iu , iuplo , j , j1 , j2 , jcol ,  &
     &        jsize , jtype , kd , lgn , liwedc , lrwedc , lwedc , m ,  &
     &        m2 , m3 , mtypes , n , nerrs , nmats , nmax , ntest ,     &
     &        ntestt
      REAL abstol , aninv , anorm , cond , ovfl , rtovfl , rtunfl ,     &
     &     temp1 , temp2 , temp3 , ulp , ulpinv , unfl , vl , vu
!     ..
!     .. Local Arrays ..
      INTEGER idumma(1) , ioldsd(4) , iseed2(4) , iseed3(4) ,           &
     &        kmagn(MAXTYP) , kmode(MAXTYP) , ktype(MAXTYP)
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLARND , SSXT1
      EXTERNAL SLAMCH , SLARND , SSXT1
!     ..
!     .. External Subroutines ..
      EXTERNAL ALASVM , CHBEV , CHBEVD , CHBEVX , CHEEV , CHEEVD ,      &
     &         CHEEVR , CHEEVX , CHET21 , CHET22 , CHPEV , CHPEVD ,     &
     &         CHPEVX , CLACPY , CLASET , CLATMR , CLATMS , SLABAD ,    &
     &         SLAFTS , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , INT , LOG , MAX , MIN , REAL , SQRT
!     ..
!     .. Data statements ..
      DATA ktype/1 , 2 , 5*4 , 5*5 , 3*8 , 3*9/
      DATA kmagn/2*1 , 1 , 1 , 1 , 2 , 3 , 1 , 1 , 1 , 2 , 3 , 1 , 2 ,  &
     &     3 , 1 , 2 , 3/
      DATA kmode/2*0 , 4 , 3 , 1 , 4 , 4 , 4 , 3 , 1 , 4 , 4 , 0 , 0 ,  &
     &     0 , 4 , 4 , 4/
!     ..
!     .. Executable Statements ..
!
!     1)      Check for errors
!
      ntestt = 0
      Info = 0
!
      badnn = .FALSE.
      nmax = 1
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
      ELSEIF ( Lda<nmax ) THEN
         Info = -9
      ELSEIF ( Ldu<nmax ) THEN
         Info = -16
      ELSEIF ( 2*MAX(2,nmax)**2>Lwork ) THEN
         Info = -22
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CDRVST',-Info)
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
      ovfl = SLAMCH('Overflow')
      CALL SLABAD(unfl,ovfl)
      ulp = SLAMCH('Epsilon')*SLAMCH('Base')
      ulpinv = ONE/ulp
      rtunfl = SQRT(unfl)
      rtovfl = SQRT(ovfl)
!
!     Loop over sizes, types
!
      DO i = 1 , 4
         iseed2(i) = Iseed(i)
         iseed3(i) = Iseed(i)
      ENDDO
!
      nerrs = 0
      nmats = 0
!
      DO jsize = 1 , Nsizes
         n = Nn(jsize)
         IF ( n>0 ) THEN
            lgn = INT(LOG(REAL(n))/LOG(TWO))
            IF ( 2**lgn<n ) lgn = lgn + 1
            IF ( 2**lgn<n ) lgn = lgn + 1
            lwedc = MAX(2*n+n*n,2*n*n)
            lrwedc = 1 + 4*n + 2*n*lgn + 3*n**2
            liwedc = 3 + 5*n
         ELSE
            lwedc = 2
            lrwedc = 8
            liwedc = 8
         ENDIF
         aninv = ONE/REAL(MAX(1,n))
!
         IF ( Nsizes/=1 ) THEN
            mtypes = MIN(MAXTYP,Ntypes)
         ELSE
            mtypes = MIN(MAXTYP+1,Ntypes)
         ENDIF
!
         DO jtype = 1 , mtypes
            IF ( Dotype(jtype) ) THEN
               nmats = nmats + 1
               ntest = 0
!
               DO j = 1 , 4
                  ioldsd(j) = Iseed(j)
               ENDDO
!
!           2)      Compute "A"
!
!                   Control parameters:
!
!               KMAGN  KMODE        KTYPE
!           =1  O(1)   clustered 1  zero
!           =2  large  clustered 2  identity
!           =3  small  exponential  (none)
!           =4         arithmetic   diagonal, (w/ eigenvalues)
!           =5         random log   Hermitian, w/ eigenvalues
!           =6         random       (none)
!           =7                      random diagonal
!           =8                      random Hermitian
!           =9                      band Hermitian, w/ eigenvalues
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
                     anorm = (rtovfl*ulp)*aninv
                  ELSEIF ( kmagn(jtype)==3 ) THEN
!
                     anorm = rtunfl*n*ulpinv
                  ELSE
!
                     anorm = ONE
                  ENDIF
!
!
                  CALL CLASET('Full',Lda,n,CZERO,CZERO,A,Lda)
                  iinfo = 0
                  cond = ulpinv
!
!           Special Matrices -- Identity & Jordan block
!
!                   Zero
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
                  ELSEIF ( itype==4 ) THEN
!
!              Diagonal Matrix, [Eigen]values Specified
!
                     CALL CLATMS(n,n,'S',Iseed,'H',Rwork,imode,cond,    &
     &                           anorm,0,0,'N',A,Lda,Work,iinfo)
!
                  ELSEIF ( itype==5 ) THEN
!
!              Hermitian, eigenvalues specified
!
                     CALL CLATMS(n,n,'S',Iseed,'H',Rwork,imode,cond,    &
     &                           anorm,n,n,'N',A,Lda,Work,iinfo)
!
                  ELSEIF ( itype==7 ) THEN
!
!              Diagonal, random eigenvalues
!
                     CALL CLATMR(n,n,'S',Iseed,'H',Work,6,ONE,CONE,'T', &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,0,0,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==8 ) THEN
!
!              Hermitian, random eigenvalues
!
                     CALL CLATMR(n,n,'S',Iseed,'H',Work,6,ONE,CONE,'T', &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,n,n,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==9 ) THEN
!
!              Hermitian banded, eigenvalues specified
!
                     ihbw = INT((n-1)*SLARND(1,iseed3))
                     CALL CLATMS(n,n,'S',Iseed,'H',Rwork,imode,cond,    &
     &                           anorm,ihbw,ihbw,'Z',U,Ldu,Work,iinfo)
!
!              Store as dense matrix for most routines.
!
                     CALL CLASET('Full',Lda,n,CZERO,CZERO,A,Lda)
                     DO idiag = -ihbw , ihbw
                        irow = ihbw - idiag + 1
                        j1 = MAX(1,idiag+1)
                        j2 = MIN(n,n+idiag)
                        DO j = j1 , j2
                           i = j - idiag
                           A(i,j) = U(irow,j)
                        ENDDO
                     ENDDO
                  ELSE
                     iinfo = 1
                  ENDIF
!
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'Generator' , iinfo , n , &
     &                      jtype , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
               ENDIF
!
!
               abstol = unfl + unfl
               IF ( n<=1 ) THEN
                  il = 1
                  iu = n
               ELSE
                  il = 1 + INT((n-1)*SLARND(1,iseed2))
                  iu = 1 + INT((n-1)*SLARND(1,iseed2))
                  IF ( il>iu ) THEN
                     itemp = il
                     il = iu
                     iu = itemp
                  ENDIF
               ENDIF
!
!           Perform tests storing upper or lower triangular
!           part of matrix.
!
               DO iuplo = 0 , 1
                  IF ( iuplo==0 ) THEN
                     uplo = 'L'
                  ELSE
                     uplo = 'U'
                  ENDIF
!
!              Call CHEEVD and CHEEVX.
!
                  CALL CLACPY(' ',n,n,A,Lda,V,Ldu)
!
                  ntest = ntest + 1
                  CALL CHEEVD('V',uplo,n,A,Ldu,D1,Work,lwedc,Rwork,     &
     &                        lrwedc,Iwork,liwedc,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVD(V,'//uplo//')' ,  &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 5
                     ENDIF
                  ENDIF
!
!              Do tests 1 and 2.
!
                  CALL CHET21(1,uplo,n,0,V,Ldu,D1,D2,A,Ldu,Z,Ldu,Tau,   &
     &                        Work,Rwork,Result(ntest))
!
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  ntest = ntest + 2
                  CALL CHEEVD('N',uplo,n,A,Ldu,D3,Work,lwedc,Rwork,     &
     &                        lrwedc,Iwork,liwedc,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVD(N,'//uplo//')' ,  &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 5
                     ENDIF
                  ENDIF
!
!              Do test 3.
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(D1(j)),ABS(D3(j)))
                     temp2 = MAX(temp2,ABS(D1(j)-D3(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
 5                CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  ntest = ntest + 1
!
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(D1(1)),ABS(D1(n)))
                     IF ( il/=1 ) THEN
                        vl = D1(il)                                     &
     &                       - MAX(HALF*(D1(il)-D1(il-1)),TEN*ulp*temp3,&
     &                       TEN*rtunfl)
                     ELSEIF ( n>0 ) THEN
                        vl = D1(1)                                      &
     &                       - MAX(HALF*(D1(n)-D1(1)),TEN*ulp*temp3,    &
     &                       TEN*rtunfl)
                     ENDIF
                     IF ( iu/=n ) THEN
                        vu = D1(iu)                                     &
     &                       + MAX(HALF*(D1(iu+1)-D1(iu)),TEN*ulp*temp3,&
     &                       TEN*rtunfl)
                     ELSEIF ( n>0 ) THEN
                        vu = D1(n)                                      &
     &                       + MAX(HALF*(D1(n)-D1(1)),TEN*ulp*temp3,    &
     &                       TEN*rtunfl)
                     ENDIF
                  ELSE
                     temp3 = ZERO
                     vl = ZERO
                     vu = ONE
                  ENDIF
!
                  CALL CHEEVX('V','A',uplo,n,A,Ldu,vl,vu,il,iu,abstol,m,&
     &                        Wa1,Z,Ldu,Work,Lwork,Rwork,Iwork,         &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVX(V,A,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 10
                     ENDIF
                  ENDIF
!
!              Do tests 4 and 5.
!
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  CALL CHET21(1,uplo,n,0,A,Ldu,Wa1,D2,Z,Ldu,V,Ldu,Tau,  &
     &                        Work,Rwork,Result(ntest))
!
                  ntest = ntest + 2
                  CALL CHEEVX('N','A',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m2,Wa2,Z,Ldu,Work,Lwork,Rwork,Iwork,      &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVX(N,A,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 10
                     ENDIF
                  ENDIF
!
!              Do test 6.
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(Wa1(j)),ABS(Wa2(j)))
                     temp2 = MAX(temp2,ABS(Wa1(j)-Wa2(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
 10               CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  ntest = ntest + 1
!
                  CALL CHEEVX('V','I',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m2,Wa2,Z,Ldu,Work,Lwork,Rwork,Iwork,      &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVX(V,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 15
                     ENDIF
                  ENDIF
!
!              Do tests 7 and 8.
!
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  CALL CHET22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Rwork,Result(ntest))
!
                  ntest = ntest + 2
!
                  CALL CHEEVX('N','I',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m3,Wa3,Z,Ldu,Work,Lwork,Rwork,Iwork,      &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVX(N,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 15
                     ENDIF
                  ENDIF
!
!              Do test 9.
!
                  temp1 = SSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = SSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
 15               CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  ntest = ntest + 1
!
                  CALL CHEEVX('V','V',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m2,Wa2,Z,Ldu,Work,Lwork,Rwork,Iwork,      &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVX(V,V,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 20
                     ENDIF
                  ENDIF
!
!              Do tests 10 and 11.
!
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  CALL CHET22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Rwork,Result(ntest))
!
                  ntest = ntest + 2
!
                  CALL CHEEVX('N','V',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m3,Wa3,Z,Ldu,Work,Lwork,Rwork,Iwork,      &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVX(N,V,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 20
                     ENDIF
                  ENDIF
!
                  IF ( m3==0 .AND. n>0 ) THEN
                     Result(ntest) = ulpinv
                     GOTO 20
                  ENDIF
!
!              Do test 12.
!
                  temp1 = SSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = SSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
!
!              Call CHPEVD and CHPEVX.
!
 20               CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
!              Load array WORK with the upper or lower triangular
!              part of the matrix in packed form.
!
                  IF ( iuplo==1 ) THEN
                     indx = 1
                     DO j = 1 , n
                        DO i = 1 , j
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ELSE
                     indx = 1
                     DO j = 1 , n
                        DO i = j , n
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ENDIF
!
                  ntest = ntest + 1
                  indwrk = n*(n+1)/2 + 1
                  CALL CHPEVD('V',uplo,n,Work,D1,Z,Ldu,Work(indwrk),    &
     &                        lwedc,Rwork,lrwedc,Iwork,liwedc,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHPEVD(V,'//uplo//')' ,  &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 25
                     ENDIF
                  ENDIF
!
!              Do tests 13 and 14.
!
                  CALL CHET21(1,uplo,n,0,A,Lda,D1,D2,Z,Ldu,V,Ldu,Tau,   &
     &                        Work,Rwork,Result(ntest))
!
                  IF ( iuplo==1 ) THEN
                     indx = 1
                     DO j = 1 , n
                        DO i = 1 , j
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ELSE
                     indx = 1
                     DO j = 1 , n
                        DO i = j , n
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ENDIF
!
                  ntest = ntest + 2
                  indwrk = n*(n+1)/2 + 1
                  CALL CHPEVD('N',uplo,n,Work,D3,Z,Ldu,Work(indwrk),    &
     &                        lwedc,Rwork,lrwedc,Iwork,liwedc,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHPEVD(N,'//uplo//')' ,  &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 25
                     ENDIF
                  ENDIF
!
!              Do test 15.
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(D1(j)),ABS(D3(j)))
                     temp2 = MAX(temp2,ABS(D1(j)-D3(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!              Load array WORK with the upper or lower triangular part
!              of the matrix in packed form.
!
 25               IF ( iuplo==1 ) THEN
                     indx = 1
                     DO j = 1 , n
                        DO i = 1 , j
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ELSE
                     indx = 1
                     DO j = 1 , n
                        DO i = j , n
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ENDIF
!
                  ntest = ntest + 1
!
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(D1(1)),ABS(D1(n)))
                     IF ( il/=1 ) THEN
                        vl = D1(il)                                     &
     &                       - MAX(HALF*(D1(il)-D1(il-1)),TEN*ulp*temp3,&
     &                       TEN*rtunfl)
                     ELSEIF ( n>0 ) THEN
                        vl = D1(1)                                      &
     &                       - MAX(HALF*(D1(n)-D1(1)),TEN*ulp*temp3,    &
     &                       TEN*rtunfl)
                     ENDIF
                     IF ( iu/=n ) THEN
                        vu = D1(iu)                                     &
     &                       + MAX(HALF*(D1(iu+1)-D1(iu)),TEN*ulp*temp3,&
     &                       TEN*rtunfl)
                     ELSEIF ( n>0 ) THEN
                        vu = D1(n)                                      &
     &                       + MAX(HALF*(D1(n)-D1(1)),TEN*ulp*temp3,    &
     &                       TEN*rtunfl)
                     ENDIF
                  ELSE
                     temp3 = ZERO
                     vl = ZERO
                     vu = ONE
                  ENDIF
!
                  CALL CHPEVX('V','A',uplo,n,Work,vl,vu,il,iu,abstol,m, &
     &                        Wa1,Z,Ldu,V,Rwork,Iwork,Iwork(5*n+1),     &
     &                        iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHPEVX(V,A,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 30
                     ENDIF
                  ENDIF
!
!              Do tests 16 and 17.
!
                  CALL CHET21(1,uplo,n,0,A,Ldu,Wa1,D2,Z,Ldu,V,Ldu,Tau,  &
     &                        Work,Rwork,Result(ntest))
!
                  ntest = ntest + 2
!
                  IF ( iuplo==1 ) THEN
                     indx = 1
                     DO j = 1 , n
                        DO i = 1 , j
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ELSE
                     indx = 1
                     DO j = 1 , n
                        DO i = j , n
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ENDIF
!
                  CALL CHPEVX('N','A',uplo,n,Work,vl,vu,il,iu,abstol,m2,&
     &                        Wa2,Z,Ldu,V,Rwork,Iwork,Iwork(5*n+1),     &
     &                        iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHPEVX(N,A,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 30
                     ENDIF
                  ENDIF
!
!              Do test 18.
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(Wa1(j)),ABS(Wa2(j)))
                     temp2 = MAX(temp2,ABS(Wa1(j)-Wa2(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
 30               ntest = ntest + 1
                  IF ( iuplo==1 ) THEN
                     indx = 1
                     DO j = 1 , n
                        DO i = 1 , j
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ELSE
                     indx = 1
                     DO j = 1 , n
                        DO i = j , n
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ENDIF
!
                  CALL CHPEVX('V','I',uplo,n,Work,vl,vu,il,iu,abstol,m2,&
     &                        Wa2,Z,Ldu,V,Rwork,Iwork,Iwork(5*n+1),     &
     &                        iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHPEVX(V,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 35
                     ENDIF
                  ENDIF
!
!              Do tests 19 and 20.
!
                  CALL CHET22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Rwork,Result(ntest))
!
                  ntest = ntest + 2
!
                  IF ( iuplo==1 ) THEN
                     indx = 1
                     DO j = 1 , n
                        DO i = 1 , j
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ELSE
                     indx = 1
                     DO j = 1 , n
                        DO i = j , n
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ENDIF
!
                  CALL CHPEVX('N','I',uplo,n,Work,vl,vu,il,iu,abstol,m3,&
     &                        Wa3,Z,Ldu,V,Rwork,Iwork,Iwork(5*n+1),     &
     &                        iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHPEVX(N,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 35
                     ENDIF
                  ENDIF
!
!              Do test 21.
!
                  temp1 = SSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = SSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
 35               ntest = ntest + 1
                  IF ( iuplo==1 ) THEN
                     indx = 1
                     DO j = 1 , n
                        DO i = 1 , j
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ELSE
                     indx = 1
                     DO j = 1 , n
                        DO i = j , n
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ENDIF
!
                  CALL CHPEVX('V','V',uplo,n,Work,vl,vu,il,iu,abstol,m2,&
     &                        Wa2,Z,Ldu,V,Rwork,Iwork,Iwork(5*n+1),     &
     &                        iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHPEVX(V,V,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 40
                     ENDIF
                  ENDIF
!
!              Do tests 22 and 23.
!
                  CALL CHET22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Rwork,Result(ntest))
!
                  ntest = ntest + 2
!
                  IF ( iuplo==1 ) THEN
                     indx = 1
                     DO j = 1 , n
                        DO i = 1 , j
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ELSE
                     indx = 1
                     DO j = 1 , n
                        DO i = j , n
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ENDIF
!
                  CALL CHPEVX('N','V',uplo,n,Work,vl,vu,il,iu,abstol,m3,&
     &                        Wa3,Z,Ldu,V,Rwork,Iwork,Iwork(5*n+1),     &
     &                        iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHPEVX(N,V,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 40
                     ENDIF
                  ENDIF
!
                  IF ( m3==0 .AND. n>0 ) THEN
                     Result(ntest) = ulpinv
                     GOTO 40
                  ENDIF
!
!              Do test 24.
!
                  temp1 = SSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = SSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
!
!              Call CHBEVD and CHBEVX.
!
 40               IF ( jtype<=7 ) THEN
                     kd = 0
                  ELSEIF ( jtype>=8 .AND. jtype<=15 ) THEN
                     kd = MAX(n-1,0)
                  ELSE
                     kd = ihbw
                  ENDIF
!
!              Load array V with the upper or lower triangular part
!              of the matrix in band form.
!
                  IF ( iuplo==1 ) THEN
                     DO j = 1 , n
                        DO i = MAX(1,j-kd) , j
                           V(kd+1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ELSE
                     DO j = 1 , n
                        DO i = j , MIN(n,j+kd)
                           V(1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ENDIF
!
                  ntest = ntest + 1
                  CALL CHBEVD('V',uplo,n,kd,V,Ldu,D1,Z,Ldu,Work,lwedc,  &
     &                        Rwork,lrwedc,Iwork,liwedc,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99002) 'CHBEVD(V,'//uplo//')' ,  &
     &                      iinfo , n , kd , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 45
                     ENDIF
                  ENDIF
!
!              Do tests 25 and 26.
!
                  CALL CHET21(1,uplo,n,0,A,Lda,D1,D2,Z,Ldu,V,Ldu,Tau,   &
     &                        Work,Rwork,Result(ntest))
!
                  IF ( iuplo==1 ) THEN
                     DO j = 1 , n
                        DO i = MAX(1,j-kd) , j
                           V(kd+1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ELSE
                     DO j = 1 , n
                        DO i = j , MIN(n,j+kd)
                           V(1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ENDIF
!
                  ntest = ntest + 2
                  CALL CHBEVD('N',uplo,n,kd,V,Ldu,D3,Z,Ldu,Work,lwedc,  &
     &                        Rwork,lrwedc,Iwork,liwedc,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99002) 'CHBEVD(N,'//uplo//')' ,  &
     &                      iinfo , n , kd , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 45
                     ENDIF
                  ENDIF
!
!              Do test 27.
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(D1(j)),ABS(D3(j)))
                     temp2 = MAX(temp2,ABS(D1(j)-D3(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!              Load array V with the upper or lower triangular part
!              of the matrix in band form.
!
 45               IF ( iuplo==1 ) THEN
                     DO j = 1 , n
                        DO i = MAX(1,j-kd) , j
                           V(kd+1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ELSE
                     DO j = 1 , n
                        DO i = j , MIN(n,j+kd)
                           V(1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ENDIF
!
                  ntest = ntest + 1
                  CALL CHBEVX('V','A',uplo,n,kd,V,Ldu,U,Ldu,vl,vu,il,iu,&
     &                        abstol,m,Wa1,Z,Ldu,Work,Rwork,Iwork,      &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHBEVX(V,A,'//uplo//')' ,&
     &                      iinfo , n , kd , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 50
                     ENDIF
                  ENDIF
!
!              Do tests 28 and 29.
!
                  CALL CHET21(1,uplo,n,0,A,Ldu,Wa1,D2,Z,Ldu,V,Ldu,Tau,  &
     &                        Work,Rwork,Result(ntest))
!
                  ntest = ntest + 2
!
                  IF ( iuplo==1 ) THEN
                     DO j = 1 , n
                        DO i = MAX(1,j-kd) , j
                           V(kd+1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ELSE
                     DO j = 1 , n
                        DO i = j , MIN(n,j+kd)
                           V(1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ENDIF
!
                  CALL CHBEVX('N','A',uplo,n,kd,V,Ldu,U,Ldu,vl,vu,il,iu,&
     &                        abstol,m2,Wa2,Z,Ldu,Work,Rwork,Iwork,     &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99002) 'CHBEVX(N,A,'//uplo//')' ,&
     &                      iinfo , n , kd , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 50
                     ENDIF
                  ENDIF
!
!              Do test 30.
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(Wa1(j)),ABS(Wa2(j)))
                     temp2 = MAX(temp2,ABS(Wa1(j)-Wa2(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!              Load array V with the upper or lower triangular part
!              of the matrix in band form.
!
 50               ntest = ntest + 1
                  IF ( iuplo==1 ) THEN
                     DO j = 1 , n
                        DO i = MAX(1,j-kd) , j
                           V(kd+1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ELSE
                     DO j = 1 , n
                        DO i = j , MIN(n,j+kd)
                           V(1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ENDIF
!
                  CALL CHBEVX('V','I',uplo,n,kd,V,Ldu,U,Ldu,vl,vu,il,iu,&
     &                        abstol,m2,Wa2,Z,Ldu,Work,Rwork,Iwork,     &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99002) 'CHBEVX(V,I,'//uplo//')' ,&
     &                      iinfo , n , kd , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 55
                     ENDIF
                  ENDIF
!
!              Do tests 31 and 32.
!
                  CALL CHET22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Rwork,Result(ntest))
!
                  ntest = ntest + 2
!
                  IF ( iuplo==1 ) THEN
                     DO j = 1 , n
                        DO i = MAX(1,j-kd) , j
                           V(kd+1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ELSE
                     DO j = 1 , n
                        DO i = j , MIN(n,j+kd)
                           V(1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ENDIF
                  CALL CHBEVX('N','I',uplo,n,kd,V,Ldu,U,Ldu,vl,vu,il,iu,&
     &                        abstol,m3,Wa3,Z,Ldu,Work,Rwork,Iwork,     &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99002) 'CHBEVX(N,I,'//uplo//')' ,&
     &                      iinfo , n , kd , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 55
                     ENDIF
                  ENDIF
!
!              Do test 33.
!
                  temp1 = SSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = SSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
!              Load array V with the upper or lower triangular part
!              of the matrix in band form.
!
 55               ntest = ntest + 1
                  IF ( iuplo==1 ) THEN
                     DO j = 1 , n
                        DO i = MAX(1,j-kd) , j
                           V(kd+1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ELSE
                     DO j = 1 , n
                        DO i = j , MIN(n,j+kd)
                           V(1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ENDIF
                  CALL CHBEVX('V','V',uplo,n,kd,V,Ldu,U,Ldu,vl,vu,il,iu,&
     &                        abstol,m2,Wa2,Z,Ldu,Work,Rwork,Iwork,     &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99002) 'CHBEVX(V,V,'//uplo//')' ,&
     &                      iinfo , n , kd , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 60
                     ENDIF
                  ENDIF
!
!              Do tests 34 and 35.
!
                  CALL CHET22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Rwork,Result(ntest))
!
                  ntest = ntest + 2
!
                  IF ( iuplo==1 ) THEN
                     DO j = 1 , n
                        DO i = MAX(1,j-kd) , j
                           V(kd+1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ELSE
                     DO j = 1 , n
                        DO i = j , MIN(n,j+kd)
                           V(1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ENDIF
                  CALL CHBEVX('N','V',uplo,n,kd,V,Ldu,U,Ldu,vl,vu,il,iu,&
     &                        abstol,m3,Wa3,Z,Ldu,Work,Rwork,Iwork,     &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99002) 'CHBEVX(N,V,'//uplo//')' ,&
     &                      iinfo , n , kd , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 60
                     ENDIF
                  ENDIF
!
                  IF ( m3==0 .AND. n>0 ) THEN
                     Result(ntest) = ulpinv
                     GOTO 60
                  ENDIF
!
!              Do test 36.
!
                  temp1 = SSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = SSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
!
!              Call CHEEV
!
 60               CALL CLACPY(' ',n,n,A,Lda,V,Ldu)
!
                  ntest = ntest + 1
                  CALL CHEEV('V',uplo,n,A,Ldu,D1,Work,Lwork,Rwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEV(V,'//uplo//')' ,   &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 65
                     ENDIF
                  ENDIF
!
!              Do tests 37 and 38
!
                  CALL CHET21(1,uplo,n,0,V,Ldu,D1,D2,A,Ldu,Z,Ldu,Tau,   &
     &                        Work,Rwork,Result(ntest))
!
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  ntest = ntest + 2
                  CALL CHEEV('N',uplo,n,A,Ldu,D3,Work,Lwork,Rwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEV(N,'//uplo//')' ,   &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 65
                     ENDIF
                  ENDIF
!
!              Do test 39
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(D1(j)),ABS(D3(j)))
                     temp2 = MAX(temp2,ABS(D1(j)-D3(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!
 65               CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
!              Call CHPEV
!
!              Load array WORK with the upper or lower triangular
!              part of the matrix in packed form.
!
                  IF ( iuplo==1 ) THEN
                     indx = 1
                     DO j = 1 , n
                        DO i = 1 , j
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ELSE
                     indx = 1
                     DO j = 1 , n
                        DO i = j , n
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ENDIF
!
                  ntest = ntest + 1
                  indwrk = n*(n+1)/2 + 1
                  CALL CHPEV('V',uplo,n,Work,D1,Z,Ldu,Work(indwrk),     &
     &                       Rwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHPEV(V,'//uplo//')' ,   &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 70
                     ENDIF
                  ENDIF
!
!              Do tests 40 and 41.
!
                  CALL CHET21(1,uplo,n,0,A,Lda,D1,D2,Z,Ldu,V,Ldu,Tau,   &
     &                        Work,Rwork,Result(ntest))
!
                  IF ( iuplo==1 ) THEN
                     indx = 1
                     DO j = 1 , n
                        DO i = 1 , j
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ELSE
                     indx = 1
                     DO j = 1 , n
                        DO i = j , n
                           Work(indx) = A(i,j)
                           indx = indx + 1
                        ENDDO
                     ENDDO
                  ENDIF
!
                  ntest = ntest + 2
                  indwrk = n*(n+1)/2 + 1
                  CALL CHPEV('N',uplo,n,Work,D3,Z,Ldu,Work(indwrk),     &
     &                       Rwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHPEV(N,'//uplo//')' ,   &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 70
                     ENDIF
                  ENDIF
!
!              Do test 42
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(D1(j)),ABS(D3(j)))
                     temp2 = MAX(temp2,ABS(D1(j)-D3(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!
!              Call CHBEV
!
 70               IF ( jtype<=7 ) THEN
                     kd = 0
                  ELSEIF ( jtype>=8 .AND. jtype<=15 ) THEN
                     kd = MAX(n-1,0)
                  ELSE
                     kd = ihbw
                  ENDIF
!
!              Load array V with the upper or lower triangular part
!              of the matrix in band form.
!
                  IF ( iuplo==1 ) THEN
                     DO j = 1 , n
                        DO i = MAX(1,j-kd) , j
                           V(kd+1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ELSE
                     DO j = 1 , n
                        DO i = j , MIN(n,j+kd)
                           V(1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ENDIF
!
                  ntest = ntest + 1
                  CALL CHBEV('V',uplo,n,kd,V,Ldu,D1,Z,Ldu,Work,Rwork,   &
     &                       iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99002) 'CHBEV(V,'//uplo//')' ,   &
     &                      iinfo , n , kd , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 75
                     ENDIF
                  ENDIF
!
!              Do tests 43 and 44.
!
                  CALL CHET21(1,uplo,n,0,A,Lda,D1,D2,Z,Ldu,V,Ldu,Tau,   &
     &                        Work,Rwork,Result(ntest))
!
                  IF ( iuplo==1 ) THEN
                     DO j = 1 , n
                        DO i = MAX(1,j-kd) , j
                           V(kd+1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ELSE
                     DO j = 1 , n
                        DO i = j , MIN(n,j+kd)
                           V(1+i-j,j) = A(i,j)
                        ENDDO
                     ENDDO
                  ENDIF
!
                  ntest = ntest + 2
                  CALL CHBEV('N',uplo,n,kd,V,Ldu,D3,Z,Ldu,Work,Rwork,   &
     &                       iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99002) 'CHBEV(N,'//uplo//')' ,   &
     &                      iinfo , n , kd , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                     ENDIF
                  ENDIF
!
!
!              Do test 45.
!
 75               temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(D1(j)),ABS(D3(j)))
                     temp2 = MAX(temp2,ABS(D1(j)-D3(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
                  CALL CLACPY(' ',n,n,A,Lda,V,Ldu)
                  ntest = ntest + 1
                  CALL CHEEVR('V','A',uplo,n,A,Ldu,vl,vu,il,iu,abstol,m,&
     &                        Wa1,Z,Ldu,Iwork,Work,Lwork,Rwork,Lrwork,  &
     &                        Iwork(2*n+1),Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVR(V,A,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 80
                     ENDIF
                  ENDIF
!
!              Do tests 45 and 46 (or ... )
!
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  CALL CHET21(1,uplo,n,0,A,Ldu,Wa1,D2,Z,Ldu,V,Ldu,Tau,  &
     &                        Work,Rwork,Result(ntest))
!
                  ntest = ntest + 2
                  CALL CHEEVR('N','A',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m2,Wa2,Z,Ldu,Iwork,Work,Lwork,Rwork,      &
     &                        Lrwork,Iwork(2*n+1),Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVR(N,A,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 80
                     ENDIF
                  ENDIF
!
!              Do test 47 (or ... )
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(Wa1(j)),ABS(Wa2(j)))
                     temp2 = MAX(temp2,ABS(Wa1(j)-Wa2(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!
 80               ntest = ntest + 1
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
                  CALL CHEEVR('V','I',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m2,Wa2,Z,Ldu,Iwork,Work,Lwork,Rwork,      &
     &                        Lrwork,Iwork(2*n+1),Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVR(V,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 85
                     ENDIF
                  ENDIF
!
!              Do tests 48 and 49 (or +??)
!
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  CALL CHET22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Rwork,Result(ntest))
!
                  ntest = ntest + 2
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
                  CALL CHEEVR('N','I',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m3,Wa3,Z,Ldu,Iwork,Work,Lwork,Rwork,      &
     &                        Lrwork,Iwork(2*n+1),Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVR(N,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 85
                     ENDIF
                  ENDIF
!
!              Do test 50 (or +??)
!
                  temp1 = SSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = SSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  Result(ntest) = (temp1+temp2)/MAX(unfl,ulp*temp3)
!
 85               ntest = ntest + 1
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
                  CALL CHEEVR('V','V',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m2,Wa2,Z,Ldu,Iwork,Work,Lwork,Rwork,      &
     &                        Lrwork,Iwork(2*n+1),Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVR(V,V,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        CYCLE
                     ENDIF
                  ENDIF
!
!              Do tests 51 and 52 (or +??)
!
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  CALL CHET22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Rwork,Result(ntest))
!
                  ntest = ntest + 2
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
                  CALL CHEEVR('N','V',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m3,Wa3,Z,Ldu,Iwork,Work,Lwork,Rwork,      &
     &                        Lrwork,Iwork(2*n+1),Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CHEEVR(N,V,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        CYCLE
                     ENDIF
                  ENDIF
!
                  IF ( m3==0 .AND. n>0 ) THEN
                     Result(ntest) = ulpinv
                     CYCLE
                  ENDIF
!
!              Do test 52 (or +??)
!
                  temp1 = SSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = SSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
                  CALL CLACPY(' ',n,n,V,Ldu,A,Lda)
!
!
!
!
!              Load array V with the upper or lower triangular part
!              of the matrix in band form.
!
!
               ENDDO
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
               ntestt = ntestt + ntest
               CALL SLAFTS('CST',n,n,jtype,ntest,Result,ioldsd,Thresh,  &
     &                     Nounit,nerrs)
            ENDIF
!
         ENDDO
      ENDDO
!
!     Summary
!
      CALL ALASVM('CST',Nounit,nerrs,ntestt,0)
!
99001 FORMAT (' CDRVST: ',A,' returned INFO=',I6,/9X,'N=',I6,', JTYPE=',&
     &        I6,', ISEED=(',3(I5,','),I5,')')
99002 FORMAT (' CDRVST: ',A,' returned INFO=',I6,/9X,'N=',I6,', KD=',I6,&
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
!
!     End of CDRVST
!
      END SUBROUTINE CDRVST
