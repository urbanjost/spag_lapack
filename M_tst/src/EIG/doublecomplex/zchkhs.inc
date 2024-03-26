!*==zchkhs.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b zchkhs
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zchkhs( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          NOUNIT, A, LDA, H, T1, T2, U, LDU, Z, UZ, W1,
!                          W3, EVECTL, EVECTR, EVECTY, EVECTX, UU, TAU,
!                          WORK, NWORK, RWORK, IWORK, SELECT, RESULT,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDU, NOUNIT, NSIZES, NTYPES, NWORK
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * ), SELECT( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
!       DOUBLE PRECISION   RESULT( 14 ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), EVECTL( LDU, * ),
!      $                   EVECTR( LDU, * ), EVECTX( LDU, * ),
!      $                   EVECTY( LDU, * ), H( LDA, * ), T1( LDA, * ),
!      $                   T2( LDA, * ), TAU( * ), U( LDU, * ),
!      $                   UU( LDU, * ), UZ( LDU, * ), W1( * ), W3( * ),
!      $                   WORK( * ), Z( LDU, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    zchkhs  checks the nonsymmetric eigenvalue problem routines.
!>
!>            ZGEHRD factors A as  U H U' , where ' means conjugate
!>            transpose, H is hessenberg, and U is unitary.
!>
!>            ZUNGHR generates the unitary matrix U.
!>
!>            ZUNMHR multiplies a matrix by the unitary matrix U.
!>
!>            ZHSEQR factors H as  Z T Z' , where Z is unitary and T
!>            is upper triangular.  It also computes the eigenvalues,
!>            w(1), ..., w(n); we define a diagonal matrix W whose
!>            (diagonal) entries are the eigenvalues.
!>
!>            ZTREVC computes the left eigenvector matrix L and the
!>            right eigenvector matrix R for the matrix T.  The
!>            columns of L are the complex conjugates of the left
!>            eigenvectors of T.  The columns of R are the right
!>            eigenvectors of T.  L is lower triangular, and R is
!>            upper triangular.
!>
!>            ZHSEIN computes the left eigenvector matrix Y and the
!>            right eigenvector matrix X for the matrix H.  The
!>            columns of Y are the complex conjugates of the left
!>            eigenvectors of H.  The columns of X are the right
!>            eigenvectors of H.  Y is lower triangular, and X is
!>            upper triangular.
!>
!>    When zchkhs is called, a number of matrix "sizes" ("n's") and a
!>    number of matrix "types" are specified.  For each size ("n")
!>    and each type of matrix, one matrix will be generated and used
!>    to test the nonsymmetric eigenroutines.  For each matrix, 14
!>    tests will be performed:
!>
!>    (1)     | A - U H U**H | / ( |A| n ulp )
!>
!>    (2)     | I - UU**H | / ( n ulp )
!>
!>    (3)     | H - Z T Z**H | / ( |H| n ulp )
!>
!>    (4)     | I - ZZ**H | / ( n ulp )
!>
!>    (5)     | A - UZ H (UZ)**H | / ( |A| n ulp )
!>
!>    (6)     | I - UZ (UZ)**H | / ( n ulp )
!>
!>    (7)     | T(Z computed) - T(Z not computed) | / ( |T| ulp )
!>
!>    (8)     | W(Z computed) - W(Z not computed) | / ( |W| ulp )
!>
!>    (9)     | TR - RW | / ( |T| |R| ulp )
!>
!>    (10)    | L**H T - W**H L | / ( |T| |L| ulp )
!>
!>    (11)    | HX - XW | / ( |H| |X| ulp )
!>
!>    (12)    | Y**H H - W**H Y | / ( |H| |Y| ulp )
!>
!>    (13)    | AX - XW | / ( |A| |X| ulp )
!>
!>    (14)    | Y**H A - W**H Y | / ( |A| |Y| ulp )
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
!>         1, ..., ULP  and random complex angles.
!>         (ULP = (first number larger than 1) - 1 )
!>    (5)  A diagonal matrix with geometrically spaced entries
!>         1, ..., ULP  and random complex angles.
!>    (6)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
!>         and random complex angles.
!>
!>    (7)  Same as (4), but multiplied by SQRT( overflow threshold )
!>    (8)  Same as (4), but multiplied by SQRT( underflow threshold )
!>
!>    (9)  A matrix of the form  U' T U, where U is unitary and
!>         T has evenly spaced entries 1, ..., ULP with random complex
!>         angles on the diagonal and random O(1) entries in the upper
!>         triangle.
!>
!>    (10) A matrix of the form  U' T U, where U is unitary and
!>         T has geometrically spaced entries 1, ..., ULP with random
!>         complex angles on the diagonal and random O(1) entries in
!>         the upper triangle.
!>
!>    (11) A matrix of the form  U' T U, where U is unitary and
!>         T has "clustered" entries 1, ULP,..., ULP with random
!>         complex angles on the diagonal and random O(1) entries in
!>         the upper triangle.
!>
!>    (12) A matrix of the form  U' T U, where U is unitary and
!>         T has complex eigenvalues randomly chosen from
!>         ULP < |z| < 1   and random O(1) entries in the upper
!>         triangle.
!>
!>    (13) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has evenly spaced entries 1, ..., ULP
!>         with random complex angles on the diagonal and random O(1)
!>         entries in the upper triangle.
!>
!>    (14) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has geometrically spaced entries
!>         1, ..., ULP with random complex angles on the diagonal
!>         and random O(1) entries in the upper triangle.
!>
!>    (15) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has "clustered" entries 1, ULP,..., ULP
!>         with random complex angles on the diagonal and random O(1)
!>         entries in the upper triangle.
!>
!>    (16) A matrix of the form  X' T X, where X has condition
!>         SQRT( ULP ) and T has complex eigenvalues randomly chosen
!>         from   ULP < |z| < 1   and random O(1) entries in the upper
!>         triangle.
!>
!>    (17) Same as (16), but multiplied by SQRT( overflow threshold )
!>    (18) Same as (16), but multiplied by SQRT( underflow threshold )
!>
!>    (19) Nonsymmetric matrix with random entries chosen from |z| < 1
!>    (20) Same as (19), but multiplied by SQRT( overflow threshold )
!>    (21) Same as (19), but multiplied by SQRT( underflow threshold )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  NSIZES - INTEGER
!>           The number of sizes of matrices to use.  If it is zero,
!>           zchkhs does nothing.  It must be at least zero.
!>           Not modified.
!>
!>  NN     - INTEGER array, dimension (NSIZES)
!>           An array containing the sizes to be used for the matrices.
!>           Zero values will be skipped.  The values must be at least
!>           zero.
!>           Not modified.
!>
!>  NTYPES - INTEGER
!>           The number of elements in DOTYPE.   If it is zero, zchkhs
!>           does nothing.  It must be at least zero.  If it is MAXTYP+1
!>           and NSIZES is 1, then an additional type, MAXTYP+1 is
!>           defined, which is to use whatever matrix is in A.  This
!>           is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>           DOTYPE(MAXTYP+1) is .TRUE. .
!>           Not modified.
!>
!>  DOTYPE - LOGICAL array, dimension (NTYPES)
!>           If DOTYPE(j) is .TRUE., then for each size in NN a
!>           matrix of that size and of type j will be generated.
!>           If NTYPES is smaller than the maximum number of types
!>           defined (PARAMETER MAXTYP), then types NTYPES+1 through
!>           MAXTYP will not be generated.  If NTYPES is larger
!>           than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
!>           will be ignored.
!>           Not modified.
!>
!>  ISEED  - INTEGER array, dimension (4)
!>           On entry ISEED specifies the seed of the random number
!>           generator. The array elements should be between 0 and 4095;
!>           if not they will be reduced mod 4096.  Also, ISEED(4) must
!>           be odd.  The random number generator uses a linear
!>           congruential sequence limited to small integers, and so
!>           should produce machine independent random numbers. The
!>           values of ISEED are changed on exit, and can be used in the
!>           next call to zchkhs to continue the same random number
!>           sequence.
!>           Modified.
!>
!>  THRESH - DOUBLE PRECISION
!>           A test will count as "failed" if the "error", computed as
!>           described above, exceeds THRESH.  Note that the error
!>           is scaled to be O(1), so THRESH should be a reasonably
!>           small multiple of 1, e.g., 10 or 100.  In particular,
!>           it should not depend on the precision (single vs. double)
!>           or the size of the matrix.  It must be at least zero.
!>           Not modified.
!>
!>  NOUNIT - INTEGER
!>           The FORTRAN unit number for printing out error messages
!>           (e.g., if a routine returns IINFO not equal to 0.)
!>           Not modified.
!>
!>  A      - COMPLEX*16 array, dimension (LDA,max(NN))
!>           Used to hold the matrix whose eigenvalues are to be
!>           computed.  On exit, A contains the last matrix actually
!>           used.
!>           Modified.
!>
!>  LDA    - INTEGER
!>           The leading dimension of A, H, T1 and T2.  It must be at
!>           least 1 and at least max( NN ).
!>           Not modified.
!>
!>  H      - COMPLEX*16 array, dimension (LDA,max(NN))
!>           The upper hessenberg matrix computed by ZGEHRD.  On exit,
!>           H contains the Hessenberg form of the matrix in A.
!>           Modified.
!>
!>  T1     - COMPLEX*16 array, dimension (LDA,max(NN))
!>           The Schur (="quasi-triangular") matrix computed by ZHSEQR
!>           if Z is computed.  On exit, T1 contains the Schur form of
!>           the matrix in A.
!>           Modified.
!>
!>  T2     - COMPLEX*16 array, dimension (LDA,max(NN))
!>           The Schur matrix computed by ZHSEQR when Z is not computed.
!>           This should be identical to T1.
!>           Modified.
!>
!>  LDU    - INTEGER
!>           The leading dimension of U, Z, UZ and UU.  It must be at
!>           least 1 and at least max( NN ).
!>           Not modified.
!>
!>  U      - COMPLEX*16 array, dimension (LDU,max(NN))
!>           The unitary matrix computed by ZGEHRD.
!>           Modified.
!>
!>  Z      - COMPLEX*16 array, dimension (LDU,max(NN))
!>           The unitary matrix computed by ZHSEQR.
!>           Modified.
!>
!>  UZ     - COMPLEX*16 array, dimension (LDU,max(NN))
!>           The product of U times Z.
!>           Modified.
!>
!>  W1     - COMPLEX*16 array, dimension (max(NN))
!>           The eigenvalues of A, as computed by a full Schur
!>           decomposition H = Z T Z'.  On exit, W1 contains the
!>           eigenvalues of the matrix in A.
!>           Modified.
!>
!>  W3     - COMPLEX*16 array, dimension (max(NN))
!>           The eigenvalues of A, as computed by a partial Schur
!>           decomposition (Z not computed, T only computed as much
!>           as is necessary for determining eigenvalues).  On exit,
!>           W3 contains the eigenvalues of the matrix in A, possibly
!>           perturbed by ZHSEIN.
!>           Modified.
!>
!>  EVECTL - COMPLEX*16 array, dimension (LDU,max(NN))
!>           The conjugate transpose of the (upper triangular) left
!>           eigenvector matrix for the matrix in T1.
!>           Modified.
!>
!>  EVEZTR - COMPLEX*16 array, dimension (LDU,max(NN))
!>           The (upper triangular) right eigenvector matrix for the
!>           matrix in T1.
!>           Modified.
!>
!>  EVECTY - COMPLEX*16 array, dimension (LDU,max(NN))
!>           The conjugate transpose of the left eigenvector matrix
!>           for the matrix in H.
!>           Modified.
!>
!>  EVECTX - COMPLEX*16 array, dimension (LDU,max(NN))
!>           The right eigenvector matrix for the matrix in H.
!>           Modified.
!>
!>  UU     - COMPLEX*16 array, dimension (LDU,max(NN))
!>           Details of the unitary matrix computed by ZGEHRD.
!>           Modified.
!>
!>  TAU    - COMPLEX*16 array, dimension (max(NN))
!>           Further details of the unitary matrix computed by ZGEHRD.
!>           Modified.
!>
!>  WORK   - COMPLEX*16 array, dimension (NWORK)
!>           Workspace.
!>           Modified.
!>
!>  NWORK  - INTEGER
!>           The number of entries in WORK.  NWORK >= 4*NN(j)*NN(j) + 2.
!>
!>  RWORK  - DOUBLE PRECISION array, dimension (max(NN))
!>           Workspace.  Could be equivalenced to IWORK, but not SELECT.
!>           Modified.
!>
!>  IWORK  - INTEGER array, dimension (max(NN))
!>           Workspace.
!>           Modified.
!>
!>  SELECT - LOGICAL array, dimension (max(NN))
!>           Workspace.  Could be equivalenced to IWORK, but not RWORK.
!>           Modified.
!>
!>  RESULT - DOUBLE PRECISION array, dimension (14)
!>           The values computed by the fourteen tests described above.
!>           The values are currently limited to 1/ulp, to avoid
!>           overflow.
!>           Modified.
!>
!>  INFO   - INTEGER
!>           If 0, then everything ran OK.
!>            -1: NSIZES < 0
!>            -2: Some NN(j) < 0
!>            -3: NTYPES < 0
!>            -6: THRESH < 0
!>            -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
!>           -14: LDU < 1 or LDU < NMAX.
!>           -26: NWORK too small.
!>           If  zlatmr, clatms, or clatme returns an error code, the
!>               absolute value of it is returned.
!>           If 1, then ZHSEQR could not find all the shifts.
!>           If 2, then the EISPACK code (for small blocks) failed.
!>           If >2, then 30*N iterations were not enough to find an
!>               eigenvalue or to decompose the problem.
!>           Modified.
!>
!>-----------------------------------------------------------------------
!>
!>     Some Local Variables and Parameters:
!>     ---- ----- --------- --- ----------
!>
!>     ZERO, ONE       Real 0 and 1.
!>     MAXTYP          The number of types defined.
!>     MTEST           The number of tests defined: care must be taken
!>                     that (1) the size of RESULT, (2) the number of
!>                     tests actually performed, and (3) MTEST agree.
!>     NTEST           The number of tests performed on this matrix
!>                     so far.  This should be less than MTEST, and
!>                     equal to it by the last test.  It will be less
!>                     if any of the routines being tested indicates
!>                     that it could not compute the matrices that
!>                     would be tested.
!>     NMAX            Largest value in NN.
!>     NMATS           The number of matrices generated so far.
!>     NERRS           The number of tests which have exceeded THRESH
!>                     so far (computed by dlafts).
!>     COND, CONDS,
!>     IMODE           Values to be passed to the matrix generators.
!>     ANORM           Norm of A; passed to matrix generators.
!>
!>     OVFL, UNFL      Overflow and underflow thresholds.
!>     ULP, ULPINV     Finest relative precision and its inverse.
!>     RTOVFL, RTUNFL,
!>     RTULP, RTULPI   Square roots of the previous 4 values.
!>
!>             The following four arrays decode JTYPE:
!>     KTYPE(j)        The general type (1-10) for type "j".
!>     KMODE(j)        The MODE value to be passed to the matrix
!>                     generator for type "j".
!>     KMAGN(j)        The order of magnitude ( O(1),
!>                     O(overflow^(1/2) ), O(underflow^(1/2) )
!>     KCONDS(j)       Selects whether CONDS is to be 1 or
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
!> \date December 2016
!
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE zchkhs(Nsizes,Nn,Ntypes,Dotype,Iseed,Thresh,Nounit,A,  &
     &                  Lda,H,T1,T2,U,Ldu,Z,Uz,W1,W3,Evectl,Evectr,     &
     &                  Evecty,Evectx,Uu,Tau,Work,Nwork,Rwork,Iwork,    &
     &                  Select,Result,Info)
      IMPLICIT NONE
!*--ZCHKHS415
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldu , Nounit , Nsizes , Ntypes , Nwork
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*) , Select(*)
      INTEGER Iseed(4) , Iwork(*) , Nn(*)
      DOUBLE PRECISION Result(14) , Rwork(*)
      COMPLEX*16 A(Lda,*) , Evectl(Ldu,*) , Evectr(Ldu,*) ,             &
     &           Evectx(Ldu,*) , Evecty(Ldu,*) , H(Lda,*) , T1(Lda,*) , &
     &           T2(Lda,*) , Tau(*) , U(Ldu,*) , Uu(Ldu,*) , Uz(Ldu,*) ,&
     &           W1(*) , W3(*) , Work(*) , Z(Ldu,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
      INTEGER MAXTYP
      PARAMETER (MAXTYP=21)
!     ..
!     .. Local Scalars ..
      LOGICAL badnn , match
      INTEGER i , ihi , iinfo , ilo , imode , in , itype , j , jcol ,   &
     &        jj , jsize , jtype , k , mtypes , n , n1 , nerrs , nmats ,&
     &        nmax , ntest , ntestt
      DOUBLE PRECISION aninv , anorm , cond , conds , ovfl , rtovfl ,   &
     &                 rtulp , rtulpi , rtunfl , temp1 , temp2 , ulp ,  &
     &                 ulpinv , unfl
!     ..
!     .. Local Arrays ..
      INTEGER idumma(1) , ioldsd(4) , kconds(MAXTYP) , kmagn(MAXTYP) ,  &
     &        kmode(MAXTYP) , ktype(MAXTYP)
      DOUBLE PRECISION dumma(4)
      COMPLEX*16 cdumma(4)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD , ZCOPY , ZGEHRD 
     EXTERNAL         ZGEMM , ZHSEIN , ZHSEQR 
     EXTERNAL         ZLACPY , ZLASET , ZTREVC 
     EXTERNAL         ZUNGHR , ZUNMHR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN , SQRT
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
!     Check for errors
!
      ntestt = 0
      Info = 0
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
      ELSEIF ( Lda<=1 .OR. Lda<nmax ) THEN
         Info = -9
      ELSEIF ( Ldu<=1 .OR. Ldu<nmax ) THEN
         Info = -14
      ELSEIF ( 4*nmax*nmax+2>Nwork ) THEN
         Info = -26
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL xerbla('zchkhs',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( Nsizes==0 .OR. Ntypes==0 ) RETURN
!
!     More important constants
!
      unfl = DLAMCH('Safe minimum')
      ovfl = DLAMCH('Overflow')
      CALL DLABAD(unfl,ovfl)
      ulp = DLAMCH('Epsilon')*DLAMCH('Base')
      ulpinv = ONE/ulp
      rtunfl = SQRT(unfl)
      rtovfl = SQRT(ovfl)
      rtulp = SQRT(ulp)
      rtulpi = ONE/rtulp
!
!     Loop over sizes, types
!
      nerrs = 0
      nmats = 0
!
      DO jsize = 1 , Nsizes
         n = Nn(jsize)
         IF ( n/=0 ) THEN
            n1 = MAX(1,n)
            aninv = ONE/DBLE(n1)
!
            IF ( Nsizes/=1 ) THEN
               mtypes = MIN(MAXTYP,Ntypes)
            ELSE
               mtypes = MIN(MAXTYP+1,Ntypes)
            ENDIF
!
            DO jtype = 1 , mtypes
               IF ( .NOT.Dotype(jtype) ) CYCLE
               nmats = nmats + 1
               ntest = 0
!
!           Save ISEED in case of an error.
!
               DO j = 1 , 4
                  ioldsd(j) = Iseed(j)
               ENDDO
!
!           Initialize RESULT
!
               DO j = 1 , 14
                  Result(j) = ZERO
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
!       =5                 random log   hermitian, w/ eigenvalues
!       =6                 random       general, w/ eigenvalues
!       =7                              random diagonal
!       =8                              random hermitian
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
                  CALL ZLASET('Full',Lda,n,CZERO,CZERO,A,Lda)
                  iinfo = 0
                  cond = ulpinv
!
!           Special Matrices
!
                  IF ( itype==1 ) THEN
!
!              Zero
!
                     iinfo = 0
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
                     CALL zlatmr(n,n,'D',Iseed,'N',Work,imode,cond,CONE,&
     &                           'T','N',Work(n+1),1,ONE,Work(2*n+1),1, &
     &                           ONE,'N',idumma,0,0,ZERO,anorm,'NO',A,  &
     &                           Lda,Iwork,iinfo)
!
                  ELSEIF ( itype==5 ) THEN
!
!              Hermitian, eigenvalues specified
!
                     CALL zlatms(n,n,'D',Iseed,'H',Rwork,imode,cond,    &
     &                           anorm,n,n,'N',A,Lda,Work,iinfo)
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
                     CALL zlatme(n,'D',Iseed,Work,imode,cond,CONE,'T',  &
     &                           'T','T',Rwork,4,conds,n,n,anorm,A,Lda, &
     &                           Work(n+1),iinfo)
!
                  ELSEIF ( itype==7 ) THEN
!
!              Diagonal, random eigenvalues
!
                     CALL zlatmr(n,n,'D',Iseed,'N',Work,6,ONE,CONE,'T', &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,0,0,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==8 ) THEN
!
!              Hermitian, random eigenvalues
!
                     CALL zlatmr(n,n,'D',Iseed,'H',Work,6,ONE,CONE,'T', &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,n,n,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==9 ) THEN
!
!              General, random eigenvalues
!
                     CALL zlatmr(n,n,'D',Iseed,'N',Work,6,ONE,CONE,'T', &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,n,n,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==10 ) THEN
!
!              Triangular, random eigenvalues
!
                     CALL zlatmr(n,n,'D',Iseed,'N',Work,6,ONE,CONE,'T', &
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
                     WRITE (Nounit,FMT=99001) 'Generator' , iinfo , n , &
     &                      jtype , ioldsd
                     Info = ABS(iinfo)
                     RETURN
                  ENDIF
               ENDIF
!
!
!           Call ZGEHRD to compute H and U, do tests.
!
               CALL ZLACPY(' ',n,n,A,Lda,H,Lda)
               ntest = 1
!
               ilo = 1
               ihi = n
!
               CALL ZGEHRD(n,ilo,ihi,H,Lda,Work,Work(n+1),Nwork-n,iinfo)
!
               IF ( iinfo/=0 ) THEN
                  Result(1) = ulpinv
                  WRITE (Nounit,FMT=99001) 'ZGEHRD' , iinfo , n ,       &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  GOTO 30
               ENDIF
!
               DO j = 1 , n - 1
                  Uu(j+1,j) = CZERO
                  DO i = j + 2 , n
                     U(i,j) = H(i,j)
                     Uu(i,j) = H(i,j)
                     H(i,j) = CZERO
                  ENDDO
               ENDDO
               CALL ZCOPY(n-1,Work,1,Tau,1)
               CALL ZUNGHR(n,ilo,ihi,U,Ldu,Work,Work(n+1),Nwork-n,iinfo)
               ntest = 2
!
               CALL zhst01(n,ilo,ihi,A,Lda,H,Lda,U,Ldu,Work,Nwork,Rwork,&
     &                     Result(1))
!
!           Call ZHSEQR to compute T1, T2 and Z, do tests.
!
!           Eigenvalues only (W3)
!
               CALL ZLACPY(' ',n,n,H,Lda,T2,Lda)
               ntest = 3
               Result(3) = ulpinv
!
               CALL ZHSEQR('E','N',n,ilo,ihi,T2,Lda,W3,Uz,Ldu,Work,     &
     &                     Nwork,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'ZHSEQR(E)' , iinfo , n ,    &
     &                   jtype , ioldsd
                  IF ( iinfo<=n+2 ) THEN
                     Info = ABS(iinfo)
                     GOTO 30
                  ENDIF
               ENDIF
!
!           Eigenvalues (W1) and Full Schur Form (T2)
!
               CALL ZLACPY(' ',n,n,H,Lda,T2,Lda)
!
               CALL ZHSEQR('S','N',n,ilo,ihi,T2,Lda,W1,Uz,Ldu,Work,     &
     &                     Nwork,iinfo)
               IF ( iinfo/=0 .AND. iinfo<=n+2 ) THEN
                  WRITE (Nounit,FMT=99001) 'ZHSEQR(S)' , iinfo , n ,    &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  GOTO 30
               ENDIF
!
!           Eigenvalues (W1), Schur Form (T1), and Schur Vectors (UZ)
!
               CALL ZLACPY(' ',n,n,H,Lda,T1,Lda)
               CALL ZLACPY(' ',n,n,U,Ldu,Uz,Ldu)
!
               CALL ZHSEQR('S','V',n,ilo,ihi,T1,Lda,W1,Uz,Ldu,Work,     &
     &                     Nwork,iinfo)
               IF ( iinfo/=0 .AND. iinfo<=n+2 ) THEN
                  WRITE (Nounit,FMT=99001) 'ZHSEQR(V)' , iinfo , n ,    &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  GOTO 30
               ENDIF
!
!           Compute Z = U' UZ
!
               CALL ZGEMM('C','N',n,n,n,CONE,U,Ldu,Uz,Ldu,CZERO,Z,Ldu)
               ntest = 8
!
!           Do Tests 3: | H - Z T Z' | / ( |H| n ulp )
!                and 4: | I - Z Z' | / ( n ulp )
!
               CALL zhst01(n,ilo,ihi,H,Lda,T1,Lda,Z,Ldu,Work,Nwork,     &
     &                     Rwork,Result(3))
!
!           Do Tests 5: | A - UZ T (UZ)' | / ( |A| n ulp )
!                and 6: | I - UZ (UZ)' | / ( n ulp )
!
               CALL zhst01(n,ilo,ihi,A,Lda,T1,Lda,Uz,Ldu,Work,Nwork,    &
     &                     Rwork,Result(5))
!
!           Do Test 7: | T2 - T1 | / ( |T| n ulp )
!
               CALL zget10(n,n,T2,Lda,T1,Lda,Work,Rwork,Result(7))
!
!           Do Test 8: | W3 - W1 | / ( max(|W1|,|W3|) ulp )
!
               temp1 = ZERO
               temp2 = ZERO
               DO j = 1 , n
                  temp1 = MAX(temp1,ABS(W1(j)),ABS(W3(j)))
                  temp2 = MAX(temp2,ABS(W1(j)-W3(j)))
               ENDDO
!
               Result(8) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!           Compute the Left and Right Eigenvectors of T
!
!           Compute the Right eigenvector Matrix:
!
               ntest = 9
               Result(9) = ulpinv
!
!           Select every other eigenvector
!
               DO j = 1 , n
                  Select(j) = .FALSE.
               ENDDO
               DO j = 1 , n , 2
                  Select(j) = .TRUE.
               ENDDO
               CALL ZTREVC('Right','All',Select,n,T1,Lda,cdumma,Ldu,    &
     &                     Evectr,Ldu,n,in,Work,Rwork,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'ZTREVC(R,A)' , iinfo , n ,  &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  GOTO 30
               ENDIF
!
!           Test 9:  | TR - RW | / ( |T| |R| ulp )
!
               CALL zget22('N','N','N',n,T1,Lda,Evectr,Ldu,W1,Work,     &
     &                     Rwork,dumma(1))
               Result(9) = dumma(1)
               IF ( dumma(2)>Thresh ) WRITE (Nounit,FMT=99002) 'Right' ,&
     &              'ZTREVC' , dumma(2) , n , jtype , ioldsd
!
!           Compute selected right eigenvectors and confirm that
!           they agree with previous right eigenvectors
!
               CALL ZTREVC('Right','Some',Select,n,T1,Lda,cdumma,Ldu,   &
     &                     Evectl,Ldu,n,in,Work,Rwork,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'ZTREVC(R,S)' , iinfo , n ,  &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  GOTO 30
               ENDIF
!
               k = 1
               match = .TRUE.
               DO j = 1 , n
                  IF ( Select(j) ) THEN
                     DO jj = 1 , n
                        IF ( Evectr(jj,j)/=Evectl(jj,k) ) THEN
                           match = .FALSE.
                           GOTO 10
                        ENDIF
                     ENDDO
                     k = k + 1
                  ENDIF
               ENDDO
 10            IF ( .NOT.match ) WRITE (Nounit,FMT=99003) 'Right' ,     &
     &                                  'ZTREVC' , n , jtype , ioldsd
!
!           Compute the Left eigenvector Matrix:
!
               ntest = 10
               Result(10) = ulpinv
               CALL ZTREVC('Left','All',Select,n,T1,Lda,Evectl,Ldu,     &
     &                     cdumma,Ldu,n,in,Work,Rwork,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'ZTREVC(L,A)' , iinfo , n ,  &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  GOTO 30
               ENDIF
!
!           Test 10:  | LT - WL | / ( |T| |L| ulp )
!
               CALL zget22('C','N','C',n,T1,Lda,Evectl,Ldu,W1,Work,     &
     &                     Rwork,dumma(3))
               Result(10) = dumma(3)
               IF ( dumma(4)>Thresh ) WRITE (Nounit,FMT=99002) 'Left' , &
     &              'ZTREVC' , dumma(4) , n , jtype , ioldsd
!
!           Compute selected left eigenvectors and confirm that
!           they agree with previous left eigenvectors
!
               CALL ZTREVC('Left','Some',Select,n,T1,Lda,Evectr,Ldu,    &
     &                     cdumma,Ldu,n,in,Work,Rwork,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'ZTREVC(L,S)' , iinfo , n ,  &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  GOTO 30
               ENDIF
!
               k = 1
               match = .TRUE.
               DO j = 1 , n
                  IF ( Select(j) ) THEN
                     DO jj = 1 , n
                        IF ( Evectl(jj,j)/=Evectr(jj,k) ) THEN
                           match = .FALSE.
                           GOTO 20
                        ENDIF
                     ENDDO
                     k = k + 1
                  ENDIF
               ENDDO
 20            IF ( .NOT.match ) WRITE (Nounit,FMT=99003) 'Left' ,      &
     &                                  'ZTREVC' , n , jtype , ioldsd
!
!           Call ZHSEIN for Right eigenvectors of H, do test 11
!
               ntest = 11
               Result(11) = ulpinv
               DO j = 1 , n
                  Select(j) = .TRUE.
               ENDDO
!
               CALL ZHSEIN('Right','Qr','Ninitv',Select,n,H,Lda,W3,     &
     &                     cdumma,Ldu,Evectx,Ldu,n1,in,Work,Rwork,Iwork,&
     &                     Iwork,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'ZHSEIN(R)' , iinfo , n ,    &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  IF ( iinfo<0 ) GOTO 30
               ELSE
!
!              Test 11:  | HX - XW | / ( |H| |X| ulp )
!
!                        (from inverse iteration)
!
                  CALL zget22('N','N','N',n,H,Lda,Evectx,Ldu,W3,Work,   &
     &                        Rwork,dumma(1))
                  IF ( dumma(1)<ulpinv ) Result(11) = dumma(1)*aninv
                  IF ( dumma(2)>Thresh ) WRITE (Nounit,FMT=99002)       &
     &                  'Right' , 'ZHSEIN' , dumma(2) , n , jtype ,     &
     &                 ioldsd
               ENDIF
!
!           Call ZHSEIN for Left eigenvectors of H, do test 12
!
               ntest = 12
               Result(12) = ulpinv
               DO j = 1 , n
                  Select(j) = .TRUE.
               ENDDO
!
               CALL ZHSEIN('Left','Qr','Ninitv',Select,n,H,Lda,W3,      &
     &                     Evecty,Ldu,cdumma,Ldu,n1,in,Work,Rwork,Iwork,&
     &                     Iwork,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'ZHSEIN(L)' , iinfo , n ,    &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  IF ( iinfo<0 ) GOTO 30
               ELSE
!
!              Test 12:  | YH - WY | / ( |H| |Y| ulp )
!
!                        (from inverse iteration)
!
                  CALL zget22('C','N','C',n,H,Lda,Evecty,Ldu,W3,Work,   &
     &                        Rwork,dumma(3))
                  IF ( dumma(3)<ulpinv ) Result(12) = dumma(3)*aninv
                  IF ( dumma(4)>Thresh ) WRITE (Nounit,FMT=99002)       &
     &                 'Left' , 'ZHSEIN' , dumma(4) , n , jtype , ioldsd
               ENDIF
!
!           Call ZUNMHR for Right eigenvectors of A, do test 13
!
               ntest = 13
               Result(13) = ulpinv
!
               CALL ZUNMHR('Left','No transpose',n,n,ilo,ihi,Uu,Ldu,Tau,&
     &                     Evectx,Ldu,Work,Nwork,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'ZUNMHR(L)' , iinfo , n ,    &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  IF ( iinfo<0 ) GOTO 30
               ELSE
!
!              Test 13:  | AX - XW | / ( |A| |X| ulp )
!
!                        (from inverse iteration)
!
                  CALL zget22('N','N','N',n,A,Lda,Evectx,Ldu,W3,Work,   &
     &                        Rwork,dumma(1))
                  IF ( dumma(1)<ulpinv ) Result(13) = dumma(1)*aninv
               ENDIF
!
!           Call ZUNMHR for Left eigenvectors of A, do test 14
!
               ntest = 14
               Result(14) = ulpinv
!
               CALL ZUNMHR('Left','No transpose',n,n,ilo,ihi,Uu,Ldu,Tau,&
     &                     Evecty,Ldu,Work,Nwork,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'ZUNMHR(L)' , iinfo , n ,    &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  IF ( iinfo<0 ) THEN
                  ENDIF
               ELSE
!
!              Test 14:  | YA - WY | / ( |A| |Y| ulp )
!
!                        (from inverse iteration)
!
                  CALL zget22('C','N','C',n,A,Lda,Evecty,Ldu,W3,Work,   &
     &                        Rwork,dumma(3))
                  IF ( dumma(3)<ulpinv ) Result(14) = dumma(3)*aninv
               ENDIF
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
!
 30            ntestt = ntestt + ntest
               CALL dlafts('ZHS',n,n,jtype,ntest,Result,ioldsd,Thresh,  &
     &                     Nounit,nerrs)
!
            ENDDO
         ENDIF
      ENDDO
!
!     Summary
!
      CALL dlasum('ZHS',Nounit,nerrs,ntestt)
!
      RETURN
!
99001 FORMAT (' zchkhs: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
99002 FORMAT (' zchkhs: ',A,' Eigenvectors from ',A,' incorrectly ',    &
     &        'normalized.',/' Bits of error=',0P,G10.3,',',9X,'N=',I6, &
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
99003 FORMAT (' zchkhs: Selected ',A,' Eigenvectors from ',A,           &
     &        ' do not match other eigenvectors ',9X,'N=',I6,', JTYPE=',&
     &        I6,', ISEED=(',3(I5,','),I5,')')
!
!     End of zchkhs
!
      END SUBROUTINE zchkhs
