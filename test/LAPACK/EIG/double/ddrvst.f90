!*==ddrvst.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DDRVST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DDRVST( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          NOUNIT, A, LDA, D1, D2, D3, D4, EVEIGS, WA1,
!                          WA2, WA3, U, LDU, V, TAU, Z, WORK, LWORK,
!                          IWORK, LIWORK, RESULT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDU, LIWORK, LWORK, NOUNIT, NSIZES,
!      $                   NTYPES
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
!       DOUBLE PRECISION   A( LDA, * ), D1( * ), D2( * ), D3( * ),
!      $                   D4( * ), EVEIGS( * ), RESULT( * ), TAU( * ),
!      $                   U( LDU, * ), V( LDU, * ), WA1( * ), WA2( * ),
!      $                   WA3( * ), WORK( * ), Z( LDU, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      DDRVST  checks the symmetric eigenvalue problem drivers.
!>
!>              DSTEV computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric tridiagonal matrix.
!>
!>              DSTEVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric tridiagonal matrix.
!>
!>              DSTEVR computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric tridiagonal matrix
!>              using the Relatively Robust Representation where it can.
!>
!>              DSYEV computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix.
!>
!>              DSYEVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix.
!>
!>              DSYEVR computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix
!>              using the Relatively Robust Representation where it can.
!>
!>              DSPEV computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix in packed
!>              storage.
!>
!>              DSPEVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix in packed
!>              storage.
!>
!>              DSBEV computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric band matrix.
!>
!>              DSBEVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric band matrix.
!>
!>              DSYEVD computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix using
!>              a divide and conquer algorithm.
!>
!>              DSPEVD computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric matrix in packed
!>              storage, using a divide and conquer algorithm.
!>
!>              DSBEVD computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric band matrix,
!>              using a divide and conquer algorithm.
!>
!>      When DDRVST is called, a number of matrix "sizes" ("n's") and a
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
!>      (3)  A diagonal matrix with evenly spaced eigenvalues
!>           1, ..., ULP  and random signs.
!>           (ULP = (first number larger than 1) - 1 )
!>      (4)  A diagonal matrix with geometrically spaced eigenvalues
!>           1, ..., ULP  and random signs.
!>      (5)  A diagonal matrix with "clustered" eigenvalues
!>           1, ULP, ..., ULP and random signs.
!>
!>      (6)  Same as (4), but multiplied by SQRT( overflow threshold )
!>      (7)  Same as (4), but multiplied by SQRT( underflow threshold )
!>
!>      (8)  A matrix of the form  U' D U, where U is orthogonal and
!>           D has evenly spaced entries 1, ..., ULP with random signs
!>           on the diagonal.
!>
!>      (9)  A matrix of the form  U' D U, where U is orthogonal and
!>           D has geometrically spaced entries 1, ..., ULP with random
!>           signs on the diagonal.
!>
!>      (10) A matrix of the form  U' D U, where U is orthogonal and
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
!>          DDRVST does nothing.  It must be at least zero.
!>          Not modified.
!>
!>  NN      INTEGER array, dimension (NSIZES)
!>          An array containing the sizes to be used for the matrices.
!>          Zero values will be skipped.  The values must be at least
!>          zero.
!>          Not modified.
!>
!>  NTYPES  INTEGER
!>          The number of elements in DOTYPE.   If it is zero, DDRVST
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
!>          next call to DDRVST to continue the same random number
!>          sequence.
!>          Modified.
!>
!>  THRESH  DOUBLE PRECISION
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
!>  A       DOUBLE PRECISION array, dimension (LDA , max(NN))
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
!>  D1      DOUBLE PRECISION array, dimension (max(NN))
!>          The eigenvalues of A, as computed by DSTEQR simlutaneously
!>          with Z.  On exit, the eigenvalues in D1 correspond with the
!>          matrix in A.
!>          Modified.
!>
!>  D2      DOUBLE PRECISION array, dimension (max(NN))
!>          The eigenvalues of A, as computed by DSTEQR if Z is not
!>          computed.  On exit, the eigenvalues in D2 correspond with
!>          the matrix in A.
!>          Modified.
!>
!>  D3      DOUBLE PRECISION array, dimension (max(NN))
!>          The eigenvalues of A, as computed by DSTERF.  On exit, the
!>          eigenvalues in D3 correspond with the matrix in A.
!>          Modified.
!>
!>  D4      DOUBLE PRECISION array, dimension
!>
!>  EVEIGS  DOUBLE PRECISION array, dimension (max(NN))
!>          The eigenvalues as computed by DSTEV('N', ... )
!>          (I reserve the right to change this to the output of
!>          whichever algorithm computes the most accurate eigenvalues).
!>
!>  WA1     DOUBLE PRECISION array, dimension
!>
!>  WA2     DOUBLE PRECISION array, dimension
!>
!>  WA3     DOUBLE PRECISION array, dimension
!>
!>  U       DOUBLE PRECISION array, dimension (LDU, max(NN))
!>          The orthogonal matrix computed by DSYTRD + DORGTR.
!>          Modified.
!>
!>  LDU     INTEGER
!>          The leading dimension of U, Z, and V.  It must be at
!>          least 1 and at least max( NN ).
!>          Not modified.
!>
!>  V       DOUBLE PRECISION array, dimension (LDU, max(NN))
!>          The Housholder vectors computed by DSYTRD in reducing A to
!>          tridiagonal form.
!>          Modified.
!>
!>  TAU     DOUBLE PRECISION array, dimension (max(NN))
!>          The Householder factors computed by DSYTRD in reducing A
!>          to tridiagonal form.
!>          Modified.
!>
!>  Z       DOUBLE PRECISION array, dimension (LDU, max(NN))
!>          The orthogonal matrix of eigenvectors computed by DSTEQR,
!>          DPTEQR, and DSTEIN.
!>          Modified.
!>
!>  WORK    DOUBLE PRECISION array, dimension (LWORK)
!>          Workspace.
!>          Modified.
!>
!>  LWORK   INTEGER
!>          The number of entries in WORK.  This must be at least
!>          1 + 4 * Nmax + 2 * Nmax * lg Nmax + 4 * Nmax**2
!>          where Nmax = max( NN(j), 2 ) and lg = log base 2.
!>          Not modified.
!>
!>  IWORK   INTEGER array,
!>             dimension (6 + 6*Nmax + 5 * Nmax * lg Nmax )
!>          where Nmax = max( NN(j), 2 ) and lg = log base 2.
!>          Workspace.
!>          Modified.
!>
!>  RESULT  DOUBLE PRECISION array, dimension (105)
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
!>          If  DLATMR, DLATMS, DSYTRD, DORGTR, DSTEQR, DSTERF,
!>              or DORMTR returns an error code, the
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
!>                       so far (computed by DLAFTS).
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
!>
!>     The tests performed are:                 Routine tested
!>    1= | A - U S U' | / ( |A| n ulp )         DSTEV('V', ... )
!>    2= | I - U U' | / ( n ulp )               DSTEV('V', ... )
!>    3= |D(with Z) - D(w/o Z)| / (|D| ulp)     DSTEV('N', ... )
!>    4= | A - U S U' | / ( |A| n ulp )         DSTEVX('V','A', ... )
!>    5= | I - U U' | / ( n ulp )               DSTEVX('V','A', ... )
!>    6= |D(with Z) - EVEIGS| / (|D| ulp)       DSTEVX('N','A', ... )
!>    7= | A - U S U' | / ( |A| n ulp )         DSTEVR('V','A', ... )
!>    8= | I - U U' | / ( n ulp )               DSTEVR('V','A', ... )
!>    9= |D(with Z) - EVEIGS| / (|D| ulp)       DSTEVR('N','A', ... )
!>    10= | A - U S U' | / ( |A| n ulp )        DSTEVX('V','I', ... )
!>    11= | I - U U' | / ( n ulp )              DSTEVX('V','I', ... )
!>    12= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSTEVX('N','I', ... )
!>    13= | A - U S U' | / ( |A| n ulp )        DSTEVX('V','V', ... )
!>    14= | I - U U' | / ( n ulp )              DSTEVX('V','V', ... )
!>    15= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSTEVX('N','V', ... )
!>    16= | A - U S U' | / ( |A| n ulp )        DSTEVD('V', ... )
!>    17= | I - U U' | / ( n ulp )              DSTEVD('V', ... )
!>    18= |D(with Z) - EVEIGS| / (|D| ulp)      DSTEVD('N', ... )
!>    19= | A - U S U' | / ( |A| n ulp )        DSTEVR('V','I', ... )
!>    20= | I - U U' | / ( n ulp )              DSTEVR('V','I', ... )
!>    21= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSTEVR('N','I', ... )
!>    22= | A - U S U' | / ( |A| n ulp )        DSTEVR('V','V', ... )
!>    23= | I - U U' | / ( n ulp )              DSTEVR('V','V', ... )
!>    24= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSTEVR('N','V', ... )
!>
!>    25= | A - U S U' | / ( |A| n ulp )        DSYEV('L','V', ... )
!>    26= | I - U U' | / ( n ulp )              DSYEV('L','V', ... )
!>    27= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEV('L','N', ... )
!>    28= | A - U S U' | / ( |A| n ulp )        DSYEVX('L','V','A', ... )
!>    29= | I - U U' | / ( n ulp )              DSYEVX('L','V','A', ... )
!>    30= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVX('L','N','A', ... )
!>    31= | A - U S U' | / ( |A| n ulp )        DSYEVX('L','V','I', ... )
!>    32= | I - U U' | / ( n ulp )              DSYEVX('L','V','I', ... )
!>    33= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVX('L','N','I', ... )
!>    34= | A - U S U' | / ( |A| n ulp )        DSYEVX('L','V','V', ... )
!>    35= | I - U U' | / ( n ulp )              DSYEVX('L','V','V', ... )
!>    36= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVX('L','N','V', ... )
!>    37= | A - U S U' | / ( |A| n ulp )        DSPEV('L','V', ... )
!>    38= | I - U U' | / ( n ulp )              DSPEV('L','V', ... )
!>    39= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEV('L','N', ... )
!>    40= | A - U S U' | / ( |A| n ulp )        DSPEVX('L','V','A', ... )
!>    41= | I - U U' | / ( n ulp )              DSPEVX('L','V','A', ... )
!>    42= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVX('L','N','A', ... )
!>    43= | A - U S U' | / ( |A| n ulp )        DSPEVX('L','V','I', ... )
!>    44= | I - U U' | / ( n ulp )              DSPEVX('L','V','I', ... )
!>    45= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVX('L','N','I', ... )
!>    46= | A - U S U' | / ( |A| n ulp )        DSPEVX('L','V','V', ... )
!>    47= | I - U U' | / ( n ulp )              DSPEVX('L','V','V', ... )
!>    48= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVX('L','N','V', ... )
!>    49= | A - U S U' | / ( |A| n ulp )        DSBEV('L','V', ... )
!>    50= | I - U U' | / ( n ulp )              DSBEV('L','V', ... )
!>    51= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEV('L','N', ... )
!>    52= | A - U S U' | / ( |A| n ulp )        DSBEVX('L','V','A', ... )
!>    53= | I - U U' | / ( n ulp )              DSBEVX('L','V','A', ... )
!>    54= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVX('L','N','A', ... )
!>    55= | A - U S U' | / ( |A| n ulp )        DSBEVX('L','V','I', ... )
!>    56= | I - U U' | / ( n ulp )              DSBEVX('L','V','I', ... )
!>    57= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVX('L','N','I', ... )
!>    58= | A - U S U' | / ( |A| n ulp )        DSBEVX('L','V','V', ... )
!>    59= | I - U U' | / ( n ulp )              DSBEVX('L','V','V', ... )
!>    60= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVX('L','N','V', ... )
!>    61= | A - U S U' | / ( |A| n ulp )        DSYEVD('L','V', ... )
!>    62= | I - U U' | / ( n ulp )              DSYEVD('L','V', ... )
!>    63= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVD('L','N', ... )
!>    64= | A - U S U' | / ( |A| n ulp )        DSPEVD('L','V', ... )
!>    65= | I - U U' | / ( n ulp )              DSPEVD('L','V', ... )
!>    66= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVD('L','N', ... )
!>    67= | A - U S U' | / ( |A| n ulp )        DSBEVD('L','V', ... )
!>    68= | I - U U' | / ( n ulp )              DSBEVD('L','V', ... )
!>    69= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVD('L','N', ... )
!>    70= | A - U S U' | / ( |A| n ulp )        DSYEVR('L','V','A', ... )
!>    71= | I - U U' | / ( n ulp )              DSYEVR('L','V','A', ... )
!>    72= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVR('L','N','A', ... )
!>    73= | A - U S U' | / ( |A| n ulp )        DSYEVR('L','V','I', ... )
!>    74= | I - U U' | / ( n ulp )              DSYEVR('L','V','I', ... )
!>    75= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVR('L','N','I', ... )
!>    76= | A - U S U' | / ( |A| n ulp )        DSYEVR('L','V','V', ... )
!>    77= | I - U U' | / ( n ulp )              DSYEVR('L','V','V', ... )
!>    78= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSYEVR('L','N','V', ... )
!>
!>    Tests 25 through 78 are repeated (as tests 79 through 132)
!>    with UPLO='U'
!>
!>    To be added in 1999
!>
!>    79= | A - U S U' | / ( |A| n ulp )        DSPEVR('L','V','A', ... )
!>    80= | I - U U' | / ( n ulp )              DSPEVR('L','V','A', ... )
!>    81= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVR('L','N','A', ... )
!>    82= | A - U S U' | / ( |A| n ulp )        DSPEVR('L','V','I', ... )
!>    83= | I - U U' | / ( n ulp )              DSPEVR('L','V','I', ... )
!>    84= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVR('L','N','I', ... )
!>    85= | A - U S U' | / ( |A| n ulp )        DSPEVR('L','V','V', ... )
!>    86= | I - U U' | / ( n ulp )              DSPEVR('L','V','V', ... )
!>    87= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSPEVR('L','N','V', ... )
!>    88= | A - U S U' | / ( |A| n ulp )        DSBEVR('L','V','A', ... )
!>    89= | I - U U' | / ( n ulp )              DSBEVR('L','V','A', ... )
!>    90= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVR('L','N','A', ... )
!>    91= | A - U S U' | / ( |A| n ulp )        DSBEVR('L','V','I', ... )
!>    92= | I - U U' | / ( n ulp )              DSBEVR('L','V','I', ... )
!>    93= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVR('L','N','I', ... )
!>    94= | A - U S U' | / ( |A| n ulp )        DSBEVR('L','V','V', ... )
!>    95= | I - U U' | / ( n ulp )              DSBEVR('L','V','V', ... )
!>    96= |D(with Z) - D(w/o Z)| / (|D| ulp)    DSBEVR('L','N','V', ... )
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
      SUBROUTINE DDRVST(Nsizes,Nn,Ntypes,Dotype,Iseed,Thresh,Nounit,A,  &
     &                  Lda,D1,D2,D3,D4,Eveigs,Wa1,Wa2,Wa3,U,Ldu,V,Tau, &
     &                  Z,Work,Lwork,Iwork,Liwork,Result,Info)
      IMPLICIT NONE
!*--DDRVST456
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldu , Liwork , Lwork , Nounit , Nsizes ,     &
     &        Ntypes
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iseed(4) , Iwork(*) , Nn(*)
      DOUBLE PRECISION A(Lda,*) , D1(*) , D2(*) , D3(*) , D4(*) ,       &
     &                 Eveigs(*) , Result(*) , Tau(*) , U(Ldu,*) ,      &
     &                 V(Ldu,*) , Wa1(*) , Wa2(*) , Wa3(*) , Work(*) ,  &
     &                 Z(Ldu,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TWO , TEN
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,TEN=10.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      INTEGER MAXTYP
      PARAMETER (MAXTYP=18)
!     ..
!     .. Local Scalars ..
      LOGICAL badnn
      CHARACTER uplo
      INTEGER i , idiag , ihbw , iinfo , il , imode , indx , irow ,     &
     &        itemp , itype , iu , iuplo , j , j1 , j2 , jcol , jsize , &
     &        jtype , kd , lgn , liwedc , lwedc , m , m2 , m3 , mtypes ,&
     &        n , nerrs , nmats , nmax , ntest , ntestt
      DOUBLE PRECISION abstol , aninv , anorm , cond , ovfl , rtovfl ,  &
     &                 rtunfl , temp1 , temp2 , temp3 , ulp , ulpinv ,  &
     &                 unfl , vl , vu
!     ..
!     .. Local Arrays ..
      INTEGER idumma(1) , ioldsd(4) , iseed2(4) , iseed3(4) ,           &
     &        kmagn(MAXTYP) , kmode(MAXTYP) , ktype(MAXTYP)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLARND , DSXT1
      EXTERNAL DLAMCH , DLARND , DSXT1
!     ..
!     .. External Subroutines ..
      EXTERNAL ALASVM , DLABAD , DLACPY , DLAFTS , DLASET , DLATMR ,    &
     &         DLATMS , DSBEV , DSBEVD , DSBEVX , DSPEV , DSPEVD ,      &
     &         DSPEVX , DSTEV , DSTEVD , DSTEVR , DSTEVX , DSTT21 ,     &
     &         DSTT22 , DSYEV , DSYEVD , DSYEVR , DSYEVX , DSYT21 ,     &
     &         DSYT22 , XERBLA
!     ..
!     .. Scalars in Common ..
      CHARACTER*32 SRNamt
!     ..
!     .. Common blocks ..
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , INT , LOG , MAX , MIN , SQRT
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
!     Keep ftrnchek happy
!
      vl = ZERO
      vu = ZERO
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
         Info = -21
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DDRVST',-Info)
         RETURN
      ENDIF
!
!     Quick return if nothing to do
!
      IF ( Nsizes==0 .OR. Ntypes==0 ) RETURN
!
!     More Important constants
!
      unfl = DLAMCH('Safe minimum')
      ovfl = DLAMCH('Overflow')
      CALL DLABAD(unfl,ovfl)
      ulp = DLAMCH('Epsilon')*DLAMCH('Base')
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
!
      DO jsize = 1 , Nsizes
         n = Nn(jsize)
         IF ( n>0 ) THEN
            lgn = INT(LOG(DBLE(n))/LOG(TWO))
            IF ( 2**lgn<n ) lgn = lgn + 1
            IF ( 2**lgn<n ) lgn = lgn + 1
            lwedc = 1 + 4*n + 2*n*lgn + 4*n**2
!           LIWEDC = 6 + 6*N + 5*N*LGN
            liwedc = 3 + 5*n
         ELSE
            lwedc = 9
!           LIWEDC = 12
            liwedc = 8
         ENDIF
         aninv = ONE/DBLE(MAX(1,n))
!
         IF ( Nsizes/=1 ) THEN
            mtypes = MIN(MAXTYP,Ntypes)
         ELSE
            mtypes = MIN(MAXTYP+1,Ntypes)
         ENDIF
!
         DO jtype = 1 , mtypes
!
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
!           =5         random log   symmetric, w/ eigenvalues
!           =6         random       (none)
!           =7                      random diagonal
!           =8                      random symmetric
!           =9                      band symmetric, w/ eigenvalues
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
                  CALL DLASET('Full',Lda,n,ZERO,ZERO,A,Lda)
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
                     CALL DLATMS(n,n,'S',Iseed,'S',Work,imode,cond,     &
     &                           anorm,0,0,'N',A,Lda,Work(n+1),iinfo)
!
                  ELSEIF ( itype==5 ) THEN
!
!              Symmetric, eigenvalues specified
!
                     CALL DLATMS(n,n,'S',Iseed,'S',Work,imode,cond,     &
     &                           anorm,n,n,'N',A,Lda,Work(n+1),iinfo)
!
                  ELSEIF ( itype==7 ) THEN
!
!              Diagonal, random eigenvalues
!
                     idumma(1) = 1
                     CALL DLATMR(n,n,'S',Iseed,'S',Work,6,ONE,ONE,'T',  &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,0,0,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==8 ) THEN
!
!              Symmetric, random eigenvalues
!
                     idumma(1) = 1
                     CALL DLATMR(n,n,'S',Iseed,'S',Work,6,ONE,ONE,'T',  &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,n,n,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==9 ) THEN
!
!              Symmetric banded, eigenvalues specified
!
                     ihbw = INT((n-1)*DLARND(1,iseed3))
                     CALL DLATMS(n,n,'S',Iseed,'S',Work,imode,cond,     &
     &                           anorm,ihbw,ihbw,'Z',U,Ldu,Work(n+1),   &
     &                           iinfo)
!
!              Store as dense matrix for most routines.
!
                     CALL DLASET('Full',Lda,n,ZERO,ZERO,A,Lda)
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
                  il = 1 + (n-1)*INT(DLARND(1,iseed2))
                  iu = 1 + (n-1)*INT(DLARND(1,iseed2))
                  IF ( il>iu ) THEN
                     itemp = il
                     il = iu
                     iu = itemp
                  ENDIF
               ENDIF
!
!           3)      If matrix is tridiagonal, call DSTEV and DSTEVX.
!
               IF ( jtype<=7 ) THEN
                  ntest = 1
                  DO i = 1 , n
                     D1(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D2(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEV'
                  CALL DSTEV('V',n,D1,D2,Z,Ldu,Work,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEV(V)' , iinfo , n ,  &
     &                      jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(1) = ulpinv
                        Result(2) = ulpinv
                        Result(3) = ulpinv
                        GOTO 5
                     ENDIF
                  ENDIF
!
!              Do tests 1 and 2.
!
                  DO i = 1 , n
                     D3(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  CALL DSTT21(n,0,D3,D4,D1,D2,Z,Ldu,Work,Result(1))
!
                  ntest = 3
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEV'
                  CALL DSTEV('N',n,D3,D4,Z,Ldu,Work,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEV(N)' , iinfo , n ,  &
     &                      jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(3) = ulpinv
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
                  Result(3) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!
 5                ntest = 4
                  DO i = 1 , n
                     Eveigs(i) = D3(i)
                     D1(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D2(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVX'
                  CALL DSTEVX('V','A',n,D1,D2,vl,vu,il,iu,abstol,m,Wa1, &
     &                        Z,Ldu,Work,Iwork,Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVX(V,A)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(4) = ulpinv
                        Result(5) = ulpinv
                        Result(6) = ulpinv
                        GOTO 10
                     ENDIF
                  ENDIF
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
!
!              Do tests 4 and 5.
!
                  DO i = 1 , n
                     D3(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  CALL DSTT21(n,0,D3,D4,Wa1,D2,Z,Ldu,Work,Result(4))
!
                  ntest = 6
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVX'
                  CALL DSTEVX('N','A',n,D3,D4,vl,vu,il,iu,abstol,m2,Wa2,&
     &                        Z,Ldu,Work,Iwork,Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVX(N,A)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(6) = ulpinv
                        GOTO 10
                     ENDIF
                  ENDIF
!
!              Do test 6.
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(Wa2(j)),ABS(Eveigs(j)))
                     temp2 = MAX(temp2,ABS(Wa2(j)-Eveigs(j)))
                  ENDDO
                  Result(6) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!
 10               ntest = 7
                  DO i = 1 , n
                     D1(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D2(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVR'
                  CALL DSTEVR('V','A',n,D1,D2,vl,vu,il,iu,abstol,m,Wa1, &
     &                        Z,Ldu,Iwork,Work,Lwork,Iwork(2*n+1),      &
     &                        Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVR(V,A)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(7) = ulpinv
                        Result(8) = ulpinv
                        GOTO 15
                     ENDIF
                  ENDIF
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
!
!              Do tests 7 and 8.
!
                  DO i = 1 , n
                     D3(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  CALL DSTT21(n,0,D3,D4,Wa1,D2,Z,Ldu,Work,Result(7))
!
                  ntest = 9
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVR'
                  CALL DSTEVR('N','A',n,D3,D4,vl,vu,il,iu,abstol,m2,Wa2,&
     &                        Z,Ldu,Iwork,Work,Lwork,Iwork(2*n+1),      &
     &                        Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVR(N,A)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(9) = ulpinv
                        GOTO 15
                     ENDIF
                  ENDIF
!
!              Do test 9.
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(Wa2(j)),ABS(Eveigs(j)))
                     temp2 = MAX(temp2,ABS(Wa2(j)-Eveigs(j)))
                  ENDDO
                  Result(9) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!
!
 15               ntest = 10
                  DO i = 1 , n
                     D1(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D2(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVX'
                  CALL DSTEVX('V','I',n,D1,D2,vl,vu,il,iu,abstol,m2,Wa2,&
     &                        Z,Ldu,Work,Iwork,Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVX(V,I)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(10) = ulpinv
                        Result(11) = ulpinv
                        Result(12) = ulpinv
                        GOTO 20
                     ENDIF
                  ENDIF
!
!              Do tests 10 and 11.
!
                  DO i = 1 , n
                     D3(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  CALL DSTT22(n,m2,0,D3,D4,Wa2,D2,Z,Ldu,Work,MAX(1,m2), &
     &                        Result(10))
!
!
                  ntest = 12
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVX'
                  CALL DSTEVX('N','I',n,D3,D4,vl,vu,il,iu,abstol,m3,Wa3,&
     &                        Z,Ldu,Work,Iwork,Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVX(N,I)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(12) = ulpinv
                        GOTO 20
                     ENDIF
                  ENDIF
!
!              Do test 12.
!
                  temp1 = DSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = DSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  Result(12) = (temp1+temp2)/MAX(unfl,ulp*temp3)
!
!
 20               ntest = 12
                  IF ( n>0 ) THEN
                     IF ( il/=1 ) THEN
                        vl = Wa1(il)                                    &
     &                       - MAX(HALF*(Wa1(il)-Wa1(il-1)),TEN*ulp*    &
     &                       temp3,TEN*rtunfl)
                     ELSE
                        vl = Wa1(1)                                     &
     &                       - MAX(HALF*(Wa1(n)-Wa1(1)),TEN*ulp*temp3,  &
     &                       TEN*rtunfl)
                     ENDIF
                     IF ( iu/=n ) THEN
                        vu = Wa1(iu)                                    &
     &                       + MAX(HALF*(Wa1(iu+1)-Wa1(iu)),TEN*ulp*    &
     &                       temp3,TEN*rtunfl)
                     ELSE
                        vu = Wa1(n)                                     &
     &                       + MAX(HALF*(Wa1(n)-Wa1(1)),TEN*ulp*temp3,  &
     &                       TEN*rtunfl)
                     ENDIF
                  ELSE
                     vl = ZERO
                     vu = ONE
                  ENDIF
!
                  DO i = 1 , n
                     D1(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D2(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVX'
                  CALL DSTEVX('V','V',n,D1,D2,vl,vu,il,iu,abstol,m2,Wa2,&
     &                        Z,Ldu,Work,Iwork,Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVX(V,V)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(13) = ulpinv
                        Result(14) = ulpinv
                        Result(15) = ulpinv
                        GOTO 25
                     ENDIF
                  ENDIF
!
                  IF ( m2==0 .AND. n>0 ) THEN
                     Result(13) = ulpinv
                     Result(14) = ulpinv
                     Result(15) = ulpinv
                     GOTO 25
                  ENDIF
!
!              Do tests 13 and 14.
!
                  DO i = 1 , n
                     D3(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  CALL DSTT22(n,m2,0,D3,D4,Wa2,D2,Z,Ldu,Work,MAX(1,m2), &
     &                        Result(13))
!
                  ntest = 15
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVX'
                  CALL DSTEVX('N','V',n,D3,D4,vl,vu,il,iu,abstol,m3,Wa3,&
     &                        Z,Ldu,Work,Iwork,Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVX(N,V)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(15) = ulpinv
                        GOTO 25
                     ENDIF
                  ENDIF
!
!              Do test 15.
!
                  temp1 = DSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = DSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  Result(15) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
!
 25               ntest = 16
                  DO i = 1 , n
                     D1(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D2(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVD'
                  CALL DSTEVD('V',n,D1,D2,Z,Ldu,Work,lwedc,Iwork,liwedc,&
     &                        iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVD(V)' , iinfo , n , &
     &                      jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(16) = ulpinv
                        Result(17) = ulpinv
                        Result(18) = ulpinv
                        GOTO 30
                     ENDIF
                  ENDIF
!
!              Do tests 16 and 17.
!
                  DO i = 1 , n
                     D3(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  CALL DSTT21(n,0,D3,D4,D1,D2,Z,Ldu,Work,Result(16))
!
                  ntest = 18
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVD'
                  CALL DSTEVD('N',n,D3,D4,Z,Ldu,Work,lwedc,Iwork,liwedc,&
     &                        iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVD(N)' , iinfo , n , &
     &                      jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(18) = ulpinv
                        GOTO 30
                     ENDIF
                  ENDIF
!
!              Do test 18.
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(Eveigs(j)),ABS(D3(j)))
                     temp2 = MAX(temp2,ABS(Eveigs(j)-D3(j)))
                  ENDDO
                  Result(18) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!
 30               ntest = 19
                  DO i = 1 , n
                     D1(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D2(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVR'
                  CALL DSTEVR('V','I',n,D1,D2,vl,vu,il,iu,abstol,m2,Wa2,&
     &                        Z,Ldu,Iwork,Work,Lwork,Iwork(2*n+1),      &
     &                        Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVR(V,I)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(19) = ulpinv
                        Result(20) = ulpinv
                        Result(21) = ulpinv
                        GOTO 35
                     ENDIF
                  ENDIF
!
!              DO tests 19 and 20.
!
                  DO i = 1 , n
                     D3(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  CALL DSTT22(n,m2,0,D3,D4,Wa2,D2,Z,Ldu,Work,MAX(1,m2), &
     &                        Result(19))
!
!
                  ntest = 21
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVR'
                  CALL DSTEVR('N','I',n,D3,D4,vl,vu,il,iu,abstol,m3,Wa3,&
     &                        Z,Ldu,Iwork,Work,Lwork,Iwork(2*n+1),      &
     &                        Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVR(N,I)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(21) = ulpinv
                        GOTO 35
                     ENDIF
                  ENDIF
!
!              Do test 21.
!
                  temp1 = DSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = DSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  Result(21) = (temp1+temp2)/MAX(unfl,ulp*temp3)
!
!
 35               ntest = 21
                  IF ( n>0 ) THEN
                     IF ( il/=1 ) THEN
                        vl = Wa1(il)                                    &
     &                       - MAX(HALF*(Wa1(il)-Wa1(il-1)),TEN*ulp*    &
     &                       temp3,TEN*rtunfl)
                     ELSE
                        vl = Wa1(1)                                     &
     &                       - MAX(HALF*(Wa1(n)-Wa1(1)),TEN*ulp*temp3,  &
     &                       TEN*rtunfl)
                     ENDIF
                     IF ( iu/=n ) THEN
                        vu = Wa1(iu)                                    &
     &                       + MAX(HALF*(Wa1(iu+1)-Wa1(iu)),TEN*ulp*    &
     &                       temp3,TEN*rtunfl)
                     ELSE
                        vu = Wa1(n)                                     &
     &                       + MAX(HALF*(Wa1(n)-Wa1(1)),TEN*ulp*temp3,  &
     &                       TEN*rtunfl)
                     ENDIF
                  ELSE
                     vl = ZERO
                     vu = ONE
                  ENDIF
!
                  DO i = 1 , n
                     D1(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D2(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVR'
                  CALL DSTEVR('V','V',n,D1,D2,vl,vu,il,iu,abstol,m2,Wa2,&
     &                        Z,Ldu,Iwork,Work,Lwork,Iwork(2*n+1),      &
     &                        Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVR(V,V)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(22) = ulpinv
                        Result(23) = ulpinv
                        Result(24) = ulpinv
                        GOTO 40
                     ENDIF
                  ENDIF
!
                  IF ( m2==0 .AND. n>0 ) THEN
                     Result(22) = ulpinv
                     Result(23) = ulpinv
                     Result(24) = ulpinv
                     GOTO 40
                  ENDIF
!
!              Do tests 22 and 23.
!
                  DO i = 1 , n
                     D3(i) = DBLE(A(i,i))
                  ENDDO
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  CALL DSTT22(n,m2,0,D3,D4,Wa2,D2,Z,Ldu,Work,MAX(1,m2), &
     &                        Result(22))
!
                  ntest = 24
                  DO i = 1 , n - 1
                     D4(i) = DBLE(A(i+1,i))
                  ENDDO
                  SRNamt = 'DSTEVR'
                  CALL DSTEVR('N','V',n,D3,D4,vl,vu,il,iu,abstol,m3,Wa3,&
     &                        Z,Ldu,Iwork,Work,Lwork,Iwork(2*n+1),      &
     &                        Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSTEVR(N,V)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(24) = ulpinv
                        GOTO 40
                     ENDIF
                  ENDIF
!
!              Do test 24.
!
                  temp1 = DSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = DSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  Result(24) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
!
!
!
               ELSE
!
                  DO i = 1 , 24
                     Result(i) = ZERO
                  ENDDO
                  ntest = 24
               ENDIF
!
!           Perform remaining tests storing upper or lower triangular
!           part of matrix.
!
 40            DO iuplo = 0 , 1
                  IF ( iuplo==0 ) THEN
                     uplo = 'L'
                  ELSE
                     uplo = 'U'
                  ENDIF
!
!              4)      Call DSYEV and DSYEVX.
!
                  CALL DLACPY(' ',n,n,A,Lda,V,Ldu)
!
                  ntest = ntest + 1
                  SRNamt = 'DSYEV'
                  CALL DSYEV('V',uplo,n,A,Ldu,D1,Work,Lwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEV(V,'//uplo//')' ,   &
     &                      iinfo , n , jtype , ioldsd
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
!              Do tests 25 and 26 (or +54)
!
                  CALL DSYT21(1,uplo,n,0,V,Ldu,D1,D2,A,Ldu,Z,Ldu,Tau,   &
     &                        Work,Result(ntest))
!
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  ntest = ntest + 2
                  SRNamt = 'DSYEV'
                  CALL DSYEV('N',uplo,n,A,Ldu,D3,Work,Lwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEV(N,'//uplo//')' ,   &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 45
                     ENDIF
                  ENDIF
!
!              Do test 27 (or +54)
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(D1(j)),ABS(D3(j)))
                     temp2 = MAX(temp2,ABS(D1(j)-D3(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
 45               CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
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
                  SRNamt = 'DSYEVX'
                  CALL DSYEVX('V','A',uplo,n,A,Ldu,vl,vu,il,iu,abstol,m,&
     &                        Wa1,Z,Ldu,Work,Lwork,Iwork,Iwork(5*n+1),  &
     &                        iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVX(V,A,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
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
!              Do tests 28 and 29 (or +54)
!
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  CALL DSYT21(1,uplo,n,0,A,Ldu,D1,D2,Z,Ldu,V,Ldu,Tau,   &
     &                        Work,Result(ntest))
!
                  ntest = ntest + 2
                  SRNamt = 'DSYEVX'
                  CALL DSYEVX('N','A',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m2,Wa2,Z,Ldu,Work,Lwork,Iwork,Iwork(5*n+1)&
     &                        ,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVX(N,A,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 50
                     ENDIF
                  ENDIF
!
!              Do test 30 (or +54)
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
 50               ntest = ntest + 1
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
                  SRNamt = 'DSYEVX'
                  CALL DSYEVX('V','I',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m2,Wa2,Z,Ldu,Work,Lwork,Iwork,Iwork(5*n+1)&
     &                        ,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVX(V,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
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
!              Do tests 31 and 32 (or +54)
!
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  CALL DSYT22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Result(ntest))
!
                  ntest = ntest + 2
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
                  SRNamt = 'DSYEVX'
                  CALL DSYEVX('N','I',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m3,Wa3,Z,Ldu,Work,Lwork,Iwork,Iwork(5*n+1)&
     &                        ,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVX(N,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 55
                     ENDIF
                  ENDIF
!
!              Do test 33 (or +54)
!
                  temp1 = DSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = DSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  Result(ntest) = (temp1+temp2)/MAX(unfl,ulp*temp3)
!
 55               ntest = ntest + 1
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
                  SRNamt = 'DSYEVX'
                  CALL DSYEVX('V','V',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m2,Wa2,Z,Ldu,Work,Lwork,Iwork,Iwork(5*n+1)&
     &                        ,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVX(V,V,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
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
!              Do tests 34 and 35 (or +54)
!
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  CALL DSYT22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Result(ntest))
!
                  ntest = ntest + 2
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
                  SRNamt = 'DSYEVX'
                  CALL DSYEVX('N','V',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m3,Wa3,Z,Ldu,Work,Lwork,Iwork,Iwork(5*n+1)&
     &                        ,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVX(N,V,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
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
!              Do test 36 (or +54)
!
                  temp1 = DSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = DSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
!
!              5)      Call DSPEV and DSPEVX.
!
 60               CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
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
                  SRNamt = 'DSPEV'
                  CALL DSPEV('V',uplo,n,Work,D1,Z,Ldu,V,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSPEV(V,'//uplo//')' ,   &
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
!              Do tests 37 and 38 (or +54)
!
                  CALL DSYT21(1,uplo,n,0,A,Lda,D1,D2,Z,Ldu,V,Ldu,Tau,   &
     &                        Work,Result(ntest))
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
                  SRNamt = 'DSPEV'
                  CALL DSPEV('N',uplo,n,Work,D3,Z,Ldu,V,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSPEV(N,'//uplo//')' ,   &
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
!              Do test 39 (or +54)
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
 65               IF ( iuplo==1 ) THEN
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
                  SRNamt = 'DSPEVX'
                  CALL DSPEVX('V','A',uplo,n,Work,vl,vu,il,iu,abstol,m, &
     &                        Wa1,Z,Ldu,V,Iwork,Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSPEVX(V,A,'//uplo//')' ,&
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
!              Do tests 40 and 41 (or +54)
!
                  CALL DSYT21(1,uplo,n,0,A,Ldu,Wa1,D2,Z,Ldu,V,Ldu,Tau,  &
     &                        Work,Result(ntest))
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
                  SRNamt = 'DSPEVX'
                  CALL DSPEVX('N','A',uplo,n,Work,vl,vu,il,iu,abstol,m2,&
     &                        Wa2,Z,Ldu,V,Iwork,Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSPEVX(N,A,'//uplo//')' ,&
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
!              Do test 42 (or +54)
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(Wa1(j)),ABS(Wa2(j)))
                     temp2 = MAX(temp2,ABS(Wa1(j)-Wa2(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
 70               IF ( iuplo==1 ) THEN
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
                  SRNamt = 'DSPEVX'
                  CALL DSPEVX('V','I',uplo,n,Work,vl,vu,il,iu,abstol,m2,&
     &                        Wa2,Z,Ldu,V,Iwork,Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSPEVX(V,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
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
!              Do tests 43 and 44 (or +54)
!
                  CALL DSYT22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Result(ntest))
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
                  SRNamt = 'DSPEVX'
                  CALL DSPEVX('N','I',uplo,n,Work,vl,vu,il,iu,abstol,m3,&
     &                        Wa3,Z,Ldu,V,Iwork,Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSPEVX(N,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 75
                     ENDIF
                  ENDIF
!
                  IF ( m3==0 .AND. n>0 ) THEN
                     Result(ntest) = ulpinv
                     GOTO 75
                  ENDIF
!
!              Do test 45 (or +54)
!
                  temp1 = DSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = DSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
 75               IF ( iuplo==1 ) THEN
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
                  SRNamt = 'DSPEVX'
                  CALL DSPEVX('V','V',uplo,n,Work,vl,vu,il,iu,abstol,m2,&
     &                        Wa2,Z,Ldu,V,Iwork,Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSPEVX(V,V,'//uplo//')' ,&
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
!              Do tests 46 and 47 (or +54)
!
                  CALL DSYT22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Result(ntest))
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
                  SRNamt = 'DSPEVX'
                  CALL DSPEVX('N','V',uplo,n,Work,vl,vu,il,iu,abstol,m3,&
     &                        Wa3,Z,Ldu,V,Iwork,Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSPEVX(N,V,'//uplo//')' ,&
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
                  IF ( m3==0 .AND. n>0 ) THEN
                     Result(ntest) = ulpinv
                     GOTO 80
                  ENDIF
!
!              Do test 48 (or +54)
!
                  temp1 = DSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = DSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
!
!              6)      Call DSBEV and DSBEVX.
!
 80               IF ( jtype<=7 ) THEN
                     kd = 1
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
                  SRNamt = 'DSBEV'
                  CALL DSBEV('V',uplo,n,kd,V,Ldu,D1,Z,Ldu,Work,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSBEV(V,'//uplo//')' ,   &
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
!              Do tests 49 and 50 (or ... )
!
                  CALL DSYT21(1,uplo,n,0,A,Lda,D1,D2,Z,Ldu,V,Ldu,Tau,   &
     &                        Work,Result(ntest))
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
                  SRNamt = 'DSBEV'
                  CALL DSBEV('N',uplo,n,kd,V,Ldu,D3,Z,Ldu,Work,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSBEV(N,'//uplo//')' ,   &
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
!              Do test 51 (or +54)
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
 85               IF ( iuplo==1 ) THEN
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
                  SRNamt = 'DSBEVX'
                  CALL DSBEVX('V','A',uplo,n,kd,V,Ldu,U,Ldu,vl,vu,il,iu,&
     &                        abstol,m,Wa2,Z,Ldu,Work,Iwork,Iwork(5*n+1)&
     &                        ,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSBEVX(V,A,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 90
                     ENDIF
                  ENDIF
!
!              Do tests 52 and 53 (or +54)
!
                  CALL DSYT21(1,uplo,n,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,Tau,  &
     &                        Work,Result(ntest))
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
                  SRNamt = 'DSBEVX'
                  CALL DSBEVX('N','A',uplo,n,kd,V,Ldu,U,Ldu,vl,vu,il,iu,&
     &                        abstol,m3,Wa3,Z,Ldu,Work,Iwork,           &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSBEVX(N,A,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 90
                     ENDIF
                  ENDIF
!
!              Do test 54 (or +54)
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(Wa2(j)),ABS(Wa3(j)))
                     temp2 = MAX(temp2,ABS(Wa2(j)-Wa3(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
 90               ntest = ntest + 1
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
                  SRNamt = 'DSBEVX'
                  CALL DSBEVX('V','I',uplo,n,kd,V,Ldu,U,Ldu,vl,vu,il,iu,&
     &                        abstol,m2,Wa2,Z,Ldu,Work,Iwork,           &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSBEVX(V,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 95
                     ENDIF
                  ENDIF
!
!              Do tests 55 and 56 (or +54)
!
                  CALL DSYT22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Result(ntest))
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
                  SRNamt = 'DSBEVX'
                  CALL DSBEVX('N','I',uplo,n,kd,V,Ldu,U,Ldu,vl,vu,il,iu,&
     &                        abstol,m3,Wa3,Z,Ldu,Work,Iwork,           &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSBEVX(N,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 95
                     ENDIF
                  ENDIF
!
!              Do test 57 (or +54)
!
                  temp1 = DSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = DSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
 95               ntest = ntest + 1
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
                  SRNamt = 'DSBEVX'
                  CALL DSBEVX('V','V',uplo,n,kd,V,Ldu,U,Ldu,vl,vu,il,iu,&
     &                        abstol,m2,Wa2,Z,Ldu,Work,Iwork,           &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSBEVX(V,V,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 100
                     ENDIF
                  ENDIF
!
!              Do tests 58 and 59 (or +54)
!
                  CALL DSYT22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Result(ntest))
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
                  SRNamt = 'DSBEVX'
                  CALL DSBEVX('N','V',uplo,n,kd,V,Ldu,U,Ldu,vl,vu,il,iu,&
     &                        abstol,m3,Wa3,Z,Ldu,Work,Iwork,           &
     &                        Iwork(5*n+1),iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSBEVX(N,V,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 100
                     ENDIF
                  ENDIF
!
                  IF ( m3==0 .AND. n>0 ) THEN
                     Result(ntest) = ulpinv
                     GOTO 100
                  ENDIF
!
!              Do test 60 (or +54)
!
                  temp1 = DSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = DSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
!
!              7)      Call DSYEVD
!
 100              CALL DLACPY(' ',n,n,A,Lda,V,Ldu)
!
                  ntest = ntest + 1
                  SRNamt = 'DSYEVD'
                  CALL DSYEVD('V',uplo,n,A,Ldu,D1,Work,lwedc,Iwork,     &
     &                        liwedc,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVD(V,'//uplo//')' ,  &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 105
                     ENDIF
                  ENDIF
!
!              Do tests 61 and 62 (or +54)
!
                  CALL DSYT21(1,uplo,n,0,V,Ldu,D1,D2,A,Ldu,Z,Ldu,Tau,   &
     &                        Work,Result(ntest))
!
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  ntest = ntest + 2
                  SRNamt = 'DSYEVD'
                  CALL DSYEVD('N',uplo,n,A,Ldu,D3,Work,lwedc,Iwork,     &
     &                        liwedc,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVD(N,'//uplo//')' ,  &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 105
                     ENDIF
                  ENDIF
!
!              Do test 63 (or +54)
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
!              8)      Call DSPEVD.
!
 105              CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
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
                  SRNamt = 'DSPEVD'
                  CALL DSPEVD('V',uplo,n,Work,D1,Z,Ldu,Work(indx),      &
     &                        lwedc-indx+1,Iwork,liwedc,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSPEVD(V,'//uplo//')' ,  &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 110
                     ENDIF
                  ENDIF
!
!              Do tests 64 and 65 (or +54)
!
                  CALL DSYT21(1,uplo,n,0,A,Lda,D1,D2,Z,Ldu,V,Ldu,Tau,   &
     &                        Work,Result(ntest))
!
                  IF ( iuplo==1 ) THEN
                     indx = 1
                     DO j = 1 , n
                        DO i = 1 , j
!
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
                  SRNamt = 'DSPEVD'
                  CALL DSPEVD('N',uplo,n,Work,D3,Z,Ldu,Work(indx),      &
     &                        lwedc-indx+1,Iwork,liwedc,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSPEVD(N,'//uplo//')' ,  &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 110
                     ENDIF
                  ENDIF
!
!              Do test 66 (or +54)
!
                  temp1 = ZERO
                  temp2 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(D1(j)),ABS(D3(j)))
                     temp2 = MAX(temp2,ABS(D1(j)-D3(j)))
                  ENDDO
                  Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!              9)      Call DSBEVD.
!
 110              IF ( jtype<=7 ) THEN
                     kd = 1
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
                  SRNamt = 'DSBEVD'
                  CALL DSBEVD('V',uplo,n,kd,V,Ldu,D1,Z,Ldu,Work,lwedc,  &
     &                        Iwork,liwedc,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSBEVD(V,'//uplo//')' ,  &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 115
                     ENDIF
                  ENDIF
!
!              Do tests 67 and 68 (or +54)
!
                  CALL DSYT21(1,uplo,n,0,A,Lda,D1,D2,Z,Ldu,V,Ldu,Tau,   &
     &                        Work,Result(ntest))
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
                  SRNamt = 'DSBEVD'
                  CALL DSBEVD('N',uplo,n,kd,V,Ldu,D3,Z,Ldu,Work,lwedc,  &
     &                        Iwork,liwedc,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSBEVD(N,'//uplo//')' ,  &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 115
                     ENDIF
                  ENDIF
!
!              Do test 69 (or +54)
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
!
 115              CALL DLACPY(' ',n,n,A,Lda,V,Ldu)
                  ntest = ntest + 1
                  SRNamt = 'DSYEVR'
                  CALL DSYEVR('V','A',uplo,n,A,Ldu,vl,vu,il,iu,abstol,m,&
     &                        Wa1,Z,Ldu,Iwork,Work,Lwork,Iwork(2*n+1),  &
     &                        Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVR(V,A,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 120
                     ENDIF
                  ENDIF
!
!              Do tests 70 and 71 (or ... )
!
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  CALL DSYT21(1,uplo,n,0,A,Ldu,Wa1,D2,Z,Ldu,V,Ldu,Tau,  &
     &                        Work,Result(ntest))
!
                  ntest = ntest + 2
                  SRNamt = 'DSYEVR'
                  CALL DSYEVR('N','A',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m2,Wa2,Z,Ldu,Iwork,Work,Lwork,Iwork(2*n+1)&
     &                        ,Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVR(N,A,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 120
                     ENDIF
                  ENDIF
!
!              Do test 72 (or ... )
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
 120              ntest = ntest + 1
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
                  SRNamt = 'DSYEVR'
                  CALL DSYEVR('V','I',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m2,Wa2,Z,Ldu,Iwork,Work,Lwork,Iwork(2*n+1)&
     &                        ,Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVR(V,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        Result(ntest+1) = ulpinv
                        Result(ntest+2) = ulpinv
                        GOTO 125
                     ENDIF
                  ENDIF
!
!              Do tests 73 and 74 (or +54)
!
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  CALL DSYT22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Result(ntest))
!
                  ntest = ntest + 2
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
                  SRNamt = 'DSYEVR'
                  CALL DSYEVR('N','I',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m3,Wa3,Z,Ldu,Iwork,Work,Lwork,Iwork(2*n+1)&
     &                        ,Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVR(N,I,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(ntest) = ulpinv
                        GOTO 125
                     ENDIF
                  ENDIF
!
!              Do test 75 (or +54)
!
                  temp1 = DSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = DSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  Result(ntest) = (temp1+temp2)/MAX(unfl,ulp*temp3)
!
 125              ntest = ntest + 1
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
                  SRNamt = 'DSYEVR'
                  CALL DSYEVR('V','V',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m2,Wa2,Z,Ldu,Iwork,Work,Lwork,Iwork(2*n+1)&
     &                        ,Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVR(V,V,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
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
!              Do tests 76 and 77 (or +54)
!
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
!
                  CALL DSYT22(1,uplo,n,m2,0,A,Ldu,Wa2,D2,Z,Ldu,V,Ldu,   &
     &                        Tau,Work,Result(ntest))
!
                  ntest = ntest + 2
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
                  SRNamt = 'DSYEVR'
                  CALL DSYEVR('N','V',uplo,n,A,Ldu,vl,vu,il,iu,abstol,  &
     &                        m3,Wa3,Z,Ldu,Iwork,Work,Lwork,Iwork(2*n+1)&
     &                        ,Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'DSYEVR(N,V,'//uplo//')' ,&
     &                      iinfo , n , jtype , ioldsd
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
!              Do test 78 (or +54)
!
                  temp1 = DSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
                  temp2 = DSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
                  IF ( n>0 ) THEN
                     temp3 = MAX(ABS(Wa1(1)),ABS(Wa1(n)))
                  ELSE
                     temp3 = ZERO
                  ENDIF
                  Result(ntest) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
                  CALL DLACPY(' ',n,n,V,Ldu,A,Lda)
!
               ENDDO
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
               ntestt = ntestt + ntest
!
               CALL DLAFTS('DST',n,n,jtype,ntest,Result,ioldsd,Thresh,  &
     &                     Nounit,nerrs)
            ENDIF
!
         ENDDO
      ENDDO
!
!     Summary
!
      CALL ALASVM('DST',Nounit,nerrs,ntestt,0)
!
99001 FORMAT (' DDRVST: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
!
!     End of DDRVST
!
      END SUBROUTINE DDRVST
