!*==cchkst.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CCHKST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKST( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          NOUNIT, A, LDA, AP, SD, SE, D1, D2, D3, D4, D5,
!                          WA1, WA2, WA3, WR, U, LDU, V, VP, TAU, Z, WORK,
!                          LWORK, RWORK, LRWORK, IWORK, LIWORK, RESULT,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDU, LIWORK, LRWORK, LWORK, NOUNIT,
!      $                   NSIZES, NTYPES
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
!       REAL               D1( * ), D2( * ), D3( * ), D4( * ), D5( * ),
!      $                   RESULT( * ), RWORK( * ), SD( * ), SE( * ),
!      $                   WA1( * ), WA2( * ), WA3( * ), WR( * )
!       COMPLEX            A( LDA, * ), AP( * ), TAU( * ), U( LDU, * ),
!      $                   V( LDU, * ), VP( * ), WORK( * ), Z( LDU, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCHKST  checks the Hermitian eigenvalue problem routines.
!>
!>    CHETRD factors A as  U S U* , where * means conjugate transpose,
!>    S is real symmetric tridiagonal, and U is unitary.
!>    CHETRD can use either just the lower or just the upper triangle
!>    of A; CCHKST checks both cases.
!>    U is represented as a product of Householder
!>    transformations, whose vectors are stored in the first
!>    n-1 columns of V, and whose scale factors are in TAU.
!>
!>    CHPTRD does the same as CHETRD, except that A and V are stored
!>    in "packed" format.
!>
!>    CUNGTR constructs the matrix U from the contents of V and TAU.
!>
!>    CUPGTR constructs the matrix U from the contents of VP and TAU.
!>
!>    CSTEQR factors S as  Z D1 Z* , where Z is the unitary
!>    matrix of eigenvectors and D1 is a diagonal matrix with
!>    the eigenvalues on the diagonal.  D2 is the matrix of
!>    eigenvalues computed when Z is not computed.
!>
!>    SSTERF computes D3, the matrix of eigenvalues, by the
!>    PWK method, which does not yield eigenvectors.
!>
!>    CPTEQR factors S as  Z4 D4 Z4* , for a
!>    Hermitian positive definite tridiagonal matrix.
!>    D5 is the matrix of eigenvalues computed when Z is not
!>    computed.
!>
!>    SSTEBZ computes selected eigenvalues.  WA1, WA2, and
!>    WA3 will denote eigenvalues computed to high
!>    absolute accuracy, with different range options.
!>    WR will denote eigenvalues computed to high relative
!>    accuracy.
!>
!>    CSTEIN computes Y, the eigenvectors of S, given the
!>    eigenvalues.
!>
!>    CSTEDC factors S as Z D1 Z* , where Z is the unitary
!>    matrix of eigenvectors and D1 is a diagonal matrix with
!>    the eigenvalues on the diagonal ('I' option). It may also
!>    update an input unitary matrix, usually the output
!>    from CHETRD/CUNGTR or CHPTRD/CUPGTR ('V' option). It may
!>    also just compute eigenvalues ('N' option).
!>
!>    CSTEMR factors S as Z D1 Z* , where Z is the unitary
!>    matrix of eigenvectors and D1 is a diagonal matrix with
!>    the eigenvalues on the diagonal ('I' option).  CSTEMR
!>    uses the Relatively Robust Representation whenever possible.
!>
!> When CCHKST is called, a number of matrix "sizes" ("n's") and a
!> number of matrix "types" are specified.  For each size ("n")
!> and each type of matrix, one matrix will be generated and used
!> to test the Hermitian eigenroutines.  For each matrix, a number
!> of tests will be performed:
!>
!> (1)     | A - V S V* | / ( |A| n ulp ) CHETRD( UPLO='U', ... )
!>
!> (2)     | I - UV* | / ( n ulp )        CUNGTR( UPLO='U', ... )
!>
!> (3)     | A - V S V* | / ( |A| n ulp ) CHETRD( UPLO='L', ... )
!>
!> (4)     | I - UV* | / ( n ulp )        CUNGTR( UPLO='L', ... )
!>
!> (5-8)   Same as 1-4, but for CHPTRD and CUPGTR.
!>
!> (9)     | S - Z D Z* | / ( |S| n ulp ) CSTEQR('V',...)
!>
!> (10)    | I - ZZ* | / ( n ulp )        CSTEQR('V',...)
!>
!> (11)    | D1 - D2 | / ( |D1| ulp )        CSTEQR('N',...)
!>
!> (12)    | D1 - D3 | / ( |D1| ulp )        SSTERF
!>
!> (13)    0 if the true eigenvalues (computed by sturm count)
!>         of S are within THRESH of
!>         those in D1.  2*THRESH if they are not.  (Tested using
!>         SSTECH)
!>
!> For S positive definite,
!>
!> (14)    | S - Z4 D4 Z4* | / ( |S| n ulp ) CPTEQR('V',...)
!>
!> (15)    | I - Z4 Z4* | / ( n ulp )        CPTEQR('V',...)
!>
!> (16)    | D4 - D5 | / ( 100 |D4| ulp )       CPTEQR('N',...)
!>
!> When S is also diagonally dominant by the factor gamma < 1,
!>
!> (17)    max | D4(i) - WR(i) | / ( |D4(i)| omega ) ,
!>          i
!>         omega = 2 (2n-1) ULP (1 + 8 gamma**2) / (1 - gamma)**4
!>                                              SSTEBZ( 'A', 'E', ...)
!>
!> (18)    | WA1 - D3 | / ( |D3| ulp )          SSTEBZ( 'A', 'E', ...)
!>
!> (19)    ( max { min | WA2(i)-WA3(j) | } +
!>            i     j
!>           max { min | WA3(i)-WA2(j) | } ) / ( |D3| ulp )
!>            i     j
!>                                              SSTEBZ( 'I', 'E', ...)
!>
!> (20)    | S - Y WA1 Y* | / ( |S| n ulp )  SSTEBZ, CSTEIN
!>
!> (21)    | I - Y Y* | / ( n ulp )          SSTEBZ, CSTEIN
!>
!> (22)    | S - Z D Z* | / ( |S| n ulp )    CSTEDC('I')
!>
!> (23)    | I - ZZ* | / ( n ulp )           CSTEDC('I')
!>
!> (24)    | S - Z D Z* | / ( |S| n ulp )    CSTEDC('V')
!>
!> (25)    | I - ZZ* | / ( n ulp )           CSTEDC('V')
!>
!> (26)    | D1 - D2 | / ( |D1| ulp )           CSTEDC('V') and
!>                                              CSTEDC('N')
!>
!> Test 27 is disabled at the moment because CSTEMR does not
!> guarantee high relatvie accuracy.
!>
!> (27)    max | D6(i) - WR(i) | / ( |D6(i)| omega ) ,
!>          i
!>         omega = 2 (2n-1) ULP (1 + 8 gamma**2) / (1 - gamma)**4
!>                                              CSTEMR('V', 'A')
!>
!> (28)    max | D6(i) - WR(i) | / ( |D6(i)| omega ) ,
!>          i
!>         omega = 2 (2n-1) ULP (1 + 8 gamma**2) / (1 - gamma)**4
!>                                              CSTEMR('V', 'I')
!>
!> Tests 29 through 34 are disable at present because CSTEMR
!> does not handle partial spectrum requests.
!>
!> (29)    | S - Z D Z* | / ( |S| n ulp )    CSTEMR('V', 'I')
!>
!> (30)    | I - ZZ* | / ( n ulp )           CSTEMR('V', 'I')
!>
!> (31)    ( max { min | WA2(i)-WA3(j) | } +
!>            i     j
!>           max { min | WA3(i)-WA2(j) | } ) / ( |D3| ulp )
!>            i     j
!>         CSTEMR('N', 'I') vs. CSTEMR('V', 'I')
!>
!> (32)    | S - Z D Z* | / ( |S| n ulp )    CSTEMR('V', 'V')
!>
!> (33)    | I - ZZ* | / ( n ulp )           CSTEMR('V', 'V')
!>
!> (34)    ( max { min | WA2(i)-WA3(j) | } +
!>            i     j
!>           max { min | WA3(i)-WA2(j) | } ) / ( |D3| ulp )
!>            i     j
!>         CSTEMR('N', 'V') vs. CSTEMR('V', 'V')
!>
!> (35)    | S - Z D Z* | / ( |S| n ulp )    CSTEMR('V', 'A')
!>
!> (36)    | I - ZZ* | / ( n ulp )           CSTEMR('V', 'A')
!>
!> (37)    ( max { min | WA2(i)-WA3(j) | } +
!>            i     j
!>           max { min | WA3(i)-WA2(j) | } ) / ( |D3| ulp )
!>            i     j
!>         CSTEMR('N', 'A') vs. CSTEMR('V', 'A')
!>
!> The "sizes" are specified by an array NN(1:NSIZES); the value of
!> each element NN(j) specifies one size.
!> The "types" are specified by a logical array DOTYPE( 1:NTYPES );
!> if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!> Currently, the list of possible types is:
!>
!> (1)  The zero matrix.
!> (2)  The identity matrix.
!>
!> (3)  A diagonal matrix with evenly spaced entries
!>      1, ..., ULP  and random signs.
!>      (ULP = (first number larger than 1) - 1 )
!> (4)  A diagonal matrix with geometrically spaced entries
!>      1, ..., ULP  and random signs.
!> (5)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
!>      and random signs.
!>
!> (6)  Same as (4), but multiplied by SQRT( overflow threshold )
!> (7)  Same as (4), but multiplied by SQRT( underflow threshold )
!>
!> (8)  A matrix of the form  U* D U, where U is unitary and
!>      D has evenly spaced entries 1, ..., ULP with random signs
!>      on the diagonal.
!>
!> (9)  A matrix of the form  U* D U, where U is unitary and
!>      D has geometrically spaced entries 1, ..., ULP with random
!>      signs on the diagonal.
!>
!> (10) A matrix of the form  U* D U, where U is unitary and
!>      D has "clustered" entries 1, ULP,..., ULP with random
!>      signs on the diagonal.
!>
!> (11) Same as (8), but multiplied by SQRT( overflow threshold )
!> (12) Same as (8), but multiplied by SQRT( underflow threshold )
!>
!> (13) Hermitian matrix with random entries chosen from (-1,1).
!> (14) Same as (13), but multiplied by SQRT( overflow threshold )
!> (15) Same as (13), but multiplied by SQRT( underflow threshold )
!> (16) Same as (8), but diagonal elements are all positive.
!> (17) Same as (9), but diagonal elements are all positive.
!> (18) Same as (10), but diagonal elements are all positive.
!> (19) Same as (16), but multiplied by SQRT( overflow threshold )
!> (20) Same as (16), but multiplied by SQRT( underflow threshold )
!> (21) A diagonally dominant tridiagonal matrix with geometrically
!>      spaced diagonal entries 1, ..., ULP.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          CCHKST does nothing.  It must be at least zero.
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
!>          The number of elements in DOTYPE.   If it is zero, CCHKST
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
!>          next call to CCHKST to continue the same random number
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
!>          (e.g., if a routine returns IINFO not equal to 0.)
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array of
!>                                  dimension ( LDA , max(NN) )
!>          Used to hold the matrix whose eigenvalues are to be
!>          computed.  On exit, A contains the last matrix actually
!>          used.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at
!>          least 1 and at least max( NN ).
!> \endverbatim
!>
!> \param[out] AP
!> \verbatim
!>          AP is COMPLEX array of
!>                      dimension( max(NN)*max(NN+1)/2 )
!>          The matrix A stored in packed format.
!> \endverbatim
!>
!> \param[out] SD
!> \verbatim
!>          SD is REAL array of
!>                             dimension( max(NN) )
!>          The diagonal of the tridiagonal matrix computed by CHETRD.
!>          On exit, SD and SE contain the tridiagonal form of the
!>          matrix in A.
!> \endverbatim
!>
!> \param[out] SE
!> \verbatim
!>          SE is REAL array of
!>                             dimension( max(NN) )
!>          The off-diagonal of the tridiagonal matrix computed by
!>          CHETRD.  On exit, SD and SE contain the tridiagonal form of
!>          the matrix in A.
!> \endverbatim
!>
!> \param[out] D1
!> \verbatim
!>          D1 is REAL array of
!>                             dimension( max(NN) )
!>          The eigenvalues of A, as computed by CSTEQR simlutaneously
!>          with Z.  On exit, the eigenvalues in D1 correspond with the
!>          matrix in A.
!> \endverbatim
!>
!> \param[out] D2
!> \verbatim
!>          D2 is REAL array of
!>                             dimension( max(NN) )
!>          The eigenvalues of A, as computed by CSTEQR if Z is not
!>          computed.  On exit, the eigenvalues in D2 correspond with
!>          the matrix in A.
!> \endverbatim
!>
!> \param[out] D3
!> \verbatim
!>          D3 is REAL array of
!>                             dimension( max(NN) )
!>          The eigenvalues of A, as computed by SSTERF.  On exit, the
!>          eigenvalues in D3 correspond with the matrix in A.
!> \endverbatim
!>
!> \param[out] D4
!> \verbatim
!>          D4 is REAL array of
!>                             dimension( max(NN) )
!>          The eigenvalues of A, as computed by CPTEQR(V).
!>          ZPTEQR factors S as  Z4 D4 Z4*
!>          On exit, the eigenvalues in D4 correspond with the matrix in A.
!> \endverbatim
!>
!> \param[out] D5
!> \verbatim
!>          D5 is REAL array of
!>                             dimension( max(NN) )
!>          The eigenvalues of A, as computed by ZPTEQR(N)
!>          when Z is not computed. On exit, the
!>          eigenvalues in D4 correspond with the matrix in A.
!> \endverbatim
!>
!> \param[out] WA1
!> \verbatim
!>          WA1 is REAL array of
!>                             dimension( max(NN) )
!>          All eigenvalues of A, computed to high
!>          absolute accuracy, with different range options.
!>          as computed by SSTEBZ.
!> \endverbatim
!>
!> \param[out] WA2
!> \verbatim
!>          WA2 is REAL array of
!>                             dimension( max(NN) )
!>          Selected eigenvalues of A, computed to high
!>          absolute accuracy, with different range options.
!>          as computed by SSTEBZ.
!>          Choose random values for IL and IU, and ask for the
!>          IL-th through IU-th eigenvalues.
!> \endverbatim
!>
!> \param[out] WA3
!> \verbatim
!>          WA3 is REAL array of
!>                             dimension( max(NN) )
!>          Selected eigenvalues of A, computed to high
!>          absolute accuracy, with different range options.
!>          as computed by SSTEBZ.
!>          Determine the values VL and VU of the IL-th and IU-th
!>          eigenvalues and ask for all eigenvalues in this range.
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array of
!>                             dimension( max(NN) )
!>          All eigenvalues of A, computed to high
!>          absolute accuracy, with different options.
!>          as computed by DSTEBZ.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX array of
!>                             dimension( LDU, max(NN) ).
!>          The unitary matrix computed by CHETRD + CUNGTR.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U, Z, and V.  It must be at least 1
!>          and at least max( NN ).
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX array of
!>                             dimension( LDU, max(NN) ).
!>          The Housholder vectors computed by CHETRD in reducing A to
!>          tridiagonal form.  The vectors computed with UPLO='U' are
!>          in the upper triangle, and the vectors computed with UPLO='L'
!>          are in the lower triangle.  (As described in CHETRD, the
!>          sub- and superdiagonal are not set to 1, although the
!>          true Householder vector has a 1 in that position.  The
!>          routines that use V, such as CUNGTR, set those entries to
!>          1 before using them, and then restore them later.)
!> \endverbatim
!>
!> \param[out] VP
!> \verbatim
!>          VP is COMPLEX array of
!>                      dimension( max(NN)*max(NN+1)/2 )
!>          The matrix V stored in packed format.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array of
!>                             dimension( max(NN) )
!>          The Householder factors computed by CHETRD in reducing A
!>          to tridiagonal form.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX array of
!>                             dimension( LDU, max(NN) ).
!>          The unitary matrix of eigenvectors computed by CSTEQR,
!>          CPTEQR, and CSTEIN.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array of
!>                      dimension( LWORK )
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          1 + 4 * Nmax + 2 * Nmax * lg Nmax + 3 * Nmax**2
!>          where Nmax = max( NN(j), 2 ) and lg = log base 2.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array,
!>          Workspace.
!> \endverbatim
!>
!> \param[out] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The number of entries in IWORK.  This must be at least
!>                  6 + 6*Nmax + 5 * Nmax * lg Nmax
!>          where Nmax = max( NN(j), 2 ) and lg = log base 2.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The number of entries in LRWORK (dimension( ??? )
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (26)
!>          The values computed by the tests described above.
!>          The values are currently limited to 1/ulp, to avoid
!>          overflow.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0, then everything ran OK.
!>           -1: NSIZES < 0
!>           -2: Some NN(j) < 0
!>           -3: NTYPES < 0
!>           -5: THRESH < 0
!>           -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
!>          -23: LDU < 1 or LDU < NMAX.
!>          -29: LWORK too small.
!>          If  CLATMR, CLATMS, CHETRD, CUNGTR, CSTEQR, SSTERF,
!>              or CUNMC2 returns an error code, the
!>              absolute value of it is returned.
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
!>       NBLOCK          Blocksize as returned by ENVIR.
!>       NMAX            Largest value in NN.
!>       NMATS           The number of matrices generated so far.
!>       NERRS           The number of tests which have exceeded THRESH
!>                       so far.
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
      SUBROUTINE CCHKST(Nsizes,Nn,Ntypes,Dotype,Iseed,Thresh,Nounit,A,  &
     &                  Lda,Ap,Sd,Se,D1,D2,D3,D4,D5,Wa1,Wa2,Wa3,Wr,U,   &
     &                  Ldu,V,Vp,Tau,Z,Work,Lwork,Rwork,Lrwork,Iwork,   &
     &                  Liwork,Result,Info)
      IMPLICIT NONE
!*--CCHKST607
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
      REAL D1(*) , D2(*) , D3(*) , D4(*) , D5(*) , Result(*) , Rwork(*) &
     &     , Sd(*) , Se(*) , Wa1(*) , Wa2(*) , Wa3(*) , Wr(*)
      COMPLEX A(Lda,*) , Ap(*) , Tau(*) , U(Ldu,*) , V(Ldu,*) , Vp(*) , &
     &        Work(*) , Z(Ldu,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO , EIGHT , TEN , HUN
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,TWO=2.0E0,EIGHT=8.0E0,TEN=10.0E0, &
     &           HUN=100.0E0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
      REAL HALF
      PARAMETER (HALF=ONE/TWO)
      INTEGER MAXTYP
      PARAMETER (MAXTYP=21)
      LOGICAL CRANGE
      PARAMETER (CRANGE=.FALSE.)
      LOGICAL CREL
      PARAMETER (CREL=.FALSE.)
!     ..
!     .. Local Scalars ..
      LOGICAL badnn , tryrac
      INTEGER i , iinfo , il , imode , inde , indrwk , itemp , itype ,  &
     &        iu , j , jc , jr , jsize , jtype , lgn , liwedc , log2ui ,&
     &        lrwedc , lwedc , m , m2 , m3 , mtypes , n , nap , nblock ,&
     &        nerrs , nmats , nmax , nsplit , ntest , ntestt
      REAL abstol , aninv , anorm , cond , ovfl , rtovfl , rtunfl ,     &
     &     temp1 , temp2 , temp3 , temp4 , ulp , ulpinv , unfl , vl , vu
!     ..
!     .. Local Arrays ..
      INTEGER idumma(1) , ioldsd(4) , iseed2(4) , kmagn(MAXTYP) ,       &
     &        kmode(MAXTYP) , ktype(MAXTYP)
      REAL dumma(1)
!     ..
!     .. External Functions ..
      INTEGER ILAENV
      REAL SLAMCH , SLARND , SSXT1
      EXTERNAL ILAENV , SLAMCH , SLARND , SSXT1
!     ..
!     .. External Subroutines ..
      EXTERNAL CCOPY , CHET21 , CHETRD , CHPT21 , CHPTRD , CLACPY ,     &
     &         CLASET , CLATMR , CLATMS , CPTEQR , CSTEDC , CSTEMR ,    &
     &         CSTEIN , CSTEQR , CSTT21 , CSTT22 , CUNGTR , CUPGTR ,    &
     &         SCOPY , SLABAD , SLASUM , SSTEBZ , SSTECH , SSTERF ,     &
     &         XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CONJG , INT , LOG , MAX , MIN , REAL , SQRT
!     ..
!     .. Data statements ..
      DATA ktype/1 , 2 , 4 , 4 , 4 , 4 , 4 , 5 , 5 , 5 , 5 , 5 , 8 , 8 ,&
     &     8 , 9 , 9 , 9 , 9 , 9 , 10/
      DATA kmagn/1 , 1 , 1 , 1 , 1 , 2 , 3 , 1 , 1 , 1 , 2 , 3 , 1 , 2 ,&
     &     3 , 1 , 1 , 1 , 2 , 3 , 1/
      DATA kmode/0 , 0 , 4 , 3 , 1 , 4 , 4 , 4 , 3 , 1 , 4 , 4 , 0 , 0 ,&
     &     0 , 4 , 3 , 1 , 4 , 4 , 3/
!     ..
!     .. Executable Statements ..
!
!     Keep ftnchek happy
      idumma(1) = 1
!
!     Check for errors
!
      ntestt = 0
      Info = 0
!
!     Important constants
!
      badnn = .FALSE.
      tryrac = .TRUE.
      nmax = 1
      DO j = 1 , Nsizes
         nmax = MAX(nmax,Nn(j))
         IF ( Nn(j)<0 ) badnn = .TRUE.
      ENDDO
!
      nblock = ILAENV(1,'CHETRD','L',nmax,-1,-1,-1)
      nblock = MIN(nmax,MAX(1,nblock))
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
         Info = -23
      ELSEIF ( 2*MAX(2,nmax)**2>Lwork ) THEN
         Info = -29
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CCHKST',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( Nsizes==0 .OR. Ntypes==0 ) RETURN
!
!     More Important constants
!
      unfl = SLAMCH('Safe minimum')
      ovfl = ONE/unfl
      CALL SLABAD(unfl,ovfl)
      ulp = SLAMCH('Epsilon')*SLAMCH('Base')
      ulpinv = ONE/ulp
      log2ui = INT(LOG(ulpinv)/LOG(TWO))
      rtunfl = SQRT(unfl)
      rtovfl = SQRT(ovfl)
!
!     Loop over sizes, types
!
      DO i = 1 , 4
         iseed2(i) = Iseed(i)
      ENDDO
      nerrs = 0
      nmats = 0
!
      DO jsize = 1 , Nsizes
         n = Nn(jsize)
         IF ( n>0 ) THEN
            lgn = INT(LOG(REAL(n))/LOG(TWO))
            IF ( 2**lgn<n ) lgn = lgn + 1
            IF ( 2**lgn<n ) lgn = lgn + 1
            lwedc = 1 + 4*n + 2*n*lgn + 4*n**2
            lrwedc = 1 + 3*n + 2*n*lgn + 4*n**2
            liwedc = 6 + 6*n + 5*n*lgn
         ELSE
            lwedc = 8
            lrwedc = 7
            liwedc = 12
         ENDIF
         nap = (n*(n+1))/2
         aninv = ONE/REAL(MAX(1,n))
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
            DO j = 1 , 4
               ioldsd(j) = Iseed(j)
            ENDDO
!
!           Compute "A"
!
!           Control parameters:
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
!           =9                      positive definite
!           =10                     diagonally dominant tridiagonal
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
               IF ( jtype<=15 ) THEN
                  cond = ulpinv
               ELSE
                  cond = ulpinv*aninv/TEN
               ENDIF
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
                  DO jc = 1 , n
                     A(jc,jc) = anorm
                  ENDDO
!
               ELSEIF ( itype==4 ) THEN
!
!              Diagonal Matrix, [Eigen]values Specified
!
                  CALL CLATMS(n,n,'S',Iseed,'H',Rwork,imode,cond,anorm, &
     &                        0,0,'N',A,Lda,Work,iinfo)
!
!
               ELSEIF ( itype==5 ) THEN
!
!              Hermitian, eigenvalues specified
!
                  CALL CLATMS(n,n,'S',Iseed,'H',Rwork,imode,cond,anorm, &
     &                        n,n,'N',A,Lda,Work,iinfo)
!
               ELSEIF ( itype==7 ) THEN
!
!              Diagonal, random eigenvalues
!
                  CALL CLATMR(n,n,'S',Iseed,'H',Work,6,ONE,CONE,'T','N',&
     &                        Work(n+1),1,ONE,Work(2*n+1),1,ONE,'N',    &
     &                        idumma,0,0,ZERO,anorm,'NO',A,Lda,Iwork,   &
     &                        iinfo)
!
               ELSEIF ( itype==8 ) THEN
!
!              Hermitian, random eigenvalues
!
                  CALL CLATMR(n,n,'S',Iseed,'H',Work,6,ONE,CONE,'T','N',&
     &                        Work(n+1),1,ONE,Work(2*n+1),1,ONE,'N',    &
     &                        idumma,n,n,ZERO,anorm,'NO',A,Lda,Iwork,   &
     &                        iinfo)
!
               ELSEIF ( itype==9 ) THEN
!
!              Positive definite, eigenvalues specified.
!
                  CALL CLATMS(n,n,'S',Iseed,'P',Rwork,imode,cond,anorm, &
     &                        n,n,'N',A,Lda,Work,iinfo)
!
               ELSEIF ( itype==10 ) THEN
!
!              Positive definite tridiagonal, eigenvalues specified.
!
                  CALL CLATMS(n,n,'S',Iseed,'P',Rwork,imode,cond,anorm, &
     &                        1,1,'N',A,Lda,Work,iinfo)
                  DO i = 2 , n
                     temp1 = ABS(A(i-1,i))
                     temp2 = SQRT(ABS(A(i-1,i-1)*A(i,i)))
                     IF ( temp1>HALF*temp2 ) THEN
                        A(i-1,i) = A(i-1,i)*(HALF*temp2/(unfl+temp1))
                        A(i,i-1) = CONJG(A(i-1,i))
                     ENDIF
                  ENDDO
!
               ELSE
!
                  iinfo = 1
               ENDIF
!
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'Generator' , iinfo , n ,    &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  RETURN
               ENDIF
            ENDIF
!
!
!           Call CHETRD and CUNGTR to compute S and U from
!           upper triangle.
!
            CALL CLACPY('U',n,n,A,Lda,V,Ldu)
!
            ntest = 1
            CALL CHETRD('U',n,V,Ldu,Sd,Se,Tau,Work,Lwork,iinfo)
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CHETRD(U)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(1) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
            CALL CLACPY('U',n,n,V,Ldu,U,Ldu)
!
            ntest = 2
            CALL CUNGTR('U',n,U,Ldu,Tau,Work,Lwork,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CUNGTR(U)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(2) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Do tests 1 and 2
!
            CALL CHET21(2,'Upper',n,1,A,Lda,Sd,Se,U,Ldu,V,Ldu,Tau,Work, &
     &                  Rwork,Result(1))
            CALL CHET21(3,'Upper',n,1,A,Lda,Sd,Se,U,Ldu,V,Ldu,Tau,Work, &
     &                  Rwork,Result(2))
!
!           Call CHETRD and CUNGTR to compute S and U from
!           lower triangle, do tests.
!
            CALL CLACPY('L',n,n,A,Lda,V,Ldu)
!
            ntest = 3
            CALL CHETRD('L',n,V,Ldu,Sd,Se,Tau,Work,Lwork,iinfo)
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CHETRD(L)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(3) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
            CALL CLACPY('L',n,n,V,Ldu,U,Ldu)
!
            ntest = 4
            CALL CUNGTR('L',n,U,Ldu,Tau,Work,Lwork,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CUNGTR(L)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(4) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
            CALL CHET21(2,'Lower',n,1,A,Lda,Sd,Se,U,Ldu,V,Ldu,Tau,Work, &
     &                  Rwork,Result(3))
            CALL CHET21(3,'Lower',n,1,A,Lda,Sd,Se,U,Ldu,V,Ldu,Tau,Work, &
     &                  Rwork,Result(4))
!
!           Store the upper triangle of A in AP
!
            i = 0
            DO jc = 1 , n
               DO jr = 1 , jc
                  i = i + 1
                  Ap(i) = A(jr,jc)
               ENDDO
            ENDDO
!
!           Call CHPTRD and CUPGTR to compute S and U from AP
!
            CALL CCOPY(nap,Ap,1,Vp,1)
!
            ntest = 5
            CALL CHPTRD('U',n,Vp,Sd,Se,Tau,iinfo)
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CHPTRD(U)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(5) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
            ntest = 6
            CALL CUPGTR('U',n,Vp,Tau,U,Ldu,Work,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CUPGTR(U)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(6) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Do tests 5 and 6
!
            CALL CHPT21(2,'Upper',n,1,Ap,Sd,Se,U,Ldu,Vp,Tau,Work,Rwork, &
     &                  Result(5))
            CALL CHPT21(3,'Upper',n,1,Ap,Sd,Se,U,Ldu,Vp,Tau,Work,Rwork, &
     &                  Result(6))
!
!           Store the lower triangle of A in AP
!
            i = 0
            DO jc = 1 , n
               DO jr = jc , n
                  i = i + 1
                  Ap(i) = A(jr,jc)
               ENDDO
            ENDDO
!
!           Call CHPTRD and CUPGTR to compute S and U from AP
!
            CALL CCOPY(nap,Ap,1,Vp,1)
!
            ntest = 7
            CALL CHPTRD('L',n,Vp,Sd,Se,Tau,iinfo)
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CHPTRD(L)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(7) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
            ntest = 8
            CALL CUPGTR('L',n,Vp,Tau,U,Ldu,Work,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CUPGTR(L)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(8) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
            CALL CHPT21(2,'Lower',n,1,Ap,Sd,Se,U,Ldu,Vp,Tau,Work,Rwork, &
     &                  Result(7))
            CALL CHPT21(3,'Lower',n,1,Ap,Sd,Se,U,Ldu,Vp,Tau,Work,Rwork, &
     &                  Result(8))
!
!           Call CSTEQR to compute D1, D2, and Z, do tests.
!
!           Compute D1 and Z
!
            CALL SCOPY(n,Sd,1,D1,1)
            IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork,1)
            CALL CLASET('Full',n,n,CZERO,CONE,Z,Ldu)
!
            ntest = 9
            CALL CSTEQR('V',n,D1,Rwork,Z,Ldu,Rwork(n+1),iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CSTEQR(V)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(9) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Compute D2
!
            CALL SCOPY(n,Sd,1,D2,1)
            IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork,1)
!
            ntest = 11
            CALL CSTEQR('N',n,D2,Rwork,Work,Ldu,Rwork(n+1),iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CSTEQR(N)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(11) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Compute D3 (using PWK method)
!
            CALL SCOPY(n,Sd,1,D3,1)
            IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork,1)
!
            ntest = 12
            CALL SSTERF(n,D3,Rwork,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'SSTERF' , iinfo , n , jtype ,  &
     &                                  ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(12) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Do Tests 9 and 10
!
            CALL CSTT21(n,0,Sd,Se,D1,dumma,Z,Ldu,Work,Rwork,Result(9))
!
!           Do Tests 11 and 12
!
            temp1 = ZERO
            temp2 = ZERO
            temp3 = ZERO
            temp4 = ZERO
!
            DO j = 1 , n
               temp1 = MAX(temp1,ABS(D1(j)),ABS(D2(j)))
               temp2 = MAX(temp2,ABS(D1(j)-D2(j)))
               temp3 = MAX(temp3,ABS(D1(j)),ABS(D3(j)))
               temp4 = MAX(temp4,ABS(D1(j)-D3(j)))
            ENDDO
!
            Result(11) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
            Result(12) = temp4/MAX(unfl,ulp*MAX(temp3,temp4))
!
!           Do Test 13 -- Sturm Sequence Test of Eigenvalues
!                         Go up by factors of two until it succeeds
!
            ntest = 13
            temp1 = Thresh*(HALF-ulp)
!
            DO j = 0 , log2ui
               CALL SSTECH(n,Sd,Se,D1,temp1,Rwork,iinfo)
               IF ( iinfo==0 ) EXIT
               temp1 = temp1*TWO
            ENDDO
!
            Result(13) = temp1
!
!           For positive definite matrices ( JTYPE.GT.15 ) call CPTEQR
!           and do tests 14, 15, and 16 .
!
            IF ( jtype>15 ) THEN
!
!              Compute D4 and Z4
!
               CALL SCOPY(n,Sd,1,D4,1)
               IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork,1)
               CALL CLASET('Full',n,n,CZERO,CONE,Z,Ldu)
!
               ntest = 14
               CALL CPTEQR('V',n,D4,Rwork,Z,Ldu,Rwork(n+1),iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'CPTEQR(V)' , iinfo , n ,    &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  IF ( iinfo<0 ) THEN
                     RETURN
                  ELSE
                     Result(14) = ulpinv
                     GOTO 20
                  ENDIF
               ENDIF
!
!              Do Tests 14 and 15
!
               CALL CSTT21(n,0,Sd,Se,D4,dumma,Z,Ldu,Work,Rwork,         &
     &                     Result(14))
!
!              Compute D5
!
               CALL SCOPY(n,Sd,1,D5,1)
               IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork,1)
!
               ntest = 16
               CALL CPTEQR('N',n,D5,Rwork,Z,Ldu,Rwork(n+1),iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'CPTEQR(N)' , iinfo , n ,    &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  IF ( iinfo<0 ) THEN
                     RETURN
                  ELSE
                     Result(16) = ulpinv
                     GOTO 20
                  ENDIF
               ENDIF
!
!              Do Test 16
!
               temp1 = ZERO
               temp2 = ZERO
               DO j = 1 , n
                  temp1 = MAX(temp1,ABS(D4(j)),ABS(D5(j)))
                  temp2 = MAX(temp2,ABS(D4(j)-D5(j)))
               ENDDO
!
               Result(16) = temp2/MAX(unfl,HUN*ulp*MAX(temp1,temp2))
            ELSE
               Result(14) = ZERO
               Result(15) = ZERO
               Result(16) = ZERO
            ENDIF
!
!           Call SSTEBZ with different options and do tests 17-18.
!
!              If S is positive definite and diagonally dominant,
!              ask for all eigenvalues with high relative accuracy.
!
            vl = ZERO
            vu = ZERO
            il = 0
            iu = 0
            IF ( jtype==21 ) THEN
               ntest = 17
               abstol = unfl + unfl
               CALL SSTEBZ('A','E',n,vl,vu,il,iu,abstol,Sd,Se,m,nsplit, &
     &                     Wr,Iwork(1),Iwork(n+1),Rwork,Iwork(2*n+1),   &
     &                     iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'SSTEBZ(A,rel)' , iinfo , n ,&
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  IF ( iinfo<0 ) THEN
                     RETURN
                  ELSE
                     Result(17) = ulpinv
                     GOTO 20
                  ENDIF
               ENDIF
!
!              Do test 17
!
               temp2 = TWO*(TWO*n-ONE)*ulp*(ONE+EIGHT*HALF**2)          &
     &                 /(ONE-HALF)**4
!
               temp1 = ZERO
               DO j = 1 , n
                  temp1 = MAX(temp1,ABS(D4(j)-Wr(n-j+1))                &
     &                    /(abstol+ABS(D4(j))))
               ENDDO
!
               Result(17) = temp1/temp2
            ELSE
               Result(17) = ZERO
            ENDIF
!
!           Now ask for all eigenvalues with high absolute accuracy.
!
            ntest = 18
            abstol = unfl + unfl
            CALL SSTEBZ('A','E',n,vl,vu,il,iu,abstol,Sd,Se,m,nsplit,Wa1,&
     &                  Iwork(1),Iwork(n+1),Rwork,Iwork(2*n+1),iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'SSTEBZ(A)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(18) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Do test 18
!
            temp1 = ZERO
            temp2 = ZERO
            DO j = 1 , n
               temp1 = MAX(temp1,ABS(D3(j)),ABS(Wa1(j)))
               temp2 = MAX(temp2,ABS(D3(j)-Wa1(j)))
            ENDDO
!
            Result(18) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!           Choose random values for IL and IU, and ask for the
!           IL-th through IU-th eigenvalues.
!
            ntest = 19
            IF ( n<=1 ) THEN
               il = 1
               iu = n
            ELSE
               il = 1 + (n-1)*INT(SLARND(1,iseed2))
               iu = 1 + (n-1)*INT(SLARND(1,iseed2))
               IF ( iu<il ) THEN
                  itemp = iu
                  iu = il
                  il = itemp
               ENDIF
            ENDIF
!
            CALL SSTEBZ('I','E',n,vl,vu,il,iu,abstol,Sd,Se,m2,nsplit,   &
     &                  Wa2,Iwork(1),Iwork(n+1),Rwork,Iwork(2*n+1),     &
     &                  iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'SSTEBZ(I)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(19) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Determine the values VL and VU of the IL-th and IU-th
!           eigenvalues and ask for all eigenvalues in this range.
!
            IF ( n>0 ) THEN
               IF ( il/=1 ) THEN
                  vl = Wa1(il) - MAX(HALF*(Wa1(il)-Wa1(il-1)),ulp*anorm,&
     &                 TWO*rtunfl)
               ELSE
                  vl = Wa1(1) - MAX(HALF*(Wa1(n)-Wa1(1)),ulp*anorm,     &
     &                 TWO*rtunfl)
               ENDIF
               IF ( iu/=n ) THEN
                  vu = Wa1(iu) + MAX(HALF*(Wa1(iu+1)-Wa1(iu)),ulp*anorm,&
     &                 TWO*rtunfl)
               ELSE
                  vu = Wa1(n) + MAX(HALF*(Wa1(n)-Wa1(1)),ulp*anorm,     &
     &                 TWO*rtunfl)
               ENDIF
            ELSE
               vl = ZERO
               vu = ONE
            ENDIF
!
            CALL SSTEBZ('V','E',n,vl,vu,il,iu,abstol,Sd,Se,m3,nsplit,   &
     &                  Wa3,Iwork(1),Iwork(n+1),Rwork,Iwork(2*n+1),     &
     &                  iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'SSTEBZ(V)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(19) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
            IF ( m3==0 .AND. n/=0 ) THEN
               Result(19) = ulpinv
               GOTO 20
            ENDIF
!
!           Do test 19
!
            temp1 = SSXT1(1,Wa2,m2,Wa3,m3,abstol,ulp,unfl)
            temp2 = SSXT1(1,Wa3,m3,Wa2,m2,abstol,ulp,unfl)
            IF ( n>0 ) THEN
               temp3 = MAX(ABS(Wa1(n)),ABS(Wa1(1)))
            ELSE
               temp3 = ZERO
            ENDIF
!
            Result(19) = (temp1+temp2)/MAX(unfl,temp3*ulp)
!
!           Call CSTEIN to compute eigenvectors corresponding to
!           eigenvalues in WA1.  (First call SSTEBZ again, to make sure
!           it returns these eigenvalues in the correct order.)
!
            ntest = 21
            CALL SSTEBZ('A','B',n,vl,vu,il,iu,abstol,Sd,Se,m,nsplit,Wa1,&
     &                  Iwork(1),Iwork(n+1),Rwork,Iwork(2*n+1),iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'SSTEBZ(A,B)' , iinfo , n ,     &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(20) = ulpinv
                  Result(21) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
            CALL CSTEIN(n,Sd,Se,m,Wa1,Iwork(1),Iwork(n+1),Z,Ldu,Rwork,  &
     &                  Iwork(2*n+1),Iwork(3*n+1),iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CSTEIN' , iinfo , n , jtype ,  &
     &                                  ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(20) = ulpinv
                  Result(21) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Do tests 20 and 21
!
            CALL CSTT21(n,0,Sd,Se,Wa1,dumma,Z,Ldu,Work,Rwork,Result(20))
!
!           Call CSTEDC(I) to compute D1 and Z, do tests.
!
!           Compute D1 and Z
!
            inde = 1
            indrwk = inde + n
            CALL SCOPY(n,Sd,1,D1,1)
            IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork(inde),1)
            CALL CLASET('Full',n,n,CZERO,CONE,Z,Ldu)
!
            ntest = 22
            CALL CSTEDC('I',n,D1,Rwork(inde),Z,Ldu,Work,lwedc,          &
     &                  Rwork(indrwk),lrwedc,Iwork,liwedc,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CSTEDC(I)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(22) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Do Tests 22 and 23
!
            CALL CSTT21(n,0,Sd,Se,D1,dumma,Z,Ldu,Work,Rwork,Result(22))
!
!           Call CSTEDC(V) to compute D1 and Z, do tests.
!
!           Compute D1 and Z
!
            CALL SCOPY(n,Sd,1,D1,1)
            IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork(inde),1)
            CALL CLASET('Full',n,n,CZERO,CONE,Z,Ldu)
!
            ntest = 24
            CALL CSTEDC('V',n,D1,Rwork(inde),Z,Ldu,Work,lwedc,          &
     &                  Rwork(indrwk),lrwedc,Iwork,liwedc,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CSTEDC(V)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(24) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Do Tests 24 and 25
!
            CALL CSTT21(n,0,Sd,Se,D1,dumma,Z,Ldu,Work,Rwork,Result(24))
!
!           Call CSTEDC(N) to compute D2, do tests.
!
!           Compute D2
!
            CALL SCOPY(n,Sd,1,D2,1)
            IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork(inde),1)
            CALL CLASET('Full',n,n,CZERO,CONE,Z,Ldu)
!
            ntest = 26
            CALL CSTEDC('N',n,D2,Rwork(inde),Z,Ldu,Work,lwedc,          &
     &                  Rwork(indrwk),lrwedc,Iwork,liwedc,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'CSTEDC(N)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  Result(26) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Do Test 26
!
            temp1 = ZERO
            temp2 = ZERO
!
            DO j = 1 , n
               temp1 = MAX(temp1,ABS(D1(j)),ABS(D2(j)))
               temp2 = MAX(temp2,ABS(D1(j)-D2(j)))
            ENDDO
!
            Result(26) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!           Only test CSTEMR if IEEE compliant
!
            IF ( ILAENV(10,'CSTEMR','VA',1,0,0,0)==1 .AND.              &
     &           ILAENV(11,'CSTEMR','VA',1,0,0,0)==1 ) THEN
!
!           Call CSTEMR, do test 27 (relative eigenvalue accuracy)
!
!              If S is positive definite and diagonally dominant,
!              ask for all eigenvalues with high relative accuracy.
!
               vl = ZERO
               vu = ZERO
               il = 0
               iu = 0
               IF ( jtype==21 .AND. CREL ) THEN
                  ntest = 27
                  abstol = unfl + unfl
                  CALL CSTEMR('V','A',n,Sd,Se,vl,vu,il,iu,m,Wr,Z,Ldu,n, &
     &                        Iwork(1),tryrac,Rwork,Lrwork,Iwork(2*n+1),&
     &                        Lwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CSTEMR(V,A,rel)' ,       &
     &                      iinfo , n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(27) = ulpinv
                        GOTO 20
                     ENDIF
                  ENDIF
!
!              Do test 27
!
                  temp2 = TWO*(TWO*n-ONE)*ulp*(ONE+EIGHT*HALF**2)       &
     &                    /(ONE-HALF)**4
!
                  temp1 = ZERO
                  DO j = 1 , n
                     temp1 = MAX(temp1,ABS(D4(j)-Wr(n-j+1))             &
     &                       /(abstol+ABS(D4(j))))
                  ENDDO
!
                  Result(27) = temp1/temp2
!
                  il = 1 + (n-1)*INT(SLARND(1,iseed2))
                  iu = 1 + (n-1)*INT(SLARND(1,iseed2))
                  IF ( iu<il ) THEN
                     itemp = iu
                     iu = il
                     il = itemp
                  ENDIF
!
                  IF ( CRANGE ) THEN
                     ntest = 28
                     abstol = unfl + unfl
                     CALL CSTEMR('V','I',n,Sd,Se,vl,vu,il,iu,m,Wr,Z,Ldu,&
     &                           n,Iwork(1),tryrac,Rwork,Lrwork,        &
     &                           Iwork(2*n+1),Lwork-2*n,iinfo)
!
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'CSTEMR(V,I,rel)' ,    &
     &                         iinfo , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(28) = ulpinv
                           GOTO 20
                        ENDIF
                     ENDIF
!
!
!                 Do test 28
!
                     temp2 = TWO*(TWO*n-ONE)*ulp*(ONE+EIGHT*HALF**2)    &
     &                       /(ONE-HALF)**4
!
                     temp1 = ZERO
                     DO j = il , iu
                        temp1 = MAX(temp1,ABS(Wr(j-il+1)-D4(n-j+1))     &
     &                          /(abstol+ABS(Wr(j-il+1))))
                     ENDDO
!
                     Result(28) = temp1/temp2
                  ELSE
                     Result(28) = ZERO
                  ENDIF
               ELSE
                  Result(27) = ZERO
                  Result(28) = ZERO
               ENDIF
!
!           Call CSTEMR(V,I) to compute D1 and Z, do tests.
!
!           Compute D1 and Z
!
               CALL SCOPY(n,Sd,1,D5,1)
               IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork,1)
               CALL CLASET('Full',n,n,CZERO,CONE,Z,Ldu)
!
               IF ( CRANGE ) THEN
                  ntest = 29
                  il = 1 + (n-1)*INT(SLARND(1,iseed2))
                  iu = 1 + (n-1)*INT(SLARND(1,iseed2))
                  IF ( iu<il ) THEN
                     itemp = iu
                     iu = il
                     il = itemp
                  ENDIF
                  CALL CSTEMR('V','I',n,D5,Rwork,vl,vu,il,iu,m,D1,Z,Ldu,&
     &                        n,Iwork(1),tryrac,Rwork(n+1),Lrwork-n,    &
     &                        Iwork(2*n+1),Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CSTEMR(V,I)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(29) = ulpinv
                        GOTO 20
                     ENDIF
                  ENDIF
!
!           Do Tests 29 and 30
!
!
!           Call CSTEMR to compute D2, do tests.
!
!           Compute D2
!
                  CALL SCOPY(n,Sd,1,D5,1)
                  IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork,1)
!
                  ntest = 31
                  CALL CSTEMR('N','I',n,D5,Rwork,vl,vu,il,iu,m,D2,Z,Ldu,&
     &                        n,Iwork(1),tryrac,Rwork(n+1),Lrwork-n,    &
     &                        Iwork(2*n+1),Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CSTEMR(N,I)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(31) = ulpinv
                        GOTO 20
                     ENDIF
                  ENDIF
!
!           Do Test 31
!
                  temp1 = ZERO
                  temp2 = ZERO
!
                  DO j = 1 , iu - il + 1
                     temp1 = MAX(temp1,ABS(D1(j)),ABS(D2(j)))
                     temp2 = MAX(temp2,ABS(D1(j)-D2(j)))
                  ENDDO
!
                  Result(31) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
!
!
!           Call CSTEMR(V,V) to compute D1 and Z, do tests.
!
!           Compute D1 and Z
!
                  CALL SCOPY(n,Sd,1,D5,1)
                  IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork,1)
                  CALL CLASET('Full',n,n,CZERO,CONE,Z,Ldu)
!
                  ntest = 32
!
                  IF ( n>0 ) THEN
                     IF ( il/=1 ) THEN
                        vl = D2(il)                                     &
     &                       - MAX(HALF*(D2(il)-D2(il-1)),ulp*anorm,    &
     &                       TWO*rtunfl)
                     ELSE
                        vl = D2(1) - MAX(HALF*(D2(n)-D2(1)),ulp*anorm,  &
     &                       TWO*rtunfl)
                     ENDIF
                     IF ( iu/=n ) THEN
                        vu = D2(iu)                                     &
     &                       + MAX(HALF*(D2(iu+1)-D2(iu)),ulp*anorm,    &
     &                       TWO*rtunfl)
                     ELSE
                        vu = D2(n) + MAX(HALF*(D2(n)-D2(1)),ulp*anorm,  &
     &                       TWO*rtunfl)
                     ENDIF
                  ELSE
                     vl = ZERO
                     vu = ONE
                  ENDIF
!
                  CALL CSTEMR('V','V',n,D5,Rwork,vl,vu,il,iu,m,D1,Z,Ldu,&
     &                        n,Iwork(1),tryrac,Rwork(n+1),Lrwork-n,    &
     &                        Iwork(2*n+1),Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CSTEMR(V,V)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(32) = ulpinv
                        GOTO 20
                     ENDIF
                  ENDIF
!
!           Do Tests 32 and 33
!
                  CALL CSTT22(n,m,0,Sd,Se,D1,dumma,Z,Ldu,Work,m,Rwork,  &
     &                        Result(32))
!
!           Call CSTEMR to compute D2, do tests.
!
!           Compute D2
!
                  CALL SCOPY(n,Sd,1,D5,1)
                  IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork,1)
!
                  ntest = 34
                  CALL CSTEMR('N','V',n,D5,Rwork,vl,vu,il,iu,m,D2,Z,Ldu,&
     &                        n,Iwork(1),tryrac,Rwork(n+1),Lrwork-n,    &
     &                        Iwork(2*n+1),Liwork-2*n,iinfo)
                  IF ( iinfo/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'CSTEMR(N,V)' , iinfo ,   &
     &                      n , jtype , ioldsd
                     Info = ABS(iinfo)
                     IF ( iinfo<0 ) THEN
                        RETURN
                     ELSE
                        Result(34) = ulpinv
                        GOTO 20
                     ENDIF
                  ENDIF
!
!           Do Test 34
!
                  temp1 = ZERO
                  temp2 = ZERO
!
                  DO j = 1 , iu - il + 1
                     temp1 = MAX(temp1,ABS(D1(j)),ABS(D2(j)))
                     temp2 = MAX(temp2,ABS(D1(j)-D2(j)))
                  ENDDO
!
                  Result(34) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
               ELSE
                  Result(29) = ZERO
                  Result(30) = ZERO
                  Result(31) = ZERO
                  Result(32) = ZERO
                  Result(33) = ZERO
                  Result(34) = ZERO
               ENDIF
!
!
!           Call CSTEMR(V,A) to compute D1 and Z, do tests.
!
!           Compute D1 and Z
!
               CALL SCOPY(n,Sd,1,D5,1)
               IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork,1)
!
               ntest = 35
!
               CALL CSTEMR('V','A',n,D5,Rwork,vl,vu,il,iu,m,D1,Z,Ldu,n, &
     &                     Iwork(1),tryrac,Rwork(n+1),Lrwork-n,         &
     &                     Iwork(2*n+1),Liwork-2*n,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'CSTEMR(V,A)' , iinfo , n ,  &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  IF ( iinfo<0 ) THEN
                     RETURN
                  ELSE
                     Result(35) = ulpinv
                     GOTO 20
                  ENDIF
               ENDIF
!
!           Do Tests 35 and 36
!
               CALL CSTT22(n,m,0,Sd,Se,D1,dumma,Z,Ldu,Work,m,Rwork,     &
     &                     Result(35))
!
!           Call CSTEMR to compute D2, do tests.
!
!           Compute D2
!
               CALL SCOPY(n,Sd,1,D5,1)
               IF ( n>0 ) CALL SCOPY(n-1,Se,1,Rwork,1)
!
               ntest = 37
               CALL CSTEMR('N','A',n,D5,Rwork,vl,vu,il,iu,m,D2,Z,Ldu,n, &
     &                     Iwork(1),tryrac,Rwork(n+1),Lrwork-n,         &
     &                     Iwork(2*n+1),Liwork-2*n,iinfo)
               IF ( iinfo/=0 ) THEN
                  WRITE (Nounit,FMT=99001) 'CSTEMR(N,A)' , iinfo , n ,  &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  IF ( iinfo<0 ) THEN
                     RETURN
                  ELSE
                     Result(37) = ulpinv
                     GOTO 20
                  ENDIF
               ENDIF
!
!           Do Test 34
!
               temp1 = ZERO
               temp2 = ZERO
!
               DO j = 1 , n
                  temp1 = MAX(temp1,ABS(D1(j)),ABS(D2(j)))
                  temp2 = MAX(temp2,ABS(D1(j)-D2(j)))
               ENDDO
!
               Result(37) = temp2/MAX(unfl,ulp*MAX(temp1,temp2))
            ENDIF
 20         ntestt = ntestt + ntest
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
!
!           Print out tests which fail.
!
            DO jr = 1 , ntest
               IF ( Result(jr)>=Thresh ) THEN
!
!                 If this is the first test to fail,
!                 print a header to the data file.
!
                  IF ( nerrs==0 ) THEN
                     WRITE (Nounit,FMT=99002) 'CST'
                     WRITE (Nounit,FMT=99003)
                     WRITE (Nounit,FMT=99004)
                     WRITE (Nounit,FMT=99005) 'Hermitian'
                     WRITE (Nounit,FMT=99006)
!
!                    Tests performed
!
                     WRITE (Nounit,FMT=99009)
                  ENDIF
                  nerrs = nerrs + 1
                  IF ( Result(jr)<10000.0E0 ) THEN
                     WRITE (Nounit,FMT=99007) n , jtype , ioldsd , jr , &
     &                      Result(jr)
                  ELSE
                     WRITE (Nounit,FMT=99008) n , jtype , ioldsd , jr , &
     &                      Result(jr)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
!     Summary
!
      CALL SLASUM('CST',Nounit,nerrs,ntestt)
      RETURN
!
99001 FORMAT (' CCHKST: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
99002 FORMAT (/1X,A3,' -- Complex Hermitian eigenvalue problem')
99003 FORMAT (' Matrix types (see CCHKST for details): ')
!
99004 FORMAT (/' Special Matrices:',                                    &
     &        /'  1=Zero matrix.                        ',              &
     &        '  5=Diagonal: clustered entries.',                       &
     &        /'  2=Identity matrix.                    ',              &
     &        '  6=Diagonal: large, evenly spaced.',                    &
     &        /'  3=Diagonal: evenly spaced entries.    ',              &
     &        '  7=Diagonal: small, evenly spaced.',                    &
     &        /'  4=Diagonal: geometr. spaced entries.')
99005 FORMAT (' Dense ',A,' Matrices:',                                 &
     &        /'  8=Evenly spaced eigenvals.            ',              &
     &        ' 12=Small, evenly spaced eigenvals.',                    &
     &        /'  9=Geometrically spaced eigenvals.     ',              &
     &        ' 13=Matrix with random O(1) entries.',                   &
     &        /' 10=Clustered eigenvalues.              ',              &
     &        ' 14=Matrix with large random entries.',                  &
     &        /' 11=Large, evenly spaced eigenvals.     ',              &
     &        ' 15=Matrix with small random entries.')
99006 FORMAT (' 16=Positive definite, evenly spaced eigenvalues',       &
     &        /' 17=Positive definite, geometrically spaced eigenvlaues'&
     &        ,/' 18=Positive definite, clustered eigenvalues',         &
     &        /' 19=Positive definite, small evenly spaced eigenvalues',&
     &        /' 20=Positive definite, large evenly spaced eigenvalues',&
     &        /' 21=Diagonally dominant tridiagonal, geometrically',    &
     &        ' spaced eigenvalues')
!
99007 FORMAT (' Matrix order=',I5,', type=',I2,', seed=',4(I4,','),     &
     &        ' result ',I3,' is',0P,F8.2)
99008 FORMAT (' Matrix order=',I5,', type=',I2,', seed=',4(I4,','),     &
     &        ' result ',I3,' is',1P,E10.3)
!
99009 FORMAT (/'Test performed:  see CCHKST for details.',/)
!     End of CCHKST
!
      END SUBROUTINE CCHKST