!*==sdrvvx.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b sdrvvx
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sdrvvx( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          NIUNIT, NOUNIT, A, LDA, H, WR, WI, WR1, WI1,
!                          VL, LDVL, VR, LDVR, LRE, LDLRE, RCONDV, RCNDV1,
!                          RCDVIN, RCONDE, RCNDE1, RCDEIN, SCALE, SCALE1,
!                          RESULT, WORK, NWORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDLRE, LDVL, LDVR, NIUNIT, NOUNIT,
!      $                   NSIZES, NTYPES, NWORK
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
!       REAL               A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ),
!      $                   RCDEIN( * ), RCDVIN( * ), RCNDE1( * ),
!      $                   RCNDV1( * ), RCONDE( * ), RCONDV( * ),
!      $                   RESULT( 11 ), SCALE( * ), SCALE1( * ),
!      $                   VL( LDVL, * ), VR( LDVR, * ), WI( * ),
!      $                   WI1( * ), WORK( * ), WR( * ), WR1( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    sdrvvx  checks the nonsymmetric eigenvalue problem expert driver
!>    SGEEVX.
!>
!>    sdrvvx uses both test matrices generated randomly depending on
!>    data supplied in the calling sequence, as well as on data
!>    read from an input file and including precomputed condition
!>    numbers to which it compares the ones it computes.
!>
!>    When sdrvvx is called, a number of matrix "sizes" ("n's") and a
!>    number of matrix "types" are specified in the calling sequence.
!>    For each size ("n") and each type of matrix, one matrix will be
!>    generated and used to test the nonsymmetric eigenroutines.  For
!>    each matrix, 9 tests will be performed:
!>
!>    (1)     | A * VR - VR * W | / ( n |A| ulp )
!>
!>      Here VR is the matrix of unit right eigenvectors.
!>      W is a block diagonal matrix, with a 1x1 block for each
!>      real eigenvalue and a 2x2 block for each complex conjugate
!>      pair.  If eigenvalues j and j+1 are a complex conjugate pair,
!>      so WR(j) = WR(j+1) = wr and WI(j) = - WI(j+1) = wi, then the
!>      2 x 2 block corresponding to the pair will be:
!>
!>              (  wr  wi  )
!>              ( -wi  wr  )
!>
!>      Such a block multiplying an n x 2 matrix  ( ur ui ) on the
!>      right will be the same as multiplying  ur + i*ui  by  wr + i*wi.
!>
!>    (2)     | A**H * VL - VL * W**H | / ( n |A| ulp )
!>
!>      Here VL is the matrix of unit left eigenvectors, A**H is the
!>      conjugate transpose of A, and W is as above.
!>
!>    (3)     | |VR(i)| - 1 | / ulp and largest component real
!>
!>      VR(i) denotes the i-th column of VR.
!>
!>    (4)     | |VL(i)| - 1 | / ulp and largest component real
!>
!>      VL(i) denotes the i-th column of VL.
!>
!>    (5)     W(full) = W(partial)
!>
!>      W(full) denotes the eigenvalues computed when VR, VL, RCONDV
!>      and RCONDE are also computed, and W(partial) denotes the
!>      eigenvalues computed when only some of VR, VL, RCONDV, and
!>      RCONDE are computed.
!>
!>    (6)     VR(full) = VR(partial)
!>
!>      VR(full) denotes the right eigenvectors computed when VL, RCONDV
!>      and RCONDE are computed, and VR(partial) denotes the result
!>      when only some of VL and RCONDV are computed.
!>
!>    (7)     VL(full) = VL(partial)
!>
!>      VL(full) denotes the left eigenvectors computed when VR, RCONDV
!>      and RCONDE are computed, and VL(partial) denotes the result
!>      when only some of VR and RCONDV are computed.
!>
!>    (8)     0 if SCALE, ILO, IHI, ABNRM (full) =
!>                 SCALE, ILO, IHI, ABNRM (partial)
!>            1/ulp otherwise
!>
!>      SCALE, ILO, IHI and ABNRM describe how the matrix is balanced.
!>      (full) is when VR, VL, RCONDE and RCONDV are also computed, and
!>      (partial) is when some are not computed.
!>
!>    (9)     RCONDV(full) = RCONDV(partial)
!>
!>      RCONDV(full) denotes the reciprocal condition numbers of the
!>      right eigenvectors computed when VR, VL and RCONDE are also
!>      computed. RCONDV(partial) denotes the reciprocal condition
!>      numbers when only some of VR, VL and RCONDE are computed.
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
!>
!>    In addition, an input file will be read from logical unit number
!>    NIUNIT. The file contains matrices along with precomputed
!>    eigenvalues and reciprocal condition numbers for the eigenvalues
!>    and right eigenvectors. For these matrices, in addition to tests
!>    (1) to (9) we will compute the following two tests:
!>
!>   (10)  |RCONDV - RCDVIN| / cond(RCONDV)
!>
!>      RCONDV is the reciprocal right eigenvector condition number
!>      computed by SGEEVX and RCDVIN (the precomputed true value)
!>      is supplied as input. cond(RCONDV) is the condition number of
!>      RCONDV, and takes errors in computing RCONDV into account, so
!>      that the resulting quantity should be O(ULP). cond(RCONDV) is
!>      essentially given by norm(A)/RCONDE.
!>
!>   (11)  |RCONDE - RCDEIN| / cond(RCONDE)
!>
!>      RCONDE is the reciprocal eigenvalue condition number
!>      computed by SGEEVX and RCDEIN (the precomputed true value)
!>      is supplied as input.  cond(RCONDE) is the condition number
!>      of RCONDE, and takes errors in computing RCONDE into account,
!>      so that the resulting quantity should be O(ULP). cond(RCONDE)
!>      is essentially given by norm(A)/RCONDV.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of sizes of matrices to use.  NSIZES must be at
!>          least zero. If it is zero, no randomly generated matrices
!>          are tested, but any test matrices read from NIUNIT will be
!>          tested.
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
!>          The number of elements in DOTYPE. NTYPES must be at least
!>          zero. If it is zero, no randomly generated test matrices
!>          are tested, but and test matrices read from NIUNIT will be
!>          tested. If it is MAXTYP+1 and NSIZES is 1, then an
!>          additional type, MAXTYP+1 is defined, which is to use
!>          whatever matrix is in A.  This is only useful if
!>          DOTYPE(1:MAXTYP) is .FALSE. and DOTYPE(MAXTYP+1) is .TRUE. .
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
!>          next call to sdrvvx to continue the same random number
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
!> \param[in] NIUNIT
!> \verbatim
!>          NIUNIT is INTEGER
!>          The FORTRAN unit number for reading in the data file of
!>          problems to solve.
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
!>          A is REAL array, dimension
!>                      (LDA, max(NN,12))
!>          Used to hold the matrix whose eigenvalues are to be
!>          computed.  On exit, A contains the last matrix actually used.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A and H.
!>          LDA >= max(NN,12), since 12 is the dimension of the largest
!>          matrix in the precomputed input file.
!> \endverbatim
!>
!> \param[out] H
!> \verbatim
!>          H is REAL array, dimension
!>                      (LDA, max(NN,12))
!>          Another copy of the test matrix A, modified by SGEEVX.
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
!>          The real and imaginary parts of the eigenvalues of A.
!>          On exit, WR + WI*i are the eigenvalues of the matrix in A.
!> \endverbatim
!>
!> \param[out] WR1
!> \verbatim
!>          WR1 is REAL array, dimension (max(NN,12))
!> \endverbatim
!>
!> \param[out] WI1
!> \verbatim
!>          WI1 is REAL array, dimension (max(NN,12))
!>
!>          Like WR, WI, these arrays contain the eigenvalues of A,
!>          but those computed when SGEEVX only computes a partial
!>          eigendecomposition, i.e. not the eigenvalues and left
!>          and right eigenvectors.
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is REAL array, dimension
!>                      (LDVL, max(NN,12))
!>          VL holds the computed left eigenvectors.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          Leading dimension of VL. Must be at least max(1,max(NN,12)).
!> \endverbatim
!>
!> \param[out] VR
!> \verbatim
!>          VR is REAL array, dimension
!>                      (LDVR, max(NN,12))
!>          VR holds the computed right eigenvectors.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          Leading dimension of VR. Must be at least max(1,max(NN,12)).
!> \endverbatim
!>
!> \param[out] LRE
!> \verbatim
!>          LRE is REAL array, dimension
!>                      (LDLRE, max(NN,12))
!>          LRE holds the computed right or left eigenvectors.
!> \endverbatim
!>
!> \param[in] LDLRE
!> \verbatim
!>          LDLRE is INTEGER
!>          Leading dimension of LRE. Must be at least max(1,max(NN,12))
!> \endverbatim
!>
!> \param[out] RCONDV
!> \verbatim
!>          RCONDV is REAL array, dimension (N)
!>          RCONDV holds the computed reciprocal condition numbers
!>          for eigenvectors.
!> \endverbatim
!>
!> \param[out] RCNDV1
!> \verbatim
!>          RCNDV1 is REAL array, dimension (N)
!>          RCNDV1 holds more computed reciprocal condition numbers
!>          for eigenvectors.
!> \endverbatim
!>
!> \param[out] RCDVIN
!> \verbatim
!>          RCDVIN is REAL array, dimension (N)
!>          When COMP = .TRUE. RCDVIN holds the precomputed reciprocal
!>          condition numbers for eigenvectors to be compared with
!>          RCONDV.
!> \endverbatim
!>
!> \param[out] RCONDE
!> \verbatim
!>          RCONDE is REAL array, dimension (N)
!>          RCONDE holds the computed reciprocal condition numbers
!>          for eigenvalues.
!> \endverbatim
!>
!> \param[out] RCNDE1
!> \verbatim
!>          RCNDE1 is REAL array, dimension (N)
!>          RCNDE1 holds more computed reciprocal condition numbers
!>          for eigenvalues.
!> \endverbatim
!>
!> \param[out] RCDEIN
!> \verbatim
!>          RCDEIN is REAL array, dimension (N)
!>          When COMP = .TRUE. RCDEIN holds the precomputed reciprocal
!>          condition numbers for eigenvalues to be compared with
!>          RCONDE.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is REAL array, dimension (N)
!>          Holds information describing balancing of matrix.
!> \endverbatim
!>
!> \param[out] SCALE1
!> \verbatim
!>          SCALE1 is REAL array, dimension (N)
!>          Holds information describing balancing of matrix.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (11)
!>          The values computed by the seven tests described above.
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
!>          max(6*12+2*12**2,6*NN(j)+2*NN(j)**2) =
!>          max(    360     ,6*NN(j)+2*NN(j)**2)    for all j.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (2*max(NN,12))
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0,  then successful exit.
!>          If <0, then input parameter -INFO is incorrect.
!>          If >0, slatmr, slatms, slatme or sget23 returned an error
!>                 code, and INFO is its absolute value.
!>
!>-----------------------------------------------------------------------
!>
!>     Some Local Variables and Parameters:
!>     ---- ----- --------- --- ----------
!>
!>     ZERO, ONE       Real 0 and 1.
!>     MAXTYP          The number of types defined.
!>     NMAX            Largest value in NN or 12.
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
      SUBROUTINE sdrvvx(Nsizes,Nn,Ntypes,Dotype,Iseed,Thresh,Niunit,    &
     &                  Nounit,A,Lda,H,Wr,Wi,Wr1,Wi1,Vl,Ldvl,Vr,Ldvr,   &
     &                  Lre,Ldlre,Rcondv,Rcndv1,Rcdvin,Rconde,Rcnde1,   &
     &                  Rcdein,Scale,Scale1,Result,Work,Nwork,Iwork,    &
     &                  Info)
      IMPLICIT NONE
!*--SDRVVX524
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldlre , Ldvl , Ldvr , Niunit , Nounit ,      &
     &        Nsizes , Ntypes , Nwork
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iseed(4) , Iwork(*) , Nn(*)
      REAL A(Lda,*) , H(Lda,*) , Lre(Ldlre,*) , Rcdein(*) , Rcdvin(*) , &
     &     Rcnde1(*) , Rcndv1(*) , Rconde(*) , Rcondv(*) , Result(11) , &
     &     Scale(*) , Scale1(*) , Vl(Ldvl,*) , Vr(Ldvr,*) , Wi(*) ,     &
     &     Wi1(*) , Work(*) , Wr(*) , Wr1(*)
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
      CHARACTER balanc
      CHARACTER*3 path
      INTEGER i , ibal , iinfo , imode , itype , iwk , j , jcol ,       &
     &        jsize , jtype , mtypes , n , nerrs , nfail , nmax ,       &
     &        nnwork , ntest , ntestf , ntestt
      REAL anorm , cond , conds , ovfl , rtulp , rtulpi , ulp , ulpinv ,&
     &     unfl
!     ..
!     .. Local Arrays ..
      CHARACTER adumma(1) , bal(4)
      INTEGER idumma(1) , ioldsd(4) , kconds(MAXTYP) , kmagn(MAXTYP) ,  &
     &        kmode(MAXTYP) , ktype(MAXTYP)
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SLABAD , SLASET 
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , SQRT
!     ..
!     .. Data statements ..
      DATA ktype/1 , 2 , 3 , 5*4 , 4*6 , 6*6 , 3*9/
      DATA kmagn/3*1 , 1 , 1 , 1 , 2 , 3 , 4*1 , 1 , 1 , 1 , 1 , 2 , 3 ,&
     &     1 , 2 , 3/
      DATA kmode/3*0 , 4 , 3 , 1 , 4 , 4 , 4 , 3 , 1 , 5 , 4 , 3 , 1 ,  &
     &     5 , 5 , 5 , 4 , 3 , 1/
      DATA kconds/3*0 , 5*0 , 4*1 , 6*2 , 3*0/
      DATA bal/'N' , 'P' , 'S' , 'B'/
!     ..
!     .. Executable Statements ..
!
      path(1:1) = 'Single precision'
      path(2:3) = 'VX'
!
!     Check for errors
!
      ntestt = 0
      ntestf = 0
      Info = 0
!
!     Important constants
!
      badnn = .FALSE.
!
!     12 is the largest dimension in the input file of precomputed
!     problems
!
      nmax = 12
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
      ELSEIF ( Lda<1 .OR. Lda<nmax ) THEN
         Info = -10
      ELSEIF ( Ldvl<1 .OR. Ldvl<nmax ) THEN
         Info = -17
      ELSEIF ( Ldvr<1 .OR. Ldvr<nmax ) THEN
         Info = -19
      ELSEIF ( Ldlre<1 .OR. Ldlre<nmax ) THEN
         Info = -21
      ELSEIF ( 6*nmax+2*nmax**2>Nwork ) THEN
         Info = -32
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL xerbla('sdrvvx',-Info)
         RETURN
      ENDIF
!
!     If nothing to do check on NIUNIT
!
      IF ( Nsizes/=0 .AND. Ntypes/=0 ) THEN
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
            IF ( Nsizes/=1 ) THEN
               mtypes = MIN(MAXTYP,Ntypes)
            ELSE
               mtypes = MIN(MAXTYP+1,Ntypes)
            ENDIF
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
                        CALL slatms(n,n,'S',Iseed,'S',Work,imode,cond,  &
     &                              anorm,0,0,'N',A,Lda,Work(n+1),iinfo)
!
                     ELSEIF ( itype==5 ) THEN
!
!              Symmetric, eigenvalues specified
!
                        CALL slatms(n,n,'S',Iseed,'S',Work,imode,cond,  &
     &                              anorm,n,n,'N',A,Lda,Work(n+1),iinfo)
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
                        CALL slatme(n,'S',Iseed,Work,imode,cond,ONE,    &
     &                              adumma,'T','T','T',Work(n+1),4,     &
     &                              conds,n,n,anorm,A,Lda,Work(2*n+1),  &
     &                              iinfo)
!
                     ELSEIF ( itype==7 ) THEN
!
!              Diagonal, random eigenvalues
!
                        CALL slatmr(n,n,'S',Iseed,'S',Work,6,ONE,ONE,   &
     &                              'T','N',Work(n+1),1,ONE,Work(2*n+1),&
     &                              1,ONE,'N',idumma,0,0,ZERO,anorm,    &
     &                              'NO',A,Lda,Iwork,iinfo)
!
                     ELSEIF ( itype==8 ) THEN
!
!              Symmetric, random eigenvalues
!
                        CALL slatmr(n,n,'S',Iseed,'S',Work,6,ONE,ONE,   &
     &                              'T','N',Work(n+1),1,ONE,Work(2*n+1),&
     &                              1,ONE,'N',idumma,n,n,ZERO,anorm,    &
     &                              'NO',A,Lda,Iwork,iinfo)
!
                     ELSEIF ( itype==9 ) THEN
!
!              General, random eigenvalues
!
                        CALL slatmr(n,n,'S',Iseed,'N',Work,6,ONE,ONE,   &
     &                              'T','N',Work(n+1),1,ONE,Work(2*n+1),&
     &                              1,ONE,'N',idumma,n,n,ZERO,anorm,    &
     &                              'NO',A,Lda,Iwork,iinfo)
                        IF ( n>=4 ) THEN
                           CALL SLASET('Full',2,n,ZERO,ZERO,A,Lda)
                           CALL SLASET('Full',n-3,1,ZERO,ZERO,A(3,1),   &
     &                                 Lda)
                           CALL SLASET('Full',n-3,2,ZERO,ZERO,A(3,n-1), &
     &                                 Lda)
                           CALL SLASET('Full',1,n,ZERO,ZERO,A(n,1),Lda)
                        ENDIF
!
                     ELSEIF ( itype==10 ) THEN
!
!              Triangular, random eigenvalues
!
                        CALL slatmr(n,n,'S',Iseed,'N',Work,6,ONE,ONE,   &
     &                              'T','N',Work(n+1),1,ONE,Work(2*n+1),&
     &                              1,ONE,'N',idumma,n,0,ZERO,anorm,    &
     &                              'NO',A,Lda,Iwork,iinfo)
!
                     ELSE
!
                        iinfo = 1
                     ENDIF
!
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99008) 'Generator' , iinfo ,  &
     &                         n , jtype , ioldsd
                        Info = ABS(iinfo)
                        RETURN
                     ENDIF
                  ENDIF
!
!
!           Test for minimal and generous workspace
!
                  DO iwk = 1 , 3
                     IF ( iwk==1 ) THEN
                        nnwork = 3*n
                     ELSEIF ( iwk==2 ) THEN
                        nnwork = 6*n + n**2
                     ELSE
                        nnwork = 6*n + 2*n**2
                     ENDIF
                     nnwork = MAX(nnwork,1)
!
!              Test for all balancing options
!
                     DO ibal = 1 , 4
                        balanc = bal(ibal)
!
!                 Perform tests
!
                        CALL sget23(.FALSE.,balanc,jtype,Thresh,ioldsd, &
     &                              Nounit,n,A,Lda,H,Wr,Wi,Wr1,Wi1,Vl,  &
     &                              Ldvl,Vr,Ldvr,Lre,Ldlre,Rcondv,      &
     &                              Rcndv1,Rcdvin,Rconde,Rcnde1,Rcdein, &
     &                              Scale,Scale1,Result,Work,nnwork,    &
     &                              Iwork,Info)
!
!                 Check for RESULT(j) > THRESH
!
                        ntest = 0
                        nfail = 0
                        DO j = 1 , 9
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
                           ntestf = 2
                        ENDIF
!
                        DO j = 1 , 9
                           IF ( Result(j)>=Thresh )                     &
     &                          WRITE (Nounit,FMT=99006) balanc , n ,   &
     &                          iwk , ioldsd , jtype , j , Result(j)
                        ENDDO
!
                        nerrs = nerrs + nfail
                        ntestt = ntestt + ntest
!
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
!
!     Read in data from file to check accuracy of condition estimation.
!     Assume input eigenvalues are sorted lexicographically (increasing
!     by real part, then decreasing by imaginary part)
!
      jtype = 0
      DO
         READ (Niunit,FMT=*,END=100) n
!
!     Read input data until N=0
!
         IF ( n==0 ) EXIT
         jtype = jtype + 1
         Iseed(1) = jtype
         DO i = 1 , n
            READ (Niunit,FMT=*) (A(i,j),j=1,n)
         ENDDO
         DO i = 1 , n
            READ (Niunit,FMT=*) Wr1(i) , Wi1(i) , Rcdein(i) , Rcdvin(i)
         ENDDO
         CALL sget23(.TRUE.,'N',22,Thresh,Iseed,Nounit,n,A,Lda,H,Wr,Wi, &
     &               Wr1,Wi1,Vl,Ldvl,Vr,Ldvr,Lre,Ldlre,Rcondv,Rcndv1,   &
     &               Rcdvin,Rconde,Rcnde1,Rcdein,Scale,Scale1,Result,   &
     &               Work,6*n+2*n**2,Iwork,Info)
!
!     Check for RESULT(j) > THRESH
!
         ntest = 0
         nfail = 0
         DO j = 1 , 11
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
            ntestf = 2
         ENDIF
!
         DO j = 1 , 11
            IF ( Result(j)>=Thresh ) WRITE (Nounit,FMT=99007) n ,       &
     &           jtype , j , Result(j)
         ENDDO
!
         nerrs = nerrs + nfail
         ntestt = ntestt + ntest
      ENDDO
!
!     Summary
!
 100  CALL slasum(path,Nounit,nerrs,ntestt)
!
99001 FORMAT (/1X,A3,' -- Real Eigenvalue-Eigenvector Decomposition',   &
     &        ' Expert Driver',                                         &
     &        /' Matrix types (see sdrvvx for details): ')
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
     &        /' 12=Well-cond., random complex ','         ',           &
     &        ' 17=Ill-cond., large rand. complx ',/' 13=Ill-condi',    &
     &        'tioned, evenly spaced.     ',                            &
     &        ' 18=Ill-cond., small rand.',' complx ')
99004 FORMAT (' 19=Matrix with random O(1) entries.    ',' 21=Matrix ', &
     &        'with small random entries.',/' 20=Matrix with large ran',&
     &        'dom entries.   ',' 22=Matrix read from input file',/)
99005 FORMAT (' Tests performed with test threshold =',F8.2,            &
     &        //' 1 = | A VR - VR W | / ( n |A| ulp ) ',                &
     &        /' 2 = | transpose(A) VL - VL W | / ( n |A| ulp ) ',      &
     &        /' 3 = | |VR(i)| - 1 | / ulp ',                           &
     &        /' 4 = | |VL(i)| - 1 | / ulp ',                           &
     &        /' 5 = 0 if W same no matter if VR or VL computed,',      &
     &        ' 1/ulp otherwise',                                       &
     &        /' 6 = 0 if VR same no matter what else computed,',       &
     &        '  1/ulp otherwise',                                      &
     &        /' 7 = 0 if VL same no matter what else computed,',       &
     &        '  1/ulp otherwise',                                      &
     &        /' 8 = 0 if RCONDV same no matter what else computed,',   &
     &        '  1/ulp otherwise',/                                     &
     &       ' 9 = 0 if SCALE, ILO, IHI, ABNRM same no matter what else'&
     &       ,' computed,  1/ulp otherwise',                            &
     &       /' 10 = | RCONDV - RCONDV(precomputed) | / cond(RCONDV),', &
     &       /' 11 = | RCONDE - RCONDE(precomputed) | / cond(RCONDE),')
99006 FORMAT (' BALANC=''',A1,''',N=',I4,',IWK=',I1,', seed=',4(I4,','),&
     &        ' type ',I2,', test(',I2,')=',G10.3)
99007 FORMAT (' N=',I5,', input example =',I3,',  test(',I2,')=',G10.3)
99008 FORMAT (' sdrvvx: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
!
!     End of sdrvvx
!
      END SUBROUTINE sdrvvx
