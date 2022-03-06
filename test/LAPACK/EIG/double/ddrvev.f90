!*==ddrvev.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DDRVEV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DDRVEV( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          NOUNIT, A, LDA, H, WR, WI, WR1, WI1, VL, LDVL,
!                          VR, LDVR, LRE, LDLRE, RESULT, WORK, NWORK,
!                          IWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDLRE, LDVL, LDVR, NOUNIT, NSIZES,
!      $                   NTYPES, NWORK
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
!       DOUBLE PRECISION   A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ),
!      $                   RESULT( 7 ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WI( * ), WI1( * ), WORK( * ), WR( * ), WR1( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DDRVEV  checks the nonsymmetric eigenvalue problem driver DGEEV.
!>
!>    When DDRVEV is called, a number of matrix "sizes" ("n's") and a
!>    number of matrix "types" are specified.  For each size ("n")
!>    and each type of matrix, one matrix will be generated and used
!>    to test the nonsymmetric eigenroutines.  For each matrix, 7
!>    tests will be performed:
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
!>    (3)     | |VR(i)| - 1 | / ulp and whether largest component real
!>
!>      VR(i) denotes the i-th column of VR.
!>
!>    (4)     | |VL(i)| - 1 | / ulp and whether largest component real
!>
!>      VL(i) denotes the i-th column of VL.
!>
!>    (5)     W(full) = W(partial)
!>
!>      W(full) denotes the eigenvalues computed when both VR and VL
!>      are also computed, and W(partial) denotes the eigenvalues
!>      computed when only W, only W and VR, or only W and VL are
!>      computed.
!>
!>    (6)     VR(full) = VR(partial)
!>
!>      VR(full) denotes the right eigenvectors computed when both VR
!>      and VL are computed, and VR(partial) denotes the result
!>      when only VR is computed.
!>
!>     (7)     VL(full) = VL(partial)
!>
!>      VL(full) denotes the left eigenvectors computed when both VR
!>      and VL are also computed, and VL(partial) denotes the result
!>      when only VL is computed.
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
!>          DDRVEV does nothing.  It must be at least zero.
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
!>          The number of elements in DOTYPE.   If it is zero, DDRVEV
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
!>          next call to DDRVEV to continue the same random number
!>          sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
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
!>          A is DOUBLE PRECISION array, dimension (LDA, max(NN))
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
!>          H is DOUBLE PRECISION array, dimension (LDA, max(NN))
!>          Another copy of the test matrix A, modified by DGEEV.
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (max(NN))
!>
!>          The real and imaginary parts of the eigenvalues of A.
!>          On exit, WR + WI*i are the eigenvalues of the matrix in A.
!> \endverbatim
!>
!> \param[out] WR1
!> \verbatim
!>          WR1 is DOUBLE PRECISION array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] WI1
!> \verbatim
!>          WI1 is DOUBLE PRECISION array, dimension (max(NN))
!>
!>          Like WR, WI, these arrays contain the eigenvalues of A,
!>          but those computed when DGEEV only computes a partial
!>          eigendecomposition, i.e. not the eigenvalues and left
!>          and right eigenvectors.
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is DOUBLE PRECISION array, dimension (LDVL, max(NN))
!>          VL holds the computed left eigenvectors.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          Leading dimension of VL. Must be at least max(1,max(NN)).
!> \endverbatim
!>
!> \param[out] VR
!> \verbatim
!>          VR is DOUBLE PRECISION array, dimension (LDVR, max(NN))
!>          VR holds the computed right eigenvectors.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          Leading dimension of VR. Must be at least max(1,max(NN)).
!> \endverbatim
!>
!> \param[out] LRE
!> \verbatim
!>          LRE is DOUBLE PRECISION array, dimension (LDLRE,max(NN))
!>          LRE holds the computed right or left eigenvectors.
!> \endverbatim
!>
!> \param[in] LDLRE
!> \verbatim
!>          LDLRE is INTEGER
!>          Leading dimension of LRE. Must be at least max(1,max(NN)).
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (7)
!>          The values computed by the seven tests described above.
!>          The values are currently limited to 1/ulp, to avoid overflow.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (NWORK)
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
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0, then everything ran OK.
!>           -1: NSIZES < 0
!>           -2: Some NN(j) < 0
!>           -3: NTYPES < 0
!>           -6: THRESH < 0
!>           -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
!>          -16: LDVL < 1 or LDVL < NMAX, where NMAX is max( NN(j) ).
!>          -18: LDVR < 1 or LDVR < NMAX, where NMAX is max( NN(j) ).
!>          -20: LDLRE < 1 or LDLRE < NMAX, where NMAX is max( NN(j) ).
!>          -23: NWORK too small.
!>          If  DLATMR, SLATMS, SLATME or DGEEV returns an error code,
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
!> \date December 2016
!
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DDRVEV(Nsizes,Nn,Ntypes,Dotype,Iseed,Thresh,Nounit,A,  &
     &                  Lda,H,Wr,Wi,Wr1,Wi1,Vl,Ldvl,Vr,Ldvr,Lre,Ldlre,  &
     &                  Result,Work,Nwork,Iwork,Info)
      IMPLICIT NONE
!*--DDRVEV409
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldlre , Ldvl , Ldvr , Nounit , Nsizes ,      &
     &        Ntypes , Nwork
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iseed(4) , Iwork(*) , Nn(*)
      DOUBLE PRECISION A(Lda,*) , H(Lda,*) , Lre(Ldlre,*) , Result(7) , &
     &                 Vl(Ldvl,*) , Vr(Ldvr,*) , Wi(*) , Wi1(*) ,       &
     &                 Work(*) , Wr(*) , Wr1(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
      INTEGER MAXTYP
      PARAMETER (MAXTYP=21)
!     ..
!     .. Local Scalars ..
      LOGICAL badnn
      CHARACTER*3 path
      INTEGER iinfo , imode , itype , iwk , j , jcol , jj , jsize ,     &
     &        jtype , mtypes , n , nerrs , nfail , nmax , nnwork ,      &
     &        ntest , ntestf , ntestt
      DOUBLE PRECISION anorm , cond , conds , ovfl , rtulp , rtulpi ,   &
     &                 tnrm , ulp , ulpinv , unfl , vmx , vrmx , vtst
!     ..
!     .. Local Arrays ..
      CHARACTER adumma(1)
      INTEGER idumma(1) , ioldsd(4) , kconds(MAXTYP) , kmagn(MAXTYP) ,  &
     &        kmode(MAXTYP) , ktype(MAXTYP)
      DOUBLE PRECISION dum(1) , res(2)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLAPY2 , DNRM2
      EXTERNAL DLAMCH , DLAPY2 , DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEEV , DGET22 , DLABAD , DLACPY , DLASET , DLASUM ,     &
     &         DLATME , DLATMR , DLATMS , XERBLA
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
!     ..
!     .. Executable Statements ..
!
      path(1:1) = 'Double precision'
      path(2:3) = 'EV'
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
      ELSEIF ( Ldvl<1 .OR. Ldvl<nmax ) THEN
         Info = -16
      ELSEIF ( Ldvr<1 .OR. Ldvr<nmax ) THEN
         Info = -18
      ELSEIF ( Ldlre<1 .OR. Ldlre<nmax ) THEN
         Info = -20
      ELSEIF ( 5*nmax+2*nmax**2>Nwork ) THEN
         Info = -23
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DDRVEV',-Info)
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
      ovfl = ONE/unfl
      CALL DLABAD(unfl,ovfl)
      ulp = DLAMCH('Precision')
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
                  CALL DLASET('Full',Lda,n,ZERO,ZERO,A,Lda)
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
                     CALL DLATME(n,'S',Iseed,Work,imode,cond,ONE,adumma,&
     &                           'T','T','T',Work(n+1),4,conds,n,n,     &
     &                           anorm,A,Lda,Work(2*n+1),iinfo)
!
                  ELSEIF ( itype==7 ) THEN
!
!              Diagonal, random eigenvalues
!
                     CALL DLATMR(n,n,'S',Iseed,'S',Work,6,ONE,ONE,'T',  &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,0,0,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==8 ) THEN
!
!              Symmetric, random eigenvalues
!
                     CALL DLATMR(n,n,'S',Iseed,'S',Work,6,ONE,ONE,'T',  &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,n,n,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==9 ) THEN
!
!              General, random eigenvalues
!
                     CALL DLATMR(n,n,'S',Iseed,'N',Work,6,ONE,ONE,'T',  &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,n,n,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
                     IF ( n>=4 ) THEN
                        CALL DLASET('Full',2,n,ZERO,ZERO,A,Lda)
                        CALL DLASET('Full',n-3,1,ZERO,ZERO,A(3,1),Lda)
                        CALL DLASET('Full',n-3,2,ZERO,ZERO,A(3,n-1),Lda)
                        CALL DLASET('Full',1,n,ZERO,ZERO,A(n,1),Lda)
                     ENDIF
!
                  ELSEIF ( itype==10 ) THEN
!
!              Triangular, random eigenvalues
!
                     CALL DLATMR(n,n,'S',Iseed,'N',Work,6,ONE,ONE,'T',  &
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
                     WRITE (Nounit,FMT=99007) 'Generator' , iinfo , n , &
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
                     nnwork = 4*n
                  ELSE
                     nnwork = 5*n + 2*n**2
                  ENDIF
                  nnwork = MAX(nnwork,1)
!
!              Initialize RESULT
!
                  DO j = 1 , 7
                     Result(j) = -ONE
                  ENDDO
!
!              Compute eigenvalues and eigenvectors, and test them
!
                  CALL DLACPY('F',n,n,A,Lda,H,Lda)
                  CALL DGEEV('V','V',n,H,Lda,Wr,Wi,Vl,Ldvl,Vr,Ldvr,Work,&
     &                       nnwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     Result(1) = ulpinv
                     WRITE (Nounit,FMT=99007) 'DGEEV1' , iinfo , n ,    &
     &                      jtype , ioldsd
                     Info = ABS(iinfo)
                     GOTO 5
                  ENDIF
!
!              Do Test (1)
!
                  CALL DGET22('N','N','N',n,A,Lda,Vr,Ldvr,Wr,Wi,Work,   &
     &                        res)
                  Result(1) = res(1)
!
!              Do Test (2)
!
                  CALL DGET22('T','N','T',n,A,Lda,Vl,Ldvl,Wr,Wi,Work,   &
     &                        res)
                  Result(2) = res(1)
!
!              Do Test (3)
!
                  DO j = 1 , n
                     tnrm = ONE
                     IF ( Wi(j)==ZERO ) THEN
                        tnrm = DNRM2(n,Vr(1,j),1)
                     ELSEIF ( Wi(j)>ZERO ) THEN
                        tnrm = DLAPY2(DNRM2(n,Vr(1,j),1),               &
     &                         DNRM2(n,Vr(1,j+1),1))
                     ENDIF
                     Result(3) = MAX(Result(3),MIN(ulpinv,ABS(tnrm-ONE)/&
     &                           ulp))
                     IF ( Wi(j)>ZERO ) THEN
                        vmx = ZERO
                        vrmx = ZERO
                        DO jj = 1 , n
                           vtst = DLAPY2(Vr(jj,j),Vr(jj,j+1))
                           IF ( vtst>vmx ) vmx = vtst
                           IF ( Vr(jj,j+1)==ZERO .AND. ABS(Vr(jj,j))    &
     &                          >vrmx ) vrmx = ABS(Vr(jj,j))
                        ENDDO
                        IF ( vrmx/vmx<ONE-TWO*ulp ) Result(3) = ulpinv
                     ENDIF
                  ENDDO
!
!              Do Test (4)
!
                  DO j = 1 , n
                     tnrm = ONE
                     IF ( Wi(j)==ZERO ) THEN
                        tnrm = DNRM2(n,Vl(1,j),1)
                     ELSEIF ( Wi(j)>ZERO ) THEN
                        tnrm = DLAPY2(DNRM2(n,Vl(1,j),1),               &
     &                         DNRM2(n,Vl(1,j+1),1))
                     ENDIF
                     Result(4) = MAX(Result(4),MIN(ulpinv,ABS(tnrm-ONE)/&
     &                           ulp))
                     IF ( Wi(j)>ZERO ) THEN
                        vmx = ZERO
                        vrmx = ZERO
                        DO jj = 1 , n
                           vtst = DLAPY2(Vl(jj,j),Vl(jj,j+1))
                           IF ( vtst>vmx ) vmx = vtst
                           IF ( Vl(jj,j+1)==ZERO .AND. ABS(Vl(jj,j))    &
     &                          >vrmx ) vrmx = ABS(Vl(jj,j))
                        ENDDO
                        IF ( vrmx/vmx<ONE-TWO*ulp ) Result(4) = ulpinv
                     ENDIF
                  ENDDO
!
!              Compute eigenvalues only, and test them
!
                  CALL DLACPY('F',n,n,A,Lda,H,Lda)
                  CALL DGEEV('N','N',n,H,Lda,Wr1,Wi1,dum,1,dum,1,Work,  &
     &                       nnwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     Result(1) = ulpinv
                     WRITE (Nounit,FMT=99007) 'DGEEV2' , iinfo , n ,    &
     &                      jtype , ioldsd
                     Info = ABS(iinfo)
                     GOTO 5
                  ENDIF
!
!              Do Test (5)
!
                  DO j = 1 , n
                     IF ( Wr(j)/=Wr1(j) .OR. Wi(j)/=Wi1(j) ) Result(5)  &
     &                    = ulpinv
                  ENDDO
!
!              Compute eigenvalues and right eigenvectors, and test them
!
                  CALL DLACPY('F',n,n,A,Lda,H,Lda)
                  CALL DGEEV('N','V',n,H,Lda,Wr1,Wi1,dum,1,Lre,Ldlre,   &
     &                       Work,nnwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     Result(1) = ulpinv
                     WRITE (Nounit,FMT=99007) 'DGEEV3' , iinfo , n ,    &
     &                      jtype , ioldsd
                     Info = ABS(iinfo)
                     GOTO 5
                  ENDIF
!
!              Do Test (5) again
!
                  DO j = 1 , n
                     IF ( Wr(j)/=Wr1(j) .OR. Wi(j)/=Wi1(j) ) Result(5)  &
     &                    = ulpinv
                  ENDDO
!
!              Do Test (6)
!
                  DO j = 1 , n
                     DO jj = 1 , n
                        IF ( Vr(j,jj)/=Lre(j,jj) ) Result(6) = ulpinv
                     ENDDO
                  ENDDO
!
!              Compute eigenvalues and left eigenvectors, and test them
!
                  CALL DLACPY('F',n,n,A,Lda,H,Lda)
                  CALL DGEEV('V','N',n,H,Lda,Wr1,Wi1,Lre,Ldlre,dum,1,   &
     &                       Work,nnwork,iinfo)
                  IF ( iinfo/=0 ) THEN
                     Result(1) = ulpinv
                     WRITE (Nounit,FMT=99007) 'DGEEV4' , iinfo , n ,    &
     &                      jtype , ioldsd
                     Info = ABS(iinfo)
                     GOTO 5
                  ENDIF
!
!              Do Test (5) again
!
                  DO j = 1 , n
                     IF ( Wr(j)/=Wr1(j) .OR. Wi(j)/=Wi1(j) ) Result(5)  &
     &                    = ulpinv
                  ENDDO
!
!              Do Test (7)
!
                  DO j = 1 , n
                     DO jj = 1 , n
                        IF ( Vl(j,jj)/=Lre(j,jj) ) Result(7) = ulpinv
                     ENDDO
                  ENDDO
!
!              End of Loop -- Check for RESULT(j) > THRESH
!
!
 5                ntest = 0
                  nfail = 0
                  DO j = 1 , 7
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
                  DO j = 1 , 7
                     IF ( Result(j)>=Thresh ) WRITE (Nounit,FMT=99006)  &
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
      CALL DLASUM(path,Nounit,nerrs,ntestt)
!
99001 FORMAT (/1X,A3,' -- Real Eigenvalue-Eigenvector Decomposition',   &
     &        ' Driver',/' Matrix types (see DDRVEV for details): ')
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
     &        //' 1 = | A VR - VR W | / ( n |A| ulp ) ',                &
     &        /' 2 = | transpose(A) VL - VL W | / ( n |A| ulp ) ',      &
     &        /' 3 = | |VR(i)| - 1 | / ulp ',                           &
     &        /' 4 = | |VL(i)| - 1 | / ulp ',                           &
     &        /' 5 = 0 if W same no matter if VR or VL computed,',      &
     &        ' 1/ulp otherwise',                                       &
     &        /' 6 = 0 if VR same no matter if VL computed,',           &
     &        '  1/ulp otherwise',                                      &
     &        /' 7 = 0 if VL same no matter if VR computed,',           &
     &        '  1/ulp otherwise',/)
99006 FORMAT (' N=',I5,', IWK=',I2,', seed=',4(I4,','),' type ',I2,     &
     &        ', test(',I2,')=',G10.3)
99007 FORMAT (' DDRVEV: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
!
!     End of DDRVEV
!
      END SUBROUTINE DDRVEV
