!*==ddrvsg2stg.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b ddrvsg2stg
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ddrvsg2stg( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                              NOUNIT, A, LDA, B, LDB, D, D2, Z, LDZ, AB,
!                              BB, AP, BP, WORK, NWORK, IWORK, LIWORK,
!                              RESULT, INFO )
!
!       IMPLICIT NONE
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LDZ, LIWORK, NOUNIT, NSIZES,
!      $                   NTYPES, NWORK
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), NN( * )
!       DOUBLE PRECISION   A( LDA, * ), AB( LDA, * ), AP( * ),
!      $                   B( LDB, * ), BB( LDB, * ), BP( * ), D( * ),
!      $                   RESULT( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      ddrvsg2stg checks the real symmetric generalized eigenproblem
!>      drivers.
!>
!>              DSYGV computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric-definite generalized
!>              eigenproblem.
!>
!>              DSYGVD computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric-definite generalized
!>              eigenproblem using a divide and conquer algorithm.
!>
!>              DSYGVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric-definite generalized
!>              eigenproblem.
!>
!>              DSPGV computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric-definite generalized
!>              eigenproblem in packed storage.
!>
!>              DSPGVD computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric-definite generalized
!>              eigenproblem in packed storage using a divide and
!>              conquer algorithm.
!>
!>              DSPGVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric-definite generalized
!>              eigenproblem in packed storage.
!>
!>              DSBGV computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric-definite banded
!>              generalized eigenproblem.
!>
!>              DSBGVD computes all eigenvalues and, optionally,
!>              eigenvectors of a real symmetric-definite banded
!>              generalized eigenproblem using a divide and conquer
!>              algorithm.
!>
!>              DSBGVX computes selected eigenvalues and, optionally,
!>              eigenvectors of a real symmetric-definite banded
!>              generalized eigenproblem.
!>
!>      When ddrvsg2stg is called, a number of matrix "sizes" ("n's") and a
!>      number of matrix "types" are specified.  For each size ("n")
!>      and each type of matrix, one matrix A of the given type will be
!>      generated; a random well-conditioned matrix B is also generated
!>      and the pair (A,B) is used to test the drivers.
!>
!>      For each pair (A,B), the following tests are performed:
!>
!>      (1) DSYGV with ITYPE = 1 and UPLO ='U':
!>
!>              | A Z - B Z D | / ( |A| |Z| n ulp )
!>              | D - D2 | / ( |D| ulp )   where D is computed by
!>                                         DSYGV and  D2 is computed by
!>                                         DSYGV_2STAGE. This test is
!>                                         only performed for DSYGV
!>
!>      (2) as (1) but calling DSPGV
!>      (3) as (1) but calling DSBGV
!>      (4) as (1) but with UPLO = 'L'
!>      (5) as (4) but calling DSPGV
!>      (6) as (4) but calling DSBGV
!>
!>      (7) DSYGV with ITYPE = 2 and UPLO ='U':
!>
!>              | A B Z - Z D | / ( |A| |Z| n ulp )
!>
!>      (8) as (7) but calling DSPGV
!>      (9) as (7) but with UPLO = 'L'
!>      (10) as (9) but calling DSPGV
!>
!>      (11) DSYGV with ITYPE = 3 and UPLO ='U':
!>
!>              | B A Z - Z D | / ( |A| |Z| n ulp )
!>
!>      (12) as (11) but calling DSPGV
!>      (13) as (11) but with UPLO = 'L'
!>      (14) as (13) but calling DSPGV
!>
!>      DSYGVD, DSPGVD and DSBGVD performed the same 14 tests.
!>
!>      DSYGVX, DSPGVX and DSBGVX performed the above 14 tests with
!>      the parameter RANGE = 'A', 'N' and 'I', respectively.
!>
!>      The "sizes" are specified by an array NN(1:NSIZES); the value
!>      of each element NN(j) specifies one size.
!>      The "types" are specified by a logical array DOTYPE( 1:NTYPES );
!>      if DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!>      This type is used for the matrix A which has half-bandwidth KA.
!>      B is generated as a well-conditioned positive definite matrix
!>      with half-bandwidth KB (<= KA).
!>      Currently, the list of possible types for A is:
!>
!>      (1)  The zero matrix.
!>      (2)  The identity matrix.
!>
!>      (3)  A diagonal matrix with evenly spaced entries
!>           1, ..., ULP  and random signs.
!>           (ULP = (first number larger than 1) - 1 )
!>      (4)  A diagonal matrix with geometrically spaced entries
!>           1, ..., ULP  and random signs.
!>      (5)  A diagonal matrix with "clustered" entries
!>           1, ULP, ..., ULP and random signs.
!>
!>      (6)  Same as (4), but multiplied by SQRT( overflow threshold )
!>      (7)  Same as (4), but multiplied by SQRT( underflow threshold )
!>
!>      (8)  A matrix of the form  U* D U, where U is orthogonal and
!>           D has evenly spaced entries 1, ..., ULP with random signs
!>           on the diagonal.
!>
!>      (9)  A matrix of the form  U* D U, where U is orthogonal and
!>           D has geometrically spaced entries 1, ..., ULP with random
!>           signs on the diagonal.
!>
!>      (10) A matrix of the form  U* D U, where U is orthogonal and
!>           D has "clustered" entries 1, ULP,..., ULP with random
!>           signs on the diagonal.
!>
!>      (11) Same as (8), but multiplied by SQRT( overflow threshold )
!>      (12) Same as (8), but multiplied by SQRT( underflow threshold )
!>
!>      (13) symmetric matrix with random entries chosen from (-1,1).
!>      (14) Same as (13), but multiplied by SQRT( overflow threshold )
!>      (15) Same as (13), but multiplied by SQRT( underflow threshold)
!>
!>      (16) Same as (8), but with KA = 1 and KB = 1
!>      (17) Same as (8), but with KA = 2 and KB = 1
!>      (18) Same as (8), but with KA = 2 and KB = 2
!>      (19) Same as (8), but with KA = 3 and KB = 1
!>      (20) Same as (8), but with KA = 3 and KB = 2
!>      (21) Same as (8), but with KA = 3 and KB = 3
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  NSIZES  INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          ddrvsg2stg does nothing.  It must be at least zero.
!>          Not modified.
!>
!>  NN      INTEGER array, dimension (NSIZES)
!>          An array containing the sizes to be used for the matrices.
!>          Zero values will be skipped.  The values must be at least
!>          zero.
!>          Not modified.
!>
!>  NTYPES  INTEGER
!>          The number of elements in DOTYPE.   If it is zero, ddrvsg2stg
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
!>          next call to ddrvsg2stg to continue the same random number
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
!>          The leading dimension of A and AB.  It must be at
!>          least 1 and at least max( NN ).
!>          Not modified.
!>
!>  B       DOUBLE PRECISION array, dimension (LDB , max(NN))
!>          Used to hold the symmetric positive definite matrix for
!>          the generailzed problem.
!>          On exit, B contains the last matrix actually
!>          used.
!>          Modified.
!>
!>  LDB     INTEGER
!>          The leading dimension of B and BB.  It must be at
!>          least 1 and at least max( NN ).
!>          Not modified.
!>
!>  D       DOUBLE PRECISION array, dimension (max(NN))
!>          The eigenvalues of A. On exit, the eigenvalues in D
!>          correspond with the matrix in A.
!>          Modified.
!>
!>  Z       DOUBLE PRECISION array, dimension (LDZ, max(NN))
!>          The matrix of eigenvectors.
!>          Modified.
!>
!>  LDZ     INTEGER
!>          The leading dimension of Z.  It must be at least 1 and
!>          at least max( NN ).
!>          Not modified.
!>
!>  AB      DOUBLE PRECISION array, dimension (LDA, max(NN))
!>          Workspace.
!>          Modified.
!>
!>  BB      DOUBLE PRECISION array, dimension (LDB, max(NN))
!>          Workspace.
!>          Modified.
!>
!>  AP      DOUBLE PRECISION array, dimension (max(NN)**2)
!>          Workspace.
!>          Modified.
!>
!>  BP      DOUBLE PRECISION array, dimension (max(NN)**2)
!>          Workspace.
!>          Modified.
!>
!>  WORK    DOUBLE PRECISION array, dimension (NWORK)
!>          Workspace.
!>          Modified.
!>
!>  NWORK   INTEGER
!>          The number of entries in WORK.  This must be at least
!>          1+5*N+2*N*lg(N)+3*N**2 where N = max( NN(j) ) and
!>          lg( N ) = smallest integer k such that 2**k >= N.
!>          Not modified.
!>
!>  IWORK   INTEGER array, dimension (LIWORK)
!>          Workspace.
!>          Modified.
!>
!>  LIWORK  INTEGER
!>          The number of entries in WORK.  This must be at least 6*N.
!>          Not modified.
!>
!>  RESULT  DOUBLE PRECISION array, dimension (70)
!>          The values computed by the 70 tests described above.
!>          Modified.
!>
!>  INFO    INTEGER
!>          If 0, then everything ran OK.
!>           -1: NSIZES < 0
!>           -2: Some NN(j) < 0
!>           -3: NTYPES < 0
!>           -5: THRESH < 0
!>           -9: LDA < 1 or LDA < NMAX, where NMAX is max( NN(j) ).
!>          -16: LDZ < 1 or LDZ < NMAX.
!>          -21: NWORK too small.
!>          -23: LIWORK too small.
!>          If  dlatmr, slatms, DSYGV, DSPGV, DSBGV, SSYGVD, SSPGVD,
!>              DSBGVD, DSYGVX, DSPGVX or SSBGVX returns an error code,
!>              the absolute value of it is returned.
!>          Modified.
!>
!> ----------------------------------------------------------------------
!>
!>       Some Local Variables and Parameters:
!>       ---- ----- --------- --- ----------
!>       ZERO, ONE       Real 0 and 1.
!>       MAXTYP          The number of types defined.
!>       NTEST           The number of tests that have been run
!>                       on this matrix.
!>       NTESTT          The total number of tests for this call.
!>       NMAX            Largest value in NN.
!>       NMATS           The number of matrices generated so far.
!>       NERRS           The number of tests which have exceeded THRESH
!>                       so far (computed by dlafts).
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE ddrvsg2stg(Nsizes,Nn,Ntypes,Dotype,Iseed,Thresh,Nounit,&
     &                      A,Lda,B,Ldb,D,D2,Z,Ldz,Ab,Bb,Ap,Bp,Work,    &
     &                      Nwork,Iwork,Liwork,Result,Info)
!
      IMPLICIT NONE
!*--DDRVSG2STG366
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldb , Ldz , Liwork , Nounit , Nsizes ,       &
     &        Ntypes , Nwork
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iseed(4) , Iwork(*) , Nn(*)
      DOUBLE PRECISION A(Lda,*) , Ab(Lda,*) , Ap(*) , B(Ldb,*) ,        &
     &                 Bb(Ldb,*) , Bp(*) , D(*) , D2(*) , Result(*) ,   &
     &                 Work(*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TEN
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TEN=10.0D0)
      INTEGER MAXTYP
      PARAMETER (MAXTYP=21)
!     ..
!     .. Local Scalars ..
      LOGICAL badnn
      CHARACTER uplo
      INTEGER i , ibtype , ibuplo , iinfo , ij , il , imode , itemp ,   &
     &        itype , iu , j , jcol , jsize , jtype , ka , ka9 , kb ,   &
     &        kb9 , m , mtypes , n , nerrs , nmats , nmax , ntest ,     &
     &        ntestt
      DOUBLE PRECISION abstol , aninv , anorm , cond , ovfl , rtovfl ,  &
     &                 rtunfl , ulp , ulpinv , unfl , vl , vu , temp1 , &
     &                 temp2
!     ..
!     .. Local Arrays ..
      INTEGER idumma(1) , ioldsd(4) , iseed2(4) , kmagn(MAXTYP) ,       &
     &        kmode(MAXTYP) , ktype(MAXTYP)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH 
      EXTERNAL LSAME , DLAMCH 
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD , DLACPY , DLASET 
     EXTERNAL         DSBGV , DSBGVD , DSBGVX , DSPGV 
     EXTERNAL         DSPGVD , DSPGVX , DSYGV , DSYGVD , DSYGVX 
     EXTERNAL         DSYGV_2STAGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN , SQRT
!     ..
!     .. Data statements ..
      DATA ktype/1 , 2 , 5*4 , 5*5 , 3*8 , 6*9/
      DATA kmagn/2*1 , 1 , 1 , 1 , 2 , 3 , 1 , 1 , 1 , 2 , 3 , 1 , 2 ,  &
     &     3 , 6*1/
      DATA kmode/2*0 , 4 , 3 , 1 , 4 , 4 , 4 , 3 , 1 , 4 , 4 , 0 , 0 ,  &
     &     0 , 6*4/
!     ..
!     .. Executable Statements ..
!
!     1)      Check for errors
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
      ELSEIF ( Lda<=1 .OR. Lda<nmax ) THEN
         Info = -9
      ELSEIF ( Ldz<=1 .OR. Ldz<nmax ) THEN
         Info = -16
      ELSEIF ( 2*MAX(nmax,3)**2>Nwork ) THEN
         Info = -21
      ELSEIF ( 2*MAX(nmax,3)**2>Liwork ) THEN
         Info = -23
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL xerbla('ddrvsg2stg',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
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
      DO i = 1 , 4
         iseed2(i) = Iseed(i)
      ENDDO
!
!     Loop over sizes, types
!
      nerrs = 0
      nmats = 0
!
      DO jsize = 1 , Nsizes
         n = Nn(jsize)
         aninv = ONE/DBLE(MAX(1,n))
!
         IF ( Nsizes/=1 ) THEN
            mtypes = MIN(MAXTYP,Ntypes)
         ELSE
            mtypes = MIN(MAXTYP+1,Ntypes)
         ENDIF
!
         ka9 = 0
         kb9 = 0
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
!           =4         arithmetic   diagonal, w/ eigenvalues
!           =5         random log   hermitian, w/ eigenvalues
!           =6         random       (none)
!           =7                      random diagonal
!           =8                      random hermitian
!           =9                      banded, w/ eigenvalues
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
                  iinfo = 0
                  cond = ulpinv
!
!           Special Matrices -- Identity & Jordan block
!
                  IF ( itype==1 ) THEN
!
!              Zero
!
                     ka = 0
                     kb = 0
                     CALL DLASET('Full',Lda,n,ZERO,ZERO,A,Lda)
!
                  ELSEIF ( itype==2 ) THEN
!
!              Identity
!
                     ka = 0
                     kb = 0
                     CALL DLASET('Full',Lda,n,ZERO,ZERO,A,Lda)
                     DO jcol = 1 , n
                        A(jcol,jcol) = anorm
                     ENDDO
!
                  ELSEIF ( itype==4 ) THEN
!
!              Diagonal Matrix, [Eigen]values Specified
!
                     ka = 0
                     kb = 0
                     CALL dlatms(n,n,'S',Iseed,'S',Work,imode,cond,     &
     &                           anorm,0,0,'N',A,Lda,Work(n+1),iinfo)
!
                  ELSEIF ( itype==5 ) THEN
!
!              symmetric, eigenvalues specified
!
                     ka = MAX(0,n-1)
                     kb = ka
                     CALL dlatms(n,n,'S',Iseed,'S',Work,imode,cond,     &
     &                           anorm,n,n,'N',A,Lda,Work(n+1),iinfo)
!
                  ELSEIF ( itype==7 ) THEN
!
!              Diagonal, random eigenvalues
!
                     ka = 0
                     kb = 0
                     CALL dlatmr(n,n,'S',Iseed,'S',Work,6,ONE,ONE,'T',  &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,0,0,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==8 ) THEN
!
!              symmetric, random eigenvalues
!
                     ka = MAX(0,n-1)
                     kb = ka
                     CALL dlatmr(n,n,'S',Iseed,'H',Work,6,ONE,ONE,'T',  &
     &                           'N',Work(n+1),1,ONE,Work(2*n+1),1,ONE, &
     &                           'N',idumma,n,n,ZERO,anorm,'NO',A,Lda,  &
     &                           Iwork,iinfo)
!
                  ELSEIF ( itype==9 ) THEN
!
!              symmetric banded, eigenvalues specified
!
!              The following values are used for the half-bandwidths:
!
!                ka = 1   kb = 1
!                ka = 2   kb = 1
!                ka = 2   kb = 2
!                ka = 3   kb = 1
!                ka = 3   kb = 2
!                ka = 3   kb = 3
!
                     kb9 = kb9 + 1
                     IF ( kb9>ka9 ) THEN
                        ka9 = ka9 + 1
                        kb9 = 1
                     ENDIF
                     ka = MAX(0,MIN(n-1,ka9))
                     kb = MAX(0,MIN(n-1,kb9))
                     CALL dlatms(n,n,'S',Iseed,'S',Work,imode,cond,     &
     &                           anorm,ka,ka,'N',A,Lda,Work(n+1),iinfo)
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
               abstol = unfl + unfl
               IF ( n<=1 ) THEN
                  il = 1
                  iu = n
               ELSE
                  il = 1 + INT((n-1)*dlarnd(1,iseed2))
                  iu = 1 + INT((n-1)*dlarnd(1,iseed2))
                  IF ( il>iu ) THEN
                     itemp = il
                     il = iu
                     iu = itemp
                  ENDIF
               ENDIF
!
!           3) Call DSYGV, DSPGV, DSBGV, SSYGVD, SSPGVD, SSBGVD,
!              DSYGVX, DSPGVX, and DSBGVX, do tests.
!
!           loop over the three generalized problems
!                 IBTYPE = 1: A*x = (lambda)*B*x
!                 IBTYPE = 2: A*B*x = (lambda)*x
!                 IBTYPE = 3: B*A*x = (lambda)*x
!
               DO ibtype = 1 , 3
!
!              loop over the setting UPLO
!
                  DO ibuplo = 1 , 2
                     IF ( ibuplo==1 ) uplo = 'U'
                     IF ( ibuplo==2 ) uplo = 'L'
!
!                 Generate random well-conditioned positive definite
!                 matrix B, of bandwidth not greater than that of A.
!
                     CALL dlatms(n,n,'U',Iseed,'P',Work,5,TEN,ONE,kb,kb,&
     &                           uplo,B,Ldb,Work(n+1),iinfo)
!
!                 Test DSYGV
!
                     ntest = ntest + 1
!
                     CALL DLACPY(' ',n,n,A,Lda,Z,Ldz)
                     CALL DLACPY(uplo,n,n,B,Ldb,Bb,Ldb)
!
                     CALL DSYGV(ibtype,'V',uplo,n,Z,Ldz,Bb,Ldb,D,Work,  &
     &                          Nwork,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'DSYGV(V,'//uplo//')' ,&
     &                         iinfo , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(ntest) = ulpinv
                           GOTO 2
                        ENDIF
                     ENDIF
!
!                 Do Test
!
                     CALL dsgt01(ibtype,uplo,n,n,A,Lda,B,Ldb,Z,Ldz,D,   &
     &                           Work,Result(ntest))
!
!                 Test DSYGV_2STAGE
!
                     ntest = ntest + 1
!
                     CALL DLACPY(' ',n,n,A,Lda,Z,Ldz)
                     CALL DLACPY(uplo,n,n,B,Ldb,Bb,Ldb)
!
                     CALL DSYGV_2STAGE(ibtype,'N',uplo,n,Z,Ldz,Bb,Ldb,  &
     &                                 D2,Work,Nwork,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'DSYGV_2STAGE(V,'//    &
     &                         uplo//')' , iinfo , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(ntest) = ulpinv
                           GOTO 2
                        ENDIF
                     ENDIF
!
!                 Do Test
!
!                  CALL dsgt01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z,
!     $                         LDZ, D, WORK, RESULT( NTEST ) )
!
!                 Do Tests | D1 - D2 | / ( |D1| ulp )
!                 D1 computed using the standard 1-stage reduction as reference
!                 D2 computed using the 2-stage reduction
!
                     temp1 = ZERO
                     temp2 = ZERO
                     DO j = 1 , n
                        temp1 = MAX(temp1,ABS(D(j)),ABS(D2(j)))
                        temp2 = MAX(temp2,ABS(D(j)-D2(j)))
                     ENDDO
!
                     Result(ntest) = temp2/MAX(unfl,ulp*MAX(temp1,temp2)&
     &                               )
!
!                 Test DSYGVD
!
                     ntest = ntest + 1
!
                     CALL DLACPY(' ',n,n,A,Lda,Z,Ldz)
                     CALL DLACPY(uplo,n,n,B,Ldb,Bb,Ldb)
!
                     CALL DSYGVD(ibtype,'V',uplo,n,Z,Ldz,Bb,Ldb,D,Work, &
     &                           Nwork,Iwork,Liwork,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'DSYGVD(V,'//uplo//    &
     &                         ')' , iinfo , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(ntest) = ulpinv
                           GOTO 2
                        ENDIF
                     ENDIF
!
!                 Do Test
!
                     CALL dsgt01(ibtype,uplo,n,n,A,Lda,B,Ldb,Z,Ldz,D,   &
     &                           Work,Result(ntest))
!
!                 Test DSYGVX
!
                     ntest = ntest + 1
!
                     CALL DLACPY(' ',n,n,A,Lda,Ab,Lda)
                     CALL DLACPY(uplo,n,n,B,Ldb,Bb,Ldb)
!
                     CALL DSYGVX(ibtype,'V','A',uplo,n,Ab,Lda,Bb,Ldb,vl,&
     &                           vu,il,iu,abstol,m,D,Z,Ldz,Work,Nwork,  &
     &                           Iwork(n+1),Iwork,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'DSYGVX(V,A'//uplo//   &
     &                         ')' , iinfo , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(ntest) = ulpinv
                           GOTO 2
                        ENDIF
                     ENDIF
!
!                 Do Test
!
                     CALL dsgt01(ibtype,uplo,n,n,A,Lda,B,Ldb,Z,Ldz,D,   &
     &                           Work,Result(ntest))
!
                     ntest = ntest + 1
!
                     CALL DLACPY(' ',n,n,A,Lda,Ab,Lda)
                     CALL DLACPY(uplo,n,n,B,Ldb,Bb,Ldb)
!
!                 since we do not know the exact eigenvalues of this
!                 eigenpair, we just set VL and VU as constants.
!                 It is quite possible that there are no eigenvalues
!                 in this interval.
!
                     vl = ZERO
                     vu = anorm
                     CALL DSYGVX(ibtype,'V','V',uplo,n,Ab,Lda,Bb,Ldb,vl,&
     &                           vu,il,iu,abstol,m,D,Z,Ldz,Work,Nwork,  &
     &                           Iwork(n+1),Iwork,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'DSYGVX(V,V,'//uplo//  &
     &                         ')' , iinfo , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(ntest) = ulpinv
                           GOTO 2
                        ENDIF
                     ENDIF
!
!                 Do Test
!
                     CALL dsgt01(ibtype,uplo,n,m,A,Lda,B,Ldb,Z,Ldz,D,   &
     &                           Work,Result(ntest))
!
                     ntest = ntest + 1
!
                     CALL DLACPY(' ',n,n,A,Lda,Ab,Lda)
                     CALL DLACPY(uplo,n,n,B,Ldb,Bb,Ldb)
!
                     CALL DSYGVX(ibtype,'V','I',uplo,n,Ab,Lda,Bb,Ldb,vl,&
     &                           vu,il,iu,abstol,m,D,Z,Ldz,Work,Nwork,  &
     &                           Iwork(n+1),Iwork,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'DSYGVX(V,I,'//uplo//  &
     &                         ')' , iinfo , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(ntest) = ulpinv
                           GOTO 2
                        ENDIF
                     ENDIF
!
!                 Do Test
!
                     CALL dsgt01(ibtype,uplo,n,m,A,Lda,B,Ldb,Z,Ldz,D,   &
     &                           Work,Result(ntest))
!
!
!                 Test DSPGV
!
 2                   ntest = ntest + 1
!
!                 Copy the matrices into packed storage.
!
                     IF ( LSAME(uplo,'U') ) THEN
                        ij = 1
                        DO j = 1 , n
                           DO i = 1 , j
                              Ap(ij) = A(i,j)
                              Bp(ij) = B(i,j)
                              ij = ij + 1
                           ENDDO
                        ENDDO
                     ELSE
                        ij = 1
                        DO j = 1 , n
                           DO i = j , n
                              Ap(ij) = A(i,j)
                              Bp(ij) = B(i,j)
                              ij = ij + 1
                           ENDDO
                        ENDDO
                     ENDIF
!
                     CALL DSPGV(ibtype,'V',uplo,n,Ap,Bp,D,Z,Ldz,Work,   &
     &                          iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'DSPGV(V,'//uplo//')' ,&
     &                         iinfo , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(ntest) = ulpinv
                           GOTO 4
                        ENDIF
                     ENDIF
!
!                 Do Test
!
                     CALL dsgt01(ibtype,uplo,n,n,A,Lda,B,Ldb,Z,Ldz,D,   &
     &                           Work,Result(ntest))
!
!                 Test DSPGVD
!
                     ntest = ntest + 1
!
!                 Copy the matrices into packed storage.
!
                     IF ( LSAME(uplo,'U') ) THEN
                        ij = 1
                        DO j = 1 , n
                           DO i = 1 , j
                              Ap(ij) = A(i,j)
                              Bp(ij) = B(i,j)
                              ij = ij + 1
                           ENDDO
                        ENDDO
                     ELSE
                        ij = 1
                        DO j = 1 , n
                           DO i = j , n
                              Ap(ij) = A(i,j)
                              Bp(ij) = B(i,j)
                              ij = ij + 1
                           ENDDO
                        ENDDO
                     ENDIF
!
                     CALL DSPGVD(ibtype,'V',uplo,n,Ap,Bp,D,Z,Ldz,Work,  &
     &                           Nwork,Iwork,Liwork,iinfo)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'DSPGVD(V,'//uplo//    &
     &                         ')' , iinfo , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(ntest) = ulpinv
                           GOTO 4
                        ENDIF
                     ENDIF
!
!                 Do Test
!
                     CALL dsgt01(ibtype,uplo,n,n,A,Lda,B,Ldb,Z,Ldz,D,   &
     &                           Work,Result(ntest))
!
!                 Test DSPGVX
!
                     ntest = ntest + 1
!
!                 Copy the matrices into packed storage.
!
                     IF ( LSAME(uplo,'U') ) THEN
                        ij = 1
                        DO j = 1 , n
                           DO i = 1 , j
                              Ap(ij) = A(i,j)
                              Bp(ij) = B(i,j)
                              ij = ij + 1
                           ENDDO
                        ENDDO
                     ELSE
                        ij = 1
                        DO j = 1 , n
                           DO i = j , n
                              Ap(ij) = A(i,j)
                              Bp(ij) = B(i,j)
                              ij = ij + 1
                           ENDDO
                        ENDDO
                     ENDIF
!
                     CALL DSPGVX(ibtype,'V','A',uplo,n,Ap,Bp,vl,vu,il,  &
     &                           iu,abstol,m,D,Z,Ldz,Work,Iwork(n+1),   &
     &                           Iwork,Info)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'DSPGVX(V,A'//uplo//   &
     &                         ')' , iinfo , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(ntest) = ulpinv
                           GOTO 4
                        ENDIF
                     ENDIF
!
!                 Do Test
!
                     CALL dsgt01(ibtype,uplo,n,m,A,Lda,B,Ldb,Z,Ldz,D,   &
     &                           Work,Result(ntest))
!
                     ntest = ntest + 1
!
!                 Copy the matrices into packed storage.
!
                     IF ( LSAME(uplo,'U') ) THEN
                        ij = 1
                        DO j = 1 , n
                           DO i = 1 , j
                              Ap(ij) = A(i,j)
                              Bp(ij) = B(i,j)
                              ij = ij + 1
                           ENDDO
                        ENDDO
                     ELSE
                        ij = 1
                        DO j = 1 , n
                           DO i = j , n
                              Ap(ij) = A(i,j)
                              Bp(ij) = B(i,j)
                              ij = ij + 1
                           ENDDO
                        ENDDO
                     ENDIF
!
                     vl = ZERO
                     vu = anorm
                     CALL DSPGVX(ibtype,'V','V',uplo,n,Ap,Bp,vl,vu,il,  &
     &                           iu,abstol,m,D,Z,Ldz,Work,Iwork(n+1),   &
     &                           Iwork,Info)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'DSPGVX(V,V'//uplo//   &
     &                         ')' , iinfo , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(ntest) = ulpinv
                           GOTO 4
                        ENDIF
                     ENDIF
!
!                 Do Test
!
                     CALL dsgt01(ibtype,uplo,n,m,A,Lda,B,Ldb,Z,Ldz,D,   &
     &                           Work,Result(ntest))
!
                     ntest = ntest + 1
!
!                 Copy the matrices into packed storage.
!
                     IF ( LSAME(uplo,'U') ) THEN
                        ij = 1
                        DO j = 1 , n
                           DO i = 1 , j
                              Ap(ij) = A(i,j)
                              Bp(ij) = B(i,j)
                              ij = ij + 1
                           ENDDO
                        ENDDO
                     ELSE
                        ij = 1
                        DO j = 1 , n
                           DO i = j , n
                              Ap(ij) = A(i,j)
                              Bp(ij) = B(i,j)
                              ij = ij + 1
                           ENDDO
                        ENDDO
                     ENDIF
!
                     CALL DSPGVX(ibtype,'V','I',uplo,n,Ap,Bp,vl,vu,il,  &
     &                           iu,abstol,m,D,Z,Ldz,Work,Iwork(n+1),   &
     &                           Iwork,Info)
                     IF ( iinfo/=0 ) THEN
                        WRITE (Nounit,FMT=99001) 'DSPGVX(V,I'//uplo//   &
     &                         ')' , iinfo , n , jtype , ioldsd
                        Info = ABS(iinfo)
                        IF ( iinfo<0 ) THEN
                           RETURN
                        ELSE
                           Result(ntest) = ulpinv
                           GOTO 4
                        ENDIF
                     ENDIF
!
!                 Do Test
!
                     CALL dsgt01(ibtype,uplo,n,m,A,Lda,B,Ldb,Z,Ldz,D,   &
     &                           Work,Result(ntest))
!
!
 4                   IF ( ibtype==1 ) THEN
!
!                    TEST DSBGV
!
                        ntest = ntest + 1
!
!                    Copy the matrices into band storage.
!
                        IF ( LSAME(uplo,'U') ) THEN
                           DO j = 1 , n
                              DO i = MAX(1,j-ka) , j
                                 Ab(ka+1+i-j,j) = A(i,j)
                              ENDDO
                              DO i = MAX(1,j-kb) , j
                                 Bb(kb+1+i-j,j) = B(i,j)
                              ENDDO
                           ENDDO
                        ELSE
                           DO j = 1 , n
                              DO i = j , MIN(n,j+ka)
                                 Ab(1+i-j,j) = A(i,j)
                              ENDDO
                              DO i = j , MIN(n,j+kb)
                                 Bb(1+i-j,j) = B(i,j)
                              ENDDO
                           ENDDO
                        ENDIF
!
                        CALL DSBGV('V',uplo,n,ka,kb,Ab,Lda,Bb,Ldb,D,Z,  &
     &                             Ldz,Work,iinfo)
                        IF ( iinfo/=0 ) THEN
                           WRITE (Nounit,FMT=99001) 'DSBGV(V,'//uplo//  &
     &                            ')' , iinfo , n , jtype , ioldsd
                           Info = ABS(iinfo)
                           IF ( iinfo<0 ) THEN
                              RETURN
                           ELSE
                              Result(ntest) = ulpinv
                              CYCLE
                           ENDIF
                        ENDIF
!
!                    Do Test
!
                        CALL dsgt01(ibtype,uplo,n,n,A,Lda,B,Ldb,Z,Ldz,D,&
     &                              Work,Result(ntest))
!
!                    TEST DSBGVD
!
                        ntest = ntest + 1
!
!                    Copy the matrices into band storage.
!
                        IF ( LSAME(uplo,'U') ) THEN
                           DO j = 1 , n
                              DO i = MAX(1,j-ka) , j
                                 Ab(ka+1+i-j,j) = A(i,j)
                              ENDDO
                              DO i = MAX(1,j-kb) , j
                                 Bb(kb+1+i-j,j) = B(i,j)
                              ENDDO
                           ENDDO
                        ELSE
                           DO j = 1 , n
                              DO i = j , MIN(n,j+ka)
                                 Ab(1+i-j,j) = A(i,j)
                              ENDDO
                              DO i = j , MIN(n,j+kb)
                                 Bb(1+i-j,j) = B(i,j)
                              ENDDO
                           ENDDO
                        ENDIF
!
                        CALL DSBGVD('V',uplo,n,ka,kb,Ab,Lda,Bb,Ldb,D,Z, &
     &                              Ldz,Work,Nwork,Iwork,Liwork,iinfo)
                        IF ( iinfo/=0 ) THEN
                           WRITE (Nounit,FMT=99001) 'DSBGVD(V,'//uplo// &
     &                            ')' , iinfo , n , jtype , ioldsd
                           Info = ABS(iinfo)
                           IF ( iinfo<0 ) THEN
                              RETURN
                           ELSE
                              Result(ntest) = ulpinv
                              CYCLE
                           ENDIF
                        ENDIF
!
!                    Do Test
!
                        CALL dsgt01(ibtype,uplo,n,n,A,Lda,B,Ldb,Z,Ldz,D,&
     &                              Work,Result(ntest))
!
!                    Test DSBGVX
!
                        ntest = ntest + 1
!
!                    Copy the matrices into band storage.
!
                        IF ( LSAME(uplo,'U') ) THEN
                           DO j = 1 , n
                              DO i = MAX(1,j-ka) , j
                                 Ab(ka+1+i-j,j) = A(i,j)
                              ENDDO
                              DO i = MAX(1,j-kb) , j
                                 Bb(kb+1+i-j,j) = B(i,j)
                              ENDDO
                           ENDDO
                        ELSE
                           DO j = 1 , n
                              DO i = j , MIN(n,j+ka)
                                 Ab(1+i-j,j) = A(i,j)
                              ENDDO
                              DO i = j , MIN(n,j+kb)
                                 Bb(1+i-j,j) = B(i,j)
                              ENDDO
                           ENDDO
                        ENDIF
!
                        CALL DSBGVX('V','A',uplo,n,ka,kb,Ab,Lda,Bb,Ldb, &
     &                              Bp,MAX(1,n),vl,vu,il,iu,abstol,m,D, &
     &                              Z,Ldz,Work,Iwork(n+1),Iwork,iinfo)
                        IF ( iinfo/=0 ) THEN
                           WRITE (Nounit,FMT=99001) 'DSBGVX(V,A'//uplo//&
     &                            ')' , iinfo , n , jtype , ioldsd
                           Info = ABS(iinfo)
                           IF ( iinfo<0 ) THEN
                              RETURN
                           ELSE
                              Result(ntest) = ulpinv
                              CYCLE
                           ENDIF
                        ENDIF
!
!                    Do Test
!
                        CALL dsgt01(ibtype,uplo,n,m,A,Lda,B,Ldb,Z,Ldz,D,&
     &                              Work,Result(ntest))
!
!
                        ntest = ntest + 1
!
!                    Copy the matrices into band storage.
!
                        IF ( LSAME(uplo,'U') ) THEN
                           DO j = 1 , n
                              DO i = MAX(1,j-ka) , j
                                 Ab(ka+1+i-j,j) = A(i,j)
                              ENDDO
                              DO i = MAX(1,j-kb) , j
                                 Bb(kb+1+i-j,j) = B(i,j)
                              ENDDO
                           ENDDO
                        ELSE
                           DO j = 1 , n
                              DO i = j , MIN(n,j+ka)
                                 Ab(1+i-j,j) = A(i,j)
                              ENDDO
                              DO i = j , MIN(n,j+kb)
                                 Bb(1+i-j,j) = B(i,j)
                              ENDDO
                           ENDDO
                        ENDIF
!
                        vl = ZERO
                        vu = anorm
                        CALL DSBGVX('V','V',uplo,n,ka,kb,Ab,Lda,Bb,Ldb, &
     &                              Bp,MAX(1,n),vl,vu,il,iu,abstol,m,D, &
     &                              Z,Ldz,Work,Iwork(n+1),Iwork,iinfo)
                        IF ( iinfo/=0 ) THEN
                           WRITE (Nounit,FMT=99001) 'DSBGVX(V,V'//uplo//&
     &                            ')' , iinfo , n , jtype , ioldsd
                           Info = ABS(iinfo)
                           IF ( iinfo<0 ) THEN
                              RETURN
                           ELSE
                              Result(ntest) = ulpinv
                              CYCLE
                           ENDIF
                        ENDIF
!
!                    Do Test
!
                        CALL dsgt01(ibtype,uplo,n,m,A,Lda,B,Ldb,Z,Ldz,D,&
     &                              Work,Result(ntest))
!
                        ntest = ntest + 1
!
!                    Copy the matrices into band storage.
!
                        IF ( LSAME(uplo,'U') ) THEN
                           DO j = 1 , n
                              DO i = MAX(1,j-ka) , j
                                 Ab(ka+1+i-j,j) = A(i,j)
                              ENDDO
                              DO i = MAX(1,j-kb) , j
                                 Bb(kb+1+i-j,j) = B(i,j)
                              ENDDO
                           ENDDO
                        ELSE
                           DO j = 1 , n
                              DO i = j , MIN(n,j+ka)
                                 Ab(1+i-j,j) = A(i,j)
                              ENDDO
                              DO i = j , MIN(n,j+kb)
                                 Bb(1+i-j,j) = B(i,j)
                              ENDDO
                           ENDDO
                        ENDIF
!
                        CALL DSBGVX('V','I',uplo,n,ka,kb,Ab,Lda,Bb,Ldb, &
     &                              Bp,MAX(1,n),vl,vu,il,iu,abstol,m,D, &
     &                              Z,Ldz,Work,Iwork(n+1),Iwork,iinfo)
                        IF ( iinfo/=0 ) THEN
                           WRITE (Nounit,FMT=99001) 'DSBGVX(V,I'//uplo//&
     &                            ')' , iinfo , n , jtype , ioldsd
                           Info = ABS(iinfo)
                           IF ( iinfo<0 ) THEN
                              RETURN
                           ELSE
                              Result(ntest) = ulpinv
                              CYCLE
                           ENDIF
                        ENDIF
!
!                    Do Test
!
                        CALL dsgt01(ibtype,uplo,n,m,A,Lda,B,Ldb,Z,Ldz,D,&
     &                              Work,Result(ntest))
!
                     ENDIF
!
                  ENDDO
               ENDDO
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
               ntestt = ntestt + ntest
               CALL dlafts('DSG',n,n,jtype,ntest,Result,ioldsd,Thresh,  &
     &                     Nounit,nerrs)
            ENDIF
         ENDDO
      ENDDO
!
!     Summary
!
      CALL dlasum('DSG',Nounit,nerrs,ntestt)
!
      RETURN
!
!     End of ddrvsg2stg
!
99001 FORMAT (' ddrvsg2stg: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,   &
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
      END SUBROUTINE ddrvsg2stg
