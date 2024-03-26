!*==cdrges.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b cdrges
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cdrges( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          NOUNIT, A, LDA, B, S, T, Q, LDQ, Z, ALPHA,
!                          BETA, WORK, LWORK, RWORK, RESULT, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDQ, LWORK, NOUNIT, NSIZES, NTYPES
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * ), DOTYPE( * )
!       INTEGER            ISEED( 4 ), NN( * )
!       REAL               RESULT( 13 ), RWORK( * )
!       COMPLEX            A( LDA, * ), ALPHA( * ), B( LDA, * ),
!      $                   BETA( * ), Q( LDQ, * ), S( LDA, * ),
!      $                   T( LDA, * ), WORK( * ), Z( LDQ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cdrges checks the nonsymmetric generalized eigenvalue (Schur form)
!> problem driver CGGES.
!>
!> CGGES factors A and B as Q*S*Z'  and Q*T*Z' , where ' means conjugate
!> transpose, S and T are  upper triangular (i.e., in generalized Schur
!> form), and Q and Z are unitary. It also computes the generalized
!> eigenvalues (alpha(j),beta(j)), j=1,...,n.  Thus,
!> w(j) = alpha(j)/beta(j) is a root of the characteristic equation
!>
!>                 det( A - w(j) B ) = 0
!>
!> Optionally it also reorder the eigenvalues so that a selected
!> cluster of eigenvalues appears in the leading diagonal block of the
!> Schur forms.
!>
!> When cdrges is called, a number of matrix "sizes" ("N's") and a
!> number of matrix "TYPES" are specified.  For each size ("N")
!> and each TYPE of matrix, a pair of matrices (A, B) will be generated
!> and used for testing. For each matrix pair, the following 13 tests
!> will be performed and compared with the threshold THRESH except
!> the tests (5), (11) and (13).
!>
!>
!> (1)   | A - Q S Z' | / ( |A| n ulp ) (no sorting of eigenvalues)
!>
!>
!> (2)   | B - Q T Z' | / ( |B| n ulp ) (no sorting of eigenvalues)
!>
!>
!> (3)   | I - QQ' | / ( n ulp ) (no sorting of eigenvalues)
!>
!>
!> (4)   | I - ZZ' | / ( n ulp ) (no sorting of eigenvalues)
!>
!> (5)   if A is in Schur form (i.e. triangular form) (no sorting of
!>       eigenvalues)
!>
!> (6)   if eigenvalues = diagonal elements of the Schur form (S, T),
!>       i.e., test the maximum over j of D(j)  where:
!>
!>                     |alpha(j) - S(j,j)|        |beta(j) - T(j,j)|
!>           D(j) = ------------------------ + -----------------------
!>                  max(|alpha(j)|,|S(j,j)|)   max(|beta(j)|,|T(j,j)|)
!>
!>       (no sorting of eigenvalues)
!>
!> (7)   | (A,B) - Q (S,T) Z' | / ( |(A,B)| n ulp )
!>       (with sorting of eigenvalues).
!>
!> (8)   | I - QQ' | / ( n ulp ) (with sorting of eigenvalues).
!>
!> (9)   | I - ZZ' | / ( n ulp ) (with sorting of eigenvalues).
!>
!> (10)  if A is in Schur form (i.e. quasi-triangular form)
!>       (with sorting of eigenvalues).
!>
!> (11)  if eigenvalues = diagonal elements of the Schur form (S, T),
!>       i.e. test the maximum over j of D(j)  where:
!>
!>                     |alpha(j) - S(j,j)|        |beta(j) - T(j,j)|
!>           D(j) = ------------------------ + -----------------------
!>                  max(|alpha(j)|,|S(j,j)|)   max(|beta(j)|,|T(j,j)|)
!>
!>       (with sorting of eigenvalues).
!>
!> (12)  if sorting worked and SDIM is the number of eigenvalues
!>       which were CELECTed.
!>
!> Test Matrices
!> =============
!>
!> The sizes of the test matrices are specified by an array
!> NN(1:NSIZES); the value of each element NN(j) specifies one size.
!> The "types" are specified by a logical array DOTYPE( 1:NTYPES ); if
!> DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
!> Currently, the list of possible types is:
!>
!> (1)  ( 0, 0 )         (a pair of zero matrices)
!>
!> (2)  ( I, 0 )         (an identity and a zero matrix)
!>
!> (3)  ( 0, I )         (an identity and a zero matrix)
!>
!> (4)  ( I, I )         (a pair of identity matrices)
!>
!>         t   t
!> (5)  ( J , J  )       (a pair of transposed Jordan blocks)
!>
!>                                     t                ( I   0  )
!> (6)  ( X, Y )         where  X = ( J   0  )  and Y = (      t )
!>                                  ( 0   I  )          ( 0   J  )
!>                       and I is a k x k identity and J a (k+1)x(k+1)
!>                       Jordan block; k=(N-1)/2
!>
!> (7)  ( D, I )         where D is diag( 0, 1,..., N-1 ) (a diagonal
!>                       matrix with those diagonal entries.)
!> (8)  ( I, D )
!>
!> (9)  ( big*D, small*I ) where "big" is near overflow and small=1/big
!>
!> (10) ( small*D, big*I )
!>
!> (11) ( big*I, small*D )
!>
!> (12) ( small*I, big*D )
!>
!> (13) ( big*D, big*I )
!>
!> (14) ( small*D, small*I )
!>
!> (15) ( D1, D2 )        where D1 is diag( 0, 0, 1, ..., N-3, 0 ) and
!>                        D2 is diag( 0, N-3, N-4,..., 1, 0, 0 )
!>           t   t
!> (16) Q ( J , J ) Z     where Q and Z are random orthogonal matrices.
!>
!> (17) Q ( T1, T2 ) Z    where T1 and T2 are upper triangular matrices
!>                        with random O(1) entries above the diagonal
!>                        and diagonal entries diag(T1) =
!>                        ( 0, 0, 1, ..., N-3, 0 ) and diag(T2) =
!>                        ( 0, N-3, N-4,..., 1, 0, 0 )
!>
!> (18) Q ( T1, T2 ) Z    diag(T1) = ( 0, 0, 1, 1, s, ..., s, 0 )
!>                        diag(T2) = ( 0, 1, 0, 1,..., 1, 0 )
!>                        s = machine precision.
!>
!> (19) Q ( T1, T2 ) Z    diag(T1)=( 0,0,1,1, 1-d, ..., 1-(N-5)*d=s, 0 )
!>                        diag(T2) = ( 0, 1, 0, 1, ..., 1, 0 )
!>
!>                                                        N-5
!> (20) Q ( T1, T2 ) Z    diag(T1)=( 0, 0, 1, 1, a, ..., a   =s, 0 )
!>                        diag(T2) = ( 0, 1, 0, 1, ..., 1, 0, 0 )
!>
!> (21) Q ( T1, T2 ) Z    diag(T1)=( 0, 0, 1, r1, r2, ..., r(N-4), 0 )
!>                        diag(T2) = ( 0, 1, 0, 1, ..., 1, 0, 0 )
!>                        where r1,..., r(N-4) are random.
!>
!> (22) Q ( big*T1, small*T2 ) Z    diag(T1) = ( 0, 0, 1, ..., N-3, 0 )
!>                                  diag(T2) = ( 0, 1, ..., 1, 0, 0 )
!>
!> (23) Q ( small*T1, big*T2 ) Z    diag(T1) = ( 0, 0, 1, ..., N-3, 0 )
!>                                  diag(T2) = ( 0, 1, ..., 1, 0, 0 )
!>
!> (24) Q ( small*T1, small*T2 ) Z  diag(T1) = ( 0, 0, 1, ..., N-3, 0 )
!>                                  diag(T2) = ( 0, 1, ..., 1, 0, 0 )
!>
!> (25) Q ( big*T1, big*T2 ) Z      diag(T1) = ( 0, 0, 1, ..., N-3, 0 )
!>                                  diag(T2) = ( 0, 1, ..., 1, 0, 0 )
!>
!> (26) Q ( T1, T2 ) Z     where T1 and T2 are random upper-triangular
!>                         matrices.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          sdrges does nothing.  NSIZES >= 0.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER array, dimension (NSIZES)
!>          An array containing the sizes to be used for the matrices.
!>          Zero values will be skipped.  NN >= 0.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, sdrges
!>          does nothing.  It must be at least zero.  If it is MAXTYP+1
!>          and NSIZES is 1, then an additional type, MAXTYP+1 is
!>          defined, which is to use whatever matrix is in A on input.
!>          This is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
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
!>          MAXTYP will not be generated. If NTYPES is larger
!>          than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
!>          will be ignored.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator. The array elements should be between 0 and 4095;
!>          if not they will be reduced mod 4096. Also, ISEED(4) must
!>          be odd.  The random number generator uses a linear
!>          congruential sequence limited to small integers, and so
!>          should produce machine independent random numbers. The
!>          values of ISEED are changed on exit, and can be used in the
!>          next call to sdrges to continue the same random number
!>          sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          A test will count as "failed" if the "error", computed as
!>          described above, exceeds THRESH.  Note that the error is
!>          scaled to be O(1), so THRESH should be a reasonably small
!>          multiple of 1, e.g., 10 or 100.  In particular, it should
!>          not depend on the precision (single vs. double) or the size
!>          of the matrix.  THRESH >= 0.
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
!>          A is COMPLEX array, dimension(LDA, max(NN))
!>          Used to hold the original A matrix.  Used as input only
!>          if NTYPES=MAXTYP+1, DOTYPE(1:MAXTYP)=.FALSE., and
!>          DOTYPE(MAXTYP+1)=.TRUE.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A, B, S, and T.
!>          It must be at least 1 and at least max( NN ).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension(LDA, max(NN))
!>          Used to hold the original B matrix.  Used as input only
!>          if NTYPES=MAXTYP+1, DOTYPE(1:MAXTYP)=.FALSE., and
!>          DOTYPE(MAXTYP+1)=.TRUE.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is COMPLEX array, dimension (LDA, max(NN))
!>          The Schur form matrix computed from A by CGGES.  On exit, S
!>          contains the Schur form matrix corresponding to the matrix
!>          in A.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDA, max(NN))
!>          The upper triangular matrix computed from B by CGGES.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDQ, max(NN))
!>          The (left) orthogonal matrix computed by CGGES.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of Q and Z. It must
!>          be at least 1 and at least max( NN ).
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension( LDQ, max(NN) )
!>          The (right) orthogonal matrix computed by CGGES.
!> \endverbatim
!>
!> \param[out] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is COMPLEX array, dimension (max(NN))
!>
!>          The generalized eigenvalues of (A,B) computed by CGGES.
!>          ALPHA(k) / BETA(k) is the k-th generalized eigenvalue of A
!>          and B.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= 3*N*N.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension ( 8*N )
!>          Real workspace.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (15)
!>          The values computed by the tests described above.
!>          The values are currently limited to 1/ulp, to avoid overflow.
!> \endverbatim
!>
!> \param[out] BWORK
!> \verbatim
!>          BWORK is LOGICAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  A routine returned an error code.  INFO is the
!>                absolute value of the INFO value returned.
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE cdrges(Nsizes,Nn,Ntypes,Dotype,Iseed,Thresh,Nounit,A,  &
     &                  Lda,B,S,T,Q,Ldq,Z,Alpha,Beta,Work,Lwork,Rwork,  &
     &                  Result,Bwork,Info)
      IMPLICIT NONE
!*--CDRGES385
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldq , Lwork , Nounit , Nsizes , Ntypes
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Bwork(*) , Dotype(*)
      INTEGER Iseed(4) , Nn(*)
      REAL Result(13) , Rwork(*)
      COMPLEX A(Lda,*) , Alpha(*) , B(Lda,*) , Beta(*) , Q(Ldq,*) ,     &
     &        S(Lda,*) , T(Lda,*) , Work(*) , Z(Ldq,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
      INTEGER MAXTYP
      PARAMETER (MAXTYP=26)
!     ..
!     .. Local Scalars ..
      LOGICAL badnn , ilabad
      CHARACTER sort
      INTEGER i , iadd , iinfo , in , isort , j , jc , jr , jsize ,     &
     &        jtype , knteig , maxwrk , minwrk , mtypes , n , n1 , nb , &
     &        nerrs , nmats , nmax , ntest , ntestt , rsub , sdim
      REAL safmax , safmin , temp1 , temp2 , ulp , ulpinv
      COMPLEX ctemp , x
!     ..
!     .. Local Arrays ..
      LOGICAL lasign(MAXTYP) , lbsign(MAXTYP)
      INTEGER ioldsd(4) , kadd(6) , kamagn(MAXTYP) , katype(MAXTYP) ,   &
     &        kazero(MAXTYP) , kbmagn(MAXTYP) , kbtype(MAXTYP) ,        &
     &        kbzero(MAXTYP) , kclass(MAXTYP) , ktrian(MAXTYP) , kz1(6) &
     &        , kz2(6)
      REAL rmagn(0:3)
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH 
!     ..
!     .. External Subroutines ..
      EXTERNAL CGGES , CLACPY , CLARFG ,     CLASET , CUNM2R , SLABAD
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , CONJG , MAX , MIN , REAL , SIGN
!     ..
!     .. Statement Functions ..
      REAL ABS1
!     ..
!     .. Statement Function definitions ..
      ABS1(x) = ABS(REAL(x)) + ABS(AIMAG(x))
!     ..
!     .. Data statements ..
      DATA kclass/15*1 , 10*2 , 1*3/
      DATA kz1/0 , 1 , 2 , 1 , 3 , 3/
      DATA kz2/0 , 0 , 1 , 2 , 1 , 1/
      DATA kadd/0 , 0 , 0 , 0 , 3 , 2/
      DATA katype/0 , 1 , 0 , 1 , 2 , 3 , 4 , 1 , 4 , 4 , 1 , 1 , 4 ,   &
     &     4 , 4 , 2 , 4 , 5 , 8 , 7 , 9 , 4*4 , 0/
      DATA kbtype/0 , 0 , 1 , 1 , 2 , -3 , 1 , 4 , 1 , 1 , 4 , 4 , 1 ,  &
     &     1 , -4 , 2 , -4 , 8*8 , 0/
      DATA kazero/6*1 , 2 , 1 , 2*2 , 2*1 , 2*2 , 3 , 1 , 3 , 4*5 ,     &
     &     4*3 , 1/
      DATA kbzero/6*1 , 1 , 2 , 2*1 , 2*2 , 2*1 , 4 , 1 , 4 , 4*6 ,     &
     &     4*4 , 1/
      DATA kamagn/8*1 , 2 , 3 , 2 , 3 , 2 , 3 , 7*1 , 2 , 3 , 3 , 2 , 1/
      DATA kbmagn/8*1 , 3 , 2 , 3 , 2 , 2 , 3 , 7*1 , 3 , 2 , 3 , 2 , 1/
      DATA ktrian/16*0 , 10*1/
      DATA lasign/6*.FALSE. , .TRUE. , .FALSE. , 2*.TRUE. , 2*.FALSE. , &
     &     3*.TRUE. , .FALSE. , .TRUE. , 3*.FALSE. , 5*.TRUE. , .FALSE./
      DATA lbsign/7*.FALSE. , .TRUE. , 2*.FALSE. , 2*.TRUE. ,           &
     &     2*.FALSE. , .TRUE. , .FALSE. , .TRUE. , 9*.FALSE./
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
      Info = 0
!
      badnn = .FALSE.
      nmax = 1
      DO j = 1 , Nsizes
         nmax = MAX(nmax,Nn(j))
         IF ( Nn(j)<0 ) badnn = .TRUE.
      ENDDO
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
      ELSEIF ( Ldq<=1 .OR. Ldq<nmax ) THEN
         Info = -14
      ENDIF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ilaenv.
!
      minwrk = 1
      IF ( Info==0 .AND. Lwork>=1 ) THEN
         minwrk = 3*nmax*nmax
         nb = MAX(1,ilaenv(1,'CGEQRF',' ',nmax,nmax,-1,-1),             &
     &        ilaenv(1,'CUNMQR','LC',nmax,nmax,nmax,-1),                &
     &        ilaenv(1,'CUNGQR',' ',nmax,nmax,nmax,-1))
         maxwrk = MAX(nmax+nmax*nb,3*nmax*nmax)
         Work(1) = maxwrk
      ENDIF
!
      IF ( Lwork<minwrk ) Info = -19
!
      IF ( Info/=0 ) THEN
         CALL xerbla('cdrges',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( Nsizes==0 .OR. Ntypes==0 ) RETURN
!
      ulp = SLAMCH('Precision')
      safmin = SLAMCH('Safe minimum')
      safmin = safmin/ulp
      safmax = ONE/safmin
      CALL SLABAD(safmin,safmax)
      ulpinv = ONE/ulp
!
!     The values RMAGN(2:3) depend on N, see below.
!
      rmagn(0) = ZERO
      rmagn(1) = ONE
!
!     Loop over matrix sizes
!
      ntestt = 0
      nerrs = 0
      nmats = 0
!
      DO jsize = 1 , Nsizes
         n = Nn(jsize)
         n1 = MAX(1,n)
         rmagn(2) = safmax*ulp/REAL(n1)
         rmagn(3) = safmin*ulpinv*REAL(n1)
!
         IF ( Nsizes/=1 ) THEN
            mtypes = MIN(MAXTYP,Ntypes)
         ELSE
            mtypes = MIN(MAXTYP+1,Ntypes)
         ENDIF
!
!        Loop over matrix types
!
         DO jtype = 1 , mtypes
            IF ( Dotype(jtype) ) THEN
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
               DO j = 1 , 13
                  Result(j) = ZERO
               ENDDO
!
!           Generate test matrices A and B
!
!           Description of control parameters:
!
!           KCLASS: =1 means w/o rotation, =2 means w/ rotation,
!                   =3 means random.
!           KATYPE: the "type" to be passed to clatm4 for computing A.
!           KAZERO: the pattern of zeros on the diagonal for A:
!                   =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
!                   =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
!                   =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
!                   non-zero entries.)
!           KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
!                   =2: large, =3: small.
!           LASIGN: .TRUE. if the diagonal elements of A are to be
!                   multiplied by a random magnitude 1 number.
!           KBTYPE, KBZERO, KBMAGN, LBSIGN: the same, but for B.
!           KTRIAN: =0: don't fill in the upper triangle, =1: do.
!           KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
!           RMAGN: used to implement KAMAGN and KBMAGN.
!
               IF ( mtypes<=MAXTYP ) THEN
                  iinfo = 0
                  IF ( kclass(jtype)<3 ) THEN
!
!              Generate A (w/o rotation)
!
                     IF ( ABS(katype(jtype))==3 ) THEN
                        in = 2*((n-1)/2) + 1
                        IF ( in/=n ) CALL CLASET('Full',n,n,CZERO,CZERO,&
     &                       A,Lda)
                     ELSE
                        in = n
                     ENDIF
                     CALL clatm4(katype(jtype),in,kz1(kazero(jtype)),   &
     &                           kz2(kazero(jtype)),lasign(jtype),      &
     &                           rmagn(kamagn(jtype)),ulp,              &
     &                           rmagn(ktrian(jtype)*kamagn(jtype)),2,  &
     &                           Iseed,A,Lda)
                     iadd = kadd(kazero(jtype))
                     IF ( iadd>0 .AND. iadd<=n ) A(iadd,iadd)           &
     &                    = rmagn(kamagn(jtype))
!
!              Generate B (w/o rotation)
!
                     IF ( ABS(kbtype(jtype))==3 ) THEN
                        in = 2*((n-1)/2) + 1
                        IF ( in/=n ) CALL CLASET('Full',n,n,CZERO,CZERO,&
     &                       B,Lda)
                     ELSE
                        in = n
                     ENDIF
                     CALL clatm4(kbtype(jtype),in,kz1(kbzero(jtype)),   &
     &                           kz2(kbzero(jtype)),lbsign(jtype),      &
     &                           rmagn(kbmagn(jtype)),ONE,              &
     &                           rmagn(ktrian(jtype)*kbmagn(jtype)),2,  &
     &                           Iseed,B,Lda)
                     iadd = kadd(kbzero(jtype))
                     IF ( iadd/=0 .AND. iadd<=n ) B(iadd,iadd)          &
     &                    = rmagn(kbmagn(jtype))
!
                     IF ( kclass(jtype)==2 .AND. n>0 ) THEN
!
!                 Include rotations
!
!                 Generate Q, Z as Householder transformations times
!                 a diagonal matrix.
!
                        DO jc = 1 , n - 1
                           DO jr = jc , n
                              Q(jr,jc) = clarnd(3,Iseed)
                              Z(jr,jc) = clarnd(3,Iseed)
                           ENDDO
                           CALL CLARFG(n+1-jc,Q(jc,jc),Q(jc+1,jc),1,    &
     &                                 Work(jc))
                           Work(2*n+jc) = SIGN(ONE,REAL(Q(jc,jc)))
                           Q(jc,jc) = CONE
                           CALL CLARFG(n+1-jc,Z(jc,jc),Z(jc+1,jc),1,    &
     &                                 Work(n+jc))
                           Work(3*n+jc) = SIGN(ONE,REAL(Z(jc,jc)))
                           Z(jc,jc) = CONE
                        ENDDO
                        ctemp = clarnd(3,Iseed)
                        Q(n,n) = CONE
                        Work(n) = CZERO
                        Work(3*n) = ctemp/ABS(ctemp)
                        ctemp = clarnd(3,Iseed)
                        Z(n,n) = CONE
                        Work(2*n) = CZERO
                        Work(4*n) = ctemp/ABS(ctemp)
!
!                 Apply the diagonal matrices
!
                        DO jc = 1 , n
                           DO jr = 1 , n
                              A(jr,jc) = Work(2*n+jr)                   &
     &                           *CONJG(Work(3*n+jc))*A(jr,jc)
                              B(jr,jc) = Work(2*n+jr)                   &
     &                           *CONJG(Work(3*n+jc))*B(jr,jc)
                           ENDDO
                        ENDDO
                        CALL CUNM2R('L','N',n,n,n-1,Q,Ldq,Work,A,Lda,   &
     &                              Work(2*n+1),iinfo)
                        IF ( iinfo==0 ) THEN
                           CALL CUNM2R('R','C',n,n,n-1,Z,Ldq,Work(n+1), &
     &                                 A,Lda,Work(2*n+1),iinfo)
                           IF ( iinfo==0 ) THEN
                              CALL CUNM2R('L','N',n,n,n-1,Q,Ldq,Work,B, &
     &                           Lda,Work(2*n+1),iinfo)
                              IF ( iinfo==0 ) THEN
                                 CALL CUNM2R('R','C',n,n,n-1,Z,Ldq,     &
     &                              Work(n+1),B,Lda,Work(2*n+1),iinfo)
                                 IF ( iinfo/=0 ) THEN
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ELSE
!
!              Random matrices
!
                     DO jc = 1 , n
                        DO jr = 1 , n
                           A(jr,jc) = rmagn(kamagn(jtype))              &
     &                                *clarnd(4,Iseed)
                           B(jr,jc) = rmagn(kbmagn(jtype))              &
     &                                *clarnd(4,Iseed)
                        ENDDO
                     ENDDO
                  ENDIF
!
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
               DO i = 1 , 13
                  Result(i) = -ONE
               ENDDO
!
!           Test with and without sorting of eigenvalues
!
               DO isort = 0 , 1
                  IF ( isort==0 ) THEN
                     sort = 'N'
                     rsub = 0
                  ELSE
                     sort = 'S'
                     rsub = 5
                  ENDIF
!
!              Call CGGES to compute H, T, Q, Z, alpha, and beta.
!
                  CALL CLACPY('Full',n,n,A,Lda,S,Lda)
                  CALL CLACPY('Full',n,n,B,Lda,T,Lda)
                  ntest = 1 + rsub + isort
                  Result(1+rsub+isort) = ulpinv
                  CALL CGGES('V','V',sort,clctes,n,S,Lda,T,Lda,sdim,    &
     &                       Alpha,Beta,Q,Ldq,Z,Ldq,Work,Lwork,Rwork,   &
     &                       Bwork,iinfo)
                  IF ( iinfo/=0 .AND. iinfo/=n+2 ) THEN
                     Result(1+rsub+isort) = ulpinv
                     WRITE (Nounit,FMT=99001) 'CGGES' , iinfo , n ,     &
     &                      jtype , ioldsd
                     Info = ABS(iinfo)
                     EXIT
                  ENDIF
!
                  ntest = 4 + rsub
!
!              Do tests 1--4 (or tests 7--9 when reordering )
!
                  IF ( isort==0 ) THEN
                     CALL cget51(1,n,A,Lda,S,Lda,Q,Ldq,Z,Ldq,Work,Rwork,&
     &                           Result(1))
                     CALL cget51(1,n,B,Lda,T,Lda,Q,Ldq,Z,Ldq,Work,Rwork,&
     &                           Result(2))
                  ELSE
                     CALL cget54(n,A,Lda,B,Lda,S,Lda,T,Lda,Q,Ldq,Z,Ldq, &
     &                           Work,Result(2+rsub))
                  ENDIF
!
                  CALL cget51(3,n,B,Lda,T,Lda,Q,Ldq,Q,Ldq,Work,Rwork,   &
     &                        Result(3+rsub))
                  CALL cget51(3,n,B,Lda,T,Lda,Z,Ldq,Z,Ldq,Work,Rwork,   &
     &                        Result(4+rsub))
!
!              Do test 5 and 6 (or Tests 10 and 11 when reordering):
!              check Schur form of A and compare eigenvalues with
!              diagonals.
!
                  ntest = 6 + rsub
                  temp1 = ZERO
!
                  DO j = 1 , n
                     ilabad = .FALSE.
                     temp2 = (ABS1(Alpha(j)-S(j,j))                     &
     &                       /MAX(safmin,ABS1(Alpha(j)),ABS1(S(j,j)))   &
     &                       +ABS1(Beta(j)-T(j,j))                      &
     &                       /MAX(safmin,ABS1(Beta(j)),ABS1(T(j,j))))   &
     &                       /ulp
!
                     IF ( j<n ) THEN
                        IF ( S(j+1,j)/=ZERO ) THEN
                           ilabad = .TRUE.
                           Result(5+rsub) = ulpinv
                        ENDIF
                     ENDIF
                     IF ( j>1 ) THEN
                        IF ( S(j,j-1)/=ZERO ) THEN
                           ilabad = .TRUE.
                           Result(5+rsub) = ulpinv
                        ENDIF
                     ENDIF
                     temp1 = MAX(temp1,temp2)
                     IF ( ilabad ) WRITE (Nounit,FMT=99002) j , n ,     &
     &                    jtype , ioldsd
                  ENDDO
                  Result(6+rsub) = temp1
!
                  IF ( isort>=1 ) THEN
!
!                 Do test 12
!
                     ntest = 12
                     Result(12) = ZERO
                     knteig = 0
                     DO i = 1 , n
                        IF ( clctes(Alpha(i),Beta(i)) )                 &
     &                       knteig = knteig + 1
                     ENDDO
                     IF ( sdim/=knteig ) Result(13) = ulpinv
                  ENDIF
!
               ENDDO
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
!
               ntestt = ntestt + ntest
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
                        WRITE (Nounit,FMT=99003) 'CGS'
!
!                    Matrix types
!
                        WRITE (Nounit,FMT=99004)
                        WRITE (Nounit,FMT=99005)
                        WRITE (Nounit,FMT=99006) 'Unitary'
!
!                    Tests performed
!
                        WRITE (Nounit,FMT=99007) 'unitary' , '''' ,     &
     &                         'transpose' , ('''',j=1,8)
!
                     ENDIF
                     nerrs = nerrs + 1
                     IF ( Result(jr)<10000.0 ) THEN
                        WRITE (Nounit,FMT=99008) n , jtype , ioldsd ,   &
     &                         jr , Result(jr)
                     ELSE
                        WRITE (Nounit,FMT=99009) n , jtype , ioldsd ,   &
     &                         jr , Result(jr)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
!
         ENDDO
      ENDDO
!
!     Summary
!
      CALL alasvm('CGS',Nounit,nerrs,ntestt,0)
!
      Work(1) = maxwrk
!
      RETURN
!
99001 FORMAT (' cdrges: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,', ISEED=(',4(I4,','),I5,')')
!
99002 FORMAT (' cdrges: S not in Schur form at eigenvalue ',I6,'.',/9X, &
     &        'N=',I6,', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
99003 FORMAT (/1X,A3,' -- Complex Generalized Schur from problem ',     &
     &        'driver')
!
99004 FORMAT (' Matrix types (see cdrges for details): ')
!
99005 FORMAT (' Special Matrices:',23X,'(J''=transposed Jordan block)', &
     &        /'   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ',  &
     &        '6=(diag(J'',I), diag(I,J''))',/' Diagonal Matrices:  ( ',&
     &        'D=diag(0,1,2,...) )',/'   7=(D,I)   9=(large*D, small*I',&
     &        ')  11=(large*I, small*D)  13=(large*D, large*I)',/       &
     &       '   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) '&
     &       ,' 14=(small*D, small*I)',/'  15=(D, reversed D)')
99006 FORMAT (' Matrices Rotated by Random ',A,' Matrices U, V:',       &
     &        /'  16=Transposed Jordan Blocks             19=geometric '&
     &        ,'alpha, beta=0,1',                                       &
     &        /'  17=arithm. alpha&beta             ',                  &
     &        '      20=arithmetic alpha, beta=0,1',/'  18=clustered ', &
     &        'alpha, beta=0,1            21=random alpha, beta=0,1',   &
     &        /' Large & Small Matrices:',/'  22=(large, small)   ',    &
     &        '23=(small,large)    24=(small,small)    25=(large,large)'&
     &        ,/'  26=random O(1) matrices.')
!
99007 FORMAT (/' Tests performed:  (S is Schur, T is triangular, ',     &
     &        'Q and Z are ',A,',',/19X,                                &
     &        'l and r are the appropriate left and right',/19X,        &
     &        'eigenvectors, resp., a is alpha, b is beta, and',/19X,A, &
     &        ' means ',A,'.)',/' Without ordering: ',                  &
     &        /'  1 = | A - Q S Z',A,                                   &
     &        ' | / ( |A| n ulp )      2 = | B - Q T Z',A,              &
     &        ' | / ( |B| n ulp )',/'  3 = | I - QQ',A,                 &
     &        ' | / ( n ulp )             4 = | I - ZZ',A,              &
     &        ' | / ( n ulp )',/'  5 = A is in Schur form S',           &
     &        /'  6 = difference between (alpha,beta)',                 &
     &        ' and diagonals of (S,T)',/' With ordering: ',            &
     &        /'  7 = | (A,B) - Q (S,T) Z',A,' | / ( |(A,B)| n ulp )',  &
     &        /'  8 = | I - QQ',A,                                      &
     &        ' | / ( n ulp )             9 = | I - ZZ',A,              &
     &        ' | / ( n ulp )',/' 10 = A is in Schur form S',           &
     &        /' 11 = difference between (alpha,beta) and diagonals',   &
     &        ' of (S,T)',/' 12 = SDIM is the correct number of ',      &
     &        'selected eigenvalues',/)
99008 FORMAT (' Matrix order=',I5,', type=',I2,', seed=',4(I4,','),     &
     &        ' result ',I2,' is',0P,F8.2)
99009 FORMAT (' Matrix order=',I5,', type=',I2,', seed=',4(I4,','),     &
     &        ' result ',I2,' is',1P,E10.3)
!
!     End of cdrges
!
      END SUBROUTINE cdrges
