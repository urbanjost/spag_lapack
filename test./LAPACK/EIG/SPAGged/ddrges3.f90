!*==ddrges3.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
 
!> \brief \b DDRGES3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DDRGES3( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                           NOUNIT, A, LDA, B, S, T, Q, LDQ, Z, ALPHAR,
!                           ALPHAI, BETA, WORK, LWORK, RESULT, BWORK,
!                           INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDQ, LWORK, NOUNIT, NSIZES, NTYPES
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * ), DOTYPE( * )
!       INTEGER            ISEED( 4 ), NN( * )
!       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
!      $                   B( LDA, * ), BETA( * ), Q( LDQ, * ),
!      $                   RESULT( 13 ), S( LDA, * ), T( LDA, * ),
!      $                   WORK( * ), Z( LDQ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DDRGES3 checks the nonsymmetric generalized eigenvalue (Schur form)
!> problem driver DGGES3.
!>
!> DGGES3 factors A and B as Q S Z'  and Q T Z' , where ' means
!> transpose, T is upper triangular, S is in generalized Schur form
!> (block upper triangular, with 1x1 and 2x2 blocks on the diagonal,
!> the 2x2 blocks corresponding to complex conjugate pairs of
!> generalized eigenvalues), and Q and Z are orthogonal. It also
!> computes the generalized eigenvalues (alpha(j),beta(j)), j=1,...,n,
!> Thus, w(j) = alpha(j)/beta(j) is a root of the characteristic
!> equation
!>                 det( A - w(j) B ) = 0
!> Optionally it also reorder the eigenvalues so that a selected
!> cluster of eigenvalues appears in the leading diagonal block of the
!> Schur forms.
!>
!> When DDRGES3 is called, a number of matrix "sizes" ("N's") and a
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
!> (5)   if A is in Schur form (i.e. quasi-triangular form)
!>       (no sorting of eigenvalues)
!>
!> (6)   if eigenvalues = diagonal blocks of the Schur form (S, T),
!>       i.e., test the maximum over j of D(j)  where:
!>
!>       if alpha(j) is real:
!>                     |alpha(j) - S(j,j)|        |beta(j) - T(j,j)|
!>           D(j) = ------------------------ + -----------------------
!>                  max(|alpha(j)|,|S(j,j)|)   max(|beta(j)|,|T(j,j)|)
!>
!>       if alpha(j) is complex:
!>                                 | det( s S - w T ) |
!>           D(j) = ---------------------------------------------------
!>                  ulp max( s norm(S), |w| norm(T) )*norm( s S - w T )
!>
!>       and S and T are here the 2 x 2 diagonal blocks of S and T
!>       corresponding to the j-th and j+1-th eigenvalues.
!>       (no sorting of eigenvalues)
!>
!> (7)   | (A,B) - Q (S,T) Z' | / ( | (A,B) | n ulp )
!>            (with sorting of eigenvalues).
!>
!> (8)   | I - QQ' | / ( n ulp ) (with sorting of eigenvalues).
!>
!> (9)   | I - ZZ' | / ( n ulp ) (with sorting of eigenvalues).
!>
!> (10)  if A is in Schur form (i.e. quasi-triangular form)
!>       (with sorting of eigenvalues).
!>
!> (11)  if eigenvalues = diagonal blocks of the Schur form (S, T),
!>       i.e. test the maximum over j of D(j)  where:
!>
!>       if alpha(j) is real:
!>                     |alpha(j) - S(j,j)|        |beta(j) - T(j,j)|
!>           D(j) = ------------------------ + -----------------------
!>                  max(|alpha(j)|,|S(j,j)|)   max(|beta(j)|,|T(j,j)|)
!>
!>       if alpha(j) is complex:
!>                                 | det( s S - w T ) |
!>           D(j) = ---------------------------------------------------
!>                  ulp max( s norm(S), |w| norm(T) )*norm( s S - w T )
!>
!>       and S and T are here the 2 x 2 diagonal blocks of S and T
!>       corresponding to the j-th and j+1-th eigenvalues.
!>       (with sorting of eigenvalues).
!>
!> (12)  if sorting worked and SDIM is the number of eigenvalues
!>       which were SELECTed.
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
!>          DDRGES3 does nothing.  NSIZES >= 0.
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
!>          The number of elements in DOTYPE.   If it is zero, DDRGES3
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
!>          next call to DDRGES3 to continue the same random number
!>          sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
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
!>          A is DOUBLE PRECISION array,
!>                                       dimension(LDA, max(NN))
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
!>          B is DOUBLE PRECISION array,
!>                                       dimension(LDA, max(NN))
!>          Used to hold the original B matrix.  Used as input only
!>          if NTYPES=MAXTYP+1, DOTYPE(1:MAXTYP)=.FALSE., and
!>          DOTYPE(MAXTYP+1)=.TRUE.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (LDA, max(NN))
!>          The Schur form matrix computed from A by DGGES3.  On exit, S
!>          contains the Schur form matrix corresponding to the matrix
!>          in A.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDA, max(NN))
!>          The upper triangular matrix computed from B by DGGES3.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ, max(NN))
!>          The (left) orthogonal matrix computed by DGGES3.
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
!>          Z is DOUBLE PRECISION array, dimension( LDQ, max(NN) )
!>          The (right) orthogonal matrix computed by DGGES3.
!> \endverbatim
!>
!> \param[out] ALPHAR
!> \verbatim
!>          ALPHAR is DOUBLE PRECISION array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] ALPHAI
!> \verbatim
!>          ALPHAI is DOUBLE PRECISION array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION array, dimension (max(NN))
!>
!>          The generalized eigenvalues of (A,B) computed by DGGES3.
!>          ( ALPHAR(k)+ALPHAI(k)*i ) / BETA(k) is the k-th
!>          generalized eigenvalue of A and B.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          LWORK >= MAX( 10*(N+1), 3*N*N ), where N is the largest
!>          matrix dimension.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (15)
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
!> \date February 2015
!
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DDRGES3(Nsizes,Nn,Ntypes,Dotype,Iseed,Thresh,Nounit,A, &
     &                   Lda,B,S,T,Q,Ldq,Z,Alphar,Alphai,Beta,Work,     &
     &                   Lwork,Result,Bwork,Info)
      IMPLICIT NONE
!*--DDRGES3407
!
!  -- LAPACK test routine (version 3.6.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     February 2015
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldq , Lwork , Nounit , Nsizes , Ntypes
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Bwork(*) , Dotype(*)
      INTEGER Iseed(4) , Nn(*)
      DOUBLE PRECISION A(Lda,*) , Alphai(*) , Alphar(*) , B(Lda,*) ,    &
     &                 Beta(*) , Q(Ldq,*) , Result(13) , S(Lda,*) ,     &
     &                 T(Lda,*) , Work(*) , Z(Ldq,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      INTEGER MAXTYP
      PARAMETER (MAXTYP=26)
!     ..
!     .. Local Scalars ..
      LOGICAL badnn , ilabad
      CHARACTER sort
      INTEGER i , i1 , iadd , ierr , iinfo , in , isort , j , jc , jr , &
     &        jsize , jtype , knteig , maxwrk , minwrk , mtypes , n ,   &
     &        n1 , nb , nerrs , nmats , nmax , ntest , ntestt , rsub ,  &
     &        sdim
      DOUBLE PRECISION safmax , safmin , temp1 , temp2 , ulp , ulpinv
!     ..
!     .. Local Arrays ..
      INTEGER iasign(MAXTYP) , ibsign(MAXTYP) , ioldsd(4) , kadd(6) ,   &
     &        kamagn(MAXTYP) , katype(MAXTYP) , kazero(MAXTYP) ,        &
     &        kbmagn(MAXTYP) , kbtype(MAXTYP) , kbzero(MAXTYP) ,        &
     &        kclass(MAXTYP) , ktrian(MAXTYP) , kz1(6) , kz2(6)
      DOUBLE PRECISION rmagn(0:3)
!     ..
!     .. External Functions ..
      LOGICAL DLCTES
      INTEGER ILAENV
      DOUBLE PRECISION DLAMCH , DLARND
      EXTERNAL DLCTES , ILAENV , DLAMCH , DLARND
!     ..
!     .. External Subroutines ..
      EXTERNAL ALASVM , DGET51 , DGET53 , DGET54 , DGGES3 , DLABAD ,    &
     &         DLACPY , DLARFG , DLASET , DLATM4 , DORM2R , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN , SIGN
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
      DATA iasign/6*0 , 2 , 0 , 2*2 , 2*0 , 3*2 , 0 , 2 , 3*0 , 5*2 , 0/
      DATA ibsign/7*0 , 2 , 2*0 , 2*2 , 2*0 , 2 , 0 , 2 , 9*0/
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
!       following subroutine, as returned by ILAENV.
!
      minwrk = 1
      IF ( Info==0 .AND. Lwork>=1 ) THEN
         minwrk = MAX(10*(nmax+1),3*nmax*nmax)
         nb = MAX(1,ILAENV(1,'DGEQRF',' ',nmax,nmax,-1,-1),             &
     &        ILAENV(1,'DORMQR','LT',nmax,nmax,nmax,-1),                &
     &        ILAENV(1,'DORGQR',' ',nmax,nmax,nmax,-1))
         maxwrk = MAX(10*(nmax+1),2*nmax+nmax*nb,3*nmax*nmax)
         Work(1) = maxwrk
      ENDIF
!
      IF ( Lwork<minwrk ) Info = -20
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DDRGES3',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( Nsizes==0 .OR. Ntypes==0 ) RETURN
!
      safmin = DLAMCH('Safe minimum')
      ulp = DLAMCH('Epsilon')*DLAMCH('Base')
      safmin = safmin/ulp
      safmax = ONE/safmin
      CALL DLABAD(safmin,safmax)
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
         rmagn(2) = safmax*ulp/DBLE(n1)
         rmagn(3) = safmin*ulpinv*DBLE(n1)
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
!           KZLASS: =1 means w/o rotation, =2 means w/ rotation,
!                   =3 means random.
!           KATYPE: the "type" to be passed to DLATM4 for computing A.
!           KAZERO: the pattern of zeros on the diagonal for A:
!                   =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
!                   =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
!                   =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
!                   non-zero entries.)
!           KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
!                   =2: large, =3: small.
!           IASIGN: 1 if the diagonal elements of A are to be
!                   multiplied by a random magnitude 1 number, =2 if
!                   randomly chosen diagonal blocks are to be rotated
!                   to form 2x2 blocks.
!           KBTYPE, KBZERO, KBMAGN, IBSIGN: the same, but for B.
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
                        IF ( in/=n ) CALL DLASET('Full',n,n,ZERO,ZERO,A,&
     &                       Lda)
                     ELSE
                        in = n
                     ENDIF
                     CALL DLATM4(katype(jtype),in,kz1(kazero(jtype)),   &
     &                           kz2(kazero(jtype)),iasign(jtype),      &
     &                           rmagn(kamagn(jtype)),ulp,              &
     &                           rmagn(ktrian(jtype)*kamagn(jtype)),2,  &
     &                           Iseed,A,Lda)
                     iadd = kadd(kazero(jtype))
                     IF ( iadd>0 .AND. iadd<=n ) A(iadd,iadd) = ONE
!
!              Generate B (w/o rotation)
!
                     IF ( ABS(kbtype(jtype))==3 ) THEN
                        in = 2*((n-1)/2) + 1
                        IF ( in/=n ) CALL DLASET('Full',n,n,ZERO,ZERO,B,&
     &                       Lda)
                     ELSE
                        in = n
                     ENDIF
                     CALL DLATM4(kbtype(jtype),in,kz1(kbzero(jtype)),   &
     &                           kz2(kbzero(jtype)),ibsign(jtype),      &
     &                           rmagn(kbmagn(jtype)),ONE,              &
     &                           rmagn(ktrian(jtype)*kbmagn(jtype)),2,  &
     &                           Iseed,B,Lda)
                     iadd = kadd(kbzero(jtype))
                     IF ( iadd/=0 .AND. iadd<=n ) B(iadd,iadd) = ONE
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
                              Q(jr,jc) = DLARND(3,Iseed)
                              Z(jr,jc) = DLARND(3,Iseed)
                           ENDDO
                           CALL DLARFG(n+1-jc,Q(jc,jc),Q(jc+1,jc),1,    &
     &                                 Work(jc))
                           Work(2*n+jc) = SIGN(ONE,Q(jc,jc))
                           Q(jc,jc) = ONE
                           CALL DLARFG(n+1-jc,Z(jc,jc),Z(jc+1,jc),1,    &
     &                                 Work(n+jc))
                           Work(3*n+jc) = SIGN(ONE,Z(jc,jc))
                           Z(jc,jc) = ONE
                        ENDDO
                        Q(n,n) = ONE
                        Work(n) = ZERO
                        Work(3*n) = SIGN(ONE,DLARND(2,Iseed))
                        Z(n,n) = ONE
                        Work(2*n) = ZERO
                        Work(4*n) = SIGN(ONE,DLARND(2,Iseed))
!
!                 Apply the diagonal matrices
!
                        DO jc = 1 , n
                           DO jr = 1 , n
                              A(jr,jc) = Work(2*n+jr)*Work(3*n+jc)      &
     &                           *A(jr,jc)
                              B(jr,jc) = Work(2*n+jr)*Work(3*n+jc)      &
     &                           *B(jr,jc)
                           ENDDO
                        ENDDO
                        CALL DORM2R('L','N',n,n,n-1,Q,Ldq,Work,A,Lda,   &
     &                              Work(2*n+1),iinfo)
                        IF ( iinfo==0 ) THEN
                           CALL DORM2R('R','T',n,n,n-1,Z,Ldq,Work(n+1), &
     &                                 A,Lda,Work(2*n+1),iinfo)
                           IF ( iinfo==0 ) THEN
                              CALL DORM2R('L','N',n,n,n-1,Q,Ldq,Work,B, &
     &                           Lda,Work(2*n+1),iinfo)
                              IF ( iinfo==0 ) THEN
                                 CALL DORM2R('R','T',n,n,n-1,Z,Ldq,     &
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
     &                                *DLARND(2,Iseed)
                           B(jr,jc) = rmagn(kbmagn(jtype))              &
     &                                *DLARND(2,Iseed)
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
!              Call DGGES3 to compute H, T, Q, Z, alpha, and beta.
!
                  CALL DLACPY('Full',n,n,A,Lda,S,Lda)
                  CALL DLACPY('Full',n,n,B,Lda,T,Lda)
                  ntest = 1 + rsub + isort
                  Result(1+rsub+isort) = ulpinv
                  CALL DGGES3('V','V',sort,DLCTES,n,S,Lda,T,Lda,sdim,   &
     &                        Alphar,Alphai,Beta,Q,Ldq,Z,Ldq,Work,Lwork,&
     &                        Bwork,iinfo)
                  IF ( iinfo/=0 .AND. iinfo/=n+2 ) THEN
                     Result(1+rsub+isort) = ulpinv
                     WRITE (Nounit,FMT=99001) 'DGGES3' , iinfo , n ,    &
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
                     CALL DGET51(1,n,A,Lda,S,Lda,Q,Ldq,Z,Ldq,Work,      &
     &                           Result(1))
                     CALL DGET51(1,n,B,Lda,T,Lda,Q,Ldq,Z,Ldq,Work,      &
     &                           Result(2))
                  ELSE
                     CALL DGET54(n,A,Lda,B,Lda,S,Lda,T,Lda,Q,Ldq,Z,Ldq, &
     &                           Work,Result(7))
                  ENDIF
                  CALL DGET51(3,n,A,Lda,T,Lda,Q,Ldq,Q,Ldq,Work,         &
     &                        Result(3+rsub))
                  CALL DGET51(3,n,B,Lda,T,Lda,Z,Ldq,Z,Ldq,Work,         &
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
                     IF ( Alphai(j)==ZERO ) THEN
                        temp2 = (ABS(Alphar(j)-S(j,j))                  &
     &                          /MAX(safmin,ABS(Alphar(j)),ABS(S(j,j))) &
     &                          +ABS(Beta(j)-T(j,j))                    &
     &                          /MAX(safmin,ABS(Beta(j)),ABS(T(j,j))))  &
     &                          /ulp
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
!
                     ELSE
                        IF ( Alphai(j)>ZERO ) THEN
                           i1 = j
                        ELSE
                           i1 = j - 1
                        ENDIF
                        IF ( i1<=0 .OR. i1>=n ) THEN
                           ilabad = .TRUE.
                        ELSEIF ( i1<n-1 ) THEN
                           IF ( S(i1+2,i1+1)/=ZERO ) THEN
                              ilabad = .TRUE.
                              Result(5+rsub) = ulpinv
                           ENDIF
                        ELSEIF ( i1>1 ) THEN
                           IF ( S(i1,i1-1)/=ZERO ) THEN
                              ilabad = .TRUE.
                              Result(5+rsub) = ulpinv
                           ENDIF
                        ENDIF
                        IF ( .NOT.ilabad ) THEN
                           CALL DGET53(S(i1,i1),Lda,T(i1,i1),Lda,Beta(j)&
     &                                 ,Alphar(j),Alphai(j),temp2,ierr)
                           IF ( ierr>=3 ) THEN
                              WRITE (Nounit,FMT=99002) ierr , j , n ,   &
     &                               jtype , ioldsd
                              Info = ABS(ierr)
                           ENDIF
                        ELSE
                           temp2 = ulpinv
                        ENDIF
!
                     ENDIF
                     temp1 = MAX(temp1,temp2)
                     IF ( ilabad ) WRITE (Nounit,FMT=99003) j , n ,     &
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
                        IF ( DLCTES(Alphar(i),Alphai(i),Beta(i)) .OR.   &
     &                       DLCTES(Alphar(i),-Alphai(i),Beta(i)) )     &
     &                       knteig = knteig + 1
                        IF ( i<n ) THEN
                           IF ( (DLCTES(Alphar(i+1),Alphai(i+1),Beta(i+1&
     &                          )) .OR.                                 &
     &                          DLCTES(Alphar(i+1),-Alphai(i+1),Beta    &
     &                          (i+1))) .AND.                           &
     &                          (.NOT.(DLCTES(Alphar(i),Alphai(i),      &
     &                          Beta(i)) .OR.                           &
     &                          DLCTES(Alphar(i),-Alphai(i),Beta(i))))  &
     &                          .AND. iinfo/=n+2 ) Result(12) = ulpinv
                        ENDIF
                     ENDDO
                     IF ( sdim/=knteig ) Result(12) = ulpinv
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
                        WRITE (Nounit,FMT=99004) 'DGS'
!
!                    Matrix types
!
                        WRITE (Nounit,FMT=99005)
                        WRITE (Nounit,FMT=99006)
                        WRITE (Nounit,FMT=99007) 'Orthogonal'
!
!                    Tests performed
!
                        WRITE (Nounit,FMT=99008) 'orthogonal' , '''' ,  &
     &                         'transpose' , ('''',j=1,8)
!
                     ENDIF
                     nerrs = nerrs + 1
                     IF ( Result(jr)<10000.0D0 ) THEN
                        WRITE (Nounit,FMT=99009) n , jtype , ioldsd ,   &
     &                         jr , Result(jr)
                     ELSE
                        WRITE (Nounit,FMT=99010) n , jtype , ioldsd ,   &
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
      CALL ALASVM('DGS',Nounit,nerrs,ntestt,0)
!
      Work(1) = maxwrk
!
      RETURN
!
99001 FORMAT (' DDRGES3: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,      &
     &        ', JTYPE=',I6,', ISEED=(',4(I4,','),I5,')')
!
99002 FORMAT (' DDRGES3: DGET53 returned INFO=',I1,' for eigenvalue ',  &
     &        I6,'.',/9X,'N=',I6,', JTYPE=',I6,', ISEED=(',4(I4,','),I5,&
     &        ')')
!
99003 FORMAT (' DDRGES3: S not in Schur form at eigenvalue ',I6,'.',/9X,&
     &        'N=',I6,', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
99004 FORMAT (/1X,A3,' -- Real Generalized Schur form driver')
!
99005 FORMAT (' Matrix types (see DDRGES3 for details): ')
!
99006 FORMAT (' Special Matrices:',23X,'(J''=transposed Jordan block)', &
     &        /'   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ',  &
     &        '6=(diag(J'',I), diag(I,J''))',/' Diagonal Matrices:  ( ',&
     &        'D=diag(0,1,2,...) )',/'   7=(D,I)   9=(large*D, small*I',&
     &        ')  11=(large*I, small*D)  13=(large*D, large*I)',/       &
     &       '   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) '&
     &       ,' 14=(small*D, small*I)',/'  15=(D, reversed D)')
99007 FORMAT (' Matrices Rotated by Random ',A,' Matrices U, V:',       &
     &        /'  16=Transposed Jordan Blocks             19=geometric '&
     &        ,'alpha, beta=0,1',                                       &
     &        /'  17=arithm. alpha&beta             ',                  &
     &        '      20=arithmetic alpha, beta=0,1',/'  18=clustered ', &
     &        'alpha, beta=0,1            21=random alpha, beta=0,1',   &
     &        /' Large & Small Matrices:',/'  22=(large, small)   ',    &
     &        '23=(small,large)    24=(small,small)    25=(large,large)'&
     &        ,/'  26=random O(1) matrices.')
!
99008 FORMAT (/' Tests performed:  (S is Schur, T is triangular, ',     &
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
     &        /'  7 = | (A,B) - Q (S,T) Z',A,' | / ( |(A,B)| n ulp )  ',&
     &        /'  8 = | I - QQ',A,                                      &
     &        ' | / ( n ulp )            9 = | I - ZZ',A,               &
     &        ' | / ( n ulp )',/' 10 = A is in Schur form S',           &
     &        /' 11 = difference between (alpha,beta) and diagonals',   &
     &        ' of (S,T)',/' 12 = SDIM is the correct number of ',      &
     &        'selected eigenvalues',/)
99009 FORMAT (' Matrix order=',I5,', type=',I2,', seed=',4(I4,','),     &
     &        ' result ',I2,' is',0P,F8.2)
99010 FORMAT (' Matrix order=',I5,', type=',I2,', seed=',4(I4,','),     &
     &        ' result ',I2,' is',1P,D10.3)
!
!     End of DDRGES3
!
      END SUBROUTINE DDRGES3
