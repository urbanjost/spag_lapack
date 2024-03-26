!*==sdrgev3.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b sdrgev3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sdrgev3( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                           NOUNIT, A, LDA, B, S, T, Q, LDQ, Z, QE, LDQE,
!                           ALPHAR, ALPHAI, BETA, ALPHR1, ALPHI1, BETA1,
!                           WORK, LWORK, RESULT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDQ, LDQE, LWORK, NOUNIT, NSIZES,
!      $                   NTYPES
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), NN( * )
!       REAL               A( LDA, * ), ALPHAI( * ), ALPHI1( * ),
!      $                   ALPHAR( * ), ALPHR1( * ), B( LDA, * ),
!      $                   BETA( * ), BETA1( * ), Q( LDQ, * ),
!      $                   QE( LDQE, * ), RESULT( * ), S( LDA, * ),
!      $                   T( LDA, * ), WORK( * ), Z( LDQ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sdrgev3 checks the nonsymmetric generalized eigenvalue problem driver
!> routine SGGEV3.
!>
!> SGGEV3 computes for a pair of n-by-n nonsymmetric matrices (A,B) the
!> generalized eigenvalues and, optionally, the left and right
!> eigenvectors.
!>
!> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
!> or a ratio  alpha/beta = w, such that A - w*B is singular.  It is
!> usually represented as the pair (alpha,beta), as there is reasonable
!> interpretation for beta=0, and even for both being zero.
!>
!> A right generalized eigenvector corresponding to a generalized
!> eigenvalue  w  for a pair of matrices (A,B) is a vector r  such that
!> (A - wB) * r = 0.  A left generalized eigenvector is a vector l such
!> that l**H * (A - wB) = 0, where l**H is the conjugate-transpose of l.
!>
!> When sdrgev3 is called, a number of matrix "sizes" ("n's") and a
!> number of matrix "types" are specified.  For each size ("n")
!> and each type of matrix, a pair of matrices (A, B) will be generated
!> and used for testing.  For each matrix pair, the following tests
!> will be performed and compared with the threshold THRESH.
!>
!> Results from SGGEV3:
!>
!> (1)  max over all left eigenvalue/-vector pairs (alpha/beta,l) of
!>
!>      | VL**H * (beta A - alpha B) |/( ulp max(|beta A|, |alpha B|) )
!>
!>      where VL**H is the conjugate-transpose of VL.
!>
!> (2)  | |VL(i)| - 1 | / ulp and whether largest component real
!>
!>      VL(i) denotes the i-th column of VL.
!>
!> (3)  max over all left eigenvalue/-vector pairs (alpha/beta,r) of
!>
!>      | (beta A - alpha B) * VR | / ( ulp max(|beta A|, |alpha B|) )
!>
!> (4)  | |VR(i)| - 1 | / ulp and whether largest component real
!>
!>      VR(i) denotes the i-th column of VR.
!>
!> (5)  W(full) = W(partial)
!>      W(full) denotes the eigenvalues computed when both l and r
!>      are also computed, and W(partial) denotes the eigenvalues
!>      computed when only W, only W and r, or only W and l are
!>      computed.
!>
!> (6)  VL(full) = VL(partial)
!>      VL(full) denotes the left eigenvectors computed when both l
!>      and r are computed, and VL(partial) denotes the result
!>      when only l is computed.
!>
!> (7)  VR(full) = VR(partial)
!>      VR(full) denotes the right eigenvectors computed when both l
!>      and r are also computed, and VR(partial) denotes the result
!>      when only l is computed.
!>
!>
!> Test Matrices
!> ---- --------
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
!>          sdrgev3 does nothing.  NSIZES >= 0.
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
!>          The number of elements in DOTYPE.   If it is zero, sdrgev3
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
!>          next call to sdrgev3 to continue the same random number
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
!>          of the matrix.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] NOUNIT
!> \verbatim
!>          NOUNIT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns IERR not equal to 0.)
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array,
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
!>          B is REAL array,
!>                                       dimension(LDA, max(NN))
!>          Used to hold the original B matrix.  Used as input only
!>          if NTYPES=MAXTYP+1, DOTYPE(1:MAXTYP)=.FALSE., and
!>          DOTYPE(MAXTYP+1)=.TRUE.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array,
!>                                 dimension (LDA, max(NN))
!>          The Schur form matrix computed from A by SGGEV3.  On exit, S
!>          contains the Schur form matrix corresponding to the matrix
!>          in A.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is REAL array,
!>                                 dimension (LDA, max(NN))
!>          The upper triangular matrix computed from B by SGGEV3.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array,
!>                                 dimension (LDQ, max(NN))
!>          The (left) eigenvectors matrix computed by SGGEV3.
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
!>          Z is REAL array, dimension( LDQ, max(NN) )
!>          The (right) orthogonal matrix computed by SGGEV3.
!> \endverbatim
!>
!> \param[out] QE
!> \verbatim
!>          QE is REAL array, dimension( LDQ, max(NN) )
!>          QE holds the computed right or left eigenvectors.
!> \endverbatim
!>
!> \param[in] LDQE
!> \verbatim
!>          LDQE is INTEGER
!>          The leading dimension of QE. LDQE >= max(1,max(NN)).
!> \endverbatim
!>
!> \param[out] ALPHAR
!> \verbatim
!>          ALPHAR is REAL array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] ALPHAI
!> \verbatim
!>          ALPHAI is REAL array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is REAL array, dimension (max(NN))
!> \verbatim
!>          The generalized eigenvalues of (A,B) computed by SGGEV3.
!>          ( ALPHAR(k)+ALPHAI(k)*i ) / BETA(k) is the k-th
!>          generalized eigenvalue of A and B.
!> \endverbatim
!>
!> \param[out] ALPHR1
!> \verbatim
!>          ALPHR1 is REAL array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] ALPHI1
!> \verbatim
!>          ALPHI1 is REAL array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] BETA1
!> \verbatim
!>          BETA1 is REAL array, dimension (max(NN))
!>
!>          Like ALPHAR, ALPHAI, BETA, these arrays contain the
!>          eigenvalues of A and B, but those computed when SGGEV3 only
!>          computes a partial eigendecomposition, i.e. not the
!>          eigenvalues and left and right eigenvectors.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK.  LWORK >= MAX( 8*N, N*(N+1) ).
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
!>          The values computed by the tests described above.
!>          The values are currently limited to 1/ulp, to avoid overflow.
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE sdrgev3(Nsizes,Nn,Ntypes,Dotype,Iseed,Thresh,Nounit,A, &
     &                   Lda,B,S,T,Q,Ldq,Z,Qe,Ldqe,Alphar,Alphai,Beta,  &
     &                   Alphr1,Alphi1,Beta1,Work,Lwork,Result,Info)
      IMPLICIT NONE
!*--SDRGEV3411
!
!  -- LAPACK test routine (version 3.6.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     February 2015
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldq , Ldqe , Lwork , Nounit , Nsizes , Ntypes
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iseed(4) , Nn(*)
      REAL A(Lda,*) , Alphai(*) , Alphi1(*) , Alphar(*) , Alphr1(*) ,   &
     &     B(Lda,*) , Beta(*) , Beta1(*) , Q(Ldq,*) , Qe(Ldqe,*) ,      &
     &     Result(*) , S(Lda,*) , T(Lda,*) , Work(*) , Z(Ldq,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      INTEGER MAXTYP
      PARAMETER (MAXTYP=26)
!     ..
!     .. Local Scalars ..
      LOGICAL badnn
      INTEGER i , iadd , ierr , in , j , jc , jr , jsize , jtype ,      &
     &        maxwrk , minwrk , mtypes , n , n1 , nerrs , nmats , nmax ,&
     &        ntestt
      REAL safmax , safmin , ulp , ulpinv
!     ..
!     .. Local Arrays ..
      INTEGER iasign(MAXTYP) , ibsign(MAXTYP) , ioldsd(4) , kadd(6) ,   &
     &        kamagn(MAXTYP) , katype(MAXTYP) , kazero(MAXTYP) ,        &
     &        kbmagn(MAXTYP) , kbtype(MAXTYP) , kbzero(MAXTYP) ,        &
     &        kclass(MAXTYP) , ktrian(MAXTYP) , kz1(6) , kz2(6)
      REAL rmagn(0:3)
!     ..
!     .. External Functions ..
      REAL SLAMCH 
      EXTERNAL SLAMCH 
!     ..
!     .. External Subroutines ..
      EXTERNAL SGGEV3 , SLABAD , SLACPY , SLARFG ,    SLASET , SORM2R 
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , REAL , SIGN
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
      ELSEIF ( Ldqe<=1 .OR. Ldqe<nmax ) THEN
         Info = -17
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
         minwrk = MAX(1,8*nmax,nmax*(nmax+1))
         maxwrk = 7*nmax + nmax*ilaenv(1,'SGEQRF',' ',nmax,1,nmax,0)
         maxwrk = MAX(maxwrk,nmax*(nmax+1))
         Work(1) = maxwrk
      ENDIF
!
      IF ( Lwork<minwrk ) Info = -25
!
      IF ( Info/=0 ) THEN
         CALL xerbla('sdrgev3',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( Nsizes==0 .OR. Ntypes==0 ) RETURN
!
      safmin = SLAMCH('Safe minimum')
      ulp = SLAMCH('Epsilon')*SLAMCH('Base')
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
!     Loop over sizes, types
!
      ntestt = 0
      nerrs = 0
      nmats = 0
!
      DO jsize = 1 , Nsizes
         n = Nn(jsize)
         n1 = MAX(1,n)
         rmagn(2) = safmax*ulp/REAL(n1)
         rmagn(3) = safmin*ulpinv*n1
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
!
!           Save ISEED in case of an error.
!
               DO j = 1 , 4
                  ioldsd(j) = Iseed(j)
               ENDDO
!
!           Generate test matrices A and B
!
!           Description of control parameters:
!
!           KCLASS: =1 means w/o rotation, =2 means w/ rotation,
!                   =3 means random.
!           KATYPE: the "type" to be passed to slatm4 for computing A.
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
                  ierr = 0
                  IF ( kclass(jtype)<3 ) THEN
!
!              Generate A (w/o rotation)
!
                     IF ( ABS(katype(jtype))==3 ) THEN
                        in = 2*((n-1)/2) + 1
                        IF ( in/=n ) CALL SLASET('Full',n,n,ZERO,ZERO,A,&
     &                       Lda)
                     ELSE
                        in = n
                     ENDIF
                     CALL slatm4(katype(jtype),in,kz1(kazero(jtype)),   &
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
                        IF ( in/=n ) CALL SLASET('Full',n,n,ZERO,ZERO,B,&
     &                       Lda)
                     ELSE
                        in = n
                     ENDIF
                     CALL slatm4(kbtype(jtype),in,kz1(kbzero(jtype)),   &
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
                              Q(jr,jc) = slarnd(3,Iseed)
                              Z(jr,jc) = slarnd(3,Iseed)
                           ENDDO
                           CALL SLARFG(n+1-jc,Q(jc,jc),Q(jc+1,jc),1,    &
     &                                 Work(jc))
                           Work(2*n+jc) = SIGN(ONE,Q(jc,jc))
                           Q(jc,jc) = ONE
                           CALL SLARFG(n+1-jc,Z(jc,jc),Z(jc+1,jc),1,    &
     &                                 Work(n+jc))
                           Work(3*n+jc) = SIGN(ONE,Z(jc,jc))
                           Z(jc,jc) = ONE
                        ENDDO
                        Q(n,n) = ONE
                        Work(n) = ZERO
                        Work(3*n) = SIGN(ONE,slarnd(2,Iseed))
                        Z(n,n) = ONE
                        Work(2*n) = ZERO
                        Work(4*n) = SIGN(ONE,slarnd(2,Iseed))
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
                        CALL SORM2R('L','N',n,n,n-1,Q,Ldq,Work,A,Lda,   &
     &                              Work(2*n+1),ierr)
                        IF ( ierr==0 ) THEN
                           CALL SORM2R('R','T',n,n,n-1,Z,Ldq,Work(n+1), &
     &                                 A,Lda,Work(2*n+1),ierr)
                           IF ( ierr==0 ) THEN
                              CALL SORM2R('L','N',n,n,n-1,Q,Ldq,Work,B, &
     &                           Lda,Work(2*n+1),ierr)
                              IF ( ierr==0 ) THEN
                                 CALL SORM2R('R','T',n,n,n-1,Z,Ldq,     &
     &                              Work(n+1),B,Lda,Work(2*n+1),ierr)
                                 IF ( ierr/=0 ) THEN
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
     &                                *slarnd(2,Iseed)
                           B(jr,jc) = rmagn(kbmagn(jtype))              &
     &                                *slarnd(2,Iseed)
                        ENDDO
                     ENDDO
                  ENDIF
!
!
                  IF ( ierr/=0 ) THEN
                     WRITE (Nounit,FMT=99001) 'Generator' , ierr , n ,  &
     &                      jtype , ioldsd
                     Info = ABS(ierr)
                     RETURN
                  ENDIF
               ENDIF
!
!
               DO i = 1 , 7
                  Result(i) = -ONE
               ENDDO
!
!           Call SGGEV3 to compute eigenvalues and eigenvectors.
!
               CALL SLACPY(' ',n,n,A,Lda,S,Lda)
               CALL SLACPY(' ',n,n,B,Lda,T,Lda)
               CALL SGGEV3('V','V',n,S,Lda,T,Lda,Alphar,Alphai,Beta,Q,  &
     &                     Ldq,Z,Ldq,Work,Lwork,ierr)
               IF ( ierr/=0 .AND. ierr/=n+1 ) THEN
                  Result(1) = ulpinv
                  WRITE (Nounit,FMT=99001) 'SGGEV31' , ierr , n ,       &
     &                   jtype , ioldsd
                  Info = ABS(ierr)
                  GOTO 10
               ENDIF
!
!           Do the tests (1) and (2)
!
               CALL sget52(.TRUE.,n,A,Lda,B,Lda,Q,Ldq,Alphar,Alphai,    &
     &                     Beta,Work,Result(1))
               IF ( Result(2)>Thresh ) WRITE (Nounit,FMT=99002) 'Left' ,&
     &              'SGGEV31' , Result(2) , n , jtype , ioldsd
!
!           Do the tests (3) and (4)
!
               CALL sget52(.FALSE.,n,A,Lda,B,Lda,Z,Ldq,Alphar,Alphai,   &
     &                     Beta,Work,Result(3))
               IF ( Result(4)>Thresh ) WRITE (Nounit,FMT=99002)         &
     &              'Right' , 'SGGEV31' , Result(4) , n , jtype , ioldsd
!
!           Do the test (5)
!
               CALL SLACPY(' ',n,n,A,Lda,S,Lda)
               CALL SLACPY(' ',n,n,B,Lda,T,Lda)
               CALL SGGEV3('N','N',n,S,Lda,T,Lda,Alphr1,Alphi1,Beta1,Q, &
     &                     Ldq,Z,Ldq,Work,Lwork,ierr)
               IF ( ierr/=0 .AND. ierr/=n+1 ) THEN
                  Result(1) = ulpinv
                  WRITE (Nounit,FMT=99001) 'SGGEV32' , ierr , n ,       &
     &                   jtype , ioldsd
                  Info = ABS(ierr)
                  GOTO 10
               ENDIF
!
               DO j = 1 , n
                  IF ( Alphar(j)/=Alphr1(j) .OR. Beta(j)/=Beta1(j) )    &
     &                 Result(5) = ulpinv
               ENDDO
!
!           Do the test (6): Compute eigenvalues and left eigenvectors,
!           and test them
!
               CALL SLACPY(' ',n,n,A,Lda,S,Lda)
               CALL SLACPY(' ',n,n,B,Lda,T,Lda)
               CALL SGGEV3('V','N',n,S,Lda,T,Lda,Alphr1,Alphi1,Beta1,Qe,&
     &                     Ldqe,Z,Ldq,Work,Lwork,ierr)
               IF ( ierr/=0 .AND. ierr/=n+1 ) THEN
                  Result(1) = ulpinv
                  WRITE (Nounit,FMT=99001) 'SGGEV33' , ierr , n ,       &
     &                   jtype , ioldsd
                  Info = ABS(ierr)
                  GOTO 10
               ENDIF
!
               DO j = 1 , n
                  IF ( Alphar(j)/=Alphr1(j) .OR. Alphai(j)/=Alphi1(j)   &
     &                 .OR. Beta(j)/=Beta1(j) ) Result(6) = ulpinv
               ENDDO
!
               DO j = 1 , n
                  DO jc = 1 , n
                     IF ( Q(j,jc)/=Qe(j,jc) ) Result(6) = ulpinv
                  ENDDO
               ENDDO
!
!           DO the test (7): Compute eigenvalues and right eigenvectors,
!           and test them
!
               CALL SLACPY(' ',n,n,A,Lda,S,Lda)
               CALL SLACPY(' ',n,n,B,Lda,T,Lda)
               CALL SGGEV3('N','V',n,S,Lda,T,Lda,Alphr1,Alphi1,Beta1,Q, &
     &                     Ldq,Qe,Ldqe,Work,Lwork,ierr)
               IF ( ierr/=0 .AND. ierr/=n+1 ) THEN
                  Result(1) = ulpinv
                  WRITE (Nounit,FMT=99001) 'SGGEV34' , ierr , n ,       &
     &                   jtype , ioldsd
                  Info = ABS(ierr)
                  GOTO 10
               ENDIF
!
               DO j = 1 , n
                  IF ( Alphar(j)/=Alphr1(j) .OR. Alphai(j)/=Alphi1(j)   &
     &                 .OR. Beta(j)/=Beta1(j) ) Result(7) = ulpinv
               ENDDO
!
               DO j = 1 , n
                  DO jc = 1 , n
                     IF ( Z(j,jc)/=Qe(j,jc) ) Result(7) = ulpinv
                  ENDDO
               ENDDO
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
!
 10            ntestt = ntestt + 7
!
!           Print out tests which fail.
!
               DO jr = 1 , 7
                  IF ( Result(jr)>=Thresh ) THEN
!
!                 If this is the first test to fail,
!                 print a header to the data file.
!
                     IF ( nerrs==0 ) THEN
                        WRITE (Nounit,FMT=99003) 'SGV'
!
!                    Matrix types
!
                        WRITE (Nounit,FMT=99004)
                        WRITE (Nounit,FMT=99005)
                        WRITE (Nounit,FMT=99006) 'Orthogonal'
!
!                    Tests performed
!
                        WRITE (Nounit,FMT=99007)
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
      CALL alasvm('SGV',Nounit,nerrs,ntestt,0)
!
      Work(1) = maxwrk
!
      RETURN
!
99001 FORMAT (' sdrgev3: ',A,' returned INFO=',I6,'.',/3X,'N=',I6,      &
     &        ', JTYPE=',I6,', ISEED=(',4(I4,','),I5,')')
!
99002 FORMAT (' sdrgev3: ',A,' Eigenvectors from ',A,                   &
     &        ' incorrectly normalized.',/' Bits of error=',0P,G10.3,   &
     &        ',',3X,'N=',I4,', JTYPE=',I3,', ISEED=(',4(I4,','),I5,')')
!
99003 FORMAT (/1X,A3,' -- Real Generalized eigenvalue problem driver')
!
99004 FORMAT (' Matrix types (see sdrgev3 for details): ')
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
99007 FORMAT (/' Tests performed:    ',                                 &
     &        /' 1 = max | ( b A - a B )''*l | / const.,',              &
     &        /' 2 = | |VR(i)| - 1 | / ulp,',                           &
     &        /' 3 = max | ( b A - a B )*r | / const.',                 &
     &        /' 4 = | |VL(i)| - 1 | / ulp,',                           &
     &        /' 5 = 0 if W same no matter if r or l computed,',        &
     &        /' 6 = 0 if l same no matter if l computed,',             &
     &        /' 7 = 0 if r same no matter if r computed,',/1X)
99008 FORMAT (' Matrix order=',I5,', type=',I2,', seed=',4(I4,','),     &
     &        ' result ',I2,' is',0P,F8.2)
99009 FORMAT (' Matrix order=',I5,', type=',I2,', seed=',4(I4,','),     &
     &        ' result ',I2,' is',1P,E10.3)
!
!     End of sdrgev3
!
      END SUBROUTINE sdrgev3
