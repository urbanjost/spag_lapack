!*==zchkgg.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b ZCHKGG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKGG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
!                          TSTDIF, THRSHN, NOUNIT, A, LDA, B, H, T, S1,
!                          S2, P1, P2, U, LDU, V, Q, Z, ALPHA1, BETA1,
!                          ALPHA3, BETA3, EVECTL, EVECTR, WORK, LWORK,
!                          RWORK, LLWORK, RESULT, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTDIF
!       INTEGER            INFO, LDA, LDU, LWORK, NOUNIT, NSIZES, NTYPES
!       DOUBLE PRECISION   THRESH, THRSHN
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * ), LLWORK( * )
!       INTEGER            ISEED( 4 ), NN( * )
!       DOUBLE PRECISION   RESULT( 15 ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), ALPHA1( * ), ALPHA3( * ),
!      $                   B( LDA, * ), BETA1( * ), BETA3( * ),
!      $                   EVECTL( LDU, * ), EVECTR( LDU, * ),
!      $                   H( LDA, * ), P1( LDA, * ), P2( LDA, * ),
!      $                   Q( LDU, * ), S1( LDA, * ), S2( LDA, * ),
!      $                   T( LDA, * ), U( LDU, * ), V( LDU, * ),
!      $                   WORK( * ), Z( LDU, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCHKGG  checks the nonsymmetric generalized eigenvalue problem
!> routines.
!>                                H          H        H
!> ZGGHRD factors A and B as U H V  and U T V , where   means conjugate
!> transpose, H is hessenberg, T is triangular and U and V are unitary.
!>
!>                                 H          H
!> ZHGEQZ factors H and T as  Q S Z  and Q P Z , where P and S are upper
!> triangular and Q and Z are unitary.  It also computes the generalized
!> eigenvalues (alpha(1),beta(1)),...,(alpha(n),beta(n)), where
!> alpha(j)=S(j,j) and beta(j)=P(j,j) -- thus, w(j) = alpha(j)/beta(j)
!> is a root of the generalized eigenvalue problem
!>
!>     det( A - w(j) B ) = 0
!>
!> and m(j) = beta(j)/alpha(j) is a root of the essentially equivalent
!> problem
!>
!>     det( m(j) A - B ) = 0
!>
!> ZTGEVC computes the matrix L of left eigenvectors and the matrix R
!> of right eigenvectors for the matrix pair ( S, P ).  In the
!> description below,  l and r are left and right eigenvectors
!> corresponding to the generalized eigenvalues (alpha,beta).
!>
!> When ZCHKGG is called, a number of matrix "sizes" ("n's") and a
!> number of matrix "types" are specified.  For each size ("n")
!> and each type of matrix, one matrix will be generated and used
!> to test the nonsymmetric eigenroutines.  For each matrix, 13
!> tests will be performed.  The first twelve "test ratios" should be
!> small -- O(1).  They will be compared with the threshold THRESH:
!>
!>                  H
!> (1)   | A - U H V  | / ( |A| n ulp )
!>
!>                  H
!> (2)   | B - U T V  | / ( |B| n ulp )
!>
!>               H
!> (3)   | I - UU  | / ( n ulp )
!>
!>               H
!> (4)   | I - VV  | / ( n ulp )
!>
!>                  H
!> (5)   | H - Q S Z  | / ( |H| n ulp )
!>
!>                  H
!> (6)   | T - Q P Z  | / ( |T| n ulp )
!>
!>               H
!> (7)   | I - QQ  | / ( n ulp )
!>
!>               H
!> (8)   | I - ZZ  | / ( n ulp )
!>
!> (9)   max over all left eigenvalue/-vector pairs (beta/alpha,l) of
!>                           H
!>       | (beta A - alpha B) l | / ( ulp max( |beta A|, |alpha B| ) )
!>
!> (10)  max over all left eigenvalue/-vector pairs (beta/alpha,l') of
!>                           H
!>       | (beta H - alpha T) l' | / ( ulp max( |beta H|, |alpha T| ) )
!>
!>       where the eigenvectors l' are the result of passing Q to
!>       DTGEVC and back transforming (JOB='B').
!>
!> (11)  max over all right eigenvalue/-vector pairs (beta/alpha,r) of
!>
!>       | (beta A - alpha B) r | / ( ulp max( |beta A|, |alpha B| ) )
!>
!> (12)  max over all right eigenvalue/-vector pairs (beta/alpha,r') of
!>
!>       | (beta H - alpha T) r' | / ( ulp max( |beta H|, |alpha T| ) )
!>
!>       where the eigenvectors r' are the result of passing Z to
!>       DTGEVC and back transforming (JOB='B').
!>
!> The last three test ratios will usually be small, but there is no
!> mathematical requirement that they be so.  They are therefore
!> compared with THRESH only if TSTDIF is .TRUE.
!>
!> (13)  | S(Q,Z computed) - S(Q,Z not computed) | / ( |S| ulp )
!>
!> (14)  | P(Q,Z computed) - P(Q,Z not computed) | / ( |P| ulp )
!>
!> (15)  max( |alpha(Q,Z computed) - alpha(Q,Z not computed)|/|S| ,
!>            |beta(Q,Z computed) - beta(Q,Z not computed)|/|P| ) / ulp
!>
!> In addition, the normalization of L and R are checked, and compared
!> with the threshold THRSHN.
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
!> (7)  ( D, I )         where D is P*D1, P is a random unitary diagonal
!>                       matrix (i.e., with random magnitude 1 entries
!>                       on the diagonal), and D1=diag( 0, 1,..., N-1 )
!>                       (i.e., a diagonal matrix with D1(1,1)=0,
!>                       D1(2,2)=1, ..., D1(N,N)=N-1.)
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
!> (15) ( D1, D2 )        where D1=P*diag( 0, 0, 1, ..., N-3, 0 ) and
!>                        D2=Q*diag( 0, N-3, N-4,..., 1, 0, 0 ), and
!>                        P and Q are random unitary diagonal matrices.
!>           t   t
!> (16) U ( J , J ) V     where U and V are random unitary matrices.
!>
!> (17) U ( T1, T2 ) V    where T1 and T2 are upper triangular matrices
!>                        with random O(1) entries above the diagonal
!>                        and diagonal entries diag(T1) =
!>                        P*( 0, 0, 1, ..., N-3, 0 ) and diag(T2) =
!>                        Q*( 0, N-3, N-4,..., 1, 0, 0 )
!>
!> (18) U ( T1, T2 ) V    diag(T1) = ( 0, 0, 1, 1, s, ..., s, 0 )
!>                        diag(T2) = ( 0, 1, 0, 1,..., 1, 0 )
!>                        s = machine precision.
!>
!> (19) U ( T1, T2 ) V    diag(T1)=( 0,0,1,1, 1-d, ..., 1-(N-5)*d=s, 0 )
!>                        diag(T2) = ( 0, 1, 0, 1, ..., 1, 0 )
!>
!>                                                        N-5
!> (20) U ( T1, T2 ) V    diag(T1)=( 0, 0, 1, 1, a, ..., a   =s, 0 )
!>                        diag(T2) = ( 0, 1, 0, 1, ..., 1, 0, 0 )
!>
!> (21) U ( T1, T2 ) V    diag(T1)=( 0, 0, 1, r1, r2, ..., r(N-4), 0 )
!>                        diag(T2) = ( 0, 1, 0, 1, ..., 1, 0, 0 )
!>                        where r1,..., r(N-4) are random.
!>
!> (22) U ( big*T1, small*T2 ) V   diag(T1) = P*( 0, 0, 1, ..., N-3, 0 )
!>                                 diag(T2) = ( 0, 1, ..., 1, 0, 0 )
!>
!> (23) U ( small*T1, big*T2 ) V   diag(T1) = P*( 0, 0, 1, ..., N-3, 0 )
!>                                 diag(T2) = ( 0, 1, ..., 1, 0, 0 )
!>
!> (24) U ( small*T1, small*T2 ) V diag(T1) = P*( 0, 0, 1, ..., N-3, 0 )
!>                                 diag(T2) = ( 0, 1, ..., 1, 0, 0 )
!>
!> (25) U ( big*T1, big*T2 ) V     diag(T1) = P*( 0, 0, 1, ..., N-3, 0 )
!>                                 diag(T2) = ( 0, 1, ..., 1, 0, 0 )
!>
!> (26) U ( T1, T2 ) V     where T1 and T2 are random upper-triangular
!>                         matrices.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of sizes of matrices to use.  If it is zero,
!>          ZCHKGG does nothing.  It must be at least zero.
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
!>          The number of elements in DOTYPE.   If it is zero, ZCHKGG
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
!>          next call to ZCHKGG to continue the same random number
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
!> \param[in] TSTDIF
!> \verbatim
!>          TSTDIF is LOGICAL
!>          Specifies whether test ratios 13-15 will be computed and
!>          compared with THRESH.
!>          = .FALSE.: Only test ratios 1-12 will be computed and tested.
!>                     Ratios 13-15 will be set to zero.
!>          = .TRUE.:  All the test ratios 1-15 will be computed and
!>                     tested.
!> \endverbatim
!>
!> \param[in] THRSHN
!> \verbatim
!>          THRSHN is DOUBLE PRECISION
!>          Threshold for reporting eigenvector normalization error.
!>          If the normalization of any eigenvector differs from 1 by
!>          more than THRSHN*ulp, then a special error message will be
!>          printed.  (This is handled separately from the other tests,
!>          since only a compiler or programming error should cause an
!>          error message, at least if THRSHN is at least 5--10.)
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
!>          A is COMPLEX*16 array, dimension (LDA, max(NN))
!>          Used to hold the original A matrix.  Used as input only
!>          if NTYPES=MAXTYP+1, DOTYPE(1:MAXTYP)=.FALSE., and
!>          DOTYPE(MAXTYP+1)=.TRUE.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A, B, H, T, S1, P1, S2, and P2.
!>          It must be at least 1 and at least max( NN ).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDA, max(NN))
!>          Used to hold the original B matrix.  Used as input only
!>          if NTYPES=MAXTYP+1, DOTYPE(1:MAXTYP)=.FALSE., and
!>          DOTYPE(MAXTYP+1)=.TRUE.
!> \endverbatim
!>
!> \param[out] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDA, max(NN))
!>          The upper Hessenberg matrix computed from A by ZGGHRD.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDA, max(NN))
!>          The upper triangular matrix computed from B by ZGGHRD.
!> \endverbatim
!>
!> \param[out] S1
!> \verbatim
!>          S1 is COMPLEX*16 array, dimension (LDA, max(NN))
!>          The Schur (upper triangular) matrix computed from H by ZHGEQZ
!>          when Q and Z are also computed.
!> \endverbatim
!>
!> \param[out] S2
!> \verbatim
!>          S2 is COMPLEX*16 array, dimension (LDA, max(NN))
!>          The Schur (upper triangular) matrix computed from H by ZHGEQZ
!>          when Q and Z are not computed.
!> \endverbatim
!>
!> \param[out] P1
!> \verbatim
!>          P1 is COMPLEX*16 array, dimension (LDA, max(NN))
!>          The upper triangular matrix computed from T by ZHGEQZ
!>          when Q and Z are also computed.
!> \endverbatim
!>
!> \param[out] P2
!> \verbatim
!>          P2 is COMPLEX*16 array, dimension (LDA, max(NN))
!>          The upper triangular matrix computed from T by ZHGEQZ
!>          when Q and Z are not computed.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension (LDU, max(NN))
!>          The (left) unitary matrix computed by ZGGHRD.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U, V, Q, Z, EVECTL, and EVEZTR.  It
!>          must be at least 1 and at least max( NN ).
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (LDU, max(NN))
!>          The (right) unitary matrix computed by ZGGHRD.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDU, max(NN))
!>          The (left) unitary matrix computed by ZHGEQZ.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDU, max(NN))
!>          The (left) unitary matrix computed by ZHGEQZ.
!> \endverbatim
!>
!> \param[out] ALPHA1
!> \verbatim
!>          ALPHA1 is COMPLEX*16 array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] BETA1
!> \verbatim
!>          BETA1 is COMPLEX*16 array, dimension (max(NN))
!>          The generalized eigenvalues of (A,B) computed by ZHGEQZ
!>          when Q, Z, and the full Schur matrices are computed.
!> \endverbatim
!>
!> \param[out] ALPHA3
!> \verbatim
!>          ALPHA3 is COMPLEX*16 array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] BETA3
!> \verbatim
!>          BETA3 is COMPLEX*16 array, dimension (max(NN))
!>          The generalized eigenvalues of (A,B) computed by ZHGEQZ
!>          when neither Q, Z, nor the Schur matrices are computed.
!> \endverbatim
!>
!> \param[out] EVECTL
!> \verbatim
!>          EVECTL is COMPLEX*16 array, dimension (LDU, max(NN))
!>          The (lower triangular) left eigenvector matrix for the
!>          matrices in S1 and P1.
!> \endverbatim
!>
!> \param[out] EVECTR
!> \verbatim
!>          EVECTR is COMPLEX*16 array, dimension (LDU, max(NN))
!>          The (upper triangular) right eigenvector matrix for the
!>          matrices in S1 and P1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The number of entries in WORK.  This must be at least
!>          max( 4*N, 2 * N**2, 1 ), for all N=NN(j).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (2*max(NN))
!> \endverbatim
!>
!> \param[out] LLWORK
!> \verbatim
!>          LLWORK is LOGICAL array, dimension (max(NN))
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (15)
!>          The values computed by the tests described above.
!>          The values are currently limited to 1/ulp, to avoid
!>          overflow.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZCHKGG(Nsizes,Nn,Ntypes,Dotype,Iseed,Thresh,Tstdif,    &
     &                  Thrshn,Nounit,A,Lda,B,H,T,S1,S2,P1,P2,U,Ldu,V,Q,&
     &                  Z,Alpha1,Beta1,Alpha3,Beta3,Evectl,Evectr,Work, &
     &                  Lwork,Rwork,Llwork,Result,Info)
      IMPLICIT NONE
!*--ZCHKGG506
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tstdif
      INTEGER Info , Lda , Ldu , Lwork , Nounit , Nsizes , Ntypes
      DOUBLE PRECISION Thresh , Thrshn
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*) , Llwork(*)
      INTEGER Iseed(4) , Nn(*)
      DOUBLE PRECISION Result(15) , Rwork(*)
      COMPLEX*16 A(Lda,*) , Alpha1(*) , Alpha3(*) , B(Lda,*) , Beta1(*) &
     &           , Beta3(*) , Evectl(Ldu,*) , Evectr(Ldu,*) , H(Lda,*) ,&
     &           P1(Lda,*) , P2(Lda,*) , Q(Ldu,*) , S1(Lda,*) ,         &
     &           S2(Lda,*) , T(Lda,*) , U(Ldu,*) , V(Ldu,*) , Work(*) , &
     &           Z(Ldu,*)
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
      PARAMETER (MAXTYP=26)
!     ..
!     .. Local Scalars ..
      LOGICAL badnn
      INTEGER i1 , iadd , iinfo , in , j , jc , jr , jsize , jtype ,    &
     &        lwkopt , mtypes , n , n1 , nerrs , nmats , nmax , ntest , &
     &        ntestt
      DOUBLE PRECISION anorm , bnorm , safmax , safmin , temp1 , temp2 ,&
     &                 ulp , ulpinv
      COMPLEX*16 ctemp
!     ..
!     .. Local Arrays ..
      LOGICAL lasign(MAXTYP) , lbsign(MAXTYP)
      INTEGER ioldsd(4) , kadd(6) , kamagn(MAXTYP) , katype(MAXTYP) ,   &
     &        kazero(MAXTYP) , kbmagn(MAXTYP) , kbtype(MAXTYP) ,        &
     &        kbzero(MAXTYP) , kclass(MAXTYP) , ktrian(MAXTYP) , kz1(6) &
     &        , kz2(6)
      DOUBLE PRECISION dumma(4) , rmagn(0:3)
      COMPLEX*16 cdumma(4)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE
      COMPLEX*16 ZLARND
      EXTERNAL DLAMCH , ZLANGE , ZLARND
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD , DLASUM , XERBLA , ZGEQR2 , ZGET51 , ZGET52 ,    &
     &         ZGGHRD , ZHGEQZ , ZLACPY , ZLARFG , ZLASET , ZLATM4 ,    &
     &         ZTGEVC , ZUNM2R
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DCONJG , MAX , MIN , SIGN
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
      lwkopt = MAX(2*nmax*nmax,4*nmax,1)
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
         Info = -10
      ELSEIF ( Ldu<=1 .OR. Ldu<nmax ) THEN
         Info = -19
      ELSEIF ( lwkopt>Lwork ) THEN
         Info = -30
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZCHKGG',-Info)
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
!     Loop over sizes, types
!
      ntestt = 0
      nerrs = 0
      nmats = 0
!
      DO jsize = 1 , Nsizes
         n = Nn(jsize)
         n1 = MAX(1,n)
         rmagn(2) = safmax*ulp/DBLE(n1)
         rmagn(3) = safmin*ulpinv*n1
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
            DO j = 1 , 15
               Result(j) = ZERO
            ENDDO
!
!           Compute A and B
!
!           Description of control parameters:
!
!           KZLASS: =1 means w/o rotation, =2 means w/ rotation,
!                   =3 means random.
!           KATYPE: the "type" to be passed to ZLATM4 for computing A.
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
!           RMAGN:  used to implement KAMAGN and KBMAGN.
!
            IF ( mtypes>MAXTYP ) GOTO 40
            iinfo = 0
            IF ( kclass(jtype)<3 ) THEN
!
!              Generate A (w/o rotation)
!
               IF ( ABS(katype(jtype))==3 ) THEN
                  in = 2*((n-1)/2) + 1
                  IF ( in/=n ) CALL ZLASET('Full',n,n,CZERO,CZERO,A,Lda)
               ELSE
                  in = n
               ENDIF
               CALL ZLATM4(katype(jtype),in,kz1(kazero(jtype)),         &
     &                     kz2(kazero(jtype)),lasign(jtype),            &
     &                     rmagn(kamagn(jtype)),ulp,                    &
     &                     rmagn(ktrian(jtype)*kamagn(jtype)),4,Iseed,A,&
     &                     Lda)
               iadd = kadd(kazero(jtype))
               IF ( iadd>0 .AND. iadd<=n ) A(iadd,iadd)                 &
     &              = rmagn(kamagn(jtype))
!
!              Generate B (w/o rotation)
!
               IF ( ABS(kbtype(jtype))==3 ) THEN
                  in = 2*((n-1)/2) + 1
                  IF ( in/=n ) CALL ZLASET('Full',n,n,CZERO,CZERO,B,Lda)
               ELSE
                  in = n
               ENDIF
               CALL ZLATM4(kbtype(jtype),in,kz1(kbzero(jtype)),         &
     &                     kz2(kbzero(jtype)),lbsign(jtype),            &
     &                     rmagn(kbmagn(jtype)),ONE,                    &
     &                     rmagn(ktrian(jtype)*kbmagn(jtype)),4,Iseed,B,&
     &                     Lda)
               iadd = kadd(kbzero(jtype))
               IF ( iadd/=0 ) B(iadd,iadd) = rmagn(kbmagn(jtype))
!
               IF ( kclass(jtype)==2 .AND. n>0 ) THEN
!
!                 Include rotations
!
!                 Generate U, V as Householder transformations times a
!                 diagonal matrix.  (Note that ZLARFG makes U(j,j) and
!                 V(j,j) real.)
!
                  DO jc = 1 , n - 1
                     DO jr = jc , n
                        U(jr,jc) = ZLARND(3,Iseed)
                        V(jr,jc) = ZLARND(3,Iseed)
                     ENDDO
                     CALL ZLARFG(n+1-jc,U(jc,jc),U(jc+1,jc),1,Work(jc))
                     Work(2*n+jc) = SIGN(ONE,DBLE(U(jc,jc)))
                     U(jc,jc) = CONE
                     CALL ZLARFG(n+1-jc,V(jc,jc),V(jc+1,jc),1,Work(n+jc)&
     &                           )
                     Work(3*n+jc) = SIGN(ONE,DBLE(V(jc,jc)))
                     V(jc,jc) = CONE
                  ENDDO
                  ctemp = ZLARND(3,Iseed)
                  U(n,n) = CONE
                  Work(n) = CZERO
                  Work(3*n) = ctemp/ABS(ctemp)
                  ctemp = ZLARND(3,Iseed)
                  V(n,n) = CONE
                  Work(2*n) = CZERO
                  Work(4*n) = ctemp/ABS(ctemp)
!
!                 Apply the diagonal matrices
!
                  DO jc = 1 , n
                     DO jr = 1 , n
                        A(jr,jc) = Work(2*n+jr)*DCONJG(Work(3*n+jc))    &
     &                             *A(jr,jc)
                        B(jr,jc) = Work(2*n+jr)*DCONJG(Work(3*n+jc))    &
     &                             *B(jr,jc)
                     ENDDO
                  ENDDO
                  CALL ZUNM2R('L','N',n,n,n-1,U,Ldu,Work,A,Lda,         &
     &                        Work(2*n+1),iinfo)
                  IF ( iinfo/=0 ) GOTO 20
                  CALL ZUNM2R('R','C',n,n,n-1,V,Ldu,Work(n+1),A,Lda,    &
     &                        Work(2*n+1),iinfo)
                  IF ( iinfo/=0 ) GOTO 20
                  CALL ZUNM2R('L','N',n,n,n-1,U,Ldu,Work,B,Lda,         &
     &                        Work(2*n+1),iinfo)
                  IF ( iinfo/=0 ) GOTO 20
                  CALL ZUNM2R('R','C',n,n,n-1,V,Ldu,Work(n+1),B,Lda,    &
     &                        Work(2*n+1),iinfo)
                  IF ( iinfo/=0 ) GOTO 20
               ENDIF
            ELSE
!
!              Random matrices
!
               DO jc = 1 , n
                  DO jr = 1 , n
                     A(jr,jc) = rmagn(kamagn(jtype))*ZLARND(4,Iseed)
                     B(jr,jc) = rmagn(kbmagn(jtype))*ZLARND(4,Iseed)
                  ENDDO
               ENDDO
            ENDIF
!
            anorm = ZLANGE('1',n,n,A,Lda,Rwork)
            bnorm = ZLANGE('1',n,n,B,Lda,Rwork)
!
!
 20         IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'Generator' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               RETURN
            ENDIF
!
!
!           Call ZGEQR2, ZUNM2R, and ZGGHRD to compute H, T, U, and V
!
 40         CALL ZLACPY(' ',n,n,A,Lda,H,Lda)
            CALL ZLACPY(' ',n,n,B,Lda,T,Lda)
            ntest = 1
            Result(1) = ulpinv
!
            CALL ZGEQR2(n,n,T,Lda,Work,Work(n+1),iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZGEQR2' , iinfo , n , jtype ,  &
     &                                  ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
!
            CALL ZUNM2R('L','C',n,n,n,T,Lda,Work,H,Lda,Work(n+1),iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZUNM2R' , iinfo , n , jtype ,  &
     &                                  ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
!
            CALL ZLASET('Full',n,n,CZERO,CONE,U,Ldu)
            CALL ZUNM2R('R','N',n,n,n,T,Lda,Work,U,Ldu,Work(n+1),iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZUNM2R' , iinfo , n , jtype ,  &
     &                                  ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
!
            CALL ZGGHRD('V','I',n,1,n,H,Lda,T,Lda,U,Ldu,V,Ldu,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZGGHRD' , iinfo , n , jtype ,  &
     &                                  ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
            ntest = 4
!
!           Do tests 1--4
!
            CALL ZGET51(1,n,A,Lda,H,Lda,U,Ldu,V,Ldu,Work,Rwork,Result(1)&
     &                  )
            CALL ZGET51(1,n,B,Lda,T,Lda,U,Ldu,V,Ldu,Work,Rwork,Result(2)&
     &                  )
            CALL ZGET51(3,n,B,Lda,T,Lda,U,Ldu,U,Ldu,Work,Rwork,Result(3)&
     &                  )
            CALL ZGET51(3,n,B,Lda,T,Lda,V,Ldu,V,Ldu,Work,Rwork,Result(4)&
     &                  )
!
!           Call ZHGEQZ to compute S1, P1, S2, P2, Q, and Z, do tests.
!
!           Compute T1 and UZ
!
!           Eigenvalues only
!
            CALL ZLACPY(' ',n,n,H,Lda,S2,Lda)
            CALL ZLACPY(' ',n,n,T,Lda,P2,Lda)
            ntest = 5
            Result(5) = ulpinv
!
            CALL ZHGEQZ('E','N','N',n,1,n,S2,Lda,P2,Lda,Alpha3,Beta3,Q, &
     &                  Ldu,Z,Ldu,Work,Lwork,Rwork,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZHGEQZ(E)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
!
!           Eigenvalues and Full Schur Form
!
            CALL ZLACPY(' ',n,n,H,Lda,S2,Lda)
            CALL ZLACPY(' ',n,n,T,Lda,P2,Lda)
!
            CALL ZHGEQZ('S','N','N',n,1,n,S2,Lda,P2,Lda,Alpha1,Beta1,Q, &
     &                  Ldu,Z,Ldu,Work,Lwork,Rwork,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZHGEQZ(S)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
!
!           Eigenvalues, Schur Form, and Schur Vectors
!
            CALL ZLACPY(' ',n,n,H,Lda,S1,Lda)
            CALL ZLACPY(' ',n,n,T,Lda,P1,Lda)
!
            CALL ZHGEQZ('S','I','I',n,1,n,S1,Lda,P1,Lda,Alpha1,Beta1,Q, &
     &                  Ldu,Z,Ldu,Work,Lwork,Rwork,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZHGEQZ(V)' , iinfo , n ,       &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
!
            ntest = 8
!
!           Do Tests 5--8
!
            CALL ZGET51(1,n,H,Lda,S1,Lda,Q,Ldu,Z,Ldu,Work,Rwork,        &
     &                  Result(5))
            CALL ZGET51(1,n,T,Lda,P1,Lda,Q,Ldu,Z,Ldu,Work,Rwork,        &
     &                  Result(6))
            CALL ZGET51(3,n,T,Lda,P1,Lda,Q,Ldu,Q,Ldu,Work,Rwork,        &
     &                  Result(7))
            CALL ZGET51(3,n,T,Lda,P1,Lda,Z,Ldu,Z,Ldu,Work,Rwork,        &
     &                  Result(8))
!
!           Compute the Left and Right Eigenvectors of (S1,P1)
!
!           9: Compute the left eigenvector Matrix without
!              back transforming:
!
            ntest = 9
            Result(9) = ulpinv
!
!           To test "SELECT" option, compute half of the eigenvectors
!           in one call, and half in another
!
            i1 = n/2
            DO j = 1 , i1
               Llwork(j) = .TRUE.
            ENDDO
            DO j = i1 + 1 , n
               Llwork(j) = .FALSE.
            ENDDO
!
            CALL ZTGEVC('L','S',Llwork,n,S1,Lda,P1,Lda,Evectl,Ldu,      &
     &                  cdumma,Ldu,n,in,Work,Rwork,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZTGEVC(L,S1)' , iinfo , n ,    &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
!
            i1 = in
            DO j = 1 , i1
               Llwork(j) = .FALSE.
            ENDDO
            DO j = i1 + 1 , n
               Llwork(j) = .TRUE.
            ENDDO
!
            CALL ZTGEVC('L','S',Llwork,n,S1,Lda,P1,Lda,Evectl(1,i1+1),  &
     &                  Ldu,cdumma,Ldu,n,in,Work,Rwork,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZTGEVC(L,S2)' , iinfo , n ,    &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
!
            CALL ZGET52(.TRUE.,n,S1,Lda,P1,Lda,Evectl,Ldu,Alpha1,Beta1, &
     &                  Work,Rwork,dumma(1))
            Result(9) = dumma(1)
            IF ( dumma(2)>Thrshn ) WRITE (Nounit,FMT=99002) 'Left' ,    &
     &           'ZTGEVC(HOWMNY=S)' , dumma(2) , n , jtype , ioldsd
!
!           10: Compute the left eigenvector Matrix with
!               back transforming:
!
            ntest = 10
            Result(10) = ulpinv
            CALL ZLACPY('F',n,n,Q,Ldu,Evectl,Ldu)
            CALL ZTGEVC('L','B',Llwork,n,S1,Lda,P1,Lda,Evectl,Ldu,      &
     &                  cdumma,Ldu,n,in,Work,Rwork,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZTGEVC(L,B)' , iinfo , n ,     &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
!
            CALL ZGET52(.TRUE.,n,H,Lda,T,Lda,Evectl,Ldu,Alpha1,Beta1,   &
     &                  Work,Rwork,dumma(1))
            Result(10) = dumma(1)
            IF ( dumma(2)>Thrshn ) WRITE (Nounit,FMT=99002) 'Left' ,    &
     &           'ZTGEVC(HOWMNY=B)' , dumma(2) , n , jtype , ioldsd
!
!           11: Compute the right eigenvector Matrix without
!               back transforming:
!
            ntest = 11
            Result(11) = ulpinv
!
!           To test "SELECT" option, compute half of the eigenvectors
!           in one call, and half in another
!
            i1 = n/2
            DO j = 1 , i1
               Llwork(j) = .TRUE.
            ENDDO
            DO j = i1 + 1 , n
               Llwork(j) = .FALSE.
            ENDDO
!
            CALL ZTGEVC('R','S',Llwork,n,S1,Lda,P1,Lda,cdumma,Ldu,      &
     &                  Evectr,Ldu,n,in,Work,Rwork,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZTGEVC(R,S1)' , iinfo , n ,    &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
!
            i1 = in
            DO j = 1 , i1
               Llwork(j) = .FALSE.
            ENDDO
            DO j = i1 + 1 , n
               Llwork(j) = .TRUE.
            ENDDO
!
            CALL ZTGEVC('R','S',Llwork,n,S1,Lda,P1,Lda,cdumma,Ldu,      &
     &                  Evectr(1,i1+1),Ldu,n,in,Work,Rwork,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZTGEVC(R,S2)' , iinfo , n ,    &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
!
            CALL ZGET52(.FALSE.,n,S1,Lda,P1,Lda,Evectr,Ldu,Alpha1,Beta1,&
     &                  Work,Rwork,dumma(1))
            Result(11) = dumma(1)
            IF ( dumma(2)>Thresh ) WRITE (Nounit,FMT=99002) 'Right' ,   &
     &           'ZTGEVC(HOWMNY=S)' , dumma(2) , n , jtype , ioldsd
!
!           12: Compute the right eigenvector Matrix with
!               back transforming:
!
            ntest = 12
            Result(12) = ulpinv
            CALL ZLACPY('F',n,n,Z,Ldu,Evectr,Ldu)
            CALL ZTGEVC('R','B',Llwork,n,S1,Lda,P1,Lda,cdumma,Ldu,      &
     &                  Evectr,Ldu,n,in,Work,Rwork,iinfo)
            IF ( iinfo/=0 ) THEN
               WRITE (Nounit,FMT=99001) 'ZTGEVC(R,B)' , iinfo , n ,     &
     &                                  jtype , ioldsd
               Info = ABS(iinfo)
               GOTO 60
            ENDIF
!
            CALL ZGET52(.FALSE.,n,H,Lda,T,Lda,Evectr,Ldu,Alpha1,Beta1,  &
     &                  Work,Rwork,dumma(1))
            Result(12) = dumma(1)
            IF ( dumma(2)>Thresh ) WRITE (Nounit,FMT=99002) 'Right' ,   &
     &           'ZTGEVC(HOWMNY=B)' , dumma(2) , n , jtype , ioldsd
!
!           Tests 13--15 are done only on request
!
            IF ( Tstdif ) THEN
!
!              Do Tests 13--14
!
               CALL ZGET51(2,n,S1,Lda,S2,Lda,Q,Ldu,Z,Ldu,Work,Rwork,    &
     &                     Result(13))
               CALL ZGET51(2,n,P1,Lda,P2,Lda,Q,Ldu,Z,Ldu,Work,Rwork,    &
     &                     Result(14))
!
!              Do Test 15
!
               temp1 = ZERO
               temp2 = ZERO
               DO j = 1 , n
                  temp1 = MAX(temp1,ABS(Alpha1(j)-Alpha3(j)))
                  temp2 = MAX(temp2,ABS(Beta1(j)-Beta3(j)))
               ENDDO
!
               temp1 = temp1/MAX(safmin,ulp*MAX(temp1,anorm))
               temp2 = temp2/MAX(safmin,ulp*MAX(temp2,bnorm))
               Result(15) = MAX(temp1,temp2)
               ntest = 15
            ELSE
               Result(13) = ZERO
               Result(14) = ZERO
               Result(15) = ZERO
               ntest = 12
            ENDIF
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
!
 60         ntestt = ntestt + ntest
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
                     WRITE (Nounit,FMT=99003) 'ZGG'
!
!                    Matrix types
!
                     WRITE (Nounit,FMT=99004)
                     WRITE (Nounit,FMT=99005)
                     WRITE (Nounit,FMT=99006) 'Unitary'
!
!                    Tests performed
!
                     WRITE (Nounit,FMT=99007) 'unitary' , '*' ,         &
     &                      'conjugate transpose' , ('*',j=1,10)
!
                  ENDIF
                  nerrs = nerrs + 1
                  IF ( Result(jr)<10000.0D0 ) THEN
                     WRITE (Nounit,FMT=99008) n , jtype , ioldsd , jr , &
     &                      Result(jr)
                  ELSE
                     WRITE (Nounit,FMT=99009) n , jtype , ioldsd , jr , &
     &                      Result(jr)
                  ENDIF
               ENDIF
            ENDDO
!
         ENDDO
      ENDDO
!
!     Summary
!
      CALL DLASUM('ZGG',Nounit,nerrs,ntestt)
      RETURN
!
99001 FORMAT (' ZCHKGG: ',A,' returned INFO=',I6,'.',/9X,'N=',I6,       &
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
99002 FORMAT (' ZCHKGG: ',A,' Eigenvectors from ',A,' incorrectly ',    &
     &        'normalized.',/' Bits of error=',0P,G10.3,',',9X,'N=',I6, &
     &        ', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
99003 FORMAT (1X,A3,' -- Complex Generalized eigenvalue problem')
!
99004 FORMAT (' Matrix types (see ZCHKGG for details): ')
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
99007 FORMAT (/' Tests performed:   (H is Hessenberg, S is Schur, B, ', &
     &        'T, P are triangular,',/20X,'U, V, Q, and Z are ',A,      &
     &        ', l and r are the',/20X,                                 &
     &        'appropriate left and right eigenvectors, resp., a is',   &
     &        /20X,'alpha, b is beta, and ',A,' means ',A,'.)',         &
     &        /' 1 = | A - U H V',A,                                    &
     &        ' | / ( |A| n ulp )      2 = | B - U T V',A,              &
     &        ' | / ( |B| n ulp )',/' 3 = | I - UU',A,                  &
     &        ' | / ( n ulp )             4 = | I - VV',A,              &
     &        ' | / ( n ulp )',/' 5 = | H - Q S Z',A,                   &
     &        ' | / ( |H| n ulp )',6X,'6 = | T - Q P Z',A,              &
     &        ' | / ( |T| n ulp )',/' 7 = | I - QQ',A,                  &
     &        ' | / ( n ulp )             8 = | I - ZZ',A,              &
     &        ' | / ( n ulp )',/' 9 = max | ( b S - a P )',A,           &
     &        ' l | / const.  10 = max | ( b H - a T )',A,              &
     &        ' l | / const.',/                                         &
     &        ' 11= max | ( b S - a P ) r | / const.   12 = max | ( b H'&
     &        ,' - a T ) r | / const.',/1X)
!
99008 FORMAT (' Matrix order=',I5,', type=',I2,', seed=',4(I4,','),     &
     &        ' result ',I2,' is',0P,F8.2)
99009 FORMAT (' Matrix order=',I5,', type=',I2,', seed=',4(I4,','),     &
     &        ' result ',I2,' is',1P,D10.3)
!
!     End of ZCHKGG
!
      END SUBROUTINE ZCHKGG
