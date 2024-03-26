!*==schkbd.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b schkbd
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE schkbd( NSIZES, MVAL, NVAL, NTYPES, DOTYPE, NRHS,
!                          ISEED, THRESH, A, LDA, BD, BE, S1, S2, X, LDX,
!                          Y, Z, Q, LDQ, PT, LDPT, U, VT, WORK, LWORK,
!                          IWORK, NOUT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDPT, LDQ, LDX, LWORK, NOUT, NRHS,
!      $                   NSIZES, NTYPES
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * )
!       REAL               A( LDA, * ), BD( * ), BE( * ), PT( LDPT, * ),
!      $                   Q( LDQ, * ), S1( * ), S2( * ), U( LDPT, * ),
!      $                   VT( LDPT, * ), WORK( * ), X( LDX, * ),
!      $                   Y( LDX, * ), Z( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> schkbd checks the singular value decomposition (SVD) routines.
!>
!> SGEBRD reduces a real general m by n matrix A to upper or lower
!> bidiagonal form B by an orthogonal transformation:  Q' * A * P = B
!> (or A = Q * B * P').  The matrix B is upper bidiagonal if m >= n
!> and lower bidiagonal if m < n.
!>
!> SORGBR generates the orthogonal matrices Q and P' from SGEBRD.
!> Note that Q and P are not necessarily square.
!>
!> SBDSQR computes the singular value decomposition of the bidiagonal
!> matrix B as B = U S V'.  It is called three times to compute
!>    1)  B = U S1 V', where S1 is the diagonal matrix of singular
!>        values and the columns of the matrices U and V are the left
!>        and right singular vectors, respectively, of B.
!>    2)  Same as 1), but the singular values are stored in S2 and the
!>        singular vectors are not computed.
!>    3)  A = (UQ) S (P'V'), the SVD of the original matrix A.
!> In addition, SBDSQR has an option to apply the left orthogonal matrix
!> U to a matrix X, useful in least squares applications.
!>
!> SBDSDC computes the singular value decomposition of the bidiagonal
!> matrix B as B = U S V' using divide-and-conquer. It is called twice
!> to compute
!>    1) B = U S1 V', where S1 is the diagonal matrix of singular
!>        values and the columns of the matrices U and V are the left
!>        and right singular vectors, respectively, of B.
!>    2) Same as 1), but the singular values are stored in S2 and the
!>        singular vectors are not computed.
!>
!>  SBDSVDX computes the singular value decomposition of the bidiagonal
!>  matrix B as B = U S V' using bisection and inverse iteration. It is
!>  called six times to compute
!>     1) B = U S1 V', RANGE='A', where S1 is the diagonal matrix of singular
!>         values and the columns of the matrices U and V are the left
!>         and right singular vectors, respectively, of B.
!>     2) Same as 1), but the singular values are stored in S2 and the
!>         singular vectors are not computed.
!>     3) B = U S1 V', RANGE='I', with where S1 is the diagonal matrix of singular
!>         values and the columns of the matrices U and V are the left
!>         and right singular vectors, respectively, of B
!>     4) Same as 3), but the singular values are stored in S2 and the
!>         singular vectors are not computed.
!>     5) B = U S1 V', RANGE='V', with where S1 is the diagonal matrix of singular
!>         values and the columns of the matrices U and V are the left
!>         and right singular vectors, respectively, of B
!>     6) Same as 5), but the singular values are stored in S2 and the
!>         singular vectors are not computed.
!>
!> For each pair of matrix dimensions (M,N) and each selected matrix
!> type, an M by N matrix A and an M by NRHS matrix X are generated.
!> The problem dimensions are as follows
!>    A:          M x N
!>    Q:          M x min(M,N) (but M x M if NRHS > 0)
!>    P:          min(M,N) x N
!>    B:          min(M,N) x min(M,N)
!>    U, V:       min(M,N) x min(M,N)
!>    S1, S2      diagonal, order min(M,N)
!>    X:          M x NRHS
!>
!> For each generated matrix, 14 tests are performed:
!>
!> Test SGEBRD and SORGBR
!>
!> (1)   | A - Q B PT | / ( |A| max(M,N) ulp ), PT = P'
!>
!> (2)   | I - Q' Q | / ( M ulp )
!>
!> (3)   | I - PT PT' | / ( N ulp )
!>
!> Test SBDSQR on bidiagonal matrix B
!>
!> (4)   | B - U S1 VT | / ( |B| min(M,N) ulp ), VT = V'
!>
!> (5)   | Y - U Z | / ( |Y| max(min(M,N),k) ulp ), where Y = Q' X
!>                                                  and   Z = U' Y.
!> (6)   | I - U' U | / ( min(M,N) ulp )
!>
!> (7)   | I - VT VT' | / ( min(M,N) ulp )
!>
!> (8)   S1 contains min(M,N) nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> (9)   | S1 - S2 | / ( |S1| ulp ), where S2 is computed without
!>                                   computing U and V.
!>
!> (10)  0 if the true singular values of B are within THRESH of
!>       those in S1.  2*THRESH if they are not.  (Tested using
!>       ssvdch)
!>
!> Test SBDSQR on matrix A
!>
!> (11)  | A - (QU) S (VT PT) | / ( |A| max(M,N) ulp )
!>
!> (12)  | X - (QU) Z | / ( |X| max(M,k) ulp )
!>
!> (13)  | I - (QU)'(QU) | / ( M ulp )
!>
!> (14)  | I - (VT PT) (PT'VT') | / ( N ulp )
!>
!> Test SBDSDC on bidiagonal matrix B
!>
!> (15)  | B - U S1 VT | / ( |B| min(M,N) ulp ), VT = V'
!>
!> (16)  | I - U' U | / ( min(M,N) ulp )
!>
!> (17)  | I - VT VT' | / ( min(M,N) ulp )
!>
!> (18)  S1 contains min(M,N) nonnegative values in decreasing order.
!>       (Return 0 if true, 1/ULP if false.)
!>
!> (19)  | S1 - S2 | / ( |S1| ulp ), where S2 is computed without
!>                                   computing U and V.
!>  Test SBDSVDX on bidiagonal matrix B
!>
!>  (20)  | B - U S1 VT | / ( |B| min(M,N) ulp ), VT = V'
!>
!>  (21)  | I - U' U | / ( min(M,N) ulp )
!>
!>  (22)  | I - VT VT' | / ( min(M,N) ulp )
!>
!>  (23)  S1 contains min(M,N) nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!>  (24)  | S1 - S2 | / ( |S1| ulp ), where S2 is computed without
!>                                    computing U and V.
!>
!>  (25)  | S1 - U' B VT' | / ( |S| n ulp )    SBDSVDX('V', 'I')
!>
!>  (26)  | I - U' U | / ( min(M,N) ulp )
!>
!>  (27)  | I - VT VT' | / ( min(M,N) ulp )
!>
!>  (28)  S1 contains min(M,N) nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!>  (29)  | S1 - S2 | / ( |S1| ulp ), where S2 is computed without
!>                                    computing U and V.
!>
!>  (30)  | S1 - U' B VT' | / ( |S1| n ulp )   SBDSVDX('V', 'V')
!>
!>  (31)  | I - U' U | / ( min(M,N) ulp )
!>
!>  (32)  | I - VT VT' | / ( min(M,N) ulp )
!>
!>  (33)  S1 contains min(M,N) nonnegative values in decreasing order.
!>        (Return 0 if true, 1/ULP if false.)
!>
!>  (34)  | S1 - S2 | / ( |S1| ulp ), where S2 is computed without
!>                                    computing U and V.
!>
!> The possible matrix types are
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
!> (6)  Same as (3), but multiplied by SQRT( overflow threshold )
!> (7)  Same as (3), but multiplied by SQRT( underflow threshold )
!>
!> (8)  A matrix of the form  U D V, where U and V are orthogonal and
!>      D has evenly spaced entries 1, ..., ULP with random signs
!>      on the diagonal.
!>
!> (9)  A matrix of the form  U D V, where U and V are orthogonal and
!>      D has geometrically spaced entries 1, ..., ULP with random
!>      signs on the diagonal.
!>
!> (10) A matrix of the form  U D V, where U and V are orthogonal and
!>      D has "clustered" entries 1, ULP,..., ULP with random
!>      signs on the diagonal.
!>
!> (11) Same as (8), but multiplied by SQRT( overflow threshold )
!> (12) Same as (8), but multiplied by SQRT( underflow threshold )
!>
!> (13) Rectangular matrix with random entries chosen from (-1,1).
!> (14) Same as (13), but multiplied by SQRT( overflow threshold )
!> (15) Same as (13), but multiplied by SQRT( underflow threshold )
!>
!> Special case:
!> (16) A bidiagonal matrix with random entries chosen from a
!>      logarithmic distribution on [ulp^2,ulp^(-2)]  (I.e., each
!>      entry is  e^x, where x is chosen uniformly on
!>      [ 2 log(ulp), -2 log(ulp) ] .)  For *this* type:
!>      (a) SGEBRD is not called to reduce it to bidiagonal form.
!>      (b) the bidiagonal is  min(M,N) x min(M,N); if M<N, the
!>          matrix will be lower bidiagonal, otherwise upper.
!>      (c) only tests 5--8 and 14 are performed.
!>
!> A subset of the full set of matrix types may be selected through
!> the logical array DOTYPE.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NSIZES
!> \verbatim
!>          NSIZES is INTEGER
!>          The number of values of M and N contained in the vectors
!>          MVAL and NVAL.  The matrix sizes are used in pairs (M,N).
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NM)
!>          The values of the matrix row dimension M.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NM)
!>          The values of the matrix column dimension N.
!> \endverbatim
!>
!> \param[in] NTYPES
!> \verbatim
!>          NTYPES is INTEGER
!>          The number of elements in DOTYPE.   If it is zero, schkbd
!>          does nothing.  It must be at least zero.  If it is MAXTYP+1
!>          and NSIZES is 1, then an additional type, MAXTYP+1 is
!>          defined, which is to use whatever matrices are in A and B.
!>          This is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
!>          DOTYPE(MAXTYP+1) is .TRUE. .
!> \endverbatim
!>
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          If DOTYPE(j) is .TRUE., then for each size (m,n), a matrix
!>          of type j will be generated.  If NTYPES is smaller than the
!>          maximum number of types defined (PARAMETER MAXTYP), then
!>          types NTYPES+1 through MAXTYP will not be generated.  If
!>          NTYPES is larger than MAXTYP, DOTYPE(MAXTYP+1) through
!>          DOTYPE(NTYPES) will be ignored.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns in the "right-hand side" matrices X, Y,
!>          and Z, used in testing SBDSQR.  If NRHS = 0, then the
!>          operations on the right-hand side will not be tested.
!>          NRHS must be at least 0.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator. The array elements should be between 0 and 4095;
!>          if not they will be reduced mod 4096.  Also, ISEED(4) must
!>          be odd.  The values of ISEED are changed on exit, and can be
!>          used in the next call to schkbd to continue the same random
!>          number sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  To have
!>          every test ratio printed, use THRESH = 0.  Note that the
!>          expected value of the test ratios is O(1), so THRESH should
!>          be a reasonably small multiple of 1, e.g., 10 or 100.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,NMAX)
!>          where NMAX is the maximum value of N in NVAL.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,MMAX),
!>          where MMAX is the maximum value of M in MVAL.
!> \endverbatim
!>
!> \param[out] BD
!> \verbatim
!>          BD is REAL array, dimension
!>                      (max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] BE
!> \verbatim
!>          BE is REAL array, dimension
!>                      (max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] S1
!> \verbatim
!>          S1 is REAL array, dimension
!>                      (max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] S2
!> \verbatim
!>          S2 is REAL array, dimension
!>                      (max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (LDX,NRHS)
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the arrays X, Y, and Z.
!>          LDX >= max(1,MMAX)
!> \endverbatim
!>
!> \param[out] Y
!> \verbatim
!>          Y is REAL array, dimension (LDX,NRHS)
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is REAL array, dimension (LDX,NRHS)
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ,MMAX)
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,MMAX).
!> \endverbatim
!>
!> \param[out] PT
!> \verbatim
!>          PT is REAL array, dimension (LDPT,NMAX)
!> \endverbatim
!>
!> \param[in] LDPT
!> \verbatim
!>          LDPT is INTEGER
!>          The leading dimension of the arrays PT, U, and V.
!>          LDPT >= max(1, max(min(MVAL(j),NVAL(j)))).
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is REAL array, dimension
!>                      (LDPT,max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is REAL array, dimension
!>                      (LDPT,max(min(MVAL(j),NVAL(j))))
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
!>          The number of entries in WORK.  This must be at least
!>          3(M+N) and  M(M + max(M,N,k) + 1) + N*min(M,N)  for all
!>          pairs  (M,N)=(MM(j),NN(j))
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension at least 8*min(M,N)
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The FORTRAN unit number for printing out error messages
!>          (e.g., if a routine returns IINFO not equal to 0.)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          If 0, then everything ran OK.
!>           -1: NSIZES < 0
!>           -2: Some MM(j) < 0
!>           -3: Some NN(j) < 0
!>           -4: NTYPES < 0
!>           -6: NRHS  < 0
!>           -8: THRESH < 0
!>          -11: LDA < 1 or LDA < MMAX, where MMAX is max( MM(j) ).
!>          -17: LDB < 1 or LDB < MMAX.
!>          -21: LDQ < 1 or LDQ < MMAX.
!>          -23: LDPT< 1 or LDPT< MNMAX.
!>          -27: LWORK too small.
!>          If  slatmr, slatms, SGEBRD, SORGBR, or SBDSQR,
!>              returns an error code, the
!>              absolute value of it is returned.
!>
!>-----------------------------------------------------------------------
!>
!>     Some Local Variables and Parameters:
!>     ---- ----- --------- --- ----------
!>
!>     ZERO, ONE       Real 0 and 1.
!>     MAXTYP          The number of types defined.
!>     NTEST           The number of tests performed, or which can
!>                     be performed so far, for the current matrix.
!>     MMAX            Largest value in NN.
!>     NMAX            Largest value in NN.
!>     MNMIN           min(MM(j), NN(j)) (the dimension of the bidiagonal
!>                     matrix.)
!>     MNMAX           The maximum value of MNMIN for j=1,...,NSIZES.
!>     NFAIL           The number of tests which have exceeded THRESH
!>     COND, IMODE     Values to be passed to the matrix generators.
!>     ANORM           Norm of A; passed to matrix generators.
!>
!>     OVFL, UNFL      Overflow and underflow thresholds.
!>     RTOVFL, RTUNFL  Square roots of the previous 2 values.
!>     ULP, ULPINV     Finest relative precision and its inverse.
!>
!>             The following four arrays decode JTYPE:
!>     KTYPE(j)        The general type (1-10) for type "j".
!>     KMODE(j)        The MODE value to be passed to the matrix
!>                     generator for type "j".
!>     KMAGN(j)        The order of magnitude ( O(1),
!>                     O(overflow^(1/2) ), O(underflow^(1/2) )
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
      SUBROUTINE schkbd(Nsizes,Mval,Nval,Ntypes,Dotype,Nrhs,Iseed,      &
     &                  Thresh,A,Lda,Bd,Be,S1,S2,X,Ldx,Y,Z,Q,Ldq,Pt,    &
     &                  Ldpt,U,Vt,Work,Lwork,Iwork,Nout,Info)
      IMPLICIT NONE
!*--SCHKBD496
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldpt , Ldq , Ldx , Lwork , Nout , Nrhs ,     &
     &        Nsizes , Ntypes
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iseed(4) , Iwork(*) , Mval(*) , Nval(*)
      REAL A(Lda,*) , Bd(*) , Be(*) , Pt(Ldpt,*) , Q(Ldq,*) , S1(*) ,   &
     &     S2(*) , U(Ldpt,*) , Vt(Ldpt,*) , Work(*) , X(Ldx,*) ,        &
     &     Y(Ldx,*) , Z(Ldx,*)
!     ..
!
! ======================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO , HALF
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,TWO=2.0E0,HALF=0.5E0)
      INTEGER MAXTYP
      PARAMETER (MAXTYP=16)
!     ..
!     .. Local Scalars ..
      LOGICAL badmm , badnn , bidiag
      CHARACTER uplo
      CHARACTER*3 path
      INTEGER i , iinfo , il , imode , itemp , itype , iu , iwbd ,      &
     &        iwbe , iwbs , iwbz , iwwork , j , jcol , jsize , jtype ,  &
     &        log2ui , m , minwrk , mmax , mnmax , mnmin , mnmin2 , mq ,&
     &        mtypes , n , nfail , nmax , ns1 , ns2 , ntest
      REAL abstol , amninv , anorm , cond , ovfl , rtovfl , rtunfl ,    &
     &     temp1 , temp2 , ulp , ulpinv , unfl , vl , vu
!     ..
!     .. Local Arrays ..
      INTEGER idum(1) , ioldsd(4) , iseed2(4) , kmagn(MAXTYP) ,         &
     &        kmode(MAXTYP) , ktype(MAXTYP)
      REAL dum(1) , dumma(1) , result(40)
!     ..
!     .. External Functions ..
      REAL SLAMCH 
      EXTERNAL SLAMCH 
!     ..
!     .. External Subroutines ..
      EXTERNAL SBDSDC , SBDSQR , SBDSVDX 
     EXTERNAL         SCOPY , SGEBRD , SGEMM , SLABAD 
     EXTERNAL         SLACPY , SLASET , SORGBR 
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , EXP , INT , LOG , MAX , MIN , SQRT
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NUNit
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NUNit , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Data statements ..
      DATA ktype/1 , 2 , 5*4 , 5*6 , 3*9 , 10/
      DATA kmagn/2*1 , 3*1 , 2 , 3 , 3*1 , 2 , 3 , 1 , 2 , 3 , 0/
      DATA kmode/2*0 , 4 , 3 , 1 , 4 , 4 , 4 , 3 , 1 , 4 , 4 , 0 , 0 ,  &
     &     0 , 0/
!     ..
!     .. Executable Statements ..
!
!     Check for errors
!
      Info = 0
!
      badmm = .FALSE.
      badnn = .FALSE.
      mmax = 1
      nmax = 1
      mnmax = 1
      minwrk = 1
      DO j = 1 , Nsizes
         mmax = MAX(mmax,Mval(j))
         IF ( Mval(j)<0 ) badmm = .TRUE.
         nmax = MAX(nmax,Nval(j))
         IF ( Nval(j)<0 ) badnn = .TRUE.
         mnmax = MAX(mnmax,MIN(Mval(j),Nval(j)))
         minwrk = MAX(minwrk,3*(Mval(j)+Nval(j)),Mval(j)                &
     &            *(Mval(j)+MAX(Mval(j),Nval(j),Nrhs)+1)+Nval(j)        &
     &            *MIN(Nval(j),Mval(j)))
      ENDDO
!
!     Check for errors
!
      IF ( Nsizes<0 ) THEN
         Info = -1
      ELSEIF ( badmm ) THEN
         Info = -2
      ELSEIF ( badnn ) THEN
         Info = -3
      ELSEIF ( Ntypes<0 ) THEN
         Info = -4
      ELSEIF ( Nrhs<0 ) THEN
         Info = -6
      ELSEIF ( Lda<mmax ) THEN
         Info = -11
      ELSEIF ( Ldx<mmax ) THEN
         Info = -17
      ELSEIF ( Ldq<mmax ) THEN
         Info = -21
      ELSEIF ( Ldpt<mnmax ) THEN
         Info = -23
      ELSEIF ( minwrk>Lwork ) THEN
         Info = -27
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL xerbla('schkbd',-Info)
         RETURN
      ENDIF
!
!     Initialize constants
!
      path(1:1) = 'Single precision'
      path(2:3) = 'BD'
      nfail = 0
      ntest = 0
      unfl = SLAMCH('Safe minimum')
      ovfl = SLAMCH('Overflow')
      CALL SLABAD(unfl,ovfl)
      ulp = SLAMCH('Precision')
      ulpinv = ONE/ulp
      log2ui = INT(LOG(ulpinv)/LOG(TWO))
      rtunfl = SQRT(unfl)
      rtovfl = SQRT(ovfl)
      INFot = 0
      abstol = 2*unfl
!
!     Loop over sizes, types
!
      DO jsize = 1 , Nsizes
         m = Mval(jsize)
         n = Nval(jsize)
         mnmin = MIN(m,n)
         amninv = ONE/MAX(m,n,1)
!
         IF ( Nsizes/=1 ) THEN
            mtypes = MIN(MAXTYP,Ntypes)
         ELSE
            mtypes = MIN(MAXTYP+1,Ntypes)
         ENDIF
!
         DO jtype = 1 , mtypes
            IF ( .NOT.Dotype(jtype) ) CYCLE
!
            DO j = 1 , 4
               ioldsd(j) = Iseed(j)
            ENDDO
!
            DO j = 1 , 34
               result(j) = -ONE
            ENDDO
!
            uplo = ' '
!
!           Compute "A"
!
!           Control parameters:
!
!           KMAGN  KMODE        KTYPE
!       =1  O(1)   clustered 1  zero
!       =2  large  clustered 2  identity
!       =3  small  exponential  (none)
!       =4         arithmetic   diagonal, (w/ eigenvalues)
!       =5         random       symmetric, w/ eigenvalues
!       =6                      nonsymmetric, w/ singular values
!       =7                      random diagonal
!       =8                      random symmetric
!       =9                      random nonsymmetric
!       =10                     random bidiagonal (log. distrib.)
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
                  anorm = (rtovfl*ulp)*amninv
               ELSEIF ( kmagn(jtype)==3 ) THEN
!
                  anorm = rtunfl*MAX(m,n)*ulpinv
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
               bidiag = .FALSE.
               IF ( itype==1 ) THEN
!
!              Zero matrix
!
                  iinfo = 0
!
               ELSEIF ( itype==2 ) THEN
!
!              Identity
!
                  DO jcol = 1 , mnmin
                     A(jcol,jcol) = anorm
                  ENDDO
!
               ELSEIF ( itype==4 ) THEN
!
!              Diagonal Matrix, [Eigen]values Specified
!
                  CALL slatms(mnmin,mnmin,'S',Iseed,'N',Work,imode,cond,&
     &                        anorm,0,0,'N',A,Lda,Work(mnmin+1),iinfo)
!
               ELSEIF ( itype==5 ) THEN
!
!              Symmetric, eigenvalues specified
!
                  CALL slatms(mnmin,mnmin,'S',Iseed,'S',Work,imode,cond,&
     &                        anorm,m,n,'N',A,Lda,Work(mnmin+1),iinfo)
!
               ELSEIF ( itype==6 ) THEN
!
!              Nonsymmetric, singular values specified
!
                  CALL slatms(m,n,'S',Iseed,'N',Work,imode,cond,anorm,m,&
     &                        n,'N',A,Lda,Work(mnmin+1),iinfo)
!
               ELSEIF ( itype==7 ) THEN
!
!              Diagonal, random entries
!
                  CALL slatmr(mnmin,mnmin,'S',Iseed,'N',Work,6,ONE,ONE, &
     &                        'T','N',Work(mnmin+1),1,ONE,              &
     &                        Work(2*mnmin+1),1,ONE,'N',Iwork,0,0,ZERO, &
     &                        anorm,'NO',A,Lda,Iwork,iinfo)
!
               ELSEIF ( itype==8 ) THEN
!
!              Symmetric, random entries
!
                  CALL slatmr(mnmin,mnmin,'S',Iseed,'S',Work,6,ONE,ONE, &
     &                        'T','N',Work(mnmin+1),1,ONE,              &
     &                        Work(m+mnmin+1),1,ONE,'N',Iwork,m,n,ZERO, &
     &                        anorm,'NO',A,Lda,Iwork,iinfo)
!
               ELSEIF ( itype==9 ) THEN
!
!              Nonsymmetric, random entries
!
                  CALL slatmr(m,n,'S',Iseed,'N',Work,6,ONE,ONE,'T','N', &
     &                        Work(mnmin+1),1,ONE,Work(m+mnmin+1),1,ONE,&
     &                        'N',Iwork,m,n,ZERO,anorm,'NO',A,Lda,Iwork,&
     &                        iinfo)
!
               ELSEIF ( itype==10 ) THEN
!
!              Bidiagonal, random entries
!
                  temp1 = -TWO*LOG(ulp)
                  DO j = 1 , mnmin
                     Bd(j) = EXP(temp1*slarnd(2,Iseed))
                     IF ( j<mnmin ) Be(j) = EXP(temp1*slarnd(2,Iseed))
                  ENDDO
!
                  iinfo = 0
                  bidiag = .TRUE.
                  IF ( m>=n ) THEN
                     uplo = 'U'
                  ELSE
                     uplo = 'L'
                  ENDIF
               ELSE
                  iinfo = 1
               ENDIF
!
               IF ( iinfo==0 ) THEN
!
!              Generate Right-Hand Side
!
                  IF ( bidiag ) THEN
                     CALL slatmr(mnmin,Nrhs,'S',Iseed,'N',Work,6,ONE,   &
     &                           ONE,'T','N',Work(mnmin+1),1,ONE,       &
     &                           Work(2*mnmin+1),1,ONE,'N',Iwork,mnmin, &
     &                           Nrhs,ZERO,ONE,'NO',Y,Ldx,Iwork,iinfo)
                  ELSE
                     CALL slatmr(m,Nrhs,'S',Iseed,'N',Work,6,ONE,ONE,   &
     &                           'T','N',Work(m+1),1,ONE,Work(2*m+1),1, &
     &                           ONE,'N',Iwork,m,Nrhs,ZERO,ONE,'NO',X,  &
     &                           Ldx,Iwork,iinfo)
                  ENDIF
               ENDIF
!
!           Error Exit
!
               IF ( iinfo/=0 ) THEN
                  WRITE (Nout,FMT=99002) 'Generator' , iinfo , m , n ,  &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  RETURN
               ENDIF
            ENDIF
!
!
!           Call SGEBRD and SORGBR to compute B, Q, and P, do tests.
!
            IF ( .NOT.bidiag ) THEN
!
!              Compute transformations to reduce A to bidiagonal form:
!              B := Q' * A * P.
!
               CALL SLACPY(' ',m,n,A,Lda,Q,Ldq)
               CALL SGEBRD(m,n,Q,Ldq,Bd,Be,Work,Work(mnmin+1),          &
     &                     Work(2*mnmin+1),Lwork-2*mnmin,iinfo)
!
!              Check error code from SGEBRD.
!
               IF ( iinfo/=0 ) THEN
                  WRITE (Nout,FMT=99002) 'SGEBRD' , iinfo , m , n ,     &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  RETURN
               ENDIF
!
               CALL SLACPY(' ',m,n,Q,Ldq,Pt,Ldpt)
               IF ( m>=n ) THEN
                  uplo = 'U'
               ELSE
                  uplo = 'L'
               ENDIF
!
!              Generate Q
!
               mq = m
               IF ( Nrhs<=0 ) mq = mnmin
               CALL SORGBR('Q',m,mq,n,Q,Ldq,Work,Work(2*mnmin+1),       &
     &                     Lwork-2*mnmin,iinfo)
!
!              Check error code from SORGBR.
!
               IF ( iinfo/=0 ) THEN
                  WRITE (Nout,FMT=99002) 'SORGBR(Q)' , iinfo , m , n ,  &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  RETURN
               ENDIF
!
!              Generate P'
!
               CALL SORGBR('P',mnmin,n,m,Pt,Ldpt,Work(mnmin+1),         &
     &                     Work(2*mnmin+1),Lwork-2*mnmin,iinfo)
!
!              Check error code from SORGBR.
!
               IF ( iinfo/=0 ) THEN
                  WRITE (Nout,FMT=99002) 'SORGBR(P)' , iinfo , m , n ,  &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  RETURN
               ENDIF
!
!              Apply Q' to an M by NRHS matrix X:  Y := Q' * X.
!
               CALL SGEMM('Transpose','No transpose',m,Nrhs,m,ONE,Q,Ldq,&
     &                    X,Ldx,ZERO,Y,Ldx)
!
!              Test 1:  Check the decomposition A := Q * B * PT
!                   2:  Check the orthogonality of Q
!                   3:  Check the orthogonality of PT
!
               CALL sbdt01(m,n,1,A,Lda,Q,Ldq,Bd,Be,Pt,Ldpt,Work,        &
     &                     result(1))
               CALL sort01('Columns',m,mq,Q,Ldq,Work,Lwork,result(2))
               CALL sort01('Rows',mnmin,n,Pt,Ldpt,Work,Lwork,result(3))
            ENDIF
!
!           Use SBDSQR to form the SVD of the bidiagonal matrix B:
!           B := U * S1 * VT, and compute Z = U' * Y.
!
            CALL SCOPY(mnmin,Bd,1,S1,1)
            IF ( mnmin>0 ) CALL SCOPY(mnmin-1,Be,1,Work,1)
            CALL SLACPY(' ',m,Nrhs,Y,Ldx,Z,Ldx)
            CALL SLASET('Full',mnmin,mnmin,ZERO,ONE,U,Ldpt)
            CALL SLASET('Full',mnmin,mnmin,ZERO,ONE,Vt,Ldpt)
!
            CALL SBDSQR(uplo,mnmin,mnmin,mnmin,Nrhs,S1,Work,Vt,Ldpt,U,  &
     &                  Ldpt,Z,Ldx,Work(mnmin+1),iinfo)
!
!           Check error code from SBDSQR.
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nout,FMT=99002) 'SBDSQR(vects)' , iinfo , m , n , &
     &                                jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  result(4) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Use SBDSQR to compute only the singular values of the
!           bidiagonal matrix B;  U, VT, and Z should not be modified.
!
            CALL SCOPY(mnmin,Bd,1,S2,1)
            IF ( mnmin>0 ) CALL SCOPY(mnmin-1,Be,1,Work,1)
!
            CALL SBDSQR(uplo,mnmin,0,0,0,S2,Work,Vt,Ldpt,U,Ldpt,Z,Ldx,  &
     &                  Work(mnmin+1),iinfo)
!
!           Check error code from SBDSQR.
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nout,FMT=99002) 'SBDSQR(values)' , iinfo , m , n ,&
     &                                jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  result(9) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Test 4:  Check the decomposition B := U * S1 * VT
!                5:  Check the computation Z := U' * Y
!                6:  Check the orthogonality of U
!                7:  Check the orthogonality of VT
!
            CALL sbdt03(uplo,mnmin,1,Bd,Be,U,Ldpt,S1,Vt,Ldpt,Work,      &
     &                  result(4))
            CALL sbdt02(mnmin,Nrhs,Y,Ldx,Z,Ldx,U,Ldpt,Work,result(5))
            CALL sort01('Columns',mnmin,mnmin,U,Ldpt,Work,Lwork,        &
     &                  result(6))
            CALL sort01('Rows',mnmin,mnmin,Vt,Ldpt,Work,Lwork,result(7))
!
!           Test 8:  Check that the singular values are sorted in
!                    non-increasing order and are non-negative
!
            result(8) = ZERO
            DO i = 1 , mnmin - 1
               IF ( S1(i)<S1(i+1) ) result(8) = ulpinv
               IF ( S1(i)<ZERO ) result(8) = ulpinv
            ENDDO
            IF ( mnmin>=1 ) THEN
               IF ( S1(mnmin)<ZERO ) result(8) = ulpinv
            ENDIF
!
!           Test 9:  Compare SBDSQR with and without singular vectors
!
            temp2 = ZERO
!
            DO j = 1 , mnmin
               temp1 = ABS(S1(j)-S2(j))                                 &
     &                 /MAX(SQRT(unfl)*MAX(S1(1),ONE),ulp*MAX(ABS(S1(j))&
     &                 ,ABS(S2(j))))
               temp2 = MAX(temp1,temp2)
            ENDDO
!
            result(9) = temp2
!
!           Test 10:  Sturm sequence test of singular values
!                     Go up by factors of two until it succeeds
!
            temp1 = Thresh*(HALF-ulp)
!
            DO j = 0 , log2ui
!               CALL ssvdch( MNMIN, BD, BE, S1, TEMP1, IINFO )
               IF ( iinfo==0 ) EXIT
               temp1 = temp1*TWO
            ENDDO
!
            result(10) = temp1
!
!           Use SBDSQR to form the decomposition A := (QU) S (VT PT)
!           from the bidiagonal form A := Q B PT.
!
            IF ( .NOT.bidiag ) THEN
               CALL SCOPY(mnmin,Bd,1,S2,1)
               IF ( mnmin>0 ) CALL SCOPY(mnmin-1,Be,1,Work,1)
!
               CALL SBDSQR(uplo,mnmin,n,m,Nrhs,S2,Work,Pt,Ldpt,Q,Ldq,Y, &
     &                     Ldx,Work(mnmin+1),iinfo)
!
!              Test 11:  Check the decomposition A := Q*U * S2 * VT*PT
!                   12:  Check the computation Z := U' * Q' * X
!                   13:  Check the orthogonality of Q*U
!                   14:  Check the orthogonality of VT*PT
!
               CALL sbdt01(m,n,0,A,Lda,Q,Ldq,S2,dumma,Pt,Ldpt,Work,     &
     &                     result(11))
               CALL sbdt02(m,Nrhs,X,Ldx,Y,Ldx,Q,Ldq,Work,result(12))
               CALL sort01('Columns',m,mq,Q,Ldq,Work,Lwork,result(13))
               CALL sort01('Rows',mnmin,n,Pt,Ldpt,Work,Lwork,result(14))
            ENDIF
!
!           Use SBDSDC to form the SVD of the bidiagonal matrix B:
!           B := U * S1 * VT
!
            CALL SCOPY(mnmin,Bd,1,S1,1)
            IF ( mnmin>0 ) CALL SCOPY(mnmin-1,Be,1,Work,1)
            CALL SLASET('Full',mnmin,mnmin,ZERO,ONE,U,Ldpt)
            CALL SLASET('Full',mnmin,mnmin,ZERO,ONE,Vt,Ldpt)
!
            CALL SBDSDC(uplo,'I',mnmin,S1,Work,U,Ldpt,Vt,Ldpt,dum,idum, &
     &                  Work(mnmin+1),Iwork,iinfo)
!
!           Check error code from SBDSDC.
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nout,FMT=99002) 'SBDSDC(vects)' , iinfo , m , n , &
     &                                jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  result(15) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Use SBDSDC to compute only the singular values of the
!           bidiagonal matrix B;  U and VT should not be modified.
!
            CALL SCOPY(mnmin,Bd,1,S2,1)
            IF ( mnmin>0 ) CALL SCOPY(mnmin-1,Be,1,Work,1)
!
            CALL SBDSDC(uplo,'N',mnmin,S2,Work,dum,1,dum,1,dum,idum,    &
     &                  Work(mnmin+1),Iwork,iinfo)
!
!           Check error code from SBDSDC.
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nout,FMT=99002) 'SBDSDC(values)' , iinfo , m , n ,&
     &                                jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  result(18) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Test 15:  Check the decomposition B := U * S1 * VT
!                16:  Check the orthogonality of U
!                17:  Check the orthogonality of VT
!
            CALL sbdt03(uplo,mnmin,1,Bd,Be,U,Ldpt,S1,Vt,Ldpt,Work,      &
     &                  result(15))
            CALL sort01('Columns',mnmin,mnmin,U,Ldpt,Work,Lwork,        &
     &                  result(16))
            CALL sort01('Rows',mnmin,mnmin,Vt,Ldpt,Work,Lwork,result(17)&
     &                  )
!
!           Test 18:  Check that the singular values are sorted in
!                     non-increasing order and are non-negative
!
            result(18) = ZERO
            DO i = 1 , mnmin - 1
               IF ( S1(i)<S1(i+1) ) result(18) = ulpinv
               IF ( S1(i)<ZERO ) result(18) = ulpinv
            ENDDO
            IF ( mnmin>=1 ) THEN
               IF ( S1(mnmin)<ZERO ) result(18) = ulpinv
            ENDIF
!
!           Test 19:  Compare SBDSQR with and without singular vectors
!
            temp2 = ZERO
!
            DO j = 1 , mnmin
               temp1 = ABS(S1(j)-S2(j))                                 &
     &                 /MAX(SQRT(unfl)*MAX(S1(1),ONE),ulp*MAX(ABS(S1(1))&
     &                 ,ABS(S2(1))))
               temp2 = MAX(temp1,temp2)
            ENDDO
!
            result(19) = temp2
!
!
!           Use SBDSVDX to compute the SVD of the bidiagonal matrix B:
!           B := U * S1 * VT
!
            IF ( jtype==10 .OR. jtype==16 ) THEN
!              =================================
!              Matrix types temporarily disabled
!              =================================
               result(20:34) = ZERO
               GOTO 20
            ENDIF
!
            iwbs = 1
            iwbd = iwbs + mnmin
            iwbe = iwbd + mnmin
            iwbz = iwbe + mnmin
            iwwork = iwbz + 2*mnmin*(mnmin+1)
            mnmin2 = MAX(1,mnmin*2)
!
            CALL SCOPY(mnmin,Bd,1,Work(iwbd),1)
            IF ( mnmin>0 ) CALL SCOPY(mnmin-1,Be,1,Work(iwbe),1)
!
            CALL SBDSVDX(uplo,'V','A',mnmin,Work(iwbd),Work(iwbe),ZERO, &
     &                   ZERO,0,0,ns1,S1,Work(iwbz),mnmin2,Work(iwwork),&
     &                   Iwork,iinfo)
!
!           Check error code from SBDSVDX.
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nout,FMT=99002) 'SBDSVDX(vects,A)' , iinfo , m ,  &
     &                                n , jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  result(20) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
            j = iwbz
            DO i = 1 , ns1
               CALL SCOPY(mnmin,Work(j),1,U(1,i),1)
               j = j + mnmin
               CALL SCOPY(mnmin,Work(j),1,Vt(i,1),Ldpt)
               j = j + mnmin
            ENDDO
!
!           Use SBDSVDX to compute only the singular values of the
!           bidiagonal matrix B;  U and VT should not be modified.
!
            IF ( jtype==9 ) THEN
!              =================================
!              Matrix types temporarily disabled
!              =================================
               result(24) = ZERO
               GOTO 20
            ENDIF
!
            CALL SCOPY(mnmin,Bd,1,Work(iwbd),1)
            IF ( mnmin>0 ) CALL SCOPY(mnmin-1,Be,1,Work(iwbe),1)
!
            CALL SBDSVDX(uplo,'N','A',mnmin,Work(iwbd),Work(iwbe),ZERO, &
     &                   ZERO,0,0,ns2,S2,Work(iwbz),mnmin2,Work(iwwork),&
     &                   Iwork,iinfo)
!
!           Check error code from SBDSVDX.
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nout,FMT=99002) 'SBDSVDX(values,A)' , iinfo , m , &
     &                                n , jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  result(24) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Save S1 for tests 30-34.
!
            CALL SCOPY(mnmin,S1,1,Work(iwbs),1)
!
!           Test 20:  Check the decomposition B := U * S1 * VT
!                21:  Check the orthogonality of U
!                22:  Check the orthogonality of VT
!                23:  Check that the singular values are sorted in
!                     non-increasing order and are non-negative
!                24:  Compare SBDSVDX with and without singular vectors
!
            CALL sbdt03(uplo,mnmin,1,Bd,Be,U,Ldpt,S1,Vt,Ldpt,           &
     &                  Work(iwbs+mnmin),result(20))
            CALL sort01('Columns',mnmin,mnmin,U,Ldpt,Work(iwbs+mnmin),  &
     &                  Lwork-mnmin,result(21))
            CALL sort01('Rows',mnmin,mnmin,Vt,Ldpt,Work(iwbs+mnmin),    &
     &                  Lwork-mnmin,result(22))
!
            result(23) = ZERO
            DO i = 1 , mnmin - 1
               IF ( S1(i)<S1(i+1) ) result(23) = ulpinv
               IF ( S1(i)<ZERO ) result(23) = ulpinv
            ENDDO
            IF ( mnmin>=1 ) THEN
               IF ( S1(mnmin)<ZERO ) result(23) = ulpinv
            ENDIF
!
            temp2 = ZERO
            DO j = 1 , mnmin
               temp1 = ABS(S1(j)-S2(j))                                 &
     &                 /MAX(SQRT(unfl)*MAX(S1(1),ONE),ulp*MAX(ABS(S1(1))&
     &                 ,ABS(S2(1))))
               temp2 = MAX(temp1,temp2)
            ENDDO
            result(24) = temp2
            anorm = S1(1)
!
!           Use SBDSVDX with RANGE='I': choose random values for IL and
!           IU, and ask for the IL-th through IU-th singular values
!           and corresponding vectors.
!
            DO i = 1 , 4
               iseed2(i) = Iseed(i)
            ENDDO
            IF ( mnmin<=1 ) THEN
               il = 1
               iu = mnmin
            ELSE
               il = 1 + INT((mnmin-1)*slarnd(1,iseed2))
               iu = 1 + INT((mnmin-1)*slarnd(1,iseed2))
               IF ( iu<il ) THEN
                  itemp = iu
                  iu = il
                  il = itemp
               ENDIF
            ENDIF
!
            CALL SCOPY(mnmin,Bd,1,Work(iwbd),1)
            IF ( mnmin>0 ) CALL SCOPY(mnmin-1,Be,1,Work(iwbe),1)
!
            CALL SBDSVDX(uplo,'V','I',mnmin,Work(iwbd),Work(iwbe),ZERO, &
     &                   ZERO,il,iu,ns1,S1,Work(iwbz),mnmin2,           &
     &                   Work(iwwork),Iwork,iinfo)
!
!           Check error code from SBDSVDX.
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nout,FMT=99002) 'SBDSVDX(vects,I)' , iinfo , m ,  &
     &                                n , jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  result(25) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
            j = iwbz
            DO i = 1 , ns1
               CALL SCOPY(mnmin,Work(j),1,U(1,i),1)
               j = j + mnmin
               CALL SCOPY(mnmin,Work(j),1,Vt(i,1),Ldpt)
               j = j + mnmin
            ENDDO
!
!           Use SBDSVDX to compute only the singular values of the
!           bidiagonal matrix B;  U and VT should not be modified.
!
            CALL SCOPY(mnmin,Bd,1,Work(iwbd),1)
            IF ( mnmin>0 ) CALL SCOPY(mnmin-1,Be,1,Work(iwbe),1)
!
            CALL SBDSVDX(uplo,'N','I',mnmin,Work(iwbd),Work(iwbe),ZERO, &
     &                   ZERO,il,iu,ns2,S2,Work(iwbz),mnmin2,           &
     &                   Work(iwwork),Iwork,iinfo)
!
!           Check error code from SBDSVDX.
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nout,FMT=99002) 'SBDSVDX(values,I)' , iinfo , m , &
     &                                n , jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  result(29) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Test 25:  Check S1 - U' * B * VT'
!                26:  Check the orthogonality of U
!                27:  Check the orthogonality of VT
!                28:  Check that the singular values are sorted in
!                     non-increasing order and are non-negative
!                29:  Compare SBDSVDX with and without singular vectors
!
            CALL sbdt04(uplo,mnmin,Bd,Be,S1,ns1,U,Ldpt,Vt,Ldpt,         &
     &                  Work(iwbs+mnmin),result(25))
            CALL sort01('Columns',mnmin,ns1,U,Ldpt,Work(iwbs+mnmin),    &
     &                  Lwork-mnmin,result(26))
            CALL sort01('Rows',ns1,mnmin,Vt,Ldpt,Work(iwbs+mnmin),      &
     &                  Lwork-mnmin,result(27))
!
            result(28) = ZERO
            DO i = 1 , ns1 - 1
               IF ( S1(i)<S1(i+1) ) result(28) = ulpinv
               IF ( S1(i)<ZERO ) result(28) = ulpinv
            ENDDO
            IF ( ns1>=1 ) THEN
               IF ( S1(ns1)<ZERO ) result(28) = ulpinv
            ENDIF
!
            temp2 = ZERO
            DO j = 1 , ns1
               temp1 = ABS(S1(j)-S2(j))                                 &
     &                 /MAX(SQRT(unfl)*MAX(S1(1),ONE),ulp*MAX(ABS(S1(1))&
     &                 ,ABS(S2(1))))
               temp2 = MAX(temp1,temp2)
            ENDDO
            result(29) = temp2
!
!           Use SBDSVDX with RANGE='V': determine the values VL and VU
!           of the IL-th and IU-th singular values and ask for all
!           singular values in this range.
!
            CALL SCOPY(mnmin,Work(iwbs),1,S1,1)
!
            IF ( mnmin>0 ) THEN
               IF ( il/=1 ) THEN
                  vu = S1(il) + MAX(HALF*ABS(S1(il)-S1(il-1)),ulp*anorm,&
     &                 TWO*rtunfl)
               ELSE
                  vu = S1(1) + MAX(HALF*ABS(S1(mnmin)-S1(1)),ulp*anorm, &
     &                 TWO*rtunfl)
               ENDIF
               IF ( iu/=ns1 ) THEN
                  vl = S1(iu)                                           &
     &                 - MAX(ulp*anorm,TWO*rtunfl,HALF*ABS(S1(iu+1)     &
     &                 -S1(iu)))
               ELSE
                  vl = S1(ns1)                                          &
     &                 - MAX(ulp*anorm,TWO*rtunfl,HALF*ABS(S1(mnmin)    &
     &                 -S1(1)))
               ENDIF
               vl = MAX(vl,ZERO)
               vu = MAX(vu,ZERO)
               IF ( vl>=vu ) vu = MAX(vu*2,vu+vl+HALF)
            ELSE
               vl = ZERO
               vu = ONE
            ENDIF
!
            CALL SCOPY(mnmin,Bd,1,Work(iwbd),1)
            IF ( mnmin>0 ) CALL SCOPY(mnmin-1,Be,1,Work(iwbe),1)
!
            CALL SBDSVDX(uplo,'V','V',mnmin,Work(iwbd),Work(iwbe),vl,vu,&
     &                   0,0,ns1,S1,Work(iwbz),mnmin2,Work(iwwork),     &
     &                   Iwork,iinfo)
!
!           Check error code from SBDSVDX.
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nout,FMT=99002) 'SBDSVDX(vects,V)' , iinfo , m ,  &
     &                                n , jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  result(30) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
            j = iwbz
            DO i = 1 , ns1
               CALL SCOPY(mnmin,Work(j),1,U(1,i),1)
               j = j + mnmin
               CALL SCOPY(mnmin,Work(j),1,Vt(i,1),Ldpt)
               j = j + mnmin
            ENDDO
!
!           Use SBDSVDX to compute only the singular values of the
!           bidiagonal matrix B;  U and VT should not be modified.
!
            CALL SCOPY(mnmin,Bd,1,Work(iwbd),1)
            IF ( mnmin>0 ) CALL SCOPY(mnmin-1,Be,1,Work(iwbe),1)
!
            CALL SBDSVDX(uplo,'N','V',mnmin,Work(iwbd),Work(iwbe),vl,vu,&
     &                   0,0,ns2,S2,Work(iwbz),mnmin2,Work(iwwork),     &
     &                   Iwork,iinfo)
!
!           Check error code from SBDSVDX.
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nout,FMT=99002) 'SBDSVDX(values,V)' , iinfo , m , &
     &                                n , jtype , ioldsd
               Info = ABS(iinfo)
               IF ( iinfo<0 ) THEN
                  RETURN
               ELSE
                  result(34) = ulpinv
                  GOTO 20
               ENDIF
            ENDIF
!
!           Test 30:  Check S1 - U' * B * VT'
!                31:  Check the orthogonality of U
!                32:  Check the orthogonality of VT
!                33:  Check that the singular values are sorted in
!                     non-increasing order and are non-negative
!                34:  Compare SBDSVDX with and without singular vectors
!
            CALL sbdt04(uplo,mnmin,Bd,Be,S1,ns1,U,Ldpt,Vt,Ldpt,         &
     &                  Work(iwbs+mnmin),result(30))
            CALL sort01('Columns',mnmin,ns1,U,Ldpt,Work(iwbs+mnmin),    &
     &                  Lwork-mnmin,result(31))
            CALL sort01('Rows',ns1,mnmin,Vt,Ldpt,Work(iwbs+mnmin),      &
     &                  Lwork-mnmin,result(32))
!
            result(33) = ZERO
            DO i = 1 , ns1 - 1
               IF ( S1(i)<S1(i+1) ) result(28) = ulpinv
               IF ( S1(i)<ZERO ) result(28) = ulpinv
            ENDDO
            IF ( ns1>=1 ) THEN
               IF ( S1(ns1)<ZERO ) result(28) = ulpinv
            ENDIF
!
            temp2 = ZERO
            DO j = 1 , ns1
               temp1 = ABS(S1(j)-S2(j))                                 &
     &                 /MAX(SQRT(unfl)*MAX(S1(1),ONE),ulp*MAX(ABS(S1(1))&
     &                 ,ABS(S2(1))))
               temp2 = MAX(temp1,temp2)
            ENDDO
            result(34) = temp2
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
!
 20         DO j = 1 , 34
               IF ( result(j)>=Thresh ) THEN
                  IF ( nfail==0 ) CALL slahd2(Nout,path)
                  WRITE (Nout,FMT=99001) m , n , jtype , ioldsd , j ,   &
     &                   result(j)
                  nfail = nfail + 1
               ENDIF
            ENDDO
            IF ( .NOT.bidiag ) THEN
               ntest = ntest + 34
            ELSE
               ntest = ntest + 30
            ENDIF
!
         ENDDO
      ENDDO
!
!     Summary
!
      CALL alasum(path,Nout,nfail,ntest,0)
!
      RETURN
!
!     End of schkbd
!
99001 FORMAT (' M=',I5,', N=',I5,', type ',I2,', seed=',4(I4,','),      &
     &        ' test(',I2,')=',G11.4)
99002 FORMAT (' schkbd: ',A,' returned INFO=',I6,'.',/9X,'M=',I6,', N=',&
     &        I6,', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
      END SUBROUTINE schkbd
