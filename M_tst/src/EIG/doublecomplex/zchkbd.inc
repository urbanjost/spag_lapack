!*==zchkbd.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b zchkbd
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zchkbd( NSIZES, MVAL, NVAL, NTYPES, DOTYPE, NRHS,
!                          ISEED, THRESH, A, LDA, BD, BE, S1, S2, X, LDX,
!                          Y, Z, Q, LDQ, PT, LDPT, U, VT, WORK, LWORK,
!                          RWORK, NOUT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDPT, LDQ, LDX, LWORK, NOUT, NRHS,
!      $                   NSIZES, NTYPES
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            ISEED( 4 ), MVAL( * ), NVAL( * )
!       DOUBLE PRECISION   BD( * ), BE( * ), RWORK( * ), S1( * ), S2( * )
!       COMPLEX*16         A( LDA, * ), PT( LDPT, * ), Q( LDQ, * ),
!      $                   U( LDPT, * ), VT( LDPT, * ), WORK( * ),
!      $                   X( LDX, * ), Y( LDX, * ), Z( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zchkbd checks the singular value decomposition (SVD) routines.
!>
!> ZGEBRD reduces a complex general m by n matrix A to real upper or
!> lower bidiagonal form by an orthogonal transformation: Q' * A * P = B
!> (or A = Q * B * P').  The matrix B is upper bidiagonal if m >= n
!> and lower bidiagonal if m < n.
!>
!> ZUNGBR generates the orthogonal matrices Q and P' from ZGEBRD.
!> Note that Q and P are not necessarily square.
!>
!> ZBDSQR computes the singular value decomposition of the bidiagonal
!> matrix B as B = U S V'.  It is called three times to compute
!>    1)  B = U S1 V', where S1 is the diagonal matrix of singular
!>        values and the columns of the matrices U and V are the left
!>        and right singular vectors, respectively, of B.
!>    2)  Same as 1), but the singular values are stored in S2 and the
!>        singular vectors are not computed.
!>    3)  A = (UQ) S (P'V'), the SVD of the original matrix A.
!> In addition, ZBDSQR has an option to apply the left orthogonal matrix
!> U to a matrix X, useful in least squares applications.
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
!> Test ZGEBRD and ZUNGBR
!>
!> (1)   | A - Q B PT | / ( |A| max(M,N) ulp ), PT = P'
!>
!> (2)   | I - Q' Q | / ( M ulp )
!>
!> (3)   | I - PT PT' | / ( N ulp )
!>
!> Test ZBDSQR on bidiagonal matrix B
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
!> (9)   0 if the true singular values of B are within THRESH of
!>       those in S1.  2*THRESH if they are not.  (Tested using
!>       dsvdch)
!>
!> (10)  | S1 - S2 | / ( |S1| ulp ), where S2 is computed without
!>                                   computing U and V.
!>
!> Test ZBDSQR on matrix A
!>
!> (11)  | A - (QU) S (VT PT) | / ( |A| max(M,N) ulp )
!>
!> (12)  | X - (QU) Z | / ( |X| max(M,k) ulp )
!>
!> (13)  | I - (QU)'(QU) | / ( M ulp )
!>
!> (14)  | I - (VT PT) (PT'VT') | / ( N ulp )
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
!>      (a) ZGEBRD is not called to reduce it to bidiagonal form.
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
!>          The number of elements in DOTYPE.   If it is zero, zchkbd
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
!>          and Z, used in testing ZBDSQR.  If NRHS = 0, then the
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
!>          used in the next call to zchkbd to continue the same random
!>          number sequence.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  To have
!>          every test ratio printed, use THRESH = 0.  Note that the
!>          expected value of the test ratios is O(1), so THRESH should
!>          be a reasonably small multiple of 1, e.g., 10 or 100.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,NMAX)
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
!>          BD is DOUBLE PRECISION array, dimension
!>                      (max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] BE
!> \verbatim
!>          BE is DOUBLE PRECISION array, dimension
!>                      (max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] S1
!> \verbatim
!>          S1 is DOUBLE PRECISION array, dimension
!>                      (max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] S2
!> \verbatim
!>          S2 is DOUBLE PRECISION array, dimension
!>                      (max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the arrays X, Y, and Z.
!>          LDX >= max(1,MMAX).
!> \endverbatim
!>
!> \param[out] Y
!> \verbatim
!>          Y is COMPLEX*16 array, dimension (LDX,NRHS)
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDX,NRHS)
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ,MMAX)
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
!>          PT is COMPLEX*16 array, dimension (LDPT,NMAX)
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
!>          U is COMPLEX*16 array, dimension
!>                      (LDPT,max(min(MVAL(j),NVAL(j))))
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is COMPLEX*16 array, dimension
!>                      (LDPT,max(min(MVAL(j),NVAL(j))))
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
!>          3(M+N) and  M(M + max(M,N,k) + 1) + N*min(M,N)  for all
!>          pairs  (M,N)=(MM(j),NN(j))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension
!>                      (5*max(min(M,N)))
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
!>          -23: LDP < 1 or LDP < MNMAX.
!>          -27: LWORK too small.
!>          If  zlatmr, clatms, ZGEBRD, ZUNGBR, or ZBDSQR,
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE zchkbd(Nsizes,Mval,Nval,Ntypes,Dotype,Nrhs,Iseed,      &
     &                  Thresh,A,Lda,Bd,Be,S1,S2,X,Ldx,Y,Z,Q,Ldq,Pt,    &
     &                  Ldpt,U,Vt,Work,Lwork,Rwork,Nout,Info)
      IMPLICIT NONE
!*--ZCHKBD418
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldpt , Ldq , Ldx , Lwork , Nout , Nrhs ,     &
     &        Nsizes , Ntypes
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iseed(4) , Mval(*) , Nval(*)
      DOUBLE PRECISION Bd(*) , Be(*) , Rwork(*) , S1(*) , S2(*)
      COMPLEX*16 A(Lda,*) , Pt(Ldpt,*) , Q(Ldq,*) , U(Ldpt,*) ,         &
     &           Vt(Ldpt,*) , Work(*) , X(Ldx,*) , Y(Ldx,*) , Z(Ldx,*)
!     ..
!
! ======================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TWO , HALF
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
      INTEGER MAXTYP
      PARAMETER (MAXTYP=16)
!     ..
!     .. Local Scalars ..
      LOGICAL badmm , badnn , bidiag
      CHARACTER uplo
      CHARACTER*3 path
      INTEGER i , iinfo , imode , itype , j , jcol , jsize , jtype ,    &
     &        log2ui , m , minwrk , mmax , mnmax , mnmin , mq , mtypes ,&
     &        n , nfail , nmax , ntest
      DOUBLE PRECISION amninv , anorm , cond , ovfl , rtovfl , rtunfl , &
     &                 temp1 , temp2 , ulp , ulpinv , unfl
!     ..
!     .. Local Arrays ..
      INTEGER ioldsd(4) , iwork(1) , kmagn(MAXTYP) , kmode(MAXTYP) ,    &
     &        ktype(MAXTYP)
      DOUBLE PRECISION dumma(1) , result(14)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH 
      EXTERNAL DLAMCH 
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DLABAD 
     EXTERNAL  ZBDSQR , ZGEBRD , ZGEMM 
     EXTERNAL  ZLACPY , ZLASET , ZUNGBR 
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
         CALL xerbla('zchkbd',-Info)
         RETURN
      ENDIF
!
!     Initialize constants
!
      path(1:1) = 'Zomplex precision'
      path(2:3) = 'BD'
      nfail = 0
      ntest = 0
      unfl = DLAMCH('Safe minimum')
      ovfl = DLAMCH('Overflow')
      CALL DLABAD(unfl,ovfl)
      ulp = DLAMCH('Precision')
      ulpinv = ONE/ulp
      log2ui = INT(LOG(ulpinv)/LOG(TWO))
      rtunfl = SQRT(unfl)
      rtovfl = SQRT(ovfl)
      INFot = 0
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
            DO j = 1 , 14
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
               CALL ZLASET('Full',Lda,n,CZERO,CZERO,A,Lda)
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
                  CALL zlatms(mnmin,mnmin,'S',Iseed,'N',Rwork,imode,    &
     &                        cond,anorm,0,0,'N',A,Lda,Work,iinfo)
!
               ELSEIF ( itype==5 ) THEN
!
!              Symmetric, eigenvalues specified
!
                  CALL zlatms(mnmin,mnmin,'S',Iseed,'S',Rwork,imode,    &
     &                        cond,anorm,m,n,'N',A,Lda,Work,iinfo)
!
               ELSEIF ( itype==6 ) THEN
!
!              Nonsymmetric, singular values specified
!
                  CALL zlatms(m,n,'S',Iseed,'N',Rwork,imode,cond,anorm, &
     &                        m,n,'N',A,Lda,Work,iinfo)
!
               ELSEIF ( itype==7 ) THEN
!
!              Diagonal, random entries
!
                  CALL zlatmr(mnmin,mnmin,'S',Iseed,'N',Work,6,ONE,CONE,&
     &                        'T','N',Work(mnmin+1),1,ONE,              &
     &                        Work(2*mnmin+1),1,ONE,'N',iwork,0,0,ZERO, &
     &                        anorm,'NO',A,Lda,iwork,iinfo)
!
               ELSEIF ( itype==8 ) THEN
!
!              Symmetric, random entries
!
                  CALL zlatmr(mnmin,mnmin,'S',Iseed,'S',Work,6,ONE,CONE,&
     &                        'T','N',Work(mnmin+1),1,ONE,              &
     &                        Work(m+mnmin+1),1,ONE,'N',iwork,m,n,ZERO, &
     &                        anorm,'NO',A,Lda,iwork,iinfo)
!
               ELSEIF ( itype==9 ) THEN
!
!              Nonsymmetric, random entries
!
                  CALL zlatmr(m,n,'S',Iseed,'N',Work,6,ONE,CONE,'T','N',&
     &                        Work(mnmin+1),1,ONE,Work(m+mnmin+1),1,ONE,&
     &                        'N',iwork,m,n,ZERO,anorm,'NO',A,Lda,iwork,&
     &                        iinfo)
!
               ELSEIF ( itype==10 ) THEN
!
!              Bidiagonal, random entries
!
                  temp1 = -TWO*LOG(ulp)
                  DO j = 1 , mnmin
                     Bd(j) = EXP(temp1*dlarnd(2,Iseed))
                     IF ( j<mnmin ) Be(j) = EXP(temp1*dlarnd(2,Iseed))
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
                     CALL zlatmr(mnmin,Nrhs,'S',Iseed,'N',Work,6,ONE,   &
     &                           CONE,'T','N',Work(mnmin+1),1,ONE,      &
     &                           Work(2*mnmin+1),1,ONE,'N',iwork,mnmin, &
     &                           Nrhs,ZERO,ONE,'NO',Y,Ldx,iwork,iinfo)
                  ELSE
                     CALL zlatmr(m,Nrhs,'S',Iseed,'N',Work,6,ONE,CONE,  &
     &                           'T','N',Work(m+1),1,ONE,Work(2*m+1),1, &
     &                           ONE,'N',iwork,m,Nrhs,ZERO,ONE,'NO',X,  &
     &                           Ldx,iwork,iinfo)
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
!           Call ZGEBRD and ZUNGBR to compute B, Q, and P, do tests.
!
            IF ( .NOT.bidiag ) THEN
!
!              Compute transformations to reduce A to bidiagonal form:
!              B := Q' * A * P.
!
               CALL ZLACPY(' ',m,n,A,Lda,Q,Ldq)
               CALL ZGEBRD(m,n,Q,Ldq,Bd,Be,Work,Work(mnmin+1),          &
     &                     Work(2*mnmin+1),Lwork-2*mnmin,iinfo)
!
!              Check error code from ZGEBRD.
!
               IF ( iinfo/=0 ) THEN
                  WRITE (Nout,FMT=99002) 'ZGEBRD' , iinfo , m , n ,     &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  RETURN
               ENDIF
!
               CALL ZLACPY(' ',m,n,Q,Ldq,Pt,Ldpt)
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
               CALL ZUNGBR('Q',m,mq,n,Q,Ldq,Work,Work(2*mnmin+1),       &
     &                     Lwork-2*mnmin,iinfo)
!
!              Check error code from ZUNGBR.
!
               IF ( iinfo/=0 ) THEN
                  WRITE (Nout,FMT=99002) 'ZUNGBR(Q)' , iinfo , m , n ,  &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  RETURN
               ENDIF
!
!              Generate P'
!
               CALL ZUNGBR('P',mnmin,n,m,Pt,Ldpt,Work(mnmin+1),         &
     &                     Work(2*mnmin+1),Lwork-2*mnmin,iinfo)
!
!              Check error code from ZUNGBR.
!
               IF ( iinfo/=0 ) THEN
                  WRITE (Nout,FMT=99002) 'ZUNGBR(P)' , iinfo , m , n ,  &
     &                   jtype , ioldsd
                  Info = ABS(iinfo)
                  RETURN
               ENDIF
!
!              Apply Q' to an M by NRHS matrix X:  Y := Q' * X.
!
               CALL ZGEMM('Conjugate transpose','No transpose',m,Nrhs,m,&
     &                    CONE,Q,Ldq,X,Ldx,CZERO,Y,Ldx)
!
!              Test 1:  Check the decomposition A := Q * B * PT
!                   2:  Check the orthogonality of Q
!                   3:  Check the orthogonality of PT
!
               CALL zbdt01(m,n,1,A,Lda,Q,Ldq,Bd,Be,Pt,Ldpt,Work,Rwork,  &
     &                     result(1))
               CALL zunt01('Columns',m,mq,Q,Ldq,Work,Lwork,Rwork,       &
     &                     result(2))
               CALL zunt01('Rows',mnmin,n,Pt,Ldpt,Work,Lwork,Rwork,     &
     &                     result(3))
            ENDIF
!
!           Use ZBDSQR to form the SVD of the bidiagonal matrix B:
!           B := U * S1 * VT, and compute Z = U' * Y.
!
            CALL DCOPY(mnmin,Bd,1,S1,1)
            IF ( mnmin>0 ) CALL DCOPY(mnmin-1,Be,1,Rwork,1)
            CALL ZLACPY(' ',m,Nrhs,Y,Ldx,Z,Ldx)
            CALL ZLASET('Full',mnmin,mnmin,CZERO,CONE,U,Ldpt)
            CALL ZLASET('Full',mnmin,mnmin,CZERO,CONE,Vt,Ldpt)
!
            CALL ZBDSQR(uplo,mnmin,mnmin,mnmin,Nrhs,S1,Rwork,Vt,Ldpt,U, &
     &                  Ldpt,Z,Ldx,Rwork(mnmin+1),iinfo)
!
!           Check error code from ZBDSQR.
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nout,FMT=99002) 'ZBDSQR(vects)' , iinfo , m , n , &
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
!           Use ZBDSQR to compute only the singular values of the
!           bidiagonal matrix B;  U, VT, and Z should not be modified.
!
            CALL DCOPY(mnmin,Bd,1,S2,1)
            IF ( mnmin>0 ) CALL DCOPY(mnmin-1,Be,1,Rwork,1)
!
            CALL ZBDSQR(uplo,mnmin,0,0,0,S2,Rwork,Vt,Ldpt,U,Ldpt,Z,Ldx, &
     &                  Rwork(mnmin+1),iinfo)
!
!           Check error code from ZBDSQR.
!
            IF ( iinfo/=0 ) THEN
               WRITE (Nout,FMT=99002) 'ZBDSQR(values)' , iinfo , m , n ,&
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
            CALL zbdt03(uplo,mnmin,1,Bd,Be,U,Ldpt,S1,Vt,Ldpt,Work,      &
     &                  result(4))
            CALL zbdt02(mnmin,Nrhs,Y,Ldx,Z,Ldx,U,Ldpt,Work,Rwork,       &
     &                  result(5))
            CALL zunt01('Columns',mnmin,mnmin,U,Ldpt,Work,Lwork,Rwork,  &
     &                  result(6))
            CALL zunt01('Rows',mnmin,mnmin,Vt,Ldpt,Work,Lwork,Rwork,    &
     &                  result(7))
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
!           Test 9:  Compare ZBDSQR with and without singular vectors
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
               CALL dsvdch(mnmin,Bd,Be,S1,temp1,iinfo)
               IF ( iinfo==0 ) EXIT
               temp1 = temp1*TWO
            ENDDO
!
            result(10) = temp1
!
!           Use ZBDSQR to form the decomposition A := (QU) S (VT PT)
!           from the bidiagonal form A := Q B PT.
!
            IF ( .NOT.bidiag ) THEN
               CALL DCOPY(mnmin,Bd,1,S2,1)
               IF ( mnmin>0 ) CALL DCOPY(mnmin-1,Be,1,Rwork,1)
!
               CALL ZBDSQR(uplo,mnmin,n,m,Nrhs,S2,Rwork,Pt,Ldpt,Q,Ldq,Y,&
     &                     Ldx,Rwork(mnmin+1),iinfo)
!
!              Test 11:  Check the decomposition A := Q*U * S2 * VT*PT
!                   12:  Check the computation Z := U' * Q' * X
!                   13:  Check the orthogonality of Q*U
!                   14:  Check the orthogonality of VT*PT
!
               CALL zbdt01(m,n,0,A,Lda,Q,Ldq,S2,dumma,Pt,Ldpt,Work,     &
     &                     Rwork,result(11))
               CALL zbdt02(m,Nrhs,X,Ldx,Y,Ldx,Q,Ldq,Work,Rwork,         &
     &                     result(12))
               CALL zunt01('Columns',m,mq,Q,Ldq,Work,Lwork,Rwork,       &
     &                     result(13))
               CALL zunt01('Rows',mnmin,n,Pt,Ldpt,Work,Lwork,Rwork,     &
     &                     result(14))
            ENDIF
!
!           End of Loop -- Check for RESULT(j) > THRESH
!
 20         DO j = 1 , 14
               IF ( result(j)>=Thresh ) THEN
                  IF ( nfail==0 ) CALL dlahd2(Nout,path)
                  WRITE (Nout,FMT=99001) m , n , jtype , ioldsd , j ,   &
     &                   result(j)
                  nfail = nfail + 1
               ENDIF
            ENDDO
            IF ( .NOT.bidiag ) THEN
               ntest = ntest + 14
            ELSE
               ntest = ntest + 5
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
!     End of zchkbd
!
99001 FORMAT (' M=',I5,', N=',I5,', type ',I2,', seed=',4(I4,','),      &
     &        ' test(',I2,')=',G11.4)
99002 FORMAT (' zchkbd: ',A,' returned INFO=',I6,'.',/9X,'M=',I6,', N=',&
     &        I6,', JTYPE=',I6,', ISEED=(',3(I5,','),I5,')')
!
      END SUBROUTINE zchkbd
