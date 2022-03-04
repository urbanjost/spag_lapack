!*==slatmr.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SLATMR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLATMR( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX,
!                          RSIGN, GRADE, DL, MODEL, CONDL, DR, MODER,
!                          CONDR, PIVTNG, IPIVOT, KL, KU, SPARSE, ANORM,
!                          PACK, A, LDA, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIST, GRADE, PACK, PIVTNG, RSIGN, SYM
!       INTEGER            INFO, KL, KU, LDA, M, MODE, MODEL, MODER, N
!       REAL               ANORM, COND, CONDL, CONDR, DMAX, SPARSE
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIVOT( * ), ISEED( 4 ), IWORK( * )
!       REAL               A( LDA, * ), D( * ), DL( * ), DR( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SLATMR generates random matrices of various types for testing
!>    LAPACK programs.
!>
!>    SLATMR operates by applying the following sequence of
!>    operations:
!>
!>      Generate a matrix A with random entries of distribution DIST
!>         which is symmetric if SYM='S', and nonsymmetric
!>         if SYM='N'.
!>
!>      Set the diagonal to D, where D may be input or
!>         computed according to MODE, COND, DMAX and RSIGN
!>         as described below.
!>
!>      Grade the matrix, if desired, from the left and/or right
!>         as specified by GRADE. The inputs DL, MODEL, CONDL, DR,
!>         MODER and CONDR also determine the grading as described
!>         below.
!>
!>      Permute, if desired, the rows and/or columns as specified by
!>         PIVTNG and IPIVOT.
!>
!>      Set random entries to zero, if desired, to get a random sparse
!>         matrix as specified by SPARSE.
!>
!>      Make A a band matrix, if desired, by zeroing out the matrix
!>         outside a band of lower bandwidth KL and upper bandwidth KU.
!>
!>      Scale A, if desired, to have maximum entry ANORM.
!>
!>      Pack the matrix if desired. Options specified by PACK are:
!>         no packing
!>         zero out upper half (if symmetric)
!>         zero out lower half (if symmetric)
!>         store the upper half columnwise (if symmetric or
!>             square upper triangular)
!>         store the lower half columnwise (if symmetric or
!>             square lower triangular)
!>             same as upper half rowwise if symmetric
!>         store the lower triangle in banded format (if symmetric)
!>         store the upper triangle in banded format (if symmetric)
!>         store the entire matrix in banded format
!>
!>    Note: If two calls to SLATMR differ only in the PACK parameter,
!>          they will generate mathematically equivalent matrices.
!>
!>          If two calls to SLATMR both have full bandwidth (KL = M-1
!>          and KU = N-1), and differ only in the PIVTNG and PACK
!>          parameters, then the matrices generated will differ only
!>          in the order of the rows and/or columns, and otherwise
!>          contain the same data. This consistency cannot be and
!>          is not maintained with less than full bandwidth.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           Number of rows of A. Not modified.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           Number of columns of A. Not modified.
!> \endverbatim
!>
!> \param[in] DIST
!> \verbatim
!>          DIST is CHARACTER*1
!>           On entry, DIST specifies the type of distribution to be used
!>           to generate a random matrix .
!>           'U' => UNIFORM( 0, 1 )  ( 'U' for uniform )
!>           'S' => UNIFORM( -1, 1 ) ( 'S' for symmetric )
!>           'N' => NORMAL( 0, 1 )   ( 'N' for normal )
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>           On entry ISEED specifies the seed of the random number
!>           generator. They should lie between 0 and 4095 inclusive,
!>           and ISEED(4) should be odd. The random number generator
!>           uses a linear congruential sequence limited to small
!>           integers, and so should produce machine independent
!>           random numbers. The values of ISEED are changed on
!>           exit, and can be used in the next call to SLATMR
!>           to continue the same random number sequence.
!>           Changed on exit.
!> \endverbatim
!>
!> \param[in] SYM
!> \verbatim
!>          SYM is CHARACTER*1
!>           If SYM='S' or 'H', generated matrix is symmetric.
!>           If SYM='N', generated matrix is nonsymmetric.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (min(M,N))
!>           On entry this array specifies the diagonal entries
!>           of the diagonal of A.  D may either be specified
!>           on entry, or set according to MODE and COND as described
!>           below. May be changed on exit if MODE is nonzero.
!> \endverbatim
!>
!> \param[in] MODE
!> \verbatim
!>          MODE is INTEGER
!>           On entry describes how D is to be used:
!>           MODE = 0 means use D as input
!>           MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND
!>           MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND
!>           MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))
!>           MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)
!>           MODE = 5 sets D to random numbers in the range
!>                    ( 1/COND , 1 ) such that their logarithms
!>                    are uniformly distributed.
!>           MODE = 6 set D to random numbers from same distribution
!>                    as the rest of the matrix.
!>           MODE < 0 has the same meaning as ABS(MODE), except that
!>              the order of the elements of D is reversed.
!>           Thus if MODE is positive, D has entries ranging from
!>              1 to 1/COND, if negative, from 1/COND to 1,
!>           Not modified.
!> \endverbatim
!>
!> \param[in] COND
!> \verbatim
!>          COND is REAL
!>           On entry, used as described under MODE above.
!>           If used, it must be >= 1. Not modified.
!> \endverbatim
!>
!> \param[in] DMAX
!> \verbatim
!>          DMAX is REAL
!>           If MODE neither -6, 0 nor 6, the diagonal is scaled by
!>           DMAX / max(abs(D(i))), so that maximum absolute entry
!>           of diagonal is abs(DMAX). If DMAX is negative (or zero),
!>           diagonal will be scaled by a negative number (or zero).
!> \endverbatim
!>
!> \param[in] RSIGN
!> \verbatim
!>          RSIGN is CHARACTER*1
!>           If MODE neither -6, 0 nor 6, specifies sign of diagonal
!>           as follows:
!>           'T' => diagonal entries are multiplied by 1 or -1
!>                  with probability .5
!>           'F' => diagonal unchanged
!>           Not modified.
!> \endverbatim
!>
!> \param[in] GRADE
!> \verbatim
!>          GRADE is CHARACTER*1
!>           Specifies grading of matrix as follows:
!>           'N'  => no grading
!>           'L'  => matrix premultiplied by diag( DL )
!>                   (only if matrix nonsymmetric)
!>           'R'  => matrix postmultiplied by diag( DR )
!>                   (only if matrix nonsymmetric)
!>           'B'  => matrix premultiplied by diag( DL ) and
!>                         postmultiplied by diag( DR )
!>                   (only if matrix nonsymmetric)
!>           'S' or 'H'  => matrix premultiplied by diag( DL ) and
!>                          postmultiplied by diag( DL )
!>                          ('S' for symmetric, or 'H' for Hermitian)
!>           'E'  => matrix premultiplied by diag( DL ) and
!>                         postmultiplied by inv( diag( DL ) )
!>                         ( 'E' for eigenvalue invariance)
!>                   (only if matrix nonsymmetric)
!>                   Note: if GRADE='E', then M must equal N.
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] DL
!> \verbatim
!>          DL is REAL array, dimension (M)
!>           If MODEL=0, then on entry this array specifies the diagonal
!>           entries of a diagonal matrix used as described under GRADE
!>           above. If MODEL is not zero, then DL will be set according
!>           to MODEL and CONDL, analogous to the way D is set according
!>           to MODE and COND (except there is no DMAX parameter for DL).
!>           If GRADE='E', then DL cannot have zero entries.
!>           Not referenced if GRADE = 'N' or 'R'. Changed on exit.
!> \endverbatim
!>
!> \param[in] MODEL
!> \verbatim
!>          MODEL is INTEGER
!>           This specifies how the diagonal array DL is to be computed,
!>           just as MODE specifies how D is to be computed.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] CONDL
!> \verbatim
!>          CONDL is REAL
!>           When MODEL is not zero, this specifies the condition number
!>           of the computed DL.  Not modified.
!> \endverbatim
!>
!> \param[in,out] DR
!> \verbatim
!>          DR is REAL array, dimension (N)
!>           If MODER=0, then on entry this array specifies the diagonal
!>           entries of a diagonal matrix used as described under GRADE
!>           above. If MODER is not zero, then DR will be set according
!>           to MODER and CONDR, analogous to the way D is set according
!>           to MODE and COND (except there is no DMAX parameter for DR).
!>           Not referenced if GRADE = 'N', 'L', 'H', 'S' or 'E'.
!>           Changed on exit.
!> \endverbatim
!>
!> \param[in] MODER
!> \verbatim
!>          MODER is INTEGER
!>           This specifies how the diagonal array DR is to be computed,
!>           just as MODE specifies how D is to be computed.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] CONDR
!> \verbatim
!>          CONDR is REAL
!>           When MODER is not zero, this specifies the condition number
!>           of the computed DR.  Not modified.
!> \endverbatim
!>
!> \param[in] PIVTNG
!> \verbatim
!>          PIVTNG is CHARACTER*1
!>           On entry specifies pivoting permutations as follows:
!>           'N' or ' ' => none.
!>           'L' => left or row pivoting (matrix must be nonsymmetric).
!>           'R' => right or column pivoting (matrix must be
!>                  nonsymmetric).
!>           'B' or 'F' => both or full pivoting, i.e., on both sides.
!>                         In this case, M must equal N
!>
!>           If two calls to SLATMR both have full bandwidth (KL = M-1
!>           and KU = N-1), and differ only in the PIVTNG and PACK
!>           parameters, then the matrices generated will differ only
!>           in the order of the rows and/or columns, and otherwise
!>           contain the same data. This consistency cannot be
!>           maintained with less than full bandwidth.
!> \endverbatim
!>
!> \param[in] IPIVOT
!> \verbatim
!>          IPIVOT is INTEGER array, dimension (N or M)
!>           This array specifies the permutation used.  After the
!>           basic matrix is generated, the rows, columns, or both
!>           are permuted.   If, say, row pivoting is selected, SLATMR
!>           starts with the *last* row and interchanges the M-th and
!>           IPIVOT(M)-th rows, then moves to the next-to-last row,
!>           interchanging the (M-1)-th and the IPIVOT(M-1)-th rows,
!>           and so on.  In terms of "2-cycles", the permutation is
!>           (1 IPIVOT(1)) (2 IPIVOT(2)) ... (M IPIVOT(M))
!>           where the rightmost cycle is applied first.  This is the
!>           *inverse* of the effect of pivoting in LINPACK.  The idea
!>           is that factoring (with pivoting) an identity matrix
!>           which has been inverse-pivoted in this way should
!>           result in a pivot vector identical to IPIVOT.
!>           Not referenced if PIVTNG = 'N'. Not modified.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>           On entry specifies the lower bandwidth of the  matrix. For
!>           example, KL=0 implies upper triangular, KL=1 implies upper
!>           Hessenberg, and KL at least M-1 implies the matrix is not
!>           banded. Must equal KU if matrix is symmetric.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>           On entry specifies the upper bandwidth of the  matrix. For
!>           example, KU=0 implies lower triangular, KU=1 implies lower
!>           Hessenberg, and KU at least N-1 implies the matrix is not
!>           banded. Must equal KL if matrix is symmetric.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] SPARSE
!> \verbatim
!>          SPARSE is REAL
!>           On entry specifies the sparsity of the matrix if a sparse
!>           matrix is to be generated. SPARSE should lie between
!>           0 and 1. To generate a sparse matrix, for each matrix entry
!>           a uniform ( 0, 1 ) random number x is generated and
!>           compared to SPARSE; if x is larger the matrix entry
!>           is unchanged and if x is smaller the entry is set
!>           to zero. Thus on the average a fraction SPARSE of the
!>           entries will be set to zero.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is REAL
!>           On entry specifies maximum entry of output matrix
!>           (output matrix will by multiplied by a constant so that
!>           its largest absolute entry equal ANORM)
!>           if ANORM is nonnegative. If ANORM is negative no scaling
!>           is done. Not modified.
!> \endverbatim
!>
!> \param[in] PACK
!> \verbatim
!>          PACK is CHARACTER*1
!>           On entry specifies packing of matrix as follows:
!>           'N' => no packing
!>           'U' => zero out all subdiagonal entries (if symmetric)
!>           'L' => zero out all superdiagonal entries (if symmetric)
!>           'C' => store the upper triangle columnwise
!>                  (only if matrix symmetric or square upper triangular)
!>           'R' => store the lower triangle columnwise
!>                  (only if matrix symmetric or square lower triangular)
!>                  (same as upper half rowwise if symmetric)
!>           'B' => store the lower triangle in band storage scheme
!>                  (only if matrix symmetric)
!>           'Q' => store the upper triangle in band storage scheme
!>                  (only if matrix symmetric)
!>           'Z' => store the entire matrix in band storage scheme
!>                      (pivoting can be provided for by using this
!>                      option to store A in the trailing rows of
!>                      the allocated storage)
!>
!>           Using these options, the various LAPACK packed and banded
!>           storage schemes can be obtained:
!>           GB               - use 'Z'
!>           PB, SB or TB     - use 'B' or 'Q'
!>           PP, SP or TP     - use 'C' or 'R'
!>
!>           If two calls to SLATMR differ only in the PACK parameter,
!>           they will generate mathematically equivalent matrices.
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>           On exit A is the desired test matrix. Only those
!>           entries of A which are significant on output
!>           will be referenced (even if A is in packed or band
!>           storage format). The 'unoccupied corners' of A in
!>           band format will be zeroed out.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           on entry LDA specifies the first dimension of A as
!>           declared in the calling program.
!>           If PACK='N', 'U' or 'L', LDA must be at least max ( 1, M ).
!>           If PACK='C' or 'R', LDA must be at least 1.
!>           If PACK='B', or 'Q', LDA must be MIN ( KU+1, N )
!>           If PACK='Z', LDA must be at least KUU+KLL+1, where
!>           KUU = MIN ( KU, N-1 ) and KLL = MIN ( KL, M-1 )
!>           Not modified.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension ( N or M)
!>           Workspace. Not referenced if PIVTNG = 'N'. Changed on exit.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           Error parameter on exit:
!>             0 => normal return
!>            -1 => M negative or unequal to N and SYM='S' or 'H'
!>            -2 => N negative
!>            -3 => DIST illegal string
!>            -5 => SYM illegal string
!>            -7 => MODE not in range -6 to 6
!>            -8 => COND less than 1.0, and MODE neither -6, 0 nor 6
!>           -10 => MODE neither -6, 0 nor 6 and RSIGN illegal string
!>           -11 => GRADE illegal string, or GRADE='E' and
!>                  M not equal to N, or GRADE='L', 'R', 'B' or 'E' and
!>                  SYM = 'S' or 'H'
!>           -12 => GRADE = 'E' and DL contains zero
!>           -13 => MODEL not in range -6 to 6 and GRADE= 'L', 'B', 'H',
!>                  'S' or 'E'
!>           -14 => CONDL less than 1.0, GRADE='L', 'B', 'H', 'S' or 'E',
!>                  and MODEL neither -6, 0 nor 6
!>           -16 => MODER not in range -6 to 6 and GRADE= 'R' or 'B'
!>           -17 => CONDR less than 1.0, GRADE='R' or 'B', and
!>                  MODER neither -6, 0 nor 6
!>           -18 => PIVTNG illegal string, or PIVTNG='B' or 'F' and
!>                  M not equal to N, or PIVTNG='L' or 'R' and SYM='S'
!>                  or 'H'
!>           -19 => IPIVOT contains out of range number and
!>                  PIVTNG not equal to 'N'
!>           -20 => KL negative
!>           -21 => KU negative, or SYM='S' or 'H' and KU not equal to KL
!>           -22 => SPARSE not in range 0. to 1.
!>           -24 => PACK illegal string, or PACK='U', 'L', 'B' or 'Q'
!>                  and SYM='N', or PACK='C' and SYM='N' and either KL
!>                  not equal to 0 or N not equal to M, or PACK='R' and
!>                  SYM='N', and either KU not equal to 0 or N not equal
!>                  to M
!>           -26 => LDA too small
!>             1 => Error return from SLATM1 (computing D)
!>             2 => Cannot scale diagonal to DMAX (max. entry is 0)
!>             3 => Error return from SLATM1 (computing DL)
!>             4 => Error return from SLATM1 (computing DR)
!>             5 => ANORM is positive, but matrix constructed prior to
!>                  attempting to scale it to have norm ANORM, is zero
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
!> \ingroup real_matgen
!
!  =====================================================================
      SUBROUTINE SLATMR(M,N,Dist,Iseed,Sym,D,Mode,Cond,Dmax,Rsign,Grade,&
     &                  Dl,Model,Condl,Dr,Moder,Condr,Pivtng,Ipivot,Kl, &
     &                  Ku,Sparse,Anorm,Pack,A,Lda,Iwork,Info)
      IMPLICIT NONE
!*--SLATMR474
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Dist , Grade , Pack , Pivtng , Rsign , Sym
      INTEGER Info , Kl , Ku , Lda , M , Mode , Model , Moder , N
      REAL Anorm , Cond , Condl , Condr , Dmax , Sparse
!     ..
!     .. Array Arguments ..
      INTEGER Ipivot(*) , Iseed(4) , Iwork(*)
      REAL A(Lda,*) , D(*) , Dl(*) , Dr(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
      REAL ONE
      PARAMETER (ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      LOGICAL badpvt , dzero , fulbnd
      INTEGER i , idist , igrade , iisub , ipack , ipvtng , irsign ,    &
     &        isub , isym , j , jjsub , jsub , k , kll , kuu , mnmin ,  &
     &        mnsub , mxsub , npvts
      REAL alpha , onorm , temp
!     ..
!     .. Local Arrays ..
      REAL tempa(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLANGB , SLANGE , SLANSB , SLANSP , SLANSY , SLATM2 , SLATM3
      EXTERNAL LSAME , SLANGB , SLANGE , SLANSB , SLANSP , SLANSY ,     &
     &         SLATM2 , SLATM3
!     ..
!     .. External Subroutines ..
      EXTERNAL SLATM1 , SSCAL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , MOD
!     ..
!     .. Executable Statements ..
!
!     1)      Decode and Test the input parameters.
!             Initialize flags & seed.
!
      Info = 0
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) RETURN
!
!     Decode DIST
!
      IF ( LSAME(Dist,'U') ) THEN
         idist = 1
      ELSEIF ( LSAME(Dist,'S') ) THEN
         idist = 2
      ELSEIF ( LSAME(Dist,'N') ) THEN
         idist = 3
      ELSE
         idist = -1
      ENDIF
!
!     Decode SYM
!
      IF ( LSAME(Sym,'S') ) THEN
         isym = 0
      ELSEIF ( LSAME(Sym,'N') ) THEN
         isym = 1
      ELSEIF ( LSAME(Sym,'H') ) THEN
         isym = 0
      ELSE
         isym = -1
      ENDIF
!
!     Decode RSIGN
!
      IF ( LSAME(Rsign,'F') ) THEN
         irsign = 0
      ELSEIF ( LSAME(Rsign,'T') ) THEN
         irsign = 1
      ELSE
         irsign = -1
      ENDIF
!
!     Decode PIVTNG
!
      IF ( LSAME(Pivtng,'N') ) THEN
         ipvtng = 0
      ELSEIF ( LSAME(Pivtng,' ') ) THEN
         ipvtng = 0
      ELSEIF ( LSAME(Pivtng,'L') ) THEN
         ipvtng = 1
         npvts = M
      ELSEIF ( LSAME(Pivtng,'R') ) THEN
         ipvtng = 2
         npvts = N
      ELSEIF ( LSAME(Pivtng,'B') ) THEN
         ipvtng = 3
         npvts = MIN(N,M)
      ELSEIF ( LSAME(Pivtng,'F') ) THEN
         ipvtng = 3
         npvts = MIN(N,M)
      ELSE
         ipvtng = -1
      ENDIF
!
!     Decode GRADE
!
      IF ( LSAME(Grade,'N') ) THEN
         igrade = 0
      ELSEIF ( LSAME(Grade,'L') ) THEN
         igrade = 1
      ELSEIF ( LSAME(Grade,'R') ) THEN
         igrade = 2
      ELSEIF ( LSAME(Grade,'B') ) THEN
         igrade = 3
      ELSEIF ( LSAME(Grade,'E') ) THEN
         igrade = 4
      ELSEIF ( LSAME(Grade,'H') .OR. LSAME(Grade,'S') ) THEN
         igrade = 5
      ELSE
         igrade = -1
      ENDIF
!
!     Decode PACK
!
      IF ( LSAME(Pack,'N') ) THEN
         ipack = 0
      ELSEIF ( LSAME(Pack,'U') ) THEN
         ipack = 1
      ELSEIF ( LSAME(Pack,'L') ) THEN
         ipack = 2
      ELSEIF ( LSAME(Pack,'C') ) THEN
         ipack = 3
      ELSEIF ( LSAME(Pack,'R') ) THEN
         ipack = 4
      ELSEIF ( LSAME(Pack,'B') ) THEN
         ipack = 5
      ELSEIF ( LSAME(Pack,'Q') ) THEN
         ipack = 6
      ELSEIF ( LSAME(Pack,'Z') ) THEN
         ipack = 7
      ELSE
         ipack = -1
      ENDIF
!
!     Set certain internal parameters
!
      mnmin = MIN(M,N)
      kll = MIN(Kl,M-1)
      kuu = MIN(Ku,N-1)
!
!     If inv(DL) is used, check to see if DL has a zero entry.
!
      dzero = .FALSE.
      IF ( igrade==4 .AND. Model==0 ) THEN
         DO i = 1 , M
            IF ( Dl(i)==ZERO ) dzero = .TRUE.
         ENDDO
      ENDIF
!
!     Check values in IPIVOT
!
      badpvt = .FALSE.
      IF ( ipvtng>0 ) THEN
         DO j = 1 , npvts
            IF ( Ipivot(j)<=0 .OR. Ipivot(j)>npvts ) badpvt = .TRUE.
         ENDDO
      ENDIF
!
!     Set INFO if an error
!
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( M/=N .AND. isym==0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( idist==-1 ) THEN
         Info = -3
      ELSEIF ( isym==-1 ) THEN
         Info = -5
      ELSEIF ( Mode<-6 .OR. Mode>6 ) THEN
         Info = -7
      ELSEIF ( (Mode/=-6 .AND. Mode/=0 .AND. Mode/=6) .AND. Cond<ONE )  &
     &         THEN
         Info = -8
      ELSEIF ( (Mode/=-6 .AND. Mode/=0 .AND. Mode/=6) .AND. irsign==-1 )&
     &         THEN
         Info = -10
      ELSEIF ( igrade==-1 .OR. (igrade==4 .AND. M/=N) .OR.              &
     &         ((igrade>=1 .AND. igrade<=4) .AND. isym==0) ) THEN
         Info = -11
      ELSEIF ( igrade==4 .AND. dzero ) THEN
         Info = -12
      ELSEIF ( (igrade==1 .OR. igrade==3 .OR. igrade==4 .OR. igrade==5) &
     &         .AND. (Model<-6 .OR. Model>6) ) THEN
         Info = -13
      ELSEIF ( (igrade==1 .OR. igrade==3 .OR. igrade==4 .OR. igrade==5) &
     &         .AND. (Model/=-6 .AND. Model/=0 .AND. Model/=6) .AND.    &
     &         Condl<ONE ) THEN
         Info = -14
      ELSEIF ( (igrade==2 .OR. igrade==3) .AND. (Moder<-6 .OR. Moder>6) &
     &         ) THEN
         Info = -16
      ELSEIF ( (igrade==2 .OR. igrade==3) .AND.                         &
     &         (Moder/=-6 .AND. Moder/=0 .AND. Moder/=6) .AND.          &
     &         Condr<ONE ) THEN
         Info = -17
      ELSEIF ( ipvtng==-1 .OR. (ipvtng==3 .AND. M/=N) .OR.              &
     &         ((ipvtng==1 .OR. ipvtng==2) .AND. isym==0) ) THEN
         Info = -18
      ELSEIF ( ipvtng/=0 .AND. badpvt ) THEN
         Info = -19
      ELSEIF ( Kl<0 ) THEN
         Info = -20
      ELSEIF ( Ku<0 .OR. (isym==0 .AND. Kl/=Ku) ) THEN
         Info = -21
      ELSEIF ( Sparse<ZERO .OR. Sparse>ONE ) THEN
         Info = -22
      ELSEIF ( ipack==-1 .OR.                                           &
     &         ((ipack==1 .OR. ipack==2 .OR. ipack==5 .OR. ipack==6)    &
     &         .AND. isym==1) .OR.                                      &
     &         (ipack==3 .AND. isym==1 .AND. (Kl/=0 .OR. M/=N)) .OR.    &
     &         (ipack==4 .AND. isym==1 .AND. (Ku/=0 .OR. M/=N)) ) THEN
         Info = -24
      ELSEIF ( ((ipack==0 .OR. ipack==1 .OR. ipack==2) .AND.            &
     &         Lda<MAX(1,M)) .OR. ((ipack==3 .OR. ipack==4) .AND. Lda<1)&
     &         .OR. ((ipack==5 .OR. ipack==6) .AND. Lda<kuu+1) .OR.     &
     &         (ipack==7 .AND. Lda<kll+kuu+1) ) THEN
         Info = -26
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLATMR',-Info)
         RETURN
      ENDIF
!
!     Decide if we can pivot consistently
!
      fulbnd = .FALSE.
      IF ( kuu==N-1 .AND. kll==M-1 ) fulbnd = .TRUE.
!
!     Initialize random number generator
!
      DO i = 1 , 4
         Iseed(i) = MOD(ABS(Iseed(i)),4096)
      ENDDO
!
      Iseed(4) = 2*(Iseed(4)/2) + 1
!
!     2)      Set up D, DL, and DR, if indicated.
!
!             Compute D according to COND and MODE
!
      CALL SLATM1(Mode,Cond,irsign,idist,Iseed,D,mnmin,Info)
      IF ( Info/=0 ) THEN
         Info = 1
         RETURN
      ENDIF
      IF ( Mode/=0 .AND. Mode/=-6 .AND. Mode/=6 ) THEN
!
!        Scale by DMAX
!
         temp = ABS(D(1))
         DO i = 2 , mnmin
            temp = MAX(temp,ABS(D(i)))
         ENDDO
         IF ( temp==ZERO .AND. Dmax/=ZERO ) THEN
            Info = 2
            RETURN
         ENDIF
         IF ( temp/=ZERO ) THEN
            alpha = Dmax/temp
         ELSE
            alpha = ONE
         ENDIF
         DO i = 1 , mnmin
            D(i) = alpha*D(i)
         ENDDO
!
      ENDIF
!
!     Compute DL if grading set
!
      IF ( igrade==1 .OR. igrade==3 .OR. igrade==4 .OR. igrade==5 ) THEN
         CALL SLATM1(Model,Condl,0,idist,Iseed,Dl,M,Info)
         IF ( Info/=0 ) THEN
            Info = 3
            RETURN
         ENDIF
      ENDIF
!
!     Compute DR if grading set
!
      IF ( igrade==2 .OR. igrade==3 ) THEN
         CALL SLATM1(Moder,Condr,0,idist,Iseed,Dr,N,Info)
         IF ( Info/=0 ) THEN
            Info = 4
            RETURN
         ENDIF
      ENDIF
!
!     3)     Generate IWORK if pivoting
!
      IF ( ipvtng>0 ) THEN
         DO i = 1 , npvts
            Iwork(i) = i
         ENDDO
         IF ( fulbnd ) THEN
            DO i = 1 , npvts
               k = Ipivot(i)
               j = Iwork(i)
               Iwork(i) = Iwork(k)
               Iwork(k) = j
            ENDDO
         ELSE
            DO i = npvts , 1 , -1
               k = Ipivot(i)
               j = Iwork(i)
               Iwork(i) = Iwork(k)
               Iwork(k) = j
            ENDDO
         ENDIF
      ENDIF
!
!     4)      Generate matrices for each kind of PACKing
!             Always sweep matrix columnwise (if symmetric, upper
!             half only) so that matrix generated does not depend
!             on PACK
!
      IF ( fulbnd ) THEN
!
!        Use SLATM3 so matrices generated with differing PIVOTing only
!        differ only in the order of their rows and/or columns.
!
         IF ( ipack==0 ) THEN
            IF ( isym==0 ) THEN
               DO j = 1 , N
                  DO i = 1 , j
                     temp = SLATM3(M,N,i,j,isub,jsub,Kl,Ku,idist,Iseed, &
     &                      D,igrade,Dl,Dr,ipvtng,Iwork,Sparse)
                     A(isub,jsub) = temp
                     A(jsub,isub) = temp
                  ENDDO
               ENDDO
            ELSEIF ( isym==1 ) THEN
               DO j = 1 , N
                  DO i = 1 , M
                     temp = SLATM3(M,N,i,j,isub,jsub,Kl,Ku,idist,Iseed, &
     &                      D,igrade,Dl,Dr,ipvtng,Iwork,Sparse)
                     A(isub,jsub) = temp
                  ENDDO
               ENDDO
            ENDIF
!
         ELSEIF ( ipack==1 ) THEN
!
            DO j = 1 , N
               DO i = 1 , j
                  temp = SLATM3(M,N,i,j,isub,jsub,Kl,Ku,idist,Iseed,D,  &
     &                   igrade,Dl,Dr,ipvtng,Iwork,Sparse)
                  mnsub = MIN(isub,jsub)
                  mxsub = MAX(isub,jsub)
                  A(mnsub,mxsub) = temp
                  IF ( mnsub/=mxsub ) A(mxsub,mnsub) = ZERO
               ENDDO
            ENDDO
!
         ELSEIF ( ipack==2 ) THEN
!
            DO j = 1 , N
               DO i = 1 , j
                  temp = SLATM3(M,N,i,j,isub,jsub,Kl,Ku,idist,Iseed,D,  &
     &                   igrade,Dl,Dr,ipvtng,Iwork,Sparse)
                  mnsub = MIN(isub,jsub)
                  mxsub = MAX(isub,jsub)
                  A(mxsub,mnsub) = temp
                  IF ( mnsub/=mxsub ) A(mnsub,mxsub) = ZERO
               ENDDO
            ENDDO
!
         ELSEIF ( ipack==3 ) THEN
!
            DO j = 1 , N
               DO i = 1 , j
                  temp = SLATM3(M,N,i,j,isub,jsub,Kl,Ku,idist,Iseed,D,  &
     &                   igrade,Dl,Dr,ipvtng,Iwork,Sparse)
!
!                 Compute K = location of (ISUB,JSUB) entry in packed
!                 array
!
                  mnsub = MIN(isub,jsub)
                  mxsub = MAX(isub,jsub)
                  k = mxsub*(mxsub-1)/2 + mnsub
!
!                 Convert K to (IISUB,JJSUB) location
!
                  jjsub = (k-1)/Lda + 1
                  iisub = k - Lda*(jjsub-1)
!
                  A(iisub,jjsub) = temp
               ENDDO
            ENDDO
!
         ELSEIF ( ipack==4 ) THEN
!
            DO j = 1 , N
               DO i = 1 , j
                  temp = SLATM3(M,N,i,j,isub,jsub,Kl,Ku,idist,Iseed,D,  &
     &                   igrade,Dl,Dr,ipvtng,Iwork,Sparse)
!
!                 Compute K = location of (I,J) entry in packed array
!
                  mnsub = MIN(isub,jsub)
                  mxsub = MAX(isub,jsub)
                  IF ( mnsub==1 ) THEN
                     k = mxsub
                  ELSE
                     k = N*(N+1)/2 - (N-mnsub+1)*(N-mnsub+2)/2 + mxsub -&
     &                   mnsub + 1
                  ENDIF
!
!                 Convert K to (IISUB,JJSUB) location
!
                  jjsub = (k-1)/Lda + 1
                  iisub = k - Lda*(jjsub-1)
!
                  A(iisub,jjsub) = temp
               ENDDO
            ENDDO
!
         ELSEIF ( ipack==5 ) THEN
!
            DO j = 1 , N
               DO i = j - kuu , j
                  IF ( i<1 ) THEN
                     A(j-i+1,i+N) = ZERO
                  ELSE
                     temp = SLATM3(M,N,i,j,isub,jsub,Kl,Ku,idist,Iseed, &
     &                      D,igrade,Dl,Dr,ipvtng,Iwork,Sparse)
                     mnsub = MIN(isub,jsub)
                     mxsub = MAX(isub,jsub)
                     A(mxsub-mnsub+1,mnsub) = temp
                  ENDIF
               ENDDO
            ENDDO
!
         ELSEIF ( ipack==6 ) THEN
!
            DO j = 1 , N
               DO i = j - kuu , j
                  temp = SLATM3(M,N,i,j,isub,jsub,Kl,Ku,idist,Iseed,D,  &
     &                   igrade,Dl,Dr,ipvtng,Iwork,Sparse)
                  mnsub = MIN(isub,jsub)
                  mxsub = MAX(isub,jsub)
                  A(mnsub-mxsub+kuu+1,mxsub) = temp
               ENDDO
            ENDDO
!
         ELSEIF ( ipack==7 ) THEN
!
            IF ( isym==0 ) THEN
               DO j = 1 , N
                  DO i = j - kuu , j
                     temp = SLATM3(M,N,i,j,isub,jsub,Kl,Ku,idist,Iseed, &
     &                      D,igrade,Dl,Dr,ipvtng,Iwork,Sparse)
                     mnsub = MIN(isub,jsub)
                     mxsub = MAX(isub,jsub)
                     A(mnsub-mxsub+kuu+1,mxsub) = temp
                     IF ( i<1 ) A(j-i+1+kuu,i+N) = ZERO
                     IF ( i>=1 .AND. mnsub/=mxsub )                     &
     &                    A(mxsub-mnsub+1+kuu,mnsub) = temp
                  ENDDO
               ENDDO
            ELSEIF ( isym==1 ) THEN
               DO j = 1 , N
                  DO i = j - kuu , j + kll
                     temp = SLATM3(M,N,i,j,isub,jsub,Kl,Ku,idist,Iseed, &
     &                      D,igrade,Dl,Dr,ipvtng,Iwork,Sparse)
                     A(isub-jsub+kuu+1,jsub) = temp
                  ENDDO
               ENDDO
            ENDIF
!
         ENDIF
!
!
!        Use SLATM2
!
      ELSEIF ( ipack==0 ) THEN
         IF ( isym==0 ) THEN
            DO j = 1 , N
               DO i = 1 , j
                  A(i,j) = SLATM2(M,N,i,j,Kl,Ku,idist,Iseed,D,igrade,Dl,&
     &                     Dr,ipvtng,Iwork,Sparse)
                  A(j,i) = A(i,j)
               ENDDO
            ENDDO
         ELSEIF ( isym==1 ) THEN
            DO j = 1 , N
               DO i = 1 , M
                  A(i,j) = SLATM2(M,N,i,j,Kl,Ku,idist,Iseed,D,igrade,Dl,&
     &                     Dr,ipvtng,Iwork,Sparse)
               ENDDO
            ENDDO
         ENDIF
!
      ELSEIF ( ipack==1 ) THEN
!
         DO j = 1 , N
            DO i = 1 , j
               A(i,j) = SLATM2(M,N,i,j,Kl,Ku,idist,Iseed,D,igrade,Dl,Dr,&
     &                  ipvtng,Iwork,Sparse)
               IF ( i/=j ) A(j,i) = ZERO
            ENDDO
         ENDDO
!
      ELSEIF ( ipack==2 ) THEN
!
         DO j = 1 , N
            DO i = 1 , j
               A(j,i) = SLATM2(M,N,i,j,Kl,Ku,idist,Iseed,D,igrade,Dl,Dr,&
     &                  ipvtng,Iwork,Sparse)
               IF ( i/=j ) A(i,j) = ZERO
            ENDDO
         ENDDO
!
      ELSEIF ( ipack==3 ) THEN
!
         isub = 0
         jsub = 1
         DO j = 1 , N
            DO i = 1 , j
               isub = isub + 1
               IF ( isub>Lda ) THEN
                  isub = 1
                  jsub = jsub + 1
               ENDIF
               A(isub,jsub) = SLATM2(M,N,i,j,Kl,Ku,idist,Iseed,D,igrade,&
     &                        Dl,Dr,ipvtng,Iwork,Sparse)
            ENDDO
         ENDDO
!
      ELSEIF ( ipack==4 ) THEN
!
         IF ( isym==0 ) THEN
            DO j = 1 , N
               DO i = 1 , j
!
!                    Compute K = location of (I,J) entry in packed array
!
                  IF ( i==1 ) THEN
                     k = j
                  ELSE
                     k = N*(N+1)/2 - (N-i+1)*(N-i+2)/2 + j - i + 1
                  ENDIF
!
!                    Convert K to (ISUB,JSUB) location
!
                  jsub = (k-1)/Lda + 1
                  isub = k - Lda*(jsub-1)
!
                  A(isub,jsub) = SLATM2(M,N,i,j,Kl,Ku,idist,Iseed,D,    &
     &                           igrade,Dl,Dr,ipvtng,Iwork,Sparse)
               ENDDO
            ENDDO
         ELSE
            isub = 0
            jsub = 1
            DO j = 1 , N
               DO i = j , M
                  isub = isub + 1
                  IF ( isub>Lda ) THEN
                     isub = 1
                     jsub = jsub + 1
                  ENDIF
                  A(isub,jsub) = SLATM2(M,N,i,j,Kl,Ku,idist,Iseed,D,    &
     &                           igrade,Dl,Dr,ipvtng,Iwork,Sparse)
               ENDDO
            ENDDO
         ENDIF
!
      ELSEIF ( ipack==5 ) THEN
!
         DO j = 1 , N
            DO i = j - kuu , j
               IF ( i<1 ) THEN
                  A(j-i+1,i+N) = ZERO
               ELSE
                  A(j-i+1,i) = SLATM2(M,N,i,j,Kl,Ku,idist,Iseed,D,      &
     &                         igrade,Dl,Dr,ipvtng,Iwork,Sparse)
               ENDIF
            ENDDO
         ENDDO
!
      ELSEIF ( ipack==6 ) THEN
!
         DO j = 1 , N
            DO i = j - kuu , j
               A(i-j+kuu+1,j) = SLATM2(M,N,i,j,Kl,Ku,idist,Iseed,D,     &
     &                          igrade,Dl,Dr,ipvtng,Iwork,Sparse)
            ENDDO
         ENDDO
!
      ELSEIF ( ipack==7 ) THEN
!
         IF ( isym==0 ) THEN
            DO j = 1 , N
               DO i = j - kuu , j
                  A(i-j+kuu+1,j) = SLATM2(M,N,i,j,Kl,Ku,idist,Iseed,D,  &
     &                             igrade,Dl,Dr,ipvtng,Iwork,Sparse)
                  IF ( i<1 ) A(j-i+1+kuu,i+N) = ZERO
                  IF ( i>=1 .AND. i/=j ) A(j-i+1+kuu,i) = A(i-j+kuu+1,j)
               ENDDO
            ENDDO
         ELSEIF ( isym==1 ) THEN
            DO j = 1 , N
               DO i = j - kuu , j + kll
                  A(i-j+kuu+1,j) = SLATM2(M,N,i,j,Kl,Ku,idist,Iseed,D,  &
     &                             igrade,Dl,Dr,ipvtng,Iwork,Sparse)
               ENDDO
            ENDDO
         ENDIF
!
!
      ENDIF
!
!     5)      Scaling the norm
!
      IF ( ipack==0 ) THEN
         onorm = SLANGE('M',M,N,A,Lda,tempa)
      ELSEIF ( ipack==1 ) THEN
         onorm = SLANSY('M','U',N,A,Lda,tempa)
      ELSEIF ( ipack==2 ) THEN
         onorm = SLANSY('M','L',N,A,Lda,tempa)
      ELSEIF ( ipack==3 ) THEN
         onorm = SLANSP('M','U',N,A,tempa)
      ELSEIF ( ipack==4 ) THEN
         onorm = SLANSP('M','L',N,A,tempa)
      ELSEIF ( ipack==5 ) THEN
         onorm = SLANSB('M','L',N,kll,A,Lda,tempa)
      ELSEIF ( ipack==6 ) THEN
         onorm = SLANSB('M','U',N,kuu,A,Lda,tempa)
      ELSEIF ( ipack==7 ) THEN
         onorm = SLANGB('M',N,kll,kuu,A,Lda,tempa)
      ENDIF
!
      IF ( Anorm>=ZERO ) THEN
!
         IF ( Anorm>ZERO .AND. onorm==ZERO ) THEN
!
!           Desired scaling impossible
!
            Info = 5
            RETURN
!
         ELSEIF ( (Anorm>ONE .AND. onorm<ONE) .OR.                      &
     &            (Anorm<ONE .AND. onorm>ONE) ) THEN
!
!           Scale carefully to avoid over / underflow
!
            IF ( ipack<=2 ) THEN
               DO j = 1 , N
                  CALL SSCAL(M,ONE/onorm,A(1,j),1)
                  CALL SSCAL(M,Anorm,A(1,j),1)
               ENDDO
!
            ELSEIF ( ipack==3 .OR. ipack==4 ) THEN
!
               CALL SSCAL(N*(N+1)/2,ONE/onorm,A,1)
               CALL SSCAL(N*(N+1)/2,Anorm,A,1)
!
            ELSEIF ( ipack>=5 ) THEN
!
               DO j = 1 , N
                  CALL SSCAL(kll+kuu+1,ONE/onorm,A(1,j),1)
                  CALL SSCAL(kll+kuu+1,Anorm,A(1,j),1)
               ENDDO
!
            ENDIF
!
!
!           Scale straightforwardly
!
         ELSEIF ( ipack<=2 ) THEN
            DO j = 1 , N
               CALL SSCAL(M,Anorm/onorm,A(1,j),1)
            ENDDO
!
         ELSEIF ( ipack==3 .OR. ipack==4 ) THEN
!
            CALL SSCAL(N*(N+1)/2,Anorm/onorm,A,1)
!
         ELSEIF ( ipack>=5 ) THEN
!
            DO j = 1 , N
               CALL SSCAL(kll+kuu+1,Anorm/onorm,A(1,j),1)
            ENDDO
!
         ENDIF
!
      ENDIF
!
!     End of SLATMR
!
      END SUBROUTINE SLATMR
