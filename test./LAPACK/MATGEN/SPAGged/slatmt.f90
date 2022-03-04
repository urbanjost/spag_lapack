!*==slatmt.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SLATMT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLATMT( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX,
!                          RANK, KL, KU, PACK, A, LDA, WORK, INFO )
!
!       .. Scalar Arguments ..
!       REAL               COND, DMAX
!       INTEGER            INFO, KL, KU, LDA, M, MODE, N, RANK
!       CHARACTER          DIST, PACK, SYM
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), D( * ), WORK( * )
!       INTEGER            ISEED( 4 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SLATMT generates random matrices with specified singular values
!>    (or symmetric/hermitian with specified eigenvalues)
!>    for testing LAPACK programs.
!>
!>    SLATMT operates by applying the following sequence of
!>    operations:
!>
!>      Set the diagonal to D, where D may be input or
!>         computed according to MODE, COND, DMAX, and SYM
!>         as described below.
!>
!>      Generate a matrix with the appropriate band structure, by one
!>         of two methods:
!>
!>      Method A:
!>          Generate a dense M x N matrix by multiplying D on the left
!>              and the right by random unitary matrices, then:
!>
!>          Reduce the bandwidth according to KL and KU, using
!>          Householder transformations.
!>
!>      Method B:
!>          Convert the bandwidth-0 (i.e., diagonal) matrix to a
!>              bandwidth-1 matrix using Givens rotations, "chasing"
!>              out-of-band elements back, much as in QR; then
!>              convert the bandwidth-1 to a bandwidth-2 matrix, etc.
!>              Note that for reasonably small bandwidths (relative to
!>              M and N) this requires less storage, as a dense matrix
!>              is not generated.  Also, for symmetric matrices, only
!>              one triangle is generated.
!>
!>      Method A is chosen if the bandwidth is a large fraction of the
!>          order of the matrix, and LDA is at least M (so a dense
!>          matrix can be stored.)  Method B is chosen if the bandwidth
!>          is small (< 1/2 N for symmetric, < .3 N+M for
!>          non-symmetric), or LDA is less than M and not less than the
!>          bandwidth.
!>
!>      Pack the matrix if desired. Options specified by PACK are:
!>         no packing
!>         zero out upper half (if symmetric)
!>         zero out lower half (if symmetric)
!>         store the upper half columnwise (if symmetric or upper
!>               triangular)
!>         store the lower half columnwise (if symmetric or lower
!>               triangular)
!>         store the lower triangle in banded format (if symmetric
!>               or lower triangular)
!>         store the upper triangle in banded format (if symmetric
!>               or upper triangular)
!>         store the entire matrix in banded format
!>      If Method B is chosen, and band format is specified, then the
!>         matrix will be generated in the band format, so no repacking
!>         will be necessary.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           The number of rows of A. Not modified.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The number of columns of A. Not modified.
!> \endverbatim
!>
!> \param[in] DIST
!> \verbatim
!>          DIST is CHARACTER*1
!>           On entry, DIST specifies the type of distribution to be used
!>           to generate the random eigen-/singular values.
!>           'U' => UNIFORM( 0, 1 )  ( 'U' for uniform )
!>           'S' => UNIFORM( -1, 1 ) ( 'S' for symmetric )
!>           'N' => NORMAL( 0, 1 )   ( 'N' for normal )
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension ( 4 )
!>           On entry ISEED specifies the seed of the random number
!>           generator. They should lie between 0 and 4095 inclusive,
!>           and ISEED(4) should be odd. The random number generator
!>           uses a linear congruential sequence limited to small
!>           integers, and so should produce machine independent
!>           random numbers. The values of ISEED are changed on
!>           exit, and can be used in the next call to SLATMT
!>           to continue the same random number sequence.
!>           Changed on exit.
!> \endverbatim
!>
!> \param[in] SYM
!> \verbatim
!>          SYM is CHARACTER*1
!>           If SYM='S' or 'H', the generated matrix is symmetric, with
!>             eigenvalues specified by D, COND, MODE, and DMAX; they
!>             may be positive, negative, or zero.
!>           If SYM='P', the generated matrix is symmetric, with
!>             eigenvalues (= singular values) specified by D, COND,
!>             MODE, and DMAX; they will not be negative.
!>           If SYM='N', the generated matrix is nonsymmetric, with
!>             singular values specified by D, COND, MODE, and DMAX;
!>             they will not be negative.
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension ( MIN( M , N ) )
!>           This array is used to specify the singular values or
!>           eigenvalues of A (see SYM, above.)  If MODE=0, then D is
!>           assumed to contain the singular/eigenvalues, otherwise
!>           they will be computed according to MODE, COND, and DMAX,
!>           and placed in D.
!>           Modified if MODE is nonzero.
!> \endverbatim
!>
!> \param[in] MODE
!> \verbatim
!>          MODE is INTEGER
!>           On entry this describes how the singular/eigenvalues are to
!>           be specified:
!>           MODE = 0 means use D as input
!>
!>           MODE = 1 sets D(1)=1 and D(2:RANK)=1.0/COND
!>           MODE = 2 sets D(1:RANK-1)=1 and D(RANK)=1.0/COND
!>           MODE = 3 sets D(I)=COND**(-(I-1)/(RANK-1))
!>
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
!>           If SYM='S' or 'H', and MODE is neither 0, 6, nor -6, then
!>              the elements of D will also be multiplied by a random
!>              sign (i.e., +1 or -1.)
!>           Not modified.
!> \endverbatim
!>
!> \param[in] COND
!> \verbatim
!>          COND is REAL
!>           On entry, this is used as described under MODE above.
!>           If used, it must be >= 1. Not modified.
!> \endverbatim
!>
!> \param[in] DMAX
!> \verbatim
!>          DMAX is REAL
!>           If MODE is neither -6, 0 nor 6, the contents of D, as
!>           computed according to MODE and COND, will be scaled by
!>           DMAX / max(abs(D(i))); thus, the maximum absolute eigen- or
!>           singular value (which is to say the norm) will be abs(DMAX).
!>           Note that DMAX need not be positive: if DMAX is negative
!>           (or zero), D will be scaled by a negative number (or zero).
!>           Not modified.
!> \endverbatim
!>
!> \param[in] RANK
!> \verbatim
!>          RANK is INTEGER
!>           The rank of matrix to be generated for modes 1,2,3 only.
!>           D( RANK+1:N ) = 0.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>           This specifies the lower bandwidth of the  matrix. For
!>           example, KL=0 implies upper triangular, KL=1 implies upper
!>           Hessenberg, and KL being at least M-1 means that the matrix
!>           has full lower bandwidth.  KL must equal KU if the matrix
!>           is symmetric.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>           This specifies the upper bandwidth of the  matrix. For
!>           example, KU=0 implies lower triangular, KU=1 implies lower
!>           Hessenberg, and KU being at least N-1 means that the matrix
!>           has full upper bandwidth.  KL must equal KU if the matrix
!>           is symmetric.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] PACK
!> \verbatim
!>          PACK is CHARACTER*1
!>           This specifies packing of matrix as follows:
!>           'N' => no packing
!>           'U' => zero out all subdiagonal entries (if symmetric)
!>           'L' => zero out all superdiagonal entries (if symmetric)
!>           'C' => store the upper triangle columnwise
!>                  (only if the matrix is symmetric or upper triangular)
!>           'R' => store the lower triangle columnwise
!>                  (only if the matrix is symmetric or lower triangular)
!>           'B' => store the lower triangle in band storage scheme
!>                  (only if matrix symmetric or lower triangular)
!>           'Q' => store the upper triangle in band storage scheme
!>                  (only if matrix symmetric or upper triangular)
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
!>           If two calls to SLATMT differ only in the PACK parameter,
!>           they will generate mathematically equivalent matrices.
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension ( LDA, N )
!>           On exit A is the desired test matrix.  A is first generated
!>           in full (unpacked) form, and then packed, if so specified
!>           by PACK.  Thus, the first M elements of the first N
!>           columns will always be modified.  If PACK specifies a
!>           packed or banded storage scheme, all LDA elements of the
!>           first N columns will be modified; the elements of the
!>           array which do not correspond to elements of the generated
!>           matrix are set to zero.
!>           Modified.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           LDA specifies the first dimension of A as declared in the
!>           calling program.  If PACK='N', 'U', 'L', 'C', or 'R', then
!>           LDA must be at least M.  If PACK='B' or 'Q', then LDA must
!>           be at least MIN( KL, M-1) (which is equal to MIN(KU,N-1)).
!>           If PACK='Z', LDA must be large enough to hold the packed
!>           array: MIN( KU, N-1) + MIN( KL, M-1) + 1.
!>           Not modified.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension ( 3*MAX( N , M ) )
!>           Workspace.
!>           Modified.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           Error code.  On exit, INFO will be set to one of the
!>           following values:
!>             0 => normal return
!>            -1 => M negative or unequal to N and SYM='S', 'H', or 'P'
!>            -2 => N negative
!>            -3 => DIST illegal string
!>            -5 => SYM illegal string
!>            -7 => MODE not in range -6 to 6
!>            -8 => COND less than 1.0, and MODE neither -6, 0 nor 6
!>           -10 => KL negative
!>           -11 => KU negative, or SYM='S' or 'H' and KU not equal to KL
!>           -12 => PACK illegal string, or PACK='U' or 'L', and SYM='N';
!>                  or PACK='C' or 'Q' and SYM='N' and KL is not zero;
!>                  or PACK='R' or 'B' and SYM='N' and KU is not zero;
!>                  or PACK='U', 'L', 'C', 'R', 'B', or 'Q', and M is not
!>                  N.
!>           -14 => LDA is less than M, or PACK='Z' and LDA is less than
!>                  MIN(KU,N-1) + MIN(KL,M-1) + 1.
!>            1  => Error return from SLATM7
!>            2  => Cannot scale to DMAX (max. sing. value is 0)
!>            3  => Error return from SLAGGE or SLAGSY
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
      SUBROUTINE SLATMT(M,N,Dist,Iseed,Sym,D,Mode,Cond,Dmax,Rank,Kl,Ku, &
     &                  Pack,A,Lda,Work,Info)
      IMPLICIT NONE
!*--SLATMT335
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      REAL Cond , Dmax
      INTEGER Info , Kl , Ku , Lda , M , Mode , N , Rank
      CHARACTER Dist , Pack , Sym
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , D(*) , Work(*)
      INTEGER Iseed(4)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
      REAL ONE
      PARAMETER (ONE=1.0E0)
      REAL TWOPI
      PARAMETER (TWOPI=6.28318530717958647692528676655900576839E+0)
!     ..
!     .. Local Scalars ..
      REAL alpha , angle , c , dummy , extra , s , temp
      INTEGER i , ic , icol , idist , iendch , iinfo , il , ilda ,      &
     &        ioffg , ioffst , ipack , ipackg , ir , ir1 , ir2 , irow , &
     &        irsign , iskew , isym , isympk , j , jc , jch , jkl ,     &
     &        jku , jr , k , llb , minlda , mnmin , mr , nc , uub
      LOGICAL givens , ilextr , iltemp , topdwn
!     ..
!     .. External Functions ..
      REAL SLARND
      LOGICAL LSAME
      EXTERNAL SLARND , LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL SLATM7 , SCOPY , SLAGGE , SLAGSY , SLAROT , SLARTG ,     &
     &         SLASET , SSCAL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , COS , MAX , MIN , MOD , REAL , SIN
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
      IF ( LSAME(Sym,'N') ) THEN
         isym = 1
         irsign = 0
      ELSEIF ( LSAME(Sym,'P') ) THEN
         isym = 2
         irsign = 0
      ELSEIF ( LSAME(Sym,'S') ) THEN
         isym = 2
         irsign = 1
      ELSEIF ( LSAME(Sym,'H') ) THEN
         isym = 2
         irsign = 1
      ELSE
         isym = -1
      ENDIF
!
!     Decode PACK
!
      isympk = 0
      IF ( LSAME(Pack,'N') ) THEN
         ipack = 0
      ELSEIF ( LSAME(Pack,'U') ) THEN
         ipack = 1
         isympk = 1
      ELSEIF ( LSAME(Pack,'L') ) THEN
         ipack = 2
         isympk = 1
      ELSEIF ( LSAME(Pack,'C') ) THEN
         ipack = 3
         isympk = 2
      ELSEIF ( LSAME(Pack,'R') ) THEN
         ipack = 4
         isympk = 3
      ELSEIF ( LSAME(Pack,'B') ) THEN
         ipack = 5
         isympk = 3
      ELSEIF ( LSAME(Pack,'Q') ) THEN
         ipack = 6
         isympk = 2
      ELSEIF ( LSAME(Pack,'Z') ) THEN
         ipack = 7
      ELSE
         ipack = -1
      ENDIF
!
!     Set certain internal parameters
!
      mnmin = MIN(M,N)
      llb = MIN(Kl,M-1)
      uub = MIN(Ku,N-1)
      mr = MIN(M,N+llb)
      nc = MIN(N,M+uub)
!
      IF ( ipack==5 .OR. ipack==6 ) THEN
         minlda = uub + 1
      ELSEIF ( ipack==7 ) THEN
         minlda = llb + uub + 1
      ELSE
         minlda = M
      ENDIF
!
!     Use Givens rotation method if bandwidth small enough,
!     or if LDA is too small to store the matrix unpacked.
!
      givens = .FALSE.
      IF ( isym==1 ) THEN
         IF ( REAL(llb+uub)<0.3*REAL(MAX(1,mr+nc)) ) givens = .TRUE.
      ELSE
         IF ( 2*llb<M ) givens = .TRUE.
      ENDIF
      IF ( Lda<M .AND. Lda>=minlda ) givens = .TRUE.
!
!     Set INFO if an error
!
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( M/=N .AND. isym/=1 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( idist==-1 ) THEN
         Info = -3
      ELSEIF ( isym==-1 ) THEN
         Info = -5
      ELSEIF ( ABS(Mode)>6 ) THEN
         Info = -7
      ELSEIF ( (Mode/=0 .AND. ABS(Mode)/=6) .AND. Cond<ONE ) THEN
         Info = -8
      ELSEIF ( Kl<0 ) THEN
         Info = -10
      ELSEIF ( Ku<0 .OR. (isym/=1 .AND. Kl/=Ku) ) THEN
         Info = -11
      ELSEIF ( ipack==-1 .OR. (isympk==1 .AND. isym==1) .OR.            &
     &         (isympk==2 .AND. isym==1 .AND. Kl>0) .OR.                &
     &         (isympk==3 .AND. isym==1 .AND. Ku>0) .OR.                &
     &         (isympk/=0 .AND. M/=N) ) THEN
         Info = -12
      ELSEIF ( Lda<MAX(1,minlda) ) THEN
         Info = -14
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLATMT',-Info)
         RETURN
      ENDIF
!
!     Initialize random number generator
!
      DO i = 1 , 4
         Iseed(i) = MOD(ABS(Iseed(i)),4096)
      ENDDO
!
      IF ( MOD(Iseed(4),2)/=1 ) Iseed(4) = Iseed(4) + 1
!
!     2)      Set up D  if indicated.
!
!             Compute D according to COND and MODE
!
      CALL SLATM7(Mode,Cond,irsign,idist,Iseed,D,mnmin,Rank,iinfo)
      IF ( iinfo/=0 ) THEN
         Info = 1
         RETURN
      ENDIF
!
!     Choose Top-Down if D is (apparently) increasing,
!     Bottom-Up if D is (apparently) decreasing.
!
      IF ( ABS(D(1))<=ABS(D(Rank)) ) THEN
         topdwn = .TRUE.
      ELSE
         topdwn = .FALSE.
      ENDIF
!
      IF ( Mode/=0 .AND. ABS(Mode)/=6 ) THEN
!
!        Scale by DMAX
!
         temp = ABS(D(1))
         DO i = 2 , Rank
            temp = MAX(temp,ABS(D(i)))
         ENDDO
!
         IF ( temp>ZERO ) THEN
            alpha = Dmax/temp
         ELSE
            Info = 2
            RETURN
         ENDIF
!
         CALL SSCAL(Rank,alpha,D,1)
!
      ENDIF
!
!     3)      Generate Banded Matrix using Givens rotations.
!             Also the special case of UUB=LLB=0
!
!               Compute Addressing constants to cover all
!               storage formats.  Whether GE, SY, GB, or SB,
!               upper or lower triangle or both,
!               the (i,j)-th element is in
!               A( i - ISKEW*j + IOFFST, j )
!
      IF ( ipack>4 ) THEN
         ilda = Lda - 1
         iskew = 1
         IF ( ipack>5 ) THEN
            ioffst = uub + 1
         ELSE
            ioffst = 1
         ENDIF
      ELSE
         ilda = Lda
         iskew = 0
         ioffst = 0
      ENDIF
!
!     IPACKG is the format that the matrix is generated in. If this is
!     different from IPACK, then the matrix must be repacked at the
!     end.  It also signals how to compute the norm, for scaling.
!
      ipackg = 0
      CALL SLASET('Full',Lda,N,ZERO,ZERO,A,Lda)
!
!     Diagonal Matrix -- We are done, unless it
!     is to be stored SP/PP/TP (PACK='R' or 'C')
!
      IF ( llb==0 .AND. uub==0 ) THEN
         CALL SCOPY(mnmin,D,1,A(1-iskew+ioffst,1),ilda+1)
         IF ( ipack<=2 .OR. ipack>=5 ) ipackg = ipack
!
      ELSEIF ( givens ) THEN
!
!        Check whether to use Givens rotations,
!        Householder transformations, or nothing.
!
         IF ( isym==1 ) THEN
!
!           Non-symmetric -- A = U D V
!
            IF ( ipack>4 ) THEN
               ipackg = ipack
            ELSE
               ipackg = 0
            ENDIF
!
            CALL SCOPY(mnmin,D,1,A(1-iskew+ioffst,1),ilda+1)
!
            IF ( topdwn ) THEN
               jkl = 0
               DO jku = 1 , uub
!
!                 Transform from bandwidth JKL, JKU-1 to JKL, JKU
!
!                 Last row actually rotated is M
!                 Last column actually rotated is MIN( M+JKU, N )
!
                  DO jr = 1 , MIN(M+jku,N) + jkl - 1
                     extra = ZERO
                     angle = TWOPI*SLARND(1,Iseed)
                     c = COS(angle)
                     s = SIN(angle)
                     icol = MAX(1,jr-jkl)
                     IF ( jr<M ) THEN
                        il = MIN(N,jr+jku) + 1 - icol
                        CALL SLAROT(.TRUE.,jr>jkl,.FALSE.,il,c,s,       &
     &                              A(jr-iskew*icol+ioffst,icol),ilda,  &
     &                              extra,dummy)
                     ENDIF
!
!                    Chase "EXTRA" back up
!
                     ir = jr
                     ic = icol
                     DO jch = jr - jkl , 1 , -jkl - jku
                        IF ( ir<M )                                     &
     &                       CALL SLARTG(A(ir+1-iskew*(ic+1)+ioffst,    &
     &                       ic+1),extra,c,s,dummy)
                        irow = MAX(1,jch-jku)
                        il = ir + 2 - irow
                        temp = ZERO
                        iltemp = jch>jku
                        CALL SLAROT(.FALSE.,iltemp,.TRUE.,il,c,-s,      &
     &                              A(irow-iskew*ic+ioffst,ic),ilda,    &
     &                              temp,extra)
                        IF ( iltemp ) THEN
                           CALL SLARTG(A(irow+1-iskew*(ic+1)+ioffst,ic+1&
     &                                 ),temp,c,s,dummy)
                           icol = MAX(1,jch-jku-jkl)
                           il = ic + 2 - icol
                           extra = ZERO
                           CALL SLAROT(.TRUE.,jch>jku+jkl,.TRUE.,il,c,  &
     &                                 -s,A(irow-iskew*icol+ioffst,icol)&
     &                                 ,ilda,extra,temp)
                           ic = icol
                           ir = irow
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
!
               jku = uub
               DO jkl = 1 , llb
!
!                 Transform from bandwidth JKL-1, JKU to JKL, JKU
!
                  DO jc = 1 , MIN(N+jkl,M) + jku - 1
                     extra = ZERO
                     angle = TWOPI*SLARND(1,Iseed)
                     c = COS(angle)
                     s = SIN(angle)
                     irow = MAX(1,jc-jku)
                     IF ( jc<N ) THEN
                        il = MIN(M,jc+jkl) + 1 - irow
                        CALL SLAROT(.FALSE.,jc>jku,.FALSE.,il,c,s,      &
     &                              A(irow-iskew*jc+ioffst,jc),ilda,    &
     &                              extra,dummy)
                     ENDIF
!
!                    Chase "EXTRA" back up
!
                     ic = jc
                     ir = irow
                     DO jch = jc - jku , 1 , -jkl - jku
                        IF ( ic<N )                                     &
     &                       CALL SLARTG(A(ir+1-iskew*(ic+1)+ioffst,    &
     &                       ic+1),extra,c,s,dummy)
                        icol = MAX(1,jch-jkl)
                        il = ic + 2 - icol
                        temp = ZERO
                        iltemp = jch>jkl
                        CALL SLAROT(.TRUE.,iltemp,.TRUE.,il,c,-s,       &
     &                              A(ir-iskew*icol+ioffst,icol),ilda,  &
     &                              temp,extra)
                        IF ( iltemp ) THEN
                           CALL SLARTG(A(ir+1-iskew*(icol+1)+ioffst,icol&
     &                                 +1),temp,c,s,dummy)
                           irow = MAX(1,jch-jkl-jku)
                           il = ir + 2 - irow
                           extra = ZERO
                           CALL SLAROT(.FALSE.,jch>jkl+jku,.TRUE.,il,c, &
     &                                 -s,A(irow-iskew*icol+ioffst,icol)&
     &                                 ,ilda,extra,temp)
                           ic = icol
                           ir = irow
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
!
            ELSE
!
!              Bottom-Up -- Start at the bottom right.
!
               jkl = 0
               DO jku = 1 , uub
!
!                 Transform from bandwidth JKL, JKU-1 to JKL, JKU
!
!                 First row actually rotated is M
!                 First column actually rotated is MIN( M+JKU, N )
!
                  iendch = MIN(M,N+jkl) - 1
                  DO jc = MIN(M+jku,N) - 1 , 1 - jkl , -1
                     extra = ZERO
                     angle = TWOPI*SLARND(1,Iseed)
                     c = COS(angle)
                     s = SIN(angle)
                     irow = MAX(1,jc-jku+1)
                     IF ( jc>0 ) THEN
                        il = MIN(M,jc+jkl+1) + 1 - irow
                        CALL SLAROT(.FALSE.,.FALSE.,jc+jkl<M,il,c,s,    &
     &                              A(irow-iskew*jc+ioffst,jc),ilda,    &
     &                              dummy,extra)
                     ENDIF
!
!                    Chase "EXTRA" back down
!
                     ic = jc
                     DO jch = jc + jkl , iendch , jkl + jku
                        ilextr = ic>0
                        IF ( ilextr )                                   &
     &                       CALL SLARTG(A(jch-iskew*ic+ioffst,ic),     &
     &                       extra,c,s,dummy)
                        ic = MAX(1,ic)
                        icol = MIN(N-1,jch+jku)
                        iltemp = jch + jku<N
                        temp = ZERO
                        CALL SLAROT(.TRUE.,ilextr,iltemp,icol+2-ic,c,s, &
     &                              A(jch-iskew*ic+ioffst,ic),ilda,     &
     &                              extra,temp)
                        IF ( iltemp ) THEN
                           CALL SLARTG(A(jch-iskew*icol+ioffst,icol),   &
     &                                 temp,c,s,dummy)
                           il = MIN(iendch,jch+jkl+jku) + 2 - jch
                           extra = ZERO
                           CALL SLAROT(.FALSE.,.TRUE.,                  &
     &                                 jch+jkl+jku<=iendch,il,c,s,      &
     &                                 A(jch-iskew*icol+ioffst,icol),   &
     &                                 ilda,temp,extra)
                           ic = icol
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
!
               jku = uub
               DO jkl = 1 , llb
!
!                 Transform from bandwidth JKL-1, JKU to JKL, JKU
!
!                 First row actually rotated is MIN( N+JKL, M )
!                 First column actually rotated is N
!
                  iendch = MIN(N,M+jku) - 1
                  DO jr = MIN(N+jkl,M) - 1 , 1 - jku , -1
                     extra = ZERO
                     angle = TWOPI*SLARND(1,Iseed)
                     c = COS(angle)
                     s = SIN(angle)
                     icol = MAX(1,jr-jkl+1)
                     IF ( jr>0 ) THEN
                        il = MIN(N,jr+jku+1) + 1 - icol
                        CALL SLAROT(.TRUE.,.FALSE.,jr+jku<N,il,c,s,     &
     &                              A(jr-iskew*icol+ioffst,icol),ilda,  &
     &                              dummy,extra)
                     ENDIF
!
!                    Chase "EXTRA" back down
!
                     ir = jr
                     DO jch = jr + jku , iendch , jkl + jku
                        ilextr = ir>0
                        IF ( ilextr )                                   &
     &                       CALL SLARTG(A(ir-iskew*jch+ioffst,jch),    &
     &                       extra,c,s,dummy)
                        ir = MAX(1,ir)
                        irow = MIN(M-1,jch+jkl)
                        iltemp = jch + jkl<M
                        temp = ZERO
                        CALL SLAROT(.FALSE.,ilextr,iltemp,irow+2-ir,c,s,&
     &                              A(ir-iskew*jch+ioffst,jch),ilda,    &
     &                              extra,temp)
                        IF ( iltemp ) THEN
                           CALL SLARTG(A(irow-iskew*jch+ioffst,jch),    &
     &                                 temp,c,s,dummy)
                           il = MIN(iendch,jch+jkl+jku) + 2 - jch
                           extra = ZERO
                           CALL SLAROT(.TRUE.,.TRUE.,                   &
     &                                 jch+jkl+jku<=iendch,il,c,s,      &
     &                                 A(irow-iskew*jch+ioffst,jch),    &
     &                                 ilda,temp,extra)
                           ir = irow
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
!
         ELSE
!
!           Symmetric -- A = U D U'
!
            ipackg = ipack
            ioffg = ioffst
!
            IF ( topdwn ) THEN
!
!              Top-Down -- Generate Upper triangle only
!
               IF ( ipack>=5 ) THEN
                  ipackg = 6
                  ioffg = uub + 1
               ELSE
                  ipackg = 1
               ENDIF
               CALL SCOPY(mnmin,D,1,A(1-iskew+ioffg,1),ilda+1)
!
               DO k = 1 , uub
                  DO jc = 1 , N - 1
                     irow = MAX(1,jc-k)
                     il = MIN(jc+1,k+2)
                     extra = ZERO
                     temp = A(jc-iskew*(jc+1)+ioffg,jc+1)
                     angle = TWOPI*SLARND(1,Iseed)
                     c = COS(angle)
                     s = SIN(angle)
                     CALL SLAROT(.FALSE.,jc>k,.TRUE.,il,c,s,            &
     &                           A(irow-iskew*jc+ioffg,jc),ilda,extra,  &
     &                           temp)
                     CALL SLAROT(.TRUE.,.TRUE.,.FALSE.,MIN(k,N-jc)+1,c, &
     &                           s,A((1-iskew)*jc+ioffg,jc),ilda,temp,  &
     &                           dummy)
!
!                    Chase EXTRA back up the matrix
!
                     icol = jc
                     DO jch = jc - k , 1 , -k
                        CALL SLARTG(A(jch+1-iskew*(icol+1)+ioffg,icol+1)&
     &                              ,extra,c,s,dummy)
                        temp = A(jch-iskew*(jch+1)+ioffg,jch+1)
                        CALL SLAROT(.TRUE.,.TRUE.,.TRUE.,k+2,c,-s,      &
     &                              A((1-iskew)*jch+ioffg,jch),ilda,    &
     &                              temp,extra)
                        irow = MAX(1,jch-k)
                        il = MIN(jch+1,k+2)
                        extra = ZERO
                        CALL SLAROT(.FALSE.,jch>k,.TRUE.,il,c,-s,       &
     &                              A(irow-iskew*jch+ioffg,jch),ilda,   &
     &                              extra,temp)
                        icol = jch
                     ENDDO
                  ENDDO
               ENDDO
!
!              If we need lower triangle, copy from upper. Note that
!              the order of copying is chosen to work for 'q' -> 'b'
!
               IF ( ipack/=ipackg .AND. ipack/=3 ) THEN
                  DO jc = 1 , N
                     irow = ioffst - iskew*jc
                     DO jr = jc , MIN(N,jc+uub)
                        A(jr+irow,jc) = A(jc-iskew*jr+ioffg,jr)
                     ENDDO
                  ENDDO
                  IF ( ipack==5 ) THEN
                     DO jc = N - uub + 1 , N
                        DO jr = N + 2 - jc , uub + 1
                           A(jr,jc) = ZERO
                        ENDDO
                     ENDDO
                  ENDIF
                  IF ( ipackg==6 ) THEN
                     ipackg = ipack
                  ELSE
                     ipackg = 0
                  ENDIF
               ENDIF
            ELSE
!
!              Bottom-Up -- Generate Lower triangle only
!
               IF ( ipack>=5 ) THEN
                  ipackg = 5
                  IF ( ipack==6 ) ioffg = 1
               ELSE
                  ipackg = 2
               ENDIF
               CALL SCOPY(mnmin,D,1,A(1-iskew+ioffg,1),ilda+1)
!
               DO k = 1 , uub
                  DO jc = N - 1 , 1 , -1
                     il = MIN(N+1-jc,k+2)
                     extra = ZERO
                     temp = A(1+(1-iskew)*jc+ioffg,jc)
                     angle = TWOPI*SLARND(1,Iseed)
                     c = COS(angle)
                     s = -SIN(angle)
                     CALL SLAROT(.FALSE.,.TRUE.,N-jc>k,il,c,s,          &
     &                           A((1-iskew)*jc+ioffg,jc),ilda,temp,    &
     &                           extra)
                     icol = MAX(1,jc-k+1)
                     CALL SLAROT(.TRUE.,.FALSE.,.TRUE.,jc+2-icol,c,s,   &
     &                           A(jc-iskew*icol+ioffg,icol),ilda,dummy,&
     &                           temp)
!
!                    Chase EXTRA back down the matrix
!
                     icol = jc
                     DO jch = jc + k , N - 1 , k
                        CALL SLARTG(A(jch-iskew*icol+ioffg,icol),extra, &
     &                              c,s,dummy)
                        temp = A(1+(1-iskew)*jch+ioffg,jch)
                        CALL SLAROT(.TRUE.,.TRUE.,.TRUE.,k+2,c,s,       &
     &                              A(jch-iskew*icol+ioffg,icol),ilda,  &
     &                              extra,temp)
                        il = MIN(N+1-jch,k+2)
                        extra = ZERO
                        CALL SLAROT(.FALSE.,.TRUE.,N-jch>k,il,c,s,      &
     &                              A((1-iskew)*jch+ioffg,jch),ilda,    &
     &                              temp,extra)
                        icol = jch
                     ENDDO
                  ENDDO
               ENDDO
!
!              If we need upper triangle, copy from lower. Note that
!              the order of copying is chosen to work for 'b' -> 'q'
!
               IF ( ipack/=ipackg .AND. ipack/=4 ) THEN
                  DO jc = N , 1 , -1
                     irow = ioffst - iskew*jc
                     DO jr = jc , MAX(1,jc-uub) , -1
                        A(jr+irow,jc) = A(jc-iskew*jr+ioffg,jr)
                     ENDDO
                  ENDDO
                  IF ( ipack==6 ) THEN
                     DO jc = 1 , uub
                        DO jr = 1 , uub + 1 - jc
                           A(jr,jc) = ZERO
                        ENDDO
                     ENDDO
                  ENDIF
                  IF ( ipackg==5 ) THEN
                     ipackg = ipack
                  ELSE
                     ipackg = 0
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
!
      ELSE
!
!        4)      Generate Banded Matrix by first
!                Rotating by random Unitary matrices,
!                then reducing the bandwidth using Householder
!                transformations.
!
!                Note: we should get here only if LDA .ge. N
!
         IF ( isym==1 ) THEN
!
!           Non-symmetric -- A = U D V
!
            CALL SLAGGE(mr,nc,llb,uub,D,A,Lda,Iseed,Work,iinfo)
         ELSE
!
!           Symmetric -- A = U D U'
!
            CALL SLAGSY(M,llb,D,A,Lda,Iseed,Work,iinfo)
!
         ENDIF
         IF ( iinfo/=0 ) THEN
            Info = 3
            RETURN
         ENDIF
      ENDIF
!
!     5)      Pack the matrix
!
      IF ( ipack/=ipackg ) THEN
         IF ( ipack==1 ) THEN
!
!           'U' -- Upper triangular, not packed
!
            DO j = 1 , M
               DO i = j + 1 , M
                  A(i,j) = ZERO
               ENDDO
            ENDDO
!
         ELSEIF ( ipack==2 ) THEN
!
!           'L' -- Lower triangular, not packed
!
            DO j = 2 , M
               DO i = 1 , j - 1
                  A(i,j) = ZERO
               ENDDO
            ENDDO
!
         ELSEIF ( ipack==3 ) THEN
!
!           'C' -- Upper triangle packed Columnwise.
!
            icol = 1
            irow = 0
            DO j = 1 , M
               DO i = 1 , j
                  irow = irow + 1
                  IF ( irow>Lda ) THEN
                     irow = 1
                     icol = icol + 1
                  ENDIF
                  A(irow,icol) = A(i,j)
               ENDDO
            ENDDO
!
         ELSEIF ( ipack==4 ) THEN
!
!           'R' -- Lower triangle packed Columnwise.
!
            icol = 1
            irow = 0
            DO j = 1 , M
               DO i = j , M
                  irow = irow + 1
                  IF ( irow>Lda ) THEN
                     irow = 1
                     icol = icol + 1
                  ENDIF
                  A(irow,icol) = A(i,j)
               ENDDO
            ENDDO
!
         ELSEIF ( ipack>=5 ) THEN
!
!           'B' -- The lower triangle is packed as a band matrix.
!           'Q' -- The upper triangle is packed as a band matrix.
!           'Z' -- The whole matrix is packed as a band matrix.
!
            IF ( ipack==5 ) uub = 0
            IF ( ipack==6 ) llb = 0
!
            DO j = 1 , uub
               DO i = MIN(j+llb,M) , 1 , -1
                  A(i-j+uub+1,j) = A(i,j)
               ENDDO
            ENDDO
!
            DO j = uub + 2 , N
               DO i = j - uub , MIN(j+llb,M)
                  A(i-j+uub+1,j) = A(i,j)
               ENDDO
            ENDDO
         ENDIF
!
!        If packed, zero out extraneous elements.
!
!        Symmetric/Triangular Packed --
!        zero out everything after A(IROW,ICOL)
!
         IF ( ipack==3 .OR. ipack==4 ) THEN
            DO jc = icol , M
               DO jr = irow + 1 , Lda
                  A(jr,jc) = ZERO
               ENDDO
               irow = 0
            ENDDO
!
         ELSEIF ( ipack>=5 ) THEN
!
!           Packed Band --
!              1st row is now in A( UUB+2-j, j), zero above it
!              m-th row is now in A( M+UUB-j,j), zero below it
!              last non-zero diagonal is now in A( UUB+LLB+1,j ),
!                 zero below it, too.
!
            ir1 = uub + llb + 2
            ir2 = uub + M + 2
            DO jc = 1 , N
               DO jr = 1 , uub + 1 - jc
                  A(jr,jc) = ZERO
               ENDDO
               DO jr = MAX(1,MIN(ir1,ir2-jc)) , Lda
                  A(jr,jc) = ZERO
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!
!
!     End of SLATMT
!
      END SUBROUTINE SLATMT
