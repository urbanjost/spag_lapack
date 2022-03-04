!*==dlatme.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DLATME
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLATME( N, DIST, ISEED, D, MODE, COND, DMAX, EI,
!         RSIGN,
!                          UPPER, SIM, DS, MODES, CONDS, KL, KU, ANORM,
!         A,
!                          LDA, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIST, RSIGN, SIM, UPPER
!       INTEGER            INFO, KL, KU, LDA, MODE, MODES, N
!       DOUBLE PRECISION   ANORM, COND, CONDS, DMAX
!       ..
!       .. Array Arguments ..
!       CHARACTER          EI( * )
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   A( LDA, * ), D( * ), DS( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLATME generates random non-symmetric square matrices with
!>    specified eigenvalues for testing LAPACK programs.
!>
!>    DLATME operates by applying the following sequence of
!>    operations:
!>
!>    1. Set the diagonal to D, where D may be input or
!>         computed according to MODE, COND, DMAX, and RSIGN
!>         as described below.
!>
!>    2. If complex conjugate pairs are desired (MODE=0 and EI(1)='R',
!>         or MODE=5), certain pairs of adjacent elements of D are
!>         interpreted as the real and complex parts of a complex
!>         conjugate pair; A thus becomes block diagonal, with 1x1
!>         and 2x2 blocks.
!>
!>    3. If UPPER='T', the upper triangle of A is set to random values
!>         out of distribution DIST.
!>
!>    4. If SIM='T', A is multiplied on the left by a random matrix
!>         X, whose singular values are specified by DS, MODES, and
!>         CONDS, and on the right by X inverse.
!>
!>    5. If KL < N-1, the lower bandwidth is reduced to KL using
!>         Householder transformations.  If KU < N-1, the upper
!>         bandwidth is reduced to KU.
!>
!>    6. If ANORM is not negative, the matrix is scaled to have
!>         maximum-element-norm ANORM.
!>
!>    (Note: since the matrix cannot be reduced beyond Hessenberg form,
!>     no packing options are available.)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The number of columns (or rows) of A. Not modified.
!> \endverbatim
!>
!> \param[in] DIST
!> \verbatim
!>          DIST is CHARACTER*1
!>           On entry, DIST specifies the type of distribution to be used
!>           to generate the random eigen-/singular values, and for the
!>           upper triangle (see UPPER).
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
!>           exit, and can be used in the next call to DLATME
!>           to continue the same random number sequence.
!>           Changed on exit.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension ( N )
!>           This array is used to specify the eigenvalues of A.  If
!>           MODE=0, then D is assumed to contain the eigenvalues (but
!>           see the description of EI), otherwise they will be
!>           computed according to MODE, COND, DMAX, and RSIGN and
!>           placed in D.
!>           Modified if MODE is nonzero.
!> \endverbatim
!>
!> \param[in] MODE
!> \verbatim
!>          MODE is INTEGER
!>           On entry this describes how the eigenvalues are to
!>           be specified:
!>           MODE = 0 means use D (with EI) as input
!>           MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND
!>           MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND
!>           MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))
!>           MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)
!>           MODE = 5 sets D to random numbers in the range
!>                    ( 1/COND , 1 ) such that their logarithms
!>                    are uniformly distributed.  Each odd-even pair
!>                    of elements will be either used as two real
!>                    eigenvalues or as the real and imaginary part
!>                    of a complex conjugate pair of eigenvalues;
!>                    the choice of which is done is random, with
!>                    50-50 probability, for each pair.
!>           MODE = 6 set D to random numbers from same distribution
!>                    as the rest of the matrix.
!>           MODE < 0 has the same meaning as ABS(MODE), except that
!>              the order of the elements of D is reversed.
!>           Thus if MODE is between 1 and 4, D has entries ranging
!>              from 1 to 1/COND, if between -1 and -4, D has entries
!>              ranging from 1/COND to 1,
!>           Not modified.
!> \endverbatim
!>
!> \param[in] COND
!> \verbatim
!>          COND is DOUBLE PRECISION
!>           On entry, this is used as described under MODE above.
!>           If used, it must be >= 1. Not modified.
!> \endverbatim
!>
!> \param[in] DMAX
!> \verbatim
!>          DMAX is DOUBLE PRECISION
!>           If MODE is neither -6, 0 nor 6, the contents of D, as
!>           computed according to MODE and COND, will be scaled by
!>           DMAX / max(abs(D(i))).  Note that DMAX need not be
!>           positive: if DMAX is negative (or zero), D will be
!>           scaled by a negative number (or zero).
!>           Not modified.
!> \endverbatim
!>
!> \param[in] EI
!> \verbatim
!>          EI is CHARACTER*1 array, dimension ( N )
!>           If MODE is 0, and EI(1) is not ' ' (space character),
!>           this array specifies which elements of D (on input) are
!>           real eigenvalues and which are the real and imaginary parts
!>           of a complex conjugate pair of eigenvalues.  The elements
!>           of EI may then only have the values 'R' and 'I'.  If
!>           EI(j)='R' and EI(j+1)='I', then the j-th eigenvalue is
!>           CMPLX( D(j) , D(j+1) ), and the (j+1)-th is the complex
!>           conjugate thereof.  If EI(j)=EI(j+1)='R', then the j-th
!>           eigenvalue is D(j) (i.e., real).  EI(1) may not be 'I',
!>           nor may two adjacent elements of EI both have the value 'I'.
!>           If MODE is not 0, then EI is ignored.  If MODE is 0 and
!>           EI(1)=' ', then the eigenvalues will all be real.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] RSIGN
!> \verbatim
!>          RSIGN is CHARACTER*1
!>           If MODE is not 0, 6, or -6, and RSIGN='T', then the
!>           elements of D, as computed according to MODE and COND, will
!>           be multiplied by a random sign (+1 or -1).  If RSIGN='F',
!>           they will not be.  RSIGN may only have the values 'T' or
!>           'F'.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] UPPER
!> \verbatim
!>          UPPER is CHARACTER*1
!>           If UPPER='T', then the elements of A above the diagonal
!>           (and above the 2x2 diagonal blocks, if A has complex
!>           eigenvalues) will be set to random numbers out of DIST.
!>           If UPPER='F', they will not.  UPPER may only have the
!>           values 'T' or 'F'.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] SIM
!> \verbatim
!>          SIM is CHARACTER*1
!>           If SIM='T', then A will be operated on by a "similarity
!>           transform", i.e., multiplied on the left by a matrix X and
!>           on the right by X inverse.  X = U S V, where U and V are
!>           random unitary matrices and S is a (diagonal) matrix of
!>           singular values specified by DS, MODES, and CONDS.  If
!>           SIM='F', then A will not be transformed.
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] DS
!> \verbatim
!>          DS is DOUBLE PRECISION array, dimension ( N )
!>           This array is used to specify the singular values of X,
!>           in the same way that D specifies the eigenvalues of A.
!>           If MODE=0, the DS contains the singular values, which
!>           may not be zero.
!>           Modified if MODE is nonzero.
!> \endverbatim
!>
!> \param[in] MODES
!> \verbatim
!>          MODES is INTEGER
!> \endverbatim
!>
!> \param[in] CONDS
!> \verbatim
!>          CONDS is DOUBLE PRECISION
!>           Same as MODE and COND, but for specifying the diagonal
!>           of S.  MODES=-6 and +6 are not allowed (since they would
!>           result in randomly ill-conditioned eigenvalues.)
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>           This specifies the lower bandwidth of the  matrix.  KL=1
!>           specifies upper Hessenberg form.  If KL is at least N-1,
!>           then A will have full lower bandwidth.  KL must be at
!>           least 1.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>           This specifies the upper bandwidth of the  matrix.  KU=1
!>           specifies lower Hessenberg form.  If KU is at least N-1,
!>           then A will have full upper bandwidth; if KU and KL
!>           are both at least N-1, then A will be dense.  Only one of
!>           KU and KL may be less than N-1.  KU must be at least 1.
!>           Not modified.
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is DOUBLE PRECISION
!>           If ANORM is not negative, then A will be scaled by a non-
!>           negative real number to make the maximum-element-norm of A
!>           to be ANORM.
!>           Not modified.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension ( LDA, N )
!>           On exit A is the desired test matrix.
!>           Modified.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           LDA specifies the first dimension of A as declared in the
!>           calling program.  LDA must be at least N.
!>           Not modified.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension ( 3*N )
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
!>            -1 => N negative
!>            -2 => DIST illegal string
!>            -5 => MODE not in range -6 to 6
!>            -6 => COND less than 1.0, and MODE neither -6, 0 nor 6
!>            -8 => EI(1) is not ' ' or 'R', EI(j) is not 'R' or 'I', or
!>                  two adjacent elements of EI are 'I'.
!>            -9 => RSIGN is not 'T' or 'F'
!>           -10 => UPPER is not 'T' or 'F'
!>           -11 => SIM   is not 'T' or 'F'
!>           -12 => MODES=0 and DS has a zero singular value.
!>           -13 => MODES is not in the range -5 to 5.
!>           -14 => MODES is nonzero and CONDS is less than 1.
!>           -15 => KL is less than 1.
!>           -16 => KU is less than 1, or KL and KU are both less than
!>                  N-1.
!>           -19 => LDA is less than N.
!>            1  => Error return from DLATM1 (computing D)
!>            2  => Cannot scale to DMAX (max. eigenvalue is 0)
!>            3  => Error return from DLATM1 (computing DS)
!>            4  => Error return from DLARGE
!>            5  => Zero singular value from DLATM1.
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
!> \ingroup double_matgen
!
!  =====================================================================
      SUBROUTINE DLATME(N,Dist,Iseed,D,Mode,Cond,Dmax,Ei,Rsign,Upper,   &
     &                  Sim,Ds,Modes,Conds,Kl,Ku,Anorm,A,Lda,Work,Info)
      IMPLICIT NONE
!*--DLATME333
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Dist , Rsign , Sim , Upper
      INTEGER Info , Kl , Ku , Lda , Mode , Modes , N
      DOUBLE PRECISION Anorm , Cond , Conds , Dmax
!     ..
!     .. Array Arguments ..
      CHARACTER Ei(*)
      INTEGER Iseed(4)
      DOUBLE PRECISION A(Lda,*) , D(*) , Ds(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=1.0D0/2.0D0)
!     ..
!     .. Local Scalars ..
      LOGICAL badei , bads , useei
      INTEGER i , ic , icols , idist , iinfo , ir , irows , irsign ,    &
     &        isim , iupper , j , jc , jcr , jr
      DOUBLE PRECISION alpha , tau , temp , xnorms
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION tempa(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLANGE , DLARAN
      EXTERNAL LSAME , DLANGE , DLARAN
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DGEMV , DGER , DLARFG , DLARGE , DLARNV ,        &
     &         DLASET , DLATM1 , DSCAL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MOD
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
      IF ( N==0 ) RETURN
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
!     Check EI
!
      useei = .TRUE.
      badei = .FALSE.
      IF ( LSAME(Ei(1),' ') .OR. Mode/=0 ) THEN
         useei = .FALSE.
      ELSEIF ( LSAME(Ei(1),'R') ) THEN
         DO j = 2 , N
            IF ( LSAME(Ei(j),'I') ) THEN
               IF ( LSAME(Ei(j-1),'I') ) badei = .TRUE.
            ELSE
               IF ( .NOT.LSAME(Ei(j),'R') ) badei = .TRUE.
            ENDIF
         ENDDO
      ELSE
         badei = .TRUE.
      ENDIF
!
!     Decode RSIGN
!
      IF ( LSAME(Rsign,'T') ) THEN
         irsign = 1
      ELSEIF ( LSAME(Rsign,'F') ) THEN
         irsign = 0
      ELSE
         irsign = -1
      ENDIF
!
!     Decode UPPER
!
      IF ( LSAME(Upper,'T') ) THEN
         iupper = 1
      ELSEIF ( LSAME(Upper,'F') ) THEN
         iupper = 0
      ELSE
         iupper = -1
      ENDIF
!
!     Decode SIM
!
      IF ( LSAME(Sim,'T') ) THEN
         isim = 1
      ELSEIF ( LSAME(Sim,'F') ) THEN
         isim = 0
      ELSE
         isim = -1
      ENDIF
!
!     Check DS, if MODES=0 and ISIM=1
!
      bads = .FALSE.
      IF ( Modes==0 .AND. isim==1 ) THEN
         DO j = 1 , N
            IF ( Ds(j)==ZERO ) bads = .TRUE.
         ENDDO
      ENDIF
!
!     Set INFO if an error
!
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( idist==-1 ) THEN
         Info = -2
      ELSEIF ( ABS(Mode)>6 ) THEN
         Info = -5
      ELSEIF ( (Mode/=0 .AND. ABS(Mode)/=6) .AND. Cond<ONE ) THEN
         Info = -6
      ELSEIF ( badei ) THEN
         Info = -8
      ELSEIF ( irsign==-1 ) THEN
         Info = -9
      ELSEIF ( iupper==-1 ) THEN
         Info = -10
      ELSEIF ( isim==-1 ) THEN
         Info = -11
      ELSEIF ( bads ) THEN
         Info = -12
      ELSEIF ( isim==1 .AND. ABS(Modes)>5 ) THEN
         Info = -13
      ELSEIF ( isim==1 .AND. Modes/=0 .AND. Conds<ONE ) THEN
         Info = -14
      ELSEIF ( Kl<1 ) THEN
         Info = -15
      ELSEIF ( Ku<1 .OR. (Ku<N-1 .AND. Kl<N-1) ) THEN
         Info = -16
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -19
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DLATME',-Info)
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
!     2)      Set up diagonal of A
!
!             Compute D according to COND and MODE
!
      CALL DLATM1(Mode,Cond,irsign,idist,Iseed,D,N,iinfo)
      IF ( iinfo/=0 ) THEN
         Info = 1
         RETURN
      ENDIF
      IF ( Mode/=0 .AND. ABS(Mode)/=6 ) THEN
!
!        Scale by DMAX
!
         temp = ABS(D(1))
         DO i = 2 , N
            temp = MAX(temp,ABS(D(i)))
         ENDDO
!
         IF ( temp>ZERO ) THEN
            alpha = Dmax/temp
         ELSEIF ( Dmax/=ZERO ) THEN
            Info = 2
            RETURN
         ELSE
            alpha = ZERO
         ENDIF
!
         CALL DSCAL(N,alpha,D,1)
!
      ENDIF
!
      CALL DLASET('Full',N,N,ZERO,ZERO,A,Lda)
      CALL DCOPY(N,D,1,A,Lda+1)
!
!     Set up complex conjugate pairs
!
      IF ( Mode==0 ) THEN
         IF ( useei ) THEN
            DO j = 2 , N
               IF ( LSAME(Ei(j),'I') ) THEN
                  A(j-1,j) = A(j,j)
                  A(j,j-1) = -A(j,j)
                  A(j,j) = A(j-1,j-1)
               ENDIF
            ENDDO
         ENDIF
!
      ELSEIF ( ABS(Mode)==5 ) THEN
!
         DO j = 2 , N , 2
            IF ( DLARAN(Iseed)>HALF ) THEN
               A(j-1,j) = A(j,j)
               A(j,j-1) = -A(j,j)
               A(j,j) = A(j-1,j-1)
            ENDIF
         ENDDO
      ENDIF
!
!     3)      If UPPER='T', set upper triangle of A to random numbers.
!             (but don't modify the corners of 2x2 blocks.)
!
      IF ( iupper/=0 ) THEN
         DO jc = 2 , N
            IF ( A(jc-1,jc)/=ZERO ) THEN
               jr = jc - 2
            ELSE
               jr = jc - 1
            ENDIF
            CALL DLARNV(idist,Iseed,jr,A(1,jc))
         ENDDO
      ENDIF
!
!     4)      If SIM='T', apply similarity transformation.
!
!                                -1
!             Transform is  X A X  , where X = U S V, thus
!
!             it is  U S V A V' (1/S) U'
!
      IF ( isim/=0 ) THEN
!
!        Compute S (singular values of the eigenvector matrix)
!        according to CONDS and MODES
!
         CALL DLATM1(Modes,Conds,0,0,Iseed,Ds,N,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = 3
            RETURN
         ENDIF
!
!        Multiply by V and V'
!
         CALL DLARGE(N,A,Lda,Iseed,Work,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = 4
            RETURN
         ENDIF
!
!        Multiply by S and (1/S)
!
         DO j = 1 , N
            CALL DSCAL(N,Ds(j),A(j,1),Lda)
            IF ( Ds(j)/=ZERO ) THEN
               CALL DSCAL(N,ONE/Ds(j),A(1,j),1)
            ELSE
               Info = 5
               RETURN
            ENDIF
         ENDDO
!
!        Multiply by U and U'
!
         CALL DLARGE(N,A,Lda,Iseed,Work,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = 4
            RETURN
         ENDIF
      ENDIF
!
!     5)      Reduce the bandwidth.
!
      IF ( Kl<N-1 ) THEN
!
!        Reduce bandwidth -- kill column
!
         DO jcr = Kl + 1 , N - 1
            ic = jcr - Kl
            irows = N + 1 - jcr
            icols = N + Kl - jcr
!
            CALL DCOPY(irows,A(jcr,ic),1,Work,1)
            xnorms = Work(1)
            CALL DLARFG(irows,xnorms,Work(2),1,tau)
            Work(1) = ONE
!
            CALL DGEMV('T',irows,icols,ONE,A(jcr,ic+1),Lda,Work,1,ZERO, &
     &                 Work(irows+1),1)
            CALL DGER(irows,icols,-tau,Work,1,Work(irows+1),1,          &
     &                A(jcr,ic+1),Lda)
!
            CALL DGEMV('N',N,irows,ONE,A(1,jcr),Lda,Work,1,ZERO,        &
     &                 Work(irows+1),1)
            CALL DGER(N,irows,-tau,Work(irows+1),1,Work,1,A(1,jcr),Lda)
!
            A(jcr,ic) = xnorms
            CALL DLASET('Full',irows-1,1,ZERO,ZERO,A(jcr+1,ic),Lda)
         ENDDO
      ELSEIF ( Ku<N-1 ) THEN
!
!        Reduce upper bandwidth -- kill a row at a time.
!
         DO jcr = Ku + 1 , N - 1
            ir = jcr - Ku
            irows = N + Ku - jcr
            icols = N + 1 - jcr
!
            CALL DCOPY(icols,A(ir,jcr),Lda,Work,1)
            xnorms = Work(1)
            CALL DLARFG(icols,xnorms,Work(2),1,tau)
            Work(1) = ONE
!
            CALL DGEMV('N',irows,icols,ONE,A(ir+1,jcr),Lda,Work,1,ZERO, &
     &                 Work(icols+1),1)
            CALL DGER(irows,icols,-tau,Work(icols+1),1,Work,1,          &
     &                A(ir+1,jcr),Lda)
!
            CALL DGEMV('C',icols,N,ONE,A(jcr,1),Lda,Work,1,ZERO,        &
     &                 Work(icols+1),1)
            CALL DGER(icols,N,-tau,Work,1,Work(icols+1),1,A(jcr,1),Lda)
!
            A(ir,jcr) = xnorms
            CALL DLASET('Full',1,icols-1,ZERO,ZERO,A(ir,jcr+1),Lda)
         ENDDO
      ENDIF
!
!     Scale the matrix to have norm ANORM
!
      IF ( Anorm>=ZERO ) THEN
         temp = DLANGE('M',N,N,A,Lda,tempa)
         IF ( temp>ZERO ) THEN
            alpha = Anorm/temp
            DO j = 1 , N
               CALL DSCAL(N,alpha,A(1,j),1)
            ENDDO
         ENDIF
      ENDIF
!
!
!     End of DLATME
!
      END SUBROUTINE DLATME
