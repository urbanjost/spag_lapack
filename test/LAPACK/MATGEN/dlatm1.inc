!*==dlatm1.f90  processed by SPAG 7.51RB at 17:44 on  4 Mar 2022
!> \brief \b dlatm1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dlatm1( MODE, COND, IRSIGN, IDIST, ISEED, D, N, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IDIST, INFO, IRSIGN, MODE, N
!       DOUBLE PRECISION   COND
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   D( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    dlatm1 computes the entries of D(1..N) as specified by
!>    MODE, COND and IRSIGN. IDIST and ISEED determine the generation
!>    of random numbers. dlatm1 is called by DLATMR to generate
!>    random test matrices for LAPACK programs.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] MODE
!> \verbatim
!>          MODE is INTEGER
!>           On entry describes how D is to be computed:
!>           MODE = 0 means do not change D.
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
!>          COND is DOUBLE PRECISION
!>           On entry, used as described under MODE above.
!>           If used, it must be >= 1. Not modified.
!> \endverbatim
!>
!> \param[in] IRSIGN
!> \verbatim
!>          IRSIGN is INTEGER
!>           On entry, if MODE neither -6, 0 nor 6, determines sign of
!>           entries of D
!>           0 => leave entries of D unchanged
!>           1 => multiply each entry of D by 1 or -1 with probability .5
!> \endverbatim
!>
!> \param[in] IDIST
!> \verbatim
!>          IDIST is INTEGER
!>           On entry, IDIST specifies the type of distribution to be
!>           used to generate a random matrix .
!>           1 => UNIFORM( 0, 1 )
!>           2 => UNIFORM( -1, 1 )
!>           3 => NORMAL( 0, 1 )
!>           Not modified.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension ( 4 )
!>           On entry ISEED specifies the seed of the random number
!>           generator. The random number generator uses a
!>           linear congruential sequence limited to small
!>           integers, and so should produce machine independent
!>           random numbers. The values of ISEED are changed on
!>           exit, and can be used in the next call to dlatm1
!>           to continue the same random number sequence.
!>           Changed on exit.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension ( N )
!>           Array to be computed according to MODE, COND and IRSIGN.
!>           May be changed on exit if MODE is nonzero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           Number of entries of D. Not modified.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>            0  => normal termination
!>           -1  => if MODE not in range -6 to 6
!>           -2  => if MODE neither -6, 0 nor 6, and
!>                  IRSIGN neither 0 nor 1
!>           -3  => if MODE neither -6, 0 nor 6 and COND less than 1
!>           -4  => if MODE equals 6 or -6 and IDIST not in range 1 to 3
!>           -7  => if N negative
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
      SUBROUTINE dlatm1(Mode,Cond,Irsign,Idist,Iseed,D,N,Info)
      IMPLICIT NONE
!*--dlatm1139
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Idist , Info , Irsign , Mode , N
      DOUBLE PRECISION Cond
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      DOUBLE PRECISION D(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      DOUBLE PRECISION alpha , temp
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
      EXTERNAL DLARNV , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , EXP , LOG
!     ..
!     .. Executable Statements ..
!
!     Decode and Test the input parameters. Initialize flags & seed.
!
      Info = 0
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Set INFO if an error
!
      IF ( Mode<-6 .OR. Mode>6 ) THEN
         Info = -1
      ELSEIF ( (Mode/=-6 .AND. Mode/=0 .AND. Mode/=6) .AND.             &
     &         (Irsign/=0 .AND. Irsign/=1) ) THEN
         Info = -2
      ELSEIF ( (Mode/=-6 .AND. Mode/=0 .AND. Mode/=6) .AND. Cond<ONE )  &
     &         THEN
         Info = -3
      ELSEIF ( (Mode==6 .OR. Mode==-6) .AND. (Idist<1 .OR. Idist>3) )   &
     &         THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -7
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('dlatm1',-Info)
         RETURN
      ENDIF
!
!     Compute D according to COND and MODE
!
      IF ( Mode/=0 ) THEN
         IF ( ABS(Mode)==2 ) THEN
!
!        One small D value:
!
            DO i = 1 , N
               D(i) = ONE
            ENDDO
            D(N) = ONE/Cond
         ELSEIF ( ABS(Mode)==3 ) THEN
!
!        Exponentially distributed D values:
!
            D(1) = ONE
            IF ( N>1 ) THEN
               alpha = Cond**(-ONE/DBLE(N-1))
               DO i = 2 , N
                  D(i) = alpha**(i-1)
               ENDDO
            ENDIF
         ELSEIF ( ABS(Mode)==4 ) THEN
!
!        Arithmetically distributed D values:
!
            D(1) = ONE
            IF ( N>1 ) THEN
               temp = ONE/Cond
               alpha = (ONE-temp)/DBLE(N-1)
               DO i = 2 , N
                  D(i) = DBLE(N-i)*alpha + temp
               ENDDO
            ENDIF
         ELSEIF ( ABS(Mode)==5 ) THEN
!
!        Randomly distributed D values on ( 1/COND , 1):
!
            alpha = LOG(ONE/Cond)
            DO i = 1 , N
               D(i) = EXP(alpha*dlaran(Iseed))
            ENDDO
         ELSEIF ( ABS(Mode)==6 ) THEN
!
!        Randomly distributed D values from IDIST
!
            CALL DLARNV(Idist,Iseed,N,D)
         ELSE
!
!        One large D value:
!
            DO i = 1 , N
               D(i) = ONE/Cond
            ENDDO
            D(1) = ONE
         ENDIF
!
!
!        If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
!        random signs to D
!
         IF ( (Mode/=-6 .AND. Mode/=0 .AND. Mode/=6) .AND. Irsign==1 )  &
     &        THEN
            DO i = 1 , N
               temp = dlaran(Iseed)
               IF ( temp>HALF ) D(i) = -D(i)
            ENDDO
         ENDIF
!
!        Reverse if MODE < 0
!
         IF ( Mode<0 ) THEN
            DO i = 1 , N/2
               temp = D(i)
               D(i) = D(N+1-i)
               D(N+1-i) = temp
            ENDDO
         ENDIF
!
      ENDIF
!
!
!     End of dlatm1
!
      END SUBROUTINE dlatm1
