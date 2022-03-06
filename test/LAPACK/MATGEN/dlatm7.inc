!*==dlatm7.f90  processed by SPAG 7.51RB at 17:44 on  4 Mar 2022
!> \brief \b dlatm7
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dlatm7( MODE, COND, IRSIGN, IDIST, ISEED, D, N,
!                          RANK, INFO )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   COND
!       INTEGER            IDIST, INFO, IRSIGN, MODE, N, RANK
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * )
!       INTEGER            ISEED( 4 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    dlatm7 computes the entries of D as specified by MODE
!>    COND and IRSIGN. IDIST and ISEED determine the generation
!>    of random numbers. dlatm7 is called by DLATMT to generate
!>    random test matrices.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  MODE   - INTEGER
!>           On entry describes how D is to be computed:
!>           MODE = 0 means do not change D.
!>
!>           MODE = 1 sets D(1)=1 and D(2:RANK)=1.0/COND
!>           MODE = 2 sets D(1:RANK-1)=1 and D(RANK)=1.0/COND
!>           MODE = 3 sets D(I)=COND**(-(I-1)/(RANK-1)) I=1:RANK
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
!>           Not modified.
!>
!>  COND   - DOUBLE PRECISION
!>           On entry, used as described under MODE above.
!>           If used, it must be >= 1. Not modified.
!>
!>  IRSIGN - INTEGER
!>           On entry, if MODE neither -6, 0 nor 6, determines sign of
!>           entries of D
!>           0 => leave entries of D unchanged
!>           1 => multiply each entry of D by 1 or -1 with probability .5
!>
!>  IDIST  - CHARACTER*1
!>           On entry, IDIST specifies the type of distribution to be
!>           used to generate a random matrix .
!>           1 => UNIFORM( 0, 1 )
!>           2 => UNIFORM( -1, 1 )
!>           3 => NORMAL( 0, 1 )
!>           Not modified.
!>
!>  ISEED  - INTEGER array, dimension ( 4 )
!>           On entry ISEED specifies the seed of the random number
!>           generator. The random number generator uses a
!>           linear congruential sequence limited to small
!>           integers, and so should produce machine independent
!>           random numbers. The values of ISEED are changed on
!>           exit, and can be used in the next call to dlatm7
!>           to continue the same random number sequence.
!>           Changed on exit.
!>
!>  D      - DOUBLE PRECISION array, dimension ( MIN( M , N ) )
!>           Array to be computed according to MODE, COND and IRSIGN.
!>           May be changed on exit if MODE is nonzero.
!>
!>  N      - INTEGER
!>           Number of entries of D. Not modified.
!>
!>  RANK   - INTEGER
!>           The rank of matrix to be generated for modes 1,2,3 only.
!>           D( RANK+1:N ) = 0.
!>           Not modified.
!>
!>  INFO   - INTEGER
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
      SUBROUTINE dlatm7(Mode,Cond,Irsign,Idist,Iseed,D,N,Rank,Info)
      IMPLICIT NONE
!*--dlatm7125
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Cond
      INTEGER Idist , Info , Irsign , Mode , N , Rank
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION D(*)
      INTEGER Iseed(4)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION alpha , temp
      INTEGER i
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
         CALL XERBLA('dlatm7',-Info)
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
            DO i = 1 , Rank - 1
               D(i) = ONE
            ENDDO
            DO i = Rank + 1 , N
               D(i) = ZERO
            ENDDO
            D(Rank) = ONE/Cond
         ELSEIF ( ABS(Mode)==3 ) THEN
!
!        Exponentially distributed D values:
!
            D(1) = ONE
            IF ( N>1 .AND. Rank>1 ) THEN
               alpha = Cond**(-ONE/DBLE(Rank-1))
               DO i = 2 , Rank
                  D(i) = alpha**(i-1)
               ENDDO
               DO i = Rank + 1 , N
                  D(i) = ZERO
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
            DO i = 2 , Rank
               D(i) = ONE/Cond
            ENDDO
            DO i = Rank + 1 , N
               D(i) = ZERO
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
!     End of dlatm7
!
      END SUBROUTINE dlatm7
