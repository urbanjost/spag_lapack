!*==slatm4.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SLATM4
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLATM4( ITYPE, N, NZ1, NZ2, ISIGN, AMAGN, RCOND,
!                          TRIANG, IDIST, ISEED, A, LDA )
!
!       .. Scalar Arguments ..
!       INTEGER            IDIST, ISIGN, ITYPE, LDA, N, NZ1, NZ2
!       REAL               AMAGN, RCOND, TRIANG
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLATM4 generates basic square matrices, which may later be
!> multiplied by others in order to produce test matrices.  It is
!> intended mainly to be used to test the generalized eigenvalue
!> routines.
!>
!> It first generates the diagonal and (possibly) subdiagonal,
!> according to the value of ITYPE, NZ1, NZ2, ISIGN, AMAGN, and RCOND.
!> It then fills in the upper triangle with random numbers, if TRIANG is
!> non-zero.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          The "type" of matrix on the diagonal and sub-diagonal.
!>          If ITYPE < 0, then type abs(ITYPE) is generated and then
!>             swapped end for end (A(I,J) := A'(N-J,N-I).)  See also
!>             the description of AMAGN and ISIGN.
!>
!>          Special types:
!>          = 0:  the zero matrix.
!>          = 1:  the identity.
!>          = 2:  a transposed Jordan block.
!>          = 3:  If N is odd, then a k+1 x k+1 transposed Jordan block
!>                followed by a k x k identity block, where k=(N-1)/2.
!>                If N is even, then k=(N-2)/2, and a zero diagonal entry
!>                is tacked onto the end.
!>
!>          Diagonal types.  The diagonal consists of NZ1 zeros, then
!>             k=N-NZ1-NZ2 nonzeros.  The subdiagonal is zero.  ITYPE
!>             specifies the nonzero diagonal entries as follows:
!>          = 4:  1, ..., k
!>          = 5:  1, RCOND, ..., RCOND
!>          = 6:  1, ..., 1, RCOND
!>          = 7:  1, a, a^2, ..., a^(k-1)=RCOND
!>          = 8:  1, 1-d, 1-2*d, ..., 1-(k-1)*d=RCOND
!>          = 9:  random numbers chosen from (RCOND,1)
!>          = 10: random numbers with distribution IDIST (see SLARND.)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.
!> \endverbatim
!>
!> \param[in] NZ1
!> \verbatim
!>          NZ1 is INTEGER
!>          If abs(ITYPE) > 3, then the first NZ1 diagonal entries will
!>          be zero.
!> \endverbatim
!>
!> \param[in] NZ2
!> \verbatim
!>          NZ2 is INTEGER
!>          If abs(ITYPE) > 3, then the last NZ2 diagonal entries will
!>          be zero.
!> \endverbatim
!>
!> \param[in] ISIGN
!> \verbatim
!>          ISIGN is INTEGER
!>          = 0: The sign of the diagonal and subdiagonal entries will
!>               be left unchanged.
!>          = 1: The diagonal and subdiagonal entries will have their
!>               sign changed at random.
!>          = 2: If ITYPE is 2 or 3, then the same as ISIGN=1.
!>               Otherwise, with probability 0.5, odd-even pairs of
!>               diagonal entries A(2*j-1,2*j-1), A(2*j,2*j) will be
!>               converted to a 2x2 block by pre- and post-multiplying
!>               by distinct random orthogonal rotations.  The remaining
!>               diagonal entries will have their sign changed at random.
!> \endverbatim
!>
!> \param[in] AMAGN
!> \verbatim
!>          AMAGN is REAL
!>          The diagonal and subdiagonal entries will be multiplied by
!>          AMAGN.
!> \endverbatim
!>
!> \param[in] RCOND
!> \verbatim
!>          RCOND is REAL
!>          If abs(ITYPE) > 4, then the smallest diagonal entry will be
!>          entry will be RCOND.  RCOND must be between 0 and 1.
!> \endverbatim
!>
!> \param[in] TRIANG
!> \verbatim
!>          TRIANG is REAL
!>          The entries above the diagonal will be random numbers with
!>          magnitude bounded by TRIANG (i.e., random numbers multiplied
!>          by TRIANG.)
!> \endverbatim
!>
!> \param[in] IDIST
!> \verbatim
!>          IDIST is INTEGER
!>          Specifies the type of distribution to be used to generate a
!>          random matrix.
!>          = 1:  UNIFORM( 0, 1 )
!>          = 2:  UNIFORM( -1, 1 )
!>          = 3:  NORMAL ( 0, 1 )
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator.  The values of ISEED are changed on exit, and can
!>          be used in the next call to SLATM4 to continue the same
!>          random number sequence.
!>          Note: ISEED(4) should be odd, for the random number generator
!>          used at present.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (LDA, N)
!>          Array to be computed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          Leading dimension of A.  Must be at least 1 and at least N.
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE SLATM4(Itype,N,Nz1,Nz2,Isign,Amagn,Rcond,Triang,Idist, &
     &                  Iseed,A,Lda)
      IMPLICIT NONE
!*--SLATM4179
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Idist , Isign , Itype , Lda , N , Nz1 , Nz2
      REAL Amagn , Rcond , Triang
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      REAL A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,TWO=2.0E0)
      REAL HALF
      PARAMETER (HALF=ONE/TWO)
!     ..
!     .. Local Scalars ..
      INTEGER i , ioff , isdb , isde , jc , jd , jr , k , kbeg , kend , &
     &        klen
      REAL alpha , cl , cr , safmin , sl , sr , sv1 , sv2 , temp
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLARAN , SLARND
      EXTERNAL SLAMCH , SLARAN , SLARND
!     ..
!     .. External Subroutines ..
      EXTERNAL SLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , EXP , LOG , MAX , MIN , MOD , REAL , SQRT
!     ..
!     .. Executable Statements ..
!
      IF ( N<=0 ) RETURN
      CALL SLASET('Full',N,N,ZERO,ZERO,A,Lda)
!
!     Insure a correct ISEED
!
      IF ( MOD(Iseed(4),2)/=1 ) Iseed(4) = Iseed(4) + 1
!
!     Compute diagonal and subdiagonal according to ITYPE, NZ1, NZ2,
!     and RCOND
!
      IF ( Itype/=0 ) THEN
         IF ( ABS(Itype)>=4 ) THEN
            kbeg = MAX(1,MIN(N,Nz1+1))
            kend = MAX(kbeg,MIN(N,N-Nz2))
            klen = kend + 1 - kbeg
         ELSE
            kbeg = 1
            kend = N
            klen = N
         ENDIF
         isdb = 1
         isde = 0
         IF ( ABS(Itype)==2 ) THEN
!
!        abs(ITYPE) = 2: Transposed Jordan block
!
            DO jd = 1 , N - 1
               A(jd+1,jd) = ONE
            ENDDO
            isdb = 1
            isde = N - 1
         ELSEIF ( ABS(Itype)==3 ) THEN
!
!        abs(ITYPE) = 3: Transposed Jordan block, followed by the
!                        identity.
!
            k = (N-1)/2
            DO jd = 1 , k
               A(jd+1,jd) = ONE
            ENDDO
            isdb = 1
            isde = k
            DO jd = k + 2 , 2*k + 1
               A(jd,jd) = ONE
            ENDDO
         ELSEIF ( ABS(Itype)==4 ) THEN
!
!        abs(ITYPE) = 4: 1,...,k
!
            DO jd = kbeg , kend
               A(jd,jd) = REAL(jd-Nz1)
            ENDDO
         ELSEIF ( ABS(Itype)==5 ) THEN
!
!        abs(ITYPE) = 5: One large D value:
!
            DO jd = kbeg + 1 , kend
               A(jd,jd) = Rcond
            ENDDO
            A(kbeg,kbeg) = ONE
         ELSEIF ( ABS(Itype)==6 ) THEN
!
!        abs(ITYPE) = 6: One small D value:
!
            DO jd = kbeg , kend - 1
               A(jd,jd) = ONE
            ENDDO
            A(kend,kend) = Rcond
         ELSEIF ( ABS(Itype)==7 ) THEN
!
!        abs(ITYPE) = 7: Exponentially distributed D values:
!
            A(kbeg,kbeg) = ONE
            IF ( klen>1 ) THEN
               alpha = Rcond**(ONE/REAL(klen-1))
               DO i = 2 , klen
                  A(Nz1+i,Nz1+i) = alpha**REAL(i-1)
               ENDDO
            ENDIF
         ELSEIF ( ABS(Itype)==8 ) THEN
!
!        abs(ITYPE) = 8: Arithmetically distributed D values:
!
            A(kbeg,kbeg) = ONE
            IF ( klen>1 ) THEN
               alpha = (ONE-Rcond)/REAL(klen-1)
               DO i = 2 , klen
                  A(Nz1+i,Nz1+i) = REAL(klen-i)*alpha + Rcond
               ENDDO
            ENDIF
         ELSEIF ( ABS(Itype)==9 ) THEN
!
!        abs(ITYPE) = 9: Randomly distributed D values on ( RCOND, 1):
!
            alpha = LOG(Rcond)
            DO jd = kbeg , kend
               A(jd,jd) = EXP(alpha*SLARAN(Iseed))
            ENDDO
         ELSEIF ( ABS(Itype)==10 ) THEN
!
!        abs(ITYPE) = 10: Randomly distributed D values from DIST
!
            DO jd = kbeg , kend
               A(jd,jd) = SLARND(Idist,Iseed)
            ENDDO
         ELSE
!
!        abs(ITYPE) = 1: Identity
!
            DO jd = 1 , N
               A(jd,jd) = ONE
            ENDDO
         ENDIF
!
!
!        Scale by AMAGN
!
         DO jd = kbeg , kend
            A(jd,jd) = Amagn*REAL(A(jd,jd))
         ENDDO
         DO jd = isdb , isde
            A(jd+1,jd) = Amagn*REAL(A(jd+1,jd))
         ENDDO
!
!        If ISIGN = 1 or 2, assign random signs to diagonal and
!        subdiagonal
!
         IF ( Isign>0 ) THEN
            DO jd = kbeg , kend
               IF ( REAL(A(jd,jd))/=ZERO ) THEN
                  IF ( SLARAN(Iseed)>HALF ) A(jd,jd) = -A(jd,jd)
               ENDIF
            ENDDO
            DO jd = isdb , isde
               IF ( REAL(A(jd+1,jd))/=ZERO ) THEN
                  IF ( SLARAN(Iseed)>HALF ) A(jd+1,jd) = -A(jd+1,jd)
               ENDIF
            ENDDO
         ENDIF
!
!        Reverse if ITYPE < 0
!
         IF ( Itype<0 ) THEN
            DO jd = kbeg , (kbeg+kend-1)/2
               temp = A(jd,jd)
               A(jd,jd) = A(kbeg+kend-jd,kbeg+kend-jd)
               A(kbeg+kend-jd,kbeg+kend-jd) = temp
            ENDDO
            DO jd = 1 , (N-1)/2
               temp = A(jd+1,jd)
               A(jd+1,jd) = A(N+1-jd,N-jd)
               A(N+1-jd,N-jd) = temp
            ENDDO
         ENDIF
!
!        If ISIGN = 2, and no subdiagonals already, then apply
!        random rotations to make 2x2 blocks.
!
         IF ( Isign==2 .AND. Itype/=2 .AND. Itype/=3 ) THEN
            safmin = SLAMCH('S')
            DO jd = kbeg , kend - 1 , 2
               IF ( SLARAN(Iseed)>HALF ) THEN
!
!                 Rotation on left.
!
                  cl = TWO*SLARAN(Iseed) - ONE
                  sl = TWO*SLARAN(Iseed) - ONE
                  temp = ONE/MAX(safmin,SQRT(cl**2+sl**2))
                  cl = cl*temp
                  sl = sl*temp
!
!                 Rotation on right.
!
                  cr = TWO*SLARAN(Iseed) - ONE
                  sr = TWO*SLARAN(Iseed) - ONE
                  temp = ONE/MAX(safmin,SQRT(cr**2+sr**2))
                  cr = cr*temp
                  sr = sr*temp
!
!                 Apply
!
                  sv1 = A(jd,jd)
                  sv2 = A(jd+1,jd+1)
                  A(jd,jd) = cl*cr*sv1 + sl*sr*sv2
                  A(jd+1,jd) = -sl*cr*sv1 + cl*sr*sv2
                  A(jd,jd+1) = -cl*sr*sv1 + sl*cr*sv2
                  A(jd+1,jd+1) = sl*sr*sv1 + cl*cr*sv2
               ENDIF
            ENDDO
         ENDIF
!
      ENDIF
!
!     Fill in upper triangle (except for 2x2 blocks)
!
      IF ( Triang/=ZERO ) THEN
         IF ( Isign/=2 .OR. Itype==2 .OR. Itype==3 ) THEN
            ioff = 1
         ELSE
            ioff = 2
            DO jr = 1 , N - 1
               IF ( A(jr+1,jr)==ZERO ) A(jr,jr+1)                       &
     &              = Triang*SLARND(Idist,Iseed)
            ENDDO
         ENDIF
!
         DO jc = 2 , N
            DO jr = 1 , jc - ioff
               A(jr,jc) = Triang*SLARND(Idist,Iseed)
            ENDDO
         ENDDO
      ENDIF
!
!
!     End of SLATM4
!
      END SUBROUTINE SLATM4
