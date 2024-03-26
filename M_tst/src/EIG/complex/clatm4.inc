!*==clatm4.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b clatm4
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE clatm4( ITYPE, N, NZ1, NZ2, RSIGN, AMAGN, RCOND,
!                          TRIANG, IDIST, ISEED, A, LDA )
!
!       .. Scalar Arguments ..
!       LOGICAL            RSIGN
!       INTEGER            IDIST, ITYPE, LDA, N, NZ1, NZ2
!       REAL               AMAGN, RCOND, TRIANG
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> clatm4 generates basic square matrices, which may later be
!> multiplied by others in order to produce test matrices.  It is
!> intended mainly to be used to test the generalized eigenvalue
!> routines.
!>
!> It first generates the diagonal and (possibly) subdiagonal,
!> according to the value of ITYPE, NZ1, NZ2, RSIGN, AMAGN, and RCOND.
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
!>             the description of AMAGN and RSIGN.
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
!>          = 10: random numbers with distribution IDIST (see clarnd.)
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
!> \param[in] RSIGN
!> \verbatim
!>          RSIGN is LOGICAL
!>          = .TRUE.:  The diagonal and subdiagonal entries will be
!>                     multiplied by random numbers of magnitude 1.
!>          = .FALSE.: The diagonal and subdiagonal entries will be
!>                     left as they are (usually non-negative real.)
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
!>          RCOND.  RCOND must be between 0 and 1.
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
!>          On entry, DIST specifies the type of distribution to be used
!>          to generate a random matrix .
!>          = 1: real and imaginary parts each UNIFORM( 0, 1 )
!>          = 2: real and imaginary parts each UNIFORM( -1, 1 )
!>          = 3: real and imaginary parts each NORMAL( 0, 1 )
!>          = 4: complex number uniform in DISK( 0, 1 )
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry ISEED specifies the seed of the random number
!>          generator.  The values of ISEED are changed on exit, and can
!>          be used in the next call to clatm4 to continue the same
!>          random number sequence.
!>          Note: ISEED(4) should be odd, for the random number generator
!>          used at present.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, N)
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE clatm4(Itype,N,Nz1,Nz2,Rsign,Amagn,Rcond,Triang,Idist, &
     &                  Iseed,A,Lda)
      IMPLICIT NONE
!*--CLATM4175
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Rsign
      INTEGER Idist , Itype , Lda , N , Nz1 , Nz2
      REAL Amagn , Rcond , Triang
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      COMPLEX A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , isdb , isde , jc , jd , jr , k , kbeg , kend , klen
      REAL alpha
      COMPLEX ctemp
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
      EXTERNAL CLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CMPLX , EXP , LOG , MAX , MIN , MOD , REAL
!     ..
!     .. Executable Statements ..
!
      IF ( N<=0 ) RETURN
      CALL CLASET('Full',N,N,CZERO,CZERO,A,Lda)
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
               A(jd+1,jd) = CONE
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
               A(jd+1,jd) = CONE
            ENDDO
            isdb = 1
            isde = k
            DO jd = k + 2 , 2*k + 1
               A(jd,jd) = CONE
            ENDDO
         ELSEIF ( ABS(Itype)==4 ) THEN
!
!        abs(ITYPE) = 4: 1,...,k
!
            DO jd = kbeg , kend
               A(jd,jd) = CMPLX(jd-Nz1)
            ENDDO
         ELSEIF ( ABS(Itype)==5 ) THEN
!
!        abs(ITYPE) = 5: One large D value:
!
            DO jd = kbeg + 1 , kend
               A(jd,jd) = CMPLX(Rcond)
            ENDDO
            A(kbeg,kbeg) = CONE
         ELSEIF ( ABS(Itype)==6 ) THEN
!
!        abs(ITYPE) = 6: One small D value:
!
            DO jd = kbeg , kend - 1
               A(jd,jd) = CONE
            ENDDO
            A(kend,kend) = CMPLX(Rcond)
         ELSEIF ( ABS(Itype)==7 ) THEN
!
!        abs(ITYPE) = 7: Exponentially distributed D values:
!
            A(kbeg,kbeg) = CONE
            IF ( klen>1 ) THEN
               alpha = Rcond**(ONE/REAL(klen-1))
               DO i = 2 , klen
                  A(Nz1+i,Nz1+i) = CMPLX(alpha**REAL(i-1))
               ENDDO
            ENDIF
         ELSEIF ( ABS(Itype)==8 ) THEN
!
!        abs(ITYPE) = 8: Arithmetically distributed D values:
!
            A(kbeg,kbeg) = CONE
            IF ( klen>1 ) THEN
               alpha = (ONE-Rcond)/REAL(klen-1)
               DO i = 2 , klen
                  A(Nz1+i,Nz1+i) = CMPLX(REAL(klen-i)*alpha+Rcond)
               ENDDO
            ENDIF
         ELSEIF ( ABS(Itype)==9 ) THEN
!
!        abs(ITYPE) = 9: Randomly distributed D values on ( RCOND, 1):
!
            alpha = LOG(Rcond)
            DO jd = kbeg , kend
               A(jd,jd) = EXP(alpha*slaran(Iseed))
            ENDDO
         ELSEIF ( ABS(Itype)==10 ) THEN
!
!        abs(ITYPE) = 10: Randomly distributed D values from DIST
!
            DO jd = kbeg , kend
               A(jd,jd) = clarnd(Idist,Iseed)
            ENDDO
         ELSE
!
!        abs(ITYPE) = 1: Identity
!
            DO jd = 1 , N
               A(jd,jd) = CONE
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
!        If RSIGN = .TRUE., assign random signs to diagonal and
!        subdiagonal
!
         IF ( Rsign ) THEN
            DO jd = kbeg , kend
               IF ( REAL(A(jd,jd))/=ZERO ) THEN
                  ctemp = clarnd(3,Iseed)
                  ctemp = ctemp/ABS(ctemp)
                  A(jd,jd) = ctemp*REAL(A(jd,jd))
               ENDIF
            ENDDO
            DO jd = isdb , isde
               IF ( REAL(A(jd+1,jd))/=ZERO ) THEN
                  ctemp = clarnd(3,Iseed)
                  ctemp = ctemp/ABS(ctemp)
                  A(jd+1,jd) = ctemp*REAL(A(jd+1,jd))
               ENDIF
            ENDDO
         ENDIF
!
!        Reverse if ITYPE < 0
!
         IF ( Itype<0 ) THEN
            DO jd = kbeg , (kbeg+kend-1)/2
               ctemp = A(jd,jd)
               A(jd,jd) = A(kbeg+kend-jd,kbeg+kend-jd)
               A(kbeg+kend-jd,kbeg+kend-jd) = ctemp
            ENDDO
            DO jd = 1 , (N-1)/2
               ctemp = A(jd+1,jd)
               A(jd+1,jd) = A(N+1-jd,N-jd)
               A(N+1-jd,N-jd) = ctemp
            ENDDO
         ENDIF
!
      ENDIF
!
!     Fill in upper triangle
!
      IF ( Triang/=ZERO ) THEN
         DO jc = 2 , N
            DO jr = 1 , jc - 1
               A(jr,jc) = Triang*clarnd(Idist,Iseed)
            ENDDO
         ENDDO
      ENDIF
!
!
!     End of clatm4
!
      END SUBROUTINE clatm4
