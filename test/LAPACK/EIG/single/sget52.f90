!*==sget52.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b SGET52
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGET52( LEFT, N, A, LDA, B, LDB, E, LDE, ALPHAR,
!                          ALPHAI, BETA, WORK, RESULT )
!
!       .. Scalar Arguments ..
!       LOGICAL            LEFT
!       INTEGER            LDA, LDB, LDE, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
!      $                   B( LDB, * ), BETA( * ), E( LDE, * ),
!      $                   RESULT( 2 ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGET52  does an eigenvector check for the generalized eigenvalue
!> problem.
!>
!> The basic test for right eigenvectors is:
!>
!>                           | b(j) A E(j) -  a(j) B E(j) |
!>         RESULT(1) = max   -------------------------------
!>                      j    n ulp max( |b(j) A|, |a(j) B| )
!>
!> using the 1-norm.  Here, a(j)/b(j) = w is the j-th generalized
!> eigenvalue of A - w B, or, equivalently, b(j)/a(j) = m is the j-th
!> generalized eigenvalue of m A - B.
!>
!> For real eigenvalues, the test is straightforward.  For complex
!> eigenvalues, E(j) and a(j) are complex, represented by
!> Er(j) + i*Ei(j) and ar(j) + i*ai(j), resp., so the test for that
!> eigenvector becomes
!>
!>                 max( |Wr|, |Wi| )
!>     --------------------------------------------
!>     n ulp max( |b(j) A|, (|ar(j)|+|ai(j)|) |B| )
!>
!> where
!>
!>     Wr = b(j) A Er(j) - ar(j) B Er(j) + ai(j) B Ei(j)
!>
!>     Wi = b(j) A Ei(j) - ai(j) B Er(j) - ar(j) B Ei(j)
!>
!>                         T   T  _
!> For left eigenvectors, A , B , a, and b  are used.
!>
!> SGET52 also tests the normalization of E.  Each eigenvector is
!> supposed to be normalized so that the maximum "absolute value"
!> of its elements is 1, where in this case, "absolute value"
!> of a complex value x is  |Re(x)| + |Im(x)| ; let us call this
!> maximum "absolute value" norm of a vector v  M(v).
!> if a(j)=b(j)=0, then the eigenvector is set to be the jth coordinate
!> vector.  The normalization test is:
!>
!>         RESULT(2) =      max       | M(v(j)) - 1 | / ( n ulp )
!>                    eigenvectors v(j)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] LEFT
!> \verbatim
!>          LEFT is LOGICAL
!>          =.TRUE.:  The eigenvectors in the columns of E are assumed
!>                    to be *left* eigenvectors.
!>          =.FALSE.: The eigenvectors in the columns of E are assumed
!>                    to be *right* eigenvectors.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrices.  If it is zero, SGET52 does
!>          nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA, N)
!>          The matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 1
!>          and at least N.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB, N)
!>          The matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of B.  It must be at least 1
!>          and at least N.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (LDE, N)
!>          The matrix of eigenvectors.  It must be O( 1 ).  Complex
!>          eigenvalues and eigenvectors always come in pairs, the
!>          eigenvalue and its conjugate being stored in adjacent
!>          elements of ALPHAR, ALPHAI, and BETA.  Thus, if a(j)/b(j)
!>          and a(j+1)/b(j+1) are a complex conjugate pair of
!>          generalized eigenvalues, then E(,j) contains the real part
!>          of the eigenvector and E(,j+1) contains the imaginary part.
!>          Note that whether E(,j) is a real eigenvector or part of a
!>          complex one is specified by whether ALPHAI(j) is zero or not.
!> \endverbatim
!>
!> \param[in] LDE
!> \verbatim
!>          LDE is INTEGER
!>          The leading dimension of E.  It must be at least 1 and at
!>          least N.
!> \endverbatim
!>
!> \param[in] ALPHAR
!> \verbatim
!>          ALPHAR is REAL array, dimension (N)
!>          The real parts of the values a(j) as described above, which,
!>          along with b(j), define the generalized eigenvalues.
!>          Complex eigenvalues always come in complex conjugate pairs
!>          a(j)/b(j) and a(j+1)/b(j+1), which are stored in adjacent
!>          elements in ALPHAR, ALPHAI, and BETA.  Thus, if the j-th
!>          and (j+1)-st eigenvalues form a pair, ALPHAR(j+1)/BETA(j+1)
!>          is assumed to be equal to ALPHAR(j)/BETA(j).
!> \endverbatim
!>
!> \param[in] ALPHAI
!> \verbatim
!>          ALPHAI is REAL array, dimension (N)
!>          The imaginary parts of the values a(j) as described above,
!>          which, along with b(j), define the generalized eigenvalues.
!>          If ALPHAI(j)=0, then the eigenvalue is real, otherwise it
!>          is part of a complex conjugate pair.  Complex eigenvalues
!>          always come in complex conjugate pairs a(j)/b(j) and
!>          a(j+1)/b(j+1), which are stored in adjacent elements in
!>          ALPHAR, ALPHAI, and BETA.  Thus, if the j-th and (j+1)-st
!>          eigenvalues form a pair, ALPHAI(j+1)/BETA(j+1) is assumed to
!>          be equal to  -ALPHAI(j)/BETA(j).  Also, nonzero values in
!>          ALPHAI are assumed to always come in adjacent pairs.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is REAL array, dimension (N)
!>          The values b(j) as described above, which, along with a(j),
!>          define the generalized eigenvalues.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N**2+N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
!>          The values computed by the test described above.  If A E or
!>          B E is likely to overflow, then RESULT(1:2) is set to
!>          10 / ulp.
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
      SUBROUTINE SGET52(Left,N,A,Lda,B,Ldb,E,Lde,Alphar,Alphai,Beta,    &
     &                  Work,Result)
      IMPLICIT NONE
!*--SGET52203
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Left
      INTEGER Lda , Ldb , Lde , N
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , Alphai(*) , Alphar(*) , B(Ldb,*) , Beta(*) ,      &
     &     E(Lde,*) , Result(2) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TEN
      PARAMETER (ZERO=0.0,ONE=1.0,TEN=10.0)
!     ..
!     .. Local Scalars ..
      LOGICAL ilcplx
      CHARACTER normab , trans
      INTEGER j , jvec
      REAL abmax , acoef , alfmax , anorm , bcoefi , bcoefr , betmax ,  &
     &     bnorm , enorm , enrmer , errnrm , safmax , safmin , salfi ,  &
     &     salfr , sbeta , scale , temp1 , ulp
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE
      EXTERNAL SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , REAL
!     ..
!     .. Executable Statements ..
!
      Result(1) = ZERO
      Result(2) = ZERO
      IF ( N<=0 ) RETURN
!
      safmin = SLAMCH('Safe minimum')
      safmax = ONE/safmin
      ulp = SLAMCH('Epsilon')*SLAMCH('Base')
!
      IF ( Left ) THEN
         trans = 'T'
         normab = 'I'
      ELSE
         trans = 'N'
         normab = 'O'
      ENDIF
!
!     Norm of A, B, and E:
!
      anorm = MAX(SLANGE(normab,N,N,A,Lda,Work),safmin)
      bnorm = MAX(SLANGE(normab,N,N,B,Ldb,Work),safmin)
      enorm = MAX(SLANGE('O',N,N,E,Lde,Work),ulp)
      alfmax = safmax/MAX(ONE,bnorm)
      betmax = safmax/MAX(ONE,anorm)
!
!     Compute error matrix.
!     Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B| |b(i) A| )
!
      ilcplx = .FALSE.
      DO jvec = 1 , N
         IF ( ilcplx ) THEN
!
!           2nd Eigenvalue/-vector of pair -- do nothing
!
            ilcplx = .FALSE.
         ELSE
            salfr = Alphar(jvec)
            salfi = Alphai(jvec)
            sbeta = Beta(jvec)
            IF ( salfi==ZERO ) THEN
!
!              Real eigenvalue and -vector
!
               abmax = MAX(ABS(salfr),ABS(sbeta))
               IF ( ABS(salfr)>alfmax .OR. ABS(sbeta)>betmax .OR.       &
     &              abmax<ONE ) THEN
                  scale = ONE/MAX(abmax,safmin)
                  salfr = scale*salfr
                  sbeta = scale*sbeta
               ENDIF
               scale = ONE/MAX(ABS(salfr)*bnorm,ABS(sbeta)*anorm,safmin)
               acoef = scale*sbeta
               bcoefr = scale*salfr
               CALL SGEMV(trans,N,N,acoef,A,Lda,E(1,jvec),1,ZERO,       &
     &                    Work(N*(jvec-1)+1),1)
               CALL SGEMV(trans,N,N,-bcoefr,B,Lda,E(1,jvec),1,ONE,      &
     &                    Work(N*(jvec-1)+1),1)
            ELSE
!
!              Complex conjugate pair
!
               ilcplx = .TRUE.
               IF ( jvec==N ) THEN
                  Result(1) = TEN/ulp
                  RETURN
               ENDIF
               abmax = MAX(ABS(salfr)+ABS(salfi),ABS(sbeta))
               IF ( ABS(salfr)+ABS(salfi)>alfmax .OR. ABS(sbeta)        &
     &              >betmax .OR. abmax<ONE ) THEN
                  scale = ONE/MAX(abmax,safmin)
                  salfr = scale*salfr
                  salfi = scale*salfi
                  sbeta = scale*sbeta
               ENDIF
               scale = ONE/MAX((ABS(salfr)+ABS(salfi))*bnorm,ABS(sbeta) &
     &                 *anorm,safmin)
               acoef = scale*sbeta
               bcoefr = scale*salfr
               bcoefi = scale*salfi
               IF ( Left ) bcoefi = -bcoefi
!
               CALL SGEMV(trans,N,N,acoef,A,Lda,E(1,jvec),1,ZERO,       &
     &                    Work(N*(jvec-1)+1),1)
               CALL SGEMV(trans,N,N,-bcoefr,B,Lda,E(1,jvec),1,ONE,      &
     &                    Work(N*(jvec-1)+1),1)
               CALL SGEMV(trans,N,N,bcoefi,B,Lda,E(1,jvec+1),1,ONE,     &
     &                    Work(N*(jvec-1)+1),1)
!
               CALL SGEMV(trans,N,N,acoef,A,Lda,E(1,jvec+1),1,ZERO,     &
     &                    Work(N*jvec+1),1)
               CALL SGEMV(trans,N,N,-bcoefi,B,Lda,E(1,jvec),1,ONE,      &
     &                    Work(N*jvec+1),1)
               CALL SGEMV(trans,N,N,-bcoefr,B,Lda,E(1,jvec+1),1,ONE,    &
     &                    Work(N*jvec+1),1)
            ENDIF
         ENDIF
      ENDDO
!
      errnrm = SLANGE('One',N,N,Work,N,Work(N**2+1))/enorm
!
!     Compute RESULT(1)
!
      Result(1) = errnrm/ulp
!
!     Normalization of E:
!
      enrmer = ZERO
      ilcplx = .FALSE.
      DO jvec = 1 , N
         IF ( ilcplx ) THEN
            ilcplx = .FALSE.
         ELSE
            temp1 = ZERO
            IF ( Alphai(jvec)==ZERO ) THEN
               DO j = 1 , N
                  temp1 = MAX(temp1,ABS(E(j,jvec)))
               ENDDO
               enrmer = MAX(enrmer,temp1-ONE)
            ELSE
               ilcplx = .TRUE.
               DO j = 1 , N
                  temp1 = MAX(temp1,ABS(E(j,jvec))+ABS(E(j,jvec+1)))
               ENDDO
               enrmer = MAX(enrmer,temp1-ONE)
            ENDIF
         ENDIF
      ENDDO
!
!     Compute RESULT(2) : the normalization error in E.
!
      Result(2) = enrmer/(REAL(N)*ulp)
!
!
!     End of SGET52
!
      END SUBROUTINE SGET52
