!*==sget22.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SGET22
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGET22( TRANSA, TRANSE, TRANSW, N, A, LDA, E, LDE, WR,
!                          WI, WORK, RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANSA, TRANSE, TRANSW
!       INTEGER            LDA, LDE, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), E( LDE, * ), RESULT( 2 ), WI( * ),
!      $                   WORK( * ), WR( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGET22 does an eigenvector check.
!>
!> The basic test is:
!>
!>    RESULT(1) = | A E  -  E W | / ( |A| |E| ulp )
!>
!> using the 1-norm.  It also tests the normalization of E:
!>
!>    RESULT(2) = max | m-norm(E(j)) - 1 | / ( n ulp )
!>                 j
!>
!> where E(j) is the j-th eigenvector, and m-norm is the max-norm of a
!> vector.  If an eigenvector is complex, as determined from WI(j)
!> nonzero, then the max-norm of the vector ( er + i*ei ) is the maximum
!> of
!>    |er(1)| + |ei(1)|, ... , |er(n)| + |ei(n)|
!>
!> W is a block diagonal matrix, with a 1 by 1 block for each real
!> eigenvalue and a 2 by 2 block for each complex conjugate pair.
!> If eigenvalues j and j+1 are a complex conjugate pair, so that
!> WR(j) = WR(j+1) = wr and WI(j) = - WI(j+1) = wi, then the 2 by 2
!> block corresponding to the pair will be:
!>
!>    (  wr  wi  )
!>    ( -wi  wr  )
!>
!> Such a block multiplying an n by 2 matrix ( ur ui ) on the right
!> will be the same as multiplying  ur + i*ui  by  wr + i*wi.
!>
!> To handle various schemes for storage of left eigenvectors, there are
!> options to use A-transpose instead of A, E-transpose instead of E,
!> and/or W-transpose instead of W.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>          Specifies whether or not A is transposed.
!>          = 'N':  No transpose
!>          = 'T':  Transpose
!>          = 'C':  Conjugate transpose (= Transpose)
!> \endverbatim
!>
!> \param[in] TRANSE
!> \verbatim
!>          TRANSE is CHARACTER*1
!>          Specifies whether or not E is transposed.
!>          = 'N':  No transpose, eigenvectors are in columns of E
!>          = 'T':  Transpose, eigenvectors are in rows of E
!>          = 'C':  Conjugate transpose (= Transpose)
!> \endverbatim
!>
!> \param[in] TRANSW
!> \verbatim
!>          TRANSW is CHARACTER*1
!>          Specifies whether or not W is transposed.
!>          = 'N':  No transpose
!>          = 'T':  Transpose, use -WI(j) instead of WI(j)
!>          = 'C':  Conjugate transpose, use -WI(j) instead of WI(j)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The matrix whose eigenvectors are in E.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (LDE,N)
!>          The matrix of eigenvectors. If TRANSE = 'N', the eigenvectors
!>          are stored in the columns of E, if TRANSE = 'T' or 'C', the
!>          eigenvectors are stored in the rows of E.
!> \endverbatim
!>
!> \param[in] LDE
!> \verbatim
!>          LDE is INTEGER
!>          The leading dimension of the array E.  LDE >= max(1,N).
!> \endverbatim
!>
!> \param[in] WR
!> \verbatim
!>          WR is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[in] WI
!> \verbatim
!>          WI is REAL array, dimension (N)
!>
!>          The real and imaginary parts of the eigenvalues of A.
!>          Purely real eigenvalues are indicated by WI(j) = 0.
!>          Complex conjugate pairs are indicated by WR(j)=WR(j+1) and
!>          WI(j) = - WI(j+1) non-zero; the real part is assumed to be
!>          stored in the j-th row/column and the imaginary part in
!>          the (j+1)-th row/column.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N*(N+1))
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
!>          RESULT(1) = | A E  -  E W | / ( |A| |E| ulp )
!>          RESULT(2) = max | m-norm(E(j)) - 1 | / ( n ulp )
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
      SUBROUTINE SGET22(Transa,Transe,Transw,N,A,Lda,E,Lde,Wr,Wi,Work,  &
     &                  Result)
      IMPLICIT NONE
!*--SGET22171
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Transa , Transe , Transw
      INTEGER Lda , Lde , N
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , E(Lde,*) , Result(2) , Wi(*) , Work(*) , Wr(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0,ONE=1.0)
!     ..
!     .. Local Scalars ..
      CHARACTER norma , norme
      INTEGER iecol , ierow , ince , ipair , itrnse , j , jcol , jvec
      REAL anorm , enorm , enrmax , enrmin , errnrm , temp1 , ulp , unfl
!     ..
!     .. Local Arrays ..
      REAL wmat(2,2)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLAMCH , SLANGE
      EXTERNAL LSAME , SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SAXPY , SGEMM , SLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
!     Initialize RESULT (in case N=0)
!
      Result(1) = ZERO
      Result(2) = ZERO
      IF ( N<=0 ) RETURN
!
      unfl = SLAMCH('Safe minimum')
      ulp = SLAMCH('Precision')
!
      itrnse = 0
      ince = 1
      norma = 'O'
      norme = 'O'
!
      IF ( LSAME(Transa,'T') .OR. LSAME(Transa,'C') ) norma = 'I'
      IF ( LSAME(Transe,'T') .OR. LSAME(Transe,'C') ) THEN
         norme = 'I'
         itrnse = 1
         ince = Lde
      ENDIF
!
!     Check normalization of E
!
      enrmin = ONE/ulp
      enrmax = ZERO
      IF ( itrnse==0 ) THEN
!
!        Eigenvectors are column vectors.
!
         ipair = 0
         DO jvec = 1 , N
            temp1 = ZERO
            IF ( ipair==0 .AND. jvec<N .AND. Wi(jvec)/=ZERO ) ipair = 1
            IF ( ipair==1 ) THEN
!
!              Complex eigenvector
!
               DO j = 1 , N
                  temp1 = MAX(temp1,ABS(E(j,jvec))+ABS(E(j,jvec+1)))
               ENDDO
               enrmin = MIN(enrmin,temp1)
               enrmax = MAX(enrmax,temp1)
               ipair = 2
            ELSEIF ( ipair==2 ) THEN
               ipair = 0
            ELSE
!
!              Real eigenvector
!
               DO j = 1 , N
                  temp1 = MAX(temp1,ABS(E(j,jvec)))
               ENDDO
               enrmin = MIN(enrmin,temp1)
               enrmax = MAX(enrmax,temp1)
               ipair = 0
            ENDIF
         ENDDO
!
      ELSE
!
!        Eigenvectors are row vectors.
!
         DO jvec = 1 , N
            Work(jvec) = ZERO
         ENDDO
!
         DO j = 1 , N
            ipair = 0
            DO jvec = 1 , N
               IF ( ipair==0 .AND. jvec<N .AND. Wi(jvec)/=ZERO )        &
     &              ipair = 1
               IF ( ipair==1 ) THEN
                  Work(jvec) = MAX(Work(jvec),ABS(E(j,jvec))+ABS(E(j,   &
     &                         jvec+1)))
                  Work(jvec+1) = Work(jvec)
               ELSEIF ( ipair==2 ) THEN
                  ipair = 0
               ELSE
                  Work(jvec) = MAX(Work(jvec),ABS(E(j,jvec)))
                  ipair = 0
               ENDIF
            ENDDO
         ENDDO
!
         DO jvec = 1 , N
            enrmin = MIN(enrmin,Work(jvec))
            enrmax = MAX(enrmax,Work(jvec))
         ENDDO
      ENDIF
!
!     Norm of A:
!
      anorm = MAX(SLANGE(norma,N,N,A,Lda,Work),unfl)
!
!     Norm of E:
!
      enorm = MAX(SLANGE(norme,N,N,E,Lde,Work),ulp)
!
!     Norm of error:
!
!     Error =  AE - EW
!
      CALL SLASET('Full',N,N,ZERO,ZERO,Work,N)
!
      ipair = 0
      ierow = 1
      iecol = 1
!
      DO jcol = 1 , N
         IF ( itrnse==1 ) THEN
            ierow = jcol
         ELSE
            iecol = jcol
         ENDIF
!
         IF ( ipair==0 .AND. Wi(jcol)/=ZERO ) ipair = 1
!
         IF ( ipair==1 ) THEN
            wmat(1,1) = Wr(jcol)
            wmat(2,1) = -Wi(jcol)
            wmat(1,2) = Wi(jcol)
            wmat(2,2) = Wr(jcol)
            CALL SGEMM(Transe,Transw,N,2,2,ONE,E(ierow,iecol),Lde,wmat, &
     &                 2,ZERO,Work(N*(jcol-1)+1),N)
            ipair = 2
         ELSEIF ( ipair==2 ) THEN
            ipair = 0
!
         ELSE
!
            CALL SAXPY(N,Wr(jcol),E(ierow,iecol),ince,Work(N*(jcol-1)+1)&
     &                 ,1)
            ipair = 0
         ENDIF
!
      ENDDO
!
      CALL SGEMM(Transa,Transe,N,N,N,ONE,A,Lda,E,Lde,-ONE,Work,N)
!
      errnrm = SLANGE('One',N,N,Work,N,Work(N*N+1))/enorm
!
!     Compute RESULT(1) (avoiding under/overflow)
!
      IF ( anorm>errnrm ) THEN
         Result(1) = (errnrm/anorm)/ulp
      ELSEIF ( anorm<ONE ) THEN
         Result(1) = (MIN(errnrm,anorm)/anorm)/ulp
      ELSE
         Result(1) = MIN(errnrm/anorm,ONE)/ulp
      ENDIF
!
!     Compute RESULT(2) : the normalization error in E.
!
      Result(2) = MAX(ABS(enrmax-ONE),ABS(enrmin-ONE))/(REAL(N)*ulp)
!
!
!     End of SGET22
!
      END SUBROUTINE SGET22
