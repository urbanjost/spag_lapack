!*==zget22.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZGET22
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGET22( TRANSA, TRANSE, TRANSW, N, A, LDA, E, LDE, W,
!                          WORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANSA, TRANSE, TRANSW
!       INTEGER            LDA, LDE, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RESULT( 2 ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), E( LDE, * ), W( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGET22 does an eigenvector check.
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
!> vector.  The max-norm of a complex n-vector x in this case is the
!> maximum of |re(x(i)| + |im(x(i)| over i = 1, ..., n.
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
!>          = 'C':  Conjugate transpose
!> \endverbatim
!>
!> \param[in] TRANSE
!> \verbatim
!>          TRANSE is CHARACTER*1
!>          Specifies whether or not E is transposed.
!>          = 'N':  No transpose, eigenvectors are in columns of E
!>          = 'T':  Transpose, eigenvectors are in rows of E
!>          = 'C':  Conjugate transpose, eigenvectors are in rows of E
!> \endverbatim
!>
!> \param[in] TRANSW
!> \verbatim
!>          TRANSW is CHARACTER*1
!>          Specifies whether or not W is transposed.
!>          = 'N':  No transpose
!>          = 'T':  Transpose, same as TRANSW = 'N'
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
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
!>          E is COMPLEX*16 array, dimension (LDE,N)
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
!> \param[in] W
!> \verbatim
!>          W is COMPLEX*16 array, dimension (N)
!>          The eigenvalues of A.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (2)
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZGET22(Transa,Transe,Transw,N,A,Lda,E,Lde,W,Work,Rwork,&
     &                  Result)
      IMPLICIT NONE
!*--ZGET22147
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
      DOUBLE PRECISION Result(2) , Rwork(*)
      COMPLEX*16 A(Lda,*) , E(Lde,*) , W(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      CHARACTER norma , norme
      INTEGER itrnse , itrnsw , j , jcol , joff , jrow , jvec
      DOUBLE PRECISION anorm , enorm , enrmax , enrmin , errnrm ,       &
     &                 temp1 , ulp , unfl
      COMPLEX*16 wtemp
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL LSAME , DLAMCH , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMM , ZLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DCONJG , DIMAG , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Initialize RESULT (in case N=0)
!
      Result(1) = ZERO
      Result(2) = ZERO
      IF ( N<=0 ) RETURN
!
      unfl = DLAMCH('Safe minimum')
      ulp = DLAMCH('Precision')
!
      itrnse = 0
      itrnsw = 0
      norma = 'O'
      norme = 'O'
!
      IF ( LSAME(Transa,'T') .OR. LSAME(Transa,'C') ) norma = 'I'
!
      IF ( LSAME(Transe,'T') ) THEN
         itrnse = 1
         norme = 'I'
      ELSEIF ( LSAME(Transe,'C') ) THEN
         itrnse = 2
         norme = 'I'
      ENDIF
!
      IF ( LSAME(Transw,'C') ) itrnsw = 1
!
!     Normalization of E:
!
      enrmin = ONE/ulp
      enrmax = ZERO
      IF ( itrnse==0 ) THEN
         DO jvec = 1 , N
            temp1 = ZERO
            DO j = 1 , N
               temp1 = MAX(temp1,ABS(DBLE(E(j,jvec)))                   &
     &                 +ABS(DIMAG(E(j,jvec))))
            ENDDO
            enrmin = MIN(enrmin,temp1)
            enrmax = MAX(enrmax,temp1)
         ENDDO
      ELSE
         DO jvec = 1 , N
            Rwork(jvec) = ZERO
         ENDDO
!
         DO j = 1 , N
            DO jvec = 1 , N
               Rwork(jvec) = MAX(Rwork(jvec),ABS(DBLE(E(jvec,j)))+ABS(  &
     &                       DIMAG(E(jvec,j))))
            ENDDO
         ENDDO
!
         DO jvec = 1 , N
            enrmin = MIN(enrmin,Rwork(jvec))
            enrmax = MAX(enrmax,Rwork(jvec))
         ENDDO
      ENDIF
!
!     Norm of A:
!
      anorm = MAX(ZLANGE(norma,N,N,A,Lda,Rwork),unfl)
!
!     Norm of E:
!
      enorm = MAX(ZLANGE(norme,N,N,E,Lde,Rwork),ulp)
!
!     Norm of error:
!
!     Error =  AE - EW
!
      CALL ZLASET('Full',N,N,CZERO,CZERO,Work,N)
!
      joff = 0
      DO jcol = 1 , N
         IF ( itrnsw==0 ) THEN
            wtemp = W(jcol)
         ELSE
            wtemp = DCONJG(W(jcol))
         ENDIF
!
         IF ( itrnse==0 ) THEN
            DO jrow = 1 , N
               Work(joff+jrow) = E(jrow,jcol)*wtemp
            ENDDO
         ELSEIF ( itrnse==1 ) THEN
            DO jrow = 1 , N
               Work(joff+jrow) = E(jcol,jrow)*wtemp
            ENDDO
         ELSE
            DO jrow = 1 , N
               Work(joff+jrow) = DCONJG(E(jcol,jrow))*wtemp
            ENDDO
         ENDIF
         joff = joff + N
      ENDDO
!
      CALL ZGEMM(Transa,Transe,N,N,N,CONE,A,Lda,E,Lde,-CONE,Work,N)
!
      errnrm = ZLANGE('One',N,N,Work,N,Rwork)/enorm
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
      Result(2) = MAX(ABS(enrmax-ONE),ABS(enrmin-ONE))/(DBLE(N)*ulp)
!
!
!     End of ZGET22
!
      END SUBROUTINE ZGET22
