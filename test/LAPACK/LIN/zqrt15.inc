!*==zqrt15.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zqrt15
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zqrt15( SCALE, RKSEL, M, N, NRHS, A, LDA, B, LDB, S,
!                          RANK, NORMA, NORMB, ISEED, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LWORK, M, N, NRHS, RANK, RKSEL, SCALE
!       DOUBLE PRECISION   NORMA, NORMB
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   S( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zqrt15 generates a matrix with full or deficient rank and of various
!> norms.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SCALE
!> \verbatim
!>          SCALE is INTEGER
!>          SCALE = 1: normally scaled matrix
!>          SCALE = 2: matrix scaled up
!>          SCALE = 3: matrix scaled down
!> \endverbatim
!>
!> \param[in] RKSEL
!> \verbatim
!>          RKSEL is INTEGER
!>          RKSEL = 1: full rank matrix
!>          RKSEL = 2: rank-deficient matrix
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of A.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of B.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB, NRHS)
!>          A matrix that is in the range space of matrix A.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension MIN(M,N)
!>          Singular values of A.
!> \endverbatim
!>
!> \param[out] RANK
!> \verbatim
!>          RANK is INTEGER
!>          number of nonzero singular values of A.
!> \endverbatim
!>
!> \param[out] NORMA
!> \verbatim
!>          NORMA is DOUBLE PRECISION
!>          one-norm norm of A.
!> \endverbatim
!>
!> \param[out] NORMB
!> \verbatim
!>          NORMB is DOUBLE PRECISION
!>          one-norm norm of B.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is integer array, dimension (4)
!>          seed for random number generator.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          length of work space required.
!>          LWORK >= MAX(M+MIN(M,N),NRHS*MIN(M,N),2*N+M)
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zqrt15(Scale,Rksel,M,N,Nrhs,A,Lda,B,Ldb,S,Rank,Norma,  &
     &                  Normb,Iseed,Work,Lwork)
      IMPLICIT NONE
!*--ZQRT15153
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldb , Lwork , M , N , Nrhs , Rank , Rksel , Scale
      DOUBLE PRECISION Norma , Normb
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      DOUBLE PRECISION S(*)
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TWO , SVMIN
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0,TWO=2.0D+0,SVMIN=0.1D+0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER info , j , mn
      DOUBLE PRECISION bignum , eps , smlnum , temp
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION dummy(1)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DASUM , DLAMCH , DLARND , DZNRM2 , ZLANGE
EXTERNAL    DASUM,   DLAMCH,   DLARND,   DZNRM2,   ZLANGE
!     ..
!     .. External Subroutines ..
EXTERNAL    DLABAD,   DLASCL,   ZDSCAL,   ZGEMM
EXTERNAL    ZLARF,   ZLARNV,   ZLAROR,   ZLASCL,   ZLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DCMPLX , MAX , MIN
!     ..
!     .. Executable Statements ..
!
      mn = MIN(M,N)
      IF ( Lwork<MAX(M+mn,mn*Nrhs,2*N+M) ) THEN
         CALL xerbla('zqrt15',16)
         RETURN
      ENDIF
!
      smlnum = DLAMCH('Safe minimum')
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
      eps = DLAMCH('Epsilon')
      smlnum = (smlnum/eps)/eps
      bignum = ONE/smlnum
!
!     Determine rank and (unscaled) singular values
!
      IF ( Rksel==1 ) THEN
         Rank = mn
      ELSEIF ( Rksel==2 ) THEN
         Rank = (3*mn)/4
         DO j = Rank + 1 , mn
            S(j) = ZERO
         ENDDO
      ELSE
         CALL xerbla('zqrt15',2)
      ENDIF
!
      IF ( Rank>0 ) THEN
!
!        Nontrivial case
!
         S(1) = ONE
         DO j = 2 , Rank
            DO
               temp = DLARND(1,Iseed)
               IF ( temp<=SVMIN ) CYCLE
               S(j) = ABS(temp)
               EXIT
            ENDDO
         ENDDO
         CALL dlaord('Decreasing',Rank,S,1)
!
!        Generate 'rank' columns of a random orthogonal matrix in A
!
         CALL ZLARNV(2,Iseed,M,Work)
         CALL ZDSCAL(M,ONE/DZNRM2(M,Work,1),Work,1)
         CALL ZLASET('Full',M,Rank,CZERO,CONE,A,Lda)
         CALL ZLARF('Left',M,Rank,Work,1,DCMPLX(TWO),A,Lda,Work(M+1))
!
!        workspace used: m+mn
!
!        Generate consistent rhs in the range space of A
!
         CALL ZLARNV(2,Iseed,Rank*Nrhs,Work)
         CALL ZGEMM('No transpose','No transpose',M,Nrhs,Rank,CONE,A,   &
     &              Lda,Work,Rank,CZERO,B,Ldb)
!
!        work space used: <= mn *nrhs
!
!        generate (unscaled) matrix A
!
         DO j = 1 , Rank
            CALL ZDSCAL(M,S(j),A(1,j),1)
         ENDDO
         IF ( Rank<N ) CALL ZLASET('Full',M,N-Rank,CZERO,CZERO,         &
     &                             A(1,Rank+1),Lda)
         CALL ZLAROR('Right','No initialization',M,N,A,Lda,Iseed,Work,  &
     &               info)
!
      ELSE
!
!        work space used 2*n+m
!
!        Generate null matrix and rhs
!
         DO j = 1 , mn
            S(j) = ZERO
         ENDDO
         CALL ZLASET('Full',M,N,CZERO,CZERO,A,Lda)
         CALL ZLASET('Full',M,Nrhs,CZERO,CZERO,B,Ldb)
!
      ENDIF
!
!     Scale the matrix
!
      IF ( Scale/=1 ) THEN
         Norma = ZLANGE('Max',M,N,A,Lda,dummy)
         IF ( Norma/=ZERO ) THEN
            IF ( Scale==2 ) THEN
!
!              matrix scaled up
!
               CALL ZLASCL('General',0,0,Norma,bignum,M,N,A,Lda,info)
               CALL DLASCL('General',0,0,Norma,bignum,mn,1,S,mn,info)
               CALL ZLASCL('General',0,0,Norma,bignum,M,Nrhs,B,Ldb,info)
            ELSEIF ( Scale==3 ) THEN
!
!              matrix scaled down
!
               CALL ZLASCL('General',0,0,Norma,smlnum,M,N,A,Lda,info)
               CALL DLASCL('General',0,0,Norma,smlnum,mn,1,S,mn,info)
               CALL ZLASCL('General',0,0,Norma,smlnum,M,Nrhs,B,Ldb,info)
            ELSE
               CALL xerbla('zqrt15',1)
               RETURN
            ENDIF
         ENDIF
      ENDIF
!
      Norma = DASUM(mn,S,1)
      Normb = ZLANGE('One-norm',M,Nrhs,B,Ldb,dummy)
!
!
!     End of zqrt15
!
      END SUBROUTINE zqrt15
