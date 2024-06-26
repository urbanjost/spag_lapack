!*==cqrt15.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cqrt15
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cqrt15( SCALE, RKSEL, M, N, NRHS, A, LDA, B, LDB, S,
!                          RANK, NORMA, NORMB, ISEED, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LWORK, M, N, NRHS, RANK, RKSEL, SCALE
!       REAL               NORMA, NORMB
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               S( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cqrt15 generates a matrix with full or deficient rank and of various
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          B is COMPLEX array, dimension (LDB, NRHS)
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
!>          S is REAL array, dimension MIN(M,N)
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
!>          NORMA is REAL
!>          one-norm norm of A.
!> \endverbatim
!>
!> \param[out] NORMB
!> \verbatim
!>          NORMB is REAL
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
!>          WORK is COMPLEX array, dimension (LWORK)
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cqrt15(Scale,Rksel,M,N,Nrhs,A,Lda,B,Ldb,S,Rank,Norma,  &
     &                  Normb,Iseed,Work,Lwork)
      IMPLICIT NONE
!*--CQRT15153
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldb , Lwork , M , N , Nrhs , Rank , Rksel , Scale
      REAL Norma , Normb
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      REAL S(*)
      COMPLEX A(Lda,*) , B(Ldb,*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO , SVMIN
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0,TWO=2.0E+0,SVMIN=0.1E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER info , j , mn
      REAL bignum , eps , smlnum , temp
!     ..
!     .. Local Arrays ..
      REAL dummy(1)
!     ..
!     .. External Functions ..
      REAL CLANGE , SASUM , SCNRM2 , SLAMCH
EXTERNAL    CLANGE,   SASUM,   SCNRM2,   SLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    CGEMM,   CLARF,   CLARNV,   CLASCL,   CLASET,   CSSCAL,   SLABAD,   SLASCL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CMPLX , MAX , MIN
!     ..
!     .. Executable Statements ..
!
      mn = MIN(M,N)
      IF ( Lwork<MAX(M+mn,mn*Nrhs,2*N+M) ) THEN
         CALL xerbla('cqrt15',16)
         RETURN
      ENDIF
!
      smlnum = SLAMCH('Safe minimum')
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
      eps = SLAMCH('Epsilon')
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
         CALL xerbla('cqrt15',2)
      ENDIF
!
      IF ( Rank>0 ) THEN
!
!        Nontrivial case
!
         S(1) = ONE
         DO j = 2 , Rank
            DO
               temp = SLARND(1,Iseed)
               IF ( temp<=SVMIN ) CYCLE
               S(j) = ABS(temp)
               EXIT
            ENDDO
         ENDDO
         CALL slaord('Decreasing',Rank,S,1)
!
!        Generate 'rank' columns of a random orthogonal matrix in A
!
         CALL CLARNV(2,Iseed,M,Work)
         CALL CSSCAL(M,ONE/SCNRM2(M,Work,1),Work,1)
         CALL CLASET('Full',M,Rank,CZERO,CONE,A,Lda)
         CALL CLARF('Left',M,Rank,Work,1,CMPLX(TWO),A,Lda,Work(M+1))
!
!        workspace used: m+mn
!
!        Generate consistent rhs in the range space of A
!
         CALL CLARNV(2,Iseed,Rank*Nrhs,Work)
         CALL CGEMM('No transpose','No transpose',M,Nrhs,Rank,CONE,A,   &
     &              Lda,Work,Rank,CZERO,B,Ldb)
!
!        work space used: <= mn *nrhs
!
!        generate (unscaled) matrix A
!
         DO j = 1 , Rank
            CALL CSSCAL(M,S(j),A(1,j),1)
         ENDDO
         IF ( Rank<N ) CALL CLASET('Full',M,N-Rank,CZERO,CZERO,         &
     &                             A(1,Rank+1),Lda)
         CALL CLAROR('Right','No initialization',M,N,A,Lda,Iseed,Work,  &
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
         CALL CLASET('Full',M,N,CZERO,CZERO,A,Lda)
         CALL CLASET('Full',M,Nrhs,CZERO,CZERO,B,Ldb)
!
      ENDIF
!
!     Scale the matrix
!
      IF ( Scale/=1 ) THEN
         Norma = CLANGE('Max',M,N,A,Lda,dummy)
         IF ( Norma/=ZERO ) THEN
            IF ( Scale==2 ) THEN
!
!              matrix scaled up
!
               CALL CLASCL('General',0,0,Norma,bignum,M,N,A,Lda,info)
               CALL SLASCL('General',0,0,Norma,bignum,mn,1,S,mn,info)
               CALL CLASCL('General',0,0,Norma,bignum,M,Nrhs,B,Ldb,info)
            ELSEIF ( Scale==3 ) THEN
!
!              matrix scaled down
!
               CALL CLASCL('General',0,0,Norma,smlnum,M,N,A,Lda,info)
               CALL SLASCL('General',0,0,Norma,smlnum,mn,1,S,mn,info)
               CALL CLASCL('General',0,0,Norma,smlnum,M,Nrhs,B,Ldb,info)
            ELSE
               CALL xerbla('cqrt15',1)
               RETURN
            ENDIF
         ENDIF
      ENDIF
!
      Norma = SASUM(mn,S,1)
      Normb = CLANGE('One-norm',M,Nrhs,B,Ldb,dummy)
!
!
!     End of cqrt15
!
      END SUBROUTINE cqrt15
