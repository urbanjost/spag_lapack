!*==dsyequb.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DSYEQUB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSYEQUB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyequb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyequb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyequb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, N
!       DOUBLE PRECISION   AMAX, SCOND
!       CHARACTER          UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), S( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSYEQUB computes row and column scalings intended to equilibrate a
!> symmetric matrix A (with respect to the Euclidean norm) and reduce
!> its condition number. The scale factors S are computed by the BIN
!> algorithm (see references) so that the scaled matrix B with elements
!> B(i,j) = S(i)*A(i,j)*S(j) has a condition number within a factor N of
!> the smallest possible condition number over all possible diagonal
!> scalings.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A. N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The N-by-N symmetric matrix whose scaling factors are to be
!>          computed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (N)
!>          If INFO = 0, S contains the scale factors for A.
!> \endverbatim
!>
!> \param[out] SCOND
!> \verbatim
!>          SCOND is DOUBLE PRECISION
!>          If INFO = 0, S contains the ratio of the smallest S(i) to
!>          the largest S(i). If SCOND >= 0.1 and AMAX is neither too
!>          large nor too small, it is not worth scaling by S.
!> \endverbatim
!>
!> \param[out] AMAX
!> \verbatim
!>          AMAX is DOUBLE PRECISION
!>          Largest absolute value of any matrix element. If AMAX is
!>          very close to overflow or very close to underflow, the
!>          matrix should be scaled.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the i-th diagonal element is nonpositive.
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
!> \date November 2017
!
!> \ingroup doubleSYcomputational
!
!> \par References:
!  ================
!>
!>  Livne, O.E. and Golub, G.H., "Scaling by Binormalization", \n
!>  Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004. \n
!>  DOI 10.1023/B:NUMA.0000016606.32820.69 \n
!>  Tech report version: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.3.1679
!>
!  =====================================================================
      SUBROUTINE DSYEQUB(Uplo,N,A,Lda,S,Scond,Amax,Work,Info)
      IMPLICIT NONE
!*--DSYEQUB135
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , N
      DOUBLE PRECISION Amax , Scond
      CHARACTER Uplo
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , S(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
      INTEGER MAX_ITER
      PARAMETER (MAX_ITER=100)
!     ..
!     .. Local Scalars ..
      INTEGER i , j , iter
      DOUBLE PRECISION avg , std , tol , c0 , c1 , c2 , t , u , si , d ,&
     &                 base , smin , smax , smlnum , bignum , scale ,   &
     &                 sumsq
      LOGICAL up
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      LOGICAL LSAME
      EXTERNAL DLAMCH , LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL DLASSQ , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , INT , LOG , MAX , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      IF ( .NOT.(LSAME(Uplo,'U') .OR. LSAME(Uplo,'L')) ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DSYEQUB',-Info)
         RETURN
      ENDIF
 
      up = LSAME(Uplo,'U')
      Amax = ZERO
!
!     Quick return if possible.
!
      IF ( N==0 ) THEN
         Scond = ONE
         RETURN
      ENDIF
 
      DO i = 1 , N
         S(i) = ZERO
      ENDDO
 
      Amax = ZERO
      IF ( up ) THEN
         DO j = 1 , N
            DO i = 1 , j - 1
               S(i) = MAX(S(i),ABS(A(i,j)))
               S(j) = MAX(S(j),ABS(A(i,j)))
               Amax = MAX(Amax,ABS(A(i,j)))
            ENDDO
            S(j) = MAX(S(j),ABS(A(j,j)))
            Amax = MAX(Amax,ABS(A(j,j)))
         ENDDO
      ELSE
         DO j = 1 , N
            S(j) = MAX(S(j),ABS(A(j,j)))
            Amax = MAX(Amax,ABS(A(j,j)))
            DO i = j + 1 , N
               S(i) = MAX(S(i),ABS(A(i,j)))
               S(j) = MAX(S(j),ABS(A(i,j)))
               Amax = MAX(Amax,ABS(A(i,j)))
            ENDDO
         ENDDO
      ENDIF
      DO j = 1 , N
         S(j) = 1.0D0/S(j)
      ENDDO
 
      tol = ONE/SQRT(2.0D0*N)
 
      DO iter = 1 , MAX_ITER
         scale = 0.0D0
         sumsq = 0.0D0
!        beta = |A|s
         DO i = 1 , N
            Work(i) = ZERO
         ENDDO
         IF ( up ) THEN
            DO j = 1 , N
               DO i = 1 , j - 1
                  Work(i) = Work(i) + ABS(A(i,j))*S(j)
                  Work(j) = Work(j) + ABS(A(i,j))*S(i)
               ENDDO
               Work(j) = Work(j) + ABS(A(j,j))*S(j)
            ENDDO
         ELSE
            DO j = 1 , N
               Work(j) = Work(j) + ABS(A(j,j))*S(j)
               DO i = j + 1 , N
                  Work(i) = Work(i) + ABS(A(i,j))*S(j)
                  Work(j) = Work(j) + ABS(A(i,j))*S(i)
               ENDDO
            ENDDO
         ENDIF
 
!        avg = s^T beta / n
         avg = 0.0D0
         DO i = 1 , N
            avg = avg + S(i)*Work(i)
         ENDDO
         avg = avg/N
 
         std = 0.0D0
         DO i = N + 1 , 2*N
            Work(i) = S(i-N)*Work(i-N) - avg
         ENDDO
         CALL DLASSQ(N,Work(N+1),1,scale,sumsq)
         std = scale*SQRT(sumsq/N)
 
         IF ( std<tol*avg ) EXIT
 
         DO i = 1 , N
            t = ABS(A(i,i))
            si = S(i)
            c2 = (N-1)*t
            c1 = (N-2)*(Work(i)-t*si)
            c0 = -(t*si)*si + 2*Work(i)*si - N*avg
            d = c1*c1 - 4*c0*c2
 
            IF ( d<=0 ) THEN
               Info = -1
               RETURN
            ENDIF
            si = -2*c0/(c1+SQRT(d))
 
            d = si - S(i)
            u = ZERO
            IF ( up ) THEN
               DO j = 1 , i
                  t = ABS(A(j,i))
                  u = u + S(j)*t
                  Work(j) = Work(j) + d*t
               ENDDO
               DO j = i + 1 , N
                  t = ABS(A(i,j))
                  u = u + S(j)*t
                  Work(j) = Work(j) + d*t
               ENDDO
            ELSE
               DO j = 1 , i
                  t = ABS(A(i,j))
                  u = u + S(j)*t
                  Work(j) = Work(j) + d*t
               ENDDO
               DO j = i + 1 , N
                  t = ABS(A(j,i))
                  u = u + S(j)*t
                  Work(j) = Work(j) + d*t
               ENDDO
            ENDIF
 
            avg = avg + (u+Work(i))*d/N
            S(i) = si
         ENDDO
      ENDDO
 
 
      smlnum = DLAMCH('SAFEMIN')
      bignum = ONE/smlnum
      smin = bignum
      smax = ZERO
      t = ONE/SQRT(avg)
      base = DLAMCH('B')
      u = ONE/LOG(base)
      DO i = 1 , N
         S(i) = base**INT(u*LOG(S(i)*t))
         smin = MIN(smin,S(i))
         smax = MAX(smax,S(i))
      ENDDO
      Scond = MAX(smin,smlnum)/MIN(smax,bignum)
!
      END SUBROUTINE DSYEQUB
