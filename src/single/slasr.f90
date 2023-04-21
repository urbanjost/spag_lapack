!*==slasr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLASR applies a sequence of plane rotations to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, PIVOT, SIDE
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), C( * ), S( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLASR applies a sequence of plane rotations to a real matrix A,
!> from either the left or the right.
!>
!> When SIDE = 'L', the transformation takes the form
!>
!>    A := P*A
!>
!> and when SIDE = 'R', the transformation takes the form
!>
!>    A := A*P**T
!>
!> where P is an orthogonal matrix consisting of a sequence of z plane
!> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
!> and P**T is the transpose of P.
!>
!> When DIRECT = 'F' (Forward sequence), then
!>
!>    P = P(z-1) * ... * P(2) * P(1)
!>
!> and when DIRECT = 'B' (Backward sequence), then
!>
!>    P = P(1) * P(2) * ... * P(z-1)
!>
!> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
!>
!>    R(k) = (  c(k)  s(k) )
!>         = ( -s(k)  c(k) ).
!>
!> When PIVOT = 'V' (Variable pivot), the rotation is performed
!> for the plane (k,k+1), i.e., P(k) has the form
!>
!>    P(k) = (  1                                            )
!>           (       ...                                     )
!>           (              1                                )
!>           (                   c(k)  s(k)                  )
!>           (                  -s(k)  c(k)                  )
!>           (                                1              )
!>           (                                     ...       )
!>           (                                            1  )
!>
!> where R(k) appears as a rank-2 modification to the identity matrix in
!> rows and columns k and k+1.
!>
!> When PIVOT = 'T' (Top pivot), the rotation is performed for the
!> plane (1,k+1), so P(k) has the form
!>
!>    P(k) = (  c(k)                    s(k)                 )
!>           (         1                                     )
!>           (              ...                              )
!>           (                     1                         )
!>           ( -s(k)                    c(k)                 )
!>           (                                 1             )
!>           (                                      ...      )
!>           (                                             1 )
!>
!> where R(k) appears in rows and columns 1 and k+1.
!>
!> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
!> performed for the plane (k,z), giving P(k) the form
!>
!>    P(k) = ( 1                                             )
!>           (      ...                                      )
!>           (             1                                 )
!>           (                  c(k)                    s(k) )
!>           (                         1                     )
!>           (                              ...              )
!>           (                                     1         )
!>           (                 -s(k)                    c(k) )
!>
!> where R(k) appears in rows and columns k and z.  The rotations are
!> performed without ever forming P(k) explicitly.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          Specifies whether the plane rotation matrix P is applied to
!>          A on the left or the right.
!>          = 'L':  Left, compute A := P*A
!>          = 'R':  Right, compute A:= A*P**T
!> \endverbatim
!>
!> \param[in] PIVOT
!> \verbatim
!>          PIVOT is CHARACTER*1
!>          Specifies the plane for which P(k) is a plane rotation
!>          matrix.
!>          = 'V':  Variable pivot, the plane (k,k+1)
!>          = 'T':  Top pivot, the plane (1,k+1)
!>          = 'B':  Bottom pivot, the plane (k,z)
!> \endverbatim
!>
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Specifies whether P is a forward or backward sequence of
!>          plane rotations.
!>          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)
!>          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  If m <= 1, an immediate
!>          return is effected.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  If n <= 1, an
!>          immediate return is effected.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension
!>                  (M-1) if SIDE = 'L'
!>                  (N-1) if SIDE = 'R'
!>          The cosines c(k) of the plane rotations.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension
!>                  (M-1) if SIDE = 'L'
!>                  (N-1) if SIDE = 'R'
!>          The sines s(k) of the plane rotations.  The 2-by-2 plane
!>          rotation part of the matrix P(k), R(k), has the form
!>          R(k) = (  c(k)  s(k) )
!>                 ( -s(k)  c(k) ).
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The M-by-N matrix A.  On exit, A is overwritten by P*A if
!>          SIDE = 'R' or by A*P**T if SIDE = 'L'.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLASR(Side,Pivot,Direct,M,N,C,S,A,Lda)
      IMPLICIT NONE
!*--SLASR203
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Direct , Pivot , Side
      INTEGER Lda , M , N
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , C(*) , S(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j
      REAL ctemp , stemp , temp
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      info = 0
      IF ( .NOT.(LSAME(Side,'L') .OR. LSAME(Side,'R')) ) THEN
         info = 1
      ELSEIF ( .NOT.(LSAME(Pivot,'V') .OR. LSAME(Pivot,'T') .OR.        &
     &         LSAME(Pivot,'B')) ) THEN
         info = 2
      ELSEIF ( .NOT.(LSAME(Direct,'F') .OR. LSAME(Direct,'B')) ) THEN
         info = 3
      ELSEIF ( M<0 ) THEN
         info = 4
      ELSEIF ( N<0 ) THEN
         info = 5
      ELSEIF ( Lda<MAX(1,M) ) THEN
         info = 9
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('SLASR ',info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( (M==0) .OR. (N==0) ) RETURN
      IF ( LSAME(Side,'L') ) THEN
!
!        Form  P * A
!
         IF ( LSAME(Pivot,'V') ) THEN
            IF ( LSAME(Direct,'F') ) THEN
               DO j = 1 , M - 1
                  ctemp = C(j)
                  stemp = S(j)
                  IF ( (ctemp/=ONE) .OR. (stemp/=ZERO) ) THEN
                     DO i = 1 , N
                        temp = A(j+1,i)
                        A(j+1,i) = ctemp*temp - stemp*A(j,i)
                        A(j,i) = stemp*temp + ctemp*A(j,i)
                     ENDDO
                  ENDIF
               ENDDO
            ELSEIF ( LSAME(Direct,'B') ) THEN
               DO j = M - 1 , 1 , -1
                  ctemp = C(j)
                  stemp = S(j)
                  IF ( (ctemp/=ONE) .OR. (stemp/=ZERO) ) THEN
                     DO i = 1 , N
                        temp = A(j+1,i)
                        A(j+1,i) = ctemp*temp - stemp*A(j,i)
                        A(j,i) = stemp*temp + ctemp*A(j,i)
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ELSEIF ( LSAME(Pivot,'T') ) THEN
            IF ( LSAME(Direct,'F') ) THEN
               DO j = 2 , M
                  ctemp = C(j-1)
                  stemp = S(j-1)
                  IF ( (ctemp/=ONE) .OR. (stemp/=ZERO) ) THEN
                     DO i = 1 , N
                        temp = A(j,i)
                        A(j,i) = ctemp*temp - stemp*A(1,i)
                        A(1,i) = stemp*temp + ctemp*A(1,i)
                     ENDDO
                  ENDIF
               ENDDO
            ELSEIF ( LSAME(Direct,'B') ) THEN
               DO j = M , 2 , -1
                  ctemp = C(j-1)
                  stemp = S(j-1)
                  IF ( (ctemp/=ONE) .OR. (stemp/=ZERO) ) THEN
                     DO i = 1 , N
                        temp = A(j,i)
                        A(j,i) = ctemp*temp - stemp*A(1,i)
                        A(1,i) = stemp*temp + ctemp*A(1,i)
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ELSEIF ( LSAME(Pivot,'B') ) THEN
            IF ( LSAME(Direct,'F') ) THEN
               DO j = 1 , M - 1
                  ctemp = C(j)
                  stemp = S(j)
                  IF ( (ctemp/=ONE) .OR. (stemp/=ZERO) ) THEN
                     DO i = 1 , N
                        temp = A(j,i)
                        A(j,i) = stemp*A(M,i) + ctemp*temp
                        A(M,i) = ctemp*A(M,i) - stemp*temp
                     ENDDO
                  ENDIF
               ENDDO
            ELSEIF ( LSAME(Direct,'B') ) THEN
               DO j = M - 1 , 1 , -1
                  ctemp = C(j)
                  stemp = S(j)
                  IF ( (ctemp/=ONE) .OR. (stemp/=ZERO) ) THEN
                     DO i = 1 , N
                        temp = A(j,i)
                        A(j,i) = stemp*A(M,i) + ctemp*temp
                        A(M,i) = ctemp*A(M,i) - stemp*temp
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ELSEIF ( LSAME(Side,'R') ) THEN
!
!        Form A * P**T
!
         IF ( LSAME(Pivot,'V') ) THEN
            IF ( LSAME(Direct,'F') ) THEN
               DO j = 1 , N - 1
                  ctemp = C(j)
                  stemp = S(j)
                  IF ( (ctemp/=ONE) .OR. (stemp/=ZERO) ) THEN
                     DO i = 1 , M
                        temp = A(i,j+1)
                        A(i,j+1) = ctemp*temp - stemp*A(i,j)
                        A(i,j) = stemp*temp + ctemp*A(i,j)
                     ENDDO
                  ENDIF
               ENDDO
            ELSEIF ( LSAME(Direct,'B') ) THEN
               DO j = N - 1 , 1 , -1
                  ctemp = C(j)
                  stemp = S(j)
                  IF ( (ctemp/=ONE) .OR. (stemp/=ZERO) ) THEN
                     DO i = 1 , M
                        temp = A(i,j+1)
                        A(i,j+1) = ctemp*temp - stemp*A(i,j)
                        A(i,j) = stemp*temp + ctemp*A(i,j)
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ELSEIF ( LSAME(Pivot,'T') ) THEN
            IF ( LSAME(Direct,'F') ) THEN
               DO j = 2 , N
                  ctemp = C(j-1)
                  stemp = S(j-1)
                  IF ( (ctemp/=ONE) .OR. (stemp/=ZERO) ) THEN
                     DO i = 1 , M
                        temp = A(i,j)
                        A(i,j) = ctemp*temp - stemp*A(i,1)
                        A(i,1) = stemp*temp + ctemp*A(i,1)
                     ENDDO
                  ENDIF
               ENDDO
            ELSEIF ( LSAME(Direct,'B') ) THEN
               DO j = N , 2 , -1
                  ctemp = C(j-1)
                  stemp = S(j-1)
                  IF ( (ctemp/=ONE) .OR. (stemp/=ZERO) ) THEN
                     DO i = 1 , M
                        temp = A(i,j)
                        A(i,j) = ctemp*temp - stemp*A(i,1)
                        A(i,1) = stemp*temp + ctemp*A(i,1)
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ELSEIF ( LSAME(Pivot,'B') ) THEN
            IF ( LSAME(Direct,'F') ) THEN
               DO j = 1 , N - 1
                  ctemp = C(j)
                  stemp = S(j)
                  IF ( (ctemp/=ONE) .OR. (stemp/=ZERO) ) THEN
                     DO i = 1 , M
                        temp = A(i,j)
                        A(i,j) = stemp*A(i,N) + ctemp*temp
                        A(i,N) = ctemp*A(i,N) - stemp*temp
                     ENDDO
                  ENDIF
               ENDDO
            ELSEIF ( LSAME(Direct,'B') ) THEN
               DO j = N - 1 , 1 , -1
                  ctemp = C(j)
                  stemp = S(j)
                  IF ( (ctemp/=ONE) .OR. (stemp/=ZERO) ) THEN
                     DO i = 1 , M
                        temp = A(i,j)
                        A(i,j) = stemp*A(i,N) + ctemp*temp
                        A(i,N) = ctemp*A(i,N) - stemp*temp
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
!
!     End of SLASR
!
      END SUBROUTINE SLASR
