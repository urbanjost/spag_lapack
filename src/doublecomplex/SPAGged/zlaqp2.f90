!*==zlaqp2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLAQP2 computes a QR factorization with column pivoting of the matrix block.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAQP2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqp2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqp2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqp2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2,
!                          WORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, M, N, OFFSET
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       DOUBLE PRECISION   VN1( * ), VN2( * )
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAQP2 computes a QR factorization with column pivoting of
!> the block A(OFFSET+1:M,1:N).
!> The block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A. N >= 0.
!> \endverbatim
!>
!> \param[in] OFFSET
!> \verbatim
!>          OFFSET is INTEGER
!>          The number of rows of the matrix A that must be pivoted
!>          but no factorized. OFFSET >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the upper triangle of block A(OFFSET+1:M,1:N) is
!>          the triangular factor obtained; the elements in block
!>          A(OFFSET+1:M,1:N) below the diagonal, together with the
!>          array TAU, represent the orthogonal matrix Q as a product of
!>          elementary reflectors. Block A(1:OFFSET,1:N) has been
!>          accordingly pivoted, but no factorized.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] JPVT
!> \verbatim
!>          JPVT is INTEGER array, dimension (N)
!>          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
!>          to the front of A*P (a leading column); if JPVT(i) = 0,
!>          the i-th column of A is a free column.
!>          On exit, if JPVT(i) = k, then the i-th column of A*P
!>          was the k-th column of A.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[in,out] VN1
!> \verbatim
!>          VN1 is DOUBLE PRECISION array, dimension (N)
!>          The vector with the partial column norms.
!> \endverbatim
!>
!> \param[in,out] VN2
!> \verbatim
!>          VN2 is DOUBLE PRECISION array, dimension (N)
!>          The vector with the exact column norms.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N)
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
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
!>    X. Sun, Computer Science Dept., Duke University, USA
!> \n
!>  Partial column norm updating strategy modified on April 2011
!>    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,
!>    University of Zagreb, Croatia.
!
!> \par References:
!  ================
!>
!> LAPACK Working Note 176
!
!> \htmlonly
!> <a href="http://www.netlib.org/lapack/lawnspdf/lawn176.pdf">[PDF]</a>
!> \endhtmlonly
!
!  =====================================================================
      SUBROUTINE ZLAQP2(M,N,Offset,A,Lda,Jpvt,Tau,Vn1,Vn2,Work)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_DZNRM2
      USE S_IDAMAX
      USE S_ZLARF
      USE S_ZLARFG
      USE S_ZSWAP
      IMPLICIT NONE
!*--ZLAQP2159
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Offset
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jpvt
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vn1
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Vn2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) :: aii
      INTEGER :: i , itemp , j , mn , offpi , pvt
      REAL(R8KIND) :: temp , temp2 , tol3z
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
      mn = MIN(M-Offset,N)
      tol3z = SQRT(DLAMCH('Epsilon'))
!
!     Compute factorization.
!
      DO i = 1 , mn
!
         offpi = Offset + i
!
!        Determine ith pivot column and swap if necessary.
!
         pvt = (i-1) + IDAMAX(N-i+1,Vn1(i),1)
!
         IF ( pvt/=i ) THEN
            CALL ZSWAP(M,A(1,pvt),1,A(1,i),1)
            itemp = Jpvt(pvt)
            Jpvt(pvt) = Jpvt(i)
            Jpvt(i) = itemp
            Vn1(pvt) = Vn1(i)
            Vn2(pvt) = Vn2(i)
         ENDIF
!
!        Generate elementary reflector H(i).
!
         IF ( offpi<M ) THEN
            CALL ZLARFG(M-offpi+1,A(offpi,i),A(offpi+1,i),1,Tau(i))
         ELSE
            CALL ZLARFG(1,A(M,i),A(M,i),1,Tau(i))
         ENDIF
!
         IF ( i<N ) THEN
!
!           Apply H(i)**H to A(offset+i:m,i+1:n) from the left.
!
            aii = A(offpi,i)
            A(offpi,i) = CONE
            CALL ZLARF('Left',M-offpi+1,N-i,A(offpi,i),1,DCONJG(Tau(i)),&
     &                 A(offpi,i+1),Lda,Work(1))
            A(offpi,i) = aii
         ENDIF
!
!        Update partial column norms.
!
         DO j = i + 1 , N
            IF ( Vn1(j)/=ZERO ) THEN
!
!              NOTE: The following 4 lines follow from the analysis in
!              Lapack Working Note 176.
!
               temp = ONE - (ABS(A(offpi,j))/Vn1(j))**2
               temp = MAX(temp,ZERO)
               temp2 = temp*(Vn1(j)/Vn2(j))**2
               IF ( temp2>tol3z ) THEN
                  Vn1(j) = Vn1(j)*SQRT(temp)
               ELSEIF ( offpi<M ) THEN
                  Vn1(j) = DZNRM2(M-offpi,A(offpi+1,j),1)
                  Vn2(j) = Vn1(j)
               ELSE
                  Vn1(j) = ZERO
                  Vn2(j) = ZERO
               ENDIF
            ENDIF
         ENDDO
!
      ENDDO
!
!
!     End of ZLAQP2
!
      END SUBROUTINE ZLAQP2
