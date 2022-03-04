!*==zlaswp.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLASWP performs a series of row interchanges on a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLASWP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaswp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaswp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaswp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, K1, K2, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLASWP performs a series of row interchanges on the matrix A.
!> One row interchange is initiated for each of rows K1 through K2 of A.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the matrix of column dimension N to which the row
!>          interchanges will be applied.
!>          On exit, the permuted matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[in] K1
!> \verbatim
!>          K1 is INTEGER
!>          The first element of IPIV for which a row interchange will
!>          be done.
!> \endverbatim
!>
!> \param[in] K2
!> \verbatim
!>          K2 is INTEGER
!>          (K2-K1+1) is the number of elements of IPIV for which a row
!>          interchange will be done.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (K1+(K2-K1)*abs(INCX))
!>          The vector of pivot indices. Only the elements in positions
!>          K1 through K1+(K2-K1)*abs(INCX) of IPIV are accessed.
!>          IPIV(K1+(K-K1)*abs(INCX)) = L implies rows K and L are to be
!>          interchanged.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of IPIV. If INCX
!>          is negative, the pivots are applied in reverse order.
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
!> \date June 2017
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Modified by
!>   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZLASWP(N,A,Lda,K1,K2,Ipiv,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!*--ZLASWP120
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      INTEGER , INTENT(IN) :: K1
      INTEGER , INTENT(IN) :: K2
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , i1 , i2 , inc , ip , ix , ix0 , j , k , n32
      COMPLEX(CX16KIND) :: temp
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
!     K1 through K2.
!
      IF ( Incx>0 ) THEN
         ix0 = K1
         i1 = K1
         i2 = K2
         inc = 1
      ELSEIF ( Incx<0 ) THEN
         ix0 = K1 + (K1-K2)*Incx
         i1 = K2
         i2 = K1
         inc = -1
      ELSE
         RETURN
      ENDIF
!
      n32 = (N/32)*32
      IF ( n32/=0 ) THEN
         DO j = 1 , n32 , 32
            ix = ix0
            DO i = i1 , i2 , inc
               ip = Ipiv(ix)
               IF ( ip/=i ) THEN
                  DO k = j , j + 31
                     temp = A(i,k)
                     A(i,k) = A(ip,k)
                     A(ip,k) = temp
                  ENDDO
               ENDIF
               ix = ix + Incx
            ENDDO
         ENDDO
      ENDIF
      IF ( n32/=N ) THEN
         n32 = n32 + 1
         ix = ix0
         DO i = i1 , i2 , inc
            ip = Ipiv(ix)
            IF ( ip/=i ) THEN
               DO k = n32 , N
                  temp = A(i,k)
                  A(i,k) = A(ip,k)
                  A(ip,k) = temp
               ENDDO
            ENDIF
            ix = ix + Incx
         ENDDO
      ENDIF
!
!
!     End of ZLASWP
!
      END SUBROUTINE ZLASWP
