!*==zlaqgb.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLAQGB scales a general band matrix, using row and column scaling factors computed by sgbequ.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAQGB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqgb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqgb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqgb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAQGB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,
!                          AMAX, EQUED )
!
!       .. Scalar Arguments ..
!       CHARACTER          EQUED
!       INTEGER            KL, KU, LDAB, M, N
!       DOUBLE PRECISION   AMAX, COLCND, ROWCND
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( * ), R( * )
!       COMPLEX*16         AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAQGB equilibrates a general M by N band matrix A with KL
!> subdiagonals and KU superdiagonals using the row and scaling factors
!> in the vectors R and C.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!>
!>          On exit, the equilibrated matrix, in the same storage format
!>          as A.  See EQUED for the form of the equilibrated matrix.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDA >= KL+KU+1.
!> \endverbatim
!>
!> \param[in] R
!> \verbatim
!>          R is DOUBLE PRECISION array, dimension (M)
!>          The row scale factors for A.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (N)
!>          The column scale factors for A.
!> \endverbatim
!>
!> \param[in] ROWCND
!> \verbatim
!>          ROWCND is DOUBLE PRECISION
!>          Ratio of the smallest R(i) to the largest R(i).
!> \endverbatim
!>
!> \param[in] COLCND
!> \verbatim
!>          COLCND is DOUBLE PRECISION
!>          Ratio of the smallest C(i) to the largest C(i).
!> \endverbatim
!>
!> \param[in] AMAX
!> \verbatim
!>          AMAX is DOUBLE PRECISION
!>          Absolute value of largest matrix entry.
!> \endverbatim
!>
!> \param[out] EQUED
!> \verbatim
!>          EQUED is CHARACTER*1
!>          Specifies the form of equilibration that was done.
!>          = 'N':  No equilibration
!>          = 'R':  Row equilibration, i.e., A has been premultiplied by
!>                  diag(R).
!>          = 'C':  Column equilibration, i.e., A has been postmultiplied
!>                  by diag(C).
!>          = 'B':  Both row and column equilibration, i.e., A has been
!>                  replaced by diag(R) * A * diag(C).
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  THRESH is a threshold value used to decide if row or column scaling
!>  should be done based on the ratio of the row or column scaling
!>  factors.  If ROWCND < THRESH, row scaling is done, and if
!>  COLCND < THRESH, column scaling is done.
!>
!>  LARGE and SMALL are threshold values used to decide if row scaling
!>  should be done based on the absolute size of the largest matrix
!>  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.
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
!> \ingroup complex16GBauxiliary
!
!  =====================================================================
      SUBROUTINE ZLAQGB(M,N,Kl,Ku,Ab,Ldab,R,C,Rowcnd,Colcnd,Amax,Equed)
      USE F77KINDS                        
      USE S_DLAMCH
      IMPLICIT NONE
!*--ZLAQGB165
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , THRESH = 0.1D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(IN) :: Rowcnd
      REAL(R8KIND) , INTENT(IN) :: Colcnd
      REAL(R8KIND) , INTENT(IN) :: Amax
      CHARACTER , INTENT(OUT) :: Equed
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: cj , large , small
      INTEGER :: i , j
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
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( M<=0 .OR. N<=0 ) THEN
         Equed = 'N'
         RETURN
      ENDIF
!
!     Initialize LARGE and SMALL.
!
      small = DLAMCH('Safe minimum')/DLAMCH('Precision')
      large = ONE/small
!
      IF ( Rowcnd>=THRESH .AND. Amax>=small .AND. Amax<=large ) THEN
!
!        No row scaling
!
         IF ( Colcnd>=THRESH ) THEN
!
!           No column scaling
!
            Equed = 'N'
         ELSE
!
!           Column scaling
!
            DO j = 1 , N
               cj = C(j)
               DO i = MAX(1,j-Ku) , MIN(M,j+Kl)
                  Ab(Ku+1+i-j,j) = cj*Ab(Ku+1+i-j,j)
               ENDDO
            ENDDO
            Equed = 'C'
         ENDIF
      ELSEIF ( Colcnd>=THRESH ) THEN
!
!        Row scaling, no column scaling
!
         DO j = 1 , N
            DO i = MAX(1,j-Ku) , MIN(M,j+Kl)
               Ab(Ku+1+i-j,j) = R(i)*Ab(Ku+1+i-j,j)
            ENDDO
         ENDDO
         Equed = 'R'
      ELSE
!
!        Row and column scaling
!
         DO j = 1 , N
            cj = C(j)
            DO i = MAX(1,j-Ku) , MIN(M,j+Kl)
               Ab(Ku+1+i-j,j) = cj*R(i)*Ab(Ku+1+i-j,j)
            ENDDO
         ENDDO
         Equed = 'B'
      ENDIF
!
!
!     End of ZLAQGB
!
      END SUBROUTINE ZLAQGB
