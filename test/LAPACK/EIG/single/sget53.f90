!*==sget53.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b SGET53
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGET53( A, LDA, B, LDB, SCALE, WR, WI, RESULT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB
!       REAL               RESULT, SCALE, WI, WR
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGET53  checks the generalized eigenvalues computed by SLAG2.
!>
!> The basic test for an eigenvalue is:
!>
!>                              | det( s A - w B ) |
!>     RESULT =  ---------------------------------------------------
!>               ulp max( s norm(A), |w| norm(B) )*norm( s A - w B )
!>
!> Two "safety checks" are performed:
!>
!> (1)  ulp*max( s*norm(A), |w|*norm(B) )  must be at least
!>      safe_minimum.  This insures that the test performed is
!>      not essentially  det(0*A + 0*B)=0.
!>
!> (2)  s*norm(A) + |w|*norm(B) must be less than 1/safe_minimum.
!>      This insures that  s*A - w*B  will not overflow.
!>
!> If these tests are not passed, then  s  and  w  are scaled and
!> tested anyway, if this is possible.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA, 2)
!>          The 2x2 matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 2.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB, N)
!>          The 2x2 upper-triangular matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of B.  It must be at least 2.
!> \endverbatim
!>
!> \param[in] SCALE
!> \verbatim
!>          SCALE is REAL
!>          The "scale factor" s in the formula  s A - w B .  It is
!>          assumed to be non-negative.
!> \endverbatim
!>
!> \param[in] WR
!> \verbatim
!>          WR is REAL
!>          The real part of the eigenvalue  w  in the formula
!>          s A - w B .
!> \endverbatim
!>
!> \param[in] WI
!> \verbatim
!>          WI is REAL
!>          The imaginary part of the eigenvalue  w  in the formula
!>          s A - w B .
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL
!>          If INFO is 2 or less, the value computed by the test
!>             described above.
!>          If INFO=3, this will just be 1/ulp.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          =0:  The input data pass the "safety checks".
!>          =1:  s*norm(A) + |w|*norm(B) > 1/safe_minimum.
!>          =2:  ulp*max( s*norm(A), |w|*norm(B) ) < safe_minimum
!>          =3:  same as INFO=2, but  s  and  w  could not be scaled so
!>               as to compute the test.
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
      SUBROUTINE SGET53(A,Lda,B,Ldb,Scale,Wr,Wi,Result,Info)
      IMPLICIT NONE
!*--SGET53130
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldb
      REAL Result , Scale , Wi , Wr
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , B(Ldb,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0,ONE=1.0)
!     ..
!     .. Local Scalars ..
      REAL absw , anorm , bnorm , ci11 , ci12 , ci22 , cnorm , cr11 ,   &
     &     cr12 , cr21 , cr22 , cscale , deti , detr , s1 , safmin ,    &
     &     scales , sigmin , temp , ulp , wis , wrs
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SQRT
!     ..
!     .. Executable Statements ..
!
!     Initialize
!
      Info = 0
      Result = ZERO
      scales = Scale
      wrs = Wr
      wis = Wi
!
!     Machine constants and norms
!
      safmin = SLAMCH('Safe minimum')
      ulp = SLAMCH('Epsilon')*SLAMCH('Base')
      absw = ABS(wrs) + ABS(wis)
      anorm = MAX(ABS(A(1,1))+ABS(A(2,1)),ABS(A(1,2))+ABS(A(2,2)),      &
     &        safmin)
      bnorm = MAX(ABS(B(1,1)),ABS(B(1,2))+ABS(B(2,2)),safmin)
!
!     Check for possible overflow.
!
      temp = (safmin*bnorm)*absw + (safmin*anorm)*scales
      IF ( temp>=ONE ) THEN
!
!        Scale down to avoid overflow
!
         Info = 1
         temp = ONE/temp
         scales = scales*temp
         wrs = wrs*temp
         wis = wis*temp
         absw = ABS(wrs) + ABS(wis)
      ENDIF
      s1 = MAX(ulp*MAX(scales*anorm,absw*bnorm),safmin*MAX(scales,absw))
!
!     Check for W and SCALE essentially zero.
!
      IF ( s1<safmin ) THEN
         Info = 2
         IF ( scales<safmin .AND. absw<safmin ) THEN
            Info = 3
            Result = ONE/ulp
            RETURN
         ENDIF
!
!        Scale up to avoid underflow
!
         temp = ONE/MAX(scales*anorm+absw*bnorm,safmin)
         scales = scales*temp
         wrs = wrs*temp
         wis = wis*temp
         absw = ABS(wrs) + ABS(wis)
         s1 = MAX(ulp*MAX(scales*anorm,absw*bnorm),                     &
     &        safmin*MAX(scales,absw))
         IF ( s1<safmin ) THEN
            Info = 3
            Result = ONE/ulp
            RETURN
         ENDIF
      ENDIF
!
!     Compute C = s A - w B
!
      cr11 = scales*A(1,1) - wrs*B(1,1)
      ci11 = -wis*B(1,1)
      cr21 = scales*A(2,1)
      cr12 = scales*A(1,2) - wrs*B(1,2)
      ci12 = -wis*B(1,2)
      cr22 = scales*A(2,2) - wrs*B(2,2)
      ci22 = -wis*B(2,2)
!
!     Compute the smallest singular value of s A - w B:
!
!                 |det( s A - w B )|
!     sigma_min = ------------------
!                 norm( s A - w B )
!
      cnorm = MAX(ABS(cr11)+ABS(ci11)+ABS(cr21),ABS(cr12)+ABS(ci12)     &
     &        +ABS(cr22)+ABS(ci22),safmin)
      cscale = ONE/SQRT(cnorm)
      detr = (cscale*cr11)*(cscale*cr22) - (cscale*ci11)*(cscale*ci22)  &
     &       - (cscale*cr12)*(cscale*cr21)
      deti = (cscale*cr11)*(cscale*ci22) + (cscale*ci11)*(cscale*cr22)  &
     &       - (cscale*ci12)*(cscale*cr21)
      sigmin = ABS(detr) + ABS(deti)
      Result = sigmin/s1
!
!     End of SGET53
!
      END SUBROUTINE SGET53
