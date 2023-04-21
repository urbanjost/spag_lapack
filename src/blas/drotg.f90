!*==drotg.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DROTG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DROTG(DA,DB,C,S)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION C,DA,DB,S
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DROTG construct givens plane rotation.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] DA
!> \verbatim
!>          DA is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] DB
!> \verbatim
!>          DB is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION
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
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DROTG(Da,Db,C,S)
      IMPLICIT NONE
!*--DROTG73
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION C , Da , Db , S
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION r , roe , scale , z
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DABS , DSIGN , DSQRT
!     ..
      scale = DABS(Da) + DABS(Db)
      IF ( scale==0.0D0 ) THEN
         C = 1.0D0
         S = 0.0D0
         r = 0.0D0
         z = 0.0D0
      ELSE
         roe = Db
         IF ( DABS(Da)>DABS(Db) ) roe = Da
         r = scale*DSQRT((Da/scale)**2+(Db/scale)**2)
         r = DSIGN(1.0D0,roe)*r
         C = Da/r
         S = Db/r
         z = 1.0D0
         IF ( DABS(Da)>DABS(Db) ) z = S
         IF ( DABS(Db)>=DABS(Da) .AND. C/=0.0D0 ) z = 1.0D0/C
      ENDIF
      Da = r
      Db = z
      END SUBROUTINE DROTG
