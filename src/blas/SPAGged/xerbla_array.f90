!*==xerbla_array.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b XERBLA_ARRAY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE XERBLA_ARRAY(SRNAME_ARRAY, SRNAME_LEN, INFO)
!
!       .. Scalar Arguments ..
!       INTEGER SRNAME_LEN, INFO
!       ..
!       .. Array Arguments ..
!       CHARACTER(1) SRNAME_ARRAY(SRNAME_LEN)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> XERBLA_ARRAY assists other languages in calling XERBLA, the LAPACK
!> and BLAS error handler.  Rather than taking a Fortran string argument
!> as the function's name, XERBLA_ARRAY takes an array of single
!> characters along with the array's length.  XERBLA_ARRAY then copies
!> up to 32 characters of that array into a Fortran string and passes
!> that to XERBLA.  If called with a non-positive SRNAME_LEN,
!> XERBLA_ARRAY will call XERBLA with a string of all blank characters.
!>
!> Say some macro or other device makes XERBLA_ARRAY available to C99
!> by a name lapack_xerbla and with a common Fortran calling convention.
!> Then a C99 program could invoke XERBLA via:
!>    {
!>      int flen = strlen(__func__);
!>      lapack_xerbla(__func__, &flen, &info);
!>    }
!>
!> Providing XERBLA_ARRAY is not necessary for intercepting LAPACK
!> errors.  XERBLA_ARRAY calls XERBLA.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SRNAME_ARRAY
!> \verbatim
!>          SRNAME_ARRAY is CHARACTER(1) array, dimension (SRNAME_LEN)
!>          The name of the routine which called XERBLA_ARRAY.
!> \endverbatim
!>
!> \param[in] SRNAME_LEN
!> \verbatim
!>          SRNAME_LEN is INTEGER
!>          The length of the name in SRNAME_ARRAY.
!> \endverbatim
!>
!> \param[in] INFO
!> \verbatim
!>          INFO is INTEGER
!>          The position of the invalid parameter in the parameter list
!>          of the calling routine.
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
!> \ingroup aux_blas
!
!  =====================================================================
      SUBROUTINE XERBLA_ARRAY(Srname_array,Srname_len,Info)
      USE S_XERBLA
      IMPLICIT NONE
!*--XERBLA_ARRAY85
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER(1) , INTENT(IN) , DIMENSION(Srname_len) :: Srname_array
      INTEGER , INTENT(IN) :: Srname_len
      INTEGER :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i
      CHARACTER(32) :: srname
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
!
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
      srname = ''
      DO i = 1 , MIN(Srname_len,LEN(srname))
         srname(i:i) = Srname_array(i)
      ENDDO
 
      CALL XERBLA(srname,Info)
 
      END SUBROUTINE XERBLA_ARRAY
