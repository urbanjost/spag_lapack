!*==zpttrs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZPTTRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZPTTRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpttrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpttrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpttrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * )
!       COMPLEX*16         B( LDB, * ), E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPTTRS solves a tridiagonal system of the form
!>    A * X = B
!> using the factorization A = U**H *D* U or A = L*D*L**H computed by ZPTTRF.
!> D is a diagonal matrix specified in the vector D, U (or L) is a unit
!> bidiagonal matrix whose superdiagonal (subdiagonal) is specified in
!> the vector E, and X and B are N by NRHS matrices.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the form of the factorization and whether the
!>          vector E is the superdiagonal of the upper bidiagonal factor
!>          U or the subdiagonal of the lower bidiagonal factor L.
!>          = 'U':  A = U**H *D*U, E is the superdiagonal of U
!>          = 'L':  A = L*D*L**H, E is the subdiagonal of L
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the tridiagonal matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The n diagonal elements of the diagonal matrix D from the
!>          factorization A = U**H *D*U or A = L*D*L**H.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX*16 array, dimension (N-1)
!>          If UPLO = 'U', the (n-1) superdiagonal elements of the unit
!>          bidiagonal factor U from the factorization A = U**H*D*U.
!>          If UPLO = 'L', the (n-1) subdiagonal elements of the unit
!>          bidiagonal factor L from the factorization A = L*D*L**H.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the right hand side vectors B for the system of
!>          linear equations.
!>          On exit, the solution vectors, X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
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
!> \date June 2016
!
!> \ingroup complex16PTcomputational
!
!  =====================================================================
      SUBROUTINE ZPTTRS(Uplo,N,Nrhs,D,E,B,Ldb,Info)
      USE F77KINDS                        
      USE S_ILAENV
      USE S_XERBLA
      USE S_ZPTTS2
      IMPLICIT NONE
!*--ZPTTRS129
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL(R8KIND) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: iuplo , j , jb , nb
      LOGICAL :: upper
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      Info = 0
      upper = (Uplo=='U' .OR. Uplo=='u')
      IF ( .NOT.upper .AND. .NOT.(Uplo=='L' .OR. Uplo=='l') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZPTTRS',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
!     Determine the number of right-hand sides to solve at a time.
!
      IF ( Nrhs==1 ) THEN
         nb = 1
      ELSE
         nb = MAX(1,ILAENV(1,'ZPTTRS',Uplo,N,Nrhs,-1,-1))
      ENDIF
!
!     Decode UPLO
!
      IF ( upper ) THEN
         iuplo = 1
      ELSE
         iuplo = 0
      ENDIF
!
      IF ( nb>=Nrhs ) THEN
         CALL ZPTTS2(iuplo,N,Nrhs,D,E,B,Ldb)
      ELSE
         DO j = 1 , Nrhs , nb
            jb = MIN(Nrhs-j+1,nb)
            CALL ZPTTS2(iuplo,N,jb,D,E,B(1,j),Ldb)
         ENDDO
      ENDIF
!
!
!     End of ZPTTRS
!
      END SUBROUTINE ZPTTRS
