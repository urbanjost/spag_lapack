!*==zhetrs_aa_2stage.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZHETRS_AA_2STAGE
!
! @generated from SRC/dsytrs_aa_2stage.f, fortran d -> c, Mon Oct 30 11:59:02 2017
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETRS_AA_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs_aa_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs_aa_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs_aa_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!      SUBROUTINE ZHETRS_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, IPIV,
!                                   IPIV2, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N, NRHS, LDA, LTB, LDB, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IPIV2( * )
!       COMPLEX*16         A( LDA, * ), TB( * ), B( LDB, * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHETRS_AA_2STAGE solves a system of linear equations A*X = B with a
!> hermitian matrix A using the factorization A = U**H*T*U or
!> A = L*T*L**H computed by ZHETRF_AA_2STAGE.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are stored
!>          as an upper or lower triangular matrix.
!>          = 'U':  Upper triangular, form is A = U**H*T*U;
!>          = 'L':  Lower triangular, form is A = L*T*L**H.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          Details of factors computed by ZHETRF_AA_2STAGE.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TB
!> \verbatim
!>          TB is COMPLEX*16 array, dimension (LTB)
!>          Details of factors computed by ZHETRF_AA_2STAGE.
!> \endverbatim
!>
!> \param[in] LTB
!> \verbatim
!>          LTB is INTEGER
!>          The size of the array TB. LTB >= 4*N.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges as computed by
!>          ZHETRF_AA_2STAGE.
!> \endverbatim
!>
!> \param[in] IPIV2
!> \verbatim
!>          IPIV2 is INTEGER array, dimension (N)
!>          Details of the interchanges as computed by
!>          ZHETRF_AA_2STAGE.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the right hand side matrix B.
!>          On exit, the solution matrix X.
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
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complex16SYcomputational
!
!  =====================================================================
      SUBROUTINE ZHETRS_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B,&
     &                            Ldb,Info)
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZGBTRS
      USE S_ZLASWP
      USE S_ZTRSM
      IMPLICIT NONE
!*--ZHETRS_AA_2STAGE157
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tb
      INTEGER , INTENT(IN) :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: ldtb , nb
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
!     ..
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
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ltb<(4*N) ) THEN
         Info = -7
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -11
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHETRS_AA_2STAGE',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
!     Read NB and compute LDTB
!
      nb = INT(Tb(1))
      ldtb = Ltb/N
!
      IF ( upper ) THEN
!
!        Solve A*X = B, where A = U**H*T*U.
!
         IF ( N>nb ) THEN
!
!           Pivot, P**T * B -> B
!
            CALL ZLASWP(Nrhs,B,Ldb,nb+1,N,Ipiv,1)
!
!           Compute (U**H \ B) -> B    [ (U**H \P**T * B) ]
!
            CALL ZTRSM('L','U','C','U',N-nb,Nrhs,ONE,A(1,nb+1),Lda,     &
     &                 B(nb+1,1),Ldb)
!
         ENDIF
!
!        Compute T \ B -> B   [ T \ (U**H \P**T * B) ]
!
         CALL ZGBTRS('N',N,nb,nb,Nrhs,Tb,ldtb,Ipiv2,B,Ldb,Info)
         IF ( N>nb ) THEN
!
!           Compute (U \ B) -> B   [ U \ (T \ (U**H \P**T * B) ) ]
!
            CALL ZTRSM('L','U','N','U',N-nb,Nrhs,ONE,A(1,nb+1),Lda,     &
     &                 B(nb+1,1),Ldb)
!
!           Pivot, P * B -> B  [ P * (U \ (T \ (U**H \P**T * B) )) ]
!
            CALL ZLASWP(Nrhs,B,Ldb,nb+1,N,Ipiv,-1)
!
         ENDIF
!
      ELSE
!
!        Solve A*X = B, where A = L*T*L**H.
!
         IF ( N>nb ) THEN
!
!           Pivot, P**T * B -> B
!
            CALL ZLASWP(Nrhs,B,Ldb,nb+1,N,Ipiv,1)
!
!           Compute (L \ B) -> B    [ (L \P**T * B) ]
!
            CALL ZTRSM('L','L','N','U',N-nb,Nrhs,ONE,A(nb+1,1),Lda,     &
     &                 B(nb+1,1),Ldb)
!
         ENDIF
!
!        Compute T \ B -> B   [ T \ (L \P**T * B) ]
!
         CALL ZGBTRS('N',N,nb,nb,Nrhs,Tb,ldtb,Ipiv2,B,Ldb,Info)
         IF ( N>nb ) THEN
!
!           Compute (L**H \ B) -> B   [ L**H \ (T \ (L \P**T * B) ) ]
!
            CALL ZTRSM('L','L','C','U',N-nb,Nrhs,ONE,A(nb+1,1),Lda,     &
     &                 B(nb+1,1),Ldb)
!
!           Pivot, P * B -> B  [ P * (L**H \ (T \ (L \P**T * B) )) ]
!
            CALL ZLASWP(Nrhs,B,Ldb,nb+1,N,Ipiv,-1)
!
         ENDIF
      ENDIF
!
!
!     End of ZHETRS_AA_2STAGE
!
      END SUBROUTINE ZHETRS_AA_2STAGE
