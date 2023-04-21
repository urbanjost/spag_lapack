!*==dpftrs.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DPFTRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DPFTRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpftrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpftrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpftrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DPFTRS( TRANSR, UPLO, N, NRHS, A, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANSR, UPLO
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( 0: * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DPFTRS solves a system of linear equations A*X = B with a symmetric
!> positive definite matrix A using the Cholesky factorization
!> A = U**T*U or A = L*L**T computed by DPFTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSR
!> \verbatim
!>          TRANSR is CHARACTER*1
!>          = 'N':  The Normal TRANSR of RFP A is stored;
!>          = 'T':  The Transpose TRANSR of RFP A is stored.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of RFP A is stored;
!>          = 'L':  Lower triangle of RFP A is stored.
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
!>          A is DOUBLE PRECISION array, dimension ( N*(N+1)/2 ).
!>          The triangular factor U or L from the Cholesky factorization
!>          of RFP A = U**T*U or RFP A = L*L**T, as computed by DPFTRF.
!>          See note below for more details about RFP A.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
!> \date December 2016
!
!> \ingroup doubleOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  We first consider Rectangular Full Packed (RFP) Format when N is
!>  even. We give an example where N = 6.
!>
!>      AP is Upper             AP is Lower
!>
!>   00 01 02 03 04 05       00
!>      11 12 13 14 15       10 11
!>         22 23 24 25       20 21 22
!>            33 34 35       30 31 32 33
!>               44 45       40 41 42 43 44
!>                  55       50 51 52 53 54 55
!>
!>
!>  Let TRANSR = 'N'. RFP holds AP as follows:
!>  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last
!>  three columns of AP upper. The lower triangle A(4:6,0:2) consists of
!>  the transpose of the first three columns of AP upper.
!>  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first
!>  three columns of AP lower. The upper triangle A(0:2,0:2) consists of
!>  the transpose of the last three columns of AP lower.
!>  This covers the case N even and TRANSR = 'N'.
!>
!>         RFP A                   RFP A
!>
!>        03 04 05                33 43 53
!>        13 14 15                00 44 54
!>        23 24 25                10 11 55
!>        33 34 35                20 21 22
!>        00 44 45                30 31 32
!>        01 11 55                40 41 42
!>        02 12 22                50 51 52
!>
!>  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the
!>  transpose of RFP A above. One therefore gets:
!>
!>
!>           RFP A                   RFP A
!>
!>     03 13 23 33 00 01 02    33 00 10 20 30 40 50
!>     04 14 24 34 44 11 12    43 44 11 21 31 41 51
!>     05 15 25 35 45 55 22    53 54 55 22 32 42 52
!>
!>
!>  We then consider Rectangular Full Packed (RFP) Format when N is
!>  odd. We give an example where N = 5.
!>
!>     AP is Upper                 AP is Lower
!>
!>   00 01 02 03 04              00
!>      11 12 13 14              10 11
!>         22 23 24              20 21 22
!>            33 34              30 31 32 33
!>               44              40 41 42 43 44
!>
!>
!>  Let TRANSR = 'N'. RFP holds AP as follows:
!>  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last
!>  three columns of AP upper. The lower triangle A(3:4,0:1) consists of
!>  the transpose of the first two columns of AP upper.
!>  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first
!>  three columns of AP lower. The upper triangle A(0:1,1:2) consists of
!>  the transpose of the last two columns of AP lower.
!>  This covers the case N odd and TRANSR = 'N'.
!>
!>         RFP A                   RFP A
!>
!>        02 03 04                00 33 43
!>        12 13 14                10 11 44
!>        22 23 24                20 21 22
!>        00 33 34                30 31 32
!>        01 11 44                40 41 42
!>
!>  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the
!>  transpose of RFP A above. One therefore gets:
!>
!>           RFP A                   RFP A
!>
!>     02 12 22 00 01             00 10 20 30 40 50
!>     03 13 23 33 11             33 11 21 31 41 51
!>     04 14 24 34 44             43 44 22 32 42 52
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DPFTRS(Transr,Uplo,N,Nrhs,A,B,Ldb,Info)
      IMPLICIT NONE
!*--DPFTRS203
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Transr , Uplo
      INTEGER Info , Ldb , N , Nrhs
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(0:*) , B(Ldb,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lower , normaltransr
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , DTFSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      normaltransr = LSAME(Transr,'N')
      lower = LSAME(Uplo,'L')
      IF ( .NOT.normaltransr .AND. .NOT.LSAME(Transr,'T') ) THEN
         Info = -1
      ELSEIF ( .NOT.lower .AND. .NOT.LSAME(Uplo,'U') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Nrhs<0 ) THEN
         Info = -4
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DPFTRS',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 ) RETURN
!
!     start execution: there are two triangular solves
!
      IF ( lower ) THEN
         CALL DTFSM(Transr,'L',Uplo,'N','N',N,Nrhs,ONE,A,B,Ldb)
         CALL DTFSM(Transr,'L',Uplo,'T','N',N,Nrhs,ONE,A,B,Ldb)
      ELSE
         CALL DTFSM(Transr,'L',Uplo,'T','N',N,Nrhs,ONE,A,B,Ldb)
         CALL DTFSM(Transr,'L',Uplo,'N','N',N,Nrhs,ONE,A,B,Ldb)
      ENDIF
!
!
!     End of DPFTRS
!
      END SUBROUTINE DPFTRS
