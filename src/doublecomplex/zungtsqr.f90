!*==zungtsqr.f90  processed by SPAG 7.51RB at 20:09 on  3 Mar 2022
!> \brief \b ZUNGTSQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNGTSQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zuntsqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungtsqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungtsqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNGTSQR( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK,
!      $                     INFO )
!
!       .. Scalar Arguments ..
!       INTEGER           INFO, LDA, LDT, LWORK, M, N, MB, NB
!       ..
!       .. Array Arguments ..
!       COMPLEX*16        A( LDA, * ), T( LDT, * ), WORK( * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNGTSQR generates an M-by-N complex matrix Q_out with orthonormal
!> columns, which are the first N columns of a product of comlpex unitary
!> matrices of order M which are returned by ZLATSQR
!>
!>      Q_out = first_N_columns_of( Q(1)_in * Q(2)_in * ... * Q(k)_in ).
!>
!> See the documentation for ZLATSQR.
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
!>          The number of columns of the matrix A. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] MB
!> \verbatim
!>          MB is INTEGER
!>          The row block size used by DLATSQR to return
!>          arrays A and T. MB > N.
!>          (Note that if MB > M, then M is used instead of MB
!>          as the row block size).
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The column block size used by ZLATSQR to return
!>          arrays A and T. NB >= 1.
!>          (Note that if NB > N, then N is used instead of NB
!>          as the column block size).
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>
!>          On entry:
!>
!>             The elements on and above the diagonal are not accessed.
!>             The elements below the diagonal represent the unit
!>             lower-trapezoidal blocked matrix V computed by ZLATSQR
!>             that defines the input matrices Q_in(k) (ones on the
!>             diagonal are not stored) (same format as the output A
!>             below the diagonal in ZLATSQR).
!>
!>          On exit:
!>
!>             The array A contains an M-by-N orthonormal matrix Q_out,
!>             i.e the columns of A are orthogonal unit vectors.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is COMPLEX*16 array,
!>          dimension (LDT, N * NIRB)
!>          where NIRB = Number_of_input_row_blocks
!>                     = MAX( 1, CEIL((M-N)/(MB-N)) )
!>          Let NICB = Number_of_input_col_blocks
!>                   = CEIL(N/NB)
!>
!>          The upper-triangular block reflectors used to define the
!>          input matrices Q_in(k), k=(1:NIRB*NICB). The block
!>          reflectors are stored in compact form in NIRB block
!>          reflector sequences. Each of NIRB block reflector sequences
!>          is stored in a larger NB-by-N column block of T and consists
!>          of NICB smaller NB-by-NB upper-triangular column blocks.
!>          (same format as the output T in ZLATSQR).
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.
!>          LDT >= max(1,min(NB1,N)).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          (workspace) COMPLEX*16 array, dimension (MAX(2,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          The dimension of the array WORK.  LWORK >= (M+NB)*N.
!>          If LWORK = -1, then a workspace query is assumed.
!>          The routine only calculates the optimal size of the WORK
!>          array, returns this value as the first entry of the WORK
!>          array, and no error message related to LWORK is issued
!>          by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!>
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2019
!
!> \ingroup complex16OTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!> November 2019, Igor Kozachenko,
!>                Computer Science Division,
!>                University of California, Berkeley
!>
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE ZUNGTSQR(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
!*--ZUNGTSQR178
!
!  -- LAPACK computational routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldt , Lwork , M , N , Mb , Nb
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , T(Ldt,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 CONE , CZERO
      PARAMETER (CONE=(1.0D+0,0.0D+0),CZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL lquery
      INTEGER iinfo , ldc , lworkopt , lc , lw , nblocal , j
!     ..
!     .. External Subroutines ..
      EXTERNAL ZCOPY , ZLAMTSQR , ZLASET , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      lquery = Lwork== - 1
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 .OR. M<N ) THEN
         Info = -2
      ELSEIF ( Mb<=N ) THEN
         Info = -3
      ELSEIF ( Nb<1 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -6
      ELSEIF ( Ldt<MAX(1,MIN(Nb,N)) ) THEN
         Info = -8
!
!        Test the input LWORK for the dimension of the array WORK.
!        This workspace is used to store array C(LDC, N) and WORK(LWORK)
!        in the call to ZLAMTSQR. See the documentation for ZLAMTSQR.
!
      ELSEIF ( Lwork<2 .AND. (.NOT.lquery) ) THEN
         Info = -10
      ELSE
!
!           Set block size for column blocks
!
         nblocal = MIN(Nb,N)
!
!           LWORK = -1, then set the size for the array C(LDC,N)
!           in ZLAMTSQR call and set the optimal size of the work array
!           WORK(LWORK) in ZLAMTSQR call.
!
         ldc = M
         lc = ldc*N
         lw = N*nblocal
!
         lworkopt = lc + lw
!
         IF ( (Lwork<MAX(1,lworkopt)) .AND. (.NOT.lquery) ) Info = -10
!
      ENDIF
!
!     Handle error in the input parameters and return workspace query.
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZUNGTSQR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         Work(1) = DCMPLX(lworkopt)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N)==0 ) THEN
         Work(1) = DCMPLX(lworkopt)
         RETURN
      ENDIF
!
!     (1) Form explicitly the tall-skinny M-by-N left submatrix Q1_in
!     of M-by-M orthogonal matrix Q_in, which is implicitly stored in
!     the subdiagonal part of input array A and in the input array T.
!     Perform by the following operation using the routine ZLAMTSQR.
!
!         Q1_in = Q_in * ( I ), where I is a N-by-N identity matrix,
!                        ( 0 )        0 is a (M-N)-by-N zero matrix.
!
!     (1a) Form M-by-N matrix in the array WORK(1:LDC*N) with ones
!     on the diagonal and zeros elsewhere.
!
      CALL ZLASET('F',M,N,CZERO,CONE,Work,ldc)
!
!     (1b)  On input, WORK(1:LDC*N) stores ( I );
!                                          ( 0 )
!
!           On output, WORK(1:LDC*N) stores Q1_in.
!
      CALL ZLAMTSQR('L','N',M,N,N,Mb,nblocal,A,Lda,T,Ldt,Work,ldc,      &
     &              Work(lc+1),lw,iinfo)
!
!     (2) Copy the result from the part of the work array (1:M,1:N)
!     with the leading dimension LDC that starts at WORK(1) into
!     the output array A(1:M,1:N) column-by-column.
!
      DO j = 1 , N
         CALL ZCOPY(M,Work((j-1)*ldc+1),1,A(1,j),1)
      ENDDO
!
      Work(1) = DCMPLX(lworkopt)
!
!     End of ZUNGTSQR
!
      END SUBROUTINE ZUNGTSQR