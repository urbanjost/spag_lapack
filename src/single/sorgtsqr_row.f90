!*==sorgtsqr_row.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SORGTSQR_ROW
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SORGTSQR_ROW + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgtsqr_row.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgtsqr_row.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgtsqr_row.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SORGTSQR_ROW( M, N, MB, NB, A, LDA, T, LDT, WORK,
!      $                         LWORK, INFO )
!       IMPLICIT NONE
!
!       .. Scalar Arguments ..
!       INTEGER           INFO, LDA, LDT, LWORK, M, N, MB, NB
!       ..
!       .. Array Arguments ..
!       REAL              A( LDA, * ), T( LDT, * ), WORK( * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SORGTSQR_ROW generates an M-by-N real matrix Q_out with
!> orthonormal columns from the output of SLATSQR. These N orthonormal
!> columns are the first N columns of a product of complex unitary
!> matrices Q(k)_in of order M, which are returned by SLATSQR in
!> a special format.
!>
!>      Q_out = first_N_columns_of( Q(1)_in * Q(2)_in * ... * Q(k)_in ).
!>
!> The input matrices Q(k)_in are stored in row and column blocks in A.
!> See the documentation of SLATSQR for more details on the format of
!> Q(k)_in, where each Q(k)_in is represented by block Householder
!> transformations. This routine calls an auxiliary routine SLARFB_GETT,
!> where the computation is performed on each individual block. The
!> algorithm first sweeps NB-sized column blocks from the right to left
!> starting in the bottom row block and continues to the top row block
!> (hence _ROW in the routine name). This sweep is in reverse order of
!> the order in which SLATSQR generates the output blocks.
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
!>          The row block size used by SLATSQR to return
!>          arrays A and T. MB > N.
!>          (Note that if MB > M, then M is used instead of MB
!>          as the row block size).
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The column block size used by SLATSQR to return
!>          arrays A and T. NB >= 1.
!>          (Note that if NB > N, then N is used instead of NB
!>          as the column block size).
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>
!>          On entry:
!>
!>             The elements on and above the diagonal are not used as
!>             input. The elements below the diagonal represent the unit
!>             lower-trapezoidal blocked matrix V computed by SLATSQR
!>             that defines the input matrices Q_in(k) (ones on the
!>             diagonal are not stored). See SLATSQR for more details.
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
!>          T is REAL array,
!>          dimension (LDT, N * NIRB)
!>          where NIRB = Number_of_input_row_blocks
!>                     = MAX( 1, CEIL((M-N)/(MB-N)) )
!>          Let NICB = Number_of_input_col_blocks
!>                   = CEIL(N/NB)
!>
!>          The upper-triangular block reflectors used to define the
!>          input matrices Q_in(k), k=(1:NIRB*NICB). The block
!>          reflectors are stored in compact form in NIRB block
!>          reflector sequences. Each of the NIRB block reflector
!>          sequences is stored in a larger NB-by-N column block of T
!>          and consists of NICB smaller NB-by-NB upper-triangular
!>          column blocks. See SLATSQR for more details on the format
!>          of T.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.
!>          LDT >= max(1,min(NB,N)).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          (workspace) REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          The dimension of the array WORK.
!>          LWORK >= NBLOCAL * MAX(NBLOCAL,(N-NBLOCAL)),
!>          where NBLOCAL=MIN(NB,N).
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
!> \date November 2020
!
!> \ingroup sigleOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!> November 2020, Igor Kozachenko,
!>                Computer Science Division,
!>                University of California, Berkeley
!>
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SORGTSQR_ROW(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
!*--SORGTSQR_ROW191
!
!  -- LAPACK computational routine (version 3.10.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2020
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldt , Lwork , M , N , Mb , Nb
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , T(Ldt,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery
      INTEGER nblocal , mb2 , m_plus_one , itmp , ib_bottom , lworkopt ,&
     &        num_all_row_blocks , jb_t , ib , imb , kb , kb_last ,     &
     &        knb , mb1
!     ..
!     .. Local Arrays ..
      REAL dummy(1,1)
!     ..
!     .. External Subroutines ..
      EXTERNAL SLARFB_GETT , SLASET , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC REAL , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      Info = 0
      lquery = Lwork== - 1
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
      ELSEIF ( Lwork<1 .AND. .NOT.lquery ) THEN
         Info = -10
      ENDIF
!
      nblocal = MIN(Nb,N)
!
!     Determine the workspace size.
!
      IF ( Info==0 ) lworkopt = nblocal*MAX(nblocal,(N-nblocal))
!
!     Handle error in the input parameters and handle the workspace query.
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SORGTSQR_ROW',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         Work(1) = REAL(lworkopt)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N)==0 ) THEN
         Work(1) = REAL(lworkopt)
         RETURN
      ENDIF
!
!     (0) Set the upper-triangular part of the matrix A to zero and
!     its diagonal elements to one.
!
      CALL SLASET('U',M,N,ZERO,ONE,A,Lda)
!
!     KB_LAST is the column index of the last column block reflector
!     in the matrices T and V.
!
      kb_last = ((N-1)/nblocal)*nblocal + 1
!
!
!     (1) Bottom-up loop over row blocks of A, except the top row block.
!     NOTE: If MB>=M, then the loop is never executed.
!
      IF ( Mb<M ) THEN
!
!        MB2 is the row blocking size for the row blocks before the
!        first top row block in the matrix A. IB is the row index for
!        the row blocks in the matrix A before the first top row block.
!        IB_BOTTOM is the row index for the last bottom row block
!        in the matrix A. JB_T is the column index of the corresponding
!        column block in the matrix T.
!
!        Initialize variables.
!
!        NUM_ALL_ROW_BLOCKS is the number of row blocks in the matrix A
!        including the first row block.
!
         mb2 = Mb - N
         m_plus_one = M + 1
         itmp = (M-Mb-1)/mb2
         ib_bottom = itmp*mb2 + Mb + 1
         num_all_row_blocks = itmp + 2
         jb_t = num_all_row_blocks*N + 1
!
         DO ib = ib_bottom , Mb + 1 , -mb2
!
!           Determine the block size IMB for the current row block
!           in the matrix A.
!
            imb = MIN(m_plus_one-ib,mb2)
!
!           Determine the column index JB_T for the current column block
!           in the matrix T.
!
            jb_t = jb_t - N
!
!           Apply column blocks of H in the row block from right to left.
!
!           KB is the column index of the current column block reflector
!           in the matrices T and V.
!
            DO kb = kb_last , 1 , -nblocal
!
!              Determine the size of the current column block KNB in
!              the matrices T and V.
!
               knb = MIN(nblocal,N-kb+1)
!
               CALL SLARFB_GETT('I',imb,N-kb+1,knb,T(1,jb_t+kb-1),Ldt,  &
     &                          A(kb,kb),Lda,A(ib,kb),Lda,Work,knb)
!
            ENDDO
!
         ENDDO
!
      ENDIF
!
!     (2) Top row block of A.
!     NOTE: If MB>=M, then we have only one row block of A of size M
!     and we work on the entire matrix A.
!
      mb1 = MIN(Mb,M)
!
!     Apply column blocks of H in the top row block from right to left.
!
!     KB is the column index of the current block reflector in
!     the matrices T and V.
!
      DO kb = kb_last , 1 , -nblocal
!
!        Determine the size of the current column block KNB in
!        the matrices T and V.
!
         knb = MIN(nblocal,N-kb+1)
!
         IF ( mb1-kb-knb+1==0 ) THEN
!
!           In SLARFB_GETT parameters, when M=0, then the matrix B
!           does not exist, hence we need to pass a dummy array
!           reference DUMMY(1,1) to B with LDDUMMY=1.
!
            CALL SLARFB_GETT('N',0,N-kb+1,knb,T(1,kb),Ldt,A(kb,kb),Lda, &
     &                       dummy(1,1),1,Work,knb)
         ELSE
            CALL SLARFB_GETT('N',mb1-kb-knb+1,N-kb+1,knb,T(1,kb),Ldt,   &
     &                       A(kb,kb),Lda,A(kb+knb,kb),Lda,Work,knb)
 
         ENDIF
!
      ENDDO
!
      Work(1) = REAL(lworkopt)
!
!     End of SORGTSQR_ROW
!
      END SUBROUTINE SORGTSQR_ROW
