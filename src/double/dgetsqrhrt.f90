!*==dgetsqrhrt.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DGETSQRHRT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGETSQRHRT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetsqrhrt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetsqrhrt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetsqrhrt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGETSQRHRT( M, N, MB1, NB1, NB2, A, LDA, T, LDT, WORK,
!      $                       LWORK, INFO )
!       IMPLICIT NONE
!
!       .. Scalar Arguments ..
!       INTEGER           INFO, LDA, LDT, LWORK, M, N, NB1, NB2, MB1
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION  A( LDA, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGETSQRHRT computes a NB2-sized column blocked QR-factorization
!> of a real M-by-N matrix A with M >= N,
!>
!>    A = Q * R.
!>
!> The routine uses internally a NB1-sized column blocked and MB1-sized
!> row blocked TSQR-factorization and perfors the reconstruction
!> of the Householder vectors from the TSQR output. The routine also
!> converts the R_tsqr factor from the TSQR-factorization output into
!> the R factor that corresponds to the Householder QR-factorization,
!>
!>    A = Q_tsqr * R_tsqr = Q * R.
!>
!> The output Q and R factors are stored in the same format as in DGEQRT
!> (Q is in blocked compact WY-representation). See the documentation
!> of DGEQRT for more details on the format.
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
!> \param[in] MB1
!> \verbatim
!>          MB1 is INTEGER
!>          The row block size to be used in the blocked TSQR.
!>          MB1 > N.
!> \endverbatim
!>
!> \param[in] NB1
!> \verbatim
!>          NB1 is INTEGER
!>          The column block size to be used in the blocked TSQR.
!>          N >= NB1 >= 1.
!> \endverbatim
!>
!> \param[in] NB2
!> \verbatim
!>          NB2 is INTEGER
!>          The block size to be used in the blocked QR that is
!>          output. NB2 >= 1.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>
!>          On entry: an M-by-N matrix A.
!>
!>          On exit:
!>           a) the elements on and above the diagonal
!>              of the array contain the N-by-N upper-triangular
!>              matrix R corresponding to the Householder QR;
!>           b) the elements below the diagonal represent Q by
!>              the columns of blocked V (compact WY-representation).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,N))
!>          The upper triangular block reflectors stored in compact form
!>          as a sequence of upper triangular blocks.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= NB2.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          The dimension of the array WORK.
!>          LWORK >= MAX( LWT + LW1, MAX( LWT+N*N+LW2, LWT+N*N+N ) ),
!>          where
!>             NUM_ALL_ROW_BLOCKS = CEIL((M-N)/(MB1-N)),
!>             NB1LOCAL = MIN(NB1,N).
!>             LWT = NUM_ALL_ROW_BLOCKS * N * NB1LOCAL,
!>             LW1 = NB1LOCAL * N,
!>             LW2 = NB1LOCAL * MAX( NB1LOCAL, ( N - NB1LOCAL ) ),
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
!
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
!> \ingroup doubleOTHERcomputational
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
      SUBROUTINE DGETSQRHRT(M,N,Mb1,Nb1,Nb2,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
!*--DGETSQRHRT182
!
!  -- LAPACK computational routine (version 3.10.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2020
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldt , Lwork , M , N , Nb1 , Nb2 , Mb1
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , T(Ldt,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery
      INTEGER i , iinfo , j , lw1 , lw2 , lwt , ldwt , lworkopt ,       &
     &        nb1local , nb2local , num_all_row_blocks
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DLATSQR , DORGTSQR_ROW , DORHR_COL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CEILING , DBLE , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      lquery = Lwork== - 1
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 .OR. M<N ) THEN
         Info = -2
      ELSEIF ( Mb1<=N ) THEN
         Info = -3
      ELSEIF ( Nb1<1 ) THEN
         Info = -4
      ELSEIF ( Nb2<1 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -7
      ELSEIF ( Ldt<MAX(1,MIN(Nb2,N)) ) THEN
         Info = -9
!
!        Test the input LWORK for the dimension of the array WORK.
!        This workspace is used to store array:
!        a) Matrix T and WORK for DLATSQR;
!        b) N-by-N upper-triangular factor R_tsqr;
!        c) Matrix T and array WORK for DORGTSQR_ROW;
!        d) Diagonal D for DORHR_COL.
!
      ELSEIF ( Lwork<N*N+1 .AND. .NOT.lquery ) THEN
         Info = -11
      ELSE
!
!           Set block size for column blocks
!
         nb1local = MIN(Nb1,N)
!
         num_all_row_blocks = MAX(1,CEILING(DBLE(M-N)/DBLE(Mb1-N)))
!
!           Length and leading dimension of WORK array to place
!           T array in TSQR.
!
         lwt = num_all_row_blocks*N*nb1local
 
         ldwt = nb1local
!
!           Length of TSQR work array
!
         lw1 = nb1local*N
!
!           Length of DORGTSQR_ROW work array.
!
         lw2 = nb1local*MAX(nb1local,(N-nb1local))
!
         lworkopt = MAX(lwt+lw1,MAX(lwt+N*N+lw2,lwt+N*N+N))
!
         IF ( (Lwork<MAX(1,lworkopt)) .AND. (.NOT.lquery) ) Info = -11
!
      ENDIF
!
!     Handle error in the input parameters and return workspace query.
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGETSQRHRT',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         Work(1) = DBLE(lworkopt)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N)==0 ) THEN
         Work(1) = DBLE(lworkopt)
         RETURN
      ENDIF
!
      nb2local = MIN(Nb2,N)
!
!
!     (1) Perform TSQR-factorization of the M-by-N matrix A.
!
      CALL DLATSQR(M,N,Mb1,nb1local,A,Lda,Work,ldwt,Work(lwt+1),lw1,    &
     &             iinfo)
!
!     (2) Copy the factor R_tsqr stored in the upper-triangular part
!         of A into the square matrix in the work array
!         WORK(LWT+1:LWT+N*N) column-by-column.
!
      DO j = 1 , N
         CALL DCOPY(j,A(1,j),1,Work(lwt+N*(j-1)+1),1)
      ENDDO
!
!     (3) Generate a M-by-N matrix Q with orthonormal columns from
!     the result stored below the diagonal in the array A in place.
!
 
      CALL DORGTSQR_ROW(M,N,Mb1,nb1local,A,Lda,Work,ldwt,Work(lwt+N*N+1)&
     &                  ,lw2,iinfo)
!
!     (4) Perform the reconstruction of Householder vectors from
!     the matrix Q (stored in A) in place.
!
      CALL DORHR_COL(M,N,nb2local,A,Lda,T,Ldt,Work(lwt+N*N+1),iinfo)
!
!     (5) Copy the factor R_tsqr stored in the square matrix in the
!     work array WORK(LWT+1:LWT+N*N) into the upper-triangular
!     part of A.
!
!     (6) Compute from R_tsqr the factor R_hr corresponding to
!     the reconstructed Householder vectors, i.e. R_hr = S * R_tsqr.
!     This multiplication by the sign matrix S on the left means
!     changing the sign of I-th row of the matrix R_tsqr according
!     to sign of the I-th diagonal element DIAG(I) of the matrix S.
!     DIAG is stored in WORK( LWT+N*N+1 ) from the DORHR_COL output.
!
!     (5) and (6) can be combined in a single loop, so the rows in A
!     are accessed only once.
!
      DO i = 1 , N
         IF ( Work(lwt+N*N+i)==-ONE ) THEN
            DO j = i , N
               A(i,j) = -ONE*Work(lwt+N*(j-1)+i)
            ENDDO
         ELSE
            CALL DCOPY(N-i+1,Work(lwt+N*(i-1)+i),N,A(i,i),Lda)
         ENDIF
      ENDDO
!
      Work(1) = DBLE(lworkopt)
!
!     End of DGETSQRHRT
!
      END SUBROUTINE DGETSQRHRT
