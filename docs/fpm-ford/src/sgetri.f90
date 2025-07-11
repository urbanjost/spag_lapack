!*==sgetri.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SGETRI
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGETRI + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgetri.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgetri.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgetri.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGETRI computes the inverse of a matrix using the LU factorization
!> computed by SGETRF.
!>
!> This method inverts U and then computes inv(A) by solving the system
!> inv(A)*L = inv(U) for inv(A).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the factors L and U from the factorization
!>          A = P*L*U as computed by SGETRF.
!>          On exit, if INFO = 0, the inverse of the original matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from SGETRF; for 1<=i<=N, row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,N).
!>          For optimal performance LWORK >= N*NB, where NB is
!>          the optimal blocksize returned by ILAENV.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!>                singular and its inverse could not be computed.
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
!> \ingroup realGEcomputational
!
!  =====================================================================
      SUBROUTINE SGETRI(N,A,Lda,Ipiv,Work,Lwork,Info)
      IMPLICIT NONE
!*--SGETRI118
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Lwork , N
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      REAL A(Lda,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery
      INTEGER i , iws , j , jb , jj , jp , ldwork , lwkopt , nb ,       &
     &        nbmin , nn
!     ..
!     .. External Functions ..
      INTEGER ILAENV
      EXTERNAL ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEMM , SGEMV , SSWAP , STRSM , STRTRI , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      nb = ILAENV(1,'SGETRI',' ',N,-1,-1,-1)
      lwkopt = N*nb
      Work(1) = lwkopt
      lquery = (Lwork==-1)
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -3
      ELSEIF ( Lwork<MAX(1,N) .AND. .NOT.lquery ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGETRI',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Form inv(U).  If INFO > 0 from STRTRI, then U is singular,
!     and the inverse is not computed.
!
      CALL STRTRI('Upper','Non-unit',N,A,Lda,Info)
      IF ( Info>0 ) RETURN
!
      nbmin = 2
      ldwork = N
      IF ( nb>1 .AND. nb<N ) THEN
         iws = MAX(ldwork*nb,1)
         IF ( Lwork<iws ) THEN
            nb = Lwork/ldwork
            nbmin = MAX(2,ILAENV(2,'SGETRI',' ',N,-1,-1,-1))
         ENDIF
      ELSE
         iws = N
      ENDIF
!
!     Solve the equation inv(A)*L = inv(U) for inv(A).
!
      IF ( nb<nbmin .OR. nb>=N ) THEN
!
!        Use unblocked code.
!
         DO j = N , 1 , -1
!
!           Copy current column of L to WORK and replace with zeros.
!
            DO i = j + 1 , N
               Work(i) = A(i,j)
               A(i,j) = ZERO
            ENDDO
!
!           Compute current column of inv(A).
!
            IF ( j<N ) CALL SGEMV('No transpose',N,N-j,-ONE,A(1,j+1),   &
     &                            Lda,Work(j+1),1,ONE,A(1,j),1)
         ENDDO
      ELSE
!
!        Use blocked code.
!
         nn = ((N-1)/nb)*nb + 1
         DO j = nn , 1 , -nb
            jb = MIN(nb,N-j+1)
!
!           Copy current block column of L to WORK and replace with
!           zeros.
!
            DO jj = j , j + jb - 1
               DO i = jj + 1 , N
                  Work(i+(jj-j)*ldwork) = A(i,jj)
                  A(i,jj) = ZERO
               ENDDO
            ENDDO
!
!           Compute current block column of inv(A).
!
            IF ( j+jb<=N ) CALL SGEMM('No transpose','No transpose',N,  &
     &                                jb,N-j-jb+1,-ONE,A(1,j+jb),Lda,   &
     &                                Work(j+jb),ldwork,ONE,A(1,j),Lda)
            CALL STRSM('Right','Lower','No transpose','Unit',N,jb,ONE,  &
     &                 Work(j),ldwork,A(1,j),Lda)
         ENDDO
      ENDIF
!
!     Apply column interchanges.
!
      DO j = N - 1 , 1 , -1
         jp = Ipiv(j)
         IF ( jp/=j ) CALL SSWAP(N,A(1,j),1,A(1,jp),1)
      ENDDO
!
      Work(1) = iws
!
!     End of SGETRI
!
      END SUBROUTINE SGETRI
