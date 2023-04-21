!*==clahef_aa.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLAHEF_AA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAHEF_AA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahef_aa.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahef_aa.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahef_aa.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAHEF_AA( UPLO, J1, M, NB, A, LDA, IPIV,
!                             H, LDH, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER    UPLO
!       INTEGER      J1, M, NB, LDA, LDH
!       ..
!       .. Array Arguments ..
!       INTEGER      IPIV( * )
!       COMPLEX      A( LDA, * ), H( LDH, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAHEF_AA factorizes a panel of a complex hermitian matrix A using
!> the Aasen's algorithm. The panel consists of a set of NB rows of A
!> when UPLO is U, or a set of NB columns when UPLO is L.
!>
!> In order to factorize the panel, the Aasen's algorithm requires the
!> last row, or column, of the previous panel. The first row, or column,
!> of A is set to be the first row, or column, of an identity matrix,
!> which is used to factorize the first panel.
!>
!> The resulting J-th row of U, or J-th column of L, is stored in the
!> (J-1)-th row, or column, of A (without the unit diagonals), while
!> the diagonal and subdiagonal of A are overwritten by those of T.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] J1
!> \verbatim
!>          J1 is INTEGER
!>          The location of the first row, or column, of the panel
!>          within the submatrix of A, passed to this routine, e.g.,
!>          when called by CHETRF_AA, for the first panel, J1 is 1,
!>          while for the remaining panels, J1 is 2.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The dimension of the submatrix. M >= 0.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The dimension of the panel to be facotorized.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,M) for
!>          the first panel, while dimension (LDA,M+1) for the
!>          remaining panels.
!>
!>          On entry, A contains the last row, or column, of
!>          the previous panel, and the trailing submatrix of A
!>          to be factorized, except for the first panel, only
!>          the panel is passed.
!>
!>          On exit, the leading panel is factorized.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the row and column interchanges,
!>          the row and column k were interchanged with the row and
!>          column IPIV(k).
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX workspace, dimension (LDH,NB).
!>
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the workspace H. LDH >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX workspace, dimension (M).
!> \endverbatim
!>
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
!> \ingroup complexSYcomputational
!
!  =====================================================================
      SUBROUTINE CLAHEF_AA(Uplo,J1,M,Nb,A,Lda,Ipiv,H,Ldh,Work)
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
      IMPLICIT NONE
!*--CLAHEF_AA153
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER M , Nb , J1 , Lda , Ldh
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX A(Lda,*) , H(Ldh,*) , Work(*)
!     ..
!
!  =====================================================================
!     .. Parameters ..
      COMPLEX ZERO , ONE
      PARAMETER (ZERO=(0.0E+0,0.0E+0),ONE=(1.0E+0,0.0E+0))
!
!     .. Local Scalars ..
      INTEGER j , k , k1 , i1 , i2 , mj
      COMPLEX piv , alpha
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ICAMAX , ILAENV
      EXTERNAL LSAME , ILAENV , ICAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL CLACGV , CGEMV , CSCAL , CAXPY , CCOPY , CSWAP , CLASET ,&
     &         XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC REAL , CONJG , MAX
!     ..
!     .. Executable Statements ..
!
      j = 1
!
!     K1 is the first column of the panel to be factorized
!     i.e.,  K1 is 2 for the first block column, and 1 for the rest of the blocks
!
      k1 = (2-J1) + 1
!
      IF ( LSAME(Uplo,'U') ) THEN
!
!        .....................................................
!        Factorize A as U**T*D*U using the upper triangle of A
!        .....................................................
!
         DO WHILE ( j<=MIN(M,Nb) )
!
!        K is the column to be factorized
!         when being called from CHETRF_AA,
!         > for the first block column, J1 is 1, hence J1+J-1 is J,
!         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1,
!
            k = J1 + j - 1
            IF ( j==M ) THEN
!
!            Only need to compute T(J, J)
!
               mj = 1
            ELSE
               mj = M - j + 1
            ENDIF
!
!        H(J:N, J) := A(J, J:N) - H(J:N, 1:(J-1)) * L(J1:(J-1), J),
!         where H(J:N, J) has been initialized to be A(J, J:N)
!
            IF ( k>2 ) THEN
!
!        K is the column to be factorized
!         > for the first block column, K is J, skipping the first two
!           columns
!         > for the rest of the columns, K is J+1, skipping only the
!           first column
!
               CALL CLACGV(j-k1,A(1,j),1)
               CALL CGEMV('No transpose',mj,j-k1,-ONE,H(j,k1),Ldh,A(1,j)&
     &                    ,1,ONE,H(j,j),1)
               CALL CLACGV(j-k1,A(1,j),1)
            ENDIF
!
!        Copy H(i:n, i) into WORK
!
            CALL CCOPY(mj,H(j,j),1,Work(1),1)
!
            IF ( j>k1 ) THEN
!
!           Compute WORK := WORK - L(J-1, J:N) * T(J-1,J),
!            where A(J-1, J) stores T(J-1, J) and A(J-2, J:N) stores U(J-1, J:N)
!
               alpha = -CONJG(A(k-1,j))
               CALL CAXPY(mj,alpha,A(k-2,j),Lda,Work(1),1)
            ENDIF
!
!        Set A(J, J) = T(J, J)
!
            A(k,j) = REAL(Work(1))
!
            IF ( j<M ) THEN
!
!           Compute WORK(2:N) = T(J, J) L(J, (J+1):N)
!            where A(J, J) stores T(J, J) and A(J-1, (J+1):N) stores U(J, (J+1):N)
!
               IF ( k>1 ) THEN
                  alpha = -A(k,j)
                  CALL CAXPY(M-j,alpha,A(k-1,j+1),Lda,Work(2),1)
               ENDIF
!
!           Find max(|WORK(2:n)|)
!
               i2 = ICAMAX(M-j,Work(2),1) + 1
               piv = Work(i2)
!
!           Apply hermitian pivot
!
               IF ( (i2/=2) .AND. (piv/=0) ) THEN
!
!              Swap WORK(I1) and WORK(I2)
!
                  i1 = 2
                  Work(i2) = Work(i1)
                  Work(i1) = piv
!
!              Swap A(I1, I1+1:N) with A(I1+1:N, I2)
!
                  i1 = i1 + j - 1
                  i2 = i2 + j - 1
                  CALL CSWAP(i2-i1-1,A(J1+i1-1,i1+1),Lda,A(J1+i1,i2),1)
                  CALL CLACGV(i2-i1,A(J1+i1-1,i1+1),Lda)
                  CALL CLACGV(i2-i1-1,A(J1+i1,i2),1)
!
!              Swap A(I1, I2+1:N) with A(I2, I2+1:N)
!
                  IF ( i2<M ) CALL CSWAP(M-i2,A(J1+i1-1,i2+1),Lda,      &
     &                 A(J1+i2-1,i2+1),Lda)
!
!              Swap A(I1, I1) with A(I2,I2)
!
                  piv = A(i1+J1-1,i1)
                  A(J1+i1-1,i1) = A(J1+i2-1,i2)
                  A(J1+i2-1,i2) = piv
!
!              Swap H(I1, 1:J1) with H(I2, 1:J1)
!
                  CALL CSWAP(i1-1,H(i1,1),Ldh,H(i2,1),Ldh)
                  Ipiv(i1) = i2
!
!
!                 Swap L(1:I1-1, I1) with L(1:I1-1, I2),
!                  skipping the first column
!
                  IF ( i1>(k1-1) ) CALL CSWAP(i1-k1+1,A(1,i1),1,A(1,i2),&
     &                 1)
               ELSE
                  Ipiv(j+1) = j + 1
               ENDIF
!
!           Set A(J, J+1) = T(J, J+1)
!
               A(k,j+1) = Work(2)
!
!
!              Copy A(J+1:N, J+1) into H(J:N, J),
!
               IF ( j<Nb ) CALL CCOPY(M-j,A(k+1,j+1),Lda,H(j+1,j+1),1)
!
!           Compute L(J+2, J+1) = WORK( 3:N ) / T(J, J+1),
!            where A(J, J+1) = T(J, J+1) and A(J+2:N, J) = L(J+2:N, J+1)
!
               IF ( j<(M-1) ) THEN
                  IF ( A(k,j+1)/=ZERO ) THEN
                     alpha = ONE/A(k,j+1)
                     CALL CCOPY(M-j-1,Work(3),1,A(k,j+2),Lda)
                     CALL CSCAL(M-j-1,alpha,A(k,j+2),Lda)
                  ELSE
                     CALL CLASET('Full',1,M-j-1,ZERO,ZERO,A(k,j+2),Lda)
                  ENDIF
               ENDIF
            ENDIF
            j = j + 1
         ENDDO
!
      ELSE
!
!        .....................................................
!        Factorize A as L*D*L**T using the lower triangle of A
!        .....................................................
!
         DO WHILE ( j<=MIN(M,Nb) )
!
!        K is the column to be factorized
!         when being called from CHETRF_AA,
!         > for the first block column, J1 is 1, hence J1+J-1 is J,
!         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1,
!
            k = J1 + j - 1
            IF ( j==M ) THEN
!
!            Only need to compute T(J, J)
!
               mj = 1
            ELSE
               mj = M - j + 1
            ENDIF
!
!        H(J:N, J) := A(J:N, J) - H(J:N, 1:(J-1)) * L(J, J1:(J-1))^T,
!         where H(J:N, J) has been initialized to be A(J:N, J)
!
            IF ( k>2 ) THEN
!
!        K is the column to be factorized
!         > for the first block column, K is J, skipping the first two
!           columns
!         > for the rest of the columns, K is J+1, skipping only the
!           first column
!
               CALL CLACGV(j-k1,A(j,1),Lda)
               CALL CGEMV('No transpose',mj,j-k1,-ONE,H(j,k1),Ldh,A(j,1)&
     &                    ,Lda,ONE,H(j,j),1)
               CALL CLACGV(j-k1,A(j,1),Lda)
            ENDIF
!
!        Copy H(J:N, J) into WORK
!
            CALL CCOPY(mj,H(j,j),1,Work(1),1)
!
            IF ( j>k1 ) THEN
!
!           Compute WORK := WORK - L(J:N, J-1) * T(J-1,J),
!            where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1)
!
               alpha = -CONJG(A(j,k-1))
               CALL CAXPY(mj,alpha,A(j,k-2),1,Work(1),1)
            ENDIF
!
!        Set A(J, J) = T(J, J)
!
            A(j,k) = REAL(Work(1))
!
            IF ( j<M ) THEN
!
!           Compute WORK(2:N) = T(J, J) L((J+1):N, J)
!            where A(J, J) = T(J, J) and A((J+1):N, J-1) = L((J+1):N, J)
!
               IF ( k>1 ) THEN
                  alpha = -A(j,k)
                  CALL CAXPY(M-j,alpha,A(j+1,k-1),1,Work(2),1)
               ENDIF
!
!           Find max(|WORK(2:n)|)
!
               i2 = ICAMAX(M-j,Work(2),1) + 1
               piv = Work(i2)
!
!           Apply hermitian pivot
!
               IF ( (i2/=2) .AND. (piv/=0) ) THEN
!
!              Swap WORK(I1) and WORK(I2)
!
                  i1 = 2
                  Work(i2) = Work(i1)
                  Work(i1) = piv
!
!              Swap A(I1+1:N, I1) with A(I2, I1+1:N)
!
                  i1 = i1 + j - 1
                  i2 = i2 + j - 1
                  CALL CSWAP(i2-i1-1,A(i1+1,J1+i1-1),1,A(i2,J1+i1),Lda)
                  CALL CLACGV(i2-i1,A(i1+1,J1+i1-1),1)
                  CALL CLACGV(i2-i1-1,A(i2,J1+i1),Lda)
!
!              Swap A(I2+1:N, I1) with A(I2+1:N, I2)
!
                  IF ( i2<M )                                           &
     &                 CALL CSWAP(M-i2,A(i2+1,J1+i1-1),1,A(i2+1,J1+i2-1)&
     &                 ,1)
!
!              Swap A(I1, I1) with A(I2, I2)
!
                  piv = A(i1,J1+i1-1)
                  A(i1,J1+i1-1) = A(i2,J1+i2-1)
                  A(i2,J1+i2-1) = piv
!
!              Swap H(I1, I1:J1) with H(I2, I2:J1)
!
                  CALL CSWAP(i1-1,H(i1,1),Ldh,H(i2,1),Ldh)
                  Ipiv(i1) = i2
!
!
!                 Swap L(1:I1-1, I1) with L(1:I1-1, I2),
!                  skipping the first column
!
                  IF ( i1>(k1-1) )                                      &
     &                 CALL CSWAP(i1-k1+1,A(i1,1),Lda,A(i2,1),Lda)
               ELSE
                  Ipiv(j+1) = j + 1
               ENDIF
!
!           Set A(J+1, J) = T(J+1, J)
!
               A(j+1,k) = Work(2)
!
!
!              Copy A(J+1:N, J+1) into H(J+1:N, J),
!
               IF ( j<Nb ) CALL CCOPY(M-j,A(j+1,k+1),1,H(j+1,j+1),1)
!
!           Compute L(J+2, J+1) = WORK( 3:N ) / T(J, J+1),
!            where A(J, J+1) = T(J, J+1) and A(J+2:N, J) = L(J+2:N, J+1)
!
               IF ( j<(M-1) ) THEN
                  IF ( A(j+1,k)/=ZERO ) THEN
                     alpha = ONE/A(j+1,k)
                     CALL CCOPY(M-j-1,Work(3),1,A(j+2,k),1)
                     CALL CSCAL(M-j-1,alpha,A(j+2,k),1)
                  ELSE
                     CALL CLASET('Full',M-j-1,1,ZERO,ZERO,A(j+2,k),Lda)
                  ENDIF
               ENDIF
            ENDIF
            j = j + 1
         ENDDO
      ENDIF
!
!     End of CLAHEF_AA
!
      END SUBROUTINE CLAHEF_AA
