!*==slasyf_aa.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLASYF_AA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASYF_AA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasyf_aa.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasyf_aa.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasyf_aa.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASYF_AA( UPLO, J1, M, NB, A, LDA, IPIV,
!                             H, LDH, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            J1, M, NB, LDA, LDH
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * ), H( LDH, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLATRF_AA factorizes a panel of a real symmetric matrix A using
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
!>          when called by SSYTRF_AA, for the first panel, J1 is 1,
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
!>          A is REAL array, dimension (LDA,M) for
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
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (M)
!>          Details of the row and column interchanges,
!>          the row and column k were interchanged with the row and
!>          column IPIV(k).
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is REAL workspace, dimension (LDH,NB).
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
!>          WORK is REAL workspace, dimension (M).
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
!> \ingroup realSYcomputational
!
!  =====================================================================
      SUBROUTINE SLASYF_AA(Uplo,J1,M,Nb,A,Lda,Ipiv,H,Ldh,Work)
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
      USE S_ILAENV
      USE S_ISAMAX
      USE S_LSAME
      USE S_SAXPY
      USE S_SCOPY
      USE S_SGEMV
      USE S_SLASET
      USE S_SSCAL
      USE S_SSWAP
      USE S_XERBLA
      IMPLICIT NONE
!*--SLASYF_AA163
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: J1
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: Nb
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      REAL :: alpha , piv
      INTEGER :: i1 , i2 , j , k , k1 , mj
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!     .. Parameters ..
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
!         when being called from SSYTRF_AA,
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
!        H(J:M, J) := A(J, J:M) - H(J:M, 1:(J-1)) * L(J1:(J-1), J),
!         where H(J:M, J) has been initialized to be A(J, J:M)
!
!
!        K is the column to be factorized
!         > for the first block column, K is J, skipping the first two
!           columns
!         > for the rest of the columns, K is J+1, skipping only the
!           first column
!
            IF ( k>2 ) CALL SGEMV('No transpose',mj,j-k1,-ONE,H(j,k1),  &
     &                            Ldh,A(1,j),1,ONE,H(j,j),1)
!
!        Copy H(i:M, i) into WORK
!
            CALL SCOPY(mj,H(j,j),1,Work(1),1)
!
            IF ( j>k1 ) THEN
!
!           Compute WORK := WORK - L(J-1, J:M) * T(J-1,J),
!            where A(J-1, J) stores T(J-1, J) and A(J-2, J:M) stores U(J-1, J:M)
!
               alpha = -A(k-1,j)
               CALL SAXPY(mj,alpha,A(k-2,j),Lda,Work(1),1)
            ENDIF
!
!        Set A(J, J) = T(J, J)
!
            A(k,j) = Work(1)
!
            IF ( j<M ) THEN
!
!           Compute WORK(2:M) = T(J, J) L(J, (J+1):M)
!            where A(J, J) stores T(J, J) and A(J-1, (J+1):M) stores U(J, (J+1):M)
!
               IF ( k>1 ) THEN
                  alpha = -A(k,j)
                  CALL SAXPY(M-j,alpha,A(k-1,j+1),Lda,Work(2),1)
               ENDIF
!
!           Find max(|WORK(2:M)|)
!
               i2 = ISAMAX(M-j,Work(2),1) + 1
               piv = Work(i2)
!
!           Apply symmetric pivot
!
               IF ( (i2/=2) .AND. (piv/=0) ) THEN
!
!              Swap WORK(I1) and WORK(I2)
!
                  i1 = 2
                  Work(i2) = Work(i1)
                  Work(i1) = piv
!
!              Swap A(I1, I1+1:M) with A(I1+1:M, I2)
!
                  i1 = i1 + j - 1
                  i2 = i2 + j - 1
                  CALL SSWAP(i2-i1-1,A(J1+i1-1,i1+1),Lda,A(J1+i1,i2),1)
!
!              Swap A(I1, I2+1:M) with A(I2, I2+1:M)
!
                  IF ( i2<M ) CALL SSWAP(M-i2,A(J1+i1-1,i2+1),Lda,      &
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
                  CALL SSWAP(i1-1,H(i1,1),Ldh,H(i2,1),Ldh)
                  Ipiv(i1) = i2
!
!
!                 Swap L(1:I1-1, I1) with L(1:I1-1, I2),
!                  skipping the first column
!
                  IF ( i1>(k1-1) ) CALL SSWAP(i1-k1+1,A(1,i1),1,A(1,i2),&
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
!              Copy A(J+1:M, J+1) into H(J:M, J),
!
               IF ( j<Nb ) CALL SCOPY(M-j,A(k+1,j+1),Lda,H(j+1,j+1),1)
!
!           Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1),
!            where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1)
!
               IF ( j<(M-1) ) THEN
                  IF ( A(k,j+1)/=ZERO ) THEN
                     alpha = ONE/A(k,j+1)
                     CALL SCOPY(M-j-1,Work(3),1,A(k,j+2),Lda)
                     CALL SSCAL(M-j-1,alpha,A(k,j+2),Lda)
                  ELSE
                     CALL SLASET('Full',1,M-j-1,ZERO,ZERO,A(k,j+2),Lda)
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
!         when being called from SSYTRF_AA,
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
!        H(J:M, J) := A(J:M, J) - H(J:M, 1:(J-1)) * L(J, J1:(J-1))^T,
!         where H(J:M, J) has been initialized to be A(J:M, J)
!
!
!        K is the column to be factorized
!         > for the first block column, K is J, skipping the first two
!           columns
!         > for the rest of the columns, K is J+1, skipping only the
!           first column
!
            IF ( k>2 ) CALL SGEMV('No transpose',mj,j-k1,-ONE,H(j,k1),  &
     &                            Ldh,A(j,1),Lda,ONE,H(j,j),1)
!
!        Copy H(J:M, J) into WORK
!
            CALL SCOPY(mj,H(j,j),1,Work(1),1)
!
            IF ( j>k1 ) THEN
!
!           Compute WORK := WORK - L(J:M, J-1) * T(J-1,J),
!            where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1)
!
               alpha = -A(j,k-1)
               CALL SAXPY(mj,alpha,A(j,k-2),1,Work(1),1)
            ENDIF
!
!        Set A(J, J) = T(J, J)
!
            A(j,k) = Work(1)
!
            IF ( j<M ) THEN
!
!           Compute WORK(2:M) = T(J, J) L((J+1):M, J)
!            where A(J, J) = T(J, J) and A((J+1):M, J-1) = L((J+1):M, J)
!
               IF ( k>1 ) THEN
                  alpha = -A(j,k)
                  CALL SAXPY(M-j,alpha,A(j+1,k-1),1,Work(2),1)
               ENDIF
!
!           Find max(|WORK(2:M)|)
!
               i2 = ISAMAX(M-j,Work(2),1) + 1
               piv = Work(i2)
!
!           Apply symmetric pivot
!
               IF ( (i2/=2) .AND. (piv/=0) ) THEN
!
!              Swap WORK(I1) and WORK(I2)
!
                  i1 = 2
                  Work(i2) = Work(i1)
                  Work(i1) = piv
!
!              Swap A(I1+1:M, I1) with A(I2, I1+1:M)
!
                  i1 = i1 + j - 1
                  i2 = i2 + j - 1
                  CALL SSWAP(i2-i1-1,A(i1+1,J1+i1-1),1,A(i2,J1+i1),Lda)
!
!              Swap A(I2+1:M, I1) with A(I2+1:M, I2)
!
                  IF ( i2<M )                                           &
     &                 CALL SSWAP(M-i2,A(i2+1,J1+i1-1),1,A(i2+1,J1+i2-1)&
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
                  CALL SSWAP(i1-1,H(i1,1),Ldh,H(i2,1),Ldh)
                  Ipiv(i1) = i2
!
!
!                 Swap L(1:I1-1, I1) with L(1:I1-1, I2),
!                  skipping the first column
!
                  IF ( i1>(k1-1) )                                      &
     &                 CALL SSWAP(i1-k1+1,A(i1,1),Lda,A(i2,1),Lda)
               ELSE
                  Ipiv(j+1) = j + 1
               ENDIF
!
!           Set A(J+1, J) = T(J+1, J)
!
               A(j+1,k) = Work(2)
!
!
!              Copy A(J+1:M, J+1) into H(J+1:M, J),
!
               IF ( j<Nb ) CALL SCOPY(M-j,A(j+1,k+1),1,H(j+1,j+1),1)
!
!           Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1),
!            where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1)
!
               IF ( j<(M-1) ) THEN
                  IF ( A(j+1,k)/=ZERO ) THEN
                     alpha = ONE/A(j+1,k)
                     CALL SCOPY(M-j-1,Work(3),1,A(j+2,k),1)
                     CALL SSCAL(M-j-1,alpha,A(j+2,k),1)
                  ELSE
                     CALL SLASET('Full',M-j-1,1,ZERO,ZERO,A(j+2,k),Lda)
                  ENDIF
               ENDIF
            ENDIF
            j = j + 1
         ENDDO
      ENDIF
!
!     End of SLASYF_AA
!
      END SUBROUTINE SLASYF_AA
