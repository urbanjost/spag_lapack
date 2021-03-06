!*==spbtrf.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SPBTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SPBTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spbtrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spbtrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spbtrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SPBTRF( UPLO, N, KD, AB, LDAB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KD, LDAB, N
!       ..
!       .. Array Arguments ..
!       REAL               AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SPBTRF computes the Cholesky factorization of a real symmetric
!> positive definite band matrix A.
!>
!> The factorization has the form
!>    A = U**T * U,  if UPLO = 'U', or
!>    A = L  * L**T,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular.
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
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is REAL array, dimension (LDAB,N)
!>          On entry, the upper or lower triangle of the symmetric band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>
!>          On exit, if INFO = 0, the triangular factor U or L from the
!>          Cholesky factorization A = U**T*U or A = L*L**T of the band
!>          matrix A, in the same storage format as A.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the leading minor of order i is not
!>                positive definite, and the factorization could not be
!>                completed.
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
!> \ingroup realOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The band storage scheme is illustrated by the following example, when
!>  N = 6, KD = 2, and UPLO = 'U':
!>
!>  On entry:                       On exit:
!>
!>      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!>
!>  Similarly, if UPLO = 'L' the format of A is as follows:
!>
!>  On entry:                       On exit:
!>
!>     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!>     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!>     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!>
!>  Array elements marked * are not used by the routine.
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
!
!  =====================================================================
      SUBROUTINE SPBTRF(Uplo,N,Kd,Ab,Ldab,Info)
      IMPLICIT NONE
!*--SPBTRF146
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Kd , Ldab , N
!     ..
!     .. Array Arguments ..
      REAL Ab(Ldab,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      INTEGER NBMAX , LDWORK
      PARAMETER (NBMAX=32,LDWORK=NBMAX+1)
!     ..
!     .. Local Scalars ..
      INTEGER i , i2 , i3 , ib , ii , j , jj , nb
!     ..
!     .. Local Arrays ..
      REAL work(LDWORK,NBMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      EXTERNAL LSAME , ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEMM , SPBTF2 , SPOTF2 , SSYRK , STRSM , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      IF ( (.NOT.LSAME(Uplo,'U')) .AND. (.NOT.LSAME(Uplo,'L')) ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Kd<0 ) THEN
         Info = -3
      ELSEIF ( Ldab<Kd+1 ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SPBTRF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Determine the block size for this environment
!
      nb = ILAENV(1,'SPBTRF',Uplo,N,Kd,-1,-1)
!
!     The block size must not exceed the semi-bandwidth KD, and must not
!     exceed the limit set by the size of the local array WORK.
!
      nb = MIN(nb,NBMAX)
!
      IF ( nb<=1 .OR. nb>Kd ) THEN
!
!        Use unblocked code
!
         CALL SPBTF2(Uplo,N,Kd,Ab,Ldab,Info)
!
!        Use blocked code
!
      ELSEIF ( LSAME(Uplo,'U') ) THEN
!
!           Compute the Cholesky factorization of a symmetric band
!           matrix, given the upper triangle of the matrix in band
!           storage.
!
!           Zero the upper triangle of the work array.
!
         DO j = 1 , nb
            DO i = 1 , j - 1
               work(i,j) = ZERO
            ENDDO
         ENDDO
!
!           Process the band matrix one diagonal block at a time.
!
         DO i = 1 , N , nb
            ib = MIN(nb,N-i+1)
!
!              Factorize the diagonal block
!
            CALL SPOTF2(Uplo,ib,Ab(Kd+1,i),Ldab-1,ii)
            IF ( ii/=0 ) THEN
               Info = i + ii - 1
               GOTO 99999
            ENDIF
            IF ( i+ib<=N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11   A12   A13
!                          A22   A23
!                                A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A12, A22 and
!                 A23 are empty if IB = KD. The upper triangle of A13
!                 lies outside the band.
!
               i2 = MIN(Kd-ib,N-i-ib+1)
               i3 = MIN(ib,N-i-Kd+1)
!
               IF ( i2>0 ) THEN
!
!                    Update A12
!
                  CALL STRSM('Left','Upper','Transpose','Non-unit',ib,  &
     &                       i2,ONE,Ab(Kd+1,i),Ldab-1,Ab(Kd+1-ib,i+ib), &
     &                       Ldab-1)
!
!                    Update A22
!
                  CALL SSYRK('Upper','Transpose',i2,ib,-ONE,            &
     &                       Ab(Kd+1-ib,i+ib),Ldab-1,ONE,Ab(Kd+1,i+ib), &
     &                       Ldab-1)
               ENDIF
!
               IF ( i3>0 ) THEN
!
!                    Copy the lower triangle of A13 into the work array.
!
                  DO jj = 1 , i3
                     DO ii = jj , ib
                        work(ii,jj) = Ab(ii-jj+1,jj+i+Kd-1)
                     ENDDO
                  ENDDO
!
!                    Update A13 (in the work array).
!
                  CALL STRSM('Left','Upper','Transpose','Non-unit',ib,  &
     &                       i3,ONE,Ab(Kd+1,i),Ldab-1,work,LDWORK)
!
!                    Update A23
!
                  IF ( i2>0 ) CALL SGEMM('Transpose','No Transpose',i2, &
     &                 i3,ib,-ONE,Ab(Kd+1-ib,i+ib),Ldab-1,work,LDWORK,  &
     &                 ONE,Ab(1+ib,i+Kd),Ldab-1)
!
!                    Update A33
!
                  CALL SSYRK('Upper','Transpose',i3,ib,-ONE,work,LDWORK,&
     &                       ONE,Ab(Kd+1,i+Kd),Ldab-1)
!
!                    Copy the lower triangle of A13 back into place.
!
                  DO jj = 1 , i3
                     DO ii = jj , ib
                        Ab(ii-jj+1,jj+i+Kd-1) = work(ii,jj)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ELSE
!
!           Compute the Cholesky factorization of a symmetric band
!           matrix, given the lower triangle of the matrix in band
!           storage.
!
!           Zero the lower triangle of the work array.
!
         DO j = 1 , nb
            DO i = j + 1 , nb
               work(i,j) = ZERO
            ENDDO
         ENDDO
!
!           Process the band matrix one diagonal block at a time.
!
         DO i = 1 , N , nb
            ib = MIN(nb,N-i+1)
!
!              Factorize the diagonal block
!
            CALL SPOTF2(Uplo,ib,Ab(1,i),Ldab-1,ii)
            IF ( ii/=0 ) THEN
               Info = i + ii - 1
               GOTO 99999
            ENDIF
            IF ( i+ib<=N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11
!                    A21   A22
!                    A31   A32   A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A21, A22 and
!                 A32 are empty if IB = KD. The lower triangle of A31
!                 lies outside the band.
!
               i2 = MIN(Kd-ib,N-i-ib+1)
               i3 = MIN(ib,N-i-Kd+1)
!
               IF ( i2>0 ) THEN
!
!                    Update A21
!
                  CALL STRSM('Right','Lower','Transpose','Non-unit',i2, &
     &                       ib,ONE,Ab(1,i),Ldab-1,Ab(1+ib,i),Ldab-1)
!
!                    Update A22
!
                  CALL SSYRK('Lower','No Transpose',i2,ib,-ONE,         &
     &                       Ab(1+ib,i),Ldab-1,ONE,Ab(1,i+ib),Ldab-1)
               ENDIF
!
               IF ( i3>0 ) THEN
!
!                    Copy the upper triangle of A31 into the work array.
!
                  DO jj = 1 , ib
                     DO ii = 1 , MIN(jj,i3)
                        work(ii,jj) = Ab(Kd+1-jj+ii,jj+i-1)
                     ENDDO
                  ENDDO
!
!                    Update A31 (in the work array).
!
                  CALL STRSM('Right','Lower','Transpose','Non-unit',i3, &
     &                       ib,ONE,Ab(1,i),Ldab-1,work,LDWORK)
!
!                    Update A32
!
                  IF ( i2>0 ) CALL SGEMM('No transpose','Transpose',i3, &
     &                 i2,ib,-ONE,work,LDWORK,Ab(1+ib,i),Ldab-1,ONE,    &
     &                 Ab(1+Kd-ib,i+ib),Ldab-1)
!
!                    Update A33
!
                  CALL SSYRK('Lower','No Transpose',i3,ib,-ONE,work,    &
     &                       LDWORK,ONE,Ab(1,i+Kd),Ldab-1)
!
!                    Copy the upper triangle of A31 back into place.
!
                  DO jj = 1 , ib
                     DO ii = 1 , MIN(jj,i3)
                        Ab(Kd+1-jj+ii,jj+i-1) = work(ii,jj)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      RETURN
!
!
!     End of SPBTRF
!
99999 END SUBROUTINE SPBTRF
