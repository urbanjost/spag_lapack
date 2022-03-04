!*==zpbtrf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZPBTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZPBTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpbtrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpbtrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpbtrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPBTRF( UPLO, N, KD, AB, LDAB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KD, LDAB, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPBTRF computes the Cholesky factorization of a complex Hermitian
!> positive definite band matrix A.
!>
!> The factorization has the form
!>    A = U**H * U,  if UPLO = 'U', or
!>    A = L  * L**H,  if UPLO = 'L',
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
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          On entry, the upper or lower triangle of the Hermitian band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>
!>          On exit, if INFO = 0, the triangular factor U or L from the
!>          Cholesky factorization A = U**H*U or A = L*L**H of the band
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
!> \ingroup complex16OTHERcomputational
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
      SUBROUTINE ZPBTRF(Uplo,N,Kd,Ab,Ldab,Info)
      USE F77KINDS                        
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      USE S_ZGEMM
      USE S_ZHERK
      USE S_ZPBTF2
      USE S_ZPOTF2
      USE S_ZTRSM
      IMPLICIT NONE
!*--ZPBTRF155
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
      INTEGER , PARAMETER  ::  NBMAX = 32 , LDWORK = NBMAX + 1
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , i2 , i3 , ib , ii , j , jj , nb
      COMPLEX(CX16KIND) , DIMENSION(LDWORK,NBMAX) :: work
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
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
         CALL XERBLA('ZPBTRF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Determine the block size for this environment
!
      nb = ILAENV(1,'ZPBTRF',Uplo,N,Kd,-1,-1)
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
         CALL ZPBTF2(Uplo,N,Kd,Ab,Ldab,Info)
!
!        Use blocked code
!
      ELSEIF ( LSAME(Uplo,'U') ) THEN
!
!           Compute the Cholesky factorization of a Hermitian band
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
            CALL ZPOTF2(Uplo,ib,Ab(Kd+1,i),Ldab-1,ii)
            IF ( ii/=0 ) THEN
               Info = i + ii - 1
               EXIT
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
                  CALL ZTRSM('Left','Upper','Conjugate transpose',      &
     &                       'Non-unit',ib,i2,CONE,Ab(Kd+1,i),Ldab-1,   &
     &                       Ab(Kd+1-ib,i+ib),Ldab-1)
!
!                    Update A22
!
                  CALL ZHERK('Upper','Conjugate transpose',i2,ib,-ONE,  &
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
                  CALL ZTRSM('Left','Upper','Conjugate transpose',      &
     &                       'Non-unit',ib,i3,CONE,Ab(Kd+1,i),Ldab-1,   &
     &                       work,LDWORK)
!
!                    Update A23
!
                  IF ( i2>0 )                                           &
     &                  CALL ZGEMM('Conjugate transpose','No transpose',&
     &                 i2,i3,ib,-CONE,Ab(Kd+1-ib,i+ib),Ldab-1,work,     &
     &                 LDWORK,CONE,Ab(1+ib,i+Kd),Ldab-1)
!
!                    Update A33
!
                  CALL ZHERK('Upper','Conjugate transpose',i3,ib,-ONE,  &
     &                       work,LDWORK,ONE,Ab(Kd+1,i+Kd),Ldab-1)
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
!           Compute the Cholesky factorization of a Hermitian band
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
            CALL ZPOTF2(Uplo,ib,Ab(1,i),Ldab-1,ii)
            IF ( ii/=0 ) THEN
               Info = i + ii - 1
               EXIT
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
                  CALL ZTRSM('Right','Lower','Conjugate transpose',     &
     &                       'Non-unit',i2,ib,CONE,Ab(1,i),Ldab-1,      &
     &                       Ab(1+ib,i),Ldab-1)
!
!                    Update A22
!
                  CALL ZHERK('Lower','No transpose',i2,ib,-ONE,         &
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
                  CALL ZTRSM('Right','Lower','Conjugate transpose',     &
     &                       'Non-unit',i3,ib,CONE,Ab(1,i),Ldab-1,work, &
     &                       LDWORK)
!
!                    Update A32
!
                  IF ( i2>0 )                                           &
     &                  CALL ZGEMM('No transpose','Conjugate transpose',&
     &                 i3,i2,ib,-CONE,work,LDWORK,Ab(1+ib,i),Ldab-1,    &
     &                 CONE,Ab(1+Kd-ib,i+ib),Ldab-1)
!
!                    Update A33
!
                  CALL ZHERK('Lower','No transpose',i3,ib,-ONE,work,    &
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
!
!
!     End of ZPBTRF
!
      END SUBROUTINE ZPBTRF
