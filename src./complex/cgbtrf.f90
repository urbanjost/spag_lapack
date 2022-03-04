!*==cgbtrf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGBTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGBTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbtrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbtrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbtrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, KL, KU, LDAB, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            AB( LDAB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGBTRF computes an LU factorization of a complex m-by-n band matrix A
!> using partial pivoting with row interchanges.
!>
!> This is the blocked version of the algorithm, calling Level 3 BLAS.
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
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB,N)
!>          On entry, the matrix A in band storage, in rows KL+1 to
!>          2*KL+KU+1; rows 1 to KL of the array need not be set.
!>          The j-th column of A is stored in the j-th column of the
!>          array AB as follows:
!>          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!>
!>          On exit, details of the factorization: U is stored as an
!>          upper triangular band matrix with KL+KU superdiagonals in
!>          rows 1 to KL+KU+1, and the multipliers used during the
!>          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!>          See below for further details.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!>               has been completed, but the factor U is exactly
!>               singular, and division by zero will occur if it is used
!>               to solve a system of equations.
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
!> \ingroup complexGBcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The band storage scheme is illustrated by the following example, when
!>  M = N = 6, KL = 2, KU = 1:
!>
!>  On entry:                       On exit:
!>
!>      *    *    *    +    +    +       *    *    *   u14  u25  u36
!>      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!>     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!>     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!>
!>  Array elements marked * are not used by the routine; elements marked
!>  + need not be set on entry, but are required by the routine to store
!>  elements of U because of fill-in resulting from the row interchanges.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CGBTRF(M,N,Kl,Ku,Ab,Ldab,Ipiv,Info)
      USE S_CCOPY
      USE S_CGBTF2
      USE S_CGEMM
      USE S_CGERU
      USE S_CLASWP
      USE S_CSCAL
      USE S_CSWAP
      USE S_CTRSM
      USE S_ICAMAX
      USE S_ILAENV
      USE S_XERBLA
      IMPLICIT NONE
!*--CGBTRF159
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDWORK = NBMAX + 1
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , i2 , i3 , ii , ip , j , j2 , j3 , jb , jj , jm ,   &
     &           jp , ju , k2 , km , kv , nb , nw
      COMPLEX :: temp
      COMPLEX , DIMENSION(LDWORK,NBMAX) :: work13 , work31
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
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
      kv = Ku + Kl
!
!     Test the input parameters.
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Kl<0 ) THEN
         Info = -3
      ELSEIF ( Ku<0 ) THEN
         Info = -4
      ELSEIF ( Ldab<Kl+kv+1 ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGBTRF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) RETURN
!
!     Determine the block size for this environment
!
      nb = ILAENV(1,'CGBTRF',' ',M,N,Kl,Ku)
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
      nb = MIN(nb,NBMAX)
!
      IF ( nb<=1 .OR. nb>Kl ) THEN
!
!        Use unblocked code
!
         CALL CGBTF2(M,N,Kl,Ku,Ab,Ldab,Ipiv,Info)
      ELSE
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
         DO j = 1 , nb
            DO i = 1 , j - 1
               work13(i,j) = ZERO
            ENDDO
         ENDDO
!
!        Zero the subdiagonal elements of the work array WORK31
!
         DO j = 1 , nb
            DO i = j + 1 , nb
               work31(i,j) = ZERO
            ENDDO
         ENDDO
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
         DO j = Ku + 2 , MIN(kv,N)
            DO i = kv - j + 2 , Kl
               Ab(i,j) = ZERO
            ENDDO
         ENDDO
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
         ju = 1
!
         DO j = 1 , MIN(M,N) , nb
            jb = MIN(nb,MIN(M,N)-j+1)
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
            i2 = MIN(Kl-jb,M-j-jb+1)
            i3 = MIN(jb,M-j-Kl+1)
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
            DO jj = j , j + jb - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
               IF ( jj+kv<=N ) THEN
                  DO i = 1 , Kl
                     Ab(i,jj+kv) = ZERO
                  ENDDO
               ENDIF
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
               km = MIN(Kl,M-jj)
               jp = ICAMAX(km+1,Ab(kv+1,jj),1)
               Ipiv(jj) = jp + jj - j
               IF ( Ab(kv+jp,jj)/=ZERO ) THEN
                  ju = MAX(ju,MIN(jj+Ku+jp-1,N))
                  IF ( jp/=1 ) THEN
!
!                    Apply interchange to columns J to J+JB-1
!
                     IF ( jp+jj-1<j+Kl ) THEN
!
                        CALL CSWAP(jb,Ab(kv+1+jj-j,j),Ldab-1,           &
     &                             Ab(kv+jp+jj-j,j),Ldab-1)
                     ELSE
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                        CALL CSWAP(jj-j,Ab(kv+1+jj-j,j),Ldab-1,         &
     &                             work31(jp+jj-j-Kl,1),LDWORK)
                        CALL CSWAP(j+jb-jj,Ab(kv+1,jj),Ldab-1,          &
     &                             Ab(kv+jp,jj),Ldab-1)
                     ENDIF
                  ENDIF
!
!                 Compute multipliers
!
                  CALL CSCAL(km,ONE/Ab(kv+1,jj),Ab(kv+2,jj),1)
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
                  jm = MIN(ju,j+jb-1)
                  IF ( jm>jj ) CALL CGERU(km,jm-jj,-ONE,Ab(kv+2,jj),1,  &
     &                 Ab(kv,jj+1),Ldab-1,Ab(kv+1,jj+1),Ldab-1)
               ELSE
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
                  IF ( Info==0 ) Info = jj
               ENDIF
!
!              Copy current column of A31 into the work array WORK31
!
               nw = MIN(jj-j+1,i3)
               IF ( nw>0 ) CALL CCOPY(nw,Ab(kv+Kl+1-jj+j,jj),1,         &
     &                                work31(1,jj-j+1),1)
            ENDDO
            IF ( j+jb<=N ) THEN
!
!              Apply the row interchanges to the other blocks.
!
               j2 = MIN(ju-j+1,kv) - jb
               j3 = MAX(0,ju-j-kv+1)
!
!              Use CLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
               CALL CLASWP(j2,Ab(kv+1-jb,j+jb),Ldab-1,1,jb,Ipiv(j),1)
!
!              Adjust the pivot indices.
!
               DO i = j , j + jb - 1
                  Ipiv(i) = Ipiv(i) + j - 1
               ENDDO
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
               k2 = j - 1 + jb + j2
               DO i = 1 , j3
                  jj = k2 + i
                  DO ii = j + i - 1 , j + jb - 1
                     ip = Ipiv(ii)
                     IF ( ip/=ii ) THEN
                        temp = Ab(kv+1+ii-jj,jj)
                        Ab(kv+1+ii-jj,jj) = Ab(kv+1+ip-jj,jj)
                        Ab(kv+1+ip-jj,jj) = temp
                     ENDIF
                  ENDDO
               ENDDO
!
!              Update the relevant part of the trailing submatrix
!
               IF ( j2>0 ) THEN
!
!                 Update A12
!
                  CALL CTRSM('Left','Lower','No transpose','Unit',jb,j2,&
     &                       ONE,Ab(kv+1,j),Ldab-1,Ab(kv+1-jb,j+jb),    &
     &                       Ldab-1)
!
!
!                    Update A22
!
                  IF ( i2>0 ) CALL CGEMM('No transpose','No transpose', &
     &                 i2,j2,jb,-ONE,Ab(kv+1+jb,j),Ldab-1,              &
     &                 Ab(kv+1-jb,j+jb),Ldab-1,ONE,Ab(kv+1,j+jb),Ldab-1)
!
!
!                    Update A32
!
                  IF ( i3>0 ) CALL CGEMM('No transpose','No transpose', &
     &                 i3,j2,jb,-ONE,work31,LDWORK,Ab(kv+1-jb,j+jb),    &
     &                 Ldab-1,ONE,Ab(kv+Kl+1-jb,j+jb),Ldab-1)
               ENDIF
!
               IF ( j3>0 ) THEN
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
                  DO jj = 1 , j3
                     DO ii = jj , jb
                        work13(ii,jj) = Ab(ii-jj+1,jj+j+kv-1)
                     ENDDO
                  ENDDO
!
!                 Update A13 in the work array
!
                  CALL CTRSM('Left','Lower','No transpose','Unit',jb,j3,&
     &                       ONE,Ab(kv+1,j),Ldab-1,work13,LDWORK)
!
!
!                    Update A23
!
                  IF ( i2>0 ) CALL CGEMM('No transpose','No transpose', &
     &                 i2,j3,jb,-ONE,Ab(kv+1+jb,j),Ldab-1,work13,LDWORK,&
     &                 ONE,Ab(1+jb,j+kv),Ldab-1)
!
!
!                    Update A33
!
                  IF ( i3>0 ) CALL CGEMM('No transpose','No transpose', &
     &                 i3,j3,jb,-ONE,work31,LDWORK,work13,LDWORK,ONE,   &
     &                 Ab(1+Kl,j+kv),Ldab-1)
!
!                 Copy the lower triangle of A13 back into place
!
                  DO jj = 1 , j3
                     DO ii = jj , jb
                        Ab(ii-jj+1,jj+j+kv-1) = work13(ii,jj)
                     ENDDO
                  ENDDO
               ENDIF
            ELSE
!
!              Adjust the pivot indices.
!
               DO i = j , j + jb - 1
                  Ipiv(i) = Ipiv(i) + j - 1
               ENDDO
            ENDIF
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
            DO jj = j + jb - 1 , j , -1
               jp = Ipiv(jj) - jj + 1
               IF ( jp/=1 ) THEN
!
!                 Apply interchange to columns J to JJ-1
!
                  IF ( jp+jj-1<j+Kl ) THEN
!
!                    The interchange does not affect A31
!
                     CALL CSWAP(jj-j,Ab(kv+1+jj-j,j),Ldab-1,            &
     &                          Ab(kv+jp+jj-j,j),Ldab-1)
                  ELSE
!
!                    The interchange does affect A31
!
                     CALL CSWAP(jj-j,Ab(kv+1+jj-j,j),Ldab-1,            &
     &                          work31(jp+jj-j-Kl,1),LDWORK)
                  ENDIF
               ENDIF
!
!              Copy the current column of A31 back into place
!
               nw = MIN(i3,jj-j+1)
               IF ( nw>0 ) CALL CCOPY(nw,work31(1,jj-j+1),1,Ab(kv+Kl+1- &
     &                                jj+j,jj),1)
            ENDDO
         ENDDO
      ENDIF
!
!
!     End of CGBTRF
!
      END SUBROUTINE CGBTRF
