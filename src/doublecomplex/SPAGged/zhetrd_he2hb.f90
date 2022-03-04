!*==zhetrd_he2hb.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZHETRD_HE2HB
!
!  @precisions fortran z -> s d c
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETRD_HE2HB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETRD_HE2HB( UPLO, N, KD, A, LDA, AB, LDAB, TAU,
!                              WORK, LWORK, INFO )
!
!       IMPLICIT NONE
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDAB, LWORK, N, KD
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), AB( LDAB, * ),
!                          TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHETRD_HE2HB reduces a complex Hermitian matrix A to complex Hermitian
!> band-diagonal form AB by a unitary similarity transformation:
!> Q**H * A * Q = AB.
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
!>          The number of superdiagonals of the reduced matrix if UPLO = 'U',
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!>          The reduced matrix is stored in the array AB.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!>          of A are overwritten by the corresponding elements of the
!>          tridiagonal matrix T, and the elements above the first
!>          superdiagonal, with the array TAU, represent the unitary
!>          matrix Q as a product of elementary reflectors; if UPLO
!>          = 'L', the diagonal and first subdiagonal of A are over-
!>          written by the corresponding elements of the tridiagonal
!>          matrix T, and the elements below the first subdiagonal, with
!>          the array TAU, represent the unitary matrix Q as a product
!>          of elementary reflectors. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] AB
!> \verbatim
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          On exit, the upper or lower triangle of the Hermitian band
!>          matrix A, stored in the first KD+1 rows of the array.  The
!>          j-th column of A is stored in the j-th column of the array AB
!>          as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (N-KD)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!>          On exit, if INFO = 0, or if LWORK=-1,
!>          WORK(1) returns the size of LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK which should be calculated
!>          by a workspace query. LWORK = MAX(1, LWORK_QUERY)
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!>          LWORK_QUERY = N*KD + N*max(KD,FACTOPTNB) + 2*KD*KD
!>          where FACTOPTNB is the blocking used by the QR or LQ
!>          algorithm, usually FACTOPTNB=128 is a good choice otherwise
!>          putting LWORK=-1 will provide the size of WORK.
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
!> \date November 2017
!
!> \ingroup complex16HEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Implemented by Azzam Haidar.
!>
!>  All details are available on technical report, SC11, SC13 papers.
!>
!>  Azzam Haidar, Hatem Ltaief, and Jack Dongarra.
!>  Parallel reduction to condensed forms for symmetric eigenvalue problems
!>  using aggregated fine-grained and memory-aware kernels. In Proceedings
!>  of 2011 International Conference for High Performance Computing,
!>  Networking, Storage and Analysis (SC '11), New York, NY, USA,
!>  Article 8 , 11 pages.
!>  http://doi.acm.org/10.1145/2063384.2063394
!>
!>  A. Haidar, J. Kurzak, P. Luszczek, 2013.
!>  An improved parallel singular value algorithm and its implementation
!>  for multicore hardware, In Proceedings of 2013 International Conference
!>  for High Performance Computing, Networking, Storage and Analysis (SC '13).
!>  Denver, Colorado, USA, 2013.
!>  Article 90, 12 pages.
!>  http://doi.acm.org/10.1145/2503210.2503292
!>
!>  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra.
!>  A novel hybrid CPU-GPU generalized eigensolver for electronic structure
!>  calculations based on fine-grained memory aware tasks.
!>  International Journal of High Performance Computing Applications.
!>  Volume 28 Issue 2, Pages 196-209, May 2014.
!>  http://hpc.sagepub.com/content/28/2/196
!>
!> \endverbatim
!>
!> \verbatim
!>
!>  If UPLO = 'U', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(k)**H . . . H(2)**H H(1)**H, where k = n-kd.
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(1:i+kd-1) = 0 and v(i+kd) = 1; conjg(v(i+kd+1:n)) is stored on exit in
!>  A(i,i+kd+1:n), and tau in TAU(i).
!>
!>  If UPLO = 'L', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(1) H(2) . . . H(k), where k = n-kd.
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(kd+1:i) = 0 and v(i+kd+1) = 1; v(i+kd+2:n) is stored on exit in
!>  A(i+kd+2:n,i), and tau in TAU(i).
!>
!>  The contents of A on exit are illustrated by the following examples
!>  with n = 5:
!>
!>  if UPLO = 'U':                       if UPLO = 'L':
!>
!>    (  ab  ab/v1  v1      v1     v1    )              (  ab                            )
!>    (      ab     ab/v2   v2     v2    )              (  ab/v1  ab                     )
!>    (             ab      ab/v3  v3    )              (  v1     ab/v2  ab              )
!>    (                     ab     ab/v4 )              (  v1     v2     ab/v3  ab       )
!>    (                            ab    )              (  v1     v2     v3     ab/v4 ab )
!>
!>  where d and e denote diagonal and off-diagonal elements of T, and vi
!>  denotes an element of the vector defining H(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZHETRD_HE2HB(Uplo,N,Kd,A,Lda,Ab,Ldab,Tau,Work,Lwork,   &
     &                        Info)
!
      USE F77KINDS                        
      USE S_ILAENV2STAGE
      USE S_LSAME
      USE S_XERBLA
      USE S_ZCOPY
      USE S_ZGELQF
      USE S_ZGEMM
      USE S_ZGEQRF
      USE S_ZHEMM
      USE S_ZHER2K
      USE S_ZLARFT
      USE S_ZLASET
      IMPLICIT NONE
!*--ZHETRD_HE2HB260
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  RONE = 1.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0) ,       &
     &                 ONE = (1.0D+0,0.0D+0) , HALF = (0.5D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , iinfo , j , lds1 , lds2 , ldt , ldw , lk , ls1 ,   &
     &           ls2 , lt , lw , lwmin , pk , pn , s1pos , s2pos ,      &
     &           tpos , wpos
      LOGICAL :: lquery , upper
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
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
!     Determine the minimal workspace size required
!     and test the input parameters
!
      Info = 0
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1)
      lwmin = ILAENV2STAGE(4,'ZHETRD_HE2HB','',N,Kd,-1,-1)
 
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Kd<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldab<MAX(1,Kd+1) ) THEN
         Info = -7
      ELSEIF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
         Info = -10
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHETRD_HE2HB',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         Work(1) = lwmin
         RETURN
      ENDIF
!
!     Quick return if possible
!     Copy the upper/lower portion of A into AB
!
      IF ( N<=Kd+1 ) THEN
         IF ( upper ) THEN
            DO i = 1 , N
               lk = MIN(Kd+1,i)
               CALL ZCOPY(lk,A(i-lk+1,i),1,Ab(Kd+1-lk+1,i),1)
            ENDDO
         ELSE
            DO i = 1 , N
               lk = MIN(Kd+1,N-i+1)
               CALL ZCOPY(lk,A(i,i),1,Ab(1,i),1)
            ENDDO
         ENDIF
         Work(1) = 1
         RETURN
      ENDIF
!
!     Determine the pointer position for the workspace
!
      ldt = Kd
      lds1 = Kd
      lt = ldt*Kd
      lw = N*Kd
      ls1 = lds1*Kd
      ls2 = lwmin - lt - lw - ls1
!      LS2 = N*MAX(KD,FACTOPTNB)
      tpos = 1
      wpos = tpos + lt
      s1pos = wpos + lw
      s2pos = s1pos + ls1
      IF ( upper ) THEN
         ldw = Kd
         lds2 = Kd
      ELSE
         ldw = N
         lds2 = N
      ENDIF
!
!
!     Set the workspace of the triangular matrix T to zero once such a
!     way every time T is generated the upper/lower portion will be always zero
!
      CALL ZLASET("A",ldt,Kd,ZERO,ZERO,Work(tpos),ldt)
!
      IF ( upper ) THEN
         DO i = 1 , N - Kd , Kd
            pn = N - i - Kd + 1
            pk = MIN(N-i-Kd+1,Kd)
!
!            Compute the LQ factorization of the current block
!
            CALL ZGELQF(Kd,pn,A(i,i+Kd),Lda,Tau(i),Work(s2pos),ls2,     &
     &                  iinfo)
!
!            Copy the upper portion of A into AB
!
            DO j = i , i + pk - 1
               lk = MIN(Kd,N-j) + 1
               CALL ZCOPY(lk,A(j,j),Lda,Ab(Kd+1,j),Ldab-1)
            ENDDO
!
            CALL ZLASET('Lower',pk,pk,ZERO,ONE,A(i,i+Kd),Lda)
!
!            Form the matrix T
!
            CALL ZLARFT('Forward','Rowwise',pn,pk,A(i,i+Kd),Lda,Tau(i), &
     &                  Work(tpos),ldt)
!
!            Compute W:
!
            CALL ZGEMM('Conjugate','No transpose',pk,pn,pk,ONE,         &
     &                 Work(tpos),ldt,A(i,i+Kd),Lda,ZERO,Work(s2pos),   &
     &                 lds2)
!
            CALL ZHEMM('Right',Uplo,pk,pn,ONE,A(i+Kd,i+Kd),Lda,         &
     &                 Work(s2pos),lds2,ZERO,Work(wpos),ldw)
!
            CALL ZGEMM('No transpose','Conjugate',pk,pk,pn,ONE,         &
     &                 Work(wpos),ldw,Work(s2pos),lds2,ZERO,Work(s1pos),&
     &                 lds1)
!
            CALL ZGEMM('No transpose','No transpose',pk,pn,pk,-HALF,    &
     &                 Work(s1pos),lds1,A(i,i+Kd),Lda,ONE,Work(wpos),   &
     &                 ldw)
!
!
!            Update the unreduced submatrix A(i+kd:n,i+kd:n), using
!            an update of the form:  A := A - V'*W - W'*V
!
            CALL ZHER2K(Uplo,'Conjugate',pn,pk,-ONE,A(i,i+Kd),Lda,      &
     &                  Work(wpos),ldw,RONE,A(i+Kd,i+Kd),Lda)
         ENDDO
!
!        Copy the upper band to AB which is the band storage matrix
!
         DO j = N - Kd + 1 , N
            lk = MIN(Kd,N-j) + 1
            CALL ZCOPY(lk,A(j,j),Lda,Ab(Kd+1,j),Ldab-1)
         ENDDO
!
      ELSE
!
!         Reduce the lower triangle of A to lower band matrix
!
         DO i = 1 , N - Kd , Kd
            pn = N - i - Kd + 1
            pk = MIN(N-i-Kd+1,Kd)
!
!            Compute the QR factorization of the current block
!
            CALL ZGEQRF(pn,Kd,A(i+Kd,i),Lda,Tau(i),Work(s2pos),ls2,     &
     &                  iinfo)
!
!            Copy the upper portion of A into AB
!
            DO j = i , i + pk - 1
               lk = MIN(Kd,N-j) + 1
               CALL ZCOPY(lk,A(j,j),1,Ab(1,j),1)
            ENDDO
!
            CALL ZLASET('Upper',pk,pk,ZERO,ONE,A(i+Kd,i),Lda)
!
!            Form the matrix T
!
            CALL ZLARFT('Forward','Columnwise',pn,pk,A(i+Kd,i),Lda,     &
     &                  Tau(i),Work(tpos),ldt)
!
!            Compute W:
!
            CALL ZGEMM('No transpose','No transpose',pn,pk,pk,ONE,      &
     &                 A(i+Kd,i),Lda,Work(tpos),ldt,ZERO,Work(s2pos),   &
     &                 lds2)
!
            CALL ZHEMM('Left',Uplo,pn,pk,ONE,A(i+Kd,i+Kd),Lda,          &
     &                 Work(s2pos),lds2,ZERO,Work(wpos),ldw)
!
            CALL ZGEMM('Conjugate','No transpose',pk,pk,pn,ONE,         &
     &                 Work(s2pos),lds2,Work(wpos),ldw,ZERO,Work(s1pos),&
     &                 lds1)
!
            CALL ZGEMM('No transpose','No transpose',pn,pk,pk,-HALF,    &
     &                 A(i+Kd,i),Lda,Work(s1pos),lds1,ONE,Work(wpos),   &
     &                 ldw)
!
!
!            Update the unreduced submatrix A(i+kd:n,i+kd:n), using
!            an update of the form:  A := A - V*W' - W*V'
!
            CALL ZHER2K(Uplo,'No transpose',pn,pk,-ONE,A(i+Kd,i),Lda,   &
     &                  Work(wpos),ldw,RONE,A(i+Kd,i+Kd),Lda)
!            ==================================================================
!            RESTORE A FOR COMPARISON AND CHECKING TO BE REMOVED
!             DO 45 J = I, I+PK-1
!                LK = MIN( KD, N-J ) + 1
!                CALL ZCOPY( LK, AB( 1, J ), 1, A( J, J ), 1 )
!   45        CONTINUE
!            ==================================================================
         ENDDO
!
!        Copy the lower band to AB which is the band storage matrix
!
         DO j = N - Kd + 1 , N
            lk = MIN(Kd,N-j) + 1
            CALL ZCOPY(lk,A(j,j),1,Ab(1,j),1)
         ENDDO
 
      ENDIF
!
      Work(1) = lwmin
!
!     End of ZHETRD_HE2HB
!
      END SUBROUTINE ZHETRD_HE2HB
