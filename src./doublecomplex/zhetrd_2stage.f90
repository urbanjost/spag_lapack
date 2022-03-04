!*==zhetrd_2stage.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZHETRD_2STAGE
!
!  @precisions fortran z -> s d c
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHETRD_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrd_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrd_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrd_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETRD_2STAGE( VECT, UPLO, N, A, LDA, D, E, TAU,
!                                 HOUS2, LHOUS2, WORK, LWORK, INFO )
!
!       IMPLICIT NONE
!
!      .. Scalar Arguments ..
!       CHARACTER          VECT, UPLO
!       INTEGER            N, LDA, LWORK, LHOUS2, INFO
!      ..
!      .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), E( * )
!       COMPLEX*16         A( LDA, * ), TAU( * ),
!                          HOUS2( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHETRD_2STAGE reduces a complex Hermitian matrix A to real symmetric
!> tridiagonal form T by a unitary similarity transformation:
!> Q1**H Q2**H* A * Q2 * Q1 = T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] VECT
!> \verbatim
!>          VECT is CHARACTER*1
!>          = 'N':  No need for the Housholder representation,
!>                  in particular for the second stage (Band to
!>                  tridiagonal) and thus LHOUS2 is of size max(1, 4*N);
!>          = 'V':  the Householder representation is needed to
!>                  either generate Q1 Q2 or to apply Q1 Q2,
!>                  then LHOUS2 is to be queried and computed.
!>                  (NOT AVAILABLE IN THIS RELEASE).
!> \endverbatim
!>
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
!>          On exit, if UPLO = 'U', the band superdiagonal
!>          of A are overwritten by the corresponding elements of the
!>          internal band-diagonal matrix AB, and the elements above
!>          the KD superdiagonal, with the array TAU, represent the unitary
!>          matrix Q1 as a product of elementary reflectors; if UPLO
!>          = 'L', the diagonal and band subdiagonal of A are over-
!>          written by the corresponding elements of the internal band-diagonal
!>          matrix AB, and the elements below the KD subdiagonal, with
!>          the array TAU, represent the unitary matrix Q1 as a product
!>          of elementary reflectors. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The off-diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (N-KD)
!>          The scalar factors of the elementary reflectors of
!>          the first stage (see Further Details).
!> \endverbatim
!>
!> \param[out] HOUS2
!> \verbatim
!>          HOUS2 is COMPLEX*16 array, dimension (LHOUS2)
!>          Stores the Householder representation of the stage2
!>          band to tridiagonal.
!> \endverbatim
!>
!> \param[in] LHOUS2
!> \verbatim
!>          LHOUS2 is INTEGER
!>          The dimension of the array HOUS2.
!>          If LWORK = -1, or LHOUS2 = -1,
!>          then a query is assumed; the routine
!>          only calculates the optimal size of the HOUS2 array, returns
!>          this value as the first entry of the HOUS2 array, and no error
!>          message related to LHOUS2 is issued by XERBLA.
!>          If VECT='N', LHOUS2 = max(1, 4*n);
!>          if VECT='V', option not yet available.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK = MAX(1, dimension)
!>          If LWORK = -1, or LHOUS2=-1,
!>          then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!>          LWORK = MAX(1, dimension) where
!>          dimension   = max(stage1,stage2) + (KD+1)*N
!>                      = N*KD + N*max(KD+1,FACTOPTNB)
!>                        + max(2*KD*KD, KD*NTHREADS)
!>                        + (KD+1)*N
!>          where KD is the blocking size of the reduction,
!>          FACTOPTNB is the blocking used by the QR or LQ
!>          algorithm, usually FACTOPTNB=128 is a good choice
!>          NTHREADS is the number of threads used when
!>          openMP compilation is enabled, otherwise =1.
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
!  =====================================================================
      SUBROUTINE ZHETRD_2STAGE(Vect,Uplo,N,A,Lda,D,E,Tau,Hous2,Lhous2,  &
     &                         Work,Lwork,Info)
!
      USE F77KINDS                        
      USE S_ILAENV2STAGE
      USE S_LSAME
      USE S_XERBLA
      USE S_ZHETRD_HB2ST
      USE S_ZHETRD_HE2HB
      IMPLICIT NONE
!*--ZHETRD_2STAGE235
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Vect
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(*) :: Hous2
      INTEGER :: Lhous2
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: abpos , ib , kd , ldab , lhmin , lwmin , lwrk , wpos
      LOGICAL :: lquery , upper , wantq
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      Info = 0
      wantq = LSAME(Vect,'V')
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1) .OR. (Lhous2==-1)
!
!     Determine the block size, the workspace size and the hous size.
!
      kd = ILAENV2STAGE(1,'ZHETRD_2STAGE',Vect,N,-1,-1,-1)
      ib = ILAENV2STAGE(2,'ZHETRD_2STAGE',Vect,N,kd,-1,-1)
      lhmin = ILAENV2STAGE(3,'ZHETRD_2STAGE',Vect,N,kd,ib,-1)
      lwmin = ILAENV2STAGE(4,'ZHETRD_2STAGE',Vect,N,kd,ib,-1)
!      WRITE(*,*),'ZHETRD_2STAGE N KD UPLO LHMIN LWMIN ',N, KD, UPLO,
!     $            LHMIN, LWMIN
!
      IF ( .NOT.LSAME(Vect,'N') ) THEN
         Info = -1
      ELSEIF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Lhous2<lhmin .AND. .NOT.lquery ) THEN
         Info = -10
      ELSEIF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
         Info = -12
      ENDIF
!
      IF ( Info==0 ) THEN
         Hous2(1) = lhmin
         Work(1) = lwmin
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHETRD_2STAGE',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
!     Determine pointer position
!
      ldab = kd + 1
      lwrk = Lwork - ldab*N
      abpos = 1
      wpos = abpos + ldab*N
      CALL ZHETRD_HE2HB(Uplo,N,kd,A,Lda,Work(abpos),ldab,Tau,Work(wpos),&
     &                  lwrk,Info)
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHETRD_HE2HB',-Info)
         RETURN
      ENDIF
      CALL ZHETRD_HB2ST('Y',Vect,Uplo,N,kd,Work(abpos),ldab,D,E,Hous2,  &
     &                  Lhous2,Work(wpos),lwrk,Info)
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHETRD_HB2ST',-Info)
         RETURN
      ENDIF
!
!
      Hous2(1) = lhmin
      Work(1) = lwmin
!
!     End of ZHETRD_2STAGE
!
      END SUBROUTINE ZHETRD_2STAGE