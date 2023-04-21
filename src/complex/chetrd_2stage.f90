!*==chetrd_2stage.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CHETRD_2STAGE
!
!  @generated from zhetrd_2stage.f, fortran z -> c, Sun Nov  6 19:34:06 2016
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHETRD_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrd_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrd_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrd_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHETRD_2STAGE( VECT, UPLO, N, A, LDA, D, E, TAU,
!                                 HOUS2, LHOUS2, WORK, LWORK, INFO )
!
!       IMPLICIT NONE
!
!      .. Scalar Arguments ..
!       CHARACTER          VECT, UPLO
!       INTEGER            N, LDA, LWORK, LHOUS2, INFO
!      ..
!      .. Array Arguments ..
!       REAL               D( * ), E( * )
!       COMPLEX            A( LDA, * ), TAU( * ),
!                          HOUS2( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHETRD_2STAGE reduces a complex Hermitian matrix A to real symmetric
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          D is REAL array, dimension (N)
!>          The diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The off-diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (N-KD)
!>          The scalar factors of the elementary reflectors of
!>          the first stage (see Further Details).
!> \endverbatim
!>
!> \param[out] HOUS2
!> \verbatim
!>          HOUS2 is COMPLEX array, dimension (LHOUS2)
!>          Stores the Householder representation of the stage2
!>          band to tridiagonal.
!> \endverbatim
!>
!> \param[in] LHOUS2
!> \verbatim
!>          LHOUS2 is INTEGER
!>          The dimension of the array HOUS2.
!>          If LWORK = -1, or LHOUS2=-1,
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
!>          WORK is COMPLEX array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK = MAX(1, dimension)
!>          If LWORK = -1, or LHOUS2 = -1,
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
!> \ingroup complexHEcomputational
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
      SUBROUTINE CHETRD_2STAGE(Vect,Uplo,N,A,Lda,D,E,Tau,Hous2,Lhous2,  &
     &                         Work,Lwork,Info)
!
      IMPLICIT NONE
!*--CHETRD_2STAGE229
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      CHARACTER Vect , Uplo
      INTEGER N , Lda , Lwork , Lhous2 , Info
!     ..
!     .. Array Arguments ..
      REAL D(*) , E(*)
      COMPLEX A(Lda,*) , Tau(*) , Hous2(*) , Work(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , upper , wantq
      INTEGER kd , ib , lwmin , lhmin , lwrk , ldab , wpos , abpos
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , CHETRD_HE2HB , CHETRD_HB2ST
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV2STAGE
      EXTERNAL LSAME , ILAENV2STAGE
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
      kd = ILAENV2STAGE(1,'CHETRD_2STAGE',Vect,N,-1,-1,-1)
      ib = ILAENV2STAGE(2,'CHETRD_2STAGE',Vect,N,kd,-1,-1)
      lhmin = ILAENV2STAGE(3,'CHETRD_2STAGE',Vect,N,kd,ib,-1)
      lwmin = ILAENV2STAGE(4,'CHETRD_2STAGE',Vect,N,kd,ib,-1)
!      WRITE(*,*),'CHETRD_2STAGE N KD UPLO LHMIN LWMIN ',N, KD, UPLO,
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
         CALL XERBLA('CHETRD_2STAGE',-Info)
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
      CALL CHETRD_HE2HB(Uplo,N,kd,A,Lda,Work(abpos),ldab,Tau,Work(wpos),&
     &                  lwrk,Info)
      IF ( Info/=0 ) THEN
         CALL XERBLA('CHETRD_HE2HB',-Info)
         RETURN
      ENDIF
      CALL CHETRD_HB2ST('Y',Vect,Uplo,N,kd,Work(abpos),ldab,D,E,Hous2,  &
     &                  Lhous2,Work(wpos),lwrk,Info)
      IF ( Info/=0 ) THEN
         CALL XERBLA('CHETRD_HB2ST',-Info)
         RETURN
      ENDIF
!
!
      Hous2(1) = lhmin
      Work(1) = lwmin
!
!     End of CHETRD_2STAGE
!
      END SUBROUTINE CHETRD_2STAGE
