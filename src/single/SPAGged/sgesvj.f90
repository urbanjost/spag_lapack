!*==sgesvj.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SGESVJ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGESVJ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesvj.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesvj.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesvj.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V,
!                          LDV, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N
!       CHARACTER*1        JOBA, JOBU, JOBV
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), SVA( N ), V( LDV, * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGESVJ computes the singular value decomposition (SVD) of a real
!> M-by-N matrix A, where M >= N. The SVD of A is written as
!>                                    [++]   [xx]   [x0]   [xx]
!>              A = U * SIGMA * V^t,  [++] = [xx] * [ox] * [xx]
!>                                    [++]   [xx]
!> where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal
!> matrix, and V is an N-by-N orthogonal matrix. The diagonal elements
!> of SIGMA are the singular values of A. The columns of U and V are the
!> left and the right singular vectors of A, respectively.
!> SGESVJ can sometimes compute tiny singular values and their singular vectors much
!> more accurately than other SVD routines, see below under Further Details.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBA
!> \verbatim
!>          JOBA is CHARACTER*1
!>          Specifies the structure of A.
!>          = 'L': The input matrix A is lower triangular;
!>          = 'U': The input matrix A is upper triangular;
!>          = 'G': The input matrix A is general M-by-N matrix, M >= N.
!> \endverbatim
!>
!> \param[in] JOBU
!> \verbatim
!>          JOBU is CHARACTER*1
!>          Specifies whether to compute the left singular vectors
!>          (columns of U):
!>          = 'U': The left singular vectors corresponding to the nonzero
!>                 singular values are computed and returned in the leading
!>                 columns of A. See more details in the description of A.
!>                 The default numerical orthogonality threshold is set to
!>                 approximately TOL=CTOL*EPS, CTOL=SQRT(M), EPS=SLAMCH('E').
!>          = 'C': Analogous to JOBU='U', except that user can control the
!>                 level of numerical orthogonality of the computed left
!>                 singular vectors. TOL can be set to TOL = CTOL*EPS, where
!>                 CTOL is given on input in the array WORK.
!>                 No CTOL smaller than ONE is allowed. CTOL greater
!>                 than 1 / EPS is meaningless. The option 'C'
!>                 can be used if M*EPS is satisfactory orthogonality
!>                 of the computed left singular vectors, so CTOL=M could
!>                 save few sweeps of Jacobi rotations.
!>                 See the descriptions of A and WORK(1).
!>          = 'N': The matrix U is not computed. However, see the
!>                 description of A.
!> \endverbatim
!>
!> \param[in] JOBV
!> \verbatim
!>          JOBV is CHARACTER*1
!>          Specifies whether to compute the right singular vectors, that
!>          is, the matrix V:
!>          = 'V':  the matrix V is computed and returned in the array V
!>          = 'A':  the Jacobi rotations are applied to the MV-by-N
!>                  array V. In other words, the right singular vector
!>                  matrix V is not computed explicitly; instead it is
!>                  applied to an MV-by-N matrix initially stored in the
!>                  first MV rows of V.
!>          = 'N':  the matrix V is not computed and the array V is not
!>                  referenced
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the input matrix A. 1/SLAMCH('E') > M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the input matrix A.
!>          M >= N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit,
!>          If JOBU = 'U' .OR. JOBU = 'C':
!>                 If INFO = 0:
!>                 RANKA orthonormal columns of U are returned in the
!>                 leading RANKA columns of the array A. Here RANKA <= N
!>                 is the number of computed singular values of A that are
!>                 above the underflow threshold SLAMCH('S'). The singular
!>                 vectors corresponding to underflowed or zero singular
!>                 values are not computed. The value of RANKA is returned
!>                 in the array WORK as RANKA=NINT(WORK(2)). Also see the
!>                 descriptions of SVA and WORK. The computed columns of U
!>                 are mutually numerically orthogonal up to approximately
!>                 TOL=SQRT(M)*EPS (default); or TOL=CTOL*EPS (JOBU = 'C'),
!>                 see the description of JOBU.
!>                 If INFO > 0,
!>                 the procedure SGESVJ did not converge in the given number
!>                 of iterations (sweeps). In that case, the computed
!>                 columns of U may not be orthogonal up to TOL. The output
!>                 U (stored in A), SIGMA (given by the computed singular
!>                 values in SVA(1:N)) and V is still a decomposition of the
!>                 input matrix A in the sense that the residual
!>                 ||A-SCALE*U*SIGMA*V^T||_2 / ||A||_2 is small.
!>          If JOBU = 'N':
!>                 If INFO = 0:
!>                 Note that the left singular vectors are 'for free' in the
!>                 one-sided Jacobi SVD algorithm. However, if only the
!>                 singular values are needed, the level of numerical
!>                 orthogonality of U is not an issue and iterations are
!>                 stopped when the columns of the iterated matrix are
!>                 numerically orthogonal up to approximately M*EPS. Thus,
!>                 on exit, A contains the columns of U scaled with the
!>                 corresponding singular values.
!>                 If INFO > 0:
!>                 the procedure SGESVJ did not converge in the given number
!>                 of iterations (sweeps).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] SVA
!> \verbatim
!>          SVA is REAL array, dimension (N)
!>          On exit,
!>          If INFO = 0 :
!>          depending on the value SCALE = WORK(1), we have:
!>                 If SCALE = ONE:
!>                 SVA(1:N) contains the computed singular values of A.
!>                 During the computation SVA contains the Euclidean column
!>                 norms of the iterated matrices in the array A.
!>                 If SCALE .NE. ONE:
!>                 The singular values of A are SCALE*SVA(1:N), and this
!>                 factored representation is due to the fact that some of the
!>                 singular values of A might underflow or overflow.
!>
!>          If INFO > 0 :
!>          the procedure SGESVJ did not converge in the given number of
!>          iterations (sweeps) and SCALE*SVA(1:N) may not be accurate.
!> \endverbatim
!>
!> \param[in] MV
!> \verbatim
!>          MV is INTEGER
!>          If JOBV = 'A', then the product of Jacobi rotations in SGESVJ
!>          is applied to the first MV rows of V. See the description of JOBV.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is REAL array, dimension (LDV,N)
!>          If JOBV = 'V', then V contains on exit the N-by-N matrix of
!>                         the right singular vectors;
!>          If JOBV = 'A', then V contains the product of the computed right
!>                         singular vector matrix and the initial matrix in
!>                         the array V.
!>          If JOBV = 'N', then V is not referenced.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V, LDV >= 1.
!>          If JOBV = 'V', then LDV >= max(1,N).
!>          If JOBV = 'A', then LDV >= max(1,MV) .
!> \endverbatim
!>
!> \param[in,out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
!>          On entry,
!>          If JOBU = 'C' :
!>          WORK(1) = CTOL, where CTOL defines the threshold for convergence.
!>                    The process stops if all columns of A are mutually
!>                    orthogonal up to CTOL*EPS, EPS=SLAMCH('E').
!>                    It is required that CTOL >= ONE, i.e. it is not
!>                    allowed to force the routine to obtain orthogonality
!>                    below EPSILON.
!>          On exit,
!>          WORK(1) = SCALE is the scaling factor such that SCALE*SVA(1:N)
!>                    are the computed singular vcalues of A.
!>                    (See description of SVA().)
!>          WORK(2) = NINT(WORK(2)) is the number of the computed nonzero
!>                    singular values.
!>          WORK(3) = NINT(WORK(3)) is the number of the computed singular
!>                    values that are larger than the underflow threshold.
!>          WORK(4) = NINT(WORK(4)) is the number of sweeps of Jacobi
!>                    rotations needed for numerical convergence.
!>          WORK(5) = max_{i.NE.j} |COS(A(:,i),A(:,j))| in the last sweep.
!>                    This is useful information in cases when SGESVJ did
!>                    not converge, as it can be used to estimate whether
!>                    the output is still useful and for post festum analysis.
!>          WORK(6) = the largest absolute value over all sines of the
!>                    Jacobi rotation angles in the last sweep. It can be
!>                    useful for a post festum analysis.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>         length of WORK, WORK >= MAX(6,M+N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, then the i-th argument had an illegal value
!>          > 0:  SGESVJ did not converge in the maximal allowed number (30)
!>                of sweeps. The output may still be useful. See the
!>                description of WORK.
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
!> \date June 2017
!
!> \ingroup realGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> The orthogonal N-by-N matrix V is obtained as a product of Jacobi plane
!> rotations. The rotations are implemented as fast scaled rotations of
!> Anda and Park [1]. In the case of underflow of the Jacobi angle, a
!> modified Jacobi transformation of Drmac [4] is used. Pivot strategy uses
!> column interchanges of de Rijk [2]. The relative accuracy of the computed
!> singular values and the accuracy of the computed singular vectors (in
!> angle metric) is as guaranteed by the theory of Demmel and Veselic [3].
!> The condition number that determines the accuracy in the full rank case
!> is essentially min_{D=diag} kappa(A*D), where kappa(.) is the
!> spectral condition number. The best performance of this Jacobi SVD
!> procedure is achieved if used in an  accelerated version of Drmac and
!> Veselic [5,6], and it is the kernel routine in the SIGMA library [7].
!> Some tuning parameters (marked with [TP]) are available for the
!> implementer. \n
!> The computational range for the nonzero singular values is the  machine
!> number interval ( UNDERFLOW , OVERFLOW ). In extreme cases, even
!> denormalized singular values can be computed with the corresponding
!> gradual loss of accurate digits.
!>
!> \par Contributors:
!  ==================
!>
!> Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)
!>
!> \par References:
!  ================
!>
!> [1] A. A. Anda and H. Park: Fast plane rotations with dynamic scaling. \n
!>    SIAM J. matrix Anal. Appl., Vol. 15 (1994), pp. 162-174. \n\n
!> [2] P. P. M. De Rijk: A one-sided Jacobi algorithm for computing the
!>    singular value decomposition on a vector computer. \n
!>    SIAM J. Sci. Stat. Comp., Vol. 10 (1998), pp. 359-371. \n\n
!> [3] J. Demmel and K. Veselic: Jacobi method is more accurate than QR. \n
!> [4] Z. Drmac: Implementation of Jacobi rotations for accurate singular
!>    value computation in floating point arithmetic. \n
!>    SIAM J. Sci. Comp., Vol. 18 (1997), pp. 1200-1222. \n\n
!> [5] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I. \n
!>    SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342. \n
!>    LAPACK Working note 169. \n\n
!> [6] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II. \n
!>    SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362. \n
!>    LAPACK Working note 170. \n\n
!> [7] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV,
!>    QSVD, (H,K)-SVD computations.\n
!>    Department of Mathematics, University of Zagreb, 2008.
!>
!> \par Bugs, Examples and Comments:
!  =================================
!>
!> Please report all bugs and send interesting test examples and comments to
!> drmac@math.hr. Thank you.
!
!  =====================================================================
      SUBROUTINE SGESVJ(Joba,Jobu,Jobv,M,N,A,Lda,Sva,Mv,V,Ldv,Work,     &
     &                  Lwork,Info)
      USE S_ISAMAX
      USE S_LSAME
      USE S_SAXPY
      USE S_SCOPY
      USE S_SDOT
      USE S_SGSVJ0
      USE S_SGSVJ1
      USE S_SLAMCH
      USE S_SLASCL
      USE S_SLASET
      USE S_SLASSQ
      USE S_SNRM2
      USE S_SROTM
      USE S_SSCAL
      USE S_SSWAP
      USE S_XERBLA
      IMPLICIT NONE
!*--SGESVJ343
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , HALF = 0.5E0 , ONE = 1.0E0
      INTEGER , PARAMETER  ::  NSWEEP = 30
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER(1) :: Joba
      CHARACTER(1) :: Jobu
      CHARACTER(1) :: Jobv
      INTEGER :: M
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(N) :: Sva
      INTEGER , INTENT(IN) :: Mv
      REAL , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL , INTENT(INOUT) , DIMENSION(Lwork) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: aapp , aapp0 , aapq , aaqq , apoaq , aqoap , big ,        &
     &        bigtheta , cs , ctol , epsln , large , mxaapq , mxsinj ,  &
     &        rootbig , rooteps , rootsfmin , roottol , sfmin , skl ,   &
     &        small , sn , t , temp1 , theta , thsign , tol
      LOGICAL :: applv , goscale , lower , lsvec , noscale , rotok ,    &
     &           rsvec , uctol , upper
      INTEGER :: blskip , emptsw , i , ibr , ierr , igl , ijblsk , ir1 ,&
     &           iswrot , jbc , jgl , kbl , lkahead , mvl , n2 , n34 ,  &
     &           n4 , nbl , notrot , p , pskipped , q , rowskip , swband
      REAL , DIMENSION(5) :: fastr
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     from BLAS
!     from LAPACK
!     ..
!     .. External Subroutines ..
!     ..
!     from BLAS
!     from LAPACK
!
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      lsvec = LSAME(Jobu,'U')
      uctol = LSAME(Jobu,'C')
      rsvec = LSAME(Jobv,'V')
      applv = LSAME(Jobv,'A')
      upper = LSAME(Joba,'U')
      lower = LSAME(Joba,'L')
!
      IF ( .NOT.(upper .OR. lower .OR. LSAME(Joba,'G')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(lsvec .OR. uctol .OR. LSAME(Jobu,'N')) ) THEN
         Info = -2
      ELSEIF ( .NOT.(rsvec .OR. applv .OR. LSAME(Jobv,'N')) ) THEN
         Info = -3
      ELSEIF ( M<0 ) THEN
         Info = -4
      ELSEIF ( (N<0) .OR. (N>M) ) THEN
         Info = -5
      ELSEIF ( Lda<M ) THEN
         Info = -7
      ELSEIF ( Mv<0 ) THEN
         Info = -9
      ELSEIF ( (rsvec .AND. (Ldv<N)) .OR. (applv .AND. (Ldv<Mv)) ) THEN
         Info = -11
      ELSEIF ( uctol .AND. (Work(1)<=ONE) ) THEN
         Info = -12
      ELSEIF ( Lwork<MAX(M+N,6) ) THEN
         Info = -13
      ELSE
         Info = 0
      ENDIF
!
!     #:(
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGESVJ',-Info)
         RETURN
      ENDIF
!
! #:) Quick return for void matrix
!
      IF ( (M==0) .OR. (N==0) ) RETURN
!
!     Set numerical parameters
!     The stopping criterion for Jacobi rotations is
!
!     max_{i<>j}|A(:,i)^T * A(:,j)|/(||A(:,i)||*||A(:,j)||) < CTOL*EPS
!
!     where EPS is the round-off and CTOL is defined as follows:
!
      IF ( uctol ) THEN
!        ... user controlled
         ctol = Work(1)
!        ... default
      ELSEIF ( lsvec .OR. rsvec .OR. applv ) THEN
         ctol = SQRT(FLOAT(M))
      ELSE
         ctol = FLOAT(M)
      ENDIF
!     ... and the machine dependent parameters are
![!]  (Make sure that SLAMCH() works properly on the target machine.)
!
      epsln = SLAMCH('Epsilon')
      rooteps = SQRT(epsln)
      sfmin = SLAMCH('SafeMinimum')
      rootsfmin = SQRT(sfmin)
      small = sfmin/epsln
      big = SLAMCH('Overflow')
!     BIG         = ONE    / SFMIN
      rootbig = ONE/rootsfmin
      large = big/SQRT(FLOAT(M*N))
      bigtheta = ONE/rooteps
!
      tol = ctol*epsln
      roottol = SQRT(tol)
!
      IF ( FLOAT(M)*epsln>=ONE ) THEN
         Info = -4
         CALL XERBLA('SGESVJ',-Info)
         RETURN
      ENDIF
!
!     Initialize the right singular vector matrix.
!
      IF ( rsvec ) THEN
         mvl = N
         CALL SLASET('A',mvl,N,ZERO,ONE,V,Ldv)
      ELSEIF ( applv ) THEN
         mvl = Mv
      ENDIF
      rsvec = rsvec .OR. applv
!
!     Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N )
!(!)  If necessary, scale A to protect the largest singular value
!     from overflow. It is possible that saving the largest singular
!     value destroys the information about the small ones.
!     This initial scaling is almost minimal in the sense that the
!     goal is to make sure that no column norm overflows, and that
!     SQRT(N)*max_i SVA(i) does not overflow. If INFinite entries
!     in A are detected, the procedure returns with INFO=-6.
!
      skl = ONE/SQRT(FLOAT(M)*FLOAT(N))
      noscale = .TRUE.
      goscale = .TRUE.
!
      IF ( lower ) THEN
!        the input matrix is M-by-N lower triangular (trapezoidal)
         DO p = 1 , N
            aapp = ZERO
            aaqq = ONE
            CALL SLASSQ(M-p+1,A(p,p),1,aapp,aaqq)
            IF ( aapp>big ) THEN
               Info = -6
               CALL XERBLA('SGESVJ',-Info)
               RETURN
            ENDIF
            aaqq = SQRT(aaqq)
            IF ( (aapp<(big/aaqq)) .AND. noscale ) THEN
               Sva(p) = aapp*aaqq
            ELSE
               noscale = .FALSE.
               Sva(p) = aapp*(aaqq*skl)
               IF ( goscale ) THEN
                  goscale = .FALSE.
                  DO q = 1 , p - 1
                     Sva(q) = Sva(q)*skl
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ELSEIF ( upper ) THEN
!        the input matrix is M-by-N upper triangular (trapezoidal)
         DO p = 1 , N
            aapp = ZERO
            aaqq = ONE
            CALL SLASSQ(p,A(1,p),1,aapp,aaqq)
            IF ( aapp>big ) THEN
               Info = -6
               CALL XERBLA('SGESVJ',-Info)
               RETURN
            ENDIF
            aaqq = SQRT(aaqq)
            IF ( (aapp<(big/aaqq)) .AND. noscale ) THEN
               Sva(p) = aapp*aaqq
            ELSE
               noscale = .FALSE.
               Sva(p) = aapp*(aaqq*skl)
               IF ( goscale ) THEN
                  goscale = .FALSE.
                  DO q = 1 , p - 1
                     Sva(q) = Sva(q)*skl
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ELSE
!        the input matrix is M-by-N general dense
         DO p = 1 , N
            aapp = ZERO
            aaqq = ONE
            CALL SLASSQ(M,A(1,p),1,aapp,aaqq)
            IF ( aapp>big ) THEN
               Info = -6
               CALL XERBLA('SGESVJ',-Info)
               RETURN
            ENDIF
            aaqq = SQRT(aaqq)
            IF ( (aapp<(big/aaqq)) .AND. noscale ) THEN
               Sva(p) = aapp*aaqq
            ELSE
               noscale = .FALSE.
               Sva(p) = aapp*(aaqq*skl)
               IF ( goscale ) THEN
                  goscale = .FALSE.
                  DO q = 1 , p - 1
                     Sva(q) = Sva(q)*skl
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
      IF ( noscale ) skl = ONE
!
!     Move the smaller part of the spectrum from the underflow threshold
!(!)  Start by determining the position of the nonzero entries of the
!     array SVA() relative to ( SFMIN, BIG ).
!
      aapp = ZERO
      aaqq = big
      DO p = 1 , N
         IF ( Sva(p)/=ZERO ) aaqq = MIN(aaqq,Sva(p))
         aapp = MAX(aapp,Sva(p))
      ENDDO
!
! #:) Quick return for zero matrix
!
      IF ( aapp==ZERO ) THEN
         IF ( lsvec ) CALL SLASET('G',M,N,ZERO,ONE,A,Lda)
         Work(1) = ONE
         Work(2) = ZERO
         Work(3) = ZERO
         Work(4) = ZERO
         Work(5) = ZERO
         Work(6) = ZERO
         RETURN
      ENDIF
!
! #:) Quick return for one-column matrix
!
      IF ( N==1 ) THEN
         IF ( lsvec ) CALL SLASCL('G',0,0,Sva(1),skl,M,1,A(1,1),Lda,    &
     &                            ierr)
         Work(1) = ONE/skl
         IF ( Sva(1)>=sfmin ) THEN
            Work(2) = ONE
         ELSE
            Work(2) = ZERO
         ENDIF
         Work(3) = ZERO
         Work(4) = ZERO
         Work(5) = ZERO
         Work(6) = ZERO
         RETURN
      ENDIF
!
!     Protect small singular values from underflow, and try to
!     avoid underflows/overflows in computing Jacobi rotations.
!
      sn = SQRT(sfmin/epsln)
      temp1 = SQRT(big/FLOAT(N))
      IF ( (aapp<=sn) .OR. (aaqq>=temp1) .OR.                           &
     &     ((sn<=aaqq) .AND. (aapp<=temp1)) ) THEN
         temp1 = MIN(big,temp1/aapp)
!         AAQQ  = AAQQ*TEMP1
!         AAPP  = AAPP*TEMP1
      ELSEIF ( (aaqq<=sn) .AND. (aapp<=temp1) ) THEN
         temp1 = MIN(sn/aaqq,big/(aapp*SQRT(FLOAT(N))))
!         AAQQ  = AAQQ*TEMP1
!         AAPP  = AAPP*TEMP1
      ELSEIF ( (aaqq>=sn) .AND. (aapp>=temp1) ) THEN
         temp1 = MAX(sn/aaqq,temp1/aapp)
!         AAQQ  = AAQQ*TEMP1
!         AAPP  = AAPP*TEMP1
      ELSEIF ( (aaqq<=sn) .AND. (aapp>=temp1) ) THEN
         temp1 = MIN(sn/aaqq,big/(SQRT(FLOAT(N))*aapp))
!         AAQQ  = AAQQ*TEMP1
!         AAPP  = AAPP*TEMP1
      ELSE
         temp1 = ONE
      ENDIF
!
!     Scale, if necessary
!
      IF ( temp1/=ONE ) CALL SLASCL('G',0,0,ONE,temp1,N,1,Sva,N,ierr)
      skl = temp1*skl
      IF ( skl/=ONE ) THEN
         CALL SLASCL(Joba,0,0,ONE,skl,M,N,A,Lda,ierr)
         skl = ONE/skl
      ENDIF
!
!     Row-cyclic Jacobi SVD algorithm with column pivoting
!
      emptsw = (N*(N-1))/2
      notrot = 0
      fastr(1) = ZERO
!
!     A is represented in factored form A = A * diag(WORK), where diag(WORK)
!     is initialized to identity. WORK is updated during fast scaled
!     rotations.
!
      DO q = 1 , N
         Work(q) = ONE
      ENDDO
!
!
      swband = 3
![TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective
!     if SGESVJ is used as a computational routine in the preconditioned
!     Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure
!     works on pivots inside a band-like region around the diagonal.
!     The boundaries are determined dynamically, based on the number of
!     pivots above a threshold.
!
      kbl = MIN(8,N)
![TP] KBL is a tuning parameter that defines the tile size in the
!     tiling of the p-q loops of pivot pairs. In general, an optimal
!     value of KBL depends on the matrix dimensions and on the
!     parameters of the computer's memory.
!
      nbl = N/kbl
      IF ( (nbl*kbl)/=N ) nbl = nbl + 1
!
      blskip = kbl**2
![TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.
!
      rowskip = MIN(5,kbl)
![TP] ROWSKIP is a tuning parameter.
!
      lkahead = 1
![TP] LKAHEAD is a tuning parameter.
!
!     Quasi block transformations, using the lower (upper) triangular
!     structure of the input matrix. The quasi-block-cycling usually
!     invokes cubic convergence. Big part of this cycle is done inside
!     canonical subspaces of dimensions less than M.
!
      IF ( (lower .OR. upper) .AND. (N>MAX(64,4*kbl)) ) THEN
![TP] The number of partition levels and the actual partition are
!     tuning parameters.
         n4 = N/4
         n2 = N/2
         n34 = 3*n4
         IF ( applv ) THEN
            q = 0
         ELSE
            q = 1
         ENDIF
!
         IF ( lower ) THEN
!
!     This works very well on lower triangular matrices, in particular
!     in the framework of the preconditioned Jacobi SVD (xGEJSV).
!     The idea is simple:
!     [+ 0 0 0]   Note that Jacobi transformations of [0 0]
!     [+ + 0 0]                                       [0 0]
!     [+ + x 0]   actually work on [x 0]              [x 0]
!     [+ + x x]                    [x x].             [x x]
!
            CALL SGSVJ0(Jobv,M-n34,N-n34,A(n34+1,n34+1),Lda,Work(n34+1),&
     &                  Sva(n34+1),mvl,V(n34*q+1,n34+1),Ldv,epsln,sfmin,&
     &                  tol,2,Work(N+1),Lwork-N,ierr)
!
            CALL SGSVJ0(Jobv,M-n2,n34-n2,A(n2+1,n2+1),Lda,Work(n2+1),   &
     &                  Sva(n2+1),mvl,V(n2*q+1,n2+1),Ldv,epsln,sfmin,   &
     &                  tol,2,Work(N+1),Lwork-N,ierr)
!
            CALL SGSVJ1(Jobv,M-n2,N-n2,n4,A(n2+1,n2+1),Lda,Work(n2+1),  &
     &                  Sva(n2+1),mvl,V(n2*q+1,n2+1),Ldv,epsln,sfmin,   &
     &                  tol,1,Work(N+1),Lwork-N,ierr)
!
            CALL SGSVJ0(Jobv,M-n4,n2-n4,A(n4+1,n4+1),Lda,Work(n4+1),    &
     &                  Sva(n4+1),mvl,V(n4*q+1,n4+1),Ldv,epsln,sfmin,   &
     &                  tol,1,Work(N+1),Lwork-N,ierr)
!
            CALL SGSVJ0(Jobv,M,n4,A,Lda,Work,Sva,mvl,V,Ldv,epsln,sfmin, &
     &                  tol,1,Work(N+1),Lwork-N,ierr)
!
            CALL SGSVJ1(Jobv,M,n2,n4,A,Lda,Work,Sva,mvl,V,Ldv,epsln,    &
     &                  sfmin,tol,1,Work(N+1),Lwork-N,ierr)
!
!
         ELSEIF ( upper ) THEN
!
!
            CALL SGSVJ0(Jobv,n4,n4,A,Lda,Work,Sva,mvl,V,Ldv,epsln,sfmin,&
     &                  tol,2,Work(N+1),Lwork-N,ierr)
!
            CALL SGSVJ0(Jobv,n2,n4,A(1,n4+1),Lda,Work(n4+1),Sva(n4+1),  &
     &                  mvl,V(n4*q+1,n4+1),Ldv,epsln,sfmin,tol,1,       &
     &                  Work(N+1),Lwork-N,ierr)
!
            CALL SGSVJ1(Jobv,n2,n2,n4,A,Lda,Work,Sva,mvl,V,Ldv,epsln,   &
     &                  sfmin,tol,1,Work(N+1),Lwork-N,ierr)
!
            CALL SGSVJ0(Jobv,n2+n4,n4,A(1,n2+1),Lda,Work(n2+1),Sva(n2+1)&
     &                  ,mvl,V(n2*q+1,n2+1),Ldv,epsln,sfmin,tol,1,      &
     &                  Work(N+1),Lwork-N,ierr)
 
         ENDIF
!
      ENDIF
!
!     .. Row-cyclic pivot strategy with de Rijk's pivoting ..
!
      DO i = 1 , NSWEEP
!
!     .. go go go ...
!
         mxaapq = ZERO
         mxsinj = ZERO
         iswrot = 0
!
         notrot = 0
         pskipped = 0
!
!     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs
!     1 <= p < q <= N. This is the first step toward a blocked implementation
!     of the rotations. New implementation, based on block transformations,
!     is under development.
!
         DO ibr = 1 , nbl
!
            igl = (ibr-1)*kbl + 1
!
            DO ir1 = 0 , MIN(lkahead,nbl-ibr)
!
               igl = igl + ir1*kbl
!
               DO p = igl , MIN(igl+kbl-1,N-1)
!
!     .. de Rijk's pivoting
!
                  q = ISAMAX(N-p+1,Sva(p),1) + p - 1
                  IF ( p/=q ) THEN
                     CALL SSWAP(M,A(1,p),1,A(1,q),1)
                     IF ( rsvec ) CALL SSWAP(mvl,V(1,p),1,V(1,q),1)
                     temp1 = Sva(p)
                     Sva(p) = Sva(q)
                     Sva(q) = temp1
                     temp1 = Work(p)
                     Work(p) = Work(q)
                     Work(q) = temp1
                  ENDIF
!
                  IF ( ir1==0 ) THEN
!
!        Column norms are periodically updated by explicit
!        norm computation.
!        Caveat:
!        Unfortunately, some BLAS implementations compute SNRM2(M,A(1,p),1)
!        as SQRT(SDOT(M,A(1,p),1,A(1,p),1)), which may cause the result to
!        overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to
!        underflow for ||A(:,p)||_2 < SQRT(underflow_threshold).
!        Hence, SNRM2 cannot be trusted, not even in the case when
!        the true norm is far from the under(over)flow boundaries.
!        If properly implemented SNRM2 is available, the IF-THEN-ELSE
!        below should read "AAPP = SNRM2( M, A(1,p), 1 ) * WORK(p)".
!
                     IF ( (Sva(p)<rootbig) .AND. (Sva(p)>rootsfmin) )   &
     &                    THEN
                        Sva(p) = SNRM2(M,A(1,p),1)*Work(p)
                     ELSE
                        temp1 = ZERO
                        aapp = ONE
                        CALL SLASSQ(M,A(1,p),1,temp1,aapp)
                        Sva(p) = temp1*SQRT(aapp)*Work(p)
                     ENDIF
                     aapp = Sva(p)
                  ELSE
                     aapp = Sva(p)
                  ENDIF
!
                  IF ( aapp>ZERO ) THEN
!
                     pskipped = 0
!
                     DO q = p + 1 , MIN(igl+kbl-1,N)
!
                        aaqq = Sva(q)
!
                        IF ( aaqq>ZERO ) THEN
!
                           aapp0 = aapp
                           IF ( aaqq>=ONE ) THEN
                              rotok = (small*aapp)<=aaqq
                              IF ( aapp<(big/aaqq) ) THEN
                                 aapq = (SDOT(M,A(1,p),1,A(1,q),1)      &
     &                                  *Work(p)*Work(q)/aaqq)/aapp
                              ELSE
                                 CALL SCOPY(M,A(1,p),1,Work(N+1),1)
                                 CALL SLASCL('G',0,0,aapp,Work(p),M,1,  &
     &                              Work(N+1),Lda,ierr)
                                 aapq = SDOT(M,Work(N+1),1,A(1,q),1)    &
     &                                  *Work(q)/aaqq
                              ENDIF
                           ELSE
                              rotok = aapp<=(aaqq/small)
                              IF ( aapp>(small/aaqq) ) THEN
                                 aapq = (SDOT(M,A(1,p),1,A(1,q),1)      &
     &                                  *Work(p)*Work(q)/aaqq)/aapp
                              ELSE
                                 CALL SCOPY(M,A(1,q),1,Work(N+1),1)
                                 CALL SLASCL('G',0,0,aaqq,Work(q),M,1,  &
     &                              Work(N+1),Lda,ierr)
                                 aapq = SDOT(M,Work(N+1),1,A(1,p),1)    &
     &                                  *Work(p)/aapp
                              ENDIF
                           ENDIF
!
                           mxaapq = MAX(mxaapq,ABS(aapq))
!
!        TO rotate or NOT to rotate, THAT is the question ...
!
                           IF ( ABS(aapq)>tol ) THEN
!
!           .. rotate
![RTD]      ROTATED = ROTATED + ONE
!
                              IF ( ir1==0 ) THEN
                                 notrot = 0
                                 pskipped = 0
                                 iswrot = iswrot + 1
                              ENDIF
!
                              IF ( rotok ) THEN
!
                                 aqoap = aaqq/aapp
                                 apoaq = aapp/aaqq
                                 theta = -HALF*ABS(aqoap-apoaq)/aapq
!
                                 IF ( ABS(theta)>bigtheta ) THEN
!
                                    t = HALF/theta
                                    fastr(3) = t*Work(p)/Work(q)
                                    fastr(4) = -t*Work(q)/Work(p)
                                    CALL SROTM(M,A(1,p),1,A(1,q),1,     &
     &                                 fastr)
                                    IF ( rsvec )                        &
     &                                 CALL SROTM(mvl,V(1,p),1,V(1,q),1,&
     &                                 fastr)
                                    Sva(q)                              &
     &                                 = aaqq*SQRT(MAX(ZERO,ONE+t*apoaq*&
     &                                 aapq))
                                    aapp = aapp*SQRT                    &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq))
                                    mxsinj = MAX(mxsinj,ABS(t))
!
                                 ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                    thsign = -SIGN(ONE,aapq)
                                    t = ONE/                            &
     &                                  (theta+thsign*SQRT(ONE+theta*   &
     &                                  theta))
                                    cs = SQRT(ONE/(ONE+t*t))
                                    sn = t*cs
!
                                    mxsinj = MAX(mxsinj,ABS(sn))
                                    Sva(q)                              &
     &                                 = aaqq*SQRT(MAX(ZERO,ONE+t*apoaq*&
     &                                 aapq))
                                    aapp = aapp*SQRT                    &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq))
!
                                    apoaq = Work(p)/Work(q)
                                    aqoap = Work(q)/Work(p)
                                    IF ( Work(p)>=ONE ) THEN
                                       IF ( Work(q)>=ONE ) THEN
                                         fastr(3) = t*apoaq
                                         fastr(4) = -t*aqoap
                                         Work(p) = Work(p)*cs
                                         Work(q) = Work(q)*cs
                                         CALL SROTM(M,A(1,p),1,A(1,q),1,&
     &                                      fastr)
                                         IF ( rsvec )                   &
     &                                      CALL SROTM(mvl,V(1,p),1,    &
     &                                      V(1,q),1,fastr)
                                       ELSE
                                         CALL SAXPY(M,-t*aqoap,A(1,q),1,&
     &                                      A(1,p),1)
                                         CALL SAXPY(M,cs*sn*apoaq,A(1,p)&
     &                                      ,1,A(1,q),1)
                                         Work(p) = Work(p)*cs
                                         Work(q) = Work(q)/cs
                                         IF ( rsvec ) THEN
                                         CALL SAXPY(mvl,-t*aqoap,V(1,q),&
     &                                      1,V(1,p),1)
                                         CALL SAXPY(mvl,cs*sn*apoaq,    &
     &                                      V(1,p),1,V(1,q),1)
                                         ENDIF
                                       ENDIF
                                    ELSEIF ( Work(q)>=ONE ) THEN
                                       CALL SAXPY(M,t*apoaq,A(1,p),1,   &
     &                                    A(1,q),1)
                                       CALL SAXPY(M,-cs*sn*aqoap,A(1,q),&
     &                                    1,A(1,p),1)
                                       Work(p) = Work(p)/cs
                                       Work(q) = Work(q)*cs
                                       IF ( rsvec ) THEN
                                         CALL SAXPY(mvl,t*apoaq,V(1,p), &
     &                                      1,V(1,q),1)
                                         CALL SAXPY(mvl,-cs*sn*aqoap,   &
     &                                      V(1,q),1,V(1,p),1)
                                       ENDIF
                                    ELSEIF ( Work(p)>=Work(q) ) THEN
                                       CALL SAXPY(M,-t*aqoap,A(1,q),1,  &
     &                                    A(1,p),1)
                                       CALL SAXPY(M,cs*sn*apoaq,A(1,p), &
     &                                    1,A(1,q),1)
                                       Work(p) = Work(p)*cs
                                       Work(q) = Work(q)/cs
                                       IF ( rsvec ) THEN
                                         CALL SAXPY(mvl,-t*aqoap,V(1,q),&
     &                                      1,V(1,p),1)
                                         CALL SAXPY(mvl,cs*sn*apoaq,    &
     &                                      V(1,p),1,V(1,q),1)
                                       ENDIF
                                    ELSE
                                       CALL SAXPY(M,t*apoaq,A(1,p),1,   &
     &                                    A(1,q),1)
                                       CALL SAXPY(M,-cs*sn*aqoap,A(1,q),&
     &                                    1,A(1,p),1)
                                       Work(p) = Work(p)/cs
                                       Work(q) = Work(q)*cs
                                       IF ( rsvec ) THEN
                                         CALL SAXPY(mvl,t*apoaq,V(1,p), &
     &                                      1,V(1,q),1)
                                         CALL SAXPY(mvl,-cs*sn*aqoap,   &
     &                                      V(1,q),1,V(1,p),1)
                                       ENDIF
                                    ENDIF
                                 ENDIF
!
                              ELSE
!              .. have to use modified Gram-Schmidt like transformation
                                 CALL SCOPY(M,A(1,p),1,Work(N+1),1)
                                 CALL SLASCL('G',0,0,aapp,ONE,M,1,      &
     &                              Work(N+1),Lda,ierr)
                                 CALL SLASCL('G',0,0,aaqq,ONE,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 temp1 = -aapq*Work(p)/Work(q)
                                 CALL SAXPY(M,temp1,Work(N+1),1,A(1,q), &
     &                              1)
                                 CALL SLASCL('G',0,0,ONE,aaqq,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 Sva(q)                                 &
     &                              = aaqq*SQRT(MAX(ZERO,ONE-aapq*aapq))
                                 mxsinj = MAX(mxsinj,sfmin)
                              ENDIF
!           END IF ROTOK THEN ... ELSE
!
!           In the case of cancellation in updating SVA(q), SVA(p)
!           recompute SVA(q), SVA(p).
!
                              IF ( (Sva(q)/aaqq)**2<=rooteps ) THEN
                                 IF ( (aaqq<rootbig) .AND.              &
     &                                (aaqq>rootsfmin) ) THEN
                                    Sva(q) = SNRM2(M,A(1,q),1)*Work(q)
                                 ELSE
                                    t = ZERO
                                    aaqq = ONE
                                    CALL SLASSQ(M,A(1,q),1,t,aaqq)
                                    Sva(q) = t*SQRT(aaqq)*Work(q)
                                 ENDIF
                              ENDIF
                              IF ( (aapp/aapp0)<=rooteps ) THEN
                                 IF ( (aapp<rootbig) .AND.              &
     &                                (aapp>rootsfmin) ) THEN
                                    aapp = SNRM2(M,A(1,p),1)*Work(p)
                                 ELSE
                                    t = ZERO
                                    aapp = ONE
                                    CALL SLASSQ(M,A(1,p),1,t,aapp)
                                    aapp = t*SQRT(aapp)*Work(p)
                                 ENDIF
                                 Sva(p) = aapp
                              ENDIF
!
                           ELSE
!        A(:,p) and A(:,q) already numerically orthogonal
                              IF ( ir1==0 ) notrot = notrot + 1
![RTD]      SKIPPED  = SKIPPED  + 1
                              pskipped = pskipped + 1
                           ENDIF
                        ELSE
!        A(:,q) is zero column
                           IF ( ir1==0 ) notrot = notrot + 1
                           pskipped = pskipped + 1
                        ENDIF
!
                        IF ( (i<=swband) .AND. (pskipped>rowskip) ) THEN
                           IF ( ir1==0 ) aapp = -aapp
                           notrot = 0
                           EXIT
                        ENDIF
!
                     ENDDO
!     END q-LOOP
!
!     bailed out of q-loop
!
                     Sva(p) = aapp
!
                  ELSE
                     Sva(p) = aapp
                     IF ( (ir1==0) .AND. (aapp==ZERO) )                 &
     &                    notrot = notrot + MIN(igl+kbl-1,N) - p
                  ENDIF
!
               ENDDO
!     end of the p-loop
!     end of doing the block ( ibr, ibr )
            ENDDO
!     end of ir1-loop
!
! ... go to the off diagonal blocks
!
            igl = (ibr-1)*kbl + 1
!
            DO jbc = ibr + 1 , nbl
!
               jgl = (jbc-1)*kbl + 1
!
!        doing the block at ( ibr, jbc )
!
               ijblsk = 0
               DO p = igl , MIN(igl+kbl-1,N)
!
                  aapp = Sva(p)
                  IF ( aapp>ZERO ) THEN
!
                     pskipped = 0
!
                     DO q = jgl , MIN(jgl+kbl-1,N)
!
                        aaqq = Sva(q)
                        IF ( aaqq>ZERO ) THEN
                           aapp0 = aapp
!
!     .. M x 2 Jacobi SVD ..
!
!        Safe Gram matrix computation
!
                           IF ( aaqq>=ONE ) THEN
                              IF ( aapp>=aaqq ) THEN
                                 rotok = (small*aapp)<=aaqq
                              ELSE
                                 rotok = (small*aaqq)<=aapp
                              ENDIF
                              IF ( aapp<(big/aaqq) ) THEN
                                 aapq = (SDOT(M,A(1,p),1,A(1,q),1)      &
     &                                  *Work(p)*Work(q)/aaqq)/aapp
                              ELSE
                                 CALL SCOPY(M,A(1,p),1,Work(N+1),1)
                                 CALL SLASCL('G',0,0,aapp,Work(p),M,1,  &
     &                              Work(N+1),Lda,ierr)
                                 aapq = SDOT(M,Work(N+1),1,A(1,q),1)    &
     &                                  *Work(q)/aaqq
                              ENDIF
                           ELSE
                              IF ( aapp>=aaqq ) THEN
                                 rotok = aapp<=(aaqq/small)
                              ELSE
                                 rotok = aaqq<=(aapp/small)
                              ENDIF
                              IF ( aapp>(small/aaqq) ) THEN
                                 aapq = (SDOT(M,A(1,p),1,A(1,q),1)      &
     &                                  *Work(p)*Work(q)/aaqq)/aapp
                              ELSE
                                 CALL SCOPY(M,A(1,q),1,Work(N+1),1)
                                 CALL SLASCL('G',0,0,aaqq,Work(q),M,1,  &
     &                              Work(N+1),Lda,ierr)
                                 aapq = SDOT(M,Work(N+1),1,A(1,p),1)    &
     &                                  *Work(p)/aapp
                              ENDIF
                           ENDIF
!
                           mxaapq = MAX(mxaapq,ABS(aapq))
!
!        TO rotate or NOT to rotate, THAT is the question ...
!
                           IF ( ABS(aapq)>tol ) THEN
                              notrot = 0
![RTD]      ROTATED  = ROTATED + 1
                              pskipped = 0
                              iswrot = iswrot + 1
!
                              IF ( rotok ) THEN
!
                                 aqoap = aaqq/aapp
                                 apoaq = aapp/aaqq
                                 theta = -HALF*ABS(aqoap-apoaq)/aapq
                                 IF ( aaqq>aapp0 ) theta = -theta
!
                                 IF ( ABS(theta)>bigtheta ) THEN
                                    t = HALF/theta
                                    fastr(3) = t*Work(p)/Work(q)
                                    fastr(4) = -t*Work(q)/Work(p)
                                    CALL SROTM(M,A(1,p),1,A(1,q),1,     &
     &                                 fastr)
                                    IF ( rsvec )                        &
     &                                 CALL SROTM(mvl,V(1,p),1,V(1,q),1,&
     &                                 fastr)
                                    Sva(q)                              &
     &                                 = aaqq*SQRT(MAX(ZERO,ONE+t*apoaq*&
     &                                 aapq))
                                    aapp = aapp*SQRT                    &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq))
                                    mxsinj = MAX(mxsinj,ABS(t))
                                 ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                    thsign = -SIGN(ONE,aapq)
                                    IF ( aaqq>aapp0 ) thsign = -thsign
                                    t = ONE/                            &
     &                                  (theta+thsign*SQRT(ONE+theta*   &
     &                                  theta))
                                    cs = SQRT(ONE/(ONE+t*t))
                                    sn = t*cs
                                    mxsinj = MAX(mxsinj,ABS(sn))
                                    Sva(q)                              &
     &                                 = aaqq*SQRT(MAX(ZERO,ONE+t*apoaq*&
     &                                 aapq))
                                    aapp = aapp*SQRT                    &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq))
!
                                    apoaq = Work(p)/Work(q)
                                    aqoap = Work(q)/Work(p)
                                    IF ( Work(p)>=ONE ) THEN
!
                                       IF ( Work(q)>=ONE ) THEN
                                         fastr(3) = t*apoaq
                                         fastr(4) = -t*aqoap
                                         Work(p) = Work(p)*cs
                                         Work(q) = Work(q)*cs
                                         CALL SROTM(M,A(1,p),1,A(1,q),1,&
     &                                      fastr)
                                         IF ( rsvec )                   &
     &                                      CALL SROTM(mvl,V(1,p),1,    &
     &                                      V(1,q),1,fastr)
                                       ELSE
                                         CALL SAXPY(M,-t*aqoap,A(1,q),1,&
     &                                      A(1,p),1)
                                         CALL SAXPY(M,cs*sn*apoaq,A(1,p)&
     &                                      ,1,A(1,q),1)
                                         IF ( rsvec ) THEN
                                         CALL SAXPY(mvl,-t*aqoap,V(1,q),&
     &                                      1,V(1,p),1)
                                         CALL SAXPY(mvl,cs*sn*apoaq,    &
     &                                      V(1,p),1,V(1,q),1)
                                         ENDIF
                                         Work(p) = Work(p)*cs
                                         Work(q) = Work(q)/cs
                                       ENDIF
                                    ELSEIF ( Work(q)>=ONE ) THEN
                                       CALL SAXPY(M,t*apoaq,A(1,p),1,   &
     &                                    A(1,q),1)
                                       CALL SAXPY(M,-cs*sn*aqoap,A(1,q),&
     &                                    1,A(1,p),1)
                                       IF ( rsvec ) THEN
                                         CALL SAXPY(mvl,t*apoaq,V(1,p), &
     &                                      1,V(1,q),1)
                                         CALL SAXPY(mvl,-cs*sn*aqoap,   &
     &                                      V(1,q),1,V(1,p),1)
                                       ENDIF
                                       Work(p) = Work(p)/cs
                                       Work(q) = Work(q)*cs
                                    ELSEIF ( Work(p)>=Work(q) ) THEN
                                       CALL SAXPY(M,-t*aqoap,A(1,q),1,  &
     &                                    A(1,p),1)
                                       CALL SAXPY(M,cs*sn*apoaq,A(1,p), &
     &                                    1,A(1,q),1)
                                       Work(p) = Work(p)*cs
                                       Work(q) = Work(q)/cs
                                       IF ( rsvec ) THEN
                                         CALL SAXPY(mvl,-t*aqoap,V(1,q),&
     &                                      1,V(1,p),1)
                                         CALL SAXPY(mvl,cs*sn*apoaq,    &
     &                                      V(1,p),1,V(1,q),1)
                                       ENDIF
                                    ELSE
                                       CALL SAXPY(M,t*apoaq,A(1,p),1,   &
     &                                    A(1,q),1)
                                       CALL SAXPY(M,-cs*sn*aqoap,A(1,q),&
     &                                    1,A(1,p),1)
                                       Work(p) = Work(p)/cs
                                       Work(q) = Work(q)*cs
                                       IF ( rsvec ) THEN
                                         CALL SAXPY(mvl,t*apoaq,V(1,p), &
     &                                      1,V(1,q),1)
                                         CALL SAXPY(mvl,-cs*sn*aqoap,   &
     &                                      V(1,q),1,V(1,p),1)
                                       ENDIF
                                    ENDIF
                                 ENDIF
!
                              ELSEIF ( aapp>aaqq ) THEN
                                 CALL SCOPY(M,A(1,p),1,Work(N+1),1)
                                 CALL SLASCL('G',0,0,aapp,ONE,M,1,      &
     &                              Work(N+1),Lda,ierr)
                                 CALL SLASCL('G',0,0,aaqq,ONE,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 temp1 = -aapq*Work(p)/Work(q)
                                 CALL SAXPY(M,temp1,Work(N+1),1,A(1,q), &
     &                              1)
                                 CALL SLASCL('G',0,0,ONE,aaqq,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 Sva(q)                                 &
     &                              = aaqq*SQRT(MAX(ZERO,ONE-aapq*aapq))
                                 mxsinj = MAX(mxsinj,sfmin)
                              ELSE
                                 CALL SCOPY(M,A(1,q),1,Work(N+1),1)
                                 CALL SLASCL('G',0,0,aaqq,ONE,M,1,      &
     &                              Work(N+1),Lda,ierr)
                                 CALL SLASCL('G',0,0,aapp,ONE,M,1,A(1,p)&
     &                              ,Lda,ierr)
                                 temp1 = -aapq*Work(q)/Work(p)
                                 CALL SAXPY(M,temp1,Work(N+1),1,A(1,p), &
     &                              1)
                                 CALL SLASCL('G',0,0,ONE,aapp,M,1,A(1,p)&
     &                              ,Lda,ierr)
                                 Sva(p)                                 &
     &                              = aapp*SQRT(MAX(ZERO,ONE-aapq*aapq))
                                 mxsinj = MAX(mxsinj,sfmin)
                              ENDIF
!           END IF ROTOK THEN ... ELSE
!
!           In the case of cancellation in updating SVA(q)
!           .. recompute SVA(q)
                              IF ( (Sva(q)/aaqq)**2<=rooteps ) THEN
                                 IF ( (aaqq<rootbig) .AND.              &
     &                                (aaqq>rootsfmin) ) THEN
                                    Sva(q) = SNRM2(M,A(1,q),1)*Work(q)
                                 ELSE
                                    t = ZERO
                                    aaqq = ONE
                                    CALL SLASSQ(M,A(1,q),1,t,aaqq)
                                    Sva(q) = t*SQRT(aaqq)*Work(q)
                                 ENDIF
                              ENDIF
                              IF ( (aapp/aapp0)**2<=rooteps ) THEN
                                 IF ( (aapp<rootbig) .AND.              &
     &                                (aapp>rootsfmin) ) THEN
                                    aapp = SNRM2(M,A(1,p),1)*Work(p)
                                 ELSE
                                    t = ZERO
                                    aapp = ONE
                                    CALL SLASSQ(M,A(1,p),1,t,aapp)
                                    aapp = t*SQRT(aapp)*Work(p)
                                 ENDIF
                                 Sva(p) = aapp
                              ENDIF
!              end of OK rotation
                           ELSE
                              notrot = notrot + 1
![RTD]      SKIPPED  = SKIPPED  + 1
                              pskipped = pskipped + 1
                              ijblsk = ijblsk + 1
                           ENDIF
                        ELSE
                           notrot = notrot + 1
                           pskipped = pskipped + 1
                           ijblsk = ijblsk + 1
                        ENDIF
!
                        IF ( (i<=swband) .AND. (ijblsk>=blskip) ) THEN
                           Sva(p) = aapp
                           notrot = 0
                           GOTO 20
                        ENDIF
                        IF ( (i<=swband) .AND. (pskipped>rowskip) ) THEN
                           aapp = -aapp
                           notrot = 0
                           EXIT
                        ENDIF
!
                     ENDDO
!        end of the q-loop
!
                     Sva(p) = aapp
!
                  ELSE
!
                     IF ( aapp==ZERO ) notrot = notrot +                &
     &                    MIN(jgl+kbl-1,N) - jgl + 1
                     IF ( aapp<ZERO ) notrot = 0
!
                  ENDIF
!
               ENDDO
!     end of the p-loop
            ENDDO
!     end of the jbc-loop
!2011 bailed out of the jbc-loop
 20         DO p = igl , MIN(igl+kbl-1,N)
               Sva(p) = ABS(Sva(p))
            ENDDO
!**
         ENDDO
!2000 :: end of the ibr-loop
!
!     .. update SVA(N)
         IF ( (Sva(N)<rootbig) .AND. (Sva(N)>rootsfmin) ) THEN
            Sva(N) = SNRM2(M,A(1,N),1)*Work(N)
         ELSE
            t = ZERO
            aapp = ONE
            CALL SLASSQ(M,A(1,N),1,t,aapp)
            Sva(N) = t*SQRT(aapp)*Work(N)
         ENDIF
!
!     Additional steering devices
!
         IF ( (i<swband) .AND. ((mxaapq<=roottol) .OR. (iswrot<=N)) )   &
     &        swband = i
!
         IF ( (i>swband+1) .AND. (mxaapq<SQRT(FLOAT(N))*tol) .AND.      &
     &        (FLOAT(N)*mxaapq*mxsinj<tol) ) GOTO 100
!
         IF ( notrot>=emptsw ) GOTO 100
!
      ENDDO
!     end i=1:NSWEEP loop
!
! #:( Reaching this point means that the procedure has not converged.
      Info = NSWEEP - 1
      GOTO 200
!
! #:) Reaching this point means numerical convergence after the i-th
!     sweep.
!
 100  Info = 0
! #:) INFO = 0 confirms successful iterations.
!
!     Sort the singular values and find how many are above
!     the underflow threshold.
!
 200  n2 = 0
      n4 = 0
      DO p = 1 , N - 1
         q = ISAMAX(N-p+1,Sva(p),1) + p - 1
         IF ( p/=q ) THEN
            temp1 = Sva(p)
            Sva(p) = Sva(q)
            Sva(q) = temp1
            temp1 = Work(p)
            Work(p) = Work(q)
            Work(q) = temp1
            CALL SSWAP(M,A(1,p),1,A(1,q),1)
            IF ( rsvec ) CALL SSWAP(mvl,V(1,p),1,V(1,q),1)
         ENDIF
         IF ( Sva(p)/=ZERO ) THEN
            n4 = n4 + 1
            IF ( Sva(p)*skl>sfmin ) n2 = n2 + 1
         ENDIF
      ENDDO
      IF ( Sva(N)/=ZERO ) THEN
         n4 = n4 + 1
         IF ( Sva(N)*skl>sfmin ) n2 = n2 + 1
      ENDIF
!
!     Normalize the left singular vectors.
!
      IF ( lsvec .OR. uctol ) THEN
         DO p = 1 , n2
            CALL SSCAL(M,Work(p)/Sva(p),A(1,p),1)
         ENDDO
      ENDIF
!
!     Scale the product of Jacobi rotations (assemble the fast rotations).
!
      IF ( rsvec ) THEN
         IF ( applv ) THEN
            DO p = 1 , N
               CALL SSCAL(mvl,Work(p),V(1,p),1)
            ENDDO
         ELSE
            DO p = 1 , N
               temp1 = ONE/SNRM2(mvl,V(1,p),1)
               CALL SSCAL(mvl,temp1,V(1,p),1)
            ENDDO
         ENDIF
      ENDIF
!
!     Undo scaling, if necessary (and possible).
      IF ( ((skl>ONE) .AND. (Sva(1)<(big/skl))) .OR.                    &
     &     ((skl<ONE) .AND. (Sva(MAX(n2,1))>(sfmin/skl))) ) THEN
         DO p = 1 , N
            Sva(p) = skl*Sva(p)
         ENDDO
         skl = ONE
      ENDIF
!
      Work(1) = skl
!     The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE
!     then some of the singular values may overflow or underflow and
!     the spectrum is given in this factored representation.
!
      Work(2) = FLOAT(n4)
!     N4 is the number of computed nonzero singular values of A.
!
      Work(3) = FLOAT(n2)
!     N2 is the number of singular values of A greater than SFMIN.
!     If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers
!     that may carry some information.
!
      Work(4) = FLOAT(i)
!     i is the index of the last sweep before declaring convergence.
!
      Work(5) = mxaapq
!     MXAAPQ is the largest absolute value of scaled pivots in the
!     last sweep
!
      Work(6) = mxsinj
!     MXSINJ is the largest absolute value of the sines of Jacobi angles
!     in the last sweep
!
!     ..
!     .. END OF SGESVJ
!     ..
      END SUBROUTINE SGESVJ
