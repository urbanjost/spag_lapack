!*==zgesvj.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> ZGESVJ </b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGESVJ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgesvj.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgesvj.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgesvj.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V,
!                          LDV, CWORK, LWORK, RWORK, LRWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDV, LWORK, LRWORK, M, MV, N
!       CHARACTER*1        JOBA, JOBU, JOBV
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ),  V( LDV, * ), CWORK( LWORK )
!       DOUBLE PRECISION   RWORK( LRWORK ),  SVA( N )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGESVJ computes the singular value decomposition (SVD) of a complex
!> M-by-N matrix A, where M >= N. The SVD of A is written as
!>                                    [++]   [xx]   [x0]   [xx]
!>              A = U * SIGMA * V^*,  [++] = [xx] * [ox] * [xx]
!>                                    [++]   [xx]
!> where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal
!> matrix, and V is an N-by-N unitary matrix. The diagonal elements
!> of SIGMA are the singular values of A. The columns of U and V are the
!> left and the right singular vectors of A, respectively.
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
!>          = 'U' or 'F': The left singular vectors corresponding to the nonzero
!>                 singular values are computed and returned in the leading
!>                 columns of A. See more details in the description of A.
!>                 The default numerical orthogonality threshold is set to
!>                 approximately TOL=CTOL*EPS, CTOL=SQRT(M), EPS=DLAMCH('E').
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
!>          = 'V' or 'J': the matrix V is computed and returned in the array V
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
!>          The number of rows of the input matrix A. 1/DLAMCH('E') > M >= 0.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit,
!>          If JOBU = 'U' .OR. JOBU = 'C':
!>                 If INFO = 0 :
!>                 RANKA orthonormal columns of U are returned in the
!>                 leading RANKA columns of the array A. Here RANKA <= N
!>                 is the number of computed singular values of A that are
!>                 above the underflow threshold DLAMCH('S'). The singular
!>                 vectors corresponding to underflowed or zero singular
!>                 values are not computed. The value of RANKA is returned
!>                 in the array RWORK as RANKA=NINT(RWORK(2)). Also see the
!>                 descriptions of SVA and RWORK. The computed columns of U
!>                 are mutually numerically orthogonal up to approximately
!>                 TOL=SQRT(M)*EPS (default); or TOL=CTOL*EPS (JOBU = 'C'),
!>                 see the description of JOBU.
!>                 If INFO > 0,
!>                 the procedure ZGESVJ did not converge in the given number
!>                 of iterations (sweeps). In that case, the computed
!>                 columns of U may not be orthogonal up to TOL. The output
!>                 U (stored in A), SIGMA (given by the computed singular
!>                 values in SVA(1:N)) and V is still a decomposition of the
!>                 input matrix A in the sense that the residual
!>                 || A - SCALE * U * SIGMA * V^* ||_2 / ||A||_2 is small.
!>          If JOBU = 'N':
!>                 If INFO = 0 :
!>                 Note that the left singular vectors are 'for free' in the
!>                 one-sided Jacobi SVD algorithm. However, if only the
!>                 singular values are needed, the level of numerical
!>                 orthogonality of U is not an issue and iterations are
!>                 stopped when the columns of the iterated matrix are
!>                 numerically orthogonal up to approximately M*EPS. Thus,
!>                 on exit, A contains the columns of U scaled with the
!>                 corresponding singular values.
!>                 If INFO > 0:
!>                 the procedure ZGESVJ did not converge in the given number
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
!>          SVA is DOUBLE PRECISION array, dimension (N)
!>          On exit,
!>          If INFO = 0 :
!>          depending on the value SCALE = RWORK(1), we have:
!>                 If SCALE = ONE:
!>                 SVA(1:N) contains the computed singular values of A.
!>                 During the computation SVA contains the Euclidean column
!>                 norms of the iterated matrices in the array A.
!>                 If SCALE .NE. ONE:
!>                 The singular values of A are SCALE*SVA(1:N), and this
!>                 factored representation is due to the fact that some of the
!>                 singular values of A might underflow or overflow.
!>
!>          If INFO > 0:
!>          the procedure ZGESVJ did not converge in the given number of
!>          iterations (sweeps) and SCALE*SVA(1:N) may not be accurate.
!> \endverbatim
!>
!> \param[in] MV
!> \verbatim
!>          MV is INTEGER
!>          If JOBV = 'A', then the product of Jacobi rotations in ZGESVJ
!>          is applied to the first MV rows of V. See the description of JOBV.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (LDV,N)
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
!> \param[in,out] CWORK
!> \verbatim
!>          CWORK is COMPLEX*16 array, dimension (max(1,LWORK))
!>          Used as workspace.
!>          If on entry LWORK = -1, then a workspace query is assumed and
!>          no computation is done; CWORK(1) is set to the minial (and optimal)
!>          length of CWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER.
!>          Length of CWORK, LWORK >= M+N.
!> \endverbatim
!>
!> \param[in,out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (max(6,LRWORK))
!>          On entry,
!>          If JOBU = 'C' :
!>          RWORK(1) = CTOL, where CTOL defines the threshold for convergence.
!>                    The process stops if all columns of A are mutually
!>                    orthogonal up to CTOL*EPS, EPS=DLAMCH('E').
!>                    It is required that CTOL >= ONE, i.e. it is not
!>                    allowed to force the routine to obtain orthogonality
!>                    below EPSILON.
!>          On exit,
!>          RWORK(1) = SCALE is the scaling factor such that SCALE*SVA(1:N)
!>                    are the computed singular values of A.
!>                    (See description of SVA().)
!>          RWORK(2) = NINT(RWORK(2)) is the number of the computed nonzero
!>                    singular values.
!>          RWORK(3) = NINT(RWORK(3)) is the number of the computed singular
!>                    values that are larger than the underflow threshold.
!>          RWORK(4) = NINT(RWORK(4)) is the number of sweeps of Jacobi
!>                    rotations needed for numerical convergence.
!>          RWORK(5) = max_{i.NE.j} |COS(A(:,i),A(:,j))| in the last sweep.
!>                    This is useful information in cases when ZGESVJ did
!>                    not converge, as it can be used to estimate whether
!>                    the output is still useful and for post festum analysis.
!>          RWORK(6) = the largest absolute value over all sines of the
!>                    Jacobi rotation angles in the last sweep. It can be
!>                    useful for a post festum analysis.
!>         If on entry LRWORK = -1, then a workspace query is assumed and
!>         no computation is done; RWORK(1) is set to the minial (and optimal)
!>         length of RWORK.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>         LRWORK is INTEGER
!>         Length of RWORK, LRWORK >= MAX(6,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, then the i-th argument had an illegal value
!>          > 0:  ZGESVJ did not converge in the maximal allowed number
!>                (NSWEEP=30) of sweeps. The output may still be useful.
!>                See the description of RWORK.
!> \endverbatim
!>
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date June 2016
!
!> \ingroup complex16GEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!> The orthogonal N-by-N matrix V is obtained as a product of Jacobi plane
!> rotations. In the case of underflow of the tangent of the Jacobi angle, a
!> modified Jacobi transformation of Drmac [3] is used. Pivot strategy uses
!> column interchanges of de Rijk [1]. The relative accuracy of the computed
!> singular values and the accuracy of the computed singular vectors (in
!> angle metric) is as guaranteed by the theory of Demmel and Veselic [2].
!> The condition number that determines the accuracy in the full rank case
!> is essentially min_{D=diag} kappa(A*D), where kappa(.) is the
!> spectral condition number. The best performance of this Jacobi SVD
!> procedure is achieved if used in an  accelerated version of Drmac and
!> Veselic [4,5], and it is the kernel routine in the SIGMA library [6].
!> Some tuning parameters (marked with [TP]) are available for the
!> implementer.
!> The computational range for the nonzero singular values is the  machine
!> number interval ( UNDERFLOW , OVERFLOW ). In extreme cases, even
!> denormalized singular values can be computed with the corresponding
!> gradual loss of accurate digits.
!> \endverbatim
!
!> \par Contributor:
!  ==================
!>
!> \verbatim
!>
!>  ============
!>
!>  Zlatko Drmac (Zagreb, Croatia)
!>
!> \endverbatim
!
!> \par References:
!  ================
!>
!> \verbatim
!>
!> [1] P. P. M. De Rijk: A one-sided Jacobi algorithm for computing the
!>    singular value decomposition on a vector computer.
!>    SIAM J. Sci. Stat. Comp., Vol. 10 (1998), pp. 359-371.
!> [2] J. Demmel and K. Veselic: Jacobi method is more accurate than QR.
!> [3] Z. Drmac: Implementation of Jacobi rotations for accurate singular
!>    value computation in floating point arithmetic.
!>    SIAM J. Sci. Comp., Vol. 18 (1997), pp. 1200-1222.
!> [4] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I.
!>    SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342.
!>    LAPACK Working note 169.
!> [5] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II.
!>    SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362.
!>    LAPACK Working note 170.
!> [6] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV,
!>    QSVD, (H,K)-SVD computations.
!>    Department of Mathematics, University of Zagreb, 2008, 2015.
!> \endverbatim
!
!> \par Bugs, examples and comments:
!  =================================
!>
!> \verbatim
!>  ===========================
!>  Please report all bugs and send interesting test examples and comments to
!>  drmac@math.hr. Thank you.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGESVJ(Joba,Jobu,Jobv,M,N,A,Lda,Sva,Mv,V,Ldv,Cwork,    &
     &                  Lwork,Rwork,Lrwork,Info)
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
      IMPLICIT NONE
!*--ZGESVJ361
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldv , Lwork , Lrwork , M , Mv , N
      CHARACTER*1 Joba , Jobu , Jobv
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , V(Ldv,*) , Cwork(Lwork)
      DOUBLE PRECISION Rwork(Lrwork) , Sva(N)
!     ..
!
!  =====================================================================
!
!     .. Local Parameters ..
      DOUBLE PRECISION ZERO , HALF , ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0))
      INTEGER NSWEEP
      PARAMETER (NSWEEP=30)
!     ..
!     .. Local Scalars ..
      COMPLEX*16 aapq , ompq
      DOUBLE PRECISION aapp , aapp0 , aapq1 , aaqq , apoaq , aqoap ,    &
     &                 big , bigtheta , cs , ctol , epsln , mxaapq ,    &
     &                 mxsinj , rootbig , rooteps , rootsfmin ,         &
     &                 roottol , skl , sfmin , small , sn , t , temp1 , &
     &                 theta , thsign , tol
      INTEGER blskip , emptsw , i , ibr , ierr , igl , ijblsk , ir1 ,   &
     &        iswrot , jbc , jgl , kbl , lkahead , mvl , n2 , n34 , n4 ,&
     &        nbl , notrot , p , pskipped , q , rowskip , swband
      LOGICAL applv , goscale , lower , lquery , lsvec , noscale ,      &
     &        rotok , rsvec , uctol , upper
!     ..
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , CONJG , DBLE , SIGN , SQRT
!     ..
!     .. External Functions ..
!     ..
!     from BLAS
      DOUBLE PRECISION DZNRM2
      COMPLEX*16 ZDOTC
      EXTERNAL ZDOTC , DZNRM2
      INTEGER IDAMAX
      EXTERNAL IDAMAX
!     from LAPACK
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!     ..
!     from BLAS
      EXTERNAL ZCOPY , ZROT , ZDSCAL , ZSWAP , ZAXPY
!     from LAPACK
      EXTERNAL DLASCL , ZLASCL , ZLASET , ZLASSQ , XERBLA
      EXTERNAL ZGSVJ0 , ZGSVJ1
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      lsvec = LSAME(Jobu,'U') .OR. LSAME(Jobu,'F')
      uctol = LSAME(Jobu,'C')
      rsvec = LSAME(Jobv,'V') .OR. LSAME(Jobv,'J')
      applv = LSAME(Jobv,'A')
      upper = LSAME(Joba,'U')
      lower = LSAME(Joba,'L')
!
      lquery = (Lwork==-1) .OR. (Lrwork==-1)
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
      ELSEIF ( uctol .AND. (Rwork(1)<=ONE) ) THEN
         Info = -12
      ELSEIF ( (Lwork<(M+N)) .AND. (.NOT.lquery) ) THEN
         Info = -13
      ELSEIF ( (Lrwork<MAX(N,6)) .AND. (.NOT.lquery) ) THEN
         Info = -15
      ELSE
         Info = 0
      ENDIF
!
!     #:(
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGESVJ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         Cwork(1) = M + N
         Rwork(1) = MAX(N,6)
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
!     max_{i<>j}|A(:,i)^* * A(:,j)| / (||A(:,i)||*||A(:,j)||) < CTOL*EPS
!
!     where EPS is the round-off and CTOL is defined as follows:
!
      IF ( uctol ) THEN
!        ... user controlled
         ctol = Rwork(1)
!        ... default
      ELSEIF ( lsvec .OR. rsvec .OR. applv ) THEN
         ctol = SQRT(DBLE(M))
      ELSE
         ctol = DBLE(M)
      ENDIF
!     ... and the machine dependent parameters are
![!]  (Make sure that SLAMCH() works properly on the target machine.)
!
      epsln = DLAMCH('Epsilon')
      rooteps = SQRT(epsln)
      sfmin = DLAMCH('SafeMinimum')
      rootsfmin = SQRT(sfmin)
      small = sfmin/epsln
      big = DLAMCH('Overflow')
!     BIG         = ONE    / SFMIN
      rootbig = ONE/rootsfmin
!      LARGE = BIG / SQRT( DBLE( M*N ) )
      bigtheta = ONE/rooteps
!
      tol = ctol*epsln
      roottol = SQRT(tol)
!
      IF ( DBLE(M)*epsln>=ONE ) THEN
         Info = -4
         CALL XERBLA('ZGESVJ',-Info)
         RETURN
      ENDIF
!
!     Initialize the right singular vector matrix.
!
      IF ( rsvec ) THEN
         mvl = N
         CALL ZLASET('A',mvl,N,CZERO,CONE,V,Ldv)
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
      skl = ONE/SQRT(DBLE(M)*DBLE(N))
      noscale = .TRUE.
      goscale = .TRUE.
!
      IF ( lower ) THEN
!        the input matrix is M-by-N lower triangular (trapezoidal)
         DO p = 1 , N
            aapp = ZERO
            aaqq = ONE
            CALL ZLASSQ(M-p+1,A(p,p),1,aapp,aaqq)
            IF ( aapp>big ) THEN
               Info = -6
               CALL XERBLA('ZGESVJ',-Info)
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
            CALL ZLASSQ(p,A(1,p),1,aapp,aaqq)
            IF ( aapp>big ) THEN
               Info = -6
               CALL XERBLA('ZGESVJ',-Info)
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
            CALL ZLASSQ(M,A(1,p),1,aapp,aaqq)
            IF ( aapp>big ) THEN
               Info = -6
               CALL XERBLA('ZGESVJ',-Info)
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
         IF ( lsvec ) CALL ZLASET('G',M,N,CZERO,CONE,A,Lda)
         Rwork(1) = ONE
         Rwork(2) = ZERO
         Rwork(3) = ZERO
         Rwork(4) = ZERO
         Rwork(5) = ZERO
         Rwork(6) = ZERO
         RETURN
      ENDIF
!
! #:) Quick return for one-column matrix
!
      IF ( N==1 ) THEN
         IF ( lsvec ) CALL ZLASCL('G',0,0,Sva(1),skl,M,1,A(1,1),Lda,    &
     &                            ierr)
         Rwork(1) = ONE/skl
         IF ( Sva(1)>=sfmin ) THEN
            Rwork(2) = ONE
         ELSE
            Rwork(2) = ZERO
         ENDIF
         Rwork(3) = ZERO
         Rwork(4) = ZERO
         Rwork(5) = ZERO
         Rwork(6) = ZERO
         RETURN
      ENDIF
!
!     Protect small singular values from underflow, and try to
!     avoid underflows/overflows in computing Jacobi rotations.
!
      sn = SQRT(sfmin/epsln)
      temp1 = SQRT(big/DBLE(N))
      IF ( (aapp<=sn) .OR. (aaqq>=temp1) .OR.                           &
     &     ((sn<=aaqq) .AND. (aapp<=temp1)) ) THEN
         temp1 = MIN(big,temp1/aapp)
!         AAQQ  = AAQQ*TEMP1
!         AAPP  = AAPP*TEMP1
      ELSEIF ( (aaqq<=sn) .AND. (aapp<=temp1) ) THEN
         temp1 = MIN(sn/aaqq,big/(aapp*SQRT(DBLE(N))))
!         AAQQ  = AAQQ*TEMP1
!         AAPP  = AAPP*TEMP1
      ELSEIF ( (aaqq>=sn) .AND. (aapp>=temp1) ) THEN
         temp1 = MAX(sn/aaqq,temp1/aapp)
!         AAQQ  = AAQQ*TEMP1
!         AAPP  = AAPP*TEMP1
      ELSEIF ( (aaqq<=sn) .AND. (aapp>=temp1) ) THEN
         temp1 = MIN(sn/aaqq,big/(SQRT(DBLE(N))*aapp))
!         AAQQ  = AAQQ*TEMP1
!         AAPP  = AAPP*TEMP1
      ELSE
         temp1 = ONE
      ENDIF
!
!     Scale, if necessary
!
      IF ( temp1/=ONE ) CALL DLASCL('G',0,0,ONE,temp1,N,1,Sva,N,ierr)
      skl = temp1*skl
      IF ( skl/=ONE ) THEN
         CALL ZLASCL(Joba,0,0,ONE,skl,M,N,A,Lda,ierr)
         skl = ONE/skl
      ENDIF
!
!     Row-cyclic Jacobi SVD algorithm with column pivoting
!
      emptsw = (N*(N-1))/2
      notrot = 0
 
      DO q = 1 , N
         Cwork(q) = CONE
      ENDDO
!
!
!
      swband = 3
![TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective
!     if ZGESVJ is used as a computational routine in the preconditioned
!     Jacobi SVD algorithm ZGEJSV. For sweeps i=1:SWBAND the procedure
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
            CALL ZGSVJ0(Jobv,M-n34,N-n34,A(n34+1,n34+1),Lda,Cwork(n34+1)&
     &                  ,Sva(n34+1),mvl,V(n34*q+1,n34+1),Ldv,epsln,     &
     &                  sfmin,tol,2,Cwork(N+1),Lwork-N,ierr)
 
            CALL ZGSVJ0(Jobv,M-n2,n34-n2,A(n2+1,n2+1),Lda,Cwork(n2+1),  &
     &                  Sva(n2+1),mvl,V(n2*q+1,n2+1),Ldv,epsln,sfmin,   &
     &                  tol,2,Cwork(N+1),Lwork-N,ierr)
 
            CALL ZGSVJ1(Jobv,M-n2,N-n2,n4,A(n2+1,n2+1),Lda,Cwork(n2+1), &
     &                  Sva(n2+1),mvl,V(n2*q+1,n2+1),Ldv,epsln,sfmin,   &
     &                  tol,1,Cwork(N+1),Lwork-N,ierr)
 
            CALL ZGSVJ0(Jobv,M-n4,n2-n4,A(n4+1,n4+1),Lda,Cwork(n4+1),   &
     &                  Sva(n4+1),mvl,V(n4*q+1,n4+1),Ldv,epsln,sfmin,   &
     &                  tol,1,Cwork(N+1),Lwork-N,ierr)
!
            CALL ZGSVJ0(Jobv,M,n4,A,Lda,Cwork,Sva,mvl,V,Ldv,epsln,sfmin,&
     &                  tol,1,Cwork(N+1),Lwork-N,ierr)
!
            CALL ZGSVJ1(Jobv,M,n2,n4,A,Lda,Cwork,Sva,mvl,V,Ldv,epsln,   &
     &                  sfmin,tol,1,Cwork(N+1),Lwork-N,ierr)
!
!
         ELSEIF ( upper ) THEN
!
!
            CALL ZGSVJ0(Jobv,n4,n4,A,Lda,Cwork,Sva,mvl,V,Ldv,epsln,     &
     &                  sfmin,tol,2,Cwork(N+1),Lwork-N,ierr)
!
            CALL ZGSVJ0(Jobv,n2,n4,A(1,n4+1),Lda,Cwork(n4+1),Sva(n4+1), &
     &                  mvl,V(n4*q+1,n4+1),Ldv,epsln,sfmin,tol,1,       &
     &                  Cwork(N+1),Lwork-N,ierr)
!
            CALL ZGSVJ1(Jobv,n2,n2,n4,A,Lda,Cwork,Sva,mvl,V,Ldv,epsln,  &
     &                  sfmin,tol,1,Cwork(N+1),Lwork-N,ierr)
!
            CALL ZGSVJ0(Jobv,n2+n4,n4,A(1,n2+1),Lda,Cwork(n2+1),        &
     &                  Sva(n2+1),mvl,V(n2*q+1,n2+1),Ldv,epsln,sfmin,   &
     &                  tol,1,Cwork(N+1),Lwork-N,ierr)
 
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
                  q = IDAMAX(N-p+1,Sva(p),1) + p - 1
                  IF ( p/=q ) THEN
                     CALL ZSWAP(M,A(1,p),1,A(1,q),1)
                     IF ( rsvec ) CALL ZSWAP(mvl,V(1,p),1,V(1,q),1)
                     temp1 = Sva(p)
                     Sva(p) = Sva(q)
                     Sva(q) = temp1
                     aapq = Cwork(p)
                     Cwork(p) = Cwork(q)
                     Cwork(q) = aapq
                  ENDIF
!
                  IF ( ir1==0 ) THEN
!
!        Column norms are periodically updated by explicit
!        norm computation.
![!]     Caveat:
!        Unfortunately, some BLAS implementations compute DZNRM2(M,A(1,p),1)
!        as SQRT(S=CDOTC(M,A(1,p),1,A(1,p),1)), which may cause the result to
!        overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to
!        underflow for ||A(:,p)||_2 < SQRT(underflow_threshold).
!        Hence, DZNRM2 cannot be trusted, not even in the case when
!        the true norm is far from the under(over)flow boundaries.
!        If properly implemented SCNRM2 is available, the IF-THEN-ELSE-END IF
!        below should be replaced with "AAPP = DZNRM2( M, A(1,p), 1 )".
!
                     IF ( (Sva(p)<rootbig) .AND. (Sva(p)>rootsfmin) )   &
     &                    THEN
                        Sva(p) = DZNRM2(M,A(1,p),1)
                     ELSE
                        temp1 = ZERO
                        aapp = ONE
                        CALL ZLASSQ(M,A(1,p),1,temp1,aapp)
                        Sva(p) = temp1*SQRT(aapp)
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
                                 aapq = (ZDOTC(M,A(1,p),1,A(1,q),1)     &
     &                                  /aaqq)/aapp
                              ELSE
                                 CALL ZCOPY(M,A(1,p),1,Cwork(N+1),1)
                                 CALL ZLASCL('G',0,0,aapp,ONE,M,1,      &
     &                              Cwork(N+1),Lda,ierr)
                                 aapq = ZDOTC(M,Cwork(N+1),1,A(1,q),1)  &
     &                                  /aaqq
                              ENDIF
                           ELSE
                              rotok = aapp<=(aaqq/small)
                              IF ( aapp>(small/aaqq) ) THEN
                                 aapq = (ZDOTC(M,A(1,p),1,A(1,q),1)     &
     &                                  /aapp)/aaqq
                              ELSE
                                 CALL ZCOPY(M,A(1,q),1,Cwork(N+1),1)
                                 CALL ZLASCL('G',0,0,aaqq,ONE,M,1,      &
     &                              Cwork(N+1),Lda,ierr)
                                 aapq = ZDOTC(M,A(1,p),1,Cwork(N+1),1)  &
     &                                  /aapp
                              ENDIF
                           ENDIF
!
 
!                           AAPQ = AAPQ * CONJG( CWORK(p) ) * CWORK(q)
                           aapq1 = -ABS(aapq)
                           mxaapq = MAX(mxaapq,-aapq1)
!
!        TO rotate or NOT to rotate, THAT is the question ...
!
                           IF ( ABS(aapq1)>tol ) THEN
                              ompq = aapq/ABS(aapq)
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
                                 theta = -HALF*ABS(aqoap-apoaq)/aapq1
!
                                 IF ( ABS(theta)>bigtheta ) THEN
!
                                    t = HALF/theta
                                    cs = ONE
 
                                    CALL ZROT(M,A(1,p),1,A(1,q),1,cs,   &
     &                                 CONJG(ompq)*t)
                                    IF ( rsvec )                        &
     &                                 CALL ZROT(mvl,V(1,p),1,V(1,q),1, &
     &                                 cs,CONJG(ompq)*t)
 
                                    Sva(q)                              &
     &                                 = aaqq*SQRT(MAX(ZERO,ONE+t*apoaq*&
     &                                 aapq1))
                                    aapp = aapp*SQRT                    &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq1))
                                    mxsinj = MAX(mxsinj,ABS(t))
!
                                 ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                    thsign = -SIGN(ONE,aapq1)
                                    t = ONE/                            &
     &                                  (theta+thsign*SQRT(ONE+theta*   &
     &                                  theta))
                                    cs = SQRT(ONE/(ONE+t*t))
                                    sn = t*cs
!
                                    mxsinj = MAX(mxsinj,ABS(sn))
                                    Sva(q)                              &
     &                                 = aaqq*SQRT(MAX(ZERO,ONE+t*apoaq*&
     &                                 aapq1))
                                    aapp = aapp*SQRT                    &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq1))
!
                                    CALL ZROT(M,A(1,p),1,A(1,q),1,cs,   &
     &                                 CONJG(ompq)*sn)
                                    IF ( rsvec )                        &
     &                                 CALL ZROT(mvl,V(1,p),1,V(1,q),1, &
     &                                 cs,CONJG(ompq)*sn)
                                 ENDIF
                                 Cwork(p) = -Cwork(q)*ompq
!
                              ELSE
!              .. have to use modified Gram-Schmidt like transformation
                                 CALL ZCOPY(M,A(1,p),1,Cwork(N+1),1)
                                 CALL ZLASCL('G',0,0,aapp,ONE,M,1,      &
     &                              Cwork(N+1),Lda,ierr)
                                 CALL ZLASCL('G',0,0,aaqq,ONE,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 CALL ZAXPY(M,-aapq,Cwork(N+1),1,A(1,q),&
     &                              1)
                                 CALL ZLASCL('G',0,0,ONE,aaqq,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 Sva(q) = aaqq*SQRT(MAX(ZERO,ONE-aapq1* &
     &                              aapq1))
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
                                    Sva(q) = DZNRM2(M,A(1,q),1)
                                 ELSE
                                    t = ZERO
                                    aaqq = ONE
                                    CALL ZLASSQ(M,A(1,q),1,t,aaqq)
                                    Sva(q) = t*SQRT(aaqq)
                                 ENDIF
                              ENDIF
                              IF ( (aapp/aapp0)<=rooteps ) THEN
                                 IF ( (aapp<rootbig) .AND.              &
     &                                (aapp>rootsfmin) ) THEN
                                    aapp = DZNRM2(M,A(1,p),1)
                                 ELSE
                                    t = ZERO
                                    aapp = ONE
                                    CALL ZLASSQ(M,A(1,p),1,t,aapp)
                                    aapp = t*SQRT(aapp)
                                 ENDIF
                                 Sva(p) = aapp
                              ENDIF
!
                           ELSE
!                             A(:,p) and A(:,q) already numerically orthogonal
                              IF ( ir1==0 ) notrot = notrot + 1
![RTD]      SKIPPED  = SKIPPED + 1
                              pskipped = pskipped + 1
                           ENDIF
                        ELSE
!                          A(:,q) is zero column
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
                                 aapq = (ZDOTC(M,A(1,p),1,A(1,q),1)     &
     &                                  /aaqq)/aapp
                              ELSE
                                 CALL ZCOPY(M,A(1,p),1,Cwork(N+1),1)
                                 CALL ZLASCL('G',0,0,aapp,ONE,M,1,      &
     &                              Cwork(N+1),Lda,ierr)
                                 aapq = ZDOTC(M,Cwork(N+1),1,A(1,q),1)  &
     &                                  /aaqq
                              ENDIF
                           ELSE
                              IF ( aapp>=aaqq ) THEN
                                 rotok = aapp<=(aaqq/small)
                              ELSE
                                 rotok = aaqq<=(aapp/small)
                              ENDIF
                              IF ( aapp>(small/aaqq) ) THEN
                                 aapq = (ZDOTC(M,A(1,p),1,A(1,q),1)     &
     &                                  /MAX(aaqq,aapp))/MIN(aaqq,aapp)
                              ELSE
                                 CALL ZCOPY(M,A(1,q),1,Cwork(N+1),1)
                                 CALL ZLASCL('G',0,0,aaqq,ONE,M,1,      &
     &                              Cwork(N+1),Lda,ierr)
                                 aapq = ZDOTC(M,A(1,p),1,Cwork(N+1),1)  &
     &                                  /aapp
                              ENDIF
                           ENDIF
!
 
!                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q)
                           aapq1 = -ABS(aapq)
                           mxaapq = MAX(mxaapq,-aapq1)
!
!        TO rotate or NOT to rotate, THAT is the question ...
!
                           IF ( ABS(aapq1)>tol ) THEN
                              ompq = aapq/ABS(aapq)
                              notrot = 0
![RTD]      ROTATED  = ROTATED + 1
                              pskipped = 0
                              iswrot = iswrot + 1
!
                              IF ( rotok ) THEN
!
                                 aqoap = aaqq/aapp
                                 apoaq = aapp/aaqq
                                 theta = -HALF*ABS(aqoap-apoaq)/aapq1
                                 IF ( aaqq>aapp0 ) theta = -theta
!
                                 IF ( ABS(theta)>bigtheta ) THEN
                                    t = HALF/theta
                                    cs = ONE
                                    CALL ZROT(M,A(1,p),1,A(1,q),1,cs,   &
     &                                 CONJG(ompq)*t)
                                    IF ( rsvec )                        &
     &                                 CALL ZROT(mvl,V(1,p),1,V(1,q),1, &
     &                                 cs,CONJG(ompq)*t)
                                    Sva(q)                              &
     &                                 = aaqq*SQRT(MAX(ZERO,ONE+t*apoaq*&
     &                                 aapq1))
                                    aapp = aapp*SQRT                    &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq1))
                                    mxsinj = MAX(mxsinj,ABS(t))
                                 ELSE
!
!                 .. choose correct signum for THETA and rotate
!
                                    thsign = -SIGN(ONE,aapq1)
                                    IF ( aaqq>aapp0 ) thsign = -thsign
                                    t = ONE/                            &
     &                                  (theta+thsign*SQRT(ONE+theta*   &
     &                                  theta))
                                    cs = SQRT(ONE/(ONE+t*t))
                                    sn = t*cs
                                    mxsinj = MAX(mxsinj,ABS(sn))
                                    Sva(q)                              &
     &                                 = aaqq*SQRT(MAX(ZERO,ONE+t*apoaq*&
     &                                 aapq1))
                                    aapp = aapp*SQRT                    &
     &                                 (MAX(ZERO,ONE-t*aqoap*aapq1))
!
                                    CALL ZROT(M,A(1,p),1,A(1,q),1,cs,   &
     &                                 CONJG(ompq)*sn)
                                    IF ( rsvec )                        &
     &                                 CALL ZROT(mvl,V(1,p),1,V(1,q),1, &
     &                                 cs,CONJG(ompq)*sn)
                                 ENDIF
                                 Cwork(p) = -Cwork(q)*ompq
!
!              .. have to use modified Gram-Schmidt like transformation
                              ELSEIF ( aapp>aaqq ) THEN
                                 CALL ZCOPY(M,A(1,p),1,Cwork(N+1),1)
                                 CALL ZLASCL('G',0,0,aapp,ONE,M,1,      &
     &                              Cwork(N+1),Lda,ierr)
                                 CALL ZLASCL('G',0,0,aaqq,ONE,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 CALL ZAXPY(M,-aapq,Cwork(N+1),1,A(1,q),&
     &                              1)
                                 CALL ZLASCL('G',0,0,ONE,aaqq,M,1,A(1,q)&
     &                              ,Lda,ierr)
                                 Sva(q) = aaqq*SQRT(MAX(ZERO,ONE-aapq1* &
     &                              aapq1))
                                 mxsinj = MAX(mxsinj,sfmin)
                              ELSE
                                 CALL ZCOPY(M,A(1,q),1,Cwork(N+1),1)
                                 CALL ZLASCL('G',0,0,aaqq,ONE,M,1,      &
     &                              Cwork(N+1),Lda,ierr)
                                 CALL ZLASCL('G',0,0,aapp,ONE,M,1,A(1,p)&
     &                              ,Lda,ierr)
                                 CALL ZAXPY(M,-CONJG(aapq),Cwork(N+1),1,&
     &                              A(1,p),1)
                                 CALL ZLASCL('G',0,0,ONE,aapp,M,1,A(1,p)&
     &                              ,Lda,ierr)
                                 Sva(p) = aapp*SQRT(MAX(ZERO,ONE-aapq1* &
     &                              aapq1))
                                 mxsinj = MAX(mxsinj,sfmin)
                              ENDIF
!           END IF ROTOK THEN ... ELSE
!
!           In the case of cancellation in updating SVA(q), SVA(p)
!           .. recompute SVA(q), SVA(p)
                              IF ( (Sva(q)/aaqq)**2<=rooteps ) THEN
                                 IF ( (aaqq<rootbig) .AND.              &
     &                                (aaqq>rootsfmin) ) THEN
                                    Sva(q) = DZNRM2(M,A(1,q),1)
                                 ELSE
                                    t = ZERO
                                    aaqq = ONE
                                    CALL ZLASSQ(M,A(1,q),1,t,aaqq)
                                    Sva(q) = t*SQRT(aaqq)
                                 ENDIF
                              ENDIF
                              IF ( (aapp/aapp0)**2<=rooteps ) THEN
                                 IF ( (aapp<rootbig) .AND.              &
     &                                (aapp>rootsfmin) ) THEN
                                    aapp = DZNRM2(M,A(1,p),1)
                                 ELSE
                                    t = ZERO
                                    aapp = ONE
                                    CALL ZLASSQ(M,A(1,p),1,t,aapp)
                                    aapp = t*SQRT(aapp)
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
            Sva(N) = DZNRM2(M,A(1,N),1)
         ELSE
            t = ZERO
            aapp = ONE
            CALL ZLASSQ(M,A(1,N),1,t,aapp)
            Sva(N) = t*SQRT(aapp)
         ENDIF
!
!     Additional steering devices
!
         IF ( (i<swband) .AND. ((mxaapq<=roottol) .OR. (iswrot<=N)) )   &
     &        swband = i
!
         IF ( (i>swband+1) .AND. (mxaapq<SQRT(DBLE(N))*tol) .AND.       &
     &        (DBLE(N)*mxaapq*mxsinj<tol) ) GOTO 100
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
         q = IDAMAX(N-p+1,Sva(p),1) + p - 1
         IF ( p/=q ) THEN
            temp1 = Sva(p)
            Sva(p) = Sva(q)
            Sva(q) = temp1
            CALL ZSWAP(M,A(1,p),1,A(1,q),1)
            IF ( rsvec ) CALL ZSWAP(mvl,V(1,p),1,V(1,q),1)
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
         DO p = 1 , n4
!            CALL ZDSCAL( M, ONE / SVA( p ), A( 1, p ), 1 )
            CALL ZLASCL('G',0,0,Sva(p),ONE,M,1,A(1,p),M,ierr)
         ENDDO
      ENDIF
!
!     Scale the product of Jacobi rotations.
!
      IF ( rsvec ) THEN
         DO p = 1 , N
            temp1 = ONE/DZNRM2(mvl,V(1,p),1)
            CALL ZDSCAL(mvl,temp1,V(1,p),1)
         ENDDO
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
      Rwork(1) = skl
!     The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE
!     then some of the singular values may overflow or underflow and
!     the spectrum is given in this factored representation.
!
      Rwork(2) = DBLE(n4)
!     N4 is the number of computed nonzero singular values of A.
!
      Rwork(3) = DBLE(n2)
!     N2 is the number of singular values of A greater than SFMIN.
!     If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers
!     that may carry some information.
!
      Rwork(4) = DBLE(i)
!     i is the index of the last sweep before declaring convergence.
!
      Rwork(5) = mxaapq
!     MXAAPQ is the largest absolute value of scaled pivots in the
!     last sweep
!
      Rwork(6) = mxsinj
!     MXSINJ is the largest absolute value of the sines of Jacobi angles
!     in the last sweep
!
!     ..
!     .. END OF ZGESVJ
!     ..
      END SUBROUTINE ZGESVJ
