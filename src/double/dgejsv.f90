!*==dgejsv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DGEJSV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGEJSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgejsv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgejsv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgejsv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP,
!                          M, N, A, LDA, SVA, U, LDU, V, LDV,
!                          WORK, LWORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       IMPLICIT    NONE
!       INTEGER     INFO, LDA, LDU, LDV, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A( LDA, * ), SVA( N ), U( LDU, * ), V( LDV, * ),
!      $            WORK( LWORK )
!       INTEGER     IWORK( * )
!       CHARACTER*1 JOBA, JOBP, JOBR, JOBT, JOBU, JOBV
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEJSV computes the singular value decomposition (SVD) of a real M-by-N
!> matrix [A], where M >= N. The SVD of [A] is written as
!>
!>              [A] = [U] * [SIGMA] * [V]^t,
!>
!> where [SIGMA] is an N-by-N (M-by-N) matrix which is zero except for its N
!> diagonal elements, [U] is an M-by-N (or M-by-M) orthonormal matrix, and
!> [V] is an N-by-N orthogonal matrix. The diagonal elements of [SIGMA] are
!> the singular values of [A]. The columns of [U] and [V] are the left and
!> the right singular vectors of [A], respectively. The matrices [U] and [V]
!> are computed and stored in the arrays U and V, respectively. The diagonal
!> of [SIGMA] is computed and stored in the array SVA.
!> DGEJSV can sometimes compute tiny singular values and their singular vectors much
!> more accurately than other SVD routines, see below under Further Details.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBA
!> \verbatim
!>          JOBA is CHARACTER*1
!>        Specifies the level of accuracy:
!>       = 'C': This option works well (high relative accuracy) if A = B * D,
!>             with well-conditioned B and arbitrary diagonal matrix D.
!>             The accuracy cannot be spoiled by COLUMN scaling. The
!>             accuracy of the computed output depends on the condition of
!>             B, and the procedure aims at the best theoretical accuracy.
!>             The relative error max_{i=1:N}|d sigma_i| / sigma_i is
!>             bounded by f(M,N)*epsilon* cond(B), independent of D.
!>             The input matrix is preprocessed with the QRF with column
!>             pivoting. This initial preprocessing and preconditioning by
!>             a rank revealing QR factorization is common for all values of
!>             JOBA. Additional actions are specified as follows:
!>       = 'E': Computation as with 'C' with an additional estimate of the
!>             condition number of B. It provides a realistic error bound.
!>       = 'F': If A = D1 * C * D2 with ill-conditioned diagonal scalings
!>             D1, D2, and well-conditioned matrix C, this option gives
!>             higher accuracy than the 'C' option. If the structure of the
!>             input matrix is not known, and relative accuracy is
!>             desirable, then this option is advisable. The input matrix A
!>             is preprocessed with QR factorization with FULL (row and
!>             column) pivoting.
!>       = 'G': Computation as with 'F' with an additional estimate of the
!>             condition number of B, where A=D*B. If A has heavily weighted
!>             rows, then using this condition number gives too pessimistic
!>             error bound.
!>       = 'A': Small singular values are the noise and the matrix is treated
!>             as numerically rank deficient. The error in the computed
!>             singular values is bounded by f(m,n)*epsilon*||A||.
!>             The computed SVD A = U * S * V^t restores A up to
!>             f(m,n)*epsilon*||A||.
!>             This gives the procedure the licence to discard (set to zero)
!>             all singular values below N*epsilon*||A||.
!>       = 'R': Similar as in 'A'. Rank revealing property of the initial
!>             QR factorization is used do reveal (using triangular factor)
!>             a gap sigma_{r+1} < epsilon * sigma_r in which case the
!>             numerical RANK is declared to be r. The SVD is computed with
!>             absolute error bounds, but more accurately than with 'A'.
!> \endverbatim
!>
!> \param[in] JOBU
!> \verbatim
!>          JOBU is CHARACTER*1
!>        Specifies whether to compute the columns of U:
!>       = 'U': N columns of U are returned in the array U.
!>       = 'F': full set of M left sing. vectors is returned in the array U.
!>       = 'W': U may be used as workspace of length M*N. See the description
!>             of U.
!>       = 'N': U is not computed.
!> \endverbatim
!>
!> \param[in] JOBV
!> \verbatim
!>          JOBV is CHARACTER*1
!>        Specifies whether to compute the matrix V:
!>       = 'V': N columns of V are returned in the array V; Jacobi rotations
!>             are not explicitly accumulated.
!>       = 'J': N columns of V are returned in the array V, but they are
!>             computed as the product of Jacobi rotations. This option is
!>             allowed only if JOBU .NE. 'N', i.e. in computing the full SVD.
!>       = 'W': V may be used as workspace of length N*N. See the description
!>             of V.
!>       = 'N': V is not computed.
!> \endverbatim
!>
!> \param[in] JOBR
!> \verbatim
!>          JOBR is CHARACTER*1
!>        Specifies the RANGE for the singular values. Issues the licence to
!>        set to zero small positive singular values if they are outside
!>        specified range. If A .NE. 0 is scaled so that the largest singular
!>        value of c*A is around DSQRT(BIG), BIG=SLAMCH('O'), then JOBR issues
!>        the licence to kill columns of A whose norm in c*A is less than
!>        DSQRT(SFMIN) (for JOBR = 'R'), or less than SMALL=SFMIN/EPSLN,
!>        where SFMIN=SLAMCH('S'), EPSLN=SLAMCH('E').
!>       = 'N': Do not kill small columns of c*A. This option assumes that
!>             BLAS and QR factorizations and triangular solvers are
!>             implemented to work in that range. If the condition of A
!>             is greater than BIG, use DGESVJ.
!>       = 'R': RESTRICTED range for sigma(c*A) is [DSQRT(SFMIN), DSQRT(BIG)]
!>             (roughly, as described above). This option is recommended.
!>                                            ~~~~~~~~~~~~~~~~~~~~~~~~~~~
!>        For computing the singular values in the FULL range [SFMIN,BIG]
!>        use DGESVJ.
!> \endverbatim
!>
!> \param[in] JOBT
!> \verbatim
!>          JOBT is CHARACTER*1
!>        If the matrix is square then the procedure may determine to use
!>        transposed A if A^t seems to be better with respect to convergence.
!>        If the matrix is not square, JOBT is ignored. This is subject to
!>        changes in the future.
!>        The decision is based on two values of entropy over the adjoint
!>        orbit of A^t * A. See the descriptions of WORK(6) and WORK(7).
!>       = 'T': transpose if entropy test indicates possibly faster
!>        convergence of Jacobi process if A^t is taken as input. If A is
!>        replaced with A^t, then the row pivoting is included automatically.
!>       = 'N': do not speculate.
!>        This option can be used to compute only the singular values, or the
!>        full SVD (U, SIGMA and V). For only one set of singular vectors
!>        (U or V), the caller should provide both U and V, as one of the
!>        matrices is used as workspace if the matrix A is transposed.
!>        The implementer can easily remove this constraint and make the
!>        code more complicated. See the descriptions of U and V.
!> \endverbatim
!>
!> \param[in] JOBP
!> \verbatim
!>          JOBP is CHARACTER*1
!>        Issues the licence to introduce structured perturbations to drown
!>        denormalized numbers. This licence should be active if the
!>        denormals are poorly implemented, causing slow computation,
!>        especially in cases of fast convergence (!). For details see [1,2].
!>        For the sake of simplicity, this perturbations are included only
!>        when the full SVD or only the singular values are requested. The
!>        implementer/user can easily add the perturbation for the cases of
!>        computing one set of singular vectors.
!>       = 'P': introduce perturbation
!>       = 'N': do not perturb
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>         The number of rows of the input matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The number of columns of the input matrix A. M >= N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
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
!>          - For WORK(1)/WORK(2) = ONE: The singular values of A. During the
!>            computation SVA contains Euclidean column norms of the
!>            iterated matrices in the array A.
!>          - For WORK(1) .NE. WORK(2): The singular values of A are
!>            (WORK(1)/WORK(2)) * SVA(1:N). This factored form is used if
!>            sigma_max(A) overflows or if small singular values have been
!>            saved from underflow by scaling the input matrix A.
!>          - If JOBR='R' then some of the singular values may be returned
!>            as exact zeros obtained by "set to zero" because they are
!>            below the numerical rank threshold or are denormalized numbers.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension ( LDU, N )
!>          If JOBU = 'U', then U contains on exit the M-by-N matrix of
!>                         the left singular vectors.
!>          If JOBU = 'F', then U contains on exit the M-by-M matrix of
!>                         the left singular vectors, including an ONB
!>                         of the orthogonal complement of the Range(A).
!>          If JOBU = 'W'  .AND. (JOBV = 'V' .AND. JOBT = 'T' .AND. M = N),
!>                         then U is used as workspace if the procedure
!>                         replaces A with A^t. In that case, [V] is computed
!>                         in U as left singular vectors of A^t and then
!>                         copied back to the V array. This 'W' option is just
!>                         a reminder to the caller that in this case U is
!>                         reserved as workspace of length N*N.
!>          If JOBU = 'N'  U is not referenced, unless JOBT='T'.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U,  LDU >= 1.
!>          IF  JOBU = 'U' or 'F' or 'W',  then LDU >= M.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension ( LDV, N )
!>          If JOBV = 'V', 'J' then V contains on exit the N-by-N matrix of
!>                         the right singular vectors;
!>          If JOBV = 'W', AND (JOBU = 'U' AND JOBT = 'T' AND M = N),
!>                         then V is used as workspace if the pprocedure
!>                         replaces A with A^t. In that case, [U] is computed
!>                         in V as right singular vectors of A^t and then
!>                         copied back to the U array. This 'W' option is just
!>                         a reminder to the caller that in this case V is
!>                         reserved as workspace of length N*N.
!>          If JOBV = 'N'  V is not referenced, unless JOBT='T'.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V,  LDV >= 1.
!>          If JOBV = 'V' or 'J' or 'W', then LDV >= N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!>          On exit, if N > 0 .AND. M > 0 (else not referenced),
!>          WORK(1) = SCALE = WORK(2) / WORK(1) is the scaling factor such
!>                    that SCALE*SVA(1:N) are the computed singular values
!>                    of A. (See the description of SVA().)
!>          WORK(2) = See the description of WORK(1).
!>          WORK(3) = SCONDA is an estimate for the condition number of
!>                    column equilibrated A. (If JOBA = 'E' or 'G')
!>                    SCONDA is an estimate of DSQRT(||(R^t * R)^(-1)||_1).
!>                    It is computed using DPOCON. It holds
!>                    N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
!>                    where R is the triangular factor from the QRF of A.
!>                    However, if R is truncated and the numerical rank is
!>                    determined to be strictly smaller than N, SCONDA is
!>                    returned as -1, thus indicating that the smallest
!>                    singular values might be lost.
!>
!>          If full SVD is needed, the following two condition numbers are
!>          useful for the analysis of the algorithm. They are provided for
!>          a developer/implementer who is familiar with the details of
!>          the method.
!>
!>          WORK(4) = an estimate of the scaled condition number of the
!>                    triangular factor in the first QR factorization.
!>          WORK(5) = an estimate of the scaled condition number of the
!>                    triangular factor in the second QR factorization.
!>          The following two parameters are computed if JOBT = 'T'.
!>          They are provided for a developer/implementer who is familiar
!>          with the details of the method.
!>
!>          WORK(6) = the entropy of A^t*A :: this is the Shannon entropy
!>                    of diag(A^t*A) / Trace(A^t*A) taken as point in the
!>                    probability simplex.
!>          WORK(7) = the entropy of A*A^t.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          Length of WORK to confirm proper allocation of work space.
!>          LWORK depends on the job:
!>
!>          If only SIGMA is needed (JOBU = 'N', JOBV = 'N') and
!>            -> .. no scaled condition estimate required (JOBE = 'N'):
!>               LWORK >= max(2*M+N,4*N+1,7). This is the minimal requirement.
!>               ->> For optimal performance (blocked code) the optimal value
!>               is LWORK >= max(2*M+N,3*N+(N+1)*NB,7). Here NB is the optimal
!>               block size for DGEQP3 and DGEQRF.
!>               In general, optimal LWORK is computed as
!>               LWORK >= max(2*M+N,N+LWORK(DGEQP3),N+LWORK(DGEQRF), 7).
!>            -> .. an estimate of the scaled condition number of A is
!>               required (JOBA='E', 'G'). In this case, LWORK is the maximum
!>               of the above and N*N+4*N, i.e. LWORK >= max(2*M+N,N*N+4*N,7).
!>               ->> For optimal performance (blocked code) the optimal value
!>               is LWORK >= max(2*M+N,3*N+(N+1)*NB, N*N+4*N, 7).
!>               In general, the optimal length LWORK is computed as
!>               LWORK >= max(2*M+N,N+LWORK(DGEQP3),N+LWORK(DGEQRF),
!>                                                     N+N*N+LWORK(DPOCON),7).
!>
!>          If SIGMA and the right singular vectors are needed (JOBV = 'V'),
!>            -> the minimal requirement is LWORK >= max(2*M+N,4*N+1,7).
!>            -> For optimal performance, LWORK >= max(2*M+N,3*N+(N+1)*NB,7),
!>               where NB is the optimal block size for DGEQP3, DGEQRF, DGELQF,
!>               DORMLQ. In general, the optimal length LWORK is computed as
!>               LWORK >= max(2*M+N,N+LWORK(DGEQP3), N+LWORK(DPOCON),
!>                       N+LWORK(DGELQF), 2*N+LWORK(DGEQRF), N+LWORK(DORMLQ)).
!>
!>          If SIGMA and the left singular vectors are needed
!>            -> the minimal requirement is LWORK >= max(2*M+N,4*N+1,7).
!>            -> For optimal performance:
!>               if JOBU = 'U' :: LWORK >= max(2*M+N,3*N+(N+1)*NB,7),
!>               if JOBU = 'F' :: LWORK >= max(2*M+N,3*N+(N+1)*NB,N+M*NB,7),
!>               where NB is the optimal block size for DGEQP3, DGEQRF, DORMQR.
!>               In general, the optimal length LWORK is computed as
!>               LWORK >= max(2*M+N,N+LWORK(DGEQP3),N+LWORK(DPOCON),
!>                        2*N+LWORK(DGEQRF), N+LWORK(DORMQR)).
!>               Here LWORK(DORMQR) equals N*NB (for JOBU = 'U') or
!>               M*NB (for JOBU = 'F').
!>
!>          If the full SVD is needed: (JOBU = 'U' or JOBU = 'F') and
!>            -> if JOBV = 'V'
!>               the minimal requirement is LWORK >= max(2*M+N,6*N+2*N*N).
!>            -> if JOBV = 'J' the minimal requirement is
!>               LWORK >= max(2*M+N, 4*N+N*N,2*N+N*N+6).
!>            -> For optimal performance, LWORK should be additionally
!>               larger than N+M*NB, where NB is the optimal block size
!>               for DORMQR.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (M+3*N).
!>          On exit,
!>          IWORK(1) = the numerical rank determined after the initial
!>                     QR factorization with pivoting. See the descriptions
!>                     of JOBA and JOBR.
!>          IWORK(2) = the number of the computed nonzero singular values
!>          IWORK(3) = if nonzero, a warning message:
!>                     If IWORK(3) = 1 then some of the column norms of A
!>                     were denormalized floats. The requested high accuracy
!>                     is not warranted by the data.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           < 0:  if INFO = -i, then the i-th argument had an illegal value.
!>           = 0:  successful exit;
!>           > 0:  DGEJSV  did not converge in the maximal allowed number
!>                 of sweeps. The computed values may be inaccurate.
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
!> \date June 2016
!
!> \ingroup doubleGEsing
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  DGEJSV implements a preconditioned Jacobi SVD algorithm. It uses DGEQP3,
!>  DGEQRF, and DGELQF as preprocessors and preconditioners. Optionally, an
!>  additional row pivoting can be used as a preprocessor, which in some
!>  cases results in much higher accuracy. An example is matrix A with the
!>  structure A = D1 * C * D2, where D1, D2 are arbitrarily ill-conditioned
!>  diagonal matrices and C is well-conditioned matrix. In that case, complete
!>  pivoting in the first QR factorizations provides accuracy dependent on the
!>  condition number of C, and independent of D1, D2. Such higher accuracy is
!>  not completely understood theoretically, but it works well in practice.
!>  Further, if A can be written as A = B*D, with well-conditioned B and some
!>  diagonal D, then the high accuracy is guaranteed, both theoretically and
!>  in software, independent of D. For more details see [1], [2].
!>     The computational range for the singular values can be the full range
!>  ( UNDERFLOW,OVERFLOW ), provided that the machine arithmetic and the BLAS
!>  & LAPACK routines called by DGEJSV are implemented to work in that range.
!>  If that is not the case, then the restriction for safe computation with
!>  the singular values in the range of normalized IEEE numbers is that the
!>  spectral condition number kappa(A)=sigma_max(A)/sigma_min(A) does not
!>  overflow. This code (DGEJSV) is best used in this restricted range,
!>  meaning that singular values of magnitude below ||A||_2 / DLAMCH('O') are
!>  returned as zeros. See JOBR for details on this.
!>     Further, this implementation is somewhat slower than the one described
!>  in [1,2] due to replacement of some non-LAPACK components, and because
!>  the choice of some tuning parameters in the iterative part (DGESVJ) is
!>  left to the implementer on a particular machine.
!>     The rank revealing QR factorization (in this code: DGEQP3) should be
!>  implemented as in [3]. We have a new version of DGEQP3 under development
!>  that is more robust than the current one in LAPACK, with a cleaner cut in
!>  rank deficient cases. It will be available in the SIGMA library [4].
!>  If M is much larger than N, it is obvious that the initial QRF with
!>  column pivoting can be preprocessed by the QRF without pivoting. That
!>  well known trick is not used in DGEJSV because in some cases heavy row
!>  weighting can be treated with complete pivoting. The overhead in cases
!>  M much larger than N is then only due to pivoting, but the benefits in
!>  terms of accuracy have prevailed. The implementer/user can incorporate
!>  this extra QRF step easily. The implementer can also improve data movement
!>  (matrix transpose, matrix copy, matrix transposed copy) - this
!>  implementation of DGEJSV uses only the simplest, naive data movement.
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>  Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)
!
!> \par References:
!  ================
!>
!> \verbatim
!>
!> [1] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I.
!>     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342.
!>     LAPACK Working note 169.
!> [2] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II.
!>     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362.
!>     LAPACK Working note 170.
!> [3] Z. Drmac and Z. Bujanovic: On the failure of rank-revealing QR
!>     factorization software - a case study.
!>     ACM Trans. Math. Softw. Vol. 35, No 2 (2008), pp. 1-28.
!>     LAPACK Working note 176.
!> [4] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV,
!>     QSVD, (H,K)-SVD computations.
!>     Department of Mathematics, University of Zagreb, 2008.
!> \endverbatim
!
!>  \par Bugs, examples and comments:
!   =================================
!>
!>  Please report all bugs and send interesting examples and/or comments to
!>  drmac@math.hr. Thank you.
!>
!  =====================================================================
      SUBROUTINE DGEJSV(Joba,Jobu,Jobv,Jobr,Jobt,Jobp,M,N,A,Lda,Sva,U,  &
     &                  Ldu,V,Ldv,Work,Lwork,Iwork,Info)
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      IMPLICIT NONE
!*--DGEJSV486
      INTEGER Info , Lda , Ldu , Ldv , Lwork , M , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , Sva(N) , U(Ldu,*) , V(Ldv,*) ,        &
     &                 Work(Lwork)
      INTEGER Iwork(*)
      CHARACTER*1 Joba , Jobp , Jobr , Jobt , Jobu , Jobv
!     ..
!
!  ===========================================================================
!
!     .. Local Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION aapp , aaqq , aatmax , aatmin , big , big1 ,     &
     &                 cond_ok , condr1 , condr2 , entra , entrat ,     &
     &                 epsln , maxprj , scalem , sconda , sfmin ,       &
     &                 small , temp1 , uscal1 , uscal2 , xsc
      INTEGER ierr , n1 , nr , numrank , p , q , warning
      LOGICAL almort , defr , errest , goscal , jracc , kill , lsvec ,  &
     &        l2aber , l2kill , l2pert , l2rank , l2tran , noscal ,     &
     &        rowpiv , rsvec , transp
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DABS , DLOG , MAX , MIN , DBLE , IDNINT , DSIGN , DSQRT
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DNRM2
      INTEGER IDAMAX
      LOGICAL LSAME
      EXTERNAL IDAMAX , LSAME , DLAMCH , DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DGELQF , DGEQP3 , DGEQRF , DLACPY , DLASCL ,     &
     &         DLASET , DLASSQ , DLASWP , DORGQR , DORMLQ , DORMQR ,    &
     &         DPOCON , DSCAL , DSWAP , DTRSM , XERBLA
!
      EXTERNAL DGESVJ
!     ..
!
!     Test the input arguments
!
      lsvec = LSAME(Jobu,'U') .OR. LSAME(Jobu,'F')
      jracc = LSAME(Jobv,'J')
      rsvec = LSAME(Jobv,'V') .OR. jracc
      rowpiv = LSAME(Joba,'F') .OR. LSAME(Joba,'G')
      l2rank = LSAME(Joba,'R')
      l2aber = LSAME(Joba,'A')
      errest = LSAME(Joba,'E') .OR. LSAME(Joba,'G')
      l2tran = LSAME(Jobt,'T')
      l2kill = LSAME(Jobr,'R')
      defr = LSAME(Jobr,'N')
      l2pert = LSAME(Jobp,'P')
!
      IF ( .NOT.(rowpiv .OR. l2rank .OR. l2aber .OR. errest .OR.        &
     &     LSAME(Joba,'C')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(lsvec .OR. LSAME(Jobu,'N') .OR. LSAME(Jobu,'W')) ) &
     &         THEN
         Info = -2
      ELSEIF ( .NOT.(rsvec .OR. LSAME(Jobv,'N') .OR. LSAME(Jobv,'W'))   &
     &         .OR. (jracc .AND. (.NOT.lsvec)) ) THEN
         Info = -3
      ELSEIF ( .NOT.(l2kill .OR. defr) ) THEN
         Info = -4
      ELSEIF ( .NOT.(l2tran .OR. LSAME(Jobt,'N')) ) THEN
         Info = -5
      ELSEIF ( .NOT.(l2pert .OR. LSAME(Jobp,'N')) ) THEN
         Info = -6
      ELSEIF ( M<0 ) THEN
         Info = -7
      ELSEIF ( (N<0) .OR. (N>M) ) THEN
         Info = -8
      ELSEIF ( Lda<M ) THEN
         Info = -10
      ELSEIF ( lsvec .AND. (Ldu<M) ) THEN
         Info = -13
      ELSEIF ( rsvec .AND. (Ldv<N) ) THEN
         Info = -15
      ELSEIF ( (.NOT.(lsvec .OR. rsvec .OR. errest) .AND.               &
     &         (Lwork<MAX(7,4*N+1,2*M+N))) .OR.                         &
     &         (.NOT.(lsvec .OR. rsvec) .AND. errest .AND.              &
     &         (Lwork<MAX(7,4*N+N*N,2*M+N))) .OR.                       &
     &         (lsvec .AND. (.NOT.rsvec) .AND.                          &
     &         (Lwork<MAX(7,2*M+N,4*N+1))) .OR.                         &
     &         (rsvec .AND. (.NOT.lsvec) .AND.                          &
     &         (Lwork<MAX(7,2*M+N,4*N+1))) .OR.                         &
     &         (lsvec .AND. rsvec .AND. (.NOT.jracc) .AND.              &
     &         (Lwork<MAX(2*M+N,6*N+2*N*N))) .OR.                       &
     &         (lsvec .AND. rsvec .AND. jracc .AND.                     &
     &         Lwork<MAX(2*M+N,4*N+N*N,2*N+N*N+6)) ) THEN
         Info = -17
      ELSE
!        #:)
         Info = 0
      ENDIF
!
      IF ( Info/=0 ) THEN
!       #:(
         CALL XERBLA('DGEJSV',-Info)
         RETURN
      ENDIF
!
!     Quick return for void matrix (Y3K safe)
! #:)
      IF ( (M==0) .OR. (N==0) ) THEN
         Iwork(1:3) = 0
         Work(1:7) = 0
         RETURN
      ENDIF
!
!     Determine whether the matrix U should be M x N or M x M
!
      IF ( lsvec ) THEN
         n1 = N
         IF ( LSAME(Jobu,'F') ) n1 = M
      ENDIF
!
!     Set numerical parameters
!
!!    NOTE: Make sure DLAMCH() does not fail on the target architecture.
!
      epsln = DLAMCH('Epsilon')
      sfmin = DLAMCH('SafeMinimum')
      small = sfmin/epsln
      big = DLAMCH('O')
!     BIG   = ONE / SFMIN
!
!     Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N
!
!(!)  If necessary, scale SVA() to protect the largest norm from
!     overflow. It is possible that this scaling pushes the smallest
!     column norm left from the underflow threshold (extreme case).
!
      scalem = ONE/DSQRT(DBLE(M)*DBLE(N))
      noscal = .TRUE.
      goscal = .TRUE.
      DO p = 1 , N
         aapp = ZERO
         aaqq = ONE
         CALL DLASSQ(M,A(1,p),1,aapp,aaqq)
         IF ( aapp>big ) THEN
            Info = -9
            CALL XERBLA('DGEJSV',-Info)
            RETURN
         ENDIF
         aaqq = DSQRT(aaqq)
         IF ( (aapp<(big/aaqq)) .AND. noscal ) THEN
            Sva(p) = aapp*aaqq
         ELSE
            noscal = .FALSE.
            Sva(p) = aapp*(aaqq*scalem)
            IF ( goscal ) THEN
               goscal = .FALSE.
               CALL DSCAL(p-1,scalem,Sva,1)
            ENDIF
         ENDIF
      ENDDO
!
      IF ( noscal ) scalem = ONE
!
      aapp = ZERO
      aaqq = big
      DO p = 1 , N
         aapp = MAX(aapp,Sva(p))
         IF ( Sva(p)/=ZERO ) aaqq = MIN(aaqq,Sva(p))
      ENDDO
!
!     Quick return for zero M x N matrix
! #:)
      IF ( aapp==ZERO ) THEN
         IF ( lsvec ) CALL DLASET('G',M,n1,ZERO,ONE,U,Ldu)
         IF ( rsvec ) CALL DLASET('G',N,N,ZERO,ONE,V,Ldv)
         Work(1) = ONE
         Work(2) = ONE
         IF ( errest ) Work(3) = ONE
         IF ( lsvec .AND. rsvec ) THEN
            Work(4) = ONE
            Work(5) = ONE
         ENDIF
         IF ( l2tran ) THEN
            Work(6) = ZERO
            Work(7) = ZERO
         ENDIF
         Iwork(1) = 0
         Iwork(2) = 0
         Iwork(3) = 0
         RETURN
      ENDIF
!
!     Issue warning if denormalized column norms detected. Override the
!     high relative accuracy request. Issue licence to kill columns
!     (set them to zero) whose norm is less than sigma_max / BIG (roughly).
! #:(
      warning = 0
      IF ( aaqq<=sfmin ) THEN
         l2rank = .TRUE.
         l2kill = .TRUE.
         warning = 1
      ENDIF
!
!     Quick return for one-column matrix
! #:)
      IF ( N==1 ) THEN
!
         IF ( lsvec ) THEN
            CALL DLASCL('G',0,0,Sva(1),scalem,M,1,A(1,1),Lda,ierr)
            CALL DLACPY('A',M,1,A,Lda,U,Ldu)
!           computing all M left singular vectors of the M x 1 matrix
            IF ( n1/=N ) THEN
               CALL DGEQRF(M,N,U,Ldu,Work,Work(N+1),Lwork-N,ierr)
               CALL DORGQR(M,n1,1,U,Ldu,Work,Work(N+1),Lwork-N,ierr)
               CALL DCOPY(M,A(1,1),1,U(1,1),1)
            ENDIF
         ENDIF
         IF ( rsvec ) V(1,1) = ONE
         IF ( Sva(1)<(big*scalem) ) THEN
            Sva(1) = Sva(1)/scalem
            scalem = ONE
         ENDIF
         Work(1) = ONE/scalem
         Work(2) = ONE
         IF ( Sva(1)/=ZERO ) THEN
            Iwork(1) = 1
            IF ( (Sva(1)/scalem)>=sfmin ) THEN
               Iwork(2) = 1
            ELSE
               Iwork(2) = 0
            ENDIF
         ELSE
            Iwork(1) = 0
            Iwork(2) = 0
         ENDIF
         Iwork(3) = 0
         IF ( errest ) Work(3) = ONE
         IF ( lsvec .AND. rsvec ) THEN
            Work(4) = ONE
            Work(5) = ONE
         ENDIF
         IF ( l2tran ) THEN
            Work(6) = ZERO
            Work(7) = ZERO
         ENDIF
         RETURN
!
      ENDIF
!
      transp = .FALSE.
      l2tran = l2tran .AND. (M==N)
!
      aatmax = -ONE
      aatmin = big
      IF ( rowpiv .OR. l2tran ) THEN
!
!     Compute the row norms, needed to determine row pivoting sequence
!     (in the case of heavily row weighted A, row pivoting is strongly
!     advised) and to collect information needed to compare the
!     structures of A * A^t and A^t * A (in the case L2TRAN.EQ..TRUE.).
!
         IF ( l2tran ) THEN
            DO p = 1 , M
               xsc = ZERO
               temp1 = ONE
               CALL DLASSQ(N,A(p,1),Lda,xsc,temp1)
!              DLASSQ gets both the ell_2 and the ell_infinity norm
!              in one pass through the vector
               Work(M+N+p) = xsc*scalem
               Work(N+p) = xsc*(scalem*DSQRT(temp1))
               aatmax = MAX(aatmax,Work(N+p))
               IF ( Work(N+p)/=ZERO ) aatmin = MIN(aatmin,Work(N+p))
            ENDDO
         ELSE
            DO p = 1 , M
               Work(M+N+p) = scalem*DABS(A(p,IDAMAX(N,A(p,1),Lda)))
               aatmax = MAX(aatmax,Work(M+N+p))
               aatmin = MIN(aatmin,Work(M+N+p))
            ENDDO
         ENDIF
!
      ENDIF
!
!     For square matrix A try to determine whether A^t  would be  better
!     input for the preconditioned Jacobi SVD, with faster convergence.
!     The decision is based on an O(N) function of the vector of column
!     and row norms of A, based on the Shannon entropy. This should give
!     the right choice in most cases when the difference actually matters.
!     It may fail and pick the slower converging side.
!
      entra = ZERO
      entrat = ZERO
      IF ( l2tran ) THEN
!
         xsc = ZERO
         temp1 = ONE
         CALL DLASSQ(N,Sva,1,xsc,temp1)
         temp1 = ONE/temp1
!
         entra = ZERO
         DO p = 1 , N
            big1 = ((Sva(p)/xsc)**2)*temp1
            IF ( big1/=ZERO ) entra = entra + big1*DLOG(big1)
         ENDDO
         entra = -entra/DLOG(DBLE(N))
!
!        Now, SVA().^2/Trace(A^t * A) is a point in the probability simplex.
!        It is derived from the diagonal of  A^t * A.  Do the same with the
!        diagonal of A * A^t, compute the entropy of the corresponding
!        probability distribution. Note that A * A^t and A^t * A have the
!        same trace.
!
         entrat = ZERO
         DO p = N + 1 , N + M
            big1 = ((Work(p)/xsc)**2)*temp1
            IF ( big1/=ZERO ) entrat = entrat + big1*DLOG(big1)
         ENDDO
         entrat = -entrat/DLOG(DBLE(M))
!
!        Analyze the entropies and decide A or A^t. Smaller entropy
!        usually means better input for the algorithm.
!
         transp = (entrat<entra)
!
!        If A^t is better than A, transpose A.
!
         IF ( transp ) THEN
!           In an optimal implementation, this trivial transpose
!           should be replaced with faster transpose.
            DO p = 1 , N - 1
               DO q = p + 1 , N
                  temp1 = A(q,p)
                  A(q,p) = A(p,q)
                  A(p,q) = temp1
               ENDDO
            ENDDO
            DO p = 1 , N
               Work(M+N+p) = Sva(p)
               Sva(p) = Work(N+p)
            ENDDO
            temp1 = aapp
            aapp = aatmax
            aatmax = temp1
            temp1 = aaqq
            aaqq = aatmin
            aatmin = temp1
            kill = lsvec
            lsvec = rsvec
            rsvec = kill
            IF ( lsvec ) n1 = N
!
            rowpiv = .TRUE.
         ENDIF
!
      ENDIF
!     END IF L2TRAN
!
!     Scale the matrix so that its maximal singular value remains less
!     than DSQRT(BIG) -- the matrix is scaled so that its maximal column
!     has Euclidean norm equal to DSQRT(BIG/N). The only reason to keep
!     DSQRT(BIG) instead of BIG is the fact that DGEJSV uses LAPACK and
!     BLAS routines that, in some implementations, are not capable of
!     working in the full interval [SFMIN,BIG] and that they may provoke
!     overflows in the intermediate results. If the singular values spread
!     from SFMIN to BIG, then DGESVJ will compute them. So, in that case,
!     one should use DGESVJ instead of DGEJSV.
!
      big1 = DSQRT(big)
      temp1 = DSQRT(big/DBLE(N))
!
      CALL DLASCL('G',0,0,aapp,temp1,N,1,Sva,N,ierr)
      IF ( aaqq>(aapp*sfmin) ) THEN
         aaqq = (aaqq/aapp)*temp1
      ELSE
         aaqq = (aaqq*temp1)/aapp
      ENDIF
      temp1 = temp1*scalem
      CALL DLASCL('G',0,0,aapp,temp1,M,N,A,Lda,ierr)
!
!     To undo scaling at the end of this procedure, multiply the
!     computed singular values with USCAL2 / USCAL1.
!
      uscal1 = temp1
      uscal2 = aapp
!
      IF ( l2kill ) THEN
!        L2KILL enforces computation of nonzero singular values in
!        the restricted range of condition number of the initial A,
!        sigma_max(A) / sigma_min(A) approx. DSQRT(BIG)/DSQRT(SFMIN).
         xsc = DSQRT(sfmin)
      ELSE
         xsc = small
!
!        Now, if the condition number of A is too big,
!        sigma_max(A) / sigma_min(A) .GT. DSQRT(BIG/N) * EPSLN / SFMIN,
!        as a precaution measure, the full SVD is computed using DGESVJ
!        with accumulated Jacobi rotations. This provides numerically
!        more robust computation, at the cost of slightly increased run
!        time. Depending on the concrete implementation of BLAS and LAPACK
!        (i.e. how they behave in presence of extreme ill-conditioning) the
!        implementor may decide to remove this switch.
         IF ( (aaqq<DSQRT(sfmin)) .AND. lsvec .AND. rsvec )             &
     &        jracc = .TRUE.
!
      ENDIF
      IF ( aaqq<xsc ) THEN
         DO p = 1 , N
            IF ( Sva(p)<xsc ) THEN
               CALL DLASET('A',M,1,ZERO,ZERO,A(1,p),Lda)
               Sva(p) = ZERO
            ENDIF
         ENDDO
      ENDIF
!
!     Preconditioning using QR factorization with pivoting
!
      IF ( rowpiv ) THEN
!        Optional row permutation (Bjoerck row pivoting):
!        A result by Cox and Higham shows that the Bjoerck's
!        row pivoting combined with standard column pivoting
!        has similar effect as Powell-Reid complete pivoting.
!        The ell-infinity norms of A are made nonincreasing.
         DO p = 1 , M - 1
            q = IDAMAX(M-p+1,Work(M+N+p),1) + p - 1
            Iwork(2*N+p) = q
            IF ( p/=q ) THEN
               temp1 = Work(M+N+p)
               Work(M+N+p) = Work(M+N+q)
               Work(M+N+q) = temp1
            ENDIF
         ENDDO
         CALL DLASWP(N,A,Lda,1,M-1,Iwork(2*N+1),1)
      ENDIF
!
!     End of the preparation phase (scaling, optional sorting and
!     transposing, optional flushing of small columns).
!
!     Preconditioning
!
!     If the full SVD is needed, the right singular vectors are computed
!     from a matrix equation, and for that we need theoretical analysis
!     of the Businger-Golub pivoting. So we use DGEQP3 as the first RR QRF.
!     In all other cases the first RR QRF can be chosen by other criteria
!     (eg speed by replacing global with restricted window pivoting, such
!     as in SGEQPX from TOMS # 782). Good results will be obtained using
!     SGEQPX with properly (!) chosen numerical parameters.
!     Any improvement of DGEQP3 improves overall performance of DGEJSV.
!
!     A * P1 = Q1 * [ R1^t 0]^t:
      DO p = 1 , N
!        .. all columns are free columns
         Iwork(p) = 0
      ENDDO
      CALL DGEQP3(M,N,A,Lda,Iwork,Work,Work(N+1),Lwork-N,ierr)
!
!     The upper triangular matrix R1 from the first QRF is inspected for
!     rank deficiency and possibilities for deflation, or possible
!     ill-conditioning. Depending on the user specified flag L2RANK,
!     the procedure explores possibilities to reduce the numerical
!     rank by inspecting the computed upper triangular factor. If
!     L2RANK or L2ABER are up, then DGEJSV will compute the SVD of
!     A + dA, where ||dA|| <= f(M,N)*EPSLN.
!
      nr = 1
      IF ( l2aber ) THEN
!        Standard absolute error bound suffices. All sigma_i with
!        sigma_i < N*EPSLN*||A|| are flushed to zero. This is an
!        aggressive enforcement of lower numerical rank by introducing a
!        backward error of the order of N*EPSLN*||A||.
         temp1 = DSQRT(DBLE(N))*epsln
         DO p = 2 , N
            IF ( DABS(A(p,p))<(temp1*DABS(A(1,1))) ) EXIT
            nr = nr + 1
         ENDDO
      ELSEIF ( l2rank ) THEN
!        .. similarly as above, only slightly more gentle (less aggressive).
!        Sudden drop on the diagonal of R1 is used as the criterion for
!        close-to-rank-deficient.
         temp1 = DSQRT(sfmin)
         DO p = 2 , N
            IF ( (DABS(A(p,p))<(epsln*DABS(A(p-1,p-1)))) .OR.           &
     &           (DABS(A(p,p))<small) .OR.                              &
     &           (l2kill .AND. (DABS(A(p,p))<temp1)) ) EXIT
            nr = nr + 1
         ENDDO
!
      ELSE
!        The goal is high relative accuracy. However, if the matrix
!        has high scaled condition number the relative accuracy is in
!        general not feasible. Later on, a condition number estimator
!        will be deployed to estimate the scaled condition number.
!        Here we just remove the underflowed part of the triangular
!        factor. This prevents the situation in which the code is
!        working hard to get the accuracy not warranted by the data.
         temp1 = DSQRT(sfmin)
         DO p = 2 , N
            IF ( (DABS(A(p,p))<small) .OR.                              &
     &           (l2kill .AND. (DABS(A(p,p))<temp1)) ) EXIT
            nr = nr + 1
         ENDDO
!
      ENDIF
!
      almort = .FALSE.
      IF ( nr==N ) THEN
         maxprj = ONE
         DO p = 2 , N
            temp1 = DABS(A(p,p))/Sva(Iwork(p))
            maxprj = MIN(maxprj,temp1)
         ENDDO
         IF ( maxprj**2>=ONE-DBLE(N)*epsln ) almort = .TRUE.
      ENDIF
!
!
      sconda = -ONE
      condr1 = -ONE
      condr2 = -ONE
!
      IF ( errest ) THEN
         IF ( N==nr ) THEN
            IF ( rsvec ) THEN
!              .. V is available as workspace
               CALL DLACPY('U',N,N,A,Lda,V,Ldv)
               DO p = 1 , N
                  temp1 = Sva(Iwork(p))
                  CALL DSCAL(p,ONE/temp1,V(1,p),1)
               ENDDO
               CALL DPOCON('U',N,V,Ldv,ONE,temp1,Work(N+1),             &
     &                     Iwork(2*N+M+1),ierr)
            ELSEIF ( lsvec ) THEN
!              .. U is available as workspace
               CALL DLACPY('U',N,N,A,Lda,U,Ldu)
               DO p = 1 , N
                  temp1 = Sva(Iwork(p))
                  CALL DSCAL(p,ONE/temp1,U(1,p),1)
               ENDDO
               CALL DPOCON('U',N,U,Ldu,ONE,temp1,Work(N+1),             &
     &                     Iwork(2*N+M+1),ierr)
            ELSE
               CALL DLACPY('U',N,N,A,Lda,Work(N+1),N)
               DO p = 1 , N
                  temp1 = Sva(Iwork(p))
                  CALL DSCAL(p,ONE/temp1,Work(N+(p-1)*N+1),1)
               ENDDO
!           .. the columns of R are scaled to have unit Euclidean lengths.
               CALL DPOCON('U',N,Work(N+1),N,ONE,temp1,Work(N+N*N+1),   &
     &                     Iwork(2*N+M+1),ierr)
            ENDIF
            sconda = ONE/DSQRT(temp1)
!           SCONDA is an estimate of DSQRT(||(R^t * R)^(-1)||_1).
!           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
         ELSE
            sconda = -ONE
         ENDIF
      ENDIF
!
      l2pert = l2pert .AND. (DABS(A(1,1)/A(nr,nr))>DSQRT(big1))
!     If there is no violent scaling, artificial perturbation is not needed.
!
!     Phase 3:
!
      IF ( .NOT.(rsvec .OR. lsvec) ) THEN
!
!         Singular Values only
!
!         .. transpose A(1:NR,1:N)
         DO p = 1 , MIN(N-1,nr)
            CALL DCOPY(N-p,A(p,p+1),Lda,A(p+1,p),1)
         ENDDO
!
!        The following two DO-loops introduce small relative perturbation
!        into the strict upper triangle of the lower triangular matrix.
!        Small entries below the main diagonal are also changed.
!        This modification is useful if the computing environment does not
!        provide/allow FLUSH TO ZERO underflow, for it prevents many
!        annoying denormalized numbers in case of strongly scaled matrices.
!        The perturbation is structured so that it does not introduce any
!        new perturbation of the singular values, and it does not destroy
!        the job done by the preconditioner.
!        The licence for this perturbation is in the variable L2PERT, which
!        should be .FALSE. if FLUSH TO ZERO underflow is active.
!
         IF ( .NOT.almort ) THEN
!
            IF ( l2pert ) THEN
!              XSC = DSQRT(SMALL)
               xsc = epsln/DBLE(N)
               DO q = 1 , nr
                  temp1 = xsc*DABS(A(q,q))
                  DO p = 1 , N
                     IF ( ((p>q) .AND. (DABS(A(p,q))<=temp1)) .OR.      &
     &                    (p<q) ) A(p,q) = DSIGN(temp1,A(p,q))
                  ENDDO
               ENDDO
            ELSE
               CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,A(1,2),Lda)
            ENDIF
!
!            .. second preconditioning using the QR factorization
!
            CALL DGEQRF(N,nr,A,Lda,Work,Work(N+1),Lwork-N,ierr)
!
!           .. and transpose upper to lower triangular
            DO p = 1 , nr - 1
               CALL DCOPY(nr-p,A(p,p+1),Lda,A(p+1,p),1)
            ENDDO
!
         ENDIF
!
!           Row-cyclic Jacobi SVD algorithm with column pivoting
!
!           .. again some perturbation (a "background noise") is added
!           to drown denormals
         IF ( l2pert ) THEN
!              XSC = DSQRT(SMALL)
            xsc = epsln/DBLE(N)
            DO q = 1 , nr
               temp1 = xsc*DABS(A(q,q))
               DO p = 1 , nr
                  IF ( ((p>q) .AND. (DABS(A(p,q))<=temp1)) .OR. (p<q) ) &
     &                 A(p,q) = DSIGN(temp1,A(p,q))
               ENDDO
            ENDDO
         ELSE
            CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,A(1,2),Lda)
         ENDIF
!
!           .. and one-sided Jacobi rotations are started on a lower
!           triangular matrix (plus perturbation which is ignored in
!           the part which destroys triangular form (confusing?!))
!
         CALL DGESVJ('L','NoU','NoV',nr,nr,A,Lda,Sva,N,V,Ldv,Work,Lwork,&
     &               Info)
!
         scalem = Work(1)
         numrank = IDNINT(Work(2))
!
!
      ELSEIF ( rsvec .AND. (.NOT.lsvec) ) THEN
!
!        -> Singular Values and Right Singular Vectors <-
!
         IF ( almort ) THEN
!
!           .. in this case NR equals N
            DO p = 1 , nr
               CALL DCOPY(N-p+1,A(p,p),Lda,V(p,p),1)
            ENDDO
            CALL DLASET('Upper',nr-1,nr-1,ZERO,ZERO,V(1,2),Ldv)
!
            CALL DGESVJ('L','U','N',N,nr,V,Ldv,Sva,nr,A,Lda,Work,Lwork, &
     &                  Info)
            scalem = Work(1)
            numrank = IDNINT(Work(2))
 
         ELSE
!
!        .. two more QR factorizations ( one QRF is not enough, two require
!        accumulated product of Jacobi rotations, three are perfect )
!
            CALL DLASET('Lower',nr-1,nr-1,ZERO,ZERO,A(2,1),Lda)
            CALL DGELQF(nr,N,A,Lda,Work,Work(N+1),Lwork-N,ierr)
            CALL DLACPY('Lower',nr,nr,A,Lda,V,Ldv)
            CALL DLASET('Upper',nr-1,nr-1,ZERO,ZERO,V(1,2),Ldv)
            CALL DGEQRF(nr,nr,V,Ldv,Work(N+1),Work(2*N+1),Lwork-2*N,    &
     &                  ierr)
            DO p = 1 , nr
               CALL DCOPY(nr-p+1,V(p,p),Ldv,V(p,p),1)
            ENDDO
            CALL DLASET('Upper',nr-1,nr-1,ZERO,ZERO,V(1,2),Ldv)
!
            CALL DGESVJ('Lower','U','N',nr,nr,V,Ldv,Sva,nr,U,Ldu,       &
     &                  Work(N+1),Lwork,Info)
            scalem = Work(N+1)
            numrank = IDNINT(Work(N+2))
            IF ( nr<N ) THEN
               CALL DLASET('A',N-nr,nr,ZERO,ZERO,V(nr+1,1),Ldv)
               CALL DLASET('A',nr,N-nr,ZERO,ZERO,V(1,nr+1),Ldv)
               CALL DLASET('A',N-nr,N-nr,ZERO,ONE,V(nr+1,nr+1),Ldv)
            ENDIF
!
            CALL DORMLQ('Left','Transpose',N,N,nr,A,Lda,Work,V,Ldv,     &
     &                  Work(N+1),Lwork-N,ierr)
!
         ENDIF
!
         DO p = 1 , N
            CALL DCOPY(N,V(p,1),Ldv,A(Iwork(p),1),Lda)
         ENDDO
         CALL DLACPY('All',N,N,A,Lda,V,Ldv)
!
         IF ( transp ) CALL DLACPY('All',N,N,V,Ldv,U,Ldu)
!
      ELSEIF ( lsvec .AND. (.NOT.rsvec) ) THEN
!
!        .. Singular Values and Left Singular Vectors                 ..
!
!        .. second preconditioning step to avoid need to accumulate
!        Jacobi rotations in the Jacobi iterations.
         DO p = 1 , nr
            CALL DCOPY(N-p+1,A(p,p),Lda,U(p,p),1)
         ENDDO
         CALL DLASET('Upper',nr-1,nr-1,ZERO,ZERO,U(1,2),Ldu)
!
         CALL DGEQRF(N,nr,U,Ldu,Work(N+1),Work(2*N+1),Lwork-2*N,ierr)
!
         DO p = 1 , nr - 1
            CALL DCOPY(nr-p,U(p,p+1),Ldu,U(p+1,p),1)
         ENDDO
         CALL DLASET('Upper',nr-1,nr-1,ZERO,ZERO,U(1,2),Ldu)
!
         CALL DGESVJ('Lower','U','N',nr,nr,U,Ldu,Sva,nr,A,Lda,Work(N+1),&
     &               Lwork-N,Info)
         scalem = Work(N+1)
         numrank = IDNINT(Work(N+2))
!
         IF ( nr<M ) THEN
            CALL DLASET('A',M-nr,nr,ZERO,ZERO,U(nr+1,1),Ldu)
            IF ( nr<n1 ) THEN
               CALL DLASET('A',nr,n1-nr,ZERO,ZERO,U(1,nr+1),Ldu)
               CALL DLASET('A',M-nr,n1-nr,ZERO,ONE,U(nr+1,nr+1),Ldu)
            ENDIF
         ENDIF
!
         CALL DORMQR('Left','No Tr',M,n1,N,A,Lda,Work,U,Ldu,Work(N+1),  &
     &               Lwork-N,ierr)
!
         IF ( rowpiv ) CALL DLASWP(n1,U,Ldu,1,M-1,Iwork(2*N+1),-1)
!
         DO p = 1 , n1
            xsc = ONE/DNRM2(M,U(1,p),1)
            CALL DSCAL(M,xsc,U(1,p),1)
         ENDDO
!
         IF ( transp ) CALL DLACPY('All',N,N,U,Ldu,V,Ldv)
!
      ELSE
!
!        .. Full SVD ..
!
         IF ( jracc ) THEN
!
!        This branch deploys a preconditioned Jacobi SVD with explicitly
!        accumulated rotations. It is included as optional, mainly for
!        experimental purposes. It does perform well, and can also be used.
!        In this implementation, this branch will be automatically activated
!        if the  condition number sigma_max(A) / sigma_min(A) is predicted
!        to be greater than the overflow threshold. This is because the
!        a posteriori computation of the singular vectors assumes robust
!        implementation of BLAS and some LAPACK procedures, capable of working
!        in presence of extreme values. Since that is not always the case, ...
!
            DO p = 1 , nr
               CALL DCOPY(N-p+1,A(p,p),Lda,V(p,p),1)
            ENDDO
!
            IF ( l2pert ) THEN
               xsc = DSQRT(small/epsln)
               DO q = 1 , nr
                  temp1 = xsc*DABS(V(q,q))
                  DO p = 1 , N
                     IF ( (p>q) .AND. (DABS(V(p,q))<=temp1) .OR. (p<q) )&
     &                    V(p,q) = DSIGN(temp1,V(p,q))
                     IF ( p<q ) V(p,q) = -V(p,q)
                  ENDDO
               ENDDO
            ELSE
               CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,V(1,2),Ldv)
            ENDIF
 
            CALL DGEQRF(N,nr,V,Ldv,Work(N+1),Work(2*N+1),Lwork-2*N,ierr)
            CALL DLACPY('L',N,nr,V,Ldv,Work(2*N+1),N)
!
            DO p = 1 , nr
               CALL DCOPY(nr-p+1,V(p,p),Ldv,U(p,p),1)
            ENDDO
 
            IF ( l2pert ) THEN
               xsc = DSQRT(small/epsln)
               DO q = 2 , nr
                  DO p = 1 , q - 1
                     temp1 = xsc*MIN(DABS(U(p,p)),DABS(U(q,q)))
                     U(p,q) = -DSIGN(temp1,U(q,p))
                  ENDDO
               ENDDO
            ELSE
               CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,U(1,2),Ldu)
            ENDIF
 
            CALL DGESVJ('G','U','V',nr,nr,U,Ldu,Sva,N,V,Ldv,            &
     &                  Work(2*N+N*nr+1),Lwork-2*N-N*nr,Info)
            scalem = Work(2*N+N*nr+1)
            numrank = IDNINT(Work(2*N+N*nr+2))
 
            IF ( nr<N ) THEN
               CALL DLASET('A',N-nr,nr,ZERO,ZERO,V(nr+1,1),Ldv)
               CALL DLASET('A',nr,N-nr,ZERO,ZERO,V(1,nr+1),Ldv)
               CALL DLASET('A',N-nr,N-nr,ZERO,ONE,V(nr+1,nr+1),Ldv)
            ENDIF
 
            CALL DORMQR('L','N',N,N,nr,Work(2*N+1),N,Work(N+1),V,Ldv,   &
     &                  Work(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,ierr)
!
!           Permute the rows of V using the (column) permutation from the
!           first QRF. Also, scale the columns to make them unit in
!           Euclidean norm. This applies to all cases.
!
            temp1 = DSQRT(DBLE(N))*epsln
            DO q = 1 , N
               DO p = 1 , N
                  Work(2*N+N*nr+nr+Iwork(p)) = V(p,q)
               ENDDO
               DO p = 1 , N
                  V(p,q) = Work(2*N+N*nr+nr+p)
               ENDDO
               xsc = ONE/DNRM2(N,V(1,q),1)
               IF ( (xsc<(ONE-temp1)) .OR. (xsc>(ONE+temp1)) )          &
     &              CALL DSCAL(N,xsc,V(1,q),1)
            ENDDO
!
!           At this moment, V contains the right singular vectors of A.
!           Next, assemble the left singular vector matrix U (M x N).
!
            IF ( nr<M ) THEN
               CALL DLASET('A',M-nr,nr,ZERO,ZERO,U(nr+1,1),Ldu)
               IF ( nr<n1 ) THEN
                  CALL DLASET('A',nr,n1-nr,ZERO,ZERO,U(1,nr+1),Ldu)
                  CALL DLASET('A',M-nr,n1-nr,ZERO,ONE,U(nr+1,nr+1),Ldu)
               ENDIF
            ENDIF
!
            CALL DORMQR('Left','No Tr',M,n1,N,A,Lda,Work,U,Ldu,Work(N+1)&
     &                  ,Lwork-N,ierr)
!
            IF ( rowpiv ) CALL DLASWP(n1,U,Ldu,1,M-1,Iwork(2*N+1),-1)
!
         ELSEIF ( .NOT.almort ) THEN
!
!           Second Preconditioning Step (QRF [with pivoting])
!           Note that the composition of TRANSPOSE, QRF and TRANSPOSE is
!           equivalent to an LQF CALL. Since in many libraries the QRF
!           seems to be better optimized than the LQF, we do explicit
!           transpose and use the QRF. This is subject to changes in an
!           optimized implementation of DGEJSV.
!
            DO p = 1 , nr
               CALL DCOPY(N-p+1,A(p,p),Lda,V(p,p),1)
            ENDDO
!
!           .. the following two loops perturb small entries to avoid
!           denormals in the second QR factorization, where they are
!           as good as zeros. This is done to avoid painfully slow
!           computation with denormals. The relative size of the perturbation
!           is a parameter that can be changed by the implementer.
!           This perturbation device will be obsolete on machines with
!           properly implemented arithmetic.
!           To switch it off, set L2PERT=.FALSE. To remove it from  the
!           code, remove the action under L2PERT=.TRUE., leave the ELSE part.
!           The following two loops should be blocked and fused with the
!           transposed copy above.
!
            IF ( l2pert ) THEN
               xsc = DSQRT(small)
               DO q = 1 , nr
                  temp1 = xsc*DABS(V(q,q))
                  DO p = 1 , N
                     IF ( (p>q) .AND. (DABS(V(p,q))<=temp1) .OR. (p<q) )&
     &                    V(p,q) = DSIGN(temp1,V(p,q))
                     IF ( p<q ) V(p,q) = -V(p,q)
                  ENDDO
               ENDDO
            ELSE
               CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,V(1,2),Ldv)
            ENDIF
!
!           Estimate the row scaled condition number of R1
!           (If R1 is rectangular, N > NR, then the condition number
!           of the leading NR x NR submatrix is estimated.)
!
            CALL DLACPY('L',nr,nr,V,Ldv,Work(2*N+1),nr)
            DO p = 1 , nr
               temp1 = DNRM2(nr-p+1,Work(2*N+(p-1)*nr+p),1)
               CALL DSCAL(nr-p+1,ONE/temp1,Work(2*N+(p-1)*nr+p),1)
            ENDDO
            CALL DPOCON('Lower',nr,Work(2*N+1),nr,ONE,temp1,            &
     &                  Work(2*N+nr*nr+1),Iwork(M+2*N+1),ierr)
            condr1 = ONE/DSQRT(temp1)
!           .. here need a second opinion on the condition number
!           .. then assume worst case scenario
!           R1 is OK for inverse <=> CONDR1 .LT. DBLE(N)
!           more conservative    <=> CONDR1 .LT. DSQRT(DBLE(N))
!
            cond_ok = DSQRT(DBLE(nr))
![TP]       COND_OK is a tuning parameter.
 
            IF ( condr1<cond_ok ) THEN
!              .. the second QRF without pivoting. Note: in an optimized
!              implementation, this QRF should be implemented as the QRF
!              of a lower triangular matrix.
!              R1^t = Q2 * R2
               CALL DGEQRF(N,nr,V,Ldv,Work(N+1),Work(2*N+1),Lwork-2*N,  &
     &                     ierr)
!
               IF ( l2pert ) THEN
                  xsc = DSQRT(small)/epsln
                  DO p = 2 , nr
                     DO q = 1 , p - 1
                        temp1 = xsc*MIN(DABS(V(p,p)),DABS(V(q,q)))
                        IF ( DABS(V(q,p))<=temp1 ) V(q,p)               &
     &                       = DSIGN(temp1,V(q,p))
                     ENDDO
                  ENDDO
               ENDIF
!
               IF ( nr/=N ) CALL DLACPY('A',N,nr,V,Ldv,Work(2*N+1),N)
!              .. save ...
!
!           .. this transposed copy should be better than naive
               DO p = 1 , nr - 1
                  CALL DCOPY(nr-p,V(p,p+1),Ldv,V(p+1,p),1)
               ENDDO
!
               condr2 = condr1
!
            ELSE
!
!              .. ill-conditioned case: second QRF with pivoting
!              Note that windowed pivoting would be equally good
!              numerically, and more run-time efficient. So, in
!              an optimal implementation, the next call to DGEQP3
!              should be replaced with eg. CALL SGEQPX (ACM TOMS #782)
!              with properly (carefully) chosen parameters.
!
!              R1^t * P2 = Q2 * R2
               DO p = 1 , nr
                  Iwork(N+p) = 0
               ENDDO
               CALL DGEQP3(N,nr,V,Ldv,Iwork(N+1),Work(N+1),Work(2*N+1), &
     &                     Lwork-2*N,ierr)
!*               CALL DGEQRF( N, NR, V, LDV, WORK(N+1), WORK(2*N+1),
!*     $              LWORK-2*N, IERR )
               IF ( l2pert ) THEN
                  xsc = DSQRT(small)
                  DO p = 2 , nr
                     DO q = 1 , p - 1
                        temp1 = xsc*MIN(DABS(V(p,p)),DABS(V(q,q)))
                        IF ( DABS(V(q,p))<=temp1 ) V(q,p)               &
     &                       = DSIGN(temp1,V(q,p))
                     ENDDO
                  ENDDO
               ENDIF
!
               CALL DLACPY('A',N,nr,V,Ldv,Work(2*N+1),N)
!
               IF ( l2pert ) THEN
                  xsc = DSQRT(small)
                  DO p = 2 , nr
                     DO q = 1 , p - 1
                        temp1 = xsc*MIN(DABS(V(p,p)),DABS(V(q,q)))
                        V(p,q) = -DSIGN(temp1,V(q,p))
                     ENDDO
                  ENDDO
               ELSE
                  CALL DLASET('L',nr-1,nr-1,ZERO,ZERO,V(2,1),Ldv)
               ENDIF
!              Now, compute R2 = L3 * Q3, the LQ factorization.
               CALL DGELQF(nr,nr,V,Ldv,Work(2*N+N*nr+1),                &
     &                     Work(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,ierr)
!              .. and estimate the condition number
               CALL DLACPY('L',nr,nr,V,Ldv,Work(2*N+N*nr+nr+1),nr)
               DO p = 1 , nr
                  temp1 = DNRM2(p,Work(2*N+N*nr+nr+p),nr)
                  CALL DSCAL(p,ONE/temp1,Work(2*N+N*nr+nr+p),nr)
               ENDDO
               CALL DPOCON('L',nr,Work(2*N+N*nr+nr+1),nr,ONE,temp1,     &
     &                     Work(2*N+N*nr+nr+nr*nr+1),Iwork(M+2*N+1),    &
     &                     ierr)
               condr2 = ONE/DSQRT(temp1)
!
!                 .. save the Householder vectors used for Q3
!                 (this overwrites the copy of R2, as it will not be
!                 needed in this branch, but it does not overwritte the
!                 Huseholder vectors of Q2.).
!                 .. and the rest of the information on Q3 is in
!                 WORK(2*N+N*NR+1:2*N+N*NR+N)
               IF ( condr2>=cond_ok )                                   &
     &              CALL DLACPY('U',nr,nr,V,Ldv,Work(2*N+1),N)
!
            ENDIF
!
            IF ( l2pert ) THEN
               xsc = DSQRT(small)
               DO q = 2 , nr
                  temp1 = xsc*V(q,q)
                  DO p = 1 , q - 1
!                    V(p,q) = - DSIGN( TEMP1, V(q,p) )
                     V(p,q) = -DSIGN(temp1,V(p,q))
                  ENDDO
               ENDDO
            ELSE
               CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,V(1,2),Ldv)
            ENDIF
!
!        Second preconditioning finished; continue with Jacobi SVD
!        The input matrix is lower trinagular.
!
!        Recover the right singular vectors as solution of a well
!        conditioned triangular matrix equation.
!
            IF ( condr1<cond_ok ) THEN
!
               CALL DGESVJ('L','U','N',nr,nr,V,Ldv,Sva,nr,U,Ldu,        &
     &                     Work(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,Info)
               scalem = Work(2*N+N*nr+nr+1)
               numrank = IDNINT(Work(2*N+N*nr+nr+2))
               DO p = 1 , nr
                  CALL DCOPY(nr,V(1,p),1,U(1,p),1)
                  CALL DSCAL(nr,Sva(p),V(1,p),1)
               ENDDO
 
!        .. pick the right matrix equation and solve it
!
               IF ( nr==N ) THEN
! :))             .. best case, R1 is inverted. The solution of this matrix
!                 equation is Q2*V2 = the product of the Jacobi rotations
!                 used in DGESVJ, premultiplied with the orthogonal matrix
!                 from the second QR factorization.
                  CALL DTRSM('L','U','N','N',nr,nr,ONE,A,Lda,V,Ldv)
               ELSE
!                 .. R1 is well conditioned, but non-square. Transpose(R2)
!                 is inverted to get the product of the Jacobi rotations
!                 used in DGESVJ. The Q-factor from the second QR
!                 factorization is then built in explicitly.
                  CALL DTRSM('L','U','T','N',nr,nr,ONE,Work(2*N+1),N,V, &
     &                       Ldv)
                  IF ( nr<N ) THEN
                     CALL DLASET('A',N-nr,nr,ZERO,ZERO,V(nr+1,1),Ldv)
                     CALL DLASET('A',nr,N-nr,ZERO,ZERO,V(1,nr+1),Ldv)
                     CALL DLASET('A',N-nr,N-nr,ZERO,ONE,V(nr+1,nr+1),   &
     &                           Ldv)
                  ENDIF
                  CALL DORMQR('L','N',N,N,nr,Work(2*N+1),N,Work(N+1),V, &
     &                        Ldv,Work(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,&
     &                        ierr)
               ENDIF
!
            ELSEIF ( condr2<cond_ok ) THEN
!
! :)           .. the input matrix A is very likely a relative of
!              the Kahan matrix :)
!              The matrix R2 is inverted. The solution of the matrix equation
!              is Q3^T*V3 = the product of the Jacobi rotations (appplied to
!              the lower triangular L3 from the LQ factorization of
!              R2=L3*Q3), pre-multiplied with the transposed Q3.
               CALL DGESVJ('L','U','N',nr,nr,V,Ldv,Sva,nr,U,Ldu,        &
     &                     Work(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,Info)
               scalem = Work(2*N+N*nr+nr+1)
               numrank = IDNINT(Work(2*N+N*nr+nr+2))
               DO p = 1 , nr
                  CALL DCOPY(nr,V(1,p),1,U(1,p),1)
                  CALL DSCAL(nr,Sva(p),U(1,p),1)
               ENDDO
               CALL DTRSM('L','U','N','N',nr,nr,ONE,Work(2*N+1),N,U,Ldu)
!              .. apply the permutation from the second QR factorization
               DO q = 1 , nr
                  DO p = 1 , nr
                     Work(2*N+N*nr+nr+Iwork(N+p)) = U(p,q)
                  ENDDO
                  DO p = 1 , nr
                     U(p,q) = Work(2*N+N*nr+nr+p)
                  ENDDO
               ENDDO
               IF ( nr<N ) THEN
                  CALL DLASET('A',N-nr,nr,ZERO,ZERO,V(nr+1,1),Ldv)
                  CALL DLASET('A',nr,N-nr,ZERO,ZERO,V(1,nr+1),Ldv)
                  CALL DLASET('A',N-nr,N-nr,ZERO,ONE,V(nr+1,nr+1),Ldv)
               ENDIF
               CALL DORMQR('L','N',N,N,nr,Work(2*N+1),N,Work(N+1),V,Ldv,&
     &                     Work(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,ierr)
            ELSE
!              Last line of defense.
! #:(          This is a rather pathological case: no scaled condition
!              improvement after two pivoted QR factorizations. Other
!              possibility is that the rank revealing QR factorization
!              or the condition estimator has failed, or the COND_OK
!              is set very close to ONE (which is unnecessary). Normally,
!              this branch should never be executed, but in rare cases of
!              failure of the RRQR or condition estimator, the last line of
!              defense ensures that DGEJSV completes the task.
!              Compute the full SVD of L3 using DGESVJ with explicit
!              accumulation of Jacobi rotations.
               CALL DGESVJ('L','U','V',nr,nr,V,Ldv,Sva,nr,U,Ldu,        &
     &                     Work(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,Info)
               scalem = Work(2*N+N*nr+nr+1)
               numrank = IDNINT(Work(2*N+N*nr+nr+2))
               IF ( nr<N ) THEN
                  CALL DLASET('A',N-nr,nr,ZERO,ZERO,V(nr+1,1),Ldv)
                  CALL DLASET('A',nr,N-nr,ZERO,ZERO,V(1,nr+1),Ldv)
                  CALL DLASET('A',N-nr,N-nr,ZERO,ONE,V(nr+1,nr+1),Ldv)
               ENDIF
               CALL DORMQR('L','N',N,N,nr,Work(2*N+1),N,Work(N+1),V,Ldv,&
     &                     Work(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,ierr)
!
               CALL DORMLQ('L','T',nr,nr,nr,Work(2*N+1),N,              &
     &                     Work(2*N+N*nr+1),U,Ldu,Work(2*N+N*nr+nr+1),  &
     &                     Lwork-2*N-N*nr-nr,ierr)
               DO q = 1 , nr
                  DO p = 1 , nr
                     Work(2*N+N*nr+nr+Iwork(N+p)) = U(p,q)
                  ENDDO
                  DO p = 1 , nr
                     U(p,q) = Work(2*N+N*nr+nr+p)
                  ENDDO
               ENDDO
!
            ENDIF
!
!           Permute the rows of V using the (column) permutation from the
!           first QRF. Also, scale the columns to make them unit in
!           Euclidean norm. This applies to all cases.
!
            temp1 = DSQRT(DBLE(N))*epsln
            DO q = 1 , N
               DO p = 1 , N
                  Work(2*N+N*nr+nr+Iwork(p)) = V(p,q)
               ENDDO
               DO p = 1 , N
                  V(p,q) = Work(2*N+N*nr+nr+p)
               ENDDO
               xsc = ONE/DNRM2(N,V(1,q),1)
               IF ( (xsc<(ONE-temp1)) .OR. (xsc>(ONE+temp1)) )          &
     &              CALL DSCAL(N,xsc,V(1,q),1)
            ENDDO
!           At this moment, V contains the right singular vectors of A.
!           Next, assemble the left singular vector matrix U (M x N).
            IF ( nr<M ) THEN
               CALL DLASET('A',M-nr,nr,ZERO,ZERO,U(nr+1,1),Ldu)
               IF ( nr<n1 ) THEN
                  CALL DLASET('A',nr,n1-nr,ZERO,ZERO,U(1,nr+1),Ldu)
                  CALL DLASET('A',M-nr,n1-nr,ZERO,ONE,U(nr+1,nr+1),Ldu)
               ENDIF
            ENDIF
!
!           The Q matrix from the first QRF is built into the left singular
!           matrix U. This applies to all cases.
!
            CALL DORMQR('Left','No_Tr',M,n1,N,A,Lda,Work,U,Ldu,Work(N+1)&
     &                  ,Lwork-N,ierr)
 
!           The columns of U are normalized. The cost is O(M*N) flops.
            temp1 = DSQRT(DBLE(M))*epsln
            DO p = 1 , nr
               xsc = ONE/DNRM2(M,U(1,p),1)
               IF ( (xsc<(ONE-temp1)) .OR. (xsc>(ONE+temp1)) )          &
     &              CALL DSCAL(M,xsc,U(1,p),1)
            ENDDO
!
!           If the initial QRF is computed with row pivoting, the left
!           singular vectors must be adjusted.
!
            IF ( rowpiv ) CALL DLASWP(n1,U,Ldu,1,M-1,Iwork(2*N+1),-1)
!
         ELSE
!
!        .. the initial matrix A has almost orthogonal columns and
!        the second QRF is not needed
!
            CALL DLACPY('Upper',N,N,A,Lda,Work(N+1),N)
            IF ( l2pert ) THEN
               xsc = DSQRT(small)
               DO p = 2 , N
                  temp1 = xsc*Work(N+(p-1)*N+p)
                  DO q = 1 , p - 1
                     Work(N+(q-1)*N+p) = -DSIGN(temp1,Work(N+(p-1)*N+q))
                  ENDDO
               ENDDO
            ELSE
               CALL DLASET('Lower',N-1,N-1,ZERO,ZERO,Work(N+2),N)
            ENDIF
!
            CALL DGESVJ('Upper','U','N',N,N,Work(N+1),N,Sva,N,U,Ldu,    &
     &                  Work(N+N*N+1),Lwork-N-N*N,Info)
!
            scalem = Work(N+N*N+1)
            numrank = IDNINT(Work(N+N*N+2))
            DO p = 1 , N
               CALL DCOPY(N,Work(N+(p-1)*N+1),1,U(1,p),1)
               CALL DSCAL(N,Sva(p),Work(N+(p-1)*N+1),1)
            ENDDO
!
            CALL DTRSM('Left','Upper','NoTrans','No UD',N,N,ONE,A,Lda,  &
     &                 Work(N+1),N)
            DO p = 1 , N
               CALL DCOPY(N,Work(N+p),N,V(Iwork(p),1),Ldv)
            ENDDO
            temp1 = DSQRT(DBLE(N))*epsln
            DO p = 1 , N
               xsc = ONE/DNRM2(N,V(1,p),1)
               IF ( (xsc<(ONE-temp1)) .OR. (xsc>(ONE+temp1)) )          &
     &              CALL DSCAL(N,xsc,V(1,p),1)
            ENDDO
!
!           Assemble the left singular vector matrix U (M x N).
!
            IF ( N<M ) THEN
               CALL DLASET('A',M-N,N,ZERO,ZERO,U(N+1,1),Ldu)
               IF ( N<n1 ) THEN
                  CALL DLASET('A',N,n1-N,ZERO,ZERO,U(1,N+1),Ldu)
                  CALL DLASET('A',M-N,n1-N,ZERO,ONE,U(N+1,N+1),Ldu)
               ENDIF
            ENDIF
            CALL DORMQR('Left','No Tr',M,n1,N,A,Lda,Work,U,Ldu,Work(N+1)&
     &                  ,Lwork-N,ierr)
            temp1 = DSQRT(DBLE(M))*epsln
            DO p = 1 , n1
               xsc = ONE/DNRM2(M,U(1,p),1)
               IF ( (xsc<(ONE-temp1)) .OR. (xsc>(ONE+temp1)) )          &
     &              CALL DSCAL(M,xsc,U(1,p),1)
            ENDDO
!
            IF ( rowpiv ) CALL DLASWP(n1,U,Ldu,1,M-1,Iwork(2*N+1),-1)
!
!
         ENDIF
!
!
!        end of the  >> almost orthogonal case <<  in the full SVD
!
         IF ( transp ) THEN
!           .. swap U and V because the procedure worked on A^t
            DO p = 1 , N
               CALL DSWAP(N,U(1,p),1,V(1,p),1)
            ENDDO
         ENDIF
!
      ENDIF
!     end of the full SVD
!
!     Undo scaling, if necessary (and possible)
!
      IF ( uscal2<=(big/Sva(1))*uscal1 ) THEN
         CALL DLASCL('G',0,0,uscal1,uscal2,nr,1,Sva,N,ierr)
         uscal1 = ONE
         uscal2 = ONE
      ENDIF
!
      IF ( nr<N ) THEN
         DO p = nr + 1 , N
            Sva(p) = ZERO
         ENDDO
      ENDIF
!
      Work(1) = uscal2*scalem
      Work(2) = uscal1
      IF ( errest ) Work(3) = sconda
      IF ( lsvec .AND. rsvec ) THEN
         Work(4) = condr1
         Work(5) = condr2
      ENDIF
      IF ( l2tran ) THEN
         Work(6) = entra
         Work(7) = entrat
      ENDIF
!
      Iwork(1) = nr
      Iwork(2) = numrank
      Iwork(3) = warning
!
!     ..
!     .. END OF DGEJSV
!     ..
      END SUBROUTINE DGEJSV
