!*==zgejsv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZGEJSV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGEJSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgejsv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgejsv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgejsv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!     SUBROUTINE ZGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP,
!                         M, N, A, LDA, SVA, U, LDU, V, LDV,
!                         CWORK, LWORK, RWORK, LRWORK, IWORK, INFO )
!
!     .. Scalar Arguments ..
!     IMPLICIT    NONE
!     INTEGER     INFO, LDA, LDU, LDV, LWORK, M, N
!     ..
!     .. Array Arguments ..
!     COMPLEX*16     A( LDA, * ),  U( LDU, * ), V( LDV, * ), CWORK( LWORK )
!     DOUBLE PRECISION   SVA( N ), RWORK( LRWORK )
!     INTEGER     IWORK( * )
!     CHARACTER*1 JOBA, JOBP, JOBR, JOBT, JOBU, JOBV
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGEJSV computes the singular value decomposition (SVD) of a complex M-by-N
!> matrix [A], where M >= N. The SVD of [A] is written as
!>
!>              [A] = [U] * [SIGMA] * [V]^*,
!>
!> where [SIGMA] is an N-by-N (M-by-N) matrix which is zero except for its N
!> diagonal elements, [U] is an M-by-N (or M-by-M) unitary matrix, and
!> [V] is an N-by-N unitary matrix. The diagonal elements of [SIGMA] are
!> the singular values of [A]. The columns of [U] and [V] are the left and
!> the right singular vectors of [A], respectively. The matrices [U] and [V]
!> are computed and stored in the arrays U and V, respectively. The diagonal
!> of [SIGMA] is computed and stored in the array SVA.
!> \endverbatim
!>
!>  Arguments:
!>  ==========
!>
!> \param[in] JOBA
!> \verbatim
!>          JOBA is CHARACTER*1
!>         Specifies the level of accuracy:
!>       = 'C': This option works well (high relative accuracy) if A = B * D,
!>              with well-conditioned B and arbitrary diagonal matrix D.
!>              The accuracy cannot be spoiled by COLUMN scaling. The
!>              accuracy of the computed output depends on the condition of
!>              B, and the procedure aims at the best theoretical accuracy.
!>              The relative error max_{i=1:N}|d sigma_i| / sigma_i is
!>              bounded by f(M,N)*epsilon* cond(B), independent of D.
!>              The input matrix is preprocessed with the QRF with column
!>              pivoting. This initial preprocessing and preconditioning by
!>              a rank revealing QR factorization is common for all values of
!>              JOBA. Additional actions are specified as follows:
!>       = 'E': Computation as with 'C' with an additional estimate of the
!>              condition number of B. It provides a realistic error bound.
!>       = 'F': If A = D1 * C * D2 with ill-conditioned diagonal scalings
!>              D1, D2, and well-conditioned matrix C, this option gives
!>              higher accuracy than the 'C' option. If the structure of the
!>              input matrix is not known, and relative accuracy is
!>              desirable, then this option is advisable. The input matrix A
!>              is preprocessed with QR factorization with FULL (row and
!>              column) pivoting.
!>       = 'G': Computation as with 'F' with an additional estimate of the
!>              condition number of B, where A=B*D. If A has heavily weighted
!>              rows, then using this condition number gives too pessimistic
!>              error bound.
!>       = 'A': Small singular values are not well determined by the data
!>              and are considered as noisy; the matrix is treated as
!>              numerically rank deficient. The error in the computed
!>              singular values is bounded by f(m,n)*epsilon*||A||.
!>              The computed SVD A = U * S * V^* restores A up to
!>              f(m,n)*epsilon*||A||.
!>              This gives the procedure the licence to discard (set to zero)
!>              all singular values below N*epsilon*||A||.
!>       = 'R': Similar as in 'A'. Rank revealing property of the initial
!>              QR factorization is used do reveal (using triangular factor)
!>              a gap sigma_{r+1} < epsilon * sigma_r in which case the
!>              numerical RANK is declared to be r. The SVD is computed with
!>              absolute error bounds, but more accurately than with 'A'.
!> \endverbatim
!>
!> \param[in] JOBU
!> \verbatim
!>          JOBU is CHARACTER*1
!>         Specifies whether to compute the columns of U:
!>       = 'U': N columns of U are returned in the array U.
!>       = 'F': full set of M left sing. vectors is returned in the array U.
!>       = 'W': U may be used as workspace of length M*N. See the description
!>              of U.
!>       = 'N': U is not computed.
!> \endverbatim
!>
!> \param[in] JOBV
!> \verbatim
!>          JOBV is CHARACTER*1
!>         Specifies whether to compute the matrix V:
!>       = 'V': N columns of V are returned in the array V; Jacobi rotations
!>              are not explicitly accumulated.
!>       = 'J': N columns of V are returned in the array V, but they are
!>              computed as the product of Jacobi rotations, if JOBT = 'N'.
!>       = 'W': V may be used as workspace of length N*N. See the description
!>              of V.
!>       = 'N': V is not computed.
!> \endverbatim
!>
!> \param[in] JOBR
!> \verbatim
!>          JOBR is CHARACTER*1
!>         Specifies the RANGE for the singular values. Issues the licence to
!>         set to zero small positive singular values if they are outside
!>         specified range. If A .NE. 0 is scaled so that the largest singular
!>         value of c*A is around SQRT(BIG), BIG=DLAMCH('O'), then JOBR issues
!>         the licence to kill columns of A whose norm in c*A is less than
!>         SQRT(SFMIN) (for JOBR = 'R'), or less than SMALL=SFMIN/EPSLN,
!>         where SFMIN=DLAMCH('S'), EPSLN=DLAMCH('E').
!>       = 'N': Do not kill small columns of c*A. This option assumes that
!>              BLAS and QR factorizations and triangular solvers are
!>              implemented to work in that range. If the condition of A
!>              is greater than BIG, use ZGESVJ.
!>       = 'R': RESTRICTED range for sigma(c*A) is [SQRT(SFMIN), SQRT(BIG)]
!>              (roughly, as described above). This option is recommended.
!>                                             ===========================
!>         For computing the singular values in the FULL range [SFMIN,BIG]
!>         use ZGESVJ.
!> \endverbatim
!>
!> \param[in] JOBT
!> \verbatim
!>          JOBT is CHARACTER*1
!>         If the matrix is square then the procedure may determine to use
!>         transposed A if A^* seems to be better with respect to convergence.
!>         If the matrix is not square, JOBT is ignored.
!>         The decision is based on two values of entropy over the adjoint
!>         orbit of A^* * A. See the descriptions of WORK(6) and WORK(7).
!>       = 'T': transpose if entropy test indicates possibly faster
!>         convergence of Jacobi process if A^* is taken as input. If A is
!>         replaced with A^*, then the row pivoting is included automatically.
!>       = 'N': do not speculate.
!>         The option 'T' can be used to compute only the singular values, or
!>         the full SVD (U, SIGMA and V). For only one set of singular vectors
!>         (U or V), the caller should provide both U and V, as one of the
!>         matrices is used as workspace if the matrix A is transposed.
!>         The implementer can easily remove this constraint and make the
!>         code more complicated. See the descriptions of U and V.
!>         In general, this option is considered experimental, and 'N'; should
!>         be preferred. This is subject to changes in the future.
!> \endverbatim
!>
!> \param[in] JOBP
!> \verbatim
!>          JOBP is CHARACTER*1
!>         Issues the licence to introduce structured perturbations to drown
!>         denormalized numbers. This licence should be active if the
!>         denormals are poorly implemented, causing slow computation,
!>         especially in cases of fast convergence (!). For details see [1,2].
!>         For the sake of simplicity, this perturbations are included only
!>         when the full SVD or only the singular values are requested. The
!>         implementer/user can easily add the perturbation for the cases of
!>         computing one set of singular vectors.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
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
!>          U is COMPLEX*16 array, dimension ( LDU, N )
!>          If JOBU = 'U', then U contains on exit the M-by-N matrix of
!>                         the left singular vectors.
!>          If JOBU = 'F', then U contains on exit the M-by-M matrix of
!>                         the left singular vectors, including an ONB
!>                         of the orthogonal complement of the Range(A).
!>          If JOBU = 'W'  .AND. (JOBV = 'V' .AND. JOBT = 'T' .AND. M = N),
!>                         then U is used as workspace if the procedure
!>                         replaces A with A^*. In that case, [V] is computed
!>                         in U as left singular vectors of A^* and then
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
!>          V is COMPLEX*16 array, dimension ( LDV, N )
!>          If JOBV = 'V', 'J' then V contains on exit the N-by-N matrix of
!>                         the right singular vectors;
!>          If JOBV = 'W', AND (JOBU = 'U' AND JOBT = 'T' AND M = N),
!>                         then V is used as workspace if the pprocedure
!>                         replaces A with A^*. In that case, [U] is computed
!>                         in V as right singular vectors of A^* and then
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
!> \param[out] CWORK
!> \verbatim
!>          CWORK is COMPLEX*16 array, dimension (MAX(2,LWORK))
!>          If the call to ZGEJSV is a workspace query (indicated by LWORK=-1 or
!>          LRWORK=-1), then on exit CWORK(1) contains the required length of
!>          CWORK for the job parameters used in the call.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          Length of CWORK to confirm proper allocation of workspace.
!>          LWORK depends on the job:
!>
!>          1. If only SIGMA is needed ( JOBU = 'N', JOBV = 'N' ) and
!>            1.1 .. no scaled condition estimate required (JOBA.NE.'E'.AND.JOBA.NE.'G'):
!>               LWORK >= 2*N+1. This is the minimal requirement.
!>               ->> For optimal performance (blocked code) the optimal value
!>               is LWORK >= N + (N+1)*NB. Here NB is the optimal
!>               block size for ZGEQP3 and ZGEQRF.
!>               In general, optimal LWORK is computed as
!>               LWORK >= max(N+LWORK(ZGEQP3),N+LWORK(ZGEQRF), LWORK(ZGESVJ)).
!>            1.2. .. an estimate of the scaled condition number of A is
!>               required (JOBA='E', or 'G'). In this case, LWORK the minimal
!>               requirement is LWORK >= N*N + 2*N.
!>               ->> For optimal performance (blocked code) the optimal value
!>               is LWORK >= max(N+(N+1)*NB, N*N+2*N)=N**2+2*N.
!>               In general, the optimal length LWORK is computed as
!>               LWORK >= max(N+LWORK(ZGEQP3),N+LWORK(ZGEQRF), LWORK(ZGESVJ),
!>                            N*N+LWORK(ZPOCON)).
!>          2. If SIGMA and the right singular vectors are needed (JOBV = 'V'),
!>             (JOBU = 'N')
!>            2.1   .. no scaled condition estimate requested (JOBE = 'N'):
!>            -> the minimal requirement is LWORK >= 3*N.
!>            -> For optimal performance,
!>               LWORK >= max(N+(N+1)*NB, 2*N+N*NB)=2*N+N*NB,
!>               where NB is the optimal block size for ZGEQP3, ZGEQRF, ZGELQ,
!>               ZUNMLQ. In general, the optimal length LWORK is computed as
!>               LWORK >= max(N+LWORK(ZGEQP3), N+LWORK(ZGESVJ),
!>                       N+LWORK(ZGELQF), 2*N+LWORK(ZGEQRF), N+LWORK(ZUNMLQ)).
!>            2.2 .. an estimate of the scaled condition number of A is
!>               required (JOBA='E', or 'G').
!>            -> the minimal requirement is LWORK >= 3*N.
!>            -> For optimal performance,
!>               LWORK >= max(N+(N+1)*NB, 2*N,2*N+N*NB)=2*N+N*NB,
!>               where NB is the optimal block size for ZGEQP3, ZGEQRF, ZGELQ,
!>               ZUNMLQ. In general, the optimal length LWORK is computed as
!>               LWORK >= max(N+LWORK(ZGEQP3), LWORK(ZPOCON), N+LWORK(ZGESVJ),
!>                       N+LWORK(ZGELQF), 2*N+LWORK(ZGEQRF), N+LWORK(ZUNMLQ)).
!>          3. If SIGMA and the left singular vectors are needed
!>            3.1  .. no scaled condition estimate requested (JOBE = 'N'):
!>            -> the minimal requirement is LWORK >= 3*N.
!>            -> For optimal performance:
!>               if JOBU = 'U' :: LWORK >= max(3*N, N+(N+1)*NB, 2*N+N*NB)=2*N+N*NB,
!>               where NB is the optimal block size for ZGEQP3, ZGEQRF, ZUNMQR.
!>               In general, the optimal length LWORK is computed as
!>               LWORK >= max(N+LWORK(ZGEQP3), 2*N+LWORK(ZGEQRF), N+LWORK(ZUNMQR)).
!>            3.2  .. an estimate of the scaled condition number of A is
!>               required (JOBA='E', or 'G').
!>            -> the minimal requirement is LWORK >= 3*N.
!>            -> For optimal performance:
!>               if JOBU = 'U' :: LWORK >= max(3*N, N+(N+1)*NB, 2*N+N*NB)=2*N+N*NB,
!>               where NB is the optimal block size for ZGEQP3, ZGEQRF, ZUNMQR.
!>               In general, the optimal length LWORK is computed as
!>               LWORK >= max(N+LWORK(ZGEQP3),N+LWORK(ZPOCON),
!>                        2*N+LWORK(ZGEQRF), N+LWORK(ZUNMQR)).
!>          4. If the full SVD is needed: (JOBU = 'U' or JOBU = 'F') and
!>            4.1. if JOBV = 'V'
!>               the minimal requirement is LWORK >= 5*N+2*N*N.
!>            4.2. if JOBV = 'J' the minimal requirement is
!>               LWORK >= 4*N+N*N.
!>            In both cases, the allocated CWORK can accommodate blocked runs
!>            of ZGEQP3, ZGEQRF, ZGELQF, SUNMQR, ZUNMLQ.
!>
!>          If the call to ZGEJSV is a workspace query (indicated by LWORK=-1 or
!>          LRWORK=-1), then on exit CWORK(1) contains the optimal and CWORK(2) contains the
!>          minimal length of CWORK for the job parameters used in the call.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (MAX(7,LWORK))
!>          On exit,
!>          RWORK(1) = Determines the scaling factor SCALE = RWORK(2) / RWORK(1)
!>                    such that SCALE*SVA(1:N) are the computed singular values
!>                    of A. (See the description of SVA().)
!>          RWORK(2) = See the description of RWORK(1).
!>          RWORK(3) = SCONDA is an estimate for the condition number of
!>                    column equilibrated A. (If JOBA = 'E' or 'G')
!>                    SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1).
!>                    It is computed using SPOCON. It holds
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
!>          RWORK(4) = an estimate of the scaled condition number of the
!>                    triangular factor in the first QR factorization.
!>          RWORK(5) = an estimate of the scaled condition number of the
!>                    triangular factor in the second QR factorization.
!>          The following two parameters are computed if JOBT = 'T'.
!>          They are provided for a developer/implementer who is familiar
!>          with the details of the method.
!>          RWORK(6) = the entropy of A^* * A :: this is the Shannon entropy
!>                    of diag(A^* * A) / Trace(A^* * A) taken as point in the
!>                    probability simplex.
!>          RWORK(7) = the entropy of A * A^*. (See the description of RWORK(6).)
!>          If the call to ZGEJSV is a workspace query (indicated by LWORK=-1 or
!>          LRWORK=-1), then on exit RWORK(1) contains the required length of
!>          RWORK for the job parameters used in the call.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          Length of RWORK to confirm proper allocation of workspace.
!>          LRWORK depends on the job:
!>
!>       1. If only the singular values are requested i.e. if
!>          LSAME(JOBU,'N') .AND. LSAME(JOBV,'N')
!>          then:
!>          1.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'),
!>               then: LRWORK = max( 7, 2 * M ).
!>          1.2. Otherwise, LRWORK  = max( 7,  N ).
!>       2. If singular values with the right singular vectors are requested
!>          i.e. if
!>          (LSAME(JOBV,'V').OR.LSAME(JOBV,'J')) .AND.
!>          .NOT.(LSAME(JOBU,'U').OR.LSAME(JOBU,'F'))
!>          then:
!>          2.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'),
!>          then LRWORK = max( 7, 2 * M ).
!>          2.2. Otherwise, LRWORK  = max( 7,  N ).
!>       3. If singular values with the left singular vectors are requested, i.e. if
!>          (LSAME(JOBU,'U').OR.LSAME(JOBU,'F')) .AND.
!>          .NOT.(LSAME(JOBV,'V').OR.LSAME(JOBV,'J'))
!>          then:
!>          3.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'),
!>          then LRWORK = max( 7, 2 * M ).
!>          3.2. Otherwise, LRWORK  = max( 7,  N ).
!>       4. If singular values with both the left and the right singular vectors
!>          are requested, i.e. if
!>          (LSAME(JOBU,'U').OR.LSAME(JOBU,'F')) .AND.
!>          (LSAME(JOBV,'V').OR.LSAME(JOBV,'J'))
!>          then:
!>          4.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'),
!>          then LRWORK = max( 7, 2 * M ).
!>          4.2. Otherwise, LRWORK  = max( 7, N ).
!>
!>          If, on entry, LRWORK = -1 or LWORK=-1, a workspace query is assumed and
!>          the length of RWORK is returned in RWORK(1).
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, of dimension at least 4, that further depends
!>          on the job:
!>
!>          1. If only the singular values are requested then:
!>             If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
!>             then the length of IWORK is N+M; otherwise the length of IWORK is N.
!>          2. If the singular values and the right singular vectors are requested then:
!>             If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
!>             then the length of IWORK is N+M; otherwise the length of IWORK is N.
!>          3. If the singular values and the left singular vectors are requested then:
!>             If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
!>             then the length of IWORK is N+M; otherwise the length of IWORK is N.
!>          4. If the singular values with both the left and the right singular vectors
!>             are requested, then:
!>             4.1. If LSAME(JOBV,'J') the length of IWORK is determined as follows:
!>                  If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
!>                  then the length of IWORK is N+M; otherwise the length of IWORK is N.
!>             4.2. If LSAME(JOBV,'V') the length of IWORK is determined as follows:
!>                  If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
!>                  then the length of IWORK is 2*N+M; otherwise the length of IWORK is 2*N.
!>
!>          On exit,
!>          IWORK(1) = the numerical rank determined after the initial
!>                     QR factorization with pivoting. See the descriptions
!>                     of JOBA and JOBR.
!>          IWORK(2) = the number of the computed nonzero singular values
!>          IWORK(3) = if nonzero, a warning message:
!>                     If IWORK(3) = 1 then some of the column norms of A
!>                     were denormalized floats. The requested high accuracy
!>                     is not warranted by the data.
!>          IWORK(4) = 1 or -1. If IWORK(4) = 1, then the procedure used A^* to
!>                     do the job as specified by the JOB parameters.
!>          If the call to ZGEJSV is a workspace query (indicated by LWORK = -1 or
!>          LRWORK = -1), then on exit IWORK(1) contains the required length of
!>          IWORK for the job parameters used in the call.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           < 0:  if INFO = -i, then the i-th argument had an illegal value.
!>           = 0:  successful exit;
!>           > 0:  ZGEJSV  did not converge in the maximal allowed number
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
!> \ingroup complex16GEsing
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  ZGEJSV implements a preconditioned Jacobi SVD algorithm. It uses ZGEQP3,
!>  ZGEQRF, and ZGELQF as preprocessors and preconditioners. Optionally, an
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
!>  & LAPACK routines called by ZGEJSV are implemented to work in that range.
!>  If that is not the case, then the restriction for safe computation with
!>  the singular values in the range of normalized IEEE numbers is that the
!>  spectral condition number kappa(A)=sigma_max(A)/sigma_min(A) does not
!>  overflow. This code (ZGEJSV) is best used in this restricted range,
!>  meaning that singular values of magnitude below ||A||_2 / DLAMCH('O') are
!>  returned as zeros. See JOBR for details on this.
!>     Further, this implementation is somewhat slower than the one described
!>  in [1,2] due to replacement of some non-LAPACK components, and because
!>  the choice of some tuning parameters in the iterative part (ZGESVJ) is
!>  left to the implementer on a particular machine.
!>     The rank revealing QR factorization (in this code: ZGEQP3) should be
!>  implemented as in [3]. We have a new version of ZGEQP3 under development
!>  that is more robust than the current one in LAPACK, with a cleaner cut in
!>  rank deficient cases. It will be available in the SIGMA library [4].
!>  If M is much larger than N, it is obvious that the initial QRF with
!>  column pivoting can be preprocessed by the QRF without pivoting. That
!>  well known trick is not used in ZGEJSV because in some cases heavy row
!>  weighting can be treated with complete pivoting. The overhead in cases
!>  M much larger than N is then only due to pivoting, but the benefits in
!>  terms of accuracy have prevailed. The implementer/user can incorporate
!>  this extra QRF step easily. The implementer can also improve data movement
!>  (matrix transpose, matrix copy, matrix transposed copy) - this
!>  implementation of ZGEJSV uses only the simplest, naive data movement.
!> \endverbatim
!
!> \par Contributor:
!  ==================
!>
!>  Zlatko Drmac, Department of Mathematics, Faculty of Science,
!>  University of Zagreb (Zagreb, Croatia); drmac@math.hr
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
!>     Department of Mathematics, University of Zagreb, 2008, 2016.
!> \endverbatim
!
!>  \par Bugs, examples and comments:
!   =================================
!>
!>  Please report all bugs and send interesting examples and/or comments to
!>  drmac@math.hr. Thank you.
!>
!  =====================================================================
      SUBROUTINE ZGEJSV(Joba,Jobu,Jobv,Jobr,Jobt,Jobp,M,N,A,Lda,Sva,U,  &
     &                  Ldu,V,Ldv,Cwork,Lwork,Rwork,Lrwork,Iwork,Info)
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      IMPLICIT NONE
!*--ZGEJSV579
      INTEGER Info , Lda , Ldu , Ldv , Lwork , Lrwork , M , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , U(Ldu,*) , V(Ldv,*) , Cwork(Lwork)
      DOUBLE PRECISION Sva(N) , Rwork(Lrwork)
      INTEGER Iwork(*)
      CHARACTER*1 Joba , Jobp , Jobr , Jobt , Jobu , Jobv
!     ..
!
!  ===========================================================================
!
!     .. Local Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0))
!     ..
!     .. Local Scalars ..
      COMPLEX*16 ctemp
      DOUBLE PRECISION aapp , aaqq , aatmax , aatmin , big , big1 ,     &
     &                 cond_ok , condr1 , condr2 , entra , entrat ,     &
     &                 epsln , maxprj , scalem , sconda , sfmin ,       &
     &                 small , temp1 , uscal1 , uscal2 , xsc
      INTEGER ierr , n1 , nr , numrank , p , q , warning
      LOGICAL almort , defr , errest , goscal , jracc , kill , lquery , &
     &        lsvec , l2aber , l2kill , l2pert , l2rank , l2tran ,      &
     &        noscal , rowpiv , rsvec , transp
!
      INTEGER optwrk , minwrk , minrwrk , miniwrk
      INTEGER lwcon , lwlqf , lwqp3 , lwqrf , lwunmlq , lwunmqr ,       &
     &        lwunmqrm , lwsvdj , lwsvdjv , lrwqp3 , lrwcon , lrwsvdj , &
     &        iwoff
      INTEGER lwrk_zgelqf , lwrk_zgeqp3 , lwrk_zgeqp3n , lwrk_zgeqrf ,  &
     &        lwrk_zgesvj , lwrk_zgesvjv , lwrk_zgesvju , lwrk_zunmlq , &
     &        lwrk_zunmqr , lwrk_zunmqrm
!     ..
!     .. Local Arrays
      COMPLEX*16 cdummy(1)
      DOUBLE PRECISION rdummy(1)
!
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DCMPLX , CONJG , DLOG , MAX , MIN , DBLE , NINT , &
     &          SQRT
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DZNRM2
      INTEGER IDAMAX , IZAMAX
      LOGICAL LSAME
      EXTERNAL IDAMAX , IZAMAX , LSAME , DLAMCH , DZNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL DLASSQ , ZCOPY , ZGELQF , ZGEQP3 , ZGEQRF , ZLACPY ,     &
     &         ZLAPMR , ZLASCL , DLASCL , ZLASET , ZLASSQ , ZLASWP ,    &
     &         ZUNGQR , ZUNMLQ , ZUNMQR , ZPOCON , DSCAL , ZDSCAL ,     &
     &         ZSWAP , ZTRSM , ZLACGV , XERBLA
!
      EXTERNAL ZGESVJ
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
      l2tran = LSAME(Jobt,'T') .AND. (M==N)
      l2kill = LSAME(Jobr,'R')
      defr = LSAME(Jobr,'N')
      l2pert = LSAME(Jobp,'P')
!
      lquery = (Lwork==-1) .OR. (Lrwork==-1)
!
      IF ( .NOT.(rowpiv .OR. l2rank .OR. l2aber .OR. errest .OR.        &
     &     LSAME(Joba,'C')) ) THEN
         Info = -1
      ELSEIF ( .NOT.(lsvec .OR. LSAME(Jobu,'N') .OR. (LSAME(Jobu,'W')   &
     &         .AND. rsvec .AND. l2tran)) ) THEN
         Info = -2
      ELSEIF ( .NOT.(rsvec .OR. LSAME(Jobv,'N') .OR. (LSAME(Jobv,'W')   &
     &         .AND. lsvec .AND. l2tran)) ) THEN
         Info = -3
      ELSEIF ( .NOT.(l2kill .OR. defr) ) THEN
         Info = -4
      ELSEIF ( .NOT.(LSAME(Jobt,'T') .OR. LSAME(Jobt,'N')) ) THEN
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
      ELSE
!        #:)
         Info = 0
      ENDIF
!
      IF ( Info==0 ) THEN
!         .. compute the minimal and the optimal workspace lengths
!         [[The expressions for computing the minimal and the optimal
!         values of LCWORK, LRWORK are written with a lot of redundancy and
!         can be simplified. However, this verbose form is useful for
!         maintenance and modifications of the code.]]
!
!        .. minimal workspace length for ZGEQP3 of an M x N matrix,
!         ZGEQRF of an N x N matrix, ZGELQF of an N x N matrix,
!         ZUNMLQ for computing N x N matrix, ZUNMQR for computing N x N
!         matrix, ZUNMQR for computing M x N matrix, respectively.
         lwqp3 = N + 1
         lwqrf = MAX(1,N)
         lwlqf = MAX(1,N)
         lwunmlq = MAX(1,N)
         lwunmqr = MAX(1,N)
         lwunmqrm = MAX(1,M)
!        .. minimal workspace length for ZPOCON of an N x N matrix
         lwcon = 2*N
!        .. minimal workspace length for ZGESVJ of an N x N matrix,
!         without and with explicit accumulation of Jacobi rotations
         lwsvdj = MAX(2*N,1)
         lwsvdjv = MAX(2*N,1)
!         .. minimal REAL workspace length for ZGEQP3, ZPOCON, ZGESVJ
         lrwqp3 = 2*N
         lrwcon = N
         lrwsvdj = N
         IF ( lquery ) THEN
            CALL ZGEQP3(M,N,A,Lda,Iwork,cdummy,cdummy,-1,rdummy,ierr)
            lwrk_zgeqp3 = cdummy(1)
            CALL ZGEQRF(N,N,A,Lda,cdummy,cdummy,-1,ierr)
            lwrk_zgeqrf = cdummy(1)
            CALL ZGELQF(N,N,A,Lda,cdummy,cdummy,-1,ierr)
            lwrk_zgelqf = cdummy(1)
         ENDIF
         minwrk = 2
         optwrk = 2
         miniwrk = N
         IF ( .NOT.(lsvec .OR. rsvec) ) THEN
!             .. minimal and optimal sizes of the complex workspace if
!             only the singular values are requested
            IF ( errest ) THEN
               minwrk = MAX(N+lwqp3,N**2+lwcon,N+lwqrf,lwsvdj)
            ELSE
               minwrk = MAX(N+lwqp3,N+lwqrf,lwsvdj)
            ENDIF
            IF ( lquery ) THEN
               CALL ZGESVJ('L','N','N',N,N,A,Lda,Sva,N,V,Ldv,cdummy,-1, &
     &                     rdummy,-1,ierr)
               lwrk_zgesvj = cdummy(1)
               IF ( errest ) THEN
                  optwrk = MAX(N+lwrk_zgeqp3,N**2+lwcon,N+lwrk_zgeqrf,  &
     &                     lwrk_zgesvj)
               ELSE
                  optwrk = MAX(N+lwrk_zgeqp3,N+lwrk_zgeqrf,lwrk_zgesvj)
               ENDIF
            ENDIF
            IF ( l2tran .OR. rowpiv ) THEN
               IF ( errest ) THEN
                  minrwrk = MAX(7,2*M,lrwqp3,lrwcon,lrwsvdj)
               ELSE
                  minrwrk = MAX(7,2*M,lrwqp3,lrwsvdj)
               ENDIF
            ELSEIF ( errest ) THEN
               minrwrk = MAX(7,lrwqp3,lrwcon,lrwsvdj)
            ELSE
               minrwrk = MAX(7,lrwqp3,lrwsvdj)
            ENDIF
            IF ( rowpiv .OR. l2tran ) miniwrk = miniwrk + M
         ELSEIF ( rsvec .AND. (.NOT.lsvec) ) THEN
!            .. minimal and optimal sizes of the complex workspace if the
!            singular values and the right singular vectors are requested
            IF ( errest ) THEN
               minwrk = MAX(N+lwqp3,lwcon,lwsvdj,N+lwlqf,2*N+lwqrf,     &
     &                  N+lwsvdj,N+lwunmlq)
            ELSE
               minwrk = MAX(N+lwqp3,lwsvdj,N+lwlqf,2*N+lwqrf,N+lwsvdj,  &
     &                  N+lwunmlq)
            ENDIF
            IF ( lquery ) THEN
               CALL ZGESVJ('L','U','N',N,N,U,Ldu,Sva,N,A,Lda,cdummy,-1, &
     &                     rdummy,-1,ierr)
               lwrk_zgesvj = cdummy(1)
               CALL ZUNMLQ('L','C',N,N,N,A,Lda,cdummy,V,Ldv,cdummy,-1,  &
     &                     ierr)
               lwrk_zunmlq = cdummy(1)
               IF ( errest ) THEN
                  optwrk = MAX(N+lwrk_zgeqp3,lwcon,lwrk_zgesvj,         &
     &                     N+lwrk_zgelqf,2*N+lwrk_zgeqrf,N+lwrk_zgesvj, &
     &                     N+lwrk_zunmlq)
               ELSE
                  optwrk = MAX(N+lwrk_zgeqp3,lwrk_zgesvj,N+lwrk_zgelqf, &
     &                     2*N+lwrk_zgeqrf,N+lwrk_zgesvj,N+lwrk_zunmlq)
               ENDIF
            ENDIF
            IF ( l2tran .OR. rowpiv ) THEN
               IF ( errest ) THEN
                  minrwrk = MAX(7,2*M,lrwqp3,lrwsvdj,lrwcon)
               ELSE
                  minrwrk = MAX(7,2*M,lrwqp3,lrwsvdj)
               ENDIF
            ELSEIF ( errest ) THEN
               minrwrk = MAX(7,lrwqp3,lrwsvdj,lrwcon)
            ELSE
               minrwrk = MAX(7,lrwqp3,lrwsvdj)
            ENDIF
            IF ( rowpiv .OR. l2tran ) miniwrk = miniwrk + M
         ELSEIF ( lsvec .AND. (.NOT.rsvec) ) THEN
!            .. minimal and optimal sizes of the complex workspace if the
!            singular values and the left singular vectors are requested
            IF ( errest ) THEN
               minwrk = N + MAX(lwqp3,lwcon,N+lwqrf,lwsvdj,lwunmqrm)
            ELSE
               minwrk = N + MAX(lwqp3,N+lwqrf,lwsvdj,lwunmqrm)
            ENDIF
            IF ( lquery ) THEN
               CALL ZGESVJ('L','U','N',N,N,U,Ldu,Sva,N,A,Lda,cdummy,-1, &
     &                     rdummy,-1,ierr)
               lwrk_zgesvj = cdummy(1)
               CALL ZUNMQR('L','N',M,N,N,A,Lda,cdummy,U,Ldu,cdummy,-1,  &
     &                     ierr)
               lwrk_zunmqrm = cdummy(1)
               IF ( errest ) THEN
                  optwrk = N + MAX(lwrk_zgeqp3,lwcon,N+lwrk_zgeqrf,     &
     &                     lwrk_zgesvj,lwrk_zunmqrm)
               ELSE
                  optwrk = N + MAX(lwrk_zgeqp3,N+lwrk_zgeqrf,           &
     &                     lwrk_zgesvj,lwrk_zunmqrm)
               ENDIF
            ENDIF
            IF ( l2tran .OR. rowpiv ) THEN
               IF ( errest ) THEN
                  minrwrk = MAX(7,2*M,lrwqp3,lrwsvdj,lrwcon)
               ELSE
                  minrwrk = MAX(7,2*M,lrwqp3,lrwsvdj)
               ENDIF
            ELSEIF ( errest ) THEN
               minrwrk = MAX(7,lrwqp3,lrwsvdj,lrwcon)
            ELSE
               minrwrk = MAX(7,lrwqp3,lrwsvdj)
            ENDIF
            IF ( rowpiv .OR. l2tran ) miniwrk = miniwrk + M
         ELSE
!            .. minimal and optimal sizes of the complex workspace if the
!            full SVD is requested
            IF ( .NOT.jracc ) THEN
               IF ( errest ) THEN
                  minwrk = MAX(N+lwqp3,N+lwcon,2*N+N**2+lwcon,2*N+lwqrf,&
     &                     2*N+lwqp3,2*N+N**2+N+lwlqf,2*N+N**2+N+N**2+  &
     &                     lwcon,2*N+N**2+N+lwsvdj,2*N+N**2+N+lwsvdjv,  &
     &                     2*N+N**2+N+lwunmqr,2*N+N**2+N+lwunmlq,       &
     &                     N+N**2+lwsvdj,N+lwunmqrm)
               ELSE
                  minwrk = MAX(N+lwqp3,2*N+N**2+lwcon,2*N+lwqrf,        &
     &                     2*N+lwqp3,2*N+N**2+N+lwlqf,2*N+N**2+N+N**2+  &
     &                     lwcon,2*N+N**2+N+lwsvdj,2*N+N**2+N+lwsvdjv,  &
     &                     2*N+N**2+N+lwunmqr,2*N+N**2+N+lwunmlq,       &
     &                     N+N**2+lwsvdj,N+lwunmqrm)
               ENDIF
               miniwrk = miniwrk + N
               IF ( rowpiv .OR. l2tran ) miniwrk = miniwrk + M
            ELSE
               IF ( errest ) THEN
                  minwrk = MAX(N+lwqp3,N+lwcon,2*N+lwqrf,               &
     &                     2*N+N**2+lwsvdjv,2*N+N**2+N+lwunmqr,         &
     &                     N+lwunmqrm)
               ELSE
                  minwrk = MAX(N+lwqp3,2*N+lwqrf,2*N+N**2+lwsvdjv,      &
     &                     2*N+N**2+N+lwunmqr,N+lwunmqrm)
               ENDIF
               IF ( rowpiv .OR. l2tran ) miniwrk = miniwrk + M
            ENDIF
            IF ( lquery ) THEN
               CALL ZUNMQR('L','N',M,N,N,A,Lda,cdummy,U,Ldu,cdummy,-1,  &
     &                     ierr)
               lwrk_zunmqrm = cdummy(1)
               CALL ZUNMQR('L','N',N,N,N,A,Lda,cdummy,U,Ldu,cdummy,-1,  &
     &                     ierr)
               lwrk_zunmqr = cdummy(1)
               IF ( .NOT.jracc ) THEN
                  CALL ZGEQP3(N,N,A,Lda,Iwork,cdummy,cdummy,-1,rdummy,  &
     &                        ierr)
                  lwrk_zgeqp3n = cdummy(1)
                  CALL ZGESVJ('L','U','N',N,N,U,Ldu,Sva,N,V,Ldv,cdummy, &
     &                        -1,rdummy,-1,ierr)
                  lwrk_zgesvj = cdummy(1)
                  CALL ZGESVJ('U','U','N',N,N,U,Ldu,Sva,N,V,Ldv,cdummy, &
     &                        -1,rdummy,-1,ierr)
                  lwrk_zgesvju = cdummy(1)
                  CALL ZGESVJ('L','U','V',N,N,U,Ldu,Sva,N,V,Ldv,cdummy, &
     &                        -1,rdummy,-1,ierr)
                  lwrk_zgesvjv = cdummy(1)
                  CALL ZUNMLQ('L','C',N,N,N,A,Lda,cdummy,V,Ldv,cdummy,  &
     &                        -1,ierr)
                  lwrk_zunmlq = cdummy(1)
                  IF ( errest ) THEN
                     optwrk = MAX(N+lwrk_zgeqp3,N+lwcon,2*N+N**2+lwcon, &
     &                        2*N+lwrk_zgeqrf,2*N+lwrk_zgeqp3n,         &
     &                        2*N+N**2+N+lwrk_zgelqf,2*N+N**2+N+N**2+   &
     &                        lwcon,2*N+N**2+N+lwrk_zgesvj,             &
     &                        2*N+N**2+N+lwrk_zgesvjv,                  &
     &                        2*N+N**2+N+lwrk_zunmqr,                   &
     &                        2*N+N**2+N+lwrk_zunmlq,                   &
     &                        N+N**2+lwrk_zgesvju,N+lwrk_zunmqrm)
                  ELSE
                     optwrk = MAX(N+lwrk_zgeqp3,2*N+N**2+lwcon,         &
     &                        2*N+lwrk_zgeqrf,2*N+lwrk_zgeqp3n,         &
     &                        2*N+N**2+N+lwrk_zgelqf,2*N+N**2+N+N**2+   &
     &                        lwcon,2*N+N**2+N+lwrk_zgesvj,             &
     &                        2*N+N**2+N+lwrk_zgesvjv,                  &
     &                        2*N+N**2+N+lwrk_zunmqr,                   &
     &                        2*N+N**2+N+lwrk_zunmlq,                   &
     &                        N+N**2+lwrk_zgesvju,N+lwrk_zunmqrm)
                  ENDIF
               ELSE
                  CALL ZGESVJ('L','U','V',N,N,U,Ldu,Sva,N,V,Ldv,cdummy, &
     &                        -1,rdummy,-1,ierr)
                  lwrk_zgesvjv = cdummy(1)
                  CALL ZUNMQR('L','N',N,N,N,cdummy,N,cdummy,V,Ldv,      &
     &                        cdummy,-1,ierr)
                  lwrk_zunmqr = cdummy(1)
                  CALL ZUNMQR('L','N',M,N,N,A,Lda,cdummy,U,Ldu,cdummy,  &
     &                        -1,ierr)
                  lwrk_zunmqrm = cdummy(1)
                  IF ( errest ) THEN
                     optwrk = MAX(N+lwrk_zgeqp3,N+lwcon,2*N+lwrk_zgeqrf,&
     &                        2*N+N**2,2*N+N**2+lwrk_zgesvjv,           &
     &                        2*N+N**2+N+lwrk_zunmqr,N+lwrk_zunmqrm)
                  ELSE
                     optwrk = MAX(N+lwrk_zgeqp3,2*N+lwrk_zgeqrf,        &
     &                        2*N+N**2,2*N+N**2+lwrk_zgesvjv,           &
     &                        2*N+N**2+N+lwrk_zunmqr,N+lwrk_zunmqrm)
                  ENDIF
               ENDIF
            ENDIF
            IF ( l2tran .OR. rowpiv ) THEN
               minrwrk = MAX(7,2*M,lrwqp3,lrwsvdj,lrwcon)
            ELSE
               minrwrk = MAX(7,lrwqp3,lrwsvdj,lrwcon)
            ENDIF
         ENDIF
         minwrk = MAX(2,minwrk)
         optwrk = MAX(minwrk,optwrk)
         IF ( Lwork<minwrk .AND. (.NOT.lquery) ) Info = -17
         IF ( Lrwork<minrwrk .AND. (.NOT.lquery) ) Info = -19
      ENDIF
!
      IF ( Info/=0 ) THEN
!       #:(
         CALL XERBLA('ZGEJSV',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         Cwork(1) = optwrk
         Cwork(2) = minwrk
         Rwork(1) = minrwrk
         Iwork(1) = MAX(4,miniwrk)
         RETURN
      ENDIF
!
!     Quick return for void matrix (Y3K safe)
! #:)
      IF ( (M==0) .OR. (N==0) ) THEN
         Iwork(1:4) = 0
         Rwork(1:7) = 0
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
      scalem = ONE/SQRT(DBLE(M)*DBLE(N))
      noscal = .TRUE.
      goscal = .TRUE.
      DO p = 1 , N
         aapp = ZERO
         aaqq = ONE
         CALL ZLASSQ(M,A(1,p),1,aapp,aaqq)
         IF ( aapp>big ) THEN
            Info = -9
            CALL XERBLA('ZGEJSV',-Info)
            RETURN
         ENDIF
         aaqq = SQRT(aaqq)
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
         IF ( lsvec ) CALL ZLASET('G',M,n1,CZERO,CONE,U,Ldu)
         IF ( rsvec ) CALL ZLASET('G',N,N,CZERO,CONE,V,Ldv)
         Rwork(1) = ONE
         Rwork(2) = ONE
         IF ( errest ) Rwork(3) = ONE
         IF ( lsvec .AND. rsvec ) THEN
            Rwork(4) = ONE
            Rwork(5) = ONE
         ENDIF
         IF ( l2tran ) THEN
            Rwork(6) = ZERO
            Rwork(7) = ZERO
         ENDIF
         Iwork(1) = 0
         Iwork(2) = 0
         Iwork(3) = 0
         Iwork(4) = -1
         RETURN
      ENDIF
!
!     Issue warning if denormalized column norms detected. Override the
!     high relative accuracy request. Issue licence to kill nonzero columns
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
            CALL ZLASCL('G',0,0,Sva(1),scalem,M,1,A(1,1),Lda,ierr)
            CALL ZLACPY('A',M,1,A,Lda,U,Ldu)
!           computing all M left singular vectors of the M x 1 matrix
            IF ( n1/=N ) THEN
               CALL ZGEQRF(M,N,U,Ldu,Cwork,Cwork(N+1),Lwork-N,ierr)
               CALL ZUNGQR(M,n1,1,U,Ldu,Cwork,Cwork(N+1),Lwork-N,ierr)
               CALL ZCOPY(M,A(1,1),1,U(1,1),1)
            ENDIF
         ENDIF
         IF ( rsvec ) V(1,1) = CONE
         IF ( Sva(1)<(big*scalem) ) THEN
            Sva(1) = Sva(1)/scalem
            scalem = ONE
         ENDIF
         Rwork(1) = ONE/scalem
         Rwork(2) = ONE
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
         Iwork(4) = -1
         IF ( errest ) Rwork(3) = ONE
         IF ( lsvec .AND. rsvec ) THEN
            Rwork(4) = ONE
            Rwork(5) = ONE
         ENDIF
         IF ( l2tran ) THEN
            Rwork(6) = ZERO
            Rwork(7) = ZERO
         ENDIF
         RETURN
!
      ENDIF
!
      transp = .FALSE.
!
      aatmax = -ONE
      aatmin = big
      IF ( rowpiv .OR. l2tran ) THEN
!
!     Compute the row norms, needed to determine row pivoting sequence
!     (in the case of heavily row weighted A, row pivoting is strongly
!     advised) and to collect information needed to compare the
!     structures of A * A^* and A^* * A (in the case L2TRAN.EQ..TRUE.).
!
         IF ( l2tran ) THEN
            DO p = 1 , M
               xsc = ZERO
               temp1 = ONE
               CALL ZLASSQ(N,A(p,1),Lda,xsc,temp1)
!              ZLASSQ gets both the ell_2 and the ell_infinity norm
!              in one pass through the vector
               Rwork(M+p) = xsc*scalem
               Rwork(p) = xsc*(scalem*SQRT(temp1))
               aatmax = MAX(aatmax,Rwork(p))
               IF ( Rwork(p)/=ZERO ) aatmin = MIN(aatmin,Rwork(p))
            ENDDO
         ELSE
            DO p = 1 , M
               Rwork(M+p) = scalem*ABS(A(p,IZAMAX(N,A(p,1),Lda)))
               aatmax = MAX(aatmax,Rwork(M+p))
               aatmin = MIN(aatmin,Rwork(M+p))
            ENDDO
         ENDIF
!
      ENDIF
!
!     For square matrix A try to determine whether A^*  would be better
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
!        Now, SVA().^2/Trace(A^* * A) is a point in the probability simplex.
!        It is derived from the diagonal of  A^* * A.  Do the same with the
!        diagonal of A * A^*, compute the entropy of the corresponding
!        probability distribution. Note that A * A^* and A^* * A have the
!        same trace.
!
         entrat = ZERO
         DO p = 1 , M
            big1 = ((Rwork(p)/xsc)**2)*temp1
            IF ( big1/=ZERO ) entrat = entrat + big1*DLOG(big1)
         ENDDO
         entrat = -entrat/DLOG(DBLE(M))
!
!        Analyze the entropies and decide A or A^*. Smaller entropy
!        usually means better input for the algorithm.
!
         transp = (entrat<entra)
!
!        If A^* is better than A, take the adjoint of A. This is allowed
!        only for square matrices, M=N.
         IF ( transp ) THEN
!           In an optimal implementation, this trivial transpose
!           should be replaced with faster transpose.
            DO p = 1 , N - 1
               A(p,p) = CONJG(A(p,p))
               DO q = p + 1 , N
                  ctemp = CONJG(A(q,p))
                  A(q,p) = CONJG(A(p,q))
                  A(p,q) = ctemp
               ENDDO
            ENDDO
            A(N,N) = CONJG(A(N,N))
            DO p = 1 , N
               Rwork(M+p) = Sva(p)
               Sva(p) = Rwork(p)
!              previously computed row 2-norms are now column 2-norms
!              of the transposed matrix
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
!     than SQRT(BIG) -- the matrix is scaled so that its maximal column
!     has Euclidean norm equal to SQRT(BIG/N). The only reason to keep
!     SQRT(BIG) instead of BIG is the fact that ZGEJSV uses LAPACK and
!     BLAS routines that, in some implementations, are not capable of
!     working in the full interval [SFMIN,BIG] and that they may provoke
!     overflows in the intermediate results. If the singular values spread
!     from SFMIN to BIG, then ZGESVJ will compute them. So, in that case,
!     one should use ZGESVJ instead of ZGEJSV.
!     >> change in the April 2016 update: allow bigger range, i.e. the
!     largest column is allowed up to BIG/N and ZGESVJ will do the rest.
      big1 = SQRT(big)
      temp1 = SQRT(big/DBLE(N))
!      TEMP1  = BIG/DBLE(N)
!
      CALL DLASCL('G',0,0,aapp,temp1,N,1,Sva,N,ierr)
      IF ( aaqq>(aapp*sfmin) ) THEN
         aaqq = (aaqq/aapp)*temp1
      ELSE
         aaqq = (aaqq*temp1)/aapp
      ENDIF
      temp1 = temp1*scalem
      CALL ZLASCL('G',0,0,aapp,temp1,M,N,A,Lda,ierr)
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
!        sigma_max(A) / sigma_min(A) approx. SQRT(BIG)/SQRT(SFMIN).
         xsc = SQRT(sfmin)
      ELSE
         xsc = small
!
!        Now, if the condition number of A is too big,
!        sigma_max(A) / sigma_min(A) .GT. SQRT(BIG/N) * EPSLN / SFMIN,
!        as a precaution measure, the full SVD is computed using ZGESVJ
!        with accumulated Jacobi rotations. This provides numerically
!        more robust computation, at the cost of slightly increased run
!        time. Depending on the concrete implementation of BLAS and LAPACK
!        (i.e. how they behave in presence of extreme ill-conditioning) the
!        implementor may decide to remove this switch.
         IF ( (aaqq<SQRT(sfmin)) .AND. lsvec .AND. rsvec )              &
     &        jracc = .TRUE.
!
      ENDIF
      IF ( aaqq<xsc ) THEN
         DO p = 1 , N
            IF ( Sva(p)<xsc ) THEN
               CALL ZLASET('A',M,1,CZERO,CZERO,A(1,p),Lda)
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
         IF ( (lsvec .AND. rsvec) .AND. .NOT.(jracc) ) THEN
            iwoff = 2*N
         ELSE
            iwoff = N
         ENDIF
         DO p = 1 , M - 1
            q = IDAMAX(M-p+1,Rwork(M+p),1) + p - 1
            Iwork(iwoff+p) = q
            IF ( p/=q ) THEN
               temp1 = Rwork(M+p)
               Rwork(M+p) = Rwork(M+q)
               Rwork(M+q) = temp1
            ENDIF
         ENDDO
         CALL ZLASWP(N,A,Lda,1,M-1,Iwork(iwoff+1),1)
      ENDIF
!
!     End of the preparation phase (scaling, optional sorting and
!     transposing, optional flushing of small columns).
!
!     Preconditioning
!
!     If the full SVD is needed, the right singular vectors are computed
!     from a matrix equation, and for that we need theoretical analysis
!     of the Businger-Golub pivoting. So we use ZGEQP3 as the first RR QRF.
!     In all other cases the first RR QRF can be chosen by other criteria
!     (eg speed by replacing global with restricted window pivoting, such
!     as in xGEQPX from TOMS # 782). Good results will be obtained using
!     xGEQPX with properly (!) chosen numerical parameters.
!     Any improvement of ZGEQP3 improves overall performance of ZGEJSV.
!
!     A * P1 = Q1 * [ R1^* 0]^*:
      DO p = 1 , N
!        .. all columns are free columns
         Iwork(p) = 0
      ENDDO
      CALL ZGEQP3(M,N,A,Lda,Iwork,Cwork,Cwork(N+1),Lwork-N,Rwork,ierr)
!
!     The upper triangular matrix R1 from the first QRF is inspected for
!     rank deficiency and possibilities for deflation, or possible
!     ill-conditioning. Depending on the user specified flag L2RANK,
!     the procedure explores possibilities to reduce the numerical
!     rank by inspecting the computed upper triangular factor. If
!     L2RANK or L2ABER are up, then ZGEJSV will compute the SVD of
!     A + dA, where ||dA|| <= f(M,N)*EPSLN.
!
      nr = 1
      IF ( l2aber ) THEN
!        Standard absolute error bound suffices. All sigma_i with
!        sigma_i < N*EPSLN*||A|| are flushed to zero. This is an
!        aggressive enforcement of lower numerical rank by introducing a
!        backward error of the order of N*EPSLN*||A||.
         temp1 = SQRT(DBLE(N))*epsln
         DO p = 2 , N
            IF ( ABS(A(p,p))<(temp1*ABS(A(1,1))) ) EXIT
            nr = nr + 1
         ENDDO
      ELSEIF ( l2rank ) THEN
!        .. similarly as above, only slightly more gentle (less aggressive).
!        Sudden drop on the diagonal of R1 is used as the criterion for
!        close-to-rank-deficient.
         temp1 = SQRT(sfmin)
         DO p = 2 , N
            IF ( (ABS(A(p,p))<(epsln*ABS(A(p-1,p-1)))) .OR.             &
     &           (ABS(A(p,p))<small) .OR.                               &
     &           (l2kill .AND. (ABS(A(p,p))<temp1)) ) EXIT
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
         temp1 = SQRT(sfmin)
         DO p = 2 , N
            IF ( (ABS(A(p,p))<small) .OR.                               &
     &           (l2kill .AND. (ABS(A(p,p))<temp1)) ) EXIT
            nr = nr + 1
         ENDDO
!
      ENDIF
!
      almort = .FALSE.
      IF ( nr==N ) THEN
         maxprj = ONE
         DO p = 2 , N
            temp1 = ABS(A(p,p))/Sva(Iwork(p))
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
               CALL ZLACPY('U',N,N,A,Lda,V,Ldv)
               DO p = 1 , N
                  temp1 = Sva(Iwork(p))
                  CALL ZDSCAL(p,ONE/temp1,V(1,p),1)
               ENDDO
               IF ( lsvec ) THEN
                  CALL ZPOCON('U',N,V,Ldv,ONE,temp1,Cwork(N+1),Rwork,   &
     &                        ierr)
               ELSE
                  CALL ZPOCON('U',N,V,Ldv,ONE,temp1,Cwork,Rwork,ierr)
               ENDIF
!
            ELSEIF ( lsvec ) THEN
!              .. U is available as workspace
               CALL ZLACPY('U',N,N,A,Lda,U,Ldu)
               DO p = 1 , N
                  temp1 = Sva(Iwork(p))
                  CALL ZDSCAL(p,ONE/temp1,U(1,p),1)
               ENDDO
               CALL ZPOCON('U',N,U,Ldu,ONE,temp1,Cwork(N+1),Rwork,ierr)
            ELSE
               CALL ZLACPY('U',N,N,A,Lda,Cwork,N)
![]            CALL ZLACPY( 'U', N, N, A, LDA, CWORK(N+1), N )
!              Change: here index shifted by N to the left, CWORK(1:N)
!              not needed for SIGMA only computation
               DO p = 1 , N
                  temp1 = Sva(Iwork(p))
![]               CALL ZDSCAL( p, ONE/TEMP1, CWORK(N+(p-1)*N+1), 1 )
                  CALL ZDSCAL(p,ONE/temp1,Cwork((p-1)*N+1),1)
               ENDDO
!           .. the columns of R are scaled to have unit Euclidean lengths.
![]               CALL ZPOCON( 'U', N, CWORK(N+1), N, ONE, TEMP1,
![]     $              CWORK(N+N*N+1), RWORK, IERR )
               CALL ZPOCON('U',N,Cwork,N,ONE,temp1,Cwork(N*N+1),Rwork,  &
     &                     ierr)
!
            ENDIF
            IF ( temp1/=ZERO ) THEN
               sconda = ONE/SQRT(temp1)
            ELSE
               sconda = -ONE
            ENDIF
!           SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1).
!           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
         ELSE
            sconda = -ONE
         ENDIF
      ENDIF
!
      l2pert = l2pert .AND. (ABS(A(1,1)/A(nr,nr))>SQRT(big1))
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
            CALL ZCOPY(N-p,A(p,p+1),Lda,A(p+1,p),1)
            CALL ZLACGV(N-p+1,A(p,p),1)
         ENDDO
         IF ( nr==N ) A(N,N) = CONJG(A(N,N))
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
!              XSC = SQRT(SMALL)
               xsc = epsln/DBLE(N)
               DO q = 1 , nr
                  ctemp = DCMPLX(xsc*ABS(A(q,q)),ZERO)
                  DO p = 1 , N
!     $                     A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) )
                     IF ( ((p>q) .AND. (ABS(A(p,q))<=temp1)) .OR.       &
     &                    (p<q) ) A(p,q) = ctemp
                  ENDDO
               ENDDO
            ELSE
               CALL ZLASET('U',nr-1,nr-1,CZERO,CZERO,A(1,2),Lda)
            ENDIF
!
!            .. second preconditioning using the QR factorization
!
            CALL ZGEQRF(N,nr,A,Lda,Cwork,Cwork(N+1),Lwork-N,ierr)
!
!           .. and transpose upper to lower triangular
            DO p = 1 , nr - 1
               CALL ZCOPY(nr-p,A(p,p+1),Lda,A(p+1,p),1)
               CALL ZLACGV(nr-p+1,A(p,p),1)
            ENDDO
!
         ENDIF
!
!           Row-cyclic Jacobi SVD algorithm with column pivoting
!
!           .. again some perturbation (a "background noise") is added
!           to drown denormals
         IF ( l2pert ) THEN
!              XSC = SQRT(SMALL)
            xsc = epsln/DBLE(N)
            DO q = 1 , nr
               ctemp = DCMPLX(xsc*ABS(A(q,q)),ZERO)
               DO p = 1 , nr
!     $                   A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) )
                  IF ( ((p>q) .AND. (ABS(A(p,q))<=temp1)) .OR. (p<q) )  &
     &                 A(p,q) = ctemp
               ENDDO
            ENDDO
         ELSE
            CALL ZLASET('U',nr-1,nr-1,CZERO,CZERO,A(1,2),Lda)
         ENDIF
!
!           .. and one-sided Jacobi rotations are started on a lower
!           triangular matrix (plus perturbation which is ignored in
!           the part which destroys triangular form (confusing?!))
!
         CALL ZGESVJ('L','N','N',nr,nr,A,Lda,Sva,N,V,Ldv,Cwork,Lwork,   &
     &               Rwork,Lrwork,Info)
!
         scalem = Rwork(1)
         numrank = NINT(Rwork(2))
!
!
      ELSEIF ( (rsvec .AND. (.NOT.lsvec) .AND. (.NOT.jracc)) .OR.       &
     &         (jracc .AND. (.NOT.lsvec) .AND. (nr/=N)) ) THEN
!
!        -> Singular Values and Right Singular Vectors <-
!
         IF ( almort ) THEN
!
!           .. in this case NR equals N
            DO p = 1 , nr
               CALL ZCOPY(N-p+1,A(p,p),Lda,V(p,p),1)
               CALL ZLACGV(N-p+1,V(p,p),1)
            ENDDO
            CALL ZLASET('U',nr-1,nr-1,CZERO,CZERO,V(1,2),Ldv)
!
            CALL ZGESVJ('L','U','N',N,nr,V,Ldv,Sva,nr,A,Lda,Cwork,Lwork,&
     &                  Rwork,Lrwork,Info)
            scalem = Rwork(1)
            numrank = NINT(Rwork(2))
 
         ELSE
!
!        .. two more QR factorizations ( one QRF is not enough, two require
!        accumulated product of Jacobi rotations, three are perfect )
!
            CALL ZLASET('L',nr-1,nr-1,CZERO,CZERO,A(2,1),Lda)
            CALL ZGELQF(nr,N,A,Lda,Cwork,Cwork(N+1),Lwork-N,ierr)
            CALL ZLACPY('L',nr,nr,A,Lda,V,Ldv)
            CALL ZLASET('U',nr-1,nr-1,CZERO,CZERO,V(1,2),Ldv)
            CALL ZGEQRF(nr,nr,V,Ldv,Cwork(N+1),Cwork(2*N+1),Lwork-2*N,  &
     &                  ierr)
            DO p = 1 , nr
               CALL ZCOPY(nr-p+1,V(p,p),Ldv,V(p,p),1)
               CALL ZLACGV(nr-p+1,V(p,p),1)
            ENDDO
            CALL ZLASET('U',nr-1,nr-1,CZERO,CZERO,V(1,2),Ldv)
!
            CALL ZGESVJ('L','U','N',nr,nr,V,Ldv,Sva,nr,U,Ldu,Cwork(N+1),&
     &                  Lwork-N,Rwork,Lrwork,Info)
            scalem = Rwork(1)
            numrank = NINT(Rwork(2))
            IF ( nr<N ) THEN
               CALL ZLASET('A',N-nr,nr,CZERO,CZERO,V(nr+1,1),Ldv)
               CALL ZLASET('A',nr,N-nr,CZERO,CZERO,V(1,nr+1),Ldv)
               CALL ZLASET('A',N-nr,N-nr,CZERO,CONE,V(nr+1,nr+1),Ldv)
            ENDIF
!
            CALL ZUNMLQ('L','C',N,N,nr,A,Lda,Cwork,V,Ldv,Cwork(N+1),    &
     &                  Lwork-N,ierr)
!
         ENDIF
!         .. permute the rows of V
!         DO 8991 p = 1, N
!            CALL ZCOPY( N, V(p,1), LDV, A(IWORK(p),1), LDA )
! 8991    CONTINUE
!         CALL ZLACPY( 'All', N, N, A, LDA, V, LDV )
         CALL ZLAPMR(.FALSE.,N,N,V,Ldv,Iwork)
!
         IF ( transp ) CALL ZLACPY('A',N,N,V,Ldv,U,Ldu)
!
      ELSEIF ( jracc .AND. (.NOT.lsvec) .AND. (nr==N) ) THEN
!
         CALL ZLASET('L',N-1,N-1,CZERO,CZERO,A(2,1),Lda)
!
         CALL ZGESVJ('U','N','V',N,N,A,Lda,Sva,N,V,Ldv,Cwork,Lwork,     &
     &               Rwork,Lrwork,Info)
         scalem = Rwork(1)
         numrank = NINT(Rwork(2))
         CALL ZLAPMR(.FALSE.,N,N,V,Ldv,Iwork)
!
      ELSEIF ( lsvec .AND. (.NOT.rsvec) ) THEN
!
!        .. Singular Values and Left Singular Vectors                 ..
!
!        .. second preconditioning step to avoid need to accumulate
!        Jacobi rotations in the Jacobi iterations.
         DO p = 1 , nr
            CALL ZCOPY(N-p+1,A(p,p),Lda,U(p,p),1)
            CALL ZLACGV(N-p+1,U(p,p),1)
         ENDDO
         CALL ZLASET('U',nr-1,nr-1,CZERO,CZERO,U(1,2),Ldu)
!
         CALL ZGEQRF(N,nr,U,Ldu,Cwork(N+1),Cwork(2*N+1),Lwork-2*N,ierr)
!
         DO p = 1 , nr - 1
            CALL ZCOPY(nr-p,U(p,p+1),Ldu,U(p+1,p),1)
            CALL ZLACGV(N-p+1,U(p,p),1)
         ENDDO
         CALL ZLASET('U',nr-1,nr-1,CZERO,CZERO,U(1,2),Ldu)
!
         CALL ZGESVJ('L','U','N',nr,nr,U,Ldu,Sva,nr,A,Lda,Cwork(N+1),   &
     &               Lwork-N,Rwork,Lrwork,Info)
         scalem = Rwork(1)
         numrank = NINT(Rwork(2))
!
         IF ( nr<M ) THEN
            CALL ZLASET('A',M-nr,nr,CZERO,CZERO,U(nr+1,1),Ldu)
            IF ( nr<n1 ) THEN
               CALL ZLASET('A',nr,n1-nr,CZERO,CZERO,U(1,nr+1),Ldu)
               CALL ZLASET('A',M-nr,n1-nr,CZERO,CONE,U(nr+1,nr+1),Ldu)
            ENDIF
         ENDIF
!
         CALL ZUNMQR('L','N',M,n1,N,A,Lda,Cwork,U,Ldu,Cwork(N+1),       &
     &               Lwork-N,ierr)
!
         IF ( rowpiv ) CALL ZLASWP(n1,U,Ldu,1,M-1,Iwork(iwoff+1),-1)
!
         DO p = 1 , n1
            xsc = ONE/DZNRM2(M,U(1,p),1)
            CALL ZDSCAL(M,xsc,U(1,p),1)
         ENDDO
!
         IF ( transp ) CALL ZLACPY('A',N,N,U,Ldu,V,Ldv)
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
!        in presence of extreme values, e.g. when the singular values spread from
!        the underflow to the overflow threshold.
!
            DO p = 1 , nr
               CALL ZCOPY(N-p+1,A(p,p),Lda,V(p,p),1)
               CALL ZLACGV(N-p+1,V(p,p),1)
            ENDDO
!
            IF ( l2pert ) THEN
               xsc = SQRT(small/epsln)
               DO q = 1 , nr
                  ctemp = DCMPLX(xsc*ABS(V(q,q)),ZERO)
                  DO p = 1 , N
!     $                V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) )
                     IF ( (p>q) .AND. (ABS(V(p,q))<=temp1) .OR. (p<q) ) &
     &                    V(p,q) = ctemp
                     IF ( p<q ) V(p,q) = -V(p,q)
                  ENDDO
               ENDDO
            ELSE
               CALL ZLASET('U',nr-1,nr-1,CZERO,CZERO,V(1,2),Ldv)
            ENDIF
 
            CALL ZGEQRF(N,nr,V,Ldv,Cwork(N+1),Cwork(2*N+1),Lwork-2*N,   &
     &                  ierr)
            CALL ZLACPY('L',N,nr,V,Ldv,Cwork(2*N+1),N)
!
            DO p = 1 , nr
               CALL ZCOPY(nr-p+1,V(p,p),Ldv,U(p,p),1)
               CALL ZLACGV(nr-p+1,U(p,p),1)
            ENDDO
 
            IF ( l2pert ) THEN
               xsc = SQRT(small/epsln)
               DO q = 2 , nr
                  DO p = 1 , q - 1
                     ctemp = DCMPLX(xsc*MIN(ABS(U(p,p)),ABS(U(q,q))),   &
     &                       ZERO)
!                  U(p,q) = - TEMP1 * ( U(q,p) / ABS(U(q,p)) )
                     U(p,q) = -ctemp
                  ENDDO
               ENDDO
            ELSE
               CALL ZLASET('U',nr-1,nr-1,CZERO,CZERO,U(1,2),Ldu)
            ENDIF
 
            CALL ZGESVJ('L','U','V',nr,nr,U,Ldu,Sva,N,V,Ldv,            &
     &                  Cwork(2*N+N*nr+1),Lwork-2*N-N*nr,Rwork,Lrwork,  &
     &                  Info)
            scalem = Rwork(1)
            numrank = NINT(Rwork(2))
 
            IF ( nr<N ) THEN
               CALL ZLASET('A',N-nr,nr,CZERO,CZERO,V(nr+1,1),Ldv)
               CALL ZLASET('A',nr,N-nr,CZERO,CZERO,V(1,nr+1),Ldv)
               CALL ZLASET('A',N-nr,N-nr,CZERO,CONE,V(nr+1,nr+1),Ldv)
            ENDIF
 
            CALL ZUNMQR('L','N',N,N,nr,Cwork(2*N+1),N,Cwork(N+1),V,Ldv, &
     &                  Cwork(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,ierr)
!
!           Permute the rows of V using the (column) permutation from the
!           first QRF. Also, scale the columns to make them unit in
!           Euclidean norm. This applies to all cases.
!
            temp1 = SQRT(DBLE(N))*epsln
            DO q = 1 , N
               DO p = 1 , N
                  Cwork(2*N+N*nr+nr+Iwork(p)) = V(p,q)
               ENDDO
               DO p = 1 , N
                  V(p,q) = Cwork(2*N+N*nr+nr+p)
               ENDDO
               xsc = ONE/DZNRM2(N,V(1,q),1)
               IF ( (xsc<(ONE-temp1)) .OR. (xsc>(ONE+temp1)) )          &
     &              CALL ZDSCAL(N,xsc,V(1,q),1)
            ENDDO
!
!           At this moment, V contains the right singular vectors of A.
!           Next, assemble the left singular vector matrix U (M x N).
!
            IF ( nr<M ) THEN
               CALL ZLASET('A',M-nr,nr,CZERO,CZERO,U(nr+1,1),Ldu)
               IF ( nr<n1 ) THEN
                  CALL ZLASET('A',nr,n1-nr,CZERO,CZERO,U(1,nr+1),Ldu)
                  CALL ZLASET('A',M-nr,n1-nr,CZERO,CONE,U(nr+1,nr+1),   &
     &                        Ldu)
               ENDIF
            ENDIF
!
            CALL ZUNMQR('L','N',M,n1,N,A,Lda,Cwork,U,Ldu,Cwork(N+1),    &
     &                  Lwork-N,ierr)
!
            IF ( rowpiv ) CALL ZLASWP(n1,U,Ldu,1,M-1,Iwork(iwoff+1),-1)
!
         ELSEIF ( .NOT.almort ) THEN
!
!           Second Preconditioning Step (QRF [with pivoting])
!           Note that the composition of TRANSPOSE, QRF and TRANSPOSE is
!           equivalent to an LQF CALL. Since in many libraries the QRF
!           seems to be better optimized than the LQF, we do explicit
!           transpose and use the QRF. This is subject to changes in an
!           optimized implementation of ZGEJSV.
!
            DO p = 1 , nr
               CALL ZCOPY(N-p+1,A(p,p),Lda,V(p,p),1)
               CALL ZLACGV(N-p+1,V(p,p),1)
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
               xsc = SQRT(small)
               DO q = 1 , nr
                  ctemp = DCMPLX(xsc*ABS(V(q,q)),ZERO)
                  DO p = 1 , N
!     $                   V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) )
                     IF ( (p>q) .AND. (ABS(V(p,q))<=temp1) .OR. (p<q) ) &
     &                    V(p,q) = ctemp
                     IF ( p<q ) V(p,q) = -V(p,q)
                  ENDDO
               ENDDO
            ELSE
               CALL ZLASET('U',nr-1,nr-1,CZERO,CZERO,V(1,2),Ldv)
            ENDIF
!
!           Estimate the row scaled condition number of R1
!           (If R1 is rectangular, N > NR, then the condition number
!           of the leading NR x NR submatrix is estimated.)
!
            CALL ZLACPY('L',nr,nr,V,Ldv,Cwork(2*N+1),nr)
            DO p = 1 , nr
               temp1 = DZNRM2(nr-p+1,Cwork(2*N+(p-1)*nr+p),1)
               CALL ZDSCAL(nr-p+1,ONE/temp1,Cwork(2*N+(p-1)*nr+p),1)
            ENDDO
            CALL ZPOCON('L',nr,Cwork(2*N+1),nr,ONE,temp1,               &
     &                  Cwork(2*N+nr*nr+1),Rwork,ierr)
            condr1 = ONE/SQRT(temp1)
!           .. here need a second opinion on the condition number
!           .. then assume worst case scenario
!           R1 is OK for inverse <=> CONDR1 .LT. DBLE(N)
!           more conservative    <=> CONDR1 .LT. SQRT(DBLE(N))
!
            cond_ok = SQRT(SQRT(DBLE(nr)))
![TP]       COND_OK is a tuning parameter.
!
            IF ( condr1<cond_ok ) THEN
!              .. the second QRF without pivoting. Note: in an optimized
!              implementation, this QRF should be implemented as the QRF
!              of a lower triangular matrix.
!              R1^* = Q2 * R2
               CALL ZGEQRF(N,nr,V,Ldv,Cwork(N+1),Cwork(2*N+1),Lwork-2*N,&
     &                     ierr)
!
               IF ( l2pert ) THEN
                  xsc = SQRT(small)/epsln
                  DO p = 2 , nr
                     DO q = 1 , p - 1
                        ctemp = DCMPLX(xsc*MIN(ABS(V(p,p)),ABS(V(q,q))),&
     &                          ZERO)
!     $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) )
                        IF ( ABS(V(q,p))<=temp1 ) V(q,p) = ctemp
                     ENDDO
                  ENDDO
               ENDIF
!
               IF ( nr/=N ) CALL ZLACPY('A',N,nr,V,Ldv,Cwork(2*N+1),N)
!              .. save ...
!
!           .. this transposed copy should be better than naive
               DO p = 1 , nr - 1
                  CALL ZCOPY(nr-p,V(p,p+1),Ldv,V(p+1,p),1)
                  CALL ZLACGV(nr-p+1,V(p,p),1)
               ENDDO
               V(nr,nr) = CONJG(V(nr,nr))
!
               condr2 = condr1
!
            ELSE
!
!              .. ill-conditioned case: second QRF with pivoting
!              Note that windowed pivoting would be equally good
!              numerically, and more run-time efficient. So, in
!              an optimal implementation, the next call to ZGEQP3
!              should be replaced with eg. CALL ZGEQPX (ACM TOMS #782)
!              with properly (carefully) chosen parameters.
!
!              R1^* * P2 = Q2 * R2
               DO p = 1 , nr
                  Iwork(N+p) = 0
               ENDDO
               CALL ZGEQP3(N,nr,V,Ldv,Iwork(N+1),Cwork(N+1),Cwork(2*N+1)&
     &                     ,Lwork-2*N,Rwork,ierr)
!*               CALL ZGEQRF( N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1),
!*     $              LWORK-2*N, IERR )
               IF ( l2pert ) THEN
                  xsc = SQRT(small)
                  DO p = 2 , nr
                     DO q = 1 , p - 1
                        ctemp = DCMPLX(xsc*MIN(ABS(V(p,p)),ABS(V(q,q))),&
     &                          ZERO)
!     $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) )
                        IF ( ABS(V(q,p))<=temp1 ) V(q,p) = ctemp
                     ENDDO
                  ENDDO
               ENDIF
!
               CALL ZLACPY('A',N,nr,V,Ldv,Cwork(2*N+1),N)
!
               IF ( l2pert ) THEN
                  xsc = SQRT(small)
                  DO p = 2 , nr
                     DO q = 1 , p - 1
                        ctemp = DCMPLX(xsc*MIN(ABS(V(p,p)),ABS(V(q,q))),&
     &                          ZERO)
!                        V(p,q) = - TEMP1*( V(q,p) / ABS(V(q,p)) )
                        V(p,q) = -ctemp
                     ENDDO
                  ENDDO
               ELSE
                  CALL ZLASET('L',nr-1,nr-1,CZERO,CZERO,V(2,1),Ldv)
               ENDIF
!              Now, compute R2 = L3 * Q3, the LQ factorization.
               CALL ZGELQF(nr,nr,V,Ldv,Cwork(2*N+N*nr+1),               &
     &                     Cwork(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,ierr)
!              .. and estimate the condition number
               CALL ZLACPY('L',nr,nr,V,Ldv,Cwork(2*N+N*nr+nr+1),nr)
               DO p = 1 , nr
                  temp1 = DZNRM2(p,Cwork(2*N+N*nr+nr+p),nr)
                  CALL ZDSCAL(p,ONE/temp1,Cwork(2*N+N*nr+nr+p),nr)
               ENDDO
               CALL ZPOCON('L',nr,Cwork(2*N+N*nr+nr+1),nr,ONE,temp1,    &
     &                     Cwork(2*N+N*nr+nr+nr*nr+1),Rwork,ierr)
               condr2 = ONE/SQRT(temp1)
!
!
!                 .. save the Householder vectors used for Q3
!                 (this overwrites the copy of R2, as it will not be
!                 needed in this branch, but it does not overwritte the
!                 Huseholder vectors of Q2.).
!                 .. and the rest of the information on Q3 is in
!                 WORK(2*N+N*NR+1:2*N+N*NR+N)
               IF ( condr2>=cond_ok )                                   &
     &              CALL ZLACPY('U',nr,nr,V,Ldv,Cwork(2*N+1),N)
!
            ENDIF
!
            IF ( l2pert ) THEN
               xsc = SQRT(small)
               DO q = 2 , nr
                  ctemp = xsc*V(q,q)
                  DO p = 1 , q - 1
!                     V(p,q) = - TEMP1*( V(p,q) / ABS(V(p,q)) )
                     V(p,q) = -ctemp
                  ENDDO
               ENDDO
            ELSE
               CALL ZLASET('U',nr-1,nr-1,CZERO,CZERO,V(1,2),Ldv)
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
               CALL ZGESVJ('L','U','N',nr,nr,V,Ldv,Sva,nr,U,Ldu,        &
     &                     Cwork(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,Rwork,&
     &                     Lrwork,Info)
               scalem = Rwork(1)
               numrank = NINT(Rwork(2))
               DO p = 1 , nr
                  CALL ZCOPY(nr,V(1,p),1,U(1,p),1)
                  CALL ZDSCAL(nr,Sva(p),V(1,p),1)
               ENDDO
 
!        .. pick the right matrix equation and solve it
!
               IF ( nr==N ) THEN
! :))             .. best case, R1 is inverted. The solution of this matrix
!                 equation is Q2*V2 = the product of the Jacobi rotations
!                 used in ZGESVJ, premultiplied with the orthogonal matrix
!                 from the second QR factorization.
                  CALL ZTRSM('L','U','N','N',nr,nr,CONE,A,Lda,V,Ldv)
               ELSE
!                 .. R1 is well conditioned, but non-square. Adjoint of R2
!                 is inverted to get the product of the Jacobi rotations
!                 used in ZGESVJ. The Q-factor from the second QR
!                 factorization is then built in explicitly.
                  CALL ZTRSM('L','U','C','N',nr,nr,CONE,Cwork(2*N+1),N, &
     &                       V,Ldv)
                  IF ( nr<N ) THEN
                     CALL ZLASET('A',N-nr,nr,CZERO,CZERO,V(nr+1,1),Ldv)
                     CALL ZLASET('A',nr,N-nr,CZERO,CZERO,V(1,nr+1),Ldv)
                     CALL ZLASET('A',N-nr,N-nr,CZERO,CONE,V(nr+1,nr+1), &
     &                           Ldv)
                  ENDIF
                  CALL ZUNMQR('L','N',N,N,nr,Cwork(2*N+1),N,Cwork(N+1), &
     &                        V,Ldv,Cwork(2*N+N*nr+nr+1),               &
     &                        Lwork-2*N-N*nr-nr,ierr)
               ENDIF
!
            ELSEIF ( condr2<cond_ok ) THEN
!
!              The matrix R2 is inverted. The solution of the matrix equation
!              is Q3^* * V3 = the product of the Jacobi rotations (appplied to
!              the lower triangular L3 from the LQ factorization of
!              R2=L3*Q3), pre-multiplied with the transposed Q3.
               CALL ZGESVJ('L','U','N',nr,nr,V,Ldv,Sva,nr,U,Ldu,        &
     &                     Cwork(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,Rwork,&
     &                     Lrwork,Info)
               scalem = Rwork(1)
               numrank = NINT(Rwork(2))
               DO p = 1 , nr
                  CALL ZCOPY(nr,V(1,p),1,U(1,p),1)
                  CALL ZDSCAL(nr,Sva(p),U(1,p),1)
               ENDDO
               CALL ZTRSM('L','U','N','N',nr,nr,CONE,Cwork(2*N+1),N,U,  &
     &                    Ldu)
!              .. apply the permutation from the second QR factorization
               DO q = 1 , nr
                  DO p = 1 , nr
                     Cwork(2*N+N*nr+nr+Iwork(N+p)) = U(p,q)
                  ENDDO
                  DO p = 1 , nr
                     U(p,q) = Cwork(2*N+N*nr+nr+p)
                  ENDDO
               ENDDO
               IF ( nr<N ) THEN
                  CALL ZLASET('A',N-nr,nr,CZERO,CZERO,V(nr+1,1),Ldv)
                  CALL ZLASET('A',nr,N-nr,CZERO,CZERO,V(1,nr+1),Ldv)
                  CALL ZLASET('A',N-nr,N-nr,CZERO,CONE,V(nr+1,nr+1),Ldv)
               ENDIF
               CALL ZUNMQR('L','N',N,N,nr,Cwork(2*N+1),N,Cwork(N+1),V,  &
     &                     Ldv,Cwork(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,  &
     &                     ierr)
            ELSE
!              Last line of defense.
! #:(          This is a rather pathological case: no scaled condition
!              improvement after two pivoted QR factorizations. Other
!              possibility is that the rank revealing QR factorization
!              or the condition estimator has failed, or the COND_OK
!              is set very close to ONE (which is unnecessary). Normally,
!              this branch should never be executed, but in rare cases of
!              failure of the RRQR or condition estimator, the last line of
!              defense ensures that ZGEJSV completes the task.
!              Compute the full SVD of L3 using ZGESVJ with explicit
!              accumulation of Jacobi rotations.
               CALL ZGESVJ('L','U','V',nr,nr,V,Ldv,Sva,nr,U,Ldu,        &
     &                     Cwork(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,Rwork,&
     &                     Lrwork,Info)
               scalem = Rwork(1)
               numrank = NINT(Rwork(2))
               IF ( nr<N ) THEN
                  CALL ZLASET('A',N-nr,nr,CZERO,CZERO,V(nr+1,1),Ldv)
                  CALL ZLASET('A',nr,N-nr,CZERO,CZERO,V(1,nr+1),Ldv)
                  CALL ZLASET('A',N-nr,N-nr,CZERO,CONE,V(nr+1,nr+1),Ldv)
               ENDIF
               CALL ZUNMQR('L','N',N,N,nr,Cwork(2*N+1),N,Cwork(N+1),V,  &
     &                     Ldv,Cwork(2*N+N*nr+nr+1),Lwork-2*N-N*nr-nr,  &
     &                     ierr)
!
               CALL ZUNMLQ('L','C',nr,nr,nr,Cwork(2*N+1),N,             &
     &                     Cwork(2*N+N*nr+1),U,Ldu,Cwork(2*N+N*nr+nr+1),&
     &                     Lwork-2*N-N*nr-nr,ierr)
               DO q = 1 , nr
                  DO p = 1 , nr
                     Cwork(2*N+N*nr+nr+Iwork(N+p)) = U(p,q)
                  ENDDO
                  DO p = 1 , nr
                     U(p,q) = Cwork(2*N+N*nr+nr+p)
                  ENDDO
               ENDDO
!
            ENDIF
!
!           Permute the rows of V using the (column) permutation from the
!           first QRF. Also, scale the columns to make them unit in
!           Euclidean norm. This applies to all cases.
!
            temp1 = SQRT(DBLE(N))*epsln
            DO q = 1 , N
               DO p = 1 , N
                  Cwork(2*N+N*nr+nr+Iwork(p)) = V(p,q)
               ENDDO
               DO p = 1 , N
                  V(p,q) = Cwork(2*N+N*nr+nr+p)
               ENDDO
               xsc = ONE/DZNRM2(N,V(1,q),1)
               IF ( (xsc<(ONE-temp1)) .OR. (xsc>(ONE+temp1)) )          &
     &              CALL ZDSCAL(N,xsc,V(1,q),1)
            ENDDO
!           At this moment, V contains the right singular vectors of A.
!           Next, assemble the left singular vector matrix U (M x N).
            IF ( nr<M ) THEN
               CALL ZLASET('A',M-nr,nr,CZERO,CZERO,U(nr+1,1),Ldu)
               IF ( nr<n1 ) THEN
                  CALL ZLASET('A',nr,n1-nr,CZERO,CZERO,U(1,nr+1),Ldu)
                  CALL ZLASET('A',M-nr,n1-nr,CZERO,CONE,U(nr+1,nr+1),   &
     &                        Ldu)
               ENDIF
            ENDIF
!
!           The Q matrix from the first QRF is built into the left singular
!           matrix U. This applies to all cases.
!
            CALL ZUNMQR('L','N',M,n1,N,A,Lda,Cwork,U,Ldu,Cwork(N+1),    &
     &                  Lwork-N,ierr)
 
!           The columns of U are normalized. The cost is O(M*N) flops.
            temp1 = SQRT(DBLE(M))*epsln
            DO p = 1 , nr
               xsc = ONE/DZNRM2(M,U(1,p),1)
               IF ( (xsc<(ONE-temp1)) .OR. (xsc>(ONE+temp1)) )          &
     &              CALL ZDSCAL(M,xsc,U(1,p),1)
            ENDDO
!
!           If the initial QRF is computed with row pivoting, the left
!           singular vectors must be adjusted.
!
            IF ( rowpiv ) CALL ZLASWP(n1,U,Ldu,1,M-1,Iwork(iwoff+1),-1)
!
         ELSE
!
!        .. the initial matrix A has almost orthogonal columns and
!        the second QRF is not needed
!
            CALL ZLACPY('U',N,N,A,Lda,Cwork(N+1),N)
            IF ( l2pert ) THEN
               xsc = SQRT(small)
               DO p = 2 , N
                  ctemp = xsc*Cwork(N+(p-1)*N+p)
                  DO q = 1 , p - 1
!                     CWORK(N+(q-1)*N+p)=-TEMP1 * ( CWORK(N+(p-1)*N+q) /
!     $                                        ABS(CWORK(N+(p-1)*N+q)) )
                     Cwork(N+(q-1)*N+p) = -ctemp
                  ENDDO
               ENDDO
            ELSE
               CALL ZLASET('L',N-1,N-1,CZERO,CZERO,Cwork(N+2),N)
            ENDIF
!
            CALL ZGESVJ('U','U','N',N,N,Cwork(N+1),N,Sva,N,U,Ldu,       &
     &                  Cwork(N+N*N+1),Lwork-N-N*N,Rwork,Lrwork,Info)
!
            scalem = Rwork(1)
            numrank = NINT(Rwork(2))
            DO p = 1 , N
               CALL ZCOPY(N,Cwork(N+(p-1)*N+1),1,U(1,p),1)
               CALL ZDSCAL(N,Sva(p),Cwork(N+(p-1)*N+1),1)
            ENDDO
!
            CALL ZTRSM('L','U','N','N',N,N,CONE,A,Lda,Cwork(N+1),N)
            DO p = 1 , N
               CALL ZCOPY(N,Cwork(N+p),N,V(Iwork(p),1),Ldv)
            ENDDO
            temp1 = SQRT(DBLE(N))*epsln
            DO p = 1 , N
               xsc = ONE/DZNRM2(N,V(1,p),1)
               IF ( (xsc<(ONE-temp1)) .OR. (xsc>(ONE+temp1)) )          &
     &              CALL ZDSCAL(N,xsc,V(1,p),1)
            ENDDO
!
!           Assemble the left singular vector matrix U (M x N).
!
            IF ( N<M ) THEN
               CALL ZLASET('A',M-N,N,CZERO,CZERO,U(N+1,1),Ldu)
               IF ( N<n1 ) THEN
                  CALL ZLASET('A',N,n1-N,CZERO,CZERO,U(1,N+1),Ldu)
                  CALL ZLASET('A',M-N,n1-N,CZERO,CONE,U(N+1,N+1),Ldu)
               ENDIF
            ENDIF
            CALL ZUNMQR('L','N',M,n1,N,A,Lda,Cwork,U,Ldu,Cwork(N+1),    &
     &                  Lwork-N,ierr)
            temp1 = SQRT(DBLE(M))*epsln
            DO p = 1 , n1
               xsc = ONE/DZNRM2(M,U(1,p),1)
               IF ( (xsc<(ONE-temp1)) .OR. (xsc>(ONE+temp1)) )          &
     &              CALL ZDSCAL(M,xsc,U(1,p),1)
            ENDDO
!
            IF ( rowpiv ) CALL ZLASWP(n1,U,Ldu,1,M-1,Iwork(iwoff+1),-1)
!
!
         ENDIF
!
!
!        end of the  >> almost orthogonal case <<  in the full SVD
!
         IF ( transp ) THEN
!           .. swap U and V because the procedure worked on A^*
            DO p = 1 , N
               CALL ZSWAP(N,U(1,p),1,V(1,p),1)
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
      Rwork(1) = uscal2*scalem
      Rwork(2) = uscal1
      IF ( errest ) Rwork(3) = sconda
      IF ( lsvec .AND. rsvec ) THEN
         Rwork(4) = condr1
         Rwork(5) = condr2
      ENDIF
      IF ( l2tran ) THEN
         Rwork(6) = entra
         Rwork(7) = entrat
      ENDIF
!
      Iwork(1) = nr
      Iwork(2) = numrank
      Iwork(3) = warning
      IF ( transp ) THEN
         Iwork(4) = 1
      ELSE
         Iwork(4) = -1
      ENDIF
 
!
!     ..
!     .. END OF ZGEJSV
!     ..
      END SUBROUTINE ZGEJSV
